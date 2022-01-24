"""
Module containing the necessary classes for the representation of a Waverider.
"""
import numpy as np
import matplotlib.pyplot as plt
import plotly.offline
import plotly.graph_objects as go
import triangle as tr
import Bezier as bzr
import Streamline
import aerodynamics
import config as cfg
from typing import List


class WRInputs:
    """
    Object with the inputs of a waverider
    Properties:
        b: Half cone shockwave angle
        s: Semispan of WR
        l: Length of WR
        per_l: Percentage of line segment of shochwave's base
        yle: Bezier control points of leading edge (y-axis)
        zle: Bezier control points of leading edge (z-axis)
        Cones: Dictionary of containing some necessary values calculated by the constrains check
        LE: The Bezier of the leading edge
        SW: The Bezier of the shochwave
    """
    b: float
    s: float
    l: float
    per_l: float
    yle: np.ndarray
    zle: np.ndarray
    LE: bzr.bezier
    SW: List[bzr.bezier]
    Cones: dict

    def __init__(self, b, s, l, per_l, yle, zle) -> None:
        self.b = float(b)
        self.s = float(s)
        self.l = float(l)
        self.per_l = float(per_l)
        self.yle = np.array(yle)
        self.zle = np.array(zle)
        self.Cones = None

        per_sw = np.array([0., 0.2, 0.5, 0.65, 1.])
        sz = len(per_sw)
        ysw = np.empty(sz)
        zsw = np.empty(sz)

        zsw[0:-1] = [np.tan(self.b*np.pi/180.)*self.l]*(sz-1)
        zsw[-1] = 0.7*zsw[0]
        for i, per in enumerate(per_sw):
            ysw[i] = self.s*(self.per_l + per*(1. - self.per_l))

        yle = [0., *self.yle, ysw[-1]]
        zle = [0., 0., *self.zle, zsw[-1]]

        self.SW = [bzr.bezier([0.0, ysw[0]], [zsw[0]] * 2, ['y', 'z']),
                bzr.bezier(ysw, zsw, ['y', 'z'])]

        self.LE = bzr.bezier(yle, zle, ['y', 'z'])

    def check(self) -> bool:
        # Performs the constraints check (True if geometry is unfeasible)

        def is_sorted(a):
            # Nested function to find if an array is sorted
            for i in range(a.size-1):
                if a[i+1] < a[i]:
                    return False
            return True

        minb = np.arcsin(1/cfg.MINF)*180/np.pi
        maxb = 30
        mins_l = 0.2
        maxs_l = 0.5
        minper_l = 0.2
        maxper_l = 1

        # b and s/l constraint
        if self.b < minb or self.b > maxb:
            return True
        elif self.s/self.l < mins_l or self.s/self.l > maxs_l:
            return True
        elif self.per_l < minper_l or self.per_l >= maxper_l:
            return True

        # First Leading Edge Point constraint
        if self.yle[0] < 0.1:
            return True

        # Checks if Normals and Leading-Edge intersect Shockwave
        yplane = np.linspace(0, self.s, cfg.PLANES)
        yplane_c = yplane[yplane > self.per_l*self.s]
        PLANES_C = len(yplane_c)
        PLANES_L = cfg.PLANES - PLANES_C
        tsw_c = self.SW[1].normal_inter(yplane_c, 'y')
        tsw_c[-1] = 1.
        
        tint = np.empty(cfg.PLANES)
        RC = np.empty(cfg.PLANES)
        dint = np.empty(cfg.PLANES)

        RC[:PLANES_L] = np.inf*np.ones([PLANES_L])
        RC[PLANES_L:] = self.SW[1].radius_curv(tsw_c)
        
        if any(self.LE.der1('y', np.linspace(0, 1, 100)) < 0):
            return True

        for i in range(cfg.PLANES-1):
            if i <= PLANES_L-1:
                tint[i] = self.LE.normal_inter(yplane[i], 'y')
                dint[i] = self.SW[0].z[0] - self.LE.curve('z', tint[i])
            else:
                j = i-PLANES_L
                k = np.tan(self.SW[1].normal_phi(tsw_c[j]))
                a = self.SW[1].normal_b(tsw_c[j])
                temp = self.LE.line_inter(k, a)
                # No intersection between normal and L.E. 
                if temp.size == 0:
                    return True
                # Intersection between L.E. and Shockwave curve
                if (self.LE.curve('z', temp) >= self.SW[1].curve('z', tsw_c[j]) or
                self.LE.curve('y', temp) > self.SW[1].curve('y', tsw_c[j])):
                    return True

                tint[i] = temp
                dint[i] = np.sqrt((self.LE.curve('y', tint[i]) - self.SW[1].curve('y', tsw_c[j])) ** 2 +\
                                (self.LE.curve('z', tint[i]) - self.SW[1].curve('z', tsw_c[j])) ** 2)
            
        tint[-1] = 1.   
        dint[-1] = 0.
        
        # Local cone vertex under LE constraint
        if any(dint > RC) or not(is_sorted(np.flip(dint))):
            return True

        # Minimum thickness constraint
        Y = 0.98*self.s
        tup = self.LE.normal_inter(Y, 'y')
        tdown = self.SW[1].normal_inter(Y, 'y')
        if tup.size == 0  or tdown.size == 0:
            return True
        thickness = np.abs(self.LE.curve('z', tup) - self.SW[1].curve('z', tdown))
        if thickness < 0.015*self.s:
            return True
        
        # Streamline start too close to cone vertex
        if any((RC - dint) < 0.05):
            return True

        self.Cones = {
            'tint': tint, 'dint': dint, 'RC': RC, 'yplane': yplane, 'yplane_c': yplane_c,
            'PLANES_C': PLANES_C, 'PLANES_L': PLANES_L, 'tsw_c': tsw_c
        }
        return False

    def reInit(self) -> 'WRInputs':
        return WRInputs(self.b, self.s, self.l, self.per_l, self.yle, self.zle)
    
    def __repr__(self) -> str:
        return f'WRInputs(b={self.b}, s={self.s}, l={self.l}, per_l={self.per_l}, yle={self.yle}, zle={self.zle})'

    def __str__(self) -> str:
        string = f'Shockwave Angle = {self.b:0.4} deg\nSemispan = {self.s:0.4} m\nLength = {self.l:0.4} m\n'
        string += f'Percentage of Line Segment = {self.per_l:0.4}\n'
        for i, y, z in zip(range(len(self.LE.y)), self.LE.y, self.LE.z):
            string += f'LE Point{i} = ({y:0.4},{z:0.4})\n'
        return string

    def __len__(self):
        return 4 + len(self.yle) + len(self.zle)

    def __iter__(self):
        return iter((self.b, self.s, self.l, self.per_l, *self.yle, *self.zle))

    def __getitem__(self, key):
        if isinstance(key, int):
            if 0 <= key < len(self):
                for i, x in enumerate(self):
                    if i == key:
                        return x
            elif -len(self) <= key < 0:
                for i, x in enumerate(self):
                    if i - len(self) == key:
                        return x
            else:
                raise IndexError('Index not in range.')
        elif isinstance(key, str):
            return getattr(self, key)
        else:
            raise IndexError('Index not integer or string.')

    def __setitem__(self, key, value):
        if isinstance(key, int):
            if key == 0:
                self.b = float(value)
            elif key == 1:
                self.s = float(value)
            elif key == 2:
                self.l = float(value)
            elif key == 3:
                self.per_l = float(value)
            elif 4 <= key < 4 + len(self.yle):
                self.yle[key - 4] = float(value)
            elif 4 + len(self.yle) <= key < len(self):
                self.zle[key - 4 - len(self.yle)] = float(value)
            elif -len(self) <= key < 0:
                self[len(self) + key] = float(value)
            else:
                raise IndexError('Index not in range.')
        elif isinstance(key, str):
            if key in vars(self):
                setattr(self, key, float(value))
            else:
                raise IndexError('No attribute with that name.')
        else:
            raise IndexError('Index not integer or string.')

    def __add__(self, other):
        if isinstance(other, WRInputs):
            b = self.b + other.b
            s = self.s + other.s
            l = self.l + other.l
            per_l = self.per_l + other.per_l
            yle = self.yle + other.yle
            zle = self.zle + other.zle
        elif isinstance(other, int):
            b = self.b + float(other) 
            s = self.s + float(other) 
            l = self.l + float(other)
            per_l = self.per_l + float(other)
            yle = self.yle + float(other) 
            zle = self.zle + float(other) 
        elif isinstance(other, float):
            b = self.b + other
            s = self.s + other
            l = self.l + other
            per_l = self.per_l + other
            yle = self.yle + other
            zle = self.zle + other
        else:
            raise TypeError('Not supported type. Must be same of instance int, float or WRInputs.')
        return WRInputs(b, s, l, per_l, yle, zle)

    def __sub__(self, other):
        if isinstance(other, WRInputs):
            b = self.b - other.b
            s = self.s - other.s
            l = self.l - other.l
            per_l = self.per_l - other.per_l
            yle = self.yle - other.yle
            zle = self.zle - other.zle
        elif isinstance(other, int):
            b = self.b - float(other) 
            s = self.s - float(other) 
            l = self.l - float(other) 
            per_l = self.per_l - float(other)
            yle = self.yle - float(other) 
            zle = self.zle - float(other) 
        elif isinstance(other, float):
            b = self.b - other
            s = self.s - other
            l = self.l - other
            per_l = self.per_l - other
            yle = self.yle - other
            zle = self.zle - other
        else:
            raise TypeError('Not supported type. Must be same of instance int, float or WRInputs.')
        return WRInputs(b, s, l, per_l, yle, zle)

    def __rsub__(self, other):
        if isinstance(other, WRInputs):
            b = other.b - self.b
            s = other.s - self.s
            l = other.l - self.l
            per_l = other.per_l - self.per_l
            yle = other.yle - self.yle
            zle = other.zle - self.zle
        elif isinstance(other, int):
            b = float(other) - self.b
            s = float(other) - self.s
            l = float(other) - self.l
            per_l = float(other) - self.per_l
            yle = float(other) - self.yle
            zle = float(other) - self.zle
        elif isinstance(other, float):
            b = other - self.b
            s = other - self.s
            l = other - self.l
            per_l = other - self.per_l
            yle = other - self.yle
            zle = other - self.zle
        else:
            raise TypeError('Not supported type. Must be same of instance int, float or WRInputs.')
        return WRInputs(b, s, l, per_l, yle, zle)

    def __mul__(self, other):
        if isinstance(other, WRInputs):
            b = self.b * other.b
            s = self.s * other.s
            l = self.l * other.l
            per_l = self.per_l * other.per_l
            yle = self.yle * other.yle
            zle = self.zle * other.zle
        elif isinstance(other, int):
            b = self.b * float(other) 
            s = self.s * float(other) 
            l = self.l * float(other)
            per_l = self.per_l * float(other)
            yle = self.yle * float(other) 
            zle = self.zle * float(other) 
        elif isinstance(other, float):
            b = self.b * other
            s = self.s * other
            l = self.l * other
            per_l = self.per_l * other
            yle = self.yle * other
            zle = self.zle * other
        else:
            raise TypeError('Not supported type. Must be same of instance int, float or WRInputs.')
        return WRInputs(b, s, l, per_l, yle, zle)

    def __truediv__(self, other):
        if isinstance(other, WRInputs):
            b = self.b / other.b
            s = self.s / other.s
            l = self.l / other.l
            per_l = self.per_l / other.per_l
            yle = self.yle / other.yle
            zle = self.zle / other.zle
        elif isinstance(other, int):
            b = self.b / float(other) 
            s = self.s / float(other) 
            l = self.l / float(other)
            per_l = self.per_l / float(other)
            yle = self.yle / float(other) 
            zle = self.zle / float(other) 
        elif isinstance(other, float):
            b = self.b / other
            s = self.s / other
            l = self.l / other
            per_l = self.per_l / other
            yle = self.yle / other
            zle = self.zle / other
        else:
            raise TypeError('Not supported type. Must be same of instance int, float or WRInputs.')
        return WRInputs(b, s, l, per_l, yle, zle)

    def __rtruediv__(self, other):
        if isinstance(other, WRInputs):
            b = other.b / self.b
            s = other.s / self.s
            l = other.l / self.l
            per_l = other.per_l / self.per_l
            yle = other.yle / self.yle
            zle = other.zle / self.zle
        elif isinstance(other, int):
            b = float(other) / self.b
            s = float(other)  / self.s
            l = float(other)  / self.l
            per_l = float(other) / self.per_l
            yle = float(other)  / self.yle
            zle = float(other)  / self.zle
        elif isinstance(other, float):
            b = other / self.b
            s = other / self.s
            l = other / self.l
            per_l = other / self.per_l
            yle = other / self.yle
            zle = other / self.zle
        else:
            raise TypeError('Not supported type. Must be same of instance int, float or WRInputs.')
        return WRInputs(b, s, l, per_l, yle, zle)
    __radd__ = __add__
    __iadd__ = __add__
    __isub__ = __sub__
    __rmul__ = __mul__
    __imul__ = __mul__
    __itruediv__ = __truediv__


class Waverider:
    """
    Object representing a waverider.
    The input is a WRInput object (see docstring of WRInput class)
    Additional properties:
        s_l: The semispan to length ratio
        nVar: The number of inputs
    """
    nVar: int
    b: float
    s: float
    l: float
    per_l: float
    yle: np.ndarray
    zle: np.ndarray
    LE: bzr.bezier
    SW: List[bzr.bezier]
    
    def __init__(self, wr_in):
        # set input atributes and curves definition
        self.nVar = len(wr_in)
        self.b = wr_in.b
        self.s = wr_in.s
        self.l = wr_in.l
        self.per_l = wr_in.per_l
        self.yle = wr_in.yle
        self.zle = wr_in.zle
        self.s_l = wr_in.s / wr_in.l
        self.LE = wr_in.LE
        self.SW = wr_in.SW
         
        # Distance if intersection and radius of curviture calculation
        if wr_in.Cones is None:
            yplane = np.linspace(0, self.s, cfg.PLANES)
            yplane_c = yplane[yplane > self.per_l*self.s]
            PLANES_C = len(yplane_c)
            PLANES_L = cfg.PLANES - PLANES_C

            tsw_c = self.SW[1].normal_inter(yplane_c, 'y')
            tsw_c[-1] = 1.
            tint = np.empty(cfg.PLANES)
            RC = np.empty(cfg.PLANES)
            dint = np.empty(cfg.PLANES)

            RC[:PLANES_L] = np.inf*np.ones([PLANES_L])
            RC[PLANES_L:] = self.SW[1].radius_curv(tsw_c)
            
            for i in range(cfg.PLANES-1):
                if i <= PLANES_L-1:
                    tint[i] = self.LE.normal_inter(yplane[i], 'y')
                    dint[i] = self.SW[0].z[0] - self.LE.curve('z', tint[i])
                else:
                    j = i-PLANES_L
                    k = np.tan(self.SW[1].normal_phi(tsw_c[j]))
                    b = self.SW[1].normal_b(tsw_c[j])
                    tint[i] = self.LE.line_inter(k, b)
                    dint[i] = np.sqrt((self.LE.curve('y', tint[i]) - self.SW[1].curve('y', tsw_c[j])) ** 2 +\
                                    (self.LE.curve('z', tint[i]) - self.SW[1].curve('z', tsw_c[j])) ** 2)
                
            tint[-1] = 1.
            dint[-1] = 0.
        else:
            yplane = wr_in.Cones['yplane']
            yplane_c = wr_in.Cones['yplane_c']
            PLANES_C = wr_in.Cones['PLANES_C']
            PLANES_L = wr_in.Cones['PLANES_L']
            tsw_c = wr_in.Cones['tsw_c']
            tint = wr_in.Cones['tint']
            RC = wr_in.Cones['RC']
            dint = wr_in.Cones['dint']

        # Taylor-Maccoll
        Taylor = Streamline.TaylorMaccoll(cfg.MINF, self.b)

        # Geometry Generation
        Dt = (dint[PLANES_L] / np.tan(self.b * np.pi / 180))/ \
            (np.sqrt(2 * cfg.CP * cfg.ATM.T +
             (cfg.MINF * cfg.ATM.v_sonic) ** 2) * Taylor.V[-1]) / cfg.N_POINTS
        
        self.LS = np.empty(cfg.PLANES, dtype=object)
        self.US = np.empty(cfg.PLANES, dtype=object)
        self.BS = np.empty(cfg.PLANES, dtype=object)

        for i in range(cfg.PLANES):
            if i < PLANES_L:
                self.LS[i] = Streamline.Wedge(cfg.MINF, self.b, dint[i], cfg.N_POINTS)
                ysw = yplane[i]
                zsw = self.SW[0].z[0]
                self.LS[i].coor_change(ysw, zsw)

            elif i >= PLANES_L and i != cfg.PLANES-1:
                self.LS[i] = Streamline.Cone(cfg.MINF, self.b, RC[i], dint[i], Taylor, Dt)
                ysw = yplane[i]
                zsw = self.SW[1].curve('z', tsw_c[i-PLANES_L])
                phi = self.SW[1].normal_phi(tsw_c[i-PLANES_L])
                self.LS[i].coor_change(ysw, zsw, phi)
                
            else:
                self.LS[i] = Streamline.Edge(Taylor, self.SW[1].y[-1], self.SW[1].z[-1])
            
            self.US[i] = Streamline.Upper(self.LS[i])
            self.BS[i] = Streamline.Base(self.LS[i])

        # Aerodynamics Run
        phi = (np.pi/2) * np.ones(cfg.PLANES)
        phi[PLANES_L:] = self.SW[1].normal_phi(tsw_c)
        aero = aerodynamics.Aero3d(self.LS, phi)
        self.L = aero.L
        self.D = aero.D
        self.V = aero.V
        self.S = aero.S
        self.CL = aero.L / \
            (0.5 * cfg.ATM.rho * aero.S * (cfg.MINF * cfg.ATM.v_sonic) ** 2)
        self.CD = aero.D / \
            (0.5 * cfg.ATM.rho * aero.S * (cfg.MINF * cfg.ATM.v_sonic) ** 2)
        self.volEff = aero.V / (aero.S ** 1.5)
        self.L_D = aero.L / aero.D

    def inputs(self) -> WRInputs:
        # returns an object of the inputs of the waverider
        return WRInputs(self.b, self.s, self.l, self.per_l, self.yle, self.zle)

    def __repr__(self) -> str:
        return f'Waverider({repr(self.inputs())})'

    def __str__(self) -> str:
        string = 'Waverider Object with:\n'
        # Design Parameters
        string += ' Design Parameters\n'
        string += f'\tMach Number: {cfg.MINF}\n'
        string += f'\tAltitude: {cfg.H} km\n'
        string += f'\tNumber of Planes: {cfg.PLANES}\n'
        string += f'\tNumber of Points: {cfg.N_POINTS}\n'
        # string += f'\tPercentage of Line Segment: {cfg.PER_L}\n'

        # Inputs
        string += ' Inputs:\n\t'
        string += '\t'.join(str(self.inputs()).splitlines(True))

        # Outputs
        string += ' Outputs:\n'
        string += f'\tLift: {self.L*1e-3:0.5} kN\n'
        string += f'\tDrag: {self.D*1e-3:0.5} kN\n'
        string += f'\tLift Coeff.: {self.CL:0.5}\n'
        string += f'\tDrag Coeff.: {self.CD:0.5}\n'
        string += f'\tL/D: {self.L_D:0.5}\n'
        string += f'\tVolume: {self.V:0.5} m^3\n'
        string += f'\tRef. Surface: {self.S:0.5} m^2\n'
        string += f'\tVolumetric Efficiency: {self.volEff:0.5}'

        return string

    def check(self):
        # Performs the constraints check (True if mission constraints are not met)
        minvolEff = 0.08
        maxvolEff = 0.25
        # minS = 0
        # maxS = 1e6
        minLength = 50
        maxLength = 80
        minL = 330585*cfg.ATM.g

        # Volumetric Efficiency Constraint
        if self.volEff > maxvolEff or self.volEff < minvolEff:
            return True
        
        # # Surface Constraint
        # if self.S > maxS or self.S < minS:
        #     return True
        
        # Length Constraint
        if self.l > maxLength or self.l < minLength:
            return True
        
        # Lift Constraint
        if self.L < minL:
            return True
        return False

    def plot3d(self, plotOffline=False):
        xls = []
        yls = []
        zls = []
        xus = []
        yus = []
        zus = []
        for ls, us in zip(self.LS, self.US):
            if ls.N != 1:
                xls.extend(list(ls.x))
                yls.extend(list(ls.y))
                zls.extend(list(ls.z))
                xus.extend(list(us.x))
                yus.extend(list(us.y))
                zus.extend(list(us.z))
            else:
                xls.append(float(ls.x))
                yls.append(float(ls.y))
                zls.append(float(ls.z))
                xus.append(float(us.x))
                yus.append(float(us.y))
                zus.append(float(us.z))

        seg = []
        for i in range(len(self.LS)-1):
            if self.LS[i+1].N == 1:
                seg.append([xls.index(self.LS[i].x[0]), xls.index(self.LS[i].x[0])+1])
            else:
                seg.append([xls.index(self.LS[i].x[0]), xls.index(self.LS[i+1].x[0])])
        for i in range(len(self.LS)-2, 0, -1):
            last = seg[-1][-1]
            seg.append([last, last-self.LS[i].N])
        seg.append([seg[-1][-1], 0])

        points_ls = np.stack((xls, yls), axis=1)
        A_ls = dict(vertices=points_ls, segments=seg, holes=[[self.LS[int(len(self.LS)/2)].x[0]*1.05, self.LS[int(len(self.LS)/2)].y[0]]])
        tri_ls = tr.triangulate(A_ls, 'p')
        points_us = np.stack((xus, yus), axis=1)
        A_us = dict(vertices=points_us, segments=seg, holes=[[self.LS[int(len(self.LS)/2)].x[0]*1.05, self.LS[int(len(self.LS)/2)].y[0]]])
        tri_us = tr.triangulate(A_us, 'p')

        data0 = go.Mesh3d(x=xls, y=yls, z=-np.array(zls), i=[i[0] for i in tri_ls['triangles']], 
                                                        j=[i[1] for i in tri_ls['triangles']], 
                                                        k=[i[2] for i in tri_ls['triangles']], color='red')
        data1 = go.Mesh3d(x=xus, y=yus, z=-np.array(zus), i=[i[0] for i in tri_us['triangles']], 
                                                        j=[i[1] for i in tri_us['triangles']], 
                                                        k=[i[2] for i in tri_us['triangles']], color='blue')
        data2 = go.Mesh3d(x=xls, y=-np.array(yls), z=-np.array(zls), i=[i[0] for i in tri_ls['triangles']], 
                                                        j=[i[1] for i in tri_ls['triangles']], 
                                                        k=[i[2] for i in tri_ls['triangles']], color='red')
        data3 = go.Mesh3d(x=xus, y=-np.array(yus), z=-np.array(zus), i=[i[0] for i in tri_us['triangles']], 
                                                        j=[i[1] for i in tri_us['triangles']], 
                                                        k=[i[2] for i in tri_us['triangles']], color='blue')
        
        if plotOffline:
        # Offline Plot (creating HTML file)
            fig = plotly.offline.plot({
                "data": [data0, data1, data2, data3], "layout": go.Layout(scene=dict(yaxis = dict(range=[-self.l/2, self.l/2]),
                                        zaxis = dict(range=[-self.l/2, self.l/2]),
                                        xaxis_visible=False, yaxis_visible=False, zaxis_visible=False),
                                        scene_aspectmode='cube')
            },
            filename= cfg.CONFIG['Viscous Method'] + ".html"
            )
        else:
        # Online Plot
            fig = go.Figure(data=[data0, data1, data2, data3])
            fig.update_layout(scene=dict(yaxis = dict(range=[-self.l/2, self.l/2]),
                                        zaxis = dict(range=[-self.l/2, self.l/2]),
                                        xaxis_visible=False, yaxis_visible=False, zaxis_visible=False),
                                        scene_aspectmode='cube')
            fig.show()

    def plotRadius(self, ax=None):
        if ax is None:
            _, ax = plt.subplots()
        
        self.LE.plot(ax, False ,'b')
        self.SW[0].plot(ax, False, 'k')
        self.SW[1].plot(ax, False, 'k', label=None)
        ax.axis('equal')
        ax.set_autoscale_on(False)

        yplane = np.linspace(0, self.s, len(self.LS))
        for y in yplane[[2*i for i in range(int(len(self.LS)/2))]]:
            if y < self.SW[1].y[0]:
                self.SW[0].plot_radius(ax, y=y)
            else:
                self.SW[1].plot_radius(ax, y=y)
        self.SW[1].plot_radius(ax, y=yplane[-1]) 

        m = np.max([self.s, self.SW[0].z[0]])*1.3
        ax.legend(['Leading Edge', 'Shockwave', 'Local Cone Radius'])
        ax.set_title('Radius of Curviture in Base Plane')
        ax.set_xlabel('y [m]')
        ax.set_ylabel('z [m]')
        ax.set(xlim=(-m/2, 3*m/2), ylim=(-m, m))
        # ax.set_xlim(-m, m)
        # ax.set_ylim(-m/2, 3*m/2)
        ax.invert_yaxis()
        ax.grid()
            
    def plotBase(self, ax=None):
        if ax is None:
            _, ax = plt.subplots()

        self.LE.plot(ax, False ,'b', label='Upper Surface')
        self.SW[0].plot(ax, False, '--k', label='Shockwave')
        self.SW[1].plot(ax, False ,'--k', label=None)

        self.LE.plot(ax, True ,'b', label=None)
        self.SW[0].plot(ax, True, '--k', label=None)
        self.SW[1].plot(ax, True ,'--k', label=None)

        ax.plot([i.y[-1] if i.N != 1 else i.y for i in self.LS],
                [i.z[-1] if i.N != 1 else i.z for i in self.LS ],
                'r', label='Lower Surface')
        ax.plot([-i.y[-1] if i.N != 1 else -i.y for i in self.LS],
                [i.z[-1] if i.N != 1 else i.z for i in self.LS ],
                'r', label=None)
        
        m = np.max([self.s, self.SW[0].z[0]])*1.3
        ax.legend()
        ax.set_title('Base Plane Projection')
        ax.set_xlabel('y [m]')
        ax.set_ylabel('z [m]')
        ax.axis('equal')
        ax.set_autoscaley_on(False)
        ax.set(xlim=(-m, m), ylim=(-m, m))
        ax.invert_yaxis()
        ax.grid()

    def plotTop(self, ax=None):
        if ax is None:
            _, ax = plt.subplots()
        
        y = [i.y[0] if i.N != 1 else i.y for i in self.LS]
        y_ = [-i.y[0] if i.N != 1 else -i.y for i in self.LS]
        x = [i.x[0] if i.N != 1 else i.x for i in self.LS]
        m = self.l*1.02

        ax.plot(y, x, 'k')
        ax.plot(y_, x, 'k')
        ax.set_title('Top View')
        ax.fill_between(y, x, color='darkgray')
        ax.fill_between(y_, x, color='darkgray')
        ax.axis('equal')
        # ax.set_ylim(0, m)
        # ax.set_xlim(-m/2, m/2)
        ax.set_ylabel('x [m]')
        ax.set_xlabel('y [m]')
        ax.grid()

    def plotAll(self, plot_3D=True, offline3D=False):
        # Creating Axes
        _, (axBase, axRad) = plt.subplots(1, 2, sharey=True, figsize=(11, 5))
        mananger1 = plt.get_current_fig_manager()
        _, axTop = plt.subplots(figsize=(5.5, 5))
        mananger2 = plt.get_current_fig_manager()
        mananger1.set_window_title(cfg.CONFIG['Viscous Method'])
        mananger2.set_window_title(cfg.CONFIG['Viscous Method'])

        self.plotTop(axTop)
        self.plotBase(axBase)
        self.plotRadius(axRad)
        if plot_3D:
            self.plot3d(offline3D)  

        axRad.set_ylabel(None)

    def todict(self) -> dict:
        """Method that returns a dictionary with the inputs of the waverider (for saving purposes)"""
        return {
            'b': self.b, 's': self.s, 'l': self.l, 'per_l': self.per_l,
            'yle': self.yle.tolist(), 'zle': self.zle.tolist()
                }


if __name__ == "__main__":
    def testing() -> None:
        import pickle
        from config_modify import config_create_window
        import time

        config_create_window()
        NWR = 1
        tic = time.time()
        for _ in range(NWR):
            b = 21.48
            s = 33.73
            l = 72.4
            per_l = 0.3719
            yle = np.array([4.258, 6.793, 20.39])
            zle = np.array([19.99, 17.37])
            guess = WRInputs(b, s, l, per_l, yle, zle)
            guess = Waverider(guess)
        toc = time.time()

        print(f'Time per Waverider = {(toc-tic)/NWR:0.4} sec')
        print(guess)
        print(repr(guess))

        guess.plotAll(plot_3D=False)
        plt.show()

        # filename = "WaveriderGuesses.pickle"
        # with open(filename, 'rb') as file:
        #     guesses = pickle.load(file)
        # print(guesses[49])

    testing()
