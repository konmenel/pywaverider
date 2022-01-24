"""Module containing the necessary classes for the representation of a Waverider."""
from functools import cached_property
import numpy as np
import matplotlib.pyplot as plt
import plotly.offline
import plotly.graph_objects as go
import triangle as tr
import bezier as bzr
import streamline as sln
import aerodynamics as aerod
import config as cfg
from typing import List, Union



class WRInputs:
    """Object with the inputs of a waverider
    Properties:
        b: Half cone shockwave angle [deg]
        s: Semispan of WR [m]
        l: Length of WR [m]
        per_l: Percentage of line segment of shochwave's base
        yle: Bezier control points of leading edge (y-axis) [m]
        zle: Bezier control points of leading edge (z-axis) [m]
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
    LE: bzr.Bezier
    SW: List[bzr.Bezier]
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

        zsw[0:-1] = [np.tan(self.b*np.pi/180.) * self.l] * (sz-1)
        zsw[-1] = 0.7 * zsw[0]
        for i, per in enumerate(per_sw):
            ysw[i] = self.s * (self.per_l + per * (1. - self.per_l))

        yle = [0., *self.yle, ysw[-1]]
        zle = [0., 0., *self.zle, zsw[-1]]

        self.SW = [bzr.Bezier([0.0, ysw[0]], [zsw[0]] * 2, ['y', 'z']),
                bzr.Bezier(ysw, zsw, ['y', 'z'])]

        self.LE = bzr.Bezier(yle, zle, ['y', 'z'])

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
        minLength = 50
        maxLength = 80

        # Shockwave angle Constraint
        if self.b < minb or self.b > maxb:
            return True
        
        # Semispan to Length Ratio Constraint
        if self.s/self.l < mins_l or self.s/self.l > maxs_l:
            return True
        
        # Shockwave Line Segment Contraint
        if self.per_l < minper_l or self.per_l >= maxper_l:
            return True

        # Length Constraint
        if self.l > maxLength or self.l < minLength:
            return True

        # First Leading Edge Point Constraint
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
        string = f'Shockwave Angle = {self.b:0.4} deg\n' + \
            f'Semispan = {self.s:0.4} m\n' + \
            f'Length = {self.l:0.4} m\n' + \
            f's/l = {self.s/self.l:0.4}\n' + \
            f'Percentage of Line Segment = {self.per_l:0.4}\n'
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
    """Object representing a waverider.
    The input is a WRInput object (see docstring of WRInput class)
    Additional properties:
        s_l: The semispan to length ratio
        nVar: The number of inputs
        LS: Numpy array of objects containing the points of the lower surface.
        US: Numpy array of objects containing the points of the upper surface.
        BS: Numpy array of objects containing the points of the base surface.
        L: The Lift of the waverider [N]
        D: The Drag of the waverider [N]
        V: The Volume of the waverider [m3]
        S: The reference Surface of the waverider [m2]
        CL: The coefficient of Lift of the waverider (L / (q*S))
        CD: The coefficient of Drag of the waverider (D / (q*S))
        volEff: The Volumetric Efficiency of the waverider (V / S^(1.5))
        L_D: The Lift to Drag ratio of the waverider (L / D)

    Method properties:
        L_loss: Waverider lift loss due to viscous interactions [N]
        D_gain: Waverider drag gain due to viscous interactions [N]
    """
    # Input related properties
    nVar: int
    b: float
    s: float
    l: float
    s_l: float
    per_l: float
    yle: np.ndarray
    zle: np.ndarray
    LE: bzr.Bezier
    SW: List[bzr.Bezier]
    Cones_cache: dict

    # Output related properties
    LS: np.ndarray
    US: np.ndarray
    BS: np.ndarray
    L: float
    D: float
    V: float
    S: float
    CL: float
    CD: float
    volEff: float
    L_D: float

    # Method properties
    L_loss: float
    D_gain: float
    Cp_base: float
    P_base: float
    S_base: float
    D_base: float

    def __init__(self, wr_in: Union[WRInputs, dict]):
        # if the inputs are given in dictionary format unpack them and
        # create the WRInputs object.
        if isinstance(wr_in, dict):
            wr_in = WRInputs(**wr_in)

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
            yplane_c = yplane[yplane > self.per_l * self.s]
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
            self.Cones_cache = {
                'tint': tint, 'dint': dint, 'RC': RC, 'yplane': yplane, 'yplane_c': yplane_c,
                'PLANES_C': PLANES_C, 'PLANES_L': PLANES_L, 'tsw_c': tsw_c
            }
        else:
            yplane = wr_in.Cones['yplane']
            yplane_c = wr_in.Cones['yplane_c']
            PLANES_C = wr_in.Cones['PLANES_C']
            PLANES_L = wr_in.Cones['PLANES_L']
            tsw_c = wr_in.Cones['tsw_c']
            tint = wr_in.Cones['tint']
            RC = wr_in.Cones['RC']
            dint = wr_in.Cones['dint']
            self.Cones_cache = wr_in.Cones

        # Taylor-Maccoll
        Taylor = sln.TaylorMaccoll(cfg.MINF, self.b)

        # Geometry Generation
        Dt = (dint[PLANES_L] / np.tan(self.b * np.pi / 180))/ \
            (np.sqrt(2 * cfg.CP * cfg.ATM.T +
             (cfg.MINF * cfg.ATM.v_sonic) ** 2) * Taylor.V[-1]) / cfg.N_POINTS
        
        self.LS = np.empty(cfg.PLANES, dtype=object)
        self.US = np.empty(cfg.PLANES, dtype=object)
        self.BS = np.empty(cfg.PLANES, dtype=object)

        for i in range(cfg.PLANES):
            if i < PLANES_L:
                self.LS[i] = sln.Wedge(cfg.MINF, self.b, dint[i], cfg.N_POINTS)
                ysw = yplane[i]
                zsw = self.SW[0].z[0]
                self.LS[i].coor_change(ysw, zsw)

            elif i >= PLANES_L and i != cfg.PLANES-1:
                self.LS[i] = sln.Cone(cfg.MINF, self.b, RC[i], dint[i], Taylor, Dt)
                ysw = yplane[i]
                zsw = self.SW[1].curve('z', tsw_c[i-PLANES_L])
                phi = self.SW[1].normal_phi(tsw_c[i-PLANES_L])
                self.LS[i].coor_change(ysw, zsw, phi)
                
            else:
                self.LS[i] = sln.Edge(Taylor, self.SW[1].y[-1], self.SW[1].z[-1])
            
            self.US[i] = sln.Upper(self.LS[i])
            self.BS[i] = sln.Base(self.LS[i])

        # Aerodynamics Run
        phi = (np.pi/2) * np.ones(cfg.PLANES)
        phi[PLANES_L:] = self.SW[1].normal_phi(tsw_c)
        aero = aerod.Aero3d(self, phi)
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

    @cached_property
    def L_loss(self) -> float:
        """Function to calculate the Lift loss."""
        # Code to calculate the angle of each plane (phi)
        yplane = np.linspace(0, self.s, cfg.PLANES)
        yplane_c = yplane[yplane > self.per_l * self.s]
        PLANES_L = cfg.PLANES - len(yplane_c)
        tsw_c = self.SW[1].normal_inter(yplane_c, 'y')
        tsw_c[-1] = 1.

        phi = (np.pi/2) * np.ones(cfg.PLANES)
        phi[PLANES_L:] = self.SW[1].normal_phi(tsw_c)
        
        return aerod.Aero3d.lift_loss(self, phi)

    @cached_property
    def D_gain(self) -> float:
        """Function to calculate the Drag gain."""        
        return aerod.Aero3d.drag_gain(self)

    @cached_property
    def Cp_base(self) -> float:
        """Return the coefficient of pressure (Cp) on the base of waverider."""
        cp = (2. / (cfg.GAM * cfg.MINF ** 2)) * (((2. / (cfg.GAM + 1)) ** 1.4) * \
                ((1 / cfg.MINF) ** 2.8) * ((2 * cfg.GAM * cfg.MINF ** 2 - (cfg.GAM - 1)) / \
                (cfg.GAM + 1)) - 1)
        return cp

    @cached_property
    def P_base(self) -> float:
        """Return the pressure on the base of waverider."""
        P, _ = aerod.base_surf_press(self.LS)
        return P

    @cached_property
    def S_base(self) -> float:
        """Return the surface of the base of waverider."""
        _, S = aerod.base_surf_press(self.LS)
        return S 

    @cached_property
    def D_base(self) -> float:
        """Return the base drag ((Pinf - Pbase) * surface) of waverider."""
        p_base, s_base = aerod.base_surf_press(self.LS)
        return (cfg.ATM.P - p_base) * s_base

    @cached_property
    def _boundaryIndeces(self) -> List[int]:
        seg = [(0, self.LS[0].N)]
        # Leading edge loop
        for i in range(1, len(self.LS) - 1):
            seg.append((seg[-1][1], seg[-1][1] + self.LS[i].N))
        # Tailing edge loop
        for i in range(len(self.LS) - 1, 0, -1):
            last = seg[-1][1]
            seg.append((last, last - self.LS[i].N))\
        # Symmetry Plane
        seg.extend([(i, i - 1) for i in range(seg[-1][1], 0, -1)])

        return seg

    def inputs(self) -> WRInputs:
        # returns an object of the inputs of the waverider
        return WRInputs(self.b, self.s, self.l, self.per_l, self.yle, self.zle)

    def __repr__(self) -> str:
        return f'Waverider({repr(self.inputs())})'

    def __str__(self) -> str:
        # Design Parameters
        string = '  Design Parameters\n' + \
            f'\tMach Number: {cfg.MINF}\n' + \
            f'\tAltitude: {cfg.H} km\n' + \
            f'\tMethod for Viscous Effects: {cfg.CONFIG["Viscous Method"]}\n' + \
            f'\tNumber of Planes: {cfg.PLANES}\n' + \
            f'\tNumber of Points: {cfg.N_POINTS}\n'
            
        # string += f'\tPercentage of Line Segment: {cfg.PER_L}\n'

        # Inputs
        string += ' Inputs:\n\t' + "\t".join(str(self.inputs()).splitlines(True))

        # Outputs
        string += '  Outputs:\n' + \
            f'\tLift: {self.L*1e-3:0.5} kN\n' + \
            f'\tDrag: {self.D*1e-3:0.5} kN\n' + \
            f'\tLift Coeff.: {self.CL:0.5}\n' + \
            f'\tDrag Coeff.: {self.CD:0.5}\n' + \
            f'\tL/D: {self.L_D:0.5}\n' + \
            f'\tVolume: {self.V:0.5} m^3\n' + \
            f'\tRef. Surface: {self.S:0.5} m^2\n' + \
            f'\tVolumetric Efficiency: {self.volEff:0.5}'

        return string

    def check(self):
        # Performs the constraints check (True if mission constraints are not met)
        minvolEff = 0.08
        maxvolEff = 0.25
        # minS = 0
        # maxS = 1e6
        minL = 330585*cfg.ATM.g

        # Volumetric Efficiency Constraint
        if self.volEff > maxvolEff or self.volEff < minvolEff:
            return True
        
        # # Surface Constraint
        # if self.S > maxS or self.S < minS:
        #     return True
        
        # Lift Constraint
        if self.L < minL:
            return True
        return False

    def getAllPlaneDataInLine(self, surface: str, attr: str) -> List[float]:
        """Method to get the data of all planes in one continuous list. 
        (i.e. LS.x of all planes)

        Parameters
        ----------
        surface : str
            [description]
        atr : str
            [description]

        Returns
        -------
        List
            [description]
        """
        surface_map = {'upper': self.US, 'lower': self.LS, 'base': self.BS}
        data = []
        for surf in  surface_map[surface]:
            if surf.N != 1:
                data.extend(list(getattr(surf, attr)))
            else:
                data.append(float(getattr(surf, attr)))

        return data

    def getAllPlaneShearInLine(self, surface: str) -> List[float]:
        """Method to get the data of all planes in one continuous list. 
        (i.e. LS.x of all planes)

        Parameters
        ----------
        surface : str
            [description]

        Returns
        -------
        List
            [description]
        """
        surface_map = {'lower': 0, 'upper': 1}
        tw = []
        for i, plane in enumerate(self.LS):
            if i < cfg.PLANES - 1:
                length = np.abs(plane.x[0] - plane.x[-1])
                x_plane = length - plane.x
            else:
                x_plane = plane.x
            
            if cfg.METHOD == 1:
                tw_for_both = aerod.van_driest_method(plane, x_plane)
            elif cfg.METHOD == 2:
                tw_for_both = aerod.ref_temp_method(plane, x_plane)

            tw.extend(tw_for_both[surface_map[surface]])

            tw[-1] = tw[-self.LS[-2].N]

        return tw

    def plotTempOnSymmetry(self, ax=None, ref=False):
        """Plot the temperature of the waverider on the symmetry.
        Arguments:
            ax: the matplotlib.axes.Axes object to plot on. If not given the
            function will create one and return it
            ref: boolean. If true plot the reference temperature used for the viscous
            calculations (default=False)
        """
        return_flag = False
        if ax == None:
            _, ax = plt.subplots()
            return_flag = True

        T = self.LS[0].T

        if ref:
            T0 = T * (1. + ((cfg.GAM - 1) / 2) * self.LS[0].M ** 2)
            Tw = T + 0.88 * (T0 - T)
            T = T * (1. + 0.032 * self.LS[0].M ** 2 + 0.58 * (Tw / T - 1.))
        
        ax.plot(self.LS[0].x , T)

        title = 'Temperature on Symmetry'
        if ref:
            title = f'Reference {title}'

        ax.set_title(title)
        ax.set_ylabel('Temperature [K]')
        ax.set_xlabel('x [m]')
        ax.grid()

        if return_flag:
            return ax

    def plotStress(self, surface='lower', ax=None, plot3d=False):
        """Plots the wall shear stress on top on the surface.
        Arguments:
            surface: a string for to specify the surface (`upper` or `lower`)
            ax: the matplotlib.axes.Axes object to plot on. If not given the
            function will create one and return it
        """
        if not cfg.VISCOUS:
            print('No shear stress exist for Inviscid Case.')
            return

        # For plotting
        levels = 100
        colormap = 'jet'
        return_flag = False

        # Indeces of Waverider bountaries
        seg = self._boundaryIndeces
        
        # 3D plot
        if plot3d:
            xls = self.getAllPlaneDataInLine('lower', 'x')
            yls = self.getAllPlaneDataInLine('lower', 'y')
            zls = self.getAllPlaneDataInLine('lower', 'z')
            tls = self.getAllPlaneShearInLine('lower')
            xus = self.getAllPlaneDataInLine('upper', 'x')
            yus = self.getAllPlaneDataInLine('upper', 'y')
            zus = self.getAllPlaneDataInLine('upper', 'z')
            tus = self.getAllPlaneShearInLine('upper')

            points_ls = np.stack((xls, yls), axis=1)
            A_ls = dict(
                vertices=points_ls,
                segments=seg,
                holes=[(self.LS[int(len(self.LS)/2)].x[0]*1.05, self.LS[int(len(self.LS)/2)].y[0])]
            )
            tri_ls = tr.triangulate(A_ls, 'p')
            points_us = np.stack((xus, yus), axis=1)
            A_us = dict(
                vertices=points_us,
                segments=seg,
                holes=[(self.LS[int(len(self.LS)/2)].x[0]*1.05, self.LS[int(len(self.LS)/2)].y[0])]
            )
            tri_us = tr.triangulate(A_us, 'p')

            data0 = go.Mesh3d(x=xls, y=yls, z=-np.array(zls),
                            i=[i[0] for i in tri_ls['triangles']],
                            j=[i[1] for i in tri_ls['triangles']],
                            k=[i[2] for i in tri_ls['triangles']],
                            intensity=tls, colorscale=colormap)
            data1 = go.Mesh3d(x=xus, y=yus, z=-np.array(zus),
                            i=[i[0] for i in tri_us['triangles']],
                            j=[i[1] for i in tri_us['triangles']],
                            k=[i[2] for i in tri_us['triangles']],
                            intensity=tus, colorscale=colormap)
            data2 = go.Mesh3d(x=xls, y=-np.array(yls), z=-np.array(zls),
                            i=[i[0] for i in tri_ls['triangles']],
                            j=[i[1] for i in tri_ls['triangles']],
                            k=[i[2] for i in tri_ls['triangles']],
                            intensity=tls, colorscale=colormap)
            data3 = go.Mesh3d(x=xus, y=-np.array(yus), z=-np.array(zus),
                            i=[i[0] for i in tri_us['triangles']],
                            j=[i[1] for i in tri_us['triangles']],
                            k=[i[2] for i in tri_us['triangles']],
                            intensity=tus, colorscale=colormap)

            fig = go.Figure(data=[data0, data1, data2, data3])
            fig.update_layout(scene=dict(yaxis = dict(range=[-self.l/2, self.l/2]),
                                        zaxis = dict(range=[-self.l/2, self.l/2]),
                                        xaxis_visible=False, yaxis_visible=False, zaxis_visible=False),
                                        scene_aspectmode='cube')
            fig.show()
            return

        # 2D plot
        if ax == None:
            _, ax = plt.subplots()
            return_flag = True

        x = self.getAllPlaneDataInLine(surface, 'x')
        y = self.getAllPlaneDataInLine(surface, 'y')

        tw = self.getAllPlaneShearInLine(surface)
        # for i, plane in enumerate(self.LS):
        #     if i < cfg.PLANES - 1:
        #         length = np.abs(plane.x[0] - plane.x[-1])
        #         x_plane = length - plane.x
        #     else:
        #         x_plane = plane.x
            
        #     if cfg.METHOD == 1:
        #         tw_ls, tw_us = aerod.van_driest_method(plane, x_plane)
        #     elif cfg.METHOD == 2:
        #         tw_ls, tw_us = aerod.ref_temp_method(plane, x_plane)

        #     if surface == 'lower':
        #         tw.extend(list(tw_ls))
        #     elif surface == 'upper':
        #         tw.extend(list(tw_us))

        #     tw[-1] = tw[-self.LS[-2].N]

        # Plotting
        vert = np.stack((x, y), axis=1)
        tri = dict(
            vertices=vert,
            segments=seg,
            holes=[(self.LS[int(len(self.LS)/2)].x[0]*1.05, self.LS[int(len(self.LS)/2)].y[0])]
            
        )
        tri = tr.triangulate(tri, 'p')
        ax.tricontourf(y, x, tri['triangles'], tw, levels, cmap=colormap)
        cntr = ax.tricontourf([-i for i in y], x, tri['triangles'], tw, levels, cmap=colormap)
        
        cbar = plt.gcf().colorbar(cntr, ax=ax)
        cbar.set_label('Shear Stress [N/m^2]')
        ax.set_title(f'{surface.capitalize()} Surface Shear Stress Contour')
        ax.axis('equal')
        ax.set_ylabel('x [m]')
        ax.set_xlabel('y [m]')
        ax.grid()

        if return_flag:
            return ax

    def plotPressure(self, surface='lower', ax=None, plot3d=False):
        """Plots the pressure on top on the surface.
        Arguments:
            surface: a string for to specify the surface (`upper` or `lower`)
            ax: the matplotlib.axes.Axes object to plot on. If not given the
            function will create one and return it
        """

        # For Plotting
        levels = 100
        colormap = 'jet'
        return_flag = False

        # Indeces of Waverider bountaries
        seg = self._boundaryIndeces

        # 3D plot
        if plot3d:
            xls = self.getAllPlaneDataInLine('lower', 'x')
            yls = self.getAllPlaneDataInLine('lower', 'y')
            zls = self.getAllPlaneDataInLine('lower', 'z')
            Pls = self.getAllPlaneDataInLine('lower', 'P')
            xus = self.getAllPlaneDataInLine('upper', 'x')
            yus = self.getAllPlaneDataInLine('upper', 'y')
            zus = self.getAllPlaneDataInLine('upper', 'z')
            Pus = self.getAllPlaneDataInLine('upper', 'P')

            points_ls = np.stack((xls, yls), axis=1)
            A_ls = dict(
                vertices=points_ls,
                segments=seg,
                holes=[(self.LS[int(len(self.LS)/2)].x[0]*1.05, self.LS[int(len(self.LS)/2)].y[0])]
            )
            tri_ls = tr.triangulate(A_ls, 'p')
            points_us = np.stack((xus, yus), axis=1)
            A_us = dict(
                vertices=points_us,
                segments=seg,
                holes=[(self.LS[int(len(self.LS)/2)].x[0]*1.05, self.LS[int(len(self.LS)/2)].y[0])]
            )
            tri_us = tr.triangulate(A_us, 'p')

            data0 = go.Mesh3d(x=xls, y=yls, z=-np.array(zls),
                            i=[i[0] for i in tri_ls['triangles']],
                            j=[i[1] for i in tri_ls['triangles']],
                            k=[i[2] for i in tri_ls['triangles']],
                            intensity=Pls, colorscale=colormap)
            data1 = go.Mesh3d(x=xus, y=yus, z=-np.array(zus),
                            i=[i[0] for i in tri_us['triangles']],
                            j=[i[1] for i in tri_us['triangles']],
                            k=[i[2] for i in tri_us['triangles']],
                            intensity=Pus, colorscale=colormap)
            data2 = go.Mesh3d(x=xls, y=-np.array(yls), z=-np.array(zls),
                            i=[i[0] for i in tri_ls['triangles']],
                            j=[i[1] for i in tri_ls['triangles']],
                            k=[i[2] for i in tri_ls['triangles']],
                            intensity=Pls, colorscale=colormap)
            data3 = go.Mesh3d(x=xus, y=-np.array(yus), z=-np.array(zus),
                            i=[i[0] for i in tri_us['triangles']],
                            j=[i[1] for i in tri_us['triangles']],
                            k=[i[2] for i in tri_us['triangles']],
                            intensity=Pus, colorscale=colormap)

            fig = go.Figure(data=[data0, data1, data2, data3])
            fig.update_layout(scene=dict(yaxis = dict(range=[-self.l/2, self.l/2]),
                                        zaxis = dict(range=[-self.l/2, self.l/2]),
                                        xaxis_visible=False, yaxis_visible=False, zaxis_visible=False),
                                        scene_aspectmode='cube')
            fig.show()
            return

        # 2D plot
        if ax == None:
            _, ax = plt.subplots()
            return_flag = True

        x = self.getAllPlaneDataInLine(surface, 'x')
        y = self.getAllPlaneDataInLine(surface, 'y')
        pressure = self.getAllPlaneDataInLine(surface, 'P')
        # if surface == 'lower':
        #     for plane in self.LS:
        #         x.extend(plane.x)
        #         y.extend(plane.y)
        #         pressure.extend(plane.P)
        # elif surface == 'upper':
        #     for plane in self.US:
        #         x.extend(plane.x)
        #         y.extend(plane.y)
        #         pressure.extend(plane.P)
        # else:
        #     raise ValueError(f'Unknown Surface "{surface}". Surface variable should be either "upper" or "lower".')
        


        vert = np.stack((x, y), axis=1)
        tri = dict(
            vertices=vert,
            segments=seg,
            holes=[(self.LS[int(len(self.LS)/2)].x[0]*1.05, self.LS[int(len(self.LS)/2)].y[0])]
            
        )
        tri = tr.triangulate(tri, 'p')

        ax.tricontourf(y, x, tri['triangles'], pressure, levels, cmap=colormap)
        cntr = ax.tricontourf([-i for i in y], x, tri['triangles'], pressure, levels, cmap=colormap)
        
        cbar = plt.gcf().colorbar(cntr, ax=ax)
        cbar.set_label('Pressure [Pa]')
        ax.set_title(f'{surface.capitalize()} Surface Pressure Contour')
        ax.axis('equal')
        ax.set_ylabel('x [m]')
        ax.set_xlabel('y [m]')
        ax.grid()

        if return_flag:
            return ax

    def plotTemperature(self, surface='lower', ax=None, ref=False, plot3d=False):
        """Plots the pressure on top on the surface.
        Arguments:
            surface: a string for to specify the surface (`upper` or `lower`)
            ax: the matplotlib.axes.Axes object to plot on. If not given the
            function will create one and return it.
            ref: boolean. If true plot the reference temperature used for the viscous
            calculations (default=False).
        """
        # For plotting
        levels = 100
        colormap = 'jet'
        return_flag = False
        
        # Indeces of Waverider bountaries
        seg = self._boundaryIndeces
        
        # 3D plot
        if plot3d:
            xls = self.getAllPlaneDataInLine('lower', 'x')
            yls = self.getAllPlaneDataInLine('lower', 'y')
            zls = self.getAllPlaneDataInLine('lower', 'z')
            Tls = self.getAllPlaneDataInLine('lower', 'T')
            xus = self.getAllPlaneDataInLine('upper', 'x')
            yus = self.getAllPlaneDataInLine('upper', 'y')
            zus = self.getAllPlaneDataInLine('upper', 'z')
            Tus = self.getAllPlaneDataInLine('upper', 'T')

            if ref:
                Mls = self.getAllPlaneDataInLine('lower', 'M')
                Tls = np.array(Tls)
                Mls = np.array(Mls)
                T0 = Tls * (1. + ((cfg.GAM - 1) / 2) * Mls ** 2)
                Tw = Tls + 0.88 * (T0 - Tls)
                Tls = Tls * (1. + 0.032 * Mls ** 2 + 0.58 * (Tw / Tls - 1.))
                
                Mus = self.getAllPlaneDataInLine('upper', 'M')
                Tus = np.array(Tus)
                Mus = np.array(Mus)
                T0 = Tus * (1. + ((cfg.GAM - 1) / 2) * Mus ** 2)
                Tw = Tus + 0.88 * (T0 - Tus)
                Tus = Tus * (1. + 0.032 * Mus ** 2 + 0.58 * (Tw / Tus - 1.))

            points_ls = np.stack((xls, yls), axis=1)
            A_ls = dict(
                vertices=points_ls,
                segments=seg,
                holes=[(self.LS[int(len(self.LS)/2)].x[0]*1.05, self.LS[int(len(self.LS)/2)].y[0])]
            )
            tri_ls = tr.triangulate(A_ls, 'p')
            points_us = np.stack((xus, yus), axis=1)
            A_us = dict(
                vertices=points_us,
                segments=seg,
                holes=[(self.LS[int(len(self.LS)/2)].x[0]*1.05, self.LS[int(len(self.LS)/2)].y[0])]
            )
            tri_us = tr.triangulate(A_us, 'p')

            data0 = go.Mesh3d(x=xls, y=yls, z=-np.array(zls),
                            i=[i[0] for i in tri_ls['triangles']],
                            j=[i[1] for i in tri_ls['triangles']],
                            k=[i[2] for i in tri_ls['triangles']],
                            intensity=Tls, colorscale=colormap)
            data1 = go.Mesh3d(x=xus, y=yus, z=-np.array(zus),
                            i=[i[0] for i in tri_us['triangles']],
                            j=[i[1] for i in tri_us['triangles']],
                            k=[i[2] for i in tri_us['triangles']],
                            intensity=Tus, colorscale=colormap)
            data2 = go.Mesh3d(x=xls, y=-np.array(yls), z=-np.array(zls),
                            i=[i[0] for i in tri_ls['triangles']],
                            j=[i[1] for i in tri_ls['triangles']],
                            k=[i[2] for i in tri_ls['triangles']],
                            intensity=Tls, colorscale=colormap)
            data3 = go.Mesh3d(x=xus, y=-np.array(yus), z=-np.array(zus),
                            i=[i[0] for i in tri_us['triangles']],
                            j=[i[1] for i in tri_us['triangles']],
                            k=[i[2] for i in tri_us['triangles']],
                            intensity=Tus, colorscale=colormap)

            fig = go.Figure(data=[data0, data1, data2, data3])
            fig.update_layout(scene=dict(yaxis = dict(range=[-self.l/2, self.l/2]),
                                        zaxis = dict(range=[-self.l/2, self.l/2]),
                                        xaxis_visible=False, yaxis_visible=False, zaxis_visible=False),
                                        scene_aspectmode='cube')
            fig.show()
            return

        # 2D plot
        if ax == None:
            _, ax = plt.subplots()
            return_flag = True


        x = self.getAllPlaneDataInLine(surface, 'x')
        y = self.getAllPlaneDataInLine(surface, 'y')
        T = self.getAllPlaneDataInLine(surface, 'T')
        
        if ref:
            M = self.getAllPlaneDataInLine(surface, 'M')
            T = np.array(T)
            M = np.array(M)
            T0 = T * (1. + ((cfg.GAM - 1) / 2) * M ** 2)
            Tw = T + 0.88 * (T0 - T)
            T = T * (1. + 0.032 * M ** 2 + 0.58 * (Tw / T - 1.))

        vert = np.stack((x, y), axis=1)
        tri = dict(
            vertices=vert,
            segments=seg,
            holes=[(self.LS[int(len(self.LS)/2)].x[0]*1.05, self.LS[int(len(self.LS)/2)].y[0])]
            
        )
        tri = tr.triangulate(tri, 'p')

        # Plotting
        ax.tricontourf(y, x, tri['triangles'], T, levels, cmap=colormap)
        cntr = ax.tricontourf([-i for i in y], x, tri['triangles'], T, levels, cmap=colormap)
        
        cbar = plt.gcf().colorbar(cntr, ax=ax)
        cbar.set_label('Temperature [K]')
        ax.set_title(f'{surface.capitalize()} Surface Temperature Contour')
        ax.axis('equal')
        ax.set_ylabel('x [m]')
        ax.set_xlabel('y [m]')
        ax.grid()

        if return_flag:
            return ax

    def plot3d(self, plotOffline=False):
        xls = self.getAllPlaneDataInLine('lower', 'x')
        yls = self.getAllPlaneDataInLine('lower', 'y')
        zls = self.getAllPlaneDataInLine('lower', 'z')
        xus = self.getAllPlaneDataInLine('upper', 'x')
        yus = self.getAllPlaneDataInLine('upper', 'y')
        zus = self.getAllPlaneDataInLine('upper', 'z')

        # Indeces of Waverider boundaries
        seg = self._boundaryIndeces

        points_ls = np.stack((xls, yls), axis=1)
        A_ls = dict(
            vertices=points_ls,
            segments=seg,
            holes=[(self.LS[int(len(self.LS)/2)].x[0]*1.05, self.LS[int(len(self.LS)/2)].y[0])]
        )
        tri_ls = tr.triangulate(A_ls, 'p')
        points_us = np.stack((xus, yus), axis=1)
        A_us = dict(
            vertices=points_us,
            segments=seg,
            holes=[(self.LS[int(len(self.LS)/2)].x[0]*1.05, self.LS[int(len(self.LS)/2)].y[0])]
        )
        tri_us = tr.triangulate(A_us, 'p')

        data0 = go.Mesh3d(x=xls, y=yls, z=-np.array(zls),
                          i=[i[0] for i in tri_ls['triangles']],
                          j=[i[1] for i in tri_ls['triangles']],
                          k=[i[2] for i in tri_ls['triangles']],
                          color='red')
        data1 = go.Mesh3d(x=xus, y=yus, z=-np.array(zus),
                          i=[i[0] for i in tri_us['triangles']],
                          j=[i[1] for i in tri_us['triangles']],
                          k=[i[2] for i in tri_us['triangles']],
                          color='blue')
        data2 = go.Mesh3d(x=xls, y=-np.array(yls), z=-np.array(zls),
                          i=[i[0] for i in tri_ls['triangles']],
                          j=[i[1] for i in tri_ls['triangles']],
                          k=[i[2] for i in tri_ls['triangles']],
                          color='red')
        data3 = go.Mesh3d(x=xus, y=-np.array(yus), z=-np.array(zus),
                          i=[i[0] for i in tri_us['triangles']],
                          j=[i[1] for i in tri_us['triangles']],
                          k=[i[2] for i in tri_us['triangles']],
                          color='blue')
        
        if plotOffline:
            # Offline Plot (creating HTML file)
            fig = plotly.offline.plot(
                {"data": [data0, data1, data2, data3], 
                "layout": go.Layout(scene=dict(yaxis = dict(range=[-self.l/2, self.l/2]),
                                    zaxis = dict(range=[-self.l/2, self.l/2]),
                                    xaxis_visible=False, yaxis_visible=False, zaxis_visible=False),
                                    scene_aspectmode='cube')},
                filename= cfg.CONFIG['Viscous Method'] + ".html")
        else:
            # Online Plot
            fig = go.Figure(data=[data0, data1, data2, data3])
            fig.update_layout(scene=dict(yaxis = dict(range=[-self.l/2, self.l/2]),
                                        zaxis = dict(range=[-self.l/2, self.l/2]),
                                        xaxis_visible=False, yaxis_visible=False, zaxis_visible=False),
                                        scene_aspectmode='cube')
            fig.show()

    def plotRadius(self, ax=None):
        return_flag = False
        if ax is None:
            _, ax = plt.subplots()
            return_flag = True
        
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

        if return_flag:
            return ax
           
    def plotBase(self, ax=None):
        return_flag = False
        if ax is None:
            _, ax = plt.subplots()
            return_flag = True

        self.LE.plot(ax, False ,'b', label='Upper Surface')
        self.SW[0].plot(ax, False, '--k', label='Shockwave')
        self.SW[1].plot(ax, False ,'--k', label=None)

        self.LE.plot(ax, True ,'b', label=None)
        self.SW[0].plot(ax, True, '--k', label=None)
        self.SW[1].plot(ax, True ,'--k', label=None)

        ax.plot([i.y[-1] for i in self.LS],
                [i.z[-1] for i in self.LS],
                'r', label='Lower Surface')
        ax.plot([-i.y[-1] for i in self.LS],
                [i.z[-1] for i in self.LS ],
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

        if return_flag:
            return ax

    def plotTop(self, ax=None):
        return_flag = False
        if ax is None:
            _, ax = plt.subplots()
            return_flag = True
        
        y = [i.y[0] for i in self.LS]
        y_ = [-i.y[0] for i in self.LS]
        x = [i.x[0] for i in self.LS]
        # m = self.l*1.02

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

        if return_flag:
            return ax

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
    from config_modify import results_load_json, modify_config

    results = results_load_json()

    test_wr = Waverider(results['WR'])

    # test_wr.plot3d()

    test_wr.plotPressure(plot3d=True)
    test_wr.plotStress(plot3d=True)
    test_wr.plotTemperature(plot3d=True)
    test_wr.plotTemperature(plot3d=True, ref=True)

    # fig_temp, ax = plt.subplots()
    # test_wr.plotTemperature(ax=ax, surface='lower')
    # fig_temp_ref, ax = plt.subplots()
    # test_wr.plotTemperature(ax=ax, surface='lower', ref=True)

    # fig_press, ax = plt.subplots()
    # test_wr.plotPressure('lower', ax=ax)
    # fig_stress_u, ax = plt.subplots()
    # test_wr.plotStress('upper', ax=ax)
    # fig_stress_l, ax = plt.subplots()
    # test_wr.plotStress('lower', ax=ax)

    # test_wr.plotAll()
    # from pathlib import Path
    # savepath = Path(__file__).parents[1].absolute()
    # savepath /= 'results\\plots'
    # fig_temp.savefig(
    #     savepath / f'temperature_contour_{cfg.CONFIG["Viscous Method"]}.svg',
    #     format='svg'
    # )
    # fig_temp_ref.savefig(
    #     savepath / f'ref_temperature_contour_{cfg.CONFIG["Viscous Method"]}.svg',
    #     format='svg'
    # )
    # fig_press.savefig(
    #     savepath / f'pressure_contour_{cfg.CONFIG["Viscous Method"]}.svg',
    #     format='svg'
    # )
    # fig_stress_u.savefig(
    #     savepath / f'stress_upper_contour_{cfg.CONFIG["Viscous Method"]}.svg',
    #     format='svg'
    # )
    # fig_stress_l.savefig(
    #     savepath / f'stress_lower_contour_{cfg.CONFIG["Viscous Method"]}.svg',
    #     format='svg'
    # )

    # plt.show()
    