"""
The module responsible for the evaluation of the aerodynamic perfomance of a Waverider.
"""
import numpy as np
import scipy.optimize as froot
import config as cfg
from typing import Tuple


def cf_solver(re_x, mach_d) -> float:
    """
    The solver of the coefficient of friction for given Reynold and Mach numbers
    Inputs:
        re_x: Reynolds at location x
        mach_d: Mach number at the edge of bountry layer
    """
    def cf_fun(c_f) -> float:
        """
        The equation for the coefficient of friction. For solution left must equal right (return 0).
        """
        left = (0.242/np.sqrt(0.2*c_f*mach_d**2))*(np.arcsin(a) + np.arcsin(a))
        right = 0.41 + np.log10(re_x*c_f)
        return left - right

    A = 0.2*mach_d**2
    a = A/np.sqrt(A**2 + 4*A)

    cfi = 0.46/(np.log10(re_x)**2.6)
    c_f = froot.brentq(cf_fun, 0.01*cfi, cfi, xtol=1e-5, rtol=1e-7)
    return c_f


def van_driest_method(l_obj, x) -> Tuple[np.ndarray, np.ndarray]:
    """
    Function that returns the wall shear stress of the plane using the Van Driest method.
    Inputs:
        l_obj: Object representing the lower surface of a crossection
               of the waverider at a given plane (instance of _Streamline)
        x: numpy array or list of x position of each point (from the leading edge)
    """
    tw_ls = np.zeros(l_obj.N)
    tw_us = np.zeros(l_obj.N)

    rhod = l_obj.P / (cfg.R_GAS * l_obj.T / cfg.MB_AIR)
    for i in range(l_obj.N):
        #Lower Surface
        ad = np.sqrt(cfg.GAM * cfg.R_GAS * l_obj.T[i] / cfg.MB_AIR)
        mid = (1.458e-6 * l_obj.T[i] ** 1.5) / (l_obj.T[i] + 110.4)
        re_x = (rhod[i] * l_obj.M[i] * ad*x[i] / mid)
        if re_x != 0:
            c_f = cf_solver(re_x, l_obj.M[i])
        else:
            c_f = 0.
        tw_ls[i] = 0.5*rhod[i]*c_f*(l_obj.M[i]*ad) ** 2

        #Upper Surface
        re_x = (cfg.ATM.rho * cfg.MINF * cfg.ATM.v_sonic * x[i] / cfg.ATM.mu)
        if re_x != 0:
            c_f = cf_solver(re_x, cfg.MINF)
        else:
            c_f = 0.
        tw_us[i] = 0.5 * cfg.ATM.rho * c_f * (cfg.MINF * cfg.ATM.v_sonic) ** 2

    return tw_ls, tw_us


def ref_temp_method(l_obj, x) -> Tuple[np.ndarray, np.ndarray]:
    """
    Function that returns the wall shear stress of the plane using the Reference Temperature method.
    Inputs:
        l_obj: Object representing the lower surface of a crossection
               of the waverider at a given plane (instance of _Streamline)
        x: numpy array or list of x position of each point (from the leading edge)
    """
    tw_ls = np.zeros(l_obj.N)
    tw_us = np.zeros(l_obj.N)
    r = 0.88 #Recovery Factor for Turbulant Flow, Pr_T=0.86 & Pr_L=0.715
    
    T0_us = cfg.ATM.T*(1. + ((cfg.GAM - 1)/2)*cfg.MINF**2)
    Tw_us = cfg.ATM.T + r*(T0_us - cfg.ATM.T)
    Tref_us = cfg.ATM.T*(1. + 0.032*cfg.MINF**2 + 0.58*(Tw_us/cfg.ATM.T - 1.))
    miref_us = (1.458e-6*Tref_us**1.5)/(Tref_us + 110.4)
    rhoref_us = cfg.ATM.P/(cfg.R_GAS*Tref_us/cfg.MB_AIR)

    for i in range(l_obj.N):
        # Lower Surface
        ad = np.sqrt(cfg.GAM*cfg.R_GAS*l_obj.T[i]/cfg.MB_AIR)
        T0 = l_obj.T[i]*(1. + ((cfg.GAM - 1)/2)*l_obj.M[i]**2)
        Tw = l_obj.T[i] + r*(T0 - l_obj.T[i])
        Tref = l_obj.T[i]*(1. + 0.032*l_obj.M[i]**2 + 0.58*(Tw/l_obj.T[i] - 1.))
        miref = (1.458e-6*Tref**1.5)/(Tref + 110.4)
        rhoref = l_obj.P[i]/(cfg.R_GAS*Tref/cfg.MB_AIR)
        re_x = (rhoref*l_obj.M[i]*ad*x[i]/miref)
        if re_x != 0:
            c_f = 0.0592/(re_x**0.2)
        else:
            c_f = 0.
        tw_ls[i] = 0.5*rhoref*c_f*(l_obj.M[i]*ad)**2

        # Upper Surface
        re_x = (rhoref_us*cfg.MINF*cfg.ATM.v_sonic*x[i]/miref_us)
        if re_x != 0:
            c_f = 0.0592/(re_x**0.2)
        else:
            c_f = 0.
        tw_us[i] = 0.5*rhoref_us*c_f*(cfg.MINF*cfg.ATM.v_sonic)**2

    return tw_ls, tw_us


def plane_area(l_obj) -> float:
    """
    Function that return the area of the crossection
    Inputs:
        l_obj: Object representing the lower surface of a crossection
               of the waverider at a given plane (instance of _Streamline)
    """
    if l_obj.N == 1:
        A = 0
    else:
        A = np.trapz(np.array([l_obj.s[0], l_obj.s[0]]), x=np.array([l_obj.x[0], l_obj.x[-1]])) - \
            np.trapz(l_obj.s, x=l_obj.x)
    return A


def base_surf_press(l_obj_list) -> Tuple[float, float]:
    """
    Function that return the base pressure and surface of the waverider
    Inputs:
        l_obj_list: List of objects representing the lower surface of the crossection
                    for each plane (instance of _Streamline)
    """
    n = len(l_obj_list)
    yls = np.empty(n)
    zls = np.empty(n)
    yus = np.empty(n)
    zus = np.empty(n)

    cp = (2./(cfg.GAM*cfg.MINF**2))*(((2./(cfg.GAM+1))**1.4)* \
        ((1/cfg.MINF)**2.8)*((2*cfg.GAM*cfg.MINF**2 - (cfg.GAM-1))/(cfg.GAM+1)) - 1)
    P = cfg.ATM.P + 0.5*cfg.ATM.rho*cp*(cfg.MINF*cfg.ATM.v_sonic)**2
     
    for i, plane in enumerate(l_obj_list):
        if plane.N == 1:
            yls[i] = plane.y
            zls[i] = plane.z
            yus[i] = plane.y
            zus[i] = plane.z
        else:
            yls[i] = plane.y[-1]
            zls[i] = plane.z[-1]
            yus[i] = plane.y[0]
            zus[i] = plane.z[0]
    S = np.trapz(zls, x=yls) - np.trapz(zus, x=yus)
    return P, S


def base_drag(l_obj_list) -> float:
    """
    Function that returns the base drag
    Inputs:
        l_obj_list: List of objects representing the lower surface of the crossection
                    for each plane (instance of _Streamline)
    """
    P, S = base_surf_press(l_obj_list)
    return P * S


def pressure_forces(l_obj) -> Tuple[float, float]:
    """
    Function that returns the 2D pressure forces for a given crossection (Lift, Drag) [N/m]
    Inputs:
        l_obj: Object representing the lower surface of a crossection
               of the waverider at a given plane (instance of _Streamline)
    """
    N = l_obj.N
    length = np.abs(l_obj.x[0] - l_obj.x[-1])

    Fx = np.empty(N-1)
    Fz = np.empty(N-1)

    for i in range(N-1):
        ds = np.sqrt((l_obj.x[i+1] - l_obj.x[i])**2 + \
            (l_obj.s[i+1] - l_obj.s[i])**2)
        delt = np.abs(np.arctan((l_obj.s[i+1] - l_obj.s[i])/ \
            (l_obj.x[i+1] - l_obj.x[i])))
        F = (l_obj.P[i+1] + l_obj.P[i])*ds/2.
        Fx[i] = F*np.sin(delt)
        Fz[i] = F*np.cos(delt)

    D = np.sum(Fx)
    L = np.sum(Fz) - cfg.ATM.P*length
    return L, D


def viscous_forces(l_obj) -> Tuple[float, float]:
    """
    Function that returns the 2D viscous forces for a given crossection (Lift, Drag) [N/m]
    Inputs:
        l_obj: Object representing the lower surface of a crossection
               of the waverider at a given plane (instance of _Streamline)
    """
    length = np.abs(l_obj.x[0] - l_obj.x[-1])
    x = length*np.ones(l_obj.N) - l_obj.x

    # Pre-allocation
    Fx = np.zeros(l_obj.N-1)
    Fz = np.zeros(l_obj.N-1)

    # Stress Calculation
    if cfg.METHOD == 1: #Van Driest Method
        tw_ls, tw_us = van_driest_method(l_obj, x)
    
    elif cfg.METHOD == 2: #Reference Tempurature
        tw_ls, tw_us = ref_temp_method(l_obj, x)

    else:
        raise Exception('Unknown Viscous Method (Cf Calulation)')
    
    # Forces Calculation
    for i in range(l_obj.N-1):
        ds = np.sqrt((x[i+1] - x[i])**2 + \
            (l_obj.s[i+1] - l_obj.s[i])**2)
        delt = np.abs(np.arctan((l_obj.s[i+1] - l_obj.s[i])/ \
            (x[i+1] - x[i])))
        F_ls = (tw_ls[i+1] + tw_ls[i])*ds/2.
        F_us = (tw_us[i+1] + tw_us[i])*ds/2.
        Fx[i] = F_ls*np.cos(delt) + F_us
        Fz[i] = F_ls*np.sin(delt)

    D = np.sum(Fx)
    L = np.sum(Fz)
    return L, D


class AeroPlane:
    """
    Class that calculates the Invicid 2D forces and the area of a plane
    Inputs:
      l_obj: Object representing the lower surface of a crossection
               of the waverider at a given plane (instance of _Streamline)

    Proberties:
        L: Lift [N/m]
        D: Drag [N/m]
        A: Area [m^2]
    """
    L: float 
    D: float
    A: float

    def __init__(self, l_obj) -> None:
        if l_obj.N == 1:
            self.A = 0.
            self.L = 0.
            self.D = 0.

        else:
            # Area of Plane
            self.A = plane_area(l_obj)

            # Pressure forces
            L_p, D_p = pressure_forces(l_obj)

            # Viscous forces
            if cfg.VISCOUS:
                L_v, D_v = viscous_forces(l_obj)
            else:
                L_v = 0.
                D_v = 0.

            self.L = L_p - L_v
            self.D = D_p + D_v


class Aero3d:
    """
    Class that calculates the Invicid 3D forces, volume and reference surface of a waverider
    Inputs:
      l_obj_list: List of bbjects representing the lower surface of the crossection
                  for each plane (subclass of _Streamline)
    
    Properties
        L: Lift [N]
        D: Drag [N]
        V: Volume [m^3]
        S: Reference Surface [m^2]
    """
    L: float 
    D: float
    V: float
    S: float

    def __init__(self, l_obj_list, phi):
        n = len(l_obj_list)
        L2d = np.empty(n)
        D2d = np.empty(n)
        A2d = np.empty(n)
        ym = np.empty(n)
        ycg = np.empty(n)
        yle = np.empty(n)
        xle = np.empty(n)

        for i in range(n):
            if i == n-1:
                yle[i] = l_obj_list[i].y
                yte = l_obj_list[i].y
                xle[i] = l_obj_list[i].x
            else:
                yle[i] = l_obj_list[i].y[0]
                yte = l_obj_list[i].y[-1]
                xle[i] = l_obj_list[i].x[0]
            
            aero2d = AeroPlane(l_obj_list[i])
            ym[i] = yle[i] + np.abs(yte - yle[i]) / 3
            ycg[i] = yle[i] + np.abs(yte - yle[i]) / 3
            L2d[i] = aero2d.L * np.sin(phi[i])
            D2d[i] = aero2d.D
            A2d[i] = aero2d.A

        self.L = 2*np.trapz(L2d, x=ycg)
        self.D = 2*np.trapz(D2d, x=ycg) - base_drag(l_obj_list)
        self.V = 2*np.trapz(A2d, x=ym)
        self.S = 2*np.trapz(xle, x=yle)
        
    def __str__(self) -> str:
        return (
            f'Lift = {self.L * 1e-3} kN\n'
            f'Drag = {self.D * 1e-3} kN\n'
            f'Volume = {self.V} m3\n'
            f'Reference Surface = {self.S} m2'
            )


if __name__ == "__main__":
    def testing() -> None:
        pass

    testing()
