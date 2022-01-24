import numpy as np
import matplotlib.pyplot as plt
import config as cfg
from typing import Tuple


class TaylorMaccoll:
    """
    Taylor-Maccoll solver
    Calculation of flow properties using the Taylor-Maccoll equation (Non-dimensional)
    Inputs:
        b: Oblique shockwave angle [deg]
        Minf: Free flow Mach number

    Properties:
        P: Pressure ratio (P/Pinf)
        V: Velocity ratio (V/Vmax, Vmax=(2*h0)^0.5)
        Vr: Radial velocity ratio (Vr/Vmax)
        Vtheta: Tangental velocity ratio (Vtheta/Vmax)
        M: Local Mach number
        T: Temperature ratio (T/Tinf)
        theta: Angles in inside the shochwave. Each angle has a corresponding set of properties.
               (limits: theta_cone <= theta <= b)
    """
    def __init__(self, Minf, b) -> None:
        def TM_eq(x, y, z) -> Tuple[float, float]:
            """
            Function for Taylor-Maccoll Equation
            Inputs:
               x: theta angle
               y = Vr
               z = Vtheta

            Output:
               dydx = dVr/dtheta = Vtheta
               dzdx = dVtheta/dtheta = d^2Vr/dtheta^2
               """
            dydx = z
            dzdx = (z**2*y - (cfg.GAM - 1)/2*(1 - y**2 - z**2)*(2*y + z*cot(x))) / \
                ((cfg.GAM - 1)/2*(1 - y**2 - z**2) - z**2)
            return dydx, dzdx

        def cot(rad) -> float:
            """Funtion for the cotangent"""
            return 1/(np.tan(rad))

        # Constants
        Dtheta = np.deg2rad(-1e-2) # Angle Step [rad]

        b = b*np.pi/180

        # Pre-allocation
        theta = np.arange(b, 0 + Dtheta, Dtheta)
        theta[-1] = 0.0
        N = len(theta)
        Vtheta = np.empty(N)
        Vr = np.empty(N)

        i = 0

        if b <= np.arcsin(1/Minf):
            Run = False
            d = 0.
            theta[i] = 0.
            Vr[i] = 0.
            Vtheta[i] = 0.
            P2P1 = 0.
            Po2P2 = 0.
            T2T1 = 0.
            To2T2 = 0.

        else:
            Run = True
            Mn1 = Minf*np.sin(b) # Normal Vector of Mach Before
            Mn2 = np.sqrt((Mn1**2 + (2/(cfg.GAM - 1))) / \
                (2*cfg.GAM/(cfg.GAM - 1)*Mn1**2 - 1)) # Normal Vector of Mach After
            d = np.arctan(2*cot(b)*((Minf*np.sin(b))**2 - 1) / \
                (Minf**2*(cfg.GAM + np.cos(2*b)) + 2)) # Deflection Angle [rad]
            M2 = Mn2/(np.sin(b - d)) # Mach Number After Shockwave
            V2 = (1 + 2/((cfg.GAM -1)*M2**2))**(-0.5) # Velocity After Shockwave
            P2P1 = 1 + (2*cfg.GAM)/(cfg.GAM + 1)*(Mn1**2 - 1) # Shockwave Relation
            Po2P2 = (1 + ((cfg.GAM - 1)/2)*M2**2)**(cfg.GAM/(cfg.GAM - 1)) # Isentropic Relation
            T2T1 = (1 + (2*cfg.GAM/(cfg.GAM + 1))*(Mn1**2 - 1))*((2 + (cfg.GAM - 1)*Mn1**2) / \
                ((cfg.GAM + 1)*Mn1**2)) # Shockwave Relation
            To2T2 = 1 + ((cfg.GAM - 1)/2)*M2**2 # Isentropic Relation
            theta[i] = b
            Vr[i] = V2*np.cos(b - d)
            Vtheta[i] = -V2*np.sin(b - d)

        # Runge-Kutta
        while Vtheta[i] < 0 and theta[i] > 0 and Run:
            x0 = theta[i]
            y0 = Vr[i]
            z0 = Vtheta[i]

            k1, l1 = TM_eq(x0, y0, z0)
            k1 *= Dtheta
            l1 *= Dtheta

            k2, l2 = TM_eq(x0 + Dtheta/2, y0 + k1/2, z0 + l1/2)
            k2 *= Dtheta
            l2 *= Dtheta

            k3, l3 = TM_eq(x0 + Dtheta/2, y0 + k2/2, z0 + l2/2)
            k3 *= Dtheta
            l3 *= Dtheta

            k4, l4 = TM_eq(x0 + Dtheta, y0 + k3, z0 + l3)
            k4 *= Dtheta
            l4 *= Dtheta

            i += 1
            Vr[i] = Vr[i-1] + (k1 + 2.*k2 + 2.*k3 + k4)/6
            Vtheta[i] = Vtheta[i-1] + (l1 + 2.*l2 + 2.*l3 + l4)/6

        # Interpolation for last angle
        theta[i] = np.interp(0, [Vtheta[i-1], Vtheta[i]], [theta[i-1], theta[i]])
        Vr[i] = np.interp(0, [Vtheta[i-1], Vtheta[i]], [Vr[i-1], Vr[i]])
        Vtheta[i] = 0.0

        self.theta = np.flip(theta[:i+1]*180/np.pi)
        self.Vr = np.flip(Vr[:i+1])
        self.Vtheta = np.flip(Vtheta[:i+1])
        self.V = np.sqrt(self.Vr**2 + self.Vtheta**2)
        self.M = np.sqrt(2./((cfg.GAM - 1)*(self.V**(-2) - 1)))
        self.P = P2P1*Po2P2*(1 + ((cfg.GAM - 1)/2)*self.M**2)**(-cfg.GAM/(cfg.GAM - 1))
        self.T = T2T1*To2T2*(1 + ((cfg.GAM - 1)/2)*self.M**2)**(-1)


class Edge:
    """
    Object for the edge of the waverider (the last plane)
    Properties:
        b: Oblique shockwave angle [rad]
        x: x coordinate of each point in the crossection's reference system (at the end of the shockwave) [m]
        s: z coordinate of each point in the crossection's reference system (at the end of the shockwave) [m]
        P: Non-dimensional local pressure for each point (P/Pinf)
        T: Non-dimensional local temperature for each point (T/Tinf)
        M: Local Mach number for each point
        N: Number of points
    """
    b: float
    x: np.ndarray
    s: np.ndarray
    P: np.ndarray
    T: np.ndarray
    M: np.ndarray
    N: int = 1

    def __init__(self, Taylor, ysw, zsw) -> None:
        self.x = np.array(0.)
        self.s = np.array(0.)
        self.y = np.array(ysw)
        self.z = np.array(zsw)
        self.P = np.array(Taylor.P[0])
        self.T = np.array(Taylor.T[0])
        self.M = np.array(Taylor.M[0])


class _Streamline:
    """
    USED ONLY GOT INHERETENCE REASONS!
    Internal object that represents a streamline.
    Properties:
        b: Oblique shockwave angle [rad]
        x_ref: x coordinate of each point in the crossection's frame of reference (at the end of the shockwave) [m]
        s: z coordinate of each point in the crossection's frame of reference (at the end of the shockwave) [m]
        y: y coordinate of each point (waverider frame of reference) [m]
        z: z coordinate of each point (waverider frame of reference) [m]
        P: Non-dimensional local pressure for each point (P/Pinf)
        T: Non-dimensional local temperature for each point (T/Tinf)
        M: Local Mach number for each point
        N: Number of points
    """
    b: float
    x: np.ndarray
    s: np.ndarray
    y: np.ndarray
    z: np.ndarray
    P: np.ndarray
    T: np.ndarray
    M: np.ndarray
    N: int

    def coor_change(self, ysw, zsw, phi=np.pi/2) -> None:
        """
        Method for change the coordinate system to that of the waverider
        """
        self.y = ysw + self.s * np.cos(phi)
        self.z = zsw + self.s * np.sin(phi)

    def stream_plot(self) -> None:
        """
        Method for plotting the streamline
        """
        x = self.x[0] - self.x
        s = self.s
        plt.plot(x, s, c='blue', label='Streamline')


class Wedge(_Streamline):
    """
    Object that represents a streamline in supersonic wedge flow
    Inputs:
        Minf: Free flow Mach number
        b: Oblique shockwave angle [deg]
        dle: distance from leading edge on the z-axis in the crossection's frame of reference [m]
        N: number of points
    """
    def __init__(self, Minf, b, dle, N) -> None:
        def cot(rad) -> float:
            """Funtion for the cotangent"""
            return 1/(np.tan(rad))

        # deg to rad
        b = b * np.pi / 180

        # Deflection Angle
        tand = 2 * cot(b) * ((Minf * np.sin(b)) ** 2 - 1) / \
            (Minf ** 2 * (cfg.GAM + np.cos(2 * b)) + 2) 
        d = np.arctan(tand)

        # Normal Mach 1
        Mn1 = Minf * np.sin(b)
        
        # Normal Mach 2
        Mn2 = np.sqrt((Mn1 ** 2 + (2 / (cfg.GAM - 1))) /
            (2 * cfg.GAM / (cfg.GAM - 1) * Mn1 ** 2 - 1))

        # Mach 2
        M2 = Mn2 / (np.sin(b - d))

        # From Shochwave Relations
        # Pressure ratio (P2/P1)
        P2P1 = 1 + (2 * cfg.GAM) / (cfg.GAM + 1) * (Mn1 ** 2 - 1)
        # Temperature ratio (T2/T1)
        T2T1 = (1 + (2 * cfg.GAM / (cfg.GAM + 1)) * (Mn1 ** 2 - 1)) * \
            ((2 + (cfg.GAM - 1) * Mn1 ** 2) / ((cfg.GAM + 1) * Mn1 ** 2))

        # Length of crossection
        L = dle/np.tan(b)

        # Set properties
        self.b = b*180/np.pi
        self.x = np.linspace(L, 0, N)
        self.s = np.linspace(-dle, L*tand-dle, N)
        self.P = P2P1*np.ones(N)*cfg.ATM.P
        self.T = T2T1*np.ones(N)*cfg.ATM.T
        self.M = M2*np.ones(N)
        self.N = N


class Cone(_Streamline):
    """
    Object that represents a streamline in supersonic conical flow
    Inputs:
        b: Shockwave cone angle [deg]
        R: Base radius of Mach cone [m]
        dle: distance from leading edge on the z-axis in the crossection's frame of reference [m]
        Taylor: From Taylor-Maccoll (TaylorMaccoll class)
        Dt: Time step [s]

    Extra properties:
        R: Base radius of Mach cone [m]
        theta_c: Half cone angle [rad]
    """
    def __init__(self, Minf, b, R, dle, Taylor, Dt) -> None:
        # Runge-Kutta Equations
        # In Vector Form: ds/dt = V
        def drdt(x) -> float:
            """
            dr/dt = Vr
            x is in rads
            """
            return np.interp(x*180./np.pi, Taylor.theta, Vr)

        def dthetadt(x, r) -> float:
            """
            dtheta/dt = Vtheta/r
            x is in rads
            """
            return np.interp(x*180./np.pi, Taylor.theta, Vtheta) / r

        b = b*np.pi/180.
        z0 = R - dle
        L = R / np.tan(b)
        h0 = cfg.CP*cfg.ATM.T + ((Minf*cfg.ATM.v_sonic)**2)/2
        Vmax = np.sqrt(2*h0)
        Vr = Taylor.Vr*Vmax
        Vtheta = Taylor.Vtheta*Vmax

        # Runge-Kutte
        i = 0
        theta_s = np.array([b])
        r_s = np.array([z0/np.sin(b)])

        while r_s[i]*np.cos(theta_s[i]) <= L:
            if i >= int(1.5*cfg.N_POINTS):
                raise Exception("Too many iterations")

            k1 = Dt*drdt(theta_s[i])
            l1 = Dt*dthetadt(theta_s[i], r_s[i])

            k2 = Dt*drdt(theta_s[i] + l1/2.)
            l2 = Dt*dthetadt(theta_s[i] + l1/2., r_s[i] + k1/2.)

            k3 = Dt*drdt(theta_s[i] + l2/2.)
            l3 = Dt*dthetadt(theta_s[i] + l2/2., r_s[i] + k2/2.)

            k4 = Dt*drdt(theta_s[i] + l3)
            l4 = Dt*dthetadt(theta_s[i] + l3, r_s[i] + k3)

            i += 1
            r_s = np.append(r_s, r_s[i-1] + (k1 + 2.*k2 + 2.*k3 + k4)/6.)
            theta_s = np.append(theta_s, theta_s[i-1] + (l1 + 2.*l2 + 2.*l3 + l4)/6.)

        x = r_s*np.cos(theta_s)
        z = r_s*np.sin(theta_s)
        P = np.interp(theta_s*180/np.pi, Taylor.theta, Taylor.P)
        T = np.interp(theta_s*180/np.pi, Taylor.theta, Taylor.T)
        M = np.interp(theta_s*180/np.pi, Taylor.theta, Taylor.M)
        V = np.interp(theta_s*180/np.pi, Taylor.theta, Taylor.V)

        # Sistima sintetagmenon sti vasi konou
        x = L - x
        z = z - R

        # Last Point interpolation
        P[-1] = np.interp(0., [x[-1], x[-2]], [P[-1], P[-2]])
        T[-1] = np.interp(0., [x[-1], x[-2]], [T[-1], T[-2]])
        M[-1] = np.interp(0., [x[-1], x[-2]], [M[-1], M[-2]])
        V[-1] = np.interp(0., [x[-1], x[-2]], [V[-1], V[-2]])
        z[-1] = np.interp(0., [x[-1], x[-2]], [z[-1], z[-2]])
        x[-1] = 0.

        # Set Properties
        # Cone
        self.b = b
        self.R = R
        self.theta_c = Taylor.theta[0]
        
        # Streamline
        self.x = x
        self.s = z
        self.M = M
        self.V = V * Vmax
        self.P = P * cfg.ATM.P
        self.T = T * cfg.ATM.T
        self.N = len(self.x)


class Upper:
    """Object representing the upper surface of each crossection"""
    def __init__(self, lower_obj) -> None:
        self.N = lower_obj.N
        self.x = lower_obj.x

        if lower_obj.y.ndim == 0:
            self.y = lower_obj.y * np.ones(self.N)
            self.z = lower_obj.z * np.ones(self.N)

        else:
            self.y = lower_obj.y[0] * np.ones(self.N)
            self.z = lower_obj.z[0] * np.ones(self.N)

        self.P = cfg.ATM.P * np.ones(self.N)
        self.T = cfg.ATM.T * np.ones(self.N)
        self.M = cfg.MINF * np.ones(self.N)


class Base:
    """Object representing the base surface of each crossection"""
    def __init__(self, lower_obj):
        self.x = np.zeros(2)

        if lower_obj.y.ndim == 0:
            self.y = np.array([lower_obj.y, lower_obj.y])
            self.z = np.array([lower_obj.z, lower_obj.z])
        
        else:
            self.y = np.array([lower_obj.y[0], lower_obj.y[-1]])
            self.z = np.array([lower_obj.z[0], lower_obj.z[-1]])


if __name__ == '__main__':
    def testing():
        import config as cfg
        from configModify import configCreateWindow
        configCreateWindow()

        b = 13.
        dint = 1.2
        T = TaylorMaccoll(cfg.MINF, b)

        R = 2
        Dt = (dint/np.tan(b*np.pi/180)) / (T.V[-1]*np.sqrt(2*1000*cfg.ATM.T + ((cfg.MINF*cfg.ATM.v_sonic)**2))) / cfg.N_POINTS
        LS = Cone(cfg.MINF, b, R, dint, T, Dt)
        LS.coor_change(0, 0, 0)

        BS = Base(LS)
        US = Upper(LS)
        
        L = R / np.tan(b*np.pi/180)
        x = L - LS.x
        z = LS.s + R
        d = T.theta[0]

        plt.plot(x, z, 'b')
        plt.plot([0, L], [0, R], 'k')
        plt.plot([0, L], [0, L*np.tan(d*np.pi/180)], 'r')
        plt.gca().invert_yaxis()
        plt.gca().axis('equal')
        plt.grid()
        plt.show()

    testing()
