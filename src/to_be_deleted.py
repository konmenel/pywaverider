"""Runge-Kutta for Taylor-Maccoll using scipy"""
import time
from typing import Tuple
import numpy as np
from scipy.integrate import solve_ivp
import config as cfg
from config_modify import config_create_window
from streamline import TaylorMaccoll


# scipy Runge-Kutta faster for Taylor-Maccoll, slower for streamline. WHY??????????
# ANSWER: use t_eval argument instead of max_step make the class much faster. Now both
#       are faster than my implementations.
#
# Implemented on the original. Speed improved a lot



def my_runge_kutta(Minf, b) -> None:
    def TM_eq(x, y, z) -> Tuple[float, float]:
        """Function for Taylor-Maccoll Equation
        Inputs:
            x: theta angle
            y = Vr
            z = Vtheta

        Output:
            dydx = dVr/dtheta = Vtheta
            dzdx = dVtheta/dtheta = d^2Vr/dtheta^2
            """
        dydx = z
        dzdx = (z ** 2 * y - 0.2 * (1 - y ** 2 - z ** 2) * (2 * y + z * cot(x))) / \
            (0.2 * (1 - y ** 2 - z ** 2) - z ** 2)
        return dydx, dzdx

    def cot(rad) -> float:
        """Funtion for the cotangent"""
        return 1/(np.tan(rad))

    # Constants
    b = b * np.pi/180

    # Pre-allocation
    theta = np.arange(b, 0 + Dtheta, Dtheta)
    theta[-1] = 0.0
    N = len(theta)
    Vtheta = np.empty(N)
    Vr = np.empty(N)

    i = 0

    
    Run = True
    Mn1 = Minf*np.sin(b) # Normal Vector of Mach Before
    Mn2 = np.sqrt((Mn1 ** 2 + 2 / 0.4) / \
        (2.8 / 0.4 * Mn1 ** 2 - 1)) # Normal Vector of Mach After
    d = np.arctan(2 * cot(b) * ((Minf * np.sin(b)) ** 2 - 1) / \
        (Minf ** 2 * (1.4 + np.cos(2*b)) + 2)) # Deflection Angle [rad]
    M2 = Mn2 / (np.sin(b - d)) # Mach Number After Shockwave
    V2 = (1 + 2 / (0.4 * M2 ** 2)) ** (-0.5) # Velocity After Shockwave
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
    Vtheta[i] = 0.

    return theta[:i], Vr[:i], Vtheta[:i]


def other_runge_kutta(Minf, b) -> None:
    def TM_eq(x, y) -> Tuple[float, float]:
        """Function for Taylor-Maccoll Equation
        Inputs:
            x = theta angle
            y[0] = y = Vr
            y[1] = z = Vtheta

        Output:
            dydx = dVr/dtheta = Vtheta
            dzdx = dVtheta/dtheta = d^2Vr/dtheta^2
            """
        dydx = y[1]
        dzdx = (y[1] ** 2 * y[0] - 0.2 * (1 - y[0] ** 2 - y[1] ** 2) * (2 * y[0] + y[1] * cot(x))) / \
            (0.2 * (1 - y[0] ** 2 - y[1] ** 2) - y[1] ** 2)
        return dydx, dzdx
    
    def cone_surf(x, y):
        return y[1]

    def cot(rad) -> float:
        """Funtion for the cotangent"""
        return 1/(np.tan(rad))

    # Constants
    b = b * np.pi/180

    # Pre-allocation
    theta = np.arange(b, 0 + Dtheta, Dtheta)
    theta[-1] = 0.0
    N = len(theta)
    Vtheta = np.empty(N)
    Vr = np.empty(N)

    i = 0

    Mn1 = Minf*np.sin(b) # Normal Vector of Mach Before
    Mn2 = np.sqrt((Mn1 ** 2 + 2 / 0.4) / \
        (2.8 / 0.4 * Mn1 ** 2 - 1)) # Normal Vector of Mach After
    d = np.arctan(2 * cot(b) * ((Minf * np.sin(b)) ** 2 - 1) / \
        (Minf ** 2 * (1.4 + np.cos(2*b)) + 2)) # Deflection Angle [rad]
    M2 = Mn2 / (np.sin(b - d)) # Mach Number After Shockwave
    V2 = (1 + 2 / (0.4 * M2 ** 2)) ** (-0.5) # Velocity After Shockwave
    theta[i] = b
    Vr[i] = V2*np.cos(b - d)
    Vtheta[i] = -V2*np.sin(b - d)

    # Runge-Kutta
    cone_surf.terminal = True
    sol = solve_ivp(
        TM_eq, (b, 0), [Vr[0], Vtheta[0]], t_eval=theta, events=cone_surf, 
        rtol=1e-12, atol=1e-12
    )

    theta = sol.t
    Vr = sol.y[0]
    Vtheta = sol.y[1]

    theta = np.append(theta, sol.t_events[0])
    Vr = np.append(Vr, sol.y_events[0][0, 0])
    Vtheta = np.append(Vtheta, 0.)    # sol.y_events[0][0, 1] is basically zero
    return theta, Vr, Vtheta


def my_streamline(Minf, b, R, dle, Taylor, Dt) -> None:
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

    # Sistima sintetagmenon sti vasi konou
    x = L - x
    z = z - R

    z[-1] = np.interp(0., [x[-1], x[-2]], [z[-1], z[-2]])
    x[-1] = 0.
    
    return x, z


def other_streamline(Minf, b, R, dle, Taylor, Dt) -> None:
    # Runge-Kutta Equations
    # In Vector Form: ds/dt = V
    def strm_func(t, y, Taylor, Vr, Vtheta, L) -> float:
        """
        dr/dt = Vr
        dtheta/dt = Vtheta/r
        y[0] = r
        y[1] = theta
        """
        drdt = np.interp(y[1]*180./np.pi, Taylor.theta, Vr)
        dthetadt = np.interp(y[1]*180./np.pi, Taylor.theta, Vtheta) / y[0]
        return drdt, dthetadt

    def stop_event(t, y, Taylor, Vr, Vtheta, L):
        return y[0] * np.cos(y[1]) - L


    time = np.arange(0, Dt*2*cfg.N_POINTS, Dt)

    b = b*np.pi/180.
    z0 = R - dle
    L = R / np.tan(b)
    h0 = cfg.CP*cfg.ATM.T + ((Minf*cfg.ATM.v_sonic)**2)/2
    Vmax = np.sqrt(2*h0)

    Vr = Taylor.Vr * Vmax
    Vtheta = Taylor.Vtheta * Vmax

    # Runge-Kutte
    stop_event.terminal = True
    sol = solve_ivp(
        fun=strm_func, t_span=(time[0], time[-1]), y0=[z0/np.sin(b), b],
        events=stop_event, t_eval=time, args=(Taylor, Vr, Vtheta, L),
        rtol=1e-8, atol=1e-8
    )

    x = sol.y[0]*np.cos(sol.y[1])
    z = sol.y[0]*np.sin(sol.y[1])

    x = np.append(x, L)
    z = np.append(z, sol.y_events[0][0,0]*np.sin(sol.y_events[0][0,1]))
    # Sistima sintetagmenon sti vasi konou
    x = L - x
    z = z - R

    # z[-1] = np.interp(0., [x[-1], x[-2]], [z[-1], z[-2]])
    # x[-1] = 0.
    
    return x, z



config_create_window()

B = 20
Dtheta = np.deg2rad(-1e-3)
taylor = TaylorMaccoll(Minf=cfg.MINF, b=B)

dint = 1.2

R = 2
Dt = (dint/np.tan(B*np.pi/180)) / \
    (taylor.V[-1]*np.sqrt(2*1000*cfg.ATM.T + ((cfg.MINF*cfg.ATM.v_sonic)**2))) / \
    cfg.N_POINTS

dint /= 10

tic = time.time()
mysol = my_runge_kutta(cfg.MINF, B)
# mysol = my_streamline(cfg.MINF, B, R, dint, taylor, Dt)
toc_my = time.time() - tic

tic = time.time()
othersol = other_runge_kutta(cfg.MINF, B)
# othersol = other_streamline(cfg.MINF, B, R, dint, taylor, Dt)
toc_other = time.time() - tic

print(f'{toc_my=} {toc_other=}')

print('---------------------------')

print(len(mysol[0]), len(othersol[0]))
print(np.max(np.abs(mysol[0] - othersol[0][:-1])))
print(othersol[-1][-1])
