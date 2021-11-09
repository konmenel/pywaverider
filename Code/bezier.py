"""
The module that contains the bezier curve class.
"""
import numpy as np
from scipy.special import comb
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from typing import Any, Union, Iterable


class Bezier:
    """
    Object that represents a 2D Bezier curve.
    Arguments:  (must be same length):
        indVar: Coordinates of point for the indepented variable
        depVar: Coordinates of point for the depented variable
        axes: axes of points (list of string ex.: ['x', 'y'])
    """
    def __init__(self, indVar: Iterable[float], depVar: Iterable[float], axes: Iterable[str]) -> None:
        self.axes = list(axes)
        self.vars = {axes[0]: np.array(indVar),
                    axes[1]: np.array(depVar)}
        self.N = len(indVar)
        self.P = np.array([indVar, depVar])

        # for i, axis in enumerate(axes):
        #     if i == 0:
        #         setattr(self, axis, indVar)
        #     else:
        #         setattr(self, axis, depVar)

    def __getattr__(self, name: str) -> Any:
        try:
            return self.vars[name]
        except KeyError:
            raise AttributeError(name)

    def curve(self, var: str, t: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Method that calculates the value of a point on the curve for a given "t".
        Arguments:
            var: the axis of the output ('x', 'y' or 'z')
            t: Bezier parameter for the point (between 0 and 1)
        """
        if isinstance(t, (list, np.ndarray)):
            c = np.zeros(len(t))
        else:
            c = 0
         
        if isinstance(t, list):
            t = np.array(t)
        
        # P = getattr(self, var)
        P = self.vars[var]

        for i in range(self.N):
            coeff = comb(self.N-1, i)
            c += coeff * (t ** i) * ((1 - t) ** (self.N - (i + 1))) * P[i]
        return c

    def der1(self, var: str, t: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Method that calculates the first derivative of a point on the curve for a given "t".
        Arguments:
            var: the axis of the output ('x', 'y' or 'z')
            t: Bezier parameter for the point (between 0 and 1)
        """
        if isinstance(t, (list, np.ndarray)):
            c = np.zeros(len(t))
        else:
            c = np.zeros(1)
        t = np.array(t)

        # P = getattr(self, var)
        P = self.vars[var]

        for i in range(self.N - 1):
            coeff = comb(self.N-2, i)
            c += coeff * (t ** i) * ((1 - t) ** (self.N - (i + 2))) * (P[i+1] - P[i])
        c *= self.N - 1
        return c
    
    def der2(self, var: str, t: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Method that calculates the second derivative of a point on the curve for a given "t".
        Arguments:
            var: the axis of the output ('x', 'y' or 'z')
            t: Bezier parameter for the point (between 0 and 1)
        """
        if isinstance(t, (list, np.ndarray)):
            c = np.zeros(len(t))
        else:
            c = np.zeros(1)
        t = np.array(t)

        # P = getattr(self, var)
        P = self.vars[var]

        for i in range(self.N - 2):
            coeff = comb(self.N-3, i)
            c += coeff * (t ** i) * ((1 - t) ** (self.N - (i + 3))) * (P[i+2] - 2 * P[i+1] + P[i])
        c *= (self.N - 1) * (self.N - 2)
        return c

    def line_inter(self, k: float, b: float) -> float:
        """
        Method that calculates the parameter 't' at the intersection of a line the Bezier curve
        Arguments:
            k: slope of line
            b: Intersection of the line with y-axis
        """
        xpoly = np.zeros(self.N)
        ypoly = np.zeros(self.N)

        for i in range(1, self.N+1):
            for j in range(i):
                xpoly[self.N-i] += comb(i-1, j) * self.P[0][j] * (-1) ** j
                ypoly[self.N-i] += comb(i-1, j) * self.P[1][j] * (-1) ** j
            xpoly[self.N-i] *= comb(self.N-1, i-1) * (-1) ** (i - 1)
            ypoly[self.N-i] *= comb(self.N-1, i-1) * (-1) ** (i - 1)

        inter_poly = ypoly - k * xpoly
        inter_poly[-1] -= b
        roots = np.roots(inter_poly)
        t = roots.real[np.abs(roots.imag) <= 1e-5]
        t = t[t >= 0]
        t = t[t <= 1]
        if t.size == 0:
            return t
        else:
            return np.max(t)

    def normal_inter(self, b: float, var: str) -> float:
        """
        Method that calculates the parameter 't' at the intersection of a normal line the
        Bezier curve.
        Arguments:
            b: value of the normal line (ex: y = b)
            var: Axis normal to line ('x', 'y' or 'z')
        """

        # P = np.array(getattr(self, var))
        P = self.vars[var]
        
        xpoly = np.zeros(self.N)

        for i in range(1, self.N+1):
            for j in range(i):
                xpoly[self.N-i] += comb(i-1, j) * P[j] * (-1) ** j
            xpoly[self.N-i] *= comb(self.N-1, i-1) * (-1) ** (i-1)

        if isinstance(b, (list, np.ndarray)):
            t = np.empty(len(b))
            lastTerm = xpoly[-1]
            for i, a in enumerate(b):
                xpoly[-1] = lastTerm - a
                roots = np.roots(xpoly)
                t_ = roots.real[np.abs(roots.imag) <= 1e-5]
                t_ = t_[t_ >= 0]
                t_ = t_[t_ <= 1]
                if np.all(t_ == 0.):
                    t[i] = 0.
                else:
                    t[i] = np.max(t_)
        else:
            xpoly[-1] -= b
            roots = np.roots(xpoly)
            t = roots.real[np.abs(roots.imag) <= 1e-5]
            t = t[t >= 0]
            t = t[t <= 1]
            if t. size != 0:
                t = np.max(t)
        return t

    def normal_phi(self, t: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Method that calculates the angle of the normal of the curve at a given point
        from the horizontal axis. (ex: arctan(k) in y = k*x + b)
        Arguments:
            t: Bezier parameter for the point (between 0 and 1)
        """
        dzdy = self.der1('z', t) / self.der1('y', t)
        if isinstance(t, (list, np.ndarray)):
            return np.pi / 2 * np.ones(len(t)) - np.abs(np.arctan(dzdy))

        return np.pi / 2 - np.abs(np.arctan(dzdy))

    def normal_b(self, t: float) -> float:
        """
        Method that calculates the intersection of the normal of the curve at a given point
        from the veritcal axis. (ex: b in y = k*x + b)
        Arguments:
            t: Bezier parameter for the point (between 0 and 1)
        """
        return self.curve('z', t) - np.tan(self.normal_phi(t)) * self.curve('y', t)
 
    def radius_curv(self, t: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Method that calculates the radius of curvature of the curve at a given point
        Arguments:
            t: Bezier parameter for the point (between 0 and 1)
        """
        if isinstance(t, (list, np.ndarray)):
            dens = np.abs(self.der1('y', t) * self.der2('z', t) - \
                self.der2('y', t) * self.der1('z', t))
            rc = np.empty(len(dens))
            for i, den in enumerate(dens):
                if den == 0:
                    rc[i] = float(1e4)
                else:
                    rc[i] = ((self.der1('z', t[i]) ** 2 + self.der1('y', t[i]) ** 2) ** (1.5)) / den
        else:
            den = np.abs(self.der1('y', t) * self.der2('z', t) - self.der2('y', t) * self.der1('z', t))
            if den == 0:
                rc = float(1e4)
            else:
                rc = ((self.der1('z', t) ** 2 + self.der1('y', t) ** 2) ** (1.5)) / den
        return rc

    def plot_radius(self, ax: Axes=None, **param) -> None:
        """
        Method that plots the radius of curvature depending on the parameters given
        Arguments:
            ax: Matplotlib Axes object to be ploted on. In not given a new one will be created
            t: if t is specified the radius will be plotted for the given t (can be either float or list-like)
            x, y or z: if one of those is specified the radius on the given point. If 2 or more points of
                       the curve have the same value the one with the maximum t will be plotted.
        """
        if 't' in param:
            t = param['t']
        
        elif self.axes[0] in param:
            if param[self.axes[0]] == self.vars[self.axes[0]][0]:
                t = 0.
            
            elif param[self.axes[0]] == self.vars[self.axes[0]][-1]:
                t = 1.
            
            else:
                t = self.normal_inter(param[self.axes[0]], self.axes[0])

        elif self.axes[1] in param:
            if param[self.axes[1]] == self.vars[self.axes[1]][0]:
                t = 0.

            elif param[self.axes[1]] == self.vars[self.axes[1]][-1]:
                t = 1.

            else:
                t = self.normal_inter(param[self.axes[1]], self.axes[1])

        r_c = self.radius_curv(t)
        if np.isinf(r_c):
            r_c = 1e3

        x = np.empty(2)
        y = np.empty(2)

        x[0] = self.curve(self.axes[0], t)
        y[0] = self.curve(self.axes[1], t)
        x[1] = x[0] - np.cos(self.normal_phi(t)) * r_c
        y[1] = y[0] - np.sin(self.normal_phi(t)) * r_c

        if ax is None:
            _, ax = plt.subplots()

        ax.plot(x, y, '--xk', linewidth=0.5)

    def plot(self, ax: Axes=None, mirrored: bool=False, *args, **kwargs) -> None:
        """
        Method that plots the Bezier curve.
        Arguments:
            ax: Matplotlib Axes object to be ploted on. In not given a new one will be created
            mirrored: if True the plot will be mirrored on the vertical axis. By default False
            *args, **kwargs: arguments and key arguments of matplotlib plots.
        """
        t = np.linspace(0, 1, 100)

        if ax is None:
            _, ax = plt.subplots()

        if mirrored:
            ax.plot(-self.curve(self.axes[0], t), self.curve(self.axes[1], t), *args, **kwargs)
        else:
            ax.plot(self.curve(self.axes[0], t), self.curve(self.axes[1], t), *args, **kwargs)
        
        ax.set_xlabel(self.axes[0])
        ax.set_ylabel(self.axes[1])


if __name__ == '__main__':
    def testing() -> None:
        t = np.linspace(0, 1, 100)
        B = Bezier((10, 30, 2, 14), (0, 2, 10, 5), ('y', 'z'))

        RC = B.radius_curv(np.linspace(0, 1, 20))
        print(B.y)
        print(B.vars)
        print(RC)
        print(B.der1('y', 0.5))

        # normal x = 15
        normal = 15
        y_n = np.array([0, 10])
        x_n = np.array([15, 15])
        
        tint = B.normal_inter(normal, 'y')
        yint = B.curve('y', tint)
        zint = B.curve('z', tint)
        print(f't = {tint}, y = {yint}, z = {zint}')

        #plt.plot(x, y)
        B.plot()
        plt.scatter(B.y, B.z)
        plt.grid()

        plt.figure()
        plt.plot(B.der1('z', t))
        plt.grid()

        plt.show()

    testing()
