# Copyright (C) 2023 Constantinos Menelaou <https://github.com/konmenel>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# @author: Constantinos Menelaou
# @github: https://github.com/konmenel
# @year: 2023
"""Configuration file.

Configuration Parameters:
    MINF: Free flow Mach number
    H: Altitude [m]
    N_POINTS: Number of points of first plane. The rest of the planes.
              are scaled according to the length of the plane.
    PLANES: Number of planes
    VISCOUS: Boolean that determines whether or not viscous effects will be included.
    METHOD: Either 1 or 2 representing the method that will be used for the viscous forces.
            (1: Van Driest Method, 2: Reference Temperature)
    ATM: ATMOSPHERE_1976 object containing the properties of the air the the given altitude. (From fluids module)
    optVar: The name of the variable that is being optimazed.
    minMax: Either 1 or -1 representing whether it is a maximization (-1) or minimization (1) problem.
    CONFIG: Dictionary that contains all of the above for saving purposes. (Check save module for more info.)

Constants:
    CP: Specific Heat Capacity at constant pressure of air. [J/(kg*K)]
    GAM: Heat Capacity ratio of air.
    MB_AIR: Molecular weight of air. [kg/kmol]
    R_GAS: Universal Gas Constant. [J/(kmol*K)]
"""
from fluids.atmosphere import ATMOSPHERE_1976


# Init Configuration
MINF: float = None
H: float = None
N_POINTS: int = None
PLANES: int = None
# PER_L: float = None
VISCOUS: bool = None
METHOD: int = None
ATM: ATMOSPHERE_1976 = None
optVar: str = None
minMax: int = None
CONFIG: dict = {}

# Constants
CP: float = 1000.
GAM: float = 1.4
MB_AIR: float = 28.97
R_GAS: float = 8314.5
