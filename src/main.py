#!/usr/bin/env python3

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
import numpy as np
import waverider as wrdr
import optimization as opti
import config as cfg
from config_modify import config_create_window
from save import save_results
from scipy.optimize import minimize

# TODO: make it command-line program

def print_header() -> None:
    """Just print the header copyright license when the program starts"""

    HEADER = """=======================================================================
PyWaverider
Copyright (C) 2023 Constantinos Menelaou <https://github.com/konmenel>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

@author: Constantinos Menelaou
@github: https://github.com/konmenel
@year: 2023
=======================================================================
    """
    print(HEADER)

def main() -> None:
    """The main function of the tool"""
    # Setup configuration
    config_create_window()
    print_header()

    # Initial guess
    b = 21.48
    s = 33.73
    l = 72.4
    per_l = 0.3719
    yle = np.array([4.258, 6.793, 20.39])
    zle = np.array([19.99, 17.37])
    guess = wrdr.WRInputs(b, s, l, per_l, yle, zle)
    guess = wrdr.Waverider(guess)

    # Optimization Run
    results = opti.NelderMead(guess, cfg.optVar, cfg.minMax)

    save_results(results)

    # Print and Plot Results
    opti.results_to_screen(results)


def other_main() -> None:
    """Just experimenting with scipy's Nelder-Mead"""
    # Setup configuration
    config_create_window()

    # Initial guess
    b = 21.48
    s = 33.73
    l = 72.4
    per_l = 0.3719
    yle = [4.258, 6.793, 20.39]
    zle = [19.99, 17.37]
    x0 = [b, s, l, per_l, *yle, *zle]

    # Optimization Run
    options={
        'maxiter': 10000,
        'return_all': True,
        'xatol': 1e-8,
        'fatol': 1e-8,
        'adaptive': True,
        'disp': True
    }

    results = minimize(
        opti.optimization_function, x0, method='Nelder-Mead',
        options=options
    )

    # save_results(results)

    # # Print and Plot Results
    # opti.results_to_screen(results)



if __name__ == "__main__":
    main()
