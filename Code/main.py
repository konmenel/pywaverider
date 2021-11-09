"""The main module of the tool. Running it will perform the calculation of an
optimal Waverider for a given mission.
"""
from os import system
import numpy as np
import waverider as wrdr
import optimization as opti
import config as cfg
from config_modify import config_create_window
from save import save_results
from scipy.optimize import minimize



def main() -> None:
    """The main function of the tool"""
    # Setup configuration
    config_create_window()
    system("title "+ cfg.CONFIG['Viscous Method'])

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
    """The main function of the tool"""
    # Setup configuration
    config_create_window()
    system("title "+ cfg.CONFIG['Viscous Method'])

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
    other_main()
