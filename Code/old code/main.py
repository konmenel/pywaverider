"""
The main module of the tool. Running it will perform the calculation of an
optimal Waverider for a given mission.
"""
import pathlib
import pickle
import json
import datetime
from os import system
import numpy as np
import Waverider as wrdr
import optimization as opti
import config as cfg
from config_modify import config_create_window
from save import old_save


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

    old_save(results)

    # # Date and Time
    # DATETIME = datetime.datetime.now()
    # DATETIME = DATETIME.strftime("%Y.%m.%d-%H.%M")

    # # Path
    # PATH = ((pathlib.Path(__file__)).parents[1]).joinpath('Results')

    # # Save results
    # filename = str(PATH) + "\\" + DATETIME  + "-" + cfg.CONFIG['Viscous Method'] + " results.pickle"
    # with open(filename, 'wb') as output:  # Overwrites any existing file.
    #     pickle.dump(results, output, protocol=4)

    # # Save Configuration
    # filename = str(PATH) + "\\" + DATETIME  + "-" + cfg.CONFIG['Viscous Method'] +  " config.json"
    # with open(filename, 'w') as output:
    #     json.dump(cfg.CONFIG, output, sort_keys=True, indent=4)

    # Print and Plot Results
    opti.results_to_screen(results)

if __name__ == "__main__":
    main()
