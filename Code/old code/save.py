"""
Module containing the functions related to saving the results and configuration.

About the saves:
The JSON file contains a dictionary with 2 fields ('results' and 'CONFIG'). 

The "results" field is a dictionary that is the output of the Nelder-Mead optimization. 
Check optimization module for more info.

The "CONFIG" field is a dictionary with the following fields:
    "Mach": The free flow Mach number.
    "Altitude [km]": The altitude of the mission in km.
    "Planes": The number of planes.
    "Number of Points": Number of points of first plane. The rest of the planes.
                        are scaled according to the length of the plane.
    "Viscous Method": The viscous effects model. Either "Inviscid", "Ref. Temperature" and "Van Driest"
    "Optimization Parameter": The name of the parameter that is being optimized.
    "Type of Problem": Either "Maximization" or "Minimization".
    "Pressure [Pa]": The pressure at that altitude in Pa.
    "Temperature [K]": The temperature at that altitude in K.
    "Density [kg/m^3]": The density at that altitude in kg/m^3.
    "Dynamic Viscosity [Ns/m^2]": The dynamic viscosity at that altitude in Ns/m^2.
"""
import config as cfg
import json
import datetime
import pathlib
import pickle
from config_modify import results_load_pickle


def save_results(results: dict) -> None:
    # Date and Time
    DATETIME = datetime.datetime.now()
    DATETIME = DATETIME.strftime("%Y.%m.%d-%H.%M")

    # Path
    PATH = ((pathlib.Path(__file__)).parents[1]).joinpath('Results')

    filename = f"{PATH}\\{DATETIME}-{cfg.CONFIG['Viscous Method']}_results.json"

    json_results = {'config': cfg.CONFIG, 'results': results}

    # Dump
    with open(filename, 'w') as output:
        json.dump(json_results, output, indent=4)


def old_save(results: dict) -> None:

    # Date and Time
    DATETIME = datetime.datetime.now()
    DATETIME = DATETIME.strftime("%Y.%m.%d-%H.%M")

    # Path
    PATH = ((pathlib.Path(__file__)).parents[1]).joinpath('Results')

    # Save results
    filename = str(PATH) + "\\" + DATETIME  + "-" + cfg.CONFIG['Viscous Method'] + " results.pickle"
    with open(filename, 'wb') as output:  # Overwrites any existing file.
        pickle.dump(results, output, protocol=4)

    # Save Configuration
    filename = str(PATH) + "\\" + DATETIME  + "-" + cfg.CONFIG['Viscous Method'] +  " config.json"
    with open(filename, 'w') as output:
        json.dump(cfg.CONFIG, output, sort_keys=True, indent=4)


def pickle_to_json_results() -> None:
    """A way to transform pickle saves (old) to json saves (new)"""
    results, filename = results_load_pickle(second_output=True)

    results['WR'] = results['WR'].todict()

    json_results = {'CONFIG': cfg.CONFIG, 'results': results}

    filename += "_results.json"
    with open(filename, 'w') as output:
        json.dump(json_results, output, indent=4)
   


if __name__ == "__main__":
    pickle_to_json_results()
