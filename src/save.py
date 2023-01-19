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
"""Module containing the functions related to saving the results and configuration.

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
import json
import datetime
import pathlib
import xlsxwriter as xls
import config as cfg
import waverider as wrdr
from config_modify import results_load_pickle, results_load_json, modify_config


def save_results(results: dict) -> None:
    """Function that saves the results dictionary in a json file."""
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
    """The old save function. Should not be used. Saved results in pickle object."""
    import pickle

    # Date and Time
    datetime_now = datetime.datetime.now()
    datetime_now = datetime_now.strftime("%Y.%m.%d-%H.%M")

    # Path
    path = ((pathlib.Path(__file__)).parents[1]).joinpath('Results')

    # Save results
    filename = str(path) + "\\" + datetime_now  + "-" + cfg.CONFIG['Viscous Method'] + " results.pickle"
    with open(filename, 'wb') as output:  # Overwrites any existing file.
        pickle.dump(results, output, protocol=4)

    # Save Configuration
    filename = str(path) + "\\" + datetime_now  + "-" + cfg.CONFIG['Viscous Method'] +  " config.json"
    with open(filename, 'w') as output:
        json.dump(cfg.CONFIG, output, sort_keys=True, indent=4)


def pickle_to_json_results() -> None:
    """A way to transform pickle saves (old) to json saves (new)."""
    results, filename = results_load_pickle(second_output=True)

    results['WR'] = results['WR'].todict()

    json_results = {'CONFIG': cfg.CONFIG, 'results': results}

    filename += "_results.json"
    with open(filename, 'w') as output:
        json.dump(json_results, output, indent=4)
   

def export_planes(wr: wrdr.Waverider=None, filename: str=None) -> None:
    """Function that exports the waverider geometry to excel file. 
    One file per plane.
    """
    if wr is None:
        # Load results
        results, filename = results_load_json(second_output=True)
        wr = wrdr.Waverider(wrdr.WRInputs(**results['WR']))

    if wr is not None and filename is None:
        raise Exception('Filename not given.')

    # Save Path
    PATH = (pathlib.Path(__file__)).parents[1] / 'Results' / filename
    PATH.mkdir(parents=True, exist_ok=True)

    # Create Excel files
    for i, ls in enumerate(wr.LS):
        filename = PATH.joinpath(f"Plane {i}.xlsx")
        with xls.Workbook(filename) as workbook:
            worksheet = workbook.add_worksheet()
            if ls != wr.LS[-1]:
                worksheet.write_column(0, 0, ls.x)
                worksheet.write_column(0, 1, ls.y)
                worksheet.write_column(0, 2, ls.z)
            else:
                worksheet.write(0, 0, ls.x)
                worksheet.write(0, 1, ls.y)
                worksheet.write(0, 2, ls.z)


def export_performance() -> None:
    """Runs the aerodynamic performace of a geometry using all methods
    (Inviscid, Van Driest, Ref. Temperature) and stores the performance to
    a file.
    """
    # All the methods
    methods = ['Inviscid', 'Van Driest', 'Ref. Temperature']

    # Laoding the result
    results, filename = results_load_json(second_output=True)
    wr_inputs = wrdr.WRInputs(**results['WR'])
    with xls.Workbook(f'{filename}.xlsx') as workbook:
        worksheet = workbook.add_worksheet()
        worksheet.write(0, 0, 'L_D')
        worksheet.write(1, 0, 'L')
        worksheet.write(2, 0, 'D')
        worksheet.write(3, 0, 'CL')
        worksheet.write(4, 0, 'CD')
        worksheet.write(5, 0, 'S')
        worksheet.write(6, 0, 'V')
        worksheet.write(7, 0, 'volEff')
        worksheet.write(8, 0, 's')
        worksheet.write(9, 0, 'l')
        worksheet.write(10, 0, 's_l')
        worksheet.write(11, 0, 'Cp_base')
        worksheet.write(12, 0, 'P_base')
        worksheet.write(13, 0, 'D_base')
        worksheet.write(14, 0, 'S_base')


        for i, method in enumerate(methods, 1):
            cfg.CONFIG['Viscous Method'] = method
            modify_config()

            wr_perf = wrdr.Waverider(wr_inputs)
            
            worksheet.write(0, i, wr_perf.L_D)
            worksheet.write(1, i, wr_perf.L)
            worksheet.write(2, i, wr_perf.D)
            worksheet.write(3, i, wr_perf.CL)
            worksheet.write(4, i, wr_perf.CD)
            worksheet.write(5, i, wr_perf.S)
            worksheet.write(6, i, wr_perf.V)
            worksheet.write(7, i, wr_perf.volEff)
            worksheet.write(8, i, wr_perf.s)
            worksheet.write(9, i, wr_perf.l)
            worksheet.write(10, i, wr_perf.s_l)
            worksheet.write(11, i, wr_perf.Cp_base)
            worksheet.write(12, i, wr_perf.P_base)
            worksheet.write(13, i, wr_perf.D_base)
            worksheet.write(14, i, wr_perf.S_base)


if __name__ == "__main__":
    export_performance()
