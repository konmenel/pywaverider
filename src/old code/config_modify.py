"""
Module the contains the function to create or load a configuration.
"""
import sys
import pickle
import json
import pathlib
import tkinter as tk
from tkinter import filedialog
from fluids.atmosphere import ATMOSPHERE_1976 as atmos
import config as cfg


def config_create_window() -> None:
    """The function that creates a window to setup a new configuration."""
    # Constants
    WINDOW_HEIGHT = 355
    WINDOW_WIDTH = 465
    SIZE = str(WINDOW_WIDTH) + 'x' + str(WINDOW_HEIGHT)
    DARK_GRAY = '#222222'
    LIGHT_GRAY = '#E1E1E1'
    DEF_FONT = ('Helvetica', 10)

    # Click Functions
    def click_start():
        cfg.CONFIG['Mach'] = float(inBoxMach.get())
        cfg.CONFIG['Altitude [km]'] = float(inBoxHeight.get())
        cfg.CONFIG['Planes'] = int(inBoxPlanes.get())
        cfg.CONFIG['Number of Points'] = int(inBoxN.get())
        # cfg.CONFIG['Percentage of Line Segment'] = float(inBoxPerL.get())/100
        cfg.CONFIG['Viscous Method'] = viscousVar.get()
        cfg.CONFIG['Optimization Parameter'] = ParamVar.get()
        cfg.CONFIG['Type of Problem'] = TypeVar.get()
        modify_config()
        window.destroy()
        
    # Create the window object
    window = tk.Tk()
    window.config(bg=DARK_GRAY)
    window.geometry(SIZE)
    window.title('Mission Configuration')
    window.resizable(False, False)

    # Enter, esc keys bind
    window.bind('<Return>', lambda e: click_start())
    window.bind('<Escape>', lambda e: sys.exit())

    # Label Text
    Mach_txt = 'Select Design Mach Number:'
    Height_txt = 'Select Cruise Altitude [km]:'
    Planes_txt = 'Select Number of Planes:'
    N_txt = 'Select Number of Points for Each Plane:'
    # PerL_txt = 'Select Percentage of Straight Shockwave [%]:'
    Method_txt = 'Select Aerodynamic Analysis Method:'
    Param_txt = 'Select Optimization Parameter:'
    Type_txt = 'Select Type of Problem:'

    # Creating the Labels
    labelMach = tk.Label(
        window, text=Mach_txt, fg=LIGHT_GRAY, bg=DARK_GRAY, font=DEF_FONT,
        anchor="w", width=40, pady=10
    )
    labelHeight = tk.Label(
        window, text=Height_txt, fg=LIGHT_GRAY, bg=DARK_GRAY, font=DEF_FONT,
        anchor="w", width=40, pady=10
    )
    labelPlanes = tk.Label(
        window, text=Planes_txt, fg=LIGHT_GRAY, bg=DARK_GRAY, font=DEF_FONT,
        anchor="w", width=40, pady=10
    )
    labelN = tk.Label(
        window, text=N_txt, fg=LIGHT_GRAY, bg=DARK_GRAY, font=DEF_FONT,
        anchor="w", width=40, pady=10
    )
    # labelPerL = tk.Label(
    #     window, text=PerL_txt, fg=LIGHT_GRAY, bg=DARK_GRAY, font=DEF_FONT,
    #     anchor="w", width=40, pady=10
    # )
    labelMethod = tk.Label(
        window, text=Method_txt, fg=LIGHT_GRAY, bg=DARK_GRAY, font=DEF_FONT,
        anchor="w", width=40, pady=10
    )
    labelParam = tk.Label(
        window, text=Param_txt, fg=LIGHT_GRAY, bg=DARK_GRAY, font=DEF_FONT,
        anchor="w", width=40, pady=10
    )
    labelType = tk.Label(
        window, text=Type_txt, fg=LIGHT_GRAY, bg=DARK_GRAY, font=DEF_FONT,
        anchor="w", width=40, pady=10
    )

    # Labels on window
    labelMach.grid(row=0, column=0)
    labelHeight.grid(row=1, column=0)
    labelPlanes.grid(row=2, column=0)
    labelN.grid(row=3, column=0)
    # labelPerL.grid(row=4, column=0)
    labelMethod.grid(row=4, column=0)
    labelParam.grid(row=5, column=0)
    labelType.grid(row=6, column=0)

    # Creating Input Boxes
    inBoxMach = tk.Entry(window, width=14)
    inBoxMach.insert(0, 5)
    inBoxHeight = tk.Entry(window, width=14)
    inBoxHeight.insert(0, 33.4)
    inBoxPlanes = tk.Entry(window, width=14)
    inBoxPlanes.insert(0, 80)
    inBoxN = tk.Entry(window, width=14)
    inBoxN.insert(0, 800)
    # inBoxPerL = tk.Entry(window, width=14)
    # inBoxPerL.insert(0, 40)

    # Input Boxes on window
    inBoxMach.grid(row=0, column=1)
    inBoxHeight.grid(row=1, column=1)
    inBoxPlanes.grid(row=2, column=1)
    inBoxN.grid(row=3, column=1)
    # inBoxPerL.grid(row=4, column=1)

    # Creating Drop Down Menus
    viscousVar = tk.StringVar(window)
    viscousOpt = ['Inviscid', 'Van Driest', 'Ref. Temperature']
    viscousVar.set(viscousOpt[0]) # default value
    viscousMenu = tk.OptionMenu(window, viscousVar, *viscousOpt)
    viscousMenu.config(width=14)
    ParamVar = tk.StringVar(window)
    ParamOpt = ['L/D', 'L', 'D', 'CL', 'CD', 'S', 'V', 'Volumetric Eff.']
    ParamVar.set(ParamOpt[0]) # default value
    ParamMenu = tk.OptionMenu(window, ParamVar, *ParamOpt)
    ParamMenu.config(width=14)
    TypeVar = tk.StringVar(window)
    TypeOpt = ['Maximization', 'Minimization']
    TypeVar.set(TypeOpt[0]) # default value
    TypeMenu = tk.OptionMenu(window, TypeVar, *TypeOpt)
    TypeMenu.config(width=14)

    # Drop Down Menus on window
    viscousMenu.grid(row=4, column=1)
    ParamMenu.grid(row=5, column=1)
    TypeMenu.grid(row=6, column=1)

    # Start, Cancel and Load Buttons buttons
    startButton = tk.Button(window, text='Start', command=click_start, width=10)
    cancelButton = tk.Button(window, text='Cancel', command=sys.exit, width=10)
    loadButton = tk.Button(window, text='Load', command=lambda: config_load_window(window), width=10)

    # Buttons on window
    startButton.place(relx=0.3, rely=0.91)
    cancelButton.place(relx=0.5, rely=0.91)
    loadButton.place(relx=0.8, rely=0.91)

    # Main loop for window
    window.mainloop()


def config_load_window(window: tk.Tk) -> None:
    """The function that load an existing configuration from a json"""
    window.filename = tk.filedialog.askopenfilename(title="Select Configuration",
                    filetypes=(('JSON files', '*.json'), ('All files', '*.*')))
    
    with open(window.filename, 'r') as file:
        results_dict = json.load(file)

    cfg.CONFIG = results_dict["CONFIG"]


def modify_config() -> None:
    """The function that modifies the configuration file."""
    # Modify config.py
    cfg.MINF = cfg.CONFIG['Mach']
    cfg.H = cfg.CONFIG['Altitude [km]']
    cfg.PLANES = cfg.CONFIG['Planes']
    cfg.N_POINTS = cfg.CONFIG['Number of Points']
    # cfg.PER_L = cfg.CONFIG['Percentage of Line Segment']

    if cfg.CONFIG['Viscous Method'] == 'Inviscid':
        cfg.VISCOUS = False
        cfg.METHOD = 0
    elif cfg.CONFIG['Viscous Method'] == 'Van Driest':
        cfg.VISCOUS = True
        cfg.METHOD = 1
    elif cfg.CONFIG['Viscous Method'] == 'Ref. Temperature':
        cfg.VISCOUS = True
        cfg.METHOD = 2

    if cfg.CONFIG['Optimization Parameter'] == 'L/D':
        cfg.optVar = 'L_D'
    elif cfg.CONFIG['Optimization Parameter'] == 'Volumetric Eff.':
        cfg.optVar = 'volEff'
    else:
        cfg.optVar = cfg.CONFIG['Optimization Parameter']

    if cfg.CONFIG['Type of Problem'] == 'Maximization':
        cfg.minMax = -1
    else:
        cfg.minMax = 1

    cfg.ATM = atmos(cfg.H*1e3)

    cfg.CONFIG['Pressure [Pa]'] = cfg.ATM.P
    cfg.CONFIG['Temperature [K]'] = cfg.ATM.T
    cfg.CONFIG['Density [kg/m^3]'] = cfg.ATM.rho
    cfg.CONFIG['Dynamic Viscosity [Ns/m^2]'] = cfg.ATM.mu


def results_load_pickle(second_output: bool=False) -> None:
    """The function that load a saved pickle result."""
    window = tk.Tk()
    window.withdraw()
    window.filename = filedialog.askopenfilename(title="Select Results",
                    initialdir=((pathlib.Path(__file__)).parents[1]).joinpath('Results'),
                    filetypes=(('pickle files', '*.pickle'), ('All files', '*.*')))
    # Load results
    with open(window.filename, 'rb') as file:
        results = pickle.load(file)
    
    # Load Config
    with open((window.filename).replace('results.pickle', 'config.json'), 'rb') as file:
        cfg.CONFIG = json.load(file)
    modify_config()

    if second_output:
        return results, window.filename.replace(' results.pickle', '')

    return results


def results_load_json(second_output: bool=False) -> None:
    """The function that load a saved json result."""
    window = tk.Tk()
    window.withdraw()
    window.filename = filedialog.askopenfilename(title="Select Results",
                    initialdir=((pathlib.Path(__file__)).parents[1]).joinpath('Results'),
                    filetypes=(('json files', '*.json'), ('All files', '*.*')))
    
    # Load json
    with open(window.filename, 'r') as file:
        results = json.load(file)

    cfg.CONFIG = results['CONFIG']
    modify_config()

    if second_output:
        return results['results'], window.filename.replace('_results.json', '')

    return results['results']



if __name__ == "__main__":
    pass
    