import matplotlib.pyplot as plt

import waverider as wrdr
import config as cfg
from config_modify import modify_config



cfg.CONFIG = {
    "Altitude [km]": 33.4,
    "Density [kg/m^3]": 0.010864279396722926,
    "Dynamic Viscosity [Ns/m^2]": 1.5051532347594733e-05,
    "Mach": 5.0,
    "Number of Points": 800,
    "Optimization Parameter": "L/D",
    "Planes": 80,
    "Pressure [Pa]": 723.7741536850451,
    "Temperature [K]": 232.0811914150587,
    "Type of Problem": "Maximization",
    "Viscous Method": "Ref. Temperature"
}
modify_config()

wr_in = {
    "b": 18.056262690116643,
    "s": 33.195353512213124,
    "l": 79.99994760609384,
    "per_l": 0.20000004898783563,
    "yle": [
        5.577090436069006,
        7.153491258236458,
        28.47187820394068
    ],
    "zle": [
        17.96887835524319,
        14.446383541177944
    ]
}

wr = wrdr.Waverider(wrdr.WRInputs(**wr_in))

fig, (ax1, ax2) = plt.subplots(1, 2)
