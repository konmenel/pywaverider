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
from os import system
import config as cfg
from config_modify import results_load_json
from optimization import results_to_screen


def main() -> None:
    """The main function of the file."""
    # User Input
    results = results_load_json()
    system("title "+ cfg.CONFIG['Viscous Method'])

    # Printing and Plotting results
    ans = input("Plot 3D? (type y/n)\n")
    while not (ans == "y" or ans == "n"):
        ans = input("Type y or n\n")

    if ans == "n":
        results_to_screen(results, plot_3D=False)

    elif ans == "y":
        ans = input("Online plot? (y/n)\n")
        while not (ans == "y" or ans == "n"):
            ans = input("Type y or n\n")
        
        if ans == "y":
            results_to_screen(results)
        elif ans == "n":
            results_to_screen(results, htmlExport=True)


if __name__ == "__main__":
    main()
