"""
A file that when run, plots and prints the results of an optimization to screen.
"""
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
