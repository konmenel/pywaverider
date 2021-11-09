import config as cfg
from os import system
from config_modify import results_load_window
from optimization import results_to_screen


def main() -> None:
    # User Input
    results = results_load_window()
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
