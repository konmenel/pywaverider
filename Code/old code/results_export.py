import xlsxwriter as xls
import pathlib
from config_modify import results_load_window


def main() -> None:
    # Load results
    results, results_filename = results_load_window(True)
    WR = results['WR']

    # Save Path
    PATH = (((pathlib.Path(__file__)).parents[1]).joinpath('Results')).joinpath(results_filename)
    PATH.mkdir(parents=True, exist_ok=True)

    # Create Excel files
    for i, ls in enumerate(WR.LS):
        filename = PATH.joinpath("Plane " + str(i) + ".xlsx")
        with xls.Workbook(filename) as workbook:
            worksheet = workbook.add_worksheet()
            if ls != WR.LS[-1]:
                worksheet.write_column(0, 0, ls.x)
                worksheet.write_column(0, 1, ls.y)
                worksheet.write_column(0, 2, ls.z)
            else:
                worksheet.write(0, 0, ls.x)
                worksheet.write(0, 1, ls.y)
                worksheet.write(0, 2, ls.z)


if __name__ == "__mai__":
    main()
