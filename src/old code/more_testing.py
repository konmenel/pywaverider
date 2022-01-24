"""Just a testing file"""
import config as cfg
import Waverider as wrdr
from config_modify import results_load_json


def testing() -> None:
    """The main function of the file."""
    results = results_load_json()
    print(cfg.CONFIG)
    print(cfg.MINF)
    print(cfg.PLANES)
    
    my_wr = wrdr.Waverider(wrdr.WRInputs(**results['WR']))
    print("--------------------------------")
    print(my_wr.inputs())
    print("-----------------------------")
    print(my_wr)

if __name__ == "__main__":
    testing()
