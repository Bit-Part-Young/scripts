import os
import stat
import sys
from os.path import basename


def print_help(
    argv_dict: dict = {"test": "description of this argv"},
    description: str = "This is the description of python srcipt",
) -> None:
    """
    print sys argv parse help message.
    """

    usage_str = f"\nusage: {sys.argv[0]} [-h --help]"

    # check if the script is executable
    script_path = os.path.abspath(sys.argv[0])
    if os.access(script_path, os.X_OK) and stat.S_ISREG(os.stat(script_path).st_mode):
        usage_str = f"\nusage: {basename(sys.argv[0])} [-h --help]"

    for key in argv_dict.keys():
        usage_str += f" [{key.upper()}]"
    print(f"{usage_str}\n")

    print(f"{description}\n")

    print("optional arguments:")

    print(f"  {'-h, --help'.ljust(22)}  show this help message and exit")

    for key, value in argv_dict.items():
        print(f"  {key.upper().ljust(22)}  {value}")

    print("\nAuthor: YSL.")

    sys.exit()


# test
# argv_dict = {
#     "poscar_file": "POSCAR file, supported format: *POSCAR*, *.vasp, *.poscar",
#     "type": "cord type, supported format: frac, cart",
# }
# description = "POSCAR file coordinates format conversion calling pymatgen."
# print_help(argv_dict=argv_dict, description=description)
