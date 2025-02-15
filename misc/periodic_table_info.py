"""打印元素周期表中元素的基本性质"""

from pymatgen.core.periodic_table import Element

elements = [
    "Ag",
    "Au",
    "Cu",
    "Mg",
    "Ni",
    "Pb",
    "Pd",
    "Pt",
    "Ta",
    "Ti",
    "Al",
    "Cr",
    "V",
    "W",
    "Mo",
    "Zr",
    "Nb",
    "Ta",
    "Hf",
]

for element in elements:
    ele = Element(element)
    # 元素 相对原子质量 原子半径（无统一值，不一定正确） 熔点
    print(
        f"{element:>2} {round(float(str(ele.atomic_mass)[:-4]), 3):<7}  {str(ele.atomic_radius)[:-4]:<5} {str(ele.melting_point)[:-2]:<6}"
    )
