import os


# TODO: 根据计算类型，删除不同的文件
def clean_vasp():
    """删除 VASP 计算过程中生成的输出文件"""

    files = [
        "AECCAR0",
        "AECCAR1",
        "AECCAR2",
        "WAVECAR",
        "CHG",
        "CHGCAR",
        "REPORT",
        "PROCAR",
        "LOCPOT",
        "DOSCAR",
        "EIGENVAL",
        "IBZKPT",
        "CONTCAR",
        "XDATCAR",
        "OSZICAR",
        "OUTCAR",
        "vasprun.xml",
    ]
    for f in files:
        try:
            os.remove(f)
        except OSError:
            pass
