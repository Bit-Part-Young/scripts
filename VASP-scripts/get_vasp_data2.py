#!/usr/bin/env python3

"""直接解析 VASP 输出，获取能量、原子数、离子步数、收敛状态和实际耗时。

仅依赖 pandas（表格输出），无需 pymatgen。energy 单位为 eV，energy_pa
为 energy / natoms，单位为 eV/atom。ionic_steps 为已完成离子步数。
time_cost 优先读取 OUTCAR 的 Elapsed time；未结束时累加 LOOP+ 的 real time。
未写完的 XML 使用最后一个完整离子步，converged=False；第一步未完成时
能量和 time_cost 为 N/A。
"""

import argparse
import glob
import re
import sys
import warnings
import xml.etree.ElementTree as ET
from math import fsum
from pathlib import Path
from xml.parsers import expat

import pandas as pd

PARAMETERS = {"NELM", "NSW", "IBRION", "EDIFFG", "ALGO", "LEPSILON", "ML_LMLFF"}


def _float(text):
    """同时支持 XML/OUTCAR 中的 E 和 Fortran D 指数。"""
    return float(text.replace("D", "E").replace("d", "e"))


def _parameters(element):
    values = {}
    for child in element:
        name = child.get("name", "").strip()
        if child.tag == "i" and name in PARAMETERS:
            value = (child.text or "").strip()
            kind = child.get("type")
            if kind == "int":
                value = int(value)
            elif kind == "logical":
                value = value.strip(".").upper() == "T"
            elif kind != "string":
                value = _float(value)
            values[name] = value
        elif child.tag not in {"i", "v"}:
            nested = _parameters(child)
            # 与 pymatgen 一致，response functions 中的默认值不覆盖主参数。
            if name == "response functions":
                nested = {k: v for k, v in nested.items() if k not in values}
            values.update(nested)
    return values


def _energies(element):
    return {
        child.get("name"): _float(child.text.strip())
        for child in element.findall("energy/i")
    }


def _converged(incar, parameters, electronic_steps, nsteps):
    """沿用 Vasprun 的电子步/离子步收敛判据。"""
    nelm = parameters["NELM"]
    algo = incar.get("ALGO", "")
    if incar.get("ML_LMLFF"):
        electronic = True
    elif incar.get("LEPSILON"):
        index = 1
        energy_keys = {"e_wo_entrp", "e_fr_energy", "e_0_energy"}
        while (
            index < len(electronic_steps)
            and set(electronic_steps[index]) == energy_keys
        ):
            index += 1
        electronic = index + 1 != nelm
    elif algo == "Exact" and incar.get("NELM") == 1:
        electronic = True
    else:
        electronic = (0 if algo.lower() == "chi" else len(electronic_steps)) < nelm

    nsw = parameters.get("NSW", 0)
    ibrion = parameters.get("IBRION", -1 if nsw in (-1, 0) else 0)
    if ibrion == 0 or (ibrion in {1, 2} and parameters.get("EDIFFG", 1) == 0):
        ionic = nsw <= 1 or nsteps == nsw
    else:
        ionic = nsw <= 1 or nsteps < nsw
    return electronic and ionic


def _read_vasprun(filename):
    incar, parameters = {}, {}
    natoms = None
    nsteps = 0
    energy = None
    electronic_steps = []
    complete = False
    stack = []
    try:
        with filename.open("rb") as handle:
            for event, element in ET.iterparse(handle, events=("start", "end")):
                if event == "start":
                    stack.append(element)
                    continue
                # 完成一个顶层节点后立即释放，避免积累整个轨迹或 DOS。
                if len(stack) == 2:
                    if element.tag == "incar":
                        incar = _parameters(element)
                    elif element.tag == "parameters":
                        parameters = _parameters(element)
                    elif element.tag == "atominfo":
                        natoms = int(element.findtext("atoms"))
                    elif element.tag == "calculation":
                        ionic_energy = _energies(element)
                        electronic_steps = [
                            _energies(step) for step in element.findall("scstep")
                        ]
                        energy = ionic_energy.get("e_0_energy")
                        if electronic_steps:
                            final = electronic_steps[-1]
                            # 某些 VASP 版本将熵项误写进离子步 e_0_energy。
                            corrected = round(
                                final["e_0_energy"]
                                - final["e_fr_energy"]
                                + ionic_energy["e_fr_energy"],
                                8,
                            )
                            if energy is None or abs(energy - corrected) > 1e-7:
                                energy = corrected
                        nsteps += 1
                    stack[0].remove(element)
                    element.clear()
                elif len(stack) == 1:
                    complete = element.tag == "modeling"
                stack.pop()
    except ET.ParseError as error:
        # 只容忍文件末尾未完成；标签错配等损坏仍然报错。
        eof_errors = {
            expat.errors.codes[message]
            for message in (
                expat.errors.XML_ERROR_NO_ELEMENTS,
                expat.errors.XML_ERROR_UNCLOSED_TOKEN,
                expat.errors.XML_ERROR_PARTIAL_CHAR,
                expat.errors.XML_ERROR_UNCLOSED_CDATA_SECTION,
            )
        }
        if error.code not in eof_errors:
            raise ValueError(f"{filename}: {error}") from error
        warnings.warn(
            f"{filename}: XML 未写完，已完成 {nsteps} 个离子步；" "converged=False",
            stacklevel=2,
        )

    if not complete and nsteps == 0:
        return float("nan"), natoms, False, 0
    if natoms is None or natoms <= 0 or energy is None or nsteps == 0:
        raise ValueError(f"{filename}: 缺少原子数或完整离子步的能量数据")
    if "NELM" not in parameters:
        raise ValueError(f"{filename}: 缺少 NELM 参数，无法判定收敛")
    converged = complete and _converged(incar, parameters, electronic_steps, nsteps)
    return energy, natoms, converged, nsteps


def _read_timing(filename):
    """返回最终墙钟时间和已完成离子步的实际时间（秒）。"""
    elapsed = None
    step_times = []
    number = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?"
    elapsed_pattern = re.compile(rf"^\s*Elapsed time \(sec\):\s*({number})\s*$")
    loop_pattern = re.compile(
        rf"^\s*LOOP\+:\s*cpu time\s*{number}:\s*real time\s*({number})\s*$"
    )
    with filename.open(errors="replace") as handle:
        for line in handle:
            if not line.endswith("\n"):
                # 文件正在写入时，末尾可能只有半行计时记录。
                continue
            match = elapsed_pattern.match(line)
            if match:
                elapsed = _float(match.group(1))
            match = loop_pattern.match(line)
            if match:
                step_times.append(_float(match.group(1)))
    return elapsed, step_times


def seconds_to_hms(seconds: float):
    """将秒数转换为时分秒格式"""

    hours, remainder = divmod(seconds, 3600)
    minutes, seconds_new = divmod(remainder, 60)

    return f"{int(hours)}h-{int(minutes)}min-{int(seconds_new)}s"


def get_vasp_data(path: str = "."):
    """返回单行 DataFrame；总能量和每原子能量均保留五位小数。"""
    directory = Path(path).expanduser()
    energy, natoms, converged, xml_steps = _read_vasprun(directory / "vasprun.xml")
    elapsed, step_times = _read_timing(directory / "OUTCAR")
    # LOOP+ 在离子步结束时写出；XML 可能尚未刷出对应的 </calculation>。
    ionic_steps = max(xml_steps, len(step_times))
    if ionic_steps == 0:
        time_cost = "N/A"
    elif elapsed is not None:
        time_cost = seconds_to_hms(elapsed)
    elif step_times:
        time_cost = seconds_to_hms(fsum(step_times))
        warnings.warn(
            f"{directory / 'OUTCAR'}: 尚无最终 Elapsed time，"
            f"time_cost 使用 {len(step_times)} 个完整离子步的 real time 之和",
            stacklevel=2,
        )
    else:
        time_cost = "N/A"
        warnings.warn(
            f"{directory / 'OUTCAR'}: 缺少 Elapsed time 和完整 LOOP+ 计时，"
            "time_cost=N/A",
            stacklevel=2,
        )

    data_dict = {
        "energy": energy,
        "energy_pa": energy / natoms if natoms else float("nan"),
        "natoms": natoms if natoms is not None else pd.NA,
        "ionic_steps": ionic_steps,
        "converged": converged,
        "time_cost": time_cost,
    }

    df = pd.DataFrame(data_dict, index=[0]).round(5)

    return df


def expand_paths(paths, seen=None):
    """展开目录通配符并去重；普通文件（如 shell 展开 * 后的文件）跳过。"""
    if seen is None:
        seen = set()
    for argument in paths:
        argument = str(Path(argument).expanduser())
        matches = (
            sorted(glob.glob(argument)) if glob.has_magic(argument) else [argument]
        )
        if not matches:
            raise ValueError(f"没有匹配的目录: {argument}")
        for match in matches:
            directory = Path(match)
            if directory.is_file():
                continue
            key = directory.resolve()
            if key not in seen:
                seen.add(key)
                yield directory


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Get VASP data directly from vasprun.xml and OUTCAR.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="Examples: %(prog)s . | %(prog)s Mo Nb | %(prog)s *",
    )

    parser.add_argument(
        "path",
        type=str,
        default=["."],
        nargs="*",
        help="VASP directories; . = current directory, '*' = immediate subdirectories",
    )

    args = parser.parse_args(argv)
    frames = []
    failed = False
    seen = set()
    for argument in args.path:
        try:
            for directory in expand_paths([argument], seen):
                try:
                    with warnings.catch_warnings(record=True) as messages:
                        frame = get_vasp_data(directory)
                    for message in messages:
                        print(f"Warning: {message.message}", file=sys.stderr)
                except (OSError, ValueError, KeyError, TypeError) as error:
                    print(f"Error: {directory}: {error}", file=sys.stderr)
                    failed = True
                    continue
                frame.insert(0, "path", str(directory))
                frames.append(frame)
        except (OSError, ValueError) as error:
            print(f"Error: {error}", file=sys.stderr)
            failed = True
    if frames:
        print(pd.concat(frames, ignore_index=True).to_string(index=False, na_rep="N/A"))
    else:
        print("Error: 没有可读取的 VASP 计算目录", file=sys.stderr)
        failed = True
    return int(failed)


if __name__ == "__main__":
    sys.exit(main())
