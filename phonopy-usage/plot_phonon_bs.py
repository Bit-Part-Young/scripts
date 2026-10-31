#!/usr/bin/env python3

"""
声子谱绘制：多组数据（DFT / NEP / EAM ...）对比 + 实验散点，支持多子图

k-path 预设（与文献常用路径一致）：
    bcc: Γ-H-P-Γ-N        fcc: Γ-X-U|K-Γ-L（U≡K，标为 K）    hcp: Γ-K-M-Γ-A
    *_full: Setyawan-Curtarolo 完整路径

高对称点以晶体学坐标定义（立方：笛卡尔坐标，单位 2π/a；六方：120° 设定下的倒格矢分数坐标），
再按每个数据自身的原胞换算为分数坐标，因此不依赖原胞的取法（vaspkit / phonopy / SC 约定均可）。

用法：
    # 单图
    plot_phonon_bs.py -i phonopy_Nb_dft.yaml phonopy_Nb_nep.yaml -l DFT NEP -k bcc \\
        -e Nb-bcc_phonon_exp.csv -t BCC-Nb -o Nb_phonon.png

    # 多子图（配置文件，路径相对于配置文件所在目录）
    plot_phonon_bs.py -c phonon_bands.yaml

配置文件示例：
    style: mlip.mplstyle          # 可选，matplotlib style 文件或名称（可为列表）
    output: ../figures/phonon_bands.png
    figsize: [7.2, 2.6]
    ncols: 3                      # 可选，每行子图数
    fig_legend: {ncols: 8}        # 可选，图上方公共图例（fig.legend 参数）；此时子图默认不画图例
    panels:
      - title: FCC-Al
        kpath: fcc
        ylim: [0, 12]
        exp: ../../2-Exp-data/phonon/Al-fcc_phonon.csv   # 或列表：[{file: a.csv, label: Expt. A}, b.csv]
        datasets:                 # 先列出的画在上层，图例按此顺序
          - {label: NEP, file: ../../3-LAMMPS/runs/phonon/NEP/phonopy_Al.yaml, color: "#e8291c"}
          - {label: DFT, file: ../../1-DFT-data/Phonon/0-yaml/phonopy_Al_dft.yaml, color: "#2b72bc"}

实验数据 CSV（'#' 开头为注释），两种格式任选：
    segment,xi,freq       # segment 如 G-X、K-G（用 k-path 中的点名，G 即 Γ），xi ∈ [0, 1] 为段内位置
    x_norm,freq           # x_norm ∈ [0, 1] 为整条路径上的归一化位置（直接从文献图中数字化时用）
"""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from matplotlib.axes import Axes
from phonopy import load
from phonopy.api_phonopy import Phonopy

# 未提供 style 时的默认参数（期刊单 / 双栏实际尺寸下的字号与线宽）
DEFAULT_RC = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "mathtext.fontset": "custom",
    "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic",
    "mathtext.cal": "Arial",
    "font.size": 9,
    "axes.linewidth": 0.8,
    "lines.linewidth": 1.2,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.major.size": 0,
    "ytick.major.size": 3,
    "ytick.minor.size": 1.5,
    "ytick.minor.visible": True,
    "ytick.right": True,
    "legend.frameon": False,
    "legend.handlelength": 1.6,
    "savefig.dpi": 600,
    "savefig.bbox": "tight",
}

# 未指定颜色时依次使用：蓝（DFT）、红（NEP）、灰（其他势）、绿、紫
DEFAULT_COLORS = ["#2b72bc", "#e8291c", "#b0b0b0", "#2ca02c", "#9467bd"]
# 多组实验数据依次使用的样式：橙色空心圆、黑色空心三角、紫色空心方块
EXP_STYLES = [
    dict(marker="o", mec="#f28e1c"),
    dict(marker="^", mec="black"),
    dict(marker="s", mec="#9467bd"),
]
EXP_COMMON = dict(ms=4, mfc="none", mew=0.8, ls="none")

# 立方：笛卡尔坐标，单位 2π/a（a 为惯用晶胞晶格常数）
CUBIC_POINTS = {
    "bcc": {"G": [0, 0, 0], "H": [0, 1, 0], "P": [1 / 2, 1 / 2, 1 / 2], "N": [1 / 2, 1 / 2, 0]},
    "fcc": {
        "G": [0, 0, 0],
        "X": [0, 1, 0],
        "W": [1 / 2, 1, 0],
        "K": [3 / 4, 3 / 4, 0],
        "U": [1 / 4, 1, 1 / 4],
        "L": [1 / 2, 1 / 2, 1 / 2],
    },
}
# 原胞体积 = a^3 / n
CUBIC_ATOMS_PER_CONV = {"bcc": 2, "fcc": 4}
# 六方：倒格矢分数坐标（a、b 夹角 120°）
HEX_POINTS = {
    "G": [0, 0, 0],
    "A": [0, 0, 1 / 2],
    "M": [1 / 2, 0, 0],
    "L": [1 / 2, 0, 1 / 2],
    "K": [1 / 3, 1 / 3, 0],
    "H": [1 / 3, 1 / 3, 1 / 2],
}

# path: 连续段列表；alias: 断点处合并标签的替换（U 与 K 等价，文献中只标 K）
KPATHS = {
    "bcc": {"lattice": "bcc", "path": [["G", "H", "P", "G", "N"]]},
    "fcc": {"lattice": "fcc", "path": [["G", "X", "U"], ["K", "G", "L"]], "alias": {"U|K": "K"}},
    "hcp": {"lattice": "hcp", "path": [["G", "K", "M", "G", "A"]]},
    "bcc_full": {"lattice": "bcc", "path": [["G", "H", "N", "G", "P", "H"], ["P", "N"]]},
    "fcc_full": {
        "lattice": "fcc",
        "path": [["G", "X", "W", "K", "G", "L", "U", "W", "L", "K"], ["U", "X"]],
    },
    "hcp_full": {
        "lattice": "hcp",
        "path": [["G", "M", "K", "G", "A", "L", "H", "A"], ["L", "M"], ["K", "H"]],
    },
}


def tex_label(name: str) -> str:
    return r"$\Gamma$" if name in ("G", "Γ", "GAMMA") else name


def load_phonopy(
    phonopy_yaml: str | Path,
    force_sets: str | Path | None = None,
    force_constants: str | Path | None = None,
) -> Phonopy:
    """读取 phonopy yaml；力 / 力常数可在 yaml 内，也可另外指定 FORCE_SETS / FORCE_CONSTANTS"""

    return load(
        phonopy_yaml=phonopy_yaml,
        force_sets_filename=force_sets,
        force_constants_filename=force_constants,
        is_nac=False,
    )


def get_qpoint_coords(phonon: Phonopy, kpath: str) -> dict[str, np.ndarray]:
    """高对称点 -> 当前原胞倒格矢下的分数坐标"""

    lattice = KPATHS[kpath]["lattice"]
    cell = np.array(phonon.primitive.cell)

    if lattice in CUBIC_POINTS:
        a = (abs(np.linalg.det(cell)) * CUBIC_ATOMS_PER_CONV[lattice]) ** (1 / 3)
        # q_frac_i = a_i · k_cart，k_cart 不含 2π（与 phonopy 约定一致）
        return {k: cell @ np.array(u) / a for k, u in CUBIC_POINTS[lattice].items()}

    a1, a2 = cell[0], cell[1]
    gamma = np.degrees(np.arccos(a1 @ a2 / np.linalg.norm(a1) / np.linalg.norm(a2)))
    if not np.isclose(gamma, 120, atol=0.5):
        raise ValueError(f"hcp k-path assumes gamma = 120 deg, got {gamma:.2f} deg")
    return {k: np.array(v, dtype=float) for k, v in HEX_POINTS.items()}


def get_xtick_labels(kpath: str) -> list[str]:
    """不连续段在断点处合并为 'X|Y'（phonopy 中前后两段共享同一 distance）"""

    preset = KPATHS[kpath]
    labels = list(preset["path"][0])
    for segment in preset["path"][1:]:
        labels[-1] = f"{labels[-1]}|{segment[0]}"
        labels.extend(segment[1:])

    alias = preset.get("alias", {})
    labels = [alias.get(label, label) for label in labels]
    return ["|".join(tex_label(s) for s in label.split("|")) for label in labels]


def run_band_structure(phonon: Phonopy, kpath: str, npoints: int = 101) -> dict:
    """返回 distances / frequencies（按段的列表）及每段的起止点名"""

    coords = get_qpoint_coords(phonon, kpath)
    bands, seg_names = [], []
    for segment in KPATHS[kpath]["path"]:
        for start, end in zip(segment[:-1], segment[1:]):
            bands.append(np.linspace(coords[start], coords[end], npoints))
            seg_names.append((start, end))

    phonon.run_band_structure(bands, is_band_connection=True)
    band = phonon.band_structure
    return {
        "distances": [np.asarray(d) for d in band.distances],
        "frequencies": [np.asarray(f) for f in band.frequencies],
        "segments": seg_names,
    }


def normalize_exp(exp) -> list[dict]:
    """exp 可为 None、单个文件，或文件 / {file, label, style} 的列表"""

    if exp is None:
        return []
    items = exp if isinstance(exp, list) else [exp]
    return [dict(item) if isinstance(item, dict) else {"file": item} for item in items]


def read_exp_data(exp_fn: str | Path, ref: dict) -> tuple[np.ndarray, np.ndarray]:
    """将实验数据映射到参考数据（第一组）的 x 坐标上"""

    df = pd.read_csv(exp_fn, comment="#", skipinitialspace=True)
    seg_bounds = {seg: (d[0], d[-1]) for seg, d in zip(ref["segments"], ref["distances"])}

    if "segment" in df.columns:
        x = []
        for seg, xi in zip(df["segment"], df["xi"]):
            start, end = (s.strip() for s in seg.replace("Γ", "G").split("-"))
            if (start, end) in seg_bounds:
                x0, x1 = seg_bounds[(start, end)]
            elif (end, start) in seg_bounds:
                x1, x0 = seg_bounds[(end, start)]
            else:
                raise ValueError(f"{exp_fn}: segment '{seg}' not in k-path")
            x.append(x0 + xi * (x1 - x0))
        x = np.array(x)
    elif "x_norm" in df.columns:
        x = df["x_norm"].to_numpy() * ref["distances"][-1][-1]
    else:
        raise ValueError(f"{exp_fn}: need columns 'segment,xi,freq' or 'x_norm,freq'")

    return x, df["freq"].to_numpy()


def plot_panel(
    ax: Axes,
    datasets: list[dict],
    kpath: str,
    exp: str | Path | list | None = None,
    title: str | None = None,
    ylim: tuple[float, float] | None = None,
    npoints: int = 101,
    legend_loc: str = "best",
    legend: bool = True,
):
    """绘制一个子图；datasets 中每项：label, file, [color, ls, lw, zorder, force_sets, force_constants]"""

    ref = None
    for i, ds in enumerate(datasets):
        phonon = load_phonopy(ds["file"], ds.get("force_sets"), ds.get("force_constants"))
        bs = run_band_structure(phonon, kpath, npoints)
        # 各数据晶格常数略有差异，统一使用第一组的 x 坐标，保证高对称点对齐
        ref = ref or bs
        color = ds.get("color", DEFAULT_COLORS[i % len(DEFAULT_COLORS)])
        for j, (x, freq) in enumerate(zip(ref["distances"], bs["frequencies"])):
            lines = ax.plot(
                x,
                freq,
                color=color,
                ls=ds.get("ls", "-"),
                lw=ds.get("lw"),
                zorder=ds.get("zorder", 2 + len(datasets) - i),
            )
            if j == 0:
                lines[0].set_label(ds["label"])

    for k, item in enumerate(normalize_exp(exp)):
        x, freq = read_exp_data(item["file"], ref)
        style = {**EXP_COMMON, **EXP_STYLES[k % len(EXP_STYLES)], **item.get("style", {})}
        ax.plot(x, freq, label=item.get("label", "Expt."), zorder=10 + k, **style)

    xticks = [ref["distances"][0][0]] + [d[-1] for d in ref["distances"]]
    xticks = list(dict.fromkeys(np.round(xticks, 10)))
    for xt in xticks[1:-1]:
        ax.axvline(xt, color="0.5", ls=(0, (3, 3)), lw=0.6, zorder=1)
    ax.set_xticks(xticks, get_xtick_labels(kpath))
    ax.set_xlim(xticks[0], xticks[-1])
    ax.set_ylim(ylim if ylim is not None else (0, None))
    if ax.get_ylim()[0] < 0:  # 有虚频时画零频线
        ax.axhline(0, color="0.5", lw=0.6, zorder=1)

    if title:
        ax.text(0.04, 0.95, title, transform=ax.transAxes, ha="left", va="top")
    if legend:
        ax.legend(loc=legend_loc)


def plot_figure(
    panels: list[dict],
    output: str | Path,
    figsize: tuple[float, float] | None = None,
    ylabel: str = "Frequency (THz)",
    panel_labels: bool = True,
    npoints: int = 101,
    ncols: int | None = None,
    fig_legend: dict | None = None,
):
    """fig_legend：各子图图例合并为一个图上方的公共图例（按标签去重），值为 fig.legend 的参数"""

    n = len(panels)
    ncols = ncols or n
    nrows = -(-n // ncols)
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=figsize or (2.6 * ncols, 2.4 * nrows),
        layout="constrained",
        squeeze=False,
    )
    for ax in axes.flat[n:]:
        ax.set_visible(False)

    for i, (ax, panel) in enumerate(zip(axes.flat, panels)):
        plot_panel(
            ax,
            datasets=panel["datasets"],
            kpath=panel["kpath"],
            exp=panel.get("exp"),
            title=panel.get("title"),
            ylim=panel.get("ylim"),
            npoints=npoints,
            legend_loc=panel.get("legend_loc", "best"),
            legend=panel.get("legend", fig_legend is None),
        )
        if i % ncols == 0:
            ax.set_ylabel(ylabel)
        if panel_labels and n > 1:
            ax.set_title(f"({chr(97 + i)})", loc="left", x=-0.12, fontweight="bold")

    if fig_legend is not None:
        handles = {}
        for ax in axes.flat[:n]:
            for h, label in zip(*ax.get_legend_handles_labels()):
                handles.setdefault(label, h)
        # 实验散点放最后
        labels = sorted(handles, key=lambda s: s.startswith("Expt"))
        kw = {"loc": "outside upper center", "ncols": len(labels), **fig_legend}
        fig.legend([handles[s] for s in labels], labels, **kw)

    fig.savefig(output)
    print(f"Phonon band structure figure saved to {output}")


def resolve_config(config_fn: str | Path) -> dict:
    """读取配置；文件路径相对于配置文件目录，缺失的数据文件给出警告并跳过"""

    config_fn = Path(config_fn).resolve()
    root = config_fn.parent
    with open(config_fn) as f:
        config = yaml.safe_load(f)

    def resolve(p):
        return None if p is None else (root / p).resolve()

    config["output"] = resolve(config.get("output", "phonon_bands.png"))
    if "style" in config:
        styles = config["style"] if isinstance(config["style"], list) else [config["style"]]
        config["style"] = [str(resolve(s)) if (root / s).exists() else s for s in styles]

    for panel in config["panels"]:
        datasets = []
        for ds in panel["datasets"]:
            for key in ("file", "force_sets", "force_constants"):
                if ds.get(key) is not None:
                    ds[key] = resolve(ds[key])
            if not ds["file"].exists():
                print(f"[{panel.get('title')}] skip {ds['label']}: {ds['file']} not found", file=sys.stderr)
                continue
            datasets.append(ds)
        if not datasets:
            raise FileNotFoundError(f"[{panel.get('title')}] no dataset found")
        panel["datasets"] = datasets

        exps = []
        for item in normalize_exp(panel.get("exp")):
            item["file"] = resolve(item["file"])
            if not item["file"].exists():
                print(f"[{panel.get('title')}] skip Expt.: {item['file']} not found", file=sys.stderr)
                continue
            exps.append(item)
        panel["exp"] = exps

    return config


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Phonon band structure comparison plot with phonopy.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("-c", "--config", metavar="YAML", help="multi-panel figure config")
    parser.add_argument(
        "-i", "--input_fn", nargs="+", metavar="FILE", help="phonopy yaml file(s) with forces"
    )
    parser.add_argument("-l", "--labels", nargs="+", help="legend labels (default: file stems)")
    parser.add_argument("-k", "--kpath", default="bcc", choices=list(KPATHS), help="k-path preset")
    parser.add_argument("-e", "--exp", metavar="CSV", help="experimental data csv")
    parser.add_argument("-t", "--title", help="panel title, e.g. BCC-Nb")
    parser.add_argument("--ylim", nargs=2, type=float, metavar=("YMIN", "YMAX"))
    parser.add_argument("-s", "--style", nargs="+", help="matplotlib style file(s) or name(s)")
    parser.add_argument("-n", "--npoints", type=int, default=101, help="q-points per segment")
    parser.add_argument("-o", "--output_fn", default="phonon_bs.png", help="output figure")
    args = parser.parse_args()

    if args.config:
        config = resolve_config(args.config)
        style = config.get("style")
        panels = config["panels"]
        output = config["output"]
        figsize = config.get("figsize")
        ylabel = config.get("ylabel", "Frequency (THz)")
        npoints = config.get("npoints", args.npoints)
        ncols = config.get("ncols")
        fig_legend = config.get("fig_legend")
    elif args.input_fn:
        labels = args.labels or [Path(fn).stem for fn in args.input_fn]
        if len(labels) != len(args.input_fn):
            parser.error("number of --labels must match --input_fn")
        style = args.style
        panels = [
            {
                "kpath": args.kpath,
                "title": args.title,
                "ylim": args.ylim,
                "exp": args.exp,
                "datasets": [{"label": l, "file": f} for l, f in zip(labels, args.input_fn)],
            }
        ]
        output, figsize, npoints = args.output_fn, (3.4, 2.8), args.npoints
        ylabel, ncols, fig_legend = "Frequency (THz)", None, None
    else:
        parser.error("either --config or --input_fn is required")

    with plt.rc_context(DEFAULT_RC), plt.style.context(style or []):
        plot_figure(
            panels, output, figsize=figsize, ylabel=ylabel, npoints=npoints, ncols=ncols,
            fig_legend=fig_legend,
        )
