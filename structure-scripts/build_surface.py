#!/usr/bin/env python3
"""BCC / FCC / HCP 表面模型构建（基于 ASE，取向构建方式与 atomsk ``--create ... orient`` 一致）

构建流程
--------
1. 由取向表（或 ``--orient``）给出 x、y、z 三个晶向（立方用 [uvw]，HCP 用四指数 [uvtw]）；
   每个晶向取该方向上的最短点阵平移矢量作为单胞基矢（与 atomsk 相同），
   因此单胞基矢、原子数与 atomsk 一致；所有取向均为右手系（自定义左手系会自动翻转 y）。
2. 标准化取向：a 沿 +x，b 在 xy 平面内，表面法向为 z（与 atomsk 相同）。
3. 终止面：将一个周期内的原子层按 z 排序，``--term K`` 表示在第 K 层下方的层间隙处切开
   （第 K 层成为底层）；``--term auto``（默认）在最大层间隙处切开；``--list-terms`` 列出所有选项。
4. 扩胞 ``-d nx ny nz``（与 atomsk ``-dup`` 一致，nz 为第三基矢的周期数，而非原子层数）。
5. 真空层：顶层与底层（经周期边界）间的法向间隙严格等于 ``-vac``，板层居中。
   ``--vac-mode z``（默认）：第三基矢改为沿 z 垂直方向 (0, 0, h)，便于 LAMMPS 与偶极修正；
   ``--vac-mode c``：沿原第三基矢方向（可能倾斜）伸长，与 atomsk ``-cell add`` 方式一致。

示例
----
下表为各表面一个周期（-d 1 1 1）的层数/原子数（Nb a=3.307、Al a=4.041、Ti a=2.928 c=4.640），
nz 为第三基矢周期数，原子层数 = 每周期层数 × nz：

    体系 表面        每周期层数  每周期原子数  |a| x |b| (Å)
    BCC  100         2          2            3.307 x 3.307
    BCC  110         2          4            4.677 x 3.307
    BCC  111         3          6            4.677 x 8.100
    BCC  112         6          6            2.864 x 4.677
    BCC  123         14         42           2.864 x 21.432
    FCC  100         2          2            2.857 x 2.857
    FCC  110         2          2            2.857 x 4.041
    FCC  111         3          6            2.857 x 4.949
    HCP  basal       2          4            2.928 x 5.071
    HCP  prismatic   4          4            2.928 x 4.640
    HCP  pyrI        6          6            2.928 x 5.487 (gamma=105.48)
    HCP  pyrII       10         20           5.487 x 5.071

(1) 面内 1x1，原子层数 >= 12（能被 12 整除时恰为 12 层）：

    build_surface.py bcc 100 -e Nb -a 3.307 -d 1 1 6                 # 12 层, 12 原子
    build_surface.py bcc 110 -e Nb -a 3.307 -d 1 1 6                 # 12 层, 24 原子
    build_surface.py bcc 111 -e Nb -a 3.307 -d 1 1 4                 # 12 层, 24 原子
    build_surface.py bcc 112 -e Nb -a 3.307 -d 1 1 2                 # 12 层, 12 原子
    build_surface.py bcc 123 -e Nb -a 3.307 -d 1 1 1                 # 14 层, 42 原子
    build_surface.py fcc 100 -e Al -a 4.041 -d 1 1 6                 # 12 层, 12 原子
    build_surface.py fcc 110 -e Al -a 4.041 -d 1 1 6                 # 12 层, 12 原子
    build_surface.py fcc 111 -e Al -a 4.041 -d 1 1 4                 # 12 层, 24 原子
    build_surface.py hcp basal -e Ti -a 2.928 -c 4.640 -d 1 1 6      # 12 层, 24 原子
    build_surface.py hcp prismatic -e Ti -a 2.928 -c 4.640 -d 1 1 3  # 12 层, 12 原子
    build_surface.py hcp pyrI -e Ti -a 2.928 -c 4.640 -d 1 1 2       # 12 层, 12 原子
    build_surface.py hcp pyrII -e Ti -a 2.928 -c 4.640 -d 1 1 2      # 20 层, 40 原子

(2) 原子层数 12（bcc 123 为 14 层、hcp pyrII 为 20 层），原子数 30~70，面内 x、y 尽量等长：

    build_surface.py bcc 100 -e Nb -a 3.307 -d 2 2 6                 # 48 原子, 6.61 x 6.61
    build_surface.py bcc 110 -e Nb -a 3.307 -d 1 2 6                 # 48 原子, 4.68 x 6.61
    build_surface.py bcc 111 -e Nb -a 3.307 -d 2 1 4                 # 48 原子, 9.35 x 8.10
    build_surface.py bcc 112 -e Nb -a 3.307 -d 2 2 2                 # 48 原子, 5.73 x 9.35
    build_surface.py bcc 123 -e Nb -a 3.307 -d 1 1 1                 # 42 原子, 2.86 x 21.43
    build_surface.py fcc 100 -e Al -a 4.041 -d 2 2 6                 # 48 原子, 5.71 x 5.71
    build_surface.py fcc 110 -e Al -a 4.041 -d 2 2 6                 # 48 原子, 5.71 x 8.08
    build_surface.py fcc 111 -e Al -a 4.041 -d 2 1 4                 # 48 原子, 5.71 x 4.95
    build_surface.py hcp basal -e Ti -a 2.928 -c 4.640 -d 2 1 6      # 48 原子, 5.86 x 5.07
    build_surface.py hcp prismatic -e Ti -a 2.928 -c 4.640 -d 2 2 3  # 48 原子, 5.86 x 9.28
    build_surface.py hcp pyrI -e Ti -a 2.928 -c 4.640 -d 3 1 2       # 36 原子, 8.78 x 5.49
    build_surface.py hcp pyrII -e Ti -a 2.928 -c 4.640 -d 1 1 2      # 40 原子, 5.49 x 5.07

其他用法：

    build_surface.py fcc all -e Al -a 4.041 -f vasp lammps
    build_surface.py hcp prismatic --list-terms -e Ti -a 2.928 -c 4.640
    build_surface.py hcp pyrI -e Ti -a 2.928 -c 4.640 -d 1 1 2 --term 1
    build_surface.py fcc custom --orient "[1-10]" "[11-2]" "[111]" -e Cu -a 3.615
    build_surface.py hcp all -e Ti -a 2.928 -c 4.640 --check-atomsk
"""

import argparse
import re
import shutil
import subprocess
import sys
import tempfile
from fractions import Fraction
from functools import reduce
from math import gcd
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.build import bulk, make_supercell
from ase.geometry import cellpar_to_cell
from ase.io import read, write

# 取向表：face -> (x, y, z, 默认 nz, 说明)；均为右手系 (x × y ∥ +z 方向)
# 与原 atomsk 脚本相比，bcc/fcc 110 改为 [-110]（原 [1-10] 为左手系），
# bcc 111 改为 [1-10][11-2]，hcp prismatic 改为 [-1-120]，hcp pyrII 改为 [1-100]
ORIENTS = {
    "bcc": {
        "100": ("[010]", "[001]", "[100]", 6, "(100)"),
        "110": ("[-110]", "[001]", "[110]", 6, "(110)"),
        "111": ("[1-10]", "[11-2]", "[111]", 4, "(111)"),
        "112": ("[11-1]", "[-110]", "[112]", 2, "(112)"),
        "123": ("[11-1]", "[-54-1]", "[123]", 1, "(123)"),
    },
    "fcc": {
        "100": ("[01-1]", "[011]", "[100]", 6, "(100)"),
        "110": ("[-110]", "[001]", "[110]", 6, "(110)"),
        "111": ("[1-10]", "[11-2]", "[111]", 4, "(111)"),
    },
    "hcp": {
        "basal": ("[11-20]", "[-1100]", "[0001]", 6, "Basal (0001)"),
        "prismatic": ("[-1-120]", "[0001]", "[-1100]", 3, "Prismatic I (-1100)"),
        "pyrI": ("[-2110]", "[1-213]", "[01-11]", 2, "Pyramidal I (01-11)"),
        "pyrII": ("[-1-123]", "[1-100]", "[11-22]", 1, "Pyramidal II (11-22)"),
    },
}
ALIASES = {"pyramidalI": "pyrI", "pyramidalII": "pyrII", "prism": "prismatic"}


# ----------------------------------------------------------------------------
# 晶向解析与最短点阵矢量
# ----------------------------------------------------------------------------
def parse_miller(s: str) -> list[int]:
    """解析 "[1-10]"、"[-1-123]" 或 "1,-1,0" / "5 -4 -1" 形式的晶向指数"""
    body = s.strip().strip("[]<>()")
    if re.search(r"[,\s]", body):
        idx = [int(t) for t in re.split(r"[,\s]+", body) if t]
    else:
        idx = [int(t) for t in re.findall(r"-?\d", body)]
    if len(idx) not in (3, 4):
        raise ValueError(f"无法解析晶向 {s!r}")
    return idx


def direction_to_cart(idx: list[int], lattice: str, a: float, c: float | None):
    """常规胞晶向 -> 笛卡尔矢量（HCP 四指数 [uvtw] -> 三指数 [u-t, v-t, w]）"""
    if lattice == "hcp":
        if len(idx) == 4:
            u, v, t, w = idx
            if u + v + t != 0:
                raise ValueError(f"四指数晶向 {idx} 不满足 u+v+t=0")
            idx = [u - t, v - t, w]
        a1 = np.array([a, 0.0, 0.0])
        a2 = np.array([-a / 2, a * np.sqrt(3) / 2, 0.0])
        a3 = np.array([0.0, 0.0, c])
        return idx[0] * a1 + idx[1] * a2 + idx[2] * a3
    if len(idx) != 3:
        raise ValueError(f"立方晶向需要三指数: {idx}")
    return a * np.array(idx, dtype=float)


def shortest_lattice_vector(cart, prim_cell) -> np.ndarray:
    """返回沿 cart 方向的最短点阵平移矢量（以原胞基矢的整数系数表示）"""
    frac = np.linalg.solve(prim_cell.T, cart)
    fr = [Fraction(x).limit_denominator(1000) for x in frac]
    lcm = reduce(lambda p, q: p * q // gcd(p, q), (f.denominator for f in fr), 1)
    ints = [int(f * lcm) for f in fr]
    g = reduce(gcd, (abs(i) for i in ints))
    return np.array(ints) // g


# ----------------------------------------------------------------------------
# 取向体相单胞
# ----------------------------------------------------------------------------
def oriented_cell(lattice, symbol, a, c, orient) -> Atoms:
    """构建以 orient 三个晶向为基矢的取向单胞（a 沿 x，b 在 xy 面内，右手系）"""
    prim = bulk(symbol, lattice, a=a, c=c) if lattice == "hcp" else bulk(symbol, lattice, a=a)
    P = np.array(
        [
            shortest_lattice_vector(direction_to_cart(parse_miller(o), lattice, a, c), prim.cell[:])
            for o in orient
        ]
    )
    if abs(np.linalg.det(P)) < 1e-8:
        raise ValueError(f"晶向 {orient} 共面，无法构成单胞")
    if np.linalg.det(P) < 0:
        print(f"[warn] 取向 {orient} 为左手系，已将 y 取反以保持右手系", file=sys.stderr)
        P[1] *= -1

    atoms = make_supercell(prim, P, wrap=True)
    # 旋转到标准取向；右手系保证这是纯旋转
    atoms.set_cell(cellpar_to_cell(atoms.cell.cellpar()), scale_atoms=True)
    atoms.wrap(eps=1e-6)
    return atoms


def find_layers(atoms: Atoms, tol: float):
    """一个周期内的原子层：返回 [(z, 原子索引列表)]，按 z 升序"""
    z = atoms.positions[:, 2]
    order = np.argsort(z)
    layers = []
    for i in order:
        if layers and z[i] - layers[-1][0] < tol:
            layers[-1][1].append(i)
        else:
            layers.append([z[i], [i]])
    return [(float(np.mean(z[ids])), ids) for _, ids in layers]


def layer_gaps(layers, period_z):
    """gap[k] = 第 k 层与其下方一层（周期性）之间的法向间距"""
    zs = np.array([z for z, _ in layers])
    return np.diff(np.concatenate([[zs[-1] - period_z], zs]))


def resolve_term(term, gaps, tol=1e-3) -> int:
    if term == "auto":
        return int(np.flatnonzero(gaps > gaps.max() - tol)[0])
    k = int(term)
    if not 0 <= k < len(gaps):
        raise ValueError(f"--term 需在 0..{len(gaps) - 1} 之间")
    return k


def set_termination(atoms: Atoms, layers, k: int) -> Atoms:
    """沿第三基矢平移，使第 k 层成为底层"""
    atoms = atoms.copy()
    s = atoms.get_scaled_positions(wrap=False)
    s3_k = np.mean(s[layers[k][1], 2])
    s[:, 2] = (s[:, 2] - s3_k) % 1.0
    s[s[:, 2] > 1 - 1e-6, 2] = 0.0
    atoms.set_scaled_positions(s)
    return atoms


def add_vacuum(slab: Atoms, vacuum: float, mode: str) -> Atoms:
    """板层法向间隙 = vacuum，板层居中"""
    slab = slab.copy()
    z = slab.positions[:, 2]
    thick = z.max() - z.min()
    h = thick + vacuum
    cell = slab.cell[:].copy()
    if mode == "z":
        cell[2] = [0.0, 0.0, h]
    else:  # 沿原第三基矢方向伸长，使其法向高度为 h
        cell[2] = cell[2] * h / cell[2, 2]
    slab.positions[:, 2] += vacuum / 2 - z.min()
    slab.set_cell(cell, scale_atoms=False)
    slab.pbc = (True, True, False)
    slab.wrap(eps=1e-6)
    slab.pbc = True
    return slab


# ----------------------------------------------------------------------------
# atomsk 对比
# ----------------------------------------------------------------------------
def atomsk_oriented(lattice, symbol, a, c, orient) -> Atoms | None:
    exe = shutil.which("atomsk")
    if exe is None:
        print("[warn] 未找到 atomsk，跳过对比", file=sys.stderr)
        return None
    lat_args = [str(a), str(c)] if lattice == "hcp" else [str(a)]
    with tempfile.TemporaryDirectory() as tmp:
        cmd = [exe, "--create", lattice, *lat_args, symbol, "orient", *orient, "out.xsf"]
        r = subprocess.run(cmd, cwd=tmp, input="n\n", capture_output=True, text=True)
        out = Path(tmp) / "out.xsf"
        if not out.exists():
            print(f"[warn] atomsk 运行失败:\n{r.stdout[-500:]}", file=sys.stderr)
            return None
        return read(out)


def same_slab(A: Atoms, B: Atoms, tol=1e-3) -> bool:
    """判断两个板层（同一面内基矢，z 非周期）是否只差一个整体平移"""
    if len(A) != len(B) or not np.allclose(A.cell[:2], B.cell[:2], atol=tol):
        return False
    inv2 = np.linalg.inv(A.cell[:2, :2])

    def rel(X, i):
        d = X.positions - X.positions[i]
        f = d[:, :2] @ inv2
        f -= np.floor(f + tol)
        return np.column_stack([f @ A.cell[:2, :2], d[:, 2]])

    bottom = lambda X: np.flatnonzero(X.positions[:, 2] < X.positions[:, 2].min() + tol)
    ra = rel(A, bottom(A)[0])
    for j in bottom(B):
        rb = rel(B, j)
        diff = ra[:, None, :] - rb[None, :, :]
        f = diff[..., :2] @ inv2
        f -= np.round(f)
        diff[..., :2] = f @ A.cell[:2, :2]
        if np.all(np.linalg.norm(diff, axis=2).min(axis=1) < 10 * tol):
            return True
    return False


def same_structure(A: Atoms, B: Atoms, tol=1e-3) -> bool:
    """同一完整胞（含真空层）的两个板层是否只差一个整体平移"""
    return np.allclose(A.cell[:], B.cell[:], atol=tol) and same_slab(A, B, tol)


def check_atomsk(unit, layers, lattice, symbol, a, c, orient):
    ref = atomsk_oriented(lattice, symbol, a, c, orient)
    if ref is None:
        return
    ok_cell = np.allclose(unit.cell[:], ref.cell[:], atol=1e-4)
    ref_slab = ref.repeat((1, 1, 2))
    match = [
        k for k in range(len(layers))
        if same_slab(set_termination(unit, layers, k).repeat((1, 1, 2)), ref_slab)
    ]
    print(
        f"    atomsk 对比: 原子数 {len(unit)} vs {len(ref)}, 基矢一致={ok_cell}, "
        f"atomsk 默认终止面对应 term={match if match else '无匹配'}"
    )


# ----------------------------------------------------------------------------
# 主流程
# ----------------------------------------------------------------------------
def build_face(args, face):
    if face == "custom":
        if not args.orient:
            raise SystemExit("custom 需要 --orient X Y Z")
        x, y, z = args.orient
        nz, label = 1, f"custom {x}{y}{z}"
    else:
        x, y, z, nz, label = ORIENTS[args.lattice][face]
        if args.orient:
            raise SystemExit("--orient 只能与 custom 一起使用")
    orient = (x, y, z)
    nx, ny, nz = args.dup if args.dup else (1, 1, nz)

    unit = oriented_cell(args.lattice, args.element, args.a, args.c, orient)
    layers = find_layers(unit, args.layer_tol)
    gaps = layer_gaps(layers, unit.cell[2, 2])

    print(f"\n{args.lattice.upper()} {args.element} {label}  orient x={x} y={y} z={z}")
    print(
        f"    单胞: {len(unit)} 原子, {len(layers)} 层/周期, "
        f"|a|={unit.cell.lengths()[0]:.4f} |b|={unit.cell.lengths()[1]:.4f} "
        f"gamma={unit.cell.angles()[2]:.2f}, 周期法向高度={unit.cell[2, 2]:.4f}"
    )
    if args.check_atomsk:
        check_atomsk(unit, layers, args.lattice, args.element, args.a, args.c, orient)

    if args.list_terms:
        auto = resolve_term("auto", gaps)
        print("    term   z(Å)    原子数  切开的层间隙(Å)")
        for k, ((zk, ids), g) in enumerate(zip(layers, gaps)):
            mark = "  <- auto" if k == auto else ""
            print(f"    {k:4d} {zk:8.4f} {len(ids):6d} {g:12.4f}{mark}")
        return

    k = resolve_term(args.term, gaps)
    slab = set_termination(unit, layers, k).repeat((nx, ny, nz))
    slab = add_vacuum(slab, args.vacuum, args.vac_mode)
    slab = slab[np.lexsort((slab.positions[:, 0], slab.positions[:, 1], slab.positions[:, 2]))]

    outdir = Path(args.outdir or f"{args.lattice}-surfaces")
    outdir.mkdir(parents=True, exist_ok=True)
    stem = f"{args.element}_{face if face != 'custom' else 'custom'}"
    written = []
    for fmt in args.format:
        if fmt == "vasp":
            fn = outdir / f"{stem}.vasp"
            write(fn, slab, format="vasp", direct=True, sort=True)
        else:
            fn = outdir / f"{stem}.lmp"
            write(
                fn,
                slab,
                format="lammps-data",
                units="metal",
                atom_style="atomic",
                specorder=[args.element],
                masses=True,
            )
        written.append(str(fn))

    zs = slab.positions[:, 2]
    print(
        f"    板层: dup={nx}x{ny}x{nz}, {len(slab)} 原子, {len(layers) * nz} 层, "
        f"term={k} (切开间隙 {gaps[k]:.4f} Å), 厚度={zs.max() - zs.min():.3f} Å, "
        f"真空={args.vacuum} Å ({args.vac_mode}), 胞高={slab.cell[2, 2]:.3f} Å"
    )
    print(f"    -> {', '.join(written)}")


def main():
    p = argparse.ArgumentParser(
        description="BCC/FCC/HCP 表面模型构建（ASE 实现，取向与 atomsk 一致）",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="可选表面:\n"
        + "\n".join(f"  {lat}: {' '.join(fs)} all custom" for lat, fs in ORIENTS.items()),
    )
    p.add_argument("lattice", choices=list(ORIENTS))
    p.add_argument("faces", nargs="+", help="表面类型，all 表示全部，custom 配合 --orient")
    p.add_argument("-e", "--element", required=True, help="元素符号")
    p.add_argument("-a", type=float, required=True, help="晶格常数 a (Å)")
    p.add_argument("-c", type=float, help="HCP 晶格常数 c (Å)")
    p.add_argument(
        "-d",
        "--dup",
        type=int,
        nargs=3,
        metavar=("NX", "NY", "NZ"),
        help="扩胞（同 atomsk -dup；默认 1 1 NZ，NZ 见取向表）",
    )
    p.add_argument("-vac", "--vacuum", type=float, default=15.0, help="真空层法向厚度 (默认 15 Å)")
    p.add_argument(
        "--vac-mode",
        choices=["z", "c"],
        default="z",
        help="z: 第三基矢改为垂直 z（默认）；c: 沿原第三基矢方向伸长（同 atomsk）",
    )
    p.add_argument("--term", default="auto", help="终止面：auto（最大层间隙处切开）或层编号 K")
    p.add_argument("--list-terms", action="store_true", help="列出各终止面选项后退出")
    p.add_argument("--layer-tol", type=float, default=0.01, help="同一原子层 z 容差 (默认 0.01 Å)")
    p.add_argument("--orient", nargs=3, metavar=("X", "Y", "Z"), help="自定义取向（face=custom）")
    p.add_argument("-f", "--format", nargs="+", choices=["vasp", "lammps"], default=["vasp"])
    p.add_argument("-o", "--outdir", help="输出目录（默认 <lattice>-surfaces）")
    p.add_argument("--check-atomsk", action="store_true", help="与 atomsk 生成的取向单胞对比")
    args = p.parse_args()

    if args.lattice == "hcp" and args.c is None:
        p.error("hcp 需要 -c")
    faces = []
    for f in args.faces:
        f = ALIASES.get(f, f)
        if f == "all":
            faces += list(ORIENTS[args.lattice])
        elif f == "custom" or f in ORIENTS[args.lattice]:
            faces.append(f)
        else:
            p.error(
                f"{args.lattice} 不支持表面 {f!r}，可选: {', '.join(ORIENTS[args.lattice])}, all, custom"
            )

    for face in dict.fromkeys(faces):
        build_face(args, face)


if __name__ == "__main__":
    main()
