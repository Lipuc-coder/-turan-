"""
双部图 Turán 数计算器 — 批量模式
================================================
给定固定的禁图 H，自动为指定范围内的所有 (n,m) 计算 ex(n,m;H)，
为每个案例生成一张可视化图，并打印带有结构注释的摘要表格。

每个极值图的结构分析：
  - 连通分量数量
  - 每个分量：是否为 K_{a,b}？（完全二部图）
  - 若非完全二部图：报告边密度

用法
-----
运行脚本并跟随提示，或导入并调用 batch_compute()。
"""

import os
import sys
import time
import pulp
from itertools import permutations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import networkx as nx

plt.rcParams["font.family"] = "DejaVu Sans"


# ══════════════════════════════════════════════════════════════════
#  共享辅助函数  (与单案例版本相同)
# ══════════════════════════════════════════════════════════════════

def _perm(n, r):
    result = 1
    for i in range(r):
        result *= (n - i)
    return result

def _transpose(H):
    p, q = len(H), len(H[0])
    return [[H[u][v] for u in range(p)] for v in range(q)]

def _get_edges(H):
    p, q = len(H), len(H[0])
    return [(u, v) for u in range(p) for v in range(q) if H[u][v]], p, q

def _orientations(H, n, m):
    H_edges, p, q = _get_edges(H)
    result = []
    if n >= p and m >= q:
        result.append((H, p, q))
    HT = _transpose(H)
    if (q != p or HT != H) and n >= q and m >= p:
        result.append((HT, q, p))
    return result


# ══════════════════════════════════════════════════════════════════
#  ILP 求解器
# ══════════════════════════════════════════════════════════════════

def compute_ex(n, m, H, time_limit=120):
    """""
    计算 ex(n,m;H)。返回 (max_edges, matrix) 或抛出 RuntimeError。
    检查 H 的两种二部定向。
    """

    orientations = _orientations(H, n, m)
    if not orientations:
        return n * m, [[1] * m for _ in range(n)]

    prob = pulp.LpProblem("BT", pulp.LpMaximize)
    x = [[pulp.LpVariable(f"x{i}_{j}", cat="Binary") for j in range(m)]
         for i in range(n)]
    prob += pulp.lpSum(x[i][j] for i in range(n) for j in range(m))

    seen_sigs = set()
    for (H_mat, p, q) in orientations:
        H_edges, _, _ = _get_edges(H_mat)
        eH = len(H_edges)
        for lm in permutations(range(n), p):
            for rm in permutations(range(m), q):
                sig = frozenset((lm[u], rm[v]) for u, v in H_edges)
                if sig not in seen_sigs:
                    seen_sigs.add(sig)
                    prob += (pulp.lpSum(x[lm[u]][rm[v]] for u, v in H_edges)
                             <= eH - 1)

    def _attempt(make_solver):
        try:
            prob.solve(make_solver())
            return pulp.LpStatus[prob.status]
        except Exception:
            return None

    status = (
        _attempt(lambda: pulp.HiGHS(msg=0, timeLimit=time_limit))
        or _attempt(lambda: pulp.PULP_CBC_CMD(msg=0, timeLimit=time_limit))
    )
    if status != "Optimal":
        raise RuntimeError(f"ILP failed: {status}")

    matrix = [[int(round(pulp.value(x[i][j]) or 0)) for j in range(m)]
              for i in range(n)]
    return sum(matrix[i][j] for i in range(n) for j in range(m)), matrix


def verify_H_free(matrix, H):
    """蛮力检查矩阵是否无 H（两种定向）。返回 bool。"""
    n, m = len(matrix), len(matrix[0])
    for (H_mat, p, q) in _orientations(H, n, m):
        H_edges, _, _ = _get_edges(H_mat)
        for lm in permutations(range(n), p):
            for rm in permutations(range(m), q):
                if all(matrix[lm[u]][rm[v]] for u, v in H_edges):
                    return False
    return True


# ══════════════════════════════════════════════════════════════════
#  结构分析
# ══════════════════════════════════════════════════════════════════

def analyze_structure(matrix):
    """
    分析以 n×m 矩阵表示的极值图的结构。

    返回 dict：
      n_components : int
      components   : 分量描述字符串列表
      summary      : 可读的一行摘要
    """
    n, m = len(matrix), len(matrix[0])
    x_nodes = {f"x{i}" for i in range(n)}
    y_nodes = {f"y{j}" for j in range(m)}

    G = nx.Graph()
    for i in range(n): G.add_node(f"x{i}")
    for j in range(m): G.add_node(f"y{j}")
    for i in range(n):
        for j in range(m):
            if matrix[i][j]:
                G.add_edge(f"x{i}", f"y{j}")

    descs = []
    for comp_nodes in nx.connected_components(G):
        sub  = G.subgraph(comp_nodes)
        X    = sorted([v for v in comp_nodes if v in x_nodes],
                      key=lambda s: int(s[1:]))
        Y    = sorted([v for v in comp_nodes if v in y_nodes],
                      key=lambda s: int(s[1:]))
        a, b = len(X), len(Y)

        if a == 0:
            descs.append(f"K_{{0,{b}}} (isolated Y-vertices)")
            continue
        if b == 0:
            descs.append(f"K_{{{a},0}} (isolated X-vertices)")
            continue

        actual   = sub.number_of_edges()
        complete = a * b
        if actual == complete:
            descs.append(f"K_{{{a},{b}}}")
        else:
            # Report which edges are missing for small components
            density = actual / complete
            if complete <= 30:
                missing = [(xi, yj) for xi in X for yj in Y
                           if not sub.has_edge(xi, yj)]
                miss_str = ", ".join(f"({v[1:]},{w[1:]})" for v, w in missing[:5])
                if len(missing) > 5:
                    miss_str += f", ...({len(missing)} total)"
                descs.append(
                    f"partial {a}x{b} [{actual}/{complete} edges, "
                    f"missing: {miss_str}]"
                )
            else:
                descs.append(
                    f"partial {a}x{b} [{actual}/{complete} edges, "
                    f"density={density:.2f}]"
                )

    # Build summary string
    # Count repeated K_{a,b} patterns
    from collections import Counter
    counts = Counter(descs)
    parts = []
    for desc, cnt in sorted(counts.items(),
                            key=lambda kv: -kv[1]):
        parts.append(f"{cnt}×{desc}" if cnt > 1 else desc)
    summary = " ⊔ ".join(parts)   # disjoint union symbol

    return {
        "n_components": nx.number_connected_components(G),
        "components":   descs,
        "summary":      summary,
    }


# ══════════════════════════════════════════════════════════════════
#  可视化 (单案例)
# ══════════════════════════════════════════════════════════════════

_CX0, _CX1 = "#C8D8FF", "#1A3A8C"
_CY0, _CY1 = "#FFD8C8", "#8C2A1A"
_CE        = "#546E7A"

def _lerp(t, c0, c1):
    r0,g0,b0 = int(c0[1:3],16),int(c0[3:5],16),int(c0[5:7],16)
    r1,g1,b1 = int(c1[1:3],16),int(c1[3:5],16),int(c1[5:7],16)
    return "#{:02x}{:02x}{:02x}".format(
        int(r0+t*(r1-r0)), int(g0+t*(g1-g0)), int(b0+t*(b1-b0)))

def _layout(n, m, x_col=0.0, y_col=3.0):
    half = (max(n, m) - 1) / 2
    pos = {}
    for i in range(n):
        pos[f"x{i}"] = (x_col, half-i*(2*half/max(n-1,1)) if n>1 else 0.0)
    for j in range(m):
        pos[f"y{j}"] = (y_col, half-j*(2*half/max(m-1,1)) if m>1 else 0.0)
    return pos

def visualize_single(n, m, max_edges, matrix, H, struct, verified,
                     save_path):
    """绘制带结构注释的单个极值图。"""
    p_orig, q_orig = len(H), len(H[0])
    H_edges_orig, _, _ = _get_edges(H)

    G = nx.Graph()
    for i in range(n): G.add_node(f"x{i}")
    for j in range(m): G.add_node(f"y{j}")
    for i in range(n):
        for j in range(m):
            if matrix[i][j]: G.add_edge(f"x{i}", f"y{j}")

    pos    = _layout(n, m)
    x_deg  = [sum(matrix[i])                      for i in range(n)]
    y_deg  = [sum(matrix[i][j] for i in range(n)) for j in range(m)]
    mx     = max(x_deg+[1]); my = max(y_deg+[1])
    x_col  = [_lerp(x_deg[i]/mx, _CX0, _CX1) for i in range(n)]
    y_col  = [_lerp(y_deg[j]/my, _CY0, _CY1) for j in range(m)]

    fig_h = max(7.0, max(n, m)*1.0)
    fig = plt.figure(figsize=(15, fig_h), facecolor="#F0F2F5")
    gs  = fig.add_gridspec(1, 2, width_ratios=[3.2, 1],
                           left=0.04, right=0.97,
                           top=0.85, bottom=0.16, wspace=0.05)
    ax  = fig.add_subplot(gs[0])
    axH = fig.add_subplot(gs[1])
    ax.set_facecolor("#FAFBFF")
    axH.set_facecolor("#FFFDE7")
    for sp in ax.spines.values():  sp.set_visible(False)
    for sp in axH.spines.values(): sp.set_color("#E0E0E0")

    # Component backgrounds
    _bg = ["#EEF0FF","#FFF4EC","#EDFFF0","#FFEDF5","#EDFFFC","#F5EDFF"]
    for ci, comp in enumerate(nx.connected_components(G)):
        xs_ = [pos[v][0] for v in comp]; ys_ = [pos[v][1] for v in comp]
        pad = 0.42
        ax.add_patch(mpatches.FancyBboxPatch(
            (min(xs_)-pad, min(ys_)-pad),
            max(xs_)-min(xs_)+2*pad, max(ys_)-min(ys_)+2*pad,
            boxstyle="round,pad=0.12", lw=1.3, ls="--",
            ec="#AAAACC", fc=_bg[ci%len(_bg)], alpha=0.5, zorder=0))
        ax.text((min(xs_)+max(xs_))/2, max(ys_)+pad+0.18,
                f"C{ci+1}", fontsize=8, color="#888", ha="center",
                va="bottom", zorder=1)

    nx.draw_networkx_edges(G, pos, edge_color=_CE,
                           width=2.2, alpha=0.55, ax=ax)
    nx.draw_networkx_nodes(G, pos,
        nodelist=[f"x{i}" for i in range(n)],
        node_color=x_col, node_size=820, edgecolors=_CX1,
        linewidths=2.4, ax=ax)
    nx.draw_networkx_nodes(G, pos,
        nodelist=[f"y{j}" for j in range(m)],
        node_color=y_col, node_size=820, edgecolors=_CY1,
        linewidths=2.4, ax=ax)

    labels = {f"x{i}": f"$x_{i}$\nd={x_deg[i]}" for i in range(n)}
    labels.update({f"y{j}": f"$y_{j}$\nd={y_deg[j]}" for j in range(m)})
    nx.draw_networkx_labels(G, pos, labels=labels, font_size=9,
                            font_weight="bold", font_color="#111", ax=ax)

    half = (max(n, m)-1)/2
    ax.text(-0.50, half, "X", fontsize=13, color=_CX1,
            fontweight="bold", va="center", ha="center")
    ax.text( 3.50, half, "Y", fontsize=13, color=_CY1,
            fontweight="bold", va="center", ha="center")

    """绘制带结构注释的单个极值图。"""
    n_comp = struct["n_components"]
    struct_lines = "\n".join(
        f"  C{i+1}: {d}" for i, d in enumerate(struct["components"])
    )
    info = (f"deg X: [{', '.join(str(d) for d in x_deg)}]\n"
            f"deg Y: [{', '.join(str(d) for d in y_deg)}]\n"
            f"Components ({n_comp}): {struct['summary']}\n"
            f"{struct_lines}")
    ax.text(0.5, -0.10, info,
            transform=ax.transAxes, ha="center", fontsize=8.5, color="#333",
            family="monospace",
            bbox=dict(boxstyle="round,pad=0.45", fc="white",
                      ec="#CCCCCC", alpha=0.92))

    vmark  = "H-free: verified" if verified else "WARNING: H found!"
    vcolor = "#1B5E20"          if verified else "#B71C1C"
    ax.text(0.02, 0.98, vmark,
            transform=ax.transAxes, ha="left", va="top",
            fontsize=9.5, color=vcolor, fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.3", fc="white",
                      ec=vcolor, alpha=0.85))

    ax.set_title(
        f"Extremal graph   ex({n}, {m}; H) = {max_edges}",
        fontsize=14, fontweight="bold", color="#1A1A2E", pad=8)
    ax.axis("off")

    # 侧面板：禁图 H
    GH = nx.Graph()
    for i in range(p_orig): GH.add_node(f"u{i}")
    for j in range(q_orig): GH.add_node(f"v{j}")
    for u, v in H_edges_orig: GH.add_edge(f"u{u}", f"v{v}")
    posH_raw = _layout(p_orig, q_orig)
    posH = {f"u{i}": posH_raw[f"x{i}"] for i in range(p_orig)}
    posH.update({f"v{j}": posH_raw[f"y{j}"] for j in range(q_orig)})
    hd_u = [sum(H[i]) for i in range(p_orig)]
    hd_v = [sum(H[i][j] for i in range(p_orig)) for j in range(q_orig)]
    nx.draw_networkx_edges(GH, posH, edge_color="#78909C",
                           width=2.0, alpha=0.8, ax=axH)
    nx.draw_networkx_nodes(GH, posH,
        nodelist=[f"u{i}" for i in range(p_orig)],
        node_color=[_lerp(hd_u[i]/max(hd_u+[1]), _CX0, _CX1)
                    for i in range(p_orig)],
        node_size=560, edgecolors=_CX1, linewidths=2, ax=axH)
    nx.draw_networkx_nodes(GH, posH,
        nodelist=[f"v{j}" for j in range(q_orig)],
        node_color=[_lerp(hd_v[j]/max(hd_v+[1]), _CY0, _CY1)
                    for j in range(q_orig)],
        node_size=560, edgecolors=_CY1, linewidths=2, ax=axH)
    H_lbl = {f"u{i}": f"u{i}(d={hd_u[i]})" for i in range(p_orig)}
    H_lbl.update({f"v{j}": f"v{j}(d={hd_v[j]})" for j in range(q_orig)})
    nx.draw_networkx_labels(GH, posH, labels=H_lbl, font_size=8,
                            font_weight="bold", ax=axH)
    half_H = (max(p_orig, q_orig)-1)/2
    axH.text(-0.22, half_H, "L", fontsize=10, color=_CX1,
             fontweight="bold", va="center", transform=axH.transData)
    axH.text( 3.22, half_H, "R", fontsize=10, color=_CY1,
             fontweight="bold", va="center", transform=axH.transData)
    axH.set_title(f"Forbidden H\n({p_orig}+{q_orig} vertices, "
                  f"{len(H_edges_orig)} edges)",
                  fontsize=9.5, fontweight="bold", color="#37474F", pad=8)
    axH.axis("off")

    legend_elems = [
        mpatches.Patch(fc=_CX0, ec=_CX1, lw=1.5, label="X (low deg)"),
        mpatches.Patch(fc=_CX1, ec=_CX1, lw=1.5, label="X (high deg)"),
        mpatches.Patch(fc=_CY0, ec=_CY1, lw=1.5, label="Y (low deg)"),
        mpatches.Patch(fc=_CY1, ec=_CY1, lw=1.5, label="Y (high deg)"),
    ]
    fig.legend(handles=legend_elems, loc="lower center", ncol=4,
               fontsize=9, frameon=True, fancybox=True,
               edgecolor="#CCCCCC", bbox_to_anchor=(0.44, 0.01))
    fig.suptitle(f"Bipartite Turan number  n={n}, m={m}",
                 fontsize=13, color="#333", y=0.98)

    plt.savefig(save_path, dpi=150, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    plt.close(fig)


# ══════════════════════════════════════════════════════════════════
#  摘要表格图片
# ══════════════════════════════════════════════════════════════════

def save_summary_table(results, H, save_path):
    """
    Save a text-table image summarising all computed cases.
    results: list of dicts with keys n, m, ex, verified, struct, elapsed
    """
    p, q = len(H), len(H[0])
    H_edges, _, _ = _get_edges(H)

    rows = [["n", "m", "ex(n,m;H)", "Verified", "Components", "Structure", "Time(s)"]]
    for r in results:
        rows.append([
            str(r["n"]),
            str(r["m"]),
            str(r["ex"]) if r["ex"] >= 0 else "FAIL",
            "YES" if r["verified"] else ("NO" if r["ex"] >= 0 else "-"),
            str(r["struct"]["n_components"]) if r["ex"] >= 0 else "-",
            r["struct"]["summary"] if r["ex"] >= 0 else r.get("error",""),
            f"{r['elapsed']:.1f}",
        ])

    col_widths = [max(len(rows[ri][ci]) for ri in range(len(rows)))
                  for ci in range(len(rows[0]))]

    # 根据行数缩放图高
    fig_h = max(3.0, 0.38 * len(rows) + 1.2)
    fig, ax = plt.subplots(figsize=(max(12, sum(col_widths)*0.13), fig_h),
                           facecolor="#F8F9FA")
    ax.set_facecolor("#F8F9FA")
    ax.axis("off")

    # 使用 matplotlib 表格绘制
    cell_text = [r for r in rows[1:]]
    col_labels = rows[0]

    tbl = ax.table(
        cellText=cell_text,
        colLabels=col_labels,
        loc="center",
        cellLoc="left",
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9.5)
    tbl.auto_set_column_width(list(range(len(col_labels))))

    # 表头样式
    for j in range(len(col_labels)):
        tbl[0, j].set_facecolor("#1A3A8C")
        tbl[0, j].set_text_props(color="white", fontweight="bold")

    # 数据行样式（交替颜色，验证列颜色编码）
    for i, r in enumerate(results):
        row_bg = "#FFFFFF" if i % 2 == 0 else "#F0F4FF"
        for j in range(len(col_labels)):
            tbl[i+1, j].set_facecolor(row_bg)
        # 验证列
        if r["ex"] >= 0:
            vcolor = "#E8F5E9" if r["verified"] else "#FFEBEE"
            tbl[i+1, 3].set_facecolor(vcolor)
        # 结构列：含 K_{a,b} 模式时绿色高亮
        if r["ex"] >= 0 and "K_" in r["struct"]["summary"]:
            tbl[i+1, 5].set_facecolor("#E8F5E9")

    ax.set_title(
        f"Batch results: ex(n, m; H)   H has {p}+{q} vertices, "
        f"{len(H_edges)} edges",
        fontsize=12, fontweight="bold", color="#1A1A2E",
        pad=14, loc="center",
    )

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"Summary table saved: {save_path}")


# ══════════════════════════════════════════════════════════════════
#  批量计算
# ══════════════════════════════════════════════════════════════════

def batch_compute(H, nm_pairs, out_dir="batch_output", time_limit=120):
    """
    为 nm_pairs 中的每个 (n,m) 计算 ex(n,m;H)。

    参数
    ----------
    H          : 禁图的邻接矩阵
    nm_pairs   : (n, m) 元组列表（n <= m）
    out_dir    : 保存 PNG 文件的目录
    time_limit : 每个案例 ILP 时间限制（秒）

    返回
    -------
    results : 结果 dict 列表
    """
    os.makedirs(out_dir, exist_ok=True)
    results = []
    total = len(nm_pairs)

    print(f"\nBatch: {total} cases, output -> '{out_dir}/'")
    print("=" * 72)
    print(f"{'n':>4} {'m':>4} {'ex':>6}  {'ok':>4}  {'comp':>5}  structure")
    print("-" * 72)

    for idx, (n, m) in enumerate(nm_pairs):
        t0 = time.time()
        result = {"n": n, "m": m, "ex": -1, "verified": False,
                  "struct": {"n_components": 0, "components": [],
                             "summary": ""}, "elapsed": 0.0}
        try:
            ex, matrix  = compute_ex(n, m, H, time_limit=time_limit)
            verified     = verify_H_free(matrix, H)
            struct       = analyze_structure(matrix)
            elapsed      = time.time() - t0

            result.update({"ex": ex, "verified": verified,
                           "struct": struct, "elapsed": elapsed})

            # 保存单个可视化图
            img_path = os.path.join(out_dir, f"ex_{n}x{m}.png")
            visualize_single(n, m, ex, matrix, H, struct, verified, img_path)

            vmark = "YES" if verified else "NO!"
            print(f"{n:>4} {m:>4} {ex:>6}  {vmark:>4}  "
                  f"{struct['n_components']:>5}  {struct['summary']}")

        except Exception as e:
            elapsed = time.time() - t0
            result["elapsed"] = elapsed
            result["error"]   = str(e)
            print(f"{n:>4} {m:>4}  {'FAIL':>6}  {'?':>4}  {'?':>5}  {e}")

        results.append(result)

    print("=" * 72)
    print(f"Done. {sum(1 for r in results if r['ex']>=0)}/{total} solved.\n")

    # 保存摘要表格图片
    save_summary_table(results, H,
                       os.path.join(out_dir, "_summary_table.png"))
    return results


# ══════════════════════════════════════════════════════════════════
#  输入辅助函数
# ══════════════════════════════════════════════════════════════════

def parse_nm_input(s):
    """
    解析用户输入的 n,m 规范。接受（请输入引号内的格式，不要带引号）：
      "3 8"          -> [(3,8)]
      "3-6 5-8"      -> 所有 3<=n<=6, 5<=m<=8 且 n<=m 的对
      "3,4,5 6,7,8"  -> 列出的值的笛卡尔积，仅保留 n<=m
      "3-6 *"        -> n 在 3..6, m = n..n+4（自动上三角）
    """
    parts = s.strip().split()
    if len(parts) != 2:
        raise ValueError("需输入两个空格分隔的规范")

    def parse_one(spec, default_max=10):
        if spec == "*":
            return None   # 特殊处理
        if "-" in spec:
            lo, hi = map(int, spec.split("-"))
            return list(range(lo, hi+1))
        elif "," in spec:
            return list(map(int, spec.split(",")))
        else:
            v = int(spec)
            return [v]

    n_vals = parse_one(parts[0])
    m_vals = parse_one(parts[1])

    pairs = []
    if m_vals is None:
        # "n_spec *": 每个 n，m 从 n 到 n+delta
        delta = 4
        for n in n_vals:
            for m in range(n, n+delta+1):
                pairs.append((n, m))
    else:
        for n in n_vals:
            for m in m_vals:
                if n <= m:
                    pairs.append((n, m))

    if not pairs:
        raise ValueError("No valid (n<=m) pairs found")
    return sorted(set(pairs))


# ══════════════════════════════════════════════════════════════════
#  主程序
# ══════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("=" * 60)
    print("  Bipartite Turan Calculator — BATCH MODE")
    print("  (both bipartition orientations checked)")
    print("=" * 60)

    # --- Input H ---
    p, q = map(int, input("\n请输入禁图 H 的行数 p 和列数 q:  ").split())
    print(f"请输入 H 邻接矩阵（{p} 行，每行 {q} 个值）：") 
    H = [list(map(int, input().split())) for _ in range(p)]

    H_edges, _, _ = _get_edges(H)
    print(f"H: {p}+{q} 个顶点, {len(H_edges)} 条边。检查的定向: "
          f"{len(_orientations(H, max(p,q)+1, max(p,q)+1))}")

    # --- Input (n,m) range ---
    print("\n指定 (n, m) 范围。示例：")
    print("  单案例     : '4 6'")
    print("  范围       : '3-6 4-8'   （所有 n<=m 的对）")
    print("  列表       : '3,4,5 5,6,7'")
    print("  自动三角   : '3-7 *'     （每个 n，m 从 n 到 n+4）")
    nm_spec = input("请输入 n m 规范（请输入引号内的格式，不要带引号）: ").strip()
    nm_pairs = parse_nm_input(nm_spec)
    print(f"-> {len(nm_pairs)} cases: {nm_pairs[:8]}"
          + (" ..." if len(nm_pairs) > 8 else ""))

    # --- Output directory ---
    out_dir = input("\n输出目录 [batch_output]:  ").strip() or "batch_output"

    # --- Time limit ---
    tl_str = input("每个案例 ILP 求解时间限制（秒）[60]:  ").strip()
    time_limit = int(tl_str) if tl_str else 60

    # --- Run ---
    batch_compute(H, nm_pairs, out_dir=out_dir, time_limit=time_limit)
    print(f"\n所有图片已保存至：{out_dir}/")
    print("  单个图表 ：ex_NxM.png")
    print("  摘要表格 ：_summary_table.png")
