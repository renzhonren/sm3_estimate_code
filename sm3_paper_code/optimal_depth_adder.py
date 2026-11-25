from __future__ import annotations

from projectq import MainEngine
from projectq.ops import H, CNOT, Measure, Toffoli, X, All, T, Tdag, S, Swap, Sdag
from projectq.backends import CircuitDrawer, ResourceCounter, CommandPrinter, ClassicalSimulator
from projectq.meta import Loop, Compute, Uncompute, Control, Dagger
import random



from projectq.ops import CNOT, H, Measure, S, Sdag, T, Tdag, Toffoli, X
# Toffoli 实现模式，与原资源统计代码保持一致：
#   TD = 0：直接使用 Toffoli 门，不进行分解。
#   TD = 1：使用无需辅助位的 Clifford+T 分解。
#   TD = 2：使用 quantum_and / quantum_and_dag 分解。
TD = 2


def quantum_and(eng, a, b, c):
    ancilla = eng.allocate_qubit()
    H | c
    CNOT | (b, ancilla)
    CNOT | (c, a)
    CNOT | (c, b)
    CNOT | (a, ancilla)
    Tdag | a
    Tdag | b
    T | c
    T | ancilla
    CNOT | (a, ancilla)
    CNOT | (c, b)
    CNOT | (c, a)
    CNOT | (b, ancilla)
    H | c
    S | c


def quantum_and_dag(eng, a, b, c):
    H | c
    Measure | c
    if int(c):
        S | a
        S | b
        X | c
        CNOT | (a, b)
        Sdag | b
        CNOT | (a, b)


def _toffoli_via_quantum_and_any_target(eng, a, b, c):
    """在 TD=2 模式下对任意初态的目标位 c 实现安全的 Toffoli。

    先计算 tmp=AND(a,b)，再执行 CNOT(tmp,c)，最后对 tmp 反算清理。
    """
    tmp = eng.allocate_qubit()[0]
    quantum_and(eng, a, b, tmp)
    CNOT | (tmp, c)
    quantum_and_dag(eng, a, b, tmp)


def toffoli_gate(eng, a, b, c, mode=True):
    """传播信号 P 路径使用的默认 Toffoli 实现。

    TD=2 时沿用逻辑与门形式；mode=True 要求目标位初始为零。
    """
    if TD == 1:
        Tdag | a
        Tdag | b
        H | c
        CNOT | (c, a)
        T | a
        CNOT | (b, c)
        CNOT | (b, a)
        T | c
        Tdag | a
        CNOT | (b, c)
        CNOT | (c, a)
        T | a
        Tdag | c
        CNOT | (b, a)
        H | c
    elif TD == 2:
        if mode:
            quantum_and(eng, a, b, c)
        else:
            quantum_and_dag(eng, a, b, c)
    else:
        Toffoli | (a, b, c)


def toffoli_gate_g(eng, a, b, c, mode=True):
    """生成信号 G 传播路径使用的 Toffoli 实现。

    TD=2 时先将与项写入临时位，再异或到目标位并清理临时位，
    因此不要求目标位初始为零。
    """
    if TD == 2:
        _toffoli_via_quantum_and_any_target(eng, a, b, c)
    else:
        toffoli_gate(eng, a, b, c, mode=mode)


def _alloc_qubit(eng):
    return eng.allocate_qureg(1)[0]


def _copy_qubit(eng, src):
    dst = _alloc_qubit(eng)
    CNOT | (src, dst)
    return dst


def _uncopy_qubit(src, dst):
    CNOT | (src, dst)


def _sklansky_tree(bitwidth: int):
    # 保持与参考实现相同的前缀树结构。
    pow2 = lambda x: 2**x
    log2 = lambda x: x.bit_length()

    bitwidth_i = bitwidth - 1
    matrix = [[None for _ in range(bitwidth_i)] for _ in range(bitwidth_i)]
    for row_index, row in enumerate(matrix):
        shift = pow2(row_index + 1)
        repeats = pow2(row_index)
        starting_index = repeats - 1
        last_index = bitwidth - 1
        for i in range(starting_index, last_index, shift):
            for j in range(repeats):
                try:
                    row[j + i] = i
                except IndexError:
                    pass
    return matrix

def _uncompute_propagation(eng, px, py, ancilla_p):
    toffoli_gate(eng, px, py, ancilla_p, mode=False)


def _step3_uncompute_p(
    eng,
    p_location,
    copy_p_location,
    ancilla_location,
    n: int,
):
    sklansky_matrix = _sklansky_tree(n)
    size = n - 1
    ancilla_index = len(ancilla_location)

    for i in range(size - 1, -1, -1):
        for j in range(size - 1, -1, -1):
            base = sklansky_matrix[i][j]
            if base is None:
                continue

            index_now = j + 1
            uses_copy = j != size - 1 and base == sklansky_matrix[i][j + 1]

            ancilla_index -= 1
            anc = ancilla_location[ancilla_index]

            if not uses_copy:
                _uncompute_propagation(eng, p_location[i][index_now], p_location[i][base], anc)
            else:
                _uncompute_propagation(eng, p_location[i][index_now], copy_p_location[i][base][-1], anc)
                _uncopy_qubit(p_location[i][base], copy_p_location[i][base][-1])
                copy_p_location[i][base].pop()

def _step3_recompute_p(
    eng,
    p_location,
    copy_p_location,
    ancilla_location,
    n: int,
):
    """精确撤销 _step3_uncompute_p 的操作。

    重新生成复制的 P 控制位和 P 辅助位。
    """
    sklansky_matrix = _sklansky_tree(n)
    size = n - 1
    ancilla_index = 0

    for i in range(size):
        for j in range(size):
            base = sklansky_matrix[i][j]
            if base is None:
                continue

            index_now = j + 1
            uses_copy = (j != size - 1 and base == sklansky_matrix[i][j + 1])

            anc = ancilla_location[ancilla_index]
            ancilla_index += 1

            if not uses_copy:
                # 撤销对应的 _uncompute_propagation 操作。
                toffoli_gate(
                    eng,
                    p_location[i][index_now],
                    p_location[i][base],
                    anc,
                    mode=True,
                )
            else:
                # 按反向顺序执行：先重新生成复制的控制位，
                copied_p = _copy_qubit(eng, p_location[i][base])
                copy_p_location[i][base].append(copied_p)
                # 再利用复制的控制位重新计算辅助位。
                toffoli_gate(
                    eng,
                    p_location[i][index_now],
                    copied_p,
                    anc,
                    mode=True,
                )


def _reverse_sklansky_prefix(
    eng,
    p_location,
    g_location,
    copy_p_location,
    copy_g_location,
    ancilla_location,
    n: int,
):
    """精确撤销 _sklansky_prefix，按节点的逆序遍历。

    正算节点依次可能复制 G、传播 G、复制 P、传播 P；
    反算则严格反向执行：撤销 P 传播、清除 P 副本、
    撤销 G 传播、清除 G 副本。
    """
    sklansky_matrix = _sklansky_tree(n)
    size = n - 1
    ancilla_index = len(ancilla_location)

    for i in range(size - 1, -1, -1):
        for j in range(size - 1, -1, -1):
            base = sklansky_matrix[i][j]
            if base is None:
                continue

            index_now = j + 1
            uses_copy = j != size - 1 and base == sklansky_matrix[i][j + 1]

            ancilla_index -= 1
            anc = ancilla_location[ancilla_index]

            # 撤销当前节点的 P 传播。
            if not uses_copy:
                p_ctrl = p_location[i][base]
            else:
                p_ctrl = copy_p_location[i][base][-1]
            toffoli_gate(
                eng,
                p_ctrl,
                p_location[i][index_now],
                anc,
                mode=False,
            )

            # 若曾复制 P 控制位，则清除该副本。
            if uses_copy:
                _uncopy_qubit(p_location[i][base], copy_p_location[i][base][-1])
                copy_p_location[i][base].pop()

            # 撤销当前节点的 G 传播。
            if not uses_copy:
                g_ctrl = g_location[i][base]
            else:
                # 反算时临时重新生成 G 副本，以对应正算节点的
                # 复制、使用和释放顺序。
                copied_g = _copy_qubit(eng, g_location[i][base])
                copy_g_location[i][base].append(copied_g)
                g_ctrl = copied_g
            toffoli_gate_g(
                eng,
                g_ctrl,
                p_location[i][index_now],
                g_location[i][index_now],
                mode=False,
            )

            # 若曾复制 G 控制位，则清除该副本。
            if uses_copy:
                _uncopy_qubit(g_location[i][base], copy_g_location[i][base][-1])
                copy_g_location[i][base].pop()







def _sklansky_prefix(eng, p, g, n: int):
    sklansky_matrix = _sklansky_tree(n)
    size = n - 1

    # 两个位置表分别按“层号、位号”索引 P 和 G 信号。
    p_location = [[None for _ in range(n)] for _ in range(n)]
    g_location = [[None for _ in range(n)] for _ in range(n)]
    for b in range(n):
        p_location[0][b] = p[b]
        g_location[0][b] = g[b]

    copy_p_location = [[[] for _ in range(n)] for _ in range(n)]
    copy_g_location = [[[] for _ in range(n)] for _ in range(n)]
    ancilla_location = []

    for i in range(size):
        for j in range(size):
            base = sklansky_matrix[i][j]
            if base is None:
                continue

            index_now = j + 1
            uses_copy = j != size - 1 and base == sklansky_matrix[i][j + 1]

            if uses_copy:
                copy_g_location[i][base].append(_copy_qubit(eng, g_location[i][base]))
                g_ctrl = copy_g_location[i][base][-1]
            else:
                g_ctrl = g_location[i][base]

            # 更新生成信号：g' = gy XOR (gx AND py)。
            toffoli_gate_g(
                eng,
                g_ctrl,
                p_location[i][index_now],
                g_location[i][index_now],
            )
            g_location[i + 1][index_now] = g_location[i][index_now]

            if uses_copy:
                _uncopy_qubit(g_location[i][base], copy_g_location[i][base][-1])
                copy_g_location[i][base].pop()

            if uses_copy:
                copy_p_location[i][base].append(_copy_qubit(eng, p_location[i][base]))
                p_ctrl = copy_p_location[i][base][-1]
            else:
                p_ctrl = p_location[i][base]

            anc_p = _alloc_qubit(eng)
            ancilla_location.append(anc_p)
            # 更新传播信号：p' = px AND py。
            toffoli_gate(eng, p_ctrl, p_location[i][index_now], anc_p)
            p_location[i + 1][index_now] = anc_p

        # 未更新的信号直接传递到下一层。
        for k in range(n):
            if g_location[i + 1][k] is None:
                g_location[i + 1][k] = g_location[i][k]
            if p_location[i + 1][k] is None:
                p_location[i + 1][k] = p_location[i][k]

    return p_location, g_location, copy_p_location, copy_g_location, ancilla_location


def _run_section3(eng, p, g, n: int):
    """执行论文中的第 3 阶段：计算前缀进位，再清理 P 辅助位。

    该阶段内部完成临时 G 副本的使用和清理，
    因此后续可按逆序撤销整个阶段。
    """
    p_location, g_location, copy_p_location, copy_g_location, ancilla_location = _sklansky_prefix(
        eng, p, g, n
    )
    _step3_uncompute_p(eng, p_location, copy_p_location, ancilla_location, n)
    return p_location, g_location, copy_p_location, copy_g_location, ancilla_location


def optimal_depth_add_in_place(eng, a, b, n: int | None = None, td: int | None = None):
    """原地模加法器：保持 a 不变，将 b 更新为 a+b（模 2^n）。"""
    global TD
    if td is not None:
        if td not in (0, 1, 2):
            raise ValueError("td must be 0, 1, or 2")
        TD = td

    if n is None:
        n = len(a)
    if n != len(a) or n != len(b):
        raise ValueError("Bit-width mismatch for optimal_depth_add_in_place")
    if n == 1:
        CNOT | (a[0], b[0])
        return b

    # 只分配进位 c_1 至 c_{n-1}，不分配占位寄存器。
    z = eng.allocate_qureg(n - 1)

    # 第一步：计算初始生成信号和传播信号。
    for i in range(n - 1):
        toffoli_gate(eng, a[i], b[i], z[i])
    for i in range(n):
        CNOT | (a[i], b[i])  # b 中此时保存传播信号 p。


    section3_state = _run_section3(eng, b, z, n - 1)

    for i in range(1, n):
        CNOT | (z[i - 1], b[i])
    for i in range(n - 1):
        X | b[i]
    for i in range(1, n - 1):
        CNOT | (a[i], b[i])


    p_location, g_location, copy_p_location, copy_g_location, ancilla_location = section3_state
    _step3_recompute_p(eng, p_location, copy_p_location, ancilla_location, n - 1)
    _reverse_sklansky_prefix(
        eng,
        p_location,
        g_location,
        copy_p_location,
        copy_g_location,
        ancilla_location,
        n - 1,
    )

    for i in range(1, n - 1):
        CNOT | (a[i], b[i])
    for i in range(n - 1):
        toffoli_gate(eng, a[i], b[i], z[i],False)
    for i in range(n - 1):
        X | b[i]



    return b
