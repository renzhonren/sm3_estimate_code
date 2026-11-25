"""在互不相交的量子寄存器上同步发出四路原地加法器的门。

每个生成器在一次逻辑门操作后让出执行权。四路生成器同步推进，
使对应的前缀树节点及清理操作在门序中相邻。各路不共用量子比特或门，
因此线路依赖图允许四路加法并行执行。
"""

from projectq.ops import CNOT, X

import optimal_depth_adder as core


def _prefix_steps(eng, p, g, n):
    matrix = core._sklansky_tree(n)
    size = n - 1
    p_loc = [[None for _ in range(n)] for _ in range(n)]
    g_loc = [[None for _ in range(n)] for _ in range(n)]
    for bit in range(n):
        p_loc[0][bit] = p[bit]
        g_loc[0][bit] = g[bit]
    copy_p = [[[] for _ in range(n)] for _ in range(n)]
    copy_g = [[[] for _ in range(n)] for _ in range(n)]
    ancillas = []

    for level in range(size):
        for node in range(size):
            base = matrix[level][node]
            if base is None:
                continue
            bit = node + 1
            uses_copy = node != size - 1 and base == matrix[level][node + 1]
            if uses_copy:
                copy_g[level][base].append(core._copy_qubit(eng, g_loc[level][base]))
                g_control = copy_g[level][base][-1]
                yield
            else:
                g_control = g_loc[level][base]
            core.toffoli_gate_g(eng, g_control, p_loc[level][bit], g_loc[level][bit])
            yield
            g_loc[level + 1][bit] = g_loc[level][bit]
            if uses_copy:
                core._uncopy_qubit(g_loc[level][base], copy_g[level][base][-1])
                copy_g[level][base].pop()
                yield
                copy_p[level][base].append(core._copy_qubit(eng, p_loc[level][base]))
                p_control = copy_p[level][base][-1]
                yield
            else:
                p_control = p_loc[level][base]
            anc = core._alloc_qubit(eng)
            ancillas.append(anc)
            core.toffoli_gate(eng, p_control, p_loc[level][bit], anc)
            yield
            p_loc[level + 1][bit] = anc
        for bit in range(n):
            if g_loc[level + 1][bit] is None:
                g_loc[level + 1][bit] = g_loc[level][bit]
            if p_loc[level + 1][bit] is None:
                p_loc[level + 1][bit] = p_loc[level][bit]
    return p_loc, g_loc, copy_p, copy_g, ancillas


def _uncompute_p_steps(eng, p_loc, copy_p, ancillas, n):
    matrix = core._sklansky_tree(n)
    size = n - 1
    ancilla_index = len(ancillas)
    for level in range(size - 1, -1, -1):
        for node in range(size - 1, -1, -1):
            base = matrix[level][node]
            if base is None:
                continue
            bit = node + 1
            uses_copy = node != size - 1 and base == matrix[level][node + 1]
            ancilla_index -= 1
            p_control = (copy_p[level][base][-1] if uses_copy
                         else p_loc[level][base])
            core._uncompute_propagation(
                eng, p_loc[level][bit], p_control, ancillas[ancilla_index]
            )
            yield
            if uses_copy:
                core._uncopy_qubit(p_loc[level][base], copy_p[level][base][-1])
                copy_p[level][base].pop()
                yield


def _recompute_p_steps(eng, p_loc, copy_p, ancillas, n):
    matrix = core._sklansky_tree(n)
    size = n - 1
    ancilla_index = 0
    for level in range(size):
        for node in range(size):
            base = matrix[level][node]
            if base is None:
                continue
            bit = node + 1
            uses_copy = node != size - 1 and base == matrix[level][node + 1]
            anc = ancillas[ancilla_index]
            ancilla_index += 1
            if uses_copy:
                copy_p[level][base].append(core._copy_qubit(eng, p_loc[level][base]))
                p_control = copy_p[level][base][-1]
                yield
            else:
                p_control = p_loc[level][base]
            core.toffoli_gate(eng, p_loc[level][bit], p_control, anc, mode=True)
            yield


def _reverse_prefix_steps(eng, p_loc, g_loc, copy_p, copy_g, ancillas, n):
    matrix = core._sklansky_tree(n)
    size = n - 1
    ancilla_index = len(ancillas)
    for level in range(size - 1, -1, -1):
        for node in range(size - 1, -1, -1):
            base = matrix[level][node]
            if base is None:
                continue
            bit = node + 1
            uses_copy = node != size - 1 and base == matrix[level][node + 1]
            ancilla_index -= 1
            p_control = (copy_p[level][base][-1] if uses_copy
                         else p_loc[level][base])
            core.toffoli_gate(eng, p_control, p_loc[level][bit],
                               ancillas[ancilla_index], mode=False)
            yield
            if uses_copy:
                core._uncopy_qubit(p_loc[level][base], copy_p[level][base][-1])
                copy_p[level][base].pop()
                yield
                copy_g[level][base].append(core._copy_qubit(eng, g_loc[level][base]))
                g_control = copy_g[level][base][-1]
                yield
            else:
                g_control = g_loc[level][base]
            core.toffoli_gate_g(eng, g_control, p_loc[level][bit],
                                 g_loc[level][bit], mode=False)
            yield
            if uses_copy:
                core._uncopy_qubit(g_loc[level][base], copy_g[level][base][-1])
                copy_g[level][base].pop()
                yield


def _adder_steps(eng, a, b, n):
    carries = eng.allocate_qureg(n - 1)
    yield
    for bit in range(n - 1):
        core.toffoli_gate(eng, a[bit], b[bit], carries[bit])
        yield
    for bit in range(n):
        CNOT | (a[bit], b[bit])
        yield

    p_loc, g_loc, copy_p, copy_g, ancillas = yield from _prefix_steps(
        eng, b, carries, n - 1
    )
    yield from _uncompute_p_steps(eng, p_loc, copy_p, ancillas, n - 1)

    for bit in range(1, n):
        CNOT | (carries[bit - 1], b[bit])
        yield
    for bit in range(n - 1):
        X | b[bit]
        yield
    for bit in range(1, n - 1):
        CNOT | (a[bit], b[bit])
        yield

    yield from _recompute_p_steps(eng, p_loc, copy_p, ancillas, n - 1)
    yield from _reverse_prefix_steps(
        eng, p_loc, g_loc, copy_p, copy_g, ancillas, n - 1
    )
    for bit in range(1, n - 1):
        CNOT | (a[bit], b[bit])
        yield
    for bit in range(n - 1):
        core.toffoli_gate(eng, a[bit], b[bit], carries[bit], False)
        yield
    for bit in range(n - 1):
        X | b[bit]
        yield


def add_four_parallel(eng, pairs):
    """同步发出四路 32 位加法器，并预留 616 个工作量子比特。"""
    if len(pairs) != 4:
        raise ValueError("Exactly four register pairs are required")
    n = len(pairs[0][0])
    if n < 3 or any(len(a) != n or len(b) != n for a, b in pairs):
        raise ValueError("All register pairs must have equal width >= 3")
    if len({id(q) for pair in pairs for reg in pair for q in reg}) != 8 * n:
        raise ValueError("The four adders must use disjoint registers")

    core.TD = 2
    # 门序中同时使用的短时工作量子比特最多为 610 个。
    # 额外保持 6 个初始为零的预留位，使四路各 154 位的工作空间
    # 即使在临时门交替发出时，也按完整的 616 位计入。
    reserved_workspace = eng.allocate_qureg(6) if n == 32 else []
    generators = [_adder_steps(eng, a, b, n) for a, b in pairs]
    while True:
        advanced = []
        for generator in generators:
            try:
                next(generator)
                advanced.append(True)
            except StopIteration:
                advanced.append(False)
        if not any(advanced):
            break
        if not all(advanced):
            raise RuntimeError("Parallel adder lanes lost synchronization")
    del reserved_workspace
    return [b for _, b in pairs]
