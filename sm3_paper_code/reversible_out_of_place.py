"""带显式清理路径的 32 位异地加法器。

原始正算函数只返回输出字；压缩线路还需反算临时前缀树。
本模块沿用相同的 Sklansky 结构，并使用原地加法器提供的清理函数，
保存反向遍历所需的状态。
"""

from dataclasses import dataclass

from projectq.ops import CNOT

import optimal_depth_adder as core


@dataclass
class AdderState:
    a: object
    b: object
    output: object
    p_location: object
    g_location: object
    copy_p_location: object
    copy_g_location: object
    ancilla_location: object
    released_p_positions: object
    released_p_count: int


def _release_clean_p_ancillas(p_loc, ancillas):
    """释放已测量并复位的 P 工作位，直到反算需要时再分配。"""
    index_by_id = {id(qubit): index for index, qubit in enumerate(ancillas)}
    positions = []
    for level, row in enumerate(p_loc):
        for bit, qubit in enumerate(row):
            index = index_by_id.get(id(qubit))
            if index is not None:
                positions.append((level, bit, index))
                row[bit] = None
    count = len(ancillas)
    ancillas[:] = [None] * count
    return positions, count


def _restore_clean_p_ancillas(eng, state):
    """为反向前缀计算重新分配初始为零的工作量子比特。"""
    ancillas = [eng.allocate_qubit()[0] for _ in range(state.released_p_count)]
    for level, bit, index in state.released_p_positions:
        state.p_location[level][bit] = ancillas[index]
    state.ancilla_location[:] = ancillas


def compute(eng, a, b):
    n = len(a)
    if n != len(b) or n < 3:
        raise ValueError("Equal input widths of at least 3 are required")
    core.TD = 2
    output = eng.allocate_qureg(n)
    for bit in range(n - 1):
        core.toffoli_gate(eng, a[bit], b[bit], output[bit + 1])
    for bit in range(n):
        CNOT | (a[bit], b[bit])

    p_loc, g_loc, copy_p, copy_g, ancillas = core._run_section3(
        eng, b, output[1:], n - 1
    )

    for bit in range(n):
        CNOT | (b[bit], output[bit])
        CNOT | (a[bit], b[bit])
    positions, count = _release_clean_p_ancillas(p_loc, ancillas)
    return AdderState(a, b, output, p_loc, g_loc, copy_p, copy_g,
                      ancillas, positions, count)


def cleanup(eng, state):
    """将输出字和所有临时前缀辅助位恢复为零。"""
    a, b, output = state.a, state.b, state.output
    n = len(a)
    _restore_clean_p_ancillas(eng, state)
    for bit in range(n):
        CNOT | (a[bit], b[bit])
        CNOT | (b[bit], output[bit])

    core._step3_recompute_p(
        eng, state.p_location, state.copy_p_location, state.ancilla_location, n - 1
    )
    core._reverse_sklansky_prefix(
        eng,
        state.p_location,
        state.g_location,
        state.copy_p_location,
        state.copy_g_location,
        state.ancilla_location,
        n - 1,
    )

    for bit in range(n):
        CNOT | (a[bit], b[bit])
    for bit in range(n - 1):
        core.toffoli_gate(eng, a[bit], b[bit], output[bit + 1], mode=False)
