"""论文中 SM3 压缩线路的单比特门级部件。

各函数分别提供正算与反算清理阶段。与运算采用原地加法器代码中的分解。
"""

from projectq.ops import CNOT, H, Measure, X

from optimal_depth_adder import quantum_and


def quantum_and_dag_figure(eng, a, b, target):
    """图 3(b) 的测量式反算，并按测量结果复位目标位。"""
    H | target
    H | b
    Measure | target
    if int(target):
        CNOT | (a, b)
        X | target
    H | b


def qfa_forward(eng, a, b, c, carry):
    """单比特 QFA；carry 为 None 时省略无用的最高位进位。"""
    CNOT | (a, b)
    CNOT | (a, c)
    if carry is not None:
        quantum_and(eng, b, c, carry)
    CNOT | (a, b)
    if carry is not None:
        CNOT | (a, carry)
    CNOT | (b, c)


def qfa_cleanup(eng, a, b, c, carry):
    """反算 QFA；若曾计算进位，也清理相应进位位。"""
    CNOT | (b, c)
    if carry is not None:
        CNOT | (a, carry)
    CNOT | (a, b)
    if carry is not None:
        quantum_and_dag_figure(eng, b, c, carry)
    CNOT | (a, c)
    CNOT | (a, b)


def ff_forward(eng, round_index, x, y, z, out):
    """将 FF 结果写入初始为零的 out；后 48 轮会临时修改 y 和 z。"""
    if round_index < 16:
        CNOT | (x, out)
        CNOT | (y, out)
        CNOT | (z, out)
    else:
        # 多数函数：majority(x,y,z) = x XOR ((y XOR x) AND (z XOR x))。
        CNOT | (x, y)
        CNOT | (x, z)
        quantum_and(eng, y, z, out)
        CNOT | (x, out)


def ff_cleanup(eng, round_index, x, y, z, out):
    if round_index < 16:
        CNOT | (z, out)
        CNOT | (y, out)
        CNOT | (x, out)
    else:
        CNOT | (x, out)
        quantum_and_dag_figure(eng, y, z, out)
        CNOT | (x, z)
        CNOT | (x, y)


def gg_forward(eng, round_index, x, y, z, out):
    """将 GG 结果写入初始为零的 out；后 48 轮会临时修改 y。"""
    if round_index < 16:
        CNOT | (x, out)
        CNOT | (y, out)
        CNOT | (z, out)
    else:
        # 选择函数：choose(x,y,z) = z XOR (x AND (y XOR z))。
        CNOT | (z, y)
        quantum_and(eng, x, y, out)
        CNOT | (z, out)


def gg_cleanup(eng, round_index, x, y, z, out):
    if round_index < 16:
        CNOT | (z, out)
        CNOT | (y, out)
        CNOT | (x, out)
    else:
        CNOT | (z, out)
        quantum_and_dag_figure(eng, x, y, out)
        CNOT | (z, y)


def qfa_words_forward(eng, a, b, c, carries):
    """32 位模加 QFA；不计算最高位向外输出的进位。"""
    if any(len(word) != 32 for word in (a, b, c, carries)):
        raise ValueError("Expected 32-bit words")
    for bit in range(32):
        target = carries[bit + 1] if bit < 31 else None
        qfa_forward(eng, a[bit], b[bit], c[bit], target)


def qfa_words_cleanup(eng, a, b, c, carries):
    for bit in range(31, -1, -1):
        target = carries[bit + 1] if bit < 31 else None
        qfa_cleanup(eng, a[bit], b[bit], c[bit], target)
