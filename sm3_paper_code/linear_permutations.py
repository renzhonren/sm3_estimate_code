"""用户提供的 SM3 P0、P1 原地线路，各包含 68 个 CNOT 门。

每对索引依次表示控制位和目标位。字采用小端位序，位 0 为最低位。
门序保持用户提供时的顺序。
"""

from projectq.ops import CNOT


P1_CNOTS = (
    (1,25), (3,11), (27,19), (9,17), (18,2), (1,9), (18,1), (27,18),
    (4,27), (13,4), (21,4), (22,13), (3,18), (20,3), (5,20), (22,5),
    (29,20), (21,29), (6,21), (15,6), (24,15), (31,22), (7,22), (23,6),
    (24,7), (23,31), (8,23), (9,24), (9,8), (26,9), (11,26), (1,9),
    (10,1), (19,10), (11,19), (12,11), (12,3), (12,27), (28,11), (29,12),
    (21,29), (30,21), (14,29), (13,28), (30,13), (5,28), (14,5), (15,30),
    (31,14), (23,31), (7,30), (0,23), (0,15), (17,0), (16,7), (16,31),
    (25,16), (25,17), (2,25), (18,2), (19,2), (11,19), (27,19), (3,11),
    (17,8), (25,17), (1,25), (9,17),
)

P0_CNOTS = (
    (0,8), (24,16), (18,10), (9,25), (26,2), (18,26), (9,18), (0,9),
    (23,0), (24,9), (7,24), (14,23), (5,14), (22,7), (6,23), (30,7),
    (5,22), (21,30), (20,5), (3,20), (21,6), (4,21), (19,4), (28,5),
    (12,21), (3,12), (26,3), (19,28), (26,19), (17,26), (8,17), (15,24),
    (18,26), (1,18), (16,1), (8,16), (15,8), (31,8), (15,0), (30,15),
    (13,30), (22,31), (13,22), (4,13), (28,13), (6,15), (14,31), (29,14),
    (29,6), (12,29), (27,4), (27,12), (10,27), (20,29), (11,20), (11,28),
    (2,11), (2,19), (10,19), (25,10), (25,2), (16,25), (24,16), (8,16),
    (26,2), (18,10), (9,25), (0,8),
)

assert len(P0_CNOTS) == len(P1_CNOTS) == 68


def apply_linear(word, network, inverse=False):
    """将 CNOT 网络原地作用于一个 32 量子比特字。"""
    if len(word) != 32:
        raise ValueError("P0 and P1 require 32-bit words")
    for source, target in reversed(network) if inverse else network:
        CNOT | (word[source], word[target])


def p0(word, inverse=False):
    apply_linear(word, P0_CNOTS, inverse)


def p1(word, inverse=False):
    apply_linear(word, P1_CNOTS, inverse)


def apply_classical(value, network):
    """在 Python 整数上执行同一网络，以独立验证位序。"""
    for source, target in network:
        value ^= ((value >> source) & 1) << target
    return value


def rotl32(value, amount):
    amount %= 32
    return ((value << amount) | (value >> (32 - amount))) & 0xFFFFFFFF


def self_check():
    from random import Random

    rng = Random(20261123)
    for network, rotations in ((P0_CNOTS, (9, 17)), (P1_CNOTS, (15, 23))):
            # 32 个基向量足以完全确定这个 32 位线性变换。
        for value in [0, 0xFFFFFFFF, 0x12345678] + [1 << bit for bit in range(32)] + [rng.getrandbits(32) for _ in range(1000)]:
            expected = value ^ rotl32(value, rotations[0]) ^ rotl32(value, rotations[1])
            actual = apply_classical(value, network)
            if actual != expected:
                raise AssertionError((rotations, hex(value), hex(actual), hex(expected)))


if __name__ == "__main__":
    self_check()
    print("P0 and P1: 68 CNOTs each; all validation vectors matched.")
