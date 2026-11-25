"""根据 SM3 公式验证完整线路的布尔接线。

为使 32 位模拟可运行，将每个 QAND/QAND† 临时替换为等价的逻辑 Toffoli 门。
此验证只检查位连接和可逆清理，不用于统计 T 门数。
"""

import argparse
from random import Random

from projectq import MainEngine
from projectq.backends import ClassicalSimulator
from projectq.ops import All, Measure, Toffoli

import optimal_depth_adder as core
import sm3_gate_blocks as blocks
from linear_permutations import rotl32
from sm3_full_gate_level import IV, allocate_block_inputs, compress_block, set_constant


MASK = 0xFFFFFFFF


def _logical_and(eng, a, b, target):
    Toffoli | (a, b, target)


core.quantum_and = _logical_and
core.quantum_and_dag = _logical_and
blocks.quantum_and = _logical_and
blocks.quantum_and_dag_figure = _logical_and


def sm3_reference(block_words, rounds):
    w = block_words[:]
    for j in range(16, 68):
        u = w[j - 16] ^ w[j - 9] ^ rotl32(w[j - 3], 15)
        w.append((u ^ rotl32(u, 15) ^ rotl32(u, 23)
                  ^ rotl32(w[j - 13], 7) ^ w[j - 6]) & MASK)
    a, b, c, d, e, f, g, h = IV
    for j in range(rounds):
        constant = 0x79CC4519 if j < 16 else 0x7A879D8A
        ss1 = rotl32((rotl32(a, 12) + e + rotl32(constant, j % 32)) & MASK, 7)
        ss2 = ss1 ^ rotl32(a, 12)
        ff = a ^ b ^ c if j < 16 else (a & b) | (a & c) | (b & c)
        gg = e ^ f ^ g if j < 16 else (e & f) | ((~e & MASK) & g)
        tt1 = (ff + d + ss2 + (w[j] ^ w[j + 4])) & MASK
        tt2 = (gg + h + ss1 + w[j]) & MASK
        p0 = tt2 ^ rotl32(tt2, 9) ^ rotl32(tt2, 17)
        a, b, c, d, e, f, g, h = (
            tt1, a, rotl32(b, 9), c, p0, e, rotl32(f, 19), g
        )
    return (a, b, c, d, e, f, g, h)


def read_word(word):
    return sum(int(qubit) << bit for bit, qubit in enumerate(word))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rounds", type=int, default=1)
    args = parser.parse_args()
    rng = Random(20261123)
    block = [rng.getrandbits(32) for _ in range(16)]
    eng = MainEngine(ClassicalSimulator())
    window, lines = allocate_block_inputs(eng)
    for word, value in zip(window, block):
        set_constant(word, value)
    for line in lines:
        for word, value in zip(line.as_tuple(), IV):
            set_constant(word, value)
    compress_block(eng, window, lines, args.rounds)
    for line in lines:
        for word in line.as_tuple():
            All(Measure) | word
    eng.flush()
    expected = sm3_reference(block, args.rounds)
    for number, line in enumerate(lines, 1):
        actual = tuple(read_word(word) for word in line.as_tuple())
        if actual != expected:
            raise AssertionError((number, actual, expected))
    print(f"Both state lines matched the SM3 equations for {args.rounds} rounds.")


if __name__ == "__main__":
    main()
