"""基于 16 字循环寄存器的 SM3 消息扩展门级实现。

生成 W16 至 W67 共 52 个消息字；每次执行四组 32 位异或和 68 个 CNOT 的 P1。
将 W_r 异或到 W_{r+4}，临时得到 W'_r；使用后再次调用同一函数，
即可还原 W_{r+4}。
"""

from projectq.ops import CNOT

from linear_permutations import P1_CNOTS, apply_classical, p1, rotl32


def xor_word(source, target, rotate_left=0):
    if len(source) != 32 or len(target) != 32:
        raise ValueError("SM3 words must have 32 qubits")
    for bit in range(32):
        CNOT | (source[(bit - rotate_left) % 32], target[bit])


def expand_next(window, word_index):
    """在 16 字循环寄存器中用 W_j 覆盖 W_{j-16}。"""
    if len(window) != 16 or not 16 <= word_index <= 67:
        raise ValueError("Need 16 words and j in 16..67")
    j = word_index
    target = window[j % 16]
    xor_word(window[(j - 9) % 16], target)
    xor_word(window[(j - 3) % 16], target, rotate_left=15)
    p1(target)
    xor_word(window[(j - 13) % 16], target, rotate_left=7)
    xor_word(window[(j - 6) % 16], target)
    return target


def derive_w_prime_in_place(window, round_index):
    """执行 W_{r+4} ^= W_r；使用 W'_r 后再次调用可还原原值。"""
    if len(window) != 16 or not 0 <= round_index <= 63:
        raise ValueError("Need 16 words and round in 0..63")
    source = window[round_index % 16]
    target = window[(round_index + 4) % 16]
    xor_word(source, target)
    return target


def classical_expand(input_words):
    """循环线路的经典对应实现，仅用于核对位序。"""
    if len(input_words) != 16:
        raise ValueError("Need 16 input words")
    window = input_words[:]
    words = input_words[:]
    for j in range(16, 68):
        value = window[j % 16]
        value ^= window[(j - 9) % 16]
        value ^= rotl32(window[(j - 3) % 16], 15)
        value = apply_classical(value, P1_CNOTS)
        value ^= rotl32(window[(j - 13) % 16], 7)
        value ^= window[(j - 6) % 16]
        value &= 0xFFFFFFFF
        window[j % 16] = value
        words.append(value)
    return words


def self_check():
    from random import Random

    rng = Random(20261123)
    for _ in range(100):
        inputs = [rng.getrandbits(32) for _ in range(16)]
        expected = inputs[:]
        for j in range(16, 68):
            value = expected[j - 16] ^ expected[j - 9] ^ rotl32(expected[j - 3], 15)
            value = value ^ rotl32(value, 15) ^ rotl32(value, 23)
            value ^= rotl32(expected[j - 13], 7) ^ expected[j - 6]
            expected.append(value & 0xFFFFFFFF)
        actual = classical_expand(inputs)
        assert actual == expected


if __name__ == "__main__":
    self_check()
    print("Rolling W16..W67 expansion matched 100 independent vectors.")
    print("Each update: 196 CNOTs; all 52 updates: 10192 CNOTs.")
