"""使用 ProjectQ 构造一个消息分组的 64 轮 SM3 门级线路。

本模块用于核对论文第 2.3、3.3 和 3.4 节的实现，采用原地与异地加法器、
各含 68 个 CNOT 门的 P0/P1 网络，以及论文中的单比特 QFA/FF/GG 线路。
两条状态线路交替运行：先构造下一轮 COMP，再清理上一轮 COMP。

"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass

from projectq import MainEngine
from projectq.backends import ResourceCounter
from projectq.ops import CNOT, X

from fixed_branch_accounting import QandDaggerCalls, add_budget
from linear_permutations import p0
from message_schedule import derive_w_prime_in_place, expand_next, xor_word
from reversible_out_of_place import AdderState, cleanup as adder_cleanup, compute as adder_compute
from sm3_gate_blocks import (
    ff_cleanup,
    ff_forward,
    gg_cleanup,
    gg_forward,
    qfa_words_cleanup,
    qfa_words_forward,
)
from optimal_depth_adder import optimal_depth_add_in_place
from parallel_in_place_adders import add_four_parallel


WORD_BITS = 32
IV = (0x7380166F, 0x4914B2B9, 0x172442D7, 0xDA8A0600,
      0xA96F30BC, 0x163138AA, 0xE38DEE4D, 0xB0FB0E4E)


def rotl_view(word, distance):
    """与论文一致，循环移位只改变逻辑位编号。"""
    distance %= WORD_BITS
    return [word[(bit - distance) % WORD_BITS] for bit in range(WORD_BITS)]


def set_constant(word, value):
    for bit in range(WORD_BITS):
        if (value >> bit) & 1:
            X | word[bit]


@dataclass
class StateWords:
    A: object
    B: object
    C: object
    D: object
    E: object
    F: object
    G: object
    H: object

    def as_tuple(self):
        return (self.A, self.B, self.C, self.D, self.E, self.F, self.G, self.H)

    def next_round(self):
        """图 9 的 SWAP 只重排寄存器引用，不产生实际量子门。"""
        return StateWords(
            self.D, self.A, rotl_view(self.B, 9), self.C,
            self.H, self.E, rotl_view(self.F, 19), self.G,
        )


@dataclass
class CompressionContext:
    round_index: int
    state: StateWords
    window: object
    t_constant: object
    t_rotated: object
    a_rotated: object
    w: object
    w_prime: object
    k1: object
    k2: object
    k3: object
    ff: object
    gg: object
    ss1_adder: AdderState
    upper_adder: AdderState
    lower_adder: AdderState
    ss1: object


def _apply_boolean_words(eng, round_index, state, ff, gg):
    for bit in range(WORD_BITS):
        ff_forward(eng, round_index, state.A[bit], state.B[bit], state.C[bit], ff[bit])
        gg_forward(eng, round_index, state.E[bit], state.F[bit], state.G[bit], gg[bit])


def _cleanup_boolean_words(eng, round_index, state, ff, gg):
    for bit in range(WORD_BITS - 1, -1, -1):
        gg_cleanup(eng, round_index, state.E[bit], state.F[bit], state.G[bit], gg[bit])
        ff_cleanup(eng, round_index, state.A[bit], state.B[bit], state.C[bit], ff[bit])


def comp_forward(eng, state, window, round_index, t_constant=None):
    """执行图 9 中 COMP 的前四层，并保留两路 TT 加法所需的输出。"""
    r = round_index
    if r >= 12:
        expand_next(window, r + 4)
    w = window[r % 16]
    w_prime = derive_w_prime_in_place(window, r)
    if t_constant is None:
        t_constant = eng.allocate_qureg(WORD_BITS)
    set_constant(t_constant, 0x79CC4519 if r < 16 else 0x7A879D8A)
    t_rotated = rotl_view(t_constant, r % 32)
    a_rotated = rotl_view(state.A, 12)

    k1 = eng.allocate_qureg(WORD_BITS)
    k2 = eng.allocate_qureg(WORD_BITS)
    k3 = eng.allocate_qureg(WORD_BITS)
    ff = eng.allocate_qureg(WORD_BITS)
    gg = eng.allocate_qureg(WORD_BITS)

    qfa_words_forward(eng, a_rotated, state.E, t_rotated, k1)
    ss1_adder = adder_compute(eng, k1, t_rotated)
    ss1 = rotl_view(ss1_adder.output, 7)

    # 在 A 和 E 仍保存原始状态值时计算 FF 与 GG。
    _apply_boolean_words(eng, r, state, ff, gg)
    xor_word(ss1, a_rotated)  # 按循环移位后的逻辑视图，A 此时保存 SS2。

    qfa_words_forward(eng, ss1, gg, w, k2)
    qfa_words_forward(eng, a_rotated, ff, w_prime, k3)
    upper_adder = adder_compute(eng, k2, w)
    lower_adder = adder_compute(eng, k3, w_prime)
    return CompressionContext(
        r, state, window, t_constant, t_rotated, a_rotated, w, w_prime,
        k1, k2, k3, ff, gg,
        ss1_adder, upper_adder, lower_adder, ss1,
    )


def add_tt_to_both_lines(eng, context, first_line, second_line):
    """按图 11 复制两路部分和，并分别加到两条状态线路的 D/H。"""
    upper = context.upper_adder.output
    lower = context.lower_adder.output
    upper_copy = eng.allocate_qureg(WORD_BITS)
    lower_copy = eng.allocate_qureg(WORD_BITS)
    xor_word(upper, upper_copy)
    xor_word(lower, lower_copy)

    add_four_parallel(eng, [
        (lower, first_line.D),
        (upper, first_line.H),
        (lower_copy, second_line.D),
        (upper_copy, second_line.H),
    ])

    xor_word(lower, lower_copy)
    xor_word(upper, upper_copy)
    p0(first_line.H)
    p0(second_line.H)
    # 下一轮开始前，显式释放已清零的两个副本。
    del upper_copy, lower_copy


def comp_cleanup(eng, context):
    """严格按数据依赖的逆序执行图 9 中的 COMP†。"""
    r = context.round_index
    s = context.state
    adder_cleanup(eng, context.lower_adder)
    adder_cleanup(eng, context.upper_adder)
    qfa_words_cleanup(eng, context.a_rotated, context.ff,
                      context.w_prime, context.k3)
    qfa_words_cleanup(eng, context.ss1, context.gg,
                      context.w, context.k2)
    xor_word(context.ss1, context.a_rotated)
    _cleanup_boolean_words(eng, r, s, context.ff, context.gg)
    adder_cleanup(eng, context.ss1_adder)
    qfa_words_cleanup(eng, context.a_rotated, s.E,
                      context.t_rotated, context.k1)
    set_constant(context.t_constant, 0x79CC4519 if r < 16 else 0x7A879D8A)
    derive_w_prime_in_place(context.window, r)


def allocate_block_inputs(eng):
    """分配符号化的 512 位消息分组和两条各 256 位的状态线路。"""
    window = [eng.allocate_qureg(WORD_BITS) for _ in range(16)]
    line_1 = StateWords(*(eng.allocate_qureg(WORD_BITS) for _ in range(8)))
    line_2 = StateWords(*(eng.allocate_qureg(WORD_BITS) for _ in range(8)))
    for source, target in zip(line_1.as_tuple(), line_2.as_tuple()):
        xor_word(source, target)
    return window, [line_1, line_2]


def compress_block(eng, window, lines, rounds=64):
    """运行至多 64 轮，并为每条状态线路保留一个常驻 T 寄存器。"""
    if not 1 <= rounds <= 64:
        raise ValueError("rounds must be 1..64")
    t_registers = [eng.allocate_qureg(WORD_BITS) for _ in range(2)]
    pending = None
    for r in range(rounds):
        active = r % 2
        passive = 1 - active
        context = comp_forward(eng, lines[active], window, r,
                               t_constant=t_registers[active])
        if pending is not None:
            old_line = pending.state
            comp_cleanup(eng, pending)
            lines[passive] = old_line.next_round()
            pending = None
        add_tt_to_both_lines(eng, context, lines[active], lines[passive])
        lines[passive] = lines[passive].next_round()
        pending = context
    if pending is not None:
        old_line = pending.state
        comp_cleanup(eng, pending)
        lines[(rounds - 1) % 2] = old_line.next_round()
    return lines


def summarize(resource):
    counts = {"T": 0, "CNOT": 0, "one_qubit_clifford": 0,
              "measurement": 0, "other": {}}
    for (gate, controls), number in resource.gate_counts.items():
        name = str(gate)
        if name in ("Allocate", "Deallocate", "FlushGate"):
            continue
        if name == "T" or name.startswith("T^"):
            counts["T"] += number
        elif name == "X" and controls == 1:
            counts["CNOT"] += number
        elif name in ("H", "S", "S^\\dagger", "X") and controls == 0:
            counts["one_qubit_clifford"] += number
        elif name == "Measure":
            counts["measurement"] += number
        else:
            counts["other"][f"{name} controls={controls}"] = number
    counts["max_width"] = resource.max_width
    return counts


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rounds", type=int, default=1)
    args = parser.parse_args()
    resource = ResourceCounter()
    # 直接统计发出的门：默认局部优化器会删除空闲的预留量子比特，
    # 从而使四路加法器的工作空间预算在统计中消失。
    eng = MainEngine(resource, engine_list=[])
    with QandDaggerCalls() as qand_calls:
        window, lines = allocate_block_inputs(eng)
        compress_block(eng, window, lines, args.rounds)
        eng.flush()
    raw = summarize(resource)
    correction = qand_calls.correction_budget(include_standard=False)
    report = add_budget(raw, correction)
    report.update({
        "rounds": args.rounds,
        "accounting": "paper_fixed_branch",
        "qand_dagger_calls": {
            "standard_adder": qand_calls.standard,
            "figure_3b_qfa_ff_gg": qand_calls.figure,
        },
        "fixed_correction_budget_added": correction,
    })
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
