"""使用真值表验证各单比特门级部件。"""

from itertools import product

from projectq import MainEngine
from projectq.backends import Simulator
from projectq.ops import All, H, Measure, X

from sm3_gate_blocks import (
    ff_cleanup,
    ff_forward,
    gg_cleanup,
    gg_forward,
    qfa_cleanup,
    qfa_forward,
    quantum_and_dag_figure,
)
from optimal_depth_adder import quantum_and


def run_case(kind, values, round_index=0):
    eng = MainEngine(Simulator())
    qubits = eng.allocate_qureg(4)
    for index, value in enumerate(values):
        if value:
            X | qubits[index]
    a, b, c, out = qubits
    if kind == "qfa":
        qfa_forward(eng, a, b, c, out)
        expected_output = (values[0] & values[1]) ^ (values[0] & values[2]) ^ (values[1] & values[2])
        expected_c = values[0] ^ values[1] ^ values[2]
        Measure | out
        Measure | c
        eng.flush()
        assert int(out) == expected_output and int(c) == expected_c
        qfa_cleanup(eng, a, b, c, out)
    elif kind == "ff":
        ff_forward(eng, round_index, a, b, c, out)
        expected_output = (values[0] ^ values[1] ^ values[2]) if round_index < 16 else int(sum(values) >= 2)
        expected_c = values[2]
        Measure | out
        eng.flush()
        assert int(out) == expected_output
        ff_cleanup(eng, round_index, a, b, c, out)
    else:
        gg_forward(eng, round_index, a, b, c, out)
        expected_output = (values[0] ^ values[1] ^ values[2]) if round_index < 16 else (values[1] if values[0] else values[2])
        expected_c = values[2]
        Measure | out
        eng.flush()
        assert int(out) == expected_output
        gg_cleanup(eng, round_index, a, b, c, out)
    All(Measure) | qubits
    eng.flush()
    actual = tuple(int(qubit) for qubit in qubits)
    assert actual == (*values, 0), (kind, round_index, values, actual)
    return expected_output, expected_c


def main():
    for left, right in product((0, 1), repeat=2):
        eng = MainEngine(Simulator())
        qubits = eng.allocate_qureg(3)
        if left:
            X | qubits[0]
        if right:
            X | qubits[1]
        quantum_and(eng, qubits[0], qubits[1], qubits[2])
        quantum_and_dag_figure(eng, qubits[0], qubits[1], qubits[2])
        All(Measure) | qubits
        eng.flush()
        assert tuple(int(qubit) for qubit in qubits) == (left, right, 0)
    for values in product((0, 1), repeat=3):
        for kind, round_index in (("qfa", 0), ("ff", 0), ("ff", 16), ("gg", 0), ("gg", 16)):
            run_case(kind, values, round_index)
    # 清理线路还必须恢复叠加态各分量之间的相对相位。
    for kind, forward, reverse in (
        ("qfa", qfa_forward, qfa_cleanup),
        ("ff", lambda e, a, b, c, o: ff_forward(e, 16, a, b, c, o),
         lambda e, a, b, c, o: ff_cleanup(e, 16, a, b, c, o)),
        ("gg", lambda e, a, b, c, o: gg_forward(e, 16, a, b, c, o),
         lambda e, a, b, c, o: gg_cleanup(e, 16, a, b, c, o)),
    ):
        eng = MainEngine(Simulator())
        qubits = eng.allocate_qureg(4)
        for qubit in qubits[:3]:
            H | qubit
        forward(eng, *qubits)
        reverse(eng, *qubits)
        eng.flush()
        mapping, amplitudes = eng.backend.cheat()
        reference = amplitudes[0]
        valid_indices = set()
        for bits in product((0, 1), repeat=3):
            index = sum(bit << mapping[qubits[j].id] for j, bit in enumerate(bits))
            valid_indices.add(index)
            amplitude = amplitudes[index]
            assert abs(amplitude - reference) < 1e-8, (kind, bits, amplitude, reference)
        assert sum(abs(value) ** 2 for i, value in enumerate(amplitudes)
                   if i not in valid_indices) < 1e-8
        All(Measure) | qubits
        eng.flush()
    print("44 forward/cleanup cases returned all inputs and outputs correctly.")


if __name__ == "__main__":
    main()
