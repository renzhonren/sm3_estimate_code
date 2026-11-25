"""在较小位宽下验证异地加法器的正算与反算清理。"""

from projectq import MainEngine
from projectq.backends import Simulator
from projectq.ops import All, H, Measure, X

from reversible_out_of_place import cleanup, compute


def read_word(word):
    return sum(int(qubit) << bit for bit, qubit in enumerate(word))


def run_case(n, left, right):
    eng = MainEngine(Simulator())
    a = eng.allocate_qureg(n)
    b = eng.allocate_qureg(n)
    for bit in range(n):
        if (left >> bit) & 1:
            X | a[bit]
        if (right >> bit) & 1:
            X | b[bit]
    state = compute(eng, a, b)
    All(Measure) | state.output
    eng.flush()
    assert read_word(state.output) == (left + right) % (1 << n)
    cleanup(eng, state)
    All(Measure) | a
    All(Measure) | b
    All(Measure) | state.output
    All(Measure) | state.ancilla_location
    eng.flush()
    assert read_word(a) == left
    assert read_word(b) == right
    assert read_word(state.output) == 0
    assert all(int(qubit) == 0 for qubit in state.ancilla_location)


def main():
    for left, right in ((0, 0), (1, 2), (3, 7), (7, 1)):
        run_case(3, left, right)
    for left, right in ((0, 0), (1, 10), (15, 15)):
        run_case(4, left, right)
    eng = MainEngine(Simulator())
    a = eng.allocate_qureg(3)
    b = eng.allocate_qureg(3)
    H | a[0]
    H | b[1]
    state = compute(eng, a, b)
    cleanup(eng, state)
    eng.flush()
    mapping, amplitudes = eng.backend.cheat()
    valid = {
        (left << mapping[a[0].id]) | (right << mapping[b[1].id])
        for left in (0, 1) for right in (0, 1)
    }
    reference = amplitudes[0]
    assert all(abs(amplitudes[index] - reference) < 1e-8 for index in valid)
    assert sum(abs(amplitude) ** 2 for index, amplitude in enumerate(amplitudes)
               if index not in valid) < 1e-8
    All(Measure) | a
    All(Measure) | b
    All(Measure) | state.output
    All(Measure) | state.ancilla_location
    eng.flush()
    print("Seven forward/cleanup adder cases passed at widths 3 and 4.")


if __name__ == "__main__":
    main()
