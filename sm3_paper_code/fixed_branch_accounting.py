"""为基于测量的 QAND 反算提供与测量结果无关的资源预算。

ResourceCounter 只执行确定性的测量零分支。包装函数记录每次 QAND† 调用，
但不向实际线路添加门。随后将条件修正的固定预算加入资源报告，
使统计结果不随测量结果为 0 或 1 而变化。
"""

from dataclasses import dataclass

import optimal_depth_adder as adder
import sm3_gate_blocks as blocks


@dataclass
class QandDaggerCalls:
    standard: int = 0
    figure: int = 0

    def __enter__(self):
        self._original_standard = adder.quantum_and_dag
        self._original_figure = blocks.quantum_and_dag_figure

        def standard(eng, a, b, target):
            self.standard += 1
            return self._original_standard(eng, a, b, target)

        def figure(eng, a, b, target):
            self.figure += 1
            return self._original_figure(eng, a, b, target)

        adder.quantum_and_dag = standard
        blocks.quantum_and_dag_figure = figure
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        adder.quantum_and_dag = self._original_standard
        blocks.quantum_and_dag_figure = self._original_figure

    def correction_budget(self, *, include_standard=True, include_figure=True):
        """计算每次反算的条件修正分支所需门数。

        原地加法器中的标准 QAND†：2 个 CNOT 门和 4 个单比特
        Clifford 门（S、S、X、S†）。
        图 3(b) 中 QFA/FF/GG 使用的 QAND†：1 个 CNOT 门和 1 个 X 门。
        两类条件修正分支均不使用 T 门。
        """
        standard = self.standard if include_standard else 0
        figure = self.figure if include_figure else 0
        return {
            "T": 0,
            "CNOT": 2 * standard + figure,
            "one_qubit_clifford": 4 * standard + figure,
        }


def add_budget(raw, budget):
    result = dict(raw)
    for key in ("T", "CNOT", "one_qubit_clifford"):
        result[key] = raw[key] + budget[key]
    return result
