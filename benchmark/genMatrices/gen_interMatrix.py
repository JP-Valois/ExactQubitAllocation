from qiskit import QuantumCircuit
from utils import * # type: ignore
import numpy as np


## Benchmark circuit names
'''
Benchmarks = ["10_qft",
              "10_sqn",
              "10_sym9",
              "11_sym9",
              "11_wim",
              "11_z4",
              "12_cycle10",
              "12_rd84",
              "12_sym9",
              "13_dist",
              "13_radd",
              "13_root",
              "14_clip",
              "14_cm42a",
              "14_cm85a",
              "15_co14",
              "15_misex1",
              "15_sqrt7",
              "16_inc",
              "16_ising",
              "16_mlp4",
              "rd_11",
              "rd_12",
              "rd_13",
              "rd_14",
              "rd_15",
              "rd_16",
              "rd_17",
              "rd_18",
              "rd_19",
              "rd_20",
              "rd_21",
              "rd_22",
              "rd_23",
              "rd_24",
              "rd_25",
              "rd_26",
              "rd_27"]
'''

# TEST
'''
Benchmarks = ["10_dj",
              "10_grover",
              "10_qaoa",
              "10_qfte",
              "10_wstate",
              "16_alu2",
              "17_5xp1",
              "17_ryy6",
              "17_t481",
              "18_alu3",
              "18_sqrt6",
              "19_add6",
              "20_cmb",
              "20_ex1010",
              "21_decod",
              "21_dk17",
              "22_apla",
              "22_cm105a",
              "23_63mod4096",
              "24_arb8",
              "24_qpe",
              "25_cu",
              "25_rd73",
              "26_in0",
              "27_grover",
              "27_sym9"]
'''

Benchmarks = ["10_vqe",
              "24_pe",
              "26_ising",
              "27_grover",
              "27_wstate"]


def InteractionMatrix (qasm_file_path):
    # Load the circuit
    qc = QuantumCircuit.from_qasm_file(qasm_file_path)
    num_qubits = qc.num_qubits

    # Initialize the interaction matrix F
    F = np.zeros((num_qubits, num_qubits), dtype=int)

    # Count two-qubit gates and populate the interaction matrix
    two_qubit_gate_count = 0

    for instruction in qc.data:

        qubits = instruction.qubits

        if len(qubits) == 2:
            qubit_indices = [qubit._index for qubit in qubits]
            i, j = qubit_indices
            if (i != j):
                two_qubit_gate_count += 1  
                F[i, j] += 1
                F[j, i] += 1

    print(f"Loaded circuit '{qasm_file_path}' with {num_qubits} qubits and {two_qubit_gate_count} two-qubit gates.")

    return F


k = 0

for bench in Benchmarks:

    input_file_path = "../circuits/" + bench + ".qasm"

    output_file_path = "../inter/" + bench + ".csv"
    #output_file_path = "../inter/MQT/" + bench + ".csv"

    SaveMatrix_csv(InteractionMatrix(input_file_path), output_file_path)

    k += 1

print(f"Succesfully loaded {k} circuits")
