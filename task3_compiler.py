#!/usr/bin/env python

"""
Solves Task 3. Simple compiler – program, which translates one quantum circuit
into another, using a restricted set of gates.

Basic gates for the input circuit, such as (I, H, X, Y, Z, RX, RY, RZ, CNOT,
CZ), are considered.
The output circuit should consists only from gates of the restricted set
(RX, RZ, CZ) only.

Analyzing the overhead:
The overhead depends on the basis set to which one we compare it,
as every gate of the basis set needs to be rebuild by the new set.
If we compare to (I, H, X, Y, Z, RX, RY, RZ, CNOT, CZ), we have:

CX = 7 gates
I  = 4 gates
X  = 4 gates
Y  = 4 gates
H  = 3 gates
RY = 3 gates
Z  = 2 gates

So if all these gates are used, we have an overhead of
(were the name stands for the number of the corresponding gates):
(CX ** 7) * (I ** 4) * (X ** 4) * (Y ** 4) * (H ** 4) * (RY ** 3) * (Z ** 4)

When optimizing the circuited I can be removed (set to I=0 in above formula)
"""

import numpy as np
import qiskit as q
from qiskit import *
from qiskit.tools.visualization import plot_histogram
from qiskit import QuantumCircuit
from cvxopt import matrix, printing
# %matplotlib inline

__author__ = "Daniel Scheiermann"
__copyright__ = ""
__credits__ = [""]
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Daniel Scheiermann"
__email__ = "daniel.scheiermann@stud.uni-hannover.de"
__status__ = "Production"

CXGATES = 7
HGATES = 3
XGATES = 4
YGATES = 4

def do_job_simulate(circuit, shots=8000, unitary_sim=False):
    if unitary_sim:
        unitary_result = q.execute(circuit,
            backend=q.Aer.get_backend('unitary_simulator')).result()
        unitary = unitary_result.get_unitary(circuit, decimals=3)

        return unitary

    else:
        statevec_result = q.execute(circuit,
            backend=q.Aer.get_backend("statevector_simulator")).result()
        state_vec = statevec_result.get_statevector(circuit)
        n_qubits = circuit.n_qubits
        circuit.measure([i for i in range(n_qubits)],
            [i for i in range(len(circuit.clbits))])
        qasm_job = q.execute(circuit,
            backend=q.Aer.get_backend("qasm_simulator"), shots=shots).result()
        counts = qasm_job.get_counts()

        return state_vec, counts


def print_matrix(unitary_simulated, decimals=2):
    printing.options['width'] = -1
    printing.options['dformat'] = "%." + str(decimals) + "f"
    unitary_simulated_print = matrix(unitary_simulated,
        unitary_simulated.shape)
    print(unitary_simulated_print)


def my_test(circ_test, circ_ideal):
    print("Tested circuit")
    tested_matrix = do_job_simulate(circ_test, unitary_sim=True)
    print_matrix(tested_matrix)
    print("Ideal circuit")
    ideal_matrix = do_job_simulate(circ_ideal, unitary_sim=True)
    print_matrix(ideal_matrix)
    if np.allclose(ideal_matrix, tested_matrix):
        print("Gate correctly implemented!\n")
        
        return True
    else:
        print("Gate WRONG!\n")
        return False
        


def my_I(circuit, qubit):
    # No gate is also id
    circuit.rz(np.pi/2, qubit)
    circuit.rz(np.pi/2, qubit)
    circuit.rz(np.pi/2, qubit)
    circuit.rz(np.pi/2, qubit)

    return circuit


def my_h(circuit, qubit):
    circuit.rz(np.pi/2, qubit)
    circuit.rx(np.pi/2, qubit)
    circuit.rz(np.pi/2, qubit)

    return circuit


def my_x(circuit, qubit):
    # circuit.h(qubit)
    # circuit.z(qubit)
    # circuit.h(qubit)
    circuit.rz(np.pi/2, qubit)
    circuit.rx(np.pi/2, qubit)
    circuit.rx(np.pi/2, qubit)
    circuit.rz(np.pi/2, qubit)

    return circuit


def my_cx(circuit, control_qubit=0, target_qubit=1):
    my_h(circuit, target_qubit)
    circuit.cz(control_qubit, target_qubit)
    my_h(circuit, target_qubit)

    return circuit


def my_y(circuit, qubit):
    circuit.rx(np.pi/2, qubit)
    circuit.rx(np.pi/2, qubit)
    circuit.rz(np.pi/2, qubit)
    circuit.rz(np.pi/2, qubit)

    return circuit


def my_z(circuit, qubit):
    circuit.rz(np.pi/2, qubit)
    circuit.rz(np.pi/2, qubit)

    return circuit

def my_ry(circuit, qubit):
    circuit.rz(np.pi/2, qubit)
    circuit.rz(np.pi/2, qubit)
    my_h(circuit, qubit)

    return circuit


def ideal_i(circuit, qubit):
    circuit.i(qubit)

    return circuit


def ideal_rx(circuit, qubit):
    circuit.rx(np.pi/2, qubit)

    return circuit


def ideal_ry(circuit, qubit):
    circuit.ry(np.pi/2, qubit)

    return circuit


def ideal_rz(circuit, qubit):
    circuit.rz(np.pi/2, qubit)

    return circuit

def ideal_x(circuit, qubit):
    circuit.x(qubit)

    return circuit

def ideal_y(circuit, qubit):
    circuit.y(qubit)

    return circuit

def ideal_z(circuit, qubit):
    circuit.z(qubit)

    return circuit


def ideal_cz(circuit, target_qubit, control_qubit=1):
    """
    target_qubit: is manipulated if control_qubit has spin up
    """
    circuit.cz(control_qubit, target_qubit)

    return circuit


def ideal_cx(circuit, target_qubit, control_qubit=1):
    """
    target_qubit: is manipulated if control_qubit has spin up
    """
    circuit.cx(control_qubit, target_qubit)

    return circuit


def brute_force_one_qubit_gate2(basis_set, target_set):
    """
    target_set is QuantumCircuit, which should be build
    by 3 components of the basis_set
    All possible combination for 2 gates of the basis_set are tested
    """
    # TODO: merge with brute_force_one_qubit_gate3()
    print("Used basis set:", basis_set)
    target_qubit = 0
    for i, first in enumerate(basis_set):
        for j, second in enumerate(basis_set):
                test_circ = QuantumCircuit(1)
                first(test_circ, target_qubit)
                second(test_circ, target_qubit)
                tested_matrix = do_job_simulate(test_circ, unitary_sim=True)
                for ideal in target_set:
                    ideal_circ = QuantumCircuit(1)
                    ideal_matrix = do_job_simulate(ideal(ideal_circ,
                        target_qubit), unitary_sim=True) 
                    if np.allclose(ideal_matrix, tested_matrix):
                        print(first, second)
                        print("Gate correctly implemented!\n")


def brute_force_one_qubit_gate3(basis_set, target_set):
    """
    target_set is QuantumCircuit, which should be build
    by 3 components of the basis_set
    All possible combination for 3 gates of the basis_set are tested
    """
    print("Used basis set:", basis_set)
    target_qubit = 0
    for i, first in enumerate(basis_set):
        for j, second in enumerate(basis_set):
            for k, third in enumerate(basis_set):
                test_circ = QuantumCircuit(1)
                first(test_circ, target_qubit)
                second(test_circ, target_qubit)
                third(test_circ, target_qubit)
                tested_matrix = do_job_simulate(test_circ, unitary_sim=True)
                for ideal in target_set:
                    ideal_circ = QuantumCircuit(1)
                    ideal_matrix = do_job_simulate(ideal(ideal_circ,
                        target_qubit), unitary_sim=True) 
                    if np.allclose(ideal_matrix, tested_matrix):
                        print(first, second, third)
                        print("Gate correctly implemented!\n")


def brute_force_one_qubit_gate4(basis_set, target_set):
    """
    target_set is QuantumCircuit, which should be build
    by 4 components of the basis_set.
    All possible combination for 4 gates of the basis_set are tested
    """
    # TODO: merge with brute_force_one_qubit_gate3()
    print("Used basis set:", basis_set)
    target_qubit = 0
    for i, first in enumerate(basis_set):
        for j, second in enumerate(basis_set):
            for k, third in enumerate(basis_set):
                for m, fourth in enumerate(basis_set):
                    test_circ = QuantumCircuit(1)
                    first(test_circ, target_qubit)
                    second(test_circ, target_qubit)
                    third(test_circ, target_qubit)
                    fourth(test_circ, target_qubit)
                    tested_matrix = do_job_simulate(test_circ,
                        unitary_sim=True)
                    for ideal in target_set:
                        ideal_circ = QuantumCircuit(1)
                        ideal_matrix = do_job_simulate(ideal(ideal_circ,
                            target_qubit), unitary_sim=True) 
                        if np.allclose(ideal_matrix, tested_matrix):
                            print(first, second, third, fourth)
                            print("Gate correctly implemented!\n")


def brute_force_two_qubit_gate6(basis_set, target_set):
    """
    target_set is QuantumCircuit, which should be build
    by 6 components of the basis_set.
    All possible combination for 6 gates of the basis_set are tested
    This way a flipped CX was found
    """
    print("Used basis set:", basis_set)
    target_qubit = 0
    for i, first in enumerate(basis_set):
        for j, second in enumerate(basis_set):
            for k, third in enumerate(basis_set):
                for m, fourth in enumerate(basis_set):
                    for n, fifth in enumerate(basis_set):
                        for o, sixth in enumerate(basis_set):
                            test_circ = QuantumCircuit(2)
                            first(test_circ, target_qubit)
                            second(test_circ, target_qubit)
                            third(test_circ, target_qubit)
                            fourth(test_circ, target_qubit)
                            fifth(test_circ, target_qubit)
                            sixth(test_circ, target_qubit)
                            tested_matrix = do_job_simulate(test_circ,
                                unitary_sim=True)
                            for ideal in target_set:
                                ideal_circ = QuantumCircuit(2)
                                ideal_matrix = do_job_simulate(
                                    ideal(ideal_circ, target_qubit),
                                    unitary_sim=True) 
                                if np.allclose(ideal_matrix, tested_matrix):
                                    print(first, second, third, fourth,
                                          fifth, sixth)
                                    print("Gate correctly implemented!\n")


def brute_force_two_qubit_gate3(basis_set, target_set):
    """
    target_set is QuantumCircuit, which should be build
    by 3 components of the basis_set.
    All possible combination for 3 gates of the basis_set are tested
    This way a flipped CX was found
    """
    print("Used basis set:", basis_set)
    target_qubit = 0
    for i, first in enumerate(basis_set):
        for j, second in enumerate(basis_set):
            for k, third in enumerate(basis_set):
                test_circ = QuantumCircuit(2)
                first(test_circ, target_qubit)
                second(test_circ, target_qubit)
                third(test_circ, target_qubit)
                tested_matrix = do_job_simulate(test_circ, unitary_sim=True)
                for ideal in target_set:
                    ideal_circ = QuantumCircuit(2)
                    ideal_matrix = do_job_simulate(ideal(ideal_circ,
                        target_qubit), unitary_sim=True) 
                    if np.allclose(ideal_matrix, tested_matrix):
                        print(first, second, third)
                        print("Gate correctly implemented!\n")


def testcases():
    number_testcases = 0
    passed_list = []

    i_ideal = QuantumCircuit(1)
    i = QuantumCircuit(1)
    my_I(i, 0)
    i_ideal.id(0)
    print("I-Gate")
    number_testcases += 1
    passed_list.append(my_test(i, i_ideal))

    h_ideal = QuantumCircuit(1)
    h = QuantumCircuit(1)
    my_h(h, 0)
    h_ideal.h(0)
    print("H-Gate")
    number_testcases += 1
    passed_list.append(my_test(h, h_ideal))

    x_ideal = QuantumCircuit(1)
    x = QuantumCircuit(1)
    my_x(x, 0)
    x_ideal.x(0)
    print("X-Gate")
    number_testcases += 1
    passed_list.append(my_test(x, x_ideal))

    y_ideal = QuantumCircuit(1)
    y = QuantumCircuit(1)
    my_y(y, 0)
    y_ideal.y(0)
    print("Y-Gate")
    number_testcases += 1
    passed_list.append(my_test(y, y_ideal))

    z_ideal = QuantumCircuit(1)
    z = QuantumCircuit(1)
    my_z(z, 0)
    z_ideal.z(0)
    print("Z-Gate")
    number_testcases += 1
    passed_list.append(my_test(z, z_ideal))

    ry_ideal = QuantumCircuit(1)
    ry = QuantumCircuit(1)
    ry_ideal.ry(np.pi/2, 0)
    my_ry(ry, 0)
    print("RY-Gate")
    number_testcases += 1
    passed_list.append(my_test(ry, ry_ideal))

    cx_ideal = QuantumCircuit(2)
    cx = QuantumCircuit(2)
    my_cx(cx, 0, 1)
    cx_ideal.cx(0, 1)
    print("CX-Gate")
    number_testcases += 1
    passed_list.append(my_test(cx, cx_ideal))

    circ = QuantumCircuit(2)
    circ.rz(np.pi/2, 0)
    circ.x(0)
    circ.cx(0, 1)
    print("Original circuit:")
    print(circ)
    print("Transpiled circuit:")
    transpiled_circ = my_transpile(circ)
    print(transpiled_circ)
    number_testcases += 1
    passed_list.append(my_test(circ, transpiled_circ))
    print("Optimzed circuit:")
    optimized_circ = my_transpile_optimized(circ)
    print(optimized_circ)
    passed_list.append(my_test(circ, optimized_circ))
    number_testcases += 1

    # IMPORTANT: As qiskit has its own transpiler it can intervene here
    # Especially this happens for CX Y, so barriers are needed
    circ = QuantumCircuit(2)
    circ.h(0)
    circ.barrier()
    circ.h(1)
    circ.barrier()
    circ.cx(0, 1)
    circ.barrier()
    circ.z(1)
    circ.id(1)
    circ.barrier()
    circ.cx(0, 1)
    circ.barrier()
    circ.y(1)
    circ.barrier()
    circ.ry(np.pi/2, 0)
    circ.barrier()
    circ.x(1)
    print("Original circuit:")
    print(circ)
    print("Transpiled circuit:")
    transpiled_circ = my_transpile(circ)
    print(transpiled_circ)
    passed_list.append(my_test(circ, transpiled_circ))
    number_testcases += 1
    print("Optimzed circuit:")
    optimized_circ = my_transpile_optimized(circ)
    print(optimized_circ)
    passed_list.append(my_test(circ, optimized_circ))
    number_testcases += 1

    # IMPORTANT: As qiskit has its own transpiler it can intervene here
    # Especially this happens for CX Y, so barriers are needed
    circ = QuantumCircuit(2)
    circ.cz(0, 1)
    circ.barrier()
    circ.cz(0, 1)
    circ.barrier()
    print("Original circuit:")
    print(circ)
    print("Transpiled circuit:")
    transpiled_circ = my_transpile(circ)
    print(transpiled_circ)
    passed_list.append(my_test(circ, transpiled_circ))
    number_testcases += 1
    print("Optimzed circuit:")
    optimized_circ = my_transpile_optimized(circ)
    print(optimized_circ)
    passed_list.append(my_test(circ, optimized_circ))
    number_testcases += 1

    # IMPORTANT: As qiskit has its own transpiler it can intervene here
    # Especially this happens for CX Y, so barriers are needed
    circ = QuantumCircuit(2)
    circ.cz(0, 1)
    circ.barrier()
    circ.cz(1, 0)
    circ.barrier()
    print("Original circuit:")
    print(circ)
    print("Transpiled circuit:")
    transpiled_circ = my_transpile(circ)
    print(transpiled_circ)
    passed_list.append(my_test(circ, transpiled_circ))
    number_testcases += 1
    print("Optimzed circuit:")
    optimized_circ = my_transpile_optimized(circ)
    print(optimized_circ)
    passed_list.append(my_test(circ, optimized_circ))
    number_testcases += 1

    # IMPORTANT: As qiskit has its own transpiler it can intervene here
    # Especially this happens for CX Y, so barriers are needed
    circ = QuantumCircuit(2)
    # circ.cx(0, 1)
    circ.x(1)
    circ.barrier()
    circ.x(1)
    # circ.cx(0, 1)
    circ.barrier()
    print("Original circuit:")
    print(circ)
    print("Transpiled circuit:")
    transpiled_circ = my_transpile(circ)
    print(transpiled_circ)
    passed_list.append(my_test(circ, transpiled_circ))
    number_testcases += 1
    print("Optimzed circuit:")
    optimized_circ = my_transpile_optimized(circ)
    print(optimized_circ)
    passed_list.append(my_test(circ, optimized_circ))
    number_testcases += 1

    passed = sum(passed_list)
    print("Passed:", passed, "/", number_testcases)
    return passed


def my_transpile(old_circuit):
    new_circuit = QuantumCircuit(*old_circuit.qregs,
        *old_circuit.cregs, name=old_circuit.name + '_transpiled')

    for inst, qargs, cargs in old_circuit._data:
        # print(inst.name)
        if inst.name == "id":
            # print("my_I removes id")
            my_I(new_circuit, qargs)
        elif inst.name == "x":
            # print("my_x used")
            my_x(new_circuit, qargs)
        elif inst.name == "y":
            # print("my_y used")
            my_y(new_circuit, qargs)
        elif inst.name == "z":
            # print("my_z used")
            my_z(new_circuit, qargs)
        elif inst.name == "h":
            # print("my_h used")
            my_h(new_circuit, qargs)
        elif inst.name == "cx":
            # print("my_cx used")
            my_cx(new_circuit, *qargs)
        elif inst.name == "ry":
            # print("my_ry used")
            my_ry(new_circuit, qargs)
        elif inst.name == "barrier":
            # print("barrier kept")
            new_circuit._append(inst, qargs, cargs)
        else:
            new_circuit._append(inst, qargs, cargs)

    return new_circuit

def append_index_to_but_reset_other(index, target_list, list_set, qubit):
    list_set.discard(target_list)
    print("list_set: ", list_set)
    # cx_list = [[] for i in old_circuit._qubits]
    # x_list = [[] for i in old_circuit._qubits]
    # y_list = [[] for i in old_circuit._qubits]

    target_list[qubit].append(index)


def my_transpile_optimized(old_circuit):
    """
    Implements following optimizations:

    1. Skip I
    2. Deletes RX, RX, RX, RX on same qubit
    3. Deletes RZ, RZ, RZ, RZ on same qubit 
    4. Deletes CZ, CZ on same qubits
    5. Deletes CX, CX on same qubits
    6. Deletes H, H on same qubit
    7. Deletes X, X on same qubit
    8. Deletes Y, Y on same qubit

    More identities could probably used to optimized further...
    """
    optimized_circuit = QuantumCircuit(*old_circuit.qregs,
        *old_circuit.cregs, name=old_circuit.name + '_optimized')
    optimized_circuit_new_set = QuantumCircuit(*old_circuit.qregs,
        *old_circuit.cregs, name=old_circuit.name + '_optimized_new_set')

    # Implement rules 5, 6, 7, 8, as these are on the old set

    # list counts consecutively gates by number of elements for every qubit
    # (qubitwise because gates on different qubits do not cancel:
    # e.g: q1: Rx Rx, q2: Rx Rx, do not cancel
    cx_list = [[] for i in old_circuit._qubits]
    h_list = [[] for i in old_circuit._qubits]
    x_list = [[] for i in old_circuit._qubits]
    y_list = [[] for i in old_circuit._qubits]
    delete_list = [[] for i in old_circuit._qubits]
    # gate index
    i = 0 
    # counts deleted rx, rz, cz gates
    number_deleted_gates_new_set = 0
    for inst, qargs, cargs in old_circuit._data:
        print(inst, qargs)
        print("cx_list:", cx_list)
        # print("qargs:", qargs[0].index)
        if inst.name == "h":
            # clear other gates list as its not consecutive anymore
            cx_list = [[] for i in old_circuit._qubits]
            x_list = [[] for i in old_circuit._qubits]
            y_list = [[] for i in old_circuit._qubits]

            h_list[qargs[0].index].append(i)

        elif inst.name == "x":
            # clear other gates list as its not consecutive anymore
            cx_list = [[] for i in old_circuit._qubits]
            h_list = [[] for i in old_circuit._qubits]
            y_list = [[] for i in old_circuit._qubits]

            x_list[qargs[0].index].append(i)

        elif inst.name == "y":
            # clear other gates list as its not consecutive anymore
            cx_list = [[] for i in old_circuit._qubits]
            h_list = [[] for i in old_circuit._qubits]
            x_list = [[] for i in old_circuit._qubits]

            y_list[qargs[0].index].append(i)

        elif inst.name == "z":
            # WARNING: No z_list as Z Z will be removed trough step 2
            # clear other gates list as its not consecutive anymore
            cx_list = [[] for i in old_circuit._qubits]
            h_list = [[] for i in old_circuit._qubits]
            x_list = [[] for i in old_circuit._qubits]
            y_list = [[] for i in old_circuit._qubits]

#         elif inst.name == "i":
#             # WARNING: No z_list as Z Z will be removed trough step 2
#             # clear other gates list as its not consecutive anymore
#             cx_list = [[] for i in old_circuit._qubits]
#             h_list = [[] for i in old_circuit._qubits]
#             x_list = [[] for i in old_circuit._qubits]
#             y_list = [[] for i in old_circuit._qubits]

        elif inst.name == "cx":
            # clear other gates list as its not consecutive anymore
            h_list = [[] for i in old_circuit._qubits]
            x_list = [[] for i in old_circuit._qubits]
            y_list = [[] for i in old_circuit._qubits]

            # As CZ is a symmetric qubit gate both cz(0, 1)
            # and cz(1, 0) need to be count
            cx_list[qargs[0].index].append(i)
            cx_list[qargs[1].index].append(i)

        # Rule 5: delete if 2 consecutive Cx
        if len(cx_list[qargs[0].index]) == 2:
            delete_list[qargs[0].index].append(cx_list[qargs[0].index])
            number_deleted_gates_new_set += 2 * CXGATES

            # reset to empty list for that qubit
            cx_list = [[] for i in old_circuit._qubits]
         
        # Rule 6: delete if 2 consecutive H
        if len(h_list[qargs[0].index]) == 2:
            delete_list[qargs[0].index].append(h_list[qargs[0].index])
            number_deleted_gates_new_set += 2 * HGATES
            
            # reset to empty list for that qubit
            h_list = [[] for i in old_circuit._qubits]

        # Rule 7: delete if 2 consecutive x
        if len(x_list[qargs[0].index]) == 2:
            delete_list[qargs[0].index].append(x_list[qargs[0].index])
            number_deleted_gates_new_set += 2 * XGATES

            # reset to empty list for that qubit
            x_list = [[] for i in old_circuit._qubits]

        # Rule 8: delete if 2 consecutive y
        if len(y_list[qargs[0].index]) == 2:
            delete_list[qargs[0].index].append(y_list[qargs[0].index])
            number_deleted_gates_new_set += 2 * YGATES

            # reset to empty list for that qubit
            y_list = [[] for i in old_circuit._qubits]

        i += 1
    # flatten the lists for different qubits
    for x, _ in enumerate(old_circuit._qubits):
        delete_list[x] = [element for line in delete_list[x] \
            for element in line]
        print("delete[" + str(x) +"]: ", delete_list[x])

    # flatten the lists for all qubits, as seperation is not needed anymore
    delete_list_flatten = [element for line in delete_list \
        for element in line]
    print("delete_list_flatten: ", delete_list_flatten)

    counter = 0
    for inst, qargs, cargs in old_circuit._data:
        # print(inst.name)
        # if gate should be deleted skip it
        # (do not insert it into optimized circ)
        if counter in delete_list_flatten:
            pass
        elif inst.name == "id":
            print("my_I removes id")
            # my_I(optimized_circuit, qargs)
        elif inst.name == "x":
            # print("my_x used")
            optimized_circuit.x(qargs)
        elif inst.name == "y":
            # print("my_y used")
            optimized_circuit.y(qargs)
        elif inst.name == "z":
            # print("my_z used")
            optimized_circuit.z(qargs)
        elif inst.name == "h":
            # print("my_h used")
            optimized_circuit.h(qargs)
        elif inst.name == "cx":
            # print("my_cx used")
            optimized_circuit.cx(*qargs)
        elif inst.name == "ry":
            # print("my_ry used")
            optimized_circuit.ry(np.pi/2, qargs)
        elif inst.name == "barrier":
            # print("barrier kept")
            optimized_circuit._append(inst, qargs, cargs)
        else:
            optimized_circuit._append(inst, qargs, cargs)

        counter += 1

    number_deleted_gates = sum([len(qubit_list) \
        for qubit_list in delete_list])
    print("Optimization deleted " + str(number_deleted_gates)
          + " gates in first step (gates from old basis set)"
          + " and " + str(number_deleted_gates_new_set)
          + " gates in terms of new set.")

    # Implement rules 1, 2, 3, 4, as these are on the new set

    # list counts consecutively gates by number of elements for every qubit
    # (qubitwise because gates on different qubits do not cancel:
    # e.g: q1: Rx Rx, q2: Rx Rx, do not cancel
    rx_list = [[] for i in old_circuit._qubits]
    rz_list = [[] for i in old_circuit._qubits]
    cz_list = [[] for i in old_circuit._qubits]
    delete_list = [[] for i in old_circuit._qubits]
    # gate index
    i = 0 
    for inst, qargs, cargs in old_circuit._data:
        # print(inst, qargs)
        # print("qargs:", qargs[0].index)
        if inst.name == "rx":
            # clear other gates list as its not consecutive anymore
            rz_list = [[] for i in old_circuit._qubits]
            cz_list = [[] for i in old_circuit._qubits]

            rx_list[qargs[0].index].append(i)

        elif inst.name == "rz":
            # clear other gates list as its not consecutive anymore
            rx_list = [[] for i in old_circuit._qubits]
            cz_list = [[] for i in old_circuit._qubits]

            rz_list[qargs[0].index].append(i)

        elif inst.name == "cz":
            # clear other gates list as its not consecutive anymore
            rx_list = [[] for i in old_circuit._qubits]
            rz_list = [[] for i in old_circuit._qubits]

            # As CZ is a symmetric qubit gate both cz(0, 1)
            # and cz(1, 0) need to be count
            cz_list[qargs[0].index].append(i)
            cz_list[qargs[1].index].append(i)
         
        # Rule 2: delete if 4 consecutive Rx
        if len(rx_list[qargs[0].index]) == 4:
            delete_list[qargs[0].index].append(rx_list[qargs[0].index])
            
            # reset to empty list for that qubit
            rx_list = [[] for i in old_circuit._qubits]

        # Rule 3: delete if 4 consecutive Rz
        if len(rz_list[qargs[0].index]) == 4:
            delete_list[qargs[0].index].append(rz_list[qargs[0].index])

            # reset to empty list for that qubit
            rz_list = [[] for i in old_circuit._qubits]

        # Rule 4: delete if 2 consecutive Cz
        if len(cz_list[qargs[0].index]) == 2:
            delete_list[qargs[0].index].append(cz_list[qargs[0].index])

            # reset to empty list for that qubit
            cz_list = [[] for i in old_circuit._qubits]

        i += 1

    # flatten the lists for different qubits
    for x, _ in enumerate(old_circuit._qubits):
        delete_list[x] = [element for line in delete_list[x] \
            for element in line]
        print("delete[" + str(x) +"]: ", delete_list[x])

    # flatten the lists for all qubits, as seperation is not needed anymore
    delete_list_flatten = [element for line in delete_list \
        for element in line]
    print("delete_list_flatten: ", delete_list_flatten)

    counter = 0
    for inst, qargs, cargs in optimized_circuit._data:
        # print(inst.name)
        # if gate should be deleted skip it
        # (do not insert it into optimized circ)
        if counter in delete_list_flatten:
            pass
        elif inst.name == "id":
            print("my_I removes id")
            # my_I(optimized_circuit_new_set, qargs)
        elif inst.name == "x":
            # print("my_x used")
            my_x(optimized_circuit_new_set, qargs)
        elif inst.name == "y":
            # print("my_y used")
            my_y(optimized_circuit_new_set, qargs)
        elif inst.name == "z":
            # print("my_z used")
            my_z(optimized_circuit_new_set, qargs)
        elif inst.name == "h":
            # print("my_h used")
            my_h(optimized_circuit_new_set, qargs)
        elif inst.name == "cx":
            # print("my_cx used")
            my_cx(optimized_circuit_new_set, *qargs)
        elif inst.name == "ry":
            # print("my_ry used")
            my_ry(optimized_circuit_new_set, qargs)
        elif inst.name == "barrier":
            # print("barrier kept")
            optimized_circuit_new_set._append(inst, qargs, cargs)
        else:
            optimized_circuit_new_set._append(inst, qargs, cargs)

        counter += 1

    number_deleted_gates = sum([len(qubit_list) \
        for qubit_list in delete_list])
    print("Optimization deleted " + str(number_deleted_gates)
          + " gates in second step (gates from new basis set).")

    return optimized_circuit_new_set


# Use flags to turn on/off parts of calculations
print_basis = True
search_i = False
search_x = False
search_y = True
search_z = False
search_ry = False
search_cx = False

if print_basis:
    # Quick look at the given basis (print matrix)
    rx_ideal = QuantumCircuit(1)
    rx = QuantumCircuit(1)
    rx_ideal.rx(np.pi/2, 0)
    print("RX-Gate")
    print(do_job_simulate(rx_ideal, unitary_sim=True))

    rz_ideal = QuantumCircuit(1)
    rz = QuantumCircuit(1)
    rz_ideal.rz(np.pi/2, 0)
    print("RZ-Gate")
    print(do_job_simulate(rz_ideal, unitary_sim=True))

    cz_ideal = QuantumCircuit(2)
    cz = QuantumCircuit(2)
    cz_ideal.cz(0, 1)
    print("CZ-Gate")
    print(do_job_simulate(cz_ideal, unitary_sim=True))

if search_i:
    # example how to get the Z gate
    basis_set = set([ideal_rx, ideal_rz])
    target_set = set([ideal_i])
    brute_force_one_qubit_gate4(basis_set, target_set)

if search_x:
    # example how to get the Z gate
    basis_set = set([ideal_rx, ideal_rz])
    target_set = set([ideal_x])
    brute_force_one_qubit_gate4(basis_set, target_set)

if search_y:
    # example how to get the Z gate
    basis_set = set([ideal_rx, ideal_rz])
    target_set = set([ideal_y])
    brute_force_one_qubit_gate4(basis_set, target_set)

if search_z:
    # example how to get the Z gate
    basis_set = set([my_h, ideal_rx, ideal_rz])
    basis_set = set([ideal_rx, ideal_rz])
    target_set = set([ideal_z])
    brute_force_one_qubit_gate2(basis_set, target_set)

if search_ry:
    # example how to get the Z gate
    # as h is build by 3 gates we use it to brute force for 5 gates
    basis_set = set([my_h, ideal_z, ideal_rx, ideal_rz])
    target_set = set([ideal_ry])
    brute_force_one_qubit_gate3(basis_set, target_set)

if search_cx:
    # get the two qubit gate cx
    basis_set = set([my_h, ideal_cz, ideal_rx, ideal_rz])
    target_set = set([ideal_cx])
    brute_force_two_qubit_gate3(basis_set, target_set)

    # all qubit until 6 components were tried. As h has 3,
    # the solution H CZ H has 7 and is the smallst
    # basis_set = set([ideal_cz, ideal_rx, ideal_rz])
    # target_set = set([ideal_cx])
    # brute_force_two_qubit_gate6(basis_set, target_set)

testcases()

