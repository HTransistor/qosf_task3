# qosf_task3

##Task 3  
Please write a simple compiler – program, which translates one quantum  
circuit into another, using a restricted set of gates.  

You need to consider just the basic gates for the input circuit,  
such as (I, H, X, Y, Z, RX, RY, RZ, CNOT, CZ).  

The output circuit should consist only from the following gates:  
RX, RZ, CZ. In other words, each gate in the original circuit must be  
replaced by an equivalent combination of gates coming from the restricted  
set (RX, RZ, CZ) only.  

For example, a Hadamard gate after compilation looks like this:  
RZ(pi/2)  
RX(pi/2)  
RZ(pi/2)  

Analyze what’s the overhead of the compiled program compared to the original  
one and propose how to improve it. What we mean by overhead is the  
following: by replacing all the initial gates with the restricted set  
of gates given in the problem, you will see that the resulting circuit  
is much more involved than the original one. This is what we called  
the overhead, and you may think about how to treat this problem,  
i.e. you could try to simplify as much as possible the resulting circuit.  


## Results  
Analyzing the overhead:  
The overhead depends on the basis set to which one we compare it,  
as every gate of the basis set needs to be rebuilt by the new set.  
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
(CX ^ 7) * (I ^ 4) * (X ^ 4) * (Y ^ 4) * (H ^ 4) * (RY ^ 3) * (Z ^ 4)  
When optimizing the circuited I can be removed (set to I=0 in above formula)  

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

In the optimized version the overhead is hardly depending on the circuit  
and the order of the used gates. There are a lot of combination that reduce  
the gate number, some are:  

CxCx = -14  
XX = -8, YY = -8, XZX = -8  
HH = -6  
I = -4, ZZ = -4, CxRy = -4, HZH = -4, HZX = -4, XZH = -4, XRy = -4  
HRy = -4, XRy = -4, YRy = -4, ZRy = -4, YZ = -4  
RxRxRxRx = -4,  RzRzRzRz = -4  
