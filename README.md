# Final Year Project (Bachelor's of Science in Physics) Thesis: "Using Higher Order Product Formulas in the Double-Bracket Iteration Quantum Algorithm"

This repository is dedicated to my honor's year thesis at Nanyang Technological University Singapore. My supervisors Dr. Marek Gluza and Dr. Koh Teck Seng have helped me make this project come into fruition. 
This work is a direct extension to my supervisor's original work on the double-bracket iteration quantum algorithm: https://arxiv.org/abs/2206.11772

### How It's Made:

- Code is written in Python following PEP8 standards, with emphasis on modular functions and readability.

Double-bracket iteration quantum algorithms (DBI) are a class of quantum algorithms used primarily for tasks like Hamiltonian diagonalization and ground-state preparation. They are inspired by the classical double-bracket flow, a mathematical method for diagonalizing matrices. 

In this project, I analyzed adding higher-order group commutator terms into the DBI algorithm. We did this to test whether incorporating these terms increases the diagonalization accuracy and efficiency of the algorithm. What we found is that the accuracy significantly improves, with the error dropping from $\approx9.98\%$ for the second-order formula (S2) to $\approx2.23\%$ for the third-order formula (S3) after four repetitions. This is important because it allows us to achieve a more precise diagonalization using fewer recursion steps, reducing overall circuit depth and making the algorithm more suitable for NISQ-era quantum devices.

The evolution generator is defined as follows: 

$\hat{W}(0) = [\hat{D}(0),\hat{H}(0)]$

and the resulting Hamiltonian is given as follows: 

$\hat{H}(s) = e^{s\hat{W}}\hat{H}(0)e^{-s\hat{W}}$


The cost function is defined as follows:

$\|\sigma(\hat{H}(s))\|^2_{HS} = tr(\sigma^2(\hat{H}(s)))$ 

We investigate the optimisation of this cost function over an arbitrarily chosen flow duration (give by $s$). $\sigma$ is the off-diagonal terms of the Hamiltonian.



### Results

We compared the performance of second-order (S2) and third-order (S3) product formulas within the Double-Bracket Iteration (DBI) framework.

- **Error Scaling:**  
  - S2 achieves error scaling of $O(s^3)$.  
  - S3 improves this to $\mathcal{O}(s^4)$, reducing approximation error by one order.  

- **Numerical Performance (Transverse-Field Ising Model, L=5):**  
  - After 14 recursion steps:  
    - DBI baseline: $\|\sigma\|$ $\approx$ 2.10  
    - S2: $\|\sigma\|$ $\approx$ 10.71  
    - S3: $\|\sigma\|$ $\approx$ 5.89

 
S3 clearly outperforms S2 in diagonalization efficiency.  

- **Repetition Improves Accuracy:**  
  - Repeating S3 four times approximates the exact commutator exponential within **2.23% error**.  
  - By contrast, four repetitions of S2 yield ~9.98% error.  

- **Trade-off:**  
  S3 requires six exponentials (vs. four for S2), leading to a higher gate cost. However, the improved accuracy and diagonalization power make S3 more favorable for Noisy-Intermediate Scale Quantum (NISQ) devices.
  
### Additional: Quantum Gates and Circuit Optimization

The Double-Bracket Iteration (DBI) algorithm and its higher-order product formula variants (S2, S3) can be implemented using standard quantum gate sets.

- **Gate Set Used**  
  - **Single-qubit rotations**: $(R_x(\theta), R_y(\theta), R_z(\theta)$
  - **Two-qubit entangling gates**: CNOT or CZ (to capture nearest-neighbor interactions, as in the Ising model)  
  - **Phase gates** ($S$, $T$) as needed for universality and decomposition  

Each product formula step requires exponentials of the Hamiltonian $e^{-iHt}$ and diagonal terms $e^{-iDt}$, which are decomposed into sequences of these native gates.

- **Circuit Optimization Benefits**  
  - **Fewer steps needed**: Higher-order product formulas improve error scaling  
    - S2 scales as $\mathcal{O}(s^3)$  
    - S3 improves this to $\mathcal{O}(s^4)$  
  
This means fewer recursion steps (and gates) are required for the same accuracy.  


  - **Shallower circuits**: Although S3 uses six exponentials per step (vs. four for S2), its improved accuracy reduces total depth, which is critical for NISQ-era devices.  
  - **Hardware efficiency**: The Hamiltonian terms (e.g., $X_jX_{j+1}, Z_j$) map directly to short gate sequences of rotations + CNOTs, reducing compilation overhead.  
  - **Structured ansatz**: Unlike variational methods with random parameters, DBI uses commutator flows, which are less prone to barren plateaus and more resource-efficient.  

In short, the DBI framework with higher-order product formulas provides a **trade-off between accuracy and gate cost**, allowing for shallower, more hardware-friendly quantum circuits.





### Future Directions

- Optimizing Product Formula Implementation: Further improvements could focus on optimizing the S3 (third-order) product formula for speed and hardware efficiency. For example, refactoring the implementation could reduce runtime by combining commuting exponentials and minimizing redundant operations.
- Scaling to Larger Systems: Testing DBI with higher-order product formulas on larger qubit systems would help understand performance trade-offs in NISQ devices and evaluate how error scaling behaves with system size.
- Hardware-Specific Compilation: Tailoring the transpilation to specific hardware architectures could reduce gate count and circuit depth, further enhancing practicality for real quantum devices.
- Hybrid Approaches: Combining DBI with variational techniques or other ansatzes might improve diagonalization efficiency while retaining algorithmic structure.

### Lessons Learned:

- Impact of Higher-Order Terms: Adding higher-order product formulas (S2 → S3) significantly improves diagonalization accuracy, as evidenced by error reduction from ~9.98% to 2.23% after four repetitions. This demonstrates that carefully chosen algorithmic refinements can outweigh brute-force or naive approaches.
- Trade-Off Between Accuracy and Gate Cost: While S3 uses more exponentials per step (6 vs. 4 for S2), the improved error scaling allows for shallower circuits overall. This highlights the importance of balancing theoretical accuracy with practical resource constraints.
- Algorithmic vs Variational Approaches: DBI is a deterministic, algorithmic flow; unlike variational circuits, it does not rely on parameter optimization. This makes it tolerant to overfitting on small systems and predictable in performance as system size grows
