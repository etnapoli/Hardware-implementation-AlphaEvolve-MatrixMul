# Scalar Version of the AlphaEvolve Algorithm for 4x4 Complex-Valued Matrix Multiplication  

This directory contains the scalar implementation of the AlphaEvolve algorithm as described in the paper **CITE THE PAPER TITLE**.

## Overview

The AlphaEvolve algorithm performs 4x4 complex matrix multiplication using an optimized combination of Frobenius expansions and inner products. The scalar implementation explicitly computes all intermediate complex terms, which are later used in the hardware implementation.  

The algorithm can be divided into two main steps:

1. **Frobenius Expansion**  
   Computes the intermediate complex terms a0 ... a47 and b0 ... b47 as linear combinations of the input matrices A and B.  

2. **Inner Product Combination**  
   Uses the results of the multiplications m0 ... m47 to compute the final output matrix C.  

All the terms are complex numbers, and the detailed equations can be found in:
- [AlphaEvolve Colab](#)  
- \cite{phialsbasement_alphaevolve_matrixmul_verification_2025}

## Implementation Details

- **Doubling Optimization:**  
  All equations are doubled to remove the factor 1/2. A division step at the end recovers the correct scaling.  

- **Real and Imaginary Expansion:**  
  Each complex equation is expanded into separate real and imaginary parts, leading to a total of 224 multi-operand addition/subtraction expressions.  

- **Examples:**  
  The real and imaginary parts of a0 are explicitly expanded in the code. The complete set of expanded expressions is available in \cite{Napoli_alphaevolve_matrixmul_circuit_2025}.

## Contents

- `expressionsA.txt` — Text file. Contains the equations for the real and imaginary parts of the a0 ... a47 variables obtained by the Frobenius product:
a_r = <u^(r), A>

- `expressionsB.txt` — Text file. Contains the equations for the real and imaginary parts of the b0 ... b47 variables obtained by the Frobenius product:
b_r = <v^(r), B>
- `expressionsC.txt` — Text file. Contains the equations for the real and imaginary parts of the C1 ... C16 variables obtained by the inner product:
C = m1*w^(1) + m2*w^(2) + ... + m48*w^(48)
 - This version is not yet optimized for hardware implementation. In particular, common subexpressions have not been isolated to reduce the hardware footprint.

## License

This code is released under the **MIT License**. See the `LICENSE` file for details.  

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

- The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.  
- THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


