# Scalar Version of the AlphaEvolve Algorithm for 4x4 Complex-Valued Matrix Multiplication  

This directory contains the scalar implementation of the AlphaEvolve algorithm as described in the paper and the Colab **[1], [2]**.  
The MATLAB version is derived from the Python version published by F. Bakreski **[3]**.

## Overview

The AlphaEvolve algorithm performs 4x4 complex matrix multiplication using an optimized combination of Frobenius expansions and inner products.  
The scalar implementation explicitly computes all intermediate complex terms, which are later used in the hardware implementation.  

The algorithm can be divided into two main steps:

1. **Frobenius Expansion**  
   Computes the intermediate complex terms a0 ... a47 and b0 ... b47 as linear combinations of the input matrices A and B.  

2. **Inner Product Combination**  
   Uses the results of the multiplications m0 ... m47 to compute the final output matrix C.  

All the terms are complex numbers, and the detailed equations can be found in **[2], [3]**:  
- [AlphaEvolve Colab](https://arxiv.org/abs/2506.13131)  
- [PhialsBasement AlphaEvolve MatrixMul Verification](https://github.com/PhialsBasement/AlphaEvolve-MatrixMul-Verification)

## Implementation Details

- **Doubling Optimization:**  
  All equations are doubled to remove the factor 1/2. A division step at the end recovers the correct scaling.  

- **Real and Imaginary Expansion:**  
  Each complex equation is expanded into separate real and imaginary parts, leading to a total of 224 multi-operand addition/subtraction expressions.  

- **Examples:**  
  The real and imaginary parts of a0 are explicitly expanded in the code. The complete set of expanded expressions is available in Napoli’s reference on AlphaEvolve circuit design.

## Contents

- `expressionsA.txt` — Contains the equations for the real and imaginary parts of a0 ... a47 obtained by the Frobenius product:  
  a_r = <u^(r), A>  

- `expressionsB.txt` — Contains the equations for the real and imaginary parts of b0 ... b47 obtained by the Frobenius product:  
  b_r = <v^(r), B>  

- `expressionsC.txt` — Contains the equations for the real and imaginary parts of C1 ... C16 obtained by the inner product:  
  C = m1·w^(1) + m2·w^(2) + ... + m48·w^(48)  

- `README.md` — this file.  

- This version is not yet optimized for hardware implementation. In particular, common subexpressions have not been isolated to reduce the hardware footprint.

## References

**[1]** A. Novikov, N. Vũ, M. Eisenberger, E. Dupont, P.-S. Huang, A. Z. Wagner, S. Shirobokov, B. Kozlovskii, F. J. R. Ruiz, A. Mehrabian, M. P. Kumar, A. See, S. Chaudhuri, G. Holland, A. Davies, S. Nowozin, P. Kohli, and M. Balog, “AlphaEvolve: A coding agent for scientific and algorithmic discovery,” 2025. [Online]. Available: https://arxiv.org/abs/2506.13131  

**[2]** Google DeepMind, “AlphaEvolve: Mathematical results,” 2025. [Online]. Available: https://colab.research.google.com/github/google-deepmind/alphaevolve_results/blob/master/mathematical_results.ipynb  

**[3]** F. Bakreski, “AlphaEvolve-MatrixMul-Verification,” GitHub repository, 2025. Available: https://github.com/PhialsBasement/AlphaEvolve-MatrixMul-Verification  

## License

This code is released under the **MIT License**. See the `LICENSE` file for details.  

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

- The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.  
- THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

