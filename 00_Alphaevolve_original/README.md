# Matlab Version of the AlphaEvolve Algorithm for 4x4 Complex-Valued Matrix Multiplication  

This directory contains the matlab implementation of the AlphaEvolve algorithm as described in the paper and the Colab **[1], [2]**.  
The MATLAB version is derived from the Python version published by F. Bakreski **[3]**.

## Overview

The AlphaEvolve algorithm performs 4x4 complex matrix multiplication using an optimized combination of Frobenius expansions and inner products.  
The matlab code below implements the algorithm using complex arithmetic and provides a test bench for the test.

All the terms are complex numbers, and the detailed equations can be found in **[2], [3]**:  
- [AlphaEvolve Colab](https://arxiv.org/abs/2506.13131)  
- [PhialsBasement AlphaEvolve MatrixMul Verification](https://github.com/PhialsBasement/AlphaEvolve-MatrixMul-Verification)

## Contents

- `alphaevolve_4x4_complex.m`  
  The matlab function that implements the algorithm.
  
- `test_alphaevolve.m`  
  The matlab scritp that tests the algorithm feeding random inputs.

- `README.md` — this file.  

## Running in matlab
- **Modify the script:**  
  Set the tolerance on the numerical pricions nad the number of random tests.

- **Run the script**  
  At the matlab command prompt, while being the in directory where both alphaevolve_4x4_complex.m and test_alphaevolve.m files are located run:
  test_alphaevolve


## References

**[1]** A. Novikov, N. Vũ, M. Eisenberger, E. Dupont, P.-S. Huang, A. Z. Wagner, S. Shirobokov, B. Kozlovskii, F. J. R. Ruiz, A. Mehrabian, M. P. Kumar, A. See, S. Chaudhuri, G. Holland, A. Davies, S. Nowozin, P. Kohli, and M. Balog, “AlphaEvolve: A coding agent for scientific and algorithmic discovery,” 2025. [Online]. Available: https://arxiv.org/abs/2506.13131  

**[2]** Google DeepMind, “AlphaEvolve: Mathematical results,” 2025. [Online]. Available: https://colab.research.google.com/github/google-deepmind/alphaevolve_results/blob/master/mathematical_results.ipynb  

**[3]** F. Bakreski, “AlphaEvolve-MatrixMul-Verification,” GitHub repository, 2025. Available: https://github.com/PhialsBasement/AlphaEvolve-MatrixMul-Verification  

## License

This code is released under the **MIT License**. See the `LICENSE` file for details.  

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

- The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.  
- THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

