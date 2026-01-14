# HDL Directory

This folder contains the SystemVerilog source code for the implementation of the AlphaEvolve  algorithm for 4x4 complex valued matrix multiplication.

## Contents
- `matrix_mult_4x4_complex_alphaevolve.sv` — SystemVerilog description of the proposed implementation of the AlphaEvolve algorithm for 4×4 complex-valued matrix multiplication.  
  The code is parameterized on **w**, which is the number of bits used for the 2's complement representation of the input matrix elements.  

- `tb_matrix_mult_4x4_complex.sv` — testbench for simulation.  
  Tests have been run for **w** values of 4, 8, 16, 24, and 32.  
  Please note that the provided testbench expects to find the test vector files **two directories above its location** (`../../`).  
  If this is not the case, correct line 50 of the testbench accordingly.     

- `README.md` — this file.  

## Notes
- Testbenches have been verified with:
  - **Questa Intel Starter FPGA Edition-64 2021.2** (bundled with **Quartus Prime 22.1 Lite**)  
  - **Cadence Xcelium Simulator 20.09-s001** (bundled with Cadence)

## License

This code is released under the **MIT License**. See the `LICENSE` file for details.  

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

- The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.  
- THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
