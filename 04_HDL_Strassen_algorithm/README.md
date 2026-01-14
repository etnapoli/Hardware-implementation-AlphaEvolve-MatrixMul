# HDL Directory

This folder contains the SystemVerilog source code for the implementation of the Strassen algorithm for 4x4 complex valued matrix multiplication.

## Contents
- `matrix_mult_4x4_complex_strassen.sv` — SystemVerilog description of the reference implementation of the Strassen algorithm for 4×4 complex-valued matrix multiplication.  
  The code is parameterized on **W**, which is the number of bits used for the 2's complement representation of the input matrix elements.  

- `strassen_2x2_complex.sv` — SystemVerilog description of 2x2 complex valued matrix multiplication using the Strassen algorithm. 
  The circuit is  exploited by `matrix_mult_4x4_complex_strassen.sv` to implement the 4x4 matrix multiplication. 
  The code is parameterized on **W**, which is the number of bits used for the 2's complement representation of the input matrix elements.  

- `mult_complex.sv` — SystemVerilog description of a binary multiplier with complex operands. 
  The circuit is exploited by `strassen_2x2_complex.sv`. 
  The code is parameterized on **W**, which is the number of bits used for the 2's complement representation of the input matrix elements.  

- `tb_matrix_mult_4x4_complex.sv` — testbench for simulation.  
  Tests have been run for **W** values of 4, 8, 16, 24, 32, and 64.  
  Please note that the provided testbench expects to find the test vector files **two directories above its location** (`../../`).  
  If this is not the case, correct line 50 of the testbench accordingly.     

- `Test_Vectors_files/` — directory containing the test vector files for the circuit.  
  Five files are provided, corresponding to different input bit-widths:  
  - `golden_values_04bit_ver_3300_test_vectors.txt` — 4-bit input elements  
  - `golden_values_08bit_ver_3300_test_vectors.txt` — 8-bit input elements  
  - `golden_values_16bit_ver_3300_test_vectors.txt` — 16-bit input elements  
  - `golden_values_24bit_ver_3300_test_vectors.txt` — 24-bit input elements  
  - `golden_values_32bit_ver_3300_test_vectors.txt` — 32-bit input elements  
  - `golden_values_64bit_ver_3300_test_vectors.txt` — 64-bit input elements  

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
