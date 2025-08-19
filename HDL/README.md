# HDL Directory

This folder contains the Verilog/SystemVerilog source code for the project.

## Contents
- `matrix_mult_4x4_complex_alphaevolve.sv` — SystemVerilog description of the proposed implementation of the AlphaEvolve algorithm for 4×4 complex-valued matrix multiplication.  
  The code is parameterized on **w**, which is the number of bits used for the 2's complement representation of the input matrix elements.  

- `tb_matrix_mult_4x4_complex.sv` — testbench for simulation.  
  Tests have been run for **w** values of 4, 8, 16, 24, and 32.  
  Please note that the provided testbench expects to find the test vector files **two directories above its location** (`../../`).  
  If this is not the case, correct line 50 of the testbench accordingly.     

- `Test_Vectors_files/` — directory containing the test vector files for the circuit.  
  Five files are provided, corresponding to different input bit-widths:  
  - `golden_values_04bit_ver_3300_test_vectors.txt` — 4-bit input elements  
  - `golden_values_08bit_ver_3300_test_vectors.txt` — 8-bit input elements  
  - `golden_values_16bit_ver_3300_test_vectors.txt` — 16-bit input elements  
  - `golden_values_24bit_ver_3300_test_vectors.txt` — 24-bit input elements  
  - `golden_values_32bit_ver_3300_test_vectors.txt` — 32-bit input elements  

## Notes
- Testbenches have been verified with:
  - **Questa Intel Starter FPGA Edition-64 2021.2** (bundled with **Quartus Prime 22.1 Lite**)  
  - **TBD** (Cadence tool version)
