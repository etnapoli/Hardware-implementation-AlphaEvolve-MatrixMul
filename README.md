# Hardware-Implementation-AlphaEvolve-MatrixMul

Verilog code and MATLAB scripts for the hardware implementation of the **4×4 Rank-48 AlphaEvolve matrix multiplication over the 0.5C algorithm**.  

## Contents

- `00_Alphaevolve_original/` — MATLAB implementation of the AlphaEvolve 4×4 complex-valued matrix multiplication algorithm.  

- `01_Alphaevolve_scalar_version/` — Text files with the scalar derivations of the equations reported in `00_Alphaevolve_original`.  

- `02_Common_Subexpression_Elimination/` — Hardware-optimized equations derived from `01_Alphaevolve_scalar_version` using a **Common Subexpression Elimination (CSE)** algorithm.  

- `03_HDL_AlpohaEvolve_algorithm/` — SystemVerilog source files and testbenches for the complete matrix multiplication circuit.  

- `04_HDL_naive_algorithm/` — SystemVerilog source files and testbenches for the implementation of the naive algorithm for matrix multiplication used as a base test for the AlphaEvolve implementation.  

- `LICENSE` — License file for the entire repository.  

- `README.md` — This file.  

## License

This project is licensed under the **CERN Open Hardware Licence v2 – Non-commercial (CERN-OHL-NC-2.0)**.  
You may use, modify, and share this work for personal, academic, or research purposes.  
Commercial use is prohibited without explicit permission from the author.  

- THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## How to Cite

If you use this repository in your work, please cite it as:  

- Napoli, E. *Hardware Implementation of the AlphaEvolve 4×4 Matrix Multiplication Algorithm*. GitHub, 2025. Available at: [https://github.com/etnapoli/Hardware-implementation-AlphaEvolve-MatrixMul](https://github.com/etnapoli/Hardware-implementation-AlphaEvolve-MatrixMul)  

