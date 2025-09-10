# AlphaEvolve — Common Subexpression Elimination (CSE) Results

This directory contains the optimized expressions for the scalar implementation of the **AlphaEvolve** algorithm after applying **Common Subexpression Elimination (CSE)**.  

## Overview

The original scalar equations expand the real and imaginary parts of all intermediate terms (\(a_0, \ldots, a_{47}\)) into **224 multi-operand addition/subtraction expressions**.  
While these equations are directly suitable for hardware mapping, they contain a large number of redundant subexpressions (identical operand combinations repeated across multiple equations).  

To reduce this redundancy, a **CSE algorithm** has been applied. The goal is to minimize the number of additions and subtractions required for hardware implementation, thereby reducing area and power without changing functionality.

## Common Subexpression Elimination (CSE)

CSE identifies repeated subexpressions across the equations and replaces them with new temporary variables.  
The process is iterative:

1. Start from the initial set of multi-operand expressions.  
2. Search for operand pairs or groups that appear multiple times.  
3. Select one of the most frequently occurring subexpressions.  
4. Define a new variable for this subexpression.  
5. Substitute the variable back into the original equations.  
6. Repeat until no subexpression occurs more than once.  

The detailed algorithm is described in the paper and summarized in Algorithm 1 (*Common Subexpression Elimination*).  

## Contents

- `expressionsA_CSE_num_op_412.txt` — equations after CSE has been applied to the expressions for the real and the imaginary parts of 'a0, ..., a47'.  The number 412 indicates that the version implements the equations using a total of 412 add/sub operations.  

- `expressionsB_CSE_num_op_291.txt` — equations after CSE has been applied to the expressions for the real and the imaginary parts of 'b0, ..., b47'.  The number 291 indicates that the version implements the equations using a total of 291 add/sub operations.  

- `expressionsC _CSE_num_op_493.txt` — equations after CSE has been applied to the expressions for the real and the imaginary parts of 'C1, ..., C16'.  The number 493 indicates that the version implements the equations using a total of 493 add/sub operations.  

- `README.md` — this file.  

## Notes

- The optimized expressions are **functionally equivalent** to the scalar version but require fewer additions/subtractions.  
- These results are intended for **hardware implementation** (e.g., Verilog/VHDL or HLS C) where efficiency is critical.  

## License

This code is released under the **MIT License**. See the `LICENSE` file for details.  

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

- The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.  
- THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
