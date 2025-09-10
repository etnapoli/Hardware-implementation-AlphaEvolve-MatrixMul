// SPDX-License-Identifier: CERN-OHL-NC-2.0
// Author: Ettore Napoli
// Affiliation: University of Salerno
// August 2025
// Description: AlphaEvolve algorithm for matrix multiplication

module matrix_mult_4x4_complex_alphaevolve #(parameter w = 16)(
    input  signed [w-1:0] A_real [0:3][0:3],
    input  signed [w-1:0] A_imag [0:3][0:3],
    input  signed [w-1:0] B_real [0:3][0:3],
    input  signed [w-1:0] B_imag [0:3][0:3],
    output signed [2*w+2:0] C_real [0:3][0:3], //2*w +3 bits
    output signed [2*w+2:0] C_imag [0:3][0:3]  //2*w +3 bits
);


/*==================================================================
The content of this file referes to a submitted scientific paper
The actual code will be released once the paper has been reviewed
=======================================================/*


// assign real part of the outputs after dividing internal signal by 8
assign C_real[0][0]  = (C1_r>>>3);
assign C_real[0][1]  = (C2_r>>>3);
assign C_real[0][2]  = (C3_r>>>3);
assign C_real[0][3]  = (C4_r>>>3);
assign C_real[1][0]  = (C5_r>>>3);
assign C_real[1][1]  = (C6_r>>>3);
assign C_real[1][2]  = (C7_r>>>3);
assign C_real[1][3]  = (C8_r>>>3);
assign C_real[2][0]  = (C9_r>>>3);
assign C_real[2][1]  = (C10_r>>>3);
assign C_real[2][2]  = (C11_r>>>3);
assign C_real[2][3]  = (C12_r>>>3);
assign C_real[3][0]  = (C13_r>>>3);
assign C_real[3][1]  = (C14_r>>>3);
assign C_real[3][2]  = (C15_r>>>3);
assign C_real[3][3]  = (C16_r>>>3);

// assign real part of the outputs after dividing internal signal by 8
assign C_imag[0][0]  = (C1_i>>>3);
assign C_imag[0][1]  = (C2_i>>>3);
assign C_imag[0][2]  = (C3_i>>>3);
assign C_imag[0][3]  = (C4_i>>>3);
assign C_imag[1][0]  = (C5_i>>>3);
assign C_imag[1][1]  = (C6_i>>>3);
assign C_imag[1][2]  = (C7_i>>>3);
assign C_imag[1][3]  = (C8_i>>>3);
assign C_imag[2][0]  = (C9_i>>>3);
assign C_imag[2][1]  = (C10_i>>>3);
assign C_imag[2][2]  = (C11_i>>>3);
assign C_imag[2][3]  = (C12_i>>>3);
assign C_imag[3][0]  = (C13_i>>>3);
assign C_imag[3][1]  = (C14_i>>>3);
assign C_imag[3][2]  = (C15_i>>>3);
assign C_imag[3][3]  = (C16_i>>>3);


endmodule

