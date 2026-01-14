// SPDX-License-Identifier: CERN-OHL-NC-2.0
// Author: Ettore Napoli
// Affiliation: University of Salerno
// january 2026
// Description: Strassen algorithm for 4x4 complex matrix multiplication

module matrix_mult_4x4_complex_strassen #(parameter w = 32)(
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

    // Combine result into 4x4 matrix C
	 // real part
	 assign C_real[0][0] = C11r[0];
	 assign C_real[0][1] = C11r[1];
	 assign C_real[1][0] = C11r[2];
	 assign C_real[1][1] = C11r[3];

	 assign C_real[0][2] = C12r[0];
	 assign C_real[0][3] = C12r[1];
	 assign C_real[1][2] = C12r[2];
	 assign C_real[1][3] = C12r[3];

	 assign C_real[2][0] = C21r[0];
	 assign C_real[2][1] = C21r[1];
	 assign C_real[3][0] = C21r[2];
	 assign C_real[3][1] = C21r[3];

	 assign C_real[2][2] = C22r[0];
	 assign C_real[2][3] = C22r[1];
	 assign C_real[3][2] = C22r[2];
	 assign C_real[3][3] = C22r[3];

	 // imag part
	 assign C_imag[0][0] = C11i[0];
	 assign C_imag[0][1] = C11i[1];
	 assign C_imag[1][0] = C11i[2];
	 assign C_imag[1][1] = C11i[3];

	 assign C_imag[0][2] = C12i[0];
	 assign C_imag[0][3] = C12i[1];
	 assign C_imag[1][2] = C12i[2];
	 assign C_imag[1][3] = C12i[3];

	 assign C_imag[2][0] = C21i[0];
	 assign C_imag[2][1] = C21i[1];
	 assign C_imag[3][0] = C21i[2];
	 assign C_imag[3][1] = C21i[3];

	 assign C_imag[2][2] = C22i[0];
	 assign C_imag[2][3] = C22i[1];
	 assign C_imag[3][2] = C22i[2];
	 assign C_imag[3][3] = C22i[3];
endmodule



