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

	// 2x2 submatrices 
   wire signed [w-1:0] A11r[0:3], A12r[0:3], A21r[0:3], A22r[0:3];
   wire signed [w-1:0] B11r[0:3], B12r[0:3], B21r[0:3], B22r[0:3];
   wire signed [w-1:0] A11i[0:3], A12i[0:3], A21i[0:3], A22i[0:3];
   wire signed [w-1:0] B11i[0:3], B12i[0:3], B21i[0:3], B22i[0:3];
	 
	// 2x2 submatrices extended to W+1 bits
	wire signed [w:0] B11r_ext[0:3], A11r_ext[0:3], A22r_ext[0:3], B22r_ext[0:3];
	wire signed [w:0] B11i_ext[0:3], A11i_ext[0:3], A22i_ext[0:3], B22i_ext[0:3];
   
	assign {A11r[0], A11r[1], A11r[2], A11r[3]} = {A_real[0][0], A_real[0][1], A_real[1][0], A_real[1][1]};
   assign {A12r[0], A12r[1], A12r[2], A12r[3]} = {A_real[0][2], A_real[0][3], A_real[1][2], A_real[1][3]};
   assign {A21r[0], A21r[1], A21r[2], A21r[3]} = {A_real[2][0], A_real[2][1], A_real[3][0], A_real[3][1]};
   assign {A22r[0], A22r[1], A22r[2], A22r[3]} = {A_real[2][2], A_real[2][3], A_real[3][2], A_real[3][3]};

	assign {B11r[0], B11r[1], B11r[2], B11r[3]} = {B_real[0][0], B_real[0][1], B_real[1][0], B_real[1][1]};
   assign {B12r[0], B12r[1], B12r[2], B12r[3]} = {B_real[0][2], B_real[0][3], B_real[1][2], B_real[1][3]};
   assign {B21r[0], B21r[1], B21r[2], B21r[3]} = {B_real[2][0], B_real[2][1], B_real[3][0], B_real[3][1]};
   assign {B22r[0], B22r[1], B22r[2], B22r[3]} = {B_real[2][2], B_real[2][3], B_real[3][2], B_real[3][3]};
	 
	assign {A11i[0], A11i[1], A11i[2], A11i[3]} = {A_imag[0][0], A_imag[0][1], A_imag[1][0], A_imag[1][1]};
   assign {A12i[0], A12i[1], A12i[2], A12i[3]} = {A_imag[0][2], A_imag[0][3], A_imag[1][2], A_imag[1][3]};
   assign {A21i[0], A21i[1], A21i[2], A21i[3]} = {A_imag[2][0], A_imag[2][1], A_imag[3][0], A_imag[3][1]};
   assign {A22i[0], A22i[1], A22i[2], A22i[3]} = {A_imag[2][2], A_imag[2][3], A_imag[3][2], A_imag[3][3]};

   assign {B11i[0], B11i[1], B11i[2], B11i[3]} = {B_imag[0][0], B_imag[0][1], B_imag[1][0], B_imag[1][1]};
   assign {B12i[0], B12i[1], B12i[2], B12i[3]} = {B_imag[0][2], B_imag[0][3], B_imag[1][2], B_imag[1][3]};
   assign {B21i[0], B21i[1], B21i[2], B21i[3]} = {B_imag[2][0], B_imag[2][1], B_imag[3][0], B_imag[3][1]};
   assign {B22i[0], B22i[1], B22i[2], B22i[3]} = {B_imag[2][2], B_imag[2][3], B_imag[3][2], B_imag[3][3]};

   // Temporary results for adds/subs    w+1 bits
   wire signed [w:0] S1r[0:3], S2r[0:3], S3r[0:3], S4r[0:3], S5r[0:3];
   wire signed [w:0] S1i[0:3], S2i[0:3], S3i[0:3], S4i[0:3], S5i[0:3];
   wire signed [w:0] T1r[0:3], T2r[0:3], T3r[0:3], T4r[0:3], T5r[0:3];
   wire signed [w:0] T1i[0:3], T2i[0:3], T3i[0:3], T4i[0:3], T5i[0:3];

   // Intermediate results for M1-M7
   wire signed [2*w+3:0] M1r[0:3], M1i[0:3];    //  2w +4 bits
   wire signed [2*w+3:0] M2r[0:3], M2i[0:3];
   wire signed [2*w+3:0] M3r[0:3], M3i[0:3];
   wire signed [2*w+3:0] M4r[0:3], M4i[0:3];
   wire signed [2*w+3:0] M5r[0:3], M5i[0:3];
	wire signed [2*w+3:0] M6r[0:3], M6i[0:3];
	wire signed [2*w+3:0] M7r[0:3], M7i[0:3];

   // M1 = (A11 + A22) * (B11 + B22)
   matrix_add_2x2_complex #(w)		// results on w+1 bits
		add1 (.Ar(A11r), .Br(A22r), .Cr(S1r),  
				.Ai(A11i), .Bi(A22i), .Ci(S1i));  
	matrix_add_2x2_complex #(w)		// results on w+1 bits
		add2 (.Ar(B11r), .Br(B22r), .Cr(T1r),
				.Ai(B11i), .Bi(B22i), .Ci(T1i));  
   strassen_2x2_complex #(w+1) mul1 //   2w+3
		(.Ar(S1r), .Br(T1r), .Cr(M1r),
		 .Ai(S1i), .Bi(T1i), .Ci(M1i) );  
		
    // M2 = (A21 + A22) * B11
    matrix_add_2x2_complex #(w)		// results on w+1 bits
		add3 (.Ar(A21r), .Br(A22r), .Cr(S2r),
				.Ai(A21i), .Bi(A22i), .Ci(S2i));  
    assign B11r_ext[0] = $signed(B11r[0]);
	 assign B11r_ext[1] = $signed(B11r[1]);
	 assign B11r_ext[2] = $signed(B11r[2]);
	 assign B11r_ext[3] = $signed(B11r[3]);																			 
    assign B11i_ext[0] = $signed(B11i[0]);
	 assign B11i_ext[1] = $signed(B11i[1]);
	 assign B11i_ext[2] = $signed(B11i[2]);
	 assign B11i_ext[3] = $signed(B11i[3]);																			 
    strassen_2x2_complex #(w+1) mul2 //   2w+3
		(.Ar(S2r), .Br(B11r_ext), .Cr(M2r),
		 .Ai(S2i), .Bi(B11i_ext), .Ci(M2i) );  
			 
    // M3 = A11 * (B12 - B22)
    matrix_sub_2x2_complex #(w)		// results on w+1 bits
		sub1 (.Ar(B12r), .Br(B22r), .Cr(T2r),
				.Ai(B12i), .Bi(B22i), .Ci(T2i));  
	 assign A11r_ext[0] = $signed(A11r[0]);
	 assign A11r_ext[1] = $signed(A11r[1]);
	 assign A11r_ext[2] = $signed(A11r[2]);
	 assign A11r_ext[3] = $signed(A11r[3]);																			 
	 assign A11i_ext[0] = $signed(A11i[0]);
	 assign A11i_ext[1] = $signed(A11i[1]);
	 assign A11i_ext[2] = $signed(A11i[2]);
	 assign A11i_ext[3] = $signed(A11i[3]);																			 
    strassen_2x2_complex #(w+1) mul3 //   2w+3
		(.Ar(A11r_ext), .Br(T2r), .Cr(M3r),
		 .Ai(A11i_ext), .Bi(T2i), .Ci(M3i) );  

    // M4 = A22 * (B21 - B11)
    matrix_sub_2x2_complex #(w)		// results on w+1 bits
		sub2 (.Ar(B21r), .Br(B11r), .Cr(T3r),
				.Ai(B21i), .Bi(B11i), .Ci(T3i));  
	 assign A22r_ext[0] = $signed(A22r[0]);
	 assign A22r_ext[1] = $signed(A22r[1]);
	 assign A22r_ext[2] = $signed(A22r[2]);
	 assign A22r_ext[3] = $signed(A22r[3]);																			 
	 assign A22i_ext[0] = $signed(A22i[0]);
	 assign A22i_ext[1] = $signed(A22i[1]);
	 assign A22i_ext[2] = $signed(A22i[2]);
	 assign A22i_ext[3] = $signed(A22i[3]);																			 
    strassen_2x2_complex #(w+1) mul4 //   2w+3
		(.Ar(A22r_ext), .Br(T3r), .Cr(M4r),
		 .Ai(A22i_ext), .Bi(T3i), .Ci(M4i) );  

    // M5 = (A11 + A12) * B22
    matrix_add_2x2_complex #(w)		// results on w+1 bits
		add4 (.Ar(A11r), .Br(A12r), .Cr(S3r),
				.Ai(A11i), .Bi(A12i), .Ci(S3i));  
	 assign B22r_ext[0] = $signed(B22r[0]);
	 assign B22r_ext[1] = $signed(B22r[1]);
	 assign B22r_ext[2] = $signed(B22r[2]);
	 assign B22r_ext[3] = $signed(B22r[3]);																			 
	 assign B22i_ext[0] = $signed(B22i[0]);
	 assign B22i_ext[1] = $signed(B22i[1]);
	 assign B22i_ext[2] = $signed(B22i[2]);
	 assign B22i_ext[3] = $signed(B22i[3]);
    strassen_2x2_complex #(w+1) mul5 //   2w+3
		(.Ar(S3r), .Br(B22r_ext), .Cr(M5r),
		 .Ai(S3i), .Bi(B22i_ext), .Ci(M5i) );  

    // M6 = (A21 - A11) * (B11 + B12)
    matrix_sub_2x2_complex #(w)		// results on w+1 bits
		sub3 (.Ar(A21r), .Br(A11r), .Cr(S4r),
				.Ai(A21i), .Bi(A11i), .Ci(S4i));  
	 
    matrix_add_2x2_complex #(w)		// results on w+1 bits
		add5 (.Ar(B11r), .Br(B12r), .Cr(T4r),
				.Ai(B11i), .Bi(B12i), .Ci(T4i));  
    strassen_2x2_complex #(w+1) mul6 //   2w+3
		(.Ar(S4r), .Br(T4r), .Cr(M6r),
		 .Ai(S4i), .Bi(T4i), .Ci(M6i) );  

    // M7 = (A12 - A22) * (B21 + B22)
    matrix_sub_2x2_complex #(w)		// results on w+1 bits
		sub4 (.Ar(A12r), .Br(A22r), .Cr(S5r),
				.Ai(A12i), .Bi(A22i), .Ci(S5i));  
    matrix_add_2x2_complex #(w)		// results on w+1 bits
		add6 (.Ar(B21r), .Br(B22r), .Cr(T5r),
				.Ai(B21i), .Bi(B22i), .Ci(T5i));  
    strassen_2x2_complex #(w+1) mul7 //   2w+3
		(.Ar(S5r), .Br(T5r), .Cr(M7r),
		 .Ai(S5i), .Bi(T5i), .Ci(M7i) );  


	 // extended C values
    wire signed [2*w+6:0] C11r_ext[0:3], C12r_ext[0:3],  C21r_ext[0:3],  C22r_ext[0:3];  // 2w+6 bits
    wire signed [2*w+6:0] C11i_ext[0:3], C12i_ext[0:3],  C21i_ext[0:3],  C22i_ext[0:3];  // 2w+6 bits
	 // Compute final submatrices
    wire signed [2*w+2:0] C11r[0:3], C12r[0:3], C21r[0:3], C22r[0:3];  // 2w+3 bits
    wire signed [2*w+2:0] C11i[0:3], C12i[0:3], C21i[0:3], C22i[0:3];  // 2w+3 bits

    // C11 = M1 + M4 - M5 + M7
	 assign C11r_ext[0] = M1r[0] + M4r[0] - M5r[0] + M7r[0];  
	 assign C11r_ext[1] = M1r[1] + M4r[1] - M5r[1] + M7r[1];
	 assign C11r_ext[2] = M1r[2] + M4r[2] - M5r[2] + M7r[2];
	 assign C11r_ext[3] = M1r[3] + M4r[3] - M5r[3] + M7r[3];
	 assign C11i_ext[0] = M1i[0] + M4i[0] - M5i[0] + M7i[0];  
	 assign C11i_ext[1] = M1i[1] + M4i[1] - M5i[1] + M7i[1];
	 assign C11i_ext[2] = M1i[2] + M4i[2] - M5i[2] + M7i[2];
	 assign C11i_ext[3] = M1i[3] + M4i[3] - M5i[3] + M7i[3];
 	 assign C11r[0] = C11r_ext[0][2*w+2:0];   
 	 assign C11r[1] = C11r_ext[1][2*w+2:0];   
 	 assign C11r[2] = C11r_ext[2][2*w+2:0];   
 	 assign C11r[3] = C11r_ext[3][2*w+2:0];   
 	 assign C11i[0] = C11i_ext[0][2*w+2:0];   
 	 assign C11i[1] = C11i_ext[1][2*w+2:0];   
 	 assign C11i[2] = C11i_ext[2][2*w+2:0];   
 	 assign C11i[3] = C11i_ext[3][2*w+2:0];   

	 // C12 = M3 + M5
	 assign C12r_ext[0] = M3r[0] + M5r[0]; 
	 assign C12r_ext[1] = M3r[1] + M5r[1];
	 assign C12r_ext[2] = M3r[2] + M5r[2];
	 assign C12r_ext[3] = M3r[3] + M5r[3];
	 assign C12i_ext[0] = M3i[0] + M5i[0]; 
	 assign C12i_ext[1] = M3i[1] + M5i[1];
	 assign C12i_ext[2] = M3i[2] + M5i[2];
	 assign C12i_ext[3] = M3i[3] + M5i[3];
 	 assign C12r[0] = C12r_ext[0][2*w+2:0];   
 	 assign C12r[1] = C12r_ext[1][2*w+2:0];   
 	 assign C12r[2] = C12r_ext[2][2*w+2:0];   
 	 assign C12r[3] = C12r_ext[3][2*w+2:0];   
 	 assign C12i[0] = C12i_ext[0][2*w+2:0];   
 	 assign C12i[1] = C12i_ext[1][2*w+2:0];   
 	 assign C12i[2] = C12i_ext[2][2*w+2:0];   
 	 assign C12i[3] = C12i_ext[3][2*w+2:0];   

    // C21 = M2 + M4
	 assign C21r_ext[0] = M2r[0] + M4r[0]; 
	 assign C21r_ext[1] = M2r[1] + M4r[1];
	 assign C21r_ext[2] = M2r[2] + M4r[2];
	 assign C21r_ext[3] = M2r[3] + M4r[3];
	 assign C21i_ext[0] = M2i[0] + M4i[0]; 
	 assign C21i_ext[1] = M2i[1] + M4i[1];
	 assign C21i_ext[2] = M2i[2] + M4i[2];
	 assign C21i_ext[3] = M2i[3] + M4i[3];
	 assign C21r[0] = C21r_ext[0][2*w+2:0]; 
	 assign C21r[1] = C21r_ext[1][2*w+2:0]; 
	 assign C21r[2] = C21r_ext[2][2*w+2:0]; 
	 assign C21r[3] = C21r_ext[3][2*w+2:0]; 
	 assign C21i[0] = C21i_ext[0][2*w+2:0]; 
	 assign C21i[1] = C21i_ext[1][2*w+2:0]; 
	 assign C21i[2] = C21i_ext[2][2*w+2:0]; 
	 assign C21i[3] = C21i_ext[3][2*w+2:0]; 

    // C22 = M1 - M2 + M3 + M6
	 assign C22r_ext[0] = M1r[0] - M2r[0] + M3r[0] + M6r[0];
	 assign C22r_ext[1] = M1r[1] - M2r[1] + M3r[1] + M6r[1];
	 assign C22r_ext[2] = M1r[2] - M2r[2] + M3r[2] + M6r[2];
	 assign C22r_ext[3] = M1r[3] - M2r[3] + M3r[3] + M6r[3];
	 assign C22i_ext[0] = M1i[0] - M2i[0] + M3i[0] + M6i[0];
	 assign C22i_ext[1] = M1i[1] - M2i[1] + M3i[1] + M6i[1];
	 assign C22i_ext[2] = M1i[2] - M2i[2] + M3i[2] + M6i[2];
	 assign C22i_ext[3] = M1i[3] - M2i[3] + M3i[3] + M6i[3];
 	 assign C22r[0] = C22r_ext[0][2*w+2:0];   
 	 assign C22r[1] = C22r_ext[1][2*w+2:0];   
 	 assign C22r[2] = C22r_ext[2][2*w+2:0];   
 	 assign C22r[3] = C22r_ext[3][2*w+2:0];   
 	 assign C22i[0] = C22i_ext[0][2*w+2:0];   
 	 assign C22i[1] = C22i_ext[1][2*w+2:0];   
 	 assign C22i[2] = C22i_ext[2][2*w+2:0];   
 	 assign C22i[3] = C22i_ext[3][2*w+2:0];   

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



module matrix_add_2x2_complex #(
    parameter w = 8
)(
    input  wire signed [w-1:0] Ar [0:3],
    input  wire signed [w-1:0] Ai [0:3],
    input  wire signed [w-1:0] Br [0:3],
    input  wire signed [w-1:0] Bi [0:3],
    output wire signed [w:0]   Cr [0:3],  // Result width is w+1 to accommodate overflow
    output wire signed [w:0]   Ci [0:3]   // Result width is w+1 to accommodate overflow
);
    assign Cr[0] = Ar[0] + Br[0];
    assign Cr[1] = Ar[1] + Br[1];
    assign Cr[2] = Ar[2] + Br[2];
    assign Cr[3] = Ar[3] + Br[3];
    assign Ci[0] = Ai[0] + Bi[0];
    assign Ci[1] = Ai[1] + Bi[1];
    assign Ci[2] = Ai[2] + Bi[2];
    assign Ci[3] = Ai[3] + Bi[3];
endmodule


module matrix_sub_2x2_complex #(
    parameter w = 8
)(
    input  wire signed [w-1:0] Ar [0:3],
    input  wire signed [w-1:0] Ai [0:3],
    input  wire signed [w-1:0] Br [0:3],
    input  wire signed [w-1:0] Bi [0:3],
    output wire signed [w:0]   Cr [0:3],  // Result width is w+1 to accommodate overflow
    output wire signed [w:0]   Ci [0:3]  // Result width is w+1 to accommodate overflow
);
    assign Cr[0] = Ar[0] - Br[0];
    assign Cr[1] = Ar[1] - Br[1];
    assign Cr[2] = Ar[2] - Br[2];
    assign Cr[3] = Ar[3] - Br[3];
    assign Ci[0] = Ai[0] - Bi[0];
    assign Ci[1] = Ai[1] - Bi[1];
    assign Ci[2] = Ai[2] - Bi[2];
    assign Ci[3] = Ai[3] - Bi[3];
endmodule

