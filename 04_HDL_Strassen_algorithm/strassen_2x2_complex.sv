// SPDX-License-Identifier: CERN-OHL-NC-2.0
// Author: Ettore Napoli
// Affiliation: University of Salerno
// january 2026
// Description: Strassen algorithm for 2x2 complex matrix multiplication

module strassen_2x2_complex #(
    parameter W = 8
)(
    input  wire signed [W-1:0] Ar[0:3],
    input  wire signed [W-1:0] Br[0:3],
    input  wire signed [W-1:0] Ai[0:3],
    input  wire signed [W-1:0] Bi[0:3],
    output wire signed [2*W+1:0] Cr[0:3], // allow for full overflow range 2*W+1 bits. Use 2*W+2 bits for coherency with imag part
    output wire signed [2*W+1:0] Ci[0:3]  // allow for full overflow range. Output on 2*W+2 bits
);
    // Local declarations for matrix elements
	 // REAL part
    wire signed [W-1:0] ar11 = Ar[0];
    wire signed [W-1:0] ar12 = Ar[1];
    wire signed [W-1:0] ar21 = Ar[2];
    wire signed [W-1:0] ar22 = Ar[3];

    wire signed [W-1:0] br11 = Br[0];
    wire signed [W-1:0] br12 = Br[1];
    wire signed [W-1:0] br21 = Br[2];
    wire signed [W-1:0] br22 = Br[3];
	 // IMAG part
    wire signed [W-1:0] ai11 = Ai[0];
    wire signed [W-1:0] ai12 = Ai[1];
    wire signed [W-1:0] ai21 = Ai[2];
    wire signed [W-1:0] ai22 = Ai[3];

    wire signed [W-1:0] bi11 = Bi[0];
    wire signed [W-1:0] bi12 = Bi[1];
    wire signed [W-1:0] bi21 = Bi[2];
    wire signed [W-1:0] bi22 = Bi[3];
	 
	 // intermediate sums (W + 1 bits)
	 wire signed  [W:0] ar11par22, ai11pai22;  // M1
	 wire signed  [W:0] br11pbr22, bi11pbi22;  // M1  
	 wire signed  [W:0] ar21par22, ai21pai22;  // M2    
	 wire signed  [W:0] br12mbr22, bi12mbi22;  // M3  
	 wire signed  [W:0] br21mbr11, bi21mbi11;  // M4  
 	 wire signed  [W:0] ar11par12, ai11pai12;  // M5
 	 wire signed  [W:0] ar21mar11, ai21mai11;  // M6  
	 wire signed  [W:0] br11pbr12, bi11pbi12;  // M6  
	 wire signed  [W:0] ar12mar22, ai12mai22;	 // M7
	 wire signed  [W:0] br21pbr22, bi21pbi22;	 // M7
 	 
	 
    // Intermediate products (2*W + 3 bits)
    wire signed [2*W+2:0] M1r; // w=3 inputs [-4,3] output max: (-4-4)*(-4-4) - (3+3)*(-4-4) = 64 + 48 = 112    min (-4-4)*(3+3) - (-4-4)*(-4-4)= -48-64 = -112 =>  8 bit  
    wire signed [2*W+2:0] M1i; // w=3 inputs [-4,3] output max: (-4-4)*(-4-4) + (-4-4)*(-4-4)= 64 + 64 = 128    min (-4-4)*(3+3) + (-4-4)*(3+3)= -48-48 = -96 =>  9 bit
    wire signed [2*W+2:0] M2r, M2i, M3r, M3i, M4r, M4i, M5r, M5i, M6r, M6i, M7r, M7i;

	 //  Cr and Ci extended 
	 wire signed [2*W+4:0] Cr_ext[0:3];  // 2*W+5 bits
	 wire signed [2*W+4:0] Ci_ext[0:3];  // 2*W+5 bits
	 
    // assign M1 = (a11 + a22) * (b11 + b22);
	 // intermediate sums
	 assign ar11par22 = ar11+ar22;  // W+1 bits
	 assign ai11pai22 = ai11+ai22;  // W+1 bits  
	 assign br11pbr22 = br11+br22;  // W+1 bits  
	 assign bi11pbi22 = bi11+bi22;  // W+1 bits  
	 mult_complex #(.W(W+1)) M1_mult (      //  on 2*(W+1)+1 = 2*W+3 bits 
		.Ar(ar11par22), .Ai(ai11pai22), .Br(br11pbr22), .Bi(bi11pbi22), .Cr(M1r), .Ci(M1i));

	 //	assign M2 = (a21 + a22) * b11;
	 // intermediate sums
	 assign ar21par22 = ar21+ar22;    // W+1 bits  
	 assign ai21pai22 = ai21+ai22;    // W+1 bits
	 mult_complex #(.W(W+1)) M2_mult (         //  on 2*(W+1)+1 = 2*W+3 bits 
		.Ar(ar21par22), .Ai(ai21pai22), .Br({br11[W-1],br11}), .Bi({bi11[W-1],bi11}), .Cr(M2r), .Ci(M2i));

	 // assign M3 = a11 * (b12 - b22);
	 // intermediate sums
	 assign br12mbr22 = br12-br22;  // W+1 bits  
	 assign bi12mbi22 = bi12-bi22;  // W+1 bits  
	 mult_complex #(.W(W+1)) M3_mult (         //  on 2*(W+1)+1 = 2*W+3 bits 
		.Ar({ar11[W-1],ar11}), .Ai({ai11[W-1],ai11}), .Br(br12mbr22), .Bi(bi12mbi22), .Cr(M3r), .Ci(M3i));
   
	 // assign M4 = a22 * (b21 - b11);
    // intermediate sums
	 assign br21mbr11 = br21-br11;  // W+1 bits  
	 assign bi21mbi11 = bi21-bi11;  // W+1 bits  
	 mult_complex #(.W(W+1)) M4_mult (         //  on 2*(W+1)+1 = 2*W+3 bits 
		.Ar({ar22[W-1],ar22}), .Ai({ai22[W-1],ai22}), .Br(br21mbr11), .Bi(bi21mbi11), .Cr(M4r), .Ci(M4i));

    // assign M5 = (a11 + a12) * b22;
    // intermediate sums
	 assign ar11par12 = ar11+ar12;  // W+1 bits  
	 assign ai11pai12 = ai11+ai12;  // W+1 bits  
	 mult_complex #(.W(W+1)) M5_mult (         //  on 2*(W+1)+1 = 2*W+3 bits 
		.Ar(ar11par12), .Ai(ai11pai12), .Br({br22[W-1],br22}), .Bi({bi22[W-1],bi22}), .Cr(M5r), .Ci(M5i));

    // assign M6 = (a21 - a11) * (b11 + b12);
    // intermediate sums
	 assign ar21mar11 = ar21-ar11;  // W+1 bits  
	 assign ai21mai11 = ai21-ai11;  // W+1 bits  
	 assign br11pbr12 = br11+br12;  // W+1 bits  
	 assign bi11pbi12 = bi11+bi12;  // W+1 bits  
	 mult_complex #(.W(W+1)) M6_mult (         //  on 2*(W+1)+1 = 2*W+3 bits 
		.Ar(ar21mar11), .Ai(ai21mai11), .Br(br11pbr12), .Bi(bi11pbi12), .Cr(M6r), .Ci(M6i));

    // assign M7 = (a12 - a22) * (b21 + b22);
    // intermediate sums
	 assign ar12mar22 = ar12-ar22;  // W+1 bits  
	 assign ai12mai22 = ai12-ai22;  // W+1 bits  
	 assign br21pbr22 = br21+br22;  // W+1 bits  
	 assign bi21pbi22 = bi21+bi22;  // W+1 bits  
	 mult_complex  #(.W(W+1)) M7_mult (         //  on 2*(W+1)+1 = 2*W+3 bits 
		.Ar(ar12mar22), .Ai(ai12mai22), .Br(br21pbr22), .Bi(bi21pbi22), .Cr(M7r), .Ci(M7i));

    // Final result computation (2W+2 bits for the real and imag parts)
	 // Cr_ext is extended allowing the maximum result of the sum
	 // It is then truncated to get the final result
    assign Cr_ext[0] = M1r + M4r - M5r + M7r;  
    assign Cr_ext[1] = M3r + M5r;
    assign Cr_ext[2] = M2r + M4r;
    assign Cr_ext[3] = M1r - M2r + M3r + M6r;
    assign Ci_ext[0] = M1i + M4i - M5i + M7i;  
    assign Ci_ext[1] = M3i + M5i;
    assign Ci_ext[2] = M2i + M4i;
    assign Ci_ext[3] = M1i - M2i + M3i + M6i;


    assign Cr[0] = Cr_ext[0][2*W+1:0];  
    assign Cr[1] = Cr_ext[1][2*W+1:0];  
    assign Cr[2] = Cr_ext[2][2*W+1:0];  
    assign Cr[3] = Cr_ext[3][2*W+1:0];  
    assign Ci[0] = Ci_ext[0][2*W+1:0];  
    assign Ci[1] = Ci_ext[1][2*W+1:0];  
    assign Ci[2] = Ci_ext[2][2*W+1:0];  
    assign Ci[3] = Ci_ext[3][2*W+1:0];  

endmodule
