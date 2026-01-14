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



/*==================================================================
The content of this file referes to a submitted scientific paper
The actual code will be released once the paper has been reviewed
=======================================================/*

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
