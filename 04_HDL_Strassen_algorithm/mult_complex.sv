// SPDX-License-Identifier: CERN-OHL-NC-2.0
// Author: Ettore Napoli
// Affiliation: University of Salerno
// january 2026
// Description: complex multiplier

module mult_complex #(
    parameter W = 16   // bit-width of real and imaginary parts
)(
    input  wire signed [W-1:0] Ar,
    input  wire signed [W-1:0] Ai,
    input  wire signed [W-1:0] Br,
    input  wire signed [W-1:0] Bi,

    output wire signed [2*W:0] Cr,  // 2*W + 1 bits 
    output wire signed [2*W:0] Ci   // 2*W + 1 bits 
);

    // Intermediate products (2W bits)
    wire signed [2*W-1:0] ArBr = Ar * Br;
    wire signed [2*W-1:0] AiBi = Ai * Bi;
    wire signed [2*W-1:0] ArBi = Ar * Bi;
    wire signed [2*W-1:0] AiBr = Ai * Br;

    // Complex multiplication result
    assign Cr = ArBr - AiBi;
    assign Ci = ArBi + AiBr;

endmodule
