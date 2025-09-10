// SPDX-License-Identifier: CERN-OHL-NC-2.0
// Author: Ettore Napoli
// Affiliation: University of Salerno
// August 2025
// Description: Naive algorithm for matrix multiplication

module matrix_mult_4x4_complex #(parameter w = 4)(
    input  signed [w-1:0] A_real [0:3][0:3],
    input  signed [w-1:0] A_imag [0:3][0:3],
    input  signed [w-1:0] B_real [0:3][0:3],
    input  signed [w-1:0] B_imag [0:3][0:3],
    output signed [2*w+2:0] C_real [0:3][0:3], //2*w +3 bits
    output signed [2*w+2:0] C_imag [0:3][0:3]  //2*w +3 bits
);

/*  Number of bits for the results.
Assume the input is on 3 bits [-4,3]
When two complex terms are multiplied (a+ib)+(c+id) = (ac-bd) +i(ad+bc)
The maximum of the real part is a=-4, c=-4, b=-4, d=+3 giving: 16+12=28.
When multiplying a row and a column the term can be summed four times getting a maximum value for the real part: 112.
The minimum of the real part is a=-4, c=+3, b=-4, d=-4 giving: -12-16=-28.
When multiplying a row and a column the term can be summed four times getting a minimum value for the real part: -112.
The number of bits to handle the range is 8 [-128,+127] that is 2*n+2


The maximum of the imaginary part is a=-4, d=-4, c=-4  b=-4  giving: 16+16=32.
When multiplying a row and a column the term can be summed four times getting a maximum value for the imaginary part: 128.
The minimum of the imaginary part is a=-4, d=+3, c=-4  b=+3  giving: -12-12=-24.
When multiplying a row and a column the term can be summed four times getting a minimum  value for the imaginary part: -96.

The number of bits to handle the range is 9 [-256,+255] that is 2*n+3   (the additional bit is just for the 128 case)

To keep coeherency, the output is kept to 2*n+3 bits for both real and imaginary parts.
*/

    integer i, j, k;
    reg signed [2*w+2:0] temp_real [0:3][0:3];
    reg signed [2*w+2:0] temp_imag [0:3][0:3];

    always @(*) begin
        for (i = 0; i < 4; i = i + 1) begin
            for (j = 0; j < 4; j = j + 1) begin
                temp_real[i][j] = 0;
                temp_imag[i][j] = 0;
                for (k = 0; k < 4; k = k + 1) begin
                    // Compute real part: (a*c - b*d)
                    temp_real[i][j] = temp_real[i][j] +
                        A_real[i][k] * B_real[k][j] - A_imag[i][k] * B_imag[k][j];
                    // Compute imag part: (a*d + b*c)
                    temp_imag[i][j] = temp_imag[i][j] +
                        A_real[i][k] * B_imag[k][j] + A_imag[i][k] * B_real[k][j];
                end
            end
        end
    end

    // Assign to outputs
    genvar x, y;
    generate
        for (x = 0; x < 4; x = x + 1) begin : row
            for (y = 0; y < 4; y = y + 1) begin : col
                assign C_real[x][y] = temp_real[x][y];
                assign C_imag[x][y] = temp_imag[x][y];
            end
        end
    endgenerate

endmodule
