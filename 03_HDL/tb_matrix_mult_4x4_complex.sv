// SPDX-License-Identifier: CERN-OHL-NC-2.0
// Author: Ettore Napoli
// Affiliation: University of Salerno
// August 2025
// Description: AlphaEvolve algorithm for matrix multiplication

`timescale 1ns / 1ps

module tb_matrix_mult_4x4_complex();

    parameter w = 16;
    parameter WIDTH_OUT = 2*w + 3;

    logic signed [w-1:0] A_real [0:3][0:3];  // rows, columns. Matlab notation
    logic signed [w-1:0] B_real [0:3][0:3];  // rows, columns. Matlab notation
    logic signed [w-1:0] A_imag [0:3][0:3];  // rows, columns. Matlab notation
    logic signed [w-1:0] B_imag [0:3][0:3];  // rows, columns. Matlab notation
    logic signed [WIDTH_OUT-1:0] C_real [0:3][0:3];
    logic signed [WIDTH_OUT-1:0] C_imag [0:3][0:3];

	 // 128 bit values read from file
    logic signed [127:0] A128r [0:3][0:3], B128r [0:3][0:3],golden_C128r [0:3][0:3];
    logic signed [127:0] A128i [0:3][0:3], B128i [0:3][0:3],golden_C128i [0:3][0:3];
    logic signed [WIDTH_OUT-1:0] golden_Cr [0:3][0:3], golden_Ci [0:3][0:3];
    logic signed [127:0] tmp128;
	 
	 
// DUT instantiation
matrix_mult_4x4_complex_alphaevolve #(.w(w)) dut (
    .A_real(A_real),
    .A_imag(A_imag),
    .B_real(B_real),
    .B_imag(B_imag),
    .C_real(C_real),
    .C_imag(C_imag)
);

 

int file;
string line;
string substrAr,substrBr,substrCr;
string substrAi,substrBi,substrCi;
int i,j,k;
int error_count;
int test;

initial begin
	$dumpfile("dump.vcd");
	file = $fopen("../../golden_values_16bit_ver_3300_test_vectors.txt", "r");
	if (file == 0)
	begin
		$display("ERROR: Cannot open file.");
		$finish;
	end

	error_count = 0;
	test = 0;
	// Loop until EOF
	while ($fgets(line, file))
	begin
		test ++;
		$display("Test n.: %0d", test);

		// extract the 128 bit and 'w' bit values for A,B and golden_C
		for (i = 0; i < 4; i = i + 1)
			for (j = 0; j < 4; j = j + 1)
			begin
				substrAr = line.substr(39+(i*1344)+(j*336), 39+(i*1344)+(j*336)+127);
				substrAi = line.substr(207+(i*1344)+(j*336), 207+(i*1344)+(j*336)+127);
				substrBr = line.substr(5415+(i*1344)+(j*336), 5415+(i*1344)+(j*336)+127);
				substrBi = line.substr(5583+(i*1344)+(j*336), 5583+(i*1344)+(j*336)+127);
				substrCr = line.substr(10791+(i*1344)+(j*336), 10791+(i*1344)+(j*336)+127);
				substrCi = line.substr(10959+(i*1344)+(j*336), 10959+(i*1344)+(j*336)+127);

				for (k = 0; k < 128; k++)
				begin
					if (substrAr[k] == "1") 
						tmp128[127 - k] = 1'b1;   // Store MSB at left (big endian)
					else 
						if (substrAr[k] == "0")	tmp128[127 - k] = 1'b0;
						else	tmp128[127 - k] = 1'bx;
				end
				A128r[i][j]=tmp128;	
				A_real[i][j]=A128r[i][j][w-1:0];	
				
				for (k = 0; k < 128; k++)
				begin
					if (substrAi[k] == "1") 
						tmp128[127 - k] = 1'b1;   // Store MSB at left (big endian)
					else 
						if (substrAi[k] == "0")	tmp128[127 - k] = 1'b0;
						else	tmp128[127 - k] = 1'bx;
				end
				A128i[i][j]=tmp128;	
				A_imag[i][j]=A128i[i][j][w-1:0];	
							
				for (k = 0; k < 128; k++)
				begin
					if (substrBr[k] == "1") 
						tmp128[127 - k] = 1'b1;   // Store MSB at left (big endian)
					else 
						if (substrBr[k] == "0")	tmp128[127 - k] = 1'b0;
						else	tmp128[127 - k] = 1'bx;
				end
				B128r[i][j]=tmp128;	
				B_real[i][j]=B128r[i][j][w-1:0];	
				
				for (k = 0; k < 128; k++)
				begin
					if (substrBi[k] == "1") 
						tmp128[127 - k] = 1'b1;   // Store MSB at left (big endian)
					else 
						if (substrBi[k] == "0")	tmp128[127 - k] = 1'b0;
						else	tmp128[127 - k] = 1'bx;
				end
				B128i[i][j]=tmp128;	
				B_imag[i][j]=B128i[i][j][w-1:0];	
				

				for (k = 0; k < 128; k++)
				begin
					if (substrCr[k] == "1") 
						tmp128[127 - k] = 1'b1;   // Store MSB at left (big endian)
					else 
						if (substrCr[k] == "0")	tmp128[127 - k] = 1'b0;
						else	tmp128[127 - k] = 1'bx;
				end
				golden_C128r[i][j]=tmp128;	
				golden_Cr[i][j]=golden_C128r[i][j][WIDTH_OUT-1:0];	

				for (k = 0; k < 128; k++)
				begin
					if (substrCi[k] == "1") 
						tmp128[127 - k] = 1'b1;   // Store MSB at left (big endian)
					else 
						if (substrCi[k] == "0")	tmp128[127 - k] = 1'b0;
						else	tmp128[127 - k] = 1'bx;
				end
				golden_C128i[i][j]=tmp128;	
				golden_Ci[i][j]=golden_C128i[i][j][WIDTH_OUT-1:0];	
			end	
		//  simulate
		#10;
		// Compare DUT output with golden_C
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
				if (C_real[i][j] !== golden_Cr[i][j])
				begin
					$display("Mismatch in test %0d at Cr[%0d][%0d]: DUT=%0d, GOLD=%0d", test, i, j, C_real[i][j], golden_Cr[i][j]);
					error_count++;
			   end
		
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
				if (C_imag[i][j] !== golden_Ci[i][j])
				begin
					$display("Mismatch in test %0d at Ci[%0d][%0d]: DUT=%0d, GOLD=%0d", test, i, j, C_imag[i][j], golden_Ci[i][j]);
					error_count++;
			   end
/*		//=========================================
		if (test==2303)
		begin
			$display("Dumping signals");
			$dumpvars(0,tb_matrix_mult_4x4_complex.dut);
		end
		if (test==2305)
		begin
			$dumpoff;
			$display("sto dumping signals");
		end
			$display("line: %s",line);
			for (int i = 0; i < 4; i++)
						$display("A[%0d][1:3]: %0d + %0di, %0d + %0di, %0d + %0di, %0d + %0di", i,
						A_real[i][0],A_imag[i][0],A_real[i][1],A_imag[i][1],
						A_real[i][2],A_imag[i][2],A_real[i][3],A_imag[i][3]);
			for (int i = 0; i < 4; i++)
						$display("B[%0d][1:3]: %0d + %0di, %0d + %0di, %0d + %0di, %0d + %0di", i,
						B_real[i][0],B_imag[i][0],B_real[i][1],B_imag[i][1],
						B_real[i][2],B_imag[i][2],B_real[i][3],B_imag[i][3]);
			for (int i = 0; i < 4; i++)
						$display("C[%0d][1:3]: %0d + %0di, %0d + %0di, %0d + %0di, %0d + %0di", i,
						C_real[i][0],C_imag[i][0],C_real[i][1],C_imag[i][1],
						C_real[i][2],C_imag[i][2],C_real[i][3],C_imag[i][3]);
			for (int i = 0; i < 4; i++)
						$display("golden_C[%0d][1:3]: %0d + %0di, %0d + %0di, %0d + %0di, %0d + %0di", i,
						golden_Cr[i][0],golden_Ci[i][0],golden_Cr[i][1],golden_Ci[i][1],
						golden_Cr[i][2],golden_Ci[i][2],golden_Cr[i][3],golden_Ci[i][3]);
		end  */
   end
	
	$fclose(file);
	if (error_count == 0)
		$display("All tests passed!");
	else
		$display("%0d mismatches found.", error_count);

	$stop;
end

endmodule
