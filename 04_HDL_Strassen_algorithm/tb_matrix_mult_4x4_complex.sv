// SPDX-License-Identifier: CERN-OHL-NC-2.0
// Author: Ettore Napoli
// Affiliation: University of Salerno
// January 2026
// Description: Test bench - Strassen algorithm for matrix multiplication

`timescale 1ns / 1ps

module tb_matrix_mult_4x4_complex();

    parameter w = 48;
    parameter WIDTH_OUT = 2*w + 3;

    logic signed [w-1:0] A_real [0:3][0:3];  // rows, columns. Matlab notation
    logic signed [w-1:0] B_real [0:3][0:3];  // rows, columns. Matlab notation
    logic signed [w-1:0] A_imag [0:3][0:3];  // rows, columns. Matlab notation
    logic signed [w-1:0] B_imag [0:3][0:3];  // rows, columns. Matlab notation
    logic signed [WIDTH_OUT-1:0] C_real [0:3][0:3];
    logic signed [WIDTH_OUT-1:0] C_imag [0:3][0:3];

	 // 550 bit values read from file
    logic signed [549:0] A550r [0:3][0:3], B550r [0:3][0:3],golden_C550r [0:3][0:3];
    logic signed [549:0] A550i [0:3][0:3], B550i [0:3][0:3],golden_C550i [0:3][0:3];
    logic signed [WIDTH_OUT-1:0] golden_Cr [0:3][0:3], golden_Ci [0:3][0:3];
    logic signed [549:0] tmp550;
	 
	 
// DUT instantiation
matrix_mult_4x4_complex_strassen #(.w(w)) dut (
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
	file = $fopen("../../golden_values_32bit_ver_3300_test_vectors.txt", "r");
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
		test = test +1;
		$display("Test n.: %0d", test);

		// extract the 550 bit and 'w' bit values for A,B and golden_C
		for (i = 0; i < 4; i = i + 1)
			for (j = 0; j < 4; j = j + 1)
			begin
				substrAr = line.substr(209+(i*6080)+(j*1520), 209+(i*6080)+(j*1520)+549);
				substrAi = line.substr(969+(i*6080)+(j*1520), 969+(i*6080)+(j*1520)+549);
				substrBr = line.substr(24529+(i*6080)+(j*1520), 24529+(i*6080)+(j*1520)+549);
				substrBi = line.substr(25289+(i*6080)+(j*1520), 25289+(i*6080)+(j*1520)+549);
				substrCr = line.substr(48849+(i*6080)+(j*1520), 48849+(i*6080)+(j*1520)+549);
				substrCi = line.substr(49609+(i*6080)+(j*1520), 49609+(i*6080)+(j*1520)+549);
//$display("test Ar: %s",substrAr);
//$display("test Ai: %s",substrAi);
//$display("test Br: %s",substrBr);
//$display("test Bi: %s",substrBi);
//$display("test Cr: %s",substrCr);
//$display("test Ci: %s",substrCi);
				for (k = 0; k < 550; k++)
				begin
					if (substrAr[k] == "1") 
						tmp550[549 - k] = 1'b1;   // Store MSB at left (big endian)
					else 
						if (substrAr[k] == "0")	tmp550[549 - k] = 1'b0;
						else	tmp550[549 - k] = 1'bx;
				end
//$display("test tmp550: %b",tmp550);
				
				A550r[i][j]=tmp550;	
				A_real[i][j]=A550r[i][j][w-1:0];
/*				if ((i==0)&(j==0)&(test==2510))
					begin
					$display("substrAr=%s",substrAr);
					$display("A550r[0][0]=%b",A550r[i][j]);
					$display("A_real[0][0]=%b",A_real[i][j]);
					end */
				for (k = 0; k < 550; k++)
				begin
					if (substrAi[k] == "1") 
						tmp550[549 - k] = 1'b1;   // Store MSB at left (big endian)
					else 
						if (substrAi[k] == "0")	tmp550[549 - k] = 1'b0;
						else	tmp550[549 - k] = 1'bx;
				end
				A550i[i][j]=tmp550;	
				A_imag[i][j]=A550i[i][j][w-1:0];	
							
				for (k = 0; k < 550; k++)
				begin
					if (substrBr[k] == "1") 
						tmp550[549 - k] = 1'b1;   // Store MSB at left (big endian)
					else 
						if (substrBr[k] == "0")	tmp550[549 - k] = 1'b0;
						else	tmp550[549 - k] = 1'bx;
				end
				B550r[i][j]=tmp550;	
				B_real[i][j]=B550r[i][j][w-1:0];	
				
				for (k = 0; k < 550; k++)
				begin
					if (substrBi[k] == "1") 
						tmp550[549 - k] = 1'b1;   // Store MSB at left (big endian)
					else 
						if (substrBi[k] == "0")	tmp550[549 - k] = 1'b0;
						else	tmp550[549 - k] = 1'bx;
				end
				B550i[i][j]=tmp550;	
				B_imag[i][j]=B550i[i][j][w-1:0];	
				

				for (k = 0; k < 550; k++)
				begin
					if (substrCr[k] == "1") 
						tmp550[549 - k] = 1'b1;   // Store MSB at left (big endian)
					else 
						if (substrCr[k] == "0")	tmp550[549 - k] = 1'b0;
						else	tmp550[549 - k] = 1'bx;
				end
				golden_C550r[i][j]=tmp550;	
				golden_Cr[i][j]=golden_C550r[i][j][WIDTH_OUT-1:0];	

				for (k = 0; k < 550; k++)
				begin
					if (substrCi[k] == "1") 
						tmp550[549 - k] = 1'b1;   // Store MSB at left (big endian)
					else 
						if (substrCi[k] == "0")	tmp550[549 - k] = 1'b0;
						else	tmp550[549 - k] = 1'bx;
				end
				golden_C550i[i][j]=tmp550;	
				golden_Ci[i][j]=golden_C550i[i][j][WIDTH_OUT-1:0];	
//$display("test Areal: %b",A_real[i][j]);
//$display("test Aimag: %b",A_imag[i][j]);
//$display("test Brreal: %b",B_real[i][j]);
//$display("test Bimag: %b",B_imag[i][j]);
//$display("test Crreal: %b",golden_Cr[i][j]);
//$display("test Cimag: %b",golden_Ci[i][j]);
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
		/*=========================================
		if ((test>4000)|(test==4000))
		begin
			$display("line: %s",line);
			for (int i = 0; i < 4; i++)
					begin
						$display("A[%0d][1:3]: %0d + %0di, %0d + %0di, %0d + %0di, %0d + %0di", i,
						A_real[i][0],A_imag[i][0],A_real[i][1],A_imag[i][1],
						A_real[i][2],A_imag[i][2],A_real[i][3],A_imag[i][3]);
					end
			for (int i = 0; i < 4; i++)
						$display("B[%0d][1:3]: %0d + %0di, %0d + %0di, %0d + %0di, %0d + %0di", i,
						B_real[i][0],B_imag[i][0],B_real[i][1],B_imag[i][1],
						B_real[i][2],B_imag[i][2],B_real[i][3],B_imag[i][3]);
			for (int i = 0; i < 4; i++)
					begin
						$display("C[%0d][1:3]: %0d + %0di, %0d + %0di, %0d + %0di, %0d + %0di", i,
						C_real[i][0],C_imag[i][0],C_real[i][1],C_imag[i][1],
						C_real[i][2],C_imag[i][2],C_real[i][3],C_imag[i][3]);
					end					
			for (int i = 0; i < 4; i++)
					begin
						$display("golden_C[%0d][1:3]: %0d + %0di, %0d + %0di, %0d + %0di, %0d + %0di", i,
						golden_Cr[i][0],golden_Ci[i][0],golden_Cr[i][1],golden_Ci[i][1],
						golden_Cr[i][2],golden_Ci[i][2],golden_Cr[i][3],golden_Ci[i][3]);
					end
		end  */
//		if (test==2510) $stop;		
				
   end
	
	$fclose(file);
	if (error_count == 0)
		$display("All tests passed!");
	else
		$display("%0d mismatches found.", error_count);

	$stop;
end

endmodule
