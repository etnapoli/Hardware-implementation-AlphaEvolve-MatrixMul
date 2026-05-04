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
wire signed [ 2 * w - 1 + 6 : 0] C10_i;  // signal num: 1
wire signed [ 2 * w - 1 + 5 : 0] C10_r;  // signal num: 2
wire signed [ 2 * w - 1 + 6 : 0] C11_i;  // signal num: 3
wire signed [ 2 * w - 1 + 5 : 0] C11_r;  // signal num: 4
wire signed [ 2 * w - 1 + 6 : 0] C12_i;  // signal num: 5
wire signed [ 2 * w - 1 + 5 : 0] C12_r;  // signal num: 6
wire signed [ 2 * w - 1 + 6 : 0] C13_i;  // signal num: 7
wire signed [ 2 * w - 1 + 5 : 0] C13_r;  // signal num: 8
wire signed [ 2 * w - 1 + 6 : 0] C14_i;  // signal num: 9
wire signed [ 2 * w - 1 + 5 : 0] C14_r;  // signal num: 10
wire signed [ 2 * w - 1 + 6 : 0] C15_i;  // signal num: 11
wire signed [ 2 * w - 1 + 5 : 0] C15_r;  // signal num: 12
wire signed [ 2 * w - 1 + 6 : 0] C16_i;  // signal num: 13
wire signed [ 2 * w - 1 + 5 : 0] C16_r;  // signal num: 14
wire signed [ 2 * w - 1 + 6 : 0] C1_i;  // signal num: 15
wire signed [ 2 * w - 1 + 5 : 0] C1_r;  // signal num: 16
wire signed [ 2 * w - 1 + 6 : 0] C2_i;  // signal num: 17
wire signed [ 2 * w - 1 + 5 : 0] C2_r;  // signal num: 18
wire signed [ 2 * w - 1 + 6 : 0] C3_i;  // signal num: 19
wire signed [ 2 * w - 1 + 5 : 0] C3_r;  // signal num: 20
wire signed [ 2 * w - 1 + 6 : 0] C4_i;  // signal num: 21
wire signed [ 2 * w - 1 + 5 : 0] C4_r;  // signal num: 22
wire signed [ 2 * w - 1 + 6 : 0] C5_i;  // signal num: 23
wire signed [ 2 * w - 1 + 5 : 0] C5_r;  // signal num: 24
wire signed [ 2 * w - 1 + 6 : 0] C6_i;  // signal num: 25
wire signed [ 2 * w - 1 + 5 : 0] C6_r;  // signal num: 26
wire signed [ 2 * w - 1 + 6 : 0] C7_i;  // signal num: 27
wire signed [ 2 * w - 1 + 5 : 0] C7_r;  // signal num: 28
wire signed [ 2 * w - 1 + 6 : 0] C8_i;  // signal num: 29
wire signed [ 2 * w - 1 + 5 : 0] C8_r;  // signal num: 30
wire signed [ 2 * w - 1 + 6 : 0] C9_i;  // signal num: 31
wire signed [ 2 * w - 1 + 5 : 0] C9_r;  // signal num: 32
wire signed [ 1 * w - 1 + 4 : 0] a0_i;  // signal num: 33
wire signed [ 1 * w - 1 + 4 : 0] a0_r;  // signal num: 34
wire signed [ 1 * w - 1 + 4 : 0] a10_i;  // signal num: 35
wire signed [ 1 * w - 1 + 4 : 0] a10_r;  // signal num: 36
wire signed [ 1 * w - 1 + 4 : 0] a11_i;  // signal num: 37
wire signed [ 1 * w - 1 + 4 : 0] a11_r;  // signal num: 38
wire signed [ 1 * w - 1 + 4 : 0] a12_i;  // signal num: 39
wire signed [ 1 * w - 1 + 4 : 0] a12_r;  // signal num: 40
wire signed [ 1 * w - 1 + 4 : 0] a13_i;  // signal num: 41
wire signed [ 1 * w - 1 + 4 : 0] a13_r;  // signal num: 42
wire signed [ 1 * w - 1 + 4 : 0] a14_i;  // signal num: 43
wire signed [ 1 * w - 1 + 4 : 0] a14_r;  // signal num: 44
wire signed [ 1 * w - 1 + 4 : 0] a15_i;  // signal num: 45
wire signed [ 1 * w - 1 + 4 : 0] a15_r;  // signal num: 46
wire signed [ 1 * w - 1 + 4 : 0] a16_i;  // signal num: 47
wire signed [ 1 * w - 1 + 4 : 0] a16_r;  // signal num: 48
wire signed [ 1 * w - 1 + 4 : 0] a17_i;  // signal num: 49
wire signed [ 1 * w - 1 + 4 : 0] a17_r;  // signal num: 50
wire signed [ 1 * w - 1 + 4 : 0] a18_i;  // signal num: 51
wire signed [ 1 * w - 1 + 4 : 0] a18_r;  // signal num: 52
wire signed [ 1 * w - 1 + 4 : 0] a19_i;  // signal num: 53
wire signed [ 1 * w - 1 + 4 : 0] a19_r;  // signal num: 54
wire signed [ 1 * w - 1 + 4 : 0] a1_i;  // signal num: 55
wire signed [ 1 * w - 1 + 4 : 0] a1_r;  // signal num: 56
wire signed [ 1 * w - 1 + 4 : 0] a20_i;  // signal num: 57
wire signed [ 1 * w - 1 + 4 : 0] a20_r;  // signal num: 58
wire signed [ 1 * w - 1 + 4 : 0] a21_i;  // signal num: 59
wire signed [ 1 * w - 1 + 4 : 0] a21_r;  // signal num: 60
wire signed [ 1 * w - 1 + 4 : 0] a22_i;  // signal num: 61
wire signed [ 1 * w - 1 + 4 : 0] a22_r;  // signal num: 62
wire signed [ 1 * w - 1 + 4 : 0] a23_i;  // signal num: 63
wire signed [ 1 * w - 1 + 4 : 0] a23_r;  // signal num: 64
wire signed [ 1 * w - 1 + 4 : 0] a24_i;  // signal num: 65
wire signed [ 1 * w - 1 + 4 : 0] a24_r;  // signal num: 66
wire signed [ 1 * w - 1 + 4 : 0] a25_i;  // signal num: 67
wire signed [ 1 * w - 1 + 4 : 0] a25_r;  // signal num: 68
wire signed [ 1 * w - 1 + 4 : 0] a26_i;  // signal num: 69
wire signed [ 1 * w - 1 + 4 : 0] a26_r;  // signal num: 70
wire signed [ 1 * w - 1 + 4 : 0] a27_i;  // signal num: 71
wire signed [ 1 * w - 1 + 4 : 0] a27_r;  // signal num: 72
wire signed [ 1 * w - 1 + 4 : 0] a28_i;  // signal num: 73
wire signed [ 1 * w - 1 + 4 : 0] a28_r;  // signal num: 74
wire signed [ 1 * w - 1 + 4 : 0] a29_i;  // signal num: 75
wire signed [ 1 * w - 1 + 4 : 0] a29_r;  // signal num: 76
wire signed [ 1 * w - 1 + 3 : 0] a2_i;  // signal num: 77
wire signed [ 1 * w - 1 + 3 : 0] a2_r;  // signal num: 78
wire signed [ 1 * w - 1 + 4 : 0] a30_i;  // signal num: 79
wire signed [ 1 * w - 1 + 4 : 0] a30_r;  // signal num: 80
wire signed [ 1 * w - 1 + 4 : 0] a31_i;  // signal num: 81
wire signed [ 1 * w - 1 + 4 : 0] a31_r;  // signal num: 82
wire signed [ 1 * w - 1 + 4 : 0] a32_i;  // signal num: 83
wire signed [ 1 * w - 1 + 4 : 0] a32_r;  // signal num: 84
wire signed [ 1 * w - 1 + 4 : 0] a33_i;  // signal num: 85
wire signed [ 1 * w - 1 + 4 : 0] a33_r;  // signal num: 86
wire signed [ 1 * w - 1 + 4 : 0] a34_i;  // signal num: 87
wire signed [ 1 * w - 1 + 4 : 0] a34_r;  // signal num: 88
wire signed [ 1 * w - 1 + 4 : 0] a35_i;  // signal num: 89
wire signed [ 1 * w - 1 + 4 : 0] a35_r;  // signal num: 90
wire signed [ 1 * w - 1 + 4 : 0] a36_i;  // signal num: 91
wire signed [ 1 * w - 1 + 4 : 0] a36_r;  // signal num: 92
wire signed [ 1 * w - 1 + 4 : 0] a37_i;  // signal num: 93
wire signed [ 1 * w - 1 + 4 : 0] a37_r;  // signal num: 94
wire signed [ 1 * w - 1 + 4 : 0] a38_i;  // signal num: 95
wire signed [ 1 * w - 1 + 4 : 0] a38_r;  // signal num: 96
wire signed [ 1 * w - 1 + 4 : 0] a39_i;  // signal num: 97
wire signed [ 1 * w - 1 + 4 : 0] a39_r;  // signal num: 98
wire signed [ 1 * w - 1 + 4 : 0] a3_i;  // signal num: 99
wire signed [ 1 * w - 1 + 4 : 0] a3_r;  // signal num: 100
wire signed [ 1 * w - 1 + 4 : 0] a40_i;  // signal num: 101
wire signed [ 1 * w - 1 + 4 : 0] a40_r;  // signal num: 102
wire signed [ 1 * w - 1 + 4 : 0] a41_i;  // signal num: 103
wire signed [ 1 * w - 1 + 4 : 0] a41_r;  // signal num: 104
wire signed [ 1 * w - 1 + 4 : 0] a42_i;  // signal num: 105
wire signed [ 1 * w - 1 + 4 : 0] a42_r;  // signal num: 106
wire signed [ 1 * w - 1 + 4 : 0] a43_i;  // signal num: 107
wire signed [ 1 * w - 1 + 4 : 0] a43_r;  // signal num: 108
wire signed [ 1 * w - 1 + 4 : 0] a44_i;  // signal num: 109
wire signed [ 1 * w - 1 + 4 : 0] a44_r;  // signal num: 110
wire signed [ 1 * w - 1 + 4 : 0] a45_i;  // signal num: 111
wire signed [ 1 * w - 1 + 4 : 0] a45_r;  // signal num: 112
wire signed [ 1 * w - 1 + 4 : 0] a46_i;  // signal num: 113
wire signed [ 1 * w - 1 + 4 : 0] a46_r;  // signal num: 114
wire signed [ 1 * w - 1 + 4 : 0] a47_i;  // signal num: 115
wire signed [ 1 * w - 1 + 4 : 0] a47_r;  // signal num: 116
wire signed [ 1 * w - 1 + 4 : 0] a4_i;  // signal num: 117
wire signed [ 1 * w - 1 + 4 : 0] a4_r;  // signal num: 118
wire signed [ 1 * w - 1 + 4 : 0] a5_i;  // signal num: 119
wire signed [ 1 * w - 1 + 4 : 0] a5_r;  // signal num: 120
wire signed [ 1 * w - 1 + 3 : 0] a6_i;  // signal num: 121
wire signed [ 1 * w - 1 + 3 : 0] a6_r;  // signal num: 122
wire signed [ 1 * w - 1 + 4 : 0] a7_i;  // signal num: 123
wire signed [ 1 * w - 1 + 4 : 0] a7_r;  // signal num: 124
wire signed [ 1 * w - 1 + 4 : 0] a8_i;  // signal num: 125
wire signed [ 1 * w - 1 + 4 : 0] a8_r;  // signal num: 126
wire signed [ 1 * w - 1 + 4 : 0] a9_i;  // signal num: 127
wire signed [ 1 * w - 1 + 4 : 0] a9_r;  // signal num: 128
wire signed [ 1 * w - 1 + 2 : 0] b0_i;  // signal num: 129
wire signed [ 1 * w - 1 + 2 : 0] b0_r;  // signal num: 130
wire signed [ 1 * w - 1 + 3 : 0] b10_i;  // signal num: 131
wire signed [ 1 * w - 1 + 3 : 0] b10_r;  // signal num: 132
wire signed [ 1 * w - 1 + 3 : 0] b11_i;  // signal num: 133
wire signed [ 1 * w - 1 + 3 : 0] b11_r;  // signal num: 134
wire signed [ 1 * w - 1 + 4 : 0] b12_i;  // signal num: 135
wire signed [ 1 * w - 1 + 4 : 0] b12_r;  // signal num: 136
wire signed [ 1 * w - 1 + 3 : 0] b13_i;  // signal num: 137
wire signed [ 1 * w - 1 + 3 : 0] b13_r;  // signal num: 138
wire signed [ 1 * w - 1 + 2 : 0] b14_i;  // signal num: 139
wire signed [ 1 * w - 1 + 2 : 0] b14_r;  // signal num: 140
wire signed [ 1 * w - 1 + 3 : 0] b15_i;  // signal num: 141
wire signed [ 1 * w - 1 + 3 : 0] b15_r;  // signal num: 142
wire signed [ 1 * w - 1 + 3 : 0] b16_i;  // signal num: 143
wire signed [ 1 * w - 1 + 3 : 0] b16_r;  // signal num: 144
wire signed [ 1 * w - 1 + 3 : 0] b17_i;  // signal num: 145
wire signed [ 1 * w - 1 + 3 : 0] b17_r;  // signal num: 146
wire signed [ 1 * w - 1 + 3 : 0] b18_i;  // signal num: 147
wire signed [ 1 * w - 1 + 3 : 0] b18_r;  // signal num: 148
wire signed [ 1 * w - 1 + 3 : 0] b19_i;  // signal num: 149
wire signed [ 1 * w - 1 + 3 : 0] b19_r;  // signal num: 150
wire signed [ 1 * w - 1 + 3 : 0] b1_i;  // signal num: 151
wire signed [ 1 * w - 1 + 3 : 0] b1_r;  // signal num: 152
wire signed [ 1 * w - 1 + 3 : 0] b20_i;  // signal num: 153
wire signed [ 1 * w - 1 + 3 : 0] b20_r;  // signal num: 154
wire signed [ 1 * w - 1 + 3 : 0] b21_i;  // signal num: 155
wire signed [ 1 * w - 1 + 3 : 0] b21_r;  // signal num: 156
wire signed [ 1 * w - 1 + 4 : 0] b22_i;  // signal num: 157
wire signed [ 1 * w - 1 + 4 : 0] b22_r;  // signal num: 158
wire signed [ 1 * w - 1 + 4 : 0] b23_i;  // signal num: 159
wire signed [ 1 * w - 1 + 4 : 0] b23_r;  // signal num: 160
wire signed [ 1 * w - 1 + 3 : 0] b24_i;  // signal num: 161
wire signed [ 1 * w - 1 + 3 : 0] b24_r;  // signal num: 162
wire signed [ 1 * w - 1 + 4 : 0] b25_i;  // signal num: 163
wire signed [ 1 * w - 1 + 4 : 0] b25_r;  // signal num: 164
wire signed [ 1 * w - 1 + 3 : 0] b26_i;  // signal num: 165
wire signed [ 1 * w - 1 + 3 : 0] b26_r;  // signal num: 166
wire signed [ 1 * w - 1 + 3 : 0] b27_i;  // signal num: 167
wire signed [ 1 * w - 1 + 3 : 0] b27_r;  // signal num: 168
wire signed [ 1 * w - 1 + 2 : 0] b28_i;  // signal num: 169
wire signed [ 1 * w - 1 + 2 : 0] b28_r;  // signal num: 170
wire signed [ 1 * w - 1 + 3 : 0] b29_i;  // signal num: 171
wire signed [ 1 * w - 1 + 3 : 0] b29_r;  // signal num: 172
wire signed [ 1 * w - 1 + 3 : 0] b2_i;  // signal num: 173
wire signed [ 1 * w - 1 + 3 : 0] b2_r;  // signal num: 174
wire signed [ 1 * w - 1 + 3 : 0] b30_i;  // signal num: 175
wire signed [ 1 * w - 1 + 3 : 0] b30_r;  // signal num: 176
wire signed [ 1 * w - 1 + 3 : 0] b31_i;  // signal num: 177
wire signed [ 1 * w - 1 + 3 : 0] b31_r;  // signal num: 178
wire signed [ 1 * w - 1 + 2 : 0] b32_i;  // signal num: 179
wire signed [ 1 * w - 1 + 2 : 0] b32_r;  // signal num: 180
wire signed [ 1 * w - 1 + 3 : 0] b33_i;  // signal num: 181
wire signed [ 1 * w - 1 + 3 : 0] b33_r;  // signal num: 182
wire signed [ 1 * w - 1 + 3 : 0] b34_i;  // signal num: 183
wire signed [ 1 * w - 1 + 3 : 0] b34_r;  // signal num: 184
wire signed [ 1 * w - 1 + 3 : 0] b35_i;  // signal num: 185
wire signed [ 1 * w - 1 + 3 : 0] b35_r;  // signal num: 186
wire signed [ 1 * w - 1 + 4 : 0] b36_i;  // signal num: 187
wire signed [ 1 * w - 1 + 4 : 0] b36_r;  // signal num: 188
wire signed [ 1 * w - 1 + 4 : 0] b37_i;  // signal num: 189
wire signed [ 1 * w - 1 + 4 : 0] b37_r;  // signal num: 190
wire signed [ 1 * w - 1 + 2 : 0] b38_i;  // signal num: 191
wire signed [ 1 * w - 1 + 2 : 0] b38_r;  // signal num: 192
wire signed [ 1 * w - 1 + 3 : 0] b39_i;  // signal num: 193
wire signed [ 1 * w - 1 + 3 : 0] b39_r;  // signal num: 194
wire signed [ 1 * w - 1 + 3 : 0] b3_i;  // signal num: 195
wire signed [ 1 * w - 1 + 3 : 0] b3_r;  // signal num: 196
wire signed [ 1 * w - 1 + 3 : 0] b40_i;  // signal num: 197
wire signed [ 1 * w - 1 + 3 : 0] b40_r;  // signal num: 198
wire signed [ 1 * w - 1 + 3 : 0] b41_i;  // signal num: 199
wire signed [ 1 * w - 1 + 3 : 0] b41_r;  // signal num: 200
wire signed [ 1 * w - 1 + 2 : 0] b42_i;  // signal num: 201
wire signed [ 1 * w - 1 + 2 : 0] b42_r;  // signal num: 202
wire signed [ 1 * w - 1 + 4 : 0] b43_i;  // signal num: 203
wire signed [ 1 * w - 1 + 4 : 0] b43_r;  // signal num: 204
wire signed [ 1 * w - 1 + 2 : 0] b44_i;  // signal num: 205
wire signed [ 1 * w - 1 + 2 : 0] b44_r;  // signal num: 206
wire signed [ 1 * w - 1 + 4 : 0] b45_i;  // signal num: 207
wire signed [ 1 * w - 1 + 4 : 0] b45_r;  // signal num: 208
wire signed [ 1 * w - 1 + 3 : 0] b46_i;  // signal num: 209
wire signed [ 1 * w - 1 + 3 : 0] b46_r;  // signal num: 210
wire signed [ 1 * w - 1 + 2 : 0] b47_i;  // signal num: 211
wire signed [ 1 * w - 1 + 2 : 0] b47_r;  // signal num: 212
wire signed [ 1 * w - 1 + 4 : 0] b4_i;  // signal num: 213
wire signed [ 1 * w - 1 + 4 : 0] b4_r;  // signal num: 214
wire signed [ 1 * w - 1 + 3 : 0] b5_i;  // signal num: 215
wire signed [ 1 * w - 1 + 3 : 0] b5_r;  // signal num: 216
wire signed [ 1 * w - 1 + 3 : 0] b6_i;  // signal num: 217
wire signed [ 1 * w - 1 + 3 : 0] b6_r;  // signal num: 218
wire signed [ 1 * w - 1 + 3 : 0] b7_i;  // signal num: 219
wire signed [ 1 * w - 1 + 3 : 0] b7_r;  // signal num: 220
wire signed [ 1 * w - 1 + 3 : 0] b8_i;  // signal num: 221
wire signed [ 1 * w - 1 + 3 : 0] b8_r;  // signal num: 222
wire signed [ 1 * w - 1 + 4 : 0] b9_i;  // signal num: 223
wire signed [ 1 * w - 1 + 4 : 0] b9_r;  // signal num: 224
wire signed [ 2 * w - 1 + 6 : 0] k00;  // signal num: 225
wire signed [ 2 * w - 1 + 7 : 0] k01;  // signal num: 226
wire signed [ 2 * w - 1 + 6 : 0] k02;  // signal num: 227
wire signed [ 2 * w - 1 + 7 : 0] k03;  // signal num: 228
wire signed [ 2 * w - 1 + 7 : 0] k04;  // signal num: 229
wire signed [ 2 * w - 1 + 7 : 0] k05;  // signal num: 230
wire signed [ 2 * w - 1 + 6 : 0] k06;  // signal num: 231
wire signed [ 2 * w - 1 + 7 : 0] k07;  // signal num: 232
wire signed [ 2 * w - 1 + 7 : 0] k08;  // signal num: 233
wire signed [ 2 * w - 1 + 7 : 0] k09;  // signal num: 234
wire signed [ 2 * w - 1 + 7 : 0] k10;  // signal num: 235
wire signed [ 2 * w - 1 + 7 : 0] k11;  // signal num: 236
wire signed [ 2 * w - 1 + 7 : 0] k12;  // signal num: 237
wire signed [ 2 * w - 1 + 7 : 0] k13;  // signal num: 238
wire signed [ 2 * w - 1 + 6 : 0] k14;  // signal num: 239
wire signed [ 2 * w - 1 + 7 : 0] k15;  // signal num: 240
wire signed [ 2 * w - 1 + 7 : 0] k16;  // signal num: 241
wire signed [ 2 * w - 1 + 7 : 0] k17;  // signal num: 242
wire signed [ 2 * w - 1 + 7 : 0] k18;  // signal num: 243
wire signed [ 2 * w - 1 + 7 : 0] k19;  // signal num: 244
wire signed [ 2 * w - 1 + 7 : 0] k20;  // signal num: 245
wire signed [ 2 * w - 1 + 7 : 0] k21;  // signal num: 246
wire signed [ 2 * w - 1 + 7 : 0] k22;  // signal num: 247
wire signed [ 2 * w - 1 + 7 : 0] k23;  // signal num: 248
wire signed [ 2 * w - 1 + 7 : 0] k26;  // signal num: 249
wire signed [ 2 * w - 1 + 7 : 0] k27;  // signal num: 250
wire signed [ 2 * w - 1 + 6 : 0] k28;  // signal num: 251
wire signed [ 2 * w - 1 + 7 : 0] k29;  // signal num: 252
wire signed [ 2 * w - 1 + 7 : 0] k30;  // signal num: 253
wire signed [ 2 * w - 1 + 7 : 0] k31;  // signal num: 254
wire signed [ 2 * w - 1 + 6 : 0] k32;  // signal num: 255
wire signed [ 2 * w - 1 + 7 : 0] k35;  // signal num: 256
wire signed [ 2 * w - 1 + 7 : 0] k36;  // signal num: 257
wire signed [ 2 * w - 1 + 6 : 0] k38;  // signal num: 258
wire signed [ 2 * w - 1 + 7 : 0] k39;  // signal num: 259
wire signed [ 2 * w - 1 + 7 : 0] k40;  // signal num: 260
wire signed [ 2 * w - 1 + 7 : 0] k41;  // signal num: 261
wire signed [ 2 * w - 1 + 6 : 0] k42;  // signal num: 262
wire signed [ 2 * w - 1 + 6 : 0] k44;  // signal num: 263
wire signed [ 2 * w - 1 + 7 : 0] k45;  // signal num: 264
wire signed [ 2 * w - 1 + 7 : 0] k46;  // signal num: 265
wire signed [ 1 * w - 1 + 2 : 0] mAi13mAr13mAi1;  // signal num: 266
wire signed [ 1 * w - 1 + 3 : 0] mAi13mAr13mAi1mAr6pAr1pAr10;  // signal num: 267
wire signed [ 1 * w - 1 + 2 : 0] mAi5mAi9;  // signal num: 268
wire signed [ 1 * w - 1 + 2 : 0] mAi5mAi9pAr5pAr9;  // signal num: 269
wire signed [ 1 * w - 1 + 2 : 0] mAi6mAi10;  // signal num: 270
wire signed [ 1 * w - 1 + 3 : 0] mAi6mAi10mAr1mAi13;  // signal num: 271
wire signed [ 1 * w - 1 + 3 : 0] mAi6mAi10mAr6mAr10;  // signal num: 272
wire signed [ 1 * w - 1 + 3 : 0] mAi6mAi10pAr1pAr10pAr6mAi13;  // signal num: 273
wire signed [ 1 * w - 1 + 2 : 0] mAi7mAi11;  // signal num: 274
wire signed [ 1 * w - 1 + 2 : 0] mAi8mAi12;  // signal num: 275
wire signed [ 1 * w - 1 + 3 : 0] mAi8mAi12mAr8mAr12;  // signal num: 276
wire signed [ 1 * w - 1 + 2 : 0] mAi8mAi12pAr1mAi13;  // signal num: 277
wire signed [ 1 * w - 1 + 2 : 0] mAr10pAr6mAi13;  // signal num: 278
wire signed [ 1 * w - 1 + 2 : 0] mAr13mAi1;  // signal num: 279
wire signed [ 1 * w - 1 + 3 : 0] mAr13mAi1mAi6mAi10pAr1pAr10pAr6mAi13;  // signal num: 280
wire signed [ 1 * w - 1 + 2 : 0] mAr13mAi1pAr10mAr6;  // signal num: 281
wire signed [ 1 * w - 1 + 3 : 0] mAr13mAi1pAr10mAr6pAi6mAi10mAr1pAi13;  // signal num: 282
wire signed [ 1 * w - 1 + 1 : 0] mAr13pAi1;  // signal num: 283
wire signed [ 1 * w - 1 + 2 : 0] mAr13pAi1pAi10mAi6;  // signal num: 284
wire signed [ 1 * w - 1 + 2 : 0] mAr13pAi1pAr1pAi13;  // signal num: 285
wire signed [ 1 * w - 1 + 3 : 0] mAr13pAi1pAr1pAi13pAi12mAi8pAr8mAr12;  // signal num: 286
wire signed [ 1 * w - 1 + 3 : 0] mAr13pAi1pAr8mAr12pAr3mAi15;  // signal num: 287
wire signed [ 1 * w - 1 + 2 : 0] mAr14mAi2;  // signal num: 288
wire signed [ 1 * w - 1 + 2 : 0] mAr14mAi2pAr2mAi14;  // signal num: 289
wire signed [ 1 * w - 1 + 1 : 0] mAr14pAi2;  // signal num: 290
wire signed [ 1 * w - 1 + 2 : 0] mAr14pAi2mAr16mAi4;  // signal num: 291
wire signed [ 1 * w - 1 + 2 : 0] mAr14pAi2mAr2mAi14;  // signal num: 292
wire signed [ 1 * w - 1 + 2 : 0] mAr15mAi3;  // signal num: 293
wire signed [ 1 * w - 1 + 1 : 0] mAr15pAi3;  // signal num: 294
wire signed [ 1 * w - 1 + 2 : 0] mAr15pAi3pAr3pAi15;  // signal num: 295
wire signed [ 1 * w - 1 + 3 : 0] mAr15pAi3pAr3pAi15pAr12mAr8pAi8mAi12;  // signal num: 296
wire signed [ 1 * w - 1 + 2 : 0] mAr16mAi4;  // signal num: 297
wire signed [ 1 * w - 1 + 1 : 0] mAr16pAi4;  // signal num: 298
wire signed [ 1 * w - 1 + 2 : 0] mAr16pAi4pAr4pAi16;  // signal num: 299
wire signed [ 1 * w - 1 + 2 : 0] mAr1mAi13;  // signal num: 300
wire signed [ 1 * w - 1 + 2 : 0] mAr1mAr10pAr6mAi13;  // signal num: 301
wire signed [ 1 * w - 1 + 1 : 0] mAr1pAi13;  // signal num: 302
wire signed [ 1 * w - 1 + 2 : 0] mAr2mAi14;  // signal num: 303
wire signed [ 1 * w - 1 + 1 : 0] mAr2pAi14;  // signal num: 304
wire signed [ 1 * w - 1 + 2 : 0] mAr2pAi14mAr14mAi2;  // signal num: 305
wire signed [ 1 * w - 1 + 2 : 0] mAr2pAi14pAr14pAi2;  // signal num: 306
wire signed [ 1 * w - 1 + 2 : 0] mAr3mAi15;  // signal num: 307
wire signed [ 1 * w - 1 + 3 : 0] mAr3mAi15mAi6mAi10;  // signal num: 308
wire signed [ 1 * w - 1 + 2 : 0] mAr3mAi15mAr15pAi3;  // signal num: 309
wire signed [ 1 * w - 1 + 3 : 0] mAr3mAi15mAr15pAi3mAi6mAi10mAr6mAr10;  // signal num: 310
wire signed [ 1 * w - 1 + 3 : 0] mAr3mAi15mAr15pAi3pAr11mAr7pAi7mAi11;  // signal num: 311
wire signed [ 1 * w - 1 + 3 : 0] mAr3mAi15mAr15pAi3pAr8pAr12mAi8mAi12;  // signal num: 312
wire signed [ 1 * w - 1 + 2 : 0] mAr3mAi15pAr6mAr10;  // signal num: 313
wire signed [ 1 * w - 1 + 1 : 0] mAr3pAi15;  // signal num: 314
wire signed [ 1 * w - 1 + 2 : 0] mAr3pAi15mAr15mAi3;  // signal num: 315
wire signed [ 1 * w - 1 + 2 : 0] mAr3pAi15pAr12mAr8;  // signal num: 316
wire signed [ 1 * w - 1 + 2 : 0] mAr4mAi16;  // signal num: 317
wire signed [ 1 * w - 1 + 2 : 0] mAr4mAi16mAr16pAi4;  // signal num: 318
wire signed [ 1 * w - 1 + 2 : 0] mAr4mAi16mAr2pAi14;  // signal num: 319
wire signed [ 1 * w - 1 + 3 : 0] mAr4mAi16mAr2pAi14pAi11mAi7pAr5pAr9;  // signal num: 320
wire signed [ 1 * w - 1 + 3 : 0] mAr4mAi16mAr2pAi14pAi9mAi5pAr7pAr11;  // signal num: 321
wire signed [ 1 * w - 1 + 2 : 0] mAr4mAi16pAr16mAi4;  // signal num: 322
wire signed [ 1 * w - 1 + 3 : 0] mAr4mAi16pAr16mAi4pAi12mAi8pAr12mAr8;  // signal num: 323
wire signed [ 1 * w - 1 + 1 : 0] mAr4pAi16;  // signal num: 324
wire signed [ 1 * w - 1 + 2 : 0] mAr4pAi16mAr16mAi4;  // signal num: 325
wire signed [ 1 * w - 1 + 2 : 0] mAr5mAr9;  // signal num: 326
wire signed [ 1 * w - 1 + 3 : 0] mAr5mAr9mAi5mAi9;  // signal num: 327
wire signed [ 1 * w - 1 + 2 : 0] mAr5mAr9pAi5pAi9;  // signal num: 328
wire signed [ 1 * w - 1 + 3 : 0] mAr5mAr9pAi5pAi9pAr16mAi4pAr4pAi16;  // signal num: 329
wire signed [ 1 * w - 1 + 2 : 0] mAr5mAr9pAi7mAi11;  // signal num: 330
wire signed [ 1 * w - 1 + 3 : 0] mAr5mAr9pAi7mAi11mAr14pAi2mAr16mAi4;  // signal num: 331
wire signed [ 1 * w - 1 + 3 : 0] mAr5mAr9pAi7mAi11pAr16pAi4pAr14mAi2;  // signal num: 332
wire signed [ 1 * w - 1 + 2 : 0] mAr6mAr10;  // signal num: 333
wire signed [ 1 * w - 1 + 3 : 0] mAr6mAr10mAr13pAi1pAr8mAr12pAr3mAi15;  // signal num: 334
wire signed [ 1 * w - 1 + 2 : 0] mAr6pAi13pAr1pAr10;  // signal num: 335
wire signed [ 1 * w - 1 + 2 : 0] mAr6pAr1pAr10;  // signal num: 336
wire signed [ 1 * w - 1 + 2 : 0] mAr7mAr11;  // signal num: 337
wire signed [ 1 * w - 1 + 3 : 0] mAr7mAr11mAi7mAi11;  // signal num: 338
wire signed [ 1 * w - 1 + 3 : 0] mAr7mAr11mAi7mAi11pAr14pAi2pAr2mAi14;  // signal num: 339
wire signed [ 1 * w - 1 + 2 : 0] mAr7mAr11pAi5mAi9;  // signal num: 340
wire signed [ 1 * w - 1 + 2 : 0] mAr8mAr12;  // signal num: 341
wire signed [ 1 * w - 1 + 2 : 0] mBi10mBi12;  // signal num: 342
wire signed [ 1 * w - 1 + 2 : 0] mBi10mBr15mBi11;  // signal num: 343
wire signed [ 1 * w - 1 + 2 : 0] mBi12pBi9mBr15mBi11;  // signal num: 344
wire signed [ 1 * w - 1 + 2 : 0] mBi13pBr9mBr11pBi15;  // signal num: 345
wire signed [ 1 * w - 1 + 2 : 0] mBi14mBi16;  // signal num: 346
wire signed [ 1 * w - 1 + 2 : 0] mBi14mBi16pBr1mBr5;  // signal num: 347
wire signed [ 1 * w - 1 + 2 : 0] mBi16pBi13pBr11mBi15;  // signal num: 348
wire signed [ 1 * w - 1 + 2 : 0] mBi1mBi5;  // signal num: 349
wire signed [ 1 * w - 1 + 2 : 0] mBi2mBi6;  // signal num: 350
wire signed [ 1 * w - 1 + 2 : 0] mBi2mBr16mBr14pBi5;  // signal num: 351
wire signed [ 1 * w - 1 + 2 : 0] mBi3mBi7;  // signal num: 352
wire signed [ 1 * w - 1 + 2 : 0] mBi4mBi8;  // signal num: 353
wire signed [ 1 * w - 1 + 3 : 0] mBi4mBi8mBi3mBi7;  // signal num: 354
wire signed [ 1 * w - 1 + 3 : 0] mBi4mBi8mBi3mBi7mBi2mBi6;  // signal num: 355
wire signed [ 1 * w - 1 + 3 : 0] mBi4mBi8mBi3mBi7pBi1pBi5;  // signal num: 356
wire signed [ 1 * w - 1 + 2 : 0] mBi4mBi8pBi1pBi5;  // signal num: 357
wire signed [ 1 * w - 1 + 2 : 0] mBi5pBi7mBi3;  // signal num: 358
wire signed [ 1 * w - 1 + 2 : 0] mBi6pBr16mBr13;  // signal num: 359
wire signed [ 1 * w - 1 + 2 : 0] mBr10mBr12;  // signal num: 360
wire signed [ 1 * w - 1 + 3 : 0] mBr10mBr12pBr2mBi14mBi16;  // signal num: 361
wire signed [ 1 * w - 1 + 2 : 0] mBr10pBr2pBi14;  // signal num: 362
wire signed [ 1 * w - 1 + 2 : 0] mBr11mBi15;  // signal num: 363
wire signed [ 1 * w - 1 + 1 : 0] mBr11pBi15;  // signal num: 364
wire signed [ 1 * w - 1 + 2 : 0] mBr11pBi15mBr10mBr12;  // signal num: 365
wire signed [ 1 * w - 1 + 3 : 0] mBr11pBi15mBr10pBr2pBi14;  // signal num: 366
wire signed [ 1 * w - 1 + 3 : 0] mBr12pBi16mBi13pBr9mBr11pBi15;  // signal num: 367
wire signed [ 1 * w - 1 + 1 : 0] mBr13pBi9;  // signal num: 368
wire signed [ 1 * w - 1 + 3 : 0] mBr14mBi10mBr15mBi11;  // signal num: 369
wire signed [ 1 * w - 1 + 3 : 0] mBr14mBi10mBr15mBi11mBr16mBi12;  // signal num: 370
wire signed [ 1 * w - 1 + 1 : 0] mBr14pBi10;  // signal num: 371
wire signed [ 1 * w - 1 + 2 : 0] mBr14pBi10mBr15pBi11;  // signal num: 372
wire signed [ 1 * w - 1 + 1 : 0] mBr14pBi5;  // signal num: 373
wire signed [ 1 * w - 1 + 2 : 0] mBr15mBi11;  // signal num: 374
wire signed [ 1 * w - 1 + 1 : 0] mBr15pBi11;  // signal num: 375
wire signed [ 1 * w - 1 + 2 : 0] mBr16mBi12;  // signal num: 376
wire signed [ 1 * w - 1 + 2 : 0] mBr16mBr14pBi5;  // signal num: 377
wire signed [ 1 * w - 1 + 1 : 0] mBr16pBi12;  // signal num: 378
wire signed [ 1 * w - 1 + 2 : 0] mBr1mBr5;  // signal num: 379
wire signed [ 1 * w - 1 + 3 : 0] mBr1mBr5pBr4pBr8pBr3pBr7;  // signal num: 380
wire signed [ 1 * w - 1 + 2 : 0] mBr2mBi14;  // signal num: 381
wire signed [ 1 * w - 1 + 1 : 0] mBr2pBi14;  // signal num: 382
wire signed [ 1 * w - 1 + 2 : 0] mBr3mBr7;  // signal num: 383
wire signed [ 1 * w - 1 + 2 : 0] mBr4mBr8;  // signal num: 384
wire signed [ 1 * w - 1 + 3 : 0] mBr4mBr8pBr12pBr10mBr6;  // signal num: 385
wire signed [ 1 * w - 1 + 2 : 0] mBr6mBr4mBr8;  // signal num: 386
wire signed [ 1 * w - 1 + 2 : 0] mBr9mBi13;  // signal num: 387
wire signed [ 1 * w - 1 + 2 : 0] mBr9pBi13pBr11mBi15;  // signal num: 388
wire signed [ 2 * w - 1 + 6 : 0] mt00_ipt34_rmt11_r;  // signal num: 389
wire signed [ 2 * w - 1 + 7 : 0] mt03_imt12_i;  // signal num: 390
wire signed [ 2 * w - 1 + 6 : 0] mt03_rpt13_rpt44_i;  // signal num: 391
wire signed [ 2 * w - 1 + 7 : 0] mt04_ipt45_imt29_i;  // signal num: 392
wire signed [ 2 * w - 1 + 6 : 0] mt06_imt10_i;  // signal num: 393
wire signed [ 2 * w - 1 + 6 : 0] mt06_rmt10_r;  // signal num: 394
wire signed [ 2 * w - 1 + 6 : 0] mt08_imt47_r;  // signal num: 395
wire signed [ 2 * w - 1 + 7 : 0] mt08_imt47_rmt03_imt12_i;  // signal num: 396
wire signed [ 2 * w - 1 + 7 : 0] mt08_imt47_rpt24_imt11_ipt37_ipt39_ipt14_ipt27_rpt03_rmt31_r;  // signal num: 397
wire signed [ 2 * w - 1 + 6 : 0] mt13_rmt44_i;  // signal num: 398
wire signed [ 2 * w - 1 + 7 : 0] mt14_imt43_ipt34_rmt11_r;  // signal num: 399
wire signed [ 2 * w - 1 + 7 : 0] mt16_imt18_ipt21_rpt34_i;  // signal num: 400
wire signed [ 2 * w - 1 + 7 : 0] mt16_rmt18_rmt24_rmt37_rpt32_imt26_i;  // signal num: 401
wire signed [ 2 * w - 1 + 7 : 0] mt17_rmt22_i;  // signal num: 402
wire signed [ 2 * w - 1 + 7 : 0] mt17_rmt22_ipt03_imt14_rmt34_ipt18_imt43_r;  // signal num: 403
wire signed [ 2 * w - 1 + 6 : 0] mt18_rmt27_i;  // signal num: 404
wire signed [ 2 * w - 1 + 6 : 0] mt19_imt35_i;  // signal num: 405
wire signed [ 2 * w - 1 + 7 : 0] mt22_rpt17_ipt18_r;  // signal num: 406
wire signed [ 2 * w - 1 + 7 : 0] mt22_rpt17_ipt18_rmt37_rpt27_imt33_r;  // signal num: 407
wire signed [ 2 * w - 1 + 7 : 0] mt24_rmt37_r;  // signal num: 408
wire signed [ 2 * w - 1 + 6 : 0] mt28_rmt40_i;  // signal num: 409
wire signed [ 2 * w - 1 + 7 : 0] mt29_rmt03_rpt13_rpt44_imt18_rmt27_ipt33_rpt37_r;  // signal num: 410
wire signed [ 2 * w - 1 + 7 : 0] mt31_imt24_rmt37_r;  // signal num: 411
wire signed [ 2 * w - 1 + 6 : 0] mt33_imt47_i;  // signal num: 412
wire signed [ 2 * w - 1 + 7 : 0] mt34_imt38_rpt36_rmt05_r;  // signal num: 413
wire signed [ 2 * w - 1 + 6 : 0] mt34_imt47_i;  // signal num: 414
wire signed [ 2 * w - 1 + 7 : 0] mt34_imt47_ipt03_rmt31_rpt11_imt08_rpt43_rmt39_i;  // signal num: 415
wire signed [ 2 * w - 1 + 6 : 0] mt34_rmt47_r;  // signal num: 416
wire signed [ 2 * w - 1 + 7 : 0] mt36_ipt30_rmt38_i;  // signal num: 417
wire signed [ 2 * w - 1 + 7 : 0] mt36_ipt30_rmt38_ipt32_rpt05_ipt18_rpt16_rmt26_r;  // signal num: 418
wire signed [ 2 * w - 1 + 7 : 0] mt37_imt24_imt16_imt18_i;  // signal num: 419
wire signed [ 2 * w - 1 + 7 : 0] mt37_rpt27_imt33_r;  // signal num: 420
wire signed [ 2 * w - 1 + 7 : 0] mt38_rpt36_rmt05_r;  // signal num: 421
wire signed [ 2 * w - 1 + 6 : 0] mt41_ipt00_imt14_i;  // signal num: 422
wire signed [ 2 * w - 1 + 6 : 0] mt42_imt46_r;  // signal num: 423
wire signed [ 2 * w - 1 + 5 : 0] p00;  // signal num: 424
wire signed [ 2 * w - 1 + 6 : 0] p01;  // signal num: 425
wire signed [ 2 * w - 1 + 5 : 0] p02;  // signal num: 426
wire signed [ 2 * w - 1 + 6 : 0] p03;  // signal num: 427
wire signed [ 2 * w - 1 + 6 : 0] p05;  // signal num: 428
wire signed [ 2 * w - 1 + 5 : 0] p06;  // signal num: 429
wire signed [ 2 * w - 1 + 6 : 0] p07;  // signal num: 430
wire signed [ 2 * w - 1 + 6 : 0] p08;  // signal num: 431
wire signed [ 2 * w - 1 + 6 : 0] p09;  // signal num: 432
wire signed [ 2 * w - 1 + 6 : 0] p10;  // signal num: 433
wire signed [ 2 * w - 1 + 6 : 0] p11;  // signal num: 434
wire signed [ 2 * w - 1 + 6 : 0] p13;  // signal num: 435
wire signed [ 2 * w - 1 + 5 : 0] p14;  // signal num: 436
wire signed [ 2 * w - 1 + 6 : 0] p15;  // signal num: 437
wire signed [ 2 * w - 1 + 6 : 0] p16;  // signal num: 438
wire signed [ 2 * w - 1 + 6 : 0] p17;  // signal num: 439
wire signed [ 2 * w - 1 + 6 : 0] p18;  // signal num: 440
wire signed [ 2 * w - 1 + 6 : 0] p19;  // signal num: 441
wire signed [ 2 * w - 1 + 6 : 0] p20;  // signal num: 442
wire signed [ 2 * w - 1 + 6 : 0] p21;  // signal num: 443
wire signed [ 2 * w - 1 + 6 : 0] p24;  // signal num: 444
wire signed [ 2 * w - 1 + 6 : 0] p26;  // signal num: 445
wire signed [ 2 * w - 1 + 6 : 0] p27;  // signal num: 446
wire signed [ 2 * w - 1 + 5 : 0] p28;  // signal num: 447
wire signed [ 2 * w - 1 + 6 : 0] p29;  // signal num: 448
wire signed [ 2 * w - 1 + 6 : 0] p30;  // signal num: 449
wire signed [ 2 * w - 1 + 6 : 0] p31;  // signal num: 450
wire signed [ 2 * w - 1 + 5 : 0] p32;  // signal num: 451
wire signed [ 2 * w - 1 + 6 : 0] p33;  // signal num: 452
wire signed [ 2 * w - 1 + 6 : 0] p34;  // signal num: 453
wire signed [ 2 * w - 1 + 6 : 0] p35;  // signal num: 454
wire signed [ 2 * w - 1 + 5 : 0] p38;  // signal num: 455
wire signed [ 2 * w - 1 + 6 : 0] p39;  // signal num: 456
wire signed [ 2 * w - 1 + 6 : 0] p40;  // signal num: 457
wire signed [ 2 * w - 1 + 6 : 0] p41;  // signal num: 458
wire signed [ 2 * w - 1 + 5 : 0] p42;  // signal num: 459
wire signed [ 2 * w - 1 + 6 : 0] p43;  // signal num: 460
wire signed [ 2 * w - 1 + 5 : 0] p44;  // signal num: 461
wire signed [ 2 * w - 1 + 6 : 0] p46;  // signal num: 462
wire signed [ 2 * w - 1 + 5 : 0] p47;  // signal num: 463
wire signed [ 1 * w - 1 + 1 : 0] pAi10mAi6;  // signal num: 464
wire signed [ 1 * w - 1 + 3 : 0] pAi10mAi6mAi13mAr13mAi1mAr6pAr1pAr10;  // signal num: 465
wire signed [ 1 * w - 1 + 3 : 0] pAi10mAi6pAr11mAr7pAi11mAi7;  // signal num: 466
wire signed [ 1 * w - 1 + 2 : 0] pAi10mAi6pAr6mAr10;  // signal num: 467
wire signed [ 1 * w - 1 + 1 : 0] pAi11mAi7;  // signal num: 468
wire signed [ 1 * w - 1 + 2 : 0] pAi11mAi7pAr5pAr9;  // signal num: 469
wire signed [ 1 * w - 1 + 2 : 0] pAi11mAi7pAr7mAr11;  // signal num: 470
wire signed [ 1 * w - 1 + 1 : 0] pAi12mAi8;  // signal num: 471
wire signed [ 1 * w - 1 + 2 : 0] pAi12mAi8mAr13mAi1;  // signal num: 472
wire signed [ 1 * w - 1 + 3 : 0] pAi12mAi8mAr13mAi1mAr3mAi15mAi6mAi10;  // signal num: 473
wire signed [ 1 * w - 1 + 2 : 0] pAi12mAi8pAr12mAr8;  // signal num: 474
wire signed [ 1 * w - 1 + 3 : 0] pAi12mAi8pAr12mAr8mAr16pAi4pAr4pAi16;  // signal num: 475
wire signed [ 1 * w - 1 + 2 : 0] pAi12mAi8pAr8mAr12;  // signal num: 476
wire signed [ 1 * w - 1 + 2 : 0] pAi13pAr1pAr10;  // signal num: 477
wire signed [ 1 * w - 1 + 3 : 0] pAi13pAr1pAr10pAr13mAi1pAi6pAi10;  // signal num: 478
wire signed [ 1 * w - 1 + 1 : 0] pAi5mAi9;  // signal num: 479
wire signed [ 1 * w - 1 + 2 : 0] pAi5mAi9pAr9mAr5;  // signal num: 480
wire signed [ 1 * w - 1 + 1 : 0] pAi5pAi9;  // signal num: 481
wire signed [ 1 * w - 1 + 2 : 0] pAi5pAi9pAr5pAr9;  // signal num: 482
wire signed [ 1 * w - 1 + 3 : 0] pAi5pAi9pAr5pAr9mAr4mAi16pAr16mAi4;  // signal num: 483
wire signed [ 1 * w - 1 + 2 : 0] pAi5pAi9pAr7mAr11;  // signal num: 484
wire signed [ 1 * w - 1 + 1 : 0] pAi6mAi10;  // signal num: 485
wire signed [ 1 * w - 1 + 2 : 0] pAi6mAi10mAr1pAi13;  // signal num: 486
wire signed [ 1 * w - 1 + 2 : 0] pAi6mAi10pAr15mAi3;  // signal num: 487
wire signed [ 1 * w - 1 + 1 : 0] pAi6pAi10;  // signal num: 488
wire signed [ 1 * w - 1 + 2 : 0] pAi6pAi10mAr6mAr10;  // signal num: 489
wire signed [ 1 * w - 1 + 3 : 0] pAi6pAi10mAr6mAr10mAr13pAi1pAr1pAi13;  // signal num: 490
wire signed [ 1 * w - 1 + 3 : 0] pAi6pAi10mAr6mAr10pAr3mAi15mAr15mAi3;  // signal num: 491
wire signed [ 1 * w - 1 + 2 : 0] pAi6pAi10pAr13pAi1;  // signal num: 492
wire signed [ 1 * w - 1 + 1 : 0] pAi7mAi11;  // signal num: 493
wire signed [ 1 * w - 1 + 1 : 0] pAi7pAi11;  // signal num: 494
wire signed [ 1 * w - 1 + 2 : 0] pAi7pAi11mAr7mAr11;  // signal num: 495
wire signed [ 1 * w - 1 + 3 : 0] pAi7pAi11mAr7mAr11mAi8mAi12mAr8mAr12;  // signal num: 496
wire signed [ 1 * w - 1 + 2 : 0] pAi7pAi11pAr5mAr9;  // signal num: 497
wire signed [ 1 * w - 1 + 2 : 0] pAi7pAi11pAr7pAr11;  // signal num: 498
wire signed [ 1 * w - 1 + 3 : 0] pAi7pAi11pAr7pAr11pAr4mAi16mAr16mAi4;  // signal num: 499
wire signed [ 1 * w - 1 + 1 : 0] pAi8mAi12;  // signal num: 500
wire signed [ 1 * w - 1 + 1 : 0] pAi8pAi12;  // signal num: 501
wire signed [ 1 * w - 1 + 2 : 0] pAi8pAi12mAr1pAi13;  // signal num: 502
wire signed [ 1 * w - 1 + 3 : 0] pAi8pAi12mAr1pAi13pAi6mAi10pAr15mAi3;  // signal num: 503
wire signed [ 1 * w - 1 + 2 : 0] pAi8pAi12mAr3pAi15;  // signal num: 504
wire signed [ 1 * w - 1 + 3 : 0] pAi8pAi12mAr3pAi15mAr13pAi1pAi10mAi6;  // signal num: 505
wire signed [ 1 * w - 1 + 1 : 0] pAi9mAi5;  // signal num: 506
wire signed [ 1 * w - 1 + 2 : 0] pAi9mAi5pAr7pAr11;  // signal num: 507
wire signed [ 1 * w - 1 + 2 : 0] pAi9mAi5pAr9mAr5;  // signal num: 508
wire signed [ 1 * w - 1 + 1 : 0] pAr10mAr6;  // signal num: 509
wire signed [ 1 * w - 1 + 1 : 0] pAr11mAr7;  // signal num: 510
wire signed [ 1 * w - 1 + 2 : 0] pAr11mAr7mAi5mAi9;  // signal num: 511
wire signed [ 1 * w - 1 + 2 : 0] pAr11mAr7pAi11mAi7;  // signal num: 512
wire signed [ 1 * w - 1 + 3 : 0] pAr11mAr7pAi11mAi7pAr16pAi4pAr4mAi16;  // signal num: 513
wire signed [ 1 * w - 1 + 2 : 0] pAr11mAr7pAi7mAi11;  // signal num: 514
wire signed [ 1 * w - 1 + 1 : 0] pAr12mAr8;  // signal num: 515
wire signed [ 1 * w - 1 + 2 : 0] pAr12mAr8pAi8mAi12;  // signal num: 516
wire signed [ 1 * w - 1 + 1 : 0] pAr13mAi1;  // signal num: 517
wire signed [ 1 * w - 1 + 2 : 0] pAr13mAi1pAi6mAi10;  // signal num: 518
wire signed [ 1 * w - 1 + 3 : 0] pAr13mAi1pAi6mAi10mAr1mAr10pAr6mAi13;  // signal num: 519
wire signed [ 1 * w - 1 + 2 : 0] pAr13mAi1pAi6pAi10;  // signal num: 520
wire signed [ 1 * w - 1 + 2 : 0] pAr13mAi1pAr1pAi13;  // signal num: 521
wire signed [ 1 * w - 1 + 3 : 0] pAr13mAi1pAr6pAr10mAr3pAi15pAr12mAr8;  // signal num: 522
wire signed [ 1 * w - 1 + 1 : 0] pAr13pAi1;  // signal num: 523
wire signed [ 1 * w - 1 + 1 : 0] pAr14mAi2;  // signal num: 524
wire signed [ 1 * w - 1 + 2 : 0] pAr14mAi2mAr2mAi14;  // signal num: 525
wire signed [ 1 * w - 1 + 1 : 0] pAr14pAi2;  // signal num: 526
wire signed [ 1 * w - 1 + 2 : 0] pAr14pAi2mAr16pAi4;  // signal num: 527
wire signed [ 1 * w - 1 + 2 : 0] pAr14pAi2pAr2mAi14;  // signal num: 528
wire signed [ 1 * w - 1 + 1 : 0] pAr15mAi3;  // signal num: 529
wire signed [ 1 * w - 1 + 3 : 0] pAr15mAi3pAi7pAi11mAr7mAr11;  // signal num: 530
wire signed [ 1 * w - 1 + 3 : 0] pAr15mAi3pAr8mAr12pAr1pAr10pAr6mAi13;  // signal num: 531
wire signed [ 1 * w - 1 + 1 : 0] pAr15pAi3;  // signal num: 532
wire signed [ 1 * w - 1 + 2 : 0] pAr15pAi3pAi8mAi12;  // signal num: 533
wire signed [ 1 * w - 1 + 3 : 0] pAr15pAi3pAi8mAi12mAi6mAi10mAr1mAi13;  // signal num: 534
wire signed [ 1 * w - 1 + 3 : 0] pAr15pAi3pAi8mAi12pAr8mAr12pAr3mAi15;  // signal num: 535
wire signed [ 1 * w - 1 + 2 : 0] pAr15pAi3pAr10mAr6;  // signal num: 536
wire signed [ 1 * w - 1 + 1 : 0] pAr16mAi4;  // signal num: 537
wire signed [ 1 * w - 1 + 2 : 0] pAr16mAi4mAr14mAi2;  // signal num: 538
wire signed [ 1 * w - 1 + 2 : 0] pAr16mAi4pAr4pAi16;  // signal num: 539
wire signed [ 1 * w - 1 + 1 : 0] pAr16pAi4;  // signal num: 540
wire signed [ 1 * w - 1 + 2 : 0] pAr16pAi4mAr4pAi16;  // signal num: 541
wire signed [ 1 * w - 1 + 2 : 0] pAr16pAi4pAr14mAi2;  // signal num: 542
wire signed [ 1 * w - 1 + 2 : 0] pAr16pAi4pAr4mAi16;  // signal num: 543
wire signed [ 1 * w - 1 + 1 : 0] pAr1mAi13;  // signal num: 544
wire signed [ 1 * w - 1 + 1 : 0] pAr1pAi13;  // signal num: 545
wire signed [ 1 * w - 1 + 1 : 0] pAr1pAr10;  // signal num: 546
wire signed [ 1 * w - 1 + 2 : 0] pAr1pAr10pAr6mAi13;  // signal num: 547
wire signed [ 1 * w - 1 + 3 : 0] pAr1pAr10pAr6mAi13pAi6pAi10pAr13pAi1;  // signal num: 548
wire signed [ 1 * w - 1 + 1 : 0] pAr2mAi14;  // signal num: 549
wire signed [ 1 * w - 1 + 2 : 0] pAr2mAi14pAr4pAi16;  // signal num: 550
wire signed [ 1 * w - 1 + 1 : 0] pAr2pAi14;  // signal num: 551
wire signed [ 1 * w - 1 + 2 : 0] pAr2pAi14mAr14pAi2;  // signal num: 552
wire signed [ 1 * w - 1 + 3 : 0] pAr2pAi14mAr14pAi2mAr15pAi3pAr3pAi15;  // signal num: 553
wire signed [ 1 * w - 1 + 2 : 0] pAr2pAi14mAr4pAi16;  // signal num: 554
wire signed [ 1 * w - 1 + 3 : 0] pAr2pAi14mAr4pAi16pAi7pAi11pAr5mAr9;  // signal num: 555
wire signed [ 1 * w - 1 + 3 : 0] pAr2pAi14mAr4pAi16pAr11mAr7mAi5mAi9;  // signal num: 556
wire signed [ 1 * w - 1 + 2 : 0] pAr2pAi14pAr14mAi2;  // signal num: 557
wire signed [ 1 * w - 1 + 1 : 0] pAr3mAi15;  // signal num: 558
wire signed [ 1 * w - 1 + 2 : 0] pAr3mAi15mAi8mAi12;  // signal num: 559
wire signed [ 1 * w - 1 + 2 : 0] pAr3mAi15mAr15mAi3;  // signal num: 560
wire signed [ 1 * w - 1 + 1 : 0] pAr3pAi15;  // signal num: 561
wire signed [ 1 * w - 1 + 1 : 0] pAr4mAi16;  // signal num: 562
wire signed [ 1 * w - 1 + 2 : 0] pAr4mAi16mAr16mAi4;  // signal num: 563
wire signed [ 1 * w - 1 + 2 : 0] pAr4mAi16mAr2mAi14;  // signal num: 564
wire signed [ 1 * w - 1 + 1 : 0] pAr4pAi16;  // signal num: 565
wire signed [ 1 * w - 1 + 1 : 0] pAr5mAr9;  // signal num: 566
wire signed [ 1 * w - 1 + 2 : 0] pAr5mAr9pAi5mAi9;  // signal num: 567
wire signed [ 1 * w - 1 + 2 : 0] pAr5mAr9pAi9mAi5;  // signal num: 568
wire signed [ 1 * w - 1 + 3 : 0] pAr5mAr9pAi9mAi5mAr14pAi2mAr2mAi14;  // signal num: 569
wire signed [ 1 * w - 1 + 1 : 0] pAr5pAr9;  // signal num: 570
wire signed [ 1 * w - 1 + 1 : 0] pAr6mAi13;  // signal num: 571
wire signed [ 1 * w - 1 + 1 : 0] pAr6mAr10;  // signal num: 572
wire signed [ 1 * w - 1 + 3 : 0] pAr6pAi13pAr1pAr10pAr13mAi1pAi6pAi10;  // signal num: 573
wire signed [ 1 * w - 1 + 1 : 0] pAr6pAr10;  // signal num: 574
wire signed [ 1 * w - 1 + 3 : 0] pAr6pAr10mAr3pAi15pAr12mAr8;  // signal num: 575
wire signed [ 1 * w - 1 + 1 : 0] pAr7mAr11;  // signal num: 576
wire signed [ 1 * w - 1 + 1 : 0] pAr7pAr11;  // signal num: 577
wire signed [ 1 * w - 1 + 1 : 0] pAr8mAr12;  // signal num: 578
wire signed [ 1 * w - 1 + 2 : 0] pAr8mAr12pAr1mAi13;  // signal num: 579
wire signed [ 1 * w - 1 + 3 : 0] pAr8mAr12pAr1pAr10pAr6mAi13;  // signal num: 580
wire signed [ 1 * w - 1 + 2 : 0] pAr8mAr12pAr3mAi15;  // signal num: 581
wire signed [ 1 * w - 1 + 1 : 0] pAr8pAr12;  // signal num: 582
wire signed [ 1 * w - 1 + 2 : 0] pAr8pAr12mAi8mAi12;  // signal num: 583
wire signed [ 1 * w - 1 + 2 : 0] pAr8pAr12pAi8pAi12;  // signal num: 584
wire signed [ 1 * w - 1 + 2 : 0] pAr8pAr12pAr13pAi1;  // signal num: 585
wire signed [ 1 * w - 1 + 3 : 0] pAr8pAr12pAr13pAi1mAi8mAi12pAr1mAi13;  // signal num: 586
wire signed [ 1 * w - 1 + 3 : 0] pAr8pAr12pAr13pAi1mAr3mAi15pAr6mAr10;  // signal num: 587
wire signed [ 1 * w - 1 + 3 : 0] pAr8pAr12pAr13pAi1pAi8pAi12mAr1pAi13;  // signal num: 588
wire signed [ 1 * w - 1 + 2 : 0] pAr8pAr12pAr15pAi3;  // signal num: 589
wire signed [ 1 * w - 1 + 3 : 0] pAr8pAr12pAr15pAi3mAr6pAi13pAr1pAr10;  // signal num: 590
wire signed [ 1 * w - 1 + 1 : 0] pAr9mAr5;  // signal num: 591
wire signed [ 1 * w - 1 + 2 : 0] pAr9mAr5mAi7mAi11;  // signal num: 592
wire signed [ 1 * w - 1 + 2 : 0] pBi10pBr15pBi11;  // signal num: 593
wire signed [ 1 * w - 1 + 1 : 0] pBi12mBi9;  // signal num: 594
wire signed [ 1 * w - 1 + 2 : 0] pBi12pBi10pBr15pBi11;  // signal num: 595
wire signed [ 1 * w - 1 + 2 : 0] pBi13pBr11mBi15;  // signal num: 596
wire signed [ 1 * w - 1 + 1 : 0] pBi16mBi13;  // signal num: 597
wire signed [ 1 * w - 1 + 3 : 0] pBi16mBi13pBr9mBr11pBi15;  // signal num: 598
wire signed [ 1 * w - 1 + 2 : 0] pBi16mBr2pBi14;  // signal num: 599
wire signed [ 1 * w - 1 + 3 : 0] pBi16mBr2pBi14mBr11pBi15mBr10mBr12;  // signal num: 600
wire signed [ 1 * w - 1 + 2 : 0] pBi16pBr2pBi14;  // signal num: 601
wire signed [ 1 * w - 1 + 2 : 0] pBi1mBi6pBr16mBr13;  // signal num: 602
wire signed [ 1 * w - 1 + 1 : 0] pBi1pBi5;  // signal num: 603
wire signed [ 1 * w - 1 + 1 : 0] pBi2mBi6;  // signal num: 604
wire signed [ 1 * w - 1 + 1 : 0] pBi2pBi6;  // signal num: 605
wire signed [ 1 * w - 1 + 1 : 0] pBi3mBi7;  // signal num: 606
wire signed [ 1 * w - 1 + 2 : 0] pBi3mBi7pBi2mBi6;  // signal num: 607
wire signed [ 1 * w - 1 + 1 : 0] pBi3pBi7;  // signal num: 608
wire signed [ 1 * w - 1 + 1 : 0] pBi4mBi8;  // signal num: 609
wire signed [ 1 * w - 1 + 1 : 0] pBi4pBi8;  // signal num: 610
wire signed [ 1 * w - 1 + 2 : 0] pBi4pBi8mBi1mBi5;  // signal num: 611
wire signed [ 1 * w - 1 + 2 : 0] pBi4pBi8pBi2pBi6;  // signal num: 612
wire signed [ 1 * w - 1 + 1 : 0] pBi5mBi1;  // signal num: 613
wire signed [ 1 * w - 1 + 2 : 0] pBi5mBi1pBi3mBi7;  // signal num: 614
wire signed [ 1 * w - 1 + 1 : 0] pBi6mBi2;  // signal num: 615
wire signed [ 1 * w - 1 + 2 : 0] pBi6mBi2mBr14pBi10;  // signal num: 616
wire signed [ 1 * w - 1 + 2 : 0] pBi6mBi2pBi7mBi3;  // signal num: 617
wire signed [ 1 * w - 1 + 2 : 0] pBi6pBr13pBi1;  // signal num: 618
wire signed [ 1 * w - 1 + 1 : 0] pBi7mBi3;  // signal num: 619
wire signed [ 1 * w - 1 + 1 : 0] pBi9mBi5;  // signal num: 620
wire signed [ 1 * w - 1 + 2 : 0] pBi9mBr15mBi11;  // signal num: 621
wire signed [ 1 * w - 1 + 1 : 0] pBr10mBr6;  // signal num: 622
wire signed [ 1 * w - 1 + 2 : 0] pBr10pBr11mBi15;  // signal num: 623
wire signed [ 1 * w - 1 + 1 : 0] pBr11mBi15;  // signal num: 624
wire signed [ 1 * w - 1 + 1 : 0] pBr12mBr9;  // signal num: 625
wire signed [ 1 * w - 1 + 3 : 0] pBr12mBr9mBi16pBi13pBr11mBi15;  // signal num: 626
wire signed [ 1 * w - 1 + 2 : 0] pBr12mBr9pBi16mBi13;  // signal num: 627
wire signed [ 1 * w - 1 + 2 : 0] pBr12pBr10mBr6;  // signal num: 628
wire signed [ 1 * w - 1 + 3 : 0] pBr12pBr4pBr8pBr6pBr10;  // signal num: 629
wire signed [ 1 * w - 1 + 1 : 0] pBr13mBr16;  // signal num: 630
wire signed [ 1 * w - 1 + 3 : 0] pBr13mBr16mBi12pBi9mBr15mBi11;  // signal num: 631
wire signed [ 1 * w - 1 + 2 : 0] pBr13mBr16pBi12mBi9;  // signal num: 632
wire signed [ 1 * w - 1 + 3 : 0] pBr13mBr16pBi12mBi9pBi4pBi8mBi1mBi5;  // signal num: 633
wire signed [ 1 * w - 1 + 1 : 0] pBr13pBi1;  // signal num: 634
wire signed [ 1 * w - 1 + 2 : 0] pBr13pBi9mBr15mBi11;  // signal num: 635
wire signed [ 1 * w - 1 + 2 : 0] pBr14pBi6mBi2;  // signal num: 636
wire signed [ 1 * w - 1 + 1 : 0] pBr14pBr16;  // signal num: 637
wire signed [ 1 * w - 1 + 2 : 0] pBr14pBr16mBi10mBi12;  // signal num: 638
wire signed [ 1 * w - 1 + 3 : 0] pBr14pBr16pBi12pBi10pBr15pBi11;  // signal num: 639
wire signed [ 1 * w - 1 + 1 : 0] pBr15pBi11;  // signal num: 640
wire signed [ 1 * w - 1 + 1 : 0] pBr16mBr13;  // signal num: 641
wire signed [ 1 * w - 1 + 1 : 0] pBr1mBr5;  // signal num: 642
wire signed [ 1 * w - 1 + 2 : 0] pBr1pBr5mBr4mBr8;  // signal num: 643
wire signed [ 1 * w - 1 + 2 : 0] pBr2mBi14mBi16;  // signal num: 644
wire signed [ 1 * w - 1 + 1 : 0] pBr2pBi14;  // signal num: 645
wire signed [ 1 * w - 1 + 2 : 0] pBr2pBi14pBr10mBr6;  // signal num: 646
wire signed [ 1 * w - 1 + 1 : 0] pBr3mBr7;  // signal num: 647
wire signed [ 1 * w - 1 + 2 : 0] pBr3mBr7pBr10mBr6;  // signal num: 648
wire signed [ 1 * w - 1 + 2 : 0] pBr3mBr7pBr5mBr1;  // signal num: 649
wire signed [ 1 * w - 1 + 1 : 0] pBr3pBr7;  // signal num: 650
wire signed [ 1 * w - 1 + 1 : 0] pBr4mBr8;  // signal num: 651
wire signed [ 1 * w - 1 + 1 : 0] pBr4pBr8;  // signal num: 652
wire signed [ 1 * w - 1 + 2 : 0] pBr4pBr8mBr1mBr5;  // signal num: 653
wire signed [ 1 * w - 1 + 2 : 0] pBr4pBr8pBr3pBr7;  // signal num: 654
wire signed [ 1 * w - 1 + 2 : 0] pBr4pBr8pBr6pBr10;  // signal num: 655
wire signed [ 1 * w - 1 + 1 : 0] pBr5mBr1;  // signal num: 656
wire signed [ 1 * w - 1 + 2 : 0] pBr5mBr4mBr8;  // signal num: 657
wire signed [ 1 * w - 1 + 1 : 0] pBr6pBr10;  // signal num: 658
wire signed [ 1 * w - 1 + 2 : 0] pBr6pBr10mBr2pBi14;  // signal num: 659
wire signed [ 1 * w - 1 + 3 : 0] pBr6pBr4pBr8pBr3pBr7;  // signal num: 660
wire signed [ 1 * w - 1 + 2 : 0] pBr6pBr7mBr3;  // signal num: 661
wire signed [ 1 * w - 1 + 1 : 0] pBr7mBr3;  // signal num: 662
wire signed [ 1 * w - 1 + 2 : 0] pBr9mBr11pBi15;  // signal num: 663
wire signed [ 1 * w - 1 + 1 : 0] pBr9pBi13;  // signal num: 664
wire signed [ 2 * w - 1 + 6 : 0] pt00_imt14_i;  // signal num: 665
wire signed [ 2 * w - 1 + 6 : 0] pt00_rmt41_r;  // signal num: 666
wire signed [ 2 * w - 1 + 6 : 0] pt01_ipt28_r;  // signal num: 667
wire signed [ 2 * w - 1 + 7 : 0] pt01_rmt16_rpt09_imt00_rmt40_imt46_rpt14_rmt28_i;  // signal num: 668
wire signed [ 2 * w - 1 + 6 : 0] pt02_imt35_r;  // signal num: 669
wire signed [ 2 * w - 1 + 7 : 0] pt03_imt14_rmt34_ipt18_imt43_r;  // signal num: 670
wire signed [ 2 * w - 1 + 7 : 0] pt03_ipt43_imt34_rmt41_ipt00_imt14_i;  // signal num: 671
wire signed [ 2 * w - 1 + 7 : 0] pt03_rmt31_rpt11_imt08_rpt43_rmt39_i;  // signal num: 672
wire signed [ 2 * w - 1 + 7 : 0] pt05_rmt36_r;  // signal num: 673
wire signed [ 2 * w - 1 + 6 : 0] pt05_rpt06_r;  // signal num: 674
wire signed [ 2 * w - 1 + 7 : 0] pt06_imt17_rmt22_i;  // signal num: 675
wire signed [ 2 * w - 1 + 6 : 0] pt06_imt24_i;  // signal num: 676
wire signed [ 2 * w - 1 + 6 : 0] pt06_rmt05_r;  // signal num: 677
wire signed [ 2 * w - 1 + 7 : 0] pt08_rpt16_ipt27_ipt11_rmt39_r;  // signal num: 678
wire signed [ 2 * w - 1 + 7 : 0] pt09_imt00_rmt40_imt46_r;  // signal num: 679
wire signed [ 2 * w - 1 + 7 : 0] pt09_imt00_rmt40_imt46_rpt14_ipt18_i;  // signal num: 680
wire signed [ 2 * w - 1 + 7 : 0] pt09_imt00_rmt40_imt46_rpt14_rmt28_i;  // signal num: 681
wire signed [ 2 * w - 1 + 7 : 0] pt09_rmt40_r;  // signal num: 682
wire signed [ 2 * w - 1 + 7 : 0] pt11_imt08_rpt43_rmt39_imt37_imt24_imt16_imt18_i;  // signal num: 683
wire signed [ 2 * w - 1 + 7 : 0] pt12_imt15_rmt02_ipt33_ipt47_i;  // signal num: 684
wire signed [ 2 * w - 1 + 7 : 0] pt12_imt15_rpt02_imt35_r;  // signal num: 685
wire signed [ 2 * w - 1 + 7 : 0] pt12_rpt15_ipt32_rpt40_i;  // signal num: 686
wire signed [ 2 * w - 1 + 6 : 0] pt13_imt44_r;  // signal num: 687
wire signed [ 2 * w - 1 + 6 : 0] pt13_rpt44_i;  // signal num: 688
wire signed [ 2 * w - 1 + 6 : 0] pt14_ipt18_i;  // signal num: 689
wire signed [ 2 * w - 1 + 6 : 0] pt14_ipt27_r;  // signal num: 690
wire signed [ 2 * w - 1 + 7 : 0] pt14_rmt27_ipt08_rpt16_i;  // signal num: 691
wire signed [ 2 * w - 1 + 7 : 0] pt14_rmt27_ipt21_imt33_r;  // signal num: 692
wire signed [ 2 * w - 1 + 6 : 0] pt14_rmt28_i;  // signal num: 693
wire signed [ 2 * w - 1 + 6 : 0] pt15_imt30_i;  // signal num: 694
wire signed [ 2 * w - 1 + 7 : 0] pt16_ipt09_rmt40_r;  // signal num: 695
wire signed [ 2 * w - 1 + 6 : 0] pt16_rmt46_i;  // signal num: 696
wire signed [ 2 * w - 1 + 7 : 0] pt16_rmt46_imt09_rmt14_r;  // signal num: 697
wire signed [ 2 * w - 1 + 7 : 0] pt17_ipt18_rpt39_rmt38_rpt05_rmt36_r;  // signal num: 698
wire signed [ 2 * w - 1 + 7 : 0] pt17_rmt15_rpt34_imt11_ipt11_rmt39_rpt08_imt43_i;  // signal num: 699
wire signed [ 2 * w - 1 + 7 : 0] pt17_rpt22_i;  // signal num: 700
wire signed [ 2 * w - 1 + 7 : 0] pt17_rpt22_imt28_rpt42_imt01_i;  // signal num: 701
wire signed [ 2 * w - 1 + 7 : 0] pt18_imt43_rpt24_imt11_ipt37_ipt39_i;  // signal num: 702
wire signed [ 2 * w - 1 + 6 : 0] pt19_rmt23_i;  // signal num: 703
wire signed [ 2 * w - 1 + 7 : 0] pt22_rmt17_i;  // signal num: 704
wire signed [ 2 * w - 1 + 7 : 0] pt22_rmt17_ipt06_rmt05_r;  // signal num: 705
wire signed [ 2 * w - 1 + 7 : 0] pt24_imt11_ipt37_ipt39_i;  // signal num: 706
wire signed [ 2 * w - 1 + 7 : 0] pt24_imt11_ipt37_ipt39_ipt14_ipt27_rpt03_rmt31_rpt34_rmt21_i;  // signal num: 707
wire signed [ 2 * w - 1 + 6 : 0] pt26_imt32_i;  // signal num: 708
wire signed [ 2 * w - 1 + 6 : 0] pt26_rmt32_r;  // signal num: 709
wire signed [ 2 * w - 1 + 7 : 0] pt26_rmt32_rmt16_rmt18_rmt24_rmt37_r;  // signal num: 710
wire signed [ 2 * w - 1 + 7 : 0] pt27_imt33_r;  // signal num: 711
wire signed [ 2 * w - 1 + 6 : 0] pt27_rmt18_i;  // signal num: 712
wire signed [ 2 * w - 1 + 7 : 0] pt27_rmt18_ipt33_ipt37_i;  // signal num: 713
wire signed [ 2 * w - 1 + 6 : 0] pt28_imt01_r;  // signal num: 714
wire signed [ 2 * w - 1 + 7 : 0] pt28_imt01_rpt11_rmt39_rpt08_imt43_i;  // signal num: 715
wire signed [ 2 * w - 1 + 6 : 0] pt29_ipt44_r;  // signal num: 716
wire signed [ 2 * w - 1 + 7 : 0] pt29_ipt44_rmt03_rmt13_imt20_rmt33_imt37_i;  // signal num: 717
wire signed [ 2 * w - 1 + 6 : 0] pt29_rmt45_rpt28_imt01_rpt09_rmt40_rmt06_imt10_i;  // signal num: 718
wire signed [ 2 * w - 1 + 7 : 0] pt30_ipt36_imt05_i;  // signal num: 719
wire signed [ 2 * w - 1 + 7 : 0] pt30_ipt36_imt05_imt28_rpt42_imt01_i;  // signal num: 720
wire signed [ 2 * w - 1 + 7 : 0] pt30_ipt36_imt05_ipt23_rpt13_imt44_r;  // signal num: 721
wire signed [ 2 * w - 1 + 6 : 0] pt30_rmt38_i;  // signal num: 722
wire signed [ 2 * w - 1 + 7 : 0] pt31_imt03_ipt22_rmt17_i;  // signal num: 723
wire signed [ 2 * w - 1 + 6 : 0] pt32_imt26_i;  // signal num: 724
wire signed [ 2 * w - 1 + 6 : 0] pt32_rpt40_i;  // signal num: 725
wire signed [ 2 * w - 1 + 6 : 0] pt33_ipt47_i;  // signal num: 726
wire signed [ 2 * w - 1 + 7 : 0] pt33_rpt37_r;  // signal num: 727
wire signed [ 2 * w - 1 + 6 : 0] pt34_imt11_i;  // signal num: 728
wire signed [ 2 * w - 1 + 7 : 0] pt34_imt11_ipt11_rmt39_rpt08_imt43_i;  // signal num: 729
wire signed [ 2 * w - 1 + 6 : 0] pt35_imt19_i;  // signal num: 730
wire signed [ 2 * w - 1 + 7 : 0] pt36_imt05_i;  // signal num: 731
wire signed [ 2 * w - 1 + 7 : 0] pt36_imt05_ipt30_rmt38_i;  // signal num: 732
wire signed [ 2 * w - 1 + 7 : 0] pt36_rmt05_r;  // signal num: 733
wire signed [ 2 * w - 1 + 7 : 0] pt37_ipt39_i;  // signal num: 734
wire signed [ 2 * w - 1 + 6 : 0] pt38_imt47_i;  // signal num: 735
wire signed [ 2 * w - 1 + 7 : 0] pt38_imt47_imt16_rmt18_rmt24_rmt37_rpt32_imt26_i;  // signal num: 736
wire signed [ 2 * w - 1 + 7 : 0] pt39_ipt47_ipt26_rmt32_rmt16_rmt18_rmt24_rmt37_r;  // signal num: 737
wire signed [ 2 * w - 1 + 7 : 0] pt39_ipt47_ipt26_rmt32_rmt16_rmt18_rmt24_rmt37_rpt30_rmt15_i;  // signal num: 738
wire signed [ 2 * w - 1 + 7 : 0] pt39_rmt38_rpt05_rmt36_r;  // signal num: 739
wire signed [ 2 * w - 1 + 6 : 0] pt41_imt00_i;  // signal num: 740
wire signed [ 2 * w - 1 + 6 : 0] pt41_rmt00_r;  // signal num: 741
wire signed [ 2 * w - 1 + 6 : 0] pt42_imt01_i;  // signal num: 742
wire signed [ 2 * w - 1 + 6 : 0] pt42_ipt46_r;  // signal num: 743
wire signed [ 2 * w - 1 + 7 : 0] pt42_ipt46_rpt08_rpt16_ipt27_ipt11_rmt39_r;  // signal num: 744
wire signed [ 2 * w - 1 + 7 : 0] pt42_rmt16_ipt28_imt01_rpt11_rmt39_rpt08_imt43_i;  // signal num: 745
wire signed [ 2 * w - 1 + 7 : 0] pt43_imt34_rpt21_imt33_r;  // signal num: 746
wire signed [ 2 * w - 1 + 7 : 0] pt43_rmt08_rmt16_imt18_ipt21_rpt34_i;  // signal num: 747
wire signed [ 2 * w - 1 + 7 : 0] pt45_imt29_ipt04_ipt46_i;  // signal num: 748
wire signed [ 2 * w - 1 + 6 : 0] pt47_rmt30_r;  // signal num: 749
wire signed [ 2 * w - 1 + 7 : 0] pvery_long_op_0001;  // signal num: 750
wire signed [ 2 * w - 1 + 7 : 0] pvery_long_op_0002;  // signal num: 751
wire signed [ 2 * w - 1 + 7 : 0] pvery_long_op_0003;  // signal num: 752
wire signed [ 2 * w - 1 + 7 : 0] pvery_long_op_0004;  // signal num: 753
wire signed [ 2 * w - 1 + 7 : 0] pvery_long_op_0005;  // signal num: 754
wire signed [ 2 * w - 1 + 7 : 0] pvery_long_op_0006;  // signal num: 755
wire signed [ 2 * w - 1 + 5 : 0] q00;  // signal num: 756
wire signed [ 2 * w - 1 + 6 : 0] q01;  // signal num: 757
wire signed [ 2 * w - 1 + 5 : 0] q02;  // signal num: 758
wire signed [ 2 * w - 1 + 6 : 0] q03;  // signal num: 759
wire signed [ 2 * w - 1 + 6 : 0] q04;  // signal num: 760
wire signed [ 2 * w - 1 + 6 : 0] q05;  // signal num: 761
wire signed [ 2 * w - 1 + 5 : 0] q06;  // signal num: 762
wire signed [ 2 * w - 1 + 6 : 0] q07;  // signal num: 763
wire signed [ 2 * w - 1 + 6 : 0] q08;  // signal num: 764
wire signed [ 2 * w - 1 + 6 : 0] q10;  // signal num: 765
wire signed [ 2 * w - 1 + 6 : 0] q11;  // signal num: 766
wire signed [ 2 * w - 1 + 6 : 0] q13;  // signal num: 767
wire signed [ 2 * w - 1 + 5 : 0] q14;  // signal num: 768
wire signed [ 2 * w - 1 + 6 : 0] q15;  // signal num: 769
wire signed [ 2 * w - 1 + 6 : 0] q16;  // signal num: 770
wire signed [ 2 * w - 1 + 6 : 0] q17;  // signal num: 771
wire signed [ 2 * w - 1 + 6 : 0] q18;  // signal num: 772
wire signed [ 2 * w - 1 + 6 : 0] q19;  // signal num: 773
wire signed [ 2 * w - 1 + 6 : 0] q20;  // signal num: 774
wire signed [ 2 * w - 1 + 6 : 0] q21;  // signal num: 775
wire signed [ 2 * w - 1 + 6 : 0] q23;  // signal num: 776
wire signed [ 2 * w - 1 + 6 : 0] q24;  // signal num: 777
wire signed [ 2 * w - 1 + 6 : 0] q26;  // signal num: 778
wire signed [ 2 * w - 1 + 6 : 0] q27;  // signal num: 779
wire signed [ 2 * w - 1 + 5 : 0] q28;  // signal num: 780
wire signed [ 2 * w - 1 + 6 : 0] q29;  // signal num: 781
wire signed [ 2 * w - 1 + 6 : 0] q30;  // signal num: 782
wire signed [ 2 * w - 1 + 6 : 0] q31;  // signal num: 783
wire signed [ 2 * w - 1 + 5 : 0] q32;  // signal num: 784
wire signed [ 2 * w - 1 + 6 : 0] q33;  // signal num: 785
wire signed [ 2 * w - 1 + 6 : 0] q34;  // signal num: 786
wire signed [ 2 * w - 1 + 6 : 0] q35;  // signal num: 787
wire signed [ 2 * w - 1 + 6 : 0] q37;  // signal num: 788
wire signed [ 2 * w - 1 + 5 : 0] q38;  // signal num: 789
wire signed [ 2 * w - 1 + 6 : 0] q39;  // signal num: 790
wire signed [ 2 * w - 1 + 6 : 0] q40;  // signal num: 791
wire signed [ 2 * w - 1 + 6 : 0] q41;  // signal num: 792
wire signed [ 2 * w - 1 + 5 : 0] q42;  // signal num: 793
wire signed [ 2 * w - 1 + 6 : 0] q43;  // signal num: 794
wire signed [ 2 * w - 1 + 5 : 0] q44;  // signal num: 795
wire signed [ 2 * w - 1 + 6 : 0] q46;  // signal num: 796
wire signed [ 2 * w - 1 + 5 : 0] q47;  // signal num: 797
wire signed [ 2 * w - 1 + 5 : 0] t00_i;  // signal num: 798
wire signed [ 2 * w - 1 + 5 : 0] t00_r;  // signal num: 799
wire signed [ 2 * w - 1 + 6 : 0] t01_i;  // signal num: 800
wire signed [ 2 * w - 1 + 6 : 0] t01_r;  // signal num: 801
wire signed [ 2 * w - 1 + 5 : 0] t02_i;  // signal num: 802
wire signed [ 2 * w - 1 + 5 : 0] t02_r;  // signal num: 803
wire signed [ 2 * w - 1 + 6 : 0] t03_i;  // signal num: 804
wire signed [ 2 * w - 1 + 6 : 0] t03_r;  // signal num: 805
wire signed [ 2 * w - 1 + 6 : 0] t05_i;  // signal num: 806
wire signed [ 2 * w - 1 + 6 : 0] t05_r;  // signal num: 807
wire signed [ 2 * w - 1 + 5 : 0] t06_i;  // signal num: 808
wire signed [ 2 * w - 1 + 5 : 0] t06_r;  // signal num: 809
wire signed [ 2 * w - 1 + 6 : 0] t07_i;  // signal num: 810
wire signed [ 2 * w - 1 + 6 : 0] t07_r;  // signal num: 811
wire signed [ 2 * w - 1 + 6 : 0] t08_i;  // signal num: 812
wire signed [ 2 * w - 1 + 6 : 0] t10_i;  // signal num: 813
wire signed [ 2 * w - 1 + 6 : 0] t10_r;  // signal num: 814
wire signed [ 2 * w - 1 + 6 : 0] t11_i;  // signal num: 815
wire signed [ 2 * w - 1 + 6 : 0] t11_r;  // signal num: 816
wire signed [ 2 * w - 1 + 6 : 0] t13_i;  // signal num: 817
wire signed [ 2 * w - 1 + 6 : 0] t13_r;  // signal num: 818
wire signed [ 2 * w - 1 + 6 : 0] t15_i;  // signal num: 819
wire signed [ 2 * w - 1 + 6 : 0] t15_r;  // signal num: 820
wire signed [ 2 * w - 1 + 6 : 0] t16_i;  // signal num: 821
wire signed [ 2 * w - 1 + 6 : 0] t16_r;  // signal num: 822
wire signed [ 2 * w - 1 + 6 : 0] t17_i;  // signal num: 823
wire signed [ 2 * w - 1 + 6 : 0] t17_r;  // signal num: 824
wire signed [ 2 * w - 1 + 6 : 0] t18_i;  // signal num: 825
wire signed [ 2 * w - 1 + 6 : 0] t18_r;  // signal num: 826
wire signed [ 2 * w - 1 + 6 : 0] t19_i;  // signal num: 827
wire signed [ 2 * w - 1 + 6 : 0] t19_r;  // signal num: 828
wire signed [ 2 * w - 1 + 6 : 0] t20_i;  // signal num: 829
wire signed [ 2 * w - 1 + 6 : 0] t20_r;  // signal num: 830
wire signed [ 2 * w - 1 + 6 : 0] t21_r;  // signal num: 831
wire signed [ 2 * w - 1 + 6 : 0] t24_r;  // signal num: 832
wire signed [ 2 * w - 1 + 6 : 0] t26_i;  // signal num: 833
wire signed [ 2 * w - 1 + 6 : 0] t26_r;  // signal num: 834
wire signed [ 2 * w - 1 + 6 : 0] t27_i;  // signal num: 835
wire signed [ 2 * w - 1 + 6 : 0] t27_r;  // signal num: 836
wire signed [ 2 * w - 1 + 5 : 0] t28_i;  // signal num: 837
wire signed [ 2 * w - 1 + 5 : 0] t28_r;  // signal num: 838
wire signed [ 2 * w - 1 + 6 : 0] t29_i;  // signal num: 839
wire signed [ 2 * w - 1 + 6 : 0] t29_r;  // signal num: 840
wire signed [ 2 * w - 1 + 6 : 0] t30_i;  // signal num: 841
wire signed [ 2 * w - 1 + 6 : 0] t30_r;  // signal num: 842
wire signed [ 2 * w - 1 + 6 : 0] t31_i;  // signal num: 843
wire signed [ 2 * w - 1 + 6 : 0] t31_r;  // signal num: 844
wire signed [ 2 * w - 1 + 5 : 0] t32_i;  // signal num: 845
wire signed [ 2 * w - 1 + 5 : 0] t32_r;  // signal num: 846
wire signed [ 2 * w - 1 + 6 : 0] t33_r;  // signal num: 847
wire signed [ 2 * w - 1 + 6 : 0] t34_r;  // signal num: 848
wire signed [ 2 * w - 1 + 6 : 0] t35_i;  // signal num: 849
wire signed [ 2 * w - 1 + 6 : 0] t35_r;  // signal num: 850
wire signed [ 2 * w - 1 + 7 : 0] t36_r;  // signal num: 851
wire signed [ 2 * w - 1 + 7 : 0] t37_r;  // signal num: 852
wire signed [ 2 * w - 1 + 5 : 0] t38_i;  // signal num: 853
wire signed [ 2 * w - 1 + 5 : 0] t38_r;  // signal num: 854
wire signed [ 2 * w - 1 + 6 : 0] t39_r;  // signal num: 855
wire signed [ 2 * w - 1 + 6 : 0] t40_i;  // signal num: 856
wire signed [ 2 * w - 1 + 6 : 0] t40_r;  // signal num: 857
wire signed [ 2 * w - 1 + 6 : 0] t41_i;  // signal num: 858
wire signed [ 2 * w - 1 + 6 : 0] t41_r;  // signal num: 859
wire signed [ 2 * w - 1 + 5 : 0] t42_i;  // signal num: 860
wire signed [ 2 * w - 1 + 5 : 0] t42_r;  // signal num: 861
wire signed [ 2 * w - 1 + 5 : 0] t44_i;  // signal num: 862
wire signed [ 2 * w - 1 + 5 : 0] t44_r;  // signal num: 863
wire signed [ 2 * w - 1 + 6 : 0] t46_i;  // signal num: 864
wire signed [ 2 * w - 1 + 6 : 0] t46_r;  // signal num: 865
wire signed [ 2 * w - 1 + 6 : 0] t47_i;  // signal num: 866
wire signed [ 2 * w - 1 + 6 : 0] t47_r;  // signal num: 867
wire signed [ 2 * w - 1 + 8 : 0] k24;  // signal num: 868
wire signed [ 2 * w - 1 + 7 : 0] k25;  // signal num: 869
wire signed [ 2 * w - 1 + 8 : 0] k33;  // signal num: 870
wire signed [ 2 * w - 1 + 8 : 0] k34;  // signal num: 871
wire signed [ 2 * w - 1 + 8 : 0] k37;  // signal num: 872
wire signed [ 2 * w - 1 + 7 : 0] k43;  // signal num: 873
wire signed [ 2 * w - 1 + 7 : 0] k47;  // signal num: 874
wire signed [ 2 * w - 1 + 7 : 0] mt02_ipt33_ipt47_i;  // signal num: 875
wire signed [ 2 * w - 1 + 7 : 0] mt02_rpt33_rpt47_r;  // signal num: 876
wire signed [ 2 * w - 1 + 7 : 0] mt03_rmt13_i;  // signal num: 877
wire signed [ 2 * w - 1 + 7 : 0] mt06_rmt10_rpt07_rmt42_r;  // signal num: 878
wire signed [ 2 * w - 1 + 7 : 0] mt07_imt21_r;  // signal num: 879
wire signed [ 2 * w - 1 + 7 : 0] mt07_rmt06_rmt10_r;  // signal num: 880
wire signed [ 2 * w - 1 + 7 : 0] mt08_ipt24_rmt11_r;  // signal num: 881
wire signed [ 2 * w - 1 + 7 : 0] mt09_rmt14_r;  // signal num: 882
wire signed [ 2 * w - 1 + 7 : 0] mt12_rpt11_rmt15_i;  // signal num: 883
wire signed [ 2 * w - 1 + 7 : 0] mt16_imt18_i;  // signal num: 884
wire signed [ 2 * w - 1 + 7 : 0] mt16_rmt18_r;  // signal num: 885
wire signed [ 2 * w - 1 + 7 : 0] mt16_rmt18_rmt24_rmt37_r;  // signal num: 886
wire signed [ 2 * w - 1 + 7 : 0] mt16_rmt18_rpt08_imt43_i;  // signal num: 887
wire signed [ 2 * w - 1 + 7 : 0] mt18_ipt11_imt08_r;  // signal num: 888
wire signed [ 2 * w - 1 + 7 : 0] mt18_rmt27_ipt33_rpt37_r;  // signal num: 889
wire signed [ 2 * w - 1 + 7 : 0] mt20_ipt00_rmt41_r;  // signal num: 890
wire signed [ 2 * w - 1 + 7 : 0] mt20_rmt33_imt37_i;  // signal num: 891
wire signed [ 2 * w - 1 + 7 : 0] mt21_rpt14_ipt27_r;  // signal num: 892
wire signed [ 2 * w - 1 + 7 : 0] mt23_rmt30_i;  // signal num: 893
wire signed [ 2 * w - 1 + 7 : 0] mt23_rmt30_ipt12_rpt15_ipt13_imt44_r;  // signal num: 894
wire signed [ 2 * w - 1 + 7 : 0] mt24_imt16_imt18_i;  // signal num: 895
wire signed [ 2 * w - 1 + 7 : 0] mt24_imt16_imt18_imt21_rpt14_ipt27_r;  // signal num: 896
wire signed [ 2 * w - 1 + 7 : 0] mt24_rpt41_imt00_i;  // signal num: 897
wire signed [ 2 * w - 1 + 7 : 0] mt28_imt20_ipt00_rmt41_r;  // signal num: 898
wire signed [ 2 * w - 1 + 7 : 0] mt28_rpt42_imt01_i;  // signal num: 899
wire signed [ 2 * w - 1 + 7 : 0] mt29_rmt03_rpt13_rpt44_i;  // signal num: 900
wire signed [ 2 * w - 1 + 7 : 0] mt33_imt37_i;  // signal num: 901
wire signed [ 2 * w - 1 + 7 : 0] mt34_ipt18_imt43_r;  // signal num: 902
wire signed [ 2 * w - 1 + 7 : 0] mt38_imt39_r;  // signal num: 903
wire signed [ 2 * w - 1 + 7 : 0] mt40_imt46_r;  // signal num: 904
wire signed [ 2 * w - 1 + 7 : 0] mt43_ipt34_rmt11_r;  // signal num: 905
wire signed [ 2 * w - 1 + 7 : 0] mt44_ipt28_imt01_rpt09_rmt40_r;  // signal num: 906
wire signed [ 2 * w - 1 + 7 : 0] mt44_rpt26_imt32_i;  // signal num: 907
wire signed [ 2 * w - 1 + 7 : 0] p04;  // signal num: 908
wire signed [ 2 * w - 1 + 7 : 0] p12;  // signal num: 909
wire signed [ 2 * w - 1 + 7 : 0] p22;  // signal num: 910
wire signed [ 2 * w - 1 + 7 : 0] p23;  // signal num: 911
wire signed [ 2 * w - 1 + 7 : 0] p25;  // signal num: 912
wire signed [ 2 * w - 1 + 7 : 0] p36;  // signal num: 913
wire signed [ 2 * w - 1 + 7 : 0] p37;  // signal num: 914
wire signed [ 2 * w - 1 + 7 : 0] p45;  // signal num: 915
wire signed [ 2 * w - 1 + 6 : 0] pt00_ipt42_r;  // signal num: 916
wire signed [ 2 * w - 1 + 7 : 0] pt00_ipt42_rpt16_rmt46_ipt01_ipt28_rpt40_imt09_i;  // signal num: 917
wire signed [ 2 * w - 1 + 7 : 0] pt01_ipt28_rpt40_imt09_i;  // signal num: 918
wire signed [ 2 * w - 1 + 7 : 0] pt01_ipt28_rpt40_imt09_ipt14_rmt27_ipt08_rpt16_i;  // signal num: 919
wire signed [ 2 * w - 1 + 7 : 0] pt01_rmt16_r;  // signal num: 920
wire signed [ 2 * w - 1 + 7 : 0] pt02_rpt33_ipt47_i;  // signal num: 921
wire signed [ 2 * w - 1 + 7 : 0] pt03_imt14_r;  // signal num: 922
wire signed [ 2 * w - 1 + 7 : 0] pt03_ipt43_imt34_r;  // signal num: 923
wire signed [ 2 * w - 1 + 7 : 0] pt03_rmt31_r;  // signal num: 924
wire signed [ 2 * w - 1 + 7 : 0] pt03_rpt34_imt11_i;  // signal num: 925
wire signed [ 2 * w - 1 + 7 : 0] pt04_ipt46_i;  // signal num: 926
wire signed [ 2 * w - 1 + 6 : 0] pt05_ipt18_r;  // signal num: 927
wire signed [ 2 * w - 1 + 7 : 0] pt07_rmt42_r;  // signal num: 928
wire signed [ 2 * w - 1 + 7 : 0] pt08_imt12_rpt11_rmt15_i;  // signal num: 929
wire signed [ 2 * w - 1 + 7 : 0] pt08_imt43_i;  // signal num: 930
wire signed [ 2 * w - 1 + 7 : 0] pt08_rpt16_i;  // signal num: 931
wire signed [ 2 * w - 1 + 7 : 0] pt08_rpt18_i;  // signal num: 932
wire signed [ 2 * w - 1 + 7 : 0] pt09_imt00_r;  // signal num: 933
wire signed [ 2 * w - 1 + 7 : 0] pt11_imt08_r;  // signal num: 934
wire signed [ 2 * w - 1 + 7 : 0] pt11_imt08_rpt43_rmt39_i;  // signal num: 935
wire signed [ 2 * w - 1 + 7 : 0] pt11_imt27_r;  // signal num: 936
wire signed [ 2 * w - 1 + 7 : 0] pt11_rmt15_i;  // signal num: 937
wire signed [ 2 * w - 1 + 7 : 0] pt11_rmt39_r;  // signal num: 938
wire signed [ 2 * w - 1 + 7 : 0] pt11_rmt39_rpt08_imt43_i;  // signal num: 939
wire signed [ 2 * w - 1 + 7 : 0] pt12_imt15_r;  // signal num: 940
wire signed [ 2 * w - 1 + 7 : 0] pt12_rpt15_i;  // signal num: 941
wire signed [ 2 * w - 1 + 7 : 0] pt12_rpt15_ipt13_imt44_r;  // signal num: 942
wire signed [ 2 * w - 1 + 7 : 0] pt14_ipt27_rmt16_imt18_i;  // signal num: 943
wire signed [ 2 * w - 1 + 7 : 0] pt14_ipt27_rpt03_rmt31_r;  // signal num: 944
wire signed [ 2 * w - 1 + 7 : 0] pt14_rmt27_i;  // signal num: 945
wire signed [ 2 * w - 1 + 7 : 0] pt15_rpt40_r;  // signal num: 946
wire signed [ 2 * w - 1 + 7 : 0] pt15_rpt40_rpt11_imt27_r;  // signal num: 947
wire signed [ 2 * w - 1 + 7 : 0] pt16_rmt26_r;  // signal num: 948
wire signed [ 2 * w - 1 + 7 : 0] pt16_rmt46_ipt01_ipt28_rpt40_imt09_i;  // signal num: 949
wire signed [ 2 * w - 1 + 7 : 0] pt17_ipt18_r;  // signal num: 950
wire signed [ 2 * w - 1 + 7 : 0] pt17_rmt15_r;  // signal num: 951
wire signed [ 2 * w - 1 + 7 : 0] pt18_imt43_r;  // signal num: 952
wire signed [ 2 * w - 1 + 7 : 0] pt20_ipt25_i;  // signal num: 953
wire signed [ 2 * w - 1 + 7 : 0] pt20_rpt25_r;  // signal num: 954
wire signed [ 2 * w - 1 + 7 : 0] pt20_rpt29_ipt44_rmt03_rmt13_i;  // signal num: 955
wire signed [ 2 * w - 1 + 7 : 0] pt21_imt33_r;  // signal num: 956
wire signed [ 2 * w - 1 + 7 : 0] pt21_rpt34_i;  // signal num: 957
wire signed [ 2 * w - 1 + 7 : 0] pt23_rpt13_imt44_r;  // signal num: 958
wire signed [ 2 * w - 1 + 7 : 0] pt24_imt11_i;  // signal num: 959
wire signed [ 2 * w - 1 + 7 : 0] pt24_imt11_ipt08_rpt18_i;  // signal num: 960
wire signed [ 2 * w - 1 + 7 : 0] pt24_imt11_ipt37_ipt39_ipt14_ipt27_rpt03_rmt31_r;  // signal num: 961
wire signed [ 2 * w - 1 + 7 : 0] pt24_rmt11_r;  // signal num: 962
wire signed [ 2 * w - 1 + 7 : 0] pt26_rmt32_rmt13_rmt44_i;  // signal num: 963
wire signed [ 2 * w - 1 + 7 : 0] pt27_ipt11_rmt39_r;  // signal num: 964
wire signed [ 2 * w - 1 + 7 : 0] pt28_imt01_rpt09_rmt40_r;  // signal num: 965
wire signed [ 2 * w - 1 + 7 : 0] pt28_imt01_rpt09_rmt40_rmt06_imt10_i;  // signal num: 966
wire signed [ 2 * w - 1 + 6 : 0] pt28_rpt05_ipt18_r;  // signal num: 967
wire signed [ 2 * w - 1 + 7 : 0] pt29_ipt44_rmt03_rmt13_i;  // signal num: 968
wire signed [ 2 * w - 1 + 7 : 0] pt29_rmt45_r;  // signal num: 969
wire signed [ 2 * w - 1 + 7 : 0] pt30_rmt15_i;  // signal num: 970
wire signed [ 2 * w - 1 + 7 : 0] pt31_imt03_i;  // signal num: 971
wire signed [ 2 * w - 1 + 7 : 0] pt32_rpt05_ipt18_r;  // signal num: 972
wire signed [ 2 * w - 1 + 7 : 0] pt32_rpt05_ipt18_rpt16_rmt26_r;  // signal num: 973
wire signed [ 2 * w - 1 + 7 : 0] pt33_ipt37_i;  // signal num: 974
wire signed [ 2 * w - 1 + 7 : 0] pt33_rpt47_r;  // signal num: 975
wire signed [ 2 * w - 1 + 7 : 0] pt34_rmt11_r;  // signal num: 976
wire signed [ 2 * w - 1 + 7 : 0] pt34_rmt21_i;  // signal num: 977
wire signed [ 2 * w - 1 + 7 : 0] pt37_rpt24_rmt11_r;  // signal num: 978
wire signed [ 2 * w - 1 + 7 : 0] pt39_ipt47_i;  // signal num: 979
wire signed [ 2 * w - 1 + 7 : 0] pt39_rmt38_r;  // signal num: 980
wire signed [ 2 * w - 1 + 7 : 0] pt40_imt09_i;  // signal num: 981
wire signed [ 2 * w - 1 + 7 : 0] pt41_rmt04_r;  // signal num: 982
wire signed [ 2 * w - 1 + 6 : 0] pt42_rmt16_i;  // signal num: 983
wire signed [ 2 * w - 1 + 7 : 0] pt43_imt34_r;  // signal num: 984
wire signed [ 2 * w - 1 + 7 : 0] pt43_rmt08_r;  // signal num: 985
wire signed [ 2 * w - 1 + 7 : 0] pt43_rmt39_i;  // signal num: 986
wire signed [ 2 * w - 1 + 7 : 0] pt45_imt29_i;  // signal num: 987
wire signed [ 2 * w - 1 + 7 : 0] pt46_imt44_r;  // signal num: 988
wire signed [ 2 * w - 1 + 7 : 0] q09;  // signal num: 989
wire signed [ 2 * w - 1 + 7 : 0] q12;  // signal num: 990
wire signed [ 2 * w - 1 + 7 : 0] q22;  // signal num: 991
wire signed [ 2 * w - 1 + 7 : 0] q25;  // signal num: 992
wire signed [ 2 * w - 1 + 7 : 0] q36;  // signal num: 993
wire signed [ 2 * w - 1 + 7 : 0] q45;  // signal num: 994
wire signed [ 2 * w - 1 + 7 : 0] t04_i;  // signal num: 995
wire signed [ 2 * w - 1 + 7 : 0] t04_r;  // signal num: 996
wire signed [ 2 * w - 1 + 7 : 0] t08_r;  // signal num: 997
wire signed [ 2 * w - 1 + 7 : 0] t09_i;  // signal num: 998
wire signed [ 2 * w - 1 + 7 : 0] t09_r;  // signal num: 999
wire signed [ 2 * w - 1 + 7 : 0] t12_i;  // signal num: 1000
wire signed [ 2 * w - 1 + 7 : 0] t12_r;  // signal num: 1001
wire signed [ 2 * w - 1 + 6 : 0] t14_i;  // signal num: 1002
wire signed [ 2 * w - 1 + 6 : 0] t14_r;  // signal num: 1003
wire signed [ 2 * w - 1 + 7 : 0] t21_i;  // signal num: 1004
wire signed [ 2 * w - 1 + 7 : 0] t22_i;  // signal num: 1005
wire signed [ 2 * w - 1 + 7 : 0] t22_r;  // signal num: 1006
wire signed [ 2 * w - 1 + 7 : 0] t23_i;  // signal num: 1007
wire signed [ 2 * w - 1 + 7 : 0] t23_r;  // signal num: 1008
wire signed [ 2 * w - 1 + 7 : 0] t24_i;  // signal num: 1009
wire signed [ 2 * w - 1 + 7 : 0] t25_i;  // signal num: 1010
wire signed [ 2 * w - 1 + 7 : 0] t25_r;  // signal num: 1011
wire signed [ 2 * w - 1 + 7 : 0] t33_i;  // signal num: 1012
wire signed [ 2 * w - 1 + 7 : 0] t34_i;  // signal num: 1013
wire signed [ 2 * w - 1 + 7 : 0] t36_i;  // signal num: 1014
wire signed [ 2 * w - 1 + 7 : 0] t37_i;  // signal num: 1015
wire signed [ 2 * w - 1 + 6 : 0] t39_i;  // signal num: 1016
wire signed [ 2 * w - 1 + 7 : 0] t43_i;  // signal num: 1017
wire signed [ 2 * w - 1 + 7 : 0] t43_r;  // signal num: 1018
wire signed [ 2 * w - 1 + 7 : 0] t45_i;  // signal num: 1019
wire signed [ 2 * w - 1 + 7 : 0] t45_r;  // signal num: 1020
//  Number of sub/add operations: 412
assign pAr6pAr10 =  + A_real[1][1] + A_real[2][1];
assign pAr1mAi13 =  + A_real[0][0] - A_imag[3][0];
assign mAr1mAi13 =  - A_real[0][0] - A_imag[3][0];
assign pAr6mAr10 =  + A_real[1][1] - A_real[2][1];
assign mAr8mAr12 =  - A_real[1][3] - A_real[2][3];
assign mAr15mAi3 =  - A_real[3][2] - A_imag[0][2];
assign pAr10mAr6 =  + A_real[2][1] - A_real[1][1];
assign pAr1pAi13 =  + A_real[0][0] + A_imag[3][0];
assign mAi7mAi11 =  - A_imag[1][2] - A_imag[2][2];
assign pAr7mAr11 =  + A_real[1][2] - A_real[2][2];
assign pAr3pAi15 =  + A_real[0][2] + A_imag[3][2];
assign mAr2mAi14 =  - A_real[0][1] - A_imag[3][1];
assign mAr1pAi13 =  - A_real[0][0] + A_imag[3][0];
assign pAr9mAr5 =  + A_real[2][0] - A_real[1][0];
assign pAr9mAr5mAi7mAi11 = pAr9mAr5 + mAi7mAi11;
assign pAr4pAi16 =  + A_real[0][3] + A_imag[3][3];
assign pAi5mAi9 =  + A_imag[1][0] - A_imag[2][0];
assign pAi5mAi9pAr9mAr5 = pAi5mAi9 + pAr9mAr5;
assign pAr2mAi14 =  + A_real[0][1] - A_imag[3][1];
assign pAr2mAi14pAr4pAi16 = pAr2mAi14 + pAr4pAi16;
assign pAr14mAi2 =  + A_real[3][1] - A_imag[0][1];
assign pAr14mAi2mAr2mAi14 = pAr14mAi2 + mAr2mAi14;
assign pAi8mAi12 =  + A_imag[1][3] - A_imag[2][3];
assign mAr6mAr10 =  - A_real[1][1] - A_real[2][1];
assign pAi7mAi11 =  + A_imag[1][2] - A_imag[2][2];
assign pAr12mAr8 =  + A_real[2][3] - A_real[1][3];
assign pAr12mAr8pAi8mAi12 = pAr12mAr8 + pAi8mAi12;
assign pAr7pAr11 =  + A_real[1][2] + A_real[2][2];
assign pAr13pAi1 =  + A_real[3][0] + A_imag[0][0];
assign mAr16mAi4 =  - A_real[3][3] - A_imag[0][3];
assign pAr6mAi13 =  + A_real[1][1] - A_imag[3][0];
assign mAr10pAr6mAi13 =  - A_real[2][1] + pAr6mAi13;
assign mAr1mAr10pAr6mAi13 =  - A_real[0][0] + mAr10pAr6mAi13;
assign pAi6pAi10 =  + A_imag[1][1] + A_imag[2][1];
assign pAi6pAi10pAr13pAi1 = pAi6pAi10 + pAr13pAi1;
assign pAi6pAi10mAr6mAr10 = pAi6pAi10 + mAr6mAr10;
assign pAr15mAi3 =  + A_real[3][2] - A_imag[0][2];
assign mAr4pAi16 =  - A_real[0][3] + A_imag[3][3];
assign mAr4pAi16mAr16mAi4 = mAr4pAi16 + mAr16mAi4;
assign mAr14mAi2 =  - A_real[3][1] - A_imag[0][1];
assign mAr14mAi2pAr2mAi14 = mAr14mAi2 + pAr2mAi14;
assign pAi10mAi6 =  + A_imag[2][1] - A_imag[1][1];
assign pAi10mAi6pAr6mAr10 = pAi10mAi6 + pAr6mAr10;
assign pAr5pAr9 =  + A_real[1][0] + A_real[2][0];
assign mAr16pAi4 =  - A_real[3][3] + A_imag[0][3];
assign mAr16pAi4pAr4pAi16 = mAr16pAi4 + pAr4pAi16;
assign mAr3pAi15 =  - A_real[0][2] + A_imag[3][2];
assign mAr3pAi15mAr15mAi3 = mAr3pAi15 + mAr15mAi3;
assign mAr3pAi15pAr12mAr8 = mAr3pAi15 + pAr12mAr8;
assign pAr6pAr10mAr3pAi15pAr12mAr8 = pAr6pAr10 + mAr3pAi15pAr12mAr8;
assign mAi5mAi9 =  - A_imag[1][0] - A_imag[2][0];
assign mAi5mAi9pAr5pAr9 = mAi5mAi9 + pAr5pAr9;
assign mAi8mAi12 =  - A_imag[1][3] - A_imag[2][3];
assign mAi8mAi12mAr8mAr12 = mAi8mAi12 + mAr8mAr12;
assign mAi8mAi12pAr1mAi13 = mAi8mAi12 + pAr1mAi13;
assign pAr3mAi15 =  + A_real[0][2] - A_imag[3][2];
assign pAr3mAi15mAi8mAi12 = pAr3mAi15 + mAi8mAi12;
assign pAr3mAi15mAr15mAi3 = pAr3mAi15 + mAr15mAi3;
assign pAi6pAi10mAr6mAr10pAr3mAi15mAr15mAi3 = pAi6pAi10mAr6mAr10 + pAr3mAi15mAr15mAi3;
assign pAr14pAi2 =  + A_real[3][1] + A_imag[0][1];
assign pAr14pAi2pAr2mAi14 = pAr14pAi2 + pAr2mAi14;
assign pAr14pAi2mAr16pAi4 = pAr14pAi2 + mAr16pAi4;
assign pAr4mAi16 =  + A_real[0][3] - A_imag[3][3];
assign pAr4mAi16mAr2mAi14 = pAr4mAi16 + mAr2mAi14;
assign pAr4mAi16mAr16mAi4 = pAr4mAi16 + mAr16mAi4;
assign mAr13pAi1 =  - A_real[3][0] + A_imag[0][0];
assign mAr13pAi1pAi10mAi6 = mAr13pAi1 + pAi10mAi6;
assign mAr13pAi1pAr1pAi13 = mAr13pAi1 + pAr1pAi13;
assign pAi6pAi10mAr6mAr10mAr13pAi1pAr1pAi13 = pAi6pAi10mAr6mAr10 + mAr13pAi1pAr1pAi13;
assign pAi6mAi10 =  + A_imag[1][1] - A_imag[2][1];
assign pAi6mAi10mAr1pAi13 = pAi6mAi10 + mAr1pAi13;
assign pAi6mAi10pAr15mAi3 = pAi6mAi10 + pAr15mAi3;
assign pAr13mAi1 =  + A_real[3][0] - A_imag[0][0];
assign pAr13mAi1pAr6pAr10mAr3pAi15pAr12mAr8 = pAr13mAi1 + pAr6pAr10mAr3pAi15pAr12mAr8;
assign pAr13mAi1pAi6pAi10 = pAr13mAi1 + pAi6pAi10;
assign pAr13mAi1pAr1pAi13 = pAr13mAi1 + pAr1pAi13;
assign pAr13mAi1pAi6mAi10 = pAr13mAi1 + pAi6mAi10;
assign pAr13mAi1pAi6mAi10mAr1mAr10pAr6mAi13 = pAr13mAi1pAi6mAi10 + mAr1mAr10pAr6mAi13;
assign pAr16pAi4 =  + A_real[3][3] + A_imag[0][3];
assign pAr16pAi4mAr4pAi16 = pAr16pAi4 + mAr4pAi16;
assign pAr16pAi4pAr14mAi2 = pAr16pAi4 + pAr14mAi2;
assign pAr16pAi4pAr4mAi16 = pAr16pAi4 + pAr4mAi16;
assign pAi5pAi9 =  + A_imag[1][0] + A_imag[2][0];
assign pAi5pAi9pAr7mAr11 = pAi5pAi9 + pAr7mAr11;
assign pAi5pAi9pAr5pAr9 = pAi5pAi9 + pAr5pAr9;
assign pAi8pAi12 =  + A_imag[1][3] + A_imag[2][3];
assign pAi8pAi12mAr1pAi13 = pAi8pAi12 + mAr1pAi13;
assign pAi8pAi12mAr1pAi13pAi6mAi10pAr15mAi3 = pAi8pAi12mAr1pAi13 + pAi6mAi10pAr15mAi3;
assign pAi8pAi12mAr3pAi15 = pAi8pAi12 + mAr3pAi15;
assign pAi8pAi12mAr3pAi15mAr13pAi1pAi10mAi6 = pAi8pAi12mAr3pAi15 + mAr13pAi1pAi10mAi6;
assign mAi6mAi10 =  - A_imag[1][1] - A_imag[2][1];
assign mAi6mAi10mAr6mAr10 = mAi6mAi10 + mAr6mAr10;
assign mAi6mAi10mAr1mAi13 = mAi6mAi10 + mAr1mAi13;
assign mAr13mAi1 =  - A_real[3][0] - A_imag[0][0];
assign mAi13mAr13mAi1 =  - A_imag[3][0] + mAr13mAi1;
assign mAr13mAi1pAr10mAr6 = mAr13mAi1 + pAr10mAr6;
assign mAr13mAi1pAr10mAr6pAi6mAi10mAr1pAi13 = mAr13mAi1pAr10mAr6 + pAi6mAi10mAr1pAi13;
assign pAr8mAr12 =  + A_real[1][3] - A_real[2][3];
assign pAr8mAr12pAr1mAi13 = pAr8mAr12 + pAr1mAi13;
assign pAr8mAr12pAr3mAi15 = pAr8mAr12 + pAr3mAi15;
assign mAr13pAi1pAr8mAr12pAr3mAi15 = mAr13pAi1 + pAr8mAr12pAr3mAi15;
assign mAr6mAr10mAr13pAi1pAr8mAr12pAr3mAi15 = mAr6mAr10 + mAr13pAi1pAr8mAr12pAr3mAi15;
assign pAi12mAi8 =  + A_imag[2][3] - A_imag[1][3];
assign pAi12mAi8pAr8mAr12 = pAi12mAi8 + pAr8mAr12;
assign mAr13pAi1pAr1pAi13pAi12mAi8pAr8mAr12 = mAr13pAi1pAr1pAi13 + pAi12mAi8pAr8mAr12;
assign pAi12mAi8pAr12mAr8 = pAi12mAi8 + pAr12mAr8;
assign pAi12mAi8pAr12mAr8mAr16pAi4pAr4pAi16 = pAi12mAi8pAr12mAr8 + mAr16pAi4pAr4pAi16;
assign pAi12mAi8mAr13mAi1 = pAi12mAi8 + mAr13mAi1;
assign pAi11mAi7 =  + A_imag[2][2] - A_imag[1][2];
assign pAi11mAi7pAr7mAr11 = pAi11mAi7 + pAr7mAr11;
assign pAi11mAi7pAr5pAr9 = pAi11mAi7 + pAr5pAr9;
assign mAr5mAr9 =  - A_real[1][0] - A_real[2][0];
assign mAr5mAr9mAi5mAi9 = mAr5mAr9 + mAi5mAi9;
assign mAr5mAr9pAi7mAi11 = mAr5mAr9 + pAi7mAi11;
assign mAr5mAr9pAi7mAi11pAr16pAi4pAr14mAi2 = mAr5mAr9pAi7mAi11 + pAr16pAi4pAr14mAi2;
assign mAr5mAr9pAi5pAi9 = mAr5mAr9 + pAi5pAi9;
assign mAr15pAi3 =  - A_real[3][2] + A_imag[0][2];
assign mAr15pAi3pAr3pAi15 = mAr15pAi3 + pAr3pAi15;
assign mAr15pAi3pAr3pAi15pAr12mAr8pAi8mAi12 = mAr15pAi3pAr3pAi15 + pAr12mAr8pAi8mAi12;
assign pAr16mAi4 =  + A_real[3][3] - A_imag[0][3];
assign pAr16mAi4mAr14mAi2 = pAr16mAi4 + mAr14mAi2;
assign pAr16mAi4pAr4pAi16 = pAr16mAi4 + pAr4pAi16;
assign mAr5mAr9pAi5pAi9pAr16mAi4pAr4pAi16 = mAr5mAr9pAi5pAi9 + pAr16mAi4pAr4pAi16;
assign mAr14pAi2 =  - A_real[3][1] + A_imag[0][1];
assign mAr14pAi2mAr16mAi4 = mAr14pAi2 + mAr16mAi4;
assign mAr5mAr9pAi7mAi11mAr14pAi2mAr16mAi4 = mAr5mAr9pAi7mAi11 + mAr14pAi2mAr16mAi4;
assign mAr14pAi2mAr2mAi14 = mAr14pAi2 + mAr2mAi14;
assign pAi9mAi5 =  + A_imag[2][0] - A_imag[1][0];
assign pAi9mAi5pAr9mAr5 = pAi9mAi5 + pAr9mAr5;
assign pAi9mAi5pAr7pAr11 = pAi9mAi5 + pAr7pAr11;
assign mAr2pAi14 =  - A_real[0][1] + A_imag[3][1];
assign mAr2pAi14mAr14mAi2 = mAr2pAi14 + mAr14mAi2;
assign mAr2pAi14pAr14pAi2 = mAr2pAi14 + pAr14pAi2;
assign pAr5mAr9 =  + A_real[1][0] - A_real[2][0];
assign pAr5mAr9pAi5mAi9 = pAr5mAr9 + pAi5mAi9;
assign pAr5mAr9pAi9mAi5 = pAr5mAr9 + pAi9mAi5;
assign pAr5mAr9pAi9mAi5mAr14pAi2mAr2mAi14 = pAr5mAr9pAi9mAi5 + mAr14pAi2mAr2mAi14;
assign mAr7mAr11 =  - A_real[1][2] - A_real[2][2];
assign mAr7mAr11pAi5mAi9 = mAr7mAr11 + pAi5mAi9;
assign mAr7mAr11mAi7mAi11 = mAr7mAr11 + mAi7mAi11;
assign mAr7mAr11mAi7mAi11pAr14pAi2pAr2mAi14 = mAr7mAr11mAi7mAi11 + pAr14pAi2pAr2mAi14;
assign pAr11mAr7 =  + A_real[2][2] - A_real[1][2];
assign pAr11mAr7pAi7mAi11 = pAr11mAr7 + pAi7mAi11;
assign pAr11mAr7mAi5mAi9 = pAr11mAr7 + mAi5mAi9;
assign pAr11mAr7pAi11mAi7 = pAr11mAr7 + pAi11mAi7;
assign pAi10mAi6pAr11mAr7pAi11mAi7 = pAi10mAi6 + pAr11mAr7pAi11mAi7;
assign pAr11mAr7pAi11mAi7pAr16pAi4pAr4mAi16 = pAr11mAr7pAi11mAi7 + pAr16pAi4pAr4mAi16;
assign mAr4mAi16 =  - A_real[0][3] - A_imag[3][3];
assign mAr4mAi16mAr16pAi4 = mAr4mAi16 + mAr16pAi4;
assign mAr4mAi16pAr16mAi4 = mAr4mAi16 + pAr16mAi4;
assign mAr4mAi16pAr16mAi4pAi12mAi8pAr12mAr8 = mAr4mAi16pAr16mAi4 + pAi12mAi8pAr12mAr8;
assign pAi5pAi9pAr5pAr9mAr4mAi16pAr16mAi4 = pAi5pAi9pAr5pAr9 + mAr4mAi16pAr16mAi4;
assign mAr4mAi16mAr2pAi14 = mAr4mAi16 + mAr2pAi14;
assign mAr4mAi16mAr2pAi14pAi9mAi5pAr7pAr11 = mAr4mAi16mAr2pAi14 + pAi9mAi5pAr7pAr11;
assign mAr4mAi16mAr2pAi14pAi11mAi7pAr5pAr9 = mAr4mAi16mAr2pAi14 + pAi11mAi7pAr5pAr9;
assign pAi7pAi11 =  + A_imag[1][2] + A_imag[2][2];
assign pAi7pAi11pAr7pAr11 = pAi7pAi11 + pAr7pAr11;
assign pAi7pAi11pAr7pAr11pAr4mAi16mAr16mAi4 = pAi7pAi11pAr7pAr11 + pAr4mAi16mAr16mAi4;
assign pAi7pAi11pAr5mAr9 = pAi7pAi11 + pAr5mAr9;
assign pAi7pAi11mAr7mAr11 = pAi7pAi11 + mAr7mAr11;
assign pAr15mAi3pAi7pAi11mAr7mAr11 = pAr15mAi3 + pAi7pAi11mAr7mAr11;
assign pAi7pAi11mAr7mAr11mAi8mAi12mAr8mAr12 = pAi7pAi11mAr7mAr11 + mAi8mAi12mAr8mAr12;
assign pAr15pAi3 =  + A_real[3][2] + A_imag[0][2];
assign pAr15pAi3pAr10mAr6 = pAr15pAi3 + pAr10mAr6;
assign pAr15pAi3pAi8mAi12 = pAr15pAi3 + pAi8mAi12;
assign pAr15pAi3pAi8mAi12pAr8mAr12pAr3mAi15 = pAr15pAi3pAi8mAi12 + pAr8mAr12pAr3mAi15;
assign pAr15pAi3pAi8mAi12mAi6mAi10mAr1mAi13 = pAr15pAi3pAi8mAi12 + mAi6mAi10mAr1mAi13;
assign mAr3mAi15 =  - A_real[0][2] - A_imag[3][2];
assign mAr3mAi15mAi6mAi10 = mAr3mAi15 + mAi6mAi10;
assign pAi12mAi8mAr13mAi1mAr3mAi15mAi6mAi10 = pAi12mAi8mAr13mAi1 + mAr3mAi15mAi6mAi10;
assign mAr3mAi15pAr6mAr10 = mAr3mAi15 + pAr6mAr10;
assign mAr3mAi15mAr15pAi3 = mAr3mAi15 + mAr15pAi3;
assign mAr3mAi15mAr15pAi3pAr11mAr7pAi7mAi11 = mAr3mAi15mAr15pAi3 + pAr11mAr7pAi7mAi11;
assign mAr3mAi15mAr15pAi3mAi6mAi10mAr6mAr10 = mAr3mAi15mAr15pAi3 + mAi6mAi10mAr6mAr10;
assign pAr2pAi14 =  + A_real[0][1] + A_imag[3][1];
assign pAr2pAi14pAr14mAi2 = pAr2pAi14 + pAr14mAi2;
assign pAr2pAi14mAr14pAi2 = pAr2pAi14 + mAr14pAi2;
assign pAr2pAi14mAr14pAi2mAr15pAi3pAr3pAi15 = pAr2pAi14mAr14pAi2 + mAr15pAi3pAr3pAi15;
assign pAr2pAi14mAr4pAi16 = pAr2pAi14 + mAr4pAi16;
assign pAr2pAi14mAr4pAi16pAi7pAi11pAr5mAr9 = pAr2pAi14mAr4pAi16 + pAi7pAi11pAr5mAr9;
assign pAr2pAi14mAr4pAi16pAr11mAr7mAi5mAi9 = pAr2pAi14mAr4pAi16 + pAr11mAr7mAi5mAi9;
assign pAr1pAr10 =  + A_real[0][0] + A_real[2][1];
assign mAr6pAr1pAr10 =  - A_real[1][1] + pAr1pAr10;
assign mAi13mAr13mAi1mAr6pAr1pAr10 = mAi13mAr13mAi1 + mAr6pAr1pAr10;
assign pAi10mAi6mAi13mAr13mAi1mAr6pAr1pAr10 = pAi10mAi6 + mAi13mAr13mAi1mAr6pAr1pAr10;
assign pAi13pAr1pAr10 =  + A_imag[3][0] + pAr1pAr10;
assign pAi13pAr1pAr10pAr13mAi1pAi6pAi10 = pAi13pAr1pAr10 + pAr13mAi1pAi6pAi10;
assign pAr6pAi13pAr1pAr10pAr13mAi1pAi6pAi10 =  + A_real[1][1] + pAi13pAr1pAr10pAr13mAi1pAi6pAi10;
assign mAr6pAi13pAr1pAr10 =  - A_real[1][1] + pAi13pAr1pAr10;
assign pAr1pAr10pAr6mAi13 = pAr1pAr10 + pAr6mAi13;
assign pAr1pAr10pAr6mAi13pAi6pAi10pAr13pAi1 = pAr1pAr10pAr6mAi13 + pAi6pAi10pAr13pAi1;
assign mAi6mAi10pAr1pAr10pAr6mAi13 = mAi6mAi10 + pAr1pAr10pAr6mAi13;
assign mAr13mAi1mAi6mAi10pAr1pAr10pAr6mAi13 = mAr13mAi1 + mAi6mAi10pAr1pAr10pAr6mAi13;
assign pAr8mAr12pAr1pAr10pAr6mAi13 = pAr8mAr12 + pAr1pAr10pAr6mAi13;
assign pAr15mAi3pAr8mAr12pAr1pAr10pAr6mAi13 = pAr15mAi3 + pAr8mAr12pAr1pAr10pAr6mAi13;
assign pAr8pAr12 =  + A_real[1][3] + A_real[2][3];
assign pAr8pAr12pAi8pAi12 = pAr8pAr12 + pAi8pAi12;
assign pAr8pAr12mAi8mAi12 = pAr8pAr12 + mAi8mAi12;
assign mAr3mAi15mAr15pAi3pAr8pAr12mAi8mAi12 = mAr3mAi15mAr15pAi3 + pAr8pAr12mAi8mAi12;
assign pAr8pAr12pAr15pAi3 = pAr8pAr12 + pAr15pAi3;
assign pAr8pAr12pAr15pAi3mAr6pAi13pAr1pAr10 = pAr8pAr12pAr15pAi3 + mAr6pAi13pAr1pAr10;
assign pAr8pAr12pAr13pAi1 = pAr8pAr12 + pAr13pAi1;
assign pAr8pAr12pAr13pAi1mAi8mAi12pAr1mAi13 = pAr8pAr12pAr13pAi1 + mAi8mAi12pAr1mAi13;
assign pAr8pAr12pAr13pAi1pAi8pAi12mAr1pAi13 = pAr8pAr12pAr13pAi1 + pAi8pAi12mAr1pAi13;
assign pAr8pAr12pAr13pAi1mAr3mAi15pAr6mAr10 = pAr8pAr12pAr13pAi1 + mAr3mAi15pAr6mAr10;
assign a0_r =  + pAi5pAi9pAr5pAr9 + pAr2pAi14pAr14mAi2 + pAr6pAi13pAr1pAr10pAr13mAi1pAi6pAi10;
assign a0_i = ( + pAr2pAi14mAr14pAi2 + mAr5mAr9pAi5pAi9 + pAi6pAi10mAr6mAr10mAr13pAi1pAr1pAi13);
assign a1_r =  + pAr5mAr9pAi9mAi5 + pAr13mAi1pAr1pAi13 + mAr4mAi16pAr16mAi4pAi12mAi8pAr12mAr8;
assign a1_i = ( + pAr5mAr9pAi5mAi9 + pAr16mAi4pAr4pAi16 + mAr13pAi1pAr1pAi13pAi12mAi8pAr8mAr12);
assign a2_r =  + mAr2pAi14 + pAi11mAi7 + pAi6mAi10 + pAr3mAi15;
assign a2_i = ( + mAr14mAi2 + pAr7mAr11 + pAr15pAi3pAr10mAr6);
assign a3_r =  + mAr4mAi16mAr2pAi14pAi9mAi5pAr7pAr11 + mAr6mAr10mAr13pAi1pAr8mAr12pAr3mAi15;
assign a3_i = ( + pAi7pAi11pAr5mAr9 + pAr15pAi3pAi8mAi12mAi6mAi10mAr1mAi13 + pAr16mAi4mAr14mAi2);
assign a4_r =  + mAr14pAi2mAr2mAi14 + mAr5mAr9mAi5mAi9 + pAr6pAi13pAr1pAr10pAr13mAi1pAi6pAi10;
assign a4_i = ( + mAi5mAi9pAr5pAr9 + pAi6pAi10mAr6mAr10mAr13pAi1pAr1pAi13 + pAr14mAi2mAr2mAi14);
assign a5_r =  + pAi7mAi11 + pAr7mAr11 + mAr4mAi16mAr16pAi4 + mAr15pAi3pAr3pAi15pAr12mAr8pAi8mAi12;
assign a5_i = ( + mAr3mAi15mAr15pAi3pAr11mAr7pAi7mAi11 + mAr4mAi16pAr16mAi4pAi12mAi8pAr12mAr8);
assign a6_r =  + pAr4mAi16 + pAr9mAr5 + pAi12mAi8mAr13mAi1;
assign a6_i = ( + pAi9mAi5 + pAr16pAi4 + pAr8mAr12pAr1mAi13);
assign a7_r =  + mAr5mAr9pAi5pAi9 + mAr2pAi14pAr14pAi2 + mAr13mAi1mAi6mAi10pAr1pAr10pAr6mAi13;
assign a7_i = ( + mAr2pAi14mAr14mAi2 + mAr5mAr9mAi5mAi9 + pAr1pAr10pAr6mAi13pAi6pAi10pAr13pAi1);
assign a8_r =  + pAi7pAi11pAr5mAr9 + pAr14pAi2mAr16pAi4 + pAr8pAr12pAr13pAi1mAr3mAi15pAr6mAr10;
assign a8_i = ( + mAr4mAi16mAr2pAi14 + mAr7mAr11pAi5mAi9 + pAi8pAi12mAr1pAi13pAi6mAi10pAr15mAi3);
assign a9_r =  + pAr12mAr8 + mAr1pAi13 + pAi12mAi8mAr13mAi1 + pAr5mAr9pAi9mAi5 + pAr16pAi4mAr4pAi16;
assign a9_i = ( + pAi12mAi8mAr13mAi1 + pAr5mAr9pAi5mAi9 + mAr4pAi16mAr16mAi4 + pAr8mAr12pAr1mAi13);
assign a10_r =  + pAr2pAi14mAr14pAi2 + pAi9mAi5pAr9mAr5 + pAr13mAi1pAi6mAi10mAr1mAr10pAr6mAi13;
assign a10_i = ( + mAr6pAi13pAr1pAr10 + pAr13mAi1pAi6mAi10 + pAr5mAr9pAi9mAi5mAr14pAi2mAr2mAi14);
assign a11_r =  + pAr2pAi14mAr4pAi16 + pAr8pAr12pAr15pAi3mAr6pAi13pAr1pAr10 + pAr9mAr5mAi7mAi11;
assign a11_i = ( + pAi9mAi5pAr7pAr11 + mAr14pAi2mAr16mAi4 + pAi8pAi12mAr3pAi15mAr13pAi1pAi10mAi6);
assign a12_r =  + pAi7pAi11pAr7pAr11 + pAr2pAi14pAr14mAi2 + mAr3mAi15mAr15pAi3mAi6mAi10mAr6mAr10;
assign a12_i = ( + pAr2pAi14mAr14pAi2 + pAr6pAr10 + mAr3mAi15mAi6mAi10 + pAr15mAi3pAi7pAi11mAr7mAr11);
assign a13_r =  + pAi6pAi10 + pAr6pAr10 + mAr3pAi15mAr15mAi3 + mAr7mAr11mAi7mAi11pAr14pAi2pAr2mAi14;
assign a13_i = ( + pAr7pAr11 + mAi7mAi11 + mAr2pAi14pAr14pAi2 + pAi6pAi10mAr6mAr10pAr3mAi15mAr15mAi3);
assign a14_r =  + pAr13mAi1pAi6mAi10 + mAr4mAi16mAr2pAi14pAi11mAi7pAr5pAr9 + pAr3mAi15mAi8mAi12;
assign a14_i = ( + pAr8pAr12pAr15pAi3mAr6pAi13pAr1pAr10 + pAr16mAi4mAr14mAi2 + pAi5pAi9pAr7mAr11);
assign a15_r =  + mAr15mAi3 + mAr8mAr12 + pAi8pAi12mAr3pAi15 + pAi7pAi11pAr7pAr11pAr4mAi16mAr16mAi4;
assign a15_i = ( + pAr16pAi4pAr4mAi16 + pAr3mAi15mAr15mAi3 + pAi7pAi11mAr7mAr11mAi8mAi12mAr8mAr12);
assign a16_r =  + pAi8pAi12mAr1pAi13pAi6mAi10pAr15mAi3 + mAr5mAr9pAi7mAi11pAr16pAi4pAr14mAi2;
assign a16_i = ( + pAr3pAi15 + mAr8mAr12 + mAr13mAi1pAr10mAr6 + pAr2pAi14mAr4pAi16pAr11mAr7mAi5mAi9);
assign a17_r =  + mAi5mAi9pAr5pAr9 + mAr14mAi2pAr2mAi14 + mAr13mAi1mAi6mAi10pAr1pAr10pAr6mAi13;
assign a17_i = ( + pAi5pAi9pAr5pAr9 + pAr14pAi2pAr2mAi14 + pAr1pAr10pAr6mAi13pAi6pAi10pAr13pAi1);
assign a18_r =  + pAr11mAr7mAi5mAi9 + pAr16mAi4mAr14mAi2 + pAi12mAi8mAr13mAi1mAr3mAi15mAi6mAi10;
assign a18_i = ( + pAr15mAi3pAr8mAr12pAr1pAr10pAr6mAi13 + pAi11mAi7pAr5pAr9 + pAr2mAi14pAr4pAi16);
assign a19_r =  + pAr8pAr12pAr15pAi3 + pAi7pAi11pAr7pAr11pAr4mAi16mAr16mAi4 + pAr3mAi15mAi8mAi12;
assign a19_i = ( + pAi7pAi11mAr7mAr11 + pAr8pAr12pAr15pAi3 + pAi8pAi12mAr3pAi15 + pAr16pAi4pAr4mAi16);
assign a20_r =  + pAr2pAi14pAr14mAi2 + pAi10mAi6pAr6mAr10 + mAr3mAi15mAr15pAi3pAr11mAr7pAi7mAi11;
assign a20_i = ( + pAr11mAr7pAi11mAi7 + pAr2pAi14mAr14pAi2 + pAi6mAi10pAr15mAi3 + mAr3mAi15pAr6mAr10);
assign a21_r =  + pAr14pAi2mAr16pAi4 + pAi12mAi8mAr13mAi1mAr3mAi15mAi6mAi10 + pAi5pAi9pAr7mAr11;
assign a21_i = ( + mAr4mAi16mAr2pAi14 + mAr5mAr9pAi7mAi11 + pAr15mAi3pAr8mAr12pAr1pAr10pAr6mAi13);
assign a22_r =  + pAi8pAi12 + mAr13pAi1 + mAr8mAr12 + mAr1mAi13 + pAi5pAi9pAr5pAr9mAr4mAi16pAr16mAi4;
assign a22_i = ( + pAr13mAi1 + mAr1mAi13 + mAi8mAi12mAr8mAr12 + mAr5mAr9pAi5pAi9pAr16mAi4pAr4pAi16);
assign a23_r =  + pAi7pAi11pAr7pAr11 + mAr16pAi4pAr4pAi16 + mAr3mAi15mAr15pAi3pAr8pAr12mAi8mAi12;
assign a23_i = ( + mAr3mAi15 + mAr4mAi16mAr16pAi4 + pAr15mAi3pAi7pAi11mAr7mAr11 + pAr8pAr12pAi8pAi12);
assign a24_r =  + pAr2pAi14mAr4pAi16 + pAr15pAi3pAi8mAi12mAi6mAi10mAr1mAi13 + pAi5pAi9pAr7mAr11;
assign a24_i = ( + mAr5mAr9pAi7mAi11mAr14pAi2mAr16mAi4 + pAr13mAi1pAr6pAr10mAr3pAi15pAr12mAr8);
assign a25_r =  + pAr11mAr7pAi7mAi11 + pAr4mAi16mAr16mAi4 + pAr15pAi3pAi8mAi12pAr8mAr12pAr3mAi15;
assign a25_i = ( + pAr15pAi3pAi8mAi12 + mAr3pAi15pAr12mAr8 + pAr11mAr7pAi11mAi7pAr16pAi4pAr4mAi16);
assign a26_r =  + pAr3pAi15 + pAr10mAr6 + pAi6mAi10pAr15mAi3 + pAr11mAr7pAi7mAi11 + pAr2pAi14pAr14mAi2;
assign a26_i = ( + pAr10mAr6 + pAi10mAi6pAr11mAr7pAi11mAi7 + pAr2pAi14mAr14pAi2mAr15pAi3pAr3pAi15);
assign a27_r =  + mAr5mAr9pAi7mAi11mAr14pAi2mAr16mAi4 + mAr6mAr10mAr13pAi1pAr8mAr12pAr3mAi15;
assign a27_i = ( + pAr11mAr7mAi5mAi9 + pAr15pAi3pAi8mAi12mAi6mAi10mAr1mAi13 + pAr4mAi16mAr2mAi14);
assign a28_r =  + mAr2pAi14mAr14mAi2 + pAi5mAi9pAr9mAr5 + mAr13mAi1pAr10mAr6pAi6mAi10mAr1pAi13;
assign a28_i = ( + pAi9mAi5pAr9mAr5 + mAr14mAi2pAr2mAi14 + pAi10mAi6mAi13mAr13mAi1mAr6pAr1pAr10);
assign a29_r =  + pAr13mAi1pAr1pAi13 + pAi12mAi8pAr12mAr8mAr16pAi4pAr4pAi16 + pAi5mAi9pAr9mAr5;
assign a29_i = ( + mAr4mAi16mAr16pAi4 + pAi9mAi5pAr9mAr5 + mAr13pAi1pAr1pAi13pAi12mAi8pAr8mAr12);
assign a30_r =  + pAi7pAi11mAr7mAr11 + mAr14mAi2pAr2mAi14 + pAi6pAi10mAr6mAr10pAr3mAi15mAr15mAi3;
assign a30_i = ( + pAr15pAi3 + pAr3mAi15 + mAi6mAi10mAr6mAr10 + mAr7mAr11mAi7mAi11pAr14pAi2pAr2mAi14);
assign a31_r =  + pAi7pAi11pAr5mAr9 + pAr8pAr12pAr15pAi3mAr6pAi13pAr1pAr10 + pAr4mAi16mAr2mAi14;
assign a31_i = ( + pAr16pAi4pAr14mAi2 + mAr7mAr11pAi5mAi9 + pAi8pAi12mAr3pAi15mAr13pAi1pAi10mAi6);
assign a32_r =  + pAi12mAi8 + mAr15mAi3 + pAr8mAr12pAr3mAi15 + pAr11mAr7pAi11mAi7pAr16pAi4pAr4mAi16;
assign a32_i = ( + pAr16pAi4mAr4pAi16 + pAi11mAi7pAr7mAr11 + pAr15pAi3pAi8mAi12pAr8mAr12pAr3mAi15);
assign a33_r =  + mAr15pAi3 + pAi10mAi6 + mAi8mAi12pAr1mAi13 + mAr5mAr9pAi7mAi11pAr16pAi4pAr14mAi2;
assign a33_i = ( + pAr8pAr12pAr13pAi1mAr3mAi15pAr6mAr10 + pAr2pAi14mAr4pAi16pAr11mAr7mAi5mAi9);
assign a34_r =  + pAr16mAi4mAr14mAi2 + pAr8pAr12pAr13pAi1mAr3mAi15pAr6mAr10 + pAr9mAr5mAi7mAi11;
assign a34_i = ( + pAi9mAi5pAr7pAr11 + pAi8pAi12mAr1pAi13pAi6mAi10pAr15mAi3 + pAr2mAi14pAr4pAi16);
assign a35_r =  + pAr11mAr7pAi11mAi7 + pAr16mAi4pAr4pAi16 + mAr15pAi3pAr3pAi15pAr12mAr8pAi8mAi12;
assign a35_i = ( + mAr3mAi15mAr15pAi3 + pAi11mAi7pAr7mAr11 + pAi12mAi8pAr12mAr8mAr16pAi4pAr4pAi16);
assign a36_r =  + mAr3pAi15 + mAr2pAi14pAr14pAi2 + pAi10mAi6pAr11mAr7pAi11mAi7 + pAr15pAi3pAr10mAr6;
assign a36_i = ( + mAr2pAi14mAr14mAi2 + pAi11mAi7pAr7mAr11 + mAr3pAi15mAr15mAi3 + pAi10mAi6pAr6mAr10);
assign a37_r =  + pAr11mAr7mAi5mAi9 + pAr8pAr12pAr15pAi3mAr6pAi13pAr1pAr10 + pAr14pAi2mAr16pAi4;
assign a37_i = ( + pAi8pAi12mAr3pAi15mAr13pAi1pAi10mAi6 + mAr4mAi16mAr2pAi14pAi11mAi7pAr5pAr9);
assign a38_r =  + pAi7pAi11mAr7mAr11 + pAi6pAi10mAr6mAr10 + pAr2pAi14mAr14pAi2mAr15pAi3pAr3pAi15;
assign a38_i = ( + mAr7mAr11mAi7mAi11 + mAr14pAi2mAr2mAi14 + mAr3mAi15mAr15pAi3mAi6mAi10mAr6mAr10);
assign a39_r =  + pAr14pAi2mAr16pAi4 + pAr15pAi3pAi8mAi12mAi6mAi10mAr1mAi13 + pAr9mAr5mAi7mAi11;
assign a39_i = ( + mAr4mAi16mAr2pAi14pAi9mAi5pAr7pAr11 + pAr13mAi1pAr6pAr10mAr3pAi15pAr12mAr8);
assign a40_r =  + mAr13pAi1pAi10mAi6 + mAr1mAr10pAr6mAi13 + pAr5mAr9pAi9mAi5mAr14pAi2mAr2mAi14;
assign a40_i = ( + pAr5mAr9pAi5mAi9 + pAr14mAi2mAr2mAi14 + pAr13mAi1pAi6mAi10mAr1mAr10pAr6mAi13);
assign a41_r =  + mAr5mAr9mAi5mAi9 + pAr16pAi4mAr4pAi16 + pAr8pAr12pAr13pAi1mAi8mAi12pAr1mAi13;
assign a41_i = ( + mAi5mAi9pAr5pAr9 + pAr8pAr12pAr13pAi1pAi8pAi12mAr1pAi13 + mAr4pAi16mAr16mAi4);
assign a42_r =  + pAr13mAi1pAr1pAi13 + pAr8pAr12mAi8mAi12 + pAi5pAi9pAr5pAr9mAr4mAi16pAr16mAi4;
assign a42_i = ( + mAr13pAi1pAr1pAi13 + mAr5mAr9pAi5pAi9pAr16mAi4pAr4pAi16 + pAr8pAr12pAi8pAi12);
assign a43_r =  + pAi12mAi8mAr13mAi1mAr3mAi15mAi6mAi10 + pAr2pAi14mAr4pAi16pAi7pAi11pAr5mAr9;
assign a43_i = ( + pAr15mAi3pAr8mAr12pAr1pAr10pAr6mAi13 + mAr14pAi2mAr16mAi4 + mAr7mAr11pAi5mAi9);
assign a44_r =  + mAr15pAi3pAr3pAi15 + mAr4mAi16mAr16pAi4 + pAi7pAi11mAr7mAr11mAi8mAi12mAr8mAr12;
assign a44_i = ( + mAr4mAi16pAr16mAi4 + mAr7mAr11mAi7mAi11 + mAr3mAi15mAr15pAi3pAr8pAr12mAi8mAi12);
assign a45_r =  + pAr5mAr9pAi9mAi5 + pAr14pAi2pAr2mAi14 + mAr13mAi1pAr10mAr6pAi6mAi10mAr1pAi13;
assign a45_i = ( + pAr5mAr9pAi5mAi9 + mAr2pAi14pAr14pAi2 + pAi10mAi6mAi13mAr13mAi1mAr6pAr1pAr10);
assign a46_r =  + pAi5pAi9pAr5pAr9 + pAr4mAi16mAr16mAi4 + pAr8pAr12pAr13pAi1mAi8mAi12pAr1mAi13;
assign a46_i = ( + mAr5mAr9pAi5pAi9 + pAr16pAi4pAr4mAi16 + pAr8pAr12pAr13pAi1pAi8pAi12mAr1pAi13);
assign a47_r =  + pAi9mAi5pAr7pAr11 + pAr16pAi4pAr14mAi2 + pAr15mAi3pAr8mAr12pAr1pAr10pAr6mAi13;
assign a47_i = ( + pAi8mAi12 + pAr3pAi15 + pAi6pAi10pAr13pAi1 + pAr2pAi14mAr4pAi16pAi7pAi11pAr5mAr9);
//  Linear combinations of elements from B;
//  Number of sub/add operations: 291
assign mBr2mBi14 =  - B_real[0][1] - B_imag[3][1];
assign mBr16mBi12 =  - B_real[3][3] - B_imag[2][3];
assign pBr16mBr13 =  + B_real[3][3] - B_real[3][0];
assign mBi6pBr16mBr13 =  - B_imag[1][1] + pBr16mBr13;
assign pBi1mBi6pBr16mBr13 =  + B_imag[0][0] + mBi6pBr16mBr13;
assign mBr16pBi12 =  - B_real[3][3] + B_imag[2][3];
assign pBi9mBi5 =  + B_imag[2][0] - B_imag[1][0];
assign mBr9mBi13 =  - B_real[2][0] - B_imag[3][0];
assign pBi3pBi7 =  + B_imag[0][2] + B_imag[1][2];
assign mBr13pBi9 =  - B_real[3][0] + B_imag[2][0];
assign pBi4mBi8 =  + B_imag[0][3] - B_imag[1][3];
assign pBr4mBr8 =  + B_real[0][3] - B_real[1][3];
assign mBr14pBi5 =  - B_real[3][1] + B_imag[1][0];
assign mBr16mBr14pBi5 =  - B_real[3][3] + mBr14pBi5;
assign mBi2mBr16mBr14pBi5 =  - B_imag[0][1] + mBr16mBr14pBi5;
assign mBr11mBi15 =  - B_real[2][2] - B_imag[3][2];
assign pBr9pBi13 =  + B_real[2][0] + B_imag[3][0];
assign mBr3mBr7 =  - B_real[0][2] - B_real[1][2];
assign mBi10mBi12 =  - B_imag[2][1] - B_imag[2][3];
assign pBi2pBi6 =  + B_imag[0][1] + B_imag[1][1];
assign pBr13pBi1 =  + B_real[3][0] + B_imag[0][0];
assign pBi6pBr13pBi1 =  + B_imag[1][1] + pBr13pBi1;
assign mBr15pBi11 =  - B_real[3][2] + B_imag[2][2];
assign mBi1mBi5 =  - B_imag[0][0] - B_imag[1][0];
assign pBi1pBi5 =  + B_imag[0][0] + B_imag[1][0];
assign pBr7mBr3 =  + B_real[1][2] - B_real[0][2];
assign pBr6pBr7mBr3 =  + B_real[1][1] + pBr7mBr3;
assign pBr5mBr1 =  + B_real[1][0] - B_real[0][0];
assign pBr1mBr5 =  + B_real[0][0] - B_real[1][0];
assign pBr15pBi11 =  + B_real[3][2] + B_imag[2][2];
assign pBi10pBr15pBi11 =  + B_imag[2][1] + pBr15pBi11;
assign pBi12pBi10pBr15pBi11 =  + B_imag[2][3] + pBi10pBr15pBi11;
assign pBi7mBi3 =  + B_imag[1][2] - B_imag[0][2];
assign mBi5pBi7mBi3 =  - B_imag[1][0] + pBi7mBi3;
assign pBi2mBi6 =  + B_imag[0][1] - B_imag[1][1];
assign mBi14mBi16 =  - B_imag[3][1] - B_imag[3][3];
assign mBi14mBi16pBr1mBr5 = mBi14mBi16 + pBr1mBr5;
assign pBr2mBi14mBi16 =  + B_real[0][1] + mBi14mBi16;
assign pBr14pBr16 =  + B_real[3][1] + B_real[3][3];
assign pBr14pBr16pBi12pBi10pBr15pBi11 = pBr14pBr16 + pBi12pBi10pBr15pBi11;
assign pBr14pBr16mBi10mBi12 = pBr14pBr16 + mBi10mBi12;
assign mBr10mBr12 =  - B_real[2][1] - B_real[2][3];
assign mBr10mBr12pBr2mBi14mBi16 = mBr10mBr12 + pBr2mBi14mBi16;
assign pBr10mBr6 =  + B_real[2][1] - B_real[1][1];
assign pBr12pBr10mBr6 =  + B_real[2][3] + pBr10mBr6;
assign mBr2pBi14 =  - B_real[0][1] + B_imag[3][1];
assign pBi16mBr2pBi14 =  + B_imag[3][3] + mBr2pBi14;
assign pBi12mBi9 =  + B_imag[2][3] - B_imag[2][0];
assign mBi2mBi6 =  - B_imag[0][1] - B_imag[1][1];
assign pBi3mBi7 =  + B_imag[0][2] - B_imag[1][2];
assign pBi3mBi7pBi2mBi6 = pBi3mBi7 + pBi2mBi6;
assign pBi5mBi1 =  + B_imag[1][0] - B_imag[0][0];
assign pBi5mBi1pBi3mBi7 = pBi5mBi1 + pBi3mBi7;
assign pBr6pBr10 =  + B_real[1][1] + B_real[2][1];
assign pBr6pBr10mBr2pBi14 = pBr6pBr10 + mBr2pBi14;
assign pBr3mBr7 =  + B_real[0][2] - B_real[1][2];
assign pBr3mBr7pBr10mBr6 = pBr3mBr7 + pBr10mBr6;
assign pBr3mBr7pBr5mBr1 = pBr3mBr7 + pBr5mBr1;
assign pBi4pBi8 =  + B_imag[0][3] + B_imag[1][3];
assign pBi4pBi8pBi2pBi6 = pBi4pBi8 + pBi2pBi6;
assign pBi4pBi8mBi1mBi5 = pBi4pBi8 + mBi1mBi5;
assign mBr1mBr5 =  - B_real[0][0] - B_real[1][0];
assign mBr14pBi10 =  - B_real[3][1] + B_imag[2][1];
assign mBr14pBi10mBr15pBi11 = mBr14pBi10 + mBr15pBi11;
assign pBi6mBi2 =  + B_imag[1][1] - B_imag[0][1];
assign pBi6mBi2pBi7mBi3 = pBi6mBi2 + pBi7mBi3;
assign pBi6mBi2mBr14pBi10 = pBi6mBi2 + mBr14pBi10;
assign pBr14pBi6mBi2 =  + B_real[3][1] + pBi6mBi2;
assign mBr4mBr8 =  - B_real[0][3] - B_real[1][3];
assign mBr4mBr8pBr12pBr10mBr6 = mBr4mBr8 + pBr12pBr10mBr6;
assign mBr6mBr4mBr8 =  - B_real[1][1] + mBr4mBr8;
assign pBr5mBr4mBr8 =  + B_real[1][0] + mBr4mBr8;
assign pBr1pBr5mBr4mBr8 =  + B_real[0][0] + pBr5mBr4mBr8;
assign pBi16mBi13 =  + B_imag[3][3] - B_imag[3][0];
assign pBr2pBi14 =  + B_real[0][1] + B_imag[3][1];
assign pBi16pBr2pBi14 =  + B_imag[3][3] + pBr2pBi14;
assign pBr2pBi14pBr10mBr6 = pBr2pBi14 + pBr10mBr6;
assign mBr10pBr2pBi14 =  - B_real[2][1] + pBr2pBi14;
assign pBr3pBr7 =  + B_real[0][2] + B_real[1][2];
assign pBr13mBr16 =  + B_real[3][0] - B_real[3][3];
assign pBr13mBr16pBi12mBi9 = pBr13mBr16 + pBi12mBi9;
assign pBr13mBr16pBi12mBi9pBi4pBi8mBi1mBi5 = pBr13mBr16pBi12mBi9 + pBi4pBi8mBi1mBi5;
assign pBr11mBi15 =  + B_real[2][2] - B_imag[3][2];
assign pBr10pBr11mBi15 =  + B_real[2][1] + pBr11mBi15;
assign pBi13pBr11mBi15 =  + B_imag[3][0] + pBr11mBi15;
assign mBr9pBi13pBr11mBi15 =  - B_real[2][0] + pBi13pBr11mBi15;
assign mBi16pBi13pBr11mBi15 =  - B_imag[3][3] + pBi13pBr11mBi15;
assign mBi3mBi7 =  - B_imag[0][2] - B_imag[1][2];
assign mBr11pBi15 =  - B_real[2][2] + B_imag[3][2];
assign mBr11pBi15mBr10pBr2pBi14 = mBr11pBi15 + mBr10pBr2pBi14;
assign mBr11pBi15mBr10mBr12 = mBr11pBi15 + mBr10mBr12;
assign pBi16mBr2pBi14mBr11pBi15mBr10mBr12 = pBi16mBr2pBi14 + mBr11pBi15mBr10mBr12;
assign pBr9mBr11pBi15 =  + B_real[2][0] + mBr11pBi15;
assign mBi13pBr9mBr11pBi15 =  - B_imag[3][0] + pBr9mBr11pBi15;
assign pBi16mBi13pBr9mBr11pBi15 = pBi16mBi13 + pBr9mBr11pBi15;
assign mBr12pBi16mBi13pBr9mBr11pBi15 =  - B_real[2][3] + pBi16mBi13pBr9mBr11pBi15;
assign pBr12mBr9 =  + B_real[2][3] - B_real[2][0];
assign pBr12mBr9mBi16pBi13pBr11mBi15 = pBr12mBr9 + mBi16pBi13pBr11mBi15;
assign pBr12mBr9pBi16mBi13 = pBr12mBr9 + pBi16mBi13;
assign pBr4pBr8 =  + B_real[0][3] + B_real[1][3];
assign pBr4pBr8mBr1mBr5 = pBr4pBr8 + mBr1mBr5;
assign pBr4pBr8pBr6pBr10 = pBr4pBr8 + pBr6pBr10;
assign pBr12pBr4pBr8pBr6pBr10 =  + B_real[2][3] + pBr4pBr8pBr6pBr10;
assign pBr4pBr8pBr3pBr7 = pBr4pBr8 + pBr3pBr7;
assign mBr1mBr5pBr4pBr8pBr3pBr7 = mBr1mBr5 + pBr4pBr8pBr3pBr7;
assign pBr6pBr4pBr8pBr3pBr7 =  + B_real[1][1] + pBr4pBr8pBr3pBr7;
assign mBi4mBi8 =  - B_imag[0][3] - B_imag[1][3];
assign mBi4mBi8pBi1pBi5 = mBi4mBi8 + pBi1pBi5;
assign mBi4mBi8mBi3mBi7 = mBi4mBi8 + mBi3mBi7;
assign mBi4mBi8mBi3mBi7pBi1pBi5 = mBi4mBi8mBi3mBi7 + pBi1pBi5;
assign mBi4mBi8mBi3mBi7mBi2mBi6 = mBi4mBi8mBi3mBi7 + mBi2mBi6;
assign mBr15mBi11 =  - B_real[3][2] - B_imag[2][2];
assign mBi10mBr15mBi11 =  - B_imag[2][1] + mBr15mBi11;
assign mBr14mBi10mBr15mBi11 =  - B_real[3][1] + mBi10mBr15mBi11;
assign mBr14mBi10mBr15mBi11mBr16mBi12 = mBr14mBi10mBr15mBi11 + mBr16mBi12;
assign pBi9mBr15mBi11 =  + B_imag[2][0] + mBr15mBi11;
assign pBr13pBi9mBr15mBi11 =  + B_real[3][0] + pBi9mBr15mBi11;
assign mBi12pBi9mBr15mBi11 =  - B_imag[2][3] + pBi9mBr15mBi11;
assign pBr13mBr16mBi12pBi9mBr15mBi11 = pBr13mBr16 + mBi12pBi9mBr15mBi11;
assign b0_r =  + mBr1mBr5 + pBr9pBi13;
assign b0_i = ( + mBi1mBi5 + mBr13pBi9);
assign b1_r =  + mBi4mBi8 + mBi2mBi6 + pBr14pBr16mBi10mBi12;
assign b1_i = ( + pBr12pBr4pBr8pBr6pBr10 + pBi16pBr2pBi14);
assign b2_r =  - B_imag[2][1] + pBr2pBi14pBr10mBr6 + pBr14pBi6mBi2;
assign b2_i = ( + mBr14pBi10 + pBi2mBi6 + pBr2pBi14pBr10mBr6);
assign b3_r =  + mBi10mBr15mBi11 + pBi7mBi3 + pBi6pBr13pBi1;
assign b3_i = ( - B_real[0][0] + pBi13pBr11mBi15 + pBr3mBr7pBr10mBr6);
assign b4_r =  + pBr3mBr7pBr5mBr1 + pBr4mBr8 + mBr12pBi16mBi13pBr9mBr11pBi15;
assign b4_i = ( + pBr13mBr16mBi12pBi9mBr15mBi11 + pBi4mBi8 + pBi5mBi1pBi3mBi7);
assign b5_r =  + pBr2mBi14mBi16 + pBr12pBr4pBr8pBr6pBr10;
assign b5_i = ( + B_imag[2][1] + B_imag[2][3] + pBr14pBr16 + pBi4pBi8pBi2pBi6);
assign b6_r =  - B_imag[2][1] + B_real[3][1] + pBi2mBi6 + pBr6pBr10mBr2pBi14;
assign b6_i = ( + pBr6pBr10mBr2pBi14 + pBi6mBi2mBr14pBi10);
assign b7_r =  + pBr5mBr1 + pBr12mBr9pBi16mBi13 + pBr4mBr8;
assign b7_i = ( + pBi5mBi1 + pBr13mBr16pBi12mBi9 + pBi4mBi8);
assign b8_r =  + B_real[2][1] + B_imag[3][1] + mBr3mBr7 + pBr1pBr5mBr4mBr8;
assign b8_i = ( + mBr14pBi10 + mBi4mBi8mBi3mBi7pBi1pBi5);
assign b9_r =  + mBi4mBi8mBi3mBi7mBi2mBi6 + pBr14pBr16pBi12pBi10pBr15pBi11;
assign b9_i = ( + mBr11pBi15mBr10mBr12 + pBr6pBr4pBr8pBr3pBr7 + pBi16pBr2pBi14);
assign b10_r =  - B_imag[0][3] + B_imag[1][3] + mBr16pBi12 + pBi6mBi2mBr14pBi10;
assign b10_i = ( - B_real[1][1] + pBr4mBr8 + mBr10mBr12pBr2mBi14mBi16);
assign b11_r =  + mBr14mBi10mBr15mBi11 + mBi4mBi8pBi1pBi5;
assign b11_i = ( - B_imag[3][1] + pBr10pBr11mBi15 + pBr4pBr8mBr1mBr5);
assign b12_r =  + mBr1mBr5pBr4pBr8pBr3pBr7 + mBr12pBi16mBi13pBr9mBr11pBi15;
assign b12_i = ( + pBr13mBr16mBi12pBi9mBr15mBi11 + pBi4pBi8mBi1mBi5 + pBi3pBi7);
assign b13_r =  - B_imag[2][0] - B_real[3][0] + pBr15pBi11 + pBi5mBi1pBi3mBi7;
assign b13_i = ( + pBr1mBr5 + pBr7mBr3 + mBi13pBr9mBr11pBi15);
assign b14_r =  - B_real[1][0] + B_real[2][0] + mBr2mBi14;
assign b14_i = ( - B_imag[0][1] + B_real[3][1] + pBi9mBi5);
assign b15_r =  + B_imag[2][0] - B_imag[2][3] + pBr13mBr16 + pBi4pBi8mBi1mBi5;
assign b15_i = ( + B_imag[3][0] - B_imag[3][3] + pBr12mBr9 + pBr1pBr5mBr4mBr8);
assign b16_r =  + B_real[1][0] + pBr2pBi14 + pBr3mBr7 + pBr9mBr11pBi15;
assign b16_i = ( + B_imag[0][1] + pBi9mBr15mBi11 + pBi3mBi7 + mBr14pBi5);
assign b17_r =  + mBi3mBi7 + pBi1pBi5 + pBr13pBi9mBr15mBi11;
assign b17_i = ( + pBr3pBr7 + mBr1mBr5 + mBr9pBi13pBr11mBi15);
assign b18_r =  + pBr13pBi9mBr15mBi11 + pBi4pBi8pBi2pBi6;
assign b18_i = ( - B_real[0][1] + mBr6mBr4mBr8 + mBr9pBi13pBr11mBi15);
assign b19_r =  - B_imag[2][0] + mBr15pBi11 + pBr13pBi1 + mBi5pBi7mBi3;
assign b19_i = ( + pBr3mBr7pBr5mBr1 + pBr9pBi13 + mBr11mBi15);
assign b20_r =  + pBr14pBr16mBi10mBi12 + pBi4pBi8pBi2pBi6;
assign b20_i = ( + pBi16mBr2pBi14 + mBr4mBr8pBr12pBr10mBr6);
assign b21_r =  - B_real[0][1] + pBr12mBr9pBi16mBi13 + pBr6pBr7mBr3;
assign b21_i = ( + pBr13mBr16pBi12mBi9 + pBi6mBi2pBi7mBi3);
assign b22_r =  + pBr13mBr16mBi12pBi9mBr15mBi11 + mBi4mBi8mBi3mBi7pBi1pBi5;
assign b22_i = ( + mBr1mBr5pBr4pBr8pBr3pBr7 + pBr12mBr9mBi16pBi13pBr11mBi15);
assign b23_r =  + B_real[2][2] + B_imag[3][2] + pBr12mBr9pBi16mBi13 + mBr1mBr5pBr4pBr8pBr3pBr7;
assign b23_i = ( + mBr15pBi11 + pBi3pBi7 + pBr13mBr16pBi12mBi9pBi4pBi8mBi1mBi5);
assign b24_r =  + pBi6mBi2 + pBr13mBr16mBi12pBi9mBr15mBi11;
assign b24_i = ( + B_real[0][1] - B_real[1][1] + pBr12mBr9mBi16pBi13pBr11mBi15);
assign b25_r =  + mBi4mBi8mBi3mBi7mBi2mBi6 + mBr16pBi12 + mBr14pBi10mBr15pBi11;
assign b25_i = ( + pBr6pBr4pBr8pBr3pBr7 + mBr11mBi15 + mBr10mBr12pBr2mBi14mBi16);
assign b26_r =  - B_real[1][1] + pBr3mBr7 + mBr11pBi15mBr10pBr2pBi14;
assign b26_i = ( + mBr14mBi10mBr15mBi11 + pBi3mBi7pBi2mBi6);
assign b27_r =  + mBi4mBi8mBi3mBi7mBi2mBi6 + mBr13pBi9;
assign b27_i = ( + B_real[0][1] + pBr6pBr4pBr8pBr3pBr7 + mBr9mBi13);
assign b28_r =  + pBr2pBi14 + pBr6pBr10;
assign b28_i = ( + mBr14pBi10 + pBi2pBi6);
assign b29_r =  + mBr14mBi10mBr15mBi11 + pBi6mBi2pBi7mBi3;
assign b29_i = ( + B_real[0][1] - B_imag[3][1] + pBr11mBi15 + pBr3mBr7pBr10mBr6);
assign b30_r =  + pBr12mBr9pBi16mBi13 + pBr4pBr8mBr1mBr5;
assign b30_i = ( + pBr13mBr16pBi12mBi9pBi4pBi8mBi1mBi5);
assign b31_r =  + mBr10mBr12 + pBr7mBr3 + mBi14mBi16pBr1mBr5;
assign b31_i = ( + B_imag[0][0] + pBr14pBr16mBi10mBi12 + mBi5pBi7mBi3);
assign b32_r =  + B_imag[2][1] + pBr14pBi6mBi2;
assign b32_i = ( - B_real[1][1] + mBr10pBr2pBi14);
assign b33_r =  + pBr12mBr9 + pBr5mBr4mBr8 + pBi16mBr2pBi14;
assign b33_i = ( + mBi4mBi8 + pBi12mBi9 + mBi2mBr16mBr14pBi5);
assign b34_r =  + pBi5mBi1 + mBr14mBi10mBr15mBi11mBr16mBi12;
assign b34_i = ( + B_real[2][3] + pBr10pBr11mBi15 + mBi14mBi16pBr1mBr5);
assign b35_r =  + pBi3mBi7pBi2mBi6 + mBr14pBi10mBr15pBi11;
assign b35_i = ( - B_real[2][1] + mBr11mBi15 + pBr6pBr7mBr3 + mBr2mBi14);
assign b36_r =  + mBr3mBr7 + mBr6mBr4mBr8 + pBi16mBr2pBi14mBr11pBi15mBr10mBr12;
assign b36_i = ( + mBi4mBi8mBi3mBi7mBi2mBi6 + mBr14mBi10mBr15mBi11mBr16mBi12);
assign b37_r =  + mBi4mBi8mBi3mBi7 + mBi12pBi9mBr15mBi11 + mBi2mBr16mBr14pBi5;
assign b37_i = ( - B_real[1][0] + pBr12mBr9 + pBr11mBi15 + pBr4pBr8pBr3pBr7 + pBr2mBi14mBi16);
assign b38_r =  + pBi5mBi1 + mBr13pBi9;
assign b38_i = ( + pBr1mBr5 + mBr9mBi13);
assign b39_r =  + mBi4mBi8 + mBi10mBi12 + pBi1mBi6pBr16mBr13;
assign b39_i = ( - B_real[0][0] + pBi16mBi13 + pBr12pBr4pBr8pBr6pBr10);
assign b40_r =  + B_real[3][1] + mBi3mBi7 + mBi2mBi6 + pBi10pBr15pBi11;
assign b40_i = ( + B_real[1][1] + pBr3pBr7 + mBr11pBi15mBr10pBr2pBi14);
assign b41_r =  + pBr12mBr9pBi16mBi13 + pBr1pBr5mBr4mBr8;
assign b41_i = ( + pBr13mBr16pBi12mBi9 + mBi4mBi8pBi1pBi5);
assign b42_r =  - B_imag[2][0] + B_real[3][0] + pBi5mBi1;
assign b42_i = ( + pBr1mBr5 + pBr9pBi13);
assign b43_r =  + B_real[0][0] + mBr3mBr7 + mBi16pBi13pBr11mBi15 + mBr4mBr8pBr12pBr10mBr6;
assign b43_i = ( + mBi4mBi8mBi3mBi7 + pBi12pBi10pBr15pBi11 + pBi1mBi6pBr16mBr13);
assign b44_r =  + pBr13pBi1 + pBi9mBi5;
assign b44_i = ( - B_real[2][0] + B_imag[3][0] + pBr5mBr1);
assign b45_r =  + pBi3mBi7pBi2mBi6 + pBi4mBi8 + pBr14pBr16pBi12pBi10pBr15pBi11;
assign b45_i = ( - B_real[0][3] + B_real[1][3] + pBr6pBr7mBr3 + pBi16mBr2pBi14mBr11pBi15mBr10mBr12);
assign b46_r =  + pBr3mBr7pBr5mBr1 + mBi13pBr9mBr11pBi15;
assign b46_i = ( + pBr13pBi9mBr15mBi11 + pBi5mBi1pBi3mBi7);
assign b47_r =  + B_real[0][0] - B_imag[3][0] + pBr6pBr10;
assign b47_i = ( + B_imag[2][1] + pBi6pBr13pBi1);
//  Perform the 48 multiplications;
//  Karatsuba-style complex multiplication
assign p00 = a0_r  *  b0_r;
assign q00 = a0_i  *  b0_i;
assign k00 = (a0_r  +  a0_i)  *  (b0_r  +  b0_i);
assign t00_r = p00  -  q00;
assign t00_i = k00  -  p00  -  q00;
assign p01 = a1_r  *  b1_r;
assign q01 = a1_i  *  b1_i;
assign k01 = (a1_r  +  a1_i)  *  (b1_r  +  b1_i);
assign t01_r = p01  -  q01;
assign t01_i = k01  -  p01  -  q01;
assign p02 = a2_r  *  b2_r;
assign q02 = a2_i  *  b2_i;
assign k02 = (a2_r  +  a2_i)  *  (b2_r  +  b2_i);
assign t02_r = p02  -  q02;
assign t02_i = k02  -  p02  -  q02;
assign p03 = a3_r  *  b3_r;
assign q03 = a3_i  *  b3_i;
assign k03 = (a3_r  +  a3_i)  *  (b3_r  +  b3_i);
assign t03_r = p03  -  q03;
assign t03_i = k03  -  p03  -  q03;
assign p04 = a4_r  *  b4_r;
assign q04 = a4_i  *  b4_i;
assign k04 = (a4_r  +  a4_i)  *  (b4_r  +  b4_i);
assign t04_r = p04  -  q04;
assign t04_i = k04  -  p04  -  q04;
assign p05 = a5_r  *  b5_r;
assign q05 = a5_i  *  b5_i;
assign k05 = (a5_r  +  a5_i)  *  (b5_r  +  b5_i);
assign t05_r = p05  -  q05;
assign t05_i = k05  -  p05  -  q05;
assign p06 = a6_r  *  b6_r;
assign q06 = a6_i  *  b6_i;
assign k06 = (a6_r  +  a6_i)  *  (b6_r  +  b6_i);
assign t06_r = p06  -  q06;
assign t06_i = k06  -  p06  -  q06;
assign p07 = a7_r  *  b7_r;
assign q07 = a7_i  *  b7_i;
assign k07 = (a7_r  +  a7_i)  *  (b7_r  +  b7_i);
assign t07_r = p07  -  q07;
assign t07_i = k07  -  p07  -  q07;
assign p08 = a8_r  *  b8_r;
assign q08 = a8_i  *  b8_i;
assign k08 = (a8_r  +  a8_i)  *  (b8_r  +  b8_i);
assign t08_r = p08  -  q08;
assign t08_i = k08  -  p08  -  q08;
assign p09 = a9_r  *  b9_r;
assign q09 = a9_i  *  b9_i;
assign k09 = (a9_r  +  a9_i)  *  (b9_r  +  b9_i);
assign t09_r = p09  -  q09;
assign t09_i = k09  -  p09  -  q09;
assign p10 = a10_r  *  b10_r;
assign q10 = a10_i  *  b10_i;
assign k10 = (a10_r  +  a10_i)  *  (b10_r  +  b10_i);
assign t10_r = p10  -  q10;
assign t10_i = k10  -  p10  -  q10;
assign p11 = a11_r  *  b11_r;
assign q11 = a11_i  *  b11_i;
assign k11 = (a11_r  +  a11_i)  *  (b11_r  +  b11_i);
assign t11_r = p11  -  q11;
assign t11_i = k11  -  p11  -  q11;
assign p12 = a12_r  *  b12_r;
assign q12 = a12_i  *  b12_i;
assign k12 = (a12_r  +  a12_i)  *  (b12_r  +  b12_i);
assign t12_r = p12  -  q12;
assign t12_i = k12  -  p12  -  q12;
assign p13 = a13_r  *  b13_r;
assign q13 = a13_i  *  b13_i;
assign k13 = (a13_r  +  a13_i)  *  (b13_r  +  b13_i);
assign t13_r = p13  -  q13;
assign t13_i = k13  -  p13  -  q13;
assign p14 = a14_r  *  b14_r;
assign q14 = a14_i  *  b14_i;
assign k14 = (a14_r  +  a14_i)  *  (b14_r  +  b14_i);
assign t14_r = p14  -  q14;
assign t14_i = k14  -  p14  -  q14;
assign p15 = a15_r  *  b15_r;
assign q15 = a15_i  *  b15_i;
assign k15 = (a15_r  +  a15_i)  *  (b15_r  +  b15_i);
assign t15_r = p15  -  q15;
assign t15_i = k15  -  p15  -  q15;
assign p16 = a16_r  *  b16_r;
assign q16 = a16_i  *  b16_i;
assign k16 = (a16_r  +  a16_i)  *  (b16_r  +  b16_i);
assign t16_r = p16  -  q16;
assign t16_i = k16  -  p16  -  q16;
assign p17 = a17_r  *  b17_r;
assign q17 = a17_i  *  b17_i;
assign k17 = (a17_r  +  a17_i)  *  (b17_r  +  b17_i);
assign t17_r = p17  -  q17;
assign t17_i = k17  -  p17  -  q17;
assign p18 = a18_r  *  b18_r;
assign q18 = a18_i  *  b18_i;
assign k18 = (a18_r  +  a18_i)  *  (b18_r  +  b18_i);
assign t18_r = p18  -  q18;
assign t18_i = k18  -  p18  -  q18;
assign p19 = a19_r  *  b19_r;
assign q19 = a19_i  *  b19_i;
assign k19 = (a19_r  +  a19_i)  *  (b19_r  +  b19_i);
assign t19_r = p19  -  q19;
assign t19_i = k19  -  p19  -  q19;
assign p20 = a20_r  *  b20_r;
assign q20 = a20_i  *  b20_i;
assign k20 = (a20_r  +  a20_i)  *  (b20_r  +  b20_i);
assign t20_r = p20  -  q20;
assign t20_i = k20  -  p20  -  q20;
assign p21 = a21_r  *  b21_r;
assign q21 = a21_i  *  b21_i;
assign k21 = (a21_r  +  a21_i)  *  (b21_r  +  b21_i);
assign t21_r = p21  -  q21;
assign t21_i = k21  -  p21  -  q21;
assign p22 = a22_r  *  b22_r;
assign q22 = a22_i  *  b22_i;
assign k22 = (a22_r  +  a22_i)  *  (b22_r  +  b22_i);
assign t22_r = p22  -  q22;
assign t22_i = k22  -  p22  -  q22;
assign p23 = a23_r  *  b23_r;
assign q23 = a23_i  *  b23_i;
assign k23 = (a23_r  +  a23_i)  *  (b23_r  +  b23_i);
assign t23_r = p23  -  q23;
assign t23_i = k23  -  p23  -  q23;
assign p24 = a24_r  *  b24_r;
assign q24 = a24_i  *  b24_i;
assign k24 = (a24_r  +  a24_i)  *  (b24_r  +  b24_i);
assign t24_r = p24  -  q24;
assign t24_i = k24  -  p24  -  q24;
assign p25 = a25_r  *  b25_r;
assign q25 = a25_i  *  b25_i;
assign k25 = (a25_r  +  a25_i)  *  (b25_r  +  b25_i);
assign t25_r = p25  -  q25;
assign t25_i = k25  -  p25  -  q25;
assign p26 = a26_r  *  b26_r;
assign q26 = a26_i  *  b26_i;
assign k26 = (a26_r  +  a26_i)  *  (b26_r  +  b26_i);
assign t26_r = p26  -  q26;
assign t26_i = k26  -  p26  -  q26;
assign p27 = a27_r  *  b27_r;
assign q27 = a27_i  *  b27_i;
assign k27 = (a27_r  +  a27_i)  *  (b27_r  +  b27_i);
assign t27_r = p27  -  q27;
assign t27_i = k27  -  p27  -  q27;
assign p28 = a28_r  *  b28_r;
assign q28 = a28_i  *  b28_i;
assign k28 = (a28_r  +  a28_i)  *  (b28_r  +  b28_i);
assign t28_r = p28  -  q28;
assign t28_i = k28  -  p28  -  q28;
assign p29 = a29_r  *  b29_r;
assign q29 = a29_i  *  b29_i;
assign k29 = (a29_r  +  a29_i)  *  (b29_r  +  b29_i);
assign t29_r = p29  -  q29;
assign t29_i = k29  -  p29  -  q29;
assign p30 = a30_r  *  b30_r;
assign q30 = a30_i  *  b30_i;
assign k30 = (a30_r  +  a30_i)  *  (b30_r  +  b30_i);
assign t30_r = p30  -  q30;
assign t30_i = k30  -  p30  -  q30;
assign p31 = a31_r  *  b31_r;
assign q31 = a31_i  *  b31_i;
assign k31 = (a31_r  +  a31_i)  *  (b31_r  +  b31_i);
assign t31_r = p31  -  q31;
assign t31_i = k31  -  p31  -  q31;
assign p32 = a32_r  *  b32_r;
assign q32 = a32_i  *  b32_i;
assign k32 = (a32_r  +  a32_i)  *  (b32_r  +  b32_i);
assign t32_r = p32  -  q32;
assign t32_i = k32  -  p32  -  q32;
assign p33 = a33_r  *  b33_r;
assign q33 = a33_i  *  b33_i;
assign k33 = (a33_r  +  a33_i)  *  (b33_r  +  b33_i);
assign t33_r = p33  -  q33;
assign t33_i = k33  -  p33  -  q33;
assign p34 = a34_r  *  b34_r;
assign q34 = a34_i  *  b34_i;
assign k34 = (a34_r  +  a34_i)  *  (b34_r  +  b34_i);
assign t34_r = p34  -  q34;
assign t34_i = k34  -  p34  -  q34;
assign p35 = a35_r  *  b35_r;
assign q35 = a35_i  *  b35_i;
assign k35 = (a35_r  +  a35_i)  *  (b35_r  +  b35_i);
assign t35_r = p35  -  q35;
assign t35_i = k35  -  p35  -  q35;
assign p36 = a36_r  *  b36_r;
assign q36 = a36_i  *  b36_i;
assign k36 = (a36_r  +  a36_i)  *  (b36_r  +  b36_i);
assign t36_r = p36  -  q36;
assign t36_i = k36  -  p36  -  q36;
assign p37 = a37_r  *  b37_r;
assign q37 = a37_i  *  b37_i;
assign k37 = (a37_r  +  a37_i)  *  (b37_r  +  b37_i);
assign t37_r = p37  -  q37;
assign t37_i = k37  -  p37  -  q37;
assign p38 = a38_r  *  b38_r;
assign q38 = a38_i  *  b38_i;
assign k38 = (a38_r  +  a38_i)  *  (b38_r  +  b38_i);
assign t38_r = p38  -  q38;
assign t38_i = k38  -  p38  -  q38;
assign p39 = a39_r  *  b39_r;
assign q39 = a39_i  *  b39_i;
assign k39 = (a39_r  +  a39_i)  *  (b39_r  +  b39_i);
assign t39_r = p39  -  q39;
assign t39_i = k39  -  p39  -  q39;
assign p40 = a40_r  *  b40_r;
assign q40 = a40_i  *  b40_i;
assign k40 = (a40_r  +  a40_i)  *  (b40_r  +  b40_i);
assign t40_r = p40  -  q40;
assign t40_i = k40  -  p40  -  q40;
assign p41 = a41_r  *  b41_r;
assign q41 = a41_i  *  b41_i;
assign k41 = (a41_r  +  a41_i)  *  (b41_r  +  b41_i);
assign t41_r = p41  -  q41;
assign t41_i = k41  -  p41  -  q41;
assign p42 = a42_r  *  b42_r;
assign q42 = a42_i  *  b42_i;
assign k42 = (a42_r  +  a42_i)  *  (b42_r  +  b42_i);
assign t42_r = p42  -  q42;
assign t42_i = k42  -  p42  -  q42;
assign p43 = a43_r  *  b43_r;
assign q43 = a43_i  *  b43_i;
assign k43 = (a43_r  +  a43_i)  *  (b43_r  +  b43_i);
assign t43_r = p43  -  q43;
assign t43_i = k43  -  p43  -  q43;
assign p44 = a44_r  *  b44_r;
assign q44 = a44_i  *  b44_i;
assign k44 = (a44_r  +  a44_i)  *  (b44_r  +  b44_i);
assign t44_r = p44  -  q44;
assign t44_i = k44  -  p44  -  q44;
assign p45 = a45_r  *  b45_r;
assign q45 = a45_i  *  b45_i;
assign k45 = (a45_r  +  a45_i)  *  (b45_r  +  b45_i);
assign t45_r = p45  -  q45;
assign t45_i = k45  -  p45  -  q45;
assign p46 = a46_r  *  b46_r;
assign q46 = a46_i  *  b46_i;
assign k46 = (a46_r  +  a46_i)  *  (b46_r  +  b46_i);
assign t46_r = p46  -  q46;
assign t46_i = k46  -  p46  -  q46;
assign p47 = a47_r  *  b47_r;
assign q47 = a47_i  *  b47_i;
assign k47 = (a47_r  +  a47_i)  *  (b47_r  +  b47_i);
assign t47_r = p47  -  q47;
assign t47_i = k47  -  p47  -  q47;
//  Construct the result matrix;
//  Number of sub/add operations: 505
assign mt34_imt47_i =  - t34_i - t47_i;
assign mt33_imt47_i =  - t33_i - t47_i;
assign pt38_imt47_i =  + t38_i - t47_i;
assign mt42_imt46_r =  - t42_i - t46_r;
assign pt46_imt44_r =  + t46_i - t44_r;
assign pt07_rmt42_r =  + t07_r - t42_r;
assign mt28_rmt40_i =  - t28_r - t40_i;
assign mt38_imt39_r =  - t38_i - t39_r;
assign mt19_imt35_i =  - t19_i - t35_i;
assign mt23_rmt30_i =  - t23_r - t30_i;
assign pt15_imt30_i =  + t15_i - t30_i;
assign pt47_rmt30_r =  + t47_r - t30_r;
assign pt16_rmt26_r =  + t16_r - t26_r;
assign pt06_imt24_i =  + t06_i - t24_i;
assign pt19_rmt23_i =  + t19_r - t23_i;
assign pt34_rmt21_i =  + t34_r - t21_i;
assign mt07_imt21_r =  - t07_i - t21_r;
assign pt35_imt19_i =  + t35_i - t19_i;
assign pt42_rmt16_i =  + t42_r - t16_i;
assign pt01_rmt16_r =  + t01_r - t16_r;
assign pt30_rmt15_i =  + t30_r - t15_i;
assign pt17_rmt15_r =  + t17_r - t15_r;
assign mt09_rmt14_r =  - t09_r - t14_r;
assign pt43_rmt08_r =  + t43_r - t08_r;
assign pt06_rmt05_r =  + t06_r - t05_r;
assign pt41_rmt04_r =  + t41_r - t04_r;
assign pt41_rmt00_r =  + t41_r - t00_r;
assign pt04_ipt46_i =  + t04_i + t46_i;
assign pt42_ipt46_r =  + t42_i + t46_r;
assign pt00_ipt42_r =  + t00_i + t42_r;
assign pt32_rpt40_i =  + t32_r + t40_i;
assign pt33_ipt37_i =  + t33_i + t37_i;
assign pt21_rpt34_i =  + t21_r + t34_i;
assign pt20_rpt25_r =  + t20_r + t25_r;
assign pt14_ipt18_i =  + t14_i + t18_i;
assign pt08_rpt18_i =  + t08_r + t18_i;
assign pt05_rpt06_r =  + t05_r + t06_r;
assign mt34_rmt47_r =  - t34_r - t47_r;
assign pt29_rmt45_r =  + t29_r - t45_r;
assign pt02_imt35_r =  + t02_i - t35_r;
assign pt27_imt33_r =  + t27_i - t33_r;
assign mt37_rpt27_imt33_r =  - t37_r + pt27_imt33_r;
assign pt14_rmt28_i =  + t14_r - t28_i;
assign pt27_rmt18_i =  + t27_r - t18_i;
assign pt27_rmt18_ipt33_ipt37_i = pt27_rmt18_i + pt33_ipt37_i;
assign mt03_imt12_i =  - t03_i - t12_i;
assign mt06_imt10_i =  - t06_i - t10_i;
assign pt33_rpt37_r =  + t33_r + t37_r;
assign pt20_ipt25_i =  + t20_i + t25_i;
assign mt08_imt47_r =  - t08_i - t47_r;
assign mt08_imt47_rmt03_imt12_i = mt08_imt47_r + mt03_imt12_i;
assign pt16_rmt46_i =  + t16_r - t46_i;
assign pt16_rmt46_imt09_rmt14_r = pt16_rmt46_i + mt09_rmt14_r;
assign mt40_imt46_r =  - t40_i - t46_r;
assign mt13_rmt44_i =  - t13_r - t44_i;
assign pt13_imt44_r =  + t13_i - t44_r;
assign pt23_rpt13_imt44_r =  + t23_r + pt13_imt44_r;
assign pt00_rmt41_r =  + t00_r - t41_r;
assign mt20_ipt00_rmt41_r =  - t20_i + pt00_rmt41_r;
assign mt28_imt20_ipt00_rmt41_r =  - t28_i + mt20_ipt00_rmt41_r;
assign pt30_rmt38_i =  + t30_r - t38_i;
assign mt36_ipt30_rmt38_i =  - t36_i + pt30_rmt38_i;
assign pt05_rmt36_r =  + t05_r - t36_r;
assign pt21_imt33_r =  + t21_i - t33_r;
assign pt26_imt32_i =  + t26_i - t32_i;
assign mt44_rpt26_imt32_i =  - t44_r + pt26_imt32_i;
assign pt45_imt29_i =  + t45_i - t29_i;
assign pt45_imt29_ipt04_ipt46_i = pt45_imt29_i + pt04_ipt46_i;
assign mt04_ipt45_imt29_i =  - t04_i + pt45_imt29_i;
assign mt18_rmt27_i =  - t18_r - t27_i;
assign mt18_rmt27_ipt33_rpt37_r = mt18_rmt27_i + pt33_rpt37_r;
assign pt11_imt27_r =  + t11_i - t27_r;
assign pt32_imt26_i =  + t32_i - t26_i;
assign mt17_rmt22_i =  - t17_r - t22_i;
assign pt06_imt17_rmt22_i =  + t06_i + mt17_rmt22_i;
assign pt22_rmt17_i =  + t22_r - t17_i;
assign pt22_rmt17_ipt06_rmt05_r = pt22_rmt17_i + pt06_rmt05_r;
assign pt11_rmt15_i =  + t11_r - t15_i;
assign mt12_rpt11_rmt15_i =  - t12_r + pt11_rmt15_i;
assign pt08_imt12_rpt11_rmt15_i =  + t08_i + mt12_rpt11_rmt15_i;
assign pt12_imt15_r =  + t12_i - t15_r;
assign pt12_imt15_rpt02_imt35_r = pt12_imt15_r + pt02_imt35_r;
assign pt00_imt14_i =  + t00_i - t14_i;
assign mt41_ipt00_imt14_i =  - t41_i + pt00_imt14_i;
assign mt03_rmt13_i =  - t03_r - t13_i;
assign mt06_rmt10_r =  - t06_r - t10_r;
assign mt06_rmt10_rpt07_rmt42_r = mt06_rmt10_r + pt07_rmt42_r;
assign mt07_rmt06_rmt10_r =  - t07_r + mt06_rmt10_r;
assign pt40_imt09_i =  + t40_i - t09_i;
assign pt36_rmt05_r =  + t36_r - t05_r;
assign mt38_rpt36_rmt05_r =  - t38_r + pt36_rmt05_r;
assign mt34_imt38_rpt36_rmt05_r =  - t34_i + mt38_rpt36_rmt05_r;
assign pt31_imt03_i =  + t31_i - t03_i;
assign pt31_imt03_ipt22_rmt17_i = pt31_imt03_i + pt22_rmt17_i;
assign pt42_imt01_i =  + t42_i - t01_i;
assign mt28_rpt42_imt01_i =  - t28_r + pt42_imt01_i;
assign pt41_imt00_i =  + t41_i - t00_i;
assign mt24_rpt41_imt00_i =  - t24_r + pt41_imt00_i;
assign pt09_imt00_r =  + t09_i - t00_r;
assign pt09_imt00_rmt40_imt46_r = pt09_imt00_r + mt40_imt46_r;
assign pt09_imt00_rmt40_imt46_rpt14_ipt18_i = pt09_imt00_rmt40_imt46_r + pt14_ipt18_i;
assign pt09_imt00_rmt40_imt46_rpt14_rmt28_i = pt09_imt00_rmt40_imt46_r + pt14_rmt28_i;
assign pt01_rmt16_rpt09_imt00_rmt40_imt46_rpt14_rmt28_i = pt01_rmt16_r + pt09_imt00_rmt40_imt46_rpt14_rmt28_i;
assign pt39_ipt47_i =  + t39_i + t47_i;
assign pt13_rpt44_i =  + t13_r + t44_i;
assign mt03_rpt13_rpt44_i =  - t03_r + pt13_rpt44_i;
assign mt29_rmt03_rpt13_rpt44_i =  - t29_r + mt03_rpt13_rpt44_i;
assign mt29_rmt03_rpt13_rpt44_imt18_rmt27_ipt33_rpt37_r = mt29_rmt03_rpt13_rpt44_i + mt18_rmt27_ipt33_rpt37_r;
assign pt29_ipt44_r =  + t29_i + t44_r;
assign pt29_ipt44_rmt03_rmt13_i = pt29_ipt44_r + mt03_rmt13_i;
assign pt20_rpt29_ipt44_rmt03_rmt13_i =  + t20_r + pt29_ipt44_rmt03_rmt13_i;
assign pt15_rpt40_r =  + t15_r + t40_r;
assign pt15_rpt40_rpt11_imt27_r = pt15_rpt40_r + pt11_imt27_r;
assign pt01_ipt28_r =  + t01_i + t28_r;
assign pt01_ipt28_rpt40_imt09_i = pt01_ipt28_r + pt40_imt09_i;
assign pt16_rmt46_ipt01_ipt28_rpt40_imt09_i = pt16_rmt46_i + pt01_ipt28_rpt40_imt09_i;
assign pt00_ipt42_rpt16_rmt46_ipt01_ipt28_rpt40_imt09_i = pt00_ipt42_r + pt16_rmt46_ipt01_ipt28_rpt40_imt09_i;
assign pt17_rpt22_i =  + t17_r + t22_i;
assign pt17_rpt22_imt28_rpt42_imt01_i = pt17_rpt22_i + mt28_rpt42_imt01_i;
assign pt08_rpt16_i =  + t08_r + t16_i;
assign pt12_rpt15_i =  + t12_r + t15_i;
assign pt12_rpt15_ipt32_rpt40_i = pt12_rpt15_i + pt32_rpt40_i;
assign pt12_rpt15_ipt13_imt44_r = pt12_rpt15_i + pt13_imt44_r;
assign mt23_rmt30_ipt12_rpt15_ipt13_imt44_r = mt23_rmt30_i + pt12_rpt15_ipt13_imt44_r;
assign pt43_rmt39_i =  + t43_r - t39_i;
assign pt39_rmt38_r =  + t39_r - t38_r;
assign pt39_rmt38_rpt05_rmt36_r = pt39_rmt38_r + pt05_rmt36_r;
assign pt43_imt34_r =  + t43_i - t34_r;
assign pt43_imt34_rpt21_imt33_r = pt43_imt34_r + pt21_imt33_r;
assign pt03_ipt43_imt34_r =  + t03_i + pt43_imt34_r;
assign pt03_ipt43_imt34_rmt41_ipt00_imt14_i = pt03_ipt43_imt34_r + mt41_ipt00_imt14_i;
assign pvery_long_op_0002 = pt20_rpt29_ipt44_rmt03_rmt13_i + pt03_ipt43_imt34_rmt41_ipt00_imt14_i;
assign pt14_rmt27_i =  + t14_r - t27_i;
assign pt14_rmt27_ipt21_imt33_r = pt14_rmt27_i + pt21_imt33_r;
assign pt14_rmt27_ipt08_rpt16_i = pt14_rmt27_i + pt08_rpt16_i;
assign pt01_ipt28_rpt40_imt09_ipt14_rmt27_ipt08_rpt16_i = pt01_ipt28_rpt40_imt09_i + pt14_rmt27_ipt08_rpt16_i;
assign pt03_imt14_r =  + t03_i - t14_r;
assign pt34_imt11_i =  + t34_i - t11_i;
assign pt03_rpt34_imt11_i =  + t03_r + pt34_imt11_i;
assign pt34_rmt11_r =  + t34_r - t11_r;
assign mt00_ipt34_rmt11_r =  - t00_i + pt34_rmt11_r;
assign mt43_ipt34_rmt11_r =  - t43_i + pt34_rmt11_r;
assign mt14_imt43_ipt34_rmt11_r =  - t14_i + mt43_ipt34_rmt11_r;
assign pvery_long_op_0003 = mt14_imt43_ipt34_rmt11_r + mt29_rmt03_rpt13_rpt44_imt18_rmt27_ipt33_rpt37_r;
assign pt33_rpt47_r =  + t33_r + t47_r;
assign mt02_rpt33_rpt47_r =  - t02_r + pt33_rpt47_r;
assign pt05_ipt18_r =  + t05_i + t18_r;
assign pt28_rpt05_ipt18_r =  + t28_r + pt05_ipt18_r;
assign pt32_rpt05_ipt18_r =  + t32_r + pt05_ipt18_r;
assign pt32_rpt05_ipt18_rpt16_rmt26_r = pt32_rpt05_ipt18_r + pt16_rmt26_r;
assign mt36_ipt30_rmt38_ipt32_rpt05_ipt18_rpt16_rmt26_r = mt36_ipt30_rmt38_i + pt32_rpt05_ipt18_rpt16_rmt26_r;
assign pt09_rmt40_r =  + t09_r - t40_r;
assign pt16_ipt09_rmt40_r =  + t16_i + pt09_rmt40_r;
assign mt33_imt37_i =  - t33_i - t37_i;
assign mt20_rmt33_imt37_i =  - t20_r + mt33_imt37_i;
assign pt29_ipt44_rmt03_rmt13_imt20_rmt33_imt37_i = pt29_ipt44_rmt03_rmt13_i + mt20_rmt33_imt37_i;
assign pt26_rmt32_r =  + t26_r - t32_r;
assign pt26_rmt32_rmt13_rmt44_i = pt26_rmt32_r + mt13_rmt44_i;
assign pt03_rmt31_r =  + t03_r - t31_r;
assign pt24_rmt11_r =  + t24_r - t11_r;
assign mt08_ipt24_rmt11_r =  - t08_i + pt24_rmt11_r;
assign pt37_rpt24_rmt11_r =  + t37_r + pt24_rmt11_r;
assign pt36_imt05_i =  + t36_i - t05_i;
assign pt36_imt05_ipt30_rmt38_i = pt36_imt05_i + pt30_rmt38_i;
assign pt30_ipt36_imt05_i =  + t30_i + pt36_imt05_i;
assign pt30_ipt36_imt05_ipt23_rpt13_imt44_r = pt30_ipt36_imt05_i + pt23_rpt13_imt44_r;
assign pt30_ipt36_imt05_imt28_rpt42_imt01_i = pt30_ipt36_imt05_i + mt28_rpt42_imt01_i;
assign pt28_imt01_r =  + t28_i - t01_r;
assign pt28_imt01_rpt09_rmt40_r = pt28_imt01_r + pt09_rmt40_r;
assign pt28_imt01_rpt09_rmt40_rmt06_imt10_i = pt28_imt01_rpt09_rmt40_r + mt06_imt10_i;
assign pt29_rmt45_rpt28_imt01_rpt09_rmt40_rmt06_imt10_i = pt29_rmt45_r + pt28_imt01_rpt09_rmt40_rmt06_imt10_i;
assign mt44_ipt28_imt01_rpt09_rmt40_r =  - t44_i + pt28_imt01_rpt09_rmt40_r;
assign pt33_ipt47_i =  + t33_i + t47_i;
assign mt02_ipt33_ipt47_i =  - t02_i + pt33_ipt47_i;
assign pt12_imt15_rmt02_ipt33_ipt47_i = pt12_imt15_r + mt02_ipt33_ipt47_i;
assign pt02_rpt33_ipt47_i =  + t02_r + pt33_ipt47_i;
assign pt37_ipt39_i =  + t37_i + t39_i;
assign pt17_ipt18_r =  + t17_i + t18_r;
assign pt17_ipt18_rpt39_rmt38_rpt05_rmt36_r = pt17_ipt18_r + pt39_rmt38_rpt05_rmt36_r;
assign mt22_rpt17_ipt18_r =  - t22_r + pt17_ipt18_r;
assign mt22_rpt17_ipt18_rmt37_rpt27_imt33_r = mt22_rpt17_ipt18_r + mt37_rpt27_imt33_r;
assign pt08_imt43_i =  + t08_i - t43_i;
assign pt18_imt43_r =  + t18_i - t43_r;
assign mt34_ipt18_imt43_r =  - t34_i + pt18_imt43_r;
assign pt03_imt14_rmt34_ipt18_imt43_r = pt03_imt14_r + mt34_ipt18_imt43_r;
assign mt17_rmt22_ipt03_imt14_rmt34_ipt18_imt43_r = mt17_rmt22_i + pt03_imt14_rmt34_ipt18_imt43_r;
assign mt24_rmt37_r =  - t24_r - t37_r;
assign mt31_imt24_rmt37_r =  - t31_i + mt24_rmt37_r;
assign pt11_imt08_r =  + t11_i - t08_r;
assign mt18_ipt11_imt08_r =  - t18_i + pt11_imt08_r;
assign pt11_imt08_rpt43_rmt39_i = pt11_imt08_r + pt43_rmt39_i;
assign pt03_rmt31_rpt11_imt08_rpt43_rmt39_i = pt03_rmt31_r + pt11_imt08_rpt43_rmt39_i;
assign mt34_imt47_ipt03_rmt31_rpt11_imt08_rpt43_rmt39_i = mt34_imt47_i + pt03_rmt31_rpt11_imt08_rpt43_rmt39_i;
assign pt11_rmt39_r =  + t11_r - t39_r;
assign pt27_ipt11_rmt39_r =  + t27_i + pt11_rmt39_r;
assign pt08_rpt16_ipt27_ipt11_rmt39_r = pt08_rpt16_i + pt27_ipt11_rmt39_r;
assign pt42_ipt46_rpt08_rpt16_ipt27_ipt11_rmt39_r = pt42_ipt46_r + pt08_rpt16_ipt27_ipt11_rmt39_r;
assign pt11_rmt39_rpt08_imt43_i = pt11_rmt39_r + pt08_imt43_i;
assign pt34_imt11_ipt11_rmt39_rpt08_imt43_i = pt34_imt11_i + pt11_rmt39_rpt08_imt43_i;
assign pt17_rmt15_rpt34_imt11_ipt11_rmt39_rpt08_imt43_i = pt17_rmt15_r + pt34_imt11_ipt11_rmt39_rpt08_imt43_i;
assign pt28_imt01_rpt11_rmt39_rpt08_imt43_i = pt28_imt01_r + pt11_rmt39_rpt08_imt43_i;
assign pt42_rmt16_ipt28_imt01_rpt11_rmt39_rpt08_imt43_i = pt42_rmt16_i + pt28_imt01_rpt11_rmt39_rpt08_imt43_i;
assign pvery_long_op_0006 = pt09_imt00_rmt40_imt46_rpt14_ipt18_i + pt42_rmt16_ipt28_imt01_rpt11_rmt39_rpt08_imt43_i;
assign mt16_imt18_i =  - t16_i - t18_i;
assign mt16_imt18_ipt21_rpt34_i = mt16_imt18_i + pt21_rpt34_i;
assign pt43_rmt08_rmt16_imt18_ipt21_rpt34_i = pt43_rmt08_r + mt16_imt18_ipt21_rpt34_i;
assign mt24_imt16_imt18_i =  - t24_i + mt16_imt18_i;
assign mt37_imt24_imt16_imt18_i =  - t37_i + mt24_imt16_imt18_i;
assign pt11_imt08_rpt43_rmt39_imt37_imt24_imt16_imt18_i = pt11_imt08_rpt43_rmt39_i + mt37_imt24_imt16_imt18_i;
assign mt16_rmt18_r =  - t16_r - t18_r;
assign mt16_rmt18_rpt08_imt43_i = mt16_rmt18_r + pt08_imt43_i;
assign mt16_rmt18_rmt24_rmt37_r = mt16_rmt18_r + mt24_rmt37_r;
assign mt16_rmt18_rmt24_rmt37_rpt32_imt26_i = mt16_rmt18_rmt24_rmt37_r + pt32_imt26_i;
assign pt38_imt47_imt16_rmt18_rmt24_rmt37_rpt32_imt26_i = pt38_imt47_i + mt16_rmt18_rmt24_rmt37_rpt32_imt26_i;
assign pt26_rmt32_rmt16_rmt18_rmt24_rmt37_r = pt26_rmt32_r + mt16_rmt18_rmt24_rmt37_r;
assign pt39_ipt47_ipt26_rmt32_rmt16_rmt18_rmt24_rmt37_r = pt39_ipt47_i + pt26_rmt32_rmt16_rmt18_rmt24_rmt37_r;
assign pt39_ipt47_ipt26_rmt32_rmt16_rmt18_rmt24_rmt37_rpt30_rmt15_i = pt39_ipt47_ipt26_rmt32_rmt16_rmt18_rmt24_rmt37_r + pt30_rmt15_i;
assign pvery_long_op_0004 = mt34_imt38_rpt36_rmt05_r + pt39_ipt47_ipt26_rmt32_rmt16_rmt18_rmt24_rmt37_rpt30_rmt15_i;
assign pt24_imt11_i =  + t24_i - t11_i;
assign pt24_imt11_ipt08_rpt18_i = pt24_imt11_i + pt08_rpt18_i;
assign pt24_imt11_ipt37_ipt39_i = pt24_imt11_i + pt37_ipt39_i;
assign pt18_imt43_rpt24_imt11_ipt37_ipt39_i = pt18_imt43_r + pt24_imt11_ipt37_ipt39_i;
assign pt14_ipt27_r =  + t14_i + t27_r;
assign pt14_ipt27_rmt16_imt18_i = pt14_ipt27_r + mt16_imt18_i;
assign mt21_rpt14_ipt27_r =  - t21_r + pt14_ipt27_r;
assign mt24_imt16_imt18_imt21_rpt14_ipt27_r = mt24_imt16_imt18_i + mt21_rpt14_ipt27_r;
assign pt14_ipt27_rpt03_rmt31_r = pt14_ipt27_r + pt03_rmt31_r;
assign pt24_imt11_ipt37_ipt39_ipt14_ipt27_rpt03_rmt31_r = pt24_imt11_ipt37_ipt39_i + pt14_ipt27_rpt03_rmt31_r;
assign pt24_imt11_ipt37_ipt39_ipt14_ipt27_rpt03_rmt31_rpt34_rmt21_i = pt24_imt11_ipt37_ipt39_ipt14_ipt27_rpt03_rmt31_r + pt34_rmt21_i;
assign pvery_long_op_0005 = mt16_rmt18_rpt08_imt43_i + pt24_imt11_ipt37_ipt39_ipt14_ipt27_rpt03_rmt31_rpt34_rmt21_i;
assign mt08_imt47_rpt24_imt11_ipt37_ipt39_ipt14_ipt27_rpt03_rmt31_r = mt08_imt47_r + pt24_imt11_ipt37_ipt39_ipt14_ipt27_rpt03_rmt31_r;
assign pvery_long_op_0001 = pt43_imt34_rpt21_imt33_r + mt08_imt47_rpt24_imt11_ipt37_ipt39_ipt14_ipt27_rpt03_rmt31_r;
assign C1_r =  - t42_r + t47_r - t17_i + mt16_rmt18_r + pt39_rmt38_r + pt36_rmt05_r + pt32_imt26_i + pt15_imt30_i + pt46_imt44_r + pt18_imt43_rpt24_imt11_ipt37_ipt39_i + mt00_ipt34_rmt11_r + pt01_ipt28_rpt40_imt09_ipt14_rmt27_ipt08_rpt16_i;
assign C1_i = ( + t00_r + pt39_ipt47_ipt26_rmt32_rmt16_rmt18_rmt24_rmt37_r + mt42_imt46_r + pt14_ipt27_rmt16_imt18_i + pt36_imt05_ipt30_rmt38_i + mt44_ipt28_imt01_rpt09_rmt40_r + pt17_rmt15_rpt34_imt11_ipt11_rmt39_rpt08_imt43_i);
assign C2_r =  - t24_i - t32_i + pt02_rpt33_ipt47_i + mt12_rpt11_rmt15_i + pt05_rpt06_r + mt28_rmt40_i + mt18_ipt11_imt08_r + mt22_rpt17_ipt18_rmt37_rpt27_imt33_r + pvery_long_op_0002;
assign C2_i = ( + t02_i + t20_i - t33_r + pt24_rmt11_r + mt33_imt37_i + mt34_ipt18_imt43_r + pt14_rmt28_i + pt32_rpt05_ipt18_r + mt29_rmt03_rpt13_rpt44_i + pt41_rmt00_r + pt06_imt17_rmt22_i + pt15_rpt40_rpt11_imt27_r + mt08_imt47_rmt03_imt12_i);
assign C3_r =  - t20_r - t25_r + pt32_imt26_i + mt02_rpt33_rpt47_r + mt38_rpt36_rmt05_r + mt19_imt35_i + mt23_rmt30_ipt12_rpt15_ipt13_imt44_r + pvery_long_op_0005;
assign C3_i = ( - t20_i - t25_i + t35_r + pt03_imt14_r + pt27_ipt11_rmt39_r + pt19_rmt23_i + mt31_imt24_rmt37_r + pt36_imt05_ipt30_rmt38_i + pt26_rmt32_rmt13_rmt44_i + pt12_imt15_rmt02_ipt33_ipt47_i + pt43_rmt08_rmt16_imt18_ipt21_rpt34_i);
assign C4_r =  - t21_r + t39_r + mt34_ipt18_imt43_r + pt41_imt00_i + pt37_rpt24_rmt11_r + mt33_imt47_i + pt31_imt03_ipt22_rmt17_i + pt45_imt29_ipt04_ipt46_i + mt06_rmt10_rpt07_rmt42_r + pt01_ipt28_rpt40_imt09_ipt14_rmt27_ipt08_rpt16_i;
assign C4_i = ( + t07_i - t04_r + pt33_rpt47_r + pt17_rpt22_i + pt00_rmt41_r + mt42_imt46_r + pt29_rmt45_rpt28_imt01_rpt09_rmt40_rmt06_imt10_i + pvery_long_op_0005);
assign C5_r =  - t17_i + t44_i + pt11_imt27_r + pvery_long_op_0004 + pvery_long_op_0006;
assign C5_i = ( + t17_r - t47_r + pt15_rpt40_r + mt18_rmt27_i + mt38_imt39_r + mt00_ipt34_rmt11_r + mt44_rpt26_imt32_i + pt30_ipt36_imt05_imt28_rpt42_imt01_i + pt16_rmt46_imt09_rmt14_r + pt11_imt08_rpt43_rmt39_imt37_imt24_imt16_imt18_i);
assign C6_r =  + t29_r + t43_r + pt03_imt14_r + mt13_rmt44_i + mt22_rpt17_ipt18_r + mt02_rpt33_rpt47_r + pt05_rpt06_r + mt08_ipt24_rmt11_r + pt03_rpt34_imt11_i + mt28_imt20_ipt00_rmt41_r + pt27_rmt18_ipt33_ipt37_i + pt12_rpt15_ipt32_rpt40_i;
assign C6_i = ( + t32_i + t11_r - t40_r + mt37_rpt27_imt33_r + pt28_rpt05_ipt18_r + pt06_imt17_rmt22_i + pt12_imt15_rmt02_ipt33_ipt47_i + pt24_imt11_ipt08_rpt18_i + pvery_long_op_0002);
assign C7_r =  + t02_r - t12_r - t31_r - t43_r - t23_i + pt11_imt08_r + mt03_rpt13_rpt44_i + pt20_ipt25_i + pt35_imt19_i + pt14_rmt27_ipt21_imt33_r + pvery_long_op_0004;
assign C7_i = ( - t31_i + t15_r + t19_r - t25_r + pt26_imt32_i + pt02_imt35_r + mt20_rmt33_imt37_i + mt43_ipt34_rmt11_r + mt38_imt39_r + mt08_imt47_rmt03_imt12_i + mt24_imt16_imt18_imt21_rpt14_ipt27_r + pt30_ipt36_imt05_ipt23_rpt13_imt44_r);
assign C8_r =  + t21_r - t27_r - t29_r + t45_r + t24_i + mt33_imt37_i + mt34_rmt47_r + pt41_rmt04_r + mt07_rmt06_rmt10_r + pt31_imt03_ipt22_rmt17_i + pvery_long_op_0006;
assign C8_i = ( - t07_i + t21_i + t40_r + mt06_imt10_i + mt18_rmt27_ipt33_rpt37_r + mt24_rpt41_imt00_i + mt04_ipt45_imt29_i + pt17_rpt22_imt28_rpt42_imt01_i + pt16_rmt46_imt09_rmt14_r + mt34_imt47_ipt03_rmt31_rpt11_imt08_rpt43_rmt39_i);
assign C9_r =  - t39_i + pt00_imt14_i + pt27_rmt18_i + pt46_imt44_r + pt16_ipt09_rmt40_r + pt30_ipt36_imt05_imt28_rpt42_imt01_i + pt17_rmt15_rpt34_imt11_ipt11_rmt39_rpt08_imt43_i + pt38_imt47_imt16_rmt18_rmt24_rmt37_rpt32_imt26_i;
assign C9_i = ( + t27_i - t44_i - t34_r - t42_r + pt26_rmt32_r + pt11_rmt15_i + pt47_rmt30_r + pt17_ipt18_rpt39_rmt38_rpt05_rmt36_r + pt11_imt08_rpt43_rmt39_imt37_imt24_imt16_imt18_i + pt01_rmt16_rpt09_imt00_rmt40_imt46_rpt14_rmt28_i);
assign C10_r =  - t02_i - t06_i - t12_i - t32_i + pt33_rpt47_r + pt41_imt00_i + mt08_ipt24_rmt11_r + pt28_rpt05_ipt18_r + pt15_rpt40_rpt11_imt27_r + pt29_ipt44_rmt03_rmt13_imt20_rmt33_imt37_i + mt17_rmt22_ipt03_imt14_rmt34_ipt18_imt43_r;
assign C10_i = ( - t03_i + t28_i + pt02_rpt33_ipt47_i + mt20_ipt00_rmt41_r + pt24_imt11_ipt08_rpt18_i + pt12_rpt15_ipt32_rpt40_i + pt22_rmt17_ipt06_rmt05_r + pvery_long_op_0003);
assign C11_r =  + t08_r - t19_r + t31_r + pt43_rmt39_i + pt20_rpt25_r + pt03_rpt34_imt11_i + pt14_rmt27_ipt21_imt33_r + pt12_imt15_rpt02_imt35_r + pt38_imt47_imt16_rmt18_rmt24_rmt37_rpt32_imt26_i + pt30_ipt36_imt05_ipt23_rpt13_imt44_r;
assign C11_i = ( + t23_i + t31_i - t02_r + mt33_imt37_i + pt39_rmt38_rpt05_rmt36_r + pt20_ipt25_i + pt03_ipt43_imt34_r + pt47_rmt30_r + mt19_imt35_i + pt26_rmt32_rmt13_rmt44_i + pt08_imt12_rpt11_rmt15_i + mt24_imt16_imt18_imt21_rpt14_ipt27_r);
assign C12_r =  + t10_i + pt11_rmt39_rpt08_imt43_i + pt31_imt03_i + mt34_rmt47_r + mt41_ipt00_imt14_i + mt07_imt21_r + pt06_imt24_i + pt16_ipt09_rmt40_r + pt17_rpt22_imt28_rpt42_imt01_i + pt27_rmt18_ipt33_ipt37_i + pt45_imt29_ipt04_ipt46_i;
assign C12_i = ( - t21_i + t24_r + pt29_rmt45_r + pt41_rmt04_r + mt22_rpt17_ipt18_rmt37_rpt27_imt33_r + mt06_rmt10_rpt07_rmt42_r + mt34_imt47_ipt03_rmt31_rpt11_imt08_rpt43_rmt39_i + pt01_rmt16_rpt09_imt00_rmt40_imt46_rpt14_rmt28_i);
assign C13_r =  - t08_i + t43_i + pt34_imt11_i + pt39_ipt47_i + pt37_rpt24_rmt11_r + pt15_imt30_i + pt14_ipt27_rmt16_imt18_i + mt44_rpt26_imt32_i + pt17_ipt18_rpt39_rmt38_rpt05_rmt36_r + pt00_ipt42_rpt16_rmt46_ipt01_ipt28_rpt40_imt09_i;
assign C13_i = ( - t00_r - t14_r - t15_r - t17_r + mt34_rmt47_r + pt18_imt43_rpt24_imt11_ipt37_ipt39_i + mt44_ipt28_imt01_rpt09_rmt40_r + pt42_ipt46_rpt08_rpt16_ipt27_ipt11_rmt39_r + mt36_ipt30_rmt38_ipt32_rpt05_ipt18_rpt16_rmt26_r);
assign C14_r =  - t18_r + t32_i + pt11_imt27_r + mt02_rpt33_rpt47_r + pt03_imt14_rmt34_ipt18_imt43_r + mt28_rmt40_i + mt24_rpt41_imt00_i + pt29_ipt44_rmt03_rmt13_imt20_rmt33_imt37_i + pt08_imt12_rpt11_rmt15_i + pt22_rmt17_ipt06_rmt05_r;
assign C14_i = ( - t05_i - t32_r + pt17_rpt22_i + pt15_rpt40_r + mt03_imt12_i + mt02_ipt33_ipt47_i + pt06_imt24_i + mt18_ipt11_imt08_r + mt28_imt20_ipt00_rmt41_r + pvery_long_op_0003);
assign C15_r =  + pt14_rmt27_i + pt31_imt03_i + pt26_imt32_i + pt39_rmt38_rpt05_rmt36_r + pt02_rpt33_ipt47_i + pt37_rpt24_rmt11_r + pt20_rpt25_r + pt35_imt19_i + mt23_rmt30_ipt12_rpt15_ipt13_imt44_r + pt43_rmt08_rmt16_imt18_ipt21_rpt34_i;
assign C15_i = ( + mt13_rmt44_i + pt20_ipt25_i + pt19_rmt23_i + pt12_imt15_rpt02_imt35_r + pvery_long_op_0001 + mt36_ipt30_rmt38_ipt32_rpt05_ipt18_rpt16_rmt26_r);
assign C16_r =  - t41_i + mt22_rpt17_ipt18_r + mt07_rmt06_rmt10_r + mt04_ipt45_imt29_i + pt00_ipt42_rpt16_rmt46_ipt01_ipt28_rpt40_imt09_i + pvery_long_op_0001;
assign C16_i = ( + t04_r + pt41_rmt00_r + mt07_imt21_r + mt33_imt47_i + mt31_imt24_rmt37_r + mt17_rmt22_ipt03_imt14_rmt34_ipt18_imt43_r + pt29_rmt45_rpt28_imt01_rpt09_rmt40_rmt06_imt10_i + pt42_ipt46_rpt08_rpt16_ipt27_ipt11_rmt39_r);


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

