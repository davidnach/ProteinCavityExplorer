/****************************************************************************
 * VASP: Volumetric Analysis of Surface Properties
 * Copyright (c) 2014 Brian Y. Chen
 * All rights reserved.
 *
 * This file is part of VASP
 * 
 * VASP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * VASP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with VASP.  If not, see <http://www.gnu.org/licenses/>.  
 *
 * File: CubeTable.cpp
 *       Contains data encoding triangles for marching cubes
 *
 * Written by 
 *       Brian Y. Chen <chen@lehigh.edu>
 *
 * WWW URL: http://cse.lehigh.edu/~chen/
 * Email: chen@lehigh.edu
 * Documentation can be found here: 
 * http://www.ploscompbiol.org/article/info:doi/10.1371/journal.pcbi.1000881
 *
 ***************************************************************************/

#include "CubeTable.h"


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


CubeTable::CubeTable()
{
	int i = 0;
	size  = 256;
	sizes = new int[size];
	
	sizes[0] = 0;		sizes[1] = 1;		sizes[2] = 1;		sizes[3] = 2;		sizes[4] = 1;
	sizes[5] = 2;		sizes[6] = 4;		sizes[7] = 3;		sizes[8] = 1;		sizes[9] = 4;
	sizes[10] = 2;		sizes[11] = 3;		sizes[12] = 2;		sizes[13] = 3;		sizes[14] = 3;
	sizes[15] = 2;		sizes[16] = 1;		sizes[17] = 2;		sizes[18] = 4;		sizes[19] = 3;
	sizes[20] = 4;		sizes[21] = 3;		sizes[22] = 5;		sizes[23] = 4;		sizes[24] = 2;
	sizes[25] = 5;		sizes[26] = 5;		sizes[27] = 4;		sizes[28] = 5;		sizes[29] = 4;
	sizes[30] = 4;		sizes[31] = 3;		sizes[32] = 1;		sizes[33] = 4;		sizes[34] = 2;
	sizes[35] = 3;		sizes[36] = 2;		sizes[37] = 5;		sizes[38] = 5;		sizes[39] = 4;
	sizes[40] = 4;		sizes[41] = 5;		sizes[42] = 3;		sizes[43] = 4;		sizes[44] = 5;
	sizes[45] = 4;		sizes[46] = 4;		sizes[47] = 3;		sizes[48] = 2;		sizes[49] = 3;
	sizes[50] = 3;		sizes[51] = 2;		sizes[52] = 5;		sizes[53] = 4;		sizes[54] = 4;
	sizes[55] = 3;		sizes[56] = 5;		sizes[57] = 4;		sizes[58] = 4;		sizes[59] = 3;
	sizes[60] = 4;		sizes[61] = 3;		sizes[62] = 3;		sizes[63] = 2;		sizes[64] = 1;
	sizes[65] = 4;		sizes[66] = 2;		sizes[67] = 5;		sizes[68] = 2;		sizes[69] = 3;
	sizes[70] = 5;		sizes[71] = 4;		sizes[72] = 4;		sizes[73] = 5;		sizes[74] = 5;
	sizes[75] = 4;		sizes[76] = 3;		sizes[77] = 4;		sizes[78] = 4;		sizes[79] = 3;
	sizes[80] = 2;		sizes[81] = 3;		sizes[82] = 5;		sizes[83] = 4;		sizes[84] = 3;
	sizes[85] = 2;		sizes[86] = 4;		sizes[87] = 3;		sizes[88] = 5;		sizes[89] = 4;
	sizes[90] = 4;		sizes[91] = 3;		sizes[92] = 4;		sizes[93] = 3;		sizes[94] = 3;
	sizes[95] = 2;		sizes[96] = 4;		sizes[97] = 5;		sizes[98] = 5;		sizes[99] = 4;
	sizes[100] = 5;	sizes[101] = 4;	sizes[102] = 4;	sizes[103] = 3;	sizes[104] = 5;
	sizes[105] = 4;	sizes[106] = 4;	sizes[107] = 3;	sizes[108] = 4;	sizes[109] = 3;
	sizes[110] = 3;	sizes[111] = 2;	sizes[112] = 3;	sizes[113] = 4;	sizes[114] = 4;
	sizes[115] = 3;	sizes[116] = 4;	sizes[117] = 3;	sizes[118] = 3;	sizes[119] = 2;
	sizes[120] = 4;	sizes[121] = 3;	sizes[122] = 3;	sizes[123] = 2;	sizes[124] = 3;
	sizes[125] = 2;	sizes[126] = 2;	sizes[127] = 1;	sizes[128] = 1;	sizes[129] = 2;
	sizes[130] = 4;	sizes[131] = 5;	sizes[132] = 4;	sizes[133] = 5;	sizes[134] = 5;
	sizes[135] = 4;	sizes[136] = 2;	sizes[137] = 5;	sizes[138] = 3;	sizes[139] = 4;
	sizes[140] = 3;	sizes[141] = 4;	sizes[142] = 4;	sizes[143] = 3;	sizes[144] = 4;
	sizes[145] = 5;	sizes[146] = 5;	sizes[147] = 4;	sizes[148] = 5;	sizes[149] = 4;
	sizes[150] = 4;	sizes[151] = 3;	sizes[152] = 5;	sizes[153] = 4;	sizes[154] = 4;
	sizes[155] = 3;	sizes[156] = 4;	sizes[157] = 3;	sizes[158] = 3;	sizes[159] = 2;
	sizes[160] = 2;	sizes[161] = 5;	sizes[162] = 3;	sizes[163] = 4;	sizes[164] = 5;
	sizes[165] = 4;	sizes[166] = 4;	sizes[167] = 3;	sizes[168] = 3;	sizes[169] = 4;
	sizes[170] = 2;	sizes[171] = 3;	sizes[172] = 4;	sizes[173] = 3;	sizes[174] = 3;
	sizes[175] = 2;	sizes[176] = 3;	sizes[177] = 4;	sizes[178] = 4;	sizes[179] = 3;
	sizes[180] = 4;	sizes[181] = 3;	sizes[182] = 3;	sizes[183] = 2;	sizes[184] = 4;
	sizes[185] = 3;	sizes[186] = 3;	sizes[187] = 2;	sizes[188] = 3;	sizes[189] = 2;
	sizes[190] = 2;	sizes[191] = 1;	sizes[192] = 2;	sizes[193] = 5;	sizes[194] = 5;
	sizes[195] = 4;	sizes[196] = 3;	sizes[197] = 4;	sizes[198] = 4;	sizes[199] = 3;
	sizes[200] = 3;	sizes[201] = 4;	sizes[202] = 4;	sizes[203] = 3;	sizes[204] = 2;
	sizes[205] = 3;	sizes[206] = 3;	sizes[207] = 2;	sizes[208] = 3;	sizes[209] = 4;
	sizes[210] = 4;	sizes[211] = 3;	sizes[212] = 4;	sizes[213] = 3;	sizes[214] = 3;
	sizes[215] = 2;	sizes[216] = 4;	sizes[217] = 3;	sizes[218] = 3;	sizes[219] = 2;
	sizes[220] = 3;	sizes[221] = 2;	sizes[222] = 2;	sizes[223] = 1;	sizes[224] = 3;
	sizes[225] = 4;	sizes[226] = 4;	sizes[227] = 3;	sizes[228] = 4;	sizes[229] = 3;
	sizes[230] = 3;	sizes[231] = 2;	sizes[232] = 4;	sizes[233] = 3;	sizes[234] = 3;
	sizes[235] = 2;	sizes[236] = 3;	sizes[237] = 2;	sizes[238] = 2;	sizes[239] = 1;
	sizes[240] = 2;	sizes[241] = 3;	sizes[242] = 3;	sizes[243] = 2;	sizes[244] = 3;
	sizes[245] = 2;	sizes[246] = 2;	sizes[247] = 1;	sizes[248] = 3;	sizes[249] = 2;
	sizes[250] = 2;	sizes[251] = 1;	sizes[252] = 2;	sizes[253] = 1;	sizes[254] = 1;
	sizes[255] = 0;
		
	polys = new int*[size];
	for(i = 0; i<size; i++){
		if(sizes[i]==0){ polys[i] = new int[1]; }
		else{ polys[i] = new int[3*sizes[i]]; }
	}
//	printf("polys[3][3] = %i, polys[3][4] = %i, polys[3][5] = %i\n", polys[3][3], polys[3][4], polys[3][5]);

	polys[1][0] = 0;				polys[1][1] = 1;				polys[1][2] = 4;
	polys[2][0] = 2;				polys[2][1] = 0;				polys[2][2] = 5;
	polys[3][0] = 2;				polys[3][1] = 1;				polys[3][2] = 4;				polys[3][3] = 4;				polys[3][4] = 5;				polys[3][5] = 2;
	polys[4][0] = 1;				polys[4][1] = 3;				polys[4][2] = 6;
	polys[5][0] = 0;				polys[5][1] = 3;				polys[5][2] = 6;				polys[5][3] = 6;				polys[5][4] = 4;				polys[5][5] = 0;
	polys[6][0] = 2;				polys[6][1] = 3;				polys[6][2] = 6;				polys[6][3] = 6;				polys[6][4] = 1;				polys[6][5] = 0;				polys[6][6] = 0;				polys[6][7] = 5;				polys[6][8] = 2;				polys[6][9] = 2;				polys[6][10] = 6;				polys[6][11] = 0;
	polys[7][0] = 2;				polys[7][1] = 3;				polys[7][2] = 6;				polys[7][3] = 6;				polys[7][4] = 4;				polys[7][5] = 5;				polys[7][6] = 5;				polys[7][7] = 2;				polys[7][8] = 6;
	polys[8][0] = 3;				polys[8][1] = 2;				polys[8][2] = 7;
	polys[9][0] = 0;				polys[9][1] = 2;				polys[9][2] = 7;				polys[9][3] = 7;				polys[9][4] = 3;				polys[9][5] = 1;				polys[9][6] = 1;				polys[9][7] = 4;				polys[9][8] = 0;				polys[9][9] = 0;				polys[9][10] = 7;				polys[9][11] = 1;
	polys[10][0] = 3;				polys[10][1] = 0;				polys[10][2] = 5;				polys[10][3] = 5;				polys[10][4] = 7;				polys[10][5] = 3;
	polys[11][0] = 3;				polys[11][1] = 1;				polys[11][2] = 4;				polys[11][3] = 4;				polys[11][4] = 5;				polys[11][5] = 7;				polys[11][6] = 7;				polys[11][7] = 3;				polys[11][8] = 4;
	polys[12][0] = 1;				polys[12][1] = 2;				polys[12][2] = 7;				polys[12][3] = 7;				polys[12][4] = 6;				polys[12][5] = 1;
	polys[13][0] = 0;				polys[13][1] = 2;				polys[13][2] = 7;				polys[13][3] = 7;				polys[13][4] = 6;				polys[13][5] = 4;				polys[13][6] = 4;				polys[13][7] = 0;				polys[13][8] = 7;
	polys[14][0] = 1;				polys[14][1] = 0;				polys[14][2] = 5;				polys[14][3] = 5;				polys[14][4] = 7;				polys[14][5] = 6;				polys[14][6] = 6;				polys[14][7] = 1;				polys[14][8] = 5;
	polys[15][0] = 4;				polys[15][1] = 5;				polys[15][2] = 7;				polys[15][3] = 7;				polys[15][4] = 6;				polys[15][5] = 4;
	polys[16][0] = 9;				polys[16][1] = 8;				polys[16][2] = 4;
	polys[17][0] = 0;				polys[17][1] = 1;				polys[17][2] = 9;				polys[17][3] = 9;				polys[17][4] = 8;				polys[17][5] = 0;
	polys[18][0] = 2;				polys[18][1] = 0;				polys[18][2] = 4;				polys[18][3] = 4;				polys[18][4] = 9;				polys[18][5] = 8;				polys[18][6] = 8;				polys[18][7] = 5;				polys[18][8] = 2;				polys[18][9] = 2;				polys[18][10] = 4;				polys[18][11] = 8;
	polys[19][0] = 2;				polys[19][1] = 1;				polys[19][2] = 9;				polys[19][3] = 9;				polys[19][4] = 8;				polys[19][5] = 5;				polys[19][6] = 5;				polys[19][7] = 2;				polys[19][8] = 9;
	polys[20][0] = 1;				polys[20][1] = 3;				polys[20][2] = 6;				polys[20][3] = 6;				polys[20][4] = 9;				polys[20][5] = 8;				polys[20][6] = 8;				polys[20][7] = 4;				polys[20][8] = 1;				polys[20][9] = 1;				polys[20][10] = 6;				polys[20][11] = 8;
	polys[21][0] = 0;				polys[21][1] = 3;				polys[21][2] = 6;				polys[21][3] = 6;				polys[21][4] = 9;				polys[21][5] = 8;				polys[21][6] = 8;				polys[21][7] = 0;				polys[21][8] = 6;
	polys[22][0] = 2;				polys[22][1] = 3;				polys[22][2] = 6;				polys[22][3] = 6;				polys[22][4] = 9;				polys[22][5] = 8;				polys[22][6] = 8;				polys[22][7] = 5;				polys[22][8] = 2;				polys[22][9] = 2;				polys[22][10] = 6;				polys[22][11] = 8;				polys[22][12] = 1;				polys[22][13] = 0;				polys[22][14] = 4;
	polys[23][0] = 2;				polys[23][1] = 3;				polys[23][2] = 6;				polys[23][3] = 6;				polys[23][4] = 9;				polys[23][5] = 8;				polys[23][6] = 8;				polys[23][7] = 5;				polys[23][8] = 2;				polys[23][9] = 2;				polys[23][10] = 6;				polys[23][11] = 8;
	polys[24][0] = 3;				polys[24][1] = 2;				polys[24][2] = 7;				polys[24][3] = 9;				polys[24][4] = 8;				polys[24][5] = 4;
	polys[25][0] = 0;				polys[25][1] = 2;				polys[25][2] = 7;				polys[25][3] = 7;				polys[25][4] = 3;				polys[25][5] = 1;				polys[25][6] = 1;				polys[25][7] = 9;				polys[25][8] = 8;				polys[25][9] = 8;				polys[25][10] = 0;				polys[25][11] = 7;				polys[25][12] = 7;				polys[25][13] = 1;				polys[25][14] = 8;
	polys[26][0] = 3;				polys[26][1] = 0;				polys[26][2] = 4;				polys[26][3] = 4;				polys[26][4] = 9;				polys[26][5] = 8;				polys[26][6] = 8;				polys[26][7] = 5;				polys[26][8] = 7;				polys[26][9] = 7;				polys[26][10] = 3;				polys[26][11] = 4;				polys[26][12] = 4;				polys[26][13] = 8;				polys[26][14] = 7;
	polys[27][0] = 3;				polys[27][1] = 1;				polys[27][2] = 9;				polys[27][3] = 9;				polys[27][4] = 8;				polys[27][5] = 5;				polys[27][6] = 5;				polys[27][7] = 7;				polys[27][8] = 3;				polys[27][9] = 3;				polys[27][10] = 9;				polys[27][11] = 5;
	polys[28][0] = 1;				polys[28][1] = 2;				polys[28][2] = 7;				polys[28][3] = 7;				polys[28][4] = 6;				polys[28][5] = 9;				polys[28][6] = 9;				polys[28][7] = 8;				polys[28][8] = 4;				polys[28][9] = 4;				polys[28][10] = 1;				polys[28][11] = 7;				polys[28][12] = 7;				polys[28][13] = 9;				polys[28][14] = 4;
	polys[29][0] = 0;				polys[29][1] = 2;				polys[29][2] = 7;				polys[29][3] = 7;				polys[29][4] = 6;				polys[29][5] = 9;				polys[29][6] = 9;				polys[29][7] = 8;				polys[29][8] = 0;				polys[29][9] = 0;				polys[29][10] = 7;				polys[29][11] = 9;
	polys[30][0] = 1;				polys[30][1] = 0;				polys[30][2] = 4;				polys[30][3] = 9;				polys[30][4] = 8;				polys[30][5] = 5;				polys[30][6] = 5;				polys[30][7] = 7;				polys[30][8] = 6;				polys[30][9] = 6;				polys[30][10] = 9;				polys[30][11] = 5;
	polys[31][0] = 9;				polys[31][1] = 8;				polys[31][2] = 5;				polys[31][3] = 5;				polys[31][4] = 7;				polys[31][5] = 6;				polys[31][6] = 6;				polys[31][7] = 9;				polys[31][8] = 5;
	polys[32][0] = 8;				polys[32][1] = 10;				polys[32][2] = 5;
	polys[33][0] = 0;				polys[33][1] = 1;				polys[33][2] = 4;				polys[33][3] = 4;				polys[33][4] = 8;				polys[33][5] = 10;				polys[33][6] = 10;				polys[33][7] = 5;				polys[33][8] = 0;				polys[33][9] = 0;				polys[33][10] = 4;				polys[33][11] = 10;
	polys[34][0] = 2;				polys[34][1] = 0;				polys[34][2] = 8;				polys[34][3] = 8;				polys[34][4] = 10;				polys[34][5] = 2;
	polys[35][0] = 2;				polys[35][1] = 1;				polys[35][2] = 4;				polys[35][3] = 4;				polys[35][4] = 8;				polys[35][5] = 10;				polys[35][6] = 10;				polys[35][7] = 2;				polys[35][8] = 4;
	polys[36][0] = 1;				polys[36][1] = 3;				polys[36][2] = 6;				polys[36][3] = 8;				polys[36][4] = 10;				polys[36][5] = 5;
	polys[37][0] = 0;				polys[37][1] = 3;				polys[37][2] = 6;				polys[37][3] = 6;				polys[37][4] = 4;				polys[37][5] = 8;				polys[37][6] = 8;				polys[37][7] = 10;				polys[37][8] = 5;				polys[37][9] = 5;				polys[37][10] = 0;				polys[37][11] = 6;				polys[37][12] = 6;				polys[37][13] = 8;				polys[37][14] = 5;
	polys[38][0] = 2;				polys[38][1] = 3;				polys[38][2] = 6;				polys[38][3] = 6;				polys[38][4] = 1;				polys[38][5] = 0;				polys[38][6] = 0;				polys[38][7] = 8;				polys[38][8] = 10;				polys[38][9] = 10;				polys[38][10] = 2;				polys[38][11] = 6;				polys[38][12] = 6;				polys[38][13] = 0;				polys[38][14] = 10;
	polys[39][0] = 2;				polys[39][1] = 3;				polys[39][2] = 6;				polys[39][3] = 6;				polys[39][4] = 4;				polys[39][5] = 8;				polys[39][6] = 8;				polys[39][7] = 10;				polys[39][8] = 2;				polys[39][9] = 2;				polys[39][10] = 6;				polys[39][11] = 8;
	polys[40][0] = 3;				polys[40][1] = 2;				polys[40][2] = 5;				polys[40][3] = 5;				polys[40][4] = 8;				polys[40][5] = 10;				polys[40][6] = 10;				polys[40][7] = 7;				polys[40][8] = 3;				polys[40][9] = 3;				polys[40][10] = 5;				polys[40][11] = 10;
	polys[41][0] = 0;				polys[41][1] = 2;				polys[41][2] = 5;				polys[41][3] = 3;				polys[41][4] = 1;				polys[41][5] = 4;				polys[41][6] = 4;				polys[41][7] = 8;				polys[41][8] = 10;				polys[41][9] = 10;				polys[41][10] = 7;				polys[41][11] = 3;				polys[41][12] = 3;				polys[41][13] = 4;				polys[41][14] = 10;
	polys[42][0] = 3;				polys[42][1] = 0;				polys[42][2] = 8;				polys[42][3] = 8;				polys[42][4] = 10;				polys[42][5] = 7;				polys[42][6] = 7;				polys[42][7] = 3;				polys[42][8] = 8;
	polys[43][0] = 3;				polys[43][1] = 1;				polys[43][2] = 4;				polys[43][3] = 4;				polys[43][4] = 8;				polys[43][5] = 10;				polys[43][6] = 10;				polys[43][7] = 7;				polys[43][8] = 3;				polys[43][9] = 3;				polys[43][10] = 4;				polys[43][11] = 10;
	polys[44][0] = 1;				polys[44][1] = 2;				polys[44][2] = 5;				polys[44][3] = 5;				polys[44][4] = 8;				polys[44][5] = 10;				polys[44][6] = 10;				polys[44][7] = 7;				polys[44][8] = 6;				polys[44][9] = 6;				polys[44][10] = 1;				polys[44][11] = 5;				polys[44][12] = 5;				polys[44][13] = 10;				polys[44][14] = 6;
	polys[45][0] = 0;				polys[45][1] = 2;				polys[45][2] = 5;				polys[45][3] = 8;				polys[45][4] = 10;				polys[45][5] = 7;				polys[45][6] = 7;				polys[45][7] = 6;				polys[45][8] = 4;				polys[45][9] = 4;				polys[45][10] = 8;				polys[45][11] = 7;
	polys[46][0] = 1;				polys[46][1] = 0;				polys[46][2] = 8;				polys[46][3] = 8;				polys[46][4] = 10;				polys[46][5] = 7;				polys[46][6] = 7;				polys[46][7] = 6;				polys[46][8] = 1;				polys[46][9] = 1;				polys[46][10] = 8;				polys[46][11] = 7;
	polys[47][0] = 8;				polys[47][1] = 10;				polys[47][2] = 7;				polys[47][3] = 7;				polys[47][4] = 6;				polys[47][5] = 4;				polys[47][6] = 4;				polys[47][7] = 8;				polys[47][8] = 7;
	polys[48][0] = 9;				polys[48][1] = 10;				polys[48][2] = 5;				polys[48][3] = 5;				polys[48][4] = 4;				polys[48][5] = 9;
	polys[49][0] = 0;				polys[49][1] = 1;				polys[49][2] = 9;				polys[49][3] = 9;				polys[49][4] = 10;				polys[49][5] = 5;				polys[49][6] = 5;				polys[49][7] = 0;				polys[49][8] = 9;
	polys[50][0] = 2;				polys[50][1] = 0;				polys[50][2] = 4;				polys[50][3] = 4;				polys[50][4] = 9;				polys[50][5] = 10;				polys[50][6] = 10;				polys[50][7] = 2;				polys[50][8] = 4;
	polys[51][0] = 2;				polys[51][1] = 1;				polys[51][2] = 9;				polys[51][3] = 9;				polys[51][4] = 10;				polys[51][5] = 2;
	polys[52][0] = 1;				polys[52][1] = 3;				polys[52][2] = 6;				polys[52][3] = 6;				polys[52][4] = 9;				polys[52][5] = 10;				polys[52][6] = 10;				polys[52][7] = 5;				polys[52][8] = 4;				polys[52][9] = 4;				polys[52][10] = 1;				polys[52][11] = 6;				polys[52][12] = 6;				polys[52][13] = 10;				polys[52][14] = 4;
	polys[53][0] = 0;				polys[53][1] = 3;				polys[53][2] = 6;				polys[53][3] = 6;				polys[53][4] = 9;				polys[53][5] = 10;				polys[53][6] = 10;				polys[53][7] = 5;				polys[53][8] = 0;				polys[53][9] = 0;				polys[53][10] = 6;				polys[53][11] = 10;
	polys[54][0] = 2;				polys[54][1] = 3;				polys[54][2] = 6;				polys[54][3] = 6;				polys[54][4] = 9;				polys[54][5] = 10;				polys[54][6] = 10;				polys[54][7] = 2;				polys[54][8] = 6;				polys[54][9] = 1;				polys[54][10] = 0;				polys[54][11] = 4;
	polys[55][0] = 2;				polys[55][1] = 3;				polys[55][2] = 6;				polys[55][3] = 6;				polys[55][4] = 9;				polys[55][5] = 10;				polys[55][6] = 10;				polys[55][7] = 2;				polys[55][8] = 6;
	polys[56][0] = 3;				polys[56][1] = 2;				polys[56][2] = 5;				polys[56][3] = 5;				polys[56][4] = 4;				polys[56][5] = 9;				polys[56][6] = 9;				polys[56][7] = 10;				polys[56][8] = 7;				polys[56][9] = 7;				polys[56][10] = 3;				polys[56][11] = 5;				polys[56][12] = 5;				polys[56][13] = 9;				polys[56][14] = 7;
	polys[57][0] = 0;				polys[57][1] = 2;				polys[57][2] = 5;				polys[57][3] = 3;				polys[57][4] = 1;				polys[57][5] = 9;				polys[57][6] = 9;				polys[57][7] = 10;				polys[57][8] = 7;				polys[57][9] = 7;				polys[57][10] = 3;				polys[57][11] = 9;
	polys[58][0] = 3;				polys[58][1] = 0;				polys[58][2] = 4;				polys[58][3] = 4;				polys[58][4] = 9;				polys[58][5] = 10;				polys[58][6] = 10;				polys[58][7] = 7;				polys[58][8] = 3;				polys[58][9] = 3;				polys[58][10] = 4;				polys[58][11] = 10;
	polys[59][0] = 3;				polys[59][1] = 1;				polys[59][2] = 9;				polys[59][3] = 9;				polys[59][4] = 10;				polys[59][5] = 7;				polys[59][6] = 7;				polys[59][7] = 3;				polys[59][8] = 9;
	polys[60][0] = 1;				polys[60][1] = 2;				polys[60][2] = 5;				polys[60][3] = 5;				polys[60][4] = 4;				polys[60][5] = 1;				polys[60][6] = 9;				polys[60][7] = 10;				polys[60][8] = 7;				polys[60][9] = 7;				polys[60][10] = 6;				polys[60][11] = 9;
	polys[61][0] = 0;				polys[61][1] = 2;				polys[61][2] = 5;				polys[61][3] = 9;				polys[61][4] = 10;				polys[61][5] = 7;				polys[61][6] = 7;				polys[61][7] = 6;				polys[61][8] = 9;
	polys[62][0] = 1;				polys[62][1] = 0;				polys[62][2] = 4;				polys[62][3] = 9;				polys[62][4] = 10;				polys[62][5] = 7;				polys[62][6] = 7;				polys[62][7] = 6;				polys[62][8] = 9;
	polys[63][0] = 9;				polys[63][1] = 10;				polys[63][2] = 7;				polys[63][3] = 7;				polys[63][4] = 6;				polys[63][5] = 9;
	polys[64][0] = 11;				polys[64][1] = 9;				polys[64][2] = 6;
	polys[65][0] = 0;				polys[65][1] = 1;				polys[65][2] = 6;				polys[65][3] = 6;				polys[65][4] = 11;				polys[65][5] = 9;				polys[65][6] = 9;				polys[65][7] = 4;				polys[65][8] = 0;				polys[65][9] = 0;				polys[65][10] = 6;				polys[65][11] = 9;
	polys[66][0] = 2;				polys[66][1] = 0;				polys[66][2] = 5;				polys[66][3] = 11;				polys[66][4] = 9;				polys[66][5] = 6;
	polys[67][0] = 2;				polys[67][1] = 1;				polys[67][2] = 6;				polys[67][3] = 6;				polys[67][4] = 11;				polys[67][5] = 9;				polys[67][6] = 9;				polys[67][7] = 4;				polys[67][8] = 5;				polys[67][9] = 5;				polys[67][10] = 2;				polys[67][11] = 6;				polys[67][12] = 6;				polys[67][13] = 9;				polys[67][14] = 5;
	polys[68][0] = 1;				polys[68][1] = 3;				polys[68][2] = 11;				polys[68][3] = 11;				polys[68][4] = 9;				polys[68][5] = 1;
	polys[69][0] = 0;				polys[69][1] = 3;				polys[69][2] = 11;				polys[69][3] = 11;				polys[69][4] = 9;				polys[69][5] = 4;				polys[69][6] = 4;				polys[69][7] = 0;				polys[69][8] = 11;
	polys[70][0] = 2;				polys[70][1] = 3;				polys[70][2] = 11;				polys[70][3] = 11;				polys[70][4] = 9;				polys[70][5] = 1;				polys[70][6] = 1;				polys[70][7] = 0;				polys[70][8] = 5;				polys[70][9] = 5;				polys[70][10] = 2;				polys[70][11] = 11;				polys[70][12] = 11;				polys[70][13] = 1;				polys[70][14] = 5;
	polys[71][0] = 2;				polys[71][1] = 3;				polys[71][2] = 11;				polys[71][3] = 11;				polys[71][4] = 9;				polys[71][5] = 4;				polys[71][6] = 4;				polys[71][7] = 5;				polys[71][8] = 2;				polys[71][9] = 2;				polys[71][10] = 11;				polys[71][11] = 4;
	polys[72][0] = 3;				polys[72][1] = 2;				polys[72][2] = 7;				polys[72][3] = 7;				polys[72][4] = 11;				polys[72][5] = 9;				polys[72][6] = 9;				polys[72][7] = 6;				polys[72][8] = 3;				polys[72][9] = 3;				polys[72][10] = 7;				polys[72][11] = 9;
	polys[73][0] = 0;				polys[73][1] = 2;				polys[73][2] = 7;				polys[73][3] = 7;				polys[73][4] = 11;				polys[73][5] = 9;				polys[73][6] = 9;				polys[73][7] = 4;				polys[73][8] = 0;				polys[73][9] = 0;				polys[73][10] = 7;				polys[73][11] = 9;				polys[73][12] = 3;				polys[73][13] = 1;				polys[73][14] = 6;
	polys[74][0] = 3;				polys[74][1] = 0;				polys[74][2] = 5;				polys[74][3] = 5;				polys[74][4] = 7;				polys[74][5] = 11;				polys[74][6] = 11;				polys[74][7] = 9;				polys[74][8] = 6;				polys[74][9] = 6;				polys[74][10] = 3;				polys[74][11] = 5;				polys[74][12] = 5;				polys[74][13] = 11;				polys[74][14] = 6;
	polys[75][0] = 3;				polys[75][1] = 1;				polys[75][2] = 6;				polys[75][3] = 11;				polys[75][4] = 9;				polys[75][5] = 4;				polys[75][6] = 4;				polys[75][7] = 5;				polys[75][8] = 7;				polys[75][9] = 7;				polys[75][10] = 11;				polys[75][11] = 4;
	polys[76][0] = 1;				polys[76][1] = 2;				polys[76][2] = 7;				polys[76][3] = 7;				polys[76][4] = 11;				polys[76][5] = 9;				polys[76][6] = 9;				polys[76][7] = 1;				polys[76][8] = 7;
	polys[77][0] = 0;				polys[77][1] = 2;				polys[77][2] = 7;				polys[77][3] = 7;				polys[77][4] = 11;				polys[77][5] = 9;				polys[77][6] = 9;				polys[77][7] = 4;				polys[77][8] = 0;				polys[77][9] = 0;				polys[77][10] = 7;				polys[77][11] = 9;
	polys[78][0] = 1;				polys[78][1] = 0;				polys[78][2] = 5;				polys[78][3] = 5;				polys[78][4] = 7;				polys[78][5] = 11;				polys[78][6] = 11;				polys[78][7] = 9;				polys[78][8] = 1;				polys[78][9] = 1;				polys[78][10] = 5;				polys[78][11] = 11;
	polys[79][0] = 11;				polys[79][1] = 9;				polys[79][2] = 4;				polys[79][3] = 4;				polys[79][4] = 5;				polys[79][5] = 7;				polys[79][6] = 7;				polys[79][7] = 11;				polys[79][8] = 4;
	polys[80][0] = 11;				polys[80][1] = 8;				polys[80][2] = 4;				polys[80][3] = 4;				polys[80][4] = 6;				polys[80][5] = 11;
	polys[81][0] = 0;				polys[81][1] = 1;				polys[81][2] = 6;				polys[81][3] = 6;				polys[81][4] = 11;				polys[81][5] = 8;				polys[81][6] = 8;				polys[81][7] = 0;				polys[81][8] = 6;
	polys[82][0] = 2;				polys[82][1] = 0;				polys[82][2] = 4;				polys[82][3] = 4;				polys[82][4] = 6;				polys[82][5] = 11;				polys[82][6] = 11;				polys[82][7] = 8;				polys[82][8] = 5;				polys[82][9] = 5;				polys[82][10] = 2;				polys[82][11] = 4;				polys[82][12] = 4;				polys[82][13] = 11;				polys[82][14] = 5;
	polys[83][0] = 2;				polys[83][1] = 1;				polys[83][2] = 6;				polys[83][3] = 6;				polys[83][4] = 11;				polys[83][5] = 8;				polys[83][6] = 8;				polys[83][7] = 5;				polys[83][8] = 2;				polys[83][9] = 2;				polys[83][10] = 6;				polys[83][11] = 8;
	polys[84][0] = 1;				polys[84][1] = 3;				polys[84][2] = 11;				polys[84][3] = 11;				polys[84][4] = 8;				polys[84][5] = 4;				polys[84][6] = 4;				polys[84][7] = 1;				polys[84][8] = 11;
	polys[85][0] = 0;				polys[85][1] = 3;				polys[85][2] = 11;				polys[85][3] = 11;				polys[85][4] = 8;				polys[85][5] = 0;
	polys[86][0] = 2;				polys[86][1] = 3;				polys[86][2] = 11;				polys[86][3] = 11;				polys[86][4] = 8;				polys[86][5] = 5;				polys[86][6] = 5;				polys[86][7] = 2;				polys[86][8] = 11;				polys[86][9] = 1;				polys[86][10] = 0;				polys[86][11] = 4;
	polys[87][0] = 2;				polys[87][1] = 3;				polys[87][2] = 11;				polys[87][3] = 11;				polys[87][4] = 8;				polys[87][5] = 5;				polys[87][6] = 5;				polys[87][7] = 2;				polys[87][8] = 11;
	polys[88][0] = 3;				polys[88][1] = 2;				polys[88][2] = 7;				polys[88][3] = 7;				polys[88][4] = 11;				polys[88][5] = 8;				polys[88][6] = 8;				polys[88][7] = 4;				polys[88][8] = 6;				polys[88][9] = 6;				polys[88][10] = 3;				polys[88][11] = 7;				polys[88][12] = 7;				polys[88][13] = 8;				polys[88][14] = 6;
	polys[89][0] = 0;				polys[89][1] = 2;				polys[89][2] = 7;				polys[89][3] = 7;				polys[89][4] = 11;				polys[89][5] = 8;				polys[89][6] = 8;				polys[89][7] = 0;				polys[89][8] = 7;				polys[89][9] = 3;				polys[89][10] = 1;				polys[89][11] = 6;
	polys[90][0] = 3;				polys[90][1] = 0;				polys[90][2] = 4;				polys[90][3] = 4;				polys[90][4] = 6;				polys[90][5] = 3;				polys[90][6] = 11;				polys[90][7] = 8;				polys[90][8] = 5;				polys[90][9] = 5;				polys[90][10] = 7;				polys[90][11] = 11;
	polys[91][0] = 3;				polys[91][1] = 1;				polys[91][2] = 6;				polys[91][3] = 11;				polys[91][4] = 8;				polys[91][5] = 5;				polys[91][6] = 5;				polys[91][7] = 7;				polys[91][8] = 11;
	polys[92][0] = 1;				polys[92][1] = 2;				polys[92][2] = 7;				polys[92][3] = 7;				polys[92][4] = 11;				polys[92][5] = 8;				polys[92][6] = 8;				polys[92][7] = 4;				polys[92][8] = 1;				polys[92][9] = 1;				polys[92][10] = 7;				polys[92][11] = 8;
	polys[93][0] = 0;				polys[93][1] = 2;				polys[93][2] = 7;				polys[93][3] = 7;				polys[93][4] = 11;				polys[93][5] = 8;				polys[93][6] = 8;				polys[93][7] = 0;				polys[93][8] = 7;
	polys[94][0] = 1;				polys[94][1] = 0;				polys[94][2] = 4;				polys[94][3] = 11;				polys[94][4] = 8;				polys[94][5] = 5;				polys[94][6] = 5;				polys[94][7] = 7;				polys[94][8] = 11;
	polys[95][0] = 11;				polys[95][1] = 8;				polys[95][2] = 5;				polys[95][3] = 5;				polys[95][4] = 7;				polys[95][5] = 11;
	polys[96][0] = 11;				polys[96][1] = 10;				polys[96][2] = 5;				polys[96][3] = 5;				polys[96][4] = 8;				polys[96][5] = 9;				polys[96][6] = 9;				polys[96][7] = 6;				polys[96][8] = 11;				polys[96][9] = 11;				polys[96][10] = 5;				polys[96][11] = 9;
	polys[97][0] = 0;				polys[97][1] = 1;				polys[97][2] = 6;				polys[97][3] = 6;				polys[97][4] = 11;				polys[97][5] = 10;				polys[97][6] = 10;				polys[97][7] = 5;				polys[97][8] = 0;				polys[97][9] = 0;				polys[97][10] = 6;				polys[97][11] = 10;				polys[97][12] = 8;				polys[97][13] = 9;				polys[97][14] = 4;
	polys[98][0] = 2;				polys[98][1] = 0;				polys[98][2] = 8;				polys[98][3] = 8;				polys[98][4] = 9;				polys[98][5] = 6;				polys[98][6] = 6;				polys[98][7] = 11;				polys[98][8] = 10;				polys[98][9] = 10;				polys[98][10] = 2;				polys[98][11] = 8;				polys[98][12] = 8;				polys[98][13] = 6;				polys[98][14] = 10;
	polys[99][0] = 2;				polys[99][1] = 1;				polys[99][2] = 6;				polys[99][3] = 6;				polys[99][4] = 11;				polys[99][5] = 10;				polys[99][6] = 10;				polys[99][7] = 2;				polys[99][8] = 6;				polys[99][9] = 8;				polys[99][10] = 9;				polys[99][11] = 4;
	polys[100][0] = 1;				polys[100][1] = 3;				polys[100][2] = 11;				polys[100][3] = 11;				polys[100][4] = 10;				polys[100][5] = 5;				polys[100][6] = 5;				polys[100][7] = 8;				polys[100][8] = 9;				polys[100][9] = 9;				polys[100][10] = 1;				polys[100][11] = 11;			polys[100][12] = 11;			polys[100][13] = 5;				polys[100][14] = 9;
	polys[101][0] = 0;				polys[101][1] = 3;				polys[101][2] = 11;				polys[101][3] = 11;				polys[101][4] = 10;				polys[101][5] = 5;				polys[101][6] = 5;				polys[101][7] = 0;				polys[101][8] = 11;				polys[101][9] = 8;				polys[101][10] = 9;				polys[101][11] = 4;
	polys[102][0] = 2;				polys[102][1] = 3;				polys[102][2] = 11;				polys[102][3] = 11;				polys[102][4] = 10;				polys[102][5] = 2;				polys[102][6] = 1;				polys[102][7] = 0;				polys[102][8] = 8;				polys[102][9] = 8;				polys[102][10] = 9;				polys[102][11] = 1;
	polys[103][0] = 2;				polys[103][1] = 3;				polys[103][2] = 11;				polys[103][3] = 11;				polys[103][4] = 10;				polys[103][5] = 2;				polys[103][6] = 8;				polys[103][7] = 9;				polys[103][8] = 4;
	polys[104][0] = 3;				polys[104][1] = 2;				polys[104][2] = 5;				polys[104][3] = 5;				polys[104][4] = 8;				polys[104][5] = 9;				polys[104][6] = 9;				polys[104][7] = 6;				polys[104][8] = 3;				polys[104][9] = 3;				polys[104][10] = 5;				polys[104][11] = 9;				polys[104][12] = 11;			polys[104][13] = 10;			polys[104][14] = 7;
	polys[105][0] = 0;				polys[105][1] = 2;				polys[105][2] = 5;				polys[105][3] = 3;				polys[105][4] = 1;				polys[105][5] = 6;				polys[105][6] = 11;				polys[105][7] = 10;				polys[105][8] = 7;				polys[105][9] = 8;				polys[105][10] = 9;				polys[105][11] = 4;
	polys[106][0] = 3;				polys[106][1] = 0;				polys[106][2] = 8;				polys[106][3] = 8;				polys[106][4] = 9;				polys[106][5] = 6;				polys[106][6] = 6;				polys[106][7] = 3;				polys[106][8] = 8;				polys[106][9] = 11;				polys[106][10] = 10;			polys[106][11] = 7;
	polys[107][0] = 3;				polys[107][1] = 1;				polys[107][2] = 6;				polys[107][3] = 11;				polys[107][4] = 10;				polys[107][5] = 7;				polys[107][6] = 8;				polys[107][7] = 9;				polys[107][8] = 4;
	polys[108][0] = 1;				polys[108][1] = 2;				polys[108][2] = 5;				polys[108][3] = 5;				polys[108][4] = 8;				polys[108][5] = 9;				polys[108][6] = 9;				polys[108][7] = 1;				polys[108][8] = 5;				polys[108][9] = 11;				polys[108][10] = 10;			polys[108][11] = 7;
	polys[109][0] = 0;				polys[109][1] = 2;				polys[109][2] = 5;				polys[109][3] = 11;				polys[109][4] = 10;				polys[109][5] = 7;				polys[109][6] = 8;				polys[109][7] = 9;				polys[109][8] = 4;
	polys[110][0] = 1;				polys[110][1] = 0;				polys[110][2] = 8;				polys[110][3] = 8;				polys[110][4] = 9;				polys[110][5] = 1;				polys[110][6] = 11;				polys[110][7] = 10;				polys[110][8] = 7;
	polys[111][0] = 11;				polys[111][1] = 10;				polys[111][2] = 7;				polys[111][3] = 8;				polys[111][4] = 9;				polys[111][5] = 4;
	polys[112][0] = 11;				polys[112][1] = 10;				polys[112][2] = 5;				polys[112][3] = 5;				polys[112][4] = 4;				polys[112][5] = 6;				polys[112][6] = 6;				polys[112][7] = 11;				polys[112][8] = 5;
	polys[113][0] = 0;				polys[113][1] = 1;				polys[113][2] = 6;				polys[113][3] = 6;				polys[113][4] = 11;				polys[113][5] = 10;				polys[113][6] = 10;				polys[113][7] = 5;				polys[113][8] = 0;				polys[113][9] = 0;				polys[113][10] = 6;				polys[113][11] = 10;
	polys[114][0] = 2;				polys[114][1] = 0;				polys[114][2] = 4;				polys[114][3] = 4;				polys[114][4] = 6;				polys[114][5] = 11;				polys[114][6] = 11;				polys[114][7] = 10;				polys[114][8] = 2;				polys[114][9] = 2;				polys[114][10] = 4;				polys[114][11] = 11;
	polys[115][0] = 2;				polys[115][1] = 1;				polys[115][2] = 6;				polys[115][3] = 6;				polys[115][4] = 11;				polys[115][5] = 10;				polys[115][6] = 10;				polys[115][7] = 2;				polys[115][8] = 6;
	polys[116][0] = 1;				polys[116][1] = 3;				polys[116][2] = 11;				polys[116][3] = 11;				polys[116][4] = 10;				polys[116][5] = 5;				polys[116][6] = 5;				polys[116][7] = 4;				polys[116][8] = 1;				polys[116][9] = 1;				polys[116][10] = 11;			polys[116][11] = 5;
	polys[117][0] = 0;				polys[117][1] = 3;				polys[117][2] = 11;				polys[117][3] = 11;				polys[117][4] = 10;				polys[117][5] = 5;				polys[117][6] = 5;				polys[117][7] = 0;				polys[117][8] = 11;
	polys[118][0] = 2;				polys[118][1] = 3;				polys[118][2] = 11;				polys[118][3] = 11;				polys[118][4] = 10;				polys[118][5] = 2;				polys[118][6] = 1;				polys[118][7] = 0;				polys[118][8] = 4;
	polys[119][0] = 2;				polys[119][1] = 3;				polys[119][2] = 11;				polys[119][3] = 11;				polys[119][4] = 10;				polys[119][5] = 2;
	polys[120][0] = 3;				polys[120][1] = 2;				polys[120][2] = 5;				polys[120][3] = 5;				polys[120][4] = 4;				polys[120][5] = 6;				polys[120][6] = 6;				polys[120][7] = 3;				polys[120][8] = 5;				polys[120][9] = 11;				polys[120][10] = 10;			polys[120][11] = 7;
	polys[121][0] = 0;				polys[121][1] = 2;				polys[121][2] = 5;				polys[121][3] = 3;				polys[121][4] = 1;				polys[121][5] = 6;				polys[121][6] = 11;				polys[121][7] = 10;				polys[121][8] = 7;
	polys[122][0] = 3;				polys[122][1] = 0;				polys[122][2] = 4;				polys[122][3] = 4;				polys[122][4] = 6;				polys[122][5] = 3;				polys[122][6] = 11;				polys[122][7] = 10;				polys[122][8] = 7;
	polys[123][0] = 3;				polys[123][1] = 1;				polys[123][2] = 6;				polys[123][3] = 11;				polys[123][4] = 10;				polys[123][5] = 7;
	polys[124][0] = 1;				polys[124][1] = 2;				polys[124][2] = 5;				polys[124][3] = 5;				polys[124][4] = 4;				polys[124][5] = 1;				polys[124][6] = 11;				polys[124][7] = 10;				polys[124][8] = 7;
	polys[125][0] = 0;				polys[125][1] = 2;				polys[125][2] = 5;				polys[125][3] = 11;				polys[125][4] = 10;				polys[125][5] = 7;
	polys[126][0] = 1;				polys[126][1] = 0;				polys[126][2] = 4;				polys[126][3] = 11;				polys[126][4] = 10;				polys[126][5] = 7;
	polys[127][0] = 11;				polys[127][1] = 10;				polys[127][2] = 7;
	polys[128][0] = 10;				polys[128][1] = 11;				polys[128][2] = 7;
	polys[129][0] = 0;				polys[129][1] = 1;				polys[129][2] = 4;				polys[129][3] = 10;				polys[129][4] = 11;				polys[129][5] = 7;
	polys[130][0] = 2;				polys[130][1] = 0;				polys[130][2] = 5;				polys[130][3] = 5;				polys[130][4] = 10;				polys[130][5] = 11;				polys[130][6] = 11;				polys[130][7] = 7;				polys[130][8] = 2;				polys[130][9] = 2;				polys[130][10] = 5;				polys[130][11] = 11;
	polys[131][0] = 2;				polys[131][1] = 1;				polys[131][2] = 4;				polys[131][3] = 4;				polys[131][4] = 5;				polys[131][5] = 10;				polys[131][6] = 10;				polys[131][7] = 11;				polys[131][8] = 7;				polys[131][9] = 7;				polys[131][10] = 2;				polys[131][11] = 4;				polys[131][12] = 4;				polys[131][13] = 10;			polys[131][14] = 7;
	polys[132][0] = 1;				polys[132][1] = 3;				polys[132][2] = 7;				polys[132][3] = 7;				polys[132][4] = 10;				polys[132][5] = 11;				polys[132][6] = 11;				polys[132][7] = 6;				polys[132][8] = 1;				polys[132][9] = 1;				polys[132][10] = 7;				polys[132][11] = 11;
	polys[133][0] = 0;				polys[133][1] = 3;				polys[133][2] = 7;				polys[133][3] = 7;				polys[133][4] = 10;				polys[133][5] = 11;				polys[133][6] = 11;				polys[133][7] = 6;				polys[133][8] = 4;				polys[133][9] = 4;				polys[133][10] = 0;				polys[133][11] = 7;				polys[133][12] = 7;				polys[133][13] = 11;			polys[133][14] = 4;
	polys[134][0] = 2;				polys[134][1] = 3;				polys[134][2] = 7;				polys[134][3] = 1;				polys[134][4] = 0;				polys[134][5] = 5;				polys[134][6] = 5;				polys[134][7] = 10;				polys[134][8] = 11;				polys[134][9] = 11;				polys[134][10] = 6;				polys[134][11] = 1;				polys[134][12] = 1;				polys[134][13] = 5;				polys[134][14] = 11;
	polys[135][0] = 2;				polys[135][1] = 3;				polys[135][2] = 7;				polys[135][3] = 10;				polys[135][4] = 11;				polys[135][5] = 6;				polys[135][6] = 6;				polys[135][7] = 4;				polys[135][8] = 5;				polys[135][9] = 5;				polys[135][10] = 10;			polys[135][11] = 6;
	polys[136][0] = 3;				polys[136][1] = 2;				polys[136][2] = 10;				polys[136][3] = 10;				polys[136][4] = 11;				polys[136][5] = 3;
	polys[137][0] = 0;				polys[137][1] = 2;				polys[137][2] = 10;				polys[137][3] = 10;				polys[137][4] = 11;				polys[137][5] = 3;				polys[137][6] = 3;				polys[137][7] = 1;				polys[137][8] = 4;				polys[137][9] = 4;				polys[137][10] = 0;				polys[137][11] = 10;			polys[137][12] = 10;			polys[137][13] = 3;				polys[137][14] = 4;
	polys[138][0] = 3;				polys[138][1] = 0;				polys[138][2] = 5;				polys[138][3] = 5;				polys[138][4] = 10;				polys[138][5] = 11;				polys[138][6] = 11;				polys[138][7] = 3;				polys[138][8] = 5;
	polys[139][0] = 3;				polys[139][1] = 1;				polys[139][2] = 4;				polys[139][3] = 4;				polys[139][4] = 5;				polys[139][5] = 10;				polys[139][6] = 10;				polys[139][7] = 11;				polys[139][8] = 3;				polys[139][9] = 3;				polys[139][10] = 4;				polys[139][11] = 10;
	polys[140][0] = 1;				polys[140][1] = 2;				polys[140][2] = 10;				polys[140][3] = 10;				polys[140][4] = 11;				polys[140][5] = 6;				polys[140][6] = 6;				polys[140][7] = 1;				polys[140][8] = 10;
	polys[141][0] = 0;				polys[141][1] = 2;				polys[141][2] = 10;				polys[141][3] = 10;				polys[141][4] = 11;				polys[141][5] = 6;				polys[141][6] = 6;				polys[141][7] = 4;				polys[141][8] = 0;				polys[141][9] = 0;				polys[141][10] = 10;			polys[141][11] = 6;
	polys[142][0] = 1;				polys[142][1] = 0;				polys[142][2] = 5;				polys[142][3] = 5;				polys[142][4] = 10;				polys[142][5] = 11;				polys[142][6] = 11;				polys[142][7] = 6;				polys[142][8] = 1;				polys[142][9] = 1;				polys[142][10] = 5;				polys[142][11] = 11;
	polys[143][0] = 10;				polys[143][1] = 11;				polys[143][2] = 6;				polys[143][3] = 6;				polys[143][4] = 4;				polys[143][5] = 5;				polys[143][6] = 5;				polys[143][7] = 10;				polys[143][8] = 6;
	polys[144][0] = 9;				polys[144][1] = 11;				polys[144][2] = 7;				polys[144][3] = 7;				polys[144][4] = 10;				polys[144][5] = 8;				polys[144][6] = 8;				polys[144][7] = 4;				polys[144][8] = 9;				polys[144][9] = 9;				polys[144][10] = 7;				polys[144][11] = 8;
	polys[145][0] = 0;				polys[145][1] = 1;				polys[145][2] = 9;				polys[145][3] = 9;				polys[145][4] = 11;				polys[145][5] = 7;				polys[145][6] = 7;				polys[145][7] = 10;				polys[145][8] = 8;				polys[145][9] = 8;				polys[145][10] = 0;				polys[145][11] = 9;				polys[145][12] = 9;				polys[145][13] = 7;				polys[145][14] = 8;
	polys[146][0] = 2;				polys[146][1] = 0;				polys[146][2] = 4;				polys[146][3] = 4;				polys[146][4] = 9;				polys[146][5] = 11;				polys[146][6] = 11;				polys[146][7] = 7;				polys[146][8] = 2;				polys[146][9] = 2;				polys[146][10] = 4;				polys[146][11] = 11;			polys[146][12] = 10;			polys[146][13] = 8;				polys[146][14] = 5;
	polys[147][0] = 2;				polys[147][1] = 1;				polys[147][2] = 9;				polys[147][3] = 9;				polys[147][4] = 11;				polys[147][5] = 7;				polys[147][6] = 7;				polys[147][7] = 2;				polys[147][8] = 9;				polys[147][9] = 10;				polys[147][10] = 8;				polys[147][11] = 5;
	polys[148][0] = 1;				polys[148][1] = 3;				polys[148][2] = 7;				polys[148][3] = 7;				polys[148][4] = 10;				polys[148][5] = 8;				polys[148][6] = 8;				polys[148][7] = 4;				polys[148][8] = 1;				polys[148][9] = 1;				polys[148][10] = 7;				polys[148][11] = 8;				polys[148][12] = 9;				polys[148][13] = 11;			polys[148][14] = 6;
	polys[149][0] = 0;				polys[149][1] = 3;				polys[149][2] = 7;				polys[149][3] = 7;				polys[149][4] = 10;				polys[149][5] = 8;				polys[149][6] = 8;				polys[149][7] = 0;				polys[149][8] = 7;				polys[149][9] = 9;				polys[149][10] = 11;			polys[149][11] = 6;
	polys[150][0] = 2;				polys[150][1] = 3;				polys[150][2] = 7;				polys[150][3] = 1;				polys[150][4] = 0;				polys[150][5] = 4;				polys[150][6] = 9;				polys[150][7] = 11;				polys[150][8] = 6;				polys[150][9] = 10;				polys[150][10] = 8;				polys[150][11] = 5;
	polys[151][0] = 2;				polys[151][1] = 3;				polys[151][2] = 7;				polys[151][3] = 9;				polys[151][4] = 11;				polys[151][5] = 6;				polys[151][6] = 10;				polys[151][7] = 8;				polys[151][8] = 5;
	polys[152][0] = 3;				polys[152][1] = 2;				polys[152][2] = 10;				polys[152][3] = 10;				polys[152][4] = 8;				polys[152][5] = 4;				polys[152][6] = 4;				polys[152][7] = 9;				polys[152][8] = 11;				polys[152][9] = 11;				polys[152][10] = 3;				polys[152][11] = 10;			polys[152][12] = 10;			polys[152][13] = 4;				polys[152][14] = 11;
	polys[153][0] = 0;				polys[153][1] = 2;				polys[153][2] = 10;				polys[153][3] = 10;				polys[153][4] = 8;				polys[153][5] = 0;				polys[153][6] = 3;				polys[153][7] = 1;				polys[153][8] = 9;				polys[153][9] = 9;				polys[153][10] = 11;			polys[153][11] = 3;
	polys[154][0] = 3;				polys[154][1] = 0;				polys[154][2] = 4;				polys[154][3] = 4;				polys[154][4] = 9;				polys[154][5] = 11;				polys[154][6] = 11;				polys[154][7] = 3;				polys[154][8] = 4;				polys[154][9] = 10;				polys[154][10] = 8;				polys[154][11] = 5;
	polys[155][0] = 3;				polys[155][1] = 1;				polys[155][2] = 9;				polys[155][3] = 9;				polys[155][4] = 11;				polys[155][5] = 3;				polys[155][6] = 10;				polys[155][7] = 8;				polys[155][8] = 5;
	polys[156][0] = 1;				polys[156][1] = 2;				polys[156][2] = 10;				polys[156][3] = 10;				polys[156][4] = 8;				polys[156][5] = 4;				polys[156][6] = 4;				polys[156][7] = 1;				polys[156][8] = 10;				polys[156][9] = 9;				polys[156][10] = 11;			polys[156][11] = 6;
	polys[157][0] = 0;				polys[157][1] = 2;				polys[157][2] = 10;				polys[157][3] = 10;				polys[157][4] = 8;				polys[157][5] = 0;				polys[157][6] = 9;				polys[157][7] = 11;				polys[157][8] = 6;
	polys[158][0] = 1;				polys[158][1] = 0;				polys[158][2] = 4;				polys[158][3] = 9;				polys[158][4] = 11;				polys[158][5] = 6;				polys[158][6] = 10;				polys[158][7] = 8;				polys[158][8] = 5;
	polys[159][0] = 9;				polys[159][1] = 11;				polys[159][2] = 6;				polys[159][3] = 10;				polys[159][4] = 8;				polys[159][5] = 5;
	polys[160][0] = 8;				polys[160][1] = 11;				polys[160][2] = 7;				polys[160][3] = 7;				polys[160][4] = 5;				polys[160][5] = 8;
	polys[161][0] = 0;				polys[161][1] = 1;				polys[161][2] = 4;				polys[161][3] = 4;				polys[161][4] = 8;				polys[161][5] = 11;				polys[161][6] = 11;				polys[161][7] = 7;				polys[161][8] = 5;				polys[161][9] = 5;				polys[161][10] = 0;				polys[161][11] = 4;				polys[161][12] = 4;				polys[161][13] = 11;			polys[161][14] = 5;
	polys[162][0] = 2;				polys[162][1] = 0;				polys[162][2] = 8;				polys[162][3] = 8;				polys[162][4] = 11;				polys[162][5] = 7;				polys[162][6] = 7;				polys[162][7] = 2;				polys[162][8] = 8;
	polys[163][0] = 2;				polys[163][1] = 1;				polys[163][2] = 4;				polys[163][3] = 4;				polys[163][4] = 8;				polys[163][5] = 11;				polys[163][6] = 11;				polys[163][7] = 7;				polys[163][8] = 2;				polys[163][9] = 2;				polys[163][10] = 4;				polys[163][11] = 11;
	polys[164][0] = 1;				polys[164][1] = 3;				polys[164][2] = 7;				polys[164][3] = 7;				polys[164][4] = 5;				polys[164][5] = 8;				polys[164][6] = 8;				polys[164][7] = 11;				polys[164][8] = 6;				polys[164][9] = 6;				polys[164][10] = 1;				polys[164][11] = 7;				polys[164][12] = 7;				polys[164][13] = 8;				polys[164][14] = 6;
	polys[165][0] = 0;				polys[165][1] = 3;				polys[165][2] = 7;				polys[165][3] = 7;				polys[165][4] = 5;				polys[165][5] = 0;				polys[165][6] = 8;				polys[165][7] = 11;				polys[165][8] = 6;				polys[165][9] = 6;				polys[165][10] = 4;				polys[165][11] = 8;
	polys[166][0] = 2;				polys[166][1] = 3;				polys[166][2] = 7;				polys[166][3] = 1;				polys[166][4] = 0;				polys[166][5] = 8;				polys[166][6] = 8;				polys[166][7] = 11;				polys[166][8] = 6;				polys[166][9] = 6;				polys[166][10] = 1;				polys[166][11] = 8;
	polys[167][0] = 2;				polys[167][1] = 3;				polys[167][2] = 7;				polys[167][3] = 8;				polys[167][4] = 11;				polys[167][5] = 6;				polys[167][6] = 6;				polys[167][7] = 4;				polys[167][8] = 8;
	polys[168][0] = 3;				polys[168][1] = 2;				polys[168][2] = 5;				polys[168][3] = 5;				polys[168][4] = 8;				polys[168][5] = 11;				polys[168][6] = 11;				polys[168][7] = 3;				polys[168][8] = 5;
	polys[169][0] = 0;				polys[169][1] = 2;				polys[169][2] = 5;				polys[169][3] = 3;				polys[169][4] = 1;				polys[169][5] = 4;				polys[169][6] = 4;				polys[169][7] = 8;				polys[169][8] = 11;				polys[169][9] = 11;				polys[169][10] = 3;				polys[169][11] = 4;
	polys[170][0] = 3;				polys[170][1] = 0;				polys[170][2] = 8;				polys[170][3] = 8;				polys[170][4] = 11;				polys[170][5] = 3;
	polys[171][0] = 3;				polys[171][1] = 1;				polys[171][2] = 4;				polys[171][3] = 4;				polys[171][4] = 8;				polys[171][5] = 11;				polys[171][6] = 11;				polys[171][7] = 3;				polys[171][8] = 4;
	polys[172][0] = 1;				polys[172][1] = 2;				polys[172][2] = 5;				polys[172][3] = 5;				polys[172][4] = 8;				polys[172][5] = 11;				polys[172][6] = 11;				polys[172][7] = 6;				polys[172][8] = 1;				polys[172][9] = 1;				polys[172][10] = 5;				polys[172][11] = 11;
	polys[173][0] = 0;				polys[173][1] = 2;				polys[173][2] = 5;				polys[173][3] = 8;				polys[173][4] = 11;				polys[173][5] = 6;				polys[173][6] = 6;				polys[173][7] = 4;				polys[173][8] = 8;
	polys[174][0] = 1;				polys[174][1] = 0;				polys[174][2] = 8;				polys[174][3] = 8;				polys[174][4] = 11;				polys[174][5] = 6;				polys[174][6] = 6;				polys[174][7] = 1;				polys[174][8] = 8;
	polys[175][0] = 8;				polys[175][1] = 11;				polys[175][2] = 6;				polys[175][3] = 6;				polys[175][4] = 4;				polys[175][5] = 8;
	polys[176][0] = 9;				polys[176][1] = 11;				polys[176][2] = 7;				polys[176][3] = 7;				polys[176][4] = 5;				polys[176][5] = 4;				polys[176][6] = 4;				polys[176][7] = 9;				polys[176][8] = 7;
	polys[177][0] = 0;				polys[177][1] = 1;				polys[177][2] = 9;				polys[177][3] = 9;				polys[177][4] = 11;				polys[177][5] = 7;				polys[177][6] = 7;				polys[177][7] = 5;				polys[177][8] = 0;				polys[177][9] = 0;				polys[177][10] = 9;				polys[177][11] = 7;
	polys[178][0] = 2;				polys[178][1] = 0;				polys[178][2] = 4;				polys[178][3] = 4;				polys[178][4] = 9;				polys[178][5] = 11;				polys[178][6] = 11;				polys[178][7] = 7;				polys[178][8] = 2;				polys[178][9] = 2;				polys[178][10] = 4;				polys[178][11] = 11;
	polys[179][0] = 2;				polys[179][1] = 1;				polys[179][2] = 9;				polys[179][3] = 9;				polys[179][4] = 11;				polys[179][5] = 7;				polys[179][6] = 7;				polys[179][7] = 2;				polys[179][8] = 9;
	polys[180][0] = 1;				polys[180][1] = 3;				polys[180][2] = 7;				polys[180][3] = 7;				polys[180][4] = 5;				polys[180][5] = 4;				polys[180][6] = 4;				polys[180][7] = 1;				polys[180][8] = 7;				polys[180][9] = 9;				polys[180][10] = 11;			polys[180][11] = 6;
	polys[181][0] = 0;				polys[181][1] = 3;				polys[181][2] = 7;				polys[181][3] = 7;				polys[181][4] = 5;				polys[181][5] = 0;				polys[181][6] = 9;				polys[181][7] = 11;				polys[181][8] = 6;
	polys[182][0] = 2;				polys[182][1] = 3;				polys[182][2] = 7;				polys[182][3] = 1;				polys[182][4] = 0;				polys[182][5] = 4;				polys[182][6] = 9;				polys[182][7] = 11;				polys[182][8] = 6;
	polys[183][0] = 2;				polys[183][1] = 3;				polys[183][2] = 7;				polys[183][3] = 9;				polys[183][4] = 11;				polys[183][5] = 6;
	polys[184][0] = 3;				polys[184][1] = 2;				polys[184][2] = 5;				polys[184][3] = 5;				polys[184][4] = 4;				polys[184][5] = 9;				polys[184][6] = 9;				polys[184][7] = 11;				polys[184][8] = 3;				polys[184][9] = 3;				polys[184][10] = 5;				polys[184][11] = 9;
	polys[185][0] = 0;				polys[185][1] = 2;				polys[185][2] = 5;				polys[185][3] = 3;				polys[185][4] = 1;				polys[185][5] = 9;				polys[185][6] = 9;				polys[185][7] = 11;				polys[185][8] = 3;
	polys[186][0] = 3;				polys[186][1] = 0;				polys[186][2] = 4;				polys[186][3] = 4;				polys[186][4] = 9;				polys[186][5] = 11;				polys[186][6] = 11;				polys[186][7] = 3;				polys[186][8] = 4;
	polys[187][0] = 3;				polys[187][1] = 1;				polys[187][2] = 9;				polys[187][3] = 9;				polys[187][4] = 11;				polys[187][5] = 3;
	polys[188][0] = 1;				polys[188][1] = 2;				polys[188][2] = 5;				polys[188][3] = 5;				polys[188][4] = 4;				polys[188][5] = 1;				polys[188][6] = 9;				polys[188][7] = 11;				polys[188][8] = 6;
	polys[189][0] = 0;				polys[189][1] = 2;				polys[189][2] = 5;				polys[189][3] = 9;				polys[189][4] = 11;				polys[189][5] = 6;
	polys[190][0] = 1;				polys[190][1] = 0;				polys[190][2] = 4;				polys[190][3] = 9;				polys[190][4] = 11;				polys[190][5] = 6;
	polys[191][0] = 9;				polys[191][1] = 11;				polys[191][2] = 6;
	polys[192][0] = 10;				polys[192][1] = 9;				polys[192][2] = 6;				polys[192][3] = 6;				polys[192][4] = 7;				polys[192][5] = 10;
	polys[193][0] = 0;				polys[193][1] = 1;				polys[193][2] = 6;				polys[193][3] = 6;				polys[193][4] = 7;				polys[193][5] = 10;				polys[193][6] = 10;				polys[193][7] = 9;				polys[193][8] = 4;				polys[193][9] = 4;				polys[193][10] = 0;				polys[193][11] = 6;				polys[193][12] = 6;				polys[193][13] = 10;			polys[193][14] = 4;
	polys[194][0] = 2;				polys[194][1] = 0;				polys[194][2] = 5;				polys[194][3] = 5;				polys[194][4] = 10;				polys[194][5] = 9;				polys[194][6] = 9;				polys[194][7] = 6;				polys[194][8] = 7;				polys[194][9] = 7;				polys[194][10] = 2;				polys[194][11] = 5;				polys[194][12] = 5;				polys[194][13] = 9;				polys[194][14] = 7;
	polys[195][0] = 2;				polys[195][1] = 1;				polys[195][2] = 6;				polys[195][3] = 6;				polys[195][4] = 7;				polys[195][5] = 2;				polys[195][6] = 10;				polys[195][7] = 9;				polys[195][8] = 4;				polys[195][9] = 4;				polys[195][10] = 5;				polys[195][11] = 10;
	polys[196][0] = 1;				polys[196][1] = 3;				polys[196][2] = 7;				polys[196][3] = 7;				polys[196][4] = 10;				polys[196][5] = 9;				polys[196][6] = 9;				polys[196][7] = 1;				polys[196][8] = 7;
	polys[197][0] = 0;				polys[197][1] = 3;				polys[197][2] = 7;				polys[197][3] = 7;				polys[197][4] = 10;				polys[197][5] = 9;				polys[197][6] = 9;				polys[197][7] = 4;				polys[197][8] = 0;				polys[197][9] = 0;				polys[197][10] = 7;				polys[197][11] = 9;
	polys[198][0] = 2;				polys[198][1] = 3;				polys[198][2] = 7;				polys[198][3] = 1;				polys[198][4] = 0;				polys[198][5] = 5;				polys[198][6] = 5;				polys[198][7] = 10;				polys[198][8] = 9;				polys[198][9] = 9;				polys[198][10] = 1;				polys[198][11] = 5;
	polys[199][0] = 2;				polys[199][1] = 3;				polys[199][2] = 7;				polys[199][3] = 10;				polys[199][4] = 9;				polys[199][5] = 4;				polys[199][6] = 4;				polys[199][7] = 5;				polys[199][8] = 10;
	polys[200][0] = 3;				polys[200][1] = 2;				polys[200][2] = 10;				polys[200][3] = 10;				polys[200][4] = 9;				polys[200][5] = 6;				polys[200][6] = 6;				polys[200][7] = 3;				polys[200][8] = 10;
	polys[201][0] = 0;				polys[201][1] = 2;				polys[201][2] = 10;				polys[201][3] = 10;				polys[201][4] = 9;				polys[201][5] = 4;				polys[201][6] = 4;				polys[201][7] = 0;				polys[201][8] = 10;				polys[201][9] = 3;				polys[201][10] = 1;				polys[201][11] = 6;
	polys[202][0] = 3;				polys[202][1] = 0;				polys[202][2] = 5;				polys[202][3] = 5;				polys[202][4] = 10;				polys[202][5] = 9;				polys[202][6] = 9;				polys[202][7] = 6;				polys[202][8] = 3;				polys[202][9] = 3;				polys[202][10] = 5;				polys[202][11] = 9;
	polys[203][0] = 3;				polys[203][1] = 1;				polys[203][2] = 6;				polys[203][3] = 10;				polys[203][4] = 9;				polys[203][5] = 4;				polys[203][6] = 4;				polys[203][7] = 5;				polys[203][8] = 10;
	polys[204][0] = 1;				polys[204][1] = 2;				polys[204][2] = 10;				polys[204][3] = 10;				polys[204][4] = 9;				polys[204][5] = 1;
	polys[205][0] = 0;				polys[205][1] = 2;				polys[205][2] = 10;				polys[205][3] = 10;				polys[205][4] = 9;				polys[205][5] = 4;				polys[205][6] = 4;				polys[205][7] = 0;				polys[205][8] = 10;
	polys[206][0] = 1;				polys[206][1] = 0;				polys[206][2] = 5;				polys[206][3] = 5;				polys[206][4] = 10;				polys[206][5] = 9;				polys[206][6] = 9;				polys[206][7] = 1;				polys[206][8] = 5;
	polys[207][0] = 10;				polys[207][1] = 9;				polys[207][2] = 4;				polys[207][3] = 4;				polys[207][4] = 5;				polys[207][5] = 10;
	polys[208][0] = 10;				polys[208][1] = 8;				polys[208][2] = 4;				polys[208][3] = 4;				polys[208][4] = 6;				polys[208][5] = 7;				polys[208][6] = 7;				polys[208][7] = 10;				polys[208][8] = 4;
	polys[209][0] = 0;				polys[209][1] = 1;				polys[209][2] = 6;				polys[209][3] = 6;				polys[209][4] = 7;				polys[209][5] = 10;				polys[209][6] = 10;				polys[209][7] = 8;				polys[209][8] = 0;				polys[209][9] = 0;				polys[209][10] = 6;				polys[209][11] = 10;
	polys[210][0] = 2;				polys[210][1] = 0;				polys[210][2] = 4;				polys[210][3] = 4;				polys[210][4] = 6;				polys[210][5] = 7;				polys[210][6] = 7;				polys[210][7] = 2;				polys[210][8] = 4;				polys[210][9] = 10;				polys[210][10] = 8;				polys[210][11] = 5;
	polys[211][0] = 2;				polys[211][1] = 1;				polys[211][2] = 6;				polys[211][3] = 6;				polys[211][4] = 7;				polys[211][5] = 2;				polys[211][6] = 10;				polys[211][7] = 8;				polys[211][8] = 5;
	polys[212][0] = 1;				polys[212][1] = 3;				polys[212][2] = 7;				polys[212][3] = 7;				polys[212][4] = 10;				polys[212][5] = 8;				polys[212][6] = 8;				polys[212][7] = 4;				polys[212][8] = 1;				polys[212][9] = 1;				polys[212][10] = 7;				polys[212][11] = 8;
	polys[213][0] = 0;				polys[213][1] = 3;				polys[213][2] = 7;				polys[213][3] = 7;				polys[213][4] = 10;				polys[213][5] = 8;				polys[213][6] = 8;				polys[213][7] = 0;				polys[213][8] = 7;
	polys[214][0] = 2;				polys[214][1] = 3;				polys[214][2] = 7;				polys[214][3] = 1;				polys[214][4] = 0;				polys[214][5] = 4;				polys[214][6] = 10;				polys[214][7] = 8;				polys[214][8] = 5;
	polys[215][0] = 2;				polys[215][1] = 3;				polys[215][2] = 7;				polys[215][3] = 10;				polys[215][4] = 8;				polys[215][5] = 5;
	polys[216][0] = 3;				polys[216][1] = 2;				polys[216][2] = 10;				polys[216][3] = 10;				polys[216][4] = 8;				polys[216][5] = 4;				polys[216][6] = 4;				polys[216][7] = 6;				polys[216][8] = 3;				polys[216][9] = 3;				polys[216][10] = 10;			polys[216][11] = 4;
	polys[217][0] = 0;				polys[217][1] = 2;				polys[217][2] = 10;				polys[217][3] = 10;				polys[217][4] = 8;				polys[217][5] = 0;				polys[217][6] = 3;				polys[217][7] = 1;				polys[217][8] = 6;
	polys[218][0] = 3;				polys[218][1] = 0;				polys[218][2] = 4;				polys[218][3] = 4;				polys[218][4] = 6;				polys[218][5] = 3;				polys[218][6] = 10;				polys[218][7] = 8;				polys[218][8] = 5;
	polys[219][0] = 3;				polys[219][1] = 1;				polys[219][2] = 6;				polys[219][3] = 10;				polys[219][4] = 8;				polys[219][5] = 5;
	polys[220][0] = 1;				polys[220][1] = 2;				polys[220][2] = 10;				polys[220][3] = 10;				polys[220][4] = 8;				polys[220][5] = 4;				polys[220][6] = 4;				polys[220][7] = 1;				polys[220][8] = 10;
	polys[221][0] = 0;				polys[221][1] = 2;				polys[221][2] = 10;				polys[221][3] = 10;				polys[221][4] = 8;				polys[221][5] = 0;
	polys[222][0] = 1;				polys[222][1] = 0;				polys[222][2] = 4;				polys[222][3] = 10;				polys[222][4] = 8;				polys[222][5] = 5;
	polys[223][0] = 10;				polys[223][1] = 8;				polys[223][2] = 5;
	polys[224][0] = 8;				polys[224][1] = 9;				polys[224][2] = 6;				polys[224][3] = 6;				polys[224][4] = 7;				polys[224][5] = 5;				polys[224][6] = 5;				polys[224][7] = 8;				polys[224][8] = 6;
	polys[225][0] = 0;				polys[225][1] = 1;				polys[225][2] = 6;				polys[225][3] = 6;				polys[225][4] = 7;				polys[225][5] = 5;				polys[225][6] = 5;				polys[225][7] = 0;				polys[225][8] = 6;				polys[225][9] = 8;				polys[225][10] = 9;				polys[225][11] = 4;
	polys[226][0] = 2;				polys[226][1] = 0;				polys[226][2] = 8;				polys[226][3] = 8;				polys[226][4] = 9;				polys[226][5] = 6;				polys[226][6] = 6;				polys[226][7] = 7;				polys[226][8] = 2;				polys[226][9] = 2;				polys[226][10] = 8;				polys[226][11] = 6;
	polys[227][0] = 2;				polys[227][1] = 1;				polys[227][2] = 6;				polys[227][3] = 6;				polys[227][4] = 7;				polys[227][5] = 2;				polys[227][6] = 8;				polys[227][7] = 9;				polys[227][8] = 4;
	polys[228][0] = 1;				polys[228][1] = 3;				polys[228][2] = 7;				polys[228][3] = 7;				polys[228][4] = 5;				polys[228][5] = 8;				polys[228][6] = 8;				polys[228][7] = 9;				polys[228][8] = 1;				polys[228][9] = 1;				polys[228][10] = 7;				polys[228][11] = 8;
	polys[229][0] = 0;				polys[229][1] = 3;				polys[229][2] = 7;				polys[229][3] = 7;				polys[229][4] = 5;				polys[229][5] = 0;				polys[229][6] = 8;				polys[229][7] = 9;				polys[229][8] = 4;
	polys[230][0] = 2;				polys[230][1] = 3;				polys[230][2] = 7;				polys[230][3] = 1;				polys[230][4] = 0;				polys[230][5] = 8;				polys[230][6] = 8;				polys[230][7] = 9;				polys[230][8] = 1;
	polys[231][0] = 2;				polys[231][1] = 3;				polys[231][2] = 7;				polys[231][3] = 8;				polys[231][4] = 9;				polys[231][5] = 4;
	polys[232][0] = 3;				polys[232][1] = 2;				polys[232][2] = 5;				polys[232][3] = 5;				polys[232][4] = 8;				polys[232][5] = 9;				polys[232][6] = 9;				polys[232][7] = 6;				polys[232][8] = 3;				polys[232][9] = 3;				polys[232][10] = 5;				polys[232][11] = 9;
	polys[233][0] = 0;				polys[233][1] = 2;				polys[233][2] = 5;				polys[233][3] = 3;				polys[233][4] = 1;				polys[233][5] = 6;				polys[233][6] = 8;				polys[233][7] = 9;				polys[233][8] = 4;
	polys[234][0] = 3;				polys[234][1] = 0;				polys[234][2] = 8;				polys[234][3] = 8;				polys[234][4] = 9;				polys[234][5] = 6;				polys[234][6] = 6;				polys[234][7] = 3;				polys[234][8] = 8;
	polys[235][0] = 3;				polys[235][1] = 1;				polys[235][2] = 6;				polys[235][3] = 8;				polys[235][4] = 9;				polys[235][5] = 4;
	polys[236][0] = 1;				polys[236][1] = 2;				polys[236][2] = 5;				polys[236][3] = 5;				polys[236][4] = 8;				polys[236][5] = 9;				polys[236][6] = 9;				polys[236][7] = 1;				polys[236][8] = 5;
	polys[237][0] = 0;				polys[237][1] = 2;				polys[237][2] = 5;				polys[237][3] = 8;				polys[237][4] = 9;				polys[237][5] = 4;
	polys[238][0] = 1;				polys[238][1] = 0;				polys[238][2] = 8;				polys[238][3] = 8;				polys[238][4] = 9;				polys[238][5] = 1;
	polys[239][0] = 8;				polys[239][1] = 9;				polys[239][2] = 4;
	polys[240][0] = 5;				polys[240][1] = 4;				polys[240][2] = 6;				polys[240][3] = 6;				polys[240][4] = 7;				polys[240][5] = 5;
	polys[241][0] = 0;				polys[241][1] = 1;				polys[241][2] = 6;				polys[241][3] = 6;				polys[241][4] = 7;				polys[241][5] = 5;				polys[241][6] = 5;				polys[241][7] = 0;				polys[241][8] = 6;
	polys[242][0] = 2;				polys[242][1] = 0;				polys[242][2] = 4;				polys[242][3] = 4;				polys[242][4] = 6;				polys[242][5] = 7;				polys[242][6] = 7;				polys[242][7] = 2;				polys[242][8] = 4;
	polys[243][0] = 2;				polys[243][1] = 1;				polys[243][2] = 6;				polys[243][3] = 6;				polys[243][4] = 7;				polys[243][5] = 2;
	polys[244][0] = 1;				polys[244][1] = 3;				polys[244][2] = 7;				polys[244][3] = 7;				polys[244][4] = 5;				polys[244][5] = 4;				polys[244][6] = 4;				polys[244][7] = 1;				polys[244][8] = 7;
	polys[245][0] = 0;				polys[245][1] = 3;				polys[245][2] = 7;				polys[245][3] = 7;				polys[245][4] = 5;				polys[245][5] = 0;
	polys[246][0] = 2;				polys[246][1] = 3;				polys[246][2] = 7;				polys[246][3] = 1;				polys[246][4] = 0;				polys[246][5] = 4;
	polys[247][0] = 2;				polys[247][1] = 3;				polys[247][2] = 7;
	polys[248][0] = 3;				polys[248][1] = 2;				polys[248][2] = 5;				polys[248][3] = 5;				polys[248][4] = 4;				polys[248][5] = 6;				polys[248][6] = 6;				polys[248][7] = 3;				polys[248][8] = 5;
	polys[249][0] = 0;				polys[249][1] = 2;				polys[249][2] = 5;				polys[249][3] = 3;				polys[249][4] = 1;				polys[249][5] = 6;
	polys[250][0] = 3;				polys[250][1] = 0;				polys[250][2] = 4;				polys[250][3] = 4;				polys[250][4] = 6;				polys[250][5] = 3;
	polys[251][0] = 3;				polys[251][1] = 1;				polys[251][2] = 6;
	polys[252][0] = 1;				polys[252][1] = 2;				polys[252][2] = 5;				polys[252][3] = 5;				polys[252][4] = 4;				polys[252][5] = 1;
	polys[253][0] = 0;				polys[253][1] = 2;				polys[253][2] = 5;
	polys[254][0] = 1;				polys[254][1] = 0;				polys[254][2] = 4;
	
}




CubeTable::~CubeTable()
{
	delete[](sizes);
	int i = 0;
	for(i = 0; i<size; i++){
		delete[](polys[i]);
	}
	delete[](polys);

}

void CubeTable::toString()
{
	printf("Printing out CubeTable:\n\n");
	printf("Size: %i\n", size);
	
	int i = 0;
	int j = 0;
	
	printf("int[] sizes = {");
	for(i = 0; i<size; i++){
		printf("\tsizes[%i] = %i;\n", i, sizes[i]);
	}
	printf("\n\n");

	for(i = 0; i<size; i++){
		for(j = 0; j<3*sizes[i]; j++){
			printf("\tpolys[%i][%i] = %i;\n", i, j, polys[i][j] );
		}
	}
}



/****************************************************************************
 * VASP: Volumetric Analysis of Surface Properties
 * Copyright (c) 2014 Brian Y. Chen
 * All rights reserved.
 ***************************************************************************/


