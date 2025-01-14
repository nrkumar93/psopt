/*
 * Copyright (c) 2021, Ramkumar Natarajan
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the Carnegie Mellon University nor the names of its
 *       contributors may be used to endorse or promote products derived from
 *       this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */
/*!
 * \file   dynamics.cpp
 * \author Ramkumar Natarajan (rnataraj@cs.cmu.edu)
 * \date   6/22/21
 */

#include "dynamics.h"
#include "mcoder_math.h"

//void xdae(adouble* derivatives, adouble* path, adouble* states,
//         adouble* controls, adouble* parameters, adouble& time,
//         adouble* xad, int iphase, Workspace* workspace)
//{
//  adouble xdot, ydot, vdot;
//
//  adouble x1 = states[ 0 ];
//  adouble x2 = states[ 1 ];
//  adouble x3 = states[ 2 ];
//  adouble x4 = states[ 3 ];
//
//  adouble u1 = controls[ 0 ];
//  adouble u2 = controls[ 1 ];
//
//  adouble num1 =  sin(x3)*( (9.0/4.0)*cos(x3)*x1*x1+2*x2*x2 )
//                  + (4.0/3.0)*(u1-u2)-(3.0/2.0)*cos(x3)*u2;
//
//  adouble num2 =  -(sin(x3)*((7.0/2.0)*x1*x1+(9.0/4.0)*cos(x3)*x2*x2)
//                    -(7.0/3.0)*u2+(3.0/2.0)*cos(x3)*(u1-u2));
//
//  adouble den  =  31.0/36.0 + 9.0/4.0*pow(sin(x3),2);
//
//  derivatives[ 0 ] = num1/den;
//  derivatives[ 1 ] = num2/den;
//  derivatives[ 2 ] = x2 - x1;
//  derivatives[ 3 ] = x1;
//}

void dae(adouble* deriv, adouble* path, adouble* state,
         adouble* input, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{
  adouble xx[134];
  xx[0] = 5.0;
  xx[1] = 0.1259556365759687;
  xx[2] = xx[1] * state[1];
  xx[3] = 2.0;
  xx[4] = 0.992035875165179;
  xx[5] = 0.5;
  xx[6] = xx[5] * state[2];
  xx[7] = sin(xx[6]);
  xx[8] = xx[4] * xx[7];
  xx[9] = xx[1] * xx[7];
  xx[7] = xx[4] * state[1];
  xx[10] = xx[9] * xx[7] - xx[2] * xx[8];
  xx[11] = xx[2] + xx[3] * xx[8] * xx[10];
  xx[12] = xx[1] * state[3];
  xx[13] = xx[11] + xx[12];
  xx[14] = cos(xx[6]);
  xx[6] = xx[3] * xx[10] * xx[14];
  xx[15] = xx[7] - xx[3] * xx[9] * xx[10];
  xx[10] = xx[4] * state[3];
  xx[16] = xx[15] + xx[10];
  xx[17] = xx[13];
  xx[18] = xx[6];
  xx[19] = xx[16];
  xx[20] = 0.1521999311646401;
  xx[21] = 2.484885165017305;
  xx[22] = 2.495895028607315;
  xx[23] = xx[20] * xx[13];
  xx[24] = xx[21] * xx[6];
  xx[25] = xx[22] * xx[16];
  pm_math_Vector3_cross_ra(xx + 17, xx + 23, xx + 26);
  xx[23] = 0.9980069827323804;
  xx[24] = xx[5] * state[4];
  xx[25] = cos(xx[24]);
  xx[29] = xx[23] * xx[25];
  xx[30] = 0.06310358482218364;
  xx[31] = sin(xx[24]);
  xx[24] = xx[30] * xx[31];
  xx[32] = xx[30] * xx[25];
  xx[25] = xx[23] * xx[31];
  xx[33] = xx[29];
  xx[34] = xx[24];
  xx[35] = xx[32];
  xx[36] = xx[25];
  pm_math_Quaternion_inverseXform_ra(xx + 33, xx + 17, xx + 37);
  xx[23] = xx[39] + state[5];
  xx[40] = xx[37];
  xx[41] = xx[38];
  xx[42] = xx[23];
  xx[30] = 4.166666666666641e-3;
  xx[31] = 0.8354166666666666;
  xx[43] = xx[30] * xx[37];
  xx[44] = xx[31] * xx[38];
  xx[45] = xx[31] * xx[23];
  pm_math_Vector3_cross_ra(xx + 40, xx + 43, xx + 46);
  xx[40] = 10.0;
  xx[41] = 0.2233135416540104;
  xx[42] = xx[5] * xx[32];
  xx[43] = xx[5] * xx[25];
  xx[44] = 0.1142118812888779;
  xx[49] = xx[41] - (xx[3] * (xx[32] * xx[42] + xx[43] * xx[25]) - xx[5]);
  xx[50] = xx[3] * (xx[43] * xx[29] + xx[24] * xx[42]);
  xx[51] = xx[44] - xx[3] * (xx[42] * xx[29] - xx[24] * xx[43]);
  pm_math_Vector3_cross_ra(xx + 17, xx + 49, xx + 52);
  pm_math_Vector3_cross_ra(xx + 17, xx + 52, xx + 55);
  pm_math_Quaternion_inverseXform_ra(xx + 33, xx + 55, xx + 17);
  xx[42] = xx[40] * xx[18];
  xx[43] = 3.335416666666667;

  xx[45] = (input[2] - xx[0] * state[4] - (xx[48] + xx[5] * xx[42])) / xx[43];
  xx[52] = xx[46] + xx[30] * state[5] * xx[38];
  xx[53] = xx[47] - xx[31] * state[5] * xx[37];
  xx[54] = xx[48] + xx[31] * xx[45];
  pm_math_Quaternion_xform_ra(xx + 33, xx + 52, xx + 46);
  xx[38] = xx[5] * state[5];
  xx[52] = xx[40] * (xx[17] - xx[38] * (xx[39] + xx[23]));
  xx[53] = xx[42] + xx[0] * xx[45];
  xx[54] = xx[40] * (xx[38] * (xx[37] + xx[37]) + xx[19]);
  pm_math_Quaternion_xform_ra(xx + 33, xx + 52, xx + 17);
  pm_math_Vector3_cross_ra(xx + 49, xx + 17, xx + 37);
  xx[23] = xx[29] * xx[29];
  xx[42] = 1.0;
  xx[52] = xx[3] * (xx[23] + xx[24] * xx[24]) - xx[42];
  xx[53] = xx[32] * xx[24];
  xx[54] = xx[29] * xx[25];
  xx[55] = xx[3] * (xx[53] - xx[54]);
  xx[56] = xx[24] * xx[25];
  xx[57] = xx[32] * xx[29];
  xx[58] = xx[3] * (xx[56] + xx[57]);
  xx[59] = xx[3] * (xx[53] + xx[54]);
  xx[53] = xx[3] * (xx[23] + xx[32] * xx[32]) - xx[42];
  xx[54] = xx[32] * xx[25];
  xx[32] = xx[24] * xx[29];
  xx[24] = xx[3] * (xx[54] - xx[32]);
  xx[29] = xx[3] * (xx[56] - xx[57]);
  xx[56] = xx[3] * (xx[54] + xx[32]);
  xx[32] = xx[3] * (xx[23] + xx[25] * xx[25]) - xx[42];
  xx[60] = xx[52];
  xx[61] = xx[55];
  xx[62] = xx[58];
  xx[63] = xx[59];
  xx[64] = xx[53];
  xx[65] = xx[24];
  xx[66] = xx[29];
  xx[67] = xx[56];
  xx[68] = xx[32];
  xx[23] = xx[31] / xx[43];
  xx[25] = xx[31] - xx[31] * xx[23];
  xx[69] = xx[30] * xx[52];
  xx[70] = xx[30] * xx[59];
  xx[71] = xx[30] * xx[29];
  xx[72] = xx[31] * xx[55];
  xx[73] = xx[31] * xx[53];
  xx[74] = xx[31] * xx[56];
  xx[75] = xx[58] * xx[25];
  xx[76] = xx[24] * xx[25];
  xx[77] = xx[25] * xx[32];
  pm_math_Matrix3x3_compose_ra(xx + 60, xx + 69, xx + 78);
  xx[25] = xx[0] / xx[43];
  xx[30] = xx[31] * xx[25];
  xx[31] = xx[30] * xx[55];
  xx[43] = xx[58] * xx[31];
  xx[54] = xx[30] * xx[53];
  xx[57] = xx[58] * xx[54];
  xx[69] = xx[30] * xx[56];
  xx[30] = xx[58] * xx[69];
  xx[70] = xx[24] * xx[31];
  xx[71] = xx[24] * xx[54];
  xx[72] = xx[24] * xx[69];
  xx[73] = xx[31] * xx[32];
  xx[31] = xx[54] * xx[32];
  xx[54] = xx[69] * xx[32];
  xx[87] = - xx[43];
  xx[88] = - xx[57];
  xx[89] = - xx[30];
  xx[90] = - xx[70];
  xx[91] = - xx[71];
  xx[92] = - xx[72];
  xx[93] = - xx[73];
  xx[94] = - xx[31];
  xx[95] = - xx[54];
  pm_math_Matrix3x3_postCross_ra(xx + 87, xx + 49, xx + 96);
  xx[69] = xx[40] - xx[0] * xx[25];
  xx[87] = xx[40] * xx[52];
  xx[88] = xx[40] * xx[59];
  xx[89] = xx[40] * xx[29];
  xx[90] = xx[55] * xx[69];
  xx[91] = xx[69] * xx[53];
  xx[92] = xx[56] * xx[69];
  xx[93] = xx[40] * xx[58];
  xx[94] = xx[40] * xx[24];
  xx[95] = xx[40] * xx[32];
  pm_math_Matrix3x3_compose_ra(xx + 60, xx + 87, xx + 105);
  pm_math_Matrix3x3_postCross_ra(xx + 105, xx + 49, xx + 58);
  pm_math_Matrix3x3_preCross_ra(xx + 58, xx + 49, xx + 87);
  xx[24] = xx[20] + xx[78] - xx[96] - xx[96] - xx[87];
  xx[29] = xx[79] - xx[97] - xx[99] - xx[88];
  xx[32] = xx[80] - xx[98] - xx[102] - xx[89];
  xx[40] = xx[81] - xx[99] - xx[97] - xx[90];
  xx[52] = xx[82] - xx[100] - xx[100] - xx[91];
  xx[53] = xx[83] - xx[101] - xx[103] - xx[92];
  xx[55] = xx[84] - xx[102] - xx[98] - xx[93];
  xx[56] = xx[85] - xx[103] - xx[101] - xx[94];
  xx[67] = xx[22] + xx[86] - xx[104] - xx[104] - xx[95];
  xx[74] = xx[24];
  xx[75] = xx[29];
  xx[76] = xx[32];
  xx[77] = xx[40];
  xx[78] = xx[21] + xx[52];
  xx[79] = xx[53];
  xx[80] = xx[55];
  xx[81] = xx[56];
  xx[82] = xx[67];
  xx[68] = xx[10] * xx[6];
  xx[69] = xx[12] * xx[15] - xx[10] * xx[11];
  xx[10] = xx[12] * xx[6];
  xx[83] = xx[68];
  xx[84] = xx[69];
  xx[85] = - xx[10];
  pm_math_Matrix3x3_xform_ra(xx + 74, xx + 83, xx + 86);
  xx[6] = xx[43] + xx[58];
  xx[12] = xx[57] + xx[61];
  xx[43] = xx[30] + xx[64];
  xx[30] = xx[70] + xx[59];
  xx[57] = xx[71] + xx[62];
  xx[58] = xx[72] + xx[65];
  xx[59] = xx[73] + xx[60];
  xx[60] = xx[31] + xx[63];
  xx[31] = xx[54] + xx[66];
  xx[70] = - xx[6];
  xx[71] = - xx[12];
  xx[72] = - xx[43];
  xx[73] = - xx[30];
  xx[74] = - xx[57];
  xx[75] = - xx[58];
  xx[76] = - xx[59];
  xx[77] = - xx[60];
  xx[78] = - xx[31];
  xx[54] = 0.04176034283181076;
  xx[61] = 0.7939134608263623;
  xx[62] = xx[54] * xx[9] + xx[61] * xx[8];
  xx[63] = xx[41] - (xx[3] * xx[8] * xx[62] - xx[61]);
  xx[41] = xx[44] - (xx[54] - xx[3] * xx[9] * xx[62]);
  xx[44] = xx[7] * xx[63] - xx[2] * xx[41];
  xx[54] = xx[7] * xx[44];
  xx[61] = xx[2] * xx[44];
  xx[44] = xx[8] * xx[54] + xx[9] * xx[61];
  xx[64] = xx[3] * xx[62] * xx[14];
  xx[62] = xx[7] * xx[7] * xx[64] + xx[2] * xx[2] * xx[64];
  xx[65] = xx[8] * xx[62];
  xx[66] = 0.7928505854813077;
  xx[79] = xx[66] * state[3];
  xx[80] = xx[3] * (xx[8] * xx[44] - xx[65] * xx[14]) - xx[54] - xx[79] * (xx[15]
    + xx[16]);
  xx[15] = xx[9] * xx[62];
  xx[16] = xx[3] * (xx[44] * xx[14] + xx[8] * xx[65] + xx[9] * xx[15]) - xx[62];
  xx[54] = xx[79] * (xx[11] + xx[13]) + xx[61] + xx[3] * (xx[15] * xx[14] - xx[9]
    * xx[44]);
  xx[89] = xx[80];
  xx[90] = xx[16];
  xx[91] = xx[54];
  pm_math_Matrix3x3_xform_ra(xx + 70, xx + 89, xx + 92);
  xx[11] = xx[26] + xx[46] + xx[37] + xx[86] + xx[92];
  xx[13] = xx[1] * xx[24] + xx[4] * xx[32] - xx[66] * xx[12];
  xx[15] = xx[28] + xx[48] + xx[39] + xx[88] + xx[94];
  pm_math_Matrix3x3_transposeXform_ra(xx + 70, xx + 83, xx + 95);
  xx[44] = 24.13716694115407;
  xx[61] = xx[44] + xx[109];
  xx[70] = xx[44] + xx[105];
  xx[71] = xx[106];
  xx[72] = xx[107];
  xx[73] = xx[108];
  xx[74] = xx[61];
  xx[75] = xx[110];
  xx[76] = xx[111];
  xx[77] = xx[112];
  xx[78] = xx[44] + xx[113];
  pm_math_Matrix3x3_xform_ra(xx + 70, xx + 89, xx + 81);
  xx[62] = xx[18] + xx[96] + xx[82];
  xx[18] = xx[1] * xx[55] + xx[4] * xx[67] - xx[66] * xx[60];
  xx[65] = xx[66] * xx[61] - (xx[1] * xx[12] + xx[4] * xx[60]);
  xx[70] = xx[1] * xx[13] + xx[4] * xx[18] + xx[66] * xx[65];

  xx[71] = (input[1] - xx[0] * state[2] - (xx[1] * xx[11] + xx[4] * xx[15] + xx
             [66] * xx[62])) / xx[70];
  xx[72] = xx[11] + xx[13] * xx[71];
  xx[11] = xx[15] + xx[18] * xx[71];
  xx[15] = xx[9] * xx[11] - xx[8] * xx[72];
  xx[26] = xx[1] * xx[40] + xx[4] * xx[53] - xx[66] * xx[57];
  xx[28] = xx[27] + xx[47] + xx[38] + xx[87] + xx[93] + xx[26] * xx[71];
  xx[37] = xx[63];
  xx[38] = xx[64];
  xx[39] = xx[41];
  xx[27] = xx[66] * xx[106] - (xx[1] * xx[6] + xx[4] * xx[59]);
  xx[46] = xx[17] + xx[95] + xx[81] + xx[27] * xx[71];
  xx[17] = xx[66] * xx[112] - (xx[1] * xx[43] + xx[4] * xx[31]);
  xx[47] = xx[19] + xx[97] + xx[83] + xx[17] * xx[71];
  xx[19] = xx[9] * xx[47] - xx[8] * xx[46];
  xx[48] = xx[62] + xx[65] * xx[71];
  xx[62] = xx[8] * xx[48];
  xx[73] = xx[9] * xx[48];
  xx[74] = xx[48] - xx[3] * (xx[19] * xx[14] + xx[8] * xx[62] + xx[9] * xx[73]);
  xx[75] = xx[46] + xx[3] * (xx[8] * xx[19] - xx[62] * xx[14]);
  xx[76] = xx[74];
  xx[77] = xx[47] + xx[3] * (xx[73] * xx[14] - xx[9] * xx[19]);
  pm_math_Vector3_cross_ra(xx + 37, xx + 75, xx + 46);
  xx[19] = xx[66] * state[1];
  xx[62] = xx[2] * xx[19];
  xx[2] = xx[14] * xx[14];
  xx[73] = xx[3] * xx[8] * xx[14];
  xx[75] = xx[3] * xx[9] * xx[8];
  xx[76] = xx[3] * xx[9] * xx[14];
  xx[81] = xx[3] * (xx[2] + xx[9] * xx[9]) - xx[42];
  xx[82] = - xx[73];
  xx[83] = xx[75];
  xx[84] = xx[73];
  xx[85] = xx[3] * xx[2] - xx[42];
  xx[86] = - xx[76];
  xx[87] = xx[75];
  xx[88] = xx[76];
  xx[89] = xx[3] * (xx[2] + xx[8] * xx[8]) - xx[42];
  xx[2] = xx[27] / xx[70];
  xx[42] = xx[65] / xx[70];
  xx[73] = xx[17] / xx[70];
  xx[90] = - (xx[6] + xx[13] * xx[2]);
  xx[91] = - (xx[12] + xx[13] * xx[42]);
  xx[92] = - (xx[43] + xx[13] * xx[73]);
  xx[93] = - (xx[30] + xx[26] * xx[2]);
  xx[94] = - (xx[57] + xx[26] * xx[42]);
  xx[95] = - (xx[58] + xx[26] * xx[73]);
  xx[96] = - (xx[59] + xx[18] * xx[2]);
  xx[97] = - (xx[60] + xx[18] * xx[42]);
  xx[98] = - (xx[31] + xx[18] * xx[73]);
  pm_math_Matrix3x3_composeTranspose_ra(xx + 90, xx + 81, xx + 114);
  pm_math_Matrix3x3_compose_ra(xx + 81, xx + 114, xx + 90);
  xx[6] = xx[65] * xx[2];
  xx[12] = xx[17] * xx[2];
  xx[30] = xx[17] * xx[42];
  xx[114] = xx[105] - xx[27] * xx[2] + xx[44];
  xx[115] = xx[106] - xx[6];
  xx[116] = xx[107] - xx[12];
  xx[117] = xx[108] - xx[6];
  xx[118] = xx[61] - xx[65] * xx[42];
  xx[119] = xx[110] - xx[30];
  xx[120] = xx[111] - xx[12];
  xx[121] = xx[112] - xx[30];
  xx[122] = xx[113] - xx[17] * xx[73] + xx[44];
  pm_math_Matrix3x3_composeTranspose_ra(xx + 114, xx + 81, xx + 99);
  pm_math_Matrix3x3_compose_ra(xx + 81, xx + 99, xx + 108);
  pm_math_Matrix3x3_postCross_ra(xx + 108, xx + 37, xx + 99);
  xx[6] = xx[92] - xx[105];
  xx[12] = xx[19] * xx[7];
  xx[7] = xx[90] - xx[99];
  xx[17] = xx[98] - xx[107];
  xx[19] = xx[96] - xx[101];
  xx[27] = xx[13] / xx[70];
  xx[30] = xx[26] * xx[27];
  xx[31] = xx[18] * xx[27];
  xx[43] = xx[26] / xx[70];
  xx[47] = xx[18] * xx[43];
  xx[57] = xx[18] / xx[70];
  xx[116] = xx[24] - xx[13] * xx[27];
  xx[117] = xx[29] - xx[30];
  xx[118] = xx[32] - xx[31];
  xx[119] = xx[40] - xx[30];
  xx[120] = xx[52] - xx[26] * xx[43] + xx[21];
  xx[121] = xx[53] - xx[47];
  xx[122] = xx[55] - xx[31];
  xx[123] = xx[56] - xx[47];
  xx[124] = xx[67] - xx[18] * xx[57];
  pm_math_Matrix3x3_composeTranspose_ra(xx + 116, xx + 81, xx + 125);
  pm_math_Matrix3x3_compose_ra(xx + 81, xx + 125, xx + 116);
  pm_math_Matrix3x3_postCross_ra(xx + 90, xx + 37, xx + 81);
  pm_math_Matrix3x3_preCross_ra(xx + 99, xx + 37, xx + 125);
  xx[13] = xx[91] - xx[102];
  xx[18] = xx[97] - xx[104];
  xx[21] = xx[1] * xx[13] + xx[4] * xx[18] + xx[66] * (xx[44] + xx[112]);
  xx[24] = xx[1] * (xx[1] * (xx[20] + xx[116] - xx[81] - xx[81] - xx[125]) + xx
                    [4] * (xx[118] - xx[83] - xx[87] - xx[127]) + xx[66] * xx[13])
    + xx[4] * (xx[1] * (xx[122] - xx[87] - xx[83] - xx[131]) + xx[4] * (xx[22] +
    xx[124] - xx[89] - xx[89] - xx[133]) + xx[66] * xx[18]) + xx[66] * xx[21];

  xx[29] = (xx[1] * xx[7] + xx[4] * xx[19] + xx[66] * xx[109]) / xx[24];
  xx[30] = xx[21] / xx[24];
  xx[31] = (xx[1] * xx[6] + xx[4] * xx[17] + xx[66] * xx[115]) / xx[24];
  xx[13] = 9.81;
  xx[18] = 0.5305552837772819;
  xx[20] = xx[5] * state[0];
  xx[5] = sin(xx[20]);
  xx[21] = xx[1] * xx[5];
  xx[22] = xx[18] * xx[21];
  xx[26] = 0.4674516989550983;
  xx[32] = cos(xx[20]);
  xx[20] = xx[26] * xx[32];
  xx[37] = xx[4] * xx[5];
  xx[5] = xx[26] * xx[37];
  xx[38] = xx[22] - xx[20] + xx[5];
  xx[39] = xx[13] * xx[38];
  xx[40] = xx[18] * xx[32];
  xx[32] = xx[18] * xx[37];
  xx[18] = xx[26] * xx[21];
  xx[21] = xx[40] + xx[32] - xx[18];
  xx[26] = xx[22] + xx[20] + xx[5];
  xx[5] = xx[13] * xx[26];
  xx[20] = xx[40] - xx[32] + xx[18];
  xx[18] = xx[3] * (xx[39] * xx[21] - xx[5] * xx[20]);
  xx[22] = xx[3] * (xx[5] * xx[21] + xx[39] * xx[20]);
  xx[20] = xx[13] - xx[3] * (xx[5] * xx[26] + xx[39] * xx[38]);
  xx[37] = xx[18];
  xx[38] = xx[22];
  xx[39] = xx[20];
  xx[5] = (input[0] - xx[0] * state[0] - (xx[1] * (xx[72] + xx[3] * (xx[8] * xx
              [15] - xx[8] * xx[28] * xx[14]) + xx[46] + xx[62] * xx[6] - xx[12]
             * xx[7]) + xx[4] * (xx[11] + xx[3] * (xx[9] * xx[28] * xx[14] - xx
              [9] * xx[15]) + xx[48] + xx[62] * xx[17] - xx[12] * xx[19]) + xx
            [66] * (xx[74] + xx[62] * xx[113] - xx[12] * xx[111]))) / xx[24] -
    pm_math_Vector3_dot_ra(xx + 29, xx + 37);
  xx[28] = xx[27];
  xx[29] = xx[43];
  xx[30] = xx[57];
  xx[0] = xx[1] * xx[5];
  xx[6] = xx[4] * xx[5];
  xx[7] = xx[9] * xx[6] - xx[0] * xx[8];
  xx[11] = xx[0] + xx[3] * xx[8] * xx[7];
  xx[13] = xx[3] * xx[7] * xx[14];
  xx[15] = xx[6] - xx[3] * xx[9] * xx[7];
  xx[37] = xx[11];
  xx[38] = xx[13];
  xx[39] = xx[15];
  xx[46] = xx[2];
  xx[47] = xx[42];
  xx[48] = xx[73];
  xx[2] = xx[18] - xx[12] - xx[6] * xx[64];
  xx[7] = xx[22] + xx[66] * xx[5] + xx[6] * xx[63] - xx[0] * xx[41];
  xx[6] = xx[8] * xx[7];
  xx[12] = xx[20] + xx[62] + xx[0] * xx[64];
  xx[0] = xx[9] * xx[12] - xx[8] * xx[2];
  xx[17] = xx[2] + xx[3] * (xx[6] * xx[14] + xx[8] * xx[0]);
  xx[2] = xx[9] * xx[7];
  xx[18] = xx[7] + xx[3] * (xx[0] * xx[14] - (xx[8] * xx[6] + xx[9] * xx[2]));
  xx[6] = xx[12] - xx[3] * (xx[2] * xx[14] + xx[9] * xx[0]);
  xx[7] = xx[17];
  xx[8] = xx[18];
  xx[9] = xx[6];
  xx[0] = xx[71] - (pm_math_Vector3_dot_ra(xx + 28, xx + 37) +
                    pm_math_Vector3_dot_ra(xx + 46, xx + 7));
  xx[7] = xx[11] + xx[1] * xx[0] + xx[68];
  xx[8] = xx[13] + xx[69];
  xx[9] = xx[15] + xx[4] * xx[0] - xx[10];
  pm_math_Quaternion_inverseXform_ra(xx + 33, xx + 7, xx + 1);
  pm_math_Vector3_cross_ra(xx + 7, xx + 49, xx + 10);
  xx[7] = xx[17] + xx[80] + xx[10];
  xx[8] = xx[18] + xx[66] * xx[0] + xx[16] + xx[11];
  xx[9] = xx[6] + xx[54] + xx[12];
  pm_math_Quaternion_inverseXform_ra(xx + 33, xx + 7, xx + 10);
  deriv[0] = state[1];
  deriv[1] = xx[5];
  deriv[2] = state[3];
  deriv[3] = xx[0];
  deriv[4] = state[5];
  deriv[5] = xx[45] - (xx[23] * xx[3] + xx[25] * xx[11]);

//  adouble xx[135];
//  xx[0] = 0.1337606471309296;
//  xx[1] = xx[0] * state[1];
//  xx[2] = 2.0;
//  xx[3] = 0.9910136675541438;
//  xx[4] = 0.5;
//  xx[5] = xx[4] * state[2];
//  xx[6] = sin(xx[5]);
//  xx[7] = xx[3] * xx[6];
//  xx[8] = xx[0] * xx[6];
//  xx[6] = xx[3] * state[1];
//  xx[9] = xx[8] * xx[6] - xx[1] * xx[7];
//  xx[10] = xx[1] + xx[2] * xx[7] * xx[9];
//  xx[11] = xx[0] * state[3];
//  xx[12] = xx[10] + xx[11];
//  xx[13] = cos(xx[5]);
//  xx[5] = xx[2] * xx[9] * xx[13];
//  xx[14] = xx[6] - xx[2] * xx[8] * xx[9];
//  xx[9] = xx[3] * state[3];
//  xx[15] = xx[14] + xx[9];
//  xx[16] = xx[12];
//  xx[17] = xx[5];
//  xx[18] = xx[15];
//  xx[19] = 0.1300986516517284;
//  xx[20] = 0.4530108549021385;
//  xx[21] = 0.48237199800506;
//  xx[22] = xx[19] * xx[12];
//  xx[23] = xx[20] * xx[5];
//  xx[24] = xx[21] * xx[15];
//  pm_math_Vector3_cross_ra<adouble>(xx + 16, xx + 22, xx + 25);
//  xx[22] = 0.9977508876353215;
//  xx[23] = xx[4] * state[4];
//  xx[24] = cos(xx[23]);
//  xx[28] = xx[22] * xx[24];
//  xx[29] = 0.06703108400531897;
//  xx[30] = sin(xx[23]);
//  xx[23] = xx[29] * xx[30];
//  xx[31] = xx[29] * xx[24];
//  xx[24] = xx[22] * xx[30];
//  xx[32] = xx[28];
//  xx[33] = xx[23];
//  xx[34] = xx[31];
//  xx[35] = xx[24];
//  pm_math_Quaternion_inverseXform_ra<adouble>(xx + 32, xx + 16, xx + 36);
//  xx[22] = xx[38] + state[5];
//  xx[39] = xx[36];
//  xx[40] = xx[37];
//  xx[41] = xx[22];
//  xx[29] = 4.166666666666641e-4;
//  xx[30] = 0.08354166666666668;
//  xx[42] = xx[29] * xx[36];
//  xx[43] = xx[30] * xx[37];
//  xx[44] = xx[30] * xx[22];
//  pm_math_Vector3_cross_ra<adouble>(xx + 39, xx + 42, xx + 45);
//  xx[39] = 0.04699417175177818;
//  xx[40] = xx[4] * xx[31];
//  xx[41] = xx[4] * xx[24];
//  xx[42] = 0.1012299713559456;
//  xx[48] = xx[39] - (xx[2] * (xx[31] * xx[40] + xx[41] * xx[24]) - xx[4]);
//  xx[49] = xx[2] * (xx[41] * xx[28] + xx[23] * xx[40]);
//  xx[50] = xx[42] - xx[2] * (xx[40] * xx[28] - xx[23] * xx[41]);
//  pm_math_Vector3_cross_ra<adouble>(xx + 16, xx + 48, xx + 51);
//  pm_math_Vector3_cross_ra<adouble>(xx + 16, xx + 51, xx + 54);
//  pm_math_Quaternion_inverseXform_ra<adouble>(xx + 32, xx + 54, xx + 16);
//  xx[40] = 0.3335416666666667;
////    ii[0] = factorSymmetricPosDef(xx + 40, 1, xx + 41);
////    if (ii[0] != 0) {
////      return sm_ssci_recordRunTimeError(
////              "sm:compiler:messages:simulationErrors:DegenerateMass",
////              "'manipulator_ip/Revolute BaseArm3' has a degenerate mass distribution on its follower side.",
////              neDiagMgr);
////    }
//
//  xx[41] = (input[2] - (xx[47] + xx[4] * xx[17])) / xx[40];
//  xx[51] = xx[45] + xx[29] * state[5] * xx[37];
//  xx[52] = xx[46] - xx[30] * state[5] * xx[36];
//  xx[53] = xx[47] + xx[30] * xx[41];
//  pm_math_Quaternion_xform_ra<adouble>(xx + 32, xx + 51, xx + 43);
//  xx[37] = xx[4] * state[5];
//  xx[51] = xx[16] - xx[37] * (xx[38] + xx[22]);
//  xx[52] = xx[17] + xx[4] * xx[41];
//  xx[53] = xx[37] * (xx[36] + xx[36]) + xx[18];
//  pm_math_Quaternion_xform_ra<adouble>(xx + 32, xx + 51, xx + 16);
//  pm_math_Vector3_cross_ra<adouble>(xx + 48, xx + 16, xx + 36);
//  xx[22] = xx[28] * xx[28];
//  xx[46] = 1.0;
//  xx[47] = xx[2] * (xx[22] + xx[23] * xx[23]) - xx[46];
//  xx[51] = xx[31] * xx[23];
//  xx[52] = xx[28] * xx[24];
//  xx[53] = xx[2] * (xx[51] - xx[52]);
//  xx[54] = xx[23] * xx[24];
//  xx[55] = xx[31] * xx[28];
//  xx[56] = xx[2] * (xx[54] + xx[55]);
//  xx[57] = xx[51] + xx[52];
//  xx[51] = xx[2] * xx[57];
//  xx[52] = xx[2] * (xx[22] + xx[31] * xx[31]) - xx[46];
//  xx[58] = xx[31] * xx[24];
//  xx[31] = xx[23] * xx[28];
//  xx[23] = xx[2] * (xx[58] - xx[31]);
//  xx[28] = xx[54] - xx[55];
//  xx[54] = xx[2] * xx[28];
//  xx[55] = xx[2] * (xx[58] + xx[31]);
//  xx[31] = xx[2] * (xx[22] + xx[24] * xx[24]) - xx[46];
//  xx[58] = xx[47];
//  xx[59] = xx[53];
//  xx[60] = xx[56];
//  xx[61] = xx[51];
//  xx[62] = xx[52];
//  xx[63] = xx[23];
//  xx[64] = xx[54];
//  xx[65] = xx[55];
//  xx[66] = xx[31];
//  xx[22] = 8.333333333333282e-4;
//  xx[24] = xx[30] / xx[40];
//  xx[67] = xx[30] - xx[30] * xx[24];
//  xx[68] = xx[29] * xx[47];
//  xx[69] = xx[22] * xx[57];
//  xx[70] = xx[22] * xx[28];
//  xx[71] = xx[30] * xx[53];
//  xx[72] = xx[30] * xx[52];
//  xx[73] = xx[30] * xx[55];
//  xx[74] = xx[56] * xx[67];
//  xx[75] = xx[23] * xx[67];
//  xx[76] = xx[67] * xx[31];
//  pm_math_Matrix3x3_compose_ra<adouble>(xx + 58, xx + 68, xx + 77);
//  xx[22] = xx[4] / xx[40];
//  xx[28] = xx[30] * xx[22];
//  xx[29] = xx[28] * xx[53];
//  xx[30] = xx[56] * xx[29];
//  xx[40] = xx[28] * xx[52];
//  xx[57] = xx[56] * xx[40];
//  xx[67] = xx[28] * xx[55];
//  xx[28] = xx[56] * xx[67];
//  xx[68] = xx[23] * xx[29];
//  xx[69] = xx[23] * xx[40];
//  xx[70] = xx[23] * xx[67];
//  xx[71] = xx[29] * xx[31];
//  xx[29] = xx[40] * xx[31];
//  xx[40] = xx[67] * xx[31];
//  xx[86] = - xx[30];
//  xx[87] = - xx[57];
//  xx[88] = - xx[28];
//  xx[89] = - xx[68];
//  xx[90] = - xx[69];
//  xx[91] = - xx[70];
//  xx[92] = - xx[71];
//  xx[93] = - xx[29];
//  xx[94] = - xx[40];
//  pm_math_Matrix3x3_postCross_ra<adouble>(xx + 86, xx + 48, xx + 95);
//  xx[67] = xx[46] - xx[4] * xx[22];
//  xx[86] = xx[47];
//  xx[87] = xx[51];
//  xx[88] = xx[54];
//  xx[89] = xx[53] * xx[67];
//  xx[90] = xx[67] * xx[52];
//  xx[91] = xx[55] * xx[67];
//  xx[92] = xx[56];
//  xx[93] = xx[23];
//  xx[94] = xx[31];
//  pm_math_Matrix3x3_compose_ra<adouble>(xx + 58, xx + 86, xx + 104);
//  pm_math_Matrix3x3_postCross_ra<adouble>(xx + 104, xx + 48, xx + 58);
//  pm_math_Matrix3x3_preCross_ra<adouble>(xx + 58, xx + 48, xx + 86);
//  xx[23] = xx[19] + xx[77] - xx[95] - xx[95] - xx[86];
//  xx[31] = xx[78] - xx[96] - xx[98] - xx[87];
//  xx[47] = xx[79] - xx[97] - xx[101] - xx[88];
//  xx[51] = xx[80] - xx[98] - xx[96] - xx[89];
//  xx[52] = xx[81] - xx[99] - xx[99] - xx[90];
//  xx[53] = xx[82] - xx[100] - xx[102] - xx[91];
//  xx[54] = xx[83] - xx[101] - xx[97] - xx[92];
//  xx[55] = xx[84] - xx[102] - xx[100] - xx[93];
//  xx[56] = xx[21] + xx[85] - xx[103] - xx[103] - xx[94];
//  xx[72] = xx[23];
//  xx[73] = xx[31];
//  xx[74] = xx[47];
//  xx[75] = xx[51];
//  xx[76] = xx[20] + xx[52];
//  xx[77] = xx[53];
//  xx[78] = xx[54];
//  xx[79] = xx[55];
//  xx[80] = xx[56];
//  xx[67] = xx[9] * xx[5];
//  xx[81] = xx[11] * xx[14] - xx[9] * xx[10];
//  xx[9] = xx[11] * xx[5];
//  xx[82] = xx[67];
//  xx[83] = xx[81];
//  xx[84] = - xx[9];
//  pm_math_Matrix3x3_xform_ra<adouble>(xx + 72, xx + 82, xx + 85);
//  xx[5] = xx[30] + xx[58];
//  xx[11] = xx[57] + xx[61];
//  xx[30] = xx[28] + xx[64];
//  xx[28] = xx[68] + xx[59];
//  xx[57] = xx[69] + xx[62];
//  xx[58] = xx[70] + xx[65];
//  xx[59] = xx[71] + xx[60];
//  xx[60] = xx[29] + xx[63];
//  xx[29] = xx[40] + xx[66];
//  xx[68] = - xx[5];
//  xx[69] = - xx[11];
//  xx[70] = - xx[30];
//  xx[71] = - xx[28];
//  xx[72] = - xx[57];
//  xx[73] = - xx[58];
//  xx[74] = - xx[59];
//  xx[75] = - xx[60];
//  xx[76] = - xx[29];
//  xx[40] = 0.03678788497604644;
//  xx[61] = 0.9707716252285514;
//  xx[62] = xx[40] * xx[8] + xx[61] * xx[7];
//  xx[63] = xx[39] - (xx[2] * xx[7] * xx[62] - xx[61]);
//  xx[39] = xx[42] - (xx[40] - xx[2] * xx[8] * xx[62]);
//  xx[40] = xx[6] * xx[63] - xx[1] * xx[39];
//  xx[42] = xx[6] * xx[40];
//  xx[61] = xx[1] * xx[40];
//  xx[40] = xx[7] * xx[42] + xx[8] * xx[61];
//  xx[64] = xx[2] * xx[62] * xx[13];
//  xx[62] = xx[6] * xx[6] * xx[64] + xx[1] * xx[1] * xx[64];
//  xx[65] = xx[7] * xx[62];
//  xx[66] = 0.9669687199762177;
//  xx[77] = xx[66] * state[3];
//  xx[78] = xx[2] * (xx[7] * xx[40] - xx[65] * xx[13]) - xx[42] - xx[77] * (xx[14]
//                                                                           + xx[15]);
//  xx[14] = xx[8] * xx[62];
//  xx[15] = xx[2] * (xx[40] * xx[13] + xx[7] * xx[65] + xx[8] * xx[14]) - xx[62];
//  xx[42] = xx[77] * (xx[10] + xx[12]) + xx[61] + xx[2] * (xx[14] * xx[13] - xx[8]
//                                                                            * xx[40]);
//  xx[88] = xx[78];
//  xx[89] = xx[15];
//  xx[90] = xx[42];
//  pm_math_Matrix3x3_xform_ra<adouble>(xx + 68, xx + 88, xx + 91);
//  xx[10] = xx[25] + xx[43] + xx[36] + xx[85] + xx[91];
//  xx[12] = xx[0] * xx[23] + xx[3] * xx[47] - xx[66] * xx[11];
//  xx[14] = xx[27] + xx[45] + xx[38] + xx[87] + xx[93];
//  pm_math_Matrix3x3_transposeXform_ra<adouble>(xx + 68, xx + 82, xx + 94);
//  xx[40] = 15.13716694115407;
//  xx[61] = xx[40] + xx[108];
//  xx[68] = xx[40] + xx[104];
//  xx[69] = xx[105];
//  xx[70] = xx[106];
//  xx[71] = xx[107];
//  xx[72] = xx[61];
//  xx[73] = xx[109];
//  xx[74] = xx[110];
//  xx[75] = xx[111];
//  xx[76] = xx[40] + xx[112];
//  pm_math_Matrix3x3_xform_ra<adouble>(xx + 68, xx + 88, xx + 82);
//  xx[62] = xx[17] + xx[95] + xx[83];
//  xx[17] = xx[0] * xx[54] + xx[3] * xx[56] - xx[66] * xx[60];
//  xx[65] = xx[66] * xx[61] - (xx[0] * xx[11] + xx[3] * xx[60]);
//  xx[68] = xx[0] * xx[12] + xx[3] * xx[17] + xx[66] * xx[65];
////    ii[0] = factorSymmetricPosDef(xx + 68, 1, xx + 69);
////    if (ii[0] != 0) {
////      return sm_ssci_recordRunTimeError(
////              "sm:compiler:messages:simulationErrors:DegenerateMass",
////              "'manipulator_ip/Revolute BaseArm2' has a degenerate mass distribution on its follower side.",
////              neDiagMgr);
////    }
//
//  xx[69] = (input[1] - (xx[0] * xx[10] + xx[3] * xx[14] + xx[66] * xx[62])) /
//           xx[68];
//  xx[70] = xx[10] + xx[12] * xx[69];
//  xx[10] = xx[14] + xx[17] * xx[69];
//  xx[14] = xx[8] * xx[10] - xx[7] * xx[70];
//  xx[25] = xx[0] * xx[51] + xx[3] * xx[53] - xx[66] * xx[57];
//  xx[27] = xx[26] + xx[44] + xx[37] + xx[86] + xx[92] + xx[25] * xx[69];
//  xx[36] = xx[63];
//  xx[37] = xx[64];
//  xx[38] = xx[39];
//  xx[26] = xx[66] * xx[105] - (xx[0] * xx[5] + xx[3] * xx[59]);
//  xx[43] = xx[16] + xx[94] + xx[82] + xx[26] * xx[69];
//  xx[16] = xx[66] * xx[111] - (xx[0] * xx[30] + xx[3] * xx[29]);
//  xx[44] = xx[18] + xx[96] + xx[84] + xx[16] * xx[69];
//  xx[18] = xx[8] * xx[44] - xx[7] * xx[43];
//  xx[45] = xx[62] + xx[65] * xx[69];
//  xx[62] = xx[7] * xx[45];
//  xx[71] = xx[8] * xx[45];
//  xx[72] = xx[45] - xx[2] * (xx[18] * xx[13] + xx[7] * xx[62] + xx[8] * xx[71]);
//  xx[73] = xx[43] + xx[2] * (xx[7] * xx[18] - xx[62] * xx[13]);
//  xx[74] = xx[72];
//  xx[75] = xx[44] + xx[2] * (xx[71] * xx[13] - xx[8] * xx[18]);
//  pm_math_Vector3_cross_ra<adouble>(xx + 36, xx + 73, xx + 43);
//  xx[18] = xx[66] * state[1];
//  xx[62] = xx[1] * xx[18];
//  xx[1] = xx[13] * xx[13];
//  xx[71] = xx[2] * xx[7] * xx[13];
//  xx[73] = xx[2] * xx[8] * xx[7];
//  xx[74] = xx[2] * xx[8] * xx[13];
//  xx[82] = xx[2] * (xx[1] + xx[8] * xx[8]) - xx[46];
//  xx[83] = - xx[71];
//  xx[84] = xx[73];
//  xx[85] = xx[71];
//  xx[86] = xx[2] * xx[1] - xx[46];
//  xx[87] = - xx[74];
//  xx[88] = xx[73];
//  xx[89] = xx[74];
//  xx[90] = xx[2] * (xx[1] + xx[7] * xx[7]) - xx[46];
//  xx[1] = xx[26] / xx[68];
//  xx[46] = xx[65] / xx[68];
//  xx[71] = xx[16] / xx[68];
//  xx[91] = - (xx[5] + xx[12] * xx[1]);
//  xx[92] = - (xx[11] + xx[12] * xx[46]);
//  xx[93] = - (xx[30] + xx[12] * xx[71]);
//  xx[94] = - (xx[28] + xx[25] * xx[1]);
//  xx[95] = - (xx[57] + xx[25] * xx[46]);
//  xx[96] = - (xx[58] + xx[25] * xx[71]);
//  xx[97] = - (xx[59] + xx[17] * xx[1]);
//  xx[98] = - (xx[60] + xx[17] * xx[46]);
//  xx[99] = - (xx[29] + xx[17] * xx[71]);
//  pm_math_Matrix3x3_composeTranspose_ra<adouble>(xx + 91, xx + 82, xx + 113);
//  pm_math_Matrix3x3_compose_ra<adouble>(xx + 82, xx + 113, xx + 91);
//  xx[5] = xx[65] * xx[1];
//  xx[11] = xx[16] * xx[1];
//  xx[28] = xx[16] * xx[46];
//  xx[113] = xx[104] - xx[26] * xx[1] + xx[40];
//  xx[114] = xx[105] - xx[5];
//  xx[115] = xx[106] - xx[11];
//  xx[116] = xx[107] - xx[5];
//  xx[117] = xx[61] - xx[65] * xx[46];
//  xx[118] = xx[109] - xx[28];
//  xx[119] = xx[110] - xx[11];
//  xx[120] = xx[111] - xx[28];
//  xx[121] = xx[112] - xx[16] * xx[71] + xx[40];
//  pm_math_Matrix3x3_composeTranspose_ra<adouble>(xx + 113, xx + 82, xx + 100);
//  pm_math_Matrix3x3_compose_ra<adouble>(xx + 82, xx + 100, xx + 109);
//  pm_math_Matrix3x3_postCross_ra(xx + 109, xx + 36, xx + 100);
//  xx[5] = xx[93] - xx[106];
//  xx[11] = xx[18] * xx[6];
//  xx[6] = xx[91] - xx[100];
//  xx[16] = xx[99] - xx[108];
//  xx[18] = xx[97] - xx[102];
//  xx[26] = xx[12] / xx[68];
//  xx[28] = xx[25] * xx[26];
//  xx[29] = xx[17] * xx[26];
//  xx[30] = xx[25] / xx[68];
//  xx[44] = xx[17] * xx[30];
//  xx[57] = xx[17] / xx[68];
//  xx[117] = xx[23] - xx[12] * xx[26];
//  xx[118] = xx[31] - xx[28];
//  xx[119] = xx[47] - xx[29];
//  xx[120] = xx[51] - xx[28];
//  xx[121] = xx[52] - xx[25] * xx[30] + xx[20];
//  xx[122] = xx[53] - xx[44];
//  xx[123] = xx[54] - xx[29];
//  xx[124] = xx[55] - xx[44];
//  xx[125] = xx[56] - xx[17] * xx[57];
//  pm_math_Matrix3x3_composeTranspose_ra<adouble>(xx + 117, xx + 82, xx + 126);
//  pm_math_Matrix3x3_compose_ra<adouble>(xx + 82, xx + 126, xx + 117);
//  pm_math_Matrix3x3_postCross_ra<adouble>(xx + 91, xx + 36, xx + 82);
//  pm_math_Matrix3x3_preCross_ra<adouble>(xx + 100, xx + 36, xx + 126);
//  xx[12] = xx[92] - xx[103];
//  xx[17] = xx[98] - xx[105];
//  xx[20] = xx[0] * xx[12] + xx[3] * xx[17] + xx[66] * (xx[40] + xx[113]);
//  xx[23] = xx[0] * (xx[0] * (xx[19] + xx[117] - xx[82] - xx[82] - xx[126]) + xx
//                                                                             [3] * (xx[119] - xx[84] - xx[88] - xx[128]) + xx[66] * xx[12])
//           + xx[3] * (xx[0] * (xx[123] - xx[88] - xx[84] - xx[132]) + xx[3] * (xx[21] +
//                                                                               xx[125] - xx[90] - xx[90] - xx[134]) + xx[66] * xx[17]) + xx[66] * xx[20];
////    ii[0] = factorSymmetricPosDef(xx + 23, 1, xx + 12);
////    if (ii[0] != 0) {
////      return sm_ssci_recordRunTimeError(
////              "sm:compiler:messages:simulationErrors:DegenerateMass",
////              "'manipulator_ip/Revolute BaseArm1' has a degenerate mass distribution on its follower side.",
////              neDiagMgr);
////    }
//
//  xx[36] = (xx[0] * xx[6] + xx[3] * xx[18] + xx[66] * xx[110]) / xx[23];
//  xx[37] = xx[20] / xx[23];
//  xx[38] = (xx[0] * xx[5] + xx[3] * xx[16] + xx[66] * xx[116]) / xx[23];
//  xx[12] = 9.81;
//  xx[17] = 0.5323909858203203;
//  xx[19] = xx[4] * state[0];
//  xx[4] = sin(xx[19]);
//  xx[20] = xx[0] * xx[4];
//  xx[21] = xx[17] * xx[20];
//  xx[25] = 0.4653599018150013;
//  xx[28] = cos(xx[19]);
//  xx[19] = xx[25] * xx[28];
//  xx[29] = xx[3] * xx[4];
//  xx[4] = xx[25] * xx[29];
//  xx[31] = xx[21] - xx[19] + xx[4];
//  xx[40] = xx[12] * xx[31];
//  xx[44] = xx[17] * xx[28];
//  xx[28] = xx[17] * xx[29];
//  xx[17] = xx[25] * xx[20];
//  xx[20] = xx[44] + xx[28] - xx[17];
//  xx[25] = xx[21] + xx[19] + xx[4];
//  xx[4] = xx[12] * xx[25];
//  xx[19] = xx[44] - xx[28] + xx[17];
//  xx[17] = xx[2] * (xx[40] * xx[20] - xx[4] * xx[19]);
//  xx[21] = xx[2] * (xx[4] * xx[20] + xx[40] * xx[19]);
//  xx[19] = xx[12] - xx[2] * (xx[4] * xx[25] + xx[40] * xx[31]);
//  xx[51] = xx[17];
//  xx[52] = xx[21];
//  xx[53] = xx[19];
//  xx[4] = (input[0] - (xx[0] * (xx[70] + xx[2] * (xx[7] * xx[14] - xx[7] * xx[27]
//                                                                   * xx[13]) + xx[43] + xx[62] * xx[5] - xx[11] * xx[6]) + xx[3] *
//                                                                                                                           (xx[10] + xx[2] * (xx[8] * xx[27] * xx[13] - xx[8] * xx
//                                                                                                                           [14]) + xx[45] + xx[62] * xx[16] - xx[11] * xx[18]) + xx[66] *
//                                                                                                                                                                                 (xx[72] + xx[62] * xx[114] - xx[11] * xx[112]))) / xx[23]
//          - pm_math_Vector3_dot_ra<adouble>(xx + 36, xx + 51);
//  xx[27] = xx[26];
//  xx[28] = xx[30];
//  xx[29] = xx[57];
//  xx[5] = xx[0] * xx[4];
//  xx[6] = xx[3] * xx[4];
//  xx[10] = xx[8] * xx[6] - xx[5] * xx[7];
//  xx[12] = xx[5] + xx[2] * xx[7] * xx[10];
//  xx[14] = xx[2] * xx[10] * xx[13];
//  xx[16] = xx[6] - xx[2] * xx[8] * xx[10];
//  xx[36] = xx[12];
//  xx[37] = xx[14];
//  xx[38] = xx[16];
//  xx[43] = xx[1];
//  xx[44] = xx[46];
//  xx[45] = xx[71];
//  xx[1] = xx[17] - xx[11] - xx[6] * xx[64];
//  xx[10] = xx[21] + xx[66] * xx[4] + xx[6] * xx[63] - xx[5] * xx[39];
//  xx[6] = xx[7] * xx[10];
//  xx[11] = xx[19] + xx[62] + xx[5] * xx[64];
//  xx[5] = xx[8] * xx[11] - xx[7] * xx[1];
//  xx[17] = xx[1] + xx[2] * (xx[6] * xx[13] + xx[7] * xx[5]);
//  xx[1] = xx[8] * xx[10];
//  xx[18] = xx[10] + xx[2] * (xx[5] * xx[13] - (xx[7] * xx[6] + xx[8] * xx[1]));
//  xx[6] = xx[11] - xx[2] * (xx[1] * xx[13] + xx[8] * xx[5]);
//  xx[19] = xx[17];
//  xx[20] = xx[18];
//  xx[21] = xx[6];
//  xx[1] = xx[69] - (pm_math_Vector3_dot_ra<adouble>(xx + 27, xx + 36) +
//                    pm_math_Vector3_dot_ra<adouble>(xx + 43, xx + 19));
//  xx[19] = xx[12] + xx[0] * xx[1] + xx[67];
//  xx[20] = xx[14] + xx[81];
//  xx[21] = xx[16] + xx[3] * xx[1] - xx[9];
//  pm_math_Quaternion_inverseXform_ra<adouble>(xx + 32, xx + 19, xx + 7);
//  pm_math_Vector3_cross_ra<adouble>(xx + 19, xx + 48, xx + 10);
//  xx[19] = xx[17] + xx[78] + xx[10];
//  xx[20] = xx[18] + xx[66] * xx[1] + xx[15] + xx[11];
//  xx[21] = xx[6] + xx[42] + xx[12];
//  pm_math_Quaternion_inverseXform_ra<adouble>(xx + 32, xx + 19, xx + 5);
//  deriv[0] = state[1];
//  deriv[1] = xx[4];
//  deriv[2] = state[3];
//  deriv[3] = xx[1];
//  deriv[4] = state[5];
//  deriv[5] = xx[41] - (xx[24] * xx[9] + xx[22] * xx[6]);
}