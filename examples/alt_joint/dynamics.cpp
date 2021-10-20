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

adouble ad_fmod(adouble x, adouble y)
{
  struct u {adouble f; uint64_t i;} ux, uy;
  ux.f = x;
  uy.f = y;
  int ex = ux.i>>52 & 0x7ff;
  int ey = uy.i>>52 & 0x7ff;
  int sx = ux.i>>63;
  uint64_t i;

  /* in the followings uxi should be ux.i, but then gcc wrongly adds */
  /* float load/store to inner loops ruining performance and code size */
  uint64_t uxi = ux.i;

  if (uy.i<<1 == 0 || isnan(y.value()) || ex == 0x7ff)
    return (x*y)/(x*y);
  if (uxi<<1 <= uy.i<<1) {
    if (uxi<<1 == uy.i<<1)
      return 0*x;
    return x;
  }

  /* normalize x and y */
  if (!ex) {
    for (i = uxi<<12; i>>63 == 0; ex--, i <<= 1);
    uxi <<= -ex + 1;
  } else {
    uxi &= -1ULL >> 12;
    uxi |= 1ULL << 52;
  }
  if (!ey) {
    for (i = uy.i<<12; i>>63 == 0; ey--, i <<= 1);
    uy.i <<= -ey + 1;
  } else {
    uy.i &= -1ULL >> 12;
    uy.i |= 1ULL << 52;
  }

  /* x mod y */
  for (; ex > ey; ex--) {
    i = uxi - uy.i;
    if (i >> 63 == 0) {
      if (i == 0)
        return 0*x;
      uxi = i;
    }
    uxi <<= 1;
  }
  i = uxi - uy.i;
  if (i >> 63 == 0) {
    if (i == 0)
      return 0*x;
    uxi = i;
  }
  for (; uxi>>52 == 0; uxi <<= 1, ex--);

  /* scale result */
  if (ex > 0) {
    uxi -= 1ULL << 52;
    uxi |= (uint64_t)ex << 52;
  } else {
    uxi >>= -ex + 1;
  }
  uxi |= (uint64_t)sx << 63;
  ux.i = uxi;
  return ux.f;
}

/*!
 * \brief normalize_angle_positive
 *
 *        Normalizes the angle to be 0 to 2*M_PI
 *        It takes and returns radians.
 */
static inline adouble normalize_angle_positive(adouble angle)
{
  return ad_fmod(ad_fmod(angle, 2.0*M_PI) + 2.0*M_PI, 2.0*M_PI);
}


void dae(adouble* deriv, adouble* path, adouble* state,
         adouble* input, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{

  bool auto_deriv = (workspace->algorithm->derivatives == "automatic")? true: false;

  if (workspace->problem->mj_backend)
  {
    assert(auto_deriv == false);
    workspace->problem->mj_handle.forwardSimulate(state, input, deriv);
//    std::cout << "deriv: " <<  deriv[0] << "\t" << deriv[1] << "\t" << deriv[2] << "\t" << deriv[3] << std::endl;
    return;
  }

  adouble xx[165];
  xx[0] = 0.7071067811865476;
  xx[1] = 0.5;
  xx[2] = xx[1] * state[2];
  xx[3] = xx[0] * cos(xx[2]);
  xx[4] = xx[0] * sin(xx[2]);
  xx[2] = - xx[4];
  xx[5] = xx[3];
  xx[6] = xx[3];
  xx[7] = xx[2];
  xx[8] = xx[4];
  xx[9] = 2.0;
  xx[10] = xx[4] * state[1];
  xx[11] = xx[3] * state[1];
  xx[12] = xx[9] * (xx[3] * xx[10] + xx[4] * xx[11]);
  xx[13] = xx[3] * xx[11];
  xx[11] = xx[4] * xx[10];
  xx[10] = xx[9] * (xx[13] - xx[11]);
  xx[14] = state[1] - xx[9] * (xx[13] + xx[11]);
  xx[11] = xx[14] + state[3];
  xx[15] = xx[12];
  xx[16] = xx[10];
  xx[17] = xx[11];
  xx[13] = 4.166666666666641e-4;
  xx[18] = 0.08354166666666671;
  xx[19] = xx[13] * xx[12];
  xx[20] = xx[18] * xx[10];
  xx[21] = xx[18] * xx[11];
  pm_math_Vector3_cross_ra(xx + 15, xx + 19, xx + 22);
  xx[19] = xx[1] * state[4];
  xx[20] = xx[0] * cos(xx[19]);
  xx[21] = xx[0] * sin(xx[19]);
  xx[25] = xx[20];
  xx[26] = - xx[20];
  xx[27] = xx[21];
  xx[28] = xx[21];
  pm_math_Quaternion_inverseXform_ra(xx + 25, xx + 15, xx + 29);
  xx[19] = xx[31] + state[5];
  xx[32] = xx[29];
  xx[33] = xx[30];
  xx[34] = xx[19];
  xx[35] = xx[13] * xx[29];
  xx[36] = xx[18] * xx[30];
  xx[37] = xx[18] * xx[19];
  pm_math_Vector3_cross_ra(xx + 32, xx + 35, xx + 38);
  xx[35] = xx[1] * state[6];
  xx[36] = xx[0] * cos(xx[35]);
  xx[37] = xx[0] * sin(xx[35]);
  xx[41] = xx[36];
  xx[42] = xx[36];
  xx[43] = - xx[37];
  xx[44] = xx[37];
  pm_math_Quaternion_inverseXform_ra(xx + 41, xx + 32, xx + 45);
  xx[35] = xx[47] + state[7];
  xx[48] = xx[45];
  xx[49] = xx[46];
  xx[50] = xx[35];
  xx[51] = xx[13] * xx[45];
  xx[52] = xx[18] * xx[46];
  xx[53] = xx[18] * xx[35];
  pm_math_Vector3_cross_ra(xx + 48, xx + 51, xx + 54);
  xx[51] = xx[1] * state[8];
  xx[52] = xx[0] * cos(xx[51]);
  xx[53] = xx[0] * sin(xx[51]);
  xx[57] = xx[52];
  xx[58] = - xx[52];
  xx[59] = xx[53];
  xx[60] = xx[53];
  pm_math_Quaternion_inverseXform_ra(xx + 57, xx + 48, xx + 61);
  xx[51] = xx[63] + state[9];
  xx[64] = xx[61];
  xx[65] = xx[62];
  xx[66] = xx[51];
  xx[67] = xx[13] * xx[61];
  xx[68] = xx[18] * xx[62];
  xx[69] = xx[18] * xx[51];
  pm_math_Vector3_cross_ra(xx + 64, xx + 67, xx + 70);
  xx[64] = xx[1] * xx[53];
  xx[65] = xx[64] * xx[53];
  xx[66] = xx[64] * xx[52];
  xx[67] = xx[1] - (xx[9] * (xx[65] + xx[65]) - xx[1]);
  xx[68] = - (xx[9] * (xx[66] - xx[66]));
  xx[69] = - (xx[9] * (xx[66] + xx[66]));
  pm_math_Vector3_cross_ra(xx + 48, xx + 67, xx + 64);
  pm_math_Vector3_cross_ra(xx + 48, xx + 64, xx + 73);
  pm_math_Quaternion_inverseXform_ra(xx + 57, xx + 73, xx + 48);
  xx[64] = 0.3335416666666667;

  xx[65] = (input[4] - (xx[72] + xx[1] * xx[49])) / xx[64];
  xx[73] = xx[70] + xx[13] * state[9] * xx[62];
  xx[74] = xx[71] - xx[18] * state[9] * xx[61];
  xx[75] = xx[72] + xx[18] * xx[65];
  pm_math_Quaternion_xform_ra(xx + 57, xx + 73, xx + 70);
  xx[62] = xx[1] * state[9];
  xx[73] = xx[48] - xx[62] * (xx[63] + xx[51]);
  xx[74] = xx[49] + xx[1] * xx[65];
  xx[75] = xx[62] * (xx[61] + xx[61]) + xx[50];
  pm_math_Quaternion_xform_ra(xx + 57, xx + 73, xx + 48);
  pm_math_Vector3_cross_ra(xx + 67, xx + 48, xx + 61);
  xx[51] = state[7] * xx[46];
  xx[46] = xx[52] * xx[52];
  xx[66] = 1.0;
  xx[73] = xx[9] * (xx[46] + xx[46]) - xx[66];
  xx[74] = xx[52] * xx[53];
  xx[52] = xx[9] * (xx[74] + xx[74]);
  xx[75] = xx[9] * (xx[74] - xx[74]);
  xx[76] = xx[74] - xx[74];
  xx[77] = xx[9] * xx[76];
  xx[78] = xx[53] * xx[53];
  xx[53] = xx[9] * (xx[46] + xx[78]) - xx[66];
  xx[79] = xx[9] * (xx[78] + xx[46]);
  xx[80] = xx[74] + xx[74];
  xx[74] = - (xx[9] * xx[80]);
  xx[81] = xx[9] * (xx[78] - xx[46]);
  xx[82] = xx[73];
  xx[83] = - xx[52];
  xx[84] = xx[75];
  xx[85] = xx[77];
  xx[86] = xx[53];
  xx[87] = xx[79];
  xx[88] = xx[74];
  xx[89] = xx[81];
  xx[90] = xx[53];
  xx[46] = 8.333333333333282e-4;
  xx[78] = xx[18] / xx[64];
  xx[91] = xx[18] - xx[18] * xx[78];
  xx[92] = xx[13] * xx[73];
  xx[93] = xx[46] * xx[76];
  xx[94] = - (xx[46] * xx[80]);
  xx[95] = - (xx[18] * xx[52]);
  xx[96] = xx[18] * xx[53];
  xx[97] = xx[18] * xx[81];
  xx[98] = xx[75] * xx[91];
  xx[99] = xx[79] * xx[91];
  xx[100] = xx[91] * xx[53];
  pm_math_Matrix3x3_compose_ra(xx + 82, xx + 92, xx + 101);
  xx[46] = xx[1] / xx[64];
  xx[64] = xx[18] * xx[46];
  xx[76] = xx[64] * xx[52];
  xx[80] = xx[75] * xx[76];
  xx[91] = xx[64] * xx[53];
  xx[92] = xx[75] * xx[91];
  xx[93] = xx[64] * xx[81];
  xx[64] = xx[75] * xx[93];
  xx[94] = xx[79] * xx[76];
  xx[95] = xx[79] * xx[91];
  xx[96] = xx[79] * xx[93];
  xx[97] = xx[76] * xx[53];
  xx[76] = xx[91] * xx[53];
  xx[91] = xx[93] * xx[53];
  xx[110] = xx[80];
  xx[111] = - xx[92];
  xx[112] = - xx[64];
  xx[113] = xx[94];
  xx[114] = - xx[95];
  xx[115] = - xx[96];
  xx[116] = xx[97];
  xx[117] = - xx[76];
  xx[118] = - xx[91];
  pm_math_Matrix3x3_postCross_ra(xx + 110, xx + 67, xx + 119);
  xx[93] = xx[66] - xx[1] * xx[46];
  xx[110] = xx[73];
  xx[111] = xx[77];
  xx[112] = xx[74];
  xx[113] = - (xx[52] * xx[93]);
  xx[114] = xx[93] * xx[53];
  xx[115] = xx[81] * xx[93];
  xx[116] = xx[75];
  xx[117] = xx[79];
  xx[118] = xx[53];
  pm_math_Matrix3x3_compose_ra(xx + 82, xx + 110, xx + 128);
  pm_math_Matrix3x3_postCross_ra(xx + 128, xx + 67, xx + 81);
  pm_math_Matrix3x3_preCross_ra(xx + 81, xx + 67, xx + 110);
  xx[52] = xx[13] + xx[101] - xx[119] - xx[119] - xx[110];
  xx[53] = state[7] * xx[45];
  xx[73] = xx[102] - xx[120] - xx[122] - xx[111];
  xx[74] = xx[80] - xx[81];
  xx[75] = xx[92] + xx[84];
  xx[77] = xx[64] + xx[87];
  xx[64] = xx[82] - xx[94];
  xx[79] = xx[95] + xx[85];
  xx[80] = xx[96] + xx[88];
  xx[81] = xx[97] - xx[83];
  xx[82] = xx[76] + xx[86];
  xx[76] = xx[91] + xx[89];
  xx[83] = xx[74];
  xx[84] = - xx[75];
  xx[85] = - xx[77];
  xx[86] = - xx[64];
  xx[87] = - xx[79];
  xx[88] = - xx[80];
  xx[89] = xx[81];
  xx[90] = - xx[82];
  xx[91] = - xx[76];
  xx[92] = xx[1] * xx[37];
  xx[93] = xx[92] * xx[37];
  xx[94] = xx[92] * xx[36];
  xx[95] = xx[1] - (xx[9] * (xx[93] + xx[93]) - xx[1]);
  xx[96] = - (xx[9] * (xx[94] - xx[94]));
  xx[97] = xx[9] * (xx[94] + xx[94]);
  pm_math_Vector3_cross_ra(xx + 32, xx + 95, xx + 92);
  pm_math_Vector3_cross_ra(xx + 32, xx + 92, xx + 98);
  pm_math_Quaternion_inverseXform_ra(xx + 41, xx + 98, xx + 32);
  xx[92] = xx[1] * state[7];
  xx[93] = xx[32] - xx[92] * (xx[47] + xx[35]);
  xx[32] = xx[92] * (xx[45] + xx[45]) + xx[34];
  xx[98] = xx[93];
  xx[99] = xx[33];
  xx[100] = xx[32];
  pm_math_Matrix3x3_xform_ra(xx + 83, xx + 98, xx + 137);
  xx[34] = xx[103] - xx[121] - xx[125] - xx[112];
  xx[35] = xx[34] - xx[1] * xx[75];
  xx[45] = xx[107] - xx[125] - xx[121] - xx[116];
  xx[47] = xx[108] - xx[126] - xx[124] - xx[117];
  xx[83] = xx[56] + xx[72] + xx[63] + xx[51] * xx[45] - xx[53] * xx[47] + xx[139];
  xx[84] = xx[66] + xx[132];
  xx[140] = xx[66] + xx[128];
  xx[141] = xx[129];
  xx[142] = xx[130];
  xx[143] = xx[131];
  xx[144] = xx[84];
  xx[145] = xx[133];
  xx[146] = xx[134];
  xx[147] = xx[135];
  xx[148] = xx[66] + xx[136];
  pm_math_Matrix3x3_xform_ra(xx + 140, xx + 98, xx + 85);
  xx[88] = xx[49] + xx[86] - (xx[51] * xx[75] - xx[53] * xx[79]);
  xx[49] = xx[18] + xx[109] - xx[127] - xx[127] - xx[118];
  xx[89] = xx[49] - xx[1] * xx[82];
  xx[90] = xx[1] * xx[84] - xx[82];
  xx[91] = xx[89] + xx[1] * xx[90];

  xx[92] = (input[3] - (xx[83] + xx[1] * xx[88])) / xx[91];
  xx[56] = xx[104] - xx[122] - xx[120] - xx[113];
  xx[63] = xx[18] + xx[105] - xx[123] - xx[123] - xx[114];
  xx[72] = xx[106] - xx[124] - xx[126] - xx[115];
  xx[94] = xx[72] - xx[1] * xx[79];
  xx[98] = xx[54] + xx[70] + xx[61] + xx[51] * xx[52] - xx[53] * xx[73] + xx[137]
           + xx[35] * xx[92];
  xx[99] = xx[55] + xx[71] + xx[62] + xx[51] * xx[56] - xx[53] * xx[63] + xx[138]
           + xx[94] * xx[92];
  xx[100] = xx[83] + xx[89] * xx[92];
  pm_math_Quaternion_xform_ra(xx + 41, xx + 98, xx + 101);
  xx[54] = xx[81] + xx[1] * xx[129];
  xx[55] = xx[1] * xx[135] - xx[76];
  xx[98] = xx[48] + xx[51] * xx[74] + xx[53] * xx[64] + xx[85] + xx[54] * xx[92];
  xx[99] = xx[88] + xx[90] * xx[92];
  xx[100] = xx[50] + xx[87] - (xx[51] * xx[77] - xx[53] * xx[80]) + xx[55] * xx
  [92];
  pm_math_Quaternion_xform_ra(xx + 41, xx + 98, xx + 85);
  pm_math_Vector3_cross_ra(xx + 95, xx + 85, xx + 98);
  xx[48] = state[5] * xx[30];
  xx[30] = xx[36] * xx[36];
  xx[50] = xx[36] * xx[37];
  xx[36] = xx[9] * (xx[50] + xx[50]);
  xx[61] = xx[9] * (xx[50] - xx[50]);
  xx[50] = xx[37] * xx[37];
  xx[37] = xx[9] * (xx[30] + xx[50]) - xx[66];
  xx[104] = xx[9] * (xx[30] + xx[30]) - xx[66];
  xx[105] = - xx[36];
  xx[106] = xx[61];
  xx[107] = xx[61];
  xx[108] = xx[37];
  xx[109] = - (xx[9] * (xx[50] + xx[30]));
  xx[110] = xx[36];
  xx[111] = xx[9] * (xx[30] - xx[50]);
  xx[112] = xx[37];
  xx[30] = xx[35] / xx[91];
  xx[36] = xx[94] * xx[30];
  xx[37] = xx[89] * xx[30];
  xx[50] = xx[94] / xx[91];
  xx[61] = xx[89] * xx[50];
  xx[62] = xx[89] / xx[91];
  xx[113] = xx[52] - xx[35] * xx[30];
  xx[114] = xx[73] - xx[36];
  xx[115] = xx[34] - xx[37];
  xx[116] = xx[56] - xx[36];
  xx[117] = xx[63] - xx[94] * xx[50];
  xx[118] = xx[72] - xx[61];
  xx[119] = xx[45] - xx[37];
  xx[120] = xx[47] - xx[61];
  xx[121] = xx[49] - xx[89] * xx[62];
  pm_math_Matrix3x3_composeTranspose_ra(xx + 113, xx + 104, xx + 137);
  pm_math_Matrix3x3_compose_ra(xx + 104, xx + 137, xx + 113);
  xx[34] = xx[54] / xx[91];
  xx[36] = xx[90] / xx[91];
  xx[37] = xx[55] / xx[91];
  xx[137] = xx[74] - xx[35] * xx[34];
  xx[138] = - (xx[75] + xx[35] * xx[36]);
  xx[139] = - (xx[77] + xx[35] * xx[37]);
  xx[140] = - (xx[64] + xx[94] * xx[34]);
  xx[141] = - (xx[79] + xx[94] * xx[36]);
  xx[142] = - (xx[80] + xx[94] * xx[37]);
  xx[143] = xx[81] - xx[89] * xx[34];
  xx[144] = - (xx[82] + xx[89] * xx[36]);
  xx[145] = - (xx[76] + xx[89] * xx[37]);
  pm_math_Matrix3x3_composeTranspose_ra(xx + 137, xx + 104, xx + 146);
  pm_math_Matrix3x3_compose_ra(xx + 104, xx + 146, xx + 137);
  pm_math_Matrix3x3_postCross_ra(xx + 137, xx + 95, xx + 146);
  xx[35] = xx[90] * xx[34];
  xx[45] = xx[55] * xx[34];
  xx[47] = xx[55] * xx[36];
  xx[155] = xx[128] - xx[54] * xx[34] + xx[66];
  xx[156] = xx[129] - xx[35];
  xx[157] = xx[130] - xx[45];
  xx[158] = xx[131] - xx[35];
  xx[159] = xx[84] - xx[90] * xx[36];
  xx[160] = xx[133] - xx[47];
  xx[161] = xx[134] - xx[45];
  xx[162] = xx[135] - xx[47];
  xx[163] = xx[136] - xx[55] * xx[37] + xx[66];
  pm_math_Matrix3x3_composeTranspose_ra(xx + 155, xx + 104, xx + 122);
  pm_math_Matrix3x3_compose_ra(xx + 104, xx + 122, xx + 155);
  pm_math_Matrix3x3_postCross_ra(xx + 155, xx + 95, xx + 104);
  pm_math_Matrix3x3_preCross_ra(xx + 104, xx + 95, xx + 122);
  xx[35] = xx[13] + xx[113] - xx[146] - xx[146] - xx[122];
  xx[45] = state[5] * xx[29];
  xx[47] = xx[114] - xx[147] - xx[149] - xx[123];
  xx[49] = xx[137] - xx[104];
  xx[52] = xx[138] - xx[107];
  xx[54] = xx[139] - xx[110];
  xx[55] = xx[140] - xx[105];
  xx[56] = xx[141] - xx[108];
  xx[61] = xx[142] - xx[111];
  xx[63] = xx[143] - xx[106];
  xx[64] = xx[144] - xx[109];
  xx[70] = xx[145] - xx[112];
  xx[104] = xx[49];
  xx[105] = xx[52];
  xx[106] = xx[54];
  xx[107] = xx[55];
  xx[108] = xx[56];
  xx[109] = xx[61];
  xx[110] = xx[63];
  xx[111] = xx[64];
  xx[112] = xx[70];
  xx[71] = xx[1] * xx[21];
  xx[72] = xx[71] * xx[21];
  xx[73] = xx[71] * xx[20];
  xx[74] = xx[1] - (xx[9] * (xx[72] + xx[72]) - xx[1]);
  xx[75] = - (xx[9] * (xx[73] - xx[73]));
  xx[76] = - (xx[9] * (xx[73] + xx[73]));
  pm_math_Vector3_cross_ra(xx + 15, xx + 74, xx + 71);
  pm_math_Vector3_cross_ra(xx + 15, xx + 71, xx + 79);
  pm_math_Quaternion_inverseXform_ra(xx + 25, xx + 79, xx + 15);
  xx[71] = xx[1] * state[5];
  xx[72] = xx[15] - xx[71] * (xx[31] + xx[19]);
  xx[15] = xx[71] * (xx[29] + xx[29]) + xx[17];
  xx[79] = xx[72];
  xx[80] = xx[16];
  xx[81] = xx[15];
  pm_math_Matrix3x3_xform_ra(xx + 104, xx + 79, xx + 82);
  xx[17] = xx[115] - xx[148] - xx[152] - xx[124];
  xx[19] = xx[17] + xx[1] * xx[52];
  xx[29] = xx[119] - xx[152] - xx[148] - xx[128];
  xx[31] = xx[120] - xx[153] - xx[151] - xx[129];
  xx[71] = xx[40] + xx[103] + xx[100] + xx[48] * xx[29] - xx[45] * xx[31] + xx
  [84];
  xx[73] = xx[66] + xx[159];
  xx[104] = xx[66] + xx[155];
  xx[105] = xx[156];
  xx[106] = xx[157];
  xx[107] = xx[158];
  xx[108] = xx[73];
  xx[109] = xx[160];
  xx[110] = xx[161];
  xx[111] = xx[162];
  xx[112] = xx[66] + xx[163];
  pm_math_Matrix3x3_xform_ra(xx + 104, xx + 79, xx + 88);
  xx[77] = xx[86] + xx[48] * xx[52] - xx[45] * xx[56] + xx[89];
  xx[79] = xx[18] + xx[121] - xx[154] - xx[154] - xx[130];
  xx[80] = xx[79] + xx[1] * xx[64];
  xx[81] = xx[64] + xx[1] * xx[73];
  xx[86] = xx[80] + xx[1] * xx[81];

  xx[91] = (input[2] - (xx[71] + xx[1] * xx[77])) / xx[86];
  xx[40] = xx[116] - xx[149] - xx[147] - xx[125];
  xx[94] = xx[18] + xx[117] - xx[150] - xx[150] - xx[126];
  xx[84] = xx[118] - xx[151] - xx[153] - xx[127];
  xx[100] = xx[84] + xx[1] * xx[56];
  xx[103] = xx[38] + xx[101] + xx[98] + xx[48] * xx[35] - xx[45] * xx[47] + xx
  [82] + xx[19] * xx[91];
  xx[104] = xx[39] + xx[102] + xx[99] + xx[48] * xx[40] - xx[45] * xx[94] + xx
  [83] + xx[100] * xx[91];
  xx[105] = xx[71] + xx[80] * xx[91];
  pm_math_Quaternion_xform_ra(xx + 25, xx + 103, xx + 106);
  xx[38] = xx[63] + xx[1] * xx[156];
  xx[39] = xx[70] + xx[1] * xx[162];
  xx[101] = xx[85] + xx[48] * xx[49] - xx[45] * xx[55] + xx[88] + xx[38] * xx[91];
  xx[102] = xx[77] + xx[81] * xx[91];
  xx[103] = xx[87] + xx[48] * xx[54] - xx[45] * xx[61] + xx[90] + xx[39] * xx[91];
  pm_math_Quaternion_xform_ra(xx + 25, xx + 101, xx + 87);
  pm_math_Vector3_cross_ra(xx + 74, xx + 87, xx + 101);
  xx[71] = xx[10] * state[3];
  xx[10] = xx[20] * xx[20];
  xx[77] = xx[20] * xx[21];
  xx[20] = xx[21] * xx[21];
  xx[21] = xx[9] * (xx[10] + xx[20]) - xx[66];
  xx[109] = xx[9] * (xx[10] + xx[10]) - xx[66];
  xx[110] = - (xx[9] * (xx[77] + xx[77]));
  xx[111] = xx[9] * (xx[77] - xx[77]);
  xx[112] = xx[9] * (xx[77] - xx[77]);
  xx[113] = xx[21];
  xx[114] = xx[9] * (xx[20] + xx[10]);
  xx[115] = - (xx[9] * (xx[77] + xx[77]));
  xx[116] = xx[9] * (xx[20] - xx[10]);
  xx[117] = xx[21];
  xx[10] = xx[19] / xx[86];
  xx[20] = xx[100] * xx[10];
  xx[21] = xx[80] * xx[10];
  xx[77] = xx[100] / xx[86];
  xx[82] = xx[80] * xx[77];
  xx[83] = xx[80] / xx[86];
  xx[118] = xx[35] - xx[19] * xx[10];
  xx[119] = xx[47] - xx[20];
  xx[120] = xx[17] - xx[21];
  xx[121] = xx[40] - xx[20];
  xx[122] = xx[94] - xx[100] * xx[77];
  xx[123] = xx[84] - xx[82];
  xx[124] = xx[29] - xx[21];
  xx[125] = xx[31] - xx[82];
  xx[126] = xx[79] - xx[80] * xx[83];
  pm_math_Matrix3x3_composeTranspose_ra(xx + 118, xx + 109, xx + 127);
  pm_math_Matrix3x3_compose_ra(xx + 109, xx + 127, xx + 118);
  xx[17] = xx[38] / xx[86];
  xx[20] = xx[81] / xx[86];
  xx[21] = xx[39] / xx[86];
  xx[127] = xx[49] - xx[19] * xx[17];
  xx[128] = xx[52] - xx[19] * xx[20];
  xx[129] = xx[54] - xx[19] * xx[21];
  xx[130] = xx[55] - xx[100] * xx[17];
  xx[131] = xx[56] - xx[100] * xx[20];
  xx[132] = xx[61] - xx[100] * xx[21];
  xx[133] = xx[63] - xx[80] * xx[17];
  xx[134] = xx[64] - xx[80] * xx[20];
  xx[135] = xx[70] - xx[80] * xx[21];
  pm_math_Matrix3x3_composeTranspose_ra(xx + 127, xx + 109, xx + 136);
  pm_math_Matrix3x3_compose_ra(xx + 109, xx + 136, xx + 127);
  pm_math_Matrix3x3_postCross_ra(xx + 127, xx + 74, xx + 136);
  xx[19] = xx[81] * xx[17];
  xx[29] = xx[39] * xx[17];
  xx[31] = xx[39] * xx[20];
  xx[145] = xx[155] - xx[38] * xx[17] + xx[66];
  xx[146] = xx[156] - xx[19];
  xx[147] = xx[157] - xx[29];
  xx[148] = xx[158] - xx[19];
  xx[149] = xx[73] - xx[81] * xx[20];
  xx[150] = xx[160] - xx[31];
  xx[151] = xx[161] - xx[29];
  xx[152] = xx[162] - xx[31];
  xx[153] = xx[163] - xx[39] * xx[21] + xx[66];
  pm_math_Matrix3x3_composeTranspose_ra(xx + 145, xx + 109, xx + 154);
  pm_math_Matrix3x3_compose_ra(xx + 109, xx + 154, xx + 145);
  pm_math_Matrix3x3_postCross_ra(xx + 145, xx + 74, xx + 109);
  pm_math_Matrix3x3_preCross_ra(xx + 109, xx + 74, xx + 154);
  xx[19] = xx[13] + xx[118] - xx[136] - xx[136] - xx[154];
  xx[13] = xx[12] * state[3];
  xx[29] = xx[119] - xx[137] - xx[139] - xx[155];
  xx[31] = xx[127] - xx[109];
  xx[35] = xx[128] - xx[112];
  xx[38] = xx[129] - xx[115];
  xx[39] = xx[130] - xx[110];
  xx[40] = xx[131] - xx[113];
  xx[47] = xx[132] - xx[116];
  xx[49] = xx[133] - xx[111];
  xx[52] = xx[134] - xx[114];
  xx[54] = xx[135] - xx[117];
  xx[109] = xx[31];
  xx[110] = xx[35];
  xx[111] = xx[38];
  xx[112] = xx[39];
  xx[113] = xx[40];
  xx[114] = xx[47];
  xx[115] = xx[49];
  xx[116] = xx[52];
  xx[117] = xx[54];
  xx[55] = xx[1] * xx[4];
  xx[56] = xx[55] * xx[3];
  xx[61] = xx[9] * (xx[56] - xx[56]);
  xx[63] = xx[61] * state[1] * state[1];
  xx[64] = xx[4] * xx[63];
  xx[79] = xx[3];
  xx[80] = xx[2];
  xx[81] = xx[4];
  xx[2] = xx[55] * xx[4];
  xx[55] = xx[1] - (xx[9] * (xx[2] + xx[2]) - xx[1]);
  xx[2] = xx[55] * state[1] * state[1];
  xx[70] = xx[4] * xx[2];
  xx[73] = xx[3] * xx[63] - xx[70];
  xx[84] = - xx[64];
  xx[85] = - xx[70];
  xx[86] = xx[73];
  pm_math_Vector3_cross_ra(xx + 79, xx + 84, xx + 98);
  xx[82] = xx[1] * state[3];
  xx[84] = xx[9] * (xx[3] * xx[64] + xx[98]) - xx[2] - xx[82] * (xx[14] + xx[11]);
  xx[2] = xx[63] + xx[9] * (xx[3] * xx[70] + xx[99]);
  xx[11] = xx[82] * (xx[12] + xx[12]) + xx[9] * (xx[100] - xx[3] * xx[73]);
  xx[98] = xx[84];
  xx[99] = xx[2];
  xx[100] = xx[11];
  pm_math_Matrix3x3_xform_ra(xx + 109, xx + 98, xx + 127);
  xx[12] = xx[120] - xx[138] - xx[142] - xx[156];
  xx[14] = xx[12] + xx[1] * xx[35];
  xx[63] = xx[124] - xx[142] - xx[138] - xx[160];
  xx[64] = xx[125] - xx[143] - xx[141] - xx[161];
  xx[70] = xx[24] + xx[108] + xx[103] + xx[71] * xx[63] - xx[13] * xx[64] + xx
  [129];
  xx[73] = xx[66] + xx[149];
  xx[109] = xx[66] + xx[145];
  xx[110] = xx[146];
  xx[111] = xx[147];
  xx[112] = xx[148];
  xx[113] = xx[73];
  xx[114] = xx[150];
  xx[115] = xx[151];
  xx[116] = xx[152];
  xx[117] = xx[66] + xx[153];
  pm_math_Matrix3x3_xform_ra(xx + 109, xx + 98, xx + 130);
  xx[82] = xx[88] + xx[71] * xx[35] - xx[13] * xx[40] + xx[131];
  xx[85] = xx[18] + xx[126] - xx[144] - xx[144] - xx[162];
  xx[86] = xx[85] + xx[1] * xx[52];
  xx[88] = xx[52] + xx[1] * xx[73];
  xx[90] = xx[86] + xx[1] * xx[88];

  xx[94] = (input[1] - (xx[70] + xx[1] * xx[82])) / xx[90];
  xx[24] = xx[121] - xx[139] - xx[137] - xx[157];
  xx[98] = xx[18] + xx[122] - xx[140] - xx[140] - xx[158];
  xx[99] = xx[123] - xx[141] - xx[143] - xx[159];
  xx[100] = xx[99] + xx[1] * xx[40];
  xx[103] = xx[22] + xx[106] + xx[101] + xx[71] * xx[19] - xx[13] * xx[29] + xx
  [127] + xx[14] * xx[94];
  xx[104] = xx[23] + xx[107] + xx[102] + xx[71] * xx[24] - xx[13] * xx[98] + xx
  [128] + xx[100] * xx[94];
  xx[105] = xx[70] + xx[86] * xx[94];
  pm_math_Quaternion_xform_ra(xx + 5, xx + 103, xx + 106);
  xx[101] = xx[55];
  xx[102] = - xx[61];
  xx[103] = xx[9] * (xx[56] + xx[56]);
  xx[22] = xx[49] + xx[1] * xx[146];
  xx[23] = xx[54] + xx[1] * xx[152];
  xx[104] = xx[87] + xx[71] * xx[31] - xx[13] * xx[39] + xx[130] + xx[22] * xx
  [94];
  xx[105] = xx[82] + xx[88] * xx[94];
  xx[106] = xx[89] + xx[71] * xx[38] - xx[13] * xx[47] + xx[132] + xx[23] * xx
  [94];
  pm_math_Quaternion_xform_ra(xx + 5, xx + 104, xx + 109);
  pm_math_Vector3_cross_ra(xx + 101, xx + 109, xx + 5);
  xx[5] = xx[1] * state[1] * state[1];
  xx[6] = xx[3] * xx[3];
  xx[8] = xx[3] * xx[4];
  xx[56] = xx[9] * (xx[8] + xx[8]);
  xx[70] = xx[9] * (xx[8] - xx[8]);
  xx[8] = xx[4] * xx[4];
  xx[82] = xx[9] * (xx[6] + xx[8]) - xx[66];
  xx[111] = xx[9] * (xx[6] + xx[6]) - xx[66];
  xx[112] = - xx[56];
  xx[113] = xx[70];
  xx[114] = xx[70];
  xx[115] = xx[82];
  xx[116] = - (xx[9] * (xx[8] + xx[6]));
  xx[117] = xx[56];
  xx[118] = xx[9] * (xx[6] - xx[8]);
  xx[119] = xx[82];
  xx[6] = xx[22] / xx[90];
  xx[8] = xx[88] / xx[90];
  xx[56] = xx[23] / xx[90];
  xx[120] = xx[31] - xx[14] * xx[6];
  xx[121] = xx[35] - xx[14] * xx[8];
  xx[122] = xx[38] - xx[14] * xx[56];
  xx[123] = xx[39] - xx[100] * xx[6];
  xx[124] = xx[40] - xx[100] * xx[8];
  xx[125] = xx[47] - xx[100] * xx[56];
  xx[126] = xx[49] - xx[86] * xx[6];
  xx[127] = xx[52] - xx[86] * xx[8];
  xx[128] = xx[54] - xx[86] * xx[56];
  pm_math_Matrix3x3_composeTranspose_ra(xx + 120, xx + 111, xx + 129);
  pm_math_Matrix3x3_compose_ra(xx + 111, xx + 129, xx + 120);
  xx[31] = xx[88] * xx[6];
  xx[35] = xx[23] * xx[6];
  xx[38] = xx[23] * xx[8];
  xx[129] = xx[145] - xx[22] * xx[6] + xx[66];
  xx[130] = xx[146] - xx[31];
  xx[131] = xx[147] - xx[35];
  xx[132] = xx[148] - xx[31];
  xx[133] = xx[73] - xx[88] * xx[8];
  xx[134] = xx[150] - xx[38];
  xx[135] = xx[151] - xx[35];
  xx[136] = xx[152] - xx[38];
  xx[137] = xx[153] - xx[23] * xx[56] + xx[66];
  pm_math_Matrix3x3_composeTranspose_ra(xx + 129, xx + 111, xx + 138);
  pm_math_Matrix3x3_compose_ra(xx + 111, xx + 138, xx + 129);
  pm_math_Matrix3x3_postCross_ra(xx + 129, xx + 101, xx + 138);
  xx[22] = xx[126] - xx[140];
  xx[23] = xx[14] / xx[90];
  xx[31] = xx[100] * xx[23];
  xx[35] = xx[86] * xx[23];
  xx[38] = xx[100] / xx[90];
  xx[39] = xx[86] * xx[38];
  xx[40] = xx[86] / xx[90];
  xx[147] = xx[19] - xx[14] * xx[23];
  xx[148] = xx[29] - xx[31];
  xx[149] = xx[12] - xx[35];
  xx[150] = xx[24] - xx[31];
  xx[151] = xx[98] - xx[100] * xx[38];
  xx[152] = xx[99] - xx[39];
  xx[153] = xx[63] - xx[35];
  xx[154] = xx[64] - xx[39];
  xx[155] = xx[85] - xx[86] * xx[40];
  pm_math_Matrix3x3_composeTranspose_ra(xx + 147, xx + 111, xx + 156);
  pm_math_Matrix3x3_compose_ra(xx + 111, xx + 156, xx + 147);
  pm_math_Matrix3x3_postCross_ra(xx + 120, xx + 101, xx + 111);
  pm_math_Matrix3x3_preCross_ra(xx + 138, xx + 101, xx + 156);
  xx[12] = xx[127] - xx[143];
  xx[14] = xx[12] + xx[1] * (xx[66] + xx[133]);
  xx[19] = xx[155] - xx[119] - xx[119] - xx[164] + xx[1] * xx[12] + xx[1] * xx
  [14] + xx[18];

  xx[12] = 9.81;
  xx[18] = xx[1] * state[0];
  xx[24] = xx[0] * cos(xx[18]);
  xx[29] = xx[0] * sin(xx[18]);
  xx[0] = xx[24] - xx[29];
  xx[18] = xx[12] * xx[0];
  xx[31] = xx[12] - xx[9] * xx[18] * xx[0];
  xx[0] = xx[9] * xx[18] * (xx[24] + xx[29]);
  xx[12] = (input[0] - (xx[108] + xx[7] - xx[5] * xx[22] + xx[1] * (xx[110] -
                                                                    xx[5] * xx[132]))) / xx[19] - (xx[31] * xx[14] / xx[19] - xx[0] *
                                                                                                                              (xx[22] + xx[1] * xx[130]) / xx[19]);
  xx[85] = xx[23];
  xx[86] = xx[38];
  xx[87] = xx[40];
  xx[7] = xx[4] * xx[12];
  xx[14] = xx[3] * xx[12];
  xx[18] = xx[9] * (xx[3] * xx[7] + xx[4] * xx[14]);
  xx[19] = xx[3] * xx[14];
  xx[14] = xx[4] * xx[7];
  xx[7] = xx[9] * (xx[19] - xx[14]);
  xx[22] = xx[12] - xx[9] * (xx[19] + xx[14]);
  xx[38] = xx[18];
  xx[39] = xx[7];
  xx[40] = xx[22];
  xx[88] = xx[6];
  xx[89] = xx[8];
  xx[90] = xx[56];
  xx[6] = xx[61] * xx[12] - (xx[0] + xx[5]);
  xx[0] = xx[31] + xx[1] * xx[12] + xx[55] * xx[12];
  xx[5] = xx[4] * xx[0];
  xx[8] = xx[4] * xx[6];
  xx[4] = xx[3] * xx[0] + xx[8];
  xx[54] = - xx[5];
  xx[55] = xx[8];
  xx[56] = xx[4];
  pm_math_Vector3_cross_ra(xx + 79, xx + 54, xx + 98);
  xx[14] = xx[6] + xx[9] * (xx[3] * xx[5] + xx[98]);
  xx[5] = xx[0] + xx[9] * (xx[99] - xx[3] * xx[8]);
  xx[0] = xx[9] * (xx[100] - xx[3] * xx[4]);
  xx[54] = xx[14];
  xx[55] = xx[5];
  xx[56] = xx[0];
  xx[3] = xx[94] - (pm_math_Vector3_dot_ra(xx + 85, xx + 38) +
                    pm_math_Vector3_dot_ra(xx + 88, xx + 54));
  xx[38] = xx[10];
  xx[39] = xx[77];
  xx[40] = xx[83];
  xx[8] = xx[18] + xx[71];
  xx[9] = xx[7] - xx[13];
  xx[10] = xx[22] + xx[3];
  pm_math_Quaternion_inverseXform_ra(xx + 25, xx + 8, xx + 22);
  xx[54] = xx[17];
  xx[55] = xx[20];
  xx[56] = xx[21];
  pm_math_Vector3_cross_ra(xx + 8, xx + 74, xx + 17);
  xx[6] = xx[14] + xx[84] + xx[17];
  xx[7] = xx[5] + xx[1] * xx[3] + xx[2] + xx[18];
  xx[8] = xx[0] + xx[11] + xx[19];
  pm_math_Quaternion_inverseXform_ra(xx + 25, xx + 6, xx + 9);
  xx[0] = xx[91] - (pm_math_Vector3_dot_ra(xx + 38, xx + 22) +
                    pm_math_Vector3_dot_ra(xx + 54, xx + 9));
  xx[4] = xx[30];
  xx[5] = xx[50];
  xx[6] = xx[62];
  xx[17] = xx[22] + xx[48];
  xx[18] = xx[23] - xx[45];
  xx[19] = xx[24] + xx[0];
  pm_math_Quaternion_inverseXform_ra(xx + 41, xx + 17, xx + 20);
  xx[23] = xx[34];
  xx[24] = xx[36];
  xx[25] = xx[37];
  pm_math_Vector3_cross_ra(xx + 17, xx + 95, xx + 26);
  xx[17] = xx[9] + xx[72] + xx[26];
  xx[18] = xx[10] + xx[1] * xx[0] + xx[16] + xx[27];
  xx[19] = xx[11] + xx[15] + xx[28];
  pm_math_Quaternion_inverseXform_ra(xx + 41, xx + 17, xx + 7);
  xx[2] = xx[92] - (pm_math_Vector3_dot_ra(xx + 4, xx + 20) +
                    pm_math_Vector3_dot_ra(xx + 23, xx + 7));
  xx[4] = xx[20] + xx[51];
  xx[5] = xx[21] - xx[53];
  xx[6] = xx[22] + xx[2];
  pm_math_Quaternion_inverseXform_ra(xx + 57, xx + 4, xx + 13);
  pm_math_Vector3_cross_ra(xx + 4, xx + 67, xx + 16);
  xx[4] = xx[7] + xx[93] + xx[16];
  xx[5] = xx[8] + xx[1] * xx[2] + xx[33] + xx[17];
  xx[6] = xx[9] + xx[32] + xx[18];
  pm_math_Quaternion_inverseXform_ra(xx + 57, xx + 4, xx + 7);
  deriv[0] = state[1];
  deriv[1] = xx[12];
  deriv[2] = state[3];
  deriv[3] = xx[3];
  deriv[4] = state[5];
  deriv[5] = xx[0];
  deriv[6] = state[7];
  deriv[7] = xx[2];
  deriv[8] = state[9];
  deriv[9] = xx[65] - (xx[78] * xx[15] + xx[46] * xx[8]);

//  for (int i=0; i<7; ++i)
//  {
//    deriv[2*i] = normalize_angle_positive(deriv[2*i]);
//  }

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
}

