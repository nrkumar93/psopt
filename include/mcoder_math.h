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
 * \file   mcoder_math.h
 * \author Ramkumar Natarajan (rnataraj@cs.cmu.edu)
 * \date   6/22/21
 */

#ifndef PSOPT_MCODER_MATH_H
#define PSOPT_MCODER_MATH_H

//#include "pm_std.h"
#include "tmwtypes.h"
#include <cstring>
#include <cmath>

//template <class T> T pm_math_Vector3_dot_ra(const T*pm_math__AuaKMC5koOPbXxPi_MZvt,
//const T*pm_math_kNtPmLll5l8eiqwtk5kfJv);
//
//template <class T> void pm_math_Vector3_cross_ra(
//const T*pm_math__AuaKMC5koOPbXxPi_MZvt,const T*
//pm_math_kNtPmLll5l8eiqwtk5kfJv,T*pm_math__1Zf2IciMRCub1vvbEr1C4);
//
//template <class T> void pm_math_Vector3_compOrthogonalBasis_ra(const T*
//pm_math_VgJW5ZqpwPpuY1inYtaofQ,T*pm_math_kwrB3ZoKf7OufTHWaHJV7a,T*
//pm_math_kyp6uAyJE40UVuAQNEYzS1,T*pm_math_V2__YrimeI4E_yWnhKofpy);
//
//template <class T> void pm_math_Quaternion_compose_ra(const T*pm_math_FbCdebDtDIhMduMKWA8Khy,
//const T*pm_math_krgkQdg3ZZ_0ded42Fk_8r,T*
//pm_math__1Zf2IciMRCub1vvbEr1C4);
//
//template <class T> void pm_math_Quaternion_composeInverse_ra(
//const T*pm_math_FbCdebDtDIhMduMKWA8Khy,const T*
//pm_math_krgkQdg3ZZ_0ded42Fk_8r,T*pm_math__1Zf2IciMRCub1vvbEr1C4);
//
//template <class T> void pm_math_Quaternion_inverseCompose_ra(const T*
//pm_math_FbCdebDtDIhMduMKWA8Khy,const T*pm_math_krgkQdg3ZZ_0ded42Fk_8r,
//T*pm_math__1Zf2IciMRCub1vvbEr1C4);
//
//template <class T> void pm_math_Quaternion_xform_ra(const
//T*pm_math_VnD_HGFKVUOWdDAWvZhyEb,const T*
//pm_math_VgJW5ZqpwPpuY1inYtaofQ,T*pm_math__1Zf2IciMRCub1vvbEr1C4);
//
//template <class T> void pm_math_Quaternion_inverseXform_ra(const T*pm_math_VnD_HGFKVUOWdDAWvZhyEb
//,const T*pm_math_VgJW5ZqpwPpuY1inYtaofQ,T*
//pm_math__1Zf2IciMRCub1vvbEr1C4);
//
//template <class T> void pm_math_Quaternion_compDeriv_ra(const
//T*pm_math_VnD_HGFKVUOWdDAWvZhyEb,const T*
//pm_math_kvL2QWFrSblY_eznRgCz87,T*pm_math_FGJfoQhjnH_kdaWooDAepH);
//
//template <class T> void pm_math_Quaternion_Matrix3x3Ctor_ra(const T*
//pm_math_FqUCZrSGGNOuePgRr82o_8,T*pm_math_VnD_HGFKVUOWdDAWvZhyEb);
//
//template <class T> void pm_math_Matrix3x3_compose_ra(const T*pm_math_Fcuud3IN0odMZi54a1R_8f,const
//T*pm_math__09m2ugY6U_OXH9Il_a7Bj,T*pm_math__1Zf2IciMRCub1vvbEr1C4);
//
//template <class T> void pm_math_Matrix3x3_composeTranspose_ra(const T*
//pm_math_Fcuud3IN0odMZi54a1R_8f,const T*pm_math__09m2ugY6U_OXH9Il_a7Bj,
//T*pm_math__1Zf2IciMRCub1vvbEr1C4);
//
//template <class T> void pm_math_Matrix3x3_transposeCompose_ra(const T*
//pm_math_Fcuud3IN0odMZi54a1R_8f,const T*pm_math__09m2ugY6U_OXH9Il_a7Bj,
//T*pm_math__1Zf2IciMRCub1vvbEr1C4);
//
//template <class T> void pm_math_Matrix3x3_preCross_ra(
//const T*pm_math_F2l4p_g4sn02huHNflQjMH,const T*
//pm_math_VgJW5ZqpwPpuY1inYtaofQ,T*pm_math__1Zf2IciMRCub1vvbEr1C4);
//
//template <class T> void pm_math_Matrix3x3_postCross_ra(const T*pm_math_F2l4p_g4sn02huHNflQjMH,
//const T*pm_math_VgJW5ZqpwPpuY1inYtaofQ,T*
//pm_math__1Zf2IciMRCub1vvbEr1C4);
//
//template <class T> void pm_math_Matrix3x3_xform_ra(const T*
//pm_math_F2l4p_g4sn02huHNflQjMH,const T*pm_math_VgJW5ZqpwPpuY1inYtaofQ,
//T*pm_math__1Zf2IciMRCub1vvbEr1C4);
//
//template <class T> void pm_math_Matrix3x3_transposeXform_ra
//(const T*pm_math_F2l4p_g4sn02huHNflQjMH,const T*
//pm_math_VgJW5ZqpwPpuY1inYtaofQ,T*pm_math__1Zf2IciMRCub1vvbEr1C4);
//
//template <class T> void pm_math_Matrix3x3_minRotation_ra(const T*pm_math_VnD_HGFKVUOWdDAWvZhyEb,
//T*pm_math_FqUCZrSGGNOuePgRr82o_8,int pm_math_kbzF46WM0FtKWeE0s7WR95[3]);

template <class T> T pm_math_Vector3_dot_ra(const T*pm_math__AuaKMC5koOPbXxPi_MZvt,
const T*pm_math_kNtPmLll5l8eiqwtk5kfJv)
{
  return pm_math__AuaKMC5koOPbXxPi_MZvt[0]*pm_math_kNtPmLll5l8eiqwtk5kfJv[0]+
pm_math__AuaKMC5koOPbXxPi_MZvt[1]*pm_math_kNtPmLll5l8eiqwtk5kfJv[1]+
pm_math__AuaKMC5koOPbXxPi_MZvt[2]*pm_math_kNtPmLll5l8eiqwtk5kfJv[2];
}

template <class T> void pm_math_Vector3_cross_ra(const T*pm_math__AuaKMC5koOPbXxPi_MZvt,const
T*pm_math_kNtPmLll5l8eiqwtk5kfJv,T*pm_math__1Zf2IciMRCub1vvbEr1C4)
{
  (void)0;;*pm_math__1Zf2IciMRCub1vvbEr1C4++=pm_math__AuaKMC5koOPbXxPi_MZvt[1]*
pm_math_kNtPmLll5l8eiqwtk5kfJv[2]-pm_math__AuaKMC5koOPbXxPi_MZvt[2]*
pm_math_kNtPmLll5l8eiqwtk5kfJv[1];*pm_math__1Zf2IciMRCub1vvbEr1C4++=
pm_math__AuaKMC5koOPbXxPi_MZvt[2]*pm_math_kNtPmLll5l8eiqwtk5kfJv[0]-
pm_math__AuaKMC5koOPbXxPi_MZvt[0]*pm_math_kNtPmLll5l8eiqwtk5kfJv[2];*
pm_math__1Zf2IciMRCub1vvbEr1C4=pm_math__AuaKMC5koOPbXxPi_MZvt[0]*
pm_math_kNtPmLll5l8eiqwtk5kfJv[1]-pm_math__AuaKMC5koOPbXxPi_MZvt[1]*
pm_math_kNtPmLll5l8eiqwtk5kfJv[0];
}

template <class T> static void pm_math_FvxVpUGZsRt_iyppjQhBrA(
const T*pm_math_VgJW5ZqpwPpuY1inYtaofQ,T*
pm_math_VhMZ7Balz3hebuYIOnZcwH)
{
  const T pm_math_Fkbrfbd5DA0yb1kdOaIHtl=
sqrt(pm_math_VgJW5ZqpwPpuY1inYtaofQ[0]*pm_math_VgJW5ZqpwPpuY1inYtaofQ[0]+
pm_math_VgJW5ZqpwPpuY1inYtaofQ[1]*pm_math_VgJW5ZqpwPpuY1inYtaofQ[1]+
pm_math_VgJW5ZqpwPpuY1inYtaofQ[2]*pm_math_VgJW5ZqpwPpuY1inYtaofQ[2]);
pm_math_VhMZ7Balz3hebuYIOnZcwH[0]=pm_math_VgJW5ZqpwPpuY1inYtaofQ[0]/
pm_math_Fkbrfbd5DA0yb1kdOaIHtl;pm_math_VhMZ7Balz3hebuYIOnZcwH[1]=
pm_math_VgJW5ZqpwPpuY1inYtaofQ[1]/pm_math_Fkbrfbd5DA0yb1kdOaIHtl;
pm_math_VhMZ7Balz3hebuYIOnZcwH[2]=pm_math_VgJW5ZqpwPpuY1inYtaofQ[2]/
pm_math_Fkbrfbd5DA0yb1kdOaIHtl;
}

template <class T> void pm_math_Vector3_compOrthogonalBasis_ra(
const T*pm_math_VgJW5ZqpwPpuY1inYtaofQ,T*
pm_math_kwrB3ZoKf7OufTHWaHJV7a,T*pm_math_kyp6uAyJE40UVuAQNEYzS1,T*
pm_math_V2__YrimeI4E_yWnhKofpy)
{
  T pm_math_Fklj3odRoMCGhX9VnS2RBz,
pm_math___TAiJ__Y1OQWD7mG_iJLJ,pm_math_VCj4vmLqKyOyie7D24jNe6;T
pm_math_FEa9HpQub3K9X1Fkj0Fvl_[3]={0.0,0.0,0.0};

  pm_math_FvxVpUGZsRt_iyppjQhBrA(pm_math_VgJW5ZqpwPpuY1inYtaofQ,pm_math_V2__YrimeI4E_yWnhKofpy);
  pm_math_Fklj3odRoMCGhX9VnS2RBz=fabs(pm_math_V2__YrimeI4E_yWnhKofpy[0]);
  pm_math___TAiJ__Y1OQWD7mG_iJLJ=fabs(pm_math_V2__YrimeI4E_yWnhKofpy[1]);
  pm_math_VCj4vmLqKyOyie7D24jNe6=fabs(pm_math_V2__YrimeI4E_yWnhKofpy[2]);

{
  const int pm_math_k0_sG8vQ8u0Bc9IThf_p20=(pm_math_Fklj3odRoMCGhX9VnS2RBz>=
pm_math___TAiJ__Y1OQWD7mG_iJLJ)?((pm_math_Fklj3odRoMCGhX9VnS2RBz>
pm_math_VCj4vmLqKyOyie7D24jNe6)?1:0):((pm_math___TAiJ__Y1OQWD7mG_iJLJ>
pm_math_VCj4vmLqKyOyie7D24jNe6)?2:0);pm_math_FEa9HpQub3K9X1Fkj0Fvl_[
pm_math_k0_sG8vQ8u0Bc9IThf_p20]=1.0;
}

  pm_math_Vector3_cross_ra(pm_math_V2__YrimeI4E_yWnhKofpy,pm_math_FEa9HpQub3K9X1Fkj0Fvl_,
pm_math_kyp6uAyJE40UVuAQNEYzS1);

  pm_math_FvxVpUGZsRt_iyppjQhBrA(
pm_math_kyp6uAyJE40UVuAQNEYzS1,pm_math_kyp6uAyJE40UVuAQNEYzS1);

  pm_math_Vector3_cross_ra(pm_math_kyp6uAyJE40UVuAQNEYzS1,
pm_math_V2__YrimeI4E_yWnhKofpy,pm_math_kwrB3ZoKf7OufTHWaHJV7a);

  pm_math_FvxVpUGZsRt_iyppjQhBrA(pm_math_kwrB3ZoKf7OufTHWaHJV7a,
pm_math_kwrB3ZoKf7OufTHWaHJV7a);
}

template <class T> void pm_math_Quaternion_compose_ra(const
T*pm_math_FbCdebDtDIhMduMKWA8Khy,const T*
pm_math_krgkQdg3ZZ_0ded42Fk_8r,T*pm_math__1Zf2IciMRCub1vvbEr1C4)
{
  (void)0;
;*pm_math__1Zf2IciMRCub1vvbEr1C4++=pm_math_FbCdebDtDIhMduMKWA8Khy[0]*
pm_math_krgkQdg3ZZ_0ded42Fk_8r[0]-(pm_math_FbCdebDtDIhMduMKWA8Khy[1]*
pm_math_krgkQdg3ZZ_0ded42Fk_8r[1]+pm_math_FbCdebDtDIhMduMKWA8Khy[2]*
pm_math_krgkQdg3ZZ_0ded42Fk_8r[2]+pm_math_FbCdebDtDIhMduMKWA8Khy[3]*
pm_math_krgkQdg3ZZ_0ded42Fk_8r[3]);*pm_math__1Zf2IciMRCub1vvbEr1C4++=
pm_math_FbCdebDtDIhMduMKWA8Khy[0]*pm_math_krgkQdg3ZZ_0ded42Fk_8r[1]+
pm_math_FbCdebDtDIhMduMKWA8Khy[1]*pm_math_krgkQdg3ZZ_0ded42Fk_8r[0]+
pm_math_FbCdebDtDIhMduMKWA8Khy[2]*pm_math_krgkQdg3ZZ_0ded42Fk_8r[3]-
pm_math_FbCdebDtDIhMduMKWA8Khy[3]*pm_math_krgkQdg3ZZ_0ded42Fk_8r[2];*
pm_math__1Zf2IciMRCub1vvbEr1C4++=pm_math_FbCdebDtDIhMduMKWA8Khy[0]*
pm_math_krgkQdg3ZZ_0ded42Fk_8r[2]+pm_math_FbCdebDtDIhMduMKWA8Khy[2]*
pm_math_krgkQdg3ZZ_0ded42Fk_8r[0]+pm_math_FbCdebDtDIhMduMKWA8Khy[3]*
pm_math_krgkQdg3ZZ_0ded42Fk_8r[1]-pm_math_FbCdebDtDIhMduMKWA8Khy[1]*
pm_math_krgkQdg3ZZ_0ded42Fk_8r[3];*pm_math__1Zf2IciMRCub1vvbEr1C4=
pm_math_FbCdebDtDIhMduMKWA8Khy[0]*pm_math_krgkQdg3ZZ_0ded42Fk_8r[3]+
pm_math_FbCdebDtDIhMduMKWA8Khy[3]*pm_math_krgkQdg3ZZ_0ded42Fk_8r[0]+
pm_math_FbCdebDtDIhMduMKWA8Khy[1]*pm_math_krgkQdg3ZZ_0ded42Fk_8r[2]-
pm_math_FbCdebDtDIhMduMKWA8Khy[2]*pm_math_krgkQdg3ZZ_0ded42Fk_8r[1];
}

template <class T> void pm_math_Quaternion_composeInverse_ra(const T*
pm_math_FbCdebDtDIhMduMKWA8Khy,const T*pm_math_krgkQdg3ZZ_0ded42Fk_8r,
T*pm_math__1Zf2IciMRCub1vvbEr1C4)
{
  (void)0;;*
pm_math__1Zf2IciMRCub1vvbEr1C4++=pm_math_FbCdebDtDIhMduMKWA8Khy[0]*-
pm_math_krgkQdg3ZZ_0ded42Fk_8r[0]-(pm_math_FbCdebDtDIhMduMKWA8Khy[1]*
pm_math_krgkQdg3ZZ_0ded42Fk_8r[1]+pm_math_FbCdebDtDIhMduMKWA8Khy[2]*
pm_math_krgkQdg3ZZ_0ded42Fk_8r[2]+pm_math_FbCdebDtDIhMduMKWA8Khy[3]*
pm_math_krgkQdg3ZZ_0ded42Fk_8r[3]);*pm_math__1Zf2IciMRCub1vvbEr1C4++=
pm_math_FbCdebDtDIhMduMKWA8Khy[0]*pm_math_krgkQdg3ZZ_0ded42Fk_8r[1]-
pm_math_FbCdebDtDIhMduMKWA8Khy[1]*pm_math_krgkQdg3ZZ_0ded42Fk_8r[0]+
pm_math_FbCdebDtDIhMduMKWA8Khy[2]*pm_math_krgkQdg3ZZ_0ded42Fk_8r[3]-
pm_math_FbCdebDtDIhMduMKWA8Khy[3]*pm_math_krgkQdg3ZZ_0ded42Fk_8r[2];*
pm_math__1Zf2IciMRCub1vvbEr1C4++=pm_math_FbCdebDtDIhMduMKWA8Khy[0]*
pm_math_krgkQdg3ZZ_0ded42Fk_8r[2]-pm_math_FbCdebDtDIhMduMKWA8Khy[2]*
pm_math_krgkQdg3ZZ_0ded42Fk_8r[0]+pm_math_FbCdebDtDIhMduMKWA8Khy[3]*
pm_math_krgkQdg3ZZ_0ded42Fk_8r[1]-pm_math_FbCdebDtDIhMduMKWA8Khy[1]*
pm_math_krgkQdg3ZZ_0ded42Fk_8r[3];*pm_math__1Zf2IciMRCub1vvbEr1C4=
pm_math_FbCdebDtDIhMduMKWA8Khy[0]*pm_math_krgkQdg3ZZ_0ded42Fk_8r[3]-
pm_math_FbCdebDtDIhMduMKWA8Khy[3]*pm_math_krgkQdg3ZZ_0ded42Fk_8r[0]+
pm_math_FbCdebDtDIhMduMKWA8Khy[1]*pm_math_krgkQdg3ZZ_0ded42Fk_8r[2]-
pm_math_FbCdebDtDIhMduMKWA8Khy[2]*pm_math_krgkQdg3ZZ_0ded42Fk_8r[1];
}

template <class T> void pm_math_Quaternion_inverseCompose_ra(const T*
pm_math_FbCdebDtDIhMduMKWA8Khy,const T*pm_math_krgkQdg3ZZ_0ded42Fk_8r,
T*pm_math__1Zf2IciMRCub1vvbEr1C4)
{
  (void)0;;*
pm_math__1Zf2IciMRCub1vvbEr1C4++= -pm_math_FbCdebDtDIhMduMKWA8Khy[0]*
pm_math_krgkQdg3ZZ_0ded42Fk_8r[0]-(pm_math_FbCdebDtDIhMduMKWA8Khy[1]*
pm_math_krgkQdg3ZZ_0ded42Fk_8r[1]+pm_math_FbCdebDtDIhMduMKWA8Khy[2]*
pm_math_krgkQdg3ZZ_0ded42Fk_8r[2]+pm_math_FbCdebDtDIhMduMKWA8Khy[3]*
pm_math_krgkQdg3ZZ_0ded42Fk_8r[3]);*pm_math__1Zf2IciMRCub1vvbEr1C4++= -
pm_math_FbCdebDtDIhMduMKWA8Khy[0]*pm_math_krgkQdg3ZZ_0ded42Fk_8r[1]+
pm_math_FbCdebDtDIhMduMKWA8Khy[1]*pm_math_krgkQdg3ZZ_0ded42Fk_8r[0]+
pm_math_FbCdebDtDIhMduMKWA8Khy[2]*pm_math_krgkQdg3ZZ_0ded42Fk_8r[3]-
pm_math_FbCdebDtDIhMduMKWA8Khy[3]*pm_math_krgkQdg3ZZ_0ded42Fk_8r[2];*
pm_math__1Zf2IciMRCub1vvbEr1C4++= -pm_math_FbCdebDtDIhMduMKWA8Khy[0]*
pm_math_krgkQdg3ZZ_0ded42Fk_8r[2]+pm_math_FbCdebDtDIhMduMKWA8Khy[2]*
pm_math_krgkQdg3ZZ_0ded42Fk_8r[0]+pm_math_FbCdebDtDIhMduMKWA8Khy[3]*
pm_math_krgkQdg3ZZ_0ded42Fk_8r[1]-pm_math_FbCdebDtDIhMduMKWA8Khy[1]*
pm_math_krgkQdg3ZZ_0ded42Fk_8r[3];*pm_math__1Zf2IciMRCub1vvbEr1C4= -
pm_math_FbCdebDtDIhMduMKWA8Khy[0]*pm_math_krgkQdg3ZZ_0ded42Fk_8r[3]+
pm_math_FbCdebDtDIhMduMKWA8Khy[3]*pm_math_krgkQdg3ZZ_0ded42Fk_8r[0]+
pm_math_FbCdebDtDIhMduMKWA8Khy[1]*pm_math_krgkQdg3ZZ_0ded42Fk_8r[2]-
pm_math_FbCdebDtDIhMduMKWA8Khy[2]*pm_math_krgkQdg3ZZ_0ded42Fk_8r[1];
}

template <class T> void pm_math_Quaternion_xform_ra(const T*pm_math_VnD_HGFKVUOWdDAWvZhyEb,const
T*pm_math_VgJW5ZqpwPpuY1inYtaofQ,T*pm_math__1Zf2IciMRCub1vvbEr1C4)
{
T pm_math_FFhjXMhOoO8gZHlh8Hu3K3[3],pm_math_FEcK0zAOzXKSZq2OQse79O[3];
pm_math_Vector3_cross_ra(pm_math_VnD_HGFKVUOWdDAWvZhyEb+1,
pm_math_VgJW5ZqpwPpuY1inYtaofQ,pm_math_FFhjXMhOoO8gZHlh8Hu3K3);
pm_math_Vector3_cross_ra(pm_math_VnD_HGFKVUOWdDAWvZhyEb+1,
pm_math_FFhjXMhOoO8gZHlh8Hu3K3,pm_math_FEcK0zAOzXKSZq2OQse79O);*
pm_math__1Zf2IciMRCub1vvbEr1C4++=pm_math_VgJW5ZqpwPpuY1inYtaofQ[0]+2.0*(
pm_math_VnD_HGFKVUOWdDAWvZhyEb[0]*pm_math_FFhjXMhOoO8gZHlh8Hu3K3[0]+
pm_math_FEcK0zAOzXKSZq2OQse79O[0]);*pm_math__1Zf2IciMRCub1vvbEr1C4++=
pm_math_VgJW5ZqpwPpuY1inYtaofQ[1]+2.0*(pm_math_VnD_HGFKVUOWdDAWvZhyEb[0]*
pm_math_FFhjXMhOoO8gZHlh8Hu3K3[1]+pm_math_FEcK0zAOzXKSZq2OQse79O[1]);*
pm_math__1Zf2IciMRCub1vvbEr1C4=pm_math_VgJW5ZqpwPpuY1inYtaofQ[2]+2.0*(
pm_math_VnD_HGFKVUOWdDAWvZhyEb[0]*pm_math_FFhjXMhOoO8gZHlh8Hu3K3[2]+
pm_math_FEcK0zAOzXKSZq2OQse79O[2]);
}

template <class T> void pm_math_Quaternion_inverseXform_ra(
const T*pm_math_VnD_HGFKVUOWdDAWvZhyEb,const T*
pm_math_VgJW5ZqpwPpuY1inYtaofQ,T*pm_math__1Zf2IciMRCub1vvbEr1C4)
{
  T pm_math_FFhjXMhOoO8gZHlh8Hu3K3[3],pm_math_FEcK0zAOzXKSZq2OQse79O[3];
pm_math_Vector3_cross_ra(pm_math_VnD_HGFKVUOWdDAWvZhyEb+1,
pm_math_VgJW5ZqpwPpuY1inYtaofQ,pm_math_FFhjXMhOoO8gZHlh8Hu3K3);
pm_math_Vector3_cross_ra(pm_math_VnD_HGFKVUOWdDAWvZhyEb+1,
pm_math_FFhjXMhOoO8gZHlh8Hu3K3,pm_math_FEcK0zAOzXKSZq2OQse79O);*
pm_math__1Zf2IciMRCub1vvbEr1C4++=pm_math_VgJW5ZqpwPpuY1inYtaofQ[0]+2.0*(-
pm_math_VnD_HGFKVUOWdDAWvZhyEb[0]*pm_math_FFhjXMhOoO8gZHlh8Hu3K3[0]+
pm_math_FEcK0zAOzXKSZq2OQse79O[0]);*pm_math__1Zf2IciMRCub1vvbEr1C4++=
pm_math_VgJW5ZqpwPpuY1inYtaofQ[1]+2.0*(-pm_math_VnD_HGFKVUOWdDAWvZhyEb[0]*
pm_math_FFhjXMhOoO8gZHlh8Hu3K3[1]+pm_math_FEcK0zAOzXKSZq2OQse79O[1]);*
pm_math__1Zf2IciMRCub1vvbEr1C4=pm_math_VgJW5ZqpwPpuY1inYtaofQ[2]+2.0*(-
pm_math_VnD_HGFKVUOWdDAWvZhyEb[0]*pm_math_FFhjXMhOoO8gZHlh8Hu3K3[2]+
pm_math_FEcK0zAOzXKSZq2OQse79O[2]);
}

template <class T> void pm_math_Quaternion_compDeriv_ra(const
T*pm_math_VnD_HGFKVUOWdDAWvZhyEb,const T*
pm_math_kvL2QWFrSblY_eznRgCz87,T*pm_math_FGJfoQhjnH_kdaWooDAepH)
{
  T
pm_math_kgCuYXt6a6hngHkW2WPaRv[3];pm_math_kgCuYXt6a6hngHkW2WPaRv[0]=0.5*
pm_math_kvL2QWFrSblY_eznRgCz87[0];pm_math_kgCuYXt6a6hngHkW2WPaRv[1]=0.5*
pm_math_kvL2QWFrSblY_eznRgCz87[1];pm_math_kgCuYXt6a6hngHkW2WPaRv[2]=0.5*
pm_math_kvL2QWFrSblY_eznRgCz87[2];*pm_math_FGJfoQhjnH_kdaWooDAepH++= -
pm_math_VnD_HGFKVUOWdDAWvZhyEb[1]*pm_math_kgCuYXt6a6hngHkW2WPaRv[0]-
pm_math_VnD_HGFKVUOWdDAWvZhyEb[2]*pm_math_kgCuYXt6a6hngHkW2WPaRv[1]-
pm_math_VnD_HGFKVUOWdDAWvZhyEb[3]*pm_math_kgCuYXt6a6hngHkW2WPaRv[2];*
pm_math_FGJfoQhjnH_kdaWooDAepH++= +pm_math_VnD_HGFKVUOWdDAWvZhyEb[0]*
pm_math_kgCuYXt6a6hngHkW2WPaRv[0]-pm_math_VnD_HGFKVUOWdDAWvZhyEb[3]*
pm_math_kgCuYXt6a6hngHkW2WPaRv[1]+pm_math_VnD_HGFKVUOWdDAWvZhyEb[2]*
pm_math_kgCuYXt6a6hngHkW2WPaRv[2];*pm_math_FGJfoQhjnH_kdaWooDAepH++= +
pm_math_VnD_HGFKVUOWdDAWvZhyEb[3]*pm_math_kgCuYXt6a6hngHkW2WPaRv[0]+
pm_math_VnD_HGFKVUOWdDAWvZhyEb[0]*pm_math_kgCuYXt6a6hngHkW2WPaRv[1]-
pm_math_VnD_HGFKVUOWdDAWvZhyEb[1]*pm_math_kgCuYXt6a6hngHkW2WPaRv[2];*
pm_math_FGJfoQhjnH_kdaWooDAepH= -pm_math_VnD_HGFKVUOWdDAWvZhyEb[2]*
pm_math_kgCuYXt6a6hngHkW2WPaRv[0]+pm_math_VnD_HGFKVUOWdDAWvZhyEb[1]*
pm_math_kgCuYXt6a6hngHkW2WPaRv[1]+pm_math_VnD_HGFKVUOWdDAWvZhyEb[0]*
pm_math_kgCuYXt6a6hngHkW2WPaRv[2];
}

template <class T> void pm_math_Quaternion_Matrix3x3Ctor_ra(
const T*pm_math_FqUCZrSGGNOuePgRr82o_8,T*
pm_math_VnD_HGFKVUOWdDAWvZhyEb)
{
  const T pm_math_VmbGGnawYJ_FY5eg3Q9VOr=
0.25*(pm_math_FqUCZrSGGNOuePgRr82o_8[0]+pm_math_FqUCZrSGGNOuePgRr82o_8[4]+
pm_math_FqUCZrSGGNOuePgRr82o_8[8]+1.0);const T
pm_math_VQ_T199jNSCPYuRLcxtO6j=pm_math_VmbGGnawYJ_FY5eg3Q9VOr-0.5*(
pm_math_FqUCZrSGGNOuePgRr82o_8[4]+pm_math_FqUCZrSGGNOuePgRr82o_8[8]);const
T pm_math_VkE6m1GxptxxYDqNBGUEZw=pm_math_VmbGGnawYJ_FY5eg3Q9VOr-0.5*(
pm_math_FqUCZrSGGNOuePgRr82o_8[8]+pm_math_FqUCZrSGGNOuePgRr82o_8[0]);const
T pm_math_VcAmFCMQr_KodumW__UKW6=pm_math_VmbGGnawYJ_FY5eg3Q9VOr-0.5*(
pm_math_FqUCZrSGGNOuePgRr82o_8[0]+pm_math_FqUCZrSGGNOuePgRr82o_8[4]);T
pm_math_kEbBObcYFIxUZ5_77V3CO_=0.0;const int pm_math_kH3bO_YJW98YWa3J39Vj1I=(
pm_math_VmbGGnawYJ_FY5eg3Q9VOr>pm_math_VQ_T199jNSCPYuRLcxtO6j)?((
pm_math_VmbGGnawYJ_FY5eg3Q9VOr>pm_math_VkE6m1GxptxxYDqNBGUEZw)?((
pm_math_VmbGGnawYJ_FY5eg3Q9VOr>pm_math_VcAmFCMQr_KodumW__UKW6)?0:3):((
pm_math_VkE6m1GxptxxYDqNBGUEZw>pm_math_VcAmFCMQr_KodumW__UKW6)?2:3)):((
pm_math_VQ_T199jNSCPYuRLcxtO6j>pm_math_VkE6m1GxptxxYDqNBGUEZw)?((
pm_math_VQ_T199jNSCPYuRLcxtO6j>pm_math_VcAmFCMQr_KodumW__UKW6)?1:3):((
pm_math_VkE6m1GxptxxYDqNBGUEZw>pm_math_VcAmFCMQr_KodumW__UKW6)?2:3));switch(
pm_math_kH3bO_YJW98YWa3J39Vj1I){case 0:pm_math_VnD_HGFKVUOWdDAWvZhyEb[0]=sqrt(
pm_math_VmbGGnawYJ_FY5eg3Q9VOr);pm_math_kEbBObcYFIxUZ5_77V3CO_=0.25/
pm_math_VnD_HGFKVUOWdDAWvZhyEb[0];pm_math_VnD_HGFKVUOWdDAWvZhyEb[1]=(
pm_math_FqUCZrSGGNOuePgRr82o_8[7]-pm_math_FqUCZrSGGNOuePgRr82o_8[5])*
pm_math_kEbBObcYFIxUZ5_77V3CO_;pm_math_VnD_HGFKVUOWdDAWvZhyEb[2]=(
pm_math_FqUCZrSGGNOuePgRr82o_8[2]-pm_math_FqUCZrSGGNOuePgRr82o_8[6])*
pm_math_kEbBObcYFIxUZ5_77V3CO_;pm_math_VnD_HGFKVUOWdDAWvZhyEb[3]=(
pm_math_FqUCZrSGGNOuePgRr82o_8[3]-pm_math_FqUCZrSGGNOuePgRr82o_8[1])*
pm_math_kEbBObcYFIxUZ5_77V3CO_;break;case 1:pm_math_VnD_HGFKVUOWdDAWvZhyEb[1]=
sqrt(pm_math_VQ_T199jNSCPYuRLcxtO6j);pm_math_kEbBObcYFIxUZ5_77V3CO_=0.25/
pm_math_VnD_HGFKVUOWdDAWvZhyEb[1];pm_math_VnD_HGFKVUOWdDAWvZhyEb[0]=(
pm_math_FqUCZrSGGNOuePgRr82o_8[7]-pm_math_FqUCZrSGGNOuePgRr82o_8[5])*
pm_math_kEbBObcYFIxUZ5_77V3CO_;pm_math_VnD_HGFKVUOWdDAWvZhyEb[2]=(
pm_math_FqUCZrSGGNOuePgRr82o_8[1]+pm_math_FqUCZrSGGNOuePgRr82o_8[3])*
pm_math_kEbBObcYFIxUZ5_77V3CO_;pm_math_VnD_HGFKVUOWdDAWvZhyEb[3]=(
pm_math_FqUCZrSGGNOuePgRr82o_8[2]+pm_math_FqUCZrSGGNOuePgRr82o_8[6])*
pm_math_kEbBObcYFIxUZ5_77V3CO_;break;case 2:pm_math_VnD_HGFKVUOWdDAWvZhyEb[2]=
sqrt(pm_math_VkE6m1GxptxxYDqNBGUEZw);pm_math_kEbBObcYFIxUZ5_77V3CO_=0.25/
pm_math_VnD_HGFKVUOWdDAWvZhyEb[2];pm_math_VnD_HGFKVUOWdDAWvZhyEb[0]=(
pm_math_FqUCZrSGGNOuePgRr82o_8[2]-pm_math_FqUCZrSGGNOuePgRr82o_8[6])*
pm_math_kEbBObcYFIxUZ5_77V3CO_;pm_math_VnD_HGFKVUOWdDAWvZhyEb[3]=(
pm_math_FqUCZrSGGNOuePgRr82o_8[5]+pm_math_FqUCZrSGGNOuePgRr82o_8[7])*
pm_math_kEbBObcYFIxUZ5_77V3CO_;pm_math_VnD_HGFKVUOWdDAWvZhyEb[1]=(
pm_math_FqUCZrSGGNOuePgRr82o_8[3]+pm_math_FqUCZrSGGNOuePgRr82o_8[1])*
pm_math_kEbBObcYFIxUZ5_77V3CO_;break;case 3:pm_math_VnD_HGFKVUOWdDAWvZhyEb[3]=
sqrt(pm_math_VcAmFCMQr_KodumW__UKW6);pm_math_kEbBObcYFIxUZ5_77V3CO_=0.25/
pm_math_VnD_HGFKVUOWdDAWvZhyEb[3];pm_math_VnD_HGFKVUOWdDAWvZhyEb[0]=(
pm_math_FqUCZrSGGNOuePgRr82o_8[3]-pm_math_FqUCZrSGGNOuePgRr82o_8[1])*
pm_math_kEbBObcYFIxUZ5_77V3CO_;pm_math_VnD_HGFKVUOWdDAWvZhyEb[1]=(
pm_math_FqUCZrSGGNOuePgRr82o_8[6]+pm_math_FqUCZrSGGNOuePgRr82o_8[2])*
pm_math_kEbBObcYFIxUZ5_77V3CO_;pm_math_VnD_HGFKVUOWdDAWvZhyEb[2]=(
pm_math_FqUCZrSGGNOuePgRr82o_8[7]+pm_math_FqUCZrSGGNOuePgRr82o_8[5])*
pm_math_kEbBObcYFIxUZ5_77V3CO_;break;}pm_math_kEbBObcYFIxUZ5_77V3CO_=1.0/sqrt(
pm_math_VnD_HGFKVUOWdDAWvZhyEb[0]*pm_math_VnD_HGFKVUOWdDAWvZhyEb[0]+
pm_math_VnD_HGFKVUOWdDAWvZhyEb[1]*pm_math_VnD_HGFKVUOWdDAWvZhyEb[1]+
pm_math_VnD_HGFKVUOWdDAWvZhyEb[2]*pm_math_VnD_HGFKVUOWdDAWvZhyEb[2]+
pm_math_VnD_HGFKVUOWdDAWvZhyEb[3]*pm_math_VnD_HGFKVUOWdDAWvZhyEb[3]);
pm_math_VnD_HGFKVUOWdDAWvZhyEb[0]*=pm_math_kEbBObcYFIxUZ5_77V3CO_;
pm_math_VnD_HGFKVUOWdDAWvZhyEb[1]*=pm_math_kEbBObcYFIxUZ5_77V3CO_;
pm_math_VnD_HGFKVUOWdDAWvZhyEb[2]*=pm_math_kEbBObcYFIxUZ5_77V3CO_;
pm_math_VnD_HGFKVUOWdDAWvZhyEb[3]*=pm_math_kEbBObcYFIxUZ5_77V3CO_;
}

template <class T> void pm_math_Matrix3x3_compose_ra(const T*pm_math_Fcuud3IN0odMZi54a1R_8f,const
T*pm_math__09m2ugY6U_OXH9Il_a7Bj,T*pm_math__1Zf2IciMRCub1vvbEr1C4)
{
  (void)0;;*pm_math__1Zf2IciMRCub1vvbEr1C4++=pm_math_Fcuud3IN0odMZi54a1R_8f[0]*
pm_math__09m2ugY6U_OXH9Il_a7Bj[0]+pm_math_Fcuud3IN0odMZi54a1R_8f[1]*
pm_math__09m2ugY6U_OXH9Il_a7Bj[3]+pm_math_Fcuud3IN0odMZi54a1R_8f[2]*
pm_math__09m2ugY6U_OXH9Il_a7Bj[6];*pm_math__1Zf2IciMRCub1vvbEr1C4++=
pm_math_Fcuud3IN0odMZi54a1R_8f[0]*pm_math__09m2ugY6U_OXH9Il_a7Bj[1]+
pm_math_Fcuud3IN0odMZi54a1R_8f[1]*pm_math__09m2ugY6U_OXH9Il_a7Bj[4]+
pm_math_Fcuud3IN0odMZi54a1R_8f[2]*pm_math__09m2ugY6U_OXH9Il_a7Bj[7];*
pm_math__1Zf2IciMRCub1vvbEr1C4++=pm_math_Fcuud3IN0odMZi54a1R_8f[0]*
pm_math__09m2ugY6U_OXH9Il_a7Bj[2]+pm_math_Fcuud3IN0odMZi54a1R_8f[1]*
pm_math__09m2ugY6U_OXH9Il_a7Bj[5]+pm_math_Fcuud3IN0odMZi54a1R_8f[2]*
pm_math__09m2ugY6U_OXH9Il_a7Bj[8];*pm_math__1Zf2IciMRCub1vvbEr1C4++=
pm_math_Fcuud3IN0odMZi54a1R_8f[3]*pm_math__09m2ugY6U_OXH9Il_a7Bj[0]+
pm_math_Fcuud3IN0odMZi54a1R_8f[4]*pm_math__09m2ugY6U_OXH9Il_a7Bj[3]+
pm_math_Fcuud3IN0odMZi54a1R_8f[5]*pm_math__09m2ugY6U_OXH9Il_a7Bj[6];*
pm_math__1Zf2IciMRCub1vvbEr1C4++=pm_math_Fcuud3IN0odMZi54a1R_8f[3]*
pm_math__09m2ugY6U_OXH9Il_a7Bj[1]+pm_math_Fcuud3IN0odMZi54a1R_8f[4]*
pm_math__09m2ugY6U_OXH9Il_a7Bj[4]+pm_math_Fcuud3IN0odMZi54a1R_8f[5]*
pm_math__09m2ugY6U_OXH9Il_a7Bj[7];*pm_math__1Zf2IciMRCub1vvbEr1C4++=
pm_math_Fcuud3IN0odMZi54a1R_8f[3]*pm_math__09m2ugY6U_OXH9Il_a7Bj[2]+
pm_math_Fcuud3IN0odMZi54a1R_8f[4]*pm_math__09m2ugY6U_OXH9Il_a7Bj[5]+
pm_math_Fcuud3IN0odMZi54a1R_8f[5]*pm_math__09m2ugY6U_OXH9Il_a7Bj[8];*
pm_math__1Zf2IciMRCub1vvbEr1C4++=pm_math_Fcuud3IN0odMZi54a1R_8f[6]*
pm_math__09m2ugY6U_OXH9Il_a7Bj[0]+pm_math_Fcuud3IN0odMZi54a1R_8f[7]*
pm_math__09m2ugY6U_OXH9Il_a7Bj[3]+pm_math_Fcuud3IN0odMZi54a1R_8f[8]*
pm_math__09m2ugY6U_OXH9Il_a7Bj[6];*pm_math__1Zf2IciMRCub1vvbEr1C4++=
pm_math_Fcuud3IN0odMZi54a1R_8f[6]*pm_math__09m2ugY6U_OXH9Il_a7Bj[1]+
pm_math_Fcuud3IN0odMZi54a1R_8f[7]*pm_math__09m2ugY6U_OXH9Il_a7Bj[4]+
pm_math_Fcuud3IN0odMZi54a1R_8f[8]*pm_math__09m2ugY6U_OXH9Il_a7Bj[7];*
pm_math__1Zf2IciMRCub1vvbEr1C4=pm_math_Fcuud3IN0odMZi54a1R_8f[6]*
pm_math__09m2ugY6U_OXH9Il_a7Bj[2]+pm_math_Fcuud3IN0odMZi54a1R_8f[7]*
pm_math__09m2ugY6U_OXH9Il_a7Bj[5]+pm_math_Fcuud3IN0odMZi54a1R_8f[8]*
pm_math__09m2ugY6U_OXH9Il_a7Bj[8];
}

template <class T> void pm_math_Matrix3x3_composeTranspose_ra(
const T*pm_math_Fcuud3IN0odMZi54a1R_8f,const T*
pm_math__09m2ugY6U_OXH9Il_a7Bj,T*pm_math__1Zf2IciMRCub1vvbEr1C4)
{
  (void)0;;*pm_math__1Zf2IciMRCub1vvbEr1C4++=pm_math_Fcuud3IN0odMZi54a1R_8f[0]*
pm_math__09m2ugY6U_OXH9Il_a7Bj[0]+pm_math_Fcuud3IN0odMZi54a1R_8f[1]*
pm_math__09m2ugY6U_OXH9Il_a7Bj[1]+pm_math_Fcuud3IN0odMZi54a1R_8f[2]*
pm_math__09m2ugY6U_OXH9Il_a7Bj[2];*pm_math__1Zf2IciMRCub1vvbEr1C4++=
pm_math_Fcuud3IN0odMZi54a1R_8f[0]*pm_math__09m2ugY6U_OXH9Il_a7Bj[3]+
pm_math_Fcuud3IN0odMZi54a1R_8f[1]*pm_math__09m2ugY6U_OXH9Il_a7Bj[4]+
pm_math_Fcuud3IN0odMZi54a1R_8f[2]*pm_math__09m2ugY6U_OXH9Il_a7Bj[5];*
pm_math__1Zf2IciMRCub1vvbEr1C4++=pm_math_Fcuud3IN0odMZi54a1R_8f[0]*
pm_math__09m2ugY6U_OXH9Il_a7Bj[6]+pm_math_Fcuud3IN0odMZi54a1R_8f[1]*
pm_math__09m2ugY6U_OXH9Il_a7Bj[7]+pm_math_Fcuud3IN0odMZi54a1R_8f[2]*
pm_math__09m2ugY6U_OXH9Il_a7Bj[8];*pm_math__1Zf2IciMRCub1vvbEr1C4++=
pm_math_Fcuud3IN0odMZi54a1R_8f[3]*pm_math__09m2ugY6U_OXH9Il_a7Bj[0]+
pm_math_Fcuud3IN0odMZi54a1R_8f[4]*pm_math__09m2ugY6U_OXH9Il_a7Bj[1]+
pm_math_Fcuud3IN0odMZi54a1R_8f[5]*pm_math__09m2ugY6U_OXH9Il_a7Bj[2];*
pm_math__1Zf2IciMRCub1vvbEr1C4++=pm_math_Fcuud3IN0odMZi54a1R_8f[3]*
pm_math__09m2ugY6U_OXH9Il_a7Bj[3]+pm_math_Fcuud3IN0odMZi54a1R_8f[4]*
pm_math__09m2ugY6U_OXH9Il_a7Bj[4]+pm_math_Fcuud3IN0odMZi54a1R_8f[5]*
pm_math__09m2ugY6U_OXH9Il_a7Bj[5];*pm_math__1Zf2IciMRCub1vvbEr1C4++=
pm_math_Fcuud3IN0odMZi54a1R_8f[3]*pm_math__09m2ugY6U_OXH9Il_a7Bj[6]+
pm_math_Fcuud3IN0odMZi54a1R_8f[4]*pm_math__09m2ugY6U_OXH9Il_a7Bj[7]+
pm_math_Fcuud3IN0odMZi54a1R_8f[5]*pm_math__09m2ugY6U_OXH9Il_a7Bj[8];*
pm_math__1Zf2IciMRCub1vvbEr1C4++=pm_math_Fcuud3IN0odMZi54a1R_8f[6]*
pm_math__09m2ugY6U_OXH9Il_a7Bj[0]+pm_math_Fcuud3IN0odMZi54a1R_8f[7]*
pm_math__09m2ugY6U_OXH9Il_a7Bj[1]+pm_math_Fcuud3IN0odMZi54a1R_8f[8]*
pm_math__09m2ugY6U_OXH9Il_a7Bj[2];*pm_math__1Zf2IciMRCub1vvbEr1C4++=
pm_math_Fcuud3IN0odMZi54a1R_8f[6]*pm_math__09m2ugY6U_OXH9Il_a7Bj[3]+
pm_math_Fcuud3IN0odMZi54a1R_8f[7]*pm_math__09m2ugY6U_OXH9Il_a7Bj[4]+
pm_math_Fcuud3IN0odMZi54a1R_8f[8]*pm_math__09m2ugY6U_OXH9Il_a7Bj[5];*
pm_math__1Zf2IciMRCub1vvbEr1C4=pm_math_Fcuud3IN0odMZi54a1R_8f[6]*
pm_math__09m2ugY6U_OXH9Il_a7Bj[6]+pm_math_Fcuud3IN0odMZi54a1R_8f[7]*
pm_math__09m2ugY6U_OXH9Il_a7Bj[7]+pm_math_Fcuud3IN0odMZi54a1R_8f[8]*
pm_math__09m2ugY6U_OXH9Il_a7Bj[8];
}

template <class T> void pm_math_Matrix3x3_transposeCompose_ra(const T*pm_math_Fcuud3IN0odMZi54a1R_8f,const T*
pm_math__09m2ugY6U_OXH9Il_a7Bj,T*pm_math__1Zf2IciMRCub1vvbEr1C4)
{
  (void)0;;*pm_math__1Zf2IciMRCub1vvbEr1C4++=pm_math_Fcuud3IN0odMZi54a1R_8f[0]*
pm_math__09m2ugY6U_OXH9Il_a7Bj[0]+pm_math_Fcuud3IN0odMZi54a1R_8f[3]*
pm_math__09m2ugY6U_OXH9Il_a7Bj[3]+pm_math_Fcuud3IN0odMZi54a1R_8f[6]*
pm_math__09m2ugY6U_OXH9Il_a7Bj[6];*pm_math__1Zf2IciMRCub1vvbEr1C4++=
pm_math_Fcuud3IN0odMZi54a1R_8f[0]*pm_math__09m2ugY6U_OXH9Il_a7Bj[1]+
pm_math_Fcuud3IN0odMZi54a1R_8f[3]*pm_math__09m2ugY6U_OXH9Il_a7Bj[4]+
pm_math_Fcuud3IN0odMZi54a1R_8f[6]*pm_math__09m2ugY6U_OXH9Il_a7Bj[7];*
pm_math__1Zf2IciMRCub1vvbEr1C4++=pm_math_Fcuud3IN0odMZi54a1R_8f[0]*
pm_math__09m2ugY6U_OXH9Il_a7Bj[2]+pm_math_Fcuud3IN0odMZi54a1R_8f[3]*
pm_math__09m2ugY6U_OXH9Il_a7Bj[5]+pm_math_Fcuud3IN0odMZi54a1R_8f[6]*
pm_math__09m2ugY6U_OXH9Il_a7Bj[8];*pm_math__1Zf2IciMRCub1vvbEr1C4++=
pm_math_Fcuud3IN0odMZi54a1R_8f[1]*pm_math__09m2ugY6U_OXH9Il_a7Bj[0]+
pm_math_Fcuud3IN0odMZi54a1R_8f[4]*pm_math__09m2ugY6U_OXH9Il_a7Bj[3]+
pm_math_Fcuud3IN0odMZi54a1R_8f[7]*pm_math__09m2ugY6U_OXH9Il_a7Bj[6];*
pm_math__1Zf2IciMRCub1vvbEr1C4++=pm_math_Fcuud3IN0odMZi54a1R_8f[1]*
pm_math__09m2ugY6U_OXH9Il_a7Bj[1]+pm_math_Fcuud3IN0odMZi54a1R_8f[4]*
pm_math__09m2ugY6U_OXH9Il_a7Bj[4]+pm_math_Fcuud3IN0odMZi54a1R_8f[7]*
pm_math__09m2ugY6U_OXH9Il_a7Bj[7];*pm_math__1Zf2IciMRCub1vvbEr1C4++=
pm_math_Fcuud3IN0odMZi54a1R_8f[1]*pm_math__09m2ugY6U_OXH9Il_a7Bj[2]+
pm_math_Fcuud3IN0odMZi54a1R_8f[4]*pm_math__09m2ugY6U_OXH9Il_a7Bj[5]+
pm_math_Fcuud3IN0odMZi54a1R_8f[7]*pm_math__09m2ugY6U_OXH9Il_a7Bj[8];*
pm_math__1Zf2IciMRCub1vvbEr1C4++=pm_math_Fcuud3IN0odMZi54a1R_8f[2]*
pm_math__09m2ugY6U_OXH9Il_a7Bj[0]+pm_math_Fcuud3IN0odMZi54a1R_8f[5]*
pm_math__09m2ugY6U_OXH9Il_a7Bj[3]+pm_math_Fcuud3IN0odMZi54a1R_8f[8]*
pm_math__09m2ugY6U_OXH9Il_a7Bj[6];*pm_math__1Zf2IciMRCub1vvbEr1C4++=
pm_math_Fcuud3IN0odMZi54a1R_8f[2]*pm_math__09m2ugY6U_OXH9Il_a7Bj[1]+
pm_math_Fcuud3IN0odMZi54a1R_8f[5]*pm_math__09m2ugY6U_OXH9Il_a7Bj[4]+
pm_math_Fcuud3IN0odMZi54a1R_8f[8]*pm_math__09m2ugY6U_OXH9Il_a7Bj[7];*
pm_math__1Zf2IciMRCub1vvbEr1C4=pm_math_Fcuud3IN0odMZi54a1R_8f[2]*
pm_math__09m2ugY6U_OXH9Il_a7Bj[2]+pm_math_Fcuud3IN0odMZi54a1R_8f[5]*
pm_math__09m2ugY6U_OXH9Il_a7Bj[5]+pm_math_Fcuud3IN0odMZi54a1R_8f[8]*
pm_math__09m2ugY6U_OXH9Il_a7Bj[8];
}

template <class T> void pm_math_Matrix3x3_preCross_ra(const
T*pm_math_F2l4p_g4sn02huHNflQjMH,const T*
pm_math_VgJW5ZqpwPpuY1inYtaofQ,T*pm_math__1Zf2IciMRCub1vvbEr1C4)
{
  (void)0;;*pm_math__1Zf2IciMRCub1vvbEr1C4++=pm_math_F2l4p_g4sn02huHNflQjMH[6]*
pm_math_VgJW5ZqpwPpuY1inYtaofQ[1]-pm_math_F2l4p_g4sn02huHNflQjMH[3]*
pm_math_VgJW5ZqpwPpuY1inYtaofQ[2];*pm_math__1Zf2IciMRCub1vvbEr1C4++=
pm_math_F2l4p_g4sn02huHNflQjMH[7]*pm_math_VgJW5ZqpwPpuY1inYtaofQ[1]-
pm_math_F2l4p_g4sn02huHNflQjMH[4]*pm_math_VgJW5ZqpwPpuY1inYtaofQ[2];*
pm_math__1Zf2IciMRCub1vvbEr1C4++=pm_math_F2l4p_g4sn02huHNflQjMH[8]*
pm_math_VgJW5ZqpwPpuY1inYtaofQ[1]-pm_math_F2l4p_g4sn02huHNflQjMH[5]*
pm_math_VgJW5ZqpwPpuY1inYtaofQ[2];*pm_math__1Zf2IciMRCub1vvbEr1C4++=
pm_math_F2l4p_g4sn02huHNflQjMH[0]*pm_math_VgJW5ZqpwPpuY1inYtaofQ[2]-
pm_math_F2l4p_g4sn02huHNflQjMH[6]*pm_math_VgJW5ZqpwPpuY1inYtaofQ[0];*
pm_math__1Zf2IciMRCub1vvbEr1C4++=pm_math_F2l4p_g4sn02huHNflQjMH[1]*
pm_math_VgJW5ZqpwPpuY1inYtaofQ[2]-pm_math_F2l4p_g4sn02huHNflQjMH[7]*
pm_math_VgJW5ZqpwPpuY1inYtaofQ[0];*pm_math__1Zf2IciMRCub1vvbEr1C4++=
pm_math_F2l4p_g4sn02huHNflQjMH[2]*pm_math_VgJW5ZqpwPpuY1inYtaofQ[2]-
pm_math_F2l4p_g4sn02huHNflQjMH[8]*pm_math_VgJW5ZqpwPpuY1inYtaofQ[0];*
pm_math__1Zf2IciMRCub1vvbEr1C4++=pm_math_F2l4p_g4sn02huHNflQjMH[3]*
pm_math_VgJW5ZqpwPpuY1inYtaofQ[0]-pm_math_F2l4p_g4sn02huHNflQjMH[0]*
pm_math_VgJW5ZqpwPpuY1inYtaofQ[1];*pm_math__1Zf2IciMRCub1vvbEr1C4++=
pm_math_F2l4p_g4sn02huHNflQjMH[4]*pm_math_VgJW5ZqpwPpuY1inYtaofQ[0]-
pm_math_F2l4p_g4sn02huHNflQjMH[1]*pm_math_VgJW5ZqpwPpuY1inYtaofQ[1];*
pm_math__1Zf2IciMRCub1vvbEr1C4=pm_math_F2l4p_g4sn02huHNflQjMH[5]*
pm_math_VgJW5ZqpwPpuY1inYtaofQ[0]-pm_math_F2l4p_g4sn02huHNflQjMH[2]*
pm_math_VgJW5ZqpwPpuY1inYtaofQ[1];
}

template <class T> void pm_math_Matrix3x3_postCross_ra(const
T*pm_math_F2l4p_g4sn02huHNflQjMH,const T*
pm_math_VgJW5ZqpwPpuY1inYtaofQ,T*pm_math__1Zf2IciMRCub1vvbEr1C4)
{
  (void)0;;*pm_math__1Zf2IciMRCub1vvbEr1C4++=pm_math_F2l4p_g4sn02huHNflQjMH[1]*
pm_math_VgJW5ZqpwPpuY1inYtaofQ[2]-pm_math_F2l4p_g4sn02huHNflQjMH[2]*
pm_math_VgJW5ZqpwPpuY1inYtaofQ[1];*pm_math__1Zf2IciMRCub1vvbEr1C4++=
pm_math_F2l4p_g4sn02huHNflQjMH[2]*pm_math_VgJW5ZqpwPpuY1inYtaofQ[0]-
pm_math_F2l4p_g4sn02huHNflQjMH[0]*pm_math_VgJW5ZqpwPpuY1inYtaofQ[2];*
pm_math__1Zf2IciMRCub1vvbEr1C4++=pm_math_F2l4p_g4sn02huHNflQjMH[0]*
pm_math_VgJW5ZqpwPpuY1inYtaofQ[1]-pm_math_F2l4p_g4sn02huHNflQjMH[1]*
pm_math_VgJW5ZqpwPpuY1inYtaofQ[0];*pm_math__1Zf2IciMRCub1vvbEr1C4++=
pm_math_F2l4p_g4sn02huHNflQjMH[4]*pm_math_VgJW5ZqpwPpuY1inYtaofQ[2]-
pm_math_F2l4p_g4sn02huHNflQjMH[5]*pm_math_VgJW5ZqpwPpuY1inYtaofQ[1];*
pm_math__1Zf2IciMRCub1vvbEr1C4++=pm_math_F2l4p_g4sn02huHNflQjMH[5]*
pm_math_VgJW5ZqpwPpuY1inYtaofQ[0]-pm_math_F2l4p_g4sn02huHNflQjMH[3]*
pm_math_VgJW5ZqpwPpuY1inYtaofQ[2];*pm_math__1Zf2IciMRCub1vvbEr1C4++=
pm_math_F2l4p_g4sn02huHNflQjMH[3]*pm_math_VgJW5ZqpwPpuY1inYtaofQ[1]-
pm_math_F2l4p_g4sn02huHNflQjMH[4]*pm_math_VgJW5ZqpwPpuY1inYtaofQ[0];*
pm_math__1Zf2IciMRCub1vvbEr1C4++=pm_math_F2l4p_g4sn02huHNflQjMH[7]*
pm_math_VgJW5ZqpwPpuY1inYtaofQ[2]-pm_math_F2l4p_g4sn02huHNflQjMH[8]*
pm_math_VgJW5ZqpwPpuY1inYtaofQ[1];*pm_math__1Zf2IciMRCub1vvbEr1C4++=
pm_math_F2l4p_g4sn02huHNflQjMH[8]*pm_math_VgJW5ZqpwPpuY1inYtaofQ[0]-
pm_math_F2l4p_g4sn02huHNflQjMH[6]*pm_math_VgJW5ZqpwPpuY1inYtaofQ[2];*
pm_math__1Zf2IciMRCub1vvbEr1C4=pm_math_F2l4p_g4sn02huHNflQjMH[6]*
pm_math_VgJW5ZqpwPpuY1inYtaofQ[1]-pm_math_F2l4p_g4sn02huHNflQjMH[7]*
pm_math_VgJW5ZqpwPpuY1inYtaofQ[0];
}

template <class T> void pm_math_Matrix3x3_xform_ra(const T
*pm_math_F2l4p_g4sn02huHNflQjMH,const T*pm_math_VgJW5ZqpwPpuY1inYtaofQ,
T*pm_math__1Zf2IciMRCub1vvbEr1C4)
{
  (void)0;;*pm_math__1Zf2IciMRCub1vvbEr1C4++=pm_math_F2l4p_g4sn02huHNflQjMH[0]*
pm_math_VgJW5ZqpwPpuY1inYtaofQ[0]+pm_math_F2l4p_g4sn02huHNflQjMH[1]*
pm_math_VgJW5ZqpwPpuY1inYtaofQ[1]+pm_math_F2l4p_g4sn02huHNflQjMH[2]*
pm_math_VgJW5ZqpwPpuY1inYtaofQ[2];*pm_math__1Zf2IciMRCub1vvbEr1C4++=
pm_math_F2l4p_g4sn02huHNflQjMH[3]*pm_math_VgJW5ZqpwPpuY1inYtaofQ[0]+
pm_math_F2l4p_g4sn02huHNflQjMH[4]*pm_math_VgJW5ZqpwPpuY1inYtaofQ[1]+
pm_math_F2l4p_g4sn02huHNflQjMH[5]*pm_math_VgJW5ZqpwPpuY1inYtaofQ[2];*
pm_math__1Zf2IciMRCub1vvbEr1C4++=pm_math_F2l4p_g4sn02huHNflQjMH[6]*
pm_math_VgJW5ZqpwPpuY1inYtaofQ[0]+pm_math_F2l4p_g4sn02huHNflQjMH[7]*
pm_math_VgJW5ZqpwPpuY1inYtaofQ[1]+pm_math_F2l4p_g4sn02huHNflQjMH[8]*
pm_math_VgJW5ZqpwPpuY1inYtaofQ[2];
}

template <class T> void pm_math_Matrix3x3_transposeXform_ra(
const T*pm_math_F2l4p_g4sn02huHNflQjMH,const T*
pm_math_VgJW5ZqpwPpuY1inYtaofQ,T*pm_math__1Zf2IciMRCub1vvbEr1C4)
{
  (void)0;
;*pm_math__1Zf2IciMRCub1vvbEr1C4++=pm_math_F2l4p_g4sn02huHNflQjMH[0]*
pm_math_VgJW5ZqpwPpuY1inYtaofQ[0]+pm_math_F2l4p_g4sn02huHNflQjMH[3]*
pm_math_VgJW5ZqpwPpuY1inYtaofQ[1]+pm_math_F2l4p_g4sn02huHNflQjMH[6]*
pm_math_VgJW5ZqpwPpuY1inYtaofQ[2];*pm_math__1Zf2IciMRCub1vvbEr1C4++=
pm_math_F2l4p_g4sn02huHNflQjMH[1]*pm_math_VgJW5ZqpwPpuY1inYtaofQ[0]+
pm_math_F2l4p_g4sn02huHNflQjMH[4]*pm_math_VgJW5ZqpwPpuY1inYtaofQ[1]+
pm_math_F2l4p_g4sn02huHNflQjMH[7]*pm_math_VgJW5ZqpwPpuY1inYtaofQ[2];*
pm_math__1Zf2IciMRCub1vvbEr1C4=pm_math_F2l4p_g4sn02huHNflQjMH[2]*
pm_math_VgJW5ZqpwPpuY1inYtaofQ[0]+pm_math_F2l4p_g4sn02huHNflQjMH[5]*
pm_math_VgJW5ZqpwPpuY1inYtaofQ[1]+pm_math_F2l4p_g4sn02huHNflQjMH[8]*
pm_math_VgJW5ZqpwPpuY1inYtaofQ[2];
}

template <class T> static T pm_math_V_FoLVGO37G4eXJuucTSjm
(const T*pm_math_F2l4p_g4sn02huHNflQjMH)
{
  const T*a=pm_math_F2l4p_g4sn02huHNflQjMH;const T*b=pm_math_F2l4p_g4sn02huHNflQjMH+3
;const T*pm_math_FFZbGh27ya8eem_J_hUtAZ=pm_math_F2l4p_g4sn02huHNflQjMH+6;
T pm_math__8EEqwe2wQpGa1Zrsu6MRI[3];pm_math_Vector3_cross_ra(b,
pm_math_FFZbGh27ya8eem_J_hUtAZ,pm_math__8EEqwe2wQpGa1Zrsu6MRI);return
pm_math_Vector3_dot_ra(a,pm_math__8EEqwe2wQpGa1Zrsu6MRI);
}

template <class T> void pm_math_Matrix3x3_minRotation_ra(const T*pm_math_VnD_HGFKVUOWdDAWvZhyEb,
T*pm_math_FqUCZrSGGNOuePgRr82o_8,int pm_math_kbzF46WM0FtKWeE0s7WR95[3])
{
const boolean_T pm_math_kA1JM95JsV0RcLxI4SAZNK=pm_math_V_FoLVGO37G4eXJuucTSjm(
pm_math_VnD_HGFKVUOWdDAWvZhyEb)<0.0;const int pm_math_V5sGO9ANMS0VVLkblQ_gRK[6
][3]={{0,1,2},{1,2,0},{2,0,1},{0,2,1},{2,1,0},{1,0,2}};const int*
pm_math_VdLWATHg4Lt5gihTpf28WR=NULL;int pm_math_VR3Gzdw6XzCOZ1z6An_brA=0;
T pm_math_FoqvJs6gbxKlbP5pszn7ig= -4.0;int pm_math_kwrB3ZoKf7OufTHWaHJV7a
;int pm_math_krYzdgglpzS1gus79MQZpM;for(pm_math_krYzdgglpzS1gus79MQZpM=0;
pm_math_krYzdgglpzS1gus79MQZpM<6;++pm_math_krYzdgglpzS1gus79MQZpM){const int*
pm__lqjegyKuwStj56WZLiC_e=pm_math_V5sGO9ANMS0VVLkblQ_gRK[
pm_math_krYzdgglpzS1gus79MQZpM];const boolean_T pm_math_VDv36i3RLWK7VqXFxwtULD
=pm_math_krYzdgglpzS1gus79MQZpM>=3;int pm_math__Ff9JZhSNPK7VeK_UaPbQQ;for(
pm_math__Ff9JZhSNPK7VeK_UaPbQQ=0;pm_math__Ff9JZhSNPK7VeK_UaPbQQ<8;++
pm_math__Ff9JZhSNPK7VeK_UaPbQQ){const boolean_T pm_math_ktbvjqcPhWdOj9Vukr1STa
=((pm_math__Ff9JZhSNPK7VeK_UaPbQQ+(pm_math__Ff9JZhSNPK7VeK_UaPbQQ>>1)+(
pm_math__Ff9JZhSNPK7VeK_UaPbQQ>>2))&1)==1;T pm_math_Vgyn9_XGvLSUaauCGjY4RB
=0.0;if((pm_math_VDv36i3RLWK7VqXFxwtULD^pm_math_ktbvjqcPhWdOj9Vukr1STa)!=
pm_math_kA1JM95JsV0RcLxI4SAZNK)continue;for(pm_math_kwrB3ZoKf7OufTHWaHJV7a=0;
pm_math_kwrB3ZoKf7OufTHWaHJV7a<3;++pm_math_kwrB3ZoKf7OufTHWaHJV7a){const T
pm_math_F32Ql82vv6pW_PYIdpkFQ0=pm_math_VnD_HGFKVUOWdDAWvZhyEb[
pm__lqjegyKuwStj56WZLiC_e[pm_math_kwrB3ZoKf7OufTHWaHJV7a]+3*
pm_math_kwrB3ZoKf7OufTHWaHJV7a];pm_math_Vgyn9_XGvLSUaauCGjY4RB+=(
pm_math__Ff9JZhSNPK7VeK_UaPbQQ&(1<<pm_math_kwrB3ZoKf7OufTHWaHJV7a))==0?+
pm_math_F32Ql82vv6pW_PYIdpkFQ0:-pm_math_F32Ql82vv6pW_PYIdpkFQ0;}if(
pm_math_Vgyn9_XGvLSUaauCGjY4RB>pm_math_FoqvJs6gbxKlbP5pszn7ig){
pm_math_VdLWATHg4Lt5gihTpf28WR=pm__lqjegyKuwStj56WZLiC_e;
pm_math_VR3Gzdw6XzCOZ1z6An_brA=pm_math__Ff9JZhSNPK7VeK_UaPbQQ;
pm_math_FoqvJs6gbxKlbP5pszn7ig=pm_math_Vgyn9_XGvLSUaauCGjY4RB;}}}memcpy(
pm_math_kbzF46WM0FtKWeE0s7WR95,pm_math_VdLWATHg4Lt5gihTpf28WR,3*sizeof(int));{
T pm_math_Fk2O4u6vQUpibmbv8Kjgnn[9];if(pm_math_FqUCZrSGGNOuePgRr82o_8==
pm_math_VnD_HGFKVUOWdDAWvZhyEb){memcpy(pm_math_Fk2O4u6vQUpibmbv8Kjgnn,
pm_math_VnD_HGFKVUOWdDAWvZhyEb,9*sizeof(T));
pm_math_VnD_HGFKVUOWdDAWvZhyEb=pm_math_Fk2O4u6vQUpibmbv8Kjgnn;}
pm_math_FqUCZrSGGNOuePgRr82o_8[0]=pm_math_VnD_HGFKVUOWdDAWvZhyEb[
pm_math_kbzF46WM0FtKWeE0s7WR95[0]+0],pm_math_FqUCZrSGGNOuePgRr82o_8[3]=
pm_math_VnD_HGFKVUOWdDAWvZhyEb[pm_math_kbzF46WM0FtKWeE0s7WR95[0]+3],
pm_math_FqUCZrSGGNOuePgRr82o_8[6]=pm_math_VnD_HGFKVUOWdDAWvZhyEb[
pm_math_kbzF46WM0FtKWeE0s7WR95[0]+6];pm_math_FqUCZrSGGNOuePgRr82o_8[1]=
pm_math_VnD_HGFKVUOWdDAWvZhyEb[pm_math_kbzF46WM0FtKWeE0s7WR95[1]+0],
pm_math_FqUCZrSGGNOuePgRr82o_8[4]=pm_math_VnD_HGFKVUOWdDAWvZhyEb[
pm_math_kbzF46WM0FtKWeE0s7WR95[1]+3],pm_math_FqUCZrSGGNOuePgRr82o_8[7]=
pm_math_VnD_HGFKVUOWdDAWvZhyEb[pm_math_kbzF46WM0FtKWeE0s7WR95[1]+6];
pm_math_FqUCZrSGGNOuePgRr82o_8[2]=pm_math_VnD_HGFKVUOWdDAWvZhyEb[
pm_math_kbzF46WM0FtKWeE0s7WR95[2]+0],pm_math_FqUCZrSGGNOuePgRr82o_8[5]=
pm_math_VnD_HGFKVUOWdDAWvZhyEb[pm_math_kbzF46WM0FtKWeE0s7WR95[2]+3],
pm_math_FqUCZrSGGNOuePgRr82o_8[8]=pm_math_VnD_HGFKVUOWdDAWvZhyEb[
pm_math_kbzF46WM0FtKWeE0s7WR95[2]+6];}for(pm_math_kwrB3ZoKf7OufTHWaHJV7a=0;
pm_math_kwrB3ZoKf7OufTHWaHJV7a<3;++pm_math_kwrB3ZoKf7OufTHWaHJV7a)if((
pm_math_VR3Gzdw6XzCOZ1z6An_brA&(1<<pm_math_kwrB3ZoKf7OufTHWaHJV7a))!=0){T
*pm_math_FL1llpmubk_sXeD10Sbe3f=pm_math_FqUCZrSGGNOuePgRr82o_8+
pm_math_kwrB3ZoKf7OufTHWaHJV7a;pm_math_FL1llpmubk_sXeD10Sbe3f[0]= -
pm_math_FL1llpmubk_sXeD10Sbe3f[0];pm_math_FL1llpmubk_sXeD10Sbe3f[3]= -
pm_math_FL1llpmubk_sXeD10Sbe3f[3];pm_math_FL1llpmubk_sXeD10Sbe3f[6]= -
pm_math_FL1llpmubk_sXeD10Sbe3f[6];}
}

template <class T> void
pm_math_VWv_QvuNRm4viyJPPMYOkS(const T*pm_math__ut5UfJwzNlZ_XZC_yEgKo,
const T*b,uint32_T n,boolean_T pm_math_FEqGkTLIgMdKjaJ1lIV70b,T*x){
uint32_T pm_math_kwrB3ZoKf7OufTHWaHJV7a,pm_math_kyp6uAyJE40UVuAQNEYzS1;for(
pm_math_kwrB3ZoKf7OufTHWaHJV7a=0;pm_math_kwrB3ZoKf7OufTHWaHJV7a<n;x-=
pm_math_kwrB3ZoKf7OufTHWaHJV7a,pm_math_kwrB3ZoKf7OufTHWaHJV7a++){T
pm_math_FQferGZUKft3_i5GvYy4Oy=0.0;const T*pm_math_FL1llpmubk_sXeD10Sbe3f
=pm_math__ut5UfJwzNlZ_XZC_yEgKo++;for(pm_math_kyp6uAyJE40UVuAQNEYzS1=0;
pm_math_kyp6uAyJE40UVuAQNEYzS1<pm_math_kwrB3ZoKf7OufTHWaHJV7a;++
pm_math_kyp6uAyJE40UVuAQNEYzS1,pm_math_FL1llpmubk_sXeD10Sbe3f+=n)
pm_math_FQferGZUKft3_i5GvYy4Oy+= *pm_math_FL1llpmubk_sXeD10Sbe3f**x++;*x=
pm_math_FEqGkTLIgMdKjaJ1lIV70b?(*b++-pm_math_FQferGZUKft3_i5GvYy4Oy):((*b++-
pm_math_FQferGZUKft3_i5GvYy4Oy)/ *pm_math_FL1llpmubk_sXeD10Sbe3f);}}

template <class T> void
pm_math_VRlre4EWbUOCcXxk9GINjY(const T*pm_math_Vi4Cp0qK964NYTFMGr9Ttn,
const T*b,uint32_T n,boolean_T pm_math_FEqGkTLIgMdKjaJ1lIV70b,T*x){
uint32_T pm_math_kwrB3ZoKf7OufTHWaHJV7a,pm_math_kyp6uAyJE40UVuAQNEYzS1;
pm_math_Vi4Cp0qK964NYTFMGr9Ttn+=n*n;b+=n-1;x+=n-1;for(
pm_math_kwrB3ZoKf7OufTHWaHJV7a=0;pm_math_kwrB3ZoKf7OufTHWaHJV7a<n;x+=
pm_math_kwrB3ZoKf7OufTHWaHJV7a,pm_math_kwrB3ZoKf7OufTHWaHJV7a++){T
pm_math_FQferGZUKft3_i5GvYy4Oy=0.0;const T*pm_math_FL1llpmubk_sXeD10Sbe3f
= --pm_math_Vi4Cp0qK964NYTFMGr9Ttn;for(pm_math_kyp6uAyJE40UVuAQNEYzS1=
pm_math_kwrB3ZoKf7OufTHWaHJV7a;pm_math_kyp6uAyJE40UVuAQNEYzS1>0;
pm_math_FL1llpmubk_sXeD10Sbe3f-=n,--pm_math_kyp6uAyJE40UVuAQNEYzS1)
pm_math_FQferGZUKft3_i5GvYy4Oy+= *pm_math_FL1llpmubk_sXeD10Sbe3f**x--;*x=
pm_math_FEqGkTLIgMdKjaJ1lIV70b?(*b-- -pm_math_FQferGZUKft3_i5GvYy4Oy):((*b-- -
pm_math_FQferGZUKft3_i5GvYy4Oy)/ *pm_math_FL1llpmubk_sXeD10Sbe3f);}}

template <class T> void pm_math_lin_alg_choleskySolve(const T*pm_math_knZMtGN5npp9jax3au9uBf
,const T*b,uint32_T n,T*x,T*pm_math_VRZCD_UL_ESThy75dC9J8D){
pm_math_VWv_QvuNRm4viyJPPMYOkS(pm_math_knZMtGN5npp9jax3au9uBf,b,n,false,
pm_math_VRZCD_UL_ESThy75dC9J8D);pm_math_VRlre4EWbUOCcXxk9GINjY(
pm_math_knZMtGN5npp9jax3au9uBf,pm_math_VRZCD_UL_ESThy75dC9J8D,n,false,x);}

template <class T> void
solveSymmetricPosDef(const T*sm_F2l4p_g4sn02huHNflQjMH,const T*b,
uint32_T n,uint32_T sm_V2__YrimeI4E_yWnhKofpy,T*x,T*
sm_VRZCD_UL_ESThy75dC9J8D){if(sm_V2__YrimeI4E_yWnhKofpy==1){if(n<=1){if(n==1)*
x= *b/ *sm_F2l4p_g4sn02huHNflQjMH;}else pm_math_lin_alg_choleskySolve(
sm_F2l4p_g4sn02huHNflQjMH,b,n,x,sm_VRZCD_UL_ESThy75dC9J8D);}else{uint32_T
pm_Fr_bHKkQKFWbfi50VWd5Pw=0;if(n<=1){if(n==1)for(pm_Fr_bHKkQKFWbfi50VWd5Pw=0;
pm_Fr_bHKkQKFWbfi50VWd5Pw<sm_V2__YrimeI4E_yWnhKofpy;++
pm_Fr_bHKkQKFWbfi50VWd5Pw)*x++= *b++/ *sm_F2l4p_g4sn02huHNflQjMH;}else{for(
pm_Fr_bHKkQKFWbfi50VWd5Pw=0;pm_Fr_bHKkQKFWbfi50VWd5Pw<
sm_V2__YrimeI4E_yWnhKofpy;++pm_Fr_bHKkQKFWbfi50VWd5Pw)
pm_math_lin_alg_choleskySolve(sm_F2l4p_g4sn02huHNflQjMH,b+
pm_Fr_bHKkQKFWbfi50VWd5Pw*n,n,x+pm_Fr_bHKkQKFWbfi50VWd5Pw*n,
sm_VRZCD_UL_ESThy75dC9J8D);}}}


#endif //PSOPT_MCODER_MATH_H
