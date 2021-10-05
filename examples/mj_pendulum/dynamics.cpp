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
 * \date   7/5/21
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

void dynamics(const double* state, const double* input, double* deriv)
{
  double xx[14];
  xx[0] = 2.0;
  xx[1] = 0.125;
  xx[2] = 0.5;
  xx[3] = xx[2] * state[2];
  xx[4] = sin(xx[3]);
  xx[5] = cos(xx[3]);
  xx[3] = xx[0] * xx[5] * xx[4];
  xx[6] = xx[2] * xx[3];
  xx[7] = 0.03125000000000001;

  xx[8] = xx[1] / xx[7];
  xx[7] = 1.0;
  xx[9] = xx[0] * xx[4] * xx[4] - xx[7];
  xx[10] = (xx[2] - xx[1] * xx[8]) * xx[9];
  xx[2] = xx[0] + xx[6] * xx[3] + xx[10] * xx[9];
  memcpy(xx + 9, xx + 2, 1 * sizeof(double));

  xx[11] = 9.810000000000002;
  xx[12] = 1.77635683940025e-15;
  xx[13] = (input[0] + xx[0] * xx[1] * state[3] * state[3] * xx[4] * xx[5]) /
           xx[9] + xx[11] * (xx[6] * (xx[0] * xx[5] * xx[5] - xx[7]) + xx[3] * xx[10]) /
                   xx[9] + xx[12] * xx[2] / xx[9];
  xx[1] = xx[13] - xx[12];
  deriv[0] = state[1];
  deriv[1] = xx[13];
  deriv[2] = state[3];
  deriv[3] = (xx[1] - xx[0] * (xx[11] * xx[4] + xx[1] * xx[5]) * xx[5]) * xx[8];
}

void dae(adouble* deriv, adouble* path, adouble* state,
         adouble* input, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{
  bool auto_deriv = (workspace->algorithm->derivatives == "automatic")? true: false;

  if (auto_deriv)
  {
    adouble xx[14];
    xx[0] = 2.0;
    xx[1] = 0.125;
    xx[2] = 0.5;
    xx[3] = xx[2] * state[2];
    xx[4] = sin(xx[3]);
    xx[5] = cos(xx[3]);
    xx[3] = xx[0] * xx[5] * xx[4];
    xx[6] = xx[2] * xx[3];
    xx[7] = 0.03125000000000001;

    xx[8] = xx[1] / xx[7];
    xx[7] = 1.0;
    xx[9] = xx[0] * xx[4] * xx[4] - xx[7];
    xx[10] = (xx[2] - xx[1] * xx[8]) * xx[9];
    xx[2] = xx[0] + xx[6] * xx[3] + xx[10] * xx[9];
    memcpy(xx + 9, xx + 2, 1 * sizeof(double));

    xx[11] = 9.810000000000002;
    xx[12] = 1.77635683940025e-15;
    xx[13] = (input[0] + xx[0] * xx[1] * state[3] * state[3] * xx[4] * xx[5]) /
             xx[9] + xx[11] * (xx[6] * (xx[0] * xx[5] * xx[5] - xx[7]) + xx[3] * xx[10]) /
                     xx[9] + xx[12] * xx[2] / xx[9];
    xx[1] = xx[13] - xx[12];
    deriv[0] = state[1];
    deriv[1] = xx[13];
    deriv[2] = state[3];
    deriv[3] = (xx[1] - xx[0] * (xx[11] * xx[4] + xx[1] * xx[5]) * xx[5]) * xx[8];
  }
  else
  {
    double x[4], u[1], dx[4];
    for (int i=0; i<4; ++i)
    {
      x[i] = state[i].value();
    }
    u[0] = input[0].value();
    dynamics(x, u, dx);
    for (int i=0; i<4; ++i)
    {
      deriv[i] = dx[i];
    }
  }
}
