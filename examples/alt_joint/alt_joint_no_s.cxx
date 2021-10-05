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
 * \file   alt_joint.cxx
 * \author Ramkumar Natarajan (rnataraj@cs.cmu.edu)
 * \date   6/22/21
 */

//////////////////////////////////////////////////////////////////////////
//////////////////        alt_joint.cxx        //////////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:                 Two link arm problem      ////////////////
//////// Last modified:         04 January 2009           ////////////////
//////// Reference:             PROPT users guide         ////////////////
//////// (See PSOPT handbook for full reference)          ////////////////
//////////////////////////////////////////////////////////////////////////
////////     Copyright (c) Victor M. Becerra, 2009        ////////////////
//////////////////////////////////////////////////////////////////////////
//////// This is part of the PSOPT software library, which ///////////////
//////// is distributed under the terms of the GNU Lesser ////////////////
//////// General Public License (LGPL)                    ////////////////
//////////////////////////////////////////////////////////////////////////

#include "dynamics.h"
#include "psopt.h"
#include <Eigen/Core>

#define DOF 5


adouble qS[] = {0.0, 0.0, 0.0, 0.0, 0.0};
//adouble qS[] = {0.1,0.1,0.1,0.1,0.1};
//adouble qS[] = {0.5,0.5,0.5,0.5,0.5};
adouble qSd[] = {0.0, 0.0, 0.0, 0.0, 0.0};
//  adouble qG[] = {0.496349540849362,4.29968989868597,5.48318530717959,5.96902604182061,5.88318530717959,6.13318530717959,0.785398163397448};
//adouble qG[] = {0.496349540849362,4.29968989868597,5.48318530717959,5.96902604182061,5.88318530717959};
adouble qG[] = {0.1,0.1,0.1,0.1,0.1};
//adouble qG[] = {0.0, 0.0, 0.0, 0.0, 0.0};
adouble qGd[] = {0.0, 0.0, 0.0, 0.0, 0.0};

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the end point (Mayer) cost function //////////
//////////////////////////////////////////////////////////////////////////

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters,adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace)
{
//  return tf;
  double rho = 100.0;
//  double rho = 1.0;
  adouble diff[DOF];
  for (int i=0; i<DOF; ++i)
  {
    diff[i] = final_states[2*i] - qG[i];
  }
  return rho*dot(diff, diff, DOF);
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the integrand (Lagrange) cost function  //////
//////////////////////////////////////////////////////////////////////////

adouble integrand_cost(adouble* states, adouble* controls,
                       adouble* parameters, adouble& time, adouble* xad,
                       int iphase, Workspace* workspace)
{
  /// Minimize control
  adouble rho = 1.0;
  return rho*dot(controls, controls, DOF);

  /// Minimize weighted control
//  adouble cost = 0;
//  adouble rho[] = {1.0, 1.0, 50.0, 100.0, 200.0};
//  for (int i=0; i<DOF; ++i)
//  {
//    cost += rho[i]*controls[i]*controls[i];
//  }
//  return cost;

  /// Minimize EE dist
//  double rho = 1.0;
//  adouble diff[DOF];
//  for (int i=0; i<DOF; ++i)
//  {
//    diff[i] = states[2*i] - qG[i];
//  }
//  return rho*dot(diff, diff, DOF);
}


////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)
{
  int i;

  for (i=0; i< 2*DOF; i++ ) {
    e[ i ]     =  initial_states[ i ];
  }

  for (i=0; i< 2*DOF; i++ ) {
    e[2*DOF + i ] =  final_states[ i ];
  }

}

///////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function ///////////////////
///////////////////////////////////////////////////////////////////////////

void linkages( adouble* linkages, adouble* xad, Workspace* workspace)
{
  // No linkages as this is a single phase problem
}


////////////////////////////////////////////////////////////////////////////
///////////////////  Define the main routine ///////////////////////////////
////////////////////////////////////////////////////////////////////////////


int main(void) {

////////////////////////////////////////////////////////////////////////////
///////////////////  Declare key structures ////////////////////////////////
////////////////////////////////////////////////////////////////////////////

  Alg algorithm;
  Sol solution;
  Prob problem;

////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem name  ////////////////////////////////
////////////////////////////////////////////////////////////////////////////

  problem.name = "alt_joint";

  problem.outfilename = "alt_joint.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

  problem.nphases = 1;
  problem.nlinkages = 0;

  psopt_level1_setup(problem);


/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup /////////////
/////////////////////////////////////////////////////////////////////////////
  int nstates = 2*DOF;
  int ncontrols = DOF;
  int nevents = 4*DOF;
  int npath = 0;
  int nodes = 40;
  double t0 = 1.0;
  double tmax = 20.0;


  problem.phases(1).nstates = nstates;
  problem.phases(1).ncontrols = ncontrols;
  problem.phases(1).nevents = nevents;
  problem.phases(1).npath = npath;
  problem.phases(1).nodes << nodes;

  psopt_level2_setup(problem, algorithm);


////////////////////////////////////////////////////////////////////////////
///////////////////  Declare DMatrix objects to store results //////////////
////////////////////////////////////////////////////////////////////////////

  DMatrix x, u, t;
  DMatrix lambda, H;

////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////

  /// Bounds on state variables
  for (int i = 0; i < DOF; ++i)
  {
    problem.phases(1).bounds.lower.states(2*i) = -2*M_PI;
//    problem.phases(1).bounds.lower.states(2*i) = -std::numeric_limits<double>::infinity();
    problem.phases(1).bounds.lower.states(2*i + 1) = -100;

    problem.phases(1).bounds.upper.states(2*i) = 2*M_PI;
//    problem.phases(1).bounds.upper.states(2*i) = std::numeric_limits<double>::infinity();
    problem.phases(1).bounds.upper.states(2*i + 1) = 100;
  }

  /// Bounds on control variables
  for (int i = 0; i < DOF; ++i)
  {
    problem.phases(1).bounds.lower.controls(i) = -50;
    problem.phases(1).bounds.upper.controls(i) = 50;
  }

  /// Start specification (lower bound)
  for (int i = 0; i < DOF; ++i)
  {
    problem.phases(1).bounds.lower.events(2*i) = qS[i].value();
//    problem.phases(1).bounds.lower.events(2*i + 1) = qSd[i].value();
    problem.phases(1).bounds.lower.events(2*i + 1) = qSd[i].value()-100;
  }

  /// Goal specification (lower bound)
  for (int i = 0; i < DOF; ++i)
  {
    problem.phases(1).bounds.lower.events(2*DOF + 2*i) = qG[i].value();
//    problem.phases(1).bounds.lower.events(2*DOF + 2*i + 1) = qGd[i].value();
    problem.phases(1).bounds.lower.events(2*DOF + 2*i + 1) = qGd[i].value()-100;
  }

  /// Start specification (upper bound)
  for (int i = 0; i < DOF; ++i)
  {
    problem.phases(1).bounds.upper.events(2*i) = qS[i].value();
//    problem.phases(1).bounds.upper.events(2*i + 1) = qSd[i].value();
    problem.phases(1).bounds.upper.events(2*i + 1) = qSd[i].value()+100;
  }

  /// Goal specification (upper bound)
  for (int i = 0; i < DOF; ++i)
  {
    problem.phases(1).bounds.upper.events(2*DOF + 2*i) = qG[i].value();
//    problem.phases(1).bounds.upper.events(2*DOF + 2*i + 1) = qGd[i].value();
    problem.phases(1).bounds.upper.events(2*DOF + 2*i + 1) = qGd[i].value()+100;
  }

//  problem.phases(1).bounds.upper.events = problem.phases(1).bounds.lower.events;


  /// Start time specification
  problem.phases(1).bounds.lower.StartTime    = 0.0;
  problem.phases(1).bounds.upper.StartTime    = 0.0;

  /// Goal time specification
  problem.phases(1).bounds.lower.EndTime      = t0;
  problem.phases(1).bounds.upper.EndTime      = tmax;


////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem functions  ///////////////////////////
////////////////////////////////////////////////////////////////////////////


  problem.integrand_cost 				= &integrand_cost;
  problem.endpoint_cost 					= &endpoint_cost;
  problem.dae 								= &dae;
  problem.events 							= &events;
  problem.linkages							= &linkages;



////////////////////////////////////////////////////////////////////////////
///////////////////  Define & register initial guess ///////////////////////
////////////////////////////////////////////////////////////////////////////


  MatrixXd x0(nstates,nodes);
  for (int i=0; i<DOF; ++i)
  {
    x0.block(2*i,0,1,nodes) = linspace(qS[2*i].value(),qG[2*i].value(), nodes);
    x0.block(2*i+1,0,1,nodes) = linspace(qSd[2*i+1].value(),qGd[2*i+1].value(), nodes);
  }

  problem.phases(1).guess.controls       = zeros(ncontrols,nodes);
//  problem.phases(1).guess.controls       = GaussianRandom(ncontrols,nodes);
  problem.phases(1).guess.states         = x0;
  problem.phases(1).guess.time           = linspace(0, tmax, nodes);

////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////


  algorithm.nlp_method                  = "IPOPT";
  algorithm.scaling                     = "automatic";
  algorithm.derivatives                 = "automatic";
  algorithm.nlp_iter_max                = 60;
  algorithm.nlp_tolerance               = 1.e-2;

////////////////////////////////////////////////////////////////////////////
///////////////////  Now call PSOPT to solve the problem   /////////////////
////////////////////////////////////////////////////////////////////////////

  psopt(solution, problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////  Extract relevant variables from solution structure   //////////
////////////////////////////////////////////////////////////////////////////

  x 		= solution.get_states_in_phase(1);
  u 		= solution.get_controls_in_phase(1);
  t 		= solution.get_time_in_phase(1);

////////////////////////////////////////////////////////////////////////////
///////////  Save solution data to files if desired ////////////////////////
////////////////////////////////////////////////////////////////////////////

  Save(x,"x.dat");
  Save(u,"u.dat");
  Save(t,"t.dat");


////////////////////////////////////////////////////////////////////////////
///////////  Plot some results if desired (requires gnuplot) ///////////////
////////////////////////////////////////////////////////////////////////////

//  plot(t,x,problem.name + ": states", "time (s)", "states", "x1 x2 x3 x4 x5 x6 x7 x8 x9 x10");
  DMatrix z = zeros(DOF, nodes);
  for (int i=0; i<DOF; ++i)
  {
    z.block(i, 0, 1, nodes) = x.block(2*i, 0, 1, nodes);
  }
  plot(t,z,problem.name + ": states", "time (s)", "states", "x0 x2 x4 x6 x8");

  plot(t,u,problem.name + ": controls", "time (s)", "controls", "u1 u2 u3 u4 u5");


//  plot(t,x,problem.name + ": states", "time (s)", "states", "x1 x2 x3 x4 x5 x6 x7 x8 x9 x10",
//       "pdf", "alt_joint_states.pdf");
//
//  plot(t,u,problem.name + ": controls", "time (s)", "controls", "u1 u2 u3 u4 u5",
//       "pdf", "alt_joint_controls.pdf");

  std::cout << "+++++++++++++++++++++++ with init guess ++++++++++++++++++++++++++" << std::endl;

  problem.phases(1).guess.controls       = u;
  problem.phases(1).guess.states         = x;
  problem.phases(1).guess.time           = t;
  psopt(solution, problem, algorithm);

}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////

