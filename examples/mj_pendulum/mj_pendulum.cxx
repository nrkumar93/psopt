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
 * \file   mj_pendulum.cxx
 * \author Ramkumar Natarajan (rnataraj@cs.cmu.edu)
 * \date   7/4/21
 */

//////////////////////////////////////////////////////////////////////////
//////////////////        twolinkarm.cxx        //////////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title:                 mj_pendulum               ////////////////
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

#include "psopt.h"

#define DOF 1
#define N_X 2
#define N_U 1

adouble qS[] = {0.0};
adouble qSd[] = {0.0};
adouble qG[] = {M_PI};
adouble qGd[] = {0.0};

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
  double rho = 1.0;
  return  rho*dot(controls, controls, N_U);
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the DAE's ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{
  workspace->problem->mj_handle.forwardSimulate(states, controls, derivatives);
//  std::cout << "deriv: " <<  derivatives[0] << "\t" << derivatives[1] << std::endl;
}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)
{
  int i;

  for (i=0; i< N_X; i++ ) {
    e[ i ]     =  initial_states[ i ];
  }

  for (i=0; i< N_X; i++ ) {
    e[N_X + i ] =  final_states[ i ];
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


int main(void)
{

////////////////////////////////////////////////////////////////////////////
///////////////////  Declare key structures ////////////////////////////////
////////////////////////////////////////////////////////////////////////////

  Alg  algorithm;
  Sol  solution;
  Prob problem;

////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem name  ////////////////////////////////
////////////////////////////////////////////////////////////////////////////

  problem.name        						= "mj_pendulum";

  problem.outfilename                 	= "mj_pendulum.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

  problem.nphases   							= 1;
  problem.nlinkages                   	= 0;

  psopt_level1_setup(problem);


/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup /////////////
/////////////////////////////////////////////////////////////////////////////
  int nodes = 40;
  double t0 = 0.1;
  double tmax = 3.0;

  problem.phases(1).nstates   				= N_X;
  problem.phases(1).ncontrols 				= N_U;
  problem.phases(1).nevents   				= 2*N_X;
  problem.phases(1).npath     				= 0;
  problem.phases(1).nodes               << nodes;

  psopt_level2_setup(problem, algorithm);

////////////////////////////////////////////////////////////////////////////
////////////////////////  Setup Mujoco interface ///////////////////////////
////////////////////////////////////////////////////////////////////////////
  const char* mjcf_file = "/home/gaussian/cmu_ri_phd/phd_misc/mujoco200_linux/model/alt_joint/pendulum.xml";
  problem.mj_handle.setupFromMJCFFile(mjcf_file);
  problem.mj_backend = true;

////////////////////////////////////////////////////////////////////////////
///////////////////  Declare DMatrix objects to store results //////////////
////////////////////////////////////////////////////////////////////////////

  DMatrix x, u, t;
  DMatrix lambda, H;

////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////

  problem.phases(1).bounds.lower.states(0) = -100;
  problem.phases(1).bounds.lower.states(1) = -100;

  problem.phases(1).bounds.upper.states(0) = 100;
  problem.phases(1).bounds.upper.states(1) = 100;

  problem.phases(1).bounds.lower.controls(0) = -100;
  problem.phases(1).bounds.upper.controls(0) = 100;

  problem.phases(1).bounds.lower.events(0) = 0.0;
  problem.phases(1).bounds.lower.events(1) = 0.0;
  problem.phases(1).bounds.lower.events(2) = M_PI;
  problem.phases(1).bounds.lower.events(3) = 0.0;

  problem.phases(1).bounds.upper.events = problem.phases(1).bounds.lower.events;



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


  MatrixXd x0(N_X,nodes);

  x0 <<  linspace(0.0,M_PI, nodes),
          linspace(0.0,0.0, nodes);

  problem.phases(1).guess.controls       = zeros(N_U,nodes);
  problem.phases(1).guess.states         = x0;
  problem.phases(1).guess.time           = linspace(0.0, tmax, nodes);

////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////


  algorithm.nlp_method                  = "IPOPT";
  algorithm.scaling                     = "automatic";
//  algorithm.derivatives                 = "automatic";
  algorithm.derivatives                 = "numerical";
  algorithm.nlp_iter_max                = 1000;
  algorithm.nlp_tolerance               = 1.e-6;

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

  plot(t,x,problem.name + ": states", "time (s)", "states", "x1 x2");

  plot(t,u,problem.name + ": controls", "time (s)", "controls", "u1");


  plot(t,x,problem.name + ": states", "time (s)", "states", "x1 x2",
       "pdf", "mj_pendulum.pdf");

  plot(t,u,problem.name + ": controls", "time (s)", "controls", "u1",
       "pdf", "mj_pendulum_controls.pdf");


}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////


