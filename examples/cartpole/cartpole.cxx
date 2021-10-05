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
 * \file   cartpole.cxx
 * \author Ramkumar Natarajan (rnataraj@cs.cmu.edu)
 * \date   7/4/21
 */

//////////////////////////////////////////////////////////////////////////
//////////////////        twolinkarm.cxx        //////////////////////////
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

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the end point (Mayer) cost function //////////
//////////////////////////////////////////////////////////////////////////

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters,adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace)
{
  return tf;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the integrand (Lagrange) cost function  //////
//////////////////////////////////////////////////////////////////////////

adouble integrand_cost(adouble* states, adouble* controls,
                       adouble* parameters, adouble& time, adouble* xad,
                       int iphase, Workspace* workspace)
{
  double rho = 1.0;
  return  rho*dot(controls, controls, 1);
}


////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)
{
  int i;

  for (i=0; i< 4; i++ ) {
    e[ i ]     =  initial_states[ i ];
  }

  for (i=0; i< 4; i++ ) {
    e[4 + i ] =  final_states[ i ];
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

  problem.name        						= "cartpole";

  problem.outfilename                 	= "cartpole.txt";

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

  problem.phases(1).nstates   				= 4;
  problem.phases(1).ncontrols 				= 1;
  problem.phases(1).nevents   				= 8;
  problem.phases(1).npath     				= 0;
  problem.phases(1).nodes               << nodes;

  psopt_level2_setup(problem, algorithm);

////////////////////////////////////////////////////////////////////////////
////////////////////////  Setup Mujoco interface ///////////////////////////
////////////////////////////////////////////////////////////////////////////
  const char* mjcf_file = "/home/gaussian/cmu_ri_phd/phd_misc/mujoco200_linux/model/cartpole.xml";
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

  problem.phases(1).bounds.lower.states(0) = -1.6;
  problem.phases(1).bounds.lower.states(1) = -1000;
  problem.phases(1).bounds.lower.states(2) = -2*M_PI;
  problem.phases(1).bounds.lower.states(3) = -1000;

  problem.phases(1).bounds.upper.states(0) = 1.6;
  problem.phases(1).bounds.upper.states(1) = 1000;
  problem.phases(1).bounds.upper.states(2) = 2*M_PI;
  problem.phases(1).bounds.upper.states(3) = 1000;

  problem.phases(1).bounds.lower.controls(0) = -100;
  problem.phases(1).bounds.upper.controls(0) = 100;

  problem.phases(1).bounds.lower.events(0) = 0.0;
  problem.phases(1).bounds.lower.events(1) = 0.0;
  problem.phases(1).bounds.lower.events(2) = 0.0;
  problem.phases(1).bounds.lower.events(3) = 0.0;
  problem.phases(1).bounds.lower.events(4) = 0.8;
  problem.phases(1).bounds.lower.events(5) = 0.0;
  problem.phases(1).bounds.lower.events(6) = M_PI;
  problem.phases(1).bounds.lower.events(7) = 0.0;

  problem.phases(1).bounds.upper.events = problem.phases(1).bounds.lower.events;



  problem.phases(1).bounds.lower.StartTime    = 0.0;
  problem.phases(1).bounds.upper.StartTime    = 0.0;

  problem.phases(1).bounds.lower.EndTime      = 0.1;
  problem.phases(1).bounds.upper.EndTime      = 2.0;


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


  MatrixXd x0(4,nodes);

  x0 <<  linspace(0.0,0.8, nodes),
          linspace(0.0,0.0, nodes),
          linspace(0.0,M_PI, nodes),
          linspace(0.0,0.0, nodes);

  problem.phases(1).guess.controls       = zeros(1,nodes);
  problem.phases(1).guess.states         = x0;
  problem.phases(1).guess.time           = linspace(0.0, 2.0, nodes);

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

  plot(t,x,problem.name + ": states", "time (s)", "states", "x1 x2 x3 x4");

  plot(t,u,problem.name + ": controls", "time (s)", "controls", "u1");


  plot(t,x,problem.name + ": states", "time (s)", "states", "x1 x2 x3 x4",
       "pdf", "cartpole_states.pdf");

  plot(t,u,problem.name + ": controls", "time (s)", "controls", "u1",
       "pdf", "cartpole_controls.pdf");


}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////


