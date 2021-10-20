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
 * \file   mujoco_interface.cxx
 * \author Ramkumar Natarajan (rnataraj@cs.cmu.edu)
 * \date   9/30/21
 */

#include "mujoco_interface.h"

#include <cassert>

MujocoPsopt::MujocoPsopt() : mj_model_(nullptr), mj_data_(nullptr)
{
  /// Activate MuJoCo
  mj_activate("/home/gaussian/cmu_ri_phd/phd_misc/mujoco200_linux/bin/mjkey.txt");
}

MujocoPsopt::MujocoPsopt(const char* filename) : mjcf_filename_(filename)
{
  /// Activate MuJoCo
  mj_activate("/home/gaussian/cmu_ri_phd/phd_misc/mujoco200_linux/bin/mjkey.txt");
  mj_model_ = mj_loadXML(mjcf_filename_, NULL, NULL, 0);
  mj_data_ = mj_makeData(mj_model_);
}

void MujocoPsopt::setupFromMJCFFile(const char *filename)
{
  // If the model and data are already allocated, first delete them
  if (mj_model_)
  {
    mj_deleteModel(mj_model_);
  }
  if (mj_data_)
  {
    mj_deleteData(mj_data_);
  }

  // Allocate new mj objects
  char loadxml_err[1000];
  mj_model_ = mj_loadXML(filename, NULL, loadxml_err, 1000);
  std::cout << loadxml_err << std::endl;
  mj_data_ = mj_makeData(mj_model_);
}

void MujocoPsopt::forwardSimulate(const adouble *x, const adouble *u, adouble *dx)
{
  if (!mj_model_ || !mj_data_)
  {
    std::cerr << "[MuJoCo Error] Mujoco models not set! Aborting"  << std::endl;
    std::abort();
  }

  assert(mj_model_->nq == mj_model_->nv);
  double q[mj_model_->nq];
  double dq[mj_model_->nv];
  double a[mj_model_->nu];

  for (int i=0; i<mj_model_->nq; ++i)
  {
    q[i] = x[2*i].value();
    dq[i+1] = x[2*i+1].value();
  }
  for (int i=0; i<mj_model_->nu; ++i)
  {
    a[i] = u[i].value();
  }

  // setting current state
  mju_copy(mj_data_->qpos, q, mj_model_->nq);
  mju_copy(mj_data_->qvel, dq, mj_model_->nv);
  // setting current control
  mju_copy(mj_data_->ctrl, a, mj_model_->nu);
  // Forward simulate
  mj_forward(mj_model_, mj_data_);
  // Get the derivatives
  for (int i=0; i<mj_model_->nq; ++i)
  {
    dx[2*i] = x[2*i+1];
    dx[2*i+1] = mj_data_->qacc[i];
  }
}

void MujocoPsopt::forwardSimulate(const double *x, const double *u, double *dx)
{
  if (!mj_model_ || !mj_data_)
  {
    std::cerr << "[MuJoCo Error] Mujoco models not set! Aborting"  << std::endl;
    std::abort();
  }

  assert(mj_model_->nq == mj_model_->nv);
  double q[mj_model_->nq];
  double dq[mj_model_->nv];
  double a[mj_model_->nu];
  double z1[mj_model_->nu];
  double z2[mj_model_->nbody*6];

  for (int i=0; i<mj_model_->nq; ++i)
  {
    q[i] = x[2*i];
    dq[i+1] = x[2*i+1];
  }
  for (int i=0; i<mj_model_->nu; ++i)
  {
    a[i] = u[i];
    z1[i] = 0;
  }
  for (int i=0; i<mj_model_->nbody; ++i)
  {
    for (int j=0; j<6; ++j)
    {
      z2[i*mj_model_->nbody + j] = 0;
    }
  }

  // setting current state
  mju_copy(mj_data_->qpos, q, mj_model_->nq);
  mju_copy(mj_data_->qvel, dq, mj_model_->nv);
  // setting current control
//  mju_copy(mj_data_->ctrl, a, mj_model_->nu);
  mju_copy(mj_data_->qfrc_applied, a, mj_model_->nv);
  mju_copy(mj_data_->ctrl, z1, mj_model_->nu);
  mju_copy(mj_data_->xfrc_applied, z2, mj_model_->nbody*6);


  // Forward simulate
  mj_forward(mj_model_, mj_data_);
  // Get the derivatives
  for (int i=0; i<mj_model_->nq; ++i)
  {
    dx[2*i] = x[2*i+1];
    dx[2*i+1] = mj_data_->qacc[i];
  }
}

void MujocoPsopt::inverseSimulate(const double *x, const double *v, const double *a, double *tau)
{
  // setting current state
  mju_copy(mj_data_->qpos, x, mj_model_->nq);
  mju_copy(mj_data_->qvel, v, mj_model_->nv);
  mju_copy(mj_data_->qacc, a, mj_model_->nv);

  // Inverse simulate
  mj_inverse(mj_model_, mj_data_);
  for (int i=0; i<mj_model_->nv; ++i)
  {
    tau[i] = mj_data_->qfrc_inverse[i];
  }
}

void MujocoPsopt::step(const double *x, const double *u, double *dx)
{

}

double MujocoPsopt::interp1(const double t0, const double tf, const double ts, const double x0, const double xf)
{
  double tratio = (ts-t0)/(tf-t0);
  assert(tratio<=1);
  return (1-tratio)*x0 + tratio*xf;
}


///// TORQUE
void MujocoPsopt::visualize(double T, std::vector<double> t, std::vector<std::vector<double>> u)
{
  mjvCamera cam;                      // abstract camera
  mjvOption opt;                      // visualization options
  mjvScene scn;                       // abstract scene
  mjrContext con;                     // custom GPU context

  // init GLFW, create window, make OpenGL context current, request v-sync
  glfwInit();
  GLFWwindow* window = glfwCreateWindow(1200, 900, "Demo", NULL, NULL);
  glfwMakeContextCurrent(window);
  glfwSwapInterval(1);

  // initialize visualization data structures
  mjv_defaultCamera(&cam);
//  mjv_defaultPerturb(&pert);
  mjv_defaultOption(&opt);
  mjv_defaultScene(&scn);
  mjr_defaultContext(&con);

  // create scene and context
  assert(mj_model_);
  mjv_makeScene(mj_model_, &scn, 1000);
  mjr_makeContext(mj_model_, &con, mjFONTSCALE_100);

  // init time to 0
  mj_data_->time = 0.0;
  int idx = 0;
  double curr_u[mj_model_->nv];
  for (int i=0; i<mj_model_->nv; ++i)
  {
    curr_u[i] = u[idx][i];
  }

  // zero vars
  double z1[mj_model_->nu];
//  double z2[mj_model_->nbody*6];
  for (int i=0; i<mj_model_->nu; ++i)
  {
    z1[i] = 0;
  }
  // init states
  mju_copy(mj_data_->qpos, z1, mj_model_->nq);
  mju_copy(mj_data_->qvel, z1, mj_model_->nv);

  // run main loop, target real-time simulation and 60 fps rendering
  while( !glfwWindowShouldClose(window) )
  {
    // Update idx and the current control input
    if (mj_data_->time > t[idx])
    {
      ++idx;
      for (int i=0; i<mj_model_->nv; ++i)
      {
        curr_u[i] = u[idx][i];
      }
    }

    // advance interactive simulation for 1/60 sec
    //  Assuming MuJoCo can simulate faster than real-time, which it usually can,
    //  this loop will finish on time for the next frame to be rendered at 60 fps.
    //  Otherwise add a cpu timer and exit this loop when it is time to render.
    mjtNum simstart = mj_data_->time;
//    std::cout << "simstart: " << simstart << std::endl;
    double frame_rate = 1.0/60.0;
    // CPU time
    Clock::time_point time_now = Clock::now();
    while( mj_data_->time - simstart < frame_rate )
    {
      // assign control
      mju_copy(mj_data_->qfrc_applied, curr_u, mj_model_->nv);
      mju_copy(mj_data_->ctrl, z1, mj_model_->nu);
      // Update idx and the current control input
      if (mj_data_->time > t[idx])
      {
        ++idx;
        for (int i=0; i<mj_model_->nv; ++i)
        {
          curr_u[i] = u[idx][i];
        }
      }

      mj_step(mj_model_, mj_data_);
//      std::cout << mj_data_->time << std::endl;
    }
    double time_left = frame_rate-duration_cast<duration<double> >(Clock::now() - time_now).count();
    while(time_left > 0)
    {
      time_left = frame_rate-duration_cast<duration<double> >(Clock::now() - time_now).count();
    }

    // get framebuffer viewport
    mjrRect viewport = {0, 0, 0, 0};
    glfwGetFramebufferSize(window, &viewport.width, &viewport.height);

    // update scene and render
    mjv_updateScene(mj_model_, mj_data_, &opt, NULL, &cam, mjCAT_ALL, &scn);
    mjr_render(viewport, &scn, &con);

    // swap OpenGL buffers (blocking call due to v-sync)
    glfwSwapBuffers(window);

    // process pending GUI events, call GLFW callbacks
    glfwPollEvents();
  }

// close GLFW, free visualization storage
  glfwTerminate();
  mjv_freeScene(&scn);
  mjr_freeContext(&con);

}

//void MujocoPsopt::visualize(double T, std::vector<double> t, std::vector<std::vector<double>> x,
//                            std::vector<std::vector<double>> dx)
//{
//  mjvCamera cam;                      // abstract camera
//  mjvOption opt;                      // visualization options
//  mjvScene scn;                       // abstract scene
//  mjrContext con;                     // custom GPU context
//
//  // init GLFW, create window, make OpenGL context current, request v-sync
//  glfwInit();
//  GLFWwindow* window = glfwCreateWindow(1200, 900, "Demo", NULL, NULL);
//  glfwMakeContextCurrent(window);
//  glfwSwapInterval(1);
//
//  // initialize visualization data structures
//  mjv_defaultCamera(&cam);
////  mjv_defaultPerturb(&pert);
//  mjv_defaultOption(&opt);
//  mjv_defaultScene(&scn);
//  mjr_defaultContext(&con);
//
//  // create scene and context
//  assert(mj_model_);
//  mjv_makeScene(mj_model_, &scn, 1000);
//  mjr_makeContext(mj_model_, &con, mjFONTSCALE_100);
//
//  // init time to 0
//  mj_data_->time = 0.0;
//  int idx = 0;
//  double curr_x[mj_model_->nq];
//  double curr_dx[mj_model_->nv];
//  for (int i=0; i<mj_model_->nq; ++i)
//  {
//    curr_x[i] = 0.0;
//    curr_dx[i] = 0.0;
//  }
//
//  double frame_rate = 1.0/60.0;
//
//  // run main loop, target real-time simulation and 60 fps rendering
//  while( !glfwWindowShouldClose(window) )
//  {
//    if (mj_data_->time > T)
//    {
//      break;
//    }
//
//    // Update idx and the current control input
//    if (mj_data_->time > t[idx+1])
//    {
//      ++idx;
//    }
//
//    // advance interactive simulation for 1/60 sec
//    //  Assuming MuJoCo can simulate faster than real-time, which it usually can,
//    //  this loop will finish on time for the next frame to be rendered at 60 fps.
//    //  Otherwise add a cpu timer and exit this loop when it is time to render.
//    for (int i=0; i<mj_model_->nv; ++i)
//    {
//      curr_x[i] = interp1(t[idx], t[idx+1], mj_data_->time, x[idx][i], x[idx+1][i]);
//      curr_dx[i] = interp1(t[idx], t[idx+1], mj_data_->time, dx[idx][i], dx[idx+1][i]);
//    }
//    // setting current state
//    mju_copy(mj_data_->qpos, curr_x, mj_model_->nq);
//    mju_copy(mj_data_->qvel, curr_dx, mj_model_->nv);
//    mj_data_->time += frame_rate;
//
//    std::cout << "idx: " << idx << " t: " << mj_data_->time << std::endl;
//    std::cout << "x: ";
//    for (int i=0; i<mj_model_->nv; ++i)
//    {
//      std::cout << mj_data_->qpos[i] << "\t";
//    }
//    std::cout << std::endl;
//    std::cout << "dx: ";
//    for (int i=0; i<mj_model_->nv; ++i)
//    {
//      std::cout << mj_data_->qvel[i] << "\t";
//    }
//    std::cout << std::endl;
//
//    // get framebuffer viewport
//    mjrRect viewport = {0, 0, 0, 0};
//    glfwGetFramebufferSize(window, &viewport.width, &viewport.height);
//
//    // update scene and render
//    mjv_updateScene(mj_model_, mj_data_, &opt, NULL, &cam, mjCAT_ALL, &scn);
//    mjr_render(viewport, &scn, &con);
//
//    // swap OpenGL buffers (blocking call due to v-sync)
//    glfwSwapBuffers(window);
//
//    // process pending GUI events, call GLFW callbacks
//    glfwPollEvents();
//  }
//
//// close GLFW, free visualization storage
//  glfwTerminate();
//  mjv_freeScene(&scn);
//  mjr_freeContext(&con);
//}

void MujocoPsopt::visualize(double T, std::vector<double> t, std::vector<std::vector<double>> x,
                            std::vector<std::vector<double>> dx)
{
  mjvCamera cam;                      // abstract camera
  mjvOption opt;                      // visualization options
  mjvScene scn;                       // abstract scene
  mjrContext con;                     // custom GPU context

  // init GLFW, create window, make OpenGL context current, request v-sync
  glfwInit();
  GLFWwindow* window = glfwCreateWindow(1200, 900, "Demo", NULL, NULL);
  glfwMakeContextCurrent(window);
  glfwSwapInterval(1);

  // initialize visualization data structures
  mjv_defaultCamera(&cam);
//  mjv_defaultPerturb(&pert);
  mjv_defaultOption(&opt);
  mjv_defaultScene(&scn);
  mjr_defaultContext(&con);

  // create scene and context
  mjv_makeScene(mj_model_, &scn, 1000);
  mjr_makeContext(mj_model_, &con, mjFONTSCALE_100);

  // init time to 0
  mj_data_->time = 0.0;
  int idx = 0;
  double curr_x[mj_model_->nq];
  double curr_dx[mj_model_->nv];
  double curr_ddx[mj_model_->nv];
  for (int i=0; i<mj_model_->nv; ++i)
  {
    curr_x[i] = x[idx][i];
    curr_dx[i] = dx[idx][i];
    curr_ddx[i] = 0;
  }

  // zero vars
  double z1[mj_model_->nu];
  for (int i=0; i<mj_model_->nu; ++i)
  {
    z1[i] = 0;
  }
  // init states
  mju_copy(mj_data_->qpos, z1, mj_model_->nq);
  mju_copy(mj_data_->qvel, z1, mj_model_->nv);

  // run main loop, target real-time simulation and 60 fps rendering
//  while( !glfwWindowShouldClose(window) )
  while(mj_data_->time < T)
  {
    // Update idx and the current control input
    if (mj_data_->time > t[idx+1])
    {
      ++idx;
    }

    // advance interactive simulation for 1/60 sec
    //  Assuming MuJoCo can simulate faster than real-time, which it usually can,
    //  this loop will finish on time for the next frame to be rendered at 60 fps.
    //  Otherwise add a cpu timer and exit this loop when it is time to render.
    mjtNum simstart = mj_data_->time;
    double frame_rate = 1.0/60.0;
    // CPU time
    Clock::time_point time_now = Clock::now();
    while( mj_data_->time - simstart < frame_rate )
    {
      for (int i=0; i<mj_model_->nv; ++i)
      {
        curr_x[i] = interp1(t[idx], t[idx + 1], mj_data_->time, x[idx][i], x[idx + 1][i]);
        curr_dx[i] = interp1(t[idx], t[idx + 1], mj_data_->time, dx[idx][i], dx[idx + 1][i]);
        curr_ddx[i] = (dx[idx+1][i]-dx[idx][i])/mj_model_->opt.timestep;
      }

      // assign control
//      mju_copy(mj_data_->ctrl, curr_x, mj_model_->nv);

      // Do ID
      double tau[mj_model_->nv];
      inverseSimulate(curr_x, curr_dx, curr_ddx, tau);
      mju_copy(mj_data_->qfrc_applied, tau, mj_model_->nq);
//      mju_copy(mj_data_->ctrl, tau, mj_model_->nq);
//      std::cout << "tau: " << std::endl;
//      for (int i=0; i<mj_model_->nq; ++i)
//      {
//        std::cout << tau[i] << "\t";
//      }
//      std::cout << std::endl;

      // Update idx and the current control input
      mj_step(mj_model_, mj_data_);
      if (mj_data_->time > t[idx+1])
      {
        ++idx;
      }

    }
//    double time_left = frame_rate-duration_cast<duration<double> >(Clock::now() - time_now).count();
//    while(time_left > 0)
//    {
//      time_left = frame_rate-duration_cast<duration<double> >(Clock::now() - time_now).count();
//    }

    // get framebuffer viewport
    mjrRect viewport = {0, 0, 0, 0};
    glfwGetFramebufferSize(window, &viewport.width, &viewport.height);

    // update scene and render
    mjv_updateScene(mj_model_, mj_data_, &opt, NULL, &cam, mjCAT_ALL, &scn);
    mjr_render(viewport, &scn, &con);

    // swap OpenGL buffers (blocking call due to v-sync)
    glfwSwapBuffers(window);

    // process pending GUI events, call GLFW callbacks
    glfwPollEvents();
  }

  std::cout << "final q: ";
  for (int i=0; i<mj_model_->nq; ++i)
  {
    std::cout << mj_data_->qpos << "\t " << std::endl;
  }
  std::cout << std::endl;


// close GLFW, free visualization storage
  glfwTerminate();
  mjv_freeScene(&scn);
  mjr_freeContext(&con);

}

