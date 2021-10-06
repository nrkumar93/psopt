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

  for (int i=0; i<mj_model_->nq; ++i)
  {
    q[i] = x[2*i];
    dq[i+1] = x[2*i+1];
  }
  for (int i=0; i<mj_model_->nu; ++i)
  {
    a[i] = u[i];
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
