// This file is part of RBDyn.
//
// RBDyn is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// RBDyn is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with RBDyn.  If not, see <http://www.gnu.org/licenses/>.

// check memory allocation in some method

// includes
// std
#include <iostream>

// boost
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Statics
#include <boost/test/unit_test.hpp>
#include <boost/math/constants/constants.hpp>

// SpaceVecAlg
#include <SpaceVecAlg/SpaceVecAlg>

// RBDyn
#include "FK.h"
#include "FV.h"
#include "FD.h"
#include "IS.h"
#include "Body.h"
#include "Joint.h"
#include "MultiBody.h"
#include "MultiBodyConfig.h"
#include "MultiBodyGraph.h"

// arm
#include "XXXarm.h"

const double TOL = 0.0000001;

BOOST_AUTO_TEST_CASE(XXXArmInvStatic)
{
	rbd::MultiBody mb;
	rbd::MultiBodyConfig mbc;
	rbd::MultiBodyGraph mbg;
	std::tie(mb, mbc, mbg) = makeXXXarm();
  rbd::InverseStatics IS(mb);
	mbc.q[1][0] = 0.;
	mbc.q[2][0] = 0.;
	mbc.q[3][0] = 0.;
	forwardKinematics(mb, mbc);
	forwardVelocity(mb, mbc);
  IS.inverseStatics(mb, mbc);

  {
    //Expected internal forces
    std::vector<Eigen::Vector6d> expF;
    expF.resize(4);
    expF[0] << 0, 0, 0, 0, 39.24, 0;
    expF[1] << 0, 0, 0, 0, 29.43, 0;
    expF[2] << 0, 0, 0, 0, 19.62, 0;
    expF[3] << 0, 0, 0, 0, 9.81, 0;

    //Expected torques
    std::vector<double> expT;
    expT.push_back(0);
    expT.push_back(0);
    expT.push_back(0);

    for (size_t i = 0; i < IS.f().size(); ++i)
      for (size_t j = 0; j < 6; ++j)
      {
        double res = fabs(IS.f()[i].vector()[j] - expF[i][j]);
        BOOST_CHECK_SMALL(res, TOL);
      }

    for (size_t i = 0; i < mbc.jointTorque.size(); ++i)
      for (size_t j = 0; j < mbc.jointTorque[i].size(); ++j)
      {
        double res = fabs(mbc.jointTorque[i][j] - expT[i]);
        BOOST_CHECK_SMALL(res, TOL);
      }
  }

	mbc.q[1][0] = M_PI;
	mbc.q[2][0] = 0.;
	mbc.q[3][0] = 0.;
	forwardKinematics(mb, mbc);
	forwardVelocity(mb, mbc);
  IS.inverseStatics(mb, mbc);

  {
    //Expected internal forces
    std::vector<Eigen::Vector6d> expF;
    expF.resize(4);
    expF[0] << 0, 0, 0, 0, 39.24, 0;
    expF[1] << 0, 0, 0, 0, -29.43, 0;
    expF[2] << 0, 0, 0, 0, -19.62, 0;
    expF[3] << 0, 0, 0, 0, -9.81, 0;

    //Expected torques
    std::vector<double> expT;
    expT.push_back(0);
    expT.push_back(0);
    expT.push_back(0);

    for (size_t i = 0; i < IS.f().size(); ++i)
      for (size_t j = 0; j < 6; ++j)
      {
        double res = fabs(IS.f()[i].vector()[j] - expF[i][j]);
        BOOST_CHECK_SMALL(res, TOL);
      }

    for (size_t i = 0; i < mbc.jointTorque.size(); ++i)
      for (size_t j = 0; j < mbc.jointTorque[i].size(); ++j)
      {
        double res = fabs(mbc.jointTorque[i][j] - expT[i]);
        BOOST_CHECK_SMALL(res, TOL);
      }
  }

	mbc.q[1][0] = M_PI/2;
	mbc.q[2][0] = 0.;
	mbc.q[3][0] = 0.;
	forwardKinematics(mb, mbc);
	forwardVelocity(mb, mbc);
  IS.inverseStatics(mb, mbc);

  {
    //Expected internal forces
    std::vector<Eigen::Vector6d> expF;
    expF.resize(4);
    expF[0] << -44.145, 0, 0, 0, 39.24, 0;
    expF[1] << -44.145, 0, 0, 0, 0, -29.43;
    expF[2] << -19.62, 0, 0, 0, 0, -19.62;
    expF[3] << -4.905, 0, 0, 0, 0, -9.81;

    //Expected torques
    std::vector<double> expT;
    expT.push_back(0);
    expT.push_back(-44.145);
    expT.push_back(-19.62 );
    expT.push_back(-4.905 );

    for (size_t i = 0; i < IS.f().size(); ++i)
      for (size_t j = 0; j < 6; ++j)
      {
        double res = fabs(IS.f()[i].vector()[j] - expF[i][j]);
        BOOST_CHECK_SMALL(res, TOL);
      }

    for (size_t i = 0; i < mbc.jointTorque.size(); ++i)
      for (size_t j = 0; j < mbc.jointTorque[i].size(); ++j)
      {
        double res = fabs(mbc.jointTorque[i][j] - expT[i]);
        BOOST_CHECK_SMALL(res, TOL);
      }
  }
}
