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
#include "util.hh"
#include "fixture.hh"

// arm
#include "XXXarm.h"

namespace rbd
{
void test(boost::shared_ptr<boost::test_tools::output_test_stream> output,
          rbd::MultiBody& mb, rbd::MultiBodyConfig& mbc,
          rbd::InverseStatics& IS, double q1, double q2, double q3)
{
  mbc.q[1][0] = q1;
  mbc.q[2][0] = q2;
  mbc.q[3][0] = q3;
  forwardKinematics(mb, mbc);
  forwardVelocity(mb, mbc);
  IS.inverseStatics(mb, mbc);
  IS.computeTorqueJacobianFD(mb, mbc, 1e-10);

  (*output) << "========================================" << std::endl;
  (*output) << "Results for mbc.q =" << mbc.q << std::endl;
  (*output) << "mbc.jointTorque =\n" << mbc.jointTorque << std::endl;
  for (auto e : IS.f())
    (*output) << "IS.f().vector =\n" << e.vector() << std::endl;
  (*output) << "IS.jointTorqueJacQ() =\n" << IS.jointTorqueJacQ() << std::endl;
  (*output) << "IS.jointTorqueJacF() =\n" << IS.jointTorqueJacF() << std::endl;
}

BOOST_AUTO_TEST_CASE(XXXArmTorqueJacobian)
{
  boost::shared_ptr<boost::test_tools::output_test_stream> output =
      retrievePattern("InverseStaticsTest");

  rbd::MultiBody mb;
  rbd::MultiBodyConfig mbc;
  rbd::MultiBodyGraph mbg;
  std::tie(mb, mbc, mbg) = makeXXXarm();
  rbd::InverseStatics IS(mb);

  test(output, mb, mbc, IS, 0, 0, 0);
  test(output, mb, mbc, IS, M_PI, 0, 0);
  test(output, mb, mbc, IS, M_PI, M_PI / 2, 0);
  test(output, mb, mbc, IS, 0.4, 0.1, 0.2);

  std::cout << output->str() << std::endl;
  BOOST_CHECK(output->match_pattern());
}  // end of namespace rbd
}

