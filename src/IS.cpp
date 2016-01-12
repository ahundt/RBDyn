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
#include <iostream>

// associated header
#include "IS.h"

// includes
// RBDyn
#include "util.hh"
#include "MultiBody.h"
#include "MultiBodyConfig.h"
#include "FK.h"
#include "FV.h"

namespace rbd
{
InverseStatics::InverseStatics(const MultiBody& mb) : f_(mb.nrBodies()) {}

void InverseStatics::inverseStatics(const MultiBody& mb, MultiBodyConfig& mbc)
{
	const std::vector<Body>& bodies = mb.bodies();
	const std::vector<Joint>& joints = mb.joints();
	const std::vector<int>& pred = mb.predecessors();

	sva::MotionVecd a_0(Eigen::Vector3d::Zero(), mbc.gravity);

	for(std::size_t i = 0; i < bodies.size(); ++i)
	{
		const sva::PTransformd& X_p_i = mbc.parentToSon[i];

		if(pred[i] != -1)
			mbc.bodyAccB[i] = X_p_i*mbc.bodyAccB[pred[i]];
		else
			mbc.bodyAccB[i] = X_p_i*a_0;

		f_[i] = bodies[i].inertia()*mbc.bodyAccB[i] -
			mbc.bodyPosW[i].dualMul(mbc.force[i]);
	}


	for(int i = static_cast<int>(bodies.size()) - 1; i >= 0; --i)
	{
		for(int j = 0; j < joints[i].dof(); ++j)
		{
			mbc.jointTorque[i][j] = mbc.motionSubspace[i].col(j).transpose()*
				f_[i].vector();
		}

		if(pred[i] != -1)
		{
			const sva::PTransformd& X_p_i = mbc.parentToSon[i];
			f_[pred[i]] = f_[pred[i]] + X_p_i.transMul(f_[i]);
		}
	}
}

void InverseStatics::computeTorqueJacobian(const MultiBody& /*mb*/,
                                            MultiBodyConfig& /*mbc*/)
{
  std::cerr << "InverseStatics::computeTorqueJacobian is not implemented yet" << std::endl;
}

void printMBC(const MultiBody& mb, const MultiBodyConfig& mbc)
{
  std::cout << "mb.bodies() = " << mb.bodies() << std::endl;
  std::cout << "mb.joints() = " << mb.joints() << std::endl;
  std::cout << "mb.predecessors() = " << mb.predecessors() << std::endl;
  std::cout << "mbc.gravity = " << mbc.gravity << std::endl;
  std::cout << "mbc.parentToSon = \n" << mbc.parentToSon << std::endl;
  std::cout << "mbc.bodyPosW = \n" << mbc.bodyPosW << std::endl;
  std::cout << "mbc.force = \n" << mbc.force << std::endl;
  std::cout << "mbc.motionSubspace = \n" << mbc.motionSubspace << std::endl;
}

void InverseStatics::sInverseStatics(const MultiBody& mb, MultiBodyConfig& mbc)
{
	checkMatchAlphaD(mb, mbc);
	checkMatchForce(mb, mbc);
	checkMatchJointConf(mb, mbc);
	checkMatchJointVelocity(mb, mbc);
	checkMatchBodyPos(mb, mbc);
	checkMatchParentToSon(mb, mbc);
	checkMatchBodyVel(mb, mbc);
	checkMatchMotionSubspace(mb, mbc);

	checkMatchBodyAcc(mb, mbc);
	checkMatchJointTorque(mb, mbc);

  inverseStatics(mb, mbc);
}

} // namespace rbd
