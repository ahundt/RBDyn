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

InverseStatics::InverseStatics(const MultiBody& mb):
	f_(mb.nrBodies())
{
}

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

  //std::cout << "mb.bodies().size() = \n" << mb.bodies().size() << std::endl;
  //std::cout << "mb.joints().size() = \n" << mb.joints().size() << std::endl;
  //std::cout << "mb.predecessors().size() = \n" << mb.predecessors().size() << std::endl;

  //for (size_t i = 0; i < mbc.q.size(); ++i)
    //for (size_t j = 0; j < mbc.q[i].size(); ++j)
      //std::cout << "mbc.q[" << i << "][" << j << "]: " << mbc.q[i][j] << std::endl;

  //std::cout << "a_0 = \n" << a_0 << std::endl;
  //for (size_t i = 0; i < mbc.force.size(); ++i)
    //std::cout << "force[" << i << "]: " << mbc.force[i] << std::endl;
  //for (size_t i = 0; i < f_.size(); ++i)
    //std::cout << "f_[" << i << "]: " << f_[i] << std::endl;
  //for (size_t i = 0; i < mbc.jointTorque.size(); ++i)
    //for (size_t j = 0; j < mbc.jointTorque[i].size(); ++j)
      //std::cout << "torque[" << i << "][" << j << "]:" << mbc.jointTorque[i][j] << std::endl;
}

void InverseStatics::computeTorqueJacobians(const MultiBody& mb, MultiBodyConfig& mbc)
{
}

void InverseStatics::computeTorqueJacobiansFD(const MultiBody& mb, MultiBodyConfig& mbc, double delta)
{
  MultiBodyConfig mbcCopy(mbc);

  size_t index = 0;
  for (size_t i = 0; i < mbc.q.size(); ++i)
    for (size_t j = 0; j < mbc.q[i].size(); ++j)
    {
      mbcCopy.q[i][j] += delta;
      forwardKinematics(mb, mbcCopy);
      forwardVelocity(mb, mbcCopy);
      inverseStatics(mb, mbcCopy);
      for (size_t k = 0; k < mbc.jointTorque.size(); ++k)
        for (size_t l = 0; l < mbc.jointTorque[k].size(); ++l)
          mbc.jointTorqueJacQ[k][l][index] = (mbcCopy.jointTorque[k][l]-mbc.jointTorque[k][l])/delta;
      mbcCopy.q[i][j] -= delta;
      index++;
    }
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


const std::vector<sva::ForceVecd>& InverseStatics::f() const
{
	return f_;
}

} // namespace rbd
