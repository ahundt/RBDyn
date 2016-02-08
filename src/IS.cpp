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
// sva
#include <SpaceVecAlg/SpaceVecAlg>

// RBDyn
#include "util.hh"
#include "MultiBody.h"
#include "MultiBodyConfig.h"
#include "FK.h"
#include "FV.h"
#include "Jacobian.h"

namespace rbd
{
InverseStatics::InverseStatics(const MultiBody& mb)
    : f_(mb.nrBodies()), df_(mb.nrBodies()), jointTorqueDiff_(mb.nrJoints())
{
  for (size_t i = 0; i < static_cast<size_t>(mb.nrBodies()); ++i)
  {
    df_[i] = Eigen::MatrixXd::Zero(6, mb.nrDof());
  }
  for (size_t i = 0; i < static_cast<size_t>(mb.nrJoints()); ++i)
  {
    jointTorqueDiff_[i].resize(mb.joint(static_cast<int>(i)).dof());
    for (int j = 0; j < mb.joint(static_cast<int>(i)).dof(); ++j)
    {
      jointTorqueDiff_[i][static_cast<size_t>(j)] =
          Eigen::VectorXd::Zero(mb.nrDof());
    }
  }
}

void InverseStatics::inverseStatics(const MultiBody& mb, MultiBodyConfig& mbc)
{
  const std::vector<Body>& bodies = mb.bodies();
  const std::vector<Joint>& joints = mb.joints();
  const std::vector<int>& pred = mb.predecessors();

  sva::MotionVecd a_0(Eigen::Vector3d::Zero(), mbc.gravity);

  for (std::size_t i = 0; i < bodies.size(); ++i)
  {
    mbc.bodyAccB[i] = mbc.bodyPosW[i] * a_0;

    f_[i] = bodies[i].inertia() * mbc.bodyAccB[i] -
            mbc.bodyPosW[i].dualMul(mbc.force[i]);
  }

  for (int i = static_cast<int>(bodies.size()) - 1; i >= 0; --i)
  {
    for (int j = 0; j < joints[i].dof(); ++j)
    {
      mbc.jointTorque[i][j] =
          mbc.motionSubspace[i].col(j).transpose() * f_[i].vector();
    }

    if (pred[i] != -1)
    {
      f_[pred[i]] = f_[pred[i]] + mbc.parentToSon[i].transMul(f_[i]);
    }
  }
}

void InverseStatics::computeTorqueJacobian(const MultiBody& mb,
                                           MultiBodyConfig& mbc)
{
  const std::vector<Body>& bodies = mb.bodies();
  const std::vector<int>& pred = mb.predecessors();
  const std::vector<Joint>& joints = mb.joints();

  sva::MotionVecd a_0(Eigen::Vector3d::Zero(), mbc.gravity);

  Eigen::Matrix<double, 6, 6> M;

  for (std::size_t i = 0; i < bodies.size(); ++i)
  {
    M.setZero();
    Jacobian jacW(mb, static_cast<int>(i));
    const Eigen::MatrixXd& jac = jacW.jacobian(mb, mbc);
    Eigen::MatrixXd fullJac = Eigen::MatrixXd::Zero(6, mb.nrDof());
    // Complete the previously computed jacobian to a full jacobian
    jacW.fullJacobian(mb, jac, fullJac);

    mbc.bodyAccB[i] = mbc.bodyPosW[i] * a_0;

    Eigen::Matrix3d RW = mbc.bodyPosW[i].rotation();
    Eigen::Matrix3d mRtW =
        -mbc.bodyPosW[i].rotation() *
        Eigen::vector3ToCrossMatrix(mbc.bodyPosW[i].translation());
    Eigen::Vector3d fC = mbc.force[i].couple();
    Eigen::Vector3d fF = mbc.force[i].force();
    Eigen::Matrix3d hatFC, hatFF, hatAF, hatmRtWfF;
    Eigen::Vector3d mRtWfF = mRtW * fF;
    hatFF = Eigen::vector3ToCrossMatrix(fF);
    hatFC = Eigen::vector3ToCrossMatrix(fC);
    hatmRtWfF = Eigen::vector3ToCrossMatrix(mRtWfF);
    hatAF = Eigen::vector3ToCrossMatrix(a_0.linear());

    f_[i] = bodies[i].inertia() * mbc.bodyAccB[i] -
            mbc.bodyPosW[i].dualMul(mbc.force[i]);

    M.block(0, 0, 3, 3) =
        bodies[i].inertia().matrix().block(0, 3, 3, 3) * RW * hatAF;
    M.block(3, 0, 3, 3) =
        bodies[i].inertia().matrix().block(3, 3, 3, 3) * RW * hatAF;
    M.block(0, 0, 3, 3) -= RW * hatFC + hatmRtWfF;
    M.block(3, 0, 3, 3) -= RW * hatFF;
    M.block(0, 3, 3, 3) -= RW * hatFF;

    df_[i] = M * fullJac;
  }

  for (int i = static_cast<int>(bodies.size()) - 1; i >= 0; --i)
  {
    for (int j = 0; j < joints[i].dof(); ++j)
    {
      jointTorqueDiff_[i][j] =
          mbc.motionSubspace[i].col(j).transpose() * df_[i];
    }

    if (pred[i] != -1)
    {
      M.setZero();
      f_[pred[i]] = f_[pred[i]] + mbc.parentToSon[i].transMul(f_[i]);
      M.block(0, 0, 3, 3) = mbc.parentToSon[i].rotation().transpose();
      M.block(0, 3, 3, 3) =
          Eigen::vector3ToCrossMatrix(mbc.parentToSon[i].translation()) *
          mbc.parentToSon[i].rotation().transpose();
      M.block(3, 3, 3, 3) = mbc.parentToSon[i].rotation().transpose();
      df_[pred[i]] = df_[pred[i]] + M * df_[i];

      for (int j = 0; j < joints[i].dof(); ++j)
      {
        Eigen::Vector3d omega = mbc.motionSubspace[i].col(j).head(3);
        Eigen::Vector3d v = mbc.motionSubspace[i].col(j).tail(3);
        Eigen::Matrix3d R = mbc.parentToSon[i].rotation();
        Eigen::Vector3d t = mbc.parentToSon[i].translation();
        Eigen::VectorXd res(6);
        res.head(3) = -omega.cross(R.transpose() * f_[i].couple());
        res.head(3) += v.cross(R.transpose() * f_[i].force());
        res.head(3) -= t.cross(omega.cross(R.transpose() * f_[i].force()));
        res.tail(3) = -omega.cross(R.transpose() * f_[i].force());
        df_[pred[i]].col(mb.jointPosInDof(i) + j) -= res;
      }
    }
  }
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
