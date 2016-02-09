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

#pragma once

// includes
// std
#include <vector>

// SpaceVecAlg
#include <SpaceVecAlg/SpaceVecAlg>

namespace rbd
{
class MultiBody;
struct MultiBodyConfig;

struct ForceWithApplicationPoint
{
  public:
    sva::ForceVecd force;
    Eigen::Vector3d pos;
}

/**
	* Inverse Dynamics algorithm.
	*/
class InverseStatics
{
public:
	InverseStatics()
	{}
	/// @param mb MultiBody associated with this algorithm.
	InverseStatics(const MultiBody& mb);

  /**
    * Compute the inverse statics.
    * @param mb MultiBody used has model.
    * @param mbc Uses force, parentToSon, bodyPosW, motionSubspace
    * and gravity.
    * Fills bodyAccB and jointTorque.
    */
        void inverseStatics(const MultiBody& mb, MultiBodyConfig& mbc);

  /**
    * NOTE: THIS CANNOT WORK PROPERLY. If a body is equipped with a force defined on a frame attached to a body,
    * that force won't follow the body during the FD computation
		* Compute the derivatives of the torques calculated by the inverse statics
    * w.r.t. q and forces
    * WARNING: This computes the derivative of the torques w.r.t some fictitious forces applied on each body at the point (0,0,0) of the reference frame.
		* Compute the derivatives of the torques calculated by the inverse statics
    * w.r.t. q and forces.
		* @param mb MultiBody used has model.
		* @param mbc Uses force, parentToSon, bodyPosW, parentToSon, motionSubspace
    * and gravity.
		* Fills jointTorqueJacQ and jointTorqueJacF.
		*/
	void computeTorqueJacobianJoint(const MultiBody& mb, MultiBodyConfig& mbc);

	/**
    * NOTE: THIS CANNOT WORK PROPERLY. If a body is equipped with a force defined on a frame attached to a body,
    * that force won't follow the body during the FD computation
		* Compute the derivatives of the torques calculated by the inverse statics
    * w.r.t. q and forces using Finite Differences.
    * WARNING: This computes the derivative of the torques w.r.t some fictitious forces applied on each body at the point (0,0,0) of the reference frame.
		* @param mb MultiBody used has model.
		* @param mbc Uses force, parentToSon, bodyPosW, parentToSon, motionSubspace
    * and gravity.
		* Fills jointTorqueJacQ and jointTorqueJacF.
		*/
	//void computeTorqueJacobianFD(const MultiBody& mb, const MultiBodyConfig& mbc, double delta = 1e-8);
	//void computeTorqueJacobianJointFD(const MultiBody& mb, const MultiBodyConfig& mbc, double delta = 1e-8);
	//void computeTorqueJacobianForceFD(const MultiBody& mb, const MultiBodyConfig& mbc, double delta = 1e-8);

	// safe version for python binding

	/** safe version of @see inverseDynamics.
		* @throw std::domain_error If mb don't match mbc.
		*/
	void sInverseStatics(const MultiBody& mb, MultiBodyConfig& mbc);

	/**
		* @brief Get the internal forces.
		* @return vector of forces transmitted from body λ(i) to body i across
		* joint i.
		*/
	const std::vector<sva::ForceVecd>& f() const{return f_;};

  const std::vector<std::vector<Eigen::VectorXd>>& jointTorqueDiff() const
  {
    return jointTorqueDiff_;
  };

  /// @brief Vector with all the individual forces applied on each body
  std::vector<std::vector<ForceWithApplicationPoint>> allForces_;

private:
	/// @brief Internal forces.
	/// f_ is the vector of forces transmitted from body λ(i) to body i across
	/// joint i.
	std::vector<sva::ForceVecd> f_;
  std::vector<Eigen::MatrixXd> df_;
  std::vector<std::vector<Eigen::VectorXd>> jointTorqueDiff_;

};

} // namespace rbd
