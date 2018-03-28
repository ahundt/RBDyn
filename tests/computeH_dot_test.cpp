#include "XYZarm.h"

#include <RBDyn/FD.h>
#include <RBDyn/FK.h>
#include <RBDyn/FV.h>

#include <iostream>

int main()
{
  srand(time(NULL));

  rbd::MultiBody mb;
  rbd::MultiBodyConfig mbc;
  rbd::MultiBodyGraph mbg;

  std::tie(mb, mbc, mbg) = makeXYZarm();

  Eigen::VectorXd q = Eigen::VectorXd::Random(mb.nrParams());
  mbc.q = rbd::vectorToParam(mb, q);
  rbd::forwardKinematics(mb, mbc);

  rbd::ForwardDynamics fd(mb);

  Eigen::VectorXd qd = Eigen::VectorXd::Random(mb.nrDof());
  mbc.alpha = rbd::vectorToParam(mb, qd);
  rbd::forwardVelocity(mb, mbc);


  fd.computeH(mb, mbc);
  Eigen::MatrixXd m1 = fd.H();

  double dt = 1e-6;
  q += dt*qd;
  mbc.q = rbd::vectorToParam(mb, q);

  rbd::forwardKinematics(mb, mbc);
  rbd::forwardVelocity(mb, mbc);

  // Computes both H and H_dot
  fd.computeH(mb, mbc);

  Eigen::MatrixXd m2 = fd.H();

  // Numerical differentiation (this is also equal to C^T+C - checked)
  Eigen::MatrixXd H_dot_diff = (m2 - m1)/dt;
  std::cout << "H_dot_diff:" << std::endl;
  std::cout << H_dot_diff << std::endl;
  std::cout << "" << std::endl;

  // True derivative
  std::cout << "H_dot:" << std::endl;
  std::cout << fd.H_dot() << std::endl;
  std::cout << "" << std::endl;



  std::cout << "Error: " << (H_dot_diff-fd.H_dot()).norm() << std::endl;

  return 0;
}

