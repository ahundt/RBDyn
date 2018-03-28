#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>

#include <RBDyn/FD.h>
#include <RBDyn/FV.h>
#include <RBDyn/FK.h>
#include <RBDyn/MultiBody.h>
#include <RBDyn/MultiBodyConfig.h>
#include <mc_rbdyn/Robots.h>
#include <mc_rbdyn/RobotLoader.h>

#include <mc_rbdyn_urdf/urdf.h>

#include <Coriolis/Coriolis.h>

#include <math.h>

static std::string source_dir = CMAKE_SOURCE_DIR;

/// @return multibody from URDF
std::tuple<rbd::MultiBody, rbd::MultiBodyConfig, rbd::MultiBodyGraph>
makeArm(std::string urdf_file)
{
    using namespace rbd;
    MultiBody mb;
    MultiBodyConfig mbc;
    MultiBodyGraph mbg;
    std::vector<std::string> filteredLinks;
    std::ifstream ifs(urdf_file);
    if(ifs.is_open())
    {
      std::stringstream urdf;
      urdf << ifs.rdbuf();
      mc_rbdyn_urdf::URDFParserResult res = mc_rbdyn_urdf::rbdyn_from_urdf(urdf.str(), true, filteredLinks, true, "base_link");
      mb = res.mb;
      mbc = res.mbc;
      mbc.gravity = Eigen::Vector3d(0., 0., 9.81);
      mbg = res.mbg;
    }
    else
    {
      std::cerr << "Failed to open model " << urdf_file << std::endl;
    }
    return std::make_tuple(mb, mbc, mbg);
}


int main()
{
  rbd::MultiBody mb;
  rbd::MultiBodyConfig mbc;
  rbd::MultiBodyGraph mbg;

  std::string urdf_file = source_dir + "/models/3_link_robot.urdf";
  std::tie(mb, mbc, mbg) = makeArm(urdf_file);
  
  Eigen::Vector3d q;
  q<< M_PI/4.0, M_PI/3.0, M_PI/8.0;
  mbc.q = rbd::vectorToParam(mb, q);
  rbd::forwardKinematics(mb, mbc);

  Eigen::Vector3d qd;
  qd<<0.05, 0.08, 0.01;
  mbc.alpha = rbd::vectorToParam(mb, qd);
  rbd::forwardVelocity(mb, mbc);
 
  // Computes both H and H_dot
  rbd::ForwardDynamics fd(mb);
  fd.computeH(mb, mbc);
  fd.computeC(mb, mbc);

  std::cout << "H:" << std::endl;
  std::cout << fd.H() << std::endl;
  std::cout << "" << std::endl;

  // True derivative
  std::cout << "H_dot:" << std::endl;
  std::cout << fd.H_dot() << std::endl;
  std::cout << "" << std::endl;

  // Coriolis matrix
  coriolis::Coriolis coriolis(mb);
  Eigen::MatrixXd C = coriolis.coriolis(mb, mbc);
  std::cout << "C:" << std::endl;
  std::cout << C << std::endl;
  std::cout << "" << std::endl;

  std::cout << "C+C^T:" << std::endl;
  std::cout << C+C.transpose() << std::endl;
  std::cout << "" << std::endl;

  std::cout << "C*alpha:" << std::endl;
  std::cout << C*qd << std::endl;
  std::cout << "" << std::endl;

  std::cout << "C*alpha + g:" << std::endl;
  std::cout << fd.C() << std::endl;
  std::cout << "" << std::endl;

  return 0;
}
