#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
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
    /* Left hand fingers */
	filteredLinks.push_back("LFinger21_link");
	filteredLinks.push_back("LFinger22_link");
	filteredLinks.push_back("LFinger23_link");
	filteredLinks.push_back("LFinger11_link");
	filteredLinks.push_back("LFinger12_link");
	filteredLinks.push_back("LFinger13_link");
	filteredLinks.push_back("LFinger31_link");
	filteredLinks.push_back("LFinger32_link");
	filteredLinks.push_back("LFinger33_link");
	filteredLinks.push_back("LFinger41_link");
	filteredLinks.push_back("LFinger42_link");
	filteredLinks.push_back("LFinger43_link");
	filteredLinks.push_back("LThumb1_link");
	filteredLinks.push_back("LThumb2_link");

	/* Right hand fingers */
	filteredLinks.push_back("RFinger21_link");
	filteredLinks.push_back("RFinger22_link");
	filteredLinks.push_back("RFinger23_link");
	filteredLinks.push_back("RFinger11_link");
	filteredLinks.push_back("RFinger12_link");
	filteredLinks.push_back("RFinger13_link");
	filteredLinks.push_back("RFinger31_link");
	filteredLinks.push_back("RFinger32_link");
	filteredLinks.push_back("RFinger33_link");
	filteredLinks.push_back("RFinger41_link");
	filteredLinks.push_back("RFinger42_link");
	filteredLinks.push_back("RFinger43_link");
	filteredLinks.push_back("RThumb1_link");
	filteredLinks.push_back("RThumb2_link");

	/* Wheels */
	filteredLinks.push_back("WheelFL_link");
	filteredLinks.push_back("WheelB_link");
	filteredLinks.push_back("WheelFR_link");

    std::ifstream ifs(urdf_file);
    if(ifs.is_open())
    {
      std::stringstream urdf;
      urdf << ifs.rdbuf();
      mc_rbdyn_urdf::URDFParserResult res = mc_rbdyn_urdf::rbdyn_from_urdf(urdf.str(), true, filteredLinks, true, "base_footprint");
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
  bool closed_loop = false;	

  rbd::MultiBody mb;
  rbd::MultiBodyConfig mbc;
  rbd::MultiBodyGraph mbg;

  std::string urdf_file = source_dir + "/models/pepper.urdf";
  std::tie(mb, mbc, mbg) = makeArm(urdf_file);
  
  std::string input_log_file_name = "mc-control-PepperMove-closed-loop-contacts-Joints.txt";
  std::string log_file = source_dir + "/logs/"+input_log_file_name;
  std::ifstream log_file_strm (log_file);

  std::ofstream output_log_file;
  std::string output_log_file_path = source_dir + "/logs/output_logs/TD_"+input_log_file_name;
  output_log_file.open(output_log_file_path); 

   // start reading log file
  if (log_file_strm.is_open()){
    int line_cnt = 0;
    std::string log_line;
    // get new log line
    while ( getline (log_file_strm, log_line) ){
    	if(line_cnt!=0){ // omit header
    		std::stringstream ss(log_line);
    		int cnt = 0;
    		std::vector<double> q_vec;
			std::vector<double> dq_vec;
			std::vector<double> ddq_vec;
    		// Get line content as vector of double
    		while(ss.good()){
    			std::string substr;
			    getline( ss, substr, ';' );
			    double subnum = std::stod(substr);
			    if(closed_loop){
				    if(cnt>=20 && cnt<=36){
				    	ddq_vec.push_back(subnum);
				    }else if(cnt>=37 && cnt<=53){
				    	dq_vec.push_back(subnum);
				    }else if(cnt>=81 && cnt<=97){
				    	q_vec.push_back(subnum);
				    }
				}else{
					if(cnt>=3 && cnt<=19){
				    	ddq_vec.push_back(subnum);
				    }else if(cnt>=20 && cnt<=36){
				    	dq_vec.push_back(subnum);
				    }else if(cnt>=64 && cnt<=80){
				    	q_vec.push_back(subnum);
				    }
				}
			    cnt++;
    		}

    		Eigen::VectorXd q(mb.nrDof());
			Eigen::VectorXd dq(mb.nrDof());
			Eigen::VectorXd ddq(mb.nrDof());

			for (int i = 0; i < mb.nrDof(); ++i)
			{
				q[i] = q_vec[i];
				dq[i] = dq_vec[i];
				ddq[i] = ddq_vec[i];
			}

    		/*std::cout << q.transpose() << std::endl;
    		std::cout << dq.transpose() << std::endl;
    		std::cout << ddq.transpose() << std::endl;*/

    		mbc.q = rbd::vectorToParam(mb, q);
  			rbd::forwardKinematics(mb, mbc);
  			mbc.alpha = rbd::vectorToParam(mb, dq);
  			mbc.alphaD = rbd::vectorToParam(mb, ddq);
  			rbd::forwardVelocity(mb, mbc);

  			rbd::ForwardDynamics fd(mb);
			fd.computeH(mb, mbc);
			fd.computeC(mb, mbc);

			Eigen::VectorXd tau_d = fd.H()*ddq + fd.C();

			//std::cout << tau_d.transpose() << std::endl;

			// Write results to output log
			output_log_file << log_line << ";";
			for(int i=0;i<tau_d.rows();++i){ 
				if(i!=tau_d.rows()-1){
				  output_log_file << tau_d[i] << ";";
				}else{
				  output_log_file << tau_d[i] << std::endl;
				}
			}
			

    	}else{
    		std::stringstream ss(log_line);
    		int cnt = 0;
    		std::vector<std::string> q_svec;
			std::vector<std::string> dq_svec;
			std::vector<std::string> ddq_svec;
    		// Get line content as vector of double
    		while(ss.good()){
    			std::string substr;
			    getline( ss, substr, ';' );
			    if(closed_loop){
				    if(cnt>=20 && cnt<=36){
				    	ddq_svec.push_back(substr);
				    }else if(cnt>=37 && cnt<=53){
				    	dq_svec.push_back(substr);
				    }else if(cnt>=81 && cnt<=97){
				    	q_svec.push_back(substr);
				    }
				}else{
					if(cnt>=3 && cnt<=19){
				    	ddq_svec.push_back(substr);
				    }else if(cnt>=20 && cnt<=36){
				    	dq_svec.push_back(substr);
				    }else if(cnt>=64 && cnt<=80){
				    	q_svec.push_back(substr);
				    }
				}
			    cnt++;
    		}

    		for (unsigned int i = 0; i < q_svec.size(); ++i){std::cout << q_svec[i] << " ";}
    		std::cout << " " << std::endl;
    		for (unsigned int i = 0; i < dq_svec.size(); ++i){std::cout << dq_svec[i] << " ";}
    		std::cout << " " << std::endl;
    		for (unsigned int i = 0; i < ddq_svec.size(); ++i){std::cout << ddq_svec[i] << " ";}
    		std::cout << " " << std::endl;

    		// write header
    		output_log_file << log_line << ";";
    		for (int i = 0; i < mb.nrDof(); ++i)
			{
				if(i!=mb.nrDof()-1){
					output_log_file << "tauOut_" << i << ";";
				}else{
					output_log_file << "tauOut_" << i << std::endl;
				}
			}
    	}
    	line_cnt++;
    	/*if(line_cnt==2){
    		break;
    	}*/
    }
  }

  output_log_file.close();
  log_file_strm.close();
  std::cout << "Finished writing log file: " << output_log_file_path << std::endl;

  return 0;
}
