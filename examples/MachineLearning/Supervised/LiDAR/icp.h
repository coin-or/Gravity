//
//  icp.h
//  Gravity
//
//  Created by Smitha on 3/31/22.
//

#ifndef icp_h
#define icp_h
#include "Lidar_utils.h"
#ifdef USE_EIGEN3
#include <Eigen/Dense>
#include <Eigen/SVD>
#endif
using namespace Eigen;
vector<double> compute_L2_error(const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data_new, const vector<vector<double>>& point_cloud_data);
void icp_new(const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data, double roll_min_i, double roll_max_i, double pitch_min_i, double pitch_max_i, double yaw_min_i, double yaw_max_i, double tx_min_i, double tx_max_i, double ty_min_i, double ty_max_i, double tz_min_i, double tz_max_i, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double tx_min, double tx_max, double ty_min, double ty_max, double tz_min, double tz_max, vector<double>& rpyt){
    rpyt.resize(7,0);
    double det=0;
    int nd=point_cloud_data.size();
    int nm=point_cloud_model.size();
    Eigen::MatrixXd H(3,3), R(3,3), U(3,3), V(3,3);
    vector<double> mu_m(3,0);
    vector<int> matching(nd);
    double error_old=nd*12+10;
    double error_new=nd*12;
    double error=0;
    
    double roll=0.5*(roll_min_i+roll_max_i);
    double pitch=0.5*(pitch_min_i+pitch_max_i);
    double yaw=0.5*(yaw_min_i+yaw_max_i);
    double tx=0.5*(tx_min_i+tx_max_i);
    double ty=0.5*(ty_min_i+ty_max_i);
    double tz=0.5*(tz_min_i+tz_max_i);
    R(0,0)= cos(yaw)*cos(roll);
    R(0,1)= cos(yaw)*sin(roll)*sin(pitch) - sin(yaw)*cos(pitch);
    R(0,2)= cos(yaw)*sin(roll)*cos(pitch) + sin(yaw)*sin(pitch);
    R(1,0)= sin(yaw)*cos(roll);
    R(1,1)= sin(yaw)*sin(roll)*sin(pitch) + cos(yaw)*cos(pitch);
    R(1,2)= sin(yaw)*sin(roll)*cos(pitch) - cos(yaw)*sin(pitch);
    R(2,0)= sin(-1*roll);
    R(2,1)= cos(roll)*sin(pitch);
    R(2,2)= cos(roll)*cos(pitch);
    
    
    Eigen::Vector3d qrt,qd,qm, t;
    
    t<<tx,ty,tz;
    int iter=0;
    while(error_old-error_new>=1e-4 && iter<=1000){
        for(auto i=0;i<3;i++){
            for(auto j=0;j<3;j++){
                H(i,j)=0.0;
            }
        }
        iter++;
        error_old=error_new;
        error_new=0.0;
        mu_m[0]=0;mu_m[1]=0;mu_m[2]=0;
        for(auto i=0;i<nd;i++){
            qd<<point_cloud_data[i][0],point_cloud_data[i][1],point_cloud_data[i][2];
            qrt=R*qd+t;
            double dist_min=12.0; int j_min=0;
            for(auto j=0;j<point_cloud_model.size();j++){
                double dist=std::pow(qrt[0] - point_cloud_model.at(j).at(0),2) + std::pow(qrt[1] - point_cloud_model.at(j).at(1),2) + std::pow(qrt[2] - point_cloud_model.at(j).at(2),2);
                if(dist<=dist_min){
                    dist_min = dist;
                    j_min=j;
                }
            }
            matching[i]=j_min;
            error_new+=dist_min;
            mu_m[0]+=point_cloud_model.at(j_min).at(0);
            mu_m[1]+=point_cloud_model.at(j_min).at(1);
            mu_m[2]+=point_cloud_model.at(j_min).at(2);
        }
        mu_m[0]/=nd;mu_m[1]/=nd;mu_m[2]/=nd;
    int j=0;
    for(auto i=0;i<nd;i++){
        qd<<point_cloud_data[i][0], point_cloud_data[i][1], point_cloud_data[i][2];
        j=matching[i];
        qm<<point_cloud_model[j][0]-mu_m[0], point_cloud_model[j][1]-mu_m[1],point_cloud_model[j][2]-mu_m[2];
        H=H+qd*qm.transpose();
    }
    JacobiSVD<MatrixXd> svd( H, ComputeFullV | ComputeFullU );
     U=svd.matrixU();
     V=svd.matrixV();

     R=V*U.transpose();
     t<<mu_m[0], mu_m[1], mu_m[2];

     det=R(0,0)*(R(1,1)*R(2,2)-R(2,1)*R(1,2))-R(0,1)*(R(1,0)*R(2,2)-R(2,0)*R(1,2))+R(0,2)*(R(1,0)*R(2,1)-R(2,0)*R(1,1));

    DebugOff("det "<<det<<endl);
    DebugOff(R(0,0)<<" "<<R(0,1)<<" "<<R(0,2)<<endl);
    DebugOff(R(1,0)<<" "<<R(1,1)<<" "<<R(1,2)<<endl);
    DebugOff(R(2,0)<<" "<<R(2,1)<<" "<<R(2,2)<<endl);
    if((det)<=0.999 || (det)>=1.001){
        DebugOn("failed "<<det<<endl);
        break;
    }
        
    }
    
    DebugOff("err "<<error_new<<" "<<iter<<endl);
   
    pitch= atan2(R(2,1), R(2,2));
    roll = atan2(-R(2,0), std::sqrt(R(2,1)*R(2,1)+R(2,2)*R(2,2)));
    yaw = atan2(R(1,0),R(0,0));
    
 
    
    if(roll_min<=roll && roll<=roll_max && pitch_min<=pitch && pitch<=pitch_max && yaw_min<=yaw && yaw<=yaw_max && tx_min<=tx && tx<=tx_max && ty_min<=ty && ty<=ty_max && tz_min<=tz && tz<=tz_max && det>=0.999 && det<=1.001){
        error=error_new;
        }
    else{
        error=error_new+100;
    }
        
    
        
        rpyt[0]=error;
        rpyt[1]=roll;
        rpyt[2]=pitch;
        rpyt[3]=yaw;
        rpyt[4]=mu_m[0];
        rpyt[5]=mu_m[1];
        rpyt[6]=mu_m[2];
}
vector<double> icp(const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data, double roll_min_i, double roll_max_i, double pitch_min_i, double pitch_max_i, double yaw_min_i, double yaw_max_i, double tx_min_i, double tx_max_i, double ty_min_i, double ty_max_i, double tz_min_i, double tz_max_i, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double tx_min, double tx_max, double ty_min, double ty_max, double tz_min, double tz_max)
{

    int nd=point_cloud_data.size();
    int nm=point_cloud_model.size();

    vector<double> rpyt(7,0);

    double error_old=nd*12+10;
    double error_new=nd*12;
    
    rpyt[1]=0.5*(roll_min_i+roll_max_i);
    rpyt[2]=0.5*(pitch_min_i+pitch_max_i);
    rpyt[3]=0.5*(yaw_min_i+yaw_max_i);
    rpyt[4]=0.5*(tx_min_i+tx_max_i);
    rpyt[5]=0.5*(ty_min_i+ty_max_i);
    rpyt[6]=0.5*(tz_min_i+tz_max_i);
    int iter=0;
    while(error_old-error_new>=1e-6 && iter<=1000){
        iter++;
        error_old=error_new;
        vector<vector<double>> point_cloud_data_new=point_cloud_data;
        apply_rot_trans(rpyt[1], rpyt[2], rpyt[3], rpyt[4], rpyt[5], rpyt[6], point_cloud_data_new);
        rpyt=compute_L2_error(point_cloud_model, point_cloud_data_new, point_cloud_data);
        error_new=rpyt[0];
}
    return rpyt;
}
    
    
vector<double> compute_L2_error(const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data_new, const vector<vector<double>>& point_cloud_data){
    vector<double> rpyt(7,0);
    Eigen::MatrixXd H(3,3);
    for(auto i=0;i<3;i++){
        for(auto j=0;j<3;j++){
            H(i,j)=0.0;
        }
    }
    int nd=point_cloud_data_new.size();
    vector<double> mu_m(3,0);
    vector<int> matching(nd);
    
    double error=0;
    Eigen::Vector3d qd,qm;
    for(auto i=0;i<nd;i++){
        double dist_min=12.0; int j_min=0;
        for(auto j=0;j<point_cloud_model.size();j++){
            double dist=std::pow(point_cloud_data_new.at(i).at(0) - point_cloud_model.at(j).at(0),2) + std::pow(point_cloud_data_new.at(i).at(1) - point_cloud_model.at(j).at(1),2) + std::pow(point_cloud_data_new.at(i).at(2) - point_cloud_model.at(j).at(2),2);
            if(dist<=dist_min){
                dist_min = dist;
                j_min=j;
            }
        }
        matching[i]=j_min;
        error+=dist_min;
        mu_m[0]+=point_cloud_model.at(j_min).at(0);
        mu_m[1]+=point_cloud_model.at(j_min).at(1);
        mu_m[2]+=point_cloud_model.at(j_min).at(2);
    }
    mu_m[0]/=nd;mu_m[1]/=nd;mu_m[2]/=nd;

    int j=0;
    for(auto i=0;i<nd;i++){
        qd<<point_cloud_data[i][0], point_cloud_data[i][1],point_cloud_data[i][2];
        j=matching[i];
        qm<<point_cloud_model[j][0]-mu_m[0], point_cloud_model[j][1]-mu_m[1],point_cloud_model[j][2]-mu_m[2];
        H=H+qd*qm.transpose();
    }
    JacobiSVD<MatrixXd> svd( H, ComputeFullV | ComputeFullU );
     auto U=svd.matrixU();
     auto V=svd.matrixV();

     auto R=V*U.transpose();

    double det=R(0,0)*(R(1,1)*R(2,2)-R(2,1)*R(1,2))-R(0,1)*(R(1,0)*R(2,2)-R(2,0)*R(1,2))+R(0,2)*(R(1,0)*R(2,1)-R(2,0)*R(1,1));

    DebugOff("det "<<det<<endl);
    DebugOff(R(0,0)<<" "<<R(0,1)<<" "<<R(0,2)<<endl);
    DebugOff(R(1,0)<<" "<<R(1,1)<<" "<<R(1,2)<<endl);
    DebugOff(R(2,0)<<" "<<R(2,1)<<" "<<R(2,2)<<endl);
    if(det<=0.999 || det>=1.001){
        DebugOn("failed "<<det<<endl);
        error=error+1000;
    }
    
    
    double roll, pitch, yaw;
    pitch= atan2(R(2,1), R(2,2));
    roll = atan2(-R(2,0), std::sqrt(R(2,1)*R(2,1)+R(2,2)*R(2,2)));
    yaw = atan2(R(1,0),R(0,0));
    
    rpyt[0]=error;
    rpyt[1]=roll;
    rpyt[2]=pitch;
    rpyt[3]=yaw;
    rpyt[4]=mu_m[0];
    rpyt[5]=mu_m[1];
    rpyt[6]=mu_m[2];
    
    
    return(rpyt);
}

#endif /* icp_h */
