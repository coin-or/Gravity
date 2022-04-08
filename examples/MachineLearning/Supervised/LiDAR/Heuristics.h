//
//  Heuristics.h
//  Gravity
//
//  Created by Smitha on 3/25/22.
//

#ifndef Heuristics_h
#define Heuristics_h
#include <cstdlib>
#include "Goicp.h"
#include "icp.h"
//#include "kdutils.h"
#include "gravity/nanoflann.hpp"
using namespace Go_ICP;
void run_ub_parallel_t(const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data, const vector<double>& roll_lb,  const vector<double>& roll_ub,  const vector<double>& pitch_lb,  const vector<double>& pitch_ub,  const vector<double>& yaw_lb, const vector<double>& yaw_ub, string error_type, vector<double>& ub_node);
void run_ub_parallel(const nanoflann::KDTreeSingleIndexAdaptor<
                     nanoflann::L2_Simple_Adaptor<double, PointCloud<double>>,
                     PointCloud<double>, 3 /* dim */>& index, const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data, const vector<vector<double>>& roll_lb,  const vector<vector<double>>& roll_ub,  const vector<vector<double>>& pitch_lb,  const vector<vector<double>>& pitch_ub,  const vector<vector<double>>& yaw_lb, const vector<vector<double>>& yaw_ub, const vector<vector<double>>& tx_lb, const vector<vector<double>>& tx_ub, const vector<vector<double>>& ty_lb, const vector<vector<double>>& ty_ub, const vector<vector<double>>& tz_lb, const vector<vector<double>>& tz_ub, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double tx_min, double tx_max, double ty_min, double ty_max, double tz_min, double tz_max, string error_type, vector<double>& ub_node, int num_threads);
void evaluate_ub_icp(const nanoflann::KDTreeSingleIndexAdaptor<
                     nanoflann::L2_Simple_Adaptor<double, PointCloud<double>>,
                     PointCloud<double>, 3 /* dim */>& index, const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data, const vector<double>& roll_lb,  const vector<double>& roll_ub,  const vector<double>& pitch_lb,  const vector<double>& pitch_ub,  const vector<double>& yaw_lb, const vector<double>& yaw_ub, const vector<double>& tx_lb, const vector<double>& tx_ub, const vector<double>& ty_lb, const vector<double>& ty_ub, const vector<double>& tz_lb, const vector<double>& tz_ub, vector<vector<double>>& res, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double tx_min, double tx_max, double ty_min, double ty_max, double tz_min, double tz_max, string error_type, double best_ub);
void evaluate_upper_bound_mid(double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max,double tx_min, double tx_max, double ty_min, double ty_max, double tz_min, double tz_max, vector<double>& res, const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data, string error_type, double best_ub);
void evaluate_upper_bound_mid_t(double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, vector<double>& res, const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data, string error_type, const vector<double>& ub_node);
vector<double> ub_heuristic_disc(const nanoflann::KDTreeSingleIndexAdaptor<
                                 nanoflann::L2_Simple_Adaptor<double, PointCloud<double>>,
                                 PointCloud<double>, 3 /* dim */>& index, vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double tx_min, double tx_max, double ty_min, double ty_max, double tz_min, double tz_max, vector<double>& best_rot_trans, double& best_ub, std::string error_type, double max_time)
{
#ifdef USE_MPI
    int worker_id, nb_workers;
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
    auto err_size = MPI_Comm_size(MPI_COMM_WORLD, &nb_workers);
#endif
    vector<double> rpy(6, 0.0);
    /* INPUT BOUNDS */
    
    /* INPUT BOUNDS */
    double prep_time_total=0;
        
    
    size_t nb_threads = std::thread::hardware_concurrency();
    DebugOn("threads "<<nb_threads);
    
    vector<vector<double>> roll_lb(nb_threads), roll_ub(nb_threads), pitch_lb(nb_threads), pitch_ub(nb_threads), yaw_lb(nb_threads), yaw_ub(nb_threads), tx_lb(nb_threads), tx_ub(nb_threads), ty_lb(nb_threads), ty_ub(nb_threads), tz_lb(nb_threads), tz_ub(nb_threads);
    
    
    int nd=point_cloud_data.size();
    
    
    vector<double> ub_node(7,0);
    
    
    pair<double,double> roll_bounds_r, pitch_bounds_r, yaw_bounds_r,tx_bounds_r,ty_bounds_r,tz_bounds_r;
    
    roll_bounds_r={roll_min, roll_max};
    pitch_bounds_r={pitch_min, pitch_max};
    yaw_bounds_r={yaw_min, yaw_max};
    tx_bounds_r={tx_min, tx_max};
    ty_bounds_r={ty_min, ty_max};
    tz_bounds_r={tz_min, tz_max};
    DebugOn("Running aGS for "<<max_time<<" seconds"<<endl);
    double ts=get_wall_time();
    best_ub=100;
    int ndisc=6;
    auto error_init=computeL2error_util(point_cloud_model,point_cloud_data, best_ub);
    DebugOn("error_init "<<error_init<<endl);
    ub_node[0]=error_init;
    bool stop=false;
    int count_probs=0;
    while(!stop){
        double rb=(roll_bounds_r.second-roll_bounds_r.first)/ndisc;
        double pb=(pitch_bounds_r.second-pitch_bounds_r.first)/ndisc;
        double yb=(yaw_bounds_r.second-yaw_bounds_r.first)/ndisc;
        double tx_b=(tx_bounds_r.second-tx_bounds_r.first)/ndisc;
        double ty_b=(ty_bounds_r.second-ty_bounds_r.first)/ndisc;
        double tz_b=(tz_bounds_r.second-tz_bounds_r.first)/ndisc;
        if(rb<=1e-9 || pb<=1e-9 ||yb<=1e-9 || tx_b<=1e-9 || ty_b<=1e-9 ||tz_b<=1e-9 ){
            stop=true;
            break;
        }
        count_probs=0;
        for(auto i=0;i<ndisc;i++){
            for(auto j=0;j<ndisc;j++){
                for(auto k=0;k<ndisc;k++){
                    for(auto l=0;l<ndisc;l++){
                        for(auto m=0;m<ndisc;m++){
                            for(auto n=0;n<ndisc;n++){
                                int thread_no=count_probs%nb_threads;
                                roll_lb[thread_no].push_back(roll_bounds_r.first+i*rb);
                                roll_ub[thread_no].push_back(roll_bounds_r.first+(i+1)*rb);
                                pitch_lb[thread_no].push_back(pitch_bounds_r.first+j*pb);
                                pitch_ub[thread_no].push_back(pitch_bounds_r.first+(j+1)*pb);
                                yaw_lb[thread_no].push_back(yaw_bounds_r.first+k*yb);
                                yaw_ub[thread_no].push_back(yaw_bounds_r.first+(k+1)*yb);
                                tx_lb[thread_no].push_back(tx_bounds_r.first+l*tx_b);
                                tx_ub[thread_no].push_back(tx_bounds_r.first+(l+1)*tx_b);
                                ty_lb[thread_no].push_back(ty_bounds_r.first+m*ty_b);
                                ty_ub[thread_no].push_back(ty_bounds_r.first+(m+1)*ty_b);
                                tz_lb[thread_no].push_back(tz_bounds_r.first+n*tz_b);
                                tz_ub[thread_no].push_back(tz_bounds_r.first+(n+1)*tz_b);
                                count_probs++;
                                
                            }
                        }
                    }
                }
            }
        }
        run_ub_parallel(index, point_cloud_model, point_cloud_data,  roll_lb,   roll_ub,   pitch_lb,   pitch_ub,  yaw_lb,  yaw_ub, tx_lb,  tx_ub,  ty_lb,  ty_ub, tz_lb,  tz_ub, roll_min,  roll_max, pitch_min,  pitch_max,  yaw_min,  yaw_max, tx_min,  tx_max, ty_min,ty_max, tz_min,  tz_max,  error_type, ub_node, nb_threads);
        roll_lb.clear();
        roll_ub.clear();
        pitch_lb.clear();
        pitch_ub.clear();
        yaw_lb.clear();
        yaw_ub.clear();
        tx_lb.clear();
        tx_ub.clear();
        ty_lb.clear();
        ty_ub.clear();
        tz_lb.clear();
        tz_ub.clear();
        best_ub=ub_node[0];
        auto roll_rad1=ub_node[1];
        auto pitch_rad1=ub_node[2];
        auto yaw_rad1 = ub_node[3];
        auto tx1=ub_node[4];
        auto ty1=ub_node[5];
        auto tz1 = ub_node[6];
        roll_bounds_r.first=roll_rad1-std::abs(roll_rad1)*0.1;
        roll_bounds_r.second=roll_rad1+std::abs(roll_rad1)*0.1;
        pitch_bounds_r.first=pitch_rad1-std::abs(pitch_rad1)*0.1;
        pitch_bounds_r.second=pitch_rad1+std::abs(pitch_rad1)*0.1;
        yaw_bounds_r.first=yaw_rad1-std::abs(yaw_rad1)*0.1;
        yaw_bounds_r.second=yaw_rad1+std::abs(yaw_rad1)*0.1;
        tx_bounds_r.first=tx1-std::abs(tx1)*0.1;
        tx_bounds_r.second=tx1+std::abs(tx1)*0.1;
        ty_bounds_r.first=ty1-std::abs(ty1)*0.1;
        ty_bounds_r.second=ty1+std::abs(ty1)*0.1;
        tz_bounds_r.first=tz1-std::abs(tz1)*0.1;
        tz_bounds_r.second=tz1+std::abs(tz1)*0.1;
        //if((get_wall_time()-ts)>=max_time){
        stop=true;
        break;
        //}
    }
    best_ub=ub_node[0];
    rpy[0]=ub_node[1];
    rpy[1]=ub_node[2];
    rpy[2]=ub_node[3];
    rpy[3]=ub_node[4];
    rpy[4]=ub_node[5];
    rpy[5]=ub_node[6];
#ifdef USE_MPI
    if(worker_id==0){
#endif
        DebugOn("Final time "<<(get_wall_time()-ts)<<endl);
        DebugOn("\n*************************\n");
        DebugOn("L2 Final "<<best_ub<<endl);
        DebugOn("*************************\n");
        DebugOn("roll rad "<< rpy[0]<<endl);
        DebugOn("pitch rad "<< rpy[1]<<endl);
        DebugOn("yaw rad "<< rpy[2]<<endl);
        DebugOn("tx "<< rpy[3]<<endl);
        DebugOn("ty "<< rpy[4]<<endl);
        DebugOn("tz "<< rpy[5]<<endl);
#ifdef USE_MPI
    }
#endif
    
    return rpy;
    
}
vector<double> ub_heuristic_disc_t(vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, vector<double>& best_rot_trans, double& best_ub, std::string error_type, double max_time=100)
{
#ifdef USE_MPI
    int worker_id, nb_workers;
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
    auto err_size = MPI_Comm_size(MPI_COMM_WORLD, &nb_workers);
#endif
    vector<double> rpy(6, 0.0);
    /* INPUT BOUNDS */
    
    /* INPUT BOUNDS */
    double prep_time_total=0;
    
    double yaw_min = -90*pi/180., yaw_max = 90*pi/180., pitch_min =-90*pi/180.,pitch_max = 90*pi/180.,roll_min =-90*pi/180.,roll_max = 90*pi/180., tx_min=-0.5, tx_max=0.5, ty_min=-0.5, ty_max=0.5,tz_min=-0.5, tz_max=0.5;
    
    
    vector<double> roll_lb, roll_ub, pitch_lb, pitch_ub, yaw_lb, yaw_ub;
    
    
    int nd=point_cloud_data.size();
    
    
    vector<double> ub_node(7,0.0);
    
    size_t nb_threads = std::thread::hardware_concurrency();
    DebugOn("threads "<<nb_threads);
    
    // nb_threads = 1;
    pair<double,double> roll_bounds_r, pitch_bounds_r, yaw_bounds_r,tx_bounds_r,ty_bounds_r,tz_bounds_r;
    
    roll_bounds_r={roll_min, roll_max};
    pitch_bounds_r={pitch_min, pitch_max};
    yaw_bounds_r={yaw_min, yaw_max};
    tx_bounds_r={tx_min, tx_max};
    ty_bounds_r={ty_min, ty_max};
    tz_bounds_r={tz_min, tz_max};
    DebugOn("Running aGS for "<<max_time<<" seconds"<<endl);
    double ts=get_wall_time();
    best_ub=100;
    int ndisc=30;
    auto error_init=computeL2error_util(point_cloud_model,point_cloud_data, best_ub);
    DebugOn("error_init "<<error_init<<endl);
    ub_node[0]=error_init;
    bool stop=false;
    int count_probs=0;
    while(!stop){
        double rb=(roll_bounds_r.second-roll_bounds_r.first)/ndisc;
        double pb=(pitch_bounds_r.second-pitch_bounds_r.first)/ndisc;
        double yb=(yaw_bounds_r.second-yaw_bounds_r.first)/ndisc;
        if(rb<=1e-9 || pb<=1e-9 ||yb<=1e-9 ){
            stop=true;
            break;
        }
        
        for(auto i=0;i<ndisc;i++){
            for(auto j=0;j<ndisc;j++){
                for(auto k=0;k<ndisc;k++){
                    count_probs++;
                    roll_lb.push_back(roll_bounds_r.first+i*rb);
                    roll_ub.push_back(roll_bounds_r.first+(i+1)*rb);
                    pitch_lb.push_back(pitch_bounds_r.first+j*pb);
                    pitch_ub.push_back(pitch_bounds_r.first+(j+1)*pb);
                    yaw_lb.push_back(yaw_bounds_r.first+k*yb);
                    yaw_ub.push_back(yaw_bounds_r.first+(k+1)*yb);
                    if(count_probs>=nb_threads){
                        run_ub_parallel_t(point_cloud_model,  point_cloud_data, roll_lb,  roll_ub,  pitch_lb,   pitch_ub,  yaw_lb,  yaw_ub,error_type, ub_node);
                        best_ub=ub_node[0];
                        count_probs=0;
                        roll_lb.clear();
                        roll_ub.clear();
                        pitch_lb.clear();
                        pitch_ub.clear();
                        yaw_lb.clear();
                        yaw_ub.clear();
                    }
                    
                    if((get_wall_time()-ts)>=max_time){
                        stop=true;
                        break;
                    }
                }
                if((get_wall_time()-ts)>=max_time){
                    stop=true;
                    break;
                }
            }
            if((get_wall_time()-ts)>=max_time){
                stop=true;
                break;
            }
        }
        auto roll_rad1=ub_node[1];
        auto pitch_rad1=ub_node[2];
        auto yaw_rad1 = ub_node[3];
        auto tx1=ub_node[4];
        auto ty1=ub_node[5];
        auto tz1 = ub_node[6];
        roll_bounds_r.first=roll_rad1-std::abs(roll_rad1)*0.1;
        roll_bounds_r.second=roll_rad1+std::abs(roll_rad1)*0.1;
        pitch_bounds_r.first=pitch_rad1-std::abs(pitch_rad1)*0.1;
        pitch_bounds_r.second=pitch_rad1+std::abs(pitch_rad1)*0.1;
        yaw_bounds_r.first=yaw_rad1-std::abs(yaw_rad1)*0.1;
        yaw_bounds_r.second=yaw_rad1+std::abs(yaw_rad1)*0.1;
        DebugOn("finished round at "<<get_wall_time()-ts<<endl);
        if((get_wall_time()-ts)>=max_time){
            stop=true;
            break;
        }
    }
    best_ub=ub_node[0];
    rpy[0]=ub_node[1];
    rpy[1]=ub_node[2];
    rpy[2]=ub_node[3];
    rpy[3]=ub_node[4];
    rpy[4]=ub_node[5];
    rpy[5]=ub_node[6];
#ifdef USE_MPI
    if(worker_id==0){
#endif
        DebugOn("Final time "<<(get_wall_time()-ts)<<endl);
        DebugOn("\n*************************\n");
        DebugOn("L2 Final "<<best_ub<<endl);
        DebugOn("*************************\n");
        DebugOn("roll rad "<< rpy[0]<<endl);
        DebugOn("pitch rad "<< rpy[1]<<endl);
        DebugOn("yaw rad "<< rpy[2]<<endl);
        DebugOn("tx "<< rpy[3]<<endl);
        DebugOn("ty "<< rpy[4]<<endl);
        DebugOn("tz "<< rpy[5]<<endl);
#ifdef USE_MPI
    }
#endif
    
    return rpy;
    
}
void run_ub_parallel(const nanoflann::KDTreeSingleIndexAdaptor<
                     nanoflann::L2_Simple_Adaptor<double, PointCloud<double>>,
                     PointCloud<double>, 3 /* dim */>& index, const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data, const vector<vector<double>>& roll_lb,  const vector<vector<double>>& roll_ub,  const vector<vector<double>>& pitch_lb,  const vector<vector<double>>& pitch_ub,  const vector<vector<double>>& yaw_lb, const vector<vector<double>>& yaw_ub, const vector<vector<double>>& tx_lb, const vector<vector<double>>& tx_ub, const vector<vector<double>>& ty_lb, const vector<vector<double>>& ty_ub, const vector<vector<double>>& tz_lb, const vector<vector<double>>& tz_ub, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max,double tx_min, double tx_max, double ty_min, double ty_max, double tz_min, double tz_max, string error_type, vector<double>& ub_node, int num_threads){
    std::vector<thread> threads;
    
    int nd=point_cloud_data.size();
    int num=roll_lb.size();
    if(num==0){
        DebugOff("in run_parallel(models...), models is empty, returning");
    }
    vector<vector<vector<double>>> vec_ub;
    vec_ub.resize(num_threads);
    
    
    int count=0;
    for (auto i = 0; i < num_threads; i++) {
        vec_ub[i].resize(roll_lb[i].size());
        threads.push_back(thread(&evaluate_ub_icp, ref(index), ref(point_cloud_model), ref(point_cloud_data),ref(roll_lb[i]), ref(roll_ub[i]), ref(pitch_lb[i]), ref(pitch_ub[i]), ref(yaw_lb[i]), ref(yaw_ub[i]), ref(tx_lb[i]), ref(tx_ub[i]), ref(ty_lb[i]), ref(ty_ub[i]), ref(tz_lb[i]), ref(tz_ub[i]), ref(vec_ub[i]), roll_min,  roll_max, pitch_min,  pitch_max,  yaw_min,  yaw_max, tx_min,  tx_max, ty_min,ty_max, tz_min,  tz_max,  error_type, ub_node[0]));
    }
    for(auto &t : threads){
        t.join();
    }
    threads.clear();
    
    for(auto i=0;i<num_threads;i++){
        for(auto j=0;j<vec_ub[i].size();j++){
            if(vec_ub[i][j][0]<=ub_node[0]){
                ub_node=vec_ub[i][j];
            }
        }
    }
}
void evaluate_ub_icp(const nanoflann::KDTreeSingleIndexAdaptor<
                     nanoflann::L2_Simple_Adaptor<double, PointCloud<double>>,
                     PointCloud<double>, 3 /* dim */>& index, const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data, const vector<double>& roll_lb,  const vector<double>& roll_ub,  const vector<double>& pitch_lb,  const vector<double>& pitch_ub,  const vector<double>& yaw_lb, const vector<double>& yaw_ub, const vector<double>& tx_lb, const vector<double>& tx_ub, const vector<double>& ty_lb, const vector<double>& ty_ub, const vector<double>& tz_lb, const vector<double>& tz_ub, vector<vector<double>>& res, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double tx_min, double tx_max, double ty_min, double ty_max, double tz_min, double tz_max, string error_type, double best_ub){
    //auto goicp=initialize_ICP_only(point_cloud_model, point_cloud_data);
    for(auto i=0;i<roll_lb.size();i++){
        res[i].resize(7);
        icp_new(index, point_cloud_model, point_cloud_data,roll_lb[i],roll_ub[i], pitch_lb[i],
                           pitch_ub[i],yaw_lb[i], yaw_ub[i], tx_lb[i], tx_ub[i],  ty_lb[i], ty_ub[i],  tz_lb[i], tz_ub[i], roll_min, roll_max,  pitch_min,  pitch_max,  yaw_min,  yaw_max,tx_min, tx_max, ty_min, ty_max, tz_min, tz_max, res[i]);
        //compute_upper_boundICP(goicp, roll_lb[i], roll_ub[i], pitch_lb[i], pitch_ub[i], yaw_lb[i], yaw_ub[i], tx_lb[i], tx_ub[i], ty_lb[i], ty_ub[i], tz_lb[i], tz_ub[i], roll_min, roll_max, pitch_min, pitch_max, yaw_min, yaw_max,  tx_min, tx_max, ty_min, ty_max, tz_min, tz_max, res[i], best_ub, point_cloud_model, point_cloud_data);
    }
    
}
void evaluate_upper_bound_mid(double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double tx_min, double tx_max, double ty_min, double ty_max, double tz_min, double tz_max, vector<double>& res, const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data, string error_type, double best_ub){
    res.resize(7);
    double roll=(roll_min+roll_max)*0.5;
    double pitch=(pitch_min+pitch_max)*0.5;
    double yaw=(yaw_min+yaw_max)*0.5;
    double tx=(tx_min+tx_max)*0.5;
    double ty=(ty_min+ty_max)*0.5;
    double tz=(tz_min+tz_max)*0.5;
    
    vector<double> rot(12);
    
    rot[0]=cos(roll)*cos(yaw);
    rot[1]=cos(yaw)*sin(roll)*sin(pitch) - sin(yaw)*cos(pitch);
    rot[2]=cos(yaw)*sin(roll)*cos(pitch) + sin(yaw)*sin(pitch);
    rot[3]=sin(yaw)*cos(roll);
    rot[4]=sin(yaw)*sin(roll)*sin(pitch) + cos(yaw)*cos(pitch);
    rot[5]=sin(yaw)*sin(roll)*cos(pitch) - cos(yaw)*sin(pitch);
    rot[6]=sin(-1*roll);
    rot[7]=cos(roll)*sin(pitch);
    rot[8]=cos(roll)*cos(pitch);
    rot[9]=tx;
    rot[10]=ty;
    rot[11]=tz;
    double error;
    
    vector<vector<double>> point_cloud_data_copy=point_cloud_data;
    
    
    
    apply_rot_trans_util(rot, point_cloud_data_copy);
    error = computeL2error_util(point_cloud_model,point_cloud_data_copy,best_ub);
    DebugOn("error"<<error<<endl);
    
    res[0]=error;
    res[1]=roll;
    res[2]=pitch;
    res[3]=yaw;
    res[4]=tx;
    res[5]=ty;
    res[6]=tz;
}
void run_ub_parallel_t(const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data, const vector<double>& roll_lb,  const vector<double>& roll_ub,  const vector<double>& pitch_lb,  const vector<double>& pitch_ub,  const vector<double>& yaw_lb, const vector<double>& yaw_ub, string error_type, vector<double>& ub_node){
    std::vector<thread> threads;
    
    int nd=point_cloud_data.size();
    int num=roll_lb.size();
    if(num==0){
        DebugOff("in run_parallel(models...), models is empty, returning");
    }
    vector<vector<double>> vec_ub;
    vec_ub.resize(num);
    
    for (auto i = 0; i < num; i++) {
        threads.push_back(thread(&evaluate_upper_bound_mid_t, roll_lb[i], roll_ub[i], pitch_lb[i], pitch_ub[i], yaw_lb[i], yaw_ub[i], ref(vec_ub[i]), ref(point_cloud_model), ref(point_cloud_data), error_type, ref(ub_node)));
    }
    for(auto &t : threads){
        t.join();
    }
    threads.clear();
    
    for(auto i=0;i<num;i++){
        if(vec_ub[i][0]<=ub_node[0]){
            ub_node=vec_ub[i];
        }
    }
}
void evaluate_upper_bound_mid_t(double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, vector<double>& res, const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data, string error_type, const vector<double>& ub_node){
    res.resize(7);
    double roll=(roll_min+roll_max)*0.5;
    double pitch=(pitch_min+pitch_max)*0.5;
    double yaw=(yaw_min+yaw_max)*0.5;
    double tx=ub_node[4];
    double ty=ub_node[5];
    double tz=ub_node[6];
    
    vector<double> rot(12);
    
    rot[0]=cos(roll)*cos(yaw);
    rot[1]=cos(yaw)*sin(roll)*sin(pitch) - sin(yaw)*cos(pitch);
    rot[2]=cos(yaw)*sin(roll)*cos(pitch) + sin(yaw)*sin(pitch);
    rot[3]=sin(yaw)*cos(roll);
    rot[4]=sin(yaw)*sin(roll)*sin(pitch) + cos(yaw)*cos(pitch);
    rot[5]=sin(yaw)*sin(roll)*cos(pitch) - cos(yaw)*sin(pitch);
    rot[6]=sin(-1*roll);
    rot[7]=cos(roll)*sin(pitch);
    rot[8]=cos(roll)*cos(pitch);
    rot[9]=tx;
    rot[10]=ty;
    rot[11]=tz;
    double error;
    for(auto i=0;i<10;i++){
        vector<vector<double>> point_cloud_data_copy=point_cloud_data;
        apply_rot_trans_util(rot, point_cloud_data_copy);
        tx=0;ty=0;tz=0;
        error = computeL2error_util(point_cloud_model,point_cloud_data_copy,ub_node[0],tx,ty,tz);
        rot[9]=tx;
        rot[10]=ty;
        rot[11]=tz;
        DebugOff("error"<<error<<endl);
    }
    
    res[0]=error;
    res[1]=roll;
    res[2]=pitch;
    res[3]=yaw;
    res[4]=tx;
    res[5]=ty;
    res[6]=tz;
    
}
//vector<double> ub_heuristic_icp(vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double tx_min, double tx_max, double ty_min, double ty_max, double tz_min, double tz_max, vector<double>& best_rot_trans, double& best_ub, std::string error_type, double max_time=100)
//{
//#ifdef USE_MPI
//    int worker_id, nb_workers;
//    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
//    auto err_size = MPI_Comm_size(MPI_COMM_WORLD, &nb_workers);
//#endif
//    vector<double> rpy(6, 0.0);
//    /* INPUT BOUNDS */
//
//    /* INPUT BOUNDS */
//    double prep_time_total=0;
//
//    size_t nb_threads = std::thread::hardware_concurrency();
//    //nb_threads=1;
//    DebugOn("threads "<<nb_threads);
//    std::vector<thread> threads;
//
//    vector<double> roll_lb, roll_ub, pitch_lb, pitch_ub, yaw_lb, yaw_ub, tx_lb, tx_ub, ty_lb, ty_ub, tz_lb, tz_ub;
//    vector<vector<double>> vec_ub(nb_threads);
//
//    int nd=point_cloud_data.size();
//
//
//    vector<double> ub_node(7,0);
//
//
//    pair<double,double> roll_bounds_r, pitch_bounds_r, yaw_bounds_r,tx_bounds_r,ty_bounds_r,tz_bounds_r;
//
//    roll_bounds_r={roll_min, roll_max};
//    pitch_bounds_r={pitch_min, pitch_max};
//    yaw_bounds_r={yaw_min, yaw_max};
//    tx_bounds_r={tx_min, tx_max};
//    ty_bounds_r={ty_min, ty_max};
//    tz_bounds_r={tz_min, tz_max};
//    DebugOn("Running aGS for "<<max_time<<" seconds"<<endl);
//    double ts=get_wall_time();
//    best_ub=100;
//    int ndisc=6;
//    auto error_init=computeL2error_util(point_cloud_model,point_cloud_data, best_ub);
//    DebugOn("error_init "<<error_init<<endl);
//    ub_node[0]=error_init;
//
//    double rb=(roll_bounds_r.second-roll_bounds_r.first)/ndisc;
//    double pb=(pitch_bounds_r.second-pitch_bounds_r.first)/ndisc;
//    double yb=(yaw_bounds_r.second-yaw_bounds_r.first)/ndisc;
//    double tx_b=(tx_bounds_r.second-tx_bounds_r.first)/ndisc;
//    double ty_b=(ty_bounds_r.second-ty_bounds_r.first)/ndisc;
//    double tz_b=(tz_bounds_r.second-tz_bounds_r.first)/ndisc;
//    int count_probs=0;
//    int count_total_probs=0;
//    for(auto i=0;i<ndisc;i++){
//        for(auto j=0;j<ndisc;j++){
//            for(auto k=0;k<ndisc;k++){
//                for(auto l=0;l<ndisc;l++){
//                    for(auto m=0;m<ndisc;m++){
//                        for(auto n=0;n<ndisc;n++){
//                            count_probs++;
//                            count_total_probs++;
//                            roll_lb.push_back(roll_bounds_r.first+i*rb);
//                            roll_ub.push_back(roll_bounds_r.first+(i+1)*rb);
//                            pitch_lb.push_back(pitch_bounds_r.first+j*pb);
//                            pitch_ub.push_back(pitch_bounds_r.first+(j+1)*pb);
//                            yaw_lb.push_back(yaw_bounds_r.first+k*yb);
//                            yaw_ub.push_back(yaw_bounds_r.first+(k+1)*yb);
//                            tx_lb.push_back(tx_bounds_r.first+l*tx_b);
//                            tx_ub.push_back(tx_bounds_r.first+(l+1)*tx_b);
//                            ty_lb.push_back(ty_bounds_r.first+m*ty_b);
//                            ty_ub.push_back(ty_bounds_r.first+(m+1)*ty_b);
//                            tz_lb.push_back(tz_bounds_r.first+n*tz_b);
//                            tz_ub.push_back(tz_bounds_r.first+(n+1)*tz_b);
//
//                            if(count_probs>=nb_threads || count_total_probs>=pow(ndisc,6)){
//                                for(auto p=0;p<count_probs;p++){
//                                    vec_ub[p].resize(7);
//                                    threads.push_back(thread(&icp_new, ref(point_cloud_model), ref(point_cloud_data), roll_lb[p], roll_ub[p], pitch_lb[p], pitch_ub[p], yaw_lb[p], yaw_ub[p], tx_lb[p], tx_ub[p], ty_lb[p], ty_ub[p], tz_lb[p], tz_ub[p],  roll_min,  roll_max, pitch_min, pitch_max, yaw_min, yaw_max, tx_min,  tx_max,  ty_min, ty_max,  tz_min,  tz_max, ref(vec_ub[p])));
//                                }
//                                for(auto &t : threads){
//                                    t.join();
//                                }
//                                threads.clear();
//
//                                for(auto p=0;p<count_probs;p++){
//                                    if(vec_ub[p][0]<=ub_node[0]){
//                                        ub_node=vec_ub[p];
//                                    }
//                                }
//                                count_probs=0;
//                                roll_lb.clear(); roll_ub.clear(); pitch_lb.clear(); pitch_ub.clear();yaw_lb.clear(); yaw_ub.clear(); tx_lb.clear(); tx_ub.clear(); ty_lb.clear(); ty_ub.clear(); tz_lb.clear(); tz_ub.clear();
//
//                            }
//                        }
//                    }
//                }
//            }
//        }
//    }
//    best_ub=ub_node[0];
//    auto roll_rad1=ub_node[1];
//    auto pitch_rad1=ub_node[2];
//    auto yaw_rad1 = ub_node[3];
//    auto tx1=ub_node[4];
//    auto ty1=ub_node[5];
//    auto tz1 = ub_node[6];
//    //if((get_wall_time()-ts)>=max_time){
//
//    best_ub=ub_node[0];
//    rpy[0]=ub_node[1];
//    rpy[1]=ub_node[2];
//    rpy[2]=ub_node[3];
//    rpy[3]=ub_node[4];
//    rpy[4]=ub_node[5];
//    rpy[5]=ub_node[6];
//#ifdef USE_MPI
//    if(worker_id==0){
//#endif
//        DebugOn("Final time "<<(get_wall_time()-ts)<<endl);
//        DebugOn("\n*************************\n");
//        DebugOn("L2 Final "<<best_ub<<endl);
//        DebugOn("*************************\n");
//        DebugOn("roll rad "<< rpy[0]<<endl);
//        DebugOn("pitch rad "<< rpy[1]<<endl);
//        DebugOn("yaw rad "<< rpy[2]<<endl);
//        DebugOn("tx "<< rpy[3]<<endl);
//        DebugOn("ty "<< rpy[4]<<endl);
//        DebugOn("tz "<< rpy[5]<<endl);
//#ifdef USE_MPI
//    }
//#endif
//
//    return rpy;
//
//}
#endif /* Heuristics_h */
