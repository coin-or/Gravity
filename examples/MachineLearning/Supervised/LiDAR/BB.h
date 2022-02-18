//
//  BB.h
//  Gravity
//
//  Created by Smitha on 8/30/21.
//

#ifndef BB_h
#define BB_h

#include "Lidar_preprocess.h"
#include "Lower_Bound.h"
#include <gravity/solver.h>

#include <thread>
using namespace std;
using namespace gravity;


void run_ub_parallel(const vector<vector<double>>& input_model_cloud, const vector<vector<double>>& input_data_cloud, const vector<vector<double>>& uav_model, const vector<vector<double>>& uav_data, const vector<vector<double>>& rpy_model, const vector<vector<double>>& rpy_data, const vector<vector<double>>& input_model_offset, const vector<vector<double>>& input_data_offset, const vector<double>& roll_lb,  const vector<double>& roll_ub,  const vector<double>& pitch_lb,  const vector<double>& pitch_ub,  const vector<double>& yaw_lb, const vector<double>& yaw_ub, string error_type, vector<double>& ub_node);
bool compute_upper_bound_mid(double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, vector<double>& best_rot, double& best_ub, vector<vector<double>>& input_model_cloud, vector<vector<double>>& input_data_cloud, vector<vector<double>>& uav_model, vector<vector<double>>& uav_data, const vector<vector<double>>& rpy_model, const vector<vector<double>>& rpy_data, vector<vector<double>>& input_model_offset, vector<vector<double>>& input_data_offset, string error_type);
void evaluate_upper_bound_mid(double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, vector<double>& res, const vector<vector<double>>& input_model_cloud, const vector<vector<double>>& input_data_cloud, const vector<vector<double>>& uav_model, const vector<vector<double>>& uav_data, const vector<vector<double>>& rpy_model, const vector<vector<double>>& rpy_data, const vector<vector<double>>& input_model_offset, const vector<vector<double>>& input_data_offset, string error_type);
void run_preprocess_parallel_Align(const vector<vector<double>>& input_cloud_model, const vector<vector<double>>& input_cloud_data, const vector<vector<double>>& uav_model, const vector<vector<double>>& uav_data, const vector<vector<double>>& rpy_model, const vector<vector<double>>& rpy_data,const vector<vector<double>>& input_model_offset, const vector<vector<double>>& input_data_offset, vector<int>& pos_vec, vector<shared_ptr<Model<double>>>& models, const vector<treenode_r>& vec_node, vector<int>& m_vec,vector<double>& vec_lb, vector<indices>& valid_cells, int nb_threads, double upper_bound, double lower_bound, vector<param<double>>& dist_cost_new, int iter, std::string error_type, vector<double>& ub_i);
void run_preprocess_model_Align(const vector<vector<double>>& input_cloud_model, const vector<vector<double>>& input_cloud_data, const vector<vector<double>>& uav_model, const vector<vector<double>>& uav_data, const vector<vector<double>>& rpy_model, const vector<vector<double>>& rpy_data, const vector<vector<double>>& input_model_offset, const vector<vector<double>>& input_data_offset, treenode_r vec_node_i, int& m_vec_i,  double& vec_lb_i,  indices& valid_cells_i, param<double>& dist_cost_i, double& prep_time_i, double upper_bound, shared_ptr<Model<double>>& model_i, std::string error_type, vector<double>& ub_i);
void run_preprocess_only_parallel(const vector<vector<double>>& input_model_cloud, const vector<vector<double>>& input_data_cloud, const vector<vector<double>>& uav_model, const vector<vector<double>>& uav_data, const vector<vector<double>>& rpy_model, const vector<vector<double>>& rpy_data, const vector<vector<double>>& input_model_offset, const vector<vector<double>>& input_data_offset, vector<treenode_r>& vec_node, vector<double>& vec_lb, vector<indices>& valid_cells, int nb_threads, double upper_bound, double lower_bound, int iter, string error_type);
void send_vector_new(const vector<size_t>& limits, vector<double>& vec_full, vector<double>& vec_worker);
vector<double> ub_heuristic(vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, vector<vector<double>>& uav_model, vector<vector<double>>& uav_data, vector<vector<double>>& rpy_model, vector<vector<double>>& rpy_data, vector<double>& best_rot_trans, double best_ub, std::string error_type, const double scanner_x, const double scanner_y,const double scanner_z, const double hr, const double hp,const double hy)
{
#ifdef USE_MPI
    int worker_id, nb_workers;
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
    auto err_size = MPI_Comm_size(MPI_COMM_WORLD, &nb_workers);
#endif
    vector<double> rpy(3);
    /* INPUT BOUNDS */
    
    /* INPUT BOUNDS */
    double time_start = get_wall_time();
    double total_time_max = 90000;
    double prep_time_total=0;
    
    
    
    double yaw_min = -2*pi/180., yaw_max = 2*pi/180., pitch_min =-2*pi/180.,pitch_max = 2*pi/180.,roll_min =-2*pi/180.,roll_max = 2*pi/180.;
    
    vector<vector<double>> input_data_cloud, input_model_cloud, input_data_offset, input_model_offset;
    
    generate_inputs(point_cloud_model, uav_model, rpy_model, scanner_x, scanner_y, scanner_z,hr,hp,hy, input_model_cloud, input_model_offset);
    generate_inputs(point_cloud_data, uav_data, rpy_data, scanner_x, scanner_y, scanner_z, hr,hp,hy,input_data_cloud, input_data_offset);
    
    indices N1 = range(1,point_cloud_data.size());
    indices N2 = range(1,point_cloud_model.size());
    int nd=point_cloud_data.size();
    vector<int> new_matching(nd);
    vector<double> res(nd);
    
    
    
    double max_time = 60;
    double max_time_init=60;
    bool max_time_increase=false;
    int max_iter = 1e6;
    int models_count=0, models_new_count=0;
    int infeasible_count=0;
    vector<pair<pair<int,int>,pair<int,int>>> incompatible_pairs;
    size_t nb_threads = std::thread::hardware_concurrency();
    //int nb_threads = 1;
    
    
    // nb_threads = 1;
    pair<double,double> roll_bounds_r, pitch_bounds_r, yaw_bounds_r;
    
    roll_bounds_r={roll_min, roll_max};
    pitch_bounds_r={pitch_min, pitch_max};
    yaw_bounds_r={yaw_min, yaw_max};
    
    vector<pair<double,double>> roll_bounds, pitch_bounds, yaw_bounds;
    
    vector<indices> valid_cells;
    vector<int> pos_vec;
    vector<double> vec_lb;
    vector<treenode_r> vec_node;
    vector<int> m_vec;
    vector<vector<double>> costs_upto_vec;
    
    vector<shared_ptr<Model<>>> models, models_new;
    double lb = 0, ub = 12, ub_=-1, best_lb = 0;
    int nb_pruned = 0;
    int depth_r=0, iter=0;
    vector<int> depth_vec;
    priority_queue<treenode_r> lb_queue;
    vector<double> costs_upto_init(nd,0.0);
    
    
    double min_cost_sum=0.0;
    
    param<double> dist_cost_r("dist_cost_r");
    indices valid_cells_r;
    
    
    
    lb_queue.push(treenode_r(roll_bounds_r, pitch_bounds_r, yaw_bounds_r, lb, ub, ub_, depth_r, valid_cells_r, false, dist_cost_r));
    treenode_r topnode=lb_queue.top();
    auto ts=get_wall_time();
    while ( lb_queue.top().depth<=6) {
        
        roll_bounds.clear();
        pitch_bounds.clear();
        yaw_bounds.clear();
        valid_cells.clear();
        depth_vec.clear();
        vec_node.clear();
        vec_lb.clear();
        
        costs_upto_vec.clear();
        iter++;
        models_count=0;
        models_new_count=0;
        topnode=lb_queue.top();
        prep_time_total=0;
        int step=8;
        
        for(auto i=0;!lb_queue.empty();i+=step){
            topnode=lb_queue.top();
            lb_queue.pop();
            
            double roll_increment,  pitch_increment, yaw_increment;
            
            
            
            roll_increment = (topnode.roll.second - topnode.roll.first)/2.0;
            pitch_increment = (topnode.pitch.second - topnode.pitch.first)/2.0;
            yaw_increment = (topnode.yaw.second - topnode.yaw.first)/2.0;
            roll_bounds.push_back({topnode.roll.first, topnode.roll.first+roll_increment});
            roll_bounds.push_back({topnode.roll.first, topnode.roll.first+roll_increment});
            roll_bounds.push_back({topnode.roll.first, topnode.roll.first+roll_increment});
            roll_bounds.push_back({topnode.roll.first, topnode.roll.first+roll_increment});
            roll_bounds.push_back({topnode.roll.first+roll_increment, topnode.roll.second});
            roll_bounds.push_back({topnode.roll.first+roll_increment, topnode.roll.second});
            roll_bounds.push_back({topnode.roll.first+roll_increment, topnode.roll.second});
            roll_bounds.push_back({topnode.roll.first+roll_increment, topnode.roll.second});
            pitch_bounds.push_back({topnode.pitch.first, topnode.pitch.first+pitch_increment});
            pitch_bounds.push_back({topnode.pitch.first, topnode.pitch.first+pitch_increment});
            pitch_bounds.push_back({topnode.pitch.first+pitch_increment, topnode.pitch.second});
            pitch_bounds.push_back({topnode.pitch.first+pitch_increment, topnode.pitch.second});
            pitch_bounds.push_back({topnode.pitch.first, topnode.pitch.first+pitch_increment});
            pitch_bounds.push_back({topnode.pitch.first, topnode.pitch.first+pitch_increment});
            pitch_bounds.push_back({topnode.pitch.first+pitch_increment, topnode.pitch.second});
            pitch_bounds.push_back({topnode.pitch.first+pitch_increment, topnode.pitch.second});
            yaw_bounds.push_back({topnode.yaw.first, topnode.yaw.first+yaw_increment});
            yaw_bounds.push_back({topnode.yaw.first+yaw_increment, topnode.yaw.second});
            yaw_bounds.push_back({topnode.yaw.first, topnode.yaw.first+yaw_increment});
            yaw_bounds.push_back({topnode.yaw.first+yaw_increment, topnode.yaw.second});
            yaw_bounds.push_back({topnode.yaw.first, topnode.yaw.first+yaw_increment});
            yaw_bounds.push_back({topnode.yaw.first+yaw_increment, topnode.yaw.second});
            yaw_bounds.push_back({topnode.yaw.first, topnode.yaw.first+yaw_increment});
            yaw_bounds.push_back({topnode.yaw.first+yaw_increment, topnode.yaw.second});
            for(auto k=0;k<8;k++){
                vec_node.push_back(treenode_r(roll_bounds[i+k],  pitch_bounds[i+k], yaw_bounds[i+k], topnode.lb, best_ub, -1.0, topnode.depth+1, topnode.valid_cells, false,topnode.dist_cost_cells));
                depth_vec.push_back(topnode.depth+1);
            }
            
            
            for(auto k=0;k<8;k++){
                auto a=compute_upper_bound_mid(roll_bounds[i+k].first, roll_bounds[i+k].second, pitch_bounds[i+k].first, pitch_bounds[i+k].second, yaw_bounds[i+k].first, yaw_bounds[i+k].second, best_rot_trans, best_ub, input_model_cloud, input_data_cloud, uav_model, uav_data, rpy_model, rpy_data, input_model_offset, input_data_offset, error_type);
                if(a){
                    DebugOn("new ub "<<best_ub<<" "<<(get_wall_time()-ts)<<endl);
                }
                
            }
            
            
        }
        //run_preprocess_parallel(pos_vec, models,vec_node, m_vec, valid_cells, vec_lb);
        
        
        
        for (int j = 0; j<vec_node.size(); j++) {
            
            
            if(vec_node[j].lb-1e-4<=best_ub)
            {
                lb_queue.push(vec_node[j]);
            }
        }
        
    }
    
    auto yaw_rad = atan2(-best_rot_trans[1],best_rot_trans[0]);
    auto roll_rad=asin(best_rot_trans[2]);
    auto pitch_rad=acos(best_rot_trans[8]/cos(roll_rad));
#ifdef USE_MPI
    if(worker_id==0){
#endif
        DebugOn("roll rad "<< roll_rad<<endl);
        DebugOn("pitch rad "<< pitch_rad<<endl);
        DebugOn("yaw rad "<< yaw_rad<<endl);
#ifdef USE_MPI
    }
#endif
    
    rpy[0]=roll_rad;
    rpy[1]=pitch_rad;
    rpy[2]=yaw_rad;
    return rpy;
    
}

vector<double> ub_heuristic_disc(vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, vector<vector<double>>& uav_model, vector<vector<double>>& uav_data, vector<vector<double>>& rpy_model, vector<vector<double>>& rpy_data, vector<double>& best_rot_trans, double& best_ub, std::string error_type, const double scanner_x, const double scanner_y,const double scanner_z, const double hr, const double hp,const double hy)
{
#ifdef USE_MPI
    int worker_id, nb_workers;
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
    auto err_size = MPI_Comm_size(MPI_COMM_WORLD, &nb_workers);
#endif
    vector<double> rpy(3, 0.0);
    /* INPUT BOUNDS */
    
    /* INPUT BOUNDS */
    double max_time=100;
    double prep_time_total=0;
    
    double yaw_min = -2*pi/180., yaw_max = 2*pi/180., pitch_min =-2*pi/180.,pitch_max = 2*pi/180.,roll_min =-2*pi/180.,roll_max = 2*pi/180.;
    
    vector<vector<double>> input_data_cloud, input_model_cloud, input_data_offset, input_model_offset;
    
    vector<double> roll_lb, roll_ub, pitch_lb, pitch_ub, yaw_lb, yaw_ub;
    
    generate_inputs(point_cloud_model, uav_model, rpy_model, scanner_x, scanner_y, scanner_z,hr,hp,hy, input_model_cloud, input_model_offset);
    generate_inputs(point_cloud_data, uav_data, rpy_data, scanner_x, scanner_y, scanner_z, hr,hp,hy,input_data_cloud, input_data_offset);
    
    indices N1 = range(1,point_cloud_data.size());
    indices N2 = range(1,point_cloud_model.size());
    int nd=point_cloud_data.size();
    vector<int> new_matching(nd);
    vector<double> res(nd);
    
    vector<double> ub_node(4,0);
    ub_node[0]=best_ub;
    size_t nb_threads = std::thread::hardware_concurrency();
    
    // nb_threads = 1;
    pair<double,double> roll_bounds_r, pitch_bounds_r, yaw_bounds_r;
    
    roll_bounds_r={roll_min, roll_max};
    pitch_bounds_r={pitch_min, pitch_max};
    yaw_bounds_r={yaw_min, yaw_max};
    
    double ts=get_wall_time();
    
    int ndisc=30;
    bool stop=false;
    int count_probs=0;
    while(!stop){
        double rb=(roll_bounds_r.second-roll_bounds_r.first)/ndisc;
        double pb=(pitch_bounds_r.second-pitch_bounds_r.first)/ndisc;
        double yb=(yaw_bounds_r.second-yaw_bounds_r.first)/ndisc;
        if(rb<=1e-9 || pb<=1e-9 ||yb<=1e-9){
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
                        run_ub_parallel(input_model_cloud,  input_data_cloud, uav_model, uav_data, rpy_model, rpy_data, input_model_offset, input_data_offset,  roll_lb,  roll_ub,  pitch_lb,   pitch_ub,  yaw_lb,  yaw_ub, error_type, ub_node);
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
        roll_bounds_r.first=roll_rad1-std::abs(roll_rad1)*0.1;
        roll_bounds_r.second=roll_rad1+std::abs(roll_rad1)*0.1;
        pitch_bounds_r.first=pitch_rad1-std::abs(pitch_rad1)*0.1;
        pitch_bounds_r.second=pitch_rad1+std::abs(pitch_rad1)*0.1;
        yaw_bounds_r.first=yaw_rad1-std::abs(yaw_rad1)*0.1;
        yaw_bounds_r.second=yaw_rad1+std::abs(yaw_rad1)*0.1;
        if((get_wall_time()-ts)>=max_time){
            stop=true;
            break;
        }
    }
    best_ub=ub_node[0];
    rpy[0]=ub_node[1];
    rpy[1]=ub_node[2];
    rpy[2]=ub_node[3];
#ifdef USE_MPI
    if(worker_id==0){
#endif
        DebugOn("final time "<<(get_wall_time()-ts)<<endl);
        DebugOn("final ub "<<best_ub<<endl);
        DebugOn("roll rad "<< rpy[0]<<endl);
        DebugOn("pitch rad "<< rpy[1]<<endl);
        DebugOn("yaw rad "<< rpy[2]<<endl);
#ifdef USE_MPI
    }
#endif

    return rpy;
    
}

vector<double> BranchBound_Align(vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, vector<vector<double>>& uav_model, vector<vector<double>>& uav_data, vector<vector<double>>& rpy_model, vector<vector<double>>& rpy_data, vector<double>& best_rot, double best_ub, std::string error_type, const double scanner_x, const double scanner_y,const double scanner_z, const double hr, const double hp,const double hy)
{
#ifdef USE_MPI
    int worker_id, nb_workers;
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
    auto err_size = MPI_Comm_size(MPI_COMM_WORLD, &nb_workers);
#endif
    vector<double> rpy, rpy_rad;
        auto roll_rad=best_rot[0];
        auto pitch_rad=best_rot[1];
        auto yaw_rad = best_rot[2];
        rpy_rad.push_back(roll_rad);
        rpy_rad.push_back(pitch_rad);
        rpy_rad.push_back(yaw_rad);
    /* INPUT BOUNDS */
    
    /* INPUT BOUNDS */
    double time_start = get_wall_time();
    double total_time_max = 90000;
    double prep_time_total=0;
    
    
    
    double yaw_min = -2*pi/180., yaw_max = 2*pi/180., pitch_min =-2*pi/180.,pitch_max = 2*pi/180.,roll_min =-2*pi/180.,roll_max = 2*pi/180.;
    
    vector<vector<double>> input_data_cloud, input_model_cloud, input_data_offset, input_model_offset;
    
    generate_inputs(point_cloud_model, uav_model, rpy_model, scanner_x, scanner_y, scanner_z,hr,hp,hy, input_model_cloud, input_model_offset);
    generate_inputs(point_cloud_data, uav_data, rpy_data, scanner_x, scanner_y, scanner_z, hr,hp,hy,input_data_cloud, input_data_offset);
    
    indices N1 = range(1,point_cloud_data.size());
    indices N2 = range(1,point_cloud_model.size());
    int nd=point_cloud_data.size();
    vector<int> new_matching(nd);
    vector<double> res(nd);
    double max_time = 30;
    double max_time_init=30;
    bool max_time_increase=false;
    int max_iter = 1e6;
    int models_count=0, models_new_count=0;
    int infeasible_count=0;
    vector<pair<pair<int,int>,pair<int,int>>> incompatible_pairs;
    size_t nb_threads = std::thread::hardware_concurrency();
    pair<double,double> roll_bounds_r, pitch_bounds_r, yaw_bounds_r;
    
    roll_bounds_r={roll_min, roll_max};
    pitch_bounds_r={pitch_min, pitch_max};
    yaw_bounds_r={yaw_min, yaw_max};
    
    vector<pair<double,double>> roll_bounds, pitch_bounds, yaw_bounds;
    
    vector<indices> valid_cells;
    vector<int> pos_vec;
    vector<double> vec_lb;
    vector<treenode_r> vec_node;
    vector<int> m_vec;
    vector<vector<double>> costs_upto_vec;
    auto point_cloud_data_copy=point_cloud_data;
    DebugOn("I will be using " << nb_threads << " parallel threads" << endl);
    vector<shared_ptr<Model<>>> models, models_new;
    double lb = 0, ub = 12, ub_=-1, best_lb = 0;
    int nb_pruned = 0;
    int depth_r=0, iter=0;
    vector<int> depth_vec, depth_vec_new;
    priority_queue<treenode_r> lb_queue;
    vector<double> costs_upto_init(nd,0.0);
    
    double min_cost_sum=0.0;
    
    vector<pair<double, double>> newmm_r;
    param<double> dist_cost_r("dist_cost_r");
    indices valid_cells_r;
    vector<param<double>> dist_cost_cells;
    double  prep_time=0.0;
    min_cost_sum=0;
    
    lb_queue.push(treenode_r(roll_bounds_r, pitch_bounds_r, yaw_bounds_r, lb, ub, ub_, depth_r, valid_cells_r, false, dist_cost_r));
    treenode_r topnode=lb_queue.top();
    
    while ( lb_queue.top().depth<=4) {
        
        roll_bounds.clear();
        pitch_bounds.clear();
        yaw_bounds.clear();
        valid_cells.clear();
        depth_vec.clear();
        vec_node.clear();
        vec_lb.clear();
        dist_cost_cells.clear();
        costs_upto_vec.clear();
        iter++;
        models_count=0;
        models_new_count=0;
        topnode=lb_queue.top();
        prep_time_total=0;
        int step=8;
        for(auto i=0;!lb_queue.empty();i+=step){
            topnode=lb_queue.top();
            lb_queue.pop();
            
            double roll_increment,  pitch_increment, yaw_increment;
            roll_increment = (topnode.roll.second - topnode.roll.first)/2.0;
            pitch_increment = (topnode.pitch.second - topnode.pitch.first)/2.0;
            yaw_increment = (topnode.yaw.second - topnode.yaw.first)/2.0;
            roll_bounds.push_back({topnode.roll.first, topnode.roll.first+roll_increment});
            roll_bounds.push_back({topnode.roll.first, topnode.roll.first+roll_increment});
            roll_bounds.push_back({topnode.roll.first, topnode.roll.first+roll_increment});
            roll_bounds.push_back({topnode.roll.first, topnode.roll.first+roll_increment});
            roll_bounds.push_back({topnode.roll.first+roll_increment, topnode.roll.second});
            roll_bounds.push_back({topnode.roll.first+roll_increment, topnode.roll.second});
            roll_bounds.push_back({topnode.roll.first+roll_increment, topnode.roll.second});
            roll_bounds.push_back({topnode.roll.first+roll_increment, topnode.roll.second});
            pitch_bounds.push_back({topnode.pitch.first, topnode.pitch.first+pitch_increment});
            pitch_bounds.push_back({topnode.pitch.first, topnode.pitch.first+pitch_increment});
            pitch_bounds.push_back({topnode.pitch.first+pitch_increment, topnode.pitch.second});
            pitch_bounds.push_back({topnode.pitch.first+pitch_increment, topnode.pitch.second});
            pitch_bounds.push_back({topnode.pitch.first, topnode.pitch.first+pitch_increment});
            pitch_bounds.push_back({topnode.pitch.first, topnode.pitch.first+pitch_increment});
            pitch_bounds.push_back({topnode.pitch.first+pitch_increment, topnode.pitch.second});
            pitch_bounds.push_back({topnode.pitch.first+pitch_increment, topnode.pitch.second});
            yaw_bounds.push_back({topnode.yaw.first, topnode.yaw.first+yaw_increment});
            yaw_bounds.push_back({topnode.yaw.first+yaw_increment, topnode.yaw.second});
            yaw_bounds.push_back({topnode.yaw.first, topnode.yaw.first+yaw_increment});
            yaw_bounds.push_back({topnode.yaw.first+yaw_increment, topnode.yaw.second});
            yaw_bounds.push_back({topnode.yaw.first, topnode.yaw.first+yaw_increment});
            yaw_bounds.push_back({topnode.yaw.first+yaw_increment, topnode.yaw.second});
            yaw_bounds.push_back({topnode.yaw.first, topnode.yaw.first+yaw_increment});
            yaw_bounds.push_back({topnode.yaw.first+yaw_increment, topnode.yaw.second});
            for(auto k=0;k<8;k++){
                vec_node.push_back(treenode_r(roll_bounds[i+k],  pitch_bounds[i+k], yaw_bounds[i+k], topnode.lb, best_ub, -1.0, topnode.depth+1, topnode.valid_cells, false,topnode.dist_cost_cells));
                depth_vec.push_back(topnode.depth+1);
            }
            
        }
        
        run_preprocess_only_parallel(input_model_cloud, input_data_cloud,uav_model, uav_data,  rpy_model, rpy_data, input_model_offset, input_data_offset,vec_node, vec_lb, valid_cells, nb_threads, best_ub, best_lb, iter, error_type);
        
        for (int j = 0; j<vec_node.size(); j++) {
            
            
            if(vec_node[j].lb-1e-4<=best_ub)
            {
                lb_queue.push(vec_node[j]);
            }
        }
    }
    best_lb = lb_queue.top().lb;
    double elapsed_time = get_wall_time() - time_start;
    double opt_gap = (best_ub - best_lb)/best_ub;
    double opt_gap_abs=(best_ub - best_lb);
    double max_opt_gap = 0.01;/* 5% opt gap */
    double eps=0.001;
    int prep_count=0;
    double ut_total=0;
    while (elapsed_time < total_time_max && lb_queue.top().lb<=best_ub && !lb_queue.empty() && opt_gap > max_opt_gap && !lb_queue.top().leaf && opt_gap_abs>0.1) {
        best_lb = lb_queue.top().lb;
        opt_gap = (best_ub - best_lb)/best_ub;
        DebugOn("Best UB so far = " << to_string_with_precision(best_ub,9) << endl);
        DebugOn("Best LB so far = " << to_string_with_precision(best_lb,9) << endl);
        DebugOn("Opt gap so far = " << to_string_with_precision(opt_gap*100,6) << "%\n");
        DebugOn("Queue size = " << lb_queue.size() << "\n");
        DebugOn("Elapsed time = " << elapsed_time << "seconds\n");
        DebugOn("iter "<<iter<<endl);
        if(elapsed_time >= total_time_max || opt_gap <= max_opt_gap)
            break;
        DebugOn("Total infeasible =  " << infeasible_count << endl);
        DebugOn("Total prep_time =  " << prep_time_total << endl);
        DebugOn("Total discarded =  " << prep_count << endl);
        double max_incr=0, max_ratio=1;
        pos_vec.clear();
        models.clear();
        
        roll_bounds.clear();
        pitch_bounds.clear();
        yaw_bounds.clear();
        valid_cells.clear();
        depth_vec.clear();
        m_vec.clear();
        vec_node.clear();
        vec_lb.clear();
        dist_cost_cells.clear();
        costs_upto_vec.clear();
        iter++;
        models_count=0;
        models_new_count=0;
        topnode=lb_queue.top();
        prep_time_total=0;
        ut_total=0;
        int step=8;
        for(auto i=0;i<nb_threads;i+=step){
            if(lb_queue.top().lb<=best_ub && !lb_queue.top().leaf && !lb_queue.empty()){
                topnode=lb_queue.top();
                lb_queue.pop();
                
                double roll_increment,  pitch_increment, yaw_increment;
                
                step=8;
                
                roll_increment = (topnode.roll.second - topnode.roll.first)/2.0;
                pitch_increment = (topnode.pitch.second - topnode.pitch.first)/2.0;
                yaw_increment = (topnode.yaw.second - topnode.yaw.first)/2.0;
                roll_bounds.push_back({topnode.roll.first, topnode.roll.first+roll_increment});
                roll_bounds.push_back({topnode.roll.first, topnode.roll.first+roll_increment});
                roll_bounds.push_back({topnode.roll.first, topnode.roll.first+roll_increment});
                roll_bounds.push_back({topnode.roll.first, topnode.roll.first+roll_increment});
                roll_bounds.push_back({topnode.roll.first+roll_increment, topnode.roll.second});
                roll_bounds.push_back({topnode.roll.first+roll_increment, topnode.roll.second});
                roll_bounds.push_back({topnode.roll.first+roll_increment, topnode.roll.second});
                roll_bounds.push_back({topnode.roll.first+roll_increment, topnode.roll.second});
                pitch_bounds.push_back({topnode.pitch.first, topnode.pitch.first+pitch_increment});
                pitch_bounds.push_back({topnode.pitch.first, topnode.pitch.first+pitch_increment});
                pitch_bounds.push_back({topnode.pitch.first+pitch_increment, topnode.pitch.second});
                pitch_bounds.push_back({topnode.pitch.first+pitch_increment, topnode.pitch.second});
                pitch_bounds.push_back({topnode.pitch.first, topnode.pitch.first+pitch_increment});
                pitch_bounds.push_back({topnode.pitch.first, topnode.pitch.first+pitch_increment});
                pitch_bounds.push_back({topnode.pitch.first+pitch_increment, topnode.pitch.second});
                pitch_bounds.push_back({topnode.pitch.first+pitch_increment, topnode.pitch.second});
                yaw_bounds.push_back({topnode.yaw.first, topnode.yaw.first+yaw_increment});
                yaw_bounds.push_back({topnode.yaw.first+yaw_increment, topnode.yaw.second});
                yaw_bounds.push_back({topnode.yaw.first, topnode.yaw.first+yaw_increment});
                yaw_bounds.push_back({topnode.yaw.first+yaw_increment, topnode.yaw.second});
                yaw_bounds.push_back({topnode.yaw.first, topnode.yaw.first+yaw_increment});
                yaw_bounds.push_back({topnode.yaw.first+yaw_increment, topnode.yaw.second});
                yaw_bounds.push_back({topnode.yaw.first, topnode.yaw.first+yaw_increment});
                yaw_bounds.push_back({topnode.yaw.first+yaw_increment, topnode.yaw.second});
                for(auto k=0;k<8;k++){
                    vec_node.push_back(treenode_r(roll_bounds[i+k],  pitch_bounds[i+k], yaw_bounds[i+k], topnode.lb, best_ub, -1.0, topnode.depth+1, topnode.valid_cells, false,topnode.dist_cost_cells));
                    depth_vec.push_back(topnode.depth+1);
                }
            }
            else{
                break;
            }
            if(lb_queue.empty()){
                break;
            }
            topnode = lb_queue.top();
        }
        elapsed_time = get_wall_time() - time_start;
        DebugOn("Elapsed time = " << elapsed_time << "seconds\n");
        if(elapsed_time + max_time > total_time_max){
            DebugOn("max time "<< max_time);
            break;
        }
        vector<double> ub_all(4);
        run_preprocess_parallel_Align(input_model_cloud, input_data_cloud,uav_model, uav_data,  rpy_model, rpy_data, input_model_offset, input_data_offset, pos_vec, models, vec_node, m_vec, vec_lb, valid_cells, nb_threads, best_ub, best_lb, dist_cost_cells, iter, error_type, ub_all);
        if(ub_all[0]<=best_ub){
            best_ub=ub_all[0];
            rpy_rad[0]=ub_all[1];
            rpy_rad[1]=ub_all[2];
            rpy_rad[2]=ub_all[3];
        }
        for (int j = 0; j<m_vec.size(); j++) {
            if(m_vec[j]==0){
                prep_count++;
            }
            else if(m_vec[j]==10){
                lb_queue.push(treenode_r(roll_bounds[j], pitch_bounds[j], yaw_bounds[j], vec_lb[j], best_ub, ub_, depth_vec[j], valid_cells[j], false, dist_cost_cells[j]));
            }
        }
        DebugOn("models size "<<models.size()<<endl);
        
        run_parallel(models, gurobi, 1e-4, nb_threads, "", max_iter, max_time, (best_ub));
        for (int j = 0; j<models.size(); j++) {
            int pos=pos_vec[j];
            if(models[j]->_status==0){
                auto lb_ =  models[j]->get_rel_obj_val();
                auto leaf_node=false;
                if(true){
                    lb = std::max(models[j]->get_rel_obj_val(), vec_lb[pos]);
                }
                if(lb-1e-4<=best_ub)
                {
                    lb_queue.push(treenode_r(roll_bounds[pos], pitch_bounds[pos], yaw_bounds[pos], lb, best_ub, ub_, depth_vec[pos], valid_cells[pos], leaf_node, dist_cost_cells[pos]));
                }
                else{
                    DebugOn("Infeasible lb "<<lb<<" "<<"best_ub "<<best_ub<<endl);
                    infeasible_count++;
                }
            }
            else{
                infeasible_count++;
            }
        }
        opt_gap = (best_ub - best_lb)/best_ub;
        opt_gap_abs=best_ub-best_lb;
        elapsed_time = get_wall_time() - time_start;
        
    }
    DebugOn("UB final "<<best_ub<<endl);
    DebugOn("LB final "<<best_lb<<endl);
    DebugOn("Gap final "<<(best_ub-best_lb)/best_ub*100.0<<endl);
    DebugOn("Elapsed time "<<elapsed_time<<endl);
    DebugOn("Total iter "<<iter<<endl);
    DebugOn("Queue size = " << lb_queue.size() << "\n");
    DebugOn("lb que top = " << lb_queue.top().lb << "\n");
    roll_rad= rpy_rad[0];
    pitch_rad=rpy_rad[1];
    yaw_rad = rpy_rad[2];
    DebugOn("roll rad "<< roll_rad<<endl);
    DebugOn("pitch rad "<< pitch_rad<<endl);
    DebugOn("yaw rad "<< yaw_rad<<endl);
    while(!lb_queue.empty())
    {
        auto node = lb_queue.top();
        DebugOn("node lb "<<node.lb<<" node.leaf "<<node.leaf<<endl);
        
        DebugOn(node.roll.first<<" "<< node.roll.second<<" "<<node.pitch.first<<" "<<node.pitch.second<<" "<<node.yaw.first<<" "<<node.yaw.second<<endl);
        if(node.roll.first-1e-6<=roll_rad && roll_rad<=node.roll.second+1e-6 && node.pitch.first-1e-6<=pitch_rad && pitch_rad<=node.pitch.second+1e-6 && node.yaw.first-1e-6<=yaw_rad && yaw_rad<=node.yaw.second+1e-6){
            DebugOn("True interval contained "<<endl);
        }
        lb_queue.pop();
    }
    
    DebugOn("roll rad "<< roll_rad<<endl);
    DebugOn("pitch rad "<< pitch_rad<<endl);
    DebugOn("yaw rad "<< yaw_rad<<endl);
    
    
    auto roll=roll_rad*180/pi;
    auto pitch=pitch_rad*180/pi;
    auto yaw=yaw_rad*180/pi;
    
    DebugOn("roll deg"<< roll<<endl);
    DebugOn("pitch deg"<< pitch<<endl);
    DebugOn("yaw deg"<< yaw<<endl);
    
    rpy.push_back(roll);
    rpy.push_back(pitch);
    rpy.push_back(yaw);
    return rpy;
}
void evaluate_upper_bound_mid(double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, vector<double>& res, const vector<vector<double>>& input_model_cloud, const vector<vector<double>>& input_data_cloud, const vector<vector<double>>& uav_model, const vector<vector<double>>& uav_data, const vector<vector<double>>& rpy_model, const vector<vector<double>>& rpy_data, const vector<vector<double>>& input_model_offset, const vector<vector<double>>& input_data_offset, string error_type){
    res.resize(4);
    double roll=(roll_min+roll_max)*0.5;
    double pitch=(pitch_min+pitch_max)*0.5;
    double yaw=(yaw_min+yaw_max)*0.5;
    
    vector<double> rot(9);
    
    
    rot[0]=cos(roll)*cos(yaw);
    rot[1]=(-1)*cos(roll)*sin(yaw);
    rot[2]=sin(roll);
    rot[3]=cos(pitch)*sin(yaw)+cos(yaw)*sin(roll)*sin(pitch);
    rot[4]=cos(pitch)*cos(yaw)-sin(roll)*sin(pitch)*sin(yaw);
    rot[5]=(-1)*cos(roll)*sin(pitch);
    rot[6]=sin(pitch)*sin(yaw)-cos(pitch)*cos(yaw)*sin(roll);
    rot[7]=cos(yaw)*sin(pitch)+cos(pitch)*sin(roll)*sin(yaw);
    rot[8]=cos(roll)*cos(pitch);
    
    vector<vector<double>> point_cloud_model_copy, point_cloud_data_copy;
    
    generate_outputs_from_inputs(roll, pitch, yaw, input_model_cloud, uav_model, rpy_model, input_model_offset, point_cloud_model_copy);
    generate_outputs_from_inputs(roll, pitch, yaw, input_data_cloud, uav_data, rpy_data, input_data_offset, point_cloud_data_copy);
    
    
    
    vector<int> matching1(input_data_cloud.size());
    vector<double> err_per_point1(input_data_cloud.size());
    double error;
    
    if(error_type=="L2"){
        error = computeL2error(point_cloud_model_copy,point_cloud_data_copy,matching1,err_per_point1);
    }
    else{
        error = computeL1error(point_cloud_model_copy,point_cloud_data_copy,matching1,err_per_point1);
    }
    res[0]=error;
    res[1]=roll;
    res[2]=pitch;
    res[3]=yaw;
    
}
bool compute_upper_bound_mid(double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, vector<double>& best_rot, double& best_ub, vector<vector<double>>& input_model_cloud, vector<vector<double>>& input_data_cloud, vector<vector<double>>& uav_model, vector<vector<double>>& uav_data, const vector<vector<double>>& rpy_model, const vector<vector<double>>& rpy_data, vector<vector<double>>& input_model_offset, vector<vector<double>>& input_data_offset, string error_type){
    bool res=false;
    double roll=(roll_min+roll_max)*0.5;
    double pitch=(pitch_min+pitch_max)*0.5;
    double yaw=(yaw_min+yaw_max)*0.5;
    
    vector<double> rot(9);
    
    
    rot[0]=cos(roll)*cos(yaw);
    rot[1]=(-1)*cos(roll)*sin(yaw);
    rot[2]=sin(roll);
    rot[3]=cos(pitch)*sin(yaw)+cos(yaw)*sin(roll)*sin(pitch);
    rot[4]=cos(pitch)*cos(yaw)-sin(roll)*sin(pitch)*sin(yaw);
    rot[5]=(-1)*cos(roll)*sin(pitch);
    rot[6]=sin(pitch)*sin(yaw)-cos(pitch)*cos(yaw)*sin(roll);
    rot[7]=cos(yaw)*sin(pitch)+cos(pitch)*sin(roll)*sin(yaw);
    rot[8]=cos(roll)*cos(pitch);
    
    vector<vector<double>> point_cloud_model_copy, point_cloud_data_copy;
    
    generate_outputs_from_inputs(roll, pitch, yaw, input_model_cloud, uav_model, rpy_model, input_model_offset, point_cloud_model_copy);
    generate_outputs_from_inputs(roll, pitch, yaw, input_data_cloud, uav_data, rpy_data, input_data_offset, point_cloud_data_copy);
    
    
    
    vector<int> matching1(input_data_cloud.size());
    vector<double> err_per_point1(input_data_cloud.size());
    double error;
    
    if(error_type=="L2"){
        error = computeL2error(point_cloud_model_copy,point_cloud_data_copy,matching1,err_per_point1);
    }
    else{
        error = computeL1error(point_cloud_model_copy,point_cloud_data_copy,matching1,err_per_point1);
    }
    
    
    if(error<=best_ub){
        best_ub=error;
        best_rot=rot;
        DebugOff("compute mid new ub "<<best_ub<<endl);
        res=true;
    }
    return res;
}
void run_ub_parallel(const vector<vector<double>>& input_model_cloud, const vector<vector<double>>& input_data_cloud, const vector<vector<double>>& uav_model, const vector<vector<double>>& uav_data, const vector<vector<double>>& rpy_model, const vector<vector<double>>& rpy_data, const vector<vector<double>>& input_model_offset, const vector<vector<double>>& input_data_offset, const vector<double>& roll_lb,  const vector<double>& roll_ub,  const vector<double>& pitch_lb,  const vector<double>& pitch_ub,  const vector<double>& yaw_lb, const vector<double>& yaw_ub, string error_type, vector<double>& ub_node){
    std::vector<thread> threads;
    
    int nd=input_data_cloud.size();
    int num=roll_lb.size();
    if(num==0){
        DebugOff("in run_parallel(models...), models is empty, returning");
    }
    vector<vector<double>> vec_ub;
    vec_ub.resize(num);
    
    for (auto i = 0; i < num; i++) {
        threads.push_back(thread(&evaluate_upper_bound_mid,roll_lb[i], roll_ub[i], pitch_lb[i], pitch_ub[i], yaw_lb[i], yaw_ub[i], ref(vec_ub[i]), ref(input_model_cloud), ref(input_data_cloud), ref(uav_model), ref(uav_data), ref(rpy_model), ref(rpy_data), ref(input_model_offset), ref(input_data_offset), error_type));
    }
    for(auto &t : threads){
        t.join();
    }
    threads.clear();

    for(auto i=0;i<num;i++){
        if(vec_ub[i][0]<=ub_node[0]){
            ub_node[0]=vec_ub[i][0];
            ub_node[1]=vec_ub[i][1];
            ub_node[2]=vec_ub[i][2];
            ub_node[3]=vec_ub[i][3];
        }
    }
}
void run_preprocess_parallel_Align(const vector<vector<double>>& input_model_cloud, const vector<vector<double>>& input_data_cloud, const vector<vector<double>>& uav_model, const vector<vector<double>>& uav_data, const vector<vector<double>>& rpy_model, const vector<vector<double>>& rpy_data, const vector<vector<double>>& input_model_offset, const vector<vector<double>>& input_data_offset, vector<int>& pos_vec, vector<shared_ptr<Model<double>>>& models, const vector<treenode_r>& vec_node, vector<int>& m_vec,vector<double>& vec_lb, vector<indices>& valid_cells, int nb_threads, double upper_bound, double lower_bound, vector<param<double>>& dist_cost_new, int iter, string error_type, vector<double>& ub_node){
    vector<shared_ptr<Model<double>>> temp_models;
    std::vector<thread> threads;
    ub_node.resize(4,1000);
    int nd=input_data_cloud.size();
    int num=vec_node.size();
    if(num==0){
        DebugOff("in run_parallel(models...), models is empty, returning");
    }
    vector<vector<double>> vec_ub;
    vector<param<double>> vec_dist_cost;
    for(auto i=0;i<num;i++){
        param<double> dist_cost("dist_cost");
        vec_dist_cost.push_back(dist_cost);
    }
    
    valid_cells.resize(num);
    m_vec.resize(num, 0);
    vec_lb.resize(num, 0.0);
    temp_models.resize(num);
    vec_ub.resize(num);
    
    vector<double> vec_prep_time;
    vec_prep_time.resize(num, 0.0);
    
    for (auto i = 0; i < num; i++) {
        threads.push_back(thread(&run_preprocess_model_Align, ref(input_model_cloud), ref(input_data_cloud), ref(uav_model), ref(uav_data), ref(rpy_model), ref(rpy_data), ref(input_model_offset), ref(input_data_offset), ref(vec_node[i]), ref(m_vec[i]),ref(vec_lb[i]), ref(valid_cells[i]), ref(vec_dist_cost[i]), ref(vec_prep_time[i]), ref(upper_bound), ref(temp_models[i]), ref(error_type), ref(vec_ub[i])));
    }
    for(auto &t : threads){
        t.join();
    }
    threads.clear();
    
    for(auto i=0;i<num;i++){
        vec_lb[i]=std::max(vec_lb[i], vec_node[i].lb);
        if(m_vec[i]==1){
            models.push_back(temp_models[i]);
            pos_vec.push_back(i);
        }
    }
    ub_node[0]=vec_ub[0][0];
    for(auto i=0;i<num;i++){
        if(vec_ub[i][0]<=ub_node[0]){
            ub_node[0]=vec_ub[i][0];
            ub_node[1]=vec_ub[i][1];
            ub_node[2]=vec_ub[i][2];
            ub_node[3]=vec_ub[i][3];
        }
    }
    dist_cost_new=vec_dist_cost;
}
void run_preprocess_model_Align(const vector<vector<double>>& input_model_cloud, const vector<vector<double>>& input_data_cloud, const vector<vector<double>>& uav_model, const vector<vector<double>>& uav_data, const vector<vector<double>>& rpy_model, const vector<vector<double>>& rpy_data, const vector<vector<double>>& input_model_offset, const vector<vector<double>>& input_data_offset, treenode_r vec_node_i, int& m_vec_i,  double& vec_lb_i,  indices& valid_cells_i, param<double>& dist_cost_i, double& prep_time_i, double upper_bound, shared_ptr<Model<double>>& model_i, string error_type, vector<double>& ub_i){
    
    evaluate_upper_bound_mid(vec_node_i.roll.first, vec_node_i.roll.second, vec_node_i.pitch.first, vec_node_i.pitch.second, vec_node_i.yaw.first ,vec_node_i.yaw.second, ub_i, input_model_cloud, input_data_cloud, uav_model, uav_data, rpy_model, rpy_data, input_model_offset,input_data_offset, error_type);
#ifdef USE_GJK
    preprocess_lid(input_model_cloud, input_data_cloud, uav_model, uav_data, rpy_model, rpy_data, input_model_offset,input_data_offset, vec_node_i.valid_cells, valid_cells_i,  dist_cost_i, vec_node_i.roll.first, vec_node_i.roll.second, vec_node_i.pitch.first, vec_node_i.pitch.second, vec_node_i.yaw.first ,vec_node_i.yaw.second, upper_bound, prep_time_i, vec_lb_i, error_type);
#endif
    bool model_created=false;
    if(valid_cells_i.size()>=input_data_cloud.size() && valid_cells_i.size()<=2e4){
        if(error_type=="L2"){
            model_i=Align_L2_model_rotation_trigonometric_scanner(input_model_cloud, input_data_cloud, uav_model, uav_data, rpy_model, rpy_data, input_model_offset, input_data_offset,vec_node_i.roll.first, vec_node_i.roll.second, vec_node_i.pitch.first, vec_node_i.pitch.second, vec_node_i.yaw.first ,vec_node_i.yaw.second, valid_cells_i, dist_cost_i);
        }
        else{
            model_i=Align_L2_model_rotation_trigonometric_scanner(input_model_cloud, input_data_cloud, uav_model, uav_data, rpy_model, rpy_data, input_model_offset, input_data_offset,vec_node_i.roll.first, vec_node_i.roll.second, vec_node_i.pitch.first, vec_node_i.pitch.second, vec_node_i.yaw.first ,vec_node_i.yaw.second, valid_cells_i, dist_cost_i);
        }
        model_created=true;
        m_vec_i=1;
    }
    else if(valid_cells_i.size()> 2e4){
        m_vec_i=10;
    }
    else if(valid_cells_i.size()<input_data_cloud.size()){
        m_vec_i=0;
    }
}
void run_preprocess_only_parallel(const vector<vector<double>>& input_model_cloud, const vector<vector<double>>& input_data_cloud, const vector<vector<double>>& uav_model, const vector<vector<double>>& uav_data, const vector<vector<double>>& rpy_model, const vector<vector<double>>& rpy_data, const vector<vector<double>>& input_model_offset, const vector<vector<double>>& input_data_offset, vector<treenode_r>& vec_node, vector<double>& vec_lb, vector<indices>& valid_cells, int nb_threads, double upper_bound, double lower_bound, int iter, string error_type){
    //size_t nb_threads = std::thread::hardware_concurrency();
    std::vector<thread> threads;
    int nd=input_data_cloud.size();
    int num=vec_node.size();
    if(num==0){
        DebugOff("in run_parallel(models...), models is empty, returning");
    }
    
    vector<param<double>> vec_dist_cost;
    for(auto i=0;i<num;i++){
        param<double> dist_cost("dist_cost");
        vec_dist_cost.push_back(dist_cost);
    }
    
    valid_cells.resize(num);
    vec_lb.resize(num, 0.0);
    
    vector<double> vec_prep_time;
    vec_prep_time.resize(num, 0.0);
    int npass=num/nb_threads+1;
    
    
    for (auto j = 0; j < npass; j++) {
        
        for (auto i = j*nb_threads; i < std::min((j+1)*nb_threads, num); i++) {
#ifdef USE_GJK
            threads.push_back(thread(&preprocess_lid, ref(input_model_cloud), ref(input_data_cloud), ref(uav_model), ref(uav_data), ref(rpy_model), ref(rpy_data), ref(input_model_offset), ref(input_data_offset), ref(vec_node[i].valid_cells), ref(valid_cells[i]),  ref(vec_dist_cost[i]), vec_node[i].roll.first, vec_node[i].roll.second, vec_node[i].pitch.first, vec_node[i].pitch.second, vec_node[i].yaw.first ,vec_node[i].yaw.second, upper_bound, ref(vec_prep_time[i]), ref(vec_lb[i]), error_type));
#endif
        }
        
        for(auto &t : threads){
            t.join();
        }
        threads.clear();
    }
    for (auto i = 0; i < num; i++) {
        if(valid_cells[i].size()>=input_data_cloud.size()){
            vec_node[i].lb=vec_lb[i];
            vec_node[i].valid_cells=valid_cells[i];
            vec_node[i].dist_cost_cells=vec_dist_cost[i];
        }
        else{
            vec_node[i].lb=1000*upper_bound;
        }
    }
}
#ifdef USE_MPI
vector<double> BranchBound_MPI(vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, vector<vector<double>>& uav_model, vector<vector<double>>& uav_data, vector<vector<double>>& rpy_model, vector<vector<double>>& rpy_data, vector<double>& best_rot, double best_ub, std::string error_type, const double scanner_x, const double scanner_y,const double scanner_z, const double hr, const double hp,const double hy)
{
    int worker_id, nb_workers;
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
    auto err_size = MPI_Comm_size(MPI_COMM_WORLD, &nb_workers);

    vector<double> rpy, rpy_rad;
    auto roll_rad=best_rot[0];
    auto pitch_rad=best_rot[1];
    auto yaw_rad = best_rot[2];
    rpy_rad.push_back(roll_rad);
    rpy_rad.push_back(pitch_rad);
    rpy_rad.push_back(yaw_rad);
    /* INPUT BOUNDS */
    
    /* INPUT BOUNDS */
    double time_start = get_wall_time();
    double total_time_max = 900000;
    double prep_time_total=0;
    
    double yaw_min = -2*pi/180., yaw_max = 2*pi/180., pitch_min =-2*pi/180.,pitch_max = 2*pi/180.,roll_min =-2*pi/180.,roll_max = 2*pi/180.;
    
    vector<vector<double>> input_data_cloud, input_model_cloud, input_data_offset, input_model_offset;
    
    generate_inputs(point_cloud_model, uav_model, rpy_model, scanner_x, scanner_y, scanner_z,hr,hp,hy, input_model_cloud, input_model_offset);
    generate_inputs(point_cloud_data, uav_data, rpy_data, scanner_x, scanner_y, scanner_z, hr,hp,hy,input_data_cloud, input_data_offset);
    
    indices N1 = range(1,point_cloud_data.size());
    indices N2 = range(1,point_cloud_model.size());
    int nd=point_cloud_data.size();
    vector<int> new_matching(nd);
    vector<double> res(nd);
    
    double max_time = 30;
    double max_time_init=30;
    bool max_time_increase=false;
    int max_iter = 1e6;
    int models_count=0, models_new_count=0;
    int infeasible_count=0;
    vector<pair<pair<int,int>,pair<int,int>>> incompatible_pairs;
    size_t nb_threads = std::thread::hardware_concurrency();
    int threads_total=nb_threads*nb_workers;
    pair<double,double> roll_bounds_r, pitch_bounds_r, yaw_bounds_r;
    
    roll_bounds_r={roll_min, roll_max};
    pitch_bounds_r={pitch_min, pitch_max};
    yaw_bounds_r={yaw_min, yaw_max};
    
    vector<pair<double,double>> roll_bounds, pitch_bounds, yaw_bounds;
    
    vector<indices> valid_cells;
    vector<int> pos_vec;
    vector<double> vec_lb;
    vector<treenode_r> vec_node;
    vector<int> m_vec;
    vector<vector<double>> costs_upto_vec;
    auto point_cloud_data_copy=point_cloud_data;
    if(worker_id==0){
        DebugOn("I will be using " << nb_threads << " parallel threads" << endl);
    }
    vector<shared_ptr<Model<>>> models, models_new;
    double lb = 0, ub = 12, ub_=-1, best_lb = 0;
    int nb_pruned = 0;
    int depth_r=0, iter=0;
    vector<int> depth_vec, depth_vec_new;
    priority_queue<treenode_r> lb_queue;
    vector<double> costs_upto_init(nd,0.0);
    double min_cost_sum=0.0;
    vector<pair<double, double>> newmm_r;
    param<double> dist_cost_r("dist_cost_r");
    indices valid_cells_r;
    vector<param<double>> dist_cost_cells;
    double  prep_time=0.0;
    min_cost_sum=0;
    lb_queue.push(treenode_r(roll_bounds_r, pitch_bounds_r, yaw_bounds_r, lb, ub, ub_, depth_r, valid_cells_r, false, dist_cost_r));
    treenode_r topnode=lb_queue.top();
    int step=8;
    while ( lb_queue.top().depth<=4) {
        roll_bounds.clear();
        pitch_bounds.clear();
        yaw_bounds.clear();
        valid_cells.clear();
        depth_vec.clear();
        vec_node.clear();
        vec_lb.clear();
        dist_cost_cells.clear();
        costs_upto_vec.clear();
        iter++;
        models_count=0;
        models_new_count=0;
        topnode=lb_queue.top();
        prep_time_total=0;
        for(auto i=0;!lb_queue.empty();i+=step){
            topnode=lb_queue.top();
            lb_queue.pop();
            double roll_increment,  pitch_increment, yaw_increment;

            roll_increment = (topnode.roll.second - topnode.roll.first)/2.0;
            pitch_increment = (topnode.pitch.second - topnode.pitch.first)/2.0;
            yaw_increment = (topnode.yaw.second - topnode.yaw.first)/2.0;
            roll_bounds.push_back({topnode.roll.first, topnode.roll.first+roll_increment});
            roll_bounds.push_back({topnode.roll.first, topnode.roll.first+roll_increment});
            roll_bounds.push_back({topnode.roll.first, topnode.roll.first+roll_increment});
            roll_bounds.push_back({topnode.roll.first, topnode.roll.first+roll_increment});
            roll_bounds.push_back({topnode.roll.first+roll_increment, topnode.roll.second});
            roll_bounds.push_back({topnode.roll.first+roll_increment, topnode.roll.second});
            roll_bounds.push_back({topnode.roll.first+roll_increment, topnode.roll.second});
            roll_bounds.push_back({topnode.roll.first+roll_increment, topnode.roll.second});
            pitch_bounds.push_back({topnode.pitch.first, topnode.pitch.first+pitch_increment});
            pitch_bounds.push_back({topnode.pitch.first, topnode.pitch.first+pitch_increment});
            pitch_bounds.push_back({topnode.pitch.first+pitch_increment, topnode.pitch.second});
            pitch_bounds.push_back({topnode.pitch.first+pitch_increment, topnode.pitch.second});
            pitch_bounds.push_back({topnode.pitch.first, topnode.pitch.first+pitch_increment});
            pitch_bounds.push_back({topnode.pitch.first, topnode.pitch.first+pitch_increment});
            pitch_bounds.push_back({topnode.pitch.first+pitch_increment, topnode.pitch.second});
            pitch_bounds.push_back({topnode.pitch.first+pitch_increment, topnode.pitch.second});
            yaw_bounds.push_back({topnode.yaw.first, topnode.yaw.first+yaw_increment});
            yaw_bounds.push_back({topnode.yaw.first+yaw_increment, topnode.yaw.second});
            yaw_bounds.push_back({topnode.yaw.first, topnode.yaw.first+yaw_increment});
            yaw_bounds.push_back({topnode.yaw.first+yaw_increment, topnode.yaw.second});
            yaw_bounds.push_back({topnode.yaw.first, topnode.yaw.first+yaw_increment});
            yaw_bounds.push_back({topnode.yaw.first+yaw_increment, topnode.yaw.second});
            yaw_bounds.push_back({topnode.yaw.first, topnode.yaw.first+yaw_increment});
            yaw_bounds.push_back({topnode.yaw.first+yaw_increment, topnode.yaw.second});
            for(auto k=0;k<8;k++){
                vec_node.push_back(treenode_r(roll_bounds[i+k],  pitch_bounds[i+k], yaw_bounds[i+k], topnode.lb, best_ub, -1.0, topnode.depth+1, topnode.valid_cells, false,topnode.dist_cost_cells));
                depth_vec.push_back(topnode.depth+1);
            }
        }
       
        run_preprocess_only_parallel(input_model_cloud, input_data_cloud,uav_model, uav_data,  rpy_model, rpy_data, input_model_offset, input_data_offset,vec_node, vec_lb, valid_cells, nb_threads, best_ub, best_lb, iter, error_type);
        
        for (int j = 0; j<vec_node.size(); j++) {
            if(vec_node[j].lb-1e-4<=best_ub)
            {
                lb_queue.push(vec_node[j]);
            }
        }
    }
    best_lb = lb_queue.top().lb;
    double elapsed_time = get_wall_time() - time_start;
    double opt_gap = (best_ub - best_lb)/best_ub;
    double opt_gap_abs=(best_ub - best_lb);
    double max_opt_gap = 0.01;/* 5% opt gap */
    double opt_gap_old=opt_gap+10;
    double eps=0.001;
    int prep_count=0;
    double ut_total=0;
    step = 8;
    while (elapsed_time < total_time_max && lb_queue.top().lb<=best_ub && !lb_queue.empty() && opt_gap > max_opt_gap && !lb_queue.top().leaf && opt_gap_abs>0.1) {
        if(worker_id==0){
            DebugOn("Best UB so far = " << to_string_with_precision(best_ub,9) << endl);
            DebugOn("Best LB so far = " << to_string_with_precision(best_lb,9) << endl);
            DebugOn("Opt gap so far = " << to_string_with_precision(opt_gap*100,6) << "%\n");
            DebugOn("Queue size = " << lb_queue.size() << "\n");
            DebugOn("Elapsed time = " << elapsed_time << "seconds\n");
            DebugOn("iter "<<iter<<endl);
        }
        if(elapsed_time >= total_time_max || opt_gap <= max_opt_gap)
            break;
        double max_incr=0, max_ratio=1;
        roll_bounds.clear();
        pitch_bounds.clear();
        yaw_bounds.clear();
        depth_vec.clear();
        vec_node.clear();
        iter++;
        topnode=lb_queue.top();
        prep_time_total=0;
        ut_total=0;
        for(auto i=0;i<threads_total;i+=step){
            if(lb_queue.top().lb<=best_ub && !lb_queue.top().leaf && !lb_queue.empty()){
                topnode=lb_queue.top();
                lb_queue.pop();
                double roll_increment,  pitch_increment, yaw_increment;
                step=8;
                roll_increment = (topnode.roll.second - topnode.roll.first)/2.0;
                pitch_increment = (topnode.pitch.second - topnode.pitch.first)/2.0;
                yaw_increment = (topnode.yaw.second - topnode.yaw.first)/2.0;
                roll_bounds.push_back({topnode.roll.first, topnode.roll.first+roll_increment});
                roll_bounds.push_back({topnode.roll.first, topnode.roll.first+roll_increment});
                roll_bounds.push_back({topnode.roll.first, topnode.roll.first+roll_increment});
                roll_bounds.push_back({topnode.roll.first, topnode.roll.first+roll_increment});
                roll_bounds.push_back({topnode.roll.first+roll_increment, topnode.roll.second});
                roll_bounds.push_back({topnode.roll.first+roll_increment, topnode.roll.second});
                roll_bounds.push_back({topnode.roll.first+roll_increment, topnode.roll.second});
                roll_bounds.push_back({topnode.roll.first+roll_increment, topnode.roll.second});
                pitch_bounds.push_back({topnode.pitch.first, topnode.pitch.first+pitch_increment});
                pitch_bounds.push_back({topnode.pitch.first, topnode.pitch.first+pitch_increment});
                pitch_bounds.push_back({topnode.pitch.first+pitch_increment, topnode.pitch.second});
                pitch_bounds.push_back({topnode.pitch.first+pitch_increment, topnode.pitch.second});
                pitch_bounds.push_back({topnode.pitch.first, topnode.pitch.first+pitch_increment});
                pitch_bounds.push_back({topnode.pitch.first, topnode.pitch.first+pitch_increment});
                pitch_bounds.push_back({topnode.pitch.first+pitch_increment, topnode.pitch.second});
                pitch_bounds.push_back({topnode.pitch.first+pitch_increment, topnode.pitch.second});
                yaw_bounds.push_back({topnode.yaw.first, topnode.yaw.first+yaw_increment});
                yaw_bounds.push_back({topnode.yaw.first+yaw_increment, topnode.yaw.second});
                yaw_bounds.push_back({topnode.yaw.first, topnode.yaw.first+yaw_increment});
                yaw_bounds.push_back({topnode.yaw.first+yaw_increment, topnode.yaw.second});
                yaw_bounds.push_back({topnode.yaw.first, topnode.yaw.first+yaw_increment});
                yaw_bounds.push_back({topnode.yaw.first+yaw_increment, topnode.yaw.second});
                yaw_bounds.push_back({topnode.yaw.first, topnode.yaw.first+yaw_increment});
                yaw_bounds.push_back({topnode.yaw.first+yaw_increment, topnode.yaw.second});
                for(auto k=0;k<8;k++){
                    vec_node.push_back(treenode_r(roll_bounds[i+k],  pitch_bounds[i+k], yaw_bounds[i+k], topnode.lb, best_ub, -1.0, topnode.depth+1, topnode.valid_cells, false,topnode.dist_cost_cells));
                    depth_vec.push_back(topnode.depth+1);
                }
            }
            else{
                break;
            }
            if(lb_queue.empty()){
                break;
            }
            topnode = lb_queue.top();
        }
        elapsed_time = get_wall_time() - time_start;
        if(worker_id==0){
            DebugOn("Elapsed time = " << elapsed_time << "seconds\n");
        }
        if(elapsed_time + max_time > total_time_max){
            DebugOn("max time "<< max_time);
            break;
        }
        auto nb_workers_ = std::min((size_t)nb_workers, vec_node.size());
        auto limits=bounds(nb_workers_, vec_node.size());
        vector<treenode_r> vec_node_worker;
        vector<double> lb_vector(vec_node.size(), -1.0);
        vector<double> lb_vec_worker;
        vector<double> ub_node(4,1000);
        vector<double> ub_all_node(4*nb_workers, 1000);
        if(worker_id+1<limits.size()){
            lb_vec_worker.resize(limits[worker_id+1]-limits[worker_id], -1);
            for (auto i = limits[worker_id]; i < limits[worker_id+1]; i++) {
                vec_node_worker.push_back(vec_node[i]);
            }
            valid_cells.clear();
            vec_lb.clear();
            dist_cost_cells.clear();
            m_vec.clear();
            pos_vec.clear();
            models.clear();
            run_preprocess_parallel_Align(input_model_cloud, input_data_cloud,uav_model, uav_data,  rpy_model, rpy_data, input_model_offset, input_data_offset, pos_vec, models, vec_node_worker, m_vec, vec_lb, valid_cells, nb_threads, best_ub, best_lb, dist_cost_cells, iter, error_type, ub_node);
            best_ub=std::min(best_ub, ub_node[0]);
            run_parallel(models, gurobi, 1e-4, nb_threads, "", max_iter, max_time, (best_ub));
            int count=0;
            for (int j = 0; j<m_vec.size(); j++) {
                if(m_vec[j]==0){
                    lb_vec_worker[j]=best_ub+10;
                    prep_count++;
                }
                else if(m_vec[j]==10){
                    lb_vec_worker[j]=vec_lb[j];
                    vec_node_worker[j].valid_cells=valid_cells[j];
                    vec_node_worker[j].lb=lb_vec_worker[j];
                }
                else if(m_vec[j]==1){
                    if(models[count]->_status==0){
                        lb_vec_worker[j]=std::max(vec_lb[j], models[count]->get_rel_obj_val());
                        vec_node_worker[j].valid_cells=valid_cells[j];
                        vec_node_worker[j].lb=lb_vec_worker[j];
                    }
                    else{
                        lb_vec_worker[j]=best_ub+10;;
                    }
                    count++;
                }
            }
        }
        send_vector_new(limits, lb_vector, lb_vec_worker);
        for(auto i=0;i<lb_vector.size();i++){
            if(lb_vector[i]<=best_ub){
                if(i>=limits[worker_id] && i<limits[worker_id+1]){
                    lb_queue.push(vec_node_worker[i-limits[worker_id]]);
                }
                else{
                    vec_node[i].lb=lb_vector[i];
                    lb_queue.push(vec_node[i]);
                }
            }
        }
        std::vector<size_t> limits_ub(nb_workers+1);
        for(auto l=0;l<=nb_workers;l++){
            limits_ub[l]=l*4;
        }
        for (auto l = limits_ub[worker_id]; l < limits_ub[worker_id+1]; l++) {
            ub_all_node[l]=ub_node[l-limits_ub[worker_id]];
        }
        std::vector<int> d, counts;
        for(auto l=limits_ub.begin()+1;l!=limits_ub.end();l++){
            counts.push_back(*l-*(l-1));
            d.push_back(*(l-1));
        }
        MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                       &ub_all_node[0], &counts[0], &d[0], MPI_DOUBLE, MPI_COMM_WORLD);
        for(auto l=0;l<ub_all_node.size();l+=4){
            if(ub_all_node[l]<=best_ub){
                best_ub=ub_all_node[l];
                rpy_rad[0]=ub_all_node[l+1];
                rpy_rad[1]=ub_all_node[l+2];
                rpy_rad[2]=ub_all_node[l+3];
            }
        }
        best_lb = lb_queue.top().lb;
        opt_gap_abs=best_ub-best_lb;
        opt_gap = (best_ub - best_lb)/best_ub;
        elapsed_time = get_wall_time() - time_start;
    }
    if(worker_id==0){
        DebugOn("UB final "<<best_ub<<endl);
        DebugOn("LB final "<<best_lb<<endl);
        DebugOn("Gap final "<<(best_ub-best_lb)/best_ub*100.0<<endl);
        DebugOn("Elapsed time "<<elapsed_time<<endl);
        DebugOn("Total iter "<<iter<<endl);
        DebugOn("Queue size = " << lb_queue.size() << "\n");
        DebugOn("lb que top = " << lb_queue.top().lb << "\n");
    }
    
    roll_rad=rpy_rad[0];
    pitch_rad=rpy_rad[1];
    yaw_rad = rpy_rad[2];
    if(worker_id==0){
        DebugOn("roll rad "<< roll_rad<<endl);
        DebugOn("pitch rad "<< pitch_rad<<endl);
        DebugOn("yaw rad "<< yaw_rad<<endl);
    }
    if(worker_id==0){
        while(true && !lb_queue.empty())
        {
            auto node = lb_queue.top();
            DebugOn("node lb "<<node.lb<<" node.leaf "<<node.leaf<<endl);
            
            DebugOn(node.roll.first<<" "<< node.roll.second<<" "<<node.pitch.first<<" "<<node.pitch.second<<" "<<node.yaw.first<<" "<<node.yaw.second<<endl);

            if(node.roll.first-1e-6<=roll_rad && roll_rad<=node.roll.second+1e-6 && node.pitch.first-1e-6<=pitch_rad && pitch_rad<=node.pitch.second+1e-6 && node.yaw.first-1e-6<=yaw_rad && yaw_rad<=node.yaw.second+1e-6){
                DebugOn("True interval contained "<<endl);
            }
            lb_queue.pop();
        }
    }
    if(worker_id==0){
        DebugOn("roll rad "<< roll_rad<<endl);
        DebugOn("pitch rad "<< pitch_rad<<endl);
        DebugOn("yaw rad "<< yaw_rad<<endl);
    }
    
    
    auto roll=roll_rad*180/pi;
    auto pitch=pitch_rad*180/pi;
    auto yaw=yaw_rad*180/pi;
    if(worker_id==0){
        DebugOn("roll deg "<< roll<<endl);
        DebugOn("pitch deg "<< pitch<<endl);
        DebugOn("yaw deg "<< yaw<<endl);
    }
    
    rpy.push_back(roll);
    rpy.push_back(pitch);
    rpy.push_back(yaw);
    return rpy;
}
#endif
#ifdef USE_MPI
void send_vector_new(const vector<size_t>& limits, vector<double>& vec_full, vector<double>& vec_worker){
    int worker_id, nb_workers;
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
    auto err_size = MPI_Comm_size(MPI_COMM_WORLD, &nb_workers);
    auto nb_workers_ =  limits.size()-1;
    DebugOff("I'm worker ID: " << worker_id << ", I'm getting ready to send my status " << endl);
    if(worker_id+1<limits.size()){
        for (auto i = limits[worker_id]; i < limits[worker_id+1]; i++) {
            vec_full[i]=vec_worker[i-limits[worker_id]];
        }
    }
    std::vector<int> d, counts;
    for(auto l=limits.begin()+1;l!=limits.end();l++){
        counts.push_back(*l-*(l-1));
        d.push_back(*(l-1));
    }
    for(auto l=nb_workers_;l!=nb_workers;l++){
        counts.push_back(0);
        d.push_back(limits.back());
    }
    if(counts.size()!=nb_workers){
        DebugOn("Error in size of counts");
    }
    if(d.size()!=nb_workers){
        DebugOn("Error in size of d");
    }
    MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                   &vec_full[0], &counts[0], &d[0], MPI_DOUBLE, MPI_COMM_WORLD);
    //MPI_Barrier(MPI_COMM_WORLD);
}
#endif
#endif /* BB_h */
