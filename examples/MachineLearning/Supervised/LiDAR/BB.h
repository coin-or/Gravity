//
//  BB.h
//  Gravity
//
//  Created by Smitha on 8/30/21.
//

#ifndef BB_h
#define BB_h

#include "Lidar_preprocess.h"

#include <thread>

void compute_upper_bound_mid(double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, vector<double>& best_rot_trans, double& best_ub, vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, vector<vector<double>>& uav_model, vector<vector<double>>& uav_data, string error_type);
void run_preprocess_parallel_Align(const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data, const vector<vector<double>>& uav_model, const vector<vector<double>>& uav_data, vector<int>& pos_vec, vector<shared_ptr<Model<double>>>& models, const vector<treenode_r>& vec_node, vector<int>& m_vec,vector<double>& vec_lb, vector<indices>& valid_cells, int nb_threads, double upper_bound, double lower_bound, vector<param<double>>& dist_cost_new, int iter, std::string error_type);
void run_preprocess_model_Align(const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data, const vector<vector<double>>& uav_model, const vector<vector<double>>& uav_data, treenode_r vec_node_i, int& m_vec_i,  double& vec_lb_i,  indices& valid_cells_i, param<double>& dist_cost_i, double& prep_time_i, double upper_bound, shared_ptr<Model<double>>& model_i, std::string error_type);

vector<double> BranchBound_Align(vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, vector<vector<double>>& uav_model, vector<vector<double>>& uav_data, vector<double>& best_rot_trans, double best_ub, std::string error_type)
{
    /* INPUT BOUNDS */
    
    /* INPUT BOUNDS */
    double time_start = get_wall_time();
    double total_time_max = 90000;
    double prep_time_total=0;
    
    double yaw_min = -2*pi/180., yaw_max = 2*pi/180., pitch_min =-2*pi/180.,pitch_max = 2*pi/180.,roll_min =-2*pi/180.,roll_max = 2*pi/180.;
    
    
    indices N1 = range(1,point_cloud_data.size());
    indices N2 = range(1,point_cloud_model.size());
    int nd=point_cloud_data.size();
    
    
    auto pcm=point_cloud_model;
    auto pcd=point_cloud_data;
    auto uavm=uav_model;
    auto uavd=uav_data;
    
    
    auto resi=run_IPH(pcm, pcd, uavm, uavd);
    //
    double roll_deg_i=get<0>(resi);
    double pitch_deg_i=get<1>(resi);
    double yaw_deg_i=get<2>(resi);
    double roll_rad_i=roll_deg_i*pi/180;
    double pitch_rad_i=pitch_deg_i*pi/180;
    double yaw_rad_i=yaw_deg_i*pi/180;
    double errori;
    vector<int> new_matching(nd);
    vector<double> res(nd);
    if(roll_rad_i>=roll_min && roll_rad_i<=roll_max && pitch_rad_i>=pitch_min && pitch_rad_i<=pitch_max && yaw_rad_i>=yaw_min && yaw_rad_i<=yaw_max){
        if(error_type=="L2"){
            errori= computeL2error(pcm,pcd,new_matching,res);
        }
        else{
            errori= computeL1error(pcm,pcd,new_matching,res);
        }
        if(errori<=best_ub){
            best_ub=errori;
            best_rot_trans[0]=cos(yaw_deg_i)*cos(roll_deg_i);
            best_rot_trans[1] = cos(yaw_deg_i)*sin(roll_deg_i)*sin(pitch_deg_i) - sin(yaw_deg_i)*cos(pitch_deg_i);
            best_rot_trans[2] = cos(yaw_deg_i)*sin(roll_deg_i)*cos(pitch_deg_i) + sin(yaw_deg_i)*sin(pitch_deg_i);
            best_rot_trans[3] = sin(yaw_deg_i)*cos(roll_deg_i);
            best_rot_trans[4] = sin(yaw_deg_i)*sin(roll_deg_i)*sin(pitch_deg_i) + cos(yaw_deg_i)*cos(pitch_deg_i);
            best_rot_trans[5] = sin(yaw_deg_i)*sin(roll_deg_i)*cos(pitch_deg_i) - cos(yaw_deg_i)*sin(pitch_deg_i);
            best_rot_trans[6] = sin(-1*roll_deg_i);
            best_rot_trans[7] = cos(roll_deg_i)*sin(pitch_deg_i);
            best_rot_trans[8] = cos(roll_deg_i)*cos(pitch_deg_i);
        }
    }
    
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
    auto point_cloud_data_copy=point_cloud_data;
    DebugOn("I will be using " << nb_threads << " parallel threads" << endl);
    vector<shared_ptr<Model<>>> models, models_new;
    double lb = 0, ub = 12, ub_=-1, best_lb = 0;
    int nb_pruned = 0;
    int depth_r=0, iter=0;
    vector<int> depth_vec, depth_vec_new;
    priority_queue<treenode_r> lb_queue;
    vector<double> costs_upto_init(nd,0.0);
    auto valid_cells_ro=indices(N1,N2);
    
    double min_cost_sum=0.0;
    
    vector<pair<double, double>> newmm_r;
    param<double> dist_cost_r("dist_cost_r");
    indices valid_cells_r;
    vector<param<double>> dist_cost_cells;
    double  prep_time=0.0;
    min_cost_sum=0;
    
    
    param<double> dist_cost_old("dist_cost_old");
    for(auto key:*valid_cells_ro._keys){
        dist_cost_old.add_val(key, 0.0);
    }
    
    
    lb_queue.push(treenode_r(roll_bounds_r, pitch_bounds_r, yaw_bounds_r, lb, ub, ub_, depth_r, valid_cells_r, false, dist_cost_r));
    double max_incr=0, max_ratio=1;
    int test_ub=16;
    treenode_r topnode=lb_queue.top();
    
    int ndisc=0;
    
    int max_cell_size=150000;
    
    double elapsed_time = get_wall_time() - time_start;
    double opt_gap = (best_ub - best_lb)/best_ub;
    double max_opt_gap = 0.01;/* 5% opt gap */
    double opt_gap_old=opt_gap+10;
    double eps=0.001;
    int prep_count=0;
    double ut_total=0;
    while (elapsed_time < total_time_max && lb_queue.top().lb<=best_ub && !lb_queue.empty() && opt_gap > max_opt_gap && !lb_queue.top().leaf) {
        best_lb = lb_queue.top().lb;
        opt_gap = (best_ub - best_lb)/best_ub;
        if(opt_gap_old-opt_gap <= eps){
            //  max_time=std::min(max_time*2, 120.0);
            //  max_time_increase=true;
        }
        if(opt_gap_old-opt_gap > eps && max_time_increase){
            max_time=std::max(max_time/2.0, max_time_init);
            if(max_time==max_time_init){
                // max_time_increase=false;
            }
        }
        opt_gap_old=opt_gap;
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
        max_incr=0, max_ratio=1;
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
                
                
                auto ut1=get_wall_time();
                if (i%test_ub==0){
                    for(auto k=0;k<8;k++){
                        compute_upper_bound_mid(roll_bounds[i+k].first, roll_bounds[i+k].second, pitch_bounds[i+k].first, pitch_bounds[i+k].second, yaw_bounds[i+k].first, yaw_bounds[i+k].second, best_rot_trans, best_ub, point_cloud_model, point_cloud_data, uav_model, uav_data, error_type);
                    }
                }
                auto ut2=get_wall_time();
                ut_total+=ut2-ut1;
            }
            else{
                break;
            }
            if(lb_queue.empty()){
                break;
            }
            topnode = lb_queue.top();
            DebugOff("lb "<<topnode.lb<<" ub "<<topnode.ub<<" size "<< models.size()<<endl);
            DebugOff("models.s "<<models.size()<<endl);
            DebugOff("validc.s "<<valid_cells.size()<<endl);
            DebugOff("shiftx.s "<<shift_x_bounds.size()<<endl);
            DebugOff("i "<<i<<endl);
            DebugOff("yaw.s "<<yaw_bounds.size()<<endl);
        }
        //run_preprocess_parallel(pos_vec, models,vec_node, m_vec, valid_cells, vec_lb);
        elapsed_time = get_wall_time() - time_start;
        DebugOn("Elapsed time = " << elapsed_time << "seconds\n");
        if(elapsed_time + max_time > total_time_max){
            DebugOn("max time "<< max_time);
            break;
        }
        DebugOn("upper bound time "<<ut_total<<endl);
        
        run_preprocess_parallel_Align(point_cloud_model, point_cloud_data,uav_model, uav_data,  pos_vec, models, vec_node, m_vec, vec_lb, valid_cells, nb_threads, best_ub, best_lb, dist_cost_cells, iter, error_type);
        
        for (int j = 0; j<m_vec.size(); j++) {
            if(!m_vec[j]){
                prep_count++;
            }
        }
        DebugOn("models size "<<models.size()<<endl);
        run_parallel(models, gurobi, 1e-4, nb_threads, "", max_iter, max_time, (best_ub));
        for (int j = 0; j<models.size(); j++) {
            int pos=pos_vec[j];
            if(models[j]->_status==0){
                ub_ = models[j]->get_obj_val();
                auto lb_ =  models[j]->get_rel_obj_val();
                auto leaf_node=false;
                if(ub_>=0){
                    if(ub_<=best_ub-1e-4){
                        vector<double> rot(9);
                        bool is_rotation  = get_solution(models[j], rot, new_matching);
                        DebugOn("gurobi found ub rotation "<<is_rotation<<endl);
                        if(is_rotation){
                            point_cloud_data_copy=point_cloud_data;
                            auto point_cloud_model_copy=point_cloud_model;
                            auto pitch_rad1 = atan2(rot[7], rot[8]);
                            auto roll_rad1 = atan2(-rot[6], std::sqrt(rot[7]*rot[7]+rot[8]*rot[8]));
                            auto yaw_rad1 = atan2(rot[3],rot[0]);
                            apply_rotation(roll_rad1*180/pi, pitch_rad1*180/pi, yaw_rad1*180/pi, point_cloud_model_copy, point_cloud_data_copy, uav_model, uav_data);
                            double err=0;
                            if(error_type=="L2"){
                                err=computeL2error(point_cloud_model_copy, point_cloud_data_copy, new_matching, res);
                            }
                            else{
                                err=computeL1error(point_cloud_model_copy, point_cloud_data_copy, new_matching, res);
                            }
                            DebugOn("gurobi found ub "<<err<<endl);
                            if(err<=best_ub){
                                best_ub=err;
                                best_rot_trans=rot;
                                DebugOn("new best ub "<<best_ub<<" ub_ "<<ub_<<" lb_ "<<lb_<<endl);
                                if((ub_-lb_/ub_)<=1e-6){
                                    leaf_node=true;
                                    DebugOn("leaf lb "<<lb_<<"  "<<error_type<<" "<<err<<" ub_ "<<ub_<<endl);
                                }
                            }
                        }
                    }
                }
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
        
        elapsed_time = get_wall_time() - time_start;
        
    }
    DebugOn("UB final "<<best_ub<<endl);
    DebugOn("LB final "<<best_lb<<endl);
    DebugOn("Gap final "<<(best_ub-best_lb)/best_ub*100.0<<endl);
    DebugOn("Elapsed time "<<elapsed_time<<endl);
    DebugOn("Total iter "<<iter<<endl);
    DebugOn("Queue size = " << lb_queue.size() << "\n");
    DebugOn("lb que top = " << lb_queue.top().lb << "\n");
    auto pitch_rad = atan2(best_rot_trans[7], best_rot_trans[8]);
    auto roll_rad = atan2(-best_rot_trans[6], std::sqrt(best_rot_trans[7]*best_rot_trans[7]+best_rot_trans[8]*best_rot_trans[8]));
    auto yaw_rad = atan2(best_rot_trans[3],best_rot_trans[0]);
    DebugOn("roll rad "<< roll_rad<<endl);
    DebugOn("pitch rad "<< pitch_rad<<endl);
    DebugOn("yaw rad "<< yaw_rad<<endl);
    while(!lb_queue.empty())
    {
        auto node = lb_queue.top();
        DebugOn("node lb "<<node.lb<<" node.leaf "<<node.leaf<<endl);
        
        DebugOn(node.roll.first<<" "<< node.roll.second<<" "<<node.pitch.first<<" "<<node.pitch.second<<" "<<node.yaw.first<<" "<<node.yaw.second<<endl);
        
        auto pitchr = atan2(best_rot_trans[7], best_rot_trans[8]);
        auto rollr = atan2(-best_rot_trans[6], std::sqrt(best_rot_trans[7]*best_rot_trans[7]+best_rot_trans[8]*best_rot_trans[8]));
        auto yawr = atan2(best_rot_trans[3],best_rot_trans[0]);
        if(node.roll.first-1e-3<=rollr && rollr<=node.roll.second+1e-3 && node.pitch.first-1e-3<=pitchr && pitchr<=node.pitch.second+1e-3 && node.yaw.first-1e-3<=yawr && yawr<=node.yaw.second+1e-3){
            DebugOn("True interval contained "<<endl);
        }
        lb_queue.pop();
        //        if(node.lb-1e-6 > best_ub)
        //            break;
    }
    pitch_rad = atan2(best_rot_trans[7], best_rot_trans[8]);
    roll_rad = atan2(-best_rot_trans[6], std::sqrt(best_rot_trans[7]*best_rot_trans[7]+best_rot_trans[8]*best_rot_trans[8]));
    yaw_rad = atan2(best_rot_trans[3],best_rot_trans[0]);
    DebugOn("roll rad "<< roll_rad<<endl);
    DebugOn("pitch rad "<< pitch_rad<<endl);
    DebugOn("yaw rad "<< yaw_rad<<endl);
    
    
    auto pitch = atan2(best_rot_trans[7], best_rot_trans[8])*180/pi;
    auto roll = atan2(-best_rot_trans[6], std::sqrt(best_rot_trans[7]*best_rot_trans[7]+best_rot_trans[8]*best_rot_trans[8]))*180/pi;
    auto yaw = atan2(best_rot_trans[3],best_rot_trans[0])*180/pi;
    
    DebugOn("roll deg"<< roll<<endl);
    DebugOn("pitch deg"<< pitch<<endl);
    DebugOn("yaw deg"<< yaw<<endl);
    
    return best_rot_trans;
    
}
void compute_upper_bound_mid(double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, vector<double>& best_rot, double& best_ub, vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, vector<vector<double>>& uav_model, vector<vector<double>>& uav_data, string error_type){
    double roll=(roll_min+roll_max)*0.5;
    double pitch=(pitch_min+pitch_max)*0.5;
    double yaw=(yaw_min+yaw_max)*0.5;
    
    vector<double> rot(9);
    rot[0]=cos(yaw)*cos(roll);
    rot[1]=cos(yaw)*sin(roll)*sin(pitch) - sin(yaw)*cos(pitch);
    rot[2]=cos(yaw)*sin(roll)*cos(pitch) + sin(yaw)*sin(pitch);
    rot[3]=sin(yaw)*cos(roll);
    rot[4]=sin(yaw)*sin(roll)*sin(pitch) + cos(yaw)*cos(pitch);
    rot[5]=sin(yaw)*sin(roll)*cos(pitch) - cos(yaw)*sin(pitch);
    rot[6]=sin(-1*roll);
    rot[7]=cos(roll)*sin(pitch);
    rot[8]=cos(roll)*cos(pitch);
    
    auto point_cloud_model_copy=point_cloud_model;
    auto point_cloud_data_copy=point_cloud_data;
    apply_rotation(roll, pitch, yaw, point_cloud_model_copy, point_cloud_data_copy, uav_model, uav_data);
    
    vector<int> matching1(point_cloud_data.size());
    vector<double> err_per_point1(point_cloud_data.size());
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
        DebugOn("compute mid new ub "<<best_ub<<endl);
    }
}
void run_preprocess_parallel_Align(const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data, const vector<vector<double>>& uav_model, const vector<vector<double>>& uav_data, vector<int>& pos_vec, vector<shared_ptr<Model<double>>>& models, const vector<treenode_r>& vec_node, vector<int>& m_vec,vector<double>& vec_lb, vector<indices>& valid_cells, int nb_threads, double upper_bound, double lower_bound, vector<param<double>>& dist_cost_new, int iter, string error_type){
    vector<shared_ptr<Model<double>>> temp_models;
    std::vector<thread> threads;
    int nd=point_cloud_data.size();
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
    m_vec.resize(num, 0);
    vec_lb.resize(num, 0.0);
    temp_models.resize(num);
    
    vector<double> vec_prep_time;
    vec_prep_time.resize(num, 0.0);
    
    for (auto i = 0; i < num; i++) {
        threads.push_back(thread(&run_preprocess_model_Align, ref(point_cloud_model), ref(point_cloud_data), ref(uav_model), ref(uav_data), ref(vec_node[i]), ref(m_vec[i]),ref(vec_lb[i]), ref(valid_cells[i]), ref(vec_dist_cost[i]), ref(vec_prep_time[i]), ref(upper_bound), ref(temp_models[i]), ref(error_type)));
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
    dist_cost_new=vec_dist_cost;
    
}
void run_preprocess_model_Align(const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data, const vector<vector<double>>& uav_model, const vector<vector<double>>& uav_data, treenode_r vec_node_i, int& m_vec_i,  double& vec_lb_i,  indices& valid_cells_i, param<double>& dist_cost_i, double& prep_time_i, double upper_bound, shared_ptr<Model<double>>& model_i, string error_type){
    
    
    vec_lb_i=preprocess_lid(point_cloud_model, point_cloud_data, uav_model, uav_data,vec_node_i.valid_cells, valid_cells_i,  dist_cost_i, vec_node_i.roll.first, vec_node_i.roll.second, vec_node_i.pitch.first, vec_node_i.pitch.second, vec_node_i.yaw.first ,vec_node_i.yaw.second, upper_bound, prep_time_i, error_type);
    //
    bool model_created=false;
    if(valid_cells_i.size()>=point_cloud_data.size()){
        if(error_type=="L2"){
            model_i=Align_L2_model(point_cloud_model, point_cloud_data, uav_model, uav_data,vec_node_i.roll.first, vec_node_i.roll.second, vec_node_i.pitch.first, vec_node_i.pitch.second, vec_node_i.yaw.first ,vec_node_i.yaw.second, valid_cells_i, dist_cost_i);
        }
        else{
            model_i=Align_L1_model(point_cloud_model, point_cloud_data, uav_model, uav_data,vec_node_i.roll.first, vec_node_i.roll.second, vec_node_i.pitch.first, vec_node_i.pitch.second, vec_node_i.yaw.first ,vec_node_i.yaw.second, valid_cells_i, dist_cost_i);
        }
        model_created=true;
    }
    if(model_created){
        m_vec_i=1;
    }
    else{
        m_vec_i=0;
    }
}
#endif /* BB_h */
