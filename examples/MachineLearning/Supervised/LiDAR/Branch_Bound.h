//
//  Branch_Bound.h
//  Gravity
//
//  Created by Smitha on 8/18/21.
//

#ifndef Branch_Bound_h
#define Branch_Bound_h
//
//  Branch_Bound.cpp
//  lidar
//
//  Created by Smitha on 8/18/21.
//
#include <gravity/solver.h>
#ifdef USE_OPT_PARSER
#include <optionParser.hpp>
#endif
#include <gravity/rapidcsv.h>
#ifdef USE_MATPLOT
#include <gravity/matplotlibcpp.h>
#endif
#include <queue>
#include <DataSet.h>
#include "lasreader.hpp"
#include "laswriter.hpp"
#include "spherical.h"
#include <algorithm>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <future>
#include <thread>
#ifdef USE_EIGEN3
#include <Eigen/Dense>
#endif
//#include <tclap/CmdLine.h>
//#include <sqlite3.h>

//#define WITH_CGAL
#ifdef WITH_CGAL
#include "CGALStuff.h"
#endif

#include "convexes.h"
#include <gravity/KDTreeVectorOfVectorsAdaptor.h>
#include <time.h>
#ifdef USE_PCL
#include <pcl/features/fpfh.h>
#include <pcl/point_types.h>
#include <pcl/features/normal_3d.h>
#endif
using namespace std;
#include <gravity/jly_goicp.h>
#include <gravity/ConfigMap.hpp>
using namespace Go_ICP;
#include <stdio.h>
#include "Lidar_utils.h"
#ifdef USE_GJK
extern "C" {
#include "openGJK.h"
}
#endif

pair<double,double> min_max_euclidean_distsq_box(vector<vector<double>> coords, vector<double> center);
bool compute_vertices(vector<vector<double>> vertex_set_a, vector<vector<double>> facets_a, vector<vector<double>> facets_b, vector<vector<int>> vertex_edge_a, vector<vector<pair<int, int>>> vertex_edge_plane_a,   vector<int> feas_set_a, vector<int> infeas_set_a, vector<vector<int>> infeas_facet_set_a, std::vector<std::vector<double>>& new_vert, const std::vector<double>& model_point);
bool vertices_box_plane_reg(const vector<double>& plane_eq, const vector<vector<double>>& big_box, vector<vector<double>>& new_verts, vector<int>& infeas_set);
void get_extreme_rotation_data(vector<vector<double>>& extreme, const vector<double>& d_pt, const vector<var<double>>& theta_vec);

#ifdef USE_GJK
double distance_polytopes_gjk(vector<vector<double>>& vec1, vector<vector<double>>& vec2);
#endif
void add_vertex(std::vector<std::vector<double>>& new_vert, std::vector<double> solution, const std::vector<double>& model_point);

void preprocess_poltyope_ve_gjk_in_centroid(const vector<vector<double>>& point_cloud_data, const vector<vector<double>>& point_cloud_model, const indices& old_cells, param<double>& dist_cost_old, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double shift_min_x, double shift_max_x, double shift_min_y, double shift_max_y, double shift_min_z, double shift_max_z, const vector<vector<vector<double>>>& model_voronoi_normals, const vector<vector<double>>& model_face_intercept, const vector<vector<vector<double>>>& model_voronoi_vertices, param<double>& dist_cost, double upper_bound, double lower_bound, double& min_cost_sum, indices& new_cells,  double& new_tx_min, double& new_tx_max, double& new_ty_min, double& new_ty_max, double& new_tz_min, double& new_tz_max, double& prep_time_total, const vector<vector<pair<double, double>>>& model_voronoi_min_max, const vector<vector<vector<int>>>& model_voronoi_vertex_edge, const vector<vector<vector<pair<int,int>>>>& model_voronoi_vertex_edge_planes, const vector<double> max_vert_vert_dist_sq, const param<double>& dii, const param<double>& djj, param<double>& maxj);/* Run the ARMO model for boresight alignment */
bool vert_enum(vector<vector<double>> vertex_set_a, vector<vector<double>> facets_a, vector<vector<int>> vertex_edge_a, vector<vector<pair<int, int>>> vertex_edge_plane_a, vector<vector<double>> vertex_set_b, vector<vector<double>> facets_b,  vector<vector<int>> vertex_edge_b, vector<vector<pair<int, int>>> vertex_edge_plane_b, std::vector<std::vector<double>>& new_vert, const std::vector<double>& model_point);
shared_ptr<Model<double>> build_linobj_convex_clean(const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data, const indices& valid_cells, double new_roll_min, double new_roll_max, double new_pitch_min, double new_pitch_max, double new_yaw_min, double new_yaw_max, double new_shift_min_x, double new_shift_max_x, double new_shift_min_y, double new_shift_max_y, double new_shift_min_z, double new_shift_max_z, param<>& dist_cost, double ub, bool& status, double lb, param<double>& maxj);
vector<double> BranchBound11(GoICP& goicp, vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, vector<vector<vector<double>>> model_voronoi_normals, vector<vector<double>> model_face_intercept, vector<vector<vector<double>>> model_voronoi_vertices,  const vector<double>& model_inner_prod_min,const vector<double>& model_inner_prod_max, vector<pair<double, double>> min_max_model,  vector<double>& best_rot_trans, double best_ub, const vector<vector<pair<double, double>>>& model_voronoi_min_max, const vector<vector<vector<int>>>& model_voronoi_vertex_edge, const vector<vector<vector<pair<int,int>>>>& model_voronoi_vertex_edge_planes, const vector<double>& max_vert_vert_dist_sq, const param<double>& dii, const param<double>& djj);
void run_preprocess_parallel_new(const vector<vector<double>>& point_cloud_data, const vector<vector<double>>& point_cloud_model, const vector<vector<vector<double>>>& model_voronoi_normals, const vector<vector<double>>& model_face_intercept, const vector<vector<vector<double>>>& model_voronoi_vertices, vector<int>& pos_vec, vector<shared_ptr<Model<double>>>& models, const vector<treenode_p>& vec_node, vector<int>& m_vec,vector<double>& vec_lb, vector<indices>& valid_cells, int nb_threads, double upper_bound, double lower_bound, vector<double>& new_shift_x_min, vector<double>& new_shift_x_max, vector<double>& new_shift_y_min, vector<double>& new_shift_y_max, vector<double>& new_shift_z_min, vector<double>& new_shift_z_max, int max_cells, const vector<vector<pair<double, double>>>& model_voronoi_min_max, const vector<vector<vector<int>>>& model_voronoi_vertex_edge, const vector<vector<vector<pair<int,int>>>>& model_voronoi_vertex_edge_planes, const param<double>& dii, const param<double>& djj, vector<param<double>>& dist_cost_new, int iter, vector<vector<double>>& costs_upto_vec, const vector<double>& max_vert_vert_dist_sq);
void run_preprocess_model(const vector<vector<double>>& point_cloud_data, const vector<vector<double>>& point_cloud_model, const vector<vector<vector<double>>>& model_voronoi_normals, const vector<vector<double>>& model_face_intercept, const vector<vector<vector<double>>>& model_voronoi_vertices, shared_ptr<Model<double>>& model_i, int& m_vec_i, int& pos_vec_i, double& vec_lb_i, treenode_p vec_node_i, indices& valid_cells_i, param<double>& dist_cost_i, int nb_threads, double upper_bound, double lower_bound, double& new_shift_x_min_i, double& new_shift_x_max_i, double& new_shift_y_min_i, double& new_shift_y_max_i, double& new_shift_z_min_i, double& new_shift_z_max_i, double& prep_time_i, int max_cells, const vector<vector<pair<double, double>>>& model_voronoi_min_max, const vector<vector<vector<int>>>& model_voronoi_vertex_edge, const vector<vector<vector<pair<int,int>>>>& model_voronoi_vertex_edge_planes, const param<double>& dii, const param<double>& djj, int iter, const vector<double>& max_vert_vert_dist_sq);
void compute_upper_boundICP(GoICP& goicp, double roll_mini, double roll_maxi, double pitch_mini, double pitch_maxi, double yaw_mini, double yaw_maxi, double shift_min_xi, double shift_max_xi, double shift_min_yi, double shift_max_yi, double shift_min_zi, double shift_max_zi, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double shift_min_x, double shift_max_x, double shift_min_y, double shift_max_y, double shift_min_z, double shift_max_z, vector<double>& best_rot_trans, double& best_ub, vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data);
GoICP Initialize_BB(vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, const vector<pair<double, double>>& min_max_model, vector<pair<double, double>>& min_max_d, double shift_min_x, double shift_max_x, double shift_min_y, double shift_max_y, double shift_min_z, double shift_max_z, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double& best_ub, vector<double>& best_rot_trans);
double run_ICP_only(GoICP& goicp, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double shift_min_x, double shift_max_x, double shift_min_y, double shift_max_y, double shift_min_z, double shift_max_z, vector<double>& rot_trans_ub);
GoICP initialize_ICP_only(const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data);
double run_ICP_only(GoICP& goicp, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double shift_min_x, double shift_max_x, double shift_min_y, double shift_max_y, double shift_min_z, double shift_max_z, vector<double>& rot_trans_ub);
vector<double> BranchBound11(GoICP& goicp, vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, vector<vector<vector<double>>> model_voronoi_normals, vector<vector<double>> model_face_intercept, vector<vector<vector<double>>> model_voronoi_vertices,  const vector<double>& model_inner_prod_min,const vector<double>& model_inner_prod_max, vector<pair<double, double>> min_max_model,  vector<double>& best_rot_trans, double best_ub, const vector<vector<pair<double, double>>>& model_voronoi_min_max, const vector<vector<vector<int>>>& model_voronoi_vertex_edge, const vector<vector<vector<pair<int,int>>>>& model_voronoi_vertex_edge_planes, const vector<double>& max_vert_vert_dist_sq, const param<double>& dii, const param<double>& djj)
{
    /* INPUT BOUNDS */
    
    /* INPUT BOUNDS */
    double time_start = get_wall_time();
    double total_time_max = 900000;
    double prep_time_total=0;
    double shift_min_x = min_max_model[0].first, shift_max_x = min_max_model[0].second, shift_min_y = min_max_model[1].first,shift_max_y = min_max_model[1].second,shift_min_z = min_max_model[2].first,shift_max_z = min_max_model[2].second;
    double yaw_min = -50*pi/180., yaw_max = 50*pi/180., pitch_min =-50*pi/180.,pitch_max = 50*pi/180.,roll_min =-50*pi/180.,roll_max = 50*pi/180.;
    indices N1 = range(1,point_cloud_data.size());
    indices N2 = range(1,point_cloud_model.size());
    int nd=point_cloud_data.size();
    vector<int> new_matching(N1.size());
    bool convex = false;
    double max_time = 10;
    double max_time_init=10;
    bool max_time_increase=false;
    int max_iter = 1e6;
    int models_count=0, models_new_count=0;
    int infeasible_count=0;
    vector<pair<pair<int,int>,pair<int,int>>> incompatible_pairs;
    size_t nb_threads = std::thread::hardware_concurrency();
    //int nb_threads = 1;


    // nb_threads = 1;
    pair<double,double> shift_x_bounds_r, shift_y_bounds_r, shift_z_bounds_r, roll_bounds_r, pitch_bounds_r, yaw_bounds_r;
    shift_x_bounds_r={shift_min_x, shift_max_x};
    shift_y_bounds_r={shift_min_y, shift_max_y};
    shift_z_bounds_r={shift_min_z, shift_max_z};
    roll_bounds_r={roll_min, roll_max};
    pitch_bounds_r={pitch_min, pitch_max};
    yaw_bounds_r={yaw_min, yaw_max};
    vector<pair<double,double>> shift_x_bounds, shift_y_bounds, shift_z_bounds;
    vector<pair<double,double>> roll_bounds, pitch_bounds, yaw_bounds;
    vector<double> new_shift_x_min, new_shift_x_max, new_shift_y_min, new_shift_y_max, new_shift_z_min, new_shift_z_max;
    vector<indices> valid_cells;
    vector<vector<double>> rot_trans(nb_threads, vector<double>(12));
    vector<double> rot_trans_temp(12);
    vector<double> rot_trans_r(12);
   
    vector<double> sol_gur(12);
    vector<int> pos_vec;
    vector<double> vec_lb;
    vector<treenode_p> vec_node;
    vector<int> m_vec;
    vector<double> res(N1.size());
    vector<vector<double>> costs_upto_vec;
    auto point_cloud_data_copy=point_cloud_data;
    DebugOn("I will be using " << nb_threads << " parallel threads" << endl);
    vector<shared_ptr<Model<>>> models, models_new;
    double lb = 0, ub = 12, ub_=-1, best_lb = 0;
    int nb_pruned = 0;
    int depth_r=0, iter=0;
    vector<int> depth_vec, depth_vec_new;
    priority_queue<treenode_p> lb_queue;
    vector<double> costs_upto_init(nd,0.0);
    auto valid_cells_ro=indices(N1,N2);
    vector<double> solution_lb_ones;
    for(auto i=0;i<nd;i++){
        solution_lb_ones.push_back(1);
    }
    double min_cost_sum=0.0;
    double shift_min_x_new, shift_max_x_new, shift_min_y_new, shift_max_y_new, shift_min_z_new, shift_max_z_new;
    vector<pair<double, double>> newmm_r;
    param<double> dist_cost_r("dist_cost_r");
    indices valid_cells_r;
    vector<param<double>> dist_cost_cells;
    double  prep_time=0.0;
    min_cost_sum=0;
    //preprocess_poltyope_cdd_gjk_centroid(point_cloud_data, point_cloud_model, valid_cells_ro, roll_min, roll_max, pitch_min, pitch_max, yaw_min, yaw_max, shift_min_x, shift_max_x, shift_min_y, shift_max_y, shift_min_z, shift_max_z, min_max_each_d, model_voronoi_normals, model_face_intercept, model_voronoi_vertices, dist_cost_r, best_ub,best_lb, min_cost_sum, valid_cells_r,  shift_min_x_new, shift_max_x_new, shift_min_y_new, shift_max_y_new, shift_min_z_new, shift_max_z_new, new_min_max_each_d, prep_time, model_voronoi_min_max);
    param<double> dist_cost_old("dist_cost_old");
    for(auto key:*valid_cells_ro._keys){
        dist_cost_old.add_val(key, 0.0);
    }
    
    //preprocess_poltyope_ve_gjk_centroid(point_cloud_data, point_cloud_model, valid_cells_ro, dist_cost_old, roll_min, roll_max, pitch_min, pitch_max, yaw_min, yaw_max, shift_min_x, shift_max_x, shift_min_y, shift_max_y, shift_min_z, shift_max_z, model_voronoi_normals, model_face_intercept, model_voronoi_vertices,    dist_cost_r, best_ub,best_lb, min_cost_sum, valid_cells_r,  shift_min_x_new, shift_max_x_new, shift_min_y_new, shift_max_y_new, shift_min_z_new, shift_max_z_new, prep_time, model_voronoi_min_max, model_voronoi_vertex_edge, model_voronoi_vertex_edge_planes, dii, djj);
    
    shift_x_bounds_r={shift_min_x_new, shift_max_x_new};
    shift_y_bounds_r={shift_min_y_new, shift_max_y_new};
    shift_z_bounds_r={shift_min_z_new, shift_max_z_new};
    
    lb_queue.push(treenode_p(roll_bounds_r, pitch_bounds_r, yaw_bounds_r, shift_x_bounds_r, shift_y_bounds_r, shift_z_bounds_r, lb, ub, ub_, depth_r, valid_cells_r, false, dist_cost_r));
    double max_incr=0, max_ratio=1;
    int nb_threads_half=nb_threads/2;
    int test_ub=10;
    treenode_p topnode=lb_queue.top();

    int ndisc=0;
    if(ndisc>1){
        lb_queue.pop();
    double sx=(shift_max_x_new-shift_min_x_new)/ndisc;
    double sy=(shift_max_y_new-shift_min_y_new)/ndisc;
    double sz=(shift_max_z_new-shift_min_z_new)/ndisc;

    for (auto i=0;i<ndisc;i++){
        for(auto j=0;j<ndisc;j++){
            for(auto k=0;k<ndisc;k++){
                pair<double,double> shift_x_bounds_i={topnode.tx.first+i*sx, topnode.tx.first+(i+1)*sx};
                pair<double,double> shift_y_bounds_i={topnode.ty.first+j*sy, topnode.ty.first+(j+1)*sy};
                pair<double,double> shift_z_bounds_i={topnode.tz.first+k*sz, topnode.tz.first+(k+1)*sz};
                pair<double,double> roll_bounds_i={topnode.roll.first, topnode.roll.second};
                pair<double,double> pitch_bounds_i={topnode.pitch.first, topnode.pitch.second};
                pair<double,double> yaw_bounds_i={topnode.yaw.first, topnode.yaw.second};
                lb_queue.push(treenode_p(roll_bounds_i, pitch_bounds_i, yaw_bounds_i, shift_x_bounds_i, shift_y_bounds_i, shift_z_bounds_i, lb, ub, ub_, topnode.depth+1, valid_cells_r, false, dist_cost_r));
            }
        }
    }
    
    }
    
    int max_cell_size=150000;
    //param<double> dist_cost_rl("dist_cost_rl");
    // auto valid_cells_rl = preprocess_poltyope_cdd_gjk_bounds(point_cloud_data, point_cloud_model, valid_cells_ro, roll_min, roll_max, pitch_min, pitch_max, yaw_min, yaw_max, shift_min_x, shift_max_x, shift_min_y, shift_max_y, shift_min_z, shift_max_z, model_voronoi_normals, model_face_intercept, model_voronoi_vertices, new_model_pts, new_model_ids, dist_cost_rl, best_ub, nb_threads, min_cost_sum, min_max_d);
    // lb_queue.push(treenode_m(roll_bounds_r, pitch_bounds_r, yaw_bounds_r, shift_x_bounds_r, shift_y_bounds_r, shift_z_bounds_r, lb, best_ub, ub_, depth_r, valid_cells_rl, false));
    //auto goicp=initialize_ICP_only(point_cloud_model, point_cloud_data);
    
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
        // Discard all nodes with high lower bounds in the queue
        //        int old_size = lb_queue.size();
        //        priority_queue<treenode_n> new_lb_queue;
        //        while(!lb_queue.empty())
        //        {
        //            auto node = lb_queue.top();
        //            lb_queue.pop();
        //            if(node.lb < best_ub)
        //                new_lb_queue.push(node);
        //            else
        //                break;
        //        }
        //        lb_queue = new_lb_queue;
        //        if(old_size - lb_queue.size()>0){
        //            DebugOn("Just pruned " << old_size - lb_queue.size() << " node(s)\n");
        //            nb_pruned += old_size - lb_queue.size();
        //        }
        DebugOn("Total infeasible =  " << infeasible_count << endl);
        DebugOn("Total prep_time =  " << prep_time_total << endl);
        DebugOn("Total discarded =  " << prep_count << endl);
        max_incr=0, max_ratio=1;
        pos_vec.clear();
        models.clear();
        shift_x_bounds.clear();
        shift_y_bounds.clear();
        shift_z_bounds.clear();
        roll_bounds.clear();
        pitch_bounds.clear();
        yaw_bounds.clear();
        valid_cells.clear();
        depth_vec.clear();
        m_vec.clear();
        vec_node.clear();
        vec_lb.clear();
        new_shift_x_min.clear();
        new_shift_x_max.clear();
        new_shift_y_min.clear();
        new_shift_y_max.clear();
        new_shift_z_min.clear();
        new_shift_z_max.clear();
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
                double x_shift_increment,  y_shift_increment, z_shift_increment;
                double roll_increment,  pitch_increment, yaw_increment;
                if(false && iter<=100){
                    step=8;
                if(topnode.depth%2==0){
                    if(topnode.tx.first<=-0.1 && topnode.tx.second>=0.1){
                        x_shift_increment = topnode.tx.first*(-1);
                    }
                    else{
                     x_shift_increment = (topnode.tx.second - topnode.tx.first)/2.0;
                    }
                    if(topnode.ty.first<=-0.1 && topnode.ty.second>=0.1){
                        y_shift_increment = topnode.ty.first*(-1);
                    }
                    else{
                     y_shift_increment = (topnode.ty.second - topnode.ty.first)/2.0;
                    }
                    if(topnode.tz.first<=-0.1 && topnode.tz.second>=0.1){
                        z_shift_increment = topnode.tz.first*(-1);
                    }
                    else{
                     z_shift_increment = (topnode.tz.second - topnode.tz.first)/2.0;
                    }
                    shift_x_bounds.push_back({topnode.tx.first, topnode.tx.first+x_shift_increment});
                    shift_x_bounds.push_back({topnode.tx.first, topnode.tx.first+x_shift_increment});
                    shift_x_bounds.push_back({topnode.tx.first, topnode.tx.first+x_shift_increment});
                    shift_x_bounds.push_back({topnode.tx.first, topnode.tx.first+x_shift_increment});
                    shift_x_bounds.push_back({topnode.tx.first+x_shift_increment, topnode.tx.second});
                    shift_x_bounds.push_back({topnode.tx.first+x_shift_increment, topnode.tx.second});
                    shift_x_bounds.push_back({topnode.tx.first+x_shift_increment, topnode.tx.second});
                    shift_x_bounds.push_back({topnode.tx.first+x_shift_increment, topnode.tx.second});
                    shift_y_bounds.push_back({topnode.ty.first, topnode.ty.first+y_shift_increment});
                    shift_y_bounds.push_back({topnode.ty.first, topnode.ty.first+y_shift_increment});
                    shift_y_bounds.push_back({topnode.ty.first+y_shift_increment, topnode.ty.second});
                    shift_y_bounds.push_back({topnode.ty.first+y_shift_increment, topnode.ty.second});
                    shift_y_bounds.push_back({topnode.ty.first, topnode.ty.first+y_shift_increment});
                    shift_y_bounds.push_back({topnode.ty.first, topnode.ty.first+y_shift_increment});
                    shift_y_bounds.push_back({topnode.ty.first+y_shift_increment, topnode.ty.second});
                    shift_y_bounds.push_back({topnode.ty.first+y_shift_increment, topnode.ty.second});
                    shift_z_bounds.push_back({topnode.tz.first, topnode.tz.first+z_shift_increment});
                    shift_z_bounds.push_back({topnode.tz.first+z_shift_increment, topnode.tz.second});
                    shift_z_bounds.push_back({topnode.tz.first, topnode.tz.first+z_shift_increment});
                    shift_z_bounds.push_back({topnode.tz.first+z_shift_increment, topnode.tz.second});
                    shift_z_bounds.push_back({topnode.tz.first, topnode.tz.first+z_shift_increment});
                    shift_z_bounds.push_back({topnode.tz.first+z_shift_increment, topnode.tz.second});
                    shift_z_bounds.push_back({topnode.tz.first, topnode.tz.first+z_shift_increment});
                    shift_z_bounds.push_back({topnode.tz.first+z_shift_increment, topnode.tz.second});
                    for(auto k=0;k<8;k++){
                        roll_bounds.push_back({topnode.roll.first, topnode.roll.second});
                        pitch_bounds.push_back({topnode.pitch.first, topnode.pitch.second});
                        yaw_bounds.push_back({topnode.yaw.first, topnode.yaw.second});
                        vec_node.push_back(treenode_p(roll_bounds[i+k],  pitch_bounds[i+k], yaw_bounds[i+k], shift_x_bounds[i+k], shift_y_bounds[i+k], shift_z_bounds[i+k], topnode.lb, best_ub, -1.0, topnode.depth+1, topnode.valid_cells, false,topnode.dist_cost_cells));
                        depth_vec.push_back(topnode.depth+1);
                    }
                }
                else{
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
                        shift_x_bounds.push_back({topnode.tx.first, topnode.tx.second});
                        shift_y_bounds.push_back({topnode.ty.first, topnode.ty.second});
                        shift_z_bounds.push_back({topnode.tz.first, topnode.tz.second});
                        vec_node.push_back(treenode_p(roll_bounds[i+k],  pitch_bounds[i+k], yaw_bounds[i+k], shift_x_bounds[i+k], shift_y_bounds[i+k], shift_z_bounds[i+k], topnode.lb, best_ub, -1.0, topnode.depth+1, topnode.valid_cells, false,topnode.dist_cost_cells));
                        depth_vec.push_back(topnode.depth+1);
                    }
                }
                }
                else{
                    step=2;
                    double x_shift_increment = (topnode.tx.second - topnode.tx.first)/2.0;
                                     max_incr = x_shift_increment;
                                     double y_shift_increment = (topnode.ty.second - topnode.ty.first)/2.0;
                                     if(y_shift_increment>max_incr)
                                         max_incr = y_shift_increment;
                                     double z_shift_increment = (topnode.tz.second - topnode.tz.first)/2.0;
                                     if(z_shift_increment>max_incr)
                                         max_incr = z_shift_increment;
                                     double roll_increment = (topnode.roll.second - topnode.roll.first)/2.0;
                                     if(roll_increment>max_incr)
                                         max_incr = roll_increment;
                                     double pitch_increment = (topnode.pitch.second - topnode.pitch.first)/2.0;
                                     if(pitch_increment>max_incr)
                                         max_incr = pitch_increment;
                                     double yaw_increment = (topnode.yaw.second - topnode.yaw.first)/2.0;
                                     if(yaw_increment>max_incr)
                                         max_incr = yaw_increment;
                                     shift_x_bounds.push_back({topnode.tx.first, topnode.tx.second});
                                     shift_y_bounds.push_back({topnode.ty.first, topnode.ty.second});
                                     shift_z_bounds.push_back({topnode.tz.first, topnode.tz.second});
                                     roll_bounds.push_back({topnode.roll.first, topnode.roll.second});
                                     pitch_bounds.push_back({topnode.pitch.first, topnode.pitch.second});
                                     yaw_bounds.push_back({topnode.yaw.first, topnode.yaw.second});
                                     depth_vec.push_back(topnode.depth+1);
                                     shift_x_bounds.push_back({topnode.tx.first, topnode.tx.second});
                                     shift_y_bounds.push_back({topnode.ty.first, topnode.ty.second});
                                     shift_z_bounds.push_back({topnode.tz.first, topnode.tz.second});
                                     roll_bounds.push_back({topnode.roll.first, topnode.roll.second});
                                     pitch_bounds.push_back({topnode.pitch.first, topnode.pitch.second});
                                     yaw_bounds.push_back({topnode.yaw.first, topnode.yaw.second});
                                     depth_vec.push_back(topnode.depth+1);
                                     if(max_incr==x_shift_increment){
                                         if(topnode.tx.first<=-0.1 && topnode.tx.second>=0.1){
                                             shift_x_bounds[i] = {topnode.tx.first, 0};
                                             shift_x_bounds[i+1] = {0, topnode.tx.second};
                                         }
                                         else{
                                             shift_x_bounds[i] = {topnode.tx.first, topnode.tx.first+x_shift_increment};
                                             shift_x_bounds[i+1] = {topnode.tx.first+x_shift_increment, topnode.tx.second};
                                         }
                                         
                                     }
                                     else if(max_incr==y_shift_increment){
                                         if(topnode.ty.first<=-0.1 && topnode.ty.second>=0.1){
                                             shift_y_bounds[i] = {topnode.ty.first, 0};
                                             shift_y_bounds[i+1] = {0, topnode.ty.second};
                                         }
                                         else{
                                             shift_y_bounds[i] = {topnode.ty.first, topnode.ty.first+y_shift_increment};
                                             shift_y_bounds[i+1] = {topnode.ty.first+y_shift_increment, topnode.ty.second};
                                         }
                                         
                                     }
                                     else if(max_incr==z_shift_increment){
                                         if(topnode.tz.first<=-0.1 && topnode.tz.second>=0.1){
                                             shift_z_bounds[i] = {topnode.tz.first, 0};
                                             shift_z_bounds[i+1] = {0, topnode.tz.second};
                                         }
                                         else{
                                             shift_z_bounds[i] = {topnode.tz.first, topnode.tz.first+z_shift_increment};
                                             shift_z_bounds[i+1] = {topnode.tz.first+z_shift_increment, topnode.tz.second};
                                         }
                                        
                                     }
                                     else if(max_incr==roll_increment){
                                         roll_bounds[i] = {topnode.roll.first, topnode.roll.first+roll_increment};
                                         roll_bounds[i+1] = {topnode.roll.first+roll_increment, topnode.roll.second};
                                         
                                     }
                                     else if(max_incr==pitch_increment){
                                         pitch_bounds[i] = {topnode.pitch.first, topnode.pitch.first+pitch_increment};
                                         pitch_bounds[i+1] = {topnode.pitch.first+pitch_increment, topnode.pitch.second};
                                        
                                     }
                                     else if(max_incr==yaw_increment){
                                         yaw_bounds[i] = {topnode.yaw.first, topnode.yaw.first+yaw_increment};
                                         yaw_bounds[i+1] = {topnode.yaw.first+yaw_increment, topnode.yaw.second};
                                        
                                     }
                                     if(max_incr==0){
                                         DebugOn("Warn: Identical child nodes"<<endl);
                                     }
                                     else{
                                         DebugOff("max_incr "<<max_incr<<endl);
                                     }
                                 
                                 
                                 vec_node.push_back(treenode_p(roll_bounds[i],  pitch_bounds[i], yaw_bounds[i], shift_x_bounds[i], shift_y_bounds[i], shift_z_bounds[i], topnode.lb, best_ub, -1.0, topnode.depth+1, topnode.valid_cells, false,topnode.dist_cost_cells));
                                 vec_node.push_back(treenode_p(roll_bounds[i+1],  pitch_bounds[i+1], yaw_bounds[i+1], shift_x_bounds[i+1], shift_y_bounds[i+1], shift_z_bounds[i+1], topnode.lb, best_ub, -1.0, topnode.depth+1, topnode.valid_cells, false, topnode.dist_cost_cells));
                }
              
            
      
                auto ut1=get_wall_time();
                if (i%test_ub==0){
                    compute_upper_boundICP(goicp, roll_bounds[i].first, roll_bounds[i].second, pitch_bounds[i].first, pitch_bounds[i].second, yaw_bounds[i].first, yaw_bounds[i].second, shift_x_bounds[i].first, shift_x_bounds[i].second, shift_y_bounds[i].first, shift_y_bounds[i].second, shift_z_bounds[i].first, shift_z_bounds[i].second, roll_min, roll_max, pitch_min, pitch_max, yaw_min, yaw_max, shift_min_x, shift_max_x, shift_min_y, shift_max_y, shift_min_z, shift_max_z, best_rot_trans, best_ub, point_cloud_model, point_cloud_data);
                    compute_upper_boundICP(goicp, roll_bounds[i+1].first, roll_bounds[i+1].second, pitch_bounds[i+1].first, pitch_bounds[i+1].second, yaw_bounds[i+1].first, yaw_bounds[i+1].second, shift_x_bounds[i+1].first, shift_x_bounds[i+1].second, shift_y_bounds[i+1].first, shift_y_bounds[i+1].second, shift_z_bounds[i+1].first, shift_z_bounds[i+1].second, roll_min, roll_max, pitch_min, pitch_max, yaw_min, yaw_max, shift_min_x, shift_max_x, shift_min_y, shift_max_y, shift_min_z, shift_max_z, best_rot_trans, best_ub, point_cloud_model, point_cloud_data);
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
      
        run_preprocess_parallel_new(point_cloud_data, point_cloud_model, model_voronoi_normals, model_face_intercept, model_voronoi_vertices, pos_vec, models, vec_node, m_vec, vec_lb, valid_cells, nb_threads, best_ub, best_lb, new_shift_x_min, new_shift_x_max, new_shift_y_min, new_shift_y_max, new_shift_z_min, new_shift_z_max, max_cell_size, model_voronoi_min_max, model_voronoi_vertex_edge,  model_voronoi_vertex_edge_planes, dii,  djj, dist_cost_cells, iter, costs_upto_vec, max_vert_vert_dist_sq);

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
                if(ub_>=0 && (ub_-lb_/ub_)<=1e-6){
                    if(ub_<=best_ub-1e-4){
                        //models[j]->print_solution();
                        //models[j]->print();
                        //models[j]->print_constraints_stats(1e-4);
                        vector<double> rot_trans(12);
                        bool is_rotation  = get_solution(models[j], rot_trans, new_matching);
                        if(is_rotation){
                            point_cloud_data_copy=point_cloud_data;
                            apply_rot_trans(rot_trans, point_cloud_data_copy);
                            auto L2err=computeL2error(point_cloud_model, point_cloud_data_copy, new_matching, res);
                            if(L2err<=best_ub){
                                best_ub=L2err;
                                best_rot_trans=rot_trans;
                                DebugOn("new best ub "<<best_ub<<" ub_ "<<ub_<<" lb_ "<<lb_<<endl);
                                leaf_node=true;
                                DebugOn("leaf lb "<<lb_<<" L2 "<<L2err<<" ub_ "<<ub_<<endl);
                                //leaf_node=true;
                                //lb=L2err;
                            }
                        }
                    }
                }
                if(true){
                    lb = std::max(models[j]->get_rel_obj_val(), vec_lb[pos]);
                }
                if(lb-1e-4<=best_ub)
                {
                    shift_x_bounds_r={new_shift_x_min[pos], new_shift_x_max[pos]};
                    shift_y_bounds_r={new_shift_y_min[pos], new_shift_y_max[pos]};
                    shift_z_bounds_r={new_shift_z_min[pos], new_shift_z_max[pos]};
                    lb_queue.push(treenode_p(roll_bounds[pos], pitch_bounds[pos], yaw_bounds[pos], shift_x_bounds_r, shift_y_bounds_r, shift_z_bounds_r, lb, best_ub, ub_, depth_vec[pos], valid_cells[pos], leaf_node, dist_cost_cells[pos]));
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
    DebugOn("tx "<<best_rot_trans[9]<<endl);
    DebugOn("ty "<<best_rot_trans[10]<<endl);
    DebugOn("tz "<<best_rot_trans[11]<<endl);
    while(!lb_queue.empty())
    {
        auto node = lb_queue.top();
        DebugOn("node lb "<<node.lb<<" node.leaf "<<node.leaf<<endl);
        DebugOn(node.tx.first<<" "<< node.tx.second<<" "<<node.ty.first<<" "<<node.ty.second<<" "<<node.tz.first<<" "<<node.tz.second<<endl);
        DebugOn(node.roll.first<<" "<< node.roll.second<<" "<<node.pitch.first<<" "<<node.pitch.second<<" "<<node.yaw.first<<" "<<node.yaw.second<<endl);
        compute_upper_boundICP(goicp, node.roll.first, node.roll.second, node.pitch.first, node.pitch.second, node.yaw.first, node.yaw.second, node.tx.first, node.tx.second, node.ty.first, node.ty.second, node.tz.first, node.tz.second, node.roll.first, node.roll.second, node.pitch.first, node.pitch.second, node.yaw.first, node.yaw.second, node.tx.first, node.tx.second, node.ty.first, node.ty.second, node.tz.first, node.tz.second, best_rot_trans, best_ub, point_cloud_model, point_cloud_data);
        auto pitchr = atan2(best_rot_trans[7], best_rot_trans[8]);
        auto rollr = atan2(-best_rot_trans[6], std::sqrt(best_rot_trans[7]*best_rot_trans[7]+best_rot_trans[8]*best_rot_trans[8]));
        auto yawr = atan2(best_rot_trans[3],best_rot_trans[0]);
        if(node.roll.first-1e-3<=rollr && rollr<=node.roll.second+1e-3 && node.pitch.first-1e-3<=pitchr && pitchr<=node.pitch.second+1e-3 && node.yaw.first-1e-3<=yawr && yawr<=node.yaw.second+1e-3 && node.tx.first-1e-3<=best_rot_trans[9] && best_rot_trans[9]<=node.tx.second+1e-3 && node.ty.first-1e-3<=best_rot_trans[10] && best_rot_trans[10]<=node.ty.second+1e-3 && node.tz.first-1e-3<=best_rot_trans[11] && best_rot_trans[11]<=node.tz.second+1e-3){
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
    DebugOn("tx "<<best_rot_trans[9]<<endl);
    DebugOn("ty "<<best_rot_trans[10]<<endl);
    DebugOn("tz "<<best_rot_trans[11]<<endl);
    apply_rot_trans(best_rot_trans, point_cloud_data);
    auto L2err=computeL2error(point_cloud_model, point_cloud_data, new_matching, res);
    DebugOn("L2err "<<L2err<<endl);
    DebugOn("Matching"<<endl);
    for(auto i=0;i<new_matching.size();i++){
        DebugOn(new_matching[i]<<endl);
    }
    auto pitch = atan2(best_rot_trans[7], best_rot_trans[8])*180/pi;
    auto roll = atan2(-best_rot_trans[6], std::sqrt(best_rot_trans[7]*best_rot_trans[7]+best_rot_trans[8]*best_rot_trans[8]))*180/pi;
    auto yaw = atan2(best_rot_trans[3],best_rot_trans[0])*180/pi;
    
    DebugOn("roll "<< roll<<endl);
    DebugOn("pitch "<< pitch<<endl);
    DebugOn("yaw "<< yaw<<endl);
    DebugOn("tx "<<best_rot_trans[9]<<endl);
    DebugOn("ty "<<best_rot_trans[10]<<endl);
    DebugOn("tz "<<best_rot_trans[11]<<endl);
    return best_rot_trans;
    
}


void run_preprocess_parallel_new(const vector<vector<double>>& point_cloud_data, const vector<vector<double>>& point_cloud_model, const vector<vector<vector<double>>>& model_voronoi_normals, const vector<vector<double>>& model_face_intercept, const vector<vector<vector<double>>>& model_voronoi_vertices, vector<int>& pos_vec, vector<shared_ptr<Model<double>>>& models, const vector<treenode_p>& vec_node, vector<int>& m_vec,vector<double>& vec_lb, vector<indices>& valid_cells, int nb_threads, double upper_bound, double lower_bound, vector<double>& new_shift_x_min, vector<double>& new_shift_x_max, vector<double>& new_shift_y_min, vector<double>& new_shift_y_max, vector<double>& new_shift_z_min, vector<double>& new_shift_z_max, int max_cells, const vector<vector<pair<double, double>>>& model_voronoi_min_max, const vector<vector<vector<int>>>& model_voronoi_vertex_edge, const vector<vector<vector<pair<int,int>>>>& model_voronoi_vertex_edge_planes, const param<double>& dii, const param<double>& djj, vector<param<double>>& dist_cost_new, int iter, vector<vector<double>>& costs_upto_vec, const vector<double>& max_vert_vert_dist_sq){
    std::vector<thread> threads;

    vector<shared_ptr<Model<double>>> temp_models;
    
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
    
    vector<double> vec_prep_time;
    
    
    valid_cells.resize(num);
    m_vec.resize(num);
    vec_lb.resize(num);
    temp_models.resize(num);
    
    
    vec_prep_time.resize(num);
    new_shift_x_min.resize(num);
    new_shift_y_min.resize(num);
    new_shift_z_min.resize(num);
    new_shift_x_max.resize(num);
    new_shift_y_max.resize(num);
    new_shift_z_max.resize(num);

    for (auto i = 0; i < num; i++) {
        threads.push_back(thread(&run_preprocess_model, ref(point_cloud_data), ref(point_cloud_model), ref(model_voronoi_normals), ref(model_face_intercept), ref(model_voronoi_vertices), ref(temp_models[i]), ref(m_vec[i]), ref(pos_vec[i]), ref(vec_lb[i]), vec_node[i], ref(valid_cells[i]), ref(vec_dist_cost[i]), nb_threads, upper_bound, lower_bound, ref(new_shift_x_min[i]), ref(new_shift_x_max[i]), ref(new_shift_y_min[i]), ref(new_shift_y_max[i]), ref(new_shift_z_min[i]), ref(new_shift_z_max[i]), ref(vec_prep_time[i]), max_cells, ref(model_voronoi_min_max), ref(model_voronoi_vertex_edge), ref(model_voronoi_vertex_edge_planes), ref(dii), ref(djj), iter, ref(max_vert_vert_dist_sq)));
    }

    /* Join the threads with the main thread */
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


void run_preprocess_model(const vector<vector<double>>& point_cloud_data, const vector<vector<double>>& point_cloud_model, const vector<vector<vector<double>>>& model_voronoi_normals, const vector<vector<double>>& model_face_intercept, const vector<vector<vector<double>>>& model_voronoi_vertices, shared_ptr<Model<double>>& model_i, int& m_vec_i, int& pos_vec_i, double& vec_lb_i, treenode_p vec_node_i, indices& valid_cells_i, param<double>& dist_cost_i, int nb_threads, double upper_bound, double lower_bound, double& new_shift_x_min_i, double& new_shift_x_max_i, double& new_shift_y_min_i, double& new_shift_y_max_i, double& new_shift_z_min_i, double& new_shift_z_max_i, double& prep_time_i, int max_cells, const vector<vector<pair<double, double>>>& model_voronoi_min_max, const vector<vector<vector<int>>>& model_voronoi_vertex_edge, const vector<vector<vector<pair<int,int>>>>& model_voronoi_vertex_edge_planes, const param<double>& dii, const param<double>& djj, int iter, const vector<double>& max_vert_vert_dist_sq){
    int nd=point_cloud_data.size();
    param<double> maxj("maxj");
       preprocess_poltyope_ve_gjk_in_centroid(ref(point_cloud_data), ref(point_cloud_model), ref(vec_node_i.valid_cells),ref(vec_node_i.dist_cost_cells), vec_node_i.roll.first,vec_node_i.roll.second, vec_node_i.pitch.first,vec_node_i.pitch.second, vec_node_i.yaw.first,vec_node_i.yaw.second, vec_node_i.tx.first,vec_node_i.tx.second,vec_node_i.ty.first,vec_node_i.ty.second,vec_node_i.tz.first,vec_node_i.tz.second, ref(model_voronoi_normals), ref(model_face_intercept), ref(model_voronoi_vertices),   ref(dist_cost_i), upper_bound, lower_bound,  ref(vec_lb_i), ref(valid_cells_i),ref(new_shift_x_min_i), ref(new_shift_x_max_i),  ref(new_shift_y_min_i),  ref(new_shift_y_max_i), ref(new_shift_z_min_i), ref(new_shift_z_max_i), ref(prep_time_i), ref(model_voronoi_min_max), ref(model_voronoi_vertex_edge), ref(model_voronoi_vertex_edge_planes),ref(max_vert_vert_dist_sq), ref(dii), ref(djj), ref(maxj));

        DebugOn("v size "<<valid_cells_i.size()<<endl);
//    if(costs_upto.size()<nd && vec_node_i.depth>0 && vec_node_i.depth%1==0){
//        costs_upto.push_back(costs_upto.back());
//    }
 
    bool model_created=false;
        if(valid_cells_i.size()>=nd){
            model_i = build_linobj_convex_clean(point_cloud_model, point_cloud_data, valid_cells_i, vec_node_i.roll.first, vec_node_i.roll.second, vec_node_i.pitch.first, vec_node_i.pitch.second, vec_node_i.yaw.first, vec_node_i.yaw.second, new_shift_x_min_i, new_shift_x_max_i, new_shift_y_min_i, new_shift_y_max_i, new_shift_z_min_i, new_shift_z_max_i, dist_cost_i, upper_bound, model_created, vec_node_i.lb, ref(maxj));
            if(!model_created){
                indices valid_cells_empty("valid_cells_empty");
                valid_cells_i=valid_cells_empty;
            }
        }
    if(model_created){
        m_vec_i=1;
    }
    else{
        m_vec_i=0;
    }
}
void compute_upper_boundICP(GoICP& goicp, double roll_mini, double roll_maxi, double pitch_mini, double pitch_maxi, double yaw_mini, double yaw_maxi, double shift_min_xi, double shift_max_xi, double shift_min_yi, double shift_max_yi, double shift_min_zi, double shift_max_zi, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double shift_min_x, double shift_max_x, double shift_min_y, double shift_max_y, double shift_min_z, double shift_max_z, vector<double>& best_rot_trans, double& best_ub, vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data){
    vector<double> rot_trans_ub;
    vector<int> new_matching(point_cloud_data.size());
    vector<double> res(point_cloud_data.size());
    double ub;
    rot_trans_ub.resize(12);
    vector<vector<double>> point_cloud_data_copy;
    ub= run_ICP_only(goicp, roll_mini, roll_maxi,  pitch_mini, pitch_maxi, yaw_mini, yaw_maxi, shift_min_xi, shift_max_xi, shift_min_yi, shift_max_yi,shift_min_zi, shift_max_zi, rot_trans_ub);
    if(ub<=best_ub+1e-2){
        DebugOff("best ub "<<best_ub<<" ub "<<ub<<endl);
        DebugOff(rot_trans_ub[0]<<" "<<rot_trans_ub[1]<<" "<<rot_trans_ub[2]<<endl);
        DebugOff(rot_trans_ub[3]<<" "<<rot_trans_ub[4]<<" "<<rot_trans_ub[5]<<endl);
        DebugOff(rot_trans_ub[6]<<" "<<rot_trans_ub[7]<<" "<<rot_trans_ub[8]<<endl);
        DebugOff(rot_trans_ub[9]<<" "<<rot_trans_ub[10]<<" "<<rot_trans_ub[11]<<endl);
        auto deti=rot_trans_ub[0]*(rot_trans_ub[4]*rot_trans_ub[8]-rot_trans_ub[7]*rot_trans_ub[5]);
        deti-=rot_trans_ub[1]*(rot_trans_ub[3]*rot_trans_ub[8]-rot_trans_ub[6]*rot_trans_ub[5]);
        deti+=rot_trans_ub[2]*(rot_trans_ub[3]*rot_trans_ub[7]-rot_trans_ub[6]*rot_trans_ub[4]);
        DebugOff("det "<<deti<<endl);
        auto pitch_rad2 = atan2(rot_trans_ub[7], rot_trans_ub[8]);
        auto roll_rad2 = atan2(rot_trans_ub[6], std::sqrt(rot_trans_ub[7]*rot_trans_ub[7]+rot_trans_ub[8]*rot_trans_ub[8]));
        auto yaw_rad2 = atan2(rot_trans_ub[3],rot_trans_ub[0]);
        if(pitch_rad2>=pitch_min && pitch_rad2<=pitch_max && roll_rad2>=roll_min && roll_rad2<=roll_max && yaw_rad2>=yaw_min && yaw_rad2<=yaw_max){
            if(rot_trans_ub[9]>=shift_min_x && rot_trans_ub[9]<=shift_max_x && rot_trans_ub[10]>=shift_min_y && rot_trans_ub[10]<=shift_max_y && rot_trans_ub[11]>=shift_min_z && rot_trans_ub[11]<=shift_max_z ){
                point_cloud_data_copy=point_cloud_data;
                apply_rot_trans(rot_trans_ub, point_cloud_data_copy);
                auto L2erri=computeL2error(point_cloud_model, point_cloud_data_copy, new_matching, res);
                if(L2erri<=best_ub){
                    best_ub=L2erri;
                    best_rot_trans=rot_trans_ub;
                    DebugOn("best ub new "<<best_ub<<" ub "<<ub<<endl);
                    DebugOn("roll "<<roll_rad2<<endl);
                    DebugOn("pitch "<<pitch_rad2<<endl);
                    DebugOn("yaw "<<yaw_rad2<<endl);
                    DebugOn("tx "<<rot_trans_ub[9]<<endl);
                    DebugOn("ty "<<rot_trans_ub[10]<<endl);
                    DebugOn("tz "<<rot_trans_ub[11]<<endl);
                   
                }
            }
        }
        
        //bool is_rotation  = get_solution(models[pos], sol_gur, new_matching);
    }
}

GoICP Initialize_BB(vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, const vector<pair<double, double>>& min_max_model, vector<pair<double, double>>& min_max_d, double shift_min_x, double shift_max_x, double shift_min_y, double shift_max_y, double shift_min_z, double shift_max_z, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double& best_ub, vector<double>& best_rot_trans)
{
    /* INPUT BOUNDS */
    double time_start = get_wall_time();
    double total_time_max = 900000;
    double prep_time_total=0;
    //    double shift_min_x =  -0.5, shift_max_x = 0.5, shift_min_y = -0.5,shift_max_y = 0.5,shift_min_z = -0.5,shift_max_z = 0.5;
    //    double yaw_min = -50*pi/180., yaw_max = 50*pi/180., pitch_min =-50*pi/180.,pitch_max = 50*pi/180.,roll_min =-50*pi/180.,roll_max = 50*pi/180.;
    indices N1 = range(1,point_cloud_data.size());
    indices N2 = range(1,point_cloud_model.size());
    int nd=point_cloud_data.size();
    vector<int> new_matching(N1.size());
    bool convex = false;
    double max_time = 10;
    double max_time_init=10;
    bool max_time_increase=false;
    int max_iter = 1e6;
    auto nb_threads = std::thread::hardware_concurrency();
    min_max_d.resize(3);
    min_max_d[0].first=-1;
    min_max_d[0].second=1;
    min_max_d[1].first=-1;
    min_max_d[1].second=1;
    min_max_d[2].first=-1;
    min_max_d[2].second=1;
    // nb_threads = 1;
    pair<double,double> shift_x_bounds_r, shift_y_bounds_r, shift_z_bounds_r, roll_bounds_r, pitch_bounds_r, yaw_bounds_r;
    shift_x_bounds_r={shift_min_x, shift_max_x};
    shift_y_bounds_r={shift_min_y, shift_max_y};
    shift_z_bounds_r={shift_min_z, shift_max_z};
    roll_bounds_r={roll_min, roll_max};
    pitch_bounds_r={pitch_min, pitch_max};
    yaw_bounds_r={yaw_min, yaw_max};
    vector<pair<double,double>> shift_x_bounds, shift_y_bounds, shift_z_bounds;
    vector<pair<double,double>> roll_bounds, pitch_bounds, yaw_bounds;
    vector<indices> valid_cells;
    vector<vector<double>> rot_trans(nb_threads, vector<double>(12));
    vector<double> rot_trans_temp(12);
    vector<double> rot_trans_r(12);
    best_rot_trans.resize(12);
    vector<double> sol_gur(12);
    vector<int> pos_vec(nb_threads), pos_vec_new;
    vector<double> res(N1.size());
    auto point_cloud_data_copy=point_cloud_data;
    DebugOn("I will be using " << nb_threads << " parallel threads" << endl);
    vector<shared_ptr<Model<>>> models, models_new;
    double lb = 0, ub = 12, ub_=-1, best_lb = 0;
    best_ub = 12;
    int nb_pruned = 0;
    int depth_r=0, iter=0;
    vector<int> depth_vec, depth_vec_new;
    priority_queue<treenode_m> lb_queue;
    auto valid_cells_ro=indices(N1,N2);
    vector<double> solution_lb_ones;
    for(auto i=0;i<nd;i++){
        solution_lb_ones.push_back(1);
    }
    auto goicp=initialize_ICP_only(point_cloud_model, point_cloud_data);
    double min_cost_sum=0.0;
    param<double> dist_cost_r("dist_cost_r");
    
    lb_queue.push(treenode_m(roll_bounds_r, pitch_bounds_r, yaw_bounds_r, shift_x_bounds_r, shift_y_bounds_r, shift_z_bounds_r, lb, ub, ub_, depth_r, valid_cells_ro, false));
    double max_incr=0, max_ratio=1;
    int nb_threads_half=nb_threads/2;
    int imax=nb_threads_half-2;
    int test_ub=1;
    treenode_m topnode=lb_queue.top();
    for(auto i=0;i<10000;i+=2){
        DebugOff("entered loop "<<i<<endl);
        topnode=lb_queue.top();
        lb_queue.pop();
        double x_shift_increment = (topnode.tx.second - topnode.tx.first)/2.0;
        max_incr = x_shift_increment;
        double y_shift_increment = (topnode.ty.second - topnode.ty.first)/2.0;
        if(y_shift_increment>max_incr)
            max_incr = y_shift_increment;
        double z_shift_increment = (topnode.tz.second - topnode.tz.first)/2.0;
        if(z_shift_increment>max_incr)
            max_incr = z_shift_increment;
        double roll_increment = (topnode.roll.second - topnode.roll.first)/2.0;
        if(roll_increment>max_incr)
            max_incr = roll_increment;
        double pitch_increment = (topnode.pitch.second - topnode.pitch.first)/2.0;
        if(pitch_increment>max_incr)
            max_incr = pitch_increment;
        double yaw_increment = (topnode.yaw.second - topnode.yaw.first)/2.0;
        if(yaw_increment>max_incr)
            max_incr = yaw_increment;
        shift_x_bounds.push_back({topnode.tx.first, topnode.tx.second});
        shift_y_bounds.push_back({topnode.ty.first, topnode.ty.second});
        shift_z_bounds.push_back({topnode.tz.first, topnode.tz.second});
        roll_bounds.push_back({topnode.roll.first, topnode.roll.second});
        pitch_bounds.push_back({topnode.pitch.first, topnode.pitch.second});
        yaw_bounds.push_back({topnode.yaw.first, topnode.yaw.second});
        shift_x_bounds.push_back({topnode.tx.first, topnode.tx.second});
        shift_y_bounds.push_back({topnode.ty.first, topnode.ty.second});
        shift_z_bounds.push_back({topnode.tz.first, topnode.tz.second});
        roll_bounds.push_back({topnode.roll.first, topnode.roll.second});
        pitch_bounds.push_back({topnode.pitch.first, topnode.pitch.second});
        yaw_bounds.push_back({topnode.yaw.first, topnode.yaw.second});
        if(max_incr==x_shift_increment){
            shift_x_bounds[i] = {topnode.tx.first, topnode.tx.first+x_shift_increment};
            shift_x_bounds[i+1] = {topnode.tx.first+x_shift_increment, topnode.tx.second};
        }
        else if(max_incr==y_shift_increment){
            shift_y_bounds[i] = {topnode.ty.first, topnode.ty.first+y_shift_increment};
            shift_y_bounds[i+1] = {topnode.ty.first+y_shift_increment, topnode.ty.second};
        }
        else if(max_incr==z_shift_increment){
            shift_z_bounds[i] = {topnode.tz.first, topnode.tz.first+z_shift_increment};
            shift_z_bounds[i+1] = {topnode.tz.first+z_shift_increment, topnode.tz.second};
        }
        else if(max_incr==roll_increment){
            roll_bounds[i] = {topnode.roll.first, topnode.roll.first+roll_increment};
            roll_bounds[i+1] = {topnode.roll.first+roll_increment, topnode.roll.second};
        }
        else if(max_incr==pitch_increment){
            pitch_bounds[i] = {topnode.pitch.first, topnode.pitch.first+pitch_increment};
            pitch_bounds[i+1] = {topnode.pitch.first+pitch_increment, topnode.pitch.second};
        }
        else if(max_incr==yaw_increment){
            yaw_bounds[i] = {topnode.yaw.first, topnode.yaw.first+yaw_increment};
            yaw_bounds[i+1] = {topnode.yaw.first+yaw_increment, topnode.yaw.second};
        }
        if(max_incr==0){
            throw invalid_argument("Warn: Identical child nodes");
        }
        else{
            DebugOff("max_incr "<<max_incr<<endl);
        }
        lb_queue.push(treenode_m(roll_bounds[i], pitch_bounds[i], yaw_bounds[i], shift_x_bounds[i], shift_y_bounds[i], shift_z_bounds[i], lb, ub, ub_, topnode.depth+1, valid_cells_ro, false));
        lb_queue.push(treenode_m(roll_bounds[i+1], pitch_bounds[i+1], yaw_bounds[i+1], shift_x_bounds[i+1], shift_y_bounds[i+1], shift_z_bounds[i+1], lb, ub, ub_, topnode.depth+1, valid_cells_ro, false));
        compute_upper_boundICP(goicp, roll_bounds[i].first, roll_bounds[i].second, pitch_bounds[i].first, pitch_bounds[i].second, yaw_bounds[i].first, yaw_bounds[i].second, shift_x_bounds[i].first, shift_x_bounds[i].second, shift_y_bounds[i].first, shift_y_bounds[i].second, shift_z_bounds[i].first, shift_z_bounds[i].second, roll_min, roll_max, pitch_min, pitch_max, yaw_min, yaw_max, shift_min_x, shift_max_x, shift_min_y, shift_max_y, shift_min_z, shift_max_z, best_rot_trans, best_ub, point_cloud_model, point_cloud_data);
        compute_upper_boundICP(goicp, roll_bounds[i+1].first, roll_bounds[i+1].second, pitch_bounds[i+1].first, pitch_bounds[i+1].second, yaw_bounds[i+1].first, yaw_bounds[i+1].second, shift_x_bounds[i+1].first, shift_x_bounds[i+1].second, shift_y_bounds[i+1].first, shift_y_bounds[i+1].second, shift_z_bounds[i+1].first, shift_z_bounds[i+1].second, roll_min, roll_max, pitch_min, pitch_max, yaw_min, yaw_max, shift_min_x, shift_max_x, shift_min_y, shift_max_y, shift_min_z, shift_max_z, best_rot_trans, best_ub, point_cloud_model, point_cloud_data);
    }
    return goicp;
}

/* Run ICP on point clouds */
double run_ICP_only(GoICP& goicp, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double shift_min_x, double shift_max_x, double shift_min_y, double shift_max_y, double shift_min_z, double shift_max_z, vector<double>& rot_trans_ub){

    double yaw=(yaw_min+yaw_max)*0.5;
    double pitch=(pitch_min+pitch_max)*0.5;
    double roll=(roll_min+roll_max)*0.5;
    
    func<> r11 = cos(yaw)*cos(roll);
    func<> r12 = cos(yaw)*sin(roll)*sin(pitch) - sin(yaw)*cos(pitch);
    func<> r13 = cos(yaw)*sin(roll)*cos(pitch) + sin(yaw)*sin(pitch);
    func<> r21 = sin(yaw)*cos(roll);
    func<> r22 = sin(yaw)*sin(roll)*sin(pitch) + cos(yaw)*cos(pitch);
    func<> r23 = sin(yaw)*sin(roll)*cos(pitch) - cos(yaw)*sin(pitch);
    func<> r31 = sin(-1*roll);
    func<> r32 = cos(roll)*sin(pitch);
    func<> r33 = cos(roll)*cos(pitch);

    
    goicp.R_init.val[0][0]=cos(yaw)*cos(roll);
    goicp.R_init.val[0][1]=cos(yaw)*sin(roll)*sin(pitch) - sin(yaw)*cos(pitch);
    goicp.R_init.val[0][2]=cos(yaw)*sin(roll)*cos(pitch) + sin(yaw)*sin(pitch);
    goicp.R_init.val[1][0]=sin(yaw)*cos(roll);
    goicp.R_init.val[1][1]=sin(yaw)*sin(roll)*sin(pitch) + cos(yaw)*cos(pitch);
    goicp.R_init.val[1][2]=sin(yaw)*sin(roll)*cos(pitch) - cos(yaw)*sin(pitch);
    goicp.R_init.val[2][0]=sin(-1*roll);
    goicp.R_init.val[2][1]=cos(roll)*sin(pitch);
    goicp.R_init.val[2][2]=cos(roll)*cos(pitch);
    goicp.T_init.val[0][0]=(shift_min_x+shift_max_x)/2.0;
    goicp.T_init.val[1][0]=(shift_min_y+shift_max_y)/2.0;
    goicp.T_init.val[2][0]=(shift_min_z+shift_max_z)/2.0;
    
    auto error=goicp.run_ICP();
    
    auto pitch_sol = atan2(goicp.R_init.val[2][1], goicp.R_init.val[2][2])*180/pi;
    auto roll_sol = atan2(-goicp.R_init.val[2][0], std::sqrt(goicp.R_init.val[2][1]*goicp.R_init.val[2][1]+goicp.R_init.val[2][2]*goicp.R_init.val[2][2]))*180/pi;
    auto yaw_sol = atan2(goicp.R_init.val[1][0],goicp.R_init.val[0][0])*180/pi;
    DebugOff("Roll (degrees) = " << to_string_with_precision(roll,12) << endl);
    DebugOff("Pitch (degrees) = " << to_string_with_precision(pitch,12) << endl);
    DebugOff("Yaw (degrees) = " << to_string_with_precision(yaw,12) << endl);
    auto tx = goicp.T_init.val[0][0];
    auto ty = goicp.T_init.val[1][0];
    auto tz = goicp.T_init.val[2][0];
    DebugOff("tx = " << tx << endl);
    DebugOff("ty = " << ty << endl);
    DebugOff("tz = " << tz << endl);
    DebugOff("ICP error "<<error<<endl);
    rot_trans_ub[0]=goicp.R_init.val[0][0];
    rot_trans_ub[1]=goicp.R_init.val[0][1];
    rot_trans_ub[2]=goicp.R_init.val[0][2];
    rot_trans_ub[3]=goicp.R_init.val[1][0];
    rot_trans_ub[4]=goicp.R_init.val[1][1];
    rot_trans_ub[5]=goicp.R_init.val[1][2];
    rot_trans_ub[6]=goicp.R_init.val[2][0];
    rot_trans_ub[7]=goicp.R_init.val[2][1];
    rot_trans_ub[8]=goicp.R_init.val[2][2];
    rot_trans_ub[9] = goicp.T_init.val[0][0];
    rot_trans_ub[10] = goicp.T_init.val[1][0];
    rot_trans_ub[11] = goicp.T_init.val[2][0];
    
    return error;
}
/* Run ICP on point clouds */
GoICP initialize_ICP_only(const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data){
    using namespace Go_ICP;
    
    int Nm = point_cloud_model.size(), Nd = point_cloud_data.size(), NdDownsampled = 0;
    clock_t  clockBegin, clockEnd;
    string modelFName, dataFName, configFName, outputFname;
    POINT3D * pModel, * pData, * pFullData;
    GoICP goicp;
    set_GoICP_options(goicp);
    // Load model and data point clouds
    pModel = (POINT3D *)malloc(sizeof(POINT3D) * Nm);
    double avg_x = 0, avg_y = 0, avg_z = 0;
    double max_x = numeric_limits<double>::lowest(), max_y = numeric_limits<double>::lowest(), max_z = numeric_limits<double>::lowest();
    double min_x = numeric_limits<double>::max(), min_y = numeric_limits<double>::max(), min_z = numeric_limits<double>::max();
    for(int i = 0; i < Nm; i++)
    {
        pModel[i].x  = point_cloud_model[i][0];
        avg_x += pModel[i].x;
        pModel[i].y  = point_cloud_model[i][1];
        avg_y += pModel[i].y;
        pModel[i].z = point_cloud_model[i][2];
        avg_z += pModel[i].z;
    }
    avg_x /= Nm;avg_y /= Nm;avg_z /= Nm;
    //            centralize(Nm, &pModel, avg_x, avg_y, avg_z);
    avg_x = 0;avg_y = 0;avg_z = 0;
    pData = (POINT3D *)malloc(sizeof(POINT3D) * Nd);
    for(int i = 0; i < Nd; i++)
    {
        pData[i].x  = point_cloud_data[i][0];
        avg_x += pData[i].x;
        pData[i].y  = point_cloud_data[i][1];
        avg_y += pData[i].y;
        pData[i].z = point_cloud_data[i][2];
        avg_z += pData[i].z;
    }
    avg_x /= Nd;avg_y /= Nd;avg_z /= Nd;
    
    bool plot_GoICP = false;
    
    goicp.pModel = pModel;
    goicp.Nm = Nm;
    goicp.pData = pData;
    goicp.Nd = Nd;
    // Build Distance Transform
    cout << "Building Distance Transform..." << flush;
    clockBegin = clock();
    goicp.BuildDT();
    clockEnd = clock();
    cout << (double)(clockEnd - clockBegin)/CLOCKS_PER_SEC << "s (CPU)" << endl;
    
    // Run GO-ICP
    if(NdDownsampled > 0)
    {
        goicp.Nd = NdDownsampled; // Only use first NdDownsampled data points (assumes data points are randomly ordered)
    }
    cout << "Model ID: " << modelFName << " (" << goicp.Nm << "), Data ID: " << dataFName << " (" << goicp.Nd << ")" << endl;
    goicp.R_init=Go_ICP::Matrix::eye(3);
    goicp.T_init=Go_ICP::Matrix::ones(3,1)*0;
    goicp.Initialize();
    
    return goicp;
}

shared_ptr<Model<double>> build_linobj_convex_clean(const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data, const indices& valid_cells, double new_roll_min, double new_roll_max, double new_pitch_min, double new_pitch_max, double new_yaw_min, double new_yaw_max, double new_shift_min_x, double new_shift_max_x, double new_shift_min_y, double new_shift_max_y, double new_shift_min_z, double new_shift_max_z, param<>& dist_cost, double ub, bool& found_all, double lb, param<double>& maxj)
{
    
    int nb_pairs = 0, min_nb_pairs = numeric_limits<int>::max(), max_nb_pairs = 0, av_nb_pairs = 0;
    size_t nm = point_cloud_model.size(), nd = point_cloud_data.size();
    vector<double> zeros = {0,0,0};
    param<> x1("x1"), x2("x2"), y1("y1"), y2("y2"), z1("z1"), z2("z2");
    int m = av_nb_pairs;
    string i_str, j_str;
    double xm_max = numeric_limits<double>::lowest(), ym_max = numeric_limits<double>::lowest(), zm_max = numeric_limits<double>::lowest();\
    double xm_min = numeric_limits<double>::max(), ym_min = numeric_limits<double>::max(), zm_min = numeric_limits<double>::max();
    
    
    for (auto i = 0; i<nd; i++){
        i_str = to_string(i+1);
        x1.add_val(i_str,point_cloud_data.at(i).at(0));
        y1.add_val(i_str,point_cloud_data.at(i).at(1));
        z1.add_val(i_str,point_cloud_data.at(i).at(2));
    }
    for (auto j = 0; j<nm; j++) {
        j_str = to_string(j+1);
        x2.add_val(j_str,point_cloud_model.at(j).at(0));
        if(point_cloud_model.at(j).at(0) > xm_max)
            xm_max = point_cloud_model.at(j).at(0);
        if(point_cloud_model.at(j).at(0) < xm_min)
            xm_min = point_cloud_model.at(j).at(0);
        y2.add_val(j_str,point_cloud_model.at(j).at(1));
        if(point_cloud_model.at(j).at(1) > ym_max)
            ym_max = point_cloud_model.at(j).at(1);
        if(point_cloud_model.at(j).at(1) < ym_min)
            ym_min = point_cloud_model.at(j).at(1);
        z2.add_val(j_str,point_cloud_model.at(j).at(2));
        if(point_cloud_model.at(j).at(2) > zm_max)
            zm_max = point_cloud_model.at(j).at(2);
        if(point_cloud_model.at(j).at(2) < zm_min)
            zm_min = point_cloud_model.at(j).at(2);
    }
    
    
    indices Pairs("Pairs");
    int idx1 = 0;
    int idx2 = 0;
    indices N1("N1"),N2("N2");
    DebugOff("nd = " << nd << endl);
    Debug("nm = " << nm << endl);
    
    N1 = range(1,nd);
    N2 = range(1,nm);
    // auto cells = indices(N1,N2);
    indices cells("cells");
    string name="TU_MIP";
    
    auto Reg=make_shared<Model<>>(name);
    for(auto i=0;i<nd;i++){
        for(auto j=1;j<=nm;j++){
            if(valid_cells.has_key(to_string(i+1)+","+to_string(j)))
                cells.insert(to_string(i+1)+","+to_string(j));
        }
    }
    
    var<int> bin("bin",0,1);
    Reg->add(bin.in(cells));
    
    var<> x_shift("x_shift", new_shift_min_x, new_shift_max_x), y_shift("y_shift", new_shift_min_y, new_shift_max_y), z_shift("z_shift", new_shift_min_z, new_shift_max_z);
    //var<> x_shift("x_shift"), y_shift("y_shift"), z_shift("z_shift");
    
    double shift_mag_max=std::max(pow(new_shift_max_x,2),pow(new_shift_min_x,2))+std::max(pow(new_shift_max_y,2),pow(new_shift_min_y,2))+std::max(pow(new_shift_max_z,2),pow(new_shift_min_z,2));
    double shift_mag_min=0;
    if(new_shift_min_x<=0 && new_shift_max_x>=0){
        shift_mag_min+=0;
    }
    else{
        shift_mag_min+=std::min(pow(new_shift_max_x,2),pow(new_shift_min_x,2));
    }
    if(new_shift_min_y<=0 && new_shift_max_y>=0){
        shift_mag_min+=0;
    }
    else{
        shift_mag_min+=std::min(pow(new_shift_max_y,2),pow(new_shift_min_y,2));
    }
    if(new_shift_min_z<=0 && new_shift_max_z>=0){
        shift_mag_min+=0;
    }
    else{
        shift_mag_min+=std::min(pow(new_shift_max_z,2),pow(new_shift_min_z,2));
    }
    
    DebugOff("Added " << cells.size() << " binary variables" << endl);
    double angle_max = 25.*pi/180.;
    var<> yaw("yaw", new_yaw_min, new_yaw_max), pitch("pitch", new_pitch_min, new_pitch_max), roll("roll", new_roll_min, new_roll_max);
    yaw.in(R(1)); pitch.in(R(1));roll.in(R(1));
    func<> r11 = cos(yaw)*cos(roll);r11.eval_all();
    func<> r12 = cos(yaw)*sin(roll)*sin(pitch) - sin(yaw)*cos(pitch);r12.eval_all();
    func<> r13 = cos(yaw)*sin(roll)*cos(pitch) + sin(yaw)*sin(pitch);r13.eval_all();
    func<> r21 = sin(yaw)*cos(roll);r21.eval_all();
    func<> r22 = sin(yaw)*sin(roll)*sin(pitch) + cos(yaw)*cos(pitch);r22.eval_all();
    func<> r23 = sin(yaw)*sin(roll)*cos(pitch) - cos(yaw)*sin(pitch);r23.eval_all();
    func<> r31 = sin(-1*roll);r31.eval_all();
    func<> r32 = cos(roll)*sin(pitch);r32.eval_all();
    func<> r33 = cos(roll)*cos(pitch);r33.eval_all();
    
    
    var<> theta11("theta11",  std::max(-1.,r11._range->first), std::min(1.,r11._range->second)), theta12("theta12", std::max(-1.,r12._range->first), std::min(1.,r12._range->second)), theta13("theta13", std::max(-1.,r13._range->first), std::min(1.,r13._range->second));
    var<> theta21("theta21", std::max(-1.,r21._range->first), std::min(1.,r21._range->second)), theta22("theta22", std::max(-1.,r22._range->first), std::min(1.,r22._range->second)), theta23("theta23", std::max(-1.,r23._range->first), std::min(1.,r23._range->second));
    var<> theta31("theta31", std::max(-1.,r31._range->first), std::min(1.,r31._range->second)), theta32("theta32", std::max(-1.,r32._range->first), std::min(1.,r32._range->second)), theta33("theta33", std::max(-1.,r33._range->first), std::min(1.,r33._range->second));
    Reg->add(theta11.in(R(1)),theta12.in(R(1)),theta13.in(R(1)));
    Reg->add(theta21.in(R(1)),theta22.in(R(1)),theta23.in(R(1)));
    Reg->add(theta31.in(R(1)),theta32.in(R(1)),theta33.in(R(1)));
    
    
    
    //var<> new_xm("new_xm", xm_min, xm_max), new_ym("new_ym", ym_min, ym_max), new_zm("new_zm", zm_min, zm_max);
    
    
    
//    double tx_max=std::max(pow(new_shift_min_x,2), pow(new_shift_max_x,2));
//    double ty_max=std::max(pow(new_shift_min_y,2), pow(new_shift_max_y,2));
//    double tz_max=std::max(pow(new_shift_min_z,2), pow(new_shift_max_z,2));
//    double tx_min, ty_min, tz_min;
//    double T11_min, T11_max, T12_min,T12_max, T13_min, T13_max;
//    double T21_min, T21_max, T22_min,T22_max, T23_min, T23_max;
//    double T31_min, T31_max, T32_min,T32_max, T33_min, T33_max, Tmin, Tmax;
//    if(new_shift_min_x<=0 && new_shift_max_x>=0)
//        tx_min=0;
//    else
//        tx_min=std::min(pow(new_shift_min_x,2), pow(new_shift_max_x,2));
//    if(new_shift_min_y<=0 && new_shift_max_y>=0)
//        ty_min=0;
//    else
//        ty_min=std::min(pow(new_shift_min_y,2), pow(new_shift_max_y,2));
//    if(new_shift_min_z<=0 && new_shift_max_z>=0)
//        tz_min=0;
//    else
//        tz_min=std::min(pow(new_shift_min_z,2), pow(new_shift_max_z,2));
//
//    Tmin=theta11.get_lb("0");
//    Tmax=theta11.get_ub("0");
//    if(Tmin<=0 && Tmax>=0)
//        T11_min=0;
//    else
//        T11_min=std::min(pow(Tmin,2), pow(Tmax,2));
//    T11_max=std::max(pow(Tmin,2), pow(Tmax,2));
//    Tmin=theta12.get_lb("0");
//    Tmax=theta12.get_ub("0");
//    if(Tmin<=0 && Tmax>=0)
//        T12_min=0;
//    else
//        T12_min=std::min(pow(Tmin,2), pow(Tmax,2));
//    T12_max=std::max(pow(Tmin,2), pow(Tmax,2));
//    Tmin=theta13.get_lb("0");
//    Tmax=theta13.get_ub("0");
//    if(Tmin<=0 && Tmax>=0)
//        T13_min=0;
//    else
//        T13_min=std::min(pow(Tmin,2), pow(Tmax,2));
//    T13_max=std::max(pow(Tmin,2), pow(Tmax,2));
//    Tmin=theta21.get_lb("0");
//    Tmax=theta21.get_ub("0");
//    if(Tmin<=0 && Tmax>=0)
//        T21_min=0;
//    else
//        T21_min=std::min(pow(Tmin,2), pow(Tmax,2));
//    T21_max=std::max(pow(Tmin,2), pow(Tmax,2));
//    Tmin=theta22.get_lb("0");
//    Tmax=theta22.get_ub("0");
//    if(Tmin<=0 && Tmax>=0)
//        T22_min=0;
//    else
//        T22_min=std::min(pow(Tmin,2), pow(Tmax,2));
//    T22_max=std::max(pow(Tmin,2), pow(Tmax,2));
//    Tmin=theta23.get_lb("0");
//    Tmax=theta23.get_ub("0");
//    if(Tmin<=0 && Tmax>=0)
//        T23_min=0;
//    else
//        T23_min=std::min(pow(Tmin,2), pow(Tmax,2));
//    T23_max=std::max(pow(Tmin,2), pow(Tmax,2));
//    Tmin=theta31.get_lb("0");
//    Tmax=theta31.get_ub("0");
//    if(Tmin<=0 && Tmax>=0)
//        T31_min=0;
//    else
//        T31_min=std::min(pow(Tmin,2), pow(Tmax,2));
//    T31_max=std::max(pow(Tmin,2), pow(Tmax,2));
//    Tmin=theta32.get_lb("0");
//    Tmax=theta32.get_ub("0");
//    if(Tmin<=0 && Tmax>=0)
//        T32_min=0;
//    else
//        T32_min=std::min(pow(Tmin,2), pow(Tmax,2));
//    T32_max=std::max(pow(Tmin,2), pow(Tmax,2));
//    Tmin=theta33.get_lb("0");
//    Tmax=theta33.get_ub("0");
//    if(Tmin<=0 && Tmax>=0)
//        T33_min=0;
//    else
//        T33_min=std::min(pow(Tmin,2), pow(Tmax,2));
//    T33_max=std::max(pow(Tmin,2), pow(Tmax,2));
    
    
    var<> tx("tx"), ty("ty"), tz("tz");
    var<> T11("T11"), T12("T12"), T13("T13");
    var<> T21("T21"), T22("T22"), T23("T23");
    var<> T31("T31"), T32("T32"), T33("T33");
    
    // Reg->add(tx.in(R(1)),ty.in(R(1)),tz.in(R(1)));
    Reg->add(x_shift.in(R(1)),y_shift.in(R(1)),z_shift.in(R(1)));
    Reg->add(T11.in(R(1)),T12.in(R(1)),T13.in(R(1)));
    Reg->add(T21.in(R(1)),T22.in(R(1)),T23.in(R(1)));
    Reg->add(T31.in(R(1)),T32.in(R(1)),T33.in(R(1)));
    double ub_root=sqrt(ub);
    //var<> new_xm("new_xm", ub_root*(-1), ub_root), new_ym("new_ym", ub_root*(-1), ub_root), new_zm("new_zm", ub_root*(-1), ub_root);
    var<> new_xm("new_xm"), new_ym("new_ym"), new_zm("new_zm");
    
    
    Constraint<> OneBin("OneBin");
    OneBin = bin.in_matrix(1, 1);
    Reg->add(OneBin.in(N1)==1);
    
    Constraint<> T11sq("T11sq");
    T11sq = pow(theta11,2) - T11;
    //txsq.add_to_callback();
    Reg->add(T11sq.in(range(0,0))==0);
    
    Constraint<> T12sq("T12sq");
    T12sq = pow(theta12,2) - T12;
    //txsq.add_to_callback();
    Reg->add(T12sq.in(range(0,0))==0);
    
    Constraint<> T13sq("T13sq");
    T13sq = pow(theta13,2) - T13;
    //txsq.add_to_callback();
    Reg->add(T13sq.in(range(0,0))==0);
    
    Constraint<> T21sq("T21sq");
    T21sq = pow(theta21,2) - T21;
    //txsq.add_to_callback();
    Reg->add(T21sq.in(range(0,0))==0);
    
    Constraint<> T22sq("T22sq");
    T22sq = pow(theta22,2) - T22;
    //txsq.add_to_callback();
    Reg->add(T22sq.in(range(0,0))==0);
    
    Constraint<> T23sq("T23sq");
    T23sq = pow(theta23,2) - T23;
    //txsq.add_to_callback();
    Reg->add(T23sq.in(range(0,0))==0);
    
    Constraint<> T31sq("T31sq");
    T31sq = pow(theta31,2) - T31;
    //txsq.add_to_callback();
    Reg->add(T31sq.in(range(0,0))==0);
    
    Constraint<> T32sq("T32sq");
    T32sq = pow(theta32,2) - T32;
    //txsq.add_to_callback();
    Reg->add(T32sq.in(range(0,0))==0);
    
    Constraint<> T33sq("T33sq");
    T33sq = pow(theta33,2) - T33;
    //txsq.add_to_callback();
    Reg->add(T33sq.in(range(0,0))==0);
    
    
    Constraint<> diag_1("diag_1");
    diag_1=1-theta11-theta22+theta33;
    Reg->add(diag_1.in(range(0,0))>=0);
    Constraint<> diag_2("diag_2");
    diag_2=1+theta11-theta22-theta33;
    Reg->add(diag_2.in(range(0,0))>=0);
    Constraint<> diag_3("diag_3");
    diag_3=1+theta11+theta22+theta33;
    Reg->add(diag_3.in(range(0,0))>=0);
    Constraint<> diag_4("diag_4");
    diag_4=1-theta11+theta22-theta33;
    Reg->add(diag_4.in(range(0,0))>=0);
    
    Constraint<> soc_12("soc_12");
    soc_12 = pow(theta13+theta31,2)-(1-theta11-theta22+theta33)*(1+theta11-theta22-theta33);
    soc_12.add_to_callback();
    Reg->add(soc_12.in(range(0,0))<=0);
    
    Constraint<> soc_13("soc_13");
    soc_13 = pow(theta12-theta21,2)-(1-theta11-theta22+theta33)*(1+theta11+theta22+theta33);
    soc_13.add_to_callback();
    Reg->add(soc_13.in(range(0,0))<=0);
    
    Constraint<> soc_14("soc_14");
    soc_14 = pow(theta23+theta32,2)-(1-theta11-theta22+theta33)*(1-theta11+theta22-theta33);
    soc_14.add_to_callback();
    Reg->add(soc_14.in(range(0,0))<=0);
    
    Constraint<> soc_23("soc_23");
    soc_23 = pow(theta23-theta32,2)-(1+theta11-theta22-theta33)*(1+theta11+theta22+theta33);
    soc_23.add_to_callback();
    Reg->add(soc_23.in(range(0,0))<=0);
    
    Constraint<> soc_24("soc_24");
    soc_24 = pow(theta12+theta21,2)-(1+theta11-theta22-theta33)*(1-theta11+theta22-theta33);
    soc_24.add_to_callback();
    Reg->add(soc_24.in(range(0,0))<=0);
    
    Constraint<> soc_34("soc_34");
    soc_34 = pow(theta31-theta13,2)-(1+theta11+theta22+theta33)*(1-theta11+theta22-theta33);
    soc_34.add_to_callback();
    Reg->add(soc_34.in(range(0,0))<=0);
    
    Constraint<> det_123("det_123");
    det_123+=(theta13+theta31)*((theta13+theta31)*(1+theta11+theta22+theta33)-(theta23-theta32)*(theta12-theta21));
    det_123-=(1-theta11-theta22+theta33)*((1+theta11-theta22-theta33)*(1+theta11+theta22+theta33)-pow(theta23-theta32,2));
    det_123-=(theta12-theta21)*((theta13+theta31)*(theta23-theta32)-(theta12-theta21)*(1+theta11-theta22-theta33));
    det_123.add_to_callback();
    Reg->add(det_123.in(range(0,0))<=0);
    
    Constraint<> det_124("det_124");
    det_124+=(theta13+theta31)*((theta13+theta31)*(1-theta11+theta22-theta33)-(theta23+theta32)*(theta12+theta21));
    det_124-=(1-theta11-theta22+theta33)*((1+theta11-theta22-theta33)*(1-theta11+theta22-theta33)-pow(theta12+theta21,2));
    det_124-=(theta23+theta32)*((theta13+theta31)*(theta12+theta21)-(theta23+theta32)*(1+theta11-theta22-theta33));
    det_124.add_to_callback();
    Reg->add(det_124.in(range(0,0))<=0);
    
    Constraint<> det_134("det_134");
    det_134+=(theta12-theta21)*((theta12-theta21)*(1-theta11+theta22-theta33)-(theta23+theta32)*(theta31-theta13));
    det_134-=(1-theta11-theta22+theta33)*((1+theta11+theta22+theta33)*(1-theta11+theta22-theta33)-pow(theta31-theta13,2));
    det_134-=(theta23+theta32)*((theta12-theta21)*(theta31-theta13)-(theta23+theta32)*(1+theta11+theta22+theta33));
    det_134.add_to_callback();
    Reg->add(det_134.in(range(0,0))<=0);
    
    Constraint<> det_234("det_234");
    det_234+=(theta23-theta32)*((theta23-theta32)*(1-theta11+theta22-theta33)-(theta12+theta21)*(theta31-theta13));
    det_234-=(1+theta11-theta22-theta33)*((1+theta11+theta22+theta33)*(1-theta11+theta22-theta33)-pow(theta31-theta13,2));
    det_234-=(theta12+theta21)*((theta23-theta32)*(theta31-theta13)-(theta12+theta21)*(1+theta11+theta22+theta33));
    det_234.add_to_callback();
    Reg->add(det_234.in(range(0,0))<=0);
    
    Constraint<> row1("row1");
    row1 = T11+T12+T13;
    // row1.add_to_callback();
    Reg->add(row1.in(range(0,0))==1);
    Constraint<> row2("row2");
    row2 = T21+T22+T23;
    //  row2.add_to_callback();
    Reg->add(row2.in(range(0,0))==1);
    Constraint<> row3("row3");
    row3 = T31+T32+T33;
    // row3.add_to_callback();
    Reg->add(row3.in(range(0,0))==1);
    Constraint<> col1("col1");
    col1 = T11+T21+T31;
    // col1.add_to_callback();
    Reg->add(col1.in(range(0,0))==1);
    Constraint<> col2("col2");
    col2 = T21+T22+T23;
    //col2.add_to_callback();
    Reg->add(col2.in(range(0,0))==1);
    Constraint<> col3("col3");
    col3 = T31+T32+T33;
    //col3.add_to_callback();
    Reg->add(col3.in(range(0,0))==1);
    
    param<> dm("dm");
    for(auto i=0;i<nm;i++){
        auto dmd=pow(point_cloud_model.at(i)[0],2)+pow(point_cloud_model.at(i)[1],2)+pow(point_cloud_model.at(i)[2],2);
        dm.add_val(to_string(i+1), dmd);
    }
 
    double shift_mag_min_root=sqrt(shift_mag_min);
    double shift_mag_max_root=sqrt(shift_mag_max);

    
    indices idsij = indices("idsij");
    idsij.add_empty_row();
    param<> xm_min_i("xm_min_i");
    param<> ym_min_i("ym_min_i");
    param<> zm_min_i("zm_min_i");
    param<> xm_max_i("xm_max_i");
    param<> ym_max_i("ym_max_i");
    param<> zm_max_i("zm_max_i");
    param<> rot_xt_min("rot_xt_min");
    param<> rot_xt_max("rot_xt_max");
    param<> rot_yt_min("rot_yt_min");
    param<> rot_yt_max("rot_yt_max");
    param<> rot_zt_min("rot_zt_min");
    param<> rot_zt_max("rot_zt_max");
    param<> sphere_inner_sq("spher_inner_sq");
    shared_ptr<pair<double,double>> x1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> y1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> z1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> new_x1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> new_y1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> new_z1_bounds = make_shared<pair<double,double>>();
    vector<double> x_lb, x_ub, y_lb, y_ub, z_lb, z_ub;
    found_all=true;
    map<int, int> new_model_pts;
    for(auto i=0;i<nd;i++){
        double siq=0;
        auto d_root=sqrt(pow(point_cloud_data.at(i)[0],2)+pow(point_cloud_data.at(i)[1],2)+pow(point_cloud_data.at(i)[2],2));
        siq=0.0;
        if(shift_mag_min_root<=d_root && shift_mag_max_root>=d_root){
            siq=0.0;
        }
        else{
            siq=std::min(pow((shift_mag_min_root-d_root),2),pow((shift_mag_max_root-d_root),2));
        }
        sphere_inner_sq.add_val(to_string(i+1), siq);
        x1_bounds->first = point_cloud_data.at(i)[0];
        x1_bounds->second = point_cloud_data.at(i)[0];
        y1_bounds->first = point_cloud_data.at(i)[1];
        y1_bounds->second = point_cloud_data.at(i)[1];
        z1_bounds->first = point_cloud_data.at(i)[2];
        z1_bounds->second = point_cloud_data.at(i)[2];
        auto x_range  = get_product_range(x1_bounds, theta11._range);
        auto y_range  = get_product_range(y1_bounds, theta12._range);
        auto z_range  = get_product_range(z1_bounds, theta13._range);
        *new_x1_bounds = {std::max(x_range->first + y_range->first + z_range->first, d_root*(-1)) + new_shift_min_x, std::min(x_range->second + y_range->second + z_range->second, d_root)+ new_shift_max_x};
        if(new_x1_bounds->first>=new_x1_bounds->second-1e-8){
            DebugOn("model infeasible in creation "<<endl);
            found_all=false;
            break;
        }
        x_lb.push_back(new_x1_bounds->first);
        x_ub.push_back(new_x1_bounds->second);
        x_range  = get_product_range(x1_bounds, theta21._range);
        y_range  = get_product_range(y1_bounds, theta22._range);
        z_range  = get_product_range(z1_bounds, theta23._range);
        *new_y1_bounds = {std::max(x_range->first + y_range->first + z_range->first, d_root*(-1)) + new_shift_min_y, std::min(x_range->second + y_range->second + z_range->second, d_root)+ new_shift_max_y};
        if(new_y1_bounds->first>=new_y1_bounds->second-1e-8){
            DebugOn("model infeasible in creation "<<endl);
            found_all=false;
            break;
        }
        y_lb.push_back(new_y1_bounds->first);
        y_ub.push_back(new_y1_bounds->second);
        x_range  = get_product_range(x1_bounds, theta31._range);
        y_range  = get_product_range(y1_bounds, theta32._range);
        z_range  = get_product_range(z1_bounds, theta33._range);
        *new_z1_bounds = {std::max(x_range->first + y_range->first + z_range->first, d_root*(-1)) + new_shift_min_z, std::min(x_range->second + y_range->second + z_range->second, d_root)+ new_shift_max_z};
        if(new_z1_bounds->first>=new_z1_bounds->second-1e-8){
            found_all=false;
            DebugOn("model infeasible in creation "<<endl);
            break;
        }
        z_lb.push_back(new_z1_bounds->first);
        z_ub.push_back(new_z1_bounds->second);
        auto x_min_val=100.0, x_max_val=-999.0;
        auto y_min_val=100.0, y_max_val=-999.0;
        auto z_min_val=100.0, z_max_val=-999.0;
        for(auto j=1;j<=nm;j++){
            if(cells.has_key(to_string(i+1)+","+to_string(j))){
                idsij.add_in_row(i, to_string(i+1)+","+to_string(j));
                if(new_model_pts.find(j)==new_model_pts.end())
                    new_model_pts.insert(std::pair<int,int>(j,1));
                else
                    new_model_pts[j]++;
                auto x=point_cloud_model.at(j-1)[0];
                auto y=point_cloud_model.at(j-1)[1];
                auto z=point_cloud_model.at(j-1)[2];
                if(x>=x_max_val){
                    x_max_val=x;
                }
                if(y>=y_max_val){
                    y_max_val=y;
                }
                if(z>=z_max_val){
                    z_max_val=z;
                }
                if(x<=x_min_val){
                    x_min_val=x;
                }
                if(y<=y_min_val){
                    y_min_val=y;
                }
                if(z<=z_min_val){
                    z_min_val=z;
                }
            }
        }
        rot_xt_min.add_val(to_string(i+1), x_lb[i]);
        rot_xt_max.add_val(to_string(i+1), x_ub[i]);
        rot_yt_min.add_val(to_string(i+1), y_lb[i]);
        rot_yt_max.add_val(to_string(i+1), y_ub[i]);
        rot_zt_min.add_val(to_string(i+1), z_lb[i]);
        rot_zt_max.add_val(to_string(i+1), z_ub[i]);
        xm_min_i.add_val(to_string(i+1), x_min_val);
        xm_max_i.add_val(to_string(i+1), x_max_val);
        ym_min_i.add_val(to_string(i+1), y_min_val);
        ym_max_i.add_val(to_string(i+1), y_max_val);
        zm_min_i.add_val(to_string(i+1), z_min_val);
        zm_max_i.add_val(to_string(i+1), z_max_val);
    }
    if(found_all){
    var<> new_x1("new_x1"),new_y1("new_y1"),new_z1("new_z1");
    var<> deltax("deltax"), deltay("deltay"), deltaz("deltaz");
    //var<> rotxt("rotxt", rot_xt_min, rot_xt_max), rotyt("rotyt", rot_yt_min, rot_yt_max), rotzt("rotzt", rot_zt_min, rot_zt_max);
    var<> rotxt("rotxt"), rotyt("rotyt"), rotzt("rotzt");
    
    Reg->add(rotxt.in(N1));
    Reg->add(rotyt.in(N1));
    Reg->add(rotzt.in(N1));
    Reg->add(new_xm.in(N1));
    Reg->add(new_ym.in(N1));
    Reg->add(new_zm.in(N1));
    Reg->add(deltax.in(N1));
    Reg->add(deltay.in(N1));
    Reg->add(deltaz.in(N1));
    indices ids = indices("in_x");
    ids.add_empty_row();
    
    for(auto i=0;i<nd;i++){
        for(auto j=1;j<=nm;j++){
            if(cells.has_key(to_string(i+1)+","+to_string(j)))
                ids.add_in_row(i, to_string(j));
        }
    }
    
    indices dist_cells=indices("dist_cells");
    vector<int> dist_matched;
    for(auto i=0;i<nd;i++){
        double dmaxij=-999.0;
        int dpos=0;
        auto x1=point_cloud_data.at(i)[0];
        auto y1=point_cloud_data.at(i)[1];
        auto z1=point_cloud_data.at(i)[2];
        if (std::find (dist_matched.begin(), dist_matched.end(), i)==dist_matched.end()){
            for(auto j=i+1;j<nd;j++){
                auto x2=point_cloud_data.at(j)[0];
                auto y2=point_cloud_data.at(j)[1];
                auto z2=point_cloud_data.at(j)[2];
                auto dij=pow(x2-x1,2)+pow(y2-y1,2)+pow(z2-z1,2);
                if(dij>=dmaxij){
                    dmaxij=dij;
                    dpos=j;
                }
            }
            dist_matched.push_back(dpos);
            dist_cells.insert(to_string(i+1)+","+to_string(dpos+1));
        }
    }
//
//    indices angle_cells=indices("angle_cells");
//    vector<int> angle_matched;
//    for(auto i=0;i<nd;i++){
//        double amaxij=99999.0;
//        int apos=0;
//        auto x1p=point_cloud_data.at(i)[0];
//        auto y1p=point_cloud_data.at(i)[1];
//        auto z1p=point_cloud_data.at(i)[2];
//        if (std::find (angle_matched.begin(), angle_matched.end(), i)==angle_matched.end()){
//            for(auto j=i+1;j<nd;j++){
//                auto x2p=point_cloud_data.at(j)[0];
//                auto y2p=point_cloud_data.at(j)[1];
//                auto z2p=point_cloud_data.at(j)[2];
//                auto aij=x1p*x2p+y1p*y2p+z1p*z2p;
//                if(aij<=amaxij){
//                    amaxij=aij;
//                    apos=j;
//                }
//            }
//            angle_matched.push_back(apos);
//            angle_cells.insert(to_string(i+1)+","+to_string(apos+1));
//        }
//    }
    
    auto ids1 = theta11.repeat_id(N1.size());
    Constraint<> Def_newxm("Def_newxm");
    Def_newxm = new_xm-product(x2.in(ids),bin.in_matrix(1, 1));
    Def_newxm += rotxt;
    Reg->add(Def_newxm.in(N1)==0);
    
    Constraint<> Def_newym("Def_newym");
    Def_newym = new_ym-product(y2.in(ids),bin.in_matrix(1, 1));
    Def_newym += rotyt;
    Reg->add(Def_newym.in(N1)==0);
    
    Constraint<> Def_newzm("Def_newzm");
    Def_newzm = new_zm-product(z2.in(ids),bin.in_matrix(1, 1));
    Def_newzm += rotzt;
    Reg->add(Def_newzm.in(N1)==0);
    
    Constraint<> sum_newxm("sum_newxm");
    sum_newxm = sum(new_xm.in(N1));
    //Reg->add(sum_newxm==0);
    
    Constraint<> sum_newym("sum_newym");
    sum_newym = sum(new_ym.in(N1));
    //Reg->add(sum_newym==0);
    
    Constraint<> sum_newzm("sum_newzm");
    sum_newzm = sum(new_zm.in(N1));
    //Reg->add(sum_newzm==0);
    
    

        
//        param<double> xc("xc"), yc("yc"), zc("zc");
//        for(auto i=0;i<nd;i++){
//            for(auto j=1;j<=nm;j++){
//                if(cells.has_key(to_string(i+1)+","+to_string(j)))
//                    xc.add_val(to_string(i+1)+","+ to_string(j), x2.eval(to_string(j)));
//                    yc.add_val(to_string(i+1)+","+ to_string(j), y2.eval(to_string(j)));
//                    zc.add_val(to_string(i+1)+","+ to_string(j), z2.eval(to_string(j)));
//            }
//        }
      //  xc.print();
    auto ids2 = x_shift.repeat_id(N1.size());
    Constraint<> x_rot1("x_rot1");
    x_rot1 += rotxt;
    x_rot1 -= x1.in(N1)*theta11.in(ids1) + y1.in(N1)*theta12.in(ids1) + z1.in(N1)*theta13.in(ids1)+x_shift.in(ids2);
    Reg->add(x_rot1.in(N1)==0);
    
    Constraint<> y_rot1("y_rot1");
    y_rot1 += rotyt;
    y_rot1 -= x1.in(N1)*theta21.in(ids1) + y1.in(N1)*theta22.in(ids1) + z1.in(N1)*theta23.in(ids1)+y_shift.in(ids2);
    Reg->add(y_rot1.in(N1)==0);
    
    Constraint<> z_rot1("z_rot1");
    z_rot1 += rotzt;
    z_rot1 -= x1.in(N1)*theta31.in(ids1) + y1.in(N1)*theta32.in(ids1) + z1.in(N1)*theta33.in(ids1)+z_shift.in(ids2);
    Reg->add(z_rot1.in(N1)==0);
    
//        Constraint<> sum_rotx("sum_rotx");
//        sum_rotx = sum(bin*xc);
//      //  Reg->add(sum_rotx==nd*x_shift);
//
//        Constraint<> sum_roty("sum_roty");
//        sum_roty = sum(bin*yc);
//      //  Reg->add(sum_roty==nd*y_shift);
//
//        Constraint<> sum_rotz("sum_rotz");
//        sum_rotz = sum(bin*zc);
//      //  Reg->add(sum_rotz==nd*z_shift);
        
    Constraint<> Def_deltax("Def_deltax");
    Def_deltax=pow(new_xm, 2)-deltax;
    Reg->add(Def_deltax.in(N1)<=0);
    
    Constraint<> Def_deltay("Def_deltay");
    Def_deltay=pow(new_ym, 2)-deltay;
    Reg->add(Def_deltay.in(N1)<=0);
    
    Constraint<> Def_deltaz("Def_deltaz");
    Def_deltaz=pow(new_zm, 2)-deltaz;
    Reg->add(Def_deltaz.in(N1)<=0);
    
    Constraint<> Def_delta("Def_delta");
    Def_delta=sum(deltax)+sum(deltay)+sum(deltaz);
    //Reg->add(Def_delta>=lb);
        param<double> maxjf("maxjf");
        indices N("N");
        for(auto it=new_model_pts.begin();it!=new_model_pts.end();it++){
            if(maxj.eval(to_string(it->first))<it->second){
                N.insert(to_string(it->first));
                maxjf.add_val(to_string(it->first), maxj.eval(to_string(it->first)));
            }
        }
   
    Constraint<> OneBin2("OneBin2");
    OneBin2 =  bin.in_matrix(0, 1)-maxjf;
    Reg->add(OneBin2.in(N)<=0);
    
    
    
    Constraint<> delta_lower("delta_lower");
    //delta_lower=(x1*x1)+(y1*y1)+(z1*z1)+(product(dm_root.in(idsij), bin.in_matrix(1, 1)))-delta;
    //delta_lower=sum(x1*x1)+sum(y1*y1)+sum(z1*z1)+nd*(x_shift*x_shift+y_shift*y_shift+z_shift*z_shift)+sum(dm_root*bin)-delta;
    //delta_lower=sum(x1*x1)+sum(y1*y1)+sum(z1*z1)+nd*(tx+ty+tz)+sum(dm_cells*bin)-delta;
    //Reg->add(delta_lower<=0);
    
    
    if(dist_cost._indices->_keys->size()!=0){
        Constraint<> delta_cost("delta_cost");
        delta_cost=product(dist_cost.in(idsij), bin.in_matrix(1,1))-deltax-deltay-deltaz;
        Reg->add(delta_cost.in(N1)<=0);
    }
    
    Constraint<> dist_rott("dist_rott");
    dist_rott=rot_xt_min*rot_xt_max+rot_yt_min*rot_yt_max+rot_zt_min*rot_zt_max+sphere_inner_sq-rotxt*(rot_xt_min+rot_xt_max)-rotyt*(rot_yt_min+rot_yt_max)-rotzt*(rot_zt_min+rot_zt_max);
    //Reg->add(dist_rott.in(N1)<=0);
        
        
    
    
    Constraint<> distij_rot("distij_rot");
    distij_rot=(rot_xt_min.from(dist_cells)-rot_xt_max.to(dist_cells))*(rot_xt_max.from(dist_cells)-rot_xt_min.to(dist_cells));
    distij_rot+=(rot_yt_min.from(dist_cells)-rot_yt_max.to(dist_cells))*(rot_yt_max.from(dist_cells)-rot_yt_min.to(dist_cells));
    distij_rot+=(rot_zt_min.from(dist_cells)-rot_zt_max.to(dist_cells))*(rot_zt_max.from(dist_cells)-rot_zt_min.to(dist_cells));
    distij_rot+=pow(x1.from(dist_cells)-x1.to(dist_cells),2)+pow(y1.from(dist_cells)-y1.to(dist_cells),2)+pow(z1.from(dist_cells)-z1.to(dist_cells),2);
    distij_rot-=(rotxt.from(dist_cells)-rotxt.to(dist_cells))*(rot_xt_min.from(dist_cells)-rot_xt_max.to(dist_cells)+rot_xt_max.from(dist_cells)-rot_xt_min.to(dist_cells));
    distij_rot-=(rotyt.from(dist_cells)-rotyt.to(dist_cells))*(rot_yt_min.from(dist_cells)-rot_yt_max.to(dist_cells)+rot_yt_max.from(dist_cells)-rot_yt_min.to(dist_cells));
    distij_rot-=(rotzt.from(dist_cells)-rotzt.to(dist_cells))*(rot_zt_min.from(dist_cells)-rot_zt_max.to(dist_cells)+rot_zt_max.from(dist_cells)-rot_zt_min.to(dist_cells));
     //Reg->add(distij_rot.in(dist_cells)<=0);
    
//    Constraint<> angleij_rot("angleij_rot");
//    angleij_rot=(rot_x_min.from(angle_cells)+rot_x_min.to(angle_cells))*(rot_x_max.from(angle_cells)+rot_x_max.to(angle_cells));
//    angleij_rot+=(rot_y_min.from(angle_cells)+rot_y_min.to(angle_cells))*(rot_y_max.from(angle_cells)+rot_y_max.to(angle_cells));
//    angleij_rot+=(rot_z_min.from(angle_cells)+rot_z_min.to(angle_cells))*(rot_z_max.from(angle_cells)+rot_z_max.to(angle_cells));
//    angleij_rot+=pow((x1.from(angle_cells)+x1.to(angle_cells)),2)+pow((y1.from(angle_cells)+y1.to(angle_cells)),2)+pow((z1.from(angle_cells)+z1.to(angle_cells)),2);
//    angleij_rot-=(rotx.from(angle_cells)+rotx.to(angle_cells))*(rot_x_min.from(angle_cells)+rot_x_min.to(angle_cells)+rot_x_max.from(angle_cells)+rot_x_max.to(angle_cells));
//    angleij_rot-=(roty.from(angle_cells)+roty.to(angle_cells))*(rot_y_min.from(angle_cells)+rot_y_min.to(angle_cells)+rot_y_max.from(angle_cells)+rot_y_max.to(angle_cells));
//    angleij_rot-=(rotz.from(angle_cells)+rotz.to(angle_cells))*(rot_z_min.from(angle_cells)+rot_z_min.to(angle_cells)+rot_z_max.from(angle_cells)+rot_z_max.to(angle_cells));
    //Reg->add(angleij_rot.in(angle_cells)<=0);
    
    
    Constraint<> dist_model("dist_model");
    dist_model=xm_min_i*xm_max_i+ym_min_i*ym_max_i+zm_min_i*zm_max_i+product(dm.in(ids),bin.in_matrix(1, 1))-(new_xm+rotxt)*(xm_min_i+xm_max_i)-(new_ym+rotyt)*(ym_min_i+ym_max_i)-(new_zm+rotzt)*(zm_min_i+zm_max_i);
   // Reg->add(dist_model.in(N1)<=0);
    
//
//        if(!incompatibles.empty()){
//            indices pairs1("pairs1"), pairs2("pairs2");
//            pairs1 = cells;
//            pairs2 = cells;
//            for (const auto &inc_pair : incompatibles) {
//                string key1 = to_string(inc_pair.first.first)+","+to_string(inc_pair.first.second);
//                string key2 = to_string(inc_pair.second.first)+","+to_string(inc_pair.second.second);
//                if(cells.has_key(key1) && cells.has_key(key2)){
//                    pairs1.add_ref(key1);
//                    pairs2.add_ref(key2);
//                }
//            }
//
//            if(pairs1.is_indexed()){
//                DebugOn("Number of incompatible pairs constraints = " << pairs1.size() << endl);
//                Constraint<> incomp_pairs("incomp_pairs");
//                incomp_pairs = bin.in(pairs1) + bin.in(pairs2);
//                Reg->add(incomp_pairs.in(range(1,pairs1.size()))<=1);
//                //incomp_pairs.print();
//                DebugOn("Number of incompatible pairs constraints = " << pairs1.size() << endl);
//
//            }
//            //        incomp_pairs.print();
//        }
    /* Objective function */
    
    func<> obj = 0;
    obj+=sum(deltax)+sum(deltay)+sum(deltaz);
    //obj+=(delta);
    
    
    Reg->min(obj);
    }
    
    //
   // Reg->print();
    return(Reg);
}


#define INTER_MAX_ITER 1000000
unsigned long planeStatPerPair;
//  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -

template< typename Convex >
bool sphericalDisjoint(const Convex & a, const Convex & b) {
    planeStatPerPair=0;
    SphericalPolygon positiveBound, tempPoly;
    int vA, vB;
    Vec3<float> dir(b.vertex(0) - a.vertex(0));
    if( ! a.differenceCoversZeroInDir(b, vA, vB, dir)) return true;
    positiveBound.clear();
    positiveBound.emplace_back(dir);
    positiveBound.clip(b.vertex(vB) - a.vertex(vA), tempPoly); positiveBound.swap(tempPoly);
    if( positiveBound.empty() ) return false;
    do {
        if( ! a.differenceCoversZeroInDir(b,vA, vB, positiveBound.averageDirection()) ) return true;
        if(planeStatPerPair >= INTER_MAX_ITER) {
            //++nbFails;
            DebugOn("Failed in intersection detection "<<endl);
            return false;
        }
        positiveBound.clip(b.vertex(vB) - a.vertex(vA), tempPoly); positiveBound.swap(tempPoly);
        if( positiveBound.empty() ) return false;
    } while( true );
}
template< typename Convex >
bool test_general(vector<vector<double>> vec_a, vector<vector<double>> vec_b){
    vector<Convex> convexes;
    convexes.reserve(2);
    convexes.emplace_back(vec_a);
    convexes.emplace_back(vec_b);
    auto res=sphericalDisjoint<VerticesOnly>(convexes[0], convexes[1]);
    return res;
    
}

void preprocess_poltyope_ve_gjk_in_centroid(const vector<vector<double>>& point_cloud_data, const vector<vector<double>>& point_cloud_model, const indices& old_cells, param<double>& dist_cost_old, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double shift_min_x, double shift_max_x, double shift_min_y, double shift_max_y, double shift_min_z, double shift_max_z, const vector<vector<vector<double>>>& model_voronoi_normals, const vector<vector<double>>& model_face_intercept, const vector<vector<vector<double>>>& model_voronoi_vertices, param<double>& dist_cost, double upper_bound, double lower_bound, double& min_cost_sum, indices& new_cells,  double& new_tx_min, double& new_tx_max, double& new_ty_min, double& new_ty_max, double& new_tz_min, double& new_tz_max, double& prep_time_total, const vector<vector<pair<double, double>>>& model_voronoi_min_max, const vector<vector<vector<int>>>& model_voronoi_vertex_edge, const vector<vector<vector<pair<int,int>>>>& model_voronoi_vertex_edge_planes, const vector<double> max_vert_vert_dist_sq, const param<double>& dii, const param<double>& djj, param<double>& maxj)
{
    bool option_cost_new=true;
    indices valid_cells("valid_cells");
    indices valid_cells_new("valid_cells_new");
    param<double> dist_cost_first("dist_cost_first");
    param<double> dist_cost_max("dist_cost_max");
    param<double> dist_cost_second("dist_cost_second");
    indices valid_cells_empty("valid_cells_empty");
    double time_start = get_wall_time();
    int new_test=0;
    planeStatPerPair = 0;
    shared_ptr<pair<double,double>> new_x1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> new_y1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> new_z1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> x1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> y1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> z1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> rot_x1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> rot_y1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> rot_z1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> tx_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> ty_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> tz_bounds = make_shared<pair<double,double>>();
    tx_bounds->first=shift_min_x;
    tx_bounds->second=shift_max_x;
    ty_bounds->first=shift_min_y;
    ty_bounds->second=shift_max_y;
    tz_bounds->first=shift_min_z;
    tz_bounds->second=shift_max_z;
    var<> yaw("yaw", yaw_min, yaw_max), pitch("pitch", pitch_min, pitch_max), roll("roll", roll_min, roll_max);
    yaw.in(R(1)); pitch.in(R(1));roll.in(R(1));
    func<> r11 = cos(yaw)*cos(roll);r11.eval_all();
    func<> r12 = cos(yaw)*sin(roll)*sin(pitch) - sin(yaw)*cos(pitch);r12.eval_all();
    func<> r13 = cos(yaw)*sin(roll)*cos(pitch) + sin(yaw)*sin(pitch);r13.eval_all();
    func<> r21 = sin(yaw)*cos(roll);r21.eval_all();
    func<> r22 = sin(yaw)*sin(roll)*sin(pitch) + cos(yaw)*cos(pitch);r22.eval_all();
    func<> r23 = sin(yaw)*sin(roll)*cos(pitch) - cos(yaw)*sin(pitch);r23.eval_all();
    func<> r31 = sin(-1*roll);r31.eval_all();
    func<> r32 = cos(roll)*sin(pitch);r32.eval_all();
    func<> r33 = cos(roll)*cos(pitch);r33.eval_all();
    
    
    var<> theta11("theta11",  std::max(-1.,r11._range->first), std::min(1.,r11._range->second)), theta12("theta12", std::max(-1.,r12._range->first), std::min(1.,r12._range->second)), theta13("theta13", std::max(-1.,r13._range->first), std::min(1.,r13._range->second));
    var<> theta21("theta21", std::max(-1.,r21._range->first), std::min(1.,r21._range->second)), theta22("theta22", std::max(-1.,r22._range->first), std::min(1.,r22._range->second)), theta23("theta23", std::max(-1.,r23._range->first), std::min(1.,r23._range->second));
    var<> theta31("theta31", std::max(-1.,r31._range->first), std::min(1.,r31._range->second)), theta32("theta32", std::max(-1.,r32._range->first), std::min(1.,r32._range->second)), theta33("theta33", std::max(-1.,r33._range->first), std::min(1.,r33._range->second));
    vector<var<double>> T1;
    T1.push_back(theta11);
    T1.push_back(theta12);
    T1.push_back(theta13);
    T1.push_back(theta21);
    T1.push_back(theta22);
    T1.push_back(theta23);
    T1.push_back(theta31);
    T1.push_back(theta32);
    T1.push_back(theta33);
    
    
    min_cost_sum=0.0;
    
    
    bool found_all=true;
    new_tx_min=0;
    new_tx_max=0;
    new_ty_min=0;
    new_ty_max=0;
    new_tz_min=0;
    new_tz_max=0;
    
    var<> x_shift("x_shift", shift_min_x, shift_max_x), y_shift("y_shift", shift_min_y, shift_max_y), z_shift("z_shift", shift_min_z, shift_max_z);
    
    map<int, map<double, int, std::greater<double>>> valid_cells_map;
    vector<double> dist_cost_max_min_i;
    map<int, map<double, int, std::greater<double>>> dist_root_novoro_map_j;
    map<int, bool> new_model_pts;
    double dist_alt_cost_sum=0;
    int nd=point_cloud_data.size();
    int nm=point_cloud_model.size();
    vector<double> dist_i;
    vector<double> sphere_inner_sq, sphere_outer_sq;
    vector<double> zeros = {0,0,0};
    vector<shared_ptr<Model<double>>> batch_models;
    int missed=0;
    vector<double> x_lb, x_ub, y_lb, y_ub, z_lb, z_ub;
    vector<vector<vector<double>>> box;
    vector<vector<vector<double>>> vertex_new;
    vector<vector<double>> box_i;
    vector<vector<double>> box_new_i;
    vector<double> coord_i;
    coord_i.resize(3);
    set<int> unique_model_pts;
    min_cost_sum=0.0;
    double ub_root=sqrt(upper_bound);
    vector<vector<int>> vertex_edge(8);
    vector<vector<std::pair<int, int>>> vertex_edge_plane(8);
    vector<vector<double>> planes;
    planes.reserve(6);
    vector<int> vf={0,1,2,3,4,5,6,7,0,1,2,3};
    vector<int> vs={1,2,3,0,5,6,7,4,4,5,6,7};
    vertex_edge[0]={1, 3, 4};
    vertex_edge[1]={0, 2, 5};
    vertex_edge[2]={1, 3, 6};
    vertex_edge[3]={0, 2, 7};
    vertex_edge[4]={5, 7, 0};
    vertex_edge[5]={4, 6, 1};
    vertex_edge[6]={5, 7, 2};
    vertex_edge[7]={4, 6, 3};
    vertex_edge_plane[0]={std::make_pair(1,2), std::make_pair(0,2), std::make_pair(0,1)};
    vertex_edge_plane[1]={std::make_pair(1,2), std::make_pair(2,3), std::make_pair(1,3)};
    vertex_edge_plane[2]={std::make_pair(2,3), std::make_pair(2,4), std::make_pair(3,4)};
    vertex_edge_plane[3]={std::make_pair(0,2), std::make_pair(2,4), std::make_pair(0,4)};
    vertex_edge_plane[4]={std::make_pair(1,5), std::make_pair(0,5), std::make_pair(0,1)};
    vertex_edge_plane[5]={std::make_pair(1,5), std::make_pair(3,5), std::make_pair(1,3)};
    vertex_edge_plane[6]={std::make_pair(3,5), std::make_pair(4,5), std::make_pair(3,4)};
    vertex_edge_plane[7]={std::make_pair(0,5), std::make_pair(4,5), std::make_pair(0,4)};
    double shift_mag_max=std::max(pow(shift_min_x,2),pow(shift_max_x,2))+std::max(pow(shift_min_y,2),pow(shift_max_y,2))+std::max(pow(shift_min_z,2),pow(shift_max_z,2));
    double shift_mag_max_root=sqrt(shift_mag_max);
    double shift_mag_min=0.0;
    if(shift_min_x<=0 && shift_max_x>=0){
        shift_mag_min+=0;
    }
    else{
        shift_mag_min+=std::min(pow(shift_min_x,2),pow(shift_max_x,2));
    }
    if(shift_min_y<=0 && shift_max_y>=0){
        shift_mag_min+=0;
    }
    else{
        shift_mag_min+=std::min(pow(shift_min_y,2),pow(shift_max_y,2));
    }
    if(shift_min_z<=0 && shift_max_z>=0){
        shift_mag_min+=0;
    }
    else{
        shift_mag_min+=std::min(pow(shift_min_z,2),pow(shift_max_z,2));
    }
    double shift_mag_min_root=sqrt(shift_mag_min);
    double siq=0;
    double min_sum_di_sq=0;
    for(auto i=0;i<nd;i++){
        siq=0.0;
        auto d_root=sqrt(pow(point_cloud_data.at(i)[0],2)+pow(point_cloud_data.at(i)[1],2)+pow(point_cloud_data.at(i)[2],2));
        if(shift_mag_min_root<=d_root && shift_mag_max_root>=d_root){
            siq=0.0;
        }
        else{
            siq=std::min(pow((shift_mag_min_root-d_root),2),pow((shift_mag_max_root-d_root),2));
        }
        sphere_inner_sq.push_back(siq);
        sphere_outer_sq.push_back(pow((d_root+shift_mag_max_root),2));
        x1_bounds->first = point_cloud_data.at(i)[0];
        x1_bounds->second = point_cloud_data.at(i)[0];
        y1_bounds->first = point_cloud_data.at(i)[1];
        y1_bounds->second = point_cloud_data.at(i)[1];
        z1_bounds->first = point_cloud_data.at(i)[2];
        z1_bounds->second = point_cloud_data.at(i)[2];
        auto x_range  = get_product_range(x1_bounds, theta11._range);
        auto y_range  = get_product_range(y1_bounds, theta12._range);
        auto z_range  = get_product_range(z1_bounds, theta13._range);
        *new_x1_bounds = {std::max(x_range->first + y_range->first + z_range->first, d_root*(-1)) + shift_min_x, std::min(x_range->second + y_range->second + z_range->second, d_root)+ shift_max_x};
        rot_x1_bounds->first=std::max(x_range->first + y_range->first + z_range->first, d_root*(-1));
        rot_x1_bounds->second=std::min(x_range->second + y_range->second + z_range->second, d_root);
        if(new_x1_bounds->first>=new_x1_bounds->second-1e-8){
            found_all=false;
            break;
        }
        x_lb.push_back(new_x1_bounds->first);
        x_ub.push_back(new_x1_bounds->second);
        x_range  = get_product_range(x1_bounds, theta21._range);
        y_range  = get_product_range(y1_bounds, theta22._range);
        z_range  = get_product_range(z1_bounds, theta23._range);
        *new_y1_bounds = {std::max(x_range->first + y_range->first + z_range->first, d_root*(-1)) + shift_min_y, std::min(x_range->second + y_range->second + z_range->second, d_root)+ shift_max_y};
        rot_y1_bounds->first=std::max(x_range->first + y_range->first + z_range->first, d_root*(-1));
        rot_y1_bounds->second=std::min(x_range->second + y_range->second + z_range->second, d_root);
        if(new_y1_bounds->first>=new_y1_bounds->second-1e-8){
            found_all=false;
            break;
        }
        y_lb.push_back(new_y1_bounds->first);
        y_ub.push_back(new_y1_bounds->second);
        x_range  = get_product_range(x1_bounds, theta31._range);
        y_range  = get_product_range(y1_bounds, theta32._range);
        z_range  = get_product_range(z1_bounds, theta33._range);
        *new_z1_bounds = {std::max(x_range->first + y_range->first + z_range->first, d_root*(-1)) + shift_min_z, std::min(x_range->second + y_range->second + z_range->second, d_root)+ shift_max_z};
        rot_z1_bounds->first=std::max(x_range->first + y_range->first + z_range->first, d_root*(-1));
        rot_z1_bounds->second=std::min(x_range->second + y_range->second + z_range->second, d_root);
        if(new_z1_bounds->first>=new_z1_bounds->second-1e-8){
            found_all=false;
            break;
        }
        z_lb.push_back(new_z1_bounds->first);
        z_ub.push_back(new_z1_bounds->second);
        x_range  = get_product_range(rot_x1_bounds, tx_bounds);
        y_range  = get_product_range(rot_y1_bounds, ty_bounds);
        z_range  = get_product_range(rot_z1_bounds, tz_bounds);
        double dist_cost_alt_i=pow(point_cloud_data[i][0],2)+pow(point_cloud_data[i][1],2)+pow(point_cloud_data[i][2],2)+shift_mag_min;
        double cost_alt=dist_cost_alt_i+2*(x_range->first+y_range->first+z_range->first);
        coord_i[0]=x_lb[i];
        coord_i[1]=y_lb[i];
        coord_i[2]=z_lb[i];
        box_i.push_back(coord_i);
        coord_i[0]=x_ub[i];
        coord_i[1]=y_lb[i];
        coord_i[2]=z_lb[i];
        box_i.push_back(coord_i);
        coord_i[0]=x_ub[i];
        coord_i[1]=y_ub[i];
        coord_i[2]=z_lb[i];
        box_i.push_back(coord_i);
        coord_i[0]=x_lb[i];
        coord_i[1]=y_ub[i];
        coord_i[2]=z_lb[i];
        box_i.push_back(coord_i);
        coord_i[0]=x_lb[i];
        coord_i[1]=y_lb[i];
        coord_i[2]=z_ub[i];
        box_i.push_back(coord_i);
        coord_i[0]=x_ub[i];
        coord_i[1]=y_lb[i];
        coord_i[2]=z_ub[i];
        box_i.push_back(coord_i);
        coord_i[0]=x_ub[i];
        coord_i[1]=y_ub[i];
        coord_i[2]=z_ub[i];
        box_i.push_back(coord_i);
        coord_i[0]=x_lb[i];
        coord_i[1]=y_ub[i];
        coord_i[2]=z_ub[i];
        box_i.push_back(coord_i);
        box.push_back(box_i);
        vector<double> pl(4,0);
        pl[0]=-1;
        pl[3]=x_lb[i];
        planes.push_back(pl);
        pl.clear();
        pl.resize(4,0);
        pl[1]=-1;
        pl[3]=y_lb[i];
        planes.push_back(pl);
        pl.clear();
        pl.resize(4,0);
        pl[2]=-1;
        pl[3]=z_lb[i];
        planes.push_back(pl);
        pl.clear();
        pl.resize(4,0);
        pl[0]=1;
        pl[3]=-x_ub[i];
        planes.push_back(pl);
        pl.clear();
        pl.resize(4,0);
        pl[1]=1;
        pl[3]=-y_ub[i];
        planes.push_back(pl);
        pl.clear();
        pl.resize(4,0);
        pl[2]=1;
        pl[3]=-z_ub[i];
        planes.push_back(pl);
        pl.clear();
        vector<vector<double>> rot_polytope;
        get_extreme_rotation_data(rot_polytope, point_cloud_data.at(i), T1);
        
        double dist_cost_max_min=9999, cost_min=9999;
        for (int j = 0; j<nm; j++) {
            auto vertices=model_voronoi_vertices.at(j);
            if(!old_cells.has_key(to_string(i+1)+","+to_string(j+1))){
                DebugOff("continued");
                continue;
            }
            bool dist_calculated=false;
            bool res=true;
            double dist=0;
            double dist_max=0;
            auto xm= point_cloud_model.at(j)[0];
            auto ym= point_cloud_model.at(j)[1];
            auto zm= point_cloud_model.at(j)[2];
            auto min_max_voro_j=model_voronoi_min_max[j];
            double cost_alt_j=cost_alt;
            double dc_ij=pow(xm,2)+pow(ym,2)+pow(zm,2);
            cost_alt_j+=dc_ij;
            double max_m_box=-999, max_m_ve=-999;
            for(auto k=0;k<box[i].size();k++){
                auto xb=box[i][k][0];
                auto yb=box[i][k][1];
                auto zb=box[i][k][2];
                auto ip=xb*xm+yb*ym+zb*zm;
                if(ip>=max_m_box){
                    max_m_box=ip;
                }
            }
            double dist_novoro=0;
            if(xm>=x_lb[i]-1e-6 && ym>=y_lb[i]-1e-6 && zm>=z_lb[i]-1e-6 && xm<=x_ub[i]-1e-6 && ym<=y_ub[i]-1e-6 && zm<=z_ub[i]-1e-6){
                res=false;
                dist_calculated=true;
            }
            auto resd=min_max_euclidean_distsq_box(box_i,point_cloud_model.at(j));
            dist=resd.first;
            dist_max=resd.second;
            dist=std::max(dist, cost_alt_j-2.0*max_m_box);
            dist_novoro=sqrt(std::max(resd.first, 0.0));
            vector<vector<double>> box_m_t;
            vector<double> m_t(3);
            m_t[0]=xm-shift_max_x;
            m_t[1]=ym-shift_max_y;
            m_t[2]=zm-shift_max_z;
            box_m_t.push_back(m_t);
            m_t[0]=xm-shift_max_x;
            m_t[1]=ym-shift_max_y;
            m_t[2]=zm-shift_min_z;
            box_m_t.push_back(m_t);
            m_t[0]=xm-shift_max_x;
            m_t[1]=ym-shift_min_y;
            m_t[2]=zm-shift_max_z;
            box_m_t.push_back(m_t);
            m_t[0]=xm-shift_max_x;
            m_t[1]=ym-shift_min_y;
            m_t[2]=zm-shift_min_z;
            box_m_t.push_back(m_t);
            m_t[0]=xm-shift_min_x;
            m_t[1]=ym-shift_max_y;
            m_t[2]=zm-shift_max_z;
            box_m_t.push_back(m_t);
            m_t[0]=xm-shift_min_x;
            m_t[1]=ym-shift_max_y;
            m_t[2]=zm-shift_min_z;
            box_m_t.push_back(m_t);
            m_t[0]=xm-shift_min_x;
            m_t[1]=ym-shift_min_y;
            m_t[2]=zm-shift_max_z;
            box_m_t.push_back(m_t);
            m_t[0]=xm-shift_min_x;
            m_t[1]=ym-shift_min_y;
            m_t[2]=zm-shift_min_z;
            box_m_t.push_back(m_t);
            
            auto resd2=distance_polytopes_gjk(rot_polytope,box_m_t);
            
            dist=std::max(resd2, dist);
//            dist_max=std::min(resd2.second, dist_max);
           
            dist_novoro=sqrt(std::max(resd.first, resd2));

         
            if(!dist_calculated){
                res=test_general<VerticesOnly>(box[i],vertices);
                DebugOff("i "<<i <<" j "<<j<<" "<<res<<endl);
            }
            if(!res){
                dist=std::max(dist,dist_cost_old.eval(to_string(i+1)+","+to_string(j+1)));
                double dist_min_v=dist;
                double dist_max_v=dist_max;
                if(dist<=upper_bound*(nd-1)/nd && dist<=dist_cost_max_min){
                    vector<double> vec1;
                    vector<vector<double>> model_halfspaces;
                    for(auto k=0;k<model_voronoi_normals[j].size();k++){
                        vec1.clear();
                        vec1.resize(4);
                        vec1[0]=model_voronoi_normals[j][k][0];
                        vec1[1]=model_voronoi_normals[j][k][1];
                        vec1[2]=model_voronoi_normals[j][k][2];
                        vec1[3]=model_face_intercept[j][k];
                        model_halfspaces.push_back(vec1);
                    }
                    
                    bool status1=true;
                    vector<vector<double>> vec_vertex1;
                    status1=vert_enum(box_i, planes, vertex_edge, vertex_edge_plane,  model_voronoi_vertices[j], model_halfspaces,  model_voronoi_vertex_edge[j], model_voronoi_vertex_edge_planes[j], vec_vertex1, point_cloud_model.at(j));
                    if(status1){
                        vector<vector<double>> mpoint;
                        mpoint.push_back(point_cloud_model.at(j));
                        auto dist1=distance_polytopes_gjk(vec_vertex1, mpoint);
                        double max_m_ve=-999;
                        double dist_mv=-999;
                        for(auto k=0;k<vec_vertex1.size();k++){
                            auto xbv=vec_vertex1[k][0];
                            auto ybv=vec_vertex1[k][1];
                            auto zbv=vec_vertex1[k][2];
                            auto ip=xbv*xm+ybv*ym+zbv*zm;
                            if(ip>=max_m_ve){
                                max_m_ve=ip;
                            }
                            auto ipm=xbv*xbv+ybv*ybv+zbv*zbv;
                            if(ipm>=dist_mv){
                                dist_mv=ipm;
                            }
                            
                        }
                        dist1=std::max(dist1-1e-5, cost_alt_j-2.0*max_m_ve-1e-6);
                        dist_min_v=std::max(dist,dist1);
                        dist_max_v=dist_mv;
                    }
                    if(dist_min_v<=upper_bound*(nd-1)/nd && dist_min_v<=dist_cost_max_min){
                        new_model_pts.insert ( std::pair<int,bool>(j,true) );
                        auto c=std::max(dist_min_v,0.0);
                        if(c<=cost_min){
                            cost_min=c;
                        }
                        valid_cells_map[i].insert(pair<double, int>(dist_novoro, j));
                        dist_root_novoro_map_j[j].insert(pair<double, int>(dist_novoro, i));
                        auto key=to_string(i+1)+","+to_string(j+1);
                        dist_cost_first.add_val(key, c);
                        dist_cost_max.add_val(key, dist_max);
                        valid_cells.insert(key);
                        if(dist_max_v<=dist_cost_max_min){
                            dist_cost_max_min=dist_max_v;
                        }
                    }
                    else{
                        DebugOff("distij "<<dist_min_v<<" upper_bound "<<upper_bound<<endl);
                    }
                }
                else{
                    DebugOff("distij "<<dist<<" upper_bound "<<upper_bound<<endl);
                }
            }
        }
        dist_cost_max_min_i.push_back(dist_cost_max_min);
        if(valid_cells_map[i].size()==0){
            found_all=false;
            DebugOff("i "<<i<<" vmap size "<<valid_cells_map[i].size()<<endl);
            break;
        }
        min_cost_sum+=cost_min;
        if(min_cost_sum>=upper_bound+1e-6){
            DebugOn("min cost_sum "<<min_cost_sum<<" upper bound "<<upper_bound<<endl);
            found_all=false;
            break;
        }
        box_i.clear();
        planes.clear();
    }
 
    if(found_all){
        for (auto j = 0; j<nm; j++) {
            if(new_model_pts.find(j)!=new_model_pts.end()){
                for(auto k=0;k<dist_root_novoro_map_j[j].size();k++){
                    auto it=std::next(dist_root_novoro_map_j[j].begin(),k);
                    auto dkj=it->first;
                    auto max_i=it->second;
                    auto dkj_u=sqrt(dist_cost_max.eval(to_string(max_i+1)+","+ to_string(j+1)));
                    if(dkj>=1e-6){
                        for(auto l=0;l<dist_root_novoro_map_j[j].size();l++){
                            if(l==k){
                                DebugOff("continued");
                                continue;
                            }
                            auto itl=std::next(dist_root_novoro_map_j[j].begin(),l);
                            auto dij=itl->first;
                            auto i=itl->second;
                            auto dik=dii.eval(to_string(i+1)+","+ to_string(max_i+1));
                            if(dkj>=dik){
                                auto dij_root=dkj-dik;
                                if(dij_root>=dij){
                                    dist_cost_first.set_val(to_string(i+1)+","+ to_string(j+1), dij_root*dij_root);
                                    DebugOff("better cost "<<dij<<" "<< dij_root*dij_root<<endl);
                                }
                            }
                            else if(dik>=dkj_u){
                                auto dij_root=dik-dkj_u;
                                if(dij_root>=dij){
                                   // DebugOn("dij_root*dij_root "<<dij_root*dij_root<<"dij "<<dij*dij<<endl);
                                    dist_cost_first.set_val(to_string(i+1)+","+ to_string(j+1), dij_root*dij_root);
                                    DebugOff("better cost "<<dij<<" "<< dij_root*dij_root<<endl);
                                }
                            }
                        }
                    }
                }
            }
        }
        for (auto i = 0; i<nd; i++) {
            for(auto k=0;k<valid_cells_map[i].size();k++){
                auto it=std::next(valid_cells_map[i].begin(),k);
                auto dik=it->first;
                auto max_j=it->second;
                auto dik_u=sqrt(dist_cost_max.eval(to_string(i+1)+","+ to_string(max_j+1)));
                if(dik>=1e-6){
                    for(auto l=0;l<valid_cells_map[i].size();l++){
                        if(l!=k){
                            auto itl=std::next(valid_cells_map[i].begin(),l);
                            auto dij=itl->first;
                            auto j=itl->second;
                            {
                                auto djk=djj.eval(to_string(j+1)+","+ to_string(max_j+1));
                                if(dik>=djk){
                                    auto dij_root=dik-djk;
                                    if(dij_root>=dij){
                                        dist_cost_first.set_val(to_string(i+1)+","+ to_string(j+1), dij_root*dij_root);
                                        DebugOff("better cost "<<dij<<" "<< dij_root*dij_root<<endl);
                                        
                                    }
                                }
                                else if(djk>=dik_u){
                                    auto dij_root=djk-dik_u;
                                    if(dij_root>=dij){
                                        dist_cost_first.set_val(to_string(i+1)+","+ to_string(j+1), dij_root*dij_root);
                                       // DebugOn("djk "<<djk<<"dik_u "<<dik_u<<endl);
                                        DebugOff("better cost "<<dij*dij<<" "<< dij_root*dij_root<<endl);
                                        
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        double min_i=0;
        double min_cost_sum_new=0;
        vector<double> v_min_i;
        for(auto i=0;i<nd;i++)
        {
            min_i=999;
            found_all=false;
            for(auto j=0;j<nm;j++){
                if(valid_cells.has_key(to_string(i+1)+","+to_string(j+1))){
                    auto dij=dist_cost_first.eval(to_string(i+1)+","+ to_string(j+1));
                    if(dij<=upper_bound*(nd-1)/nd && dij<=dist_cost_max_min_i[i]){
                        if(dij<=min_i){
                            min_i=dij;
                        }
                        found_all=true;
                    }
                }
            }
            v_min_i.push_back(min_i);
            if(!found_all){
                DebugOn("found_all_false");
                break;
            }
        }
       
        for(auto i=0;i<nd && found_all;i++){
            min_i=999;
            found_all=false;
            double sum_other_i=0;
            double xm_min=9999.0, ym_min=9999.0, zm_min=9999.0;
            double xm_max=-9999.0, ym_max=-9999.0, zm_max=-9999.0;
            for(auto l=0;l<nd;l++){
                if(l!=i){
                    sum_other_i+=v_min_i[l];
                }
            }
            for(auto j=0;j<nm;j++){
                if(valid_cells.has_key(to_string(i+1)+","+to_string(j+1))){
                    auto xm= point_cloud_model.at(j)[0];
                    auto ym= point_cloud_model.at(j)[1];
                    auto zm= point_cloud_model.at(j)[2];
                    auto dij=dist_cost_first.eval(to_string(i+1)+","+ to_string(j+1));
                    if(dij<=upper_bound*(nd-1)/nd && dij<=dist_cost_max_min_i[i] && dij<=(upper_bound-sum_other_i)){
                        valid_cells_new.insert(to_string(i+1)+","+to_string(j+1));
                        dist_cost_second.add_val(to_string(i+1)+","+to_string(j+1), dij);
                        if(xm<=xm_min){
                            xm_min=xm;
                        }
                        if(xm>=xm_max){
                            xm_max=xm;
                        }
                        if(ym<=ym_min){
                            ym_min=ym;
                        }
                        if(ym>=ym_max){
                            ym_max=ym;
                        }
                        if(zm<=zm_min){
                            zm_min=zm;
                        }
                        if(zm>=zm_max){
                            zm_max=zm;
                        }
                        if(dij<=min_i){
                            min_i=dij;
                        }
                        found_all=true;
                    }
                }
            }
            min_cost_sum_new+=min_i;
            if(min_cost_sum_new>=upper_bound+1e-6){
                DebugOn("min cost sum exceeds ub");
                found_all=false;
            }
            if(!found_all){
                DebugOn("found_all_false");
                break;
            }
            new_tx_min+=xm_min;
            new_tx_max+=xm_max;
            new_ty_min+=ym_min;
            new_ty_max+=ym_max;
            new_tz_min+=zm_min;
            new_tz_max+=zm_max;
        }
        min_cost_sum=min_cost_sum_new;
    }
    if(found_all){
        
        for(auto it=new_model_pts.begin();it!=new_model_pts.end();it++){
            int j=it->first;
             vector<int> numi(nd,0);
             for(auto i=0;i<nd;i++){
                 if(valid_cells_new.has_key(to_string(i+1)+","+to_string(j+1))){
                     numi[i]++;
                     for(auto k=0;k<nd;k++){
                         if(k!=i){
                         if(valid_cells_new.has_key(to_string(k+1)+","+to_string(j+1))){
                             if(pow(dii.eval(to_string(i+1)+","+ to_string(k+1)),2)<=max_vert_vert_dist_sq[j]){
                                 numi[i]++;
                             }
                         }
                         }
                     }
                 }
             }
             maxj.add_val(to_string(j+1), *max_element(numi.begin(), numi.end()));
         }
        double ndd=point_cloud_data.size()*1.0;
        new_tx_min/=ndd;
        new_tx_max/=ndd;
        new_ty_min/=ndd;
        new_ty_max/=ndd;
        new_tz_min/=ndd;
        new_tz_max/=ndd;
        new_tx_min=std::max(new_tx_min-1e-6, shift_min_x);
        new_tx_max=std::min(new_tx_max+1e-6, shift_max_x);
        new_ty_min=std::max(new_ty_min-1e-6, shift_min_y);
        new_ty_max=std::min(new_ty_max+1e-6, shift_max_y);
        new_tz_min=std::max(new_tz_min-1e-6, shift_min_z);
        new_tz_max=std::min(new_tz_max+1e-6, shift_max_z);
        if(new_tx_min>=new_tx_max+1e-10 || new_tx_max<=new_tx_min-1e-10){
            found_all=false;
            DebugOn("new tx lb ub bounds cross "<<new_tx_min<<" "<<new_tx_max<<endl);
        }
        if(new_ty_min>=new_ty_max+1e-10 || new_ty_max<=new_ty_min-1e-10){
            found_all=false;
            DebugOn("new ty lb ub bounds cross "<<new_ty_min<<" "<<new_ty_max<<endl);
        }
        if(new_tz_min>=new_tz_max+1e-10 || new_tz_max<=new_tz_min-1e-10){
            found_all=false;
            DebugOn("new tz lb ub bounds cross "<<new_tz_min<<" "<<new_tz_max<<endl);
        }
        if(new_tx_min>shift_min_x+1e-3){
            DebugOn("new tx lb "<<new_tx_min<<endl);
        }
        if(new_ty_min>shift_min_y+1e-3){
            DebugOn("new ty lb "<<new_ty_min<<endl);
        }
        if(new_tz_min>shift_min_z+1e-3){
            DebugOn("new tz lb "<<new_tz_min<<endl);
        }
        if(new_tx_max<shift_max_x-1e-3){
            DebugOn("new tx ub "<<new_tx_max<<endl);
        }
        if(new_ty_max<shift_max_y-1e-3){
            DebugOn("new ty ub "<<new_ty_max<<endl);
        }
        if(new_tz_max<shift_max_z-1e-3){
            DebugOn("new tz ub "<<new_tz_max<<endl);
        }
    }
    auto time_end = get_wall_time();
    auto prep_time = time_end - time_start;
    prep_time_total+=prep_time;
    DebugOn("Voronoi preprocessing time = " << prep_time << endl);
    if(found_all){
        double cr=((old_cells.size()-valid_cells_new.size())*1.0)/(old_cells.size()*1.0)*100.0;
        DebugOn("cells removed "<<cr<<endl);
    }
    else{
        DebugOn("cells removed "<<100<<endl);
    }
    
    if(found_all){
        min_cost_sum=std::max(min_cost_sum,0.0);
        DebugOn("min_cost_sum "<<min_cost_sum<<endl);
        new_cells=valid_cells_new;
        dist_cost=dist_cost_second;
    }
    else{
        new_cells=valid_cells_empty;
    }
}

bool vert_enum(vector<vector<double>> vertex_set_a, vector<vector<double>> facets_a, vector<vector<int>> vertex_edge_a, vector<vector<pair<int, int>>> vertex_edge_plane_a, vector<vector<double>> vertex_set_b, vector<vector<double>> facets_b,  vector<vector<int>> vertex_edge_b, vector<vector<pair<int, int>>> vertex_edge_plane_b, std::vector<std::vector<double>>& new_vert, const std::vector<double>& model_point){
    bool res=true;
    vector<int> feas_set_a, infeas_set_a,  feas_set_b, infeas_set_b;
    vector<vector<int>> infeas_facet_set_a(vertex_set_a.size()), infeas_facet_set_b(vertex_set_b.size());
    /*k is added to feas_set_a if feasible wrt all facets of B*/
    for(auto k=0;k<vertex_set_a.size();k++){
        auto xk=vertex_set_a[k][0];
        auto yk=vertex_set_a[k][1];
        auto zk=vertex_set_a[k][2];
        bool feas=true;
        vector<int> inf_facet;
        for(auto l=0;l<facets_b.size();l++){
            auto al=facets_b[l][0];
            auto bl=facets_b[l][1];
            auto cl=facets_b[l][2];
            auto dl=facets_b[l][3];
            auto val=al*xk+bl*yk+cl*zk+dl;
            DebugOff("val a"<< al*xk+bl*yk+cl*zk+dl<<endl);
            if(al*xk+bl*yk+cl*zk+dl>=0){
                feas=false;
                inf_facet.push_back(l);
            }
        }
        if(feas){
            feas_set_a.push_back(k);
        }
        else{
            infeas_set_a.push_back(k);
            infeas_facet_set_a[k]=inf_facet;
        }
    }
    /*k is added to feas_set_b if feasible wrt all facets of A*/
    for(auto k=0;k<vertex_set_b.size();k++){
        auto xk=vertex_set_b[k][0];
        auto yk=vertex_set_b[k][1];
        auto zk=vertex_set_b[k][2];
        bool feas=true;
        vector<int> inf_facet;
        for(auto l=0;l<facets_a.size();l++){
            auto al=facets_a[l][0];
            auto bl=facets_a[l][1];
            auto cl=facets_a[l][2];
            auto dl=facets_a[l][3];
            auto valb=al*xk+bl*yk+cl*zk+dl;
            DebugOff("val b"<< al*xk+bl*yk+cl*zk+dl<<endl);
            if(al*xk+bl*yk+cl*zk+dl>=0){
                feas=false;
                inf_facet.push_back(l);
            }
        }
        if(feas){
            feas_set_b.push_back(k);
        }
        else{
            infeas_set_b.push_back(k);
            infeas_facet_set_b[k]=inf_facet;
        }
    }
    new_vert.clear();
//    if(feas_set_a.size()>0 && feas_set_b.size()>0){
//        DebugOn("Both vertices intersect");
//        return false;
//    }
    for(auto i=0;i<feas_set_a.size();i++){
         add_vertex(new_vert, vertex_set_a[feas_set_a[i]], model_point);
    }
    for(auto i=0;i<feas_set_b.size();i++){
     add_vertex(new_vert, vertex_set_b[feas_set_b[i]], model_point);
    }
    
   if(feas_set_a.size()>=0 && feas_set_a.size()<vertex_set_a.size()){
        res=compute_vertices(vertex_set_a, facets_a, facets_b, vertex_edge_a, vertex_edge_plane_a,   feas_set_a, infeas_set_a, infeas_facet_set_a, new_vert, model_point);
    }
    if(feas_set_b.size()>=0 && feas_set_b.size()<vertex_set_b.size()){
        res=compute_vertices(vertex_set_b, facets_b, facets_a, vertex_edge_b, vertex_edge_plane_b,   feas_set_b, infeas_set_b, infeas_facet_set_b, new_vert, model_point);
    }
//    else if(feas_set_a.size()==0 && feas_set_b.size()==0){
//        DebugOn("no vertex intersection found "<<endl);
//        return false;
//    }
    return res;
}
bool compute_vertices(vector<vector<double>> vertex_set_a, vector<vector<double>> facets_a, vector<vector<double>> facets_b, vector<vector<int>> vertex_edge_a, vector<vector<pair<int, int>>> vertex_edge_plane_a,   vector<int> feas_set_a, vector<int> infeas_set_a, vector<vector<int>> infeas_facet_set_a, std::vector<std::vector<double>>& new_vert, const std::vector<double>& model_point){
    bool vertex_found=true, vertex_found_i=true;
    /*Loops over all vertices of A that are  inside B*/
    for(auto i=0;i<feas_set_a.size();i++){
        auto vi=feas_set_a[i];
        auto ve_i=vertex_edge_a[vi];
        auto vep_i=vertex_edge_plane_a[vi];
        auto xa=vertex_set_a[vi][0];
        auto ya=vertex_set_a[vi][1];
        auto za=vertex_set_a[vi][2];
        /*Loops over all vertices that have edge with i*/
        for(auto k=0;k<ve_i.size();k++){
            auto vj=ve_i[k];
            auto xb=vertex_set_a[vj][0];
            auto yb=vertex_set_a[vj][1];
            auto zb=vertex_set_a[vj][2];
            /* finds vertex vj that is infeasible to set B*/
            if(std::find(infeas_set_a.begin(),infeas_set_a.end(), vj)!=infeas_set_a.end()){
                auto edge_pair=vep_i[k];
                vertex_found_i=false;
                auto inf_facet_set_j=infeas_facet_set_a[vj];
//                for(auto m=0;m<infeas_facet_set_a[vi].size();m++){
//                    inf_facet_set_j.push_back(infeas_facet_set_a[vi][m]);
//                }
                /* finds plane of B (inf_facet_set_j[l]) that separates vertices i and vj*/
                for(auto l=0;l<inf_facet_set_j.size();l++){
                    vector<double> plane1, plane2, plane3;
                    plane1=facets_a[edge_pair.first];
                    plane2=facets_a[edge_pair.second];
                    plane3=facets_b[inf_facet_set_j[l]];
                    Eigen::Matrix3f A;
                    Eigen::Vector3f b;
                    A << plane1[0], plane1[1], plane1[2],  plane2[0], plane2[1], plane2[2],  plane3[0], plane3[1], plane3[2];
                    b << plane1[3]*(-1), plane2[3]*(-1), plane3[3]*(-1);
                    Eigen::Vector3f xs = A.fullPivLu().solve(b);
                    double error=(A*xs - b).norm();
                    DebugOff("err "<<error<<endl);
                    if(error<=1e-6){
                        vector<double> solution(3);
                        solution[0]=xs[0];
                        solution[1]=xs[1];
                        solution[2]=xs[2];
                                bool feas1=true, feas2=true;
                        for(auto j=0;j<facets_a.size();j++){
                                    auto al=facets_a[j][0];
                                    auto bl=facets_a[j][1];
                                    auto cl=facets_a[j][2];
                                    auto dl=facets_a[j][3];
                                    if(al*solution[0]+bl*solution[1]+cl*solution[2]+dl>=1e-7){
                                        feas1=false;
                                        if(al*solution[0]+bl*solution[1]+cl*solution[2]+dl<=1e-3){
                                        DebugOff("false "<<al*solution[0]+bl*solution[1]+cl*solution[2]+dl<<endl);
                                        DebugOff("al  "<<al<<" "<<bl<<" "<<cl<<" "<<dl<<endl);
                                        DebugOff("sol  "<<solution[0]<<" "<<solution[1]<<" "<<solution[2] <<endl);
                                        DebugOff("end"<<endl);
                                        }
                                    }
                                }
                                for(auto j=0;j<facets_b.size();j++){
                                    auto al=facets_b[j][0];
                                    auto bl=facets_b[j][1];
                                    auto cl=facets_b[j][2];
                                    auto dl=facets_b[j][3];
                                    if(al*solution[0]+bl*solution[1]+cl*solution[2]+dl>=1e-7){
                                        feas2=false;
                                        if(al*solution[0]+bl*solution[1]+cl*solution[2]+dl<=1e-3){
                                        DebugOff("false "<<al*solution[0]+bl*solution[1]+cl*solution[2]+dl<<endl);
                                        DebugOff("al  "<<al<<" "<<bl<<" "<<cl<<" "<<dl<<endl);
                                        DebugOff("sol  "<<solution[0]<<" "<<solution[1]<<" "<<solution[2] <<endl);
                                        DebugOff("end"<<endl);
                                        }
                                    }
                                }
                                if(feas1 && feas2){
                                    vertex_found_i=true;
                                    add_vertex(new_vert, solution, model_point);
                                }
                            }
                    else{
                        DebugOn("solve failed");
                    }
                }
            }
        }
        if(vertex_found && vertex_found_i){
            vertex_found=true;
        }
        else{
            vertex_found=false;
        }
    }
    vertex_found_i=true;
    for(auto i=0;i<infeas_set_a.size();i++)
    {
        auto vi=infeas_set_a[i];
        auto ve_i=vertex_edge_a[vi];
        auto vep_i=vertex_edge_plane_a[vi];
        auto xa=vertex_set_a[vi][0];
        auto ya=vertex_set_a[vi][1];
        auto za=vertex_set_a[vi][2];
        auto inf_facet_set_i=infeas_facet_set_a[vi];
        /*Loops over all vertices that have edge with i*/
        for(auto k=0;k<ve_i.size();k++){
            auto vj=ve_i[k];
            auto xb=vertex_set_a[vj][0];
            auto yb=vertex_set_a[vj][1];
            auto zb=vertex_set_a[vj][2];
            /* finds vertex vj that is infeasible to set B*/
            if(std::find(infeas_set_a.begin(),infeas_set_a.end(), vj)!=infeas_set_a.end()){
                auto edge_pair=vep_i[k];
                auto inf_facet_set_j=infeas_facet_set_a[vj];
                auto inf_facet_set_u=inf_facet_set_j;
                for(auto m=0;m<infeas_facet_set_a[vi].size();m++){
                    inf_facet_set_u.push_back(infeas_facet_set_a[vi][m]);
                }
                /* finds plane of B (inf_facet_set_j[l]) that separates vertices i and vj*/
                for(auto l=0;l<inf_facet_set_u.size();l++){
                    if(std::find(inf_facet_set_i.begin(), inf_facet_set_i.end(), inf_facet_set_u[l])==inf_facet_set_i.end()||std::find(inf_facet_set_j.begin(), inf_facet_set_j.end(), inf_facet_set_u[l])==inf_facet_set_j.end()){
                        vertex_found_i=false;
                    vector<double> plane1, plane2, plane3;
                    plane1=facets_a[edge_pair.first];
                    plane2=facets_a[edge_pair.second];
                    plane3=facets_b[inf_facet_set_u[l]];
                    Eigen::Matrix3f A;
                    Eigen::Vector3f b;
                    A << plane1[0], plane1[1], plane1[2],  plane2[0], plane2[1], plane2[2],  plane3[0], plane3[1], plane3[2];
                    b << plane1[3]*(-1), plane2[3]*(-1), plane3[3]*(-1);
                    Eigen::Vector3f xs = A.fullPivLu().solve(b);
                    double error=(A*xs - b).norm();
                    if(error<=1e-6){
                        vertex_found_i=true;
                        vector<double> solution(3);
                        solution[0]=xs[0];
                        solution[1]=xs[1];
                        solution[2]=xs[2];
                        bool feas1=true, feas2=true;
                        for(auto j=0;j<facets_a.size();j++){
                                    auto al=facets_a[j][0];
                                    auto bl=facets_a[j][1];
                                    auto cl=facets_a[j][2];
                                    auto dl=facets_a[j][3];
                                    if(al*solution[0]+bl*solution[1]+cl*solution[2]+dl>=1e-7){
                                        feas1=false;
                                        if(al*solution[0]+bl*solution[1]+cl*solution[2]+dl<=1e-3){
                                        DebugOff("false "<<al*solution[0]+bl*solution[1]+cl*solution[2]+dl<<endl);
                                        DebugOff("al  "<<al<<" "<<bl<<" "<<cl<<" "<<dl<<endl);
                                        DebugOff("sol  "<<solution[0]<<" "<<solution[1]<<" "<<solution[2] <<endl);
                                        DebugOff("end"<<endl);
                                    }
                                  }
                                }
                                for(auto j=0;j<facets_b.size();j++){
                                    auto al=facets_b[j][0];
                                    auto bl=facets_b[j][1];
                                    auto cl=facets_b[j][2];
                                    auto dl=facets_b[j][3];
                                    if(al*solution[0]+bl*solution[1]+cl*solution[2]+dl>=1e-7){
                                        feas2=false;
                                        if(al*solution[0]+bl*solution[1]+cl*solution[2]+dl<=1e-3){
                                            DebugOff("false "<<al*solution[0]+bl*solution[1]+cl*solution[2]+dl<<endl);
                                        DebugOff("al  "<<al<<" "<<bl<<" "<<cl<<" "<<dl<<endl);
                                        DebugOff("sol  "<<solution[0]<<" "<<solution[1]<<" "<<solution[2] <<endl);
                                        DebugOff("end"<<endl);
                                        }
                                    }
                                }
                                if(feas1 && feas2){
                                    vertex_found_i=true;
                                    add_vertex(new_vert, solution, model_point);
                                }
                            }
                    else{
                        DebugOn("solve failed");
                    }
                    }
                }
            }
        }
        if(vertex_found && vertex_found_i){
            vertex_found=true;
        }
        else{
            vertex_found=false;
        }
    }
    
    return vertex_found;
}

void add_vertex(std::vector<std::vector<double>>& new_vert, std::vector<double> solution, const std::vector<double>& model_point){
    auto xs=solution[0];
    auto ys=solution[1];
    auto zs=solution[2];
    auto xm=model_point[0];
    auto ym=model_point[1];
    auto zm=model_point[2];
    auto ds=pow(xs-xm,2)+pow(ys-ym,2)+pow(zs-zm,2);
    bool added=false;
    for(auto i=0;i<new_vert.size();i++){
        auto x=new_vert[i][0];
        auto y=new_vert[i][1];
        auto z=new_vert[i][2];
        if(std::abs(x-xs)<=1e-6 && std::abs(y-ys)<=1e-6 && std::abs(z-zs)<=1e-6){
            auto dm=pow(x-xm,2)+pow(y-ym,2)+pow(z-zm,2);
            if(ds<=dm){
                new_vert[i][0]=xs;
                new_vert[i][1]=ys;
                new_vert[i][2]=zs;
            }
            added=true;
            DebugOff("identical vertex found"<<endl);
        }
    }
    if(!added){
        new_vert.push_back(solution);
    }
}


#ifdef USE_GJK
double distance_polytopes_gjk(vector<vector<double>>& vec1, vector<vector<double>>& vec2){
    
    /* Squared distance computed by openGJK.                                 */
    double dd;
    /* Structure of simplex used by openGJK.                                 */
    struct simplex  s;
    /* Number of vertices defining body 1 and body 2, respectively.          */
    int    nvrtx1, nvrtx2;
    
    double **arr1 = (double **)malloc(vec1.size() * sizeof(double *));
    for (int i=0; i<vec1.size(); i++){
        arr1[i] = (double *)malloc(3 * sizeof(double));
    }
    
    /* Read and store vertices' coordinates. */
    for (auto i = 0; i<vec1.size(); i++)
    {
        arr1[i][0]=vec1[i][0];
        arr1[i][1]=vec1[i][1];
        arr1[i][2]=vec1[i][2];
    }
    
    double **arr2 = (double **)malloc(vec2.size() * sizeof(double *));
    for (int i=0; i<vec2.size(); i++){
        arr2[i] = (double *)malloc(3 * sizeof(double));
    }
    
    /* Read and store vertices' coordinates. */
    for (auto i = 0; i<vec2.size(); i++)
    {
        arr2[i][0]=vec2[i][0];
        arr2[i][1]=vec2[i][1];
        arr2[i][2]=vec2[i][2];
    }
    
    nvrtx1=vec1.size();
    nvrtx2=vec2.size();
    
    /* Structures of body 1 and body 2, respectively.                        */
    struct bd       bd1;
    struct bd       bd2;
    
    bd1.coord = arr1;
    bd1.numpoints = nvrtx1;
    
    /* Import coordinates of object 2. */
    
    bd2.coord = arr2;
    bd2.numpoints = nvrtx2;
    
    /* Initialise simplex as empty */
    s.nvrtx = 0;
    
    /* For importing openGJK this is Step 3: invoke the GJK procedure. */
    /* Compute squared distance using GJK algorithm. */
    dd = gjk (bd1, bd2, &s);
    
    for (int i=0; i<bd1.numpoints; i++){
        free(bd1.coord[i]);
    }
    free(bd1.coord);
    for (int i=0; i<bd2.numpoints; i++){
        free(bd2.coord[i]);
    }
    free(bd2.coord);
    
    dd*=dd;
    DebugOff("Distance sqrd "<<dd<<endl);
    return dd;
}

#endif
/*Min and maximum distance squared to given centre ((x-cx)^2+(y-cy)^2+(z-cz)^2) is computed over a polytope Y. Bounding box X, where each vertex corresponds to bounds (not a general polytope). x,y,z in X: xL \le x \le xu, yL \le yU, zL \le z \le zU. Max distance is computed by a general polytope Y given the vertices of the polytope. Minimum distance is calculated over a bounding box X which oversetimates Y
 @coords: vector of extreme points of Y
 @center: point to which distance is computed. center is  avector of coordinates cx,cy,cz*/
pair<double,double> min_max_euclidean_distsq_box(vector<vector<double>> coords, vector<double> center){
    double max_dist=-999.0, min_dist=0.0;
    double xl=99999.0, xu=-999.0, yl=99999.0, yu=-999.0, zl=99999.0, zu=-999.0;
    const double tol=1e-10;
    auto cx=center[0];
    auto cy=center[1];
    auto cz=center[2];
    for(auto i=0;i<coords.size();i++){
        auto x=coords[i][0];
        auto y=coords[i][1];
        auto z=coords[i][2];
        auto d=pow(x-cx,2)+pow(y-cy,2)+pow(z-cz,2);
        if(d>=max_dist){
            max_dist=d;
        }
        if(x<=xl){
            xl=x;
        }
        if(y<=yl){
            yl=y;
        }
        if(z<=zl){
            zl=z;
        }
        if(x>=xu){
            xu=x;
        }
        if(y>=yu){
            yu=y;
        }
        if(z>=zu){
            zu=z;
        }
    }
    if(cx-xl>=-tol && xu-cx>=-tol){
        min_dist+=0.0;
    }
    else if(xu-cx<=-tol){
        min_dist+=pow(xu-cx,2);
    }
    else if(cx-xl<=-tol){
        min_dist+=pow(xl-cx,2);
    }
    
    if(cy-yl>=-tol && yu-cy>=-tol){
        min_dist+=0.0;
    }
    else if(yu-cy<=-tol){
        min_dist+=pow(yu-cy,2);
    }
    else if(cy-yl<=-tol){
        min_dist+=pow(yl-cy,2);
    }
    
    if(cz-zl>=-tol && zu-cz>=-tol){
        min_dist+=0.0;
    }
    else if(zu-cz<=-tol){
        min_dist+=pow(zu-cz,2);
    }
    else if(cz-zl<=-tol){
        min_dist+=pow(zl-cz,2);
    }
    std::pair<double,double> res;
    res={min_dist, max_dist};
    return res;
}
/* Extreme points of the feasible region in which a rotated data point can lie.
 @extreme: Vertices of extreme points of the feasible region: [R[d_pt]]
 @d_pt: Data point coordinate
 @theta_vec: vector of variables defining rotation with roll, pitch, yaw
 SECANT underestimator
 Tangent overestimator
 */
void get_extreme_rotation_data(vector<vector<double>>& extreme, const vector<double>& d_pt, const vector<var<double>>& theta_vec){
    vector<double> d(3);
    
    vector<vector<double>> new_vert_a, new_vert_b;
    vector<int> infeas_set_a, infeas_set_b, infeas_set;
    bool vertex_found_a=false, vertex_found_b=false;
    
    d[0]=d_pt[0];
    d[1]=d_pt[1];
    d[2]=d_pt[2];
    
    double d_mag=pow(d[0],2)+pow(d[1],2)+pow(d[2],2);
    double d_root=sqrt(d_mag);
    
    auto theta11=theta_vec[0];
    auto theta12=theta_vec[1];
    auto theta13=theta_vec[2];
    auto theta21=theta_vec[3];
    auto theta22=theta_vec[4];
    auto theta23=theta_vec[5];
    auto theta31=theta_vec[6];
    auto theta32=theta_vec[7];
    auto theta33=theta_vec[8];
    

    vector<vector<double>> box_big;
    vector<double> coord_i;
    coord_i.resize(3);
    
    shared_ptr<pair<double,double>> x1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> y1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> z1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> rot_x1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> rot_y1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> rot_z1_bounds = make_shared<pair<double,double>>();
    
    x1_bounds->first = d[0];
    x1_bounds->second = d[0];
    y1_bounds->first = d[1];
    y1_bounds->second = d[1];
    z1_bounds->first = d[2];
    z1_bounds->second = d[2];
    auto x_range  = get_product_range(x1_bounds, theta11._range);
    auto y_range  = get_product_range(y1_bounds, theta12._range);
    auto z_range  = get_product_range(z1_bounds, theta13._range);
    
    rot_x1_bounds->first=std::max(x_range->first + y_range->first + z_range->first, d_root*(-1));
    rot_x1_bounds->second=std::min(x_range->second + y_range->second + z_range->second, d_root);
    
    double x_lb=rot_x1_bounds->first;
    double x_ub=rot_x1_bounds->second;
    
    x_range  = get_product_range(x1_bounds, theta21._range);
    y_range  = get_product_range(y1_bounds, theta22._range);
    z_range  = get_product_range(z1_bounds, theta23._range);
    
    rot_y1_bounds->first=std::max(x_range->first + y_range->first + z_range->first, d_root*(-1));
    rot_y1_bounds->second=std::min(x_range->second + y_range->second + z_range->second, d_root);
    
    double y_lb=rot_y1_bounds->first;
    double y_ub=rot_y1_bounds->second;
    
    
    x_range  = get_product_range(x1_bounds, theta31._range);
    y_range  = get_product_range(y1_bounds, theta32._range);
    z_range  = get_product_range(z1_bounds, theta33._range);
    
    rot_z1_bounds->first=std::max(x_range->first + y_range->first + z_range->first, d_root*(-1));
    rot_z1_bounds->second=std::min(x_range->second + y_range->second + z_range->second, d_root);
    
    double z_lb=rot_z1_bounds->first;
    double z_ub=rot_z1_bounds->second;
    
    coord_i[0]=x_lb;
    coord_i[1]=y_lb;
    coord_i[2]=z_lb;
    box_big.push_back(coord_i);
    coord_i[0]=x_ub;
    coord_i[1]=y_lb;
    coord_i[2]=z_lb;
    box_big.push_back(coord_i);
    coord_i[0]=x_ub;
    coord_i[1]=y_ub;
    coord_i[2]=z_lb;
    box_big.push_back(coord_i);
    coord_i[0]=x_lb;
    coord_i[1]=y_ub;
    coord_i[2]=z_lb;
    box_big.push_back(coord_i);
    coord_i[0]=x_lb;
    coord_i[1]=y_lb;
    coord_i[2]=z_ub;
    box_big.push_back(coord_i);
    coord_i[0]=x_ub;
    coord_i[1]=y_lb;
    coord_i[2]=z_ub;
    box_big.push_back(coord_i);
    coord_i[0]=x_ub;
    coord_i[1]=y_ub;
    coord_i[2]=z_ub;
    box_big.push_back(coord_i);
    coord_i[0]=x_lb;
    coord_i[1]=y_ub;
    coord_i[2]=z_ub;
    box_big.push_back(coord_i);
    
    /*Equation of secan which underestimates feasible region*/
    vector<double> secant(4,0.0);
    secant[0]=(x_lb+x_ub)*(-1);
    secant[1]=(y_lb+y_ub)*(-1);
    secant[2]=(z_lb+z_ub)*(-1);
    secant[3]=d_mag+x_lb*x_ub+y_lb*y_ub+z_lb*z_ub;
    
    vertex_found_a=vertices_box_plane_reg(secant, box_big, new_vert_a, infeas_set_a);
    
    /*Find coordinates of the point to define tangent, Tangent is found parallel to secant, if this fails midpoint is used and checked to ensure does not intersect secant in this range*/
   
    double xt, yt, zt;
    bool plane_eq_found=false;
    vector<double> tangent(4);
    double temp_x=sqrt(d_mag/(1+((y_lb+y_ub)*(y_lb+y_ub)+(z_lb+z_ub)*(z_lb+z_ub))/((x_lb+x_ub)*(x_lb+x_ub))));

    double lamda=2.0*temp_x/((x_lb+x_ub)*(-1));
    
    double temp_y=lamda*(y_lb+y_ub)*(-1)/2.0;
    
    double temp_z=lamda*(z_lb+z_ub)*(-1)/2.0;
    
    if(temp_x>=x_lb && temp_x<=x_ub && temp_y>=y_lb && temp_y<=y_ub && temp_z>=z_lb && temp_z<=z_ub){
        xt=temp_x;
        yt=temp_y;
        zt=temp_z;
        plane_eq_found=true;
    }
    else if (temp_x*(-1)>=x_lb && temp_x*(-1)<=x_ub && temp_y*(-1)>=y_lb && temp_y*(-1)<=y_ub && temp_z*(-1)>=z_lb && temp_z*(-1)<=z_ub){
        xt=temp_x*(-1);
        yt=temp_y*(-1);
        zt=temp_z*(-1);
        plane_eq_found=true;
    }
    if(!plane_eq_found){
        xt=(x_lb+x_ub)*0.5;
        yt=(y_lb+y_ub)*0.5;
        auto zttp=sqrt(d_mag-pow(xt,2)-pow(yt,2));
        if( zttp <=z_ub && zttp>=z_lb){
            zt=zttp;
        }
        else{
            zt=zttp*(-1);
        }
        auto dt=-2*(xt*xt+yt*yt+zt*zt);
        plane_eq_found=true;
        if(vertex_found_a){
            for(auto v: new_vert_a){
                auto x=v[0];
                auto y=v[1];
                auto z=v[2];
                if(2*xt*x+2*yt*y+2*zt*z+dt>=1e-9){
                    plane_eq_found=false;
                }
            }
        }
    }
    if(plane_eq_found){
        tangent[0]=2*xt;
        tangent[1]=2*yt;
        tangent[2]=2*zt;
        tangent[3]=-2*(xt*xt+yt*yt+zt*zt);
        DebugOff("Plane Eq found in extreme data "<<endl);
        vertex_found_b=vertices_box_plane_reg(tangent, box_big, new_vert_b, infeas_set_b);
    }
    else{
        DebugOn("Plane Eq not found in extreme data "<<endl);
    }
    if(vertex_found_a){
        for(auto v: new_vert_a){
            extreme.push_back(v);
        }
        for(auto i: infeas_set_a){
            infeas_set.push_back(i);
        }
    }
    else{
        DebugOn("Vertex A no found "<<endl);
    }
    if(vertex_found_b){
        for(auto v: new_vert_b){
            extreme.push_back(v);
        }
        for(auto i: infeas_set_b){
            infeas_set.push_back(i);
        }
    }
    else{
        DebugOn("Vertex B no found "<<plane_eq_found<<endl);
    }
    for(auto i=0;i<box_big.size();i++){
        if(std::find (infeas_set.begin(), infeas_set.end(), i) ==infeas_set.end()){
            extreme.push_back(box_big[i]);
        }
    }

}
/* Extreme points of a bounded box (defined x_lb, x_ub, y_lb. y_ub z_lb, z_ub) when intersected with a plane. Returns true when new vertices found
 @plane_eq: equation of plane intersecting the box
 @big_box: Vector of coordinates of the original bounded box
 @new_verts:New vertices created when plane intersects box
 @infeas_set: Set of indices of old coordinates that are infeasible to plane eq
 */
bool vertices_box_plane_reg(const vector<double>& plane_eq, const vector<vector<double>>& big_box, vector<vector<double>>& new_verts, vector<int>& infeas_set){
    
    const double feas_tol=1e-6;
    
    double a=plane_eq[0];double b=plane_eq[1];double c=plane_eq[2];double d=plane_eq[3];
    
    double x_lb=big_box[0][0];
    double y_lb=big_box[0][1];
    double z_lb=big_box[0][2];
    double x_ub=big_box[6][0];
    double y_ub=big_box[6][1];
    double z_ub=big_box[6][2];
    /*vertex vf[i] and vs[i] are connected by an edge*/
    vector<int> vf={0,1,2,3,4,5,6,7,0,1,2,3};
    vector<int> vs={1,2,3,0,5,6,7,4,4,5,6,7};
    vector<vector<int>> vert_edge(8);
    /*for each vertex i,vert_edge[i] has all the edges that are incident upon the vertex*/
    vert_edge[0]={0, 3, 8};
    vert_edge[1]={0, 1, 9};
    vert_edge[2]={1, 2, 10};
    vert_edge[3]={2, 3, 11};
    vert_edge[4]={4, 7, 8};
    vert_edge[5]={4, 5, 9};
    vert_edge[6]={5, 6, 10};
    vert_edge[7]={6, 7, 11};
    /*for each edge i,plane_x[i] shows whether:
     0: no x plane defines the edge
     -1: plane x=x_lb defines the edge
     1: plane x=x_ub defines the edge
     */
    vector<int> plane_x={0,1,0,-1,0,1,0,-1,-1,1,1,-1};
    vector<int> plane_y={-1,0,1,0,-1,0,1,0,-1,-1,1,1};
    vector<int> plane_z={-1,-1,-1,-1,1,1,1,1,0,0,0,0};
    vector<int> feas_set;
    for(auto k=0;k<8;k++){
        auto xk=big_box[k][0];
        auto yk=big_box[k][1];
        auto zk=big_box[k][2];
        if(a*xk+b*yk+c*zk+d<=feas_tol){
            feas_set.push_back(k);
        }
        else{
            infeas_set.push_back(k);
        }
    }
    bool vertex_found=true;
    for(auto k=0;k<infeas_set.size() && vertex_found;k++){
        bool vertex_found_k=true;
        auto edge_set_k=vert_edge[infeas_set[k]];
        for(auto l=0;l<edge_set_k.size();l++){
            auto el=edge_set_k[l];
            auto elf=vf[el];
            auto els=vs[el];
            /* New vertex must lie on an edge which connects infeasible point inf[k] and some feasible point */
            if(std::find (infeas_set.begin(), infeas_set.end(), elf)==infeas_set.end()||std::find (infeas_set.begin(), infeas_set.end(), els)==infeas_set.end()){
                bool vertex_found_kl=false;
                double xel=100000.0, yel=100000.0, zel=100000.0;
                if(plane_x[el]==-1){
                    xel=x_lb;
                }
                if(plane_x[el]==1){
                    xel=x_ub;
                }
                if(plane_y[el]==-1){
                    yel=y_lb;
                }
                if(plane_y[el]==1){
                    yel=y_ub;
                }
                if(plane_z[el]==-1){
                    zel=z_lb;
                }
                if(plane_z[el]==1){
                    zel=z_ub;
                }
                if(plane_x[el]==0 && abs(a)>1e-10){
                    xel=(-d-c*zel-b*yel)/a;
                    if(xel<=x_ub+1e-9 && xel>=x_lb-1e-9){
                        vertex_found_kl=true;
                    }
                }
                if(plane_y[el]==0 && abs(b)>1e-10){
                    yel=(-d-c*zel-a*xel)/b;
                    if(yel<=y_ub+1e-9 && yel>=y_lb-1e-9){
                        vertex_found_kl=true;
                    }
                }
                if(plane_z[el]==0 && abs(c)>1e-10){
                    zel=(-d-a*xel-b*yel)/c;
                    if(zel<=z_ub+1e-9 && zel>=z_lb-1e-9){
                        vertex_found_kl=true;
                    }
                }
                if(vertex_found_kl){
                    vector<double> v(3);
                    v[0]=xel;
                    v[1]=yel;
                    v[2]=zel;
                    new_verts.push_back(v);
                }
                else{
                    DebugOn(vertex_found_kl<<endl);
                    DebugOn("Failed to find vertex in prep "<<endl);
                }
                if(vertex_found_kl && vertex_found_k){
                    vertex_found_k=true;
                }
                else{
                    vertex_found_k=false;
                    break;
                }
            }
        }
        if(vertex_found_k && vertex_found){
            vertex_found=true;
        }
        else{
            vertex_found=false;
            break;
        }
    }
    return vertex_found;
}

#endif /* Branch_Bound_h */
