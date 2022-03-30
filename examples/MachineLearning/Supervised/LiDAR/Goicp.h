//
//  Goicp.h
//  Gravity
//
//  Created by Smitha on 3/25/22.
//

#ifndef Goicp_h
#define Goicp_h
#include <gravity/jly_goicp.h>
#include <gravity/ConfigMap.hpp>
#include "Lidar_utils.h"
#include <DataSet.h>
using namespace Go_ICP;
void set_GoICP_options(GoICP& goicp, double mse=1e-3, int dt=300);
double run_ICP_only(GoICP& goicp, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double shift_min_x, double shift_max_x, double shift_min_y, double shift_max_y, double shift_min_z, double shift_max_z, vector<double>& rot_trans_ub);
void compute_upper_boundICP(GoICP& goicp, double roll_mini, double roll_maxi, double pitch_mini, double pitch_maxi, double yaw_mini, double yaw_maxi, double shift_min_xi, double shift_max_xi, double shift_min_yi, double shift_max_yi, double shift_min_zi, double shift_max_zi, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double shift_min_x, double shift_max_x, double shift_min_y, double shift_max_y, double shift_min_z, double shift_max_z, vector<double>& best_rot_trans, double& best_ub, const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data);
GoICP Initialize_BB(vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, const vector<pair<double, double>>& min_max_model, vector<pair<double, double>>& min_max_d, double shift_min_x, double shift_max_x, double shift_min_y, double shift_max_y, double shift_min_z, double shift_max_z, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double& best_ub, vector<double>& best_rot_trans);
GoICP initialize_ICP_only(const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data);
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
/* Run ICP on point clouds */
double run_ICP_only(GoICP& goicp, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double shift_min_x, double shift_max_x, double shift_min_y, double shift_max_y, double shift_min_z, double shift_max_z, vector<double>& rot_trans_ub){
 
    double yaw=(yaw_min+yaw_max)*0.5;
    double pitch=(pitch_min+pitch_max)*0.5;
    double roll=(roll_min+roll_max)*0.5;
    
    
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

tuple<double,double,double,double,double,double,double> run_GoICP(const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data, double mse, int dt){
    using namespace Go_ICP;
    
    int Nm = point_cloud_model.size(), Nd = point_cloud_data.size(), NdDownsampled = 0;
    clock_t  clockBegin, clockEnd;
    string modelFName, dataFName, configFName, outputFname;
    POINT3D * pModel, * pData, * pFullData;
    GoICP goicp;
    set_GoICP_options(goicp, mse, dt);
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
    cout << "Registering..." << endl;
    clockBegin = clock();
    goicp.Register();
    clockEnd = clock();
    double time = (double)(clockEnd - clockBegin)/CLOCKS_PER_SEC;
    cout << "Optimal Rotation Matrix:" << endl;
    cout << goicp.optR << endl;
    cout << "Optimal Translation Vector:" << endl;
    cout << goicp.optT << endl;
    cout << "Finished in " << time << endl;
    auto pitch = atan2(goicp.optR.val[2][1], goicp.optR.val[2][2])*180/pi;
    auto roll = atan2(-goicp.optR.val[2][0], std::sqrt(goicp.optR.val[2][1]*goicp.optR.val[2][1]+goicp.optR.val[2][2]*goicp.optR.val[2][2]))*180/pi;
    auto yaw = atan2(goicp.optR.val[1][0],goicp.optR.val[0][0])*180/pi;
    DebugOn("Roll (degrees) = " << to_string_with_precision(roll,12) << endl);
    DebugOn("Pitch (degrees) = " << to_string_with_precision(pitch,12) << endl);
    DebugOn("Yaw (degrees) = " << to_string_with_precision(yaw,12) << endl);
    auto tx = goicp.optT.val[0][0];
    auto ty = goicp.optT.val[1][0];
    auto tz = goicp.optT.val[2][0];
    DebugOn("tx = " << tx << endl);
    DebugOn("ty = " << ty << endl);
    DebugOn("tz = " << tz << endl);
    DebugOn("err "<<goicp.optError<<endl);
    return {roll,pitch,yaw,tx,ty,tz,goicp.optError};
}
void compute_upper_boundICP(GoICP& goicp, double roll_mini, double roll_maxi, double pitch_mini, double pitch_maxi, double yaw_mini, double yaw_maxi, double shift_min_xi, double shift_max_xi, double shift_min_yi, double shift_max_yi, double shift_min_zi, double shift_max_zi, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double shift_min_x, double shift_max_x, double shift_min_y, double shift_max_y, double shift_min_z, double shift_max_z, vector<double>& best_rot_trans, double& best_ub, const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data){
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
//                    best_ub=L2erri;
//                    best_rot_trans=rot_trans_ub;
                    DebugOn("best ub new "<<L2erri<<" ub "<<ub<<endl);
                    DebugOn("roll rad"<<roll_rad2<<endl);
                    DebugOn("pitch rad"<<pitch_rad2<<endl);
                    DebugOn("yaw rad"<<yaw_rad2<<endl);
                    DebugOn("roll deg "<<roll_rad2*180/pi<<endl);
                    DebugOn("pitch deg "<<pitch_rad2*180/pi<<endl);
                    DebugOn("yaw deg "<<yaw_rad2*180/pi<<endl);
                    DebugOn("tx "<<rot_trans_ub[9]<<endl);
                    DebugOn("ty "<<rot_trans_ub[10]<<endl);
                    DebugOn("tz "<<rot_trans_ub[11]<<endl);
                    best_rot_trans[0]=L2erri;
                    best_rot_trans[1]=roll_rad2;
                    best_rot_trans[2]=pitch_rad2;
                    best_rot_trans[3]=yaw_rad2;
                    best_rot_trans[4]=rot_trans_ub[9];
                    best_rot_trans[5]=rot_trans_ub[10];
                    best_rot_trans[6]=rot_trans_ub[11];
                   
                }
            }
        }
        
        //bool is_rotation  = get_solution(models[pos], sol_gur, new_matching);
    }
}
void set_GoICP_options(GoICP& goicp, double mse, int dt){
    goicp.MSEThresh = mse;
    goicp.initNodeRot.a = -pi/2;
    goicp.initNodeRot.b = -pi/2;
    goicp.initNodeRot.c = -pi/2;
    goicp.initNodeRot.w = pi;
    goicp.initNodeTrans.x = -0.5;
    goicp.initNodeTrans.y = -0.5;
    goicp.initNodeTrans.z = -0.5;
    goicp.initNodeTrans.w = 1.0;
    goicp.trimFraction = 0;
    //    goicp.optError = 11;
    // If < 0.1% trimming specified, do no trimming
    if(goicp.trimFraction < 0.001)
    {
        goicp.doTrim = false;
    }
    goicp.dt.SIZE = dt;
    goicp.dt.expandFactor = 2.0;
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
#endif /* Goicp_h */
