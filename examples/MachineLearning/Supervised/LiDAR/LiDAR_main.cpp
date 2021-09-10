//
//  LiDAR_main.cpp
//  Gravity
//
//  Created by Hassan Hijazi on 3 April 2020.
//
//
#include <stdio.h>
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
#include "Lidar_utils.h"
#include "Branch_Bound.h"
#include "IPH.h"
#include "Lower_Bound.h"
#include "Lidar_preprocess.h"
# include "BB.h"
#ifdef USE_EIGEN3
//#include </Users/smitha/Utils/eigen-3.3.9/Eigen/Dense>
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
#ifdef USE_QHULL
#include "libqhullcpp/RboxPoints.h"
#include "libqhullcpp/QhullError.h"
#include "libqhullcpp/QhullQh.h"
#include "libqhullcpp/QhullFacet.h"
#include "libqhullcpp/QhullFacetList.h"
#include "libqhullcpp/QhullFacetSet.h"
#include "libqhullcpp/QhullLinkedList.h"
#include "libqhullcpp/QhullPoint.h"
#include "libqhullcpp/QhullUser.h"
#include "libqhullcpp/QhullVertex.h"
#include "libqhullcpp/QhullVertexSet.h"
#include "libqhullcpp/Qhull.h"

using orgQhull::Qhull;
using orgQhull::QhullError;
using orgQhull::QhullFacet;
using orgQhull::QhullFacetList;
using orgQhull::QhullFacetListIterator;
using orgQhull::QhullFacetSet;
using orgQhull::QhullFacetSetIterator;
using orgQhull::QhullPoint;
using orgQhull::QhullPoints;
using orgQhull::QhullPointsIterator;
using orgQhull::QhullQh;
using orgQhull::QhullUser;
using orgQhull::QhullVertex;
using orgQhull::QhullVertexList;
using orgQhull::QhullVertexListIterator;
using orgQhull::QhullVertexSet;
using orgQhull::QhullVertexSetIterator;
using orgQhull::RboxPoints;


double largest_inscribed_sphere_centre(double x0, double y0, double z0, const vector<double>& point_cloud_data, double& limax, Qhull& qt);
void compute_voronoi(Qhull& qhull, vector<vector<double>>& voronoiVertices, vector<vector<int>>& voronoiRegions);
using namespace orgQhull;
#endif

#ifdef USE_VORO
#include "voro++.hh"
using namespace voro;
#endif

#ifdef USE_CDD
extern "C" {
#include "setoper.h"
#include "cdd.h"
}
#endif

#ifdef USE_GJK
extern "C" {
#include "openGJK.h"
}
#endif

#define DEFAULT_OUTPUT_FNAME "output.txt"
#define DEFAULT_CONFIG_FNAME "config.txt"
#define DEFAULT_MODEL_FNAME "model.txt"
#define DEFAULT_DATA_FNAME "data.txt"


/* Read input files */


/* Set the different options for GoICP  */
void set_GoICP_options(GoICP & goicp);




double get_GoICP_dist(double radius_r, double radius_t, const vector<double>& p, bool L1norm);




/* Plot two point clouds */
void plot(const vector<vector<double>>& point_cloud1, const vector<vector<double>>& point_cloud2, double point_thick = 0.1);
void plot(const vector<vector<double>>& point_cloud1, const vector<vector<double>>& point_cloud2,const vector<vector<double>>& point_cloud3, double point_thick = 0.1);
void plot(const vector<vector<double>>& ext_model, const vector<vector<double>>& ext_data, const vector<vector<double>>& ext_data1,const vector<vector<double>>& ext_data2, double point_thick=0.1);

/* Run Go-ICP on point clouds, return best solution */
tuple<double,double,double,double,double,double,double> run_GoICP(const vector<vector<double>>& point_cloud1, const vector<vector<double>>& point_cloud2);










#ifdef USE_PCL
/* PCL functions */
/* Compute fpfh for each point in ext_model */
pair<pcl::PointCloud<pcl::PointNormal>::Ptr,pcl::PointCloud<pcl::FPFHSignature33>::Ptr> compute_features(const vector<vector<double>>& ext_model);

/* Save the features to a file */
void save_feature_file(const string& filename, const pcl::PointCloud<pcl::PointNormal>::Ptr& augmented_cloud, const pcl::PointCloud<pcl::FPFHSignature33>::Ptr& cloud_features);

#endif

const bool pts_compare(const pair<double,int>& c1, const pair<double,int>& c2) {
    return c1.first > c2.first;
}





int main (int argc, char * argv[])
{
    //    read_laz("/Users/l297598/Downloads/Ta51_powerlines_3__2020_12_18_combined.laz");
    //    return 0;
    string prob_type = "Reg";
    if(argc>1){
        prob_type = argv[1];
    }
#ifdef USE_MPI
    auto err_init = MPI_Init(nullptr,nullptr);
    int worker_id, nb_workers;
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
    auto err_size = MPI_Comm_size(MPI_COMM_WORLD, &nb_workers);
#endif
    
    bool Registration = prob_type=="Reg";/* Solve the Registration problem */
    bool skip_first_line = true; /* First line in Go-ICP input files can be ignored */
    
    /* Boresight Alignment Problem */
    vector<vector<double>> full_point_cloud_model, full_point_cloud_data;
    vector<vector<double>> point_cloud_model, point_cloud_data;
    vector<vector<double>> full_uav_model, full_uav_data;
    vector<vector<double>> uav_model, uav_data;
//    string Model_file = string(prj_dir)+"/data_sets/LiDAR/point_cloud1.txt";
//    string Data_file = string(prj_dir)+"/data_sets/LiDAR/point_cloud1.txt";
//    string red_Model_file = string(prj_dir)+"/data_sets/LiDAR/red_point_cloud1.txt";
//    string red_Data_file = string(prj_dir)+"/data_sets/LiDAR/red_point_cloud1.txt";
//    string algo = "ARMO", global_str = "local";
//
//    if(argc>2){
//        Model_file = argv[2];
//    }
//    if(argc>3){
//        Data_file = argv[3];
//    }
//    if(argc>4){
//        red_Model_file = argv[4];
//        rapidcsv::Document  red_Model_doc(red_Model_file, rapidcsv::LabelParams(0, -1),rapidcsv::SeparatorParams(' '));
//        read_data(red_Model_doc, point_cloud_model, uav_model);
//    }
//    if(argc>5){
//        red_Data_file = argv[5];
//        rapidcsv::Document  red_Data_doc(red_Data_file, rapidcsv::LabelParams(0, -1),rapidcsv::SeparatorParams(' '));
//        read_data(red_Data_doc, point_cloud_data, uav_data);
//
//    }
//    rapidcsv::Document  Model_doc(Model_file, rapidcsv::LabelParams(0, -1),rapidcsv::SeparatorParams(' '));
//    rapidcsv::Document  Data_doc(Data_file, rapidcsv::LabelParams(0, -1),rapidcsv::SeparatorParams(' '));
//
//
    
       //string file_u="/Users/smitha/Desktop/LiDAR_data/Ta51_powerlines_1__2020_12_18_combined.laz";
    string file_u="/Users/smitha/Downloads/Ta51_powerlines_3__2020_12_18_combined.laz";
    //string file_u="/Users/smitha/Desktop/LiDAR_data/DAG4_L_2__2019_06_20_18_combined_RPY_000_frames_701-763_1181-1276.laz";
       auto uav_cloud_u=read_laz(file_u);
       vector<vector<double>> empty_vec, uav_xy;

    double uxmin=numeric_limits<double>::max(), uxmax=numeric_limits<double>::min(), uymin=numeric_limits<double>::max(), uymax=numeric_limits<double>::min(), uzmin=numeric_limits<double>::max(), uzmax=numeric_limits<double>::min();
    vector<vector<double>> uav_coords;
    vector<int> pos(6);
    for(auto i=0;i<uav_cloud_u.size();i++){
        auto x=uav_cloud_u.at(i)[0]-uav_cloud_u.at(0)[0];
        auto y=uav_cloud_u.at(i)[1]-uav_cloud_u.at(0)[1];
        auto z=uav_cloud_u.at(i)[2]-uav_cloud_u.at(0)[2];
        if(x<=uxmin){
            uxmin=x;
            pos[0]=i;
        }
        if(y<=uymin){
            uymin=y;
            pos[1]=i;
        }
        if(z<=uzmin){
            uzmin=z;
            pos[2]=i;
        }
        if(x>=uxmax){
            uxmax=x;
            pos[3]=i;
        }
        if(y>=uymax){
            uymax=y;
            pos[4]=i;
        }
        if(z>=uzmax){
            uzmax=z;
            pos[5]=i;
        }
        uav_xy.push_back({x,y,z});
    }
    for(auto i=0;i<6;i++){
        uav_coords.push_back(uav_cloud_u.at(pos[i]));
    }
    DebugOn("umin max in x "<<uxmin<<" "<<uxmax<<endl);
    DebugOn("umin max in y "<<uymin<<" "<<uymax<<endl);
    DebugOn("umin max in z "<<uzmin<<" "<<uzmax<<endl);
    
    //empty_vec.push_back({uxmin, uymin, 0});
       //plot(uav_xy, empty_vec);
    auto uav_xyz=filter_z_slope(uav_xy);
    
//    for(auto i=0;i<uav_xyz.size();i++){
//        auto x=uav_xyz.at(i)[0];
//        auto y=uav_xyz.at(i)[1];
//        auto z=uav_xyz.at(i)[2]+uav_cloud_u.at(0)[2];
//        uav_xy.push_back({x,y,z});
//    }
    auto slices=extract_slices(uav_xyz);
    for(auto i=0;i<slices.size();i++){
        if(slices[i].size()>10){
        empty_vec.push_back(slices[i][0]);
        
        auto ulist=reg_slope_lines(slices[i]);
            vector<vector<double>> slice_plot, u_plot;
            for(auto j=0;j<slices[i].size();j++){
                auto x=slices[i].at(j)[0]+uav_cloud_u.at(0)[0];
                auto y=slices[i].at(j)[1]+uav_cloud_u.at(0)[1];
                auto z=slices[i].at(j)[2]+uav_cloud_u.at(0)[2];
                slice_plot.push_back({x,y,z});
            }
            plot(slice_plot, uav_coords, 3);
            for(auto j=0;j<ulist.size();j++){
                auto x=ulist.at(j)[0]+uav_cloud_u.at(0)[0];
                auto y=ulist.at(j)[1]+uav_cloud_u.at(0)[1];
                auto z=ulist.at(j)[2]+uav_cloud_u.at(0)[2];
                u_plot.push_back({x,y,z});
            }
        plot(u_plot, slice_plot, uav_coords,3);
        empty_vec.clear();
        }
    }
    auto ulist=reg_slope_lines(slices[1]);
    plot(ulist, uav_xyz, 3);
    
    plot(ulist, empty_vec, 3);
    
  //      plot(ulist, empty_vec);
    //plot(ulist, empty_vec, 10);
       
//       auto resm=read_data(Model_doc, full_point_cloud_model, full_uav_model);
//       auto resd=read_data(Data_doc, full_point_cloud_data, full_uav_data);
//
//       double zl=std::min(resm[2].first, resd[2].first);
//       double zu=std::max(resm[2].second, resd[2].second);
//       double zm_range=(resm[2].second-resm[2].first);
//       double zd_range=(resd[2].second-resd[2].first);
//
//
//       //    read_laz("/Users/smitha/Downloads/Ta51_powerlines_1__2020_12_18_combined.laz");
//       //    ///Users/smitha/Downloads/DAG4_L_2__2019_06_20_18_combined_RPY_000_frames_701-763_1181-1276.laz
//       vector<vector<double>> down_point_cloud_data, down_point_cloud_model, down_uav_data, down_uav_model;
//       vector<vector<double>> down_point_cloud_data1, down_point_cloud_model1, down_uav_data1, down_uav_model1;
//       int xmin=point_cloud_data.at(0)[0];
//       int ymin=point_cloud_data.at(0)[1];
//       int zmin=point_cloud_data.at(0)[2];
//       bool scale=true;
//       if(!scale){
//           xmin=0;
//           ymin=0;
//           zmin=0;
//       }
//
//
//
//        for(auto i=0;i<full_point_cloud_data.size();i++){
//                vector<double> pc1(3);
//                pc1[0]=full_uav_data.at(i)[0]-xmin;
//                pc1[1]=full_uav_data.at(i)[1]-ymin;
//                pc1[2]=full_uav_data.at(i)[2]-zmin;
//                down_uav_data1.push_back(pc1);
//
//                vector<double> pc2(3);
//                pc2[0]=full_point_cloud_data.at(i)[0]-xmin;
//                pc2[1]=full_point_cloud_data.at(i)[1]-ymin;
//                pc2[2]=full_point_cloud_data.at(i)[2]-zmin;
//                down_point_cloud_data1.push_back(pc2);
//
//
//        }
//
//       double yl=numeric_limits<double>::max(), yu=numeric_limits<double>::min();
//
//        for(auto i=0;i<full_point_cloud_model.size();i++){
//                vector<double> pc1(3);
//                pc1[0]=full_point_cloud_model.at(i)[0]-xmin;
//                pc1[1]=full_point_cloud_model.at(i)[1]-ymin;
//                pc1[2]=full_point_cloud_model.at(i)[2]-zmin;
//                down_point_cloud_model1.push_back(pc1);
//
//
//                vector<double> pc2(3);
//                pc2[0]=full_uav_model.at(i)[0]-xmin;
//                pc2[1]=full_uav_model.at(i)[1]-ymin;
//                pc2[2]=full_uav_model.at(i)[2]-zmin;
//                down_uav_model1.push_back(pc2);
//            if(pc2[1]>=yu){
//                yu=pc2[1];
//            }
//            if(pc2[1]<=yl){
//                yl=pc2[1];
//            }
//            }
//       double yrange=(yu-yl);
//       for(auto i=0;i<down_point_cloud_model1.size();i++){
//           double y=down_uav_model1.at(i)[1];
//           double frac=(y-yl)/(yu-yl);
//           if(frac>=0.5 && frac<=0.51){
//           vector<double> pc1(3);
//           pc1[0]=down_point_cloud_model1.at(i)[0];
//           pc1[1]=down_point_cloud_model1.at(i)[1];
//           pc1[2]=down_point_cloud_model1.at(i)[2];
//           down_point_cloud_model.push_back(pc1);
//           vector<double> pc2(3);
//           pc2[0]=down_uav_model1.at(i)[0];
//           pc2[1]=down_uav_model1.at(i)[1];
//           pc2[2]=down_uav_model1.at(i)[2];
//           down_uav_model.push_back(pc2);
//           }
//       }
//
//       for(auto i=0;i<down_point_cloud_data1.size();i++){
//           double y=down_uav_data1.at(i)[1];
//           double frac=(y-yl)/(yu-yl);
//           if(frac>=0.5 && frac<=0.51){
//           vector<double> pc1(3);
//           pc1[0]=down_point_cloud_data1.at(i)[0];
//           pc1[1]=down_point_cloud_data1.at(i)[1];
//           pc1[2]=down_point_cloud_data1.at(i)[2];
//           down_point_cloud_data.push_back(pc1);
//           vector<double> pc2(3);
//           pc2[0]=down_uav_data1.at(i)[0];
//           pc2[1]=down_uav_data1.at(i)[1];
//           pc2[2]=down_uav_data1.at(i)[2];
//           down_uav_data.push_back(pc2);
//           }
//       }
//
//   #ifdef USE_MATPLOT
////       plot(down_point_cloud_data,  down_point_cloud_model);
////
////       plot(down_point_cloud_data1, down_point_cloud_model1);
//   #endif
//
//
//
//    bool run_goICP = false;
//    if(run_goICP){/* Run GoICP inline */
//        run_GoICP(point_cloud_data, point_cloud_model);
//    }
//    param<double> dist_cells("dist_cells");
//
//    int n1=down_point_cloud_data1.size();
//    int n2=down_point_cloud_model1.size();
////    indices N1("N1"),N2("N2");
////    N1 = range(1,n1);
////    N2 = range(1,n2);
////    auto valid_cells_old = indices(N1,N2);
//    indices new_cells("new_cells");
//
//    double roll_min=1.77279-0.1;
//    double roll_max=1.77279+0.1;
//    double pitch_min=-0.03;
//    double pitch_max=-0.01;
//    double yaw_min=-0.4;
//    double yaw_max=-0.3;
//
//
//
//    auto down_point_cloud_model_copy= down_point_cloud_model;
//    auto down_point_cloud_data_copy= down_point_cloud_data;
//    auto down_uav_model_copy= down_uav_model;
//    auto down_uav_data_copy= down_uav_data;
//
//    vector<int> matching(down_point_cloud_model1.size());
//    vector<double> err_per_point(down_point_cloud_model1.size());
//
////
////    auto L2error_init_down = computeL2error(down_point_cloud_model1,down_point_cloud_data1,matching,err_per_point);
////    auto L1error_init_down = computeL1error(down_point_cloud_model1,down_point_cloud_data1,matching,err_per_point);
////
////    DebugOn("L2 error init on down model (20) "<<L2error_init_down<<endl);
////    DebugOn("L1 error init on down model (20) "<<L1error_init_down<<endl);
////
////    vector<int> matching1(point_cloud_data.size());
////    vector<double> err_per_point1(point_cloud_data.size());
////
////    auto L2error_init = computeL2error(point_cloud_model,point_cloud_data,matching1,err_per_point1);
////    auto L1error_init = computeL1error(point_cloud_model,point_cloud_data,matching1,err_per_point1);
////
////
////    DebugOn("L2 error init on down model in paper "<<L2error_init<<endl);
////    DebugOn("L1 error init on down model in paper "<<L1error_init<<endl);
//
//    double final_roll = 0, final_pitch = 0, final_yaw = 0;
//
//
//
//    /*preprocess_lid(down_point_cloud_model1, down_point_cloud_data1, down_uav_model1, down_uav_data1,valid_cells_old,new_cells,  dist_cells, roll_min*pi/180, roll_max*pi/180, pitch_min*pi/180, pitch_max*pi/180, yaw_min*pi/180, yaw_max*pi/180, L2error_init_down);
//
//    auto A_M=Align_model(down_point_cloud_model1, down_point_cloud_data1, down_uav_model1, down_uav_data1,roll_min*pi/180, roll_max*pi/180, pitch_min*pi/180, pitch_max*pi/180, yaw_min*pi/180, yaw_max*pi/180, new_cells, dist_cells);*/
//
////    vector<double> rot(9);
////    bool is_rotation = get_solution(A_M, rot, matching);
////    auto pitch_rad = atan2(rot[7], rot[8]);
////    auto roll_rad = atan2(-rot[6], std::sqrt(rot[7]*rot[7]+rot[8]*rot[8]));
////    auto yaw_rad = atan2(rot[3],rot[0]);
//
//
//
////    DebugOn("Angle in radians roll "<<roll_rad<<endl);
////    DebugOn("Angle in radians pitch "<<pitch_rad<<endl);
////    DebugOn("Angle in radians yaw "<<yaw_rad<<endl);
////
////    Angle in deg roll 4.65224
////    Angle in deg pitch -5.01149
////    Angle in deg yaw 5.01657
//
//    /*auto res=run_IPH(down_point_cloud_model1, down_point_cloud_data1, down_uav_model1, down_uav_data1);
//
//    double roll_deg=get<0>(res);
//    double pitch_deg=get<1>(res);
//    double yaw_deg=get<2>(res);*/
////
////    double roll_deg=roll_rad*180/pi;
////    double pitch_deg=pitch_rad*180/pi;
////    double yaw_deg=yaw_rad*180/pi;
//
//    double roll_deg=-0.25;
//    double pitch_deg=0.5;
//    double yaw_deg=-0.7;
//
//    DebugOn("Angle in deg roll "<<roll_deg<<endl);
//    DebugOn("Angle in deg pitch "<<pitch_deg<<endl);
//    DebugOn("Angle in deg yaw "<<yaw_deg<<endl);
//
//
//
////    apply_rotation(roll_deg, pitch_deg, yaw_deg, down_point_cloud_model, down_point_cloud_data, down_uav_model, down_uav_data);
////
////    apply_rotation(roll_deg, pitch_deg, yaw_deg, down_point_cloud_model1, down_point_cloud_data1, down_uav_model1, down_uav_data1);
//
////    auto L2error = computeL2error(down_point_cloud_model1,down_point_cloud_data1,matching,err_per_point);
////    auto L1error = computeL1error(down_point_cloud_model1,down_point_cloud_data1,matching,err_per_point);
//
////    DebugOn("L2 error final on down model (20) "<<L2error<<endl);
////    DebugOn("L1 error final on down model (20) "<<L1error<<endl);
//////
//   // apply_rotation(roll_deg, pitch_deg, yaw_deg, point_cloud_model, point_cloud_data, uav_model, uav_data);
//
////    L2error = computeL2error(point_cloud_model,point_cloud_data,matching1,err_per_point1);
////    L1error = computeL1error(point_cloud_model,point_cloud_data,matching1,err_per_point1);
//
////    DebugOn("L2 error final on down model in paper "<<L2error<<endl);
////    DebugOn("L1 error final on down model in paper "<<L1error<<endl);
//#ifdef USE_MATPLOT
////    plot(down_point_cloud_data, down_point_cloud_model);
////    plot(down_point_cloud_data1,  down_point_cloud_model1);
////    plot(down_point_cloud_data1, down_point_cloud_model1, down_uav_data1, down_uav_model1);
//    //  plot(point_cloud_model, point_cloud_data);
//#endif
//
//    double total_time =0, time_start = 0, time_end = 0;
//
//
//
//    apply_rotation(roll_deg, pitch_deg, yaw_deg, full_point_cloud_model, full_point_cloud_data, full_uav_model, full_uav_data);
//
//    DebugOff("Initial L2 error = " << L2error_init << endl);
//    DebugOff("Initial L1 error = " << L1error_init << endl);
//    time_start = get_wall_time();
//
//    final_roll=roll_deg;
//    final_pitch=pitch_deg;
//    final_yaw=yaw_deg;
//
//    /*DebugOn("n1 = " << down_point_cloud_model.size() << endl);
//     DebugOn("n2 = " << down_point_cloud_data.size() << endl);
//     DebugOn("Initial L2 error = " << L2error_init << endl);
//     DebugOn("Final L2 error = " << L2error << endl);
//     DebugOn("Relative L2 gap = " << 100*(L2error_init-L2error)/L2error_init << "%" << endl);
//     DebugOn("Initial L1 error = " << L1error_init << endl);
//     DebugOn("Final L1 error = " << L1error << endl);
//     DebugOn("Relative L1 gap = " << 100*(L1error_init-L1error)/L1error_init << "%" << endl);
//     time_end = get_wall_time();
//     DebugOn("Total wall clock time = " << time_end - time_start << endl);
//     double shifted_x, shifted_y, shifted_z;
//     auto tot_pts = full_point_cloud_model.size()+full_point_cloud_data.size();
//     apply_rotation(final_roll, final_pitch, final_yaw, full_point_cloud_data, full_point_cloud_model, full_uav_data, full_uav_model);*/
//
//    bool save_file = true;
//    if(save_file){
//        auto name = Model_file.substr(0,Model_file.find('.'));
//        auto fname = name+"_ARMO_RPY_"+to_string(final_roll)+"_"+to_string(final_pitch)+"_"+to_string(final_yaw)+".laz";
//        save_laz(fname,full_point_cloud_data, full_point_cloud_model);
//        fname = name+"_ARMO_RPY_"+to_string(final_roll)+"_"+to_string(final_pitch)+"_"+to_string(final_yaw)+"_sub.laz";
//        save_laz(fname,down_point_cloud_data, down_point_cloud_model);
//    }
    return 0;
}





void set_GoICP_options(GoICP& goicp){
    goicp.MSEThresh = 1e-3;
    goicp.initNodeRot.a = -2;
    goicp.initNodeRot.b = -2;
    goicp.initNodeRot.c = -2;
    goicp.initNodeRot.w = 2;
    goicp.initNodeTrans.x = -0.5;
    goicp.initNodeTrans.y = -0.5;
    goicp.initNodeTrans.z = -0.5;
    goicp.initNodeTrans.w = 2;
    goicp.trimFraction = 0;
    //    goicp.optError = 11;
    // If < 0.1% trimming specified, do no trimming
    if(goicp.trimFraction < 0.001)
    {
        goicp.doTrim = false;
    }
    goicp.dt.SIZE = 300;
    goicp.dt.expandFactor = 2.0;
}

tuple<double,double,double,double,double,double,double> run_GoICP(const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data){
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
