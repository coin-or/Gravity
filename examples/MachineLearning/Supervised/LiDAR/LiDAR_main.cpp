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
#include <gravity/solver.h>
//#include "Branch_Bound.h"
#include "IPH.h"
//#include "Lower_Bound.h"
//#include "Lidar_preprocess.h"
#include "BB.h"
//#ifdef USE_CGAL
//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Polyhedron_3.h>
//#include <CGAL/convex_hull_3.h>
//#endif
//#ifndef CGAL_HAS_THREADS
//#define CGAL_HAS_THREADS
//#endif

#ifdef USE_EIGEN3
//#include </Users/smitha/Utils/eigen-3.3.9/Eigen/Dense>
#include <Eigen/Dense>
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
#include <iostream>
#include <algorithm>

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
#ifdef USE_MPI
    auto err_init = MPI_Init(nullptr,nullptr);
    int worker_id, nb_workers;
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
    auto err_size = MPI_Comm_size(MPI_COMM_WORLD, &nb_workers);
#endif
    string prob_type = "Reg";
    if(argc>1){
        prob_type = argv[1];
    }
    string Model_file = string(prj_dir)+"/data_sets/LiDAR/point_cloud1.txt";
    string Data_file = string(prj_dir)+"/data_sets/LiDAR/point_cloud1.txt";
    string red_Model_file = string(prj_dir)+"/data_sets/LiDAR/red_point_cloud1.txt";
    string red_Data_file = string(prj_dir)+"/data_sets/LiDAR/red_point_cloud1.txt";
    string algo = "ARMO", global_str = "local";
    
    if(argc>2){
        Model_file = argv[2];
    }
    if(argc>3){
        Data_file = argv[3];
    }
    if(argc>4){
        red_Model_file = argv[4];
        
    }
    if(argc>5){
        red_Data_file = argv[5];
        
    }
#ifdef USE_MPI
    auto err_init = MPI_Init(nullptr,nullptr);
    int worker_id, nb_workers;
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
    auto err_size = MPI_Comm_size(MPI_COMM_WORLD, &nb_workers);
#endif
    
    bool Registration = prob_type=="Reg";/* Solve the Registration problem */
    bool skip_first_line = true; /* First line in Go-ICP input files can be ignored */
    bool if_plot=false;
    
    vector<vector<double>> full_point_cloud_model, full_point_cloud_data;
    vector<vector<double>> point_cloud_model, point_cloud_data;
    vector<vector<double>> point_cloud_model1, point_cloud_data1;
    vector<vector<double>> full_uav_model, full_uav_data;
    vector<vector<double>> uav_model, uav_data;
    vector<vector<double>> uav_model1, uav_data1;
    vector<vector<double>> full_rpy_model, full_rpy_data;
    vector<vector<double>> rpy_model, rpy_data;
    vector<vector<double>> rpy_model1, rpy_data1;

    
    if(if_plot){
        rapidcsv::Document  red_Model_doc(red_Model_file, rapidcsv::LabelParams(0, -1),rapidcsv::SeparatorParams(' '));
        read_data(red_Model_doc, point_cloud_model, uav_model);
        rapidcsv::Document  red_Data_doc(red_Data_file, rapidcsv::LabelParams(0, -1),rapidcsv::SeparatorParams(' '));
        read_data(red_Data_doc, point_cloud_data, uav_data);
        rapidcsv::Document  Model_doc(Model_file, rapidcsv::LabelParams(0, -1),rapidcsv::SeparatorParams(' '));
        rapidcsv::Document  Data_doc(Data_file, rapidcsv::LabelParams(0, -1),rapidcsv::SeparatorParams(' '));
        auto resm=read_data(Model_doc, full_point_cloud_model, full_uav_model);
        auto resd=read_data(Data_doc, full_point_cloud_data, full_uav_data);
        
        
        
        double roll_deg=-1.56249980053;
        double pitch_deg=-0.688109578051;
        double yaw_deg=1.12496586489;
        
        DebugOn("Angle in deg roll "<<roll_deg<<endl);
        DebugOn("Angle in deg pitch "<<pitch_deg<<endl);
        DebugOn("Angle in deg yaw "<<yaw_deg<<endl);
        
        
        
        apply_rotation(roll_deg, pitch_deg, yaw_deg, full_point_cloud_model, full_point_cloud_data, full_uav_model, full_uav_data);
        
        DebugOff("Initial L2 error = " << L2error_init << endl);
        DebugOff("Initial L1 error = " << L1error_init << endl);
        
        double final_roll=roll_deg;
        double final_pitch=pitch_deg;
        double final_yaw=yaw_deg;
        
        
        
        bool save_file = true;
        if(save_file){
            auto name = Model_file.substr(0,Model_file.find('.'));
            auto fname = name+"_ARMO_RPY_"+to_string(final_roll)+"_"+to_string(final_pitch)+"_"+to_string(final_yaw)+".laz";
            save_laz(fname,full_point_cloud_data, full_point_cloud_model);
        }
        return 0;
    }
    
    /* Boresight Alignment Problem */
    
    //string file_u="/Users/smitha/Desktop/LiDAR_data/DAG4_L_2__2019_06_20_18_combined_frames701_761_0_original_2fold.laz";
    // string file_u="/Users/smitha/Desktop/LiDAR_data/14.ABSOLUTE.DAG4_L_2__2019_06_20_18_combine.20211013.RPY.1.45.0.9.-0.25.frames.712-761.1189-1267.fCP.LE.adc.laz";
    //string file_u="/Users/smitha/Desktop/LiDAR_data/Truck_only.laz";
    //string file_u="/Users/smitha/Desktop/LiDAR_data/DAG4_L_2__2019_06_20_18_combined_3line_1.45_0.9_-0.25.las";
    //string file_u1="/Users/smitha/Desktop/LiDAR_data/DAG4_L_2__2019_06_20_18_combined_frames_701_763_1.45_0.9_-0.25.las";
    // string file_u1="/Users/smitha/Desktop/LiDAR_data/DAG4_L_2__2019_06_20_18_combined_frames701_761_0.0.0.las";
    //string file_u2="/Users/smitha/Desktop/LiDAR_data/DAG4_L_2__2019_06_20_18_combined_frames_1191_1281_0.0.0.las";
    string file_u1="/Users/smitha/Desktop/LiDAR_data/5b.ABSOLUTE.DAG4_L_2__2019_06_20_18_combined.20211117.RPY.000.frames.712-761.1189-1267.fCP.LE.adc.laz";
    //string file_u1="/Users/smitha/Downloads/Centennial_USA_2015_02_19_14_00_36_frames_10_to_end_2594.laz";
    //string file_u1="/Users/smitha/Downloads/pound_farm_flight4_all_frames_1_to_5672.laz";
    //string file_u1="/Users/smitha/Desktop/LiDAR_data/14b.ABSOLUTE.DAG4_L_2__2019_06_20_18_combined.20211117.RPY.1.45.0.9.-0.25.frames.712-761.1189-1267.fCP.LE.adc.laz";
   // string file_u1="/Users/smitha/Desktop/LiDAR_data/5.ABSOLUTE.DAG4_L_2__2019_06_20_18_combined.20211008.RPY.000.frames.712-761.1189-1267.fCP.LE.adc.laz";
    //string file_u1="/Users/smitha/Downloads/21.SDM21Fig8Cars.ABS.TA51__2020_09_18_18_45_45.20211118.RPY.000.frames.1251-1416.1620-1778.fCP.LE.adc.laz";
    //string file_u1="/Users/smitha/Downloads/22.SDM21Fig11Tent.ABS.TL_Shrubz_6__2018_07_15_combined.20211118.RPY.000.frames.1643-1689.2616-2662.fCP.LE.adc.laz";
    //string file_u1="/Users/smitha/Desktop/LiDAR_data/Ta51_powerlines_1__2020_12_18_combined.laz";
    //string file_u="/Users/smitha/Desktop/LiDAR_data/17.ABSOLUTE.DAG4_L_2__2019_06_20_18_combined.20211107.RPY.1.45.0.9.-0.25.frames.712-761.1189-1267.fCP.LE.adc.laz";
    //string file_u2="/Users/smitha/Desktop/LiDAR_data/16.ABSOLUTE.DAG4_L_2__2019_06_20_18_combined.20211107.RPY.000.frames.712-761.1189-1267.fCP.LE.adc.laz";
    //string file_u="/Users/smitha/Desktop/LiDAR_data/7.ABSOLUTE.DAG4_L_2__2019_06_20_18_combined.20211008.RPY.000.frames.712-761.1189-1267.fLE.adc.laz";
    //    string file_u1="/Users/smitha/Downloads/14.ABSOLUTE.DAG4_L_2__2019_06_20_18_combined.20211013.RPY.1.45.0.9.-0.25.frames.712-761.1189-1267.fCP.LE.adc.laz";
    // string file_u1="/Users/smitha/Downloads/15.ABSOLUTE.DAG4_L_2__2019_06_20_18_combined.20211020.RPY_1.460181_0.431533_-0.06725.frames.712-761.1189-1267.fCP.LE.adc.laz";
    //string file_u1="/Users/smitha/Downloads/14.ABSOLUTE.DAG4_L_2__2019_06_20_18_combined.20211013.RPY.1.45.0.9.-0.25.frames.712-761.1189-1267.fCP.LE.adc.laz";
    vector<vector<double>> lidar_point_cloud1, lidar_point_cloud2;
    vector<vector<double>> roll_pitch_yaw1, roll_pitch_yaw2;
    //read_laz_new(file_u1, file_u2);
    
    //   auto ucloud=::read_laz(file_u1, lidar_point_cloud1, roll_pitch_yaw1);
    auto uav_cloud_u1=::read_laz(file_u1, lidar_point_cloud1, roll_pitch_yaw1);
    // auto uav_cloud_u2=::read_laz(file_u2, lidar_point_cloud2, roll_pitch_yaw2);
    
    vector<vector<double>> em4;
    int n=roll_pitch_yaw1.size();
    double scanner_x=0.0, scanner_y=0.160999998450279, scanner_z=0.016;
    const double hr=0.000288477458525449, hp=0.000563974492251873, hy=0.00173089664895087;
//    auto roll_deg=-1.44922*pi/180;
//    auto pitch_deg=0.542969*pi/180;
//    auto yaw_deg=0.00390625*pi/180;
    
    //    auto roll_deg =-1.45*pi/180;
    //    auto pitch_deg = 0.899999999*pi/180;
    //    auto yaw_deg = -0.25*pi/180;
    //    auto roll_deg= -0.0254859253764153;
    //    auto pitch_deg =0.0153557369485497;
    //    auto yaw_deg =-0.00552612589672208;
//    auto roll_deg=-0.0255034808069468;
//    auto pitch_deg = 0.0152935786172748;
//    auto yaw_deg =-0.00570988887920976;
            auto roll_deg =-1.42578*pi/180;
            auto pitch_deg = 0.621094*pi/180;
            auto yaw_deg = 0.113281*pi/180;
    //apply_transform_new_order_Test(roll_deg, pitch_deg, yaw_deg, lidar_point_cloud, uav_cloud_u, roll_pitch_yaw, scanner_x, scanner_y, scanner_z);
    //apply_transform_new_order(roll_deg, pitch_deg, yaw_deg, lidar_point_cloud1, uav_cloud_u1, roll_pitch_yaw1, scanner_x, scanner_y, scanner_z, hr, hp, hy);
    // apply_transform_new_order_Test(roll_deg, pitch_deg, yaw_deg, lidar_point_cloud1, uav_cloud_u1, roll_pitch_yaw1, scanner_x, scanner_y, scanner_z);
    
    //save_laz(file_u1.substr(0,Model_file.find('.'))+to_string(scanner_x)+"_"+to_string(scanner_y)+"_"+to_string(scanner_z)+"_manual_lv_5_.laz", lidar_point_cloud2, em4);
    //save_laz(file_u1.substr(0,Model_file.find('.'))+to_string(roll_deg)+"_"+to_string(pitch_deg)+"_"+to_string(yaw_deg)+"_pf_scanner.laz", lidar_point_cloud1, em4);
    
    
    
    auto uav_cloud_u=uav_cloud_u1;
    auto lidar_point_cloud=lidar_point_cloud1;
    auto roll_pitch_yaw=roll_pitch_yaw1;
    auto file_u=file_u1;
    
    double uxmin=numeric_limits<double>::max(), uxmax=-uxmin, uymin=numeric_limits<double>::max(), uymax=-uymin, uzmin=numeric_limits<double>::max(), uzmax=-uzmin;
    for(auto i=0;i<uav_cloud_u.size();i++)
    {
        auto x=uav_cloud_u.at(i)[0];
        auto y=uav_cloud_u.at(i)[1];
        auto z=uav_cloud_u.at(i)[2];
        if(x<=uxmin){
            uxmin=x;
            
        }
        if(y<=uymin){
            uymin=y;
            
        }
        if(z<=uzmin){
            uzmin=z;
            
        }
        if(x>=uxmax){
            uxmax=x;
           
        }
        if(y>=uymax){
            uymax=y;
            
        }
        if(z>=uzmax){
            uzmax=z;
            
        }
        
    }
        DebugOn("umin max in x "<<(uxmin)<<" "<<(uxmax)<<endl);
        double diff_x=uxmax-uxmin;
        DebugOn("umin max in y "<<uymin<<" "<<uymax<<endl);
        double diff_y=uymax-uymin;
        DebugOn("umin max in z "<<uzmin<<" "<<uzmax<<endl);
        double diff_z=uzmax-uzmin;
        
        DebugOn(diff_x<<" "<<diff_y<<" "<<diff_z<<endl);
    vector<vector<double>> em3;
    
    uxmin=numeric_limits<double>::max(); uxmax=-uxmin; uymin=numeric_limits<double>::max(); uymax=-uymin; uzmin=numeric_limits<double>::max(); uzmax=-uzmin;
    vector<vector<double>> empty_vec, uav_xy;
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
   
    
    uxmin=numeric_limits<double>::max(); uxmax=-uxmin; uymin=numeric_limits<double>::max(); uymax=-uymin; uzmin=numeric_limits<double>::max(); uzmax=-uzmin;
    vector<vector<double>> rough_pc;
    auto index_set=filter_z_slope(uav_xy);
    empty_vec.push_back({uxmin, uymin, 0});
    //plot(uav_xy, empty_vec);
    vector<int> turns(0);
    auto slices_indices=extract_slices(uav_xy, index_set, turns);
    vector<vector<vector<double>>> uplot_array;
    vector<vector<vector<double>>> slice_array;
    vector<int> turn_array;
    vector<vector<int>> ulist_array;
    multimap<double, int, greater<double>> rank_map;
    multimap<double, int, greater<double>> z_map_model;
    multimap<double, int, greater<double>> z_map_data;
    for(auto i=0;i<slices_indices.size();i++){
        if(slices_indices[i].second-slices_indices[i].first>10){
            auto ulist=reg_slope_lines(uav_xy, slices_indices[i], turns[i]);
            vector<vector<double>> slice_plot, u_plot;
            for(auto j=slices_indices[i].first;j<=slices_indices[i].second;j++){
                auto x=uav_xy.at(j)[0]+uav_cloud_u.at(0)[0];
                auto y=uav_xy.at(j)[1]+uav_cloud_u.at(0)[1];
                auto z=uav_xy.at(j)[2]+uav_cloud_u.at(0)[2];
                slice_plot.push_back({x,y,z});
            }
            //plot(slice_plot, uav_coords, 3);
            if(ulist.size()>=4){
                DebugOn("ulist "<<endl);
                for(auto j=0;j<ulist.size();j++){
                    auto x=uav_xy[ulist[j]][0]+uav_cloud_u.at(0)[0];
                    auto y=uav_xy[ulist[j]][1]+uav_cloud_u.at(0)[1];
                    auto z=uav_xy[ulist[j]][2]+uav_cloud_u.at(0)[2];
                    u_plot.push_back({x,y,z});
                    DebugOn(x<<" "<<y<<" "<<z<<endl);
                }
                auto x1=uav_xy[ulist[0]][0];
                auto y1=uav_xy[ulist[0]][1];
                auto z1=uav_xy[ulist[0]][2];
                auto x2=uav_xy[ulist[1]][0];
                auto y2=uav_xy[ulist[1]][1];
                auto z2=uav_xy[ulist[1]][2];
                auto x3=uav_xy[ulist[2]][0];
                auto y3=uav_xy[ulist[2]][1];
                auto z3=uav_xy[ulist[2]][2];
                auto x4=uav_xy[ulist[3]][0];
                auto y4=uav_xy[ulist[3]][1];
                auto z4=uav_xy[ulist[3]][2];
                double s1=(y2-y1)/(x2-x1);
                double s2=(y4-y3)/(x4-x3);
                double d1=pow(y2-y1,2)+pow(x2-x1,2);
                double d2=pow(y4-y3,2)+pow(x4-x3,2);
                DebugOn("d1 "<<d1<<" d2 "<<d2<<endl);
                DebugOn("s1 "<<s1<<" s2 "<<s2<<endl);
                uplot_array.push_back(u_plot);
                ulist_array.push_back(ulist);
                slice_array.push_back(slice_plot);
                turn_array.push_back(turns[i]);
                double score=(d1+d2)/(1e5)-abs(s1-s2)/(abs(s1+s2));
                DebugOn("score "<<score<<endl);
                rank_map.insert(pair<double, int>(score, slice_array.size()-1));
                //plot(u_plot, slice_plot, uav_coords,3);
            }
            empty_vec.clear();
        }
    }
    vector<vector<double>> cloud1, cloud2, uav1, uav2, pc1, pc2, u1, u2, rpy1, rpy2;
    vector<vector<double>> em;
    if(rank_map.size()>=1){
        auto it=rank_map.begin();
        int pos=it->second;
        //plot(uplot_array[pos], uav_coords,3);
        for(auto j=0;j<uplot_array[pos].size();j++){
            auto x=uplot_array[pos].at(j)[0];
            auto y=uplot_array[pos].at(j)[1];
            auto z=uplot_array[pos].at(j)[2];
            DebugOff(x<<" "<<y<<" "<<z<<endl);
        }
        auto uplot=uplot_array[pos];
        auto x1=uplot.at(0)[0];
        auto y1=uplot.at(0)[1];
        auto z1=uplot.at(0)[2];
        auto x2=uplot.at(1)[0];
        auto y2=uplot.at(1)[1];
        auto z2=uplot.at(1)[2];
        auto x3=uplot.at(2)[0];
        auto y3=uplot.at(2)[1];
        auto z3=uplot.at(2)[2];
        auto x4=uplot.at(3)[0];
        auto y4=uplot.at(3)[1];
        auto z4=uplot.at(3)[2];
        double s1=(y2-y1)/(x2-x1);
        double s2=(y4-y3)/(x4-x3);
        double d1=pow(y2-y1,2)+pow(x2-x1,2);
        double d2=pow(y4-y3,2)+pow(x4-x3,2);
        double score=(d1+d2)/100-abs(s1-s2);
        DebugOff("d1 "<<d1<<" d2 "<<d2<<endl);
        DebugOff("s1 "<<s1<<" s2 "<<s2<<endl);
        auto ulist=ulist_array[pos];
        vector<int> frame1(0), frame2(0);
        
        for(auto i=ulist[0];i<ulist[1];i++){
            cloud1.push_back(lidar_point_cloud.at(i));
            uav1.push_back(uav_cloud_u.at(i));
            rpy1.push_back(roll_pitch_yaw.at(i));
        }
        for(auto i=ulist[2];i<ulist[3];i++){
            cloud2.push_back(lidar_point_cloud.at(i));
            uav2.push_back(uav_cloud_u.at(i));
            rpy2.push_back(roll_pitch_yaw.at(i));
        }
        DebugOn("cloud1.size() "<<cloud1.size()<<endl);
        DebugOn("cloud2.size() "<<cloud2.size()<<endl);
        auto name_u=file_u.substr(0,Model_file.find('.'))+"_split.laz";
        save_laz(name_u,cloud1, cloud2);
        auto name_u_uav=file_u.substr(0,Model_file.find('.'))+"_uav.laz";
        save_laz(name_u_uav,cloud1, uav1);
        //plot(cloud1,cloud2);
        if(cloud1.size()>=cloud2.size()){
            full_point_cloud_model=cloud1;
            full_point_cloud_data=cloud2;
            full_uav_model=uav1;
            full_uav_data=uav2;
            full_rpy_model=rpy1;
            full_rpy_data=rpy2;
        }
        else{
            full_point_cloud_model=cloud2;
            full_point_cloud_data=cloud1;
            full_uav_model=uav2;
            full_uav_data=uav1;
            full_rpy_model=rpy2;
            full_rpy_data=rpy1;
        }
        
        uxmin=numeric_limits<double>::max(); uxmax=-uxmin; uymin=numeric_limits<double>::max(); uymax=-uymin; uzmin=numeric_limits<double>::max(); uzmax=-uzmin;
        for(auto i=0;i<full_point_cloud_model.size();i++){
            auto x=full_point_cloud_model.at(i)[0];
            auto y=full_point_cloud_model.at(i)[1];
            auto z=full_point_cloud_model.at(i)[2];
            if(x<=uxmin){
                uxmin=x;
            }
            if(y<=uymin){
                uymin=y;
            }
            if(z<=uzmin){
                uzmin=z;
            }
            if(x>=uxmax){
                uxmax=x;
            }
            if(y>=uymax){
                uymax=y;
            }
            if(z>=uzmax){
                uzmax=z;
            }
        }
        for(auto i=0;i<full_point_cloud_model.size();i++){
            auto x=full_point_cloud_model.at(i)[0];
            auto y=full_point_cloud_model.at(i)[1];
            auto z=full_point_cloud_model.at(i)[2];
            double lamdax=(x-uxmin)/(uxmax-uxmin);
            double lamday=(y-uymin)/(uymax-uymin);
            if(lamdax>=0.28 && lamdax<=0.325){//0.3
                if(lamday>=0.45 && lamday<=0.52){
                    if(z>=1261){
                        if((lamdax>0.3 && z>=1261.6) || lamdax<=0.3){
                    point_cloud_model1.push_back(full_point_cloud_model.at(i));
                    uav_model1.push_back(full_uav_model.at(i));
                    rpy_model1.push_back(full_rpy_model.at(i));
                        }
                    }
                }
            }
        }
        uxmin=numeric_limits<double>::max(); uxmax=-uxmin; uymin=numeric_limits<double>::max(); uymax=-uymin; uzmin=numeric_limits<double>::max(); uzmax=-uzmin;
        for(auto i=0;i<full_point_cloud_data.size();i++){
            auto x=full_point_cloud_data.at(i)[0];
            auto y=full_point_cloud_data.at(i)[1];
            auto z=full_point_cloud_data.at(i)[2];
            if(x<=uxmin){
                uxmin=x;
            }
            if(y<=uymin){
                uymin=y;
            }
            if(z<=uzmin){
                uzmin=z;
            }
            if(x>=uxmax){
                uxmax=x;
            }
            if(y>=uymax){
                uymax=y;
            }
            if(z>=uzmax){
                uzmax=z;
            }
        }
        for(auto i=0;i<full_point_cloud_data.size();i++){
            auto x=full_point_cloud_data.at(i)[0];
            auto y=full_point_cloud_data.at(i)[1];
            auto z=full_point_cloud_data.at(i)[2];
            double lamdax=(x-uxmin)/(uxmax-uxmin);
            double lamday=(y-uymin)/(uymax-uymin);
            if(lamdax>=0.45 && lamdax<=0.5){//0.5
                if(lamday>=0.4 && lamday<=0.46){
                    if(z>=1262.5){
                        if((lamdax>0.48 && z>=1262.7) || lamdax<=0.48){
                    point_cloud_data1.push_back(full_point_cloud_data.at(i));
                    uav_data1.push_back(full_uav_data.at(i));
                    rpy_data1.push_back(full_rpy_data.at(i));
                        }
                    }
                }
            }
        }
        DebugOn(point_cloud_model1.size()<<endl);
        DebugOn(point_cloud_data1.size()<<endl);
#ifdef USE_MATPLOT
        plot(point_cloud_model1, point_cloud_data1);
#endif
         save_laz(file_u.substr(0,Model_file.find('.'))+"_model1.laz", point_cloud_model1, em);
         save_laz(file_u.substr(0,Model_file.find('.'))+"_data1.laz", point_cloud_data1, em);
        int count=0;
        vector<vector<double>> point_cloud_model_temp, uav_model_temp, rpy_model_temp;
       
        for(auto i=0;i<point_cloud_model1.size();i+=1){
            point_cloud_data.push_back(point_cloud_model1.at(i));
            uav_data.push_back(uav_model1.at(i));
            rpy_data.push_back(rpy_model1.at(i));
            
        }
        vector<vector<double>> em, point_cloud_data_temp, uav_data_temp, rpy_data_temp;
        // save_laz(file_u.substr(0,Model_file.find('.'))+"_model.laz", point_cloud_model, em);
        count=0;
     
        for(auto i=0;i<point_cloud_data1.size();i+=2){
            point_cloud_model.push_back(point_cloud_data1.at(i));
            uav_model.push_back(uav_data1.at(i));
            rpy_model.push_back(rpy_data1.at(i));
            
        }
        
        double x_scale=0, y_scale=0, z_scale=0;
        double xmin=point_cloud_model.at(0)[0], ymin=point_cloud_model.at(0)[1], zmin=point_cloud_model.at(0)[2];
        bool scale=true;
        if(scale){
            x_scale=uxmin;
            y_scale=uymin;
            z_scale=uzmin;
            for(auto i=0;i<point_cloud_model.size();i++){
                point_cloud_model.at(i)[0]-=uav_cloud_u.at(0)[0];
                point_cloud_model.at(i)[1]-=uav_cloud_u.at(0)[1];
                point_cloud_model.at(i)[2]-=uav_cloud_u.at(0)[2];
            }
            for(auto i=0;i<point_cloud_data.size();i++){
                point_cloud_data.at(i)[0]-=uav_cloud_u.at(0)[0];
                point_cloud_data.at(i)[1]-=uav_cloud_u.at(0)[1];
                point_cloud_data.at(i)[2]-=uav_cloud_u.at(0)[2];
            }
            for(auto i=0;i<uav_model.size();i++){
                uav_model.at(i)[0]-=uav_cloud_u.at(0)[0];
                uav_model.at(i)[1]-=uav_cloud_u.at(0)[1];
                uav_model.at(i)[2]-=uav_cloud_u.at(0)[2];
            }
            for(auto i=0;i<uav_data.size();i++){
                uav_data.at(i)[0]-=uav_cloud_u.at(0)[0];
                uav_data.at(i)[1]-=uav_cloud_u.at(0)[1];
                uav_data.at(i)[2]-=uav_cloud_u.at(0)[2];
            }
        }
        DebugOn(point_cloud_model.size()<<endl);
        DebugOn(point_cloud_data.size()<<endl);
#ifdef USE_MATPLOT
        plot(point_cloud_model, point_cloud_data, 1);
#endif
        /*plot(uav_model, uav_data, 1);
        plot(rpy_model, rpy_data, 1);*/
        
        indices N1 = range(1,point_cloud_data.size());
        indices N2 = range(1,point_cloud_model.size());
        //auto valid_cells_old = indices(N1,N2);
        indices valid_cells_old("valid_cells_old");
        DebugOn("valid cells old size "<<valid_cells_old.size()<<endl);
        indices new_cells("new_cells");
        indices new_cells1("new_cells1");
        param<double> dist_cells("dist_cells");
        param<double> dist_cells1("dist_cells1");
        double upper_bound=3e4;
        double prep_time=0;
        
        double roll_min=-2*pi/180;
        double roll_max=2*pi/180;
        double pitch_min=-2*pi/180;
        double pitch_max=2*pi/180;
        double yaw_min=-2*pi/180;
        double yaw_max=2*pi/180;
        
        auto point_cloud_model_copy=point_cloud_model;
        auto point_cloud_data_copy=point_cloud_data;
        
        vector<int> matching(point_cloud_data.size());
        vector<double> err_per_point(point_cloud_data.size());
        
        auto L2init=computeL2error(point_cloud_model,point_cloud_data,matching,err_per_point);
        auto L1init=computeL1error(point_cloud_model,point_cloud_data,matching,err_per_point);
        
        DebugOn("L2init  "<<L2init<<endl);
        DebugOn("L1init  "<<L1init<<endl);
        
        string error_type="L2";
        bool algo_IPH=false;
        
        vector<double> best_rot_trans(9,0.0);
        double best_ub=1e5;
        if(error_type=="L2"){
            best_ub=L2init;
        }
        else{
            best_ub=L1init;
        }
        if(algo_IPH){
            auto pcm=point_cloud_model;
            auto pcd=point_cloud_data;
            auto uavm=uav_model;
            auto uavd=uav_data;
            
            auto pcm1=point_cloud_model;
            auto pcd1=point_cloud_data;
            auto uavm1=uav_model;
            auto uavd1=uav_data;
            
            auto resi=run_IPH(pcm, pcd, uavm, uavd, rpy_model, rpy_data, scanner_x, scanner_y, scanner_z, hr, hp, hy);
            
            double roll_deg_i=get<0>(resi);
            double pitch_deg_i=get<1>(resi);
            double yaw_deg_i=get<2>(resi);
            double roll_rad_i=roll_deg_i*pi/180;
            double pitch_rad_i=pitch_deg_i*pi/180;
            double yaw_rad_i=yaw_deg_i*pi/180;
            double errori;
            
            apply_transform_new_order(roll_rad_i, pitch_rad_i, yaw_rad_i, pcm1, uavm1, rpy_model, scanner_x,scanner_y,scanner_z ,hr,hp,hy);
            apply_transform_new_order(roll_rad_i, pitch_rad_i, yaw_rad_i, pcd1, uavd1, rpy_data, scanner_x,scanner_y,scanner_z,hr,hp,hy);
            apply_transform_new_order(roll_rad_i, pitch_rad_i, yaw_rad_i, lidar_point_cloud, uav_cloud_u, roll_pitch_yaw, scanner_x,scanner_y,scanner_z,hr,hp,hy);
            
            auto errorl2= computeL2error(pcm1,pcd1,matching,err_per_point);
            errori= computeL1error(pcm1,pcd1,matching,err_per_point);
            DebugOn("error L1 "<<errori<<endl);
            DebugOn("error L2 "<<errorl2<<endl);
            DebugOn("Percentage improved L1 "<<(L1init-errori)/L1init*100.0<<endl);
            save_laz(file_u.substr(0,Model_file.find('.'))+"_"+to_string(roll_deg_i)+"_"+to_string(pitch_deg_i)+"_"+to_string(yaw_deg_i)+"_full_set.laz", lidar_point_cloud, em);
        }
        else{
            auto rot= BranchBound_Align(point_cloud_model, point_cloud_data, uav_model, uav_data, rpy_model, rpy_data, best_rot_trans, best_ub, error_type, scanner_x, scanner_y, scanner_z, hr, hp, hy);
    
           
            auto roll_deg_bb = rot[0];
            auto pitch_deg_bb = rot[1];
            auto yaw_deg_bb = rot[2];
            
            apply_transform_new_order(roll_deg_bb*pi/180, pitch_deg_bb*pi/180, yaw_deg_bb*pi/180, point_cloud_model, uav_model, rpy_model, scanner_x,scanner_y,scanner_z,hr,hp,hy);
            apply_transform_new_order(roll_deg_bb*pi/180, pitch_deg_bb*pi/180, yaw_deg_bb*pi/180, point_cloud_data, uav_data, rpy_data, scanner_x,scanner_y,scanner_z,hr,hp,hy);
            save_laz(file_u.substr(0,Model_file.find('.'))+"_"+to_string(roll_deg_bb)+"_"+to_string(pitch_deg_bb)+"_"+to_string(yaw_deg_bb)+"_opt_set.laz", point_cloud_model, point_cloud_data);
            apply_transform_new_order(roll_deg_bb*pi/180, pitch_deg_bb*pi/180, yaw_deg_bb*pi/180, lidar_point_cloud, uav_cloud_u, roll_pitch_yaw, scanner_x,scanner_y,scanner_z,hr,hp,hy);
            save_laz(file_u.substr(0,Model_file.find('.'))+"_"+to_string(roll_deg_bb)+"_"+to_string(pitch_deg_bb)+"_"+to_string(yaw_deg_bb)+"_full_set.laz", lidar_point_cloud, em);
            
            auto L2=computeL2error(point_cloud_model,point_cloud_data,matching,err_per_point);
            auto L1=computeL1error(point_cloud_model,point_cloud_data,matching,err_per_point);
            
            DebugOn("L2  "<<L2<<endl);
            DebugOn("L1  "<<L1<<endl);
            
            DebugOn("Percentage improved L2 "<<(L2init-L2)/L2init*100.0<<endl);
            DebugOn("Percentage improved L1 "<<(L1init-L1)/L1init*100.0<<endl);
        }
#ifdef USE_MATPLOT
        plot(point_cloud_model, point_cloud_data);
#endif
    }
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
