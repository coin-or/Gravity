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
#include <gravity/matplotlibcpp.h>
#include <DataSet.h>
#include "lasreader.hpp"
#include "laswriter.hpp"
#include <gravity/KDTreeVectorOfVectorsAdaptor.h>
#include <time.h>
using namespace std;
#include <gravity/jly_goicp.h>
#include <gravity/ConfigMap.hpp>

#define DEFAULT_OUTPUT_FNAME "output.txt"
#define DEFAULT_CONFIG_FNAME "config.txt"
#define DEFAULT_MODEL_FNAME "model.txt"
#define DEFAULT_DATA_FNAME "data.txt"

/* Set the different options for GoICP  */
void set_GoICP_options(GoICP & goicp);

/* Centralize point cloud around origin */
void centralize(int n, POINT3D **  p, double avg_x, double avg_y, double avg_z);

/* Scale point clouds to [-1,1] */
void unit_scale(int n1, POINT3D **  p1, int n2, POINT3D **  p2);

/* Scale point cloud using provided max values */
void scale_all(int n1, POINT3D **  p1, double max_x, double max_y, double max_z, double min_x, double min_y, double min_z);

/* Computes the interpolation coefficient based on time elapsed  */
double get_interpolation_coef(const double& lidar_time, UAVPoint* p1, UAVPoint* p2);

/* Return the min-max values for x, y and z  for all possible rotations of p with angle +- angle*/
vector<pair<double,double>> get_min_max(double angle, const vector<double>& p, const vector<double>& ref);

/* Return true if two cubes intersect
 The cube is stored using a vector of size 3: {x,y,z}, where each entry is [min,max] on the corresponding axis
 */
bool intersect(const vector<pair<double,double>>& a, const vector<pair<double,double>>& b);

/* Returns the coordinates of the cube center
 The cube is stored using a vector of size 3: {x,y,z}, where each entry is [min,max] on the corresponding axis
 */
tuple<double,double,double> get_center(const vector<pair<double,double>>& cube);

/* Apply rotation + translation on input data (Registration) */
void apply_rot_trans(double roll, double pitch, double yaw, double x_shift, double y_shift, double z_shift, vector<vector<double>>& point_cloud);

/* Return vector of 6 extreme points (2 per axis) from point cloud */
vector<vector<double>> get_extreme_points(const vector<vector<double>>& point_cloud);

/* Return vector of n extreme points from point cloud */
vector<vector<double>> get_n_extreme_points(int n, const vector<vector<double>>& point_cloud);

/* Run iterative heuristic */
tuple<double,double,double,double,double,double> run_heuristic(const vector<vector<double>>& ext_model, vector<vector<double>>& ext_data, vector<vector<double>>& point_cloud_data);

/* Compute the L2 error */
double computeL2error(const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data);

/* Return central point from point cloud */
vector<double> get_center(const vector<vector<double>>& point_cloud);


/* Apply rotation on input data (Boresight alignment) */
void apply_rotation(double roll, double pitch, double yaw, vector<double>& x_vec1, vector<double>& y_vec1, vector<double>& z_vec1, vector<double>& x_vec2, vector<double>& y_vec2, vector<double>& z_vec2, const vector<double>& x_shift1, const vector<double>& y_shift1, const vector<double>& z_shift1, const vector<double>& x_shift2, const vector<double>& y_shift2, const vector<double>& z_shift2);

/* Run the ARMO model for boresight alignment */
tuple<double,double,double> run_ARMO(string axis, const vector<double>& x_vec1, const vector<double>& y_vec1, const vector<double>& z_vec1, const vector<double>& x_vec2, const vector<double>& y_vec2, const vector<double>& z_vec2, const vector<double>& x_shift1, const vector<double>& y_shift1, const vector<double>& z_shift1, const vector<double>& x_shift2, const vector<double>& y_shift2, const vector<double>& z_shift2, const vector<vector<double>>& point_cloud1, const vector<vector<double>>& point_cloud2, const vector<vector<double>>& uav1, const vector<vector<double>>& uav2);

/* Plot two point clouds */
void plot(const vector<vector<double>>& point_cloud1, const vector<vector<double>>& point_cloud2, double point_thick = 0.1);

/* Run Go-ICP on point clouds, return best solution */
tuple<double,double,double,double,double,double> run_GoICP(const vector<vector<double>>& point_cloud1, const vector<vector<double>>& point_cloud2);

/* Run the ARMO model for registration */
tuple<double,double,double,double,double,double> run_ARMO(bool bypass, string axis, const vector<vector<double>>& point_cloud1, const vector<vector<double>>& point_cloud2);

/* Run the Global ARMO model for registration */
tuple<double,double,double,double,double,double> run_ARMO_Global(bool bypass, string axis, const vector<vector<double>>& point_cloud1, const vector<vector<double>>& point_cloud2);

/* Run the MINLP ARMO model for registration */
tuple<double,double,double,double,double,double> run_ARMO_MINLP(bool bypass, string axis, const vector<vector<double>>& point_cloud1, const vector<vector<double>>& point_cloud2);



int main (int argc, char * argv[])
{
    bool test_global = false;
    if (test_global) {
        vector<vector<double>> point_cloud_model, point_cloud_data;
        point_cloud_model.resize(6);
        for(int i = 0; i<6; i++){
            point_cloud_model[i].resize(3);
        }
        
        point_cloud_model[0][0] = 0.1;
        point_cloud_model[0][1] = 0.1;
        point_cloud_model[0][2] = 0.1;
        
        point_cloud_model[1][0] = 0.1;
        point_cloud_model[1][1] = 0.9;
        point_cloud_model[1][2] = 0.1;
        
        point_cloud_model[2][0] = 0.9;
        point_cloud_model[2][1] = 0.1;
        point_cloud_model[2][2] = 0.1;
        
        point_cloud_model[3][0] = 0.9;
        point_cloud_model[3][1] = 0.9;
        point_cloud_model[3][2] = 0.1;
        
        point_cloud_model[4][0] = 0.5;
        point_cloud_model[4][1] = 0.5;
        point_cloud_model[4][2] = 0.9;
        
        point_cloud_model[5][0] = 0.5;
        point_cloud_model[5][1] = 0.5;
        point_cloud_model[5][2] = 0.8;
        
        point_cloud_data = point_cloud_model;
        apply_rot_trans(40, -25, 25, 0.1, 0.04, 0.2, point_cloud_data);
        plot(point_cloud_model,point_cloud_data,2);
//        auto res = run_GoICP(point_cloud_model,point_cloud_data);
//        auto res = run_ARMO_Global(false, "full", point_cloud_model, point_cloud_data);
        int nb_iter = 0, max_nb_iter = 10;
        double roll = 0, pitch = 0, yaw = 0, x_shift = 0, y_shift = 0, z_shift = 1;
        while(nb_iter < max_nb_iter && std::abs(roll)+std::abs(pitch)+std::abs(yaw)+std::abs(x_shift)+std::abs(y_shift)+std::abs(z_shift)>1e-1) {
            auto L2error = computeL2error(point_cloud_model,point_cloud_data);
            DebugOn("L2 error before = " << L2error << endl);
            auto res = run_ARMO(false, "full", point_cloud_model, point_cloud_data);
            roll = get<0>(res); pitch = get<1>(res); yaw = get<2>(res); x_shift = get<3>(res); y_shift = get<4>(res); z_shift = get<5>(res);
            apply_rot_trans(roll, pitch, yaw, x_shift, y_shift, z_shift, point_cloud_data);
            L2error = computeL2error(point_cloud_model,point_cloud_data);
            DebugOn("L2 error after = " << L2error << endl);
            res = run_ARMO(false, "z", point_cloud_model, point_cloud_data);
            roll = get<0>(res); pitch = get<1>(res); yaw = get<2>(res); x_shift = get<3>(res); y_shift = get<4>(res); z_shift = get<5>(res);
            apply_rot_trans(roll, pitch, yaw, x_shift, y_shift, z_shift, point_cloud_data);
            L2error = computeL2error(point_cloud_model,point_cloud_data);
            DebugOn("L2 error after = " << L2error << endl);
            res = run_ARMO(false, "x", point_cloud_model, point_cloud_data);
            roll = get<0>(res); pitch = get<1>(res); yaw = get<2>(res); x_shift = get<3>(res); y_shift = get<4>(res); z_shift = get<5>(res);
            apply_rot_trans(roll, pitch, yaw, x_shift, y_shift, z_shift, point_cloud_data);
            L2error = computeL2error(point_cloud_model,point_cloud_data);
            DebugOn("L2 error after = " << L2error << endl);
            res = run_ARMO(false, "y", point_cloud_model, point_cloud_data);
            roll = get<0>(res); pitch = get<1>(res); yaw = get<2>(res); x_shift = get<3>(res); y_shift = get<4>(res); z_shift = get<5>(res);
            apply_rot_trans(roll, pitch, yaw, x_shift, y_shift, z_shift, point_cloud_data);
            L2error = computeL2error(point_cloud_model,point_cloud_data);
            DebugOn("L2 error after = " << L2error << endl);
            DebugOn("ITERATION " << nb_iter << endl);
            nb_iter++;
        }
//        auto res = run_ARMO_MINLP(false, "full", point_cloud_model, point_cloud_data);
        
        plot(point_cloud_model,point_cloud_data,2);
        return 0;
    }
    bool Registration = true;/* Solve the Registration problem */
    bool skip_first_line = true; /* First line in Go-ICP input files can be ignored */
    if(Registration){
        vector<double> x_vec0, y_vec0, z_vec0, x_vec1, y_vec1, z_vec1;
        vector<vector<double>> point_cloud_model, point_cloud_data;
        string Model_file = string(prj_dir)+"/data_sets/LiDAR/model.txt";
        string Data_file = string(prj_dir)+"/data_sets/LiDAR/data.txt";
        string algo = "ARMO", global_str = "local";
        if(argc>1){
            Model_file = argv[1];
        }
        if(argc>2){
            Data_file = argv[2];
        }
        if(argc>3){
            algo = argv[3];
        }
        if(argc>4){
            global_str = argv[4];
        }
        rapidcsv::Document  Model_doc(Model_file, rapidcsv::LabelParams(0, -1),rapidcsv::SeparatorParams(' '));
        rapidcsv::Document  Data_doc(Data_file, rapidcsv::LabelParams(0, -1),rapidcsv::SeparatorParams(' '));
        int model_nb_rows = Model_doc.GetRowCount();
        int data_nb_rows = Data_doc.GetRowCount();
        if(model_nb_rows<3){
            throw invalid_argument("Model file with less than 2 points");
            return 0;
        }
        if(data_nb_rows<3){
            throw invalid_argument("Data file with less than 2 points");
            return 0;
        }
        DebugOn("Model file has " << model_nb_rows-1 << " rows" << endl);
        DebugOn("Data file has " << data_nb_rows-1 << " rows" << endl);
        int row0 = 0;
        if(skip_first_line)
            row0 = 1;
        point_cloud_model.resize(model_nb_rows-2);
        for (int i = row0; i< model_nb_rows-1; i++) { // Input iterator
            auto x = Model_doc.GetCell<double>(0, i);
            auto y = Model_doc.GetCell<double>(1, i);
            auto z = Model_doc.GetCell<double>(2, i);
            x_vec0.push_back(x);
            y_vec0.push_back(y);
            z_vec0.push_back(z);
            point_cloud_model[i-1].resize(3);
            point_cloud_model[i-1][0] = x;
            point_cloud_model[i-1][1] = y;
            point_cloud_model[i-1][2] = z;
        }
        point_cloud_data.resize(data_nb_rows-2);
        for (int i = row0; i< data_nb_rows-1; i++) { // Input iterator
            auto x = Data_doc.GetCell<double>(0, i);
            auto y = Data_doc.GetCell<double>(1, i);
            auto z = Data_doc.GetCell<double>(2, i);
            x_vec1.push_back(x);
            y_vec1.push_back(y);
            z_vec1.push_back(z);
            point_cloud_data[i-1].resize(3);
            point_cloud_data[i-1][0] = x;
            point_cloud_data[i-1][1] = y;
            point_cloud_data[i-1][2] = z;
        }
        auto old_point_cloud = point_cloud_data;
        int nb_ext = 10;
        bool global = global_str=="global";
        bool filter_extremes = (algo!="GoICP");
        auto ext_model = point_cloud_model;
        auto ext_data = point_cloud_data;
        if (filter_extremes) {
            ext_model = get_n_extreme_points(nb_ext, point_cloud_model);
            ext_data = get_n_extreme_points(nb_ext, point_cloud_data);
        }
        
        /* Rotating and translating */
        auto thetax = atan2(-0.0081669, -0.0084357)/2;
        auto thetay = atan2(0.9999311, std::sqrt(0.0081669*0.0081669+0.0084357*0.0084357))/2;
        auto thetaz = atan2(-0.0081462,-0.0084556)/2;
        
//        apply_rot_trans(thetay*180/pi, thetax*180/pi, thetaz*180/pi, 0.2163900/2, -0.1497952/2, 0.0745708/2, ext_data);
//        apply_rot_trans(thetay*180/pi, thetax*180/pi, thetaz*180/pi, 0.2163900/2, -0.1497952/2, 0.0745708/2, point_cloud_data);

        
        bool show_extremes = false;
        vector<double> x_vec_model(ext_model.size()), y_vec_model(ext_model.size()), z_vec_model(ext_model.size());
        vector<double> x_vec_data(ext_data.size()), y_vec_data(ext_data.size()), z_vec_data(ext_data.size());
        if (show_extremes) {
            plot(ext_model,ext_data,1);
        }
        
        bool show_full_set = false;
        if (show_full_set) {
            plot(point_cloud_model,point_cloud_data);
        }
        tuple<double,double,double,double,double,double> res, res1, res2;
        bool run_goICP = (algo=="GoICP");
        if(run_goICP){/* Run GoICP inline */
            res = run_GoICP(point_cloud_model, point_cloud_data);
            auto roll = get<0>(res);auto pitch = get<1>(res);auto yaw = get<2>(res);auto x_shift = get<3>(res);auto y_shift = get<4>(res);auto z_shift = get<5>(res);
            apply_rot_trans(roll, pitch, yaw, x_shift, y_shift, z_shift, point_cloud_data);
            auto L2error = computeL2error(point_cloud_model, point_cloud_data);
            DebugOn("L2Error = " << L2error << endl);
        }
        else {
            res1 = run_heuristic(ext_model,ext_data,point_cloud_data);
            auto L2error = computeL2error(ext_model,ext_data);
            DebugOn("L2 error with exterme set after = " << L2error << endl);
            L2error = computeL2error(point_cloud_model,point_cloud_data);
            DebugOn("L2 error with full set after = " << L2error << endl);
            ext_data = get_n_extreme_points(100, point_cloud_data);
            res2 = run_heuristic(point_cloud_model,ext_data,point_cloud_data);
            L2error = computeL2error(point_cloud_model,ext_data);
            DebugOn("L2 error with exterme set after = " << L2error << endl);
            L2error = computeL2error(point_cloud_model,point_cloud_data);
            DebugOn("L2 error with full set after = " << L2error << endl);
            res = res1;
            get<0>(res) += get<0>(res2);
            get<1>(res) += get<1>(res2);
            get<2>(res) += get<2>(res2);
            get<3>(res) += get<3>(res2);
            get<4>(res) += get<4>(res2);
            get<5>(res) += get<5>(res2);
            auto roll = get<0>(res);auto pitch = get<1>(res);auto yaw = get<2>(res);auto x_shift = get<3>(res);auto y_shift = get<4>(res);auto z_shift = get<5>(res);
            DebugOn("Final Roll = " << roll << endl);
            DebugOn("Final Pitch = " << pitch << endl);
            DebugOn("Final Yaw = " << yaw << endl);
            DebugOn("Final tx = " << x_shift << endl);
            DebugOn("Final ty = " << y_shift << endl);
            DebugOn("Final tz = " << z_shift << endl);
//            plot(point_cloud_model,old_point_cloud);
            apply_rot_trans(roll, pitch, yaw, x_shift, y_shift, z_shift, old_point_cloud);
            plot(point_cloud_model,old_point_cloud);
        }
        
        bool show_after = false;
        if (show_after) {
            plot(point_cloud_model,ext_data,0.5);
            plot(point_cloud_model,point_cloud_data);
        }
        bool save_file = false;
        if(save_file){
            DebugOn("Saving new las file");
            //    LASreadOpener lasreadopener_final;
            //    lasreadopener_final.set_file_name(LiDAR_file1.c_str());
            //    lasreadopener_final.set_populate_header(TRUE);
            //    LASreader* lasreader = lasreadopener_final.open();
            LASheader lasheader;
            lasheader.global_encoding = 1;
            lasheader.x_scale_factor = 0.01;
            lasheader.y_scale_factor = 0.01;
            lasheader.z_scale_factor = 0.01;
            lasheader.x_offset =  500000.0;
            lasheader.y_offset = 4100000.0;
            lasheader.z_offset = 0.0;
            lasheader.point_data_format = 1;
            lasheader.point_data_record_length = 28;
            
            LASwriteOpener laswriteopener;
            auto name = Model_file.substr(0,Model_file.find('.'));
            auto fname = name+".laz";
            laswriteopener.set_file_name(fname.c_str());
            LASwriter* laswriter = laswriteopener.open(&lasheader);
            LASpoint laspoint;
            laspoint.init(&lasheader, lasheader.point_data_format, lasheader.point_data_record_length, 0);
            for (auto i = 0; i< x_vec0.size(); i++) {
                laspoint.set_x(x_vec0[i]);
                laspoint.set_y(y_vec0[i]);
                laspoint.set_z(z_vec0[i]);
                laswriter->write_point(&laspoint);
                laswriter->update_inventory(&laspoint);
            }
            for (auto i = 0; i< x_vec1.size(); i++) {
                laspoint.set_x(x_vec1[i]);
                laspoint.set_y(y_vec1[i]);
                laspoint.set_z(z_vec1[i]);
                laswriter->write_point(&laspoint);
                laswriter->update_inventory(&laspoint);
            }
            laswriter->update_header(&lasheader, TRUE);
            laswriter->close();
            delete laswriter;
        }
        
        return 0;
        
    }
    
    bool show_umbrella = false;
    if(show_umbrella){
        vector<double> x_vec0, y_vec0, z_vec0, x_vec, y_vec, z_vec;
        map<string, tuple<int,int,int>> mapping;
        int x_min_idx, x_max_idx, y_min_idx, y_max_idx, z_min_idx, z_max_idx;
        double x1 = 10;
        double y1 = 10;
        double z1 = 10;
        x_vec0.push_back(x1);
        y_vec0.push_back(y1);
        z_vec0.push_back(z1);
        for (int a = 0; a < 10; a++){
            x_vec0.push_back(a);
            y_vec0.push_back(a);
            z_vec0.push_back(a);
        }
        double x_rot1, y_rot1, z_rot1, x_min = numeric_limits<double>::max(), x_max = numeric_limits<double>::lowest(), y_min = numeric_limits<double>::max(), y_max = numeric_limits<double>::lowest(), z_min = numeric_limits<double>::max(), z_max = numeric_limits<double>::lowest();
        double angles[] = {0, 0.5, -0.5};
        //      list<double> mylist (angles,angles+3);
        for (int a = 0; a < 3; a++){
            for (int b = 0; b <3; b++){
                for (int c = 0; c <3; c++){
                    
                    x_rot1 = x1*cos(angles[a])*cos(angles[b]) + y1*(cos(angles[a])*sin(angles[b])*sin(angles[c]) - sin(angles[a])*cos(angles[c])) + z1*(cos(angles[a])*sin(angles[b])*cos(angles[c]) + sin(angles[a])*sin(angles[c]));
                    
                    y_rot1 = x1*sin(angles[a])*cos(angles[b]) + y1*(sin(angles[a])*sin(angles[b])*sin(angles[c]) + cos(angles[a])*cos(angles[c])) + z1*(sin(angles[a])*sin(angles[b])*cos(angles[c]) - cos(angles[a])*sin(angles[c]));
                    
                    z_rot1 = x1*sin(-1*angles[b]) + y1*(cos(angles[b])*sin(angles[c])) + z1*(cos(angles[b])*cos(angles[c]));
                    
                    x_vec.push_back(x_rot1);
                    y_vec.push_back(y_rot1);
                    z_vec.push_back(z_rot1);
                    if(x_min>x_rot1){
                        x_min = x_rot1;
                        x_min_idx = x_vec.size()-1;
                        mapping["x_min"] = {a,b,c};
                    }
                    if(y_min>y_rot1){
                        y_min = y_rot1;
                        y_min_idx = x_vec.size()-1;
                        mapping["y_min"] = {a,b,c};
                    }
                    if(z_min>z_rot1){
                        z_min = z_rot1;
                        z_min_idx = x_vec.size()-1;
                        mapping["z_min"] = {a,b,c};
                    }
                    if(x_max<x_rot1){
                        x_max = x_rot1;
                        x_max_idx = x_vec.size()-1;
                        mapping["x_max"] = {a,b,c};
                    }
                    if(y_max<y_rot1){
                        y_max = y_rot1;
                        y_max_idx = x_vec.size()-1;
                        mapping["y_max"] = {a,b,c};
                    }
                    if(z_max<z_rot1){
                        z_max = z_rot1;
                        z_max_idx = x_vec.size()-1;
                        mapping["z_max"] = {a,b,c};
                    }
                }
            }
        }
        x_vec0.push_back(x_vec.at(x_min_idx));
        y_vec0.push_back(y_vec.at(x_min_idx));
        z_vec0.push_back(z_vec.at(x_min_idx));
        x_vec0.push_back(x_vec.at(x_max_idx));
        y_vec0.push_back(y_vec.at(x_max_idx));
        z_vec0.push_back(z_vec.at(x_max_idx));
        x_vec0.push_back(x_vec.at(y_min_idx));
        y_vec0.push_back(y_vec.at(y_min_idx));
        z_vec0.push_back(z_vec.at(y_min_idx));
        x_vec0.push_back(x_vec.at(y_max_idx));
        y_vec0.push_back(y_vec.at(y_max_idx));
        z_vec0.push_back(z_vec.at(y_max_idx));
        x_vec0.push_back(x_vec.at(z_min_idx));
        y_vec0.push_back(y_vec.at(z_min_idx));
        z_vec0.push_back(z_vec.at(z_min_idx));
        x_vec0.push_back(x_vec.at(z_max_idx));
        y_vec0.push_back(y_vec.at(z_max_idx));
        z_vec0.push_back(z_vec.at(z_max_idx));
        DebugOn("Extreme points mapping:" << endl);
        for (const auto &triplets: mapping) {
            DebugOn("(" << angles[get<0>(triplets.second)] << "," << angles[get<1>(triplets.second)] << "," << angles[get<2>(triplets.second)] << ")" << endl);
        }
        
        
        angles[1] = 0.1;
        angles[2] = -0.1;
        //      list<double> mylist (angles,angles+3);
        for (int a = 0; a < 3; a++){
            for (int b = 0; b <3; b++){
                for (int c = 0; c <3; c++){
                    
                    x_rot1 = x1*cos(angles[a])*cos(angles[b]) + y1*(cos(angles[a])*sin(angles[b])*sin(angles[c]) - sin(angles[a])*cos(angles[c])) + z1*(cos(angles[a])*sin(angles[b])*cos(angles[c]) + sin(angles[a])*sin(angles[c]));
                    
                    y_rot1 = x1*sin(angles[a])*cos(angles[b]) + y1*(sin(angles[a])*sin(angles[b])*sin(angles[c]) + cos(angles[a])*cos(angles[c])) + z1*(sin(angles[a])*sin(angles[b])*cos(angles[c]) - cos(angles[a])*sin(angles[c]));
                    
                    z_rot1 = x1*sin(-1*angles[b]) + y1*(cos(angles[b])*sin(angles[c])) + z1*(cos(angles[b])*cos(angles[c]));
                    
                    x_vec.push_back(x_rot1);
                    y_vec.push_back(y_rot1);
                    z_vec.push_back(z_rot1);
                }}}
        
        angles[1] = 0.25;
        angles[2] = -0.25;
        //      list<double> mylist (angles,angles+3);
        for (int a = 0; a < 3; a++){
            for (int b = 0; b <3; b++){
                for (int c = 0; c <3; c++){
                    
                    x_rot1 = x1*cos(angles[a])*cos(angles[b]) + y1*(cos(angles[a])*sin(angles[b])*sin(angles[c]) - sin(angles[a])*cos(angles[c])) + z1*(cos(angles[a])*sin(angles[b])*cos(angles[c]) + sin(angles[a])*sin(angles[c]));
                    
                    y_rot1 = x1*sin(angles[a])*cos(angles[b]) + y1*(sin(angles[a])*sin(angles[b])*sin(angles[c]) + cos(angles[a])*cos(angles[c])) + z1*(sin(angles[a])*sin(angles[b])*cos(angles[c]) - cos(angles[a])*sin(angles[c]));
                    
                    z_rot1 = x1*sin(-1*angles[b]) + y1*(cos(angles[b])*sin(angles[c])) + z1*(cos(angles[b])*cos(angles[c]));
                    
                    x_vec.push_back(x_rot1);
                    y_vec.push_back(y_rot1);
                    z_vec.push_back(z_rot1);
                }}}
        
        
        
        
        namespace plt = matplotlibcpp;
        
        std::map<std::string, std::string> keywords, keywords2;
        keywords["marker"] = "s";
        keywords["linestyle"] = "None";
        keywords["ms"] = "2";
        plt::plot3(x_vec0, y_vec0, z_vec0, x_vec, y_vec, z_vec, keywords);
        
        //    keywords2["marker"] = "s";
        //    keywords2["ms"] = "0.1";
        
        plt::show();
        return 0;
    }
    
    
    
    int output = 0;
    double tol = 1e-6;
    bool optimize = true;
    double solver_time_end, total_time_end, solve_time, total_time;
    int grid_step = 1;
    map<pair<double,double>,tuple<int,double,double,double,UAVPoint*>> grid1, grid2, grid; /* Grid with all cells where <x,y> is the key and <nb_measurements, min_z, max_z, av_z> is the cell data */
    int ground_z = 0;
    
    bool one_csv_file = true;
    
    if(one_csv_file){
        string CSV_file = string(prj_dir)+"/data_sets/LiDAR/CSV_Export.csv";
        if(argc>1){
            CSV_file = argv[1];
        }
        rapidcsv::Document  CSV_data(CSV_file, rapidcsv::LabelParams(0, -1));
        int nb_rows = CSV_data.GetRowCount();
        if(nb_rows<2){
            throw invalid_argument("csv file with less than 2 points");
            return 0;
        }
        DebugOn("csv file has " << nb_rows << " rows" << endl);
        DebugOn("csv file has " << CSV_data.GetColumnCount() << " columns" << endl);
        /* Values below are used to identify u-turns in drone flight */
        bool neg_x = false;/* x is decreasing */
        bool neg_y = false;/* y is decreasing */
        
        vector<UAVPoint*> UAVPoints;
        vector<LidarPoint*> LidarPoints;
        map<int,shared_ptr<Frame>> frames;
        map<int,shared_ptr<Frame>> frames1, frames2;
        vector<double> uav_x, uav_y, uav_z;
        vector<double> uav_x1, uav_y1, uav_z1;
        vector<double> old_x_vec1,old_y_vec1,old_z_vec1;
        vector<double> old_x_vec2,old_y_vec2,old_z_vec2;
        vector<double> x_vec1,y_vec1,z_vec1,zmin_vec1,zmax_vec1;
        vector<double> x_vec2,y_vec2,z_vec2,zmin_vec2,zmax_vec2;
        vector<double> x_vec1_r,y_vec1_r,z_vec1_r;
        vector<double> x_vec2_r,y_vec2_r,z_vec2_r;
        vector<double> x_shift1_r,y_shift1_r,z_shift1_r;
        vector<double> x_shift2_r,y_shift2_r,z_shift2_r;
        vector<double> x_shift,y_shift,z_shift;
        vector<double> x_shift1,y_shift1,z_shift1;
        vector<double> x_shift2,y_shift2,z_shift2;
        vector<double> x_shift_all1,y_shift_all1,z_shift_all1;
        vector<double> x_shift_all2,y_shift_all2,z_shift_all2;
        vector<double> uav_roll1,uav_pitch1,uav_yaw1;
        vector<double> uav_roll2,uav_pitch2,uav_yaw2;
        vector<double> x_combined,y_combined,z_combined,zmin_combined,zmax_combined;
        set<double> timestamps;
        size_t uav_id = 0;
        bool new_uav = true, u_turn = false, frame1 = true;
        double unix_time, delta_x = 0, delta_y = 0;
        double scale = 1e-2;
        pair<map<int,shared_ptr<Frame>>::iterator,bool> frame_ptr;
        for (int i = 0; i< nb_rows-1; i++) { // Input iterator
            auto laser_id = CSV_data.GetCell<int>("LaserID", i);
//            if(laser_id!=15){/* Only keep points from Nadir laser */
//                continue;
//            }
            
            
            auto frame_id = CSV_data.GetCell<int>(0, i);
            unix_time = CSV_data.GetCell<double>("Time", i);
            new_uav = (uav_id==0) || (UAVPoints[uav_id-1]->_frame_id != frame_id);
            if(new_uav){
                auto uav_x1 = CSV_data.GetCell<double>("Track_UTM_E", i)*scale;
                auto uav_y1 = CSV_data.GetCell<double>("Track_UTM_N", i)*scale;
                if(UAVPoints.size()==2){
                    auto uav_x0 = UAVPoints.back()->_x;
                    auto uav_y0 = UAVPoints.back()->_y;
                    neg_x = (uav_x1 - uav_x0) < 0;/* x is decreasing */
                    neg_y = (uav_y1 - uav_y0) < 0;/* y is decreasing */
                }
                else if(UAVPoints.size()>2){
                    auto uav_x0 = UAVPoints.back()->_x;
                    auto uav_y0 = UAVPoints.back()->_y;
                    bool neg_x_new = (uav_x1 - uav_x0) < 0;/* x is decreasing */
                    bool neg_y_new = (uav_y1 - uav_y0) < 0;/* y is decreasing */
                    if(neg_x_new!=neg_x || neg_y_new!=neg_y){/* A U turn is being detected */
                        u_turn = true;
                        frame1 = false;
                        neg_x = neg_x_new;
                        neg_y = neg_y_new;
                    }
                    else {
                        u_turn = false;
                    }
                }
                UAVPoints.push_back(new UAVPoint());
                UAVPoints[uav_id]->_frame_id = frame_id;
                UAVPoints[uav_id]->_x = uav_x1;
                UAVPoints[uav_id]->_y = uav_y1;
                UAVPoints[uav_id]->_height = CSV_data.GetCell<double>("Track_UTM_Height", i)*scale;
                UAVPoints[uav_id]->set_unix_time(unix_time);
                if(uav_id>0){
                    UAVPoints[uav_id]->_prev = UAVPoints[uav_id-1];
                    UAVPoints[uav_id-1]->_next = UAVPoints[uav_id];
                }
                uav_x.push_back(UAVPoints[uav_id]->_x);
                uav_y.push_back(UAVPoints[uav_id]->_y);
                uav_z.push_back(UAVPoints[uav_id]->_height);
                frame_ptr = frames.insert(make_pair(UAVPoints[uav_id]->_frame_id, make_shared<Frame>(UAVPoints[uav_id]->_frame_id, UAVPoints[uav_id]->_unix_time)));
                frame_ptr.first->second->add_UAV_point(UAVPoints[uav_id]);
                if(frame1){/* Has not performed a u-turn yet, keep adding to frames1 */
                    frames1.insert(make_pair(frame_ptr.first->second->_id, frame_ptr.first->second));
                    DebugOff("Frame " << frame_ptr.first->first << " in flight line 1" << endl);
                }
                else{/* Already turned, keep adding to frames2 */
                    frames2.insert(make_pair(frame_ptr.first->second->_id, frame_ptr.first->second));
                    DebugOff("Frame " << frame_ptr.first->first << " in flight line 2" << endl);
                }
                if(u_turn){
                    DebugOn("Detected a Uturn at frame " << frame_ptr.first->first << endl);
                }
                uav_id++;
            }
            
            auto xpos = CSV_data.GetCell<double>("UTM_E", i)*scale;
            auto ypos = CSV_data.GetCell<double>("UTM_N", i)*scale;
            auto zpos = CSV_data.GetCell<double>("UTM_Height", i)*scale;
            LidarPoints.push_back(new LidarPoint(laser_id,unix_time,xpos,ypos,zpos));
            frame_ptr.first->second->add_lidar_point(LidarPoints.back());
            LidarPoints.back()->_uav_pt = frame_ptr.first->second->_uav_point;
            
            //                uav_x1.push_back((frame_ptr.first->second._uav_points.front()->_longitude+582690.8242)*1e-5);
            //                uav_y1.push_back((frame_ptr.first->second._uav_points.front()->_latitude+4107963.58)*1e-5);
            //                uav_z1.push_back(frame_ptr.first->second._uav_points.front()->_height*100);
        }
        DebugOn("Read " << uav_id << " frames" << endl);
        DebugOn(frames1.size() << " frames in flight line 1" << endl);
        DebugOn(frames2.size() << " frames in flight line 2" << endl);
        DebugOn(LidarPoints.size() << " lidar points read" << endl);
        
        vector<vector<double>> point_cloud1, point_cloud2, uav1, uav2;
        
        int nb_pts_per_frame1 = 0, nb_pts_per_frame2 = 0, idx1 = 1, idx2 = 1;
        double coef = 0;
        for (const auto &frame: frames1) {
            if (frame.second->_id == frames1.rbegin()->first) {/* Ignore last frame */
                continue;
            }
            nb_pts_per_frame1 += frame.second->_lidar_points->size();
            int i = 0;
            for (const auto &p: *frame.second->_lidar_points) {
//                coef = get_interpolation_coef(p->_unix_time, frame.second->_uav_point, frame.second->_uav_point->_next);
//                if(i%10==0 && p->_laser_id==15){
                if(i%100==0 && p->_laser_id>6 && p->_laser_id<23){
//                if(p->_z>12.6){
                    x_vec1.push_back(p->_x);
                    x_shift1.push_back(frame.second->_uav_point->_x+coef*(frame.second->_uav_point->_next->_x - frame.second->_uav_point->_x));
                    y_vec1.push_back(p->_y);
                    y_shift1.push_back(frame.second->_uav_point->_y+coef*(frame.second->_uav_point->_next->_y - frame.second->_uav_point->_y));
                    z_vec1.push_back(p->_z);
                    z_shift1.push_back(frame.second->_uav_point->_height+coef*(frame.second->_uav_point->_next->_height - frame.second->_uav_point->_height));
                    point_cloud1.push_back({p->_x, p->_y, p->_z});
                    uav1.push_back({x_shift1.back(), y_shift1.back(), z_shift1.back()});
                }
                old_x_vec1.push_back(p->_x*1/scale);
                x_shift_all1.push_back((frame.second->_uav_point->_x+coef*(frame.second->_uav_point->_next->_x - frame.second->_uav_point->_x))*1/scale);
                old_y_vec1.push_back(p->_y*1/scale);
                y_shift_all1.push_back((frame.second->_uav_point->_y+coef*(frame.second->_uav_point->_next->_y - frame.second->_uav_point->_y))*1/scale);
                old_z_vec1.push_back(p->_z*1/scale);
                z_shift_all1.push_back((frame.second->_uav_point->_height+coef*(frame.second->_uav_point->_next->_height - frame.second->_uav_point->_height))*1/scale);
                i++;
            }
            
        }
        DebugOff("last frame id = " << frames2.rbegin()->first << endl);
        for (const auto &frame: frames2) {
            if (frame.second->_id == frames2.rbegin()->first) {/* Ignore last frame */
                continue;
            }
            nb_pts_per_frame2 += frame.second->_lidar_points->size();
            int i = 0;
            for (auto const &p: *frame.second->_lidar_points) {
//                coef = get_interpolation_coef(p->_unix_time, frame.second->_uav_point, frame.second->_uav_point->_next);
//                if(i%10==0 && p->_laser_id==15){
                if(i%100==0 && p->_laser_id>6 && p->_laser_id<23){
//                if(p->_z>12.6){
                    x_vec2.push_back(p->_x);
                    x_shift2.push_back(frame.second->_uav_point->_x+coef*(frame.second->_uav_point->_next->_x - frame.second->_uav_point->_x));
                    y_vec2.push_back(p->_y);
                    y_shift2.push_back(frame.second->_uav_point->_y+coef*(frame.second->_uav_point->_next->_y - frame.second->_uav_point->_y));
                    z_vec2.push_back(p->_z);
                    z_shift2.push_back(frame.second->_uav_point->_height+coef*(frame.second->_uav_point->_next->_height - frame.second->_uav_point->_height));
                    point_cloud2.push_back({p->_x, p->_y, p->_z});
                    uav2.push_back({x_shift2.back(), y_shift2.back(), z_shift2.back()});
                }
                old_x_vec2.push_back(p->_x*1/scale);
                x_shift_all2.push_back((frame.second->_uav_point->_x+coef*(frame.second->_uav_point->_next->_x - frame.second->_uav_point->_x))*1/scale);
                old_y_vec2.push_back(p->_y*1/scale);
                y_shift_all2.push_back((frame.second->_uav_point->_y+coef*(frame.second->_uav_point->_next->_y - frame.second->_uav_point->_y))*1/scale);
                old_z_vec2.push_back(p->_z*1/scale);
                z_shift_all2.push_back((frame.second->_uav_point->_height+coef*(frame.second->_uav_point->_next->_height - frame.second->_uav_point->_height))*1/scale);
                i++;
            }
        }
        if(frames1.size()!=0)
            DebugOn("Average number of points per frame in flight line 1 = " << nb_pts_per_frame1/frames1.size() << endl);
        if(frames2.size()!=0)
            DebugOn("Average number of points per frame in flight line 2 = " << nb_pts_per_frame2/frames2.size() << endl);
        DebugOn("Number of points in flight line 1 = " << x_vec1.size() << endl);
        DebugOn("Number of points in flight line 2 = " << x_vec2.size() << endl);
        
        
        bool run_goICP = false;
        if(run_goICP){/* Run GoICP inline */
            using namespace Go_ICP;
            int Nm = old_x_vec1.size(), Nd = x_vec2.size(), NdDownsampled = 0;
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
                pModel[i].x  = old_x_vec1[i];
                avg_x += pModel[i].x;
                pModel[i].y  = old_y_vec1[i];
                avg_y += pModel[i].y;
                pModel[i].z = old_z_vec1[i];
                avg_z += pModel[i].z;
            }
            avg_x /= Nm;avg_y /= Nm;avg_z /= Nm;
            centralize(Nm, &pModel, avg_x, avg_y, avg_z);
            avg_x = 0;avg_y = 0;avg_z = 0;
            pData = (POINT3D *)malloc(sizeof(POINT3D) * Nd);
            for(int i = 0; i < Nd; i++)
            {
                pData[i].x  = x_vec2[i];
                avg_x += pData[i].x;
                pData[i].y  = y_vec2[i];
                avg_y += pData[i].y;
                pData[i].z = z_vec2[i];
                avg_z += pData[i].z;
            }
            avg_x /= Nd;avg_y /= Nd;avg_z /= Nd;
            centralize(Nd, &pData, avg_x, avg_y, avg_z);
            for(int i = 0; i < Nd; i++)
            {
                if(max_x<pData[i].x){
                    max_x = pData[i].x;
                }
                if(min_x>pData[i].x){
                    min_x = pData[i].x;
                }
                if(max_y<pData[i].y){
                    max_y = pData[i].y;
                }
                if(min_y>pData[i].y){
                    min_y = pData[i].y;
                }
                if(max_z<pData[i].z){
                    max_z = pData[i].z;
                }
                if(min_z>pData[i].z){
                    min_z = pData[i].z;
                }
            }
            for(int i = 0; i < Nm; i++)
            {
                if(max_x<pModel[i].x){
                    max_x = pModel[i].x;
                }
                if(min_x>pModel[i].x){
                    min_x = pModel[i].x;
                }
                if(max_y<pModel[i].y){
                    max_y = pModel[i].y;
                }
                if(min_y>pModel[i].y){
                    min_y = pModel[i].y;
                }
                if(max_z<pModel[i].z){
                    max_z = pModel[i].z;
                }
                if(min_z>pModel[i].z){
                    min_z = pModel[i].z;
                }
            }
            scale_all(Nm, &pModel, max_x, max_y, max_z, min_x, min_y, min_z);
            scale_all(Nd, &pData, max_x, max_y, max_z, min_x, min_y, min_z);
            bool plot_GoICP = true;
            if (plot_GoICP) {
                vector<double> x_model, y_model, z_model, x_data, y_data, z_data;
                x_model.resize(Nm); y_model.resize(Nm); z_model.resize(Nm);
                x_data.resize(Nd); y_data.resize(Nd); z_data.resize(Nd);
                for(int i = 0; i < Nm; i++)
                {
                    x_model[i] = pModel[i].x;
                    y_model[i] = pModel[i].y;
                    z_model[i] = pModel[i].z;
                }
                for(int i = 0; i < Nd; i++)
                {
                    x_data[i] = pData[i].x;
                    y_data[i] = pData[i].y;
                    z_data[i] = pData[i].z;
                }
                namespace plt = matplotlibcpp;
                
                std::map<std::string, std::string> keywords, keywords2;
                keywords["marker"] = "s";
                keywords["linestyle"] = "None";
                keywords["ms"] = "0.05";
                //    plt::plot3(x_combined, y_combined, zmax_combined, keywords);
                plt::plot3(x_model, y_model, z_model, x_data, y_data, z_data, keywords);
//                plt::plot3(x_data, y_data, z_data, keywords);
                plt::show();
            }
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
            std::ofstream outputFileFinal("./optimal_pointcloud2.txt");
            bool plot_unscaled = false;
            if (true||plot_GoICP) {
                for(int i = 0; i < goicp.Nd; i++)
                {
                    POINT3D& p = goicp.pData[i];
                    goicp.pData[i].x = goicp.optR.val[0][0]*p.x + goicp.optR.val[0][1]*p.y + goicp.optR.val[0][2]*p.z + goicp.optT.val[0][0];
                    goicp.pData[i].y = goicp.optR.val[1][0]*p.x + goicp.optR.val[1][1]*p.y + goicp.optR.val[1][2]*p.z + goicp.optT.val[1][0];
                    goicp.pData[i].z = goicp.optR.val[2][0]*p.x + goicp.optR.val[2][1]*p.y + goicp.optR.val[2][2]*p.z + goicp.optT.val[2][0];
                }
                vector<double> x_model, y_model, z_model, x_data, y_data, z_data;
                x_model.resize(Nm); y_model.resize(Nm); z_model.resize(Nm);
                x_data.resize(Nd); y_data.resize(Nd); z_data.resize(Nd);
                for(int i = 0; i < Nm; i++)
                {
                    if(plot_unscaled){
                        x_model[i] = old_x_vec1[i];
                        y_model[i] = old_y_vec1[i];
                        z_model[i] = old_z_vec1[i];

                    }
                    else{
                        x_model[i] = pModel[i].x;
                        y_model[i] = pModel[i].y;
                        z_model[i] = pModel[i].z;
                    }
                }
                for(int i = 0; i < Nd; i++)
                {
                    if(plot_unscaled){
                        x_data[i] = 0.5*(pData[i].x + 1) * (max_x - min_x) + min_x;
                        y_data[i] = 0.5*(pData[i].y + 1) * (max_y - min_y) + min_y;
                        z_data[i] = 0.5*(pData[i].z + 1) * (max_z - min_z) + min_z;
                        x_data[i] += avg_x;
                        y_data[i] += avg_y;
                        z_data[i] += avg_z;
                    }
                    else {
                        x_data[i] = pData[i].x;
                        y_data[i] = pData[i].y;
                        z_data[i] = pData[i].z;
                    }
                }
                namespace plt = matplotlibcpp;
                
                std::map<std::string, std::string> keywords, keywords2;
                keywords["marker"] = "s";
                keywords["linestyle"] = "None";
                keywords["ms"] = "0.05";
                //    plt::plot3(x_combined, y_combined, zmax_combined, keywords);
                plt::plot3(x_model, y_model, z_model, x_data, y_data, z_data, keywords);
//                plt::plot3(x_data, y_data, z_data, keywords);
                plt::show();
            }
//            for(int i = 0; i < goicp.Nd; i++)
//            {
//                POINT3D& p = goicp.pData[i];
//                goicp.pData[i].x = goicp.optR.val[0][0]*p.x + goicp.optR.val[0][1]*p.y + goicp.optR.val[0][2]*p.z + goicp.optT.val[0][0];
//                goicp.pData[i].y = goicp.optR.val[1][0]*p.x + goicp.optR.val[1][1]*p.y + goicp.optR.val[1][2]*p.z + goicp.optT.val[1][0];
//                goicp.pData[i].z = goicp.optR.val[2][0]*p.x + goicp.optR.val[2][1]*p.y + goicp.optR.val[2][2]*p.z + goicp.optT.val[2][0];
////                outputFileFinal << goicp.pData[i].x << " " << goicp.pData[i].y << " " << goicp.pData[i].z << std::endl;
//            }
            /* Copy model points unchanged */
            int n1 = old_x_vec1.size(), n2 = old_x_vec2.size();
            
            auto tot_pts = n1+n2;
            x_combined.resize(tot_pts);
            y_combined.resize(tot_pts);
            z_combined.resize(tot_pts);
            
            
            for (auto i = 0; i< n1; i++) {
                x_combined[i] = old_x_vec1[i];
                y_combined[i] = old_y_vec1[i];
                z_combined[i] = old_z_vec1[i];
            }
//            for (auto i = 0; i< n2; i++) {
//                x_combined[n1+i] = old_x_vec2[i];
//                y_combined[n1+i] = old_y_vec2[i];
//                z_combined[n1+i] = old_z_vec2[i];
//            }
            pFullData = (POINT3D *)malloc(sizeof(POINT3D) * n2);
            for(int i = 0; i < n2; i++)
            {
                pFullData[i].x  = old_x_vec2[i];
                pFullData[i].y  = old_y_vec2[i];
                pFullData[i].z = old_z_vec2[i];
            }
            centralize(n2, &pFullData, avg_x, avg_y, avg_z);/* We're using the averages of x_vec2 */
            scale_all(n2, &pFullData, max_x, max_y, max_z, min_x, min_y, min_z);/* We're using the min/max values of x_vec2 and x_vec1 */
            vector<double> new_data_x, new_data_y, new_data_z;
            for(int i = 0; i < n2; i++)
            {
                double px = pFullData[i].x;double py = pFullData[i].y;double pz = pFullData[i].z;
//                POINT3D& p = pFullData[i];
                pFullData[i].x = goicp.optR.val[0][0]*px + goicp.optR.val[0][1]*py + goicp.optR.val[0][2]*pz + goicp.optT.val[0][0];
                pFullData[i].y = goicp.optR.val[1][0]*px + goicp.optR.val[1][1]*py + goicp.optR.val[1][2]*pz + goicp.optT.val[1][0];
                pFullData[i].z = goicp.optR.val[2][0]*px + goicp.optR.val[2][1]*py + goicp.optR.val[2][2]*pz + goicp.optT.val[2][0];
                /* Rescaling to original coordinates */
                x_combined[n1+i] = 0.5*(pFullData[i].x + 1) * (max_x - min_x) + min_x;
                y_combined[n1+i] = 0.5*(pFullData[i].y + 1) * (max_y - min_y) + min_y;
                z_combined[n1+i] = 0.5*(pFullData[i].z + 1) * (max_z - min_z) + min_z;
                x_combined[n1+i] += avg_x;
                y_combined[n1+i] += avg_y;
                z_combined[n1+i] += avg_z;
                if(i<10)
                    outputFileFinal << x_combined[n1+i] << " " <<  y_combined[n1+i] << " " << z_combined[n1+i] << std::endl;
            }
            bool save_file = true;
            if(save_file){
                DebugOn("Saving new las file");
                //    LASreadOpener lasreadopener_final;
                //    lasreadopener_final.set_file_name(LiDAR_file1.c_str());
                //    lasreadopener_final.set_populate_header(TRUE);
                //    LASreader* lasreader = lasreadopener_final.open();
                LASheader lasheader;
                lasheader.global_encoding = 1;
                lasheader.x_scale_factor = 0.01;
                lasheader.y_scale_factor = 0.01;
                lasheader.z_scale_factor = 0.01;
                lasheader.x_offset =  500000.0;
                lasheader.y_offset = 4100000.0;
                lasheader.z_offset = 0.0;
                lasheader.point_data_format = 1;
                lasheader.point_data_record_length = 28;
                
                LASwriteOpener laswriteopener;
                auto fname = "GoICPmod30NoTrans.laz";
                laswriteopener.set_file_name(fname);
                LASwriter* laswriter = laswriteopener.open(&lasheader);
                LASpoint laspoint;
                laspoint.init(&lasheader, lasheader.point_data_format, lasheader.point_data_record_length, 0);
                for (auto i = 0; i< x_combined.size(); i++) {
                    laspoint.set_x(x_combined[i]);
                    laspoint.set_y(y_combined[i]);
                    laspoint.set_z(z_combined[i]);
                    laswriter->write_point(&laspoint);
                    laswriter->update_inventory(&laspoint);
                }
                laswriter->update_header(&lasheader, TRUE);
                laswriter->close();
                delete laswriter;
                namespace plt = matplotlibcpp;
                
                std::map<std::string, std::string> keywords, keywords2;
                keywords["marker"] = "s";
                keywords["linestyle"] = "None";
                keywords["ms"] = "0.05";
                //    plt::plot3(x_combined, y_combined, zmax_combined, keywords);
                plt::plot3(x_combined, y_combined, z_combined, keywords);
                plt::show();
            }
            delete(pModel);
            delete(pData);
            return 0;
        }
        
        bool scale_data = false;
        if(scale_data){
            auto Nm = x_vec1.size();
            auto Nd = x_vec2.size();
            POINT3D * pModel, * pData, * pFullData;
            // Load model and data point clouds
            pModel = (POINT3D *)malloc(sizeof(POINT3D) * Nm);
            double avg_x = 0, avg_y = 0, avg_z = 0;
            double max_x = numeric_limits<double>::lowest(), max_y = numeric_limits<double>::lowest(), max_z = numeric_limits<double>::lowest();
            double min_x = numeric_limits<double>::max(), min_y = numeric_limits<double>::max(), min_z = numeric_limits<double>::max();
            for(int i = 0; i < Nm; i++)
            {
                pModel[i].x  = old_x_vec1[i];
                avg_x += pModel[i].x;
                pModel[i].y  = old_y_vec1[i];
                avg_y += pModel[i].y;
                pModel[i].z = old_z_vec1[i];
                avg_z += pModel[i].z;
            }
            avg_x /= Nm;avg_y /= Nm;avg_z /= Nm;
            centralize(Nm, &pModel, avg_x, avg_y, avg_z);
            avg_x = 0;avg_y = 0;avg_z = 0;
            pData = (POINT3D *)malloc(sizeof(POINT3D) * Nd);
            for(int i = 0; i < Nd; i++)
            {
                pData[i].x  = x_vec2[i];
                avg_x += pData[i].x;
                pData[i].y  = y_vec2[i];
                avg_y += pData[i].y;
                pData[i].z = z_vec2[i];
                avg_z += pData[i].z;
            }
            avg_x /= Nd;avg_y /= Nd;avg_z /= Nd;
            centralize(Nd, &pData, avg_x, avg_y, avg_z);
            unit_scale(Nm, &pModel, Nd, &pData);
            /* What about the drone points?*/
        }
        bool bypass = true;
        double roll1 = 0, pitch1 = 0, yaw1 = 0, roll2 = 0, pitch2 = 0, yaw2 = 0, roll3 = 0, pitch3 = 0, yaw3 = 0;
        if(!bypass){

            tuple<double,double,double> opt_sol = run_ARMO("z", x_vec1, y_vec1, z_vec1, x_vec2, y_vec2, z_vec2, x_shift1, y_shift1, z_shift1, x_shift2, y_shift2, z_shift2, point_cloud1, point_cloud2, uav1, uav2);
            
    //
            roll1 = get<0>(opt_sol);
            pitch1 = get<1>(opt_sol);
            yaw1 = get<2>(opt_sol);
            apply_rotation(roll1, pitch1, yaw1, x_vec1, y_vec1, z_vec1, x_vec2, y_vec2, z_vec2, x_shift1, y_shift1, z_shift1, x_shift2, y_shift2, z_shift2);
            apply_rotation(roll1, pitch1, yaw1, old_x_vec1, old_y_vec1, old_z_vec1, old_x_vec2, old_y_vec2, old_z_vec2, x_shift_all1, y_shift_all1, z_shift_all1, x_shift_all2, y_shift_all2, z_shift_all2);
    //        auto projected1 = project("z", point_cloud1);/* Project on (x,y) plane, i.e., set z to zero */
    //        vector<double> zeros(0,x_vec2.size());
            opt_sol = run_ARMO("full", x_vec1, y_vec1, z_vec1, x_vec2, y_vec2, z_vec2, x_shift1, y_shift1, z_shift1, x_shift2, y_shift2, z_shift2, point_cloud1, point_cloud2, uav1, uav2);
//            roll2 = get<0>(opt_sol);
//            pitch2 = get<1>(opt_sol);
//            yaw2 = get<2>(opt_sol);
//            apply_rotation(roll2, pitch2, yaw2, x_vec1, y_vec1, z_vec1, x_vec2, y_vec2, z_vec2, x_shift1, y_shift1, z_shift1, x_shift2, y_shift2, z_shift2);
//            apply_rotation(roll2, pitch2, yaw2, old_x_vec1, old_y_vec1, old_z_vec1, old_x_vec2, old_y_vec2, old_z_vec2, x_shift_all1, y_shift_all1, z_shift_all1, x_shift_all2, y_shift_all2, z_shift_all2);
//            opt_sol = run_ARMO("y", x_vec1, y_vec1, z_vec1, x_vec2, y_vec2, z_vec2, x_shift1, y_shift1, z_shift1, x_shift2, y_shift2, z_shift2, point_cloud1, point_cloud2, uav1, uav2);
            roll3 = get<0>(opt_sol);
            pitch3 = get<1>(opt_sol);
            yaw3 = get<2>(opt_sol);
            apply_rotation(roll3, pitch3, yaw3, x_vec1, y_vec1, z_vec1, x_vec2, y_vec2, z_vec2, x_shift1, y_shift1, z_shift1, x_shift2, y_shift2, z_shift2);
        }
        double shifted_x, shifted_y, shifted_z;
        auto tot_pts = old_x_vec1.size()+old_x_vec2.size();
        x_combined.resize(tot_pts);
        y_combined.resize(tot_pts);
        z_combined.resize(tot_pts);
        
        roll3 = -0.25;pitch3=0.5;yaw3=-0.7;/* Manual values*/
        double beta = roll3*pi/180;// roll in radians
        double gamma = pitch3*pi/180; // pitch in radians
        double alpha = yaw3*pi/180; // yaw in radians

        for (auto i = 0; i< old_x_vec1.size(); i++) {
            shifted_x = old_x_vec1[i] - x_shift_all1[i];
            shifted_y = old_y_vec1[i] - y_shift_all1[i];
            shifted_z = old_z_vec1[i] - z_shift_all1[i];
            x_combined[i] = shifted_x*cos(alpha)*cos(beta) + shifted_y*(cos(alpha)*sin(beta)*sin(gamma) - sin(alpha)*cos(gamma)) + shifted_z*(cos(alpha)*sin(beta)*cos(gamma) + sin(alpha)*sin(gamma));
            y_combined[i] = shifted_x*sin(alpha)*cos(beta) + shifted_y*(sin(alpha)*sin(beta)*sin(gamma) + cos(alpha)*cos(gamma)) + shifted_z*(sin(alpha)*sin(beta)*cos(gamma) - cos(alpha)*sin(gamma));
            z_combined[i] = shifted_x*(-sin(beta)) + shifted_y*(cos(beta)*sin(gamma)) + shifted_z*(cos(beta)*cos(gamma));
            x_combined[i] += x_shift_all1[i];
            y_combined[i] += y_shift_all1[i];
            z_combined[i] += z_shift_all1[i];
        }
        beta *= -1;
        alpha *= -1;
        for (auto i = 0; i< old_x_vec2.size(); i++) {
            shifted_x = old_x_vec2[i] - x_shift_all2[i];
            shifted_y = old_y_vec2[i] - y_shift_all2[i];
            shifted_z = old_z_vec2[i] - z_shift_all2[i];
            x_combined[old_x_vec1.size()+i] = shifted_x*cos(alpha)*cos(beta) + shifted_y*(cos(alpha)*sin(beta)*sin(gamma) - sin(alpha)*cos(gamma)) + shifted_z*(cos(alpha)*sin(beta)*cos(gamma) + sin(alpha)*sin(gamma));
            y_combined[old_x_vec1.size()+i] = shifted_x*sin(alpha)*cos(beta) + shifted_y*(sin(alpha)*sin(beta)*sin(gamma) + cos(alpha)*cos(gamma)) + shifted_z*(sin(alpha)*sin(beta)*cos(gamma) - cos(alpha)*sin(gamma));
            z_combined[old_x_vec1.size()+i] = shifted_x*(-sin(beta)) + shifted_y*(cos(beta)*sin(gamma)) + shifted_z*(cos(beta)*cos(gamma));
            x_combined[old_x_vec1.size()+i] += x_shift_all2[i];
            y_combined[old_x_vec1.size()+i] += y_shift_all2[i];
            z_combined[old_x_vec1.size()+i] += z_shift_all2[i];
        }
        apply_rotation(roll3, pitch3, yaw3, old_x_vec1, old_y_vec1, old_z_vec1, old_x_vec2, old_y_vec2, old_z_vec2, x_shift_all1, y_shift_all1, z_shift_all1, x_shift_all2, y_shift_all2,z_shift_all2);
        bool plot_data = true;
        if(plot_data){
            namespace plt = matplotlibcpp;
            
            std::map<std::string, std::string> keywords, keywords2;
            keywords["marker"] = "s";
            keywords["linestyle"] = "None";
            keywords["ms"] = "0.1";
            //    plt::plot3(x_combined, y_combined, zmax_combined, keywords);
            //            plt::plot3(x_vec1, y_vec1, z_vec1, x_vec2, y_vec2, z_vec2, keywords);
//            plt::plot3(x_combined, y_combined, z_combined, keywords);
//            plt::plot3(uav_x, uav_y, uav_z, x_combined, y_combined, z_combined, keywords);
            keywords2["marker"] = "s";
            keywords2["ms"] = "0.05";
            //    plt::plot3(uav_x1, uav_y1, uav_z1, uav_x, uav_y, uav_z, keywords2);
            //    plt::plot3(uav_x1, uav_y1, uav_z1, keywords);
                plt::plot3(x_vec1, y_vec1, z_vec1, x_vec2, y_vec2, z_vec2, keywords);
//            plt::plot3(old_x_vec1, old_y_vec1, old_z_vec1, old_x_vec2, old_y_vec2, old_z_vec2, keywords2);
            //    plt::colorbar();
            // Enable legend.
            //    plt::legend();
            plt::show();
        }
        for (auto p : LidarPoints)
        {
            delete p;
        }
        LidarPoints.clear();
        for (auto p : UAVPoints)
        {
            delete p;
        }
        UAVPoints.clear();
        
        bool save_file = true;
        if(save_file){
            DebugOn("Saving new las file\n");
            //    LASreadOpener lasreadopener_final;
            //    lasreadopener_final.set_file_name(LiDAR_file1.c_str());
            //    lasreadopener_final.set_populate_header(TRUE);
            //    LASreader* lasreader = lasreadopener_final.open();
            LASheader lasheader;
            lasheader.global_encoding = 1;
            lasheader.x_scale_factor = 0.01;
            lasheader.y_scale_factor = 0.01;
            lasheader.z_scale_factor = 0.01;
            lasheader.x_offset =  500000.0;
            lasheader.y_offset = 4100000.0;
            lasheader.z_offset = 0.0;
            lasheader.point_data_format = 1;
            lasheader.point_data_record_length = 28;
            
            LASwriteOpener laswriteopener;
            auto roll = roll1 + roll2 + roll3;
            auto pitch = pitch1 + pitch2 + pitch3;
            auto yaw = yaw1 + yaw2 + yaw3;
            DebugOn("Final Roll = " << roll << ", final Pitch = " << pitch << ", final Yaw = " << yaw << endl);
            auto name = CSV_file.substr(0,CSV_file.find('.'));
            auto fname = name+"_ARMO_RPY_"+to_string(roll)+"_"+to_string(pitch)+"_"+to_string(yaw)+".laz";
            laswriteopener.set_file_name(fname.c_str());
            LASwriter* laswriter = laswriteopener.open(&lasheader);
            LASpoint laspoint;
            laspoint.init(&lasheader, lasheader.point_data_format, lasheader.point_data_record_length, 0);
            for (auto i = 0; i< x_combined.size(); i++) {
                laspoint.set_x(x_combined[i]);
                laspoint.set_y(y_combined[i]);
                laspoint.set_z(z_combined[i]);
                laswriter->write_point(&laspoint);
                laswriter->update_inventory(&laspoint);
            }
            laswriter->update_header(&lasheader, TRUE);
            laswriter->close();
            delete laswriter;
        }
        
    }
    return 0;
}

/* Centralize point cloud around origin */
void centralize(int n, POINT3D **  p, double avg_x, double avg_y, double avg_z){
    for(int i = 0; i < n; i++)
    {
        (*p)[i].x  -= avg_x;
        (*p)[i].y  -= avg_y;
        (*p)[i].z -= avg_z;
    }
}

/* Scale point clouds to [-1,1] */
void unit_scale(int n1, POINT3D **  p1, int n2, POINT3D **  p2){
    double max_x = numeric_limits<double>::lowest(), max_y = numeric_limits<double>::lowest(), max_z = numeric_limits<double>::lowest();
    double min_x = numeric_limits<double>::max(), min_y = numeric_limits<double>::max(), min_z = numeric_limits<double>::max();
    for(int i = 0; i < n1; i++)
    {
        if(max_x<(*p1)[i].x){
            max_x = (*p1)[i].x;
        }
        if(min_x>(*p1)[i].x){
            min_x = (*p1)[i].x;
        }
        if(max_y<(*p1)[i].y){
            max_y = (*p1)[i].y;
        }
        if(min_y>(*p1)[i].y){
            min_y = (*p1)[i].y;
        }
        if(max_z<(*p1)[i].z){
            max_z = (*p1)[i].z;
        }
        if(min_z>(*p1)[i].z){
            min_z = (*p1)[i].z;
        }
    }
    for(int i = 0; i < n2; i++)
    {
        if(max_x<(*p2)[i].x){
            max_x = (*p2)[i].x;
        }
        if(min_x>(*p2)[i].x){
            min_x = (*p2)[i].x;
        }
        if(max_y<(*p2)[i].y){
            max_y = (*p2)[i].y;
        }
        if(min_y>(*p2)[i].y){
            min_y = (*p2)[i].y;
        }
        if(max_z<(*p2)[i].z){
            max_z = (*p2)[i].z;
        }
        if(min_z>(*p2)[i].z){
            min_z = (*p2)[i].z;
        }
    }
    scale_all(n1, p1, max_x, max_y, max_z, min_x, min_y, min_z);
    scale_all(n2, p2, max_x, max_y, max_z, min_x, min_y, min_z);
}

/* Scale point cloud using privded max values */
void scale_all(int n1, POINT3D **  p1, double max_x, double max_y, double max_z, double min_x, double min_y, double min_z){
    for(int i = 0; i < n1; i++)
    {
        (*p1)[i].x = 2*(((*p1)[i].x - min_x)/(max_x - min_x)) - 1;
        (*p1)[i].y = 2*(((*p1)[i].y - min_y)/(max_y - min_y)) - 1;
        (*p1)[i].z = 2*(((*p1)[i].z - min_z)/(max_z - min_z)) - 1;
    }
}
    
void set_GoICP_options(GoICP& goicp){
    goicp.MSEThresh = 0.001;
    goicp.initNodeRot.a = -3.1416;
    goicp.initNodeRot.b = -3.1416;
    goicp.initNodeRot.c = -3.1416;
    goicp.initNodeRot.w = 6.2832;
    goicp.initNodeTrans.x = -0.25;
    goicp.initNodeTrans.y = -0.25;
    goicp.initNodeTrans.z = -0.25;
    goicp.initNodeTrans.w = 1;
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

/* Computes the interpolation coefficient based on time elapsed  */
double get_interpolation_coef(const double& lidar_time, UAVPoint* p1, UAVPoint* p2) {
    double t1 = p1->_unix_time;
    double t2 = p2->_unix_time;
    double tot_t = t2 - t1;
    double curr_t = lidar_time - t1;
    if(tot_t<0 || curr_t<0){
        throw invalid_argument("something is wrong with time records");
    }
    double perc_time = curr_t/tot_t;
    return perc_time;
}
    /* Return true if two cubes intersect */
    bool intersect(const vector<pair<double,double>>& a, const vector<pair<double,double>>& b) {
        return (a[0].first <= b[0].second && a[0].second >= b[0].first) &&
        (a[1].first <= b[1].second && a[1].second >= b[1].first) &&
        (a[2].first <= b[2].second && a[2].second >= b[2].first);
    }
    
    
    
    /* Return the coordinates of the cube center*/
    tuple<double,double,double> get_center(const vector<pair<double,double>>& cube){
        double x = 0, y = 0, z = 0;
        for (int a = 0; a < 2; a++){
            for (int b = 0; b <2; b++){
                for (int c = 0; c <2; c++){
                    if(a==1)
                        x += cube.at(0).first;
                    else
                        x += cube.at(0).second;
                    if(b==1)
                        y += cube.at(1).first;
                    else
                        y += cube.at(1).second;
                    if(c==1)
                        z += cube.at(2).first;
                    else
                        z += cube.at(2).second;
                }
            }
        }
        tuple<double,double,double> res;
        get<0>(res) = x/8.;
        get<1>(res) = y/8.;
        get<2>(res) = z/8.;
        return res;
    }
    
    /* Return the min-max values for x, y and z  for all possible rotations of p with angle +- angle*/
    vector<pair<double,double>> get_min_max(double angle, const vector<double>& p, const vector<double>& ref){
        double x1 = p[0], y1 = p[1], z1 = p[2], shifted_x, shifted_y, shifted_z, alpha, beta, gamma;
        double x_ref = ref[0], y_ref = ref[1], z_ref = ref[2];
        double x_rot1, y_rot1, z_rot1, x_min = numeric_limits<double>::max(), x_max = numeric_limits<double>::lowest(), y_min = numeric_limits<double>::max(), y_max = numeric_limits<double>::lowest(), z_min = numeric_limits<double>::max(), z_max = numeric_limits<double>::lowest();
        double angles[] = {0, angle, -angle};
        vector<pair<double,double>> min_max;
        for (int a = 0; a < 3; a++){
            for (int b = 0; b <3; b++){
                for (int c = 0; c <3; c++){
                    shifted_x = x1 - x_ref;
                    shifted_y = y1 - y_ref;
                    shifted_z = z1 - z_ref;
                    alpha = angles[a];
                    beta = angles[b];
                    gamma = angles[c];
                    x_rot1 = shifted_x*cos(alpha)*cos(beta) + shifted_y*(cos(alpha)*sin(beta)*sin(gamma) - sin(alpha)*cos(gamma)) + shifted_z*(cos(alpha)*sin(beta)*cos(gamma) + sin(alpha)*sin(gamma));
                    y_rot1 = shifted_x*sin(alpha)*cos(beta) + shifted_y*(sin(alpha)*sin(beta)*sin(gamma) + cos(alpha)*cos(gamma)) + shifted_z*(sin(alpha)*sin(beta)*cos(gamma) - cos(alpha)*sin(gamma));
                    z_rot1 = shifted_x*(-sin(beta)) + shifted_y*(cos(beta)*sin(gamma)) + shifted_z*(cos(beta)*cos(gamma));
                    x_rot1 += x_ref;
                    y_rot1 += y_ref;
                    z_rot1 += z_ref;
                    
                    if(x_min>x_rot1){
                        x_min = x_rot1;
                    }
                    if(y_min>y_rot1){
                        y_min = y_rot1;
                    }
                    if(z_min>z_rot1){
                        z_min = z_rot1;
                    }
                    if(x_max<x_rot1){
                        x_max = x_rot1;
                    }
                    if(y_max<y_rot1){
                        y_max = y_rot1;
                    }
                    if(z_max<z_rot1){
                        z_max = z_rot1;
                    }
                }
            }
        }
        min_max.push_back({x_min,x_max});
        min_max.push_back({y_min,y_max});
        min_max.push_back({z_min,z_max});
        return min_max;
    }

/* Run the MINLP ARMO model for registration */
tuple<double,double,double,double,double,double> run_ARMO_MINLP(bool bypass, string axis, const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data){
    double angle_max = 3, shift_max = 0.25;
    double roll_1 = 0, yaw_1 = 0, pitch_1 = 0;
    int nb_pairs = 0, min_nb_pairs = numeric_limits<int>::max(), max_nb_pairs = 0, av_nb_pairs = 0;
    size_t nm = point_cloud_model.size(), nd = point_cloud_data.size();
    vector<pair<double,double>> min_max_data;
    vector<vector<pair<double,double>>> min_max_model(nm);
    vector<int> nb_neighbors(nd);
    vector<vector<int>> neighbors(nd);
    vector<double> zeros = {0,0,0};
    
    param<> x1("x1"), x2("x2"), y1("y1"), y2("y2"), z1("z1"), z2("z2");
    param<> x_uav1("x_uav1"), y_uav1("y_uav1"), z_uav1("z_uav1");
    param<> x_uav2("x_uav2"), y_uav2("y_uav2"), z_uav2("z_uav2");
    //        return 0;
    int m = av_nb_pairs;
    //            int m = 1;
    vector<double> min_dist(nd,numeric_limits<double>::max());
    vector<int> nearest(nd);
    vector<string> nearest_id(nd);
    string i_str, j_str;
    indices Pairs("Pairs"), cells("cells");
    map<int,int> n2_map;
    int idx1 = 0;
    int idx2 = 0;
    int nb_max_neigh = 1;
    double dist_sq = 0;
    /* Compute nearest points in data point cloud */
    for (auto i = 0; i<nd; i++) {
        i_str = to_string(i+1);
        x1.add_val(i_str,point_cloud_data.at(i).at(0));
        y1.add_val(i_str,point_cloud_data.at(i).at(1));
        z1.add_val(i_str,point_cloud_data.at(i).at(2));
    }
    for (auto j = 0; j<nm; j++) {
        j_str = to_string(j+1);
        x2.add_val(j_str,point_cloud_model.at(j).at(0));
        y2.add_val(j_str,point_cloud_model.at(j).at(1));
        z2.add_val(j_str,point_cloud_model.at(j).at(2));
    }
    for (auto i = 0; i< nd; i++) {
        double min_dist = numeric_limits<double>::max();
        for (auto j = 0; j< nm; j++) {
            dist_sq = std::pow(point_cloud_data.at(i).at(0) - point_cloud_model.at(j).at(0),2) + std::pow(point_cloud_data.at(i).at(1) - point_cloud_model.at(j).at(1),2) + std::pow(point_cloud_data.at(i).at(2) - point_cloud_model.at(j).at(2),2);
            if(min_dist>dist_sq){
                min_dist = dist_sq;
                j_str = to_string(j+1);
                nearest_id[i] = j_str;
            }
        }
    }
    idx1 = 0;
    indices N1("N1"),N2("N2");
    DebugOn("nd = " << nd << endl);
    DebugOn("nm = " << nm << endl);
    
    N1 = range(1,nd);
    N2 = range(1,nm);
    cells = indices(N1,N2);
    //    cells.print();
    //    for (auto i = 0; i<nd; i++) {
    //        i_str = to_string(i+1);
    //        for (auto j = 0; j<nm; j++) {
    //            j_str = to_string(j+1);
    //            cells.add(i_str+","+j_str);
    //        }
    //    }
    Model<> Reg("Reg");
    var<> new_x1("new_x1"), new_y1("new_y1"), new_z1("new_z1");
    var<> yaw("yaw", -angle_max, angle_max), pitch("pitch", -angle_max, angle_max), roll("roll", -angle_max, angle_max);
    var<> x_shift("x_shift", -shift_max, shift_max), y_shift("y_shift", -shift_max, shift_max), z_shift("z_shift", -shift_max, shift_max);
    //                    var<> yaw("yaw", -1e-6, 1e-6), pitch("pitch",-1e-6, 1e-6), roll("roll", -1e-6, 1e-6);
    //                    var<> x_shift("x_shift", 0, 0), y_shift("y_shift", 0, 0), z_shift("z_shift", 0, 0);
    var<> delta("delta", pos_);
    var<> delta_min("delta_min", pos_);
    var<int> bin("bin",0,1);
    
    Reg.add(bin.in(cells));
    Reg.add(delta.in(cells), delta_min.in(N1));
    Reg.add(yaw.in(R(1)),pitch.in(R(1)),roll.in(R(1)));
    Reg.add(x_shift.in(R(1)),y_shift.in(R(1)),z_shift.in(R(1)));
    Reg.add(new_x1.in(N1), new_y1.in(N1), new_z1.in(N1));
    //        Reg.add(x_diff.in(cells), y_diff.in(cells), z_diff.in(cells));
    //                Reg.add(z_diff.in(cells));
    DebugOn("There are " << cells.size() << " cells" << endl);
//    for (int i = 0; i<6; i++) {
//        bin(to_string(i+1)+","+to_string(i+1))=1;
//    }
    
    
    Constraint<> DeltaMin("DeltaMin");
    DeltaMin = delta_min;
    DeltaMin -= bin.in_matrix(1, 1)*delta.in_matrix(1, 1);
    Reg.add(DeltaMin.in(N1)>=0);
    
    
    Constraint<> OneBin("OneBin");
    OneBin = bin.in_matrix(1, 1);
    Reg.add(OneBin.in(N1)==1);

    
    Constraint<> Norm2("Norm2");
    Norm2 += delta - pow(new_x1.from(cells) - x2.to(cells),2) - pow(new_y1.from(cells) - y2.to(cells),2) - pow(new_z1.from(cells) - z2.to(cells),2);
    Reg.add(Norm2.in(cells)>=0);
    
    
    
    auto ids1 = yaw.repeat_id(cells.size());
    
    Constraint<> x_rot1("x_rot1");
    x_rot1 += new_x1 - x_shift.in(ids1);
    x_rot1 -= (x1.in(N1))*cos(yaw.in(ids1))*cos(roll.in(ids1)) + (y1.in(N1))*(cos(yaw.in(ids1))*sin(roll.in(ids1))*sin(pitch.in(ids1)) - sin(yaw.in(ids1))*cos(pitch.in(ids1))) + (z1.in(N1))*(cos(yaw.in(ids1))*sin(roll.in(ids1))*cos(pitch.in(ids1)) + sin(yaw.in(ids1))*sin(pitch.in(ids1)));
    Reg.add(x_rot1.in(N1)==0);
    
    
    
    Constraint<> y_rot1("y_rot1");
    y_rot1 += new_y1 - y_shift.in(ids1);
    y_rot1 -= (x1.in(N1))*sin(yaw.in(ids1))*cos(roll.in(ids1)) + (y1.in(N1))*(sin(yaw.in(ids1))*sin(roll.in(ids1))*sin(pitch.in(ids1)) + cos(yaw.in(ids1))*cos(pitch.in(ids1))) + (z1.in(N1))*(sin(yaw.in(ids1))*sin(roll.in(ids1))*cos(pitch.in(ids1)) - cos(yaw.in(ids1))*sin(pitch.in(ids1)));
    Reg.add(y_rot1.in(N1)==0);
    
    Constraint<> z_rot1("z_rot1");
    z_rot1 += new_z1 - z_shift.in(ids1);
    z_rot1 -= (x1.in(N1))*sin(-1*roll.in(ids1)) + (y1.in(N1))*(cos(roll.in(ids1))*sin(pitch.in(ids1))) + (z1.in(N1))*(cos(roll.in(ids1))*cos(pitch.in(ids1)));
    Reg.add(z_rot1.in(N1)==0);
    Reg.min(sum(delta_min));
    
//    Reg.print();
    
    solver<> S(Reg,bonmin);
    S.run();
//    Reg.print_solution();
    //        S.run(0, 1e-10, 1000);
    
    
    //        for (int i = 0; i<500; i++) {
    //            pre_x.add_val(x_rot1.eval(i));
    //            pre_y.add_val(y_rot1.eval(i));
    //            pre_z.add_val(z_rot1.eval(i));
    //            x_uav.add_val(x_uav1.eval(i));
    //            y_uav.add_val(y_uav1.eval(i));
    //            z_uav.add_val(z_uav1.eval(i));
    //        }
    //        for (int i = 0; i<500; i++) {
    //            pre_x.add_val(x_rot2.eval(i));
    //            pre_y.add_val(y_rot2.eval(i));
    //            pre_z.add_val(z_rot2.eval(i));
    //            x_uav.add_val(x_uav2.eval(i));
    //            y_uav.add_val(y_uav2.eval(i));
    //            z_uav.add_val(z_uav2.eval(i));
    //        }
    //    M.print_solution();
    
    DebugOn("Pitch (degrees) = " << pitch.eval()*180/pi << endl);
    DebugOn("Roll (degrees) = " << roll.eval()*180/pi << endl);
    DebugOn("Yaw (degrees) = " << yaw.eval()*180/pi << endl);
    DebugOn("Pitch = " << pitch.eval() << endl);
    DebugOn("Roll = " << roll.eval() << endl);
    DebugOn("Yaw = " << yaw.eval() << endl);
    DebugOn("x shift = " << x_shift.eval() << endl);
    DebugOn("y shift = " << y_shift.eval() << endl);
    DebugOn("z shift = " << z_shift.eval() << endl);
    roll_1 = roll.eval()*180/pi;
    pitch_1 = pitch.eval()*180/pi;
    yaw_1 = yaw.eval()*180/pi;
    return {roll_1, pitch_1, yaw_1, x_shift.eval(), y_shift.eval(), z_shift.eval()};
}


/* Run the Global ARMO model for registration */
tuple<double,double,double,double,double,double> run_ARMO_Global(bool bypass, string axis, const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data){
    double angle_max = 1, shift_max = 0.2;
    double roll_1 = 0, yaw_1 = 0, pitch_1 = 0;
    int nb_pairs = 0, min_nb_pairs = numeric_limits<int>::max(), max_nb_pairs = 0, av_nb_pairs = 0;
    size_t nm = point_cloud_model.size(), nd = point_cloud_data.size();
    vector<pair<double,double>> min_max_data;
    vector<vector<pair<double,double>>> min_max_model(nm);
    vector<int> nb_neighbors(nd);
    vector<vector<int>> neighbors(nd);
    vector<double> zeros = {0,0,0};
    
    param<> x1("x1"), x2("x2"), y1("y1"), y2("y2"), z1("z1"), z2("z2");
    param<> x_uav1("x_uav1"), y_uav1("y_uav1"), z_uav1("z_uav1");
    param<> x_uav2("x_uav2"), y_uav2("y_uav2"), z_uav2("z_uav2");
    //        return 0;
    int m = av_nb_pairs;
    //            int m = 1;
    vector<double> min_dist(nd,numeric_limits<double>::max());
    vector<int> nearest(nd);
    vector<string> nearest_id(nd);
    string i_str, j_str;
    indices Pairs("Pairs"), cells("cells");
    map<int,int> n2_map;
    int idx1 = 0;
    int idx2 = 0;
    int nb_max_neigh = 1;
    double dist_sq = 0;
    /* Compute nearest points in data point cloud */
    for (auto i = 0; i<nd; i++) {
        i_str = to_string(i+1);
        x1.add_val(i_str,point_cloud_data.at(i).at(0));
        y1.add_val(i_str,point_cloud_data.at(i).at(1));
        z1.add_val(i_str,point_cloud_data.at(i).at(2));
    }
    for (auto j = 0; j<nm; j++) {
        j_str = to_string(j+1);
        x2.add_val(j_str,point_cloud_model.at(j).at(0));
        y2.add_val(j_str,point_cloud_model.at(j).at(1));
        z2.add_val(j_str,point_cloud_model.at(j).at(2));
    }
    for (auto i = 0; i< nd; i++) {
        double min_dist = numeric_limits<double>::max();
        for (auto j = 0; j< nm; j++) {
            dist_sq = std::pow(point_cloud_data.at(i).at(0) - point_cloud_model.at(j).at(0),2) + std::pow(point_cloud_data.at(i).at(1) - point_cloud_model.at(j).at(1),2) + std::pow(point_cloud_data.at(i).at(2) - point_cloud_model.at(j).at(2),2);
            if(min_dist>dist_sq){
                min_dist = dist_sq;
                j_str = to_string(j+1);
                nearest_id[i] = j_str;
            }
        }
    }
    idx1 = 0;
    indices N1("N1"),N2("N2");
    DebugOn("nd = " << nd << endl);
    DebugOn("nm = " << nm << endl);
    
    N1 = range(1,nd);
    N2 = range(1,nm);
    cells = indices(N1,N2);
//    cells.print();
//    for (auto i = 0; i<nd; i++) {
//        i_str = to_string(i+1);
//        for (auto j = 0; j<nm; j++) {
//            j_str = to_string(j+1);
//            cells.add(i_str+","+j_str);
//        }
//    }
    Model<> Reg("Reg");
    var<> new_x1("new_x1"), new_y1("new_y1"), new_z1("new_z1");
    var<> x_diff("x_diff", pos_), y_diff("y_diff", pos_), z_diff("z_diff", pos_);
    
    //            var<> yaw("yaw", thetaz, thetaz), pitch("pitch", thetax, thetax), roll("roll", thetay, thetay);
    //            var<> x_shift("x_shift", 0.2163900, 0.2163900), y_shift("y_shift", -0.1497952, -0.1497952), z_shift("z_shift", 0.0745708, 0.0745708);
    var<> cosr("cosr",  std::cos(angle_max), 1), sinr("sinr", -std::sin(angle_max), std::sin(angle_max));
    var<> cosp("cosp",  std::cos(angle_max), 1), sinp("sinp", -std::sin(angle_max), std::sin(angle_max));
    var<> cosy("cosy",  std::cos(angle_max), 1), siny("siny", -std::sin(angle_max), std::sin(angle_max));
    var<> cosy_sinr("cosy_sinr", -std::sin(angle_max), 1), siny_sinr("siny_sinr", -std::sin(angle_max)*std::sin(angle_max), std::sin(angle_max)*std::sin(angle_max));
    var<> yaw("yaw", -angle_max, angle_max), pitch("pitch", -angle_max, angle_max), roll("roll", -angle_max, angle_max);
    var<> x_shift("x_shift", -shift_max, shift_max), y_shift("y_shift", -shift_max, shift_max), z_shift("z_shift", -shift_max, shift_max);
//                    var<> yaw("yaw", -1e-6, 1e-6), pitch("pitch",-1e-6, 1e-6), roll("roll", -1e-6, 1e-6);
//                    var<> x_shift("x_shift", 0, 0), y_shift("y_shift", 0, 0), z_shift("z_shift", 0, 0);
    var<> delta("delta", pos_);
    var<> delta_min("delta_min", pos_);
    var<int> bin("bin",0,1);
    Reg.add(bin.in(cells));
    Reg.add(delta.in(cells), delta_min.in(N1));
    Reg.add(yaw.in(R(1)),pitch.in(R(1)),roll.in(R(1)));
    Reg.add(cosr.in(R(1)),cosp.in(R(1)),cosy.in(R(1)));
    Reg.add(sinr.in(R(1)),sinp.in(R(1)),siny.in(R(1)));
    Reg.add(cosy_sinr.in(R(1)),siny_sinr.in(R(1)));
    Reg.add(x_shift.in(R(1)),y_shift.in(R(1)),z_shift.in(R(1)));
    Reg.add(new_x1.in(N1), new_y1.in(N1), new_z1.in(N1));
//        Reg.add(x_diff.in(cells), y_diff.in(cells), z_diff.in(cells));
    //                Reg.add(z_diff.in(cells));
    DebugOn("There are " << cells.size() << " cells" << endl);
    
//    for (int i = 0; i<nd; i++) {
//        vector<var<>> delta_vec(nm);
//        for (int j = 0; j<nm; j++) {
//            delta_vec[j] = delta(to_string(i+1)+","+to_string(j+1));
//        }
//        Constraint<> DeltaMin("DeltaMin_"+to_string(i));
//        DeltaMin += delta_min[i+1] - min(delta_vec);
//        Reg.add(DeltaMin==0);
//    }

    Constraint<> DeltaMin("DeltaMin");
    DeltaMin = delta_min;
    DeltaMin -= bin.in_matrix(1, 1)*delta.in_matrix(1, 1);
    Reg.add(DeltaMin.in(N1)==0);
    
    
    Constraint<> OneBin("OneBin");
    OneBin = bin.in_matrix(1, 1);
    Reg.add(OneBin.in(N1)==1);
    
    Constraint<> Norm2("Norm2");
    Norm2 += delta - pow(new_x1.from(cells) - x2.to(cells),2) - pow(new_y1.from(cells) - y2.to(cells),2) - pow(new_z1.from(cells) - z2.to(cells),2);
    Reg.add(Norm2.in(cells)>=0);
    
    Constraint<> trigR("trigR");
    trigR = pow(cosr,2) + pow(sinr,2);
    Reg.add(trigR==1);

    Constraint<> trigP("trigP");
    trigP = pow(cosp,2) + pow(sinp,2);
    Reg.add(trigP==1);

    Constraint<> trigY("trigY");
    trigY = pow(cosy,2) + pow(siny,2);
    Reg.add(trigY==1);
    
//    Constraint<> cos_roll("cos_roll");
//    cos_roll = cosr - cos(roll);
//    Reg.add(cos_roll==0);
//
//    Constraint<> sin_roll("sin_roll");
//    sin_roll = sinr - sin(roll);
//    Reg.add(sin_roll==0);
//
//    Constraint<> cos_pitch("cos_pitch");
//    cos_pitch = cosp - cos(pitch);
//    Reg.add(cos_pitch==0);
//
//    Constraint<> sin_pitch("sin_pitch");
//    sin_pitch = sinp - sin(pitch);
//    Reg.add(sin_pitch==0);
//
//    Constraint<> cos_yaw("cos_yaw");
//    cos_yaw = cosy - cos(yaw);
//    Reg.add(cos_yaw==0);
//
//    Constraint<> sin_yaw("sin_yaw");
//    sin_yaw = siny - sin(yaw);
//    Reg.add(sin_yaw==0);
    
    Constraint<> cosy_sinr_prod("cosy_sinr");
    cosy_sinr_prod = cosy_sinr - cosy*sinr;
    Reg.add(cosy_sinr_prod==0);
    
    Constraint<> siny_sinr_prod("siny_sinr");
    siny_sinr_prod = siny_sinr - siny*sinr;
    Reg.add(siny_sinr_prod==0);
    
    auto ids1 = cosy.repeat_id(cells.size());
    
    /* alpha = yaw_, beta = pitch_ and gamma = roll_ */
    Constraint<> x_rot1("x_rot1");
    x_rot1 += new_x1 - x_shift.in(ids1);
    x_rot1 -= (x1.in(N1))*cosy.in(ids1)*cosr.in(ids1) + (y1.in(N1))*(cosy_sinr.in(ids1)*sinp.in(ids1) - siny.in(ids1)*cosp.in(ids1)) + (z1.in(N1))*(cosy_sinr.in(ids1)*cosp.in(ids1) + siny.in(ids1)*sinp.in(ids1));
    Reg.add(x_rot1.in(N1)==0);
    
    
    Constraint<> y_rot1("y_rot1");
    y_rot1 += new_y1 - y_shift.in(ids1);
    y_rot1 -= (x1.in(N1))*siny.in(ids1)*cosr.in(ids1) + (y1.in(N1))*(siny_sinr.in(ids1)*sinp.in(ids1) + cosy.in(ids1)*cosp.in(ids1)) + (z1.in(N1))*(siny_sinr.in(ids1)*cosp.in(ids1) - cosy.in(ids1)*sinp.in(ids1));
    Reg.add(y_rot1.in(N1)==0);
    
    Constraint<> z_rot1("z_rot1");
    z_rot1 += new_z1 - z_shift.in(ids1);
    z_rot1 -= (x1.in(N1))*-1*sinr.in(ids1) + (y1.in(N1))*(cosr.in(ids1)*sinp.in(ids1)) + (z1.in(N1))*(cosr.in(ids1)*cosp.in(ids1));
    Reg.add(z_rot1.in(N1)==0);
    
    
    //    M.min(sum(z_diff)/nb_overlap);
    
    //        M.min(sum(z_diff));
    if(axis == "full")
        //            Reg.min(sum(x_diff) + sum(y_diff) + sum(z_diff));
        //                Reg.min(sum(x_diff)/cells.size() + sum(y_diff)/cells.size() + sum(z_diff)/cells.size());
        Reg.min(sum(delta_min));
    else if(axis == "x")
        Reg.min(sum(x_diff)/cells.size());
    else if (axis == "y")
        Reg.min(sum(y_diff)/cells.size());
    else
        Reg.min(sum(z_diff)/cells.size());
    
    //                Reg.min(sum(x_diff)/cells.size() + sum(y_diff)/cells.size() + sum(z_diff)/cells.size());
    
        Reg.print();
    
    solver<> S(Reg,gurobi);
    S.run();
    Reg.print_solution();
//        S.run(0, 1e-10, 1000);
    
    
    //        for (int i = 0; i<500; i++) {
    //            pre_x.add_val(x_rot1.eval(i));
    //            pre_y.add_val(y_rot1.eval(i));
    //            pre_z.add_val(z_rot1.eval(i));
    //            x_uav.add_val(x_uav1.eval(i));
    //            y_uav.add_val(y_uav1.eval(i));
    //            z_uav.add_val(z_uav1.eval(i));
    //        }
    //        for (int i = 0; i<500; i++) {
    //            pre_x.add_val(x_rot2.eval(i));
    //            pre_y.add_val(y_rot2.eval(i));
    //            pre_z.add_val(z_rot2.eval(i));
    //            x_uav.add_val(x_uav2.eval(i));
    //            y_uav.add_val(y_uav2.eval(i));
    //            z_uav.add_val(z_uav2.eval(i));
    //        }
    //    M.print_solution();
//    roll.in(R(1));pitch.in(R(1));yaw.in(R(1));
    pitch = std::asin(sinp.eval());
    roll = std::asin(sinr.eval());
    yaw = std::asin(siny.eval());
    DebugOn("Pitch (degrees) = " << pitch.eval()*180/pi << endl);
    DebugOn("Roll (degrees) = " << roll.eval()*180/pi << endl);
    DebugOn("Yaw (degrees) = " << yaw.eval()*180/pi << endl);
    DebugOn("Pitch = " << pitch.eval() << endl);
    DebugOn("Roll = " << roll.eval() << endl);
    DebugOn("Yaw = " << yaw.eval() << endl);
    DebugOn("x shift = " << x_shift.eval() << endl);
    DebugOn("y shift = " << y_shift.eval() << endl);
    DebugOn("z shift = " << z_shift.eval() << endl);
    roll_1 = roll.eval()*180/pi;
    pitch_1 = pitch.eval()*180/pi;
    yaw_1 = yaw.eval()*180/pi;
    return {roll_1, pitch_1, yaw_1, x_shift.eval(), y_shift.eval(), z_shift.eval()};
}

/* Run the ARMO model for registration */
tuple<double,double,double,double,double,double> run_ARMO(bool bypass, string axis, const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data){
    auto thetax = atan2(-0.0081669, -0.0084357)/2;
    auto thetay = atan2(0.9999311, std::sqrt(0.0081669*0.0081669+0.0084357*0.0084357))/2;
    auto thetaz = atan2(-0.0081462,-0.0084556)/2;
    DebugOff("thetax = " << thetax << endl);
    DebugOff("thetay = " << thetay << endl);
    DebugOff("thetaz = " << thetaz << endl);
    
    if(!bypass){
        double angle_max = 3, shift_max = 0.25;
        double roll_1 = 0, yaw_1 = 0, pitch_1 = 0;
        int nb_pairs = 0, min_nb_pairs = numeric_limits<int>::max(), max_nb_pairs = 0, av_nb_pairs = 0;
        size_t nm = point_cloud_model.size(), nd = point_cloud_data.size();
        vector<pair<double,double>> min_max_data;
        vector<vector<pair<double,double>>> min_max_model(nm);
        vector<int> nb_neighbors(nd);
        vector<vector<int>> neighbors(nd);
        vector<double> zeros = {0,0,0};
        
        param<> x1("x1"), x2("x2"), y1("y1"), y2("y2"), z1("z1"), z2("z2");
        param<> x_uav1("x_uav1"), y_uav1("y_uav1"), z_uav1("z_uav1");
        param<> x_uav2("x_uav2"), y_uav2("y_uav2"), z_uav2("z_uav2");
        //        return 0;
        bool solve_lidar_cube = false, solve_lidar_iter = !solve_lidar_cube;
        int m = av_nb_pairs;
        //            int m = 1;
        vector<double> min_dist(nd,numeric_limits<double>::max());
        vector<int> nearest(nd);
        vector<string> nearest_id(nd);
        string i_str, j_str;
        indices Pairs("Pairs"), cells("cells");
        map<int,int> n2_map;
        int idx1 = 0;
        int idx2 = 0;
        int nb_max_neigh = 1;
        double dist_sq = 0;
        if(solve_lidar_cube)
            nb_max_neigh = m;
        /* Compute nearest points in data point cloud */
        for (auto i = 0; i<nd; i++) {
            i_str = to_string(i+1);
            x1.add_val(i_str,point_cloud_data.at(i).at(0));
            y1.add_val(i_str,point_cloud_data.at(i).at(1));
            z1.add_val(i_str,point_cloud_data.at(i).at(2));
        }
        for (auto j = 0; j<nm; j++) {
            j_str = to_string(j+1);
            x2.add_val(j_str,point_cloud_model.at(j).at(0));
            y2.add_val(j_str,point_cloud_model.at(j).at(1));
            z2.add_val(j_str,point_cloud_model.at(j).at(2));
        }
        for (auto i = 0; i< nd; i++) {
            double min_dist = numeric_limits<double>::max();
            for (auto j = 0; j< nm; j++) {
                if(axis=="full")
                    dist_sq = std::pow(point_cloud_data.at(i).at(0) - point_cloud_model.at(j).at(0),2) + std::pow(point_cloud_data.at(i).at(1) - point_cloud_model.at(j).at(1),2) + std::pow(point_cloud_data.at(i).at(2) - point_cloud_model.at(j).at(2),2);
                else if (axis == "z")
                    dist_sq = std::pow(point_cloud_data.at(i).at(0) - point_cloud_model.at(j).at(0),2) + std::pow(point_cloud_data.at(i).at(1) - point_cloud_model.at(j).at(1),2);
                else if(axis=="y")
                    dist_sq = std::pow(point_cloud_data.at(i).at(0) - point_cloud_model.at(j).at(0),2) + std::pow(point_cloud_data.at(i).at(2) - point_cloud_model.at(j).at(2),2);
                else
                    dist_sq = std::pow(point_cloud_data.at(i).at(1) - point_cloud_model.at(j).at(1),2) + std::pow(point_cloud_data.at(i).at(2) - point_cloud_model.at(j).at(2),2);
                
                if(min_dist>dist_sq){
                    min_dist = dist_sq;
                    j_str = to_string(j+1);
                    nearest_id[i] = j_str;
                }
            }
        }
        idx1 = 0;
        indices N1("N1"),N2("N2");
        DebugOn("nd = " << nd << endl);
        DebugOn("nm = " << nm << endl);
        
        N1 = range(1,nd);
        N2 = range(1,nm);
        if(solve_lidar_iter){
            for (auto i = 0; i<nd; i++) {
                i_str = to_string(i+1);
                j_str = nearest_id[i];
                cells.add(i_str+","+j_str);
            }
            Model<> Reg("Reg");
            var<> new_x1("new_x1"), new_y1("new_y1"), new_z1("new_z1");
            var<> x_diff("x_diff", pos_), y_diff("y_diff", pos_), z_diff("z_diff", pos_);
            
//            var<> yaw("yaw", thetaz, thetaz), pitch("pitch", thetax, thetax), roll("roll", thetay, thetay);
//            var<> x_shift("x_shift", 0.2163900, 0.2163900), y_shift("y_shift", -0.1497952, -0.1497952), z_shift("z_shift", 0.0745708, 0.0745708);
            var<> yaw("yaw", -angle_max, angle_max), pitch("pitch", -angle_max, angle_max), roll("roll", -angle_max, angle_max);
            var<> x_shift("x_shift", -shift_max, shift_max), y_shift("y_shift", -shift_max, shift_max), z_shift("z_shift", -shift_max, shift_max);
//            var<> yaw("yaw", 0, 0), pitch("pitch", 0, 0), roll("roll", 0, 0);
//            var<> x_shift("x_shift", 0, 0), y_shift("y_shift", 0, 0), z_shift("z_shift", 0, 0);
            var<> delta("delta", pos_);
            Reg.add(delta.in(cells));
            Reg.add(yaw.in(R(1)),pitch.in(R(1)),roll.in(R(1)));
            Reg.add(x_shift.in(R(1)),y_shift.in(R(1)),z_shift.in(R(1)));
            Reg.add(new_x1.in(N1), new_y1.in(N1), new_z1.in(N1));
            Reg.add(x_diff.in(cells), y_diff.in(cells), z_diff.in(cells));
            //                Reg.add(z_diff.in(cells));
            DebugOn("There are " << cells.size() << " cells" << endl);
            
            if(axis == "full"){
                Constraint<> Norm2("Norm2");
                Norm2 += delta - pow(new_x1.from(cells) - x2.to(cells),2) - pow(new_y1.from(cells) - y2.to(cells),2) - pow(new_z1.from(cells) - z2.to(cells),2);
                Reg.add(Norm2.in(cells)>=0);
            }
            else if(axis == "x"){
                Constraint<> Norm2("Norm2");
                Norm2 += delta - pow(new_y1.from(cells) - y2.to(cells),2) - pow(new_z1.from(cells) - z2.to(cells),2);
                Reg.add(Norm2.in(cells)>=0);
            }
            else if (axis == "y"){
                Constraint<> Norm2("Norm2");
                Norm2 += delta - pow(new_x1.from(cells) - x2.to(cells),2) - pow(new_z1.from(cells) - z2.to(cells),2);
                Reg.add(Norm2.in(cells)>=0);
            }
            else {
                Constraint<> Norm2("Norm2");
                Norm2 += delta - pow(new_x1.from(cells) - x2.to(cells),2) - pow(new_y1.from(cells) - y2.to(cells),2);
                Reg.add(Norm2.in(cells)>=0);
            }

            
//            Constraint<> x_abs1("x_abs1");
//            x_abs1 += x_diff - (new_x1.from(cells) - x2.to(cells));
//            Reg.add(x_abs1.in(cells)>=0);
//
//            Constraint<> x_abs2("x_abs2");
//            x_abs2 += x_diff - (x2.to(cells) - new_x1.from(cells));
//            Reg.add(x_abs2.in(cells)>=0);
//
//            Constraint<> y_abs1("y_abs1");
//            y_abs1 += y_diff - (new_y1.from(cells) - y2.to(cells));
//            Reg.add(y_abs1.in(cells)>=0);
//
//            Constraint<> y_abs2("y_abs2");
//            y_abs2 += y_diff - (y2.to(cells) - new_y1.from(cells));
//            Reg.add(y_abs2.in(cells)>=0);
//
//            Constraint<> z_abs1("z_abs1");
//            z_abs1 += z_diff - (new_z1.from(cells) - z2.to(cells));
//            Reg.add(z_abs1.in(cells)>=0);
//
//            Constraint<> z_abs2("z_abs2");
//            z_abs2 += z_diff - (z2.to(cells) - new_z1.from(cells));
//            Reg.add(z_abs2.in(cells)>=0);
            
            auto ids1 = yaw.repeat_id(cells.size());
            
            /* alpha = yaw_, beta = pitch_ and gamma = roll_ */
            Constraint<> x_rot1("x_rot1");
            x_rot1 += new_x1 - x_shift.in(ids1);
            x_rot1 -= (x1.in(N1))*cos(yaw.in(ids1))*cos(roll.in(ids1)) + (y1.in(N1))*(cos(yaw.in(ids1))*sin(roll.in(ids1))*sin(pitch.in(ids1)) - sin(yaw.in(ids1))*cos(pitch.in(ids1))) + (z1.in(N1))*(cos(yaw.in(ids1))*sin(roll.in(ids1))*cos(pitch.in(ids1)) + sin(yaw.in(ids1))*sin(pitch.in(ids1)));
            Reg.add(x_rot1.in(N1)==0);
            
            
            
            Constraint<> y_rot1("y_rot1");
            y_rot1 += new_y1 - y_shift.in(ids1);
            y_rot1 -= (x1.in(N1))*sin(yaw.in(ids1))*cos(roll.in(ids1)) + (y1.in(N1))*(sin(yaw.in(ids1))*sin(roll.in(ids1))*sin(pitch.in(ids1)) + cos(yaw.in(ids1))*cos(pitch.in(ids1))) + (z1.in(N1))*(sin(yaw.in(ids1))*sin(roll.in(ids1))*cos(pitch.in(ids1)) - cos(yaw.in(ids1))*sin(pitch.in(ids1)));
            Reg.add(y_rot1.in(N1)==0);
            
            Constraint<> z_rot1("z_rot1");
            z_rot1 += new_z1 - z_shift.in(ids1);
            z_rot1 -= (x1.in(N1))*sin(-1*roll.in(ids1)) + (y1.in(N1))*(cos(roll.in(ids1))*sin(pitch.in(ids1))) + (z1.in(N1))*(cos(roll.in(ids1))*cos(pitch.in(ids1)));
            Reg.add(z_rot1.in(N1)==0);
            
            
            //    M.min(sum(z_diff)/nb_overlap);
            
            //        M.min(sum(z_diff));
//            if(axis == "full")
                //            Reg.min(sum(x_diff) + sum(y_diff) + sum(z_diff));
//                Reg.min(sum(x_diff)/cells.size() + sum(y_diff)/cells.size() + sum(z_diff)/cells.size());
                        Reg.min(sum(delta));
//            else if(axis == "x")
//                Reg.min(sum(x_diff)/cells.size());
//            else if (axis == "y")
//                Reg.min(sum(y_diff)/cells.size());
//            else
//                Reg.min(sum(z_diff)/cells.size());
            
            //                Reg.min(sum(x_diff)/cells.size() + sum(y_diff)/cells.size() + sum(z_diff)/cells.size());
            
            //    M.print();
            
            solver<> S(Reg,ipopt);
            S.run(0, 1e-10, 1000);
            DebugOn("Objective = " << Reg.get_obj_val() << endl);
            
            
            //        for (int i = 0; i<500; i++) {
            //            pre_x.add_val(x_rot1.eval(i));
            //            pre_y.add_val(y_rot1.eval(i));
            //            pre_z.add_val(z_rot1.eval(i));
            //            x_uav.add_val(x_uav1.eval(i));
            //            y_uav.add_val(y_uav1.eval(i));
            //            z_uav.add_val(z_uav1.eval(i));
            //        }
            //        for (int i = 0; i<500; i++) {
            //            pre_x.add_val(x_rot2.eval(i));
            //            pre_y.add_val(y_rot2.eval(i));
            //            pre_z.add_val(z_rot2.eval(i));
            //            x_uav.add_val(x_uav2.eval(i));
            //            y_uav.add_val(y_uav2.eval(i));
            //            z_uav.add_val(z_uav2.eval(i));
            //        }
            //    M.print_solution();
            
            DebugOn("Pitch (degrees) = " << pitch.eval()*180/pi << endl);
            DebugOn("Roll (degrees) = " << roll.eval()*180/pi << endl);
            DebugOn("Yaw (degrees) = " << yaw.eval()*180/pi << endl);
            DebugOn("Pitch = " << pitch.eval() << endl);
            DebugOn("Roll = " << roll.eval() << endl);
            DebugOn("Yaw = " << yaw.eval() << endl);
            DebugOn("x shift = " << x_shift.eval() << endl);
            DebugOn("y shift = " << y_shift.eval() << endl);
            DebugOn("z shift = " << z_shift.eval() << endl);
            roll_1 = roll.eval()*180/pi;
            pitch_1 = pitch.eval()*180/pi;
            yaw_1 = yaw.eval()*180/pi;
            return {roll_1, pitch_1, yaw_1, x_shift.eval(), y_shift.eval(), z_shift.eval()};
        }
    }
    else
        return {thetay*180/pi, thetax*180/pi, thetaz*180/pi, 0.2163900/2, -0.1497952/2, 0.0745708/2};
}

tuple<double,double,double> run_ARMO(string axis, const vector<double>& x_vec1, const vector<double>& y_vec1, const vector<double>& z_vec1, const vector<double>& x_vec2, const vector<double>& y_vec2, const vector<double>& z_vec2, const vector<double>& x_shift1, const vector<double>& y_shift1, const vector<double>& z_shift1, const vector<double>& x_shift2, const vector<double>& y_shift2, const vector<double>& z_shift2, const vector<vector<double>>& point_cloud1, const vector<vector<double>>& point_cloud2, const vector<vector<double>>& uav1, const vector<vector<double>>& uav2){
    double angle_max = 0.05;
    int nb_pairs = 0, min_nb_pairs = numeric_limits<int>::max(), max_nb_pairs = 0, av_nb_pairs = 0;
    vector<pair<double,double>> min_max1;
    vector<vector<pair<double,double>>> min_max2(x_vec2.size());
    vector<int> nb_neighbors(x_vec1.size());
    vector<vector<int>> neighbors;
    /* Compute cube for all points in point cloud 2 */
    for (auto i = 0; i<x_vec2.size(); i++) {
        min_max2[i] = get_min_max(angle_max, point_cloud2.at(i), uav2.at(i));
    }
    double roll_1 = 0, yaw_1 = 0, pitch_1 = 0;
    bool bypass = false;
    if(!bypass){
        /* Check if cubes intersect */
        neighbors.resize(x_vec1.size());
        for (auto i = 0; i<x_vec1.size(); i++) {
            nb_pairs = 0;
            min_max1 = get_min_max(angle_max, point_cloud1.at(i), uav1.at(i));
            DebugOff("For point (" << point_cloud1.at(i).at(0) << "," <<  point_cloud1.at(i).at(1) << "," << point_cloud1.at(i).at(2) << "): ");
            DebugOff("\n neighbors in umbrella : \n");
            for (size_t j = 0; j < x_vec2.size(); j++){
                if(intersect(min_max1, min_max2[j])){ /* point is in umbrella */
                    nb_pairs++;
                    neighbors[i].push_back(j);
                    DebugOff("(" << point_cloud2.at(j).at(0) << "," <<  point_cloud2.at(j).at(1) << "," << point_cloud2.at(j).at(2) << ")\n");
                }
            }
            
            DebugOff("nb points in umbrella = " << nb_pairs << endl);
            if(nb_pairs>max_nb_pairs)
                max_nb_pairs = nb_pairs;
                if(nb_pairs<min_nb_pairs)
                    min_nb_pairs = nb_pairs;
                    av_nb_pairs += nb_pairs;
                    
                    //        std::cout << "For point (" << point_cloud1.at(i).at(0) << "," <<  point_cloud1.at(i).at(1) << "," << point_cloud1.at(i).at(2) << ")"<< " knnSearch(n="<<m<<"): \n";
                    //        for (size_t k = 0; k < m; k++)
                    //            std::cout << "ret_index["<<k<<"]=" << ret_indexes[k] << " out_dist_sqr=" << out_dists_sqr[k] << " point = (" << point_cloud2.at(ret_indexes[k]).at(0) << "," <<  point_cloud2.at(ret_indexes[k]).at(1) << "," << point_cloud2.at(ret_indexes[k]).at(2) << ")" << std::endl;
                    nb_neighbors[i] = nb_pairs;
                    }
        av_nb_pairs /= x_vec1.size();
        DebugOn("Min nb of Pairs = " << min_nb_pairs << endl);
        DebugOn("Max nb of Pairs = " << max_nb_pairs << endl);
        DebugOn("Average nb of Pairs = " << av_nb_pairs << endl);
        param<> x1("x1"), x2("x2"), y1("y1"), y2("y2"), z1("z1"), z2("z2");
        param<> x_uav1("x_uav1"), y_uav1("y_uav1"), z_uav1("z_uav1");
        param<> x_uav2("x_uav2"), y_uav2("y_uav2"), z_uav2("z_uav2");
        //        return 0;
        bool solve_lidar_cube = false, solve_lidar_iter = !solve_lidar_cube;
        int m = av_nb_pairs;
        //            int m = 1;
        vector<double> min_dist(x_vec1.size(),numeric_limits<double>::max());
        vector<int> nearest(x_vec1.size());
        vector<string> nearest_id(x_vec1.size());
        string i_str, j_str;
        indices Pairs("Pairs"), cells("cells");
        map<int,int> n2_map;
        int idx1 = 0;
        int idx2 = 0;
        int nb_max_neigh = 1;
        double dist_sq = 0;
        if(solve_lidar_cube)
            nb_max_neigh = m;
        /* Keep points with neighbors >= m */
        for (auto i = 0; i<x_vec1.size(); i++) {
            if(solve_lidar_iter)
                nb_max_neigh = 1;
            else
                nb_max_neigh = m;
            if(nb_neighbors[i]>=nb_max_neigh){
                i_str = to_string(idx1+1);
                x_uav1.add_val(i_str,x_shift1.at(i));
                x1.add_val(i_str,x_vec1.at(i));
                y_uav1.add_val(i_str,y_shift1.at(i));
                y1.add_val(i_str,y_vec1.at(i));
                z1.add_val(i_str,z_vec1.at(i));
                z_uav1.add_val(i_str,z_shift1.at(i));
                if(solve_lidar_iter){
                    nb_max_neigh = nb_neighbors[i];
                }
                for (auto j = 0; j<nb_max_neigh; j++) {
                    auto k = neighbors[i].at(j);
                    auto res = n2_map.find(k);
                    if(res==n2_map.end()){
                        n2_map[k] = idx2;
                        j_str = to_string(idx2+1);
                        x_uav2.add_val(j_str,x_shift2.at(k));
                        x2.add_val(j_str,x_vec2.at(k));
                        y_uav2.add_val(j_str,y_shift2.at(k));
                        y2.add_val(j_str,y_vec2.at(k));
                        z_uav2.add_val(j_str,z_shift2.at(k));
                        z2.add_val(j_str,z_vec2.at(k));
                        idx2++;
                    }
                    else {
                        j_str = to_string(res->second+1);
                    }
//                    if(axis=="x")
                        dist_sq = std::pow(x_vec1.at(i) - x_vec2.at(k),2) + std::pow(y_vec1.at(i) - y_vec2.at(k),2) + std::pow(z_vec1.at(i) - z_vec2.at(k),2);
//                    else
//                        dist_sq = std::pow(x_vec1.at(i) - x_vec2.at(k),2) + std::pow(y_vec1.at(i) - y_vec2.at(k),2);
                    if(min_dist[i]>dist_sq){
                        min_dist[i] = dist_sq;
                        nearest[i] = k;
                        nearest_id[i] = j_str;
                    }
                    
                    if(solve_lidar_cube)
                        Pairs.add(i_str+","+j_str);
                }
                idx1++;
            }
        }
        idx1 = 0;
        indices N1("N1"),N2("N2");
        if(solve_lidar_iter){
            for (auto i = 0; i<x_vec1.size(); i++) {
                if(nb_neighbors[i]>=1){
                    i_str = to_string(idx1+1);
                    j_str = nearest_id[i];
                    cells.add(i_str+","+j_str);
                    if(!N2.has_key(j_str))
                        N2.add(j_str);
                        idx1++;
                }
            }
        }
        
        
        int n1 = x1.get_dim();
        int n2 = x2.get_dim();
        DebugOn("n1 = " << n1 << endl);
        DebugOn("n2 = " << n2 << endl);
        
        N1 = range(1,n1);
        if(solve_lidar_cube)
            N2 = range(1,n2);
            indices M("M");
            M = range(1,m);
            
            DebugOn("Total size of Pairs = " << Pairs.size() << endl);
            
            
            indices NM("NM");
            NM = indices(N1,M);
            indices S1("S1"), S2("S2"), Sm1("Sm1"), S2m2("S2m2"), S3m1("S3m1"), Sm("Sm"), K("K");
            S1 = indices(N1, range(1,1));
            S2 = indices(N1, range(2,2));
            Sm = indices(N1, range(m,m));
            Sm1 = indices(N1, range(m-1,m-1));
            S2m2 = indices(N1, range(2,m-2));
            S3m1 = indices(N1, range(3,m-1));
            K = indices(N1,range(2,m-1));
            
            
            
            
            
            if (solve_lidar_cube) {
                Model<> Lidar("Lidar");
                var<> new_x1("new_x1"), new_y1("new_y1"), new_z1("new_z1");
                var<> new_x2("new_x2"), new_y2("new_y2"), new_z2("new_z2");
                //            var<> x_diff("x_diff", pos_), y_diff("y_diff", pos_), z_diff("z_diff", pos_);
                var<> mu("mu", pos_), mu_k("mu_k", pos_), delta("delta", pos_);
                //                            var<> yaw1("yaw1", -0.5*pi/180, -0.5*pi/180), pitch1("pitch1", 0.9*pi/180, 0.9*pi/180), roll1("roll1", 1.375*pi/180, 1.375*pi/180);
                var<> yaw1("yaw1", -0.1, 0.1), pitch1("pitch1", -0.1, 0.1), roll1("roll1", -0.1, 0.1);
                //                var<> yaw1("yaw1", 0.25*pi/180, 0.25*pi/180), pitch1("pitch1", 0.9*pi/180, 0.9*pi/180), roll1("roll1", -1.45*pi/180, -1.45*pi/180);
                //                var<> yaw1("yaw1", 0, 0), pitch1("pitch1", -0.5778*pi/180, -0.5778*pi/180), roll1("roll1", -1.44581*pi/180, -1.44581*pi/180);
                //                var<> yaw1("yaw1", 0, 0), pitch1("pitch1", -0.573231*pi/180, -0.573231*pi/180), roll1("roll1", -1.45338*pi/180, -1.45338*pi/180);
                //                var<> yaw1("yaw1", 0.0249847*pi/180, 0.0249847*pi/180), pitch1("pitch1", -0.507086*pi/180, -0.507086*pi/180), roll1("roll1", -1.3698*pi/180, -1.3698*pi/180);
                //                var<> yaw1("yaw1", 0, 0), pitch1("pitch1", 0, 0), roll1("roll1", 0, 0);
                var<> yaw2("yaw2", -0.1, 0.1), pitch2("pitch2", -0.1, 0.1), roll2("roll2", -0.1, 0.1);
                //            yaw1 = -0.5*pi/180;
                //            pitch1 = 0.9*pi/180;
                //            roll1 = 1.375*pi/180;
                Lidar.add(yaw1.in(R(1)),pitch1.in(R(1)),roll1.in(R(1)));
                Lidar.add(yaw2.in(R(1)),pitch2.in(R(1)),roll2.in(R(1)));
                Lidar.add(new_x1.in(N1), new_y1.in(N1), new_z1.in(N1));
                Lidar.add(new_x2.in(N2), new_y2.in(N2), new_z2.in(N2));
                //            Lidar.add(x_diff.in(NM), y_diff.in(NM), z_diff.in(NM));
                Lidar.add(delta.in(NM));
                Lidar.add(mu.in(N1));
                Lidar.add(mu_k.in(K));
                
                Constraint<> Equal_pitch("Equal_pitch");
                Equal_pitch += pitch1 - pitch2;
                Lidar.add(Equal_pitch==0);
                
                Constraint<> Opp_roll("Opp_roll");
                Opp_roll += roll1 + roll2;
                Lidar.add(Opp_roll==0);
                
                Constraint<> Opp_yaw("Opp_yaw");
                Opp_yaw += yaw1 + yaw2;
                Lidar.add(Opp_yaw==0);
                
                Constraint<> Mu_2("Mu_2");
                Mu_2 += mu_k.in(S2) - gravity::min(delta.in(S1), delta.in(S2));
                Lidar.add(Mu_2.in(S2)==0);
                
                Constraint<> Mu_k("Mu_k");
                Mu_k += mu_k.in(S3m1) - gravity::min(mu_k.in(S2m2), delta.in(S3m1));
                Lidar.add(Mu_k.in(S3m1)==0);
                
                Constraint<> Mu("Mu");
                Mu += mu - gravity::min(mu_k.in(Sm1), delta.in(Sm));
                Lidar.add(Mu.in(N1)==0);
                
                //                            Constraint<> Norm2("Norm2");
                //                            Norm2 += delta - pow(new_x1.in_ignore_ith(1, 1, Pairs) - new_x2.in_ignore_ith(0, 1, Pairs),2) - pow(new_y1.in_ignore_ith(1, 1, Pairs) - new_y2.in_ignore_ith(0, 1, Pairs),2) - pow(new_z1.in_ignore_ith(1, 1, Pairs) - new_z2.in_ignore_ith(0, 1, Pairs),2);
                //                            Lidar.add(Norm2.in(Pairs)==0);
                
                Constraint<> Norm1("Norm1");
                Norm1 += delta - abs(new_x1.in_ignore_ith(1, 1, Pairs) - new_x2.in_ignore_ith(0, 1, Pairs)) - abs(new_y1.in_ignore_ith(1, 1, Pairs) - new_y2.in_ignore_ith(0, 1, Pairs)) - abs(new_z1.in_ignore_ith(1, 1, Pairs) - new_z2.in_ignore_ith(0, 1, Pairs));
                //                Norm1 += delta - abs(new_z1.in_ignore_ith(1, 1, Pairs) - new_z2.in_ignore_ith(0, 1, Pairs));
                Lidar.add(Norm1.in(Pairs)==0);
                
                
                //            Constraint<> z_abs1("z_abs1");
                //            z_abs1 += z_diff - (new_z1.in_ignore_ith(1, 1, Pairs) - new_z2.in_ignore_ith(0, 1, Pairs));
                //            Lidar.add(z_abs1.in(Pairs)>=0);
                //
                //            Constraint<> z_abs2("z_abs2");
                //            z_abs2 += z_diff - (new_z2.in_ignore_ith(1, 0, Pairs) - new_z1.in_ignore_ith(1, 1, Pairs));
                //            Lidar.add(z_abs2.in(Pairs)>=0);
                //
                //            Constraint<> x_abs1("x_abs1");
                //            x_abs1 += x_diff - (new_x1.in_ignore_ith(1, 1, Pairs) - new_x2.in_ignore_ith(0, 1, Pairs));
                //            Lidar.add(x_abs1.in(Pairs)>=0);
                //
                //            Constraint<> x_abs2("x_abs2");
                //            x_abs2 += x_diff - (new_x2.in_ignore_ith(0, 1, Pairs) - new_x1.in_ignore_ith(1, 1, Pairs));
                //            Lidar.add(x_abs2.in(Pairs)>=0);
                //
                //
                //            Constraint<> y_abs1("y_abs1");
                //            y_abs1 += y_diff - (new_y1.in_ignore_ith(1, 1, Pairs) - new_y2.in_ignore_ith(0, 1, Pairs));
                //            Lidar.add(y_abs1.in(Pairs)>=0);
                //
                //            Constraint<> y_abs2("y_abs2");
                //            y_abs2 += y_diff - (new_y2.in_ignore_ith(0, 1, Pairs) - new_y1.in_ignore_ith(1, 1, Pairs));
                //            Lidar.add(y_abs2.in(Pairs)>=0);
                
                auto ids1 = yaw1.repeat_id(n1);
                
                /* alpha = yaw_, beta = pitch_ and gamma = roll_ */
                Constraint<> x_rot1("x_rot1");
                x_rot1 += new_x1 - x_uav1.in(N1);
                x_rot1 -= (x1.in(N1)-x_uav1.in(N1))*cos(yaw1.in(ids1))*cos(roll1.in(ids1)) + (y1.in(N1)-y_uav1.in(N1))*(cos(yaw1.in(ids1))*sin(roll1.in(ids1))*sin(pitch1.in(ids1)) - sin(yaw1.in(ids1))*cos(pitch1.in(ids1))) + (z1.in(N1)-z_uav1.in(N1))*(cos(yaw1.in(ids1))*sin(roll1.in(ids1))*cos(pitch1.in(ids1)) + sin(yaw1.in(ids1))*sin(pitch1.in(ids1)));
                Lidar.add(x_rot1.in(N1)==0);
                
                auto ids2 = yaw2.repeat_id(n2);
                
                Constraint<> x_rot2("x_rot2");
                x_rot2 += new_x2 - x_uav2.in(N2);
                x_rot2 -= (x2.in(N2)-x_uav2.in(N2))*cos(yaw2.in(ids2))*cos(roll2.in(ids2)) + (y2.in(N2)-y_uav2.in(N2))*(cos(yaw2.in(ids2))*sin(roll2.in(ids2))*sin(pitch2.in(ids2)) - sin(yaw2.in(ids2))*cos(pitch2.in(ids2))) + (z2.in(N2)-z_uav2.in(N2))*(cos(yaw2.in(ids2))*sin(roll2.in(ids2))*cos(pitch2.in(ids2)) + sin(yaw2.in(ids2))*sin(pitch2.in(ids2)));
                Lidar.add(x_rot2.in(N2)==0);
                
                
                Constraint<> y_rot1("y_rot1");
                y_rot1 += new_y1 - y_uav1.in(N1);
                y_rot1 -= (x1.in(N1)-x_uav1.in(N1))*sin(yaw1.in(ids1))*cos(roll1.in(ids1)) + (y1.in(N1)-y_uav1.in(N1))*(sin(yaw1.in(ids1))*sin(roll1.in(ids1))*sin(pitch1.in(ids1)) + cos(yaw1.in(ids1))*cos(pitch1.in(ids1))) + (z1.in(N1)-z_uav1.in(N1))*(sin(yaw1.in(ids1))*sin(roll1.in(ids1))*cos(pitch1.in(ids1)) - cos(yaw1.in(ids1))*sin(pitch1.in(ids1)));
                Lidar.add(y_rot1.in(N1)==0);
                
                Constraint<> y_rot2("y_rot2");
                y_rot2 += new_y2 - y_uav2.in(N2);
                y_rot2 -= (x2.in(N2)-x_uav2.in(N2))*sin(yaw2.in(ids2))*cos(roll2.in(ids2)) + (y2.in(N2)-y_uav2.in(N2))*(sin(yaw2.in(ids2))*sin(roll2.in(ids2))*sin(pitch2.in(ids2)) + cos(yaw2.in(ids2))*cos(pitch2.in(ids2))) + (z2.in(N2)-z_uav2.in(N2))*(sin(yaw2.in(ids2))*sin(roll2.in(ids2))*cos(pitch2.in(ids2)) - cos(yaw2.in(ids2))*sin(pitch2.in(ids2)));
                Lidar.add(y_rot2.in(N2)==0);
                
                Constraint<> z_rot1("z_rot1");
                z_rot1 += new_z1 - z_uav1.in(N1);
                z_rot1 -= (x1.in(N1)-x_uav1.in(N1))*sin(-1*roll1.in(ids1)) + (y1.in(N1)-y_uav1.in(N1))*(cos(roll1.in(ids1))*sin(pitch1.in(ids1))) + (z1.in(N1)-z_uav1.in(N1))*(cos(roll1.in(ids1))*cos(pitch1.in(ids1)));
                Lidar.add(z_rot1.in(N1)==0);
                
                
                Constraint<> z_rot2("z_rot2");
                z_rot2 += new_z2 - z_uav2.in(N2);
                z_rot2 -= (x2.in(N2)-x_uav2.in(N2))*sin(-1*roll2.in(ids2)) + (y2.in(N2)-y_uav2.in(N2))*(cos(roll2.in(ids2))*sin(pitch2.in(ids2))) + (z2.in(N2)-z_uav2.in(N2))*(cos(roll2.in(ids2))*cos(pitch2.in(ids2)));
                Lidar.add(z_rot2.in(N2)==0);
                
                //            Lidar.min(sum(mu) + 1e2*pow(yaw1,2) + 1e2*pow(roll1,2) + 1e2*pow(pitch1,2));
                //            Lidar.min(1e3*sum(mu) + (pow(yaw1,2) + pow(roll1,2) + pow(pitch1,2)));
                Lidar.min(sum(mu));
                
                //            Lidar.print();
                //            Lidar.initialize_zero();
                //            return 0;
                solver<> S(Lidar,ipopt);
                //            S.set_option("tol", 1e-10);
                //            S.run(5,1e-10);
                S.run();
                DebugOn("Pitch1 = " << pitch1.eval()*180/pi << endl);
                DebugOn("Roll1 = " << roll1.eval()*180/pi << endl);
                DebugOn("Yaw1 = " << yaw1.eval()*180/pi << endl);
                DebugOn("Pitch2 = " << pitch2.eval()*180/pi << endl);
                DebugOn("Roll2 = " << roll2.eval()*180/pi << endl);
                DebugOn("Yaw2 = " << yaw2.eval()*180/pi << endl);
                roll_1 = roll1.eval()*180/pi;
                pitch_1 = pitch1.eval()*180/pi;
                yaw_1 = yaw1.eval()*180/pi;
                //            return 0;
            }
        if (solve_lidar_iter) {
            Model<> Lidar("Lidar");
            var<> new_x1("new_x1"), new_y1("new_y1"), new_z1("new_z1");
            var<> new_x2("new_x2"), new_y2("new_y2"), new_z2("new_z2");
            var<> x_diff("x_diff", pos_), y_diff("y_diff", pos_), z_diff("z_diff", pos_);
//            var<> yaw1("yaw1", 0.25*pi/180, 0.25*pi/180), pitch1("pitch1", 0.5*pi/180, 0.5*pi/180), roll1("roll1", 0.7*pi/180, 0.7*pi/180);
//            var<> yaw1("yaw1", 0.25*pi/180, 0.25*pi/180), pitch1("pitch1", 0.9*pi/180, 0.9*pi/180), roll1("roll1", -1.45*pi/180, -1.45*pi/180);
            //                var<> yaw1("yaw1", 0, 0), pitch1("pitch1", 0, 0), roll1("roll1", 0, 0);
            //                var<> yaw1("yaw1", -0.5*pi/180, -0.5*pi/180), pitch1("pitch1", 0.9*pi/180, 0.9*pi/180), roll1("roll1", 1.375*pi/180, 1.375*pi/180);
            var<> yaw1("yaw1", -0.1, 0.1), pitch1("pitch1", -0.1, 0.1), roll1("roll1", -0.1, 0.1);
            var<> yaw2("yaw2", -0.1, 0.1), pitch2("pitch2", -0.1, 0.1), roll2("roll2", -0.1, 0.1);
            
            Lidar.add(yaw1.in(R(1)),pitch1.in(R(1)),roll1.in(R(1)));
            Lidar.add(yaw2.in(R(1)),pitch2.in(R(1)),roll2.in(R(1)));
            Lidar.add(new_x1.in(N1), new_y1.in(N1), new_z1.in(N1));
            Lidar.add(new_x2.in(N2), new_y2.in(N2), new_z2.in(N2));
            Lidar.add(x_diff.in(cells), y_diff.in(cells), z_diff.in(cells));
            //                Lidar.add(z_diff.in(cells));
            
            Constraint<> Equal_pitch("Equal_pitch");
            Equal_pitch += pitch1 - pitch2;
            Lidar.add(Equal_pitch==0);
            
            Constraint<> Opp_roll("Opp_roll");
            Opp_roll += roll1 + roll2;
            Lidar.add(Opp_roll==0);
            
            Constraint<> Opp_yaw("Opp_yaw");
            Opp_yaw += yaw1 + yaw2;
            Lidar.add(Opp_yaw==0);
            
            
            Constraint<> x_abs1("x_abs1");
            x_abs1 += x_diff - (new_x1.from(cells) - new_x2.to(cells));
            Lidar.add(x_abs1.in(cells)>=0);
            
            Constraint<> x_abs2("x_abs2");
            x_abs2 += x_diff - (new_x2.to(cells) - new_x1.from(cells));
            Lidar.add(x_abs2.in(cells)>=0);
            
            Constraint<> y_abs1("y_abs1");
            y_abs1 += y_diff - (new_y1.from(cells) - new_y2.to(cells));
            Lidar.add(y_abs1.in(cells)>=0);
            
            Constraint<> y_abs2("y_abs2");
            y_abs2 += y_diff - (new_y2.to(cells) - new_y1.from(cells));
            Lidar.add(y_abs2.in(cells)>=0);
            
            Constraint<> z_abs1("z_abs1");
            z_abs1 += z_diff - (new_z1.from(cells) - new_z2.to(cells));
            Lidar.add(z_abs1.in(cells)>=0);
            
            Constraint<> z_abs2("z_abs2");
            z_abs2 += z_diff - (new_z2.to(cells) - new_z1.from(cells));
            Lidar.add(z_abs2.in(cells)>=0);
            
            auto ids1 = yaw1.repeat_id(cells.size());
            
            /* alpha = yaw_, beta = pitch_ and gamma = roll_ */
            Constraint<> x_rot1("x_rot1");
            x_rot1 += new_x1 - x_uav1.in(N1);
            x_rot1 -= (x1.in(N1)-x_uav1.in(N1))*cos(yaw1.in(ids1))*cos(roll1.in(ids1)) + (y1.in(N1)-y_uav1.in(N1))*(cos(yaw1.in(ids1))*sin(roll1.in(ids1))*sin(pitch1.in(ids1)) - sin(yaw1.in(ids1))*cos(pitch1.in(ids1))) + (z1.in(N1)-z_uav1.in(N1))*(cos(yaw1.in(ids1))*sin(roll1.in(ids1))*cos(pitch1.in(ids1)) + sin(yaw1.in(ids1))*sin(pitch1.in(ids1)));
            Lidar.add(x_rot1.in(N1)==0);
            
            auto ids2 = yaw2.repeat_id(cells.size());
            
            Constraint<> x_rot2("x_rot2");
            x_rot2 += new_x2 - x_uav2.in(N2);
            x_rot2 -= (x2.in(N2)-x_uav2.in(N2))*cos(yaw2.in(ids2))*cos(roll2.in(ids2)) + (y2.in(N2)-y_uav2.in(N2))*(cos(yaw2.in(ids2))*sin(roll2.in(ids2))*sin(pitch2.in(ids2)) - sin(yaw2.in(ids2))*cos(pitch2.in(ids2))) + (z2.in(N2)-z_uav2.in(N2))*(cos(yaw2.in(ids2))*sin(roll2.in(ids2))*cos(pitch2.in(ids2)) + sin(yaw2.in(ids2))*sin(pitch2.in(ids2)));
            Lidar.add(x_rot2.in(N2)==0);
            
            
            Constraint<> y_rot1("y_rot1");
            y_rot1 += new_y1 - y_uav1.in(N1);
            y_rot1 -= (x1.in(N1)-x_uav1.in(N1))*sin(yaw1.in(ids1))*cos(roll1.in(ids1)) + (y1.in(N1)-y_uav1.in(N1))*(sin(yaw1.in(ids1))*sin(roll1.in(ids1))*sin(pitch1.in(ids1)) + cos(yaw1.in(ids1))*cos(pitch1.in(ids1))) + (z1.in(N1)-z_uav1.in(N1))*(sin(yaw1.in(ids1))*sin(roll1.in(ids1))*cos(pitch1.in(ids1)) - cos(yaw1.in(ids1))*sin(pitch1.in(ids1)));
            Lidar.add(y_rot1.in(N1)==0);
            
            Constraint<> y_rot2("y_rot2");
            y_rot2 += new_y2 - y_uav2.in(N2);
            y_rot2 -= (x2.in(N2)-x_uav2.in(N2))*sin(yaw2.in(ids2))*cos(roll2.in(ids2)) + (y2.in(N2)-y_uav2.in(N2))*(sin(yaw2.in(ids2))*sin(roll2.in(ids2))*sin(pitch2.in(ids2)) + cos(yaw2.in(ids2))*cos(pitch2.in(ids2))) + (z2.in(N2)-z_uav2.in(N2))*(sin(yaw2.in(ids2))*sin(roll2.in(ids2))*cos(pitch2.in(ids2)) - cos(yaw2.in(ids2))*sin(pitch2.in(ids2)));
            Lidar.add(y_rot2.in(N2)==0);
            
            Constraint<> z_rot1("z_rot1");
            z_rot1 += new_z1 - z_uav1.in(N1);
            z_rot1 -= (x1.in(N1)-x_uav1.in(N1))*sin(-1*roll1.in(ids1)) + (y1.in(N1)-y_uav1.in(N1))*(cos(roll1.in(ids1))*sin(pitch1.in(ids1))) + (z1.in(N1)-z_uav1.in(N1))*(cos(roll1.in(ids1))*cos(pitch1.in(ids1)));
            Lidar.add(z_rot1.in(N1)==0);
            
            
            Constraint<> z_rot2("z_rot2");
            z_rot2 += new_z2 - z_uav2.in(N2);
            z_rot2 -= (x2.in(N2)-x_uav2.in(N2))*sin(-1*roll2.in(ids2)) + (y2.in(N2)-y_uav2.in(N2))*(cos(roll2.in(ids2))*sin(pitch2.in(ids2))) + (z2.in(N2)-z_uav2.in(N2))*(cos(roll2.in(ids2))*cos(pitch2.in(ids2)));
            Lidar.add(z_rot2.in(N2)==0);
            
            //    M.min(sum(z_diff)/nb_overlap);
            
            //        M.min(sum(z_diff));
            if(axis == "full")
                Lidar.min(sum(x_diff)/cells.size() + sum(y_diff)/cells.size() + sum(z_diff)/cells.size());
            else if(axis == "x")
                Lidar.min(sum(x_diff)/cells.size());
            else if (axis == "y")
                Lidar.min(sum(y_diff)/cells.size());
            else
                Lidar.min(sum(z_diff)/cells.size());
            
            //                Lidar.min(sum(x_diff)/cells.size() + sum(y_diff)/cells.size() + sum(z_diff)/cells.size());
            
            //    M.print();
            
            solver<> S(Lidar,ipopt);
            S.run();
            
            
            //        for (int i = 0; i<500; i++) {
            //            pre_x.add_val(x_rot1.eval(i));
            //            pre_y.add_val(y_rot1.eval(i));
            //            pre_z.add_val(z_rot1.eval(i));
            //            x_uav.add_val(x_uav1.eval(i));
            //            y_uav.add_val(y_uav1.eval(i));
            //            z_uav.add_val(z_uav1.eval(i));
            //        }
            //        for (int i = 0; i<500; i++) {
            //            pre_x.add_val(x_rot2.eval(i));
            //            pre_y.add_val(y_rot2.eval(i));
            //            pre_z.add_val(z_rot2.eval(i));
            //            x_uav.add_val(x_uav2.eval(i));
            //            y_uav.add_val(y_uav2.eval(i));
            //            z_uav.add_val(z_uav2.eval(i));
            //        }
            //    M.print_solution();
            
            DebugOn("Pitch1 = " << pitch1.eval()*180/pi << endl);
            DebugOn("Roll1 = " << roll1.eval()*180/pi << endl);
            DebugOn("Yaw1 = " << yaw1.eval()*180/pi << endl);
            DebugOn("Pitch2 = " << pitch2.eval()*180/pi << endl);
            DebugOn("Roll2 = " << roll2.eval()*180/pi << endl);
            DebugOn("Yaw2 = " << yaw2.eval()*180/pi << endl);
            roll_1 = roll1.eval()*180/pi;
            pitch_1 = pitch1.eval()*180/pi;
            yaw_1 = yaw1.eval()*180/pi;
        }
    }
    return {roll_1, pitch_1, yaw_1};
}


void apply_rotation(double roll, double pitch, double yaw, vector<double>& x_vec1, vector<double>& y_vec1, vector<double>& z_vec1, vector<double>& x_vec2, vector<double>& y_vec2, vector<double>& z_vec2, const vector<double>& x_shift1, const vector<double>& y_shift1, const vector<double>& z_shift1, const vector<double>& x_shift2, const vector<double>& y_shift2, const vector<double>& z_shift2){
    double beta = roll*pi/180;// roll in radians
    double gamma = pitch*pi/180; // pitch in radians
    double alpha = yaw*pi/180; // yaw in radians
    double shifted_x, shifted_y, shifted_z;
    /* Apply rotation */
    for (auto i = 0; i< x_vec1.size(); i++) {
        shifted_x = x_vec1[i] - x_shift1[i];
        shifted_y = y_vec1[i] - y_shift1[i];
        shifted_z = z_vec1[i] - z_shift1[i];
        x_vec1[i] = shifted_x*cos(alpha)*cos(beta) + shifted_y*(cos(alpha)*sin(beta)*sin(gamma) - sin(alpha)*cos(gamma)) + shifted_z*(cos(alpha)*sin(beta)*cos(gamma) + sin(alpha)*sin(gamma));
        y_vec1[i] = shifted_x*sin(alpha)*cos(beta) + shifted_y*(sin(alpha)*sin(beta)*sin(gamma) + cos(alpha)*cos(gamma)) + shifted_z*(sin(alpha)*sin(beta)*cos(gamma) - cos(alpha)*sin(gamma));
        z_vec1[i] = shifted_x*(-sin(beta)) + shifted_y*(cos(beta)*sin(gamma)) + shifted_z*(cos(beta)*cos(gamma));
        x_vec1[i] += x_shift1[i];
        y_vec1[i] += y_shift1[i];
        z_vec1[i] += z_shift1[i];
    }
    beta *= -1;
    alpha *= -1;
    for (auto i = 0; i< x_vec2.size(); i++) {
        shifted_x = x_vec2[i] - x_shift2[i];
        shifted_y = y_vec2[i] - y_shift2[i];
        shifted_z = z_vec2[i] - z_shift2[i];
        x_vec2[i] = shifted_x*cos(alpha)*cos(beta) + shifted_y*(cos(alpha)*sin(beta)*sin(gamma) - sin(alpha)*cos(gamma)) + shifted_z*(cos(alpha)*sin(beta)*cos(gamma) + sin(alpha)*sin(gamma));
        y_vec2[i] = shifted_x*sin(alpha)*cos(beta) + shifted_y*(sin(alpha)*sin(beta)*sin(gamma) + cos(alpha)*cos(gamma)) + shifted_z*(sin(alpha)*sin(beta)*cos(gamma) - cos(alpha)*sin(gamma));
        z_vec2[i] = shifted_x*(-sin(beta)) + shifted_y*(cos(beta)*sin(gamma)) + shifted_z*(cos(beta)*cos(gamma));
        x_vec2[i] += x_shift2[i];
        y_vec2[i] += y_shift2[i];
        z_vec2[i] += z_shift2[i];
    }
}


void apply_rot_trans(double roll, double pitch, double yaw, double x_shift, double y_shift, double z_shift, vector<vector<double>>& point_cloud){
    double beta = roll*pi/180;// roll in radians
    double gamma = pitch*pi/180; // pitch in radians
    double alpha = yaw*pi/180; // yaw in radians
    double shifted_x, shifted_y, shifted_z;
    size_t n = point_cloud.size();
    /* Apply rotation */
    for (auto i = 0; i< n; i++) {
        shifted_x = point_cloud[i][0];
        shifted_y = point_cloud[i][1];
        shifted_z = point_cloud[i][2];
        point_cloud[i][0] = shifted_x*cos(alpha)*cos(beta) + shifted_y*(cos(alpha)*sin(beta)*sin(gamma) - sin(alpha)*cos(gamma)) + shifted_z*(cos(alpha)*sin(beta)*cos(gamma) + sin(alpha)*sin(gamma));
        point_cloud[i][1] = shifted_x*sin(alpha)*cos(beta) + shifted_y*(sin(alpha)*sin(beta)*sin(gamma) + cos(alpha)*cos(gamma)) + shifted_z*(sin(alpha)*sin(beta)*cos(gamma) - cos(alpha)*sin(gamma));
        point_cloud[i][2] = shifted_x*(-sin(beta)) + shifted_y*(cos(beta)*sin(gamma)) + shifted_z*(cos(beta)*cos(gamma));
        point_cloud[i][0] += x_shift;
        point_cloud[i][1] += y_shift;
        point_cloud[i][2] += z_shift;
    }
}

/* Return vector of extreme points from point cloud */
vector<vector<double>> get_n_extreme_points(int n, const vector<vector<double>>& point_cloud){
    vector<vector<double>> ext;
    set<int> added_indices;
    int nb = 0;
    while (nb<n) {
        double x_min_val = numeric_limits<double>::max(), y_min_val = numeric_limits<double>::max(), z_min_val = numeric_limits<double>::max();
        double x_max_val = numeric_limits<double>::lowest(), y_max_val = numeric_limits<double>::lowest(), z_max_val = numeric_limits<double>::lowest();
        vector<double> x_min(3), y_min(3), z_min(3);
        vector<double> x_max(3), y_max(3), z_max(3);
        size_t n = point_cloud.size();
        for (auto i = 0; i< n; i++) {
            if(added_indices.count(i)!=0)
                continue;
            if(point_cloud[i][0] < x_min_val){
                added_indices.insert(i);
                x_min_val = point_cloud[i][0];
                x_min[0] = point_cloud[i][0];
                x_min[1] = point_cloud[i][1];
                x_min[2] = point_cloud[i][2];
                if(nb<n){
                    ext.push_back(x_min);
                    nb++;
                }
                continue;
            }
            if(point_cloud[i][0] > x_max_val){
                added_indices.insert(i);
                x_max_val = point_cloud[i][0];
                x_max[0] = point_cloud[i][0];
                x_max[1] = point_cloud[i][1];
                x_max[2] = point_cloud[i][2];
                if(nb<n){
                    ext.push_back(x_max);
                    nb++;
                }
                continue;
            }
            if(point_cloud[i][1] < y_min_val){
                added_indices.insert(i);
                y_min_val = point_cloud[i][1];
                y_min[0] = point_cloud[i][0];
                y_min[1] = point_cloud[i][1];
                y_min[2] = point_cloud[i][2];
                if(nb<n){
                    ext.push_back(y_min);
                    nb++;
                }
                continue;
            }
            if(point_cloud[i][1] > y_max_val){
                added_indices.insert(i);
                y_max_val = point_cloud[i][1];
                y_max[0] = point_cloud[i][0];
                y_max[1] = point_cloud[i][1];
                y_max[2] = point_cloud[i][2];
                if(nb<n){
                    ext.push_back(y_max);
                    nb++;
                }
                continue;
            }
            if(point_cloud[i][2] < z_min_val){
                added_indices.insert(i);
                z_min_val = point_cloud[i][2];
                z_min[0] = point_cloud[i][0];
                z_min[1] = point_cloud[i][1];
                z_min[2] = point_cloud[i][2];
                if(nb<n){
                    ext.push_back(z_min);
                    nb++;
                }
                continue;
            }
            if(point_cloud[i][2] > z_max_val){
                added_indices.insert(i);
                z_max_val = point_cloud[i][2];
                z_max[0] = point_cloud[i][0];
                z_max[1] = point_cloud[i][1];
                z_max[2] = point_cloud[i][2];
                if(nb<n){
                    ext.push_back(z_max);
                    nb++;
                }
                continue;
            }
        }
    }
    return ext;
}


/* Return vector of extreme points from point cloud */
vector<vector<double>> get_extreme_points(const vector<vector<double>>& point_cloud){
    return get_n_extreme_points(6, point_cloud);
}

/* Return central point from point cloud */
vector<double> get_center(const vector<vector<double>>& point_cloud){
    
}


/* Compute the L2 error */
double computeL2error(const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data){
    size_t n = point_cloud_data.size();
    size_t m = point_cloud_model.size();
    double dist_sq = 0, err = 0;
    for (auto i = 0; i< n; i++) {
        double min_dist = numeric_limits<double>::max();
        for (auto j = 0; j< m; j++) {
            dist_sq = std::pow(point_cloud_data.at(i).at(0) - point_cloud_model.at(j).at(0),2) + std::pow(point_cloud_data.at(i).at(1) - point_cloud_model.at(j).at(1),2) + std::pow(point_cloud_data.at(i).at(2) - point_cloud_model.at(j).at(2),2);
            if(min_dist>dist_sq){
                min_dist = dist_sq;
            }
        }
        err += min_dist;
    }
    return err;
}


/* Plot two point clouds */
void plot(const vector<vector<double>>& ext_model, const vector<vector<double>>& ext_data, double point_thick){
    namespace plt = matplotlibcpp;
    vector<double> x_vec_model(ext_model.size()), y_vec_model(ext_model.size()), z_vec_model(ext_model.size());
    vector<double> x_vec_data(ext_data.size()), y_vec_data(ext_data.size()), z_vec_data(ext_data.size());
    for (int i = 0; i<ext_model.size(); i++) {
        x_vec_model[i] = ext_model[i][0];
        y_vec_model[i] = ext_model[i][1];
        z_vec_model[i] = ext_model[i][2];
    }
    for (int i = 0; i<ext_data.size(); i++) {
        x_vec_data[i] = ext_data[i][0];
        y_vec_data[i] = ext_data[i][1];
        z_vec_data[i] = ext_data[i][2];
    }
    std::map<std::string, std::string> keywords;
    keywords["marker"] = "s";
    keywords["linestyle"] = "None";
    keywords["ms"] = to_string(point_thick);
    plt::plot3(x_vec_model, y_vec_model, z_vec_model,x_vec_data, y_vec_data, z_vec_data, keywords);
    plt::show();
}


/* Run Go-ICP on point clouds */
tuple<double,double,double,double,double,double> run_GoICP(const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data){
    using namespace Go_ICP;
    
    int Nm = point_cloud_model.size(), Nd = point_cloud_data.size(), NdDownsampled = 1000;
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
    DebugOn("Roll = " << roll << endl);
    DebugOn("Pitch = " << pitch << endl);
    DebugOn("Yaw = " << yaw << endl);
    auto tx = goicp.optT.val[0][0];
    auto ty = goicp.optT.val[1][0];
    auto tz = goicp.optT.val[2][0];
    DebugOn("tx = " << tx << endl);
    DebugOn("ty = " << ty << endl);
    DebugOn("tz = " << tz << endl);
    return {roll,pitch,yaw,tx,ty,tz};
}

tuple<double,double,double,double,double,double> run_heuristic(const vector<vector<double>>& ext_model, vector<vector<double>>& ext_data, vector<vector<double>>& point_cloud_data){
    double roll = 0, pitch = 0, yaw = 0, x_shift = 0, y_shift = 0, z_shift = 1; /* at least one nonzero to enter the while loop */
    double final_roll = 0, final_pitch = 0, final_yaw = 0, final_x_shift = 0, final_y_shift = 0, final_z_shift = 0;
    int nb_iter = 0, max_nb_iter = 100;
    tuple<double,double,double,double,double,double> res;
    while(nb_iter < max_nb_iter && std::abs(roll)+std::abs(pitch)+std::abs(yaw)+std::abs(x_shift)+std::abs(y_shift)+std::abs(z_shift)>1e-1) {
        auto L2error = computeL2error(ext_model,ext_data);
        DebugOn("L2 error with exterme set before = " << L2error << endl);
        res = run_ARMO(false, "full", ext_model, ext_data);
        roll = get<0>(res);pitch = get<1>(res);yaw = get<2>(res);x_shift = get<3>(res);y_shift = get<4>(res);z_shift = get<5>(res);
        final_roll += roll;final_pitch += pitch;final_yaw += yaw;final_x_shift += x_shift;final_y_shift += y_shift;final_z_shift += z_shift;
        apply_rot_trans(roll, pitch, yaw, x_shift, y_shift, z_shift, ext_data);
        apply_rot_trans(roll, pitch, yaw, x_shift, y_shift, z_shift, point_cloud_data);
        L2error = computeL2error(ext_model,ext_data);
        DebugOn("L2 error with exterme set after full = " << L2error << endl);
        nb_iter++;
        DebugOn("ITERATION " << nb_iter << endl);
    }
    nb_iter = 0;z_shift=1;max_nb_iter = 100;
    while(nb_iter < max_nb_iter && std::abs(roll)+std::abs(pitch)+std::abs(yaw)+std::abs(x_shift)+std::abs(y_shift)+std::abs(z_shift)>1e-1) {
        res = run_ARMO(false, "z", ext_model, ext_data);
        roll = get<0>(res);pitch = get<1>(res);yaw = get<2>(res);x_shift = get<3>(res);y_shift = get<4>(res);z_shift = get<5>(res);
        final_roll += roll;final_pitch += pitch;final_yaw += yaw;final_x_shift += x_shift;final_y_shift += y_shift;final_z_shift += z_shift;
        apply_rot_trans(roll, pitch, yaw, x_shift, y_shift, z_shift, ext_data);
        apply_rot_trans(roll, pitch, yaw, x_shift, y_shift, z_shift, point_cloud_data);
        auto L2error = computeL2error(ext_model,ext_data);
        DebugOn("L2 error with exterme set after z = " << L2error << endl);
        nb_iter++;
        DebugOn("ITERATION " << nb_iter << endl);
    }
//    DebugOn("Plotting after z" << endl);
//    plot(ext_model,ext_data,1);
    nb_iter = 0;z_shift=1;max_nb_iter = 100;
    while(nb_iter < max_nb_iter && std::abs(roll)+std::abs(pitch)+std::abs(yaw)+std::abs(x_shift)+std::abs(y_shift)+std::abs(z_shift)>1e-1) {
        res = run_ARMO(false, "y", ext_model, ext_data);
        roll = get<0>(res);pitch = get<1>(res);yaw = get<2>(res);x_shift = get<3>(res);y_shift = get<4>(res);z_shift = get<5>(res);
        final_roll += roll;final_pitch += pitch;final_yaw += yaw;final_x_shift += x_shift;final_y_shift += y_shift;final_z_shift += z_shift;
        apply_rot_trans(roll, pitch, yaw, x_shift, y_shift, z_shift, ext_data);
        apply_rot_trans(roll, pitch, yaw, x_shift, y_shift, z_shift, point_cloud_data);
        auto L2error = computeL2error(ext_model,ext_data);
        DebugOn("L2 error with exterme set after y = " << L2error << endl);
        nb_iter++;
        DebugOn("ITERATION " << nb_iter << endl);
    }
//    DebugOn("Plotting after y" << endl);
//    plot(ext_model,ext_data,1);
    nb_iter = 0;z_shift=1;max_nb_iter = 100;
    while(nb_iter < max_nb_iter && std::abs(roll)+std::abs(pitch)+std::abs(yaw)+std::abs(x_shift)+std::abs(y_shift)+std::abs(z_shift)>1e-1) {
        res = run_ARMO(false, "x", ext_model, ext_data);
        roll = get<0>(res);pitch = get<1>(res);yaw = get<2>(res);x_shift = get<3>(res);y_shift = get<4>(res);z_shift = get<5>(res);
        final_roll += roll;final_pitch += pitch;final_yaw += yaw;final_x_shift += x_shift;final_y_shift += y_shift;final_z_shift += z_shift;
        apply_rot_trans(roll, pitch, yaw, x_shift, y_shift, z_shift, ext_data);
        apply_rot_trans(roll, pitch, yaw, x_shift, y_shift, z_shift, point_cloud_data);
        auto L2error = computeL2error(ext_model,ext_data);
        DebugOn("L2 error with exterme set after x = " << L2error << endl);
        nb_iter++;
        DebugOn("ITERATION " << nb_iter << endl);
    }
//    DebugOn("Plotting after x" << endl);
//    plot(ext_model,ext_data,1);
    nb_iter = 0;z_shift=1;max_nb_iter = 100;
    while(nb_iter < max_nb_iter && std::abs(roll)+std::abs(pitch)+std::abs(yaw)+std::abs(x_shift)+std::abs(y_shift)+std::abs(z_shift)>1e-1) {
        res = run_ARMO(false, "full", ext_model, ext_data);
        roll = get<0>(res);pitch = get<1>(res);yaw = get<2>(res);x_shift = get<3>(res);y_shift = get<4>(res);z_shift = get<5>(res);
        final_roll += roll;final_pitch += pitch;final_yaw += yaw;final_x_shift += x_shift;final_y_shift += y_shift;final_z_shift += z_shift;
        apply_rot_trans(roll, pitch, yaw, x_shift, y_shift, z_shift, ext_data);
        apply_rot_trans(roll, pitch, yaw, x_shift, y_shift, z_shift, point_cloud_data);
        auto L2error = computeL2error(ext_model,ext_data);
        DebugOn("L2 error with exterme set after full = " << L2error << endl);
        nb_iter++;
        DebugOn("ITERATION " << nb_iter << endl);
    }
    DebugOn("Final Roll (degrees) = " << final_roll << endl);
    DebugOn("Final Pitch (degrees) = " << final_pitch << endl);
    DebugOn("Final Yaw (degrees) = " << final_yaw << endl);
    DebugOn("Final Roll (radians)= " << final_roll *pi/180 << endl);
    DebugOn("Final Pitch (radians) = " << final_pitch *pi/180 << endl);
    DebugOn("Final Yaw (radians) = " << final_yaw *pi/180 << endl);
    DebugOn("Final x shift = " << final_x_shift << endl);
    DebugOn("Final y shift = " << final_y_shift << endl);
    DebugOn("Final z shift = " << final_z_shift << endl);
    return {final_roll,final_pitch,final_yaw,final_x_shift,final_y_shift,final_z_shift};
}
