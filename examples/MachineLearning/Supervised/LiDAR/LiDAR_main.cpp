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
#include "Lidar_utils.h"
#include "Heuristics.h"
#include "Goicp.h"
#include "icp.h"
#include "BB.h"
#include <algorithm>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <thread>
#ifdef USE_EIGEN3
#include <Eigen/Dense>
#include <Eigen/SVD>
#endif
#include <gravity/KDTreeVectorOfVectorsAdaptor.h>
#include <time.h>
using namespace std;
#include <gravity/jly_goicp.h>
#include <gravity/ConfigMap.hpp>
#include "gravity/nanoflann.hpp"
using namespace Go_ICP;
#ifdef USE_VORO
#include "voro++.hh"
using namespace voro;
#endif
#ifdef USE_GJK
extern "C" {
#include "openGJK.h"
}
#endif


/* Read input files */
void read_data(const rapidcsv::Document& doc,vector<vector<double>>& point_cloud, vector<vector<double>>& uav);

/* Return central point from point cloud */
vector<double> get_center(const vector<vector<double>>& point_cloud);

vector<pair<double,double>> center_point_cloud(vector<vector<double>>& point_cloud);


int main (int argc, char * argv[])
{
//    string file_u="/Users/smitha/Downloads/42.Kelley_MLOPT1_2021_09_03_COMBINED_BUSH.laz";
//    vector<vector<double>> uav_model, uav_data,uav_model1, uav_data1,cloud1, cloud2, uava1, uava2, rpy1, rpy2;
//    vector<vector<double>> full_rpy_model, full_rpy_data, rpy_model, rpy_data,rpy_model1, rpy_data1;
//    vector<vector<double>> lidar_point_cloud, roll_pitch_yaw, em;
//    auto uav_cloud_u=::read_laz1(file_u, lidar_point_cloud, roll_pitch_yaw);
//    int mid_i=0;
//    for(auto i=1;i<uav_cloud_u.size();i++)
//            {
//                auto x=uav_cloud_u.at(i)[0];
//                auto y=uav_cloud_u.at(i)[1];
//                auto z=uav_cloud_u.at(i)[2];
//                auto x_prev=uav_cloud_u.at(i-1)[0];
//                auto y_prev=uav_cloud_u.at(i-1)[1];
//                auto z_prev=uav_cloud_u.at(i-1)[2];
//                  if((abs(x-x_prev)>=1 && abs(y-y_prev)>=1)){
//    if(mid_i!=0){
//        DebugOn("More than Two flight lines are detected "<<mid_i<<" "<<endl);
//        DebugOn("Invalid flight line selection!!!!!!!!!!!!!!!!!!!!!!!!");
//    }
//    mid_i=i;
//    DebugOn("Two flight lines are detected "<<mid_i<<endl);
//
//}
//            }
//            /*If two flight lines identified2*/
//            if(mid_i==0){
//                invalid_argument("Two flight lines are not detected!");
//            }
//            for(auto i=0;i<mid_i;i++){
//                cloud1.push_back(lidar_point_cloud.at(i));
//                uava1.push_back(uav_cloud_u.at(i));
//                rpy1.push_back(roll_pitch_yaw.at(i));
//            }
//            for(auto i=mid_i;i<lidar_point_cloud.size();i++){
//                cloud2.push_back(lidar_point_cloud.at(i));
//                uava2.push_back(uav_cloud_u.at(i));
//                rpy2.push_back(roll_pitch_yaw.at(i));
//            }
//    vector<vector<double>> empty_vec;
//    empty_vec.push_back(uav_cloud_u[0]);
//#ifdef USE_MATPLOT
//     plot(uav_cloud_u,empty_vec, 0.1);
//#endif
//    empty_vec.clear();
//    empty_vec.push_back(lidar_point_cloud[0]);
//#ifdef USE_MATPLOT
//     plot(lidar_point_cloud,empty_vec, 0.1);
//#endif
    //    return 0
    double mse=1e-3;
    int dt=300;
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
        vector<vector<double>> point_cloud_model, point_cloud_data;
        string Model_file = string(prj_dir)+"/data_sets/LiDAR/toy_model.txt";
        string Data_file = string(prj_dir)+"/data_sets/LiDAR/toy_data.txt";
        string algo = "nsbb";
        if(argc>2){
            Model_file = argv[2];
        }
        if(argc>3){
            Data_file = argv[3];
        }
        if(argc>4){
            algo = argv[4];
        }
        bool run_goICP = (algo=="GoICP");
        bool run_gurobi = (algo=="GRB");
        if(run_goICP && argc==7){
            mse=stod(argv[5]);
            dt=stoi(argv[6]);
        }
        rapidcsv::Document  Model_doc(Model_file, rapidcsv::LabelParams(0, -1),rapidcsv::SeparatorParams('\t'));
        rapidcsv::Document  Data_doc(Data_file, rapidcsv::LabelParams(0, -1),rapidcsv::SeparatorParams('\t'));
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
        
        DebugOn("Model file has " << model_nb_rows << " rows" << endl);
        DebugOn("Data file has " << data_nb_rows << " rows" << endl);
        int row0 = 0;
        for (int i = row0; i< model_nb_rows; i++) { // Input iterator
            auto x = Model_doc.GetCell<double>(0, i);
            auto y = Model_doc.GetCell<double>(1, i);
            auto z = Model_doc.GetCell<double>(2, i);
            vector<double> xyz;
            xyz.push_back(x);
            xyz.push_back(y);
            xyz.push_back(z);
            point_cloud_model.push_back(xyz);
        }
        //point_cloud_data.resize(data_nb_rows);
        for (int i = row0; i< data_nb_rows; i++) { // Input iterator
            auto x = Data_doc.GetCell<double>(0, i);
            auto y = Data_doc.GetCell<double>(1, i);
            auto z = Data_doc.GetCell<double>(2, i);
            vector<double> xyz;
            xyz.push_back(x);
            xyz.push_back(y);
            xyz.push_back(z);
            point_cloud_data.push_back(xyz);
        }
       
        auto min_max_data=center_point_cloud(point_cloud_data);
        auto min_max_model=center_point_cloud(point_cloud_model);
       
        
        double best_ub=100;
        vector<double> best_rot_trans;
        //input
        double shift_min_x =  std::max(min_max_model[0].first,-0.5), shift_max_x = std::min(min_max_model[0].second,0.5), shift_min_y = std::max(min_max_model[1].first,-0.5),shift_max_y = std::min(min_max_model[1].second,0.5),shift_min_z = std::max(min_max_model[2].first, -0.5),shift_max_z = std::min(min_max_model[2].second,0.5);
        vector<pair<double, double>> min_max_d;
        double yaw_min = -90*pi/180., yaw_max = 90*pi/180., pitch_min =-90*pi/180.,pitch_max = 90*pi/180.,roll_min =-90*pi/180.,roll_max = 90*pi/180.;
        vector<vector<vector<double>>> model_voronoi_vertices(model_nb_rows);/* Store the normal vector of each facet of the voronoi cell of each point */
        
    PointCloud<double> cloud;
        
        cloud.pts.resize(point_cloud_model.size());
        for (size_t i = 0; i < point_cloud_model.size(); i++)
        {
            cloud.pts[i].x = point_cloud_model[i][0];
            cloud.pts[i].y = point_cloud_model[i][1];
            cloud.pts[i].z = point_cloud_model[i][2];
        }
        nanoflann::KDTreeSingleIndexAdaptor<
            nanoflann::L2_Simple_Adaptor<double, PointCloud<double>>,
            PointCloud<double>, 3 /* dim */> index(3 /*dim*/, cloud, {10 /* max leaf */});
        index.buildIndex();
        //vector<double> rpyt_H;
        auto rpyt_H=ub_heuristic_disc(index, point_cloud_model, point_cloud_data, roll_min, roll_max, pitch_min, pitch_max, yaw_min, yaw_max,shift_min_x, shift_max_x, shift_min_y, shift_max_y, shift_min_z, shift_max_z, best_ub, "L2", 100);
  
        double ub_sq_root=sqrt(best_ub);
        
        
        map<int, vector<int>> incomp;
        
        double dmax=0;
        for(auto i=0;i<point_cloud_data.size()-1;i++){
            vector<int> red;
            for(auto j=i+1;j<point_cloud_data.size();j++){
                auto d=pow(point_cloud_data.at(i)[0]-point_cloud_data.at(j)[0],2)+
                pow(point_cloud_data.at(i)[1]-point_cloud_data.at(j)[1],2)+
                pow(point_cloud_data.at(i)[2]-point_cloud_data.at(j)[2],2);
                if(d>=3*best_ub+1e-9){
                    red.push_back(j);
                    DebugOff("cannot");
                }
                if(d>=dmax){
                    dmax=d;
                }
            }
            if(red.size()>=1){
                incomp[i]=red;
            }
        }
        
        
        map<int, vector<int>> incompj;
        map<int, vector<int>> incomp_dmax;
        
        for(auto i=0;i<point_cloud_model.size()-1;i++){
            vector<int> red;
            for(auto j=i+1;j<point_cloud_model.size();j++){
                auto d=pow(point_cloud_model.at(i)[0]-point_cloud_model.at(j)[0],2)+
                pow(point_cloud_model.at(i)[1]-point_cloud_model.at(j)[1],2)+
                pow(point_cloud_model.at(i)[2]-point_cloud_model.at(j)[2],2);
                if(d>=3*best_ub+1e-9){
                    //red.push_back(j);
                    //DebugOn("nonmatch"<<endl);
                }
                if(d>=pow(sqrt(3*best_ub)+sqrt(dmax),2)){
                    red.push_back(j);
                    DebugOff("nonmatch"<<endl);
                }
            }
            if(red.size()>=1){
                incompj[i]=red;
            }
        }
        DebugOff("screened");
        
        param<double> dist_jj("dist_jj");
        param<double> dist_ii("dist_ii");
        for(auto j=0;j<point_cloud_model.size()-1;j++){
            for(auto i=j+1;i<point_cloud_model.size();i++){
                auto d=pow(point_cloud_model.at(i)[0]-point_cloud_model.at(j)[0],2)+
                pow(point_cloud_model.at(i)[1]-point_cloud_model.at(j)[1],2)+
                pow(point_cloud_model.at(i)[2]-point_cloud_model.at(j)[2],2);
                auto d_sq=sqrt(d);
                dist_jj.add_val(to_string(i+1)+","+to_string(j+1), d_sq);
                dist_jj.add_val(to_string(j+1)+","+to_string(i+1), d_sq);
            }
        }
        for(auto j=0;j<point_cloud_data.size()-1;j++){
            for(auto i=j+1;i<point_cloud_data.size();i++){
                auto d=pow(point_cloud_data.at(i)[0]-point_cloud_data.at(j)[0],2)+
                pow(point_cloud_data.at(i)[1]-point_cloud_data.at(j)[1],2)+
                pow(point_cloud_data.at(i)[2]-point_cloud_data.at(j)[2],2);
                auto d_sq=sqrt(d);
                dist_ii.add_val(to_string(i+1)+","+to_string(j+1), d_sq);
                dist_ii.add_val(to_string(j+1)+","+to_string(i+1), d_sq);
            }
        }
        

#ifdef USE_VORO
        //container model_con(min_max_d[0].first-1e-4,min_max_d[0].second+1e-4,min_max_d[1].first-1e-4,min_max_d[1].second+1e-4,min_max_d[2].first-1e-4,min_max_d[2].second+1e-4,10,10,10,false,false,false,8);
        container model_con(-1,1,-1,1,-1,1,10,10,10,false,false,false,8);
        //container model_con(min_max_d[0].first,min_max_d[0].second,min_max_d[1].first,min_max_d[1].second,min_max_d[2].first,min_max_d[2].second,10,10,10,false,false,false,8);
        //container model_con(x_min,x_max,y_min, y_max, z_min, z_max,10,10,10,false,false,false,8);
        for (int i = row0; i< model_nb_rows; i++) { // Input iterator
            model_con.put(i, point_cloud_model[i][0], point_cloud_model[i][1], point_cloud_model[i][2]);
        }
        /* Compute the facets of the Voronoi cells of all model points */
        c_loop_all cl(model_con);
        int idx,nx,ny,nz;
        double x,y,z,x1,y1,z1;
        voronoicell c;
        vector<int> face_vert;
        vector<double> v;
        vector<vector<vector<double>>> model_voronoi_normals(model_nb_rows);/* Store the normal vector of each facet of the voronoi cell of each point */
        vector<vector<pair<double, double>>> model_voronoi_min_max(model_nb_rows);/* Store the normal vector of each facet of the voronoi cell of each point */
        vector<vector<vector<double>>> model_face_pts(model_nb_rows);/* Store a point from each facet of the voronoi cell of each point */
        vector<vector<vector<int>>> model_voronoi_vertex_edge(model_nb_rows);/* For each cell, for each vertex, store list of other vertices it is linked to */
        vector<vector<vector<pair<int, int>>>> model_voronoi_vertex_edge_planes(model_nb_rows);/* For each vertex, corresponding to each edge containing the vertex (Vertex_edge), the two planes that make the edge are given*/
        vector<vector<vector<int>>> facet_vertex(model_nb_rows);/* For each cell, for each facet, store list of vertices on the facet */
        vector<vector<double>> model_face_intercept(model_nb_rows);/* Store the constant part (intercept) in the equation of the voronoi cell of each point */
        vector<double> model_voronoi_in_radius(model_nb_rows);/* Store the radius of the largest ball contained IN the voronoi cell of each point */
        vector<double> model_voronoi_max_seg_dist_sq(model_nb_rows);/* Store the radius of the largest ball contained IN the voronoi cell of each point */
        vector<double> model_voronoi_out_radius(model_nb_rows);/* Store the radius of the smallest ball enclosing the voronoi cell of each point */
        vector<double> model_voronoi_out_radius_sq(model_nb_rows);/* Store the radius of the smallest ball enclosing the voronoi cell of each point */
        vector<double> model_inner_prod_max(model_nb_rows);
        vector<double> model_inner_prod_min(model_nb_rows);
        vector<double> max_vert_vert_dist_sq(model_nb_rows);
        param<> norm_x("norm_x"), norm_y("norm_y"), norm_z("norm_z"), intercept("intercept");
        param<> model_radius("model_radius");
        indices m_facets("m_facets");
        vector<pair<double,int>> volume(model_nb_rows);/*volume of each voronoi cell <volume,point_id>*/
        int total_nb_faces = 0;
        if(cl.start()) do if(model_con.compute_cell(c,cl)) {
            cl.pos(x,y,z);
            idx=cl.pid();
            c.vertices(x,y,z,v);
            auto ce=c.ed;
            auto cnu=c.nu;
            
            int nb_vertices = v.size()/3;
            vector<vector<int>> ve(nb_vertices);
            DebugOn("cell "<<idx<<"with "<<nb_vertices<<" vertices"<<endl);
            for (auto k=0;k<nb_vertices;k++){
                for(auto l=0;l<cnu[k];l++){
                    DebugOff(ce[k][l]<<" ");
                    ve[k].push_back(ce[k][l]);
                }
                DebugOff(endl);
            }
            model_voronoi_vertex_edge[idx]=ve;
            volume[idx] = {c.volume(),idx};
            int v_idx = 0;
            double max_dist = 0;
            model_voronoi_vertices[idx].resize(nb_vertices);
            vector<double> vertices;
            vertices.resize(3);
            vector<int> v_order;
            Debug("Cell has " << c.number_of_edges() << " edges \n");
            double model_ip_max=-999.0;
            double model_ip_min=1000.0;
            double xmin=100, xmax=-100, ymin=100, ymax=-100, zmin=100, zmax=-100;
            for (int i = 0; i<nb_vertices; i++) {
                double dist_sq = std::pow(x - v[v_idx],2) + std::pow(y - v[v_idx+1],2) + std::pow(z - v[v_idx+2],2);
                if(dist_sq>max_dist)
                    max_dist = dist_sq;
                vertices[0]=v[v_idx];
                vertices[1]=v[v_idx+1];
                vertices[2]=v[v_idx+2];
                auto ip=vertices[0]*x+vertices[1]*y+vertices[2]*z;
                if(ip>=model_ip_max){
                    model_ip_max=ip;
                }
                if(ip<=model_ip_min){
                    model_ip_min=ip;
                }
                if(vertices[0]<=xmin){
                    xmin=vertices[0];
                }
                if(vertices[0]>=xmax){
                    xmax=vertices[0];
                }
                if(vertices[1]<=ymin){
                    ymin=vertices[1];
                }
                if(vertices[1]>=ymax){
                    ymax=vertices[1];
                }
                if(vertices[2]<=zmin){
                    zmin=vertices[2];
                }
                if(vertices[2]>=zmax){
                    zmax=vertices[2];
                }
                model_voronoi_vertices[idx][i]=vertices;
                v_idx += 3;
            }
            double dist_max_sq=0;
            for (int i = 0; i<nb_vertices; i++) {
                auto vi=model_voronoi_vertices[idx][i];
                for (int j = i+1; j<nb_vertices; j++) {
                    auto vj=model_voronoi_vertices[idx][j];
                    double d=pow(vi[0]-vj[0],2)+pow(vi[1]-vj[1],2)+pow(vi[2]-vj[2],2);
                    if(d>=dist_max_sq){
                        dist_max_sq=d;
                    }
                }
            }
            max_vert_vert_dist_sq[idx]=dist_max_sq;
            vector<pair<double, double>> min_max;
            min_max.resize(3);
            min_max[0].first=xmin;
            min_max[0].second=xmax;
            min_max[1].first=ymin;
            min_max[1].second=ymax;
            min_max[2].first=zmin;
            min_max[2].second=zmax;
            model_voronoi_max_seg_dist_sq[idx]=pow(xmax-xmin, 2)+pow(ymax-ymin, 2)+pow(zmax-zmin, 2);
            model_voronoi_min_max[idx]=min_max;
            model_inner_prod_max[idx]=model_ip_max;
            model_inner_prod_min[idx]=model_ip_min;
            model_voronoi_out_radius[idx] = std::sqrt(max_dist);/* the radius of the smallest ball enclosing the voronoi cell is the distance from the center to the farthest vertex */
            model_voronoi_out_radius_sq[idx] = (max_dist);
            model_radius.add_val(to_string(idx+1), model_voronoi_out_radius[idx]);
            vector<double> normals;
            c.normals(normals);
            int nb_normals = normals.size()/3;
            model_voronoi_normals[idx].resize(nb_normals);
            int normal_idx = 0;
            for (int i = 0; i<nb_normals; i++) {
                model_voronoi_normals[idx][i] = {normals[normal_idx],normals[normal_idx+1],normals[normal_idx+2]};
                normal_idx += 3;
            }
            c.face_vertices(face_vert);
            
            int nb_faces = c.number_of_faces();
            total_nb_faces += nb_faces;
            DebugOn("The Voronoi cell of point : (" << x << "," << y << ","<< z << ") has " << nb_faces << " facets.\n");
            model_face_pts[idx].resize(nb_faces);
            model_face_intercept[idx].resize(nb_faces);
            int face_id = 0;
            double min_dist = numeric_limits<double>::max();
            for (int f_id = 0; f_id<nb_faces; f_id++) {
                int nb_v = face_vert[face_id];
                Debug("Facet " << f_id << " has " << nb_v << " vertices.\n");
                face_id++;
                x1 = v[3*face_vert[face_id]];y1 = v[3*face_vert[face_id]+1]; z1 = v[3*face_vert[face_id]+2];
                face_id += nb_v;
                model_face_pts[idx][f_id] = {x1,y1,z1};
                /* Getting the plane equation of each face in standard form: ax + by + cz +d = 0 */
                /* The plane's equation is given by: n_x(x - x1) + n_y(y - y1) + n_z(z  - z1) = 0 */
                double a = model_voronoi_normals[idx][f_id][0];/* n_x */
                double b = model_voronoi_normals[idx][f_id][1];/* n_y */
                double c = model_voronoi_normals[idx][f_id][2];/* n_z */
                double d = -1*a*x1 - b*y1 - c*z1; /* -n_x*x1 - n_y*y1 - n_z*z1 */
                model_face_intercept[idx][f_id] = d;
                /* check which side of the plane (x,y,z) is on */
                double eq = a*x + b*y + c*z + d;
                /* make sure that the model point satisfies eq < 0, i.e. points insides the voronoi cell have to satisfy eq <= 0 */
                if (eq==0) {
                    throw invalid_argument("model point cannot be on the voronoi face!");
                }
                if (eq>0) {
                    model_voronoi_normals[idx][f_id][0] *= -1;
                    model_voronoi_normals[idx][f_id][1] *= -1;
                    model_voronoi_normals[idx][f_id][2] *= -1;
                    model_face_intercept[idx][f_id] *= -1;
                }
                string key = to_string(idx+1)+","+to_string(f_id+1);
                norm_x.add_val(key, model_voronoi_normals[idx][f_id][0]);
                norm_y.add_val(key, model_voronoi_normals[idx][f_id][1]);
                norm_z.add_val(key, model_voronoi_normals[idx][f_id][2]);
                intercept.add_val(key, model_face_intercept[idx][f_id]);
                m_facets.insert(key);
                Debug("for key: " << key << endl);
                Debug("one point on this facet: (" << x1 << "," << y1 << ","<< z1 << ").\n");
                Debug("Facet " << f_id << " equation: " << model_voronoi_normals[idx][f_id][0] << "x + " << model_voronoi_normals[idx][f_id][1] << "y + " << model_voronoi_normals[idx][f_id][2] << "z + " << model_face_intercept[idx][f_id] << " = 0\n");
                
                //                double dist = std::abs(data_voronoi_normals[idx][f_id][0]*x + data_voronoi_normals[idx][f_id][1]*y + data_voronoi_normals[idx][f_id][2]*z + data_voronoi_normals[idx][f_id][0]*x1 + data_voronoi_normals[idx][f_id][1]*y1 + data_voronoi_normals[idx][f_id][2]*z1)/std::sqrt(std::pow(data_voronoi_normals[idx][f_id][0],2) + std::pow(data_voronoi_normals[idx][f_id][1],2) + std::pow(data_voronoi_normals[idx][f_id][2],2));
                //                if(dist<min_dist){
                //                    min_dist = dist;
                //                }
            }
            model_voronoi_in_radius[idx] = min_dist;/* the radius of the largest ball in the voronoi cell is the distance from the center to the nearest facet */
            face_id=0;
            /* fv stores for each plane, all the vertices that are on the plane*/
            vector<vector<int>> fv(nb_faces);
            for (int f_id = 0; f_id<nb_faces; f_id++) {
                int nb_v = face_vert[face_id];
                DebugOn("f_id "<<f_id<<"\t");
                for(auto k=1;k<=nb_v;k++){
                    ++face_id;
                    fv[f_id].push_back(face_vert[face_id]);
                    DebugOn(face_vert[face_id]<<" ");
                }
                face_id++;
                DebugOn(endl);
            }
            facet_vertex[idx]=fv;
            /* fv stores for each plane, all the vertices that are on the plane*/
            vector<vector<pair<int, int>>> vep(nb_vertices);
            for(auto i=0;i<nb_vertices;i++){
                for(auto j=0;j<ve[i].size();j++){
                    vector<int> p;
                    int vj=ve[i][j];/*vj loops over all vertices that are connected to i by edge*/
                    for(auto k=0;k<fv.size();k++){
                        if(std::find(fv[k].begin(),fv[k].end(),i)!=fv[k].end() && std::find(fv[k].begin(),fv[k].end(),vj)!=fv[k].end()){
                            p.push_back(k);
                        }
                        if(p.size()==2){ /* Each edge is defined by  two planes*/
                            break;
                        }
                    }
                    if(p.size()!=2){
                        throw invalid_argument("no two planes found defining an edge");
                    }
                    pair<int,int> pair_p={p[0], p[1]};
                    vep[i].push_back(pair_p);/* for each vertices i and vj, for each edge, ve[i][j], vep[i][j] shows pair of planes defining that edge */
                }
            }
            model_voronoi_vertex_edge_planes[idx]=vep;
            v.clear();
            face_vert.clear();
        } while (cl.inc());
#endif

        if(run_goICP){/* Run GoICP inline */
            vector<int> matching(point_cloud_data.size());
            vector<double> err_per_point(point_cloud_data.size());
//            auto res_icp = run_GoICP(point_cloud_model, point_cloud_data, mse, dt);
//            auto roll = get<0>(res_icp);auto pitch = get<1>(res_icp);auto yaw = get<2>(res_icp);auto x_shift = get<3>(res_icp);auto y_shift = get<4>(res_icp);auto z_shift = get<5>(res_icp);
//            double roll=0.161051;
//            double pitch=1.72447;
//            double yaw=-1.42412;
//            double x_shift=-0.128593;
//            double y_shift=-0.143401;
//            double z_shift=0.118928;
            double roll=-0.0949683;
            double pitch=0.0425285;
            double yaw=-0.199007;
            double x_shift=0.218107;
            double y_shift=-0.250144;
            double z_shift=0.0160762;
//            double roll=-1.04111;
//            double pitch=1.19823;
//            double yaw=0.558101;
//            double x_shift=0.0568501;
//            double y_shift=-0.129687;
//            double z_shift=0.0191642;
            apply_rot_trans(roll, pitch, yaw, x_shift, y_shift, z_shift, point_cloud_data);
            //plot(point_cloud_model,point_cloud_data,1);
            auto err=computeL2error(point_cloud_model, point_cloud_data, matching, err_per_point);
            DebugOn("Go ICP error "<<err<<endl);
        }
        else if(run_gurobi){
            double error=0;
            indices valid_cells_old("valid_cells_old");
            indices new_cells("new_cells");
            param<double> dist_cells("dist_cells");
            double prep_time=0;
            double min_cost_sum=0;
            preprocess_lid(point_cloud_model, point_cloud_data, model_voronoi_vertices, valid_cells_old,new_cells, dist_cells, roll_min, roll_max, pitch_min, pitch_max, yaw_min, yaw_max,shift_min_x, shift_max_x, shift_min_y, shift_max_y, shift_min_z, shift_max_z,  best_ub, prep_time, min_cost_sum, "L2",  dist_ii, dist_jj, model_voronoi_out_radius_sq);
            map<int, vector<int>> incomp;
            auto model=Reg_L2_model_rotation_trigonometric(point_cloud_model, point_cloud_data, roll_min, roll_max, pitch_min, pitch_max, yaw_min, yaw_max,shift_min_x, shift_max_x, shift_min_y, shift_max_y, shift_min_z, shift_max_z, new_cells, dist_cells, incomp);
            model->print();
            solver<> S(model,gurobi);
            //            S.use_callback();
            //           // R->replace_integers();
        S.run(5,1e-6, 1e9, 1, false, {false,""}, 172700,best_ub, 72);
//#ifdef USE_MATPLOT
//            apply_rot_trans(rpyt[0], rpyt[1], rpyt[2], rpyt[3], rpyt[4], rpyt[5], point_cloud_data);
//            //plot(point_cloud_model,point_cloud_data,1);
//#endif
        }
        else{
            double error=0;
            indices valid_cells_old("valid_cells_old");
            indices new_cells("new_cells");
            param<double> dist_cells("dist_cells");
            double prep_time=0;
            double min_cost_sum=0;
            preprocess_lid(point_cloud_model, point_cloud_data, model_voronoi_vertices, valid_cells_old,new_cells, dist_cells, roll_min, roll_max, pitch_min, pitch_max, yaw_min, yaw_max,shift_min_x, shift_max_x, shift_min_y, shift_max_y, shift_min_z, shift_max_z,  best_ub, prep_time, min_cost_sum, "L2", dist_ii, dist_jj, model_voronoi_out_radius_sq);
            vector<double> rpyt;
#ifdef USE_MATPLOT
            rpyt=BranchBound_Align(point_cloud_model, point_cloud_data, model_voronoi_vertices, rpyt_H, best_ub, shift_min_x, shift_max_x, shift_min_y, shift_max_y, shift_min_z, shift_max_z, "L2",  index, dist_ii, dist_jj, model_voronoi_out_radius_sq);
#else
            rpyt=BranchBound_MPI(point_cloud_model, point_cloud_data, model_voronoi_vertices, rpyt_H, best_ub, shift_min_x, shift_max_x, shift_min_y, shift_max_y, shift_min_z, shift_max_z, "L2",  index, dist_ii, dist_jj, model_voronoi_out_radius_sq);
#endif
#ifdef USE_MATPLOT
            apply_rot_trans(rpyt[0], rpyt[1], rpyt[2], rpyt[3], rpyt[4], rpyt[5], point_cloud_data);
            //plot(point_cloud_model,point_cloud_data,1);
#endif
            
          
        }
    return 0;
}


/* Return central point from point cloud */
vector<double> get_center(const vector<vector<double>>& point_cloud){
    int n=point_cloud.size();
    vector<double> res(3);
    double cx=0, cy=0, cz=0;
    for(auto i=0;i<n;i++)
    {
        auto x=point_cloud.at(i)[0];
        auto y=point_cloud.at(i)[1];
        auto z=point_cloud.at(i)[2];
        cx+=x;
        cy+=y;
        cz+=z;
    }
    cx=cx/n;
    cy=cy/n;
    cz=cz/n;
    res[0]=cx;
    res[1]=cy;
    res[2]=cz;
    return(res);
}
/*Assumption all points are in between -1,1*/


vector<pair<double,double>> center_point_cloud(vector<vector<double>>& point_cloud){
    double cx=0,cy=0,cz=0;
    int n=point_cloud.size();
    double min_x=10, max_x=-1, min_y=10, max_y=-1,min_z=10, max_z=-1;
    vector<pair<double,double>> res;
    for(auto i=0;i<n;i++)
    {
        auto x=point_cloud.at(i)[0];
        auto y=point_cloud.at(i)[1];
        auto z=point_cloud.at(i)[2];
        cx+=x;
        cy+=y;
        cz+=z;
    }
    cx=cx/n;
    cy=cy/n;
    cz=cz/n;
    double dmin=100, dmax=-1,d;
    for(auto i=0;i<point_cloud.size();i++)
    {
        point_cloud.at(i)[0]-=cx;
        point_cloud.at(i)[1]-=cy;
        point_cloud.at(i)[2]-=cz;
        if(point_cloud.at(i)[0]<=min_x){
            min_x=point_cloud.at(i)[0];
        }
        if(point_cloud.at(i)[0]>=max_x){
            max_x=point_cloud.at(i)[0];
        }
        if(point_cloud.at(i)[1]<=min_y){
            min_y=point_cloud.at(i)[1];
        }
        if(point_cloud.at(i)[1]>=max_y){
            max_y=point_cloud.at(i)[1];
        }
        if(point_cloud.at(i)[2]<=min_z){
            min_z=point_cloud.at(i)[2];
        }
        if(point_cloud.at(i)[2]>=max_z){
            max_z=point_cloud.at(i)[2];
        }
        d=pow(point_cloud.at(i)[0],2)+pow(point_cloud.at(i)[1],2)+pow(point_cloud.at(i)[2],2);
        if(d<=dmin){
            dmin=d;
        }
        if(d>=dmax){
            dmax=d;
        }
        
    }
    res.push_back(make_pair(min_x, max_x));
    res.push_back(make_pair(min_y, max_y));
    res.push_back(make_pair(min_z, max_z));
    res.push_back(make_pair(dmin, dmax));
    res.push_back(make_pair(sqrt(dmin), sqrt(dmax)));
    return res;
}

#ifdef USE_MATPLOT
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
    //    keywords["aspect"] = "equal";
    
    plt::plot3(x_vec_model, y_vec_model, z_vec_model,x_vec_data, y_vec_data, z_vec_data, keywords);
    
    plt::show();
}

void plot(const vector<vector<double>>& ext_model, const vector<vector<double>>& ext_data, const vector<vector<double>>& ext_data1, double point_thick)
{
    namespace plt = matplotlibcpp;
    vector<double> x_vec_model(ext_model.size()), y_vec_model(ext_model.size()), z_vec_model(ext_model.size());
    vector<double> x_vec_data(ext_data.size()), y_vec_data(ext_data.size()), z_vec_data(ext_data.size());
    vector<double> a_vec_data(ext_data1.size()), b_vec_data(ext_data1.size()), c_vec_data(ext_data1.size());
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
    for (int i = 0; i<ext_data1.size(); i++) {
        a_vec_data[i] = ext_data1[i][0];
        b_vec_data[i] = ext_data1[i][1];
        c_vec_data[i] = ext_data1[i][2];
    }
    DebugOn(a_vec_data.size()<<" "<<b_vec_data.size()<<" "<<c_vec_data.size());
    std::map<std::string, std::string> keywords;
    keywords["marker"] = "s";
    keywords["linestyle"] = "None";
    keywords["ms"] = to_string(point_thick);
    //    keywords["aspect"] = "equal";
    
    plt::plot3(x_vec_model, y_vec_model, z_vec_model,x_vec_data, y_vec_data, z_vec_data,a_vec_data, b_vec_data, c_vec_data, keywords);
    
    plt::show();
}
void plot(const vector<vector<double>>& ext_model, const vector<vector<double>>& ext_data, const vector<vector<double>>& ext_data1,const vector<vector<double>>& ext_data2, double point_thick)
{
    namespace plt = matplotlibcpp;
    vector<double> x_vec_model(ext_model.size()), y_vec_model(ext_model.size()), z_vec_model(ext_model.size());
    vector<double> x_vec_data(ext_data.size()), y_vec_data(ext_data.size()), z_vec_data(ext_data.size());
    vector<double> a_vec_data(ext_data1.size()), b_vec_data(ext_data1.size()), c_vec_data(ext_data1.size());
    vector<double> e_vec_data(ext_data2.size()), f_vec_data(ext_data2.size()), g_vec_data(ext_data2.size());
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
    for (int i = 0; i<ext_data1.size(); i++) {
        a_vec_data[i] = ext_data1[i][0];
        b_vec_data[i] = ext_data1[i][1];
        c_vec_data[i] = ext_data1[i][2];
    }
    for (int i = 0; i<ext_data2.size(); i++) {
        e_vec_data[i] = ext_data2[i][0];
        f_vec_data[i] = ext_data2[i][1];
        g_vec_data[i] = ext_data2[i][2];
    }
    DebugOn(a_vec_data.size()<<" "<<b_vec_data.size()<<" "<<c_vec_data.size());
    std::map<std::string, std::string> keywords;
    keywords["marker"] = "s";
    keywords["linestyle"] = "None";
    keywords["ms"] = to_string(point_thick);
    //    keywords["aspect"] = "equal";
    
    plt::plot3(x_vec_model, y_vec_model, z_vec_model,x_vec_data, y_vec_data, z_vec_data,a_vec_data, b_vec_data, c_vec_data,e_vec_data, f_vec_data, g_vec_data, keywords);
    
    plt::show();
}
#endif
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
