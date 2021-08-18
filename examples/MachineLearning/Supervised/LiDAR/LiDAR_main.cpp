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
void read_data(const rapidcsv::Document& doc,vector<vector<double>>& point_cloud, vector<vector<double>>& uav);

/* Read LAZ files */
void read_laz(const string& fname);

/* Save LAZ files */
void save_laz(const string& fname, const vector<vector<double>>& point_cloud1, const vector<vector<double>>& point_cloud2);

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
vector<pair<double,double>> get_min_max(double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, const vector<double>& p, const vector<double>& ref);

double get_GoICP_dist(double radius_r, double radius_t, const vector<double>& p, bool L1norm);

double get_max_dist(double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double tx_min, double tx_max, double ty_min, double ty_max, double tz_min, double tz_max, const vector<double>& p, const vector<double>& ref, bool L1norm = false);

bool compute_vertices(vector<vector<double>> vertex_set_a, vector<vector<double>> facets_a, vector<vector<double>> facets_b, vector<vector<int>> vertex_edge_a, vector<vector<pair<int, int>>> vertex_edge_plane_a,   vector<int> feas_set_a, vector<int> infeas_set_a, vector<vector<int>> infeas_facet_set_a, std::vector<std::vector<double>>& new_vert, const std::vector<double>& model_point);

/* Return true if two cubes intersect
 The cube is stored using a vector of size 3: {x,y,z}, where each entry is [min,max] on the corresponding axis
 */
bool intersect(const vector<pair<double,double>>& a, const vector<pair<double,double>>& b);

/* Returns the coordinates of the cube center
 The cube is stored using a vector of size 3: {x,y,z}, where each entry is [min,max] on the corresponding axis
 */
tuple<double,double,double> get_center(const vector<pair<double,double>>& cube);

/* Apply rotation + translation on input data (using rotation +translation matrix) */
void apply_rot_trans(const vector<double>& theta_matrix, vector<vector<double>>& point_cloud);

/* Apply rotation + translation on input data (using 3 angles) */
void apply_rot_trans(double roll, double pitch, double yaw, double x_shift, double y_shift, double z_shift, vector<vector<double>>& point_cloud);

/* Return vector of 6 extreme points (2 per axis) from point cloud */
vector<vector<double>> get_extreme_points(const vector<vector<double>>& point_cloud);


/* Return vector of n extreme points from point cloud */
vector<vector<double>> get_n_extreme_points(int n, const vector<vector<double>>& point_cloud);


/* Return vector of n extreme points from point cloud */
pair<vector<vector<double>>,vector<vector<double>>> get_n_extreme_points(int n, const vector<vector<double>>& point_cloud, const vector<vector<double>>& uav);


/* Run the iterative projection heuristic on the Registration problem */
tuple<double,double,double,double,double,double> run_IPH(const vector<vector<double>>& ext_model, vector<vector<double>>& ext_data, vector<vector<double>>& point_cloud_data);

/* Run the iterative projection heuristic on the Boresight Alignement problem */
tuple<double,double,double> run_IPH(vector<vector<double>>& ext_model, vector<vector<double>>& ext_data, const vector<vector<double>>& uav1, const vector<vector<double>>& uav2);

/* Compute the L2 error */
double computeL2error(const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data, vector<int>& matching, vector<double>& err_per_point);

double computeL2error(const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data, indices valid_cells, vector<int>& matching, vector<double>& err_per_point);

/* Compute the L1 error */
double computeL1error(const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data, vector<int>& matching, vector<double>& err_per_point);

/* Return central point from point cloud */
vector<double> get_center(const vector<vector<double>>& point_cloud);


/* Apply rotation on input data (Boresight alignment) */
void apply_rotation(double roll, double pitch, double yaw, vector<vector<double>>& point_cloud1, vector<vector<double>>& point_cloud2, const vector<vector<double>>& uav1, const vector<vector<double>>& uav2);

/* Update point cloud coordinates */
void update_xyz(vector<vector<double>>& point_cloud1, vector<double>& x_vec1, vector<double>& y_vec1, vector<double>& z_vec1);
void preprocess_poltyope_ve_gjk_in_centroid(const vector<vector<double>>& point_cloud_data, const vector<vector<double>>& point_cloud_model, const indices& old_cells, param<double>& dist_cost_old, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double shift_min_x, double shift_max_x, double shift_min_y, double shift_max_y, double shift_min_z, double shift_max_z, const vector<vector<vector<double>>>& model_voronoi_normals, const vector<vector<double>>& model_face_intercept, const vector<vector<vector<double>>>& model_voronoi_vertices, param<double>& dist_cost, double upper_bound, double lower_bound, double& min_cost_sum, indices& new_cells,  double& new_tx_min, double& new_tx_max, double& new_ty_min, double& new_ty_max, double& new_tz_min, double& new_tz_max, double& prep_time_total, const vector<vector<pair<double, double>>>& model_voronoi_min_max, const vector<vector<vector<int>>>& model_voronoi_vertex_edge, const vector<vector<vector<pair<int,int>>>>& model_voronoi_vertex_edge_planes, const vector<double> max_vert_vert_dist_sq, const param<double>& dii, const param<double>& djj, param<double>& maxj);/* Run the ARMO model for boresight alignment */
void run_preprocess_parallel_new(const vector<vector<double>>& point_cloud_data, const vector<vector<double>>& point_cloud_model, const vector<vector<vector<double>>>& model_voronoi_normals, const vector<vector<double>>& model_face_intercept, const vector<vector<vector<double>>>& model_voronoi_vertices, vector<int>& pos_vec, vector<shared_ptr<Model<double>>>& models, const vector<treenode_p>& vec_node, vector<int>& m_vec,vector<double>& vec_lb, vector<indices>& valid_cells, int nb_threads, double upper_bound, double lower_bound, vector<double>& new_shift_x_min, vector<double>& new_shift_x_max, vector<double>& new_shift_y_min, vector<double>& new_shift_y_max, vector<double>& new_shift_z_min, vector<double>& new_shift_z_max, int max_cells, const vector<vector<pair<double, double>>>& model_voronoi_min_max, const vector<vector<vector<int>>>& model_voronoi_vertex_edge, const vector<vector<vector<pair<int,int>>>>& model_voronoi_vertex_edge_planes, const param<double>& dii, const param<double>& djj, vector<param<double>>& dist_cost_new, int iter, vector<vector<double>>& costs_upto_vec, const vector<double>& max_vert_vert_dist_sq);
void run_preprocess_model(const vector<vector<double>>& point_cloud_data, const vector<vector<double>>& point_cloud_model, const vector<vector<vector<double>>>& model_voronoi_normals, const vector<vector<double>>& model_face_intercept, const vector<vector<vector<double>>>& model_voronoi_vertices, shared_ptr<Model<double>>& model_i, int& m_vec_i, int& pos_vec_i, double& vec_lb_i, treenode_p vec_node_i, indices& valid_cells_i, param<double>& dist_cost_i, int nb_threads, double upper_bound, double lower_bound, double& new_shift_x_min_i, double& new_shift_x_max_i, double& new_shift_y_min_i, double& new_shift_y_max_i, double& new_shift_z_min_i, double& new_shift_z_max_i, double& prep_time_i, int max_cells, const vector<vector<pair<double, double>>>& model_voronoi_min_max, const vector<vector<vector<int>>>& model_voronoi_vertex_edge, const vector<vector<vector<pair<int,int>>>>& model_voronoi_vertex_edge_planes, const param<double>& dii, const param<double>& djj, int iter, const vector<double>& max_vert_vert_dist_sq);
bool vert_enum(vector<vector<double>> vertex_set_a, vector<vector<double>> facets_a, vector<vector<int>> vertex_edge_a, vector<vector<pair<int, int>>> vertex_edge_plane_a, vector<vector<double>> vertex_set_b, vector<vector<double>> facets_b,  vector<vector<int>> vertex_edge_b, vector<vector<pair<int, int>>> vertex_edge_plane_b, std::vector<std::vector<double>>& new_vert, const std::vector<double>& model_point);
tuple<double,double,double> run_ARMO(string axis, const vector<vector<double>>& point_cloud1, const vector<vector<double>>& point_cloud2, const vector<vector<double>>& uav1, const vector<vector<double>>& uav2);

/* Plot two point clouds */
void plot(const vector<vector<double>>& point_cloud1, const vector<vector<double>>& point_cloud2, double point_thick = 0.1);
void plot(const vector<vector<double>>& point_cloud1, const vector<vector<double>>& point_cloud2,const vector<vector<double>>& point_cloud3, double point_thick = 0.1);
void plot(const vector<vector<double>>& ext_model, const vector<vector<double>>& ext_data, const vector<vector<double>>& ext_data1,const vector<vector<double>>& ext_data2, double point_thick=0.1);

/* Run Go-ICP on point clouds, return best solution */
tuple<double,double,double,double,double,double,double> run_GoICP(const vector<vector<double>>& point_cloud1, const vector<vector<double>>& point_cloud2);

/* Run the ARMO model for registration */
tuple<double,double,double,double,double,double> run_ARMO(bool bypass, string axis, const vector<vector<double>>& point_cloud1, const vector<vector<double>>& point_cloud2);

/* Run the Global ARMO model for registration */
vector<double> run_ARMO_Global(bool bypass, string axis, const vector<vector<double>>& point_cloud1, const vector<vector<double>>& point_cloud2, bool norm1 = false);



bool get_solution(const shared_ptr<Model<double>>& M, vector<double>& rot_trans, vector<int>& new_matching);
void get_rotation_transl_matrix(const shared_ptr<Model<double>>& M, vector<double>& rot_trans);
void get_angle_rotation_transl_matrix(const shared_ptr<Model<double>>& M, vector<double>& rot_trans);
void update_matching(shared_ptr<Model<double>>& M, vector<int>& new_matching);
void round_bin(shared_ptr<Model<double>>& M, int nd, int nm);



void initialize_model_from_parent(shared_ptr<Model<double>>parent, shared_ptr<Model<double>> child);

shared_ptr<Model<double>> build_norm2_SOC_MIQCP(vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, const indices& valid_cells, const indices& new_model_ids, const param<>& dist_cost, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double shift_min_x, double shift_max_x, double shift_min_y, double shift_max_y, double shift_min_z, double shift_max_z, vector<double>& rot_trans, bool convex, const vector<pair<pair<int,int>,pair<int,int>>>& incompatibles, param<>& norm_x, param<>& norm_y, param<>& norm_z,  param<>& intercept, const vector<int>& init_matching, const vector<double>& error_per_point, param<>& model_radius, bool relax_ints, bool relax_sdp = false, bool rigid_transf = true, double perc_outliers = 0);


shared_ptr<Model<double>> build_linobj_convex_clean(const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data, const indices& valid_cells, double new_roll_min, double new_roll_max, double new_pitch_min, double new_pitch_max, double new_yaw_min, double new_yaw_max, double new_shift_min_x, double new_shift_max_x, double new_shift_min_y, double new_shift_max_y, double new_shift_min_z, double new_shift_max_z, param<>& dist_cost, double ub, bool& status, double lb, param<double>& maxj);




vector<double> BranchBound11(GoICP& goicp, vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, vector<vector<vector<double>>> model_voronoi_normals, vector<vector<double>> model_face_intercept, vector<vector<vector<double>>> model_voronoi_vertices,  const vector<double>& model_inner_prod_min,const vector<double>& model_inner_prod_max, vector<pair<double, double>> min_max_model,  vector<double>& best_rot_trans, double best_ub, const vector<vector<pair<double, double>>>& model_voronoi_min_max, const vector<vector<vector<int>>>& model_voronoi_vertex_edge, const vector<vector<vector<pair<int,int>>>>& model_voronoi_vertex_edge_planes, const vector<double>& dm,  const param<double>& dii, const param<double>& djj);

//indices get_valid_pairs(vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double shift_min_x, double shift_max_x, double shift_min_y, double shift_max_y, double shift_min_z, double shift_max_z, vector<vector<vector<double>>> model_voronoi_normals, vector<vector<double>> model_face_intercept, const vector<double>& model_voronoi_out_radius, vector<vector<vector<double>>> model_voronoi_vertices, bool norm1);


pair<double,double> min_max_euclidean_distsq_box(vector<vector<double>> coords, vector<double> center);

pair<double,double> min_max_euclidean_distsq_box_plane(vector<vector<double>> coords, vector<vector<double>> new_coords, vector<double> center, vector<double> plane, double xl, double xu, double yl, double yu, double zl, double zu);



vector<pair<double,double>> center_point_cloud(vector<vector<double>>& point_cloud);



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

double build_upperbound(vector<double>& solutionlb,vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, const indices& valid_cells, double new_roll_min, double new_roll_max, double new_pitch_min, double new_pitch_max, double new_yaw_min, double new_yaw_max, double new_shift_min_x, double new_shift_max_x, double new_shift_min_y, double new_shift_max_y, double new_shift_min_z, double new_shift_max_z, vector<double>& rot_trans);
GoICP initialize_ICP_only(const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data);
double run_ICP_only(GoICP& goicp, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double shift_min_x, double shift_max_x, double shift_min_y, double shift_max_y, double shift_min_z, double shift_max_z, vector<double>& rot_trans_ub);

void compute_upper_boundICP(GoICP& goicp, double roll_mini, double roll_maxi, double pitch_mini, double pitch_maxi, double yaw_mini, double yaw_maxi, double shift_min_xi, double shift_max_xi, double shift_min_yi, double shift_max_yi, double shift_min_zi, double shift_max_zi, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double shift_min_x, double shift_max_x, double shift_min_y, double shift_max_y, double shift_min_z, double shift_max_z, vector<double>& best_rot_trans, double& best_ub, vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data);
GoICP Initialize_BB(vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, const vector<pair<double, double>>& min_max_model, vector<pair<double, double>>& min_max_data,double shift_min_x, double shift_max_x, double shift_min_y, double shift_max_y, double shift_min_z, double shift_max_z, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double& best_ub, vector<double>& best_rot_trans);

#ifdef USE_GJK
double distance_polytopes_gjk(vector<vector<double>>& vec1, vector<vector<double>>& vec2);
#endif
void add_vertex(std::vector<std::vector<double>>& new_vert, std::vector<double> solution, const std::vector<double>& model_point);

/* Return the min-max values for x, y and z  for all possible rotations of p with angle +- angle*/
vector<pair<double,double>> get_min_max(double angle, const vector<double>& p, const vector<double>& ref);

/* Run the ARMO model for boresight alignment */
tuple<double,double,double> run_ARMO(string axis, const vector<vector<double>>& point_cloud1, const vector<vector<double>>& point_cloud2, const vector<vector<double>>& uav1, const vector<vector<double>>& uav2);

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
    vector<vector<double>> full_point_cloud1, full_point_cloud2;
    vector<vector<double>> point_cloud1, point_cloud2;
    vector<vector<double>> full_uav1, full_uav2;
    vector<vector<double>> uav1, uav2;
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
        rapidcsv::Document  red_Model_doc(red_Model_file, rapidcsv::LabelParams(0, -1),rapidcsv::SeparatorParams(' '));
        read_data(red_Model_doc, point_cloud1, uav1);
    }
    if(argc>5){
        red_Data_file = argv[5];
        rapidcsv::Document  red_Data_doc(red_Data_file, rapidcsv::LabelParams(0, -1),rapidcsv::SeparatorParams(' '));
        read_data(red_Data_doc, point_cloud2, uav2);
        
    }
    rapidcsv::Document  Model_doc(Model_file, rapidcsv::LabelParams(0, -1),rapidcsv::SeparatorParams(' '));
    rapidcsv::Document  Data_doc(Data_file, rapidcsv::LabelParams(0, -1),rapidcsv::SeparatorParams(' '));
    
    
    read_data(Model_doc, full_point_cloud1, full_uav1);
    read_data(Data_doc, full_point_cloud2, full_uav2);
    
    
    
    bool run_goICP = false;
    if(run_goICP){/* Run GoICP inline */
        run_GoICP(point_cloud1, point_cloud2);
    }
    double final_roll = 0, final_pitch = 0, final_yaw = 0;
    double total_time =0, time_start = 0, time_end = 0;
    double L2error_init = 0, L1error_init = 0;
    vector<int> matching(point_cloud2.size());
    vector<double> err_per_point(point_cloud2.size());
    L2error_init = computeL2error(point_cloud1,point_cloud2,matching,err_per_point);
    L1error_init = computeL1error(point_cloud1,point_cloud2,matching,err_per_point);
    DebugOn("Initial L2 error = " << L2error_init << endl);
    DebugOn("Initial L1 error = " << L1error_init << endl);
    time_start = get_wall_time();
    auto res = run_IPH(point_cloud1, point_cloud2, uav1, uav2);
    final_roll = get<0>(res);final_pitch = get<1>(res); final_yaw = get<2>(res);
    auto L2error = computeL2error(point_cloud1,point_cloud2,matching,err_per_point);
    auto L1error = computeL1error(point_cloud1,point_cloud2,matching,err_per_point);
    DebugOn("n1 = " << point_cloud1.size() << endl);
    DebugOn("n2 = " << point_cloud2.size() << endl);
    DebugOn("Initial L2 error = " << L2error_init << endl);
    DebugOn("Final L2 error = " << L2error << endl);
    DebugOn("Relative L2 gap = " << 100*(L2error_init-L2error)/L2error_init << "%" << endl);
    DebugOn("Initial L1 error = " << L1error_init << endl);
    DebugOn("Final L1 error = " << L1error << endl);
    DebugOn("Relative L1 gap = " << 100*(L1error_init-L1error)/L1error_init << "%" << endl);
    time_end = get_wall_time();
    DebugOn("Total wall clock time = " << time_end - time_start << endl);
    double shifted_x, shifted_y, shifted_z;
    auto tot_pts = full_point_cloud1.size()+full_point_cloud2.size();
    apply_rotation(final_roll, final_pitch, final_yaw, full_point_cloud1, full_point_cloud2, full_uav1, full_uav2);
    
    bool save_file = true;
    if(save_file){
        auto name = Model_file.substr(0,Model_file.find('.'));
        auto fname = name+"_ARMO_RPY_"+to_string(final_roll)+"_"+to_string(final_pitch)+"_"+to_string(final_yaw)+".laz";
        save_laz(fname,full_point_cloud1, full_point_cloud2);
        fname = name+"_ARMO_RPY_"+to_string(final_roll)+"_"+to_string(final_pitch)+"_"+to_string(final_yaw)+"_sub.laz";
        save_laz(fname,point_cloud1, point_cloud2);
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
vector<pair<double,double>> get_min_max(double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, const vector<double>& p, const vector<double>& ref){
    double x1 = p[0], y1 = p[1], z1 = p[2], shifted_x, shifted_y, shifted_z, alpha, beta, gamma;
    double x_ref = ref[0], y_ref = ref[1], z_ref = ref[2];
    double x_rot1, y_rot1, z_rot1, x_min = numeric_limits<double>::max(), x_max = numeric_limits<double>::lowest(), y_min = numeric_limits<double>::max(), y_max = numeric_limits<double>::lowest(), z_min = numeric_limits<double>::max(), z_max = numeric_limits<double>::lowest();
    double angles_alpha[] = {0, yaw_min, yaw_max};
    double angles_beta[] = {0, roll_min, roll_max};
    double angles_gamma[] = {0, pitch_min, pitch_max};
    vector<pair<double,double>> min_max;
    shifted_x = x1 - x_ref;
    shifted_y = y1 - y_ref;
    shifted_z = z1 - z_ref;
    for (int a = 0; a < 3; a++){
        alpha = angles_alpha[a];
        for (int b = 0; b <3; b++){
            beta = angles_beta[b];
            for (int c = 0; c <3; c++){
                gamma = angles_gamma[c];
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

double get_GoICP_dist(double radius_r, double radius_t, const vector<double>& p, bool L1norm){
    //    radius_r *= sqrt(2);
    double t_radius = std::sqrt(3)*radius_t;
    DebugOff("GoICP t radius = " << to_string_with_precision(t_radius, 6) << endl);
    DebugOff("GoICP r radius = " << to_string_with_precision(2*std::sin(std::min(sqrt(3)*radius_r/2.,pi/2.))*std::sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]), 6) << endl);
    if(L1norm)
        return 2*std::sin(std::min(3*radius_r/2,pi/2))*(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]) + 3*radius_t;
    return 2*std::sin(std::min(sqrt(3)*radius_r/2.,pi/2.))*std::sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]) + std::sqrt(3)*radius_t;
}

/* Return the min-max values for x, y and z  for all possible rotations of p with angle +- angle assuming we're in [-pi,pi]*/
double get_max_dist(double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double tx_min, double tx_max, double ty_min, double ty_max, double tz_min, double tz_max, const vector<double>& p, const vector<double>& ref, bool L1norm){
    double x1 = p[0], y1 = p[1], z1 = p[2], shifted_x, shifted_y, shifted_z, yaw, roll, pitch;
    double x_ref = ref[0], y_ref = ref[1], z_ref = ref[2];
    double x_rot1, y_rot1, z_rot1, x_min = numeric_limits<double>::max(), x_max = numeric_limits<double>::lowest(), y_min = numeric_limits<double>::max(), y_max = numeric_limits<double>::lowest(), z_min = numeric_limits<double>::max(), z_max = numeric_limits<double>::lowest();
    vector<double> angles_yaw = {-pi, -pi/2, 0, pi/2, pi, yaw_min, yaw_max};
    vector<double> angles_roll = {-pi, -pi/2, 0, pi/2, pi, roll_min, roll_max};
    vector<double> angles_pitch = {-pi, -pi/2, 0, pi/2, pi, pitch_min, pitch_max};
    if((roll_min >=-pi && roll_max <= -pi/2) || (roll_min >=-pi/2 && roll_max <= 0) || (roll_min >=0 && roll_max <= pi/2)  || (roll_min >= pi/2 && roll_max <= pi)){
        angles_roll = {roll_min, roll_max};
    }
    else if(roll_min >=-pi && roll_max <= 0){
        angles_roll = {roll_min, roll_max, -pi/2};
    }
    else if(roll_min >=-pi/2 && roll_max <= pi/2){
        angles_roll = {roll_min, roll_max, 0};
    }
    else if(roll_min >=0 && roll_max <= pi){
        angles_roll = {roll_min, roll_max, pi/2};
    }
    else if(roll_min >=0){
        angles_roll = {roll_min, roll_max, pi/2, pi};
    }
    else if(roll_max <=0){
        angles_roll = {roll_min, roll_max, -pi/2, -pi};
    }
    if((pitch_min >=-pi && pitch_max <= -pi/2) || (pitch_min >=-pi/2 && pitch_max <= 0) || (pitch_min >=0 && pitch_max <= pi/2)  || (pitch_min >= pi/2 && pitch_max <= pi)){
        angles_pitch = {pitch_min, pitch_max};
    }
    else if(pitch_min >=-pi && pitch_max <= 0){
        angles_pitch = {pitch_min, pitch_max, -pi/2};
    }
    else if(pitch_min >=-pi/2 && pitch_max <= pi/2){
        angles_pitch = {pitch_min, pitch_max, 0};
    }
    else if(pitch_min >=0 && pitch_max <= pi){
        angles_pitch = {pitch_min, pitch_max, pi/2};
    }
    else if(pitch_min >=0){
        angles_pitch = {pitch_min, pitch_max, pi/2, pi};
    }
    else if(pitch_max <=0){
        angles_pitch = {pitch_min, pitch_max, -pi/2, -pi};
    }
    if((yaw_min >=-pi && yaw_max <= -pi/2) || (yaw_min >=-pi/2 && yaw_max <= 0) || (yaw_min >=0 && yaw_max <= pi/2)  || (yaw_min >= pi/2 && yaw_max <= pi)){
        angles_yaw = {yaw_min, yaw_max};
    }
    else if(yaw_min >=-pi && yaw_max <= 0){
        angles_yaw = {yaw_min, yaw_max, -pi/2};
    }
    else if(yaw_min >=-pi/2 && yaw_max <= pi/2){
        angles_yaw = {yaw_min, yaw_max, 0};
    }
    else if(yaw_min >=0 && yaw_max <= pi){
        angles_yaw = {yaw_min, yaw_max, pi/2};
    }
    else if(yaw_min >=0){
        angles_yaw = {yaw_min, yaw_max, pi/2, pi};
    }
    else if(yaw_max <=0){
        angles_yaw = {yaw_min, yaw_max, -pi/2, -pi};
    }
    double tx[] = {tx_min, tx_max};
    double ty[] = {ty_min, ty_max};
    double tz[] = {tz_min, tz_max};
    double max_dist = 0, dist = 0;
    shifted_x = x1 - x_ref;
    shifted_y = y1 - y_ref;
    shifted_z = z1 - z_ref;
    for (int a = 0; a < angles_yaw.size(); a++){
        yaw = angles_yaw[a];
        for (int b = 0; b <angles_roll.size(); b++){
            roll = angles_roll[b];
            for (int c = 0; c <angles_pitch.size(); c++){
                pitch = angles_pitch[c];
                x_rot1 = shifted_x*cos(yaw)*cos(roll) + shifted_y*(cos(yaw)*sin(roll)*sin(pitch) - sin(yaw)*cos(pitch)) + shifted_z*(cos(yaw)*sin(roll)*cos(pitch) + sin(yaw)*sin(pitch));
                y_rot1 = shifted_x*sin(yaw)*cos(roll) + shifted_y*(sin(yaw)*sin(roll)*sin(pitch) + cos(yaw)*cos(pitch)) + shifted_z*(sin(yaw)*sin(roll)*cos(pitch) - cos(yaw)*sin(pitch));
                z_rot1 = shifted_x*(-sin(roll)) + shifted_y*(cos(roll)*sin(pitch)) + shifted_z*(cos(roll)*cos(pitch));
                x_rot1 += x_ref;
                y_rot1 += y_ref;
                z_rot1 += z_ref;
                if(L1norm){
                    dist = std::abs(x1 - x_rot1) + std::abs(y1 - y_rot1) + std::abs(z1 - z_rot1);
                }
                else{
                    dist = std::pow(x1 - x_rot1,2) + std::pow(y1 - y_rot1,2) + std::pow(z1 - z_rot1,2);
                }
                if(dist>max_dist){
                    max_dist = dist;
                }
            }
        }
    }
    double t_radius = std::sqrt(std::pow(std::max(std::abs(tx_min),std::abs(tx_max)),2) + std::pow(std::max(std::abs(ty_min),std::abs(ty_max)),2) + std::pow(std::max(std::abs(tz_min),std::abs(tz_max)),2));
    DebugOff("Our t radius = " << to_string_with_precision(t_radius, 6) << endl);
    DebugOff("Our r radius = " << to_string_with_precision(std::sqrt(max_dist), 6) << endl);
    if(!L1norm)
        return std::sqrt(max_dist) + t_radius;
    return max_dist + std::max(std::abs(tx_min),std::abs(tx_max)) + std::max(std::abs(ty_min),std::abs(ty_max)) + std::max(std::abs(tz_min),std::abs(tz_max));
}






void min_max_dist_box(double x0, double y0, double z0,  double xl, double xu, double yl, double yu, double zl, double zu, double& dbox_min, double& dbox_max){
    
    vector<double> d;
    
    d.push_back(pow(xu-x0,2)+pow(yu-y0,2)+pow(zu-z0,2));
    d.push_back(pow(xl-x0,2)+pow(yu-y0,2)+pow(zu-z0,2));
    d.push_back(pow(xl-x0,2)+pow(yl-y0,2)+pow(zu-z0,2));
    d.push_back(pow(xu-x0,2)+pow(yl-y0,2)+pow(zu-z0,2));
    d.push_back(pow(xu-x0,2)+pow(yu-y0,2)+pow(zl-z0,2));
    d.push_back(pow(xl-x0,2)+pow(yu-y0,2)+pow(zl-z0,2));
    d.push_back(pow(xl-x0,2)+pow(yl-y0,2)+pow(zl-z0,2));
    d.push_back(pow(xu-x0,2)+pow(yl-y0,2)+pow(zl-z0,2));
    
    auto a=std::minmax_element(d.begin(), d.end());
    dbox_min=*a.first;
    dbox_max=*a.second;
    
    if((xl<=x0) && (xu>=x0) && (yl<=y0) && (yu>=y0) && (zl<=z0) && (zu>=z0)){
        dbox_min=0;
    }
    
}




void update_theta_bounds(shared_ptr<Model<double>>& m, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max){
    var<> yaw("yaw", yaw_min, yaw_max), pitch("pitch", pitch_min, pitch_max), roll("roll", roll_min, roll_max);
    yaw.in(R(1)); pitch.in(R(1));roll.in(R(1));
    func<> r11 = cos(yaw)*cos(roll);
    func<> r12 = cos(yaw)*sin(roll)*sin(pitch) - sin(yaw)*cos(pitch);
    func<> r13 = cos(yaw)*sin(roll)*cos(pitch) + sin(yaw)*sin(pitch);
    func<> r21 = sin(yaw)*cos(roll);
    func<> r22 = sin(yaw)*sin(roll)*sin(pitch) + cos(yaw)*cos(pitch);
    func<> r23 = sin(yaw)*sin(roll)*cos(pitch) - cos(yaw)*sin(pitch);
    func<> r31 = sin(-1*roll);
    func<> r32 = cos(roll)*sin(pitch);
    func<> r33 = cos(roll)*cos(pitch);
    auto theta11 = m->get_ptr_var<double>("theta11"); auto theta12 = m->get_ptr_var<double>("theta12"); auto theta13 = m->get_ptr_var<double>("theta13");
    auto theta21 = m->get_ptr_var<double>("theta21"); auto theta22 = m->get_ptr_var<double>("theta22"); auto theta23 = m->get_ptr_var<double>("theta23");
    auto theta31 = m->get_ptr_var<double>("theta31"); auto theta32 = m->get_ptr_var<double>("theta32"); auto theta33 = m->get_ptr_var<double>("theta33");
    func<> cosr_f = cos(roll);
    func<> sinr_f = sin(roll);
    func<> cosp_f = cos(pitch);
    func<> sinp_f = sin(pitch);
    func<> cosy_f = cos(yaw);
    func<> siny_f = sin(yaw);
    auto cosr = m->get_ptr_var<double>("cosr");
    cosr->set_bounds(cosr_f._range->first, cosr_f._range->second);
    auto sinr = m->get_ptr_var<double>("sinr");
    sinr->set_bounds(sinr_f._range->first, sinr_f._range->second);
    auto cosp = m->get_ptr_var<double>("cosp");
    cosp->set_bounds(cosp_f._range->first, cosp_f._range->second);
    auto sinp = m->get_ptr_var<double>("sinp");
    sinp->set_bounds(sinp_f._range->first, sinp_f._range->second);
    auto cosy = m->get_ptr_var<double>("cosy");
    cosy->set_bounds(cosy_f._range->first, cosy_f._range->second);
    auto siny = m->get_ptr_var<double>("siny");
    siny->set_bounds(siny_f._range->first, siny_f._range->second);
    auto cosy_sinr_range = get_product_range(cosy->_range, sinr->_range);
    auto siny_sinr_range = get_product_range(siny->_range, sinr->_range);
    auto cosy_sinr = m->get_ptr_var<double>("cosy_sinr");
    cosy_sinr->set_bounds(cosy_sinr_range->first, cosy_sinr_range->second);
    auto siny_sinr = m->get_ptr_var<double>("siny_sinr");
    siny_sinr->set_bounds(siny_sinr_range->first, siny_sinr_range->second);
    
    theta11->set_bounds(std::max(-1.,r11._range->first), std::min(1.,r11._range->second));
    theta12->set_bounds(std::max(-1.,r12._range->first), std::min(1.,r12._range->second));
    theta13->set_bounds(std::max(-1.,r13._range->first), std::min(1.,r13._range->second));
    theta21->set_bounds(std::max(-1.,r21._range->first), std::min(1.,r21._range->second));
    theta21->set_bounds(std::max(-1.,r22._range->first), std::min(1.,r22._range->second));
    theta23->set_bounds(std::max(-1.,r23._range->first), std::min(1.,r23._range->second));
    theta31->set_bounds(std::max(-1.,r31._range->first), std::min(1.,r31._range->second));
    theta32->set_bounds(std::max(-1.,r32._range->first), std::min(1.,r32._range->second));
    theta33->set_bounds(std::max(-1.,r33._range->first), std::min(1.,r33._range->second));
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

void round_bin(shared_ptr<Model<double>>& M, int nd, int nm){
    auto bin = M->get_ptr_var<double>("bin");
    shared_ptr<vector<double>> cont_vals = bin->_val;
    int idx = 0;
    for (int i = 0; i<nd; i++) {
        double max = 0;
        int max_i = 0;
        for (int j = 0; j<nm; j++) {
            if(max<cont_vals->at(idx)){
                max = cont_vals->at(idx);
                max_i = idx;
            }
            idx++;
        }
        idx -= nm;
        for (int j = 0; j<nm; j++) {
            cont_vals->at(idx++) = 0;
        }
        cont_vals->at(max_i) = 1;
    }
    bin->fix();
}

void update_matching(shared_ptr<Model<double>>& M, vector<int>& new_matching){
    auto bin = M->get_var_ptr("bin");
    shared_ptr<vector<int>> bin_vals = nullptr;
    shared_ptr<vector<double>> cont_vals = nullptr;
    if(bin->is_integer()){
        bin_vals = static_pointer_cast<var<int>>(bin)->_val;
    }
    else{/* integers were replaced by continuous vars*/
        bin_vals = static_pointer_cast<var<int>>(M->get_int_var(bin->get_id()))->_val;
        cont_vals = static_pointer_cast<var<double>>(bin)->_val;
    }
    int nd = new_matching.size();
    int nm = bin_vals->size()/nd;
    int idx = 0;
    for (int i = 0; i<nd; i++) {
        for (int j = 0; j<nm; j++) {
            if(new_matching[i]==j){
                bin_vals->at(idx)=1;
                if(cont_vals)
                    cont_vals->at(idx)=1;
            }
            else{
                bin_vals->at(idx)=0;
                if(cont_vals)
                    cont_vals->at(idx)=0;
            }
            idx++;
        }
    }
}

void get_rotation_transl_matrix(const shared_ptr<Model<double>>& M, vector<double>& rot_trans){
    auto theta11 = M->get_var<double>("theta11");auto theta12 = M->get_var<double>("theta12");auto theta13 = M->get_var<double>("theta13");
    auto theta21 = M->get_var<double>("theta21");auto theta22 = M->get_var<double>("theta22");auto theta23 = M->get_var<double>("theta23");
    auto theta31 = M->get_var<double>("theta31");auto theta32 = M->get_var<double>("theta32");auto theta33 = M->get_var<double>("theta33");
    auto x_shift = M->get_var<double>("x_shift");auto y_shift = M->get_var<double>("y_shift");auto z_shift = M->get_var<double>("z_shift");
    
    Debug("Theta matrix = " << endl);
    Debug("|" << theta11.eval() << " " << theta12.eval() << " " << theta13.eval() << "|" << endl);
    Debug("|" << theta21.eval() << " " << theta22.eval() << " " << theta23.eval() << "|" << endl);
    Debug("|" << theta31.eval() << " " << theta32.eval() << " " << theta33.eval() << "|" << endl);
    
    rot_trans[0]=theta11.eval();
    rot_trans[1]=theta12.eval();
    rot_trans[2]=theta13.eval();;
    rot_trans[3]=theta21.eval();
    rot_trans[4]=theta22.eval();
    rot_trans[5]=theta23.eval();
    rot_trans[6]=theta31.eval();
    rot_trans[7]=theta32.eval();
    rot_trans[8]=theta33.eval();
    rot_trans[9]=x_shift.eval();
    rot_trans[10]=y_shift.eval();
    rot_trans[11]=z_shift.eval();
}

void get_angle_rotation_transl_matrix(const shared_ptr<Model<double>>& M, vector<double>& rot_trans){
    auto theta11 = M->get_var<double>("theta11");auto theta12 = M->get_var<double>("theta12");auto theta13 = M->get_var<double>("theta13");
    auto theta21 = M->get_var<double>("theta21");auto theta22 = M->get_var<double>("theta22");auto theta23 = M->get_var<double>("theta23");
    auto theta31 = M->get_var<double>("theta31");auto theta32 = M->get_var<double>("theta32");auto theta33 = M->get_var<double>("theta33");
    auto x_shift = M->get_var<double>("x_shift");auto y_shift = M->get_var<double>("y_shift");auto z_shift = M->get_var<double>("z_shift");
    
    Debug("Theta matrix = " << endl);
    Debug("|" << theta11.eval() << " " << theta12.eval() << " " << theta13.eval() << "|" << endl);
    Debug("|" << theta21.eval() << " " << theta22.eval() << " " << theta23.eval() << "|" << endl);
    Debug("|" << theta31.eval() << " " << theta32.eval() << " " << theta33.eval() << "|" << endl);
    
    auto pitch_val = std::atan2(theta32.eval(), theta33.eval())*180/pi;
    auto roll_val = std::atan2(-1*theta31.eval(), std::sqrt(theta32.eval()*theta32.eval()+theta33.eval()*theta33.eval()))*180/pi;
    auto yaw_val = std::atan2(theta21.eval(),theta11.eval())*180/pi;
    
    rot_trans[0]=roll_val;
    rot_trans[1]=pitch_val;
    rot_trans[2]=yaw_val;
    rot_trans[3]=x_shift.eval();
    rot_trans[4]=y_shift.eval();
    rot_trans[5]=z_shift.eval();
}

bool get_solution(const shared_ptr<Model<double>>& M, vector<double>& rot_trans, vector<int>& new_matching){
    auto theta11 = M->get_var<double>("theta11");auto theta12 = M->get_var<double>("theta12");auto theta13 = M->get_var<double>("theta13");
    auto theta21 = M->get_var<double>("theta21");auto theta22 = M->get_var<double>("theta22");auto theta23 = M->get_var<double>("theta23");
    auto theta31 = M->get_var<double>("theta31");auto theta32 = M->get_var<double>("theta32");auto theta33 = M->get_var<double>("theta33");
    auto x_shift = M->get_var<double>("x_shift");auto y_shift = M->get_var<double>("y_shift");auto z_shift = M->get_var<double>("z_shift");
    auto bin = M->get_var_ptr("bin");
    shared_ptr<vector<int>> bin_vals;
    if(bin->is_integer()){
        bin_vals = static_pointer_cast<var<int>>(bin)->_val;
    }
    else{/* integers were replaced by continuous vars*/
        bin_vals = static_pointer_cast<var<int>>(M->get_int_var(bin->get_id()))->_val;
    }
    int nd = new_matching.size();
    int nm = bin_vals->size()/nd;
    int idx = 0;
    for (int i = 0; i<nd; i++) {
        for (int j = 0; j<nm; j++) {
            if(bin_vals->at(idx++)==1){
                new_matching[i]=j;
            }
        }
    }
    DebugOff("Theta matrix = " << endl);
    DebugOff("|" << theta11.eval() << " " << theta12.eval() << " " << theta13.eval() << "|" << endl);
    DebugOff("|" << theta21.eval() << " " << theta22.eval() << " " << theta23.eval() << "|" << endl);
    DebugOff("|" << theta31.eval() << " " << theta32.eval() << " " << theta33.eval() << "|" << endl);
    constant<> row1 = pow(theta11.eval(),2)+pow(theta12.eval(),2)+pow(theta13.eval(),2);
    constant<> row2 = pow(theta21.eval(),2)+pow(theta22.eval(),2)+pow(theta23.eval(),2);
    constant<> row3 = pow(theta31.eval(),2)+pow(theta32.eval(),2)+pow(theta33.eval(),2);
    constant<> col1 = pow(theta11.eval(),2)+pow(theta21.eval(),2)+pow(theta31.eval(),2);
    constant<> col2 = pow(theta12.eval(),2)+pow(theta22.eval(),2)+pow(theta32.eval(),2);
    constant<> col3 = pow(theta13.eval(),2)+pow(theta23.eval(),2)+pow(theta33.eval(),2);
    DebugOff("row 1 " << row1.eval() << endl);
    DebugOff("row 2 " << row2.eval() << endl);
    DebugOff("row 3 " << row3.eval() << endl);
    DebugOff("col 1 " << col1.eval() << endl);
    DebugOff("col 2 " << col2.eval() << endl);
    DebugOff("col 3 " << col3.eval() << endl);
    constant<> det=theta11.eval()*(theta22.eval()*theta33.eval()-theta32.eval()*theta23.eval())
    -theta12.eval()*(theta21.eval()*theta33.eval()-theta31.eval()*theta23.eval())+theta13.eval()*(theta21.eval()*theta32.eval()-theta31.eval()*theta22.eval());
    constant<> row12 = (theta11.eval()*theta21.eval())+(theta12.eval()*theta22.eval())+(theta13.eval()*theta23.eval());
    constant<> row13 = (theta11.eval()*theta31.eval())+(theta12.eval()*theta32.eval())+(theta13.eval()*theta33.eval());
    constant<> row23 = (theta21.eval()*theta31.eval())+(theta22.eval()*theta32.eval())+(theta23.eval()*theta33.eval());
    DebugOff("row 12 " << row12.eval() << endl);
    DebugOff("row 13 " << row13.eval() << endl);
    DebugOff("row 23 " << row23.eval() << endl);
    
    DebugOff("Determinant "<<det.eval()<<endl);
    
    bool is_rotation = row1.is_approx(1) && row2.is_approx(1) && row3.is_approx(1) && col1.is_approx(1) && col2.is_approx(1) && col3.is_approx(1) && det.is_approx(1) && row12.is_approx(0) && row13.is_approx(0) && row23.is_approx(0);
    
    rot_trans[0]=theta11.eval();
    rot_trans[1]=theta12.eval();
    rot_trans[2]=theta13.eval();;
    rot_trans[3]=theta21.eval();
    rot_trans[4]=theta22.eval();
    rot_trans[5]=theta23.eval();
    rot_trans[6]=theta31.eval();
    rot_trans[7]=theta32.eval();
    rot_trans[8]=theta33.eval();
    rot_trans[9]=x_shift.eval();
    rot_trans[10]=y_shift.eval();
    rot_trans[11]=z_shift.eval();
    if(rot_trans.size()>12){
        auto scale_x = M->get_var<double>("scale_x");auto scale_y = M->get_var<double>("scale_y");auto scale_z = M->get_var<double>("scale_z");
        rot_trans[12]=scale_x.eval();
        rot_trans[13]=scale_y.eval();
        rot_trans[14]=scale_z.eval();
        DebugOn("scale_x = " << scale_x.eval() << endl);
        DebugOn("scale_y = " << scale_y.eval() << endl);
        DebugOn("scale_z = " << scale_z.eval() << endl);
    }
    
    auto pitch_val = std::atan2(theta32.eval(), theta33.eval())*180/pi;
    auto roll_val = std::atan2(-1*theta31.eval(), std::sqrt(theta32.eval()*theta32.eval()+theta33.eval()*theta33.eval()))*180/pi;
    auto yaw_val = std::atan2(theta21.eval(),theta11.eval())*180/pi;
    DebugOff("Roll (degrees) = " << to_string_with_precision(roll_val,12) << endl);
    DebugOff("Pitch (degrees) = " << to_string_with_precision(pitch_val,12) << endl);
    DebugOff("Yaw (degrees) = " << to_string_with_precision(yaw_val,12) << endl);
    DebugOff("x shift = " << x_shift.eval() << endl);
    DebugOff("y shift = " << y_shift.eval() << endl);
    DebugOff("z shift = " << z_shift.eval() << endl);
    if(!is_rotation){
        DebugOn("WARNING, returned matrix is not a Rotation!\n");
    }
    return is_rotation;
}

tuple<double,double,double> run_ARMO(string axis, const vector<vector<double>>& point_cloud1, const vector<vector<double>>& point_cloud2, const vector<vector<double>>& uav1, const vector<vector<double>>& uav2){
    double angle_max = 0.05;
    int nb_pairs = 0, min_nb_pairs = numeric_limits<int>::max(), max_nb_pairs = 0, av_nb_pairs = 0;
    vector<pair<double,double>> min_max1;
    vector<vector<pair<double,double>>> min_max2(point_cloud2.size());
    vector<int> nb_neighbors(point_cloud1.size());
    vector<vector<int>> neighbors;
    /* Compute cube for all points in point cloud 2 */
    for (auto i = 0; i<point_cloud2.size(); i++) {
        min_max2[i] = get_min_max(angle_max, point_cloud2.at(i), uav2.at(i));
    }
    double roll_1 = 0, yaw_1 = 0, pitch_1 = 0;
    bool bypass = false;
    if(!bypass){
        /* Check if cubes intersect */
        neighbors.resize(point_cloud1.size());
        for (auto i = 0; i<point_cloud1.size(); i++) {
            nb_pairs = 0;
            min_max1 = get_min_max(angle_max, point_cloud1.at(i), uav1.at(i));
            DebugOff("For point (" << point_cloud1.at(i).at(0) << "," <<  point_cloud1.at(i).at(1) << "," << point_cloud1.at(i).at(2) << "): ");
            DebugOff("\n neighbors in umbrella : \n");
            for (size_t j = 0; j < point_cloud2.size(); j++){
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
        av_nb_pairs /= point_cloud1.size();
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
        vector<double> min_dist(point_cloud1.size(),numeric_limits<double>::max());
        vector<int> nearest(point_cloud1.size());
        vector<string> nearest_id(point_cloud1.size());
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
        for (auto i = 0; i<point_cloud1.size(); i++) {
            if(solve_lidar_iter)
                nb_max_neigh = 1;
            else
                nb_max_neigh = m;
            if(nb_neighbors[i]>=nb_max_neigh){
                i_str = to_string(idx1+1);
                x_uav1.add_val(i_str,uav1.at(i)[0]);
                x1.add_val(i_str,point_cloud1.at(i)[0]);
                y_uav1.add_val(i_str,uav1.at(i)[1]);
                y1.add_val(i_str,point_cloud1.at(i)[1]);
                z_uav1.add_val(i_str,uav1.at(i)[2]);
                z1.add_val(i_str,point_cloud1.at(i)[2]);
                if(solve_lidar_iter){
                    nb_max_neigh = nb_neighbors[i];
                }
                for (auto j = 0; j<nb_max_neigh; j++) {
                    auto k = neighbors[i].at(j);
                    auto res = n2_map.find(k);
                    if(res==n2_map.end()){
                        n2_map[k] = idx2;
                        j_str = to_string(idx2+1);
                        x_uav2.add_val(j_str,uav2.at(k)[0]);
                        x2.add_val(j_str,point_cloud2.at(k)[0]);
                        y_uav2.add_val(j_str,uav2.at(k)[1]);
                        y2.add_val(j_str,point_cloud2.at(k)[1]);
                        z_uav2.add_val(j_str,uav2.at(k)[2]);
                        z2.add_val(j_str,point_cloud2.at(k)[2]);
                        idx2++;
                    }
                    else {
                        j_str = to_string(res->second+1);
                    }
                    if(axis=="x")
                        dist_sq = std::pow(point_cloud1.at(i)[1] - point_cloud2.at(k)[1],2) + std::pow(point_cloud1.at(i)[2] - point_cloud2.at(k)[2],2);
                    else if(axis=="y")
                        dist_sq = std::pow(point_cloud1.at(i)[0] - point_cloud2.at(k)[0],2) + std::pow(point_cloud1.at(i)[2] - point_cloud2.at(k)[2],2);
                    else if(axis=="z")
                        dist_sq = std::pow(point_cloud1.at(i)[0] - point_cloud2.at(k)[0],2) + std::pow(point_cloud1.at(i)[1] - point_cloud2.at(k)[1],2);
                    else
                        dist_sq = std::pow(point_cloud1.at(i)[0] - point_cloud2.at(k)[0],2) + std::pow(point_cloud1.at(i)[1] - point_cloud2.at(k)[1],2) + std::pow(point_cloud1.at(i)[2] - point_cloud2.at(k)[2],2);
                    
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
            for (auto i = 0; i<point_cloud1.size(); i++) {
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
//            Sm1 = indices(N1, range(m-1,m-1));
//            S2m2 = indices(N1, range(2,m-2));
//            S3m1 = indices(N1, range(3,m-1));
//            K = indices(N1,range(2,m-1));
        
            
        
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
            bool L2Norm = false;
            bool L1Norm = !L2Norm;
            
            if(L2Norm){
                Constraint<> XNorm2("XNorm2");
                XNorm2 += x_diff - pow(new_x1.from(cells) - new_x2.to(cells),2);
                Lidar.add(XNorm2.in(cells)>=0);
                
                Constraint<> YNorm2("YNorm2");
                YNorm2 += y_diff - pow(new_y1.from(cells) - new_y2.to(cells),2);
                Lidar.add(YNorm2.in(cells)>=0);
                
                Constraint<> ZNorm2("ZNorm2");
                ZNorm2 += z_diff - pow(new_z1.from(cells) - new_z2.to(cells),2);
                Lidar.add(ZNorm2.in(cells)>=0);
            }
            if(L1Norm){
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
            }
            
            auto ids1 = yaw1.repeat_id(cells.size());
            auto ids2 = yaw2.repeat_id(cells.size());
            
            /* alpha = yaw_, beta = pitch_ and gamma = roll_ */
            if(axis!="x"){
                Constraint<> x_rot1("x_rot1");
                x_rot1 += new_x1 - x_uav1.in(N1);
                x_rot1 -= (x1.in(N1)-x_uav1.in(N1))*cos(yaw1.in(ids1))*cos(roll1.in(ids1)) + (y1.in(N1)-y_uav1.in(N1))*(cos(yaw1.in(ids1))*sin(roll1.in(ids1))*sin(pitch1.in(ids1)) - sin(yaw1.in(ids1))*cos(pitch1.in(ids1))) + (z1.in(N1)-z_uav1.in(N1))*(cos(yaw1.in(ids1))*sin(roll1.in(ids1))*cos(pitch1.in(ids1)) + sin(yaw1.in(ids1))*sin(pitch1.in(ids1)));
                Lidar.add(x_rot1.in(N1)==0);
                
            
                Constraint<> x_rot2("x_rot2");
                x_rot2 += new_x2 - x_uav2.in(N2);
                x_rot2 -= (x2.in(N2)-x_uav2.in(N2))*cos(yaw2.in(ids2))*cos(roll2.in(ids2)) + (y2.in(N2)-y_uav2.in(N2))*(cos(yaw2.in(ids2))*sin(roll2.in(ids2))*sin(pitch2.in(ids2)) - sin(yaw2.in(ids2))*cos(pitch2.in(ids2))) + (z2.in(N2)-z_uav2.in(N2))*(cos(yaw2.in(ids2))*sin(roll2.in(ids2))*cos(pitch2.in(ids2)) + sin(yaw2.in(ids2))*sin(pitch2.in(ids2)));
                Lidar.add(x_rot2.in(N2)==0);
            }
            
            if(axis!="y"){
                Constraint<> y_rot1("y_rot1");
                y_rot1 += new_y1 - y_uav1.in(N1);
                y_rot1 -= (x1.in(N1)-x_uav1.in(N1))*sin(yaw1.in(ids1))*cos(roll1.in(ids1)) + (y1.in(N1)-y_uav1.in(N1))*(sin(yaw1.in(ids1))*sin(roll1.in(ids1))*sin(pitch1.in(ids1)) + cos(yaw1.in(ids1))*cos(pitch1.in(ids1))) + (z1.in(N1)-z_uav1.in(N1))*(sin(yaw1.in(ids1))*sin(roll1.in(ids1))*cos(pitch1.in(ids1)) - cos(yaw1.in(ids1))*sin(pitch1.in(ids1)));
                Lidar.add(y_rot1.in(N1)==0);
                
                Constraint<> y_rot2("y_rot2");
                y_rot2 += new_y2 - y_uav2.in(N2);
                y_rot2 -= (x2.in(N2)-x_uav2.in(N2))*sin(yaw2.in(ids2))*cos(roll2.in(ids2)) + (y2.in(N2)-y_uav2.in(N2))*(sin(yaw2.in(ids2))*sin(roll2.in(ids2))*sin(pitch2.in(ids2)) + cos(yaw2.in(ids2))*cos(pitch2.in(ids2))) + (z2.in(N2)-z_uav2.in(N2))*(sin(yaw2.in(ids2))*sin(roll2.in(ids2))*cos(pitch2.in(ids2)) - cos(yaw2.in(ids2))*sin(pitch2.in(ids2)));
                Lidar.add(y_rot2.in(N2)==0);
            }
            
            if(axis!="z"){
                Constraint<> z_rot1("z_rot1");
                z_rot1 += new_z1 - z_uav1.in(N1);
                z_rot1 -= (x1.in(N1)-x_uav1.in(N1))*sin(-1*roll1.in(ids1)) + (y1.in(N1)-y_uav1.in(N1))*(cos(roll1.in(ids1))*sin(pitch1.in(ids1))) + (z1.in(N1)-z_uav1.in(N1))*(cos(roll1.in(ids1))*cos(pitch1.in(ids1)));
                Lidar.add(z_rot1.in(N1)==0);
                
                
                Constraint<> z_rot2("z_rot2");
                z_rot2 += new_z2 - z_uav2.in(N2);
                z_rot2 -= (x2.in(N2)-x_uav2.in(N2))*sin(-1*roll2.in(ids2)) + (y2.in(N2)-y_uav2.in(N2))*(cos(roll2.in(ids2))*sin(pitch2.in(ids2))) + (z2.in(N2)-z_uav2.in(N2))*(cos(roll2.in(ids2))*cos(pitch2.in(ids2)));
                Lidar.add(z_rot2.in(N2)==0);
            }
            
            //    M.min(sum(z_diff)/nb_overlap);
            
            //        M.min(sum(z_diff));
            if(axis == "full")
                Lidar.min(sum(x_diff)/cells.size() + sum(y_diff)/cells.size() + sum(z_diff)/cells.size());
            else if(axis == "x"){
                Lidar.min(sum(y_diff)/cells.size() + sum(z_diff)/cells.size());
                roll1.set_lb(0);
                roll1.set_ub(0);
                yaw1.set_lb(0);
                yaw1.set_ub(0);
//                x1.set_val(0);
//                x2.set_val(0);
            }
            else if (axis == "y") {
                Lidar.min(sum(x_diff)/cells.size() + sum(z_diff)/cells.size());
                yaw1.set_lb(0);
                yaw1.set_ub(0);
                pitch1.set_lb(0);
                pitch1.set_ub(0);
//                y1.set_val(0);
//                y2.set_val(0);
            }
            else if (axis == "z") {
                Lidar.min(sum(x_diff)/cells.size() + sum(y_diff)/cells.size());
                pitch1.set_lb(0);
                pitch1.set_ub(0);
                roll1.set_lb(0);
                roll1.set_ub(0);
//                z1.set_val(0);
//                z2.set_val(0);
            }
            else if (axis == "only_x")
                Lidar.min(sum(x_diff)/cells.size());
            else if (axis == "only_y")
                Lidar.min(sum(y_diff)/cells.size());
            else if (axis == "only_z")
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
        else if (solve_lidar_cube) {
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
    }
    return {roll_1, pitch_1, yaw_1};
}

tuple<double,double,double,double,double,double> run_ARMO(bool bypass, string axis, const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data){
    auto thetax = atan2(-0.0081669, -0.0084357)/2;
    auto thetay = atan2(0.9999311, std::sqrt(0.0081669*0.0081669+0.0084357*0.0084357))/2;
    auto thetaz = atan2(-0.0081462,-0.0084556)/2;
    DebugOff("thetax = " << thetax << endl);
    DebugOff("thetay = " << thetay << endl);
    DebugOff("thetaz = " << thetaz << endl);
    
    if(!bypass){
        double angle_max = 1, shift_max = 0.25;
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
//            S.run();
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
/* Update point cloud coordinates */
void update_xyz(vector<vector<double>>& point_cloud1, vector<double>& x_vec1, vector<double>& y_vec1, vector<double>& z_vec1){
    for (auto i = 0; i< point_cloud1.size(); i++) {
        point_cloud1[i][0] = x_vec1[i];
        point_cloud1[i][1] = y_vec1[i];
        point_cloud1[i][2] = z_vec1[i];
    }
}

void apply_rotation(double roll, double pitch, double yaw, vector<vector<double>>& point_cloud1, vector<vector<double>>& point_cloud2, const vector<vector<double>>& uav1, const vector<vector<double>>& uav2){
    double beta = roll*pi/180;// roll in radians
    double gamma = pitch*pi/180; // pitch in radians
    double alpha = yaw*pi/180; // yaw in radians
    double shifted_x, shifted_y, shifted_z;
    /* Apply rotation */
    for (auto i = 0; i< point_cloud1.size(); i++) {
        shifted_x = point_cloud1[i][0] - uav1[i][0];
        shifted_y = point_cloud1[i][1] - uav1[i][1];
        shifted_z = point_cloud1[i][2] - uav1[i][2];
        point_cloud1[i][0] = shifted_x*cos(alpha)*cos(beta) + shifted_y*(cos(alpha)*sin(beta)*sin(gamma) - sin(alpha)*cos(gamma)) + shifted_z*(cos(alpha)*sin(beta)*cos(gamma) + sin(alpha)*sin(gamma));
        point_cloud1[i][1] = shifted_x*sin(alpha)*cos(beta) + shifted_y*(sin(alpha)*sin(beta)*sin(gamma) + cos(alpha)*cos(gamma)) + shifted_z*(sin(alpha)*sin(beta)*cos(gamma) - cos(alpha)*sin(gamma));
        point_cloud1[i][2] = shifted_x*(-sin(beta)) + shifted_y*(cos(beta)*sin(gamma)) + shifted_z*(cos(beta)*cos(gamma));
        point_cloud1[i][0] += uav1[i][0];
        point_cloud1[i][1] += uav1[i][1];
        point_cloud1[i][2] += uav1[i][2];
    }
    beta *= -1;
    alpha *= -1;
    for (auto i = 0; i< point_cloud2.size(); i++) {
        shifted_x = point_cloud2[i][0] - uav2[i][0];
        shifted_y = point_cloud2[i][1] - uav2[i][1];
        shifted_z = point_cloud2[i][2] - uav2[i][2];
        point_cloud2[i][0] = shifted_x*cos(alpha)*cos(beta) + shifted_y*(cos(alpha)*sin(beta)*sin(gamma) - sin(alpha)*cos(gamma)) + shifted_z*(cos(alpha)*sin(beta)*cos(gamma) + sin(alpha)*sin(gamma));
        point_cloud2[i][1] = shifted_x*sin(alpha)*cos(beta) + shifted_y*(sin(alpha)*sin(beta)*sin(gamma) + cos(alpha)*cos(gamma)) + shifted_z*(sin(alpha)*sin(beta)*cos(gamma) - cos(alpha)*sin(gamma));
        point_cloud2[i][2] = shifted_x*(-sin(beta)) + shifted_y*(cos(beta)*sin(gamma)) + shifted_z*(cos(beta)*cos(gamma));
        point_cloud2[i][0] += uav2[i][0];
        point_cloud2[i][1] += uav2[i][1];
        point_cloud2[i][2] += uav2[i][2];
    }
}

void apply_rot_trans(const vector<double>& theta_matrix, vector<vector<double>>& point_cloud){
    double shifted_x, shifted_y, shifted_z;
    size_t n = point_cloud.size();
    /* Apply rotation */
    for (auto i = 0; i< n; i++) {
        shifted_x = point_cloud[i][0];
        shifted_y = point_cloud[i][1];
        shifted_z = point_cloud[i][2];
        point_cloud[i][0] = shifted_x*theta_matrix[0] + shifted_y*theta_matrix[1] + shifted_z*theta_matrix[2];
        point_cloud[i][1] = shifted_x*theta_matrix[3] + shifted_y*theta_matrix[4] + shifted_z*theta_matrix[5];
        point_cloud[i][2] = shifted_x*theta_matrix[6] + shifted_y*theta_matrix[7] + shifted_z*theta_matrix[8];
        if (theta_matrix.size()>=13) {
            point_cloud[i][0] *= theta_matrix[12];
            point_cloud[i][1] *= theta_matrix[13];
            point_cloud[i][2] *= theta_matrix[14];
        }
        point_cloud[i][0] += theta_matrix[9];
        point_cloud[i][1] += theta_matrix[10];
        point_cloud[i][2] += theta_matrix[11];
    }
}

void apply_rot_trans(double roll, double pitch, double yaw, double x_shift, double y_shift, double z_shift, vector<vector<double>>& point_cloud){
    double beta = roll*pi/180;// roll in radians
    double gamma = pitch*pi/180; // pitch in radians
    double alpha = yaw*pi/180; // yaw in radians
    double shifted_x, shifted_y, shifted_z;
    size_t n = point_cloud.size();
    DebugOn(roll<<" "<<pitch<<" "<<yaw<<endl);
    DebugOn(beta<<" "<<gamma<<" "<<alpha<<endl);
    DebugOn(cos(beta)<<endl<<sin(beta)<<endl<<cos(gamma)<<endl<<sin(gamma)<<endl<<cos(alpha)<<endl<<sin(alpha)<<endl);
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
vector<vector<double>> get_n_extreme_points(int max_n, const vector<vector<double>>& point_cloud){
    vector<vector<double>> ext;
    set<int> added_indices;
    while (ext.size()<max_n) {
        double x_min_val = numeric_limits<double>::max(), y_min_val = numeric_limits<double>::max(), z_min_val = numeric_limits<double>::max();
        double x_max_val = numeric_limits<double>::lowest(), y_max_val = numeric_limits<double>::lowest(), z_max_val = numeric_limits<double>::lowest();
        int x_min, y_min, z_min;
        int x_max, y_max, z_max;
        size_t n = point_cloud.size();
        for (auto i = 0; i< n; i++) {
            if(added_indices.count(i)!=0)
                continue;
            if(point_cloud[i][0] < x_min_val){
                x_min_val = point_cloud[i][0];
                x_min = i;
            }
            if(point_cloud[i][0] > x_max_val){
                x_max_val = point_cloud[i][0];
                x_max = i;
            }
            if(point_cloud[i][1] < y_min_val){
                y_min_val = point_cloud[i][1];
                y_min = i;
            }
            if(point_cloud[i][1] > y_max_val){
                y_max_val = point_cloud[i][1];
                y_max = i;
            }
            if(point_cloud[i][2] < z_min_val){
                z_min_val = point_cloud[i][2];
                z_min = i;
            }
            if(point_cloud[i][2] > z_max_val){
                z_max_val = point_cloud[i][2];
                z_max = i;
            }
        }
        if(ext.size()<max_n && added_indices.insert(x_min).second)
            ext.push_back(point_cloud[x_min]);
        if(ext.size()<max_n && added_indices.insert(x_max).second)
            ext.push_back(point_cloud[x_max]);
        if(ext.size()<max_n && added_indices.insert(y_min).second)
            ext.push_back(point_cloud[y_min]);
        if(ext.size()<max_n && added_indices.insert(y_max).second)
            ext.push_back(point_cloud[y_max]);
        if(ext.size()<max_n && added_indices.insert(z_min).second)
            ext.push_back(point_cloud[z_min]);
        if(ext.size()<max_n && added_indices.insert(z_max).second)
            ext.push_back(point_cloud[z_max]);
    }
    return ext;
}


/* Return vector of extreme points from point cloud */
pair<vector<vector<double>>,vector<vector<double>>> get_n_extreme_points(int n, const vector<vector<double>>& point_cloud, const vector<vector<double>>& uav){
    vector<vector<double>> ext, ext_uav;
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
                    ext_uav.push_back(uav[i]);
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
                    ext_uav.push_back(uav[i]);
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
                    ext_uav.push_back(uav[i]);
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
                    ext_uav.push_back(uav[i]);
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
                    ext_uav.push_back(uav[i]);
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
                    ext_uav.push_back(uav[i]);
                    nb++;
                }
                continue;
            }
        }
    }
    return {ext,ext_uav};
}


/* Return vector of extreme points from point cloud */
vector<vector<double>> get_extreme_points(const vector<vector<double>>& point_cloud){
    return get_n_extreme_points(6, point_cloud);
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


/* Compute the L1 error for model and data sets
 @param[in] point_cloud_model, Model point cloud
 @param[in] point_cloud_data, Data point cloud
 @param[in] matching, vecor used to store the optimal match
 @param[in] err_per_point, vector used to store the minimum L1 error per data point
 */
double computeL2error(const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data, vector<int>& matching, vector<double>& err_per_point){
    size_t n = point_cloud_data.size();
    size_t m = point_cloud_model.size();
    double dist_sq = 0, err = 0;
    for (auto i = 0; i< n; i++) {
        double min_dist = numeric_limits<double>::max();
        for (auto j = 0; j< m; j++) {
            dist_sq = std::pow(point_cloud_data.at(i).at(0) - point_cloud_model.at(j).at(0),2) + std::pow(point_cloud_data.at(i).at(1) - point_cloud_model.at(j).at(1),2) + std::pow(point_cloud_data.at(i).at(2) - point_cloud_model.at(j).at(2),2);
            if(min_dist>dist_sq){
                min_dist = dist_sq;
                matching[i] = j;
            }
        }
        DebugOff("error(" << i+1 << ") = " << to_string_with_precision(min_dist,12) << endl);
        DebugOff("matching(" << i+1 << ") = " << matching[i]+1 << endl);
        err_per_point[i] = min_dist;
        err += min_dist;
    }
    return err;
}



/* Compute the L1 error for model and data sets
 @param[in] point_cloud_model, Model point cloud
 @param[in] point_cloud_data, Data point cloud
 @param[in] matching, vecor used to store the optimal match
 @param[in] err_per_point, vector used to store the minimum L1 error per data point
 */
double computeL1error(const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data, vector<int>& matching, vector<double>& err_per_point){
    size_t n = point_cloud_data.size();
    size_t m = point_cloud_model.size();
    double dist_abs = 0, err = 0;
    for (auto i = 0; i< n; i++) {
        double min_dist = numeric_limits<double>::max();
        for (auto j = 0; j< m; j++) {
            dist_abs = std::abs(point_cloud_data.at(i).at(0) - point_cloud_model.at(j).at(0)) + std::abs(point_cloud_data.at(i).at(1) - point_cloud_model.at(j).at(1)) + std::abs(point_cloud_data.at(i).at(2) - point_cloud_model.at(j).at(2));
            if(min_dist>dist_abs){
                min_dist = dist_abs;
                matching[i] = j;
            }
        }
        DebugOn("error(" << i+1 << ") = " << to_string_with_precision(min_dist,12) << endl);
        DebugOn("matching(" << i+1 << ") = " << matching[i]+1 << endl);
        err_per_point[i] = min_dist;
        err += min_dist;
    }
    return err;
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

/* Run ICP on point clouds */
double run_ICP_only(GoICP& goicp, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double shift_min_x, double shift_max_x, double shift_min_y, double shift_max_y, double shift_min_z, double shift_max_z, vector<double>& rot_trans_ub){
    //    var<> yaw("yaw", yaw_min, yaw_max), pitch("pitch", pitch_min, pitch_max), roll("roll", roll_min, roll_max);
    //    func<> r11 = cos(yaw)*cos(roll);
    //    func<> r12 = cos(yaw)*sin(roll)*sin(pitch) - sin(yaw)*cos(pitch);
    //    func<> r13 = cos(yaw)*sin(roll)*cos(pitch) + sin(yaw)*sin(pitch);
    //    func<> r21 = sin(yaw)*cos(roll);
    //    func<> r22 = sin(yaw)*sin(roll)*sin(pitch) + cos(yaw)*cos(pitch);
    //    func<> r23 = sin(yaw)*sin(roll)*cos(pitch) - cos(yaw)*sin(pitch);
    //    func<> r31 = sin(-1*roll);
    //    func<> r32 = cos(roll)*sin(pitch);
    //    func<> r33 = cos(roll)*cos(pitch);
    //
    //
    //    var<> theta11("theta11",  std::max(-1.,r11._range->first), std::min(1.,r11._range->second)), theta12("theta12", std::max(-1.,r12._range->first), std::min(1.,r12._range->second)), theta13("theta13", std::max(-1.,r13._range->first), std::min(1.,r13._range->second));
    //    var<> theta21("theta21", std::max(-1.,r21._range->first), std::min(1.,r21._range->second)), theta22("theta22", std::max(-1.,r22._range->first), std::min(1.,r22._range->second)), theta23("theta23", std::max(-1.,r23._range->first), std::min(1.,r23._range->second));
    //    var<> theta31("theta31", std::max(-1.,r31._range->first), std::min(1.,r31._range->second)), theta32("theta32", std::max(-1.,r32._range->first), std::min(1.,r32._range->second)), theta33("theta33", std::max(-1.,r33._range->first), std::min(1.,r33._range->second));
    //
    //    goicp.R_init.val[0][0]=(std::max(-1.,r11._range->first)+ std::min(1.,r11._range->second))/2.0;
    //    goicp.R_init.val[0][1]=(std::max(-1.,r12._range->first)+ std::min(1.,r12._range->second))/2.0;
    //    goicp.R_init.val[0][2]=(std::max(-1.,r13._range->first)+ std::min(1.,r13._range->second))/2.0;
    //    goicp.R_init.val[1][0]=(std::max(-1.,r21._range->first)+ std::min(1.,r21._range->second))/2.0;
    //    goicp.R_init.val[1][1]=(std::max(-1.,r22._range->first)+ std::min(1.,r22._range->second))/2.0;
    //    goicp.R_init.val[1][2]=(std::max(-1.,r23._range->first)+ std::min(1.,r23._range->second))/2.0;
    //    goicp.R_init.val[2][0]=(std::max(-1.,r31._range->first)+ std::min(1.,r31._range->second))/2.0;
    //    goicp.R_init.val[2][1]=(std::max(-1.,r32._range->first)+ std::min(1.,r32._range->second))/2.0;
    //    goicp.R_init.val[2][2]=(std::max(-1.,r33._range->first)+ std::min(1.,r33._range->second))/2.0;
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
    
    //
    //    var<> theta11("theta11",  std::max(-1.,r11._range->first), std::min(1.,r11._range->second)), theta12("theta12", std::max(-1.,r12._range->first), std::min(1.,r12._range->second)), theta13("theta13", std::max(-1.,r13._range->first), std::min(1.,r13._range->second));
    //    var<> theta21("theta21", std::max(-1.,r21._range->first), std::min(1.,r21._range->second)), theta22("theta22", std::max(-1.,r22._range->first), std::min(1.,r22._range->second)), theta23("theta23", std::max(-1.,r23._range->first), std::min(1.,r23._range->second));
    //    var<> theta31("theta31", std::max(-1.,r31._range->first), std::min(1.,r31._range->second)), theta32("theta32", std::max(-1.,r32._range->first), std::min(1.,r32._range->second)), theta33("theta33", std::max(-1.,r33._range->first), std::min(1.,r33._range->second));
    
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

tuple<double,double,double> run_IPH(vector<vector<double>>& ext_model, vector<vector<double>>& ext_data, const vector<vector<double>>& uav1, const vector<vector<double>>& uav2){
    double roll = 0, pitch = 0, yaw = 1; /* at least one nonzero to enter the while loop */
    double final_roll = 0, final_pitch = 0, final_yaw = 0;
    int nb_iter = 0, max_nb_iter = 100;
    tuple<double,double,double> res;
    while(nb_iter < max_nb_iter && std::abs(roll)+std::abs(pitch)+std::abs(yaw)>1e-1) {
        //        auto L2error = computeL1error(ext_model,ext_data);
        //        DebugOn("L2 error with exterme set before = " << L2error << endl);
        res = run_ARMO("full", ext_model, ext_data, uav1, uav2);
        roll = get<0>(res);pitch = get<1>(res);yaw = get<2>(res);
        final_roll += roll;final_pitch += pitch;final_yaw += yaw;
        apply_rotation(roll, pitch, yaw, ext_model, ext_data, uav1, uav2);
        //        L2error = computeL1error(ext_model,ext_data);
        //        DebugOn("L2 error with exterme set after full = " << L2error << endl);
        nb_iter++;
        DebugOn("No projection, ITERATION " << nb_iter << endl);
    }
    nb_iter = 0;yaw=1;max_nb_iter = 100;
    //    while(nb_iter < max_nb_iter && std::abs(yaw)>1e-1) {
    while(nb_iter < max_nb_iter && std::abs(roll)+std::abs(pitch)+std::abs(yaw)>1e-1) {
        res = run_ARMO("z", ext_model, ext_data, uav1, uav2);
        roll = get<0>(res);pitch = get<1>(res);yaw = get<2>(res);
        final_roll += roll;final_pitch += pitch;final_yaw += yaw;
        apply_rotation(roll, pitch, yaw, ext_model, ext_data, uav1, uav2);
        //        auto L2error = computeL1error(ext_model,ext_data);
        //        DebugOn("L2 error with exterme set after z = " << L2error << endl);
        nb_iter++;
        DebugOn("Projceting out z axis, ITERATION " << nb_iter << endl);
    }
    //    DebugOn("Plotting after z" << endl);
    //
    nb_iter = 0;yaw=1;max_nb_iter = 100;
    //    while(nb_iter < max_nb_iter && std::abs(roll)>1e-1) {
    while(nb_iter < max_nb_iter && std::abs(roll)+std::abs(pitch)+std::abs(yaw)>1e-1) {
        res = run_ARMO("y", ext_model, ext_data, uav1, uav2);
        roll = get<0>(res);pitch = get<1>(res);yaw = get<2>(res);
        final_roll += roll;final_pitch += pitch;final_yaw += yaw;
        apply_rotation(roll, pitch, yaw, ext_model, ext_data, uav1, uav2);
        //        auto L2error = computeL1error(ext_model,ext_data);
        //        DebugOn("L2 error with exterme set after y = " << L2error << endl);
        nb_iter++;
        DebugOn("Projceting out y axis, ITERATION " << nb_iter << endl);
    }
    //    DebugOn("Plotting after y" << endl);
    //    plot(ext_model,ext_data,1);
    nb_iter = 0;yaw=1;max_nb_iter = 100;
    //    while(nb_iter < max_nb_iter && std::abs(pitch)>1e-1) {
    while(nb_iter < max_nb_iter && std::abs(roll)+std::abs(pitch)+std::abs(yaw)>1e-1) {
        res = run_ARMO("x", ext_model, ext_data, uav1, uav2);
        roll = get<0>(res);pitch = get<1>(res);yaw = get<2>(res);
        final_roll += roll;final_pitch += pitch;final_yaw += yaw;
        apply_rotation(roll, pitch, yaw, ext_model, ext_data, uav1, uav2);
        //        auto L2error = computeL1error(ext_model,ext_data);
        //        DebugOn("L2 error with exterme set after x = " << L2error << endl);
        nb_iter++;
        DebugOn("Projceting out x axis, ITERATION " << nb_iter << endl);
    }
    //    DebugOn("Plotting after x" << endl);
    //    plot(ext_model,ext_data,1);
    nb_iter = 0;yaw=1;max_nb_iter = 100;
    while(nb_iter < max_nb_iter && std::abs(roll)+std::abs(pitch)+std::abs(yaw)>1e-1) {
        res = run_ARMO("full", ext_model, ext_data, uav1, uav2);
        roll = get<0>(res);pitch = get<1>(res);yaw = get<2>(res);
        final_roll += roll;final_pitch += pitch;final_yaw += yaw;
        apply_rotation(roll, pitch, yaw, ext_model, ext_data, uav1, uav2);
        //        auto L2error = computeL1error(ext_model,ext_data);
        //        DebugOn("L2 error with exterme set after full = " << L2error << endl);
        nb_iter++;
        DebugOn("No projection, ITERATION " << nb_iter << endl);
    }
    DebugOn("Final Roll (degrees) = " << final_roll << endl);
    DebugOn("Final Pitch (degrees) = " << final_pitch << endl);
    DebugOn("Final Yaw (degrees) = " << final_yaw << endl);
    DebugOn("Final Roll (radians)= " << final_roll *pi/180 << endl);
    DebugOn("Final Pitch (radians) = " << final_pitch *pi/180 << endl);
    DebugOn("Final Yaw (radians) = " << final_yaw *pi/180 << endl);
    return {final_roll,final_pitch,final_yaw};
}


tuple<double,double,double,double,double,double> run_IPH(const vector<vector<double>>& ext_model, vector<vector<double>>& ext_data, vector<vector<double>>& point_cloud_data){
    double roll = 0, pitch = 0, yaw = 0, x_shift = 0, y_shift = 0, z_shift = 1; /* at least one nonzero to enter the while loop */
    double final_roll = 0, final_pitch = 0, final_yaw = 0, final_x_shift = 0, final_y_shift = 0, final_z_shift = 0;
    int nb_iter = 0, max_nb_iter = 100;
    tuple<double,double,double,double,double,double> res;
    while(nb_iter < max_nb_iter && std::abs(roll)+std::abs(pitch)+std::abs(yaw)+std::abs(x_shift)+std::abs(y_shift)+std::abs(z_shift)>1e-1) {
        //        auto L2error = computeL1error(ext_model,ext_data);
        //        DebugOn("L2 error with exterme set before = " << L2error << endl);
        res = run_ARMO(false, "full", ext_model, ext_data);
        roll = get<0>(res);pitch = get<1>(res);yaw = get<2>(res);x_shift = get<3>(res);y_shift = get<4>(res);z_shift = get<5>(res);
        final_roll += roll;final_pitch += pitch;final_yaw += yaw;final_x_shift += x_shift;final_y_shift += y_shift;final_z_shift += z_shift;
        apply_rot_trans(roll, pitch, yaw, x_shift, y_shift, z_shift, ext_data);
        apply_rot_trans(roll, pitch, yaw, x_shift, y_shift, z_shift, point_cloud_data);
        //        L2error = computeL1error(ext_model,ext_data);
        //        DebugOn("L2 error with exterme set after full = " << L2error << endl);
        nb_iter++;
        DebugOn("ITERATION " << nb_iter << endl);
    }
    nb_iter = 0;z_shift=1;max_nb_iter = 100;
    while(nb_iter < max_nb_iter && std::abs(roll)+std::abs(pitch)+std::abs(yaw)+std::abs(x_shift)+std::abs(y_shift)+std::abs(z_shift)>1) {
        res = run_ARMO(false, "z", ext_model, ext_data);
        roll = get<0>(res);pitch = get<1>(res);yaw = get<2>(res);x_shift = get<3>(res);y_shift = get<4>(res);z_shift = get<5>(res);
        final_roll += roll;final_pitch += pitch;final_yaw += yaw;final_x_shift += x_shift;final_y_shift += y_shift;final_z_shift += z_shift;
        apply_rot_trans(roll, pitch, yaw, x_shift, y_shift, z_shift, ext_data);
        apply_rot_trans(roll, pitch, yaw, x_shift, y_shift, z_shift, point_cloud_data);
        //        auto L2error = computeL1error(ext_model,ext_data);
        //        DebugOn("L2 error with exterme set after z = " << L2error << endl);
        nb_iter++;
        DebugOn("ITERATION " << nb_iter << endl);
    }
    //    DebugOn("Plotting after z" << endl);
    //    plot(ext_model,ext_data,1);
    nb_iter = 0;z_shift=1;max_nb_iter = 100;
    while(nb_iter < max_nb_iter && std::abs(roll)+std::abs(pitch)+std::abs(yaw)+std::abs(x_shift)+std::abs(y_shift)+std::abs(z_shift)>1) {
        res = run_ARMO(false, "y", ext_model, ext_data);
        roll = get<0>(res);pitch = get<1>(res);yaw = get<2>(res);x_shift = get<3>(res);y_shift = get<4>(res);z_shift = get<5>(res);
        final_roll += roll;final_pitch += pitch;final_yaw += yaw;final_x_shift += x_shift;final_y_shift += y_shift;final_z_shift += z_shift;
        apply_rot_trans(roll, pitch, yaw, x_shift, y_shift, z_shift, ext_data);
        apply_rot_trans(roll, pitch, yaw, x_shift, y_shift, z_shift, point_cloud_data);
        //        auto L2error = computeL1error(ext_model,ext_data);
        //        DebugOn("L2 error with exterme set after y = " << L2error << endl);
        nb_iter++;
        DebugOn("ITERATION " << nb_iter << endl);
    }
    //    DebugOn("Plotting after y" << endl);
    //    plot(ext_model,ext_data,1);
    nb_iter = 0;z_shift=1;max_nb_iter = 100;
    while(nb_iter < max_nb_iter && std::abs(roll)+std::abs(pitch)+std::abs(yaw)+std::abs(x_shift)+std::abs(y_shift)+std::abs(z_shift)>1) {
        res = run_ARMO(false, "x", ext_model, ext_data);
        roll = get<0>(res);pitch = get<1>(res);yaw = get<2>(res);x_shift = get<3>(res);y_shift = get<4>(res);z_shift = get<5>(res);
        final_roll += roll;final_pitch += pitch;final_yaw += yaw;final_x_shift += x_shift;final_y_shift += y_shift;final_z_shift += z_shift;
        apply_rot_trans(roll, pitch, yaw, x_shift, y_shift, z_shift, ext_data);
        apply_rot_trans(roll, pitch, yaw, x_shift, y_shift, z_shift, point_cloud_data);
        //        auto L2error = computeL1error(ext_model,ext_data);
        //        DebugOn("L2 error with exterme set after x = " << L2error << endl);
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
        //        auto L2error = computeL1error(ext_model,ext_data);
        //        DebugOn("L2 error with exterme set after full = " << L2error << endl);
        nb_iter++;
        DebugOn("ITERATION " << nb_iter << endl);
    }
    Debug("Final Roll (degrees) = " << final_roll << endl);
    Debug("Final Pitch (degrees) = " << final_pitch << endl);
    Debug("Final Yaw (degrees) = " << final_yaw << endl);
    Debug("Final Roll (radians)= " << final_roll *pi/180 << endl);
    Debug("Final Pitch (radians) = " << final_pitch *pi/180 << endl);
    Debug("Final Yaw (radians) = " << final_yaw *pi/180 << endl);
    Debug("Final x shift = " << final_x_shift << endl);
    Debug("Final y shift = " << final_y_shift << endl);
    Debug("Final z shift = " << final_z_shift << endl);
    return {final_roll,final_pitch,final_yaw,final_x_shift,final_y_shift,final_z_shift};
}



/* Read input files */
void read_data(const rapidcsv::Document& Model_doc,vector<vector<double>>& point_cloud, vector<vector<double>>& uav){
    int model_nb_rows = Model_doc.GetRowCount();
    if(model_nb_rows<3){
        throw invalid_argument("Input file with less than 2 points");
    }
    DebugOn("Input file has " << model_nb_rows << " rows" << endl);
    point_cloud.resize(model_nb_rows);
    uav.resize(model_nb_rows);
    for (int i = 0; i< model_nb_rows; i++) {
        auto laser_id = Model_doc.GetCell<int>(0, i);
        auto x = Model_doc.GetCell<double>(1, i);
        auto y = Model_doc.GetCell<double>(2, i);
        auto z = Model_doc.GetCell<double>(3, i);
        auto uav_x = Model_doc.GetCell<double>(4, i);
        auto uav_y = Model_doc.GetCell<double>(5, i);
        auto uav_z = Model_doc.GetCell<double>(6, i);
        point_cloud[i].resize(3);
        point_cloud[i][0] = x;
        point_cloud[i][1] = y;
        point_cloud[i][2] = z;
        uav[i].resize(3);
        uav[i][0] = uav_x;
        uav[i][1] = uav_y;
        uav[i][2] = uav_z;
    }
}

/* Read Laz files */
void read_laz(const string& fname){
    LASreadOpener lasreadopener;
    lasreadopener.set_file_name(fname.c_str());
    lasreadopener.set_populate_header(TRUE);
    param<> x1("x1"), y1("x1"), z1("x1");
    int xdim1=0, ydim1=0, zdim1=0;
    if (!lasreadopener.active())
    {
        throw invalid_argument("ERROR: no input specified\n");
    }
    vector<double> x_vec1,y_vec1,z_vec1,zmin_vec1,zmax_vec1;
    vector<double> x_shift,y_shift,z_shift;
    vector<double> x_combined,y_combined,z_combined,zmin_combined,zmax_combined;
    set<double> timestamps;
    while (lasreadopener.active())
    {
        LASreader* lasreader = lasreadopener.open();
        if (lasreader == 0)
        {
            throw invalid_argument("ERROR: could not open lasreader\n");
        }
        
        DebugOn("Number of points = " << lasreader->npoints << endl);
        DebugOn("min x axis = " << lasreader->header.min_x << endl);
        DebugOn("max x axis = " << lasreader->header.max_x << endl);
        DebugOn("min y axis = " << lasreader->header.min_y << endl);
        DebugOn("max y axis = " << lasreader->header.max_y << endl);
        DebugOn("min z axis = " << lasreader->header.min_z << endl);
        DebugOn("max z axis = " << lasreader->header.max_z << endl);
        
        int nb_dots; /* Number of measurements inside cell */
        int xpos, ypos;
        double z, min_z, max_z, av_z;
        pair<int,int> pos;
        size_t nb_pts = 0;
        tuple<double,double,double,double,UAVPoint*> cell; /* <min_z,max_z,av_z> */
        /* Now get rid of the first points
         lasreader->read_point();
         auto gps_time = lasreader->point.get_gps_time();
         while (lasreader->point.get_gps_time()==gps_time)
         {
         lasreader->read_point();
         }
         */
        
        //        /* Values below are used to identify u-turns in drone flight */
        bool neg_x = false;/* x is decreasing */
        bool neg_y = false;/* y is decreasing */
        //
        vector<UAVPoint*> UAVPoints;
        vector<LidarPoint*> LidarPoints;
        map<int,shared_ptr<Frame>> frames;
        map<int,shared_ptr<Frame>> frames1, frames2;
        vector<double> uav_x, uav_y, uav_z;
        vector<double> uav_x1, uav_y1, uav_z1;
        vector<double> x_vec1,y_vec1,z_vec1,zmin_vec1,zmax_vec1;
        vector<double> x_vec2,y_vec2,z_vec2,zmin_vec2,zmax_vec2;
        vector<double> x_shift1,y_shift1,z_shift1;
        vector<double> x_shift2,y_shift2,z_shift2;
        vector<double> x_shift,y_shift,z_shift;
        vector<double> uav_roll1,uav_pitch1,uav_yaw1;
        vector<double> uav_roll2,uav_pitch2,uav_yaw2;
        vector<double> x_combined,y_combined,z_combined,zmin_combined,zmax_combined;
        set<double> timestamps;
        set<int> xvals;
        size_t uav_id = 0;
        bool new_uav = true, u_turn = false, frame1 = true, u_turn_2=false;
        double unix_time, delta_x = 0, delta_y = 0;
        pair<map<int,shared_ptr<Frame>>::iterator,bool> frame_ptr;
        bool exit = false;
        vector<vector<double>> point_cloud1, point_cloud2;
        while (lasreader->read_point() && LidarPoints.size()!=10e6)
        {
            //            if(nb_pts++<1e4)
            //                continue;
            //            if(nb_pts==0){
            //                DebugOn(to_string_with_precision(10.*(lasreader->point.get_gps_time()+315964800. - 18.),24) << ": (" << to_string_with_precision(lasreader->point.get_x(),10) <<"," << to_string_with_precision(lasreader->point.get_y(),10) << ","<< to_string_with_precision(lasreader->point.get_z(),10) <<")"<<endl);
            //                    //                return 0;
            //            }
            auto laser_id = lasreader->point.get_point_source_ID();
            //            if(laser_id!=15){/* Only keep points from Nadir laser */
            //                continue;
            //            }
            auto unix_time = lasreader->point.get_gps_time();
            auto x = lasreader->point.get_x();
            auto y = lasreader->point.get_y();
            auto z = lasreader->point.get_z();
            LidarPoints.push_back(new LidarPoint(laser_id,unix_time,x,y,z));
            point_cloud1.push_back({x,y,z});
            //            if(!xvals.insert(x*100).second){/* A U turn is being detected */
            //                u_turn = true;
            //                DebugOn("Detected a Uturn at point " << LidarPoints.size() << endl);
            //                if(u_turn) {
            //                    DebugOn("This is the second Uturn! " << endl);
            //                    u_turn_2 = true;
            //                }
            //                frame1 = false;
            //            }
            //            if(frame1){
            //                point_cloud1.push_back({x,y,z});
            //            }
            //            else{
            //                point_cloud2.push_back({x,y,z});
            //            }
        }
        //        plot(point_cloud1, point_cloud2);
        DebugOn("Read " << LidarPoints.size() << " points" << endl);
        DebugOn(point_cloud1.size() << " points in flight line 1" << endl);
        DebugOn(point_cloud2.size() << " points in flight line 2" << endl);
        
        save_laz("flight3.laz", point_cloud1, point_cloud2);
        //                auto frame_id = CSV_data.GetCell<int>(0, i);
        //                new_uav = (uav_id==0) || (UAVPoints[uav_id-1]->_frame_id != frame_id);
        //                if(new_uav){
        //                    auto uav_x1 = CSV_data.GetCell<double>("Track_UTM_E", i);
        //                    auto uav_y1 = CSV_data.GetCell<double>("Track_UTM_N", i);
        //                    if(UAVPoints.size()==2){
        //                        auto uav_x0 = UAVPoints.back()->_x;
        //                        auto uav_y0 = UAVPoints.back()->_y;
        //                        neg_x = (uav_x1 - uav_x0) < 0;/* x is decreasing */
        //                        neg_y = (uav_y1 - uav_y0) < 0;/* y is decreasing */
        //                    }
        //                    else if(UAVPoints.size()>2){
        //                        auto uav_x0 = UAVPoints.back()->_x;
        //                        auto uav_y0 = UAVPoints.back()->_y;
        //                        bool neg_x_new = (uav_x1 - uav_x0) < 0;/* x is decreasing */
        //                        bool neg_y_new = (uav_y1 - uav_y0) < 0;/* y is decreasing */
        //                        if(neg_x_new!=neg_x || neg_y_new!=neg_y){/* A U turn is being detected */
        //                            u_turn = true;
        //                            frame1 = false;
        //                            neg_x = neg_x_new;
        //                            neg_y = neg_y_new;
        //                        }
        //                        else {
        //                            u_turn = false;
        //                        }
        //                    }
        //                    UAVPoints.push_back(new UAVPoint());
        //                    UAVPoints[uav_id]->_frame_id = frame_id;
        //                    UAVPoints[uav_id]->_x = uav_x1;
        //                    UAVPoints[uav_id]->_y = uav_y1;
        //                    UAVPoints[uav_id]->_height = CSV_data.GetCell<double>("Track_UTM_Height", i);
        //                    unix_time = CSV_data.GetCell<double>("Time", i);
        //                    UAVPoints[uav_id]->set_unix_time(unix_time);
        //                    uav_x.push_back(UAVPoints[uav_id]->_x);
        //                    uav_y.push_back(UAVPoints[uav_id]->_y);
        //                    uav_z.push_back(UAVPoints[uav_id]->_height);
        //                    frame_ptr = frames.insert(make_pair(UAVPoints[uav_id]->_frame_id, make_shared<Frame>(UAVPoints[uav_id]->_frame_id, UAVPoints[uav_id]->_unix_time)));
        //                    frame_ptr.first->second->add_UAV_point(UAVPoints[uav_id]);
        //                    if(frame1){/* Has not performed a u-turn yet, keep adding to frames1 */
        //                        frames1.insert(make_pair(frame_ptr.first->second->_id, frame_ptr.first->second));
        //                    }
        //                    else{/* Already turned, keep adding to frames2 */
        //                        frames2.insert(make_pair(frame_ptr.first->second->_id, frame_ptr.first->second));
        //                    }
        //                    if(u_turn){
        //                        DebugOn("Detected a Uturn at frame " << frame_ptr.first->first << endl);
        //                    }
        //                    uav_id++;
        //                }
        //
        //                auto xpos = CSV_data.GetCell<double>("UTM_E", i);
        //                auto ypos = CSV_data.GetCell<double>("UTM_N", i);
        //                auto zpos = CSV_data.GetCell<double>("UTM_Height", i);
        //                LidarPoints.push_back(new LidarPoint(laser_id,unix_time,xpos,ypos,zpos));
        //                frame_ptr.first->second->add_lidar_point(LidarPoints.back());
        //                LidarPoints.back()->_uav_pt = frame_ptr.first->second->_uav_point;
        //
        //                //                uav_x1.push_back((frame_ptr.first->second._uav_points.front()->_longitude+582690.8242)*1e-5);
        //                //                uav_y1.push_back((frame_ptr.first->second._uav_points.front()->_latitude+4107963.58)*1e-5);
        //                //                uav_z1.push_back(frame_ptr.first->second._uav_points.front()->_height*100);
        //            }
        //            DebugOn("Read " << uav_id << " frames" << endl);
        //            DebugOn(frames1.size() << " frames in flight line 1" << endl);
        //            DebugOn(frames2.size() << " frames in flight line 2" << endl);
        //            DebugOn(LidarPoints.size() << " lidar points read" << endl);
        //            int nb_pts_per_frame1 = 0, nb_pts_per_frame2 = 0;
        //            for (const auto &frame: frames1) {
        //                nb_pts_per_frame1 += frame.second->_lidar_points->size();
        //                int i = 0;
        //                for (const auto &p: *frame.second->_lidar_points) {
        //                    if(i%10==0){
        //                        x_vec1.push_back(p->_x);
        //                        x_shift1.push_back(frame.second->_uav_point->_x);
        //                        y_vec1.push_back(p->_y);
        //                        y_shift1.push_back(frame.second->_uav_point->_y);
        //                        z_vec1.push_back(p->_z);
        //                        z_shift1.push_back(frame.second->_uav_point->_height);
        //                    }
        //                    i++;
        //                }
        //
        //            }
        //            for (const auto &frame: frames2) {
        //                nb_pts_per_frame2 += frame.second->_lidar_points->size();
        //                int i = 0;
        //                for (auto const &p: *frame.second->_lidar_points) {
        //                    if(i%10==0){
        //                        x_vec2.push_back(p->_x);
        //                        x_shift2.push_back(frame.second->_uav_point->_x);
        //                        y_vec2.push_back(p->_y);
        //                        y_shift2.push_back(frame.second->_uav_point->_y);
        //                        z_vec2.push_back(p->_z);
        //                        z_shift2.push_back(frame.second->_uav_point->_height);
        //                    }
        //                    i++;
        //                }
        //            }
        //            if(frames1.size()!=0)
        //                DebugOn("Average number of points per frame in flight line 1 = " << nb_pts_per_frame1/frames1.size() << endl);
        //            if(frames2.size()!=0)
        //                DebugOn("Average number of points per frame in flight line 2 = " << nb_pts_per_frame2/frames2.size() << endl);
        //            bool plot_data = false;
        //            if(plot_data){
        //            }
        //        while (lasreader->read_point())
        //        {
        //            if(nb_pts==0){
        //                DebugOn(to_string_with_precision(10.*(lasreader->point.get_gps_time()+315964800. - 18.),24) << ": (" << to_string_with_precision(lasreader->point.get_x(),10) <<"," << to_string_with_precision(lasreader->point.get_y(),10) << ","<< to_string_with_precision(lasreader->point.get_z(),10) <<")"<<endl);
        ////                return 0;
        //            }
        //
        //        }
    }
}

vector<double> projection(vector<double> normal, double intercept, vector<double> point){
    vector<double> res;
    res.resize(3);
    auto xi=point[0];
    auto yi=point[1];
    auto zi=point[2];
    auto a=normal[0]/(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
    auto b=normal[1]/(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
    auto c=normal[2]/(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
    auto d=intercept/(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
    auto dist=(a*xi+b*yi+c*zi+d);
    res[0]=xi-dist*a;
    res[1]=yi-dist*b;
    res[2]=zi-dist*c;
    return res;
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
/* Save LAZ files */
void save_laz(const string& fname, const vector<vector<double>>& point_cloud1, const vector<vector<double>>& point_cloud2){
    DebugOn("Saving new las file\n");
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
    
    auto n1 = point_cloud1.size();
    auto n2 = point_cloud2.size();
    LASwriteOpener laswriteopener;
    laswriteopener.set_file_name(fname.c_str());
    LASwriter* laswriter = laswriteopener.open(&lasheader);
    LASpoint laspoint;
    laspoint.init(&lasheader, lasheader.point_data_format, lasheader.point_data_record_length, 0);
    for (auto i = 0; i< n1; i++) {
        laspoint.set_x(point_cloud1[i][0]*1e2);
        laspoint.set_y(point_cloud1[i][1]*1e2);
        laspoint.set_z(point_cloud1[i][2]*1e2);
        laswriter->write_point(&laspoint);
        laswriter->update_inventory(&laspoint);
    }
    for (auto i = 0; i< n2; i++) {
        laspoint.set_x(point_cloud2[i][0]*1e2);
        laspoint.set_y(point_cloud2[i][1]*1e2);
        laspoint.set_z(point_cloud2[i][2]*1e2);
        laswriter->write_point(&laspoint);
        laswriter->update_inventory(&laspoint);
    }
    laswriter->update_header(&lasheader, TRUE);
    laswriter->close();
    delete laswriter;
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
