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

/* Compute the L1 error */
double computeL1error(const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data, vector<int>& matching, vector<double>& err_per_point);

/* Return central point from point cloud */
vector<double> get_center(const vector<vector<double>>& point_cloud);


/* Apply rotation on input data (Boresight alignment) */
void apply_rotation(double roll, double pitch, double yaw, vector<vector<double>>& point_cloud1, vector<vector<double>>& point_cloud2, const vector<vector<double>>& uav1, const vector<vector<double>>& uav2);

/* Update point cloud coordinates */
void update_xyz(vector<vector<double>>& point_cloud1, vector<double>& x_vec1, vector<double>& y_vec1, vector<double>& z_vec1);

/* Run the ARMO model for boresight alignment */
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

/* Run the Global ARMO model for registration */
tuple<double,double,double,double,double,double> run_ARMO_Global_reform(bool bypass, string axis, const vector<vector<double>>& point_cloud1, const vector<vector<double>>& point_cloud2);

shared_ptr<Model<double>> model_Global_reform(bool bypass, string axis, vector<vector<double>>& point_cloud1, vector<vector<double>>& point_cloud2, vector<double>& rot_trans, bool norm1 = false);

bool get_solution(const shared_ptr<Model<double>>& M, vector<double>& rot_trans, vector<int>& new_matching);
void get_rotation_transl_matrix(const shared_ptr<Model<double>>& M, vector<double>& rot_trans);
void get_angle_rotation_transl_matrix(const shared_ptr<Model<double>>& M, vector<double>& rot_trans);
void update_matching(shared_ptr<Model<double>>& M, vector<int>& new_matching);
void round_bin(shared_ptr<Model<double>>& M, int nd, int nm);

shared_ptr<Model<double>> build_TU_MIP(vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, vector<double>& rot_trans, const vector<pair<pair<int,int>,pair<int,int>>>& incompatibles);


shared_ptr<Model<double>> build_SDP(vector<double>& point, vector<double>& rot_mat);

void update_theta_bounds(shared_ptr<Model<double>>& mod, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max);

void initialize_model_from_parent(shared_ptr<Model<double>>parent, shared_ptr<Model<double>> child);

vector<double> BranchBound(vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, param<>& norm_x, param<>& norm_y, param<>& norm_z,  param<>& intercept, const vector<int>& init_matching, const vector<double>& init_err_per_point, param<>& model_radius, vector<vector<vector<double>>> model_voronoi_normals, vector<vector<double>> model_face_intercept, vector<vector<vector<double>>> model_voronoi_vertices, vector<int>& new_model_pts, indices& new_model_ids, param<>& dist_cost, bool relax_ints, bool relax_sdp, bool rigid_transf);


vector<double> BranchBound2(vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, param<>& norm_x, param<>& norm_y, param<>& norm_z,  param<>& intercept, const vector<int>& init_matching, const vector<double>& init_err_per_point, param<>& model_radius, vector<vector<vector<double>>> model_voronoi_normals, vector<vector<double>> model_face_intercept, vector<vector<vector<double>>> model_voronoi_vertices, vector<int>& new_model_pts, indices& new_model_ids, param<>& dist_cost, bool relax_ints, bool relax_sdp, bool rigid_transf);

vector<double> BranchBound3(vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, param<>& norm_x, param<>& norm_y, param<>& norm_z,  param<>& intercept, const vector<int>& init_matching, const vector<double>& init_err_per_point, param<>& model_radius, vector<vector<vector<double>>> model_voronoi_normals, vector<vector<double>> model_face_intercept, vector<vector<vector<double>>> model_voronoi_vertices, vector<int>& new_model_pts, indices& new_model_ids, param<>& dist_cost, bool relax_ints, bool relax_sdp, bool rigid_transf);

shared_ptr<Model<double>> build_norm2_SOC_MIQCP(vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, const indices& valid_cells, const indices& new_model_ids, const param<>& dist_cost, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double shift_min_x, double shift_max_x, double shift_min_y, double shift_max_y, double shift_min_z, double shift_max_z, vector<double>& rot_trans, bool convex, const vector<pair<pair<int,int>,pair<int,int>>>& incompatibles, param<>& norm_x, param<>& norm_y, param<>& norm_z,  param<>& intercept, const vector<int>& init_matching, const vector<double>& error_per_point, param<>& model_radius, bool relax_ints, bool relax_sdp = false, bool rigid_transf = true, double perc_outliers = 0);

shared_ptr<Model<double>> build_linobj_convex(vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, const indices& valid_cells, double new_roll_min, double new_roll_max, double new_pitch_min, double new_pitch_max, double new_yaw_min, double new_yaw_max, double new_shift_min_x, double new_shift_max_x, double new_shift_min_y, double new_shift_max_y, double new_shift_min_z, double new_shift_max_z, vector<double>& rot_trans, bool separate, const vector<pair<pair<int,int>,pair<int,int>>>& incompatibles,  param<>& norm_x,  param<>& norm_y,  param<>& norm_z,  param<>& intercept,const vector<int>& init_matching, const vector<double>& error_per_point, bool relax_inits, double ub);

shared_ptr<Model<double>> build_linobj_convex_OLD(vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, const indices& valid_cells, double new_roll_min, double new_roll_max, double new_pitch_min, double new_pitch_max, double new_yaw_min, double new_yaw_max, double new_shift_min_x, double new_shift_max_x, double new_shift_min_y, double new_shift_max_y, double new_shift_min_z, double new_shift_max_z, vector<double>& rot_trans, bool separate, const vector<pair<pair<int,int>,pair<int,int>>>& incompatibles,  param<>& norm_x,  param<>& norm_y,  param<>& norm_z,  param<>& intercept,const vector<int>& init_matching, const vector<double>& error_per_point, bool relax_inits);


indices get_valid_pairs(vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double shift_min_x, double shift_max_x, double shift_min_y, double shift_max_y, double shift_min_z, double shift_max_z, param<>& norm_x, param<>& norm_y, param<>& norm_z,  param<>& intercept, const vector<double>& model_voronoi_out_radius, vector<int>& new_model_pts, indices& new_model_ids, bool norm1);

//indices get_valid_pairs(vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double shift_min_x, double shift_max_x, double shift_min_y, double shift_max_y, double shift_min_z, double shift_max_z, vector<vector<vector<double>>> model_voronoi_normals, vector<vector<double>> model_face_intercept, const vector<double>& model_voronoi_out_radius, vector<vector<vector<double>>> model_voronoi_vertices, bool norm1);

indices preprocess_QP(const vector<vector<double>>& point_cloud_data, const vector<vector<double>>& point_cloud_model, const indices& old_cells, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double shift_min_x, double shift_max_x, double shift_min_y, double shift_max_y, double shift_min_z, double shift_max_z, const vector<vector<vector<double>>>& model_voronoi_normals, const vector<vector<double>>& model_face_intercept, const vector<vector<vector<double>>>& model_voronoi_vertices, vector<int>& new_model_pts, indices& new_model_ids, param<>& dist_cost, double upper_bound, int nb_total_threads);


shared_ptr<Model<double>> build_SOC_MIQCP(vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, vector<double>& rot_trans, bool convex, const vector<pair<pair<int,int>,pair<int,int>>>& incompatibles, param<>& norm_x, param<>& norm_y, param<>& norm_z, param<>& intercept, const vector<int>& init_matching);

shared_ptr<Model<double>> build_new_SOC_MIQCP(vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, vector<double>& rot_trans, bool convex, const vector<pair<pair<int,int>,pair<int,int>>>& incompatibles, param<>& norm_x, param<>& norm_y, param<>& norm_z, param<>& intercept, const vector<int>& init_matching);

shared_ptr<Model<double>> build_projected_SOC_MIQCP(vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, vector<double>& rot_trans, bool convex, const vector<pair<pair<int,int>,pair<int,int>>>& incompatibles, param<>& norm_x, param<>& norm_y, param<>& norm_z, param<>& intercept, const vector<int>& init_matching);



shared_ptr<Model<double>> build_norm1_SOC_MIQCP(vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, const indices& valid_cells, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double shift_min_x, double shift_max_x, double shift_min_y, double shift_max_y, double shift_min_z, double shift_max_z, vector<double>& rot_trans, bool convex, const vector<pair<pair<int,int>,pair<int,int>>>& incompatibles, param<>& norm_x, param<>& norm_y, param<>& norm_z,  param<>& intercept, const vector<int>& init_matching, const vector<double>& error_per_point, bool relax_ints);

/* Run the MINLP ARMO model for registration */
tuple<double,double,double,double,double,double> run_ARMO_MINLP(bool bypass, string axis, const vector<vector<double>>& point_cloud1, const vector<vector<double>>& point_cloud2);

shared_ptr<Model<double>> three_point_model(vector<vector<double>> d, vector<vector<double>>m, bool convex);


#ifdef USE_QHULL
vector<pair<double, double>> get_translation_bounds(vector<pair<double, double>>min_max_model, const vector<double>& flat_point_cloud_data);
double  maximum_inscribed_sphere_all_points(vector<vector<double>> point_cloud_data,  Qhull& qt);
#endif

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

indices preprocess(vector<vector<double>> point_cloud_data, vector<vector<double>> point_cloud_model, double angle_max_deg, double shift_min_x,double shift_max_x,double shift_min_y,double shift_max_y,double shift_min_z,double shift_max_z, vector<vector<vector<double>>> model_voronoi_normals, vector<vector<double>> model_face_intercept);
indices preprocess_new(const vector<vector<double>>& point_cloud_data, const vector<vector<double>>& point_cloud_model, const indices& old_cells, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double shift_min_x, double shift_max_x, double shift_min_y, double shift_max_y, double shift_min_z, double shift_max_z, const vector<vector<vector<double>>>& model_voronoi_normals, const vector<vector<double>>& model_face_intercept, const vector<vector<vector<double>>>& model_voronoi_vertices, vector<int>& new_model_pts, indices& new_model_ids, param<>& dist_cost, double upper_bound, int nb_total_threads);

indices preprocess_poltyope_intersect(const vector<vector<double>>& point_cloud_data, const vector<vector<double>>& point_cloud_model, const indices& old_cells, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double shift_min_x, double shift_max_x, double shift_min_y, double shift_max_y, double shift_min_z, double shift_max_z, const vector<vector<vector<double>>>& model_voronoi_normals, const vector<vector<double>>& model_face_intercept, const vector<vector<vector<double>>>& model_voronoi_vertices, vector<int>& new_model_pts, indices& new_model_ids, param<>& dist_cost, double upper_bound, int nb_total_threads);

shared_ptr<Model<double>> build_linear_lowerbound(vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, const indices& valid_cells, double new_roll_min, double new_roll_max, double new_pitch_min, double new_pitch_max, double new_yaw_min, double new_yaw_max, double new_shift_min_x, double new_shift_max_x, double new_shift_min_y, double new_shift_max_y, double new_shift_min_z, double new_shift_max_z, vector<double>& rot_trans, bool separate, const vector<pair<pair<int,int>,pair<int,int>>>& incompatibles,  param<>& norm_x,  param<>& norm_y,  param<>& norm_z,  param<>& intercept,const vector<int>& init_matching, const vector<double>& error_per_point, bool relax_inits);

shared_ptr<Model<double>> build_polyhedral(vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, const indices& valid_cells, double new_roll_min, double new_roll_max, double new_pitch_min, double new_pitch_max, double new_yaw_min, double new_yaw_max, double new_shift_min_x, double new_shift_max_x, double new_shift_min_y, double new_shift_max_y, double new_shift_min_z, double new_shift_max_z, vector<double>& rot_trans, bool separate, const vector<pair<pair<int,int>,pair<int,int>>>& incompatibles,  param<>& norm_x,  param<>& norm_y,  param<>& norm_z,  param<>& intercept,const vector<int>& init_matching, const vector<double>& error_per_point, const vector<vector<vector<double>>>& model_voronoi_normals, const vector<vector<double>>& model_face_intercept,  bool relax_inits);
double build_upperbound(vector<double>& solutionlb,vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, const indices& valid_cells, double new_roll_min, double new_roll_max, double new_pitch_min, double new_pitch_max, double new_yaw_min, double new_yaw_max, double new_shift_min_x, double new_shift_max_x, double new_shift_min_y, double new_shift_max_y, double new_shift_min_z, double new_shift_max_z, vector<double>& rot_trans);
GoICP initialize_ICP_only(const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data);
double run_ICP_only(GoICP& goicp, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double shift_min_x, double shift_max_x, double shift_min_y, double shift_max_y, double shift_min_z, double shift_max_z, vector<double>& rot_trans_ub);
shared_ptr<Model<double>> build_polyhedral_ipopt(vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, const indices& valid_cells, double new_roll_min, double new_roll_max, double new_pitch_min, double new_pitch_max, double new_yaw_min, double new_yaw_max, double new_shift_min_x, double new_shift_max_x, double new_shift_min_y, double new_shift_max_y, double new_shift_min_z, double new_shift_max_z, vector<double>& rot_trans, bool separate, const vector<pair<pair<int,int>,pair<int,int>>>& incompatibles,  param<>& norm_x,  param<>& norm_y,  param<>& norm_z,  param<>& intercept,const vector<int>& init_matching, const vector<double>& error_per_point, const vector<vector<vector<double>>>& model_voronoi_normals, const vector<vector<double>>& model_face_intercept, bool relax_inits);
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
    if(Registration){
        double x_min = -1, x_max = 1, y_min = -1, y_max = 1, z_min = -1, z_max = 1;
        vector<double> x_vec0, y_vec0, z_vec0, x_vec1, y_vec1, z_vec1;
        vector<vector<double>> initial_point_cloud_model, initial_point_cloud_data;
        vector<vector<double>> point_cloud_model, point_cloud_data;
        string Model_file = string(prj_dir)+"/data_sets/LiDAR/toy_model.txt";
        string Data_file = string(prj_dir)+"/data_sets/LiDAR/toy_data.txt";
        string algo = "MP", global_str = "global", convex_str = "nonconvex", reform_str="no", obbt_str="yes", norm_str="norm2";
        if(argc>2){
            Model_file = argv[2];
        }
        if(argc>3){
            Data_file = argv[3];
        }
        if(argc>4){
            algo = argv[4];
        }
        if(argc>5){
            global_str = argv[5];
        }
        if(argc>6){
            convex_str = argv[6];
        }
        if(argc>7){
            reform_str = argv[7];
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
        vector<double> flat_point_cloud_data, flat_point_cloud_model;
        DebugOn("Model file has " << model_nb_rows << " rows" << endl);
        DebugOn("Data file has " << data_nb_rows << " rows" << endl);
        int row0 = 0;
        //point_cloud_model.resize(model_nb_rows);
        int fwdm=1;
        int fwdd=1;
        bool downsample=true;
        if(downsample && data_nb_rows>10){
            fwdm=4;
            fwdd=2;
            
        }
        for (int i = row0; i< model_nb_rows; i++) { // Input iterator
            auto x = Model_doc.GetCell<double>(0, i);
            auto y = Model_doc.GetCell<double>(1, i);
            auto z = Model_doc.GetCell<double>(2, i);
            vector<double> xyz;
            xyz.push_back(x);
            xyz.push_back(y);
            xyz.push_back(z);
            initial_point_cloud_model.push_back(xyz);
            if(i%fwdm!=0)
                continue;
            x_vec0.push_back(x);
            y_vec0.push_back(y);
            z_vec0.push_back(z);
            //            model_con.put(i, x, y, z);
            flat_point_cloud_model.push_back(x);
            flat_point_cloud_model.push_back(y);
            flat_point_cloud_model.push_back(z);
            //            point_cloud_model[model_nb_rows-i-1].resize(3);
            //            point_cloud_model[model_nb_rows-i-1][0] = x;
            //            point_cloud_model[model_nb_rows-i-1][1] = y;
            //            point_cloud_model[model_nb_rows-i-1][2] = z;
            
            
            point_cloud_model.push_back(xyz);
            
        }
        //point_cloud_data.resize(data_nb_rows);
        for (int i = row0; i< data_nb_rows; i++) { // Input iterator
            if(i%fwdm!=0)
                continue;
            auto x = Data_doc.GetCell<double>(0, i);
            auto y = Data_doc.GetCell<double>(1, i);
            auto z = Data_doc.GetCell<double>(2, i);
            vector<double> xyz;
            xyz.push_back(x);
            xyz.push_back(y);
            xyz.push_back(z);
            initial_point_cloud_data.push_back(xyz);
            x_vec1.push_back(x);
            y_vec1.push_back(y);
            z_vec1.push_back(z);
            flat_point_cloud_data.push_back(x);
            flat_point_cloud_data.push_back(y);
            flat_point_cloud_data.push_back(z);
            point_cloud_data.push_back(xyz);
        }
        data_nb_rows=point_cloud_data.size();
        model_nb_rows=point_cloud_model.size();
        auto min_max_data=center_point_cloud(point_cloud_data);
        auto min_max_model=center_point_cloud(point_cloud_model);
        vector<double> pcd,pcm;
        for (int i = 0; i< data_nb_rows; i++) { // Input iterator
            pcd.push_back(point_cloud_data[i][0]);
            pcd.push_back(point_cloud_data[i][1]);
            pcd.push_back(point_cloud_data[i][2]);
        }
        vector<pair<double, double>> min_max_t;
#ifdef USE_QHULL
        min_max_t=get_translation_bounds(min_max_model, pcd);
#endif
        
        
        int reduced_nb_data = 50;
        int reduced_nb_model = 200;
        bool subsample = false;
        if (subsample) {
            model_nb_rows = reduced_nb_model;
            data_nb_rows = reduced_nb_data;
            point_cloud_model = get_n_extreme_points(reduced_nb_model, point_cloud_model);
            point_cloud_data = get_n_extreme_points(reduced_nb_data, point_cloud_data);
        }
        
#ifdef USE_VORO
        container model_con(x_min,x_max,y_min,y_max,z_min,z_max,10,10,10,false,false,false,8);
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
        vector<vector<vector<double>>> model_voronoi_vertices(model_nb_rows);/* Store the normal vector of each facet of the voronoi cell of each point */
        vector<vector<vector<double>>> model_face_pts(model_nb_rows);/* Store a point from each facet of the voronoi cell of each point */
        vector<vector<double>> model_face_intercept(model_nb_rows);/* Store the constant part (intercept) in the equation of the voronoi cell of each point */
        vector<double> model_voronoi_in_radius(model_nb_rows);/* Store the radius of the largest ball contained IN the voronoi cell of each point */
        vector<double> model_voronoi_out_radius(model_nb_rows);/* Store the radius of the smallest ball enclosing the voronoi cell of each point */
        param<> norm_x("norm_x"), norm_y("norm_y"), norm_z("norm_z"), intercept("intercept");
        param<> model_radius("model_radius");
        indices m_facets("m_facets");
        vector<pair<double,int>> volume(model_nb_rows);/*volume of each voronoi cell <volume,point_id>*/
        int total_nb_faces = 0;
        if(cl.start()) do if(model_con.compute_cell(c,cl)) {
            cl.pos(x,y,z);
            idx=cl.pid();
            c.vertices(x,y,z,v);
            int nb_vertices = v.size()/3;
            volume[idx] = {c.volume(),idx};
            int v_idx = 0;
            double max_dist = 0;
            model_voronoi_vertices[idx].resize(nb_vertices);
            vector<double> vertices;
            vertices.resize(3);
            vector<int> v_order;
            Debug("Cell has " << c.number_of_edges() << " edges \n");
            for (int i = 0; i<nb_vertices; i++) {
                double dist_sq = std::pow(x - v[v_idx],2) + std::pow(y - v[v_idx+1],2) + std::pow(z - v[v_idx+2],2);
                if(dist_sq>max_dist)
                    max_dist = dist_sq;
                vertices[0]=v[v_idx];
                vertices[1]=v[v_idx+1];
                vertices[2]=v[v_idx+2];
                model_voronoi_vertices[idx][i]=vertices;
                v_idx += 3;
            }
            model_voronoi_out_radius[idx] = std::sqrt(max_dist);/* the radius of the smallest ball enclosing the voronoi cell is the distance from the center to the farthest vertex */
            model_radius.add_val(to_string(idx+1), model_voronoi_out_radius[idx]);
            vector
            <double> normals;
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
            v.clear();
            face_vert.clear();
        } while (cl.inc());
        DebugOn("The total number of faces = " << total_nb_faces << endl);
        
        sort(volume.begin(), volume.end(), pts_compare);
        int nb_reduced_model = std::min(model_nb_rows,150);
        vector<vector<double>> red_point_cloud_model(nb_reduced_model);
        for (int i =0; i<nb_reduced_model; i++) {
            red_point_cloud_model[i] = point_cloud_model[volume[i].second];
        }
#ifdef USE_MATPLOT
        //                  plot(red_point_cloud_model,point_cloud_model,2);
#endif
        // Output the particle positions in gnuplot format
        //        model_con.draw_particles("my_points_p.gnu");
        // Output the Voronoi cells in gnuplot format
        //        model_con.draw_cells_gnuplot("my_points_v.gnu");
        
        
        //#ifdef USE_MATPLOT
        //        bool plot_voronoi = false;
        //        if(plot_voronoi){
        //            vector<vector<double>> data_voronoiVertices, test;
        //            test.push_back(point_cloud_data.back());
        //            auto voronoi_faces = data_voronoi_normals.back();
        //            for (int i = 0; i<voronoi_faces.size(); i++) {
        //                data_voronoiVertices.push_back({voronoi_faces[i][0],voronoi_faces[i][1],voronoi_faces[i][3]});
        //            }
        //            plot(test,data_voronoiVertices, 1);
        //        }
        //#endif
        //auto cells=preprocess(point_cloud_data, point_cloud_model, 25, 0.23,0.24,-0.24,-0.23,-0.02,-0.01,  model_voronoi_normals, model_face_intercept);
        double shift_min_xa = 0.23, shift_max_xa = 0.24, shift_min_ya = -0.24,shift_max_ya= -0.23,shift_min_za = -0.02,shift_max_za = -0.01;
        double yaw_mina = -25*pi/180., yaw_maxa = 25*pi/180., pitch_mina = -25*pi/180.,pitch_maxa = 25.*pi/180.,roll_mina = -25*pi/180.,roll_maxa = 25*pi/180.;
#endif
        
        auto old_point_cloud = point_cloud_data;
        
        bool global = global_str=="global";
        bool norm1 = norm_str=="norm1";
        bool obbt=obbt_str=="yes";
        auto ext_model = point_cloud_model;
        auto ext_data = point_cloud_data;
#ifdef USE_MATPLOT
        //                   plot(ext_model,ext_data,1);
#endif
        //        ext_data = get_n_extreme_points(nb_ext, point_cloud_data);
        vector<double> x_vec_model(ext_model.size()), y_vec_model(ext_model.size()), z_vec_model(ext_model.size());
        vector<double> x_vec_data(ext_data.size()), y_vec_data(ext_data.size()), z_vec_data(ext_data.size());
        tuple<double,double,double,double,double,double,double> res_icp;
        tuple<double,double,double,double,double,double> res, res1, res2;
        vector<int> L2matching(data_nb_rows), L1matching(data_nb_rows);
        vector<double> L2err_per_point(data_nb_rows), L1err_per_point(data_nb_rows);
        auto L2error_init = computeL2error(ext_model,ext_data,L2matching,L2err_per_point);
        DebugOn("L2 on reduced set = " << L2error_init << endl);
        // run_ARMO_Global(false, "full", ext_model, ext_data);
#ifdef USE_QHULL
        Qhull qt_data;//("data", 4, point_cloud_data.size(), flat_point_cloud_data.data(), "");
        qt_data.runQhull("data", 3, point_cloud_data.size(), flat_point_cloud_data.data(), "");
        vector<vector<double>> data_facet_centers;
        for(auto it = qt_data.facetList().begin();it!=qt_data.facetList().end();it++){
            auto v = it->getCenter().toStdVector();
            data_facet_centers.push_back(v);
        }
#ifdef USE_MATPLOT
        // plot(ext_data,data_facet_centers, 2);
#endif
        DebugOn("Computing Voronoi regions for data point cloud\n");
        vector<vector<double>> data_voronoiVertices;
        vector<vector<int>> data_voronoiRegions;
        compute_voronoi(qt_data, data_voronoiVertices, data_voronoiRegions);
#ifdef USE_MATPLOT
        // plot(ext_data,data_voronoiVertices, 2);
#endif
        Qhull qt_model("model", 4, point_cloud_model.size(), flat_point_cloud_model.data(), "");
        DebugOn("Computing Voronoi regions for model point cloud\n");
        vector<vector<double>> model_voronoiVertices;
        vector<vector<int>> model_voronoiRegions;
        compute_voronoi(qt_model, model_voronoiVertices, model_voronoiRegions);
#endif
        bool run_goICP = (algo=="GoICP");
        bool IPH = (algo=="IPH");
        bool ARMO = (algo=="ARMO");
        
        if(ARMO || IPH || run_goICP){
            if(IPH){/* Run run_IPH */
                res1 = run_IPH(ext_model,ext_data,point_cloud_data);
                auto roll = get<0>(res1);auto pitch = get<1>(res1);auto yaw = get<2>(res1);auto x_shift = get<3>(res1);auto y_shift = get<4>(res1);auto z_shift = get<5>(res1);
                apply_rot_trans(roll, pitch, yaw, x_shift, y_shift, z_shift, point_cloud_data);
            }
            else if(run_goICP){/* Run GoICP inline */
                res_icp = run_GoICP(ext_model, ext_data);
                auto roll = get<0>(res_icp);auto pitch = get<1>(res_icp);auto yaw = get<2>(res_icp);auto x_shift = get<3>(res_icp);auto y_shift = get<4>(res_icp);auto z_shift = get<5>(res_icp);
                apply_rot_trans(roll, pitch, yaw, x_shift, y_shift, z_shift, point_cloud_data);
            }
            else {
                if(global){
                    bool convex = convex_str=="convex";
                    bool reform = reform_str=="yes";
                    if(reform){
                        res=run_ARMO_Global_reform(convex, "full", ext_model, ext_data);
                        auto roll = get<0>(res);auto pitch = get<1>(res);auto yaw = get<2>(res);auto x_shift = get<3>(res);auto y_shift = get<4>(res);auto z_shift = get<5>(res);
                        apply_rot_trans(roll, pitch, yaw, x_shift, y_shift, z_shift, point_cloud_data);
                    }
                    else{
                        auto res = run_ARMO_Global(convex, "full", ext_model, ext_data, norm1);
                        if(!norm1){
                            auto roll = res[0];auto pitch = res[1];auto yaw = res[2];auto x_shift = res[3];auto y_shift = res[4];auto z_shift = res[5];
                            apply_rot_trans(roll, pitch, yaw, x_shift, y_shift, z_shift, point_cloud_data);
                        }
                        else{
                            apply_rot_trans(res, point_cloud_data);
                        }
                    }
                }
                else{
                    res1 = run_IPH(ext_model,ext_data,point_cloud_data);
                }
            }
        }
        else{
            
            //res_icp = run_GoICP(ext_model, ext_data);
            // auto roll = get<0>(res_icp);auto pitch = get<1>(res_icp);auto yaw = get<2>(res_icp);auto x_shift = get<3>(res_icp);auto y_shift = get<4>(res_icp);auto z_shift = get<5>(res_icp);
            //auto upper_bound=get<6>(res_icp);
            vector<double> rot_trans;
            bool norm1 = true;
            if(norm1){
                rot_trans.resize(12,0.0);
            }
            else{
                rot_trans.resize(6,0.0);
            }
            //            auto Reg_nc=model_Global_reform(false, "full", point_cloud_model, point_cloud_data, rot_trans, norm1);
            bool separate=true;
            bool linearize=false;
            //input
            double shift_min_x =  -0.1, shift_max_x = 0.1, shift_min_y = -0.1,shift_max_y = 0.1,shift_min_z = -0.1,shift_max_z = 0.1;
            double yaw_min = -10*pi/180., yaw_max = 10*pi/180., pitch_min =-10*pi/180.,pitch_max = 10*pi/180.,roll_min =-10*pi/180.,roll_max = 10*pi/180.;
            
            double roll_mid = (roll_max + roll_min)/2.;
            double pitch_mid = (pitch_max + pitch_min)/2.;
            double yaw_mid = (yaw_max + yaw_min)/2.;
            double shift_mid_x = (shift_max_x + shift_min_x)/2.;
            double shift_mid_y = (shift_max_y + shift_min_y)/2.;
            double shift_mid_z = (shift_max_z + shift_min_z)/2.;
            bool convex = false;
            double new_roll_min = roll_min;
            double new_roll_max = roll_max;
            double new_pitch_min = pitch_min;
            double new_pitch_max = pitch_max;
            double new_yaw_min = yaw_min;
            double new_yaw_max = yaw_max;
            double new_shift_min_x = shift_min_x;
            double new_shift_max_x = shift_max_x;
            double new_shift_min_y = shift_min_y;
            double new_shift_max_y = shift_max_y;
            double new_shift_min_z = shift_min_z;
            double new_shift_max_z = shift_max_z;
            bool go_to_midpoint = false;
            if(go_to_midpoint){
                apply_rot_trans(roll_mid, pitch_mid, yaw_mid, shift_mid_x, shift_mid_y, shift_mid_z, point_cloud_data);
                auto L2error_init = computeL2error(point_cloud_model,point_cloud_data,L2matching,L2err_per_point);
                auto L1error_init = computeL1error(point_cloud_model,point_cloud_data,L1matching,L1err_per_point);
                DebugOn("Initial L2 with midpoint = " << L2error_init << endl);
                DebugOn("Initial L1 with midpoint = " << L1error_init << endl);
                
                new_roll_min = -(roll_max - roll_min)/2.;
                new_roll_max = (roll_max - roll_min)/2.;
                new_pitch_min = -(pitch_max - pitch_min)/2.;
                new_pitch_max = (pitch_max - pitch_min)/2.;
                new_yaw_min = -(yaw_max - yaw_min)/2.;
                new_yaw_max = (yaw_max - yaw_min)/2.;
                new_shift_min_x = -(shift_max_x - shift_min_x)/2.;
                new_shift_max_x = (shift_max_x - shift_min_x)/2.;
                new_shift_min_y = -(shift_max_y - shift_min_y)/2.;
                new_shift_max_y = (shift_max_y - shift_min_y)/2.;
                new_shift_min_z = -(shift_max_z - shift_min_z)/2.;
                new_shift_max_z = (shift_max_z - shift_min_z)/2.;
                
                
            }
#ifdef USE_MATPLOT
            //            plot(point_cloud_model,point_cloud_data,1);
#endif
            
            //            auto TU_MIP = build_TU_MIP(point_cloud_model, point_cloud_data, rot_trans, incompatibles);
            
            // auto NC_SOC_MIQCP = build_projected_SOC_MIQCP(point_cloud_model, point_cloud_data, rot_trans, convex, incompatibles, norm_x, norm_y, norm_z, intercept, matching);
            
            
            
            double time_start = get_wall_time();
            vector<int> new_model_pts;
            param<> dist_cost("dist_cost");
            indices new_model_ids("new_model_ids");
            double upper_bound=L2error_init;
            int nb_total_threads=8;
            indices N1 = range(1,point_cloud_data.size());
            indices N2 = range(1,point_cloud_model.size());
            indices valid_cells("valid_cells");
            
            bool use_features = false;
            if(use_features){
#ifdef USE_PCL
                auto model_features = compute_features(point_cloud_model);
                auto data_features = compute_features(point_cloud_data);
                
                pcl::search::KdTree<pcl::FPFHSignature33>::Ptr hist_kdtree (new pcl::search::KdTree<pcl::FPFHSignature33> ());
                hist_kdtree->setInputCloud(model_features.second);
                // K nearest neighbor search
                
                int K = 50;
                
                std::vector<int> pointIdxNKNSearch(K);
                std::vector<float> pointNKNSquaredDistance(K);
                map<int,set<int>> valid_cells_map;
                for (int i = 0; i<N1.size(); i++) {
                    
                    Debug("K nearest neighbor search at features of data point " << i << endl);
                    
                    if ( hist_kdtree->nearestKSearch (data_features.second->at(i), K, pointIdxNKNSearch, pointNKNSquaredDistance) > 0 )
                    {
                        for (std::size_t j = 0; j < pointIdxNKNSearch.size (); ++j){
                            Debug("Neighbor id = " << pointIdxNKNSearch[j] << " (squared distance: " << pointNKNSquaredDistance[j] << ")" << std::endl);
                            valid_cells_map[i].insert(pointIdxNKNSearch[j]);
                        }
                    }
                }
                
                for (const auto &cells: valid_cells_map) {
                    for (const auto m_id: cells.second) {
                        valid_cells.insert(to_string(cells.first+1)+","+to_string(m_id+1));
                    }
                }
                valid_cells.print();
#endif
            }
            else {
                valid_cells = indices(N1,N2);
            }
            bool preprocess = true;
            
            //            valid_cells = get_valid_pairs(point_cloud_data, point_cloud_model, new_roll_min, new_roll_max, new_pitch_min, new_pitch_max, new_yaw_min, new_yaw_max, new_shift_min_x, new_shift_max_x, new_shift_min_y, new_shift_max_y, new_shift_min_z, new_shift_max_z, norm_x, norm_y, norm_z, intercept, model_voronoi_out_radius, new_model_pts, new_model_ids,false);
            /* Build the index set of incompatible pairs */
            vector<pair<pair<int,int>,pair<int,int>>> incompatibles;
            if(preprocess){
            
                valid_cells= preprocess_poltyope_intersect(point_cloud_data,  point_cloud_model, valid_cells, roll_min, roll_max, pitch_min, pitch_max,yaw_min, yaw_max, shift_min_x,  shift_max_x,shift_min_y, shift_max_y, shift_min_z,  shift_max_z,  model_voronoi_normals,  model_face_intercept,  model_voronoi_vertices, new_model_pts,  new_model_ids, dist_cost, upper_bound, nb_total_threads);
                valid_cells.print();
                
                double time_end = get_wall_time();
                auto prep_time = time_end - time_start;
                /* Terminal output */
                DebugOn("Preprocessing time = " << prep_time << endl);
                
                
                int reduced_nm = new_model_pts.size();
                /* Compute pair-wise distance between model points */
                vector<vector<double>> model_pairwise_dist(reduced_nm-1);
                /* Compute pair-wise min and max distance between Voronoi cells of model point pairs */
                vector<vector<double>> voronoi_min_dist(reduced_nm-1), voronoi_max_dist(reduced_nm-1);
                
                for (int i = 0; i<reduced_nm-1; i++) {
                    auto pi = point_cloud_model[new_model_pts[i]];
                    model_pairwise_dist[i].resize(reduced_nm-i-1);
                    voronoi_min_dist[i].resize(reduced_nm-i-1);
                    voronoi_max_dist[i].resize(reduced_nm-i-1);
                    int pj_id = 0;
                    for (int j = i+1; j<reduced_nm; j++) {
                        double min_dist = numeric_limits<double>::max(), max_dist = 0;
                        auto pj = point_cloud_model[new_model_pts[j]];
                        model_pairwise_dist[i][pj_id] = std::sqrt(std::pow(pi[0] - pj[0],2) + std::pow(pi[1] - pj[1],2) + std::pow(pi[2] - pj[2],2));
                        voronoi_min_dist[i][pj_id] = model_pairwise_dist[i][pj_id] - model_voronoi_out_radius[i] - model_voronoi_out_radius[j];
                        voronoi_max_dist[i][pj_id] = model_pairwise_dist[i][pj_id] + model_voronoi_out_radius[i] + model_voronoi_out_radius[j];
                        pj_id++;
                    }
                }
                
                /* Compute pair-wise distance between data points */
                vector<vector<double>> data_pairwise_dist(data_nb_rows-1);
                for (int i = 0; i<data_nb_rows-1; i++) {
                    auto pi = point_cloud_data[i];
                    data_pairwise_dist[i].resize(data_nb_rows-i-1);
                    int pj_id = 0;
                    for (int j = i+1; j<data_nb_rows; j++) {
                        auto pj = point_cloud_data[j];
                        data_pairwise_dist[i][pj_id] = std::sqrt(std::pow(pi[0] - pj[0],2) + std::pow(pi[1] - pj[1],2) + std::pow(pi[2] - pj[2],2));
                        pj_id++;
                    }
                }
                
                int nb_pairs = 0, nb_pairs_max = 10000;
                
                for (int i = 0; i<data_nb_rows-1; i++) {
                    int pj_id = 0;
                    for (int j = i+1; j<data_nb_rows; j++) {
                        double dp = data_pairwise_dist[i][pj_id];
                        for (int k = 0; k<reduced_nm-1; k++) {
                            int pk_id = 0;
                            for (int l = k+1; l<reduced_nm; l++) {
                                double dv_min = voronoi_min_dist[k][pk_id];
                                double dv_max = voronoi_max_dist[k][pk_id];
                                if(dp < dv_min || dp > dv_max){
                                    incompatibles.push_back({{i,j},{new_model_pts[k],new_model_pts[l]}});
                                    Debug("Imcompatible pairs: data(" << i << "," << j << ") with model(" << k << "," << l << ")" << endl);
                                    nb_pairs++;
                                }
                                pk_id++;
                            }
                        }
                        pj_id++;
                    }
                }
                DebugOn("number of incompatible pairs = " << incompatibles.size() << endl);
                if(incompatibles.size() > nb_pairs_max)
                    incompatibles.resize(nb_pairs_max);
            }
        
            bool relax_integers = false, relax_sdp = false, rigid_transf = true;
            vector<int> init_matching;
            vector<double> error_per_point;
//            GoICP go= initialize_ICP_only(point_cloud_model, point_cloud_data);
//            vector<double> rt(12);
//            auto err=run_ICP_only(go, new_roll_min, new_roll_max, new_pitch_min, new_pitch_max, new_yaw_min, new_yaw_max, new_shift_min_x, new_shift_max_x, new_shift_min_y, new_shift_max_y, new_shift_min_z, new_shift_max_z, rt);
//            go.Clear();
            auto R=build_polyhedral_ipopt(point_cloud_model, point_cloud_data, valid_cells, new_roll_min, new_roll_max,new_pitch_min,new_pitch_max,  new_yaw_min, new_yaw_max,  new_shift_min_x, new_shift_max_x, new_shift_min_y, new_shift_max_y, new_shift_min_z, new_shift_max_z,  rot_trans, false,  incompatibles,  norm_x, norm_y,  norm_z,  intercept, init_matching,  error_per_point,  model_voronoi_normals,  model_face_intercept,  relax_integers);
                        solver<> S(R,ipopt);
//            S.use_callback();
//           // R->replace_integers();
                        S.run(5,1e-6);
            R->print();
//
            R->print_solution();
            rot_trans = BranchBound3(point_cloud_model, point_cloud_data, norm_x, norm_y, norm_z, intercept, L2matching, L2err_per_point, model_radius, model_voronoi_normals, model_face_intercept, model_voronoi_vertices, new_model_pts, new_model_ids, dist_cost, relax_integers, relax_sdp, rigid_transf);
                      //  auto NC_SOC_MIQCP = build_norm2_SOC_MIQCP(point_cloud_model, point_cloud_data, valid_cells, new_model_ids, dist_cost, new_roll_min, new_roll_max,  new_pitch_min, new_pitch_max, new_yaw_min, new_yaw_max, new_shift_min_x, new_shift_max_x, new_shift_min_y, new_shift_max_y, new_shift_min_z, new_shift_max_z, rot_trans, convex, incompatibles, norm_x, norm_y, norm_z, intercept, L2matching, L2err_per_point, model_radius, relax_integers, relax_sdp, rigid_transf);
//                        solver<> S1(NC_SOC_MIQCP,gurobi);
//                       S1.use_callback();
//                        S1.run(5,1e-6,9000,100000);
            // get_solution(NC_SOC_MIQCP, rot_trans, L2matching);
            double ub_solver_tol=1e-6, lb_solver_tol=1e-8, range_tol=1e-3, opt_rel_tol=1e-2, opt_abs_tol=1e6;
            unsigned max_iter=1e3, max_time=3000;
            int nb_threads=1;
            SolverType ub_solver_type = ipopt, lb_solver_type = ipopt;
            //            NC_SOC_MIQCP->replace_integers();
            //            SOC_MIQCP->replace_integers();
            //            auto res=NC_SOC_MIQCP->run_obbt(SOC_MIQCP, max_time, max_iter, opt_rel_tol, opt_abs_tol, nb_threads, ub_solver_type, lb_solver_type, ub_solver_tol, lb_solver_tol, range_tol);
            //            vector<int> new_matching(point_cloud_model.size());
            //            vector<int> matching(point_cloud_model.size());
            //
             // auto SOC_MIP = build_linobj_convex(point_cloud_model, point_cloud_data, valid_cells, roll_min, roll_max, pitch_min, pitch_max, yaw_min, yaw_max, shift_min_x, shift_max_x, shift_min_y, shift_max_y, shift_min_z, shift_max_z, rot_trans, separate=false, incompatibles, norm_x, norm_y, norm_z, intercept,L2matching, L2err_per_point, false, 12.0);
            
            
            //auto SOC_MIPlin = build_polyhedral(point_cloud_model, point_cloud_data, valid_cells, roll_min, roll_max, pitch_min, pitch_max, yaw_min, yaw_max, shift_min_x, shift_max_x, shift_min_y, shift_max_y, shift_min_z, shift_max_z, rot_trans, separate=false, incompatibles, norm_x, norm_y, norm_z, intercept,L2matching, L2err_per_point,  model_voronoi_normals,  model_face_intercept, false);
//            apply_rot_trans(rot_trans, point_cloud_data);
//
//            apply_rot_trans(rot_trans, point_cloud_data);
//
//
//            apply_rot_trans(rot_trans, initial_point_cloud_data);
            
            //SOC_MIP->print();
            //  solver<> S(SOC_MIP,gurobi);
            // S.run();
            //            SOC_MIP->print_int_solution();
            //            SOC_MIP->print_solution();
            
            //            vector<int> new_matching(point_cloud_data.size());
            //            auto SOC_MIP = build_linobj_convex(ext_model, ext_data, rot_trans,new_matching, separate=true, incompatibles, norm_x, norm_y, norm_z, intercept,min_max_model, L2matching,L2err_per_point,true);
            ////            SOC_MIP->print();
            //            solver<> S(SOC_MIP,ipopt);
            //            S.run();
            //            SOC_MIP->print_int_solution();
            //            SOC_MIP->print_solution();
            //            bool is_rotation = get_solution(SOC_MIP, rot_trans,new_matching);
            //            apply_rot_trans(rot_trans, point_cloud_data);
            //            auto L2error_it1 = computeL2error(point_cloud_model,point_cloud_data,matching);
            //            DebugOn("L2 at root node = " << to_string_with_precision(L2error_it1,12) << endl);
            //            SOC_MIP->round_and_fix();
            //            SOC_MIP->reset_constrs();
            //            S.run();
            //            SOC_MIP->print_solution();
            //            is_rotation = get_solution(SOC_MIP, rot_trans,new_matching);
            //            SOC_MIP->print();
            //            vector<double>x(SOC_MIP->get_nb_vars());
            //            auto c = SOC_MIP->get_var<double>("c");
            //            auto bin = SOC_MIP->get_var<int>("bin");
            //            for (int i = 0; i<ext_data.size(); i++) {
            //                for (int j = 0; j<ext_model.size(); j++) {
            //                    string key = to_string(i+1)+","+to_string(j+1);
            //                    c.param<double>::set_val(key, c.eval(key)*bin.eval(key));
            //                }
            //            }
            //            SOC_MIP->get_solution(x);
            //            SOC_MIP->print_solution();
            //            auto SOC_MIP_sep = build_linobj_convex(ext_model, ext_data, rot_trans,separate=true, incompatibles, norm_x, norm_y, norm_z, intercept,min_max_model, matching, true);
            //            SOC_MIP_sep->print();
            //            SOC_MIP_sep->set_solution(x);
            //            bool feas = SOC_MIP_sep->is_feasible(1e-6);
            
            ////            SOC_MIP->print();
            //            if(linearize){
            //                int constr_viol=1;
            //                Model<> interior_model;
            //                auto lin=SOC_MIP->buildOA();
            //                interior_model=lin->add_outer_app_solution(*SOC_MIP);
            //                    //                auto SOC_MIP = build_linobj_convex(point_cloud_model, point_cloud_data, rot_trans,separate);
            //                    //            auto lin=SOC_MIP->outer_approximate_continuous_relaxation(10,constr_viol);
            //                    // lin->print();
            //                constr_viol=1;
            //                int oacuts=0;
            //                vector<double> solution(lin->_nb_vars);
            //                int nb_count=0;
            //                double obj_old=-15, obj_new=-15;
            //                while(constr_viol==1){
            //                    solver<> S(lin,gurobi);
            //                    S.run();
            //                    lin->print_int_solution();
            //                    lin->get_solution(solution);
            //                    obj_new=lin->get_obj_val();
            //                    if(obj_new<=obj_old){
            //                        break;
            //                    }
            //
            //                    constr_viol=SOC_MIP->add_iterative(interior_model, solution, lin, "allvar", oacuts, 1e-8);
            //                    nb_count++;
            //                    obj_old=obj_new;
            //                }
            //                DebugOn("nb count "<<nb_count);
            //            }
            //           SOC_MIP->print_solution();
            //            if(norm1){
            //                apply_rot_trans(rot_trans, point_cloud_data);
            //            }
            //            else{
            //                apply_rot_trans(rot_trans[0], rot_trans[1], rot_trans[2], rot_trans[3], rot_trans[4], rot_trans[5], point_cloud_data);
            //            }
        }
#ifdef USE_MATPLOT
        plot(initial_point_cloud_model,initial_point_cloud_data,1);
        plot(point_cloud_model,point_cloud_data,1);
#endif
        DebugOn("L2 before optimization = " << to_string_with_precision(L2error_init,12) << endl);
        auto L2error_final = computeL2error(point_cloud_model,point_cloud_data,L2matching,L2err_per_point);
        auto L1error_final = computeL1error(point_cloud_model,point_cloud_data,L1matching,L1err_per_point);
        DebugOn("L2 after optimization on downsampled = " << to_string_with_precision(L2error_final,12) << endl);
        DebugOn("L1 after optimization on downsampled = " << to_string_with_precision(L1error_final,12) << endl);
        L2error_final = computeL2error(initial_point_cloud_model,initial_point_cloud_data,L2matching,L2err_per_point);
        L1error_final = computeL1error(initial_point_cloud_model,initial_point_cloud_data,L1matching,L1err_per_point);
        DebugOn("L2 after optimization on full set = " << to_string_with_precision(L2error_final,12) << endl);
        DebugOn("L1 after optimization on full set = " << to_string_with_precision(L1error_final,12) << endl);
        return 0;
    }
    
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

/* Run the MINLP ARMO model for registration */
tuple<double,double,double,double,double,double> run_ARMO_MINLP(bool bypass, string axis, const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data){
    double shift_min_x = 0.125, shift_max_x = 0.25, shift_min_y = -0.25,shift_max_y = -0.125,shift_min_z = -0.125,shift_max_z = 0;
    double yaw_min = -12.5*pi/180., yaw_max = 0, pitch_min = 12.5*pi/180.,pitch_max = 25.*pi/180.,roll_min = -12.5*pi/180.,roll_max = 0;
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
    Debug("nd = " << nd << endl);
    Debug("nm = " << nm << endl);
    
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
    var<> yaw("yaw", yaw_min, yaw_max), pitch("pitch", pitch_min, pitch_max), roll("roll", roll_min, roll_max);
    var<> x_shift("x_shift", shift_min_x, shift_max_x), y_shift("y_shift", shift_min_y, shift_max_y), z_shift("z_shift", shift_min_z, shift_max_z);
    //                    var<> yaw("yaw", -1e-6, 1e-6), pitch("pitch",-1e-6, 1e-6), roll("roll", -1e-6, 1e-6);
    //                    var<> x_shift("x_shift", 0, 0), y_shift("y_shift", 0, 0), z_shift("z_shift", 0, 0);
    var<> delta("delta", pos_);
    var<> delta_min("delta_min", pos_);
    var<> bin("bin",0,1);
    
    Reg.add(bin.in(cells));
    Reg.add(delta.in(cells), delta_min.in(N1));
    Reg.add(yaw.in(R(1)),pitch.in(R(1)),roll.in(R(1)));
    Reg.add(x_shift.in(R(1)),y_shift.in(R(1)),z_shift.in(R(1)));
    Reg.add(new_x1.in(N1), new_y1.in(N1), new_z1.in(N1));
    //        Reg.add(x_diff.in(cells), y_diff.in(cells), z_diff.in(cells));
    //                Reg.add(z_diff.in(cells));
    DebugOn("There are " << cells.size() << " cells" << endl);
    for (int i = 0; i<N1.size(); i++) {
        bin(to_string(i+1)+","+to_string(i+1)).set_lb(1);//=1;
    }
    
    
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
    
    solver<> S(Reg,ipopt);
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
vector<double> run_ARMO_Global(bool convex, string axis, const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data, bool norm1){};
//{
//    double yaw_min = -5*pi/180., yaw_max = 5, pitch_min = -5.*pi/180.,pitch_max = 5.*pi/180.,roll_min = -5*pi/180.,roll_max = 5*pi/180.;
//    double roll_1 = 0, yaw_1 = 0, pitch_1 = 0;
//    int nb_pairs = 0, min_nb_pairs = numeric_limits<int>::max(), max_nb_pairs = 0, av_nb_pairs = 0;
//    size_t nm = point_cloud_model.size(), nd = point_cloud_data.size();
//    vector<pair<double,double>> min_max_data;
//    vector<vector<pair<double,double>>> min_max_model(nm);
//    vector<int> nb_neighbors(nd);
//    vector<vector<int>> neighbors(nd);
//    vector<double> zeros = {0,0,0};
//
//    var<> theta11("theta11",  0, 1), theta12("theta12", -2, 2), theta13("theta13", -2, 2);
//    var<> theta21("theta21",  -1, 1), theta22("theta22", -2, 2), theta23("theta23", -2, 2);
//    var<> theta31("theta31",  -1, 1), theta32("theta32", -1, 1), theta33("theta33", 0, 1);
//
//    param<> x1("x1"), x2("x2"), y1("y1"), y2("y2"), z1("z1"), z2("z2");
//    param<> x_uav1("x_uav1"), y_uav1("y_uav1"), z_uav1("z_uav1");
//    param<> x_uav2("x_uav2"), y_uav2("y_uav2"), z_uav2("z_uav2");
//        //        return 0;
//    int m = av_nb_pairs;
//        //            int m = 1;
//    vector<double> min_dist(nd,numeric_limits<double>::max());
//    vector<int> nearest(nd);
//    vector<string> nearest_id(nd);
//    string i_str, j_str;
//    indices Pairs("Pairs"), cells("cells");
//    map<int,int> n2_map;
//    int idx1 = 0;
//    int idx2 = 0;
//    int nb_max_neigh = 1;
//    double dist_sq = 0;
//    /* Compute nearest points in data point cloud */
//    for (auto i = 0; i<nd; i++) {
//        i_str = to_string(i+1);
//        x1.add_val(i_str,point_cloud_data.at(i).at(0));
//        y1.add_val(i_str,point_cloud_data.at(i).at(1));
//        z1.add_val(i_str,point_cloud_data.at(i).at(2));
//    }
//    for (auto j = 0; j<nm; j++) {
//        j_str = to_string(j+1);
//        x2.add_val(j_str,point_cloud_model.at(j).at(0));
//        y2.add_val(j_str,point_cloud_model.at(j).at(1));
//        z2.add_val(j_str,point_cloud_model.at(j).at(2));
//    }
//    for (auto i = 0; i< nd; i++) {
//        double min_dist = numeric_limits<double>::max();
//        for (auto j = 0; j< nm; j++) {
//            dist_sq = std::pow(point_cloud_data.at(i).at(0) - point_cloud_model.at(j).at(0),2) + std::pow(point_cloud_data.at(i).at(1) - point_cloud_model.at(j).at(1),2) + std::pow(point_cloud_data.at(i).at(2) - point_cloud_model.at(j).at(2),2);
//            if(min_dist>dist_sq){
//                min_dist = dist_sq;
//                j_str = to_string(j+1);
//                nearest_id[i] = j_str;
//            }
//        }
//    }
//    idx1 = 0;
//    indices N1("N1"),N2("N2");
//    DebugOn("nd = " << nd << endl);
//    DebugOn("nm = " << nm << endl);
//
//    N1 = range(1,nd);
//    N2 = range(1,nm);
//    cells = indices(N1,N2);
//        //    cells.print();
//        //    for (auto i = 0; i<nd; i++) {
//        //        i_str = to_string(i+1);
//        //        for (auto j = 0; j<nm; j++) {
//        //            j_str = to_string(j+1);
//        //            cells.add(i_str+","+j_str);
//        //        }
//        //    }
//    Model<> Reg("Reg");
//    var<> new_x1("new_x1", -1,1), new_y1("new_y1", -1, 1), new_z1("new_z1", -1,1);
//    var<> x_diff("x_diff", pos_), y_diff("y_diff", pos_), z_diff("z_diff", pos_);
//
//        //            var<> yaw("yaw", thetaz, thetaz), pitch("pitch", thetax, thetax), roll("roll", thetay, thetay);
//        //            var<> x_shift("x_shift", 0.2163900, 0.2163900), y_shift("y_shift", -0.1497952, -0.1497952), z_shift("z_shift", 0.0745708, 0.0745708);
//    var<> cosr("cosr",  std::min(cos(roll_min),cos(roll_max)), 1), sinr("sinr", std::sin(roll_min), std::sin(roll_max));
//    var<> cosp("cosp",  std::min(cos(pitch_min),cos(pitch_max)), 1), sinp("sinp", std::sin(pitch_min), std::sin(pitch_max));
//    var<> cosy("cosy",  std::cos(std::min(yaw_min,yaw_max)), 1), siny("siny", std::sin(yaw_min), std::sin(yaw_max));
//    var<> cosy_sinr("cosy_sinr", std::sin(roll_min), std::sin(roll_max)), siny_sinr("siny_sinr", std::sin(yaw_min)*std::sin(roll_min), std::sin(yaw_max)*std::sin(roll_max));
//        //    x_rot1 -= (x1.in(N1))*cosy.in(ids1)*cosr.in(ids1) + (y1.in(N1))*(cosy_sinr.in(ids1)*sinp.in(ids1) - siny.in(ids1)*cosp.in(ids1)) + (z1.in(N1))*(cosy_sinr.in(ids1)*cosp.in(ids1) + siny.in(ids1)*sinp.in(ids1));
//        //    y_rot1 -= (x1.in(N1))*siny.in(ids1)*cosr.in(ids1) + (y1.in(N1))*(siny_sinr.in(ids1)*sinp.in(ids1) + cosy.in(ids1)*cosp.in(ids1)) + (z1.in(N1))*(siny_sinr.in(ids1)*cosp.in(ids1) - cosy.in(ids1)*sinp.in(ids1));
//        //    z_rot1 -= (x1.in(N1))*-1*sinr.in(ids1) + (y1.in(N1))*(cosr.in(ids1)*sinp.in(ids1)) + (z1.in(N1))*(cosr.in(ids1)*cosp.in(ids1));
//    var<> siny_sinp("siny_sinp", std::sin(yaw_min)*std::sin(pitch_min), std::sin(yaw_max)*std::sin(pitch_max));
//    var<> cosy_sinp("cosy_sinp", std::sin(pitch_min), std::sin(pitch_max));
//    var<> cosy_cosr("cosy_cosr", std::cos(std::max(abs(roll_min),abs(roll_max))), 1), cosy_sinr_sinp("cosy_sinr_sinp", -std::sin(angle_max)*std::sin(angle_max), std::sin(angle_max)*std::sin(angle_max));
//    var<> cosy_cosp("cosy_cosp", std::cos(angle_max), 1);
//    var<> siny_cosp("siny_cosp", -std::sin(angle_max), std::sin(angle_max)), cosy_sinr_cosp("cosy_sinr_cosp", -std::sin(angle_max), std::sin(angle_max));
//    var<> siny_cosr("siny_cosr", -std::sin(angle_max), std::sin(angle_max)), siny_sinr_sinp("siny_sinr_sinp", -std::pow(std::sin(angle_max),3), std::pow(std::sin(angle_max),3)), siny_sinr_cosp("siny_sinr_cosp", -std::sin(angle_max)*std::sin(angle_max), std::sin(angle_max)*std::sin(angle_max));
//    var<> cosr_sinp("cosr_sinp", -std::sin(angle_max), std::sin(angle_max)), cosr_cosp("cosr_cosp", std::cos(angle_max), 1);
//    var<> yaw("yaw", yaw_min, yaw_max), pitch("pitch", pitch_min, pitch_max), roll("roll", roll_min, roll_max);
//    var<> x_shift("x_shift", shift_min_x, shift_max_x), y_shift("y_shift", shift_min_y, shift_max_y), z_shift("z_shift", shift_min_z, shift_max_z);
//        //                    var<> yaw("yaw", -1e-6, 1e-6), pitch("pitch",-1e-6, 1e-6), roll("roll", -1e-6, 1e-6);
//        //                    var<> x_shift("x_shift", 0, 0), y_shift("y_shift", 0, 0), z_shift("z_shift", 0, 0);
//    var<> delta("delta", 0, 12);
//    var<> delta_min("delta_min", pos_);
//    var<int> bin("bin",0,1);
//    Reg.add(bin.in(cells));
//    Reg.add(delta.in(cells), delta_min.in(N1));
//    if(norm1){
//        Reg.add(theta11.in(R(1)),theta12.in(R(1)),theta13.in(R(1)));
//        Reg.add(theta21.in(R(1)),theta22.in(R(1)),theta23.in(R(1)));
//        Reg.add(theta31.in(R(1)),theta32.in(R(1)),theta33.in(R(1)));
//    }
//    else{
//        Reg.add(yaw.in(R(1)),pitch.in(R(1)),roll.in(R(1)));
//        Reg.add(cosr.in(R(1)),cosp.in(R(1)),cosy.in(R(1)));
//        Reg.add(sinr.in(R(1)),sinp.in(R(1)),siny.in(R(1)));
//        Reg.add(cosy_sinr.in(R(1)),siny_sinr.in(R(1)));
//    }
//    if(convex && !norm1){
//        Reg.add(siny_sinp.in(R(1)),cosy_sinp.in(R(1)));
//        Reg.add(cosy_cosr.in(R(1)), cosy_cosp.in(R(1)), cosy_sinr_sinp.in(R(1)));
//        Reg.add(siny_cosp.in(R(1)), cosy_sinr_cosp.in(R(1)));
//        Reg.add(siny_cosr.in(R(1)), siny_sinr_sinp.in(R(1)), siny_sinr_cosp.in(R(1)));
//        Reg.add(cosr_sinp.in(R(1)), cosr_cosp.in(R(1)));
//    }
//    Reg.add(x_shift.in(R(1)),y_shift.in(R(1)),z_shift.in(R(1)));
//    Reg.add(new_x1.in(N1), new_y1.in(N1), new_z1.in(N1));
//    Reg.add(x_diff.in(cells), y_diff.in(cells), z_diff.in(cells));
//        //        Reg.add(x_diff.in(cells), y_diff.in(cells), z_diff.in(cells));
//        //                Reg.add(z_diff.in(cells));
//    DebugOn("There are " << cells.size() << " cells" << endl);
//        //        for (int i = 0; i<nd; i++) {
//        //            bin(to_string(i+1)+","+to_string(i+1)).set_lb(1);
//        //        }
//        //    for (int i = 0; i<nd; i++) {
//        //        vector<var<>> delta_vec(nm);
//        //        for (int j = 0; j<nm; j++) {
//        //            delta_vec[j] = delta(to_string(i+1)+","+to_string(j+1));
//        //        }
//        //        Constraint<> DeltaMin("DeltaMin_"+to_string(i));
//        //        DeltaMin += delta_min[i+1] - min(delta_vec);
//        //        Reg.add(DeltaMin==0);
//        //    }
//
//    Constraint<> DeltaMin("DeltaMin");
//    DeltaMin = delta_min;
//    DeltaMin -= bin.in_matrix(1, 1)*delta.in_matrix(1, 1);
//    Reg.add(DeltaMin.in(N1)==0);
//
//
//    Constraint<> OneBin("OneBin");
//    OneBin = bin.in_matrix(1, 1);
//    Reg.add(OneBin.in(N1)==1);
//
//
//    Constraint<> OneBin2("OneBin2");
//    OneBin2 = bin.in_matrix(0, 1);
//    Reg.add(OneBin2.in(N2)<=1);
//
//    if(norm1){
//        Constraint<> x_abs1("x_abs1");
//        x_abs1 += x_diff - (new_x1.from(cells) - x2.to(cells));
//        Reg.add(x_abs1.in(cells)>=0);
//
//        Constraint<> x_abs2("x_abs2");
//        x_abs2 += x_diff + (new_x1.from(cells) - x2.to(cells));
//        Reg.add(x_abs2.in(cells)>=0);
//
//        Constraint<> y_abs1("y_abs1");
//        y_abs1 += y_diff - (new_y1.from(cells) - y2.to(cells));
//        Reg.add(y_abs1.in(cells)>=0);
//
//        Constraint<> y_abs2("y_abs2");
//        y_abs2 += y_diff + (new_y1.from(cells) - y2.to(cells));
//        Reg.add(y_abs2.in(cells)>=0);
//
//        Constraint<> z_abs1("z_abs1");
//        z_abs1 += z_diff - (new_z1.from(cells) - z2.to(cells));
//        Reg.add(z_abs1.in(cells)>=0);
//
//        Constraint<> z_abs2("z_abs2");
//        z_abs2 += z_diff + (new_z1.from(cells) - z2.to(cells));
//        Reg.add(z_abs2.in(cells)>=0);
//
//        Constraint<> Norm1("Norm1");
//        Norm1 += delta - (x_diff + y_diff + z_diff);
//        Reg.add(Norm1.in(cells)>=0);
//
//        auto ids1 = theta11.repeat_id(cells.size());
//        Constraint<> x_rot1("x_rot1");
//        x_rot1 += new_x1 -x_shift;
//        x_rot1 -= x1.in(N1)*theta11.in(ids1) + y1.in(N1)*theta12.in(ids1) + z1.in(N1)*theta13.in(ids1);
//        Reg.add(x_rot1.in(N1)==0);
//
//        Constraint<> y_rot1("y_rot1");
//        y_rot1 += new_y1 - y_shift;
//        y_rot1 -= x1.in(N1)*theta21.in(ids1) + y1.in(N1)*theta22.in(ids1) + z1.in(N1)*theta23.in(ids1);
//        Reg.add(y_rot1.in(N1)==0);
//
//        Constraint<> z_rot1("z_rot1");
//        z_rot1 += new_z1 -z_shift;
//        z_rot1 -= x1.in(N1)*theta31.in(ids1) + y1.in(N1)*theta32.in(ids1) + z1.in(N1)*theta33.in(ids1);
//        Reg.add(z_rot1.in(N1)==0);
//    }
//    else{
//        Constraint<> Norm2("Norm2");
//        Norm2 += delta - pow(new_x1.from(cells) - x2.to(cells),2) - pow(new_y1.from(cells) - y2.to(cells),2) - pow(new_z1.from(cells) - z2.to(cells),2);
//        Reg.add(Norm2.in(cells)>=0);
//
//
//        Constraint<> trigR("trigR");
//        trigR = pow(cosr,2) + pow(sinr,2);
//        if(!convex)
//            Reg.add(trigR==1);
//        else
//            Reg.add(trigR<=1);
//
//        Constraint<> trigP("trigP");
//        trigP = pow(cosp,2) + pow(sinp,2);
//        if(!convex)
//            Reg.add(trigP==1);
//        else
//            Reg.add(trigP<=1);
//
//        Constraint<> trigY("trigY");
//        trigY = pow(cosy,2) + pow(siny,2);
//        if(!convex)
//            Reg.add(trigY==1);
//        else
//            Reg.add(trigY<=1);
//
//            //    Constraint<> cos_roll("cos_roll");
//            //    cos_roll = cosr - cos(roll);
//            //    Reg.add(cos_roll==0);
//            //
//            //    Constraint<> sin_roll("sin_roll");
//            //    sin_roll = sinr - sin(roll);
//            //    Reg.add(sin_roll==0);
//            //
//            //    Constraint<> cos_pitch("cos_pitch");
//            //    cos_pitch = cosp - cos(pitch);
//            //    Reg.add(cos_pitch==0);
//            //
//            //    Constraint<> sin_pitch("sin_pitch");
//            //    sin_pitch = sinp - sin(pitch);
//            //    Reg.add(sin_pitch==0);
//            //
//            //    Constraint<> cos_yaw("cos_yaw");
//            //    cos_yaw = cosy - cos(yaw);
//            //    Reg.add(cos_yaw==0);
//            //
//            //    Constraint<> sin_yaw("sin_yaw");
//            //    sin_yaw = siny - sin(yaw);
//            //    Reg.add(sin_yaw==0);
//
//
//        if(!convex){
//            Constraint<> cosy_sinr_prod("cosy_sinr");
//            cosy_sinr_prod = cosy_sinr - cosy*sinr;
//            Reg.add(cosy_sinr_prod==0);
//        }
//        else{
//            Reg.add_McCormick("cosy_sinr", cosy_sinr, cosy, sinr);
//        }
//
//        if(!convex){
//            Constraint<> siny_sinr_prod("siny_sinr");
//            siny_sinr_prod = siny_sinr - siny*sinr;
//            Reg.add(siny_sinr_prod==0);
//        }
//        else {
//            Reg.add_McCormick("siny_sinr", siny_sinr, siny, sinr);
//        }
//
//        auto ids1 = cosy.repeat_id(cells.size());
//
//        if(!convex){
//            /* alpha = yaw_, beta = roll and gamma = pitch */
//            Constraint<> x_rot1("x_rot1");
//            x_rot1 += new_x1 - x_shift.in(ids1);
//            x_rot1 -= (x1.in(N1))*cosy.in(ids1)*cosr.in(ids1) + (y1.in(N1))*(cosy_sinr.in(ids1)*sinp.in(ids1) - siny.in(ids1)*cosp.in(ids1)) + (z1.in(N1))*(cosy_sinr.in(ids1)*cosp.in(ids1) + siny.in(ids1)*sinp.in(ids1));
//            Reg.add(x_rot1.in(N1)==0);
//
//
//            Constraint<> y_rot1("y_rot1");
//            y_rot1 += new_y1 - y_shift.in(ids1);
//            y_rot1 -= (x1.in(N1))*siny.in(ids1)*cosr.in(ids1) + (y1.in(N1))*(siny_sinr.in(ids1)*sinp.in(ids1) + cosy.in(ids1)*cosp.in(ids1)) + (z1.in(N1))*(siny_sinr.in(ids1)*cosp.in(ids1) - cosy.in(ids1)*sinp.in(ids1));
//            Reg.add(y_rot1.in(N1)==0);
//
//            Constraint<> z_rot1("z_rot1");
//            z_rot1 += new_z1 - z_shift.in(ids1);
//            z_rot1 -= (x1.in(N1))*-1*sinr.in(ids1) + (y1.in(N1))*(cosr.in(ids1)*sinp.in(ids1)) + (z1.in(N1))*(cosr.in(ids1)*cosp.in(ids1));
//            Reg.add(z_rot1.in(N1)==0);
//        }
//        else {
//                //        Reg.add_McCormick("cosy_sinr", cosy_sinr, cosy, sinr);
//                //        Reg.add_McCormick("siny_sinr", siny_sinr, siny, sinr);
//            Reg.add_McCormick("cosy_cosr", cosy_cosr, cosy, cosr);
//            Reg.add_McCormick("cosy_sinr_sinp", cosy_sinr_sinp, cosy_sinr, sinp);
//            Reg.add_McCormick("siny_cosp", siny_cosp, siny, cosp);
//            Reg.add_McCormick("siny_sinp", siny_sinp, siny, sinp);
//            Reg.add_McCormick("siny_cosr", siny_cosr, siny, cosr);
//            Reg.add_McCormick("siny_sinr_sinp", siny_sinr_sinp, siny_sinr, sinp);
//            Reg.add_McCormick("cosy_cosp", cosy_cosp, cosy, cosp);
//            Reg.add_McCormick("siny_sinr_cosp", siny_sinr_cosp, siny_sinr, cosp);
//            Reg.add_McCormick("cosy_sinp", cosy_sinp, cosy, sinp);
//            Reg.add_McCormick("cosr_sinp", cosr_sinp, cosr, sinp);
//            Reg.add_McCormick("cosr_cosp", cosr_cosp, cosr, cosp);
//
//                //        point_cloud[i][0] = shifted_x*cos(alpha)*cos(beta) + shifted_y*(cos(alpha)*sin(beta)*sin(gamma) - sin(alpha)*cos(gamma)) + shifted_z*(cos(alpha)*sin(beta)*cos(gamma) + sin(alpha)*sin(gamma));
//                //        point_cloud[i][1] = shifted_x*sin(alpha)*cos(beta) + shifted_y*(sin(alpha)*sin(beta)*sin(gamma) + cos(alpha)*cos(gamma)) + shifted_z*(sin(alpha)*sin(beta)*cos(gamma) - cos(alpha)*sin(gamma));
//                //        point_cloud[i][2] = shifted_x*(-sin(beta)) + shifted_y*(cos(beta)*sin(gamma)) + shifted_z*(cos(beta)*cos(gamma));
//                //        double beta = roll*pi/180;// roll in radians
//                //        double gamma = pitch*pi/180; // pitch in radians
//                //        double alpha = yaw*pi/180; // yaw in radians
//            Constraint<> x_rot1("x_rot1");
//            x_rot1 += new_x1 - x_shift.in(ids1);
//            x_rot1 -= (x1.in(N1))*cosy_cosr.in(ids1) + (y1.in(N1))*(cosy_sinr_sinp.in(ids1) - siny_cosp.in(ids1)) + (z1.in(N1))*(cosy_sinr_cosp.in(ids1) + siny_sinp.in(ids1));
//            Reg.add(x_rot1.in(N1)==0);
//
//
//            Constraint<> y_rot1("y_rot1");
//            y_rot1 += new_y1 - y_shift.in(ids1);
//            y_rot1 -= (x1.in(N1))*siny_cosr.in(ids1) + (y1.in(N1))*(siny_sinr_sinp.in(ids1) + cosy_cosp.in(ids1)) + (z1.in(N1))*(siny_sinr_cosp.in(ids1) - cosy_sinp.in(ids1));
//            Reg.add(y_rot1.in(N1)==0);
//
//            Constraint<> z_rot1("z_rot1");
//            z_rot1 += new_z1 - z_shift.in(ids1);
//            z_rot1 -= (x1.in(N1))*-1*sinr.in(ids1) + (y1.in(N1))*(cosr_sinp.in(ids1)) + (z1.in(N1))*(cosr_cosp.in(ids1));
//            Reg.add(z_rot1.in(N1)==0);
//        }
//    }
//
//        //    M.min(sum(z_diff)/nb_overlap);
//
//        //        M.min(sum(z_diff));
//    if(axis == "full")
//            //            Reg.min(sum(x_diff) + sum(y_diff) + sum(z_diff));
//            //                Reg.min(sum(x_diff)/cells.size() + sum(y_diff)/cells.size() + sum(z_diff)/cells.size());
//        Reg.min(sum(delta_min));
//    else if(axis == "x")
//        Reg.min(sum(x_diff)/cells.size());
//    else if (axis == "y")
//        Reg.min(sum(y_diff)/cells.size());
//    else
//        Reg.min(sum(z_diff)/cells.size());
//
//        //                Reg.min(sum(x_diff)/cells.size() + sum(y_diff)/cells.size() + sum(z_diff)/cells.size());
//
//    Reg.print();
//
//    if(convex){
//            //        Reg.replace_integers();
//            //        auto Rel = Reg.relax();
//        solver<> S(Reg,gurobi);
//        S.run();
//    }
//    else {
//        solver<> S(Reg,gurobi);
//        S.run();
//    }
//    Reg.print_solution();
//        //        S.run(0, 1e-10, 1000);
//
//
//        //        for (int i = 0; i<500; i++) {
//        //            pre_x.add_val(x_rot1.eval(i));
//        //            pre_y.add_val(y_rot1.eval(i));
//        //            pre_z.add_val(z_rot1.eval(i));
//        //            x_uav.add_val(x_uav1.eval(i));
//        //            y_uav.add_val(y_uav1.eval(i));
//        //            z_uav.add_val(z_uav1.eval(i));
//        //        }
//        //        for (int i = 0; i<500; i++) {
//        //            pre_x.add_val(x_rot2.eval(i));
//        //            pre_y.add_val(y_rot2.eval(i));
//        //            pre_z.add_val(z_rot2.eval(i));
//        //            x_uav.add_val(x_uav2.eval(i));
//        //            y_uav.add_val(y_uav2.eval(i));
//        //            z_uav.add_val(z_uav2.eval(i));
//        //        }
//        //    M.print_solution();
//        //    roll.in(R(1));pitch.in(R(1));yaw.in(R(1));
//    vector<double> rot_trans;
//    if(norm1){
//        rot_trans.resize(12);
//        DebugOn("Theta matrix = " << endl);
//        DebugOn("|" << theta11.eval() << " " << theta12.eval() << " " << theta13.eval() << "|" << endl);
//        DebugOn("|" << theta21.eval() << " " << theta22.eval() << " " << theta23.eval() << "|" << endl);
//        DebugOn("|" << theta31.eval() << " " << theta32.eval() << " " << theta33.eval() << "|" << endl);
//        rot_trans[0]=theta11.eval();
//        rot_trans[1]=theta12.eval();
//        rot_trans[2]=theta13.eval();;
//        rot_trans[3]=theta21.eval();
//        rot_trans[4]=theta22.eval();
//        rot_trans[5]=theta23.eval();
//        rot_trans[6]=theta31.eval();
//        rot_trans[7]=theta32.eval();
//        rot_trans[8]=theta33.eval();
//        rot_trans[9]=x_shift.eval();
//        rot_trans[10]=y_shift.eval();
//        rot_trans[11]=z_shift.eval();
//        DebugOn("x shift = " << x_shift.eval() << endl);
//        DebugOn("y shift = " << y_shift.eval() << endl);
//        DebugOn("z shift = " << z_shift.eval() << endl);
//    }
//    else{
//        rot_trans.resize(6);
//        pitch = std::asin(sinp.eval());
//        roll = std::asin(sinr.eval());
//        yaw = std::asin(siny.eval());
//        DebugOn("Roll (degrees) = " << to_string_with_precision(roll.eval()*180/pi,12) << endl);
//        DebugOn("Pitch (degrees) = " << to_string_with_precision(pitch.eval()*180/pi,12) << endl);
//        DebugOn("Yaw (degrees) = " << to_string_with_precision(yaw.eval()*180/pi,12) << endl);
//        DebugOn("x shift = " << x_shift.eval() << endl);
//        DebugOn("y shift = " << y_shift.eval() << endl);
//        DebugOn("z shift = " << z_shift.eval() << endl);
//        roll_1 = roll.eval()*180/pi;
//        pitch_1 = pitch.eval()*180/pi;
//        yaw_1 = yaw.eval()*180/pi;
//        rot_trans[0] = roll_1;
//        rot_trans[1] = pitch_1;
//        rot_trans[2] = yaw_1;
//        rot_trans[3] = x_shift.eval();
//        rot_trans[4] = y_shift.eval();
//        rot_trans[5] = z_shift.eval();
//    }
//    return rot_trans;
//}

/* Run the Reformulated Global ARMO model for registration, given in ARMO_inequalities */
tuple<double,double,double,double,double,double> run_ARMO_Global_reform(bool convex, string axis, const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data){
    double angle_max = 3.14, shift_max = 0.5;
    double roll_1 = 0, yaw_1 = 0, pitch_1 = 0;
    int nb_pairs = 0, min_nb_pairs = numeric_limits<int>::max(), max_nb_pairs = 0, av_nb_pairs = 0;
    size_t nm = point_cloud_model.size(), nd = point_cloud_data.size();
    
    vector<double> zeros = {0,0,0};
    
    param<> x1("x1"), x2("x2"), y1("y1"), y2("y2"), z1("z1"), z2("z2"), nx2("nx2"), ny2("ny2"), nz2("nz2");
    param<> x_uav1("x_uav1"), y_uav1("y_uav1"), z_uav1("z_uav1");
    param<> x_uav2("x_uav2"), y_uav2("y_uav2"), z_uav2("z_uav2");
    //        return 0;
    int m = av_nb_pairs;
    
    string i_str, j_str;
    indices Pairs("Pairs"), cells("cells");
    //map<int,int> n2_map;
    int idx1 = 0;
    int idx2 = 0;
    
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
    
    double x,y,z,nx,ny,nz,kx,ky,kz, d=10000, dist_min=100000;
    for (auto j = 0; j<nm; j++) {
        j_str = to_string(j+1);
        x=x2.eval(j_str);
        y=y2.eval(j_str);
        z=z2.eval(j_str);
        dist_min=15;
        for (auto k = 0; k<nm; k++) {
            if(k!=j){
                auto k_str = to_string(k+1);
                kx=x2.eval(k_str);
                ky=y2.eval(k_str);
                kz=z2.eval(k_str);
                d=std::pow(kx-x,2)+std::pow(ky-y,2)+std::pow(kz-z,2);
                if(d<dist_min){
                    nx=kx;
                    ny=ky;
                    nz=kz;
                    dist_min=d;
                }
            }
            
        }
        nx2.add_val(j_str,nx);
        ny2.add_val(j_str,ny);
        nz2.add_val(j_str,nz);
    }
    
    
    idx1 = 0;
    indices N1("N1"),N2("N2");
    Debug("nd = " << nd << endl);
    Debug("nm = " << nm << endl);
    
    N1 = range(1,nd);
    N2 = range(1,nm);
    cells = indices(N1,N2);
    
    Model<> Reg("Reg");
    var<> new_x1("new_x1", -1, 1), new_y1("new_y1", -1, 1), new_z1("new_z1", -1, 1);
    var<> new_nx("new_nx", -1, 1), new_ny("new_ny", -1, 1), new_nz("new_nz", -1, 1);
    var<> new_xm("new_xm", -1, 1), new_ym("new_ym", -1, 1), new_zm("new_zm", -1, 1);
    var<> x_diff("x_diff", pos_), y_diff("y_diff", pos_), z_diff("z_diff", pos_);
    
    
    /*  bool bounded=true;
     if(!bounded){
     var<> cosr("cosr",  -1, 1), sinr("sinr", -1, 1);
     var<> cosp("cosp",   0, 1), sinp("sinp", -1, 1);
     var<> cosy("cosy",  -1, 1), siny("siny", -1, 1);
     var<> cosy_sinr("cosy_sinr", -1, 1), siny_sinr("siny_sinr", -1, 1);
     var<> siny_sinp("siny_sinp", -1, 1);
     var<> cosy_sinp("cosy_sinp", -1, 1);
     var<> cosy_cosr("cosy_cosr", -1, 1), cosy_sinr_sinp("cosy_sinr_sinp", -1, 1);
     var<> cosy_cosp("cosy_cosp", -1, 1);
     var<> siny_cosp("siny_cosp", -1, 1), cosy_sinr_cosp("cosy_sinr_cosp", -1, 1);
     var<> siny_cosr("siny_cosr", -1, 1), siny_sinr_sinp("siny_sinr_sinp", -1, 1), siny_sinr_cosp("siny_sinr_cosp", -1,1);
     var<> cosr_sinp("cosr_sinp", -1,1), cosr_cosp("cosr_cosp", -1, 1);
     }*/
    
    angle_max=1;
    var<> cosr("cosr",  std::cos(angle_max), 1), sinr("sinr", -std::sin(angle_max), std::sin(angle_max));
    var<> cosp("cosp",  std::cos(angle_max), 1), sinp("sinp", -std::sin(angle_max), std::sin(angle_max));
    var<> cosy("cosy",  std::cos(angle_max), 1), siny("siny", -std::sin(angle_max), std::sin(angle_max));
    var<> cosy_sinr("cosy_sinr", -std::sin(angle_max), std::sin(angle_max)), siny_sinr("siny_sinr", -std::sin(angle_max)*std::sin(angle_max), std::sin(angle_max)*std::sin(angle_max));
    
    var<> siny_sinp("siny_sinp", -std::sin(angle_max)*std::sin(angle_max), std::sin(angle_max)*std::sin(angle_max));
    var<> cosy_sinp("cosy_sinp", -std::sin(angle_max), std::sin(angle_max));
    var<> cosy_cosr("cosy_cosr", std::cos(angle_max), 1), cosy_sinr_sinp("cosy_sinr_sinp", -std::sin(angle_max)*std::sin(angle_max), std::sin(angle_max)*std::sin(angle_max));
    var<> cosy_cosp("cosy_cosp", std::cos(angle_max), 1);
    var<> siny_cosp("siny_cosp", -std::sin(angle_max), std::sin(angle_max)), cosy_sinr_cosp("cosy_sinr_cosp", -std::sin(angle_max), std::sin(angle_max));
    var<> siny_cosr("siny_cosr", -std::sin(angle_max), std::sin(angle_max)), siny_sinr_sinp("siny_sinr_sinp", -std::pow(std::sin(angle_max),3), std::pow(std::sin(angle_max),3)), siny_sinr_cosp("siny_sinr_cosp", -std::sin(angle_max)*std::sin(angle_max), std::sin(angle_max)*std::sin(angle_max));
    var<> cosr_sinp("cosr_sinp", -std::sin(angle_max), std::sin(angle_max)), cosr_cosp("cosr_cosp", std::cos(angle_max), 1);
    
    //var<> yaw("yaw", -angle_max, angle_max), pitch("pitch", -angle_max, angle_max), roll("roll", -angle_max, angle_max);
    var<> x_shift("x_shift", -shift_max, shift_max), y_shift("y_shift", -shift_max, shift_max), z_shift("z_shift", -shift_max, shift_max);
    
    var<> delta("delta", 0,12);
    
    var<> bin("bin",0,1);
    Reg.add(bin.in(cells));
    DebugOn("Added binary variables" << endl);
    Reg.add(delta.in(N1));
    Reg.add(cosr.in(R(1)),cosp.in(R(1)),cosy.in(R(1)));
    Reg.add(sinr.in(R(1)),sinp.in(R(1)),siny.in(R(1)));
    Reg.add(cosy_sinr.in(R(1)),siny_sinr.in(R(1)));
    if(convex){
        Reg.add(siny_sinp.in(R(1)),cosy_sinp.in(R(1)));
        Reg.add(cosy_cosr.in(R(1)), cosy_cosp.in(R(1)), cosy_sinr_sinp.in(R(1)));
        Reg.add(siny_cosp.in(R(1)), cosy_sinr_cosp.in(R(1)));
        Reg.add(siny_cosr.in(R(1)), siny_sinr_sinp.in(R(1)), siny_sinr_cosp.in(R(1)));
        Reg.add(cosr_sinp.in(R(1)), cosr_cosp.in(R(1)));
    }
    Reg.add(x_shift.in(R(1)),y_shift.in(R(1)),z_shift.in(R(1)));
    Reg.add(new_x1.in(N1), new_y1.in(N1), new_z1.in(N1));
    Reg.add(new_xm.in(N1), new_ym.in(N1), new_zm.in(N1));
    DebugOn("There are " << cells.size() << " cells" << endl);
    
    indices ids = indices("in_x");
    ids.add_empty_row();
    
    for(auto i=0;i<nd;i++){
        for(auto j=1;j<=nm;j++){
            if(cells.has_key(to_string(i+1)+","+to_string(j)))
                ids.add_in_row(i, to_string(j));
        }
    }
    
    
    Constraint<> Def_newxm("Def_newxm");
    Def_newxm = new_xm-product(x2.in(ids),bin.in_matrix(1, 1));
    Reg.add(Def_newxm.in(N1)==0);
    
    Constraint<> Def_newym("Def_newym");
    Def_newym = new_ym-product(y2.in(ids),bin.in_matrix(1, 1));
    Reg.add(Def_newym.in(N1)==0);
    
    Constraint<> Def_newzm("Def_newzm");
    Def_newzm = new_zm-product(z2.in(ids),bin.in_matrix(1, 1));
    Reg.add(Def_newzm.in(N1)==0);
    
    Constraint<> OneBin("OneBin");
    OneBin = bin.in_matrix(1, 1);
    Reg.add(OneBin.in(N1)>=1);
    Constraint<> OneBin2("OneBin2");
    OneBin2 = bin.in_matrix(0, 1);
    Reg.add(OneBin2.in(N2)<=1);
    
    //Can also try hull relaxation of the big-M here
    bool vi_M=false;
    if(vi_M){
        Constraint<> VI_M("VI_M");
        VI_M = 2*((x2.to(cells)-nx2.to(cells))*new_x1.from(cells)+(y2.to(cells)-ny2.to(cells))*new_y1.from(cells)+(z2.to(cells)-nz2.to(cells))*new_z1.from(cells))+ ((pow(nx2.to(cells),2)+pow(ny2.to(cells),2)+pow(nz2.to(cells),2))-(pow(x2.to(cells),2)+pow(y2.to(cells),2)+pow(z2.to(cells),2)))*bin.in(cells)+(3)*(1-bin.in(cells));
        Reg.add(VI_M.in(cells)>=0);
    }
    bool vi_reform=false;
    if(vi_reform){
        // VI.print_symbolic();
        bool vi_nonconvex=true;
        Reg.add(new_nx.in(N1), new_ny.in(N1), new_nz.in(N1));
        
        Constraint<> Def_newnx("Def_newnx");
        Def_newnx = new_nx-product(nx2.in(ids),bin.in_matrix(1, 1));
        Reg.add(Def_newnx.in(N1)==0);
        
        Constraint<> Def_newny("Def_newny");
        Def_newny = new_ny-product(ny2.in(ids),bin.in_matrix(1, 1));
        Reg.add(Def_newny.in(N1)==0);
        
        Constraint<> Def_newnz("Def_newnz");
        Def_newnz = new_nz-product(nz2.in(ids),bin.in_matrix(1, 1));
        Reg.add(Def_newnz.in(N1)==0);
        
        
        if(vi_nonconvex){
            
            Constraint<> VI_nonconvex("VI_nonconvex");
            VI_nonconvex = 2*((new_xm-new_nx)*new_x1+(new_ym-new_ny)*new_y1+(new_zm-new_nz)*new_z1)+ ((pow(new_nx,2)+pow(new_ny,2)+pow(new_nz,2))-(pow(new_xm,2)+pow(new_ym,2)+pow(new_zm,2)));
            Reg.add(VI_nonconvex.in(N1)>=0);
        }
        else{
            
            
            var<> px("px", -1, 1), py("py", -1, 1), pz("pz", -1, 1), nlift("nlift", 0, 3);
            Reg.add(px.in(N1),py.in(N1),pz.in(N1));
            
            Constraint<> Def_px_U("Def_px_U");
            Def_px_U = px.from(cells)-new_x1.from(cells)*(x2.to(cells)-nx2.to(cells))+bin.in(cells)-1;
            Reg.add(Def_px_U.in(cells)<=0);
            
            Constraint<> Def_px_L("Def_px_L");
            Def_px_L = px.from(cells)-new_x1.from(cells)*(x2.to(cells)-nx2.to(cells))-bin.in(cells)+1;
            Reg.add(Def_px_L.in(cells)>=0);
            
            Constraint<> Def_py_U("Def_py_U");
            Def_py_U = py.from(cells)-new_y1.from(cells)*(y2.to(cells)-ny2.to(cells))+bin.in(cells)-1;
            Reg.add(Def_py_U.in(cells)<=0);
            
            Constraint<> Def_py_L("Def_py_L");
            Def_py_L = py.from(cells)-new_y1.from(cells)*(y2.to(cells)-ny2.to(cells))-bin.in(cells)+1;
            Reg.add(Def_py_L.in(cells)>=0);
            
            Constraint<> Def_pz_U("Def_pz_U");
            Def_pz_U = pz.from(cells)-new_z1.from(cells)*(z2.to(cells)-nz2.to(cells))+bin.in(cells)-1;
            Reg.add(Def_pz_U.in(cells)<=0);
            
            Constraint<> Def_pz_L("Def_pz_L");
            Def_pz_L = pz.from(cells)-new_z1.from(cells)*(z2.to(cells)-nz2.to(cells))-bin.in(cells)+1;
            Reg.add(Def_pz_L.in(cells)>=0);
            
            if(convex){
                Reg.add(nlift.in(N1));
                
                Constraint<> Def_nlift("Def_nlift");
                Def_nlift = nlift-pow(new_nx,2)-pow(new_ny,2)-pow(new_nz,2);
                Reg.add(Def_nlift.in(N1)>=0);
                
                
                Constraint<> VI_convex("VI_convex");
                VI_convex = 2*(px+py+pz)+nlift-(pow(new_xm,2)+pow(new_ym,2)+pow(new_zm,2));
                Reg.add(VI_convex.in(N1)>=0);
            }
            else{
                Constraint<> VI_nc("VI_nc");
                VI_nc = 2*(px+py+pz)+(pow(new_nx,2)+pow(new_ny,2)+pow(new_nz,2))-(pow(new_xm,2)+pow(new_ym,2)+pow(new_zm,2));
                Reg.add(VI_nc.in(N1)>=0);
            }
        }
        
        
        
    }
    
    
    Constraint<> Norm2("Norm2");
    Norm2 += delta - pow(new_x1 - new_xm,2) - pow(new_y1 - new_ym,2) - pow(new_z1 - new_zm,2);
    Reg.add(Norm2.in(N1)>=0);
    
    Constraint<> trigR("trigR");
    trigR = pow(cosr,2) + pow(sinr,2);
    if(!convex)
        Reg.add(trigR==1);
    else
        Reg.add(trigR<=1);
    
    Constraint<> trigP("trigP");
    trigP = pow(cosp,2) + pow(sinp,2);
    if(!convex)
        Reg.add(trigP==1);
    else
        Reg.add(trigP<=1);
    
    Constraint<> trigY("trigY");
    trigY = pow(cosy,2) + pow(siny,2);
    if(!convex)
        Reg.add(trigY==1);
    else
        Reg.add(trigY<=1);
    
    
    if(!convex){
        Constraint<> cosy_sinr_prod("cosy_sinr");
        cosy_sinr_prod = cosy_sinr - cosy*sinr;
        Reg.add(cosy_sinr_prod==0);
    }
    else{
        Reg.add_McCormick("cosy_sinr", cosy_sinr, cosy, sinr);
    }
    
    if(!convex){
        Constraint<> siny_sinr_prod("siny_sinr");
        siny_sinr_prod = siny_sinr - siny*sinr;
        Reg.add(siny_sinr_prod==0);
    }
    else {
        Reg.add_McCormick("siny_sinr", siny_sinr, siny, sinr);
    }
    
    auto ids1 = cosy.repeat_id(cells.size());
    
    if(!convex){
        /* alpha = yaw_, beta = roll and gamma = pitch */
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
    }
    else {
        //        Reg.add_McCormick("cosy_sinr", cosy_sinr, cosy, sinr);
        //        Reg.add_McCormick("siny_sinr", siny_sinr, siny, sinr);
        Reg.add_McCormick("cosy_cosr", cosy_cosr, cosy, cosr);
        Reg.add_McCormick("cosy_sinr_sinp", cosy_sinr_sinp, cosy_sinr, sinp);
        Reg.add_McCormick("siny_cosp", siny_cosp, siny, cosp);
        Reg.add_McCormick("siny_sinp", siny_sinp, siny, sinp);
        Reg.add_McCormick("siny_cosr", siny_cosr, siny, cosr);
        Reg.add_McCormick("siny_sinr_sinp", siny_sinr_sinp, siny_sinr, sinp);
        Reg.add_McCormick("cosy_cosp", cosy_cosp, cosy, cosp);
        Reg.add_McCormick("siny_sinr_cosp", siny_sinr_cosp, siny_sinr, cosp);
        Reg.add_McCormick("cosy_sinp", cosy_sinp, cosy, sinp);
        Reg.add_McCormick("cosr_sinp", cosr_sinp, cosr, sinp);
        Reg.add_McCormick("cosr_cosp", cosr_cosp, cosr, cosp);
        
        
        Constraint<> x_rot1("x_rot1");
        x_rot1 += new_x1 - x_shift.in(ids1);
        x_rot1 -= (x1.in(N1))*cosy_cosr.in(ids1) + (y1.in(N1))*(cosy_sinr_sinp.in(ids1) - siny_cosp.in(ids1)) + (z1.in(N1))*(cosy_sinr_cosp.in(ids1) + siny_sinp.in(ids1));
        Reg.add(x_rot1.in(N1)==0);
        
        
        Constraint<> y_rot1("y_rot1");
        y_rot1 += new_y1 - y_shift.in(ids1);
        y_rot1 -= (x1.in(N1))*siny_cosr.in(ids1) + (y1.in(N1))*(siny_sinr_sinp.in(ids1) + cosy_cosp.in(ids1)) + (z1.in(N1))*(siny_sinr_cosp.in(ids1) - cosy_sinp.in(ids1));
        Reg.add(y_rot1.in(N1)==0);
        
        Constraint<> z_rot1("z_rot1");
        z_rot1 += new_z1 - z_shift.in(ids1);
        z_rot1 -= (x1.in(N1))*-1*sinr.in(ids1) + (y1.in(N1))*(cosr_sinp.in(ids1)) + (z1.in(N1))*(cosr_cosp.in(ids1));
        Reg.add(z_rot1.in(N1)==0);
    }
    
    
    
    if(axis == "full")
        
        Reg.min(sum(delta));
    else if(axis == "x")
        Reg.min(sum(x_diff)/cells.size());
    else if (axis == "y")
        Reg.min(sum(y_diff)/cells.size());
    else
        Reg.min(sum(z_diff)/cells.size());
    
    Reg.print();
    
    if(convex){
        //        Reg.replace_integers();
        //        auto Rel = Reg.relax();
        solver<> S(Reg,ipopt);
        S.run();
    }
    else {
        solver<> S(Reg,ipopt);
        S.run();
    }
    Reg.print_solution();
    //        S.run(0, 1e-10, 1000);
    
    
    auto pitch = std::atan2(sinp.eval(), cosp.eval());
    auto roll = std::atan2(sinr.eval(),cosr.eval());
    auto yaw = std::atan2(siny.eval(),cosy.eval());
    DebugOn("Roll (degrees) = " << to_string_with_precision(roll*180/pi,12) << endl);
    DebugOn("Pitch (degrees) = " << to_string_with_precision(pitch*180/pi,12) << endl);
    DebugOn("Yaw (degrees) = " << to_string_with_precision(yaw*180/pi,12) << endl);
    DebugOn("x shift = " << x_shift.eval() << endl);
    DebugOn("y shift = " << y_shift.eval() << endl);
    DebugOn("z shift = " << z_shift.eval() << endl);
    roll_1 = roll*180/pi;
    pitch_1 = pitch*180/pi;
    yaw_1 = yaw*180/pi;
    return {roll_1, pitch_1, yaw_1, x_shift.eval(), y_shift.eval(), z_shift.eval()};
}

#ifdef USE_QHULL


vector<pair<double, double>> get_translation_bounds(vector<pair<double, double>>min_max_model, const vector<double>& flat_point_cloud_data){
    double lmin=0, limax=0;
    double txmin, txmax, tymin, tymax, tzmin, tzmax;
    vector<pair<double, double>> res;
    Qhull qt;
    lmin=largest_inscribed_sphere_centre(0, 0, 0, flat_point_cloud_data, limax, qt);
    txmin=min_max_model[0].first+lmin;
    txmax=min_max_model[0].second-lmin;
    tymin=min_max_model[1].first+lmin;
    tymax=min_max_model[1].second-lmin;
    tzmin=min_max_model[2].first+lmin;
    tzmax=min_max_model[2].second-lmin;
    res.push_back(make_pair(txmin, txmax));
    res.push_back(make_pair(tymin, tymax));
    res.push_back(make_pair(tzmin, tzmax));
    return res;
}
double largest_inscribed_sphere_centre(double x0, double y0, double z0, const vector<double>& flat_point_cloud_data, double& limax, Qhull& qt){
    double lmin=999;
    vector<double> p;
    p.push_back(x0);
    p.push_back(y0);
    p.push_back(z0);
    
    
    qt.runQhull("obj", 3, flat_point_cloud_data.size()/3, flat_point_cloud_data.data(), "");
    //  cout << qt.facetList();
    for(auto it = qt.facetList().begin();it!=qt.facetList().end();it++){
        //        auto v = it->getCenter().toStdVector();
        auto d = it->distance(p.data());
        // cout<<d<<endl;
        if(abs(d)<=lmin){
            lmin=abs(d);
        }
        if(d>=0){
            DebugOn("centre outside convex hull");
            lmin=0;
            break;
        }
    }
    
    DebugOn("minimum radius found "<<lmin<<endl);
    
    vector<vector<double>> sphere;
    vector<double> radius;
    
    radius.resize(3);
    double r=lmin/10.0;
    for (auto i=-10;i<10;i++){
        auto xs=i*r;
        for(auto j=-10;j<10;j++){
            auto ys=j*r;
            auto zr=pow(lmin,2)-pow(ys,2)-pow(xs,2);
            if(zr>0){
                auto zs=sqrt(zr);
                DebugOff(xs<<" "<<ys<<" "<<zs<<endl);
                radius[0]=(xs);
                radius[1]=(ys);
                radius[2]=(zs);
                sphere.push_back(radius);
                radius[0]=(xs);
                radius[1]=ys;
                radius[2]=(zs*(-1));
                DebugOff(xs<<" "<<ys<<" "<<zs*(-1)<<endl);
                sphere.push_back(radius);
            }
        }
    }
    for (auto i=-10;i<10;i++){
        auto zs=i*r;
        for(auto j=-10;j<10;j++){
            auto ys=j*r;
            auto zr=pow(lmin,2)-pow(ys,2)-pow(zs,2);
            if(zr>0){
                auto xs=sqrt(zr);
                DebugOff(xs<<" "<<ys<<" "<<zs<<endl);
                radius[0]=(xs);
                radius[1]=(ys);
                radius[2]=(zs);
                sphere.push_back(radius);
                radius[0]=(xs*(-1));
                radius[1]=ys;
                radius[2]=(zs);
                DebugOff(xs<<" "<<ys<<" "<<zs*(-1)<<endl);
                sphere.push_back(radius);
            }
        }
    }
    
    
    DebugOn(sphere.size()<<endl);
    
    vector<vector<double>> point_cloud_data;
    for(auto i=0;i<flat_point_cloud_data.size();i+=3){
        vector<double> r;
        r.push_back(flat_point_cloud_data[i]);
        r.push_back(flat_point_cloud_data[i+1]);
        r.push_back(flat_point_cloud_data[i+2]);
        point_cloud_data.push_back(r);
    }
    limax= maximum_inscribed_sphere_all_points( point_cloud_data,  qt);
#ifdef USE_MATPLOT
    //plot(sphere,point_cloud_data,5);
#endif
    return lmin;
    
}

double  maximum_inscribed_sphere_all_points(vector<vector<double>> point_cloud_data,  Qhull& qt){
    double lmax=-1.0, lmini;
    vector<double> p;
    
    for(auto i=0;i<point_cloud_data.size();i++){
        lmini=999;
        p=point_cloud_data.at(i);
        for(auto it = qt.facetList().begin();it!=qt.facetList().end();it++){
            auto d = it->distance(p.data());
            // cout<<d<<endl;
            if(abs(d)<=lmini){
                lmini=abs(d);
            }
            if(d>=1e-6){
                DebugOn("centre outside convex hull");
                lmini=0;
                break;
            }
        }
        if(lmini>=lmax){
            lmax=lmini;
        }
    }
    return lmax;
}
#endif


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

shared_ptr<Model<double>> build_projected_SOC_MIQCP(vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, vector<double>& rot_trans, bool convex, const vector<pair<pair<int,int>,pair<int,int>>>& incompatibles, param<>& norm_x, param<>& norm_y, param<>& norm_z,  param<>& intercept, const vector<int>& init_matching){}


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

shared_ptr<Model<double>> build_norm2_SOC_MIQCP(vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, const indices& valid_cells, const indices& new_model_ids, const param<>& dist_cost, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double shift_min_x, double shift_max_x, double shift_min_y, double shift_max_y, double shift_min_z, double shift_max_z, vector<double>& rot_trans, bool convex, const vector<pair<pair<int,int>,pair<int,int>>>& incompatibles, param<>& norm_x, param<>& norm_y, param<>& norm_z,  param<>& intercept, const vector<int>& init_matching, const vector<double>& error_per_point, param<>& model_radius, bool relax_ints, bool relax_sdp, bool rigid_transf, double perc_outliers){
    double roll_1 = 0, yaw_1 = 0, pitch_1 = 0;
    int nb_pairs = 0, min_nb_pairs = numeric_limits<int>::max(), max_nb_pairs = 0, av_nb_pairs = 0;
    size_t nm = point_cloud_model.size(), nd = point_cloud_data.size();
    double xm_max = numeric_limits<double>::lowest(), ym_max = numeric_limits<double>::lowest(), zm_max = numeric_limits<double>::lowest();
    double xm_min = numeric_limits<double>::max(), ym_min = numeric_limits<double>::max(), zm_min = numeric_limits<double>::max();
    vector<double> zeros = {0,0,0};
    
    param<> x1("x1"), x2("x2"), y1("y1"), y2("y2"), z1("z1"), z2("z2");
    int m = av_nb_pairs;
    string i_str, j_str;
    for (auto i = 0; i<nd; i++) {
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
    //    for (const auto &model_key: *new_model_ids._keys) {
    //        int j = stoi(model_key);
    //        x2.add_val(model_key,point_cloud_model.at(j-1).at(0));
    //        if(point_cloud_model.at(j-1).at(0) > xm_max)
    //            xm_max = point_cloud_model.at(j-1).at(0);
    //        if(point_cloud_model.at(j-1).at(0) < xm_min)
    //            xm_min = point_cloud_model.at(j-1).at(0);
    //        y2.add_val(model_key,point_cloud_model.at(j-1).at(1));
    //        if(point_cloud_model.at(j-1).at(1) > ym_max)
    //            ym_max = point_cloud_model.at(j-1).at(1);
    //        if(point_cloud_model.at(j-1).at(1) < ym_min)
    //            ym_min = point_cloud_model.at(j-1).at(1);
    //        z2.add_val(model_key,point_cloud_model.at(j-1).at(2));
    //        if(point_cloud_model.at(j-1).at(2) > zm_max)
    //            zm_max = point_cloud_model.at(j-1).at(2);
    //        if(point_cloud_model.at(j-1).at(2) < zm_min)
    //            zm_min = point_cloud_model.at(j-1).at(2);
    //    }
    
    
    indices Pairs("Pairs"), cells("cells");
    int idx1 = 0;
    int idx2 = 0;
    indices N1("N1"),N2("N2");
    Debug("nd = " << nd << endl);
    Debug("nm = " << nm << endl);
    
    N1 = range(1,nd);
    N2 = range(1,nm);
    cells = valid_cells;
    
    //    cells = indices(N1,N2);
    string name="Norm2_MISDP";
    
    auto Reg=make_shared<Model<>>(name);
    
    var<> scale_x("scale_x", 0.6, 1.4);
    var<> scale_y("scale_y", 0.6, 1.4);
    var<> scale_z("scale_z", 0.6, 1.4);
    if(!rigid_transf){
        Reg->add(scale_x.in(R(1)),scale_y.in(R(1)),scale_z.in(R(1)));
        rot_trans.resize(15);
    }
    else{
        rot_trans.resize(12);
    }
    
    Reg->add_param(model_radius);
    Reg->add_param(x1);Reg->add_param(y1);Reg->add_param(z1);
    Reg->add_param(x2);Reg->add_param(y2);Reg->add_param(z2);
    var<> x_shift("x_shift", shift_min_x, shift_max_x), y_shift("y_shift", shift_min_y, shift_max_y), z_shift("z_shift", shift_min_z, shift_max_z);
    
    Reg->add(x_shift.in(R(1)),y_shift.in(R(1)),z_shift.in(R(1)));
    
    var<int> bin("bin",0,1);
    Reg->add(bin.in(cells));
    DebugOff("Added " << cells.size() << " binary variables" << endl);
    
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
    
    
    Reg->add(theta11.in(R(1)),theta12.in(R(1)),theta13.in(R(1)));
    Reg->add(theta21.in(R(1)),theta22.in(R(1)),theta23.in(R(1)));
    Reg->add(theta31.in(R(1)),theta32.in(R(1)),theta33.in(R(1)));
    
    param<> x_new_lb("x_new_lb");
    x_new_lb.in(N1);
    param<> x_new_ub("x_new_ub");
    x_new_ub.in(N1);
    param<> y_new_lb("y_new_lb");
    y_new_lb.in(N1);
    param<> y_new_ub("y_new_ub");
    y_new_ub.in(N1);
    param<> z_new_lb("z_new_lb");
    z_new_lb.in(N1);
    param<> z_new_ub("z_new_ub");
    z_new_ub.in(N1);
    
    param<> mid_point_lb("mid_point_lb");
    mid_point_lb.in(N1);
    
    double lower_bound = 0, go_icp_lb = 0, SDP_lb = 0;
    double x_lb = 0, y_lb = 0, z_lb = 0, x1_i = 0, y1_i = 0, z1_i = 0;
    shared_ptr<pair<double,double>> x1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> y1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> z1_bounds = make_shared<pair<double,double>>();
    for (int i = 0; i<nd; i++) {
        x1_bounds->first = x1.eval(i);
        x1_bounds->second = x1.eval(i);
        y1_bounds->first = y1.eval(i);
        y1_bounds->second = y1.eval(i);
        z1_bounds->first = z1.eval(i);
        z1_bounds->second = z1.eval(i);
        auto x_range  = get_product_range(x1_bounds, theta11._range);
        auto y_range  = get_product_range(y1_bounds, theta12._range);
        auto z_range  = get_product_range(z1_bounds, theta13._range);
        auto bounds = get_min_max(roll_min, roll_max, pitch_min, pitch_max, yaw_min, yaw_max, point_cloud_data[i], zeros);
        
        //        double new_roll_max1 = acos(sqrt(cos(sqrt(3)*roll_max)*(2./3.) + 1./3.));
        //        double new_roll_max = asin(sin(sqrt(3.)*roll_max)/sqrt(3.) + 1./3.*(1 - cos(sqrt(3.)*roll_max)));
        
        //        auto go_icp_max_dist = get_GoICP_dist(roll_max, shift_max_x, point_cloud_data[i], false);
        
        //        double v1 = 1./std::sqrt(3.), v2 = 1./std::sqrt(3.), v3 = 1./std::sqrt(3.), ct = std::cos(std::sqrt(3.)*roll_max), ct2 = (1.-std::cos(std::sqrt(3.)*roll_max)), st = std::sin(std::sqrt(3.)*roll_max);
        //        double tmp231 = v2*v3*ct2, tmp232 = v1*st;
        //        auto tmp121 = v2*v2*ct2;
        //        auto tmp122 = v3*st;
        //        auto tmp131 = v1*v3*ct2;
        //        auto tmp132 = v2*st;
        //        auto R11 = ct + v1*v1*ct2;
        //        auto R12 = tmp121 - tmp122;
        //        auto R13 = tmp131 + tmp132;
        //        auto R21 = tmp121 + tmp122;
        //        auto R22 = ct + v2*v2*ct2;
        //        auto R23 = tmp231 - tmp232;
        //        double R31 = tmp131 - tmp132;
        //        double R32 = tmp231 + tmp232;
        //        double R33 = ct + v3*v3*ct2;
        //        vector<double> rot_mat = {R11, R12, R13, R21, R22, R23, R31, R32, R33};
        //
        ////        auto SDP = build_SDP(point_cloud_data[i], rot_mat);
        ////        SDP->is_feasible(1e-6);
        ////        double max_sdp = std::sqrt(SDP->get_obj_val());
        //
        //        double new_pitch_max = std::atan2(R32, R33);
        //        auto new_yaw_max = std::atan2(R21,R11);
        //        auto new_roll_max = std::atan2(-1*R31, std::sqrt(R32*R32+R33*R33));
        
        auto go_icp_max_dist = get_GoICP_dist(roll_max, shift_max_x, point_cloud_data[i], false);
        //        DebugOn("SDP max dist = " << to_string_with_precision(max_sdp, 6) << endl);
        //        auto max_dist = get_max_dist(-new_roll_max, new_roll_max, -new_pitch_max, new_pitch_max, -new_yaw_max, new_yaw_max, shift_min_x, shift_max_x, shift_min_y, shift_max_y, shift_min_z, shift_max_z, point_cloud_data[i], zeros, false);
        
        //        lower_bound += std::max(0.,error_per_point[i] - max_dist);
        go_icp_lb += std::max(0.,error_per_point[i] - go_icp_max_dist);
        //        SDP_lb += std::max(0.,error_per_point[i] - max_sdp);
        //        mid_point_lb.set_val(i,std::max(0.,error_per_point[i] - max_dist));
        x_range->first = x_range->first + y_range->first + z_range->first;
        x_range->second = x_range->second + y_range->second + z_range->second;
        //        x_new_lb.set_val(i,bounds[0].first + x_shift.get_lb().eval());
        //        x_new_ub.set_val(i,bounds[0].second + x_shift.get_ub().eval());
        auto new_range = get_product_range(x_range, scale_x._range);
        x_new_lb.set_val(i, new_range->first);
        x_new_ub.set_val(i, new_range->second);
        x_range  = get_product_range(x1_bounds, theta21._range);
        y_range  = get_product_range(y1_bounds, theta22._range);
        z_range  = get_product_range(z1_bounds, theta23._range);
        y_range->first = x_range->first + y_range->first + z_range->first;
        y_range->second = x_range->second + y_range->second + z_range->second;
        
        new_range = get_product_range(y_range, scale_y._range);
        //        auto ylb = x_range->first + y_range->first + z_range->first + y_shift.get_lb().eval();
        //        auto yub = x_range->second + y_range->second + z_range->second + y_shift.get_ub().eval();
        //        y_new_lb.set_val(i,bounds[1].first + y_shift.get_lb().eval());
        //        y_new_ub.set_val(i,bounds[1].second + y_shift.get_ub().eval());
        y_new_lb.set_val(i, new_range->first);
        y_new_ub.set_val(i, new_range->second);
        x_range  = get_product_range(x1_bounds, theta31._range);
        y_range  = get_product_range(y1_bounds, theta32._range);
        z_range  = get_product_range(z1_bounds, theta33._range);
        z_range->first = x_range->first + y_range->first + z_range->first;
        z_range->second = x_range->second + y_range->second + z_range->second;
        
        new_range = get_product_range(z_range, scale_z._range);
        //        auto zlb = x_range->first + y_range->first + z_range->first + z_shift.get_lb().eval();
        //        auto zub = x_range->second + y_range->second + z_range->second + z_shift.get_ub().eval();
        //        z_new_lb.set_val(i,bounds[2].first + z_shift.get_lb().eval());
        //        z_new_ub.set_val(i,bounds[2].second + z_shift.get_ub().eval());
        z_new_lb.set_val(i, new_range->first);
        z_new_ub.set_val(i, new_range->second);
    }
    //    DebugOn("Lower bound = " << to_string_with_precision(lower_bound,6) << endl);
    DebugOff("GoICP lower bound = " << to_string_with_precision(go_icp_lb,6) << endl);
    //    DebugOn("SDP lower bound = " << to_string_with_precision(SDP_lb,6) << endl);
    
    var<> new_xm("new_xm", xm_min, xm_max), new_ym("new_ym", ym_min, ym_max), new_zm("new_zm", zm_min, zm_max);
    var<> new_x1("new_x1", x_new_lb, x_new_ub), new_y1("new_y1", y_new_lb, y_new_ub), new_z1("new_z1", z_new_lb, z_new_ub);
    //            var<> new_x1("new_x1", -1, 1), new_y1("new_y1", -1, 1), new_z1("new_z1", -1, 1);
    Reg->add(new_xm.in(N1), new_ym.in(N1), new_zm.in(N1));
    if(!rigid_transf){
        Reg->add(new_x1.in(N1), new_y1.in(N1), new_z1.in(N1));
    }
    
    //    var<> new_x1_ij("new_x1_ij", -1, 1), new_y1_ij("new_y1_ij", -1, 1), new_z1_ij("new_z1_ij", -1, 1);
    //    Reg->add(new_x1_ij.in(cells), new_y1_ij.in(cells), new_z1_ij.in(cells));
    
    
    var<> delta("delta", 0,12);
    //    Reg->add(delta.in(N1));
    
    
    //    var<> d1("d1", 0,4),d2("d2", 0,4),d3("d3", 0,4),d4("d4", 0,4);
    //    var<> l12("l12", -2,2),l13("l13", -2,2),l14("l14", -2,2),l23("l23", -2,2),l24("l24", -2,2),l34("l34", -2,2);
    //    Reg->add(l12.in(R(1)),l13.in(R(1)),l14.in(R(1)));
    //    Reg->add(l23.in(R(1)),l24.in(R(1)),l34.in(R(1)));
    //    Reg->add(d1.in(R(1)),d2.in(R(1)),d3.in(R(1)),d4.in(R(1)));
    
    
    
    
    indices ids = indices("in_x");
    ids = N2;
    ids.add_empty_row();
    for(auto i=0;i<nd;i++){
        for(auto j=1;j<=nm;j++){
            string model_key = to_string(j);
            if(cells.has_key(to_string(i+1)+","+model_key)){
                ids.add_in_row(i, model_key);
            }
        }
        if(ids._ids->at(i).size()==0){
            Reg->_status = -1;
            return Reg;
        }
    }
    
    
    
    bool add_voronoi = false;
    
    if(add_voronoi){
        indices voronoi_ids("voronoi_ids");
        for (const auto &key: *norm_x._indices->_keys) {
            auto m_id = key.substr(0, key.find_first_of(","));
            if(new_model_ids.has_key(m_id)){
                for(auto i=0;i<3;i++){
                    if(cells.has_key(to_string(i+1)+","+m_id)){
                        voronoi_ids.insert(to_string(i+1)+","+key);
                    }
                }
            }
        }
        auto voronoi_ids_coefs = voronoi_ids.ignore_ith(0, 1);
        auto voronoi_ids_m = voronoi_ids.from_ith(0, 1);
        auto voronoi_ids_data = voronoi_ids.ignore_ith(1, 2);
        auto voronoi_ids_bin = voronoi_ids.ignore_ith(2, 1);
        auto ids1 = theta11.repeat_id(voronoi_ids.size());
        Constraint<> Voronoi("Voronoi");
        Voronoi = norm_x.in(voronoi_ids_coefs)*(x1.in(voronoi_ids_data)*theta11.in(ids1) + y1.in(voronoi_ids_data)*theta12.in(ids1) + z1.in(voronoi_ids_data)*theta13.in(ids1)+x_shift.in(ids1)) + norm_y.in(voronoi_ids_coefs)*(x1.in(voronoi_ids_data)*theta21.in(ids1) + y1.in(voronoi_ids_data)*theta22.in(ids1) + z1.in(voronoi_ids_data)*theta23.in(ids1)+y_shift.in(ids1)) + norm_z.in(voronoi_ids_coefs)*(x1.in(voronoi_ids_data)*theta31.in(ids1) + y1.in(voronoi_ids_data)*theta32.in(ids1) + z1.in(voronoi_ids_data)*theta33.in(ids1)+z_shift.in(ids1)) + intercept.in(voronoi_ids_coefs);
        Reg->add_on_off_multivariate_refined(Voronoi.in(voronoi_ids)<=0, bin.in(voronoi_ids_bin), true);
        //        Reg->print();
    }
    
    
    theta11.initialize_all(1);
    theta22.initialize_all(1);
    theta33.initialize_all(1);
    /* REF */
    bool spatial_branching = false;
    if(spatial_branching){
        var<> yaw("yaw", yaw_min, yaw_max), pitch("pitch", pitch_min, pitch_max), roll("roll", roll_min, roll_max);
        yaw.in(R(1)); pitch.in(R(1));roll.in(R(1));
        func<> cosr_f = cos(roll);
        func<> sinr_f = sin(roll);
        func<> cosp_f = cos(pitch);
        func<> sinp_f = sin(pitch);
        func<> cosy_f = cos(yaw);
        func<> siny_f = sin(yaw);
        var<> cosr("cosr", cosr_f._range->first, cosr_f._range->second), sinr("sinr", sinr_f._range->first, sinr_f._range->second);
        var<> cosp("cosp",  cosp_f._range->first, cosp_f._range->second), sinp("sinp", sinp_f._range->first, sinp_f._range->second);
        var<> cosy("cosy",  cosy_f._range->first, cosy_f._range->second), siny("siny", siny_f._range->first, siny_f._range->second);
        auto cosy_sinr_range = get_product_range(cosy._range, sinr._range);
        auto siny_sinr_range = get_product_range(siny._range, sinr._range);
        var<> cosy_sinr("cosy_sinr", cosy_sinr_range->first, cosy_sinr_range->second), siny_sinr("siny_sinr", siny_sinr_range->first, siny_sinr_range->second);
        
        Reg->add(cosr.in(R(1)),cosp.in(R(1)),cosy.in(R(1)));
        Reg->add(sinr.in(R(1)),sinp.in(R(1)),siny.in(R(1)));
        Reg->add(cosy_sinr.in(R(1)),siny_sinr.in(R(1)));
        
        Constraint<> cosy_sinr_prod("cosy_sinr");
        cosy_sinr_prod = cosy_sinr - cosy*sinr;
        Reg->add(cosy_sinr_prod==0);
        
        Constraint<> siny_sinr_prod("siny_sinr");
        siny_sinr_prod = siny_sinr - siny*sinr;
        Reg->add(siny_sinr_prod==0);
        
        
        Constraint<> R11("R11");
        R11 += theta11 - cosy*cosr;
        Reg->add(R11==0);
        
        Constraint<> R12("R12");
        R12 += theta12 - (cosy_sinr*sinp - siny*cosp);
        Reg->add(R12==0);
        
        Constraint<> R13("R13");
        R13 += theta13 - (cosy_sinr*cosp + siny*sinp);
        Reg->add(R13==0);
        
        Constraint<> R21("R21");
        R21 += theta21 - siny*cosr;
        Reg->add(R21==0);
        
        Constraint<> R22("R22");
        R22 += theta22 - (siny_sinr*sinp + cosy*cosp);
        Reg->add(R22==0);
        
        Constraint<> R23("R23");
        R23 += theta23 - (siny_sinr*cosp - cosy*sinp);
        Reg->add(R23==0);
        
        Constraint<> R31("R31");
        R31 += theta31 + sinr;
        Reg->add(R31==0);
        
        Constraint<> R32("R32");
        R32 += theta32 - cosr*sinp;
        Reg->add(R32==0);
        
        Constraint<> R33("R33");
        R33 += theta33 - cosr*cosp;
        Reg->add(R33==0);
        
        /*
         func<> r11 = cos(yaw)*cos(roll);r11.eval_all();
         func<> r12 = cos(yaw)*sin(roll)*sin(pitch) - sin(yaw)*cos(pitch);r12.eval_all();
         func<> r13 = cos(yaw)*sin(roll)*cos(pitch) + sin(yaw)*sin(pitch);r13.eval_all();
         func<> r21 = sin(yaw)*cos(roll);r21.eval_all();
         func<> r22 = sin(yaw)*sin(roll)*sin(pitch) + cos(yaw)*cos(pitch);r22.eval_all();
         func<> r23 = sin(yaw)*sin(roll)*cos(pitch) - cos(yaw)*sin(pitch);r23.eval_all();
         func<> r31 = sin(-1*roll);r31.eval_all();
         func<> r32 = cos(roll)*sin(pitch);r32.eval_all();
         func<> r33 = cos(roll)*cos(pitch);r33.eval_all();
         
         */
        
        bool add_spatial_bins = true;
        
        if(add_spatial_bins){
            /* Spatial branching vars */
            int nb_pieces = 10; // Divide each axis into nb_pieces
            indices spatial_ids("spatial_ids");
            spatial_ids = range(1,nb_pieces);
            indices theta_ids("theta_ids");
            theta_ids = indices(range(1,3),range(1,3));
            indices shift_ids("shift_ids");
            shift_ids.insert({"x", "y", "z"});
            indices theta_spatial("theta_spatial");
            theta_spatial = indices(theta_ids, spatial_ids);
            
            indices shift_spatial("shift_spatial");
            shift_spatial = indices(shift_ids, spatial_ids);
            
            //    var<int> sbin_shift("sbin_shift", 0, 1);
            //    var<int> sbin_theta("sbin_theta", 0, 1);
            
            //    Reg->add(sbin_theta.in(theta_spatial),sbin_shift.in(shift_spatial));
            
            var<int> sbin_tx("sbin_tx", 0, 1), sbin_ty("sbin_ty", 0, 1), sbin_tz("sbin_tz", 0, 1);
            var<int> sbin_roll("sbin_roll", 0, 1), sbin_pitch("sbin_pitch", 0, 1), sbin_yaw("sbin_yaw", 0, 1);
            sbin_roll._priority = 1e6;sbin_pitch._priority = 1e6;sbin_yaw._priority = 1e6;
            sbin_tx._priority = 2e6;sbin_ty._priority = 2e6;sbin_tz._priority = 2e6;
            Reg->add(sbin_roll.in(spatial_ids),sbin_pitch.in(spatial_ids),sbin_yaw.in(spatial_ids));
            Reg->add(sbin_tz.in(spatial_ids), sbin_tx.in(spatial_ids),sbin_ty.in(spatial_ids));
            /* Spatial branching constraints */
            Constraint<> OneBinRoll("OneBinRoll");
            OneBinRoll = sum(sbin_roll);
            Reg->add(OneBinRoll==1);
            
            Constraint<> OneBinPitch("OneBinPitch");
            OneBinPitch = sum(sbin_pitch);
            Reg->add(OneBinPitch==1);
            
            Constraint<> OneBinYaw("OneBinYaw");
            OneBinYaw = sum(sbin_yaw);
            Reg->add(OneBinYaw==1);
            
            
            Constraint<> OneBinShiftSpatialx("OneBinShiftSpatialx");
            OneBinShiftSpatialx = sum(sbin_tx);
            Reg->add(OneBinShiftSpatialx==1);
            
            Constraint<> OneBinShiftSpatialy("OneBinShiftSpatialy");
            OneBinShiftSpatialy = sum(sbin_ty);
            Reg->add(OneBinShiftSpatialy==1);
            
            Constraint<> OneBinShiftSpatialz("OneBinShiftSpatialz");
            OneBinShiftSpatialz = sum(sbin_tz);
            Reg->add(OneBinShiftSpatialz==1);
            
            //    Reg->print();
            double angle_min = std::min(std::min(roll_min,pitch_min),yaw_min), angle_max = std::max(std::max(roll_max,pitch_max),yaw_max);
            double t_min = std::min(std::min(shift_min_x,shift_min_y),shift_min_z), t_max = std::max(std::max(shift_max_x,shift_max_y),shift_max_z);
            double angle_increment = (angle_max - angle_min)/nb_pieces;/* Angles are defined in [-2,2] in radians (-120,120) in degrees */
            double shift_increment = (t_max - t_min)/nb_pieces;/* Shifts are defined in [-0.5,0.5] */
            
            auto spatial_ids_n = range(1,nb_pieces-1);
            auto spatial_ids_1 = range(2,nb_pieces);
            param<> angle_lb("angle_lb"), angle_ub("angle_ub");
            param<> cos_lb("cos_lb"), cos_ub("cos_ub");
            param<> sin_lb("sin_lb"), sin_ub("sin_ub");
            param<> t_lb("t_lb"), t_ub("t_ub");
            angle_ub.in(spatial_ids);
            angle_lb.in(spatial_ids);
            cos_ub.in(spatial_ids);
            cos_lb.in(spatial_ids);
            sin_ub.in(spatial_ids);
            sin_lb.in(spatial_ids);
            t_ub.in(spatial_ids);
            t_lb.in(spatial_ids);
            double lb = 0, ub = 0;
            for (int i = 0; i<nb_pieces; i++) {
                lb = angle_min + i*angle_increment;
                ub = angle_min+(i+1)*angle_increment;
                angle_lb.set_val(i,lb);
                angle_ub.set_val(i,ub);
                cos_lb.set_val(i,std::min(cos(lb), cos(ub)));
                cos_ub.set_val(i,std::max(cos(lb), cos(ub)));
                if(lb < 0 && ub > 0){/* zero is in the domain, i.e., cos max is 1 */
                    cos_ub.set_val(i,1);
                }
                if((lb < -pi && ub > -pi) || (lb < pi && ub > pi)){/* -pi or pi is in the domain, i.e., cos min is -1 */
                    cos_lb.set_val(i,-1);
                }
                sin_lb.set_val(i,std::min(sin(lb), sin(ub)));
                sin_ub.set_val(i,std::max(sin(lb), sin(ub)));
                if((lb < -3*pi/2 && ub > -3*pi/2) || (lb < pi/2 && ub > pi/2)){/* -3pi/2 or pi/2 is in the domain, i.e., sin max is 1 */
                    sin_ub.set_val(i,1);
                }
                if((lb < 3*pi/2 && ub > 3*pi/2) || (lb < -pi/2 && ub > -pi/2)){/* 3pi/2 or -pi/2 is in the domain, i.e., sin min is -1 */
                    sin_lb.set_val(i,-1);
                }
                t_lb.set_val(i,t_min + i*shift_increment);
                t_ub.set_val(i,t_min + (i+1)*shift_increment);
            }
            Reg->add_param(angle_lb);
            Reg->add_param(angle_ub);
            Reg->add_param(cos_lb);
            Reg->add_param(cos_ub);
            Reg->add_param(sin_lb);
            Reg->add_param(sin_ub);
            Reg->add_param(t_lb);
            Reg->add_param(t_ub);
            
            auto ids_repeat = theta11.repeat_id(nb_pieces-1);
            
            
            Constraint<> CosRoll_UB("CosRoll_UB");
            CosRoll_UB = cosr.in(ids_repeat) - cos_ub.in(spatial_ids_n);
            Reg->add_on_off_multivariate_refined(CosRoll_UB.in(spatial_ids_n)<=0, sbin_roll.in(spatial_ids_n), true);
            
            Constraint<> CosRoll_LB("CosRoll_LB");
            CosRoll_LB = cos_lb.in(spatial_ids_1) - cosr.in(ids_repeat);
            Reg->add_on_off_multivariate_refined(CosRoll_LB.in(spatial_ids_1) <= 0, sbin_roll.in(spatial_ids_1), true);
            
            Constraint<> SinRoll_UB("SinRoll_UB");
            SinRoll_UB = sinr.in(ids_repeat) - sin_ub.in(spatial_ids_n);
            Reg->add_on_off_multivariate_refined(SinRoll_UB.in(spatial_ids_n)<=0, sbin_roll.in(spatial_ids_n), true);
            
            Constraint<> SinRoll_LB("SinRoll_LB");
            SinRoll_LB = sin_lb.in(spatial_ids_1) - sinr.in(ids_repeat);
            Reg->add_on_off_multivariate_refined(SinRoll_LB.in(spatial_ids_1) <= 0, sbin_roll.in(spatial_ids_1), true);
            
            Constraint<> CosPitch_UB("CosPitch_UB");
            CosPitch_UB = cosp.in(ids_repeat) - cos_ub.in(spatial_ids_n);
            Reg->add_on_off_multivariate_refined(CosPitch_UB.in(spatial_ids_n)<=0, sbin_pitch.in(spatial_ids_n), true);
            
            Constraint<> CosPitch_LB("CosPitch_LB");
            CosPitch_LB = cos_lb.in(spatial_ids_1) - cosp.in(ids_repeat);
            Reg->add_on_off_multivariate_refined(CosPitch_LB.in(spatial_ids_1) <= 0, sbin_pitch.in(spatial_ids_1), true);
            
            Constraint<> SinPitch_UB("SinPitch_UB");
            SinPitch_UB = sinp.in(ids_repeat) - sin_ub.in(spatial_ids_n);
            Reg->add_on_off_multivariate_refined(SinPitch_UB.in(spatial_ids_n)<=0, sbin_pitch.in(spatial_ids_n), true);
            
            Constraint<> SinPitch_LB("SinPitch_LB");
            SinPitch_LB = sin_lb.in(spatial_ids_1) - sinp.in(ids_repeat);
            Reg->add_on_off_multivariate_refined(SinPitch_LB.in(spatial_ids_1) <= 0, sbin_pitch.in(spatial_ids_1), true);
            
            
            Constraint<> CosYaw_UB("CosYaw_UB");
            CosYaw_UB = cosy.in(ids_repeat) - cos_ub.in(spatial_ids_n);
            Reg->add_on_off_multivariate_refined(CosYaw_UB.in(spatial_ids_n)<=0, sbin_yaw.in(spatial_ids_n), true);
            
            Constraint<> CosYaw_LB("CosYaw_LB");
            CosYaw_LB = cos_lb.in(spatial_ids_1) - cosy.in(ids_repeat);
            Reg->add_on_off_multivariate_refined(CosYaw_LB.in(spatial_ids_1) <= 0, sbin_yaw.in(spatial_ids_1), true);
            
            Constraint<> SinYaw_UB("SinYaw_UB");
            SinYaw_UB = siny.in(ids_repeat) - sin_ub.in(spatial_ids_n);
            Reg->add_on_off_multivariate_refined(SinYaw_UB.in(spatial_ids_n)<=0, sbin_yaw.in(spatial_ids_n), true);
            
            Constraint<> SinYaw_LB("SinYaw_LB");
            SinYaw_LB = sin_lb.in(spatial_ids_1) - siny.in(ids_repeat);
            Reg->add_on_off_multivariate_refined(SinYaw_LB.in(spatial_ids_1) <= 0, sbin_yaw.in(spatial_ids_1), true);
            
            Constraint<> xshift_Spatial_UB("xshift_Spatial_UB");
            xshift_Spatial_UB = x_shift.in(ids_repeat) - t_ub.in(spatial_ids_n);
            Reg->add_on_off_multivariate_refined(xshift_Spatial_UB.in(spatial_ids_n) <= 0, sbin_tx.in(spatial_ids_n), true);
            
            Constraint<> xshift_Spatial_LB("xshift_Spatial_LB");
            xshift_Spatial_LB = t_lb.in(spatial_ids_1) - x_shift.in(ids_repeat);
            Reg->add_on_off_multivariate_refined(xshift_Spatial_LB.in(spatial_ids_1) <= 0, sbin_tx.in(spatial_ids_1), true);
            
            Constraint<> yshift_Spatial_UB("yshift_Spatial_UB");
            yshift_Spatial_UB = y_shift.in(ids_repeat) - t_ub.in(spatial_ids_n);
            Reg->add_on_off_multivariate_refined(yshift_Spatial_UB.in(spatial_ids_n) <= 0, sbin_ty.in(spatial_ids_n), true);
            
            Constraint<> yshift_Spatial_LB("yshift_Spatial_LB");
            yshift_Spatial_LB = t_lb.in(spatial_ids_1) - y_shift.in(ids_repeat);
            Reg->add_on_off_multivariate_refined(yshift_Spatial_LB.in(spatial_ids_1) <= 0, sbin_ty.in(spatial_ids_1), true);
            
            Constraint<> zshift_Spatial_UB("zshift_Spatial_UB");
            zshift_Spatial_UB = z_shift.in(ids_repeat) - t_ub.in(spatial_ids_n);
            Reg->add_on_off_multivariate_refined(zshift_Spatial_UB.in(spatial_ids_n) <= 0, sbin_tz.in(spatial_ids_n), true);
            
            Constraint<> zshift_Spatial_LB("zshift_Spatial_LB");
            zshift_Spatial_LB = t_lb.in(spatial_ids_1) - z_shift.in(ids_repeat);
            Reg->add_on_off_multivariate_refined(zshift_Spatial_LB.in(spatial_ids_1) <= 0, sbin_tz.in(spatial_ids_1), true);
        }
        //
        Constraint<> trigR("trigR");
        trigR = pow(cosr,2) + pow(sinr,2);
        trigR.add_to_callback();
        Reg->add(trigR.in(range(1,1))<=1);
        
        Constraint<> trigP("trigP");
        trigP = pow(cosp,2) + pow(sinp,2);
        trigP.add_to_callback();
        Reg->add(trigP.in(range(1,1))<=1);
        
        Constraint<> trigY("trigY");
        trigY = pow(cosy,2) + pow(siny,2);
        trigY.add_to_callback();
        Reg->add(trigY.in(range(1,1))<=1);
        
        
        Constraint<> trigR_NC("trigR_NC");
        trigR_NC = pow(cosr,2) + pow(sinr,2);
        Reg->add(trigR_NC.in(range(1,1))>=1);
        
        Constraint<> trigP_NC("trigP_NC");
        trigP_NC = pow(cosp,2) + pow(sinp,2);
        Reg->add(trigP_NC.in(range(1,1))>=1);
        
        Constraint<> trigY_NC("trigY_NC");
        trigY_NC = pow(cosy,2) + pow(siny,2);
        Reg->add(trigY_NC.in(range(1,1))>=1);
        //        Reg->print();
    }
    
    
    Constraint<> Def_newxm("Def_newxm");
    Def_newxm = new_xm-product(x2.in(ids),bin.in_matrix(1, 1));
    Reg->add(Def_newxm.in(N1)==0);
    
    Constraint<> Def_newym("Def_newym");
    Def_newym = new_ym-product(y2.in(ids),bin.in_matrix(1, 1));
    Reg->add(Def_newym.in(N1)==0);
    
    Constraint<> Def_newzm("Def_newzm");
    Def_newzm = new_zm-product(z2.in(ids),bin.in_matrix(1, 1));
    Reg->add(Def_newzm.in(N1)==0);
    
    var<int> keep("keep",0,1);
    
    if(perc_outliers==0){
        Constraint<> OneBin("OneBin");
        OneBin = bin.in_matrix(1, 1);
        Reg->add(OneBin.in(N1)==1);
    }
    else{
        Reg->add(keep.in(N1));
        Constraint<> OneBin("OneBin");
        OneBin = bin.in_matrix(1, 1) - keep;
        Reg->add(OneBin.in(N1)==0);
        
        Constraint<> OneBinMin("OneBinMin");
        OneBinMin = nd*(100-perc_outliers)/100. - sum(keep);
        Reg->add(OneBinMin<=0);
        //        Reg->print();
    }
    
    //    Constraint<> OneBin2("OneBin2");
    //    OneBin2 = bin.in_matrix(0, 1);
    //    Reg->add(OneBin2.in(N2)<=1);
    
    if(!incompatibles.empty()){
        indices pairs1("pairs1"), pairs2("pairs2");
        pairs1 = cells;
        pairs2 = cells;
        for (const auto &inc_pair : incompatibles) {
            string key1 = to_string(inc_pair.first.first+1)+","+to_string(inc_pair.second.first+1);
            string key2 = to_string(inc_pair.first.second+1)+","+to_string(inc_pair.second.second+1);
            if(cells.has_key(key1) && cells.has_key(key2)){
                pairs1.add_ref(key1);
                pairs2.add_ref(key2);
            }
        }
        
        if(pairs1.is_indexed()){
            DebugOn("Number of incompatible pairs constraints = " << pairs1.size() << endl);
            Constraint<> incomp_pairs("incomp_pairs");
            incomp_pairs = bin.in(pairs1) + bin.in(pairs2);
            Reg->add(incomp_pairs.in(range(1,pairs1.size()))<=1);
        }
        //        incomp_pairs.print();
    }
    
    
    
    var<> x_diff("x_diff", pos_), y_diff("y_diff", pos_), z_diff("z_diff", pos_);
    Reg->add(x_diff.in(N1), y_diff.in(N1), z_diff.in(N1));
    
    auto ids1 = theta11.repeat_id(N1.size());
    if(!rigid_transf){
        Constraint<> x_rot1("x_rot1");
        x_rot1 += new_x1;
        x_rot1 -= x1.in(N1)*scale_x*theta11.in(ids1) + y1.in(N1)*scale_x*theta12.in(ids1) + z1.in(N1)*scale_x*theta13.in(ids1);
        Reg->add(x_rot1.in(N1)==0);
        
        Constraint<> y_rot1("y_rot1");
        y_rot1 += new_y1;
        y_rot1 -= x1.in(N1)*scale_y*theta21.in(ids1) + y1.in(N1)*scale_y*theta22.in(ids1) + z1.in(N1)*scale_y*theta23.in(ids1);
        Reg->add(y_rot1.in(N1)==0);
        
        Constraint<> z_rot1("z_rot1");
        z_rot1 += new_z1;
        z_rot1 -= x1.in(N1)*scale_z*theta31.in(ids1) + y1.in(N1)*scale_z*theta32.in(ids1) + z1.in(N1)*theta33.in(ids1);
        Reg->add(z_rot1.in(N1)==0);
    }
    Constraint<> x_norm("x_norm");
    if(rigid_transf)
        x_norm += pow((x1.in(N1)*theta11.in(ids1) + y1.in(N1)*theta12.in(ids1) + z1.in(N1)*theta13.in(ids1) + x_shift.in(ids1)) - new_xm, 2) - x_diff;
    else
        x_norm += pow((new_x1.in(N1) + x_shift.in(ids1)) - new_xm, 2) - x_diff;
    Reg->add(x_norm.in(N1)<=0);
    
    
    Constraint<> y_norm("y_norm");
    if(rigid_transf)
        y_norm += pow((x1.in(N1)*theta21.in(ids1) + y1.in(N1)*theta22.in(ids1) + z1.in(N1)*theta23.in(ids1) + y_shift.in(ids1)) - new_ym, 2) - y_diff;
    else
        y_norm += pow((new_y1.in(N1) + y_shift.in(ids1)) - new_ym, 2) - y_diff;
    Reg->add(y_norm.in(N1)<=0);
    
    Constraint<> z_norm("z_norm");
    if(rigid_transf)
        z_norm += pow((x1.in(N1)*theta31.in(ids1) + y1.in(N1)*theta32.in(ids1) + z1.in(N1)*theta33.in(ids1) + z_shift.in(ids1)) - new_zm, 2) - z_diff;
    else
        z_norm += pow((new_z1.in(N1) + z_shift.in(ids1)) - new_zm, 2) - z_diff;
    Reg->add(z_norm.in(N1)<=0);
    
    bool add_delta_ij = false;
    if(add_delta_ij) {
        var<> delta_ij("delta_ij", 0, 12);
        Reg->add(delta_ij.in(cells));
        
        Constraint<> Norm2("Norm2");
        Norm2 -= delta_ij - pow(new_x1.from(cells) - x2.to(cells),2) - pow(new_y1.from(cells) - y2.to(cells),2) - pow(new_z1.from(cells) - z2.to(cells),2);
        Reg->add(Norm2.in(cells)<=0);
        
        
        Constraint<> DeltaMin("DeltaMin");
        DeltaMin -= delta.from(cells);
        DeltaMin += delta_ij;
        Reg->add_on_off_multivariate_refined(DeltaMin.in(cells)<=0, bin, true);
    }
    
    
    //    Reg->print();
    
    
    if(!relax_sdp){
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
        if(convex){
            Constraint<> row1("row1");
            row1 = pow(theta11,2)+pow(theta12,2)+pow(theta13,2);
            Reg->add(row1.in(range(0,0))<=1);
            Constraint<> row2("row2");
            row2 = pow(theta21,2)+pow(theta22,2)+pow(theta23,2);
            Reg->add(row2.in(range(0,0))<=1);
            Constraint<> row3("row3");
            row3 = pow(theta31,2)+pow(theta32,2)+pow(theta33,2);
            Reg->add(row3.in(range(0,0))<=1);
            Constraint<> col1("col1");
            col1 = pow(theta11,2)+pow(theta21,2)+pow(theta31,2);
            Reg->add(col1.in(range(0,0))<=1);
            Constraint<> col2("col2");
            col2 = pow(theta12,2)+pow(theta22,2)+pow(theta32,2);
            Reg->add(col2.in(range(0,0))<=1);
            Constraint<> col3("col3");
            col3 = pow(theta13,2)+pow(theta23,2)+pow(theta33,2);
            Reg->add(col3.in(range(0,0))<=1);
        }
        else {
            Constraint<> row1("row1");
            row1 = pow(theta11,2)+pow(theta12,2)+pow(theta13,2);
            Reg->add(row1==1);
            Constraint<> row2("row2");
            row2 = pow(theta21,2)+pow(theta22,2)+pow(theta23,2);
            Reg->add(row2==1);
            Constraint<> row3("row3");
            row3 = pow(theta31,2)+pow(theta32,2)+pow(theta33,2);
            Reg->add(row3==1);
            Constraint<> col1("col1");
            col1 = pow(theta11,2)+pow(theta21,2)+pow(theta31,2);
            Reg->add(col1<=1);
            Constraint<> col2("col2");
            col2 = pow(theta12,2)+pow(theta22,2)+pow(theta32,2);
            Reg->add(col2<=1);
            Constraint<> col3("col3");
            col3 = pow(theta13,2)+pow(theta23,2)+pow(theta33,2);
            Reg->add(col3<=1);
        }
        
    }
    
    
    
    //    for (int i = 1; i<=nd; i++) {
    //        string key = to_string(i)+","+to_string(init_matching.at(i-1)+1);
    //        bin._val->at(bin._indices->_keys_map->at(key)) = 1;
    //    }
    
    /* Objective function */
    if(false && mid_point_lb._range->second>0){
        Constraint<> MidPointLB("MidPointLB");
        MidPointLB = x_diff + y_diff + z_diff - mid_point_lb;
        Reg->add(MidPointLB.in(N1) >= 0);
    }
    
    //    Constraint<> min_cost("min_cost");
    //    min_cost=dist_cost*(bin)-(x_diff.from(cells)+y_diff.from(cells)+z_diff.from(cells));
    //    Reg->add(min_cost.in(cells)<=0);
    //    Constraint<> objub("objub");
    //    objub += sum(x_diff + y_diff + z_diff);
    //    Reg->add(objub <= 0.256);
    param<> ones("ones");
    ones.in(N1);
    ones = 1;
    if(perc_outliers==0)
        Reg->min(sum(x_diff + y_diff + z_diff));
    else
        Reg->min(ones.tr()*(keep*x_diff) + ones.tr()*(keep*y_diff) + ones.tr()*(keep*z_diff));
    //    bin._val->at(bin._indices->_keys_map->at("1,1")) = 1;
    //    bin.set_lb("1,1",1);
    //    bin.set_lb("2,2",1);
    //    bin.set_lb("3,3",1);
    //    for (int i = 0; i<nd; i++) {
    //        string key = to_string(i+1)+","+to_string(i+1);
    //        bin._val->at(bin._indices->_keys_map->at(key)) = 1;
    //    }
    //    Reg->print();
    //    solver<> S1(Reg,ipopt);
    //    for(int i = 1; i<= nd; i++){
    //        bin.param<int>::set_val(to_string(i)+","+to_string(i), 1);
    //    }
    //    Reg->print();
    Reg->replace_integers();
    //    if(relax_ints){
    //        solver<> S(Reg,ipopt);
    //        S.run();
    //        Reg->round_solution();
    //    }
    //    else {
    //        solver<> S(Reg,gurobi);
    //        if(!relax_sdp){
    //            S.use_callback();
    //        }
    //        S.run();
    //    }
    //    Reg->print();
    //    Reg->print_int_solution();
    
    //    Reg->print_solution();
    //    DebugOn("Obj before reading solution = " << to_string_with_precision(Reg->get_obj_val(),10) << endl);
    //    Reg->read_solution(Reg->get_name()+"_solution.txt");
    ////    Reg->print_solution();
    //    Reg->reset_constrs();
    //    Reg->is_feasible(1e-6);
    //    Reg->_obj->reeval();
    //    DebugOn("Obj after reading solution = " << to_string_with_precision(Reg->get_obj_val(),10) << endl);
    
    //    {
    //        indices voronoi_ids("voronoi_ids");
    //        voronoi_ids = indices(N1, *norm_x._indices);
    //        auto voronoi_ids_coefs = voronoi_ids.ignore_ith(0, 1);
    ////        auto voronoi_ids_m = voronoi_ids.from_ith(0, 1);
    //        auto voronoi_ids_data = voronoi_ids.ignore_ith(1, 2);
    //        auto voronoi_ids_bin = voronoi_ids.ignore_ith(2, 1);
    ////        Constraint<> Voronoi_model("Voronoi_model");
    ////        Voronoi_model = norm_x.in(voronoi_ids_coefs)*new_xm.in(voronoi_ids_m) + norm_y.in(voronoi_ids_coefs)*new_ym.in(voronoi_ids_m) + norm_z.in(voronoi_ids_coefs)*new_zm.in(voronoi_ids_m) + intercept.in(voronoi_ids_coefs);
    ////        Reg->add_on_off_multivariate_refined(Voronoi_model.in(voronoi_ids)<=0, bin.in(voronoi_ids_bin), true);
    //
    //
    //        auto func = bin.in(voronoi_ids_bin)*(norm_x.in(voronoi_ids_coefs)*new_x1.in(voronoi_ids_data) + norm_y.in(voronoi_ids_coefs)*new_y1.in(voronoi_ids_data) + norm_z.in(voronoi_ids_coefs)*new_z1.in(voronoi_ids_data) + intercept.in(voronoi_ids_coefs));
    //
    //        func.allocate_mem();
    //        func.eval_all();
    //        for (int i = 0; i<func.get_nb_inst(); i++) {
    //            if(func._val->at(i)>1e-6){
    //                DebugOn("instance " <<  i << " is violated \n");
    //                func.print(i,10);
    //                DebugOn(" | violation = " <<  func._val->at(i) << endl);
    //            }
    //        }
    //
    //    }
    //    Debug("row 1 " << pow(theta11.eval(),2)+pow(theta12.eval(),2)+pow(theta13.eval(),2)
    //            << endl);
    //    Debug("row 2 " << pow(theta21.eval(),2)+pow(theta22.eval(),2)+pow(theta23.eval(),2)
    //            << endl);
    //    Debug("row 3 " << pow(theta31.eval(),2)+pow(theta32.eval(),2)+pow(theta33.eval(),2)
    //            << endl);
    //    Debug("col 1 " << pow(theta11.eval(),2)+pow(theta21.eval(),2)+pow(theta31.eval(),2)
    //            << endl);
    //    Debug("col 2 " << pow(theta12.eval(),2)+pow(theta22.eval(),2)+pow(theta32.eval(),2)
    //            << endl);
    //    Debug("col 3 " << pow(theta13.eval(),2)+pow(theta23.eval(),2)+pow(theta33.eval(),2)
    //            << endl);
    //
    //    Debug("row 12 " << (theta11.eval()*theta21.eval())+(theta12.eval()*theta22.eval())+(theta13.eval()*theta23.eval())
    //            << endl);
    //    Debug("row 13 " << (theta11.eval()*theta31.eval())+(theta12.eval()*theta32.eval())+(theta13.eval()*theta33.eval())
    //            << endl);
    //    Debug("row 23 " << (theta21.eval()*theta31.eval())+(theta22.eval()*theta32.eval())+(theta23.eval()*theta33.eval())
    //            << endl);
    //    Debug("Theta matrix = " << endl);
    //    Debug("|" << theta11.eval() << " " << theta12.eval() << " " << theta13.eval() << "|" << endl);
    //    Debug("|" << theta21.eval() << " " << theta22.eval() << " " << theta23.eval() << "|" << endl);
    //    Debug("|" << theta31.eval() << " " << theta32.eval() << " " << theta33.eval() << "|" << endl);
    //    rot_trans[0]=theta11.eval();
    //    rot_trans[1]=theta12.eval();
    //    rot_trans[2]=theta13.eval();;
    //    rot_trans[3]=theta21.eval();
    //    rot_trans[4]=theta22.eval();
    //    rot_trans[5]=theta23.eval();
    //    rot_trans[6]=theta31.eval();
    //    rot_trans[7]=theta32.eval();
    //    rot_trans[8]=theta33.eval();
    //    rot_trans[9]=x_shift.eval();
    //    rot_trans[10]=y_shift.eval();
    //    rot_trans[11]=z_shift.eval();
    //    if(!rigid_transf){
    //        rot_trans[12]=scale_x.eval();
    //        rot_trans[13]=scale_y.eval();
    //        rot_trans[14]=scale_z.eval();
    //        Debug("scale_x = " << scale_x.eval() << endl);
    //        Debug("scale_y = " << scale_y.eval() << endl);
    //        Debug("scale_z = " << scale_z.eval() << endl);
    //    }
    //    auto pitch_val = std::atan2(theta32.eval(), theta33.eval())*180/pi;
    //    auto roll_val = std::atan2(-1*theta31.eval(), std::sqrt(theta32.eval()*theta32.eval()+theta33.eval()*theta33.eval()))*180/pi;
    //    auto yaw_val = std::atan2(theta21.eval(),theta11.eval())*180/pi;
    //    Debug("Roll (degrees) = " << to_string_with_precision(roll_val,12) << endl);
    //    Debug("Pitch (degrees) = " << to_string_with_precision(pitch_val,12) << endl);
    //    Debug("Yaw (degrees) = " << to_string_with_precision(yaw_val,12) << endl);
    //    Debug("x shift = " << x_shift.eval() << endl);
    //    Debug("y shift = " << y_shift.eval() << endl);
    //    Debug("z shift = " << z_shift.eval() << endl);
    
    //    Reg->write_solution();
    //    Reg->write();
    Reg->_status = 0;
    Reg->print();
    return(Reg);
}

shared_ptr<Model<double>> build_norm1_SOC_MIQCP(vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, const indices& valid_cells, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double shift_min_x, double shift_max_x, double shift_min_y, double shift_max_y, double shift_min_z, double shift_max_z, vector<double>& rot_trans, bool convex, const vector<pair<pair<int,int>,pair<int,int>>>& incompatibles, param<>& norm_x, param<>& norm_y, param<>& norm_z,  param<>& intercept, const vector<int>& init_matching, const vector<double>& error_per_point, bool relax_ints){
    double roll_1 = 0, yaw_1 = 0, pitch_1 = 0;
    int nb_pairs = 0, min_nb_pairs = numeric_limits<int>::max(), max_nb_pairs = 0, av_nb_pairs = 0;
    size_t nm = point_cloud_model.size(), nd = point_cloud_data.size();
    
    //    double shift_min_x = -0.25, shift_max_x = 0.25, shift_min_y = -0.25,shift_max_y = 0.25,shift_min_z = -0.25,shift_max_z = 0.25;
    //    double shift_min_x = 0.23, shift_max_x = 0.24, shift_min_y = -0.24,shift_max_y = -0.23,shift_min_z = -0.02,shift_max_z = -0.01;
    //    double shift_min_x = 0.125, shift_max_x = 0.25, shift_min_y = -0.25,shift_max_y = -0.125,shift_min_z = -0.125,shift_max_z = 0;
    //    double yaw_min = -15*pi/180., yaw_max = -10*pi/180., pitch_min = 15*pi/180.,pitch_max = 20.*pi/180.,roll_min = -10*pi/180.,roll_max = -5*pi/180.;
    
    
    vector<double> zeros = {0,0,0};
    
    param<> x1("x1"), x2("x2"), y1("y1"), y2("y2"), z1("z1"), z2("z2");
    int m = av_nb_pairs;
    string i_str, j_str;
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
    
    
    indices Pairs("Pairs"), cells("cells");
    int idx1 = 0;
    int idx2 = 0;
    indices N1("N1"),N2("N2");
    DebugOn("nd = " << nd << endl);
    DebugOn("nm = " << nm << endl);
    
    N1 = range(1,nd);
    N2 = range(1,nm);
    //    cells = valid_cells;
    cells = indices(N1,N2);
    string name="Norm1_MISDP";
    
    auto Reg=make_shared<Model<>>(name);
    
    
    //    Reg->add_param(x1);Reg->add_param(y1);Reg->add_param(z1);
    //    Reg->add_param(x2);Reg->add_param(y2);Reg->add_param(z2);
    var<> x_shift("x_shift", shift_min_x, shift_max_x), y_shift("y_shift", shift_min_y, shift_max_y), z_shift("z_shift", shift_min_z, shift_max_z);
    //    var<> x_shift("x_shift", 0.23, 0.24), y_shift("y_shift", -0.24, -0.23), z_shift("z_shift", -0.02, -0.01);
    
    Reg->add(x_shift.in(R(1)),y_shift.in(R(1)),z_shift.in(R(1)));
    
    var<> scale_x("scale_x", 0.6, 1.4);
    var<> scale_y("scale_y", 0.6, 1.4);
    var<> scale_z("scale_z", 0.6, 1.4);
    Reg->add(scale_x.in(R(1)),scale_y.in(R(1)),scale_z.in(R(1)));
    
    var<int> bin("bin",0,1);
    Reg->add(bin.in(cells));
    DebugOn("Added " << cells.size() << " binary variables" << endl);
    //    double angle_max = 25.*pi/180.;
    
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
    
    //            var<> theta11("theta11",  0.8, 1), theta12("theta12", -1, 1), theta13("theta13", -1, 1);
    //            var<> theta21("theta21",  -1, 1), theta22("theta22", -1, 1), theta23("theta23", -1, 1);
    //            var<> theta31("theta31",  -1, 1), theta32("theta32", -1, 1), theta33("theta33", 0.8, 1);
    
    //            var<> theta11("theta11",  0.96, 0.97), theta12("theta12", 0.15, 0.16), theta13("theta13", -0.22, -0.2);
    //            var<> theta21("theta21",  -0.21, -0.2), theta22("theta22", 0.95, 0.96), theta23("theta23", -0.24, 0.23);
    //            var<> theta31("theta31",  0.17, 0.18), theta32("theta32", 0.26, 0.27), theta33("theta33", 0.94, 0.95);
    Reg->add(theta11.in(R(1)),theta12.in(R(1)),theta13.in(R(1)));
    Reg->add(theta21.in(R(1)),theta22.in(R(1)),theta23.in(R(1)));
    Reg->add(theta31.in(R(1)),theta32.in(R(1)),theta33.in(R(1)));
    //    Reg->print();
    
    param<> x_new_lb("x_new_lb");
    x_new_lb.in(N1);
    param<> x_new_ub("x_new_ub");
    x_new_ub.in(N1);
    param<> y_new_lb("y_new_lb");
    y_new_lb.in(N1);
    param<> y_new_ub("y_new_ub");
    y_new_ub.in(N1);
    param<> z_new_lb("z_new_lb");
    z_new_lb.in(N1);
    param<> z_new_ub("z_new_ub");
    z_new_ub.in(N1);
    
    param<> mid_point_lb("mid_point_lb");
    mid_point_lb.in(N1);
    
    double lower_bound = 0, go_icp_lb = 0, SDP_lb = 0;
    double x_lb = 0, y_lb = 0, z_lb = 0, x1_i = 0, y1_i = 0, z1_i = 0;
    shared_ptr<pair<double,double>> x1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> y1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> z1_bounds = make_shared<pair<double,double>>();
    for (int i = 0; i<nd; i++) {
        x1_bounds->first = x1.eval(i);
        x1_bounds->second = x1.eval(i);
        y1_bounds->first = y1.eval(i);
        y1_bounds->second = y1.eval(i);
        z1_bounds->first = z1.eval(i);
        z1_bounds->second = z1.eval(i);
        auto x_range  = get_product_range(x1_bounds, theta11._range);
        auto y_range  = get_product_range(y1_bounds, theta12._range);
        auto z_range  = get_product_range(z1_bounds, theta13._range);
        auto bounds = get_min_max(roll_min, roll_max, pitch_min, pitch_max, yaw_min, yaw_max, point_cloud_data[i], zeros);
        
        //        double new_roll_max1 = acos(sqrt(cos(sqrt(3)*roll_max)*(2./3.) + 1./3.));
        //        double new_roll_max = asin(sin(sqrt(3.)*roll_max)/sqrt(3.) + 1./3.*(1 - cos(sqrt(3.)*roll_max)));
        
        //        auto go_icp_max_dist = get_GoICP_dist(roll_max, shift_max_x, point_cloud_data[i], false);
        
        //        double v1 = 1./std::sqrt(3.), v2 = 1./std::sqrt(3.), v3 = 1./std::sqrt(3.), ct = std::cos(std::sqrt(3.)*roll_max), ct2 = (1.-std::cos(std::sqrt(3.)*roll_max)), st = std::sin(std::sqrt(3.)*roll_max);
        //        double tmp231 = v2*v3*ct2, tmp232 = v1*st;
        //        auto tmp121 = v2*v2*ct2;
        //        auto tmp122 = v3*st;
        //        auto tmp131 = v1*v3*ct2;
        //        auto tmp132 = v2*st;
        //        auto R11 = ct + v1*v1*ct2;
        //        auto R12 = tmp121 - tmp122;
        //        auto R13 = tmp131 + tmp132;
        //        auto R21 = tmp121 + tmp122;
        //        auto R22 = ct + v2*v2*ct2;
        //        auto R23 = tmp231 - tmp232;
        //        double R31 = tmp131 - tmp132;
        //        double R32 = tmp231 + tmp232;
        //        double R33 = ct + v3*v3*ct2;
        //        vector<double> rot_mat = {R11, R12, R13, R21, R22, R23, R31, R32, R33};
        //
        ////        auto SDP = build_SDP(point_cloud_data[i], rot_mat);
        ////        SDP->is_feasible(1e-6);
        ////        double max_sdp = std::sqrt(SDP->get_obj_val());
        //
        //        double new_pitch_max = std::atan2(R32, R33);
        //        auto new_yaw_max = std::atan2(R21,R11);
        //        auto new_roll_max = std::atan2(-1*R31, std::sqrt(R32*R32+R33*R33));
        
        auto go_icp_max_dist = get_GoICP_dist(roll_max, shift_max_x, point_cloud_data[i], false);
        //        DebugOn("SDP max dist = " << to_string_with_precision(max_sdp, 6) << endl);
        //        auto max_dist = get_max_dist(-new_roll_max, new_roll_max, -new_pitch_max, new_pitch_max, -new_yaw_max, new_yaw_max, shift_min_x, shift_max_x, shift_min_y, shift_max_y, shift_min_z, shift_max_z, point_cloud_data[i], zeros, false);
        
        //        lower_bound += std::max(0.,error_per_point[i] - max_dist);
        go_icp_lb += std::max(0.,error_per_point[i] - go_icp_max_dist);
        //        SDP_lb += std::max(0.,error_per_point[i] - max_sdp);
        //        mid_point_lb.set_val(i,std::max(0.,error_per_point[i] - max_dist));
        auto xlb = x_range->first + y_range->first + z_range->first + x_shift.get_lb().eval();
        auto xub = x_range->second + y_range->second + z_range->second+ x_shift.get_ub().eval();
        x_new_lb.set_val(i,bounds[0].first + x_shift.get_lb().eval());
        x_new_ub.set_val(i,bounds[0].second + x_shift.get_ub().eval());
        //        x_new_lb.set_val(i, x_range->first + y_range->first + z_range->first + x_shift.get_lb().eval());
        //        x_new_ub.set_val(i, x_range->second + y_range->second + z_range->second+ x_shift.get_ub().eval());
        x_range  = get_product_range(x1_bounds, theta21._range);
        y_range  = get_product_range(y1_bounds, theta22._range);
        z_range  = get_product_range(z1_bounds, theta23._range);
        auto ylb = x_range->first + y_range->first + z_range->first + y_shift.get_lb().eval();
        auto yub = x_range->second + y_range->second + z_range->second + y_shift.get_ub().eval();
        y_new_lb.set_val(i,bounds[1].first + y_shift.get_lb().eval());
        y_new_ub.set_val(i,bounds[1].second + y_shift.get_ub().eval());
        //        y_new_lb.set_val(i, x_range->first + y_range->first + z_range->first + y_shift.get_lb().eval());
        //        y_new_ub.set_val(i, x_range->second + y_range->second + z_range->second+ y_shift.get_ub().eval());
        x_range  = get_product_range(x1_bounds, theta31._range);
        y_range  = get_product_range(y1_bounds, theta32._range);
        z_range  = get_product_range(z1_bounds, theta33._range);
        auto zlb = x_range->first + y_range->first + z_range->first + z_shift.get_lb().eval();
        auto zub = x_range->second + y_range->second + z_range->second + z_shift.get_ub().eval();
        z_new_lb.set_val(i,bounds[2].first + z_shift.get_lb().eval());
        z_new_ub.set_val(i,bounds[2].second + z_shift.get_ub().eval());
        //        z_new_lb.set_val(i, x_range->first + y_range->first + z_range->first + z_shift.get_lb().eval());
        //        z_new_ub.set_val(i, x_range->second + y_range->second + z_range->second+ z_shift.get_ub().eval());
    }
    //    DebugOn("Lower bound = " << to_string_with_precision(lower_bound,6) << endl);
    DebugOn("GoICP lower bound = " << to_string_with_precision(go_icp_lb,6) << endl);
    //    DebugOn("SDP lower bound = " << to_string_with_precision(SDP_lb,6) << endl);
    
    var<> new_xm("new_xm", -1, 1), new_ym("new_ym", -1, 1), new_zm("new_zm", -1, 1);
    var<> new_x1("new_x1", x_new_lb, x_new_ub), new_y1("new_y1", y_new_lb, y_new_ub), new_z1("new_z1", z_new_lb, z_new_ub);
    //            var<> new_x1("new_x1", -1, 1), new_y1("new_y1", -1, 1), new_z1("new_z1", -1, 1);
    Reg->add(new_xm.in(N1), new_ym.in(N1), new_zm.in(N1));
    //    Reg->add(new_x1.in(N1), new_y1.in(N1), new_z1.in(N1));
    
    //    var<> new_x1_ij("new_x1_ij", -1, 1), new_y1_ij("new_y1_ij", -1, 1), new_z1_ij("new_z1_ij", -1, 1);
    //    Reg->add(new_x1_ij.in(cells), new_y1_ij.in(cells), new_z1_ij.in(cells));
    
    
    var<> delta("delta", 0,12);
    //    Reg->add(delta.in(N1));
    
    
    //    var<> d1("d1", 0,4),d2("d2", 0,4),d3("d3", 0,4),d4("d4", 0,4);
    //    var<> l12("l12", -2,2),l13("l13", -2,2),l14("l14", -2,2),l23("l23", -2,2),l24("l24", -2,2),l34("l34", -2,2);
    //    Reg->add(l12.in(R(1)),l13.in(R(1)),l14.in(R(1)));
    //    Reg->add(l23.in(R(1)),l24.in(R(1)),l34.in(R(1)));
    //    Reg->add(d1.in(R(1)),d2.in(R(1)),d3.in(R(1)),d4.in(R(1)));
    
    
    indices ids = indices("in_x");
    ids.add_empty_row();
    for(auto i=0;i<nd;i++){
        for(auto j=1;j<=nm;j++){
            if(cells.has_key(to_string(i+1)+","+to_string(j)))
                ids.add_in_row(i, to_string(j));
        }
    }
    
    bool add_shift_cut = false;
    
    if(add_shift_cut){
        Constraint<> Centroid("Centroid");
        Centroid = sum(new_x1) + sum(new_y1) + sum(new_z1) - (sum(new_xm) + sum(new_ym) + sum(new_zm));
        Reg->add(Centroid==0);
        
        Constraint<> Centroid_tx("Centroid_tx");
        Centroid_tx = nd*x_shift - sum(new_xm);
        Reg->add(Centroid_tx==0);
        
        Constraint<> Centroid_ty("Centroid_ty");
        Centroid_ty = nd*y_shift - sum(new_ym);
        Reg->add(Centroid_ty==0);
        
        Constraint<> Centroid_tz("Centroid_tz");
        Centroid_tz = nd*z_shift - sum(new_zm);
        Reg->add(Centroid_tz==0);
    }
    
    bool add_voronoi = false;
    
    if(add_voronoi){
        indices voronoi_ids("voronoi_ids");
        voronoi_ids = indices(range(1,3), *norm_x._indices);
        auto voronoi_ids_coefs = voronoi_ids.ignore_ith(0, 1);
        auto voronoi_ids_m = voronoi_ids.from_ith(0, 1);
        auto voronoi_ids_data = voronoi_ids.ignore_ith(1, 2);
        auto voronoi_ids_bin = voronoi_ids.ignore_ith(2, 1);
        //        Constraint<> Voronoi_model("Voronoi_model");
        //        Voronoi_model = norm_x.in(voronoi_ids_coefs)*new_xm.in(voronoi_ids_m) + norm_y.in(voronoi_ids_coefs)*new_ym.in(voronoi_ids_m) + norm_z.in(voronoi_ids_coefs)*new_zm.in(voronoi_ids_m) + intercept.in(voronoi_ids_coefs);
        //        Reg->add_on_off_multivariate_refined(Voronoi_model.in(voronoi_ids)<=0, bin.in(voronoi_ids_bin), true);
        
        
        Constraint<> Voronoi("Voronoi");
        Voronoi = norm_x.in(voronoi_ids_coefs)*new_x1.in(voronoi_ids_data) + norm_y.in(voronoi_ids_coefs)*new_y1.in(voronoi_ids_data) + norm_z.in(voronoi_ids_coefs)*new_z1.in(voronoi_ids_data) + intercept.in(voronoi_ids_coefs);
        Reg->add_on_off_multivariate_refined(Voronoi.in(voronoi_ids)<=0, bin.in(voronoi_ids_bin), true);
        //        Reg->print();
    }
    
    
    //    Reg->print();
    
    
    //
    theta11.initialize_all(1);
    theta22.initialize_all(1);
    theta33.initialize_all(1);
    
    bool spatial_branching = false;
    if(spatial_branching){
        /* Spatial branching vars */
        int nb_pieces = 5; // Divide each axis into nb_pieces
        indices spatial_ids("spatial_ids");
        spatial_ids = range(1,nb_pieces);
        indices theta_ids("theta_ids");
        theta_ids = indices(range(1,3),range(1,3));
        indices shift_ids("shift_ids");
        shift_ids.insert({"x", "y", "z"});
        indices theta_spatial("theta_spatial");
        theta_spatial = indices(theta_ids, spatial_ids);
        
        indices shift_spatial("shift_spatial");
        shift_spatial = indices(shift_ids, spatial_ids);
        
        //    var<int> sbin_shift("sbin_shift", 0, 1);
        //    var<int> sbin_theta("sbin_theta", 0, 1);
        
        //    Reg->add(sbin_theta.in(theta_spatial),sbin_shift.in(shift_spatial));
        
        var<int> sbin_tx("sbin_tx", 0, 1), sbin_ty("sbin_ty", 0, 1), sbin_tz("sbin_tz", 0, 1);
        var<int> sbin_theta11("sbin_theta11", 0, 1), sbin_theta12("sbin_theta12", 0, 1), sbin_theta13("sbin_theta13", 0, 1);
        var<int> sbin_theta21("sbin_theta21", 0, 1), sbin_theta22("sbin_theta22", 0, 1), sbin_theta23("sbin_theta23", 0, 1);
        var<int> sbin_theta31("sbin_theta31", 0, 1), sbin_theta32("sbin_theta32", 0, 1), sbin_theta33("sbin_theta33", 0, 1);
        Reg->add(sbin_theta11.in(spatial_ids),sbin_theta12.in(spatial_ids),sbin_theta13.in(spatial_ids));
        Reg->add(sbin_theta21.in(spatial_ids),sbin_theta22.in(spatial_ids),sbin_theta23.in(spatial_ids));
        Reg->add(sbin_theta31.in(spatial_ids),sbin_theta32.in(spatial_ids),sbin_theta33.in(spatial_ids));
        Reg->add(sbin_tx.in(spatial_ids),sbin_ty.in(spatial_ids),sbin_tz.in(spatial_ids));
        /* Spatial branching constraints */
        Constraint<> OneBinAngleSpatial11("OneBinAngleSpatial11");
        OneBinAngleSpatial11 = sum(sbin_theta11);
        Reg->add(OneBinAngleSpatial11==1);
        
        Constraint<> OneBinAngleSpatial12("OneBinAngleSpatial12");
        OneBinAngleSpatial12 = sum(sbin_theta12);
        Reg->add(OneBinAngleSpatial12==1);
        
        Constraint<> OneBinAngleSpatial13("OneBinAngleSpatial13");
        OneBinAngleSpatial13 = sum(sbin_theta13);
        Reg->add(OneBinAngleSpatial13==1);
        
        Constraint<> OneBinAngleSpatial21("OneBinAngleSpatial21");
        OneBinAngleSpatial21 = sum(sbin_theta21);
        Reg->add(OneBinAngleSpatial21==1);
        
        Constraint<> OneBinAngleSpatial22("OneBinAngleSpatial22");
        OneBinAngleSpatial22 = sum(sbin_theta22);
        Reg->add(OneBinAngleSpatial22==1);
        
        
        Constraint<> OneBinAngleSpatial23("OneBinAngleSpatial23");
        OneBinAngleSpatial23 = sum(sbin_theta23);
        Reg->add(OneBinAngleSpatial23==1);
        
        Constraint<> OneBinAngleSpatial31("OneBinAngleSpatial31");
        OneBinAngleSpatial31 = sum(sbin_theta31);
        Reg->add(OneBinAngleSpatial31==1);
        
        Constraint<> OneBinAngleSpatial32("OneBinAngleSpatial32");
        OneBinAngleSpatial32 = sum(sbin_theta32);
        Reg->add(OneBinAngleSpatial32==1);
        
        Constraint<> OneBinAngleSpatial33("OneBinAngleSpatial33");
        OneBinAngleSpatial33 = sum(sbin_theta33);
        Reg->add(OneBinAngleSpatial33==1);
        
        Constraint<> OneBinShiftSpatialx("OneBinShiftSpatialx");
        OneBinShiftSpatialx = sum(sbin_tx);
        Reg->add(OneBinShiftSpatialx==1);
        
        Constraint<> OneBinShiftSpatialy("OneBinShiftSpatialy");
        OneBinShiftSpatialy = sum(sbin_ty);
        Reg->add(OneBinShiftSpatialy==1);
        
        Constraint<> OneBinShiftSpatialz("OneBinShiftSpatialz");
        OneBinShiftSpatialz = sum(sbin_tz);
        Reg->add(OneBinShiftSpatialz==1);
        
        //    Reg->print();
        
        double diag_increment = 1./nb_pieces;/* Diagonals are defined in [0,1] */
        double off_diag_increment = 2./nb_pieces;/* Diagonals are defined in [-1,1] */
        double shift_increment = 0.5/nb_pieces;/* Shifts are defined in [-0.25,0.25] */
        
        auto spatial_ids_n = range(1,nb_pieces-1);
        auto spatial_ids_1 = range(2,nb_pieces);
        param<> diag_lb("diag_lb"), diag_ub("diag_ub"), off_diag_lb("off_diag_lb"), off_diag_ub("off_diag_ub");
        param<> t_lb("t_lb"), t_ub("t_ub");
        diag_ub.in(spatial_ids);
        diag_lb.in(spatial_ids);
        off_diag_ub.in(spatial_ids);
        off_diag_lb.in(spatial_ids);
        t_ub.in(spatial_ids);
        t_lb.in(spatial_ids);
        for (int i = 0; i<nb_pieces; i++) {
            diag_ub.set_val(i,(i+1)*diag_increment);
            off_diag_ub.set_val(i,-1+(i+1)*off_diag_increment);
            t_ub.set_val(i,-0.25 + (i+1)*diag_increment);
            diag_lb.set_val(i,i*diag_increment);
            off_diag_lb.set_val(i,-1 + i*off_diag_increment);
            t_lb.set_val(i,-0.25 + i*diag_increment);
        }
        auto ids_repeat = theta11.repeat_id(nb_pieces-1);
        
        Constraint<> Diag_Spatial_UB11("Diag_Spatial_UB11");
        Diag_Spatial_UB11 = theta11.in(ids_repeat) - diag_ub.in(spatial_ids_n);
        Reg->add_on_off_multivariate_refined(Diag_Spatial_UB11.in(spatial_ids_n)<=0, sbin_theta11.in(spatial_ids_n), true);
        
        Constraint<> Diag_Spatial_LB11("Diag_Spatial_LB11");
        Diag_Spatial_LB11 = diag_lb.in(spatial_ids_1) - theta11.in(ids_repeat);
        Reg->add_on_off_multivariate_refined(Diag_Spatial_LB11.in(spatial_ids_1) <= 0, sbin_theta11.in(spatial_ids_1), true);
        
        Constraint<> Diag_Spatial_UB22("Diag_Spatial_UB22");
        Diag_Spatial_UB22 = theta22.in(ids_repeat) - diag_ub.in(spatial_ids_n);
        Reg->add_on_off_multivariate_refined(Diag_Spatial_UB22.in(spatial_ids_n)<=0, sbin_theta22.in(spatial_ids_n), true);
        
        Constraint<> Diag_Spatial_LB33("Diag_Spatial_LB33");
        Diag_Spatial_LB33 = diag_lb.in(spatial_ids_1) - theta33.in(ids_repeat);
        Reg->add_on_off_multivariate_refined(Diag_Spatial_LB33.in(spatial_ids_1) <= 0, sbin_theta33.in(spatial_ids_1), true);
        
        
        Constraint<> Diag_Spatial_UB12("Diag_Spatial_UB12");
        Diag_Spatial_UB12 = theta12.in(ids_repeat) - off_diag_ub.in(spatial_ids_n);
        Reg->add_on_off_multivariate_refined(Diag_Spatial_UB12.in(spatial_ids_n) <= 0, sbin_theta12.in(spatial_ids_n), true);
        
        Constraint<> Diag_Spatial_LB12("Diag_Spatial_LB12");
        Diag_Spatial_LB12 = off_diag_lb.in(spatial_ids_1) - theta12.in(ids_repeat);
        Reg->add_on_off_multivariate_refined(Diag_Spatial_LB12.in(spatial_ids_1) <= 0, sbin_theta12.in(spatial_ids_1), true);
        
        Constraint<> Diag_Spatial_UB13("Diag_Spatial_UB13");
        Diag_Spatial_UB13 = theta13.in(ids_repeat) - off_diag_ub.in(spatial_ids_n);
        Reg->add_on_off_multivariate_refined(Diag_Spatial_UB13.in(spatial_ids_n) <= 0, sbin_theta13.in(spatial_ids_n), true);
        
        Constraint<> Diag_Spatial_LB13("Diag_Spatial_LB13");
        Diag_Spatial_LB13 = off_diag_lb.in(spatial_ids_1) - theta13.in(ids_repeat);
        Reg->add_on_off_multivariate_refined(Diag_Spatial_LB13.in(spatial_ids_1) <= 0, sbin_theta13.in(spatial_ids_1), true);
        
        Constraint<> Diag_Spatial_UB21("Diag_Spatial_UB21");
        Diag_Spatial_UB21 = theta21.in(ids_repeat) - off_diag_ub.in(spatial_ids_n);
        Reg->add_on_off_multivariate_refined(Diag_Spatial_UB21.in(spatial_ids_n) <= 0, sbin_theta21.in(spatial_ids_n), true);
        
        Constraint<> Diag_Spatial_LB21("Diag_Spatial_LB21");
        Diag_Spatial_LB21 = off_diag_lb.in(spatial_ids_1) - theta21.in(ids_repeat);
        Reg->add_on_off_multivariate_refined(Diag_Spatial_LB21.in(spatial_ids_1) <= 0, sbin_theta21.in(spatial_ids_1), true);
        
        Constraint<> Diag_Spatial_UB23("Diag_Spatial_UB23");
        Diag_Spatial_UB23 = theta23.in(ids_repeat) - off_diag_ub.in(spatial_ids_n);
        Reg->add_on_off_multivariate_refined(Diag_Spatial_UB23.in(spatial_ids_n) <= 0, sbin_theta23.in(spatial_ids_n), true);
        
        Constraint<> Diag_Spatial_LB23("Diag_Spatial_LB23");
        Diag_Spatial_LB23 = off_diag_lb.in(spatial_ids_1) - theta23.in(ids_repeat);
        Reg->add_on_off_multivariate_refined(Diag_Spatial_LB23.in(spatial_ids_1) <= 0, sbin_theta23.in(spatial_ids_1), true);
        
        Constraint<> Diag_Spatial_UB31("Diag_Spatial_UB31");
        Diag_Spatial_UB31 = theta31.in(ids_repeat) - off_diag_ub.in(spatial_ids_n);
        Reg->add_on_off_multivariate_refined(Diag_Spatial_UB31.in(spatial_ids_n) <= 0, sbin_theta31.in(spatial_ids_n), true);
        
        Constraint<> Diag_Spatial_LB31("Diag_Spatial_LB31");
        Diag_Spatial_LB31 = off_diag_lb.in(spatial_ids_1) - theta31.in(ids_repeat);
        Reg->add_on_off_multivariate_refined(Diag_Spatial_LB31.in(spatial_ids_1) <= 0, sbin_theta31.in(spatial_ids_1), true);
        
        Constraint<> Diag_Spatial_UB32("Diag_Spatial_UB32");
        Diag_Spatial_UB32 = theta32.in(ids_repeat) - off_diag_ub.in(spatial_ids_n);
        Reg->add_on_off_multivariate_refined(Diag_Spatial_UB32.in(spatial_ids_n) <= 0, sbin_theta32.in(spatial_ids_n), true);
        
        Constraint<> Diag_Spatial_LB32("Diag_Spatial_LB32");
        Diag_Spatial_LB32 = off_diag_lb.in(spatial_ids_1) - theta32.in(ids_repeat);
        Reg->add_on_off_multivariate_refined(Diag_Spatial_LB32.in(spatial_ids_1) <= 0, sbin_theta32.in(spatial_ids_1), true);
        
        Constraint<> xshift_Spatial_UB("xshift_Spatial_UB");
        xshift_Spatial_UB = x_shift.in(ids_repeat) - t_ub.in(spatial_ids_n);
        Reg->add_on_off_multivariate_refined(xshift_Spatial_UB.in(spatial_ids_n) <= 0, sbin_tx.in(spatial_ids_n), true);
        
        Constraint<> xshift_Spatial_LB("xshift_Spatial_LB");
        xshift_Spatial_LB = t_lb.in(spatial_ids_1) - x_shift.in(ids_repeat);
        Reg->add_on_off_multivariate_refined(xshift_Spatial_LB.in(spatial_ids_1) <= 0, sbin_tx.in(spatial_ids_1), true);
        
        Constraint<> yshift_Spatial_UB("yshift_Spatial_UB");
        yshift_Spatial_UB = y_shift.in(ids_repeat) - t_ub.in(spatial_ids_n);
        Reg->add_on_off_multivariate_refined(yshift_Spatial_UB.in(spatial_ids_n) <= 0, sbin_ty.in(spatial_ids_n), true);
        
        Constraint<> yshift_Spatial_LB("yshift_Spatial_LB");
        yshift_Spatial_LB = t_lb.in(spatial_ids_1) - y_shift.in(ids_repeat);
        Reg->add_on_off_multivariate_refined(yshift_Spatial_LB.in(spatial_ids_1) <= 0, sbin_ty.in(spatial_ids_1), true);
        
        Constraint<> zshift_Spatial_UB("zshift_Spatial_UB");
        zshift_Spatial_UB = z_shift.in(ids_repeat) - t_ub.in(spatial_ids_n);
        Reg->add_on_off_multivariate_refined(zshift_Spatial_UB.in(spatial_ids_n) <= 0, sbin_tz.in(spatial_ids_n), true);
        
        Constraint<> zshift_Spatial_LB("zshift_Spatial_LB");
        zshift_Spatial_LB = t_lb.in(spatial_ids_1) - z_shift.in(ids_repeat);
        Reg->add_on_off_multivariate_refined(zshift_Spatial_LB.in(spatial_ids_1) <= 0, sbin_tz.in(spatial_ids_1), true);
        //        Reg->print();
    }
    
    
    Constraint<> Def_newxm("Def_newxm");
    Def_newxm = new_xm-product(x2.in(ids),bin.in_matrix(1, 1));
    Reg->add(Def_newxm.in(N1)==0);
    
    Constraint<> Def_newym("Def_newym");
    Def_newym = new_ym-product(y2.in(ids),bin.in_matrix(1, 1));
    Reg->add(Def_newym.in(N1)==0);
    
    Constraint<> Def_newzm("Def_newzm");
    Def_newzm = new_zm-product(z2.in(ids),bin.in_matrix(1, 1));
    Reg->add(Def_newzm.in(N1)==0);
    
    Constraint<> OneBin("OneBin");
    OneBin = bin.in_matrix(1, 1);
    Reg->add(OneBin.in(N1)==1);
    
    //    Constraint<> OneBin2("OneBin2");
    //    OneBin2 = bin.in_matrix(0, 1);
    //    Reg->add(OneBin2.in(N2)<=1);
    
    if(false && !incompatibles.empty()){
        indices pairs1("pairs1"), pairs2("pairs2");
        pairs1 = cells;
        pairs2 = cells;
        for (const auto &inc_pair : incompatibles) {
            pairs1.add_ref(to_string(inc_pair.first.first+1)+","+to_string(inc_pair.second.first+1));
            pairs2.add_ref(to_string(inc_pair.first.second+1)+","+to_string(inc_pair.second.second+1));
        }
        
        Constraint<> incomp_pairs("incomp_pairs");
        incomp_pairs = bin.in(pairs1) + bin.in(pairs2);
        Reg->add(incomp_pairs.in(range(1,pairs1.size()))<=1);
        //        incomp_pairs.print();
    }
    
    
    
    var<> x_diff("x_diff", pos_), y_diff("y_diff", pos_), z_diff("z_diff", pos_);
    Reg->add(x_diff.in(N1), y_diff.in(N1), z_diff.in(N1));
    
    auto ids1 = theta11.repeat_id(N1.size());
    //    Constraint<> x_rot1("x_rot1");
    //    x_rot1 += new_x1 -x_shift;
    //    x_rot1 -= x1.in(N1)*theta11.in(ids1) + y1.in(N1)*theta12.in(ids1) + z1.in(N1)*theta13.in(ids1);
    //    Reg->add(x_rot1.in(N1)==0);
    //
    //    Constraint<> y_rot1("y_rot1");
    //    y_rot1 += new_y1 - y_shift;
    //    y_rot1 -= x1.in(N1)*theta21.in(ids1) + y1.in(N1)*theta22.in(ids1) + z1.in(N1)*theta23.in(ids1);
    //    Reg->add(y_rot1.in(N1)==0);
    //
    //    Constraint<> z_rot1("z_rot1");
    //    z_rot1 += new_z1 -z_shift;
    //    z_rot1 -= x1.in(N1)*theta31.in(ids1) + y1.in(N1)*theta32.in(ids1) + z1.in(N1)*theta33.in(ids1);
    //    Reg->add(z_rot1.in(N1)==0);
    Constraint<> x_abs1("x_abs1");
    x_abs1 += x_diff - ((x1.in(N1)*scale_x.in(ids1)*theta11.in(ids1) + y1.in(N1)*scale_x.in(ids1)*theta12.in(ids1) + z1.in(N1)*scale_x.in(ids1)*theta13.in(ids1) + x_shift.in(ids1)) - new_xm);
    Reg->add(x_abs1.in(N1)>=0);
    
    Constraint<> x_abs2("x_abs2");
    x_abs2 += x_diff - (new_xm - (x1.in(N1)*scale_x.in(ids1)*theta11.in(ids1) + y1.in(N1)*scale_x.in(ids1)*theta12.in(ids1) + z1.in(N1)*scale_x.in(ids1)*theta13.in(ids1) + x_shift.in(ids1)));
    Reg->add(x_abs2.in(N1)>=0);
    
    Constraint<> y_abs1("y_abs1");
    y_abs1 += y_diff - ((x1.in(N1)*scale_y.in(ids1)*theta21.in(ids1) + y1.in(N1)*scale_y.in(ids1)*theta22.in(ids1) + z1.in(N1)*scale_y.in(ids1)*theta23.in(ids1) + y_shift.in(ids1)) - new_ym);
    Reg->add(y_abs1.in(N1)>=0);
    
    Constraint<> y_abs2("y_abs2");
    y_abs2 += y_diff - (new_ym - (x1.in(N1)*scale_y.in(ids1)*theta21.in(ids1) + y1.in(N1)*scale_y.in(ids1)*theta22.in(ids1) + z1.in(N1)*scale_y.in(ids1)*theta23.in(ids1) + y_shift.in(ids1)));
    Reg->add(y_abs2.in(N1)>=0);
    
    Constraint<> z_abs1("z_abs1");
    z_abs1 += z_diff - ((x1.in(N1)*scale_z.in(ids1)*theta31.in(ids1) + y1.in(N1)*scale_z.in(ids1)*theta32.in(ids1) + z1.in(N1)*scale_z.in(ids1)*theta33.in(ids1) + z_shift.in(ids1)) - new_zm);
    Reg->add(z_abs1.in(N1)>=0);
    
    Constraint<> z_abs2("z_abs2");
    z_abs2 += z_diff - (new_zm - (x1.in(N1)*scale_z.in(ids1)*theta31.in(ids1) + y1.in(N1)*scale_z.in(ids1)*theta32.in(ids1) + z1.in(N1)*scale_z.in(ids1)*theta33.in(ids1) + z_shift.in(ids1)));
    Reg->add(z_abs2.in(N1)>=0);
    
    bool add_delta_ij = false;
    if(add_delta_ij) {
        //        Reg->add(delta.in(N1));
        //
        //        Constraint<> Norm1("Norm1");
        //        Norm1 += delta - (x_diff + y_diff + z_diff);
        //        Reg->add(Norm1.in(N1)>=0);
        
        var<> x_diff_ij("x_diff_ij", pos_), y_diff_ij("y_diff_ij", pos_), z_diff_ij("z_diff_ij", pos_);
        Reg->add(x_diff_ij.in(cells), y_diff_ij.in(cells), z_diff_ij.in(cells));
        
        //        var<> delta_ij("delta_ij", 0, 6);
        //        Reg->add(delta_ij.in(cells));
        auto ids1 = theta11.repeat_id(cells.size());
        Constraint<> x_abs1_new("x_abs1_new");
        x_abs1_new += x_diff_ij - ((x1.from(cells)*theta11.in(ids1) + y1.from(cells)*theta12.in(ids1) + z1.from(cells)*theta13.in(ids1) + x_shift.in(ids1)) - x2.to(cells));
        Reg->add(x_abs1_new.in(cells)>=0);
        
        Constraint<> x_abs2_new("x_abs2_new");
        x_abs2_new += x_diff_ij - (x2.to(cells) - (x1.from(cells)*theta11.in(ids1) + y1.from(cells)*theta12.in(ids1) + z1.from(cells)*theta13.in(ids1) + x_shift.in(ids1)));
        Reg->add(x_abs2_new.in(cells)>=0);
        
        Constraint<> y_abs1_new("y_abs1_new");
        y_abs1_new += y_diff_ij - (y2.to(cells) - (x1.from(cells)*theta21.in(ids1) + y1.from(cells)*theta22.in(ids1) + z1.from(cells)*theta23.in(ids1) + y_shift.in(ids1)));
        Reg->add(y_abs1_new.in(cells)>=0);
        
        Constraint<> y_abs2_new("y_abs2_new");
        y_abs2_new += y_diff_ij - ((x1.from(cells)*theta21.in(ids1) + y1.from(cells)*theta22.in(ids1) + z1.from(cells)*theta23.in(ids1) + y_shift.in(ids1)) - y2.to(cells));
        Reg->add(y_abs2_new.in(cells)>=0);
        
        Constraint<> z_abs1_new("z_abs1_new");
        z_abs1_new += z_diff_ij - ((x1.from(cells)*theta31.in(ids1) + y1.from(cells)*theta32.in(ids1) + z1.from(cells)*theta33.in(ids1) + z_shift.in(ids1)) - z2.to(cells));
        Reg->add(z_abs1_new.in(cells)>=0);
        
        Constraint<> z_abs2_new("z_abs2_new");
        z_abs2_new += z_diff_ij - (z2.to(cells) - (x1.from(cells)*theta31.in(ids1) + y1.from(cells)*theta32.in(ids1) + z1.from(cells)*theta33.in(ids1) + z_shift.in(ids1)));
        Reg->add(z_abs2_new.in(cells)>=0);
        
        
        Constraint<> x_min("x_min");
        x_min += x_diff.from(cells) - x_diff_ij;
        Reg->add(x_min.in(cells)<=0);
        
        Constraint<> y_min("y_min");
        y_min += y_diff.from(cells) - y_diff_ij;
        Reg->add(y_min.in(cells)<=0);
        
        Constraint<> z_min("z_min");
        z_min += z_diff.from(cells) - z_diff_ij;
        Reg->add(z_min.in(cells)<=0);
        
        
        
        
        //        Constraint<> Norm1_ij("Norm1_ij");
        //        Norm1_ij += delta_ij - (x_diff_ij + y_diff_ij + z_diff_ij);
        //        Reg->add(Norm1_ij.in(cells)>=0);
        
        //        Constraint<> DeltaMin("DeltaMin");
        //        DeltaMin -= delta.from(cells);
        //        DeltaMin += delta_ij;
        //        Reg->add_on_off_multivariate_refined(DeltaMin.in(cells)<=0, bin, true);
    }
    
    
    //    Reg->print();
    
    bool add_sdp_rel = true;
    if(add_sdp_rel){
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
        //        soc_12.add_to_callback();
        Reg->add(soc_12.in(range(0,0))<=0);
        
        Constraint<> soc_13("soc_13");
        soc_13 = pow(theta12-theta21,2)-(1-theta11-theta22+theta33)*(1+theta11+theta22+theta33);
        //        soc_13.add_to_callback();
        Reg->add(soc_13.in(range(0,0))<=0);
        
        Constraint<> soc_14("soc_14");
        soc_14 = pow(theta23+theta32,2)-(1-theta11-theta22+theta33)*(1-theta11+theta22-theta33);
        //        soc_14.add_to_callback();
        Reg->add(soc_14.in(range(0,0))<=0);
        
        Constraint<> soc_23("soc_23");
        soc_23 = pow(theta23-theta32,2)-(1+theta11-theta22-theta33)*(1+theta11+theta22+theta33);
        //        soc_23.add_to_callback();
        Reg->add(soc_23.in(range(0,0))<=0);
        
        Constraint<> soc_24("soc_24");
        soc_24 = pow(theta12+theta21,2)-(1+theta11-theta22-theta33)*(1-theta11+theta22-theta33);
        //        soc_24.add_to_callback();
        Reg->add(soc_24.in(range(0,0))<=0);
        
        Constraint<> soc_34("soc_34");
        soc_34 = pow(theta31-theta13,2)-(1+theta11+theta22+theta33)*(1-theta11+theta22-theta33);
        //        soc_34.add_to_callback();
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
        if(convex){
            Constraint<> row1("row1");
            row1 = pow(theta11,2)+pow(theta12,2)+pow(theta13,2);
            Reg->add(row1.in(range(0,0))<=1);
            Constraint<> row2("row2");
            row2 = pow(theta21,2)+pow(theta22,2)+pow(theta23,2);
            Reg->add(row2.in(range(0,0))<=1);
            Constraint<> row3("row3");
            row3 = pow(theta31,2)+pow(theta32,2)+pow(theta33,2);
            Reg->add(row3.in(range(0,0))<=1);
            Constraint<> col1("col1");
            col1 = pow(theta11,2)+pow(theta21,2)+pow(theta31,2);
            Reg->add(col1.in(range(0,0))<=1);
            Constraint<> col2("col2");
            col2 = pow(theta12,2)+pow(theta22,2)+pow(theta32,2);
            Reg->add(col2.in(range(0,0))<=1);
            Constraint<> col3("col3");
            col3 = pow(theta13,2)+pow(theta23,2)+pow(theta33,2);
            Reg->add(col3.in(range(0,0))<=1);
        }
        else {
            Constraint<> row1("row1");
            row1 = pow(theta11,2)+pow(theta12,2)+pow(theta13,2);
            Reg->add(row1==1);
            Constraint<> row2("row2");
            row2 = pow(theta21,2)+pow(theta22,2)+pow(theta23,2);
            Reg->add(row2==1);
            Constraint<> row3("row3");
            row3 = pow(theta31,2)+pow(theta32,2)+pow(theta33,2);
            Reg->add(row3==1);
            Constraint<> col1("col1");
            col1 = pow(theta11,2)+pow(theta21,2)+pow(theta31,2);
            Reg->add(col1==1);
            Constraint<> col2("col2");
            col2 = pow(theta12,2)+pow(theta22,2)+pow(theta32,2);
            Reg->add(col2==1);
            Constraint<> col3("col3");
            col3 = pow(theta13,2)+pow(theta23,2)+pow(theta33,2);
            Reg->add(col3==1);
        }
    }
    
    
    
    
    //    for (int i = 1; i<=nd; i++) {
    //        string key = to_string(i)+","+to_string(init_matching.at(i-1)+1);
    //        bin._val->at(bin._indices->_keys_map->at(key)) = 1;
    //    }
    
    /* Objective function */
    if(false && mid_point_lb._range->second>0){
        Constraint<> MidPointLB("MidPointLB");
        MidPointLB = x_diff + y_diff + z_diff - mid_point_lb;
        Reg->add(MidPointLB.in(N1) >= 0);
    }
    
    Reg->min(sum(x_diff + y_diff + z_diff));
    //    bin._val->at(bin._indices->_keys_map->at("1,1")) = 1;
    //    bin.set_lb("1,1",1);
    //    bin.set_lb("2,2",1);
    //    bin.set_lb("3,3",1);
    //    for (int i = 0; i<nd; i++) {
    //        string key = to_string(i+1)+","+to_string(i+1);
    //        bin._val->at(bin._indices->_keys_map->at(key)) = 1;
    //    }
    //    Reg->print();
    //    solver<> S1(Reg,ipopt);
    //    for(int i = 1; i<= nd; i++){
    //        bin.param<int>::set_val(to_string(i)+","+to_string(i), 1);
    //    }
    //            Reg->print();
    if(relax_ints){
        solver<> S(Reg,ipopt);
        S.run();
        Reg->round_solution();
    }
    else {
        solver<> S(Reg,gurobi);
        S.use_callback();
        S.run();
    }
    //    Reg->print();
    //    Reg->print_int_solution();
    
    Reg->print_solution();
    //    {
    //        indices voronoi_ids("voronoi_ids");
    //        voronoi_ids = indices(N1, *norm_x._indices);
    //        auto voronoi_ids_coefs = voronoi_ids.ignore_ith(0, 1);
    ////        auto voronoi_ids_m = voronoi_ids.from_ith(0, 1);
    //        auto voronoi_ids_data = voronoi_ids.ignore_ith(1, 2);
    //        auto voronoi_ids_bin = voronoi_ids.ignore_ith(2, 1);
    ////        Constraint<> Voronoi_model("Voronoi_model");
    ////        Voronoi_model = norm_x.in(voronoi_ids_coefs)*new_xm.in(voronoi_ids_m) + norm_y.in(voronoi_ids_coefs)*new_ym.in(voronoi_ids_m) + norm_z.in(voronoi_ids_coefs)*new_zm.in(voronoi_ids_m) + intercept.in(voronoi_ids_coefs);
    ////        Reg->add_on_off_multivariate_refined(Voronoi_model.in(voronoi_ids)<=0, bin.in(voronoi_ids_bin), true);
    //
    //
    //        auto func = bin.in(voronoi_ids_bin)*(norm_x.in(voronoi_ids_coefs)*new_x1.in(voronoi_ids_data) + norm_y.in(voronoi_ids_coefs)*new_y1.in(voronoi_ids_data) + norm_z.in(voronoi_ids_coefs)*new_z1.in(voronoi_ids_data) + intercept.in(voronoi_ids_coefs));
    //
    //        func.allocate_mem();
    //        func.eval_all();
    //        for (int i = 0; i<func.get_nb_inst(); i++) {
    //            if(func._val->at(i)>1e-6){
    //                DebugOn("instance " <<  i << " is violated \n");
    //                func.print(i,10);
    //                DebugOn(" | violation = " <<  func._val->at(i) << endl);
    //            }
    //        }
    //
    //    }
    DebugOn("row 1 " << pow(theta11.eval(),2)+pow(theta12.eval(),2)+pow(theta13.eval(),2)
            << endl);
    DebugOn("row 2 " << pow(theta21.eval(),2)+pow(theta22.eval(),2)+pow(theta23.eval(),2)
            << endl);
    DebugOn("row 3 " << pow(theta31.eval(),2)+pow(theta32.eval(),2)+pow(theta33.eval(),2)
            << endl);
    DebugOn("col 1 " << pow(theta11.eval(),2)+pow(theta21.eval(),2)+pow(theta31.eval(),2)
            << endl);
    DebugOn("col 2 " << pow(theta12.eval(),2)+pow(theta22.eval(),2)+pow(theta32.eval(),2)
            << endl);
    DebugOn("col 3 " << pow(theta13.eval(),2)+pow(theta23.eval(),2)+pow(theta33.eval(),2)
            << endl);
    
    DebugOn("row 12 " << (theta11.eval()*theta21.eval())+(theta12.eval()*theta22.eval())+(theta13.eval()*theta23.eval())
            << endl);
    DebugOn("row 13 " << (theta11.eval()*theta31.eval())+(theta12.eval()*theta32.eval())+(theta13.eval()*theta33.eval())
            << endl);
    DebugOn("row 23 " << (theta21.eval()*theta31.eval())+(theta22.eval()*theta32.eval())+(theta23.eval()*theta33.eval())
            << endl);
    DebugOn("Theta matrix = " << endl);
    DebugOn("|" << theta11.eval() << " " << theta12.eval() << " " << theta13.eval() << "|" << endl);
    DebugOn("|" << theta21.eval() << " " << theta22.eval() << " " << theta23.eval() << "|" << endl);
    DebugOn("|" << theta31.eval() << " " << theta32.eval() << " " << theta33.eval() << "|" << endl);
    rot_trans.resize(15,0.0);
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
    rot_trans[12]=scale_x.eval();
    rot_trans[13]=scale_y.eval();
    rot_trans[14]=scale_z.eval();
    auto pitch_val = std::atan2(theta32.eval(), theta33.eval())*180/pi;
    auto roll_val = std::atan2(-1*theta31.eval(), std::sqrt(theta32.eval()*theta32.eval()+theta33.eval()*theta33.eval()))*180/pi;
    auto yaw_val = std::atan2(theta21.eval(),theta11.eval())*180/pi;
    DebugOn("Roll (degrees) = " << to_string_with_precision(roll_val,12) << endl);
    DebugOn("Pitch (degrees) = " << to_string_with_precision(pitch_val,12) << endl);
    DebugOn("Yaw (degrees) = " << to_string_with_precision(yaw_val,12) << endl);
    DebugOn("x shift = " << x_shift.eval() << endl);
    DebugOn("y shift = " << y_shift.eval() << endl);
    DebugOn("z shift = " << z_shift.eval() << endl);
    DebugOn("scale_x = " << scale_x.eval() << endl);
    DebugOn("scale_y = " << scale_y.eval() << endl);
    DebugOn("scale_z = " << scale_z.eval() << endl);
    return(Reg);
}
/* Find the farthest point attainable from input point by a 3d Rotation given rotation matrix bounds from rot_mat */
shared_ptr<Model<double>> build_SDP(vector<double>& point, vector<double>& rot_mat){
    
    
    
    param<> x1("x1"), y1("y1"), z1("z1");
    x1 = point[0];
    y1 = point[1];
    z1 = point[2];
    string name="SDP";
    
    auto Reg=make_shared<Model<>>(name);
    
    
    
    var<> theta11("theta11",  std::max(-rot_mat[0],rot_mat[0]), 1), theta12("theta12", std::min(-rot_mat[1],rot_mat[1]), std::max(-rot_mat[1], rot_mat[1])), theta13("theta13", std::min(-rot_mat[2],rot_mat[2]), std::max(-rot_mat[2], rot_mat[2]));
    var<> theta21("theta21", std::min(-rot_mat[3],rot_mat[3]), std::max(-rot_mat[3], rot_mat[3])), theta22("theta22", std::max(-rot_mat[4],rot_mat[4]), 1), theta23("theta23", std::min(-rot_mat[5],rot_mat[5]), std::max(-rot_mat[5], rot_mat[5]));
    var<> theta31("theta31", std::min(-rot_mat[6],rot_mat[6]), std::max(-rot_mat[6], rot_mat[6])), theta32("theta32", std::min(-rot_mat[7],rot_mat[7]), std::max(-rot_mat[7], rot_mat[7])), theta33("theta33", std::max(-rot_mat[8],rot_mat[8]), 1);
    
    
    Reg->add(theta11.in(R(1)),theta12.in(R(1)),theta13.in(R(1)));
    Reg->add(theta21.in(R(1)),theta22.in(R(1)),theta23.in(R(1)));
    Reg->add(theta31.in(R(1)),theta32.in(R(1)),theta33.in(R(1)));
    
    indices N1 = range(1,1);
    param<> x_new_lb("x_new_lb");
    x_new_lb.in(N1);
    param<> x_new_ub("x_new_ub");
    x_new_ub.in(N1);
    param<> y_new_lb("y_new_lb");
    y_new_lb.in(N1);
    param<> y_new_ub("y_new_ub");
    y_new_ub.in(N1);
    param<> z_new_lb("z_new_lb");
    z_new_lb.in(N1);
    param<> z_new_ub("z_new_ub");
    z_new_ub.in(N1);
    
    int nd = 1;
    double x_lb = 0, y_lb = 0, z_lb = 0, x1_i = 0, y1_i = 0, z1_i = 0;
    shared_ptr<pair<double,double>> x1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> y1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> z1_bounds = make_shared<pair<double,double>>();
    for (int i = 0; i<nd; i++) {
        x1_bounds->first = x1.eval(i);
        x1_bounds->second = x1.eval(i);
        y1_bounds->first = y1.eval(i);
        y1_bounds->second = y1.eval(i);
        z1_bounds->first = z1.eval(i);
        z1_bounds->second = z1.eval(i);
        auto x_range  = get_product_range(x1_bounds, theta11._range);
        auto y_range  = get_product_range(y1_bounds, theta12._range);
        auto z_range  = get_product_range(z1_bounds, theta13._range);
        x_new_lb.set_val(i, x_range->first + y_range->first + z_range->first);
        x_new_ub.set_val(i, x_range->second + y_range->second + z_range->second);
        x_range  = get_product_range(x1_bounds, theta21._range);
        y_range  = get_product_range(y1_bounds, theta22._range);
        z_range  = get_product_range(z1_bounds, theta23._range);
        y_new_lb.set_val(i, x_range->first + y_range->first + z_range->first);
        y_new_ub.set_val(i, x_range->second + y_range->second + z_range->second);
        x_range  = get_product_range(x1_bounds, theta31._range);
        y_range  = get_product_range(y1_bounds, theta32._range);
        z_range  = get_product_range(z1_bounds, theta33._range);
        z_new_lb.set_val(i, x_range->first + y_range->first + z_range->first);
        z_new_ub.set_val(i, x_range->second + y_range->second + z_range->second);
    }
    
    //    var<> new_x1("new_x1", x_new_lb, x_new_ub), new_y1("new_y1", y_new_lb, y_new_ub), new_z1("new_z1", z_new_lb, z_new_ub);
    var<> new_x1("new_x1", -1, 1), new_y1("new_y1", -1, 1), new_z1("new_z1", -1, 1);
    Reg->add(new_x1.in(N1), new_y1.in(N1), new_z1.in(N1));
    
    var<int> bin("bin", 0, 1);
    Reg->add(bin.in(N1));
    
    //    var<> delta("delta", pos_);
    //    Reg->add(delta.in(N1));
    
    
    theta11.initialize_all(1);
    theta22.initialize_all(1);
    theta33.initialize_all(1);
    
    
    
    
    
    //    Constraint<> Norm2_new("Norm2_new");
    //    Norm2_new -= delta - (pow(new_x1 - x1,2) + pow(new_y1 - y1,2) + pow(new_z1 - z1,2));
    //    Reg->add(Norm2_new.in(N1)==0);
    
    
    Constraint<> x_rot1("x_rot1");
    x_rot1 += new_x1;
    x_rot1 -= x1*theta11 + y1*theta12 + z1*theta13;
    Reg->add(x_rot1.in(N1)==0);
    
    Constraint<> y_rot1("y_rot1");
    y_rot1 += new_y1;
    y_rot1 -= x1*theta21 + y1*theta22 + z1*theta23;
    Reg->add(y_rot1.in(N1)==0);
    
    Constraint<> z_rot1("z_rot1");
    z_rot1 += new_z1;
    z_rot1 -= x1*theta31 + y1*theta32 + z1*theta33;
    Reg->add(z_rot1.in(N1)==0);
    
    
    
    /* Objective function */
    
    Reg->max((pow(new_x1 - x1,2) + pow(new_y1 - y1,2) + pow(new_z1 - z1,2)) - bin);
    
    
    bool add_sdp_rel = true;
    if(add_sdp_rel){
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
        //        soc_12.add_to_callback();
        Reg->add(soc_12.in(range(0,0))<=0);
        
        Constraint<> soc_13("soc_13");
        soc_13 = pow(theta12-theta21,2)-(1-theta11-theta22+theta33)*(1+theta11+theta22+theta33);
        //        soc_13.add_to_callback();
        Reg->add(soc_13.in(range(0,0))<=0);
        
        Constraint<> soc_14("soc_14");
        soc_14 = pow(theta23+theta32,2)-(1-theta11-theta22+theta33)*(1-theta11+theta22-theta33);
        //        soc_14.add_to_callback();
        Reg->add(soc_14.in(range(0,0))<=0);
        
        Constraint<> soc_23("soc_23");
        soc_23 = pow(theta23-theta32,2)-(1+theta11-theta22-theta33)*(1+theta11+theta22+theta33);
        //        soc_23.add_to_callback();
        Reg->add(soc_23.in(range(0,0))<=0);
        
        Constraint<> soc_24("soc_24");
        soc_24 = pow(theta12+theta21,2)-(1+theta11-theta22-theta33)*(1-theta11+theta22-theta33);
        //        soc_24.add_to_callback();
        Reg->add(soc_24.in(range(0,0))<=0);
        
        Constraint<> soc_34("soc_34");
        soc_34 = pow(theta31-theta13,2)-(1+theta11+theta22+theta33)*(1-theta11+theta22-theta33);
        //        soc_34.add_to_callback();
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
        row1 = pow(theta11,2)+pow(theta12,2)+pow(theta13,2);
        Reg->add(row1==1);
        Constraint<> row2("row2");
        row2 = pow(theta21,2)+pow(theta22,2)+pow(theta23,2);
        Reg->add(row2==1);
        Constraint<> row3("row3");
        row3 = pow(theta31,2)+pow(theta32,2)+pow(theta33,2);
        Reg->add(row3==1);
        Constraint<> col1("col1");
        col1 = pow(theta11,2)+pow(theta21,2)+pow(theta31,2);
        Reg->add(col1==1);
        Constraint<> col2("col2");
        col2 = pow(theta12,2)+pow(theta22,2)+pow(theta32,2);
        Reg->add(col2==1);
        Constraint<> col3("col3");
        col3 = pow(theta13,2)+pow(theta23,2)+pow(theta33,2);
        Reg->add(col3==1);
    }
    
    Reg->print();
    
    solver<> GS(Reg, gurobi);
    GS.run();
    Reg->print_solution();
    
    
    DebugOn("Theta matrix = " << endl);
    DebugOn("|" << theta11.eval() << " " << theta12.eval() << " " << theta13.eval() << "|" << endl);
    DebugOn("|" << theta21.eval() << " " << theta22.eval() << " " << theta23.eval() << "|" << endl);
    DebugOn("|" << theta31.eval() << " " << theta32.eval() << " " << theta33.eval() << "|" << endl);
    constant<> row1 = pow(theta11.eval(),2)+pow(theta12.eval(),2)+pow(theta13.eval(),2);
    constant<> row2 = pow(theta21.eval(),2)+pow(theta22.eval(),2)+pow(theta23.eval(),2);
    constant<> row3 = pow(theta31.eval(),2)+pow(theta32.eval(),2)+pow(theta33.eval(),2);
    constant<> col1 = pow(theta11.eval(),2)+pow(theta21.eval(),2)+pow(theta31.eval(),2);
    constant<> col2 = pow(theta12.eval(),2)+pow(theta22.eval(),2)+pow(theta32.eval(),2);
    constant<> col3 = pow(theta13.eval(),2)+pow(theta23.eval(),2)+pow(theta33.eval(),2);
    DebugOn("row 1 " << row1.eval() << endl);
    DebugOn("row 2 " << row2.eval() << endl);
    DebugOn("row 3 " << row3.eval() << endl);
    DebugOn("col 1 " << col1.eval() << endl);
    DebugOn("col 2 " << col2.eval() << endl);
    DebugOn("col 3 " << col3.eval() << endl);
    constant<> det=theta11.eval()*(theta22.eval()*theta33.eval()-theta32.eval()*theta23.eval())
    -theta12.eval()*(theta21.eval()*theta33.eval()-theta31.eval()*theta23.eval())+theta13.eval()*(theta21.eval()*theta32.eval()-theta31.eval()*theta22.eval());
    constant<> row12 = (theta11.eval()*theta21.eval())+(theta12.eval()*theta22.eval())+(theta13.eval()*theta23.eval());
    constant<> row13 = (theta11.eval()*theta31.eval())+(theta12.eval()*theta32.eval())+(theta13.eval()*theta33.eval());
    constant<> row23 = (theta21.eval()*theta31.eval())+(theta22.eval()*theta32.eval())+(theta23.eval()*theta33.eval());
    DebugOn("row 12 " << row12.eval() << endl);
    DebugOn("row 13 " << row13.eval() << endl);
    DebugOn("row 23 " << row23.eval() << endl);
    
    DebugOn("Determinant "<<det.eval()<<endl);
    
    return(Reg);
}

shared_ptr<Model<double>> build_new_SOC_MIQCP(vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, vector<double>& rot_trans, bool convex, const vector<pair<pair<int,int>,pair<int,int>>>& incompatibles, param<>& norm_x, param<>& norm_y, param<>& norm_z,  param<>& intercept, const vector<int>& init_matching){
    double roll_1 = 0, yaw_1 = 0, pitch_1 = 0;
    int nb_pairs = 0, min_nb_pairs = numeric_limits<int>::max(), max_nb_pairs = 0, av_nb_pairs = 0;
    size_t nm = point_cloud_model.size(), nd = point_cloud_data.size();
    
    double shift_min_x = 0.125, shift_max_x = 0.25, shift_min_y = -0.25,shift_max_y = -0.125,shift_min_z = -0.125,shift_max_z = 0;
    double yaw_min = -12.5*pi/180., yaw_max = 0, pitch_min = 12.5*pi/180.,pitch_max = 25.*pi/180.,roll_min = -12.5*pi/180.,roll_max = 0;
    
    vector<double> zeros = {0,0,0};
    
    param<> x1("x1"), x2("x2"), y1("y1"), y2("y2"), z1("z1"), z2("z2");
    int m = av_nb_pairs;
    string i_str, j_str;
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
    
    
    indices Pairs("Pairs"), cells("cells");
    int idx1 = 0;
    int idx2 = 0;
    indices N1("N1"),N2("N2");
    Debug("nd = " << nd << endl);
    Debug("nm = " << nm << endl);
    
    N1 = range(1,nd);
    N2 = range(1,nm);
    cells = indices(N1,N2);
    string name="SOC_MIQCP";
    
    auto Reg=make_shared<Model<>>(name);
    
    
    //    double shift_min_x = -0.25, shift_max_x = 0.25, shift_min_y = -0.25,shift_max_y = 0.25,shift_min_z = -0.25,shift_max_z = 0.25;
    //    double shift_min_x = 0.23, shift_max_x = 0.24, shift_min_y = -0.24,shift_max_y = -0.23,shift_min_z = -0.02,shift_max_z = -0.01;
    var<> x_shift("x_shift", shift_min_x, shift_max_x), y_shift("y_shift", shift_min_y, shift_max_y), z_shift("z_shift", shift_min_z, shift_max_z);
    //    var<> x_shift("x_shift", 0.23, 0.24), y_shift("y_shift", -0.24, -0.23), z_shift("z_shift", -0.011, -0.01);
    
    Reg->add(x_shift.in(R(1)),y_shift.in(R(1)),z_shift.in(R(1)));
    
    var<int> bin("bin",0,1);
    Reg->add(bin.in(cells));
    DebugOn("Added " << cells.size() << " binary variables" << endl);
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
    
    
    var<> theta11("theta11",  r11._range->first, r11._range->second), theta12("theta12", r12._range->first, r12._range->second), theta13("theta13", r13._range->first, r13._range->second);
    var<> theta21("theta21", r21._range->first, r21._range->second), theta22("theta22", r22._range->first, r22._range->second), theta23("theta23", r23._range->first, r23._range->second);
    var<> theta31("theta31", r31._range->first, r31._range->second), theta32("theta32", r32._range->first, r32._range->second), theta33("theta33", r33._range->first, r33._range->second);
    
    //            var<> theta11("theta11",  0.8, 1), theta12("theta12", -1, 1), theta13("theta13", -1, 1);
    //            var<> theta21("theta21",  -1, 1), theta22("theta22", -1, 1), theta23("theta23", -1, 1);
    //            var<> theta31("theta31",  -1, 1), theta32("theta32", -1, 1), theta33("theta33", 0.8, 1);
    
    //    var<> theta11("theta11",  0.96, 0.97), theta12("theta12", 0.15, 0.16), theta13("theta13", -0.22, -0.2);
    //    var<> theta21("theta21",  -0.21, -0.2), theta22("theta22", 0.95, 0.96), theta23("theta23", -0.24, 0.23);
    //    var<> theta31("theta31",  0.17, 0.18), theta32("theta32", 0.26, 0.27), theta33("theta33", 0.94, 0.95);
    Reg->add(theta11.in(R(1)),theta12.in(R(1)),theta13.in(R(1)));
    Reg->add(theta21.in(R(1)),theta22.in(R(1)),theta23.in(R(1)));
    Reg->add(theta31.in(R(1)),theta32.in(R(1)),theta33.in(R(1)));
    //    Reg->print();
    
    param<> x_new_lb("x_new_lb");
    x_new_lb.in(N1);
    param<> x_new_ub("x_new_ub");
    x_new_ub.in(N1);
    param<> y_new_lb("y_new_lb");
    y_new_lb.in(N1);
    param<> y_new_ub("y_new_ub");
    y_new_ub.in(N1);
    param<> z_new_lb("z_new_lb");
    z_new_lb.in(N1);
    param<> z_new_ub("z_new_ub");
    z_new_ub.in(N1);
    
    
    double x_lb = 0, y_lb = 0, z_lb = 0, x1_i = 0, y1_i = 0, z1_i = 0;
    shared_ptr<pair<double,double>> x1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> y1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> z1_bounds = make_shared<pair<double,double>>();
    for (int i = 0; i<nd; i++) {
        x1_bounds->first = x1.eval(i);
        x1_bounds->second = x1.eval(i);
        y1_bounds->first = y1.eval(i);
        y1_bounds->second = y1.eval(i);
        z1_bounds->first = z1.eval(i);
        z1_bounds->second = z1.eval(i);
        auto x_range  = get_product_range(x1_bounds, theta11._range);
        auto y_range  = get_product_range(y1_bounds, theta12._range);
        auto z_range  = get_product_range(z1_bounds, theta13._range);
        auto bounds = get_min_max(roll_min, roll_max, pitch_min, pitch_max, yaw_min, yaw_max, point_cloud_data[i], zeros);
        auto xlb = x_range->first + y_range->first + z_range->first + x_shift.get_lb().eval();
        auto xub = x_range->second + y_range->second + z_range->second+ x_shift.get_ub().eval();
        x_new_lb.set_val(i,bounds[0].first + x_shift.get_lb().eval());
        x_new_ub.set_val(i,bounds[0].second + x_shift.get_ub().eval());
        //        x_new_lb.set_val(i, x_range->first + y_range->first + z_range->first + x_shift.get_lb().eval());
        //        x_new_ub.set_val(i, x_range->second + y_range->second + z_range->second+ x_shift.get_ub().eval());
        x_range  = get_product_range(x1_bounds, theta21._range);
        y_range  = get_product_range(y1_bounds, theta22._range);
        z_range  = get_product_range(z1_bounds, theta23._range);
        auto ylb = x_range->first + y_range->first + z_range->first + y_shift.get_lb().eval();
        auto yub = x_range->second + y_range->second + z_range->second + y_shift.get_ub().eval();
        y_new_lb.set_val(i,bounds[1].first + y_shift.get_lb().eval());
        y_new_ub.set_val(i,bounds[1].second + y_shift.get_ub().eval());
        //        y_new_lb.set_val(i, x_range->first + y_range->first + z_range->first + y_shift.get_lb().eval());
        //        y_new_ub.set_val(i, x_range->second + y_range->second + z_range->second+ y_shift.get_ub().eval());
        x_range  = get_product_range(x1_bounds, theta31._range);
        y_range  = get_product_range(y1_bounds, theta32._range);
        z_range  = get_product_range(z1_bounds, theta33._range);
        auto zlb = x_range->first + y_range->first + z_range->first + z_shift.get_lb().eval();
        auto zub = x_range->second + y_range->second + z_range->second + z_shift.get_ub().eval();
        z_new_lb.set_val(i,bounds[2].first + z_shift.get_lb().eval());
        z_new_ub.set_val(i,bounds[2].second + z_shift.get_ub().eval());
        //        z_new_lb.set_val(i, x_range->first + y_range->first + z_range->first + z_shift.get_lb().eval());
        //        z_new_ub.set_val(i, x_range->second + y_range->second + z_range->second+ z_shift.get_ub().eval());
    }
    
    var<> new_xm("new_xm", -1, 1), new_ym("new_ym", -1, 1), new_zm("new_zm", -1, 1);
    var<> new_x1("new_x1", x_new_lb, x_new_ub), new_y1("new_y1", y_new_lb, y_new_ub), new_z1("new_z1", z_new_lb, z_new_ub);
    //            var<> new_x1("new_x1", -1, 1), new_y1("new_y1", -1, 1), new_z1("new_z1", -1, 1);
    Reg->add(new_xm.in(N1), new_ym.in(N1), new_zm.in(N1));
    Reg->add(new_x1.in(N1), new_y1.in(N1), new_z1.in(N1));
    
    //    var<> new_x1_ij("new_x1_ij", -1, 1), new_y1_ij("new_y1_ij", -1, 1), new_z1_ij("new_z1_ij", -1, 1);
    //    Reg->add(new_x1_ij.in(cells), new_y1_ij.in(cells), new_z1_ij.in(cells));
    
    
    var<> delta("delta", pos_);
    Reg->add(delta.in(N1));
    
    
    //    var<> d1("d1", 0,4),d2("d2", 0,4),d3("d3", 0,4),d4("d4", 0,4);
    //    var<> l12("l12", -2,2),l13("l13", -2,2),l14("l14", -2,2),l23("l23", -2,2),l24("l24", -2,2),l34("l34", -2,2);
    //    Reg->add(l12.in(R(1)),l13.in(R(1)),l14.in(R(1)));
    //    Reg->add(l23.in(R(1)),l24.in(R(1)),l34.in(R(1)));
    //    Reg->add(d1.in(R(1)),d2.in(R(1)),d3.in(R(1)),d4.in(R(1)));
    
    
    indices ids = indices("in_x");
    ids.add_empty_row();
    for(auto i=0;i<nd;i++){
        for(auto j=1;j<=nm;j++){
            if(cells.has_key(to_string(i+1)+","+to_string(j)))
                ids.add_in_row(i, to_string(j));
        }
    }
    
    bool add_shift_cut = false;
    
    if(add_shift_cut){
        Constraint<> Centroid("Centroid");
        Centroid = sum(new_x1) + sum(new_y1) + sum(new_z1) - (sum(new_xm) + sum(new_ym) + sum(new_zm));
        Reg->add(Centroid==0);
        
        Constraint<> Centroid_tx("Centroid_tx");
        Centroid_tx = nd*x_shift - sum(new_xm);
        Reg->add(Centroid_tx==0);
        
        Constraint<> Centroid_ty("Centroid_ty");
        Centroid_ty = nd*y_shift - sum(new_ym);
        Reg->add(Centroid_ty==0);
        
        Constraint<> Centroid_tz("Centroid_tz");
        Centroid_tz = nd*z_shift - sum(new_zm);
        Reg->add(Centroid_tz==0);
    }
    
    bool add_voronoi = false;
    
    if(add_voronoi){
        indices voronoi_ids("voronoi_ids");
        voronoi_ids = indices(range(1,3), *norm_x._indices);
        auto voronoi_ids_coefs = voronoi_ids.ignore_ith(0, 1);
        auto voronoi_ids_m = voronoi_ids.from_ith(0, 1);
        auto voronoi_ids_data = voronoi_ids.ignore_ith(1, 2);
        auto voronoi_ids_bin = voronoi_ids.ignore_ith(2, 1);
        //        Constraint<> Voronoi_model("Voronoi_model");
        //        Voronoi_model = norm_x.in(voronoi_ids_coefs)*new_xm.in(voronoi_ids_m) + norm_y.in(voronoi_ids_coefs)*new_ym.in(voronoi_ids_m) + norm_z.in(voronoi_ids_coefs)*new_zm.in(voronoi_ids_m) + intercept.in(voronoi_ids_coefs);
        //        Reg->add_on_off_multivariate_refined(Voronoi_model.in(voronoi_ids)<=0, bin.in(voronoi_ids_bin), true);
        
        
        Constraint<> Voronoi("Voronoi");
        Voronoi = norm_x.in(voronoi_ids_coefs)*new_x1.in(voronoi_ids_data) + norm_y.in(voronoi_ids_coefs)*new_y1.in(voronoi_ids_data) + norm_z.in(voronoi_ids_coefs)*new_z1.in(voronoi_ids_data) + intercept.in(voronoi_ids_coefs);
        Reg->add_on_off_multivariate_refined(Voronoi.in(voronoi_ids)<=0, bin.in(voronoi_ids_bin), true);
        //        Reg->print();
    }
    
    
    //    Reg->print();
    
    
    //
    theta11.initialize_all(1);
    theta22.initialize_all(1);
    theta33.initialize_all(1);
    
    bool spatial_branching = false;
    if(spatial_branching){
        /* Spatial branching vars */
        int nb_pieces = 5; // Divide each axis into nb_pieces
        indices spatial_ids("spatial_ids");
        spatial_ids = range(1,nb_pieces);
        indices theta_ids("theta_ids");
        theta_ids = indices(range(1,3),range(1,3));
        indices shift_ids("shift_ids");
        shift_ids.insert({"x", "y", "z"});
        indices theta_spatial("theta_spatial");
        theta_spatial = indices(theta_ids, spatial_ids);
        
        indices shift_spatial("shift_spatial");
        shift_spatial = indices(shift_ids, spatial_ids);
        
        //    var<int> sbin_shift("sbin_shift", 0, 1);
        //    var<int> sbin_theta("sbin_theta", 0, 1);
        
        //    Reg->add(sbin_theta.in(theta_spatial),sbin_shift.in(shift_spatial));
        
        var<int> sbin_tx("sbin_tx", 0, 1), sbin_ty("sbin_ty", 0, 1), sbin_tz("sbin_tz", 0, 1);
        var<int> sbin_theta11("sbin_theta11", 0, 1), sbin_theta12("sbin_theta12", 0, 1), sbin_theta13("sbin_theta13", 0, 1);
        var<int> sbin_theta21("sbin_theta21", 0, 1), sbin_theta22("sbin_theta22", 0, 1), sbin_theta23("sbin_theta23", 0, 1);
        var<int> sbin_theta31("sbin_theta31", 0, 1), sbin_theta32("sbin_theta32", 0, 1), sbin_theta33("sbin_theta33", 0, 1);
        Reg->add(sbin_theta11.in(spatial_ids),sbin_theta12.in(spatial_ids),sbin_theta13.in(spatial_ids));
        Reg->add(sbin_theta21.in(spatial_ids),sbin_theta22.in(spatial_ids),sbin_theta23.in(spatial_ids));
        Reg->add(sbin_theta31.in(spatial_ids),sbin_theta32.in(spatial_ids),sbin_theta33.in(spatial_ids));
        Reg->add(sbin_tx.in(spatial_ids),sbin_ty.in(spatial_ids),sbin_tz.in(spatial_ids));
        /* Spatial branching constraints */
        Constraint<> OneBinAngleSpatial11("OneBinAngleSpatial11");
        OneBinAngleSpatial11 = sum(sbin_theta11);
        Reg->add(OneBinAngleSpatial11==1);
        
        Constraint<> OneBinAngleSpatial12("OneBinAngleSpatial12");
        OneBinAngleSpatial12 = sum(sbin_theta12);
        Reg->add(OneBinAngleSpatial12==1);
        
        Constraint<> OneBinAngleSpatial13("OneBinAngleSpatial13");
        OneBinAngleSpatial13 = sum(sbin_theta13);
        Reg->add(OneBinAngleSpatial13==1);
        
        Constraint<> OneBinAngleSpatial21("OneBinAngleSpatial21");
        OneBinAngleSpatial21 = sum(sbin_theta21);
        Reg->add(OneBinAngleSpatial21==1);
        
        Constraint<> OneBinAngleSpatial22("OneBinAngleSpatial22");
        OneBinAngleSpatial22 = sum(sbin_theta22);
        Reg->add(OneBinAngleSpatial22==1);
        
        
        Constraint<> OneBinAngleSpatial23("OneBinAngleSpatial23");
        OneBinAngleSpatial23 = sum(sbin_theta23);
        Reg->add(OneBinAngleSpatial23==1);
        
        Constraint<> OneBinAngleSpatial31("OneBinAngleSpatial31");
        OneBinAngleSpatial31 = sum(sbin_theta31);
        Reg->add(OneBinAngleSpatial31==1);
        
        Constraint<> OneBinAngleSpatial32("OneBinAngleSpatial32");
        OneBinAngleSpatial32 = sum(sbin_theta32);
        Reg->add(OneBinAngleSpatial32==1);
        
        Constraint<> OneBinAngleSpatial33("OneBinAngleSpatial33");
        OneBinAngleSpatial33 = sum(sbin_theta33);
        Reg->add(OneBinAngleSpatial33==1);
        
        Constraint<> OneBinShiftSpatialx("OneBinShiftSpatialx");
        OneBinShiftSpatialx = sum(sbin_tx);
        Reg->add(OneBinShiftSpatialx==1);
        
        Constraint<> OneBinShiftSpatialy("OneBinShiftSpatialy");
        OneBinShiftSpatialy = sum(sbin_ty);
        Reg->add(OneBinShiftSpatialy==1);
        
        Constraint<> OneBinShiftSpatialz("OneBinShiftSpatialz");
        OneBinShiftSpatialz = sum(sbin_tz);
        Reg->add(OneBinShiftSpatialz==1);
        
        //    Reg->print();
        
        double diag_increment = 1./nb_pieces;/* Diagonals are defined in [0,1] */
        double off_diag_increment = 2./nb_pieces;/* Diagonals are defined in [-1,1] */
        double shift_increment = 0.5/nb_pieces;/* Shifts are defined in [-0.25,0.25] */
        
        auto spatial_ids_n = range(1,nb_pieces-1);
        auto spatial_ids_1 = range(2,nb_pieces);
        param<> diag_lb("diag_lb"), diag_ub("diag_ub"), off_diag_lb("off_diag_lb"), off_diag_ub("off_diag_ub");
        param<> t_lb("t_lb"), t_ub("t_ub");
        diag_ub.in(spatial_ids);
        diag_lb.in(spatial_ids);
        off_diag_ub.in(spatial_ids);
        off_diag_lb.in(spatial_ids);
        t_ub.in(spatial_ids);
        t_lb.in(spatial_ids);
        for (int i = 0; i<nb_pieces; i++) {
            diag_ub.set_val(i,(i+1)*diag_increment);
            off_diag_ub.set_val(i,-1+(i+1)*off_diag_increment);
            t_ub.set_val(i,-0.25 + (i+1)*diag_increment);
            diag_lb.set_val(i,i*diag_increment);
            off_diag_lb.set_val(i,-1 + i*off_diag_increment);
            t_lb.set_val(i,-0.25 + i*diag_increment);
        }
        auto ids_repeat = theta11.repeat_id(nb_pieces-1);
        
        Constraint<> Diag_Spatial_UB11("Diag_Spatial_UB11");
        Diag_Spatial_UB11 = theta11.in(ids_repeat) - diag_ub.in(spatial_ids_n);
        Reg->add_on_off_multivariate_refined(Diag_Spatial_UB11.in(spatial_ids_n)<=0, sbin_theta11.in(spatial_ids_n), true);
        
        Constraint<> Diag_Spatial_LB11("Diag_Spatial_LB11");
        Diag_Spatial_LB11 = diag_lb.in(spatial_ids_1) - theta11.in(ids_repeat);
        Reg->add_on_off_multivariate_refined(Diag_Spatial_LB11.in(spatial_ids_1) <= 0, sbin_theta11.in(spatial_ids_1), true);
        
        Constraint<> Diag_Spatial_UB22("Diag_Spatial_UB22");
        Diag_Spatial_UB22 = theta22.in(ids_repeat) - diag_ub.in(spatial_ids_n);
        Reg->add_on_off_multivariate_refined(Diag_Spatial_UB22.in(spatial_ids_n)<=0, sbin_theta22.in(spatial_ids_n), true);
        
        Constraint<> Diag_Spatial_LB33("Diag_Spatial_LB33");
        Diag_Spatial_LB33 = diag_lb.in(spatial_ids_1) - theta33.in(ids_repeat);
        Reg->add_on_off_multivariate_refined(Diag_Spatial_LB33.in(spatial_ids_1) <= 0, sbin_theta33.in(spatial_ids_1), true);
        
        
        Constraint<> Diag_Spatial_UB12("Diag_Spatial_UB12");
        Diag_Spatial_UB12 = theta12.in(ids_repeat) - off_diag_ub.in(spatial_ids_n);
        Reg->add_on_off_multivariate_refined(Diag_Spatial_UB12.in(spatial_ids_n) <= 0, sbin_theta12.in(spatial_ids_n), true);
        
        Constraint<> Diag_Spatial_LB12("Diag_Spatial_LB12");
        Diag_Spatial_LB12 = off_diag_lb.in(spatial_ids_1) - theta12.in(ids_repeat);
        Reg->add_on_off_multivariate_refined(Diag_Spatial_LB12.in(spatial_ids_1) <= 0, sbin_theta12.in(spatial_ids_1), true);
        
        Constraint<> Diag_Spatial_UB13("Diag_Spatial_UB13");
        Diag_Spatial_UB13 = theta13.in(ids_repeat) - off_diag_ub.in(spatial_ids_n);
        Reg->add_on_off_multivariate_refined(Diag_Spatial_UB13.in(spatial_ids_n) <= 0, sbin_theta13.in(spatial_ids_n), true);
        
        Constraint<> Diag_Spatial_LB13("Diag_Spatial_LB13");
        Diag_Spatial_LB13 = off_diag_lb.in(spatial_ids_1) - theta13.in(ids_repeat);
        Reg->add_on_off_multivariate_refined(Diag_Spatial_LB13.in(spatial_ids_1) <= 0, sbin_theta13.in(spatial_ids_1), true);
        
        Constraint<> Diag_Spatial_UB21("Diag_Spatial_UB21");
        Diag_Spatial_UB21 = theta21.in(ids_repeat) - off_diag_ub.in(spatial_ids_n);
        Reg->add_on_off_multivariate_refined(Diag_Spatial_UB21.in(spatial_ids_n) <= 0, sbin_theta21.in(spatial_ids_n), true);
        
        Constraint<> Diag_Spatial_LB21("Diag_Spatial_LB21");
        Diag_Spatial_LB21 = off_diag_lb.in(spatial_ids_1) - theta21.in(ids_repeat);
        Reg->add_on_off_multivariate_refined(Diag_Spatial_LB21.in(spatial_ids_1) <= 0, sbin_theta21.in(spatial_ids_1), true);
        
        Constraint<> Diag_Spatial_UB23("Diag_Spatial_UB23");
        Diag_Spatial_UB23 = theta23.in(ids_repeat) - off_diag_ub.in(spatial_ids_n);
        Reg->add_on_off_multivariate_refined(Diag_Spatial_UB23.in(spatial_ids_n) <= 0, sbin_theta23.in(spatial_ids_n), true);
        
        Constraint<> Diag_Spatial_LB23("Diag_Spatial_LB23");
        Diag_Spatial_LB23 = off_diag_lb.in(spatial_ids_1) - theta23.in(ids_repeat);
        Reg->add_on_off_multivariate_refined(Diag_Spatial_LB23.in(spatial_ids_1) <= 0, sbin_theta23.in(spatial_ids_1), true);
        
        Constraint<> Diag_Spatial_UB31("Diag_Spatial_UB31");
        Diag_Spatial_UB31 = theta31.in(ids_repeat) - off_diag_ub.in(spatial_ids_n);
        Reg->add_on_off_multivariate_refined(Diag_Spatial_UB31.in(spatial_ids_n) <= 0, sbin_theta31.in(spatial_ids_n), true);
        
        Constraint<> Diag_Spatial_LB31("Diag_Spatial_LB31");
        Diag_Spatial_LB31 = off_diag_lb.in(spatial_ids_1) - theta31.in(ids_repeat);
        Reg->add_on_off_multivariate_refined(Diag_Spatial_LB31.in(spatial_ids_1) <= 0, sbin_theta31.in(spatial_ids_1), true);
        
        Constraint<> Diag_Spatial_UB32("Diag_Spatial_UB32");
        Diag_Spatial_UB32 = theta32.in(ids_repeat) - off_diag_ub.in(spatial_ids_n);
        Reg->add_on_off_multivariate_refined(Diag_Spatial_UB32.in(spatial_ids_n) <= 0, sbin_theta32.in(spatial_ids_n), true);
        
        Constraint<> Diag_Spatial_LB32("Diag_Spatial_LB32");
        Diag_Spatial_LB32 = off_diag_lb.in(spatial_ids_1) - theta32.in(ids_repeat);
        Reg->add_on_off_multivariate_refined(Diag_Spatial_LB32.in(spatial_ids_1) <= 0, sbin_theta32.in(spatial_ids_1), true);
        
        Constraint<> xshift_Spatial_UB("xshift_Spatial_UB");
        xshift_Spatial_UB = x_shift.in(ids_repeat) - t_ub.in(spatial_ids_n);
        Reg->add_on_off_multivariate_refined(xshift_Spatial_UB.in(spatial_ids_n) <= 0, sbin_tx.in(spatial_ids_n), true);
        
        Constraint<> xshift_Spatial_LB("xshift_Spatial_LB");
        xshift_Spatial_LB = t_lb.in(spatial_ids_1) - x_shift.in(ids_repeat);
        Reg->add_on_off_multivariate_refined(xshift_Spatial_LB.in(spatial_ids_1) <= 0, sbin_tx.in(spatial_ids_1), true);
        
        Constraint<> yshift_Spatial_UB("yshift_Spatial_UB");
        yshift_Spatial_UB = y_shift.in(ids_repeat) - t_ub.in(spatial_ids_n);
        Reg->add_on_off_multivariate_refined(yshift_Spatial_UB.in(spatial_ids_n) <= 0, sbin_ty.in(spatial_ids_n), true);
        
        Constraint<> yshift_Spatial_LB("yshift_Spatial_LB");
        yshift_Spatial_LB = t_lb.in(spatial_ids_1) - y_shift.in(ids_repeat);
        Reg->add_on_off_multivariate_refined(yshift_Spatial_LB.in(spatial_ids_1) <= 0, sbin_ty.in(spatial_ids_1), true);
        
        Constraint<> zshift_Spatial_UB("zshift_Spatial_UB");
        zshift_Spatial_UB = z_shift.in(ids_repeat) - t_ub.in(spatial_ids_n);
        Reg->add_on_off_multivariate_refined(zshift_Spatial_UB.in(spatial_ids_n) <= 0, sbin_tz.in(spatial_ids_n), true);
        
        Constraint<> zshift_Spatial_LB("zshift_Spatial_LB");
        zshift_Spatial_LB = t_lb.in(spatial_ids_1) - z_shift.in(ids_repeat);
        Reg->add_on_off_multivariate_refined(zshift_Spatial_LB.in(spatial_ids_1) <= 0, sbin_tz.in(spatial_ids_1), true);
        //        Reg->print();
    }
    
    
    Constraint<> Def_newxm("Def_newxm");
    Def_newxm = new_xm-product(x2.in(ids),bin.in_matrix(1, 1));
    Reg->add(Def_newxm.in(N1)==0);
    
    Constraint<> Def_newym("Def_newym");
    Def_newym = new_ym-product(y2.in(ids),bin.in_matrix(1, 1));
    Reg->add(Def_newym.in(N1)==0);
    
    Constraint<> Def_newzm("Def_newzm");
    Def_newzm = new_zm-product(z2.in(ids),bin.in_matrix(1, 1));
    Reg->add(Def_newzm.in(N1)==0);
    
    Constraint<> OneBin("OneBin");
    OneBin = bin.in_matrix(1, 1);
    Reg->add(OneBin.in(N1)==1);
    
    //    Constraint<> OneBin2("OneBin2");
    //    OneBin2 = bin.in_matrix(0, 1);
    //    Reg->add(OneBin2.in(N2)<=1);
    
    if(false && !incompatibles.empty()){
        indices pairs1("pairs1"), pairs2("pairs2");
        pairs1 = cells;
        pairs2 = cells;
        for (const auto &inc_pair : incompatibles) {
            pairs1.add_ref(to_string(inc_pair.first.first+1)+","+to_string(inc_pair.second.first+1));
            pairs2.add_ref(to_string(inc_pair.first.second+1)+","+to_string(inc_pair.second.second+1));
        }
        
        Constraint<> incomp_pairs("incomp_pairs");
        incomp_pairs = bin.in(pairs1) + bin.in(pairs2);
        Reg->add(incomp_pairs.in(range(1,pairs1.size()))<=1);
        //        incomp_pairs.print();
    }
    
    
    
    double scaling = 1;
    Constraint<> Norm2_new("Norm2_new");
    Norm2_new -= delta - scaling*(pow(new_x1 - new_xm,2) + pow(new_y1 - new_ym,2) + pow(new_z1 - new_zm,2));
    Reg->add(Norm2_new.in(N1)<=0);
    
    bool add_delta_cut = true;
    if(add_delta_cut){
        var<> new_x1_sqr("new_x1_sqr", 0, max(pow(new_x1.get_lb(),2),pow(new_x1.get_ub(),2)));
        var<> new_y1_sqr("new_y1_sqr", 0, max(pow(new_y1.get_lb(),2),pow(new_y1.get_ub(),2)));
        var<> new_z1_sqr("new_z1_sqr", 0, max(pow(new_z1.get_lb(),2),pow(new_y1.get_ub(),2)));
        Reg->add(new_x1_sqr.in(N1),new_y1_sqr.in(N1),new_z1_sqr.in(N1));
        
        bool split = true, convexify = true;
        Constraint<> new_x1_Square("new_x1_Square");
        new_x1_Square -= new_x1_sqr - pow(new_x1,2);
        //        Reg->add(new_x1_Square.in(N1)==0,convexify,"on/off",split);/* Convexify and split nonconvex equation */
        Reg->add(new_x1_Square.in(N1)<=0);
        
        Constraint<> new_y1_Square("new_y1_Square");
        new_y1_Square -= new_y1_sqr - pow(new_y1,2);
        //        Reg->add(new_y1_Square.in(N1)==0,convexify,"on/off",split);/* Convexify and split nonconvex equation */
        Reg->add(new_y1_Square.in(N1)<=0);
        
        Constraint<> new_z1_Square("new_z1_Square");
        new_z1_Square -= new_z1_sqr - pow(new_z1,2);
        //        Reg->add(new_z1_Square.in(N1)==0,convexify,"on/off",split);/* Convexify and split nonconvex equation */
        Reg->add(new_z1_Square.in(N1)<=0);
        
        Constraint<> DeltaCut("DeltaCut");
        DeltaCut -= delta.from(cells);
        DeltaCut += scaling*(pow(x2.to(cells),2) + pow(y2.to(cells),2) + pow(z2.to(cells),2));
        DeltaCut -= scaling*(2*new_x1.from(cells)*x2.to(cells) + 2*new_y1.from(cells)*y2.to(cells) + 2*new_z1.from(cells)*z2.to(cells));
        DeltaCut += scaling*(new_x1_sqr.from(cells) + new_y1_sqr.from(cells) + new_z1_sqr.from(cells));
        Reg->add_on_off_multivariate_refined(DeltaCut.in(cells)<=0, bin, true);
    }
    bool add_delta_ij = false;
    
    if(add_delta_ij) {
        var<> delta_ij("delta_ij", 0, 12);
        Reg->add(delta_ij.in(cells));
        
        Constraint<> Norm2("Norm2");
        Norm2 -= delta_ij - pow(new_x1.from(cells) - x2.to(cells),2) - pow(new_y1.from(cells) - y2.to(cells),2) - pow(new_z1.from(cells) - z2.to(cells),2);
        Reg->add(Norm2.in(cells)<=0);
        
        
        Constraint<> DeltaMin("DeltaMin");
        DeltaMin -= delta.from(cells);
        DeltaMin += delta_ij;
        Reg->add_on_off_multivariate_refined(DeltaMin.in(cells)<=0, bin, true);
    }
    
    bool add_distance_cut = true;
    if(add_distance_cut){
        
    }
    //    Reg->print();
    
    bool add_sdp_rel = true;
    if(add_sdp_rel){
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
        if(convex){
            Constraint<> row1("row1");
            row1 = pow(theta11,2)+pow(theta12,2)+pow(theta13,2);
            Reg->add(row1.in(range(0,0))<=1);
            Constraint<> row2("row2");
            row2 = pow(theta21,2)+pow(theta22,2)+pow(theta23,2);
            Reg->add(row2.in(range(0,0))<=1);
            Constraint<> row3("row3");
            row3 = pow(theta31,2)+pow(theta32,2)+pow(theta33,2);
            Reg->add(row3.in(range(0,0))<=1);
            Constraint<> col1("col1");
            col1 = pow(theta11,2)+pow(theta21,2)+pow(theta31,2);
            Reg->add(col1.in(range(0,0))<=1);
            Constraint<> col2("col2");
            col2 = pow(theta12,2)+pow(theta22,2)+pow(theta32,2);
            Reg->add(col2.in(range(0,0))<=1);
            Constraint<> col3("col3");
            col3 = pow(theta13,2)+pow(theta23,2)+pow(theta33,2);
            Reg->add(col3.in(range(0,0))<=1);
        }
        else {
            Constraint<> row1("row1");
            row1 = pow(theta11,2)+pow(theta12,2)+pow(theta13,2);
            Reg->add(row1==1);
            Constraint<> row2("row2");
            row2 = pow(theta21,2)+pow(theta22,2)+pow(theta23,2);
            Reg->add(row2==1);
            Constraint<> row3("row3");
            row3 = pow(theta31,2)+pow(theta32,2)+pow(theta33,2);
            Reg->add(row3==1);
            Constraint<> col1("col1");
            col1 = pow(theta11,2)+pow(theta21,2)+pow(theta31,2);
            Reg->add(col1==1);
            Constraint<> col2("col2");
            col2 = pow(theta12,2)+pow(theta22,2)+pow(theta32,2);
            Reg->add(col2==1);
            Constraint<> col3("col3");
            col3 = pow(theta13,2)+pow(theta23,2)+pow(theta33,2);
            Reg->add(col3==1);
        }
    }
    auto ids1 = theta11.repeat_id(cells.size());
    Constraint<> x_rot1("x_rot1");
    x_rot1 += new_x1 -x_shift;
    x_rot1 -= x1.in(N1)*theta11.in(ids1) + y1.in(N1)*theta12.in(ids1) + z1.in(N1)*theta13.in(ids1);
    Reg->add(x_rot1.in(N1)==0);
    
    Constraint<> y_rot1("y_rot1");
    y_rot1 += new_y1 - y_shift;
    y_rot1 -= x1.in(N1)*theta21.in(ids1) + y1.in(N1)*theta22.in(ids1) + z1.in(N1)*theta23.in(ids1);
    Reg->add(y_rot1.in(N1)==0);
    
    Constraint<> z_rot1("z_rot1");
    z_rot1 += new_z1 -z_shift;
    z_rot1 -= x1.in(N1)*theta31.in(ids1) + y1.in(N1)*theta32.in(ids1) + z1.in(N1)*theta33.in(ids1);
    Reg->add(z_rot1.in(N1)==0);
    
    
    
    for (int i = 1; i<=nd; i++) {
        string key = to_string(i)+","+to_string(init_matching.at(i-1)+1);
        bin._val->at(bin._indices->_keys_map->at(key)) = 1;
    }
    
    /* Objective function */
    
    Reg->min(sum(delta));
    //    bin._val->at(bin._indices->_keys_map->at("1,1")) = 1;
    //    bin.set_lb("1,1",1);
    //    bin.set_lb("2,2",1);
    //    bin.set_lb("3,3",1);
    //    for (int i = 0; i<nd; i++) {
    //        string key = to_string(i+1)+","+to_string(i+1);
    //        bin._val->at(bin._indices->_keys_map->at(key)) = 1;
    //    }
    //    Reg->print();
    //    solver<> S1(Reg,ipopt);
    for(int i = 1; i<= nd; i++){
        bin.param<int>::set_val(to_string(i)+","+to_string(i), 1);
    }
    Reg->print();
    //    solver<> S(Reg,gurobi);
    //    S.run();
    //        //    Reg->print();
    //    Reg->print_int_solution();
    //
    //    Reg->print_solution();
    //    {
    //        indices voronoi_ids("voronoi_ids");
    //        voronoi_ids = indices(N1, *norm_x._indices);
    //        auto voronoi_ids_coefs = voronoi_ids.ignore_ith(0, 1);
    ////        auto voronoi_ids_m = voronoi_ids.from_ith(0, 1);
    //        auto voronoi_ids_data = voronoi_ids.ignore_ith(1, 2);
    //        auto voronoi_ids_bin = voronoi_ids.ignore_ith(2, 1);
    ////        Constraint<> Voronoi_model("Voronoi_model");
    ////        Voronoi_model = norm_x.in(voronoi_ids_coefs)*new_xm.in(voronoi_ids_m) + norm_y.in(voronoi_ids_coefs)*new_ym.in(voronoi_ids_m) + norm_z.in(voronoi_ids_coefs)*new_zm.in(voronoi_ids_m) + intercept.in(voronoi_ids_coefs);
    ////        Reg->add_on_off_multivariate_refined(Voronoi_model.in(voronoi_ids)<=0, bin.in(voronoi_ids_bin), true);
    //
    //
    //        auto func = bin.in(voronoi_ids_bin)*(norm_x.in(voronoi_ids_coefs)*new_x1.in(voronoi_ids_data) + norm_y.in(voronoi_ids_coefs)*new_y1.in(voronoi_ids_data) + norm_z.in(voronoi_ids_coefs)*new_z1.in(voronoi_ids_data) + intercept.in(voronoi_ids_coefs));
    //
    //        func.allocate_mem();
    //        func.eval_all();
    //        for (int i = 0; i<func.get_nb_inst(); i++) {
    //            if(func._val->at(i)>1e-6){
    //                DebugOn("instance " <<  i << " is violated \n");
    //                func.print(i,10);
    //                DebugOn(" | violation = " <<  func._val->at(i) << endl);
    //            }
    //        }
    //
    //    }
    //    DebugOn("row 1 " << pow(theta11.eval(),2)+pow(theta12.eval(),2)+pow(theta13.eval(),2)
    //            << endl);
    //    DebugOn("row 2 " << pow(theta21.eval(),2)+pow(theta22.eval(),2)+pow(theta23.eval(),2)
    //            << endl);
    //    DebugOn("row 3 " << pow(theta31.eval(),2)+pow(theta32.eval(),2)+pow(theta33.eval(),2)
    //            << endl);
    //    DebugOn("col 1 " << pow(theta11.eval(),2)+pow(theta21.eval(),2)+pow(theta31.eval(),2)
    //            << endl);
    //    DebugOn("col 2 " << pow(theta12.eval(),2)+pow(theta22.eval(),2)+pow(theta32.eval(),2)
    //            << endl);
    //    DebugOn("col 3 " << pow(theta13.eval(),2)+pow(theta23.eval(),2)+pow(theta33.eval(),2)
    //            << endl);
    //
    //    DebugOn("row 12 " << (theta11.eval()*theta21.eval())+(theta12.eval()*theta22.eval())+(theta13.eval()*theta23.eval())
    //            << endl);
    //    DebugOn("row 13 " << (theta11.eval()*theta31.eval())+(theta12.eval()*theta32.eval())+(theta13.eval()*theta33.eval())
    //            << endl);
    //    DebugOn("row 23 " << (theta21.eval()*theta31.eval())+(theta22.eval()*theta32.eval())+(theta23.eval()*theta33.eval())
    //            << endl);
    //    DebugOn("Theta matrix = " << endl);
    //    DebugOn("|" << theta11.eval() << " " << theta12.eval() << " " << theta13.eval() << "|" << endl);
    //    DebugOn("|" << theta21.eval() << " " << theta22.eval() << " " << theta23.eval() << "|" << endl);
    //    DebugOn("|" << theta31.eval() << " " << theta32.eval() << " " << theta33.eval() << "|" << endl);
    //    rot_trans[0]=theta11.eval();
    //    rot_trans[1]=theta12.eval();
    //    rot_trans[2]=theta13.eval();;
    //    rot_trans[3]=theta21.eval();
    //    rot_trans[4]=theta22.eval();
    //    rot_trans[5]=theta23.eval();
    //    rot_trans[6]=theta31.eval();
    //    rot_trans[7]=theta32.eval();
    //    rot_trans[8]=theta33.eval();
    //    rot_trans[9]=x_shift.eval();
    //    rot_trans[10]=y_shift.eval();
    //    rot_trans[11]=z_shift.eval();
    //    auto pitch_val = std::atan2(theta32.eval(), theta33.eval())*180/pi;
    //    auto roll_val = std::atan2(-1*theta31.eval(), std::sqrt(theta32.eval()*theta32.eval()+theta33.eval()*theta33.eval()))*180/pi;
    //    auto yaw_val = std::atan2(theta21.eval(),theta11.eval())*180/pi;
    //    DebugOn("Roll (degrees) = " << to_string_with_precision(roll_val,12) << endl);
    //    DebugOn("Pitch (degrees) = " << to_string_with_precision(pitch_val,12) << endl);
    //    DebugOn("Yaw (degrees) = " << to_string_with_precision(yaw_val,12) << endl);
    //    DebugOn("x shift = " << x_shift.eval() << endl);
    //    DebugOn("y shift = " << y_shift.eval() << endl);
    //    DebugOn("z shift = " << z_shift.eval() << endl);
    
    return(Reg);
}



shared_ptr<Model<double>> build_SOC_MIQCP(vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, vector<double>& rot_trans, bool convex, const vector<pair<pair<int,int>,pair<int,int>>>& incompatibles, param<>& norm_x, param<>& norm_y, param<>& norm_z,  param<>& intercept, const vector<int>& init_matching){
    double roll_1 = 0, yaw_1 = 0, pitch_1 = 0;
    int nb_pairs = 0, min_nb_pairs = numeric_limits<int>::max(), max_nb_pairs = 0, av_nb_pairs = 0;
    size_t nm = point_cloud_model.size(), nd = point_cloud_data.size();
    
    vector<double> zeros = {0,0,0};
    
    param<> x1("x1"), x2("x2"), y1("y1"), y2("y2"), z1("z1"), z2("z2");
    int m = av_nb_pairs;
    string i_str, j_str;
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
    
    
    indices Pairs("Pairs"), cells("cells");
    int idx1 = 0;
    int idx2 = 0;
    indices N1("N1"),N2("N2");
    Debug("nd = " << nd << endl);
    Debug("nm = " << nm << endl);
    
    N1 = range(1,nd);
    N2 = range(1,nm);
    cells = indices(N1,N2);
    string name="SOC_MIQCP";
    
    auto Reg=make_shared<Model<>>(name);
    
    double shift_min_x = 0.125, shift_max_x = 0.25, shift_min_y = -0.25,shift_max_y = -0.125,shift_min_z = -0.125,shift_max_z = 0;
    double yaw_min = -12.5*pi/180., yaw_max = 0, pitch_min = 12.5*pi/180.,pitch_max = 25.*pi/180.,roll_min = -12.5*pi/180.,roll_max = 0;
    //    double shift_min_x = -0.25, shift_max_x = 0.25, shift_min_y = -0.25,shift_max_y = 0.25,shift_min_z = -0.25,shift_max_z = 0.25;
    var<> x_shift("x_shift", shift_min_x, shift_max_x), y_shift("y_shift", shift_min_y, shift_max_y), z_shift("z_shift", shift_min_z, shift_max_z);
    //    var<> x_shift("x_shift", 0.23, 0.24), y_shift("y_shift", -0.24, -0.23), z_shift("z_shift", -0.011, -0.01);
    
    Reg->add(x_shift.in(R(1)),y_shift.in(R(1)),z_shift.in(R(1)));
    
    var<int> bin("bin",0,1);
    Reg->add(bin.in(cells));
    DebugOn("Added " << cells.size() << " binary variables" << endl);
    double angle_max = 25.*pi/180.;
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
    
    
    var<> theta11("theta11",  r11._range->first, r11._range->second), theta12("theta12", r12._range->first, r12._range->second), theta13("theta13", r13._range->first, r13._range->second);
    var<> theta21("theta21", r21._range->first, r21._range->second), theta22("theta22", r22._range->first, r22._range->second), theta23("theta23", r23._range->first, r23._range->second);
    var<> theta31("theta31", r31._range->first, r31._range->second), theta32("theta32", r32._range->first, r32._range->second), theta33("theta33", r33._range->first, r33._range->second);
    
    //            var<> theta11("theta11",  0.8, 1), theta12("theta12", -1, 1), theta13("theta13", -1, 1);
    //            var<> theta21("theta21",  -1, 1), theta22("theta22", -1, 1), theta23("theta23", -1, 1);
    //            var<> theta31("theta31",  -1, 1), theta32("theta32", -1, 1), theta33("theta33", 0.8, 1);
    
    //    var<> theta11("theta11",  0.96, 0.97), theta12("theta12", 0.15, 0.16), theta13("theta13", -0.22, -0.2);
    //    var<> theta21("theta21",  -0.21, -0.2), theta22("theta22", 0.95, 0.96), theta23("theta23", -0.24, 0.23);
    //    var<> theta31("theta31",  0.17, 0.18), theta32("theta32", 0.26, 0.27), theta33("theta33", 0.94, 0.95);
    Reg->add(theta11.in(R(1)),theta12.in(R(1)),theta13.in(R(1)));
    Reg->add(theta21.in(R(1)),theta22.in(R(1)),theta23.in(R(1)));
    Reg->add(theta31.in(R(1)),theta32.in(R(1)),theta33.in(R(1)));
    //    Reg->print();
    
    param<> x_new_lb("x_new_lb");
    x_new_lb.in(N1);
    param<> x_new_ub("x_new_ub");
    x_new_ub.in(N1);
    param<> y_new_lb("y_new_lb");
    y_new_lb.in(N1);
    param<> y_new_ub("y_new_ub");
    y_new_ub.in(N1);
    param<> z_new_lb("z_new_lb");
    z_new_lb.in(N1);
    param<> z_new_ub("z_new_ub");
    z_new_ub.in(N1);
    
    
    double x_lb = 0, y_lb = 0, z_lb = 0, x1_i = 0, y1_i = 0, z1_i = 0;
    shared_ptr<pair<double,double>> x1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> y1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> z1_bounds = make_shared<pair<double,double>>();
    for (int i = 0; i<nd; i++) {
        x1_bounds->first = x1.eval(i);
        x1_bounds->second = x1.eval(i);
        y1_bounds->first = y1.eval(i);
        y1_bounds->second = y1.eval(i);
        z1_bounds->first = z1.eval(i);
        z1_bounds->second = z1.eval(i);
        auto x_range  = get_product_range(x1_bounds, theta11._range);
        auto y_range  = get_product_range(y1_bounds, theta12._range);
        auto z_range  = get_product_range(z1_bounds, theta13._range);
        x_new_lb.set_val(i, x_range->first + y_range->first + z_range->first + x_shift.get_lb().eval());
        x_new_ub.set_val(i, x_range->second + y_range->second + z_range->second+ x_shift.get_ub().eval());
        x_range  = get_product_range(x1_bounds, theta21._range);
        y_range  = get_product_range(y1_bounds, theta22._range);
        z_range  = get_product_range(z1_bounds, theta23._range);
        y_new_lb.set_val(i, x_range->first + y_range->first + z_range->first + y_shift.get_lb().eval());
        y_new_ub.set_val(i, x_range->second + y_range->second + z_range->second+ y_shift.get_ub().eval());
        x_range  = get_product_range(x1_bounds, theta31._range);
        y_range  = get_product_range(y1_bounds, theta32._range);
        z_range  = get_product_range(z1_bounds, theta33._range);
        z_new_lb.set_val(i, x_range->first + y_range->first + z_range->first + z_shift.get_lb().eval());
        z_new_ub.set_val(i, x_range->second + y_range->second + z_range->second+ z_shift.get_ub().eval());
    }
    
    var<> new_xm("new_xm", -1, 1), new_ym("new_ym", -1, 1), new_zm("new_zm", -1, 1);
    var<> new_x1("new_x1", x_new_lb, x_new_ub), new_y1("new_y1", y_new_lb, y_new_ub), new_z1("new_z1", z_new_lb, z_new_ub);
    //            var<> new_x1("new_x1", -1, 1), new_y1("new_y1", -1, 1), new_z1("new_z1", -1, 1);
    Reg->add(new_xm.in(N1), new_ym.in(N1), new_zm.in(N1));
    Reg->add(new_x1.in(N1), new_y1.in(N1), new_z1.in(N1));
    
    //    var<> new_x1_ij("new_x1_ij", -1, 1), new_y1_ij("new_y1_ij", -1, 1), new_z1_ij("new_z1_ij", -1, 1);
    //    Reg->add(new_x1_ij.in(cells), new_y1_ij.in(cells), new_z1_ij.in(cells));
    
    
    var<> delta("delta", pos_);
    Reg->add(delta.in(N1));
    
    
    //    var<> d1("d1", 0,4),d2("d2", 0,4),d3("d3", 0,4),d4("d4", 0,4);
    //    var<> l12("l12", -2,2),l13("l13", -2,2),l14("l14", -2,2),l23("l23", -2,2),l24("l24", -2,2),l34("l34", -2,2);
    //    Reg->add(l12.in(R(1)),l13.in(R(1)),l14.in(R(1)));
    //    Reg->add(l23.in(R(1)),l24.in(R(1)),l34.in(R(1)));
    //    Reg->add(d1.in(R(1)),d2.in(R(1)),d3.in(R(1)),d4.in(R(1)));
    
    
    indices ids = indices("in_x");
    ids.add_empty_row();
    for(auto i=0;i<nd;i++){
        for(auto j=1;j<=nm;j++){
            if(cells.has_key(to_string(i+1)+","+to_string(j)))
                ids.add_in_row(i, to_string(j));
        }
    }
    
    bool add_shift_cut = false;
    
    if(add_shift_cut){
        Constraint<> Centroid("Centroid");
        Centroid = sum(new_x1) + sum(new_y1) + sum(new_z1) - (sum(new_xm) + sum(new_ym) + sum(new_zm));
        Reg->add(Centroid==0);
        
        Constraint<> Centroid_tx("Centroid_tx");
        Centroid_tx = nd*x_shift - sum(new_xm);
        Reg->add(Centroid_tx==0);
        
        Constraint<> Centroid_ty("Centroid_ty");
        Centroid_ty = nd*y_shift - sum(new_ym);
        Reg->add(Centroid_ty==0);
        
        Constraint<> Centroid_tz("Centroid_tz");
        Centroid_tz = nd*z_shift - sum(new_zm);
        Reg->add(Centroid_tz==0);
    }
    
    bool add_voronoi = true;
    
    if(add_voronoi){
        indices voronoi_ids("voronoi_ids");
        voronoi_ids = indices(range(1,3), *norm_x._indices);
        auto voronoi_ids_coefs = voronoi_ids.ignore_ith(0, 1);
        auto voronoi_ids_m = voronoi_ids.from_ith(0, 1);
        auto voronoi_ids_data = voronoi_ids.ignore_ith(1, 2);
        auto voronoi_ids_bin = voronoi_ids.ignore_ith(2, 1);
        //        Constraint<> Voronoi_model("Voronoi_model");
        //        Voronoi_model = norm_x.in(voronoi_ids_coefs)*new_xm.in(voronoi_ids_m) + norm_y.in(voronoi_ids_coefs)*new_ym.in(voronoi_ids_m) + norm_z.in(voronoi_ids_coefs)*new_zm.in(voronoi_ids_m) + intercept.in(voronoi_ids_coefs);
        //        Reg->add_on_off_multivariate_refined(Voronoi_model.in(voronoi_ids)<=0, bin.in(voronoi_ids_bin), true);
        
        
        Constraint<> Voronoi("Voronoi");
        Voronoi = norm_x.in(voronoi_ids_coefs)*new_x1.in(voronoi_ids_data) + norm_y.in(voronoi_ids_coefs)*new_y1.in(voronoi_ids_data) + norm_z.in(voronoi_ids_coefs)*new_z1.in(voronoi_ids_data) + intercept.in(voronoi_ids_coefs);
        Reg->add_on_off_multivariate_refined(Voronoi.in(voronoi_ids)<=0, bin.in(voronoi_ids_bin), true);
        //        Reg->print();
    }
    
    
    //    Reg->print();
    
    
    //
    theta11.initialize_all(1);
    theta22.initialize_all(1);
    theta33.initialize_all(1);
    
    bool spatial_branching = true;
    if(spatial_branching){
        /* Spatial branching vars */
        int nb_pieces = 5; // Divide each axis into nb_pieces
        indices spatial_ids("spatial_ids");
        spatial_ids = range(1,nb_pieces);
        indices theta_ids("theta_ids");
        theta_ids = indices(range(1,3),range(1,3));
        indices shift_ids("shift_ids");
        shift_ids.insert({"x", "y", "z"});
        indices theta_spatial("theta_spatial");
        theta_spatial = indices(theta_ids, spatial_ids);
        
        indices shift_spatial("shift_spatial");
        shift_spatial = indices(shift_ids, spatial_ids);
        
        //    var<int> sbin_shift("sbin_shift", 0, 1);
        //    var<int> sbin_theta("sbin_theta", 0, 1);
        
        //    Reg->add(sbin_theta.in(theta_spatial),sbin_shift.in(shift_spatial));
        
        var<int> sbin_tx("sbin_tx", 0, 1), sbin_ty("sbin_ty", 0, 1), sbin_tz("sbin_tz", 0, 1);
        var<int> sbin_theta11("sbin_theta11", 0, 1), sbin_theta12("sbin_theta12", 0, 1), sbin_theta13("sbin_theta13", 0, 1);
        var<int> sbin_theta21("sbin_theta21", 0, 1), sbin_theta22("sbin_theta22", 0, 1), sbin_theta23("sbin_theta23", 0, 1);
        var<int> sbin_theta31("sbin_theta31", 0, 1), sbin_theta32("sbin_theta32", 0, 1), sbin_theta33("sbin_theta33", 0, 1);
        Reg->add(sbin_theta11.in(spatial_ids),sbin_theta12.in(spatial_ids),sbin_theta13.in(spatial_ids));
        Reg->add(sbin_theta21.in(spatial_ids),sbin_theta22.in(spatial_ids),sbin_theta23.in(spatial_ids));
        Reg->add(sbin_theta31.in(spatial_ids),sbin_theta32.in(spatial_ids),sbin_theta33.in(spatial_ids));
        Reg->add(sbin_tx.in(spatial_ids),sbin_ty.in(spatial_ids),sbin_tz.in(spatial_ids));
        /* Spatial branching constraints */
        Constraint<> OneBinAngleSpatial11("OneBinAngleSpatial11");
        OneBinAngleSpatial11 = sum(sbin_theta11);
        Reg->add(OneBinAngleSpatial11==1);
        
        Constraint<> OneBinAngleSpatial12("OneBinAngleSpatial12");
        OneBinAngleSpatial12 = sum(sbin_theta12);
        Reg->add(OneBinAngleSpatial12==1);
        
        Constraint<> OneBinAngleSpatial13("OneBinAngleSpatial13");
        OneBinAngleSpatial13 = sum(sbin_theta13);
        Reg->add(OneBinAngleSpatial13==1);
        
        Constraint<> OneBinAngleSpatial21("OneBinAngleSpatial21");
        OneBinAngleSpatial21 = sum(sbin_theta21);
        Reg->add(OneBinAngleSpatial21==1);
        
        Constraint<> OneBinAngleSpatial22("OneBinAngleSpatial22");
        OneBinAngleSpatial22 = sum(sbin_theta22);
        Reg->add(OneBinAngleSpatial22==1);
        
        
        Constraint<> OneBinAngleSpatial23("OneBinAngleSpatial23");
        OneBinAngleSpatial23 = sum(sbin_theta23);
        Reg->add(OneBinAngleSpatial23==1);
        
        Constraint<> OneBinAngleSpatial31("OneBinAngleSpatial31");
        OneBinAngleSpatial31 = sum(sbin_theta31);
        Reg->add(OneBinAngleSpatial31==1);
        
        Constraint<> OneBinAngleSpatial32("OneBinAngleSpatial32");
        OneBinAngleSpatial32 = sum(sbin_theta32);
        Reg->add(OneBinAngleSpatial32==1);
        
        Constraint<> OneBinAngleSpatial33("OneBinAngleSpatial33");
        OneBinAngleSpatial33 = sum(sbin_theta33);
        Reg->add(OneBinAngleSpatial33==1);
        
        Constraint<> OneBinShiftSpatialx("OneBinShiftSpatialx");
        OneBinShiftSpatialx = sum(sbin_tx);
        Reg->add(OneBinShiftSpatialx==1);
        
        Constraint<> OneBinShiftSpatialy("OneBinShiftSpatialy");
        OneBinShiftSpatialy = sum(sbin_ty);
        Reg->add(OneBinShiftSpatialy==1);
        
        Constraint<> OneBinShiftSpatialz("OneBinShiftSpatialz");
        OneBinShiftSpatialz = sum(sbin_tz);
        Reg->add(OneBinShiftSpatialz==1);
        
        //    Reg->print();
        
        double diag_increment = 1./nb_pieces;/* Diagonals are defined in [0,1] */
        double off_diag_increment = 2./nb_pieces;/* Diagonals are defined in [-1,1] */
        double shift_increment = 0.5/nb_pieces;/* Shifts are defined in [-0.25,0.25] */
        
        auto spatial_ids_n = range(1,nb_pieces-1);
        auto spatial_ids_1 = range(2,nb_pieces);
        param<> diag_lb("diag_lb"), diag_ub("diag_ub"), off_diag_lb("off_diag_lb"), off_diag_ub("off_diag_ub");
        param<> t_lb("t_lb"), t_ub("t_ub");
        diag_ub.in(spatial_ids);
        diag_lb.in(spatial_ids);
        off_diag_ub.in(spatial_ids);
        off_diag_lb.in(spatial_ids);
        t_ub.in(spatial_ids);
        t_lb.in(spatial_ids);
        for (int i = 0; i<nb_pieces; i++) {
            diag_ub.set_val(i,(i+1)*diag_increment);
            off_diag_ub.set_val(i,-1+(i+1)*off_diag_increment);
            t_ub.set_val(i,-0.25 + (i+1)*diag_increment);
            diag_lb.set_val(i,i*diag_increment);
            off_diag_lb.set_val(i,-1 + i*off_diag_increment);
            t_lb.set_val(i,-0.25 + i*diag_increment);
        }
        auto ids_repeat = theta11.repeat_id(nb_pieces-1);
        
        Constraint<> Diag_Spatial_UB11("Diag_Spatial_UB11");
        Diag_Spatial_UB11 = theta11.in(ids_repeat) - diag_ub.in(spatial_ids_n);
        Reg->add_on_off_multivariate_refined(Diag_Spatial_UB11.in(spatial_ids_n)<=0, sbin_theta11.in(spatial_ids_n), true);
        
        Constraint<> Diag_Spatial_LB11("Diag_Spatial_LB11");
        Diag_Spatial_LB11 = diag_lb.in(spatial_ids_1) - theta11.in(ids_repeat);
        Reg->add_on_off_multivariate_refined(Diag_Spatial_LB11.in(spatial_ids_1) <= 0, sbin_theta11.in(spatial_ids_1), true);
        
        Constraint<> Diag_Spatial_UB22("Diag_Spatial_UB22");
        Diag_Spatial_UB22 = theta22.in(ids_repeat) - diag_ub.in(spatial_ids_n);
        Reg->add_on_off_multivariate_refined(Diag_Spatial_UB22.in(spatial_ids_n)<=0, sbin_theta22.in(spatial_ids_n), true);
        
        Constraint<> Diag_Spatial_LB33("Diag_Spatial_LB33");
        Diag_Spatial_LB33 = diag_lb.in(spatial_ids_1) - theta33.in(ids_repeat);
        Reg->add_on_off_multivariate_refined(Diag_Spatial_LB33.in(spatial_ids_1) <= 0, sbin_theta33.in(spatial_ids_1), true);
        
        
        Constraint<> Diag_Spatial_UB12("Diag_Spatial_UB12");
        Diag_Spatial_UB12 = theta12.in(ids_repeat) - off_diag_ub.in(spatial_ids_n);
        Reg->add_on_off_multivariate_refined(Diag_Spatial_UB12.in(spatial_ids_n) <= 0, sbin_theta12.in(spatial_ids_n), true);
        
        Constraint<> Diag_Spatial_LB12("Diag_Spatial_LB12");
        Diag_Spatial_LB12 = off_diag_lb.in(spatial_ids_1) - theta12.in(ids_repeat);
        Reg->add_on_off_multivariate_refined(Diag_Spatial_LB12.in(spatial_ids_1) <= 0, sbin_theta12.in(spatial_ids_1), true);
        
        Constraint<> Diag_Spatial_UB13("Diag_Spatial_UB13");
        Diag_Spatial_UB13 = theta13.in(ids_repeat) - off_diag_ub.in(spatial_ids_n);
        Reg->add_on_off_multivariate_refined(Diag_Spatial_UB13.in(spatial_ids_n) <= 0, sbin_theta13.in(spatial_ids_n), true);
        
        Constraint<> Diag_Spatial_LB13("Diag_Spatial_LB13");
        Diag_Spatial_LB13 = off_diag_lb.in(spatial_ids_1) - theta13.in(ids_repeat);
        Reg->add_on_off_multivariate_refined(Diag_Spatial_LB13.in(spatial_ids_1) <= 0, sbin_theta13.in(spatial_ids_1), true);
        
        Constraint<> Diag_Spatial_UB21("Diag_Spatial_UB21");
        Diag_Spatial_UB21 = theta21.in(ids_repeat) - off_diag_ub.in(spatial_ids_n);
        Reg->add_on_off_multivariate_refined(Diag_Spatial_UB21.in(spatial_ids_n) <= 0, sbin_theta21.in(spatial_ids_n), true);
        
        Constraint<> Diag_Spatial_LB21("Diag_Spatial_LB21");
        Diag_Spatial_LB21 = off_diag_lb.in(spatial_ids_1) - theta21.in(ids_repeat);
        Reg->add_on_off_multivariate_refined(Diag_Spatial_LB21.in(spatial_ids_1) <= 0, sbin_theta21.in(spatial_ids_1), true);
        
        Constraint<> Diag_Spatial_UB23("Diag_Spatial_UB23");
        Diag_Spatial_UB23 = theta23.in(ids_repeat) - off_diag_ub.in(spatial_ids_n);
        Reg->add_on_off_multivariate_refined(Diag_Spatial_UB23.in(spatial_ids_n) <= 0, sbin_theta23.in(spatial_ids_n), true);
        
        Constraint<> Diag_Spatial_LB23("Diag_Spatial_LB23");
        Diag_Spatial_LB23 = off_diag_lb.in(spatial_ids_1) - theta23.in(ids_repeat);
        Reg->add_on_off_multivariate_refined(Diag_Spatial_LB23.in(spatial_ids_1) <= 0, sbin_theta23.in(spatial_ids_1), true);
        
        Constraint<> Diag_Spatial_UB31("Diag_Spatial_UB31");
        Diag_Spatial_UB31 = theta31.in(ids_repeat) - off_diag_ub.in(spatial_ids_n);
        Reg->add_on_off_multivariate_refined(Diag_Spatial_UB31.in(spatial_ids_n) <= 0, sbin_theta31.in(spatial_ids_n), true);
        
        Constraint<> Diag_Spatial_LB31("Diag_Spatial_LB31");
        Diag_Spatial_LB31 = off_diag_lb.in(spatial_ids_1) - theta31.in(ids_repeat);
        Reg->add_on_off_multivariate_refined(Diag_Spatial_LB31.in(spatial_ids_1) <= 0, sbin_theta31.in(spatial_ids_1), true);
        
        Constraint<> Diag_Spatial_UB32("Diag_Spatial_UB32");
        Diag_Spatial_UB32 = theta32.in(ids_repeat) - off_diag_ub.in(spatial_ids_n);
        Reg->add_on_off_multivariate_refined(Diag_Spatial_UB32.in(spatial_ids_n) <= 0, sbin_theta32.in(spatial_ids_n), true);
        
        Constraint<> Diag_Spatial_LB32("Diag_Spatial_LB32");
        Diag_Spatial_LB32 = off_diag_lb.in(spatial_ids_1) - theta32.in(ids_repeat);
        Reg->add_on_off_multivariate_refined(Diag_Spatial_LB32.in(spatial_ids_1) <= 0, sbin_theta32.in(spatial_ids_1), true);
        
        Constraint<> xshift_Spatial_UB("xshift_Spatial_UB");
        xshift_Spatial_UB = x_shift.in(ids_repeat) - t_ub.in(spatial_ids_n);
        Reg->add_on_off_multivariate_refined(xshift_Spatial_UB.in(spatial_ids_n) <= 0, sbin_tx.in(spatial_ids_n), true);
        
        Constraint<> xshift_Spatial_LB("xshift_Spatial_LB");
        xshift_Spatial_LB = t_lb.in(spatial_ids_1) - x_shift.in(ids_repeat);
        Reg->add_on_off_multivariate_refined(xshift_Spatial_LB.in(spatial_ids_1) <= 0, sbin_tx.in(spatial_ids_1), true);
        
        Constraint<> yshift_Spatial_UB("yshift_Spatial_UB");
        yshift_Spatial_UB = y_shift.in(ids_repeat) - t_ub.in(spatial_ids_n);
        Reg->add_on_off_multivariate_refined(yshift_Spatial_UB.in(spatial_ids_n) <= 0, sbin_ty.in(spatial_ids_n), true);
        
        Constraint<> yshift_Spatial_LB("yshift_Spatial_LB");
        yshift_Spatial_LB = t_lb.in(spatial_ids_1) - y_shift.in(ids_repeat);
        Reg->add_on_off_multivariate_refined(yshift_Spatial_LB.in(spatial_ids_1) <= 0, sbin_ty.in(spatial_ids_1), true);
        
        Constraint<> zshift_Spatial_UB("zshift_Spatial_UB");
        zshift_Spatial_UB = z_shift.in(ids_repeat) - t_ub.in(spatial_ids_n);
        Reg->add_on_off_multivariate_refined(zshift_Spatial_UB.in(spatial_ids_n) <= 0, sbin_tz.in(spatial_ids_n), true);
        
        Constraint<> zshift_Spatial_LB("zshift_Spatial_LB");
        zshift_Spatial_LB = t_lb.in(spatial_ids_1) - z_shift.in(ids_repeat);
        Reg->add_on_off_multivariate_refined(zshift_Spatial_LB.in(spatial_ids_1) <= 0, sbin_tz.in(spatial_ids_1), true);
        //        Reg->print();
    }
    
    
    Constraint<> Def_newxm("Def_newxm");
    Def_newxm = new_xm-product(x2.in(ids),bin.in_matrix(1, 1));
    Reg->add(Def_newxm.in(N1)==0);
    
    Constraint<> Def_newym("Def_newym");
    Def_newym = new_ym-product(y2.in(ids),bin.in_matrix(1, 1));
    Reg->add(Def_newym.in(N1)==0);
    
    Constraint<> Def_newzm("Def_newzm");
    Def_newzm = new_zm-product(z2.in(ids),bin.in_matrix(1, 1));
    Reg->add(Def_newzm.in(N1)==0);
    
    Constraint<> OneBin("OneBin");
    OneBin = bin.in_matrix(1, 1);
    Reg->add(OneBin.in(N1)==1);
    
    //    Constraint<> OneBin2("OneBin2");
    //    OneBin2 = bin.in_matrix(0, 1);
    //    Reg->add(OneBin2.in(N2)<=1);
    
    if(!incompatibles.empty()){
        indices pairs1("pairs1"), pairs2("pairs2");
        pairs1 = cells;
        pairs2 = cells;
        for (const auto &inc_pair : incompatibles) {
            pairs1.add_ref(to_string(inc_pair.first.first+1)+","+to_string(inc_pair.second.first+1));
            pairs2.add_ref(to_string(inc_pair.first.second+1)+","+to_string(inc_pair.second.second+1));
        }
        
        Constraint<> incomp_pairs("incomp_pairs");
        incomp_pairs = bin.in(pairs1) + bin.in(pairs2);
        Reg->add(incomp_pairs.in(range(1,pairs1.size()))<=1);
        //        incomp_pairs.print();
    }
    
    
    
    double scaling = 1;
    Constraint<> Norm2_new("Norm2_new");
    Norm2_new -= delta - scaling*(pow(new_x1 - new_xm,2) + pow(new_y1 - new_ym,2) + pow(new_z1 - new_zm,2));
    Reg->add(Norm2_new.in(N1)<=0);
    
    bool add_delta_cut = false;
    if(add_delta_cut){
        var<> new_x1_sqr("new_x1_sqr", 0, max(pow(new_x1.get_lb(),2),pow(new_x1.get_ub(),2)));
        var<> new_y1_sqr("new_y1_sqr", 0, max(pow(new_y1.get_lb(),2),pow(new_y1.get_ub(),2)));
        var<> new_z1_sqr("new_z1_sqr", 0, max(pow(new_z1.get_lb(),2),pow(new_y1.get_ub(),2)));
        Reg->add(new_x1_sqr.in(N1),new_y1_sqr.in(N1),new_z1_sqr.in(N1));
        
        bool split = true, convexify = true;
        Constraint<> new_x1_Square("new_x1_Square");
        new_x1_Square -= new_x1_sqr - pow(new_x1,2);
        //        Reg->add(new_x1_Square.in(N1)==0,convexify,"on/off",split);/* Convexify and split nonconvex equation */
        Reg->add(new_x1_Square.in(N1)<=0);
        
        Constraint<> new_y1_Square("new_y1_Square");
        new_y1_Square -= new_y1_sqr - pow(new_y1,2);
        //        Reg->add(new_y1_Square.in(N1)==0,convexify,"on/off",split);/* Convexify and split nonconvex equation */
        Reg->add(new_y1_Square.in(N1)<=0);
        
        Constraint<> new_z1_Square("new_z1_Square");
        new_z1_Square -= new_z1_sqr - pow(new_z1,2);
        //        Reg->add(new_z1_Square.in(N1)==0,convexify,"on/off",split);/* Convexify and split nonconvex equation */
        Reg->add(new_z1_Square.in(N1)<=0);
        
        Constraint<> DeltaCut("DeltaCut");
        DeltaCut -= delta.from(cells);
        DeltaCut += scaling*(pow(x2.to(cells),2) + pow(y2.to(cells),2) + pow(z2.to(cells),2));
        DeltaCut -= scaling*(2*new_x1.from(cells)*x2.to(cells) + 2*new_y1.from(cells)*y2.to(cells) + 2*new_z1.from(cells)*z2.to(cells));
        DeltaCut += scaling*(new_x1_sqr.from(cells) + new_y1_sqr.from(cells) + new_z1_sqr.from(cells));
        Reg->add_on_off_multivariate_refined(DeltaCut.in(cells)<=0, bin, true);
    }
    bool add_delta_ij = false;
    
    if(add_delta_ij) {
        var<> delta_ij("delta_ij", 0, sqrt(8));
        Reg->add(delta_ij.in(cells));
        
        Constraint<> Norm2("Norm2");
        Norm2 -= delta_ij - pow(new_x1.from(cells) - x2.to(cells),2) - pow(new_y1.from(cells) - y2.to(cells),2) - pow(new_z1.from(cells) - z2.to(cells),2);
        Reg->add(Norm2.in(cells)<=0);
        
        
        Constraint<> DeltaMin("DeltaMin");
        DeltaMin -= delta;
        DeltaMin += scaling*bin.in_matrix(1, 1)*delta_ij.in_matrix(1, 1);
        Reg->add(DeltaMin.in(N1)<=0);
    }
    
    bool add_distance_cut = true;
    if(add_distance_cut){
        
    }
    //    Reg->print();
    
    bool add_sdp_rel = true;
    if(add_sdp_rel){
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
        Reg->add(soc_12.in(range(0,0))<=0);
        
        Constraint<> soc_13("soc_13");
        soc_13 = pow(theta12-theta21,2)-(1-theta11-theta22+theta33)*(1+theta11+theta22+theta33);
        Reg->add(soc_13.in(range(0,0))<=0);
        
        Constraint<> soc_14("soc_14");
        soc_14 = pow(theta23+theta32,2)-(1-theta11-theta22+theta33)*(1-theta11+theta22-theta33);
        Reg->add(soc_14.in(range(0,0))<=0);
        
        Constraint<> soc_23("soc_23");
        soc_23 = pow(theta23-theta32,2)-(1+theta11-theta22-theta33)*(1+theta11+theta22+theta33);
        Reg->add(soc_23.in(range(0,0))<=0);
        
        Constraint<> soc_24("soc_24");
        soc_24 = pow(theta12+theta21,2)-(1+theta11-theta22-theta33)*(1-theta11+theta22-theta33);
        Reg->add(soc_24.in(range(0,0))<=0);
        
        Constraint<> soc_34("soc_34");
        soc_34 = pow(theta31-theta13,2)-(1+theta11+theta22+theta33)*(1-theta11+theta22-theta33);
        Reg->add(soc_34.in(range(0,0))<=0);
        
        Constraint<> det_123("det_123");
        det_123+=(theta13+theta31)*((theta13+theta31)*(1+theta11+theta22+theta33)-(theta23-theta32)*(theta12-theta21));
        det_123-=(1-theta11-theta22+theta33)*((1+theta11-theta22-theta33)*(1+theta11+theta22+theta33)-pow(theta23-theta32,2));
        det_123-=(theta12-theta21)*((theta13+theta31)*(theta23-theta32)-(theta12-theta21)*(1+theta11-theta22-theta33));
        Reg->add(det_123.in(range(0,0))<=0);
        
        Constraint<> det_124("det_124");
        det_124+=(theta13+theta31)*((theta13+theta31)*(1-theta11+theta22-theta33)-(theta23+theta32)*(theta12+theta21));
        det_124-=(1-theta11-theta22+theta33)*((1+theta11-theta22-theta33)*(1-theta11+theta22-theta33)-pow(theta12+theta21,2));
        det_124-=(theta23+theta32)*((theta13+theta31)*(theta12+theta21)-(theta23+theta32)*(1+theta11-theta22-theta33));
        Reg->add(det_124.in(range(0,0))<=0);
        
        Constraint<> det_134("det_134");
        det_134+=(theta12-theta21)*((theta12-theta21)*(1-theta11+theta22-theta33)-(theta23+theta32)*(theta31-theta13));
        det_134-=(1-theta11-theta22+theta33)*((1+theta11+theta22+theta33)*(1-theta11+theta22-theta33)-pow(theta31-theta13,2));
        det_134-=(theta23+theta32)*((theta12-theta21)*(theta31-theta13)-(theta23+theta32)*(1+theta11+theta22+theta33));
        Reg->add(det_134.in(range(0,0))<=0);
        
        Constraint<> det_234("det_234");
        det_234+=(theta23-theta32)*((theta23-theta32)*(1-theta11+theta22-theta33)-(theta12+theta21)*(theta31-theta13));
        det_234-=(1+theta11-theta22-theta33)*((1+theta11+theta22+theta33)*(1-theta11+theta22-theta33)-pow(theta31-theta13,2));
        det_234-=(theta12+theta21)*((theta23-theta32)*(theta31-theta13)-(theta12+theta21)*(1+theta11+theta22+theta33));
        Reg->add(det_234.in(range(0,0))<=0);
        if(convex){
            Constraint<> row1("row1");
            row1 = pow(theta11,2)+pow(theta12,2)+pow(theta13,2);
            Reg->add(row1.in(range(0,0))<=1);
            Constraint<> row2("row2");
            row2 = pow(theta21,2)+pow(theta22,2)+pow(theta23,2);
            Reg->add(row2.in(range(0,0))<=1);
            Constraint<> row3("row3");
            row3 = pow(theta31,2)+pow(theta32,2)+pow(theta33,2);
            Reg->add(row3.in(range(0,0))<=1);
            Constraint<> col1("col1");
            col1 = pow(theta11,2)+pow(theta21,2)+pow(theta31,2);
            Reg->add(col1.in(range(0,0))<=1);
            Constraint<> col2("col2");
            col2 = pow(theta12,2)+pow(theta22,2)+pow(theta32,2);
            Reg->add(col2.in(range(0,0))<=1);
            Constraint<> col3("col3");
            col3 = pow(theta13,2)+pow(theta23,2)+pow(theta33,2);
            Reg->add(col3.in(range(0,0))<=1);
        }
        else {
            Constraint<> row1("row1");
            row1 = pow(theta11,2)+pow(theta12,2)+pow(theta13,2);
            Reg->add(row1==1);
            Constraint<> row2("row2");
            row2 = pow(theta21,2)+pow(theta22,2)+pow(theta23,2);
            Reg->add(row2==1);
            Constraint<> row3("row3");
            row3 = pow(theta31,2)+pow(theta32,2)+pow(theta33,2);
            Reg->add(row3==1);
            Constraint<> col1("col1");
            col1 = pow(theta11,2)+pow(theta21,2)+pow(theta31,2);
            Reg->add(col1==1);
            Constraint<> col2("col2");
            col2 = pow(theta12,2)+pow(theta22,2)+pow(theta32,2);
            Reg->add(col2==1);
            Constraint<> col3("col3");
            col3 = pow(theta13,2)+pow(theta23,2)+pow(theta33,2);
            Reg->add(col3==1);
        }
    }
    auto ids1 = theta11.repeat_id(cells.size());
    Constraint<> x_rot1("x_rot1");
    x_rot1 += new_x1 -x_shift;
    x_rot1 -= x1.in(N1)*theta11.in(ids1) + y1.in(N1)*theta12.in(ids1) + z1.in(N1)*theta13.in(ids1);
    Reg->add(x_rot1.in(N1)==0);
    
    Constraint<> y_rot1("y_rot1");
    y_rot1 += new_y1 - y_shift;
    y_rot1 -= x1.in(N1)*theta21.in(ids1) + y1.in(N1)*theta22.in(ids1) + z1.in(N1)*theta23.in(ids1);
    Reg->add(y_rot1.in(N1)==0);
    
    Constraint<> z_rot1("z_rot1");
    z_rot1 += new_z1 -z_shift;
    z_rot1 -= x1.in(N1)*theta31.in(ids1) + y1.in(N1)*theta32.in(ids1) + z1.in(N1)*theta33.in(ids1);
    Reg->add(z_rot1.in(N1)==0);
    
    
    
    for (int i = 1; i<=nd; i++) {
        string key = to_string(i)+","+to_string(init_matching.at(i-1)+1);
        bin._val->at(bin._indices->_keys_map->at(key)) = 1;
    }
    
    /* Objective function */
    
    Reg->min(sum(delta));
    //    bin._val->at(bin._indices->_keys_map->at("1,1")) = 1;
    //    bin.set_lb("1,1",1);
    //    bin.set_lb("2,2",1);
    //    bin.set_lb("3,3",1);
    //    for (int i = 0; i<nd; i++) {
    //        string key = to_string(i+1)+","+to_string(i+1);
    //        bin._val->at(bin._indices->_keys_map->at(key)) = 1;
    //    }
    //    Reg->print();
    //    solver<> S1(Reg,ipopt);
    for(int i = 1; i<= nd; i++){
        bin.param<int>::set_val(to_string(i)+","+to_string(i), 1);
    }
    solver<> S(Reg,gurobi);
    S.run();
    //    Reg->print();
    Reg->print_int_solution();
    
    Reg->print_solution();
    //    {
    //        indices voronoi_ids("voronoi_ids");
    //        voronoi_ids = indices(N1, *norm_x._indices);
    //        auto voronoi_ids_coefs = voronoi_ids.ignore_ith(0, 1);
    ////        auto voronoi_ids_m = voronoi_ids.from_ith(0, 1);
    //        auto voronoi_ids_data = voronoi_ids.ignore_ith(1, 2);
    //        auto voronoi_ids_bin = voronoi_ids.ignore_ith(2, 1);
    ////        Constraint<> Voronoi_model("Voronoi_model");
    ////        Voronoi_model = norm_x.in(voronoi_ids_coefs)*new_xm.in(voronoi_ids_m) + norm_y.in(voronoi_ids_coefs)*new_ym.in(voronoi_ids_m) + norm_z.in(voronoi_ids_coefs)*new_zm.in(voronoi_ids_m) + intercept.in(voronoi_ids_coefs);
    ////        Reg->add_on_off_multivariate_refined(Voronoi_model.in(voronoi_ids)<=0, bin.in(voronoi_ids_bin), true);
    //
    //
    //        auto func = bin.in(voronoi_ids_bin)*(norm_x.in(voronoi_ids_coefs)*new_x1.in(voronoi_ids_data) + norm_y.in(voronoi_ids_coefs)*new_y1.in(voronoi_ids_data) + norm_z.in(voronoi_ids_coefs)*new_z1.in(voronoi_ids_data) + intercept.in(voronoi_ids_coefs));
    //
    //        func.allocate_mem();
    //        func.eval_all();
    //        for (int i = 0; i<func.get_nb_inst(); i++) {
    //            if(func._val->at(i)>1e-6){
    //                DebugOn("instance " <<  i << " is violated \n");
    //                func.print(i,10);
    //                DebugOn(" | violation = " <<  func._val->at(i) << endl);
    //            }
    //        }
    //
    //    }
    DebugOn("row 1 " << pow(theta11.eval(),2)+pow(theta12.eval(),2)+pow(theta13.eval(),2)
            << endl);
    DebugOn("row 2 " << pow(theta21.eval(),2)+pow(theta22.eval(),2)+pow(theta23.eval(),2)
            << endl);
    DebugOn("row 3 " << pow(theta31.eval(),2)+pow(theta32.eval(),2)+pow(theta33.eval(),2)
            << endl);
    DebugOn("col 1 " << pow(theta11.eval(),2)+pow(theta21.eval(),2)+pow(theta31.eval(),2)
            << endl);
    DebugOn("col 2 " << pow(theta12.eval(),2)+pow(theta22.eval(),2)+pow(theta32.eval(),2)
            << endl);
    DebugOn("col 3 " << pow(theta13.eval(),2)+pow(theta23.eval(),2)+pow(theta33.eval(),2)
            << endl);
    
    DebugOn("row 12 " << (theta11.eval()*theta21.eval())+(theta12.eval()*theta22.eval())+(theta13.eval()*theta23.eval())
            << endl);
    DebugOn("row 13 " << (theta11.eval()*theta31.eval())+(theta12.eval()*theta32.eval())+(theta13.eval()*theta33.eval())
            << endl);
    DebugOn("row 23 " << (theta21.eval()*theta31.eval())+(theta22.eval()*theta32.eval())+(theta23.eval()*theta33.eval())
            << endl);
    DebugOn("Theta matrix = " << endl);
    DebugOn("|" << theta11.eval() << " " << theta12.eval() << " " << theta13.eval() << "|" << endl);
    DebugOn("|" << theta21.eval() << " " << theta22.eval() << " " << theta23.eval() << "|" << endl);
    DebugOn("|" << theta31.eval() << " " << theta32.eval() << " " << theta33.eval() << "|" << endl);
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
    auto pitch_val = std::atan2(theta32.eval(), theta33.eval())*180/pi;
    auto roll_val = std::atan2(-1*theta31.eval(), std::sqrt(theta32.eval()*theta32.eval()+theta33.eval()*theta33.eval()))*180/pi;
    auto yaw_val = std::atan2(theta21.eval(),theta11.eval())*180/pi;
    DebugOn("Roll (degrees) = " << to_string_with_precision(roll_val,12) << endl);
    DebugOn("Pitch (degrees) = " << to_string_with_precision(pitch_val,12) << endl);
    DebugOn("Yaw (degrees) = " << to_string_with_precision(yaw_val,12) << endl);
    DebugOn("x shift = " << x_shift.eval() << endl);
    DebugOn("y shift = " << y_shift.eval() << endl);
    DebugOn("z shift = " << z_shift.eval() << endl);
    
    return(Reg);
}

shared_ptr<Model<double>> build_TU_MIP(vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, vector<double>& rot_trans, const vector<pair<pair<int,int>,pair<int,int>>>& incompatibles){
    double angle_max = 1;
    double roll_1 = 0, yaw_1 = 0, pitch_1 = 0;
    int nb_pairs = 0, min_nb_pairs = numeric_limits<int>::max(), max_nb_pairs = 0, av_nb_pairs = 0;
    size_t nm = point_cloud_model.size(), nd = point_cloud_data.size();
    
    vector<double> zeros = {0,0,0};
    
    param<> x1("x1"), x2("x2"), y1("y1"), y2("y2"), z1("z1"), z2("z2");
    int m = av_nb_pairs;
    string i_str, j_str;
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
    
    
    indices Pairs("Pairs"), cells("cells");
    int idx1 = 0;
    int idx2 = 0;
    indices N1("N1"),N2("N2");
    Debug("nd = " << nd << endl);
    Debug("nm = " << nm << endl);
    
    N1 = range(1,nd);
    N2 = range(1,nm);
    cells = indices(N1,N2);
    
    
    
    string name="TU_MIP";
    
    auto Reg=make_shared<Model<>>(name);
    
    
    double shift_min_x = -1, shift_max_x = 1, shift_min_y = -1,shift_max_y = 1,shift_min_z = -1,shift_max_z = 1;
    var<> x_shift("x_shift", shift_min_x, shift_max_x), y_shift("y_shift", shift_min_y, shift_max_y), z_shift("z_shift", shift_min_z, shift_max_z);
    
    Reg->add(x_shift.in(R(1)),y_shift.in(R(1)),z_shift.in(R(1)));
    
    var<int> bin("bin",0,1);
    Reg->add(bin.in(cells));
    DebugOn("Added " << cells.size() << " binary variables" << endl);
    
    
    var<> theta11("theta11",  0, 1), theta12("theta12", -1, 1), theta13("theta13", -1, 1);
    var<> theta21("theta21",  -1, 1), theta22("theta22", -1, 1), theta23("theta23", -1, 1);
    var<> theta31("theta31",  -1, 1), theta32("theta32", -1, 1), theta33("theta33", 0, 1);
    Reg->add(theta11.in(R(1)),theta12.in(R(1)),theta13.in(R(1)));
    Reg->add(theta21.in(R(1)),theta22.in(R(1)),theta23.in(R(1)));
    Reg->add(theta31.in(R(1)),theta32.in(R(1)),theta33.in(R(1)));
    
    indices ids = indices("in_x");
    ids.add_empty_row();
    
    
    Constraint<> OneBin("OneBin");
    OneBin = bin.in_matrix(1, 1);
    Reg->add(OneBin.in(N1)==1);
    
    Constraint<> OneBin2("OneBin2");
    OneBin2 = bin.in_matrix(0, 1);
    Reg->add(OneBin2.in(N2)<=1);
    
    indices pairs1("pairs1"), pairs2("pairs2");
    pairs1 = cells;
    pairs2 = cells;
    for (const auto &inc_pair : incompatibles) {
        pairs1.add_ref(to_string(inc_pair.first.first+1)+","+to_string(inc_pair.second.first+1));
        pairs2.add_ref(to_string(inc_pair.first.second+1)+","+to_string(inc_pair.second.second+1));
    }
    
    Constraint<> incomp_pairs("incomp_pairs");
    incomp_pairs = bin.in(pairs1) + bin.in(pairs2);
    Reg->add(incomp_pairs.in(range(1,pairs1.size()))<=1);
    
    Constraint<> row1("row1");
    row1 = pow(theta11,2)+pow(theta12,2)+pow(theta13,2);
    Reg->add(row1<=1);
    Constraint<> row2("row2");
    row2 = pow(theta21,2)+pow(theta22,2)+pow(theta23,2);
    Reg->add(row2<=1);
    Constraint<> row3("row3");
    row3 = pow(theta31,2)+pow(theta32,2)+pow(theta33,2);
    Reg->add(row3<=1);
    Constraint<> col1("col1");
    col1 = pow(theta11,2)+pow(theta21,2)+pow(theta31,2);
    Reg->add(col1<=1);
    Constraint<> col2("col2");
    col2 = pow(theta12,2)+pow(theta22,2)+pow(theta32,2);
    Reg->add(col2<=1);
    Constraint<> col3("col3");
    col3 = pow(theta13,2)+pow(theta23,2)+pow(theta33,2);
    Reg->add(col3<=1);
    //    Constraint<> row1("row1");
    //    row1 = pow(theta11,2)+pow(theta12,2)+pow(theta13,2);
    //    Reg->add(row1<=1);
    //    Constraint<> row1("row1");
    //    row1 = pow(theta11,2)+pow(theta12,2)+pow(theta13,2);
    //    Reg->add(row1<=1);
    
    
    
    /* Objective function */
    func<> obj = nd*(x_shift*x_shift + y_shift*y_shift + z_shift*z_shift);
    obj += sum(x1*x1 + y1*y1 + z1*z1);
    auto ids_repeat = x_shift.repeat_id(cells.size());
    obj -= 2*sum(x2.to(cells)*x_shift.in(ids_repeat)*bin) + 2*sum(y2.to(cells)*y_shift.in(ids_repeat)*bin) + 2*sum(z2.to(cells)*z_shift.in(ids_repeat)*bin);
    
    obj += sum(x2.to(cells)*x2.to(cells)*bin) + sum(y2.to(cells)*y2.to(cells)*bin) + sum(z2.to(cells)*z2.to(cells)*bin);
    
    auto ids1 = theta11.repeat_id(cells.size());
    obj -= 2*sum(x2.to(cells)*x1.from(cells)*bin*theta11.in(ids1));
    obj -= 2*sum(x2.to(cells)*y1.from(cells)*bin*theta12.in(ids1));
    obj -= 2*sum(x2.to(cells)*z1.from(cells)*bin*theta13.in(ids1));
    obj -= 2*sum(y2.to(cells)*x1.from(cells)*bin*theta21.in(ids1));
    obj -= 2*sum(y2.to(cells)*y1.from(cells)*bin*theta22.in(ids1));
    obj -= 2*sum(y2.to(cells)*z1.from(cells)*bin*theta23.in(ids1));
    obj -= 2*sum(z2.to(cells)*x1.from(cells)*bin*theta31.in(ids1));
    obj -= 2*sum(z2.to(cells)*y1.from(cells)*bin*theta32.in(ids1));
    obj -= 2*sum(z2.to(cells)*z1.from(cells)*bin*theta33.in(ids1));
    
    
    //    func<> obj = sum(x1*x1 + y1*y1 + z1*z1);
    //    for (auto i = 0; i<nd; i++) {
    ////        obj += x1.eval(i)*x1.eval(i) + y1.eval(i)*y1.eval(i) + z1.eval(i)*z1.eval(i);
    //        for (auto j = 0; j<nm; j++) {
    //            string ij = to_string(i+1) +","+to_string(j+1);
    //            obj -= /*(x2.eval(j)*x2.eval(j) + y2.eval(j)*y2.eval(j) + z2.eval(j)*z2.eval(j) -*/ /*2*(x_shift[0]*x2.eval(j) + y_shift[0]*y2.eval(j) + z_shift[0]*z2.eval(j)) -*/ bin(ij)*(2*x2.eval(j)*(x1.eval(i)*theta11[0] + y1.eval(i)*theta12[0] + z1.eval(i)*theta13[0]) + 2*y2.eval(j)*(x1.eval(i)*theta21[0] + y1.eval(i)*theta22[0] + z1.eval(i)*theta23[0]) + 2*z2.eval(j)*(x1.eval(i)*theta31[0] + y1.eval(i)*theta32[0] + z1.eval(i)*theta33[0]) );
    //        }
    //    }
    
    //    auto ids1 = theta11.repeat_id(cells.size());
    //    obj += x1.in(N1)*theta11.in(ids1) + y1.in(N1)*theta12.in(ids1) + z1.in(N1)*theta13.in(ids1);
    //    obj.print();
    Reg->min(obj);
    
    //    Reg->print();
    //    solver<> S(Reg,ipopt);
    solver<> S(Reg,gurobi);
    S.run();
    Reg->print_int_solution();
    
    //    Reg->print_solution();
    
    DebugOn("Theta matrix = " << endl);
    DebugOn("|" << theta11.eval() << " " << theta12.eval() << " " << theta13.eval() << "|" << endl);
    DebugOn("|" << theta21.eval() << " " << theta22.eval() << " " << theta23.eval() << "|" << endl);
    DebugOn("|" << theta31.eval() << " " << theta32.eval() << " " << theta33.eval() << "|" << endl);
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
    DebugOn("x shift = " << x_shift.eval() << endl);
    DebugOn("y shift = " << y_shift.eval() << endl);
    DebugOn("z shift = " << z_shift.eval() << endl);
    
    return(Reg);
}

void run_ICR(vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, vector<double>& rot_trans, bool separate, const vector<pair<pair<int,int>,pair<int,int>>>& incompatibles, param<>& norm_x, param<>& norm_y, param<>& norm_z, param<>& intercept, const vector<pair<double,double>>& min_max_model, const vector<int>& matching){
    //    auto Reg = build_linobj_convex(point_cloud_model, point_cloud_data, rot_trans, separate, incompatibles, norm_x, norm_y, norm_z, intercept, min_max_model, matching);
    //    double obj = Reg->get_obj_val();
    //    double obj_init = obj;
    bool progress = true;
    while(progress){
        //        auto x2 = Reg->get_p
    }
    
}


shared_ptr<Model<double>> build_linobj_convex(vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, const indices& valid_cells, double new_roll_min, double new_roll_max, double new_pitch_min, double new_pitch_max, double new_yaw_min, double new_yaw_max, double new_shift_min_x, double new_shift_max_x, double new_shift_min_y, double new_shift_max_y, double new_shift_min_z, double new_shift_max_z, vector<double>& rot_trans, bool separate, const vector<pair<pair<int,int>,pair<int,int>>>& incompatibles,  param<>& norm_x,  param<>& norm_y,  param<>& norm_z,  param<>& intercept,const vector<int>& init_matching, const vector<double>& error_per_point, bool relax_inits, double ub){
    
    //    double shift_min_x = -0.25, shift_max_x = 0.25, shift_min_y = -0.25,shift_max_y = 0.25,shift_min_z = -0.25,shift_max_z = 0.25;
    //    double yaw_min = -25*pi/180., yaw_max = 25*pi/180., pitch_min = -25*pi/180.,pitch_max = 25.*pi/180.,roll_min = -25*pi/180.,roll_max = 25*pi/180.;
    //
    //
    //    double roll_1 = 0, yaw_1 = 0, pitch_1 = 0;
    int nb_pairs = 0, min_nb_pairs = numeric_limits<int>::max(), max_nb_pairs = 0, av_nb_pairs = 0;
    size_t nm = point_cloud_model.size(), nd = point_cloud_data.size();
    vector<double> zeros = {0,0,0};
    param<> x1("x1"), x2("x2"), y1("y1"), y2("y2"), z1("z1"), z2("z2");
    int m = av_nb_pairs;
    string i_str, j_str;
    double xm_max = numeric_limits<double>::lowest(), ym_max = numeric_limits<double>::lowest(), zm_max = numeric_limits<double>::lowest();
    double xm_min = numeric_limits<double>::max(), ym_min = numeric_limits<double>::max(), zm_min = numeric_limits<double>::max();
    
    
    for (auto i = 0; i<nd; i++) {
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
    Debug("nd = " << nd << endl);
    Debug("nm = " << nm << endl);
    
    N1 = range(1,nd);
    N2 = range(1,nm);
    // auto cells = indices(N1,N2);
    auto cells=valid_cells;
    string name="TU_MIP";
    
    auto Reg=make_shared<Model<>>(name);
    
    var<int> bin("bin",0,1);
    Reg->add(bin.in(cells));
    
    var<> x_shift("x_shift", new_shift_min_x, new_shift_max_x), y_shift("y_shift", new_shift_min_y, new_shift_max_y), z_shift("z_shift", new_shift_min_z, new_shift_max_z);
    
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
    
   
    
    
    DebugOn("Added " << cells.size() << " binary variables" << endl);
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
  
    
    bool hybrid=true;
    
    DebugOn("Added " << cells.size() << " binary variables" << endl);
    
    
    double tx_max=std::max(pow(new_shift_min_x,2), pow(new_shift_max_x,2));
    double ty_max=std::max(pow(new_shift_min_y,2), pow(new_shift_max_y,2));
    double tz_max=std::max(pow(new_shift_min_z,2), pow(new_shift_max_z,2));
    double tx_min, ty_min, tz_min;
    double T11_min, T11_max, T12_min,T12_max, T13_min, T13_max;
    double T21_min, T21_max, T22_min,T22_max, T23_min, T23_max;
    double T31_min, T31_max, T32_min,T32_max, T33_min, T33_max, Tmin, Tmax;
    if(new_shift_min_x<=0 && new_shift_max_x>=0)
        tx_min=0;
    else
        tx_min=std::min(pow(new_shift_min_x,2), pow(new_shift_max_x,2));
    if(new_shift_min_y<=0 && new_shift_max_y>=0)
        ty_min=0;
    else
        ty_min=std::min(pow(new_shift_min_y,2), pow(new_shift_max_y,2));
    if(new_shift_min_z<=0 && new_shift_max_z>=0)
        tz_min=0;
    else
        tz_min=std::min(pow(new_shift_min_z,2), pow(new_shift_max_z,2));
    
    Tmin=theta11.get_lb("0");
    Tmax=theta11.get_ub("0");
    if(Tmin<=0 && Tmax>=0)
        T11_min=0;
    else
        T11_min=std::min(pow(Tmin,2), pow(Tmax,2));
    T11_max=std::max(pow(Tmin,2), pow(Tmax,2));
    Tmin=theta12.get_lb("0");
    Tmax=theta12.get_ub("0");
    if(Tmin<=0 && Tmax>=0)
        T12_min=0;
    else
        T12_min=std::min(pow(Tmin,2), pow(Tmax,2));
    T12_max=std::max(pow(Tmin,2), pow(Tmax,2));
    Tmin=theta13.get_lb("0");
    Tmax=theta13.get_ub("0");
    if(Tmin<=0 && Tmax>=0)
        T13_min=0;
    else
        T13_min=std::min(pow(Tmin,2), pow(Tmax,2));
    T13_max=std::max(pow(Tmin,2), pow(Tmax,2));
    Tmin=theta21.get_lb("0");
    Tmax=theta21.get_ub("0");
    if(Tmin<=0 && Tmax>=0)
        T21_min=0;
    else
        T21_min=std::min(pow(Tmin,2), pow(Tmax,2));
    T21_max=std::max(pow(Tmin,2), pow(Tmax,2));
    Tmin=theta22.get_lb("0");
    Tmax=theta22.get_ub("0");
    if(Tmin<=0 && Tmax>=0)
        T22_min=0;
    else
        T22_min=std::min(pow(Tmin,2), pow(Tmax,2));
    T22_max=std::max(pow(Tmin,2), pow(Tmax,2));
    Tmin=theta23.get_lb("0");
    Tmax=theta23.get_ub("0");
    if(Tmin<=0 && Tmax>=0)
        T23_min=0;
    else
        T23_min=std::min(pow(Tmin,2), pow(Tmax,2));
    T23_max=std::max(pow(Tmin,2), pow(Tmax,2));
    Tmin=theta31.get_lb("0");
    Tmax=theta31.get_ub("0");
    if(Tmin<=0 && Tmax>=0)
        T31_min=0;
    else
        T31_min=std::min(pow(Tmin,2), pow(Tmax,2));
    T31_max=std::max(pow(Tmin,2), pow(Tmax,2));
    Tmin=theta32.get_lb("0");
    Tmax=theta32.get_ub("0");
    if(Tmin<=0 && Tmax>=0)
        T32_min=0;
    else
        T32_min=std::min(pow(Tmin,2), pow(Tmax,2));
    T32_max=std::max(pow(Tmin,2), pow(Tmax,2));
    Tmin=theta33.get_lb("0");
    Tmax=theta33.get_ub("0");
    if(Tmin<=0 && Tmax>=0)
        T33_min=0;
    else
        T33_min=std::min(pow(Tmin,2), pow(Tmax,2));
    T33_max=std::max(pow(Tmin,2), pow(Tmax,2));
    
//    var<> tx("tx", tx_min, tx_max), ty("ty", ty_min, ty_max), tz("tz", tz_min, tz_max);
//    var<> T11("T11", T11_min, T11_max), T12("T12", T12_min, T12_max), T13("T13", T13_min, T13_max);
//    var<> T21("T21", T21_min, T21_max), T22("T22", T22_min, T22_max), T23("T23", T23_min, T23_max);
//    var<> T31("T31", T31_min, T31_max), T32("T32", T32_min, T32_max), T33("T33", T33_min, T33_max);

    var<> tx("tx"), ty("ty"), tz("tz");
    var<> T11("T11"), T12("T12"), T13("T13");
    var<> T21("T21"), T22("T22"), T23("T23");
    var<> T31("T31"), T32("T32"), T33("T33");
    
    Reg->add(tx.in(R(1)),ty.in(R(1)),tz.in(R(1)));
    Reg->add(x_shift.in(R(1)),y_shift.in(R(1)),z_shift.in(R(1)));
    Reg->add(T11.in(R(1)),T12.in(R(1)),T13.in(R(1)));
    Reg->add(T21.in(R(1)),T22.in(R(1)),T23.in(R(1)));
    Reg->add(T31.in(R(1)),T32.in(R(1)),T33.in(R(1)));
    var<> new_xm("new_xm"), new_ym("new_ym"), new_zm("new_zm");
    var<> delta("delta");
    
    theta11.initialize_all(1);
    theta22.initialize_all(1);
    theta33.initialize_all(1);
    
    
    
    
    
    //    indices ids = indices("in_x");
    //    ids.add_empty_row();
    
    
    Constraint<> OneBin("OneBin");
    OneBin = bin.in_matrix(1, 1);
    Reg->add(OneBin.in(N1)==1);
    
    
    bool add_voronoi = false;
    if(add_voronoi){
        indices voronoi_ids("voronoi_ids");
        voronoi_ids = indices(range(1, 3), *norm_x._indices);
        auto voronoi_ids_coefs = voronoi_ids.ignore_ith(0, 1);
        //        auto voronoi_ids_m = voronoi_ids.from_ith(0, 1);
        auto voronoi_ids_data = voronoi_ids.ignore_ith(1, 2);
        auto voronoi_ids_bin = voronoi_ids.ignore_ith(2, 1);
        //        Constraint<> Voronoi_model("Voronoi_model");
        //        Voronoi_model = norm_x.in(voronoi_ids_coefs)*new_xm.in(voronoi_ids_m) + norm_y.in(voronoi_ids_coefs)*new_ym.in(voronoi_ids_m) + norm_z.in(voronoi_ids_coefs)*new_zm.in(voronoi_ids_m) + intercept.in(voronoi_ids_coefs);
        //        Reg->add_on_off_multivariate_refined(Voronoi_model.in(voronoi_ids)<=0, bin.in(voronoi_ids_bin), true);
        
        auto ids1 = theta11.repeat_id(voronoi_ids.size());
        Constraint<> Voronoi("Voronoi");
        Voronoi = norm_x.in(voronoi_ids_coefs)*(x1.in(voronoi_ids_data)*theta11.in(ids1) + y1.in(voronoi_ids_data)*theta12.in(ids1) + z1.in(voronoi_ids_data)*theta13.in(ids1)+x_shift.in(ids1)) + norm_y.in(voronoi_ids_coefs)*(x1.in(voronoi_ids_data)*theta21.in(ids1) + y1.in(voronoi_ids_data)*theta22.in(ids1) + z1.in(voronoi_ids_data)*theta23.in(ids1)+y_shift.in(ids1)) + norm_z.in(voronoi_ids_coefs)*(x1.in(voronoi_ids_data)*theta31.in(ids1) + y1.in(voronoi_ids_data)*theta32.in(ids1) + z1.in(voronoi_ids_data)*theta33.in(ids1)+z_shift.in(ids1)) + intercept.in(voronoi_ids_coefs);
        Reg->add_on_off_multivariate_refined(Voronoi.in(voronoi_ids)<=0, bin.in(voronoi_ids_bin), true);
        
    }
    
    if(false && !incompatibles.empty()){
        indices pairs1("pairs1"), pairs2("pairs2");
        pairs1 = cells;
        pairs2 = cells;
        for (const auto &inc_pair : incompatibles) {
            pairs1.add_ref(to_string(inc_pair.first.first+1)+","+to_string(inc_pair.second.first+1));
            pairs2.add_ref(to_string(inc_pair.first.second+1)+","+to_string(inc_pair.second.second+1));
        }
        
        Constraint<> incomp_pairs("incomp_pairs");
        incomp_pairs = bin.in(pairs1) + bin.in(pairs2);
        // Reg->add(incomp_pairs.in(range(1,pairs1.size()))<=1);
        //        incomp_pairs.print();
    }
    
    


    Constraint<> txsq("txsq");
    txsq = pow(x_shift,2) - tx;
    //txsq.add_to_callback();
    Reg->add(txsq.in(range(0,0))<=0);

    Constraint<> tysq("tysq");
    tysq = pow(y_shift,2) - ty;
   // tysq.add_to_callback();
    Reg->add(tysq.in(range(0,0))<=0);

    Constraint<> tzsq("tzsq");
    tzsq = pow(z_shift,2) - tz;
    //tzsq.add_to_callback();
    Reg->add(tzsq.in(range(0,0))<=0);
    
    Constraint<> T11sq("T11sq");
    T11sq = pow(theta11,2) - T11;
    //txsq.add_to_callback();
    Reg->add(T11sq.in(range(0,0))<=0);
        
    Constraint<> T12sq("T12sq");
    T12sq = pow(theta12,2) - T12;
    //txsq.add_to_callback();
    Reg->add(T12sq.in(range(0,0))<=0);
    
    Constraint<> T13sq("T13sq");
    T13sq = pow(theta13,2) - T13;
    //txsq.add_to_callback();
    Reg->add(T13sq.in(range(0,0))<=0);
    
    Constraint<> T21sq("T21sq");
    T21sq = pow(theta21,2) - T21;
    //txsq.add_to_callback();
    Reg->add(T21sq.in(range(0,0))<=0);
        
    Constraint<> T22sq("T22sq");
    T22sq = pow(theta22,2) - T22;
    //txsq.add_to_callback();
    Reg->add(T22sq.in(range(0,0))<=0);
    
    Constraint<> T23sq("T23sq");
    T23sq = pow(theta23,2) - T23;
    //txsq.add_to_callback();
    Reg->add(T23sq.in(range(0,0))<=0);
    
    Constraint<> T31sq("T31sq");
    T31sq = pow(theta31,2) - T31;
    //txsq.add_to_callback();
    Reg->add(T31sq.in(range(0,0))<=0);
        
    Constraint<> T32sq("T32sq");
    T32sq = pow(theta32,2) - T32;
    //txsq.add_to_callback();
    Reg->add(T32sq.in(range(0,0))<=0);
    
    Constraint<> T33sq("T33sq");
    T33sq = pow(theta33,2) - T33;
    //txsq.add_to_callback();
    Reg->add(T33sq.in(range(0,0))<=0);
        
    
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
    
    
    param<> deltax_lb("deltax_lb");
    // x_new_lb.in(N1);
    param<> deltax_ub("deltax_ub");
    //x_new_ub.in(N1);
    param<> deltay_lb("deltay_lb");
    //y_new_lb.in(N1);
    param<> deltay_ub("deltay_ub");
    //y_new_ub.in(N1);
    param<> deltaz_lb("deltaz_lb");
    //z_new_lb.in(N1);
    param<> deltaz_ub("deltaz_ub");
    //z_new_ub.in(N1);
    param<> x_new_lb("x_new_lb");
    // x_new_lb.in(N1);
    param<> x_new_ub("x_new_ub");
    //x_new_ub.in(N1);
    param<> y_new_lb("y_new_lb");
    //y_new_lb.in(N1);
    param<> y_new_ub("y_new_ub");
    //y_new_ub.in(N1);
    param<> z_new_lb("z_new_lb");
    //z_new_lb.in(N1);
    param<> z_new_ub("z_new_ub");
    //z_new_ub.in(N1);
    double deltax_min, deltay_min, deltaz_min;
    shared_ptr<pair<double,double>> new_x1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> new_y1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> new_z1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> x1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> y1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> z1_bounds = make_shared<pair<double,double>>();
    for(auto i=0;i<nd;i++){
        x1_bounds->first = point_cloud_data.at(i)[0];
        x1_bounds->second = point_cloud_data.at(i)[0];
        y1_bounds->first = point_cloud_data.at(i)[1];
        y1_bounds->second = point_cloud_data.at(i)[1];
        z1_bounds->first = point_cloud_data.at(i)[2];
        z1_bounds->second = point_cloud_data.at(i)[2];
        auto x_range  = get_product_range(x1_bounds, theta11._range);
        auto y_range  = get_product_range(y1_bounds, theta12._range);
        auto z_range  = get_product_range(z1_bounds, theta13._range);
        *new_x1_bounds = {x_range->first + y_range->first + z_range->first + new_shift_min_x,
            x_range->second + y_range->second + z_range->second+ new_shift_max_x};
        x_new_lb.add_val(to_string(i+1),new_x1_bounds->first);
        x_new_ub.add_val(to_string(i+1),new_x1_bounds->second);
        auto xm_lb = xm_min-new_x1_bounds->second;
        auto xm_ub = xm_max-new_x1_bounds->first;
        if(xm_lb<=0 && xm_ub>=0)
            deltax_min=0;
        else
            deltax_min=std::min(pow(xm_lb,2), pow(xm_ub,2));
        deltax_lb.add_val(to_string(i+1), deltax_min);
        deltax_ub.add_val(to_string(i+1), std::max(pow(xm_lb,2), pow(xm_ub,2)));
        x_range  = get_product_range(x1_bounds, theta21._range);
        y_range  = get_product_range(y1_bounds, theta22._range);
        z_range  = get_product_range(z1_bounds, theta23._range);
        *new_y1_bounds = {x_range->first + y_range->first + z_range->first + new_shift_min_y,
            x_range->second + y_range->second + z_range->second+ new_shift_max_y};
        y_new_lb.add_val(to_string(i+1),new_y1_bounds->first);
        y_new_ub.add_val(to_string(i+1),new_y1_bounds->second);
        auto ym_lb = ym_min-new_y1_bounds->second;
        auto ym_ub = ym_max-new_y1_bounds->first;
        if(ym_lb<=0 && ym_ub>=0)
            deltay_min=0;
        else
            deltay_min=std::min(pow(ym_lb,2), pow(ym_ub,2));
        deltay_lb.add_val(to_string(i+1), deltay_min);
        deltay_ub.add_val(to_string(i+1), std::max(pow(ym_lb,2), pow(ym_ub,2)));
        x_range  = get_product_range(x1_bounds, theta31._range);
        y_range  = get_product_range(y1_bounds, theta32._range);
        z_range  = get_product_range(z1_bounds, theta33._range);
        *new_z1_bounds = {x_range->first + y_range->first + z_range->first + new_shift_min_z,
            x_range->second + y_range->second + z_range->second+ new_shift_max_z};
        z_new_lb.add_val(to_string(i+1),new_z1_bounds->first);
        z_new_ub.add_val(to_string(i+1),new_z1_bounds->second);
        auto zm_lb = zm_min-new_z1_bounds->second;
        auto zm_ub = zm_max-new_z1_bounds->first;
        if(zm_lb<=0 && zm_ub>=0)
            deltaz_min=0;
        else
            deltaz_min=std::min(pow(zm_lb,2), pow(zm_ub,2));
        deltaz_lb.add_val(to_string(i+1), deltaz_min);
        deltaz_ub.add_val(to_string(i+1), std::max(pow(zm_lb,2), pow(zm_ub,2)));
    }
    // var<> new_x1("new_x1", x_new_lb, x_new_ub),new_y1("new_y1", y_new_lb, y_new_ub),new_z1("new_z1", z_new_lb, z_new_ub);
    
    var<> new_x1("new_x1"),new_y1("new_y1"),new_z1("new_z1");
    var<> deltax("deltax"), deltay("deltay"), deltaz("deltaz");
    if(hybrid){
        //var<> deltax("deltax", deltax_lb, deltax_ub), deltay("deltay", deltay_lb, deltay_ub), deltaz("deltaz", deltaz_lb, deltaz_ub);
        
        // double dmax=min_max_model[3].second;
//        Reg->add(new_x1.in(N1));
//        Reg->add(new_y1.in(N1));
//        Reg->add(new_z1.in(N1));
        Reg->add(new_xm.in(N1));
        Reg->add(new_ym.in(N1));
        Reg->add(new_zm.in(N1));
        Reg->add(delta.in(R(1)));
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
        
        auto ids1 = theta11.repeat_id(N1.size());
        Constraint<> Def_newxm("Def_newxm");
        Def_newxm = new_xm-product(x2.in(ids),bin.in_matrix(1, 1));
        Def_newxm += x1.in(N1)*theta11.in(ids1) + y1.in(N1)*theta12.in(ids1) + z1.in(N1)*theta13.in(ids1)+x_shift;
        Reg->add(Def_newxm.in(N1)==0);
        
        Constraint<> Def_newym("Def_newym");
        Def_newym = new_ym-product(y2.in(ids),bin.in_matrix(1, 1));
        Def_newym += x1.in(N1)*theta21.in(ids1) + y1.in(N1)*theta22.in(ids1) + z1.in(N1)*theta23.in(ids1)+y_shift;
        Reg->add(Def_newym.in(N1)==0);
        
        Constraint<> Def_newzm("Def_newzm");
        Def_newzm = new_zm-product(z2.in(ids),bin.in_matrix(1, 1));
        Def_newzm += x1.in(N1)*theta31.in(ids1) + y1.in(N1)*theta32.in(ids1) + z1.in(N1)*theta33.in(ids1)+z_shift;
        Reg->add(Def_newzm.in(N1)==0);
        
        Constraint<> sum_newxm("sum_newxm");
        sum_newxm = sum(new_xm.in(N1))-nd*x_shift;
        //Reg->add(sum_newxm==0);
        
        Constraint<> sum_newym("sum_newym");
        sum_newym = sum(new_ym.in(N1))-nd*y_shift;
        //Reg->add(sum_newym==0);
        
        Constraint<> sum_newzm("sum_newzm");
        sum_newzm = sum(new_zm.in(N1))-nd*z_shift;
        //Reg->add(sum_newzm==0);
        
        //auto ids1 = theta11.repeat_id(N1.size());
        Constraint<> x_rot1("x_rot1");
        x_rot1 += new_x1 -x_shift;
        x_rot1 -= x1.in(N1)*theta11.in(ids1) + y1.in(N1)*theta12.in(ids1) + z1.in(N1)*theta13.in(ids1);
        //Reg->add(x_rot1.in(N1)==0);
        
        Constraint<> y_rot1("y_rot1");
        y_rot1 += new_y1 - y_shift;
        y_rot1 -= x1.in(N1)*theta21.in(ids1) + y1.in(N1)*theta22.in(ids1) + z1.in(N1)*theta23.in(ids1);
        //Reg->add(y_rot1.in(N1)==0);
        
        Constraint<> z_rot1("z_rot1");
        z_rot1 += new_z1 -z_shift;
        z_rot1 -= x1.in(N1)*theta31.in(ids1) + y1.in(N1)*theta32.in(ids1) + z1.in(N1)*theta33.in(ids1);
        //Reg->add(z_rot1.in(N1)==0);
        
        param<> ub_par("ub_par");
        ub_par=ub;
        Constraint<> cut_off("cut_off");
        cut_off=delta-ub_par;
        Reg->add(cut_off.in(range(0,0))<=0);
        
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
        Def_delta=sum(deltax)+sum(deltay)+sum(deltaz)-delta;
        Reg->add(Def_delta<=0);
        
        if(hybrid){
            auto idstheta = theta11.repeat_id(N1.size());
            Constraint<> limit_neg("limit_neg");
            limit_neg=2*(new_xm*x1*theta11.in(idstheta));
            limit_neg+= 2*(new_xm*y1*theta12.in(idstheta));
            limit_neg+= 2*(new_xm*z1*theta13.in(idstheta));
            limit_neg+= 2*(new_ym*x1*theta21.in(idstheta));
            limit_neg+= 2*(new_ym*y1*theta22.in(idstheta));
            limit_neg+= 2*(new_ym*z1*theta23.in(idstheta));
            limit_neg+= 2*(new_zm*x1*theta31.in(idstheta));
            limit_neg+= 2*(new_zm*y1*theta32.in(idstheta));
            limit_neg+= 2*(new_zm*z1*theta33.in(idstheta));
            limit_neg-=pow(x1,2)+pow(y1,2)+pow(z1,2);
            limit_neg-=pow(new_xm,2)+pow(new_ym,2)+pow(new_zm,2);
            // Reg->add(limit_neg.in(N1)<=0);
            
            Constraint<> limit_pos("limit_pos");
            limit_pos-=2*(new_xm*x1*theta11.in(idstheta));
            limit_pos-= 2*(new_xm*y1*theta12.in(idstheta));
            limit_pos-= 2*(new_xm*z1*theta13.in(idstheta));
            limit_pos-= 2*(new_ym*x1*theta21.in(idstheta));
            limit_pos-= 2*(new_ym*y1*theta22.in(idstheta));
            limit_pos-= 2*(new_ym*z1*theta23.in(idstheta));
            limit_pos-= 2*(new_zm*x1*theta31.in(idstheta));
            limit_pos-= 2*(new_zm*y1*theta32.in(idstheta));
            limit_pos-= 2*(new_zm*z1*theta33.in(idstheta));
            limit_pos-=pow(x1,2)+pow(y1,2)+pow(z1,2);
            limit_pos-=pow(new_xm,2)+pow(new_ym,2)+pow(new_zm,2);
            // Reg->add(limit_pos.in(N1)<=0);
            
            var<> ex("ex");
            // Reg->add(ex.in(N1));
            param<> dm("dm");
            param<> d_sq_root("d_sq_root");
            for(auto i=0;i<nm;i++){
                auto dmd=pow(point_cloud_model.at(i)[0],2)+pow(point_cloud_model.at(i)[1],2)+pow(point_cloud_model.at(i)[2],2);
                auto dmd_r=sqrt(dmd);
                dm.add_val(to_string(i+1), dmd);
                d_sq_root.add_val(to_string(i+1), dmd_r);
            }
            
            param<> cd("cd");
            for(auto i=0;i<nd;i++){
                auto dmd=pow(point_cloud_data.at(i)[0],2)+pow(point_cloud_data.at(i)[1],2)+pow(point_cloud_data.at(i)[2],2);
                auto dmd_r=2*(sqrt(dmd)+sqrt(shift_mag_max));
                cd.add_val(to_string(i+1), dmd_r);
            }
            
            indices idsij = indices("idsij");
            idsij.add_empty_row();
            
            param<> dm_root("dm_root");
            dm_root.in(cells);
            for(auto i=0;i<nd;i++){
                auto di=pow(point_cloud_data.at(i)[0],2)+pow(point_cloud_data.at(i)[1],2)+pow(point_cloud_data.at(i)[2],2);
                auto di_r=2*(sqrt(di)+sqrt(shift_mag_max));
                for(auto j=1;j<=nm;j++){
                    if(cells.has_key(to_string(i+1)+","+to_string(j))){
                        auto dmd=pow(point_cloud_model.at(j-1)[0],2)+pow(point_cloud_model.at(j-1)[1],2)+pow(point_cloud_model.at(j-1)[2],2);
                        auto dmd_r=dmd-sqrt(dmd)*di_r;
                        dm_root.set_val(to_string(i+1)+","+to_string(j), dmd_r);
                        idsij.add_in_row(i, to_string(i+1)+","+to_string(j));
                    }
                }
            }
            
            param<> cL("cL");
            cL.in(cells);
            param<> cU("cU");
            cU.in(cells);
            param<> dm_cells("dm_cells");
            dm_cells.in(cells);
            
            vector<vector<vector<double>>> box;
            vector<vector<double>> box_i;
            vector<double> coord_i;
            coord_i.resize(3);
            
            //            shared_ptr<pair<double,double>> x1_bounds = make_shared<pair<double,double>>();
            //            shared_ptr<pair<double,double>> y1_bounds = make_shared<pair<double,double>>();
            //            shared_ptr<pair<double,double>> z1_bounds = make_shared<pair<double,double>>();
            vector<double> x_lb, x_ub, y_lb, y_ub, z_lb, z_ub;
            
            for(auto i=0;i<nd;i++){
                //  auto max_dist_i = get_max_dist(roll_min, roll_max, pitch_min, pitch_max, yaw_min, yaw_max, shift_min_x, shift_max_x, shift_min_y, shift_max_y, shift_min_z, shift_max_z, point_cloud_data[i], zeros, false);
                //dist_i.push_back(max_dist_i);
                // sphere_outer_sq.push_back(pow((d_root+shift_mag_max_root),2));
                x1_bounds->first = point_cloud_data.at(i)[0];
                x1_bounds->second = point_cloud_data.at(i)[0];
                y1_bounds->first = point_cloud_data.at(i)[1];
                y1_bounds->second = point_cloud_data.at(i)[1];
                z1_bounds->first = point_cloud_data.at(i)[2];
                z1_bounds->second = point_cloud_data.at(i)[2];
                auto x_range  = get_product_range(x1_bounds, theta11._range);
                auto y_range  = get_product_range(y1_bounds, theta12._range);
                auto z_range  = get_product_range(z1_bounds, theta13._range);
                x_lb.push_back(x_range->first + y_range->first + z_range->first + new_shift_min_x);
                x_ub.push_back(x_range->second + y_range->second + z_range->second+ new_shift_max_x);
                x_range  = get_product_range(x1_bounds, theta21._range);
                y_range  = get_product_range(y1_bounds, theta22._range);
                z_range  = get_product_range(z1_bounds, theta23._range);
                y_lb.push_back(x_range->first + y_range->first + z_range->first + new_shift_min_y);
                y_ub.push_back(x_range->second + y_range->second + z_range->second+ new_shift_max_y);
                x_range  = get_product_range(x1_bounds, theta31._range);
                y_range  = get_product_range(y1_bounds, theta32._range);
                z_range  = get_product_range(z1_bounds, theta33._range);
                z_lb.push_back(x_range->first + y_range->first + z_range->first + new_shift_min_z);
                z_ub.push_back(x_range->second + y_range->second + z_range->second+ new_shift_max_z);
                coord_i[0]=x_lb[i];
                coord_i[1]=y_lb[i];
                coord_i[2]=z_lb[i];
                box_i.push_back(coord_i);
                coord_i[0]=x_lb[i];
                coord_i[1]=y_ub[i];
                coord_i[2]=z_lb[i];
                box_i.push_back(coord_i);
                coord_i[0]=x_ub[i];
                coord_i[1]=y_ub[i];
                coord_i[2]=z_lb[i];
                box_i.push_back(coord_i);
                coord_i[0]=x_ub[i];
                coord_i[1]=y_lb[i];
                coord_i[2]=z_lb[i];
                box_i.push_back(coord_i);
                coord_i[0]=x_lb[i];
                coord_i[1]=y_lb[i];
                coord_i[2]=z_ub[i];
                box_i.push_back(coord_i);
                coord_i[0]=x_lb[i];
                coord_i[1]=y_ub[i];
                coord_i[2]=z_ub[i];
                box_i.push_back(coord_i);
                coord_i[0]=x_ub[i];
                coord_i[1]=y_ub[i];
                coord_i[2]=z_ub[i];
                box_i.push_back(coord_i);
                coord_i[0]=x_ub[i];
                coord_i[1]=y_lb[i];
                coord_i[2]=z_ub[i];
                box_i.push_back(coord_i);
                auto di=pow(point_cloud_data.at(i)[0],2)+pow(point_cloud_data.at(i)[1],2)+pow(point_cloud_data.at(i)[2],2);
                auto di_r=2*(sqrt(di)+sqrt(shift_mag_max));
                for(auto j=1;j<=nm;j++){
                    if(cells.has_key(to_string(i+1)+","+to_string(j))){
                        auto x=point_cloud_model.at(j-1)[0];
                        auto y=point_cloud_model.at(j-1)[1];
                        auto z=point_cloud_model.at(j-1)[2];
                        double c_min=100, c_max=-999;
                        for(auto k=0;k<8;k++){
                            auto inner_prod=2*(box_i[k][0]*x+box_i[k][1]*y+box_i[k][2]*z);
                            if(inner_prod>=c_max){
                                c_max=inner_prod;
                            }
                            if(inner_prod<=c_min){
                                c_min=inner_prod;
                            }
                        }
                        c_max=std::min(c_max,sqrt(pow(x,2)+pow(y,2)+pow(z,2))*di_r);
                        cL.set_val(to_string(i+1)+","+to_string(j), c_min);
                        cU.set_val(to_string(i+1)+","+to_string(j), c_max);
                        dm_cells.set_val(to_string(i+1)+","+to_string(j), pow(x,2)+pow(y,2)+pow(z,2)-c_max);
                        auto dmd_r=sqrt(pow(x,2)+pow(y,2)+pow(z,2))*di_r-c_max;
                        if(dmd_r<-0.0001){
                            DebugOn("dmdr "<<dmd_r<<endl);
                        }
                    }
                }
                box_i.clear();
            }
            
            
            Constraint<> delta_lower("delta_lower");
            //delta_lower=(x1*x1)+(y1*y1)+(z1*z1)+(product(dm_root.in(idsij), bin.in_matrix(1, 1)))-delta;
            //delta_lower=sum(x1*x1)+sum(y1*y1)+sum(z1*z1)+nd*(x_shift*x_shift+y_shift*y_shift+z_shift*z_shift)+sum(dm_root*bin)-delta;
            delta_lower=sum(x1*x1)+sum(y1*y1)+sum(z1*z1)+nd*(tx+ty+tz)+sum(dm_cells*bin)-delta;
            Reg->add(delta_lower<=0);
            
            
            
            Constraint<> dist_model("dist_model");
            dist_model=pow(new_xm,2)+pow(new_ym,2)+pow(new_zm,2)-product(dm.in(ids),bin.in_matrix(1, 1));
            // Reg->add(dist_model.in(N1)==0);
            
            
            Constraint<> dist_data1("dist_data1");
            dist_data1=new_x1*new_x1+new_y1*new_y1+new_z1*new_z1-pow((sqrt(pow(x1,2)+pow(y1,2)+pow(z1,2))+sqrt(shift_mag_max)),2);
            // Reg->add(dist_data1.in(N1)<=0);
            
            Constraint<> dist_data2("dist_data2");
            dist_data2=x_new_lb*x_new_ub+y_new_lb*y_new_ub+z_new_lb*z_new_ub+pow((sqrt(pow(x1,2)+pow(y1,2)+pow(z1,2))-sqrt(shift_mag_max)),2)-new_x1*(x_new_ub+x_new_lb)-new_y1*(y_new_ub+y_new_lb)-new_z1*(z_new_ub+z_new_lb);
           // Reg->add(dist_data2.in(N1)<=0);
            
            
            
            
            
            Constraint<> def_ex("def_ex");
            def_ex=ex-product(d_sq_root.in(ids),bin.in_matrix(1, 1));
            //Reg->add(def_ex.in(N1)==0);
            
            
            
            Constraint<> limit_neg_bin("limit_neg_bin");
            limit_neg_bin=2*(new_xm*x1*theta11.in(idstheta));
            limit_neg_bin+= 2*(new_xm*y1*theta12.in(idstheta));
            limit_neg_bin+= 2*(new_xm*z1*theta13.in(idstheta));
            limit_neg_bin+= 2*(new_ym*x1*theta21.in(idstheta));
            limit_neg_bin+= 2*(new_ym*y1*theta22.in(idstheta));
            limit_neg_bin+= 2*(new_ym*z1*theta23.in(idstheta));
            limit_neg_bin+= 2*(new_zm*x1*theta31.in(idstheta));
            limit_neg_bin+= 2*(new_zm*y1*theta32.in(idstheta));
            limit_neg_bin+= 2*(new_zm*z1*theta33.in(idstheta));
            limit_neg_bin-=pow(x1,2)+pow(y1,2)+pow(z1,2);
            limit_neg_bin-=product(dm.in(ids),bin.in_matrix(1, 1));
            // Reg->add(limit_neg_bin.in(N1)<=0);
            
        }
    }
    bool spatial_branch=false;
    if(spatial_branch)
    {
        var<> yaw("yaw", new_yaw_min, new_yaw_max), pitch("pitch", new_pitch_min, new_pitch_max), roll("roll", new_roll_min, new_roll_max);
        yaw.in(R(1)); pitch.in(R(1));roll.in(R(1));
        func<> cosr_f = cos(roll);
        func<> sinr_f = sin(roll);
        func<> cosp_f = cos(pitch);
        func<> sinp_f = sin(pitch);
        func<> cosy_f = cos(yaw);
        func<> siny_f = sin(yaw);
        var<> cosr("cosr", cosr_f._range->first, cosr_f._range->second), sinr("sinr", sinr_f._range->first, sinr_f._range->second);
        var<> cosp("cosp",  cosp_f._range->first, cosp_f._range->second), sinp("sinp", sinp_f._range->first, sinp_f._range->second);
        var<> cosy("cosy",  cosy_f._range->first, cosy_f._range->second), siny("siny", siny_f._range->first, siny_f._range->second);
        auto cosy_sinr_range = get_product_range(cosy._range, sinr._range);
        auto siny_sinr_range = get_product_range(siny._range, sinr._range);
        var<> cosy_sinr("cosy_sinr", cosy_sinr_range->first, cosy_sinr_range->second), siny_sinr("siny_sinr", siny_sinr_range->first, siny_sinr_range->second);
        
        Reg->add(cosr.in(R(1)),cosp.in(R(1)),cosy.in(R(1)));
        Reg->add(sinr.in(R(1)),sinp.in(R(1)),siny.in(R(1)));
        Reg->add(cosy_sinr.in(R(1)),siny_sinr.in(R(1)));
        
        Constraint<> cosy_sinr_prod("cosy_sinr");
        cosy_sinr_prod = cosy_sinr - cosy*sinr;
        Reg->add(cosy_sinr_prod==0);
        
        Constraint<> siny_sinr_prod("siny_sinr");
        siny_sinr_prod = siny_sinr - siny*sinr;
        Reg->add(siny_sinr_prod==0);
        
        
        Constraint<> R11("R11");
        R11 += theta11 - cosy*cosr;
        Reg->add(R11==0);
        
        Constraint<> R12("R12");
        R12 += theta12 - (cosy_sinr*sinp - siny*cosp);
        Reg->add(R12==0);
        
        Constraint<> R13("R13");
        R13 += theta13 - (cosy_sinr*cosp + siny*sinp);
        Reg->add(R13==0);
        
        Constraint<> R21("R21");
        R21 += theta21 - siny*cosr;
        Reg->add(R21==0);
        
        Constraint<> R22("R22");
        R22 += theta22 - (siny_sinr*sinp + cosy*cosp);
        Reg->add(R22==0);
        
        Constraint<> R23("R23");
        R23 += theta23 - (siny_sinr*cosp - cosy*sinp);
        Reg->add(R23==0);
        
        Constraint<> R31("R31");
        R31 += theta31 + sinr;
        Reg->add(R31==0);
        
        Constraint<> R32("R32");
        R32 += theta32 - cosr*sinp;
        Reg->add(R32==0);
        
        Constraint<> R33("R33");
        R33 += theta33 - cosr*cosp;
        Reg->add(R33==0);
        
        /*
         func<> r11 = cos(yaw)*cos(roll);r11.eval_all();
         func<> r12 = cos(yaw)*sin(roll)*sin(pitch) - sin(yaw)*cos(pitch);r12.eval_all();
         func<> r13 = cos(yaw)*sin(roll)*cos(pitch) + sin(yaw)*sin(pitch);r13.eval_all();
         func<> r21 = sin(yaw)*cos(roll);r21.eval_all();
         func<> r22 = sin(yaw)*sin(roll)*sin(pitch) + cos(yaw)*cos(pitch);r22.eval_all();
         func<> r23 = sin(yaw)*sin(roll)*cos(pitch) - cos(yaw)*sin(pitch);r23.eval_all();
         func<> r31 = sin(-1*roll);r31.eval_all();
         func<> r32 = cos(roll)*sin(pitch);r32.eval_all();
         func<> r33 = cos(roll)*cos(pitch);r33.eval_all();
         
         */
        
        bool add_spatial_bins = true;
        
        if(add_spatial_bins){
            /* Spatial branching vars */
            int nb_pieces = 10; // Divide each axis into nb_pieces
            indices spatial_ids("spatial_ids");
            spatial_ids = range(1,nb_pieces);
            indices theta_ids("theta_ids");
            theta_ids = indices(range(1,3),range(1,3));
            indices shift_ids("shift_ids");
            shift_ids.insert({"x", "y", "z"});
            indices theta_spatial("theta_spatial");
            theta_spatial = indices(theta_ids, spatial_ids);
            
            indices shift_spatial("shift_spatial");
            shift_spatial = indices(shift_ids, spatial_ids);
            
            //    var<int> sbin_shift("sbin_shift", 0, 1);
            //    var<int> sbin_theta("sbin_theta", 0, 1);
            
            //    Reg->add(sbin_theta.in(theta_spatial),sbin_shift.in(shift_spatial));
            
            var<int> sbin_tx("sbin_tx", 0, 1), sbin_ty("sbin_ty", 0, 1), sbin_tz("sbin_tz", 0, 1);
            var<int> sbin_roll("sbin_roll", 0, 1), sbin_pitch("sbin_pitch", 0, 1), sbin_yaw("sbin_yaw", 0, 1);
            sbin_roll._priority = 1e6;sbin_pitch._priority = 1e6;sbin_yaw._priority = 1e6;
            sbin_tx._priority = 2e6;sbin_ty._priority = 2e6;sbin_tz._priority = 2e6;
            Reg->add(sbin_roll.in(spatial_ids),sbin_pitch.in(spatial_ids),sbin_yaw.in(spatial_ids));
            Reg->add(sbin_tz.in(spatial_ids), sbin_tx.in(spatial_ids),sbin_ty.in(spatial_ids));
            /* Spatial branching constraints */
            Constraint<> OneBinRoll("OneBinRoll");
            OneBinRoll = sum(sbin_roll);
            Reg->add(OneBinRoll==1);
            
            Constraint<> OneBinPitch("OneBinPitch");
            OneBinPitch = sum(sbin_pitch);
            Reg->add(OneBinPitch==1);
            
            Constraint<> OneBinYaw("OneBinYaw");
            OneBinYaw = sum(sbin_yaw);
            Reg->add(OneBinYaw==1);
            
            
            Constraint<> OneBinShiftSpatialx("OneBinShiftSpatialx");
            OneBinShiftSpatialx = sum(sbin_tx);
            Reg->add(OneBinShiftSpatialx==1);
            
            Constraint<> OneBinShiftSpatialy("OneBinShiftSpatialy");
            OneBinShiftSpatialy = sum(sbin_ty);
            Reg->add(OneBinShiftSpatialy==1);
            
            Constraint<> OneBinShiftSpatialz("OneBinShiftSpatialz");
            OneBinShiftSpatialz = sum(sbin_tz);
            Reg->add(OneBinShiftSpatialz==1);
            
            //    Reg->print();
            double angle_min = std::min(std::min(new_roll_min,new_pitch_min),new_yaw_min), angle_max = std::max(std::max(new_roll_max,new_pitch_max),new_yaw_max);
            double t_min = std::min(std::min(new_shift_min_x,new_shift_min_y),new_shift_min_z), t_max = std::max(std::max(new_shift_max_x,new_shift_max_y),new_shift_max_z);
            double angle_increment = (angle_max - angle_min)/nb_pieces;/* Angles are defined in [-2,2] in radians (-120,120) in degrees */
            double shift_increment = (t_max - t_min)/nb_pieces;/* Shifts are defined in [-0.5,0.5] */
            
            auto spatial_ids_n = range(1,nb_pieces-1);
            auto spatial_ids_1 = range(2,nb_pieces);
            param<> angle_lb("angle_lb"), angle_ub("angle_ub");
            param<> cos_lb("cos_lb"), cos_ub("cos_ub");
            param<> sin_lb("sin_lb"), sin_ub("sin_ub");
            param<> t_lb("t_lb"), t_ub("t_ub");
            angle_ub.in(spatial_ids);
            angle_lb.in(spatial_ids);
            cos_ub.in(spatial_ids);
            cos_lb.in(spatial_ids);
            sin_ub.in(spatial_ids);
            sin_lb.in(spatial_ids);
            t_ub.in(spatial_ids);
            t_lb.in(spatial_ids);
            double lb = 0, ub = 0;
            for (int i = 0; i<nb_pieces; i++) {
                lb = angle_min + i*angle_increment;
                ub = angle_min+(i+1)*angle_increment;
                angle_lb.set_val(i,lb);
                angle_ub.set_val(i,ub);
                cos_lb.set_val(i,std::min(cos(lb), cos(ub)));
                cos_ub.set_val(i,std::max(cos(lb), cos(ub)));
                if(lb < 0 && ub > 0){/* zero is in the domain, i.e., cos max is 1 */
                    cos_ub.set_val(i,1);
                }
                if((lb < -pi && ub > -pi) || (lb < pi && ub > pi)){/* -pi or pi is in the domain, i.e., cos min is -1 */
                    cos_lb.set_val(i,-1);
                }
                sin_lb.set_val(i,std::min(sin(lb), sin(ub)));
                sin_ub.set_val(i,std::max(sin(lb), sin(ub)));
                if((lb < -3*pi/2 && ub > -3*pi/2) || (lb < pi/2 && ub > pi/2)){/* -3pi/2 or pi/2 is in the domain, i.e., sin max is 1 */
                    sin_ub.set_val(i,1);
                }
                if((lb < 3*pi/2 && ub > 3*pi/2) || (lb < -pi/2 && ub > -pi/2)){/* 3pi/2 or -pi/2 is in the domain, i.e., sin min is -1 */
                    sin_lb.set_val(i,-1);
                }
                t_lb.set_val(i,t_min + i*shift_increment);
                t_ub.set_val(i,t_min + (i+1)*shift_increment);
            }
            Reg->add_param(angle_lb);
            Reg->add_param(angle_ub);
            Reg->add_param(cos_lb);
            Reg->add_param(cos_ub);
            Reg->add_param(sin_lb);
            Reg->add_param(sin_ub);
            Reg->add_param(t_lb);
            Reg->add_param(t_ub);
            
            auto ids_repeat = theta11.repeat_id(nb_pieces-1);
            
            
            Constraint<> CosRoll_UB("CosRoll_UB");
            CosRoll_UB = cosr.in(ids_repeat) - cos_ub.in(spatial_ids_n);
            Reg->add_on_off_multivariate_refined(CosRoll_UB.in(spatial_ids_n)<=0, sbin_roll.in(spatial_ids_n), true);
            
            Constraint<> CosRoll_LB("CosRoll_LB");
            CosRoll_LB = cos_lb.in(spatial_ids_1) - cosr.in(ids_repeat);
            Reg->add_on_off_multivariate_refined(CosRoll_LB.in(spatial_ids_1) <= 0, sbin_roll.in(spatial_ids_1), true);
            
            Constraint<> SinRoll_UB("SinRoll_UB");
            SinRoll_UB = sinr.in(ids_repeat) - sin_ub.in(spatial_ids_n);
            Reg->add_on_off_multivariate_refined(SinRoll_UB.in(spatial_ids_n)<=0, sbin_roll.in(spatial_ids_n), true);
            
            Constraint<> SinRoll_LB("SinRoll_LB");
            SinRoll_LB = sin_lb.in(spatial_ids_1) - sinr.in(ids_repeat);
            Reg->add_on_off_multivariate_refined(SinRoll_LB.in(spatial_ids_1) <= 0, sbin_roll.in(spatial_ids_1), true);
            
            Constraint<> CosPitch_UB("CosPitch_UB");
            CosPitch_UB = cosp.in(ids_repeat) - cos_ub.in(spatial_ids_n);
            Reg->add_on_off_multivariate_refined(CosPitch_UB.in(spatial_ids_n)<=0, sbin_pitch.in(spatial_ids_n), true);
            
            Constraint<> CosPitch_LB("CosPitch_LB");
            CosPitch_LB = cos_lb.in(spatial_ids_1) - cosp.in(ids_repeat);
            Reg->add_on_off_multivariate_refined(CosPitch_LB.in(spatial_ids_1) <= 0, sbin_pitch.in(spatial_ids_1), true);
            
            Constraint<> SinPitch_UB("SinPitch_UB");
            SinPitch_UB = sinp.in(ids_repeat) - sin_ub.in(spatial_ids_n);
            Reg->add_on_off_multivariate_refined(SinPitch_UB.in(spatial_ids_n)<=0, sbin_pitch.in(spatial_ids_n), true);
            
            Constraint<> SinPitch_LB("SinPitch_LB");
            SinPitch_LB = sin_lb.in(spatial_ids_1) - sinp.in(ids_repeat);
            Reg->add_on_off_multivariate_refined(SinPitch_LB.in(spatial_ids_1) <= 0, sbin_pitch.in(spatial_ids_1), true);
            
            
            Constraint<> CosYaw_UB("CosYaw_UB");
            CosYaw_UB = cosy.in(ids_repeat) - cos_ub.in(spatial_ids_n);
            Reg->add_on_off_multivariate_refined(CosYaw_UB.in(spatial_ids_n)<=0, sbin_yaw.in(spatial_ids_n), true);
            
            Constraint<> CosYaw_LB("CosYaw_LB");
            CosYaw_LB = cos_lb.in(spatial_ids_1) - cosy.in(ids_repeat);
            Reg->add_on_off_multivariate_refined(CosYaw_LB.in(spatial_ids_1) <= 0, sbin_yaw.in(spatial_ids_1), true);
            
            Constraint<> SinYaw_UB("SinYaw_UB");
            SinYaw_UB = siny.in(ids_repeat) - sin_ub.in(spatial_ids_n);
            Reg->add_on_off_multivariate_refined(SinYaw_UB.in(spatial_ids_n)<=0, sbin_yaw.in(spatial_ids_n), true);
            
            Constraint<> SinYaw_LB("SinYaw_LB");
            SinYaw_LB = sin_lb.in(spatial_ids_1) - siny.in(ids_repeat);
            Reg->add_on_off_multivariate_refined(SinYaw_LB.in(spatial_ids_1) <= 0, sbin_yaw.in(spatial_ids_1), true);
            
            Constraint<> xshift_Spatial_UB("xshift_Spatial_UB");
            xshift_Spatial_UB = x_shift.in(ids_repeat) - t_ub.in(spatial_ids_n);
            Reg->add_on_off_multivariate_refined(xshift_Spatial_UB.in(spatial_ids_n) <= 0, sbin_tx.in(spatial_ids_n), true);
            
            Constraint<> xshift_Spatial_LB("xshift_Spatial_LB");
            xshift_Spatial_LB = t_lb.in(spatial_ids_1) - x_shift.in(ids_repeat);
            Reg->add_on_off_multivariate_refined(xshift_Spatial_LB.in(spatial_ids_1) <= 0, sbin_tx.in(spatial_ids_1), true);
            
            Constraint<> yshift_Spatial_UB("yshift_Spatial_UB");
            yshift_Spatial_UB = y_shift.in(ids_repeat) - t_ub.in(spatial_ids_n);
            Reg->add_on_off_multivariate_refined(yshift_Spatial_UB.in(spatial_ids_n) <= 0, sbin_ty.in(spatial_ids_n), true);
            
            Constraint<> yshift_Spatial_LB("yshift_Spatial_LB");
            yshift_Spatial_LB = t_lb.in(spatial_ids_1) - y_shift.in(ids_repeat);
            Reg->add_on_off_multivariate_refined(yshift_Spatial_LB.in(spatial_ids_1) <= 0, sbin_ty.in(spatial_ids_1), true);
            
            Constraint<> zshift_Spatial_UB("zshift_Spatial_UB");
            zshift_Spatial_UB = z_shift.in(ids_repeat) - t_ub.in(spatial_ids_n);
            Reg->add_on_off_multivariate_refined(zshift_Spatial_UB.in(spatial_ids_n) <= 0, sbin_tz.in(spatial_ids_n), true);
            
            Constraint<> zshift_Spatial_LB("zshift_Spatial_LB");
            zshift_Spatial_LB = t_lb.in(spatial_ids_1) - z_shift.in(ids_repeat);
            Reg->add_on_off_multivariate_refined(zshift_Spatial_LB.in(spatial_ids_1) <= 0, sbin_tz.in(spatial_ids_1), true);
        }
        //
        Constraint<> trigR("trigR");
        trigR = pow(cosr,2) + pow(sinr,2);
        trigR.add_to_callback();
        Reg->add(trigR.in(range(1,1))<=1);
        
        Constraint<> trigP("trigP");
        trigP = pow(cosp,2) + pow(sinp,2);
        trigP.add_to_callback();
        Reg->add(trigP.in(range(1,1))<=1);
        
        Constraint<> trigY("trigY");
        trigY = pow(cosy,2) + pow(siny,2);
        trigY.add_to_callback();
        Reg->add(trigY.in(range(1,1))<=1);
        
        
        Constraint<> trigR_NC("trigR_NC");
        trigR_NC = pow(cosr,2) + pow(sinr,2);
        Reg->add(trigR_NC.in(range(1,1))>=1);
        
        Constraint<> trigP_NC("trigP_NC");
        trigP_NC = pow(cosp,2) + pow(sinp,2);
        Reg->add(trigP_NC.in(range(1,1))>=1);
        
        Constraint<> trigY_NC("trigY_NC");
        trigY_NC = pow(cosy,2) + pow(siny,2);
        Reg->add(trigY_NC.in(range(1,1))>=1);
        //        Reg->print();
    }
    /* Objective function */
    
    func<> obj = 0;
    
    //obj+= sum(x1*x1 + y1*y1 + z1*z1);
    //obj += sum(x2.to(cells)*x2.to(cells)*bin) + sum(y2.to(cells)*y2.to(cells)*bin) + sum(z2.to(cells)*z2.to(cells)*bin);
    
    // obj -=nd*(x_shift*x_shift+y_shift*y_shift+z_shift*z_shift);
    
    //  obj -= 2*x_shift*sum_xm+2*y_shift*sum_ym+2*z_shift*sum_zm;
    
    //        auto ids_repeat1 = theta11.repeat_id(N1.size());
    //        obj -= 2*sum(new_xm*x1*theta11.in(ids_repeat1));
    //        obj -= 2*sum(new_xm*y1*theta12.in(ids_repeat1));
    //        obj -= 2*sum(new_xm*z1*theta13.in(ids_repeat1));
    //        obj -= 2*sum(new_ym*x1*theta21.in(ids_repeat1));
    //        obj -= 2*sum(new_ym*y1*theta22.in(ids_repeat1));
    //        obj -= 2*sum(new_ym*z1*theta23.in(ids_repeat1));
    //        obj -= 2*sum(new_zm*x1*theta31.in(ids_repeat1));
    //        obj -= 2*sum(new_zm*y1*theta32.in(ids_repeat1));
    //        obj -= 2*sum(new_zm*z1*theta33.in(ids_repeat1));
    param<> one("one");
    one.in(N1);
    one = 1;
    //obj+=sum(deltax)+sum(deltay)+sum(deltaz);
    obj+=(delta);
    
    
    Reg->min(obj);
    
//    if(model_init){
//        auto binc = model_parent->get_var<int>("bin");
//        auto v_keys_bin=binc.get_keys();
//        for(auto it_key=v_keys_bin->begin(); it_key!=v_keys_bin->end(); it_key++){
//            auto key = *it_key;
//            if(cells.has_key(key)){
//                bin._val->at(bin._indices->_keys_map->at(key)) = binc.eval(key);
//            }
//        }
//        auto xmc = model_parent->get_var<double>("new_xm");
//        auto ymc = model_parent->get_var<double>("new_ym");
//        auto zmc = model_parent->get_var<double>("new_zm");
//        auto v_keys_xm=xmc.get_keys();
//        for(auto it_key=v_keys_xm->begin(); it_key!=v_keys_xm->end(); it_key++){
//            auto key = *it_key;
//            new_xm._val->at(new_xm._indices->_keys_map->at(key)) = xmc.eval(key);
//            new_ym._val->at(new_ym._indices->_keys_map->at(key)) = ymc.eval(key);
//            new_zm._val->at(new_zm._indices->_keys_map->at(key)) = zmc.eval(key);
//        }
//        auto xshiftc = model_parent->get_var<double>("x_shift");
//        auto yshiftc = model_parent->get_var<double>("y_shift");
//        auto zshiftc = model_parent->get_var<double>("z_shift");
//        x_shift.initialize_all(xshiftc.eval());
//        y_shift.initialize_all(xshiftc.eval());
//        z_shift.initialize_all(xshiftc.eval());
//        auto txc = model_parent->get_var<double>("tx");
//        auto tyc = model_parent->get_var<double>("ty");
//        auto tzc = model_parent->get_var<double>("tz");
//        tx.initialize_all(tx.eval());
//        ty.initialize_all(ty.eval());
//        tz.initialize_all(tz.eval());
//        auto dc = model_parent->get_var<double>("delta");
//        tz.initialize_all(tz.eval());
//    }
        
    //    for (int i = 0; i<nd; i++) {
    //        string key = to_string(i+1)+","+to_string(i+1);
    //        bin._val->at(bin._indices->_keys_map->at(key)) = 1;
    //    }
    //}
    
//        Reg->print();
//        solver<> S(Reg,gurobi);
//        S.use_callback();
//        S.run();
//        Reg->print_solution();
//        DebugOn("Theta matrix = " << endl);
//        DebugOn("|" << theta11.eval() << " " << theta12.eval() << " " << theta13.eval() << "|" << endl);
//        DebugOn("|" << theta21.eval() << " " << theta22.eval() << " " << theta23.eval() << "|" << endl);
//        DebugOn("|" << theta31.eval() << " " << theta32.eval() << " " << theta33.eval() << "|" << endl);
//        DebugOn("row 1 " << pow(theta11.eval(),2)+pow(theta12.eval(),2)+pow(theta13.eval(),2)
//                << endl);
//        DebugOn("row 2 " << pow(theta21.eval(),2)+pow(theta22.eval(),2)+pow(theta23.eval(),2)
//                << endl);
//        DebugOn("row 3 " << pow(theta31.eval(),2)+pow(theta32.eval(),2)+pow(theta33.eval(),2)
//                << endl);
//        DebugOn("col 1 " << pow(theta11.eval(),2)+pow(theta21.eval(),2)+pow(theta31.eval(),2)
//                << endl);
//        DebugOn("col 2 " << pow(theta12.eval(),2)+pow(theta22.eval(),2)+pow(theta32.eval(),2)
//                << endl);
//        DebugOn("col 3 " << pow(theta13.eval(),2)+pow(theta23.eval(),2)+pow(theta33.eval(),2)
//                << endl);
//        double det=theta11.eval()*(theta22.eval()*theta33.eval()-theta32.eval()*theta23.eval())
//        -theta12.eval()*(theta21.eval()*theta33.eval()-theta31.eval()*theta23.eval())+theta13.eval()*(theta21.eval()*theta32.eval()-theta31.eval()*theta22.eval());
//
//
//        DebugOn("row 12 " << (theta11.eval()*theta21.eval())+(theta12.eval()*theta22.eval())+(theta13.eval()*theta23.eval())
//                << endl);
//        DebugOn("row 13 " << (theta11.eval()*theta31.eval())+(theta12.eval()*theta32.eval())+(theta13.eval()*theta33.eval())
//                << endl);
//        DebugOn("row 23 " << (theta21.eval()*theta31.eval())+(theta22.eval()*theta32.eval())+(theta23.eval()*theta33.eval())
//                << endl);
//
//        DebugOn("Determinant "<<det<<endl);
    //        Reg->print_solution();
    //    rot_trans[0]=theta11.eval();
    //    rot_trans[1]=theta12.eval();
    //    rot_trans[2]=theta13.eval();;
    //    rot_trans[3]=theta21.eval();
    //    rot_trans[4]=theta22.eval();
    //    rot_trans[5]=theta23.eval();
    //    rot_trans[6]=theta31.eval();
    //    rot_trans[7]=theta32.eval();
    //    rot_trans[8]=theta33.eval();
    //    rot_trans[9]=x_shift.eval();
    //    rot_trans[10]=y_shift.eval();
    //    rot_trans[11]=z_shift.eval();
    //
    //
    //
    //
    //    auto pitch_val = std::atan2(theta32.eval(), theta33.eval())*180/pi;
    //    auto roll_val = std::atan2(-1*theta31.eval(), std::sqrt(theta32.eval()*theta32.eval()+theta33.eval()*theta33.eval()))*180/pi;
    //    auto yaw_val = std::atan2(theta21.eval(),theta11.eval())*180/pi;
    //    DebugOn("Roll (degrees) = " << to_string_with_precision(roll_val,12) << endl);
    //    DebugOn("Pitch (degrees) = " << to_string_with_precision(pitch_val,12) << endl);
    //    DebugOn("Yaw (degrees) = " << to_string_with_precision(yaw_val,12) << endl);
    //    DebugOff("x shift = " << x_shift.eval() << endl);
    //    DebugOff("y shift = " << y_shift.eval() << endl);
    //    DebugOff("z shift = " << z_shift.eval() << endl);
    //
    
    return(Reg);
}
void initialize_model_from_parent(shared_ptr<Model<double>>parent, shared_ptr<Model<double>> child){
    for (auto it=parent->_vars.begin(); it!=parent->_vars.end(); it++){
        auto vname=it->second->_name;
        auto v = it->second.get();
        switch (v->get_intype()) {
            case double_: {
                auto var_par = (var<double>*)v;
                auto var_child = child->get_var<double>(vname);
                auto v_keys=var_par->get_keys();
                auto child_ind=*(var_child._indices);
                for(auto it_key=v_keys->begin(); it_key!=v_keys->end(); it_key++){
                    auto key = *it_key;
                    if(child_ind.has_key(key)){
                        var_child._val->at(var_child._indices->_keys_map->at(key)) = var_par->eval(key);
                    }
                }
                break;
            }
            case integer_: {
                auto var_par = (var<int>*)v;
                auto var_child = child->get_var<int>(vname);
                auto v_keys=var_par->get_keys();
                auto child_ind=*(var_child._indices);
                for(auto it_key=v_keys->begin(); it_key!=v_keys->end(); it_key++){
                    auto key = *it_key;
                    if(child_ind.has_key(key)){
                        var_child._val->at(var_child._indices->_keys_map->at(key)) = var_par->eval(key);
                    }
                }
                break;
            }
            case binary_: {
                auto var_par = (var<bool>*)v;
                auto var_child = child->get_var<bool>(vname);
                auto v_keys=var_par->get_keys();
                auto child_ind=*(var_child._indices);
                for(auto it_key=v_keys->begin(); it_key!=v_keys->end(); it_key++){
                    auto key = *it_key;
                    if(child_ind.has_key(key)){
                        var_child._val->at(var_child._indices->_keys_map->at(key)) = var_par->eval(key);
                    }
                }
                break;
            }
            case float_: {
                auto var_par = (var<float>*)v;
                auto var_child = child->get_var<float>(vname);
                auto v_keys=var_par->get_keys();
                auto child_ind=*(var_child._indices);
                for(auto it_key=v_keys->begin(); it_key!=v_keys->end(); it_key++){
                    auto key = *it_key;
                    if(child_ind.has_key(key)){
                        var_child._val->at(var_child._indices->_keys_map->at(key)) = var_par->eval(key);
                    }
                }
                break;
            }
            case long_: {
                auto var_par = (var<long double>*)v;
                auto var_child = child->get_var<long double>(vname);
                auto v_keys=var_par->get_keys();
                auto child_ind=*(var_child._indices);
                for(auto it_key=v_keys->begin(); it_key!=v_keys->end(); it_key++){
                    auto key = *it_key;
                    if(child_ind.has_key(key)){
                        var_child._val->at(var_child._indices->_keys_map->at(key)) = var_par->eval(key);
                    }
                }
                break;
            }
            case short_: {
                auto var_par = (var<short>*)v;
                auto var_child = child->get_var<short>(vname);
                auto v_keys=var_par->get_keys();
                auto child_ind=*(var_child._indices);
                for(auto it_key=v_keys->begin(); it_key!=v_keys->end(); it_key++){
                    auto key = *it_key;
                    if(child_ind.has_key(key)){
                        var_child._val->at(var_child._indices->_keys_map->at(key)) = var_par->eval(key);
                    }
                }
                break;
            }
        }
    }
}
double build_upperbound(vector<double>& solutionlb,vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, const indices& valid_cells, double new_roll_min, double new_roll_max, double new_pitch_min, double new_pitch_max, double new_yaw_min, double new_yaw_max, double new_shift_min_x, double new_shift_max_x, double new_shift_min_y, double new_shift_max_y, double new_shift_min_z, double new_shift_max_z, vector<double>& rot_trans){
    
    int nb_pairs = 0, min_nb_pairs = numeric_limits<int>::max(), max_nb_pairs = 0, av_nb_pairs = 0;
    size_t nm = point_cloud_model.size(), nd = point_cloud_data.size();
    vector<double> zeros = {0,0,0};
    param<> x1("x1"), x2("x2"), y1("y1"), y2("y2"), z1("z1"), z2("z2");
    int m = av_nb_pairs;
    string i_str, j_str;
    double xm_max = numeric_limits<double>::lowest(), ym_max = numeric_limits<double>::lowest(), zm_max = numeric_limits<double>::lowest();
    double xm_min = numeric_limits<double>::max(), ym_min = numeric_limits<double>::max(), zm_min = numeric_limits<double>::max();
    double res=numeric_limits<double>::max();
    param<> d1("d1"), d2("d2");
    for (auto i = 0; i<nd; i++) {
        i_str = to_string(i+1);
        x1.add_val(i_str,point_cloud_data.at(i).at(0));
        y1.add_val(i_str,point_cloud_data.at(i).at(1));
        z1.add_val(i_str,point_cloud_data.at(i).at(2));
        // d1.add_val(i_str, pow(point_cloud_data.at(i).at(0),2)+pow(point_cloud_data.at(i).at(1),2)+pow(point_cloud_data.at(i).at(2),2));
        
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
        // d2.add_val(j_str, pow(point_cloud_model.at(j).at(0),2)+pow(point_cloud_model.at(j).at(1),2)+pow(point_cloud_model.at(j).at(2),2));
    }
    
    
    indices Pairs("Pairs");
    int idx1 = 0;
    int idx2 = 0;
    indices N1("N1"),N2("N2");
    Debug("nd = " << nd << endl);
    Debug("nm = " << nm << endl);
    
    N1 = range(1,nd);
    N2 = range(1,nm);
    // auto cells = indices(N1,N2);
    auto cells=valid_cells;
    string name="lin_lb";
    
    auto Reg=make_shared<Model<>>(name);
    
    int count=0;
    param<> bin("bin");
    bin.in(cells);
    for(auto i=0;i<nd;i++){
        for(auto j=0;j<nm;j++){
            if(cells.has_key(to_string(i+1)+","+to_string(j+1))){
                bin.set_val(to_string(i+1)+","+to_string(j+1), round(solutionlb[count++]));
            }
        }
    }
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
    var<> x_shift("x_shift", new_shift_min_x, new_shift_max_x), y_shift("y_shift", new_shift_min_y, new_shift_max_y), z_shift("z_shift", new_shift_min_z, new_shift_max_z);
    Reg->add(x_shift.in(R(1)),y_shift.in(R(1)),z_shift.in(R(1)));
    
    
    DebugOn("Added " << cells.size() << " binary variables" << endl);
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
    
    DebugOff("Added " << cells.size() << " binary variables" << endl);
    

    
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
    // soc_12.add_to_callback();
    Reg->add(soc_12.in(range(0,0))<=0);
    
    Constraint<> soc_13("soc_13");
    soc_13 = pow(theta12-theta21,2)-(1-theta11-theta22+theta33)*(1+theta11+theta22+theta33);
    //soc_13.add_to_callback();
    Reg->add(soc_13.in(range(0,0))<=0);
    
    Constraint<> soc_14("soc_14");
    soc_14 = pow(theta23+theta32,2)-(1-theta11-theta22+theta33)*(1-theta11+theta22-theta33);
    //soc_14.add_to_callback();
    Reg->add(soc_14.in(range(0,0))<=0);
    
    Constraint<> soc_23("soc_23");
    soc_23 = pow(theta23-theta32,2)-(1+theta11-theta22-theta33)*(1+theta11+theta22+theta33);
    //soc_23.add_to_callback();
    Reg->add(soc_23.in(range(0,0))<=0);
    
    Constraint<> soc_24("soc_24");
    soc_24 = pow(theta12+theta21,2)-(1+theta11-theta22-theta33)*(1-theta11+theta22-theta33);
    //soc_24.add_to_callback();
    Reg->add(soc_24.in(range(0,0))<=0);
    
    Constraint<> soc_34("soc_34");
    soc_34 = pow(theta31-theta13,2)-(1+theta11+theta22+theta33)*(1-theta11+theta22-theta33);
    // soc_34.add_to_callback();
    Reg->add(soc_34.in(range(0,0))<=0);
    
    Constraint<> det_123("det_123");
    det_123+=(theta13+theta31)*((theta13+theta31)*(1+theta11+theta22+theta33)-(theta23-theta32)*(theta12-theta21));
    det_123-=(1-theta11-theta22+theta33)*((1+theta11-theta22-theta33)*(1+theta11+theta22+theta33)-pow(theta23-theta32,2));
    det_123-=(theta12-theta21)*((theta13+theta31)*(theta23-theta32)-(theta12-theta21)*(1+theta11-theta22-theta33));
    //det_123.add_to_callback();
    Reg->add(det_123.in(range(0,0))<=0);
    
    Constraint<> det_124("det_124");
    det_124+=(theta13+theta31)*((theta13+theta31)*(1-theta11+theta22-theta33)-(theta23+theta32)*(theta12+theta21));
    det_124-=(1-theta11-theta22+theta33)*((1+theta11-theta22-theta33)*(1-theta11+theta22-theta33)-pow(theta12+theta21,2));
    det_124-=(theta23+theta32)*((theta13+theta31)*(theta12+theta21)-(theta23+theta32)*(1+theta11-theta22-theta33));
    //det_124.add_to_callback();
    Reg->add(det_124.in(range(0,0))<=0);
    
    Constraint<> det_134("det_134");
    det_134+=(theta12-theta21)*((theta12-theta21)*(1-theta11+theta22-theta33)-(theta23+theta32)*(theta31-theta13));
    det_134-=(1-theta11-theta22+theta33)*((1+theta11+theta22+theta33)*(1-theta11+theta22-theta33)-pow(theta31-theta13,2));
    det_134-=(theta23+theta32)*((theta12-theta21)*(theta31-theta13)-(theta23+theta32)*(1+theta11+theta22+theta33));
    //det_134.add_to_callback();
    Reg->add(det_134.in(range(0,0))<=0);
    
    Constraint<> det_234("det_234");
    det_234+=(theta23-theta32)*((theta23-theta32)*(1-theta11+theta22-theta33)-(theta12+theta21)*(theta31-theta13));
    det_234-=(1+theta11-theta22-theta33)*((1+theta11+theta22+theta33)*(1-theta11+theta22-theta33)-pow(theta31-theta13,2));
    det_234-=(theta12+theta21)*((theta23-theta32)*(theta31-theta13)-(theta12+theta21)*(1+theta11+theta22+theta33));
    //det_234.add_to_callback();
    Reg->add(det_234.in(range(0,0))<=0);
    //
    Constraint<> row1("row1");
    row1 = pow(theta11,2)+pow(theta12,2)+pow(theta13,2);
    // row1.add_to_callback();
    Reg->add(row1.in(range(0,0))<=1);
    Constraint<> row2("row2");
    row2 = pow(theta21,2)+pow(theta22,2)+pow(theta23,2);
    //row2.add_to_callback();
    Reg->add(row2.in(range(0,0))<=1);
    Constraint<> row3("row3");
    row3 = pow(theta31,2)+pow(theta32,2)+pow(theta33,2);
    //row3.add_to_callback();
    Reg->add(row3.in(range(0,0))<=1);
    Constraint<> col1("col1");
    col1 = pow(theta11,2)+pow(theta21,2)+pow(theta31,2);
    // col1.add_to_callback();
    Reg->add(col1.in(range(0,0))<=1);
    Constraint<> col2("col2");
    col2 = pow(theta12,2)+pow(theta22,2)+pow(theta32,2);
    //col2.add_to_callback();
    Reg->add(col2.in(range(0,0))<=1);
    Constraint<> col3("col3");
    col3 = pow(theta13,2)+pow(theta23,2)+pow(theta33,2);
    //col3.add_to_callback();
    Reg->add(col3.in(range(0,0))<=1);
    
    vector<vector<vector<double>>> box;
    vector<vector<double>> box_i;
    vector<double> coord_i;
    coord_i.resize(3);
    
    param<> cL("cL");
    cL.in(cells);
    param<> cU("cU");
    cU.in(cells);
    
    shared_ptr<pair<double,double>> x1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> y1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> z1_bounds = make_shared<pair<double,double>>();
    vector<double> x_lb, x_ub, y_lb, y_ub, z_lb, z_ub;
    
    for(auto i=0;i<nd;i++){
        //  auto max_dist_i = get_max_dist(roll_min, roll_max, pitch_min, pitch_max, yaw_min, yaw_max, shift_min_x, shift_max_x, shift_min_y, shift_max_y, shift_min_z, shift_max_z, point_cloud_data[i], zeros, false);
        //dist_i.push_back(max_dist_i);
        
        
        // sphere_outer_sq.push_back(pow((d_root+shift_mag_max_root),2));
        
        x1_bounds->first = point_cloud_data.at(i)[0];
        x1_bounds->second = point_cloud_data.at(i)[0];
        y1_bounds->first = point_cloud_data.at(i)[1];
        y1_bounds->second = point_cloud_data.at(i)[1];
        z1_bounds->first = point_cloud_data.at(i)[2];
        z1_bounds->second = point_cloud_data.at(i)[2];
        auto x_range  = get_product_range(x1_bounds, theta11._range);
        auto y_range  = get_product_range(y1_bounds, theta12._range);
        auto z_range  = get_product_range(z1_bounds, theta13._range);
        x_lb.push_back(x_range->first + y_range->first + z_range->first + new_shift_min_x);
        x_ub.push_back(x_range->second + y_range->second + z_range->second+ new_shift_max_x);
        x_range  = get_product_range(x1_bounds, theta21._range);
        y_range  = get_product_range(y1_bounds, theta22._range);
        z_range  = get_product_range(z1_bounds, theta23._range);
        y_lb.push_back(x_range->first + y_range->first + z_range->first + new_shift_min_y);
        y_ub.push_back(x_range->second + y_range->second + z_range->second+ new_shift_max_y);
        x_range  = get_product_range(x1_bounds, theta31._range);
        y_range  = get_product_range(y1_bounds, theta32._range);
        z_range  = get_product_range(z1_bounds, theta33._range);
        z_lb.push_back(x_range->first + y_range->first + z_range->first + new_shift_min_z);
        z_ub.push_back(x_range->second + y_range->second + z_range->second+ new_shift_max_z);
        coord_i[0]=x_lb[i];
        coord_i[1]=y_lb[i];
        coord_i[2]=z_lb[i];
        box_i.push_back(coord_i);
        coord_i[0]=x_lb[i];
        coord_i[1]=y_ub[i];
        coord_i[2]=z_lb[i];
        box_i.push_back(coord_i);
        coord_i[0]=x_ub[i];
        coord_i[1]=y_ub[i];
        coord_i[2]=z_lb[i];
        box_i.push_back(coord_i);
        coord_i[0]=x_ub[i];
        coord_i[1]=y_lb[i];
        coord_i[2]=z_lb[i];
        box_i.push_back(coord_i);
        coord_i[0]=x_lb[i];
        coord_i[1]=y_lb[i];
        coord_i[2]=z_ub[i];
        box_i.push_back(coord_i);
        coord_i[0]=x_lb[i];
        coord_i[1]=y_ub[i];
        coord_i[2]=z_ub[i];
        box_i.push_back(coord_i);
        coord_i[0]=x_ub[i];
        coord_i[1]=y_ub[i];
        coord_i[2]=z_ub[i];
        box_i.push_back(coord_i);
        coord_i[0]=x_ub[i];
        coord_i[1]=y_lb[i];
        coord_i[2]=z_ub[i];
        box_i.push_back(coord_i);
        for(auto j=1;j<=nm;j++){
            if(cells.has_key(to_string(i+1)+","+to_string(j))){
                auto x=point_cloud_model.at(j-1)[0];
                auto y=point_cloud_model.at(j-1)[1];
                auto z=point_cloud_model.at(j-1)[2];
                double c_min=100, c_max=-999;
                for(auto k=0;k<8;k++){
                    auto inner_prod=box_i[k][0]*x+box_i[k][1]*y+box_i[k][2]*z;
                    if(inner_prod>=c_max){
                        c_max=inner_prod;
                    }
                    if(inner_prod<=c_min){
                        c_min=inner_prod;
                    }
                }
                
                cL.set_val(to_string(i+1)+","+to_string(j), c_min*bin.eval(to_string(i+1)+","+to_string(j)));
                cU.set_val(to_string(i+1)+","+to_string(j), c_max*bin.eval(to_string(i+1)+","+to_string(j)));
                
            }
        }
        box_i.clear();
    }
    var<> c("c");
    Reg->add(c.in(cells));
    
    double tx_max=std::max(pow(new_shift_min_x,2), pow(new_shift_max_x,2));
    double ty_max=std::max(pow(new_shift_min_y,2), pow(new_shift_max_y,2));
    double tz_max=std::max(pow(new_shift_min_z,2), pow(new_shift_max_z,2));
    var<> tx("tx", 0, tx_max), ty("ty", 0, ty_max), tz("tz", 0, tz_max);
    
    Reg->add(tx.in(R(1)),ty.in(R(1)),tz.in(R(1)));
    
    
    indices idsij = indices("idsij");
    idsij.add_empty_row();
    param<> dm("dm");
    dm.in(cells);
    param<> dm_pr("dm_pr");
    dm_pr.in(cells);
    param<> dd("dd");
    for(auto i=0;i<nd;i++){
        auto di=pow(point_cloud_data.at(i)[0],2)+pow(point_cloud_data.at(i)[1],2)+pow(point_cloud_data.at(i)[2],2);
        auto di_r=(sqrt(di)+sqrt(shift_mag_max));
        dd.add_val(to_string(i+1), pow(di_r,2));
        for(auto j=1;j<=nm;j++){
            if(cells.has_key(to_string(i+1)+","+to_string(j))){
                auto dmd=(pow(point_cloud_model.at(j-1)[0],2)+pow(point_cloud_model.at(j-1)[1],2)+pow(point_cloud_model.at(j-1)[2],2))*bin.eval(to_string(i+1)+","+to_string(j));
                DebugOff(bin.eval(to_string(i+1)+","+to_string(j)));
                auto dmd_r=sqrt(dmd)*di_r*(-1);
                dm.set_val(to_string(i+1)+","+to_string(j), dmd);
                dm_pr.set_val(to_string(i+1)+","+to_string(j), dmd_r);
                idsij.add_in_row(i, to_string(i+1)+","+to_string(j));
            }
        }
    }
    
    Constraint<> txsq("txsq");
    txsq = pow(x_shift,2) - tx;
    //txsq.add_to_callback();
    Reg->add(txsq.in(range(0,0))<=0);
    
    Constraint<> tysq("tysq");
    tysq = pow(y_shift,2) - ty;
    //tysq.add_to_callback();
    Reg->add(tysq.in(range(0,0))<=0);
    
    Constraint<> tzsq("tzsq");
    tzsq = pow(z_shift,2) - tz;
    // tzsq.add_to_callback();
    Reg->add(tzsq.in(range(0,0))<=0);
    
    
    
    
    auto idst = theta11.repeat_id(cells.size());
    Constraint<> def_c("def_c");
    def_c=c-((x1.from(cells)*theta11.in(idst) + y1.from(cells)*theta12.in(idst) + z1.from(cells)*theta13.in(idst)+x_shift)*x2.to(cells)+(x1.from(cells)*theta21.in(idst) + y1.from(cells)*theta22.in(idst) + z1.from(cells)*theta23.in(idst)+y_shift)*y2.to(cells)+(x1.from(cells)*theta31.in(idst) + y1.from(cells)*theta32.in(idst) + z1.from(cells)*theta33.in(idst)+z_shift)*z2.to(cells))*bin.in(cells)*2;
    Reg->add(def_c.in(cells)<=0);
    

    
    func<> obj = 0;
    
    
    obj+=(sum(x1*x1)+sum(y1*y1)+sum(z1*z1)+sum(dm)-sum(c)+nd*(tx+ty+tz));
    
    Reg->min(obj);
    
    // Reg->set_solution(solutionlb);
    //  Reg->round_and_fix();
    // Reg->print();
    solver<> S(Reg, ipopt);
    S.run(0, 1e-6);
    // Reg->print_solution();
    
    
    if(Reg->_status==0){
        res=Reg->get_obj_val();
        double det=theta11.eval()*(theta22.eval()*theta33.eval()-theta32.eval()*theta23.eval())
        -theta12.eval()*(theta21.eval()*theta33.eval()-theta31.eval()*theta23.eval())+theta13.eval()*(theta21.eval()*theta32.eval()-theta31.eval()*theta22.eval());
        
        
        DebugOff("row 12 " << (theta11.eval()*theta21.eval())+(theta12.eval()*theta22.eval())+(theta13.eval()*theta23.eval())
                << endl);
        DebugOff("row 13 " << (theta11.eval()*theta31.eval())+(theta12.eval()*theta32.eval())+(theta13.eval()*theta33.eval())
                << endl);
        DebugOff("row 23 " << (theta21.eval()*theta31.eval())+(theta22.eval()*theta32.eval())+(theta23.eval()*theta33.eval())
                << endl);
        
        DebugOff("Determinant "<<det<<endl);
        DebugOff("ub from ipopt "<<res<<endl);
    }
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
    return(res);
}
shared_ptr<Model<double>> build_linear_lowerbound(vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, const indices& valid_cells, double new_roll_min, double new_roll_max, double new_pitch_min, double new_pitch_max, double new_yaw_min, double new_yaw_max, double new_shift_min_x, double new_shift_max_x, double new_shift_min_y, double new_shift_max_y, double new_shift_min_z, double new_shift_max_z, vector<double>& rot_trans, bool separate, const vector<pair<pair<int,int>,pair<int,int>>>& incompatibles,  param<>& norm_x,  param<>& norm_y,  param<>& norm_z,  param<>& intercept,const vector<int>& init_matching, const vector<double>& error_per_point, bool relax_inits){
    
    int nb_pairs = 0, min_nb_pairs = numeric_limits<int>::max(), max_nb_pairs = 0, av_nb_pairs = 0;
    size_t nm = point_cloud_model.size(), nd = point_cloud_data.size();
    vector<double> zeros = {0,0,0};
    param<> x1("x1"), x2("x2"), y1("y1"), y2("y2"), z1("z1"), z2("z2");
    int m = av_nb_pairs;
    string i_str, j_str;
    double xm_max = numeric_limits<double>::lowest(), ym_max = numeric_limits<double>::lowest(), zm_max = numeric_limits<double>::lowest();
    double xm_min = numeric_limits<double>::max(), ym_min = numeric_limits<double>::max(), zm_min = numeric_limits<double>::max();
    
    
    for (auto i = 0; i<nd; i++) {
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
    Debug("nd = " << nd << endl);
    Debug("nm = " << nm << endl);
    
    N1 = range(1,nd);
    N2 = range(1,nm);
    // auto cells = indices(N1,N2);
    auto cells=valid_cells;
    string name="lin_lb";
    
    auto Reg=make_shared<Model<>>(name);
    
    var<int> bin("bin",0,1);
    Reg->add(bin.in(cells));
    
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
    var<> x_shift("x_shift", new_shift_min_x, new_shift_max_x), y_shift("y_shift", new_shift_min_y, new_shift_max_y), z_shift("z_shift", new_shift_min_z, new_shift_max_z);
    Reg->add(x_shift.in(R(1)),y_shift.in(R(1)),z_shift.in(R(1)));
    
    
    DebugOn("Added " << cells.size() << " binary variables" << endl);
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
    
    DebugOn("Added " << cells.size() << " binary variables" << endl);
    
    
    
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
    row1 = pow(theta11,2)+pow(theta12,2)+pow(theta13,2);
    // row1.add_to_callback();
    Reg->add(row1.in(range(0,0))==1);
    Constraint<> row2("row2");
    row2 = pow(theta21,2)+pow(theta22,2)+pow(theta23,2);
    //  row2.add_to_callback();
    Reg->add(row2.in(range(0,0))==1);
    Constraint<> row3("row3");
    row3 = pow(theta31,2)+pow(theta32,2)+pow(theta33,2);
    // row3.add_to_callback();
    Reg->add(row3.in(range(0,0))==1);
    Constraint<> col1("col1");
    col1 = pow(theta11,2)+pow(theta21,2)+pow(theta31,2);
    // col1.add_to_callback();
    Reg->add(col1.in(range(0,0))<=1);
    Constraint<> col2("col2");
    col2 = pow(theta12,2)+pow(theta22,2)+pow(theta32,2);
    //col2.add_to_callback();
    Reg->add(col2.in(range(0,0))<=1);
    Constraint<> col3("col3");
    col3 = pow(theta13,2)+pow(theta23,2)+pow(theta33,2);
    //col3.add_to_callback();
    Reg->add(col3.in(range(0,0))<=1);
    
    
    Constraint<> OneBin("OneBin");
    OneBin = bin.in_matrix(1, 1);
    Reg->add(OneBin.in(N1)==1);
    
    
    param<> x_new_lb("x_new_lb");
    //x_new_lb.in(N1);
    param<> x_new_ub("x_new_ub");
    //x_new_ub.in(N1);
    param<> y_new_lb("y_new_lb");
    //y_new_lb.in(N1);
    param<> y_new_ub("y_new_ub");
    //y_new_ub.in(N1);
    param<> z_new_lb("z_new_lb");
    //z_new_lb.in(N1);
    param<> z_new_ub("z_new_ub");
    //z_new_ub.in(N1);
    param<> cx_L("cx_L");
    param<> cx_U("cx_U");
    param<> cy_L("cy_L");
    param<> cy_U("cy_U");
    param<> cz_L("cz_L");
    param<> cz_U("cz_U");
    shared_ptr<pair<double,double>> new_x1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> new_y1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> new_z1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> x1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> y1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> z1_bounds = make_shared<pair<double,double>>();
    for(auto i=0;i<nd;i++){
        x1_bounds->first = point_cloud_data.at(i)[0];
        x1_bounds->second = point_cloud_data.at(i)[0];
        y1_bounds->first = point_cloud_data.at(i)[1];
        y1_bounds->second = point_cloud_data.at(i)[1];
        z1_bounds->first = point_cloud_data.at(i)[2];
        z1_bounds->second = point_cloud_data.at(i)[2];
        auto x_range  = get_product_range(x1_bounds, theta11._range);
        auto y_range  = get_product_range(y1_bounds, theta12._range);
        auto z_range  = get_product_range(z1_bounds, theta13._range);
        *new_x1_bounds = {x_range->first + y_range->first + z_range->first + new_shift_min_x,
            x_range->second + y_range->second + z_range->second+ new_shift_max_x};
        x_new_lb.add_val(to_string(i+1), new_x1_bounds->first);
        x_new_ub.add_val(to_string(i+1), new_x1_bounds->second);
        auto x=std::max(pow(new_x1_bounds->first, 2),pow(new_x1_bounds->second, 2))+std::max(xm_min, xm_max);
        cx_L.add_val(to_string(i+1), x*(-0.5));
        cx_U.add_val(to_string(i+1), x*(0.5));
        x_range  = get_product_range(x1_bounds, theta21._range);
        y_range  = get_product_range(y1_bounds, theta22._range);
        z_range  = get_product_range(z1_bounds, theta23._range);
        *new_y1_bounds = {x_range->first + y_range->first + z_range->first + new_shift_min_y,
            x_range->second + y_range->second + z_range->second+ new_shift_max_y};
        y_new_lb.add_val(to_string(i+1), new_y1_bounds->first);
        y_new_ub.add_val(to_string(i+1), new_y1_bounds->second);
        auto y=std::max(pow(new_y1_bounds->first, 2),pow(new_y1_bounds->second, 2))+std::max(ym_min, ym_max);
        cy_L.add_val(to_string(i+1), y*(-0.5));
        cy_U.add_val(to_string(i+1), y*(0.5));
        x_range  = get_product_range(x1_bounds, theta31._range);
        y_range  = get_product_range(y1_bounds, theta32._range);
        z_range  = get_product_range(z1_bounds, theta33._range);
        *new_z1_bounds = {x_range->first + y_range->first + z_range->first + new_shift_min_z,
            x_range->second + y_range->second + z_range->second+ new_shift_max_z};
        z_new_lb.add_val(to_string(i+1), new_z1_bounds->first);
        z_new_ub.add_val(to_string(i+1), new_z1_bounds->second);
        auto z=std::max(pow(new_z1_bounds->first, 2),pow(new_z1_bounds->second, 2))+std::max(zm_min, zm_max);
        cz_L.add_val(to_string(i+1), z*(-0.5));
        cz_U.add_val(to_string(i+1), z*(0.5));
    }
    var<> new_x1("new_x1", x_new_lb, x_new_ub),new_y1("new_y1", y_new_lb, y_new_ub),new_z1("new_z1", z_new_lb, z_new_ub);
    var<> new_xm("new_xm", xm_min, xm_max), new_ym("new_ym", ym_min, ym_max), new_zm("new_zm", zm_min, zm_max);
    //var<> new_xm("new_xm"), new_ym("new_ym"), new_zm("new_zm");
    var<> delta("delta", pos_);
    
    Reg->add(new_x1.in(N1));
    Reg->add(new_y1.in(N1));
    Reg->add(new_z1.in(N1));
    Reg->add(new_xm.in(N1));
    Reg->add(new_ym.in(N1));
    Reg->add(new_zm.in(N1));
    //  Reg->add(delta.in(N1));
    Reg->add(delta.in(R(1)));
    
    
    double tx_max=std::max(pow(new_shift_min_x,2), pow(new_shift_max_x,2));
    double ty_max=std::max(pow(new_shift_min_y,2), pow(new_shift_max_y,2));
    double tz_max=std::max(pow(new_shift_min_z,2), pow(new_shift_max_z,2));
    var<> tx("tx", 0, tx_max), ty("ty", 0, ty_max), tz("tz", 0, tz_max);
    
    // Reg->add(tx.in(R(1)),ty.in(R(1)),tz.in(R(1)));
    param<> dm_root("dm_root");
    param<> x2_sq("x2_sq");
    param<> y2_sq("y2_sq");
    param<> z2_sq("z2_sq");
    for(auto i=0;i<nm;i++){
        auto x=pow(point_cloud_model.at(i)[0],2)*(-1);
        auto y=pow(point_cloud_model.at(i)[1],2)*(-1);
        auto z=pow(point_cloud_model.at(i)[2],2)*(-1);
        x2_sq.add_val(to_string(i+1), x);
        y2_sq.add_val(to_string(i+1), y);
        z2_sq.add_val(to_string(i+1), z);
        dm_root.add_val(to_string(i+1), (x+y+z)*(-1));
    }
    
    indices idsij = indices("idsij");
    idsij.add_empty_row();
    param<> dm("dm");
    dm.in(cells);
    param<> dm_pr("dm_pr");
    dm_pr.in(cells);
    param<> dd("dd");
    for(auto i=0;i<nd;i++){
        auto di=pow(point_cloud_data.at(i)[0],2)+pow(point_cloud_data.at(i)[1],2)+pow(point_cloud_data.at(i)[2],2);
        auto di_r=(sqrt(di)+sqrt(shift_mag_max));
        dd.add_val(to_string(i+1), pow(di_r,2));
        for(auto j=1;j<=nm;j++){
            if(cells.has_key(to_string(i+1)+","+to_string(j))){
                auto dmd=pow(point_cloud_model.at(j-1)[0],2)+pow(point_cloud_model.at(j-1)[1],2)+pow(point_cloud_model.at(j-1)[2],2);
                auto dmd_r=sqrt(dmd)*di_r*(-1);
                dm.set_val(to_string(i+1)+","+to_string(j), dmd);
                dm_pr.set_val(to_string(i+1)+","+to_string(j), dmd_r);
                idsij.add_in_row(i, to_string(i+1)+","+to_string(j));
            }
        }
    }
    
    Constraint<> txsq("txsq");
    txsq = pow(x_shift,2) - tx;
    txsq.add_to_callback();
    // Reg->add(txsq.in(range(0,0))<=0);
    
    Constraint<> tysq("tysq");
    tysq = pow(y_shift,2) - ty;
    tysq.add_to_callback();
    //Reg->add(tysq.in(range(0,0))<=0);
    
    Constraint<> tzsq("tzsq");
    tzsq = pow(z_shift,2) - tz;
    tzsq.add_to_callback();
    //Reg->add(tzsq.in(range(0,0))<=0);
    
    indices ids = indices("in_x");
    ids.add_empty_row();
    for(auto i=0;i<nd;i++){
        for(auto j=1;j<=nm;j++){
            if(cells.has_key(to_string(i+1)+","+to_string(j)))
                ids.add_in_row(i, to_string(j));
        }
    }
    
    auto ids1 = theta11.repeat_id(N1.size());
    Constraint<> Def_newxm("Def_newxm");
    Def_newxm = new_xm-product(x2.in(ids),bin.in_matrix(1, 1));
    //Def_newxm += x1.in(N1)*theta11.in(ids1) + y1.in(N1)*theta12.in(ids1) + z1.in(N1)*theta13.in(ids1)+x_shift;
    Reg->add(Def_newxm.in(N1)==0);
    
    Constraint<> Def_newym("Def_newym");
    Def_newym = new_ym-product(y2.in(ids),bin.in_matrix(1, 1));
    // Def_newym += x1.in(N1)*theta21.in(ids1) + y1.in(N1)*theta22.in(ids1) + z1.in(N1)*theta23.in(ids1)+y_shift;
    Reg->add(Def_newym.in(N1)==0);
    
    Constraint<> Def_newzm("Def_newzm");
    Def_newzm = new_zm-product(z2.in(ids),bin.in_matrix(1, 1));
    // Def_newzm += x1.in(N1)*theta31.in(ids1) + y1.in(N1)*theta32.in(ids1) + z1.in(N1)*theta33.in(ids1)+z_shift;
    Reg->add(Def_newzm.in(N1)==0);
    
    
    
    //auto ids1 = theta11.repeat_id(N1.size());
    Constraint<> x_rot1("x_rot1");
    x_rot1 += new_x1 -x_shift;
    x_rot1 -= x1.in(N1)*theta11.in(ids1) + y1.in(N1)*theta12.in(ids1) + z1.in(N1)*theta13.in(ids1);
    Reg->add(x_rot1.in(N1)==0);
    
    Constraint<> y_rot1("y_rot1");
    y_rot1 += new_y1 - y_shift;
    y_rot1 -= x1.in(N1)*theta21.in(ids1) + y1.in(N1)*theta22.in(ids1) + z1.in(N1)*theta23.in(ids1);
    Reg->add(y_rot1.in(N1)==0);
    
    Constraint<> z_rot1("z_rot1");
    z_rot1 += new_z1 -z_shift;
    z_rot1 -= x1.in(N1)*theta31.in(ids1) + y1.in(N1)*theta32.in(ids1) + z1.in(N1)*theta33.in(ids1);
    Reg->add(z_rot1.in(N1)==0);
    
    var<> cx("cx",cx_L, cx_U), cy("cy",cy_L, cy_U), cz("cz",cz_L,cz_U);
    Reg->add(cx.in(N1), cy.in(N1), cz.in(N1));
    
    Constraint<> def_cx("def_cx");
    def_cx=cx-new_x1*new_xm;
    Reg->add(def_cx.in(N1)==0);
    
    Constraint<> def_cy("def_cy");
    def_cy=cy-new_y1*new_ym;
    Reg->add(def_cy.in(N1)==0);
    
    Constraint<> def_cz("def_cz");
    def_cz=cz-new_z1*new_zm;
    Reg->add(def_cz.in(N1)==0);
    
    Constraint<> cx_l("cx_l");
    cx_l=product(x2_sq.in(ids), bin.in_matrix(1, 1))-new_x1*new_x1-2*cx;
    //Reg->add(cx_l.in(N1)<=0);
    
    Constraint<> cy_l("cy_l");
    cy_l=product(y2_sq.in(ids), bin.in_matrix(1, 1))-new_y1*new_y1-2*cy;
    //Reg->add(cy_l.in(N1)<=0);
    
    Constraint<> cz_l("cz_l");
    cz_l=product(z2_sq.in(ids), bin.in_matrix(1, 1))-new_z1*new_z1-2*cz;
    // Reg->add(cz_l.in(N1)<=0);
    
    Constraint<> c_l("c_l");
    c_l=dd*(-1)-product(dm_root.in(ids), bin.in_matrix(1, 1))-2*(cx+cy+cz);
    //Reg->add(c_l.in(N1)<=0);
    param<> two("2");
    two.in(N1);
    two = 2;
    Constraint<> def_delta("def_delta");
    def_delta=(sum(new_x1*new_x1)+sum(new_y1*new_y1)+sum(new_z1*new_z1)-two.tr()*cx-two.tr()*cy-two.tr()*cz+sum(new_xm*new_xm)+sum(new_ym*new_ym)+sum(new_zm*new_zm))-delta;
    Reg->add(def_delta<=0);
    /* Objective function */
    
    func<> obj = 0;
    
    
    //obj+=(sum(x1*x1)+sum(y1*y1)+sum(z1*z1)+nd*(tx+ty+tz)-two.tr()*cx-two.tr()*cy-two.tr()*cz+sum(dm*bin));
    //obj+=(sum(new_x1*new_x1)+sum(new_y1*new_y1)+sum(new_z1*new_z1)-two.tr()*cx-two.tr()*cy-two.tr()*cz+sum(new_xm*new_xm)+sum(new_ym*new_ym)+sum(new_zm*new_zm));
    obj+=delta;
    Reg->min(obj);
    
    Reg->print();
    solver<> S(Reg,gurobi);
    S.use_callback();
    S.run();
    
    
    return(Reg);
}
shared_ptr<Model<double>> build_polyhedral(vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, const indices& valid_cells, double new_roll_min, double new_roll_max, double new_pitch_min, double new_pitch_max, double new_yaw_min, double new_yaw_max, double new_shift_min_x, double new_shift_max_x, double new_shift_min_y, double new_shift_max_y, double new_shift_min_z, double new_shift_max_z, vector<double>& rot_trans, bool separate, const vector<pair<pair<int,int>,pair<int,int>>>& incompatibles,  param<>& norm_x,  param<>& norm_y,  param<>& norm_z,  param<>& intercept,const vector<int>& init_matching, const vector<double>& error_per_point, const vector<vector<vector<double>>>& model_voronoi_normals, const vector<vector<double>>& model_face_intercept, bool relax_inits){
    
    int nb_pairs = 0, min_nb_pairs = numeric_limits<int>::max(), max_nb_pairs = 0, av_nb_pairs = 0;
    size_t nm = point_cloud_model.size(), nd = point_cloud_data.size();
    vector<double> zeros = {0,0,0};
    param<> x1("x1"), x2("x2"), y1("y1"), y2("y2"), z1("z1"), z2("z2");
    int m = av_nb_pairs;
    string i_str, j_str;
    double xm_max = numeric_limits<double>::lowest(), ym_max = numeric_limits<double>::lowest(), zm_max = numeric_limits<double>::lowest();
    double xm_min = numeric_limits<double>::max(), ym_min = numeric_limits<double>::max(), zm_min = numeric_limits<double>::max();
    
    
    for (auto i = 0; i<nd; i++) {
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
    Debug("nd = " << nd << endl);
    Debug("nm = " << nm << endl);
    
    N1 = range(1,nd);
    N2 = range(1,nm);
    // auto cells = indices(N1,N2);
    auto cells=valid_cells;
    string name="lin_lb";
    
    auto Reg=make_shared<Model<>>(name);
    
    var<int> bin("bin",0,1);
    Reg->add(bin.in(cells));
    
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
    var<> x_shift("x_shift", new_shift_min_x, new_shift_max_x), y_shift("y_shift", new_shift_min_y, new_shift_max_y), z_shift("z_shift",new_shift_min_z,new_shift_max_z);
    Reg->add(x_shift.in(R(1)),y_shift.in(R(1)),z_shift.in(R(1)));
    
    
    DebugOn("Added " << cells.size() << " binary variables" << endl);
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
    
//    Reg->add(theta11.in(R(1)),theta12.in(R(1)),theta13.in(R(1)));
//    Reg->add(theta21.in(R(1)),theta22.in(R(1)),theta23.in(R(1)));
//    Reg->add(theta31.in(R(1)),theta32.in(R(1)),theta33.in(R(1)));
    
    DebugOn("Added " << cells.size() << " binary variables" << endl);
    
    
    
//    Constraint<> diag_1("diag_1");
//    diag_1=theta11+theta22-theta33-1;
//    Reg->add(diag_1.in(range(0,0))<=0);
//    Constraint<> diag_2("diag_2");
//    diag_2=theta22+theta33-(1+theta11);
//    Reg->add(diag_2.in(range(0,0))<=0);
//    Constraint<> diag_3("diag_3");
//    diag_3-=1+theta11+theta22+theta33;
//    Reg->add(diag_3.in(range(0,0))<=0);
//    Constraint<> diag_4("diag_4");
//    diag_4=theta11-theta22+theta33-1;
//    Reg->add(diag_4.in(range(0,0))<=0);
//
//            Constraint<> soc_12("soc_12");
//            soc_12 = pow(theta13+theta31,2)-(1-theta11-theta22+theta33)*(1+theta11-theta22-theta33);
//           // soc_12.add_to_callback();
//            Reg->add(soc_12.in(range(0,0))<=0);
//
//            Constraint<> soc_13("soc_13");
//            soc_13 = pow(theta12-theta21,2)-(1-theta11-theta22+theta33)*(1+theta11+theta22+theta33);
//            //soc_13.add_to_callback();
//            Reg->add(soc_13.in(range(0,0))<=0);
//
//            Constraint<> soc_14("soc_14");
//            soc_14 = pow(theta23+theta32,2)-(1-theta11-theta22+theta33)*(1-theta11+theta22-theta33);
//            //soc_14.add_to_callback();
//            Reg->add(soc_14.in(range(0,0))<=0);
//
//            Constraint<> soc_23("soc_23");
//            soc_23 = pow(theta23-theta32,2)-(1+theta11-theta22-theta33)*(1+theta11+theta22+theta33);
//            //soc_23.add_to_callback();
//            Reg->add(soc_23.in(range(0,0))<=0);
//
//            Constraint<> soc_24("soc_24");
//            soc_24 = pow(theta12+theta21,2)-(1+theta11-theta22-theta33)*(1-theta11+theta22-theta33);
//            //soc_24.add_to_callback();
//            Reg->add(soc_24.in(range(0,0))<=0);
//
//            Constraint<> soc_34("soc_34");
//            soc_34 = pow(theta31-theta13,2)-(1+theta11+theta22+theta33)*(1-theta11+theta22-theta33);
//           // soc_34.add_to_callback();
//            Reg->add(soc_34.in(range(0,0))<=0);
//    //
//            Constraint<> det_123("det_123");
//            det_123+=(theta13+theta31)*((theta13+theta31)*(1+theta11+theta22+theta33)-(theta23-theta32)*(theta12-theta21));
//            det_123-=(1-theta11-theta22+theta33)*((1+theta11-theta22-theta33)*(1+theta11+theta22+theta33)-pow(theta23-theta32,2));
//            det_123-=(theta12-theta21)*((theta13+theta31)*(theta23-theta32)-(theta12-theta21)*(1+theta11-theta22-theta33));
//            det_123.add_to_callback();
//            Reg->add(det_123.in(range(0,0))<=0);
//
//            Constraint<> det_124("det_124");
//            det_124+=(theta13+theta31)*((theta13+theta31)*(1-theta11+theta22-theta33)-(theta23+theta32)*(theta12+theta21));
//            det_124-=(1-theta11-theta22+theta33)*((1+theta11-theta22-theta33)*(1-theta11+theta22-theta33)-pow(theta12+theta21,2));
//            det_124-=(theta23+theta32)*((theta13+theta31)*(theta12+theta21)-(theta23+theta32)*(1+theta11-theta22-theta33));
//            det_124.add_to_callback();
//            Reg->add(det_124.in(range(0,0))<=0);
//
//            Constraint<> det_134("det_134");
//            det_134+=(theta12-theta21)*((theta12-theta21)*(1-theta11+theta22-theta33)-(theta23+theta32)*(theta31-theta13));
//            det_134-=(1-theta11-theta22+theta33)*((1+theta11+theta22+theta33)*(1-theta11+theta22-theta33)-pow(theta31-theta13,2));
//            det_134-=(theta23+theta32)*((theta12-theta21)*(theta31-theta13)-(theta23+theta32)*(1+theta11+theta22+theta33));
//            det_134.add_to_callback();
//            Reg->add(det_134.in(range(0,0))<=0);
//
//            Constraint<> det_234("det_234");
//            det_234+=(theta23-theta32)*((theta23-theta32)*(1-theta11+theta22-theta33)-(theta12+theta21)*(theta31-theta13));
//            det_234-=(1+theta11-theta22-theta33)*((1+theta11+theta22+theta33)*(1-theta11+theta22-theta33)-pow(theta31-theta13,2));
//            det_234-=(theta12+theta21)*((theta23-theta32)*(theta31-theta13)-(theta12+theta21)*(1+theta11+theta22+theta33));
//            det_234.add_to_callback();
//            Reg->add(det_234.in(range(0,0))<=0);
//    //
//    Constraint<> row1("row1");
//    row1 = pow(theta11,2)+pow(theta12,2)+pow(theta13,2);
//    // row1.add_to_callback();
//    Reg->add(row1.in(range(0,0))<=1);
//    Constraint<> row2("row2");
//    row2 = pow(theta21,2)+pow(theta22,2)+pow(theta23,2);
//    //row2.add_to_callback();
//    Reg->add(row2.in(range(0,0))<=1);
//    Constraint<> row3("row3");
//    row3 = pow(theta31,2)+pow(theta32,2)+pow(theta33,2);
//    //row3.add_to_callback();
//    Reg->add(row3.in(range(0,0))<=1);
//    Constraint<> col1("col1");
//    col1 = pow(theta11,2)+pow(theta21,2)+pow(theta31,2);
//    // col1.add_to_callback();
//    Reg->add(col1.in(range(0,0))<=1);
//    Constraint<> col2("col2");
//    col2 = pow(theta12,2)+pow(theta22,2)+pow(theta32,2);
//    //col2.add_to_callback();
//    Reg->add(col2.in(range(0,0))<=1);
//    Constraint<> col3("col3");
//    col3 = pow(theta13,2)+pow(theta23,2)+pow(theta33,2);
//    //col3.add_to_callback();
//    Reg->add(col3.in(range(0,0))<=1);
    
    
    Constraint<> OneBin("OneBin");
    OneBin = bin.in_matrix(1, 1);
    Reg->add(OneBin.in(N1)==1);

    param<> xL("xL");
    xL.in(cells);
    param<> xU("xU");
    xU.in(cells);
    param<> yL("yL");
    yL.in(cells);
    param<> yU("yU");
    yU.in(cells);
    param<> zL("zL");
    zL.in(cells);
    param<> zU("zU");
    zU.in(cells);
    param<> cL("cL");
    cL.in(cells);
    param<> cU("cU");
    cU.in(cells);
    param<> xxU("xxU"),xxL("xxL"),yyU("yyU"),yyL("yyL"),zzU("zzU"),zzL("zzL");
    
    for(auto j=1;j<=nm;j++){
        auto x=point_cloud_model.at(j-1)[0];
        auto y=point_cloud_model.at(j-1)[1];
        auto z=point_cloud_model.at(j-1)[2];
        xxL.add_val(to_string(j),std::min(new_shift_min_x*x, new_shift_max_x*x));
        xxU.add_val(to_string(j),std::max(new_shift_min_x*x, new_shift_max_x*x));
        yyL.add_val(to_string(j),std::min(new_shift_min_y*y, new_shift_max_y*y));
        yyU.add_val(to_string(j),std::max(new_shift_min_y*y, new_shift_max_y*y));
        zzL.add_val(to_string(j),std::min(new_shift_min_z*z, new_shift_max_z*z));
        zzU.add_val(to_string(j),std::max(new_shift_min_z*z, new_shift_max_z*z));
    }
    
    vector<vector<vector<double>>> box;
    vector<vector<double>> box_i;
    vector<double> coord_i;
    coord_i.resize(3);
    
    shared_ptr<pair<double,double>> x1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> y1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> z1_bounds = make_shared<pair<double,double>>();
    vector<double> x_lb, x_ub, y_lb, y_ub, z_lb, z_ub;
    indices vnorm("vnorm");
    vnorm.add_empty_row();
    indices vface("vface");
    vface.add_empty_row();
    int q=0;
    int nfacets=0;
    for(auto i=0;i<nd;i++){
        for(auto j=1;j<=nm;j++){
            if(cells.has_key(to_string(i+1)+","+to_string(j))){
                nfacets+=model_face_intercept.at(j-1).size();
            }
        }
    }
    indices facets("facets");
    facets=range(1, nfacets);
    param<> la("la"), lb("lb"),lc("lc"),ld("ld");
    la.in(facets);lb.in(facets);lc.in(facets);ld.in(facets);
    for(auto i=0;i<nd;i++){
        x1_bounds->first = point_cloud_data.at(i)[0];
        x1_bounds->second = point_cloud_data.at(i)[0];
        y1_bounds->first = point_cloud_data.at(i)[1];
        y1_bounds->second = point_cloud_data.at(i)[1];
        z1_bounds->first = point_cloud_data.at(i)[2];
        z1_bounds->second = point_cloud_data.at(i)[2];
        auto x_range  = get_product_range(x1_bounds, theta11._range);
        auto y_range  = get_product_range(y1_bounds, theta12._range);
        auto z_range  = get_product_range(z1_bounds, theta13._range);
        x_lb.push_back(x_range->first + y_range->first + z_range->first);
        x_ub.push_back(x_range->second + y_range->second + z_range->second);
        x_range  = get_product_range(x1_bounds, theta21._range);
        y_range  = get_product_range(y1_bounds, theta22._range);
        z_range  = get_product_range(z1_bounds, theta23._range);
        y_lb.push_back(x_range->first + y_range->first + z_range->first);
        y_ub.push_back(x_range->second + y_range->second + z_range->second);
        x_range  = get_product_range(x1_bounds, theta31._range);
        y_range  = get_product_range(y1_bounds, theta32._range);
        z_range  = get_product_range(z1_bounds, theta33._range);
        z_lb.push_back(x_range->first + y_range->first + z_range->first );
        z_ub.push_back(x_range->second + y_range->second + z_range->second);
        coord_i[0]=x_lb[i];
        coord_i[1]=y_lb[i];
        coord_i[2]=z_lb[i];
        box_i.push_back(coord_i);
        coord_i[0]=x_lb[i];
        coord_i[1]=y_ub[i];
        coord_i[2]=z_lb[i];
        box_i.push_back(coord_i);
        coord_i[0]=x_ub[i];
        coord_i[1]=y_ub[i];
        coord_i[2]=z_lb[i];
        box_i.push_back(coord_i);
        coord_i[0]=x_ub[i];
        coord_i[1]=y_lb[i];
        coord_i[2]=z_lb[i];
        box_i.push_back(coord_i);
        coord_i[0]=x_lb[i];
        coord_i[1]=y_lb[i];
        coord_i[2]=z_ub[i];
        box_i.push_back(coord_i);
        coord_i[0]=x_lb[i];
        coord_i[1]=y_ub[i];
        coord_i[2]=z_ub[i];
        box_i.push_back(coord_i);
        coord_i[0]=x_ub[i];
        coord_i[1]=y_ub[i];
        coord_i[2]=z_ub[i];
        box_i.push_back(coord_i);
        coord_i[0]=x_ub[i];
        coord_i[1]=y_lb[i];
        coord_i[2]=z_ub[i];
        box_i.push_back(coord_i);
        auto di=pow(point_cloud_data.at(i)[0],2)+pow(point_cloud_data.at(i)[1],2)+pow(point_cloud_data.at(i)[2],2);
        auto di_r=(sqrt(di));
        for(auto j=1;j<=nm;j++){
            if(cells.has_key(to_string(i+1)+","+to_string(j))){
                auto x=point_cloud_model.at(j-1)[0];
                auto y=point_cloud_model.at(j-1)[1];
                auto z=point_cloud_model.at(j-1)[2];
                double c_min=100, c_max=-999;
                for(auto k=0;k<8;k++){
                    auto inner_prod=box_i[k][0]*x+box_i[k][1]*y+box_i[k][2]*z;
                    if(inner_prod>=c_max){
                        c_max=inner_prod;
                    }
                    if(inner_prod<=c_min){
                        c_min=inner_prod;
                    }
                }
                auto normals=model_voronoi_normals.at(j-1);
                auto intercepts=model_face_intercept.at(j-1);
                for(auto k=0;k<intercepts.size();k++){
                    la.set_val(q, normals[k][0]);
                    lb.set_val(q, normals[k][1]);
                    lc.set_val(q, normals[k][2]);
                    ld.set_val(q, intercepts[k]);
                    vnorm.add_in_row(q, to_string(i+1)+","+to_string(j));
                    vface.add_in_row(q, to_string(q+1));
                    q++;
                }
                auto dmd=pow(point_cloud_model.at(j-1)[0],2)+pow(point_cloud_model.at(j-1)[1],2)+pow(point_cloud_model.at(j-1)[2],2);
                auto cL_val=std::max(sqrt(pow(x,2)+pow(y,2)+pow(z,2))*di_r*(-1), c_min);
                auto cU_val=std::min(sqrt(pow(x,2)+pow(y,2)+pow(z,2))*di_r, c_max);
                cL.set_val(to_string(i+1)+","+to_string(j), cL_val);
                cU.set_val(to_string(i+1)+","+to_string(j), dmd-cU_val*2);
                xL.set_val(to_string(i+1)+","+to_string(j), x_lb[i]);
                xU.set_val(to_string(i+1)+","+to_string(j), x_ub[i]);
                yL.set_val(to_string(i+1)+","+to_string(j), y_lb[i]);
                yU.set_val(to_string(i+1)+","+to_string(j), y_ub[i]);
                zL.set_val(to_string(i+1)+","+to_string(j), z_lb[i]);
                zU.set_val(to_string(i+1)+","+to_string(j), z_ub[i]);
            }
        }
        box_i.clear();
    }
    //var<> c("c",cL, cU);
    //Reg->add(c.in(cells));
    
    double tx_max=std::max(pow(new_shift_min_x,2), pow(new_shift_max_x,2));
    double ty_max=std::max(pow(new_shift_min_y,2), pow(new_shift_max_y,2));
    double tz_max=std::max(pow(new_shift_min_z,2), pow(new_shift_max_z,2));
    //var<> tx("tx", pos_), ty("ty",pos_), tz("tz", pos_);
    
    //Reg->add(tx.in(R(1)),ty.in(R(1)),tz.in(R(1)));
    param<> dm_root("dm_root");
    param<> x2_sq("x2_sq");
    param<> y2_sq("y2_sq");
    param<> z2_sq("z2_sq");
    for(auto i=0;i<nm;i++){
        auto x=pow(point_cloud_model.at(i)[0],2)*(-1);
        auto y=pow(point_cloud_model.at(i)[1],2)*(-1);
        auto z=pow(point_cloud_model.at(i)[2],2)*(-1);
        x2_sq.add_val(to_string(i+1), x);
        y2_sq.add_val(to_string(i+1), y);
        z2_sq.add_val(to_string(i+1), z);
        dm_root.add_val(to_string(i+1), (x+y+z)*(-1));
    }
    
//    indices idsij = indices("idsij");
//    idsij.add_empty_row();
//    param<> dm("dm");
//    dm.in(cells);
//    param<> dm_pr("dm_pr");
//    dm_pr.in(cells);
//    param<> dd("dd");
//    for(auto i=0;i<nd;i++){
//        auto di=pow(point_cloud_data.at(i)[0],2)+pow(point_cloud_data.at(i)[1],2)+pow(point_cloud_data.at(i)[2],2);
//        auto di_r=(sqrt(di)+sqrt(shift_mag_max));
//        dd.add_val(to_string(i+1), pow(di_r,2));
//        for(auto j=1;j<=nm;j++){
//            if(cells.has_key(to_string(i+1)+","+to_string(j))){
//                auto dmd=pow(point_cloud_model.at(j-1)[0],2)+pow(point_cloud_model.at(j-1)[1],2)+pow(point_cloud_model.at(j-1)[2],2);
//                auto dmd_r=sqrt(dmd)*di_r*(-1);
//                dm.set_val(to_string(i+1)+","+to_string(j), dmd);
//                dm_pr.set_val(to_string(i+1)+","+to_string(j), dmd_r);
//                idsij.add_in_row(i, to_string(i+1)+","+to_string(j));
//            }
//        }
//    }
    
//    Constraint<> txsq("txsq");
//    txsq = pow(x_shift,2) - tx;
//    //txsq.add_to_callback();
//    Reg->add(txsq.in(range(0,0))<=0);
//
//    Constraint<> tysq("tysq");
//    tysq = pow(y_shift,2) - ty;
//    //tysq.add_to_callback();
//    Reg->add(tysq.in(range(0,0))<=0);
//
//    Constraint<> tzsq("tzsq");
//    tzsq = pow(z_shift,2) - tz;
//    // tzsq.add_to_callback();
//    Reg->add(tzsq.in(range(0,0))<=0);
    
    indices ids = indices("in_x");
    ids.add_empty_row();
    for(auto i=0;i<nd;i++){
        for(auto j=1;j<=nm;j++){
            if(cells.has_key(to_string(i+1)+","+to_string(j)))
                ids.add_in_row(i, to_string(j));
        }
    }
    
//    auto ids1 = theta11.repeat_id(N1.size());
//
//    auto idst = theta11.repeat_id(cells.size());
//    Constraint<> def_c("def_c");
//    def_c=c-(x1.from(cells)*theta11.in(idst) + y1.from(cells)*theta12.in(idst) + z1.from(cells)*theta13.in(idst))*x2.to(cells)-(x1.from(cells)*theta21.in(idst) + y1.from(cells)*theta22.in(idst) + z1.from(cells)*theta23.in(idst))*y2.to(cells)-(x1.from(cells)*theta31.in(idst) + y1.from(cells)*theta32.in(idst) + z1.from(cells)*theta33.in(idst))*z2.to(cells);
//    //Reg->add(def_c.in(cells)<=0);
//
//    Constraint<> def_c1("def_c1");
//    def_c1=c-(x1.from(cells)*theta11.in(idst) + y1.from(cells)*theta12.in(idst) + z1.from(cells)*theta13.in(idst))*x2.to(cells)-xxU.to(cells)-(x1.from(cells)*theta21.in(idst) + y1.from(cells)*theta22.in(idst) + z1.from(cells)*theta23.in(idst))*y2.to(cells)-yyU.to(cells)-(x1.from(cells)*theta31.in(idst) + y1.from(cells)*theta32.in(idst) + z1.from(cells)*theta33.in(idst))*z2.to(cells)-zzU.to(cells);
    //Reg->add(def_c1.in(cells)>=0);
    
//        var<> xij("xij", xL, xU), yij("yij", yL, yU), zij("zij", zL, zU);
//      //  Reg->add(xij.in(cells), yij.in(cells), zij.in(cells));
//    //
//        Constraint<> cl("cl");
//        cl=c-xij*x2.to(cells)-yij*y2.to(cells)-zij*z2.to(cells);
//        //Reg->add(cl.in(cells)<=0);
//
//
//        Constraint<> region("region");
//        region=la.in(vface)*xij.in(vnorm)+lb.in(vface)*yij.in(vnorm)+lc.in(vface)*zij.in(vnorm)+ld.in(facets);
//       // Reg->add(region.in(facets)<=0);
//
//        Constraint<> dist_data2("dist_data2");
//        dist_data2=xL*xU+yL*yU+zL*zU+pow((sqrt(pow(x1.from(cells),2)+pow(y1.from(cells),2)+pow(z1.from(cells),2))-sqrt(shift_mag_max)),2)-xij*(xL+xU)-yij*(yL+yU)-zij*(zL+zU);
       // Reg->add(dist_data2.in(cells)<=0);
    
    var<> new_xm("new_xm", -1, 1), new_ym("new_ym", -1, 1),new_zm("new_zm", -1, 1);
    Reg->add(new_xm.in(N1), new_ym.in(N1),new_zm.in(N1));
    Constraint<> Def_newxm("Def_newxm");
    Def_newxm = new_xm-product(x2.in(ids),bin.in_matrix(1, 1));
    Reg->add(Def_newxm.in(N1)==0);
    
    Constraint<> Def_newym("Def_newym");
    Def_newym = new_ym-product(y2.in(ids),bin.in_matrix(1, 1));
    Reg->add(Def_newym.in(N1)==0);
    
    Constraint<> Def_newzm("Def_newzm");
    Def_newzm = new_zm-product(z2.in(ids),bin.in_matrix(1, 1));
    Reg->add(Def_newzm.in(N1)==0);
    
    Constraint<> Def_tx("Def_tx");
    Def_tx = nd*x_shift-sum(new_xm);
    Reg->add(Def_tx==0);
    
    Constraint<> Def_ty("Def_ty");
    Def_ty = nd*y_shift-sum(new_ym);
    Reg->add(Def_ty==0);
    
    Constraint<> Def_tz("Def_tz");
    Def_tz = nd*z_shift-sum(new_zm);
    Reg->add(Def_tz==0);
    

    
    param<> two("2");
    two.in(cells);
    two = 2;
    
    //    Constraint<> def_delta("def_delta");
    //    def_delta=(sum(new_x1*new_x1)+sum(new_y1*new_y1)+sum(new_z1*new_z1)-two.tr()*cx-two.tr()*cy-two.tr()*cz+sum(new_xm*new_xm)+sum(new_ym*new_ym)+sum(new_zm*new_zm))-delta;
    //    Reg->add(def_delta<=0);
    /* Objective function */
    
    func<> obj = 0;
    
    
    //obj+=(sum(x1*x1)+sum(y1*y1)+sum(z1*z1)+nd*(tx+ty+tz)-two.tr()*(c*bin)+sum(dm*bin));
   // obj+=(sum(x1*x1)+sum(y1*y1)+sum(z1*z1)-nd*(x_shift*x_shift+y_shift*y_shift+z_shift*z_shift)-two.tr()*(c*bin)+sum(dm*bin));
    
   // obj+=(sum(x1*x1)+sum(y1*y1)+sum(z1*z1)-nd*(x_shift*(new_shift_min_x+new_shift_max_x)+y_shift*(new_shift_min_y+new_shift_max_y)+z_shift*(new_shift_min_z+new_shift_max_z)-new_shift_min_x*new_shift_max_x-new_shift_min_y*new_shift_max_y-new_shift_min_z*new_shift_max_z)-product(cU, bin)+sum(dm*bin));
    
    obj+=(sum(x1*x1)+sum(y1*y1)+sum(z1*z1)-nd*(x_shift*(new_shift_min_x+new_shift_max_x)+y_shift*(new_shift_min_y+new_shift_max_y)+z_shift*(new_shift_min_z+new_shift_max_z)-new_shift_min_x*new_shift_max_x-new_shift_min_y*new_shift_max_y-new_shift_min_z*new_shift_max_z)+sum(cU*bin));
    
    
    //obj+=(sum(new_x1*new_x1)+sum(new_y1*new_y1)+sum(new_z1*new_z1)-two.tr()*cx-two.tr()*cy-two.tr()*cz+sum(new_xm*new_xm)+sum(new_ym*new_ym)+sum(new_zm*new_zm));
    //obj+=delta;
    Reg->min(obj);
    
    //Reg->print();
    //solver<> S(Reg,gurobi);
    //S.use_callback();
    //S.run();
    //Reg->print_solution();
//    vector<double> solution_lb, rot_trans_ub;
//    //solution_lb.resize(Reg->_nb_vars);
//
//    rot_trans_ub.resize(12);
//    // Reg->get_solution(solution_lb);
//
//
//    double det=theta11.eval()*(theta22.eval()*theta33.eval()-theta32.eval()*theta23.eval())
//    -theta12.eval()*(theta21.eval()*theta33.eval()-theta31.eval()*theta23.eval())+theta13.eval()*(theta21.eval()*theta32.eval()-theta31.eval()*theta22.eval());
//
//
//    DebugOn("row 12 " << (theta11.eval()*theta21.eval())+(theta12.eval()*theta22.eval())+(theta13.eval()*theta23.eval())
//            << endl);
//    DebugOn("row 13 " << (theta11.eval()*theta31.eval())+(theta12.eval()*theta32.eval())+(theta13.eval()*theta33.eval())
//            << endl);
//    DebugOn("row 23 " << (theta21.eval()*theta31.eval())+(theta22.eval()*theta32.eval())+(theta23.eval()*theta33.eval())
//            << endl);
//
//    DebugOn("Determinant "<<det<<endl);
    return(Reg);
}
shared_ptr<Model<double>> build_polyhedral_ipopt(vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, const indices& valid_cells, double new_roll_min, double new_roll_max, double new_pitch_min, double new_pitch_max, double new_yaw_min, double new_yaw_max, double new_shift_min_x, double new_shift_max_x, double new_shift_min_y, double new_shift_max_y, double new_shift_min_z, double new_shift_max_z, vector<double>& rot_trans, bool separate, const vector<pair<pair<int,int>,pair<int,int>>>& incompatibles,  param<>& norm_x,  param<>& norm_y,  param<>& norm_z,  param<>& intercept,const vector<int>& init_matching, const vector<double>& error_per_point, const vector<vector<vector<double>>>& model_voronoi_normals, const vector<vector<double>>& model_face_intercept, bool relax_inits){
    
    int nb_pairs = 0, min_nb_pairs = numeric_limits<int>::max(), max_nb_pairs = 0, av_nb_pairs = 0;
    size_t nm = point_cloud_model.size(), nd = point_cloud_data.size();
    vector<double> zeros = {0,0,0};
    param<> x1("x1"), x2("x2"), y1("y1"), y2("y2"), z1("z1"), z2("z2");
    int m = av_nb_pairs;
    string i_str, j_str;
    double xm_max = numeric_limits<double>::lowest(), ym_max = numeric_limits<double>::lowest(), zm_max = numeric_limits<double>::lowest();
    double xm_min = numeric_limits<double>::max(), ym_min = numeric_limits<double>::max(), zm_min = numeric_limits<double>::max();
    
    
    for (auto i = 0; i<nd; i++) {
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
    Debug("nd = " << nd << endl);
    Debug("nm = " << nm << endl);
    
    N1 = range(1,nd);
    N2 = range(1,nm);
    // auto cells = indices(N1,N2);
    auto cells=valid_cells;
    string name="lin_lb";
    
    auto Reg=make_shared<Model<>>(name);
    
    var<int> bin("bin",0,1);
    Reg->add(bin.in(cells));
    
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
    var<> x_shift("x_shift", new_shift_min_x, new_shift_max_x), y_shift("y_shift", new_shift_min_y, new_shift_max_y), z_shift("z_shift",new_shift_min_z,new_shift_max_z);
    Reg->add(x_shift.in(R(1)),y_shift.in(R(1)),z_shift.in(R(1)));
    
    
    DebugOn("Added " << cells.size() << " binary variables" << endl);
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
    
    DebugOn("Added " << cells.size() << " binary variables" << endl);
    
    
    
    Constraint<> diag_1("diag_1");
    diag_1=theta11+theta22-theta33-1;
    Reg->add(diag_1.in(range(0,0))<=0);
    Constraint<> diag_2("diag_2");
    diag_2=theta22+theta33-(1+theta11);
    Reg->add(diag_2.in(range(0,0))<=0);
    Constraint<> diag_3("diag_3");
    diag_3-=1+theta11+theta22+theta33;
    Reg->add(diag_3.in(range(0,0))<=0);
    Constraint<> diag_4("diag_4");
    diag_4=theta11-theta22+theta33-1;
    Reg->add(diag_4.in(range(0,0))<=0);

            Constraint<> soc_12("soc_12");
            soc_12 = pow(theta13+theta31,2)-(1-theta11-theta22+theta33)*(1+theta11-theta22-theta33);
           // soc_12.add_to_callback();
            Reg->add(soc_12.in(range(0,0))<=0);

            Constraint<> soc_13("soc_13");
            soc_13 = pow(theta12-theta21,2)-(1-theta11-theta22+theta33)*(1+theta11+theta22+theta33);
            //soc_13.add_to_callback();
            Reg->add(soc_13.in(range(0,0))<=0);

            Constraint<> soc_14("soc_14");
            soc_14 = pow(theta23+theta32,2)-(1-theta11-theta22+theta33)*(1-theta11+theta22-theta33);
            //soc_14.add_to_callback();
            Reg->add(soc_14.in(range(0,0))<=0);

            Constraint<> soc_23("soc_23");
            soc_23 = pow(theta23-theta32,2)-(1+theta11-theta22-theta33)*(1+theta11+theta22+theta33);
            //soc_23.add_to_callback();
            Reg->add(soc_23.in(range(0,0))<=0);

            Constraint<> soc_24("soc_24");
            soc_24 = pow(theta12+theta21,2)-(1+theta11-theta22-theta33)*(1-theta11+theta22-theta33);
            //soc_24.add_to_callback();
            Reg->add(soc_24.in(range(0,0))<=0);

            Constraint<> soc_34("soc_34");
            soc_34 = pow(theta31-theta13,2)-(1+theta11+theta22+theta33)*(1-theta11+theta22-theta33);
           // soc_34.add_to_callback();
            Reg->add(soc_34.in(range(0,0))<=0);
    //
            Constraint<> det_123("det_123");
            det_123+=(theta13+theta31)*((theta13+theta31)*(1+theta11+theta22+theta33)-(theta23-theta32)*(theta12-theta21));
            det_123-=(1-theta11-theta22+theta33)*((1+theta11-theta22-theta33)*(1+theta11+theta22+theta33)-pow(theta23-theta32,2));
            det_123-=(theta12-theta21)*((theta13+theta31)*(theta23-theta32)-(theta12-theta21)*(1+theta11-theta22-theta33));
            //det_123.add_to_callback();
            Reg->add(det_123.in(range(0,0))<=0);

            Constraint<> det_124("det_124");
            det_124+=(theta13+theta31)*((theta13+theta31)*(1-theta11+theta22-theta33)-(theta23+theta32)*(theta12+theta21));
            det_124-=(1-theta11-theta22+theta33)*((1+theta11-theta22-theta33)*(1-theta11+theta22-theta33)-pow(theta12+theta21,2));
            det_124-=(theta23+theta32)*((theta13+theta31)*(theta12+theta21)-(theta23+theta32)*(1+theta11-theta22-theta33));
            //det_124.add_to_callback();
            Reg->add(det_124.in(range(0,0))<=0);

            Constraint<> det_134("det_134");
            det_134+=(theta12-theta21)*((theta12-theta21)*(1-theta11+theta22-theta33)-(theta23+theta32)*(theta31-theta13));
            det_134-=(1-theta11-theta22+theta33)*((1+theta11+theta22+theta33)*(1-theta11+theta22-theta33)-pow(theta31-theta13,2));
            det_134-=(theta23+theta32)*((theta12-theta21)*(theta31-theta13)-(theta23+theta32)*(1+theta11+theta22+theta33));
            //det_134.add_to_callback();
            Reg->add(det_134.in(range(0,0))<=0);

            Constraint<> det_234("det_234");
            det_234+=(theta23-theta32)*((theta23-theta32)*(1-theta11+theta22-theta33)-(theta12+theta21)*(theta31-theta13));
            det_234-=(1+theta11-theta22-theta33)*((1+theta11+theta22+theta33)*(1-theta11+theta22-theta33)-pow(theta31-theta13,2));
            det_234-=(theta12+theta21)*((theta23-theta32)*(theta31-theta13)-(theta12+theta21)*(1+theta11+theta22+theta33));
            //det_234.add_to_callback();
            Reg->add(det_234.in(range(0,0))<=0);
    //
    Constraint<> row1("row1");
    row1 = pow(theta11,2)+pow(theta12,2)+pow(theta13,2);
    // row1.add_to_callback();
    Reg->add(row1.in(range(0,0))<=1);
    Constraint<> row2("row2");
    row2 = pow(theta21,2)+pow(theta22,2)+pow(theta23,2);
    //row2.add_to_callback();
    Reg->add(row2.in(range(0,0))<=1);
    Constraint<> row3("row3");
    row3 = pow(theta31,2)+pow(theta32,2)+pow(theta33,2);
    //row3.add_to_callback();
    Reg->add(row3.in(range(0,0))<=1);
    Constraint<> col1("col1");
    col1 = pow(theta11,2)+pow(theta21,2)+pow(theta31,2);
    // col1.add_to_callback();
    Reg->add(col1.in(range(0,0))<=1);
    Constraint<> col2("col2");
    col2 = pow(theta12,2)+pow(theta22,2)+pow(theta32,2);
    //col2.add_to_callback();
    Reg->add(col2.in(range(0,0))<=1);
    Constraint<> col3("col3");
    col3 = pow(theta13,2)+pow(theta23,2)+pow(theta33,2);
    //col3.add_to_callback();
    Reg->add(col3.in(range(0,0))<=1);
    
    
    Constraint<> OneBin("OneBin");
    OneBin = bin.in_matrix(1, 1);
    Reg->add(OneBin.in(N1)==1);

    param<> xL("xL");
    xL.in(cells);
    param<> xU("xU");
    xU.in(cells);
    param<> yL("yL");
    yL.in(cells);
    param<> yU("yU");
    yU.in(cells);
    param<> zL("zL");
    zL.in(cells);
    param<> zU("zU");
    zU.in(cells);
    param<> cL("cL");
    cL.in(cells);
    param<> cU("cU");
    cU.in(cells);
    param<> xxU("xxU"),xxL("xxL"),yyU("yyU"),yyL("yyL"),zzU("zzU"),zzL("zzL");
    
    for(auto j=1;j<=nm;j++){
        auto x=point_cloud_model.at(j-1)[0];
        auto y=point_cloud_model.at(j-1)[1];
        auto z=point_cloud_model.at(j-1)[2];
        xxL.add_val(to_string(j),std::min(new_shift_min_x*x, new_shift_max_x*x));
        xxU.add_val(to_string(j),std::max(new_shift_min_x*x, new_shift_max_x*x));
        yyL.add_val(to_string(j),std::min(new_shift_min_y*y, new_shift_max_y*y));
        yyU.add_val(to_string(j),std::max(new_shift_min_y*y, new_shift_max_y*y));
        zzL.add_val(to_string(j),std::min(new_shift_min_z*z, new_shift_max_z*z));
        zzU.add_val(to_string(j),std::max(new_shift_min_z*z, new_shift_max_z*z));
    }
    
    vector<vector<vector<double>>> box;
    vector<vector<double>> box_i;
    vector<double> coord_i;
    coord_i.resize(3);
    
    shared_ptr<pair<double,double>> x1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> y1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> z1_bounds = make_shared<pair<double,double>>();
    vector<double> x_lb, x_ub, y_lb, y_ub, z_lb, z_ub;
    indices vnorm("vnorm");
    vnorm.add_empty_row();
    indices vface("vface");
    vface.add_empty_row();
    int q=0;
    int nfacets=0;
    for(auto i=0;i<nd;i++){
        for(auto j=1;j<=nm;j++){
            if(cells.has_key(to_string(i+1)+","+to_string(j))){
                nfacets+=model_face_intercept.at(j-1).size();
            }
        }
    }
    indices facets("facets");
    facets=range(1, nfacets);
    param<> la("la"), lb("lb"),lc("lc"),ld("ld");
    la.in(facets);lb.in(facets);lc.in(facets);ld.in(facets);
    for(auto i=0;i<nd;i++){
        x1_bounds->first = point_cloud_data.at(i)[0];
        x1_bounds->second = point_cloud_data.at(i)[0];
        y1_bounds->first = point_cloud_data.at(i)[1];
        y1_bounds->second = point_cloud_data.at(i)[1];
        z1_bounds->first = point_cloud_data.at(i)[2];
        z1_bounds->second = point_cloud_data.at(i)[2];
        auto x_range  = get_product_range(x1_bounds, theta11._range);
        auto y_range  = get_product_range(y1_bounds, theta12._range);
        auto z_range  = get_product_range(z1_bounds, theta13._range);
        x_lb.push_back(x_range->first + y_range->first + z_range->first);
        x_ub.push_back(x_range->second + y_range->second + z_range->second);
        x_range  = get_product_range(x1_bounds, theta21._range);
        y_range  = get_product_range(y1_bounds, theta22._range);
        z_range  = get_product_range(z1_bounds, theta23._range);
        y_lb.push_back(x_range->first + y_range->first + z_range->first);
        y_ub.push_back(x_range->second + y_range->second + z_range->second);
        x_range  = get_product_range(x1_bounds, theta31._range);
        y_range  = get_product_range(y1_bounds, theta32._range);
        z_range  = get_product_range(z1_bounds, theta33._range);
        z_lb.push_back(x_range->first + y_range->first + z_range->first);
        z_ub.push_back(x_range->second + y_range->second + z_range->second);
        coord_i[0]=x_lb[i];
        coord_i[1]=y_lb[i];
        coord_i[2]=z_lb[i];
        box_i.push_back(coord_i);
        coord_i[0]=x_lb[i];
        coord_i[1]=y_ub[i];
        coord_i[2]=z_lb[i];
        box_i.push_back(coord_i);
        coord_i[0]=x_ub[i];
        coord_i[1]=y_ub[i];
        coord_i[2]=z_lb[i];
        box_i.push_back(coord_i);
        coord_i[0]=x_ub[i];
        coord_i[1]=y_lb[i];
        coord_i[2]=z_lb[i];
        box_i.push_back(coord_i);
        coord_i[0]=x_lb[i];
        coord_i[1]=y_lb[i];
        coord_i[2]=z_ub[i];
        box_i.push_back(coord_i);
        coord_i[0]=x_lb[i];
        coord_i[1]=y_ub[i];
        coord_i[2]=z_ub[i];
        box_i.push_back(coord_i);
        coord_i[0]=x_ub[i];
        coord_i[1]=y_ub[i];
        coord_i[2]=z_ub[i];
        box_i.push_back(coord_i);
        coord_i[0]=x_ub[i];
        coord_i[1]=y_lb[i];
        coord_i[2]=z_ub[i];
        box_i.push_back(coord_i);
        auto di=pow(point_cloud_data.at(i)[0],2)+pow(point_cloud_data.at(i)[1],2)+pow(point_cloud_data.at(i)[2],2);
        auto di_r=(sqrt(di));
        for(auto j=1;j<=nm;j++){
            if(cells.has_key(to_string(i+1)+","+to_string(j))){
                auto x=point_cloud_model.at(j-1)[0];
                auto y=point_cloud_model.at(j-1)[1];
                auto z=point_cloud_model.at(j-1)[2];
                double c_min=100, c_max=-999;
                for(auto k=0;k<8;k++){
                    auto inner_prod=box_i[k][0]*x+box_i[k][1]*y+box_i[k][2]*z;
                    if(inner_prod>=c_max){
                        c_max=inner_prod;
                    }
                    if(inner_prod<=c_min){
                        c_min=inner_prod;
                    }
                }
                auto normals=model_voronoi_normals.at(j-1);
                auto intercepts=model_face_intercept.at(j-1);
                for(auto k=0;k<intercepts.size();k++){
                    la.set_val(q, normals[k][0]);
                    lb.set_val(q, normals[k][1]);
                    lc.set_val(q, normals[k][2]);
                    ld.set_val(q, intercepts[k]);
                    vnorm.add_in_row(q, to_string(i+1)+","+to_string(j));
                    vface.add_in_row(q, to_string(q+1));
                    q++;
                }
                auto dmd=pow(point_cloud_model.at(j-1)[0],2)+pow(point_cloud_model.at(j-1)[1],2)+pow(point_cloud_model.at(j-1)[2],2);
                auto cL_val=std::max(sqrt(pow(x,2)+pow(y,2)+pow(z,2))*di_r*(-1), c_min);
                auto cU_val=std::min(sqrt(pow(x,2)+pow(y,2)+pow(z,2))*di_r, c_max);
                cL.set_val(to_string(i+1)+","+to_string(j), std::min(cL_val,0.0));
                cU.set_val(to_string(i+1)+","+to_string(j), std::max(cU_val,0.0));
                xL.set_val(to_string(i+1)+","+to_string(j), x_lb[i]);
                xU.set_val(to_string(i+1)+","+to_string(j), x_ub[i]);
                yL.set_val(to_string(i+1)+","+to_string(j), y_lb[i]);
                yU.set_val(to_string(i+1)+","+to_string(j), y_ub[i]);
                zL.set_val(to_string(i+1)+","+to_string(j), z_lb[i]);
                zU.set_val(to_string(i+1)+","+to_string(j), z_ub[i]);
            }
        }
        box_i.clear();
    }
    var<> c("c",cL, cU);
    Reg->add(c.in(cells));
    
    double tx_max=std::max(pow(new_shift_min_x,2), pow(new_shift_max_x,2));
    double ty_max=std::max(pow(new_shift_min_y,2), pow(new_shift_max_y,2));
    double tz_max=std::max(pow(new_shift_min_z,2), pow(new_shift_max_z,2));
    double tx_min, ty_min, tz_min;
    if(new_shift_min_x<=0 && new_shift_max_x>=0){
        tx_min=0;
    }
    else{
        tx_min=std::min(pow(new_shift_min_x,2), pow(new_shift_max_x,2));
    }
    if(new_shift_min_y<=0 && new_shift_max_y>=0){
        ty_min=0;
    }
    else{
        ty_min=std::min(pow(new_shift_min_y,2), pow(new_shift_max_y,2));
    }
    if(new_shift_min_z<=0 && new_shift_max_z>=0){
        tz_min=0;
    }
    else{
        tz_min=std::min(pow(new_shift_min_z,2), pow(new_shift_max_z,2));
    }
    var<> tx("tx", tx_min, tx_max), ty("ty",ty_min, ty_max), tz("tz", tz_min, tz_max);
    
    Reg->add(tx.in(R(1)),ty.in(R(1)),tz.in(R(1)));
    param<> dm_root("dm_root");
    param<> x2_sq("x2_sq");
    param<> y2_sq("y2_sq");
    param<> z2_sq("z2_sq");
    for(auto i=0;i<nm;i++){
        auto x=pow(point_cloud_model.at(i)[0],2)*(-1);
        auto y=pow(point_cloud_model.at(i)[1],2)*(-1);
        auto z=pow(point_cloud_model.at(i)[2],2)*(-1);
        x2_sq.add_val(to_string(i+1), x);
        y2_sq.add_val(to_string(i+1), y);
        z2_sq.add_val(to_string(i+1), z);
        dm_root.add_val(to_string(i+1), (x+y+z)*(-1));
    }
    
    indices idsij = indices("idsij");
    idsij.add_empty_row();
    param<> dm("dm");
    dm.in(cells);
    param<> dm_pr("dm_pr");
    dm_pr.in(cells);
    param<> dd("dd");
    for(auto i=0;i<nd;i++){
        auto di=pow(point_cloud_data.at(i)[0],2)+pow(point_cloud_data.at(i)[1],2)+pow(point_cloud_data.at(i)[2],2);
        auto di_r=(sqrt(di)+sqrt(shift_mag_max));
        dd.add_val(to_string(i+1), pow(di_r,2));
        for(auto j=1;j<=nm;j++){
            if(cells.has_key(to_string(i+1)+","+to_string(j))){
                auto dmd=pow(point_cloud_model.at(j-1)[0],2)+pow(point_cloud_model.at(j-1)[1],2)+pow(point_cloud_model.at(j-1)[2],2);
                auto dmd_r=sqrt(dmd)*di_r*(-1);
                dm.set_val(to_string(i+1)+","+to_string(j), dmd);
                dm_pr.set_val(to_string(i+1)+","+to_string(j), dmd_r);
                idsij.add_in_row(i, to_string(i+1)+","+to_string(j));
            }
        }
    }
    
    Constraint<> txsq("txsq");
    txsq = pow(x_shift,2) - tx;
    //txsq.add_to_callback();
    Reg->add(txsq.in(range(0,0))<=0);

    Constraint<> tysq("tysq");
    tysq = pow(y_shift,2) - ty;
    //tysq.add_to_callback();
    Reg->add(tysq.in(range(0,0))<=0);

    Constraint<> tzsq("tzsq");
    tzsq = pow(z_shift,2) - tz;
    // tzsq.add_to_callback();
    Reg->add(tzsq.in(range(0,0))<=0);
    
    indices ids = indices("in_x");
    ids.add_empty_row();
    for(auto i=0;i<nd;i++){
        for(auto j=1;j<=nm;j++){
            if(cells.has_key(to_string(i+1)+","+to_string(j)))
                ids.add_in_row(i, to_string(j));
        }
    }
    
//    auto ids1 = theta11.repeat_id(N1.size());
//
    auto idst = theta11.repeat_id(cells.size());
    Constraint<> def_c("def_c");
    def_c=c+cL*(1-bin)-((x1.from(cells)*theta11.in(idst) + y1.from(cells)*theta12.in(idst) + z1.from(cells)*theta13.in(idst))*x2.to(cells)+(x1.from(cells)*theta21.in(idst) + y1.from(cells)*theta22.in(idst) + z1.from(cells)*theta23.in(idst))*y2.to(cells)+(x1.from(cells)*theta31.in(idst) + y1.from(cells)*theta32.in(idst) + z1.from(cells)*theta33.in(idst))*z2.to(cells));
    Reg->add(def_c.in(cells)<=0);
    
    Constraint<> def_c1("def_c1");
    def_c1=c-cU*bin;
    Reg->add(def_c1.in(cells)<=0);
    
    
//    auto idst = theta11.repeat_id(cells.size());
//    Constraint<> def_c("def_c");
//    def_c=c+cL*(1-bin)-((x1.from(cells)*theta11.in(idst) + y1.from(cells)*theta12.in(idst) + z1.from(cells)*theta13.in(idst)+x_shift)*x2.to(cells)+(x1.from(cells)*theta21.in(idst) + y1.from(cells)*theta22.in(idst) + z1.from(cells)*theta23.in(idst)+y_shift)*y2.to(cells)+(x1.from(cells)*theta31.in(idst) + y1.from(cells)*theta32.in(idst) + z1.from(cells)*theta33.in(idst)+z_shift)*z2.to(cells));
//    Reg->add(def_c.in(cells)<=0);
//
//    Constraint<> def_c1("def_c1");
//    def_c1=c-cU*bin;
//    Reg->add(def_c1.in(cells)<=0);
//
//    Constraint<> def_c1("def_c1");
//    def_c1=c-(x1.from(cells)*theta11.in(idst) + y1.from(cells)*theta12.in(idst) + z1.from(cells)*theta13.in(idst))*x2.to(cells)-xxU.to(cells)-(x1.from(cells)*theta21.in(idst) + y1.from(cells)*theta22.in(idst) + z1.from(cells)*theta23.in(idst))*y2.to(cells)-yyU.to(cells)-(x1.from(cells)*theta31.in(idst) + y1.from(cells)*theta32.in(idst) + z1.from(cells)*theta33.in(idst))*z2.to(cells)-zzU.to(cells);
    //Reg->add(def_c1.in(cells)>=0);
    
//        var<> xij("xij", xL, xU), yij("yij", yL, yU), zij("zij", zL, zU);
//      //  Reg->add(xij.in(cells), yij.in(cells), zij.in(cells));
//    //
//        Constraint<> cl("cl");
//        cl=c-xij*x2.to(cells)-yij*y2.to(cells)-zij*z2.to(cells);
//        //Reg->add(cl.in(cells)<=0);
//
//
//        Constraint<> region("region");
//        region=la.in(vface)*xij.in(vnorm)+lb.in(vface)*yij.in(vnorm)+lc.in(vface)*zij.in(vnorm)+ld.in(facets);
//       // Reg->add(region.in(facets)<=0);
//
//        Constraint<> dist_data2("dist_data2");
//        dist_data2=xL*xU+yL*yU+zL*zU+pow((sqrt(pow(x1.from(cells),2)+pow(y1.from(cells),2)+pow(z1.from(cells),2))-sqrt(shift_mag_max)),2)-xij*(xL+xU)-yij*(yL+yU)-zij*(zL+zU);
       // Reg->add(dist_data2.in(cells)<=0);
    
    var<> new_xm("new_xm", -1, 1), new_ym("new_ym", -1, 1),new_zm("new_zm", -1, 1);
    Reg->add(new_xm.in(N1), new_ym.in(N1),new_zm.in(N1));
    Constraint<> Def_newxm("Def_newxm");
    Def_newxm = new_xm-product(x2.in(ids),bin.in_matrix(1, 1));
    Reg->add(Def_newxm.in(N1)==0);
    
    Constraint<> Def_newym("Def_newym");
    Def_newym = new_ym-product(y2.in(ids),bin.in_matrix(1, 1));
    Reg->add(Def_newym.in(N1)==0);
    
    Constraint<> Def_newzm("Def_newzm");
    Def_newzm = new_zm-product(z2.in(ids),bin.in_matrix(1, 1));
    Reg->add(Def_newzm.in(N1)==0);
    
    Constraint<> Def_tx("Def_tx");
    Def_tx = nd*x_shift-sum(new_xm);
    Reg->add(Def_tx==0);
    
    Constraint<> Def_ty("Def_ty");
    Def_ty = nd*y_shift-sum(new_ym);
    Reg->add(Def_ty==0);
    
    Constraint<> Def_tz("Def_tz");
    Def_tz = nd*z_shift-sum(new_zm);
    Reg->add(Def_tz==0);
    

    
    param<> two("2");
    two.in(cells);
    two = 2;
    
    //    Constraint<> def_delta("def_delta");
    //    def_delta=(sum(new_x1*new_x1)+sum(new_y1*new_y1)+sum(new_z1*new_z1)-two.tr()*cx-two.tr()*cy-two.tr()*cz+sum(new_xm*new_xm)+sum(new_ym*new_ym)+sum(new_zm*new_zm))-delta;
    //    Reg->add(def_delta<=0);
    /* Objective function */
    
    func<> obj = 0;
    
    
    //obj+=(sum(x1*x1)+sum(y1*y1)+sum(z1*z1)+nd*(tx+ty+tz)-two.tr()*(c)+sum(dm*bin));
   // obj+=(sum(x1*x1)+sum(y1*y1)+sum(z1*z1)-nd*(x_shift*x_shift+y_shift*y_shift+z_shift*z_shift)-two.tr()*(c*bin)+sum(dm*bin));
    

    
    obj+=(sum(x1*x1)+sum(y1*y1)+sum(z1*z1)-nd*(x_shift*(new_shift_min_x+new_shift_max_x)+y_shift*(new_shift_min_y+new_shift_max_y)+z_shift*(new_shift_min_z+new_shift_max_z)-new_shift_min_x*new_shift_max_x-new_shift_min_y*new_shift_max_y-new_shift_min_z*new_shift_max_z)+sum(dm*bin)-sum(two*c));
    
    
    //obj+=(sum(new_x1*new_x1)+sum(new_y1*new_y1)+sum(new_z1*new_z1)-two.tr()*cx-two.tr()*cy-two.tr()*cz+sum(new_xm*new_xm)+sum(new_ym*new_ym)+sum(new_zm*new_zm));
    //obj+=delta;
    Reg->min(obj);
    
    //Reg->print();
    //solver<> S(Reg,gurobi);
    //S.use_callback();
    //S.run();
    //Reg->print_solution();
//    vector<double> solution_lb, rot_trans_ub;
//    //solution_lb.resize(Reg->_nb_vars);
//
//    rot_trans_ub.resize(12);
//    // Reg->get_solution(solution_lb);
//
//
//    double det=theta11.eval()*(theta22.eval()*theta33.eval()-theta32.eval()*theta23.eval())
//    -theta12.eval()*(theta21.eval()*theta33.eval()-theta31.eval()*theta23.eval())+theta13.eval()*(theta21.eval()*theta32.eval()-theta31.eval()*theta22.eval());
//
//
//    DebugOn("row 12 " << (theta11.eval()*theta21.eval())+(theta12.eval()*theta22.eval())+(theta13.eval()*theta23.eval())
//            << endl);
//    DebugOn("row 13 " << (theta11.eval()*theta31.eval())+(theta12.eval()*theta32.eval())+(theta13.eval()*theta33.eval())
//            << endl);
//    DebugOn("row 23 " << (theta21.eval()*theta31.eval())+(theta22.eval()*theta32.eval())+(theta23.eval()*theta33.eval())
//            << endl);
//
//    DebugOn("Determinant "<<det<<endl);
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
    DebugOn("Theta matrix = " << endl);
    DebugOn("|" << theta11.eval() << " " << theta12.eval() << " " << theta13.eval() << "|" << endl);
    DebugOn("|" << theta21.eval() << " " << theta22.eval() << " " << theta23.eval() << "|" << endl);
    DebugOn("|" << theta31.eval() << " " << theta32.eval() << " " << theta33.eval() << "|" << endl);
    constant<> row1 = pow(theta11.eval(),2)+pow(theta12.eval(),2)+pow(theta13.eval(),2);
    constant<> row2 = pow(theta21.eval(),2)+pow(theta22.eval(),2)+pow(theta23.eval(),2);
    constant<> row3 = pow(theta31.eval(),2)+pow(theta32.eval(),2)+pow(theta33.eval(),2);
    constant<> col1 = pow(theta11.eval(),2)+pow(theta21.eval(),2)+pow(theta31.eval(),2);
    constant<> col2 = pow(theta12.eval(),2)+pow(theta22.eval(),2)+pow(theta32.eval(),2);
    constant<> col3 = pow(theta13.eval(),2)+pow(theta23.eval(),2)+pow(theta33.eval(),2);
    DebugOn("row 1 " << row1.eval() << endl);
    DebugOn("row 2 " << row2.eval() << endl);
    DebugOn("row 3 " << row3.eval() << endl);
    DebugOn("col 1 " << col1.eval() << endl);
    DebugOn("col 2 " << col2.eval() << endl);
    DebugOn("col 3 " << col3.eval() << endl);
    constant<> det=theta11.eval()*(theta22.eval()*theta33.eval()-theta32.eval()*theta23.eval())
    -theta12.eval()*(theta21.eval()*theta33.eval()-theta31.eval()*theta23.eval())+theta13.eval()*(theta21.eval()*theta32.eval()-theta31.eval()*theta22.eval());
    constant<> row12 = (theta11.eval()*theta21.eval())+(theta12.eval()*theta22.eval())+(theta13.eval()*theta23.eval());
    constant<> row13 = (theta11.eval()*theta31.eval())+(theta12.eval()*theta32.eval())+(theta13.eval()*theta33.eval());
    constant<> row23 = (theta21.eval()*theta31.eval())+(theta22.eval()*theta32.eval())+(theta23.eval()*theta33.eval());
    DebugOn("row 12 " << row12.eval() << endl);
    DebugOn("row 13 " << row13.eval() << endl);
    DebugOn("row 23 " << row23.eval() << endl);
    
    DebugOn("Determinant "<<det.eval()<<endl);
    
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
    DebugOn("Roll (degrees) = " << to_string_with_precision(roll_val,12) << endl);
    DebugOn("Pitch (degrees) = " << to_string_with_precision(pitch_val,12) << endl);
    DebugOn("Yaw (degrees) = " << to_string_with_precision(yaw_val,12) << endl);
    DebugOn("x shift = " << x_shift.eval() << endl);
    DebugOn("y shift = " << y_shift.eval() << endl);
    DebugOn("z shift = " << z_shift.eval() << endl);
    if(!is_rotation){
        DebugOn("WARNING, returned matrix is not a Rotation!\n");
    }
    return is_rotation;
}


shared_ptr<Model<double>> three_point_model(vector<vector<double>> d, vector<vector<double>>m, bool convex){
    auto Reg=make_shared<Model<>>("bounds");
    var<> theta11("theta11",  -1, 1), theta12("theta12", -1, 1), theta13("theta13", -1, 1);
    var<> theta21("theta21",  -1, 1), theta22("theta22", -1, 1), theta23("theta23", -1, 1);
    var<> theta31("theta31",  -1, 1), theta32("theta32", -1, 1), theta33("theta33", -1, 1);
    double shift_min_x = -1, shift_max_x = 1, shift_min_y = -1,shift_max_y = 1,shift_min_z = -1,shift_max_z = 1;
    var<> x_shift("x_shift", shift_min_x, shift_max_x), y_shift("y_shift", shift_min_y, shift_max_y), z_shift("z_shift", shift_min_z, shift_max_z);
    var<> tx("tx", 0,1), ty("ty",0,1), tz("tz",0,1);
    Reg->add(theta11.in(R(1)),theta12.in(R(1)),theta13.in(R(1)));
    Reg->add(theta21.in(R(1)),theta22.in(R(1)),theta23.in(R(1)));
    Reg->add(theta31.in(R(1)),theta32.in(R(1)),theta33.in(R(1)));
    Reg->add(x_shift.in(R(1)),y_shift.in(R(1)),z_shift.in(R(1)));
    Reg->add(tx.in(R(1)),ty.in(R(1)),tz.in(R(1)));
    param<> x1("x1"), x2("x2"), y1("y1"), y2("y2"), z1("z1"), z2("z2"),dd,dm;
    for (auto i = 0; i<3; i++) {
        auto i_str = to_string(i+1);
        x1.add_val(i_str,d[i].at(0));
        y1.add_val(i_str,d[i].at(1));
        z1.add_val(i_str,d[i].at(2));
        dd.add_val(i_str, pow(d[i].at(0),2)+pow(d[i].at(1),2)+pow(d[i].at(2),2));
    }
    for (auto i = 0; i<3; i++) {
        auto i_str = to_string(i+1);
        x2.add_val(i_str,m[i].at(0));
        y2.add_val(i_str,m[i].at(1));
        z2.add_val(i_str,m[i].at(2));
        dm.add_val(i_str, pow(d[i].at(0),2)+pow(d[i].at(1),2)+pow(d[i].at(2),2));
    }
    x1.print();
    y1.print();
    z1.print();
    x2.print();
    y2.print();
    z2.print();
    dd.print();
    dm.print();
    // auto a=sum(dd);
    //auto b=sum(dm);
    // DebugOn(" "<< a.eval() <<" "<< b.eval() <<endl);
    indices N1("N1");
    N1 = range(1,3);
    auto ids_repeat = x_shift.repeat_id(N1.size());
    func<> obj =sum(x1.in(N1)*x1.in(N1) + y1.in(N1)*y1.in(N1) + z1.in(N1)*z1.in(N1));
    obj +=3*(tx +ty+tz);
    
    obj -= 2*sum(x2.in(N1)*x_shift.in(ids_repeat)) + 2*sum(y2.in(N1)*y_shift.in(ids_repeat)) + 2*sum(z2.in(N1)*z_shift.in(ids_repeat));
    
    obj += sum(x2.in(N1)*x2.in(N1)) + sum(y2.in(N1)*y2.in(N1)) + sum(z2.in(N1)*z2.in(N1));
    
    auto ids1 = theta11.repeat_id(N1.size());
    obj -= 2*sum(x2.in(N1)*x1.in(N1)*theta11.in(ids1));
    obj -= 2*sum(x2.in(N1)*y1.in(N1)*theta12.in(ids1));
    obj -= 2*sum(x2.in(N1)*z1.in(N1)*theta13.in(ids1));
    obj -= 2*sum(y2.in(N1)*x1.in(N1)*theta21.in(ids1));
    obj -= 2*sum(y2.in(N1)*y1.in(N1)*theta22.in(ids1));
    obj -= 2*sum(y2.in(N1)*z1.in(N1)*theta23.in(ids1));
    obj -= 2*sum(z2.in(N1)*x1.in(N1)*theta31.in(ids1));
    obj -= 2*sum(z2.in(N1)*y1.in(N1)*theta32.in(ids1));
    obj -= 2*sum(z2.in(N1)*z1.in(N1)*theta33.in(ids1));
    
    
    Constraint<> soc1("soc1");
    soc1 = pow((theta13+theta31),2)-(1-theta11-theta22+theta33)*(1+theta11-theta22-theta33);
    Reg->add(soc1.in(range(0,0))<=0);
    
    Constraint<> soc2("soc2");
    soc2 = pow((theta12-theta21),2)-(1-theta11-theta22+theta33)*(1+theta11+theta22+theta33);
    Reg->add(soc2.in(range(0,0))<=0);
    
    Constraint<> soc3("soc3");
    soc3 = pow((theta23+theta32),2)-(1-theta11-theta22+theta33)*(1-theta11+theta22-theta33);
    Reg->add(soc3.in(range(0,0))<=0);
    
    Constraint<> soc4("soc4");
    soc4 = pow((theta23-theta32),2)-(1+theta11-theta22-theta33)*(1+theta11+theta22+theta33);
    Reg->add(soc4.in(range(0,0))<=0);
    
    Constraint<> soc5("soc5");
    soc5 = pow((theta12+theta21),2)-(1+theta11-theta22-theta33)*(1-theta11+theta22-theta33);
    Reg->add(soc5.in(range(0,0))<=0);
    
    Constraint<> soc6("soc6");
    soc6 = pow((theta31-theta13),2)-(1+theta11+theta22+theta33)*(1-theta11+theta22-theta33);
    Reg->add(soc6.in(range(0,0))<=0);
    if(convex){
        Constraint<> row1("row1");
        row1 = pow(theta11,2)+pow(theta12,2)+pow(theta13,2);
        Reg->add(row1.in(range(0,0))<=1);
        Constraint<> row2("row2");
        row2 = pow(theta21,2)+pow(theta22,2)+pow(theta23,2);
        Reg->add(row2.in(range(0,0))<=1);
        Constraint<> row3("row3");
        row3 = pow(theta31,2)+pow(theta32,2)+pow(theta33,2);
        Reg->add(row3.in(range(0,0))<=1);
        Constraint<> col1("col1");
        col1 = pow(theta11,2)+pow(theta21,2)+pow(theta31,2);
        Reg->add(col1.in(range(0,0))<=1);
        Constraint<> col2("col2");
        col2 = pow(theta12,2)+pow(theta22,2)+pow(theta32,2);
        Reg->add(col2.in(range(0,0))<=1);
        Constraint<> col3("col3");
        col3 = pow(theta13,2)+pow(theta23,2)+pow(theta33,2);
        Reg->add(col3.in(range(0,0))<=1);
    }
    else{
        Constraint<> row1("row1");
        row1 = pow(theta11,2)+pow(theta12,2)+pow(theta13,2);
        Reg->add(row1.in(range(0,0))==1);
        Constraint<> row2("row2");
        row2 = pow(theta21,2)+pow(theta22,2)+pow(theta23,2);
        Reg->add(row2.in(range(0,0))==1);
        Constraint<> row3("row3");
        row3 = pow(theta31,2)+pow(theta32,2)+pow(theta33,2);
        Reg->add(row3.in(range(0,0))==1);
        Constraint<> col1("col1");
        col1 = pow(theta11,2)+pow(theta21,2)+pow(theta31,2);
        Reg->add(col1.in(range(0,0))==1);
        Constraint<> col2("col2");
        col2 = pow(theta12,2)+pow(theta22,2)+pow(theta32,2);
        Reg->add(col2.in(range(0,0))==1);
        Constraint<> col3("col3");
        col3 = pow(theta13,2)+pow(theta23,2)+pow(theta33,2);
        Reg->add(col3.in(range(0,0))==1);
    }
    
    //obj.print_symbolic();
    Reg->min(obj);
    Reg->print();
    return Reg;
}
/* Run the ARMO model for registration */
tuple<double,double,double,double,double,double> run_ARMO(bool bypass, string axis, const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data){
    auto thetax = atan2(-0.0081669, -0.0084357)/2;
    auto thetay = atan2(0.9999311, std::sqrt(0.0081669*0.0081669+0.0084357*0.0084357))/2;
    auto thetaz = atan2(-0.0081462,-0.0084556)/2;
    DebugOff("thetax = " << thetax << endl);
    DebugOff("thetay = " << thetay << endl);
    DebugOff("thetaz = " << thetaz << endl);
    double shift_min_x = 0.125, shift_max_x = 0.25, shift_min_y = -0.25,shift_max_y = -0.125,shift_min_z = -0.125,shift_max_z = 0;
    double yaw_min = -12.5*pi/180., yaw_max = 0, pitch_min = 12.5*pi/180.,pitch_max = 25.*pi/180.,roll_min = -12.5*pi/180.,roll_max = 0;
    
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
        Debug("nd = " << nd << endl);
        Debug("nm = " << nm << endl);
        
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
            var<> yaw("yaw", yaw_min, yaw_max), pitch("pitch", pitch_min, pitch_max), roll("roll", roll_min, roll_max);
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

tuple<double,double,double> run_ARMO(string axis, const vector<vector<double>>& point_cloud1, const vector<vector<double>>& point_cloud2, const vector<vector<double>>& uav1, const vector<vector<double>>& uav2){};
//{
//    double angle_max = 0.05;
//    int nb_pairs = 0, min_nb_pairs = numeric_limits<int>::max(), max_nb_pairs = 0, av_nb_pairs = 0;
//    vector<pair<double,double>> min_max1;
//    vector<vector<pair<double,double>>> min_max2(point_cloud2.size());
//    vector<int> nb_neighbors(point_cloud1.size());
//    vector<vector<int>> neighbors;
//    /* Compute cube for all points in point cloud 2 */
//    for (auto i = 0; i<point_cloud2.size(); i++) {
//        min_max2[i] = get_min_max(angle_max, point_cloud2.at(i), uav2.at(i));
//    }
//    double roll_1 = 0, yaw_1 = 0, pitch_1 = 0;
//    bool bypass = false;
//    if(!bypass){
//        /* Check if cubes intersect */
//        neighbors.resize(point_cloud1.size());
//        for (auto i = 0; i<point_cloud1.size(); i++) {
//            nb_pairs = 0;
//            min_max1 = get_min_max(angle_max, point_cloud1.at(i), uav1.at(i));
//            DebugOff("For point (" << point_cloud1.at(i).at(0) << "," <<  point_cloud1.at(i).at(1) << "," << point_cloud1.at(i).at(2) << "): ");
//            DebugOff("\n neighbors in umbrella : \n");
//            for (size_t j = 0; j < point_cloud2.size(); j++){
//                if(intersect(min_max1, min_max2[j])){ /* point is in umbrella */
//                    nb_pairs++;
//                    neighbors[i].push_back(j);
//                    DebugOff("(" << point_cloud2.at(j).at(0) << "," <<  point_cloud2.at(j).at(1) << "," << point_cloud2.at(j).at(2) << ")\n");
//                }
//            }
//
//            DebugOff("nb points in umbrella = " << nb_pairs << endl);
//            if(nb_pairs>max_nb_pairs)
//                max_nb_pairs = nb_pairs;
//            if(nb_pairs<min_nb_pairs)
//                min_nb_pairs = nb_pairs;
//            av_nb_pairs += nb_pairs;
//
//                //        std::cout << "For point (" << point_cloud1.at(i).at(0) << "," <<  point_cloud1.at(i).at(1) << "," << point_cloud1.at(i).at(2) << ")"<< " knnSearch(n="<<m<<"): \n";
//                //        for (size_t k = 0; k < m; k++)
//                //            std::cout << "ret_index["<<k<<"]=" << ret_indexes[k] << " out_dist_sqr=" << out_dists_sqr[k] << " point = (" << point_cloud2.at(ret_indexes[k]).at(0) << "," <<  point_cloud2.at(ret_indexes[k]).at(1) << "," << point_cloud2.at(ret_indexes[k]).at(2) << ")" << std::endl;
//            nb_neighbors[i] = nb_pairs;
//        }
//        av_nb_pairs /= point_cloud1.size();
//        DebugOn("Min nb of Pairs = " << min_nb_pairs << endl);
//        DebugOn("Max nb of Pairs = " << max_nb_pairs << endl);
//        DebugOn("Average nb of Pairs = " << av_nb_pairs << endl);
//        param<> x1("x1"), x2("x2"), y1("y1"), y2("y2"), z1("z1"), z2("z2");
//        param<> x_uav1("x_uav1"), y_uav1("y_uav1"), z_uav1("z_uav1");
//        param<> x_uav2("x_uav2"), y_uav2("y_uav2"), z_uav2("z_uav2");
//            //        return 0;
//        bool solve_lidar_cube = false, solve_lidar_iter = !solve_lidar_cube;
//        int m = av_nb_pairs;
//            //            int m = 1;
//        vector<double> min_dist(point_cloud1.size(),numeric_limits<double>::max());
//        vector<int> nearest(point_cloud1.size());
//        vector<string> nearest_id(point_cloud1.size());
//        string i_str, j_str;
//        indices Pairs("Pairs"), cells("cells");
//        map<int,int> n2_map;
//        int idx1 = 0;
//        int idx2 = 0;
//        int nb_max_neigh = 1;
//        double dist_sq = 0;
//        if(solve_lidar_cube)
//            nb_max_neigh = m;
//        /* Keep points with neighbors >= m */
//        for (auto i = 0; i<point_cloud1.size(); i++) {
//            if(solve_lidar_iter)
//                nb_max_neigh = 1;
//            else
//                nb_max_neigh = m;
//            if(nb_neighbors[i]>=nb_max_neigh){
//                i_str = to_string(idx1+1);
//                x_uav1.add_val(i_str,uav1.at(i)[0]);
//                x1.add_val(i_str,point_cloud1.at(i)[0]);
//                y_uav1.add_val(i_str,uav1.at(i)[1]);
//                y1.add_val(i_str,point_cloud1.at(i)[1]);
//                z_uav1.add_val(i_str,uav1.at(i)[2]);
//                z1.add_val(i_str,point_cloud1.at(i)[2]);
//                if(solve_lidar_iter){
//                    nb_max_neigh = nb_neighbors[i];
//                }
//                for (auto j = 0; j<nb_max_neigh; j++) {
//                    auto k = neighbors[i].at(j);
//                    auto res = n2_map.find(k);
//                    if(res==n2_map.end()){
//                        n2_map[k] = idx2;
//                        j_str = to_string(idx2+1);
//                        x_uav2.add_val(j_str,uav2.at(k)[0]);
//                        x2.add_val(j_str,point_cloud2.at(k)[0]);
//                        y_uav2.add_val(j_str,uav2.at(k)[1]);
//                        y2.add_val(j_str,point_cloud2.at(k)[1]);
//                        z_uav2.add_val(j_str,uav2.at(k)[2]);
//                        z2.add_val(j_str,point_cloud2.at(k)[2]);
//                        idx2++;
//                    }
//                    else {
//                        j_str = to_string(res->second+1);
//                    }
//                    if(axis=="x")
//                        dist_sq = std::pow(point_cloud1.at(i)[1] - point_cloud2.at(k)[1],2) + std::pow(point_cloud1.at(i)[2] - point_cloud2.at(k)[2],2);
//                    else if(axis=="y")
//                        dist_sq = std::pow(point_cloud1.at(i)[0] - point_cloud2.at(k)[0],2) + std::pow(point_cloud1.at(i)[2] - point_cloud2.at(k)[2],2);
//                    else if(axis=="z")
//                        dist_sq = std::pow(point_cloud1.at(i)[0] - point_cloud2.at(k)[0],2) + std::pow(point_cloud1.at(i)[1] - point_cloud2.at(k)[1],2);
//                    else
//                        dist_sq = std::pow(point_cloud1.at(i)[0] - point_cloud2.at(k)[0],2) + std::pow(point_cloud1.at(i)[1] - point_cloud2.at(k)[1],2) + std::pow(point_cloud1.at(i)[2] - point_cloud2.at(k)[2],2);
//
//                    if(min_dist[i]>dist_sq){
//                        min_dist[i] = dist_sq;
//                        nearest[i] = k;
//                        nearest_id[i] = j_str;
//                    }
//
//                    if(solve_lidar_cube)
//                        Pairs.add(i_str+","+j_str);
//                }
//                idx1++;
//            }
//        }
//        idx1 = 0;
//        indices N1("N1"),N2("N2");
//        if(solve_lidar_iter){
//            for (auto i = 0; i<point_cloud1.size(); i++) {
//                if(nb_neighbors[i]>=1){
//                    i_str = to_string(idx1+1);
//                    j_str = nearest_id[i];
//                    cells.add(i_str+","+j_str);
//                    if(!N2.has_key(j_str))
//                        N2.add(j_str);
//                    idx1++;
//                }
//            }
//        }
//
//
//        int n1 = x1.get_dim();
//        int n2 = x2.get_dim();
//        DebugOn("n1 = " << n1 << endl);
//        DebugOn("n2 = " << n2 << endl);
//
//        N1 = range(1,n1);
//        if(solve_lidar_cube)
//            N2 = range(1,n2);
//        indices M("M");
//        M = range(1,m);
//
//        DebugOn("Total size of Pairs = " << Pairs.size() << endl);
//
//
//        indices NM("NM");
//        NM = indices(N1,M);
//        indices S1("S1"), S2("S2"), Sm1("Sm1"), S2m2("S2m2"), S3m1("S3m1"), Sm("Sm"), K("K");
//        S1 = indices(N1, range(1,1));
//        S2 = indices(N1, range(2,2));
//        Sm = indices(N1, range(m,m));
//            //            Sm1 = indices(N1, range(m-1,m-1));
//            //            S2m2 = indices(N1, range(2,m-2));
//            //            S3m1 = indices(N1, range(3,m-1));
//            //            K = indices(N1,range(2,m-1));
//
//
//
//        if (solve_lidar_iter) {
//            Model<> Lidar("Lidar");
//            var<> new_x1("new_x1"), new_y1("new_y1"), new_z1("new_z1");
//            var<> new_x2("new_x2"), new_y2("new_y2"), new_z2("new_z2");
//            var<> x_diff("x_diff", pos_), y_diff("y_diff", pos_), z_diff("z_diff", pos_);
//                //            var<> yaw1("yaw1", 0.25*pi/180, 0.25*pi/180), pitch1("pitch1", 0.5*pi/180, 0.5*pi/180), roll1("roll1", 0.7*pi/180, 0.7*pi/180);
//                //            var<> yaw1("yaw1", 0.25*pi/180, 0.25*pi/180), pitch1("pitch1", 0.9*pi/180, 0.9*pi/180), roll1("roll1", -1.45*pi/180, -1.45*pi/180);
//                //                var<> yaw1("yaw1", 0, 0), pitch1("pitch1", 0, 0), roll1("roll1", 0, 0);
//                //                var<> yaw1("yaw1", -0.5*pi/180, -0.5*pi/180), pitch1("pitch1", 0.9*pi/180, 0.9*pi/180), roll1("roll1", 1.375*pi/180, 1.375*pi/180);
//            var<> yaw1("yaw1", -0.1, 0.1), pitch1("pitch1", -0.1, 0.1), roll1("roll1", -0.1, 0.1);
//            var<> yaw2("yaw2", -0.1, 0.1), pitch2("pitch2", -0.1, 0.1), roll2("roll2", -0.1, 0.1);
//
//            Lidar.add(yaw1.in(R(1)),pitch1.in(R(1)),roll1.in(R(1)));
//            Lidar.add(yaw2.in(R(1)),pitch2.in(R(1)),roll2.in(R(1)));
//            Lidar.add(new_x1.in(N1), new_y1.in(N1), new_z1.in(N1));
//            Lidar.add(new_x2.in(N2), new_y2.in(N2), new_z2.in(N2));
//            Lidar.add(x_diff.in(cells), y_diff.in(cells), z_diff.in(cells));
//                //                Lidar.add(z_diff.in(cells));
//
//            Constraint<> Equal_pitch("Equal_pitch");
//            Equal_pitch += pitch1 - pitch2;
//            Lidar.add(Equal_pitch==0);
//
//            Constraint<> Opp_roll("Opp_roll");
//            Opp_roll += roll1 + roll2;
//            Lidar.add(Opp_roll==0);
//
//            Constraint<> Opp_yaw("Opp_yaw");
//            Opp_yaw += yaw1 + yaw2;
//            Lidar.add(Opp_yaw==0);
//            bool L2Norm = false;
//            bool L1Norm = !L2Norm;
//
//            if(L2Norm){
//                Constraint<> XNorm2("XNorm2");
//                XNorm2 += x_diff - pow(new_x1.from(cells) - new_x2.to(cells),2);
//                Lidar.add(XNorm2.in(cells)>=0);
//
//                Constraint<> YNorm2("YNorm2");
//                YNorm2 += y_diff - pow(new_y1.from(cells) - new_y2.to(cells),2);
//                Lidar.add(YNorm2.in(cells)>=0);
//
//                Constraint<> ZNorm2("ZNorm2");
//                ZNorm2 += z_diff - pow(new_z1.from(cells) - new_z2.to(cells),2);
//                Lidar.add(ZNorm2.in(cells)>=0);
//            }
//            if(L1Norm){
//                Constraint<> x_abs1("x_abs1");
//                x_abs1 += x_diff - (new_x1.from(cells) - new_x2.to(cells));
//                Lidar.add(x_abs1.in(cells)>=0);
//
//                Constraint<> x_abs2("x_abs2");
//                x_abs2 += x_diff - (new_x2.to(cells) - new_x1.from(cells));
//                Lidar.add(x_abs2.in(cells)>=0);
//
//                Constraint<> y_abs1("y_abs1");
//                y_abs1 += y_diff - (new_y1.from(cells) - new_y2.to(cells));
//                Lidar.add(y_abs1.in(cells)>=0);
//
//                Constraint<> y_abs2("y_abs2");
//                y_abs2 += y_diff - (new_y2.to(cells) - new_y1.from(cells));
//                Lidar.add(y_abs2.in(cells)>=0);
//
//                Constraint<> z_abs1("z_abs1");
//                z_abs1 += z_diff - (new_z1.from(cells) - new_z2.to(cells));
//                Lidar.add(z_abs1.in(cells)>=0);
//
//                Constraint<> z_abs2("z_abs2");
//                z_abs2 += z_diff - (new_z2.to(cells) - new_z1.from(cells));
//                Lidar.add(z_abs2.in(cells)>=0);
//            }
//
//            auto ids1 = yaw1.repeat_id(cells.size());
//            auto ids2 = yaw2.repeat_id(cells.size());
//
//            /* alpha = yaw_, beta = pitch_ and gamma = roll_ */
//            if(axis!="x"){
//                Constraint<> x_rot1("x_rot1");
//                x_rot1 += new_x1 - x_uav1.in(N1);
//                x_rot1 -= (x1.in(N1)-x_uav1.in(N1))*cos(yaw1.in(ids1))*cos(roll1.in(ids1)) + (y1.in(N1)-y_uav1.in(N1))*(cos(yaw1.in(ids1))*sin(roll1.in(ids1))*sin(pitch1.in(ids1)) - sin(yaw1.in(ids1))*cos(pitch1.in(ids1))) + (z1.in(N1)-z_uav1.in(N1))*(cos(yaw1.in(ids1))*sin(roll1.in(ids1))*cos(pitch1.in(ids1)) + sin(yaw1.in(ids1))*sin(pitch1.in(ids1)));
//                Lidar.add(x_rot1.in(N1)==0);
//
//
//                Constraint<> x_rot2("x_rot2");
//                x_rot2 += new_x2 - x_uav2.in(N2);
//                x_rot2 -= (x2.in(N2)-x_uav2.in(N2))*cos(yaw2.in(ids2))*cos(roll2.in(ids2)) + (y2.in(N2)-y_uav2.in(N2))*(cos(yaw2.in(ids2))*sin(roll2.in(ids2))*sin(pitch2.in(ids2)) - sin(yaw2.in(ids2))*cos(pitch2.in(ids2))) + (z2.in(N2)-z_uav2.in(N2))*(cos(yaw2.in(ids2))*sin(roll2.in(ids2))*cos(pitch2.in(ids2)) + sin(yaw2.in(ids2))*sin(pitch2.in(ids2)));
//                Lidar.add(x_rot2.in(N2)==0);
//            }
//
//            if(axis!="y"){
//                Constraint<> y_rot1("y_rot1");
//                y_rot1 += new_y1 - y_uav1.in(N1);
//                y_rot1 -= (x1.in(N1)-x_uav1.in(N1))*sin(yaw1.in(ids1))*cos(roll1.in(ids1)) + (y1.in(N1)-y_uav1.in(N1))*(sin(yaw1.in(ids1))*sin(roll1.in(ids1))*sin(pitch1.in(ids1)) + cos(yaw1.in(ids1))*cos(pitch1.in(ids1))) + (z1.in(N1)-z_uav1.in(N1))*(sin(yaw1.in(ids1))*sin(roll1.in(ids1))*cos(pitch1.in(ids1)) - cos(yaw1.in(ids1))*sin(pitch1.in(ids1)));
//                Lidar.add(y_rot1.in(N1)==0);
//
//                Constraint<> y_rot2("y_rot2");
//                y_rot2 += new_y2 - y_uav2.in(N2);
//                y_rot2 -= (x2.in(N2)-x_uav2.in(N2))*sin(yaw2.in(ids2))*cos(roll2.in(ids2)) + (y2.in(N2)-y_uav2.in(N2))*(sin(yaw2.in(ids2))*sin(roll2.in(ids2))*sin(pitch2.in(ids2)) + cos(yaw2.in(ids2))*cos(pitch2.in(ids2))) + (z2.in(N2)-z_uav2.in(N2))*(sin(yaw2.in(ids2))*sin(roll2.in(ids2))*cos(pitch2.in(ids2)) - cos(yaw2.in(ids2))*sin(pitch2.in(ids2)));
//                Lidar.add(y_rot2.in(N2)==0);
//            }
//
//            if(axis!="z"){
//                Constraint<> z_rot1("z_rot1");
//                z_rot1 += new_z1 - z_uav1.in(N1);
//                z_rot1 -= (x1.in(N1)-x_uav1.in(N1))*sin(-1*roll1.in(ids1)) + (y1.in(N1)-y_uav1.in(N1))*(cos(roll1.in(ids1))*sin(pitch1.in(ids1))) + (z1.in(N1)-z_uav1.in(N1))*(cos(roll1.in(ids1))*cos(pitch1.in(ids1)));
//                Lidar.add(z_rot1.in(N1)==0);
//
//
//                Constraint<> z_rot2("z_rot2");
//                z_rot2 += new_z2 - z_uav2.in(N2);
//                z_rot2 -= (x2.in(N2)-x_uav2.in(N2))*sin(-1*roll2.in(ids2)) + (y2.in(N2)-y_uav2.in(N2))*(cos(roll2.in(ids2))*sin(pitch2.in(ids2))) + (z2.in(N2)-z_uav2.in(N2))*(cos(roll2.in(ids2))*cos(pitch2.in(ids2)));
//                Lidar.add(z_rot2.in(N2)==0);
//            }
//
//                //    M.min(sum(z_diff)/nb_overlap);
//
//                //        M.min(sum(z_diff));
//            if(axis == "full")
//                Lidar.min(sum(x_diff)/cells.size() + sum(y_diff)/cells.size() + sum(z_diff)/cells.size());
//            else if(axis == "x"){
//                Lidar.min(sum(y_diff)/cells.size() + sum(z_diff)/cells.size());
//                roll1.set_lb(0);
//                roll1.set_ub(0);
//                yaw1.set_lb(0);
//                yaw1.set_ub(0);
//                    //                x1.set_val(0);
//                    //                x2.set_val(0);
//            }
//            else if (axis == "y") {
//                Lidar.min(sum(x_diff)/cells.size() + sum(z_diff)/cells.size());
//                yaw1.set_lb(0);
//                yaw1.set_ub(0);
//                pitch1.set_lb(0);
//                pitch1.set_ub(0);
//                    //                y1.set_val(0);
//                    //                y2.set_val(0);
//            }
//            else if (axis == "z") {
//                Lidar.min(sum(x_diff)/cells.size() + sum(y_diff)/cells.size());
//                pitch1.set_lb(0);
//                pitch1.set_ub(0);
//                roll1.set_lb(0);
//                roll1.set_ub(0);
//                    //                z1.set_val(0);
//                    //                z2.set_val(0);
//            }
//            else if (axis == "only_x")
//                Lidar.min(sum(x_diff)/cells.size());
//            else if (axis == "only_y")
//                Lidar.min(sum(y_diff)/cells.size());
//            else if (axis == "only_z")
//                Lidar.min(sum(z_diff)/cells.size());
//
//                //                Lidar.min(sum(x_diff)/cells.size() + sum(y_diff)/cells.size() + sum(z_diff)/cells.size());
//
//                //    M.print();
//
//            solver<> S(Lidar,ipopt);
//            S.run();
//
//
//                //        for (int i = 0; i<500; i++) {
//                //            pre_x.add_val(x_rot1.eval(i));
//                //            pre_y.add_val(y_rot1.eval(i));
//                //            pre_z.add_val(z_rot1.eval(i));
//                //            x_uav.add_val(x_uav1.eval(i));
//                //            y_uav.add_val(y_uav1.eval(i));
//                //            z_uav.add_val(z_uav1.eval(i));
//                //        }
//                //        for (int i = 0; i<500; i++) {
//                //            pre_x.add_val(x_rot2.eval(i));
//                //            pre_y.add_val(y_rot2.eval(i));
//                //            pre_z.add_val(z_rot2.eval(i));
//                //            x_uav.add_val(x_uav2.eval(i));
//                //            y_uav.add_val(y_uav2.eval(i));
//                //            z_uav.add_val(z_uav2.eval(i));
//                //        }
//                //    M.print_solution();
//
//            DebugOn("Pitch1 = " << pitch1.eval()*180/pi << endl);
//            DebugOn("Roll1 = " << roll1.eval()*180/pi << endl);
//            DebugOn("Yaw1 = " << yaw1.eval()*180/pi << endl);
//            DebugOn("Pitch2 = " << pitch2.eval()*180/pi << endl);
//            DebugOn("Roll2 = " << roll2.eval()*180/pi << endl);
//            DebugOn("Yaw2 = " << yaw2.eval()*180/pi << endl);
//            roll_1 = roll1.eval()*180/pi;
//            pitch_1 = pitch1.eval()*180/pi;
//            yaw_1 = yaw1.eval()*180/pi;
//        }
//        else if (solve_lidar_cube) {
//            Model<> Lidar("Lidar");
//            var<> new_x1("new_x1"), new_y1("new_y1"), new_z1("new_z1");
//            var<> new_x2("new_x2"), new_y2("new_y2"), new_z2("new_z2");
//                //            var<> x_diff("x_diff", pos_), y_diff("y_diff", pos_), z_diff("z_diff", pos_);
//            var<> mu("mu", pos_), mu_k("mu_k", pos_), delta("delta", pos_);
//                //                            var<> yaw1("yaw1", -0.5*pi/180, -0.5*pi/180), pitch1("pitch1", 0.9*pi/180, 0.9*pi/180), roll1("roll1", 1.375*pi/180, 1.375*pi/180);
//            var<> yaw1("yaw1", -0.1, 0.1), pitch1("pitch1", -0.1, 0.1), roll1("roll1", -0.1, 0.1);
//                //                var<> yaw1("yaw1", 0.25*pi/180, 0.25*pi/180), pitch1("pitch1", 0.9*pi/180, 0.9*pi/180), roll1("roll1", -1.45*pi/180, -1.45*pi/180);
//                //                var<> yaw1("yaw1", 0, 0), pitch1("pitch1", -0.5778*pi/180, -0.5778*pi/180), roll1("roll1", -1.44581*pi/180, -1.44581*pi/180);
//                //                var<> yaw1("yaw1", 0, 0), pitch1("pitch1", -0.573231*pi/180, -0.573231*pi/180), roll1("roll1", -1.45338*pi/180, -1.45338*pi/180);
//                //                var<> yaw1("yaw1", 0.0249847*pi/180, 0.0249847*pi/180), pitch1("pitch1", -0.507086*pi/180, -0.507086*pi/180), roll1("roll1", -1.3698*pi/180, -1.3698*pi/180);
//                //                var<> yaw1("yaw1", 0, 0), pitch1("pitch1", 0, 0), roll1("roll1", 0, 0);
//            var<> yaw2("yaw2", -0.1, 0.1), pitch2("pitch2", -0.1, 0.1), roll2("roll2", -0.1, 0.1);
//                //            yaw1 = -0.5*pi/180;
//                //            pitch1 = 0.9*pi/180;
//                //            roll1 = 1.375*pi/180;
//            Lidar.add(yaw1.in(R(1)),pitch1.in(R(1)),roll1.in(R(1)));
//            Lidar.add(yaw2.in(R(1)),pitch2.in(R(1)),roll2.in(R(1)));
//            Lidar.add(new_x1.in(N1), new_y1.in(N1), new_z1.in(N1));
//            Lidar.add(new_x2.in(N2), new_y2.in(N2), new_z2.in(N2));
//                //            Lidar.add(x_diff.in(NM), y_diff.in(NM), z_diff.in(NM));
//            Lidar.add(delta.in(NM));
//            Lidar.add(mu.in(N1));
//            Lidar.add(mu_k.in(K));
//
//            Constraint<> Equal_pitch("Equal_pitch");
//            Equal_pitch += pitch1 - pitch2;
//            Lidar.add(Equal_pitch==0);
//
//            Constraint<> Opp_roll("Opp_roll");
//            Opp_roll += roll1 + roll2;
//            Lidar.add(Opp_roll==0);
//
//            Constraint<> Opp_yaw("Opp_yaw");
//            Opp_yaw += yaw1 + yaw2;
//            Lidar.add(Opp_yaw==0);
//
//            Constraint<> Mu_2("Mu_2");
//            Mu_2 += mu_k.in(S2) - gravity::min(delta.in(S1), delta.in(S2));
//            Lidar.add(Mu_2.in(S2)==0);
//
//            Constraint<> Mu_k("Mu_k");
//            Mu_k += mu_k.in(S3m1) - gravity::min(mu_k.in(S2m2), delta.in(S3m1));
//            Lidar.add(Mu_k.in(S3m1)==0);
//
//            Constraint<> Mu("Mu");
//            Mu += mu - gravity::min(mu_k.in(Sm1), delta.in(Sm));
//            Lidar.add(Mu.in(N1)==0);
//
//                //                            Constraint<> Norm2("Norm2");
//                //                            Norm2 += delta - pow(new_x1.in_ignore_ith(1, 1, Pairs) - new_x2.in_ignore_ith(0, 1, Pairs),2) - pow(new_y1.in_ignore_ith(1, 1, Pairs) - new_y2.in_ignore_ith(0, 1, Pairs),2) - pow(new_z1.in_ignore_ith(1, 1, Pairs) - new_z2.in_ignore_ith(0, 1, Pairs),2);
//                //                            Lidar.add(Norm2.in(Pairs)==0);
//
//            Constraint<> Norm1("Norm1");
//            Norm1 += delta - abs(new_x1.in_ignore_ith(1, 1, Pairs) - new_x2.in_ignore_ith(0, 1, Pairs)) - abs(new_y1.in_ignore_ith(1, 1, Pairs) - new_y2.in_ignore_ith(0, 1, Pairs)) - abs(new_z1.in_ignore_ith(1, 1, Pairs) - new_z2.in_ignore_ith(0, 1, Pairs));
//                //                Norm1 += delta - abs(new_z1.in_ignore_ith(1, 1, Pairs) - new_z2.in_ignore_ith(0, 1, Pairs));
//            Lidar.add(Norm1.in(Pairs)==0);
//
//
//                //            Constraint<> z_abs1("z_abs1");
//                //            z_abs1 += z_diff - (new_z1.in_ignore_ith(1, 1, Pairs) - new_z2.in_ignore_ith(0, 1, Pairs));
//                //            Lidar.add(z_abs1.in(Pairs)>=0);
//                //
//                //            Constraint<> z_abs2("z_abs2");
//                //            z_abs2 += z_diff - (new_z2.in_ignore_ith(1, 0, Pairs) - new_z1.in_ignore_ith(1, 1, Pairs));
//                //            Lidar.add(z_abs2.in(Pairs)>=0);
//                //
//                //            Constraint<> x_abs1("x_abs1");
//                //            x_abs1 += x_diff - (new_x1.in_ignore_ith(1, 1, Pairs) - new_x2.in_ignore_ith(0, 1, Pairs));
//                //            Lidar.add(x_abs1.in(Pairs)>=0);
//                //
//                //            Constraint<> x_abs2("x_abs2");
//                //            x_abs2 += x_diff - (new_x2.in_ignore_ith(0, 1, Pairs) - new_x1.in_ignore_ith(1, 1, Pairs));
//                //            Lidar.add(x_abs2.in(Pairs)>=0);
//                //
//                //
//                //            Constraint<> y_abs1("y_abs1");
//                //            y_abs1 += y_diff - (new_y1.in_ignore_ith(1, 1, Pairs) - new_y2.in_ignore_ith(0, 1, Pairs));
//                //            Lidar.add(y_abs1.in(Pairs)>=0);
//                //
//                //            Constraint<> y_abs2("y_abs2");
//                //            y_abs2 += y_diff - (new_y2.in_ignore_ith(0, 1, Pairs) - new_y1.in_ignore_ith(1, 1, Pairs));
//                //            Lidar.add(y_abs2.in(Pairs)>=0);
//
//            auto ids1 = yaw1.repeat_id(n1);
//
//            /* alpha = yaw_, beta = pitch_ and gamma = roll_ */
//            Constraint<> x_rot1("x_rot1");
//            x_rot1 += new_x1 - x_uav1.in(N1);
//            x_rot1 -= (x1.in(N1)-x_uav1.in(N1))*cos(yaw1.in(ids1))*cos(roll1.in(ids1)) + (y1.in(N1)-y_uav1.in(N1))*(cos(yaw1.in(ids1))*sin(roll1.in(ids1))*sin(pitch1.in(ids1)) - sin(yaw1.in(ids1))*cos(pitch1.in(ids1))) + (z1.in(N1)-z_uav1.in(N1))*(cos(yaw1.in(ids1))*sin(roll1.in(ids1))*cos(pitch1.in(ids1)) + sin(yaw1.in(ids1))*sin(pitch1.in(ids1)));
//            Lidar.add(x_rot1.in(N1)==0);
//
//            auto ids2 = yaw2.repeat_id(n2);
//
//            Constraint<> x_rot2("x_rot2");
//            x_rot2 += new_x2 - x_uav2.in(N2);
//            x_rot2 -= (x2.in(N2)-x_uav2.in(N2))*cos(yaw2.in(ids2))*cos(roll2.in(ids2)) + (y2.in(N2)-y_uav2.in(N2))*(cos(yaw2.in(ids2))*sin(roll2.in(ids2))*sin(pitch2.in(ids2)) - sin(yaw2.in(ids2))*cos(pitch2.in(ids2))) + (z2.in(N2)-z_uav2.in(N2))*(cos(yaw2.in(ids2))*sin(roll2.in(ids2))*cos(pitch2.in(ids2)) + sin(yaw2.in(ids2))*sin(pitch2.in(ids2)));
//            Lidar.add(x_rot2.in(N2)==0);
//
//
//            Constraint<> y_rot1("y_rot1");
//            y_rot1 += new_y1 - y_uav1.in(N1);
//            y_rot1 -= (x1.in(N1)-x_uav1.in(N1))*sin(yaw1.in(ids1))*cos(roll1.in(ids1)) + (y1.in(N1)-y_uav1.in(N1))*(sin(yaw1.in(ids1))*sin(roll1.in(ids1))*sin(pitch1.in(ids1)) + cos(yaw1.in(ids1))*cos(pitch1.in(ids1))) + (z1.in(N1)-z_uav1.in(N1))*(sin(yaw1.in(ids1))*sin(roll1.in(ids1))*cos(pitch1.in(ids1)) - cos(yaw1.in(ids1))*sin(pitch1.in(ids1)));
//            Lidar.add(y_rot1.in(N1)==0);
//
//            Constraint<> y_rot2("y_rot2");
//            y_rot2 += new_y2 - y_uav2.in(N2);
//            y_rot2 -= (x2.in(N2)-x_uav2.in(N2))*sin(yaw2.in(ids2))*cos(roll2.in(ids2)) + (y2.in(N2)-y_uav2.in(N2))*(sin(yaw2.in(ids2))*sin(roll2.in(ids2))*sin(pitch2.in(ids2)) + cos(yaw2.in(ids2))*cos(pitch2.in(ids2))) + (z2.in(N2)-z_uav2.in(N2))*(sin(yaw2.in(ids2))*sin(roll2.in(ids2))*cos(pitch2.in(ids2)) - cos(yaw2.in(ids2))*sin(pitch2.in(ids2)));
//            Lidar.add(y_rot2.in(N2)==0);
//
//            Constraint<> z_rot1("z_rot1");
//            z_rot1 += new_z1 - z_uav1.in(N1);
//            z_rot1 -= (x1.in(N1)-x_uav1.in(N1))*sin(-1*roll1.in(ids1)) + (y1.in(N1)-y_uav1.in(N1))*(cos(roll1.in(ids1))*sin(pitch1.in(ids1))) + (z1.in(N1)-z_uav1.in(N1))*(cos(roll1.in(ids1))*cos(pitch1.in(ids1)));
//            Lidar.add(z_rot1.in(N1)==0);
//
//
//            Constraint<> z_rot2("z_rot2");
//            z_rot2 += new_z2 - z_uav2.in(N2);
//            z_rot2 -= (x2.in(N2)-x_uav2.in(N2))*sin(-1*roll2.in(ids2)) + (y2.in(N2)-y_uav2.in(N2))*(cos(roll2.in(ids2))*sin(pitch2.in(ids2))) + (z2.in(N2)-z_uav2.in(N2))*(cos(roll2.in(ids2))*cos(pitch2.in(ids2)));
//            Lidar.add(z_rot2.in(N2)==0);
//
//                //            Lidar.min(sum(mu) + 1e2*pow(yaw1,2) + 1e2*pow(roll1,2) + 1e2*pow(pitch1,2));
//                //            Lidar.min(1e3*sum(mu) + (pow(yaw1,2) + pow(roll1,2) + pow(pitch1,2)));
//            Lidar.min(sum(mu));
//
//                //            Lidar.print();
//                //            Lidar.initialize_zero();
//                //            return 0;
//            solver<> S(Lidar,ipopt);
//                //            S.set_option("tol", 1e-10);
//                //            S.run(5,1e-10);
//            S.run();
//            DebugOn("Pitch1 = " << pitch1.eval()*180/pi << endl);
//            DebugOn("Roll1 = " << roll1.eval()*180/pi << endl);
//            DebugOn("Yaw1 = " << yaw1.eval()*180/pi << endl);
//            DebugOn("Pitch2 = " << pitch2.eval()*180/pi << endl);
//            DebugOn("Roll2 = " << roll2.eval()*180/pi << endl);
//            DebugOn("Yaw2 = " << yaw2.eval()*180/pi << endl);
//            roll_1 = roll1.eval()*180/pi;
//            pitch_1 = pitch1.eval()*180/pi;
//            yaw_1 = yaw1.eval()*180/pi;
//                //            return 0;
//        }
//    }
//    return {roll_1, pitch_1, yaw_1};
//}

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
    var<> yaw("yaw", yaw_min, yaw_max), pitch("pitch", pitch_min, pitch_max), roll("roll", roll_min, roll_max);
    func<> r11 = cos(yaw)*cos(roll);
    func<> r12 = cos(yaw)*sin(roll)*sin(pitch) - sin(yaw)*cos(pitch);
    func<> r13 = cos(yaw)*sin(roll)*cos(pitch) + sin(yaw)*sin(pitch);
    func<> r21 = sin(yaw)*cos(roll);
    func<> r22 = sin(yaw)*sin(roll)*sin(pitch) + cos(yaw)*cos(pitch);
    func<> r23 = sin(yaw)*sin(roll)*cos(pitch) - cos(yaw)*sin(pitch);
    func<> r31 = sin(-1*roll);
    func<> r32 = cos(roll)*sin(pitch);
    func<> r33 = cos(roll)*cos(pitch);
    
    
    var<> theta11("theta11",  std::max(-1.,r11._range->first), std::min(1.,r11._range->second)), theta12("theta12", std::max(-1.,r12._range->first), std::min(1.,r12._range->second)), theta13("theta13", std::max(-1.,r13._range->first), std::min(1.,r13._range->second));
    var<> theta21("theta21", std::max(-1.,r21._range->first), std::min(1.,r21._range->second)), theta22("theta22", std::max(-1.,r22._range->first), std::min(1.,r22._range->second)), theta23("theta23", std::max(-1.,r23._range->first), std::min(1.,r23._range->second));
    var<> theta31("theta31", std::max(-1.,r31._range->first), std::min(1.,r31._range->second)), theta32("theta32", std::max(-1.,r32._range->first), std::min(1.,r32._range->second)), theta33("theta33", std::max(-1.,r33._range->first), std::min(1.,r33._range->second));

    goicp.R_init.val[0][0]=(std::max(-1.,r11._range->first)+ std::min(1.,r11._range->second))/2.0;
    goicp.R_init.val[0][1]=(std::max(-1.,r12._range->first)+ std::min(1.,r12._range->second))/2.0;
    goicp.R_init.val[0][2]=(std::max(-1.,r13._range->first)+ std::min(1.,r13._range->second))/2.0;
    goicp.R_init.val[1][0]=(std::max(-1.,r21._range->first)+ std::min(1.,r21._range->second))/2.0;
    goicp.R_init.val[1][1]=(std::max(-1.,r22._range->first)+ std::min(1.,r22._range->second))/2.0;
    goicp.R_init.val[1][2]=(std::max(-1.,r23._range->first)+ std::min(1.,r23._range->second))/2.0;
    goicp.R_init.val[2][0]=(std::max(-1.,r31._range->first)+ std::min(1.,r31._range->second))/2.0;
    goicp.R_init.val[2][1]=(std::max(-1.,r32._range->first)+ std::min(1.,r32._range->second))/2.0;
    goicp.R_init.val[2][2]=(std::max(-1.,r33._range->first)+ std::min(1.,r33._range->second))/2.0;
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
    DebugOn("ICP error "<<error<<endl);
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

#ifdef USE_PCL
//vector<vector<double>> filter_pcl(vector<vector<double>> model){
//    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZ>);
//        //pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_filtered (new pcl::PointCloud<pcl::PointXYZ>);
//    for (int i = 0; i<ext_model.size(); i++) {
//        pcl::PointXYZ p;
//        p.x = model[i][0];
//        p.y = model[i][1];
//        p.z = model[i][2];
//        cloud->push_back(p);
//    }
//
//    return model;
//
//}
void save_feature_file(const string& filename, const pcl::PointCloud<pcl::PointNormal>::Ptr& augmented_cloud, const pcl::PointCloud<pcl::FPFHSignature33>::Ptr& cloud_features){
    string fname = filename+"_features.bin";
    FILE* fid = fopen(fname.c_str(), "wb");
    int nV = augmented_cloud->size(), nDim = 33;
    fwrite(&nV, sizeof(int), 1, fid);
    fwrite(&nDim, sizeof(int), 1, fid);
    for (int v = 0; v < nV; v++) {
        const pcl::PointNormal &pt = augmented_cloud->points[v];
        float xyz[3] = {pt.x, pt.y, pt.z};
        fwrite(xyz, sizeof(float), 3, fid);
        const pcl::FPFHSignature33 &feature = cloud_features->points[v];
        fwrite(feature.histogram, sizeof(float), 33, fid);
        if(v<10){
            DebugOn("For point: (" << pt.x << "," << pt.y << "," << pt.z << "): ");
            for (int i = 0; i<33; i++) {
                DebugOn(feature.histogram[i] << " , ");
            }
            DebugOn(endl);
        }
    }
    DebugOn(endl);
    fclose(fid);
}

pair<pcl::PointCloud<pcl::PointNormal>::Ptr,pcl::PointCloud<pcl::FPFHSignature33>::Ptr> compute_features(const vector<vector<double>>& ext_model) {
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZ>);
    
    for (int i = 0; i<ext_model.size(); i++) {
        pcl::PointXYZ p;
        p.x = ext_model[i][0];
        p.y = ext_model[i][1];
        p.z = ext_model[i][2];
        cloud->push_back(p);
    }
    
    // Create the normal estimation class, and pass the input dataset to it
    pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> ne;
    ne.setInputCloud (cloud);
    
    // Create an empty kdtree representation, and pass it to the normal estimation object.
    // Its content will be filled inside the object, based on the given input dataset (as no other search surface is given).
    pcl::search::KdTree<pcl::PointXYZ>::Ptr tree (new pcl::search::KdTree<pcl::PointXYZ> ());
    ne.setSearchMethod (tree);
    
    // Output datasets
    pcl::PointCloud<pcl::Normal>::Ptr cloud_normals (new pcl::PointCloud<pcl::Normal>);
    
    // Use all neighbors in a sphere of radius 0.1
    ne.setRadiusSearch (0.1);
    
    // Compute the normals and make sure they are finite.
    ne.compute (*cloud_normals);
    for (int i = 0; i < cloud_normals->size(); i++)
    {
        if (!pcl::isFinite<pcl::Normal>((*cloud_normals)[i]))
        {
            PCL_WARN("normals[%d] is not finite\n", i);
        }
    }
    // cloud_normals->size () should have the same size as the input cloud->size ()*
    // Assume a point cloud with normal is given as
    pcl::PointCloud<pcl::PointNormal>::Ptr augmented_cloud(new  pcl::PointCloud<pcl::PointNormal>());
    augmented_cloud->resize(cloud_normals->size());
    for(std::size_t i = 0; i<cloud_normals->size(); ++i)
    {
        (*augmented_cloud)[i].x = (*cloud)[i].x;
        (*augmented_cloud)[i].y = (*cloud)[i].y;
        (*augmented_cloud)[i].z = (*cloud)[i].z;
        (*augmented_cloud)[i].normal_x = (*cloud_normals)[i].normal_x;
        (*augmented_cloud)[i].normal_y = (*cloud_normals)[i].normal_y;
        (*augmented_cloud)[i].normal_z = (*cloud_normals)[i].normal_z;
    }
    // Create the FPFH estimation class, and pass the input dataset+normals to it
    pcl::FPFHEstimation<pcl::PointXYZ, pcl::Normal, pcl::FPFHSignature33> fpfh;
    pcl::PointCloud<pcl::FPFHSignature33>::Ptr cloud_features(new pcl::PointCloud<pcl::FPFHSignature33>());
    
    fpfh.setRadiusSearch(1);
    fpfh.setInputCloud(cloud);
    fpfh.setInputNormals(cloud_normals);
    fpfh.compute(*cloud_features);
    return {augmented_cloud,cloud_features};
}

#endif

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

indices get_valid_pairs(vector<vector<double>>& point_cloud_data, vector<vector<double>>& point_cloud_model, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double shift_min_x, double shift_max_x, double shift_min_y, double shift_max_y, double shift_min_z, double shift_max_z, param<>& norm_x, param<>& norm_y, param<>& norm_z,  param<>& intercept, const vector<double>& model_voronoi_out_radius, vector<int>& new_model_pts, indices& new_model_ids, bool norm1){
    indices valid_cells("valid_cells");
    size_t nm = point_cloud_model.size(), nd = point_cloud_data.size();
    set<int> unique_model_pts;
    vector<double> zeros = {0,0,0};
    for (int i = 0; i<nd; i++) {
        auto max_dist = get_max_dist(roll_min, roll_max, pitch_min, pitch_max, yaw_min, yaw_max, shift_min_x, shift_max_x, shift_min_y, shift_max_y, shift_min_z, shift_max_z, point_cloud_data[i], zeros, norm1);
        for(auto j=0;j<nm;j++){
            auto voronoi_radius = model_voronoi_out_radius[j];
            auto dist = std::sqrt(std::pow(point_cloud_data[i][0] - point_cloud_model[j][0],2) + std::pow(point_cloud_data[i][1] - point_cloud_model[j][1],2) + std::pow(point_cloud_data[i][2] - point_cloud_model[j][2],2));
            if(dist>max_dist+voronoi_radius){
                DebugOn("incompatible pair: (" << i+1 << "," << j+1 << ")\n");
            }
            else{
                valid_cells.insert(to_string(i+1) +"," + to_string(j+1));
                unique_model_pts.insert(j);
            }
        }
    }
    for (int m_id: unique_model_pts) {
        new_model_pts.push_back(m_id);
        new_model_ids.insert(to_string(m_id+1));
    }
    DebugOn("Number of valid cells = " << valid_cells.size() << endl);
    DebugOn("Number of discarded pairs = " << nm*nd - valid_cells.size() << endl);
    return valid_cells;
}

//indices get_valid_pairs(vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double shift_min_x, double shift_max_x, double shift_min_y, double shift_max_y, double shift_min_z, double shift_max_z, vector<vector<vector<double>>> model_voronoi_normals, vector<vector<double>> model_face_intercept, const vector<double>& model_voronoi_out_radius, vector<vector<vector<double>>> model_voronoi_vertices, bool norm1){
//    norm1=false;
//    indices valid_cells("valid_cells");
//    size_t nm = point_cloud_model.size(), nd = point_cloud_data.size();
//    vector<double> zeros = {0,0,0};
//    double max_dist_i, min_dist_ik;
//    for (int i = 0; i<nd; i++) {
//         max_dist_i = get_max_dist(roll_min, roll_max, pitch_min, pitch_max, yaw_min, yaw_max, shift_min_x, shift_max_x, shift_min_y, shift_max_y, shift_min_z, shift_max_z, point_cloud_data[i], zeros, norm1);
//        auto xi=point_cloud_data.at(i)[0];
//        auto yi=point_cloud_data.at(i)[1];
//        auto zi=point_cloud_data.at(i)[2];
//        for(auto j=0;j<nm;j++){
//            bool center_voronoi_feasible=true;
//            bool vertex_intersect_circle=false;
//            auto normals=model_voronoi_normals.at(j);
//            auto intercepts=model_face_intercept.at(j);
//            auto vertex=model_voronoi_vertices.at(j);
//            vector<bool> ij_inter;
//            vector<double> feasible, min_dist_ik;
//            feasible.resize(intercepts.size());
//            min_dist_ik.resize(intercepts.size());
//            ij_inter.resize(intercepts.size(), false);
//            int count_ik_true=0;
//            for(auto k=0;k<intercepts.size();k++){
//                auto a=normals[k][0];
//                auto b=normals[k][1];
//                auto c=normals[k][2];
//                auto d=intercepts[k];
//                feasible[k]=a*xi+b*yi+c*zi+d;
//                if(feasible[k]>=1e-9){
//                    center_voronoi_feasible=false;
//                }
//                min_dist_ik[k]=std::abs(a*xi+b*yi+c*zi+d)/(a*a+b*b+c*c);
//                if(min_dist_ik[k]<=max_dist_i){
//                    ij_inter.at(k)=true;
//                    count_ik_true++;
//                }
//            }
//            if(!center_voronoi_feasible){
//                for(auto k=0;k<vertex.size();k++){
//                    auto a=vertex[k][0];
//                    auto b=vertex[k][1];
//                    auto c=vertex[k][2];
//                    if((a-xi)*(a-xi)+(b-yi)*(b-yi)+(c-zi)*(c-zi)<=(max_dist_i*max_dist_i)){
//                        vertex_intersect_circle=true;
//                    }
//                }
//            }
//
//            if(count_ik_true>=2){
//                valid_cells.insert(to_string(i+1) +"," + to_string(j+1));
//            }
//            else{
//                DebugOn("incompatible pair: (" << i+1 << "," << j+1 << ")\n");
//            }
//        }
//    }
//    DebugOn("Number of valid cells = " << valid_cells.size() << endl);
//    DebugOn("Number of discarded pairs = " << nm*nd - valid_cells.size() << endl);
//    return valid_cells;
//}


indices preprocess(vector<vector<double>> point_cloud_data, vector<vector<double>> point_cloud_model, double angle_max_deg, double shift_min_x,double shift_max_x,double shift_min_y,double shift_max_y,double shift_min_z,double shift_max_z, vector<vector<vector<double>>> model_voronoi_normals, vector<vector<double>> model_face_intercept){
    double angle_max = angle_max_deg*pi/180.;
    var<> yaw("yaw", -angle_max, angle_max), pitch("pitch", -angle_max, angle_max), roll("roll", -angle_max, angle_max);
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
    
    
    var<> theta11("theta11",  r11._range->first, r11._range->second), theta12("theta12", r12._range->first, r12._range->second), theta13("theta13", r13._range->first, r13._range->second);
    var<> theta21("theta21", r21._range->first, r21._range->second), theta22("theta22", r22._range->first, r22._range->second), theta23("theta23", r23._range->first, r23._range->second);
    var<> theta31("theta31", r31._range->first, r31._range->second), theta32("theta32", r32._range->first, r32._range->second), theta33("theta33", r33._range->first, r33._range->second);
    var<> x_shift("x_shift", shift_min_x, shift_max_x), y_shift("y_shift", shift_min_y, shift_max_y), z_shift("z_shift", shift_min_z, shift_max_z);
    gravity::indices cells("cells");
    double x_lb = 0, y_lb = 0, z_lb = 0, x1_i = 0, y1_i = 0, z1_i = 0;
    shared_ptr<pair<double,double>> new_x1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> new_y1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> new_z1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> x1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> y1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> z1_bounds = make_shared<pair<double,double>>();
    int nd=point_cloud_data.size();
    int nm=point_cloud_model.size();
    int missed=0;
    for (int i = 0; i<nd; i++) {
        string i_str = to_string(i+1);
        //        auto bounds = get_min_max(angle_max, point_cloud_data[i], zeros);
        x1_bounds->first = point_cloud_data.at(i)[0];
        x1_bounds->second = point_cloud_data.at(i)[0];
        y1_bounds->first = point_cloud_data.at(i)[1];
        y1_bounds->second = point_cloud_data.at(i)[1];
        z1_bounds->first = point_cloud_data.at(i)[2];
        z1_bounds->second = point_cloud_data.at(i)[2];
        auto x_range  = get_product_range(x1_bounds, theta11._range);
        auto y_range  = get_product_range(y1_bounds, theta12._range);
        auto z_range  = get_product_range(z1_bounds, theta13._range);
        *new_x1_bounds = {x_range->first + y_range->first + z_range->first + shift_min_x,
            x_range->second + y_range->second + z_range->second+ shift_max_x};
        x_range  = get_product_range(x1_bounds, theta21._range);
        y_range  = get_product_range(y1_bounds, theta22._range);
        z_range  = get_product_range(z1_bounds, theta23._range);
        *new_y1_bounds = {x_range->first + y_range->first + z_range->first + shift_min_y,
            x_range->second + y_range->second + z_range->second+ shift_max_y};
        x_range  = get_product_range(x1_bounds, theta31._range);
        y_range  = get_product_range(y1_bounds, theta32._range);
        z_range  = get_product_range(z1_bounds, theta33._range);
        *new_z1_bounds = {x_range->first + y_range->first + z_range->first + shift_min_z,
            x_range->second + y_range->second + z_range->second+ shift_max_z};
        for(auto j=0;j<nm;j++){
            auto voro_model=make_shared<Model<double>>("voro_model");
            var<> x("x", new_x1_bounds->first, new_x1_bounds->second), y("y", new_y1_bounds->first, new_y1_bounds->second),z("z", new_z1_bounds->first, new_z1_bounds->second);
            voro_model->add(x.in(R(1)), y.in(R(1)), z.in(R(1)));
            bool intersect=false;
            param<> a("a"), b("b"), c("c"), d("d");
            auto normals=model_voronoi_normals.at(j);
            auto intercepts=model_face_intercept.at(j);
            for(auto k=0;k<intercepts.size();k++){
                a.add_val(to_string(k), normals[k][0]);
                b.add_val(to_string(k), normals[k][1]);
                c.add_val(to_string(k), normals[k][2]);
                d.add_val(to_string(k), intercepts[k]);
            }
            Constraint<> region("region");
            region=a*x+b*y+c*z+d;
            voro_model->add(region<=0);
            
            func<> obj=0;
            voro_model->min(obj);
            // voro_model->print();
            solver<> s(voro_model,ipopt);
            s.run(0, 1e-6);
            if(voro_model->_status!=0){
                DebugOn("Missed pair: (" << i+1 << "," << j+1 << ")\n");
                missed++;
            }
            else{
                auto key=to_string(i+1)+","+to_string(j+1);
                cells.add(key);
            }
        }
    }
    DebugOff("Num "<<missed<<endl);
    return (cells);
}
indices preprocess_QP(const vector<vector<double>>& point_cloud_data, const vector<vector<double>>& point_cloud_model, const indices& old_cells, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double shift_min_x, double shift_max_x, double shift_min_y, double shift_max_y, double shift_min_z, double shift_max_z, const vector<vector<vector<double>>>& model_voronoi_normals, const vector<vector<double>>& model_face_intercept, const vector<vector<vector<double>>>& model_voronoi_vertices, vector<int>& new_model_pts, indices& new_model_ids, param<>& dist_cost, double upper_bound, int nb_total_threads){
    nb_total_threads=4;
    double time_start = get_wall_time();
    
    shared_ptr<pair<double,double>> new_x1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> new_y1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> new_z1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> x1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> y1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> z1_bounds = make_shared<pair<double,double>>();
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
    
    
    var<> theta11("theta11",  std::max(-1.,r11._range->first), std::min(1.,r11._range->second)), theta12("theta12", std::max(-1.,r12._range->first), std::min(1.,r12._range->second)), theta13("theta13", std::max(-1.,r13._range->first), std::min(1.,r13._range->second));
    var<> theta21("theta21", std::max(-1.,r21._range->first), std::min(1.,r21._range->second)), theta22("theta22", std::max(-1.,r22._range->first), std::min(1.,r22._range->second)), theta23("theta23", std::max(-1.,r23._range->first), std::min(1.,r23._range->second));
    var<> theta31("theta31", std::max(-1.,r31._range->first), std::min(1.,r31._range->second)), theta32("theta32", std::max(-1.,r32._range->first), std::min(1.,r32._range->second)), theta33("theta33", std::max(-1.,r33._range->first), std::min(1.,r33._range->second));
    
    var<> x_shift("x_shift", shift_min_x, shift_max_x), y_shift("y_shift", shift_min_y, shift_max_y), z_shift("z_shift", shift_min_z, shift_max_z);
    indices valid_cells("valid_cells");
    map<int,vector<int>> valid_cells_map;
    map<int, vector<double>> dist_cost_map;
    int nd=point_cloud_data.size();
    int nm=point_cloud_model.size();
    vector<double> dist_i, coord_i;
    coord_i.resize(3);
    vector<vector<vector<double>>> box;
    vector<vector<double>> box_i;
    vector<double> zeros = {0,0,0};
    vector<shared_ptr<Model<double>>> batch_models;
    int missed=0;
    vector<double> x_lb, x_ub, y_lb, y_ub, z_lb, z_ub;
    for(auto i=0;i<nd;i++){
        auto max_dist_i = get_max_dist(roll_min, roll_max, pitch_min, pitch_max, yaw_min, yaw_max, shift_min_x, shift_max_x, shift_min_y, shift_max_y, shift_min_z, shift_max_z, point_cloud_data[i], zeros, false);
        dist_i.push_back(max_dist_i);
        x1_bounds->first = point_cloud_data.at(i)[0];
        x1_bounds->second = point_cloud_data.at(i)[0];
        y1_bounds->first = point_cloud_data.at(i)[1];
        y1_bounds->second = point_cloud_data.at(i)[1];
        z1_bounds->first = point_cloud_data.at(i)[2];
        z1_bounds->second = point_cloud_data.at(i)[2];
        auto x_range  = get_product_range(x1_bounds, theta11._range);
        auto y_range  = get_product_range(y1_bounds, theta12._range);
        auto z_range  = get_product_range(z1_bounds, theta13._range);
        *new_x1_bounds = {x_range->first + y_range->first + z_range->first + shift_min_x,
            x_range->second + y_range->second + z_range->second+ shift_max_x};
        x_lb.push_back(x_range->first + y_range->first + z_range->first + shift_min_x);
        x_ub.push_back(x_range->second + y_range->second + z_range->second+ shift_max_x);
        x_range  = get_product_range(x1_bounds, theta21._range);
        y_range  = get_product_range(y1_bounds, theta22._range);
        z_range  = get_product_range(z1_bounds, theta23._range);
        *new_y1_bounds = {x_range->first + y_range->first + z_range->first + shift_min_y,
            x_range->second + y_range->second + z_range->second+ shift_max_y};
        y_lb.push_back(x_range->first + y_range->first + z_range->first + shift_min_y);
        y_ub.push_back(x_range->second + y_range->second + z_range->second+ shift_max_y);
        x_range  = get_product_range(x1_bounds, theta31._range);
        y_range  = get_product_range(y1_bounds, theta32._range);
        z_range  = get_product_range(z1_bounds, theta33._range);
        *new_z1_bounds = {x_range->first + y_range->first + z_range->first + shift_min_z,
            x_range->second + y_range->second + z_range->second+ shift_max_z};
        z_lb.push_back(x_range->first + y_range->first + z_range->first + shift_min_z);
        z_ub.push_back(x_range->second + y_range->second + z_range->second+ shift_max_z);
        //        coord_i[0]=x_lb[i];
        //        coord_i[1]=y_lb[i];
        //        coord_i[2]=z_lb[i];
        //        box_i.push_back(coord_i);
        //        coord_i[0]=x_lb[i];
        //        coord_i[1]=y_ub[i];
        //        coord_i[2]=z_lb[i];
        //        box_i.push_back(coord_i);
        //        coord_i[0]=x_ub[i];
        //        coord_i[1]=y_ub[i];
        //        coord_i[2]=z_lb[i];
        //        box_i.push_back(coord_i);
        //        coord_i[0]=x_ub[i];
        //        coord_i[1]=y_lb[i];
        //        coord_i[2]=z_lb[i];
        //        box_i.push_back(coord_i);
        //        coord_i[0]=x_lb[i];
        //        coord_i[1]=y_lb[i];
        //        coord_i[2]=z_ub[i];
        //        box_i.push_back(coord_i);
        //        coord_i[0]=x_lb[i];
        //        coord_i[1]=y_ub[i];
        //        coord_i[2]=z_ub[i];
        //        box_i.push_back(coord_i);
        //        coord_i[0]=x_ub[i];
        //        coord_i[1]=y_ub[i];
        //        coord_i[2]=z_ub[i];
        //        box_i.push_back(coord_i);
        //        coord_i[0]=x_ub[i];
        //        coord_i[1]=y_lb[i];
        //        coord_i[2]=z_ub[i];
        //        box_i.push_back(coord_i);
        //        box.push_back(box_i);
        //        box_i.clear();
    }
    set<int> unique_model_pts;
    for (int j = 0; j<nm; j++) {
        auto voro_model=make_shared<Model<double>>("voro_model");
        var<> x("x"), y("y"),z("z");
        var<> lamda("lamda", 0, 1);
        var<> lamda1("lamda1", 0, 1);
        voro_model->add(x.in(R(1)), y.in(R(1)), z.in(R(1)));
        auto vertices=model_voronoi_vertices.at(j);
        int number_vert=vertices.size();
        voro_model->add(lamda.in(range(0,number_vert-1)));
        voro_model->add(lamda1.in(range(0,7)));
        param<> xi("xi"), yi("yi"), zi("zi");
        
        param<> xo("xo"), yo("yo"), zo("zo");
        param<> a("a"), b("b"), c("c"), d("d");
        param<> x_data("x_data"), y_data("y_data"), z_data("z_data");
        param<> x_model("x_model"), y_model("y_model"), z_model("z_model");
        param<> min_dist_sq("min_dist_sq");
        x_model.add_val("0", point_cloud_model.at(j)[0]);
        y_model.add_val("0", point_cloud_model.at(j)[1]);
        z_model.add_val("0", point_cloud_model.at(j)[2]);
        x_data.add_val("0", 0);
        y_data.add_val("0", 0);
        z_data.add_val("0", 0);
        min_dist_sq.add_val("0",0);
        auto normals=model_voronoi_normals.at(j);
        auto intercepts=model_face_intercept.at(j);
        indices voronoi_ids("voronoi_ids");
        voronoi_ids = range(1,intercepts.size());
        a.in(voronoi_ids);b.in(voronoi_ids);c.in(voronoi_ids);
        xi.in(range(0,number_vert-1));yi.in(range(0,number_vert-1));zi.in(range(0,number_vert-1));
        xo.in(range(0,7));yo.in(range(0,7));zo.in(range(0,7));
        for(auto k=0;k<intercepts.size();k++){
            a.set_val(k, normals[k][0]);
            b.set_val(k, normals[k][1]);
            c.set_val(k, normals[k][2]);
            d.add_val(k, intercepts[k]);
        }
        //        for(auto k=0;k<vertices.size();k++){
        //            xi.set_val(k, vertices[k][0]);
        //            yi.set_val(k, vertices[k][1]);
        //            zi.set_val(k, vertices[k][2]);
        //        }
        Constraint<> region("region");
        region=a*x+b*y+c*z+d;
        voro_model->add(region.in(voronoi_ids)<=0);
        
        //        Constraint<> def_x("def_x");
        //        def_x=x-sum(lamda*xi);
        //        voro_model->add(def_x==0);
        //
        //        Constraint<> def_y("def_y");
        //        def_y=y-sum(lamda*yi);
        //        voro_model->add(def_y==0);
        //
        //        Constraint<> def_z("def_z");
        //        def_z=z-sum(lamda*zi);
        //        voro_model->add(def_z==0);
        //
        //        Constraint<> def_lam("def_lam");
        //        def_lam=sum(lamda)-1;
        //        voro_model->add(def_lam==0);
        
        //        Constraint<> def_lam1("def_lam1");
        //        def_lam1=sum(lamda1)-1;
        //        voro_model->add(def_lam1==0);
        //
        //        Constraint<> def_x1("def_x1");
        //        def_x1=x-sum(lamda1*xo);
        //        voro_model->add(def_x1==0);
        //
        //        Constraint<> def_y1("def_y1");
        //        def_y1=y-sum(lamda1*yo);
        //        voro_model->add(def_y1==0);
        //
        //        Constraint<> def_z1("def_z1");
        //        def_z1=z-sum(lamda1*zo);
        //        voro_model->add(def_z1==0);
        
        
        //        Constraint<> radius("radius");
        //        radius=(x-x_data)*(x-x_data)+(y-y_data)*(y-y_data)+(z-z_data)*(z-z_data)-min_dist_sq;
        //        voro_model->add(radius.in(range(0,0))<=0);
        
        //        func<> obj=(x-x_model)*(x-x_model)+(y-y_model)*(y-y_model)+(z-z_model)*(z-z_model);
        //func<> obj=(x-x_data)*(x-x_data)+(y-y_data)*(y-y_data)+(z-z_data)*(z-z_data);
        func<> obj=0;
        voro_model->min(obj);
        for(auto i=0;i<nd;i++){
            string key_name=to_string(i)+","+to_string(j);
            if(!old_cells.has_key(to_string(i+1)+","+to_string(j+1)))
                continue;
            bool intersect=false;
            x_data = point_cloud_data[i][0];
            y_data = point_cloud_data[i][1];
            z_data = point_cloud_data[i][2];
            x.set_lb(x_lb[i]);
            x.set_ub(x_ub[i]);
            y.set_lb(y_lb[i]);
            y.set_ub(y_ub[i]);
            z.set_lb(z_lb[i]);
            z.set_ub(z_ub[i]);
            //            for(auto k=0;k<8;k++){
            //                xo.set_val(k, box[i][k][0]);
            //                yo.set_val(k, box[i][k][1]);
            //                zo.set_val(k, box[i][k][2]);
            //            }
            min_dist_sq=dist_i[i];
            auto modelc=voro_model->copy();
            modelc->set_name(key_name);
            //  modelc->print();
            batch_models.push_back(modelc);
        }
    }
    new_model_pts.clear();
    new_model_ids = indices();
    double time_end = get_wall_time();
    auto prep_time = time_end - time_start;
    DebugOn("Model build time = " << prep_time << endl);
    run_parallel(batch_models, ipopt, 1e-6, nb_total_threads, 1000);
    for(auto k=0;k<batch_models.size();k++){
        auto status_k=batch_models[k]->_status;
        auto obj_k=batch_models[k]->get_obj_val();
        auto keyk=batch_models[k]->get_name();
        auto ik=keyk.substr(0,keyk.find_first_of(","));
        auto jk=keyk.substr(keyk.find_first_of(",")+1);
        auto ik_int=std::atoi(ik.c_str());
        auto jk_int=std::atoi(jk.c_str());
        // if(status_k==0 && (sqrt(obj_k) - dist_i[ik_int])<=1e-6){
        if(status_k==0){
            auto d=std::max((obj_k-1e-6),0.0);
            DebugOff("Distance to voronoi cell = " << ik_int+1 << "," << jk_int+1 <<" "<< d << endl);
            if(unique_model_pts.insert(jk_int).second){
                new_model_pts.push_back(jk_int);
                new_model_ids.insert(to_string(jk_int+1));
            }
            valid_cells_map[ik_int].push_back(jk_int);
            dist_cost_map[ik_int].push_back(d);
            DebugOff("Valid pair: (" << ik_int+1 << "," << jk_int+1 << ")\n");
        }
        else{
            DebugOff("incompatible pair: (" << ik_int+1 << "," << jk_int+1 << ")\n");
        }
    }
    batch_models.clear();
    
    for (const auto &vcel:valid_cells_map) {
        auto key_data=to_string(vcel.first+1);
        auto cost_data=dist_cost_map[vcel.first];
        int count=0;
        for (auto const model_id: vcel.second) {
            auto key=key_data+","+to_string(model_id+1);
            valid_cells.insert(key);
            //            dist_cost.add_val(key, cost_data[count++]);
        }
    }
    DebugOn("Number of old pairs = " << old_cells.size() << endl);
    DebugOn("Number of valid cells = " << valid_cells.size() << endl);
    DebugOn("Number of discarded pairs = " << old_cells.size() - valid_cells.size() << endl);
    //    valid_cells.print();
    time_end = get_wall_time();
    prep_time = time_end - time_start;
    DebugOn("Voronoi preprocessing time = " << prep_time << endl);
    return valid_cells;
}
indices preprocess_new(const vector<vector<double>>& point_cloud_data, const vector<vector<double>>& point_cloud_model, const indices& old_cells, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double shift_min_x, double shift_max_x, double shift_min_y, double shift_max_y, double shift_min_z, double shift_max_z, const vector<vector<vector<double>>>& model_voronoi_normals, const vector<vector<double>>& model_face_intercept, const vector<vector<vector<double>>>& model_voronoi_vertices, vector<int>& new_model_pts, indices& new_model_ids, param<>& dist_cost, double upper_bound, int nb_total_threads){
    double time_start = get_wall_time();
    int new_test=0;
    shared_ptr<pair<double,double>> new_x1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> new_y1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> new_z1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> x1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> y1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> z1_bounds = make_shared<pair<double,double>>();
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
    
    
    var<> theta11("theta11",  std::max(-1.,r11._range->first), std::min(1.,r11._range->second)), theta12("theta12", std::max(-1.,r12._range->first), std::min(1.,r12._range->second)), theta13("theta13", std::max(-1.,r13._range->first), std::min(1.,r13._range->second));
    var<> theta21("theta21", std::max(-1.,r21._range->first), std::min(1.,r21._range->second)), theta22("theta22", std::max(-1.,r22._range->first), std::min(1.,r22._range->second)), theta23("theta23", std::max(-1.,r23._range->first), std::min(1.,r23._range->second));
    var<> theta31("theta31", std::max(-1.,r31._range->first), std::min(1.,r31._range->second)), theta32("theta32", std::max(-1.,r32._range->first), std::min(1.,r32._range->second)), theta33("theta33", std::max(-1.,r33._range->first), std::min(1.,r33._range->second));
    
    var<> x_shift("x_shift", shift_min_x, shift_max_x), y_shift("y_shift", shift_min_y, shift_max_y), z_shift("z_shift", shift_min_z, shift_max_z);
    indices valid_cells("valid_cells");
    map<int,vector<int>> valid_cells_map;
    map<int, vector<double>> dist_cost_map;
    int nd=point_cloud_data.size();
    int nm=point_cloud_model.size();
    vector<double> dist_i;
    vector<double> sphere_inner_sq, sphere_outer_sq;
    vector<double> zeros = {0,0,0};
    vector<shared_ptr<Model<double>>> batch_models;
    int missed=0;
    vector<double> x_lb, x_ub, y_lb, y_ub, z_lb, z_ub;
    double shift_mag_max=std::max(pow(shift_min_x,2),pow(shift_max_x,2))+std::max(pow(shift_min_y,2),pow(shift_max_y,2))+std::max(pow(shift_min_z,2),pow(shift_max_z,2));
    double shift_mag_max_root=sqrt(shift_mag_max);
    for(auto i=0;i<nd;i++){
        auto max_dist_i = get_max_dist(roll_min, roll_max, pitch_min, pitch_max, yaw_min, yaw_max, shift_min_x, shift_max_x, shift_min_y, shift_max_y, shift_min_z, shift_max_z, point_cloud_data[i], zeros, false);
        dist_i.push_back(max_dist_i);
        
        auto d_root=sqrt(pow(point_cloud_data.at(i)[0],2)+pow(point_cloud_data.at(i)[1],2)+pow(point_cloud_data.at(i)[2],2));
        sphere_inner_sq.push_back(pow((d_root-shift_mag_max_root),2));
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
        *new_x1_bounds = {x_range->first + y_range->first + z_range->first + shift_min_x,
            x_range->second + y_range->second + z_range->second+ shift_max_x};
        x_lb.push_back(x_range->first + y_range->first + z_range->first + shift_min_x);
        x_ub.push_back(x_range->second + y_range->second + z_range->second+ shift_max_x);
        x_range  = get_product_range(x1_bounds, theta21._range);
        y_range  = get_product_range(y1_bounds, theta22._range);
        z_range  = get_product_range(z1_bounds, theta23._range);
        *new_y1_bounds = {x_range->first + y_range->first + z_range->first + shift_min_y,
            x_range->second + y_range->second + z_range->second+ shift_max_y};
        y_lb.push_back(x_range->first + y_range->first + z_range->first + shift_min_y);
        y_ub.push_back(x_range->second + y_range->second + z_range->second+ shift_max_y);
        x_range  = get_product_range(x1_bounds, theta31._range);
        y_range  = get_product_range(y1_bounds, theta32._range);
        z_range  = get_product_range(z1_bounds, theta33._range);
        *new_z1_bounds = {x_range->first + y_range->first + z_range->first + shift_min_z,
            x_range->second + y_range->second + z_range->second+ shift_max_z};
        z_lb.push_back(x_range->first + y_range->first + z_range->first + shift_min_z);
        z_ub.push_back(x_range->second + y_range->second + z_range->second+ shift_max_z);
    }
    set<int> unique_model_pts;
    for (int j = 0; j<nm; j++) {
        auto voro_model=make_shared<Model<double>>("voro_model");
        var<> x("x",-1,1), y("y", -1,1),z("z", -1,1);
        voro_model->add(x.in(R(1)), y.in(R(1)), z.in(R(1)));
        param<> a("a"), b("b"), c("c"), d("d");
        param<> xlb("xlb"), xub("xub"),ylb("ylb"), yub("yub"),zlb("zlb"), zub("zub");
        param<> inner_radius_sq("inner_radius_sq");
        //        param<> x_data("x_data"), y_data("y_data"), z_data("z_data");
        //        param<> x_model("x_model"), y_model("y_model"), z_model("z_model");
        //        x_model.add_val("0", point_cloud_model.at(j)[0]);
        //        y_model.add_val("0", point_cloud_model.at(j)[1]);
        //        z_model.add_val("0", point_cloud_model.at(j)[2]);
        //        x_data.add_val("0", 0);
        //        y_data.add_val("0", 0);
        //        z_data.add_val("0", 0);
        xlb.add_val("0",0);
        xub.add_val("0",0);
        ylb.add_val("0",0);
        yub.add_val("0",0);
        zlb.add_val("0",0);
        zub.add_val("0",0);
        inner_radius_sq.add_val("0",0);
        auto vertices=model_voronoi_vertices.at(j);
        double dist_vj_max=-999;
        for(auto k=0;k<vertices.size();k++){
            auto v=vertices[k];
            auto dist_v=pow(v[0],2)+pow(v[1],2)+pow(v[2],2);
            if(dist_v>=dist_vj_max){
                dist_vj_max=dist_v;
            }
        }
        auto normals=model_voronoi_normals.at(j);
        auto intercepts=model_face_intercept.at(j);
        indices voronoi_ids("voronoi_ids");
        voronoi_ids = range(1,intercepts.size());
        a.in(voronoi_ids);b.in(voronoi_ids);c.in(voronoi_ids);
        for(auto k=0;k<intercepts.size();k++){
            a.set_val(k, normals[k][0]);
            b.set_val(k, normals[k][1]);
            c.set_val(k, normals[k][2]);
            d.add_val(k, intercepts[k]);
        }
        Constraint<> region("region");
        region=a*x+b*y+c*z+d;
        voro_model->add(region.in(voronoi_ids)<=0);
        
        Constraint<> inner_sphere("inner_sphere");
        inner_sphere=xlb*xub+ylb*yub+zlb*zub+inner_radius_sq-x*(xub+xlb)-y*(yub+ylb)-z*(zub+zlb);
        voro_model->add(inner_sphere.in(range(0,0))<=0);
        
        //        Constraint<> radius("radius");
        //        radius=(x-x_data)*(x-x_data)+(y-y_data)*(y-y_data)+(z-z_data)*(z-z_data)-min_dist_sq;
        //        voro_model->add(radius.in(range(0,0))<=0);
        
        //        func<> obj=(x-x_model)*(x-x_model)+(y-y_model)*(y-y_model)+(z-z_model)*(z-z_model);
        func<> obj=(x)*(x)+(y)*(y)+(z)*(z);
        voro_model->min(obj);
        for(auto i=0;i<nd;i++){
            string key_name=to_string(i)+","+to_string(j);
            if(!old_cells.has_key(to_string(i+1)+","+to_string(j+1)))
                continue;
            if(dist_vj_max<=sphere_inner_sq[i]-1e-6){
                new_test++;
                continue;
            }
            
            bool intersect=false;
            //            x_data = point_cloud_data[i][0];
            //            y_data = point_cloud_data[i][1];
            //            z_data = point_cloud_data[i][2];
            inner_radius_sq=sphere_inner_sq[i];
            xlb=x_lb[i];
            xub=x_ub[i];
            ylb=y_lb[i];
            yub=y_ub[i];
            zlb=z_lb[i];
            zub=z_ub[i];
            x.set_lb(x_lb[i]);
            x.set_ub(x_ub[i]);
            y.set_lb(y_lb[i]);
            y.set_ub(y_ub[i]);
            z.set_lb(z_lb[i]);
            z.set_ub(z_ub[i]);
            
            auto modelc=voro_model->copy();
            modelc->set_name(key_name);
            batch_models.push_back(modelc);
        }
    }
    
    new_model_pts.clear();
    new_model_ids = indices();
    double time_end = get_wall_time();
    auto prep_time = time_end - time_start;
    DebugOn("Model build time = " << prep_time << endl);
    run_parallel(batch_models, ipopt, 1e-6, nb_total_threads, 1000);
    for(auto k=0;k<batch_models.size();k++){
        auto status_k=batch_models[k]->_status;
        auto obj_k=batch_models[k]->get_obj_val();
        auto keyk=batch_models[k]->get_name();
        auto ik=keyk.substr(0,keyk.find_first_of(","));
        auto jk=keyk.substr(keyk.find_first_of(",")+1);
        auto ik_int=std::atoi(ik.c_str());
        auto jk_int=std::atoi(jk.c_str());
        if(status_k==0 && ((obj_k - sphere_outer_sq[ik_int])<=1e-6)){
            auto d=std::max((obj_k-1e-6),0.0);
            DebugOff("Distance to voronoi cell = " << ik_int+1 << "," << jk_int+1 <<" "<< d << endl);
            if(unique_model_pts.insert(jk_int).second){
                new_model_pts.push_back(jk_int);
                new_model_ids.insert(to_string(jk_int+1));
            }
            valid_cells_map[ik_int].push_back(jk_int);
            dist_cost_map[ik_int].push_back(d);
            DebugOff("Valid pair: (" << ik_int+1 << "," << jk_int+1 << ")\n");
        }
        else{
            DebugOff("incompatible pair: (" << ik_int+1 << "," << jk_int+1 << ")\n");
        }
    }
    batch_models.clear();
    
    for (const auto &vcel:valid_cells_map) {
        auto key_data=to_string(vcel.first+1);
        auto cost_data=dist_cost_map[vcel.first];
        int count=0;
        for (auto const model_id: vcel.second) {
            auto key=key_data+","+to_string(model_id+1);
            valid_cells.insert(key);
            //            dist_cost.add_val(key, cost_data[count++]);
        }
    }
    DebugOn("Number of old pairs = " << old_cells.size() << endl);
    DebugOn("Number of valid cells = " << valid_cells.size() << endl);
    DebugOn("Number of discarded pairs = " << old_cells.size() - valid_cells.size() << endl);
    DebugOn("Number of discarded pairs by inner sphere test = " << new_test << endl);
    //    valid_cells.print();
    time_end = get_wall_time();
    prep_time = time_end - time_start;
    DebugOn("Voronoi preprocessing time = " << prep_time << endl);
    return valid_cells;
}
#define INTER_MAX_ITER 100000
unsigned long planeStatPerPair;
unsigned long nbFails;
static SphericalPolygon tupTemp_;
struct ThreeUniquePoints {
    vector<Vec3f> pts_;
    int pos_ = 0;
    ThreeUniquePoints() : pts_() {}
    void add(const Vec3f & v) {
        bool same_found(false);
        for( const auto & p : pts_ ) {
            if( v == p ) {
                same_found = true;
                break;
            }
        }
        if( same_found ) return;
#define TUPSIZE 3
        if( pts_.size() < TUPSIZE ) {
            pts_.push_back(v);
            return;
        } else {
            pts_[pos_] = v;
            pos_ = (pos_ + 1) % TUPSIZE;
        }
    }
    void clip(const SphericalPolygon & src, SphericalPolygon & dest, const Vec3f & delta) {
        int i(1);
        src.clip(pts_[0] - delta, dest);
        for( ; i < pts_.size(); ) {
            dest.clip(pts_[i++] - delta, tupTemp_);
            dest.swap(tupTemp_);
        }
    }
};

//  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -

template< typename Convex >
bool sphericalDisjoint(const Convex & a, const Convex & b) {
    static SphericalPolygon positiveBound, tempPoly;
    int vA, vB;
    Vec3f dir(b.vertex(0) - a.vertex(0));
    if( ! a.differenceCoversZeroInDir(b, vA, vB, dir) ) return true;
    positiveBound.clear();
    positiveBound.emplace_back(dir);
    positiveBound.clip(b.vertex(vB) - a.vertex(vA), tempPoly); positiveBound.swap(tempPoly);
    if( positiveBound.empty() ) return false;
    do {
        if( ! a.differenceCoversZeroInDir(b, vA, vB, positiveBound.averageDirection()) ) return true;
        if(planeStatPerPair >= INTER_MAX_ITER) {
            ++nbFails;
            DebugOn("fail "<<endl);
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
    auto res=sphericalDisjoint(convexes[0], convexes[1]);
    return res;
    
}
indices preprocess_poltyope_intersect(const vector<vector<double>>& point_cloud_data, const vector<vector<double>>& point_cloud_model, const indices& old_cells, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double shift_min_x, double shift_max_x, double shift_min_y, double shift_max_y, double shift_min_z, double shift_max_z, const vector<vector<vector<double>>>& model_voronoi_normals, const vector<vector<double>>& model_face_intercept, const vector<vector<vector<double>>>& model_voronoi_vertices, vector<int>& new_model_pts, indices& new_model_ids, param<>& dist_cost, double upper_bound, int nb_total_threads){
    nbFails=0;
    planeStatPerPair = 0;
    indices valid_cells("valid_cells");
    
    double time_start = get_wall_time();
    int new_test=0;
    shared_ptr<pair<double,double>> new_x1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> new_y1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> new_z1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> x1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> y1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> z1_bounds = make_shared<pair<double,double>>();
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
    
    
    var<> theta11("theta11",  std::max(-1.,r11._range->first), std::min(1.,r11._range->second)), theta12("theta12", std::max(-1.,r12._range->first), std::min(1.,r12._range->second)), theta13("theta13", std::max(-1.,r13._range->first), std::min(1.,r13._range->second));
    var<> theta21("theta21", std::max(-1.,r21._range->first), std::min(1.,r21._range->second)), theta22("theta22", std::max(-1.,r22._range->first), std::min(1.,r22._range->second)), theta23("theta23", std::max(-1.,r23._range->first), std::min(1.,r23._range->second));
    var<> theta31("theta31", std::max(-1.,r31._range->first), std::min(1.,r31._range->second)), theta32("theta32", std::max(-1.,r32._range->first), std::min(1.,r32._range->second)), theta33("theta33", std::max(-1.,r33._range->first), std::min(1.,r33._range->second));
    
    
    func<> row1 = pow(theta11,2)+pow(theta12,2)+pow(theta13,2);
    func<> row2 = pow(theta21,2)+pow(theta22,2)+pow(theta23,2);
    func<> row3 = pow(theta31,2)+pow(theta32,2)+pow(theta33,2);
    func<> col1 = pow(theta11,2)+pow(theta21,2)+pow(theta31,2);
    func<> col2 = pow(theta12,2)+pow(theta22,2)+pow(theta32,2);
    func<> col3 = pow(theta13,2)+pow(theta23,2)+pow(theta33,2);
    func<> row12 = theta11*theta21+theta12*theta22+theta13*theta23;
    func<> row13 = theta11*theta31+theta12*theta32+theta13*theta33;
    func<> row23 = theta21*theta31+theta22*theta32+theta23*theta33;
    func<> det= theta11*(theta22*theta33-theta32*theta23)
    -theta12*(theta21*theta33-theta31*theta23)+theta13*(theta21*theta32-theta31*theta22);
    
    bool feas=true;
    
    if(row1._range->first>=(1+1e-4) || row1._range->second<=(1-1e-4)){
        DebugOn("row1 "<<row1._range->first<<" "<<row1._range->second<<endl);
        feas=false;
    }
    else if(row2._range->first>=(1+1e-4) || row2._range->second<=(1-1e-4)){
        DebugOn("row2 "<<row2._range->first<<" "<<row2._range->second<<endl);
        feas=false;
    }
    else if(row3._range->first>=(1+1e-4) || row3._range->second<=(1-1e-4)){
        DebugOn("row3 "<<row3._range->first<<" "<<row3._range->second<<endl);
        feas=false;
    }
    else if(col1._range->first>=(1+1e-4) || col1._range->second<=(1-1e-4)){
        DebugOn("col1 "<<col1._range->first<<" "<<col1._range->second<<endl);
        feas=false;
    }
    else if(col2._range->first>=(1+1e-4) || col2._range->second<=(1-1e-4)){
        DebugOn("col2 "<<col2._range->first<<" "<<col2._range->second<<endl);
        feas=false;
    }
    else if(col3._range->first>=(1+1e-4) || col3._range->second<=(1-1e-4)){
        DebugOn("col3 "<<col3._range->first<<" "<<col3._range->second<<endl);
        feas=false;
    }
    else if(row12._range->first>=1e-4 || row12._range->second<=(-1e-4)){
        DebugOn("row12 "<<row12._range->first<<" "<<row12._range->second<<endl);
        feas=false;
    }
    else if(row13._range->first>=1e-4 || row13._range->second<=(-1e-4)){
        DebugOn("row13 "<<row13._range->first<<" "<<row13._range->second<<endl);
        feas=false;
    }
    else if(row23._range->first>=1e-4 || row23._range->second<=(-1e-4)){
        DebugOn("row23 "<<row23._range->first<<" "<<row23._range->second<<endl);
        feas=false;
    }
    else if(det._range->first>=(1+1e-4) || det._range->second<=(1-1e-4)){
        DebugOn("det "<<det._range->first<<" "<<det._range->second<<endl);
        feas=false;
    }

    if(feas){
        
        var<> x_shift("x_shift", shift_min_x, shift_max_x), y_shift("y_shift", shift_min_y, shift_max_y), z_shift("z_shift", shift_min_z, shift_max_z);
        
        map<int,vector<int>> valid_cells_map;
        map<int, vector<double>> dist_cost_map;
        int nd=point_cloud_data.size();
        int nm=point_cloud_model.size();
        vector<double> dist_i;
        vector<double> sphere_inner_sq, sphere_outer_sq;
        vector<double> zeros = {0,0,0};
        vector<shared_ptr<Model<double>>> batch_models;
        int missed=0;
        vector<double> x_lb, x_ub, y_lb, y_ub, z_lb, z_ub;
        vector<vector<vector<double>>> box;
        vector<vector<double>> box_i;
        vector<double> coord_i;
        coord_i.resize(3);
        double shift_mag_max=std::max(pow(shift_min_x,2),pow(shift_max_x,2))+std::max(pow(shift_min_y,2),pow(shift_max_y,2))+std::max(pow(shift_min_z,2),pow(shift_max_z,2));
        double shift_mag_max_root=sqrt(shift_mag_max);
        for(auto i=0;i<nd;i++){
            //  auto max_dist_i = get_max_dist(roll_min, roll_max, pitch_min, pitch_max, yaw_min, yaw_max, shift_min_x, shift_max_x, shift_min_y, shift_max_y, shift_min_z, shift_max_z, point_cloud_data[i], zeros, false);
            //dist_i.push_back(max_dist_i);
            
            auto d_root=sqrt(pow(point_cloud_data.at(i)[0],2)+pow(point_cloud_data.at(i)[1],2)+pow(point_cloud_data.at(i)[2],2));
            sphere_inner_sq.push_back(pow((d_root-shift_mag_max_root),2));
            // sphere_outer_sq.push_back(pow((d_root+shift_mag_max_root),2));
            
            x1_bounds->first = point_cloud_data.at(i)[0];
            x1_bounds->second = point_cloud_data.at(i)[0];
            y1_bounds->first = point_cloud_data.at(i)[1];
            y1_bounds->second = point_cloud_data.at(i)[1];
            z1_bounds->first = point_cloud_data.at(i)[2];
            z1_bounds->second = point_cloud_data.at(i)[2];
            auto x_range  = get_product_range(x1_bounds, theta11._range);
            auto y_range  = get_product_range(y1_bounds, theta12._range);
            auto z_range  = get_product_range(z1_bounds, theta13._range);
            *new_x1_bounds = {x_range->first + y_range->first + z_range->first + shift_min_x,
                x_range->second + y_range->second + z_range->second+ shift_max_x};
            x_lb.push_back(x_range->first + y_range->first + z_range->first + shift_min_x);
            x_ub.push_back(x_range->second + y_range->second + z_range->second+ shift_max_x);
            x_range  = get_product_range(x1_bounds, theta21._range);
            y_range  = get_product_range(y1_bounds, theta22._range);
            z_range  = get_product_range(z1_bounds, theta23._range);
            *new_y1_bounds = {x_range->first + y_range->first + z_range->first + shift_min_y,
                x_range->second + y_range->second + z_range->second+ shift_max_y};
            y_lb.push_back(x_range->first + y_range->first + z_range->first + shift_min_y);
            y_ub.push_back(x_range->second + y_range->second + z_range->second+ shift_max_y);
            x_range  = get_product_range(x1_bounds, theta31._range);
            y_range  = get_product_range(y1_bounds, theta32._range);
            z_range  = get_product_range(z1_bounds, theta33._range);
            *new_z1_bounds = {x_range->first + y_range->first + z_range->first + shift_min_z,
                x_range->second + y_range->second + z_range->second+ shift_max_z};
            z_lb.push_back(x_range->first + y_range->first + z_range->first + shift_min_z);
            z_ub.push_back(x_range->second + y_range->second + z_range->second+ shift_max_z);
            coord_i[0]=x_lb[i];
            coord_i[1]=y_lb[i];
            coord_i[2]=z_lb[i];
            box_i.push_back(coord_i);
            coord_i[0]=x_lb[i];
            coord_i[1]=y_ub[i];
            coord_i[2]=z_lb[i];
            box_i.push_back(coord_i);
            coord_i[0]=x_ub[i];
            coord_i[1]=y_ub[i];
            coord_i[2]=z_lb[i];
            box_i.push_back(coord_i);
            coord_i[0]=x_ub[i];
            coord_i[1]=y_lb[i];
            coord_i[2]=z_lb[i];
            box_i.push_back(coord_i);
            coord_i[0]=x_lb[i];
            coord_i[1]=y_lb[i];
            coord_i[2]=z_ub[i];
            box_i.push_back(coord_i);
            coord_i[0]=x_lb[i];
            coord_i[1]=y_ub[i];
            coord_i[2]=z_ub[i];
            box_i.push_back(coord_i);
            coord_i[0]=x_ub[i];
            coord_i[1]=y_ub[i];
            coord_i[2]=z_ub[i];
            box_i.push_back(coord_i);
            coord_i[0]=x_ub[i];
            coord_i[1]=y_lb[i];
            coord_i[2]=z_ub[i];
            box_i.push_back(coord_i);
            box.push_back(box_i);
            box_i.clear();
        }
        set<int> unique_model_pts;
        for (int j = 0; j<nm; j++) {
            auto vertices=model_voronoi_vertices.at(j);
            for(auto i=0;i<nd;i++){
                if(!old_cells.has_key(to_string(i+1)+","+to_string(j+1))){
                    DebugOff("continued");
                    continue;
                }
                auto res=test_general<VerticesOnly>(box[i],vertices);
                DebugOff("i "<<i <<" j "<<j<<" "<<res<<endl);
                if(!res){
                    auto dist_vj_max=-999.0;
                    for(auto k=0;k<vertices.size();k++){
                        auto v=vertices[k];
                        auto dist_v=pow(v[0],2)+pow(v[1],2)+pow(v[2],2);
                        if(dist_v>=dist_vj_max){
                            dist_vj_max=dist_v;
                        }
                    }
                    if(dist_vj_max>=sphere_inner_sq[i]-1e-6){
                        if(unique_model_pts.insert(j).second){
                            new_model_pts.push_back(j);
                            new_model_ids.insert(to_string(i));
                        }
                        valid_cells_map[i].push_back(j);
                    }
                }
            }
        }
        for(auto i=0;i<nd;i++){
            DebugOff(valid_cells_map[i].size()<<endl);
        }
        bool found_all=true;
        if(valid_cells_map.size()<nd){
            found_all=false;
        }
        if(found_all){
            for (const auto &vcel:valid_cells_map) {
                auto key_data=to_string(vcel.first+1);
                //  auto cost_data=dist_cost_map[vcel.first];
                int count=0;
                for (auto const model_id: vcel.second) {
                    auto key=key_data+","+to_string(model_id+1);
                    valid_cells.insert(key);
                    count++;
                    //            dist_cost.add_val(key, cost_data[count++]);
                }
                if(count==0){
                    DebugOn("no possible match found for data point "<<key_data<<endl);
                    found_all=false;
                    break;
                }
            }
        }
        if(found_all){
            DebugOff("Number of old pairs = " << old_cells.size() << endl);
            DebugOff("Number of valid cells = " << valid_cells.size() << endl);
            DebugOff("Number of discarded pairs = " << old_cells.size() - valid_cells.size() << endl);
            
        }
        else{
            
            DebugOff("Number of valid cells = " << 0 << endl);
        }
    }
    //  DebugOn("Number of discarded pairs by inner sphere test = " << new_test << endl);
    //    valid_cells.print();
    //        time_end = get_wall_time();
    //        prep_time = time_end - time_start;
    //        DebugOn("Voronoi preprocessing time = " << prep_time << endl);
    auto time_end = get_wall_time();
    auto prep_time = time_end - time_start;
    DebugOff("Voronoi preprocessing time = " << prep_time << endl);
    return valid_cells;
}
//void bounds_new_x(const vector<vector<double>>& point_cloud_data, const vector<vector<double>>& point_cloud_model, const indices& old_cells, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double shift_min_x, double shift_max_x, double shift_min_y, double shift_max_y, double shift_min_z, double shift_max_z, const vector<vector<vector<double>>>& model_voronoi_normals, const vector<vector<double>>& model_face_intercept, const vector<vector<vector<double>>>& model_voronoi_vertices, vector<int>& new_model_pts, indices& new_model_ids, param<>& dist_cost, double upper_bound, int nb_total_threads){
//    nbFails=0;
//    planeStatPerPair = 0;
//
//
//        double time_start = get_wall_time();
//        int new_test=0;
//        shared_ptr<pair<double,double>> new_x1_bounds = make_shared<pair<double,double>>();
//        shared_ptr<pair<double,double>> new_y1_bounds = make_shared<pair<double,double>>();
//        shared_ptr<pair<double,double>> new_z1_bounds = make_shared<pair<double,double>>();
//        shared_ptr<pair<double,double>> x1_bounds = make_shared<pair<double,double>>();
//        shared_ptr<pair<double,double>> y1_bounds = make_shared<pair<double,double>>();
//        shared_ptr<pair<double,double>> z1_bounds = make_shared<pair<double,double>>();
//        var<> yaw("yaw", yaw_min, yaw_max), pitch("pitch", pitch_min, pitch_max), roll("roll", roll_min, roll_max);
//        yaw.in(R(1)); pitch.in(R(1));roll.in(R(1));
//        func<> r11 = cos(yaw)*cos(roll);
//        func<> r12 = cos(yaw)*sin(roll)*sin(pitch) - sin(yaw)*cos(pitch);
//        func<> r13 = cos(yaw)*sin(roll)*cos(pitch) + sin(yaw)*sin(pitch);
//        func<> r21 = sin(yaw)*cos(roll);
//        func<> r22 = sin(yaw)*sin(roll)*sin(pitch) + cos(yaw)*cos(pitch);
//        func<> r23 = sin(yaw)*sin(roll)*cos(pitch) - cos(yaw)*sin(pitch);
//        func<> r31 = sin(-1*roll);
//        func<> r32 = cos(roll)*sin(pitch);
//        func<> r33 = cos(roll)*cos(pitch);
//
//
//        var<> theta11("theta11",  std::max(-1.,r11._range->first), std::min(1.,r11._range->second)), theta12("theta12", std::max(-1.,r12._range->first), std::min(1.,r12._range->second)), theta13("theta13", std::max(-1.,r13._range->first), std::min(1.,r13._range->second));
//        var<> theta21("theta21", std::max(-1.,r21._range->first), std::min(1.,r21._range->second)), theta22("theta22", std::max(-1.,r22._range->first), std::min(1.,r22._range->second)), theta23("theta23", std::max(-1.,r23._range->first), std::min(1.,r23._range->second));
//        var<> theta31("theta31", std::max(-1.,r31._range->first), std::min(1.,r31._range->second)), theta32("theta32", std::max(-1.,r32._range->first), std::min(1.,r32._range->second)), theta33("theta33", std::max(-1.,r33._range->first), std::min(1.,r33._range->second));
//
//        var<> x_shift("x_shift", shift_min_x, shift_max_x), y_shift("y_shift", shift_min_y, shift_max_y), z_shift("z_shift", shift_min_z, shift_max_z);
//        indices valid_cells("valid_cells");
//        map<int,vector<int>> valid_cells_map;
//        map<int, vector<double>> dist_cost_map;
//        int nd=point_cloud_data.size();
//        int nm=point_cloud_model.size();
//        vector<double> dist_i;
//        vector<double> sphere_inner_sq, sphere_outer_sq;
//        vector<double> zeros = {0,0,0};
//        vector<shared_ptr<Model<double>>> batch_models;
//        int missed=0;
//        vector<double> x_lb, x_ub, y_lb, y_ub, z_lb, z_ub;
//        vector<vector<vector<double>>> box;
//        vector<vector<double>> box_i;
//        vector<double> coord_i;
//        coord_i.resize(3);
//        double shift_mag_max=std::max(pow(shift_min_x,2),pow(shift_max_x,2))+std::max(pow(shift_min_y,2),pow(shift_max_y,2))+std::max(pow(shift_min_z,2),pow(shift_max_z,2));
//        double shift_mag_max_root=sqrt(shift_mag_max);
//        for(auto i=0;i<nd;i++){
//          //  auto max_dist_i = get_max_dist(roll_min, roll_max, pitch_min, pitch_max, yaw_min, yaw_max, shift_min_x, shift_max_x, shift_min_y, shift_max_y, shift_min_z, shift_max_z, point_cloud_data[i], zeros, false);
//            //dist_i.push_back(max_dist_i);
//
//         auto d_root=sqrt(pow(point_cloud_data.at(i)[0],2)+pow(point_cloud_data.at(i)[1],2)+pow(point_cloud_data.at(i)[2],2));
//            sphere_inner_sq.push_back(pow((d_root-shift_mag_max_root),2));
//            sphere_outer_sq.push_back(pow((d_root+shift_mag_max_root),2));
//
//            x1_bounds->first = point_cloud_data.at(i)[0];
//            x1_bounds->second = point_cloud_data.at(i)[0];
//            y1_bounds->first = point_cloud_data.at(i)[1];
//            y1_bounds->second = point_cloud_data.at(i)[1];
//            z1_bounds->first = point_cloud_data.at(i)[2];
//            z1_bounds->second = point_cloud_data.at(i)[2];
//            auto x_range  = get_product_range(x1_bounds, theta11._range);
//            auto y_range  = get_product_range(y1_bounds, theta12._range);
//            auto z_range  = get_product_range(z1_bounds, theta13._range);
//            *new_x1_bounds = {x_range->first + y_range->first + z_range->first + shift_min_x,
//                x_range->second + y_range->second + z_range->second+ shift_max_x};
//            x_lb.push_back(x_range->first + y_range->first + z_range->first + shift_min_x);
//            x_ub.push_back(x_range->second + y_range->second + z_range->second+ shift_max_x);
//            x_range  = get_product_range(x1_bounds, theta21._range);
//            y_range  = get_product_range(y1_bounds, theta22._range);
//            z_range  = get_product_range(z1_bounds, theta23._range);
//            *new_y1_bounds = {x_range->first + y_range->first + z_range->first + shift_min_y,
//                x_range->second + y_range->second + z_range->second+ shift_max_y};
//            y_lb.push_back(x_range->first + y_range->first + z_range->first + shift_min_y);
//            y_ub.push_back(x_range->second + y_range->second + z_range->second+ shift_max_y);
//            x_range  = get_product_range(x1_bounds, theta31._range);
//            y_range  = get_product_range(y1_bounds, theta32._range);
//            z_range  = get_product_range(z1_bounds, theta33._range);
//            *new_z1_bounds = {x_range->first + y_range->first + z_range->first + shift_min_z,
//                x_range->second + y_range->second + z_range->second+ shift_max_z};
//            z_lb.push_back(x_range->first + y_range->first + z_range->first + shift_min_z);
//            z_ub.push_back(x_range->second + y_range->second + z_range->second+ shift_max_z);
//            coord_i[0]=x_lb[i];
//            coord_i[1]=y_lb[i];
//            coord_i[2]=z_lb[i];
//            box_i.push_back(coord_i);
//            coord_i[0]=x_lb[i];
//            coord_i[1]=y_ub[i];
//            coord_i[2]=z_lb[i];
//            box_i.push_back(coord_i);
//            coord_i[0]=x_ub[i];
//            coord_i[1]=y_ub[i];
//            coord_i[2]=z_lb[i];
//            box_i.push_back(coord_i);
//            coord_i[0]=x_ub[i];
//            coord_i[1]=y_lb[i];
//            coord_i[2]=z_lb[i];
//            box_i.push_back(coord_i);
//            coord_i[0]=x_lb[i];
//            coord_i[1]=y_lb[i];
//            coord_i[2]=z_ub[i];
//            box_i.push_back(coord_i);
//            coord_i[0]=x_lb[i];
//            coord_i[1]=y_ub[i];
//            coord_i[2]=z_ub[i];
//            box_i.push_back(coord_i);
//            coord_i[0]=x_ub[i];
//            coord_i[1]=y_ub[i];
//            coord_i[2]=z_ub[i];
//            box_i.push_back(coord_i);
//            coord_i[0]=x_ub[i];
//            coord_i[1]=y_lb[i];
//            coord_i[2]=z_ub[i];
//            box_i.push_back(coord_i);
//            box.push_back(box_i);
//            box_i.clear();
//        }
//        set<int> unique_model_pts;
//
//            for(auto i=0;i<nd;i++){
//                auto voro_model=make_shared<Model<double>>("voro_model");
//                int jw=10;
//                indices facets("facets");
//                facets=range(0, jw);
//                var<> xi("xi"), yi("yi"), zi("zi");
//                voro_model->add(xi.in(facets), yi.in(facets),zi.in(facets));
//                var<> lamda("lamda");
//                voro_model->add(lamda.in(facets));
//
//                var<> x("x", x_lb[i], x_ub[i]), y("y", y_lb[i], y_ub[i]), z("z",z_lb[i], z_ub[i]);
//                voro_model->add(x.in(R(1)), y.in(R(1)),z.in(R(1)));
//                param<> a("a"), b("b"), c("c"), d("d");
//
//
//                Constraint<> inner_sphere("inner_sphere");
//                inner_sphere=x_lb[i]*x_ub[i]+y_lb[i]*y_ub[i]+z_lb[i]*z_ub[i]+sphere_inner_sq[i]-x*(x_ub[i]+x_lb[i])-y*(y_ub[i]+y_lb[i])-z*(z_ub[i]+z_lb[i]);
//                voro_model->add(inner_sphere.in(range(0,0))<=0);
//                for (int j = 0; j<jw; j++) {
//                    auto vertices=model_voronoi_vertices.at(j);
//                if(!old_cells.has_key(to_string(i+1)+","+to_string(j+1))){
//                    DebugOff("continued");
//                    continue;
//                }
//                    auto normals=model_voronoi_normals.at(j);
//                    auto intercepts=model_face_intercept.at(j);
//                    indices voronoi_ids("voronoi_ids");
//                    voronoi_ids = range(1,intercepts.size());
//                    a.in(voronoi_ids);b.in(voronoi_ids);c.in(voronoi_ids);
//                    for(auto k=0;k<intercepts.size();k++){
//                        a.set_val(k, normals[k][0]);
//                        b.set_val(k, normals[k][1]);
//                        c.set_val(k, normals[k][2]);
//                        d.add_val(k, intercepts[k]);
//                    }
//                    Constraint<> region("region"+to_string(j));
//                    region=a*xi(j)+b*yi(j)+c*zi(j)+d*lamda(j);
//                    voro_model->add(region.in(voronoi_ids)<=0);
//                }
//
//
////                DebugOff("i "<<i <<" j "<<j<<" "<<res<<endl);
////                if(!res){
////                if(unique_model_pts.insert(j).second){
////                    new_model_pts.push_back(j);
////                    new_model_ids.insert(to_string(i));
////                }
////                valid_cells_map[i].push_back(j);
////                }
//
//        }
//    for(auto i=0;i<nd;i++){
//        DebugOff(valid_cells_map[i].size()<<endl);
//    }
//    bool found_all=true;
//    if(valid_cells_map.size()<nd){
//        found_all=false;
//    }
//    if(found_all){
//        for (const auto &vcel:valid_cells_map) {
//            auto key_data=to_string(vcel.first+1);
//          //  auto cost_data=dist_cost_map[vcel.first];
//            int count=0;
//            for (auto const model_id: vcel.second) {
//                auto key=key_data+","+to_string(model_id+1);
//                valid_cells.insert(key);
//                count++;
//    //            dist_cost.add_val(key, cost_data[count++]);
//            }
//            if(count==0){
//                DebugOn("no possible match found for data point "<<key_data<<endl);
//                found_all=false;
//                break;
//            }
//        }
//    }
//    if(found_all){
//        DebugOn("Number of old pairs = " << old_cells.size() << endl);
//        DebugOn("Number of valid cells = " << valid_cells.size() << endl);
//        DebugOn("Number of discarded pairs = " << old_cells.size() - valid_cells.size() << endl);
//
//    }
//    else{
//
//        DebugOn("Number of valid cells = " << 0 << endl);
//    }
//      //  DebugOn("Number of discarded pairs by inner sphere test = " << new_test << endl);
//    //    valid_cells.print();
////        time_end = get_wall_time();
////        prep_time = time_end - time_start;
////        DebugOn("Voronoi preprocessing time = " << prep_time << endl);
//    auto time_end = get_wall_time();
//    auto prep_time = time_end - time_start;
//    DebugOn("Voronoi preprocessing time = " << prep_time << endl);
//        return valid_cells;
//    }


vector<double> BranchBound(vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, param<>& norm_x, param<>& norm_y, param<>& norm_z,  param<>& intercept, const vector<int>& init_matching, const vector<double>& init_err_per_point, param<>& model_radius, vector<vector<vector<double>>> model_voronoi_normals, vector<vector<double>> model_face_intercept, vector<vector<vector<double>>> model_voronoi_vertices, vector<int>& new_model_pts, indices& new_model_ids, param<>& dist_cost, bool relax_ints, bool relax_sdp, bool rigid_transf) {
    /* INPUT BOUNDS */
    double time_start = get_wall_time();
    double total_time_max = 9000;
    double shift_min_x =  -0.1, shift_max_x = 0.1, shift_min_y = -0.1,shift_max_y = 0.1,shift_min_z = -0.1,shift_max_z = 0.1;
    double yaw_min = -10*pi/180., yaw_max = 10*pi/180., pitch_min =-10*pi/180.,pitch_max = 10*pi/180.,roll_min =-10*pi/180.,roll_max = 10*pi/180.;
    indices N1 = range(1,point_cloud_data.size());
    indices N2 = range(1,point_cloud_model.size());
    vector<int> new_matching(N1.size());
    bool convex = false;
    double max_time = 30;
    int max_iter = 1e6;
    vector<pair<pair<int,int>,pair<int,int>>> incompatible_pairs;
    auto nb_threads = std::thread::hardware_concurrency()/2;
    nb_threads = 1;
    vector<pair<double,double>> shift_x_bounds(nb_threads), shift_y_bounds(nb_threads), shift_z_bounds(nb_threads);
    vector<pair<double,double>> roll_bounds(nb_threads), pitch_bounds(nb_threads), yaw_bounds(nb_threads);
    vector<indices> valid_cells(nb_threads, indices("valid_cells"));
    vector<vector<double>> rot_trans(nb_threads, vector<double>(12));
    vector<double> best_rot_trans(12);
    DebugOn("I will be using " << nb_threads << " parallel threads" << endl);
    vector<shared_ptr<Model<>>> models;
    double lb = 0, ub = 0, best_lb = 0, best_ub = numeric_limits<double>::max();
    int nb_pruned = 0;
    priority_queue<treenode> lb_queue;
    for (int i = 0; i<nb_threads; i++) {
        valid_cells[i] = indices(N1,N2);
        double x_shift_increment = (shift_max_x - shift_min_x)/nb_threads;
        lb = shift_min_x + i*x_shift_increment;
        ub = shift_min_x+(i+1)*x_shift_increment;
        shift_x_bounds[i] = {lb, ub};
        shift_y_bounds[i] = {shift_min_y, shift_max_y};
        shift_z_bounds[i] = {shift_min_z, shift_max_z};
        roll_bounds[i] = {roll_min, roll_max};
        pitch_bounds[i] = {pitch_min, pitch_max};
        yaw_bounds[i] = {yaw_min, yaw_max};
        valid_cells[i] = preprocess_poltyope_intersect(point_cloud_data, point_cloud_model, valid_cells[i], roll_bounds[i].first, roll_bounds[i].second,  pitch_bounds[i].first, pitch_bounds[i].second, yaw_bounds[i].first, yaw_bounds[i].second, shift_x_bounds[i].first, shift_x_bounds[i].second, shift_y_bounds[i].first, shift_y_bounds[i].second, shift_z_bounds[i].first, shift_z_bounds[i].second, model_voronoi_normals, model_face_intercept, model_voronoi_vertices, new_model_pts, new_model_ids, dist_cost, 0, nb_threads);
        if(valid_cells.size()!=0){
            auto model = build_norm2_SOC_MIQCP(point_cloud_model, point_cloud_data, valid_cells[i], N2, dist_cost, roll_min, roll_max,  pitch_min, pitch_max, yaw_min, yaw_max, lb, ub, shift_min_y, shift_max_y, shift_min_z, shift_max_z, rot_trans[i], convex, incompatible_pairs, norm_x, norm_y, norm_z, intercept, init_matching, init_err_per_point, model_radius, relax_ints, relax_sdp, rigid_transf);
            model->write();
            models.push_back(model);
        }
    }
    run_parallel(models, gurobi, 1e-4, nb_threads, "", max_iter, max_time);
    for (int i = 0; i<models.size(); i++) {
        if(models[i]->_status==0){
            ub = models[i]->get_obj_val();
            lb = std::max(models[i]->get_rel_obj_val(), best_lb);
            if(lb<best_lb)
                best_lb = lb;
            if(ub<best_ub){
                best_ub = ub;
                bool is_rotation  = get_solution(models[i], best_rot_trans, new_matching);
            }
            lb_queue.push(treenode(models[i], roll_bounds[i],  pitch_bounds[i], yaw_bounds[i], shift_x_bounds[i], shift_y_bounds[i], shift_z_bounds[i], lb, ub, valid_cells[i]));
        }
        else {
            DebugOn("Infeasible model after preprocessing!\n");
            nb_pruned++;
            DebugOn("Nb pruned =  " << nb_pruned << endl);
        }
    }
    
    double elapsed_time = get_wall_time() - time_start;
    double opt_gap = (best_ub - best_lb)/best_ub;
    double max_opt_gap = 0.05;/* 5% opt gap */
    double opt_gap_old=opt_gap;
    double eps=0.01;
    while (elapsed_time < total_time_max && !lb_queue.empty() && opt_gap > max_opt_gap) {
        best_lb = lb_queue.top().lb;
        opt_gap = (best_ub - best_lb)/best_ub;
        if(opt_gap_old-opt_gap <= eps){
            max_time*=2;
        }
        opt_gap_old=opt_gap;
        DebugOn("Best UB so far = " << to_string_with_precision(best_ub,6) << endl);
        DebugOn("Best LB so far = " << to_string_with_precision(best_lb,6) << endl);
        DebugOn("Opt gap so far = " << to_string_with_precision(opt_gap*100,6) << "%\n");
        DebugOn("Queue size = " << lb_queue.size() << "\n");
        DebugOn("Elapsed time = " << elapsed_time << "seconds\n");
        if(elapsed_time >= total_time_max || opt_gap <= max_opt_gap)
            break;
        // Discard all nodes with high lower bounds in the queue
        int old_size = lb_queue.size();
        priority_queue<treenode> new_lb_queue;
        while(!lb_queue.empty())
        {
            auto node = lb_queue.top();
            lb_queue.pop();
            if(node.lb < best_ub)
                new_lb_queue.push(node);
            else
                break;
        }
        lb_queue = new_lb_queue;
        if(old_size - lb_queue.size()>0)
            DebugOn("Just pruned " << old_size - lb_queue.size() << " node(s)\n");
        nb_pruned += old_size - lb_queue.size();
        DebugOn("Total pruned =  " << nb_pruned << endl);
        DebugOn("Queue size = " << lb_queue.size() << "\n");
        treenode topnode = lb_queue.top();
        lb_queue.pop();
        models.clear();
        for (int i = 0; i<nb_threads; i++) {
            double x_shift_increment = (topnode.tx.second - topnode.tx.first)/nb_threads;
            double max_incr = x_shift_increment;
            double y_shift_increment = (topnode.ty.second - topnode.ty.first)/nb_threads;
            if(y_shift_increment>max_incr)
                max_incr = y_shift_increment;
            double z_shift_increment = (topnode.tz.second - topnode.tz.first)/nb_threads;
            if(z_shift_increment>max_incr)
                max_incr = z_shift_increment;
            double roll_increment = (topnode.roll.second - topnode.roll.first)/nb_threads;
            if(roll_increment>max_incr)
                max_incr = roll_increment;
            double pitch_increment = (topnode.pitch.second - topnode.pitch.first)/nb_threads;
            if(pitch_increment>max_incr)
                max_incr = pitch_increment;
            double yaw_increment = (topnode.yaw.second - topnode.yaw.first)/nb_threads;
            if(yaw_increment>max_incr)
                max_incr = yaw_increment;
            shift_x_bounds[i] = {topnode.tx.first, topnode.tx.second};
            shift_y_bounds[i] = {topnode.ty.first, topnode.ty.second};
            shift_z_bounds[i] = {topnode.tz.first, topnode.tz.second};
            roll_bounds[i] = {topnode.roll.first, topnode.roll.second};
            pitch_bounds[i] = {topnode.pitch.first, topnode.pitch.second};
            yaw_bounds[i] = {topnode.yaw.first, topnode.yaw.second};
            if(max_incr==x_shift_increment){
                lb = topnode.tx.first + i*x_shift_increment;
                ub = topnode.tx.first+(i+1)*x_shift_increment;
                shift_x_bounds[i] = {lb, ub};
            }
            else if(max_incr==y_shift_increment){
                lb = topnode.ty.first + i*y_shift_increment;
                ub = topnode.ty.first+(i+1)*y_shift_increment;
                shift_y_bounds[i] = {shift_min_y, shift_max_y};
            }
            else if(max_incr==z_shift_increment){
                lb = topnode.tz.first + i*z_shift_increment;
                ub = topnode.tz.first+(i+1)*z_shift_increment;
                shift_z_bounds[i] = {shift_min_z, shift_max_z};
            }
            else if(max_incr==roll_increment){
                lb = topnode.roll.first + i*roll_increment;
                ub = topnode.roll.first+(i+1)*roll_increment;
                roll_bounds[i] = {lb, ub};
            }
            else if(max_incr==pitch_increment){
                lb = topnode.pitch.first + i*pitch_increment;
                ub = topnode.pitch.first+(i+1)*pitch_increment;
                pitch_bounds[i] = {lb, ub};
            }
            else{
                lb = topnode.yaw.first + i*yaw_increment;
                ub = topnode.yaw.first+(i+1)*yaw_increment;
                yaw_bounds[i] = {lb, ub};
            }
            valid_cells[i] = preprocess_poltyope_intersect(point_cloud_data, point_cloud_model, topnode.valid_cells, roll_bounds[i].first, roll_bounds[i].second,  pitch_bounds[i].first, pitch_bounds[i].second, yaw_bounds[i].first, yaw_bounds[i].second, shift_x_bounds[i].first, shift_x_bounds[i].second, shift_y_bounds[i].first, shift_y_bounds[i].second, shift_z_bounds[i].first, shift_z_bounds[i].second, model_voronoi_normals, model_face_intercept,model_voronoi_vertices, new_model_pts, new_model_ids, dist_cost, 0, nb_threads);
            
            
            auto m = build_norm2_SOC_MIQCP(point_cloud_model, point_cloud_data, valid_cells[i], N2, dist_cost, roll_bounds[i].first, roll_bounds[i].second,  pitch_bounds[i].first, pitch_bounds[i].second, yaw_bounds[i].first, yaw_bounds[i].second, shift_x_bounds[i].first, shift_x_bounds[i].second, shift_y_bounds[i].first, shift_y_bounds[i].second, shift_z_bounds[i].first, shift_z_bounds[i].second, rot_trans[i], convex, incompatible_pairs, norm_x, norm_y, norm_z, intercept, init_matching, init_err_per_point, model_radius, relax_ints, relax_sdp, rigid_transf);
            if(m->_status==0 &&  (valid_cells[i].size()!=0))
                models.push_back(m);
            else {
                DebugOn("Infeasible model\n");
                nb_pruned++;
                DebugOn("Total pruned =  " << nb_pruned << endl);
                DebugOn("Queue size = " << lb_queue.size() << "\n");
            }
        }
        elapsed_time = get_wall_time() - time_start;
        DebugOn("Elapsed time = " << elapsed_time << "seconds\n");
        if(elapsed_time + max_time > total_time_max)
            break;
        run_parallel(models, gurobi, 1e-4, nb_threads, "", max_iter, max_time);
        for (int i = 0; i<models.size(); i++) {
            if(models[i]->_status==0){
                ub = models[i]->get_obj_val();
                lb = std::max(models[i]->get_rel_obj_val(), best_lb);
                if(ub<best_ub){
                    best_ub = ub;
                    bool is_rotation  = get_solution(models[i], best_rot_trans, new_matching);
                }
                lb_queue.push(treenode(models[i], roll_bounds[i],  pitch_bounds[i], yaw_bounds[i], shift_x_bounds[i], shift_y_bounds[i], shift_z_bounds[i], lb, ub, valid_cells[i]));
            }
        }
        elapsed_time = get_wall_time() - time_start;
        //        max_time += 1;
    }
    return best_rot_trans;
}

vector<double> BranchBound2(vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, param<>& norm_x, param<>& norm_y, param<>& norm_z,  param<>& intercept, const vector<int>& init_matching, const vector<double>& init_err_per_point, param<>& model_radius, vector<vector<vector<double>>> model_voronoi_normals, vector<vector<double>> model_face_intercept, vector<vector<vector<double>>> model_voronoi_vertices, vector<int>& new_model_pts, indices& new_model_ids, param<>& dist_cost, bool relax_ints, bool relax_sdp, bool rigid_transf) {
    /* INPUT BOUNDS */
    double time_start = get_wall_time();
    double total_time_max = 9000;
    double shift_min_x =  -0.1, shift_max_x = 0.1, shift_min_y = -0.1,shift_max_y = 0.1,shift_min_z = -0.1,shift_max_z = 0.1;
    double yaw_min = -10*pi/180., yaw_max = 10*pi/180., pitch_min =-10*pi/180.,pitch_max = 10*pi/180.,roll_min =-10*pi/180.,roll_max = 10*pi/180.;
    indices N1 = range(1,point_cloud_data.size());
    indices N2 = range(1,point_cloud_model.size());
    vector<int> new_matching(N1.size());
    bool convex = false;
    double max_time = 30;
    double max_time_init=30;
    bool max_time_increase=false;
    int max_iter = 1e6;
    vector<pair<pair<int,int>,pair<int,int>>> incompatible_pairs;
    auto nb_threads = std::thread::hardware_concurrency()/2;
    //nb_threads = 2;
    vector<pair<double,double>> shift_x_bounds(nb_threads), shift_y_bounds(nb_threads), shift_z_bounds(nb_threads);
    vector<pair<double,double>> roll_bounds(nb_threads), pitch_bounds(nb_threads), yaw_bounds(nb_threads);
    vector<indices> valid_cells(nb_threads, indices("valid_cells"));
    vector<vector<double>> rot_trans(nb_threads, vector<double>(12));
    vector<double> best_rot_trans(12);
    DebugOn("I will be using " << nb_threads << " parallel threads" << endl);
    vector<shared_ptr<Model<>>> models;
    double lb = 0, ub = 12, ub_=12, best_lb = 0, best_ub = 12;
    int nb_pruned = 0;
    priority_queue<treenode> lb_queue;
    for (int i = 0; i<nb_threads; i++) {
        valid_cells[i] = indices(N1,N2);
        double x_shift_increment = (shift_max_x - shift_min_x)/nb_threads;
        lb = shift_min_x + i*x_shift_increment;
        ub = shift_min_x+(i+1)*x_shift_increment;
        shift_x_bounds[i] = {lb, ub};
        shift_y_bounds[i] = {shift_min_y, shift_max_y};
        shift_z_bounds[i] = {shift_min_z, shift_max_z};
        roll_bounds[i] = {roll_min, roll_max};
        pitch_bounds[i] = {pitch_min, pitch_max};
        yaw_bounds[i] = {yaw_min, yaw_max};
        valid_cells[i] = preprocess_poltyope_intersect(point_cloud_data, point_cloud_model, valid_cells[i], roll_bounds[i].first, roll_bounds[i].second,  pitch_bounds[i].first, pitch_bounds[i].second, yaw_bounds[i].first, yaw_bounds[i].second, shift_x_bounds[i].first, shift_x_bounds[i].second, shift_y_bounds[i].first, shift_y_bounds[i].second, shift_z_bounds[i].first, shift_z_bounds[i].second, model_voronoi_normals, model_face_intercept, model_voronoi_vertices, new_model_pts, new_model_ids, dist_cost, 0, nb_threads);
        DebugOff("v size "<<valid_cells[i].size()<<endl);
        if(valid_cells[i].size()!=0){
            auto model = build_linobj_convex(point_cloud_model, point_cloud_data, valid_cells[i], roll_bounds[i].first, roll_bounds[i].second,  pitch_bounds[i].first, pitch_bounds[i].second, yaw_bounds[i].first, yaw_bounds[i].second, shift_x_bounds[i].first, shift_x_bounds[i].second, shift_y_bounds[i].first, shift_y_bounds[i].second, shift_z_bounds[i].first, shift_z_bounds[i].second, rot_trans[i], false, incompatible_pairs, norm_x, norm_y, norm_z, intercept,init_matching, init_err_per_point, false, best_ub);
            model->write();
            models.push_back(model);
        }
    }
    run_parallel(models, gurobi, 1e-4, nb_threads, "", max_iter, max_time);
    for (int i = 0; i<models.size(); i++) {
        if(models[i]->_status==0){
            ub_ = models[i]->get_obj_val();
            if(ub_>=0){
                vector<double> solution_lb;
                vector<double> rot_trans_ub;
                rot_trans_ub.resize(12);
                auto binc = models[i]->get_var<int>("bin");
                auto v_keys=binc.get_keys();
                for(auto it_key=v_keys->begin(); it_key!=v_keys->end(); it_key++){
                    auto key = *it_key;
                    solution_lb.push_back(binc.eval(key));
                }
                ub= build_upperbound(solution_lb,point_cloud_model, point_cloud_data, *(binc._indices), roll_min, roll_max, pitch_min, pitch_max, yaw_min, yaw_max, shift_min_x, shift_max_x, shift_min_y, shift_max_y, shift_min_z, shift_max_z, rot_trans_ub);
                if(ub<best_ub){
                    best_ub = ub;
               // bool is_rotation  = get_solution(models[i], best_rot_trans, new_matching);
                }
            }
            lb = std::max(models[i]->get_rel_obj_val(), best_lb);
            if(lb<best_lb)
                best_lb = lb;
            
            lb_queue.push(treenode(models[i], roll_bounds[i],  pitch_bounds[i], yaw_bounds[i], shift_x_bounds[i], shift_y_bounds[i], shift_z_bounds[i], lb, ub, valid_cells[i]));
        }
        else {
            DebugOn("Infeasible model after preprocessing!\n");
            nb_pruned++;
            DebugOn("Nb pruned =  " << nb_pruned << endl);
        }
    }
    
    double elapsed_time = get_wall_time() - time_start;
    double opt_gap = (best_ub - best_lb)/best_ub;
    double max_opt_gap = 0.05;/* 5% opt gap */
    double opt_gap_old=opt_gap;
    double eps=0.001;
    while (elapsed_time < total_time_max && !lb_queue.empty() && opt_gap > max_opt_gap) {
        best_lb = lb_queue.top().lb;
        opt_gap = (best_ub - best_lb)/best_ub;
        if(opt_gap_old-opt_gap <= eps){
            max_time=std::min(max_time*2, 120.0);
            max_time_increase=true;
        }
        if(opt_gap_old-opt_gap > eps && max_time_increase){
            max_time=std::max(max_time/2.0, max_time_init);
            if(max_time==max_time_init){
                max_time_increase=false;
            }
        }
        opt_gap_old=opt_gap;
        DebugOn("Best UB so far = " << to_string_with_precision(best_ub,6) << endl);
        DebugOn("Best LB so far = " << to_string_with_precision(best_lb,6) << endl);
        DebugOn("Opt gap so far = " << to_string_with_precision(opt_gap*100,6) << "%\n");
        DebugOn("Queue size = " << lb_queue.size() << "\n");
        DebugOn("Elapsed time = " << elapsed_time << "seconds\n");
        if(elapsed_time >= total_time_max || opt_gap <= max_opt_gap)
            break;
        // Discard all nodes with high lower bounds in the queue
        int old_size = lb_queue.size();
        priority_queue<treenode> new_lb_queue;
        while(!lb_queue.empty())
        {
            auto node = lb_queue.top();
            lb_queue.pop();
            if(node.lb < best_ub)
                new_lb_queue.push(node);
            else
                break;
        }
        lb_queue = new_lb_queue;
        if(old_size - lb_queue.size()>0)
            DebugOn("Just pruned " << old_size - lb_queue.size() << " node(s)\n");
        nb_pruned += old_size - lb_queue.size();
        DebugOn("Total pruned =  " << nb_pruned << endl);
        DebugOn("Queue size = " << lb_queue.size() << "\n");
        treenode topnode = lb_queue.top();
        double max_incr=0;
        lb_queue.pop();
        models.clear();
        for (int i = 0; i<nb_threads; i++) {
            double x_shift_increment = (topnode.tx.second - topnode.tx.first)/nb_threads;
            max_incr = x_shift_increment;
            double y_shift_increment = (topnode.ty.second - topnode.ty.first)/nb_threads;
            if(y_shift_increment>max_incr)
                max_incr = y_shift_increment;
            double z_shift_increment = (topnode.tz.second - topnode.tz.first)/nb_threads;
            if(z_shift_increment>max_incr)
                max_incr = z_shift_increment;
            double roll_increment = (topnode.roll.second - topnode.roll.first)/nb_threads;
            if(roll_increment>max_incr)
                max_incr = roll_increment;
            double pitch_increment = (topnode.pitch.second - topnode.pitch.first)/nb_threads;
            if(pitch_increment>max_incr)
                max_incr = pitch_increment;
            double yaw_increment = (topnode.yaw.second - topnode.yaw.first)/nb_threads;
            if(yaw_increment>max_incr)
                max_incr = yaw_increment;
            shift_x_bounds[i] = {topnode.tx.first, topnode.tx.second};
            shift_y_bounds[i] = {topnode.ty.first, topnode.ty.second};
            shift_z_bounds[i] = {topnode.tz.first, topnode.tz.second};
            roll_bounds[i] = {topnode.roll.first, topnode.roll.second};
            pitch_bounds[i] = {topnode.pitch.first, topnode.pitch.second};
            yaw_bounds[i] = {topnode.yaw.first, topnode.yaw.second};
            if(max_incr==x_shift_increment){
                lb = topnode.tx.first + i*x_shift_increment;
                ub = topnode.tx.first+(i+1)*x_shift_increment;
                shift_x_bounds[i] = {lb, ub};
            }
            else if(max_incr==y_shift_increment){
                lb = topnode.ty.first + i*y_shift_increment;
                ub = topnode.ty.first+(i+1)*y_shift_increment;
                shift_y_bounds[i] = {lb,ub};
            }
            else if(max_incr==z_shift_increment){
                lb = topnode.tz.first + i*z_shift_increment;
                ub = topnode.tz.first+(i+1)*z_shift_increment;
                shift_z_bounds[i] = {lb,ub};
            }
            else if(max_incr==roll_increment){
                lb = topnode.roll.first + i*roll_increment;
                ub = topnode.roll.first+(i+1)*roll_increment;
                roll_bounds[i] = {lb, ub};
            }
            else if(max_incr==pitch_increment){
                lb = topnode.pitch.first + i*pitch_increment;
                ub = topnode.pitch.first+(i+1)*pitch_increment;
                pitch_bounds[i] = {lb, ub};
            }
            else{
                lb = topnode.yaw.first + i*yaw_increment;
                ub = topnode.yaw.first+(i+1)*yaw_increment;
                yaw_bounds[i] = {lb, ub};
            }
            if(max_incr==0){
                DebugOn("identical"<<endl);
            }
            else{
                DebugOn("max_incr "<<max_incr<<endl);
            }
            valid_cells[i] = preprocess_poltyope_intersect(point_cloud_data, point_cloud_model, topnode.valid_cells, roll_bounds[i].first, roll_bounds[i].second,  pitch_bounds[i].first, pitch_bounds[i].second, yaw_bounds[i].first, yaw_bounds[i].second, shift_x_bounds[i].first, shift_x_bounds[i].second, shift_y_bounds[i].first, shift_y_bounds[i].second, shift_z_bounds[i].first, shift_z_bounds[i].second, model_voronoi_normals, model_face_intercept,model_voronoi_vertices, new_model_pts, new_model_ids, dist_cost, 0, nb_threads);
            
            if((valid_cells[i].size()!=0)){
                auto m = build_linobj_convex(point_cloud_model, point_cloud_data, valid_cells[i], roll_bounds[i].first, roll_bounds[i].second,  pitch_bounds[i].first, pitch_bounds[i].second, yaw_bounds[i].first, yaw_bounds[i].second, shift_x_bounds[i].first, shift_x_bounds[i].second, shift_y_bounds[i].first, shift_y_bounds[i].second, shift_z_bounds[i].first, shift_z_bounds[i].second, rot_trans[i], false, incompatible_pairs, norm_x, norm_y, norm_z, intercept,init_matching, init_err_per_point, false, best_ub);
                models.push_back(m);
            }
            else {
                DebugOn("v size "<<valid_cells[i].size()<<endl);
                DebugOn("Infeasible model\n");
                nb_pruned++;
                DebugOn("Total pruned =  " << nb_pruned << endl);
                DebugOn("Queue size = " << lb_queue.size() << "\n");
            }
        }
        elapsed_time = get_wall_time() - time_start;
        DebugOn("Elapsed time = " << elapsed_time << "seconds\n");
        if(elapsed_time + max_time > total_time_max){
            DebugOn("max time "<< max_time);
            break;
        }
        run_parallel(models, gurobi, 1e-4, nb_threads, "", max_iter, max_time);
        for (int i = 0; i<models.size(); i++) {
            if(models[i]->_status==0){
                ub_ = models[i]->get_obj_val();
                if(ub_>=0){
                    vector<double> solution_lb;
                    vector<double> rot_trans_ub;
                    rot_trans_ub.resize(12);
                    auto binc = models[i]->get_var<int>("bin");
                    auto v_keys=binc.get_keys();
                    for(auto it_key=v_keys->begin(); it_key!=v_keys->end(); it_key++){
                        auto key = *it_key;
                        solution_lb.push_back(binc.eval(key));
                    }
                    ub= build_upperbound(solution_lb,point_cloud_model, point_cloud_data, *(binc._indices), roll_min, roll_max, pitch_min, pitch_max, yaw_min, yaw_max, shift_min_x, shift_max_x, shift_min_y, shift_max_y, shift_min_z, shift_max_z, rot_trans_ub);
                    if(ub<best_ub){
                        best_ub = ub;
                   // bool is_rotation  = get_solution(models[i], best_rot_trans, new_matching);
                    }
                }
                lb = std::max(models[i]->get_rel_obj_val(), best_lb);
                if(lb<best_lb)
                    best_lb = lb;
                lb_queue.push(treenode(models[i], roll_bounds[i],  pitch_bounds[i], yaw_bounds[i], shift_x_bounds[i], shift_y_bounds[i], shift_z_bounds[i], lb, ub, valid_cells[i]));
            }
        }
        elapsed_time = get_wall_time() - time_start;
        //        max_time += 1;
    }
    return best_rot_trans;
}
vector<double> BranchBound3(vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, param<>& norm_x, param<>& norm_y, param<>& norm_z,  param<>& intercept, const vector<int>& init_matching, const vector<double>& init_err_per_point, param<>& model_radius, vector<vector<vector<double>>> model_voronoi_normals, vector<vector<double>> model_face_intercept, vector<vector<vector<double>>> model_voronoi_vertices, vector<int>& new_model_pts, indices& new_model_ids, param<>& dist_cost, bool relax_ints, bool relax_sdp, bool rigid_transf) {
    /* INPUT BOUNDS */
    double time_start = get_wall_time();
    double total_time_max = 9000;
    double shift_min_x =  -0.1, shift_max_x = 0.1, shift_min_y = -0.1,shift_max_y = 0.1,shift_min_z = -0.1,shift_max_z = 0.1;
    double yaw_min = -10*pi/180., yaw_max = 10*pi/180., pitch_min =-10*pi/180.,pitch_max = 10*pi/180.,roll_min =-10*pi/180.,roll_max = 10*pi/180.;
    indices N1 = range(1,point_cloud_data.size());
    indices N2 = range(1,point_cloud_model.size());
    vector<int> new_matching(N1.size());
    bool convex = false;
    double max_time = 30;
    double max_time_init=30;
    bool max_time_increase=false;
    int max_iter = 1e6;
    int models_count=0;
    int infeasible_count=0;
    vector<pair<pair<int,int>,pair<int,int>>> incompatible_pairs;
    auto nb_threads = std::thread::hardware_concurrency();
    //nb_threads = 1;
    pair<double,double> shift_x_bounds_r, shift_y_bounds_r, shift_z_bounds_r, roll_bounds_r, pitch_bounds_r, yaw_bounds_r;
    shift_x_bounds_r={shift_min_x, shift_max_x};
    shift_y_bounds_r={shift_min_y, shift_max_y};
    shift_z_bounds_r={shift_min_z, shift_max_z};
    roll_bounds_r={roll_min, roll_max};
    pitch_bounds_r={pitch_min, pitch_max};
    yaw_bounds_r={yaw_min, yaw_max};
    vector<pair<double,double>> shift_x_bounds(nb_threads), shift_y_bounds(nb_threads), shift_z_bounds(nb_threads);
    vector<pair<double,double>> roll_bounds(nb_threads), pitch_bounds(nb_threads), yaw_bounds(nb_threads);
    vector<indices> valid_cells(nb_threads, indices("valid_cells"));
    vector<vector<double>> rot_trans(nb_threads, vector<double>(12));
    vector<double> rot_trans_r(12);
    vector<double> best_rot_trans(12);
    vector<double> sol_gur(12);
    vector<int> pos_vec(nb_threads);
    DebugOn("I will be using " << nb_threads << " parallel threads" << endl);
    vector<shared_ptr<Model<>>> models;
    double lb = 0, ub = 12, ub_=-1, best_lb = 0, best_ub = 12;
    int nb_pruned = 0;
    int depth_r=0, iter=0;
    priority_queue<treenode_n> lb_queue;
    auto valid_cells_ro=indices(N1,N2);
    auto valid_cells_r = preprocess_poltyope_intersect(point_cloud_data, point_cloud_model, valid_cells_ro, roll_min, roll_max, pitch_min, pitch_max, yaw_min, yaw_max, shift_min_x, shift_max_x, shift_min_y, shift_max_y, shift_min_z, shift_max_z, model_voronoi_normals, model_face_intercept, model_voronoi_vertices, new_model_pts, new_model_ids, dist_cost, 0, nb_threads);
    //auto model = build_linobj_convex(point_cloud_model, point_cloud_data, valid_cells_r, roll_min, roll_max, pitch_min, pitch_max, yaw_min, yaw_max, shift_min_x, shift_max_x, shift_min_y, shift_max_y, shift_min_z, shift_max_z, rot_trans_r, false, incompatible_pairs, norm_x, norm_y, norm_z, intercept,init_matching, init_err_per_point, false, best_ub);
    auto model=build_polyhedral_ipopt(point_cloud_model, point_cloud_data, valid_cells_r, roll_min, roll_max,pitch_min,pitch_max,  yaw_min, yaw_max,  shift_min_x, shift_max_x, shift_min_y, shift_max_y, shift_min_z, shift_max_z,  rot_trans_r, false,  incompatible_pairs,  norm_x, norm_y,  norm_z,  intercept, init_matching,  init_err_per_point,  model_voronoi_normals,  model_face_intercept,  false);
 
    lb_queue.push(treenode_n(model, roll_bounds_r, pitch_bounds_r, yaw_bounds_r, shift_x_bounds_r, shift_y_bounds_r, shift_z_bounds_r, lb, ub, ub_, depth_r, valid_cells_r));
    //auto goicp=initialize_ICP_only(point_cloud_model, point_cloud_data);

    double elapsed_time = get_wall_time() - time_start;
    double opt_gap = (best_ub - best_lb)/best_ub;
    double max_opt_gap = 0.05;/* 5% opt gap */
    double opt_gap_old=opt_gap+10;
    double eps=0.01;
    while (elapsed_time < total_time_max && !lb_queue.empty() && opt_gap > max_opt_gap) {
        best_lb = lb_queue.top().lb;
        opt_gap = (best_ub - best_lb)/best_ub;
        if(opt_gap_old-opt_gap <= eps){
            max_time=std::min(max_time*2, 120.0);
            max_time_increase=true;
        }
        if(opt_gap_old-opt_gap > eps && max_time_increase){
            max_time=std::max(max_time/2.0, max_time_init);
            if(max_time==max_time_init){
                max_time_increase=false;
            }
        }
        opt_gap_old=opt_gap;
        DebugOn("Best UB so far = " << to_string_with_precision(best_ub,6) << endl);
        DebugOn("Best LB so far = " << to_string_with_precision(best_lb,6) << endl);
        DebugOn("Opt gap so far = " << to_string_with_precision(opt_gap*100,6) << "%\n");
        DebugOn("Queue size = " << lb_queue.size() << "\n");
        DebugOn("Elapsed time = " << elapsed_time << "seconds\n");
        DebugOn("iter "<<iter<<endl);
        if(elapsed_time >= total_time_max || opt_gap <= max_opt_gap)
            break;
        // Discard all nodes with high lower bounds in the queue
        int old_size = lb_queue.size();
        priority_queue<treenode_n> new_lb_queue;
        while(!lb_queue.empty())
        {
            auto node = lb_queue.top();
            lb_queue.pop();
            if(node.lb < best_ub)
                new_lb_queue.push(node);
        }
        lb_queue = new_lb_queue;
        if(old_size - lb_queue.size()>0){
            DebugOn("Just pruned " << old_size - lb_queue.size() << " node(s)\n");
            nb_pruned += old_size - lb_queue.size();
        }
        DebugOn("Total infeasible =  " << infeasible_count << endl);
        DebugOn("Total pruned =  " << nb_pruned << endl);
        DebugOn("Queue size = " << lb_queue.size() << "\n");
        
        double max_incr=0, max_ratio=1;
        pos_vec.clear();
        models.clear();
        models_count=0;
        iter++;
        for (int i = 0; i<nb_threads; i+=2) {
            treenode_n topnode = lb_queue.top();
            lb_queue.pop();
            if(topnode.ub_>0 && false){
                auto model_parent= topnode.mod;
                vector<double> sol(6);
                get_angle_rotation_transl_matrix(model_parent, sol);
                auto roll_=sol[0];auto pitch_=sol[1];auto yaw_=sol[2];
                auto tx_=sol[3];auto ty_=sol[4];auto tz_=sol[5];
                max_ratio=1;
                double x_shift_increment = std::abs(0.5-(topnode.tx.second-tx_)/(topnode.tx.second-topnode.tx.first));
                if(x_shift_increment<=0.4){
                    max_ratio = x_shift_increment;
                }
                double y_shift_increment = std::abs(0.5-(topnode.ty.second-ty_)/(topnode.ty.second-topnode.ty.first));
                if(y_shift_increment<=0.4 && y_shift_increment>max_ratio)
                    max_ratio = y_shift_increment;
                double z_shift_increment = std::abs(0.5-(topnode.tz.second-tz_)/(topnode.tz.second-topnode.tz.first));
                if(z_shift_increment<=0.4 && z_shift_increment>max_ratio)
                    max_ratio = z_shift_increment;
                double roll_increment = std::abs(0.5-(topnode.roll.second-roll_)/(topnode.roll.second-topnode.roll.first));
                if(roll_increment<=0.4 && roll_increment>max_ratio)
                    max_ratio = roll_increment;
                double pitch_increment = std::abs(0.5-(topnode.pitch.second-pitch_)/(topnode.pitch.second-topnode.pitch.first));
                if(pitch_increment<=0.4 && pitch_increment>max_ratio)
                    max_ratio = pitch_increment;
                double yaw_increment = std::abs(0.5-(topnode.yaw.second-yaw_)/(topnode.yaw.second-topnode.yaw.first));
                if(yaw_increment<=0.4 && yaw_increment>max_ratio)
                    max_ratio = yaw_increment;
                shift_x_bounds[i] = {topnode.tx.first, topnode.tx.second};
                shift_y_bounds[i] = {topnode.ty.first, topnode.ty.second};
                shift_z_bounds[i] = {topnode.tz.first, topnode.tz.second};
                roll_bounds[i] = {topnode.roll.first, topnode.roll.second};
                pitch_bounds[i] = {topnode.pitch.first, topnode.pitch.second};
                yaw_bounds[i] = {topnode.yaw.first, topnode.yaw.second};
                shift_x_bounds[i+1] = {topnode.tx.first, topnode.tx.second};
                shift_y_bounds[i+1] = {topnode.ty.first, topnode.ty.second};
                shift_z_bounds[i+1] = {topnode.tz.first, topnode.tz.second};
                roll_bounds[i+1] = {topnode.roll.first, topnode.roll.second};
                pitch_bounds[i+1] = {topnode.pitch.first, topnode.pitch.second};
                yaw_bounds[i+1] = {topnode.yaw.first, topnode.yaw.second};
                if(max_ratio==x_shift_increment){
                    shift_x_bounds[i] = {topnode.tx.first, tx_};
                    shift_x_bounds[i+1] = {tx_,topnode.tx.second};
                }
                else if(max_ratio==y_shift_increment){
                    shift_y_bounds[i] = {topnode.ty.first, ty_};
                    shift_y_bounds[i+1] = {ty_,topnode.ty.second};
                }
                else if(max_ratio==z_shift_increment){
                    shift_z_bounds[i] = {topnode.tz.first, tz_};
                    shift_z_bounds[i+1] = {tz_,topnode.tz.second};
                }
                else if(max_ratio==roll_increment){
                    roll_bounds[i] = {topnode.roll.first, roll_};
                    roll_bounds[i+1] = {roll_,topnode.roll.second};
                }
                else if(max_ratio==pitch_increment){
                    pitch_bounds[i] = {topnode.pitch.first, pitch_};
                    pitch_bounds[i+1] = {pitch_,topnode.pitch.second};
                }
                else{
                    yaw_bounds[i] = {topnode.yaw.first, yaw_};
                    yaw_bounds[i+1] = {yaw_,topnode.yaw.second};
                }
            }
            if((topnode.ub_<=0 || max_ratio==1) || true){
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
                shift_x_bounds[i] = {topnode.tx.first, topnode.tx.second};
                shift_y_bounds[i] = {topnode.ty.first, topnode.ty.second};
                shift_z_bounds[i] = {topnode.tz.first, topnode.tz.second};
                roll_bounds[i] = {topnode.roll.first, topnode.roll.second};
                pitch_bounds[i] = {topnode.pitch.first, topnode.pitch.second};
                yaw_bounds[i] = {topnode.yaw.first, topnode.yaw.second};
                shift_x_bounds[i+1] = {topnode.tx.first, topnode.tx.second};
                shift_y_bounds[i+1] = {topnode.ty.first, topnode.ty.second};
                shift_z_bounds[i+1] = {topnode.tz.first, topnode.tz.second};
                roll_bounds[i+1] = {topnode.roll.first, topnode.roll.second};
                pitch_bounds[i+1] = {topnode.pitch.first, topnode.pitch.second};
                yaw_bounds[i+1] = {topnode.yaw.first, topnode.yaw.second};
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
                else{
                    yaw_bounds[i] = {topnode.yaw.first, topnode.yaw.first+yaw_increment};
                    yaw_bounds[i+1] = {topnode.yaw.first+yaw_increment, topnode.yaw.second};
                }
                if(max_incr==0){
                    DebugOn("Warn: Identical child nodes"<<endl);
                }
                else{
                    DebugOff("max_incr "<<max_incr<<endl);
                }
            }
            
            valid_cells[i] = preprocess_poltyope_intersect(point_cloud_data, point_cloud_model, topnode.valid_cells, roll_bounds[i].first, roll_bounds[i].second,  pitch_bounds[i].first, pitch_bounds[i].second, yaw_bounds[i].first, yaw_bounds[i].second, shift_x_bounds[i].first, shift_x_bounds[i].second, shift_y_bounds[i].first, shift_y_bounds[i].second, shift_z_bounds[i].first, shift_z_bounds[i].second, model_voronoi_normals, model_face_intercept,model_voronoi_vertices, new_model_pts, new_model_ids, dist_cost, 0, nb_threads);
            
            if((valid_cells[i].size()!=0)){
                //auto m = build_linobj_convex(point_cloud_model, point_cloud_data, valid_cells[i], roll_bounds[i].first, roll_bounds[i].second,  pitch_bounds[i].first, pitch_bounds[i].second, yaw_bounds[i].first, yaw_bounds[i].second, shift_x_bounds[i].first, shift_x_bounds[i].second, shift_y_bounds[i].first, shift_y_bounds[i].second, shift_z_bounds[i].first, shift_z_bounds[i].second, rot_trans[i], false, incompatible_pairs, norm_x, norm_y, norm_z, intercept,init_matching, init_err_per_point, false, best_ub);
                auto m=build_polyhedral_ipopt(point_cloud_model, point_cloud_data, valid_cells[i], roll_bounds[i].first, roll_bounds[i].second,  pitch_bounds[i].first, pitch_bounds[i].second, yaw_bounds[i].first, yaw_bounds[i].second, shift_x_bounds[i].first, shift_x_bounds[i].second, shift_y_bounds[i].first, shift_y_bounds[i].second, shift_z_bounds[i].first, shift_z_bounds[i].second, rot_trans[i], false,  incompatible_pairs,  norm_x, norm_y,  norm_z,  intercept, init_matching,  init_err_per_point,  model_voronoi_normals,  model_face_intercept,  false);
                if(topnode.ub_>0 && max_ratio<1){
                   // initialize_model_from_parent(topnode.mod, m);
                }
                models.push_back(m);
                pos_vec.push_back(i);
            }
            else {
                DebugOff("v size "<<valid_cells[i].size()<<endl);
                DebugOff("Infeasible model\n");
                nb_pruned++;
                DebugOff("Total pruned =  " << nb_pruned << endl);
                DebugOff("Queue size = " << lb_queue.size() << "\n");
            }
            valid_cells[i+1] = preprocess_poltyope_intersect(point_cloud_data, point_cloud_model, topnode.valid_cells, roll_bounds[i+1].first, roll_bounds[i+1].second,  pitch_bounds[i+1].first, pitch_bounds[i+1].second, yaw_bounds[i+1].first, yaw_bounds[i+1].second, shift_x_bounds[i+1].first, shift_x_bounds[i+1].second, shift_y_bounds[i+1].first, shift_y_bounds[i+1].second, shift_z_bounds[i+1].first, shift_z_bounds[i+1].second, model_voronoi_normals, model_face_intercept,model_voronoi_vertices, new_model_pts, new_model_ids, dist_cost, 0, nb_threads);
            
            if((valid_cells[i+1].size()!=0)){
                //auto m = build_linobj_convex(point_cloud_model, point_cloud_data, valid_cells[i+1], roll_bounds[i+1].first, roll_bounds[i+1].second,  pitch_bounds[i+1].first, pitch_bounds[i+1].second, yaw_bounds[i+1].first, yaw_bounds[i+1].second, shift_x_bounds[i+1].first, shift_x_bounds[i+1].second, shift_y_bounds[i+1].first, shift_y_bounds[i+1].second, shift_z_bounds[i+1].first, shift_z_bounds[i+1].second, rot_trans[i+1], false, incompatible_pairs, norm_x, norm_y, norm_z, intercept,init_matching, init_err_per_point, false, best_ub);
                auto m=build_polyhedral_ipopt(point_cloud_model, point_cloud_data, valid_cells[i+1], roll_bounds[i+1].first, roll_bounds[i+1].second,  pitch_bounds[i+1].first, pitch_bounds[i+1].second, yaw_bounds[i+1].first, yaw_bounds[i+1].second, shift_x_bounds[i+1].first, shift_x_bounds[i+1].second, shift_y_bounds[i+1].first, shift_y_bounds[i+1].second, shift_z_bounds[i+1].first, shift_z_bounds[i+1].second, rot_trans[i+1], false,  incompatible_pairs,  norm_x, norm_y,  norm_z,  intercept, init_matching,  init_err_per_point,  model_voronoi_normals,  model_face_intercept,  false);
                if(topnode.ub_>0 && max_ratio<1){
                    //initialize_model_from_parent(topnode.mod, m);
                }
                models.push_back(m);
                pos_vec.push_back(i+1);
            }
            else {
                DebugOff("v size "<<valid_cells[i].size()<<endl);
                DebugOff("Infeasible model\n");
                nb_pruned++;
                DebugOff("Total pruned =  " << nb_pruned << endl);
                DebugOff("Queue size = " << lb_queue.size() << "\n");
            }
            if(lb_queue.empty()){
                break;
            }
        }
        elapsed_time = get_wall_time() - time_start;
        DebugOn("Elapsed time = " << elapsed_time << "seconds\n");
        if(elapsed_time + max_time > total_time_max){
            DebugOn("max time "<< max_time);
            break;
        }
        run_parallel(models, ipopt, 1e-4, nb_threads, "", max_iter, max_time);
        for (int i = 0; i<models.size(); i++) {
            if(models[i]->_status==0){
                ub_ = models[i]->get_obj_val();
                int pos=pos_vec[i];
                    if(ub_<=best_ub && iter%5==0){
                        vector<double> solution_lb;
                        vector<double> rot_trans_ub;
                        rot_trans_ub.resize(12);
                        auto binc = models[i]->get_var<double>("bin");
                        auto v_keys=binc.get_keys();
                        for(auto it_key=v_keys->begin(); it_key!=v_keys->end(); it_key++){
                            auto key = *it_key;
                            solution_lb.push_back(binc.eval(key));
                        }
                        ub= build_upperbound(solution_lb,point_cloud_model, point_cloud_data, *(binc._indices), roll_min, roll_max, pitch_min, pitch_max, yaw_min, yaw_max, shift_min_x, shift_max_x, shift_min_y, shift_max_y, shift_min_z, shift_max_z, rot_trans_ub);
                        if(ub<best_ub){
                            best_ub = ub;
                            best_rot_trans=rot_trans_ub;
                            //bool is_rotation  = get_solution(models[pos], sol_gur, new_matching);
                        }
//                        ub= run_ICP_only(goicp, roll_bounds[pos].first, roll_bounds[pos].second,  pitch_bounds[pos].first, pitch_bounds[pos].second, yaw_bounds[pos].first, yaw_bounds[pos].second, shift_x_bounds[pos].first, shift_x_bounds[pos].second, shift_y_bounds[pos].first, shift_y_bounds[pos].second, shift_z_bounds[pos].first, shift_z_bounds[pos].second,rot_trans_ub);
//                        if(ub<best_ub){
//                            best_ub = ub;
//                            best_rot_trans=rot_trans_ub;
//                            //bool is_rotation  = get_solution(models[pos], sol_gur, new_matching);
//                        }
                    }
               // lb = models[i]->get_rel_obj_val();
                lb = models[i]->get_obj_val();
                //lb = std::max(models[i]->get_rel_obj_val(), best_lb);
                if(lb<best_lb)
                    best_lb = lb;
                lb_queue.push(treenode_n(models[i], roll_bounds[pos],  pitch_bounds[pos], yaw_bounds[pos], shift_x_bounds[pos], shift_y_bounds[pos], shift_z_bounds[pos], lb, ub, ub_, iter, valid_cells[pos]));
            }
            else{
                infeasible_count++;
            }
        }
        best_lb = lb_queue.top().lb;
        opt_gap = (best_ub - best_lb)/best_ub;
    
        elapsed_time = get_wall_time() - time_start;
        //        max_time += 1;
    }
    DebugOn("UB final "<<best_ub<<endl);
    DebugOn("LB final "<<best_lb<<endl);
    DebugOn("Gap final "<<(best_ub-best_lb)/best_ub*100.0<<endl);

    vector<double> res(N1.size());
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

/*
 
 auto siny_sinp_range = get_product_range(siny._range, sinp._range);
 var<> siny_sinp("siny_sinp", siny_sinp_range->first, siny_sinp_range->second);
 auto cosy_sinp_range = get_product_range(cosy._range, sinp._range);
 var<> cosy_sinp("cosy_sinp", cosy_sinp_range->first, cosy_sinp_range->second);
 auto cosy_cosr_range = get_product_range(cosy._range, cosr._range);
 var<> cosy_cosr("cosy_cosr", cosy_cosr_range->first, cosy_cosr_range->second);
 auto cosy_sinr_sinp_range = get_product_range(cosy_sinr._range, sinp._range);
 var<> cosy_sinr_sinp("cosy_sinr_sinp", cosy_sinr_sinp_range->first, cosy_sinr_sinp_range->second);
 auto cosy_cosp_range = get_product_range(cosy._range, cosp._range);
 var<> cosy_cosp("cosy_cosp", cosy_cosp_range->first, cosy_cosp_range->second);
 auto siny_cosp_range = get_product_range(siny._range, cosp._range);
 var<> siny_cosp("siny_cosp", siny_cosp_range->first, siny_cosp_range->second);
 auto cosy_sinr_cosp_range = get_product_range(cosy_sinr._range, cosp._range);
 var<> cosy_sinr_cosp("cosy_sinr_cosp", cosy_sinr_cosp_range->first, cosy_sinr_cosp_range->second);
 auto siny_cosr_range = get_product_range(siny._range, cosr._range);
 var<> siny_cosr("siny_cosr", siny_cosr_range->first, siny_cosr_range->second);
 auto siny_sinr_sinp_range = get_product_range(siny_sinr._range, sinp._range);
 var<> siny_sinr_sinp("siny_sinr_sinp", siny_sinr_sinp_range->first, siny_sinr_sinp_range->second);
 auto siny_sinr_cosp_range = get_product_range(siny_sinr._range, cosp._range);
 var<> siny_sinr_cosp("siny_sinr_cosp", siny_sinr_cosp_range->first, siny_sinr_cosp_range->second);
 auto cosr_sinp_range = get_product_range(cosr._range, sinp._range);
 var<> cosr_sinp("cosr_sinp", cosr_sinp_range->first, cosr_sinp_range->second);
 auto cosr_cosp_range = get_product_range(cosr._range, cosp._range);
 var<> cosr_cosp("cosr_cosp", cosr_cosp_range->first, cosr_cosp_range->second);
 
 func<> r11 = cos(yaw)*cos(roll);r11.eval_all();
 func<> r12 = cos(yaw)*sin(roll)*sin(pitch) - sin(yaw)*cos(pitch);r12.eval_all();
 func<> r13 = cos(yaw)*sin(roll)*cos(pitch) + sin(yaw)*sin(pitch);r13.eval_all();
 func<> r21 = sin(yaw)*cos(roll);r21.eval_all();
 func<> r22 = sin(yaw)*sin(roll)*sin(pitch) + cos(yaw)*cos(pitch);r22.eval_all();
 func<> r23 = sin(yaw)*sin(roll)*cos(pitch) - cos(yaw)*sin(pitch);r23.eval_all();
 func<> r31 = sin(-1*roll);r31.eval_all();
 func<> r32 = cos(roll)*sin(pitch);r32.eval_all();
 func<> r33 = cos(roll)*cos(pitch);r33.eval_all();
 
 
 */
