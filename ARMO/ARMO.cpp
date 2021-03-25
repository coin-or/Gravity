    //  ARMO: Alignment and Registration via Mathematical Optimization
    //  ARMO.cpp
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
#include "lasreader.hpp"
#include "laswriter.hpp"
#include "ARMO.hpp"
#include <gravity/KDTreeVectorOfVectorsAdaptor.h>
#include <time.h>
using namespace std;


#ifdef USE_VORO
#include "voro++.hh"
using namespace voro;

#endif

#define DEFAULT_OUTPUT_FNAME "output.txt"
#define DEFAULT_CONFIG_FNAME "config.txt"
#define DEFAULT_MODEL_FNAME "model.txt"
#define DEFAULT_DATA_FNAME "data.txt"



int main (int argc, char * argv[])
{
    string prob_type = "Reg";/* Problem type "Reg" for registration and "Align" for boresight alignment */
    string non_prop_scaling = "noscale";/* if set to "scale", will allow non-proportional scaling */
    double perc_outliers = 0;/* Percentage of outliers */
    if(argc>1){
        prob_type = argv[1];
    }
    bool Registration = prob_type=="Reg";/* Solve the Registration problem */
    if(Registration){
        double x_min = -1, x_max = 1, y_min = -1, y_max = 1, z_min = -1, z_max = 1;
        vector<vector<double>> initial_point_cloud_model, initial_point_cloud_data;
        vector<vector<double>> point_cloud_model, point_cloud_data;
        string Model_file = string(prj_dir)+"/data_sets/Registration/toy_model.txt";
        string Data_file = string(prj_dir)+"/data_sets/Registration/toy_data.txt";
        string algo = "MIQCP", global_str = "global", convex_str = "nonconvex", reform_str="no", obbt_str="yes", norm_str="norm2";
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
            non_prop_scaling = argv[5];
        }
        if(argc>6){
            perc_outliers = atof(argv[6]);
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
        vector<pair<double,double>> bounds(3);/* Bounds on x y and z of model points */
        DebugOn("Model file has " << model_nb_rows << " rows" << endl);
        DebugOn("Data file has " << data_nb_rows << " rows" << endl);
        for (int i = 0; i< model_nb_rows; i++) { // Input iterator
            auto x = Model_doc.GetCell<double>(0, i);
            auto y = Model_doc.GetCell<double>(1, i);
            auto z = Model_doc.GetCell<double>(2, i);
            if(bounds[0].first > x)
                bounds[0].first = x;
            if(bounds[0].second < x)
                bounds[0].second = x;
            if(bounds[1].first > y)
                bounds[1].first = y;
            if(bounds[1].second < y)
                bounds[1].second = y;
            if(bounds[2].first > z)
                bounds[2].first = z;
            if(bounds[2].second < z)
                bounds[2].second = z;
            initial_point_cloud_model.push_back({x,y,z});
            point_cloud_model.push_back({x,y,z});
        }
        for (int i = 0; i< data_nb_rows; i++) { // Input iterator
            auto x = Data_doc.GetCell<double>(0, i);
            auto y = Data_doc.GetCell<double>(1, i);
            auto z = Data_doc.GetCell<double>(2, i);
            initial_point_cloud_data.push_back({x,y,z});
            point_cloud_data.push_back({x,y,z});
        }
        data_nb_rows=point_cloud_data.size();
        model_nb_rows=point_cloud_model.size();
        int reduced_nb_data = 50;
        int reduced_nb_model = 50;
        bool subsample = false;
        if (subsample) {
            model_nb_rows = reduced_nb_model;
            data_nb_rows = reduced_nb_data;
            point_cloud_model = get_n_extreme_points(reduced_nb_model, point_cloud_model);
            point_cloud_data = get_n_extreme_points(reduced_nb_data, point_cloud_data);
        }
        center_point_cloud(point_cloud_model);
        center_point_cloud(point_cloud_data);
        center_point_cloud(initial_point_cloud_model);
        center_point_cloud(initial_point_cloud_data);
#ifdef USE_VORO
        container model_con(x_min,x_max,y_min,y_max,z_min,z_max,10,10,10,false,false,false,8);
        for (int i = 0; i< model_nb_rows; i++) { // Input iterator
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
        int total_nb_faces = 0;
        if(cl.start()) do if(model_con.compute_cell(c,cl)) {
            cl.pos(x,y,z);
            idx=cl.pid();
            c.vertices(x,y,z,v);
            int nb_vertices = v.size()/3;
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
                v_idx += 3;
                vertices[0]=v[v_idx];
                vertices[1]=v[v_idx+1];
                vertices[2]=v[v_idx+2];
                model_voronoi_vertices[idx][i]=vertices;
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
                    DebugOn("WARNING: model point " << idx << " on a voronoi face!");
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
            }
            model_voronoi_in_radius[idx] = min_dist;/* the radius of the largest ball in the voronoi cell is the distance from the center to the nearest facet */
        } while (cl.inc());
        DebugOn("The total number of faces = " << total_nb_faces << endl);
        
        int nb_reduced_model = std::min(model_nb_rows,150);
#endif
        bool norm1 = norm_str=="norm1";
        tuple<double,double,double,double,double,double,double> res_icp;
        tuple<double,double,double,double,double,double> res, res1, res2;
        vector<int> L2matching(data_nb_rows), L1matching(data_nb_rows);
        vector<double> L2err_per_point(data_nb_rows), L1err_per_point(data_nb_rows);
        auto L2error_init = computeL2error(point_cloud_model,point_cloud_data,L2matching,L2err_per_point);
        DebugOn("Initial L2 = " << L2error_init << endl);
        
        bool GoICP = (algo=="GoICP");
        bool MIQCP = (algo=="MIQCP");
                
        if(GoICP){/* Run GoICP inline */
            res_icp = run_GoICP(point_cloud_model, point_cloud_data);
            auto roll = get<0>(res_icp);auto pitch = get<1>(res_icp);auto yaw = get<2>(res_icp);auto x_shift = get<3>(res_icp);auto y_shift = get<4>(res_icp);auto z_shift = get<5>(res_icp);
            apply_rot_trans(roll, pitch, yaw, x_shift, y_shift, z_shift, point_cloud_data);
            apply_rot_trans(roll, pitch, yaw, x_shift, y_shift, z_shift, initial_point_cloud_data);
        }
        else if(MIQCP){
            vector<double> rot_trans(12);
            double shift_min_x =  -0.1, shift_max_x = 0.1, shift_min_y = -0.1,shift_max_y = 0.1,shift_min_z = -0.1,shift_max_z = 0.1;
            double yaw_min = -10*pi/180., yaw_max = 10*pi/180., pitch_min =-10*pi/180.,pitch_max = 10*pi/180.,roll_min =-10*pi/180.,roll_max = 10*pi/180.;
            /* Use wider bounds for small instances */
            if(data_nb_rows<100){
                shift_min_x =  -0.25; shift_max_x = 0.25; shift_min_y = -0.25; shift_max_y = 0.25; shift_min_z = -0.25; shift_max_z = 0.25;
                yaw_min = -120*pi/180.; yaw_max = 120*pi/180.; pitch_min =-120*pi/180.; pitch_max = 120*pi/180.; roll_min =-120*pi/180.; roll_max=120*pi/180.;
            }
            /* Check if 2D data */
            bool fixed_x = bounds[0].second - bounds[0].first < 1e-6;
            if(fixed_x){
                yaw_min = 0;
                yaw_max = 0;
                pitch_min = 0;
                pitch_max = 0;
                shift_min_x = 0;
                shift_max_x = 0;
            }
            bool fixed_y = bounds[1].second - bounds[1].first < 1e-6;
            if(fixed_y){
                roll_min = 0;
                roll_max = 0;
                yaw_min = 0;
                yaw_max = 0;
                shift_min_y = 0;
                shift_max_y = 0;
            }
            bool fixed_z = bounds[2].second - bounds[2].first < 1e-6;
            if(fixed_z){
                roll_min = 0;
                roll_max = 0;
                pitch_min = 0;
                pitch_max = 0;
                shift_min_z = 0;
                shift_max_z = 0;
            }
            
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
            valid_cells = indices(N1,N2);
            bool preprocess = false, nonprop_scale = non_prop_scaling=="scale";
            if(preprocess){
                double scale = 1;
                if(nonprop_scale)
                    scale = 1.4;/* Allow for 40% scaling in each axis */
                valid_cells=preprocess_QP(point_cloud_data, point_cloud_model, valid_cells, roll_min, roll_max, pitch_min, pitch_max, yaw_min, yaw_max, shift_min_x, shift_max_x, shift_min_y, shift_max_y, shift_min_z, shift_max_z, model_voronoi_normals, model_face_intercept, new_model_pts, new_model_ids, dist_cost, upper_bound, nb_total_threads, scale);
                double time_end = get_wall_time();
                auto prep_time = time_end - time_start;
                /* Terminal output */
                DebugOn("Preprocessing time = " << prep_time << endl);
            }
            vector<pair<pair<int,int>,pair<int,int>>> incompatibles;
            bool convex = false, relax_integers = false, relax_sdp = false;
            auto NC_SOC_MIQCP = build_norm2_SOC_MIQCP(point_cloud_model, point_cloud_data, valid_cells, new_model_ids, dist_cost, roll_min, roll_max,  pitch_min, pitch_max, yaw_min, yaw_max, shift_min_x, shift_max_x, shift_min_y, shift_max_y, shift_min_z, shift_max_z, rot_trans, convex, incompatibles, norm_x, norm_y, norm_z, intercept, L2matching, L2err_per_point, model_radius, relax_integers, relax_sdp, nonprop_scale, perc_outliers);
            double time_limit = 120;
#ifdef USE_GUROBI
            solver<> S(NC_SOC_MIQCP,gurobi);
            S.use_callback();
            S.run(5,1e-6,time_limit,1000);
#else
            DebugOn("WARNING: this version was compiled without Gurobi, please install Gurobi and rerun cmake.\n");
            return 0;
#endif
            get_solution(NC_SOC_MIQCP, rot_trans, L2matching);
            
            apply_rot_trans(rot_trans, point_cloud_data);
            apply_rot_trans(rot_trans, initial_point_cloud_data);
            
#ifdef USE_MATPLOT
            plot(initial_point_cloud_model,initial_point_cloud_data,1);
#endif
            DebugOn("L2 before optimization = " << to_string_with_precision(L2error_init,12) << endl);
            auto L2error_final = computeL2error(point_cloud_model,point_cloud_data,L2matching,L2err_per_point);
            auto L1error_final = computeL1error(point_cloud_model,point_cloud_data,L1matching,L1err_per_point);
            auto L2error_final_full = computeL2error(initial_point_cloud_model,initial_point_cloud_data,L2matching,L2err_per_point);
            auto L1error_final_full = computeL1error(initial_point_cloud_model,initial_point_cloud_data,L1matching,L1err_per_point);
            DebugOn("L2 after optimization on downsampled = " << to_string_with_precision(L2error_final,12) << endl);
            DebugOn("L1 after optimization on downsampled = " << to_string_with_precision(L1error_final,12) << endl);
            DebugOn("L2 after optimization on full set = " << to_string_with_precision(L2error_final_full,12) << endl);
            DebugOn("L1 after optimization on full set = " << to_string_with_precision(L1error_final_full,12) << endl);
        }
        return 0;
    }
    
    /* Boresight Alignment Problem */
    vector<vector<double>> full_point_cloud_model, full_point_cloud_data;
    vector<vector<double>> point_cloud_model, point_cloud_data;
    vector<vector<double>> full_uav_model, full_uav_data;
    vector<vector<double>> uav_model, uav_data;
    string Model_file = string(prj_dir)+"/data_sets/LiDAR/Cars_model.txt";
    string Data_file = string(prj_dir)+"/data_sets/LiDAR/Cars_data.txt";
    string red_Model_file = string(prj_dir)+"/data_sets/LiDAR/Cars_model_sub.txt";
    string red_Data_file = string(prj_dir)+"/data_sets/LiDAR/Cars_data_sub.txt";
    string algo = "IPH";
    if(argc>2){
        Model_file = argv[2];
    }
    if(argc>3){
        Data_file = argv[3];
    }
    vector<pair<double,double>> bounds(3, {numeric_limits<double>::max(),numeric_limits<double>::lowest()});
    if(argc>4){
        red_Model_file = argv[4];
    }
    if(argc>5){
        red_Data_file = argv[5];
    }
    rapidcsv::Document  red_Model_doc(red_Model_file, rapidcsv::LabelParams(0, -1),rapidcsv::SeparatorParams(' '));
    read_data(bounds, red_Model_doc, point_cloud_model, uav_model);
    rapidcsv::Document  red_Data_doc(red_Data_file, rapidcsv::LabelParams(0, -1),rapidcsv::SeparatorParams(' '));
    read_data(bounds, red_Data_doc, point_cloud_data, uav_data);
    rapidcsv::Document  Model_doc(Model_file, rapidcsv::LabelParams(0, -1),rapidcsv::SeparatorParams(' '));
    read_data(bounds, Model_doc, full_point_cloud_model, full_uav_model);
    rapidcsv::Document  Data_doc(Data_file, rapidcsv::LabelParams(0, -1),rapidcsv::SeparatorParams(' '));
    read_data(bounds, Data_doc, full_point_cloud_data, full_uav_data);
    
    
#ifdef USE_MATPLOT
    plot(point_cloud_model,point_cloud_data,0.4);
#endif
    
    
    
    
    double final_roll = 0, final_pitch = 0, final_yaw = 0;
    double total_time =0, time_start = 0, time_end = 0;
    double L2error_init = 0, L1error_init = 0;
    vector<int> matching(point_cloud_data.size());
    vector<double> err_per_point(point_cloud_data.size());
    L2error_init = computeL2error(point_cloud_model,point_cloud_data,matching,err_per_point);
    L1error_init = computeL1error(point_cloud_model,point_cloud_data,matching,err_per_point);
    DebugOn("Initial L2 error = " << L2error_init << endl);
    DebugOn("Initial L1 error = " << L1error_init << endl);
    time_start = get_wall_time();
    auto res = run_IPH(point_cloud_model, point_cloud_data, uav_model, uav_data);
    final_roll = get<0>(res);final_pitch = get<1>(res); final_yaw = get<2>(res);
#ifdef USE_MATPLOT
    plot(point_cloud_model,point_cloud_data,0.4);
#endif
    
    auto L2error = computeL2error(point_cloud_model,point_cloud_data,matching,err_per_point);
    auto L1error = computeL1error(point_cloud_model,point_cloud_data,matching,err_per_point);
    DebugOn("nm = " << point_cloud_model.size() << endl);
    DebugOn("nd = " << point_cloud_data.size() << endl);
    DebugOn("Initial L2 error = " << L2error_init << endl);
    DebugOn("Final L2 error = " << L2error << endl);
    DebugOn("Relative L2 gap = " << 100*(L2error_init-L2error)/L2error_init << "%" << endl);
    DebugOn("Initial L1 error = " << L1error_init << endl);
    DebugOn("Final L1 error = " << L1error << endl);
    DebugOn("Relative L1 gap = " << 100*(L1error_init-L1error)/L1error_init << "%" << endl);
    time_end = get_wall_time();
    DebugOn("Total wall clock time = " << time_end - time_start << endl);
    double shifted_x, shifted_y, shifted_z;
    auto tot_pts = full_point_cloud_model.size()+full_point_cloud_data.size();
    apply_rotation(final_roll, final_pitch, final_yaw, full_point_cloud_model, full_point_cloud_data, full_uav_model, full_uav_data);
    
    bool save_file = true;
    if(save_file){
        auto name = Model_file.substr(0,Model_file.find('.'));
        auto fname = name+"_ARMO_RPY_"+to_string(final_roll)+"_"+to_string(final_pitch)+"_"+to_string(final_yaw)+".laz";
        save_laz(fname,full_point_cloud_model, full_point_cloud_data);
        fname = name+"_ARMO_RPY_"+to_string(final_roll)+"_"+to_string(final_pitch)+"_"+to_string(final_yaw)+"_sub.laz";
        save_laz(fname,point_cloud_model, point_cloud_data);
    }
    DebugOn("Finished saving the LAZ files, please go to https://plas.io/ and upload the files to visualize the results!\n");
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

/* Scale point cloud using provided max values */
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


/* Return true if two cubes intersect with tolerance tol */
bool intersect(const vector<pair<double,double>>& a, const vector<pair<double,double>>& b, double tol) {
    return (a[0].first <= b[0].second + tol && a[0].second + tol >= b[0].first) &&
    (a[1].first <= b[1].second + tol && a[1].second + tol >= b[1].first) &&
    (a[2].first <= b[2].second + tol && a[2].second + tol >= b[2].first);
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


/* Return the max distance one can reach from point p given rotation and translation bounds */
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
tuple<double,double,double,double,double,double> run_MINLP(bool bypass, string axis, const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data){
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




/* Update the bounds on the rotation matrix entries based on the bounds of roll, pitch and yaw */
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

shared_ptr<Model<double>> build_norm2_SOC_MIQCP(vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, const indices& valid_cells, const indices& new_model_ids, const param<>& dist_cost, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double shift_min_x, double shift_max_x, double shift_min_y, double shift_max_y, double shift_min_z, double shift_max_z, vector<double>& rot_trans, bool convex, const vector<pair<pair<int,int>,pair<int,int>>>& incompatibles, param<>& norm_x, param<>& norm_y, param<>& norm_z,  param<>& intercept, const vector<int>& init_matching, const vector<double>& error_per_point, param<>& model_radius, bool relax_ints, bool relax_sdp, bool nonprop_scale, double perc_outliers){
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
    if(nonprop_scale){
        Reg->add(scale_x.in(R(1)),scale_y.in(R(1)),scale_z.in(R(1)));
        rot_trans.resize(15);
    }
    Reg->add_param(model_radius);
    Reg->add_param(x1);Reg->add_param(y1);Reg->add_param(z1);
    Reg->add_param(x2);Reg->add_param(y2);Reg->add_param(z2);
    var<> x_shift("x_shift", shift_min_x, shift_max_x), y_shift("y_shift", shift_min_y, shift_max_y), z_shift("z_shift", shift_min_z, shift_max_z);
    
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
    
    
    var<> theta11("theta11",  std::max(-1.,r11._range->first), std::min(1.,r11._range->second)), theta12("theta12", std::max(-1.,r12._range->first), std::min(1.,r12._range->second)), theta13("theta13", std::max(-1.,r13._range->first), std::min(1.,r13._range->second));
    var<> theta21("theta21", std::max(-1.,r21._range->first), std::min(1.,r21._range->second)), theta22("theta22", std::max(-1.,r22._range->first), std::min(1.,r22._range->second)), theta23("theta23", std::max(-1.,r23._range->first), std::min(1.,r23._range->second));
    var<> theta31("theta31", std::max(-1.,r31._range->first), std::min(1.,r31._range->second)), theta32("theta32", std::max(-1.,r32._range->first), std::min(1.,r32._range->second)), theta33("theta33", std::max(-1.,r33._range->first), std::min(1.,r33._range->second));
    
    
    Reg->add(theta11.in(R(1)),theta12.in(R(1)),theta13.in(R(1)));
    Reg->add(theta21.in(R(1)),theta22.in(R(1)),theta23.in(R(1)));
    Reg->add(theta31.in(R(1)),theta32.in(R(1)),theta33.in(R(1)));
    
    /* Bounds on the transformed data points */
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
    
    var<> new_x1("new_x1", x_new_lb, x_new_ub), new_y1("new_y1", y_new_lb, y_new_ub), new_z1("new_z1", z_new_lb, z_new_ub);
    if(nonprop_scale){
        /* Compute bounds on transformed data points using bound propagation rules */
        shared_ptr<pair<double,double>> x_range = make_shared<pair<double,double>>();
        shared_ptr<pair<double,double>> y_range = make_shared<pair<double,double>>();
        shared_ptr<pair<double,double>> z_range = make_shared<pair<double,double>>();
        for (int i = 0; i<nd; i++) {
            auto bounds = get_min_max(roll_min, roll_max, pitch_min, pitch_max, yaw_min, yaw_max, point_cloud_data[i], zeros);
            x_range->first = bounds[0].first + x_shift.get_lb().eval();
            x_range->second = bounds[0].second + x_shift.get_ub().eval();
            auto new_range = get_product_range(x_range, scale_x._range);
            x_new_lb.set_val(i, new_range->first);
            x_new_ub.set_val(i, new_range->second);
            y_range->first = bounds[1].first + y_shift.get_lb().eval();
            y_range->second = bounds[1].second + y_shift.get_ub().eval();
            new_range = get_product_range(y_range, scale_y._range);
            y_new_lb.set_val(i, new_range->first);
            y_new_ub.set_val(i, new_range->second);
            z_range->first = bounds[2].first + z_shift.get_lb().eval();
            z_range->second = bounds[2].second + z_shift.get_ub().eval();
            new_range = get_product_range(z_range, scale_z._range);
            z_new_lb.set_val(i, new_range->first);
            z_new_ub.set_val(i, new_range->second);
        }
        Reg->add(new_x1.in(N1), new_y1.in(N1), new_z1.in(N1));
    }
    var<> new_xm("new_xm", xm_min, xm_max), new_ym("new_ym", ym_min, ym_max), new_zm("new_zm", zm_min, zm_max);
    Reg->add(new_xm.in(N1), new_ym.in(N1), new_zm.in(N1));
    
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
    
    
    /* Add voronoi cuts */
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
    /* Add spatial branching variables and constraints */
    bool spatial_branching = true;
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
    
    /* Constraints defining the new model points matched with each data point */
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
    }
    
    /* Incompatible pairs constraints */
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
    }
    
    
    
    var<> x_diff("x_diff", pos_), y_diff("y_diff", pos_), z_diff("z_diff", pos_);
    Reg->add(x_diff.in(N1), y_diff.in(N1), z_diff.in(N1));
    
    auto ids1 = theta11.repeat_id(N1.size());
    if(nonprop_scale){
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
    if(!nonprop_scale)
        x_norm += pow((x1.in(N1)*theta11.in(ids1) + y1.in(N1)*theta12.in(ids1) + z1.in(N1)*theta13.in(ids1) + x_shift.in(ids1)) - new_xm, 2) - x_diff;
    else
        x_norm += pow((new_x1.in(N1) + x_shift.in(ids1)) - new_xm, 2) - x_diff;
    Reg->add(x_norm.in(N1)<=0);
    
    
    Constraint<> y_norm("y_norm");
    if(!nonprop_scale)
        y_norm += pow((x1.in(N1)*theta21.in(ids1) + y1.in(N1)*theta22.in(ids1) + z1.in(N1)*theta23.in(ids1) + y_shift.in(ids1)) - new_ym, 2) - y_diff;
    else
        y_norm += pow((new_y1.in(N1) + y_shift.in(ids1)) - new_ym, 2) - y_diff;
    Reg->add(y_norm.in(N1)<=0);
    
    Constraint<> z_norm("z_norm");
    if(!nonprop_scale)
        z_norm += pow((x1.in(N1)*theta31.in(ids1) + y1.in(N1)*theta32.in(ids1) + z1.in(N1)*theta33.in(ids1) + z_shift.in(ids1)) - new_zm, 2) - z_diff;
    else
        z_norm += pow((new_z1.in(N1) + z_shift.in(ids1)) - new_zm, 2) - z_diff;
    Reg->add(z_norm.in(N1)<=0);
    
    
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
            Reg->add(col1==1);
            Constraint<> col2("col2");
            col2 = pow(theta12,2)+pow(theta22,2)+pow(theta32,2);
            Reg->add(col2==1);
            Constraint<> col3("col3");
            col3 = pow(theta13,2)+pow(theta23,2)+pow(theta33,2);
            Reg->add(col3==1);
        }
        
    }

    /* Objective function */
    
    if(perc_outliers==0){
        Reg->min(sum(x_diff + y_diff + z_diff));
    }
    else {
        param<> ones("ones");
        ones.in(N1);
        ones = 1;
        Reg->min(ones.tr()*(keep*x_diff) + ones.tr()*(keep*y_diff) + ones.tr()*(keep*z_diff));
    }
    Reg->_status = 0;
    return(Reg);
}

shared_ptr<Model<double>> build_norm1_SOC_MIQCP(vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, const indices& valid_cells, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double shift_min_x, double shift_max_x, double shift_min_y, double shift_max_y, double shift_min_z, double shift_max_z, vector<double>& rot_trans, bool convex, const vector<pair<pair<int,int>,pair<int,int>>>& incompatibles, param<>& norm_x, param<>& norm_y, param<>& norm_z,  param<>& intercept, const vector<int>& init_matching, const vector<double>& error_per_point, bool relax_ints){
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
        
//        auto go_icp_max_dist = get_GoICP_dist(roll_max, shift_max_x, point_cloud_data[i], false);
            //        DebugOn("SDP max dist = " << to_string_with_precision(max_sdp, 6) << endl);
            //        auto max_dist = get_max_dist(-new_roll_max, new_roll_max, -new_pitch_max, new_pitch_max, -new_yaw_max, new_yaw_max, shift_min_x, shift_max_x, shift_min_y, shift_max_y, shift_min_z, shift_max_z, point_cloud_data[i], zeros, false);
        
            //        lower_bound += std::max(0.,error_per_point[i] - max_dist);
//        go_icp_lb += std::max(0.,error_per_point[i] - go_icp_max_dist);
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



/* Retrieve the value of the rotation and translation variables from model M, then compute roll, pitch and yaw angles and return true if the rotation matrix is valid (unitary)*/
bool get_angles(const shared_ptr<Model<double>>& M, vector<double>& angles, vector<int>& new_matching){
    auto theta11 = M->get_var<double>("theta11");auto theta12 = M->get_var<double>("theta12");auto theta13 = M->get_var<double>("theta13");
    auto theta21 = M->get_var<double>("theta21");auto theta22 = M->get_var<double>("theta22");auto theta23 = M->get_var<double>("theta23");
    auto theta31 = M->get_var<double>("theta31");auto theta32 = M->get_var<double>("theta32");auto theta33 = M->get_var<double>("theta33");
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
    
    
    
    auto pitch_val = std::atan2(theta32.eval(), theta33.eval())*180/pi;
    auto roll_val = std::atan2(-1*theta31.eval(), std::sqrt(theta32.eval()*theta32.eval()+theta33.eval()*theta33.eval()))*180/pi;
    auto yaw_val = std::atan2(theta21.eval(),theta11.eval())*180/pi;
    DebugOn("Roll (degrees) = " << to_string_with_precision(roll_val,12) << endl);
    DebugOn("Pitch (degrees) = " << to_string_with_precision(pitch_val,12) << endl);
    DebugOn("Yaw (degrees) = " << to_string_with_precision(yaw_val,12) << endl);
    angles[0] = roll_val;
    angles[1] = pitch_val;
    angles[2] = yaw_val;
    
    if(!is_rotation){
        DebugOn("WARNING, returned matrix is not a Rotation!\n");
    }
    return is_rotation;
}


/* Retrieve the value of the rotation, translation and scale variables from model M, return true if the rotation matrix is valid (unitary)*/
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
    if(rot_trans.size()>9){
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


vector<double> run_MISDP(vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, const vector<vector<double>>& uav_model, const vector<vector<double>>& uav_data)
{
    double angle_max = 0.1;
    int nb_pairs = 0, min_nb_pairs = numeric_limits<int>::max(), max_nb_pairs = 0, av_nb_pairs = 0;
    vector<pair<double,double>> min_max1;
    vector<vector<pair<double,double>>> min_max2(point_cloud_model.size());
    vector<int> nb_neighbors(point_cloud_data.size());
    vector<map<double,int>> neighbors;
    /* Compute cube for all points in point cloud 2 */
    for (auto i = 0; i<point_cloud_model.size(); i++) {
        min_max2[i] = get_min_max(angle_max, point_cloud_model.at(i), uav_model.at(i));
    }
    double roll_1 = 0, yaw_1 = 0, pitch_1 = 0;
    double dist_sq = 0;
    /* Check if cubes intersect */
    neighbors.resize(point_cloud_data.size());
    bool filter_neighbors = true;
    for (auto i = 0; i<point_cloud_data.size(); i++) {
        nb_pairs = 0;
        min_max1 = get_min_max(angle_max, point_cloud_data.at(i), uav_data.at(i));
        DebugOff("For point (" << point_cloud_data.at(i).at(0) << "," <<  point_cloud_data.at(i).at(1) << "," << point_cloud_data.at(i).at(2) << "): ");
        DebugOff("\n neighbors in umbrella : \n");
        for (size_t j = 0; j < point_cloud_model.size(); j++){
            if(filter_neighbors && intersect(min_max1, min_max2[j])){ /* point is in umbrella */
                nb_pairs++;
                dist_sq = std::pow(point_cloud_data.at(i)[0] - point_cloud_model.at(j)[0],2) + std::pow(point_cloud_data.at(i)[1] - point_cloud_model.at(j)[1],2) + std::pow(point_cloud_data.at(i)[2] - point_cloud_model.at(j)[2],2);
                neighbors[i].insert({dist_sq,j});/* TODO if many neighbors with exact same distance */
                DebugOff("(" << point_cloud_model.at(j).at(0) << "," <<  point_cloud_model.at(j).at(1) << "," << point_cloud_model.at(j).at(2) << ")\n");
            }
        }
        
        DebugOff("nb points in umbrella = " << nb_pairs << endl);
        if(nb_pairs>max_nb_pairs)
            max_nb_pairs = nb_pairs;
        if(nb_pairs<min_nb_pairs)
            min_nb_pairs = nb_pairs;
        av_nb_pairs += nb_pairs;
        
            //        std::cout << "For point (" << point_cloud_data.at(i).at(0) << "," <<  point_cloud_data.at(i).at(1) << "," << point_cloud_data.at(i).at(2) << ")"<< " knnSearch(n="<<m<<"): \n";
            //        for (size_t k = 0; k < m; k++)
            //            std::cout << "ret_index["<<k<<"]=" << ret_indexes[k] << " out_dist_sqr=" << out_dists_sqr[k] << " point = (" << point_cloud_model.at(ret_indexes[k]).at(0) << "," <<  point_cloud_model.at(ret_indexes[k]).at(1) << "," << point_cloud_model.at(ret_indexes[k]).at(2) << ")" << std::endl;
        nb_neighbors[i] = nb_pairs;
    }
    av_nb_pairs /= point_cloud_data.size();
    DebugOn("Min nb of Pairs = " << min_nb_pairs << endl);
    DebugOn("Max nb of Pairs = " << max_nb_pairs << endl);
    DebugOn("Average nb of Pairs = " << av_nb_pairs << endl);
    param<> x1("x1"), x2("x2"), y1("y1"), y2("y2"), z1("z1"), z2("z2");
    param<> x_uav1("x_uav1"), y_uav1("y_uav1"), z_uav1("z_uav1");
    param<> x_uav2("x_uav2"), y_uav2("y_uav2"), z_uav2("z_uav2");
        //        return 0;
    int m = 500;
        //            int m = 1;
    double xm_max = numeric_limits<double>::lowest(), ym_max = numeric_limits<double>::lowest(), zm_max = numeric_limits<double>::lowest();
    double xm_min = numeric_limits<double>::max(), ym_min = numeric_limits<double>::max(), zm_min = numeric_limits<double>::max();
    vector<double> min_dist(point_cloud_data.size(),numeric_limits<double>::max());
    vector<int> nearest(point_cloud_data.size());
    vector<string> nearest_id(point_cloud_data.size());
    string i_str, j_str;
    indices valid_cells("valid_cells");
    map<int,int> n2_map;
    int nb_max_neigh = m;
    map<int,map<double,int>> valid_cells_map;
    
    for (auto j = 0; j<point_cloud_model.size(); j++) {
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
        x_uav2.add_val(j_str,uav_model.at(j)[0]);
        y_uav2.add_val(j_str,uav_model.at(j)[1]);
        z_uav2.add_val(j_str,uav_model.at(j)[2]);
    }
    /* Keep points with neighbors >= m */
    for (auto i = 0; i<point_cloud_data.size(); i++) {
        nb_max_neigh = m;
            //        if(nb_neighbors[i]>=1){
        i_str = to_string(i+1);
        x_uav1.add_val(i_str,uav_data.at(i)[0]);
        x1.add_val(i_str,point_cloud_data.at(i)[0]);
        y_uav1.add_val(i_str,uav_data.at(i)[1]);
        y1.add_val(i_str,point_cloud_data.at(i)[1]);
        z_uav1.add_val(i_str,uav_data.at(i)[2]);
        z1.add_val(i_str,point_cloud_data.at(i)[2]);
        auto it = neighbors[i].begin();
        for (auto j = 0; j<std::min(m,(int)neighbors[i].size()); j++) {
            auto k = it->second;
            j_str = to_string(k+1);
            valid_cells_map[i].insert({it->first,k});
            if(min_dist[i]>dist_sq){
                min_dist[i] = dist_sq;
                nearest[i] = k;
                nearest_id[i] = j_str;
            }
            it++;
        }
            //        }
    }
    int nm = x2.get_dim();
    vector<int> face_vert;
    vector<double> v;
    vector<vector<vector<double>>> model_voronoi_normals(nm);/* Store the normal vector of each facet of the voronoi cell of each point */
    vector<vector<vector<double>>> model_voronoi_vertices(nm);/* Store the normal vector of each facet of the voronoi cell of each point */
    vector<vector<vector<double>>> model_face_pts(nm);/* Store a point from each facet of the voronoi cell of each point */
    vector<vector<double>> model_face_intercept(nm);/* Store the constant part (intercept) in the equation of the voronoi cell of each point */
    vector<double> model_voronoi_in_radius(nm);/* Store the radius of the largest ball contained IN the voronoi cell of each point */
    vector<double> model_voronoi_out_radius(nm);/* Store the radius of the smallest ball enclosing the voronoi cell of each point */
    param<> norm_x("norm_x"), norm_y("norm_y"), norm_z("norm_z"), intercept("intercept");
    param<> model_radius("model_radius");
    indices m_facets("m_facets");
    
    bool compute_voronoi = false;
    if(compute_voronoi){
#ifdef USE_VORO
        
        container model_con(xm_min,xm_max,ym_min,ym_max,zm_min,zm_max,10,10,10,false,false,false,8);
        for (int i = 0; i< nm; i++) { // Input iterator
            model_con.put(i, point_cloud_model[i][0], point_cloud_model[i][1], point_cloud_model[i][2]);
        }
        /* Compute the facets of the Voronoi cells of all model points */
        c_loop_all cl(model_con);
        int idx,nx,ny,nz;
        double x,y,z,x1,y1,z1;
        voronoicell c;
        vector<pair<double,int>> volume(nm);/*volume of each voronoi cell <volume,point_id>*/
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
                v_idx += 3;
                vertices[0]=v[v_idx];
                vertices[1]=v[v_idx+1];
                vertices[2]=v[v_idx+2];
                model_voronoi_vertices[idx][i]=vertices;
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
                    DebugOn("WARNING: model point cannot be on the voronoi face!\n");
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
        } while (cl.inc());
        DebugOn("The total number of faces = " << total_nb_faces << endl);
        
        int nb_reduced_model = std::min(nm,150);
        vector<vector<double>> red_point_cloud_model(nb_reduced_model);
        for (int i =0; i<nb_reduced_model; i++) {
            red_point_cloud_model[i] = point_cloud_model[volume[i].second];
        }
#endif
    }
    indices N1("N1"),N2("N2");
    
    int n1 = x1.get_dim();
    int n2 = x2.get_dim();
    DebugOn("n1 = " << n1 << endl);
    DebugOn("n2 = " << n2 << endl);
    
    indices ids = indices("in_x");
    ids = N2;
    ids.add_empty_row();
    int row_id = 0;
    for (const auto &vcel:valid_cells_map) {
        auto key_data=to_string(vcel.first+1);
        for (auto const model_id: vcel.second) {
            auto key=key_data+","+to_string(model_id.second+1);
            ids.add_in_row(row_id, to_string(model_id.second+1));
            valid_cells.insert(key);
        }
        row_id++;
    }
    
    N1 = range(1,n1);
    N2 = range(1,n2);
    indices M("M");
    M = range(1,m);
    vector<int> new_model_pts;
    indices new_model_ids;
    param<> dist_cost;
    double ub = 100;
    int nb_threads = 8;
        //    valid_cells = get_valid_pairs(point_cloud_data, point_cloud_model, -angle_max, angle_max, -angle_max, angle_max, -angle_max, angle_max, 0, 0, 0, 0, 0, 0, norm_x, norm_y, norm_z, intercept, model_voronoi_out_radius, new_model_pts, new_model_ids, true);
        //    valid_cells = preprocess_QP(point_cloud_data, point_cloud_model, valid_cells, -angle_max, angle_max, -angle_max, angle_max, -angle_max, angle_max, 0, 0, 0, 0, 0, 0, model_voronoi_normals, model_face_intercept, new_model_pts, new_model_ids, dist_cost, ub, nb_threads);
    DebugOn("Total size of valid cells = " << valid_cells.size() << endl);
    
    
    auto cells = valid_cells;
    
    string name="Norm2_MISDP";
    
    auto Reg=make_shared<Model<>>(name);
    
    Reg->add_param(x1);Reg->add_param(y1);Reg->add_param(z1);
    Reg->add_param(x2);Reg->add_param(y2);Reg->add_param(z2);
    var<int> bin("bin",0,1);
    Reg->add(bin.in(cells));
    DebugOff("Added " << cells.size() << " binary variables" << endl);
    double yaw_min = -5*pi/180., yaw_max = 5*pi/180., pitch_min =-5*pi/180.,pitch_max = 5*pi/180.,roll_min =-5*pi/180.,roll_max = 5*pi/180.;
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
        //    var<> new_x1("new_x1"), new_y1("new_y1"), new_z1("new_z1");
        //    var<> new_x2("new_x2"), new_y2("new_y2"), new_z2("new_z2");
        //
        //    Reg->add(new_x1.in(N1), new_y1.in(N1), new_z1.in(N1));
        //    Reg->add(new_x2.in(N2), new_y2.in(N2), new_z2.in(N2));
    
        //    indices ids = indices("id");
        //    ids = N2;
        //    ids.add_empty_row();
        //    for(auto i=0;i<n1;i++){
        //        for(auto j=1;j<=n2;j++){
        //            string model_key = to_string(j);
        //            if(cells.has_key(to_string(i+1)+","+model_key)){
        //                ids.add_in_row(i, model_key);
        //            }
        //        }
        //        if(ids._ids->at(i).size()==0){
        //            Reg->_status = -1;
        //            return {0,0,0};
        //        }
        //    }
    
    
        //    for(auto i=0;i<N1.size();i++){
        //        for(auto j=1;j<=N2.size();j++){
        //            string model_key = to_string(j);
        //            if(cells.has_key(to_string(i+1)+","+model_key)){
        //                ids.add_in_row(i, model_key);
        //            }
        //        }
        ////        if(ids._ids->at(i).size()==0){
        ////            Reg->_status = -1;
        ////            throw invalid_argument("infeasible model!");
        ////        }
        //    }
    
    var<> new_xm("new_xm", xm_min, xm_max), new_ym("new_ym", ym_min, ym_max), new_zm("new_zm", zm_min, zm_max);
    Reg->add(new_xm.in(N1), new_ym.in(N1), new_zm.in(N1));
    
    
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
    
    
    var<> x_diff("x_diff", pos_), y_diff("y_diff", pos_), z_diff("z_diff", pos_);
    Reg->add(x_diff.in(N1), y_diff.in(N1), z_diff.in(N1));
    
    auto ids1 = theta11.repeat_id(N1.size());
    Constraint<> x_abs1("x_abs1");
    x_abs1 += x_diff - (((x1 - x_uav1)*theta11.in(ids1) + (y1 - y_uav1)*theta12.in(ids1) + (z1 - z_uav1)*theta13.in(ids1) + x_uav1) - new_xm);
    Reg->add(x_abs1.in(N1)>=0);
    
    Constraint<> x_abs2("x_abs2");
    x_abs2 += x_diff - (new_xm - ((x1 - x_uav1)*theta11.in(ids1) + (y1 - y_uav1)*theta12.in(ids1) + (z1 - z_uav1)*theta13.in(ids1) + x_uav1));
    Reg->add(x_abs2.in(N1)>=0);
    
    Constraint<> y_abs1("y_abs1");
    y_abs1 += y_diff - (((x1 - x_uav1)*theta21.in(ids1) + (y1 - y_uav1)*theta22.in(ids1) + (z1 - z_uav1)*theta23.in(ids1) + y_uav1) - new_ym);
    Reg->add(y_abs1.in(N1)>=0);
    
    Constraint<> y_abs2("y_abs2");
    y_abs2 += y_diff - (new_ym - ((x1 - x_uav1)*theta21.in(ids1) + (y1 - y_uav1)*theta22.in(ids1) + (z1 - z_uav1)*theta23.in(ids1) + y_uav1));
    Reg->add(y_abs2.in(N1)>=0);
    
    Constraint<> z_abs1("z_abs1");
    z_abs1 += z_diff - (((x1 - x_uav1)*theta31.in(ids1) + (y1 - y_uav1)*theta32.in(ids1) + (z1 - z_uav1)*theta33.in(ids1) + z_uav1) - new_zm);
    Reg->add(z_abs1.in(N1)>=0);
    
    Constraint<> z_abs2("z_abs2");
    z_abs2 += z_diff - (new_zm - ((x1 - x_uav1)*theta31.in(ids1) + (y1 - y_uav1)*theta32.in(ids1) + (z1 - z_uav1)*theta33.in(ids1) + z_uav1));
    Reg->add(z_abs2.in(N1)>=0);
    
    
    
    bool relax_sdp = false;
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
            //        if(convex){
            //            Constraint<> row1("row1");
            //            row1 = pow(theta11,2)+pow(theta12,2)+pow(theta13,2);
            //            Reg->add(row1.in(range(0,0))<=1);
            //            Constraint<> row2("row2");
            //            row2 = pow(theta21,2)+pow(theta22,2)+pow(theta23,2);
            //            Reg->add(row2.in(range(0,0))<=1);
            //            Constraint<> row3("row3");
            //            row3 = pow(theta31,2)+pow(theta32,2)+pow(theta33,2);
            //            Reg->add(row3.in(range(0,0))<=1);
            //            Constraint<> col1("col1");
            //            col1 = pow(theta11,2)+pow(theta21,2)+pow(theta31,2);
            //            Reg->add(col1.in(range(0,0))<=1);
            //            Constraint<> col2("col2");
            //            col2 = pow(theta12,2)+pow(theta22,2)+pow(theta32,2);
            //            Reg->add(col2.in(range(0,0))<=1);
            //            Constraint<> col3("col3");
            //            col3 = pow(theta13,2)+pow(theta23,2)+pow(theta33,2);
            //            Reg->add(col3.in(range(0,0))<=1);
            //        }
            //        else {
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
            //        }
        
    }
    bool spatial_branching = true;
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
            
            var<int> sbin_roll("sbin_roll", 0, 1), sbin_pitch("sbin_pitch", 0, 1), sbin_yaw("sbin_yaw", 0, 1);
            sbin_roll._priority = 1e6;sbin_pitch._priority = 1e6;sbin_yaw._priority = 1e6;
            Reg->add(sbin_roll.in(spatial_ids),sbin_pitch.in(spatial_ids),sbin_yaw.in(spatial_ids));
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
            
            
                //    Reg->print();
            double angle_min = std::min(std::min(roll_min,pitch_min),yaw_min), angle_max = std::max(std::max(roll_max,pitch_max),yaw_max);
            double angle_increment = (angle_max - angle_min)/nb_pieces;/* Angles are defined in [-2,2] in radians (-120,120) in degrees */
            
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
            }
            Reg->add_param(angle_lb);
            Reg->add_param(angle_ub);
            Reg->add_param(cos_lb);
            Reg->add_param(cos_ub);
            Reg->add_param(sin_lb);
            Reg->add_param(sin_ub);
            
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
    Reg->min(sum(x_diff + y_diff + z_diff));
    Reg->write();
    solver<> S(Reg,gurobi);
    S.use_callback();
    S.run(5,1e-6,300,1000);
    
    Reg->print_solution();
    vector<double> rot(9);
    vector<int> matching(n1);
//    bool is_rotation = get_solution_rot(Reg, rot, matching);
    return rot;
}


tuple<double,double,double> run_NLP(string axis, const vector<vector<double>>& point_cloud1, const vector<vector<double>>& point_cloud2, const vector<vector<double>>& uav1, const vector<vector<double>>& uav2)
{
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

/* Update point cloud coordinates */
void update_xyz(vector<vector<double>>& point_cloud1, vector<double>& x_vec1, vector<double>& y_vec1, vector<double>& z_vec1){
    for (auto i = 0; i< point_cloud1.size(); i++) {
        point_cloud1[i][0] = x_vec1[i];
        point_cloud1[i][1] = y_vec1[i];
        point_cloud1[i][2] = z_vec1[i];
    }
}


void apply_rotation(double roll, double pitch, double yaw, vector<vector<double>>& point_cloud1, const vector<vector<double>>& uav1){
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
}

void apply_rotation_mat(const vector<double>& R, vector<vector<double>>& point_cloud1, const vector<vector<double>>& uav1){
    /* Apply rotation */
    double shifted_x, shifted_y, shifted_z;
    for (auto i = 0; i< point_cloud1.size(); i++) {
        shifted_x = point_cloud1[i][0] - uav1[i][0];
        shifted_y = point_cloud1[i][1] - uav1[i][1];
        shifted_z = point_cloud1[i][2] - uav1[i][2];
        point_cloud1[i][0] = shifted_x*R[0] + shifted_y*R[1] + shifted_z*R[2];
        point_cloud1[i][1] = shifted_x*R[3] + shifted_y*R[4] + shifted_z*R[5];
        point_cloud1[i][2] = shifted_x*R[6] + shifted_y*R[7] + shifted_z*R[8];
        point_cloud1[i][0] += uav1[i][0];
        point_cloud1[i][1] += uav1[i][1];
        point_cloud1[i][2] += uav1[i][2];
    }
}

void apply_rotation_mat(const vector<double>& R, vector<vector<double>>& point_cloud1, vector<vector<double>>& point_cloud2, const vector<vector<double>>& uav1, const vector<vector<double>>& uav2){
    /* Apply rotation */
    double shifted_x, shifted_y, shifted_z;
    for (auto i = 0; i< point_cloud1.size(); i++) {
        shifted_x = point_cloud1[i][0] - uav1[i][0];
        shifted_y = point_cloud1[i][1] - uav1[i][1];
        shifted_z = point_cloud1[i][2] - uav1[i][2];
        point_cloud1[i][0] = shifted_x*R[0] + shifted_y*R[1] + shifted_z*R[2];
        point_cloud1[i][1] = shifted_x*R[3] + shifted_y*R[4] + shifted_z*R[5];
        point_cloud1[i][2] = shifted_x*R[6] + shifted_y*R[7] + shifted_z*R[8];
        point_cloud1[i][0] += uav1[i][0];
        point_cloud1[i][1] += uav1[i][1];
        point_cloud1[i][2] += uav1[i][2];
    }
    for (auto i = 0; i< point_cloud2.size(); i++) {
        shifted_x = point_cloud2[i][0] - uav2[i][0];
        shifted_y = point_cloud2[i][1] - uav2[i][1];
        shifted_z = point_cloud2[i][2] - uav2[i][2];
        point_cloud2[i][0] = shifted_x*R[0] - shifted_y*R[1] - shifted_z*R[2];
        point_cloud2[i][1] = -1*shifted_x*R[3] + shifted_y*R[4] + shifted_z*R[5];
        point_cloud2[i][2] = -1*shifted_x*R[6] + shifted_y*R[7] + shifted_z*R[8];
        point_cloud2[i][0] += uav2[i][0];
        point_cloud2[i][1] += uav2[i][1];
        point_cloud2[i][2] += uav2[i][2];
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



/* Return vector of extreme points from point cloud with corresponding INS refernce points */
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

/* Make sure the point cloud is centered around the origin */
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
        DebugOff("error(" << i+1 << ") = " << to_string_with_precision(min_dist,12) << endl);
        DebugOff("matching(" << i+1 << ") = " << matching[i]+1 << endl);
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
    
    plt::plot3(x_vec_model, y_vec_model, z_vec_model,x_vec_data, y_vec_data, z_vec_data, keywords);
    
    plt::show();
}
#endif

/* Run Go-ICP on point clouds */
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

tuple<double,double,double> run_IPH(vector<vector<double>>& ext_model, vector<vector<double>>& ext_data, const vector<vector<double>>& uav1, const vector<vector<double>>& uav2){
    double roll = 0, pitch = 0, yaw = 1; /* at least one nonzero to enter the while loop */
    double final_roll = 0, final_pitch = 0, final_yaw = 0;
    int nb_iter = 0, max_nb_iter = 100;
    tuple<double,double,double> res;
    while(nb_iter < max_nb_iter && std::abs(roll)+std::abs(pitch)+std::abs(yaw)>1e-1) {
        res = run_NLP("full", ext_model, ext_data, uav1, uav2);
        roll = get<0>(res);pitch = get<1>(res);yaw = get<2>(res);
        final_roll += roll;final_pitch += pitch;final_yaw += yaw;
        apply_rotation(roll, pitch, yaw, ext_model, ext_data, uav1, uav2);
        nb_iter++;
        DebugOn("No projection, ITERATION " << nb_iter << endl);
    }
    nb_iter = 0;yaw=1;max_nb_iter = 100;
    while(nb_iter < max_nb_iter && std::abs(roll)+std::abs(pitch)+std::abs(yaw)>1e-1) {
        res = run_NLP("z", ext_model, ext_data, uav1, uav2);
        roll = get<0>(res);pitch = get<1>(res);yaw = get<2>(res);
        final_roll += roll;final_pitch += pitch;final_yaw += yaw;
        apply_rotation(roll, pitch, yaw, ext_model, ext_data, uav1, uav2);
        nb_iter++;
        DebugOn("Projceting out z axis, ITERATION " << nb_iter << endl);
    }
    nb_iter = 0;yaw=1;max_nb_iter = 100;
    while(nb_iter < max_nb_iter && std::abs(roll)+std::abs(pitch)+std::abs(yaw)>1e-1) {
        res = run_NLP("y", ext_model, ext_data, uav1, uav2);
        roll = get<0>(res);pitch = get<1>(res);yaw = get<2>(res);
        final_roll += roll;final_pitch += pitch;final_yaw += yaw;
        apply_rotation(roll, pitch, yaw, ext_model, ext_data, uav1, uav2);
        nb_iter++;
        DebugOn("Projceting out y axis, ITERATION " << nb_iter << endl);
    }
    nb_iter = 0;yaw=1;max_nb_iter = 100;
    while(nb_iter < max_nb_iter && std::abs(roll)+std::abs(pitch)+std::abs(yaw)>1e-1) {
        res = run_NLP("x", ext_model, ext_data, uav1, uav2);
        roll = get<0>(res);pitch = get<1>(res);yaw = get<2>(res);
        final_roll += roll;final_pitch += pitch;final_yaw += yaw;
        apply_rotation(roll, pitch, yaw, ext_model, ext_data, uav1, uav2);
        nb_iter++;
        DebugOn("Projceting out x axis, ITERATION " << nb_iter << endl);
    }
    nb_iter = 0;yaw=1;max_nb_iter = 100;
    while(nb_iter < max_nb_iter && std::abs(roll)+std::abs(pitch)+std::abs(yaw)>1e-1) {
        res = run_NLP("full", ext_model, ext_data, uav1, uav2);
        roll = get<0>(res);pitch = get<1>(res);yaw = get<2>(res);
        final_roll += roll;final_pitch += pitch;final_yaw += yaw;
        apply_rotation(roll, pitch, yaw, ext_model, ext_data, uav1, uav2);
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






/* Read input files */
void read_data(vector<pair<double,double>>& bounds, const rapidcsv::Document& Model_doc,vector<vector<double>>& point_cloud, vector<vector<double>>& uav, bool only_keep_nadir){
    int model_nb_rows = Model_doc.GetRowCount();
    if(model_nb_rows<3){
        throw invalid_argument("Input file with less than 2 points");
    }
    DebugOn("Input file has " << model_nb_rows << " rows" << endl);
    if(!only_keep_nadir){
        point_cloud.resize(model_nb_rows);
        uav.resize(model_nb_rows);
    }
    for (int i = 0; i< model_nb_rows; i++) {
        auto laser_id = Model_doc.GetCell<int>(0, i);
        if(only_keep_nadir && laser_id!=15){/* Only keep points from Nadir laser */
            continue;
        }
        auto x = Model_doc.GetCell<double>(1, i);
        auto y = Model_doc.GetCell<double>(2, i);
        auto z = Model_doc.GetCell<double>(3, i);
        if(bounds[0].first > x)
            bounds[0].first = x;
        if(bounds[0].second < x)
            bounds[0].second = x;
        if(bounds[1].first > y)
            bounds[1].first = y;
        if(bounds[1].second < y)
            bounds[1].second = y;
        if(bounds[2].first > z)
            bounds[2].first = z;
        if(bounds[2].second < z)
            bounds[2].second = z;
        auto uav_x = Model_doc.GetCell<double>(4, i);
        auto uav_y = Model_doc.GetCell<double>(5, i);
        auto uav_z = Model_doc.GetCell<double>(6, i);
        if(only_keep_nadir){
            point_cloud.push_back(vector<double>(3));
            point_cloud.back()[0] = x;
            point_cloud.back()[1] = y;
            point_cloud.back()[2] = z;
            uav.push_back(vector<double>(3));
            uav.back()[0] = uav_x;
            uav.back()[1] = uav_y;
            uav.back()[2] = uav_z;
        }
        else {
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
}



indices preprocess_QP(const vector<vector<double>>& point_cloud_data, const vector<vector<double>>& point_cloud_model, const indices& old_cells, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double shift_min_x, double shift_max_x, double shift_min_y, double shift_max_y, double shift_min_z, double shift_max_z, const vector<vector<vector<double>>>& model_voronoi_normals, const vector<vector<double>>& model_face_intercept, vector<int>& new_model_pts, indices& new_model_ids, param<>& dist_cost, double upper_bound, int nb_total_threads, double scale){
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
    vector<double> dist_i;
    vector<double> zeros = {0,0,0};
    vector<shared_ptr<Model<double>>> batch_models;
    int missed=0;
    vector<double> x_lb, x_ub, y_lb, y_ub, z_lb, z_ub;
    for(auto i=0;i<nd;i++){
        auto max_dist_i = get_max_dist(roll_min, roll_max, pitch_min, pitch_max, yaw_min, yaw_max, shift_min_x, shift_max_x, shift_min_y, shift_max_y, shift_min_z, shift_max_z, point_cloud_data[i], zeros, false);
        dist_i.push_back(max_dist_i*scale);
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
        for(auto k=0;k<intercepts.size();k++){
            a.set_val(k, normals[k][0]);
            b.set_val(k, normals[k][1]);
            c.set_val(k, normals[k][2]);
            d.add_val(k, intercepts[k]);
        }
        Constraint<> region("region");
        region=a*x+b*y+c*z+d;
        voro_model->add(region.in(voronoi_ids)<=0);
        
            //        Constraint<> radius("radius");
            //        radius=(x-x_data)*(x-x_data)+(y-y_data)*(y-y_data)+(z-z_data)*(z-z_data)-min_dist_sq;
            //        voro_model->add(radius.in(range(0,0))<=0);
        
            //        func<> obj=(x-x_model)*(x-x_model)+(y-y_model)*(y-y_model)+(z-z_model)*(z-z_model);
        func<> obj=(x-x_data)*(x-x_data)+(y-y_data)*(y-y_data)+(z-z_data)*(z-z_data);
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
            min_dist_sq=dist_i[i];
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
        if(status_k==0 && (sqrt(obj_k) - dist_i[ik_int])<=1e-6){
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
