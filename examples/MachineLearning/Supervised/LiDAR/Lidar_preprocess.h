//
//  Lidar_preprocess.h
//  Gravity
//
//  Created by Smitha on 8/24/21.
//

#ifndef Lidar_preprocess_h
#define Lidar_preprocess_h
#include "Branch_Bound.h"
bool vertices_box_plane(const vector<double>& plane_eq,  const vector<vector<double>>& big_box, vector<vector<double>>& new_verts, vector<int>& infeas_set);
void get_extreme_point_data(vector<vector<double>>& extreme, const vector<double>& uav_d, const vector<double>& d_pt, const vector<var<double>>& theta_vec);
void get_extreme_point_model(vector<vector<double>>& extreme, const vector<double>& uav_d, const vector<double>& uav_m, const vector<double>& m_pt, const vector<var<double>>& theta_vec);
double max_distance_polytopes(const vector<vector<double>>& poly1,  const vector<vector<double>>& poly2);
/* New valid cells and minimum distance each cell
 @point_cloud_model: model lidar point cloud
 @point_cloud_data: data lidar point cloud
 @uav_model: model uav point cloud
 @uav_data: data uav point cloud
 @valid_cells_old: previous estimate of valid cells
 @new_cells: valid cells identified by this algo
 @dist_cells: minimum distance min_dist_ij for each cell
 @roll_min, roll_max, pitch_min, pitch_max, yaw_min, yaw_max: Bounds on rotation angles
 
 ALL DISTANCES ARE SQUARED UNLESS _ROOT
 */
void preprocess_lid(vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data, const vector<vector<double>>& uav_model, const vector<vector<double>>& uav_data, indices& valid_cells_old, indices& new_cells, param<double>& dist_cells, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double upper_bound)
{
    
    indices valid_cells("valid_cells");
    indices valid_cells_new("valid_cells_new");
    indices valid_cells_empty("valid_cells_empty");
    param<double> dist_cells_old ("dist_cells_old");
    double time_start = get_wall_time();
    int new_test=0;
    planeStatPerPair = 0;
   
    /*variables to define rotation on the positive side (model side)*/

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
    
    theta11.print();
    theta12.print();
    theta13.print();
    theta21.print();
    theta22.print();
    theta23.print();
    theta31.print();
    theta32.print();
    theta33.print();
    
    /*variables to define rotation on the reverse side (data side)*/
    
    func<> r11r = cos(yaw*(-1))*cos(roll*(-1));r11r.eval_all();
    func<> r12r = cos(yaw*(-1))*sin(roll*(-1))*sin(pitch) - sin(yaw*(-1))*cos(pitch);r12r.eval_all();
    func<> r13r = cos(yaw*(-1))*sin(roll*(-1))*cos(pitch) + sin(yaw*(-1))*sin(pitch);r13r.eval_all();
    func<> r21r = sin(yaw*(-1))*cos(roll*(-1));r21r.eval_all();
    func<> r22r = sin(yaw*(-1))*sin(roll*(-1))*sin(pitch) + cos(yaw*(-1))*cos(pitch);r22r.eval_all();
    func<> r23r = sin(yaw*(-1))*sin(roll*(-1))*cos(pitch) - cos(yaw*(-1))*sin(pitch);r23r.eval_all();
    func<> r31r = sin(-1*roll*(-1));r31r.eval_all();
    func<> r32r = cos(roll*(-1))*sin(pitch);r32r.eval_all();
    func<> r33r = cos(roll*(-1))*cos(pitch);r33r.eval_all();
    
   
    

    var<> theta11r("theta11r",  std::max(-1.,r11r._range->first), std::min(1.,r11r._range->second)), theta12r("theta12r", std::max(-1.,r12r._range->first), std::min(1.,r12r._range->second)), theta13r("theta13r", std::max(-1.,r13r._range->first), std::min(1.,r13r._range->second));
    var<> theta21r("theta21r", std::max(-1.,r21r._range->first), std::min(1.,r21r._range->second)), theta22r("theta22r", std::max(-1.,r22r._range->first), std::min(1.,r22r._range->second)), theta23r("theta23r", std::max(-1.,r23r._range->first), std::min(1.,r23r._range->second));
    var<> theta31r("theta31r", std::max(-1.,r31r._range->first), std::min(1.,r31r._range->second)), theta32r("theta32r", std::max(-1.,r32r._range->first), std::min(1.,r32r._range->second)), theta33r("theta33r", std::max(-1.,r33r._range->first), std::min(1.,r33r._range->second));
    
    vector<var<double>> T2;
    T2.push_back(theta11r);
    T2.push_back(theta12r);
    T2.push_back(theta13r);
    T2.push_back(theta21r);
    T2.push_back(theta22r);
    T2.push_back(theta23r);
    T2.push_back(theta31r);
    T2.push_back(theta32r);
    T2.push_back(theta33r);
    
    theta11r.print();
    theta12r.print();
    theta13r.print();
    theta21r.print();
    theta22r.print();
    theta23r.print();
    theta31r.print();
    theta32r.print();
    theta33r.print();
    
    double min_cost_sum=0.0;
    
    
    bool found_all=true;

    int nd=point_cloud_data.size();
    int nm=point_cloud_model.size();

    vector<map<double, int>> valid_cells_map(nd);
    
    map<int, bool> new_model_pts;

   
    for(auto i=0;i<nd;i++){
        double min_dist_ij_max=numeric_limits<double>::max();
        double min_dist_ij_min=numeric_limits<double>::max();
        vector<vector<double>> extreme_i;
        /*Feasible region R(data-uav_data)*/
        get_extreme_point_data(extreme_i, uav_data[i], point_cloud_data[i], T2);
        for (int j = 0; j<nm; j++) {
            string key= to_string(i+1)+","+to_string(j+1);
            if(!valid_cells_old.has_key(key)){
                DebugOff("continued");
                continue;
            }
            double dist_ij_min, dist_ij_max;
            vector<vector<double>> extreme_j;
            /*Feasible region R(model-uav_model)+(uav_model-uav_data)*/
            get_extreme_point_model(extreme_j, uav_data[i], uav_model[j], point_cloud_model[j], T1);
            /*Calling GJK*/
            dist_ij_min=std::max(distance_polytopes_gjk(extreme_i, extreme_j)-1e-6, 0.0);
            if(dist_ij_min<=upper_bound && dist_ij_min<=min_dist_ij_max){
                dist_ij_max=max_distance_polytopes(extreme_i, extreme_j);
                if(dist_ij_max<=min_dist_ij_max){
                    min_dist_ij_max=dist_ij_max;
                }
                if(dist_ij_min<=min_dist_ij_min){
                    min_dist_ij_min=dist_ij_min;
                }
                valid_cells_map[i].insert(pair<double, int>(dist_ij_max, j));
                dist_cells_old.add_val(key, dist_ij_min);
                valid_cells.insert(key);
            }
        }
        min_cost_sum+=min_dist_ij_min;
        if(valid_cells_map[i].empty() || min_cost_sum>=upper_bound+1e-6){
            found_all=false;
            break;
        }
    }
            
    /*Looping again to ensure all valid cells have min_dist less than min_dist_ij_max*/
    if(found_all){
        for(auto i=0;i<nd;i++){
            auto it=valid_cells_map[i].begin();
            double min_dist_ij_max=it->first;
            for (int j = 0; j<nm; j++) {
                string key= to_string(i+1)+","+to_string(j+1);
                if(!valid_cells.has_key(key)){
                    DebugOff("continued");
                    continue;
                }
                auto dij=dist_cells_old.eval(key);
                if(dij<=min_dist_ij_max){
                    valid_cells_new.insert(key);
                    dist_cells.add_val(key, dij);
                }
            }
        }
        DebugOn("min_cost_sum "<<min_cost_sum<<endl);
        double vo=valid_cells_old.size();
        double vn=valid_cells_new.size();
        double remo=(vo-vn)/vo*100.0;
        DebugOn("valid cells old size "<<vo<<endl);
        DebugOn("valid cells old size "<<vn<<endl);
        DebugOn("rem percen "<<remo<<endl);
        new_cells=valid_cells_new;
    }
    else{
        new_cells=valid_cells_empty;
        DebugOn("No valid cells found "<<endl);
    }
}
/* Extreme points of the feasible region in which a model point can lie for alignment
 @extreme: Vertices of extreme points of the feasible region: [R'[d_pt-uav_d]]
 @uav_d: uav coordinates of data point
 @d_pt: Data point lidar coordinate
 @theta_vec: vector of variables defining rotation with +roll and +yaw
 SECANT reverse trainge inequality
 Tangent Traingle inequality
 */
void get_extreme_point_data(vector<vector<double>>& extreme, const vector<double>& uav_d, const vector<double>& d_pt, const vector<var<double>>& theta_vec){
    vector<double> d(3);
    
    vector<vector<double>> new_vert_a, new_vert_b;
    vector<int> infeas_set_a, infeas_set_b, infeas_set;
    bool vertex_found_a=false, vertex_found_b=false;
    
    d[0]=d_pt[0]-uav_d[0];
    d[1]=d_pt[1]-uav_d[1];
    d[2]=d_pt[2]-uav_d[2];
    
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
    
    vertex_found_a=vertices_box_plane(secant, box_big, new_vert_a, infeas_set_a);
    
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
        zt=(z_lb+z_ub)*0.5;
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
        vertex_found_b=vertices_box_plane(tangent, box_big, new_vert_b, infeas_set_b);
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
    if(vertex_found_b){
        for(auto v: new_vert_b){
            extreme.push_back(v);
        }
        for(auto i: infeas_set_b){
            infeas_set.push_back(i);
        }
    }
    for(auto i=0;i<box_big.size();i++){
        if(std::find (infeas_set.begin(), infeas_set.end(), i) ==infeas_set.end()){
            extreme.push_back(box_big[i]);
        }
    }

}
/* Extreme points of the feasible region in which a model point can lie for alignment
 @extreme: Vertices of extreme points of the feasible region: [R[m-uav_m]+uav_m-uav_d]
 @uav_d: uav coordinates of data point
 @uav_m: uav coordinates of model point
 @m_pt: Model point lidar coordinate
 @theta_vec: vector of variables defining rotation with -roll and -yaw
 SECANT reverse trainge inequality
 Tangent Traingle inequality
 */
void get_extreme_point_model(vector<vector<double>>& extreme, const vector<double>& uav_d, const vector<double>& uav_m, const vector<double>& m_pt, const vector<var<double>>& theta_vec)
{
    vector<double> d(3), b(3);
    
    vector<vector<double>> new_vert_a, new_vert_b;
    vector<int> infeas_set_a, infeas_set_b, infeas_set;
    bool vertex_found_a=false, vertex_found_b=false;
    
    d[0]=m_pt[0]-uav_m[0];
    d[1]=m_pt[1]-uav_m[1];
    d[2]=m_pt[2]-uav_m[2];
    
    b[0]=uav_m[0]-uav_d[0];
    b[1]=uav_m[1]-uav_d[1];
    b[2]=uav_m[2]-uav_d[2];
    
    double d_mag=pow(d[0],2)+pow(d[1],2)+pow(d[2],2);
    double b_mag=pow(b[0],2)+pow(b[1],2)+pow(b[2],2);
    double d_root=sqrt(d_mag);
    double b_root=sqrt(b_mag);
    
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
    
    rot_x1_bounds->first=std::max(x_range->first + y_range->first + z_range->first, d_root*(-1))+b[0];
    rot_x1_bounds->second=std::min(x_range->second + y_range->second + z_range->second, d_root)+b[0];
    
    double x_lb=rot_x1_bounds->first;
    double x_ub=rot_x1_bounds->second;
    
    x_range  = get_product_range(x1_bounds, theta21._range);
    y_range  = get_product_range(y1_bounds, theta22._range);
    z_range  = get_product_range(z1_bounds, theta23._range);
    
    rot_y1_bounds->first=std::max(x_range->first + y_range->first + z_range->first, d_root*(-1))+b[1];
    rot_y1_bounds->second=std::min(x_range->second + y_range->second + z_range->second, d_root)+b[1];
    
    double y_lb=rot_y1_bounds->first;
    double y_ub=rot_y1_bounds->second;
    
    
    x_range  = get_product_range(x1_bounds, theta31._range);
    y_range  = get_product_range(y1_bounds, theta32._range);
    z_range  = get_product_range(z1_bounds, theta33._range);
    
    rot_z1_bounds->first=std::max(x_range->first + y_range->first + z_range->first, d_root*(-1))+b[2];
    rot_z1_bounds->second=std::min(x_range->second + y_range->second + z_range->second, d_root)+b[2];
    
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
    
    vector<double> secant(4,0.0);
    secant[0]=(x_lb+x_ub)*(-1);
    secant[1]=(y_lb+y_ub)*(-1);
    secant[2]=(z_lb+z_ub)*(-1);
    secant[3]=pow(d_root-b_root,2)+x_lb*x_ub+y_lb*y_ub+z_lb*z_ub;
    
    vertex_found_a=vertices_box_plane(secant, box_big, new_vert_a, infeas_set_a);
    

    
    /*Find coordinates of the point to define tangent, Tangent is found parallel to secant, if this fails midpoint is used and checked to ensure does not intersect secant in this range*/
    double xt, yt, zt;
    bool plane_eq_found=false;
    vector<double> tangent(4);
    double temp_x=sqrt(pow(d_root+b_root,2)/(1+((y_lb+y_ub)*(y_lb+y_ub)+(z_lb+z_ub)*(z_lb+z_ub))/((x_lb+x_ub)*(x_lb+x_ub))));

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
        zt=(z_lb+z_ub)*0.5;
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
        DebugOff("Plane Eq found in extreme model "<<endl);
        vertex_found_b=vertices_box_plane(tangent, box_big, new_vert_b, infeas_set_b);
    }
    else{
        DebugOn("Plane Eq no found in extreme model "<<plane_eq_found<<endl);
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
bool vertices_box_plane(const vector<double>& plane_eq, const vector<vector<double>>& big_box, vector<vector<double>>& new_verts, vector<int>& infeas_set){
    
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
/* Distance between two polytopes
 @poly1: Vertices of polytope1
 @poly2: Vertices of polytope2
 */
double max_distance_polytopes(const vector<vector<double>>& poly1,  const vector<vector<double>>& poly2){
    double dmax=0;
    for(auto v1: poly1){
        for(auto v2: poly2){
            auto d=pow(v1[0]-v2[0],2)+pow(v1[1]-v2[1],2)+pow(v1[2]-v2[2],2);
            if(d>=dmax){
                dmax=d;
            }
        }
    }
    return dmax;
}
#endif /* Lidar_preprocess_h */

