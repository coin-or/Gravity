//
//  Lidar_utils.hpp
//  lidar
//
//  Created by Smitha on 8/18/21.
//

#ifndef Lidar_utils_hpp
#define Lidar_utils_hpp

#include <stdio.h>
#include <DataSet.h> 
using namespace std;
#include <gravity/jly_goicp.h>
#include <gravity/ConfigMap.hpp>
using namespace Go_ICP;
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
void get_angle_rotation_transl_matrix(const shared_ptr<Model<double>>& M, vector<double>& rot_trans);
void round_bin(shared_ptr<Model<double>>& M, int nd, int nm);
void update_matching(shared_ptr<Model<double>>& M, vector<int>& new_matching);
void round_bin(shared_ptr<Model<double>>& M, int nd, int nm)
{
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
void update_matching(shared_ptr<Model<double>>& M, vector<int>& new_matching)
{
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
void get_angle_rotation_transl_matrix(const shared_ptr<Model<double>>& M, vector<double>& rot_trans)
{
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
bool get_solution(const shared_ptr<Model<double>>& M, vector<double>& rot_trans, vector<int>& new_matching)
{
    auto theta11 = M->get_var<double>("theta11");auto theta12 = M->get_var<double>("theta12");auto theta13 = M->get_var<double>("theta13");
    auto theta21 = M->get_var<double>("theta21");auto theta22 = M->get_var<double>("theta22");auto theta23 = M->get_var<double>("theta23");
    auto theta31 = M->get_var<double>("theta31");auto theta32 = M->get_var<double>("theta32");auto theta33 = M->get_var<double>("theta33");
    var<> x_shift, y_shift, z_shift;
    if(rot_trans.size()>9){
        x_shift = M->get_var<double>("x_shift");y_shift = M->get_var<double>("y_shift");z_shift = M->get_var<double>("z_shift");
    }
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
    if(rot_trans.size()>9){
        rot_trans[9]=x_shift.eval();
        rot_trans[10]=y_shift.eval();
        rot_trans[11]=z_shift.eval();
    }
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
/*Roll Pitch Yaw in degrees
 Aplly Rotation on a point cloud*/
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
/*Roll Pitch Yaw in degrees*/
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
/* Update point cloud coordinates */
void update_xyz(vector<vector<double>>& point_cloud1, vector<double>& x_vec1, vector<double>& y_vec1, vector<double>& z_vec1){
    for (auto i = 0; i< point_cloud1.size(); i++) {
        point_cloud1[i][0] = x_vec1[i];
        point_cloud1[i][1] = y_vec1[i];
        point_cloud1[i][2] = z_vec1[i];
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

/* Centralize point cloud around origin */
void centralize(int n, POINT3D **  p, double avg_x, double avg_y, double avg_z){
    for(int i = 0; i < n; i++)
    {
        (*p)[i].x  -= avg_x;
        (*p)[i].y  -= avg_y;
        (*p)[i].z -= avg_z;
    }
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


/* Compute the L2 error for model and data sets
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
/* Read input files */
vector<pair<double, double>> read_data(const rapidcsv::Document& Model_doc,vector<vector<double>>& point_cloud, vector<vector<double>>& uav){
    vector<pair<double, double>> res(3);
    int model_nb_rows = Model_doc.GetRowCount();
    if(model_nb_rows<3){
        throw invalid_argument("Input file with less than 2 points");
    }
    DebugOn("Input file has " << model_nb_rows << " rows" << endl);
    point_cloud.resize(model_nb_rows);
    uav.resize(model_nb_rows);
    double xmin=numeric_limits<double>::max(), xmax=numeric_limits<double>::min(), ymin=numeric_limits<double>::max(), ymax=numeric_limits<double>::min(), zmin=numeric_limits<double>::max(), zmax=numeric_limits<double>::min();
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
        if(x<=xmin){
            xmin=x;
        }
        if(y<=ymin){
            ymin=y;
        }
        if(z<=zmin){
            zmin=z;
        }
        if(x>=xmax){
            xmax=x;
        }
        if(y>=ymax){
            ymax=y;
        }
        if(z>=zmax){
            zmax=z;
        }
    }
    res[0]={xmin, xmax};
    res[1]={ymin, ymax};
    res[2]={zmin, zmax};
    return res;
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
/* Read Laz files */
vector<vector<double>> read_laz(const string& fname){
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
    vector<vector<double>> uav_cloud;
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
        while (lasreader->read_point() && LidarPoints.size()!=200e6)
        {
            nb_pts++;
            //            if(nb_pts++<1e4)
            //                continue;
            //            if(nb_pts==0){
            //                DebugOn(to_string_with_precision(10.*(lasreader->point.get_gps_time()+315964800. - 18.),24) << ": (" << to_string_with_precision(lasreader->point.get_x(),10) <<"," << to_string_with_precision(lasreader->point.get_y(),10) << ","<< to_string_with_precision(lasreader->point.get_z(),10) <<")"<<endl);
            //                    //                return 0;
            //            }
            auto laser_id = lasreader->point.get_point_source_ID();
            DebugOff("lid "<<laser_id<<endl);
                        if(nb_pts%10!=0){/* Only keep points from Nadir laser */
                            continue;
                        }
            auto unix_time = lasreader->point.get_gps_time();
            auto x = lasreader->point.get_x();
            auto y = lasreader->point.get_y();
            auto z = lasreader->point.get_z();
            auto uav_x = lasreader->point.get_attribute_as_float(1);
            auto uav_y = lasreader->point.get_attribute_as_float(2);
            auto uav_z = lasreader->point.get_attribute_as_float(3);
            LidarPoints.push_back(new LidarPoint(laser_id,unix_time,x,y,z));
            point_cloud1.push_back({x,y,z});
            uav_cloud.push_back({uav_x,uav_y, uav_z});
//                        if(!xvals.insert(x*100).second){/* A U turn is being detected */
//                            u_turn = true;
//                            DebugOn("Detected a Uturn at point " << LidarPoints.size() << endl);
//                            if(u_turn) {
//                                DebugOn("This is the second Uturn! " << endl);
//                                u_turn_2 = true;
//                            }
//                            frame1 = false;
//                        }
//                        if(frame1){
//                            point_cloud1.push_back({x,y,z});
//                        }
//                        else{
//                            point_cloud2.push_back({x,y,z});
//                        }
        }

        DebugOn("Read " << LidarPoints.size() << " points" << endl);
        DebugOn(point_cloud1.size() << " points in flight line 1" << endl);
        DebugOn(point_cloud2.size() << " points in flight line 2" << endl);
        vector<vector<double>> empty_vec;
#ifdef USE_MATPLOT
       // plot(uav_cloud,empty_vec, 0.1);
#endif
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
    DebugOn("finished read laz"<<endl);
    return uav_cloud;
}
/*scale uav_cloud and then call this*/
vector<vector<double>> turn_detect(vector<vector<double>> uav_cloud){
    vector<vector<double>> u_list;
    vector<double> uturn;
    const double zero_tol=1e-12;
    const double tol=1e-2;
    bool line_found=false;
    bool turn_found=false;
    for(auto i=0;i<uav_cloud.size()-2;i++){
        double p1x=uav_cloud[i][0];
        double p1y=uav_cloud[i][1];
        double p2x=uav_cloud[i+1][0];
        double p2y=uav_cloud[i+1][1];
        double p3x=uav_cloud[i+2][0];
        double p3y=uav_cloud[i+2][1];
        
        double sl1=1000, sl2=2000, sl3=3000;
        if(abs(p2x-p1x)>=zero_tol){
            sl1=(p2y-p1y)/(p2x-p1x);
        }
        if(abs(p3x-p2x)>=zero_tol){
            sl2=(p3y-p2y)/(p3x-p2x);
        }
        if(abs(p3x-p1x)>=zero_tol){
            sl3=(p3y-p1y)/(p3x-p1x);
        }
        if(abs(sl1-sl2)<=tol && abs(sl2-sl3)<=tol && abs(sl1-sl3)<=tol){
            if(!line_found && u_list.empty()){
                u_list.push_back(uav_cloud[i+1]);
            }
            line_found=true;
        }
        else{
            if(line_found){
                turn_found=true;
                DebugOn("Turn found"<<endl);
                u_list.push_back(uav_cloud[i+1]);
            }
            line_found=false;
        }
    }
    return u_list;
}

/*scale uav_cloud and then call this*/
vector<vector<double>> turn_detect_slope_map(vector<vector<double>> uav_cloud){
    vector<vector<double>> u_list;
    map<double, pair<int, int>> slope_map;
    vector<double> uturn;
    const double zero_tol=1e-12;
    const double tol=1e-2;
    bool line_found=false;
    bool turn_found=false;
    for(auto i=0;i<uav_cloud.size()-10;i++){
        double p1x=uav_cloud[i][0];
        double p1y=uav_cloud[i][1];
        double p2x=uav_cloud[i+10][0];
        double p2y=uav_cloud[i+10][1];
        
        
        double sl1=1000, sl2=2000, sl3=3000;
        if(abs(p2x-p1x)>=zero_tol){
            sl1=(p2y-p1y)/(p2x-p1x);
        }
        
        if(!slope_map.insert(pair<double, pair<int,int>>(sl1, pair<int, int>(i, i+10))).second){
            auto p1_p2=slope_map[sl1];
            auto j=p1_p2.first;
            auto k=p1_p2.second;
            double p1x_o=uav_cloud[j][0];
            double p1y_o=uav_cloud[j][1];
            double p2x_o=uav_cloud[k][0];
            double p2y_o=uav_cloud[k][1];
           
            if(((p2x_o-p1x_o)<=-zero_tol && (p2x-p1x)>=zero_tol) || ((p2x_o-p1x_o)>=zero_tol && (p2x-p1x)<=-zero_tol) || ((p2y_o-p1y_o)<=-zero_tol && (p2y-p1y)>=zero_tol) || ((p2y_o-p1y_o)>=zero_tol && (p2y-p1y)<=-zero_tol)){
                u_list.push_back(uav_cloud[i]);
                u_list.push_back(uav_cloud[i+10]);
                u_list.push_back(uav_cloud[j]);
                u_list.push_back(uav_cloud[k]);
                DebugOn("a "<<setprecision(15)<<uav_cloud[i][0]<<" "<<setprecision(15)<<uav_cloud[i][1]<<endl);
                DebugOn("b "<<setprecision(15)<<uav_cloud[i+10][0]<<" "<<setprecision(15)<<uav_cloud[i+10][1]<<endl);
                DebugOn("c "<<setprecision(15)<<uav_cloud[j][0]<<" "<<setprecision(15)<<uav_cloud[j][1]<<endl);
                DebugOn("d "<<setprecision(15)<<uav_cloud[k][0]<<" "<<setprecision(15)<<uav_cloud[k][1]<<endl);
            }
            
        }
       
        
    }
    return u_list;
}

/*scale uav_cloud and then call this*/
vector<vector<double>> reg_slope_lines(vector<vector<double>> uav_cloud){
    vector<vector<double>> lines;
    multimap<double, vector<double>, greater <double>> dist_lines;
    vector<vector<double>> ulist;
    vector<double> uturn;
    const double zero_tol=1e-12;
    const double tol=1e-2;
    bool line_found=false;
    bool turn_found=false;
    for(auto i=1;i<uav_cloud.size();i++){
        double p1x=uav_cloud[i-1][0];
        double p1y=uav_cloud[i-1][1];
        double p2x=uav_cloud[i][0];
        double p2y=uav_cloud[i][1];
        
        double sl_prev=0;
        double c_prev=0;
        double sum_xy_prev,sum_xx_prev,sum_x_prev, sum_y_prev, n_prev;
        double begp, endp;
        
        if(!lines.empty()){
            sl_prev=lines.back()[0];
            c_prev=lines.back()[1];
            sum_y_prev=lines.back()[2];
            sum_x_prev=lines.back()[3];
            sum_xx_prev=lines.back()[4];
            sum_xy_prev=lines.back()[5];
            n_prev=lines.back()[6];
            begp=lines.back()[7];
            endp=lines.back()[8];
        }
        
        if(lines.empty() || abs(p2y-sl_prev*p2x-c_prev)/abs(p2y)>=1e-1){
            n_prev=2;
            sum_y_prev=p2y+p1y;
            sum_x_prev=p2x+p1x;
            sum_xx_prev=pow(p2x,2)+pow(p1x,2);
            sum_xy_prev=p2x*p2y+p1x*p1y;
            begp=i-1;
            endp=i;
//            if(!lines.empty() && lines.back()[6]==2){
//                lines.erase(lines.end()-1, lines.end());
//            }
            vector<double> res(9);
            lines.push_back(res);
            DebugOn("error "<<abs(p2y-sl_prev*p2x-c_prev)<<endl);
        }
        else{
            sum_y_prev+=p2y;
            sum_x_prev+=p2x;
            sum_xx_prev+=pow(p2x,2);
            sum_xy_prev+=p2x*p2y;
            n_prev++;
            endp=i;
        }
        if(abs(n_prev*sum_xx_prev-sum_x_prev*sum_x_prev)>=zero_tol){
            sl_prev=(n_prev*sum_xy_prev-sum_x_prev*sum_y_prev)/(n_prev*sum_xx_prev-sum_x_prev*sum_x_prev);
            c_prev=(sum_y_prev-sl_prev*sum_x_prev)/n_prev;
        }
        else{
            c_prev=sum_x_prev/n_prev;
            sl_prev=1e5;
        }
        lines.back()[0]=sl_prev;
        lines.back()[1]=c_prev;
        lines.back()[2]=sum_y_prev;
        lines.back()[3]=sum_x_prev;
        lines.back()[4]=sum_xx_prev;
        lines.back()[5]=sum_xy_prev;
        lines.back()[6]=n_prev;
        lines.back()[7]=begp;
        lines.back()[8]=endp;
    }
    
    auto s=lines.size();
        for(auto i=0;i<s;i++){
            int p1=lines[i][7];
            int p2=lines[i][8];
            auto x1=uav_cloud[p1][0];
            auto y1=uav_cloud[p1][1];
            auto x2=uav_cloud[p2][0];
            auto y2=uav_cloud[p2][1];
            DebugOn("p1 p2 "<<p1 <<" "<<p2<<endl);
            DebugOn("x1 y1 "<<x1 <<" "<<y1<<endl);
            DebugOn("x2 y2 "<<x2 <<" "<<y2<<endl);
            double d=pow(x2-x1,2)+pow(y2-y1,2);
            dist_lines.insert(pair<double, vector<double>>(d, lines[i]));
        }
    
    auto it=dist_lines.begin();
    int p1=it->second[7];
    int p2=it->second[8];
    it++;
    int p3=it->second[7];
    int p4=it->second[8];
    ulist.push_back(uav_cloud[p1]);
    ulist.push_back(uav_cloud[p2]);
    ulist.push_back(uav_cloud[p3]);
    ulist.push_back(uav_cloud[p4]);
    it++;
    int p5=it->second[7];
    int p6=it->second[8];
//    ulist.push_back(uav_cloud[p5]);
//    ulist.push_back(uav_cloud[p6]);
    

//    auto s=lines.size();
//    for(auto i=0;i<s-1;i++){
//        for(auto j=i+1;j<s-1;j++){
//            if(abs(lines[i][0]-lines[j][0])<=1e-3){
//                int p1=lines[i][7];
//                int p2=lines[i][8];
//                int p3=lines[j][7];
//                int p4=lines[j][8];
//                ulist.push_back(uav_cloud[p1]);
//                ulist.push_back(uav_cloud[p2]);
//                ulist.push_back(uav_cloud[p3]);
//                ulist.push_back(uav_cloud[p4]);
//            }
//        }
//    }
    return ulist;
}

vector<vector<double>> filter_z_slope(vector<vector<double>> uav_cloud){
    vector<vector<double>> res;
    const double angle_max=30;
    const double tan_angle_max=tan(30*pi/180);
    
    for(auto i=1;i<uav_cloud.size();i++){
        auto x2=uav_cloud.at(i)[0];
        auto y2=uav_cloud.at(i)[1];
        auto z2=uav_cloud.at(i)[2];
        auto x1=uav_cloud.at(i-1)[0];
        auto y1=uav_cloud.at(i-1)[1];
        auto z1=uav_cloud.at(i-1)[2];
        double slx=(z2-z1)/(x2-x1);
        double sly=(z2-z1)/(y2-y1);
        if(slx<tan_angle_max && sly <tan_angle_max){
            res.push_back(uav_cloud.at(i));
        }
    }
    return res;
}
//vector<double> projection_parallel()
//{
//    vector<vector<double>> u_list;
//    vector<double> uturn;
//    const double zero_tol=1e-12;
//    const double tol=1e-3;
//    bool line_found=false;
//    bool turn_found=false;
//    for(auto i=0;i<uav_cloud.size()-2;i++){
//        double p1x=uav_cloud[i][0];
//        double p1y=uav_cloud[i][1];
//        double p2x=uav_cloud[i+1][0];
//        double p2y=uav_cloud[i+1][1];
//        double p3x=uav_cloud[i+2][0];
//        double p3y=uav_cloud[i+2][1];
//
//        double sl1=1000, sl2=2000, sl3=3000;
//        if(abs(p2x-p1x)>=zero_tol){
//            sl1=(p2y-p1y)/(p2x-p1x);
//        }
//        if(abs(p3x-p2x)>=zero_tol){
//            sl2=(p3y-p2y)/(p3x-p2x);
//        }
//        if(abs(p3x-p2x)>=zero_tol){
//            sl3=(p3y-p1y)/(p3x-p1x);
//        }
//        if(abs(sl1-sl2)<=1e-3 && abs(sl2-sl3)<=1e-3 && abs(sl1-sl3)<=1e-3){
//            line_found=true;
//        }
//        else{
//            if(line_found){
//                turn_found=true;
//                DebugOn("Turn found"<<endl);
//                u_list.push_back(uav_cloud[i+1]);
//            }
//            line_found=false;
//        }
//    }
//    return u_list;
//}
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


#endif /* Lidar_utils_hpp */
