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
using namespace gravity;

#include <gravity/ConfigMap.hpp>

/* Computes the interpolation coefficient based on time elapsed  */
double get_interpolation_coef(const double& lidar_time, UAVPoint* p1, UAVPoint* p2);

/* Return the min-max values for x, y and z  for all possible rotations of p with angle +- angle*/
vector<pair<double,double>> get_min_max(double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, const vector<double>& p, const vector<double>& ref);


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
    DebugOn("Roll (degrees) = " << to_string_with_precision(roll_val,12) << endl);
    DebugOn("Pitch (degrees) = " << to_string_with_precision(pitch_val,12) << endl);
    DebugOn("Yaw (degrees) = " << to_string_with_precision(yaw_val,12) << endl);
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

void apply_rotation_new(double roll, double pitch, double yaw, vector<vector<double>>& full_point_cloud, vector<vector<double>>& full_uav, vector<vector<double>>& roll_pitch_yaw_uav, double scanner_x, double scanner_y, double scanner_z){
    double shifted_x, shifted_y, shifted_z;
    /* Apply rotation */
    for (auto i = 0; i< full_point_cloud.size(); i++) {
        
        auto rollu=(-pi+roll_pitch_yaw_uav[i][1])*(-1);
        auto pitchu=roll_pitch_yaw_uav[i][0];
        auto yawu=(-pi/2+roll_pitch_yaw_uav[i][2])*(-1);
        
        auto tx=scanner_x*cos(pitchu)*cos(yawu)+scanner_y*(cos(rollu)*sin(yawu) +cos(yawu)*sin(rollu)*sin(pitchu))+scanner_z*(sin(rollu)*sin(yawu)-cos(rollu)*cos(yawu)*sin(pitchu));
        auto ty=-scanner_x*cos(pitchu)*sin(yawu)+scanner_y*(cos(rollu)*cos(yawu)-sin(rollu)*sin(pitchu)*sin(yawu))+scanner_z*((cos(yawu)*sin(rollu)+cos(rollu)*sin(pitchu)*sin(yawu)));
        auto tz=scanner_x*sin(pitchu)+scanner_y*(-cos(pitchu)*sin(rollu))+scanner_z*(cos(rollu)*cos(pitchu));
        
        shifted_x = full_point_cloud[i][0];// -tx;
        shifted_y = full_point_cloud[i][1];// -ty;
        shifted_z = full_point_cloud[i][2];// -tz;
       
        
        full_point_cloud[i][0] = shifted_x*cos(pitchu)*cos(yawu) - shifted_y*cos(pitchu)*sin(yawu) + shifted_z*sin(pitchu);
        full_point_cloud[i][1] = shifted_x*(cos(rollu)*sin(yawu) +cos(yawu)*sin(rollu)*sin(pitchu))+ shifted_y*(cos(rollu)*cos(yawu)-sin(rollu)*sin(pitchu)*sin(yawu))  + shifted_z*(-cos(pitchu)*sin(rollu));
        full_point_cloud[i][2] = shifted_x*(sin(rollu)*sin(yawu)-cos(rollu)*cos(yawu)*sin(pitchu)) + shifted_y*(cos(yawu)*sin(rollu)+cos(rollu)*sin(pitchu)*sin(yawu)) + shifted_z*(cos(rollu)*cos(pitchu));
        
        


        auto fx = full_point_cloud[i][0]*cos(pitch)*cos(yaw) - full_point_cloud[i][1]*cos(pitch)*sin(yaw) + full_point_cloud[i][2]*sin(pitch);
        auto fy = full_point_cloud[i][0]*(cos(roll)*sin(yaw) +cos(yaw)*sin(roll)*sin(pitch))+ full_point_cloud[i][1]*(cos(roll)*cos(yaw)-sin(roll)*sin(pitch)*sin(yaw))  + full_point_cloud[i][2]*(-cos(pitch)*sin(roll));
        auto fz = full_point_cloud[i][0]*(sin(roll)*sin(yaw)-cos(roll)*cos(yaw)*sin(pitch)) + full_point_cloud[i][1]*(cos(yaw)*sin(roll)+cos(roll)*sin(pitch)*sin(yaw)) + full_point_cloud[i][2]*(cos(roll)*cos(pitch));
        
        vector<double> pnew(3);
        
        pnew[0]=fx*cos(pitchu)*cos(yawu)+fy*(cos(rollu)*sin(yawu) +cos(yawu)*sin(rollu)*sin(pitchu))+fz*(sin(rollu)*sin(yawu)-cos(rollu)*cos(yawu)*sin(pitchu));
        pnew[1]=-fx*cos(pitchu)*sin(yawu)+fy*(cos(rollu)*cos(yawu)-sin(rollu)*sin(pitchu)*sin(yawu))+fz*((cos(yawu)*sin(rollu)+cos(rollu)*sin(pitchu)*sin(yawu)));
      
        
        
//        auto tx=scanner_x*cos(pitchu)*cos(yawu) - scanner_y*cos(pitchu)*sin(yawu) + scanner_z*sin(pitchu);
//        auto ty = scanner_x*(cos(rollu)*sin(yawu) +cos(yawu)*sin(rollu)*sin(pitchu))+ scanner_y*(cos(rollu)*cos(yawu)-sin(rollu)*sin(pitchu)*sin(yawu))  + scanner_z*(-cos(pitchu)*sin(rollu));
//        auto tz = scanner_x*(sin(rollu)*sin(yawu)-cos(rollu)*cos(yawu)*sin(pitchu)) + scanner_y*(cos(yawu)*sin(rollu)+cos(rollu)*sin(pitchu)*sin(yawu)) + scanner_z*(cos(rollu)*cos(pitchu));
        
     
        
        
        

        full_point_cloud[i][0] = pnew[0]+full_uav[i][0];//+tx;
        full_point_cloud[i][1] = pnew[1]+full_uav[i][1];//+ty;
        full_point_cloud[i][2] = pnew[2]+full_uav[i][2];//+tz;
    }
    
}
vector<double> apply_rotation_new_order(double roll, double pitch, double yaw, double x, double y, double z){
    vector<double> res(3);
    res[0] = x*cos(roll)*cos(yaw) - y*cos(roll)*sin(yaw) + z*sin(roll);
    res[1] = x*(cos(pitch)*sin(yaw) +cos(yaw)*sin(roll)*sin(pitch))+ y*(cos(pitch)*cos(yaw)-sin(roll)*sin(pitch)*sin(yaw))  + z*(-cos(roll)*sin(pitch));
    res[2] =x*(sin(pitch)*sin(yaw)-cos(pitch)*cos(yaw)*sin(roll)) + y*(cos(yaw)*sin(pitch)+cos(pitch)*sin(roll)*sin(yaw)) + z*(cos(roll)*cos(pitch));
    return res;
}

vector<double> apply_rotation_inverse_new_order(double roll, double pitch, double yaw, double x, double y, double z){
    vector<double> res(3);
    res[0] = x*cos(roll)*cos(yaw) + y*(cos(pitch)*sin(yaw) +cos(yaw)*sin(roll)*sin(pitch)) + z*(sin(pitch)*sin(yaw)-cos(pitch)*cos(yaw)*sin(roll));
    res[1]=x*cos(roll)*sin(yaw)*(-1)+y*(cos(pitch)*cos(yaw)-sin(roll)*sin(pitch)*sin(yaw))+z*(cos(yaw)*sin(pitch)+cos(pitch)*sin(roll)*sin(yaw));
    res[2]=x*sin(roll)+y*(-cos(roll)*sin(pitch))+z*(cos(roll)*cos(pitch));
    return res;
}
void apply_transform_new_order(double roll, double pitch, double yaw, vector<vector<double>>& full_point_cloud, const vector<vector<double>>& full_uav, const vector<vector<double>>& roll_pitch_yaw_uav, const double scanner_x, const double scanner_y, const double scanner_z, const double hr, const double hp, const double hy){
    double shifted_x, shifted_y, shifted_z;
    /* Apply rotation */
    for (auto i = 0; i< full_point_cloud.size(); i++) {
        
        auto rollu=roll_pitch_yaw_uav[i][0];
        auto pitchu=(-pi+roll_pitch_yaw_uav[i][1])*(-1);
        auto yawu=(-pi/2+roll_pitch_yaw_uav[i][2])*(-1);
        
        auto res_s=apply_rotation_inverse_new_order(roll_pitch_yaw_uav[i][0], (roll_pitch_yaw_uav[i][1]), (-pi/2+roll_pitch_yaw_uav[i][2])*(-1), scanner_x, scanner_y, scanner_z);
        
        shifted_x = full_point_cloud[i][0]-full_uav[i][0] -res_s[0];// -tx;
        shifted_y = full_point_cloud[i][1]-full_uav[i][1] -res_s[1];// -ty;
        shifted_z = full_point_cloud[i][2]-full_uav[i][2] -res_s[2];// -tz;
        
        auto res_inv_ins1=apply_rotation_new_order(rollu, pitchu, yawu, shifted_x, shifted_y, shifted_z);
        auto res_inv_ins=apply_rotation_inverse_new_order(hr, hp, hy, res_inv_ins1[0], res_inv_ins1[1], res_inv_ins1[2]);
        
        auto res_bore=apply_rotation_new_order(roll, pitch, yaw, res_inv_ins[0], res_inv_ins[1], res_inv_ins[2]);
        
        auto res_ins=apply_rotation_inverse_new_order(rollu, pitchu, yawu, res_bore[0], res_bore[1], res_bore[2]);
        
        full_point_cloud[i][0] = res_ins[0]+full_uav[i][0]+res_s[0];//+tx;
        full_point_cloud[i][1] = res_ins[1]+full_uav[i][1]+res_s[1];//+ty;
        full_point_cloud[i][2] = res_ins[2]+full_uav[i][2]+res_s[2];//+tz;
    }
}

void generate_inputs(vector<vector<double>>& full_point_cloud, const vector<vector<double>>& full_uav, const vector<vector<double>>& roll_pitch_yaw_uav, const double scanner_x, const double scanner_y, const double scanner_z, const double hr, const double hp, const double hy, vector<vector<double>>& input_point_cloud, vector<vector<double>>& input_rot_offset){
    double shifted_x, shifted_y, shifted_z;
    /* Apply rotation */
    for (auto i = 0; i< full_point_cloud.size(); i++) {
        
        auto rollu=roll_pitch_yaw_uav[i][0];
        auto pitchu=(-pi+roll_pitch_yaw_uav[i][1])*(-1);
        auto yawu=(-pi/2+roll_pitch_yaw_uav[i][2])*(-1);
        
        auto res_s=apply_rotation_inverse_new_order(roll_pitch_yaw_uav[i][0], (roll_pitch_yaw_uav[i][1]), (-pi/2+roll_pitch_yaw_uav[i][2])*(-1), scanner_x, scanner_y, scanner_z);
        
        shifted_x = full_point_cloud[i][0]-full_uav[i][0] -res_s[0];// -tx;
        shifted_y = full_point_cloud[i][1]-full_uav[i][1] -res_s[1];// -ty;
        shifted_z = full_point_cloud[i][2]-full_uav[i][2] -res_s[2];// -tz;
        
        auto res_inv_ins1=apply_rotation_new_order(rollu, pitchu, yawu, shifted_x, shifted_y, shifted_z);
       auto res_inv_ins=apply_rotation_inverse_new_order(hr,hp, hy, res_inv_ins1[0], res_inv_ins1[1], res_inv_ins1[2]);
        input_point_cloud.push_back(res_inv_ins);
        input_rot_offset.push_back(res_s);
    }
}
void generate_outputs_from_inputs(double roll, double pitch, double yaw, const vector<vector<double>>& input_point_cloud, const vector<vector<double>>& uav, const vector<vector<double>>& roll_pitch_yaw_uav, const vector<vector<double>>& input_offset, vector<vector<double>>& output_point_cloud){
    /* Apply rotation */
    for (auto i = 0; i< input_point_cloud.size(); i++) {
        vector<double> res(3);
        auto rollu=roll_pitch_yaw_uav[i][0];
        auto pitchu=(-pi+roll_pitch_yaw_uav[i][1])*(-1);
        auto yawu=(-pi/2+roll_pitch_yaw_uav[i][2])*(-1);
        
        
        auto res_bore=apply_rotation_new_order(roll, pitch, yaw, input_point_cloud[i][0], input_point_cloud[i][1], input_point_cloud[i][2]);
        
        auto res_ins=apply_rotation_inverse_new_order(rollu, pitchu, yawu, res_bore[0], res_bore[1], res_bore[2]);
        
        res[0] = res_ins[0]+uav[i][0]+input_offset[i][0];
        res[1] = res_ins[1]+uav[i][1]+input_offset[i][1];
        res[2] = res_ins[2]+uav[i][2]+input_offset[i][2];
        
        output_point_cloud.push_back(res);
    }
}
void apply_transform_new_order_Test(double roll, double pitch, double yaw, vector<vector<double>>& full_point_cloud, const vector<vector<double>>& full_uav, const vector<vector<double>>& roll_pitch_yaw_uav, double scanner_x, double scanner_y, double scanner_z){
    scanner_x=0;
    scanner_y=0;
    scanner_z=0;
    double shifted_x, shifted_y, shifted_z;
    /* Apply rotation */
    for (auto i = 0; i< full_point_cloud.size(); i++) {
        
        auto rollu=roll_pitch_yaw_uav[i][0];
        auto pitchu=(-pi+roll_pitch_yaw_uav[i][1])*(-1);
        auto yawu=(-pi/2+roll_pitch_yaw_uav[i][2])*(-1);
        
        auto res_s=apply_rotation_inverse_new_order(roll_pitch_yaw_uav[i][0], (roll_pitch_yaw_uav[i][1]), (-pi/2+roll_pitch_yaw_uav[i][2])*(-1), scanner_x, scanner_y, scanner_z);
        
        shifted_x = full_point_cloud[i][0]-full_uav[i][0] -res_s[0];// -tx;
        shifted_y = full_point_cloud[i][1]-full_uav[i][1] -res_s[1];// -ty;
        shifted_z = full_point_cloud[i][2]-full_uav[i][2] -res_s[2];// -tz;
        
        auto res_inv_ins=apply_rotation_new_order(rollu, pitchu, yawu, shifted_x, shifted_y, shifted_z);
        
        auto res_bore=apply_rotation_new_order(roll, pitch, yaw, res_inv_ins[0], res_inv_ins[1], res_inv_ins[2]);
        
        auto res_ins=apply_rotation_inverse_new_order(rollu, pitchu, yawu, res_bore[0], res_bore[1], res_bore[2]);
        
        full_point_cloud[i][0] = res_ins[0]+full_uav[i][0]+res_s[0];//+tx;
        full_point_cloud[i][1] = res_ins[1]+full_uav[i][1]+res_s[1];//+ty;
        full_point_cloud[i][2] = res_ins[2]+full_uav[i][2]+res_s[2];//+tz;
    }
}

void apply_transform_new_order_scanner(double roll, double pitch, double yaw, vector<vector<double>>& full_point_cloud, const vector<vector<double>>& full_uav, const vector<vector<double>>& roll_pitch_yaw_uav, double scanner_x, double scanner_y, double scanner_z){
    double shifted_x, shifted_y, shifted_z;
    /* Apply rotation */
    for (auto i = 0; i< full_point_cloud.size(); i++) {
        
        auto rollu=roll_pitch_yaw_uav[i][0];
        auto pitchu=(-pi+roll_pitch_yaw_uav[i][1])*(-1);
        auto yawu=(-pi/2+roll_pitch_yaw_uav[i][2])*(-1);
        
        shifted_x = full_point_cloud[i][0]-full_uav[i][0];//-scanner_x;
        shifted_y = full_point_cloud[i][1]-full_uav[i][1];//-scanner_y;
        shifted_z = full_point_cloud[i][2]-full_uav[i][2];//-scanner_z;
        
        auto res_inv_ins=apply_rotation_new_order(rollu, pitchu, yawu, shifted_x, shifted_y, shifted_z);
        
        auto res_bore=apply_rotation_new_order(roll, pitch, yaw, res_inv_ins[0]-scanner_x, res_inv_ins[1]-scanner_y, res_inv_ins[2]-scanner_z);
        
        //auto res_bore=apply_rotation_new_order(roll, pitch, yaw, res_inv_ins[0], res_inv_ins[1], res_inv_ins[2]);
        
        auto res_ins=apply_rotation_inverse_new_order(rollu, pitchu, yawu, res_bore[0], res_bore[1], res_bore[2]);
        
        auto res_ins_scanner=apply_rotation_inverse_new_order(rollu, pitchu, yawu, scanner_x, scanner_y, scanner_z);
        
        full_point_cloud[i][0] = res_ins[0]+full_uav[i][0]+res_ins_scanner[0];
        full_point_cloud[i][1] = res_ins[1]+full_uav[i][1]+res_ins_scanner[1];
        full_point_cloud[i][2] = res_ins[2]+full_uav[i][2]+res_ins_scanner[2];
//        full_point_cloud[i][0] = res_ins[0]+full_uav[i][0]+scanner_x;
//        full_point_cloud[i][1] = res_ins[1]+full_uav[i][1]+scanner_y;
//        full_point_cloud[i][2] = res_ins[2]+full_uav[i][2]+scanner_z;
    }
}

void apply_rotation_bore_ins_new_order(double roll, double pitch, double yaw, vector<vector<double>>& full_point_cloud, vector<vector<double>>& full_uav, vector<vector<double>>& roll_pitch_yaw_uav, double scanner_x, double scanner_y, double scanner_z){
    double shifted_x, shifted_y, shifted_z;
    /* Apply rotation */
    for (auto i = 0; i< full_point_cloud.size(); i++) {
        
        auto rollu=roll_pitch_yaw_uav[i][0];
        auto pitchu=(-pi+roll_pitch_yaw_uav[i][1])*(-1);
        auto yawu=(-pi/2+roll_pitch_yaw_uav[i][2])*(-1);
        
        
        auto res_bore=apply_rotation_new_order(roll, pitch, yaw, full_point_cloud[i][0], full_point_cloud[i][1], full_point_cloud[i][2]);
        
        auto res_ins=apply_rotation_inverse_new_order(rollu, pitchu, yawu, res_bore[0], res_bore[1], res_bore[2]);

        full_point_cloud[i][0] = res_ins[0]+full_uav[i][0];//+tx;
        full_point_cloud[i][1] = res_ins[1]+full_uav[i][1];//+ty;
        full_point_cloud[i][2] = res_ins[2]+full_uav[i][2];//+tz;
    }
}


void apply_rotation_lidarviewer(double roll, double pitch, double yaw, vector<vector<double>>& full_point_cloud, vector<vector<double>>& full_uav, vector<vector<double>>& roll_pitch_yaw_uav, double scanner_x, double scanner_y, double scanner_z){
    double shifted_x, shifted_y, shifted_z;
    /* Apply rotation */
    for (auto i = 0; i< full_point_cloud.size(); i++) {
        
        auto rollu=(-pi+roll_pitch_yaw_uav[i][1])*(-1);
        auto pitchu=roll_pitch_yaw_uav[i][0];
        auto yawu=(-pi/2+roll_pitch_yaw_uav[i][2])*(-1);
        
        auto tx=scanner_x*cos(pitchu)*cos(yawu)+scanner_y*(cos(rollu)*sin(yawu) +cos(yawu)*sin(rollu)*sin(pitchu))+scanner_z*(sin(rollu)*sin(yawu)-cos(rollu)*cos(yawu)*sin(pitchu));
        auto ty=-scanner_x*cos(pitchu)*sin(yawu)+scanner_y*(cos(rollu)*cos(yawu)-sin(rollu)*sin(pitchu)*sin(yawu))+scanner_z*((cos(yawu)*sin(rollu)+cos(rollu)*sin(pitchu)*sin(yawu)));
        auto tz=scanner_x*sin(pitchu)+scanner_y*(-cos(pitchu)*sin(rollu))+scanner_z*(cos(rollu)*cos(pitchu));
        
        shifted_x = full_point_cloud[i][0];// -tx;
        shifted_y = full_point_cloud[i][1];// -ty;
        shifted_z = full_point_cloud[i][2];// -tz;
        
        full_point_cloud[i][1]*=-1;
        full_point_cloud[i][2]*=-1;
        


        full_point_cloud[i][0] = shifted_x*cos(pitchu)*cos(yawu) - shifted_y*cos(pitchu)*sin(yawu) + shifted_z*sin(pitchu);
        full_point_cloud[i][1] = (shifted_x*(cos(rollu)*sin(yawu) +cos(yawu)*sin(rollu)*sin(pitchu))+ shifted_y*(cos(rollu)*cos(yawu)-sin(rollu)*sin(pitchu)*sin(yawu))  + shifted_z*(-cos(pitchu)*sin(rollu)));
        full_point_cloud[i][2] =(shifted_x*(sin(rollu)*sin(yawu)-cos(rollu)*cos(yawu)*sin(pitchu)) + shifted_y*(cos(yawu)*sin(rollu)+cos(rollu)*sin(pitchu)*sin(yawu)) + shifted_z*(cos(rollu)*cos(pitchu)));
//
        auto rolla=rollu-roll;
        auto pitcha=(pitchu-pitch)*(-1);
        auto yawa=(yawu-yaw)*(-1);
       
      


        auto fx = full_point_cloud[i][0]*cos(pitcha)*cos(yawa) + full_point_cloud[i][1]*(cos(rolla)*sin(yawa) +cos(yawa)*sin(rolla)*sin(pitcha)) + full_point_cloud[i][2]*(sin(rolla)*sin(yawa)-cos(rolla)*cos(yawa)*sin(pitcha));
        auto fy = -full_point_cloud[i][0]*cos(pitcha)*sin(yawa)+ full_point_cloud[i][1]*(cos(rolla)*cos(yawa)-sin(rolla)*sin(pitcha)*sin(yawa))  + full_point_cloud[i][2]*((cos(yawa)*sin(rolla)+cos(rolla)*sin(pitcha)*sin(yawa)));
        auto fz = full_point_cloud[i][0]*sin(pitcha) + full_point_cloud[i][1]*(-cos(pitcha)*sin(rolla)) + full_point_cloud[i][2]*(cos(rolla)*cos(pitcha));
        
      
        
        
//        auto tx=scanner_x*cos(pitchu)*cos(yawu) - scanner_y*cos(pitchu)*sin(yawu) + scanner_z*sin(pitchu);
//        auto ty = scanner_x*(cos(rollu)*sin(yawu) +cos(yawu)*sin(rollu)*sin(pitchu))+ scanner_y*(cos(rollu)*cos(yawu)-sin(rollu)*sin(pitchu)*sin(yawu))  + scanner_z*(-cos(pitchu)*sin(rollu));
//        auto tz = scanner_x*(sin(rollu)*sin(yawu)-cos(rollu)*cos(yawu)*sin(pitchu)) + scanner_y*(cos(yawu)*sin(rollu)+cos(rollu)*sin(pitchu)*sin(yawu)) + scanner_z*(cos(rollu)*cos(pitchu));
        
     
        
        
        

        full_point_cloud[i][0] = fx+full_uav[i][0];//+tx;
        full_point_cloud[i][1] = fy+full_uav[i][1];//+ty;
        full_point_cloud[i][2] = fz+full_uav[i][2];//+tz;
    }
    
}

vector<double> apply_rotation_transpose_new(double roll, double pitch, double yaw, double x, double y, double z){
    vector<double> pnew(3);
    
    pnew[0]=x*cos(pitch)*cos(yaw)+y*(cos(roll)*sin(yaw) +cos(yaw)*sin(roll)*sin(pitch))+z*(sin(roll)*sin(yaw)-cos(roll)*cos(yaw)*sin(pitch));
    pnew[1]=-x*cos(pitch)*sin(yaw)+y*(cos(roll)*cos(yaw)-sin(roll)*sin(pitch)*sin(yaw))+z*((cos(yaw)*sin(roll)+cos(roll)*sin(pitch)*sin(yaw)));
    pnew[2]=x*sin(pitch)+y*(-cos(pitch)*sin(roll))+z*(cos(roll)*cos(pitch));
    return pnew;

}

void apply_rotation(double roll, double pitch, double yaw, vector<vector<double>>& full_point_cloud, vector<vector<double>>& full_uav, vector<vector<double>>& roll_pitch_yaw_uav){
    double shifted_x, shifted_y, shifted_z;
    double beta = roll;// roll in radians
    double gamma = pitch; // pitch in radians
    double alpha = yaw; // yaw in radians
    /* Apply rotation */
    for (auto i = 0; i< full_point_cloud.size(); i++) {
        shifted_x = full_point_cloud[i][0];
        shifted_y = full_point_cloud[i][1];
        shifted_z = full_point_cloud[i][2] ;
        auto rollu=roll_pitch_yaw_uav[i][0];
        auto pitchu=roll_pitch_yaw_uav[i][1];
        auto yawu=roll_pitch_yaw_uav[i][2];
        
        auto px = shifted_x*cos(alpha)*cos(beta) + shifted_y*(cos(alpha)*sin(beta)*sin(gamma) - sin(alpha)*cos(gamma)) + shifted_z*(cos(alpha)*sin(beta)*cos(gamma) + sin(alpha)*sin(gamma));
       auto py = shifted_x*sin(alpha)*cos(beta) + shifted_y*(sin(alpha)*sin(beta)*sin(gamma) + cos(alpha)*cos(gamma)) + shifted_z*(sin(alpha)*sin(beta)*cos(gamma) - cos(alpha)*sin(gamma));
       auto pz = shifted_x*(-sin(beta)) + shifted_y*(cos(beta)*sin(gamma)) + shifted_z*(cos(beta)*cos(gamma));
   
        double betau = rollu;// roll in radians
        double gammau = pitchu; // pitch in radians
        double alphau = yawu; // yaw in radians

        auto fx = px*cos(alphau)*cos(betau) + py*(cos(alphau)*sin(betau)*sin(gammau) - sin(alphau)*cos(gammau)) + pz*(cos(alphau)*sin(betau)*cos(gammau) + sin(alphau)*sin(gammau));
       auto fy = px*sin(alphau)*cos(betau) + py*(sin(alphau)*sin(betau)*sin(gammau) + cos(alphau)*cos(gammau)) + pz*(sin(alphau)*sin(betau)*cos(gammau) - cos(alphau)*sin(gammau));
       auto fz = px*(-sin(betau)) + py*(cos(betau)*sin(gammau)) + pz*(cos(betau)*cos(gammau));
        
        

        full_point_cloud[i][0] = fx+full_uav[i][0];
        full_point_cloud[i][1] = fy+full_uav[i][1];
        full_point_cloud[i][2] = fz+full_uav[i][2];
    }
    
}

vector<double> apply_rotation_transpose(double roll, double pitch, double yaw, double x, double y, double z){
    
    double beta = roll;// roll in radians
    double gamma = pitch; // pitch in radians
    double alpha = yaw; // yaw in radians
    
    vector<double> pnew(3);
    
    
    pnew[0]=x*cos(alpha)*cos(beta)+y*sin(alpha)*cos(beta) +z*(-sin(beta));
    pnew[1]=x*(cos(alpha)*sin(beta)*sin(gamma) - sin(alpha)*cos(gamma))+y*(sin(alpha)*sin(beta)*sin(gamma) + cos(alpha)*cos(gamma))+z*(cos(beta)*sin(gamma));
    pnew[2]=x*(cos(alpha)*sin(beta)*cos(gamma) + sin(alpha)*sin(gamma))+y*(sin(alpha)*sin(beta)*cos(gamma) - cos(alpha)*sin(gamma))+z*(cos(beta)*cos(gamma));
    return pnew;

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
int get_sign(double y){
    int res;
    const double zero_tol=-1e-4;
    if(y>=zero_tol){
        res=1;
    }
    //    else if(y>=-zero_tol && y<=zero_tol){
    //        res=0;
    //    }
    else{
        res=-1;
    }
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
    lasheader.x_offset = 500000.0;
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
vector<vector<double>> read_laz(const string& fname, vector<vector<double>>& lidar_point_cloud, vector<vector<double>>& roll_pitch_yaw){
    string namef= fname.substr(0,fname.find('.'));
    string name=namef+"_original.laz";
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
    double time_start, time_end;
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
        DebugOn("xscale = "<<lasreader->header.x_scale_factor<<endl);
        DebugOn("yscale = "<<lasreader->header.y_scale_factor<<endl);
        DebugOn("zscale = "<<lasreader->header.z_scale_factor<<endl);
        DebugOn("xoffset = "<<lasreader->header.x_offset<<endl);
        DebugOn("yoffset = "<<lasreader->header.y_offset<<endl);
        DebugOn("zoffset = "<<lasreader->header.z_offset<<endl);
       
        
        int total_pts=lasreader->npoints;
        int skip=1;
        if(total_pts>=15e6 && total_pts<=1e8){
            skip=10;
        }
        else if(total_pts>1e8){
            skip=1000;
        }
        int nb_dots; /* Number of measurements inside cell */
        int xpos, ypos;
        double z, min_z, max_z, av_z;
        pair<int,int> pos;
        size_t nb_pts = 0;
        tuple<double,double,double,double,UAVPoint*> cell; /* <min_z,max_z,av_z> */
       
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
        
        
        lasreader->read_point();
        auto time_start = lasreader->point.get_gps_time();
        DebugOn("time_start "<<time_start<<endl);
        while (lasreader->read_point() && LidarPoints.size()!=200e6)
        {
            nb_pts++;
            time_end = lasreader->point.get_gps_time();
            auto laser_id = lasreader->point.get_point_source_ID();
            auto unix_time = lasreader->point.get_gps_time();
            DebugOff("lid "<<laser_id<<endl);
//            if(!((unix_time-time_start)>=150 && (unix_time-time_start)<=250)){
//                continue;
//            }
//            if(nb_pts%10!=0){/* Only keep points from Nadir laser */
//                continue;
//            }
            auto X = (lasreader->point.get_X());
            auto Y = (lasreader->point.get_Y());
            auto Z = (lasreader->point.get_Z());
            auto x = (lasreader->point.get_x());
            auto y = (lasreader->point.get_y());
            auto z = (lasreader->point.get_z());
            auto uav_x = (lasreader->point.get_attribute_as_float(1));
            auto uav_y = (lasreader->point.get_attribute_as_float(2));
            auto uav_z = lasreader->point.get_attribute_as_float(3);
            auto roll = lasreader->point.get_attribute_as_float(4);
            auto pitch = lasreader->point.get_attribute_as_float(5);
            auto yaw = lasreader->point.get_attribute_as_float(6);
            for(int i = 0; i < 11; i++){
                DebugOff("attribute "<< i << " = " << lasreader->point.get_attribute_name(i) << endl);
                DebugOff("attribute "<< i << " value = " << lasreader->point.get_attribute_as_float(i) << endl);
            }
            DebugOff("rpy "<<to_string_with_precision(roll,9)<<" "<<to_string_with_precision(pitch,9)<<" "<<to_string_with_precision(yaw,9)<<endl);
            DebugOff("uav "<<to_string_with_precision(uav_x,9)<<" "<<to_string_with_precision(uav_y,9)<<" "<<to_string_with_precision(uav_z,9)<<endl);
            DebugOff("lidar "<<to_string_with_precision(x,9)<<" "<<to_string_with_precision(y,9)<<" "<<to_string_with_precision(z,9)<<endl);
            if(!isnan(uav_x) && !isnan(uav_y) && !isnan(uav_z)){
                LidarPoints.push_back(new LidarPoint(laser_id,unix_time,x,y,z));
                point_cloud1.push_back({x,y,z});
                uav_cloud.push_back({uav_x,uav_y, uav_z});
                lidar_point_cloud.push_back({x,y,z});
                roll_pitch_yaw.push_back({roll,pitch,yaw});
            }
        }
        
        DebugOn("Read " << LidarPoints.size() << " points" << endl);
        DebugOn(point_cloud1.size() << " points in flight line 1" << endl);
        DebugOn(point_cloud2.size() << " points in flight line 2" << endl);
        vector<vector<double>> empty_vec;
        empty_vec.push_back(uav_cloud[0]);
#ifdef USE_MATPLOT
         plot(uav_cloud,empty_vec, 0.1);
#endif
        empty_vec.clear();
        empty_vec.push_back(lidar_point_cloud[0]);
#ifdef USE_MATPLOT
         plot(lidar_point_cloud,empty_vec, 0.1);
#endif
        save_laz(name, point_cloud1, point_cloud2);
        DebugOn("time_start "<<time_start);
        DebugOn("time_end "<<time_end-time_start);

    }
    DebugOn("finished read laz"<<endl);

    return uav_cloud;
}
void read_laz_new(const string& fname1, const string& fname2)
{
   string namef= fname1.substr(0,fname1.find('.'));
   string name=namef+"_original_2fold.laz";
   LASreadOpener lasreadopener1;
   lasreadopener1.set_file_name(fname1.c_str());
   lasreadopener1.set_populate_header(TRUE);
   param<> x1("x1"), y1("x1"), z1("x1");
   int xdim1=0, ydim1=0, zdim1=0;
   if (!lasreadopener1.active())
   {
       throw invalid_argument("ERROR: no input specified\n");
   }
 
   vector<vector<double>> point_cloud1, point_cloud2;
   while (lasreadopener1.active())
   {

       LASreader* lasreader = lasreadopener1.open();
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

       int total_pts=lasreader->npoints;
       int skip=1;
       if(total_pts>=15e6 && total_pts<=1e8){
           skip=10;
       }
       else if(total_pts>1e8){
           skip=100;
       }
       int nb_dots; /* Number of measurements inside cell */
       int xpos, ypos;
       double z, min_z, max_z, av_z;
       pair<int,int> pos;
       size_t nb_pts = 0;
       tuple<double,double,double,double,UAVPoint*> cell; /* <min_z,max_z,av_z> */
       

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
       while (lasreader->read_point() && LidarPoints.size()!=200e6)
       {
           nb_pts++;
           
           auto laser_id = lasreader->point.get_point_source_ID();
           DebugOff("lid "<<laser_id<<endl);
           if(nb_pts%skip!=0){/* Only keep points from Nadir laser */
               continue;
           }
           auto unix_time = lasreader->point.get_gps_time();
           auto x = lasreader->point.get_x();
           auto y = lasreader->point.get_y();
           auto z = lasreader->point.get_z();
           auto uav_x = lasreader->point.get_attribute_as_float(1);
           auto uav_y = lasreader->point.get_attribute_as_float(2);
           auto uav_z = lasreader->point.get_attribute_as_float(3);
           auto roll = lasreader->point.get_attribute_as_float(4);
           auto pitch = lasreader->point.get_attribute_as_float(5);
           auto yaw = lasreader->point.get_attribute_as_float(6);
           for(int i = 0; i < 11; i++){
               DebugOff("attribute "<< i << " = " << lasreader->point.get_attribute_name(i) << endl);
               DebugOff("attribute "<< i << " value = " << lasreader->point.get_attribute_as_float(i) << endl);
           }
           DebugOff("uav "<<to_string_with_precision(uav_x,9)<<" "<<to_string_with_precision(uav_y,9)<<" "<<to_string_with_precision(uav_z,9)<<endl);
           if(!isnan(uav_x) && !isnan(uav_y) && !isnan(uav_z)){
               LidarPoints.push_back(new LidarPoint(laser_id,unix_time,x,y,z));
               point_cloud1.push_back({x,y,z});
               //uav_cloud.push_back({uav_x,uav_y, uav_z});
               //lidar_point_cloud.push_back({x,y,z});
              // roll_pitch_yaw.push_back({roll,pitch,yaw});
           }
          
       }
   }
   LASreadOpener lasreadopener;
   lasreadopener.set_file_name(fname2.c_str());
   lasreadopener.set_populate_header(TRUE);
   if (!lasreadopener.active())
   {
       throw invalid_argument("ERROR: no input specified\n");
   }

   while (lasreadopener.active())
   {
          
   set<double> timestamps;
   vector<vector<double>> uav_cloud;
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

       int total_pts=lasreader->npoints;
       int skip=1;
       if(total_pts>=15e6 && total_pts<=1e8){
           skip=10;
       }
       else if(total_pts>1e8){
           skip=100;
       }
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
       set<int> xvals;
       size_t uav_id = 0;
       bool new_uav = true, u_turn = false, frame1 = true, u_turn_2=false;
       double unix_time, delta_x = 0, delta_y = 0;
       pair<map<int,shared_ptr<Frame>>::iterator,bool> frame_ptr;
       bool exit = false;
       while (lasreader->read_point() && LidarPoints.size()!=200e6)
       {
           nb_pts++;
          
           auto laser_id = lasreader->point.get_point_source_ID();
           DebugOff("lid "<<laser_id<<endl);
           if(nb_pts%skip!=0){/* Only keep points from Nadir laser */
               continue;
           }
           auto unix_time = lasreader->point.get_gps_time();
           auto x = lasreader->point.get_x();
           auto y = lasreader->point.get_y();
           auto z = lasreader->point.get_z();
           auto uav_x = lasreader->point.get_attribute_as_float(1);
           auto uav_y = lasreader->point.get_attribute_as_float(2);
           auto uav_z = lasreader->point.get_attribute_as_float(3);
           auto roll = lasreader->point.get_attribute_as_float(4);
           auto pitch = lasreader->point.get_attribute_as_float(5);
           auto yaw = lasreader->point.get_attribute_as_float(6);
           for(int i = 0; i < 11; i++){
               DebugOff("attribute "<< i << " = " << lasreader->point.get_attribute_name(i) << endl);
               DebugOff("attribute "<< i << " value = " << lasreader->point.get_attribute_as_float(i) << endl);
           }
           DebugOff("uav "<<to_string_with_precision(uav_x,9)<<" "<<to_string_with_precision(uav_y,9)<<" "<<to_string_with_precision(uav_z,9)<<endl);
           if(!isnan(uav_x) && !isnan(uav_y) && !isnan(uav_z)){
               LidarPoints.push_back(new LidarPoint(laser_id,unix_time,x,y,z));
               point_cloud1.push_back({x,y,z});
               //uav_cloud.push_back({uav_x,uav_y, uav_z});
               //lidar_point_cloud.push_back({x,y,z});
              // roll_pitch_yaw.push_back({roll,pitch,yaw});
           }
          
       }

       DebugOn("Read " << LidarPoints.size() << " points" << endl);
       DebugOn(point_cloud1.size() << " points in flight line 1" << endl);
       DebugOn(point_cloud2.size() << " points in flight line 2" << endl);

 
   }
   
   save_laz(name, point_cloud1, point_cloud2);
   DebugOn("finished read laz"<<endl);
  // return uav_cloud;
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
//vector<vector<double>> turn_detect_slope_map(vector<vector<double>> uav_cloud){
//    vector<vector<double>> u_list;
//    map<double, pair<int, int>> rank_map;
//    vector<double> uturn;
//    const double zero_tol=1e-12;
//    const double tol=1e-2;
//    bool line_found=false;
//    bool turn_found=false;
//    for(auto i=0;i<uav_cloud.size()-10;i++){
//        double p1x=uav_cloud[i][0];
//        double p1y=uav_cloud[i][1];
//        double p2x=uav_cloud[i+10][0];
//        double p2y=uav_cloud[i+10][1];
//
//
//        double sl1=1000, sl2=2000, sl3=3000;
//        if(abs(p2x-p1x)>=zero_tol){
//            sl1=(p2y-p1y)/(p2x-p1x);
//        }
//
//        if(!slope_map.insert(pair<double, pair<int,int>>(sl1, pair<int, int>(i, i+10))).second){
//            auto p1_p2=slope_map[sl1];
//            auto j=p1_p2.first;
//            auto k=p1_p2.second;
//            double p1x_o=uav_cloud[j][0];
//            double p1y_o=uav_cloud[j][1];
//            double p2x_o=uav_cloud[k][0];
//            double p2y_o=uav_cloud[k][1];
//
//            if(((p2x_o-p1x_o)<=-zero_tol && (p2x-p1x)>=zero_tol) || ((p2x_o-p1x_o)>=zero_tol && (p2x-p1x)<=-zero_tol) || ((p2y_o-p1y_o)<=-zero_tol && (p2y-p1y)>=zero_tol) || ((p2y_o-p1y_o)>=zero_tol && (p2y-p1y)<=-zero_tol)){
//                u_list.push_back(uav_cloud[i]);
//                u_list.push_back(uav_cloud[i+10]);
//                u_list.push_back(uav_cloud[j]);
//                u_list.push_back(uav_cloud[k]);
//                DebugOn("a "<<setprecision(15)<<uav_cloud[i][0]<<" "<<setprecision(15)<<uav_cloud[i][1]<<endl);
//                DebugOn("b "<<setprecision(15)<<uav_cloud[i+10][0]<<" "<<setprecision(15)<<uav_cloud[i+10][1]<<endl);
//                DebugOn("c "<<setprecision(15)<<uav_cloud[j][0]<<" "<<setprecision(15)<<uav_cloud[j][1]<<endl);
//                DebugOn("d "<<setprecision(15)<<uav_cloud[k][0]<<" "<<setprecision(15)<<uav_cloud[k][1]<<endl);
//            }
//
//        }
//
//
//    }
//    return u_list;
//}

/*scale uav_cloud and then call this*/
vector<int> reg_slope_lines(vector<vector<double>> uav_cloud, pair<int,int> slice, int turn){
    vector<vector<double>> lines;
    multimap<double, vector<double>, greater <double>> dist_lines;
    multimap<double, pair<int, int>, greater <double>> rank_map;
    vector<int> ulist;
    vector<double> uturn;
    const double zero_tol=1e-12;
    const double tol=1e-2;
    const double slope_zero=0.1;
    bool line1_found=false;
    bool line2_found=false;
    bool line_found=false;
    bool turn_found=false;
    double slope_line=0;
    int line_num=0;
    for(auto i=slice.first+1;i<=slice.second;i++){
        double p1x=uav_cloud[i-1][0];
        double p1y=uav_cloud[i-1][1];
        double p2x=uav_cloud[i][0];
        double p2y=uav_cloud[i][1];
        double p3x=uav_cloud[i+1][0];
        double p3y=uav_cloud[i+1][1];
       
        double sl_prev=0;
        double c_prev=0;
        double sum_xy_prev,sum_xx_prev,sum_x_prev, sum_y_prev, n_prev;
        double begp, endp;
        
        double sl_now=0, c_now=0;
        if(abs(3*(pow(p3x,2)+pow(p2x,2)+pow(p1x,2))-(p3x+p2x+p1x)*(p3x+p2x+p1x))>=zero_tol){
            sl_now=(3*(p3x*p3y+p2x*p2y+p1x*p1y)-(p3x+p2x+p1x)*(p3y+p2y+p1y))/(3*(pow(p3x,2)+pow(p2x,2)+pow(p1x,2))-(p3x+p2x+p1x)*(p3x+p2x+p1x));
            c_now=((p3y+p2y+p1y)-sl_now*(p3x+p2x+p1x))/3;
        }
        else{
            sl_now=1e5;
        }
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
        
        if(lines.empty() ||  (abs(p2y-sl_prev*p2x-c_prev)/abs(p2y)>=1e-1 && abs(p2y-sl_prev*p2x-c_prev)>=1e-1) || (begp<turn && i>=turn)){
            n_prev=2;
            sum_y_prev=p2y+p1y;
            sum_x_prev=p2x+p1x;
            sum_xx_prev=pow(p2x,2)+pow(p1x,2);
            sum_xy_prev=p2x*p2y+p1x*p1y;
            begp=i;
            endp=i+1;
            if(line_found){
                i+=10;
                DebugOn(abs(p2y-sl_prev*p2x-c_prev)/abs(p2y)<<" "<<abs(p2y-sl_prev*p2x-c_prev)<<endl);
            }
            //            if(!lines.empty() && lines.back()[6]==2){
            //                lines.erase(lines.end()-1, lines.end());
            //            }
            line_found=false;
            
            vector<double> res(9);
            lines.push_back(res);
            DebugOff("error "<<abs(p2y-sl_prev*p2x-c_prev)<<endl);
        }
        else{
            sum_y_prev+=p2y;
            sum_x_prev+=p2x;
            sum_xx_prev+=pow(p2x,2);
            sum_xy_prev+=p2x*p2y;
            n_prev++;
            endp=i;
            line_found=true;
            slope_line=sl_prev;
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
        auto sign_y_i=get_sign(y2-y1);
        auto sign_x_i=get_sign(x2-x1);
        double s1=(y2-y1)/(x2-x1);
        if(abs(p1-p2)>10){
            for(auto j=i+1;j<s;j++){
                int p3=lines[j][7];
                int p4=lines[j][8];
                auto x3=uav_cloud[p3][0];
                auto y3=uav_cloud[p3][1];
                auto x4=uav_cloud[p4][0];
                auto y4=uav_cloud[p4][1];
                if(abs(p3-p4)>10){
                    auto sign_y_j=get_sign(y4-y3);
                    auto sign_x_j=get_sign(x4-x3);
                    double s2=(y4-y3)/(x4-x3);
                    if((abs(s1)>=slope_zero && abs(s2)>=slope_zero && abs(sign_y_j-sign_y_i)>=1) || (abs(s1)<=slope_zero && abs(s2)<=slope_zero && abs(sign_x_j-sign_x_i)>=1)){
                        DebugOff("p1 p2 "<<p1 <<" "<<p2<<endl);
                        DebugOff("x1 y1 "<<x1 <<" "<<y1<<endl);
                        DebugOff("x2 y2 "<<x2 <<" "<<y2<<endl);
                        double d1=pow(x2-x1,2)+pow(y2-y1,2);
                        double d2=pow(x3-x4,2)+pow(y3-y4,2);
                        rank_map.insert(pair<double, pair<int,int>>(d1+d2-abs(s2-s1), pair<int,int>(i,j)));
                    }
                }
            }
        }
    }
    if(rank_map.size()>=1){
        auto it=rank_map.begin();
        int l1=it->second.first;
        int l2=it->second.second;
        int p1=lines[l1][7];
        int p2=lines[l1][8];
        int p3=lines[l2][7];
        int p4=lines[l2][8];
//        ulist.push_back(uav_cloud[p1]);
//        ulist.push_back(uav_cloud[p2]);
//        ulist.push_back(uav_cloud[p3]);
//        ulist.push_back(uav_cloud[p4]);
        ulist.push_back(p1);
        ulist.push_back(p2);
        ulist.push_back(p3);
        ulist.push_back(p4);
    }
    //    it++;
    //    int p5=it->second[7];
    //    int p6=it->second[8];
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


vector<pair<int, int>> extract_slices(vector<vector<double>> uav_cloud, vector<int> index_set, vector<int>& turns){
    vector<vector<vector<double>>> slice;
    vector<pair<int, int>> res;
    const double slope_zero=0.1;
    const double zero_tol=1e-12;
    const double angle_max=50;
    const double tan_angle_max=tan(30*pi/180);
    bool sign_change1=false;
    bool sign_change2=false;
    double sl_now=0, c_now=0, slx, sly;
    double sl_prev=0;
    double c_prev=0;
    int endp=0;
    int begp=0;
    int turn1=0, turn2=0;
    int i=0;
    begp=1;
    while((index_set.size()-endp)>=10){
        vector<double> q(0);
        vector<vector<double>> slice_i(0);
        // slice_i.push_back(q);
        slice.push_back(slice_i);
        std::pair<int,int> pair1;
        pair1.first=index_set[begp];
        pair1.second=index_set[begp];
        res.push_back(pair1);
        
        int sign_x_prev=0, sign_x_now=0;
        int sign_y_prev=0, sign_y_now=0;
        
        sign_change1=false;
        sign_change2=false;
        bool sign_compute1=false;
        bool sign_compute2=false;
        for(i=begp;i<index_set.size()-1;i++){
            double p1x=uav_cloud[index_set[i-1]][0];
            double p1y=uav_cloud[index_set[i-1]][1];
            double p1z=uav_cloud[index_set[i-1]][2];
            double p2x=uav_cloud[index_set[i]][0];
            double p2y=uav_cloud[index_set[i]][1];
            double p2z=uav_cloud[index_set[i]][2];
            double p3x=uav_cloud[index_set[i+1]][0];
            double p3y=uav_cloud[index_set[i+1]][1];
            double p3z=uav_cloud[index_set[i+1]][2];
            
            if(abs(3*(pow(p3x,2)+pow(p2x,2)+pow(p1x,2))-(p3x+p2x+p1x)*(p3x+p2x+p1x))>=zero_tol){
                sl_now=(3*(p3x*p3y+p2x*p2y+p1x*p1y)-(p3x+p2x+p1x)*(p3y+p2y+p1y))/(3*(pow(p3x,2)+pow(p2x,2)+pow(p1x,2))-(p3x+p2x+p1x)*(p3x+p2x+p1x));
                c_now=((p3y+p2y+p1y)-sl_now*(p3x+p2x+p1x))/3;
            }
            else{
                sl_now=1e5;
            }
            if((abs(p2y-p1y)<=1e-3 && abs(p2x-p1x)<=1e-3)){
                sign_x_now=sign_x_prev;
                sign_y_now=sign_y_prev;
            }
            else{
                if(abs(sl_now)>=slope_zero && abs(sl_now)<=1e5-slope_zero){
                    if(abs(p2y-p1y)<=1e-3){
                        sign_y_now=sign_y_prev;
                    }
                    else{
                        sign_y_now=get_sign(p2y-p1y);
                        if(!sign_compute2 && sign_compute1){
                            sign_compute2=true;
                        }
                        if(!sign_compute1){
                            sign_compute1=true;
                        }
                    }
                    sign_x_now=sign_y_now;
                }
                else if(abs(sl_now)<=slope_zero){
                    if(abs(p2y-p1y)<=1e-3){
                        sign_y_now=sign_y_prev;
                    }
                    else{
                        sign_y_now=get_sign(p2y-p1y);
                        if(!sign_compute2 && sign_compute1){
                            sign_compute2=true;
                        }
                        if(!sign_compute1){
                            sign_compute1=true;
                        }
                    }
                    sign_x_now=sign_y_now;
                }
                else{
                    if(abs(p2x-p1x)<=1e-3){
                        sign_x_now=sign_x_prev;
                    }
                    else{
                        sign_x_now=get_sign(p2x-p1x);
                        if(!sign_compute2 && sign_compute1){
                            sign_compute2=true;
                        }
                        if(!sign_compute1){
                            sign_compute1=true;
                        }
                    }
                    sign_y_now=sign_x_now;
                }
                
                
            }
            if(abs(sign_x_now-sign_x_prev)>=1 && abs(sign_y_now-sign_y_prev)>=1 && sign_compute2){
                if(!sign_change1){
                    sign_change1=true;
                    turn1=index_set[i];
                    turns.push_back(turn1);
                    begp=i+10;
                    DebugOff("turn1 "<<i<<" "<<sign_y_now<<" "<<sign_y_prev<<endl);
                    DebugOff(p2y<<" "<<p1y<<endl);
                    DebugOff(p2x<<" "<<p1x<<endl);
                }
                else{
                    if(i>=begp){
                        sign_change2=true;
                        turn2=index_set[i];
                        endp=i;
                        DebugOff("turn2 "<<i<<" "<<sign_y_now<<" "<<sign_y_prev<<endl);
                        DebugOff(p2y<<" "<<p1y<<endl);
                        DebugOff(p2x<<" "<<p1x<<endl);
                    }
                }
            }
            
            sign_x_prev=sign_x_now;
            sign_y_prev=sign_y_now;
            
            if(!sign_change2){
                slice.back().push_back(uav_cloud[index_set[i]]);
            }
            else
            {
                res.back().second=turn2;
                break;
            }
        }
        if(i>=index_set.size()-10){
            res.back().second=index_set[i];
            break;
        }
    }
    for(auto i=0;i<slice.size();i++){
        vector<vector<double>> empty_vec;
        empty_vec.push_back(slice[i][0]);
        
        //plot(slice[i], empty_vec, 1);
    }
    return res;
}

/*Remove points with either atan(dz/dy) or atan(dz/dx) greater than +- 30*/
vector<int> filter_z_slope(vector<vector<double>> uav_cloud){
    vector<vector<double>> res;
    vector<int> index_set;
    const double angle_max=30;
    const double zero_tol=1e-12;
    const double tan_angle_max=tan(30*pi/180);
    double slx, sly;
    for(auto i=1;i<uav_cloud.size()-1;i++){
        auto p1x=uav_cloud.at(i-1)[0];
        auto p1y=uav_cloud.at(i-1)[1];
        auto p1z=uav_cloud.at(i-1)[2];
        auto p2x=uav_cloud.at(i)[0];
        auto p2y=uav_cloud.at(i)[1];
        auto p2z=uav_cloud.at(i)[2];
        auto p3x=uav_cloud.at(i+1)[0];
        auto p3y=uav_cloud.at(i+1)[1];
        auto p3z=uav_cloud.at(i+1)[2];
        if(abs(3*(pow(p3x,2)+pow(p2x,2)+pow(p1x,2))-(p3x+p2x+p1x)*(p3x+p2x+p1x))>=zero_tol){
            slx=(3*(p3x*p3z+p2x*p2z+p1x*p1z)-(p3x+p2x+p1x)*(p3z+p2z+p1z))/(3*(pow(p3x,2)+pow(p2x,2)+pow(p1x,2))-(p3x+p2x+p1x)*(p3x+p2x+p1x));
            //c_now=((p3y+p2y+p1y)-sl_now*(p3x+p2x+p1x))/3;
        }
        else{
            slx=1e5;
        }
        if(abs(3*(pow(p3y,2)+pow(p2y,2)+pow(p1y,2))-(p3y+p2y+p1y)*(p3y+p2y+p1y))>=zero_tol){
            sly=(3*(p3y*p3z+p2y*p2z+p1y*p1z)-(p3y+p2y+p1y)*(p3z+p2z+p1z))/(3*(pow(p3y,2)+pow(p2y,2)+pow(p1y,2))-(p3y+p2y+p1y)*(p3y+p2y+p1y));
            //c_now=((p3y+p2y+p1y)-sl_now*(p3x+p2x+p1x))/3;
        }
        else{
            sly=1e5;
        }
        
        //        double slx=(z3-z1)/(x3-x1);
        //        double sly=(z3-z1)/(y3-y1);
        if((slx<tan_angle_max) && (slx>-tan_angle_max) && (sly <tan_angle_max) && (sly >-tan_angle_max)){
            res.push_back(uav_cloud.at(i));
            index_set.push_back(i);
        }
    }
    vector<vector<double>> empty_vec;
    empty_vec.push_back(res[0]);
   // plot(res, empty_vec,1);
    return index_set;
}
void fit_points_line(vector<vector<double>> uav_cloud, int start, int stop, double& sl, double& c){
    double sum_xx=0, sum_x=0, sum_y=0, sum_xy=0;
    int n=0;
    const double zero_tol=1e-12;
    for(auto i=start;i<stop; i++){
        auto x=uav_cloud.at(i)[0];
        auto y=uav_cloud.at(i)[1];
        auto z=uav_cloud.at(i)[2];
        n++;
        sum_xx+=x*x;
        sum_x+=x;
        sum_y+=y;
        sum_xy+=x*y;
    }
        if(abs(n*sum_xx-sum_x*sum_x)>=zero_tol){
            sl=(n*sum_xy-sum_x*sum_y)/(n*sum_xx-sum_x*sum_x);
            c=(sum_y-sl*sum_x)/n;
        }
        else{
            sl=1e5;
            c=sum_x/n;
    }
}

void get_frame(vector<vector<double>> uav_cloud, int start1, int stop1, int start2, int stop2, vector<int>& frame1_index, vector<int>& frame2_index){
    vector<vector<double>> frame1, frame2;
    double sum_xx=0, sum_x=0, sum_y=0, sum_xy=0;
    int n=0;
    const double zero_tol=1e-12;
    double sl1=0, c1=0, sl2=0, c2=0, rsl1=0, rsl2=0, rc1=0, rc2=0;
    double px1=uav_cloud.at(start1)[0];
    double py1=uav_cloud.at(start1)[1];
    double pz1=uav_cloud.at(start1)[2];
    double px2=uav_cloud.at(stop1)[0];
    double py2=uav_cloud.at(stop1)[1];
    double pz2=uav_cloud.at(stop1)[2];
    double px3=uav_cloud.at(start2)[0];
    double py3=uav_cloud.at(start2)[1];
    double pz3=uav_cloud.at(start2)[2];
    double px4=uav_cloud.at(stop2)[0];
    double py4=uav_cloud.at(stop2)[1];
    double pz4=uav_cloud.at(stop2)[2];
    if((pow(px2-px1, 2)+pow(py2-py1, 2))>=(pow(px3-px4, 2)+pow(py3-py4, 2))){
        auto t1=start2;
        auto t2=stop2;
        start2=start1;
        stop2=stop1;
        start1=t1;
        stop1=t2;
    
    px1=uav_cloud.at(start1)[0];
    py1=uav_cloud.at(start1)[1];
    pz1=uav_cloud.at(start1)[2];
    px2=uav_cloud.at(stop1)[0];
    py2=uav_cloud.at(stop1)[1];
    pz2=uav_cloud.at(stop1)[2];
    px3=uav_cloud.at(start2)[0];
    py3=uav_cloud.at(start2)[1];
    pz3=uav_cloud.at(start2)[2];
    px4=uav_cloud.at(stop2)[0];
    py4=uav_cloud.at(stop2)[1];
    pz4=uav_cloud.at(stop2)[2];
    }
    double mid_x1=(px1*0.5+px2*0.5);
    double mid_x2=(px1*0.4+px2*0.6);
    double mid=(mid_x1+mid_x2)*0.5;
   
    //double mid_x2=(px3+px4)*0.5;
   
  
    
    double xmin1=1e9, xmax1=-1e9, ymin1=1e9, ymax1=-1e9, zmin1=1e9, zmax1=-1e9;
    int start_line, stop_line;
    bool con_frame=false;
    for(auto i=start1;i<stop1;i++){
        auto x=uav_cloud.at(i)[0];
        auto y=uav_cloud.at(i)[1];
        auto z=uav_cloud[start1][2];
        double lamda=(x-px2)/(px1-px2);
        if(lamda>=0.4 && lamda<=0.5){
            if(!con_frame){
                start_line=i;
            }
            frame1.push_back({x,y,z});
            frame1_index.push_back(i);
            if(x<=xmin1){
                xmin1=x;
            }
            if(x>=xmax1){
                xmax1=x;
            }
            if(y<=ymin1){
                ymin1=y;
            }
            if(y>=ymax1){
                ymax1=y;
            }
            con_frame=true;
        }
        else{
            if(con_frame){
                stop_line=i;
            }
            con_frame=false;
        }
    }
    DebugOn("mid "<<mid_x1<<" "<<mid_x2<<" "<<mid<<endl);
    DebugOn("line start "<<uav_cloud[start_line][0]<<" "<<uav_cloud[stop_line][0]<<" "<<mid<<endl);
    
    fit_points_line(uav_cloud, start_line, stop_line, sl1, c1);
    vector<vector<double>> line;
    for(auto i=start_line;i<stop_line;i++)
    {
        auto x=uav_cloud.at(i)[0];
        auto y=sl1*x+c1;
        auto z=uav_cloud[start1][2];
        line.push_back({x,y,z});
    }
    double midy=sl1*mid+c1;
    rsl1=-1.0/sl1;
    double mid_y1=sl1*mid_x1+c1;
    double mid_y2=sl1*mid_x2+c1;
    rc1=mid_y1-rsl1*mid_x1;
    rc2=mid_y2-rsl1*mid_x2;
    
    
    for(auto i=-10;i<10;i++)
    {
        auto x=mid_x1-i;
        auto y=rsl1*x+rc1;
        auto z=uav_cloud[start1][2];
        line.push_back({x,y,z});
    }
    for(auto i=-10;i<10;i++)
    {
        auto x=mid_x2-i;
        auto y=rsl1*x+rc2;
        auto z=uav_cloud[start1][2];
        line.push_back({x,y,z});
    }
    vector<vector<double>> pl1, pl2;
    
    pl1.push_back(uav_cloud[start1]);
    pl1.push_back(uav_cloud[start2]);
    pl1.push_back(uav_cloud[stop1]);
    pl1.push_back(uav_cloud[stop2]);
    
    vector<double> v(3);
    v[0]=mid_x1;
    v[1]=mid_y1;
    v[2]=uav_cloud[start1][2];
    pl2.push_back(v);
    v.clear();
    v.resize(3);
    v[0]=mid_x2;
    v[1]=mid_y2;
    v[2]=uav_cloud[start1][2];
    pl2.push_back(v);
    v.clear();
    v.resize(3);
    v[0]=mid;
    v[1]=midy;
    v[2]=uav_cloud[start1][2];
    pl2.push_back(v);
    
    
    
    //plot(pl1, line, 1);
    
    auto sign_mid_line1=get_sign(midy-rsl1*mid-rc1);
    auto sign_mid_line2=get_sign(midy-rsl1*mid-rc2);
    double xmin2=1e9, xmax2=-1e9, ymin2=1e9, ymax2=-1e9;
    DebugOn("sign mid "<<sign_mid_line1<<endl);
    DebugOn("sign mid "<<sign_mid_line2<<endl);
    for(auto i=start2;i<stop2;i++){
        auto x=uav_cloud.at(i)[0];
        auto y=uav_cloud.at(i)[1];
        auto z=uav_cloud[start1][2];
      
        auto sign1=get_sign(y-rsl1*x-rc1);
        auto sign2=get_sign(y-rsl1*x-rc2);
        DebugOff("sign1 "<<sign1<<" sign2 "<<sign2<<endl);
        if(abs(sign1-sign_mid_line1)==0 && abs(sign2-sign_mid_line2)==0){
            if(x<=xmin2){
                xmin2=x;
            }
            if(x>=xmax2){
                xmax2=x;
            }
            if(y<=ymin2){
                ymin2=y;
            }
            if(y>=ymax2){
                ymax2=y;
            }
            frame2.push_back({x,y,z});
            frame2_index.push_back(i);
        }
    }
    DebugOn("frame 1 size "<<frame1.size()<<endl);
    DebugOn("xmin "<<xmin1<<" "<<xmax1<<endl);
    DebugOn("ymin "<<ymin1<<" "<<ymax1<<endl);
    DebugOn("frame 2 size "<<frame2.size()<<endl);
    DebugOn("xmin "<<xmin2<<" "<<xmax2<<endl);
    DebugOn("ymin "<<ymin2<<" "<<ymax2<<endl);
    //plot(frame1, frame2, 1);
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
