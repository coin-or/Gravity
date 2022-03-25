//
//  Lidar_utils.h
//  Gravity
//
//  Created by Smitha on 3/21/22.
//

#ifndef Lidar_utils_h
#define Lidar_utils_h
/* Apply rotation + translation on input data (using rotation +translation matrix) */
void apply_rot_trans(const vector<double>& theta_matrix, vector<vector<double>>& point_cloud);

/* Apply rotation + translation on input data (using 3 angles) */
void apply_rot_trans(double roll, double pitch, double yaw, double x_shift, double y_shift, double z_shift, vector<vector<double>>& point_cloud);

void apply_rot_trans_util(const vector<double>& theta_matrix, vector<vector<double>>& point_cloud){
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
double computeL2error_util(const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data, double best_ub);
/* Compute the L1 error for model and data sets
 @param[in] point_cloud_model, Model point cloud
 @param[in] point_cloud_data, Data point cloud
 @param[in] matching, vecor used to store the optimal match
 @param[in] err_per_point, vector used to store the minimum L1 error per data point
 */
double computeL2error_util(const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data, double best_ub){
    size_t n = point_cloud_data.size();
    size_t m = point_cloud_model.size();
    double dist_sq = 0, err = 0;
    for (auto i = 0; i< n; i++) {
        double min_dist = best_ub;
        for (auto j = 0; j< m; j++) {
            dist_sq = std::pow(point_cloud_data.at(i).at(0) - point_cloud_model.at(j).at(0),2) + std::pow(point_cloud_data.at(i).at(1) - point_cloud_model.at(j).at(1),2) + std::pow(point_cloud_data.at(i).at(2) - point_cloud_model.at(j).at(2),2);
            if(min_dist>dist_sq){
                min_dist = dist_sq;
            }
        }
        err += min_dist;
        if(err>=best_ub){
            err=best_ub+10;
            break;
        }
    }
    return err;
}
double computeL2error_util(const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data, double best_ub, double& tx,double& ty, double& tz ){
    size_t n = point_cloud_data.size();
    size_t m = point_cloud_model.size();
    double dist_sq = 0, err = 0, txa,tya,tza;
    for (auto i = 0; i< n; i++) {
        double min_dist = best_ub;
        for (auto j = 0; j< m; j++) {
            dist_sq = std::pow(point_cloud_data.at(i).at(0) - point_cloud_model.at(j).at(0),2) + std::pow(point_cloud_data.at(i).at(1) - point_cloud_model.at(j).at(1),2) + std::pow(point_cloud_data.at(i).at(2) - point_cloud_model.at(j).at(2),2);
            if(min_dist>dist_sq){
                min_dist = dist_sq;
                txa=point_cloud_model.at(j).at(0);
                tya=point_cloud_model.at(j).at(1);
                tza=point_cloud_model.at(j).at(2);
            }
        }
        err += min_dist;
        tx+=txa;
        ty+=tya;
        tz+=tza;
    }
    tx/=n;
    ty/=n;
    tz/=n;
    return err;
}
double computeL2error(const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data, vector<int>& matching, vector<double>& err_per_point);

/* Compute the L1 error */
double computeL1error(const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data, vector<int>& matching, vector<double>& err_per_point);

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


#endif /* Lidar_utils_h */
