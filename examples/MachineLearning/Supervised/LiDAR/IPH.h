//
//  IPH.h
//  Gravity
//
//  Created by Smitha on 8/18/21.
//

#ifndef IPH_h
#define IPH_h

/* Return the min-max values for x, y and z  for all possible rotations of p with angle +- angle*/
vector<pair<double,double>> get_min_max(double angle, const vector<double>& p, const vector<double>& ref);
/* Run the iterative projection heuristic on the Registration problem */
tuple<double,double,double,double,double,double> run_IPH(const vector<vector<double>>& ext_model, vector<vector<double>>& ext_data, vector<vector<double>>& point_cloud_data);
/* Run the ARMO model for registration */
tuple<double,double,double,double,double,double> run_ARMO(bool bypass, string axis, const vector<vector<double>>& point_cloud1, const vector<vector<double>>& point_cloud2);
/* Run the ARMO model for boresight alignment */
tuple<double,double,double> run_ARMO(string axis, const vector<vector<double>>& point_cloud1, const vector<vector<double>>& point_cloud2, const vector<vector<double>>& uav1, const vector<vector<double>>& uav2, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max);
tuple<double,double,double,double,double,double> run_ARMO(bool bypass, string axis, const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data)
{
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
            Reg.add(yaw.in(range(0,0)),pitch.in(range(0,0)),roll.in(range(0,0)));
            Reg.add(x_shift.in(range(0,0)),y_shift.in(range(0,0)),z_shift.in(range(0,0)));
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
tuple<double,double,double> run_NLP(string axis, const vector<vector<double>>& point_cloud1, const vector<vector<double>>& point_cloud2, const vector<vector<double>>& uav1, const vector<vector<double>>& uav2, const vector<vector<double>>& rollpitchyawins1, const vector<vector<double>>& rollpitchyawins2)
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
        param<> x1i("x1i"), x2i("x2i"), y1i("y1i"), y2i("y2i"), z1i("z1i"), z2i("z2i");
        param<> x_uav1("x_uav1"), y_uav1("y_uav1"), z_uav1("z_uav1");
        param<> x_uav2("x_uav2"), y_uav2("y_uav2"), z_uav2("z_uav2");
        param<> roll_ins1("roll_ins1"),pitch_ins1("pitch_ins1"), yaw_ins1("yaw_ins1");
        param<> roll_ins2("roll_ins2"), pitch_ins2("pitch_ins2"), yaw_ins2("yaw_ins2");
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
                roll_ins1.add_val(i_str, rollpitchyawins1.at(i)[0]);
                pitch_ins1.add_val(i_str, (-pi+rollpitchyawins1.at(i)[1])*(-1));
                yaw_ins1.add_val(i_str, (-pi/2+rollpitchyawins1.at(i)[2])*(-1));
                auto res_d=apply_rotation_new_order(rollpitchyawins1.at(i)[0], (-pi+rollpitchyawins1.at(i)[1])*(-1), (-pi/2+rollpitchyawins1.at(i)[2])*(-1), point_cloud1.at(i)[0]-uav1.at(i)[0], point_cloud1.at(i)[1]-uav1.at(i)[1], point_cloud1.at(i)[2]-uav1.at(i)[2]);
                x1i.add_val(i_str,res_d[0]);
                y1i.add_val(i_str,res_d[1]);
                z1i.add_val(i_str,res_d[2]);
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
                        roll_ins2.add_val(j_str, rollpitchyawins2.at(k)[0]);
                        pitch_ins2.add_val(j_str, (-pi+rollpitchyawins2.at(k)[1])*(-1));
                        yaw_ins2.add_val(j_str, (-pi/2+rollpitchyawins2.at(k)[2])*(-1));
                        auto res_m=apply_rotation_new_order(rollpitchyawins2.at(k)[0], (-pi+rollpitchyawins2.at(k)[1])*(-1), (-pi/2+rollpitchyawins2.at(k)[2])*(-1), point_cloud2.at(k)[0]-uav2.at(k)[0], point_cloud2.at(k)[1]-uav2.at(k)[1], point_cloud2.at(k)[2]-uav2.at(k)[2]);
                        x2i.add_val(j_str,res_m[0]);
                        y2i.add_val(j_str,res_m[1]);
                        z2i.add_val(j_str,res_m[2]);
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

        
        
        
        if (solve_lidar_iter) {
            Model<> Lidar("Lidar");
            var<> new_x1("new_x1"), new_y1("new_y1"), new_z1("new_z1");
            var<> new_x2("new_x2"), new_y2("new_y2"), new_z2("new_z2");
            var<> xb_d("xb_d"), yb_d("yb_d"), zb_d("zb_d");
            var<> xb_m("xb_m"), yb_m("yb_m"), zb_m("zb_m");
            var<> x_diff("x_diff", pos_), y_diff("y_diff", pos_), z_diff("z_diff", pos_);
    
            var<> yaw("yaw", -0.00872, 0.00872), pitch("pitch",  -0.00872, 0.00872), roll("roll",  -0.00872, 0.00872);
            
            Lidar.add(yaw.in(range(0,0)),pitch.in(range(0,0)),roll.in(range(0,0)));
            Lidar.add(new_x1.in(N1), new_y1.in(N1), new_z1.in(N1));
            Lidar.add(new_x2.in(N2), new_y2.in(N2), new_z2.in(N2));
            Lidar.add(xb_d.in(N1), yb_d.in(N1), zb_d.in(N1));
            Lidar.add(xb_m.in(N2), yb_m.in(N2), zb_m.in(N2));
            Lidar.add(x_diff.in(cells), y_diff.in(cells), z_diff.in(cells));
            //                Lidar.add(z_diff.in(cells));
            
            
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
            
            
           
            
            func<> r11 = cos(roll)*cos(yaw);r11.eval_all();
            func<> r12 = (-1)*cos(roll)*sin(yaw);r12.eval_all();
            func<> r13 =sin(roll);r13.eval_all();
            
            func<> r21 =cos(pitch)*sin(yaw)+cos(yaw)*sin(roll)*sin(pitch) ;r21.eval_all();
            func<> r22 = cos(pitch)*cos(yaw)-sin(roll)*sin(pitch)*sin(yaw) ;r22.eval_all();
            func<> r23 = (-1)*cos(roll)*sin(pitch);r23.eval_all();
            
            func<> r31 = sin(pitch)*sin(yaw)-cos(pitch)*cos(yaw)*sin(roll);r31.eval_all();
            func<> r32 = cos(yaw)*sin(pitch)+cos(pitch)*sin(roll)*sin(yaw);r32.eval_all();
            func<> r33 = cos(roll)*cos(pitch);r33.eval_all();
            
            var<> theta11("theta11",  std::max(-1.,r11._range->first), std::min(1.,r11._range->second)), theta12("theta12", std::max(-1.,r12._range->first), std::min(1.,r12._range->second)), theta13("theta13", std::max(-1.,r13._range->first), std::min(1.,r13._range->second));
            var<> theta21("theta21", std::max(-1.,r21._range->first), std::min(1.,r21._range->second)), theta22("theta22", std::max(-1.,r22._range->first), std::min(1.,r22._range->second)), theta23("theta23", std::max(-1.,r23._range->first), std::min(1.,r23._range->second));
            var<> theta31("theta31", std::max(-1.,r31._range->first), std::min(1.,r31._range->second)), theta32("theta32", std::max(-1.,r32._range->first), std::min(1.,r32._range->second)), theta33("theta33", std::max(-1.,r33._range->first), std::min(1.,r33._range->second));
            
            
            Lidar.add(theta11.in(range(0,0)),theta12.in(range(0,0)),theta13.in(range(0,0)));
            Lidar.add(theta21.in(range(0,0)),theta22.in(range(0,0)),theta23.in(range(0,0)));
            Lidar.add(theta31.in(range(0,0)),theta32.in(range(0,0)),theta33.in(range(0,0)));
            
            
            
            
            Constraint<> T11("T11");
            T11=theta11.in(range(0,0));
            T11-=cos(roll.in(range(0,0)))*cos(yaw.in(range(0,0)));
            Lidar.add(T11.in(range(0,0))==0);
            
            Constraint<> T12("T12");
            T12=theta12.in(range(0,0));
            T12-=(-1)*cos(roll.in(range(0,0)))*sin(yaw.in(range(0,0)));
            Lidar.add(T12.in(range(0,0))==0);
            
            Constraint<> T13("T13");
            T13=theta13.in(range(0,0));
            T13-=sin(roll.in(range(0,0)));
            Lidar.add(T13.in(range(0,0))==0);
            
            Constraint<> T21("T21");
            T21+=theta21.in(range(0,0));
            T21-=cos(pitch.in(range(0,0)))*sin(yaw.in(range(0,0)));
            T21-=cos(yaw.in(range(0,0)))*sin(roll.in(range(0,0)))*sin(pitch.in(range(0,0)));
            Lidar.add(T21.in(range(0,0))==0);
            
            Constraint<> T22("T22");
            T22+=theta22.in(range(0,0));
            T22-=cos(pitch.in(range(0,0)))*cos(yaw.in(range(0,0)));
            T22-=(-1)*sin(roll.in(range(0,0)))*sin(pitch.in(range(0,0)))*sin(yaw.in(range(0,0)));
            Lidar.add(T22.in(range(0,0))==0);
            
            Constraint<> T23("T23");
            T23+=theta23.in(range(0,0));
            T23-=(-1)*cos(roll.in(range(0,0)))*sin(pitch.in(range(0,0)));
            Lidar.add(T23.in(range(0,0))==0);
        
            
            Constraint<> T31("T31");
            T31+=theta31.in(range(0,0));
            T31-=sin(pitch.in(range(0,0)))*sin(yaw.in(range(0,0)));
            T31-=(-1)*cos(pitch.in(range(0,0)))*cos(yaw.in(range(0,0)))*sin(roll.in(range(0,0)));
            Lidar.add(T31.in(range(0,0))==0);
            
            Constraint<> T32("T32");
            T32+=theta32.in(range(0,0));
            T32-=cos(yaw.in(range(0,0)))*sin(pitch.in(range(0,0)));
            T32-=cos(pitch.in(range(0,0)))*sin(roll.in(range(0,0)))*sin(yaw.in(range(0,0)));
            Lidar.add(T32.in(range(0,0))==0);
            
            Constraint<> T33("T33");
            T33+=theta33.in(range(0,0));
            T33-=cos(roll.in(range(0,0)))*cos(pitch.in(range(0,0)));
            Lidar.add(T33.in(range(0,0))==0);
            
            auto ids1 = theta11.repeat_id(N1.size());
            auto ids2 = theta11.repeat_id(N2.size());
            
            
            Constraint<> xd_bore("xd_bore");
            xd_bore=x1i.in(N1)*theta11.in(ids1)+y1i.in(N1)*theta12.in(ids1)+z1i.in(N1)*theta13.in(ids1)-xb_d.in(N1);
            Lidar.add(xd_bore.in(N1)==0);
            
            Constraint<> xm_bore("xm_bore");
            xm_bore=x2i.in(N2)*theta11.in(ids2)+y2i.in(N2)*theta12.in(ids2)+z2i.in(N2)*theta13.in(ids2)-xb_m.in(N2);
            Lidar.add(xm_bore.in(N2)==0);
            
            Constraint<> yd_bore("yd_bore");
            yd_bore=x1i.in(N1)*theta21.in(ids1)+y1i.in(N1)*theta22.in(ids1)+z1i.in(N1)*theta23.in(ids1)-yb_d.in(N1);
            Lidar.add(yd_bore.in(N1)==0);
            
            Constraint<> ym_bore("ym_bore");
            ym_bore=x2i.in(N2)*theta21.in(ids2)+y2i.in(N2)*theta22.in(ids2)+z2i.in(N2)*theta23.in(ids2)-yb_m.in(N2);
            Lidar.add(ym_bore.in(N2)==0);
            
            Constraint<> zd_bore("zd_bore");
            zd_bore=x1i.in(N1)*theta31.in(ids1)+y1i.in(N1)*theta32.in(ids1)+z1i.in(N1)*theta33.in(ids1)-zb_d.in(N1);
            Lidar.add(zd_bore.in(N1)==0);
            
            Constraint<> zm_bore("zm_bore");
            zm_bore=x2i.in(N2)*theta31.in(ids2)+y2i.in(N2)*theta32.in(ids2)+z2i.in(N2)*theta33.in(ids2)-zb_m.in(N2);
            Lidar.add(zm_bore.in(N2)==0);
            
            
            
//            Constraint<> xm_bore("xm_bore");
//            xm_bore=x2*cos(roll.in(ids2))*cos(yaw.in(ids2))-y2*cos(roll.in(ids2))*sin(yaw.in(ids2))+z2*sin(roll.in(ids2))-xb_m;
//            Lidar.add(xm_bore.in(N2)==0);
//
//            Constraint<> yd_bore("yd_bore");
//            yd_bore=x1*(cos(pitch.in(ids1))*sin(yaw.in(ids1)) +cos(yaw.in(ids1))*sin(roll.in(ids1))*sin(pitch.in(ids1)))+y1*(cos(pitch.in(ids1))*cos(yaw.in(ids1))-sin(roll.in(ids1))*sin(pitch.in(ids1))*sin(yaw.in(ids1)))+z1*((-1)*cos(roll.in(ids1))*sin(pitch.in(ids1)))-yb_d;
//            Lidar.add(yd_bore.in(N1)==0);
//
//
//            Constraint<> ym_bore("ym_bore");
//            ym_bore=x2*(cos(pitch.in(ids2))*sin(yaw.in(ids2)) +cos(yaw.in(ids2))*sin(roll.in(ids2))*sin(pitch.in(ids2)))+y2*(cos(pitch.in(ids2))*cos(yaw.in(ids2))-sin(roll.in(ids2))*sin(pitch.in(ids2))*sin(yaw.in(ids2)))+z2*((-1)*cos(roll.in(ids2))*sin(pitch.in(ids2)))-yb_m;
//            Lidar.add(ym_bore.in(N2)==0);
//
//            Constraint<> zd_bore("zd_bore");
//            zd_bore=x1.in(N1)*(sin(pitch.in(ids1))*sin(yaw.in(ids1))-cos(pitch.in(ids1))*cos(yaw.in(ids1))*sin(roll).in(ids1))+y1.in(N1)*(cos(yaw.in(ids1))*sin(pitch.in(ids1))+cos(pitch.in(ids1))*sin(roll.in(ids1))*sin(yaw.in(ids1)))+z1.in(N1)*cos(roll.in(ids1))*cos(pitch.in(ids1))-zb_d;
//            Lidar.add(zd_bore.in(N1)==0);
//
//            Constraint<> zm_bore("zm_bore");
//            zm_bore=x2*(sin(pitch.in(ids2))*sin(yaw.in(ids2))-cos(pitch.in(ids2))*cos(yaw.in(ids2))*sin(roll).in(ids2))+y2*(cos(yaw.in(ids2))*sin(pitch.in(ids2))+cos(pitch.in(ids2))*sin(roll.in(ids2))*sin(yaw.in(ids2)))+z2*(cos(roll.in(ids2))*cos(pitch.in(ids2)))-zb_m;
//            Lidar.add(zm_bore.in(N2)==0);
            
            
            
           
            /* alpha = yaw_, beta = pitch_ and gamma = roll_ */
            if(axis!="x"){
                
Constraint<> x_rot1("x_rot1");
x_rot1 += new_x1 - x_uav1.in(N1);
x_rot1 -= xb_d*cos(roll_ins1)*cos(yaw_ins1)+yb_d*(cos(pitch_ins1)*sin(yaw_ins1) +cos(yaw_ins1)*sin(roll_ins1)*sin(pitch_ins1))+zb_d*(sin(pitch_ins1)*sin(yaw_ins1)-cos(pitch_ins1)*cos(yaw_ins1)*sin(roll_ins1));
Lidar.add(x_rot1.in(N1)==0);
                
Constraint<> x_rot2("x_rot2");
x_rot2 += new_x2 - x_uav2.in(N2);
x_rot2 -= xb_m*cos(roll_ins2)*cos(yaw_ins2)+yb_m*(cos(pitch_ins2)*sin(yaw_ins2) +cos(yaw_ins2)*sin(roll_ins2)*sin(pitch_ins2))+zb_m*(sin(pitch_ins2)*sin(yaw_ins2)-cos(pitch_ins2)*cos(yaw_ins2)*sin(roll_ins2));
Lidar.add(x_rot2.in(N2)==0);
            }
            
            if(axis!="y"){
                
                Constraint<> y_rot1("y_rot1");
                y_rot1 += new_y1 - y_uav1.in(N1);
                y_rot1 -= xb_d*(-1)*cos(roll_ins1)*sin(yaw_ins1)+yb_d*(cos(pitch_ins1)*cos(yaw_ins1) -sin(roll_ins1)*sin(pitch_ins1)*sin(yaw_ins1))+zb_d*(cos(yaw_ins1)*sin(pitch_ins1)+cos(pitch_ins1)*sin(roll_ins1)*sin(yaw_ins1));
                Lidar.add(y_rot1.in(N1)==0);

                
                Constraint<> y_rot2("y_rot2");
                y_rot2 += new_y2 - y_uav2.in(N2);
                y_rot2 -= xb_m*(-1)*cos(roll_ins2)*sin(yaw_ins2)+yb_m*(cos(pitch_ins2)*cos(yaw_ins2) -sin(roll_ins2)*sin(pitch_ins2)*sin(yaw_ins2))+zb_m*(cos(yaw_ins2)*sin(pitch_ins2)+cos(pitch_ins2)*sin(roll_ins2)*sin(yaw_ins2));
                Lidar.add(y_rot2.in(N2)==0);
            }
            
            if(axis!="z"){
                
                Constraint<> z_rot1("z_rot1");
                z_rot1 += new_z1 - z_uav1.in(N1);
                z_rot1 -= xb_d*sin(roll_ins1)+yb_d*(-1)*cos(roll_ins1)*sin(pitch_ins1)+zb_d*(cos(roll_ins1)*cos(pitch_ins1));
                Lidar.add(z_rot1.in(N1)==0);
                
                Constraint<> z_rot2("z_rot2");
                z_rot2 += new_z2 - z_uav2.in(N2);
                z_rot2 -= xb_m*sin(roll_ins2)+yb_m*(-1)*cos(roll_ins2)*sin(pitch_ins2)+zb_m*(cos(roll_ins2)*cos(pitch_ins2));
                Lidar.add(z_rot2.in(N2)==0);
            }
            
            //    M.min(sum(z_diff)/nb_overlap);
            
            //        M.min(sum(z_diff));
            if(axis == "full")
                Lidar.min(sum(x_diff)/cells.size() + sum(y_diff)/cells.size() + sum(z_diff)/cells.size());
            else if(axis == "x"){
                Lidar.min(sum(y_diff)/cells.size() + sum(z_diff)/cells.size());
                roll.set_lb(0);
                roll.set_ub(0);
                yaw.set_lb(0);
                yaw.set_ub(0);
                //                x1.set_val(0);
                //                x2.set_val(0);
            }
            else if (axis == "y") {
                Lidar.min(sum(x_diff)/cells.size() + sum(z_diff)/cells.size());
                yaw.set_lb(0);
                yaw.set_ub(0);
                pitch.set_lb(0);
                pitch.set_ub(0);
                //                y1.set_val(0);
                //                y2.set_val(0);
            }
            else if (axis == "z") {
                Lidar.min(sum(x_diff)/cells.size() + sum(y_diff)/cells.size());
                pitch.set_lb(0);
                pitch.set_ub(0);
                roll.set_lb(0);
                roll.set_ub(0);
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
            //Lidar.print();
            solver<> S(Lidar,ipopt);
            S.run();
           // Lidar.print_solution();
            
            
            
            DebugOn("Pitch = " << pitch.eval()*180/pi << endl);
            DebugOn("Roll = " << roll.eval()*180/pi << endl);
            DebugOn("Yaw = " << yaw.eval()*180/pi << endl);
            //            DebugOn("Pitch2 = " << pitch2.eval()*180/pi << endl);
            //            DebugOn("Roll2 = " << roll2.eval()*180/pi << endl);
            //            DebugOn("Yaw2 = " << yaw2.eval()*180/pi << endl);
            roll_1 = roll.eval()*180/pi;
            pitch_1 = pitch.eval()*180/pi;
            yaw_1 = yaw.eval()*180/pi;
            
            if(abs(theta11.eval()-cos(roll.eval())*cos(yaw.eval()))>=1e-6){
                  DebugOn("T11.eval()"<<endl);
            }
            
            
           if(abs(theta12.eval()-(-1)*cos(roll.eval())*sin(yaw.eval()))>=1e-6){
             DebugOn("T12.eval()"<<endl);
           }
            
            if(abs(theta13.eval()-sin(roll.eval()))>=1e-6){
             DebugOn("T13.eval()"<<endl);
            }
           if(abs(theta21.eval()-cos(pitch.eval())*sin(yaw.eval())-cos(yaw.eval())*sin(roll.eval())*sin(pitch.eval()))>=1e-6){
           DebugOn("T21.eval()"<<endl);
           }
            
            
           if(abs(theta22.eval()-cos(pitch.eval())*cos(yaw.eval())-(-1)*sin(roll.eval())*sin(pitch.eval())*sin(yaw.eval()))>=1e-6){
           DebugOn("T22.eval()"<<endl);
           }
               
           if(abs(theta23.eval()-(-1)*cos(roll.eval())*sin(pitch.eval()))>=1e-6){
            DebugOn("T23.eval()"<<endl);
           }
          
            
            
            if(abs(theta31.eval()-sin(pitch.eval())*sin(yaw.eval())-(-1)*cos(pitch.eval())*cos(yaw.eval())*sin(roll.eval()))>=1e-6){
                  DebugOn("T31.eval()"<<endl);
            }
            
            
            if(abs(theta32.eval()-cos(yaw.eval())*sin(pitch.eval())-cos(pitch.eval())*sin(roll.eval())*sin(yaw.eval()))>=1e-6){
             DebugOn("T32.eval()"<<endl);
            }
            
            if(abs(theta33.eval()-cos(roll.eval())*cos(pitch.eval()))>=1e-6){
                  DebugOn("T33.eval()"<<endl);
            }
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
            Lidar.add(yaw1.in(range(0,0)),pitch1.in(range(0,0)),roll1.in(range(0,0)));
            Lidar.add(yaw2.in(range(0,0)),pitch2.in(range(0,0)),roll2.in(range(0,0)));
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
            //            DebugOn("Pitch2 = " << pitch2.eval()*180/pi << endl);
            //            DebugOn("Roll2 = " << roll2.eval()*180/pi << endl);
            //            DebugOn("Yaw2 = " << yaw2.eval()*180/pi << endl);
            roll_1 = roll1.eval()*180/pi;
            pitch_1 = pitch1.eval()*180/pi;
            yaw_1 = yaw1.eval()*180/pi;
            //            return 0;
        }
    }
    return {roll_1, pitch_1, yaw_1};
}
tuple<double,double,double> run_IPH(vector<vector<double>>& ext_model, vector<vector<double>>& ext_data, const vector<vector<double>>& uav1, const vector<vector<double>>& uav2, const vector<vector<double>>& rollpitchyawins1, const vector<vector<double>>& rollpitchyawins2){
    double roll = 0, pitch = 0, yaw = 1; /* at least one nonzero to enter the while loop */
    double final_roll = 0, final_pitch = 0, final_yaw = 0;
    int nb_iter = 0, max_nb_iter = 20;
    tuple<double,double,double> res;
    while(nb_iter < max_nb_iter && std::abs(roll)+std::abs(pitch)+std::abs(yaw)>0.5) {
        res = run_NLP("full", ext_model, ext_data, uav1, uav2, rollpitchyawins1, rollpitchyawins2);
        roll = get<0>(res);pitch = get<1>(res);yaw = get<2>(res);
        final_roll += roll;final_pitch += pitch;final_yaw += yaw;
        apply_transform_new_order(roll*pi/180, pitch*pi/180, yaw*pi/180, ext_model, uav1, rollpitchyawins1, 0.0,0.0,0.0);
        apply_transform_new_order(roll*pi/180, pitch*pi/180, yaw*pi/180, ext_data, uav2, rollpitchyawins2, 0.0,0.0,0.0);
        nb_iter++;
        DebugOn("No projection, ITERATION " << nb_iter << endl);
    }
    nb_iter = 0;yaw=1;max_nb_iter = 100;
    while(nb_iter < max_nb_iter && std::abs(roll)+std::abs(pitch)+std::abs(yaw)>0.5) {
        res = run_NLP("z", ext_model, ext_data, uav1, uav2, rollpitchyawins1, rollpitchyawins2);
        roll = get<0>(res);pitch = get<1>(res);yaw = get<2>(res);
        final_roll += roll;final_pitch += pitch;final_yaw += yaw;
        apply_transform_new_order(roll*pi/180, pitch*pi/180, yaw*pi/180, ext_model, uav1, rollpitchyawins1, 0.0,0.0,0.0);
        apply_transform_new_order(roll*pi/180, pitch*pi/180, yaw*pi/180, ext_data, uav2, rollpitchyawins2, 0.0,0.0,0.0);
        nb_iter++;
        DebugOn("Projceting out z axis, ITERATION " << nb_iter << endl);
    }
    nb_iter = 0;yaw=1;max_nb_iter = 100;
    while(nb_iter < max_nb_iter && std::abs(roll)+std::abs(pitch)+std::abs(yaw)>0.5) {
        res = run_NLP("y", ext_model, ext_data, uav1, uav2,  rollpitchyawins1, rollpitchyawins2);
        roll = get<0>(res);pitch = get<1>(res);yaw = get<2>(res);
        final_roll += roll;final_pitch += pitch;final_yaw += yaw;
        apply_transform_new_order(roll*pi/180, pitch*pi/180, yaw*pi/180, ext_model, uav1, rollpitchyawins1, 0.0,0.0,0.0);
        apply_transform_new_order(roll*pi/180, pitch*pi/180, yaw*pi/180, ext_data, uav2, rollpitchyawins2, 0.0,0.0,0.0);
        nb_iter++;
        DebugOn("Projceting out y axis, ITERATION " << nb_iter << endl);
    }
    nb_iter = 0;yaw=1;max_nb_iter = 100;
    while(nb_iter < max_nb_iter && std::abs(roll)+std::abs(pitch)+std::abs(yaw)>0.5) {
        res = run_NLP("x", ext_model, ext_data, uav1, uav2,  rollpitchyawins1, rollpitchyawins2);
        roll = get<0>(res);pitch = get<1>(res);yaw = get<2>(res);
        final_roll += roll;final_pitch += pitch;final_yaw += yaw;
        apply_transform_new_order(roll*pi/180, pitch*pi/180, yaw*pi/180, ext_model, uav1, rollpitchyawins1, 0.0,0.0,0.0);
        apply_transform_new_order(roll*pi/180, pitch*pi/180, yaw*pi/180, ext_data, uav2, rollpitchyawins2, 0.0,0.0,0.0);
        nb_iter++;
        DebugOn("Projceting out x axis, ITERATION " << nb_iter << endl);
    }
    nb_iter = 0;yaw=1;max_nb_iter = 100;
    while(nb_iter < max_nb_iter && std::abs(roll)+std::abs(pitch)+std::abs(yaw)>0.5) {
        res = run_NLP("full", ext_model, ext_data, uav1, uav2,  rollpitchyawins1, rollpitchyawins2);
        roll = get<0>(res);pitch = get<1>(res);yaw = get<2>(res);
        final_roll += roll;final_pitch += pitch;final_yaw += yaw;
        apply_transform_new_order(roll*pi/180, pitch*pi/180, yaw*pi/180, ext_model, uav1, rollpitchyawins1, 0.0,0.0,0.0);
        apply_transform_new_order(roll*pi/180, pitch*pi/180, yaw*pi/180, ext_data, uav2, rollpitchyawins2, 0.0,0.0,0.0);
        nb_iter++;
        DebugOn("No projection, ITERATION " << nb_iter << endl);
    }
    DebugOn("Final Roll (degrees) = " << final_roll << endl);
    DebugOn("Final Pitch (degrees) = " << final_pitch << endl);
    DebugOn("Final Yaw (degrees) = " << final_yaw << endl);
    DebugOn("Final Roll (radians)= " << final_roll *pi/180 << endl);
    DebugOn("Final Pitch (radians) = " << final_pitch *pi/180 << endl);
    DebugOn("Final Yaw (radians) = " << final_yaw *pi/180 << endl);
//    apply_transform_new_order(roll*pi/180, pitch*pi/180, yaw*pi/180, ext_model, uav1, rollpitchyawins1, 0.0,0.0,0.0);
//    apply_transform_new_order(roll*pi/180, pitch*pi/180, yaw*pi/180, ext_data, uav2, rollpitchyawins2, 0.0,0.0,0.0);
    return {final_roll,final_pitch,final_yaw};
}
#endif /* IPH_h */
