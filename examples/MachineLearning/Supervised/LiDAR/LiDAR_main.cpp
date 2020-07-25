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
#include <time.h>
using namespace std;

int main (int argc, char * argv[])
{
    bool show_umbrella = false;
    if(show_umbrella){
        vector<double> x_vec0, y_vec0, z_vec0, x_vec, y_vec, z_vec;
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
        double x_rot1, y_rot1, z_rot1;
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
                }}}
        
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
//        return 0;
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
        vector<double> x_vec1,y_vec1,z_vec1,zmin_vec1,zmax_vec1;
        vector<double> x_vec2,y_vec2,z_vec2,zmin_vec2,zmax_vec2;
        vector<double> x_shift1,y_shift1,z_shift1;
        vector<double> x_shift2,y_shift2,z_shift2;
        vector<double> x_shift,y_shift,z_shift;
        vector<double> uav_roll1,uav_pitch1,uav_yaw1;
        vector<double> uav_roll2,uav_pitch2,uav_yaw2;
        vector<double> x_combined,y_combined,z_combined,zmin_combined,zmax_combined;
        set<double> timestamps;
        size_t uav_id = 0;
        bool new_uav = true, u_turn = false, frame1 = true;
        double unix_time, delta_x = 0, delta_y = 0;
        pair<map<int,shared_ptr<Frame>>::iterator,bool> frame_ptr;
        for (int i = 0; i< nb_rows-1; i++) { // Input iterator
            auto laser_id = CSV_data.GetCell<int>("LaserID", i);
            if(laser_id!=15){/* Only keep points from Nadir laser */
                continue;
            }
            auto frame_id = CSV_data.GetCell<int>(0, i);
            new_uav = (uav_id==0) || (UAVPoints[uav_id-1]->_frame_id != frame_id);
            if(new_uav){
                auto uav_x1 = CSV_data.GetCell<double>("Track_UTM_E", i);
                auto uav_y1 = CSV_data.GetCell<double>("Track_UTM_N", i);
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
                UAVPoints[uav_id]->_height = CSV_data.GetCell<double>("Track_UTM_Height", i);
                unix_time = CSV_data.GetCell<double>("Time", i);
                UAVPoints[uav_id]->set_unix_time(unix_time);
                uav_x.push_back(UAVPoints[uav_id]->_x);
                uav_y.push_back(UAVPoints[uav_id]->_y);
                uav_z.push_back(UAVPoints[uav_id]->_height);
                frame_ptr = frames.insert(make_pair(UAVPoints[uav_id]->_frame_id, make_shared<Frame>(UAVPoints[uav_id]->_frame_id, UAVPoints[uav_id]->_unix_time)));
                frame_ptr.first->second->add_UAV_point(UAVPoints[uav_id]);
                if(frame1){/* Has not performed a u-turn yet, keep adding to frames1 */
                    frames1.insert(make_pair(frame_ptr.first->second->_id, frame_ptr.first->second));
                }
                else{/* Already turned, keep adding to frames2 */
                    frames2.insert(make_pair(frame_ptr.first->second->_id, frame_ptr.first->second));
                }
                if(u_turn){
                    DebugOn("Detected a Uturn at frame " << frame_ptr.first->first << endl);
                }
                uav_id++;
            }
            
            auto xpos = CSV_data.GetCell<double>("UTM_E", i);
            auto ypos = CSV_data.GetCell<double>("UTM_N", i);
            auto zpos = CSV_data.GetCell<double>("UTM_Height", i);
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
        param<> x1("x1"), x2("x2"), y1("y1"), y2("y2"), z1("z1"), z2("z2");
        param<> x_uav1("x_uav1"), y_uav1("y_uav1"), z_uav1("z_uav1");
        param<> x_uav2("x_uav2"), y_uav2("y_uav2"), z_uav2("z_uav2");
        int nb_pts_per_frame1 = 0, nb_pts_per_frame2 = 0, idx1 = 1, idx2 = 1;
        for (const auto &frame: frames1) {
            nb_pts_per_frame1 += frame.second->_lidar_points->size();
            int i = 0;
            for (const auto &p: *frame.second->_lidar_points) {
                if(i%10==0){
                    x_vec1.push_back(p->_x);
                    x_shift1.push_back(frame.second->_uav_point->_x);
                    x_uav1.add_val(to_string(idx1),x_shift1.back());
                    x1.add_val(to_string(idx1),x_vec1.back());
                    y_vec1.push_back(p->_y);
                    y_shift1.push_back(frame.second->_uav_point->_y);
                    y_uav1.add_val(to_string(idx1),y_shift1.back());
                    y1.add_val(to_string(idx1),y_vec1.back());
                    z_vec1.push_back(p->_z);
                    z_shift1.push_back(frame.second->_uav_point->_height);
                    z_uav1.add_val(to_string(idx1),z_shift1.back());
                    z1.add_val(to_string(idx1),z_vec1.back());
                    idx1++;
                }
                i++;
            }
            
        }
        for (const auto &frame: frames2) {
            nb_pts_per_frame2 += frame.second->_lidar_points->size();
            int i = 0;
            for (auto const &p: *frame.second->_lidar_points) {
                if(i%100==0){
                    x_vec2.push_back(p->_x);
                    x_shift2.push_back(frame.second->_uav_point->_x);
                    x_uav2.add_val(to_string(idx2),x_shift2.back());
                    x2.add_val(to_string(idx2),x_vec2.back());
                    y_vec2.push_back(p->_y);
                    y_shift2.push_back(frame.second->_uav_point->_y);
                    y_uav2.add_val(to_string(idx2),y_shift2.back());
                    y2.add_val(to_string(idx2),y_vec2.back());
                    z_vec2.push_back(p->_z);
                    z_shift2.push_back(frame.second->_uav_point->_height);
                    z_uav2.add_val(to_string(idx2),z_shift2.back());
                    z2.add_val(to_string(idx2),z_vec2.back());
                    idx2++;
                }
                i++;
            }
        }
        if(frames1.size()!=0)
            DebugOn("Average number of points per frame in flight line 1 = " << nb_pts_per_frame1/frames1.size() << endl);
        if(frames2.size()!=0)
            DebugOn("Average number of points per frame in flight line 2 = " << nb_pts_per_frame2/frames2.size() << endl);
        
        DebugOn("Number of selected points in flight line 1 = " << x_vec1.size() << endl);
        DebugOn("Number of selected points in flight line 2 = " << x_vec2.size() << endl);
        int n = std::min(x_vec1.size(),x_vec2.size());
        indices N("N");
        N = range(1, n);
        bool solve_lidar = true;
        if (solve_lidar) {
            Model<> M("Lidar");
            var<> new_x1("new_x1"), new_y1("new_y1"), new_z1("new_z1");
            var<> new_x2("new_x2"), new_y2("new_y2"), new_z2("new_z2");
            var<> x_diff("x_diff", pos_), y_diff("y_diff", pos_), z_diff("z_diff", pos_);
            var<> yaw1("yaw1", -0.1, 0.1), pitch1("pitch1", -0.1, 0.1), roll1("roll1", -0.1, 0.1);
            var<> yaw2("yaw2", -0.1, 0.1), pitch2("pitch2", -0.1, 0.1), roll2("roll2", -0.1, 0.1);
            
            M.add(yaw1.in(R(1)),pitch1.in(R(1)),roll1.in(R(1)));
            M.add(yaw2.in(R(1)),pitch2.in(R(1)),roll2.in(R(1)));
            M.add(new_x1.in(N), new_y1.in(N), new_z1.in(N));
            M.add(new_x2.in(N), new_y2.in(N), new_z2.in(N));
            M.add(x_diff.in(N), y_diff.in(N), z_diff.in(N));

            Constraint<> Equal_pitch("Equal_pitch");
            Equal_pitch += pitch1 - pitch2;
            M.add(Equal_pitch==0);

            Constraint<> Opp_roll("Opp_roll");
            Opp_roll += roll1 + roll2;
            M.add(Opp_roll==0);

            Constraint<> Opp_yaw("Opp_yaw");
            Opp_yaw += yaw1 + yaw2;
            M.add(Opp_yaw==0);


            Constraint<> z_abs1("z_abs1");
            z_abs1 += z_diff - (new_z1 - new_z2);
            M.add(z_abs1.in(N)>=0);

            Constraint<> z_abs2("z_abs2");
            z_abs2 += z_diff - (new_z2 - new_z1);
            M.add(z_abs2.in(N)>=0);
            
            Constraint<> x_abs1("x_abs1");
            x_abs1 += x_diff - (new_x1 - new_x2);
            M.add(x_abs1.in(N)>=0);
    
            Constraint<> x_abs2("x_abs2");
            x_abs2 += x_diff - (new_x2 - new_x1);
            M.add(x_abs2.in(N)>=0);
    
    
            Constraint<> y_abs1("y_abs1");
            y_abs1 += y_diff - (new_y1 - new_y2);
            M.add(y_abs1.in(N)>=0);
    
            Constraint<> y_abs2("y_abs2");
            y_abs2 += y_diff - (new_y2 - new_y1);
            M.add(y_abs2.in(N)>=0);

            auto ids = yaw1.repeat_id(n);

            /* alpha = yaw_, beta = pitch_ and gamma = roll_ */
            Constraint<> x_rot1("x_rot1");
            x_rot1 += new_x1 - x_uav1.in(N);
            x_rot1 -= (x1.in(N)-x_uav1.in(N))*cos(yaw1.in(ids))*cos(roll1.in(ids)) + (y1.in(N)-y_uav1.in(N))*(cos(yaw1.in(ids))*sin(roll1.in(ids))*sin(pitch1.in(ids)) - sin(yaw1.in(ids))*cos(pitch1.in(ids))) + (z1.in(N)-z_uav1.in(N))*(cos(yaw1.in(ids))*sin(roll1.in(ids))*cos(pitch1.in(ids)) + sin(yaw1.in(ids))*sin(pitch1.in(ids)));
            M.add(x_rot1.in(N)==0);

            Constraint<> x_rot2("x_rot2");
            x_rot2 += new_x2 - x_uav2.in(N);
            x_rot2 -= (x2.in(N)-x_uav2.in(N))*cos(yaw2.in(ids))*cos(roll2.in(ids)) + (y2.in(N)-y_uav2.in(N))*(cos(yaw2.in(ids))*sin(roll2.in(ids))*sin(pitch2.in(ids)) - sin(yaw2.in(ids))*cos(pitch2.in(ids))) + (z2.in(N)-z_uav2.in(N))*(cos(yaw2.in(ids))*sin(roll2.in(ids))*cos(pitch2.in(ids)) + sin(yaw2.in(ids))*sin(pitch2.in(ids)));
            M.add(x_rot2.in(N)==0);


            Constraint<> y_rot1("y_rot1");
            y_rot1 += new_y1 - y_uav1.in(N);
            y_rot1 -= (x1.in(N)-x_uav1.in(N))*sin(yaw1.in(ids))*cos(roll1.in(ids)) + (y1.in(N)-y_uav1.in(N))*(sin(yaw1.in(ids))*sin(roll1.in(ids))*sin(pitch1.in(ids)) + cos(yaw1.in(ids))*cos(pitch1.in(ids))) + (z1.in(N)-z_uav1.in(N))*(sin(yaw1.in(ids))*sin(roll1.in(ids))*cos(pitch1.in(ids)) - cos(yaw1.in(ids))*sin(pitch1.in(ids)));
            M.add(y_rot1.in(N)==0);

            Constraint<> y_rot2("y_rot2");
            y_rot2 += new_y2 - y_uav2.in(N);
            y_rot2 -= (x2.in(N)-x_uav2.in(N))*sin(yaw2.in(ids))*cos(roll2.in(ids)) + (y2.in(N)-y_uav2.in(N))*(sin(yaw2.in(ids))*sin(roll2.in(ids))*sin(pitch2.in(ids)) + cos(yaw2.in(ids))*cos(pitch2.in(ids))) + (z2.in(N)-z_uav2.in(N))*(sin(yaw2.in(ids))*sin(roll2.in(ids))*cos(pitch2.in(ids)) - cos(yaw2.in(ids))*sin(pitch2.in(ids)));
            M.add(y_rot2.in(N)==0);

            Constraint<> z_rot1("z_rot1");
            z_rot1 += new_z1 - z_uav1.in(N);
            z_rot1 -= (x1.in(N)-x_uav1.in(N))*sin(-1*roll1.in(ids)) + (y1.in(N)-y_uav1.in(N))*(cos(roll1.in(ids))*sin(pitch1.in(ids))) + (z1.in(N)-z_uav1.in(N))*(cos(roll1.in(ids))*cos(pitch1.in(ids)));
            M.add(z_rot1.in(N)==0);


            Constraint<> z_rot2("z_rot2");
            z_rot2 += new_z2 - z_uav2.in(N);
            z_rot2 -= (x2.in(N)-x_uav2.in(N))*sin(-1*roll2.in(ids)) + (y2.in(N)-y_uav2.in(N))*(cos(roll2.in(ids))*sin(pitch2.in(ids))) + (z2.in(N)-z_uav2.in(N))*(cos(roll2.in(ids))*cos(pitch2.in(ids)));
            M.add(z_rot2.in(N)==0);

            M.min(sum(x_diff) + sum(y_diff) + sum(z_diff));

            M.print();

            solver<> S(M,ipopt);
            S.run();
        }
        
        /* Adjusting the data */
        double roll = 1.375, pitch = 0.9, yaw = -0.5; /* In degrees */
        double shifted_x, shifted_y, shifted_z;
        double beta = roll*pi/180;// roll in radians
        double gamma = pitch*pi/180; // pitch in radians
        double alpha = yaw*pi/180; // yaw in radians
        
        auto tot_pts = x_vec1.size()+x_vec2.size();
        x_combined.resize(tot_pts);
        y_combined.resize(tot_pts);
        z_combined.resize(tot_pts);


        for (auto i = 0; i< x_vec1.size(); i++) {
            shifted_x = x_vec1[i] - x_shift1[i];
            shifted_y = y_vec1[i] - y_shift1[i];
            shifted_z = z_vec1[i] - z_shift1[i];
            x_combined[i] = shifted_x*cos(alpha)*cos(beta) + shifted_y*(cos(alpha)*sin(beta)*sin(gamma) - sin(alpha)*cos(gamma)) + shifted_z*(cos(alpha)*sin(beta)*cos(gamma) + sin(alpha)*sin(gamma));
            y_combined[i] = shifted_x*sin(alpha)*cos(beta) + shifted_y*(sin(alpha)*sin(beta)*sin(gamma) + cos(alpha)*cos(gamma)) + shifted_z*(sin(alpha)*sin(beta)*cos(gamma) - cos(alpha)*sin(gamma));
            z_combined[i] = shifted_x*(-sin(beta)) + shifted_y*(cos(beta)*sin(gamma)) + shifted_z*(cos(beta)*cos(gamma));
            x_combined[i] += x_shift1[i];
            y_combined[i] += y_shift1[i];
            z_combined[i] += z_shift1[i];
        }
        beta *= -1;
        alpha *= -1;
        for (auto i = 0; i< x_vec2.size(); i++) {
            shifted_x = x_vec2[i] - x_shift2[i];
            shifted_y = y_vec2[i] - y_shift2[i];
            shifted_z = z_vec2[i] - z_shift2[i];
            x_combined[x_vec1.size()+i] = shifted_x*cos(alpha)*cos(beta) + shifted_y*(cos(alpha)*sin(beta)*sin(gamma) - sin(alpha)*cos(gamma)) + shifted_z*(cos(alpha)*sin(beta)*cos(gamma) + sin(alpha)*sin(gamma));
            y_combined[x_vec1.size()+i] = shifted_x*sin(alpha)*cos(beta) + shifted_y*(sin(alpha)*sin(beta)*sin(gamma) + cos(alpha)*cos(gamma)) + shifted_z*(sin(alpha)*sin(beta)*cos(gamma) - cos(alpha)*sin(gamma));
            z_combined[x_vec1.size()+i] = shifted_x*(-sin(beta)) + shifted_y*(cos(beta)*sin(gamma)) + shifted_z*(cos(beta)*cos(gamma));
            x_combined[x_vec1.size()+i] += x_shift2[i];
            y_combined[x_vec1.size()+i] += y_shift2[i];
            z_combined[x_vec1.size()+i] += z_shift2[i];
        }
        bool plot_data = true;
        if(plot_data){
            namespace plt = matplotlibcpp;
            
            std::map<std::string, std::string> keywords, keywords2;
            keywords["marker"] = "s";
            keywords["linestyle"] = "None";
            keywords["ms"] = "0.05";
            //    plt::plot3(x_combined, y_combined, zmax_combined, keywords);
            plt::plot3(x_vec1, y_vec1, z_vec1, x_vec2, y_vec2, z_vec2, keywords);
            plt::plot3(x_combined, y_combined, z_combined, keywords);
            //        plt::plot3(uav_x, uav_y, uav_z, x_combined, y_combined, zmax_combined, keywords);
            keywords2["marker"] = "s";
            keywords2["ms"] = "0.1";
            //    plt::plot3(uav_x1, uav_y1, uav_z1, uav_x, uav_y, uav_z, keywords2);
            //    plt::plot3(uav_x1, uav_y1, uav_z1, keywords);
            //    plt::plot3(x_vec2, y_vec2, zmax_vec2, keywords);
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
    }
//    else {
//
//        string log_level="0";
//        string LiDAR_file1 = string(prj_dir)+"/data_sets/LiDAR/test1.las";
//        string LiDAR_file2 = string(prj_dir)+"/data_sets/LiDAR/test2.las";
//        string GPS_file = string(prj_dir)+"/data_sets/LiDAR/RawGNSS.csv";
//
//    #ifdef USE_OPT_PARSER
//        /** create a OptionParser with options */
//        op::OptionParser opt;
//        opt.add_option("h", "help", "shows option help"); // no default value means boolean options, which default value is false
//        opt.add_option("f", "file", "Input file name (def. ../data_sets/LiDAR/*)", fname_phi );
//        opt.add_option("l", "log", "Log level (def. 0)", log_level );
//
//        /** parse the options and verify that all went well. If not, errors and help will be shown */
//        bool correct_parsing = opt.parse_options(argc, argv);
//
//        if(!correct_parsing){
//            return EXIT_FAILURE;
//        }
//
//        LiDAR_file = opt["f"];
//        output = op::str2int(opt["l"]);
//        output = 5;
//        bool has_help = op::str2bool(opt["h"]);
//        /** show help */
//        if(has_help) {
//            opt.show_help();
//            exit(0);
//        }
//    #else
//        if(argc>1){
//            LiDAR_file1 = argv[1];
//        }
//        if(argc>2){
//            LiDAR_file2 = argv[2];
//        }
//        if(argc>3){
//            GPS_file = argv[3];
//        }
//        if(argc>4){
//            char *p;
//            grid_step =  strtol(argv[4], &p, 10);
//        }
//        if(argc>5){
//            char *p;
//            optimize =  (strtol(argv[5], &p, 10)!=0);
//        }
//    #endif
//        rapidcsv::Document  GPS_data(GPS_file);
//        int n = GPS_data.GetRowCount();
//        vector<UAVPoint> UAVPoints(n);
//        vector<LidarPoint> LidarPoints;
//        map<double,Frame> frames, frames1;
//        vector<double> uav_x, uav_y, uav_z;
//        vector<double> uav_x1, uav_y1, uav_z1;
//        for (int i = 0; i< n-1; i++) { // Input iterator
//            UAVPoints[i]._x = (GPS_data.GetCell<double>("UtmPos_X", i));
//            UAVPoints[i]._y = (GPS_data.GetCell<double>("UtmPos_Y", i));
//            UAVPoints[i]._latitude = GPS_data.GetCell<double>("AbsPos_Y", i);
//            UAVPoints[i]._longitude = GPS_data.GetCell<double>("AbsPos_X", i);
//            UAVPoints[i]._height = GPS_data.GetCell<double>("AbsPos_Z", i);
//            UAVPoints[i]._roll = GPS_data.GetCell<double>("Roll", i);
//            UAVPoints[i]._pitch = GPS_data.GetCell<double>("Pitch", i);
//            UAVPoints[i]._yaw = GPS_data.GetCell<double>("Heading", i);
//            double unix_time = round(GPS_data.GetCell<double>("Time", i)*10);
//            UAVPoints[i].set_unix_time(unix_time);
//    //        uav_x.push_back(UAVPoints[i]._longitude+582690.8242);
//    //        uav_y.push_back(UAVPoints[i]._latitude+4107963.58);
//    //        uav_z.push_back(UAVPoints[i]._height);
//    //        DebugOff("Unix time = " << GPS_data.GetCell<double>("Unix Time", i) << endl);
//    //        DebugOff("Latitude = " << GPS_data.GetCell<double>("Latitude (degrees)", i) << endl);
//    //        DebugOff("Longitude = " << GPS_data.GetCell<double>("Longitude (degrees)", i) << endl);
//    //        DebugOff("Height = " << GPS_data.GetCell<double>("Height", i) << endl);
//    //        auto frame_ptr = frames.insert(make_pair(UAVPoints[i]._unix_time, Frame(UAVPoints[i]._unix_time)));
//    //        frame_ptr.first->second.add_UAV_point(UAVPoints[i]);
//    //        UAVPoints[i]._latitude = GPS_data.GetCell<double>("Latitude (degrees)", i);
//    //        UAVPoints[i]._longitude = GPS_data.GetCell<double>("Longitude (degrees)", i);
//    //        UAVPoints[i]._height = GPS_data.GetCell<double>("Height", i);
//    //        UAVPoints[i].set_unix_time(GPS_data.GetCell<double>("Unix Time", i));
//            /* If longitude = -116.06925521318, x = 582690.824176164 */
//    //        uav_x.push_back((UAVPoints[i]._longitude*582690.824176164/(-116.06925521318))*1e-5);
//            uav_x.push_back(UAVPoints[i]._x);
//            uav_y.push_back(UAVPoints[i]._y);
//            uav_z.push_back(UAVPoints[i]._height);
//            DebugOff("Unix time = " << GPS_data.GetCell<double>("Unix Time", i) << endl);
//            auto frame_ptr = frames.insert(make_pair(UAVPoints[i]._unix_time, Frame(UAVPoints[i]._unix_time)));
//            frame_ptr.first->second.add_UAV_point(UAVPoints[i]);
//        }
//        DebugOn("Read " << n << "UAV points" << endl);
//    //    LiDAR_file = "/Users/l297598/Downloads/PostPost_DAG4_lidar_LDRD_Reserve_20200326_LAS_and_ANPP_files/Roll_Pitch_Yaw_automation_2020_03_26/3.ins_postpost_SGZ_DAG4_L_2_2019_07_31_18_combined_LDRD_Reserve_20200326_Frames_2130_2222_RPY_1.375_0.9_-0.55_NoKin_adc.las";
//    //   LiDAR_file = "/Users/l297598/Downloads/PostPost_DAG4_lidar_LDRD_Reserve_20200326_LAS_and_ANPP_files/Roll_Pitch_Yaw_automation_2020_03_26/1.raw_postpost_SGZ_DAG4_L_2_2019_07_31_18_combined_LDRD_Reserve_20200326_Frames_2130_2222_RPY_0_0_0_NoKin_adc.las";
//    //    LiDAR_file = "/Users/l297598/Downloads/PostPost_DAG4_lidar_LDRD_Reserve_20200326_LAS_and_ANPP_files/Roll_Pitch_Yaw_automation_2020_03_26/6.ins_postpost_SGZ_DAG4_L_2_2019_07_31_18_combined_LDRD_Reserve_20200326_Frames_2130_2222_and_2511_2587_RPY_1.375_0.9_-0.55_NoKin_adc.las";
//    //   LiDAR_file = "/Users/l297598/Downloads/PostPost_DAG4_lidar_LDRD_Reserve_20200326_LAS_and_ANPP_files/Roll_Pitch_Yaw_automation_2020_03_26/5.raw_postpost_SGZ_DAG4_L_2_2019_07_31_18_combined_LDRD_Reserve_20200326_Frames_2130_2222_and_2511_2587_RPY_0_0_0_NoKin_adc.las";
//    //    LiDAR_file = "/Users/l297598/Downloads/PostPost_DAG4_lidar_LDRD_Reserve_20200326_LAS_and_ANPP_files/Roll_Pitch_Yaw_automation_2020_03_26/7.raw_postpost_SGZ_DAG4_L_2_2019_07_31_18_combined_LDRD_Reserve_20200326_Frames_2130_2587_RPY_0_0_0_NoKin_adc.las";
//    //    LiDAR_file = "/Users/l297598/Downloads/PostPost_DAG4_lidar_LDRD_Reserve_20200326_LAS_and_ANPP_files/Roll_Pitch_Yaw_automation_2020_03_26/8.ins_postpost_SGZ_DAG4_L_2_2019_07_31_18_combined_LDRD_Reserve_20200326_Frames_2130_2587_RPY_1.375_0.9_-0.55_NoKin_adc.las";
//    //    LiDAR_file = "/Users/l297598/Downloads/PostPost_DAG4_lidar_LDRD_Reserve_20200326_LAS_and_ANPP_files/PostPost_DAG4_lidar_LDRD_Reserve_20200326_LAS_files/raw_postpost_1g11d_DAG4_L_2_2019_07_31_18_combined_LDRD_Reserve_20200326_Frames_2892_3073_RPY_0_0_0_NoKin_adc.las";
//    //    LiDAR_file = "/Users/l297598/Downloads/PostPost_DAG4_lidar_LDRD_Reserve_20200326_LAS_and_ANPP_files/PostPost_DAG4_lidar_LDRD_Reserve_20200326_LAS_files/ins_postpost_1g11d_DAG4_L_2_2019_07_31_18_combined_LDRD_Reserve_20200326_Frames_2892_3073_RPY_1.375_0.9_-0.55_NoKin_adc.las";
//        LASreadOpener lasreadopener;
//        lasreadopener.set_file_name(LiDAR_file1.c_str());
//        lasreadopener.set_populate_header(TRUE);
//        param<> x1("x1"), x2("x2"), y1("y1"), y2("y2"), z1("z1"), z2("z2");
//        param<> x_uav1("x_uav1"), y_uav1("y_uav1"), z_uav1("z_uav1");
//        param<> x_uav2("x_uav2"), y_uav2("y_uav2"), z_uav2("z_uav2");
//        param<> x_uav("x_uav"), y_uav("y_uav"), z_uav("z_uav");
//        param<> pre_x("pre_x"), pre_y("pre_y"), pre_z("pre_z");
//        int xdim1=0, ydim1=0, zdim1=0;
//        if (!lasreadopener.active())
//        {
//            throw invalid_argument("ERROR: no input specified\n");
//        }
//        vector<double> x_vec1,y_vec1,z_vec1,zmin_vec1,zmax_vec1;
//        vector<double> x_vec2,y_vec2,z_vec2,zmin_vec2,zmax_vec2;
//        vector<double> x_shift1,y_shift1,z_shift1;
//        vector<double> x_shift2,y_shift2,z_shift2;
//        vector<double> x_shift,y_shift,z_shift;
//        vector<double> uav_roll1,uav_pitch1,uav_yaw1;
//        vector<double> uav_roll2,uav_pitch2,uav_yaw2;
//        vector<double> x_combined,y_combined,z_combined,zmin_combined,zmax_combined;
//        set<double> timestamps;
//        while (lasreadopener.active())
//        {
//            LASreader* lasreader = lasreadopener.open();
//            if (lasreader == 0)
//            {
//                throw invalid_argument("ERROR: could not open lasreader\n");
//            }
//
//            xdim1 = ceil(lasreader->header.max_x*grid_step - lasreader->header.min_x*grid_step);
//            ydim1 = ceil(lasreader->header.max_y*grid_step - lasreader->header.min_y*grid_step);
//    //        zdim1 = ceil(lasreader->header.max_z*100*grid_step - lasreader->header.min_z*100*grid_step);
//            ground_z = floor(lasreader->header.min_z);
//            size_t nb_pts_per_frame = lasreader->npoints/frames.size()+1;
//            DebugOn("Number of points = " << lasreader->npoints << endl);
//            DebugOn("Number of points per frame = " <<  nb_pts_per_frame << endl);
//            DebugOn("min x axis = " << lasreader->header.min_x << endl);
//            DebugOn("max x axis = " << lasreader->header.max_x << endl);
//            DebugOn("min y axis = " << lasreader->header.min_y << endl);
//            DebugOn("max y axis = " << lasreader->header.max_y << endl);
//    //        uav_x1.push_back(lasreader->header.min_x);
//    //        uav_y1.push_back(lasreader->header.min_y);
//    //        uav_z1.push_back(ground_z);
//            DebugOn("dimension in x axis = " << xdim1 << endl);
//            DebugOn("dimension in y axis = " << ydim1 << endl);
//    //        DebugOn("dimension in z axis = " << zdim1 << endl);
//
//
//            int nb_dots; /* Number of measurements inside cell */
//            double xpos, ypos;
//            double z, min_z, max_z, av_z;
//            pair<double,double> pos;
//            size_t nb_pts = 0;
//            tuple<int,double,double,double,UAVPoint*> cell; /* <nb_dots,min_z,max_z,av_z> */
//            /* Now get rid of the first points
//            lasreader->read_point();
//            auto gps_time = lasreader->point.get_gps_time();
//            while (lasreader->point.get_gps_time()==gps_time)
//            {
//                lasreader->read_point();
//            }
//             */
//            auto frames_it = frames.begin();
//            while (lasreader->read_point())
//            {
//                if(nb_pts==0){
//                    DebugOn(to_string_with_precision(10.*(lasreader->point.get_gps_time()+315964800. - 18.),24) << ": (" << to_string_with_precision(lasreader->point.get_x(),10) <<"," << to_string_with_precision(lasreader->point.get_y(),10) << ","<< to_string_with_precision(lasreader->point.get_z(),10) <<")"<<endl);
//    //                return 0;
//                }
//                xpos = (((lasreader->point.get_x())));
//                ypos = (((lasreader->point.get_y())));
//                if (ypos < 0)
//                {
//                    fprintf(stderr, "ERROR: ypos = %d\n", ypos);
//                }
//
//
//                //            z = lasreader->point.get_Z();
//                //            z = floor((lasreader->point.get_z()*100*grid_step - lasreader->header.min_z*100*grid_step));
//                z = lasreader->point.get_z();
//                timestamps.insert(lasreader->point.get_gps_time()+1./(nb_pts%10));
//                DebugOff("GPS time = " << to_string_with_precision(lasreader->point.get_gps_time(),24) << endl);
//                LidarPoints.push_back(LidarPoint(lasreader->point.get_gps_time(),xpos,ypos,z));
//                if(nb_pts>=nb_pts_per_frame && nb_pts%nb_pts_per_frame==0 && next(frames_it)!=frames.end()){
//                    frames_it++;
//                }
//                auto frame_ptr = frames.insert(make_pair(frames_it->first, Frame(LidarPoints.back()._unix_time*10+(nb_pts%10))));
//                frame_ptr.first->second.add_lidar_point(LidarPoints.back());
//                if(frame_ptr.second){
//                    throw invalid_argument("Frame with missing UAV point at " + to_string(LidarPoints.back()._unix_time*10+(nb_pts%10)));
//                }
//                else{
//                    LidarPoints.back()._uav_pt = frame_ptr.first->second._uav_point;
//                    frames1.insert(make_pair(LidarPoints.back()._unix_time*10+(nb_pts%10), Frame(frame_ptr.first->second)));
//    //                uav_x1.push_back((frame_ptr.first->second._uav_points.front()->_longitude+582690.8242)*1e-5);
//    //                uav_y1.push_back((frame_ptr.first->second._uav_points.front()->_latitude+4107963.58)*1e-5);
//    //                uav_z1.push_back(frame_ptr.first->second._uav_points.front()->_height*100);
//                }
//                DebugOff("Added lidar point to frame at " << LidarPoints.back()._unix_time << " : " << LidarPoints.back()._hour << ":" << LidarPoints.back()._minutes << ":" << LidarPoints.back()._seconds << endl);
//                //            xpos = floor((lasreader->point.get_x()*grid_step - lasreader->header.min_x*grid_step));
//                //            ypos = floor((lasreader->point.get_y()*grid_step - lasreader->header.min_y*grid_step));
//
//                //            z = round((lasreader->point.get_z()*100*grid_step));
//    //            xpos -= frame_ptr.first->second._uav_point->_x;
//    //            ypos -= frame_ptr.first->second._uav_point->_y;
//    //            z -= frame_ptr.first->second._uav_point->_height;
//                pos = make_pair(xpos,ypos);
//                cell = make_tuple(1,z,z,z,frame_ptr.first->second._uav_point);
//                DebugOff("new z value = " << z << endl);
//                DebugOff("fractional z value = " << lasreader->point.get_z() << endl);
//                DebugOff("z range = " << grid_maxz[pos] - grid_minz[pos] << endl << endl);
//                auto map_pair = grid1.insert(make_pair(pos,cell));
//                auto map_pair_all = grid.insert(make_pair(pos,cell));
//                if (!map_pair.second) {/* Cell already exists */
//    //                DebugOn("Found overlapping cell after " << nb_pts << endl);
//    //                DebugOn("(" << to_string_with_precision(lasreader->point.get_x(),10) <<"," << to_string_with_precision(lasreader->point.get_y(),10) << ","<< to_string_with_precision(lasreader->point.get_z(),10) <<")"<<endl);
//    //                                return 0;
//                    nb_dots = get<0>(map_pair.first->second);
//                    min_z = get<1>(map_pair.first->second);
//                    max_z = get<2>(map_pair.first->second);
//                    av_z = get<3>(map_pair.first->second);
//                    /* Update cell data */
//                    get<0>(map_pair.first->second)++;
//                    if(min_z>z){
//                        DebugOff("updating  min z at cell [" << xpos << "," << ypos << "]" << endl);
//                        DebugOff("previous min z = " << min_z << endl);
//                        DebugOff("new min z = " << z << endl);
//                        get<1>(map_pair.first->second) = z;
//                    }
//                    if(max_z<z){
//                        DebugOff("updating  max z at cell [" << xpos << "," << ypos << "]" << endl);
//                        DebugOff("previous max z = " << max_z << endl);
//                        DebugOff("new max z = " << z << endl);
//                        get<2>(map_pair.first->second) = z;
//                    }
//                    av_z *= nb_dots;/* remove previous denominator */
//                    get<3>(map_pair.first->second) = (av_z + z)/(nb_dots+1); /* Update denominator */
//                }
//                if (!map_pair_all.second) {/* Cell already exists */
//                    nb_dots = get<0>(map_pair_all.first->second);
//                    min_z = get<1>(map_pair_all.first->second);
//                    max_z = get<2>(map_pair_all.first->second);
//                    av_z = get<3>(map_pair_all.first->second);
//                    /* Update cell data */
//                    get<0>(map_pair_all.first->second)++;
//                    if(min_z>z){
//                        get<1>(map_pair_all.first->second) = z;
//                    }
//                    if(max_z<z){
//                        get<2>(map_pair_all.first->second) = z;
//                    }
//                    av_z *= nb_dots;/* remove previous denominator */
//                    get<3>(map_pair_all.first->second) = (av_z + z)/(nb_dots+1); /* Update denominator */
//                }
//                nb_pts++;
//            }
//            DebugOn("Read " << nb_pts << " points" << endl);
//            DebugOn("2D grid has " << grid1.size() << " cells " << endl);
//            DebugOn("xdim = " << xdim1 << endl);
//            DebugOn("ydim = " << ydim1 << endl);
//
//            int cell_counter = 0, max_nb_dots = 0;
//            double av_nb_dots = 0, max_z_range = 0, z_range = 0,  av_z_range = 0;
//            for (auto &cell : grid1) {
//                DebugOff("Cell num " << cell_counter++ << " at [" << cell.first.first << "," << cell.first.second << "]" << endl);
//                nb_dots = get<0>(cell.second);
//                min_z = get<1>(cell.second);
//                max_z = get<2>(cell.second);
//                av_z = get<3>(cell.second);
//                x_vec1.push_back(cell.first.first);
//                x_shift1.push_back(get<4>(cell.second)->_x);
//                x_shift.push_back(get<4>(cell.second)->_x);
//                x_uav1.add_val(to_string(cell.first.first)+","+to_string(cell.first.second),(get<4>(cell.second)->_x));
//                x1.add_val(to_string(cell.first.first)+","+to_string(cell.first.second), (cell.first.first));
//                y_vec1.push_back(cell.first.second);
//                y_shift1.push_back(get<4>(cell.second)->_y);
//                y_shift.push_back(get<4>(cell.second)->_y);
//                y_uav1.add_val(to_string(cell.first.first)+","+to_string(cell.first.second),(get<4>(cell.second)->_y));
//                y1.add_val(to_string(cell.first.first)+","+to_string(cell.first.second), (cell.first.second));
//                z_vec1.push_back(av_z);
//                zmin_vec1.push_back(min_z);
//                zmax_vec1.push_back(max_z);
//                z_shift1.push_back(get<4>(cell.second)->_height);
//                z_shift.push_back(get<4>(cell.second)->_height);
//                z_uav1.add_val(to_string(cell.first.first)+","+to_string(cell.first.second),get<4>(cell.second)->_height);
//                z1.add_val(to_string(cell.first.first)+","+to_string(cell.first.second), max_z);
//                uav_roll1.push_back(get<4>(cell.second)->_roll);
//                uav_pitch1.push_back(get<4>(cell.second)->_pitch);
//                uav_yaw1.push_back(get<4>(cell.second)->_yaw);
//                if(max_nb_dots<nb_dots)
//                    max_nb_dots = nb_dots;
//                av_nb_dots += nb_dots;
//                z_range = max_z - min_z;
//                DebugOff("z range = " << z_range << endl);
//                av_z_range += z_range;
//                if(z_range > max_z_range)
//                    max_z_range = z_range;
//                DebugOff("z min = " << min_z << endl);
//                DebugOff("z max = " << max_z << endl);
//                DebugOff("z av = " << av_z << endl);
//                DebugOff("nb points = " << nb_dots << endl << endl);
//            }
//            av_nb_dots /= grid1.size();
//            av_z_range /= grid1.size();
//            DebugOn("Max z range for a cell = " << max_z_range << endl);
//            DebugOn("Average z range for a cell = " << av_z_range << endl);
//            DebugOn("Average points per cell = " << av_nb_dots << endl);
//            DebugOn("Max points per cell = " << max_nb_dots << endl << endl);
//    //        LASwriteOpener laswriteopener;
//    //        laswriteopener.set_file_name("compressed.laz");
//    //        LASwriter* laswriter = laswriteopener.open(&lasreader->header);
//    //
//    //        while (lasreader->read_point()) laswriter->write_point(&lasreader->point);
//    //
//    //        laswriter->close();
//    //        delete laswriter;
//
//            lasreader->close();
//            delete lasreader;
//        }
//
//        int xdim2=0, ydim2=0, zdim2=0;
//
//        lasreadopener.set_file_name(LiDAR_file2.c_str());
//        lasreadopener.set_populate_header(TRUE);
//        if (!lasreadopener.active())
//        {
//            throw invalid_argument("ERROR: no input specified\n");
//        }
//        while (lasreadopener.active())
//        {
//            LASreader* lasreader = lasreadopener.open();
//            if (lasreader == 0)
//            {
//                throw invalid_argument("ERROR: could not open lasreader for second file\n");
//            }
//
//            xdim1 = ceil(lasreader->header.max_x*grid_step - lasreader->header.min_x*grid_step);
//            ydim1 = ceil(lasreader->header.max_y*grid_step - lasreader->header.min_y*grid_step);
//            ground_z = floor(lasreader->header.min_z);
//            auto nb_pts_per_frame = lasreader->npoints/frames.size();
//            DebugOn("Number of points = " << lasreader->npoints << endl);
//            DebugOn("Number of points per frame = " <<  nb_pts_per_frame << endl);
//            DebugOn("min x axis = " << lasreader->header.min_x << endl);
//            DebugOn("max x axis = " << lasreader->header.max_x << endl);
//            DebugOn("min y axis = " << lasreader->header.min_y << endl);
//            DebugOn("max y axis = " << lasreader->header.max_y << endl);
//            DebugOn("dimension in x axis = " << xdim1 << endl);
//            DebugOn("dimension in y axis = " << ydim1 << endl);
//
//
//            int nb_dots; /* Number of measurements inside cell */
//            double xpos, ypos;
//            double z, min_z, max_z, av_z;
//            pair<double,double> pos;
//            size_t nb_pts = 0;
//            tuple<int,double,double,double,UAVPoint*> cell; /* <nb_dots,min_z,max_z,av_z> */
//            /* Now get rid of the first points
//             lasreader->read_point();
//             auto gps_time = lasreader->point.get_gps_time();
//             while (lasreader->point.get_gps_time()==gps_time)
//             {
//             lasreader->read_point();
//             }
//             */
//            auto frames_it = frames.begin();
//            while (lasreader->read_point())
//            {
//                if(nb_pts==0){
//                    DebugOn(to_string_with_precision(10.*(lasreader->point.get_gps_time()+315964800. - 18.),24) << ": (" << to_string_with_precision(lasreader->point.get_x(),10) <<"," << to_string_with_precision(lasreader->point.get_y(),10) << ","<< to_string_with_precision(lasreader->point.get_z(),10) <<")"<<endl);
//                    //                return 0;
//                }
//                xpos = (((lasreader->point.get_x())));
//                ypos = (((lasreader->point.get_y())));
//                if (ypos < 0)
//                {
//                    fprintf(stderr, "ERROR: ypos = %d\n", ypos);
//                }
//
//
//                //            z = lasreader->point.get_Z();
//                //            z = floor((lasreader->point.get_z()*100*grid_step - lasreader->header.min_z*100*grid_step));
//                z = lasreader->point.get_z();
//                timestamps.insert(lasreader->point.get_gps_time()+1./(nb_pts%10));
//                DebugOff("GPS time = " << to_string_with_precision(lasreader->point.get_gps_time(),24) << endl);
//                LidarPoints.push_back(LidarPoint(lasreader->point.get_gps_time(),xpos,ypos,z));
//                if(nb_pts>=nb_pts_per_frame && nb_pts%nb_pts_per_frame==0 && next(frames_it)!=frames.end()){
//                    frames_it++;
//                }
//                auto frame_ptr = frames.insert(make_pair(frames_it->first, Frame(LidarPoints.back()._unix_time*10+(nb_pts%10))));
//    //            auto frame_ptr = frames.insert(make_pair(LidarPoints.back()._unix_time*10+(nb_pts%10), Frame(LidarPoints.back()._unix_time*10+(nb_pts%10))));
//                frame_ptr.first->second.add_lidar_point(LidarPoints.back());
//                if(frame_ptr.second){
//                    throw invalid_argument("Frame with missing UAV point at " + to_string(LidarPoints.back()._unix_time*10+(nb_pts%10)));
//                }
//                else{
//                    LidarPoints.back()._uav_pt = frame_ptr.first->second._uav_point;
//                    frames1.insert(make_pair(LidarPoints.back()._unix_time*10+(nb_pts%10), Frame(frame_ptr.first->second)));
//                    //                uav_x1.push_back((frame_ptr.first->second._uav_points.front()->_longitude+582690.8242)*1e-5);
//                    //                uav_y1.push_back((frame_ptr.first->second._uav_points.front()->_latitude+4107963.58)*1e-5);
//                    //                uav_z1.push_back(frame_ptr.first->second._uav_points.front()->_height*100);
//                }
//                DebugOff("Added lidar point to frame at " << LidarPoints.back()._unix_time << " : " << LidarPoints.back()._hour << ":" << LidarPoints.back()._minutes << ":" << LidarPoints.back()._seconds << endl);
//                //            xpos = floor((lasreader->point.get_x()*grid_step - lasreader->header.min_x*grid_step));
//                //            ypos = floor((lasreader->point.get_y()*grid_step - lasreader->header.min_y*grid_step));
//
//                //            z = round((lasreader->point.get_z()*100*grid_step));
//                //            xpos -= frame_ptr.first->second._uav_point->_x;
//                //            ypos -= frame_ptr.first->second._uav_point->_y;
//                //            z -= frame_ptr.first->second._uav_point->_height;
//                pos = make_pair(xpos,ypos);
//                cell = make_tuple(1,z,z,z,frame_ptr.first->second._uav_point);
//                DebugOff("new z value = " << z << endl);
//                DebugOff("fractional z value = " << lasreader->point.get_z() << endl);
//                DebugOff("z range = " << grid_maxz[pos] - grid_minz[pos] << endl << endl);
//                auto map_pair = grid2.insert(make_pair(pos,cell));
//                auto map_pair_all = grid.insert(make_pair(pos,cell));
//                if (!map_pair.second) {/* Cell already exists */
//                    //                DebugOn("Found overlapping cell after " << nb_pts << endl);
//                    //                DebugOn("(" << to_string_with_precision(lasreader->point.get_x(),10) <<"," << to_string_with_precision(lasreader->point.get_y(),10) << ","<< to_string_with_precision(lasreader->point.get_z(),10) <<")"<<endl);
//                    //                                return 0;
//                    nb_dots = get<0>(map_pair.first->second);
//                    min_z = get<1>(map_pair.first->second);
//                    max_z = get<2>(map_pair.first->second);
//                    av_z = get<3>(map_pair.first->second);
//                    /* Update cell data */
//                    get<0>(map_pair.first->second)++;
//                    if(min_z>z){
//                        DebugOff("updating  min z at cell [" << xpos << "," << ypos << "]" << endl);
//                        DebugOff("previous min z = " << min_z << endl);
//                        DebugOff("new min z = " << z << endl);
//                        get<1>(map_pair.first->second) = z;
//                    }
//                    if(max_z<z){
//                        DebugOff("updating  max z at cell [" << xpos << "," << ypos << "]" << endl);
//                        DebugOff("previous max z = " << max_z << endl);
//                        DebugOff("new max z = " << z << endl);
//                        get<2>(map_pair.first->second) = z;
//                    }
//                    av_z *= nb_dots;/* remove previous denominator */
//                    get<3>(map_pair.first->second) = (av_z + z)/(nb_dots+1); /* Update denominator */
//                }
//                if (!map_pair_all.second) {/* Cell already exists */
//                    nb_dots = get<0>(map_pair_all.first->second);
//                    min_z = get<1>(map_pair_all.first->second);
//                    max_z = get<2>(map_pair_all.first->second);
//                    av_z = get<3>(map_pair_all.first->second);
//                    /* Update cell data */
//                    get<0>(map_pair_all.first->second)++;
//                    if(min_z>z){
//                        get<1>(map_pair_all.first->second) = z;
//                    }
//                    if(max_z<z){
//                        get<2>(map_pair_all.first->second) = z;
//                    }
//                    av_z *= nb_dots;/* remove previous denominator */
//                    get<3>(map_pair_all.first->second) = (av_z + z)/(nb_dots+1); /* Update denominator */
//                }
//                nb_pts++;
//            }
//            DebugOn("Read " << nb_pts << " points in file 2" << endl);
//            DebugOn("2D grid2 has " << grid2.size() << " cells " << endl);
//            DebugOn("xdim2 = " << xdim2 << endl);
//            DebugOn("ydim2 = " << ydim2 << endl);
//
//            int cell_counter = 0, max_nb_dots = 0;
//            double av_nb_dots = 0, z_range = 0,  max_z_range = 0, av_z_range = 0;
//            for (auto &cell : grid2) {
//                DebugOff("Cell num " << cell_counter++ << " at [" << cell.first.first << "," << cell.first.second << "]" << endl);
//                nb_dots = get<0>(cell.second);
//                min_z = get<1>(cell.second);
//                max_z = get<2>(cell.second);
//                av_z = get<3>(cell.second);
//                x_vec2.push_back(cell.first.first);
//                x_shift2.push_back(get<4>(cell.second)->_x);
//                x_shift.push_back(get<4>(cell.second)->_x);
//                x_uav2.add_val(to_string(cell.first.first)+","+to_string(cell.first.second),(get<4>(cell.second)->_x));
//                x2.add_val(to_string(cell.first.first)+","+to_string(cell.first.second), (cell.first.first));
//                y_vec2.push_back(cell.first.second);
//                y_shift2.push_back(get<4>(cell.second)->_y);
//                y_shift.push_back(get<4>(cell.second)->_y);
//                y_uav2.add_val(to_string(cell.first.first)+","+to_string(cell.first.second),(get<4>(cell.second)->_y));
//                y2.add_val(to_string(cell.first.first)+","+to_string(cell.first.second), (cell.first.second));
//                z_vec2.push_back(av_z);
//                zmin_vec2.push_back(min_z);
//                zmax_vec2.push_back(max_z);
//                z_shift2.push_back(get<4>(cell.second)->_height);
//                z_shift.push_back(get<4>(cell.second)->_height);
//                z_uav2.add_val(to_string(cell.first.first)+","+to_string(cell.first.second),get<4>(cell.second)->_height);
//                z2.add_val(to_string(cell.first.first)+","+to_string(cell.first.second), max_z);
//                uav_roll2.push_back(get<4>(cell.second)->_roll);
//                uav_pitch2.push_back(get<4>(cell.second)->_pitch);
//                uav_yaw2.push_back(get<4>(cell.second)->_yaw);
//                if(max_nb_dots<nb_dots)
//                    max_nb_dots = nb_dots;
//                av_nb_dots += nb_dots;
//                z_range = max_z - min_z;
//                DebugOff("z range = " << z_range << endl);
//                av_z_range += z_range;
//                if(z_range > max_z_range)
//                    max_z_range = z_range;
//                DebugOff("z min = " << min_z << endl);
//                DebugOff("z max = " << max_z << endl);
//                DebugOff("z av = " << av_z << endl);
//                DebugOff("nb points = " << nb_dots << endl << endl);
//            }
//            av_nb_dots /= grid2.size();
//            av_z_range /= grid2.size();
//            DebugOn("2D grid has " << grid2.size() << " cells " << endl);
//            DebugOn("Max z range for a cell in grid2 = " << max_z_range << endl);
//            DebugOn("Average z range for a cell in grid2 = " << av_z_range << endl);
//            DebugOn("Average points per cell in grid2 = " << av_nb_dots << endl);
//            DebugOn("Max points per cell in grid2 = " << max_nb_dots << endl << endl);
//            //        LASwriteOpener laswriteopener;
//            //        laswriteopener.set_file_name("compressed.laz");
//            //        LASwriter* laswriter = laswriteopener.open(&lasreader->header);
//            //
//            //        while (lasreader->read_point()) laswriter->write_point(&lasreader->point);
//            //
//            //        laswriter->close();
//            //        delete laswriter;
//
//            lasreader->close();
//            delete lasreader;
//        }
//        int cell_counter = 0, max_nb_dots = 0, nb_dots = 0;
//        double av_nb_dots = 0, av_z_range = 0, max_z_range = 0, z_range = 0, min_z = 0, max_z = 0, av_z = 0;
//        indices cells("cells");
//        int nb_overlap = 0;
//        bool insert = true;
//        for (auto &cell : grid) {
//            DebugOff("Cell num " << cell_counter++ << " at [" << cell.first.first << "," << cell.first.second << "]" << endl);
//            if(grid1.count(cell.first)>0 && grid2.count(cell.first)){
//                nb_overlap++;
//    //            if(insert){
//                if(nb_overlap%30==1){
//                    cells.insert(to_string(cell.first.first) + "," + to_string(cell.first.second));
//                }
//                insert = !insert;
//                nb_dots = get<0>(cell.second);
//                min_z = get<1>(cell.second);
//                max_z = get<2>(cell.second);
//                av_z = get<3>(cell.second);
//                if(max_nb_dots<nb_dots)
//                    max_nb_dots = nb_dots;
//                av_nb_dots += nb_dots;
//                z_range = max_z - min_z;
//                DebugOff("z range = " << z_range << endl);
//                av_z_range += z_range;
//                if(z_range > max_z_range)
//                    max_z_range = z_range;
//                DebugOff("z min = " << min_z << endl);
//                DebugOff("z max = " << max_z << endl);
//                DebugOff("z av = " << av_z << endl);
//                DebugOff("nb points = " << nb_dots << endl << endl);
//            }
//        }
//        av_nb_dots /= nb_overlap;
//        av_z_range /= nb_overlap;
//        DebugOn("Number of overlapping points = " << nb_overlap << endl);
//        DebugOn("Max z range for a cell in combined grid = " << max_z_range << endl);
//        DebugOn("Average z range for a cell in combined grid = " << av_z_range << endl);
//        DebugOn("Average points per cell in combined grid = " << av_nb_dots << endl);
//        DebugOn("Max points per cell in combined grid = " << max_nb_dots << endl << endl);
//            int max_nb_points = 0, av_nb_points = 0, nb_empty_frames = 0;
//            for (const auto &frame: frames) {
//                int nb_pts = frame.second._lidar_points->size();
//                if(nb_pts==0)
//                    nb_empty_frames++;
//                else{
//                    max_nb_points = std::max(max_nb_points, nb_pts);
//                    av_nb_points += nb_pts;
//                }
//            }
//            av_nb_points/=(frames.size()-nb_empty_frames);
//            DebugOn("Number of empty frames = " << nb_empty_frames << endl);
//            DebugOn("Average number of points per frame = " << av_nb_points << endl);
//            DebugOn("Max number of points per frame = " << max_nb_points << endl);
//            DebugOn("Number of timestamps = " << timestamps.size() << endl);
//    //
//    //        //        LASwriteOpener laswriteopener;
//    //        //        laswriteopener.set_file_name("compressed.laz");
//    //        //        LASwriter* laswriter = laswriteopener.open(&lasreader->header);
//    //        //
//    //        //        while (lasreader->read_point()) laswriter->write_point(&lasreader->point);
//    //        //
//    //        //        laswriter->close();
//    //        //        delete laswriter;
//    //
//    //        lasreader->close();
//    //        delete lasreader;
//    //    }
//        for (const auto &frame: frames1) {
//            uav_x1.push_back((frame.second._uav_point->_x));
//            uav_y1.push_back((frame.second._uav_point->_y));
//            uav_z1.push_back(frame.second._uav_point->_height);
//        }
//
//
//
//
//        using namespace gravity;
//        double roll_1 = 0, pitch_1 = 0, yaw_1 = 0;
//        double roll_2 = 0, pitch_2 = 0, yaw_2 = 0;
//
//        bool solve_lidar = true;
//        if (solve_lidar) {
//            Model<> M("Lidar");
//            var<> new_x1("new_x1"), new_y1("new_y1"), new_z1("new_z1");
//            var<> new_x2("new_x2"), new_y2("new_y2"), new_z2("new_z2");
//            var<> x_diff("x_diff", pos_), y_diff("y_diff", pos_), z_diff("z_diff", pos_);
//            var<> yaw1("yaw1", 0, 0), pitch1("pitch1", -0.1, 0.1), roll1("roll1", 0, 0);
//            var<> yaw2("yaw2", 0, 0), pitch2("pitch2", -0.1, 0.1), roll2("roll2", 0, 0);
//
//            if(!optimize){
//    //            yaw1.set_lb(0.25*pi/180);
//    //            yaw1.set_ub(0.25*pi/180);
//    //            pitch1.set_lb(-0.9*pi/180);
//    //            pitch1.set_ub(-0.9*pi/180);
//    //            roll1.set_lb(1.45*pi/180.);
//    //            roll1.set_ub(1.45*pi/180.);
//    //            yaw2.set_lb(-0.25*pi/180);
//    //            yaw2.set_ub(-0.25*pi/180);
//    //            pitch2.set_lb(0.9*pi/180);
//    //            pitch2.set_ub(0.9*pi/180);
//    //            roll2.set_lb(-1.45*pi/180.);
//    //            roll2.set_ub(-1.45*pi/180.);
//                yaw1.set_lb(0);
//                yaw1.set_ub(0);
//                pitch1.set_lb(0);
//                pitch1.set_ub(0);
//                roll1.set_lb(0);
//                roll1.set_ub(0);
//
//                yaw2.set_lb(0);
//                yaw2.set_ub(0);
//                pitch2.set_lb(0);
//                pitch2.set_ub(0);
//                roll2.set_lb(0);
//                roll2.set_ub(0);
//            }
//        //    var<> sin_yaw1("yaw", -1, 1), sin_pitch1("pitch", -1, 1), sin_roll1("roll", -1, 1);
//        //    var<> sin_yaw2("yaw", -1, 1), sin_pitch2("pitch", -0.1, 0.1), roll2("roll", -0.1, 0.1);
//
//        //
//        //    var<> yaw1("yaw1", 0, 0), pitch1("pitch1", 0, 0), roll1("roll1", 0, 0);
//        //    var<> yaw2("yaw2", 0, 0), pitch2("pitch2", 0, 0), roll2("roll2", 0, 0);
//
//        //    var<> abs_yaw1("abs_yaw1", 0, 0.1), abs_pitch1("abs_pitch1", 0, 0.1), abs_roll1("abs_roll1", 0, 0.1);
//        //    var<> abs_yaw2("abs_yaw2", 0, 0.2), abs_pitch2("abs_pitch2", 0, 0.2), abs_roll2("abs_roll2", 0, 0.2);
//        //    var<> yaw1("yaw1"), pitch1("pitch1"), roll1("roll1");
//        //    var<> yaw2("yaw2"), pitch2("pitch2"), roll2("roll2");
//            M.add(yaw1.in(R(1)),pitch1.in(R(1)),roll1.in(R(1)));
//            M.add(yaw2.in(R(1)),pitch2.in(R(1)),roll2.in(R(1)));
//        //    M.add(abs_yaw1.in(R(1)),abs_pitch1.in(R(1)),abs_roll1.in(R(1)));
//        //    M.add(abs_yaw2.in(R(1)),abs_pitch2.in(R(1)),abs_roll2.in(R(1)));
//        //    yaw1 = 0.1; yaw2 = -0.1; pitch1 = 0.1; pitch2 = 0.1; roll1 = 0.1; roll2 = -0.1;
//        //    yaw1 = -0.01; roll1 = 0.001; pitch1 = 0.001;
//        //    pitch1 = 0.017;
//            M.add(new_x1.in(cells), new_y1.in(cells), new_z1.in(cells));
//            M.add(new_x2.in(cells), new_y2.in(cells), new_z2.in(cells));
//    //        M.add(x_diff.in(cells), y_diff.in(cells), z_diff.in(cells));
//            M.add(z_diff.in(cells));
//
//            Constraint<> Equal_pitch("Equal_pitch");
//            Equal_pitch += pitch1 - pitch2;
//            M.add(Equal_pitch==0);
//    //
//            Constraint<> Opp_roll("Opp_roll");
//            Opp_roll += roll1 + roll2;
//            M.add(Opp_roll==0);
//
//            Constraint<> Opp_yaw("Opp_yaw");
//            Opp_yaw += yaw1 + yaw2;
//            M.add(Opp_yaw==0);
//
//
//        //    Constraint<> x_norm("x_norm");
//        //    x_norm += x_diff - pow((new_x1 - new_x2),2);
//        //    M.add(x_norm.in(cells)>=0);
//        //
//        //    Constraint<> y_norm("y_norm");
//        //    y_norm += y_diff - pow((new_y1 - new_y2),2);
//        //    M.add(y_norm.in(cells)>=0);
//        //
//        //    Constraint<> z_norm("z_norm");
//        //    z_norm += z_diff - pow((new_z1 - new_z2),2);
//        //    M.add(z_norm.in(cells)>=0);
//
//        //    Constraint<> yaw_abs1("yaw_abs1");
//        //    yaw_abs1 += abs_yaw1 - yaw2;
//        //    M.add(yaw_abs1 >= 0);
//        //
//        //    Constraint<> yaw_abs2("yaw_abs2");
//        //    yaw_abs2 += abs_yaw1 + yaw2;
//        //    M.add(yaw_abs2 >= 0);
//        //
//        //
//        //    Constraint<> roll_abs1("roll_abs1");
//        //    roll_abs1 += abs_roll1 - roll2;
//        //    M.add(roll_abs1 >= 0);
//        //
//        //    Constraint<> roll_abs2("roll_abs2");
//        //    roll_abs2 += abs_roll1 + roll2;
//        //    M.add(roll_abs2 >= 0);
//        //
//        //    Constraint<> pitch_abs1("pitch_abs1");
//        //    pitch_abs1 += abs_pitch1 - pitch2;
//        //    M.add(pitch_abs1 >= 0);
//        //
//        //    Constraint<> pitch_abs2("pitch_abs2");
//        //    pitch_abs2 += abs_pitch1 + pitch2;
//        //    M.add(pitch_abs2 >= 0);
//
//
//    //        Constraint<> x_abs1("x_abs1");
//    //        x_abs1 += x_diff - (new_x1 - new_x2);
//    //        M.add(x_abs1.in(cells)>=0);
//    //
//    //        Constraint<> x_abs2("x_abs2");
//    //        x_abs2 += x_diff - (new_x2 - new_x1);
//    //        M.add(x_abs2.in(cells)>=0);
//    //
//    //
//    //        Constraint<> y_abs1("y_abs1");
//    //        y_abs1 += y_diff - (new_y1 - new_y2);
//    //        M.add(y_abs1.in(cells)>=0);
//    //
//    //        Constraint<> y_abs2("y_abs2");
//    //        y_abs2 += y_diff - (new_y2 - new_y1);
//    //        M.add(y_abs2.in(cells)>=0);
//
//            Constraint<> z_abs1("z_abs1");
//            z_abs1 += z_diff - (new_z1 - new_z2);
//            M.add(z_abs1.in(cells)>=0);
//
//            Constraint<> z_abs2("z_abs2");
//            z_abs2 += z_diff - (new_z2 - new_z1);
//            M.add(z_abs2.in(cells)>=0);
//
//            auto ids = yaw1.repeat_id(cells.size());
//
//            /* alpha = yaw_, beta = pitch_ and gamma = roll_ */
//            Constraint<> x_rot1("x_rot1");
//            x_rot1 += new_x1 - x_uav1.in(cells);
//    //        x_rot1 -= (x1.in(cells)-x_uav1.in(cells))*cos(yaw1.in(ids))*cos(pitch1.in(ids)) + (y1.in(cells)-y_uav1.in(cells))*(cos(yaw1.in(ids))*sin(pitch1.in(ids))*sin(roll1.in(ids)) - sin(yaw1.in(ids))*cos(roll1.in(ids))) + (z1.in(cells)-z_uav1.in(cells))*(cos(yaw1.in(ids))*sin(pitch1.in(ids))*cos(roll1.in(ids)) + sin(yaw1.in(ids))*sin(roll1.in(ids)));
//    //        M.add(x_rot1.in(cells)==0);
//
//            x_rot1 -= (x1.in(cells)-x_uav1.in(cells))*cos(yaw1.in(ids))*cos(roll1.in(ids)) + (y1.in(cells)-y_uav1.in(cells))*(cos(yaw1.in(ids))*sin(roll1.in(ids))*sin(pitch1.in(ids)) - sin(yaw1.in(ids))*cos(pitch1.in(ids))) + (z1.in(cells)-z_uav1.in(cells))*(cos(yaw1.in(ids))*sin(roll1.in(ids))*cos(pitch1.in(ids)) + sin(yaw1.in(ids))*sin(pitch1.in(ids)));
//            M.add(x_rot1.in(cells)==0);
//
//        //    Constraint<> x_rot2("x_rot2");
//        //    x_rot2 += new_x2 + x_uav2.in(cells);
//        //    x_rot2 -= (x2.in(cells)-x_uav2.in(cells))*cos(yaw1.in(ids))*cos(pitch1.in(ids)) + (y2.in(cells)-y_uav2.in(cells))*(cos(yaw1.in(ids))*sin(pitch1.in(ids))*sin(roll1.in(ids)) - sin(yaw1.in(ids))*cos(roll1.in(ids))) + (z2.in(cells)-z_uav2.in(cells))*(cos(yaw1.in(ids))*sin(pitch1.in(ids))*cos(roll1.in(ids)) + sin(yaw1.in(ids))*sin(roll1.in(ids)));
//        //    M.add(x_rot2.in(cells)==0);
//
//    //        Constraint<> x_rot2("x_rot2");
//    //        x_rot2 += new_x2 - x_uav2.in(cells);
//    //        x_rot2 -= (x2.in(cells)-x_uav2.in(cells))*cos(yaw2.in(ids))*cos(pitch2.in(ids)) + (y2.in(cells)-y_uav2.in(cells))*(cos(yaw2.in(ids))*sin(pitch2.in(ids))*sin(roll2.in(ids)) - sin(yaw2.in(ids))*cos(roll2.in(ids))) + (z2.in(cells)-z_uav2.in(cells))*(cos(yaw2.in(ids))*sin(pitch2.in(ids))*cos(roll2.in(ids)) + sin(yaw2.in(ids))*sin(roll2.in(ids)));
//    //        M.add(x_rot2.in(cells)==0);
//
//            Constraint<> x_rot2("x_rot2");
//            x_rot2 += new_x2 - x_uav2.in(cells);
//            x_rot2 -= (x2.in(cells)-x_uav2.in(cells))*cos(yaw2.in(ids))*cos(roll2.in(ids)) + (y2.in(cells)-y_uav2.in(cells))*(cos(yaw2.in(ids))*sin(roll2.in(ids))*sin(pitch2.in(ids)) - sin(yaw2.in(ids))*cos(pitch2.in(ids))) + (z2.in(cells)-z_uav2.in(cells))*(cos(yaw2.in(ids))*sin(roll2.in(ids))*cos(pitch2.in(ids)) + sin(yaw2.in(ids))*sin(pitch2.in(ids)));
//            M.add(x_rot2.in(cells)==0);
//
//    //        Constraint<> y_rot1("y_rot1");
//    //        y_rot1 += new_y1 - y_uav1.in(cells);
//    //        y_rot1 -= (x1.in(cells)-x_uav1.in(cells))*sin(yaw1.in(ids))*cos(pitch1.in(ids)) + (y1.in(cells)-y_uav1.in(cells))*(sin(yaw1.in(ids))*sin(pitch1.in(ids))*sin(roll1.in(ids)) + cos(yaw1.in(ids))*cos(roll1.in(ids))) + (z1.in(cells)-z_uav1.in(cells))*(sin(yaw1.in(ids))*sin(pitch1.in(ids))*cos(roll1.in(ids)) - cos(yaw1.in(ids))*sin(roll1.in(ids)));
//    //        M.add(y_rot1.in(cells)==0);
//
//            Constraint<> y_rot1("y_rot1");
//            y_rot1 += new_y1 - y_uav1.in(cells);
//            y_rot1 -= (x1.in(cells)-x_uav1.in(cells))*sin(yaw1.in(ids))*cos(roll1.in(ids)) + (y1.in(cells)-y_uav1.in(cells))*(sin(yaw1.in(ids))*sin(roll1.in(ids))*sin(pitch1.in(ids)) + cos(yaw1.in(ids))*cos(pitch1.in(ids))) + (z1.in(cells)-z_uav1.in(cells))*(sin(yaw1.in(ids))*sin(roll1.in(ids))*cos(pitch1.in(ids)) - cos(yaw1.in(ids))*sin(pitch1.in(ids)));
//            M.add(y_rot1.in(cells)==0);
//
//    //        Constraint<> y_rot2("y_rot2");
//    //        y_rot2 += new_y2 - y_uav2.in(cells);
//    //        y_rot2 -= (x2.in(cells)-x_uav2.in(cells))*sin(yaw2.in(ids))*cos(pitch2.in(ids)) + (y2.in(cells)-y_uav2.in(cells))*(sin(yaw2.in(ids))*sin(pitch2.in(ids))*sin(roll2.in(ids)) + cos(yaw2.in(ids))*cos(roll2.in(ids))) + (z2.in(cells)-z_uav2.in(cells))*(sin(yaw2.in(ids))*sin(pitch2.in(ids))*cos(roll2.in(ids)) - cos(yaw2.in(ids))*sin(roll2.in(ids)));
//    //        M.add(y_rot2.in(cells)==0);
//
//            Constraint<> y_rot2("y_rot2");
//            y_rot2 += new_y2 - y_uav2.in(cells);
//            y_rot2 -= (x2.in(cells)-x_uav2.in(cells))*sin(yaw2.in(ids))*cos(roll2.in(ids)) + (y2.in(cells)-y_uav2.in(cells))*(sin(yaw2.in(ids))*sin(roll2.in(ids))*sin(pitch2.in(ids)) + cos(yaw2.in(ids))*cos(pitch2.in(ids))) + (z2.in(cells)-z_uav2.in(cells))*(sin(yaw2.in(ids))*sin(roll2.in(ids))*cos(pitch2.in(ids)) - cos(yaw2.in(ids))*sin(pitch2.in(ids)));
//            M.add(y_rot2.in(cells)==0);
//
//        //    Constraint<> y_rot2("y_rot2");
//        //    y_rot2 += new_y2 + y_uav2.in(cells);
//        //    y_rot2 -= (x2.in(cells)-x_uav2.in(cells))*sin(yaw1.in(ids))*cos(pitch1.in(ids)) + (y2.in(cells)-y_uav2.in(cells))*(sin(yaw1.in(ids))*sin(pitch1.in(ids))*sin(roll1.in(ids)) + cos(yaw1.in(ids))*cos(roll1.in(ids))) + (z2.in(cells)-z_uav2.in(cells))*(sin(yaw1.in(ids))*sin(pitch1.in(ids))*cos(roll1.in(ids)) - cos(yaw1.in(ids))*sin(roll1.in(ids)));
//        //    M.add(y_rot2.in(cells)==0);
//
//
//    //        Constraint<> z_rot1("z_rot1");
//    //        z_rot1 += new_z1 - z_uav1.in(cells);
//    //        z_rot1 -= (x1.in(cells)-x_uav1.in(cells))*sin(-1*pitch1.in(ids)) + (y1.in(cells)-y_uav1.in(cells))*(cos(pitch1.in(ids))*sin(roll1.in(ids))) + (z1.in(cells)-z_uav1.in(cells))*(cos(pitch1.in(ids))*cos(roll1.in(ids)));
//    //        M.add(z_rot1.in(cells)==0);
//
//            Constraint<> z_rot1("z_rot1");
//            z_rot1 += new_z1 - z_uav1.in(cells);
//            z_rot1 -= (x1.in(cells)-x_uav1.in(cells))*sin(-1*roll1.in(ids)) + (y1.in(cells)-y_uav1.in(cells))*(cos(roll1.in(ids))*sin(pitch1.in(ids))) + (z1.in(cells)-z_uav1.in(cells))*(cos(roll1.in(ids))*cos(pitch1.in(ids)));
//            M.add(z_rot1.in(cells)==0);
//
//    //        Constraint<> z_rot2("z_rot2");
//    //        z_rot2 += new_z2 - z_uav2.in(cells);
//    //        z_rot2 -= (x2.in(cells)-x_uav2.in(cells))*sin(-1*pitch2.in(ids)) + (y2.in(cells)-y_uav2.in(cells))*(cos(pitch2.in(ids))*sin(roll2.in(ids))) + (z2.in(cells)-z_uav2.in(cells))*(cos(pitch2.in(ids))*cos(roll2.in(ids)));
//    //        M.add(z_rot2.in(cells)==0);
//
//            Constraint<> z_rot2("z_rot2");
//            z_rot2 += new_z2 - z_uav2.in(cells);
//            z_rot2 -= (x2.in(cells)-x_uav2.in(cells))*sin(-1*roll2.in(ids)) + (y2.in(cells)-y_uav2.in(cells))*(cos(roll2.in(ids))*sin(pitch2.in(ids))) + (z2.in(cells)-z_uav2.in(cells))*(cos(roll2.in(ids))*cos(pitch2.in(ids)));
//            M.add(z_rot2.in(cells)==0);
//
//        //    Constraint<> z_rot2("z_rot2");
//        //    z_rot2 += new_z2 + z_uav2.in(cells);
//        //    z_rot2 -= (x2.in(cells)-x_uav2.in(cells))*sin(-1*pitch1.in(ids)) + (y2.in(cells)-y_uav2.in(cells))*(cos(pitch1.in(ids))*sin(roll1.in(ids))) + (z2.in(cells)-z_uav2.in(cells))*(cos(pitch1.in(ids))*cos(roll1.in(ids)));
//        //    M.add(z_rot2.in(cells)==0);
//
//        //    M.min(sum(z_diff)/nb_overlap);
//
//    //        M.min(sum(z_diff));
//            M.min(sum(z_diff)/cells.size());
//
//    //        M.min(sum(x_diff)/cells.size() + sum(y_diff)/cells.size() + sum(z_diff)/cells.size());
//
//        //    M.print();
//
//            solver<> S(M,ipopt);
//            S.run();
//
//
//    //        for (int i = 0; i<500; i++) {
//    //            pre_x.add_val(x_rot1.eval(i));
//    //            pre_y.add_val(y_rot1.eval(i));
//    //            pre_z.add_val(z_rot1.eval(i));
//    //            x_uav.add_val(x_uav1.eval(i));
//    //            y_uav.add_val(y_uav1.eval(i));
//    //            z_uav.add_val(z_uav1.eval(i));
//    //        }
//    //        for (int i = 0; i<500; i++) {
//    //            pre_x.add_val(x_rot2.eval(i));
//    //            pre_y.add_val(y_rot2.eval(i));
//    //            pre_z.add_val(z_rot2.eval(i));
//    //            x_uav.add_val(x_uav2.eval(i));
//    //            y_uav.add_val(y_uav2.eval(i));
//    //            z_uav.add_val(z_uav2.eval(i));
//    //        }
//    //    M.print_solution();
//
//            pitch_1 = pitch1.eval();
//            roll_1 = roll1.eval();
//            yaw_1 = yaw1.eval();
//        //    pitch_2 = pitch1.eval();
//        //    roll_2 = roll1.eval();
//        //    yaw_2 = yaw1.eval();
//            pitch_2 = pitch2.eval();
//            roll_2 = roll2.eval();
//            yaw_2 = yaw2.eval();
//        }
//        DebugOn("Pitch1 = " << pitch_1*180/pi << endl);
//        DebugOn("Roll1 = " << roll_1*180/pi << endl);
//        DebugOn("Yaw1 = " << yaw_1*180/pi << endl);
//        DebugOn("Pitch2 = " << pitch_2*180/pi << endl);
//        DebugOn("Roll2 = " << roll_2*180/pi << endl);
//        DebugOn("Yaw2 = " << yaw_2*180/pi << endl);
//        //    roll_ = -0.55;
//        //    pitch_ = 0.01375;
//        //    pitch_ = 1.375*pi/180.;
//        //    roll_ = -0.55*pi/180.;
//        //    roll_ = -0.9*pi/180.;
//        //    roll_ = 1.375*pi/180.;
//        //    pitch_ = 0.9*pi/180.;
//        //    yaw_ = -0.55*pi/180.;
//        //    yaw_ = 0.9*pi/180.;
//        //    yaw_ = -0.55;
//        double alpha = 0, beta = 0, gamma = 0; /* alpha = yaw_, beta = pitch_ and gamma = roll_ */
//        //    alpha = -0.55;
//        //    beta = 0.9;
//        //    gamma = 1.375;
//        //    beta = 0.02;
//        //    gamma = 0.01375;
//        //    beta = 0.01;
//        //    gamma = 1.375;
//        //    alpha = -30.*pi/180.;
//        //    beta = 30.*pi/180.;
//        //    gamma = 30.*pi/180.;
//        auto tot_pts = x_vec1.size()+x_vec2.size();
//        x_combined.resize(tot_pts);
//    //    x_combined.resize(x_vec1.size());
//        y_combined.resize(tot_pts);
//        z_combined.resize(tot_pts);
//        zmin_combined.resize(tot_pts);
//        zmax_combined.resize(tot_pts);
//        double shifted_x, shifted_y, shifted_z;
//        beta = roll_1;
//        gamma = pitch_1;
//        alpha = yaw_1;
//        for (auto i = 0; i< x_vec1.size(); i++) {
//            //        gamma = roll_ + uav_roll1[i];
//            //        beta = pitch_ + uav_pitch1[i];
//            //        alpha = yaw_ + uav_yaw_1[i];
//            shifted_x = x_vec1[i] - x_shift1[i];
//            shifted_y = y_vec1[i] - y_shift1[i];
//            shifted_z = zmax_vec1[i] - z_shift1[i];
//            x_combined[i] = shifted_x*cos(alpha)*cos(beta) + shifted_y*(cos(alpha)*sin(beta)*sin(gamma) - sin(alpha)*cos(gamma)) + shifted_z*(cos(alpha)*sin(beta)*cos(gamma) + sin(alpha)*sin(gamma));
//            y_combined[i] = shifted_x*sin(alpha)*cos(beta) + shifted_y*(sin(alpha)*sin(beta)*sin(gamma) + cos(alpha)*cos(gamma)) + shifted_z*(sin(alpha)*sin(beta)*cos(gamma) - cos(alpha)*sin(gamma));
//            zmax_combined[i] = shifted_x*(-sin(beta)) + shifted_y*(cos(beta)*sin(gamma)) + shifted_z*(cos(beta)*cos(gamma));
//            x_combined[i] += x_shift1[i];
//            y_combined[i] += y_shift1[i];
//            zmax_combined[i] += z_shift1[i];
//
//            //        x_combined[i] = x_vec1[i]*cos(alpha)*cos(beta) + y_vec1[i]*(cos(alpha)*sin(beta)*sin(gamma) - sin(alpha)*cos(gamma)) + zmax_vec1[i]*(cos(alpha)*sin(beta)*cos(gamma) + sin(alpha)*sin(gamma));
//            //        y_combined[i] = x_vec1[i]*sin(alpha)*cos(beta) + y_vec1[i]*(sin(alpha)*sin(beta)*sin(gamma) + cos(alpha)*cos(gamma)) + zmax_vec1[i]*(sin(alpha)*sin(beta)*cos(gamma) - cos(alpha)*sin(gamma));
//            //        zmax_combined[i] = x_vec1[i]*(-sin(beta)) + y_vec1[i]*(cos(beta)*sin(gamma)) + zmax_vec1[i]*(cos(beta)*cos(gamma));
//            //        x_combined[i] = x_vec1[i]*cos(alpha)*cos(beta) + y_vec1[i]*sin(alpha)*cos(beta) + zmax_vec1[i]*(-sin(beta));
//            //        y_combined[i] = x_vec1[i]*(cos(alpha)*sin(beta)*sin(gamma) - sin(alpha)*cos(gamma)) + y_vec1[i]*(sin(alpha)*sin(beta)*sin(gamma) + cos(alpha)*cos(gamma)) + zmax_vec1[i]*(cos(beta)*sin(gamma));
//            //        zmax_combined[i] = x_vec1[i]*(cos(alpha)*sin(beta)*cos(gamma) + sin(alpha)*sin(gamma)) + y_vec1[i]*(sin(alpha)*sin(beta)*cos(gamma) - cos(alpha)*sin(gamma)) + zmax_vec1[i]*(cos(beta)*cos(gamma));
//
//            //        x_combined[i] = x_vec1[i];
//            //        y_combined[i] = y_vec1[i];
//            //        z_combined[i] = x_vec1[i]*(-sin(beta)) + y_vec1[i]*(cos(beta)*sin(gamma)) + z_vec1[i]*(cos(beta)*cos(gamma));
//            //        zmin_combined[i] = x_vec1[i]*(-sin(beta)) + y_vec1[i]*(cos(beta)*sin(gamma)) + zmin_vec1[i]*(cos(beta)*cos(gamma));
//
//        }
//        /* Moving to second flight line assumed to be in opposite direction */
//    //    roll_2 *= -1;
//    ////    pitch_ *= -1;
//    //    yaw_2 *= -1;
//        beta = roll_2;
//        gamma = pitch_2;
//        alpha = yaw_2;
//            for (auto i = 0; i< x_vec2.size(); i++) {
//        //        gamma = roll_ + uav_roll2[i];
//        //        beta = pitch_ + uav_pitch2[i];
//        //        alpha = yaw_ + uav_yaw2[i];
//                shifted_x = x_vec2[i] - x_shift2[i];
//                shifted_y = y_vec2[i] - y_shift2[i];
//                shifted_z = zmax_vec2[i] - z_shift2[i];
//                x_combined[x_vec1.size()+i] = shifted_x*cos(alpha)*cos(beta) + shifted_y*(cos(alpha)*sin(beta)*sin(gamma) - sin(alpha)*cos(gamma)) + shifted_z*(cos(alpha)*sin(beta)*cos(gamma) + sin(alpha)*sin(gamma));
//                y_combined[x_vec1.size()+i] = shifted_x*sin(alpha)*cos(beta) + shifted_y*(sin(alpha)*sin(beta)*sin(gamma) + cos(alpha)*cos(gamma)) + shifted_z*(sin(alpha)*sin(beta)*cos(gamma) - cos(alpha)*sin(gamma));
//                zmax_combined[x_vec1.size()+i] = shifted_x*(-sin(beta)) + shifted_y*(cos(beta)*sin(gamma)) + shifted_z*(cos(beta)*cos(gamma));
//                x_combined[x_vec1.size()+i] += x_shift2[i];
//                y_combined[x_vec1.size()+i] += y_shift2[i];
//                zmax_combined[x_vec1.size()+i] += z_shift2[i];
//            }
//
//        bool solve_post = false;
//        double post_roll = 0, post_pitch = 0, post_yaw = 0;
//        if(solve_post) {
//            tot_pts = pre_x.get_dim();
//            indices all_pts("all_pts");
//            all_pts = range(1, tot_pts);
//            DebugOn("all_pts.size() = " << all_pts.size()<< endl);
//            DebugOn("tot_pts = " << tot_pts << endl);
//            Model<> MPost("MPost");
//            var<> new_x("new_x"), new_y("new_y"), new_z("new_z");
//            var<> zmax("zmax"), zmin("zmin");
//
//            var<> yaw("yaw", -0.1, 0.1), pitch("pitch", -0.1, 0.1), roll("roll", -0.1, 0.1);
//            MPost.add(yaw.in(R(1)),pitch.in(R(1)),roll.in(R(1)));
//            MPost.add(new_x.in(all_pts), new_y.in(all_pts), new_z.in(all_pts));
//            MPost.add(zmax.in(R(1)), zmin.in(R(1)));
//
//            auto z_ids = zmax.repeat_id(tot_pts);
//
//    //        Constraint<> z_max("z_max");
//    //        z_max += zmax.in(z_ids) - new_z;
//    //        MPost.add(z_max.in(all_pts)>=0);
//    //
//    //        Constraint<> z_min("z_min");
//    //        z_min += zmin.in(z_ids) - new_z;
//    //        MPost.add(z_min.in(all_pts)<=0);
//
//            auto ids = yaw.repeat_id(tot_pts);
//
//            /* alpha = yaw_, beta = pitch_ and gamma = roll_ */
//            Constraint<> x_rot("x_rot");
//            x_rot += new_x - x_uav;
//            x_rot -= (pre_x - x_uav)*cos(yaw.in(ids))*cos(pitch.in(ids)) + (pre_y)*(cos(yaw.in(ids))*sin(pitch.in(ids))*sin(roll.in(ids)) - sin(yaw.in(ids))*cos(roll.in(ids))) + (pre_z - z_uav)*(cos(yaw.in(ids))*sin(pitch.in(ids))*cos(roll.in(ids)) + sin(yaw.in(ids))*sin(roll.in(ids)));
//            MPost.add(x_rot.in(all_pts)==0);
//
//            Constraint<> y_rot("y_rot");
//            y_rot += new_y - y_uav;
//            y_rot -= (pre_x - x_uav)*sin(yaw.in(ids))*cos(pitch.in(ids)) + (pre_y - y_uav)*(sin(yaw.in(ids))*sin(pitch.in(ids))*sin(roll.in(ids)) + cos(yaw.in(ids))*cos(roll.in(ids))) + (pre_z - z_uav)*(sin(yaw.in(ids))*sin(pitch.in(ids))*cos(roll.in(ids)) - cos(yaw.in(ids))*sin(roll.in(ids)));
//            MPost.add(y_rot.in(all_pts)==0);
//
//            Constraint<> z_rot("z_rot");
//            z_rot += new_z - z_uav;
//            z_rot -= (pre_x - x_uav)*sin(-1*pitch.in(ids)) + (pre_y - y_uav)*(cos(pitch.in(ids))*sin(roll.in(ids))) + (pre_z - z_uav)*(cos(pitch.in(ids))*cos(roll.in(ids)));
//            MPost.add(z_rot.in(all_pts)==0);
//
//    //        MPost.min(zmax - zmin);
//            MPost.min(norm2(new_z - 1255));
//    //        MPost.print();
//
//            solver<> S(MPost,ipopt);
//            S.run();
//            post_roll = roll.eval();
//            post_pitch = pitch.eval();
//            post_yaw = yaw.eval();
//
//        }
//        if(post_roll!=0 || post_pitch!=0 || post_yaw!=0){
//            gamma = post_roll;
//            beta = post_pitch;
//            alpha = post_yaw;
//            for (auto i = 0; i< x_combined.size(); i++) {
//                shifted_x = x_combined[i] - x_shift[i];
//                shifted_y = y_combined[i] - y_shift[i];
//                shifted_z = zmax_combined[i] - z_shift[i];
//                x_combined[i] = shifted_x*cos(alpha)*cos(beta) + shifted_y*(cos(alpha)*sin(beta)*sin(gamma) - sin(alpha)*cos(gamma)) + shifted_z*(cos(alpha)*sin(beta)*cos(gamma) + sin(alpha)*sin(gamma));
//                y_combined[i] = shifted_x*sin(alpha)*cos(beta) + shifted_y*(sin(alpha)*sin(beta)*sin(gamma) + cos(alpha)*cos(gamma)) + shifted_z*(sin(alpha)*sin(beta)*cos(gamma) - cos(alpha)*sin(gamma));
//                zmax_combined[i] = shifted_x*(-sin(beta)) + shifted_y*(cos(beta)*sin(gamma)) + shifted_z*(cos(beta)*cos(gamma));
//                x_combined[i] += x_shift[i];
//                y_combined[i] += y_shift[i];
//                zmax_combined[i] += z_shift[i];
//            }
//        }
//        DebugOn("Post-Pitch = " << post_pitch*180/pi << endl);
//        DebugOn("Post-Roll = " << post_roll*180/pi << endl);
//        DebugOn("Post-Yaw = " << post_yaw*180/pi << endl);
//
//
//
//        bool plot_data = true;
//        if(plot_data){
//            namespace plt = matplotlibcpp;
//
//            std::map<std::string, std::string> keywords, keywords2;
//            keywords["marker"] = "s";
//            keywords["linestyle"] = "None";
//            keywords["ms"] = "0.05";
//        //    plt::plot3(x_combined, y_combined, zmax_combined, keywords);
//            plt::plot3(uav_x1, uav_y1, uav_z1, x_combined, y_combined, zmax_combined, keywords);
//        //        plt::plot3(uav_x, uav_y, uav_z, x_combined, y_combined, zmax_combined, keywords);
//            keywords2["marker"] = "s";
//            keywords2["ms"] = "0.1";
//            //    plt::plot3(uav_x1, uav_y1, uav_z1, uav_x, uav_y, uav_z, keywords2);
//            //    plt::plot3(uav_x1, uav_y1, uav_z1, keywords);
//            //    plt::plot3(x_vec2, y_vec2, zmax_vec2, keywords);
//            //    plt::colorbar();
//            // Enable legend.
//            //    plt::legend();
//            plt::show();
//        }
//
//        bool save_file = false;
//        if(save_file){
//            DebugOn("Saving new las file");
//        //    LASreadOpener lasreadopener_final;
//        //    lasreadopener_final.set_file_name(LiDAR_file1.c_str());
//        //    lasreadopener_final.set_populate_header(TRUE);
//        //    LASreader* lasreader = lasreadopener_final.open();
//            LASheader lasheader;
//            lasheader.global_encoding = 1;
//            lasheader.x_scale_factor = 0.01;
//            lasheader.y_scale_factor = 0.01;
//            lasheader.z_scale_factor = 0.01;
//            lasheader.x_offset =  500000.0;
//            lasheader.y_offset = 4100000.0;
//            lasheader.z_offset = 0.0;
//            lasheader.point_data_format = 1;
//            lasheader.point_data_record_length = 28;
//
//            LASwriteOpener laswriteopener;
//            laswriteopener.set_file_name("corrected.laz");
//            LASwriter* laswriter = laswriteopener.open(&lasheader);
//            LASpoint laspoint;
//            laspoint.init(&lasheader, lasheader.point_data_format, lasheader.point_data_record_length, 0);
//            for (auto i = 0; i< x_combined.size(); i++) {
//                laspoint.set_x(x_combined[i]);
//                laspoint.set_y(y_combined[i]);
//                laspoint.set_z(zmax_combined[i]);
//                laswriter->write_point(&laspoint);
//                laswriter->update_inventory(&laspoint);
//            }
//            laswriter->update_header(&lasheader, TRUE);
//            laswriter->close();
//            delete laswriter;
//        }
//    }

    
    return 0;
}
