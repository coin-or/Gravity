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

#include "lasreader.hpp"
#include "laswriter.hpp"

using namespace std;

int main (int argc, char * argv[])
{
    int output = 0;
    double tol = 1e-6;
    double solver_time_end, total_time_end, solve_time, total_time;
    double grid_step = 1;
    map<pair<int,int>,tuple<int,int,int,int>> grid1, grid2; /* Grid with all cells where <x,y> is the key and <nb_measurements, min_z, max_z, av_z> is the cell data */
    int ground_z = 0;
    string log_level="0";
    string LiDAR_file1 = string(prj_dir)+"/data_sets/LiDAR/test1.las";
    string LiDAR_file2 = string(prj_dir)+"/data_sets/LiDAR/test2.las";
#ifdef USE_OPT_PARSER
    /** create a OptionParser with options */
    op::OptionParser opt;
    opt.add_option("h", "help", "shows option help"); // no default value means boolean options, which default value is false
    opt.add_option("f", "file", "Input file name (def. ../data_sets/LiDAR/*)", fname_phi );
    opt.add_option("l", "log", "Log level (def. 0)", log_level );
    
    /** parse the options and verify that all went well. If not, errors and help will be shown */
    bool correct_parsing = opt.parse_options(argc, argv);
    
    if(!correct_parsing){
        return EXIT_FAILURE;
    }
    
    LiDAR_file = opt["f"];
    output = op::str2int(opt["l"]);
    output = 5;
    bool has_help = op::str2bool(opt["h"]);
    /** show help */
    if(has_help) {
        opt.show_help();
        exit(0);
    }
#else
    if(argc>1){
        LiDAR_file1 = argv[1];
    }
    if(argc>2){
        LiDAR_file2 = argv[2];
    }
    if(argc>3){
        char *p;
        grid_step =  strtol(argv[3], &p, 10);
    }
#endif
    
//    LiDAR_file = "/Users/l297598/Downloads/PostPost_DAG4_lidar_LDRD_Reserve_20200326_LAS_and_ANPP_files/Roll_Pitch_Yaw_automation_2020_03_26/3.ins_postpost_SGZ_DAG4_L_2_2019_07_31_18_combined_LDRD_Reserve_20200326_Frames_2130_2222_RPY_1.375_0.9_-0.55_NoKin_adc.las";
//   LiDAR_file = "/Users/l297598/Downloads/PostPost_DAG4_lidar_LDRD_Reserve_20200326_LAS_and_ANPP_files/Roll_Pitch_Yaw_automation_2020_03_26/1.raw_postpost_SGZ_DAG4_L_2_2019_07_31_18_combined_LDRD_Reserve_20200326_Frames_2130_2222_RPY_0_0_0_NoKin_adc.las";
//    LiDAR_file = "/Users/l297598/Downloads/PostPost_DAG4_lidar_LDRD_Reserve_20200326_LAS_and_ANPP_files/Roll_Pitch_Yaw_automation_2020_03_26/6.ins_postpost_SGZ_DAG4_L_2_2019_07_31_18_combined_LDRD_Reserve_20200326_Frames_2130_2222_and_2511_2587_RPY_1.375_0.9_-0.55_NoKin_adc.las";
//   LiDAR_file = "/Users/l297598/Downloads/PostPost_DAG4_lidar_LDRD_Reserve_20200326_LAS_and_ANPP_files/Roll_Pitch_Yaw_automation_2020_03_26/5.raw_postpost_SGZ_DAG4_L_2_2019_07_31_18_combined_LDRD_Reserve_20200326_Frames_2130_2222_and_2511_2587_RPY_0_0_0_NoKin_adc.las";
//    LiDAR_file = "/Users/l297598/Downloads/PostPost_DAG4_lidar_LDRD_Reserve_20200326_LAS_and_ANPP_files/Roll_Pitch_Yaw_automation_2020_03_26/7.raw_postpost_SGZ_DAG4_L_2_2019_07_31_18_combined_LDRD_Reserve_20200326_Frames_2130_2587_RPY_0_0_0_NoKin_adc.las";
//    LiDAR_file = "/Users/l297598/Downloads/PostPost_DAG4_lidar_LDRD_Reserve_20200326_LAS_and_ANPP_files/Roll_Pitch_Yaw_automation_2020_03_26/8.ins_postpost_SGZ_DAG4_L_2_2019_07_31_18_combined_LDRD_Reserve_20200326_Frames_2130_2587_RPY_1.375_0.9_-0.55_NoKin_adc.las";
//    LiDAR_file = "/Users/l297598/Downloads/PostPost_DAG4_lidar_LDRD_Reserve_20200326_LAS_and_ANPP_files/PostPost_DAG4_lidar_LDRD_Reserve_20200326_LAS_files/raw_postpost_1g11d_DAG4_L_2_2019_07_31_18_combined_LDRD_Reserve_20200326_Frames_2892_3073_RPY_0_0_0_NoKin_adc.las";
//    LiDAR_file = "/Users/l297598/Downloads/PostPost_DAG4_lidar_LDRD_Reserve_20200326_LAS_and_ANPP_files/PostPost_DAG4_lidar_LDRD_Reserve_20200326_LAS_files/ins_postpost_1g11d_DAG4_L_2_2019_07_31_18_combined_LDRD_Reserve_20200326_Frames_2892_3073_RPY_1.375_0.9_-0.55_NoKin_adc.las";
    LASreadOpener lasreadopener;
    lasreadopener.set_file_name(LiDAR_file1.c_str());
    lasreadopener.set_populate_header(TRUE);
    int xdim1=0, ydim1=0, zdim1=0;
    if (!lasreadopener.active())
    {
        throw invalid_argument("ERROR: no input specified\n");
    }
    vector<double> x_vec1,y_vec1,z_vec1,zmin_vec1,zmax_vec1;
    vector<double> x_combined,y_combined,z_combined,zmin_combined,zmax_combined;
    
    while (lasreadopener.active())
    {
        LASreader* lasreader = lasreadopener.open();
        if (lasreader == 0)
        {
            throw invalid_argument("ERROR: could not open lasreader\n");
        }
        
        xdim1 = ceil(lasreader->header.max_x*grid_step - lasreader->header.min_x*grid_step);
        ydim1 = ceil(lasreader->header.max_y*grid_step - lasreader->header.min_y*grid_step);
        zdim1 = ceil(lasreader->header.max_z*grid_step - lasreader->header.min_z*grid_step);
        ground_z = floor(lasreader->header.min_z);
        DebugOn("Number of points = " << lasreader->npoints << endl);
        DebugOn("min x axis = " << lasreader->header.min_x << endl);
        DebugOn("max x axis = " << lasreader->header.max_x << endl);
        DebugOn("min y axis = " << lasreader->header.min_y << endl);
        DebugOn("max y axis = " << lasreader->header.max_y << endl);
        
        DebugOn("dimension in x axis = " << xdim1 << endl);
        DebugOn("dimension in y axis = " << ydim1 << endl);
        DebugOn("dimension in z axis = " << zdim1 << endl);
        
        
        int nb_dots; /* Number of measurements inside cell */
        int xpos, ypos, z, min_z, max_z, av_z;
        pair<int,int> pos;
        size_t nb_pts = 0;
        tuple<int,int,int,int> cell; /* <min_z,max_z,av_z> */
        while (lasreader->read_point())
        {
//            if(nb_pts==0){
//                DebugOn("(" << to_string_with_precision(lasreader->point.get_x(),10) <<"," << to_string_with_precision(lasreader->point.get_y(),10) << ","<< to_string_with_precision(lasreader->point.get_z(),10) <<")"<<endl);
//                return 0;
//            }
            xpos = floor((lasreader->point.get_x()*grid_step - lasreader->header.min_x*grid_step));
            ypos = floor((lasreader->point.get_y()*grid_step - lasreader->header.min_y*grid_step));
            DebugOff("(" << drone_x <<"," << drone_y << ","<<drone_z<<")"<<endl);
            if (ypos < 0)
            {
                fprintf(stderr, "ERROR: ypos = %d\n", ypos);
            }
            
            pos = make_pair(xpos,ypos);
            z = lasreader->point.get_Z();
            cell = make_tuple(1,z,z,z);
            DebugOff("new z value = " << z << endl);
            DebugOff("fractional z value = " << lasreader->point.get_z() << endl);
            DebugOff("z range = " << grid_maxz[pos] - grid_minz[pos] << endl << endl);
            auto map_pair = grid1.insert(make_pair(pos,cell));
            if (!map_pair.second) {/* Cell already exists */
//                DebugOn("Found overlapping cell after " << nb_pts << endl);
//                DebugOn("(" << to_string_with_precision(lasreader->point.get_x(),10) <<"," << to_string_with_precision(lasreader->point.get_y(),10) << ","<< to_string_with_precision(lasreader->point.get_z(),10) <<")"<<endl);
//                                return 0;
                nb_dots = get<0>(map_pair.first->second);
                min_z = get<1>(map_pair.first->second);
                max_z = get<2>(map_pair.first->second);
                av_z = get<3>(map_pair.first->second);
                /* Update cell data */
                get<0>(map_pair.first->second)++;
                if(min_z>z){
                    DebugOff("updating  min z at cell [" << xpos << "," << ypos << "]" << endl);
                    DebugOff("previous min z = " << min_z << endl);
                    DebugOff("new min z = " << z << endl);
                    get<1>(map_pair.first->second) = z;
                }
                if(max_z<z){
                    DebugOff("updating  max z at cell [" << xpos << "," << ypos << "]" << endl);
                    DebugOff("previous max z = " << max_z << endl);
                    DebugOff("new max z = " << z << endl);
                    get<2>(map_pair.first->second) = z;
                }
                av_z *= nb_dots;/* remove previous denominator */
                get<3>(map_pair.first->second) = (av_z + z)/(nb_dots+1); /* Update denominator */
            }
            nb_pts++;
        }
        DebugOn("Read " << nb_pts << " points" << endl);
        DebugOn("2D grid has " << grid1.size() << " cells " << endl);
        DebugOn("xdim = " << xdim1 << endl);
        DebugOn("ydim = " << ydim1 << endl);
        int max_z_range = 0, cell_counter = 0, z_range = 0,  max_nb_dots = 0;
        double av_nb_dots = 0, av_z_range = 0;
        for (auto &cell : grid1) {
            DebugOff("Cell num " << cell_counter++ << " at [" << cell.first.first << "," << cell.first.second << "]" << endl);
            nb_dots = get<0>(cell.second);
            min_z = get<1>(cell.second);
            max_z = get<2>(cell.second);
            av_z = get<3>(cell.second);
            x_vec1.push_back(cell.first.first);
            y_vec1.push_back(cell.first.second);
            z_vec1.push_back(av_z);
            zmin_vec1.push_back(min_z);
            zmax_vec1.push_back(max_z);
            if(max_nb_dots<nb_dots)
                max_nb_dots = nb_dots;
            av_nb_dots += nb_dots;
            z_range = max_z - min_z;
            DebugOff("z range = " << z_range << endl);
            av_z_range += z_range;
            if(z_range > max_z_range)
                max_z_range = z_range;
            DebugOff("z min = " << min_z << endl);
            DebugOff("z max = " << max_z << endl);
            DebugOff("z av = " << av_z << endl);
            DebugOff("nb points = " << nb_dots << endl << endl);
        }
        av_nb_dots /= grid1.size();
        av_z_range /= grid1.size();
        DebugOn("Max z range for a cell = " << max_z_range << endl);
        DebugOn("Average z range for a cell = " << av_z_range << endl);
        DebugOn("Average points per cell = " << av_nb_dots << endl);
        DebugOn("Max points per cell = " << max_nb_dots << endl << endl);
//        LASwriteOpener laswriteopener;
//        laswriteopener.set_file_name("compressed.laz");
//        LASwriter* laswriter = laswriteopener.open(&lasreader->header);
//
//        while (lasreader->read_point()) laswriter->write_point(&lasreader->point);
//
//        laswriter->close();
//        delete laswriter;
        
        lasreader->close();
        delete lasreader;
    }
    
    int xdim2=0, ydim2=0, zdim2=0;
    lasreadopener.set_file_name(LiDAR_file2.c_str());
    lasreadopener.set_populate_header(TRUE);
    if (!lasreadopener.active())
    {
        throw invalid_argument("ERROR: no input specified\n");
    }
    vector<double> x_vec2,y_vec2,z_vec2,zmin_vec2,zmax_vec2;
    
    while (lasreadopener.active())
    {
        LASreader* lasreader = lasreadopener.open();
        if (lasreader == 0)
        {
            throw invalid_argument("ERROR: could not open lasreader\n");
        }
        
        xdim2 = ceil(lasreader->header.max_x*grid_step - lasreader->header.min_x*grid_step);
        ydim2 = ceil(lasreader->header.max_y*grid_step - lasreader->header.min_y*grid_step);
        zdim2 = ceil(lasreader->header.max_z*grid_step - lasreader->header.min_z*grid_step);
        ground_z = std::min(ground_z,(int)floor(lasreader->header.min_z));
        DebugOn("Number of points = " << lasreader->npoints << endl);
        DebugOn("min x axis = " << lasreader->header.min_x << endl);
        DebugOn("max x axis = " << lasreader->header.max_x << endl);
        DebugOn("min y axis = " << lasreader->header.min_y << endl);
        DebugOn("max y axis = " << lasreader->header.max_y << endl);
        
        DebugOn("dimension in x axis = " << xdim2 << endl);
        DebugOn("dimension in y axis = " << ydim2 << endl);
        DebugOn("dimension in z axis = " << zdim2 << endl);
        
        
        int nb_dots; /* Number of measurements inside cell */
        int xpos, ypos, z, min_z, max_z, av_z;
        pair<int,int> pos;
        size_t nb_pts = 0;
        tuple<int,int,int,int> cell; /* <min_z,max_z,av_z> */
        while (lasreader->read_point())
        {
            //            if(nb_pts==0){
            //                DebugOn("(" << to_string_with_precision(lasreader->point.get_x(),10) <<"," << to_string_with_precision(lasreader->point.get_y(),10) << ","<< to_string_with_precision(lasreader->point.get_z(),10) <<")"<<endl);
            //                return 0;
            //            }
            xpos = floor((lasreader->point.get_x()*grid_step - lasreader->header.min_x*grid_step));
            ypos = floor((lasreader->point.get_y()*grid_step - lasreader->header.min_y*grid_step));
            DebugOff("(" << drone_x <<"," << drone_y << ","<<drone_z<<")"<<endl);
            if (ypos < 0)
            {
                fprintf(stderr, "ERROR: ypos = %d\n", ypos);
            }
            
            pos = make_pair(xpos,ypos);
            z = lasreader->point.get_Z();
            cell = make_tuple(1,z,z,z);
            DebugOff("new z value = " << z << endl);
            DebugOff("fractional z value = " << lasreader->point.get_z() << endl);
            DebugOff("z range = " << grid_maxz[pos] - grid_minz[pos] << endl << endl);
            auto map_pair = grid2.insert(make_pair(pos,cell));
            if (!map_pair.second) {/* Cell already exists */
                //                DebugOn("Found overlapping cell after " << nb_pts << endl);
                //                DebugOn("(" << to_string_with_precision(lasreader->point.get_x(),10) <<"," << to_string_with_precision(lasreader->point.get_y(),10) << ","<< to_string_with_precision(lasreader->point.get_z(),10) <<")"<<endl);
                //                                return 0;
                nb_dots = get<0>(map_pair.first->second);
                min_z = get<1>(map_pair.first->second);
                max_z = get<2>(map_pair.first->second);
                av_z = get<3>(map_pair.first->second);
                /* Update cell data */
                get<0>(map_pair.first->second)++;
                if(min_z>z){
                    DebugOff("updating  min z at cell [" << xpos << "," << ypos << "]" << endl);
                    DebugOff("previous min z = " << min_z << endl);
                    DebugOff("new min z = " << z << endl);
                    get<1>(map_pair.first->second) = z;
                }
                if(max_z<z){
                    DebugOff("updating  max z at cell [" << xpos << "," << ypos << "]" << endl);
                    DebugOff("previous max z = " << max_z << endl);
                    DebugOff("new max z = " << z << endl);
                    get<2>(map_pair.first->second) = z;
                }
                av_z *= nb_dots;/* remove previous denominator */
                get<3>(map_pair.first->second) = (av_z + z)/(nb_dots+1); /* Update denominator */
            }
            nb_pts++;
        }
        DebugOn("Read " << nb_pts << " points" << endl);
        DebugOn("2D grid has " << grid2.size() << " cells " << endl);
        DebugOn("xdim = " << xdim2 << endl);
        DebugOn("ydim = " << ydim2 << endl);
        int max_z_range = 0, cell_counter = 0, z_range = 0,  max_nb_dots = 0;
        double av_nb_dots = 0, av_z_range = 0;
        for (auto &cell : grid2) {
            DebugOff("Cell num " << cell_counter++ << " at [" << cell.first.first << "," << cell.first.second << "]" << endl);
            nb_dots = get<0>(cell.second);
            min_z = get<1>(cell.second);
            max_z = get<2>(cell.second);
            av_z = get<3>(cell.second);
            x_vec2.push_back(cell.first.first);
            y_vec2.push_back(cell.first.second);
            z_vec2.push_back(av_z);
            zmin_vec2.push_back(min_z);
            zmax_vec2.push_back(max_z);
            if(max_nb_dots<nb_dots)
                max_nb_dots = nb_dots;
            av_nb_dots += nb_dots;
            z_range = max_z - min_z;
            DebugOff("z range = " << z_range << endl);
            av_z_range += z_range;
            if(z_range > max_z_range)
                max_z_range = z_range;
            DebugOff("z min = " << min_z << endl);
            DebugOff("z max = " << max_z << endl);
            DebugOff("z av = " << av_z << endl);
            DebugOff("nb points = " << nb_dots << endl << endl);
        }
        av_nb_dots /= grid2.size();
        av_z_range /= grid2.size();
        DebugOn("Max z range for a cell = " << max_z_range << endl);
        DebugOn("Average z range for a cell = " << av_z_range << endl);
        DebugOn("Average points per cell = " << av_nb_dots << endl);
        DebugOn("Max points per cell = " << max_nb_dots << endl << endl);
        //        LASwriteOpener laswriteopener;
        //        laswriteopener.set_file_name("compressed.laz");
        //        LASwriter* laswriter = laswriteopener.open(&lasreader->header);
        //
        //        while (lasreader->read_point()) laswriter->write_point(&lasreader->point);
        //
        //        laswriter->close();
        //        delete laswriter;
        
        lasreader->close();
        delete lasreader;
    }
    double alpha = 0, beta = 0, gamma = 0; /* alpha = yaw, beta = pitch and gamma = roll */
//    alpha = -0.55;
    beta = 68*pi/180.;
//    gamma = 60*pi/180.;
    x_combined.resize(x_vec1.size()+x_vec2.size());
    y_combined.resize(x_combined.size());
    z_combined.resize(x_combined.size());
    zmin_combined.resize(x_combined.size());
    zmax_combined.resize(x_combined.size());
    for (auto i = 0; i< x_vec1.size(); i++) {
//        x_combined[i] = x_vec1[i]*cos(alpha)*cos(beta) + y_vec1[i]*(cos(alpha)*sin(beta)*sin(gamma) - sin(alpha)*cos(gamma)) + z_vec1[i]*(cos(alpha)*sin(beta)*cos(gamma) + sin(alpha)*sin(gamma));
//        y_combined[i] = x_vec1[i]*sin(alpha)*cos(beta) + y_vec1[i]*(sin(alpha)*sin(beta)*sin(gamma) + cos(alpha)*cos(gamma)) + z_vec1[i]*(sin(alpha)*sin(beta)*cos(gamma) - cos(alpha)*sin(gamma));
        x_combined[i] = x_vec1[i];
        y_combined[i] = y_vec1[i];
        z_combined[i] = x_vec1[i]*(-sin(beta)) + y_vec1[i]*(cos(beta)*sin(gamma)) + z_vec1[i]*(cos(beta)*cos(gamma));
        zmin_combined[i] = x_vec1[i]*(-sin(beta)) + y_vec1[i]*(cos(beta)*sin(gamma)) + zmin_vec1[i]*(cos(beta)*cos(gamma));
        zmax_combined[i] = x_vec1[i]*(-sin(beta)) + y_vec1[i]*(cos(beta)*sin(gamma)) + zmax_vec1[i]*(cos(beta)*cos(gamma));
    }
    /* Moving to second flight line assumed to be in opposite direction, pitch (beta) is not affected by drone direction */
//    alpha *= -1;
//    gamma *= -1;
    beta *= -1;
    for (auto i = 0; i< x_vec2.size(); i++) {
//        x_combined[x_vec1.size()+i] = x_vec2[i]*cos(alpha)*cos(beta) + y_vec2[i]*(cos(alpha)*sin(beta)*sin(gamma) - sin(alpha)*cos(gamma)) + z_vec2[i]*(cos(alpha)*sin(beta)*cos(gamma) + sin(alpha)*sin(gamma));
//        y_combined[x_vec1.size()+i] = x_vec2[i]*sin(alpha)*cos(beta) + y_vec2[i]*(sin(alpha)*sin(beta)*sin(gamma) + cos(alpha)*cos(gamma)) + z_vec2[i]*(sin(alpha)*sin(beta)*cos(gamma) - cos(alpha)*sin(gamma));
        x_combined[x_vec1.size()+i] = x_vec2[i];
        y_combined[x_vec1.size()+i] = y_vec2[i];
        z_combined[x_vec1.size()+i] = x_vec2[i]*(-sin(beta)) + y_vec2[i]*(cos(beta)*sin(gamma)) + z_vec2[i]*(cos(beta)*cos(gamma));
        zmin_combined[x_vec1.size()+i] = x_vec2[i]*(-sin(beta)) + y_vec2[i]*(cos(beta)*sin(gamma)) + zmin_vec2[i]*(cos(beta)*cos(gamma));
        zmax_combined[x_vec1.size()+i] = x_vec2[i]*(-sin(beta)) + y_vec2[i]*(cos(beta)*sin(gamma)) + zmax_vec2[i]*(cos(beta)*cos(gamma))-180;
    }
    namespace plt = matplotlibcpp;
//    std::vector<std::vector<double>> x, y, z;
//    for (int i = 0; i <= xdim;  i++) {
//        std::vector<double> x_row, y_row, z_row;
//        for (int j = 0; j <= ydim; j++) {
//            auto it = grid.find(make_pair(i, j));
//            if(it!=grid.end()){
//                x_row.push_back(i);
//                y_row.push_back(j);
//                z_row.push_back(get<3>(it->second));
//            }
//            else{
//                x_row.push_back(i);
//                y_row.push_back(j);
//                z_row.push_back(min_z);
//            }
//        }
//        x.push_back(x_row);
//        y.push_back(y_row);
//        z.push_back(z_row);
//    }
    std::map<std::string, std::string> keywords;
//    keywords["marker"] = ",";
    keywords["marker"] = "s";
    keywords["linestyle"] = "None";
    keywords["ms"] = "0.05";
//    keywords["markerfacecolor"] = "(1, 1, 0, 0.5)";
//    keywords["color"] = "orange";
//    plt::scatter(x_vec, y_vec, z_vec);
//    plt::scatter(x_vec, y_vec, zmin_vec);
//    plt::scatter(x_vec, y_vec, zmax_vec);
//    plt::scatter(x_vec, y_vec, zmax_vec, keywords);
    plt::plot3(x_combined, y_combined, zmax_combined, keywords);
//    plt::plot3(x_vec2, y_vec2, zmax_vec2, keywords);
//    plt::colorbar();
    // Enable legend.
//    plt::legend();
    plt::show();
    //using namespace gravity;
    return 0;
}
