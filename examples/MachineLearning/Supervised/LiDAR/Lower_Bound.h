//
//  Lower_Bound.h
//  Gravity
//
//  Created by Smitha on 8/20/21.
//

#ifndef Lower_Bound_h
#define Lower_Bound_h


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
#include <stdio.h>
#include "Lidar_utils.h"
#include "IPH.h"
#ifdef USE_GJK
extern "C" {
#include "openGJK.h"
}
#endif
#ifdef USE_VORO
#include "voro++.hh"
using namespace voro;
#endif
shared_ptr<Model<double>> Align_L2_model(const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data, const vector<vector<double>>& uav_model, const vector<vector<double>>& uav_data,  double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, indices& cells, param<double> dist_cost)
{
    double angle_max = 0.1;
    int nb_pairs = 0, min_nb_pairs = numeric_limits<int>::max(), max_nb_pairs = 0, av_nb_pairs = 0;

    param<> x1("x1"), x2("x2"), y1("y1"), y2("y2"), z1("z1"), z2("z2");
    param<> x_uav1("x_uav1"), y_uav1("y_uav1"), z_uav1("z_uav1");
    param<> x_uav2("x_uav2"), y_uav2("y_uav2"), z_uav2("z_uav2");
    string j_str, i_str;
    double xm_max = numeric_limits<double>::lowest(), ym_max = numeric_limits<double>::lowest(), zm_max = numeric_limits<double>::lowest();
    double xm_min = numeric_limits<double>::max(), ym_min = numeric_limits<double>::max(), zm_min = numeric_limits<double>::max();
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
    for (auto i = 0; i<point_cloud_data.size(); i++) {
        i_str = to_string(i+1);
        x_uav1.add_val(i_str,uav_data.at(i)[0]);
        x1.add_val(i_str,point_cloud_data.at(i)[0]);
        y_uav1.add_val(i_str,uav_data.at(i)[1]);
        y1.add_val(i_str,point_cloud_data.at(i)[1]);
        z_uav1.add_val(i_str,uav_data.at(i)[2]);
        z1.add_val(i_str,point_cloud_data.at(i)[2]);
    }
    

    indices N1("N1"),N2("N2");
    
    int n1 = x1.get_dim();
    int n2 = x2.get_dim();
    DebugOn("n1 = " << n1 << endl);
    DebugOn("n2 = " << n2 << endl);
    
    indices ids = indices("in_x");
    indices idsij = indices("idsij");
       idsij.add_empty_row();
    ids = N2;
    ids.add_empty_row();
    for(auto i=0;i<n1;i++){
        for(auto j=1;j<=n2;j++){
            if(cells.has_key(to_string(i+1)+","+to_string(j))){
                ids.add_in_row(i, to_string(j));
                idsij.add_in_row(i, to_string(i+1)+","+to_string(j));
            }
        }
    }
    N1 = range(1,n1);
    N2 = range(1,n2);
    //auto cells = indices(N1,N2);
    
    vector<int> new_model_pts;
    indices new_model_ids;
    string name="Norm2_Align";
    
    auto Reg=make_shared<Model<>>(name);
    
    Reg->add_param(x1);Reg->add_param(y1);Reg->add_param(z1);
    Reg->add_param(x2);Reg->add_param(y2);Reg->add_param(z2);
    var<int> bin("bin",0,1);
    Reg->add(bin.in(cells));
    DebugOff("Added " << cells.size() << " binary variables" << endl);
    
    var<> yaw("yaw", yaw_min, yaw_max), pitch("pitch", pitch_min, pitch_max), roll("roll", roll_min, roll_max);
    yaw.in(R(1)); pitch.in(R(1));roll.in(R(1));
    func<> r11 = cos(yaw)*cos(roll);r11.eval_all();
    //R11=r11
    func<> r12 = cos(yaw)*sin(roll)*sin(pitch) - sin(yaw)*cos(pitch);r12.eval_all();
    //R12=-r12
    func<> r13 = cos(yaw)*sin(roll)*cos(pitch) + sin(yaw)*sin(pitch);r13.eval_all();
    //R13=-r13
    func<> r21 = sin(yaw)*cos(roll);r21.eval_all();
    //R21=-r21
    func<> r22 = sin(yaw)*sin(roll)*sin(pitch) + cos(yaw)*cos(pitch);r22.eval_all();
    //R22=r22
    func<> r23 = sin(yaw)*sin(roll)*cos(pitch) - cos(yaw)*sin(pitch);r23.eval_all();
    //R23=r23
    func<> r31 = sin(-1*roll);r31.eval_all();
    //R31=-r31
    func<> r32 = cos(roll)*sin(pitch);r32.eval_all();
    //R32=r32
    func<> r33 = cos(roll)*cos(pitch);r33.eval_all();
    //R33=r33
    
    
    var<> theta11("theta11",  std::max(-1.,r11._range->first), std::min(1.,r11._range->second)), theta12("theta12", std::max(-1.,r12._range->first), std::min(1.,r12._range->second)), theta13("theta13", std::max(-1.,r13._range->first), std::min(1.,r13._range->second));
    var<> theta21("theta21", std::max(-1.,r21._range->first), std::min(1.,r21._range->second)), theta22("theta22", std::max(-1.,r22._range->first), std::min(1.,r22._range->second)), theta23("theta23", std::max(-1.,r23._range->first), std::min(1.,r23._range->second));
    var<> theta31("theta31", std::max(-1.,r31._range->first), std::min(1.,r31._range->second)), theta32("theta32", std::max(-1.,r32._range->first), std::min(1.,r32._range->second)), theta33("theta33", std::max(-1.,r33._range->first), std::min(1.,r33._range->second));
    
    
    Reg->add(theta11.in(R(1)),theta12.in(R(1)),theta13.in(R(1)));
    Reg->add(theta21.in(R(1)),theta22.in(R(1)),theta23.in(R(1)));
    Reg->add(theta31.in(R(1)),theta32.in(R(1)),theta33.in(R(1)));

    
    var<> new_xm("new_xm"), new_ym("new_ym"), new_zm("new_zm");
    
    var<> xm("xm"), ym("ym"), zm("zm");
    Reg->add(new_xm.in(N1), new_ym.in(N1), new_zm.in(N1));
    Reg->add(xm.in(N2), ym.in(N2), zm.in(N2));
    
    
  
    
    Constraint<> Def_newxm("Def_newxm");
    Def_newxm = new_xm-product(xm.in(ids),bin.in_matrix(1, 1));
    Reg->add(Def_newxm.in(N1)==0);
    

    Constraint<> Def_newym("Def_newym");
    Def_newym = new_ym-product(ym.in(ids),bin.in_matrix(1, 1));
    Reg->add(Def_newym.in(N1)==0);
    
    Constraint<> Def_newzm("Def_newzm");
    Def_newzm = new_zm-product(zm.in(ids),bin.in_matrix(1, 1));
    Reg->add(Def_newzm.in(N1)==0);
    
    Constraint<> OneBin("OneBin");
    OneBin = bin.in_matrix(1, 1);
    Reg->add(OneBin.in(N1)==1);
    
    
    var<> x_diff("x_diff"), y_diff("y_diff"), z_diff("z_diff");
    Reg->add(x_diff.in(N1), y_diff.in(N1), z_diff.in(N1));
    
    auto ids1 = theta11.repeat_id(N1.size());
    Constraint<> xd_trans("xd_trans");
    xd_trans += x_diff - (((x1 - x_uav1)*theta11.in(ids1) + (y1 - y_uav1)*theta12.in(ids1)*(-1) + (z1 - z_uav1)*theta13.in(ids1)*(-1) + x_uav1) - new_xm);
    Reg->add(xd_trans.in(N1)==0);
    

    Constraint<> yd_trans("yd_trans");
    yd_trans += y_diff - (((x1 - x_uav1)*theta21.in(ids1)*(-1) + (y1 - y_uav1)*theta22.in(ids1) + (z1 - z_uav1)*theta23.in(ids1) + y_uav1) - new_ym);
    Reg->add(yd_trans.in(N1)==0);
    
    
    Constraint<> zd_trans("zd_trans");
    zd_trans += z_diff - (((x1 - x_uav1)*theta31.in(ids1)*(-1) + (y1 - y_uav1)*theta32.in(ids1) + (z1 - z_uav1)*theta33.in(ids1) + z_uav1) - new_zm);
    Reg->add(zd_trans.in(N1)==0);
    
    auto ids2 = theta11.repeat_id(N2.size());
    Constraint<> xm_trans("xm_trans");
    xm_trans += xm - ((x2 - x_uav2)*theta11.in(ids2) + (y2 - y_uav2)*theta12.in(ids2) + (z2 - z_uav2)*theta13.in(ids2) + x_uav2);
    Reg->add(xm_trans.in(N2)==0);
    

    Constraint<> ym_trans("ym_trans");
    ym_trans += ym - ((x2 - x_uav2)*theta21.in(ids2) + (y2 - y_uav2)*theta22.in(ids2) + (z2 - z_uav2)*theta23.in(ids2) + y_uav2);
    Reg->add(ym_trans.in(N2)==0);
    
    
    Constraint<> zm_trans("zm_trans");
    zm_trans += zm - ((x2 - x_uav2)*theta31.in(ids2) + (y2 - y_uav2)*theta32.in(ids2) + (z2 - z_uav2)*theta33.in(ids2) + z_uav2);
    Reg->add(zm_trans.in(N2)==0);

    var<> deltax("deltax"), deltay("deltay"), deltaz("deltaz");
    Reg->add(deltax.in(N1));
    Reg->add(deltay.in(N1));
    Reg->add(deltaz.in(N1));
    
    Constraint<> Def_deltax("Def_deltax");
    Def_deltax=pow(x_diff, 2)-deltax;
    Reg->add(Def_deltax.in(N1)<=0);
    
    Constraint<> Def_deltay("Def_deltay");
    Def_deltay=pow(y_diff, 2)-deltay;
    Reg->add(Def_deltay.in(N1)<=0);
    
    Constraint<> Def_deltaz("Def_deltaz");
    Def_deltaz=pow(z_diff, 2)-deltaz;
    Reg->add(Def_deltaz.in(N1)<=0);
    
    
    if(dist_cost._indices->_keys->size()!=0){
        Constraint<> delta_cost("delta_cost");
        delta_cost=product(dist_cost.in(idsij), bin.in_matrix(1,1))-deltax-deltay-deltaz;
        Reg->add(delta_cost.in(N1)<=0);
    }
    
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
 
  //  Reg->min(sum(pow(x_diff,2) + pow(y_diff,2) + pow(z_diff,2)));
    
    Reg->min(sum(deltax) + sum(deltay)+sum(deltaz));
    //Reg->print();

    
   // solver<> S(Reg,gurobi);
    //S.use_callback();
    //S.run(5,1e-6,30000,1000);
    
    //Reg->print_solution();
//    vector<double> rot(9);
//    vector<int> matching(n1);
//    bool is_rotation = get_solution(Reg, rot, matching);
//    auto pitch_rad = atan2(rot[7], rot[8]);
//    auto roll_rad = atan2(-rot[6], std::sqrt(rot[7]*rot[7]+rot[8]*rot[8]));
//    auto yaw_rad = atan2(rot[3],rot[0]);
//    DebugOn("roll rad "<< roll_rad<<endl);
//    DebugOn("pitch rad "<< pitch_rad<<endl);
//    DebugOn("yaw rad "<< yaw_rad<<endl);
    
   
    return Reg;
}
shared_ptr<Model<double>> Align_L1_model(const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data, const vector<vector<double>>& uav_model, const vector<vector<double>>& uav_data,  double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, indices& cells, param<double> dist_cost)
{
    double angle_max = 0.1;
    int nb_pairs = 0, min_nb_pairs = numeric_limits<int>::max(), max_nb_pairs = 0, av_nb_pairs = 0;

    param<> x1("x1"), x2("x2"), y1("y1"), y2("y2"), z1("z1"), z2("z2");
    param<> x_uav1("x_uav1"), y_uav1("y_uav1"), z_uav1("z_uav1");
    param<> x_uav2("x_uav2"), y_uav2("y_uav2"), z_uav2("z_uav2");
    string j_str, i_str;
    double xm_max = numeric_limits<double>::lowest(), ym_max = numeric_limits<double>::lowest(), zm_max = numeric_limits<double>::lowest();
    double xm_min = numeric_limits<double>::max(), ym_min = numeric_limits<double>::max(), zm_min = numeric_limits<double>::max();
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
    for (auto i = 0; i<point_cloud_data.size(); i++) {
        i_str = to_string(i+1);
        x_uav1.add_val(i_str,uav_data.at(i)[0]);
        x1.add_val(i_str,point_cloud_data.at(i)[0]);
        y_uav1.add_val(i_str,uav_data.at(i)[1]);
        y1.add_val(i_str,point_cloud_data.at(i)[1]);
        z_uav1.add_val(i_str,uav_data.at(i)[2]);
        z1.add_val(i_str,point_cloud_data.at(i)[2]);
    }
    

    indices N1("N1"),N2("N2");
    
    int n1 = x1.get_dim();
    int n2 = x2.get_dim();
    DebugOn("n1 = " << n1 << endl);
    DebugOn("n2 = " << n2 << endl);
    
    indices ids = indices("in_x");
    indices idsij = indices("idsij");
       idsij.add_empty_row();
    ids = N2;
    ids.add_empty_row();
    for(auto i=0;i<n1;i++){
        for(auto j=1;j<=n2;j++){
            if(cells.has_key(to_string(i+1)+","+to_string(j))){
                ids.add_in_row(i, to_string(j));
                idsij.add_in_row(i, to_string(i+1)+","+to_string(j));
            }
        }
    }
    N1 = range(1,n1);
    N2 = range(1,n2);
    //auto cells = indices(N1,N2);
    
    vector<int> new_model_pts;
    indices new_model_ids;
    string name="Norm2_Align";
    
    auto Reg=make_shared<Model<>>(name);
    
    Reg->add_param(x1);Reg->add_param(y1);Reg->add_param(z1);
    Reg->add_param(x2);Reg->add_param(y2);Reg->add_param(z2);
    var<int> bin("bin",0,1);
    Reg->add(bin.in(cells));
    DebugOff("Added " << cells.size() << " binary variables" << endl);
    
    var<> yaw("yaw", yaw_min, yaw_max), pitch("pitch", pitch_min, pitch_max), roll("roll", roll_min, roll_max);
    yaw.in(R(1)); pitch.in(R(1));roll.in(R(1));
    func<> r11 = cos(yaw)*cos(roll);r11.eval_all();
    //R11=r11
    func<> r12 = cos(yaw)*sin(roll)*sin(pitch) - sin(yaw)*cos(pitch);r12.eval_all();
    //R12=-r12
    func<> r13 = cos(yaw)*sin(roll)*cos(pitch) + sin(yaw)*sin(pitch);r13.eval_all();
    //R13=-r13
    func<> r21 = sin(yaw)*cos(roll);r21.eval_all();
    //R21=-r21
    func<> r22 = sin(yaw)*sin(roll)*sin(pitch) + cos(yaw)*cos(pitch);r22.eval_all();
    //R22=r22
    func<> r23 = sin(yaw)*sin(roll)*cos(pitch) - cos(yaw)*sin(pitch);r23.eval_all();
    //R23=r23
    func<> r31 = sin(-1*roll);r31.eval_all();
    //R31=-r31
    func<> r32 = cos(roll)*sin(pitch);r32.eval_all();
    //R32=r32
    func<> r33 = cos(roll)*cos(pitch);r33.eval_all();
    //R33=r33
    
    
    var<> theta11("theta11",  std::max(-1.,r11._range->first), std::min(1.,r11._range->second)), theta12("theta12", std::max(-1.,r12._range->first), std::min(1.,r12._range->second)), theta13("theta13", std::max(-1.,r13._range->first), std::min(1.,r13._range->second));
    var<> theta21("theta21", std::max(-1.,r21._range->first), std::min(1.,r21._range->second)), theta22("theta22", std::max(-1.,r22._range->first), std::min(1.,r22._range->second)), theta23("theta23", std::max(-1.,r23._range->first), std::min(1.,r23._range->second));
    var<> theta31("theta31", std::max(-1.,r31._range->first), std::min(1.,r31._range->second)), theta32("theta32", std::max(-1.,r32._range->first), std::min(1.,r32._range->second)), theta33("theta33", std::max(-1.,r33._range->first), std::min(1.,r33._range->second));
    
    
    Reg->add(theta11.in(R(1)),theta12.in(R(1)),theta13.in(R(1)));
    Reg->add(theta21.in(R(1)),theta22.in(R(1)),theta23.in(R(1)));
    Reg->add(theta31.in(R(1)),theta32.in(R(1)),theta33.in(R(1)));

    
    var<> new_xm("new_xm"), new_ym("new_ym"), new_zm("new_zm");
    
    var<> xm("xm"), ym("ym"), zm("zm");
    Reg->add(new_xm.in(N1), new_ym.in(N1), new_zm.in(N1));
    Reg->add(xm.in(N2), ym.in(N2), zm.in(N2));
    
    
  
    
    Constraint<> Def_newxm("Def_newxm");
    Def_newxm = new_xm-product(xm.in(ids),bin.in_matrix(1, 1));
    Reg->add(Def_newxm.in(N1)==0);
    

    Constraint<> Def_newym("Def_newym");
    Def_newym = new_ym-product(ym.in(ids),bin.in_matrix(1, 1));
    Reg->add(Def_newym.in(N1)==0);
    
    Constraint<> Def_newzm("Def_newzm");
    Def_newzm = new_zm-product(zm.in(ids),bin.in_matrix(1, 1));
    Reg->add(Def_newzm.in(N1)==0);
    
    Constraint<> OneBin("OneBin");
    OneBin = bin.in_matrix(1, 1);
    Reg->add(OneBin.in(N1)==1);
    
    
    var<> x_diff("x_diff"), y_diff("y_diff"), z_diff("z_diff");
    Reg->add(x_diff.in(N1), y_diff.in(N1), z_diff.in(N1));
    
    auto ids1 = theta11.repeat_id(N1.size());
    Constraint<> xd_trans1("xd_trans1");
    xd_trans1 += x_diff - (((x1 - x_uav1)*theta11.in(ids1) + (y1 - y_uav1)*theta12.in(ids1)*(-1) + (z1 - z_uav1)*theta13.in(ids1)*(-1) + x_uav1) - new_xm);
    Reg->add(xd_trans1.in(N1)>=0);
    
    Constraint<> xd_trans2("xd_trans2");
    xd_trans2 += x_diff + (((x1 - x_uav1)*theta11.in(ids1) + (y1 - y_uav1)*theta12.in(ids1)*(-1) + (z1 - z_uav1)*theta13.in(ids1)*(-1) + x_uav1) - new_xm);
    Reg->add(xd_trans2.in(N1)>=0);
    

    Constraint<> yd_trans1("yd_trans1");
    yd_trans1 += y_diff - (((x1 - x_uav1)*theta21.in(ids1)*(-1) + (y1 - y_uav1)*theta22.in(ids1) + (z1 - z_uav1)*theta23.in(ids1) + y_uav1) - new_ym);
    Reg->add(yd_trans1.in(N1)>=0);
    
    Constraint<> yd_trans2("yd_trans2");
    yd_trans2 += y_diff + (((x1 - x_uav1)*theta21.in(ids1)*(-1) + (y1 - y_uav1)*theta22.in(ids1) + (z1 - z_uav1)*theta23.in(ids1) + y_uav1) - new_ym);
    Reg->add(yd_trans2.in(N1)>=0);
    
    
    Constraint<> zd_trans1("zd_trans1");
    zd_trans1 += z_diff - (((x1 - x_uav1)*theta31.in(ids1)*(-1) + (y1 - y_uav1)*theta32.in(ids1) + (z1 - z_uav1)*theta33.in(ids1) + z_uav1) - new_zm);
    Reg->add(zd_trans1.in(N1)>=0);
    
    
    Constraint<> zd_trans2("zd_trans2");
    zd_trans2 += z_diff + (((x1 - x_uav1)*theta31.in(ids1)*(-1) + (y1 - y_uav1)*theta32.in(ids1) + (z1 - z_uav1)*theta33.in(ids1) + z_uav1) - new_zm);
    Reg->add(zd_trans2.in(N1)>=0);
    
    auto ids2 = theta11.repeat_id(N2.size());
    Constraint<> xm_trans("xm_trans");
    xm_trans += xm - ((x2 - x_uav2)*theta11.in(ids2) + (y2 - y_uav2)*theta12.in(ids2) + (z2 - z_uav2)*theta13.in(ids2) + x_uav2);
    Reg->add(xm_trans.in(N2)==0);
    

    Constraint<> ym_trans("ym_trans");
    ym_trans += ym - ((x2 - x_uav2)*theta21.in(ids2) + (y2 - y_uav2)*theta22.in(ids2) + (z2 - z_uav2)*theta23.in(ids2) + y_uav2);
    Reg->add(ym_trans.in(N2)==0);
    
    
    Constraint<> zm_trans("zm_trans");
    zm_trans += zm - ((x2 - x_uav2)*theta31.in(ids2) + (y2 - y_uav2)*theta32.in(ids2) + (z2 - z_uav2)*theta33.in(ids2) + z_uav2);
    Reg->add(zm_trans.in(N2)==0);
    
    if(dist_cost._indices->_keys->size()!=0){
        Constraint<> delta_cost("delta_cost");
        delta_cost=product(dist_cost.in(idsij), bin.in_matrix(1,1))-x_diff-y_diff-z_diff;
        Reg->add(delta_cost.in(N1)<=0);
    }
    
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
 
  //  Reg->min(sum(pow(x_diff,2) + pow(y_diff,2) + pow(z_diff,2)));
    
    Reg->min(sum(x_diff) + sum(y_diff)+sum(z_diff));
//    Reg->print();
//
//
//    solver<> S(Reg,gurobi);
//    S.use_callback();
//    S.run(5,1e-6,30000,1000);
//
//    Reg->print_solution();
//    vector<double> rot(9);
//    vector<int> matching(n1);
//    bool is_rotation = get_solution(Reg, rot, matching);
//    auto pitch_rad = atan2(rot[7], rot[8]);
//    auto roll_rad = atan2(-rot[6], std::sqrt(rot[7]*rot[7]+rot[8]*rot[8]));
//    auto yaw_rad = atan2(rot[3],rot[0]);
//    DebugOn("roll rad "<< roll_rad<<endl);
//    DebugOn("pitch rad "<< pitch_rad<<endl);
//    DebugOn("yaw rad "<< yaw_rad<<endl);
    
   
    return Reg;
}
shared_ptr<Model<double>> Align_L2_model_noise(const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data, const vector<vector<double>>& uav_model, const vector<vector<double>>& uav_data,  double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, indices& cells, param<double> dist_cost)
{
    double angle_max = 0.1;
    int nb_pairs = 0, min_nb_pairs = numeric_limits<int>::max(), max_nb_pairs = 0, av_nb_pairs = 0;

    param<> x1("x1"), x2("x2"), y1("y1"), y2("y2"), z1("z1"), z2("z2");
    param<> x_uav1("x_uav1"), y_uav1("y_uav1"), z_uav1("z_uav1");
    param<> x_uav2("x_uav2"), y_uav2("y_uav2"), z_uav2("z_uav2");
    string j_str, i_str;
    double xm_max = numeric_limits<double>::lowest(), ym_max = numeric_limits<double>::lowest(), zm_max = numeric_limits<double>::lowest();
    double xm_min = numeric_limits<double>::max(), ym_min = numeric_limits<double>::max(), zm_min = numeric_limits<double>::max();
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
    for (auto i = 0; i<point_cloud_data.size(); i++) {
        i_str = to_string(i+1);
        x_uav1.add_val(i_str,uav_data.at(i)[0]);
        x1.add_val(i_str,point_cloud_data.at(i)[0]);
        y_uav1.add_val(i_str,uav_data.at(i)[1]);
        y1.add_val(i_str,point_cloud_data.at(i)[1]);
        z_uav1.add_val(i_str,uav_data.at(i)[2]);
        z1.add_val(i_str,point_cloud_data.at(i)[2]);
    }
    

    indices N1("N1"),N2("N2");
    
    int n1 = x1.get_dim();
    int n2 = x2.get_dim();
    DebugOn("n1 = " << n1 << endl);
    DebugOn("n2 = " << n2 << endl);
    
    indices ids = indices("in_x");
    indices idsij = indices("idsij");
       idsij.add_empty_row();
    ids = N2;
    ids.add_empty_row();
    for(auto i=0;i<n1;i++){
        for(auto j=1;j<=n2;j++){
            if(cells.has_key(to_string(i+1)+","+to_string(j))){
                ids.add_in_row(i, to_string(j));
                idsij.add_in_row(i, to_string(i+1)+","+to_string(j));
            }
        }
    }
    N1 = range(1,n1);
    N2 = range(1,n2);
    //auto cells = indices(N1,N2);
    
    vector<int> new_model_pts;
    indices new_model_ids;
    string name="Norm2_Align";
    
    auto Reg=make_shared<Model<>>(name);
    
    Reg->add_param(x1);Reg->add_param(y1);Reg->add_param(z1);
    Reg->add_param(x2);Reg->add_param(y2);Reg->add_param(z2);
    var<int> bin("bin",0,1);
    Reg->add(bin.in(cells));
    DebugOff("Added " << cells.size() << " binary variables" << endl);
    
    var<> yaw("yaw", yaw_min, yaw_max), pitch("pitch", pitch_min, pitch_max), roll("roll", roll_min, roll_max);
    yaw.in(R(1)); pitch.in(R(1));roll.in(R(1));
    func<> r11 = cos(yaw)*cos(roll);r11.eval_all();
    //R11=r11
    func<> r12 = cos(yaw)*sin(roll)*sin(pitch) - sin(yaw)*cos(pitch);r12.eval_all();
    //R12=-r12
    func<> r13 = cos(yaw)*sin(roll)*cos(pitch) + sin(yaw)*sin(pitch);r13.eval_all();
    //R13=-r13
    func<> r21 = sin(yaw)*cos(roll);r21.eval_all();
    //R21=-r21
    func<> r22 = sin(yaw)*sin(roll)*sin(pitch) + cos(yaw)*cos(pitch);r22.eval_all();
    //R22=r22
    func<> r23 = sin(yaw)*sin(roll)*cos(pitch) - cos(yaw)*sin(pitch);r23.eval_all();
    //R23=r23
    func<> r31 = sin(-1*roll);r31.eval_all();
    //R31=-r31
    func<> r32 = cos(roll)*sin(pitch);r32.eval_all();
    //R32=r32
    func<> r33 = cos(roll)*cos(pitch);r33.eval_all();
    //R33=r33
    
    
    var<> theta11("theta11",  std::max(-1.,r11._range->first), std::min(1.,r11._range->second)), theta12("theta12", std::max(-1.,r12._range->first), std::min(1.,r12._range->second)), theta13("theta13", std::max(-1.,r13._range->first), std::min(1.,r13._range->second));
    var<> theta21("theta21", std::max(-1.,r21._range->first), std::min(1.,r21._range->second)), theta22("theta22", std::max(-1.,r22._range->first), std::min(1.,r22._range->second)), theta23("theta23", std::max(-1.,r23._range->first), std::min(1.,r23._range->second));
    var<> theta31("theta31", std::max(-1.,r31._range->first), std::min(1.,r31._range->second)), theta32("theta32", std::max(-1.,r32._range->first), std::min(1.,r32._range->second)), theta33("theta33", std::max(-1.,r33._range->first), std::min(1.,r33._range->second));
    
    
    Reg->add(theta11.in(R(1)),theta12.in(R(1)),theta13.in(R(1)));
    Reg->add(theta21.in(R(1)),theta22.in(R(1)),theta23.in(R(1)));
    Reg->add(theta31.in(R(1)),theta32.in(R(1)),theta33.in(R(1)));

    
    var<> new_xm("new_xm"), new_ym("new_ym"), new_zm("new_zm");
    
    var<> xm("xm"), ym("ym"), zm("zm");
    Reg->add(new_xm.in(N1), new_ym.in(N1), new_zm.in(N1));
    Reg->add(xm.in(N2), ym.in(N2), zm.in(N2));
    
    
  
    
    Constraint<> Def_newxm("Def_newxm");
    Def_newxm = new_xm-product(xm.in(ids),bin.in_matrix(1, 1));
    Reg->add(Def_newxm.in(N1)==0);
    

    Constraint<> Def_newym("Def_newym");
    Def_newym = new_ym-product(ym.in(ids),bin.in_matrix(1, 1));
    Reg->add(Def_newym.in(N1)==0);
    
    Constraint<> Def_newzm("Def_newzm");
    Def_newzm = new_zm-product(zm.in(ids),bin.in_matrix(1, 1));
    Reg->add(Def_newzm.in(N1)==0);
    
    Constraint<> OneBin("OneBin");
    OneBin = bin.in_matrix(1, 1);
    Reg->add(OneBin.in(N1)==1);
    
    
    var<> x_diff("x_diff"), y_diff("y_diff"), z_diff("z_diff");
    Reg->add(x_diff.in(N1), y_diff.in(N1), z_diff.in(N1));
    
    auto ids1 = theta11.repeat_id(N1.size());
    Constraint<> xd_trans("xd_trans");
    xd_trans += x_diff - (((x1 - x_uav1)*theta11.in(ids1) + (y1 - y_uav1)*theta12.in(ids1)*(-1) + (z1 - z_uav1)*theta13.in(ids1)*(-1) + x_uav1) - new_xm);
    Reg->add(xd_trans.in(N1)==0);
    

    Constraint<> yd_trans("yd_trans");
    yd_trans += y_diff - (((x1 - x_uav1)*theta21.in(ids1)*(-1) + (y1 - y_uav1)*theta22.in(ids1) + (z1 - z_uav1)*theta23.in(ids1) + y_uav1) - new_ym);
    Reg->add(yd_trans.in(N1)==0);
    
    
    Constraint<> zd_trans("zd_trans");
    zd_trans += z_diff - (((x1 - x_uav1)*theta31.in(ids1)*(-1) + (y1 - y_uav1)*theta32.in(ids1) + (z1 - z_uav1)*theta33.in(ids1) + z_uav1) - new_zm);
    Reg->add(zd_trans.in(N1)==0);
    
    auto ids2 = theta11.repeat_id(N2.size());
    Constraint<> xm_trans("xm_trans");
    xm_trans += xm - ((x2 - x_uav2)*theta11.in(ids2) + (y2 - y_uav2)*theta12.in(ids2) + (z2 - z_uav2)*theta13.in(ids2) + x_uav2);
    Reg->add(xm_trans.in(N2)==0);
    

    Constraint<> ym_trans("ym_trans");
    ym_trans += ym - ((x2 - x_uav2)*theta21.in(ids2) + (y2 - y_uav2)*theta22.in(ids2) + (z2 - z_uav2)*theta23.in(ids2) + y_uav2);
    Reg->add(ym_trans.in(N2)==0);
    
    
    Constraint<> zm_trans("zm_trans");
    zm_trans += zm - ((x2 - x_uav2)*theta31.in(ids2) + (y2 - y_uav2)*theta32.in(ids2) + (z2 - z_uav2)*theta33.in(ids2) + z_uav2);
    Reg->add(zm_trans.in(N2)==0);

    var<> deltax("deltax"), deltay("deltay"), deltaz("deltaz");
    Reg->add(deltax.in(N1));
    Reg->add(deltay.in(N1));
    Reg->add(deltaz.in(N1));
    
    Constraint<> Def_deltax("Def_deltax");
    Def_deltax=pow(x_diff, 2)-deltax;
    Reg->add(Def_deltax.in(N1)<=0);
    
    Constraint<> Def_deltay("Def_deltay");
    Def_deltay=pow(y_diff, 2)-deltay;
    Reg->add(Def_deltay.in(N1)<=0);
    
    Constraint<> Def_deltaz("Def_deltaz");
    Def_deltaz=pow(z_diff, 2)-deltaz;
    Reg->add(Def_deltaz.in(N1)<=0);
    
    
    if(dist_cost._indices->_keys->size()!=0){
        Constraint<> delta_cost("delta_cost");
        delta_cost=product(dist_cost.in(idsij), bin.in_matrix(1,1))-deltax-deltay-deltaz;
        Reg->add(delta_cost.in(N1)<=0);
    }
    
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
 
  //  Reg->min(sum(pow(x_diff,2) + pow(y_diff,2) + pow(z_diff,2)));
    
    Reg->min(sum(deltax) + sum(deltay)+sum(deltaz));
    //Reg->print();

    
   // solver<> S(Reg,gurobi);
    //S.use_callback();
    //S.run(5,1e-6,30000,1000);
    
    //Reg->print_solution();
//    vector<double> rot(9);
//    vector<int> matching(n1);
//    bool is_rotation = get_solution(Reg, rot, matching);
//    auto pitch_rad = atan2(rot[7], rot[8]);
//    auto roll_rad = atan2(-rot[6], std::sqrt(rot[7]*rot[7]+rot[8]*rot[8]));
//    auto yaw_rad = atan2(rot[3],rot[0]);
//    DebugOn("roll rad "<< roll_rad<<endl);
//    DebugOn("pitch rad "<< pitch_rad<<endl);
//    DebugOn("yaw rad "<< yaw_rad<<endl);
    
   
    return Reg;
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

#endif /* Lower_Bound_h */
