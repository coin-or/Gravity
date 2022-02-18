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
#include <algorithm>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <future>
#include <thread>
#include <time.h>
using namespace gravity;
using namespace std;
#include <stdio.h>
#include "Lidar_utils.h"

shared_ptr<Model<double>> Align_L2_model_rotation_neworder(const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data, const vector<vector<double>>& uav_model, const vector<vector<double>>& uav_data, const vector<vector<double>>& rollpitchyawins_model, const vector<vector<double>>& rollpitchyawins_data, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, indices& cells, param<double> dist_cost)
{
    double angle_max = 0.1;
    int nb_pairs = 0, min_nb_pairs = numeric_limits<int>::max(), max_nb_pairs = 0, av_nb_pairs = 0;

    param<> x1("x1"), x2("x2"), y1("y1"), y2("y2"), z1("z1"), z2("z2");
    param<> x_uav1("x_uav1"), y_uav1("y_uav1"), z_uav1("z_uav1");
    param<> x_uav2("x_uav2"), y_uav2("y_uav2"), z_uav2("z_uav2");
    param<> x1i("x1i"), x2i("x2i"), y1i("y1i"), y2i("y2i"), z1i("z1i"), z2i("z2i");
    param<> roll_ins1("roll_ins1"),pitch_ins1("pitch_ins1"), yaw_ins1("yaw_ins1");
    param<> roll_ins2("roll_ins2"), pitch_ins2("pitch_ins2"), yaw_ins2("yaw_ins2");
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
        roll_ins2.add_val(j_str, rollpitchyawins_model.at(j)[0]);
        pitch_ins2.add_val(j_str, (-pi+rollpitchyawins_model.at(j)[1])*(-1));
        yaw_ins2.add_val(j_str, (-pi/2+rollpitchyawins_model.at(j)[2])*(-1));
        auto res_m=apply_rotation_new_order(rollpitchyawins_model.at(j)[0], (-pi+rollpitchyawins_model.at(j)[1])*(-1), (-pi/2+rollpitchyawins_model.at(j)[2])*(-1), point_cloud_model.at(j)[0]-uav_model.at(j)[0], point_cloud_model.at(j)[1]-uav_model.at(j)[1], point_cloud_model.at(j)[2]-uav_model.at(j)[2]);
        x2i.add_val(j_str,res_m[0]);
        y2i.add_val(j_str,res_m[1]);
        z2i.add_val(j_str,res_m[2]);
    }
    
    for (auto i = 0; i<point_cloud_data.size(); i++) {
        i_str = to_string(i+1);
        x_uav1.add_val(i_str,uav_data.at(i)[0]);
        x1.add_val(i_str,point_cloud_data.at(i)[0]);
        y_uav1.add_val(i_str,uav_data.at(i)[1]);
        y1.add_val(i_str,point_cloud_data.at(i)[1]);
        z_uav1.add_val(i_str,uav_data.at(i)[2]);
        z1.add_val(i_str,point_cloud_data.at(i)[2]);
        roll_ins1.add_val(i_str, rollpitchyawins_data.at(i)[0]);
        pitch_ins1.add_val(i_str, (-pi+rollpitchyawins_data.at(i)[1])*(-1));
        yaw_ins1.add_val(i_str, (-pi/2+rollpitchyawins_data.at(i)[2])*(-1));
        auto res_d=apply_rotation_new_order(rollpitchyawins_data.at(i)[0], (-pi+rollpitchyawins_data.at(i)[1])*(-1), (-pi/2+rollpitchyawins_data.at(i)[2])*(-1), point_cloud_data.at(i)[0]-uav_data.at(i)[0], point_cloud_data.at(i)[1]-uav_data.at(i)[1], point_cloud_data.at(i)[2]-uav_data.at(i)[2]);
        x1i.add_val(i_str,res_d[0]);
        y1i.add_val(i_str,res_d[1]);
        z1i.add_val(i_str,res_d[2]);
    }
    

    indices N1("N1"),N2("N2");
    
    int n1 = x1.get_dim();
    int n2 = x2.get_dim();
    DebugOn("n1 = " << n1 << endl);
    DebugOn("n2 = " << n2 << endl);
    
    N1 = range(1,n1);
    N2 = range(1,n2);
    
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
    var<> xb_d("xb_d"), yb_d("yb_d"), zb_d("zb_d");
    Reg->add(xb_d.in(N1), yb_d.in(N1), zb_d.in(N1));
    var<> xb_m("xb_m"), yb_m("yb_m"), zb_m("zb_m");
    Reg->add(xb_m.in(N2), yb_m.in(N2), zb_m.in(N2));
    
    auto ids1 = theta11.repeat_id(N1.size());
    
    Constraint<> xd_bore("xd_bore");
    xd_bore=x1i*theta11.in(ids1)+y1i*theta12.in(ids1)+z1i*theta13.in(ids1)-xb_d;
    Reg->add(xd_bore.in(N1)==0);
    
    Constraint<> yd_bore("yd_bore");
    yd_bore=x1i*theta21.in(ids1)+y1i*theta22.in(ids1)+z1i*theta23.in(ids1)-yb_d;
    Reg->add(yd_bore.in(N1)==0);

    
    Constraint<> zd_bore("zd_bore");
    zd_bore=x1i*theta31.in(ids1)+y1i*theta32.in(ids1)+z1i*theta33.in(ids1)-zb_d;
    Reg->add(zd_bore.in(N1)==0);
    
  
   
    
    Constraint<> xd_ins_inv("xd_ins_inv");
    xd_ins_inv=xb_d*cos(roll_ins1)*cos(yaw_ins1)+yb_d*(cos(pitch_ins1)*sin(yaw_ins1) +cos(yaw_ins1)*sin(roll_ins1)*sin(pitch_ins1))+zb_d*(sin(pitch_ins1)*sin(yaw_ins1)-cos(pitch_ins1)*cos(yaw_ins1)*sin(roll_ins1))+x_uav1-new_xm-x_diff;
    Reg->add(xd_ins_inv.in(N1)==0);
    
    Constraint<> yd_ins_inv("yd_ins_inv");
    yd_ins_inv=xb_d*(-1)*cos(roll_ins1)*sin(yaw_ins1)+yb_d*(cos(pitch_ins1)*cos(yaw_ins1) -sin(roll_ins1)*sin(pitch_ins1)*sin(yaw_ins1))+zb_d*(cos(yaw_ins1)*sin(pitch_ins1)+cos(pitch_ins1)*sin(roll_ins1)*sin(yaw_ins1))+y_uav1-new_ym-y_diff;
    Reg->add(yd_ins_inv.in(N1)==0);
    
    Constraint<> zd_ins_inv("zd_ins_inv");
    zd_ins_inv=xb_d*sin(roll_ins1)+yb_d*(-1)*cos(roll_ins1)*sin(pitch_ins1)+zb_d*(cos(roll_ins1)*cos(pitch_ins1))+z_uav1-new_zm-z_diff;
    Reg->add(zd_ins_inv.in(N1)==0);
 
    auto ids2 = theta11.repeat_id(N2.size());
    
    Constraint<> xm_bore("xm_bore");
    xm_bore=x2i*theta11.in(ids2)+y2i*theta12.in(ids2)+z2i*theta13.in(ids2)-xb_m;
    Reg->add(xm_bore.in(N2)==0);
    
    Constraint<> ym_bore("ym_bore");
    ym_bore=x2i*theta21.in(ids2)+y2i*theta22.in(ids2)+z2i*theta23.in(ids2)-yb_m;
    Reg->add(ym_bore.in(N2)==0);

    
    Constraint<> zm_bore("zm_bore");
    zm_bore=x2i*theta31.in(ids2)+y2i*theta32.in(ids2)+z2i*theta33.in(ids2)-zb_m;
    Reg->add(zm_bore.in(N2)==0);
    

    
    Constraint<> xm_ins_inv("xm_ins_inv");
    xm_ins_inv=xb_m*cos(roll_ins2)*cos(yaw_ins2)+yb_m*(cos(pitch_ins2)*sin(yaw_ins2) +cos(yaw_ins2)*sin(roll_ins2)*sin(pitch_ins2))+zb_m*(sin(pitch_ins2)*sin(yaw_ins2)-cos(pitch_ins2)*cos(yaw_ins2)*sin(roll_ins2))+x_uav2-xm;
    Reg->add(xm_ins_inv.in(N2)==0);
    
    
    Constraint<> ym_ins_inv("ym_ins_inv");
    ym_ins_inv=xb_m*(-1)*cos(roll_ins2)*sin(yaw_ins2)+yb_m*(cos(pitch_ins2)*cos(yaw_ins2) -sin(roll_ins2)*sin(pitch_ins2)*sin(yaw_ins2))+zb_m*(cos(yaw_ins2)*sin(pitch_ins2)+cos(pitch_ins2)*sin(roll_ins2)*sin(yaw_ins2))+y_uav2-ym;
    Reg->add(ym_ins_inv.in(N2)==0);
    
    Constraint<> zm_ins_inv("zm_ins_inv");
    zm_ins_inv=xb_m*sin(roll_ins2)+yb_m*(-1)*cos(roll_ins2)*sin(pitch_ins2)+zb_m*(cos(roll_ins2)*cos(pitch_ins2))+z_uav2-zm;
    Reg->add(zm_ins_inv.in(N2)==0);
    


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

   
    return Reg;
}
shared_ptr<Model<double>> Align_L2_model_rotation_trigonometric(const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data, const vector<vector<double>>& uav_model, const vector<vector<double>>& uav_data, const vector<vector<double>>& rollpitchyawins_model, const vector<vector<double>>& rollpitchyawins_data, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, indices& cells, param<double> dist_cost)
{
    double angle_max = 0.1;
    int nb_pairs = 0, min_nb_pairs = numeric_limits<int>::max(), max_nb_pairs = 0, av_nb_pairs = 0;

    param<> x1("x1"), x2("x2"), y1("y1"), y2("y2"), z1("z1"), z2("z2");
    param<> x_uav1("x_uav1"), y_uav1("y_uav1"), z_uav1("z_uav1");
    param<> x_uav2("x_uav2"), y_uav2("y_uav2"), z_uav2("z_uav2");
    param<> x1i("x1i"), x2i("x2i"), y1i("y1i"), y2i("y2i"), z1i("z1i"), z2i("z2i");
    param<> roll_ins1("roll_ins1"),pitch_ins1("pitch_ins1"), yaw_ins1("yaw_ins1");
    param<> roll_ins2("roll_ins2"), pitch_ins2("pitch_ins2"), yaw_ins2("yaw_ins2");
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
        roll_ins2.add_val(j_str, rollpitchyawins_model.at(j)[0]);
        pitch_ins2.add_val(j_str, (-pi+rollpitchyawins_model.at(j)[1])*(-1));
        yaw_ins2.add_val(j_str, (-pi/2+rollpitchyawins_model.at(j)[2])*(-1));
        auto res_m=apply_rotation_new_order(rollpitchyawins_model.at(j)[0], (-pi+rollpitchyawins_model.at(j)[1])*(-1), (-pi/2+rollpitchyawins_model.at(j)[2])*(-1), point_cloud_model.at(j)[0]-uav_model.at(j)[0], point_cloud_model.at(j)[1]-uav_model.at(j)[1], point_cloud_model.at(j)[2]-uav_model.at(j)[2]);
        x2i.add_val(j_str,res_m[0]);
        y2i.add_val(j_str,res_m[1]);
        z2i.add_val(j_str,res_m[2]);
    }
    
    for (auto i = 0; i<point_cloud_data.size(); i++) {
        i_str = to_string(i+1);
        x_uav1.add_val(i_str,uav_data.at(i)[0]);
        x1.add_val(i_str,point_cloud_data.at(i)[0]);
        y_uav1.add_val(i_str,uav_data.at(i)[1]);
        y1.add_val(i_str,point_cloud_data.at(i)[1]);
        z_uav1.add_val(i_str,uav_data.at(i)[2]);
        z1.add_val(i_str,point_cloud_data.at(i)[2]);
        roll_ins1.add_val(i_str, rollpitchyawins_data.at(i)[0]);
        pitch_ins1.add_val(i_str, (-pi+rollpitchyawins_data.at(i)[1])*(-1));
        yaw_ins1.add_val(i_str, (-pi/2+rollpitchyawins_data.at(i)[2])*(-1));
        auto res_d=apply_rotation_new_order(rollpitchyawins_data.at(i)[0], (-pi+rollpitchyawins_data.at(i)[1])*(-1), (-pi/2+rollpitchyawins_data.at(i)[2])*(-1), point_cloud_data.at(i)[0]-uav_data.at(i)[0], point_cloud_data.at(i)[1]-uav_data.at(i)[1], point_cloud_data.at(i)[2]-uav_data.at(i)[2]);
        x1i.add_val(i_str,res_d[0]);
        y1i.add_val(i_str,res_d[1]);
        z1i.add_val(i_str,res_d[2]);
    }
    

    indices N1("N1"),N2("N2");
    
    int n1 = x1.get_dim();
    int n2 = x2.get_dim();
    DebugOff("n1 = " << n1 << endl);
    DebugOff("n2 = " << n2 << endl);
    
    N1 = range(1,n1);
    N2 = range(1,n2);
    
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
    var<> xb_d("xb_d"), yb_d("yb_d"), zb_d("zb_d");
    Reg->add(xb_d.in(N1), yb_d.in(N1), zb_d.in(N1));
    var<> xb_m("xb_m"), yb_m("yb_m"), zb_m("zb_m");
    Reg->add(xb_m.in(N2), yb_m.in(N2), zb_m.in(N2));
    
    auto ids1 = theta11.repeat_id(N1.size());
    
    Constraint<> xd_bore("xd_bore");
    xd_bore=x1i*theta11.in(ids1)+y1i*theta12.in(ids1)+z1i*theta13.in(ids1)-xb_d;
    Reg->add(xd_bore.in(N1)==0);
    
    Constraint<> yd_bore("yd_bore");
    yd_bore=x1i*theta21.in(ids1)+y1i*theta22.in(ids1)+z1i*theta23.in(ids1)-yb_d;
    Reg->add(yd_bore.in(N1)==0);

    
    Constraint<> zd_bore("zd_bore");
    zd_bore=x1i*theta31.in(ids1)+y1i*theta32.in(ids1)+z1i*theta33.in(ids1)-zb_d;
    Reg->add(zd_bore.in(N1)==0);
    
  
   
    
    Constraint<> xd_ins_inv("xd_ins_inv");
    xd_ins_inv=xb_d*cos(roll_ins1)*cos(yaw_ins1)+yb_d*(cos(pitch_ins1)*sin(yaw_ins1) +cos(yaw_ins1)*sin(roll_ins1)*sin(pitch_ins1))+zb_d*(sin(pitch_ins1)*sin(yaw_ins1)-cos(pitch_ins1)*cos(yaw_ins1)*sin(roll_ins1))+x_uav1-new_xm-x_diff;
    Reg->add(xd_ins_inv.in(N1)==0);
    
    Constraint<> yd_ins_inv("yd_ins_inv");
    yd_ins_inv=xb_d*(-1)*cos(roll_ins1)*sin(yaw_ins1)+yb_d*(cos(pitch_ins1)*cos(yaw_ins1) -sin(roll_ins1)*sin(pitch_ins1)*sin(yaw_ins1))+zb_d*(cos(yaw_ins1)*sin(pitch_ins1)+cos(pitch_ins1)*sin(roll_ins1)*sin(yaw_ins1))+y_uav1-new_ym-y_diff;
    Reg->add(yd_ins_inv.in(N1)==0);
    
    Constraint<> zd_ins_inv("zd_ins_inv");
    zd_ins_inv=xb_d*sin(roll_ins1)+yb_d*(-1)*cos(roll_ins1)*sin(pitch_ins1)+zb_d*(cos(roll_ins1)*cos(pitch_ins1))+z_uav1-new_zm-z_diff;
    Reg->add(zd_ins_inv.in(N1)==0);
 
    auto ids2 = theta11.repeat_id(N2.size());
    
    Constraint<> xm_bore("xm_bore");
    xm_bore=x2i*theta11.in(ids2)+y2i*theta12.in(ids2)+z2i*theta13.in(ids2)-xb_m;
    Reg->add(xm_bore.in(N2)==0);
    
    Constraint<> ym_bore("ym_bore");
    ym_bore=x2i*theta21.in(ids2)+y2i*theta22.in(ids2)+z2i*theta23.in(ids2)-yb_m;
    Reg->add(ym_bore.in(N2)==0);

    
    Constraint<> zm_bore("zm_bore");
    zm_bore=x2i*theta31.in(ids2)+y2i*theta32.in(ids2)+z2i*theta33.in(ids2)-zb_m;
    Reg->add(zm_bore.in(N2)==0);
    

    
    Constraint<> xm_ins_inv("xm_ins_inv");
    xm_ins_inv=xb_m*cos(roll_ins2)*cos(yaw_ins2)+yb_m*(cos(pitch_ins2)*sin(yaw_ins2) +cos(yaw_ins2)*sin(roll_ins2)*sin(pitch_ins2))+zb_m*(sin(pitch_ins2)*sin(yaw_ins2)-cos(pitch_ins2)*cos(yaw_ins2)*sin(roll_ins2))+x_uav2-xm;
    Reg->add(xm_ins_inv.in(N2)==0);
    
    
    Constraint<> ym_ins_inv("ym_ins_inv");
    ym_ins_inv=xb_m*(-1)*cos(roll_ins2)*sin(yaw_ins2)+yb_m*(cos(pitch_ins2)*cos(yaw_ins2) -sin(roll_ins2)*sin(pitch_ins2)*sin(yaw_ins2))+zb_m*(cos(yaw_ins2)*sin(pitch_ins2)+cos(pitch_ins2)*sin(roll_ins2)*sin(yaw_ins2))+y_uav2-ym;
    Reg->add(ym_ins_inv.in(N2)==0);
    
    Constraint<> zm_ins_inv("zm_ins_inv");
    zm_ins_inv=xb_m*sin(roll_ins2)+yb_m*(-1)*cos(roll_ins2)*sin(pitch_ins2)+zb_m*(cos(roll_ins2)*cos(pitch_ins2))+z_uav2-zm;
    Reg->add(zm_ins_inv.in(N2)==0);
    


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
    
    func<> cosr_f = cos(roll);cosr_f.eval_all();
    func<> sinr_f = sin(roll);sinr_f.eval_all();
    func<> cosp_f = cos(pitch);cosp_f.eval_all();
    func<> sinp_f = sin(pitch);sinp_f.eval_all();
    func<> cosy_f = cos(yaw);cosy_f.eval_all();
    func<> siny_f = sin(yaw);siny_f.eval_all();
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
    
    Constraint<> trigR_NC("trigR_NC");
    trigR_NC = pow(cosr,2) + pow(sinr,2);
    Reg->add(trigR_NC.in(range(0,0))==1);
    
    Constraint<> trigP_NC("trigP_NC");
    trigP_NC = pow(cosp,2) + pow(sinp,2);
    Reg->add(trigP_NC.in(range(0,0))==1);
    
    Constraint<> trigY_NC("trigY_NC");
    trigY_NC = pow(cosy,2) + pow(siny,2);
    Reg->add(trigY_NC.in(range(0,0))==1);

    Constraint<> T11("T11");
    T11=theta11.in(range(0,0));
    T11-=cosr*cosy;
    Reg->add(T11.in(range(0,0))==0);
    
    Constraint<> T12("T12");
    T12=theta12.in(range(0,0));
    T12-=(-1)*cosr*siny;
    Reg->add(T12.in(range(0,0))==0);
    
    Constraint<> T13("T13");
    T13=theta13.in(range(0,0));
    T13-=sinr;
    Reg->add(T13.in(range(0,0))==0);
   
    Constraint<> T21("T21");
    T21+=theta21.in(range(0,0));
    T21-=cosp*siny;
    T21-=cosy_sinr*sinp;
    Reg->add(T21.in(range(0,0))==0);
    
    Constraint<> T22("T22");
    T22+=theta22.in(range(0,0));
    T22-=cosp*cosy;
    T22-=(-1)*sinp*siny_sinr;
    Reg->add(T22.in(range(0,0))==0);
    
    Constraint<> T23("T23");
    T23+=theta23.in(range(0,0));
    T23-=(-1)*cosr*sinp;
    Reg->add(T23.in(range(0,0))==0);
    
    Constraint<> T31("T31");
    T31+=theta31.in(range(0,0));
    T31-=sinp*siny;
    T31-=(-1)*cosp*cosy_sinr;
    Reg->add(T31.in(range(0,0))==0);
    
    Constraint<> T32("T32");
    T32+=theta32.in(range(0,0));
    T32-=cosy*sinp;
    T32-=cosp*siny_sinr;
    Reg->add(T32.in(range(0,0))==0);
    
    Constraint<> T33("T33");
    T33+=theta33.in(range(0,0));
    T33-=cosr*cosp;
    Reg->add(T33.in(range(0,0))==0);
 
  //  Reg->min(sum(pow(x_diff,2) + pow(y_diff,2) + pow(z_diff,2)));
    
    Reg->min(sum(deltax) + sum(deltay)+sum(deltaz));
    //Reg->print();

   
    return Reg;
}
shared_ptr<Model<double>> Align_L2_model_rotation_trigonometric_scanner(const vector<vector<double>>& input_model_cloud, const vector<vector<double>>& input_data_cloud, const vector<vector<double>>& uav_model, const vector<vector<double>>& uav_data, const vector<vector<double>>& rollpitchyawins_model, const vector<vector<double>>& rollpitchyawins_data, const vector<vector<double>>& input_model_offset, const vector<vector<double>>& input_data_offset, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, indices& cells, param<double> dist_cost)
{
    double angle_max = 0.1;
    int nb_pairs = 0, min_nb_pairs = numeric_limits<int>::max(), max_nb_pairs = 0, av_nb_pairs = 0;

    param<> x_uav1("x_uav1"), y_uav1("y_uav1"), z_uav1("z_uav1");
    param<> x_uav2("x_uav2"), y_uav2("y_uav2"), z_uav2("z_uav2");
    param<> x1i("x1i"), x2i("x2i"), y1i("y1i"), y2i("y2i"), z1i("z1i"), z2i("z2i");
    param<> roll_ins1("roll_ins1"),pitch_ins1("pitch_ins1"), yaw_ins1("yaw_ins1");
    param<> roll_ins2("roll_ins2"), pitch_ins2("pitch_ins2"), yaw_ins2("yaw_ins2");
    param<> rot_sc_mod_x("rot_sc_mod_x"), rot_sc_mod_y("rot_sc_mod_y"),rot_sc_mod_z("rot_sc_mod_z");
    param<> rot_sc_data_x("rot_sc_data_x"), rot_sc_data_y("rot_sc_data_y"),rot_sc_data_z("rot_sc_data_z");
    string j_str, i_str;
 
    for (auto j = 0; j<input_model_cloud.size(); j++) {
        j_str = to_string(j+1);
        x_uav2.add_val(j_str,uav_model.at(j)[0]);
        y_uav2.add_val(j_str,uav_model.at(j)[1]);
        z_uav2.add_val(j_str,uav_model.at(j)[2]);
        roll_ins2.add_val(j_str, rollpitchyawins_model.at(j)[0]);
        pitch_ins2.add_val(j_str, (-pi+rollpitchyawins_model.at(j)[1])*(-1));
        yaw_ins2.add_val(j_str, (-pi/2+rollpitchyawins_model.at(j)[2])*(-1));
        x2i.add_val(j_str,input_model_cloud.at(j)[0]);
        y2i.add_val(j_str,input_model_cloud.at(j)[1]);
        z2i.add_val(j_str,input_model_cloud.at(j)[2]);
        rot_sc_mod_x.add_val(j_str, input_model_offset.at(j)[0]);
        rot_sc_mod_y.add_val(j_str, input_model_offset.at(j)[1]);
        rot_sc_mod_z.add_val(j_str, input_model_offset.at(j)[2]);
    }
    
    for (auto i = 0; i<input_data_cloud.size(); i++) {
        i_str = to_string(i+1);
        x_uav1.add_val(i_str,uav_data.at(i)[0]);
        y_uav1.add_val(i_str,uav_data.at(i)[1]);
        z_uav1.add_val(i_str,uav_data.at(i)[2]);
        roll_ins1.add_val(i_str, rollpitchyawins_data.at(i)[0]);
        pitch_ins1.add_val(i_str, (-pi+rollpitchyawins_data.at(i)[1])*(-1));
        yaw_ins1.add_val(i_str, (-pi/2+rollpitchyawins_data.at(i)[2])*(-1));
        x1i.add_val(i_str,input_data_cloud.at(i)[0]);
        y1i.add_val(i_str,input_data_cloud.at(i)[1]);
        z1i.add_val(i_str,input_data_cloud.at(i)[2]);
        rot_sc_data_x.add_val(i_str, input_data_offset.at(i)[0]);
        rot_sc_data_y.add_val(i_str, input_data_offset.at(i)[1]);
        rot_sc_data_z.add_val(i_str, input_data_offset.at(i)[2]);
    }
    

    indices N1("N1"),N2("N2");
    
    int n1 = x1i.get_dim();
    int n2 = x2i.get_dim();
    DebugOff("n1 = " << n1 << endl);
    DebugOff("n2 = " << n2 << endl);
    
    N1 = range(1,n1);
    N2 = range(1,n2);
    
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
  
    //auto cells = indices(N1,N2);
    
    vector<int> new_model_pts;
    indices new_model_ids;
    string name="Norm2_Align";
    
    auto Reg=make_shared<Model<>>(name);
    
    Reg->add_param(x1i);Reg->add_param(y1i);Reg->add_param(z1i);
    Reg->add_param(x2i);Reg->add_param(y2i);Reg->add_param(z2i);
    var<int> bin("bin",0,1);
    Reg->add(bin.in(cells));
    DebugOff("Added " << cells.size() << " binary variables" << endl);
    
    
    var<> yaw("yaw", yaw_min, yaw_max), pitch("pitch", pitch_min, pitch_max), roll("roll", roll_min, roll_max);
    yaw.in(R(1)); pitch.in(R(1));roll.in(R(1));
    


  

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
    var<> xb_d("xb_d"), yb_d("yb_d"), zb_d("zb_d");
    Reg->add(xb_d.in(N1), yb_d.in(N1), zb_d.in(N1));
    var<> xb_m("xb_m"), yb_m("yb_m"), zb_m("zb_m");
    Reg->add(xb_m.in(N2), yb_m.in(N2), zb_m.in(N2));
    
    auto ids1 = theta11.repeat_id(N1.size());
    
    Constraint<> xd_bore("xd_bore");
    xd_bore=x1i*theta11.in(ids1)+y1i*theta12.in(ids1)+z1i*theta13.in(ids1)-xb_d;
    Reg->add(xd_bore.in(N1)==0);
    
    Constraint<> yd_bore("yd_bore");
    yd_bore=x1i*theta21.in(ids1)+y1i*theta22.in(ids1)+z1i*theta23.in(ids1)-yb_d;
    Reg->add(yd_bore.in(N1)==0);

    
    Constraint<> zd_bore("zd_bore");
    zd_bore=x1i*theta31.in(ids1)+y1i*theta32.in(ids1)+z1i*theta33.in(ids1)-zb_d;
    Reg->add(zd_bore.in(N1)==0);
    
  
   
    
    Constraint<> xd_ins_inv("xd_ins_inv");
    xd_ins_inv=xb_d*cos(roll_ins1)*cos(yaw_ins1)+yb_d*(cos(pitch_ins1)*sin(yaw_ins1) +cos(yaw_ins1)*sin(roll_ins1)*sin(pitch_ins1))+zb_d*(sin(pitch_ins1)*sin(yaw_ins1)-cos(pitch_ins1)*cos(yaw_ins1)*sin(roll_ins1))+x_uav1+rot_sc_data_x-new_xm-x_diff;
    Reg->add(xd_ins_inv.in(N1)==0);
    
    Constraint<> yd_ins_inv("yd_ins_inv");
    yd_ins_inv=xb_d*(-1)*cos(roll_ins1)*sin(yaw_ins1)+yb_d*(cos(pitch_ins1)*cos(yaw_ins1) -sin(roll_ins1)*sin(pitch_ins1)*sin(yaw_ins1))+zb_d*(cos(yaw_ins1)*sin(pitch_ins1)+cos(pitch_ins1)*sin(roll_ins1)*sin(yaw_ins1))+y_uav1+rot_sc_data_y-new_ym-y_diff;
    Reg->add(yd_ins_inv.in(N1)==0);
    
    Constraint<> zd_ins_inv("zd_ins_inv");
    zd_ins_inv=xb_d*sin(roll_ins1)+yb_d*(-1)*cos(roll_ins1)*sin(pitch_ins1)+zb_d*(cos(roll_ins1)*cos(pitch_ins1))+z_uav1+rot_sc_data_z-new_zm-z_diff;
    Reg->add(zd_ins_inv.in(N1)==0);
 
    auto ids2 = theta11.repeat_id(N2.size());
    
    Constraint<> xm_bore("xm_bore");
    xm_bore=x2i*theta11.in(ids2)+y2i*theta12.in(ids2)+z2i*theta13.in(ids2)-xb_m;
    Reg->add(xm_bore.in(N2)==0);
    
    Constraint<> ym_bore("ym_bore");
    ym_bore=x2i*theta21.in(ids2)+y2i*theta22.in(ids2)+z2i*theta23.in(ids2)-yb_m;
    Reg->add(ym_bore.in(N2)==0);

    
    Constraint<> zm_bore("zm_bore");
    zm_bore=x2i*theta31.in(ids2)+y2i*theta32.in(ids2)+z2i*theta33.in(ids2)-zb_m;
    Reg->add(zm_bore.in(N2)==0);
    

    
    Constraint<> xm_ins_inv("xm_ins_inv");
    xm_ins_inv=xb_m*cos(roll_ins2)*cos(yaw_ins2)+yb_m*(cos(pitch_ins2)*sin(yaw_ins2) +cos(yaw_ins2)*sin(roll_ins2)*sin(pitch_ins2))+zb_m*(sin(pitch_ins2)*sin(yaw_ins2)-cos(pitch_ins2)*cos(yaw_ins2)*sin(roll_ins2))+x_uav2+rot_sc_mod_x-xm;
    Reg->add(xm_ins_inv.in(N2)==0);
    
    
    Constraint<> ym_ins_inv("ym_ins_inv");
    ym_ins_inv=xb_m*(-1)*cos(roll_ins2)*sin(yaw_ins2)+yb_m*(cos(pitch_ins2)*cos(yaw_ins2) -sin(roll_ins2)*sin(pitch_ins2)*sin(yaw_ins2))+zb_m*(cos(yaw_ins2)*sin(pitch_ins2)+cos(pitch_ins2)*sin(roll_ins2)*sin(yaw_ins2))+y_uav2+rot_sc_mod_y-ym;
    Reg->add(ym_ins_inv.in(N2)==0);
    
    Constraint<> zm_ins_inv("zm_ins_inv");
    zm_ins_inv=xb_m*sin(roll_ins2)+yb_m*(-1)*cos(roll_ins2)*sin(pitch_ins2)+zb_m*(cos(roll_ins2)*cos(pitch_ins2))+z_uav2+rot_sc_mod_z-zm;
    Reg->add(zm_ins_inv.in(N2)==0);
    


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
    
    
    if(dist_cost._val->size()!=0){
        Constraint<> delta_cost("delta_cost");
        delta_cost=product(dist_cost.in(idsij), bin.in_matrix(1,1))-deltax-deltay-deltaz;
        Reg->add(delta_cost.in(N1)<=0);
    }
    
    func<> cosr_f = cos(roll);cosr_f.eval_all();
    func<> sinr_f = sin(roll);sinr_f.eval_all();
    func<> cosp_f = cos(pitch);cosp_f.eval_all();
    func<> sinp_f = sin(pitch);sinp_f.eval_all();
    func<> cosy_f = cos(yaw);cosy_f.eval_all();
    func<> siny_f = sin(yaw);siny_f.eval_all();
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
    
    Constraint<> trigR_NC("trigR_NC");
    trigR_NC = pow(cosr,2) + pow(sinr,2);
    Reg->add(trigR_NC.in(range(0,0))==1);
    
    Constraint<> trigP_NC("trigP_NC");
    trigP_NC = pow(cosp,2) + pow(sinp,2);
    Reg->add(trigP_NC.in(range(0,0))==1);
    
    Constraint<> trigY_NC("trigY_NC");
    trigY_NC = pow(cosy,2) + pow(siny,2);
    Reg->add(trigY_NC.in(range(0,0))==1);

    Constraint<> T11("T11");
    T11=theta11.in(range(0,0));
    T11-=cosr*cosy;
    Reg->add(T11.in(range(0,0))==0);
    
    Constraint<> T12("T12");
    T12=theta12.in(range(0,0));
    T12-=(-1)*cosr*siny;
    Reg->add(T12.in(range(0,0))==0);
    
    Constraint<> T13("T13");
    T13=theta13.in(range(0,0));
    T13-=sinr;
    Reg->add(T13.in(range(0,0))==0);
   
    Constraint<> T21("T21");
    T21+=theta21.in(range(0,0));
    T21-=cosp*siny;
    T21-=cosy_sinr*sinp;
    Reg->add(T21.in(range(0,0))==0);
    
    Constraint<> T22("T22");
    T22+=theta22.in(range(0,0));
    T22-=cosp*cosy;
    T22-=(-1)*sinp*siny_sinr;
    Reg->add(T22.in(range(0,0))==0);
    
    Constraint<> T23("T23");
    T23+=theta23.in(range(0,0));
    T23-=(-1)*cosr*sinp;
    Reg->add(T23.in(range(0,0))==0);
    
    Constraint<> T31("T31");
    T31+=theta31.in(range(0,0));
    T31-=sinp*siny;
    T31-=(-1)*cosp*cosy_sinr;
    Reg->add(T31.in(range(0,0))==0);
    
    Constraint<> T32("T32");
    T32+=theta32.in(range(0,0));
    T32-=cosy*sinp;
    T32-=cosp*siny_sinr;
    Reg->add(T32.in(range(0,0))==0);
    
    Constraint<> T33("T33");
    T33+=theta33.in(range(0,0));
    T33-=cosr*cosp;
    Reg->add(T33.in(range(0,0))==0);
 
  //  Reg->min(sum(pow(x_diff,2) + pow(y_diff,2) + pow(z_diff,2)));
    
    Reg->min(sum(deltax) + sum(deltay)+sum(deltaz));
    //Reg->print();

   
    return Reg;
}
shared_ptr<Model<double>> Align_L2_model_rotation_orthogonal(const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data, const vector<vector<double>>& uav_model, const vector<vector<double>>& uav_data, const vector<vector<double>>& rollpitchyawins_model, const vector<vector<double>>& rollpitchyawins_data, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, indices& cells, param<double> dist_cost)
{
    double angle_max = 0.1;
    int nb_pairs = 0, min_nb_pairs = numeric_limits<int>::max(), max_nb_pairs = 0, av_nb_pairs = 0;

    param<> x1("x1"), x2("x2"), y1("y1"), y2("y2"), z1("z1"), z2("z2");
    param<> x_uav1("x_uav1"), y_uav1("y_uav1"), z_uav1("z_uav1");
    param<> x_uav2("x_uav2"), y_uav2("y_uav2"), z_uav2("z_uav2");
    param<> x1i("x1i"), x2i("x2i"), y1i("y1i"), y2i("y2i"), z1i("z1i"), z2i("z2i");
    param<> roll_ins1("roll_ins1"),pitch_ins1("pitch_ins1"), yaw_ins1("yaw_ins1");
    param<> roll_ins2("roll_ins2"), pitch_ins2("pitch_ins2"), yaw_ins2("yaw_ins2");
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
        roll_ins2.add_val(j_str, rollpitchyawins_model.at(j)[0]);
        pitch_ins2.add_val(j_str, (-pi+rollpitchyawins_model.at(j)[1])*(-1));
        yaw_ins2.add_val(j_str, (-pi/2+rollpitchyawins_model.at(j)[2])*(-1));
        auto res_m=apply_rotation_new_order(rollpitchyawins_model.at(j)[0], (-pi+rollpitchyawins_model.at(j)[1])*(-1), (-pi/2+rollpitchyawins_model.at(j)[2])*(-1), point_cloud_model.at(j)[0]-uav_model.at(j)[0], point_cloud_model.at(j)[1]-uav_model.at(j)[1], point_cloud_model.at(j)[2]-uav_model.at(j)[2]);
        x2i.add_val(j_str,res_m[0]);
        y2i.add_val(j_str,res_m[1]);
        z2i.add_val(j_str,res_m[2]);
    }
    
    for (auto i = 0; i<point_cloud_data.size(); i++) {
        i_str = to_string(i+1);
        x_uav1.add_val(i_str,uav_data.at(i)[0]);
        x1.add_val(i_str,point_cloud_data.at(i)[0]);
        y_uav1.add_val(i_str,uav_data.at(i)[1]);
        y1.add_val(i_str,point_cloud_data.at(i)[1]);
        z_uav1.add_val(i_str,uav_data.at(i)[2]);
        z1.add_val(i_str,point_cloud_data.at(i)[2]);
        roll_ins1.add_val(i_str, rollpitchyawins_data.at(i)[0]);
        pitch_ins1.add_val(i_str, (-pi+rollpitchyawins_data.at(i)[1])*(-1));
        yaw_ins1.add_val(i_str, (-pi/2+rollpitchyawins_data.at(i)[2])*(-1));
        auto res_d=apply_rotation_new_order(rollpitchyawins_data.at(i)[0], (-pi+rollpitchyawins_data.at(i)[1])*(-1), (-pi/2+rollpitchyawins_data.at(i)[2])*(-1), point_cloud_data.at(i)[0]-uav_data.at(i)[0], point_cloud_data.at(i)[1]-uav_data.at(i)[1], point_cloud_data.at(i)[2]-uav_data.at(i)[2]);
        x1i.add_val(i_str,res_d[0]);
        y1i.add_val(i_str,res_d[1]);
        z1i.add_val(i_str,res_d[2]);
    }
    

    indices N1("N1"),N2("N2");
    
    int n1 = x1.get_dim();
    int n2 = x2.get_dim();
    DebugOn("n1 = " << n1 << endl);
    DebugOn("n2 = " << n2 << endl);
    
    N1 = range(1,n1);
    N2 = range(1,n2);
    
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
    var<> xb_d("xb_d"), yb_d("yb_d"), zb_d("zb_d");
    Reg->add(xb_d.in(N1), yb_d.in(N1), zb_d.in(N1));
    var<> xb_m("xb_m"), yb_m("yb_m"), zb_m("zb_m");
    Reg->add(xb_m.in(N2), yb_m.in(N2), zb_m.in(N2));
    
    auto ids1 = theta11.repeat_id(N1.size());
    
    Constraint<> xd_bore("xd_bore");
    xd_bore=x1i*theta11.in(ids1)+y1i*theta12.in(ids1)+z1i*theta13.in(ids1)-xb_d;
    Reg->add(xd_bore.in(N1)==0);
    
    Constraint<> yd_bore("yd_bore");
    yd_bore=x1i*theta21.in(ids1)+y1i*theta22.in(ids1)+z1i*theta23.in(ids1)-yb_d;
    Reg->add(yd_bore.in(N1)==0);

    
    Constraint<> zd_bore("zd_bore");
    zd_bore=x1i*theta31.in(ids1)+y1i*theta32.in(ids1)+z1i*theta33.in(ids1)-zb_d;
    Reg->add(zd_bore.in(N1)==0);
    
  
   
    
    Constraint<> xd_ins_inv("xd_ins_inv");
    xd_ins_inv=xb_d*cos(roll_ins1)*cos(yaw_ins1)+yb_d*(cos(pitch_ins1)*sin(yaw_ins1) +cos(yaw_ins1)*sin(roll_ins1)*sin(pitch_ins1))+zb_d*(sin(pitch_ins1)*sin(yaw_ins1)-cos(pitch_ins1)*cos(yaw_ins1)*sin(roll_ins1))+x_uav1-new_xm-x_diff;
    Reg->add(xd_ins_inv.in(N1)==0);
    
    Constraint<> yd_ins_inv("yd_ins_inv");
    yd_ins_inv=xb_d*(-1)*cos(roll_ins1)*sin(yaw_ins1)+yb_d*(cos(pitch_ins1)*cos(yaw_ins1) -sin(roll_ins1)*sin(pitch_ins1)*sin(yaw_ins1))+zb_d*(cos(yaw_ins1)*sin(pitch_ins1)+cos(pitch_ins1)*sin(roll_ins1)*sin(yaw_ins1))+y_uav1-new_ym-y_diff;
    Reg->add(yd_ins_inv.in(N1)==0);
    
    Constraint<> zd_ins_inv("zd_ins_inv");
    zd_ins_inv=xb_d*sin(roll_ins1)+yb_d*(-1)*cos(roll_ins1)*sin(pitch_ins1)+zb_d*(cos(roll_ins1)*cos(pitch_ins1))+z_uav1-new_zm-z_diff;
    Reg->add(zd_ins_inv.in(N1)==0);
 
    auto ids2 = theta11.repeat_id(N2.size());
    
    Constraint<> xm_bore("xm_bore");
    xm_bore=x2i*theta11.in(ids2)+y2i*theta12.in(ids2)+z2i*theta13.in(ids2)-xb_m;
    Reg->add(xm_bore.in(N2)==0);
    
    Constraint<> ym_bore("ym_bore");
    ym_bore=x2i*theta21.in(ids2)+y2i*theta22.in(ids2)+z2i*theta23.in(ids2)-yb_m;
    Reg->add(ym_bore.in(N2)==0);

    
    Constraint<> zm_bore("zm_bore");
    zm_bore=x2i*theta31.in(ids2)+y2i*theta32.in(ids2)+z2i*theta33.in(ids2)-zb_m;
    Reg->add(zm_bore.in(N2)==0);
    

    
    Constraint<> xm_ins_inv("xm_ins_inv");
    xm_ins_inv=xb_m*cos(roll_ins2)*cos(yaw_ins2)+yb_m*(cos(pitch_ins2)*sin(yaw_ins2) +cos(yaw_ins2)*sin(roll_ins2)*sin(pitch_ins2))+zb_m*(sin(pitch_ins2)*sin(yaw_ins2)-cos(pitch_ins2)*cos(yaw_ins2)*sin(roll_ins2))+x_uav2-xm;
    Reg->add(xm_ins_inv.in(N2)==0);
    
    
    Constraint<> ym_ins_inv("ym_ins_inv");
    ym_ins_inv=xb_m*(-1)*cos(roll_ins2)*sin(yaw_ins2)+yb_m*(cos(pitch_ins2)*cos(yaw_ins2) -sin(roll_ins2)*sin(pitch_ins2)*sin(yaw_ins2))+zb_m*(cos(yaw_ins2)*sin(pitch_ins2)+cos(pitch_ins2)*sin(roll_ins2)*sin(yaw_ins2))+y_uav2-ym;
    Reg->add(ym_ins_inv.in(N2)==0);
    
    Constraint<> zm_ins_inv("zm_ins_inv");
    zm_ins_inv=xb_m*sin(roll_ins2)+yb_m*(-1)*cos(roll_ins2)*sin(pitch_ins2)+zb_m*(cos(roll_ins2)*cos(pitch_ins2))+z_uav2-zm;
    Reg->add(zm_ins_inv.in(N2)==0);
    


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
    
    var<> minor1("minor1"), minor2("minor2"), minor3("minor3");
    Reg->add(minor1.in(R(1)), minor2.in(R(1)), minor3.in(R(1)));
    
    Constraint<> row1("row1");
    row1 = pow(theta11,2)+pow(theta12,2)+pow(theta13,2);
    Reg->add(row1==1);
    Constraint<> row2("row2");
    row2 = pow(theta21,2)+pow(theta22,2)+pow(theta23,2);
    Reg->add(row2==1);
    Constraint<> row3("row3");
    row3 = pow(theta31,2)+pow(theta32,2)+pow(theta33,2);
    Reg->add(row3==1);
    
    Constraint<> row1_col2("row1_col2");
    row1_col2 = theta11*theta21+theta12*theta22+theta13*theta23;
    Reg->add(row1_col2==0);
    
    Constraint<> row1_col3("row1_col3");
    row1_col3 = theta11*theta31+theta12*theta32+theta13*theta33;
    Reg->add(row1_col3==0);
    
    Constraint<> row2_col3("row2_col3");
    row2_col3 = theta21*theta31+theta22*theta32+theta23*theta33;
    Reg->add(row2_col3==0);
    
    Constraint<> minor1_Def("minor1_Def");
    minor1_Def=theta22*theta33-theta23*theta32-minor1;
    Reg->add(minor1_Def==0);
    
    Constraint<> minor2_Def("minor2_Def");
    minor2_Def=theta21*theta33-theta31*theta23-minor2;
    Reg->add(minor2_Def==0);
    
    Constraint<> minor3_Def("minor3_Def");
    minor3_Def=theta21*theta32-theta31*theta22-minor3;
    Reg->add(minor3_Def==0);
    
    Constraint<> Det_Def("Det_Def");
    Det_Def=theta11*minor1-theta12*minor2+theta13*minor3;
    Reg->add(Det_Def==1);
    
    

   
 
  //  Reg->min(sum(pow(x_diff,2) + pow(y_diff,2) + pow(z_diff,2)));
    
    Reg->min(sum(deltax) + sum(deltay)+sum(deltaz));
    //Reg->print();

   
    return Reg;
}
shared_ptr<Model<double>> Align_L1_model_rotation_neworder(const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data, const vector<vector<double>>& uav_model, const vector<vector<double>>& uav_data, const vector<vector<double>>& rollpitchyawins_model, const vector<vector<double>>& rollpitchyawins_data, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, indices& cells, param<double> dist_cost)
{
    double angle_max = 0.1;
    int nb_pairs = 0, min_nb_pairs = numeric_limits<int>::max(), max_nb_pairs = 0, av_nb_pairs = 0;

    param<> x1("x1"), x2("x2"), y1("y1"), y2("y2"), z1("z1"), z2("z2");
    param<> x_uav1("x_uav1"), y_uav1("y_uav1"), z_uav1("z_uav1");
    param<> x_uav2("x_uav2"), y_uav2("y_uav2"), z_uav2("z_uav2");
    param<> x1i("x1i"), x2i("x2i"), y1i("y1i"), y2i("y2i"), z1i("z1i"), z2i("z2i");
    param<> roll_ins1("roll_ins1"),pitch_ins1("pitch_ins1"), yaw_ins1("yaw_ins1");
    param<> roll_ins2("roll_ins2"), pitch_ins2("pitch_ins2"), yaw_ins2("yaw_ins2");
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
        roll_ins2.add_val(j_str, rollpitchyawins_model.at(j)[0]);
        pitch_ins2.add_val(j_str, (-pi+rollpitchyawins_model.at(j)[1])*(-1));
        yaw_ins2.add_val(j_str, (-pi/2+rollpitchyawins_model.at(j)[2])*(-1));
        auto res_m=apply_rotation_new_order(rollpitchyawins_model.at(j)[0], (-pi+rollpitchyawins_model.at(j)[1])*(-1), (-pi/2+rollpitchyawins_model.at(j)[2])*(-1), point_cloud_model.at(j)[0]-uav_model.at(j)[0], point_cloud_model.at(j)[1]-uav_model.at(j)[1], point_cloud_model.at(j)[2]-uav_model.at(j)[2]);
        x2i.add_val(j_str,res_m[0]);
        y2i.add_val(j_str,res_m[1]);
        z2i.add_val(j_str,res_m[2]);
    }
    
    for (auto i = 0; i<point_cloud_data.size(); i++) {
        i_str = to_string(i+1);
        x_uav1.add_val(i_str,uav_data.at(i)[0]);
        x1.add_val(i_str,point_cloud_data.at(i)[0]);
        y_uav1.add_val(i_str,uav_data.at(i)[1]);
        y1.add_val(i_str,point_cloud_data.at(i)[1]);
        z_uav1.add_val(i_str,uav_data.at(i)[2]);
        z1.add_val(i_str,point_cloud_data.at(i)[2]);
        roll_ins1.add_val(i_str, rollpitchyawins_data.at(i)[0]);
        pitch_ins1.add_val(i_str, (-pi+rollpitchyawins_data.at(i)[1])*(-1));
        yaw_ins1.add_val(i_str, (-pi/2+rollpitchyawins_data.at(i)[2])*(-1));
        auto res_d=apply_rotation_new_order(rollpitchyawins_data.at(i)[0], (-pi+rollpitchyawins_data.at(i)[1])*(-1), (-pi/2+rollpitchyawins_data.at(i)[2])*(-1), point_cloud_data.at(i)[0]-uav_data.at(i)[0], point_cloud_data.at(i)[1]-uav_data.at(i)[1], point_cloud_data.at(i)[2]-uav_data.at(i)[2]);
        x1i.add_val(i_str,res_d[0]);
        y1i.add_val(i_str,res_d[1]);
        z1i.add_val(i_str,res_d[2]);
    }
    

    indices N1("N1"),N2("N2");
    
    int n1 = x1.get_dim();
    int n2 = x2.get_dim();
    DebugOn("n1 = " << n1 << endl);
    DebugOn("n2 = " << n2 << endl);
    
    N1 = range(1,n1);
    N2 = range(1,n2);
    
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
    var<> xb_d("xb_d"), yb_d("yb_d"), zb_d("zb_d");
    Reg->add(xb_d.in(N1), yb_d.in(N1), zb_d.in(N1));
    var<> xb_m("xb_m"), yb_m("yb_m"), zb_m("zb_m");
    Reg->add(xb_m.in(N2), yb_m.in(N2), zb_m.in(N2));
    
    auto ids1 = theta11.repeat_id(N1.size());
    
    Constraint<> xd_bore("xd_bore");
    xd_bore=x1i*theta11.in(ids1)+y1i*theta12.in(ids1)+z1i*theta13.in(ids1)-xb_d;
    Reg->add(xd_bore.in(N1)==0);
    
    Constraint<> yd_bore("yd_bore");
    yd_bore=x1i*theta21.in(ids1)+y1i*theta22.in(ids1)+z1i*theta23.in(ids1)-yb_d;
    Reg->add(yd_bore.in(N1)==0);

    
    Constraint<> zd_bore("zd_bore");
    zd_bore=x1i*theta31.in(ids1)+y1i*theta32.in(ids1)+z1i*theta33.in(ids1)-zb_d;
    Reg->add(zd_bore.in(N1)==0);
    
  
   
    
    Constraint<> xd_ins_inv1("xd_ins_inv1");
    xd_ins_inv1=x_diff;
    xd_ins_inv1+=xb_d*cos(roll_ins1)*cos(yaw_ins1)+yb_d*(cos(pitch_ins1)*sin(yaw_ins1) +cos(yaw_ins1)*sin(roll_ins1)*sin(pitch_ins1))+zb_d*(sin(pitch_ins1)*sin(yaw_ins1)-cos(pitch_ins1)*cos(yaw_ins1)*sin(roll_ins1))+x_uav1-new_xm;
    Reg->add(xd_ins_inv1.in(N1)>=0);
    
    Constraint<> xd_ins_inv2("xd_ins_inv2");
    xd_ins_inv2=x_diff;
    xd_ins_inv2-=xb_d*cos(roll_ins1)*cos(yaw_ins1)+yb_d*(cos(pitch_ins1)*sin(yaw_ins1) +cos(yaw_ins1)*sin(roll_ins1)*sin(pitch_ins1))+zb_d*(sin(pitch_ins1)*sin(yaw_ins1)-cos(pitch_ins1)*cos(yaw_ins1)*sin(roll_ins1))+x_uav1-new_xm;
    Reg->add(xd_ins_inv2.in(N1)>=0);
    
    Constraint<> yd_ins_inv1("yd_ins_inv1");
    yd_ins_inv1=y_diff;
    yd_ins_inv1+=xb_d*(-1)*cos(roll_ins1)*sin(yaw_ins1)+yb_d*(cos(pitch_ins1)*cos(yaw_ins1) -sin(roll_ins1)*sin(pitch_ins1)*sin(yaw_ins1))+zb_d*(cos(yaw_ins1)*sin(pitch_ins1)+cos(pitch_ins1)*sin(roll_ins1)*sin(yaw_ins1))+y_uav1-new_ym;
    Reg->add(yd_ins_inv1.in(N1)>=0);
    
    Constraint<> yd_ins_inv2("yd_ins_inv2");
    yd_ins_inv2=y_diff;
    yd_ins_inv1-=xb_d*(-1)*cos(roll_ins1)*sin(yaw_ins1)+yb_d*(cos(pitch_ins1)*cos(yaw_ins1) -sin(roll_ins1)*sin(pitch_ins1)*sin(yaw_ins1))+zb_d*(cos(yaw_ins1)*sin(pitch_ins1)+cos(pitch_ins1)*sin(roll_ins1)*sin(yaw_ins1))+y_uav1-new_ym;
    Reg->add(yd_ins_inv2.in(N1)>=0);
    
    Constraint<> zd_ins_inv1("zd_ins_inv1");
    zd_ins_inv1=z_diff;
    zd_ins_inv1+=xb_d*sin(roll_ins1)+yb_d*(-1)*cos(roll_ins1)*sin(pitch_ins1)+zb_d*(cos(roll_ins1)*cos(pitch_ins1))+z_uav1-new_zm;
    Reg->add(zd_ins_inv1.in(N1)>=0);
    
    Constraint<> zd_ins_inv2("zd_ins_inv2");
    zd_ins_inv2=z_diff;
    zd_ins_inv2-=xb_d*sin(roll_ins1)+yb_d*(-1)*cos(roll_ins1)*sin(pitch_ins1)+zb_d*(cos(roll_ins1)*cos(pitch_ins1))+z_uav1-new_zm;
    Reg->add(zd_ins_inv2.in(N1)>=0);
 
    auto ids2 = theta11.repeat_id(N2.size());
    
    Constraint<> xm_bore("xm_bore");
    xm_bore=x2i*theta11.in(ids2)+y2i*theta12.in(ids2)+z2i*theta13.in(ids2)-xb_m;
    Reg->add(xm_bore.in(N2)==0);
    
    Constraint<> ym_bore("ym_bore");
    ym_bore=x2i*theta21.in(ids2)+y2i*theta22.in(ids2)+z2i*theta23.in(ids2)-yb_m;
    Reg->add(ym_bore.in(N2)==0);

    
    Constraint<> zm_bore("zm_bore");
    zm_bore=x2i*theta31.in(ids2)+y2i*theta32.in(ids2)+z2i*theta33.in(ids2)-zb_m;
    Reg->add(zm_bore.in(N2)==0);
    

    
    Constraint<> xm_ins_inv("xm_ins_inv");
    xm_ins_inv=xb_m*cos(roll_ins2)*cos(yaw_ins2)+yb_m*(cos(pitch_ins2)*sin(yaw_ins2) +cos(yaw_ins2)*sin(roll_ins2)*sin(pitch_ins2))+zb_m*(sin(pitch_ins2)*sin(yaw_ins2)-cos(pitch_ins2)*cos(yaw_ins2)*sin(roll_ins2))+x_uav2-xm;
    Reg->add(xm_ins_inv.in(N2)==0);
    
    
    Constraint<> ym_ins_inv("ym_ins_inv");
    ym_ins_inv=xb_m*(-1)*cos(roll_ins2)*sin(yaw_ins2)+yb_m*(cos(pitch_ins2)*cos(yaw_ins2) -sin(roll_ins2)*sin(pitch_ins2)*sin(yaw_ins2))+zb_m*(cos(yaw_ins2)*sin(pitch_ins2)+cos(pitch_ins2)*sin(roll_ins2)*sin(yaw_ins2))+y_uav2-ym;
    Reg->add(ym_ins_inv.in(N2)==0);
    
    Constraint<> zm_ins_inv("zm_ins_inv");
    zm_ins_inv=xb_m*sin(roll_ins2)+yb_m*(-1)*cos(roll_ins2)*sin(pitch_ins2)+zb_m*(cos(roll_ins2)*cos(pitch_ins2))+z_uav2-zm;
    Reg->add(zm_ins_inv.in(N2)==0);
    


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
#endif /* Lower_Bound_h */
