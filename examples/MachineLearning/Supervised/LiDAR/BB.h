//
//  BB.h
//  Gravity
//
//  Created by Smitha on 4/2/22.
//

#ifndef BB_h
#define BB_h
#ifdef USE_GJK
extern "C" {
#include "openGJK.h"
}
double distance_polytopes_gjk(vector<vector<double>>& vec1, vector<vector<double>>& vec2){
    
    /* Squared distance computed by openGJK.                                 */
    double dd;
    /* Structure of simplex used by openGJK.                                 */
    struct simplex  s;
    /* Number of vertices defining body 1 and body 2, respectively.          */
    int    nvrtx1, nvrtx2;
    
    double **arr1 = (double **)malloc(vec1.size() * sizeof(double *));
    for (int i=0; i<vec1.size(); i++){
        arr1[i] = (double *)malloc(3 * sizeof(double));
    }
    
    /* Read and store vertices' coordinates. */
    for (auto i = 0; i<vec1.size(); i++)
    {
        arr1[i][0]=vec1[i][0];
        arr1[i][1]=vec1[i][1];
        arr1[i][2]=vec1[i][2];
    }
    
    double **arr2 = (double **)malloc(vec2.size() * sizeof(double *));
    for (int i=0; i<vec2.size(); i++){
        arr2[i] = (double *)malloc(3 * sizeof(double));
    }
    
    /* Read and store vertices' coordinates. */
    for (auto i = 0; i<vec2.size(); i++)
    {
        arr2[i][0]=vec2[i][0];
        arr2[i][1]=vec2[i][1];
        arr2[i][2]=vec2[i][2];
    }
    
    nvrtx1=vec1.size();
    nvrtx2=vec2.size();
    
    /* Structures of body 1 and body 2, respectively.                        */
    struct bd       bd1;
    struct bd       bd2;
    
    bd1.coord = arr1;
    bd1.numpoints = nvrtx1;
    
    /* Import coordinates of object 2. */
    
    bd2.coord = arr2;
    bd2.numpoints = nvrtx2;
    
    /* Initialise simplex as empty */
    s.nvrtx = 0;
    
    /* For importing openGJK this is Step 3: invoke the GJK procedure. */
    /* Compute squared distance using GJK algorithm. */
    dd = gjk (bd1, bd2, &s);
    
    for (int i=0; i<bd1.numpoints; i++){
        free(bd1.coord[i]);
    }
    free(bd1.coord);
    for (int i=0; i<bd2.numpoints; i++){
        free(bd2.coord[i]);
    }
    free(bd2.coord);
    
    dd*=dd;
    DebugOff("Distance sqrd "<<dd<<endl);
    return dd;
}
#endif

shared_ptr<Model<double>> Reg_L2_model_rotation_trigonometric(const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max,double tx_min, double tx_max, double ty_min, double ty_max, double tz_min, double tz_max, indices& cells, param<double> dist_cost)
{
    double angle_max = 0.1;
    int nb_pairs = 0, min_nb_pairs = numeric_limits<int>::max(), max_nb_pairs = 0, av_nb_pairs = 0;
    
    param<> x1("x1"), x2("x2"), y1("y1"), y2("y2"), z1("z1"), z2("z2");
    
    string j_str, i_str;
    
    for (auto j = 0; j<point_cloud_model.size(); j++) {
        j_str = to_string(j+1);
        x2.add_val(j_str,point_cloud_model.at(j).at(0));
        y2.add_val(j_str,point_cloud_model.at(j).at(1));
        z2.add_val(j_str,point_cloud_model.at(j).at(2));
    }
    
    for (auto i = 0; i<point_cloud_data.size(); i++) {
        i_str = to_string(i+1);
        x1.add_val(i_str,point_cloud_data.at(i)[0]);
        y1.add_val(i_str,point_cloud_data.at(i)[1]);
        z1.add_val(i_str,point_cloud_data.at(i)[2]);
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
    
    string name="Norm2_Reg";
    
    auto Reg=make_shared<Model<>>(name);
    
    Reg->add_param(x1);Reg->add_param(y1);Reg->add_param(z1);
    Reg->add_param(x2);Reg->add_param(y2);Reg->add_param(z2);
    var<int> bin("bin",0,1);
    Reg->add(bin.in(cells));
    DebugOff("Added " << cells.size() << " binary variables" << endl);
    
    
    var<> yaw("yaw", yaw_min, yaw_max), pitch("pitch", pitch_min, pitch_max), roll("roll", roll_min, roll_max);
    yaw.in(R(1)); pitch.in(R(1));roll.in(R(1));
    
    
    func<> r11 = cos(roll)*cos(yaw);r11.eval_all();
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
    
    var<> tx("tx", tx_min, tx_max), ty("ty", ty_min, ty_max),tz("tz", tz_min, tz_max);
    Reg->add(tx.in(R(1)),ty.in(R(1)),tz.in(R(1)));
    
    
    Reg->add(theta11.in(R(1)),theta12.in(R(1)),theta13.in(R(1)));
    Reg->add(theta21.in(R(1)),theta22.in(R(1)),theta23.in(R(1)));
    Reg->add(theta31.in(R(1)),theta32.in(R(1)),theta33.in(R(1)));
    
    
    var<> new_xm("new_xm"), new_ym("new_ym"), new_zm("new_zm");
    
    var<> xm("xm"), ym("ym"), zm("zm");
    Reg->add(new_xm.in(N1), new_ym.in(N1), new_zm.in(N1));
    Reg->add(xm.in(N2), ym.in(N2), zm.in(N2));
    
    
    
    
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
    
    
    
    
    
    var<> x_diff("x_diff"), y_diff("y_diff"), z_diff("z_diff");
    Reg->add(x_diff.in(N1), y_diff.in(N1), z_diff.in(N1));
    
    
    auto ids1 = theta11.repeat_id(N1.size());
    
    Constraint<> xd("xd");
    xd=x1*theta11.in(ids1)+y1*theta12.in(ids1)+z1*theta13.in(ids1)+tx.in(ids1)-new_xm-x_diff;
    Reg->add(xd.in(N1)==0);
    
    Constraint<> yd("yd");
    yd=x1*theta21.in(ids1)+y1*theta22.in(ids1)+z1*theta23.in(ids1)+ty.in(ids1)-new_ym-y_diff;
    Reg->add(yd.in(N1)==0);
    
    
    Constraint<> zd("zd");
    zd=x1*theta31.in(ids1)+y1*theta32.in(ids1)+z1*theta33.in(ids1)+tz.in(ids1)-new_zm-z_diff;
    Reg->add(zd.in(N1)==0);
    
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
    T12-=cosy_sinr*sinp;
    T12-=(-1)*cosp*siny;
    Reg->add(T12.in(range(0,0))==0);
    
    Constraint<> T13("T13");
    T13=theta13.in(range(0,0));
    T13-=cosy_sinr*cosp;
    T13-=sinp*siny;
    Reg->add(T13.in(range(0,0))==0);
    
    Constraint<> T21("T21");
    T21+=theta21.in(range(0,0));
    T21-=cosr*siny;
    Reg->add(T21.in(range(0,0))==0);
    
    Constraint<> T22("T22");
    T22+=theta22.in(range(0,0));
    T22-=siny_sinr*sinp;
    T22-=cosp*cosy;
    Reg->add(T22.in(range(0,0))==0);
    
    Constraint<> T23("T23");
    T23+=theta23.in(range(0,0));
    T23-=siny_sinr*cosp;
    T23-=(-1)*sinp*cosy;
    Reg->add(T23.in(range(0,0))==0);
    
    Constraint<> T31("T31");
    T31+=theta31.in(range(0,0));
    T31-=(-1)*sinr;
    Reg->add(T31.in(range(0,0))==0);
    
    Constraint<> T32("T32");
    T32+=theta32.in(range(0,0));
    T32-=cosr*sinp;
    Reg->add(T32.in(range(0,0))==0);
    
    Constraint<> T33("T33");
    T33+=theta33.in(range(0,0));
    T33-=cosr*cosp;
    Reg->add(T33.in(range(0,0))==0);
    
    
    Reg->min(sum(deltax) + sum(deltay)+sum(deltaz));
    
    return Reg;
}
#ifdef USE_GJK
void preprocess_lid(const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data, const vector<vector<vector<double>>>& model_voronoi_vertices, indices& valid_cells_old, indices& new_cells, param<double>& dist_cells, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double yaw_max,double tx_min, double tx_max, double ty_min, double ty_max, double tz_min, double tz_max, double upper_bound, double& prep_time, double& min_cost_sum, string error_type)
{
    prep_time=0;
    indices valid_cells("valid_cells");
    indices valid_cells_new("valid_cells_new");
    indices valid_cells_empty("valid_cells_empty");
    param<double> dist_cells_old ("dist_cells_old");
    double time_start = get_wall_time();
    
    /*variables to define rotation on the positive side (model side)*/
    
    var<> yaw("yaw", yaw_min, yaw_max), pitch("pitch", pitch_min, pitch_max), roll("roll", roll_min, roll_max);
    yaw.in(R(1)); pitch.in(R(1));roll.in(R(1));
    
    func<> r11 = cos(roll)*cos(yaw);r11.eval_all();
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
    var<> tx("tx", tx_min, tx_max), ty("ty", ty_min, ty_max),tz("tz", tz_min, tz_max);
    
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
    vector<var<double>> t;
    
    min_cost_sum=0.0;
    
    bool found_all=true;
    
    int nd=input_data_cloud.size();
    int nm=input_model_cloud.size();
    
    vector<map<double, int>> valid_cells_map(nd);
    
    for(auto i=0;i<nd;i++){
        double min_dist_ij_max=numeric_limits<double>::max();
        double min_dist_ij_min=numeric_limits<double>::max();
        vector<vector<double>> extreme_i;
        /*Feasible region R'_insR(input_data_cloud)*/
        get_extreme_point(extreme_i, point_cloud_data[i], T1, );
        for (int j = 0; j<nm; j++) {
            string key= to_string(i+1)+","+to_string(j+1);
            if(valid_cells_old.size()>=input_data_cloud.size()){
                if(!valid_cells_old.has_key(key)){
                    DebugOff("continued");
                    continue;
                }
            }
            double dist_ij_min, dist_ij_max;
            double dij_voro=std::max(distance_polytopes_gjk(extreme_i, model_voronoi_vertices.at(j))-1e-6, 0.0);
            /*Calling GJK*/
            if(dij_voro<=1e-6){
                vector<vector<double>> extreme_j;
                extreme_j.push_back(point_cloud_model.at(j));
                dist_ij_min=std::max(distance_polytopes_gjk(extreme_i, extreme_j)-1e-6, 0.0);
                if(dist_ij_min<=upper_bound && dist_ij_min<=min_dist_ij_max){
                    dist_ij_max=std::min(max_distance_polytopes(extreme_i, extreme_j), model_voronoi_out_radius[j]);
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
        DebugOff("min_cost_sum "<<min_cost_sum<<endl);
        double vo=valid_cells_old.size();
        if(vo==0){
            vo=input_data_cloud.size()*input_model_cloud.size();
        }
        double vn=valid_cells_new.size();
        double remo=(vo-vn)/vo*100.0;
        DebugOff("valid cells old size "<<vo<<endl);
        DebugOff("valid cells new size "<<vn<<endl);
        DebugOff("rem percen "<<remo<<endl);
        new_cells=valid_cells_new;
    }
    else{
        new_cells=valid_cells_empty;
        DebugOff("No valid cells found "<<endl);
    }
    //return min_cost_sum;
}
#endif

#endif /* BB_h */
