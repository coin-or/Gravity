//
//  BB.h
//  Gravity
//
//  Created by Smitha on 4/2/22.
//

#ifndef BB_h
#define BB_h
#include "gravity/nanoflann.hpp"
#include "Heuristics.h"
void make_box(double xl, double xu, double yl, double yu,double zl, double zu, vector<vector<double>>& box){
    vector<double> res(3,0);
    bool repeat;
    res[0]=xl;
    res[1]=yl;
    res[2]=zl;
    repeat=false;
    for(auto k=0;k<box.size();k++){
        if(abs(box[k][0]-res[0])<=1e-9 && abs(box[k][1]-res[1])<=1e-9 && abs(box[k][2]-res[2])<=1e-9)
            repeat=true;
        DebugOff("repeat "<<box[k][0]<<" "<<box[k][1]<<" "<<box[k][2]<<endl);
    }
    if(!repeat){
        box.push_back(res);
    }
    res[0]=xu;
    res[1]=yl;
    res[2]=zl;
    repeat=false;
    for(auto k=0;k<box.size();k++){
        if(abs(box[k][0]-res[0])<=1e-9 && abs(box[k][1]-res[1])<=1e-9 && abs(box[k][2]-res[2])<=1e-9)
            repeat=true;
        DebugOff("repeat "<<box[k][0]<<" "<<box[k][1]<<" "<<box[k][2]<<endl);
    }
    if(!repeat){
        box.push_back(res);
    }
    res[0]=xl;
    res[1]=yu;
    res[2]=zl;
    repeat=false;
    for(auto k=0;k<box.size();k++){
        if(abs(box[k][0]-res[0])<=1e-9 && abs(box[k][1]-res[1])<=1e-9 && abs(box[k][2]-res[2])<=1e-9)
            repeat=true;
        DebugOff("repeat "<<box[k][0]<<" "<<box[k][1]<<" "<<box[k][2]<<endl);
    }
    if(!repeat){
        box.push_back(res);
    }
    res[0]=xu;
    res[1]=yu;
    res[2]=zl;
    repeat=false;
    for(auto k=0;k<box.size();k++){
        if(abs(box[k][0]-res[0])<=1e-9 && abs(box[k][1]-res[1])<=1e-9 && abs(box[k][2]-res[2])<=1e-9)
            repeat=true;
        DebugOff("repeat "<<box[k][0]<<" "<<box[k][1]<<" "<<box[k][2]<<endl);
    }
    if(!repeat){
        box.push_back(res);
    }
    res[0]=xl;
    res[1]=yl;
    res[2]=zu;
    repeat=false;
    for(auto k=0;k<box.size();k++){
        if(abs(box[k][0]-res[0])<=1e-9 && abs(box[k][1]-res[1])<=1e-9 && abs(box[k][2]-res[2])<=1e-9)
            repeat=true;
        DebugOff("repeat "<<box[k][0]<<" "<<box[k][1]<<" "<<box[k][2]<<endl);
    }
    if(!repeat){
        box.push_back(res);
    }
    res[0]=xu;
    res[1]=yl;
    res[2]=zu;
    repeat=false;
    for(auto k=0;k<box.size();k++){
        if(abs(box[k][0]-res[0])<=1e-9 && abs(box[k][1]-res[1])<=1e-9 && abs(box[k][2]-res[2])<=1e-9)
            repeat=true;
        DebugOff("repeat "<<box[k][0]<<" "<<box[k][1]<<" "<<box[k][2]<<endl);
    }
    if(!repeat){
        box.push_back(res);
    }
    res[0]=xl;
    res[1]=yu;
    res[2]=zu;
    repeat=false;
    for(auto k=0;k<box.size();k++){
        if(abs(box[k][0]-res[0])<=1e-9 && abs(box[k][1]-res[1])<=1e-9 && abs(box[k][2]-res[2])<=1e-9)
            repeat=true;
        DebugOff("repeat "<<box[k][0]<<" "<<box[k][1]<<" "<<box[k][2]<<endl);
    }
    if(!repeat){
        box.push_back(res);
    }
    res[0]=xu;
    res[1]=yu;
    res[2]=zu;
    repeat=false;
    for(auto k=0;k<box.size();k++){
        if(abs(box[k][0]-res[0])<=1e-9 && abs(box[k][1]-res[1])<=1e-9 && abs(box[k][2]-res[2])<=1e-9)
            repeat=true;
        DebugOff("repeat "<<box[k][0]<<" "<<box[k][1]<<" "<<box[k][2]<<endl);
    }
    if(!repeat){
        box.push_back(res);
    }
    
}
void filter_extremei(vector<vector<double>>& ex_i, double xl, double xu, double yl, double yu,double zl, double zu, vector<vector<double>>& ex_ij){
    for(auto k=0;k<8;k++){
        double xl_ik=ex_i[k*8][0];
        double yl_ik=ex_i[k*8][1];
        double zl_ik=ex_i[k*8][2];
        double xu_ik=ex_i[(k+1)*8-1][0];
        double yu_ik=ex_i[(k+1)*8-1][1];
        double zu_ik=ex_i[(k+1)*8-1][2];
        xl_ik=std::max(xl_ik, xl);
        yl_ik=std::max(yl_ik, yl);
        zl_ik=std::max(zl_ik, zl);
        xu_ik=std::min(xu_ik, xu);
        yu_ik=std::min(yu_ik, yu);
        zu_ik=std::min(zu_ik, zu);
        if(xl_ik<=xu_ik && yl_ik<=yu_ik && zl_ik<=zu_ik){
            make_box(xl_ik,xu_ik,yl_ik,yu_ik,zl_ik,zu_ik,ex_ij);
        }
    }
}
void run_preprocess_model_Align(const vector<vector<double>>& point_cloud_model,
                                const vector<vector<double>>& point_cloud_data, const vector<vector<vector<double>>>& model_voronoi_vertices, treenode_p& vec_node_i, int& m_vec_i, double& vec_lb_i, indices& valid_cells_i,param<double>& dist_cost_i, double& prep_time_i, double upper_bound, shared_ptr<Model<double>>& model_i, std::string error_type, vector<double>& ub_i, double roll_min, double roll_max,  double pitch_min,  double pitch_max,  double yaw_min,  double yaw_max, double tx_min, double tx_max, double ty_min, double ty_max, double tz_min,double tz_max,const nanoflann::KDTreeSingleIndexAdaptor<
                                nanoflann::L2_Simple_Adaptor<double, PointCloud<double>>,
                                PointCloud<double>, 3 /* dim */>& index, const map<int, vector<int>>& incomp, const param<double>& dii, const param<double>& djj, const vector<double>& model_voronoi_radius_sq);
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
                    DebugOff(vertex_found_kl<<endl);
                    DebugOff("Failed to find vertex in prep "<<endl);
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
/* Extreme points of the feasible region in which a data point can lie for Registration
 @extreme: Superset of Vertices of extreme points of the feasible region: [Rd+t]; Convex hull of this set is the feasible region.
 @d: Data point
 @theta_vec: vector of variables defining rotation
 @box_t: feasible region for the [T] bounds
 SECANT reverse trainge inequality
 Tangent Traingle inequality
 */
void get_extreme_point(vector<vector<double>>& extreme, const vector<double>& d, const vector<var<double>>& theta_vec, vector<vector<double>>& box_t){
    double box_xlb, box_xub, box_ylb, box_yub, box_zlb, box_zub;
    double tx_min=box_t[0][0],ty_min=box_t[0][1],tz_min=box_t[0][2],tx_max=box_t[7][0],ty_max=box_t[7][1],tz_max=box_t[7][2];
    
    vector<vector<double>> new_vert_a, new_vert_b;
    vector<int> infeas_set_a, infeas_set_b, infeas_set;
    bool vertex_found_a=false, vertex_found_b=false;
    
    double d_mag=pow(d[0],2)+pow(d[1],2)+pow(d[2],2);
    double d_root=sqrt(d_mag);
    double shift_mag_max_root=sqrt(std::max(pow(tx_min,2),pow(tx_max,2))+std::max(pow(ty_min,2),pow(ty_max,2))+std::max(pow(tz_min,2),pow(tz_max,2)));
    double box_max=d_root+shift_mag_max_root;
    
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
    box_xlb=std::max(x_lb+tx_min, box_max*(-1));
    box_xub=std::min(x_ub+tx_max, box_max);
    
    x_range  = get_product_range(x1_bounds, theta21._range);
    y_range  = get_product_range(y1_bounds, theta22._range);
    z_range  = get_product_range(z1_bounds, theta23._range);
    
    rot_y1_bounds->first=std::max(x_range->first + y_range->first + z_range->first, d_root*(-1));
    rot_y1_bounds->second=std::min(x_range->second + y_range->second + z_range->second, d_root);
    
    double y_lb=rot_y1_bounds->first;
    double y_ub=rot_y1_bounds->second;
    box_ylb=std::max(y_lb+ty_min, box_max*(-1));
    box_yub=std::min(y_ub+ty_max, box_max);
    
    
    x_range  = get_product_range(x1_bounds, theta31._range);
    y_range  = get_product_range(y1_bounds, theta32._range);
    z_range  = get_product_range(z1_bounds, theta33._range);
    
    rot_z1_bounds->first=std::max(x_range->first + y_range->first + z_range->first, d_root*(-1));
    rot_z1_bounds->second=std::min(x_range->second + y_range->second + z_range->second, d_root);
    
    double z_lb=rot_z1_bounds->first;
    double z_ub=rot_z1_bounds->second;
    box_zlb=std::max(z_lb+tz_min, box_max*(-1));
    box_zub=std::min(z_ub+tz_max, box_max);
    
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
        auto zttp=sqrt(d_mag-pow(xt,2)-pow(yt,2));
        if( zttp <=z_ub && zttp>=z_lb){
            zt=zttp;
        }
        else{
            zt=zttp*(-1);
        }
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
        DebugOff("Plane Eq not found in extreme data "<<endl);
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
        DebugOff("Vertex A not found "<<endl);
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
        DebugOff("Vertex B not found "<<plane_eq_found<<endl);
    }
    for(auto i=0;i<box_big.size();i++){
        if(std::find (infeas_set.begin(), infeas_set.end(), i) ==infeas_set.end()){
            extreme.push_back(box_big[i]);
        }
    }
    auto extreme_rot=extreme;
    extreme.clear();
    vector<double> res(3);
    for(auto i=0;i<extreme_rot.size();i++){
        for(auto j=0;j<box_t.size();j++){
            res[0]=extreme_rot[i][0]+box_t[j][0];
            res[0]=std::min(std::max(res[0], box_xlb), box_xub);
            res[1]=extreme_rot[i][1]+box_t[j][1];
            res[1]=std::min(std::max(res[1], box_ylb), box_yub);
            res[2]=extreme_rot[i][2]+box_t[j][2];
            res[2]=std::min(std::max(res[2], box_zlb), box_zub);
            extreme.push_back(res);
        }
    }
}
void get_extreme_point(vector<vector<double>>& extreme, const vector<double>& d, const vector<var<double>>& theta_vec, double tx_min, double tx_max, double ty_min, double ty_max, double tz_min, double tz_max)
{
    double d_mag=pow(d[0],2)+pow(d[1],2)+pow(d[2],2);
    double d_root=sqrt(d_mag);
    double shift_mag_max=std::max(pow(tx_min,2),pow(tx_max,2))+std::max(pow(ty_min,2),pow(ty_max,2))+std::max(pow(tz_min,2),pow(tz_max,2));
    double shift_mag_max_root=sqrt(shift_mag_max);
    double shift_mag_min=0.0;
    if(tx_min<=0 && tx_max>=0){
        shift_mag_min+=0;
    }
    else{
        shift_mag_min+=std::min(pow(tx_min,2),pow(tx_max,2));
    }
    if(ty_min<=0 && ty_max>=0){
        shift_mag_min+=0;
    }
    else{
        shift_mag_min+=std::min(pow(ty_min,2),pow(ty_max,2));
    }
    if(tz_min<=0 && tz_max>=0){
        shift_mag_min+=0;
    }
    else{
        shift_mag_min+=std::min(pow(tz_min,2),pow(tz_max,2));
    }
    double shift_mag_min_root=sqrt(shift_mag_min);
    double siq, soq=0;
    if(shift_mag_min_root<=d_root && shift_mag_max_root>=d_root){
        siq=0.0;
    }
    else{
        siq=std::min(pow((shift_mag_min_root-d_root),2),pow((shift_mag_max_root-d_root),2));
    }
    soq=pow((d_root+shift_mag_max_root),2);
    
    vector<vector<double>> new_vert_a, new_vert_b;
    vector<int> infeas_set_a, infeas_set_b, infeas_set;
    bool vertex_found_a=false, vertex_found_b=false;
    
    
    
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
    
    rot_x1_bounds->first=std::max(std::max(x_range->first + y_range->first + z_range->first, d_root*(-1))+tx_min,(d_root+shift_mag_max_root)*(-1));
    rot_x1_bounds->second=std::min(std::min(x_range->second + y_range->second + z_range->second, d_root)+tx_max,(d_root+shift_mag_max_root));
    
    double x_lb=rot_x1_bounds->first;
    double x_ub=rot_x1_bounds->second;
    
    x_range  = get_product_range(x1_bounds, theta21._range);
    y_range  = get_product_range(y1_bounds, theta22._range);
    z_range  = get_product_range(z1_bounds, theta23._range);
    
    rot_y1_bounds->first=std::max(std::max(x_range->first + y_range->first + z_range->first, d_root*(-1))+ty_min,(d_root+shift_mag_max_root)*(-1));
    rot_y1_bounds->second=std::min(std::min(x_range->second + y_range->second + z_range->second, d_root)+ty_max,(d_root+shift_mag_max_root));
    
    double y_lb=rot_y1_bounds->first;
    double y_ub=rot_y1_bounds->second;
    
    
    x_range  = get_product_range(x1_bounds, theta31._range);
    y_range  = get_product_range(y1_bounds, theta32._range);
    z_range  = get_product_range(z1_bounds, theta33._range);
    
    rot_z1_bounds->first=std::max(std::max(x_range->first + y_range->first + z_range->first, d_root*(-1))+tz_min,(d_root+shift_mag_max_root)*(-1));
    rot_z1_bounds->second=std::min(std::min(x_range->second + y_range->second + z_range->second, d_root)+tz_max,(d_root+shift_mag_max_root));
    
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
    secant[3]=siq+x_lb*x_ub+y_lb*y_ub+z_lb*z_ub;
    
    
    vertex_found_a=vertices_box_plane(secant, box_big, new_vert_a, infeas_set_a);
    
    /*Find coordinates of the point to define tangent, Tangent is found parallel to secant, if this fails midpoint is used and checked to ensure does not intersect secant in this range*/
    
    double xt, yt, zt;
    
    vector<double> tangent(4);
    
    xt=(x_lb+x_ub)*0.5;
    yt=(y_lb+y_ub)*0.5;
    auto ztemp=(z_lb+z_ub)*0.5;
    zt=ztemp;
    if(z_lb<=sqrt(soq-xt*xt-yt*yt)<=z_ub){
        zt=sqrt(soq-xt*xt-yt*yt);
    }
    else if(z_lb<=sqrt(soq-xt*xt-yt*yt)*(-1)<=z_ub){
        zt=sqrt(soq-xt*xt-yt*yt)*(-1);
    }
    
    bool plane_eq_found=true;
    if(vertex_found_a){
        for(auto v: new_vert_a){
            auto x=v[0];
            auto y=v[1];
            auto z=v[2];
            if(2*xt*x+2*yt*y+2*zt*z-soq-(xt*xt+yt*yt+zt*zt)>=1e-9){
                plane_eq_found=false;
            }
        }
    }
    
    if(plane_eq_found){
        tangent[0]=2*xt;
        tangent[1]=2*yt;
        tangent[2]=2*zt;
        tangent[3]=-soq-(xt*xt+yt*yt+zt*zt);
        DebugOff("Plane Eq found in extreme data "<<endl);
        vertex_found_b=vertices_box_plane(tangent, box_big, new_vert_b, infeas_set_b);
    }
    else{
        DebugOff("Plane Eq not found in extreme data "<<endl);
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
        DebugOff("Vertex A not found "<<endl);
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
        DebugOff("Vertex B not found "<<plane_eq_found<<endl);
    }
    for(auto i=0;i<box_big.size();i++){
        if(std::find (infeas_set.begin(), infeas_set.end(), i) ==infeas_set.end()){
            extreme.push_back(box_big[i]);
        }
    }
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
#ifdef USE_GJK
extern "C" {
#include "openGJK.h"
}
double distance_polytopes_gjk(const vector<vector<double>>& vec1, const vector<vector<double>>& vec2){
    
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

shared_ptr<Model<double>> Reg_L2_model_rotation_trigonometric(const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max,double tx_min, double tx_max, double ty_min, double ty_max, double tz_min, double tz_max, indices& cells, param<double> dist_cost, const map<int, vector<int>>& incomp)
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
    
    int nd = point_cloud_data.size();
    int nm = point_cloud_model.size();
    DebugOff("nd = " << nd << endl);
    DebugOff("nm = " << nm << endl);
    
    N1 = range(1,nd);
    N2 = range(1,nm);
    
    indices ids = indices("in_x");
    indices idsij = indices("idsij");
    idsij.add_empty_row();
    ids = N2;
    ids.add_empty_row();
    for(auto i=0;i<nd;i++){
        for(auto j=1;j<=nm;j++){
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
    
    //    vector<pair<pair<int,int>,pair<int,int>>> incompatibles;
    //    int count=0;
    //    for(auto i=0;i<nd;i++){
    //        if(incomp.find(i)!=incomp.end()){
    //            int count_i=0;
    //            for (int j = 0; j<nm ; j++) {
    //                string key= to_string(i+1)+","+to_string(j+1);
    //                if(!cells.has_key(key)){
    //                    DebugOff("continued");
    //                    continue;
    //                }
    //                for (int k = 0; k<incomp.at(i).size(); k++) {
    //                    int i1=incomp.at(i)[k];
    //                    string key1= to_string(i1+1)+","+to_string(j+1);
    //                    if(cells.has_key(key1)){
    //                        incompatibles.push_back({{i,j},{i1,j}});
    //                        count++;
    //                        count_i++;
    //                    }
    //                }
    //            }
    //        }
    //    }
    //    if(!incompatibles.empty()){
    //        indices pairs1("pairs1"), pairs2("pairs2");
    //        pairs1 = cells;
    //        pairs2 = cells;
    //        for (const auto &inc_pair : incompatibles) {
    //            string key1 = to_string(inc_pair.first.first+1)+","+to_string(inc_pair.first.second+1);
    //            string key2 = to_string(inc_pair.second.first+1)+","+to_string(inc_pair.second.second+1);
    //            if(cells.has_key(key1) && cells.has_key(key2)){
    //                pairs1.add_ref(key1);
    //                pairs2.add_ref(key2);
    //            }
    //        }
    //
    //        if(pairs1.is_indexed()){
    //            DebugOn("Number of incompatible pairs constraints = " << pairs1.size() << endl);
    //            Constraint<> incomp_pairs("incomp_pairs");
    //            incomp_pairs = bin.in(pairs1) + bin.in(pairs2);
    //            Reg->add(incomp_pairs.in(range(1,pairs1.size()))<=1);
    //            //incomp_pairs.print();
    //        }
    //    }
    
    
    Reg->min(sum(deltax) + sum(deltay)+sum(deltaz));
    
    return Reg;
}
shared_ptr<Model<double>> Reg_L2_model_rotation_trigonometric_small(const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max,double tx_min, double tx_max, double ty_min, double ty_max, double tz_min, double tz_max, indices& cells, param<double> dist_cost, bool add_nc)
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
    
    int n1 = point_cloud_data.size();
    int n2 = point_cloud_model.size();
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
    Reg->add(new_xm.in(N1), new_ym.in(N1), new_zm.in(N1));
    
    auto ids1 = theta11.repeat_id(N1.size());
    
    Constraint<> Def_newxm("Def_newxm");
    Def_newxm = new_xm+x1*theta11.in(ids1)+y1*theta12.in(ids1)+z1*theta13.in(ids1)+tx.in(ids1)-product(x2.in(ids),bin.in_matrix(1, 1));
    Reg->add(Def_newxm.in(N1)==0);
    
    
    Constraint<> Def_newym("Def_newym");
    Def_newym = new_ym+x1*theta21.in(ids1)+y1*theta22.in(ids1)+z1*theta23.in(ids1)+ty.in(ids1)-product(y2.in(ids),bin.in_matrix(1, 1));
    Reg->add(Def_newym.in(N1)==0);
    
    Constraint<> Def_newzm("Def_newzm");
    Def_newzm = new_zm+x1*theta31.in(ids1)+y1*theta32.in(ids1)+z1*theta33.in(ids1)+tz.in(ids1)-product(z2.in(ids),bin.in_matrix(1, 1));
    Reg->add(Def_newzm.in(N1)==0);
    
    Constraint<> OneBin("OneBin");
    OneBin = bin.in_matrix(1, 1);
    Reg->add(OneBin.in(N1)==1);
    
    var<> deltax("deltax"), deltay("deltay"), deltaz("deltaz");
    Reg->add(deltax.in(N1));
    Reg->add(deltay.in(N1));
    Reg->add(deltaz.in(N1));
    
    Constraint<> Def_deltax("Def_deltax");
    Def_deltax=pow(new_xm, 2)-deltax;
    Reg->add(Def_deltax.in(N1)<=0);
    
    Constraint<> Def_deltay("Def_deltay");
    Def_deltay=pow(new_ym, 2)-deltay;
    Reg->add(Def_deltay.in(N1)<=0);
    
    Constraint<> Def_deltaz("Def_deltaz");
    Def_deltaz=pow(new_zm, 2)-deltaz;
    Reg->add(Def_deltaz.in(N1)<=0);
    
    
    if(dist_cost._indices->_keys->size()!=0){
        Constraint<> delta_cost("delta_cost");
        delta_cost=product(dist_cost.in(idsij), bin.in_matrix(1,1))-deltax-deltay-deltaz;
        Reg->add(delta_cost.in(N1)<=0);
    }
    
    //    func<> d1_f = 1-theta11-theta22+theta33;
    //    func<> d2_f = 1+theta11-theta22-theta33;;
    //    func<> d3_f = 1+theta11+theta22+theta33;
    //    func<> d4_f = 1-theta11+theta22-theta33;
    //
    //    var<> d1("d1", 0, d1_f._range->second);
    //    var<> d2("d2", 0, d2_f._range->second);
    //    var<> d3("d3", 0, d3_f._range->second);
    //    var<> d4("d4", 0, d4_f._range->second);
    //    Reg->add(d1.in(R(1)));Reg->add(d2.in(R(1)));Reg->add(d3.in(R(1)));Reg->add(d4.in(R(1)));
    //    var<> p1("p1");
    //    var<> p2("p2");
    //    var<> p3("p3");
    //    var<> p4("p4");
    //    var<> p5("p5");
    //    var<> p6("p6");
    //    Reg->add(p1.in(R(1)));Reg->add(p2.in(R(1)));Reg->add(p3.in(R(1)));Reg->add(p4.in(R(1)));Reg->add(p5.in(R(1)));
    //    Reg->add(p6.in(R(1)));
    
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
    
    //    Constraint<> p_1("p_1");
    //    p_1=theta13+theta31-p1;
    //    Reg->add(p_1.in(range(0,0))==0);
    //
    //    Constraint<> p_2("p_2");
    //    p_2=theta12-theta21-p2;
    //    Reg->add(p_2.in(range(0,0))==0);
    //
    //    Constraint<> p_3("p_3");
    //    p_3=theta23+theta32-p3;
    //    Reg->add(p_3.in(range(0,0))==0);
    //
    //    Constraint<> p_4("p_4");
    //    p_4=theta23-theta32-p4;
    //    Reg->add(p_4.in(range(0,0))==0);
    //
    //    Constraint<> p_5("p_5");
    //    p_5=theta12+theta21-p5;
    //    Reg->add(p_5.in(range(0,0))==0);
    //
    //    Constraint<> p_6("p_6");
    //    p_6=theta31-theta13-p6;
    //    Reg->add(p_6.in(range(0,0))==0);
    //
    //
    //
    //    Constraint<> soc_12("soc_12");
    //    soc_12 = pow(p1,2)-d1*d2;
    //    Reg->add(soc_12.in(range(0,0))<=0);
    //
    //    Constraint<> soc_13("soc_13");
    //    soc_13 = pow(p2,2)-d1*d3;
    //    Reg->add(soc_13.in(range(0,0))<=0);
    //
    //    Constraint<> soc_14("soc_14");
    //    soc_14 = pow(p3,2)-d1*d4;
    //    Reg->add(soc_14.in(range(0,0))<=0);
    //
    //    Constraint<> soc_23("soc_23");
    //    soc_23 = pow(p4,2)-d2*d3;
    //    Reg->add(soc_23.in(range(0,0))<=0);
    //
    //    Constraint<> soc_24("soc_24");
    //    soc_24 = pow(p5,2)-d2*d4;
    //    Reg->add(soc_24.in(range(0,0))<=0);
    //
    //    Constraint<> soc_34("soc_34");
    //    soc_34 = pow(p6,2)-d3*d4;
    //    Reg->add(soc_34.in(range(0,0))<=0);
    
    if(add_nc){
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
    }
    
    //    Constraint<> Trow1("Trow1");
    //    Trow1=pow(theta11,2)+pow(theta12,2)+pow(theta13,2);
    //    Reg->add(Trow1.in(range(0,0))<=1);
    //
    //    Constraint<> Trow2("Trow2");
    //    Trow2=pow(theta21,2)+pow(theta22,2)+pow(theta23,2);
    //    Reg->add(Trow2.in(range(0,0))<=1);
    //
    //    Constraint<> Trow3("Trow3");
    //    Trow3=pow(theta31,2)+pow(theta32,2)+pow(theta33,2);
    //    Reg->add(Trow3.in(range(0,0))<=1);
    //
    //    Constraint<> Tcol1("Tcol1");
    //    Tcol1=pow(theta11,2)+pow(theta21,2)+pow(theta31,2);
    //    Reg->add(Tcol1.in(range(0,0))<=1);
    //
    //    Constraint<> Tcol2("Tcol2");
    //    Tcol2=pow(theta12,2)+pow(theta22,2)+pow(theta32,2);
    //    Reg->add(Tcol2.in(range(0,0))<=1);
    //
    //    Constraint<> Tcol3("Tcol3");
    //    Tcol3=pow(theta13,2)+pow(theta23,2)+pow(theta33,2);
    //    Reg->add(Tcol3.in(range(0,0))<=1);
    
    
    Constraint<> sec_row1("sec_row1");
    sec_row1=1+theta11.get_lb()*theta11.get_ub()+theta12.get_lb()*theta12.get_ub()+theta13.get_lb()*theta13.get_ub()-(theta11.get_lb()+theta11.get_ub())*theta11-(theta12.get_lb()+theta12.get_ub())*theta12-(theta13.get_lb()+theta13.get_ub())*theta13;
    Reg->add(sec_row1.in(range(0,0))<=0);
    
    Constraint<> sec_row2("sec_row2");
    sec_row2=1+theta21.get_lb()*theta21.get_ub()+theta22.get_lb()*theta22.get_ub()+theta23.get_lb()*theta23.get_ub()-(theta21.get_lb()+theta21.get_ub())*theta21-(theta22.get_lb()+theta22.get_ub())*theta22-(theta23.get_lb()+theta23.get_ub())*theta23;
    Reg->add(sec_row2.in(range(0,0))<=0);
    
    Constraint<> sec_row3("sec_row3");
    sec_row3=1+theta31.get_lb()*theta31.get_ub()+theta32.get_lb()*theta32.get_ub()+theta33.get_lb()*theta33.get_ub()-(theta31.get_lb()+theta31.get_ub())*theta31-(theta32.get_lb()+theta32.get_ub())*theta32-(theta33.get_lb()+theta33.get_ub())*theta33;
    Reg->add(sec_row3.in(range(0,0))<=0);
    
    
    Constraint<> sec_col1("sec_col1");
    sec_col1=1+theta11.get_lb()*theta11.get_ub()+theta21.get_lb()*theta21.get_ub()+theta31.get_lb()*theta31.get_ub()-(theta11.get_lb()+theta11.get_ub())*theta11-(theta21.get_lb()+theta21.get_ub())*theta21-(theta31.get_lb()+theta31.get_ub())*theta31;
    Reg->add(sec_col1.in(range(0,0))<=0);
    
    Constraint<> sec_col2("sec_col2");
    sec_col2=1+theta12.get_lb()*theta12.get_ub()+theta22.get_lb()*theta22.get_ub()+theta32.get_lb()*theta32.get_ub()-(theta12.get_lb()+theta12.get_ub())*theta12-(theta22.get_lb()+theta22.get_ub())*theta22-(theta32.get_lb()+theta32.get_ub())*theta32;
    Reg->add(sec_col2.in(range(0,0))<=0);
    
    Constraint<> sec_col3("sec_col3");
    sec_col3=1+theta13.get_lb()*theta13.get_ub()+theta23.get_lb()*theta23.get_ub()+theta33.get_lb()*theta33.get_ub()-(theta13.get_lb()+theta13.get_ub())*theta13-(theta23.get_lb()+theta23.get_ub())*theta23-(theta33.get_lb()+theta33.get_ub())*theta33;
    Reg->add(sec_col3.in(range(0,0))<=0);
    
    
    
    
    Reg->min(sum(deltax) + sum(deltay)+sum(deltaz));
    
    return Reg;
}
shared_ptr<Model<double>> Reg_L2_model_rotation_trigonometric_ve(const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max,double tx_min, double tx_max, double ty_min, double ty_max, double tz_min, double tz_max, indices& cells, param<double> dist_cost, bool add_nc)
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
    
    int n1 = point_cloud_data.size();
    int n2 = point_cloud_model.size();
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
    
    param<> rotx_min("rotx_min");
    param<> rotx_max("rotx_max");
    param<> roty_min("roty_min");
    param<> roty_max("roty_max");
    param<> rotz_min("rotz_min");
    param<> rotz_max("rotz_max");
    param<> d_mag("d_mag");
    shared_ptr<pair<double,double>> x1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> y1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> z1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> new_x1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> new_y1_bounds = make_shared<pair<double,double>>();
    shared_ptr<pair<double,double>> new_z1_bounds = make_shared<pair<double,double>>();
    for(auto i=0;i<n1;i++){
        auto d_sq=(pow(point_cloud_data.at(i)[0],2)+pow(point_cloud_data.at(i)[1],2)+pow(point_cloud_data.at(i)[2],2));
        auto d_root=sqrt(pow(point_cloud_data.at(i)[0],2)+pow(point_cloud_data.at(i)[1],2)+pow(point_cloud_data.at(i)[2],2));
        x1_bounds->first = point_cloud_data.at(i)[0];
        x1_bounds->second = point_cloud_data.at(i)[0];
        y1_bounds->first = point_cloud_data.at(i)[1];
        y1_bounds->second = point_cloud_data.at(i)[1];
        z1_bounds->first = point_cloud_data.at(i)[2];
        z1_bounds->second = point_cloud_data.at(i)[2];
        auto x_range  = get_product_range(x1_bounds, theta11._range);
        auto y_range  = get_product_range(y1_bounds, theta12._range);
        auto z_range  = get_product_range(z1_bounds, theta13._range);
        *new_x1_bounds = {std::max(x_range->first + y_range->first + z_range->first, d_root*(-1)),std::min(x_range->second + y_range->second + z_range->second, d_root)};
        rotx_min.add_val(to_string(i+1), new_x1_bounds->first);
        rotx_max.add_val(to_string(i+1), new_x1_bounds->second);
        
        x_range  = get_product_range(x1_bounds, theta21._range);
        y_range  = get_product_range(y1_bounds, theta22._range);
        z_range  = get_product_range(z1_bounds, theta23._range);
        *new_y1_bounds = {std::max(x_range->first + y_range->first + z_range->first, d_root*(-1)), std::min(x_range->second + y_range->second + z_range->second, d_root)};
        roty_min.add_val(to_string(i+1), new_y1_bounds->first);
        roty_max.add_val(to_string(i+1), new_y1_bounds->second);
        
        x_range  = get_product_range(x1_bounds, theta31._range);
        y_range  = get_product_range(y1_bounds, theta32._range);
        z_range  = get_product_range(z1_bounds, theta33._range);
        *new_z1_bounds = {std::max(x_range->first + y_range->first + z_range->first, d_root*(-1)), std::min(x_range->second + y_range->second + z_range->second, d_root)};
        rotz_min.add_val(to_string(i+1), new_z1_bounds->first);
        rotz_max.add_val(to_string(i+1), new_z1_bounds->second);
        d_mag.add_val(to_string(i+1), d_sq);
    }
    
    
    var<> new_xm("new_xm"), new_ym("new_ym"), new_zm("new_zm");
    Reg->add(new_xm.in(N1), new_ym.in(N1), new_zm.in(N1));
    
    
    var<> rotx("rotx", rotx_min, rotx_max), roty("roty",roty_min, roty_max), rotz("rotz",rotz_min, rotz_max);
    Reg->add(rotx.in(N1), roty.in(N1), rotz.in(N1));
    
    auto ids1 = theta11.repeat_id(N1.size());
    
    Constraint<> Def_newxm("Def_newxm");
    Def_newxm = new_xm-product(x2.in(ids),bin.in_matrix(1, 1))+rotx+tx.in(ids1);
    Reg->add(Def_newxm.in(N1)==0);
    
    
    Constraint<> Def_newym("Def_newym");
    Def_newym = new_ym-product(y2.in(ids),bin.in_matrix(1, 1))+roty+ty.in(ids1);
    Reg->add(Def_newym.in(N1)==0);
    
    Constraint<> Def_newzm("Def_newzm");
    Def_newzm = new_zm-product(z2.in(ids),bin.in_matrix(1, 1))+rotz+tz.in(ids1);
    Reg->add(Def_newzm.in(N1)==0);
    
    
    Constraint<> xd("xd");
    xd=x1*theta11.in(ids1)+y1*theta12.in(ids1)+z1*theta13.in(ids1)-rotx;
    Reg->add(xd.in(N1)==0);
    
    Constraint<> yd("yd");
    yd=x1*theta21.in(ids1)+y1*theta22.in(ids1)+z1*theta23.in(ids1)-roty;
    Reg->add(yd.in(N1)==0);
    
    
    Constraint<> zd("zd");
    zd=x1*theta31.in(ids1)+y1*theta32.in(ids1)+z1*theta33.in(ids1)-rotz;
    Reg->add(zd.in(N1)==0);
    
    Constraint<> sec_rot("sec_rot");
    sec_rot=d_mag+rotx_min*rotx_max+roty_min*roty_max+rotz_min*rotz_max-(rotx_min+rotx_max)*rotx-(roty_min+roty_max)*roty-(rotz_min+rotz_max)*rotz;
    Reg->add(sec_rot.in(N1)<=0);
    
    Constraint<> OneBin("OneBin");
    OneBin = bin.in_matrix(1, 1);
    Reg->add(OneBin.in(N1)==1);
    
    var<> deltax("deltax"), deltay("deltay"), deltaz("deltaz");
    Reg->add(deltax.in(N1));
    Reg->add(deltay.in(N1));
    Reg->add(deltaz.in(N1));
    
    Constraint<> Def_deltax("Def_deltax");
    Def_deltax=pow(new_xm, 2)-deltax;
    Reg->add(Def_deltax.in(N1)<=0);
    
    Constraint<> Def_deltay("Def_deltay");
    Def_deltay=pow(new_ym, 2)-deltay;
    Reg->add(Def_deltay.in(N1)<=0);
    
    Constraint<> Def_deltaz("Def_deltaz");
    Def_deltaz=pow(new_zm, 2)-deltaz;
    Reg->add(Def_deltaz.in(N1)<=0);
    
    
    if(dist_cost._indices->_keys->size()!=0){
        Constraint<> delta_cost("delta_cost");
        delta_cost=product(dist_cost.in(idsij), bin.in_matrix(1,1))-deltax-deltay-deltaz;
        Reg->add(delta_cost.in(N1)<=0);
    }
    if(add_nc){
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
    }
    
    Reg->min(sum(deltax) + sum(deltay)+sum(deltaz));
    
    return Reg;
}
#ifdef USE_GJK
void preprocess_lid(const vector<vector<double>>& point_cloud_model, const vector<vector<double>>& point_cloud_data, const vector<vector<vector<double>>>& model_voronoi_vertices, indices& valid_cells_old, indices& new_cells, param<double>& dist_cells, double roll_min, double roll_max, double pitch_min, double pitch_max, double yaw_min, double yaw_max, double& tx_min, double& tx_max, double& ty_min, double& ty_max, double& tz_min, double& tz_max, double upper_bound, double& prep_time, double& min_cost_sum, string error_type, const param<double>& dist_ii, const param<double>& dist_jj, const vector<double>& model_voronoi_radius_sq)
{
    int nd=point_cloud_data.size();
    int nm=point_cloud_model.size();
    
    multimap<double,int,greater<double>> map_nd;
    
    double ub_sq_root=sqrt(upper_bound);
    
    min_cost_sum=0.0;
    bool found_all=true;
    
    vector<vector<double>> box_t;
    vector<double> coord(3,0);
    coord[0]=tx_min;
    coord[1]=ty_min;
    coord[2]=tz_min;
    box_t.push_back(coord);
    coord[0]=tx_min;
    coord[1]=ty_max;
    coord[2]=tz_min;
    box_t.push_back(coord);
    coord[0]=tx_max;
    coord[1]=ty_min;
    coord[2]=tz_min;
    box_t.push_back(coord);
    coord[0]=tx_max;
    coord[1]=ty_max;
    coord[2]=tz_min;
    box_t.push_back(coord);
    coord[0]=tx_min;
    coord[1]=ty_min;
    coord[2]=tz_max;
    box_t.push_back(coord);
    coord[0]=tx_min;
    coord[1]=ty_max;
    coord[2]=tz_max;
    box_t.push_back(coord);
    coord[0]=tx_max;
    coord[1]=ty_min;
    coord[2]=tz_max;
    box_t.push_back(coord);
    coord[0]=tx_max;
    coord[1]=ty_max;
    coord[2]=tz_max;
    box_t.push_back(coord);
    
    vector<vector<double>> box_j;
    
    vector<vector<int>> nd_vec(nd);
    vector<vector<int>> nm_vec(nm);
    
    prep_time=0;
    indices valid_cells("valid_cells");
    indices valid_cells_new("valid_cells_new");
    indices valid_cells_empty("valid_cells_empty");
    param<double> dist_cells_old("dist_cells_old");
    param<double> dist_cells_max("dist_cells_max");
    param<double> dist_cells_max_novoro("dist_cells_max_novoro");
    double time_start = get_wall_time();
    vector<map<double, int>> valid_cells_map(nd);
    
    
    
    map<int, bool> new_model_pts;
    for(auto i=0;i<nd;i++){
        for (int j = 0; j<nm; j++) {
            string key= to_string(i+1)+","+to_string(j+1);
            if(valid_cells_old.size()>=point_cloud_data.size()){
                if(!valid_cells_old.has_key(key)){
                    continue;
                }
            }
            new_model_pts.insert ( std::pair<int,bool>(j,true) );
        }
    }
    
    for (auto it=new_model_pts.begin();it!=new_model_pts.end();it++) {
        box_j.push_back(point_cloud_model.at(it->first));
    }
    auto d_boxt_boxj=std::max(distance_polytopes_gjk(box_t, box_j)-1e-6, 0.0);
    DebugOff("d_boxt_boxj "<<d_boxt_boxj<<endl);
    if(d_boxt_boxj>=1e-6){
        DebugOn("intersection test failed "<<d_boxt_boxj<<endl);
        found_all=false;
    }
    else{
        found_all=true;
    }
    
    
    //triangle inequality?
    
    
    /*variables to define rotation on the positive side (model side)*/
    
    if(found_all){
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
        for(auto i=0;i<nd;i++){
            double min_dist_ij_max=numeric_limits<double>::max();
            double min_dist_ij_min=numeric_limits<double>::max();
            vector<vector<double>> extreme_i;
            get_extreme_point(extreme_i, point_cloud_data[i], T1, box_t);
            // get_extreme_point(extreme_i, point_cloud_data[i], T1, tx_min, tx_max, ty_min, ty_max,tz_min, tz_max);
            for (int j = 0; j<nm; j++) {
                string key= to_string(i+1)+","+to_string(j+1);
                if(valid_cells_old.size()>=point_cloud_data.size()){
                    if(!valid_cells_old.has_key(key)){
                        DebugOff("continued");
                        continue;
                    }
                }
                double dist_ij_min, dist_ij_max;
                double dij_voro=std::max(distance_polytopes_gjk(extreme_i, model_voronoi_vertices.at(j)), 0.0);
                /*Calling GJK*/
                if(dij_voro<=1e-9){
                    vector<vector<double>> extreme_j;
                    extreme_j.push_back(point_cloud_model.at(j));
                    dist_ij_min=std::max(distance_polytopes_gjk(extreme_i, extreme_j)-1e-9, 0.0);
                    if(dist_ij_min<=upper_bound){
                        dist_ij_max=max_distance_polytopes(extreme_i, extreme_j);
                        dist_cells_max_novoro.add_val(key, dist_ij_max);
                        dist_ij_max=std::min(dist_ij_max, model_voronoi_radius_sq[j]);
                        if(dist_ij_min>=dist_ij_max){
                            DebugOn("crossover"<<" "<<dist_ij_min<<" "<<dist_ij_max<<endl);
                        }
                        else{
                            if(dist_ij_min<=min_dist_ij_max){
                                if(dist_ij_min<=min_dist_ij_min){
                                    min_dist_ij_min=dist_ij_min;
                                }
                                if(dist_ij_max<=min_dist_ij_max){
                                    min_dist_ij_max=dist_ij_max;
                                }
                                valid_cells_map[i].insert(pair<double, int>(dist_ij_max, j));
                                dist_cells_old.add_val(key, dist_ij_min);
                                dist_cells_max.add_val(key, dist_ij_max);
                                valid_cells.insert(key);
                                nd_vec[i].push_back(j);
                                nm_vec[j].push_back(i);
                            }
                        }
                    }
                    else{
                        DebugOff("false "<<dist_ij_min<<endl);
                    }
                }
                
            }
            min_cost_sum+=min_dist_ij_min;
            map_nd.insert(std::pair<double,int>(min_dist_ij_min,i));
            if(valid_cells_map[i].empty() || min_cost_sum>=upper_bound+1e-6){
                found_all=false;
                break;
            }
        }
    }
    if(found_all){
        for(auto i=0;i<nd;i++){
            if(nd_vec[i].size()>=2){
                for(auto l=0;l<nd_vec[i].size()-1;l++){
                    auto j=nd_vec[i][l];
                    auto key_j=to_string(i+1)+","+to_string(j+1);
                    auto dij_min_sq=sqrt(dist_cells_old.eval(key_j));
                    auto dij_max_sq=sqrt(dist_cells_max_novoro.eval(key_j));
                    for(auto m=l+1;m<nd_vec[i].size();m++){
                        auto k=nd_vec[i][m];
                        auto key_k=to_string(i+1)+","+to_string(k+1);
                        auto dik_min_sq=sqrt(dist_cells_old.eval(key_k));
                        auto dik_max_sq=sqrt(dist_cells_max_novoro.eval(key_k));
                        auto djk=dist_jj.eval(to_string(j+1)+","+to_string(k+1));
                        if(djk<=dik_min_sq || djk>=dik_max_sq ){
                            auto temp=std::max(djk-dik_max_sq,dik_min_sq-djk);
                            if(temp>=dij_min_sq){
                                dij_min_sq=temp;
                                dist_cells_old.set_val(key_j, temp*temp);
                            }
                        }
                        if(djk<=dij_min_sq || djk>=dij_max_sq ){
                            auto temp=std::max(djk-dij_max_sq,dij_min_sq-djk);
                            if(temp>=dik_min_sq){
                                dik_min_sq=temp;
                                dist_cells_old.set_val(key_k, temp*temp);
                            }
                        }
                    }
                }
            }
        }
        for(auto j=0;j<nm;j++){
            if(nm_vec[j].size()>=2){
                for(auto l=0;l<nm_vec[j].size()-1;l++){
                    auto i=nm_vec[j][l];
                    auto key_i=to_string(i+1)+","+to_string(j+1);
                    auto dij_min_sq=sqrt(dist_cells_old.eval(key_i));
                    auto dij_max_sq=sqrt(dist_cells_max_novoro.eval(key_i));
                    for(auto m=l+1;m<nm_vec[j].size();m++){
                        auto k=nm_vec[j][m];
                        auto key_k=to_string(k+1)+","+to_string(j+1);
                        auto dkj_min_sq=sqrt(dist_cells_old.eval(key_k));
                        auto dkj_max_sq=sqrt(dist_cells_max_novoro.eval(key_k));
                        auto dik=dist_ii.eval(to_string(i+1)+","+to_string(k+1));
                        if(dik<=dkj_min_sq || dik>=dkj_max_sq ){
                            auto temp=std::max(dik-dkj_max_sq,dkj_min_sq-dik);
                            if(temp>=dij_min_sq){
                                dij_min_sq=temp;
                                dist_cells_old.set_val(key_i, temp*temp);
                            }
                        }
                        if(dik<=dij_min_sq || dik>=dij_max_sq ){
                            auto temp=std::max(dik-dij_max_sq,dij_min_sq-dik);
                            if(temp>=dkj_min_sq){
                                dkj_min_sq=temp;
                                dist_cells_old.set_val(key_k, temp*temp);
                            }
                        }
                    }
                }
            }
        }
    }
    double min_cost_nd_sum=0;
    /*Looping again to ensure all valid cells have min_dist less than min_dist_ij_max*/
    if(found_all){
        double new_tx_min=0, new_ty_min=0,new_tz_min=0,new_tx_max=0,new_ty_max=0,new_tz_max=0;
        for(auto itn=map_nd.begin();itn!=map_nd.end();itn++)
        {
            auto i=itn->second;
            min_cost_nd_sum+=itn->first;
            double xm_min=9999.0, ym_min=9999.0, zm_min=9999.0;
            double xm_max=-9999.0, ym_max=-9999.0, zm_max=-9999.0;
            auto it=valid_cells_map[i].begin();
            double min_dist_ij_max=it->first;
            found_all=false;
            for (int j = 0; j<nm; j++) {
                string key= to_string(i+1)+","+to_string(j+1);
                if(!valid_cells.has_key(key)){
                    DebugOff("continued");
                    continue;
                }
                auto dij=dist_cells_old.eval(key);
                if(dij<=min_dist_ij_max && dij<=upper_bound-min_cost_nd_sum+1e-6){
                    found_all=true;
                }
                else{
                    dist_cells_old.set_val(key,min_dist_ij_max*10);
                }
                auto xm= point_cloud_model.at(j)[0];
                auto ym= point_cloud_model.at(j)[1];
                auto zm= point_cloud_model.at(j)[2];
                if(xm<=xm_min){
                    xm_min=xm;
                }
                if(xm>=xm_max){
                    xm_max=xm;
                }
                if(ym<=ym_min){
                    ym_min=ym;
                }
                if(ym>=ym_max){
                    ym_max=ym;
                }
                if(zm<=zm_min){
                    zm_min=zm;
                }
                if(zm>=zm_max){
                    zm_max=zm;
                }
            }
            if(!found_all){
                break;
            }
            new_tx_min+=xm_min;
            new_tx_max+=xm_max;
            new_ty_min+=ym_min;
            new_ty_max+=ym_max;
            new_tz_min+=zm_min;
            new_tz_max+=zm_max;
        }
        new_tx_min/=nd;
        new_tx_max/=nd;
        new_ty_min/=nd;
        new_ty_max/=nd;
        new_tz_min/=nd;
        new_tz_max/=nd;
        tx_min=std::max(new_tx_min-1e-6, tx_min);
        tx_max=std::min(new_tx_max+1e-6, tx_max);
        ty_min=std::max(new_ty_min-1e-6, ty_min);
        ty_max=std::min(new_ty_max+1e-6, ty_max);
        tz_min=std::max(new_tz_min-1e-6, tz_min);
        tz_max=std::min(new_tz_max+1e-6, tz_max);
        DebugOff("tx "<<tx_min<<" "<<tx_max<<endl);
        DebugOff("ty "<<ty_min<<" "<<ty_max<<endl);
        DebugOff("tz "<<tz_min<<" "<<tz_max<<endl);
        if(tx_min>=tx_max+1e-9){
            DebugOn("tx cross "<<tx_min<<" "<<tx_max<<endl);
            found_all=false;
        }
        if(ty_min>=ty_max+1e-9){
            DebugOn("ty cross "<<ty_min<<" "<<ty_max<<endl);
            found_all=false;
        }
        if(tz_min>=tz_max+1e-9){
            DebugOn("tz cross "<<tz_min<<" "<<tz_max<<endl);
            found_all=false;
        }
        if(tx_max<=tx_min-1e-9){
            DebugOn("tx cross "<<tx_min<<" "<<tx_max<<endl);
            found_all=false;
        }
        if(ty_max<=ty_min-1e-9){
            DebugOn("ty cross "<<ty_min<<" "<<ty_max<<endl);
            found_all=false;
        }
        if(tz_max<=tz_min-1e-9){
            DebugOn("tz cross "<<tz_min<<" "<<tz_max<<endl);
            found_all=false;
        }
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
        }
        DebugOn("min_cost_sum "<<min_cost_sum<<endl);
        double vo=valid_cells_old.size();
        if(vo==0){
            vo=point_cloud_data.size()*point_cloud_model.size();
        }
        double vn=valid_cells_new.size();
        double remo=(vo-vn)/vo*100.0;
        DebugOn("valid cells old size "<<vo<<endl);
        DebugOn("valid cells new size "<<vn<<endl);
        DebugOn("rem percen "<<remo<<endl);
        new_cells=valid_cells_new;
    }
    if(!found_all)
    {
        new_cells=valid_cells_empty;
        DebugOn("No valid cells found "<<endl);
    }
    prep_time=get_wall_time()-time_start;
    DebugOff("prep time "<<prep_time<<" "<<roll_min<<" "<<roll_max<<" "<<pitch_min<<" "<<pitch_max<<" "<<yaw_min<<" "<<yaw_max<<" "<<tx_min<<" "<<tx_max<<" "<<ty_min<<" "<<ty_max<<" "<<tz_min<<" "<<tz_max<<endl);
    //return min_cost_sum;
}

#endif
void run_preprocess_only_parallel(const vector<vector<double>>& point_cloud_model,
                                  const vector<vector<double>>& point_cloud_data, const vector<vector<vector<double>>>& model_voronoi_vertices, vector<treenode_p>& vec_node, vector<double>& vec_lb, vector<indices>& valid_cells, int nb_threads, double upper_bound, double lower_bound, int iter, string error_type){
    //size_t nb_threads = std::thread::hardware_concurrency();
    std::vector<thread> threads;
    int nd=point_cloud_data.size();
    int num=vec_node.size();
    if(num==0){
        DebugOff("in run_parallel(models...), models is empty, returning");
    }
    
    vector<param<double>> vec_dist_cost;
    for(auto i=0;i<num;i++){
        param<double> dist_cost("dist_cost");
        vec_dist_cost.push_back(dist_cost);
    }
    
    valid_cells.resize(num);
    vec_lb.resize(num, 0.0);
    
    vector<double> vec_prep_time;
    vec_prep_time.resize(num, 0.0);
    int npass=num/nb_threads+1;
    DebugOn("npass num "<<npass<<" "<<num<<endl);
    for (auto j = 0; j < npass; j++) {
        
        for (auto i = j*nb_threads; i < std::min((j+1)*nb_threads, num); i++) {
#ifdef USE_GJK
            //threads.push_back(thread(&preprocess_lid, ref(point_cloud_model), ref(point_cloud_data), ref(model_voronoi_vertices),  ref(vec_node[i].valid_cells), ref(valid_cells[i]),  ref(vec_dist_cost[i]), vec_node[i].roll.first, vec_node[i].roll.second, vec_node[i].pitch.first, vec_node[i].pitch.second, vec_node[i].yaw.first ,vec_node[i].yaw.second, ref(vec_node[i].tx.first), ref(vec_node[i].tx.second), ref(vec_node[i].ty.first), ref(vec_node[i].ty.second), ref(vec_node[i].tz.first) ,ref(vec_node[i].tz.second), upper_bound, ref(vec_prep_time[i]), ref(vec_lb[i]), error_type));
#endif
        }
        
        for(auto &t : threads){
            t.join();
        }
        threads.clear();
        DebugOn("one pass "<<j<<" "<<num<<endl);
    }
    for (auto i = 0; i < num; i++) {
        if(valid_cells[i].size()>=point_cloud_data.size()){
            vec_node[i].lb=vec_lb[i];
            vec_node[i].valid_cells=valid_cells[i];
            vec_node[i].dist_cost_cells=vec_dist_cost[i];
        }
        else{
            vec_node[i].lb=1000*upper_bound;
        }
    }
}
void run_preprocess_model_Align(const vector<vector<double>>& point_cloud_model,
                                const vector<vector<double>>& point_cloud_data, const vector<vector<vector<double>>>& model_voronoi_vertices, treenode_p& vec_node_i, int& m_vec_i, double& vec_lb_i, indices& valid_cells_i,param<double>& dist_cost_i, double& prep_time_i, double upper_bound, shared_ptr<Model<double>>& model_i, std::string error_type, vector<double>& ub_i, double roll_min, double roll_max,  double pitch_min,  double pitch_max,  double yaw_min,  double yaw_max, double tx_min, double tx_max, double ty_min, double ty_max, double tz_min,double tz_max,const nanoflann::KDTreeSingleIndexAdaptor<
                                nanoflann::L2_Simple_Adaptor<double, PointCloud<double>>,
                                PointCloud<double>, 3 /* dim */>& index, const map<int, vector<int>>& incomp, const param<double>& dii, const param<double>& djj, const vector<double>& model_voronoi_radius_sq){
    icp_new(index, point_cloud_model, point_cloud_data,vec_node_i.roll.first, vec_node_i.roll.second, vec_node_i.pitch.first, vec_node_i.pitch.second, vec_node_i.yaw.first ,vec_node_i.yaw.second, vec_node_i.tx.first,vec_node_i.tx.second, vec_node_i.ty.first,vec_node_i.ty.second, vec_node_i.tz.first,vec_node_i.tz.second, roll_min, roll_max,  pitch_min,  pitch_max,  yaw_min,  yaw_max,tx_min, tx_max, ty_min, ty_max, tz_min, tz_max, ub_i);
    
#ifdef USE_GJK
    preprocess_lid(ref(point_cloud_model), ref(point_cloud_data), ref(model_voronoi_vertices),  ref(vec_node_i.valid_cells), ref(valid_cells_i),  ref(dist_cost_i), vec_node_i.roll.first, vec_node_i.roll.second, vec_node_i.pitch.first, vec_node_i.pitch.second, vec_node_i.yaw.first ,vec_node_i.yaw.second, ref(vec_node_i.tx.first), ref(vec_node_i.tx.second), ref(vec_node_i.ty.first), ref(vec_node_i.ty.second), ref(vec_node_i.tz.first) ,ref(vec_node_i.tz.second), upper_bound, ref(prep_time_i), ref(vec_lb_i), error_type, ref(dii), ref(djj), ref(model_voronoi_radius_sq));
#endif
    bool model_created=false;
    if(valid_cells_i.size()>=point_cloud_data.size() && valid_cells_i.size()<=1e4){
        if(error_type=="L2"){
            model_i=Reg_L2_model_rotation_trigonometric(point_cloud_model, point_cloud_data, vec_node_i.roll.first, vec_node_i.roll.second, vec_node_i.pitch.first, vec_node_i.pitch.second, vec_node_i.yaw.first ,vec_node_i.yaw.second, vec_node_i.tx.first, vec_node_i.tx.second, vec_node_i.ty.first, vec_node_i.ty.second, vec_node_i.tz.first ,vec_node_i.tz.second,valid_cells_i,dist_cost_i,incomp);
        }
        else{
            DebugOn("L1 objective not supported");
        }
        model_created=true;
        m_vec_i=1;
    }
    else if(valid_cells_i.size()> 1e4){
        m_vec_i=10;
    }
    else if(valid_cells_i.size()<point_cloud_data.size()){
        m_vec_i=0;
    }
}

void run_preprocess_parallel_Align(const vector<vector<double>>& point_cloud_model,
                                   const vector<vector<double>>& point_cloud_data, const vector<vector<vector<double>>>& model_voronoi_vertices, vector<int>& pos_vec, vector<shared_ptr<Model<double>>>& models, vector<treenode_p>& vec_node, vector<int>& m_vec,vector<double>& vec_lb, vector<indices>& valid_cells, int nb_threads, double upper_bound, double lower_bound, vector<param<double>>& dist_cost_new, int iter, string error_type, vector<double>& ub_node,double roll_min, double roll_max,  double pitch_min,  double pitch_max,  double yaw_min,  double yaw_max, double tx_min, double tx_max, double ty_min, double ty_max, double tz_min, double tz_max,const nanoflann::KDTreeSingleIndexAdaptor<
                                   nanoflann::L2_Simple_Adaptor<double, PointCloud<double>>,
                                   PointCloud<double>, 3 /* dim */>& index, const map<int, vector<int>>& incomp, const param<double>& dii, const param<double>& djj, const vector<double>& model_voronoi_radius_sq){
    vector<shared_ptr<Model<double>>> temp_models;
    std::vector<thread> threads;
    ub_node.resize(7,1000);
    int nd=point_cloud_data.size();
    int num=vec_node.size();
    if(num==0){
        DebugOff("in run_parallel(models...), models is empty, returning");
    }
    vector<vector<double>> vec_ub;
    vector<param<double>> vec_dist_cost;
    for(auto i=0;i<num;i++){
        param<double> dist_cost("dist_cost");
        vec_dist_cost.push_back(dist_cost);
    }
    
    valid_cells.resize(num);
    m_vec.resize(num, 0);
    vec_lb.resize(num, 0.0);
    temp_models.resize(num);
    vec_ub.resize(num);
    
    vector<double> vec_prep_time;
    vec_prep_time.resize(num, 0.0);
    for (auto i = 0; i < num; i++) {
        threads.push_back(thread(&run_preprocess_model_Align, ref(point_cloud_model), ref(point_cloud_data),ref(model_voronoi_vertices), ref(vec_node[i]), ref(m_vec[i]),ref(vec_lb[i]), ref(valid_cells[i]), ref(vec_dist_cost[i]), ref(vec_prep_time[i]), upper_bound, ref(temp_models[i]), error_type, ref(vec_ub[i]),roll_min, roll_max, pitch_min, pitch_max,yaw_min, yaw_max, tx_min, tx_max,ty_min, ty_max, tz_min, tz_max, ref(index), ref(incomp),ref(dii), ref(djj), ref(model_voronoi_radius_sq)));
    }
    for(auto &t : threads){
        t.join();
    }
    threads.clear();
    
    for(auto i=0;i<num;i++){
        vec_lb[i]=std::max(vec_lb[i], vec_node[i].lb);
        if(m_vec[i]==1){
            models.push_back(temp_models[i]);
            pos_vec.push_back(i);
        }
    }
    ub_node[0]=vec_ub[0][0];
    for(auto i=0;i<num;i++){
        if(vec_ub[i][0]<=ub_node[0]){
            ub_node[0]=vec_ub[i][0];
            ub_node[1]=vec_ub[i][1];
            ub_node[2]=vec_ub[i][2];
            ub_node[3]=vec_ub[i][3];
            ub_node[4]=vec_ub[i][4];
            ub_node[5]=vec_ub[i][5];
            ub_node[6]=vec_ub[i][6];
        }
    }
    dist_cost_new=vec_dist_cost;
}
vector<double> BranchBound_Align(vector<vector<double>>& point_cloud_model, vector<vector<double>>& point_cloud_data,const vector<vector<vector<double>>>& model_voronoi_vertices, vector<double>& best_rot, double best_ub, double tx_min, double tx_max, double ty_min, double ty_max, double tz_min, double tz_max, std::string error_type, const nanoflann::KDTreeSingleIndexAdaptor<
                                 nanoflann::L2_Simple_Adaptor<double, PointCloud<double>>,
                                 PointCloud<double>, 3 /* dim */>& index,const param<double>& dii, const param<double>& djj, const vector<double>& model_voronoi_radius_sq)
{
#ifdef USE_MPI
    int worker_id, nb_workers;
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
    auto err_size = MPI_Comm_size(MPI_COMM_WORLD, &nb_workers);
#endif
    vector<double> rpyt, rpy_rad, txyz;
    auto roll_rad=best_rot[0];
    auto pitch_rad=best_rot[1];
    auto yaw_rad = best_rot[2];
    rpy_rad.push_back(roll_rad);
    rpy_rad.push_back(pitch_rad);
    rpy_rad.push_back(yaw_rad);
    txyz.push_back(best_rot[3]);
    txyz.push_back(best_rot[4]);
    txyz.push_back(best_rot[5]);
    /* INPUT BOUNDS */
    
    /* INPUT BOUNDS */
    double time_start = get_wall_time();
    double total_time_max = 900000;
    double prep_time_total=0;
    
    
    
    double yaw_min = -90*pi/180., yaw_max = 90*pi/180., pitch_min =-90*pi/180.,pitch_max = 90*pi/180.,roll_min =-90*pi/180.,roll_max = 90*pi/180.;
    
    int nd=point_cloud_data.size();
    int nm=point_cloud_model.size();
    indices N1 = range(1,nd);
    indices N2 = range(1,nm);
    vector<int> new_matching(nd);
    vector<double> res(nd);
    double max_time = 5;
    double max_time_init=5;
    bool max_time_increase=false;
    int max_iter = 1e6;
    int models_count=0, models_new_count=0;
    int infeasible_count=0;
    vector<pair<pair<int,int>,pair<int,int>>> incompatible_pairs;
    size_t nb_threads = std::thread::hardware_concurrency();
    pair<double,double> roll_bounds_r, pitch_bounds_r, yaw_bounds_r,tx_bounds_r,ty_bounds_r,tz_bounds_r;
    
    roll_bounds_r={roll_min, roll_max};
    pitch_bounds_r={pitch_min, pitch_max};
    yaw_bounds_r={yaw_min, yaw_max};
    tx_bounds_r={tx_min, tx_max};
    ty_bounds_r={ty_min, ty_max};
    tz_bounds_r={tz_min, tz_max};
    
    vector<pair<double,double>> roll_bounds, pitch_bounds, yaw_bounds,tx_bounds,ty_bounds,tz_bounds;
    
    vector<indices> valid_cells;
    vector<int> pos_vec;
    vector<double> vec_lb;
    vector<treenode_p> vec_node;
    vector<int> m_vec;
    vector<vector<double>> costs_upto_vec;
    auto point_cloud_data_copy=point_cloud_data;
    DebugOn("I will be using " << nb_threads << " parallel threads" << endl);
    vector<shared_ptr<Model<>>> models, models_new;
    double lb = 0, ub = 12, ub_=-1, best_lb = 0;
    int nb_pruned = 0;
    int depth_r=0, iter=0;
    vector<int> depth_vec, depth_vec_new;
    priority_queue<treenode_p> lb_queue;
    vector<double> costs_upto_init(nd,0.0);
    
    double min_cost_sum=0.0;
    
    param<double> dist_cost_r("dist_cost_r");
    indices valid_cells_r;
    vector<param<double>> dist_cost_cells;
    double  prep_time=0.0;
    min_cost_sum=0;
    
    map<int, vector<int>> incomp;
    
    
    for(auto i=0;i<point_cloud_data.size()-1;i++){
        vector<int> red;
        for(auto j=i+1;j<point_cloud_data.size();j++){
            auto d=pow(point_cloud_data.at(i)[0]-point_cloud_data.at(j)[0],2)+
            pow(point_cloud_data.at(i)[1]-point_cloud_data.at(j)[1],2)+
            pow(point_cloud_data.at(i)[2]-point_cloud_data.at(j)[2],2);
            if(d>=3*best_ub+1e-9){
                red.push_back(j);
                DebugOff("cannot match with same j"<<endl);
            }
        }
        if(red.size()>=1){
            incomp[i]=red;
        }
    }
    
    
    
    lb_queue.push(treenode_p(roll_bounds_r, pitch_bounds_r, yaw_bounds_r, tx_bounds_r, ty_bounds_r,tz_bounds_r,lb, ub, ub_, depth_r, valid_cells_r, false, dist_cost_r));
    treenode_p topnode=lb_queue.top();
    
    best_lb = lb_queue.top().lb;
    double elapsed_time = get_wall_time() - time_start;
    double opt_gap = (best_ub - best_lb)/best_ub;
    double opt_gap_abs=(best_ub - best_lb);
    double max_opt_gap = 0.01;/* 5% opt gap */
    double eps=0.001;
    int prep_count=0;
    double ut_total=0;
    while (elapsed_time < total_time_max && lb_queue.top().lb<=best_ub && !lb_queue.empty() && opt_gap > max_opt_gap && !lb_queue.top().leaf) {
        best_lb = lb_queue.top().lb;
        opt_gap = (best_ub - best_lb)/best_ub;
        DebugOn("Best UB so far = " << to_string_with_precision(best_ub,9) << endl);
        DebugOn("Best LB so far = " << to_string_with_precision(best_lb,9) << endl);
        DebugOn("Opt gap so far = " << to_string_with_precision(opt_gap*100,6) << "%\n");
        DebugOn("Queue size = " << lb_queue.size() << "\n");
        DebugOn("Elapsed time = " << elapsed_time << "seconds\n");
        DebugOn("iter "<<iter<<endl);
        if(elapsed_time >= total_time_max || opt_gap <= max_opt_gap)
            break;
        DebugOn("Total infeasible =  " << infeasible_count << endl);
        DebugOn("Total prep_time =  " << prep_time_total << endl);
        DebugOn("Total discarded =  " << prep_count << endl);
        double max_incr=0, max_ratio=1;
        pos_vec.clear();
        models.clear();
        roll_bounds.clear();
        pitch_bounds.clear();
        yaw_bounds.clear();
        tx_bounds.clear();
        ty_bounds.clear();
        tz_bounds.clear();
        valid_cells.clear();
        depth_vec.clear();
        m_vec.clear();
        vec_node.clear();
        vec_lb.clear();
        dist_cost_cells.clear();
        costs_upto_vec.clear();
        iter++;
        models_count=0;
        models_new_count=0;
        topnode=lb_queue.top();
        prep_time_total=0;
        ut_total=0;
        int step=8;
        for(auto i=0;i<nb_threads;i+=step){
            if(lb_queue.top().lb<=best_ub && !lb_queue.top().leaf && !lb_queue.empty()){
                topnode=lb_queue.top();
                lb_queue.pop();
                if((topnode.depth%2==0)){
                    DebugOn("R branch "<<topnode.depth<<endl);
                    double roll_increment,  pitch_increment, yaw_increment;
                    roll_increment = (topnode.roll.second - topnode.roll.first)/2.0;
                    pitch_increment = (topnode.pitch.second - topnode.pitch.first)/2.0;
                    yaw_increment = (topnode.yaw.second - topnode.yaw.first)/2.0;
                    roll_bounds.push_back({topnode.roll.first, topnode.roll.first+roll_increment});
                    roll_bounds.push_back({topnode.roll.first, topnode.roll.first+roll_increment});
                    roll_bounds.push_back({topnode.roll.first, topnode.roll.first+roll_increment});
                    roll_bounds.push_back({topnode.roll.first, topnode.roll.first+roll_increment});
                    roll_bounds.push_back({topnode.roll.first+roll_increment, topnode.roll.second});
                    roll_bounds.push_back({topnode.roll.first+roll_increment, topnode.roll.second});
                    roll_bounds.push_back({topnode.roll.first+roll_increment, topnode.roll.second});
                    roll_bounds.push_back({topnode.roll.first+roll_increment, topnode.roll.second});
                    pitch_bounds.push_back({topnode.pitch.first, topnode.pitch.first+pitch_increment});
                    pitch_bounds.push_back({topnode.pitch.first, topnode.pitch.first+pitch_increment});
                    pitch_bounds.push_back({topnode.pitch.first+pitch_increment, topnode.pitch.second});
                    pitch_bounds.push_back({topnode.pitch.first+pitch_increment, topnode.pitch.second});
                    pitch_bounds.push_back({topnode.pitch.first, topnode.pitch.first+pitch_increment});
                    pitch_bounds.push_back({topnode.pitch.first, topnode.pitch.first+pitch_increment});
                    pitch_bounds.push_back({topnode.pitch.first+pitch_increment, topnode.pitch.second});
                    pitch_bounds.push_back({topnode.pitch.first+pitch_increment, topnode.pitch.second});
                    yaw_bounds.push_back({topnode.yaw.first, topnode.yaw.first+yaw_increment});
                    yaw_bounds.push_back({topnode.yaw.first+yaw_increment, topnode.yaw.second});
                    yaw_bounds.push_back({topnode.yaw.first, topnode.yaw.first+yaw_increment});
                    yaw_bounds.push_back({topnode.yaw.first+yaw_increment, topnode.yaw.second});
                    yaw_bounds.push_back({topnode.yaw.first, topnode.yaw.first+yaw_increment});
                    yaw_bounds.push_back({topnode.yaw.first+yaw_increment, topnode.yaw.second});
                    yaw_bounds.push_back({topnode.yaw.first, topnode.yaw.first+yaw_increment});
                    yaw_bounds.push_back({topnode.yaw.first+yaw_increment, topnode.yaw.second});
                    for(auto k=0;k<8;k++){
                        tx_bounds.push_back({topnode.tx.first, topnode.tx.second});
                        ty_bounds.push_back({topnode.ty.first, topnode.ty.second});
                        tz_bounds.push_back({topnode.tz.first, topnode.tz.second});
                    }
                }
                else{
                    DebugOn("t branch "<<topnode.depth<<endl);
                    double tx_increment,  ty_increment, tz_increment;
                    if(false && topnode.tx.first<=-0.001 && topnode.tx.second>=0.001){
                        tx_increment = topnode.tx.first*(-1);
                    }
                    else{
                        tx_increment = (topnode.tx.second - topnode.tx.first)/2.0;
                    }
                    if(false && topnode.ty.first<=-0.001 && topnode.ty.second>=0.001){
                        ty_increment = topnode.ty.first*(-1);
                    }
                    else{
                        ty_increment = (topnode.ty.second - topnode.ty.first)/2.0;
                    }
                    if(false && topnode.tz.first<=-0.001 && topnode.tz.second>=0.001){
                        tz_increment = topnode.tz.first*(-1);
                    }
                    else{
                        tz_increment = (topnode.tz.second - topnode.tz.first)/2.0;
                    }
                    tx_bounds.push_back({topnode.tx.first, topnode.tx.first+tx_increment});
                    tx_bounds.push_back({topnode.tx.first, topnode.tx.first+tx_increment});
                    tx_bounds.push_back({topnode.tx.first, topnode.tx.first+tx_increment});
                    tx_bounds.push_back({topnode.tx.first, topnode.tx.first+tx_increment});
                    tx_bounds.push_back({topnode.tx.first+tx_increment, topnode.tx.second});
                    tx_bounds.push_back({topnode.tx.first+tx_increment, topnode.tx.second});
                    tx_bounds.push_back({topnode.tx.first+tx_increment, topnode.tx.second});
                    tx_bounds.push_back({topnode.tx.first+tx_increment, topnode.tx.second});
                    ty_bounds.push_back({topnode.ty.first, topnode.ty.first+ty_increment});
                    ty_bounds.push_back({topnode.ty.first, topnode.ty.first+ty_increment});
                    ty_bounds.push_back({topnode.ty.first+ty_increment, topnode.ty.second});
                    ty_bounds.push_back({topnode.ty.first+ty_increment, topnode.ty.second});
                    ty_bounds.push_back({topnode.ty.first, topnode.ty.first+ty_increment});
                    ty_bounds.push_back({topnode.ty.first, topnode.ty.first+ty_increment});
                    ty_bounds.push_back({topnode.ty.first+ty_increment, topnode.ty.second});
                    ty_bounds.push_back({topnode.ty.first+ty_increment, topnode.ty.second});
                    tz_bounds.push_back({topnode.tz.first, topnode.tz.first+tz_increment});
                    tz_bounds.push_back({topnode.tz.first+tz_increment, topnode.tz.second});
                    tz_bounds.push_back({topnode.tz.first, topnode.tz.first+tz_increment});
                    tz_bounds.push_back({topnode.tz.first+tz_increment, topnode.tz.second});
                    tz_bounds.push_back({topnode.tz.first, topnode.tz.first+tz_increment});
                    tz_bounds.push_back({topnode.tz.first+tz_increment, topnode.tz.second});
                    tz_bounds.push_back({topnode.tz.first, topnode.tz.first+tz_increment});
                    tz_bounds.push_back({topnode.tz.first+tz_increment, topnode.tz.second});
                    for(auto k=0;k<8;k++){
                        roll_bounds.push_back({topnode.roll.first, topnode.roll.second});
                        pitch_bounds.push_back({topnode.pitch.first, topnode.pitch.second});
                        yaw_bounds.push_back({topnode.yaw.first, topnode.yaw.second});
                    }
                }
                for(auto k=0;k<8;k++){
                    vec_node.push_back(treenode_p(roll_bounds[i+k],  pitch_bounds[i+k], yaw_bounds[i+k], tx_bounds[i+k],  ty_bounds[i+k], tz_bounds[i+k], topnode.lb, best_ub, -1.0, topnode.depth+1, topnode.valid_cells, false,topnode.dist_cost_cells));
                    depth_vec.push_back(topnode.depth+1);
                }
            }
            else{
                break;
            }
            if(lb_queue.empty()){
                break;
            }
            topnode = lb_queue.top();
        }
        elapsed_time = get_wall_time() - time_start;
        DebugOn("Elapsed time = " << elapsed_time << "seconds\n");
        if(elapsed_time + max_time > total_time_max){
            DebugOn("max time "<< max_time);
            break;
        }
        vector<double> ub_all(7);
        run_preprocess_parallel_Align(point_cloud_model,point_cloud_data, model_voronoi_vertices, pos_vec, models, vec_node, m_vec, vec_lb, valid_cells, nb_threads, best_ub, best_lb, dist_cost_cells, iter, error_type, ub_all, roll_min, roll_max, pitch_min, pitch_max, yaw_min, yaw_max, tx_min, tx_max, ty_min, ty_max, tz_min, tz_max,index, incomp,dii, djj, model_voronoi_radius_sq);
        if(ub_all[0]<=best_ub){
            best_ub=ub_all[0];
            rpy_rad[0]=ub_all[1];
            rpy_rad[1]=ub_all[2];
            rpy_rad[2]=ub_all[3];
            txyz[0]=ub_all[4];
            txyz[1]=ub_all[5];
            txyz[2]=ub_all[6];
        }
        for (int j = 0; j<m_vec.size(); j++) {
            if(m_vec[j]==0){
                prep_count++;
            }
            else if(m_vec[j]==10){
                lb_queue.push(treenode_p(roll_bounds[j], pitch_bounds[j], yaw_bounds[j],tx_bounds[j], ty_bounds[j], tz_bounds[j], vec_lb[j], best_ub, ub_, depth_vec[j], valid_cells[j], false, dist_cost_cells[j]));
            }
        }
        DebugOn("models size "<<models.size()<<endl);
        
        run_parallel(models, gurobi, 1e-6, nb_threads, "", max_iter, max_time, (best_ub+1e-6));
        for (int j = 0; j<models.size(); j++) {
            int pos=pos_vec[j];
            if(models[j]->_status==0){
                auto lb_ =  models[j]->get_rel_obj_val();
                auto leaf_node=false;
                if(true){
                    lb = std::max(models[j]->get_rel_obj_val(), vec_lb[pos]);
                }
                if(lb-1e-4<=best_ub)
                {
                    lb_queue.push(treenode_p(roll_bounds[pos], pitch_bounds[pos], yaw_bounds[pos], tx_bounds[pos],ty_bounds[pos],tz_bounds[pos], lb, best_ub, ub_, depth_vec[pos], valid_cells[pos], leaf_node, dist_cost_cells[pos]));
                }
                else{
                    DebugOn("Infeasible lb "<<lb<<" "<<"best_ub "<<best_ub<<endl);
                    infeasible_count++;
                }
            }
            else{
                infeasible_count++;
            }
        }
        opt_gap = (best_ub - best_lb)/best_ub;
        opt_gap_abs=best_ub-best_lb;
        elapsed_time = get_wall_time() - time_start;
        
    }
    DebugOn("UB final "<<best_ub<<endl);
    DebugOn("LB final "<<best_lb<<endl);
    DebugOn("Gap final "<<(best_ub-best_lb)/best_ub*100.0<<endl);
    DebugOn("Elapsed time "<<elapsed_time<<endl);
    DebugOn("Total iter "<<iter<<endl);
    DebugOn("Queue size = " << lb_queue.size() << "\n");
    DebugOn("lb que top = " << lb_queue.top().lb << "\n");
    roll_rad= rpy_rad[0];
    pitch_rad=rpy_rad[1];
    yaw_rad = rpy_rad[2];
    auto tx=txyz[0];
    auto ty=txyz[1];
    auto tz=txyz[2];
    DebugOn("roll rad "<< roll_rad<<endl);
    DebugOn("pitch rad "<< pitch_rad<<endl);
    DebugOn("yaw rad "<< yaw_rad<<endl);
    while(!lb_queue.empty())
    {
        auto node = lb_queue.top();
        DebugOn("node lb "<<node.lb<<" node.leaf "<<node.leaf<<endl);
        
        DebugOn(node.roll.first<<" "<< node.roll.second<<" "<<node.pitch.first<<" "<<node.pitch.second<<" "<<node.yaw.first<<" "<<node.yaw.second<<endl);
        DebugOn(node.tx.first<<" "<< node.tx.second<<" "<<node.ty.first<<" "<<node.ty.second<<" "<<node.tz.first<<" "<<node.tz.second<<endl);
        if(node.roll.first-1e-6<=roll_rad && roll_rad<=node.roll.second+1e-6 && node.pitch.first-1e-6<=pitch_rad && pitch_rad<=node.pitch.second+1e-6 && node.yaw.first-1e-6<=yaw_rad && yaw_rad<=node.yaw.second+1e-6){
            if(node.tx.first-1e-6<=tx && tx<=node.tx.second+1e-6 && node.ty.first-1e-6<=ty && ty<=node.ty.second+1e-6 && node.tz.first-1e-6<=tz && tz<=node.tz.second+1e-6){
                DebugOn("True interval contained "<<endl);
            }
        }
        lb_queue.pop();
    }
    
    DebugOn("roll rad "<< roll_rad<<endl);
    DebugOn("pitch rad "<< pitch_rad<<endl);
    DebugOn("yaw rad "<< yaw_rad<<endl);
    DebugOn("tx "<< tx<<endl);
    DebugOn("ty "<< ty<<endl);
    DebugOn("tz "<< tz<<endl);
    
    auto roll=roll_rad*180/pi;
    auto pitch=pitch_rad*180/pi;
    auto yaw=yaw_rad*180/pi;
    
    
    DebugOn("roll deg"<< roll<<endl);
    DebugOn("pitch deg"<< pitch<<endl);
    DebugOn("yaw deg"<< yaw<<endl);
    DebugOn("tx "<< tx<<endl);
    DebugOn("ty "<< ty<<endl);
    DebugOn("tz "<< tz<<endl);
    
    rpyt.push_back(roll_rad);
    rpyt.push_back(pitch_rad);
    rpyt.push_back(yaw_rad);
    rpyt.push_back(tx);
    rpyt.push_back(ty);
    rpyt.push_back(tz);
    return rpyt;
}
#ifdef USE_MPI
void send_vector_new(const vector<size_t>& limits, vector<double>& vec_full, vector<double>& vec_worker){
    int worker_id, nb_workers;
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
    auto err_size = MPI_Comm_size(MPI_COMM_WORLD, &nb_workers);
    auto nb_workers_ =  limits.size()-1;
    DebugOff("I'm worker ID: " << worker_id << ", I'm getting ready to send my status " << endl);
    if(worker_id+1<limits.size()){
        for (auto i = limits[worker_id]; i < limits[worker_id+1]; i++) {
            vec_full[i]=vec_worker[i-limits[worker_id]];
        }
    }
    std::vector<int> d, counts;
    for(auto l=limits.begin()+1;l!=limits.end();l++){
        counts.push_back(*l-*(l-1));
        d.push_back(*(l-1));
    }
    for(auto l=nb_workers_;l!=nb_workers;l++){
        counts.push_back(0);
        d.push_back(limits.back());
    }
    if(counts.size()!=nb_workers){
        DebugOn("Error in size of counts");
    }
    if(d.size()!=nb_workers){
        DebugOn("Error in size of d");
    }
    MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                   &vec_full[0], &counts[0], &d[0], MPI_DOUBLE, MPI_COMM_WORLD);
    //MPI_Barrier(MPI_COMM_WORLD);
}
#endif
#endif /* BB_h */
