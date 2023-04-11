//
//  misdp.cpp
//  misdp
//
//  Created by Smitha on 7/7/22.
//

#include <stdio.h>
#include "rip_model.h"
#include "../CBF/read_misdp.h"
#include <gravity/solver.h>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <fstream>


using namespace gravity;
using namespace std;

int main(int argc, char * argv[]){
    //string fname=string(prj_dir)+"/data_sets/MISDP/2x3_3bars.cbf";
    //string fname=string(prj_dir)+"/data_sets/MISDP/rip_zeroone_153005l.txt";
    string fname=string(prj_dir)+"/data_sets/MISDP/rip_band404003l.txt";
    //string fname=string(prj_dir)+"/data_sets/MISDP/spinglass2g_33.txt";
    //string fname=string(prj_dir)+"/data_sets/MISDP/spinglass2g_1111.txt";
    //string fname=string(prj_dir)+"/data_sets/MISDP/2x7_3bars.cbf";
    //string fname=string(prj_dir)+"/data_sets/MISDP/2x4_2scen_3bars.cbf";
    //string fname=string(prj_dir)+"/data_sets/MISDP/coloncancer_1_100_5.cbf";
    //string fname=string(prj_dir)+"/data_sets/MISDP/2g_4_164_k3_5_6.cbf";
    
    

    bool root_refine = false;
    if(argc>=2){
        fname=argv[1];
    }
    if(argc>=3){
        root_refine = true;
    }
    if(argc>=4){
        /*Write sparse input file for RIP instance*/
        CBF_write_sparse_rip(fname.c_str());
    }
    else{
    auto m=make_shared<Model<double>>("rip");
    
    auto g=model_rip(fname, m);
    //m->print();
    
    auto ts=get_wall_time();
    solver<> sc(m, gurobi);
    if(root_refine)
        sc.run(true);
    sc.run(false);
    auto tf=get_wall_time();
    m->print_solution();
    m->print_constraints_stats(1e-9);
    
    auto eig_value1=m->check_PSD_bags();
    
    
    
    string out_file_name=fname;
    auto pos=out_file_name.find_last_of("/");
    out_file_name=out_file_name.substr(pos+1);
    pos=out_file_name.find_first_of(".");
    auto out_file_name1=out_file_name.substr(0,pos);
    out_file_name=string(prj_dir)+"/results_misdp/"+out_file_name1+".txt";
    ofstream fout(out_file_name.c_str());
    fout<<out_file_name1<<"&"<<m->get_obj_val()<<"&"<<tf-ts<<"&"<<eig_value1<<"&"<<m->_rel_obj_val<<"&"<<m->num_cuts[0]<<"&"<<m->num_cuts[1]<<"&"<<m->num_cuts[2]<<"&"<<m->num_cuts[3]<<"&"<<m->num_cuts[4]<<"&"<<m->num_cuts[5]<<"\n";
    fout.close();
    cout<<out_file_name1<<" obj "<<m->get_obj_val()<<" time "<<tf-ts<<" smallest eig "<<eig_value1<<" lower bound  "<<m->_rel_obj_val<<endl;
    cout<<"soc incumbent "<<m->num_cuts[0]<<"threed incumbent "<<m->num_cuts[1]<<" eig bags incumbent "<<m->num_cuts[2]<<" eig full incumbent "<<m->num_cuts[3]<<" soc mipnode "<<m->num_cuts[4]<<" eig bags mipnode "<<m->num_cuts[6]<<" eig full mipnode "<<m->num_cuts[7]<<"\n";
    }
}
