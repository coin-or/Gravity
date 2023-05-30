//
//  misdp.cpp
//  misdp
//
//  Created by Smitha on 7/7/22.
//

#include <stdio.h>
#include "mink_model.h"
#include <gravity/solver.h>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <fstream>


using namespace gravity;
using namespace std;

int main(int argc, char * argv[]){
    //string fname=string(prj_dir)+"/data_sets/MISDP/2x3_3bars.cbf";
    string fname=string(prj_dir)+"/data_sets/MISDP/band50_3.txt";
    //string fname=string(prj_dir)+"/data_sets/MISDP/Block_30_20_4_6.txt";
    //string fname=string(prj_dir)+"/data_sets/MISDP/spinglass2g_33.txt";
    //string fname=string(prj_dir)+"/data_sets/MISDP/spinglass2g_1111.txt";
    //string fname=string(prj_dir)+"/data_sets/MISDP/2x7_3bars.cbf";
    //string fname=string(prj_dir)+"/data_sets/MISDP/2x4_2scen_3bars.cbf";
    //string fname=string(prj_dir)+"/data_sets/MISDP/coloncancer_1_100_5.cbf";
    //string fname=string(prj_dir)+"/data_sets/MISDP/2g_4_164_k3_5_6.cbf";
    
    //string fname=string(prj_dir)+"/data_sets/MISDP/spinglass2g_1111.txt";

    
    bool root_refine = true, add_soc=true, add_threed=true, add_bag=true, hierarc=false;
    string root_refine_s = "false", add_soc_s="false", add_threed_s="false", add_bag_s="false", hierarc_s="false";
    if(argc>=2){
        fname=argv[1];
    }
    if(argc>=3){
        root_refine_s=argv[2];
        if(root_refine_s.compare("true")==0)
            root_refine = true;
    }
    if(argc>=4){
        add_soc_s=argv[3];
        if(add_soc_s=="true")
            add_soc = true;
    }
    if(argc>=5){
        add_threed_s=argv[4];
        if(add_threed_s=="true")
            add_threed = true;
    }
    if(argc>=6){
        write_mink_cbf(fname);
        exit(0);
    }
    auto m=make_shared<Model<double>>("mink");
    
    int num_nodes=0, num_sparse_edges=0, num_part=0,width=0;
    auto g=model_mink(fname, m, num_nodes,num_sparse_edges, num_part, width);
    m->add_soc=add_soc;
    m->add_threed=add_threed;
    m->add_bag=add_bag;
    m->add_hierarc=hierarc;
    //m->print();
    
    solver<> sc(m,gurobi);
    bool relax = false;
    auto ts=get_wall_time();
    if(root_refine){
        sc.run(relax = true);
    }
    sc.run();
    auto tf=get_wall_time();
    m->print_solution();
    m->print_constraints_stats(1e-9);
    
    auto eig_value=m->check_PSD_bags();
    
    double num_dense_edges=num_nodes*(num_nodes-1)*0.5;
    double matrix_density=num_sparse_edges/num_dense_edges;
    string out_file_name=fname;
    auto pos=out_file_name.find_last_of("/");
    out_file_name=out_file_name.substr(pos+1);
    pos=out_file_name.find_first_of(".");
    auto out_file_name1=out_file_name.substr(0,pos);
    out_file_name=string(prj_dir)+"/results_misdp/"+out_file_name1+".txt";
    ofstream fout(out_file_name.c_str());
    fout<<"Num_nodes "<<num_nodes<<" Num_sparse_edges "<<num_sparse_edges<<" Num_part "<<num_part<<" Num_dense_edges "<<num_dense_edges<<" Density "<<matrix_density<<" Width "<<width<<endl;
    fout<<out_file_name1<<"&"<<m->get_obj_val()<<"&"<<tf-ts<<"&"<<eig_value<<"&"<<m->_rel_obj_val<<"&"<<m->num_cuts[0]<<"&"<<m->num_cuts[1]<<"&"<<m->num_cuts[2]<<"&"<<m->num_cuts[3]<<"&"<<m->num_cuts[4]<<"&"<<m->num_cuts[5]<<m->num_cuts[6]<<"&"<<m->num_cuts[7]<<"&"<<m->num_cuts[8]<<"&"<<m->num_cuts[9]<<"&"<<m->num_cuts[10]<<"&"<<m->num_cuts[11]<<"\n";
    fout.close();
    cout<<"Num_nodes "<<num_nodes<<" Num_sparse_edges "<<num_sparse_edges<<" Num_part "<<num_part<<" Num_dense_edges "<<num_dense_edges<<" Density "<<matrix_density<<" max clique size "<<width<<endl;
    cout<<out_file_name1<<" obj "<<m->get_obj_val()<<" time "<<tf-ts<<" smallest eig "<<eig_value<<" lower bound  "<<m->_rel_obj_val<<endl;
    cout<<"soc incumbent "<<m->num_cuts[0]<<" threed incumbent "<<m->num_cuts[1]<<" eig bags incumbent "<<m->num_cuts[2]<<" eig full incumbent "<<m->num_cuts[3]<<" soc mipnode "<<m->num_cuts[4]<<" threed mipnode "<<m->num_cuts[5]<<" eig bags mipnode "<<m->num_cuts[6]<<" eig full mipnode "<<m->num_cuts[7]<<"\n";
}
