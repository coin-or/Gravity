//  Created by Smitha on 11/2/22.
//
#ifndef mink_model_h
#define mink_model_h
#include <gravity/solver.h>
#include <gravity/Net.h>
#include <stdio.h>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

using namespace gravity;
using namespace std;
using namespace Eigen;

double check_PSD_full_rip(shared_ptr<Model<double>>& m, int num_part){
    const double active_tol=1e-9;
    var<double> X=m->get_var<double>("X");
    var<double> Xij=m->get_var<double>("Xij");
    double neg_eig_value=1;

    int dim_full=X._indices->_keys->size();
    
    Eigen::MatrixXd mat_full(dim_full,dim_full);
    int count=0;
    vector<string> all_names;
    for(auto k:*X._indices->_keys){
        mat_full(count, count)=X.eval(k);
        all_names.push_back(k);
        count++;
    }
    
    for(auto i=0;i<all_names.size()-1;i++){
        for(auto j=i+1;j<all_names.size();j++){
            auto k=all_names[i]+","+all_names[j];
            if (Xij._indices->has_key(k)){
                mat_full(i, j)=Xij.eval(k);
                mat_full(j, i)=Xij.eval(k);
            }
            else{
                mat_full(i, j)=1;
                mat_full(j, i)=1;
            }
        }
    }
Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
es.compute(mat_full);
    DebugOn("full eig values ");
for(auto m=0;m<dim_full;m++){
    if(es.eigenvalues()[m]<=-active_tol)
        DebugOn("negative eigen value "<<endl);
    if(es.eigenvalues()[m]<=neg_eig_value)
        neg_eig_value=es.eigenvalues()[m];
    DebugOn(std::setprecision(12)<<es.eigenvalues()[m]<<"\t");
    }

return neg_eig_value;
}
/*For testing sparse primal RIP models; sparse data files are created from full data files; read_misdp may be used to test dense dual RIP model*/
Net model_rip(string fname, shared_ptr<Model<double>>& m) {
    int k=5;
Net g;
  int Num_nodes = 0;
  int Num_edges = 0;
  ifstream infile(fname);
  string sLine;

  if (infile.good()) {
    getline(infile, sLine);
    istringstream iss(sLine);
    iss >> Num_nodes;
    iss >> Num_edges;
  } else {
    fprintf(stderr, "can’t open input file %s\n", fname.c_str());
    exit(1);
  }
  string name;

  Node* node = nullptr;
    indices nodes;
    DebugOn("nodes "<<Num_nodes<<endl);
    DebugOn("edges "<<Num_edges<<endl);

  for (int i = 0; i < Num_nodes; i++) {
    name = to_string(i);
    node = new Node(name);
    g.add_node(node);
    nodes.insert(name);
  }

  // get arcs
  Arc* arc = NULL;
  Arc* arc_clone = NULL;
  Arc* arc_chordal = NULL;

  // note that src, dest are names of nodes.
  string src, dest;
  int src_a, dest_a;
  double weight;
  param<double> w("w"), v("v");
  while (getline(infile, sLine, '\n')) {
    istringstream iss(sLine);
    iss >> src_a >> dest_a >> weight;
    // cout << src  << ", " << dest << ", " << weight << endl;
      if(std::abs(weight)>=1e-9){
      if(src_a!=dest_a){
      if(src_a<dest_a){
          src=to_string(src_a);
          dest=to_string(dest_a);
      }
      else{
          src=to_string(dest_a);
          dest=to_string(src_a);
      }
          

    name = src + "," + dest;  //

    arc = new Arc(name);
    arc->_src = g.get_node(src);
    arc->_dest = g.get_node(dest);
    arc->_weight = weight;
    g.add_arc(arc);
    arc->connect();
    w.add_val(name, weight);
      }
      else{
          v.add_val(to_string(src_a), weight);
      }
      }
  }
infile.close();
    
auto node_pairs=g.get_node_pairs();

g.sdp_3d_cuts=false;



 
    


g.get_tree_decomp_bags();
auto bags_3d=g.decompose_bags_3d();
m->sdp_dual=false;

//auto node_pairs_chord = g.get_node_pairs_chord;
auto node_pairs_chord = g.get_node_pairs_chord(g._bags);
    
    vector<vector<pair<pair<size_t,size_t>, double>>> Xij_x_map(node_pairs_chord.size());
    vector<vector<pair<pair<size_t,size_t>, double>>> Xii_x_map(nodes.size());
    vector<vector<pair<pair<size_t,size_t>, double>>> Xij_y_map(node_pairs_chord.size());
    vector<vector<pair<pair<size_t,size_t>, double>>> Xii_y_map(nodes.size());
    vector<double> Xij_cons_map(node_pairs_chord.size(),0);
    vector<double> Xii_cons_map(nodes.size(),0);

var<> X("X", 0, 1);
m->add(X.in(nodes));
   // X.initialize_all(1.0/(Num_nodes*1.0));

var<> Xij("Xij", -1.0, 1.0);
m->add(Xij.in(node_pairs_chord));


var<int> y("y",0,1);
m->add(y.in(nodes));

Constraint<> sum_Xii("sum_Xii");
sum_Xii=sum(X)-1;
m->add(sum_Xii==0);
    
    Constraint<> sum_y("sum_y");
    sum_y=sum(y)-k;
    m->add(sum_y<=0);
    
    Constraint<> Xi_yi("Xi_yi");
    Xi_yi=X-y;
    m->add(Xi_yi.in(nodes)<=0);
    
    Constraint<> yi_Xi("yi_Xi");
    yi_Xi=y*(-1)-X;
   // m->add(yi_Xi.in(nodes)<=0);

//    Constraint<> Xij_yi("Xij_yi");
//    Xij_yi=Xij-y.to(node_pairs_chord);
//    m->add(Xij_yi.in(node_pairs_chord) <= 0);
//
//    Constraint<> yi_Xij("yi_Xij");
//    yi_Xij=y.to(node_pairs_chord)*(-1)-Xij;
//    m->add(yi_Xij.in(node_pairs_chord) <= 0);
//
//    Constraint<> Xij_yia("Xij_yia");
//    Xij_yia=Xij-y.from(node_pairs_chord);
//    m->add(Xij_yia.in(node_pairs_chord) <= 0);
//
//    Constraint<> yi_Xija("yi_Xija");
//    yi_Xija=y.from(node_pairs_chord)*(-1)-Xij;
//    m->add(yi_Xija.in(node_pairs_chord) <= 0);
//
    Constraint<> SOC("SOC");
    SOC = Xij*Xij - X.from(node_pairs_chord)*X.to(node_pairs_chord);
    SOC.add_to_callback();
    m->add(SOC.in(node_pairs_chord) <= 0);
    
    auto bag_size = bags_3d.size();
    auto Wij_ = Xij.pairs_in_bags(bags_3d, 3);
    auto Wii_ = X.in_bags(bags_3d, 3);
    
    Constraint<> SDP3("SDP_3D");
    
    SDP3 += (pow(Wij_[0], 2)) * Wii_[2];
    SDP3 += (pow(Wij_[1], 2)) * Wii_[0];
    SDP3 += (pow(Wij_[2], 2)) * Wii_[1];
    SDP3 -= 2 * Wij_[0] * (Wij_[1] * Wij_[2]);
    SDP3 -= Wii_[0] * Wii_[1] * Wii_[2];
    SDP3.add_to_callback();
    m->add(SDP3.in(range(0, bag_size-1)) <= 0);
    
    int count=0;
    std::vector<pair<int,std::vector<string>>> _bag_names;
    m->sdp_dual=false;
    for(auto b:g._bags){
        pair<int,vector<string>> bn;
        if(b.second.size()>=3){
            bn.first=count++;
            DebugOn("bag "<<count<<endl);
            for(auto n:b.second){
                bn.second.push_back(n->_name);
                DebugOn(n->_name <<"\t");
            }
            DebugOn(endl);
            _bag_names.push_back(bn);
        }
    }
    
    for(auto k:*node_pairs_chord._keys){
        if(!node_pairs.has_key(k))
            w.add_val(k,0);
    }
    DebugOn("true edges "<<node_pairs.size());
func<> obj=w.tr()*Xij+v.tr()*X;
m->min(obj);
    
    DebugOn("largest ub "<<obj._range->second<<endl);

m->print();
    var<double> x("x", 0, 0);
    m->add(x.in(range(0,0)));
    
    
    m->_bag_names=_bag_names;
    m->Xii_x_map=Xii_x_map;
    m->Xii_y_map=Xii_y_map;
    m->Xij_x_map=Xij_x_map;
    m->Xij_y_map=Xij_y_map;
    m->Xii_cons_map=Xii_cons_map;
    m->Xij_cons_map=Xij_cons_map;
    m->make_PSD(X, Xij);
    
return g;

}


#endif /* read_misdp_hpp */
