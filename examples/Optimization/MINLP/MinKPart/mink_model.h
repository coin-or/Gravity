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

double check_PSD_full_mink(shared_ptr<Model<double>>& m, int num_part){
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
/*Mink in primal form*/
Net model_mink(string fname, shared_ptr<Model<double>>& m) {
Net g;
  int num_part=3;
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
    name = to_string(i + 1);
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
param<double> w("w");
  while (getline(infile, sLine, '\n')) {
    istringstream iss(sLine);
    iss >> src_a >> dest_a >> weight;
    // cout << src  << ", " << dest << ", " << weight << endl;
    
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
infile.close();
    
auto node_pairs=g.get_node_pairs();

    
map<string, vector<pair<string, double>>> map_x;
map<string, vector<pair<string, double>>> map_y;
map<string, double> map_const;




g.sdp_3d_cuts=false;
g.get_tree_decomp_bags();
m->sdp_dual=false;

//auto node_pairs_chord = g.get_node_pairs_chord;
auto node_pairs_chord = g.get_node_pairs_chord(g._bags);

var<> X("X", 1, 1);
X._psd_diag=true;
m->add(X.in(nodes));
var<> Xij("Xij", -1.0/(num_part-1.0), 1);
Xij._psd_off_diag=true;
m->add(Xij.in(node_pairs_chord));
auto X_ptr = m->get_ptr_var<double>("X");
auto Xij_ptr = m->get_ptr_var<double>("Xij");
X_ptr->_psd_var=Xij_ptr;
Xij_ptr->_psd_var=X_ptr;


var<int> y("y",0,1);
m->add(y.in(node_pairs_chord));

Constraint<> def_Xij("Xij");
def_Xij=Xij-(-1.0/(num_part-1.0)+num_part*1.0/(num_part-1.0)*y);
m->add(def_Xij.in(node_pairs_chord)==0);
    
    for(auto k:*node_pairs_chord._keys){
        if(!node_pairs.has_key(k)){
            w.add_val(k,0.0);
        }
        map_y[k].push_back(make_pair(k, num_part*1.0/(num_part-1.0)));
        map_const[k]=(-1.0)/(num_part-1.0);
    }
    
    Constraint<> SOC("SOC");
    SOC = Xij*Xij - X.from(node_pairs_chord)*X.to(node_pairs_chord);
    SOC.add_to_callback();
    m->add(SOC.in(node_pairs_chord) <= 0);
    
    int count=0;
    std::vector<pair<int,std::vector<string>>> _bag_names;
    m->sdp_dual=false;
    for(auto b:g._bags){
        pair<int,vector<string>> bn;
        bn.first=count++;
        if(b.second.size()>=3){
            DebugOn("bag "<<count<<endl);
            for(auto n:b.second){
                bn.second.push_back(n->_name);
                DebugOn(n->_name <<"\t");
            }
            DebugOn(endl);
            _bag_names.push_back(bn);
        }
    }
  
func<> obj=w.tr()*y;
m->min(obj);

m->print();
    var<double> x("x", 0, 0);
    m->add(x.in(range(0,0)));
    
    m->map_x=map_x;
    m->map_y=map_y;
    m->map_const=map_const;
    m->_bag_names=_bag_names;
   
    
return g;

}
Net model_mink_dense(string fname, shared_ptr<Model<double>>& m)
{
Net g;
  int num_part=3;
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
    name = to_string(i + 1);
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
param<double> w("w");
  while (getline(infile, sLine, '\n')) {
    istringstream iss(sLine);
    iss >> src_a >> dest_a >> weight;
    // cout << src  << ", " << dest << ", " << weight << endl;
    
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
infile.close();
    
auto node_pairs=g.get_node_pairs();

    
map<string, vector<pair<string, double>>> map_x;
map<string, vector<pair<string, double>>> map_y;
map<string, double> map_const;




g.sdp_3d_cuts=false;
g.get_tree_decomp_bags();
m->sdp_dual=false;

//auto node_pairs_chord = g.get_node_pairs_chord;
auto node_pairs_chord = g.get_node_pairs_chord(g._bags);

var<> X("X", 1, 1);
X._psd_diag=true;
m->add(X.in(nodes));
var<> Xij("Xij", -1.0/(num_part-1.0), 1);
Xij._psd_off_diag=true;
m->add(Xij.in(node_pairs_chord));
auto X_ptr = m->get_ptr_var<double>("X");
auto Xij_ptr = m->get_ptr_var<double>("Xij");
X_ptr->_psd_var=Xij_ptr;
Xij_ptr->_psd_var=X_ptr;


var<int> y("y",0,1);
m->add(y.in(node_pairs_chord));

Constraint<> def_Xij("Xij");
def_Xij=Xij-(-1.0/(num_part-1.0)+num_part*1.0/(num_part-1.0)*y);
m->add(def_Xij.in(node_pairs_chord)==0);
    
    for(auto k:*node_pairs_chord._keys){
        if(!node_pairs.has_key(k)){
            w.add_val(k,0.0);
        }
        map_y[k].push_back(make_pair(k, num_part*1.0/(num_part-1.0)));
        map_const[k]=(-1.0)/(num_part-1.0);
    }
    
    Constraint<> SOC("SOC");
    SOC = Xij*Xij - X.from(node_pairs_chord)*X.to(node_pairs_chord);
    SOC.add_to_callback();
    m->add(SOC.in(node_pairs_chord) <= 0);
    
    int count=0;
    std::vector<pair<int,std::vector<string>>> _bag_names;
    m->sdp_dual=false;
    for(auto b:g._bags){
        pair<int,vector<string>> bn;
        bn.first=count++;
        if(b.second.size()>=3){
            DebugOn("bag "<<count<<endl);
            for(auto n:b.second){
                bn.second.push_back(n->_name);
                DebugOn(n->_name <<"\t");
            }
            DebugOn(endl);
            _bag_names.push_back(bn);
        }
    }
  
func<> obj=w.tr()*y;
m->min(obj);

m->print();
    var<double> x("x", 0, 0);
    m->add(x.in(range(0,0)));
    
    m->map_x=map_x;
    m->map_y=map_y;
    m->map_const=map_const;
    m->_bag_names=_bag_names;
   
    
return g;

}
Net model_mink_free(string fname, shared_ptr<Model<double>>& m) {
Net g;
  int num_part=3;
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
    
    DebugOn("Nodes "<<Num_nodes<<endl);
    DebugOn("Edges "<<Num_edges<<endl);

  Node* node = nullptr;
    indices nodes;

  for (int i = 0; i < Num_nodes; i++) {
    name = to_string(i + 1);
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
param<double> w("w");
  while (getline(infile, sLine, '\n')) {
    istringstream iss(sLine);
    iss >> src_a >> dest_a >> weight;
    // cout << src  << ", " << dest << ", " << weight << endl;
    
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
infile.close();
    
auto node_pairs=g.get_node_pairs();

    
map<string, vector<pair<string, double>>> map_x;
map<string, vector<pair<string, double>>> map_y;
map<string, double> map_const;




g.sdp_3d_cuts=false;
g.get_tree_decomp_bags();
m->sdp_dual=false;

//auto node_pairs_chord = g.get_node_pairs_chord;
auto node_pairs_chord = g.get_node_pairs_chord(g._bags);

var<> X("X", 1, 1);
X._psd=true;
m->add(X.in(nodes));
var<> Xij("Xij", -1.0/(num_part-1.0), 1);
Xij._psd=true;
m->add(Xij.in(node_pairs_chord));


var<int> y("y",0,1);
m->add(y.in(node_pairs));

Constraint<> def_Xij("Xij");
def_Xij=Xij-(-1.0/(num_part-1.0)+num_part*1.0/(num_part-1.0)*y);
m->add(def_Xij.in(node_pairs)==0);
    
    for(auto k:*node_pairs_chord._keys){
        if(node_pairs.has_key(k)){
            map_y[k].push_back(make_pair(k, num_part*1.0/(num_part-1.0)));
            map_const[k]=(-1.0)/(num_part-1.0);
        }
//        else{
//            w.add_val(k,0.0);
//        }
    }
    
    Constraint<> SOC("SOC");
    SOC = Xij*Xij - X.from(node_pairs_chord)*X.to(node_pairs_chord);
   // SOC.add_to_callback();
   // m->add(SOC.in(node_pairs_chord) <= 0);
    
    int count=0;
    std::vector<pair<int,std::vector<string>>> _bag_names;
    m->sdp_dual=false;
    for(auto b:g._bags){
        pair<int,vector<string>> bn;
        bn.first=count++;
        if(b.second.size()>=2){
            DebugOn("bag "<<count<<endl);
            for(auto n:b.second){
                bn.second.push_back(n->_name);
                DebugOn(n->_name <<"\t");
            }
            DebugOn(endl);
            _bag_names.push_back(bn);
        }
    }
  

func<> obj=w.tr()*y;
m->min(obj);

m->print();
    var<double> x("x");
    m->add(x.in(range(0,0)));
    
    m->map_x=map_x;
    m->map_y=map_y;
    m->map_const=map_const;
    m->_bag_names=_bag_names;
   
    
return g;

}


#endif /* read_misdp_hpp */
