//
//  MinkMinkmodel.cpp
//  Gravity
//
//  Created by Guanglei Wang on 19/6/17.
//
//

#include "Minkmodel.hpp"
//#define USEDEBUG
#ifdef USEDEBUG
#define Debug(x) cout << x
#else
#define Debug(x)
#endif
#define DebugOn(x) cout << x
#define DebugOff(x)

using namespace std;

//template<typename T>
//void Sort(T &a, T &b, T &c){
//    if (a > b) {
//        std::swap(a, b);
//    }
//    if (a > c) {
//        std::swap(a, c);
//    }
//    if (b > c) {
//        std::swap(b, c);
//    }
//}



Minkmodel::Minkmodel() {};

Minkmodel::~Minkmodel() {};

//Minkmodel::Minkmodel(ModelType type, Net* graph, double K):_type(type),_solver(cplex),_K(K),_graph(graph),zij("zij"),X("X"){};
Minkmodel::Minkmodel(ModelType type, Net* graph, double K):_type(type),_solver(cplex),_K(K),_graph(graph) {
    _cliqueid = make_shared<map<string,vector<unsigned>>>();
};

//Minkmodel::Minkmodel(ModelType type, Net* graph, double K,SolverType solver):_type(type),_solver(solver),_K(K),_graph(graph),zij("zij"),X("X"){};
Minkmodel::Minkmodel(ModelType type, Net* graph, double K,SolverType solver):_type(type),_solver(solver),_K(K),_graph(graph) {
    _cliqueid = make_shared<map<string,vector<unsigned>>>();
};

void Minkmodel::build() {
    switch (_type) {
    case MIP:
        add_vars_origin();
        add_triangle();
        add_clique();
        break;
    case SDP:
        add_vars_lifted();
        //add_triangle_lifted();
        //add_clique_lifted();
        break;
    case MIP_tree:
        add_vars_origin_tree();
        cliquetree_decompose();
        add_triangle_tree();
        add_clique_tree();
        break;
    case SDP_tree:
        add_vars_lifted_tree();
        //cliquetree_decompose();
        //add_triangle_lifted_tree();
        //add_clique_lifted_tree();
        break;
    default:
        break;
    }
}
void Minkmodel::reset() {};
void Minkmodel::add_vars_origin() {
    var<bool> zij("zij");
    _model.add_var(zij^(_graph->nodes.size()*(_graph->nodes.size()-1)/2));

    func_ obj_MIP;
    int i=0, j=0;
    for (auto a: _graph->arcs) {
        i = (a->src)->ID;
        j = (a->dest)->ID;
        if (i <= j)
            obj_MIP += (a->weight)*zij(i,j);
        else
            obj_MIP += (a->weight)*zij(j,i);
    }
    _model.set_objective(min(obj_MIP));
}

void Minkmodel::add_vars_origin_tree() {
    var<bool> zij("zij");
    // the number of arcs in the chordal extension
    _model.add_var(zij^((_graph->_chordalextension)->arcs.size()));//
    //cout << _graph->arcs.size() << endl;

    func_ obj_MIP;
    int i=0, j=0;
    for (auto a: _graph->arcs) {
        i = (a->src)->ID;
        j = (a->dest)->ID;
        if (i <= j)
            obj_MIP += (a->weight)*zij(i,j);
        else
            obj_MIP += (a->weight)*zij(j,i);
    }
    _model.set_objective(min(obj_MIP));
}

void Minkmodel::add_vars_lifted() {
   sdpvar<double> X("X");
    _model.add_var(X^(_graph->nodes.size())); // S^n.

    int i=0, j=0;
    func_ obj;
    for (auto a: _graph->arcs) {
        i = (a->src)->ID;
        j = (a->dest)->ID;
        if (i <= j)
            obj += a->weight*((_K-1)*X(i,j) + 1)/_K;
        else
            obj += a->weight*((_K-1)*X(j,i) + 1)/_K;
    }
    _model.set_objective(min(obj));
    
    for (i = 0; i < _graph->nodes.size(); i++){
        Constraint diag("("+to_string(i)+")");
        diag = X(i,i)-1;
        _model.add_constraint(diag=0);
    }
    
    for (i= 0; i < _graph->nodes.size(); i ++)
        for (j = i+1; j< _graph->nodes.size(); j++){
            Constraint bound("("+to_string(i)+ "," + to_string(j) +")");
            bound = X(i,j) + 1/(_K-1);
            _model.add_constraint(bound >=0);
        }
}

void Minkmodel::add_vars_lifted_tree() {
    sdpvar<double> X("X");
    _model.add_var(X^(_graph->nodes.size()));
    
    int i=0, j=0;
    func_ obj;
    for (auto a: _graph->arcs) {
        i = (a->src)->ID;
        j = (a->dest)->ID;
        if (i <= j)
            obj += a->weight*((_K-1)*X(i,j) + 1)/_K;
        else
            obj += a->weight*((_K-1)*X(j,i) + 1)/_K;
    }
    _model.set_objective(min(obj));
    
    for (i = 0; i < _graph->nodes.size(); i++){
        Constraint diag("("+to_string(i)+")");
        diag = X(i,i)-1;
        _model.add_constraint(diag=0);
    }
    
    for (auto a: _graph->_chordalextension->arcs){
        i = (a->src)->ID;
        j = (a->dest)->ID;
        Constraint bound("("+to_string(i)+ "," + to_string(j) +")");
        bound = X(i,j) + 1/(_K-1);
        _model.add_constraint(bound >=0);

    }
}

void Minkmodel::add_triangle() {
    auto n = _graph->nodes.size();
    auto zij = (*(var<bool>*)(_model.get_var("zij")));

    for (auto i=0; i<n; i++)
        for (auto h=i+1; h<n; h++)
            for (auto j=h+1; j<n; j++) {
                Constraint Triangle1("Triangle1("+to_string(i)+","+to_string(h)+ ","+to_string(j)+")");
                Triangle1 = zij(i,h)+ zij(h,j)-zij(i,j);
                Constraint Triangle2("Triangle2("+to_string(i)+","+to_string(h)+ ","+to_string(j)+")");
                Triangle2 = zij(i,h)+zij(i,j)-zij(h,j);
                Constraint Triangle3("Triangle3("+to_string(i)+","+to_string(h)+ ","+to_string(j)+")");
                Triangle3 = zij(i,j)+zij(h,j)- zij(i,h);
                _model.add_constraint(Triangle1 <=1);
                _model.add_constraint(Triangle2 <=1);
                _model.add_constraint(Triangle3 <=1);
            }
}

void Minkmodel::add_triangle_lifted() {
    auto n = _graph->nodes.size();
    auto X = (*(var<bool>*)(_model.get_var("X")));

    for (auto i=0; i<n; i++)
        for (auto h=i+1; h<n; h++)
            for (auto j=h+1; j<n; j++) {
                Constraint Triangle1("Triangle1("+to_string(i)+","+to_string(h)+ ","+to_string(j)+")");
                Triangle1 = X(i,h)+X(h,j)-X(i,j);
                Constraint Triangle2("Triangle2("+to_string(i)+","+to_string(h)+ ","+to_string(j)+")");
                Triangle2 = X(i,h)+X(i,j)-X(h,j);
                Constraint Triangle3("Triangle3("+to_string(i)+","+to_string(h)+ ","+to_string(j)+")");
                Triangle3 = X(i,j)+X(h,j)- X(i,h);
                _model.add_constraint(Triangle1<=1);
                _model.add_constraint(Triangle2<=1);
                _model.add_constraint(Triangle3<=1);
            }
}

void Minkmodel::add_clique() {
    auto n = _graph->nodes.size();
    auto zij = (*(var<bool>*)(_model.get_var("zij")));
    if (_K >2) {
        for (auto i=0; i<n; i++)
            for (auto h=i+1; h<n; h++)
                for (auto j=h+1; j<n; j++)
                    for (auto l=j+1; l<n; l++)
                    {
                        Constraint Clique("ZClique("+to_string(i)+","+to_string(h)+ ","+to_string(j)+ ", "+to_string(l)+")");
                        Clique = zij(i,h) +zij(i,j) + zij(i,l) + zij(h,j) + zij(h,l) +zij(j,l);
                        _model.add_constraint(Clique >=1);
                    }
    }
    else
    {
        for (auto i=0; i<n-1; i++)
            for (auto h=i+1; h<n; h++)
                for (auto j=h+1; j<n; j++)
                {
                    Constraint Clique("Clique("+to_string(i)+","+to_string(h)+ ","+to_string(j)+")");
                    Clique = zij(i,h) +zij(i,j) + zij(h,j);
                    _model.add_constraint(Clique >=1);
                }
    }
}

void Minkmodel::add_clique_lifted() {
    auto n = _graph->nodes.size();
    auto X = (*(var<bool>*)(_model.get_var("X")));
    if (_K >2) {
        for (auto i=0; i<n; i++)
            for (auto h=i+1; h<n; h++)
                for (auto j=h+1; j<n; j++)
                    for (auto l=j+1; l<n; l++)
                    {
                        Constraint Clique("Clique("+to_string(i)+","+to_string(h)+ ","+to_string(j)+ ", "+to_string(l)+")");
                        Clique = X(i,h) +X(i,j) + X(i,l) + X(h,j) + X(h,l) +X(j,l);
                        _model.add_constraint(Clique >=-0.5*_K);
                    }
    }
    else {
        for (auto i=0; i<n-1; i++)
            for (auto h=i+1; h<n; h++)
                for (auto j=h+1; j<n; j++) {
                    Constraint Clique("Clique("+to_string(i)+","+to_string(h)+ ","+to_string(j)+")");
                    Clique = X(i,h) + X(h,j) + X(i,j);
                    _model.add_constraint(Clique >=-0.5*_K);
                }
    }
}

// construct a temporary container;
vector<unsigned> temp;


void Minkmodel::nchoosek(int bag_id, int offset, int K){
    if (K == 0) {
        static int count = 0;
        cout << "combination no " << (++count) << ": [ ";
        for (int i = 0; i < temp.size(); ++i) { cout << temp[i] << " "; }
        cout << "] " << endl;
        std::string key;
        for (unsigned i = 0; i < temp.size(); ++i) {
            key += to_string(temp[i]);
            if(i<temp.size()-1) {
                key+=",";
            }
        }
        
        auto iter = _cliqueid->find(key);
        if (iter == _cliqueid->end()) {
            _cliqueid->insert(make_pair<>(key,temp));
        }
        return;
    }
    auto bag= _graph->_bags.at(bag_id);
    for (int i = offset; i <= bag.size() - K; i++) {
        temp.push_back(bag[i]->ID);
        //cout << "bag_id: " << bag_id << " i+1: " << i+1 << " K-1 " << K-1 << endl;
        nchoosek(bag_id, i+1, K-1);
        temp.pop_back();
    }
}


void Minkmodel::cliquetree_decompose() {
    //_graph->get_tree_decomp_bags();
    int i1,i2,i3,i4;
    int j1,j2,j3,j4;
    for (int i = 0; i < _graph->_bags.size(); i++) {
        auto bag = _graph->_bags.at(i);
        if (bag.size()<3) {
            continue;
        }
        for (int j = 0; j < bag.size()-2; j++)
            for (int h=j+1; h<bag.size()-1; h++)
                for (int l=h+1; l<bag.size(); l++) {
                    i1 = bag[j]->ID; // zero index.
                    i2 = bag[h]->ID;
                    i3 = bag[l]->ID;
                    //cout << "(i1, i2, i3) " << i1 << " " << i2 << " " << i3 << endl;
                    // Sort(i1,i2,i3);
                    if(_ids.count(make_tuple(i1, i2, i3))==0) {
                        _ids.insert(make_tuple(i1, i2, i3));
                    }
                    else {
                        continue;
                    }
                }
        if (bag.size()> _K+1 && _K > 2) {
            nchoosek(i,0,_K+1);
            for(j1 = 0; j1 < bag.size()-3; j1++)
                for(j2=j1+1; j2<bag.size()-2; j2++)
                    for(j3=j2+1; j3<bag.size()-1; j3++)
                        for(j4=j3+1; j4<bag.size(); j4++) {
                            i1 = bag[j1]->ID;
                            i2 = bag[j2]->ID;
                            i3 = bag[j3]->ID;
                            i4 = bag[j4]->ID;
                // cout << "(i1, i2, i3, i4) " << i1 << " " << i2 << " " << i3 <<" " << i4<< endl; Sort(i1,i2,i3,i4);
                            if(_ids4.count(make_tuple(i1, i2, i3,i4))==0) {
                                _ids4.insert(make_tuple(i1, i2,i3,i4));
                            }
                            else {
                                continue;
                            }
                        }
        }
    }
    cout << "size of triangle inequalties: " << _ids.size() << endl;
    cout << "size of clique inequalties: " << _cliqueid->size() << endl;
    cout << "size of clique inequalties_ids4 " << _ids4.size() << endl;

}

void Minkmodel::add_3Dcuts() {
    auto X = (*(var<bool>*)(_model.get_var("X")));
    int i1,i2,i3;
    for (auto it: _ids) {
        i1 = get<0>(it);
        i2 = get<1>(it);
        i3 = get<2>(it);
        Constraint SDP3("SDP3("+to_string(i1)+","+to_string(i2)+","+to_string(i3)+")");
        SDP3 = -2*X(i1,i2)*X(i2,i3)*X(i1,i3);
        SDP3 -= 1;
        SDP3 += power(X(i1,i2),2);
        SDP3 += power(X(i1,i3),2);
        SDP3 += power(X(i2,i3),2);
        _model.add_constraint(SDP3);
    }
};

void Minkmodel::add_triangle_tree() {
    int i1,i2,i3;
    auto zij = (*(var<bool>*)(_model.get_var("zij")));
    for (auto it: _ids) {
        i1 = get<0>(it);
        i2 = get<1>(it);
        i3 = get<2>(it);
        //cout << "(i1, i2, i3): " << i1 << ", " << i2 << ", " << i3 << endl;
        Constraint Triangle1("ZTriangle1("+to_string(i1)+","+to_string(i2)+ ","+to_string(i3)+")");
        Triangle1 = zij(i1,i2)+zij(i1,i3)-zij(i2,i3);
        Constraint Triangle2("ZTriangle2("+to_string(i1)+","+to_string(i2)+ ","+to_string(i3)+")");
        Triangle2 = zij(i1,i3)+zij(i2,i3)-zij(i1,i2);
        Constraint Triangle3("ZTriangle3("+to_string(i1)+","+to_string(i2)+ ","+to_string(i3)+")");
        Triangle3 = zij(i2,i3)+zij(i1,i2)- zij(i1,i3);
        _model.add_constraint(Triangle1 <=1);
        _model.add_constraint(Triangle2 <=1);
        _model.add_constraint(Triangle3 <=1);
    }
}

void Minkmodel::add_triangle_lifted_tree() {
    int i1,i2,i3;
    auto X = (*(var<bool>*)(_model.get_var("X")));
    for (auto it: _ids) {
        i1 = get<0>(it);
        i2 = get<1>(it);
        i3 = get<2>(it);
        //cout << "(i1, i2, i3): " << i1 << ", " << i2 <<", " << i3 <<endl;
        Constraint Triangle1("Triangle1("+to_string(i1)+","+to_string(i2)+ ","+to_string(i3)+")");
        Triangle1 = X(i1,i2)+X(i1,i3)-X(i2,i3);
        Constraint Triangle2("Triangle2("+to_string(i1)+","+to_string(i2)+ ","+to_string(i3)+")");
        Triangle2 = X(i1,i3)+X(i2,i3)-X(i1,i2);
        Constraint Triangle3("Triangle3("+to_string(i1)+","+to_string(i2)+ ","+to_string(i3)+")");
        Triangle3 = X(i2,i3)+X(i1,i2)- X(i1,i3);
        _model.add_constraint(Triangle1<=1);
        _model.add_constraint(Triangle2<=1);
        _model.add_constraint(Triangle3<=1);
    }

}

void Minkmodel::add_clique_tree() {
//    int i1,i2,i3,i4;
   auto zij = (*(var<bool>*)(_model.get_var("zij")));
    if (_K>2){
        for (auto it: (*_cliqueid)) {
            auto key = it.first;
            auto value = it.second;
            Constraint Clique("ZClique["+ key +"]");
            for (int i = 0; i < value.size()-1; i++){
                auto id1 = value[i];
                for (int j = i+1 ; j< value.size(); j++){
                    auto id2 = value[j];
                    if (id1 <= id2)
                        Clique += zij(id1,id2);
                    else
                        Clique += zij(id2,id1);
                }
            }
           _model.add_constraint(Clique >=1);
            //Clique.print();
        }
    }
    else {
        for (auto it: _ids) {
            auto i1 = get<0>(it);
            auto i2 = get<1>(it);
            auto i3 = get<2>(it);
            //cout << "(i1, i2, i3): " << i1 << ", " << i2 << ", " << i3<< endl;
            Constraint Clique("ZClique("+to_string(i1)+","+to_string(i2)+ ","+to_string(i3)+")");
            Clique = zij(i1,i2) +zij(i1,i3) + zij(i2,i3);
            _model.add_constraint(Clique >=1);
        }
    }
}

void Minkmodel::add_general_clique() {
//  t=1, q=2,...,k-1
    auto zij = (*(var<bool>*)(_model.get_var("zij")));
    for (int i = 0; i < _graph->_bags.size(); i++) {
        auto bag = _graph->_bags.at(i);
        if (bag.size() > _K+1) {
            Constraint GClique("General_Clique"+to_string(i));
            int t = floor(bag.size()/_K);
            int q = (bag.size()-t*_K);
            for (int j= 0; j< bag.size()-1; j++) {
                int idj=bag[j]->ID;
                for (int k=j+1; k< bag.size(); k++) {
                    int idk=bag[k]->ID;
                    if (idj < idk)
                        GClique +=zij(idj,idk);
                    else
                        GClique +=zij(idk,idj);

                }
            }
            _model.add_constraint(GClique >= 0.5*t*(t-1)*(t-q)+0.5*t*(t+1)*q);
        }
    }
}

void Minkmodel::add_wheel() {
// auto zij = (*(var<bool>*)(_model.get_var("zij")));
//    int q = 5;
//    for (int i = 0; i < _graph->_bags.size(); i++){
//        auto bag = _graph->_bags.at(i);
//        if (bag.size()>q+1) {
//            for (int j= 0; j< bag.size(); j++){
//                Constraint GClique("General_Clique: bag "+to_string(i)+" constr: "+to_string(j));
//                 // wheel center
//                int idfirst=bag[j]->ID;
//                // the spokes part
//                for (int k=0; k< j; k++){
//                    int idsec=bag[k]->ID;
//                    if (idfirst < idsec)
//                        GClique += zij(idfirst,idsec);
//                    else
//                        GClique += zij(idsec,idfirst);
//
//                }
//                for (int k=j+1; k< bag.size(); k++){
//                    int idsec=bag[k]->ID;
//                    GClique += zij(idfirst,idsec);
//                }
//                // the wheel part
//                // generate the wheel set
//                vector<int> wheel;
//                for (int k=0; k < bag.size(); k++)
//                    if (k!=j) wheel.push_back(k);
//
//                for (auto iter: wheel){
//
//
//                }
//
//
//            }
//        }
//    }
}

void Minkmodel::add_bicycle() {}

//void Minkmodel::add_clique_lifted_tree() {
//    auto X = (*(var<bool>*)(_model.get_var("X")));
//    int i1,i2,i3,i4;
//    if (_K>2)
//        for (auto it: _ids4) {
//            i1 = get<0>(it);
//            i2 = get<1>(it);
//            i3 = get<2>(it);
//            i4 = get<3> (it);
//            //cout << "(i1, i2, i3, i4): " << i1 << ", " << i2 << ", " << i3 << ", " << i4<< endl;
//            Constraint Clique("XClique("+to_string(i1)+","+to_string(i2)+ ","+to_string(i3)+ ", "+to_string(i4)+")");
//            Clique = X(i1,i2) +X(i1,i3) + X(i1,i4) + X(i2,i3) + X(i2,i4) +X(i3,i4);
//            _model.add_constraint(Clique >= -0.5*_K);
//        }
//    else {
//        for (auto it: _ids) {
//            i1 = get<0>(it);
//            i2 = get<1>(it);
//            i3 = get<2>(it);
//            //cout << "(i1, i2, i3): " << i1 << ", " << i2 << ", " << i3<< endl;
//            Constraint Clique("XClique("+to_string(i1)+","+to_string(i2)+ ","+to_string(i3)+")");
//            Clique = X(i1,i2) +X(i1,i3) + X(i2,i3);
//            _model.add_constraint(Clique >= -0.5*_K);
//        }
//    }
//}

int Minkmodel::solve(int output, bool relax) {
    solver s(_model,_solver);
    cout << "Running the relaxation model \n";
    s.run(output,relax);
    return 1;
}

void Minkmodel::construct_fsol() {
// construct a feasible solution from z
    auto zij = (*(var<bool>*)(_model.get_var("zij")));
    param<int> sol("sol");
    int i=0,j=0;
    if (_type==MIP) {
        cout << "The MIP solution is: " << endl;
        for (i= 0; i < _graph->nodes.size()-1; i++) {
            for (j= i+1; j < _graph->nodes.size(); j++) {
                sol(i,j)=zij(i,j).getvalue();
                cout << sol(i,j).to_str() << endl;
            }
        }
    }
   else{
        cout << "The constructed solution is: " << endl;
         for (auto a: _graph->_chordalextension->arcs) {
             i = (a->src)->ID;
             j = (a->dest)->ID;
             if (i <= j) {
                 sol(i,j)=zij(i,j).getvalue();
             }
             else {
                 // generally, will never enter
                 cerr << "something wrong with lables of chordal extension graph";
                 exit(1);
                 sol(j,i)=zij(j,i).getvalue();
             }
         }

         Node* n = nullptr;
         Node* nn = nullptr;
         Arc* arc_chordal=nullptr;
         bool allzeros=true;
         double temp=0;

         for (i=0; i < _graph->_chordalextension->nodes.size()-1; i++) {
             n=(_graph->_chordalextension->get_node(to_string(i+1)));
             for (j = i+1; j< _graph->_chordalextension->nodes.size(); j++) {
                 nn=_graph->_chordalextension->get_node(to_string(j+1));
                 //cout<< "(n, nn): " << n->ID << ", " << nn->ID << endl;
                 if (n->is_connected(nn)) {
                     cout << sol(i,j).to_str() << endl;
                 }
                 else
                 {
                     string name = to_string(_graph->_chordalextension->arcs.size()+1);
                     arc_chordal = new Arc(name);
                     arc_chordal->id = _graph->_chordalextension->arcs.size();
                     arc_chordal->src = n;
                     arc_chordal->dest = nn;
                     arc_chordal->weight = 0;
                     arc_chordal->connect();
                     // now find the maximal clique containing edge (i,j)

                     // for neigbour i intersect neighbour j
                     set<int> idi;
                     set<int> idj;
                     Debug("neighbours of " << i << ": ");

                     for (auto a: n->branches) {
                         idi.insert(a->neighbour(n)->ID);
                         Debug(a->neighbour(n)->ID << ", ");
                     }
                     Debug(endl << "neighbours of " << j << ": ");

                     for (auto a: nn->branches) {
                         idj.insert(a->neighbour(nn)->ID);
                         Debug(a->neighbour(nn)->ID << ", ");
                     }
                     Debug(endl);

                     // intersection of idi and idj.
                     std::vector<int> inter;
                     set_intersection(idi.begin(), idi.end(), idj.begin(),idj.end(),std::back_inserter(inter));

                     // how to get values of zij
                     if (inter.size()==0) {
                         sol(i,j)=1;
                     }
                     else {
                         allzeros=true;
                         for (auto h: inter) {
                             temp=0;
                             Debug("(i, j, h): " << i << " " << j << " " << h <<endl);
                             if (h <= i) {
                                 temp += sol(h,i).getvalue();
                                 // cout << sol(h,i).getvalue() << endl;
                             }
                             else {
                                 temp += sol(i,h).getvalue();
                                 //cout << sol(i, h).getvalue() << endl;
                             }
                             if (h<=j) {
                                 temp += sol(h,j).getvalue();
                                 //cout << sol(h,j).getvalue() << endl;
                             }
                             else {
                                 temp += sol(j,h).getvalue();
                                 //cout << sol(j,h).getvalue() << endl;
                             }

                             //cout << "temp: " << temp << endl;
                             //DebugOn("xhi + xhj = " << temp);

                             if (temp==2)
                             {
                                 sol(i,j)=1;
                                 allzeros=false;
                                 break;
                             }
                             else if (temp==1)
                             {
                                 sol(i,j)=0;
                                 allzeros=false;
                                 break;
                             }
                             else
                                 continue;
                         }
                         if (allzeros) {
                             // if inter.size() = _K-1, then the clique inequality is not implied in the chordal graph.. we need to enfore..
                             // if inter.size() >= _K, then some of the clique inequality has been eforced by intersect + i/j, but
                             //if (inter.size()>=_K -1)
                             sol(i,j) = 1;
                             //else
                             //  sol(i,j) = 0;
                         }
                     }
                     cout << sol(i,j).to_str() << endl;
                     (_graph->_chordalextension)->add_arc(arc_chordal);
                 }
             }
         }
        }
    }


