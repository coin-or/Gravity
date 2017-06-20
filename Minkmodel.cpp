//
//  MinkMinkmodel.cpp
//  Gravity
//
//  Created by Guanglei Wang on 19/6/17.
//
//

#include "Minkmodel.hpp"
using namespace std;


template<typename T>
void Sort(T &a, T &b, T &c){
    if (a > b) {
        std::swap(a, b);
    }
    if (a > c) {
        std::swap(a, c);
    }
    if (b > c) {
        std::swap(b, c);
    }
}



Minkmodel::Minkmodel(){};
Minkmodel::~Minkmodel(){ _graph->~Net(); _model.~Model();};

Minkmodel::Minkmodel(ModelType type, Net* graph, double K):_type(type),_graph(graph),_K(K),_solver(cplex),zij("zij"),Xij("Xij"){};


Minkmodel::Minkmodel(ModelType type, Net* graph, double K,SolverType solver):_type(type),_graph(graph),_K(K),_solver(solver),zij("zij"),Xij("Xij"){};


void Minkmodel::build(){
    switch (_type) {
        case MIP:
            add_vars_origin();
            add_triangle();
            add_clique();
            break;
        case SDP:
            add_vars_lifted();
            add_triangle_lifted();
            add_clique_lifted();
            break;
        case MIP_tree:
            add_vars_origin();
            add_triangle_tree();
           // add_clique();
        case SDP_tree:
            add_vars_lifted();
            add_triangle_lifted_tree();
            add_clique_lifted();
        default:
            break;
    }
}
void Minkmodel::reset(){};

void Minkmodel::add_vars_origin(){
    _model.add_var(zij^(_graph->nodes.size()*(_graph->nodes.size()-1)/2));
    
    func_ obj_MIP;
    int i=0, j=0;
    for (auto a: _graph->arcs){
        i = (a->src)->ID;
        j = (a->dest)->ID;
        if (i <= j)
            obj_MIP += (a->weight)*zij(i,j);
        else
            obj_MIP += (a->weight)*zij(j,i);
    }
    _model.set_objective(min(obj_MIP));

    
}

void Minkmodel::add_vars_lifted(){
    Xij.add_bounds(-1/(_K-1), 1);
    _model.add_var(Xij^(_graph->nodes.size()*(_graph->nodes.size()-1)/2));
    
    int i=0, j=0;
    func_ obj;
    for (auto a: _graph->arcs){
        i = (a->src)->ID;
        j = (a->dest)->ID;
        if (i <= j)
            obj += a->weight*((_K-1)*Xij(i,j) + 1)/_K;
        else
            obj += a->weight*((_K-1)*Xij(j,i) + 1)/_K;
    }
    _model.set_objective(min(obj));
}

void Minkmodel::add_triangle(){
    auto n = _graph->nodes.size();
    for (auto i=0; i<n-1; i++)
        for (auto h=i+1; h<n; h++)
            for (auto j=h+1; j<n; j++){
                Constraint Triangle1("Triangle1("+to_string(i)+","+to_string(h)+ ","+to_string(j)+")");
                Triangle1 = zij(i,h)+zij(h,j)-zij(i,j);
                Constraint Triangle2("Triangle2("+to_string(i)+","+to_string(h)+ ","+to_string(j)+")");
                Triangle2 = zij(i,h)+zij(i,j)-zij(h,j);
                Constraint Triangle3("Triangle3("+to_string(i)+","+to_string(h)+ ","+to_string(j)+")");
                Triangle3 = zij(i,j)+zij(h,j)- zij(i,h);
                _model.add_constraint(Triangle1 <=1);
                _model.add_constraint(Triangle2 <=1);
                _model.add_constraint(Triangle3 <=1);
                }
}

void Minkmodel::add_triangle_lifted(){
    auto n = _graph->nodes.size();
    for (auto i=0; i<n; i++)
        for (auto h=i+1; h<n; h++)
            for (auto j=h+1; j<n;j++){
                Constraint Triangle1("Triangle1("+to_string(i)+","+to_string(h)+ ","+to_string(j)+")");
                Triangle1 = Xij(i,h)+Xij(h,j)-Xij(i,j);
                Constraint Triangle2("Triangle2("+to_string(i)+","+to_string(h)+ ","+to_string(j)+")");
                Triangle2 = Xij(i,h)+Xij(i,j)-Xij(h,j);
                Constraint Triangle3("Triangle3("+to_string(i)+","+to_string(h)+ ","+to_string(j)+")");
                Triangle3 = Xij(i,j)+Xij(h,j)- Xij(i,h);
                _model.add_constraint(Triangle1<=1);
                _model.add_constraint(Triangle2<=1);
                _model.add_constraint(Triangle3<=1);
                }
}

void Minkmodel::add_clique(){
    auto n = _graph->nodes.size();
    if (_K >2) {
        for (auto i=0; i<n-1; i++)
            for (auto h=i+1; h<n; h++)
                for (auto j=h+1; j<n;j++)
                    for (auto l=j+1;l<n;l++)
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
                for (auto j=h+1; j<n;j++)
                {
                    Constraint Clique("Clique("+to_string(i)+","+to_string(h)+ ","+to_string(j)+")");
                    Clique = zij(i,h) +zij(i,j) + zij(h,j);
                    _model.add_constraint(Clique >=1);
                }
        }
}

void Minkmodel::add_clique_lifted(){
    auto n = _graph->nodes.size();
    if (_K >2) {
        for (auto i=0; i<n; i++)
            for (auto h=i+1; h<n; h++)
                for (auto j=h+1; j<n;j++)
                    for (auto l=j+1;l<n;l++)
                    {
                    Constraint Clique("Clique("+to_string(i)+","+to_string(h)+ ","+to_string(j)+ ", "+to_string(l)+")");
                    Clique = Xij(i,h) +Xij(i,j) + Xij(i,l) + Xij(h,j) + Xij(h,l) +Xij(j,l);
                    _model.add_constraint(Clique >=-0.5*_K);
                    }
            }
        else{
            for (auto i=0; i<n-1; i++)
                for (auto h=i+1; h<n; h++)
                    for (auto j=h+1; j<n;j++){
                        Constraint Clique("Clique("+to_string(i)+","+to_string(h)+ ","+to_string(j)+")");
                        Clique = Xij(i,h) + Xij(h,j) + Xij(i,j);
                        _model.add_constraint(Clique >=-0.5*_K);
                    }
        }
}

void Minkmodel::tree_decompose(){
    _graph->get_tree_decomp_bags();
    int i1,i2,i3;
    for (int i = 0; i < _graph->_bags.size(); i++){
        auto bag = _graph->_bags.at(i);
        if (bag.size()<3) {
            continue;
        }
        
        for (int j = 0; j < bag.size()-2;j++)
            for (int h=j+1; h<bag.size()-1;h++)
                for (int l=j+1; l<bag.size();l++){
                    i1 = bag[j]->ID;
                    i2 = bag[h]->ID;
                    i3 = bag[l]->ID;
                    if (i1 != i2 && i2 !=i3 &&i3!= i1){
                        Sort(i1,i2,i3);
                        if(_ids.count(make_tuple(i1, i2, i3))==0){
                            _ids.insert(make_tuple(i1, i2, i3));
                        }
                        else {
                            continue;
                        }
                    }
                }
    }
    cout << "size of ids: " << _ids.size() << endl;
}

void Minkmodel::add_3Dcuts(){
    int i1,i2,i3;
    for (auto it: _ids){
        i1 = get<0>(it);
        i2 = get<1>(it);
        i3 = get<2>(it);
        Constraint SDP3("SDP3("+to_string(i1)+","+to_string(i2)+","+to_string(i3)+")");
        SDP3 = -2*Xij(i1,i2)*Xij(i2,i3)*Xij(i1,i3);
        SDP3 -= 1;
        SDP3 += power(Xij(i1,i2),2);
        SDP3 += power(Xij(i1,i3),2);
        SDP3 += power(Xij(i2,i3),2);
        _model.add_constraint(SDP3);
    }
};

void Minkmodel::add_triangle_tree(){
    int i1,i2,i3;
    for (auto it: _ids){
        i1 = get<0>(it);
        i2 = get<1>(it);
        i3 = get<2>(it);
        cout << "(i1, i2, i3): " << i1 << ", " << i2 << ", " << i3 << endl;
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

void Minkmodel::add_triangle_lifted_tree(){
    int i1,i2,i3;
    for (auto it: _ids){
        i1 = get<0>(it);
        i2 = get<1>(it);
        i3 = get<2>(it);
        cout << "(i1, i2, i3): " << i1 << ", " << i2 <<", " << i3 <<endl;
        Constraint Triangle1("Triangle1("+to_string(i1)+","+to_string(i2)+ ","+to_string(i3)+")");
        Triangle1 = Xij(i1,i2)+Xij(i1,i3)-Xij(i2,i3);
        Constraint Triangle2("Triangle2("+to_string(i1)+","+to_string(i2)+ ","+to_string(i3)+")");
        Triangle2 = Xij(i1,i3)+Xij(i2,i3)-Xij(i1,i2);
        Constraint Triangle3("Triangle3("+to_string(i1)+","+to_string(i2)+ ","+to_string(i3)+")");
        Triangle3 = Xij(i2,i3)+Xij(i1,i2)- Xij(i1,i3);
        _model.add_constraint(Triangle1<=1);
        _model.add_constraint(Triangle2<=1);
        _model.add_constraint(Triangle3<=1);
    }

}

int Minkmodel::solve(){
    solver s(_model,_solver);
    //double wall0 = get_wall_time();
    //double cpu0  = get_cpu_time();
    cout << "Running the relaxation model \n";
    s.run();
    //double wall1 = get_wall_time();
    //double cpu1  = get_cpu_time();
    // relax.print_solution();
    return 1;
}
