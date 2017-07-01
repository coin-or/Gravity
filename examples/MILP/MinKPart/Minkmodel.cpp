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


Minkmodel::Minkmodel(){
};

Minkmodel::~Minkmodel(){
};

Minkmodel::Minkmodel(ModelType type, Net* graph, double K):_type(type),_solver(cplex),_K(K),_graph(graph),zij("zij"),Xij("Xij"){};


Minkmodel::Minkmodel(ModelType type, Net* graph, double K,SolverType solver):_type(type),_solver(solver),_K(K),_graph(graph),zij("zij"),Xij("Xij"){};


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
            add_vars_origin_tree();
            tree_decompose();
            add_triangle_tree();
            add_clique_tree();
            //add_clique();
            break;
        case SDP_tree:
            add_vars_lifted();
            tree_decompose();
            add_triangle_lifted_tree();
            //add_clique_lifted();
            add_clique_lifted_tree();
            break;
        default:
            break;
    }
}
void Minkmodel::reset(){};

void Minkmodel::add_vars_origin(){
    zij.add_bounds(0,1);
    _model.add_var(zij^(_graph->nodes.size()*(_graph->nodes.size()-1)/2));
    
    func_ obj_MIP;
    int i=0, j=0;
    for (auto a: _graph->arcs){
        i = (a->src)->ID;
        j = (a->dest)->ID;
        if (i <= j &&a->weight !=0)
            obj_MIP += (a->weight)*zij(i,j);
        else
            obj_MIP += (a->weight)*zij(j,i);
    }
    _model.set_objective(min(obj_MIP));
}

void Minkmodel::add_vars_origin_tree(){
    zij.add_bounds(0,1);
    
    _model.add_var(zij^(_graph->arcs.size()));
    
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
    //_graph->get_tree_decomp_bags();
    int i1,i2,i3,i4;
    int j1,j2,j3,j4;
    for (int i = 0; i < _graph->_bags.size(); i++){
        auto bag = _graph->_bags.at(i);
        if (bag.size()<3) {
            continue;
        }
        for (int j = 0; j < bag.size()-2;j++)
            for (int h=j+1; h<bag.size()-1;h++)
                for (int l=h+1; l<bag.size();l++){
                    i1 = bag[j]->ID;
                    i2 = bag[h]->ID;
                    i3 = bag[l]->ID;
                    //cout << "(i1, i2, i3) " << i1 << " " << i2 << " " << i3 << endl;
                       // Sort(i1,i2,i3);
                        if(_ids.count(make_tuple(i1, i2, i3))==0){
                            _ids.insert(make_tuple(i1, i2, i3));
                        }
                        else {
                            continue;
                        }
                }
        if (bag.size()>3){
            for(j1 = 0; j1 < bag.size()-3;j1++)
                for(j2=j1+1; j2<bag.size()-2;j2++)
                    for(j3=j2+1; j3<bag.size()-1;j3++)
                        for(j4=j3+1;j4<bag.size();j4++){
                            i1 = bag[j1]->ID;
                            i2 = bag[j2]->ID;
                            i3 = bag[j3]->ID;
                            i4 = bag[j4]->ID;
                       // cout << "(i1, i2, i3, i4) " << i1 << " " << i2 << " " << i3 <<" " << i4<< endl;
                         //   Sort(i1,i2,i3,i4);
                            if(_ids4.count(make_tuple(i1, i2, i3,i4))==0){
                                _ids4.insert(make_tuple(i1, i2,i3,i4));
                            }
                            else {
                                continue;
                            }
                    }
        }
    }
    cout << "size of 3d ids: " << _ids.size() << endl;
    cout << "size of 4d ids: " << _ids4.size() << endl;
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
    Constraint redudant("redudant");
    redudant = zij(2,3);
    _model.add_constraint(redudant <= 1);    
}

void Minkmodel::add_triangle_lifted_tree(){
    int i1,i2,i3;
    for (auto it: _ids){
        i1 = get<0>(it);
        i2 = get<1>(it);
        i3 = get<2>(it);
        //cout << "(i1, i2, i3): " << i1 << ", " << i2 <<", " << i3 <<endl;
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
void Minkmodel::add_clique_tree(){
    int i1,i2,i3,i4;
    if (_K>2)
        for (auto it: _ids4){
            i1 = get<0>(it);
            i2 = get<1>(it);
            i3 = get<2>(it);
            i4 = get<3> (it);
            //cout << "(i1, i2, i3, i4): " << i1 << ", " << i2 << ", " << i3 << ", " << i4<< endl;
            Constraint Clique("ZClique("+to_string(i1)+","+to_string(i2)+ ","+to_string(i3)+ ", "+to_string(i4)+")");
            Clique = zij(i1,i2) +zij(i1,i3) + zij(i1,i4) + zij(i2,i3) + zij(i2,i4) +zij(i3,i4);
            _model.add_constraint(Clique >=1);
        }
    else{
        for (auto it: _ids){
            i1 = get<0>(it);
            i2 = get<1>(it);
            i3 = get<2>(it);
            //cout << "(i1, i2, i3): " << i1 << ", " << i2 << ", " << i3<< endl;
            Constraint Clique("ZClique("+to_string(i1)+","+to_string(i2)+ ","+to_string(i3)+")");
            Clique = zij(i1,i2) +zij(i1,i3) + zij(i2,i3);
            _model.add_constraint(Clique >=1);
        }
    }
}

void Minkmodel::add_general_clique(){
//  t=1, q=2,...,k-1
    for (int i = 0; i < _graph->_bags.size(); i++){
        auto bag = _graph->_bags.at(i);
        if (bag.size() > _K+1) {
            Constraint GClique("General_Clique"+to_string(i));
            int t = floor(bag.size()/_K);
            int q = (bag.size()-t*_K);
            for (int j= 0; j< bag.size()-1; j++){
                int idj=bag[j]->ID;
                for (int k=j+1; k< bag.size(); k++){
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

void Minkmodel::add_wheel(){
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

void Minkmodel::add_bicycle(){}

void Minkmodel::add_clique_lifted_tree(){
    int i1,i2,i3,i4;
    if (_K>2)
        for (auto it: _ids4){
            i1 = get<0>(it);
            i2 = get<1>(it);
            i3 = get<2>(it);
            i4 = get<3> (it);
            //cout << "(i1, i2, i3, i4): " << i1 << ", " << i2 << ", " << i3 << ", " << i4<< endl;
            Constraint Clique("XClique("+to_string(i1)+","+to_string(i2)+ ","+to_string(i3)+ ", "+to_string(i4)+")");
            Clique = Xij(i1,i2) +Xij(i1,i3) + Xij(i1,i4) + Xij(i2,i3) + Xij(i2,i4) +Xij(i3,i4);
            _model.add_constraint(Clique >= -0.5*_K);
        }
    else{
        for (auto it: _ids){
            i1 = get<0>(it);
            i2 = get<1>(it);
            i3 = get<2>(it);
            //cout << "(i1, i2, i3): " << i1 << ", " << i2 << ", " << i3<< endl;
            Constraint Clique("XClique("+to_string(i1)+","+to_string(i2)+ ","+to_string(i3)+")");
            Clique = Xij(i1,i2) +Xij(i1,i3) + Xij(i2,i3);
            _model.add_constraint(Clique >= -0.5*_K);
        }
    }
}

int Minkmodel::solve(int output, bool relax){
    solver s(_model,_solver);
    //double wall0 = get_wall_time();
    //double cpu0  = get_cpu_time();
    cout << "Running the relaxation model \n";
    s.run(output,relax);
    //double wall1 = get_wall_time();
    //double cpu1  = get_cpu_time();
    //relax.print_solution();
    return 1;
}

//bool Minkmodel::check_eigenvalues(){
//   arma::SpMat<double> A(_graph->nodes.size(), _graph->nodes.size());
//    //for (i= 0; i < _data->nbV[0]*nbServers; i++)
//    //  for (j= 0; j < _data->nbV[0]*nbServers; j++){
//    // matrix A contain values of X - x*x^T.
//    //A(i,j) = (double)val_M[0][i][j];
//    //}
//    arma::vec eigval;
//    arma::eigs_sym(eigval,_eigvec,A, 1,"sa",0.000001);
//    if (eigval(0) < -0.0001)
//    {
//        cout << "there exists negative eigenvalue: " << eigval << endl;
//        return true;
//    }
//    else{
//        cout << "The solution x satisfy SDP constraint  X>=0" << endl;
//        return false;
//    }
//}
