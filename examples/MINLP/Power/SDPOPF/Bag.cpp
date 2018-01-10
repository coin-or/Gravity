//
// Created by kbestuzheva on 12/18/17.
//

#include <gravity/param.h>
#include "Bag.h"

Bag::Bag(int id, PowerNet* grid, vector<Node*> nodes):_id(id),_grid(grid),_nodes(nodes) {
    Arc* aij = NULL;
    _wstarp.set_name("wstar");
    _wmin.set_name("wmin");
    _wmax.set_name("wmax");
    string namewr, namewi, namew, namepair;
    bool reversed;
    for(int i = 0; i < _nodes.size()-1; i++){
        for(int j = i+1; j < _nodes.size(); j++){
            aij = _grid->get_directed_arc(_nodes[i]->_name,_nodes[j]->_name);
            if(aij==NULL) {
                aij = _grid->get_directed_arc(_nodes[j]->_name,_nodes[i]->_name);
                cout << "\nBag id = " << _id << ", reversed arc " << _nodes[i]->_name << "," << _nodes[j]->_name;
                reversed = true;
                namepair = _nodes[j]->_name+","+_nodes[i]->_name;
            }else{
                reversed = false;
                namepair = _nodes[i]->_name+","+_nodes[j]->_name;
            }
            if(aij->_imaginary) {
                _all_lines = false;
                cout << "\nMissing lines in bag, skipping";
                return;
            }

            namewr = "wr(" + _nodes[i]->_name + "," + _nodes[j]->_name + ")";
            namewi = "wi(" + _nodes[i]->_name + "," + _nodes[j]->_name + ")";

            _indices.push_back(new index_(namewr));
            _indices.push_back(new index_(namewi));

            _wmin.set_val(namewr,_grid->wr_min(namepair).eval());
            _wmin.set_val(namewi,_grid->wi_min(namepair).eval());
            _wmax.set_val(namewr,_grid->wr_max(namepair).eval());
            _wmax.set_val(namewi,_grid->wi_max(namepair).eval());

            _wstarp.set_val(namewr,((Line*)aij)->wr);
            if(!reversed) _wstarp.set_val(namewi,((Line*)aij)->wi);
            else _wstarp.set_val(namewi,-((Line*)aij)->wi);
        }
    }
    for(int i = 0; i < _nodes.size(); i++){
        namew = "w(" + _nodes[i]->_name + ")";
        _indices.push_back(new index_(namew));
        _wmin.set_val(namew,_grid->w_min(_nodes[i]->_name).eval());
        _wmax.set_val(namew,_grid->w_max(_nodes[i]->_name).eval());
        _wstarp.set_val(namew,((Bus*)_nodes[i])->w);
    }

//    cout << "\nBus pairs: ";
//    for(auto& bp: _bus_pairs._keys) cout << "\n" << bp->_name;
//    cout << "\n---------\n";

//    _wstarp.print(true); cout << "\n";
//    _wmin.print(true); cout << "\n";
//    _wmax.print(true);
};

Bag::~Bag(){
    for(gravity::index_* idx: _indices) {
        delete idx;
    }
    _indices.clear();
}

bool Bag::add_lines(){
//    if(_nodes.size() != 3 || _all_lines) return;
    Line *a12, *a13, *a23;
    Bus *n1, *n2, *n3;
    double tol = 0.00001;
    n1 = (Bus*)_grid->get_node(_nodes[0]->_name);
    n2 = (Bus*)_grid->get_node(_nodes[1]->_name);
    n3 = (Bus*)_grid->get_node(_nodes[2]->_name);
    a12 = (Line*)_grid->get_arc(n1,n2);
    a13 = (Line*)_grid->get_arc(n1,n3);
    a23 = (Line*)_grid->get_arc(n2,n3);
//    if((a12->_free && a13->_free) || (a23->_free && a13->_free) || (a12->_free && a23->_free)) //more than 2 unassigned lines
//        return;

    double wr12, wi12, wr13, wi13, wr23, wi23;
    double w1 = n1->w; double w2 = n2->w; double w3 = n3->w;

    wr12 = a12->wr; wr13 = a13->wr; wr23 = a23->wr;
    if(_grid->get_directed_arc(n1->_name,n2->_name)!=nullptr) wi12 = a12->wi;
    else wi12 = -a12->wi;
    if(_grid->get_directed_arc(n1->_name,n3->_name)!=nullptr) wi13 = a13->wi;
    else wi13 = -a13->wi;
    if(_grid->get_directed_arc(n2->_name,n3->_name)!=nullptr) wi23 = a23->wi;
    else wi23 = -a23->wi;

    string s12, s13, s23;
    s12 = a12->_src->_name+","+a12->_dest->_name;
    s13 = a13->_src->_name+","+a13->_dest->_name;
    s23 = a23->_src->_name+","+a23->_dest->_name;

    if(a13->_free) {
        cout << "Calculating values for " << a13->_name;
        a13->_free = false;

        a13->wr = (wr12 * wr23 - wi12 * wi23) / w2;
        if(_grid->get_directed_arc(n1->_name,n3->_name)!=nullptr) a13->wi = (wi12 * wr23 + wr12 * wi23) / w2;
        else a13->wi = -(wi12 * wr23 + wr12 * wi23) / w2;


        double SDP = wr12*(wr23*wr13 + wi23*wi13) + wi12*(-wi23*wr13 + wr23*wi13);
        SDP *= 2;
        SDP -= (wr12*wr12 + wi12*wi12)*w3 + (wr13*wr13 + wi13*wi13)*w2 + (wr23*wr23 + wi23*wi23)*w1;
        SDP += w1*w2*w3;
//            double R1 = wr13*wr13 + wi13*wi13 - ((wr12*wr12+wi12*wi12)*w3 + (wr32*wr32+wi32*wi32)*w1 - w1*w2*w3)/w2;
//            cout << "\nR1 = " << R1 << ", SOC = " << wr12*wr12 + wi12*wi12 - w1*w2 << ", SOC2 = " << wr32*wr32 + wi32*wi32 - w2*w3;
        cout << "\nNo a13, SDP = " << SDP;

//        if(_net->sdp_alg==1) return false;

        if (!(wr13 >= _grid->wr_min(s13).eval()-tol && wr13 <= _grid->wr_max(s13).eval()+tol
              && wi13 >= _grid->wi_min(s13).eval()-tol && wi13 <= _grid->wi_max(s13).eval()+tol)){
//            cout << "\nBounds are violated";
            return false;
        }
        if(wr13*wr13+wi13*wi13 > w1*w3+tol) {
//            cout << "\nSOCP is violated";
            return false;
        }
        return true;
    }

    if(a23->_free) {
        a23->_free = false;

//        cout << "Calculating values for " << a23->_name;
        a23->wr = (wr12 * wr13 + wi12 * wi13) / w1;
        if(_grid->get_directed_arc(n2->_name,n3->_name)!=nullptr) a23->wi = (wr12 * wi13 - wi12 * wr13) / w1;
        else a23->wi = -(wr12 * wi13 - wi12 * wr13) / w1;

        double SDP = wr12*(wr23*wr13 + wi23*wi13) + wi12*(-wi23*wr13 + wr23*wi13);//todo: check these, did I change the sign of wi23 in the prev case?
        SDP *= 2;
        SDP -= (wr12*wr12 + wi12*wi12)*w3 + (wr13*wr13 + wi13*wi13)*w2 + (wr23*wr23 + wi23*wi23)*w1;
        SDP += w1*w2*w3;
//            double R1 = (wr12*wr12*wr13*wr13 + wi12*wi12*wi13*wi13 + wi12*wi12*wr13*wr13 + wr12*wr12*wi13*wi13)/(w1*w1);
//            R1 -= ((wr13*wr13+wi13*wi13)*w2 + (wr12*wr12+wi12*wi12)*w3 - w1*w2*w3)/w1;
//            cout << "\nR1 = " << R1 << ", SOC = " << wr12*wr12 + wi12*wi12 - w1*w2 << ", SOC2 = " << wr13*wr13 + wi13*wi13 - w1*w3;
        cout << "\nNo a32, SDP = " << SDP;

//        if(_net->sdp_alg==1) return false;

        if (!(wr23 >= _grid->wr_min(s23).eval()-tol && wr23 <= _grid->wr_max(s23).eval()+tol && wi23 >= _grid->wi_min(s23).eval()-tol
              && wi23 <= _grid->wi_max(s23).eval()+tol && wr23*wr23+wi23*wi23 > w2*w3+tol)){
//            cout << "\nBounds or SOCP is violated";
            return false;
        }
        return true;
    }

    if(a12->_free) {
        a12->_free = false;

//        cout << "Calculating values for " << a12->_name;
        a12->wr = (wr23 * wr13 + wi23 * wi13) / w3;
        if(_grid->get_directed_arc(n1->_name,n2->_name)!=nullptr) a12->wi = (-wi23 * wr13 + wr23 * wi13) / w3;
        else a12->wi = -(-wi23 * wr13 + wr23 * wi13) / w3;

        double SDP = wr12*(wr23*wr13 + wi23*wi13) + wi12*(-wi23*wr13 + wr23*wi13);
        SDP *= 2;
        SDP -= (wr12*wr12 + wi12*wi12)*w3 + (wr13*wr13 + wi13*wi13)*w2 + (wr23*wr23 + wi23*wi23)*w1;
        SDP += w1*w2*w3;
//            double R1 = wr12*wr12 + wi12*wi12 - ((wr13*wr13+wi13*wi13)*w2 + (wr32*wr32+wi32*wi32)*w1 - w1*w2*w3)/w3;
//            cout << "\nR1 = " << R1;
        cout << "\nNo a12, SDP = " << SDP;

//        if(_net->sdp_alg==1) return false;

        if (!(wr12 >= _grid->wr_min(s12).eval()-tol && wr12 <= _grid->wr_max(s12).eval()+tol && wi12 >= _grid->wi_min(s12).eval()-tol
              && wi12 <= _grid->wi_max(s12).eval()+tol && wr12*wr12+wi12*wi12 > w1*w2+tol)){
//            cout << "\nBounds or SOCP is violated";
            return false;
        }
        return true;
    }
    return false;
}

bool Bag::is_PSD(){
    if(_nodes.size() != 3) return true;
    int free_lines = 0;
    for(int i = 0; i < _nodes.size()-1; i++){
        for(int j = i+1; j < _nodes.size(); j++){
            Arc *aij = _grid->get_arc(_nodes[i]->_name, _nodes[j]->_name);
            if(aij->_free) free_lines++;
        }
    }
    if(free_lines > 1) return true;
    if(free_lines == 1) return add_lines();

    //the bag has all lines

    double tol = 0.0005;
//    int n = igraph_vector_size(cl);
    int n = _nodes.size();
    Node* node;
    arma::cx_mat A(n,n);

    for(int i = 0; i < n; i++){
        for(int j = 0; j <= i; j++) {
            if(i==j) {
                A(i,j) = arma::cx_double(((Bus*)_grid->get_node(_nodes[i]->_name))->w,0);
            }else{
                double AijI;
                if(_grid->has_directed_arc(_nodes[j],_nodes[i])) {
                    AijI = ((Line*)_grid->get_arc(_nodes[j], _nodes[i]))->wi;
                }
                else {
                    AijI = -((Line*)_grid->get_arc(_nodes[j], _nodes[i]))->wi;
                }
                A(i,j) = arma::cx_double(((Line*)_grid->get_arc(_nodes[j],_nodes[i]))->wr,AijI);
                A(j,i) = arma::cx_double(((Line*)_grid->get_arc(_nodes[j],_nodes[i]))->wr,-AijI);
            }
        }
    }
//    A.print();
    arma::cx_mat R;
    arma::vec v = arma::eig_sym(A);
    cout << "\n";
    double min_eig = 0, max_eig = -1;
    for(auto eig: v) {
        if(eig < min_eig) min_eig = eig;
        if(eig > max_eig) max_eig = eig;
    }
    if(min_eig/max_eig > -tol) {cout << "\nPSD"; return true;}
    else {cout << "\nNot PSD"; return false;}
}

param<double> Bag::nfp(){
    param<double> what;

//    cout << "\nIndices:";
//    for(auto& i: _indices) cout << "\n" << i->_name;
//    cout << "\n---------\n";

    //    var<double> R_Wij("R_Wij", _grid->wr_min.in(_bus_pairs._keys), _grid->wr_max.in(_bus_pairs._keys));
//    NPP.add_var(R_Wij ^ (n*(n-1)/2));
//
//    var<double> I_Wij("I_Wij", _grid->wi_min.in(_bus_pairs._keys), _grid->wi_max.in(_bus_pairs._keys));
//    NPP.add_var(I_Wij ^ (n*(n-1)/2));
//
//    var<double> Wii("Wii", _grid->w_min.in(_nodes), _grid->w_max.in(_nodes));
//    NPP.add_var(Wii ^ (n));

    Model NPP("NPP model");
    int n = _nodes.size();

    sdpvar<double> W("W");
    NPP.add_var(W ^ (2*n));

    var<double> z("z");
    z.in_q_cone();
    NPP.add_var(z ^ (n*n+1));

    var<double> w("w", _wmin.in(_indices), _wmax.in(_indices));
    NPP.add_var(w ^ (n*n));

    func_ obj;
    obj = z("obj");
    NPP.set_objective(min(obj));

//    min t  s.t.
//    ||z|| <= t    (z(obj) = t)
//    w*-w = z
//    X is PSD,
//    L <= X <= U
//    matr structure

    /** Constraints **/

    /* w*-w = z */
    Constraint svec("svec");
    svec = _wstarp.in(_indices) - w.in(_indices) - z.in(_indices);
    cout << "\nadding svec";
    NPP.add_constraint(svec==0);

    string namew, namewr, namewi;

    /* matrix structure */
    for(int i = 0; i < n; i++){
        for(int j = i; j < n; j++){
            if(i==j){
                namew = "w(" + _nodes[i]->_name + ")";

                Constraint mstruct1("mstruct"+to_string(i+n)+to_string(i+n));
                mstruct1 = w(namew) - W(i+n,i+n);
                NPP.add_constraint(mstruct1==0);

                Constraint mstruct2("mstruct"+to_string(i)+to_string(i));
                mstruct2 = w(namew) - W(i,i);
                NPP.add_constraint(mstruct2==0);

                /* zeros */
                Constraint zero("zero"+to_string(i)+to_string(i+n));
                zero = W(i,i+n);
                NPP.add_constraint(zero==0);
                Constraint zero2("zero"+to_string(i+n)+to_string(i));
                zero2 = W(i+n,i);
                NPP.add_constraint(zero2==0);
            }
            else {
                namewr = "wr(" + _nodes[i]->_name + "," + _nodes[j]->_name + ")";

                Constraint mstruct("mstructr"+to_string(i)+to_string(j));
                mstruct = w(namewr) - W(i,j);
                NPP.add_constraint(mstruct==0);
                Constraint mstruct2("mstructr"+to_string(i+n)+to_string(j+n));
                mstruct2 = w(namewr) - W(i+n,j+n);
                NPP.add_constraint(mstruct2==0);

                namewi = "wi(" + _nodes[i]->_name + "," + _nodes[j]->_name + ")";

                Constraint mstructi("mstructi"+to_string(i)+to_string(j+n));
                mstructi = w(namewi) - W(i,j+n);
                NPP.add_constraint(mstructi==0);
                Constraint mstructi2("mstructi"+to_string(j)+to_string(i+n));
                mstructi2 = w(namewi) + W(j,i+n);
                NPP.add_constraint(mstructi2==0);
            }
        }
    }

    solver s(NPP,mosek_);
    s.run(1,0);

//    z.print(); cout << "\n";
//    W.print(true); cout << "\n";
//    w.print(); cout << "\n";

    for(int i = 0; i < n; i++){
        for(int j = i; j < n; j++){
            if(i==j){
                namew = "w(" + _nodes[i]->_name + ")";
                what.set_val(namew,w(namew).eval());
            }
            else {
                namewr = "wr(" + _nodes[i]->_name + "," + _nodes[j]->_name + ")";
                namewi = "wi(" + _nodes[i]->_name + "," + _nodes[j]->_name + ")";
                what.set_val(namewr, w(namewr).eval());
                if(_grid->get_directed_arc(_nodes[i]->_name,_nodes[j]->_name)!=nullptr)
                    what.set_val(namewi, w(namewi).eval());
                else
                    what.set_val(namewi, -w(namewi).eval());
            }
        }
    }
    what.set_name("w_hat");
    what.print(true);
    return what;
}
