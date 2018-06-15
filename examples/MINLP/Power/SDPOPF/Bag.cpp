//
// Created by kbestuzheva on 12/18/17.
//

#include <Eigen/Dense>
//#include <Eigen/Eigenvalues>
#include <gravity/param.h>
#include "Bag.h"

#define FLAT

Bag::Bag(int id, const PowerNet& grid, vector<Node*> nodes):_id(id),_grid((PowerNet*)&grid),_nodes(nodes) {
    Arc* aij = NULL;
    _wstarp.set_name("wstar");
    _W_star.set_name("W_star");
    _wmin.set_name("wmin");
    _wmax.set_name("wmax");
    string namewr, namewi, namew, namepair;
//    bool reversed;
    for(int i = 0; i < _nodes.size()-1; i++){
        for(int j = i+1; j < _nodes.size(); j++){
            aij = _grid->get_directed_arc(_nodes[i]->_name,_nodes[j]->_name);
            if(aij==NULL) {
                aij = _grid->get_directed_arc(_nodes[j]->_name,_nodes[i]->_name);
//                cout << "\nBag id = " << _id << ", reversed arc " << _nodes[i]->_name << "," << _nodes[j]->_name;
//                reversed = true;
                namepair = _nodes[j]->_name+","+_nodes[i]->_name;
            }else{
//                reversed = false;
                namepair = _nodes[i]->_name+","+_nodes[j]->_name;
            }
//            if(aij->_imaginary) {
//                _all_lines = false; //todo: do I need this?
////                cout << "\nMissing lines in bag, skipping";
////                return;
//            }

            namewr = "wr(" + _nodes[i]->_name + "," + _nodes[j]->_name + ")";
            namewi = "wi(" + _nodes[i]->_name + "," + _nodes[j]->_name + ")";

//            cout << "\n" << namewr << " = " <<  ((Line*)aij)->wr;
//            cout << "\n" << namewi << " = " <<  ((Line*)aij)->wi;

            _indices.push_back(index_(namewr));
            _indices.push_back(index_(namewi));

            _wmin.set_val(namewr,_grid->wr_min(namepair).eval());
            _wmin.set_val(namewi,_grid->wi_min(namepair).eval());
            _wmax.set_val(namewr,_grid->wr_max(namepair).eval());
            _wmax.set_val(namewi,_grid->wi_max(namepair).eval());

            _Wmin.set_val(to_string(i)+","+to_string(j),_grid->wr_min(namepair).eval());
            _Wmin.set_val(to_string(i+_nodes.size())+","+to_string(j+_nodes.size()),_grid->wr_min(namepair).eval());
            _Wmin.set_val(to_string(i)+","+to_string(j+_nodes.size()),_grid->wi_min(namepair).eval());
            _Wmin.set_val(to_string(j)+","+to_string(i+_nodes.size()),_grid->wi_min(namepair).eval());

            _Wmax.set_val(to_string(i)+","+to_string(j),_grid->wr_max(namepair).eval());
            _Wmax.set_val(to_string(i+_nodes.size())+","+to_string(j+_nodes.size()),_grid->wr_max(namepair).eval());
            _Wmax.set_val(to_string(i)+","+to_string(j+_nodes.size()),_grid->wi_max(namepair).eval());
            _Wmax.set_val(to_string(j)+","+to_string(i+_nodes.size()),_grid->wi_max(namepair).eval());

//            _wstarp.set_val(namewr,((Line*)aij)->wr);
//            if(!reversed) _wstarp.set_val(namewi,((Line*)aij)->wi);
//            else _wstarp.set_val(namewi,-((Line*)aij)->wi);

//            _W_star.set_val(to_string(i)+","+to_string(j),((Line*)aij)->wr);
//            _W_star.set_val(to_string(i+_nodes.size())+","+to_string(j+_nodes.size()),((Line*)aij)->wr);
//            _W_star.set_val(to_string(i)+","+to_string(j+_nodes.size()),((Line*)aij)->wi);
//            _W_star.set_val(to_string(j)+","+to_string(i+_nodes.size()),-((Line*)aij)->wi);
        }
    }

    for(int i = 0; i < _nodes.size(); i++){
        namew = "w(" + _nodes[i]->_name + ")";
        _indices.push_back(index_(namew));
        _wmin.set_val(namew,_grid->w_min(_nodes[i]->_name).eval());
        _wmax.set_val(namew,_grid->w_max(_nodes[i]->_name).eval());

        _Wmin.set_val(to_string(i)+","+to_string(i),_grid->w_min(_nodes[i]->_name).eval());
        _Wmin.set_val(to_string(i+_nodes.size())+","+to_string(i+_nodes.size()),_grid->w_min(_nodes[i]->_name).eval());
        _Wmin.set_val(to_string(i)+","+to_string(i+_nodes.size()),0);

        _Wmax.set_val(to_string(i)+","+to_string(i),_grid->w_max(_nodes[i]->_name).eval());
        _Wmax.set_val(to_string(i+_nodes.size())+","+to_string(i+_nodes.size()),_grid->w_max(_nodes[i]->_name).eval());
        _Wmax.set_val(to_string(i)+","+to_string(i+_nodes.size()),0);

//        _wstarp.set_val(namew,((Bus*)_nodes[i])->w);
//        cout << "\n" << namew << " = " <<  ((Bus*)_nodes[i])->w;

//        _W_star.set_val(to_string(i)+","+to_string(i),((Bus*)_nodes[i])->w);
//        _W_star.set_val(to_string(i+_nodes.size())+","+to_string(i+_nodes.size()),((Bus*)_nodes[i])->w);
//        _W_star.set_val(to_string(i)+","+to_string(i+_nodes.size()),0.0);
    }

//    cout << "\nBus pairs: ";
//    for(auto& bp: _bus_pairs._keys) cout << "\n" << bp->_name;
//    cout << "\n---------\n";

//    _wstarp.print(true); cout << "\n";
//    _wmin.print(true); cout << "\n";
//    _wmax.print(true);
//    _W_star.print(true);
};

Bag::~Bag(){
    _indices.clear();
}

param<double> Bag::fill_wstar(){
    param<double> W_star;
    W_star.set_name("W_star");
    Line* aij;
    bool reversed;
    string namewr, namewi, namew;

    for(int i = 0; i < _nodes.size()-1; i++){
        for(int j = i+1; j < _nodes.size(); j++){
            aij = (Line*)_grid->get_directed_arc(_nodes[i]->_name,_nodes[j]->_name);
            if(aij==NULL) {
                aij = (Line*)_grid->get_directed_arc(_nodes[j]->_name,_nodes[i]->_name);
//                cout << "\nBag id = " << _id << ", reversed arc " << _nodes[i]->_name << "," << _nodes[j]->_name;
                reversed = true;
            }else reversed = false;

            namewr = "wr(" + _nodes[i]->_name + "," + _nodes[j]->_name + ")";
            namewi = "wi(" + _nodes[i]->_name + "," + _nodes[j]->_name + ")";

            _wstarp.set_val(namewr,aij->wr);
            if(!reversed) _wstarp.set_val(namewi,aij->wi);
            else _wstarp.set_val(namewi,-aij->wi);

            _W_star.set_val(to_string(i)+","+to_string(j),aij->wr);
            _W_star.set_val(to_string(i+_nodes.size())+","+to_string(j+_nodes.size()),aij->wr);
            if(!reversed) {
                _W_star.set_val(to_string(i) + "," + to_string(j + _nodes.size()), aij->wi);
                _W_star.set_val(to_string(j) + "," + to_string(i + _nodes.size()), -aij->wi);
            }else{
                _W_star.set_val(to_string(i) + "," + to_string(j + _nodes.size()), -aij->wi);
                _W_star.set_val(to_string(j) + "," + to_string(i + _nodes.size()), aij->wi);
            }
        }
    }
    for(int i = 0; i < _nodes.size(); i++){
        namew = "w(" + _nodes[i]->_name + ")";
        _wstarp.set_val(namew,((Bus*)_nodes[i])->w);

        _W_star.set_val(to_string(i)+","+to_string(i),((Bus*)_nodes[i])->w);
        _W_star.set_val(to_string(i+_nodes.size())+","+to_string(i+_nodes.size()),((Bus*)_nodes[i])->w);
        _W_star.set_val(to_string(i)+","+to_string(i+_nodes.size()),0.0);
    }
//    cout << "\nwstar: "; _wstarp.print(true);
    return W_star;
}

bool Bag::add_lines(){
    // this is only called for 3d bags with 1 missing line
//    if(_nodes.size() != 3 || _all_lines) return;

    int free_lines = 0;
    for (int i = 0; i < _nodes.size() - 1; i++) {
        for (int j = i + 1; j < _nodes.size(); j++) {
            Arc *aij = _grid->get_arc(_nodes[i]->_name, _nodes[j]->_name);
            if (aij->_free) free_lines++;
        }
    }
    if (free_lines !=1) return true;

    Line *a12, *a13, *a23;
    Bus *n1, *n2, *n3;
    double tol = 0.00001;
    bool psd = true;
    n1 = (Bus*)_grid->get_node(_nodes[0]->_name);
    n2 = (Bus*)_grid->get_node(_nodes[1]->_name);
    n3 = (Bus*)_grid->get_node(_nodes[2]->_name);
    a12 = (Line*)_grid->get_arc(n1,n2);
    a13 = (Line*)_grid->get_arc(n1,n3);
    a23 = (Line*)_grid->get_arc(n2,n3);

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
        DebugOff("\nCalculating values for a13 = " << a13->_name << endl);
        a13->_free = false;

        a13->wr = (wr12 * wr23 - wi12 * wi23) / w2;
        if(_grid->get_directed_arc(n1->_name,n3->_name)!=nullptr) a13->wi = (wi12 * wr23 + wr12 * wi23) / w2;
        else a13->wi = -(wi12 * wr23 + wr12 * wi23) / w2;

//        wi13 = (wi12 * wr32 - wr12 * wi32) / w2;

//        double SDP = wr12*(wr23*a13->wr + wi23*a13->wi) + wi12*(-wi23*a13->wr + wr23*a13->wi);
//        SDP *= 2;
//        SDP -= (wr12*wr12 + wi12*wi12)*w3 + (a13->wr*a13->wr + a13->wi*a13->wi)*w2 + (wr23*wr23 + wi23*wi23)*w1;
//        SDP += w1*w2*w3;
//        double R1 = wr13*wr13 + wi13*wi13 - ((wr12*wr12+wi12*wi12)*w3 + (wr23*wr23+wi23*wi23)*w1 - w1*w2*w3)/w2;
//            cout << "\nR1 = " << R1 << ", SOC = " << wr12*wr12 + wi12*wi12 - w1*w2 << ", SOC2 = " << wr23*wr23 + wi23*wi23 - w2*w3;
//        cout << "\nNo a13, SDP = " << SDP;

//        if(_net->sdp_alg==1) return false;

        if (!(wr13 >= _grid->wr_min(s13).eval()-tol && wr13 <= _grid->wr_max(s13).eval()+tol
              && wi13 >= _grid->wi_min(s13).eval()-tol && wi13 <= _grid->wi_max(s13).eval()+tol)){
            DebugOff("\nBounds are violated");
            psd = false;
        }
        if(wr13*wr13+wi13*wi13 > w1*w3+tol) {
            DebugOff("\nSOCP is violated");
            psd = false;
        }
    }

    if(a23->_free) {
        a23->_free = false;

        DebugOff("\nCalculating values for a23 = " << a23->_name << endl);
        a23->wr = (wr12 * wr13 + wi12 * wi13) / w1;
        if(_grid->get_directed_arc(n2->_name,n3->_name)!=nullptr) a23->wi = (wr12 * wi13 - wi12 * wr13) / w1;
        else a23->wi = -(wr12 * wi13 - wi12 * wr13) / w1;

//        double SDP = wr12*(wr23*wr13 + wi23*wi13) + wi12*(-wi23*wr13 + wr23*wi13);
//        SDP *= 2;
//        SDP -= (wr12*wr12 + wi12*wi12)*w3 + (wr13*wr13 + wi13*wi13)*w2 + (wr23*wr23 + wi23*wi23)*w1;
//        SDP += w1*w2*w3;
//            double R1 = (wr12*wr12*wr13*wr13 + wi12*wi12*wi13*wi13 + wi12*wi12*wr13*wr13 + wr12*wr12*wi13*wi13)/(w1*w1);
//            R1 -= ((wr13*wr13+wi13*wi13)*w2 + (wr12*wr12+wi12*wi12)*w3 - w1*w2*w3)/w1;
//            cout << "\nR1 = " << R1 << ", SOC = " << wr12*wr12 + wi12*wi12 - w1*w2 << ", SOC2 = " << wr13*wr13 + wi13*wi13 - w1*w3;
//        cout << "\nNo a32, SDP = " << SDP;

//        if(_net->sdp_alg==1) return false;

        if (!(wr23 >= _grid->wr_min(s23).eval()-tol && wr23 <= _grid->wr_max(s23).eval()+tol && wi23 >= _grid->wi_min(s23).eval()-tol
              && wi23 <= _grid->wi_max(s23).eval()+tol && wr23*wr23+wi23*wi23 <= w2*w3+tol)){
            DebugOff("\nBounds or SOCP is violated");
            psd = false;
        }
    }

    if(a12->_free) {
        a12->_free = false;

        DebugOff("\nCalculating values for a12 = " << a12->_name << endl);
        a12->wr = (wr23 * wr13 + wi23 * wi13) / w3;
        if(_grid->get_directed_arc(n1->_name,n2->_name)!=nullptr) a12->wi = (-wi23 * wr13 + wr23 * wi13) / w3;
        else a12->wi = -(-wi23 * wr13 + wr23 * wi13) / w3;

//        double SDP = wr12*(wr23*wr13 + wi23*wi13) + wi12*(-wi23*wr13 + wr23*wi13);
//        SDP *= 2;
//        SDP -= (wr12*wr12 + wi12*wi12)*w3 + (wr13*wr13 + wi13*wi13)*w2 + (wr23*wr23 + wi23*wi23)*w1;
//        SDP += w1*w2*w3;
//            double R1 = wr12*wr12 + wi12*wi12 - ((wr13*wr13+wi13*wi13)*w2 + (wr32*wr32+wi32*wi32)*w1 - w1*w2*w3)/w3;
//            cout << "\nR1 = " << R1;
//        cout << "\nNo a12, SDP = " << SDP;

//        if(_net->sdp_alg==1) return false;

        if (!(wr12 >= _grid->wr_min(s12).eval()-tol && wr12 <= _grid->wr_max(s12).eval()+tol && wi12 >= _grid->wi_min(s12).eval()-tol
              && wi12 <= _grid->wi_max(s12).eval()+tol && wr12*wr12+wi12*wi12 <= w1*w2+tol)){
            DebugOff("\nBounds or SOCP is violated");
            psd = false;
        }
    }
    return psd;
}

bool Bag::is_PSD(){
    if(_nodes.size() == 2) return true;

    int free_lines = 0;
    for (int i = 0; i < _nodes.size() - 1; i++) {
        for (int j = i + 1; j < _nodes.size(); j++) {
            Arc *aij = _grid->get_arc(_nodes[i]->_name, _nodes[j]->_name);
            if (aij->_free) free_lines++;
        }
    }
    if (free_lines > 1) return true;
    if (free_lines == 1 && _nodes.size() == 3) return add_lines();

    if(_nodes.size()>3) {
        if(free_lines > 0) DebugOn("\nFree lines in a larger bag!");
        else DebugOff("\nNo free lines in a larger bag :)");
//        return true;
    }

    //the bag has all lines or has > 3 nodes

    double tol = 0.0001;
//    int n = igraph_vector_size(cl);
    int n = _nodes.size();
    Node* node;
    arma::cx_mat A(n,n);

    fill_wstar();

    Line *a12, *a13, *a23;
    Bus *n1, *n2, *n3;
    n1 = (Bus*)_grid->get_node(_nodes[0]->_name);
    n2 = (Bus*)_grid->get_node(_nodes[1]->_name);
    n3 = (Bus*)_grid->get_node(_nodes[2]->_name);
    a12 = (Line*)_grid->get_arc(n1,n2);
    a13 = (Line*)_grid->get_arc(n1,n3);
    a23 = (Line*)_grid->get_arc(n2,n3);

    double wr12, wi12, wr13, wi13, wr23, wi23;
    double w1 = n1->w; double w2 = n2->w; double w3 = n3->w;

    wr12 = a12->wr; wr13 = a13->wr; wr23 = a23->wr;
    if(_grid->get_directed_arc(n1->_name,n2->_name)!=nullptr) wi12 = a12->wi;
    else wi12 = -a12->wi;
    if(_grid->get_directed_arc(n1->_name,n3->_name)!=nullptr) wi13 = a13->wi;
    else wi13 = -a13->wi;
    if(_grid->get_directed_arc(n2->_name,n3->_name)!=nullptr) wi23 = a23->wi;
    else wi23 = -a23->wi;

    double SDP = wr12*(wr23*wr13 + wi23*wi13) + wi12*(-wi23*wr13 + wr23*wi13);
    SDP *= 2;
    SDP -= (wr12*wr12 + wi12*wi12)*w3 + (wr13*wr13 + wi13*wi13)*w2 + (wr23*wr23 + wi23*wi23)*w1;
    SDP += w1*w2*w3;
    DebugOff("\nSDP = " << SDP);

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
    arma::cx_mat P;
    arma::vec v;// = arma::eig_sym(A);
    arma::eig_sym(v,P,A);
    DebugOff("\n");
    double min_eig = 0, max_eig = -1;
    for(auto eig: v) {
        if(eig < min_eig) min_eig = eig;
        if(eig > max_eig) max_eig = eig;
    }
    if(min_eig/max_eig > -tol) {
        DebugOff("\nBag is PSD");
        return true;
    }
    else {
        double pos_tol = -0.00001;
        DebugOff("\nBag is not PSD\n");
        for(int i = 0; i < n; i++) {
            if(v[i] < 0) {
                v[i] = 0;
//                P.col(i).zeros();
            }
        }
        arma::mat D(n,n);
        D.zeros();
        D.diag() = v;
        arma::cx_mat W_hat = P*D*P.t();

        return false;
    }
}

void Bag::update_PSD(){
    _is_psd = is_PSD();
}

/* the indices of all elements in the upper triangular matrix of size n */
vector<index_> triang_indices(int n) {
    vector<index_> res;
    for(int i = 0; i < n; i++) {
        for(int j = i; j < n; j++) {
            res.push_back(index_(to_string(i)+","+to_string(j)));
        }
    }
    return res;
}

/* the following functions return the indices corresponding to some parts
 * of a 2n x 2n matrix that is divided into 4 submatrices of equal size: */

/* upper triangular part (including the diagonal) of the upper left submatrix */
vector<index_> ul_u(int n) {
    vector<index_> res;
    for(int i = 0; i < n; i++) {
        for(int j = i; j < n; j++) {
            res.push_back(index_(to_string(i)+","+to_string(j)));
        }
    }
    return res;
}

/* upper triangular part (including the diagonal) of the lower right submatrix */
vector<index_> lr_u(int n) {
    vector<index_> res;
    for(int i = n; i < 2*n; i++) {
        for(int j = i; j < 2*n; j++) {
            res.push_back(index_(to_string(i)+","+to_string(j)));
        }
    }
    return res;
}

/* upper triangular part (excluding the diagonal) of the upper right submatrix */
vector<index_> ur_u(int n) {
    vector<index_> res;
    for(int i = 0; i < n; i++) {
        for(int j = n+i+1; j < 2*n; j++) {
            res.push_back(index_(to_string(i)+","+to_string(j)));
        }
    }
    return res;
}

/* lower triangular part (excluding the diagonal) of the upper right submatrix */
vector<index_> ur_l(int n) {
    vector<index_> res;
    for(int i = 0; i < n; i++) {
        for(int j = n; j < n+i; j++) {
            res.push_back(index_(to_string(i)+","+to_string(j)));
        }
    }
    return res;
}

/*---------------------------------------*/
#ifdef FLAT
param<double> Bag::nfp(){
    param<double> what;
    fill_wstar();

//    cout << "\nIndices:";
//    for(auto& i: _indices) cout << "\n" << i->_name;
//    cout << "\n---------\n";

    Model NPP("NPP model");
    int n = _nodes.size();
    DebugOff("\nn = " << n);

    var<double> W("W");
    W._psd = true;
    NPP.add_var(W.in(R(2*n,2*n))); // note: number of instances in upper triangle is n*(2*n+1)
//#endif

    var<double> z("z");
    z.in_q_cone();
    auto Rn = R(n*n+1);
    NPP.add_var(z.in(Rn));

    func_ obj;
    obj = z("obj");
    NPP.min(obj);

//    min t
//    s.t.
//    ||z|| <= t    (z(obj) = t)
//    w*-w = z
//    W is PSD,
//    L <= W <= U
//    matr structure

    /** Constraints **/

    /* w*-w = z */
    Constraint svec("svec");

    auto idxs_ul_u = ul_u(n);
//    cout << "\nUpper left, upper: ";
//    for(auto& idx: idxs_ul_u) DebugOff(idx._name << "; ");

    auto idxs_lr_u = lr_u(n);
//    cout << "\nLower right, upper: ";
//    for(auto& idx: idxs_lr_u) DebugOff(idx._name << "; ");

    auto idxs_ur_u = ur_u(n);
//    cout << "\nUpper right, upper: ";
//    for(auto& idx: idxs_ur_u) DebugOff(idx._name << "; ");

    auto idxs_ur_l = ur_l(n);
//    cout << "\nUpper right, lower: ";
//    for(auto& idx: idxs_ur_l) DebugOff(idx._name << "; ");

    vector<index_> obj_idxs;
    obj_idxs.insert(obj_idxs.end(),idxs_ul_u.begin(),idxs_ul_u.end());
    obj_idxs.insert(obj_idxs.end(),idxs_ur_u.begin(),idxs_ur_u.end());
    DebugOff("\nObj: ");
    for(auto& idx: obj_idxs) DebugOn(idx._name << "; ");

    auto idxs = triang_indices(2*n);
    svec = _W_star - W - z;
    NPP.add_constraint(svec.in(obj_idxs)==0);
//    NPP.add_constraint(svec.in(idxs)==0);
//    svec.print_expanded();

    vector<index_> zero_idxs;
    for(int i = 0; i < n; i++) zero_idxs.push_back(index_(to_string(i)+","+to_string(i+n)));

    Constraint zeros("zeros");
    zeros = W.in(zero_idxs);
    NPP.add_constraint(zeros==0);
//    zeros.print_expanded();

    Constraint real_symm("real_symm");
    real_symm = W.submat(idxs_ul_u) - W.submat(idxs_lr_u);
    NPP.add_constraint(real_symm==0);
//    real_symm.print_expanded();

    Constraint imag_symm("imag_symm");
    imag_symm = W.submat(idxs_ur_u) + W.submat(idxs_ur_l);
    NPP.add_constraint(imag_symm==0);
//    imag_symm.print_expanded();

    Constraint W_lb("W_lb");
    W_lb = _Wmin - W;
    NPP.add_constraint(W_lb.in(obj_idxs) <= 0);
//    W_lb.print_expanded();

    Constraint W_ub("W_ub");
    W_ub = _Wmax - W;
    NPP.add_constraint(W_ub.in(obj_idxs) >= 0);
//    W_ub.print_expanded();

    string namew, namewr, namewi;

    solver s(NPP,Mosek);
    s.run(0,0);
//    exit(0);

//    z.print(); cout << "\n";
//    W.print(); cout << "\n";
//    w.print(); cout << "\n";

    for(int i = 0; i < n; i++){
        for(int j = i; j < n; j++){
            if(i==j){
                namew = "w(" + _nodes[i]->_name + ")";
                what.set_val(namew,W(to_string(i)+","+to_string(j)).eval());
            }
            else {
                namewr = "wr(" + _nodes[i]->_name + "," + _nodes[j]->_name + ")";
                namewi = "wi(" + _nodes[i]->_name + "," + _nodes[j]->_name + ")";
                what.set_val(namewr, W(to_string(i)+","+to_string(j)).eval());
                if(_grid->get_directed_arc(_nodes[i]->_name,_nodes[j]->_name)!=nullptr)
                    what.set_val(namewi, W(to_string(i)+","+to_string(j+_nodes.size())).eval());
                else
                    what.set_val(namewi, -W(to_string(i)+","+to_string(j+_nodes.size())).eval());
            }
        }
    }

    what.set_name("w_hat");
//    what.print(true);
    return what;
}
#else
/*--------------------------------------------------------------------------------*/
param<double> Bag::nfp(){
    param<double> what;
    fill_wstar();

    //    cout << "\nIndices:";
    //    for(auto& i: _indices) cout << "\n" << i->_name;
    //    cout << "\n---------\n";

    Model NPP("NPP model");
    int n = _nodes.size();
    DebugOff("\nn = " << n);

    //    _wstarp.print(true);

//#ifndef FLAT
    var<double> W("W");
    W._psd = true;
//    W._is_matrix = true;
    NPP.add_var(W.in(R(2*n,2*n)));
//#else
//    var<double> W("W");
//    W._psd = true;
//    NPP.add_var(W.in(R(2*n,2*n))); // note: number of instances in apper triangle is n*(2*n+1)
//#endif

    var<double> z("z");
    z.in_q_cone();
//#ifndef FLAT
    auto Rn = R(n*n+1);
//#else
//    auto Rn = R(2*n*2*n+1);
//#endif
    NPP.add_var(z.in(Rn));

    var<double> w("w", _wmin, _wmax);
    NPP.add_var(w.in(_indices));

    func_ obj;
    obj = z("obj");
    NPP.min(obj);

    //    min t
    //    s.t.
    //    ||z|| <= t    (z(obj) = t)
    //    w*-w = z
    //    W is PSD,
    //    L <= W <= U
    //    matr structure

    /** Constraints **/

    /* w*-w = z */
    Constraint svec("svec");

//#ifndef FLAT
    svec = _wstarp - w - z;
    NPP.add_constraint(svec.in(_indices)==0);
//#else
//    auto idxs = triang_indices(2*n);
//    svec = _W_star - W - z;
//    NPP.add_constraint(svec.in(idxs)==0);
//    //    svec.print_expanded();
//

//    auto idxs_ul_u = ul_u(n);
//    cout << "\nUpper left, upper: ";
//    for(auto& idx: idxs_ul_u) DebugOff(idx._name << "; ");

//    auto idxs_lr_u = lr_u(n);
//    cout << "\nLower right, upper: ";
//    for(auto& idx: idxs_lr_u) DebugOff(idx._name << "; ");

//    auto idxs_ur_u = ur_u(n);
//    cout << "\nUpper right, upper: ";
//    for(auto& idx: idxs_ur_u) DebugOff(idx._name << "; ");

//    auto idxs_ur_l = ur_l(n);

//    Constraint real_symm("real_symm");
//    real_symm = W.submat(idxs_ul_u) - W.submat(idxs_lr_u);
//    NPP.add_constraint(real_symm==0);
//    real_symm.print_expanded();

//    Constraint imag_symm("imag_symm");
//    imag_symm = W.submat(idxs_ur_u) + W.submat(idxs_ur_l);
//    NPP.add_constraint(imag_symm==0);

//    vector<index_> zero_idxs;
//    for(int i = 0; i < n; i++) {
//        zero_idxs.push_back(index_(to_string(i)+","+to_string(i+n)));
//    }
//    Constraint zeros("zeros");
//    zeros = W.in(zero_idxs);
//    NPP.add_constraint(zeros==0);
//    zeros.print_expanded();

//    Constraint W_lb("W_lb");
//    W_lb = _Wmin - W;
//    NPP.add_constraint(W_lb.in(idxs) <= 0);
//    //    W_lb.print_expanded();

//    Constraint W_ub("W_ub");
//    W_ub = _Wmax - W;
//    NPP.add_constraint(W_ub.in(idxs) >= 0);
//    //    W_ub.print_expanded();

//#endif

    string namew, namewr, namewi;

//#ifndef FLAT
    /* matrix structure */
    for(int i = 0; i < n; i++){
        for(int j = i; j < n; j++){
            if(i==j){
                namew = "w(" + _nodes[i]->_name + ")";

                Constraint mstruct1("mstruct"+to_string(i+n)+to_string(i+n));
                mstruct1 = w(namew) - W(to_string(i+n)+","+to_string(i+n));
                NPP.add_constraint(mstruct1==0);

                Constraint mstruct2("mstruct"+to_string(i)+to_string(i));
                mstruct2 = w(namew) - W(to_string(i)+","+to_string(i));
                NPP.add_constraint(mstruct2==0);

                /* zeros */
                Constraint zero("zero"+to_string(i)+to_string(i+n));
                zero = W(to_string(i)+","+to_string(i+n));
                NPP.add_constraint(zero==0);
                Constraint zero2("zero"+to_string(i+n)+to_string(i));
                zero2 = W(to_string(i+n)+","+to_string(i));
                NPP.add_constraint(zero2==0);
            }
            else {
                namewr = "wr(" + _nodes[i]->_name + "," + _nodes[j]->_name + ")";

                Constraint mstruct("mstructr"+to_string(i)+to_string(j));
                mstruct = w(namewr) - W(to_string(i)+","+to_string(j));
                NPP.add_constraint(mstruct==0);

                Constraint mstruct2("mstructr"+to_string(i+n)+to_string(j+n));
                mstruct2 = w(namewr) - W(to_string(i+n)+","+to_string(j+n));
                NPP.add_constraint(mstruct2==0);

                namewi = "wi(" + _nodes[i]->_name + "," + _nodes[j]->_name + ")";

                Constraint mstructi("mstructi"+to_string(i)+to_string(j+n));
                mstructi = w(namewi) - W(to_string(i)+","+to_string(j+n));
                NPP.add_constraint(mstructi==0);
                Constraint mstructi2("mstructi"+to_string(j)+to_string(i+n));
                mstructi2 = w(namewi) + W(to_string(j)+","+to_string(i+n));
                NPP.add_constraint(mstructi2==0);
            }
        }
    }
//#endif

//    zeros.print_expanded();

    solver s(NPP,Mosek);
    s.run(0,0);

    //    z.print(); cout << "\n";
    //    W.print(true); cout << "\n";
    //    w.print(); cout << "\n";
//#ifndef FLAT
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
//#else
//    for(int i = 0; i < n; i++){
//        for(int j = i; j < n; j++){
//            if(i==j){
//                namew = "w(" + _nodes[i]->_name + ")";
//                what.set_val(namew,W(to_string(i)+","+to_string(j)).eval());
//            }
//            else {
//                namewr = "wr(" + _nodes[i]->_name + "," + _nodes[j]->_name + ")";
//                namewi = "wi(" + _nodes[i]->_name + "," + _nodes[j]->_name + ")";
//                what.set_val(namewr, W(to_string(i)+","+to_string(j)).eval());
//                if(_grid->get_directed_arc(_nodes[i]->_name,_nodes[j]->_name)!=nullptr)
//                    what.set_val(namewi, W(to_string(i)+","+to_string(j+_nodes.size())).eval());
//                else
//                what.set_val(namewi, -W(to_string(i)+","+to_string(j+_nodes.size())).eval());
//            }
//        }
//    }
//#endif
    what.set_name("w_hat");
    //    what.print(true);
    //    exit(0);
    return what;
}
//param<double> Bag::nfp(){
//    param<double> what;
//    fill_wstar();
//
//    cout << "\nIndices:";
//    for(auto& i: _indices) cout << "\n" << i._name;
//    cout << "\n---------\n";
//
//    Model NPP("NPP model");
//    int n = _nodes.size();
//    DebugOff("\nn = " << n);
//
////    _wstarp.print(true);
//
//    sdpvar<double> W("W");
//    NPP.add_var(W^(2*n));
//
//    var<double> z("z");
//    z.in_q_cone();
//    auto Rn = R(n*n+1);
//    NPP.add_var(z.in(Rn));
//
//    var<double> w("w");//, _wmin, _wmax);
//    NPP.add_var(w.in(_indices));
//
//    func_ obj;
//    obj = z("obj");
//    NPP.min(obj);
//
////    min t   s.t.
////    ||z|| <= t    (z(obj) = t)
////    w*-w = z
////    W is PSD,
////    L <= W <= U
////    matr structure
//
//    /** Constraints **/
//
//    /* w*-w = z */
//    Constraint svec("svec");
//
//    svec = _wstarp - w - z;
//    NPP.add_constraint(svec.in(_indices)==0);
//
//    string namew, namewr, namewi;
//
//    /* matrix structure */
//    for(int i = 0; i < n; i++){
//        for(int j = i; j < n; j++){
//            if(i==j){
//                namew = "w(" + _nodes[i]->_name + ")";
//
//                Constraint mstruct1("mstruct"+to_string(i+n)+to_string(i+n));
//                mstruct1 = w(namew) - W(i+n,i+n);
//                NPP.add_constraint(mstruct1==0);
//                mstruct1.print_expanded();
//
//                Constraint mstruct2("mstruct"+to_string(i)+to_string(i));
//                mstruct2 = w(namew) - W(i,i);
//                NPP.add_constraint(mstruct2==0);
//
//                /* zeros */
//                Constraint zero("zero"+to_string(i)+to_string(i+n));
//                zero = W(i,i+n);
//                NPP.add_constraint(zero==0);
//                Constraint zero2("zero"+to_string(i+n)+to_string(i));
//                zero2 = W(i+n,i);
//                NPP.add_constraint(zero2==0);
//            }
//            else {
//                namewr = "wr(" + _nodes[i]->_name + "," + _nodes[j]->_name + ")";
//
//                Constraint mstruct("mstructr"+to_string(i)+to_string(j));
//                mstruct = w(namewr) - W(i,j);
//                NPP.add_constraint(mstruct==0);
//                Constraint mstruct2("mstructr"+to_string(i+n)+to_string(j+n));
//                mstruct2 = w(namewr) - W(i+n,j+n);
//                NPP.add_constraint(mstruct2==0);
//
//                namewi = "wi(" + _nodes[i]->_name + "," + _nodes[j]->_name + ")";
//
//                Constraint mstructi("mstructi"+to_string(i)+to_string(j+n));
//                mstructi = w(namewi) - W(i,j+n);
//                NPP.add_constraint(mstructi==0);
//                Constraint mstructi2("mstructi"+to_string(j)+to_string(i+n));
//                mstructi2 = w(namewi) + W(j,i+n);
//                NPP.add_constraint(mstructi2==0);
//            }
//        }
//    }
//
//    solver s(NPP,Mosek);
//    s.run(0,0);
//
////    z.print(); cout << "\n";
////    W.print(); cout << "\n";
////    w.print(); cout << "\n";
//
//    for(int i = 0; i < n; i++){
//        for(int j = i; j < n; j++){
//            if(i==j){
//                namew = "w(" + _nodes[i]->_name + ")";
//                what.set_val(namew,w(namew).eval());
//            }
//            else {
//                namewr = "wr(" + _nodes[i]->_name + "," + _nodes[j]->_name + ")";
//                namewi = "wi(" + _nodes[i]->_name + "," + _nodes[j]->_name + ")";
//                what.set_val(namewr, w(namewr).eval());
//                if(_grid->get_directed_arc(_nodes[i]->_name,_nodes[j]->_name)!=nullptr)
//                    what.set_val(namewi, w(namewi).eval());
//                else
//                    what.set_val(namewi, -w(namewi).eval());
//            }
//        }
//    }
//
////    double w1, w2, w3, wr12, wr13, wr23, wi12, wi13, wi23;
////    w1 = W("0,0").eval(); w2 = W("1,1").eval(); w3 = W("2,2").eval();
////    wr12 = W("0,1").eval(); wr13 = W("0,2").eval(); wr23 = W("1,2").eval();
////    wi12 = W("0,4").eval(); wi13 = W("0,5").eval(); wi23 = W("1,5").eval();
////
////    double SDP = wr12*(wr23*wr13 + wi23*wi13) + wi12*(-wi23*wr13 + wr23*wi13);
////    SDP *= 2;
////    SDP -= (wr12*wr12 + wi12*wi12)*w3 + (wr13*wr13 + wi13*wi13)*w2 + (wr23*wr23 + wi23*wi23)*w1;
////    SDP += w1*w2*w3;
////    DebugOn("\nSDP of w_hat = " << SDP);
////    DebugOn("\nSOCP12 = " << wr12*wr12 + wi12*wi12 - w1*w2);
////    DebugOn("\nSOCP13 = " << wr13*wr13 + wi13*wi13 - w1*w3);
////    DebugOn("\nSOCP23 = " << wr23*wr23 + wi23*wi23 - w3*w2);
//
//    what.set_name("w_hat");
////    what.print(true);
//    return what;
//}

#endif
