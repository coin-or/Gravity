//
//  Test.cpp
//  
//
//  Created by Hassan on 3 Jan 2016.
//
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <stdio.h>
#include <iostream>
#include <string>
#include <stdio.h>
#include <cstring>
#include <fstream>
//#include <gravity/Net.h>
#include <gravity/model.h>
#include <gravity/solver.h>
#include <gravity/func.h>
#include <stdio.h>
#include <stdlib.h>
#include <gravity/doctest.h>
//#include <PowerNet.h>
//#include <variant>


using namespace std;
using namespace gravity;

//
//TEST_CASE("testing variants") {
//    class type{
//    public:
//        variant<bool, short, int, float, double> val;
//    };
//    class cst{
//    public:
//        type _val;
//        cst& operator=(double v){
//            _val.val = v;
//            return *this;
//        }
//
//        cst& operator=(int v){
//            _val.val = v;
//            return *this;
//        }
//    };
//
//    cst c;
//    c = 1.2;
//    cst c2;
//    c2 = 3;
//    cout << "c = " << get<double>(c._val.val) << endl;
//    cout << "c2 = " << get<int>(c2._val.val) << endl;
//    cout << "c2 = " << get<int>(c2._val.val) << endl;
//
//
////    v = true; // v contains bool
////    w = 1.;
////    int pos = v.index();
//////    auto alt = variant_alternative<0,bool>(v);
////    cout << "type of v = " << typeid(get<int>(v)).name() << endl;
////    cout << "type of w = " << typeid(w).name() << endl;
////    cout << "type of v = " << typeid(v).name() << endl;
//}
//
//TEST_CASE("testing something") {
//    constant<int> ic;
//    constant<bool> bc;
//    shared_ptr<constant_> icp = make_shared<constant<int>>(ic);
//    shared_ptr<constant_> bcp = make_shared<constant<bool>>(bc);
//    cout << "type of ic = " << typeid(ic).name() << endl;
//    cout << "type of icp = " << typeid(*icp.get()).name() << endl;
//    cout << "type of bc = " << typeid(bc).name() << endl;
//    cout << "type of bcp = " << typeid(*bcp.get()).name() << endl;
//}

TEST_CASE("testing constants") {
    constant<> c0;
    CHECK(c0.is_double());
    c0 = 3.5;
    CHECK(c0==3.5);
    CHECK(c0.is_positive());
    CHECK(!c0.is_negative());
    constant<Cpx> cx0;
    cx0 = Cpx(-1,1);
    CHECK(cx0.is_complex());
    CHECK(cx0==Cpx(-1,1));
    constant<Cpx> cx1;
    cx1 = Cpx(-1,-2);
    auto cx2 = cx0 + cx1;
    cx2.print();
    CHECK(cx2.is_complex());
    CHECK(cx2.is_negative());
    auto mag0 = sqrmag(cx2);
    CHECK(abs(mag0.eval()-5)<1e-8);
    auto ang0 = angle(cx2);
    ang0.println();
    CHECK(abs(ang0.eval()-(-2.677945045))<1e-8);
    auto cx3 = conj(cx1);
    CHECK(cx3.eval().real()==cx1.eval().real());
    CHECK(cx3.eval().imag()==-1*cx1.eval().imag());
    CHECK(real(cx3) == cx3.eval().real());
    CHECK(imag(cx3) == cx3.eval().imag());
    CHECK(cx3.get_dim() == 1);
    try{
        cx3.get_dim(2);
    }
    catch(invalid_argument& arg){
        cout << "Error successfully caught: "<< endl;
        cout << arg.what() << endl;
    }
    CHECK(cx3.is_number());
}

TEST_CASE("testing parameters") {
    param<> c0("c0");
    CHECK(c0.is_double());
    c0 = 3.5;
    CHECK(c0.is_positive());
    CHECK(!c0.is_negative());
    param<Cpx> cx0("p0");
    cx0 = Cpx(-1,1);
    CHECK(cx0.is_complex());
    param<Cpx> cx1("p1");
    cx1 = Cpx(-1,-2);
    auto mag0 = sqrmag(cx1);
    auto ang0 = ang(cx1);
    mag0.print();
    ang0.print();
    auto cx3 = conj(cx1);
    CHECK(cx3.eval().real()==cx1.eval().real());
    CHECK(cx3.get_dim() == 1);
}

TEST_CASE("testing param copy operator") {
    param<int> ip("ip");
    ip.print();
    ip.add_val(2);
    ip.add_val(3);
    ip.print();
    param<int> ip2(ip);
    ip2.print();
    CHECK(ip==ip2);
}
//
TEST_CASE("testing param indexing, add_val() and set_val() functions") {
    param<> ip("ip");
    ip.print();
    ip.add_val(2);
    ip.add_val(-1.3);
    ip.print();
    ip.set_val(0, 1.5);
    ip.print();
    CHECK(ip.eval(0)==1.5);
    indices ids("index_set");
    ids = {"id1", "id2", "key3"};
    param<> dp("dp");
    dp.in(ids);
    dp.print();
    dp("id1") = 1.5;
    dp("id2") = -231.5;
    dp.print();
    CHECK(dp.eval("id1")==1.5);
    CHECK(dp.eval("id2")==-231.5);
    CHECK(dp.eval("key3")==0);
    REQUIRE_THROWS_AS(dp("unexisting_key").eval(), invalid_argument);
    auto ndp = dp.in(ids.exclude("id2"));
    ndp.print();
    CHECK(ndp.get_dim()==2);
}

TEST_CASE("testing matrix params") {
    param<> mat("mat");
    mat.set_size(3,3);
    mat.set_val(1, 2, 2.3);
    mat.print();
    CHECK(mat.eval(1,2)==2.3);
    CHECK(mat(1,2).eval()==2.3);
    auto tr_mat = mat.tr();
    tr_mat.print();
    CHECK(tr_mat.eval(2,1)==2.3);
    CHECK(tr_mat(2,1).eval()==2.3);
    /* Complex matrices */
    param<Cpx> Cmat("Cmat");
    Cmat.set_size(3,3);
    Cpx cval = Cpx(-1,2);
    Cmat.set_val(0, 1, cval);
    Cmat.print();
    CHECK(Cmat.eval(0,1)==cval);
    CHECK(Cmat(0,1).eval()==cval);
}

TEST_CASE("testing param sign functions") {
    param<int> ip("ip");
    ip.print();
    ip.add_val(2);
    ip.add_val(3);
    ip.print();
    CHECK(ip.is_non_negative());
    CHECK(ip.is_positive());
    CHECK(!ip.is_negative());
    CHECK(!ip.is_non_positive());
    param<> dp("dp");
    dp.add_val(0);
    dp.add_val(-3);
    dp.print();
    CHECK(!dp.is_positive());
    CHECK(!dp.is_non_negative());
    CHECK(!dp.is_negative());
    CHECK(dp.is_non_positive());
}
//
TEST_CASE("testing variables indexing") {
    indices ids("index_set");
    ids = {"id1", "id2", "key3"};
    var<> iv("iv",-2, 5);
    iv.in(ids);
    iv.print();
    param<Cpx> lb("lb"), ub("ub");
    lb.in(ids);
    lb = Cpx(0,-1);
    ub.in(ids);
    ub = Cpx(1,1);
    var<Cpx> cv("cv", lb, ub);
    cv.in(ids.exclude("id2"));
    cv.print();
    CHECK(cv.get_dim()==2);
    var<Cpx> cv1("cv", Cpx(0,-1),Cpx(1,1));
    cv1.in(ids.exclude("id2"));
    cv1.print();
    CHECK(cv1.get_dim()==2);
    CHECK(cv==cv1);
    CHECK(cv.get_indices()!=ids);
}



TEST_CASE("testing vector dot product"){
    param<> a("a");
    a.set_size(3);
    a.set_val(0, 1);
    a.set_val(1, -1);
    a.set_val(2, 2);
    param<int> b("b");
    b=5;
    CHECK(b.get_dim()==1);
    CHECK(b.is_integer());
    CHECK(b.eval(0)==5);
    b=0;
    b=2;
    b=5;
    b=4;
    CHECK(b.get_dim()==5);
    CHECK(b.eval(1)==0);
    b.set_val(1, 3);
    CHECK(b.eval(1)==3);
    b.print();
    var<> z("z",-1,1);
    z.in(R(4));
    CHECK(z.get_dim()==4);
    z.print();
    z.in(R(3));
    CHECK(z.get_dim()==3);
    z.print();
    param<Cpx> cp("cp");
    cp = Cpx(-1,-1);
    auto lin = cp+z;
    CHECK(lin.is_complex());
    CHECK(lin.is_linear());
    auto lin2 = cp*z;
    auto lin3 = a.tr()*z - b;
    CHECK(lin3.is_double());
    CHECK(lin3.get_dim()==5);
    lin3.print();
    var<Cpx> y("y");
    y.in(C(5));
    auto lin4 = b*y*y;
    lin4.print();
    CHECK(lin4.is_complex());
    CHECK(lin4.is_convex());
    CHECK(lin4.get_dim()==5);
    auto lin5 = lin4 - lin3;
    CHECK(lin5.is_complex());
    CHECK(lin5.get_dim()==5);
    CHECK(lin5.get_nb_vars()==2);
    CHECK(lin5.is_convex());
    lin5.print_symbolic();
    lin5.print();
    CHECK(lin5.to_str()=="(b)y² - ([a]ᵀ)z + (b)");
    auto lin6 = (a+exp(b)).tr()*z;
    lin6.print_symbolic();
    CHECK(lin6.is_linear());
    lin6.print();
    param<Cpx> cpx("cpx");
    cpx = Cpx(1,1);
    cpx = Cpx(2,2);
    cpx = Cpx(3,1);
    auto cpx_f = 2*exp(cpx)*z;
    cpx_f.print_symbolic();
    CHECK(cpx_f.to_str()=="((2,0) * exp(cpx))z");
    cpx_f.print();
}
//    auto df = lin.get_dfdx(z);
//    CHECK(df.is_constant());
//    CHECK(df.get_dim()==3);
//    df.print();
//    auto lin2 = sum(z);
//    lin2 += 2*z;
//    lin2.print_symbolic();
//    lin2.print();
//    CHECK(lin2.get_nb_instances()==3);
//    auto dfdz = lin2.get_dfdx(z);
//    dfdz.print_symbolic();
//    CHECK(dfdz==2);
//    auto dfdvecz = lin2.get_dfdx(z.vec());
//    dfdvecz.print_symbolic();
//}
//
TEST_CASE("testing complex numbers") {
    indices ids("index_set");
    ids = {"id1", "id2", "key3", "key4"};
    var<> iv("x",-2, 5);
    iv.in(ids);
    var<Cpx> cv("y", Cpx(0,-1),Cpx(1,1));
    constant<Cpx> cx(Cpx(-1,-2));
    auto cx_conj = conj(cx);
    CHECK(real(cx_conj)==real(cx));
    CHECK(imag(cx_conj).eval()==-1*imag(cx).eval());
    param<Cpx> px("px");
    px = Cpx(-1,-2);
    px.print();
    auto px_conj = conj(px);
    CHECK(real(px_conj)._is_real);
    CHECK(imag(px_conj)._is_imag);
}
//
TEST_CASE("testing complex functions") {
    indices ids("index_set");
    ids = {"id1", "id2", "key3", "key4"};
    var<> iv("x",-2, 5);
    iv.in(ids);
    var<Cpx> cv("y", Cpx(0,-1),Cpx(1,1));
    cv.in(ids);
    auto f = 2*iv;
    f.print_symbolic();
    CHECK(f.is_linear());
    CHECK(f.is_convex());
    f+= pow(iv,2);
    f+=2;
    f.print_symbolic();
    f.print();
    CHECK(f.is_quadratic());
    CHECK(f.is_convex());
//    auto dfx = f.get_dfdx(iv);
//    dfx.print();
//    CHECK(dfx.is_linear());
//    CHECK(dfx==2*iv+2);
//    f *= iv;
//    f.print_symbolic();
//    CHECK(f.is_polynomial());
//    CHECK(!f.is_complex());
//    dfx = f.get_dfdx(iv);
//    dfx.print();
//    CHECK(dfx.is_quadratic());
//    CHECK(dfx==3*power(iv,2)+4*iv);
//    f += iv*log(cv);
//    f.print_symbolic();
//    CHECK(f.is_nonlinear());
//    CHECK(f.is_complex());
//    dfx = f.get_dfdx(iv);
//    dfx.print();
//    CHECK(dfx==3*power(iv,2)+4*iv + log(cv));
//    auto dfy = f.get_dfdx(cv);
//    dfy.print();
//    CHECK(dfy==iv*(1/cv));
//    f.in(ids.exclude("id2"));
//    CHECK(f.get_nb_vars()==2);
//    CHECK(f.get_nb_instances()==3);
}
//
//
TEST_CASE("testing complex matrix product") {
    var<Cpx> X("X", Cpx(0,-1),Cpx(1,1));
    X.in(C(5,5));
    param<Cpx> A("A");
    A.set_size(3, 5);
    A.set_val(0,0,Cpx(-1,1));
    A.set_val(1,0,Cpx(2,1));
    A.set_val(2,1,Cpx(-1,-1));
    A.set_val(1,1,Cpx(1,1));
    auto f = A*X;
    f.print_symbolic();
    CHECK(f.is_linear());
    CHECK(f.is_convex());
    CHECK(f.get_dim()==15);
    var<Cpx> Y("Y", Cpx(0,-1),Cpx(1,1));
    Y.in(C(5,5));
    f+= X*Y;
    f.print_symbolic();
    f += X;
    f.print_symbolic();
    CHECK(f.is_quadratic());
    CHECK(!f.is_convex());
    CHECK(f.is_complex());
    CHECK(f.get_nb_vars()==50);
    CHECK(f.to_str()=="[XY + (A + (1,0))X]");
}
//
TEST_CASE("testing function convexity"){
    var<> dp("dp",0.5,10.);
    var<int> ip("ip",-1,1);
    auto expr = log(dp) + sqrt(ip);
    expr.print_symbolic();
    CHECK(expr.is_concave());
    var<> p("p",0,1);
    var<> q("q");
    auto cc = p*p + q*q;
    cc.print_symbolic();
    CHECK(cc.is_convex());
    auto cc1 = cc * -1;
    cc1.print_symbolic();
    CHECK(cc1.is_concave());
    cc1 += 2*p*q;
    cc1.print_symbolic();
//    CHECK(cc1.is_rotated_soc());
    param<int> aa("aa");
    aa = -1;
    aa = -3;
    auto ff = (aa)*p*p;
    ff.print_symbolic();
    CHECK(ff.is_concave());
    ff *= aa;
    ff.print_symbolic();
    CHECK(ff.is_convex());
    ff *= -1;
    ff.print_symbolic();
    CHECK(ff.is_concave());
    ff *= aa;
    ff.print_symbolic();
    CHECK(ff.is_convex());
    ff += aa*(ip + dp)*q*q;
    ff.print_symbolic();
    CHECK(!ff.is_convex());
    CHECK(!ff.is_concave());
    CHECK(ff.is_polynomial());
    CHECK(ff.get_nb_vars()==4);
    param<> b("b");
    b = 1;
    auto fn = p*p + q*q;
    fn.print_symbolic();
    CHECK(fn.is_convex());
    fn -= 2*p*q;
    fn.print_symbolic();
    fn += exp(p);
    fn.print_symbolic();
//    CHECK(fn.is_convex());
}
//
//TEST_CASE("testing nonlinear expressions"){
//    var<> x1("x1", -1, 1), x2("x2", 0, 3), x3("x3");
//    Constraint cstr("cycle");
//    cstr = cos(x1*x2);
//    CHECK(cstr.get_nb_vars()==2);
//    cstr += sin(x2*x3) + x1*expo(x2*x3) + log(x2);
//    CHECK(cstr.get_nb_vars()==3);
//}
//
//TEST_CASE("testing monomials"){
//    var<> x1("x1"), x2("x2"), x3("x3");
//    Constraint cstr("cycle");
//    cstr = x1*x2 + x2*x3 + x1*x3;
//    auto monoms = cstr.get_monomials(9);
//}
//
//TEST_CASE("testing acopf"){
//    string fname = string(prj_dir)+"/data_sets/Power/nesta_case5_pjm.m";
//    unsigned nb_threads = 2;
//    double tol = 1e-6;
//    string mehrotra = "no", log_level="0";
//    PowerNet grid1,grid2;
//    grid1.readgrid(fname);
//    auto ACOPF1 = grid1.build_ACOPF(ACRECT);
//    fname = string(prj_dir)+"/data_sets/Power/nesta_case14_ieee.m";
//    auto ACOPF2 = grid1.build_ACOPF(ACPOL);
//    auto models = {ACOPF1, ACOPF2};
//    /* run in parallel */
//    run_parallel(models, ipopt, tol = 1e-6, nb_threads=2);
//    CHECK(abs(ACOPF1->_obj_val-17551.8909275818)<tol);
//    CHECK(ACOPF1->is_feasible(tol));
//    ACOPF1->print_solution();
//    auto Mc = ACOPF1->build_McCormick();
//    auto clone = grid1.clone();
//    CHECK(clone->arcs.size()==grid1.arcs.size());
//    clone->remove_arc(clone->arcs.at(0));
//    CHECK(clone->arcs.size()==grid1.arcs.size()-1);
//}
//
//TEST_CASE("testing socopf"){
//    auto time_start = get_cpu_time();
//    string fname = string(prj_dir)+"/data_sets/Power/nesta_case5_pjm.m";
//    int output = 0;
//    bool relax = false;
//    double tol = 1e-6;
//    string mehrotra = "no", log_level="0";
//    PowerNet grid;
//    grid.readgrid(fname);
//    auto SOCOPF = grid.build_SCOPF();
//    solver OPF(*SOCOPF,ipopt);
//    OPF.run(output, relax = false, tol = 1e-6);
//    auto time_end = get_cpu_time();
//    DebugOn("Total cpu time = " << time_end - time_start << " secs" << endl);
//    CHECK(abs(SOCOPF->_obj_val-14999.715037743885)<tol);
//}
