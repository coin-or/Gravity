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
#include <functional>
#include <stdio.h>
#include <cstring>
#include <fstream>
//#include <gravity/Net.h>
#include <gravity/model.h>
#include <gravity/solver.h>
#include <gravity/constraint.h>
#include <stdio.h>
#include <stdlib.h>
#include <gravity/doctest.h>
//#include <PowerNet.h>
//#include <variant>


using namespace std;
using namespace gravity;


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
    param<Cpx> cx0("cx0");
    cx0 = Cpx(-1,1);
    CHECK(cx0.is_complex());
    param<Cpx> cx1("cx1");
    cx1 = Cpx(-1,-2);
    auto mag0 = sqrmag(cx1);
    auto ang0 = ang(cx1);
    mag0.print();
    ang0.print();
    auto cx2 = conj(cx1);
    CHECK(cx2.eval().real()==cx1.eval().real());
    CHECK(cx2.get_dim() == 1);
    param<Cpx> cx3("cx3");
    cx3.in(C(6));
    CHECK(cx3.get_dim() == 6);
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
    CHECK(ip2.is_param());
    CHECK(ip2.get_intype()==integer_);
    CHECK(ip2.is_positive());
    CHECK(ip2._range->first==2);
    CHECK(ip2._range->second==3);
}

TEST_CASE("testing param indexing, add_val() and set_val() functions") {
    param<> ip("ip");
    ip.print();
    ip.add_val(2);
    ip.add_val(-1.3);
    CHECK(ip._range->first==-1.3);
    CHECK(ip._range->second==2);
    ip.print();
    ip.set_val(1, 1.5);
    ip.print();
    CHECK(ip.eval(1)==1.5);
    CHECK(ip._range->first==1.5);
    CHECK(ip._range->second==2);
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
    CHECK(dp._range->first==-231.5);
    CHECK(dp._range->second==1.5);
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
    cp.in(C(3));
    cp.set_val(Cpx(1,1));
    auto lin = cp+z;
    CHECK(!lin.is_evaluated());
    lin.print_symbolic();
    lin.print();
    CHECK(lin._range->first==Cpx(0,1));
    CHECK(lin._range->second==Cpx(2,1));
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
    lin4.print_symbolic();
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

TEST_CASE("testing complex functions") {
    var<> z("z",-1,1);
    z.in(R(3));
    param<Cpx> cpx("cpx");
    cpx = Cpx(1,1);
    cpx = Cpx(2,2);
    cpx = Cpx(3,1);
    auto cpx_f = 2*exp(cpx)*z;
    cpx_f.print_symbolic();
    CHECK(cpx_f.to_str()=="((2,0) * exp(cpx))z");
    cpx_f.print();
}

TEST_CASE("testing ReLU") {
    var<> A("A",-1,1), B("B",-1,1);
    A.in(R(4));
    B.in(R(1));
    param<> X("X");
    X.in(R(4));
    X(0) = 0.2;
    X(1) = 0.3;
    X(2) = 0.4;
    X(3) = -0.2;
    auto f = ReLU(A.tr()*X + B);
    f.print_symbolic();
    f.print();
    auto dfdA = f.get_derivative(B);
    dfdA.print_symbolic();
    CHECK(dfdA.to_str()=="UnitStep(B + (Xᵀ)[A])");
    dfdA.print();
}

TEST_CASE("testing range propagation") {
    indices ids("index_set");
    ids = {"id1", "id2", "key3", "key4"};
    var<> x("x",-2, 5);
    x.in(ids);
    auto f = 2*x + 2;
    f.print_symbolic();
    CHECK(f.is_linear());
    CHECK(f.is_convex());
    CHECK(f._range->first==-2);
    CHECK(f._range->second==12);
    f+= pow(x,2);
    CHECK(f._range->first==-2);
    CHECK(f._range->second==37);
    f+=2;
    CHECK(f._range->first==0);
    CHECK(f._range->second==39);
    CHECK(f.to_str()=="x² + 2x + 4");
    f.print_symbolic();
    f.print();
    CHECK(f.is_quadratic());
    CHECK(f.is_convex());
    auto dfx = f.get_derivative(x);
    dfx.print();
    CHECK(dfx.is_linear());
    CHECK(dfx.get_dim()==4);
    param<> a("a");
    a.in(R(4));
    a = 2;
    a.print();
    var<> x1("x1",-2, 2);
    x1.in(ids);
    var<> x2("x2",0, 3);
    x2.in(ids);
    auto f2 = a*pow(x1,4) * pow(x2,2);
    f2.print_symbolic();
    f2.print();
    CHECK(f2.is_polynomial());
    CHECK(f2._range->first==0);
    CHECK(f2._range->second==2*pow(2,4)*pow(3,2));
    var<> x3("x3",-2, 2);
    x3.in(ids);
    f2 *= pow(x3,3);
    f2.print_symbolic();
    f2.print();
    CHECK(f2._range->first==-2*pow(2,4)*pow(3,2)*pow(2,3));
    CHECK(f2._range->second==2*pow(2,4)*pow(3,2)*pow(2,3));
    auto dfdx1 = f2.get_derivative(x1);
    dfdx1.print_symbolic();
    dfdx1.print();
    CHECK(dfdx1._range->first==-8*pow(2,3)*pow(3,2)*pow(2,3));
    CHECK(dfdx1._range->second==8*pow(2,3)*pow(3,2)*pow(2,3));
}

TEST_CASE("testing polynomial functions") {
    param<int> a("a");
    a.set_size(3);
    a.set_val(0,1);
    a.set_val(1,-1);
    a.set_val(2,2);
    param<int> b("b");
    b=5;
    b=0;
    b=2;
    var<> x("x",-1,1);
    x.in(R(3));
    var<> y("y",-1,1);
    y.in(R(3));
    var<> z("z",-1,1);
    z.in(R(3));
    auto poly = pow(x,2)*pow(y,3)*pow(z,4) + b*pow(y,2)*pow(z,3);
    CHECK(poly.to_str()=="x²y³z⁴ + (b)y²z³");
    poly.print_symbolic();
    poly.print();
    poly += a*x;
    poly.print_symbolic();
    poly.print();
    auto dfdx = poly.get_derivative(x);
    CHECK(dfdx.to_str()=="2xy³z⁴ + a");
    dfdx.print_symbolic();
    dfdx.print();
    auto dfd2x = dfdx.get_derivative(x);
    CHECK(dfd2x.to_str()=="2y³z⁴");
    dfd2x.print_symbolic();
    dfd2x.print();
}

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
    f.print();
    CHECK(f.is_linear());
    CHECK(f.is_convex());
    CHECK(f.get_dim()==15);
    var<Cpx> Y("Y", Cpx(0,-1),Cpx(1,1));
    Y.in(C(5,5));
    f+= A*Y;
    f.print_symbolic();
    f.print();
    CHECK(f.is_complex());
    CHECK(f.get_nb_vars()==50);
    auto f2 = X*X;
    f2.print_symbolic();
    f2.print();
    CHECK(f2.is_convex());
    f2 -= X*Y;
    f2.print_symbolic();
    f2.print();
    CHECK(!f2.is_convex());
    CHECK(f2.get_dim()==25);
    auto f3 = exp(A*X);
    f3.print_symbolic();
    f3.print();
}

TEST_CASE("testing function convexity"){
    var<> dp("dp",0.5,10.);
    var<int> ip("ip",3,4);
    auto expr = log(dp) + sqrt(ip);
    expr.print_symbolic();
    CHECK(expr.is_concave());
    CHECK(expr._range->first==log(0.5)+(int)sqrt(3));
    CHECK(expr._range->second==log(10)+(int)sqrt(4));
    CHECK(expr.is_positive());
    var<> p("p",0,1);
    var<> q("q",-0.5, 0.5);
    auto cc = p*p + q*q;
    cc.print_symbolic();
    CHECK(cc.is_convex());
    auto cc1 = cc * -1;
    cc1.print_symbolic();
    CHECK(cc1.is_concave());
    cc1 += expr;
    cc1.print_symbolic();
    CHECK(cc1.is_concave());
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
    CHECK(ff._range->first==-10.5);
    CHECK(ff._range->second==37.5);
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
    fn += exp(p);
    fn.print_symbolic();
    CHECK(fn.is_convex());
    auto f2 = pow(p,4);
    CHECK(f2.is_convex());
}

TEST_CASE("testing constraints"){
    var<> x1("x1", -1, 1), x2("x2", 0.1, 3), x3("x3",2,4);
    x1.in(R(2));
    x2.in(R(2));
    x3.in(R(2));
    Constraint<> cstr("cstr");
    cstr = x1 + exp(x2) - x3;
    CHECK(cstr.get_nb_vars()==3);
    CHECK(cstr.is_nonlinear());
    CHECK(cstr.is_convex());
    CHECK(cstr.get_dim()==2);
    CHECK(cstr._range->first==-1+exp(0.1)-4);
    CHECK(cstr._range->second==1+exp(3)-2);
    cstr.print_symbolic();
    cstr.print();
    auto dfdx2 = cstr.get_derivative(x2);
    dfdx2.print_symbolic();
    param<int> a("a");
    a = -1;
    a = -4;
    Constraint<> cstr1("cstr1");
    cstr1 = x3 - sqrt(x2) + ReLU(x1);
    cstr1 >= a + 3;
    cstr1.print_symbolic();
    cstr1.print();
    CHECK(cstr1.is_concave());
    CHECK(cstr1._range->first==2-sqrt(3)-3+1);
    CHECK(cstr1._range->second==4-sqrt(0.1)+2);
}

TEST_CASE("testing soc/rotated soc constraints"){
    var<> x1("x1", -1, 1), x2("x2", 0.1, 3), x3("x3",2,4);
    x1.in(R(2));
    x2.in(R(2));
    x3.in(R(2));
    param<> a("a");
    a = 1;
    a = 4;
    Constraint<> cstr("cstr");
    cstr = (x2+a)*(x3+log(a)) - pow(x1,2);
    cstr >= 0;
    cstr.print_symbolic();
    cstr.print();
    CHECK(cstr.get_nb_vars()==3);
    CHECK(cstr.is_quadratic());
    CHECK(cstr.check_rotated_soc());
    CHECK(cstr.get_dim()==2);
    CHECK(cstr._range->first==1.1*(2+log(1))-1);
    CHECK(cstr._range->second==(3+4)*(4+log(4)));
    auto dfdx2 = cstr.get_derivative(x2);
    dfdx2.print_symbolic();
    CHECK(dfdx2.to_str()=="x3 + log(a)");
    CHECK(dfdx2.is_linear());
    CHECK(dfdx2._range->first==2+log(1));
    CHECK(dfdx2._range->second==4+log(4));
    Constraint<> cstr2("cstr2");
    cstr2 = pow(x2+a,2) + pow(x3+sqrt(a),2) - pow(x1,2);
    cstr2 <= 0;
    cstr2.print_symbolic();
    cstr2.print();
    CHECK(cstr2.get_nb_vars()==3);
    CHECK(cstr2.is_quadratic());
    CHECK(cstr2.check_soc());
    CHECK(cstr2.get_dim()==2);
    CHECK(cstr2._range->first==1.1*1.1 + 3*3 - 1);
    CHECK(cstr2._range->second==7*7 + 6*6);
}

TEST_CASE("testing nonlinear expressions"){
    var<> x1("x1", -1, 1), x2("x2", 0.1, 3), x3("x3");
    x1.in(R(1));
    x2.in(R(1));
    x3.in(R(1));
    auto cstr = x1*exp(x2*x3);
    CHECK(cstr.get_nb_vars()==3);
    CHECK(cstr.is_nonlinear());
    CHECK(cstr.get_dim()==1);
    cstr.print_symbolic();
    cstr.print();
    auto dfdx2 = cstr.get_derivative(x2);
    dfdx2.print_symbolic();
}

TEST_CASE("testing monomials"){
    var<> x1("x1"), x2("x2"), x3("x3");
    auto cstr = x1*x2 + x2*x3 + x1*x3;
    auto monoms = cstr.get_monomials(5);
}

TEST_CASE("testing Model"){
    Model<> M("MyModel");
    var<> x("x", 1,2), y("y", 2, 4), z("z", -1, 2);
    M.add(x.in(R(2)),y.in(R(2)),z.in(R(2)));
    Constraint<> cstr1("cstr1");
    cstr1 = pow(x,2) + pow(y,2) - pow(z,2);
    M.add(cstr1 <= 0);
    M.max(sum(x));
    M.print_symbolic();
    M.print();
    CHECK(M.is_convex());
    CHECK(M.get_nb_vars()==6);
    CHECK(M.get_nb_cons()==2);
//    M.print();
}

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
