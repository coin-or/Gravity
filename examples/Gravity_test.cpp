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
#include <gravity/constant.h>
#include <stdio.h>
#include <stdlib.h>
#include <gravity/doctest.h>
//#include <PowerNet.h>


using namespace std;
using namespace gravity;


TEST_CASE("testing constants") {
    constant<> c0;
    CHECK(c0.is_double());
    c0 = 3.5;
    CHECK(c0==3.5);
    CHECK(c0.is_positive());
    constant<Cpx> cx0;
    cx0 = Cpx(-1,1);
    CHECK(cx0.is_complex());
    CHECK(cx0==Cpx(-1,1));
    constant<Cpx> cx1;
    cx1 = Cpx(-1,-2);
    auto cx2 = cx0 + cx1;
    CHECK(cx2.is_complex());
    CHECK(cx2.is_negative());
    cx2.print();
}

//TEST_CASE("testing param copy operator") {
//    param<int> ip("ip");
//    ip.print();
//    ip.add_val(2);
//    ip.add_val(3);
//    ip.print();
//    param<int> ip2(ip);
//    ip2.print();
//    CHECK(ip==ip2);
//}
//
//TEST_CASE("testing param indexing, add_val() and set_val() functions") {
//    param<> ip("ip");
//    ip.print();
//    ip.add_val(2);
//    ip.add_val(-1.3);
//    ip.print();
//    ip.set_val(0, 1.5);
//    ip.print();
//    CHECK(ip.eval(0)==1.5);
//    indices ids("index_set");
//    ids = {"id1", "id2", "key3"};
//    param<> dp("dp");
//    dp.in(ids);
//    dp.print();
//    dp("id1") = 1.5;
//    dp("id2") = -231.5;
//    dp.print();
//    CHECK(dp.eval("id1")==1.5);
//    CHECK(dp.eval("id2")==-231.5);
//    CHECK(dp.eval("key3")==0);
//    REQUIRE_THROWS_AS(dp("unexisting_key").eval(), invalid_argument);
//    auto ndp = dp.in(ids.exclude("id2"));
//    ndp.print();
//    CHECK(ndp.get_dim()==2);
//
//}
//TEST_CASE("testing matrix params") {
//    param<> mat("mat");
//    mat.set_size(3,3);
//    mat.set_val(1, 2, 2.3);
//    mat.print();
//    CHECK(mat.eval(1,2)==2.3);
//    CHECK(mat(1,2).eval()==2.3);
//    auto tr_mat = mat.tr();
//    tr_mat.print();
//    CHECK(tr_mat.eval(2,1)==2.3);
//    CHECK(tr_mat(2,1).eval()==2.3);
//    /* Complex matrices */
//    param<Cpx> Cmat("Cmat");
//    Cmat.set_size(3,3);
//    Cpx cval = Cpx(-1,2);
//    Cmat.set_val(0, 1, cval);
//    Cmat.print();
//    CHECK(Cmat.eval(0,1)==cval);
//    CHECK(Cmat(0,1).eval()==cval);
//}
//
//TEST_CASE("testing param sign functions") {
//
//
//    param<int> ip("ip");
//    ip.print();
//    ip.add_val(2);
//    ip.add_val(3);
//    ip.print();
//    CHECK(ip.is_non_negative());
//    CHECK(ip.is_positive());
//    CHECK(!ip.is_negative());
//    CHECK(!ip.is_non_positive());
//    param<> dp("dp");
//    dp.add_val(0);
//    dp.add_val(-3);
//    dp.print();
//    CHECK(!dp.is_positive());
//    CHECK(!dp.is_non_negative());
//    CHECK(!dp.is_negative());
//    CHECK(dp.is_non_positive());
//}
//
//TEST_CASE("testing variables indexing") {
//    indices ids("index_set");
//    ids = {"id1", "id2", "key3"};
//    var<> iv("iv",-2, 5);
//    iv.in(ids);
//    iv.print();
//    var<Cpx> cv("cv", Cpx(0,-1),Cpx(1,1));
//    cv.in(ids.exclude("id2"));
//    cv.print();
//    CHECK(cv.get_dim()==2);
//}
//TEST_CASE("testing vector dot product"){
//    var<> z("z",-1,1);
//    z.in(R(4));
//    param<> a("a");
//    a.set_size(3);
//    a.set_val(0, 1);
//    a.set_val(1, -1);
//    a.set_val(2, 2);
//    param<> b("b");
//    b=5;
//    CHECK(b.get_dim()==1);
//    CHECK(b.eval(0)==5);
//    b=-1;
//    b=2;
//    CHECK(b.get_dim()==3);
//    CHECK(b.eval(1)==-1);
//    b.set_val(1, 3);
//    CHECK(b.eval(1)==3);
//    z.in(R(3));
//    auto lin = (a+expo(b)).tr()*z;
//    lin.print_symbolic();
//    CHECK(lin.is_linear());
//    CHECK(lin.get_nb_instances()==1);
//    lin.print();
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
//TEST_CASE("testing complex numbers") {
//    indices ids("index_set");
//    ids = {"id1", "id2", "key3", "key4"};
//    var<> iv("x",-2, 5);
//    iv.in(ids);
//    var<Cpx> cv("y", Cpx(0,-1),Cpx(1,1));
//    constant<Cpx> cx(Cpx(-1,-2));
//    auto cx_conj = conj(cx);
//    CHECK(real(cx_conj)==real(cx));
//    CHECK(imag(cx_conj).eval()==-1*imag(cx).eval());
//    param<Cpx> px("px");
//    px = Cpx(-1,-2);
//    px.print();
//    auto px_conj = conj(px);
//    CHECK(real(px_conj)._is_real);
//    CHECK(imag(px_conj)._is_imag);
//}
//
//TEST_CASE("testing complex functions") {
//    indices ids("index_set");
//    ids = {"id1", "id2", "key3", "key4"};
//    var<> iv("x",-2, 5);
//    iv.in(ids);
//    var<Cpx> cv("y", Cpx(0,-1),Cpx(1,1));
//    cv.in(ids);
//    auto f = 2*iv;
//    f.print_symbolic();
//    CHECK(f.is_linear());
//    CHECK(f.is_convex());
//    f+= power(iv,2);
//    f.print_symbolic();
//    CHECK(f.is_quadratic());
//    CHECK(f.is_convex());
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
//}
//
//
//TEST_CASE("testing complex matrix product") {
//    var<Cpx> X("X", Cpx(0,-1),Cpx(1,1));
//    X.in(C(5,5));
//    param<Cpx> A("A");
//    A.set_size(3, 5);
//    A.set_val(0,0,Cpx(-1,1));
//    A.set_val(1,0,Cpx(2,1));
//    A.set_val(2,1,Cpx(-1,-1));
//    A.set_val(1,1,Cpx(1,1));
//    auto f = A*X;
//    f.print_symbolic();
//    CHECK(f.is_linear());
//    CHECK(f.is_convex());
//    CHECK(f.get_dim()==15);
//    var<Cpx> Y("Y", Cpx(0,-1),Cpx(1,1));
//    Y.in(C(5,5));
//    f+= X*Y;
//    f.print_symbolic();
//    CHECK(f.is_quadratic());
//    CHECK(!f.is_convex());
//    CHECK(f.is_complex());
//    CHECK(f.get_nb_vars()==50);
//}
//
//TEST_CASE("testing function convexity"){
//    var<> dp("dp",0.5,10.);
//    var<int> ip("ip",-1,1);
//    auto exp = log(dp) + sqrt(ip);
//    exp.print_symbolic();
//    CHECK(exp.is_concave());
//    var<> p("p");
//    var<> q("q");
//    auto cc = p*p + q*q;
//    cc.print_symbolic();
//    CHECK(cc.is_convex());
//    auto cc1 = cc * -1;
//    cc1.print_symbolic();
//    CHECK(cc1.is_concave());
//    cc1 += 2*p*q;
//    cc1.print_symbolic();
//    CHECK(cc1.is_rotated_soc());
//    param<int> aa("aa");
//    aa = -1;
//    aa = -3;
//    auto ff = (aa)*p*p;
//    ff.print_symbolic();
//    CHECK(ff.is_concave());
//    ff *= aa;
//    ff.print_symbolic();
//    CHECK(ff.is_convex());
//    ff *= -1;
//    ff.print_symbolic();
//    CHECK(ff.is_concave());
//    ff *= aa;
//    ff.print_symbolic();
//    CHECK(ff.is_convex());
//    ff += aa*(ip + dp)*q*q;
//    ff.print_symbolic();
//    CHECK(!ff.is_convex());
//    CHECK(!ff.is_concave());
//    CHECK(ff.is_polynomial());
//    CHECK(ff.get_nb_vars()==4);
//    param<> b("b");
//    b = 1;
//    auto fn = p*p + q*q;
//    fn.print_symbolic();
//    CHECK(fn.is_convex());
//    fn -= 2*p*q;
//    fn.print_symbolic();
//    fn += expo(p);
//    fn.print_symbolic();
//    CHECK(fn.is_convex());
//}
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
