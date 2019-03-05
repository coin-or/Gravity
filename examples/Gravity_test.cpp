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
#include <gravity/model.h>
#include <gravity/solver.h>
#include <gravity/constraint.h>
#include <stdio.h>
#include <stdlib.h>
#include <gravity/doctest.h>
#include <PowerNet.h>
//#include <variant>

using namespace std;
using namespace gravity;



TEST_CASE("testing numerical precision") {
    double bMVA = 100.0;
    double smax = 996.0;
    double smax2 = smax/bMVA;
    CHECK(smax2==9.9600000000000008);
}

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
    lin5.print_symbolic();
    lin5.print();
    CHECK(lin5.is_complex());
    CHECK(lin5.get_dim()==5);
    CHECK(lin5.get_nb_vars()==4);
    CHECK(lin5.is_convex());    
    CHECK(lin5.to_str()=="(b)y² - ([a]ᵀ)[z] + b");
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
    CHECK(cpx_f.to_str()=="((2,0)exp(cpx))z");
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

TEST_CASE("2d Polynomial") {
    var<> x1("x1",-1,1), x2("x2",-1,1);
    Model<> M("test");
    M.add(x1.in(R(1)),x2.in(R(1)));
    M.initialize_uniform();
    M.min(pow(x1,2) - 2*pow(x1,4) + x1*x2 - 4*pow(x2,2) + 4*pow(x2,8));
    solver<> s(M,ipopt);
    s.run(5,1e-6);
    M.print();
    M.print_solution();
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

TEST_CASE("testing bounds copy"){
    var<int> x("x",-3,3);
    var<> y("y", -1.2, 1.4);
    y.copy_bounds(x);
    y.print();
    CHECK(y.get_lb()==-3);
    CHECK(y.get_ub()==3);
}

TEST_CASE("testing quadratic function factorization"){
    var<> x1("x1",0.5,10), x2("x2",-1,1);
    x1.in(R(4));x2.in(R(4));
    auto f = 3*pow(x1,2) + 5*pow(x2,2) - 2*x1*x2; /* Convex since it can be factorized to (x1 - x2)^2 + 2x1^2 + 4x2^2 */
    f.print_symbolic();
    f.print();
    CHECK(f.is_convex());
    f = 4*x1*x2 - 2*pow(x1,2) - 3*pow(x2,2); /* Concave since it can be factorized to  -(2-sqrt(2))x1^2 - (3-sqrt(2))x2^2 - (sqrt(2)x1 - sqrt(2)x2)^2*/
    f.print_symbolic();
    f.print();
    CHECK(f.is_concave());
    /* Checking convexity in model objective */
    var<> x("x",-10,10);
    auto test = make_shared<Model<>>("test");
    test->add(x.in(R(2)));
    x.initialize_uniform();
    test->min(6*x(0)*x(0) + 4*x(1)*x(1) - 2.5*x(0)*x(1));
    Constraint<> c1("c1");
    c1 = x(0)*x(1);
    test->add(c1>=8);
    test->print();
    solver<> s(test,ipopt);
    s.run(5, 1e-6);
    test->print_solution();
    CHECK(abs(test->_obj->get_val()-58.3836718)<1e-6);
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
    CHECK(cstr.to_str()==" - x1² + x2x3 + (log(a))x2 + (a)x3 + a * log(a)");
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

TEST_CASE("testing simple model"){
    Model<> M("MyModel");
    var<> x("x", 1,2), y("y", 2, 4), z("z", -1, 2);
    M.add(x.in(R(2)),y.in(R(2)),z.in(R(2)));
    Constraint<> cstr1("cstr1");
    cstr1 = pow(x,2) + pow(y,2) - pow(z,2);
    M.add(cstr1 <= 0);
    Constraint<> cstr2("cstr2");
    cstr2 = 1;
    M.add(cstr2 >= 2);//Check wrong range.
    M.max(sum(x));
    M.print();
    M.round_solution();
    CHECK(M.is_convex());
    CHECK(M.get_nb_vars()==6);
    CHECK(M.get_nb_cons()==2);
}

TEST_CASE("testing nonlinear Model"){
    Model<> M("MyModel2");
    param<> x_lb("x_lb"), x_ub("x_ub"), y_lb("y_lb"), y_ub("y_ub");
    x_lb.in(R(4));x_ub.in(R(4));y_lb.in(R(4));y_ub.in(R(4));
    x_lb = -1; x_ub = 4; y_lb = -5; y_ub = 3;
    x_lb(2) = -12;
    var<> x("x", x_lb,x_ub), y("y", y_lb, y_ub), z("z", -1, 2);
    M.add(x.in(R(4)),y.in(R(4)),z.in(R(4)));
    CHECK(x.get_lb(2)==-12);
    Constraint<> cstr1("cstr1");
    cstr1 = pow(x,2) + pow(y,2) - pow(z,2);
    M.add(cstr1 <= 0);
    param<> a("a");
    a = 2;a = 3;a = 4;a = 5;
    Constraint<> cstr2("cstr2");
    cstr2 = a*x*y*cos(x-z);
    M.add(cstr2 == 2);
    CHECK(cstr2.to_str()==" - 2 + (a)xy * cos(x - z)");
    M.max(sum(x));
    M.print_symbolic();
    M.print();
    CHECK(!M.is_convex());
    CHECK(M.get_nb_vars()==12);
    CHECK(M.get_nb_cons()==8);
}

TEST_CASE("testing acopf"){
    string fname = string(prj_dir)+"/data_sets/Power/nesta_case5_pjm.m";
    unsigned nb_threads = 2;
    double tol = 1e-6;
    string mehrotra = "no", log_level="0";
    PowerNet grid1,grid2;
    grid1.readgrid(fname);
    auto ACOPF1 = grid1.build_ACOPF(ACRECT);
    fname = string(prj_dir)+"/data_sets/Power/nesta_case14_ieee.m";
    auto ACOPF2 = grid1.build_ACOPF(ACPOL);
    auto models = {ACOPF1, ACOPF2};
    /* run in parallel */
    run_parallel(models, ipopt, tol = 1e-6, nb_threads=2);
    CHECK(abs(ACOPF1->get_obj_val()-17551.8909275818)<tol);
    CHECK(abs(ACOPF2->get_obj_val()-17551.8909275818)<tol);
    CHECK(ACOPF1->is_feasible(tol));
    ACOPF1->print_solution();
    auto Mc = ACOPF1->build_McCormick();
    auto clone = grid1.clone();
    CHECK(clone->arcs.size()==grid1.arcs.size());
    clone->remove_arc(clone->arcs.at(0));
    CHECK(clone->arcs.size()==grid1.arcs.size()-1);
}

TEST_CASE("testing socopf"){
    string fname = string(prj_dir)+"/data_sets/Power/nesta_case5_pjm.m";
//    string fname = "/Users/hlh/Dropbox/Work/Dev/pglib-opf-18.08/pglib_opf_case2383wp_k.m";
    int output = 0;
    double tol = 1e-6;
    string mehrotra = "no", log_level="0";
    PowerNet grid;
    grid.readgrid(fname);
    auto SOCOPF = grid.build_SCOPF();
//    SOCOPF->print();
    CHECK(SOCOPF->_type==quad_m);
    solver<> OPF(SOCOPF,ipopt);
    auto time_start = get_wall_time();
    OPF.run(output=5, tol=1e-6);
    auto time_end = get_wall_time();
    DebugOn("Total wall time = " << time_end - time_start << " secs" << endl);
    CHECK(abs(SOCOPF->_obj->get_val()-14999.71503774388)<tol);
    CHECK(OPF.get_nb_iterations()==24);
}

TEST_CASE("Bug in Cplex MIQCP presolve"){
    /* This test identifies a bug in Cplex's MIQCP presolve. Turn Cplex = true and relax = false to get an infeasible status, with relax=true, Cplex returns feasible. Also turning off presolve returns feasible. */
    
    /* Start indexing from 1 */
    indices R1 = range(1,1);
    indices R2 = range(1,2);
    indices R4 = range(1,4);
    indices R10 = range(1,10);
    
    Model<> m;
    
    /* Bounds */
    param<> x_lb("x_lb");
    param<> x_ub("x_ub");
    x_lb = {-2,-2,-1,35,-4,-3,0,0,-4,0};
    x_ub = {2,2,1,37,-2,5,25,4,4,4};
    
    var<> x("x",x_lb,x_ub);
    var<int> A1("A1",0,1), A2("A2",0,1), A6("A6",0,1);
    m.add(A1.in(R1), A2.in(R1), A6.in(R1));
    var<> L7("L7",0,1), L8("L8",0,1), L9("L9",0,1), L10("L10",0,1);
    m.add(x.in(R10), L7.in(R2), L8.in(R2), L9.in(R4), L10.in(R2));
    
    
    m.min(x[3]+x[4]+x[5]);
    
    Constraint<> c1("c1");
    c1 = x[3] - x[7];
    m.add(c1==0);
    Constraint<> c2("c2");
    c2 = x[4] + 14* x[1] - 3*x[8] + 14*x[2] - 6*x[9] - 3*x[10] - 19;
    m.add(c2==0);
    Constraint<> c3("c3");
    c3 = x[5] + 32* x[1] - 12*x[8] - 48*x[2] + 36*x[9] - 27*x[10] - 18;
    m.add(c3==0);
    Constraint<> c4("c4");
    c4 = x[6] - x[2] - x[1] - 1;
    m.add(c4==0);
    Constraint<> c5("c5");
    c5 = A1[1] - 1;
    m.add(c5==0);
    Constraint<> c6("c6");
    c6 = x[1] + 2*A1[1];
    m.add(c6>=0);
    Constraint<> c7("c7");
    c7 = x[1] - 2*A1[1];
    m.add(c7<=0);
    Constraint<> c8("c8");
    c8 = L8[1] + L8[2];
    m.add(c8==1);
    Constraint<> c9("c9");
    c9 = x[8] - 4*L8[1] - 4* L8[2];
    m.add(c9<=0);
    Constraint<> c10("c10");
    c10 = L8[1] - A1[1];
    m.add(c10<=0);
    Constraint<> c11("c11");
    c11 = L8[2] - A1[1];
    m.add(c11<=0);
    Constraint<> c12("c12");
    c12 = x[1] + 2* A1[1];
    m.add(c12>=0);
    Constraint<> c13("c13");
    c13 = x[1] - 2* A1[1];
    m.add(c13<=0);
    Constraint<> c14("c14");
    c14 = x[1] + 2* L8[1] - 2* L8[2];
    m.add(c14==0);
    Constraint<> c15("c15");
    c15 = A2[1] - 1;
    m.add(c15==0);
    Constraint<> c16("c16");
    c16 = x[2] + 2* A2[1];
    m.add(c16>=0);
    Constraint<> c17("c17");
    c17 = x[2] - 2* A2[1];
    m.add(c17<=0);
    Constraint<> c18("c18");
    c18 = L9[1] + L9[2] + L9[3] + L9[4] - 1;
    m.add(c18==0);
    Constraint<> c19("c19");
    c19 = x[9] - 4* L9[1] + 4* L9[2] + 4* L9[3] - 4* L9[4];
    m.add(c19==0);
    Constraint<> c20("c20");
    c20 = L9[1] + L9[3] - A2[1];
    m.add(c20<=0);
    Constraint<> c21("c21");
    c21 = L9[2] + L9[4] - A2[1];
    m.add(c21<=0);
    Constraint<> c22("c22");
    c22 = x[2] + 2* L9[1] + 2* L9[3] - 2* L9[2] - 2* L9[4];
    m.add(c22==0);
    Constraint<> c23("c23");
    c23 = L9[1] + L9[2] - A1[1];
    m.add(c23<=0);
    Constraint<> c24("c24");
    c24 = L9[3] + L9[4] - A1[1];
    m.add(c24<=0);
    Constraint<> c25("c25");
    c25 = x[1] + 2* L9[1] + 2* L9[2] - 2* L9[3] - 2* L9[4];
    m.add(c25==0);
    Constraint<> c26("c26");
    c26 = L10[1] + L10[2];
    m.add(c26==1);
    Constraint<> c27("c27");
    c27 = x[10] - 4* L10[1] - 4* L10[2];
    m.add(c27<=0);
    Constraint<> c28("c28");
    c28 = x[2] + 2* L10[1] - 2* L10[2];
    m.add(c28==0);
    Constraint<> c29("c29");
    c29 = A6[1];
    m.add(c29==1);
    Constraint<> c30("c30");
    c30 = L7[1] + L7[2];
    m.add(c30==1);
    Constraint<> c31("c31");
    c31 = x[7] - 9* L7[1] - 25* L7[2];
    m.add(c31<=0);
    Constraint<> c32("c32");
    c32 = x[6] + 3* L7[1] - 5* L7[2];
    m.add(c32==0);
    Constraint<> c33("c33");
    c33 = -1*pow(x[1],2) + x[8];
    m.add(c33>=0);
    Constraint<> c34("c34");
    c34 = -1*pow(x[2],2) + x[10];
    m.add(c34>=0);
    Constraint<> c35("c35");
    c35 = -1*pow(x[6],2) + x[7];
    m.add(c35>=0);
    m.print();
    bool use_cplex = false, relax = true;
    if(use_cplex){
        solver<> s(m,cplex);
        s.run(relax);
    }
    else {
        solver<> s(m,ipopt);
        s.run(5,1e-6);
    }
    CHECK(abs(m.get_obj_val()-31.003139013)<1e-6);
    m.print_solution();
}

TEST_CASE("Alpine issue") {
    Model<> m;
    var<> x1("x1", -2, 2), x2("x2", -2, 2);
    var<> y1("y1"), y2("y2"), y3("y3"), y4("y4");
    m.add(x1.in(R(1)), x2.in(R(1)),y1.in(R(1)), y2.in(R(1)),y3.in(R(1)),y4.in(R(1)));
    Constraint<> c1("c1");
    c1 = pow((x1 +x2 +1),2);
    CHECK(c1.is_convex());
    y1.set_lb(c1._range->first);
    y1.set_ub(c1._range->second);
    c1 -= y1;
    m.add(c1==0);
    cout << "y1 lower bound = " << y1._range->first << endl;
    cout << "y1 upper bound = " << y1._range->second << endl;
    Constraint<> c2("c2");
    //    c2 = 19 - 14*x1 + 3*pow(x1,2) - 14*x2 + 6*x1*x2 + 3*pow(x2,2);
    c2 = 19 - 14*x1 - 14*x2 + pow(sqrt(3)*x1+sqrt(3)*x2,2);
    y2.set_lb(c2._range->first);
    y2.set_ub(c2._range->second);
    cout << "y2 lower bound = " << y2._range->first << endl;
    cout << "y2 upper bound = " << y2._range->second << endl;
    c2 -= y2;
    CHECK(c2.is_convex());
    m.add(c2==0);
    Constraint<> c3("c3");
    c3 = pow(2*x1 - 3*x2,2);
    y3.set_lb(c3._range->first);
    y3.set_ub(c3._range->second);
    c3 -= y3;
    CHECK(c3.is_convex());
    m.add(c3==0);
    cout << "y3 lower bound = " << y3._range->first << endl;
    cout << "y3 upper bound = " << y3._range->second << endl;
    Constraint<> c4("c4");
    //    c4 = 18 - 32*x1 + 12*pow(x1,2) + 48*x2 - 36*x1*x2 + 27*pow(x2,2);
    c4 = 18 - 32*x1 - 6*pow(x1,2) + 48*x2 + pow(sqrt(18)*x1 - sqrt(18)*x2,2) + 9*pow(x2,2);
    y4.set_lb(c4._range->first);
    y4.set_ub(c4._range->second);
    c4 -= y4;
    CHECK(!c4.is_convex());
    m.add(c4==0);
    cout << "y4 lower bound = " << y4._range->first << endl;
    cout << "y4 upper bound = " << y4._range->second << endl;
    
    auto obj = 30 + y3*y4 +30*y1*y2 + y1*y2*y3*y4;
    cout << "Objective lower bound = " << obj._range->first << endl;
    cout << "Objective upper bound = " << obj._range->second << endl;
    m.min(obj);
    m.print_symbolic();
    //    m.initialize_uniform();
    solver<> s(m,ipopt);
    //    s.run(5,1e-6,1000000000);
    //    m.print_solution();
    Model<> relax("Lifted Relaxation");
    relax.add(x1.in(R(1)), x2.in(R(1)),y1.in(R(1)), y2.in(R(1)),y3.in(R(1)),y4.in(R(1)));
    var<> y1234("y1234"), y12("y12"), y34("y34"), x12("(x12"), x11("x11"), x22("x22");
    relax.add(x11.in(R(1)), x22.in(R(1)),x12.in(R(1)), y12.in(R(1)),y34.in(R(1)),y1234.in(R(1)));
    y1234.set_lb((y1*y2*y3*y4)._range->first);
    y1234.set_ub((y1*y2*y3*y4)._range->second);
    y12.set_lb((y1*y2)._range->first);
    y12.set_ub((y1*y2)._range->second);
    y34.set_lb((y3*y4)._range->first);
    y34.set_ub((y3*y4)._range->second);
    x12.set_lb((x1*x2)._range->first);
    x12.set_ub((x1*x2)._range->second);
    x11.set_lb((x1*x1)._range->first);
    x11.set_ub((x1*x1)._range->second);
    x22.set_lb((x1*x1)._range->first);
    x22.set_ub((x2*x2)._range->second);
    //    func<> relax_obj = 30 + y34 +30*y12 + y1234;
    func<> relax_obj = y1234;
    relax.max(relax_obj);
    relax.add(c1<=0);
    relax.add(c2<=0);
    relax.add(c3<=0);
    
    Constraint<> obj_ub("obj_ub");
    
    obj_ub = 30 + y34 +30*y12 + y1234;
    relax.add(obj_ub<=3);
    
    Constraint<> c1_relax("c1_relax");
    c1_relax = x11 + 2*x12 + x22 + 2*x1 + 2*x2 - y1 + 1;
    relax.add(c1_relax==0);
    Constraint<> c2_relax("c2_relax");
    c2_relax = 3*x11 + 6*x12 + 3*x22 - 14*x1 - 14*x2 - y2 + 19;
    relax.add(c2_relax==0);
    Constraint<> c3_relax("c3_relax");
    c3_relax = 4*x11 - 12*x12 + 9*x22 - y3;
    relax.add(c3_relax==0);
    Constraint<> c4_relax("c4_relax");
    c4_relax = 12*x11 -36*x12 + 27*x22 - 32*x1 + 48*x2 - y4 + 18;
    relax.add(c4_relax==0);
    
    Constraint<> soc1("soc1");
    soc1 = x11 - pow(x1,2);
    relax.add(soc1 >= 0);
    Constraint<> soc2("soc2");
    soc2 = x22 - pow(x2,2);
    relax.add(soc2 >= 0);
    
    Constraint<> rot_soc("rot_soc");
    rot_soc = x11*x22 - pow(x12,2);
    relax.add(rot_soc >= 0);
    CHECK(rot_soc.is_rotated_soc());
    
    relax.add_McCormick("x12", x12, x1, x2);
    relax.add_McCormick("y12", y12, y1, y2);
    relax.add_McCormick("y34", y34, y3, y4);
    relax.add_McCormick("y1234", y1234, y12, y34);
    
    relax.print_symbolic();
    solver<> s2(relax,ipopt);
    s2.run();
    relax.print_solution();
}
