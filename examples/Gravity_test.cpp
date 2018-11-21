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
#include <gravity/Net.h>
#include <gravity/model.h>
#include <gravity/solver.h>
#include <stdio.h>
#include <stdlib.h>
#include <gravity/doctest.h>



using namespace std;
using namespace gravity;

//TEST_CASE("testing shared_ptrs") {
//    param<int> ip("ip");
//    ip.add_val(2);
//    ip.add_val(3);
//    param<> dp("dp");
//    dp.print();
//    dp.add_val(1.2);
//    dp.add_val(-1.3);
//    vector<shared_ptr<param_>> ps;
//    ps.push_back(make_shared<param<int>>(ip));
//    ps.push_back(make_shared<param<>>(dp));
//    if(ps[0]->is_integer()){
//        auto nip = (param<int>*)ps[0].get();
//        nip->print();
//        auto ndp = (param<>*)ps[1].get();
//        ndp->print();
//    }
//    
//}

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
    

//    auto cc = p*p + q*q;
////    cc.print();
//    auto cc1 = cc * -1;
//    cc1.print();//SHOULD PRINT CONCAVE
//    cc1 += 2*p*q;
////    cc1.print();
//    
//    param<int> aa("aa");
//    aa = -1;
//    aa = -3;
//    auto ff = (aa)*p*p;
////    ff.print();
//    ff *= aa;
//    ff.print();
//    ff *= -1;
//    ff.print();
////    ff.print();
//    
//    
//    ff *= aa;
////    ff.print();
////    ff += aa*(ip + dp)*q*q;
////    ff.print();
//    ff *= aa*p*q;
////    ff.print();
//    auto ppp = p*p*p;
//    ppp.print();
//    auto qqq = q*q*q;
//    qqq.print();
//    auto ss = ppp + qqq;
////    ss.print();
//    (ss.get_derivative(p)).print();
//    ss += 2*ip*ppp;
////    ss.print();
//    (ss.get_derivative(p)).print();
//    ss -= 2*ip*ppp + ppp;
////    ss.print();
//    (ss.get_derivative(q)).print();
////    auto exp = log(ff);
////    exp.print();
////    l11 *= -2;
////    l11.print();
//    auto l00 = 2*p(3,1) + q+ p(3,1);
////    l00.print();
//    (l00.get_derivative(p(3,1))).print();
//    auto f0 = 0.1*q;
////    f0.print();
//    f0 -= (0.1+1e-7)*q;
////    f0.print();
////    ip.print(true);
//    auto vec_prod = (aa+ip).tr()*v1;//fix print!
//    vec_prod.print();
//    vec_prod += ip.tr()*p;
//    vec_prod.print();
////    auto quad = (aa+ip)*(v1.tr()*v1) + q;// FIX _is_vector for v1!
//    auto quad = v1.tr()*v1;// FIX _is_vector for v1!
//    
//    quad.print(true,true);
////    int C = 10;
////    int n = 10, ni = 3;
////    var<float> alpha_ij("alpha_ij", 0, C);
////    var<float> alpha_kl("alpha_kl", 0, C);
////    param<int> y_ij("y_ij");
////    param<int> y_kl("y_kl");
////    var<>x("x"), y("y");
////    auto SVM = (3-3*x-y) + (4-3*x-y) + (4-4*x-y) + (5-4*x-y) + (y) + (y-1+x) + (y+2*x) + (y-2+3*x);
////    auto SVM = 5 - y -4*x -y + 4 -3*x - y + 3 - 3*x +4 -4*x + y + y -1 +x +y +2*x +y - 2 + 3*x;
////    SVM.print();
////    auto f1 = sum(alpha,n,ni);
////    auto f1 = sqrt(v1.tr()*v1) + p*q; auto f1 = sqrt(v1.tr()*v1) + ip + log(p) + quad;
//    
//    auto f1 = sqrt(v1.tr()*v1) + ip + log(p) + quad + (p*p*p)/(q*q*dp);
////    f1.print();
////    auto f2 = v11*sqrt(v1.tr()*v1) + ip + log(p) - p + expo(q) + cos(p+ip*q(1)) + sin(dp(2));
//    auto f2 = sin(dp);
////    f2.print();
////    f2 = v1/2 + sin((ip/dp)*p('i')) + 3;
//    f2 = sin((ip/dp)*p);
////    f2.print();
//    f2 -= 2.2;
//    f2 *= 2;
////  f2.print();
////  f2.print();
//  var<int> x("x");
//  var<int> y("y");
//  var<int> z("z");
////    
////  auto poly = -2*x - y*z + y*y*z + y*z*z + x*x*y;
////  auto poly = y*y*z + y*z*z;
////  poly.print();
////  auto dfdx = poly.get_dfdx(x);
////    
////  auto df2dx = dfdx.get_dfdx(x);
////  auto dfdy = poly.get_dfdx(y);
////  auto dfdydz = dfdy.get_dfdx(z);
////    
////  dfdy.print();
////  auto df2dy = dfdy.get_dfdx(y);
////
////    df2dx.print();
////    df2dy.print();
////    dfdydz.print();
//    Model m;
//    var<double> Xij("Xij", 0, 1);
//    var<double> Xii("Xii", 0, 1);
//    unsigned n = 2;
//    m.add_var(Xij^(n*(n-1)/2));
//    m.add_var(Xii^n);
//    
//    // new way to generate SOCP
//    ordered_pairs indices(1,n);
//    Constraint SOCP("SOCP");
//    SOCP = power(Xij.in(indices), 2) - Xii.from(indices)*Xii.to(indices);
//    m.add_constraint(SOCP <= 0);
//    constant<int> ones(1);
//    
//    constant<int> twos(2);
//    auto obj = ones.tr()*Xii + twos.tr()*Xij;
//    obj.print();
//    m.set_objective(min(obj));
//    solver s(m,ipopt);
//     s.run();
////    return 0;
////   }
////  df2dx.print();
////  df2dy.print();
////  dfdydz.print();
////  test SDPA solver.
//    
//    
//    func_ c;
//    c = 2*x;
////    test SDPA
////    auto sdpa_inst = new SdpaProgram();
////    string fname = "../data_sets/Minkcut/toy.txt";
////    std::cout << "Let me test"<< endl;
////
////    std::string file_input = "../data_sets/SDPs/example1.dat";
////    std::string file_param = "../data_sets/SDPs/param.sdpa";
////    char* f_input= new char[file_input.length()+1];
////    char* f_param= new char[file_param.length()+1];
////    strcpy(f_input,file_input.c_str());
////    strcpy(f_param,file_param.c_str());
////    sdpa_inst->read_model(f_input,f_param);
//
//    var<double> x("x");
//    var<double> y("y");
//    var<double> v1("v1");
//    var<double> v2("v2");
//    Model toy_nlp;
//    toy_nlp.add_var(x^3);
//    toy_nlp.add_var(y^3);
//    toy_nlp.add_var(v1^3);
//    toy_nlp.add_var(v2^3);
////    auto vec = constant<double>(1).tr()*power(X,2);
////    vec.print();
////
//    param<double> a("a");
//    param<double> b("b");
//    a^3;
//    b^3;
//    auto nl = cos(a*x+b*y) + sin(a*x - b*y);
//    nl.print(true);
//    auto dfdx = nl.get_derivative(x);
//    dfdx.print(true);
//    auto dfdy = nl.get_derivative(y);
//    dfdy.print(true);
//    
//    auto nl2 = v1*v2*cos(a*x+b*y) + v1*v2*sin(a*x - b*y) + expo(v1*a*x*y);
//    nl2.print(true);
//    auto dfdx2 = nl2.get_derivative(x);
//    dfdx2.print(true);
//    auto dfdy2 = nl2.get_derivative(y);
//    dfdy2.print(true);
//    auto dfdv1 = nl2.get_derivative(v1);
//    dfdv1.print(true);
//    
////    p.push_back(&b);
////    p[1]->print(true);
////    cout << "\n";
////    a = 0;
////    p[0]->print(true);
////    cout << "\n";
////    p[1]->print(true);
//    
//    func_ polynomial;
//    var<double> xvar("x", -1, 1);
//    param<double> c("c");
//    param<double> d("d");
//
//    c = 1;
//    d = -0.2;
//    x^3;
//    polynomial += xvar(1)*xvar(2)*xvar(3) + power(xvar(1), 2) + power(xvar(3),2) - 2*xvar(1)*c + c*1.2*d + c + d;
//    polynomial.print(true);
//    cout << "the constant of this polynomial function is: " << poly_eval(polynomial.get_cst()) << endl;
//    
//    unsigned n = 10;
//    param<double> ones("ones");
//    param<double> zeros("zeros");
//
//    ones.set_size(n, 1);
//    zeros.set_size(n, 0);
//    var<double> yvar("y", zeros, ones);
//    var<double> Xvar("X", -1, 1);
//
//    Xvar(1,1) =2;
//    Xvar.print(true);
//    yvar.print(true);

//    param<> p("def_p"); /* Default type is double */
//    p("bus1") = 0.3; /* Adding bus1 to the set of indices and assigning corresponding real value */
//    p("bus2") = 2.3; /* Adding bus2 to the set of indices and assigning corresponding real value */
//    p.print(true); /* Will print:  def_p = [ (bus1=0.300000) (bus2=2.300000) ]; */
//    p.set_val("bus1", 0.6); /* The function set_val can be used to set the value of the corresponding indices */
//    p.print(true);/* Will print:  def_p = [ (bus1=0.600000) (bus2=2.300000) ]; */
//    
//    param<int> pi("int_p"); /* Integer constant */
//    pi(1) = 0; /* Indexing can also be done using numerical indices */
//    pi(5) = 1;
//    pi.print(true);/* Will print:  int_p = [ (1=0) (5=1) ];*/
//    
//    param<int> lb("lb"), ub("ub");
//    lb("bus1") = -1;
//    ub("bus1") = 2;
//    var<> x("x", -0.1, 3.5);
//    var<int> y("y", lb,ub);
//    
//    x(0).print(true);
//    y("bus1").print(true);
//    
//    Constraint cstr1("cstr1");
//    cstr1 += power(x,4) - power(y, 2) + pi*x*y + 2*x;
//    cstr1 <= 2;
//    cstr1.print();
//    Constraint cstr2("cstr2");
//    cstr2 += cos(pi*x) + expo(2*y);
//    cstr2 <= 2;
//    cstr2.print();
    
//    var<> x("x"), y("y"), z("z",pos_);
//    vector<index_pair> bus_pairs;
//    bus_pairs.push_back(index_pair("bus1", "bus2"));
//    /* Second-order cone constraints */
//    Constraint SOC("SOC");
//    SOC =  power(x, 2) + power(y, 2);
//    SOC.in(bus_pairs) <= 0;
//    SOC.print();
//
//    
////    auto cc = p*p + q*q;
////    DebugOn(cc.to_str(true) << endl);
////    auto cc1 = cc * -1;
////    DebugOn(cc1.to_str(true) << endl);//SHOULD PRINT CONCAVE
////    cc1 += 2*p*q;
////    DebugOn(cc1.to_str(true) << endl);
//
//    
//    param<int> a("a");
//    a(0) = -1;
//    a(1) = -3;
//    a.print(true);
//    auto f = (a)*power(x,2);
//    f.print(true,true);
//    f *= a;
//    f.print(true,true);
//    f *= -1;
//    f.print(true,true);
//    
//    /* Declare model */
//    Model SOCP("Second-Order Cone Model");
//    SOCP.add_var(x.in(R(10))); /* Will add a vector of size 10 representing variables named x */
//    SOCP.add_constraint(SOC.in(bus_pairs) <= 0); /* Will add second-order constraints indexed by bus_pairs (see previous Code Block) */
//    SOCP.min(x+2*y); /* Declaring the objective function */
//    return 0;
   }

TEST_CASE("testing variables indexing") {
    indices ids("index_set");
    ids = {"id1", "id2", "key3"};
    var<> iv("iv",-2, 5);
    iv.in(ids);
    iv.print();
    var<Cpx> cv("cv", Cpx(0,-1),Cpx(1,1));
    cv.in(ids.exclude("id2"));
    cv.print();
    CHECK(cv.get_dim()==2);
}
TEST_CASE("testing vector dot product"){
    var<> z("z",-1,1);
    z.in(R(4));
    param<> a("a");
    a.set_size(3);
    a.set_val(0, 1);
    a.set_val(1, -1);
    a.set_val(2, 2);
    param<> b("b");
    b=5;
    CHECK(b.get_dim()==1);
    CHECK(b.eval(0)==5);
    b=-1;
    b=2;
    CHECK(b.get_dim()==3);
    CHECK(b.eval(1)==-1);
    b.set_val(1, 3);
    CHECK(b.eval(1)==3);
    z.in(R(3));
    auto lin = (a+expo(b)).tr()*z;
    lin.print_symbolic();
    CHECK(lin.is_linear());
    CHECK(lin.get_nb_instances()==1);
    lin.print();
    auto df = lin.get_dfdx(z);
    CHECK(df.is_constant());
    CHECK(df.get_dim()==3);
    df.print();
    auto lin2 = sum(z);
    lin2 += 2*z;
    lin2.print_symbolic();
    lin2.print();
    CHECK(lin2.get_nb_instances()==3);
    auto dfdz = lin2.get_dfdx(z);
    dfdz.print_symbolic();
    CHECK(dfdz==2);
    auto dfdvecz = lin2.get_dfdx(z.vec());
    dfdvecz.print_symbolic();
}

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
    f+= power(iv,2);
    f.print_symbolic();
    CHECK(f.is_quadratic());
    CHECK(f.is_convex());
    auto dfx = f.get_dfdx(iv);
    dfx.print();
    CHECK(dfx.is_linear());
    CHECK(dfx==2*iv+2);
    f *= iv;
    f.print_symbolic();
    CHECK(f.is_polynomial());
    CHECK(!f.is_complex());
    dfx = f.get_dfdx(iv);
    dfx.print();
    CHECK(dfx.is_quadratic());
    CHECK(dfx==3*power(iv,2)+4*iv);
    f += iv*log(cv);
    f.print_symbolic();
    CHECK(f.is_nonlinear());
    CHECK(f.is_complex());
    dfx = f.get_dfdx(iv);
    dfx.print();
    CHECK(dfx==3*power(iv,2)+4*iv + log(cv));
    auto dfy = f.get_dfdx(cv);
    dfy.print();
    CHECK(dfy==iv*(1/cv));
    f.in(ids.exclude("id2"));
    CHECK(f.get_nb_vars()==2);
    CHECK(f.get_nb_instances()==3);
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
    CHECK(f.is_linear());
    CHECK(f.is_convex());
    CHECK(f.get_dim()==15);
    var<Cpx> Y("Y", Cpx(0,-1),Cpx(1,1));
    Y.in(C(5,5));
    f+= X*Y;
    f.print_symbolic();
    CHECK(f.is_quadratic());
    CHECK(!f.is_convex());
    CHECK(f.is_complex());
    CHECK(f.get_nb_vars()==50);
}

TEST_CASE("testing function convexity"){
    var<> dp("dp",0.5,10.);
    var<int> ip("ip",-1,1);
    auto exp = log(dp) + sqrt(ip);
    exp.print_symbolic();
    CHECK(exp.is_concave());
    var<> p("p");
    var<> q("q");
    auto cc = p*p + q*q;
    cc.print_symbolic();
    CHECK(cc.is_convex());
    auto cc1 = cc * -1;
    cc1.print_symbolic();
    CHECK(cc1.is_concave());
    cc1 += 2*p*q;
    cc1.print_symbolic();
    CHECK(cc1.is_rotated_soc());
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
    fn += expo(p);
    fn.print_symbolic();
    CHECK(fn.is_convex());
}
