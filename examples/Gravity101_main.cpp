//
//  Test.cpp
//  
//
//  Created by Hassan on 3 Jan 2016.
//

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

using namespace std;
using namespace gravity;

int main (int argc, const char * argv[])
{

    //  Start Timers
//    std::cout << "HELLO!\n";
//    std::cout << "Understanding the numerical limits of your machine:" << endl;
//    std::cout << "type\tlowest\thighest\n";
//    
//    std::cout << "bool\t"
//    << std::numeric_limits<bool>::lowest() << '\t'
//    << std::numeric_limits<bool>::max() << '\n';
//    std::cout << "short\t"
//    << std::numeric_limits<short>::lowest() << '\t'
//    << std::numeric_limits<short>::max() << '\n';
//    std::cout << "unsigned\t"
//    << std::numeric_limits<unsigned>::lowest() << '\t'
//    << std::numeric_limits<unsigned>::max() << '\n';
//    std::cout << "int\t"
//    << std::numeric_limits<int>::lowest() << '\t'
//    << std::numeric_limits<int>::max() << '\n';
//    std::cout << "long int\t"
//    << std::numeric_limits<long int>::lowest() << '\t'
//    << std::numeric_limits<long int>::max() << '\n';
//    std::cout << "double\t"
//    << std::numeric_limits<double>::lowest() << '\t'
//    << std::numeric_limits<double>::max() << '\n';
//    std::cout << "long double\t"
//    << std::numeric_limits<long double>::lowest() << '\t' << std::numeric_limits<long double>::max() << '\n';
////    constant<int> c(2);
////    c.print();
////    constant<float> cf(1.60);
////    cf.print();
////    param<> p;
//    param<int> ip("ip");
//    std::cout << ip.to_str() << endl;
//    
//   // ip.print();
//    ip(2,3)=5;
//    
//    cout<<ip(2,3).to_str()<<endl;
//    cout<<ip(2,3).getvalue()<<endl;
//    cout <<"dim: " << ip.get_dim()<<endl;
//    ip(1,2,3)=999;
//    cout<<ip(1,2,3).to_str()<<endl;
//    cout<<ip(1,2,3).getvalue()<<endl;
//    // dim of the parameter.
//    // 0 if it is not a vector.
//    // else it is considered as a vector.
//    cout <<"dim: " << ip.get_dim()<<endl;
//
//    ip.print(true);
//    
//    
//    for (int i = 0; i<100000; i++) {
//        ip = 222;
//    }
//    param<double> dp("dp");
//    dp = 1.8;
//
////    dp.print();
//    dp = 1909092.55;
////    dp.print();
////    auto exp = log(dp) + sqrt(ip);
////    auto exp = dp*2;
////    exp.print();
//    
//    var<> v1("v1", 4, 10);
//    v1.set_size(3);
//    v1.add_bounds(-1,2);
////    v1.print(true);
//    var<double> p("p", 1, 50);
//    p.add_bounds(0.1, 20);
//    p.param<double>::print(true);
//    
//    p.add_bounds(10, 60);
//    var<float> q("q");
//    p.set_size(200);    
//    q.set_size(200);
//    
////  auto c1 = (p_ij^2)+(q_ij^2)-(dp^2);
////  func_ f(constant<>(2));
////  func_ f(2);
//    
////    p(1,2).print(true);
////    q(1,1).print(true);
//    auto f = dp*p*p;
////    f.print(true);
//    auto c1 = p + q - dp + 1;
////    c1.print(true);
//    c1 += p;
////    c1.print();
//    auto l2 = 2*p;
////    l2.print();
//    auto q1 = l2 * q;
////    q1.print();
//    auto l4 = -1*(dp-ip);
////    l4.print();
//    l4 = l4*2;
////    l4.print();
//    l4 += 2*dp;
////    l4.print();
//    l4 -= 2*ip - 2 + p;
////    l4.print();
////    constant<> zero = 0;
//    auto l5 = l4*1;
////    l5.print();
////    l3 *= -1;
////    l3.print();
//    var<short> v11("v11");
//    param<float> p11("p11");
////
//    auto l11 = p11*v11 + ip*p - dp*q;
////    l11.print();
//    l11 += v11;
////    l11.print();
//    auto l22 = p11*v11;
//    l22.print();
//    l22 += 1 - p11*v11*2;
//    l22.print();
//    l11 += dp*q + q - ip*p + p - p11*v11*2 + v11;
////    l11.print();
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
//    sdpvar<double> X("X");
//    X^5;
//    X.print();
//    auto t2 =X(3,3)+X(4, 2);
//    t2.print();
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
    
    var<> x("x"), y("y"), z("z",pos_);
    vector<index_pair> bus_pairs;
    bus_pairs.push_back(index_pair("bus1", "bus2"));
    /* Second-order cone constraints */
    Constraint SOC("SOC");
    SOC =  power(x, 2) + power(y, 2);
    SOC.in(bus_pairs) <= 0;
    SOC.print();

    
//    auto cc = p*p + q*q;
//    DebugOn(cc.to_str(true) << endl);
//    auto cc1 = cc * -1;
//    DebugOn(cc1.to_str(true) << endl);//SHOULD PRINT CONCAVE
//    cc1 += 2*p*q;
//    DebugOn(cc1.to_str(true) << endl);

    
    param<int> a("a");
    a(0) = -1;
    a(1) = -3;
    a.print(true);
    auto f = (a)*power(x,2);
    f.print(true,true);
    f *= a;
    f.print(true,true);
    f *= -1;
    f.print(true,true);
    
    /* Declare model */
    Model SOCP("Second-Order Cone Model");
    SOCP.add_var(x.in(R(10))); /* Will add a vector of size 10 representing variables named x */
    SOCP.add_constraint(SOC.in(bus_pairs) <= 0); /* Will add second-order constraints indexed by bus_pairs (see previous Code Block) */
    SOCP.min(x+2*y); /* Declaring the objective function */
    return 0;
   } 
