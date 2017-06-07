//
//  Test.cpp
//  
//
//  Created by Hassan on 3 Jan 2016.
//
//


#include <stdio.h>
#include <iostream>
#include <string>
#include <stdio.h>
#include <cstring>
#include <fstream>
#include <gravity/Net.h>
#include <gravity/func.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

//  Windows
#ifdef _WIN32
#include <Windows.h>
double get_wall_time(){
    LARGE_INTEGER time,freq;
    if (!QueryPerformanceFrequency(&freq)){
        //  Handle error
        return 0;
    }
    if (!QueryPerformanceCounter(&time)){
        //  Handle error
        return 0;
    }
    return (double)time.QuadPart / freq.QuadPart;
}
double get_cpu_time(){
    FILETIME a,b,c,d;
    if (GetProcessTimes(GetCurrentProcess(),&a,&b,&c,&d) != 0){
        //  Returns total user time.
        //  Can be tweaked to include kernel times as well.
        return
        (double)(d.dwLowDateTime |
                 ((unsigned long long)d.dwHighDateTime << 32)) * 0.0000001;
    }else{
        //  Handle error
        return 0;
    }
}

//  Posix/Linux
#else
#include <time.h>
#include <sys/time.h>
double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
double get_cpu_time(){
    return (double)clock() / CLOCKS_PER_SEC;
}
#endif

int main (int argc, const char * argv[])
{

    //  Start Timers
    std::cout << "HELLO!\n";
    std::cout << "Understanding the numerical limits of your machine:" << endl;
    std::cout << "type\tlowest\thighest\n";
    std::cout << "short\t"
    << std::numeric_limits<short>::lowest() << '\t'
    << std::numeric_limits<short>::max() << '\n';
    std::cout << "unsigned\t"
    << std::numeric_limits<unsigned>::lowest() << '\t'
    << std::numeric_limits<unsigned>::max() << '\n';
    std::cout << "int\t"
    << std::numeric_limits<int>::lowest() << '\t'
    << std::numeric_limits<int>::max() << '\n';
    std::cout << "long int\t"
    << std::numeric_limits<long int>::lowest() << '\t'
    << std::numeric_limits<long int>::max() << '\n';
    std::cout << "double\t"
    << std::numeric_limits<double>::lowest() << '\t'
    << std::numeric_limits<double>::max() << '\n';
    std::cout << "long double\t"
    << std::numeric_limits<long double>::lowest() << '\t' << std::numeric_limits<long double>::max() << '\n';
//    constant<int> c(2);
//    c.print();
//    constant<float> cf(1.60);
//    cf.print();
//    param<> p;
    param<int> ip("ip");
    ip = 2;
   ip.print();
    ip = 5;
//    ip.print();
    for (int i = 0; i<100000; i++) {
        ip = 222;
    }
    param<double> dp("dp");
    dp = 1.8;
//    dp.print();
    dp = 1909092.55;
//    dp.print();
//    auto exp = log(dp) + sqrt(ip);
//    auto exp = dp*2;
//    exp.print();
    
    var<> v1("v1", 4, 10);
    v1.set_size(3);
    v1.add_bounds(-1,2);
    v1.print(true);
    var<double> p("p", 1, 50);
    p.add_bounds(0.1, 20);
    p.add_bounds(10, 60);
    var<float> q("q");
    p.set_size(200);
    q.set_size(200);
    
//  auto c1 = (p_ij^2)+(q_ij^2)-(dp^2);
//  func_ f(constant<>(2));
//  func_ f(2);
    
    p(1).print(true);
    q.print(true);
    auto f = dp*p*p;
    f.print(true);
    auto c1 = p + q - dp + 1;
    c1.print(true);
    c1 += p;
    c1.print();
    auto l2 = 2*p;
    l2.print();
    auto q1 = l2 * q;
    q1.print();
    auto l4 = -1*(dp-ip);
    l4.print();
    l4 = l4*2;
    l4.print();
    l4 += 2*dp;
    l4.print();
    l4 -= 2*ip - 2 + p;
    l4.print();
    constant<> zero = 0;
    auto l5 = l4*1;
    l5.print();
//    l3 *= -1;
//    l3.print();
    var<short> v11("v11");
    param<float> p11("p11");
//
    auto l11 = p11*v11 + ip*p - dp*q;
    l11.print();
    l11 += v11;
    l11.print();
    auto l22 = p11*v11;
    l22.print();
    l22 += 1 - p11*v11*2;
    l22.print();
    l11 += dp*q + q - ip*p + p - p11*v11*2 + v11;
    l11.print();
    auto cc = p*p + q*q;
    cc.print();
    auto cc1 = cc * -1;
    cc1.print();//SHOULD PRINT CONCAVE
    cc1 += 2*p*q;
    cc1.print();
    param<> aa("aa");
    aa = -1;
//    aa = 0;
//    auto ff = (aa*-1)*p*p - (ip + dp)*q*q;
    auto ff = (aa*-1)*p*p;
//    auto ff = -1*aa*v11*v11;
    ff.print();
    ff *= aa;
    ff.print();
    ff *= aa*aa;
    ff.print();
    ff *= aa;
    ff.print();
//    ff += aa*(ip + dp)*q*q;
//    ff.print();
    ff *= aa*p*q;
    ff.print();
    auto ppp = p*p*p;
    ppp.print();
    auto qqq = q*q*q;
    qqq.print();
    auto ss = ppp + qqq;
    ss.print();
    (ss.get_dfdx(p)).print();
    ss += 2*ip*ppp;
    ss.print();
    (ss.get_dfdx(p)).print();
    ss -= 2*ip*ppp + ppp;
    ss.print();
    (ss.get_dfdx(q)).print();
//    auto exp = log(ff);
//    exp.print();
//    l11 *= -2;
//    l11.print();
    auto l00 = 2*p(3,1) + q(1)+ p(3);
    l00.print();
    (l00.get_dfdx(p(3,1))).print();
    auto f0 = 0.1*q;
    f0.print();
    f0 -= (0.1+1e-7)*q;
    f0.print();
//    ip.print(true);
    auto vec_prod = (aa+ip).tr()*v1;
    vec_prod.print();
    vec_prod += ip.tr()*p;
    vec_prod.print();
    auto quad = (aa+ip)*(v1.tr()*v1) + q;
    quad.print();
//    int C = 10;
//    int n = 10, ni = 3;
//    var<float> alpha_ij("alpha_ij", 0, C);
//    var<float> alpha_kl("alpha_kl", 0, C);
//    param<int> y_ij("y_ij");
//    param<int> y_kl("y_kl");
//    var<>x("x"), y("y");
//    auto SVM = (3-3*x-y) + (4-3*x-y) + (4-4*x-y) + (5-4*x-y) + (y) + (y-1+x) + (y+2*x) + (y-2+3*x);
//    auto SVM = 5 - y -4*x -y + 4 -3*x - y + 3 - 3*x +4 -4*x + y + y -1 +x +y +2*x +y - 2 + 3*x;
//    SVM.print();
//    auto f1 = sum(alpha,n,ni);
//    auto f1 = sqrt(v1.tr()*v1) + p*q;    auto f1 = sqrt(v1.tr()*v1) + ip + log(p) + quad;
    
    auto f1 = sqrt(v1.tr()*v1) + ip + log(p) + quad + (p*p*p)/(q*q*dp);
    f1.print();
//    auto f2 = v11*sqrt(v1.tr()*v1) + ip + log(p) - p + expo(q) + cos(p+ip*q(1)) + sin(dp(2));
    auto f2 = sin(dp(2));
    f2.print();
//    f2 = v1/2 + sin((ip/dp)*p('i')) + 3;
    f2 = sin((ip/dp)*p);
    f2.print();
    f2 -= 2.2;
//  f2.print();
    
//  assuming a vector of coefficients are given, try to generate the graph structure.
//  1: (a, b)  1
//  2: (a, c)  1
//  3: (b, c)  1
    
    
    Net graph;
    /*
    // add nodes
    Node* node1= nullptr;
    Node* node2 = nullptr;
    Node* node3 = nullptr;

    int id = 0;
    int id2 = 1;
    int id3 = 2;
    
    node1 = new Node(to_string(id),id);
    node2 = new Node(to_string(id2),id2);
    node3 = new Node(to_string(id3),id3);
    
    Node* nodeclone1 = new Node(to_string(id),id);
    Node* nodeclone2 = new Node(to_string(id2),id2);
    Node* nodeclone3 = new Node(to_string(id3),id3);
    
    graph.add_node(node1);
    graph.add_node(node2);
    graph.add_node(node3);
    
    Arc* arc1 = NULL;
    Arc* arc2 = NULL;
    Arc* arc3 = NULL;

    string src, dest;
    arc1 = new Arc(node1, node2);
    arc2 = new Arc(node1, node3);
    arc3 = new Arc(node2, node3);
    
    
    //arc1_clone = new Arc(nodeclone1, nodeclone2);
    //arc2_clone = new Arc(nodeclone1, nodeclone3);
    //arc3_clone = new Arc(nodeclone2, nodeclone3);
    
    //Net* clone;
    graph._clone = new Net();
    graph._clone->add_node(nodeclone1);
    graph._clone->add_node(nodeclone2);
    graph._clone->add_node(nodeclone3);
    //graph.add_arc(arc1);
    //graph.add_arc(arc2);
    //graph.add_arc(arc3);
    
    //arc_clone->src = _clone->get_node(src);
    //arc_clone->dest = _clone->get_node(dest);
    
    
    graph._clone->add_arc(arc1);
    graph._clone->add_arc(arc2);
    graph._clone->add_arc(arc3);
    */
    /*
    graph._clone = new Net();

    for (unsigned i=0; i<3; i++){
        Node* node = NULL;
        Node* node_clone = NULL;
        node = new Node(to_string(i), i);
        node_clone = new Node(to_string(i),i);
        //node->vs = vs;
        graph.add_node(node);
        graph._clone->add_node(node_clone);
    }
    
    Arc* arc = NULL;
    Arc* arc_clone = NULL;
    string src, dest;
    int id;
    for (unsigned i=0; i<3; i++){
        for (unsigned j=i+1; j<3;j++){
            src = to_string(i);
            dest = to_string(j);
            id = (int)graph.arcs.size() + 1;
            arc = new Arc(to_string(id));
            arc_clone = new Arc(to_string(id));
            arc->id = id;
            arc_clone->id = id;
            arc->src = graph.get_node(src);
            arc->dest= graph.get_node(dest);
            arc_clone->src = graph._clone->get_node(src);
            arc_clone->dest = graph._clone->get_node(dest);
            graph._clone->add_arc(arc_clone);
        }
    }
*/
    
    //get the chordal extension
    //graph.get_tree_decomp_bags();
    //string filename = "/Users/guangleiwang/ANU/nesta/nesta_case30_ieee.m";

    
    
    //graph.test(filename);
    graph.test();
    
    
    // arc1->print();
    // populate the graph:
    
    
    
//    string fname = "/Users/guangleiwang/Downloads/Stable_set_instances/p.3n150.txt";
//    graph.readFile(fname);
    
    f2 *= 2;
    f2.print();
    return 0;
//    l00.print();
//    auto q00 = l00*q(1,2,3);
//    q00.print();
//    auto t00 = q(1)*q(1,2,3) + 2;
//    t00.print();
//    
//    t00 *= 0;
//    t00.print();
//    
//    q00 -= t00;
//    q00.print();
//    int i = 1;
//    int j = 2;
//    auto q11 = p(i,j)*q(i);
//    q11.print();
//    l11 *= ip;
//    l11.print();
//    c1->print();
//    c1.concretize_all(200);
    
    
//    Constraint c1;
//    c1 += log(v1[i]);
//    c1 += v1[i];
//    auto exp2 = exp*log(v1);
//    exp2.print();
//    constant_ t;
//    auto test = log(t);
//    expr<param<>> exp(ip);
//    bool test = exp.contains(c);
//    cout << test;
//    exp = 2;
//    exp = 3.455;
//    vector<constant*> constants;
    double wall0 = get_wall_time();
    double cpu0  = get_cpu_time();
//    constant a("a");
//    constant b("b");
//    constant c("c");
//    constant tmp("tmp");
//    a = 2;
//    a = 2;
//    b = 2.1;
//    b = 2;
//    c = 3;
//    c = 2;
//    tmp = 2*(a+b)/c;
//    constant tmp2("tmp2");
//    tmp2 = tmp;
//    tmp.print();
//    tmp2.print();
//    c = 20;
//    b = 1;
//    a = 9;
//    a.print();
//    b.print();
//    c.print();
//    constant tmp3("tmp3");
//    tmp3 = 2/(a+b);
//    tmp3.print();
//    constant tmp4("tmp4");
//    tmp4 = tmp3; 
//    tmp4.print();
//    var<bool> z("z");
//    z.set_id(1);
//    z = 1;
//    z = 1;
//    var<> x("x", -2.5, 10);
//    x.set_id(2);
//    x = 1;
//    x = -2;
//    x = 12;
//    var<float> y("y");
//    y.set_id(3);
//    y = 2;
//    y = 200;
//    var<int> u("u", 1,5);
//    u.set_id(4);
//    var<float> t("t", 0,10);
//    t = (float)1.2;
//    t.set_id(5);
//    z.print();
//    x.print();
//    y.print();
//    u.print();
//    t.print();
//    Function q;
////    Function q2 = c*x;
////    q=tmp3*x;
//    q = y^3;
//    q.print();
////    q += 2*a*z;
////    q.print();
////    Function dqdz = q.get_dfdx(z);
////    dqdz.print();
////    double val = dqdz.eval();
////    cout << "val = " << val << endl;
//////    return 0;
////    q += 2*a*x;
////    q.print();
////    Function dqdx = q.get_dfdx(x);
////    dqdx.print();
//////    q+=y;
////    q += 2*a*x;
////    q.print();
////    dqdx = q.get_dfdx(x);
////    dqdx.print();
////    //    q+=y;
////    
////    q -= 2*a*x;
//////    q+=u;
//////    q+=t;
////    q.print();
////    q -= a*z;
//////    q+=x;
////    q.print();
////    q -= -2*a*x+(2*a/c)*(z*y);
//    q -= (a/c)*(z*y);
//    q.print();
////    val = q.eval_dfdx(&x, 0);
////    cout << "val = " << val << endl;    
////    q += c/(2*a)*t+(2*a/c)*(z*y);
////    q.print();
////    Function dqdy = q.get_dfdx(y);
////    dqdy.print();
////    dqdz = q.get_dfdx(z);
////    dqdz.print();
////    return 0;
//    Function f;
//    f = a*x*y*z + 2*x;
//    f.print();
//    Function g;
//    g = b*t*x + a*u + a;
//    g.print();
//    
//    Function fg = f*g*x + (2*y^4);
//    fg.print();
////    fg.compute_dfdx(y);
////    Function* dfgdy = fg.get_dfdx(y);
////    dfgdy->print();
////    fg.compute_dfdx(x);
////    Function* dfgdx = fg.get_dfdx(x);
////    dfgdx->print();
//    
////    cout << "val = " << dfgdx->eval() << endl;
////    cout << "val = " << fg.eval_dfdx(&x, 0) << endl;
//
////    double v = dfgdy->eval();
////    assert(fabs(v-503.68)<10e-6);
////    assert(v ==fg.eval_dfdx(&y, 0));
////    cout << "val = " << v << endl;
//
////    return 0;
////    fg -= a*b*x*t*(y*x);
////    fg -= (a^2)*u*y*x;
////    fg.print();
////    return 0;
//    var<> wr12("wr12");
//    wr12.set_id(13);
//    var<> wr23("wr23");
//    wr23.set_id(14);
//    var<> wr31("wr31");
//    wr31.set_id(15);
//    var<> wi12("wi12");
//    wi12.set_id(16);
//    var<> wi23("wi23");
//    wi23.set_id(17);
//    var<> wi31("wi31");
//    wi31.set_id(18);
//    var<>w1("w1");
//    w1.set_id(10);
//    var<>w2("w2");
//    w2.set_id(11);
//    var<>w3("w3");
//    w3.set_id(12);
//    Function sdp = 2*(wr12*(wr23*wr31 - wi23*wi31) - wi12*(wi23*wr31 + wr23*wi31)) - (wr12*wr12 +wi12*wi12)*w3 - (wr31*wr31 + wi31*wi31)*w2 - (wr23*wr23 +wi23*wi23)*w1 + w1*w2*w3;
//    sdp.print();
////    sdp.compute_dfdx(w1);
////    Function* dfsdpdw1 = sdp.get_dfdx(w1);
////    dfsdpdw1->print();
////    sdp.compute_dfdx(wr12);
////    Function* dfsdpdwr12 = sdp.get_dfdx(wr12);
////    dfsdpdwr12->print();
//    Function q1;
//    q1 = (2/a);
//    Function q0;
//    q0 = (2/b);
//    Function res = q1 + q0;
//    res.print();
//    
//    Function func;
//    Function func0;
//    func0.print();
//    func = func0 + x;
//    func += b/a*(z^2);
//    func.print();
//    func0 += a*y;
//    func0.print();
//    func0 = 2;
//    func0 /=x;
//    func0.print();
//    func -= func0;
//    func.print();
//    func -= func0;
//    func.print();
//    func *= log(x);
//    func.print();
//    func += expo(a*y + 2/b*(t^2)) + z;
//    func.print();
//    func *= x^2;
//    func.print();
//    func /= log(y^2);
//    func.print();
//    func += x + (y^2);
//    func.print();
//    Function f1 = (2/a)*log(x);
//    f1 *= (4/c);
//    f1.print();
//    f1 += 2/a;
//    f1.print();
//    f1 = a*y*z + 2*a*b + (3/b)*x*y*t;
//    f1.print();
//    double v = f1.eval();
//    cout << "val = " << v << endl;
//    v = f1.eval_dfdx(&x, 0);
////    cout << "dfdx val = " << v << endl;
////    f1.compute_dfdx(x);
////    Function* df1dx = f1.get_dfdx(x);
////    cout << "dfdx: ";
////    df1dx->print();
////    cout << "dfdx val = " << df1dx->eval() << endl;
//    Model m("mymodel");
////    Constraint c1 = a/b*x*y*cos(t);
////    c1 <= 1;
////    m.add_constraint(c1);
//    Constraint c2("thermal_limit");
//    c2 = a*(x^2)/z + (y^2)/z - (b^2)*z;
//    c2 <= 0;
////    Constraint c2 = (2*(y^2) + (1/c)*(z^2))/(a*t);
////    c2 <= 3;
//    c2.print();
//    m.add_constraint(c2);
//    Constraint c3("thermal_limit2");
//    c3 = b*(x^2)/y + (y^2)/x - (a^2)*t;
//    c3 <= 0;
//    //    Constraint c2 = (2*(y^2) + (1/c)*(z^2))/(a*t);
//    //    c2 <= 3;
//    c3.print();
//    m.add_constraint(c3);
//    
//    m.print_functions();
////    constant t1 = -3*a;
////    constant t2 = -2*(-1*a^2);
////    constant t1t2 = t1*t2;
////    t1t2.print();
////    c2.print();
//    double cv2 = c2.eval();
//    printf("cv2 = %f\n", cv2);
//    cout << "eval c2 = " << c2.eval() << endl;
//    cout << "eval c3 = " << c3.eval() << endl;
//    m.hess_links();
//    cout << "nb nnz in Jacobian = " << m.get_nb_nnz_j() << endl;
//    cout << "nb nnz in Hessian = " << m.get_nb_nnz_h() << endl;
    return 0;

//    q-=x;
//    q.print();
//    q-= a*z;
//    q.print();
//    q += q2;
//    q.print();
//    q -= (tmp4+c)*x;
//    q.print();
//    q += 2*(a+b)*y;
//    q.print();
//    q += 2*b*u;
//    q.print();
//    q -= (2*b + 1)*u;
//    q.print();
//    q -= y;
//    q.print();
//    q -= (a+b)*y;
//    q.print();
//    q -= (a+b)*y;
//    q.print();
//    q -= 3*t;
//    q.print();
//    q += x*y;
//    q.print();
//    q += t*y;
//    q.print();
//    q -= 2*a*(x*y);
//    q.print();
//    q -= a*z;
//    q.print();
//    q += c*(t*y);
//    q.print();

    
    
    
    //  Stop timers
    double wall1 = get_wall_time();
    double cpu1  = get_cpu_time();
    cout << "\nWall clock computing time =  " << wall1 - wall0 << "\n";
    cout << "CPU computing time =  " << cpu1 - cpu0 << "\n";

    return 0;
    
}
