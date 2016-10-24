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
#include <gravity/var.h>

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
    cout << "HELLO!\n";
//    constant<int> c(2);
//    c.print();
//    constant<float> cf(1.60);
//    cf.print();
//    param<> p;
    param<short> ip("ip");
    ip = 2;
//    ip.print();
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
    auto exp = 2./3*(log(dp) + ip)/(dp^2);
//    auto exp = dp*2;
    exp.print();
    
    var<> v1("v1", 4, 10);
    v1.reserve(3);
    v1.add_bounds(-1,2);
    v1.print(true);
    
    auto exp2 = exp*log(v1);
    exp2.print();
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
