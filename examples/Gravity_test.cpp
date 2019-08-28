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
#include <stdio.h>
#include <stdlib.h>
#include <gravity/solver.h>
#include <gravity/doctest.h>
#include <PowerNet.h>
//#include <variant>

using namespace std;
using namespace gravity;

#ifdef USE_MPI
#include <mpi.h>
auto err_init = MPI_Init(nullptr,nullptr);
#endif


TEST_CASE("testing constants") {
    int worker_id = 0;
#ifdef USE_MPI
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
#endif
    if(worker_id==0){
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
        CHECK(std::abs(mag0.eval()-5)<1e-8);
        auto ang0 = angle(cx2);
        ang0.println();
        CHECK(std::abs(ang0.eval()-(-2.677945045))<1e-8);
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
}

TEST_CASE("testing parameters") {
    int worker_id = 0;
#ifdef USE_MPI
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
#endif
    if(worker_id==0){
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
}

TEST_CASE("testing param copy operator") {
    int worker_id = 0;
#ifdef USE_MPI
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
#endif
    if(worker_id==0){
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
}

TEST_CASE("testing param indexing, add_val() and set_val() functions") {
    int worker_id = 0;
#ifdef USE_MPI
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
#endif
    if(worker_id==0){
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
        ids.add({"id1", "id2", "key3"});
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
}

TEST_CASE("testing matrix params") {
    int worker_id = 0;
#ifdef USE_MPI
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
#endif
    if(worker_id==0){
        param<> mat("M");
        mat.set_size(3,4);
        for (auto i = 0; i<3;i++) {
            for (auto j = 0; j<4;j++) {
                mat.set_val(i, j, 10*i+j);
            }
        }
        mat.print();
        CHECK(mat.eval(1,2)==12);
        CHECK(mat(1,2).eval()==12);
        auto tr_mat = mat.tr();
        tr_mat.print();
        CHECK(tr_mat.eval(2,1)==12);
        CHECK(tr_mat(2,1).eval()==12);

        var<> v("v",0,1);
        v.in(R(4));
        Constraint<> Mv("Mv");
        Mv = product(mat,v);
        Mv.print();

        /* Complex matrices */
        param<Cpx> Cmat("Cmat");
        Cmat.set_size(3,3);
        Cpx cval = Cpx(-1,2);
        Cmat.set_val(0, 1, cval);
        Cmat.print();
        CHECK(Cmat.eval(0,1)==cval);
        CHECK(Cmat(0,1).eval()==cval);
    }
}

TEST_CASE("testing param sign functions") {
    int worker_id = 0;
#ifdef USE_MPI
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
#endif
    if(worker_id==0){
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
}

TEST_CASE("testing variables indexing") {
    int worker_id = 0;
#ifdef USE_MPI
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
#endif
    if(worker_id==0){
        indices ids("index_set");
        ids.add({"id1", "id2", "key3"});
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
}



TEST_CASE("testing vector dot product"){
    int worker_id = 0;
#ifdef USE_MPI
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
#endif
    if(worker_id==0){
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
        auto lin6 = (2*a-exp(a)+1).tr()*z;
        lin6.print_symbolic();
        CHECK(lin6.to_str()=="([2a + 1 + -exp(a)]ᵀ)[z]");
        CHECK(lin6.is_linear());
        lin6.print();
    }
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
    int worker_id = 0;
#ifdef USE_MPI
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
#endif
    if(worker_id==0){
        indices ids("index_set");
        ids.add({"id1", "id2", "key3", "key4"});
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
}

TEST_CASE("testing complex functions") {
    int worker_id = 0;
#ifdef USE_MPI
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
#endif
    if(worker_id==0){
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
}

TEST_CASE("testing ReLU") {
    int worker_id = 0;
#ifdef USE_MPI
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
#endif
    if(worker_id==0){
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
}

TEST_CASE("2d Polynomial") {
    int worker_id = 0;
#ifdef USE_MPI
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
#endif
    if(worker_id==0){
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
}

TEST_CASE("testing range propagation") {
    int worker_id = 0;
#ifdef USE_MPI
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
#endif
    if(worker_id==0){
        indices ids("index_set");
        ids.add({"id1", "id2", "key3", "key4"});
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
}

TEST_CASE("testing polynomial functions") {
    int worker_id = 0;
#ifdef USE_MPI
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
#endif
    if(worker_id==0){
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
}

TEST_CASE("testing complex matrix product") {
    int worker_id = 0;
#ifdef USE_MPI
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
#endif
    if(worker_id==0){
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
}

TEST_CASE("testing function convexity"){
    int worker_id = 0;
#ifdef USE_MPI
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
#endif
    if(worker_id==0){
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
}

TEST_CASE("testing bounds copy"){
    int worker_id = 0;
#ifdef USE_MPI
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
#endif
    if(worker_id==0){
        var<int> x("x",-3,3);
        var<> y("y", -1.2, 1.4);
        y.copy_bounds(x);
        y.print();
        CHECK(y.get_lb(0)==-3);
        CHECK(y.get_ub(0)==3);
    }
}

TEST_CASE("testing bounds get"){
    int worker_id = 0;
#ifdef USE_MPI
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
#endif
    if(worker_id==0){
        indices ids("index_set");
        ids.add({"id1", "id2"});
        param<> lb("lb");
        lb.in(ids);
        lb.print();
        lb("id1") = 1.5;
        lb("id2") = -231.5;
        param<> ub("ub");
        ub.in(ids);
        ub.print();
        ub("id1") = 2.5;
        ub("id2") = 231.5;
        var<> x("x",lb,ub);
        x.in(ids);
        CHECK(x.get_lb("id1")==1.5);
        CHECK(x.get_lb("id2")==-231.5);
        CHECK(x.get_ub(0)==2.5);
        CHECK(x.get_ub(1)==231.5);
        Constraint<> Cst("Cst");
        Cst += x*x.get_lb();
        Cst <= x.get_ub();
        Cst.print();
        x("id1").set_lb(2);
        x("id1").set_ub(2.1);
        /** preferred version **/
        x.set_lb("id2",-1);
        x.set_ub("id2",1);
        Cst.print();/** Cst has new params */
    }
}

TEST_CASE("testing quadratic function factorization"){
    int worker_id = 0;
#ifdef USE_MPI
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
#endif
    if(worker_id==0){
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
        CHECK(std::abs(test->_obj->get_val()-58.3836718)<1e-6);
    }
}

TEST_CASE("testing constraints"){
    int worker_id = 0;
#ifdef USE_MPI
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
#endif
    if(worker_id==0){
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
}

TEST_CASE("testing soc/rotated soc constraints"){
    int worker_id = 0;
#ifdef USE_MPI
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
#endif
    if(worker_id==0){
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
}

TEST_CASE("testing nonlinear expressions"){
    int worker_id = 0;
#ifdef USE_MPI
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
#endif
    if(worker_id==0){
        var<> x1("x1", -1, 1), x2("x2", 0.1, 3), x3("x3");
        x1.in(R(1));
        x2.in(R(1));
        x3.in(R(1));
        param<> a("a"), b("b");
        a = 0.2;b = -0.6;
        param<> c("c"), d("d");
        c = 0.1;d = -0.5;
        auto cstr = x1*exp(x2*x3);
        CHECK(cstr.get_nb_vars()==3);
        CHECK(cstr.is_nonlinear());
        CHECK(cstr.get_dim()==1);
        cstr.print_symbolic();
        cstr.print();
        auto dfdx2 = cstr.get_derivative(x2);
        dfdx2.print_symbolic();
        auto f = 2*cos(min(acos(a/b), asin(c/d)))*x1;
        f.print_symbolic();
        CHECK(f.to_str()=="(2cos(min(acos(a/b), asin(c/d))))x1");
        f.print();
    }
}

TEST_CASE("testing monomials"){
    int worker_id = 0;
#ifdef USE_MPI
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
#endif
    if(worker_id==0){
        var<> x1("x1"), x2("x2"), x3("x3");
        auto cstr = x1*x2 + x2*x3 + x1*x3;
        auto monoms = cstr.get_monomials(5);
    }
}

TEST_CASE("testing monomials"){
    int worker_id = 0;
#ifdef USE_MPI
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
#endif
    if(worker_id==0){
        auto v = build_compositions(4, 2);
        for (auto i = 0; i< v.size(); i++) {
            for (auto &row:v[i]) {
                cout << row << " ";
            }
            cout << endl;
        }
    }
}

TEST_CASE("testing simple model"){
    int worker_id = 0;
#ifdef USE_MPI
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
#endif
    if(worker_id==0){
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
}

TEST_CASE("testing nonlinear Model"){
    int worker_id = 0;
#ifdef USE_MPI
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
#endif
    if(worker_id==0){
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
}

TEST_CASE("testing complex constraint expansion"){
    int worker_id = 0;
#ifdef USE_MPI
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
#endif
    if(worker_id==0){
        Model<> M("test");
        auto ids = indices("ids");
        ids.add({"id1", "id2"});
        var<> x("x",-1,1), y("y",-2,2), u1("u1",0,1), v1("v1",-3,1), u2("u2",2,3), v2("v2",0.5,1.4);
        M.add(x.in(ids),y.in(ids),u1.in(ids),v1.in(ids),u2.in(ids),v2.in(ids));
        var<Cpx> w1("w1"), w2("w2"), z("z");
        z.real_imag(x, y);
        w1.real_imag(u1, v1);
        w2.real_imag(u2, v2);

        param<> pr1("pr1"), pi1("pi1"),pr2("pr2");
        pr1 = {1,2};pi1 = {0,-1};pr2 = {-2,2};
        param<Cpx> p1("p1");
        p1.real_imag(pr1, pi1);

        Constraint<Cpx> C_lin1("C_lin1");
        C_lin1 = p1*z;
        M.add(C_lin1.in(ids)==0);

        param<Cpx> p2("p2");
        p2.set_real(pr2);/* zero imaginary */
        Constraint<Cpx> C_lin2("C_lin2");
        C_lin2 = p2*z;
        M.add(C_lin2.in(ids)==0);

        Constraint<Cpx> C_quad("C_quad");
        C_quad = z*w1;
        M.add(C_quad.in(ids)==0);

        Constraint<Cpx> C_norm("C_norm");
        C_norm = z*conj(z);
        M.add(C_norm.in(ids)==0);

        Constraint<Cpx> C_pol("C_pol");
        C_pol = z*w1*w2;
        M.add(C_pol.in(ids)==0);
        M.print();


    }
}

TEST_CASE("testing multithread solve"){
    int worker_id = 0;
#ifdef USE_MPI
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
#endif
    if(worker_id==0){
        string fname = string(prj_dir)+"/data_sets/Power/nesta_case5_pjm.m";
        unsigned nb_threads = 2;
        double tol = 1e-6;
        string mehrotra = "no", log_level="0";
        PowerNet grid1,grid2;
        grid1.readgrid(fname);
        auto ACOPF1 = build_ACOPF(grid1,ACRECT);
        fname = string(prj_dir)+"/data_sets/Power/nesta_case14_ieee.m";
        auto ACOPF2 = build_ACOPF(grid1,ACPOL);
        auto models = {ACOPF1, ACOPF2};
        /* run in parallel */
        run_parallel(models, ipopt, tol = 1e-6, nb_threads=2);
        CHECK(std::abs(ACOPF1->get_obj_val()-17551.89)<1e-3);
        CHECK(std::abs(ACOPF2->get_obj_val()-17551.89)<1e-3);
        CHECK(ACOPF1->is_feasible(tol));
        ACOPF1->print_solution();
        auto Mc = ACOPF1->build_McCormick();
        auto clone = grid1.clone();
        CHECK(clone->arcs.size()==grid1.arcs.size());
        clone->remove_arc(clone->arcs.at(0));
        CHECK(clone->arcs.size()==grid1.arcs.size()-1);
    }
}

TEST_CASE("testing socopf"){
    int worker_id = 0;
#ifdef USE_MPI
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
#endif
    if(worker_id==0){
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
        CHECK(std::abs(SOCOPF->_obj->get_val()-14999.715)<1e-3);
    }
}

TEST_CASE("Bug in Cplex MIQCP presolve"){
    int worker_id = 0;
#ifdef USE_MPI
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
#endif
    if(worker_id==0){
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
        CHECK(std::abs(m.get_obj_val()-31.003139013)<1e-6);
        m.print_solution();
    }
}

TEST_CASE("Alpine issue") {
    int worker_id = 0;
#ifdef USE_MPI
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
#endif
    if(worker_id==0){
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
}

TEST_CASE("testing absolute value function") {
    int worker_id = 0;
#ifdef USE_MPI
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
#endif
    if(worker_id==0){
        
        
        Model<> M("Test");
        var<> x("x", -6, 2);
        var<> y("y", -3, 2);
        var<> obj("obj", pos_);
        M.add(x.in(R(1)),y.in(R(1)),obj.in(R(1)));
        
        M.min(obj);
        Constraint<> obj_lb("obj_lb");
        obj_lb = obj;
        M.add(obj_lb == abs(x)*y + 6);
        
        M.print();
        solver<> NLP(M,ipopt);
        int output;
        double tol;
        string lin_solver;
        NLP.run(output=5,tol=1e-6);
        M.print_solution();
        M.print_symbolic();

    }
}


TEST_CASE("testing normal distributions") {
    int worker_id = 0;
#ifdef USE_MPI
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
#endif
    if(worker_id==0){
        
        param<> A("A");
        int n = 200, p = 50;
        A.in(R(n,p));
        double mean = 0, dev = 1;
        A.initialize_normal(mean, dev);
        
        Model<> M("Normal");
        var<> X("X", -1, 1), Xp("Xp", pos_);
        M.add(X.in(R(p)), Xp.in(R(p)));
        
        var<> s("slack");
        M.add(s.in(R(n)));
        
        M.min(sum(Xp) + 1e3*product(1,pow(s,2)));
        
        Constraint<> Abs_p("x_abs_p");
        Abs_p = X - Xp;
        M.add(Abs_p <= 0);
        
        Constraint<> Abs_n("x_abs_n");
        Abs_n = Xp + X;
        M.add(Abs_n >= 0);
        
        Constraint<> Norm2("norm2");
        Norm2 = product(1,pow(X,2));
        M.add(Norm2==1);
        
        Constraint<> Lin("lin");
        Lin = product(A,X);
        M.add(Lin==s);
        
        M.print_symbolic();
        //    M.print();
        solver<> NLP(M,ipopt);
        int output;
        double tol;
        string lin_solver;
        NLP.run(output=5,tol=1e-6);
        M.print_solution();
        
    }
}

TEST_CASE("testing set union unindexed") {
    indices ids1("index_set1");
    ids1 = indices(range(1,2), range(2,4));
    indices ids2("index_set2");
    ids2 = indices(range(1,5), range(2,4));
    auto union_set = union_ids(ids1, ids2);
    union_set.print();
    CHECK(union_set.size()==15);
    indices ids3("index_set3");
    ids3 = indices(range(1,5));
    REQUIRE_THROWS_AS(union_ids(ids1,ids3), invalid_argument);
}

TEST_CASE("testing set union indexed") {
    DebugOn("testing set union indexed" << endl);
    indices ids1("index_set1");
    ids1.add("1","2","3","4","5","6","7","8","9");
    var<>  v1("v1");
    v1.in(ids1);
    
    indices ids2("index_set2");
    ids2.add("5,4", "3,5", "6,7", "7,8", "2,6", "8,9");
    var<>  v2("v2");
    v2 = v1.from(ids2);
    
    var<>  v3("v3");
    v3 = v1.to(ids2);
    
    auto union_set = union_ids(*v2._indices, *v3._indices); //in this case, the keys are same, so the union should check not only keys but also ._ids in the individual index sets and add accordingly
    CHECK(union_set.size()==8); //in the _ids, the function should work in a way that union_set._ids = [4,2,5,6,1,7,3,8]
    REQUIRE_THROWS_AS(union_ids(ids1,ids2), invalid_argument);
}

TEST_CASE("testing from_ith() function") {
    indices ids("index_set");
    ids = indices(range(1,3),range(9,10), range(2,4));
    param<> dp("dp");
    dp.in(range(2,4));
    dp.print();
    dp("2") = 1.5;
    dp("4") = -231.5;
    dp.print();
    CHECK(dp.eval("2")==1.5);
    CHECK(dp.eval("4")==-231.5);
    CHECK(dp._range->first==-231.5);
    CHECK(dp._range->second==1.5);
    REQUIRE_THROWS_AS(dp("unexisting_key").eval(), invalid_argument);
    ids.print();
    auto ndp = dp.from_ith(2,ids);
    ndp.print();
    CHECK(ndp.get_dim()==ids.size());
    indices ids2("index_set2");
    ids2 = indices(range(1,3),range(9,10), range(2,4));
    var<> dv("dv");
    dv.in(range(9,10));
    int precision = 5;
    dv.print_vals(precision=5);
    auto ndv = dv.from_ith(1,ids2);
    ndv.print_vals(precision=5);
    var<> dv2("dv2");
    dv2.in(range(1,3));
    auto ndv2 = dv2.from_ith(0,ids2);
    ndv2.print_vals(precision=5);
}


TEST_CASE("testing in_ignore_ith() function") {
    indices ids("index_set");
    ids = indices(range(1,3),range(9,10), range(2,4));
    param<> dp("dp");
    dp.in(range(1,3),range(2,4));
    dp("1,2") = 1.5;
    dp("3,4") = -231.5;
    dp.print();
    auto ndp = dp.in_ignore_ith(1, 1, ids);
    ndp.print();
    CHECK(ndp.get_dim()==ids.size());
    indices idsv1("index_setv1");
    idsv1.add("id1", "id11", "id111");
    indices idsv2("index_setv2");
    idsv2.add("id2", "id22", "id222");
    var<> dv("dv");
    dv.in(idsv1,idsv2);
    int precision = 5;
    dv.print_vals(precision=5);
    indices ids_all("index_set_all");
    ids_all.add("id1,id3,id2", "id11,id33,id22", "id111,id333,id222");
    auto ndv = dv.in_ignore_ith(1, 1, ids_all);
    ndv.print_vals(precision=5);
}

TEST_CASE("testing get_matrix()") {
    auto ids = indices(range(1,3),range(8,12));
    var<> dv("dv");
    dv = dv.in(ids);
    dv.print_vals(4);
    
    Constraint<> Sum1("Sum1");
    Sum1 = sum(dv.in_matrix(0,1)); // when we say this it should start from entry 0 and get the first indices range(1,3) in the second dimension
    Sum1.print();
    CHECK(Sum1.get_nb_instances() == 5); // so when we actually sum, the number of instances should be equal to range(8,12) = 5
    
    Constraint<> Sum2("Sum2");
    Sum2 = sum(dv.in_matrix(1,1));
    Sum2.print();
    CHECK(Sum2.get_nb_instances() == 3);
    
    auto ids1 = indices(range(1,3),range(4,7),range(8,12));
    var<> dv1("dv1");
    dv1 = dv1.in(ids1);
    dv1.print_vals(4);
    
    auto dv4 = dv1.in_matrix(1,1); //this case is wrong as well, it was working before because all the entries had only 3 elements.
    Constraint<> Sum3("Sum3");
    Sum3 = sum(dv4);
    Sum3.print();
    CHECK(Sum3.get_nb_instances() == 15);
}

TEST_CASE("testing sum_ith()") {
    indices ids1("index_set2");
    ids1.add("5,4,1", "5,2,1", "7,8,4", "5,5,1", "7,6,4", "7,9,4");
    var<>  v1("v1");
    v1.in(ids1);
    Constraint<> Sum1("Sum1");
    Sum1 = sum_ith(v1, 1, 1);
    stringstream buffer;
    auto console = cout.rdbuf(buffer.rdbuf());
    Sum1.print();
    cout.rdbuf(console);
    CHECK(buffer.str()==" Sum1 (Linear) : \nSum1[0]: v1[5,4,1] + v1[5,2,1] + v1[5,5,1] <= 0;\nSum1[1]: v1[7,8,4] + v1[7,6,4] + v1[7,9,4] <= 0;\n");
    CHECK(Sum1.get_nb_instances() == 2);
}

TEST_CASE("sum over outgoing") {
    int precision = 4;
    /** Define indices */
    indices arcs("arcs");
    arcs.add("a1,1,2", "a2,1,3", "a3,1,4", "a4,3,4", "a5,2,4", "a6,1,5", "a7,1,5");
    indices nodes("nodes");
    nodes.add("1", "2", "3", "4", "5");
    /** Declare model,vars and constraints */
    Model<> M("Test");
    var<>  v1("flux", 0, 1);
    M.add(v1.in(arcs));
    Constraint<> Sum0("Sum0");
    Sum0 = v1.sum_out(nodes) + v1.sum_in(nodes);
    M.add(Sum0.in(nodes) == 0);
    M.print_symbolic();
    M.print();
    CHECK(Sum0.get_nb_instances() == nodes.size());
}

TEST_CASE("testing sum_ith() func<> version"){
    indices ids1("index set1""");
    ids1 = indices(range(1,3),range(1,4),range(1,6));
    
    indices ids2("index set2");
    ids2 = indices(range(1,3),range(1,4),range(1,5),range(1,6));
    
    param<> p1("p1");
    p1.in(ids1);
    size_t pos = 0;
    for(auto i = 1; i<= 3; i++){
        for(auto j = 1; j<= 4; j++){
            for(auto k = 1; k<= 6; k++){
                p1.set_val(pos, 100*i+10*j+k);
                pos++;
            }
        }
    }
    
    var<> v2("v2");
    v2.in(ids2);
    auto pp1 = p1.in_ignore_ith(2,1,ids2);
    pp1.print_vals(4);
    Constraint<> Sum0("Sum0");//p1.in_ignore_ith(2,1,ids2)*
//    Sum0 = sum_ith(pp1,1,1);
//    Sum0 = sum(pp1);
//    Sum0.print();
    
    Constraint<> Sum1("Sum1");//p1.in_ignore_ith(2,1,ids2)*
    Sum1 = sum_ith(v2,1,2); //supposed to ignore range(1,5) and then sum over range(1,4) and range(1,5) assuming the function has 4 entries in the index set (I believe it should)
    // I am also assuming that this function works for param<> as well
//    Sum1 = pp1.tr()*v2.in_matrix(1,2);
    Sum1.print();
    CHECK(Sum1.get_nb_instances() == 3*6);
    
}

TEST_CASE("testing Outer Approximation") {
    indices buses("buses");
    buses.insert("1", "2", "3", "4");
    indices bus_pairs("bpairs");
    bus_pairs.insert("1,2", "1,3", "3,4", "4,1");
    
    Model<> Mtest("Mtest");
    var<>  R_Wij("R_Wij", -1, 1);
    /* Imaginary part of Wij = ViVj */
    var<>  Im_Wij("Im_Wij", -1, 1);
    var<>  Wii("Wii", 0.8, 1.21);
    Mtest.add(R_Wij.in(bus_pairs), Im_Wij.in(bus_pairs), Wii.in(buses));
    Constraint<> SOC("SOC");
    SOC = pow(R_Wij, 2) + pow(Im_Wij, 2) - Wii.from(bus_pairs);
    SOC.in(bus_pairs);
    SOC.print();
    R_Wij.gravity::param<double>::set_val("1,2",0.5);
    R_Wij.gravity::param<double>::set_val("3,4",-0.5);
    Im_Wij.gravity::param<double>::set_val("1,2",-0.5);
    Im_Wij.gravity::param<double>::set_val("3,4",-0.5);
    Wii.gravity::param<double>::set_val("1",0.5);
    Wii.gravity::param<double>::set_val("3",0.5);
    SOC.eval_all();
    SOC.compute_derivatives();
    Constraint<> OA ("OA");
    OA = SOC.get_outer_app();
    Mtest.add(OA<=0);
    OA.print();
    Mtest.print();
}

#ifdef USE_MPI
TEST_CASE("testing OpenMPI") {
    DebugOn("testing OpenMPI" << endl);
    int worker_id;
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
    unsigned nb_threads = 2;
    double tol = 1e-6;
    PowerNet grid1,grid2;
    string fname = string(prj_dir)+"/data_sets/Power/nesta_case5_pjm.m";
    grid1.readgrid(fname);
    fname = string(prj_dir)+"/data_sets/Power/nesta_case14_ieee.m";
    grid2.readgrid(fname);
    auto ACOPF1 = build_ACOPF(grid1,ACRECT);
    auto SOCOPF1 = build_SDPOPF(grid1);
    auto ACOPF2 = build_ACOPF(grid2,ACPOL);
    auto SOCOPF2 = build_SDPOPF(grid2);
    auto models = {ACOPF1, SOCOPF1, ACOPF2, SOCOPF2};
    /* run in parallel */
    run_MPI(models, ipopt, tol=1e-6, nb_threads=2);
    if(worker_id==0){
    	//CHECK(std::abs(ACOPF1->get_obj_val()-17551.890927)<tol);
    	//CHECK(ACOPF1->is_feasible(tol));
    	ACOPF1->print_solution();
    	ACOPF2->print_solution();
    	SOCOPF2->print_solution();
    }
    MPI_Finalize();
}
#endif

