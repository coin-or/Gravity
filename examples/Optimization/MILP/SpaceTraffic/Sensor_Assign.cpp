//
//  main.cpp
//  bilevel_sensor
//
//  Created by Svetlana Riabova on 6/8/22.
//

#include <iostream>
#include "Sensor_Assign.hpp"
#include <gravity/rapidcsv.h>
int main(int argc, const char * argv[]) {
    
    int time_horizon = 3600;
    int nb_steps = 360;
    double scale = 1;
    double dt = time_horizon/nb_steps;
    double mu = 398600.435507; // GM in km^3/s^2
    double min_dist = 10;
    
    /* Model */
    Model<> Mopt("Mopt");
    
    /* Index sets */
    indices T("T"), T_excl_first("T_excl_first"), T_excl_last("T_excl_last");
    T = range(1,nb_steps);
    T_excl_first = range(2,nb_steps);
    T_excl_last = range(1,nb_steps-1);
    
    /* Input parameters */
    param<> x_p("x_p"),y_p("y_p"),z_p("z_p"),vx_p("vx_p"),vy_p("vy_p"),vz_p("vz_p");/*< Primary coorindates and velocity */
    x_p.in(T);y_p.in(T);z_p.in(T);
    vx_p.in(T);vy_p.in(T);vz_p.in(T);
    string fname = string(prj_dir)+"/data_sets/SpaceCraft/primary_motion.csv";
    rapidcsv::Document input_data(fname, rapidcsv::LabelParams(-1, -1),
                                  rapidcsv::SeparatorParams(','));
    auto nb_rows = input_data.GetRowCount();
    if(nb_rows<time_horizon){
        throw invalid_argument("Input file has insufficient data points, it requires " + to_string(nb_steps) + " rows\n");
    }
    for (int i = 0; i< nb_steps; i++) {
        x_p.set_val(i,input_data.GetCell<double>(1, i*dt));
        y_p.set_val(i,input_data.GetCell<double>(2, i*dt));
        z_p.set_val(i,input_data.GetCell<double>(3, i*dt));
        vx_p.set_val(i,input_data.GetCell<double>(4, i*dt));
        vy_p.set_val(i,input_data.GetCell<double>(5, i*dt));
        vz_p.set_val(i,input_data.GetCell<double>(6, i*dt));
    }
    DebugOn("Done reading input file\n");
    
    
    /* Variables */
    var<> x("x",-1e4, 1e4), y("y", -1e4, 1e4), z("z", 1e3, 1e4);/*< (x,y,z) object coordinates */
    var<> vx("vx", -1e2, 1e2), vy("vy", -1e2, 1e2), vz("vz", -1e2, 1e2);/*< (vx,vy,vz) object velocity */
    var<> ux("ux", -1e4, 1e4), uy("uy", -1e4, 1e4), uz("uz", -1e4, 1e4);/*< (ux,uy,uz) thrust applied at given coordinates */
    var<int> b("b", 0, 1);/*< if b_t = 1, thrust is applied at time step t */
//    var<> ux("ux", 0, 0), uy("uy", 0, 0), uz("uz", -1e2, 1e2);/*< (ux,uy,uz) thrust applied at given coordinates */
    /* Runge-Kutta auxiliary variables */
    var<> k1_x("k1_x"), k2_x("k2_x"), k3_x("k3_x"), k4_x("k4_x");
    var<> k1_y("k1_y"), k2_y("k2_y"), k3_y("k3_y"), k4_y("k4_y");
    var<> k1_z("k1_z"), k2_z("k2_z"), k3_z("k3_z"), k4_z("k4_z");
    var<> k1_vx("k1_vx"), k2_vx("k2_vx"), k3_vx("k3_vx"), k4_vx("k4_vx");
    var<> k1_vy("k1_vy"), k2_vy("k2_vy"), k3_vy("k3_vy"), k4_vy("k4_vy");
    var<> k1_vz("k1_vz"), k2_vz("k2_vz"), k3_vz("k3_vz"), k4_vz("k4_vz");
    var<> r1("r1",5e3,sqrt(3)*1e4), r2("r2",5e3,sqrt(3)*1e4), r3("r3",5e3,sqrt(3)*1e4), r4("r4",5e3,sqrt(3)*1e4);/*< r = sqrt(x^2 + y^2 + z^2) */
    var<> rc1("rc1",1./std::pow((sqrt(3)*1e4),3), 1./std::pow(5e3,3)), rc2("rc2",1./std::pow((sqrt(3)*1e4),3), 1./std::pow(5e3,3)), rc3("rc3",1./std::pow((sqrt(3)*1e4),3), 1./std::pow(5e3,3)), rc4("rc4",1./std::pow((sqrt(3)*1e4),3), 1./std::pow(5e3,3));/*< rci = ri^3*/
    
    Mopt.add(x.in(T), y.in(T), z.in(T), vx.in(T), vy.in(T), vz.in(T), ux.in(T), uy.in(T), uz.in(T), k1_x.in(T), k2_x.in(T), k3_x.in(T), k4_x.in(T), k1_y.in(T), k2_y.in(T), k3_y.in(T), k4_y.in(T), k1_z.in(T), k2_z.in(T), k3_z.in(T), k4_z.in(T), k1_vx.in(T), k2_vx.in(T), k3_vx.in(T), k4_vx.in(T), k1_vy.in(T), k2_vy.in(T), k3_vy.in(T), k4_vy.in(T), k1_vz.in(T), k2_vz.in(T), k3_vz.in(T), k4_vz.in(T), r1.in(T), r2.in(T), r3.in(T), r4.in(T), rc1.in(T), rc2.in(T), rc3.in(T), rc4.in(T));
    
    bool add_impulse_max = false;
    if(add_impulse_max){
        
        Mopt.add(b.in(T));
        
        Constraint<> thrust_on_off_ux("thrust_on_off_ux");
        thrust_on_off_ux = ux - 1e4*b;
        Mopt.add(thrust_on_off_ux.in(T) <= 0);
        
        Constraint<> thrust_on_off_uy("thrust_on_off_uy");
        thrust_on_off_uy = uy - 1e4*b;
        Mopt.add(thrust_on_off_uy.in(T) <= 0);
        
        Constraint<> thrust_on_off_uz("thrust_on_off_uz");
        thrust_on_off_uz = uz - 1e4*b;
        Mopt.add(thrust_on_off_uz.in(T) <= 0);
        
        Constraint<> thrust_on_off_ux_neg("thrust_on_off_ux_neg");
        thrust_on_off_ux_neg = ux + 1e4*b;
        Mopt.add(thrust_on_off_ux_neg.in(T) >= 0);
        
        Constraint<> thrust_on_off_uy_neg("thrust_on_off_uy_neg");
        thrust_on_off_uy_neg = uy + 1e4*b;
        Mopt.add(thrust_on_off_uy_neg.in(T) >= 0);
        
        Constraint<> thrust_on_off_uz_neg("thrust_on_off_uz_neg");
        thrust_on_off_uz_neg = uz + 1e4*b;
        Mopt.add(thrust_on_off_uz_neg.in(T) >= 0);
        
        Constraint<> max_impulses("max_impulses");
        max_impulses = sum(b) - 2;
        Mopt.add(max_impulses <= 0);
    }
    
    
    /* Objective */
    Mopt.min(sum(ux*ux) + sum(uy*uy) + sum(uz*uz));
    
    /* Constraints */
    
    Constraint<> min_dist_def("min_dist_def");
    min_dist_def = abs(x-x_p) + abs(y-y_p) + abs(z-z_p);
    Mopt.add(min_dist_def.in(T) >= min_dist);
//
//    Constraint<> min_dx_lb("min_dx_lb");
//    min_dx_lb = abs(x - x_p);
//    Mopt.add(min_dx_lb.in(T) >= min_dist);

//    Constraint<> min_dx_ub("min_dx_ub");
//    min_dx_ub = (x_p - x);
//    Mopt.add(min_dx_ub.in(T) >= min_dist);

//    Constraint<> min_dy_lb("min_dy_lb");
//    min_dy_lb = abs(y - y_p);
//    Mopt.add(min_dy_lb.in(T) >= min_dist);
//
//    Constraint<> min_dy_ub("min_dy_ub");
//    min_dy_ub = (y_p - y);
//    Mopt.add(min_dy_ub.in(T) >= min_dist);
//
//    Constraint<> min_dz_lb("min_dz_lb");
//    min_dz_lb = abs(z - z_p);
//    Mopt.add(min_dz_lb.in(T) >= min_dist);
//
//    Constraint<> min_dz_ub("min_dz_ub");
//    min_dz_ub = (z_p - z);
//    Mopt.add(min_dz_ub.in(T) >= min_dist);
    
    Constraint<> linking_x("linking_x");
    linking_x = x.in(T_excl_first) - (x.in(T_excl_last) + dt/6*(k1_x.in(T_excl_last) + 2*k2_x.in(T_excl_last) + 2*k3_x.in(T_excl_last) + k4_x.in(T_excl_last)));
    Mopt.add(linking_x.in(T_excl_first)==0);
    
    Constraint<> linking_y("linking_y");
    linking_y = y.in(T_excl_first) - (y.in(T_excl_last) + dt/6*(k1_y.in(T_excl_last) + 2*k2_y.in(T_excl_last) + 2*k3_y.in(T_excl_last) + k4_y.in(T_excl_last)));
    Mopt.add(linking_y.in(T_excl_first)==0);
    
    Constraint<> linking_z("linking_z");
    linking_z = z.in(T_excl_first) - (z.in(T_excl_last) + dt/6*(k1_z.in(T_excl_last) + 2*k2_z.in(T_excl_last) + 2*k3_z.in(T_excl_last) + k4_z.in(T_excl_last)));
    Mopt.add(linking_z.in(T_excl_first)==0);
    
    Constraint<> linking_vx("linking_vx");
    linking_vx = vx.in(T_excl_first) - (vx.in(T_excl_last) + dt/6*(k1_vx.in(T_excl_last) + 2*k2_vx.in(T_excl_last) + 2*k3_vx.in(T_excl_last) + k4_vx.in(T_excl_last)));
    Mopt.add(linking_vx.in(T_excl_first)==0);
    
    Constraint<> linking_vy("linking_vy");
    linking_vy = vy.in(T_excl_first) - (vy.in(T_excl_last) + dt/6*(k1_vy.in(T_excl_last) + 2*k2_vy.in(T_excl_last) + 2*k3_vy.in(T_excl_last) + k4_vy.in(T_excl_last)));
    Mopt.add(linking_vy.in(T_excl_first)==0);
    
    Constraint<> linking_vz("linking_vz");
    linking_vz = vz.in(T_excl_first) - (vz.in(T_excl_last) + dt/6*(k1_vz.in(T_excl_last) + 2*k2_vz.in(T_excl_last) + 2*k3_vz.in(T_excl_last) + k4_vz.in(T_excl_last)));
    Mopt.add(linking_vz.in(T_excl_first)==0);
    
    Constraint<> k1_x_def("k1_x_def");
    k1_x_def = k1_x - vx;
    Mopt.add(k1_x_def.in(T)==0);

    Constraint<> k1_y_def("k1_y_def");
    k1_y_def = k1_y - vy;
    Mopt.add(k1_y_def.in(T)==0);

    Constraint<> k1_z_def("k1_z_def");
    k1_z_def = k1_z - vz;
    Mopt.add(k1_z_def.in(T)==0);

    Constraint<> k1_vx_def("k1_vx_def");
    k1_vx_def = k1_vx - (ux - mu*x*rc1);
    Mopt.add(k1_vx_def.in(T)==0);

    Constraint<> k1_vy_def("k1_vy_def");
    k1_vy_def = k1_vy - (uy - mu*y*rc1);
    Mopt.add(k1_vy_def.in(T)==0);

    Constraint<> k1_vz_def("k1_vz_def");
    k1_vz_def = k1_vz - (uz - mu*z*rc1);
    Mopt.add(k1_vz_def.in(T)==0);
    
    
    Constraint<> k2_x_def("k2_x_def");
    k2_x_def = k2_x - (vx + k1_vx*dt/2.);
    Mopt.add(k2_x_def.in(T)==0);

    Constraint<> k2_y_def("k2_y_def");
    k2_y_def = k2_y - (vy + k1_vy*dt/2.);
    Mopt.add(k2_y_def.in(T)==0);

    Constraint<> k2_z_def("k2_z_def");
    k2_z_def = k2_z - (vz + k1_vz*dt/2.);
    Mopt.add(k2_z_def.in(T)==0);

    Constraint<> k2_vx_def("k2_vx_def");
    k2_vx_def = k2_vx.in(T_excl_last) - ((ux.in(T_excl_last) + ux.in(T_excl_first))/2. - mu*rc2.in(T_excl_last)*(x.in(T_excl_last) + k1_x.in(T_excl_last)*dt/2.));
    Mopt.add(k2_vx_def.in(T_excl_last)==0);

    Constraint<> k2_vy_def("k2_vy_def");
    k2_vy_def = k2_vy.in(T_excl_last) - ((uy.in(T_excl_last) + uy.in(T_excl_first))/2. - mu*rc2.in(T_excl_last)*(y.in(T_excl_last) + k1_y.in(T_excl_last)*dt/2.));
    Mopt.add(k2_vy_def.in(T_excl_last)==0);

    Constraint<> k2_vz_def("k2_vz_def");
    k2_vz_def = k2_vz.in(T_excl_last) - ((uz.in(T_excl_last) + uz.in(T_excl_first))/2. - mu*rc2.in(T_excl_last)*(z.in(T_excl_last) + k1_z.in(T_excl_last)*dt/2.));
    Mopt.add(k2_vz_def.in(T_excl_last)==0);
    
    Constraint<> k3_x_def("k3_x_def");
    k3_x_def = k3_x - (vx + k2_vx*dt/2.);
    Mopt.add(k3_x_def.in(T)==0);

    Constraint<> k3_y_def("k3_y_def");
    k3_y_def = k3_y - (vy + k2_vy*dt/2.);
    Mopt.add(k3_y_def.in(T)==0);

    Constraint<> k3_z_def("k3_z_def");
    k3_z_def = k3_z - (vz + k2_vz*dt/2.);
    Mopt.add(k3_z_def.in(T)==0);

    Constraint<> k3_vx_def("k3_vx_def");
    k3_vx_def = k3_vx.in(T_excl_last) - ((ux.in(T_excl_last) + ux.in(T_excl_first))/2. - mu*rc3.in(T_excl_last)*(x.in(T_excl_last) + k2_x.in(T_excl_last)*dt/2.));
    Mopt.add(k3_vx_def.in(T_excl_last)==0);

    Constraint<> k3_vy_def("k3_vy_def");
    k3_vy_def = k3_vy.in(T_excl_last) - ((uy.in(T_excl_last) + uy.in(T_excl_first))/2. - mu*rc3.in(T_excl_last)*(y.in(T_excl_last) + k2_y.in(T_excl_last)*dt/2.));
    Mopt.add(k3_vy_def.in(T_excl_last)==0);

    Constraint<> k3_vz_def("k3_vz_def");
    k3_vz_def = k3_vz.in(T_excl_last) - ((uz.in(T_excl_last) + uz.in(T_excl_first))/2. - mu*rc3.in(T_excl_last)*(z.in(T_excl_last) + k2_z.in(T_excl_last)*dt/2.));
    Mopt.add(k3_vz_def.in(T_excl_last)==0);
    
    
    Constraint<> k4_x_def("k4_x_def");
    k4_x_def = k4_x - (vx + k3_vx*dt);
    Mopt.add(k4_x_def.in(T)==0);

    Constraint<> k4_y_def("k4_y_def");
    k4_y_def = k4_y - (vy + k3_vy*dt);
    Mopt.add(k4_y_def.in(T)==0);

    Constraint<> k4_z_def("k4_z_def");
    k4_z_def = k4_z - (vz + k3_vz*dt);
    Mopt.add(k4_z_def.in(T)==0);

    Constraint<> k4_vx_def("k4_vx_def");
    k4_vx_def = k4_vx.in(T_excl_last) - (ux.in(T_excl_first) - mu*rc4.in(T_excl_last)*(x.in(T_excl_last) + k3_x.in(T_excl_last)*dt));
    Mopt.add(k4_vx_def.in(T_excl_last)==0);

    Constraint<> k4_vy_def("k4_vy_def");
    k4_vy_def = k4_vy.in(T_excl_last) - (uy.in(T_excl_first) - mu*rc4.in(T_excl_last)*(y.in(T_excl_last) + k3_y.in(T_excl_last)*dt));
    Mopt.add(k4_vy_def.in(T_excl_last)==0);

    Constraint<> k4_vz_def("k4_vz_def");
    k4_vz_def = k4_vz.in(T_excl_last) - (uz.in(T_excl_first) - mu*rc4.in(T_excl_last)*(z.in(T_excl_last) + k3_z.in(T_excl_last)*dt));
    Mopt.add(k4_vz_def.in(T_excl_last)==0);

    Constraint<> r1_def("r1_def");
    r1_def = r1 - sqrt(x*x + y*y + z*z);
    Mopt.add(r1_def.in(T)==0);

    Constraint<> r2_def("r2_def");
    r2_def = r2 - sqrt(pow(x+k1_x*dt/2.,2) + pow(y+k1_y*dt/2.,2) + pow(z+k1_z*dt/2.,2));
    Mopt.add(r2_def.in(T)==0);

    Constraint<> r3_def("r3_def");
    r3_def = r3 - sqrt(pow(x+k2_x*dt/2.,2) + pow(y+k2_y*dt/2.,2) + pow(z+k2_z*dt/2.,2));
    Mopt.add(r3_def.in(T)==0);

    Constraint<> r4_def("r4_def");
    r4_def = r4 - sqrt(pow(x+k3_x*dt/2.,2) + pow(y+k3_y*dt/2.,2) + pow(z+k3_z*dt/2.,2));
    Mopt.add(r4_def.in(T)==0);
    
    Constraint<> rc1_def("rc1_def");
    rc1_def = rc1*pow(r1,3);
    Mopt.add(rc1_def.in(T)==1);

    Constraint<> rc2_def("rc2_def");
    rc2_def = rc2*pow(r2,3);
    Mopt.add(rc2_def.in(T)==1);

    Constraint<> rc3_def("rc3_def");
    rc3_def = rc3*pow(r3,3);
    Mopt.add(rc3_def.in(T)==1);

    Constraint<> rc4_def("rc4_def");
    rc4_def = rc4*pow(r4,3);
    Mopt.add(rc4_def.in(T)==1);
    
    /* Initial Conditions */
    
    indices T_start("T_start");
    T_start.add(*T._keys->begin());
    
    
    Constraint<> x_init("x_init");
    x_init = x.in(T_start) - (x_p.in(T_start) + min_dist);
    Mopt.add(x_init.in(T_start)==0);

    Constraint<> y_init("y_init");
    y_init = y.in(T_start) - (y_p.in(T_start) + min_dist);
    Mopt.add(y_init.in(T_start)==0);

    Constraint<> z_init("z_init");
    z_init = z.in(T_start) - (z_p.in(T_start) + min_dist);
    Mopt.add(z_init.in(T_start)==0);

    Constraint<> vx_init("vx_init");
    vx_init = vx.in(T_start) - vx_p.in(T_start);
    Mopt.add(vx_init.in(T_start)==0);

    Constraint<> vy_init("vy_init");
    vy_init = vy.in(T_start) - vy_p.in(T_start);
    Mopt.add(vy_init.in(T_start)==0);

    Constraint<> vz_init("vz_init");
    vz_init = vz.in(T_start) - vz_p.in(T_start);
    Mopt.add(vz_init.in(T_start)==0);

//    Mopt.initialize_midpoint();
    x.initialize_all(10);
    y.initialize_all(10);
    z.initialize_all(6510);
    vx.initialize_all(4);
    vy.initialize_all(8);
    r1.initialize_all(5e3);
    r2.initialize_all(5e3);
    r3.initialize_all(5e3);
    r4.initialize_all(5e3);
    rc1.initialize_all(1./(5e3));
    rc2.initialize_all(1./(5e3));
    rc3.initialize_all(1./(5e3));
    rc4.initialize_all(1./(5e3));
    Mopt.print();
    
    solver<> S(Mopt, ipopt);
    double tol=1e-9, time_limit = 36000;
    S.run(tol=1e-6, time_limit=600);
    Mopt.print_solution();
    return 0;
}

vector<param<double>> myModel::readData(int argc, const char * argv[]){
    N = 4; M = 7; K = 3;
    int degree = 100;

    if(argc>=2){
            string fname = argv[1];
            auto dims = graph.read_pairwise_list(fname);
            N = dims.first;
            M = dims.second;
        }
        else {
            graph.generate_bipartite_random(N,M,degree);
        }

    /* Graph nodes are indexed in {0,...,n+m-1})*/
    assert(N + M == graph.nodes.size());/* Make sure we have the right number of nodes */
    
    //init ownr...
    random_device rd; // obtain a random number from hardware
    mt19937 gen(rd()); // seed the generator
    uniform_int_distribution<> distr(0, K-1); // define the range

    for(int i = 0; i < N; i++)
        owner.push_back(distr(gen)); // generate numbers
    
    /* Indexing sets */
    //sensors = range (0, N - 1);
    //objects = range (N, N + M - 1);
    arcs.add(graph.arcs);
    for (int i = 0; i < N; i++) {
        sensors.add("sensor" + to_string(i));
    }
    for (int i = N; i < N + M; i++) {
        objects.add("object" + to_string(i));
    }
    
    jk = indices(objects, range(0,K-1));
    
    //Parameters
    vector<param<double>> par (3);
    param<double> w0("w0");
    param<double> w_own("w_own");
    param<double> w_bought("w_bought");
    if (argc >= 3) {
        string fname = argv[2];
        fstream file;
        file.open(fname);
        //FILE *fp = fopen(fname.c_str(),"r");
        /*if(f == NULL)
        {
                cout << "Can’t open input file " << fname;
                exit(1);
        }*/
        string tmp;
        string tmp1;
        string tmp2;
        file >> tmp1;
        K = stoi(tmp1);
        for (int i = 0; i < N; i++) {
            file >> tmp2;
            owner[i] = stoi(tmp2);
        }
        for (int k = 0; k < K; k++) {
            for (int i = 0; i < N; i++) {
                if (owner[i] == k) {
                    own_sens.add("sensor" + to_string(i) + "," + to_string(k));
                    for (Arc* a: graph.get_node("sensor" + to_string(i))->get_out()) {
                        own_arcs.add("sensor" + to_string(i) + "," + a->_dest->_name +  "," +  to_string(k));
                        for (Arc* b: graph.get_node(a->_dest->_name)->get_in()) {
                            if ((owner[stoi(b->_src->_name.substr(6, b->_src->_name.find(",")))] == k) && (stoi(b->_src->_name.substr(6, b->_src->_name.find(","))) != i)) {
                                own_rplc.add(a->_src->_name + "," + b->_src->_name + "," + a->_dest->_name +  "," +  to_string(k));
                            }
                            else if (stoi(b->_src->_name.substr(6, b->_src->_name.find(","))) != i) {
                                oths_rplc.add("sensor" + to_string(i) + "," + b->_src->_name + "," + a->_dest->_name +  "," +  to_string(k));
                            }
                        }
                    }
                }
                else {
                    bought_sens.add("sensor" + to_string(i) + "," + to_string(k));
                    for (Arc* a: graph.get_node("sensor" + to_string(i))->get_out()) {
                        bought_arcs.add(a->_src->_name + "," + a->_dest->_name +  "," +  to_string(k));
                    }
                }
            }
        }
        w0.in(arcs);
        w_own.in(own_arcs);
        w_bought.in(bought_arcs);
        for (Arc* a: graph.arcs) {
            for (int k = 0; k < owner[stoi(a->_src->_name.substr(6, a->_src->_name.find(",")))]; k++) {
                file >> tmp;
                w_bought(a->_src->_name + "," + a->_dest->_name + "," + to_string(k)) = stoi(tmp);
            }
            file >> tmp;
            w_own(a->_src->_name + "," + a->_dest->_name + "," + to_string(owner[stoi(a->_src->_name.substr(6, a->_src->_name.find(",")))])) = stod(tmp);
            for (int k = owner[stoi(a->_src->_name.substr(6, a->_src->_name.find(",")))] + 1; k < K; k++) {
                file >> tmp;
                w_bought(a->_src->_name + "," + a->_dest->_name + "," + to_string(k)) = stod(tmp);
            }
        }
        w0.initialize_normal(2, 1);
    }
    /*else {
        w0.initialize_normal(2, 1);
        w_own.initialize_normal(2, 1);
        w_bought.initialize_normal(2, 1);
    }*/
    par[0] = w0;
    par[1] = w_own;
    par[2] = w_bought;
    
    return par;
}

void myModel::InitBilevel(param<double> w0, param<double> w_own, param<double> w_bought) {

    Model<> model("BilevelSensor");
    
    double e = 0.001;

    /*Variables*/
    var<double> p("p", pos_);
    model.add(p.in(sensors));
    var<double> y("y", pos_);
    model.add(y.in(sensors));
    var<int> s("s", 0, 1);
    model.add(s.in(own_arcs));
    var<int> sn("sn", 0, 1);
    model.add(sn.in(sensors));
    var<int> z0("z0", 0, 1);
    model.add(z0.in(arcs));
    var<int> z("z", 0, 1);
    model.add(z.in(bought_arcs));

    var<> obj_var("obj_var", 0,10000);
    model.add(obj_var.in(range(0,0)));
    
    param<> ones("ones");
    ones.in(sensors);
    ones = 1;
    
    
    /*Objective*/
    Constraint<> obj("obj_def");
    obj += product(w_own + w0.in_ignore_ith(2, 1, own_arcs), s);
    obj += ones.tr()*p*sn;
    obj += product(w_bought + w0.in_ignore_ith(2, 1, bought_arcs), z);
    obj -= ones.in_ignore_ith(1, 2, bought_arcs).tr()*p.in_ignore_ith(1, 2, bought_arcs)*z;    
    obj -= e * ones.in_ignore_ith(1, 2, bought_arcs).tr()*p.in_ignore_ith(1, 2, bought_arcs)*z;
    obj -= obj_var;    
    model.add(obj.in(range(0,0)) == 0);
    
    model.max(obj_var);
    
    /*Constraints*/
    //Upper level
    
    Constraint<> ub("Unique_Bought_Obsrvn");
    ub = sum(z.in_matrix(1, 2)) - sn;
    model.add(ub.in(sensors) == 0);

    Constraint<> luo("Leader_Unique_Object");
    luo = sum(z0.in_matrix(0, 1));
    model.add(luo.in(objects) <= 1);
    
    //Lower level
        //----Primal Feasibility----

    Constraint<> fua("Unique_Own_Assignment");
    fua = sum(s.in_matrix(1, 1)) + sn;
    model.add(fua.in(own_sens) == 1);
    
    Constraint<> fuab("Unique_Bought_Assignment");
    fuab = sum(z.in_matrix(1, 1));
    model.add(fuab.in(bought_sens) <= 1);
    
    indices z_ids = z.get_matrix_ids(0, 1);
    Constraint<> fub("Follower_Unique_Object");
    fub = sum(s.sum_over(z_ids, 0)) + sum(z.in(z_ids));
    model.add(fub.in(z_ids.ignore_ith(0, 1)) <= 1);
    
        //----Fair price----
    Constraint<> fp("FairPrice");
    fp = p.in_ignore_ith(1, 2, bought_arcs) - (product(w_bought, z) + y.in_ignore_ith(1, 2, bought_arcs))/2;
//    model.add(fp.in(bought_arcs) >= 0);

    //indices s_ids = s.get_matrix_ids(0, 1);
    //s_ids.print();
    /*Constraint<> sl1("Seller lb1");
    sl1 = y.in_ignore_ith(1, 2, own_arcs) - w_own * (1 - sum(s.sum_over(own_rplc.ignore_ith(0, 1), 0)) - sum(z.sum_over(own_rplc.ignore_ith(0, 1), 0)));
    model.add(sl1.in(own_arcs) >= 0);*/
    
    //w_own.print();
    //w_bought.print();
    for (int i = 0; i < N; i++) {
        for (Arc* b: graph.get_node("sensor" + to_string(i))->get_out()) {
            string j = b->_dest->_name;
            Constraint<> sl1("Seller lb1:" + to_string(i) + "," + b->_dest->_name + "," + to_string(owner[i]));
            sl1 += y("sensor" + to_string(i)) - w_own("sensor" + to_string(i) + "," + b->_dest->_name + "," + to_string(owner[i]));
            for (Arc* a: graph.get_node(j)->get_in()) {
                if ((owner[stoi(a->_src->_name.substr(6, a->_src->_name.find(",")))] == owner[i]) && (stoi(a->_src->_name.substr(6, a->_src->_name.find(","))) != i)) {
                    sl1 += w_own("sensor" + to_string(i) + "," + b->_dest->_name + "," + to_string(owner[i])) * s(a->_src->_name + "," + j +  "," + to_string(owner[i]));
                }
                else if (owner[stoi(a->_src->_name.substr(6, a->_src->_name.find(",")))] != owner[i]) {
                    sl1 += w_own("sensor" + to_string(i) + "," + b->_dest->_name + "," + to_string(owner[i])) * z(a->_src->_name + "," + j +  "," + to_string(owner[i]));
                }
            }
            model.add(sl1 >= 0);
        }
    }
    
    
    Constraint<> sl2("Seller lb2");
    sl2 = y.in_ignore_ith(1, 3, own_rplc) - (w_own.in_ignore_ith(1, 1, own_rplc) - w_own.in_ignore_ith(0, 1, own_rplc)) * s.in_ignore_ith(0, 1, own_rplc);
    model.add(sl2.in(own_rplc) >= 0);
    
    Constraint<> sl3("Seller lb3");
    sl3 = y.in_ignore_ith(1, 3, oths_rplc) - (w_own.in_ignore_ith(1, 1, oths_rplc) - 1)*z.in_ignore_ith(0, 1, oths_rplc);// - ones.in_ignore_ith(1, 3, oths_rplc).tr()*p.in_ignore_ith(1, 3, oths_rplc) * z.in_ignore_ith(0, 1, oths_rplc);
//    model.add(sl3.in(oths_rplc) >= 0);
    
    
    model.print();
    solver<> sol(model, ipopt);
//
    sol.run();
    //model.restructure();
    model.print();

}
