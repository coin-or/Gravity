#include <iostream>
#include <stdio.h>
#include <math.h>
#include <list>
#include <gravity/matplotlibcpp.h>
#include <gravity/solver.h>
#include <gravity/rapidcsv.h>
#include <DataSet.h>
#include <time.h>
using namespace std;

int main ()
{
    vector<double> x_vec, y_vec, z_vec;
    double x1 = 1;
    double y1 = 1;
    double z1 = 1;
    double x_rot1, y_rot1, z_rot1;
    double angles[] = {0, 0.1, -0.1};
//      list<double> mylist (angles,angles+3);
    for (int a = 0; a < 3; a++){
        for (int b = 0; b <3; b++){
            for (int c = 0; c <3; c++){

                    x_rot1 = x1*cos(angles[a])*cos(angles[b]) + y1*(cos(angles[a])*sin(angles[b])*sin(angles[c]) - sin(angles[a])*cos(angles[c])) + z1*(cos(angles[a])*sin(angles[b])*cos(angles[c]) + sin(angles[a])*sin(angles[c]));

                   y_rot1 = x1*sin(angles[a])*cos(angles[b]) + y1*(sin(angles[a])*sin(angles[b])*sin(angles[c]) + cos(angles[a])*cos(angles[c])) + z1*(sin(angles[a])*sin(angles[b])*cos(angles[c]) - cos(angles[a])*sin(angles[c]));

                   z_rot1 = x1*sin(-1*angles[b]) + y1*(cos(angles[b])*sin(angles[c])) + z1*(cos(angles[b])*cos(angles[c]));

                   x_vec.push_back(x_rot1);
                   y_vec.push_back(y_rot1);
                   z_vec.push_back(z_rot1);
          }}}

    


namespace plt = matplotlibcpp;
   
    std::map<std::string, std::string> keywords, keywords2;
    keywords["marker"] = "s";
    keywords["linestyle"] = "None";
    keywords["ms"] = "0.05";
    plt::plot3(x_vec, y_vec, z_vec, keywords);

    keywords2["marker"] = "s";
    keywords2["ms"] = "0.1";

    plt::show();

}
