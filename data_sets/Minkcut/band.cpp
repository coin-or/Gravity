/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : band.cpp

* Purpose :

* Creation Date : 15-08-2017

* Last Modified : Tue 15 Aug 09:26:36 2017

* Created By : Guanglei Wang
_._._._._._._._._._._._._._._._._._._._._.*/
#include<iostream>
#include<fstream>
#include<stdio.h>

using namespace std;

int main (int argc, const char * argv[])
{
    int dim = 0;
    if (argc < 1) {
        cerr << "input the size of the graph (vertices)! \n"; 
    }
    else {
       dim = atoi(argv[1]); 
    ofstream outfile("band"+to_string(dim)+".txt");
    if (!outfile)
        cerr << "Oops! Uable to save session data! \n";
    else {
        outfile << dim << " " << dim -1 << "\n";
        for (int i = 1; i < dim; i++)
            outfile << i << " " << i+1 << " " << 1 << "\n";
    }
    }
}
