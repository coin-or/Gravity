/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : band.cpp

* Purpose :

* Creation Date : 15-08-2017

* Last Modified : Tue 15 Aug 10:33:07 2017

* Created By : Guanglei Wang
_._._._._._._._._._._._._._._._._._._._._.*/
#include<iostream>
#include<fstream>
#include<string>
#include<string.h>
#include<stdio.h>
#include<stdlib.h>

using namespace std;

int main (int argc, const char * argv[])
{
    int dim = 0;
    if (argc < 1) {
        cerr << "input the size of the graph (vertices)! \n"; 
    }
    else {
       dim = atoi(argv[1]); 
    	ofstream outfile("band"+to_string(dim)+".txt", ios_base::out);
    if (!outfile)
        cerr << "Oops! Uable to save session data! \n";
    else {
        outfile << dim << " " << 3*dim -6 << "\n";
        int a = 1;
        for (int i = 1; i < dim - 2; i++){
            outfile << i << " " << i+1 << " " << a << "\n";
            outfile << i << " " << i+2 << " " << a << "\n";
            outfile << i << " " << i+3 << " " << a << "\n";
            a = -1*a;
        }
            outfile << dim - 2 << " " << dim - 1 << " " << a << "\n";
            outfile << dim - 2 << " " << dim  << " " << a << "\n";
            outfile << dim -1 << " " << dim  << " " << a << "\n";
    }
    outfile.close();
    }
    return 0;
}
