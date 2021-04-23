/********************************************************************
Header File for Go-ICP Class
Last modified: Apr 21, 2014

"Go-ICP: Solving 3D Registration Efficiently and Globally Optimally"
Jiaolong Yang, Hongdong Li, Yunde Jia
International Conference on Computer Vision (ICCV), 2013

Copyright (C) 2013 Jiaolong Yang (BIT and ANU)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*********************************************************************/

#ifndef JLY_GOICP_H
#define JLY_GOICP_H

#include <queue>
using namespace std;

#include "jly_icp3d.hpp"
#include "jly_3ddt.h"

#define PI 3.1415926536
#define SQRT3 1.732050808

typedef struct _POINT3D
{
	double x, y, z;
}POINT3D;

typedef struct _ROTNODE
{
	double a, b, c, w;
	double ub, lb;
	int l;
	friend bool operator < (const struct _ROTNODE & n1, const struct _ROTNODE & n2)
	{
		if(n1.lb != n2.lb)
			return n1.lb > n2.lb;
		else
			return n1.w < n2.w;
			//return n1.ub > n2.ub;
	}
	
}ROTNODE;

typedef struct _TRANSNODE
{
	double x, y, z, w;
	double ub, lb;
	friend bool operator < (const struct _TRANSNODE & n1, const struct _TRANSNODE & n2)
	{
		if(n1.lb != n2.lb)
			return n1.lb > n2.lb;
		else
			return n1.w < n2.w;
			//return n1.ub > n2.ub;
	}
}TRANSNODE;

/********************************************************/



/********************************************************/

#define MAXROTLEVEL 20

class GoICP
{
public:
	int Nm, Nd;
	POINT3D * pModel, * pData;

	ROTNODE initNodeRot;
	TRANSNODE initNodeTrans;

	DT3D dt;

	ROTNODE optNodeRot;
	TRANSNODE optNodeTrans;

	GoICP();
	double Register();
   
	void BuildDT();

	double MSEThresh;
	double SSEThresh;
	double icpThresh;

	double optError;
	Go_ICP::Matrix optR;
	Go_ICP::Matrix optT;
    
    Go_ICP::Matrix R_init;
    Go_ICP::Matrix T_init;
    
    double run_ICP();

	clock_t clockBegin;

	double trimFraction;
	int inlierNum;
	bool doTrim;

private:
	//temp variables
	double * normData;
	double * minDis;
	double** maxRotDis;
	double * maxRotDisL;
	POINT3D * pDataTemp;
	POINT3D * pDataTempICP;
	
	ICP3D<double> icp3d;
	double * M_icp;
	double * D_icp;

	double ICP(Go_ICP::Matrix& R_icp, Go_ICP::Matrix& t_icp);
	double InnerBnB(double* maxRotDisL, TRANSNODE* nodeTransOut, int& nb_nodes);
	double OuterBnB();
	void Initialize();
	void Clear();

};

/********************************************************/

#endif
