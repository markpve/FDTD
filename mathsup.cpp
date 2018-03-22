//
//

#include "mathsup.h"

#include <stdexcept>

using namespace std;

namespace FDTD
{
	double max(double x, double y)
	{
		if (x >= y)
			return x;
		else
			return y;
	}

	double min(double x, double y)
	{
		if (x <= y)
			return x;
		else
			return y;
	}

	//**************************************************************************************************************
	//Functions for class Array1D

	Array1D::Array1D(int M) : a(M, 0.0)
	{

	}

	double& Array1D::operator() (const int& m) { return a[m]; }



	//******************************************************************************************************************
	//functions for class Array1D_Vector3D

	Array1D_Vector3D::Array1D_Vector3D(int numElements, int Size1, int Size2, int Size3) : 
		a(numElements, vector<vector <vector<double>>>(Size1, vector<vector<double>>(Size2, vector<double>(Size3, 0))))
	{

	}

	vector<vector<vector<double>>>& Array1D_Vector3D::operator() (const int& m) { return a.at(m); }



	//********************************************************************************************************************
	//Functions for class Array1D_Vector1D

	Array1D_Vector1D::Array1D_Vector1D(int numElements, int Size1) : a(numElements, vector<double>(Size1, 0))
	{

	}

	vector<double>& Array1D_Vector1D::operator() (const int& m) { return a.at(m); }



	//*********************************************************************************************************************
	//Functions for class Array2D

	Array2D::Array2D(int M, int N) : a(M, vector<double>(N, 0))
	{

	}

	double& Array2D::operator() (const int& m, const int& n) { return a[m][n]; }



	//**********************************************************************************************************************
	//Functions for class Array3D

	Array3D::Array3D(int numX, int numY, int numZ) : a(numX, vector <vector<double>>(numY, vector<double>(numZ, 0)))
	{

	}

	double& Array3D::operator() (const int& m, const int& n, const int& p) { return a[m][n][p]; }

}//End namespace FDTD