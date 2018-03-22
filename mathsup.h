#ifndef _MATHSUP_H
#define _MATHSUP_H

#include <vector>

using namespace std;

namespace FDTD

{
	double max(double x, double y);

	double min(double x, double y);


	//**************************************************************************************************************************
	class Array1D
	{
	public:
		Array1D() {};
		Array1D(int M);

		double& operator() (const int&);

	private:
		vector<double> a;

	};


	//**************************************************************************************************************************
	class Array1D_Vector3D
	{
	public:
		Array1D_Vector3D(){}
		Array1D_Vector3D(int numElements, int Size1, int Size2, int Size3);

		vector <vector <vector<double>>>& operator() (const int&);

	private:
		vector<vector<vector<vector<double>>>> a;

	};


	//**************************************************************************************************************************
	class Array1D_Vector1D
	{
	public:
		Array1D_Vector1D() {}
		Array1D_Vector1D(int numElements, int Size1);

		vector<double>& operator() (const int&);

	private:
		vector<vector<double>> a;

	};


	//**************************************************************************************************************************
	class Array2D
	{
	public:
		Array2D() {};
		Array2D(int M, int N);

		double& operator() (const int&, const int&);

	private:
		vector <vector<double>> a;

	};



	//**************************************************************************************************************************
	class Array3D
	{
	public:
		Array3D() {};
		Array3D(int numX, int numY, int numZ);

		double& operator() (const int&, const int&, const int&);

	private:
		vector <double> v;
		vector <vector <vector<double>>> a;

	};


}//End namespace FDTD


#endif
