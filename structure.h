#ifndef _STRUCTURE_H
#define _STRUCTURE_H


#include <vector>
#include "constants.h"


using namespace std;

namespace FDTD
{

	class Medium
	{
	public:
		Medium() {}
		Medium(double aconductivity, double arelPermitivity, double arelPermeability);	//arelPermitivity and arel Permeability are 
																						//relative nondimensional values. aconductivity in SI units.

		double sigma();		//Returns electrical conductivity
		double ep();		//Returns permitivity
		double mu();		//Returns magnetic permeability

	private:

		double elecconductivity;
		double relPermitivity;
		double relPermeability;

		double permitivity;
		double permeability;

	};



	class RectElement
	{
	public:
		RectElement() {}

		RectElement(double axD, double ayD, double azD, double apermit, double apermea, double acond);

		//Rectangular element with dimensions axD, ayD, azD and 						
		//properties permitivity aeps, permeability aper, electrical conductivity acond.

		double sigma() const;		//Returns electrical conductivity
		double eps() const;			//Returns permitivity
		double mu() const;			//Returns permeability

		double xLen() const;				//Returns x-side length
		double yLen() const;				//Returns y-side length
		double zLen() const;				//Returns z-side length

	private:
		double xL;					//x-side length
		double yL;					//y-side length
		double zL;					//z-side length

		double permit;				//permitivity
		double permea;				//permeability
		double cond;				//conductivity
	};


	class Structure1D
	{
	public:
		Structure1D() {}
		Structure1D(int assgnumLayers, double aWidth, double assgeps[], double assgmu[], double assgsig[], double assgthickness[]);

		double th(int i) const;
		double ep(int i) const;
		double mu(int i) const;
		double sigma(int i) const;

		double min_c() const;		//Returns minimum speed of electromagnetic radiation in any of the layers.

		int numL() const;
		double width() const;	//Width of device

	private:
		int numLayers;			//Number of layers in device
		double myWidth;

		double* thickness;
		double* eps;
		double* muu;
		double* sig;

		double minSpeedLight;		//Minimum speed of electromagnetic radiation in any of the layers.

	};
	


//*****************************************************************************************************************************************************
//class for three dimensional structures. An object of this class represents a structure consisting of a rectangular latitice of elements. 
//The lattice is formed by decomposing a rectangular 3D region with edges parallel to the x,y and z directions into layers perpendicular to 
//the x-direction, layers perpendicular to the y-direction and layers perpendicular to the z-direction.
//Letting nLx be the number of layers perpendicular to the x-direction, these layers are identified by an index i = 0, 1, 2, ....nLX - 1. 
//Letting nLy be the number of layers perpendicular to the y-direction, these regions are identified by an index j = 0, 1, 2, ...nLY - 1.
//Letting nLz be the number of layers perpendicular to the y-direction, these regions are identified by an index j = 0, 1, 2, ...nLY - 1.
//Thus each element in the structure is identified by an ordered triple of indexes i,j, k.

//Each element has electrical and other properties, including permitivity, permeability and electrical conductivity. 
//Groups of elements may be composed of the same material. Thus for example, all of the elements in an entire x-direction layer may be of the 
//same material, i.e. have the same properties, so that the entire layer is of the same material. 



	class Structure3D
	{
	public:
		Structure3D() {}

		Structure3D(vector <double>& axLThick, vector <double>& ayLThick, vector <double>& azLThick, vector <vector <vector <RectElement>>>& aelement, bool uC, bool yDP, bool zDP);
					//Constructs a Structure3D  of the type described above with the vector axLThick of x-layer thicknesses, the vector ayLThick of y-layer thicknesses,
					//the vector azLThick of z-layer thicknesses and the three-dimensional vector of rectangular elements aelement. 
					//The elements of aelement specify the material properties of the elements of the structure. 
					//If argument unitCell is true, this structure is taken to be the unit cell of a periodic structure. 

		Structure3D(Structure3D& strTobeCopied);

		//Accessors
		int numLx() const;					//Returns the number of layers in the x-direction.	
		int numLy() const;					//Returns the number of layers in the y-direction.		
		int numLz() const;					//Returns the number of layers in the z-direction.

		double xLth(int i) const;			//Returns the thickness of the layer perpendicular to the x-direction with index i.
		double yLth(int j) const;			//Returns the thickness of the layer perpendicular to the y-direction with index j.
		double zLth(int k) const;			//Returns the thickness of the layer perpendicular to the z-direction with index k.

		double sigma(int i, int j, int k) const;	//Returns electrical conductivity of element i,j,k
		double eps(int i, int j, int k) const;		//Returns permitivity of element i,j,k
		double mu(int i, int j, int k) const;		//Returns permrability of element i,j,k

		bool Is_a_unitCell() const;					//Returns true if the structure is the unit cell of a larger periodic structure
		bool Is_periodic_in_yDirection() const;		//Returns true if the larger structure is periodic in the y-direction
		bool Is_periodic_in_zDirection() const;		//Returns true if the larger structure is periodic in the z-direction

		double min_c() const;		//Returns minimum speed of electromagnetic radiation in any of the layers.


	private:

		int nLx;						//Number of layers in the x-direction
		int nLy;						//Number of layers in the y-direction
		int nLz;						//Number of layers in the z-direction

		double txD;						//Total depth of structure in the x-direction
		double tyD;						//Total hieght of structure in the y-direction
		double tzD;						//Total width of structure in the z-direction



		vector <double> xLThick;		//Vector of x-layer thicknesses
		vector <double> yLThick;		//Vector of y-layer thicknesses
		vector <double> zLThick;		//Vector of z-layer thicknesses

		vector <vector <vector <RectElement>>> element;		//3D array of rectaangular elements. See RectElement class.

		bool unitCell;		//Set to true if object represents a unit cell of periodic structure
		bool yPeriodicity;	//Set to true if y-direction periodicity;
		bool zPeriodicity;	//Set to true if z-direction periodicity;

	};



}//End namespace FDTD


#endif
