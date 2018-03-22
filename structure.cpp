//
//

#include "structure.h"

#include <stdexcept>


using namespace std;


namespace FDTD
{
	//***********************************************************************************************************
	//member functions for class Medium

	Medium::Medium(double aconductivity, double arelPermitivity, double arelPermeability)
	{
		elecconductivity = aconductivity;
		relPermitivity = arelPermitivity;
		relPermeability = arelPermeability;

		permitivity = eps0*relPermitivity;
		permeability = mu0*relPermeability;

	}


	//Accessors for class Medium

	double Medium::sigma() { return elecconductivity; }
	double Medium::ep() { return permitivity; }
	double Medium::mu() { return permeability; }



	//*************************************************************************************************************
	//Member functions for class RectElement

	RectElement::RectElement(double axL, double ayL, double azL, double apermit, double apermea, double acond)
	{
		xL = axL;
		yL = ayL;
		zL = azL;

		permit = apermit;
		permea = apermea;
		cond = acond;
	}//End constructor


	 //Accessors
	double RectElement::sigma() const { return cond; }		//Returns electrical conductivity
	double RectElement::eps() const { return permit; }		//Returns permitivity
	double RectElement::mu() const { return permea; }		//Returns permeability


	double RectElement::xLen() const { return xL; }			//Returns x-side length

	double RectElement::yLen() const { return yL; }			//Returns y-side length
	double RectElement::zLen() const { return zL; }			//Returns z-side length



	//**************************************************************************************************************
	//member functions for class Structure3D

	Structure3D::Structure3D(vector <double>& axLThick, vector <double>& ayLThick, vector <double>& azLThick, vector <vector <vector <RectElement>>>& aelement, bool uC, bool yDP, bool zDP) :
		nLx(axLThick.size()), nLy(ayLThick.size()), nLz(azLThick.size()),
		xLThick(axLThick), yLThick(ayLThick), zLThick(azLThick), element(aelement)
	{
		int i, j, k;

		nLx = axLThick.size();
		nLy = ayLThick.size();
		nLz = azLThick.size();

		txD = 0;
		for (i = 0; i < nLx; i++)
			txD += xLThick[i];
		//End for

		tyD = 0;
		for (j = 0; j < nLy; j++)
			tyD += yLThick[j];
		//End for

		tzD = 0;
		for (k = 0; k < nLz; k++)
			tzD += zLThick[k];
		//End for

		unitCell = uC;
		yPeriodicity = yDP;
		zPeriodicity = zDP;

	}//End constructor

	Structure3D::Structure3D(Structure3D& strTobeCopied)
	{
		nLx = strTobeCopied.nLx;
		nLy = strTobeCopied.nLy;
		nLz = strTobeCopied.nLz;

		txD = strTobeCopied.txD;
		tyD = strTobeCopied.tyD;
		tzD = strTobeCopied.tzD;

		xLThick = strTobeCopied.xLThick;
		yLThick = strTobeCopied.yLThick;
		zLThick = strTobeCopied.zLThick;

	}//End copy construdtor


	 //Accessors

	double Structure3D::sigma(int i, int j, int k) const { return element[i][j][k].sigma(); }	//Returns electrical conductivity of element i,j,k
	double Structure3D::eps(int i, int j, int k) const { return element[i][j][k].eps(); }		//Returns permitivity of element i,j,k
	double Structure3D::mu(int i, int j, int k) const { return element[i][j][k].mu(); }			//Returns permrability of element i,j,k

	int Structure3D::numLx() const { return nLx; }			//Returns the number of layers in the x-direction.	
	int Structure3D::numLy() const { return nLy; }			//Returns the number of layers in the y-direction.		
	int Structure3D::numLz() const { return nLz; }			//Returns the number of layers in the z-direction.

	double Structure3D::xLth(int i) const { return xLThick[i]; }		//Returns the thickness of the layer perpendicular to the x-direction with index i.
	double Structure3D::yLth(int j) const { return yLThick[j]; }		//Returns the thickness of the layer perpendicular to the y-direction with index j.
	double Structure3D::zLth(int k) const { return zLThick[k]; }		//Returns the thickness of the layer perpendicular to the z-direction with index k.

	bool Structure3D::Is_a_unitCell() const { return unitCell; }					//Returns true if the structure is the unit cell of a larger periodic structure
	bool Structure3D::Is_periodic_in_yDirection() const { return yPeriodicity; }		//Returns true if the larger structure is periodic in the y-direction
	bool Structure3D::Is_periodic_in_zDirection() const { return zPeriodicity; }	//Returns true if the larger structure is periodic in the z-direction

	double Structure3D::min_c() const 
	{
		int i, j, k;
		double ci;
		double minSpeedLight;

		minSpeedLight = 1 / sqrt(element[0][0][0].eps( )*element[0][0][0].mu());

		for (i = 0; i < nLx; i++)
			for (j = 0; j < nLy; j++)
				for (k = 0; k < nLz; k++)
					{
						ci = 1 / sqrt(element[i][j][k].eps()*element[i][j][k].mu());
						if (ci < minSpeedLight)
						minSpeedLight = ci;
					}//End for

		return minSpeedLight;

	}		//Returns minimum speed of electromagnetic radiation in any of the layers.



//**************************************************************************************************************************************************************
//Member functions for class Device1D


	Structure1D::Structure1D(int assgnumLayers, double awidth, double assgeps[], double assgmu[], double assgsig[], double assgthickness[])
	{
		int i;
		double ci;

		numLayers = assgnumLayers;
		myWidth = awidth;

		thickness = new double[numLayers];
		eps = new double[numLayers];
		muu = new double[numLayers];
		sig = new double[numLayers];

		for (i = 0; i < numLayers; i++)
		{
			thickness[i] = assgthickness[i];	// *1E-6;	//Converted to meters
			eps[i] = assgeps[i];
			muu[i] = assgmu[i];
			sig[i] = assgsig[i];
		}//end for each layer

		minSpeedLight = 1 / sqrt(eps[0] * muu[0]);

		for (i = 1; i < numLayers; i++)
		{
			ci = 1 / sqrt(eps[i] * muu[i]);
			if (ci < minSpeedLight)
				minSpeedLight = ci;
		}//End for

	}//End Constructor


	 //Accessors for class Device1D
	double Structure1D::th(int i) const { return thickness[i]; }	//returns thickness of sub i layer of device

	double Structure1D::ep(int i) const { return eps[i]; }

	double Structure1D::mu(int i) const { return muu[i]; }

	double Structure1D::sigma(int i) const { return sig[i]; }

	int Structure1D::numL() const { return numLayers; }

	double Structure1D::width() const { return myWidth; }

	double Structure1D::min_c() const { return minSpeedLight; }		//Returns minimum speed of electromagnetic radiation in any of the layers.
	
	
	

}//End namespace FDTD