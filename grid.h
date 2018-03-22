#ifndef _GRID_H
#define _GRID_H

#include <vector>
#include "structure.h"
#include "mathsup.h"
#include "constants.h"
#include "report.h"


using namespace std;

namespace FDTD
{
	const int minCpt = 20;			//Minumum number of spacial steps across a layer
	const int NLambda = 60;			//Minimum number of spacial steps per wavelength

	


	class Grid3D
	{
	public:
		Grid3D() {}
		//Default constructor

		Grid3D(Structure3D& str, double Lambda, Report& myReport);

		//An object of the class represents a three-dimensional multilayer grid.  
		//The grid is generated to fit the structure represented by the object of the class Structure3D referenced by str. 
		//x-direction grid layers match structure layers except grid has 
		//one layer on top (left side of grid) of structure and one layer below (right side) the structure.
		//deltx below the structure = deltx above the structure.

		//Lambda is the spacial width of the source pulse or wavelength of a harmonic source and is used by Grid3D to deternine spacial step sizes

		//To fit the structure, the grid is divided into the same planar regions perpendicular to the x,y and z-directions as the structure, 
		//thereby forming a rectangular three-dimensional array of rectangular elements that coincide with the structure elements. The grid also has 
		//one layer on top (low-x side of grid) of structure and one layer below (high-x side) the structure. 
		//Let nLx be the number of grid layers perpendicular to the x-direction. These layers are identified by an index i = 0, 1, 2, ....nLx - 1. 
		//Let nLy be the number of grid layers perpendicular to the y-direction. These layers are identified by an index j = 0, 1, 2, ...nLy - 1. 
		//Let nLz be the number of grid layers perpendicular to the z-direction. These layers are identified by an index k = 0, 1, 2, ...nLz - 1. 
		//Thus each element in the grid is identified by an ordered triple of indexes i,j, k and the geometric properties of the elements form arrays of values.

		//The rectangular elements of the grid are further divided into Yee cells.

		//The x-direction node spacing may only vary in the x-direction from one layer to another. Thus all of the elements in a given layer perpendicular 
		//to the x-direction have the same x-direction node spacing.
		//Similarily, the y-direction node spacing may only vary in the y-direction from one planar region to another. Thus all of the elements in a given layer perpendicular 
		//to the y-direction have the same y-direction node spacing.
		//Similarily, the z-direction node spacing may only vary in the z-direction from one planar region to another. Thus all of the elements in a given layer perpendicular 
		//to the z-direction have the same z-direction node spacing.

		//The x-direction Layer above (low-x side of) the structure begins on the grid boundary, which contains Ez, Hx and Ey nodes, and ends on a plane of Ez, Hx and Ey nodes nodes 
		//at the top surface of structure. 
		//All x-layers of the grid extend between planes of Ez, Hx and Ey nodes.

		//The y-direction Layer on the low y-side of the structure begins on the grid boundary, which contains Ez, Ex and Hy nodes at y = 0, and ends on a plane of Ez, Ex and Hy nodes nodes. 
		//All y-layers of the grid extend between planes of Ez, Ex and Hy nodes. 

		//The z-direction Layer on the low z-side of the structure begins on the grid boundary, which contains Ey, Ex and Hz nodes at z = 0, and ends on a plane of Ey, Ex and Hz nodes. 
		//All z-layers extend between planes of Ey, Ex and Hz nodes nodes.

		//With the grid structure described here, each grid element i,j,k contains some number of Yee cells. The interfaces between layers contain Yee cell faces. 
		//None of the Yee cells are cut by interfaces bewteen layers.

		//The number of spacial steps across an interior layer is the maximum of (width of the layer/Lambda)*Nlambda rounded to an integer and the constant minCpt, 
		//where minCpt is the minimum number rerquired for any layer.

		//The spacial step size across an interior layer is width of the layer/ Number of spacial steps. So if (width of the layer/Lambda)*Nlambda is an integer > minCpt, the spacial step size will be Lambda/Nlambda.
		//If (width of the layer/Lambda)*Nlambda is an not an integer but its rounded value is > minCpt, step size is width of the layer/Rounded((width of the layer/Lambda)*Nlambda).
		//If the layer is thin enought that Rounded ((width of the layer/Lambda)*Nlambda) is < minCpt, then the spacial step size for the given direction in the layer will be width of the layer/minCpt.
		//So relatively very thin layers may have much smaller node spacing than their adjacent layers.
		//Thus the spacial step size in any layer of the Grid3D will be Lambda/Nlamda rounded so that the number of spacial steps across the layer is an integer, unless
		//the layer is so thin that the number of spacial steps would then be less than the constant minCpt defined above. 

		//So to fit a layer in the structure, the grid in the layer may be stretched or compressed compared to that in other layers.

		//The number of spacial steps across the x-direction grid layer above the structure will be NLambda
		//The spacial step size in the x-direction grid layer above the structure will be the same as in the first x-direction layer in the structure.
		//So the thickness of the x-direction grid layer above the structure will be Lambda if the step size in the first structure layer is Lambda/Nlamda.
		//If the number of steps across the first x-direction structure layer is minCpt, 
		//the thickness of the x-direction grid layer above the structure will be Nlamda*(thickness of the first structure layer)/minCpt.

		//The number of spacial steps across the x-direction grid layer below the structure will be NLambda.
		//The spacial step size in the x-direction grid layer below the structure will be the same as in the last x-direction layer in the structure.

		//Bottom layer may be used for a PML to terminate the grid or a vacuum layer with an ABC on the right.


		Grid3D(Grid3D& aGrid);			//Copy constructor

		int sizeX() const;				//Returns the number of Ez nodes on a line in x-direction across grid.
		int sizeY() const;				//Returns the number of Ey nodes in y-direction across grid.
		int sizeZ() const;				//Returns the number of Ex nodes in z-direction across grid.

		int numLx() const;				//Returns the number of layers normal to the x-direction.
		int numLy() const;				//Returns the number of layers in the y-direction.		
		int numLz() const;				//Returns the number of layers in the z-direction.

		int nStepsXL(int i) const;		//Returns the number of  spacial steps across x-layer with index i on a line in the x-direction. 
		int nStepsYL(int j) const;		//Returns the number of spacial steps accross y-layer with index j on a line in the y-direction.
		int nStepsZL(int k) const;		//Returns the number of spacial steps accross z-layer with index k on a line in the z-direction.

		double deltaxL(int i) const;	//deltaxL(i) returns the x-direction node spacing in micrometers for x-layer with index i.
		double deltayL(int j) const;	//deltayL(j) returns the y-direction node spacing in micrometers for y-layer with index j.
		double deltazL(int k) const;	//deltazL(i) returns the z-direction node spacing in micrometers for z-layer with index k.

		double minDelta() const;				//Returns the minimum of spacial steps (minimum Yee cell edge length)

		int numSteps_xDir_throuh(int indlast) const;		//Returns the number of steps across the grid in the x-direction beginning at x = 0 up through end of layer with index indLast.
															//If indlast < 0, returns 0

		int numSteps_yDir_throuh(int indLast) const;		//Returns the number of steps across the grid in the y-direction beginning at y = 0 up through end of layer with index indLast.
															//If indlast < 0, returns 0

		int numSteps_zDir_throuh(int indlast) const;		//Returns the number of steps across the grid in the z-direction beginning at z = 0 up through end of layer with index indLast.
															//If indlast < 0, returns 0

		double delAx(int mInd, int nInd, int pInd) const;		//Returns the x-face area of the Yee cell with indexes mInd, nInd, pInd

		bool find(int m, int n, int p, int& i, int& j, int& k, bool& interior, bool& xLower, bool& yLower, bool& zLower, bool& xBoundary, bool& yBoundary, bool& zBoundary) const;

		//If the Yee cell m,n,p is in the grid, find finds the grid element containing the Yee cell m,n,p and sets indexes iInd,jInd, kInd to the grid element's indexes. Returns true if the element found.
		//If the Yee cell sits on the high x-boundary or high y-boundary or high z-boundary (which is the case if for example the Ex node of the Yee cell is on the high z grid boundary), 
		//find sets the indexes iInd, JInd, kInd and returns true. 
		//If find returns false, then error condition.
		//If m is the first node index in the x-layer with layer index i, so that the Yee cell m,n,p is on the x = const lower boundary of the layer with index i, find sets xLower to true. 
		//If n is the first node index in the y-layer with layer index j, so that the Yee cell m,n,p is on the y = const lower boundary of the layer with index j, find sets yLower to true. 
		//If p is the first node index in the z-layer with layer index k, so that the Yee cell m,n,p is on the z = const lower boundary of the layer with index k, find sets zLower to true.
		//If m is the last x-node index  in the grid, so that the Yee cell sits on the high x-boundary, find sets xBoundary to true.
		//If n is the last y-node index  in the grid, so that the Yee cell sits on the high y-boundary, find sets yBoundary to true.
		//If p is the last z-node index  in the grid, so that the Yee cell sits on the high z-boundary, find sets zBoundary to true.


	private:

		int M;							//Number of Ez node planes across grid in x-direction including the boundary surfaces
		int N;							//Number of Ez node planes across grid in y-direction including the boundary surfaces
		int P;							//Number of Ex node planes across grid in z-direction including the boundary surfaces

		int mD;							//Device/structure starts at x = mD*deltx[0].

		int nLx;						//Number of layers in the x-direction
		int nLy;						//Number of layers in the y-direction
		int nLz;						//Number of layers in the z-direction

		vector <double> xLThick;		//Vector of x-layer thicknesses
		vector <double> yLThick;		//Vector of y-layer thicknesses
		vector <double> zLThick;		//Vector of z-layer thicknesses

		vector <int> numStepsXL;		//numStepsXL[i] is the number spacial steps across layer with index i on a line in the x-direction. 
		vector <int> numStepsYL;		//numStepsYL[j] is the number of spacial steps across layer with index j in the y-direction.
		vector <int> numStepsZL;		//numStepsZL[k] is the number of spacial steps across layer with index k in the z-direction.


		vector <double> deltxL;			//deltxL[i] is the x-direction node spacing in micrometers for x-layer with index i.
		vector <double> deltyL;			//deltyL[j] is the y-direction node spacing in micrometers for y-layer with index j.
		vector <double> deltzL;			//deltzL[k] is the z-direction node spacing in micrometers for z-layer with index k.

		double minDx;					//Minimum deltxL over all x-layers
		double minDy;					//Minimum deltyL over all y-layers		
		double minDz;					//Minimum deltzL over all z-layers

		double minDel;					//Minimum of all spacial increments


	};


	
	//*****************************************************************************************************************************************************


	//Enumerated data types used in the class Structure3DWithGrid

	enum CType_Ex{ Interior_Ex, yInterface_Ex, zInterface_Ex, yzEdge_Ex };			//A vector CExE_for has four elements: the elemnt with index Interior_Ex = 0 is the array of interior coefficients etc. 

	enum CType_Ey { Interior_Ey, xInterface_Ey, zInterface_Ey, xzEdge_Ey};

	enum CType_Ez { Interior_Ez, xInterface_Ez, yInterface_Ez, xyEdge_Ez};

	enum CType_Hz{ Interior_Hz, zInterface_Hz};

	enum CType_Hy { Interior_Hy, yInterface_Hy};

	enum CType_Hx { Interior_Hx, xInterface_Hx};

	enum CB_High_y_Type_Ex {Interior_of_yBoundary_element_Ex, zInterface_on_yBoundary_element_Ex};

	enum CB_High_z_Type_Ex {Interior_of_zBoundary_element_Ex, yInterface_on_zBoundary_element_Ex };

	enum CB_High_z_Type_Ey {Interior_of_zBoundary_element_Ey, xInterface_on_zBoundary_element_Ey};

	enum CB_High_x_Type_Ey { xBoundary_Ey, xBoundary_zInterface_Ey};

	enum CB_High_y_Type_Ez { Interior_of_yBoundary_element_Ez, xInterface_on_yBoundary_element_Ez};

	//enum CBType_Ez { xBoundary_Ez, Interior_of_yBoundary_element_Ez, xBoundary_yInterface_Ez, xInterface_on_yBoundary_element_Ez };

	//enum CBBType_Ex { yBoundary_Ex, zBoundary_Ez, zBoundary_yInterface_Ex, zBoundary_yBoundary_Ex };

	

	class Structure3DWithGrid
	{
	public:
		Structure3DWithGrid() {}
		//Default constructor

		Structure3DWithGrid(Structure3D& str, double Lambda, Report& myReport);

		//Constructs object representing the three dimensional structure str with grid fitted to it. An object of this class represents a structure consisting of a rectangular latitice of elements 
		//as described in the documentation for the class Structure3D and its grid, generated based on Lambda as described in the documentation for the class Grid3D.
		//The structure occupies a rectangular 3D region with edges parallel to the x,y and z directions.

		//Once the structure is defined or copied from a Structurre3D object, the grid is generated to fit it as described in the Grid3D documentation.
		//The grid  has x-layers above and below the structure as described in the Grid3D documentation.
		//The grid is divided into layers perpendicular to the x-direction, layers perpendicular to the y-direction and layers perpendicular to the z-direction.
		//Letting nLx be the number of layers perpendicular to the x-direction, these layers are identified by an index i = 0, 1, 2, ....nLX - 1. 
		//Letting nLy be the number of planar regions perpendicular to the y-direction, these regions are identified by an index j = 0, 1, 2, ...nLY - 1.
		//Letting nLz be the number of planar regions perpendicular to the z-direction, these regions are identified by an index k = 0, 1, 2, ...nLZ - 1.
		//Thus the grid is decomposed into elements with each element identified by an ordered triple of indexes i,j, k. Each element has a grid within it as
		//described in the Grid3D documentaion.

		//Each element has electrical and other properties, including permitivity, permeability and electrical conductivity. 
		//Groups of elements of the grid within the structure portion may be composed of the same material. Thus for example, all of the elements in an entire x-direction layer may be of the 
		//same material, i.e. have the same properties, so that the entire layer is of the same material.  The x-direction grid layers above and below the structure 
		//each are assumed to be vacuum layers and the update coefficients in those layers are computed accordingly.

		//After the grid is generated, the update coefficients are computed using 
		//delt = theGrid.minDelta()*Sc / c0 with Sc = 1 / sqrt(3.0)			//delt is in microseconds and delx in micrometers. 

		Structure3DWithGrid(vector <double>& axLThick, vector <double>& ayLThick, vector <double>& azLThick, vector <vector <vector <RectElement>>>& aelement, double Lambda,
			bool uC, bool yDP, bool zDP, Report& myReport);


		//Constructs a Structure3D  of the type described above with the vector axLThick of x-layer thicknesses, the vector ayLThick of y-layer thicknesses,
		//and the vector azLThick of z-layer thicknesses. If uC is true, the structure is taken to be the unit cell of a larger structure periodic in the y and z directions. 
		//After the structure is defined, a Grid3D is generated as described above to fit the structure. 
		
		//After the grid is generated, the update coefficients are computed using 
		//delt = theGrid.minDelta()*Sc / c0 with Sc = 1 / sqrt(3.0)			//delt is in microseconds and delx in micrometers.


		Grid3D* grid();				//Returns pointer to the Grid3D member object
		Structure3D* structure();	//Returns pointer to the structure member object

		int numLx() const;				//Returns the number of grid layers normal to the x-direction. (two more than the number of x-direction layers in the structure.)
		int numLy() const;				//Returns the number of layers in the y-direction.		
		int numLz() const;				//Returns the number of layers in the z-direction.

		double deltaxL(int i) const;	//deltaxL(i) returns the x-direction node spacing in micrometers for x-layer with index i.
		double deltayL(int j) const;	//deltayL(j) returns the y-direction node spacing in micrometers for y-layer with index j.
		double deltazL(int k) const;	//deltazL(i) returns the z-direction node spacing in micrometers for z-layer with index k.

		int nStepsXL(int i) const;		//Returns the number of  spacial steps across x-layer with index i on a line in the x-direction. 
		int nStepsYL(int j) const;		//Returns the number of spacial steps accross y-layer with index j on a line in the y-direction.
		int nStepsZL(int k) const;		//Returns the number of spacial steps accross z-layer with index k on a line in the z-direction.

		int sizeX() const;				//Returns the number of Ez nodes on a line in x-direction across grid.
		int sizeY() const;				//Returns the number of Ez nodes in y-direction across grid.
		int sizeZ() const;				//Returns the number of Ex nodes in z-direction across grid.

		double delAx(int mInd, int nInd, int pInd) const;		//Returns the x-face area of the Yee cell with indexes mInd, nInd, pInd

		int numSteps_xDir_throuh(int indlast) const;		//Returns the number of steps across the grid in the x-direction beginning at x = 0 up through end of layer with index indLast.

		int numSteps_yDir_throuh(int indlast) const;		//Returns the number of steps across the grid in the y-direction beginning at y = 0 up through end of layer with index indLast.

		int numSteps_zDir_throuh(int indlast) const;		//Returns the number of steps across the grid in the z-direction beginning at z = 0 up through end of layer with index indLast.

		double Sc3D()	const;									//Returns Sc for the simulation of the 3D field;



		//Functions that return the update coefficient values for the grid nodes
		double CEzE(int m, int n, int p);
		double CEzHx(int m, int n, int p);
		double CEzHy(int m, int n, int p);
		double CHxH(int m, int n, int p);
		double CHxEz(int m, int n, int p);
		double CHxEy(int m, int n, int p);
		double CHyH(int m, int n, int p);
		double CHyEx(int m, int n, int p);
		double CHyEz(int m, int n, int p);
		double CExE(int m, int n, int p);
		double CExHy(int m, int n, int p);
		double CExHz(int m, int n, int p);
		double CEyE(int m, int n, int p);
		double CEyHx(int m, int n, int p);
		double CEyHz(int m, int n, int p);
		double CHzH(int m, int n, int p);
		double CHzEy(int m, int n, int p);
		double CHzEx(int m, int n, int p);

		bool Is_a_unitCell() const;					//Returns true if the structure is the unit cell of a larger periodic structure
		bool Is_periodic_in_yDirection() const;		//Returns true if the larger structure is periodic in the y-direction
		bool Is_periodic_in_zDirection() const;		//Returns true if the larger structure is periodic in the z-direction

		bool find(int m, int n, int p, int& i, int& j, int& k, bool& interior, bool& xLower, bool& yLower, bool& zLower, bool& xBoundary, bool& yBoundary, bool& zBoundary) const;

		//If the Yee cell m,n,p is in the grid, find finds the grid element containing the Yee cell m,n,p and sets indexes i,j, k to the grid element's indexes. Returns true if the element found.
		//If the Yee cell sits on the high x-boundary or high y-boundary or high z-boundary (which is the case if for example the Ex node of the Yee cell is on the high z grid boundary), 
		//find sets the indexes i, j, k and returns true. 
		//If find returns false, then error condition.
		//Also
		//If m is the first node index in the x-layer with layer index i, so that the Yee cell m,n,p is on the x = const lower boundary of the layer with index i, find sets xLower to true. 
		//If n is the first node index in the y-layer with layer index j, so that the Yee cell m,n,p is on the y = const lower boundary of the layer with index j, find sets yLower to true. 
		//If p is the first node index in the z-layer with layer index k, so that the Yee cell m,n,p is on the z = const lower boundary of the layer with index k, find sets zLower to true.
		//If m is the last x-node index  in the grid, so that the Yee cell sits on the high x-boundary, find sets xBoundary to true.
		//If n is the last y-node index  in the grid, so that the Yee cell sits on the high y-boundary, find sets yBoundary to true.
		//If p is the last z-node index  in the grid, so that the Yee cell sits on the high z-boundary, find sets zBoundary to true.

	private:

		Structure3D theStructure;
		Grid3D theGrid;

		int M;		//Number of Ez nodes on a line across grid in the x-direction
		int N;		//Number of Ez nodes on a line across grid in the y-direction
		int P;		//Number of Ex nodes on a line across grid in the z-direction

		int nLx;	//Number of grid layers in the x-direction
		int nLy;	//Number of grid layers in the y-direction
		int nLz;	//Number of gridlayers in the z-direction



		const double Sc = 1 / sqrt(3.0);	//Choose delt = Sc*delMin/c0 so that local Courant number will not exceed 1/sqrt(3.0).
		double delt;

		
		vector <Array3D> CExE_for;			//A vector CExE_for has four elements: the elemnt with index Interior_Ex = 0 is the array of interior coefficients etc. 
		vector <Array3D> CExHy_for;
		vector <Array3D> CExHz_for;

		vector <Array3D> CEyE_for;
		vector <Array3D> CEyHx_for;
		vector <Array3D> CEyHz_for;

		vector <Array3D> CEzE_for;
		vector <Array3D> CEzHx_for;
		vector <Array3D> CEzHy_for;

		vector <Array3D> CHxH_for;
		vector <Array3D> CHxEy_for;
		vector <Array3D> CHxEz_for;

		vector <Array3D> CHyH_for;
		vector <Array3D> CHyEx_for;
		vector <Array3D> CHyEz_for;

		vector <Array3D> CHzH_for;
		vector <Array3D> CHzEx_for;
		vector <Array3D> CHzEy_for;

		vector <Array2D> CExE_on_high_yBoundary_for;
		vector <Array2D> CExHy_on_high_yBoundary_for;
		vector <Array2D> CExHz_on_high_yBoundary_for;

		vector <Array2D> CExE_on_high_zBoundary_for;
		vector <Array2D> CExHy_on_high_zBoundary_for;
		vector <Array2D> CExHz_on_high_zBoundary_for;

		vector <Array2D> CEyE_on_high_zBoundary_for;
		vector <Array2D> CEyHx_on_high_zBoundary_for;
		vector <Array2D> CEyHz_on_high_zBoundary_for;

		vector <Array2D> CEzE_on_high_yBoundary_for;
		vector <Array2D> CEzHx_on_high_yBoundary_for;
		vector <Array2D> CEzHy_on_high_yBoundary_for;

		Array2D CHzH_on_high_zBoundary_for;
		Array2D CHzEx_on_high_zBoundary_for;
		Array2D CHzEy_on_high_zBoundary_for;

		Array2D CHyH_on_high_yBoundary_for;
		Array2D CHyEx_on_high_yBoundary_for;
		Array2D CHyEz_on_high_yBoundary_for;

		Array2D CHxH_on_high_xBoundary_for;
		Array2D CHxEy_on_high_xBoundary_for;
		Array2D CHxEz_on_high_xBoundary_for;

		Array1D CExE_for_GridEdge;
		Array1D CExHy_for_GridEdge;
		Array1D CExHz_for_GridEdge;



		void computeUpdateCoefficientsForGridElements();
							//Computes the arrays of update coefficients for the nodes in the elements. E.g. CEzE_for[Interior](i, j, k) is an update coefficient for nodes in the interior of element i,j,k.


		double selectCEx(int i, int j, int k,bool found,  bool xLower, bool yLower, bool zLower, bool xBoundary, bool yBoundary, bool zBoundary, 
							vector <Array3D>& CExF_for, vector <Array2D>& CExF_for_high_yGridBoundary, vector <Array2D>& CExF_for_high_zGridBoundary, Array1D& CExF_for_GridEdge);
							//The bool value found indicates wether a given Yee cell is inside one of the grid elements
							//Indexes i,j,k identify the grid element containing the Yee cell if found is true or having the Yee cell on its surface if found is false. 
							//The bool values xlower, yLower, zLower indicate wether 
							//the Yee cell is on one of the low boundaaries of the element.  
							//If found is true and all are false, then the Yee cell is in the interior of the element i,j,k. 
							//Otherwise, with found true, the Yee cell is on one or more of the low boundaries of element i,j,k as indicated by xlower, yLower, zLower.
							//If found is false, then the Yee cell sits outside of the computational domain on a boundary. xBoundary, yBoundary, zBoundary indicate the boundary. 
							//This function selects the corresponding vector of coefficient arrays referenced by the parameters, 
							//either the vector CEx_for which contains the arrays for grid elements inside the computational domain,  
							//or one of the vectors for elements that sit on a computational boundary. It also then selects the proper array in it according to the values of the bools,
							//e.g. the array of CExF values for the interior of grid elements: 
							//The correct CExF coefficient value, e.g., CExF(i,j,k) for element i,j,k will then be returned.

							//If there is no CExF coefficient corresponding with the information in the bools, then an error condition is set.
							
							

		double selectCEy(int i, int j, int k, bool found, bool xLower, bool yLower, bool zLower, bool xBoundary, bool yBoundary, bool zBoundary, 
							vector <Array3D>& CEyF_for, vector <Array2D>& CExF_for_high_zGridBoundary);
							//Selects one the vectors of coefficient arrays referenced by the parameters and the element of that vector containing the CEy coefficient arrray who's value 
							//corresponding with i,j,k and xLower, yLower, zLower and xBoundary, yBoundary, zBoundary is to be returned and returns that value.
							

		double selectCEz(int i, int j, int k, bool found, bool xLower, bool yLower, bool yBoundary, bool zBoundary, vector <Array3D>& CEzF_for, vector <Array2D>& CEzF_for_high_yBoundary_for);
							//Selects one the vectors of coefficient arrays referenced by the parameters and the element of that vector containing the CEz coefficient arrray who's value 
							//corresponding with i,j,k and xLower, yLower, zLower and xBoundary, yBoundary, zBoundary is to be returned and returns that value.

		double selectCHx(int i, int j, int k, bool found, bool xLower, bool yLower, bool zLower, bool xBoundary, bool yBoundary, bool zBoundary, vector <Array3D>& CHxF_for, Array2D& CHxF_for_high_xBoundary_for);
						//Selects  the element of the vector of coefficient arrays containing the CHx coefficient arrray who's value corresponding with i,j,k 
						//and xLower, yLower, zLower and xBoundary, yBoundary, zBoundary is to be returned and returns that value.

		double selectCHy(int i, int j, int k, bool found, bool yLower, bool yBoundary, bool zBoundary, vector <Array3D>& CHyF_for, Array2D& CHyF_for_high_yBoundary_for);
						//Selects  the element of the vector of coefficient arrays containing the CHy coefficient arrray who's value corresponding with i,j,k 
						//and xLower, yLower, zLower and xBoundary, yBoundary, zBoundary is to be returned and returns that value.
						//In the case of uniform magnetic permeability = mu0, no boundary CHy arrays needed.
							

		double selectCHz(int i, int j, int k, bool found, bool xLast, bool yLast, bool zLast, bool yBoundary, bool zBoundary, vector <Array3D>& CHzF_for, Array2D& CHzF_for_high_zGridBoundary);
					//Selects  the element of the vector of coefficient arrays containing the CHz coefficient arrray who's value corresponding with i,j,k 
					//and xLower, yLower, zLower and xBoundary, yBoundary, zBoundary is to be returned and returns that value.
					//In the case of uniform magnetic permeability = mu0, no boundary CHy arrays needed.


	};
	


	//***************************************************************************************************************************************************************************

	class Grid1D
	{
	public:
		Grid1D() {}

		Grid1D(Structure1D& dev, double minLambda);	//Constructor for a Grid based on a 1D device which may have multiple layers in x-direction. 
													//Grid has a layer above the device and a layer below the device. 
													//Grid with N Ez electric field node planes begins on an Ez node at x=0 and 
													//ends on an Hy magnetic field node at (N-1/2)*delx. 
													//Each layer except the top layer begins and ends on a magnetic field node  plane.
													//Thickness in micrometers of layer sub i except top and bottom layers is 
													//thickness of the device layer: dev.th(i - 1) The number of nodes in grid layer in the device is
													//max((int)((dev.th(i - 1) / minLambda)*NLambdaMin + 0.5), minCpt).
													//Number of nodes in layer above device is minCpt. Number of nodes in layer under device is minCpt. 
													//Node spacing for grid layer sub i except top layer is delx(i) = th(i) / nNds(i). (See accessor function declarations below.)
													//Node spacing in layer above device is equal to node spacing in first device layr. 
													//deltx for bottom layer is the minimum deltx.



		Grid1D(int anNodes, double adelx);		//Two-layer one-dimensional Grid with anNodes Ez electric field node planes in top layer and minCpt Ez electric field nodes in bottom layer. 
												//Total number of Ez nodes is then N = anNodes + minCpt.
												//Grid has uniform delx assigned the value adelx
												//Grid begins on an Ez node at x=0 and ends on an Hy magnetic field node at (anNodes - 1)*adlex +  adelx/2 + minCpt*adelx. 
												//Depth is (anNodes - 1)*adlex +  adelx/2 + minCpt*adelx. 
												//The second layer may be used for a lossy layer terminating the grid i.e., Pefectly Matched Layer (PML)).
												//adelx should be in micrometers


		int size() const;						//Returns number of Ez electric field nodes
		int numL() const;						//Returns number of grid layers
		int nNds(int layer) const;				//Returns number of Ez electric field nodes in grid layer with subscript layer
		double delx(int layer) const;			//Returns node spacing in grid layer with subscript layer
		double th(int layer) const;				//Returns thickness of grid layer with subscript layer in micrometers
		double minDelx();						//Returns minimum delx over all of the layers.

	private:

		int N;				//Number of electric field nodes (= number of magnetic field nodes)
		int numLayers;		//Number of layers in the grid (= number of device layers + 2)
		int mD;				//If grid fitted to structure or device, device or structure starts at x-index mD + 1/2. Each layer ends on a magnetic field node.
		double minDeltx;	//Minimum delx over all of the layers

		int* nNodes;		//nNodes[i] is number of nodes in layer i+1 (sub i) of grid.
		double* thickness;  //thickness[i] is thickness in micrometers of layer i+1 of grid.
		double * deltx;		//deltx[i] is deltax in micrometers for layer i+1 of grid.
	};



}//End namespace FDTD


#endif
