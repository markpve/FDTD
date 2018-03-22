#include "grid.h"
//#include "constants.h"

#include <stdexcept>


using namespace std;


namespace FDTD
{
	//**********************************************************************************************************************************************
	//Functions for Structure3DWithGrid

	Structure3DWithGrid::Structure3DWithGrid(Structure3D& str, double minLambda, Report& myReport) : theStructure(str), theGrid(str, minLambda, myReport), nLx(theGrid.numLx()), nLy(theGrid.numLy()), nLz(theGrid.numLz()),

		CExE_for(4, Array3D(nLx, nLy, nLz)), CExHy_for(4, Array3D(nLx, nLy, nLz)), CExHz_for(4, Array3D(nLx, nLy, nLz)),
		CEyE_for(4, Array3D(nLx, nLy, nLz)), CEyHx_for(4, Array3D(nLx, nLy, nLz)), CEyHz_for(4, Array3D(nLx, nLy, nLz)),
		CEzE_for(4, Array3D(nLx, nLy, nLz)), CEzHx_for(4, Array3D(nLx, nLy, nLz)), CEzHy_for(4, Array3D(nLx, nLy, nLz)),

		CHxH_for(2, Array3D(nLx, nLy, nLz)), CHxEy_for(2, Array3D(nLx, nLy, nLz)), CHxEz_for(2, Array3D(nLx, nLy, nLz)),
		CHyH_for(2, Array3D(nLx, nLy, nLz)), CHyEx_for(2, Array3D(nLx, nLy, nLz)), CHyEz_for(2, Array3D(nLx, nLy, nLz)),
		CHzH_for(2, Array3D(nLx, nLy, nLz)), CHzEx_for(2, Array3D(nLx, nLy, nLz)), CHzEy_for(2, Array3D(nLx, nLy, nLz)),


		CExE_on_high_yBoundary_for(2, Array2D(nLx, nLz)), CExHy_on_high_yBoundary_for(2, Array2D(nLx, nLz)), CExHz_on_high_yBoundary_for(2, Array2D(nLx, nLz)),
		CExE_on_high_zBoundary_for(2, Array2D(nLx, nLy)), CExHy_on_high_zBoundary_for(2, Array2D(nLx, nLy)), CExHz_on_high_zBoundary_for(2, Array2D(nLx, nLy)),
		CEyE_on_high_zBoundary_for(2, Array2D(nLx, nLy)), CEyHx_on_high_zBoundary_for(2, Array2D(nLx, nLy)), CEyHz_on_high_zBoundary_for(2, Array2D(nLx, nLy)),
		CEzE_on_high_yBoundary_for(2, Array2D(nLx, nLz)), CEzHx_on_high_yBoundary_for(2, Array2D(nLx, nLz)), CEzHy_on_high_yBoundary_for(2, Array2D(nLx, nLz)),

		CHzH_on_high_zBoundary_for(nLx, nLy), CHzEx_on_high_zBoundary_for(nLx, nLy), CHzEy_on_high_zBoundary_for(nLx, nLy),

		CHxH_on_high_xBoundary_for(nLy, nLz), CHxEy_on_high_xBoundary_for(nLy, nLz), CHxEz_on_high_xBoundary_for(nLy, nLz),

		CHyH_on_high_yBoundary_for(nLx, nLz), CHyEx_on_high_yBoundary_for(nLx, nLz), CHyEz_on_high_yBoundary_for(nLx, nLz),

		CExE_for_GridEdge(nLx), CExHy_for_GridEdge(nLx), CExHz_for_GridEdge(nLx)

	{
		delt = theGrid.minDelta()*Sc / c0;	//delt is in microseconds and delx in micrometers.

		M = theGrid.sizeX();
		N = theGrid.sizeY();
		P = theGrid.sizeZ();

		computeUpdateCoefficientsForGridElements();

	}//End constructor Structure3DWithGrid



	Structure3DWithGrid::Structure3DWithGrid(vector <double>& axLThick, vector <double>& ayLThick, vector <double>& azLThick, vector <vector <vector <RectElement>>>& aelement, double minLambda,
		bool uC, bool yDP, bool zDP, Report& myReport) :
		theStructure(axLThick, ayLThick, azLThick, aelement, uC, yDP, zDP), theGrid(theStructure, minLambda, myReport), nLx(theGrid.numLx()), nLy(theGrid.numLy()), nLz(theGrid.numLz()),



		CExE_for(4, Array3D(nLx, nLy, nLz)), CExHy_for(4, Array3D(nLx, nLy, nLz)), CExHz_for(4, Array3D(nLx, nLy, nLz)),
		CEyE_for(4, Array3D(nLx, nLy, nLz)), CEyHx_for(4, Array3D(nLx, nLy, nLz)), CEyHz_for(4, Array3D(nLx, nLy, nLz)),
		CEzE_for(4, Array3D(nLx, nLy, nLz)), CEzHx_for(4, Array3D(nLx, nLy, nLz)), CEzHy_for(4, Array3D(nLx, nLy, nLz)),

		CHxH_for(2, Array3D(nLx, nLy, nLz)), CHxEy_for(2, Array3D(nLx, nLy, nLz)), CHxEz_for(2, Array3D(nLx, nLy, nLz)),
		CHyH_for(2, Array3D(nLx, nLy, nLz)), CHyEx_for(2, Array3D(nLx, nLy, nLz)), CHyEz_for(2, Array3D(nLx, nLy, nLz)),
		CHzH_for(2, Array3D(nLx, nLy, nLz)), CHzEx_for(2, Array3D(nLx, nLy, nLz)), CHzEy_for(2, Array3D(nLx, nLy, nLz)),


		CExE_on_high_yBoundary_for(2, Array2D(nLx, nLz)), CExHy_on_high_yBoundary_for(2, Array2D(nLx, nLz)), CExHz_on_high_yBoundary_for(2, Array2D(nLx, nLz)),
		CExE_on_high_zBoundary_for(2, Array2D(nLx, nLy)), CExHy_on_high_zBoundary_for(2, Array2D(nLx, nLy)), CExHz_on_high_zBoundary_for(2, Array2D(nLx, nLy)),
		CEyE_on_high_zBoundary_for(2, Array2D(nLx, nLy)), CEyHx_on_high_zBoundary_for(2, Array2D(nLx, nLy)), CEyHz_on_high_zBoundary_for(2, Array2D(nLx, nLy)),
		CEzE_on_high_yBoundary_for(2, Array2D(nLx, nLz)), CEzHx_on_high_yBoundary_for(2, Array2D(nLx, nLz)), CEzHy_on_high_yBoundary_for(2, Array2D(nLx, nLz)),

		CHzH_on_high_zBoundary_for(nLx, nLy), CHzEx_on_high_zBoundary_for(nLx, nLy), CHzEy_on_high_zBoundary_for(nLx, nLy),

		CHxH_on_high_xBoundary_for(nLy, nLz), CHxEy_on_high_xBoundary_for(nLy, nLz), CHxEz_on_high_xBoundary_for(nLy, nLz),

		CHyH_on_high_yBoundary_for(nLx, nLz), CHyEx_on_high_yBoundary_for(nLx, nLz), CHyEz_on_high_yBoundary_for(nLx, nLz),

		CExE_for_GridEdge(nLx), CExHy_for_GridEdge(nLx), CExHz_for_GridEdge(nLx)

	{
		delt = theGrid.minDelta()*Sc / c0;	//delt is in microseconds and delx in micrometers.

		M = theGrid.sizeX();
		N = theGrid.sizeY();
		P = theGrid.sizeZ();

		computeUpdateCoefficientsForGridElements();

	}//End constructor Structure3DWithGrid	



	 //********************************************************************************************************************************	
	 //******************************************************************************************************************************
	void  Structure3DWithGrid::computeUpdateCoefficientsForGridElements()
	{
		int i;
		int j;
		int k;

		double loss;
		double permit;
		double cond;
		double permea;

		//double condEdge;
		//double permitEdge;
		//double permeaEdge;
		//double lossEdge;

		double cond_xyEdge;
		double permit_xyEdge;
		double permea_xyEdge;
		double loss_xyEdge;

		double cond_xzEdge;
		double permit_xzEdge;
		double permea_xzEdge;
		double loss_xzEdge;

		double cond_yzEdge;
		double permit_yzEdge;
		double permea_yzEdge;
		double loss_yzEdge;

		double deltaxLInterface;
		double deltayLInterface;
		double deltazLInterface;


		
		//***************************************************************************************************************************
		//Compute coeficients for grid layer above structure

		i = 0;

			for (j = 0; j < nLy; j++)		//j is the subscript of the grid y-layer
			{
				for (k = 0; k < nLz; k++)
				{
					//Compute properties for elements above structure
					permit = eps0;				//For vacuum
					cond = 0;
					permea = mu0;

					loss = cond*delt*1E-6 / (2 * permit);

					CEzE_for[Interior_Ez](i, j, k) = (1 - loss) / (1 + loss);								
					CEzHx_for[Interior_Ez](i, j, k) = delt / (permit*theGrid.deltayL(j)) / (1 + loss);
					CEzHy_for[Interior_Ez](i, j, k) = delt / (permit*theGrid.deltaxL(i)) / (1 + loss);

					CHxH_for[Interior_Hx](i, j, k) = 1;
					CHxEy_for[Interior_Hx](i, j, k) = delt / (permea*theGrid.deltazL(k));		//CHx coefficients not needed for n = N - 1
					CHxEz_for[Interior_Hx](i, j, k) = delt / (permea*theGrid.deltayL(j));

					CHyH_for[Interior_Hy](i, j, k) = 1;											//CHy coefficients not needed for m = M - 1.
					CHyEx_for[Interior_Hy](i, j, k) = delt / (permea*theGrid.deltazL(k));
					CHyEz_for[Interior_Hy](i, j, k) = delt / (permea*theGrid.deltaxL(i));

					CExE_for[Interior_Ex](i, j, k) = (1 - loss) / (1 + loss);								//CEx coefficents not needed for m = M-1, n = 0, n =  N-1, p = P-1
					CExHz_for[Interior_Ex](i, j, k) = delt / (permit*theGrid.deltayL(j)) / (1 + loss);
					CExHy_for[Interior_Ex](i, j, k) = delt / (permit*theGrid.deltazL(k)) / (1 + loss);

					CEyE_for[Interior_Ey](i, j, k) = (1 - loss) / (1 + loss);
					CEyHx_for[Interior_Ey](i, j, k) = delt / (permit*theGrid.deltazL(k)) / (1 + loss);		//CEy coefficients not needed for m = 0, m = M-1, n = N-1, p = p-1
					CEyHz_for[Interior_Ey](i, j, k) = delt / (permit*theGrid.deltaxL(i)) / (1 + loss);

					CHzH_for[Interior_Hz](i, j, k) = 1;
					CHzEx_for[Interior_Hz](i, j, k) = delt / (permea*theGrid.deltayL(j));		//CHz coefficients not needed for m = M-1, n = N-1
					CHzEy_for[Interior_Hz](i, j, k) = delt / (permea*theGrid.deltaxL(i));


					//Coefficients for x-interface, xz-edge and xy-edge nodes on element i,j,k
					//Assume uniform magnetic permeability

					//For element above the structure
					deltaxLInterface = theGrid.deltaxL(i);

					CHxH_for[xInterface_Hx](i, j, k) = 1;
					CHxEz_for[xInterface_Hx](i, j, k) = CHxEz_for[Interior_Hx](i, j, k);
					CHxEy_for[xInterface_Hx](i, j, k) = CHxEy_for[Interior_Hx](i, j, k);



					 //Compute coefficients for y-interface nodes

					if (j > 0)
					{
						deltayLInterface = theGrid.deltayL(j) / 2 + theGrid.deltayL(j - 1) / 2;
							
						permit_xyEdge = eps0;
						cond_xyEdge = 0;
						permea_xyEdge = mu0;
						loss_xyEdge = cond_xyEdge*delt*1E-6 / (2 * permit_xyEdge);
							

						CEzE_for[yInterface_Ez](i, j, k) = (1 - loss) / (1 + loss);								
						CEzHx_for[yInterface_Ez](i, j, k) = delt / (permit*deltayLInterface) / (1 + loss);
						CEzHy_for[yInterface_Ez](i, j, k) = delt / (permit*theGrid.deltaxL(i)) / (1 + loss);

						CExE_for[yInterface_Ex](i, j, k) = (1 - loss) / (1 + loss);								
						CExHz_for[yInterface_Ex](i, j, k) = delt / (permit*deltayLInterface) / (1 + loss);
						CExHy_for[yInterface_Ex](i, j, k) = delt / (permit*theGrid.deltazL(k)) / (1 + loss);

						CHyH_for[yInterface_Hy](i, j, k) = 1;											//CHy coefficients not needed for m = M - 1.
						CHyEx_for[yInterface_Hy](i, j, k) = delt / (permea*theGrid.deltazL(k));
						CHyEz_for[yInterface_Hy](i, j, k) = delt / (permea*theGrid.deltaxL(i));

						CEzE_for[xyEdge_Ez](i, j, k) = (1 - loss_xyEdge) / (1 + loss_xyEdge);								//
						CEzHx_for[xyEdge_Ez](i, j, k) = delt / (permit_xyEdge*deltayLInterface) / (1 + loss_xyEdge);
						CEzHy_for[xyEdge_Ez](i, j, k) = delt / (permit_xyEdge*deltaxLInterface) / (1 + loss_xyEdge);

					}
					else
					{
						permit_xyEdge = eps0;
						cond_xyEdge = 0;
						permea_xyEdge = mu0;
						loss_xyEdge = cond_xyEdge*delt*1E-6 / (2 * permit_xyEdge);

						deltayLInterface = theGrid.deltayL(j) / 2 + theGrid.deltayL(nLy - 1) / 2;	//Assumes y - periodicity

						CEzE_for[yInterface_Ez](i, j, k) = (1 - loss) / (1 + loss);								//
						CEzHx_for[yInterface_Ez](i, j, k) = delt / (permit*deltayLInterface) / (1 + loss);
						CEzHy_for[yInterface_Ez](i, j, k) = delt / (permit*theGrid.deltaxL(i)) / (1 + loss);

						CExE_for[yInterface_Ex](i, j, k) = (1 - loss) / (1 + loss);								//CEx coefficents not needed for m = M-1, n = 0, n =  N-1, p = P-1
						CExHz_for[yInterface_Ex](i, j, k) = delt / (permit*deltayLInterface) / (1 + loss);
						CExHy_for[yInterface_Ex](i, j, k) = delt / (permit*theGrid.deltazL(k)) / (1 + loss);

						CHyH_for[yInterface_Hy](i, j, k) = 1;											//CHy coefficients not needed for m = M - 1.
						CHyEx_for[yInterface_Hy](i, j, k) = delt / (permea*theGrid.deltazL(k));
						CHyEz_for[yInterface_Hy](i, j, k) = delt / (permea*theGrid.deltaxL(i));

						CEzE_for[xyEdge_Ez](i, j, k) = (1 - loss_xyEdge) / (1 + loss_xyEdge);								//
						CEzHx_for[xyEdge_Ez](i, j, k) = delt / (permit_xyEdge*deltayLInterface) / (1 + loss_xyEdge);
						CEzHy_for[xyEdge_Ez](i, j, k) = delt / (permit_xyEdge*deltaxLInterface) / (1 + loss_xyEdge);

					}//End if j > 0 else


					 //Compute coefficients for z-interface nodes

					if (k > 0)
					{
						deltazLInterface = theGrid.deltazL(k) / 2 + theGrid.deltazL(k - 1) / 2;

						if (j > 0)
							deltayLInterface = theGrid.deltayL(j) / 2 + theGrid.deltayL(j - 1) / 2;
						else
							deltayLInterface = theGrid.deltayL(j) / 2 + theGrid.deltayL(nLy - 1) / 2;	//Assume y-direction periodicity
						//End if

						permit_xzEdge = eps0;
						cond_xzEdge = 0;
						permea_xzEdge = mu0;
						loss_xzEdge = cond_xzEdge*delt*1E-6 / (2 * permit_xzEdge);

						cond_yzEdge = 0;
						permit_yzEdge = eps0;
						permea_yzEdge = mu0;
						loss_yzEdge = cond_yzEdge*delt*1E-6 / (2 * permit_yzEdge);
							

						CExE_for[zInterface_Ex](i, j, k) = (1 - loss) / (1 + loss);								//CEx coefficents not needed for m = M-1, n = 0, n =  N-1, p = P-1
						CExHz_for[zInterface_Ex](i, j, k) = delt / (permit*theGrid.deltayL(j)) / (1 + loss);
						CExHy_for[zInterface_Ex](i, j, k) = delt / (permit*deltazLInterface) / (1 + loss);

						CEyE_for[zInterface_Ey](i, j, k) = (1 - loss) / (1 + loss);
						CEyHx_for[zInterface_Ey](i, j, k) = delt / (permit*deltazLInterface) / (1 + loss);		
						CEyHz_for[zInterface_Ey](i, j, k) = delt / (permit*theGrid.deltaxL(i)) / (1 + loss);

						CHzH_for[zInterface_Hz](i, j, k) = 1;
						CHzEx_for[zInterface_Hz](i, j, k) = delt / (permea*theGrid.deltayL(j));		//CHz coefficients not needed for m = M-1, n = N-1
						CHzEy_for[zInterface_Hz](i, j, k) = delt / (permea*theGrid.deltaxL(i));

						CEyE_for[xzEdge_Ey](i, j, k) = (1 - loss_xzEdge) / (1 + loss_xzEdge);
						CEyHx_for[xzEdge_Ey](i, j, k) = delt / (permit_xzEdge*deltazLInterface) / (1 + loss_xzEdge);		//CEy coefficients not needed for m = 0, m = M-1, n = N-1, p = p-1
						CEyHz_for[xzEdge_Ey](i, j, k) = delt / (permit_xzEdge*deltaxLInterface) / (1 + loss_xzEdge);

						CExE_for[yzEdge_Ex](i, j, k) = (1 - loss_yzEdge) / (1 + loss_yzEdge);								//CEx coefficents not needed for m = M-1, n = 0, n =  N-1, p = P-1
						CExHz_for[yzEdge_Ex](i, j, k) = delt / (permit_yzEdge*deltayLInterface) / (1 + loss_yzEdge);
						CExHy_for[yzEdge_Ex](i, j, k) = delt / (permit_yzEdge*deltazLInterface) / (1 + loss_yzEdge);

					}
					else
						//Calculations for k = 0
					{
						deltazLInterface = theGrid.deltazL(k) / 2 + theGrid.deltazL(nLz - 1) / 2;		//Assumes z periodicity

						if (j > 0)
							deltayLInterface = theGrid.deltayL(j) / 2 + theGrid.deltayL(j - 1) / 2;
						else
							//Assume  y direction periodicity
							deltayLInterface = theGrid.deltayL(j) / 2 + theGrid.deltayL(nLy - 1) / 2;
						//End if 

						permit_xzEdge = eps0;
						cond_xzEdge = 0;
						permea_xzEdge = mu0;
						loss_xzEdge = cond_xzEdge*delt*1E-6 / (2 * permit_xzEdge);

						cond_yzEdge = 0;
						permit_yzEdge = eps0;
						permea_yzEdge = mu0;
						loss_yzEdge = cond_yzEdge*delt*1E-6 / (2 * permit_yzEdge);

						CExE_for[zInterface_Ex](i, j, k) = (1 - loss) / (1 + loss);								//CEx coefficents not needed for m = M-1, n = 0, n =  N-1, p = P-1
						CExHz_for[zInterface_Ex](i, j, k) = delt / (permit*theGrid.deltayL(j)) / (1 + loss);
						CExHy_for[zInterface_Ex](i, j, k) = delt / (permit*deltazLInterface) / (1 + loss);

						CEyE_for[zInterface_Ey](i, j, k) = (1 - loss) / (1 + loss);
						CEyHx_for[zInterface_Ey](i, j, k) = delt / (permit*deltazLInterface) / (1 + loss);		//CEy coefficients not needed for m = 0, m = M-1, n = N-1, p = p-1
						CEyHz_for[zInterface_Ey](i, j, k) = delt / (permit*theGrid.deltaxL(i)) / (1 + loss);

						CHzH_for[zInterface_Hz](i, j, k) = 1;
						CHzEx_for[zInterface_Hz](i, j, k) = delt / (permea*theGrid.deltayL(j));		//CHz coefficients not needed for m = M-1, n = N-1
						CHzEy_for[zInterface_Hz](i, j, k) = delt / (permea*theGrid.deltaxL(i));

						CEyE_for[xzEdge_Ey](i, j, k) = (1 - loss_xzEdge) / (1 + loss_xzEdge);
						CEyHx_for[xzEdge_Ey](i, j, k) = delt / (permit_xzEdge*deltazLInterface) / (1 + loss_xzEdge);		//CEy coefficients not needed for m = 0, m = M-1, n = N-1, p = p-1
						CEyHz_for[xzEdge_Ey](i, j, k) = delt / (permit_xzEdge*deltaxLInterface) / (1 + loss_xzEdge);

						CExE_for[yzEdge_Ex](i, j, k) = (1 - loss_yzEdge) / (1 + loss_yzEdge);								//CEx coefficents not needed for m = M-1, n = 0, n =  N-1, p = P-1
						CExHz_for[yzEdge_Ex](i, j, k) = delt / (permit_yzEdge*deltayLInterface) / (1 + loss_yzEdge);
						CExHy_for[yzEdge_Ex](i, j, k) = delt / (permit_yzEdge*deltazLInterface) / (1 + loss_yzEdge);

					}//End if k > 0 else


				}//End for each k
			}//End for each j


			//***************************************************************************************
			//Compute coefficients for high y and z-boundaries of grid above structure


			permit = eps0;				//For vacuum above structure
			cond = 0;
			permea = mu0;

			loss = cond*delt*1E-6 / (2 * permit);


			//Compute coefficients for high y-boundary of grid

			deltayLInterface = theGrid.deltayL(nLy - 1) / 2 + theGrid.deltayL(0) / 2;

			for (k = 0; k < nLz; k++)
			{
				if(k > 0)
					deltazLInterface = theGrid.deltazL(k) / 2 + theGrid.deltazL(k - 1) / 2;		
				else
					deltazLInterface = theGrid.deltazL(k) / 2 + theGrid.deltazL(nLz - 1) / 2;		//Assumes z periodicity
				//End if

				CEzE_on_high_yBoundary_for[Interior_of_yBoundary_element_Ez](i, k) = (1 - loss) / (1 + loss);
				CEzHx_on_high_yBoundary_for[Interior_of_yBoundary_element_Ez](i, k) = delt / (permit*deltayLInterface) / (1 + loss);
				CEzHy_on_high_yBoundary_for[Interior_of_yBoundary_element_Ez](i, k) = delt / (permit*theGrid.deltaxL(i)) / (1 + loss);

				CEzE_on_high_yBoundary_for[xInterface_on_yBoundary_element_Ez](i, k) = (1 - loss) / (1 + loss);
				CEzHx_on_high_yBoundary_for[xInterface_on_yBoundary_element_Ez](i, k) = delt / (permit*deltayLInterface) / (1 + loss);
				CEzHy_on_high_yBoundary_for[xInterface_on_yBoundary_element_Ez](i, k) = delt / (permit*theGrid.deltaxL(i)) / (1 + loss);

				CExE_on_high_yBoundary_for[Interior_of_yBoundary_element_Ex](i, k) = (1 - loss) / (1 + loss);
				CExHz_on_high_yBoundary_for[Interior_of_yBoundary_element_Ex](i, k) = delt / (permit*deltayLInterface) / (1 + loss);
				CExHy_on_high_yBoundary_for[Interior_of_yBoundary_element_Ex](i, k) = delt / (permit*theGrid.deltazL(k)) / (1 + loss);

				CExE_on_high_yBoundary_for[zInterface_on_yBoundary_element_Ex](i, k) = (1 - loss) / (1 + loss);
				CExHz_on_high_yBoundary_for[zInterface_on_yBoundary_element_Ex](i, k) = delt / (permit*deltayLInterface) / (1 + loss);
				CExHy_on_high_yBoundary_for[zInterface_on_yBoundary_element_Ex](i, k) = delt / (permit*deltazLInterface) / (1 + loss);

				CHyH_on_high_yBoundary_for(i, k) = 1;
				CHyEx_on_high_yBoundary_for(i, k) = delt / (permea*theGrid.deltazL(k));
				CHyEz_on_high_yBoundary_for(i, k) = delt / (permea*theGrid.deltaxL(i));

			}//End for each k


			//Compute coefficients for high z-boundary of grid

			deltazLInterface = theGrid.deltazL(nLz - 1) / 2 + theGrid.deltazL(0) / 2;

			

			for (j = 0; j < nLy; j++)
			{
				if (j > 0)
					deltayLInterface = theGrid.deltayL(j) / 2 + theGrid.deltayL(j - 1) / 2;
				else
					deltayLInterface = theGrid.deltayL(j) / 2 + theGrid.deltayL(nLy - 1) / 2;		//Assumes y-periodicity
																									//End if

				CExE_on_high_zBoundary_for[Interior_of_zBoundary_element_Ex](i, j) = (1 - loss) / (1 + loss);
				CExHz_on_high_zBoundary_for[Interior_of_zBoundary_element_Ex](i, j) = delt / (permit*theGrid.deltayL(j)) / (1 + loss);
				CExHy_on_high_zBoundary_for[Interior_of_zBoundary_element_Ex](i, j) = delt / (permit*deltazLInterface) / (1 + loss);

				CExE_on_high_zBoundary_for[yInterface_on_zBoundary_element_Ex](i, j) = (1 - loss) / (1 + loss);
				CExHz_on_high_zBoundary_for[yInterface_on_zBoundary_element_Ex](i, j) = delt / (permit*deltayLInterface) / (1 + loss);
				CExHy_on_high_zBoundary_for[yInterface_on_zBoundary_element_Ex](i, j) = delt / (permit*deltazLInterface) / (1 + loss);

				CEyE_on_high_zBoundary_for[Interior_of_zBoundary_element_Ey](i, j) = (1 - loss) / (1 + loss);
				CEyHx_on_high_zBoundary_for[Interior_of_zBoundary_element_Ey](i, j) = delt / (permit*deltazLInterface) / (1 + loss);
				CEyHz_on_high_zBoundary_for[Interior_of_zBoundary_element_Ey](i, j) = delt / (permit*theGrid.deltaxL(i)) / (1 + loss);

				CEyE_on_high_zBoundary_for[xInterface_on_zBoundary_element_Ey](i, j) = (1 - loss) / (1 + loss);
				CEyHx_on_high_zBoundary_for[xInterface_on_zBoundary_element_Ey](i, j) = delt / (permit*deltazLInterface) / (1 + loss);
				CEyHz_on_high_zBoundary_for[xInterface_on_zBoundary_element_Ey](i, j) = delt / (permit*theGrid.deltaxL(i)) / (1 + loss);

				CHzH_on_high_zBoundary_for(i, j) = 1;
				CHzEx_on_high_zBoundary_for(i, j) = delt / (permea*theGrid.deltayL(j));
				CHzEy_on_high_zBoundary_for(i, j) = delt / (permea*theGrid.deltaxL(i));

			}//End for each j

			CExE_for_GridEdge(i) = (1 - loss) / (1 + loss);
			CExHy_for_GridEdge(i) = delt / (permit*deltazLInterface) / (1 + loss);
			CExHz_for_GridEdge(i) = delt / (permit*deltayLInterface) / (1 + loss);

		

		//End calculations for elements above structure
		//***********************************************************************************************************


		//***************************************************************************************************************************
		//Compute coefficients for grid elements in structure 

		for (i = 1; i < nLx - 1; i++)
		{
			for (j = 0; j < nLy; j++)		//j is the subscript of the grid y-layer
			{
				for (k = 0; k < nLz; k++)
				{
					
					cond = theStructure.sigma(i - 1, j, k);
					permit = theStructure.eps(i - 1, j, k);
					permea = theStructure.mu(i - 1, j, k);
					loss = cond*delt*1E-6 / (2 * permit);

					CEzE_for[Interior_Ez](i, j, k) = (1 - loss) / (1 + loss);								//
					CEzHx_for[Interior_Ez](i, j, k) = delt / (permit*theGrid.deltayL(j)) / (1 + loss);
					CEzHy_for[Interior_Ez](i, j, k) = delt / (permit*theGrid.deltaxL(i)) / (1 + loss);

					CHxH_for[Interior_Hx](i, j, k) = 1;
					CHxEy_for[Interior_Hx](i, j, k) = delt / (permea*theGrid.deltazL(k));		//CHx coefficients not needed for n = N - 1
					CHxEz_for[Interior_Hx](i, j, k) = delt / (permea*theGrid.deltayL(j));

					CHyH_for[Interior_Hy](i, j, k) = 1;											//CHy coefficients not needed for m = M - 1.
					CHyEx_for[Interior_Hy](i, j, k) = delt / (permea*theGrid.deltazL(k));
					CHyEz_for[Interior_Hy](i, j, k) = delt / (permea*theGrid.deltaxL(i));

					CExE_for[Interior_Ex](i, j, k) = (1 - loss) / (1 + loss);								//CEx coefficents not needed for m = M-1, n = 0, n =  N-1, p = P-1
					CExHz_for[Interior_Ex](i, j, k) = delt / (permit*theGrid.deltayL(j)) / (1 + loss);
					CExHy_for[Interior_Ex](i, j, k) = delt / (permit*theGrid.deltazL(k)) / (1 + loss);

					CEyE_for[Interior_Ey](i, j, k) = (1 - loss) / (1 + loss);
					CEyHx_for[Interior_Ey](i, j, k) = delt / (permit*theGrid.deltazL(k)) / (1 + loss);		//CEy coefficients not needed for m = 0, m = M-1, n = N-1, p = p-1
					CEyHz_for[Interior_Ey](i, j, k) = delt / (permit*theGrid.deltaxL(i)) / (1 + loss);

					CHzH_for[Interior_Hz](i, j, k) = 1;
					CHzEx_for[Interior_Hz](i, j, k) = delt / (permea*theGrid.deltayL(j));		//CHz coefficients not needed for m = M-1, n = N-1
					CHzEy_for[Interior_Hz](i, j, k) = delt / (permea*theGrid.deltaxL(i));


					//Coefficients for x-interface, xz-edge and xy-edge nodes on element i,j,k
					//Assume uniform magnetic permeability

					deltaxLInterface = theGrid.deltaxL(i) / 2 + theGrid.deltaxL(i - 1) / 2;

					//Compute properties for element in structure
					if (i > 1 && i < nLx - 1)
					{
						cond = theStructure.sigma(i - 1, j, k) / 2 + theStructure.sigma(i - 2, j, k) / 2;
						permit = theStructure.eps(i - 1, j, k) / 2 + theStructure.eps(i - 2, j, k) / 2;
						permea = theStructure.mu(i - 1, j, k) / 2 + theStructure.mu(i - 2, j, k) / 2;
						loss = cond*delt*1E-6 / (2 * permit);
					}
					else
						if (i == 1)
						{
							cond = theStructure.sigma(i - 1, j, k) / 2;								//Conductivity is zero in grid layer above the theStructure
							permit = theStructure.eps(i - 1, j, k) / 2 + eps0 / 2;
							permea = theStructure.mu(i - 1, j, k) / 2 + mu0 / 2;
							loss = cond*delt*1E-6 / (2 * permit);
						}//End if
					//End if
							


					CHxH_for[xInterface_Hx](i, j, k) = 1;
					CHxEz_for[xInterface_Hx](i, j, k) = CHxEz_for[Interior_Hx](i, j, k);
					CHxEy_for[xInterface_Hx](i, j, k) = CHxEy_for[Interior_Hx](i, j, k);

					CEzE_for[xInterface_Ez](i, j, k) = (1 - loss) / (1 + loss);								//
					CEzHx_for[xInterface_Ez](i, j, k) = delt / (permit*theGrid.deltayL(j)) / (1 + loss);
					CEzHy_for[xInterface_Ez](i, j, k) = delt / (permit*deltaxLInterface) / (1 + loss);

					CEyE_for[xInterface_Ey](i, j, k) = (1 - loss) / (1 + loss);
					CEyHx_for[xInterface_Ey](i, j, k) = delt / (permit*theGrid.deltazL(k)) / (1 + loss);		//CEy coefficients not needed for m = 0, m = M-1, n = N-1, p = p-1
					CEyHz_for[xInterface_Ey](i, j, k) = delt / (permit*theGrid.deltaxL(i)) / (1 + loss);


					 //Compute coefficients for y-interface nodes

					if (j > 0)
					{
						deltayLInterface = theGrid.deltayL(j) / 2 + theGrid.deltayL(j - 1) / 2;

						//Compute properties for element in structure
							cond = theStructure.sigma(i - 1, j, k) / 2 + theStructure.sigma(i - 1, j - 1, k) / 2;
							permit = theStructure.eps(i - 1, j, k) / 2 + theStructure.eps(i - 1, j - 1, k) / 2;
							permea = theStructure.mu(i - 1, j, k) / 2 + theStructure.mu(i - 1, j - 1, k) / 2;
							loss = cond*delt*1E-6 / (2 * permit);

							//Needed for x - y edge values
							deltaxLInterface = theGrid.deltaxL(i) / 2 + theGrid.deltaxL(i - 1) / 2;

							if (i > 1)
							{
								cond_xyEdge = theStructure.sigma(i - 1, j, k) / 4 + theStructure.sigma(i - 1, j - 1, k) / 4 + theStructure.sigma(i - 2, j, k) / 4 + theStructure.sigma(i - 2, j - 1, k) / 4;
								permit_xyEdge = theStructure.eps(i - 1, j, k) / 4 + theStructure.eps(i - 1, j - 1, k) / 4 + theStructure.eps(i - 2, j, k) / 4 + theStructure.eps(i - 2, j - 1, k) / 4;
								permea_xyEdge = theStructure.mu(i - 1, j, k) / 4 + theStructure.mu(i - 1, j - 1, k) / 4 + theStructure.mu(i - 2, j, k) / 4 + theStructure.mu(i - 2, j - 1, k) / 4;
								loss_xyEdge = cond_xyEdge*delt*1E-6 / (2 * permit_xyEdge);
							}
							else
							{
								cond_xyEdge = theStructure.sigma(i - 1, j, k) / 4 + theStructure.sigma(i - 1, j - 1, k) / 4;
								permit_xyEdge = theStructure.eps(i - 1, j, k) / 4 + theStructure.eps(i - 1, j - 1, k) / 4 + eps0 / 4 + eps0 / 4;
								permea_xyEdge = theStructure.mu(i - 1, j, k) / 4 + theStructure.mu(i - 1, j - 1, k) / 4 + mu0 / 4 + mu0 / 4;
								loss_xyEdge = cond_xyEdge*delt*1E-6 / (2 * permit_xyEdge);
							}//end if
						

						CEzE_for[yInterface_Ez](i, j, k) = (1 - loss) / (1 + loss);								//
						CEzHx_for[yInterface_Ez](i, j, k) = delt / (permit*deltayLInterface) / (1 + loss);
						CEzHy_for[yInterface_Ez](i, j, k) = delt / (permit*theGrid.deltaxL(i)) / (1 + loss);

						CExE_for[yInterface_Ex](i, j, k) = (1 - loss) / (1 + loss);								//CEx coefficents not needed for m = M-1, n = 0, n =  N-1, p = P-1
						CExHz_for[yInterface_Ex](i, j, k) = delt / (permit*deltayLInterface) / (1 + loss);
						CExHy_for[yInterface_Ex](i, j, k) = delt / (permit*theGrid.deltazL(k)) / (1 + loss);

						CHyH_for[yInterface_Hy](i, j, k) = 1;											//CHy coefficients not needed for m = M - 1.
						CHyEx_for[yInterface_Hy](i, j, k) = delt / (permea*theGrid.deltazL(k));
						CHyEz_for[yInterface_Hy](i, j, k) = delt / (permea*theGrid.deltaxL(i));

						CEzE_for[xyEdge_Ez](i, j, k) = (1 - loss_xyEdge) / (1 + loss_xyEdge);								//
						CEzHx_for[xyEdge_Ez](i, j, k) = delt / (permit_xyEdge*deltayLInterface) / (1 + loss_xyEdge);
						CEzHy_for[xyEdge_Ez](i, j, k) = delt / (permit_xyEdge*deltaxLInterface) / (1 + loss_xyEdge);

					}
					else
					{
						deltayLInterface = theGrid.deltayL(j) / 2 + theGrid.deltayL(nLy - 1) / 2;	//Assumes y - periodicity

						//Compute properties for element in structure
						cond = theStructure.sigma(i - 1, j, k) / 2 + theStructure.sigma(i - 1, nLy - 1, k) / 2;		//Assumes y-periodicity
						permit = theStructure.eps(i - 1, j, k) / 2 + theStructure.eps(i - 1, nLy - 1, k) / 2;			//Assumes y-periodicity
						loss = cond*delt*1E-6 / (2 * permit);

						if(i > 1)
						{
							//Assume y-periodicity
							cond_xyEdge = theStructure.sigma(i - 1, j, k) / 4 + theStructure.sigma(i - 1, nLy - 1, k) / 4 + theStructure.sigma(i - 2, j, k) / 4 + theStructure.sigma(i - 2, nLy - 1, k) / 4;
							permit_xyEdge = theStructure.eps(i - 1, j, k) / 4 + theStructure.eps(i - 1, nLy - 1, k) / 4 + theStructure.eps(i - 2, j, k) / 4 + theStructure.eps(i - 2, nLy - 1, k) / 4;
							permea_xyEdge = theStructure.mu(i - 1, j, k) / 4 + theStructure.mu(i - 1, nLy - 1, k) / 4 + theStructure.mu(i - 2, j, k) / 4 + theStructure.mu(i - 2, nLy - 1, k) / 4;
							loss_xyEdge = cond_xyEdge*delt*1E-6 / (2 * permit_xyEdge);
						}
						else
						{
							cond_xyEdge = theStructure.sigma(i - 1, j, k) / 4 + theStructure.sigma(i - 1, nLy - 1, k) / 4;
							permit_xyEdge = theStructure.eps(i - 1, j, k) / 4 + theStructure.eps(i - 1, nLy - 1, k) / 4 + eps0 / 4 + eps0 / 4;
							permea_xyEdge = theStructure.mu(i - 1, j, k) / 4 + theStructure.mu(i - 1, nLy - 1, k) / 4 + mu0 / 4 + mu0 / 4;
							loss_xyEdge = cond_xyEdge*delt*1E-6 / (2 * permit_xyEdge);
						}//end if
						

						CHyH_for[yInterface_Hy](i, j, k) = 1;											//CHy coefficients not needed for m = M - 1.
						CHyEx_for[yInterface_Hy](i, j, k) = delt / (permea*theGrid.deltazL(k));
						CHyEz_for[yInterface_Hy](i, j, k) = delt / (permea*theGrid.deltaxL(i));

						CEzE_for[yInterface_Ez](i, j, k) = (1 - loss) / (1 + loss);
						CEzHy_for[yInterface_Ez](i, j, k) = delt / (permit*theGrid.deltaxL(i)) / (1 + loss);
						CEzHx_for[yInterface_Ez](i, j, k) = delt / (permit*deltayLInterface) / (1 + loss);

						CEzE_for[xyEdge_Ez](i, j, k) = (1 - loss_xyEdge) / (1 + loss_xyEdge);								//
						CEzHx_for[xyEdge_Ez](i, j, k) = delt / (permit_xyEdge*deltayLInterface) / (1 + loss_xyEdge);
						CEzHy_for[xyEdge_Ez](i, j, k) = delt / (permit_xyEdge*deltaxLInterface) / (1 + loss_xyEdge);

						CExE_for[yInterface_Ex](i, j, k) = (1 - loss) / (1 + loss);								//CEx coefficents not needed for m = M-1, n = 0, n =  N-1, p = P-1
						CExHz_for[yInterface_Ex](i, j, k) = delt / (permit*deltayLInterface) / (1 + loss);
						CExHy_for[yInterface_Ex](i, j, k) = delt / (permit*theGrid.deltazL(k)) / (1 + loss);

						
					}//End if j > 0 else



					 //Compute coefficients for z-interface nodes

					if (k > 0)
					{
						deltazLInterface = theGrid.deltazL(k) / 2 + theGrid.deltazL(k - 1) / 2;
	
						//Compute properties for low z interface nodes not on an edge *******************************
						cond = theStructure.sigma(i - 1, j, k) / 2 + theStructure.sigma(i - 1, j, k - 1) / 2;
						permit = theStructure.eps(i - 1, j, k) / 2 + theStructure.eps(i - 1, j, k - 1) / 2;
						permea = theStructure.mu(i - 1, j, k) / 2 + theStructure.mu(i - 1, j, k - 1) / 2;
						loss = cond*delt*1E-6 / (2 * permit);

						//Compute properties for edge nodes *********************************************************

						//Compute properties for yz edge nodes
						if (j > 0)
						{
							deltayLInterface = theGrid.deltayL(j) / 2 + theGrid.deltayL(j - 1) / 2;

							cond_yzEdge = theStructure.sigma(i - 1, j, k) / 4 + theStructure.sigma(i - 1, j, k - 1) / 4 + theStructure.sigma(i - 1, j - 1, k) / 4 + theStructure.sigma(i - 1, j - 1, k - 1) / 4;
							permit_yzEdge = theStructure.eps(i - 1, j, k) / 4 + theStructure.eps(i - 1, j, k - 1) / 4 + theStructure.eps(i - 1, j - 1, k) / 4 + theStructure.eps(i - 1, j - 1, k - 1) / 4;
							permea_yzEdge = theStructure.mu(i - 1, j, k) / 4 + theStructure.mu(i - 1, j, k - 1) / 4 + theStructure.mu(i - 1, j - 1, k) / 4 + theStructure.mu(i - 1, j - 1, k - 1) / 4;
							loss_yzEdge = cond_yzEdge*delt*1E-6 / (2 * permit_yzEdge);
						}
						else
						{
							//Assume  y direction periodicity
							deltayLInterface = theGrid.deltayL(j) / 2 + theGrid.deltayL(nLy - 1) / 2;

							cond_yzEdge = theStructure.sigma(i - 1, j, k) / 4 + theStructure.sigma(i - 1, j, k - 1) / 4 + theStructure.sigma(i - 1, nLy - 1, k) / 4 + theStructure.sigma(i - 1, nLy - 1, k - 1) / 4;
							permit_yzEdge = theStructure.eps(i - 1, j, k) / 4 + theStructure.eps(i - 1, j, k - 1) / 4 + theStructure.eps(i - 1, nLy - 1, k) / 4 + theStructure.eps(i - 1, nLy - 1, k - 1) / 4;
							permea_yzEdge = theStructure.mu(i - 1, j, k) / 4 + theStructure.mu(i - 1, j, k - 1) / 4 + theStructure.mu(i - 1, nLy - 1, k) / 4 + theStructure.mu(i - 1, nLy - 1, k - 1) / 4;
							loss_yzEdge = cond_yzEdge*delt*1E-6 / (2 * permit_yzEdge);
						}//End if

						//End computation of yz edge properties

						//Compute properties for xz edge nodes

						if (i > 1)
						{
							cond_xzEdge = theStructure.sigma(i - 1, j, k) / 4 + theStructure.sigma(i - 1, j, k - 1) / 4 + theStructure.sigma(i - 2, j, k) / 4 + theStructure.sigma(i - 2, j, k - 1) / 4;
							permit_xzEdge = theStructure.eps(i - 1, j, k) / 4 + theStructure.eps(i - 1, j, k - 1) / 4 + theStructure.eps(i - 2, j, k) / 4 + theStructure.eps(i - 2, j, k - 1) / 4;
							permea_xzEdge = theStructure.mu(i - 1, j, k) / 4 + theStructure.mu(i - 1, j, k - 1) / 4 + theStructure.mu(i - 2, j, k) / 4 + theStructure.mu(i - 2, j, k - 1) / 4;
							loss_xzEdge = cond_xzEdge*delt*1E-6 / (2 * permit_xzEdge);
						}
						else
						{
							cond_xzEdge = theStructure.sigma(i - 1, j, k) / 4 + theStructure.sigma(i - 1, j, k - 1) / 4;						//Conduvtivity zero above structure
							permit_xzEdge = theStructure.eps(i - 1, j, k) / 4 + theStructure.eps(i - 1, j, k - 1) / 4 + eps0 / 4 + eps0 / 4;
							permea_xzEdge = theStructure.mu(i - 1, j, k) / 4 + theStructure.mu(i - 1, j, k - 1) / 4 + mu0 / 4 + mu0 / 4;
							loss_xzEdge = cond_xzEdge*delt*1E-6 / (2 * permit_xzEdge);
						}//end if

						CExE_for[zInterface_Ex](i, j, k) = (1 - loss) / (1 + loss);								//CEx coefficents not needed for m = M-1, n = 0, n =  N-1, p = P-1
						CExHz_for[zInterface_Ex](i, j, k) = delt / (permit*theGrid.deltayL(j)) / (1 + loss);
						CExHy_for[zInterface_Ex](i, j, k) = delt / (permit*deltazLInterface) / (1 + loss);

						CEyE_for[zInterface_Ey](i, j, k) = (1 - loss) / (1 + loss);
						CEyHx_for[zInterface_Ey](i, j, k) = delt / (permit*deltazLInterface) / (1 + loss);		//CEy coefficients not needed for m = 0, m = M-1, n = N-1, p = p-1
						CEyHz_for[zInterface_Ey](i, j, k) = delt / (permit*theGrid.deltaxL(i)) / (1 + loss);

						CHzH_for[zInterface_Hz](i, j, k) = 1;
						CHzEx_for[zInterface_Hz](i, j, k) = delt / (permea*theGrid.deltayL(j));		//CHz coefficients not needed for m = M-1, n = N-1
						CHzEy_for[zInterface_Hz](i, j, k) = delt / (permea*theGrid.deltaxL(i));

						CEyE_for[xzEdge_Ey](i, j, k) = (1 - loss_xzEdge) / (1 + loss_xzEdge);
						CEyHx_for[xzEdge_Ey](i, j, k) = delt / (permit_xzEdge*deltazLInterface) / (1 + loss_xzEdge);		//CEy coefficients not needed for m = 0, m = M-1, n = N-1, p = p-1
						CEyHz_for[xzEdge_Ey](i, j, k) = delt / (permit_xzEdge*deltaxLInterface) / (1 + loss_xzEdge);

						CExE_for[yzEdge_Ex](i, j, k) = (1 - loss_yzEdge) / (1 + loss_yzEdge);								//CEx coefficents not needed for m = M-1, n = 0, n =  N-1, p = P-1
						CExHz_for[yzEdge_Ex](i, j, k) = delt / (permit_yzEdge*deltayLInterface) / (1 + loss_yzEdge);
						CExHy_for[yzEdge_Ex](i, j, k) = delt / (permit_yzEdge*deltazLInterface) / (1 + loss_yzEdge);

					}
					else
						//Calculations for k = 0
					{
						deltazLInterface = theGrid.deltazL(k) / 2 + theGrid.deltazL(nLz - 1) / 2;		//Assumes z periodicity

						//Compute properties for element in structure

							cond = theStructure.sigma(i - 1, j, k) / 2 + theStructure.sigma(i - 1, nLy - 1, k) / 2;		//Assumes y-periodicity
							permit = theStructure.eps(i - 1, j, k) / 2 + theStructure.eps(i - 1, nLy - 1, k) / 2;			//Assumes y-periodicity
							loss = cond*delt*1E-6 / (2 * permit);

							//Compute properties for edge nodes
				
							//Compute properties for yz edge nodes
							if (j > 0)
							{
								deltayLInterface = theGrid.deltayL(j) / 2 + theGrid.deltayL(j - 1) / 2;

								//Assume z -direction periodicity
								cond_yzEdge = theStructure.sigma(i - 1, j, k) / 4 + theStructure.sigma(i - 1, j - 1, k) / 4 + theStructure.sigma(i - 1, j, nLz - 1) / 4 + theStructure.sigma(i - 1, j - 1, nLz - 1) / 4;
								permit_yzEdge = theStructure.eps(i - 1, j, k) / 4 + theStructure.eps(i - 1, j - 1, k) / 4 + theStructure.eps(i - 1, j, nLz - 1) / 4 + theStructure.eps(i - 1, j - 1, nLz - 1) / 4;
								permea_yzEdge = theStructure.mu(i - 1, j, k) / 4 + theStructure.mu(i - 1, j - 1, k) / 4 + theStructure.mu(i - 1, j, nLz - 1) / 4 + theStructure.mu(i - 1, j - 1, nLz - 1) / 4;
								loss_yzEdge = cond_yzEdge*delt*1E-6 / (2 * permit_yzEdge);
							}
							else
							{
								//Assume  y and z direction periodicity

								deltayLInterface = theGrid.deltayL(j) / 2 + theGrid.deltayL(nLy - 1) / 2;

								cond_yzEdge = theStructure.sigma(i - 1, j, k) / 4 + theStructure.sigma(i - 1, j, nLz - 1) / 4 + theStructure.sigma(i - 1, nLy - 1, k) / 4 + theStructure.sigma(i - 1, nLy - 1, nLz - 1) / 4;
								permit_yzEdge = theStructure.eps(i - 1, j, k) / 4 + theStructure.eps(i - 1, j, nLz - 1) / 4 + theStructure.eps(i - 1, nLy - 1, k) / 4 + theStructure.eps(i - 1, nLy - 1, nLz - 1) / 4;
								permea_yzEdge = theStructure.mu(i - 1, j, k) / 4 + theStructure.mu(i - 1, j, nLz - 1) / 4 + theStructure.mu(i - 1, nLy - 1, k) / 4 + theStructure.mu(i - 1, nLy - 1, nLz - 1) / 4;
								loss_yzEdge = cond_yzEdge*delt*1E-6 / (2 * permit_yzEdge);


							}//End if

							 //End computation of yz edge properties

							 //Compute xz edge properties
							if (i > 1)
							{
								//Assume z-periodicity
								cond_xzEdge = theStructure.sigma(i - 1, j, k) / 4 + theStructure.sigma(i - 2, j, k) / 4 + theStructure.sigma(i - 2, j, nLz - 1) / 4 + theStructure.sigma(i - 1, j, nLz - 1) / 4;
								permit_xzEdge = theStructure.eps(i - 1, j, k) / 4 + theStructure.eps(i - 2, j, k) / 4 + theStructure.eps(i - 2, j, nLz - 1) / 4 + theStructure.eps(i - 1, j, nLz - 1) / 4;
								permea_xzEdge = theStructure.mu(i - 1, j, k) / 4 + theStructure.mu(i - 2, j, k) / 4 + theStructure.mu(i - 2, j, nLz - 1) / 4 + theStructure.mu(i - 1, j, nLz - 1) / 4;
								loss_xzEdge = cond_xzEdge*delt*1E-6 / (2 * permit_xzEdge);
							}
							else
							{
								//Assume z-periodicity
								cond_xzEdge = theStructure.sigma(i - 1, j, k) / 4 + theStructure.sigma(i - 1, j, nLz - 1) / 4;
								permit_xzEdge = theStructure.eps(i - 1, j, k) / 4 + theStructure.eps(i - 1, j, nLz - 1) / 4 + eps0/4 + eps0/4;
								permea_xzEdge = theStructure.mu(i - 1, j, k) / 4 +  theStructure.mu(i - 1, j, nLz - 1) / 4 + mu0 / 4 + mu0 / 4;
								loss_xzEdge = cond_xzEdge*delt*1E-6 / (2 * permit_xzEdge);
							}//end if
						

						CExE_for[zInterface_Ex](i, j, k) = (1 - loss) / (1 + loss);								//CEx coefficents not needed for m = M-1, n = 0, n =  N-1, p = P-1
						CExHz_for[zInterface_Ex](i, j, k) = delt / (permit*theGrid.deltayL(j)) / (1 + loss);
						CExHy_for[zInterface_Ex](i, j, k) = delt / (permit*deltazLInterface) / (1 + loss);

						CEyE_for[zInterface_Ey](i, j, k) = (1 - loss) / (1 + loss);
						CEyHx_for[zInterface_Ey](i, j, k) = delt / (permit*deltazLInterface) / (1 + loss);		//CEy coefficients not needed for m = 0, m = M-1, n = N-1, p = p-1
						CEyHz_for[zInterface_Ey](i, j, k) = delt / (permit*theGrid.deltaxL(i)) / (1 + loss);

						CHzH_for[zInterface_Hz](i, j, k) = 1;
						CHzEx_for[zInterface_Hz](i, j, k) = delt / (permea*theGrid.deltayL(j));		//CHz coefficients not needed for m = M-1, n = N-1
						CHzEy_for[zInterface_Hz](i, j, k) = delt / (permea*theGrid.deltaxL(i));

						CEyE_for[xzEdge_Ey](i, j, k) = (1 - loss_xzEdge) / (1 + loss_xzEdge);
						CEyHx_for[xzEdge_Ey](i, j, k) = delt / (permit_xzEdge*deltazLInterface) / (1 + loss_xzEdge);		//CEy coefficients not needed for m = 0, m = M-1, n = N-1, p = p-1
						CEyHz_for[xzEdge_Ey](i, j, k) = delt / (permit_xzEdge*deltaxLInterface) / (1 + loss_xzEdge);

						CExE_for[yzEdge_Ex](i, j, k) = (1 - loss_yzEdge) / (1 + loss_yzEdge);								//CEx coefficents not needed for m = M-1, n = 0, n =  N-1, p = P-1
						CExHz_for[yzEdge_Ex](i, j, k) = delt / (permit_yzEdge*deltayLInterface) / (1 + loss_yzEdge);
						CExHy_for[yzEdge_Ex](i, j, k) = delt / (permit_yzEdge*deltazLInterface) / (1 + loss_yzEdge);

					}//End if k > 0 else

					 
				}//End for each k
			}//End for each j


			 //***************************************************************************************
			 //Compute coefficients for high y and z-boundaries of grid


			//Compute coefficients for high y-boundary of grid

			deltayLInterface = theGrid.deltayL(nLy - 1) / 2 + theGrid.deltayL(0) / 2;
			deltaxLInterface = theGrid.deltaxL(i) / 2 + theGrid.deltaxL(i - 1) / 2;

			for (k = 0; k < nLz; k++)
			{
				//Properties for nodes in the interior of y-boundary surface elements
				cond = theStructure.sigma(i - 1, nLy - 1, k) / 2 + theStructure.sigma(i - 1, 0, k) / 2;			//Assumes y-periodicity
				permit = theStructure.eps(i - 1, nLy - 1, k) / 2 + theStructure.eps(i - 1, 0, k) / 2;
				permea = theStructure.mu(i - 1, nLy - 1, k) / 2 + theStructure.mu(i - 1, 0, k) / 2;
				loss = cond*delt*1E-6 / (2 * permit);

				//For y-z edge nodes
				if (k > 0)
				{
					deltazLInterface = theGrid.deltazL(k) / 2 + theGrid.deltazL(k - 1) / 2;

					//properties for y-boundary z-interface nodes
					//Assume y-direction periodicity
					cond_yzEdge = theStructure.sigma(i - 1, nLy - 1, k) / 4 + theStructure.sigma(i - 1, nLy - 1, k - 1) / 4 + theStructure.sigma(i - 1, 0, k) / 4 + theStructure.sigma(i - 1, 0, k - 1) / 4;
					permit_yzEdge = theStructure.eps(i - 1, nLy - 1, k) / 4 + theStructure.eps(i - 1, nLy - 1, k - 1) / 4 + theStructure.eps(i - 1, 0, k) / 4 + theStructure.eps(i - 1, 0, k - 1) / 4;
					permea_yzEdge = theStructure.mu(i - 1, nLy - 1, k) / 4 + theStructure.mu(i - 1, nLy - 1, k - 1) / 4 + theStructure.mu(i - 1, 0, k) / 4 + theStructure.mu(i - 1, 0, k - 1) / 4;
					loss_yzEdge = cond_yzEdge*delt*1E-6 / (2 * permit_yzEdge);
				}
					
				else
				{
					deltazLInterface = theGrid.deltazL(k) / 2 + theGrid.deltazL(nLz - 1) / 2;		//Assumes z-periodicity
					
					//properties for y-boundary z-interface nodes
					//Assume y-direction and z-direction periodicity
					cond_yzEdge = theStructure.sigma(i - 1, nLy - 1, k) / 4 + theStructure.sigma(i - 1, nLy - 1, 0) / 4 + theStructure.sigma(i - 1, 0, k) / 4 + theStructure.sigma(i - 1, 0, 0) / 4;
					permit_yzEdge = theStructure.eps(i - 1, nLy - 1, k) / 4 + theStructure.eps(i - 1, nLy - 1, 0) / 4 + theStructure.eps(i - 1, 0, k) / 4 + theStructure.eps(i - 1, 0, 0) / 4;
					permea_yzEdge = theStructure.mu(i - 1, nLy - 1, k) / 4 + theStructure.mu(i - 1, nLy - 1, 0) / 4 + theStructure.mu(i - 1, 0, k) / 4 + theStructure.mu(i - 1, 0, 0) / 4;
					loss_yzEdge = cond_yzEdge*delt*1E-6 / (2 * permit_yzEdge);
				}	
				//End if


				//For x - y edge nodes

				if (i > 1)
				{
					//properties for y-boundary x-interface nodes
					//Assume y-direction periodicity
					cond_xyEdge = theStructure.sigma(i - 1, nLy - 1, k) / 4 + theStructure.sigma(i - 1, 0, k) / 4 + theStructure.sigma(i - 2, nLy - 1, k) / 4 + theStructure.sigma(i - 2, 0, k) / 4;
					permit_xyEdge = theStructure.eps(i - 1, nLy - 1, k) / 4 + theStructure.eps(i - 1, 0, k) / 4 + theStructure.eps(i - 2, nLy - 1, k) / 4 + theStructure.eps(i - 2, 0, k) / 4;
					permea_xyEdge = theStructure.mu(i - 1, nLy - 1, k) / 4 + theStructure.mu(i - 1, 0, k) / 4 + theStructure.mu(i - 2, nLy - 1, k) / 4 + theStructure.mu(i - 2, 0, k) / 4;
					loss_xyEdge = cond_xyEdge*delt*1E-6 / (2 * permit_xyEdge);
				}
				else
				{
					//properties for y-boundary x-interface nodes
					//Assume y-direction periodicity
					cond_xyEdge = theStructure.sigma(i - 1, nLy - 1, k) / 4 + theStructure.sigma(i - 1, 0, k) / 4;
					permit_xyEdge = theStructure.eps(i - 1, nLy - 1, k) / 4 + theStructure.eps(i - 1, 0, k) / 4 + eps0 / 4 + eps0 / 4;
					permea_xyEdge = theStructure.mu(i - 1, nLy - 1, k) / 4 + theStructure.mu(i - 1, 0, k) / 4 + mu0 / 4 + mu0 / 4;
					loss_xyEdge = cond_xyEdge*delt*1E-6 / (2 * permit_xyEdge);
				}//end if


				CEzE_on_high_yBoundary_for[Interior_of_yBoundary_element_Ez](i, k) = (1 - loss) / (1 + loss);
				CEzHx_on_high_yBoundary_for[Interior_of_yBoundary_element_Ez](i, k) = delt / (permit*deltayLInterface) / (1 + loss);
				CEzHy_on_high_yBoundary_for[Interior_of_yBoundary_element_Ez](i, k) = delt / (permit*theGrid.deltaxL(i)) / (1 + loss);

				CEzE_on_high_yBoundary_for[xInterface_on_yBoundary_element_Ez](i, k) = (1 - loss_xyEdge) / (1 + loss_xyEdge);
				CEzHx_on_high_yBoundary_for[xInterface_on_yBoundary_element_Ez](i, k) = delt / (permit_xyEdge*deltayLInterface) / (1 + loss_xyEdge);
				CEzHy_on_high_yBoundary_for[xInterface_on_yBoundary_element_Ez](i, k) = delt / (permit_xyEdge*deltaxLInterface) / (1 + loss_xyEdge);

				CExE_on_high_yBoundary_for[Interior_of_yBoundary_element_Ex](i, k) = (1 - loss) / (1 + loss);
				CExHz_on_high_yBoundary_for[Interior_of_yBoundary_element_Ex](i, k) = delt / (permit*deltayLInterface) / (1 + loss);
				CExHy_on_high_yBoundary_for[Interior_of_yBoundary_element_Ex](i, k) = delt / (permit*theGrid.deltazL(k)) / (1 + loss);

				CExE_on_high_yBoundary_for[zInterface_on_yBoundary_element_Ex](i, k) = (1 - loss_yzEdge) / (1 + loss_yzEdge);
				CExHz_on_high_yBoundary_for[zInterface_on_yBoundary_element_Ex](i, k) = delt / (permit_yzEdge*deltayLInterface) / (1 + loss_yzEdge);
				CExHy_on_high_yBoundary_for[zInterface_on_yBoundary_element_Ex](i, k) = delt / (permit_yzEdge*deltazLInterface) / (1 + loss_yzEdge);

				CHyH_on_high_yBoundary_for(i, k) = 1;
				CHyEx_on_high_yBoundary_for(i, k) = delt / (permea*theGrid.deltazL(k));
				CHyEz_on_high_yBoundary_for(i, k) = delt / (permea*theGrid.deltaxL(i));

			}//End for each k


			 //Compute coefficients for high z-boundary of grid

			deltazLInterface = theGrid.deltazL(nLz - 1) / 2 + theGrid.deltazL(0) / 2;

			for (j = 0; j < nLy; j++)
			{
				//Properties for nodes in the interior of z-boundary surface elements
				cond = theStructure.sigma(i - 1, j, nLz - 1) / 2 + theStructure.sigma(i - 1, j, 0) / 2;			//Assumes z-periodicity
				permit = theStructure.eps(i - 1, j, nLz - 1) / 2 + theStructure.eps(i - 1, j, 0) / 2;
				permea = theStructure.mu(i - 1, j, nLz - 1) / 2 + theStructure.mu(i - 1, j, 0) / 2;
				loss = cond*delt*1E-6 / (2 * permit);


			if (j > 0)
			{
				deltayLInterface = theGrid.deltayL(j) / 2 + theGrid.deltayL(j - 1) / 2;

				//properties for y-boundary z-interface nodes
				//Assume z-direction periodicity
				cond_yzEdge = theStructure.sigma(i - 1, j, nLz - 1) / 4 + theStructure.sigma(i - 1, j, 0) / 4 + theStructure.sigma(i - 1, j - 1, nLz - 1) / 4 + theStructure.sigma(i - 1, j - 1, 0) / 4;
				permit_yzEdge = theStructure.eps(i - 1, j, nLz - 1) / 4 + theStructure.eps(i - 1, j, 0) / 4 + theStructure.eps(i - 1, j - 1, nLz - 1) / 4 + theStructure.eps(i - 1, j - 1, 0) / 4;
				permea_yzEdge = theStructure.mu(i - 1, j, nLz - 1) / 4 + theStructure.mu(i - 1, j, 0) / 4 + theStructure.mu(i - 1, j - 1, nLz - 1) / 4 + theStructure.mu(i - 1, j - 1, 0) / 4;
				loss_yzEdge = cond_yzEdge*delt*1E-6 / (2 * permit_yzEdge);
			}
				
			else
			{
				deltayLInterface = theGrid.deltayL(j) / 2 + theGrid.deltayL(nLy - 1) / 2;		//Assumes y-periodicity

				//properties for y-boundary z-interface nodes
				//Assume y and z-direction periodicity
				cond_yzEdge = theStructure.sigma(i - 1, j, nLz - 1) / 4 + theStructure.sigma(i - 1, j, 0) / 4 + theStructure.sigma(i - 1, nLy - 1, nLz - 1) / 4 + theStructure.sigma(i - 1, nLy - 1, 0) / 4;
				permit_yzEdge = theStructure.eps(i - 1, j, nLz - 1) / 4 + theStructure.eps(i - 1, j, 0) / 4 + theStructure.eps(i - 1, nLy - 1, nLz - 1) / 4 + theStructure.eps(i - 1, nLy - 1, 0) / 4;
				permea_yzEdge = theStructure.mu(i - 1, j, nLz - 1) / 4 + theStructure.mu(i - 1, j, 0) / 4 + theStructure.mu(i - 1, nLy - 1, nLz - 1) / 4 + theStructure.mu(i - 1, nLy - 1, 0) / 4;
				loss_yzEdge = cond_yzEdge*delt*1E-6 / (2 * permit_yzEdge);

			}
			//End if

			//Compute xz edge properties
			if (i > 1)
			{
				//Assume z-periodicity
				cond_xzEdge = theStructure.sigma(i - 1, j, nLz - 1) / 4 + theStructure.sigma(i - 2, j, nLz - 1) / 4 + theStructure.sigma(i - 2, j, 0) / 4 + theStructure.sigma(i - 1, j, 0) / 4;
				permit_xzEdge = theStructure.eps(i - 1, j, nLz - 1) / 4 + theStructure.eps(i - 2, j, nLz - 1) / 4 + theStructure.eps(i - 2, j, 0) / 4 + theStructure.eps(i - 1, j, 0) / 4;
				permea_xzEdge = theStructure.mu(i - 1, j, nLz - 1) / 4 + theStructure.mu(i - 2, j, nLz - 1) / 4 + theStructure.mu(i - 2, j, 0) / 4 + theStructure.mu(i - 1, j, 0) / 4;
				loss_xzEdge = cond_xzEdge*delt*1E-6 / (2 * permit_xzEdge);
			}
			else
			{
				//Assume z-periodicity
				cond_xzEdge = theStructure.sigma(i - 1, j, nLz - 1) / 4 + theStructure.sigma(i - 1, j, 0) / 4;
				permit_xzEdge = theStructure.eps(i - 1, j, nLz - 1) / 4 + theStructure.eps(i - 1, j, 0) / 4 + eps0 / 4 + eps0 / 4;
				permea_xzEdge = theStructure.mu(i - 1, j, nLz - 1) / 4 + theStructure.mu(i - 1, j, 0) / 4 + mu0 / 4 + mu0 / 4;
				loss_xzEdge = cond_xzEdge*delt*1E-6 / (2 * permit_xzEdge);
			}//end if

			
				CExE_on_high_zBoundary_for[Interior_of_zBoundary_element_Ex](i, j) = (1 - loss) / (1 + loss);
				CExHz_on_high_zBoundary_for[Interior_of_zBoundary_element_Ex](i, j) = delt / (permit*theGrid.deltayL(j)) / (1 + loss);
				CExHy_on_high_zBoundary_for[Interior_of_zBoundary_element_Ex](i, j) = delt / (permit*deltazLInterface) / (1 + loss);

				CExE_on_high_zBoundary_for[yInterface_on_zBoundary_element_Ex](i, j) = (1 - loss_yzEdge) / (1 + loss_yzEdge);
				CExHz_on_high_zBoundary_for[yInterface_on_zBoundary_element_Ex](i, j) = delt / (permit_yzEdge*deltayLInterface) / (1 + loss_yzEdge);
				CExHy_on_high_zBoundary_for[yInterface_on_zBoundary_element_Ex](i, j) = delt / (permit_yzEdge*deltazLInterface) / (1 + loss_yzEdge);

				CEyE_on_high_zBoundary_for[Interior_of_zBoundary_element_Ey](i, j) = (1 - loss) / (1 + loss);
				CEyHx_on_high_zBoundary_for[Interior_of_zBoundary_element_Ey](i, j) = delt / (permit*deltazLInterface) / (1 + loss);
				CEyHz_on_high_zBoundary_for[Interior_of_zBoundary_element_Ey](i, j) = delt / (permit*theGrid.deltaxL(i)) / (1 + loss);

				CEyE_on_high_zBoundary_for[xInterface_on_zBoundary_element_Ey](i, j) = (1 - loss_xzEdge) / (1 + loss_xzEdge);
				CEyHx_on_high_zBoundary_for[xInterface_on_zBoundary_element_Ey](i, j) = delt / (permit_xzEdge*deltazLInterface) / (1 + loss_xzEdge);
				CEyHz_on_high_zBoundary_for[xInterface_on_zBoundary_element_Ey](i, j) = delt / (permit_xzEdge*deltaxLInterface) / (1 + loss_xzEdge);

				CHzH_on_high_zBoundary_for(i, j) = 1;
				CHzEx_on_high_zBoundary_for(i, j) = delt / (permea*theGrid.deltayL(j));
				CHzEy_on_high_zBoundary_for(i, j) = delt / (permea*theGrid.deltaxL(i));

			}//End for each j


			//properties on edge of grid at high y-boundary and high z-boundary
			deltayLInterface = theGrid.deltayL(nLy - 1) / 2 + theGrid.deltayL(0) / 2;

			//properties for nodes
			//Assume y and z-direction periodicity
			cond_yzEdge = theStructure.sigma(i - 1, nLy - 1, nLz - 1) / 4 + theStructure.sigma(i - 1, nLy - 1, 0) / 4 + theStructure.sigma(i - 1, 0, nLz - 1) / 4 + theStructure.sigma(i - 1, 0, 0) / 4;
			permit_yzEdge = theStructure.eps(i - 1, nLy - 1, nLz - 1) / 4 + theStructure.eps(i - 1, nLy - 1, 0) / 4 + theStructure.eps(i - 1, 0, nLz - 1) / 4 + theStructure.eps(i - 1, 0, 0) / 4;
			permea_yzEdge = theStructure.mu(i - 1, nLy - 1, nLz - 1) / 4 + theStructure.mu(i - 1, nLy - 1, 0) / 4 + theStructure.mu(i - 1, 0, nLz - 1) / 4 + theStructure.mu(i - 1, 0, 0) / 4;
			loss_yzEdge = cond_yzEdge*delt*1E-6 / (2 * permit_yzEdge);

			CExE_for_GridEdge(i) = (1 - loss_yzEdge) / (1 + loss_yzEdge);
			CExHy_for_GridEdge(i) = delt / (permit_yzEdge*deltazLInterface) / (1 + loss_yzEdge);
			CExHz_for_GridEdge(i) = delt / (permit_yzEdge*deltayLInterface) / (1 + loss_yzEdge);


		}//End for each i

	//End calculations for elements in structure


	//**********************************************************************************
	//Calculations for elements below structure

		i = nLx - 1;
		
		for (j = 0; j < nLy; j++)		//j is the subscript of the grid y-layer
		{
			for (k = 0; k < nLz; k++)
			{
				
				permit = eps0;				//For vacuum
				cond = 0;
				permea = mu0;

				loss = cond*delt*1E-6 / (2 * permit);

				CEzE_for[Interior_Ez](i, j, k) = (1 - loss) / (1 + loss);								//
				CEzHx_for[Interior_Ez](i, j, k) = delt / (permit*theGrid.deltayL(j)) / (1 + loss);
				CEzHy_for[Interior_Ez](i, j, k) = delt / (permit*theGrid.deltaxL(i)) / (1 + loss);

				CHxH_for[Interior_Hx](i, j, k) = 1;
				CHxEy_for[Interior_Hx](i, j, k) = delt / (permea*theGrid.deltazL(k));		//CHx coefficients not needed for n = N - 1
				CHxEz_for[Interior_Hx](i, j, k) = delt / (permea*theGrid.deltayL(j));

				CHyH_for[Interior_Hy](i, j, k) = 1;											//CHy coefficients not needed for m = M - 1.
				CHyEx_for[Interior_Hy](i, j, k) = delt / (permea*theGrid.deltazL(k));
				CHyEz_for[Interior_Hy](i, j, k) = delt / (permea*theGrid.deltaxL(i));

				CExE_for[Interior_Ex](i, j, k) = (1 - loss) / (1 + loss);								//CEx coefficents not needed for m = M-1, n = 0, n =  N-1, p = P-1
				CExHz_for[Interior_Ex](i, j, k) = delt / (permit*theGrid.deltayL(j)) / (1 + loss);
				CExHy_for[Interior_Ex](i, j, k) = delt / (permit*theGrid.deltazL(k)) / (1 + loss);

				CEyE_for[Interior_Ey](i, j, k) = (1 - loss) / (1 + loss);
				CEyHx_for[Interior_Ey](i, j, k) = delt / (permit*theGrid.deltazL(k)) / (1 + loss);		//CEy coefficients not needed for m = 0, m = M-1, n = N-1, p = p-1
				CEyHz_for[Interior_Ey](i, j, k) = delt / (permit*theGrid.deltaxL(i)) / (1 + loss);

				CHzH_for[Interior_Hz](i, j, k) = 1;
				CHzEx_for[Interior_Hz](i, j, k) = delt / (permea*theGrid.deltayL(j));		//CHz coefficients not needed for m = M-1, n = N-1
				CHzEy_for[Interior_Hz](i, j, k) = delt / (permea*theGrid.deltaxL(i));


				//Coefficients for x-interface, xz-edge and xy-edge nodes on element i,j,k
				//Assume uniform magnetic permeability

				deltaxLInterface = theGrid.deltaxL(i) / 2 + theGrid.deltaxL(i - 1) / 2;
			
				cond = theStructure.sigma(i - 2, j, k) / 2;								//Conductivity is zero in grid layer beneath the the structure
				permit = theStructure.eps(i - 2, j, k) / 2 + eps0 / 2;
				permea = theStructure.mu(i - 2, j, k) / 2 + mu0 / 2;
				loss = cond*delt*1E-6 / (2 * permit);
						
				CHxH_for[xInterface_Hx](i, j, k) = 1;
				CHxEz_for[xInterface_Hx](i, j, k) = CHxEz_for[Interior_Hx](i, j, k);
				CHxEy_for[xInterface_Hx](i, j, k) = CHxEy_for[Interior_Hx](i, j, k);

				CEzE_for[xInterface_Ez](i, j, k) = (1 - loss) / (1 + loss);								//
				CEzHx_for[xInterface_Ez](i, j, k) = delt / (permit*theGrid.deltayL(j)) / (1 + loss);
				CEzHy_for[xInterface_Ez](i, j, k) = delt / (permit*deltaxLInterface) / (1 + loss);

				CEyE_for[xInterface_Ey](i, j, k) = (1 - loss) / (1 + loss);
				CEyHx_for[xInterface_Ey](i, j, k) = delt / (permit*theGrid.deltazL(k)) / (1 + loss);		//CEy coefficients not needed for m = 0, m = M-1, n = N-1, p = p-1
				CEyHz_for[xInterface_Ey](i, j, k) = delt / (permit*theGrid.deltaxL(i)) / (1 + loss);


				//Compute coefficients for y-interface nodes

				if (j > 0)
				{
					deltayLInterface = theGrid.deltayL(j) / 2 + theGrid.deltayL(j - 1) / 2;

					permit = eps0;				//For vacuum
					cond = 0;
					permea = mu0;
					loss = cond*delt*1E-6 / (2 * permit);

							
					cond_xyEdge = theStructure.sigma(i - 2, j, k) / 4 + theStructure.sigma(i - 2, j - 1, k) / 4;
					permit_xyEdge = theStructure.eps(i - 2, j, k) / 4 + theStructure.eps(i - 2, j - 1, k) / 4 + eps0 / 4 + eps0 / 4;
					permea_xyEdge = theStructure.mu(i - 2, j, k) / 4 + theStructure.mu(i - 2, j - 1, k) / 4 + mu0 / 4 + mu0 / 4;
					loss_xyEdge = cond_xyEdge*delt*1E-6 / (2 * permit_xyEdge);
						

					CEzE_for[yInterface_Ez](i, j, k) = (1 - loss) / (1 + loss);								//
					CEzHx_for[yInterface_Ez](i, j, k) = delt / (permit*deltayLInterface) / (1 + loss);
					CEzHy_for[yInterface_Ez](i, j, k) = delt / (permit*theGrid.deltaxL(i)) / (1 + loss);

					CExE_for[yInterface_Ex](i, j, k) = (1 - loss) / (1 + loss);								//CEx coefficents not needed for m = M-1, n = 0, n =  N-1, p = P-1
					CExHz_for[yInterface_Ex](i, j, k) = delt / (permit*deltayLInterface) / (1 + loss);
					CExHy_for[yInterface_Ex](i, j, k) = delt / (permit*theGrid.deltazL(k)) / (1 + loss);

					CHyH_for[yInterface_Hy](i, j, k) = 1;											//CHy coefficients not needed for m = M - 1.
					CHyEx_for[yInterface_Hy](i, j, k) = delt / (permea*theGrid.deltazL(k));
					CHyEz_for[yInterface_Hy](i, j, k) = delt / (permea*theGrid.deltaxL(i));

					CEzE_for[xyEdge_Ez](i, j, k) = (1 - loss_xyEdge) / (1 + loss_xyEdge);								//
					CEzHx_for[xyEdge_Ez](i, j, k) = delt / (permit_xyEdge*deltayLInterface) / (1 + loss_xyEdge);
					CEzHy_for[xyEdge_Ez](i, j, k) = delt / (permit_xyEdge*deltaxLInterface) / (1 + loss_xyEdge);

				}
				else
				{
					deltayLInterface = theGrid.deltayL(j) / 2 + theGrid.deltayL(nLy - 1) / 2;	//Assumes y - periodicity

					permit = eps0;				//For vacuum
					cond = 0;
					permea = mu0;
					loss = cond*delt*1E-6 / (2 * permit);

							
					cond_xyEdge = theStructure.sigma(i - 2, j, k) / 4 + theStructure.sigma(i - 2, nLy - 1, k) / 4;
					permit_xyEdge = theStructure.eps(i - 2, j, k) / 4 + theStructure.eps(i - 2, nLy - 1, k) / 4 + eps0 / 4 + eps0 / 4;
					permea_xyEdge = theStructure.mu(i - 2, j, k) / 4 + theStructure.mu(i - 2, nLy - 1, k) / 4 + mu0 / 4 + mu0 / 4;
					loss_xyEdge = cond_xyEdge*delt*1E-6 / (2 * permit_xyEdge);

					CHyH_for[yInterface_Hy](i, j, k) = 1;											//CHy coefficients not needed for m = M - 1.
					CHyEx_for[yInterface_Hy](i, j, k) = delt / (permea*theGrid.deltazL(k));
					CHyEz_for[yInterface_Hy](i, j, k) = delt / (permea*theGrid.deltaxL(i));

					CEzE_for[yInterface_Ez](i, j, k) = (1 - loss) / (1 + loss);
					CEzHy_for[yInterface_Ez](i, j, k) = delt / (permit*theGrid.deltaxL(i)) / (1 + loss);
					CEzHx_for[yInterface_Ez](i, j, k) = delt / (permit*deltayLInterface) / (1 + loss);

					CEzE_for[xyEdge_Ez](i, j, k) = (1 - loss_xyEdge) / (1 + loss_xyEdge);								//
					CEzHx_for[xyEdge_Ez](i, j, k) = delt / (permit_xyEdge*deltayLInterface) / (1 + loss_xyEdge);
					CEzHy_for[xyEdge_Ez](i, j, k) = delt / (permit_xyEdge*deltaxLInterface) / (1 + loss_xyEdge);

					CExE_for[yInterface_Ex](i, j, k) = (1 - loss) / (1 + loss);								//CEx coefficents not needed for m = M-1, n = 0, n =  N-1, p = P-1
					CExHz_for[yInterface_Ex](i, j, k) = delt / (permit*deltayLInterface) / (1 + loss);
					CExHy_for[yInterface_Ex](i, j, k) = delt / (permit*theGrid.deltazL(k)) / (1 + loss);


				}//End if j > 0 else



				//Compute coefficients for z-interface nodes

				if (k > 0)
				{
					deltazLInterface = theGrid.deltazL(k) / 2 + theGrid.deltazL(k - 1) / 2;

					
					//Compute properties for element not in structure
					permit = eps0;				//For vacuum
					cond = 0;
					permea = mu0;
					loss = cond*delt*1E-6 / (2 * permit);

					if (j > 0)
						deltayLInterface = theGrid.deltayL(j) / 2 + theGrid.deltayL(j - 1) / 2;
					else
						deltayLInterface = theGrid.deltayL(j) / 2 + theGrid.deltayL(nLy - 1) / 2;	//Assume y-direction periodicity
					//End if

							
					cond_xzEdge = theStructure.sigma(i - 2, j, k) / 4 + theStructure.sigma(i - 2, j, k - 1) / 4;
					permit_xzEdge = theStructure.eps(i - 2, j, k) / 4 + theStructure.eps(i - 2, j, k - 1) / 4 + eps0 / 4 + eps0 / 4;
					permea_xzEdge = theStructure.mu(i - 2, j, k) / 4 + theStructure.mu(i - 2, j, k - 1) / 4 + mu0 / 4 + mu0 / 4;
					loss_xzEdge = cond_xzEdge*delt*1E-6 / (2 * permit_xzEdge);

					cond_yzEdge = 0;
					permit_yzEdge = eps0;
					permea_yzEdge = mu0;
					loss_yzEdge = cond_yzEdge*delt*1E-6 / (2 * permit_yzEdge);
						

					CExE_for[zInterface_Ex](i, j, k) = (1 - loss) / (1 + loss);								//CEx coefficents not needed for m = M-1, n = 0, n =  N-1, p = P-1
					CExHz_for[zInterface_Ex](i, j, k) = delt / (permit*theGrid.deltayL(j)) / (1 + loss);
					CExHy_for[zInterface_Ex](i, j, k) = delt / (permit*deltazLInterface) / (1 + loss);

					CEyE_for[zInterface_Ey](i, j, k) = (1 - loss) / (1 + loss);
					CEyHx_for[zInterface_Ey](i, j, k) = delt / (permit*deltazLInterface) / (1 + loss);		//CEy coefficients not needed for m = 0, m = M-1, n = N-1, p = p-1
					CEyHz_for[zInterface_Ey](i, j, k) = delt / (permit*theGrid.deltaxL(i)) / (1 + loss);

					CHzH_for[zInterface_Hz](i, j, k) = 1;
					CHzEx_for[zInterface_Hz](i, j, k) = delt / (permea*theGrid.deltayL(j));		//CHz coefficients not needed for m = M-1, n = N-1
					CHzEy_for[zInterface_Hz](i, j, k) = delt / (permea*theGrid.deltaxL(i));

					CEyE_for[xzEdge_Ey](i, j, k) = (1 - loss_xzEdge) / (1 + loss_xzEdge);
					CEyHx_for[xzEdge_Ey](i, j, k) = delt / (permit_xzEdge*deltazLInterface) / (1 + loss_xzEdge);		//CEy coefficients not needed for m = 0, m = M-1, n = N-1, p = p-1
					CEyHz_for[xzEdge_Ey](i, j, k) = delt / (permit_xzEdge*deltaxLInterface) / (1 + loss_xzEdge);

					CExE_for[yzEdge_Ex](i, j, k) = (1 - loss_yzEdge) / (1 + loss_yzEdge);								//CEx coefficents not needed for m = M-1, n = 0, n =  N-1, p = P-1
					CExHz_for[yzEdge_Ex](i, j, k) = delt / (permit_yzEdge*deltayLInterface) / (1 + loss_yzEdge);
					CExHy_for[yzEdge_Ex](i, j, k) = delt / (permit_yzEdge*deltazLInterface) / (1 + loss_yzEdge);

					}
					else
						//Calculations for k = 0
					{
						deltazLInterface = theGrid.deltazL(k) / 2 + theGrid.deltazL(nLz - 1) / 2;		//Assumes z periodicity

						if (j > 0)
							deltayLInterface = theGrid.deltayL(j) / 2 + theGrid.deltayL(j - 1) / 2;
						else
								//Assume  y and z direction periodicity
								deltayLInterface = theGrid.deltayL(j) / 2 + theGrid.deltayL(nLy - 1) / 2;
							//End if 

						permit = eps0;				//For vacuum
						cond = 0;
						permea = mu0;
						loss = cond*delt*1E-6 / (2 * permit);

						cond_xzEdge = theStructure.sigma(i - 2, j, k) / 4 + theStructure.sigma(i - 2, j, nLz - 1) / 4;
						permit_xzEdge = theStructure.eps(i - 2, j, k) / 4 + theStructure.eps(i - 2, j, nLz - 1) / 4 + eps0 / 4 + eps0 / 4;
						permea_xzEdge = theStructure.mu(i - 2, j, k) / 4 + theStructure.mu(i - 2, j, nLz - 1) / 4 + mu0 / 4 + mu0 / 4;
						loss_xzEdge = cond_xzEdge*delt*1E-6 / (2 * permit_xzEdge);

						cond_yzEdge = 0;
						permit_yzEdge = eps0;
						permea_yzEdge = mu0;
						loss_yzEdge = cond_yzEdge*delt*1E-6 / (2 * permit_yzEdge);

						CExE_for[zInterface_Ex](i, j, k) = (1 - loss) / (1 + loss);								//CEx coefficents not needed for m = M-1, n = 0, n =  N-1, p = P-1
						CExHz_for[zInterface_Ex](i, j, k) = delt / (permit*theGrid.deltayL(j)) / (1 + loss);
						CExHy_for[zInterface_Ex](i, j, k) = delt / (permit*deltazLInterface) / (1 + loss);

						CEyE_for[zInterface_Ey](i, j, k) = (1 - loss) / (1 + loss);
						CEyHx_for[zInterface_Ey](i, j, k) = delt / (permit*deltazLInterface) / (1 + loss);		//CEy coefficients not needed for m = 0, m = M-1, n = N-1, p = p-1
						CEyHz_for[zInterface_Ey](i, j, k) = delt / (permit*theGrid.deltaxL(i)) / (1 + loss);

						CHzH_for[zInterface_Hz](i, j, k) = 1;
						CHzEx_for[zInterface_Hz](i, j, k) = delt / (permea*theGrid.deltayL(j));		//CHz coefficients not needed for m = M-1, n = N-1
						CHzEy_for[zInterface_Hz](i, j, k) = delt / (permea*theGrid.deltaxL(i));

						CEyE_for[xzEdge_Ey](i, j, k) = (1 - loss_xzEdge) / (1 + loss_xzEdge);
						CEyHx_for[xzEdge_Ey](i, j, k) = delt / (permit_xzEdge*deltazLInterface) / (1 + loss_xzEdge);		//CEy coefficients not needed for m = 0, m = M-1, n = N-1, p = p-1
						CEyHz_for[xzEdge_Ey](i, j, k) = delt / (permit_xzEdge*deltaxLInterface) / (1 + loss_xzEdge);

						CExE_for[yzEdge_Ex](i, j, k) = (1 - loss_yzEdge) / (1 + loss_yzEdge);								//CEx coefficents not needed for m = M-1, n = 0, n =  N-1, p = P-1
						CExHz_for[yzEdge_Ex](i, j, k) = delt / (permit_yzEdge*deltayLInterface) / (1 + loss_yzEdge);
						CExHy_for[yzEdge_Ex](i, j, k) = delt / (permit_yzEdge*deltazLInterface) / (1 + loss_yzEdge);

					}//End if k > 0 else


				}//End for each k
			}//End for each j


			 //***************************************************************************************
			 //Compute coefficients for high y and z-boundaries of grid below structure


			permit = eps0;				//For vacuum above structure
			cond = 0;
			permea = mu0;

			loss = cond*delt*1E-6 / (2 * permit);


			//Compute coefficients for high y-boundary of grid

			deltayLInterface = theGrid.deltayL(nLy - 1) / 2 + theGrid.deltayL(0) / 2;

			deltaxLInterface = theGrid.deltaxL(i - 1) / 2 + theGrid.deltaxL(i) / 2;

			for (k = 0; k < nLz; k++)
			{
				if (k > 0)
					deltazLInterface = theGrid.deltazL(k) / 2 + theGrid.deltazL(k - 1) / 2;
				else
					deltazLInterface = theGrid.deltazL(k) / 2 + theGrid.deltazL(nLz - 1) / 2;		//Assumes z periodicity
				//End if

				//properties for y-boundary x-interface nodes
				//Assume y-direction periodicity
				cond_xyEdge =  theStructure.sigma(i - 2, nLy - 1, k) / 4 + theStructure.sigma(i - 2, 0, k) / 4;
				permit_xyEdge = theStructure.eps(i - 2, nLy - 1, k) / 4 + theStructure.eps(i - 2, 0, k) / 4 + eps0/4 + eps0/4;
				permea_xyEdge = theStructure.mu(i - 2, nLy - 1, k) / 4 + theStructure.mu(i - 2, 0, k) / 4 + mu0/4 + mu0/4;
				loss_xyEdge = cond_xyEdge*delt*1E-6 / (2 * permit_xyEdge);

				CEzE_on_high_yBoundary_for[Interior_of_yBoundary_element_Ez](i, k) = (1 - loss) / (1 + loss);
				CEzHx_on_high_yBoundary_for[Interior_of_yBoundary_element_Ez](i, k) = delt / (permit*deltayLInterface) / (1 + loss);
				CEzHy_on_high_yBoundary_for[Interior_of_yBoundary_element_Ez](i, k) = delt / (permit*theGrid.deltaxL(i)) / (1 + loss);

				CEzE_on_high_yBoundary_for[xInterface_on_yBoundary_element_Ez](i, k) = (1 - loss_xyEdge) / (1 + loss_xyEdge);
				CEzHx_on_high_yBoundary_for[xInterface_on_yBoundary_element_Ez](i, k) = delt / (permit_xyEdge*deltayLInterface) / (1 + loss_xyEdge);
				CEzHy_on_high_yBoundary_for[xInterface_on_yBoundary_element_Ez](i, k) = delt / (permit_xyEdge*deltaxLInterface) / (1 + loss_xyEdge);

				CExE_on_high_yBoundary_for[Interior_of_yBoundary_element_Ex](i, k) = (1 - loss) / (1 + loss);
				CExHz_on_high_yBoundary_for[Interior_of_yBoundary_element_Ex](i, k) = delt / (permit*deltayLInterface) / (1 + loss);
				CExHy_on_high_yBoundary_for[Interior_of_yBoundary_element_Ex](i, k) = delt / (permit*theGrid.deltazL(k)) / (1 + loss);

				CExE_on_high_yBoundary_for[zInterface_on_yBoundary_element_Ex](i, k) = (1 - loss) / (1 + loss);
				CExHz_on_high_yBoundary_for[zInterface_on_yBoundary_element_Ex](i, k) = delt / (permit*deltayLInterface) / (1 + loss);
				CExHy_on_high_yBoundary_for[zInterface_on_yBoundary_element_Ex](i, k) = delt / (permit*deltazLInterface) / (1 + loss);

				CHyH_on_high_yBoundary_for(i, k) = 1;
				CHyEx_on_high_yBoundary_for(i, k) = delt / (permea*theGrid.deltazL(k));
				CHyEz_on_high_yBoundary_for(i, k) = delt / (permea*theGrid.deltaxL(i));

			}//End for each k


			 //Compute coefficients for high z-boundary of grid

			deltazLInterface = theGrid.deltazL(nLz - 1) / 2 + theGrid.deltazL(0) / 2;


			for (j = 0; j < nLy; j++)
			{
				if (j > 0)
					deltayLInterface = theGrid.deltayL(j) / 2 + theGrid.deltayL(j - 1) / 2;
				else
					deltayLInterface = theGrid.deltayL(j) / 2 + theGrid.deltayL(nLy - 1) / 2;		//Assumes y-periodicity
				//End if

				CExE_on_high_zBoundary_for[Interior_of_zBoundary_element_Ex](i, j) = (1 - loss) / (1 + loss);
				CExHz_on_high_zBoundary_for[Interior_of_zBoundary_element_Ex](i, j) = delt / (permit*theGrid.deltayL(j)) / (1 + loss);
				CExHy_on_high_zBoundary_for[Interior_of_zBoundary_element_Ex](i, j) = delt / (permit*deltazLInterface) / (1 + loss);

				CExE_on_high_zBoundary_for[yInterface_on_zBoundary_element_Ex](i, j) = (1 - loss) / (1 + loss);
				CExHz_on_high_zBoundary_for[yInterface_on_zBoundary_element_Ex](i, j) = delt / (permit*deltayLInterface) / (1 + loss);
				CExHy_on_high_zBoundary_for[yInterface_on_zBoundary_element_Ex](i, j) = delt / (permit*deltazLInterface) / (1 + loss);

				CEyE_on_high_zBoundary_for[Interior_of_zBoundary_element_Ey](i, j) = (1 - loss) / (1 + loss);
				CEyHx_on_high_zBoundary_for[Interior_of_zBoundary_element_Ey](i, j) = delt / (permit*deltazLInterface) / (1 + loss);
				CEyHz_on_high_zBoundary_for[Interior_of_zBoundary_element_Ey](i, j) = delt / (permit*theGrid.deltaxL(i)) / (1 + loss);

				CEyE_on_high_zBoundary_for[xInterface_on_zBoundary_element_Ey](i, j) = (1 - loss) / (1 + loss);
				CEyHx_on_high_zBoundary_for[xInterface_on_zBoundary_element_Ey](i, j) = delt / (permit*deltazLInterface) / (1 + loss);
				CEyHz_on_high_zBoundary_for[xInterface_on_zBoundary_element_Ey](i, j) = delt / (permit*theGrid.deltaxL(i)) / (1 + loss);

				CHzH_on_high_zBoundary_for(i, j) = 1;
				CHzEx_on_high_zBoundary_for(i, j) = delt / (permea*theGrid.deltayL(j));
				CHzEy_on_high_zBoundary_for(i, j) = delt / (permea*theGrid.deltaxL(i));

			}//End for each j

			CExE_for_GridEdge(i) = (1 - loss) / (1 + loss);
			CExHy_for_GridEdge(i) = delt / (permit*deltazLInterface) / (1 + loss);
			CExHz_for_GridEdge(i) = delt / (permit*deltayLInterface) / (1 + loss);


			//Compute CHx coefficients for high x-boundary of grid

			permea = mu0;

			for (j = 0; j < nLy; j++)
				for (k = 0; k < nLz; k++)
				{
					CHxH_on_high_xBoundary_for(j, k) = CHxH_for[Interior_Hx](i, j, k);
					CHxEy_on_high_xBoundary_for(j, k) = CHxEy_for[Interior_Hx](i, j, k);		//CHx coefficients not needed for n = N - 1
					CHxEz_on_high_xBoundary_for(j, k) = CHxEz_for[Interior_Hx](i, j, k);
				}//End for
			

			//End calculations for elements below structure

	}//End function computeUpdateCoefficientsForGridElements

			
	
	//Accessors

	Grid3D* Structure3DWithGrid::grid() { return &theGrid; }				// Returns pointer to the grid member object

	Structure3D* Structure3DWithGrid::structure() { return &theStructure; } //Returns pointer to the structure member object

	int Structure3DWithGrid::numLx() const { return nLx; }					//Returns the number of layers normal to the x-direction.
	int Structure3DWithGrid::numLy() const { return nLy; }					//Returns the number of layers in the y-direction.		
	int Structure3DWithGrid::numLz() const { return nLz; }					//Returns the number of layers in the z-direction.

	int Structure3DWithGrid::sizeX() const { return M; }						//Returns the number of Ez nodes on a line in x-direction across grid.
	int Structure3DWithGrid::sizeY() const { return N; }						//Returns the number of Ez nodes in y-direction across grid.
	int Structure3DWithGrid::sizeZ() const { return P; }						//Returns the number of Ex nodes in z-direction across grid.



	double Structure3DWithGrid::deltaxL(int i) const					//deltaxL(i) returns the x-direction node spacing in micrometers for x-layer with index i.
	{
		return theGrid.deltaxL(i);
	}

	double Structure3DWithGrid::deltayL(int j) const					//deltayL(j) returns the y-direction node spacing in micrometers for y-layer with index j.
	{
		return theGrid.deltayL(j);
	}

	double Structure3DWithGrid::deltazL(int k) const					//deltazL(i) returns the z-direction node spacing in micrometers for z-layer with index k.
	{
		return theGrid.deltazL(k);
	}


	int Structure3DWithGrid::nStepsXL(int i) const		//Returns the number of  Ez nodes across x-layer with index i on a line in the x-direction. 
	{
		return theGrid.nStepsXL(i);
	}
			
	int Structure3DWithGrid::nStepsYL(int j) const		//Returns the number of Ez nodes accross y-layer with index j on a line in the y-direction.
	{
		return theGrid.nStepsYL(j);
	}

	int Structure3DWithGrid::nStepsZL(int k) const		//Returns the number of Ex nodes accross z-layer with index k on a line in the z-direction exluding interface node at lowest z-value.
	{
		return theGrid.nStepsZL(k);
	}



	//*****************************************************************************************************************************************************
	//Function fund
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
			

	bool Structure3DWithGrid::find(int m, int n, int p, int& iInd, int& jInd, int& kInd, bool& interior, bool& xLower, bool& yLower, bool& zLower, bool& xBoundary, bool& yBoundary, bool& zBoundary) const
	{
		int i, j, k;
				
		int start, nstart, pStart;
		int numSteps;

		bool found, xFound, yFound, zFound;

		interior = true;
		xLower = false;
		yLower = false;
		zLower = false;

		xBoundary = false;
		yBoundary = false;
		zBoundary = false;

		xFound = false;
		yFound = false;
		zFound = false;
		found = false;

		start = 0;
		numSteps = 0;

		for (i = 0; i < theGrid.numLx() && !xFound; i++)
		{
			numSteps += theGrid.nStepsXL(i);

			if (m < numSteps)
			{
				iInd = i;
				xFound = true;

				if (m == start)
				{
					xLower = true;
					interior = false;
				}//End if

			}//End if

			start = numSteps;
		}//End for each x-layer index
						

		start = 0;
		numSteps = 0;

		for (j = 0; j < theGrid.numLy() && !yFound; j++)
		{
			numSteps += theGrid.nStepsYL(j);

			if (n < numSteps)
			{
				jInd = j;
				yFound = true;

				if (n == start)
				{
					yLower = true;
					interior = false;
				}
				//End if

			}//End if

			start = numSteps;
		}//End for each y-layer index


		start = 0;
		numSteps = 0;

		for (k = 0; k < theGrid.numLz() && !zFound; k++)
		{
			numSteps += theGrid.nStepsZL(k);

			if (p < numSteps)
			{
				kInd = k;
				zFound = true;

				if (p == start)
				{
					zLower = true;
						interior = false;
				}
				//End if

			}//End if
		}//End for each z-layer index


		if (xFound && yFound && zFound)
			found = true;
		else
		{
			found = false;
			//Determine if m,n,p is a boundary value

			if (!xFound)
				if (m == M - 1)
				{
					iInd = theGrid.numLx();
					xBoundary = true;
				}
				else
				{
					//error condition: m too high
				}//End if

			if (!yFound)
				if (n == N - 1)
				{
					jInd = theGrid.numLy();
					yBoundary = true;
				}
				else
				{
					//error condition: n too high
				}//End if

			if (!zFound)
				if (p == P - 1)
				{
					kInd = theGrid.numLz();
					zBoundary = true;
				}
				else
				{
					//error condition: p too high
				}//End if

			}//End if

			return found;

	}//End function find of Structure3DWithGrid


			 //******************************************************************************************************************************
			double Structure3DWithGrid::selectCEx(int i, int j, int k, bool found, bool xLower, bool yLower, bool zLower, bool xboundary, bool yBoundary, bool zBoundary,
													vector <Array3D>& CExF_for, vector <Array2D>& CExF_for_high_yGridBoundary, vector <Array2D>& CExF_for_high_zGridBoundary, Array1D& CExF_for_GridEdge)
			{

				/*
				enum CType_Ex{ Interior_Ex, yInterface_Ex, zInterface_Ex, yzEdge_Ex };
				
				enum CB_High_y_Type_Ex {Interior_of_yBoundary_element_Ex, zInterface_on_yBoundary_element_Ex};

				enum CB_High_z_Type_Ex {Interior_of_zBoundary_element_Ex, yInterface_on_zBoundary_element_Ex };

				*/


				if (found)
				{
					if (!yLower && !zLower)
						return CExF_for[Interior_Ex](i, j, k);
					else if (!yLower && zLower)
						return CExF_for[zInterface_Ex](i, j, k);
					else if (yLower && !zLower)
						return CExF_for[yInterface_Ex](i, j, k);
					else if (yLower && zLower)
						return CExF_for[yzEdge_Ex](i, j, k);
					//end if
				}
				else
						if (yBoundary && zBoundary)
							return CExF_for_GridEdge(i);

						else if (yBoundary && !zBoundary)
							if (!zLower)
								return CExF_for_high_yGridBoundary[Interior_of_yBoundary_element_Ex](i, k);
							else
								return CExF_for_high_yGridBoundary[zInterface_on_yBoundary_element_Ex](i, k);
							//End if

						else if (!yBoundary && zBoundary)
							if (!yLower)
								return CExF_for_high_zGridBoundary[Interior_of_zBoundary_element_Ex](i, j);
							else
								return CExF_for_high_zGridBoundary[yInterface_on_zBoundary_element_Ex](i, j);
							//End if
						else
						{
							//error CEx function should not have been called
							return 0;
						}//End if

						return 0;
			}//End function selectCEx


			 //*******************************************************************************************************************************
			double Structure3DWithGrid::selectCEy(int i, int j, int k, bool found, bool xLower, bool yLower, bool zLower, bool xBoundary, bool yBoundary, bool zBoundary,
				vector <Array3D>& CEyF_for, vector <Array2D>& CEyF_for_high_zGridBoundary)
			{
				

				/*	enum CType_Ey { Interior_Ey, xInterface_Ey, zInterface_Ey, xzEdge_Ey};

					enum CB_High_z_Type_Ey {Interior_of_zBoundary_element_Ey, xInterface_on_zBoundary_element_Ey};

					enum CB_High_x_Type_Ey { xBoundary_Ey, xBoundary_zInterface_Ey};
					
				*/

				if (found)
					if (!xLower && !zLower)
						return CEyF_for[Interior_Ey](i, j, k);
					else if(!xLower && zLower)
						return CEyF_for[zInterface_Ey](i, j, k);
					else if (xLower && !zLower)
						return CEyF_for[xInterface_Ey](i, j, k);
					else
						return  CEyF_for[xzEdge_Ey](i, j, k);
					//End if
				else
					if(zBoundary)
						if(!xLower)
							if(!yBoundary)
								return CEyF_for_high_zGridBoundary[Interior_of_zBoundary_element_Ey](i, j);
							else
							{//Error: CEy should not have been called or m,n,p out of range
								return 0;
							}//End if
						else
							if (!yBoundary)
								return CEyF_for_high_zGridBoundary[xInterface_on_zBoundary_element_Ey](i, j);
							else
							{//Error: CEy should not have been called or m,n,p out of range
								return 0;
							}//End if
						//End if
					else
					{
						//Error: CEy should not have been called or m,n,p out of range
						return 0;
					}//End if
				//End if


			}//End function selectCEy


			 //*********************************************************************************************************************************************************
			double Structure3DWithGrid::selectCEz(int i, int j, int k, bool found, bool xLower, bool yLower, bool yBoundary, bool zBoundary, vector <Array3D>& CEzF_for, vector <Array2D>& CEzF_for_high_yBoundary_for)
			{
				/*
				enum CType_Ez { Interior_Ez, xInterface_Ez, yInterface_Ez, xyEdge_Ez};

				enum CB_High_y_Type_Ez { Interior_of_yBoundary_element_Ez, xInterface_on_yBoundary_element_Ez}; */


				if (found)
					if (!xLower && !yLower)
						return CEzF_for[Interior_Ez](i, j, k);
					else if (!xLower && yLower)
						return CEzF_for[yInterface_Ez](i, j, k);
					else if (xLower && !yLower)
						return CEzF_for[xInterface_Ez](i, j, k);
					else
						return CEzF_for[xyEdge_Ez](i, j, k);
					//End if
				else
					if (yBoundary)
						if(!xLower)
							if(!zBoundary)
								return CEzF_for_high_yBoundary_for[Interior_of_yBoundary_element_Ez](i, k);
							else
							{//Error: CEz should not have been called or m,n,p out of range
								return 0;
							}//End if
						else
							if (!zBoundary)
								return CEzF_for_high_yBoundary_for[xInterface_on_yBoundary_element_Ez](i, k);
							else
							{//Error: CEz should not have been called or m,n,p out of range
								return 0;
							}//End if
						//End if
					else
					{
						//Error: CEz should not have been called or m,n,p out of range
						return 0;
					}//End if
				//End if


			}//End function selectCEz


			 //**********************************************************************************************************************************
			double Structure3DWithGrid::selectCHx(int i, int j, int k, bool found, bool xLower, bool yLower, bool zLower, bool xBoundary, bool yBoundary, bool zBoundary, vector <Array3D>& CHxF_for, Array2D& CHxF_for_high_xBoundary_for)
			{
				/* enum CType_Hx { Interior_Hx, xInterface_Hx}; */

				if (found)
					if (!xLower)
						return CHxF_for[Interior_Hx](i, j, k);
					else
						return CHxF_for[xInterface_Hx](i, j, k);
				//End if
				else
					if (xBoundary)
						if(!yBoundary && !zBoundary)
							return CHxF_for_high_xBoundary_for(j, k);
						else
						{	//Error: CHx should not have been called or m,n,p out of range
							return 0;
						}//End if
					else
					{
						//Error: CHx should not have been found or m,n,p out of range
						return 0;
					}//End if
				//End if

			}//End function selectCHx 


			 //**********************************************************************************************************************************
			double Structure3DWithGrid::selectCHy(int i, int j, int k, bool found, bool yLower, bool yBoundary, bool zBoundary, vector <Array3D>& CHyF_for, Array2D& CHyF_for_high_yBoundary_for)
			{

				/* enum CType_Hy { Interior_Hy, yInterface_Hy};*/

				if (found)
					if (!yLower)
						return CHyF_for[Interior_Hy](i, j, k);
					else
						return CHyF_for[yInterface_Hy](i, j, k);
				//End if
				else
					if (yBoundary)
						if(!zBoundary)
							return CHyF_for_high_yBoundary_for(i, k);
						else
						{	//Error: CHy should not have been called or m,n,p out of range
							return 0;
						}//End if
					else
					{
						//Error: CHy should not have been called or m,n,p out of range 
						return 0;
					}//End if
				//End if


			}//End function selectCHy



			 //**********************************************************************************************************************************
			double Structure3DWithGrid::selectCHz(int i, int j, int k, bool found, bool xLower, bool yLower, bool zLower, bool yBoundary, bool zBoundary,
													vector <Array3D>& CHzF_for, Array2D& CHzF_for_high_zGridBoundary)
			{
				/* enum CType_Hz{ Interior_Hz, zInterface_Hz};

				enum CB_High_z_Type_Hz { Interior_of_zBoundary_element_Hz}; */

				if (found)
					if (!zLower)
						return CHzF_for[Interior_Hz](i, j, k);
					else
						return CHzF_for[zInterface_Hz](i, j, k);
					//End if
				else
					if (zBoundary)
						if(!yBoundary)
							return CHzF_for_high_zGridBoundary(i, j);
						else
						{	//Error: CHz should not have been called or m,n,p out of range
							return 0;
						}//End if
					else
					{
						//Error: CHz should not have been called or m,n,p out of range
						return 0;
					}//End if
				//End if



			}//End function selectCHz 


			 //**********************************************************************************************************************************
			double Structure3DWithGrid::CExE(int m, int n, int p)
			{
				bool interior = true;
				bool xLower = false;
				bool yLower = false;
				bool zLower = false;

				bool xBoundary = false;
				bool yBoundary = false;
				bool zBoundary = false;

				bool found;

				int i, j, k;

				double val;


				found = find(m, n, p, i, j, k, interior, xLower, yLower, zLower, xBoundary, yBoundary, zBoundary);

				if (!found)
				{/*Set error message if m,n,p  out of range*/
				}

				val = selectCEx(i, j, k, found, xLower, yLower, zLower, xBoundary, yBoundary, zBoundary, CExE_for, CExE_on_high_yBoundary_for, CExE_on_high_zBoundary_for, CExE_for_GridEdge);

				return val;

			}//End function CExE


			 //******************************************************************************************************************************
			double Structure3DWithGrid::CExHz(int m, int n, int p)
			{
				bool interior = true;
				bool xLower = false;
				bool yLower = false;
				bool zLower = false;

				bool xBoundary = false;
				bool yBoundary = false;
				bool zBoundary = false;

				bool found;

				int i, j, k;

				double val;


				found = find(m, n, p, i, j, k, interior, xLower, yLower, zLower, xBoundary, yBoundary, zBoundary);

				if (!found)
				{/*Set error message if m,n,p  out of range*/
				}

				val = selectCEx(i, j, k, found, xLower, yLower, zLower, xBoundary, yBoundary, zBoundary, CExHz_for, CExHz_on_high_yBoundary_for, CExHz_on_high_zBoundary_for, CExHz_for_GridEdge);

				return val;


			}//End function CExHy


			 //******************************************************************************************************************************
			double Structure3DWithGrid::CExHy(int m, int n, int p)
			{
				bool interior = true;
				bool xLower = false;
				bool yLower = false;
				bool zLower = false;

				bool xBoundary = false;
				bool yBoundary = false;
				bool zBoundary = false;

				bool found;

				int i, j, k;

				double val;


				found = find(m, n, p, i, j, k, interior, xLower, yLower, zLower, xBoundary, yBoundary, zBoundary);

				if (!found)
				{/*Set error message if m,n,p  out of range*/
				}

				val = selectCEx(i, j, k, found, xLower, yLower, zLower, xBoundary, yBoundary, zBoundary, CExHy_for, CExHy_on_high_yBoundary_for, CExHy_on_high_zBoundary_for, CExHy_for_GridEdge);

				return val;


			}//End function CExHy


			 //******************************************************************************************************************************
			double Structure3DWithGrid::CEyE(int m, int n, int p)
			{
				bool interior = true;
				bool xLower = false;
				bool yLower = false;
				bool zLower = false;

				bool xBoundary = false;
				bool yBoundary = false;
				bool zBoundary = false;

				bool found;

				int i, j, k;

				double val;


				found = find(m, n, p, i, j, k, interior, xLower, yLower, zLower, xBoundary, yBoundary, zBoundary);

				if (!found)
				{/*Set error message if m,n,p  out of range*/
				}

				val = selectCEy(i, j, k, found, xLower, yLower, zLower, xBoundary, yBoundary, zBoundary, CEyE_for, CEyE_on_high_zBoundary_for);

				return val;

			}//End function CEyE


			 //******************************************************************************************************************************
			double Structure3DWithGrid::CEyHx(int m, int n, int p)
			{
				bool interior = true;
				bool xLower = false;
				bool yLower = false;
				bool zLower = false;

				bool xBoundary = false;
				bool yBoundary = false;
				bool zBoundary = false;

				bool found;

				int i, j, k;

				double val;


				found = find(m, n, p, i, j, k, interior, xLower, yLower, zLower, xBoundary, yBoundary, zBoundary);

				if (!found)
				{/*Set error message if m,n,p  out of range*/
				}

				val = selectCEy(i, j, k, found, xLower, yLower, zLower, xBoundary, yBoundary, zBoundary, CEyHx_for, CEyHx_on_high_zBoundary_for);

				return val;

			}//End function CEyHx


			 //******************************************************************************************************************************
			double Structure3DWithGrid::CEyHz(int m, int n, int p)
			{
				bool interior = true;
				bool xLower = false;
				bool yLower = false;
				bool zLower = false;

				bool xBoundary = false;
				bool yBoundary = false;
				bool zBoundary = false;

				bool found;

				int i, j, k;

				double val;


				found = find(m, n, p, i, j, k, interior, xLower, yLower, zLower, xBoundary, yBoundary, zBoundary);

				if (!found)
				{/*Set error message*/
				}

				val = selectCEy(i, j, k, found, xLower, yLower, zLower, xBoundary, yBoundary, zBoundary, CEyHz_for, CEyHz_on_high_zBoundary_for);

				return val;

			}//End function CEyHz


			 //******************************************************************************************************************************
			double Structure3DWithGrid::CHzH(int m, int n, int p)
			{
				bool interior = true;
				bool xLower = false;
				bool yLower = false;
				bool zLower = false;

				bool xBoundary = false;
				bool yBoundary = false;
				bool zBoundary = false;

				bool found;

				int i, j, k;

				double val;


				found = find(m, n, p, i, j, k, interior, xLower, yLower, zLower, xBoundary, yBoundary, zBoundary);

				if (!found)
				{/*Set error message if m,np out of range*/
				}

				val = selectCHz(i, j, k, found, xLower, yLower, zLower, yBoundary, zBoundary, CHzH_for, CHzH_on_high_zBoundary_for);

				return val;

			}//End function CHzH


			 //******************************************************************************************************************************
			double Structure3DWithGrid::CHzEy(int m, int n, int p)
			{
				bool interior = true;
				bool xLower = false;
				bool yLower = false;
				bool zLower = false;

				bool xBoundary = false;
				bool yBoundary = false;
				bool zBoundary = false;

				bool found;

				int i, j, k;

				double val;


				found = find(m, n, p, i, j, k, interior, xLower, yLower, zLower, xBoundary, yBoundary, zBoundary);

				if (!found)
				{/*Set error message if m,np out of range*/
				}

				val = selectCHz(i, j, k, found, xLower, yLower, zLower, yBoundary, zBoundary, CHzEy_for, CHzEy_on_high_zBoundary_for);

				return val;

			}//End function CHzEy


			 //**********************************************************************************************************************************
			double Structure3DWithGrid::CHzEx(int m, int n, int p)
			{
				bool interior = true;
				bool xLower = false;
				bool yLower = false;
				bool zLower = false;

				bool xBoundary = false;
				bool yBoundary = false;
				bool zBoundary = false;

				bool found;

				int i, j, k;

				double val;


				found = find(m, n, p, i, j, k, interior, xLower, yLower, zLower, xBoundary, yBoundary, zBoundary);

				if (!found)
				{/*Set error message if m,np out of range*/
				}

				val = selectCHz(i, j, k, found, xLower, yLower, zLower, yBoundary, zBoundary, CHzEx_for, CHzEx_on_high_zBoundary_for);

				return val;

			}//End function CHzEx


			 //**********************************************************************************************************************************
			double Structure3DWithGrid::CHxH(int m, int n, int p)
			{
				bool interior = true;
				bool xLower = false;
				bool yLower = false;
				bool zLower = false;

				bool xBoundary = false;
				bool yBoundary = false;
				bool zBoundary = false;

				bool found;

				int i, j, k;

				double val;


				found = find(m, n, p, i, j, k, interior, xLower, yLower, zLower, xBoundary, yBoundary, zBoundary);

				if (!found)
				{/*Set error message if m,n,p out of range*/
				}

				val = selectCHx(i, j, k, found, xLower, yLower, zLower, xBoundary, yBoundary, zBoundary, CHxH_for, CHxH_on_high_xBoundary_for);

				return val;

			}//End function CHxH


			 //**********************************************************************************************************************************
			double Structure3DWithGrid::CHxEy(int m, int n, int p)
			{
				bool interior = true;
				bool xLower = false;
				bool yLower = false;
				bool zLower = false;

				bool xBoundary = false;
				bool yBoundary = false;
				bool zBoundary = false;

				bool found;

				int i, j, k;

				double val;


				found = find(m, n, p, i, j, k, interior, xLower, yLower, zLower, xBoundary, yBoundary, zBoundary);

				if (!found)
				{/*Set error message if m,n,p out of range*/
				}

				val = selectCHx(i, j, k, found, xLower, yLower, zLower, xBoundary, yBoundary, zBoundary, CHxEy_for, CHxEy_on_high_xBoundary_for);

				return val;

			}//End function CHxEy


			 //**********************************************************************************************************************************
			double Structure3DWithGrid::CHxEz(int m, int n, int p)
			{
				bool interior = true;
				bool xLower = false;
				bool yLower = false;
				bool zLower = false;

				bool xBoundary = false;
				bool yBoundary = false;
				bool zBoundary = false;

				bool found;

				int i, j, k;

				double val;


				found = find(m, n, p, i, j, k, interior, xLower, yLower, zLower, xBoundary, yBoundary, zBoundary);

				if (!found)
				{/*Set error message if m,n,p out of range*/
				}

				val = selectCHx(i, j, k, found, xLower, yLower, zLower, xBoundary, yBoundary, zBoundary, CHxEz_for, CHxEz_on_high_xBoundary_for);

				return val;

			}//End function CHxEz


			 //**********************************************************************************************************************************
			double Structure3DWithGrid::CHyH(int m, int n, int p)
			{
				bool interior = true;
				bool xLower = false;
				bool yLower = false;
				bool zLower = false;

				bool xBoundary = false;
				bool yBoundary = false;
				bool zBoundary = false;

				bool found;

				int i, j, k;

				double val;


				found = find(m, n, p, i, j, k, interior, xLower, yLower, zLower, xBoundary, yBoundary, zBoundary);

				if (!found)
				{/*Set error message if m,np out of range*/
				}

				val = selectCHy(i, j, k, found, yLower, yBoundary, zBoundary, CHyH_for, CHyH_on_high_yBoundary_for);

				return val;

			}//End function CHyH


			 //**********************************************************************************************************************************
			double Structure3DWithGrid::CHyEx(int m, int n, int p)
			{
				bool interior = true;
				bool xLower = false;
				bool yLower = false;
				bool zLower = false;

				bool xBoundary = false;
				bool yBoundary = false;
				bool zBoundary = false;

				bool found;

				int i, j, k;

				double val;


				found = find(m, n, p, i, j, k, interior, xLower, yLower, zLower, xBoundary, yBoundary, zBoundary);

				if (!found)
				{/*Set error message if m,np out of range*/
				}

				val = selectCHy(i, j, k, found, yLower, yBoundary, zBoundary, CHyEx_for, CHyEx_on_high_yBoundary_for);

				return val;

			}//End function CHyEx


			 //**********************************************************************************************************************************
			double Structure3DWithGrid::CHyEz(int m, int n, int p)
			{
				bool interior = true;
				bool xLower = false;
				bool yLower = false;
				bool zLower = false;

				bool xBoundary = false;
				bool yBoundary = false;
				bool zBoundary = false;

				bool found;

				int i, j, k;

				double val;


				found = find(m, n, p, i, j, k, interior, xLower, yLower, zLower, xBoundary, yBoundary, zBoundary);

				if (!found)
				{/*Set error message if m,np out of range*/
				}

				val = selectCHy(i, j, k, found, yLower, yBoundary, zBoundary, CHyEz_for, CHyEz_on_high_yBoundary_for);

				return val;


			}//End function CHyEz


			 //**********************************************************************************************************************************
			double Structure3DWithGrid::CEzE(int m, int n, int p)
			{
				bool interior = true;
				bool xLower = false;
				bool yLower = false;
				bool zLower = false;

				bool xBoundary = false;
				bool yBoundary = false;
				bool zBoundary = false;

				bool found;

				int i, j, k;

				double val;


				found = find(m, n, p, i, j, k, interior, xLower, yLower, zLower, xBoundary, yBoundary, zBoundary);

				if (!found)
				{/*Set error message if m,n,pout of range*/
				}

				val = selectCEz(i, j, k, found, xLower, yLower, yBoundary, zBoundary, CEzE_for, CEzE_on_high_yBoundary_for);

				return val;

			}//End function CEzE


			 //**********************************************************************************************************************************
			double Structure3DWithGrid::CEzHy(int m, int n, int p)
			{
				bool interior = true;
				bool xLower = false;
				bool yLower = false;
				bool zLower = false;

				bool xBoundary = false;
				bool yBoundary = false;
				bool zBoundary = false;

				bool found;

				int i, j, k;

				double val;


				found = find(m, n, p, i, j, k, interior, xLower, yLower, zLower, xBoundary, yBoundary, zBoundary);

				if (!found)
				{/*Set error message if m,n,p out of range*/
				}

				val = selectCEz(i, j, k, found, xLower, yLower, yBoundary, zBoundary, CEzHy_for, CEzHy_on_high_yBoundary_for);

				return val;

			}//End function CEzHy


			 //**********************************************************************************************************************************
			double Structure3DWithGrid::CEzHx(int m, int n, int p)
			{

				bool interior = true;
				bool xLower = false;
				bool yLower = false;
				bool zLower = false;

				bool xBoundary = false;
				bool yBoundary = false;
				bool zBoundary = false;

				bool found;

				int i, j, k;

				double val;


				found = find(m, n, p, i, j, k, interior, xLower, yLower, zLower, xBoundary, yBoundary, zBoundary);

				if (!found)
				{/*Set error message if m,n,p out of range*/
				}

				val = selectCEz(i, j, k, found, xLower, yLower, yBoundary, zBoundary, CEzHx_for, CEzHx_on_high_yBoundary_for);

				return val;

			}//End function CEzHx


			bool Structure3DWithGrid::Is_a_unitCell() const { return theStructure.Is_a_unitCell(); }			//Returns true if the structure is the unit cell of a larger periodic structure
			bool Structure3DWithGrid::Is_periodic_in_yDirection() const { return theStructure.Is_periodic_in_yDirection(); }			//Returns true if the larger structure is periodic in the y-direction
			bool Structure3DWithGrid::Is_periodic_in_zDirection() const { return theStructure.Is_periodic_in_zDirection(); }						//Returns true if the larger structure is periodic in the z-direction


			//****************************************************************************************************************
			double Structure3DWithGrid::delAx(int mInd, int nInd, int pInd) const		//Returns the x-face area of the Yee cell with indexes mInd, nInd, pInd
			{
				return theGrid.delAx(mInd, nInd, pInd);
			}


			int Structure3DWithGrid::numSteps_xDir_throuh(int iInd) const
			{
				return theGrid.numSteps_xDir_throuh(iInd);
			}

			int Structure3DWithGrid::numSteps_yDir_throuh(int iInd) const
			{
				return theGrid.numSteps_yDir_throuh(iInd);
			}

			int Structure3DWithGrid::numSteps_zDir_throuh(int iInd) const
			{
				return theGrid.numSteps_zDir_throuh(iInd);
			}

			double Structure3DWithGrid::Sc3D()	const
			{
				return Sc;
			}


//****************************************************************************************************************************************
	//functions for Grid3D

			Grid3D::Grid3D(Structure3D& str, double minLambda, Report& myReport) : xLThick(str.numLx() + 2), yLThick(str.numLy()), zLThick(str.numLz()), numStepsXL(str.numLx() + 2),
				numStepsYL(str.numLy()), numStepsZL(str.numLz()), deltxL(str.numLx() + 2), deltyL(str.numLy()), deltzL(str.numLz())
			{
				int i;
				int j;
				int k;
				int numYinc;

				double dy;

				if (minLambda <= 0)
				{
					myReport.setErrorFlag(true);	
					myReport.setMessage("Error: Attempt to set negative pulse width / wavelength. Using default.");
					myReport.setVal(2);

				}
				else
				{
					//Error flag should have been initialized to false. So no need to set to false
					myReport.setErrorFlag(true);
					myReport.setMessage("No error: Acceptable pulse width");
					myReport.setVal(minLambda);


				}//End if

				//Set the numbers of grid layers in each of the coordinate directions
				nLx = str.numLx() + 2;
				nLy = str.numLy();
				nLz = str.numLz();


				//Calculate x-direction grid properties

				//Set number of spacial steps across the grid layer above the structure 
				numStepsXL[0] = NLambda;

				//x-index of the Ez field nodes at the topsurface of the sructure.
				mD = numStepsXL[0];

				M = numStepsXL[0] + 1;	//Starting value of M before generating grid in structure


				//Calculate grid in the structure x-layers

				//x-direction properties
			

				for (i = 1; i < nLx - 1; i++)
				{
					xLThick[i] = str.xLth(i - 1);
					numStepsXL[i] = max((int)((xLThick[i] / minLambda)*NLambda + 0.5), minCpt);	
					deltxL[i] = xLThick[i] / numStepsXL[i];

					M = M + numStepsXL[i];

					//To find minimum deltx
					if (i == 1)
						minDx = deltxL[i];
					else
						if (deltxL[i] < minDx)
							minDx = deltxL[i];
					//End if

				}//End for each x-layer in structure

				 //Finish setting grid values for layer above structure.
				deltxL[0] = deltxL[1];
				xLThick[0] = numStepsXL[0] * deltxL[0];


				//Calculate values for grid layer beneath structure and add number of nodes to M.
				numStepsXL[nLx - 1] = NLambda;
				M = M + numStepsXL[nLx - 1];					//M becomes the number of Ez node planes
				deltxL[nLx - 1] = deltxL[nLx - 2];				//OK not to be deltx[0] ?
				xLThick[nLx - 1] = numStepsXL[nLx - 1] * deltxL[nLx - 1];

				//No need to compare minDx with deltxL[nLx - 1] since deltxL[nLx - 1] = deltxL[nLx - 2] 


				//Calculate y - direction grid properties

				minDy = deltyL[0];

				N = 1;

				for (j = 0; j < nLy; j++)
				{
					yLThick[j] = str.yLth(j);
					numStepsYL[j] = max((int)((yLThick[j] / minLambda)*NLambda + 0.5), minCpt);	
					deltyL[j] = yLThick[j] / numStepsYL[j];

					N = N + numStepsYL[j];

					//To find minimum delty
					if (j == 0)
						minDy = deltyL[j];
					else
						if (deltyL[j] < minDy)
							minDy = deltyL[j];
					//End if


				}//End for each y-layer in grid 


				//Calculate z-direction grid properties

				P = 1;

				for (k = 0; k < nLz; k++)
				{
					zLThick[k] = str.zLth(k);
					numStepsZL[k] = max((int)((zLThick[k] / minLambda)*NLambda + 0.5), minCpt);	
					deltzL[k] = zLThick[k] / numStepsZL[k];

					P = P + numStepsZL[k];

					//To find minimum deltz
					if (k == 0)
						minDz = deltzL[k];
					else
						if (deltzL[k] < minDz)
							minDz = deltzL[k];
					//End if

				}//End for each z-layer in grid 

				//Determine the minimum spacial increment
				minDel = minDx;
				if (minDy < minDel)
					minDel = minDy;
				//End if
				if (minDz < minDel)
					minDel = minDz;
				//End if

			}//End Constructor Grid3D



			Grid3D::Grid3D(Grid3D& aGrid)
			{
				M = aGrid.sizeX();							//Number of Ez nodes across grid in x-direction
				N = aGrid.sizeY();							//Number of Ez nodes across grid in y-direction
				P = aGrid.sizeZ();							//Number of Ex nodes across grid in z-direction

				int mD;							//Device/structure starts at x = mD*deltx[0].

				nLx = aGrid.nLx;						//Number of layers in the x-direction
				nLy = aGrid.nLy;						//Number of layers in the y-direction
				nLz = aGrid.nLz;						//Number of layers in the z-direction

				xLThick = aGrid.xLThick;		//Vector of x-layer thicknesses
				yLThick = aGrid.yLThick;		//Vector of y-layer thicknesses
				zLThick = aGrid.zLThick;		//Vector of z-layer thicknesses

				numStepsXL = aGrid.numStepsXL;			//numStepsXL[i] is the number spacial steps across layer with index i on a line in the x-direction. 
				numStepsYL = aGrid.numStepsYL;			//numStepsYL[j] is the number of spacial steps across layer with index j in the y-direction.
				numStepsZL = aGrid.numStepsZL;			//numStepsZL[k] is the number of spacial steps across layer with index k in the z-direction.


				deltxL = aGrid.deltxL;			//deltxL[i] is the x-direction node spacing in micrometers for x-layer with index i.
				deltyL = aGrid.deltyL;			//deltyL[j] is the y-direction node spacing in micrometers for y-layer with index j.
				deltzL = aGrid.deltzL;			//deltzL[k] is the z-direction node spacing in micrometers for z-layer with index k.

				minDx = aGrid.minDx;					//Minimum deltxL over all x-layers
				minDy = aGrid.minDy;					//Minimum deltyL over all y-layers		
				minDz = aGrid.minDz;					//Minimum deltzL over all z-layers

				minDel = aGrid.minDel;					//Minimum of all spacial increments

			}//End copy constructor for Grid3D



			 //Accessors
			int Grid3D::sizeX() const { return M; }						//Returns the number of Ez nodes on a line in x-direction across grid.
			int Grid3D::sizeY() const { return N; }						//Returns the number of Ez nodes in y-direction across grid.
			int Grid3D::sizeZ() const { return P; }						//Returns the number of Ex nodes in z-direction across grid.


			int Grid3D::numLx() const { return nLx; }						//Returns the number of layers normal to the x-direction.
			int Grid3D::numLy() const { return nLy; }					//Returns the number of layers in the y-direction.		
			int Grid3D::numLz() const { return nLz; }					//Returns the number of layers in the z-direction.

			int Grid3D::nStepsXL(int i) const { return numStepsXL[i]; }		//Returns the number of  Ez nodes across x-layer with index i on a line in the x-direction. 
			int Grid3D::nStepsYL(int j) const { return numStepsYL[j]; }		//Returns the number of Ez nodes accross y-layer with index j on a line in the y-direction exluding interface node at lowest y-value.
			int Grid3D::nStepsZL(int k) const { return numStepsZL[k]; }		//Returns the number of Ex nodes accross z-layer with index k on a line in the z-direction exluding interface node at lowest z-value.

			double Grid3D::deltaxL(int i) const { return deltxL[i]; }	//deltaxL(i) returns the x-direction node spacing in micrometers for x-layer with index i.
			double Grid3D::deltayL(int j) const { return deltyL[j]; }	//deltayL(j) returns the y-direction node spacing in micrometers for y-layer with index j.
			double Grid3D::deltazL(int k) const { return deltzL[k]; }	//deltazL(i) returns the z-direction node spacing in micrometers for z-layer with index k.

			double Grid3D::minDelta() const { return minDel; }			//Returns the minimum of spacial steps (minimum Yee cell edge length)


			int Grid3D::numSteps_xDir_throuh(int indLast) const			//Returns the number of steps across the grid in the x-direction beginning at x = 0 up through end of layer with index indLast.
			{
				int total = 0;
				int i;

				for (i = 0; i < indLast + 1; i++)
					total += nStepsXL(i);
				//End for

				return total;
			}//End function numSteps_xDir_throuh


			int Grid3D::numSteps_yDir_throuh(int indLast) const			//Returns the number of steps across the grid in the y-direction beginning at y = 0 up through end of layer with index indLast.
			{
				int total = 0;
				int j;

				for (j = 0; j < indLast + 1; j++)
					total += nStepsYL(j);
				//End for

				return total;
			}//End function numSteps_yDir_throuh


			int Grid3D::numSteps_zDir_throuh(int indLast) const			//Returns the number of steps across the grid in the z-direction beginning at z = 0 up through end of layer with index indLast.
			{
				int total = 0;
				int k;

				for (k = 0; k < indLast + 1; k++)
					total += nStepsYL(k);
				//End for

				return total;
			}//End function numSteps_zDir_throuh


			double Grid3D::delAx(int mInd, int nInd, int pInd) const		//Returns the x-face area of the Yee cell with indexes mInd, nInd, pInd
			{
				bool interior = true;
				bool xLower = false;
				bool yLower = false;
				bool zLower = false;

				bool xBoundary = false;
				bool yBoundary = false;
				bool zBoundary = false;

				bool found;

				int i, j, k;

				found = find(mInd, nInd, pInd, i, j, k, interior, xLower, yLower, zLower, xBoundary, yBoundary, zBoundary);

				if (found)
					return  deltayL(j)*deltazL(k);
				else
				{
					//Error condition: delAxshould not have been called or mInd, nInd, pInd out of range
					return 0;
				}//End if

			}//End function delAx


			bool Grid3D::find(int m, int n, int p, int& iInd, int& jInd, int& kInd, bool& interior, bool& xLower, bool& yLower, bool& zLower, bool& xBoundary, bool& yBoundary, bool& zBoundary) const
			{
				int i, j, k;

				int start, nstart, pStart;
				int numSteps;

				bool found, xFound, yFound, zFound;

				interior = true;
				xLower = false;
				yLower = false;
				zLower = false;

				xBoundary = false;
				yBoundary = false;
				zBoundary = false;

				xFound = false;
				yFound = false;
				zFound = false;
				found = false;


				start = 0;
				numSteps = 0;

				for (i = 0; i < numLx() && !xFound; i++)
				{

					numSteps += nStepsXL(i);

					if (m < numSteps)
					{
						iInd = i;
						xFound = true;

						if (m == start)
						{
							xLower = true;
							interior = false;
						}//End if

					}//End if

					start = numSteps;
				}//End for each x-layer index


				start = 0;
				numSteps = 0;

				for (j = 0; j < numLy() && !yFound; j++)
				{
					numSteps += nStepsYL(j);

					if (n < numSteps)
					{
						jInd = j;
						yFound = true;

						if (n == start)
						{
							yLower = true;
							interior = false;
						}
						//End if

					}//End if

					start = numSteps;
				}//End for each y-layer index


				start = 0;
				numSteps = 0;

				for (k = 0; k < numLz() && !zFound; k++)
				{
					numSteps += nStepsZL(k);

					if (p < numSteps)
					{
						kInd = k;
						zFound = true;

						if (p == start)
						{
							zLower = true;
							interior = false;
						}
						//End if

					}//End if
				}//End for each z-layer index


				if (xFound && yFound && zFound)
					found = true;
				else
				{
					found = false;
					//Determine if m,n,p is a boundary value

					if (!xFound)
						if (m == M - 1)
						{
							iInd = numLx();
							xBoundary = true;
						}
						else
						{
							//error condition: m too high
						}//End if

					if (!yFound)
						if (n == N - 1)
						{
							jInd = numLy();
							yBoundary = true;
						}
						else
						{
							//error condition: n too high
						}//End if

					if (!zFound)
						if (n == N - 1)
						{
							kInd = numLz();
							zBoundary = true;
						}
						else
						{
							//error condition: p too high
						}//End if

				}//End if

				return found;

			}//End function find of Grid3D



//****************************************************************************************************************************************
//Functions for Grid1D

			Grid1D::Grid1D(Structure1D& dev, double minLambda)
			{
				int i;

				numLayers = dev.numL() + 2; //Grid has a layer above the device and a layer below the device.

				nNodes = new int[numLayers];
				thickness = new double[numLayers];
				deltx = new double[numLayers];

				nNodes[0] = minCpt;
				//nNodes[numLayers - 1] = minCpt;
				mD = nNodes[0] - 1;

				N = nNodes[0];

				for (i = 1; i < numLayers - 1; i++)
				{
					thickness[i] = dev.th(i - 1);
					nNodes[i] = max((int)((thickness[i] / minLambda)*NLambda + 0.5), minCpt);
					deltx[i] = thickness[i] / nNodes[i];
					N = N + nNodes[i];
				}//End for each device layer

				deltx[0] = deltx[1];
				thickness[0] = (nNodes[0] - 1) * deltx[0] + deltx[0] / 2;

				nNodes[numLayers - 1] = minCpt;
				N = N + nNodes[numLayers - 1];

				//Find minimum delx for grid above bottom layer
				minDeltx = deltx[0];
				for (i = 0; i < numLayers - 1; i++)
					if (deltx[i] < minDeltx)
						minDeltx = deltx[i];

				deltx[numLayers - 1] = minDeltx;
				thickness[numLayers - 1] = nNodes[numLayers - 1] * deltx[numLayers - 1];

			}//End Constructor



			Grid1D::Grid1D(int anNodes, double adelx)
			{
				int i;

				numLayers = 2;

				nNodes = new int[numLayers];
				thickness = new double[numLayers];
				deltx = new double[numLayers];

				nNodes[0] = anNodes;
				deltx[0] = adelx;
				thickness[0] = (nNodes[0] - 1) * deltx[0] + deltx[0] / 2;

				deltx[numLayers - 1] = deltx[numLayers - 2];
				nNodes[numLayers - 1] = minCpt;
				thickness[numLayers - 1] = nNodes[numLayers - 1] * deltx[numLayers - 1];

				N = nNodes[0] + nNodes[numLayers - 1];

				mD = nNodes[0] - 1;

				//Minimum delx for this grid is the uniform delx
				minDeltx = deltx[0];

			}//End Constructor Grid1D


			 //Accessors

			int Grid1D::size()const { return N; }
			int Grid1D::numL()const { return numLayers; }
			int Grid1D::nNds(int layer)const { return nNodes[layer]; }
			double Grid1D::delx(int layer)const { return deltx[layer]; }	//Returns delx for layer sub i of grid 
			double Grid1D::th(int layer)const { return thickness[layer]; }	//Returns thicknes of layer sub i of grid
			double Grid1D::minDelx() { return minDeltx; }



	//**********************************************************************************************************



			 //*********************************************************************************************************************************

	//Cut outs

			/*-------------------

			void Structure3DWithGrid::computeBoundaryUpdateCoefficientsForBoundaryGridElements()
			{
				int i;
				int j;
				int k;

				double loss;
				double permit;
				double cond;
				double permea;

				double deltayYboundary;
				double deltazZboundary;


				//Compute coefficients for nodes on boundaries of region above structure if the structure is the unit cell of a larger periodic structure
				i = 0;

				if (theStructure.Is_periodic_in_yDirection())
				{
					//Compute coefficients for y-boundary nodes
					for (k = 0; k < nLz; k++)
					{
						deltayYboundary = theGrid.deltayL(nLy - 1) / 2 + theGrid.deltayL(0) / 2;		//Periodicity

						CEzE_for_Grid[yBoundary_Ez](i, k) = 1;
						CEzHy_for_Grid[yBoundary_Ez](i, k) = delt / (eps0*theGrid.deltaxL(i));
						CEzHx_for_Grid[yBoundary_Ez](i, k) = delt / (eps0*deltayYboundary);

						CEzE_for_Grid[yBoundary_zInterface_Ez](i, k) = 1;						//not used if k = nLz - 1
						CEzHy_for_Grid[yBoundary_zInterface_Ez](i, k) = CEzHy_for_Grid[yBoundary_Ez](i, k);	//not used if k = nLz - 1
						CEzHx_for_Grid[yBoundary_zInterface_Ez](i, k) = CEzHx_for_Grid[yBoundary_Ez](i, k);	//not used if k = nLz - 1

						CExE_for_Grid[yBoundary_Ex](i, k) = 1;
						CExHy_for_Grid[yBoundary_Ex](i, k) = delt / (eps0*theGrid.deltazL(k));
						CExHz_for_Grid[yBoundary_Ex](i, k) = delt / (eps0*deltayYboundary);

						cond = theStructure.sigma(i, nLy - 1, k) / 4 + theStructure.sigma(i, 0, k) / 4;						//zero conductivity above structure
						permit = eps0 / 4 + eps0 / 4 + theStructure.eps(i, nLy - 1, k) / 4 + theStructure.eps(i, 0, k) / 4;
						loss = cond*delt / (2 * permit);

						CExE_for_Grid[yBoundary_xInterface_Ex](i, k) = (1 - loss) / (1 + loss);
						CExHy_for_Grid[yBoundary_xInterface_Ex](i, k) = (delt / (permit*theGrid.deltazL(k))) / (1 + loss);
						CExHz_for_Grid[yBoundary_xInterface_Ex](i, k) = (delt / (permit*deltayYboundary)) / (1 + loss);

					}//End for

					deltazZboundary = theGrid.deltazL(nLz - 1) / 2 + theGrid.deltazL(0) / 2;		//Periodicity

					CExE_for_GridEdge[yBoundary_zBoundary_Ex](i) = 1;
					CExHy_for_GridEdge[yBoundary_zBoundary_Ex](i) = delt / (eps0*deltazZboundary);
					CExHz_for_GridEdge[yBoundary_zBoundary_Ex](i) = delt / (eps0*deltayYboundary);

					if (theStructure.Is_periodic_in_zDirection())
					{
						cond = theStructure.sigma(i, nLy - 1, nLz - 1) / 8 + theStructure.sigma(i, 0, nLz - 1) / 8 + theStructure.sigma(i, nLy - 1, 0) / 8 + theStructure.sigma(i, 0, 0) / 8;
						//Periodicity. Zero conductivity above structure
						permit = eps0 / 8 + eps0 / 8 + eps0 / 8 + eps0 / 8 + theStructure.eps(i, nLy - 1, nLz - 1) / 8 + theStructure.eps(i, 0, nLz - 1) / 8 + theStructure.eps(i, nLy - 1, 0) / 8 + theStructure.eps(i, 0, 0) / 8;
						loss = cond*delt / (2 * permit);

						CExE_for_GridEdge[yBoundary_zBoundary_xInterface_Ex](i) = (1 - loss) / (1 + loss);
						CExHy_for_GridEdge[yBoundary_zBoundary_xInterface_Ex](i) = (delt / (permit*deltazZboundary)) / (1 + loss);
						CExHz_for_GridEdge[yBoundary_zBoundary_xInterface_Ex](i) = (delt / (permit*deltayYboundary)) / (1 + loss);
					}
					else
					{
						cond = theStructure.sigma(i, nLy - 1, nLz - 1) / 8 + theStructure.sigma(i, 0, nLz - 1) / 8;						//zero conductivity above structure and outside z-boundary
						permit = eps0 / 8 + eps0 / 8 + eps0 / 8 + eps0 / 8 + theStructure.eps(i, nLy - 1, nLz - 1) / 8 + theStructure.eps(i, 0, nLz - 1) / 8;	//Permitivity outside z-boundary assumed eps0
						loss = cond*delt / (2 * permit);

						CExE_for_GridEdge[yBoundary_zBoundary_xInterface_Ex](i) = (1 - loss) / (1 + loss);
						CExHy_for_GridEdge[yBoundary_zBoundary_xInterface_Ex](i) = (delt / (permit*theGrid.deltazL(k))) / (1 + loss);
						CExHz_for_GridEdge[yBoundary_zBoundary_xInterface_Ex](i) = (delt / (permit*deltayYboundary)) / (1 + loss);
					}//End if

				}//End if

				if (theStructure.Is_periodic_in_zDirection())
				{
					//Compute coefficients for z-boundary nodes
					for (j = 0; j < nLy; j++)
					{
						deltazZboundary = theGrid.deltazL(nLz - 1) / 2 + theGrid.deltazL(0) / 2;		//Periodicity

						CExE_for_Grid[zBoundary_Ex](i, j) = 1;
						CExHz_for_Grid[zBoundary_Ex](i, j) = delt / (eps0*theGrid.deltayL(j));
						CExHy_for_Grid[zBoundary_Ex](i, j) = delt / (eps0*deltazZboundary);

						cond = theStructure.sigma(i, j, nLz - 1) / 4 + theStructure.sigma(i, j, 0) / 4;						//zero conductivity above structure
						permit = eps0 / 4 + eps0 / 4 + theStructure.eps(i, j, nLz - 1) / 4 + theStructure.eps(i, j, 0) / 4;
						loss = cond*delt / (2 * permit);

						CExE_for_Grid[zBoundary_xInterface_Ex](i, j) = (1 - loss) / (1 + loss);
						CExHz_for_Grid[zBoundary_xInterface_Ex](i, j) = delt / (permit*theGrid.deltayL(j));
						CExHy_for_Grid[zBoundary_xInterface_Ex](i, j) = delt / (permit*deltazZboundary);

						CEyE_for_Grid[zBoundary_Ey](i, j) = 1;
						CEyHx_for_Grid[zBoundary_Ey](i, j) = delt / (eps0 *deltazZboundary);
						CEyHz_for_Grid[zBoundary_Ey](i, j) = delt / (eps0 *theGrid.deltaxL(i));

						if (j < nLy - 1)
						{
							cond = theStructure.sigma(i - 1, j, nLz - 1) / 4 + theStructure.sigma(i - 1, j, 0) / 4 + theStructure.sigma(i - 1, j + 1, nLz - 1) / 4 + theStructure.sigma(i - 1, j + 1, 0) / 4;
							permit = theStructure.eps(i - 1, j, nLz - 1) / 4 + theStructure.eps(i - 1, j, 0) / 4 + theStructure.eps(i - 1, j + 1, nLz - 1) / 4 + theStructure.eps(i - 1, j + 1, 0) / 4;
							loss = cond*delt / (2 * permit);

							CEyE_for_Grid[zBoundary_yInterface_Ey](i, j) = 1;
							CEyHx_for_Grid[zBoundary_yInterface_Ey](i, j) = delt / (eps0 *deltazZboundary);
							CEyHz_for_Grid[zBoundary_yInterface_Ey](i, j) = delt / (eps0 *theGrid.deltaxL(i));
						}//End if

					}//End for


					 //Needed in case structure is not periodic in y-direction

					if (!theStructure.Is_periodic_in_yDirection())
					{
						deltayYboundary = theGrid.deltayL(nLy - 1) / 2 + theGrid.deltayL(0) / 2;		//Periodicity
						deltazZboundary = theGrid.deltazL(nLz - 1) / 2 + theGrid.deltazL(0) / 2;		//Periodicity

						CExE_for_GridEdge[yBoundary_zBoundary_Ex](i) = 1;
						CExHy_for_GridEdge[yBoundary_zBoundary_Ex](i) = delt / (eps0*deltazZboundary);
						CExHz_for_GridEdge[yBoundary_zBoundary_Ex](i) = delt / (eps0*deltayYboundary);

						cond = theStructure.sigma(i, nLy - 1, nLz - 1) / 8 + theStructure.sigma(i, nLy - 1, 0) / 8;		//zero conductivity above structure and outside y-boundary
						permit = eps0 / 8 + eps0 / 8 + eps0 / 8 + eps0 / 8 + eps0 / 8 + eps0 / 8 + theStructure.eps(i, nLy - 1, nLz - 1) / 8 + theStructure.eps(i, nLy - 1, 0) / 8;
						//permitivity outside y-boundary assumed eps0
						loss = cond*delt / (2 * permit);

						CExE_for_GridEdge[yBoundary_zBoundary_xInterface_Ex](i) = (1 - loss) / (1 + loss);
						CExHz_for_GridEdge[yBoundary_zBoundary_xInterface_Ex](i) = (delt / (permit*deltayYboundary)) / (1 + loss);
						CExHy_for_GridEdge[yBoundary_zBoundary_xInterface_Ex](i) = (delt / (permit*deltazZboundary)) / (1 + loss);
					}//End if


				}//End if


				 //Compute coefficients for nodes on boundaries of structure if the structure is the unit cell of a larger periodic structure

				if (theStructure.Is_periodic_in_yDirection())
				{
					//Compute coefficients for y-boundary nodes
					for (i = 1; i < nLx - 1; i++)
						for (k = 0; k < nLz; k++)
						{
							deltayYboundary = theGrid.deltayL(nLy - 1) / 2 + theGrid.deltayL(0) / 2;		//Periodicity

							cond = theStructure.sigma(i - 1, nLy - 1, k) / 2 + theStructure.sigma(i - 1, 0, k) / 2;		//Periodicity
							permit = theStructure.eps(i - 1, nLy - 1, k) / 2 + theStructure.eps(i - 1, 0, k) / 2;		//Periodicity
							loss = cond*delt / (2 * permit);

							CEzE_for_Grid[yBoundary_Ez](i, k) = (1 - loss) / (1 + loss);
							CEzHy_for_Grid[yBoundary_Ez](i, k) = (delt / (permit*theGrid.deltaxL(i))) / (1 + loss);
							CEzHx_for_Grid[yBoundary_Ez](i, k) = (delt / (permit*deltayYboundary)) / (1 + loss);

							CExE_for_Grid[yBoundary_Ex](i, k) = (1 - loss) / (1 + loss);
							CExHy_for_Grid[yBoundary_Ex](i, k) = (delt / (permit*theGrid.deltazL(k))) / (1 + loss);
							CExHz_for_Grid[yBoundary_Ex](i, k) = (delt / (permit*deltayYboundary)) / (1 + loss);

							if (i < nLx - 2)
							{
								cond = theStructure.sigma(i - 1, nLy - 1, k) / 4 + theStructure.sigma(i - 1, 0, k) / 4 + theStructure.sigma(i, nLy - 1, k) / 4 + theStructure.sigma(i, 0, k) / 4;						//zero conductivity above structure
								permit = theStructure.eps(i - 1, nLy - 1, k) / 4 + theStructure.eps(i - 1, 0, k) / 4 + theStructure.eps(i, nLy - 1, k) / 4 + theStructure.eps(i, 0, k) / 4;
								loss = cond*delt / (2 * permit);
							}
							else
							{
								cond = theStructure.sigma(i - 1, nLy - 1, k) / 4 + theStructure.sigma(i - 1, 0, k) / 4 + theStructure.sigma(i, nLy - 1, k) / 4 + theStructure.sigma(i, 0, k) / 4;						//zero conductivity above structure
								permit = theStructure.eps(i - 1, nLy - 1, k) / 4 + theStructure.eps(i - 1, 0, k) / 4 + theStructure.eps(i, nLy - 1, k) / 4 + theStructure.eps(i, 0, k) / 4;
								loss = cond*delt / (2 * permit);
							}//End if

							CExE_for_Grid[yBoundary_xInterface_Ex](i, k) = (1 - loss) / (1 + loss);
							CExHy_for_Grid[yBoundary_xInterface_Ex](i, k) = (delt / (permit*theGrid.deltazL(k))) / (1 + loss);
							CExHz_for_Grid[yBoundary_xInterface_Ex](i, k) = (delt / (permit*deltayYboundary)) / (1 + loss);

							if (k < nLz - 1)
							{
								cond = theStructure.sigma(i - 1, nLy - 1, k) / 4 + theStructure.sigma(i - 1, 0, k) / 4 + theStructure.sigma(i - 1, nLy - 1, k + 1) / 4 + theStructure.sigma(i - 1, 0, k + 1) / 4;
								permit = theStructure.eps(i - 1, nLy - 1, k) / 4 + theStructure.eps(i - 1, 0, k) / 4 + theStructure.eps(i - 1, nLy - 1, k + 1) / 4 + theStructure.eps(i - 1, 0, k + 1) / 4;
								loss = cond*delt / (2 * permit);

								CEzE_for_Grid[yBoundary_zInterface_Ez](i, k) = (1 - loss) / (1 + loss);
								CEzHy_for_Grid[yBoundary_zInterface_Ez](i, k) = (delt / (permit*theGrid.deltaxL(i))) / (1 + loss);
								CEzHx_for_Grid[yBoundary_zInterface_Ez](i, k) = (delt / (permit*deltayYboundary)) / (1 + loss);

							}//End if

						}//End for
						 //End for

					deltazZboundary = theGrid.deltazL(nLz - 1) / 2 + theGrid.deltazL(0) / 2;		//Periodicity

					if (theStructure.Is_periodic_in_zDirection())
					{
						cond = theStructure.sigma(i - 1, nLy - 1, nLz - 1) / 4 + theStructure.sigma(i - 1, 0, nLz - 1) / 4 + theStructure.sigma(i - 1, nLy - 1, 0) / 4 + theStructure.sigma(i - 1, 0, 0) / 4;
						permit = theStructure.eps(i - 1, nLy - 1, nLz - 1) / 4 + theStructure.eps(i - 1, 0, nLz - 1) / 4 + theStructure.eps(i - 1, nLy - 1, 0) / 4 + theStructure.eps(i - 1, 0, 0) / 4;
						loss = cond*delt / (2 * permit);

						CExE_for_GridEdge[yBoundary_zBoundary_Ex](i) = (1 - loss) / (1 + loss);
						CExHy_for_GridEdge[yBoundary_zBoundary_Ex](i) = (delt / (permit*deltazZboundary)) / (1 + loss);
						CExHz_for_GridEdge[yBoundary_zBoundary_Ex](i) = (delt / (permit*deltayYboundary)) / (1 + loss);

						if (i < nLx - 2)
						{
							cond = theStructure.sigma(i - 1, nLy - 1, nLz - 1) / 8 + theStructure.sigma(i - 1, 0, nLz - 1) / 8 + theStructure.sigma(i - 1, nLy - 1, 0) / 8 + theStructure.sigma(i - 1, 0, 0) / 8
								+ theStructure.sigma(i, nLy - 1, nLz - 1) / 8 + theStructure.sigma(i, 0, nLz - 1) / 8 + theStructure.sigma(i, nLy - 1, 0) / 8 + theStructure.sigma(i, 0, 0) / 8;
							//Periodicity. 
							permit = theStructure.eps(i - 1, nLy - 1, nLz - 1) / 8 + theStructure.eps(i - 1, 0, nLz - 1) / 8 + theStructure.eps(i - 1, nLy - 1, 0) / 8 + theStructure.eps(i - 1, 0, 0) / 8
								+ theStructure.eps(i, nLy - 1, nLz - 1) / 8 + theStructure.eps(i, 0, nLz - 1) / 8 + theStructure.eps(i, nLy - 1, 0) / 8 + theStructure.eps(i, 0, 0) / 8;
							//Periodicity. 
							loss = cond*delt / (2 * permit);
						}
						else
						{
							cond = theStructure.sigma(i - 1, nLy - 1, nLz - 1) / 8 + theStructure.sigma(i - 1, 0, nLz - 1) / 8 + theStructure.sigma(i - 1, nLy - 1, 0) / 8 + theStructure.sigma(i - 1, 0, 0) / 8;
							//Periodicity. Conductivity beneath structure is zero
							permit = theStructure.eps(i - 1, nLy - 1, nLz - 1) / 8 + theStructure.eps(i - 1, 0, nLz - 1) / 8 + theStructure.eps(i - 1, nLy - 1, 0) / 8 + theStructure.eps(i - 1, 0, 0) / 8
								+ eps0 / 8 + eps0 / 8 + eps0 / 8 + eps0 / 8;
							//Periodicity. Permitivity beneath structure is eps0
							loss = cond*delt / (2 * permit);

						}//End if

						CExE_for_GridEdge[yBoundary_zBoundary_xInterface_Ex](i) = (1 - loss) / (1 + loss);
						CExHy_for_GridEdge[yBoundary_zBoundary_xInterface_Ex](i) = (delt / (permit*deltazZboundary)) / (1 + loss);
						CExHz_for_GridEdge[yBoundary_zBoundary_xInterface_Ex](i) = (delt / (permit*deltayYboundary)) / (1 + loss);


					}
					else
					{
						cond = theStructure.sigma(i - 1, nLy - 1, nLz - 1) / 4 + theStructure.sigma(i - 1, 0, nLz - 1) / 4;			//Conductivity outside z-boundary zero
						permit = theStructure.eps(i - 1, nLy - 1, nLz - 1) / 4 + theStructure.eps(i - 1, 0, nLz - 1) / 4 + eps0 / 4 + eps0 / 4;		//Permitivity outside z-boundary assumed eps0
						loss = cond*delt / (2 * permit);

						CExE_for_GridEdge[yBoundary_zBoundary_Ex](i) = (1 - loss) / (1 + loss);
						CExHy_for_GridEdge[yBoundary_zBoundary_Ex](i) = (delt / (permit*deltazZboundary)) / (1 + loss);
						CExHz_for_GridEdge[yBoundary_zBoundary_Ex](i) = (delt / (permit*deltayYboundary)) / (1 + loss);

						if (i < nLx - 2)
						{
							cond = theStructure.sigma(i - 1, nLy - 1, nLz - 1) / 8 + theStructure.sigma(i - 1, 0, nLz - 1) / 8
								+ theStructure.sigma(i, nLy - 1, nLz - 1) / 8 + theStructure.sigma(i, 0, nLz - 1) / 8;					//Conductivity outside z-boundary zero		//Periodicity. 
							permit = theStructure.eps(i - 1, nLy - 1, nLz - 1) / 8 + theStructure.eps(i - 1, 0, nLz - 1) / 8
								+ theStructure.eps(i, nLy - 1, nLz - 1) / 8 + theStructure.eps(i, 0, nLz - 1) / 8 + eps0 / 8 + eps0 / 8;	//Permitivity outside z-boundary assumed eps0 //Periodicity. 
							loss = cond*delt / (2 * permit);
						}
						else
						{
							cond = theStructure.sigma(i - 1, nLy - 1, nLz - 1) / 8 + theStructure.sigma(i - 1, 0, nLz - 1) / 8;			//Conductivity beneath structure and outside z-boundary zero. Periodicity. 
																																		//Periodicity. Conductivity beneath structure is zero
							permit = theStructure.eps(i - 1, nLy - 1, nLz - 1) / 8 + theStructure.eps(i - 1, 0, nLz - 1) / 8 + eps0 / 8 + eps0 / 8 + eps0 / 8 + eps0 / 8 + eps0 / 8 + eps0 / 8;
							//Periodicity. Permitivity beneath structure is eps0
							loss = cond*delt / (2 * permit);

						}//End if

						CExE_for_GridEdge[yBoundary_zBoundary_xInterface_Ex](i) = (1 - loss) / (1 + loss);
						CExHy_for_GridEdge[yBoundary_zBoundary_xInterface_Ex](i) = (delt / (permit*deltazZboundary)) / (1 + loss);
						CExHz_for_GridEdge[yBoundary_zBoundary_xInterface_Ex](i) = (delt / (permit*deltayYboundary)) / (1 + loss);
					}//End if

				}//End if


				if (theStructure.Is_periodic_in_zDirection())
				{
					//Compute coefficients for z-boundary nodes
					for (i = 1; i < nLx - 1; i++)
						for (j = 0; j < nLy; j++)
						{
							deltazZboundary = theGrid.deltazL(nLz - 1) / 2 + theGrid.deltazL(0) / 2;					//Periodicity

							cond = theStructure.sigma(i - 1, j, nLz - 1) / 2 + theStructure.sigma(i - 1, j, 0) / 2;		//Periodicity
							permit = theStructure.eps(i - 1, j, nLz - 1) / 2 + theStructure.eps(i - 1, j, 0) / 2;		//Periodicity
							loss = cond*delt / (2 * permit);

							CExE_for_Grid[zBoundary_Ex](i, j) = (1 - loss) / (1 + loss);
							CExHz_for_Grid[zBoundary_Ex](i, j) = (delt / (permit*theGrid.deltayL(j))) / (1 + loss);
							CExHy_for_Grid[zBoundary_Ex](i, j) = (delt / (permit*deltazZboundary)) / (1 + loss);

							if (i < nLx - 2)
							{
								cond = theStructure.sigma(i - 1, j, nLz - 1) / 4 + theStructure.sigma(i - 1, j, 0) / 4 + theStructure.sigma(i, j, nLz - 1) / 4 + theStructure.sigma(i, j, 0) / 4;
								permit = theStructure.eps(i - 1, j, nLz - 1) / 4 + theStructure.eps(i - 1, j, 0) / 4 + theStructure.eps(i, j, nLz - 1) / 4 + theStructure.eps(i, j, 0) / 4;
								loss = cond*delt / (2 * permit);

								CExE_for_Grid[zBoundary_xInterface_Ex](i, j) = (1 - loss) / (1 + loss);
								CExHz_for_Grid[zBoundary_xInterface_Ex](i, j) = (delt / (permit*theGrid.deltayL(j))) / (1 + loss);
								CExHy_for_Grid[zBoundary_xInterface_Ex](i, j) = (delt / (permit*deltazZboundary)) / (1 + loss);
							}
							else
							{
								cond = theStructure.sigma(i - 1, j, nLz - 1) / 2 + theStructure.sigma(i - 1, j, 0) / 2;		//Periodicity. Conductivity beneath structure is zero
								permit = eps0 / 4 + eps0 / 4 + theStructure.eps(i - 1, j, nLz - 1) / 4 + theStructure.eps(i - 1, j, 0) / 4;
								loss = cond*delt / (2 * permit);

								CExE_for_Grid[zBoundary_xInterface_Ex](i, j) = (1 - loss) / (1 + loss);
								CExHz_for_Grid[zBoundary_xInterface_Ex](i, j) = (delt / (permit*theGrid.deltayL(j))) / (1 + loss);
								CExHy_for_Grid[zBoundary_xInterface_Ex](i, j) = (delt / (permit*deltazZboundary)) / (1 + loss);
							}//End if

							cond = theStructure.sigma(i - 1, j, nLz - 1) / 2 + theStructure.sigma(i - 1, j, 0) / 2;		//Periodicity
							permit = theStructure.eps(i - 1, j, nLz - 1) / 2 + theStructure.eps(i - 1, j, 0) / 2;		//Periodicity
							loss = cond*delt / (2 * permit);

							CEyE_for_Grid[zBoundary_Ey](i, j) = (1 - loss) / (1 + loss);
							CEyHx_for_Grid[zBoundary_Ey](i, j) = (delt / (permit*deltazZboundary)) / (1 + loss);
							CEyHz_for_Grid[zBoundary_Ey](i, j) = (delt / (permit *theGrid.deltaxL(i))) / (1 + loss);

							if (j < nLy - 1)
							{
								cond = theStructure.sigma(i - 1, j, nLz - 1) / 4 + theStructure.sigma(i - 1, j, 0) / 4 + theStructure.sigma(i - 1, j + 1, nLz - 1) / 4 + theStructure.sigma(i - 1, j + 1, 0) / 4;
								permit = theStructure.eps(i - 1, j, nLz - 1) / 4 + theStructure.eps(i - 1, j, 0) / 4 + theStructure.eps(i - 1, j + 1, nLz - 1) / 4 + theStructure.eps(i - 1, j + 1, 0) / 4;
								loss = cond*delt / (2 * permit);

								CEyE_for_Grid[zBoundary_yInterface_Ey](i, j) = (1 - loss) / (1 + loss);
								CEyHx_for_Grid[zBoundary_yInterface_Ey](i, j) = (delt / (permit*deltazZboundary)) / (1 + loss);
								CEyHz_for_Grid[zBoundary_yInterface_Ey](i, j) = (delt / (permit *theGrid.deltaxL(i))) / (1 + loss);
							}//End if


						}//End for
						 //End for


						 //Needed in case structure is not periodic in y-direction

					if (!theStructure.Is_periodic_in_yDirection())
					{
						deltayYboundary = theGrid.deltayL(nLy - 1) / 2 + theGrid.deltayL(0) / 2;		//Periodicity
						deltazZboundary = theGrid.deltazL(nLz - 1) / 2 + theGrid.deltazL(0) / 2;		//Periodicity

						cond = theStructure.sigma(i - 1, nLy - 1, nLz - 1) / 4 + theStructure.sigma(i - 1, nLy - 1, 0) / 4;		//zero conductivity outside y-boundary
						permit = eps0 / 4 + eps0 / 4 + theStructure.eps(i, nLy - 1, nLz - 1) / 4 + theStructure.eps(i, nLy - 1, 0) / 4;
						loss = cond*delt / (2 * permit);

						CExE_for_GridEdge[yBoundary_zBoundary_Ex](i) = (1 - loss) / (1 + loss);
						CExHy_for_GridEdge[yBoundary_zBoundary_Ex](i) = (delt / (permit*deltazZboundary)) / (1 + loss);
						CExHz_for_GridEdge[yBoundary_zBoundary_Ex](i) = (delt / (permit*deltayYboundary)) / (1 + loss);

						if (i < nLx - 2)
						{
							cond = theStructure.sigma(i - 1, nLy - 1, nLz - 1) / 8 + theStructure.sigma(i - 1, nLy - 1, 0) / 8
								+ theStructure.sigma(i, nLy - 1, nLz - 1) / 8 + theStructure.sigma(i, nLy - 1, 0) / 8;		//zero conductivity outside y-boundary
							permit = theStructure.eps(i - 1, nLy - 1, nLz - 1) / 8 + theStructure.eps(i - 1, nLy - 1, 0) / 8 +
								theStructure.sigma(i, nLy - 1, nLz - 1) / 8 + theStructure.sigma(i, nLy - 1, 0) / 8 + eps0 / 8 + eps0 / 8 + eps0 / 8 + eps0 / 8;
							//permitivity outside y-boundary assumed eps0
							loss = cond*delt / (2 * permit);
						}
						else
						{
							cond = theStructure.sigma(i, nLy - 1, nLz - 1) / 8 + theStructure.sigma(i, nLy - 1, 0) / 8;		//zero conductivity outside y-boundary
							permit = theStructure.eps(i - 1, nLy - 1, nLz - 1) / 8 + theStructure.eps(i - 1, nLy - 1, 0) / 8 + eps0 / 8 + eps0 / 8 + eps0 / 8 + eps0 / 8 + eps0 / 8 + eps0 / 8;
							//permitivity outside y-boundary assumed eps0
							loss = cond*delt / (2 * permit);
						}//End if

						CExE_for_GridEdge[yBoundary_zBoundary_xInterface_Ex](i) = (1 - loss) / (1 + loss);
						CExHz_for_GridEdge[yBoundary_zBoundary_xInterface_Ex](i) = (delt / (permit*deltayYboundary)) / (1 + loss);
						CExHy_for_GridEdge[yBoundary_zBoundary_xInterface_Ex](i) = (delt / (permit*deltazZboundary)) / (1 + loss);
					}//End if

				}//End if



				 //Compute coefficients for nodes on boundaries of region below structure if the structure is the unit cell of a larger periodic structure
				i = nLx - 1;

				if (theStructure.Is_periodic_in_yDirection())
				{
					//Compute coefficients for y-boundary nodes
					for (k = 0; k < nLz; k++)
					{
						deltayYboundary = theGrid.deltayL(nLy - 1) / 2 + theGrid.deltayL(0) / 2;		//Periodicity

						CEzE_for_Grid[yBoundary_Ez](i, k) = 1;
						CEzHy_for_Grid[yBoundary_Ez](i, k) = delt / (eps0*theGrid.deltaxL(i));
						CEzHx_for_Grid[yBoundary_Ez](i, k) = delt / (eps0*deltayYboundary);

						CEzE_for_Grid[yBoundary_zInterface_Ez](i, k) = 1;						//not used if k = nLz - 1
						CEzHy_for_Grid[yBoundary_zInterface_Ez](i, k) = CEzHy_for_Grid[yBoundary_Ez](i, k);	//not used if k = nLz - 1
						CEzHx_for_Grid[yBoundary_zInterface_Ez](i, k) = CEzHx_for_Grid[yBoundary_Ez](i, k);	//not used if k = nLz - 1

						CExE_for_Grid[yBoundary_Ex](i, k) = 1;
						CExHy_for_Grid[yBoundary_Ex](i, k) = delt / (eps0*theGrid.deltazL(k));
						CExHz_for_Grid[yBoundary_Ex](i, k) = delt / (eps0*deltayYboundary);
					}//End for

					deltazZboundary = theGrid.deltazL(nLz - 1) / 2 + theGrid.deltazL(0) / 2;		//Periodicity

					CExE_for_GridEdge[yBoundary_zBoundary_Ex](i) = 1;
					CExHy_for_GridEdge[yBoundary_zBoundary_Ex](i) = delt / (eps0*deltazZboundary);
					CExHz_for_GridEdge[yBoundary_zBoundary_Ex](i) = delt / (eps0*deltayYboundary);

				}//End if

				if (theStructure.Is_periodic_in_zDirection())
				{
					//Compute coefficients for z-boundary nodes
					for (j = 0; j < nLy; j++)
					{
						deltazZboundary = theGrid.deltazL(nLz - 1) / 2 + theGrid.deltazL(0) / 2;		//Periodicity

						CExE_for_Grid[zBoundary_Ex](i, j) = 1;
						CExHz_for_Grid[zBoundary_Ex](i, j) = delt / (eps0*theGrid.deltayL(j));
						CExHy_for_Grid[zBoundary_Ex](i, j) = delt / (eps0*deltazZboundary);

						CEyE_for_Grid[zBoundary_Ey](i, j) = 1;
						CEyHx_for_Grid[zBoundary_Ey](i, j) = delt / (eps0 *deltazZboundary);
						CEyHz_for_Grid[zBoundary_Ey](i, j) = delt / (eps0 *theGrid.deltaxL(i));

						if (j < nLy - 1)
						{
							cond = theStructure.sigma(i - 1, j, nLz - 1) / 4 + theStructure.sigma(i - 1, j, 0) / 4 + theStructure.sigma(i - 1, j + 1, nLz - 1) / 4 + theStructure.sigma(i - 1, j + 1, 0) / 4;
							permit = theStructure.eps(i - 1, j, nLz - 1) / 4 + theStructure.eps(i - 1, j, 0) / 4 + theStructure.eps(i - 1, j + 1, nLz - 1) / 4 + theStructure.eps(i - 1, j + 1, 0) / 4;
							loss = cond*delt / (2 * permit);

							CEyE_for_Grid[zBoundary_yInterface_Ey](i, j) = 1;
							CEyHx_for_Grid[zBoundary_yInterface_Ey](i, j) = delt / (eps0 *deltazZboundary);
							CEyHz_for_Grid[zBoundary_yInterface_Ey](i, j) = delt / (eps0 *theGrid.deltaxL(i));
						}//End if

					}//End for


					 //Needed in case structure is not periodic in y-direction

					if (!theStructure.Is_periodic_in_yDirection())
					{
						deltayYboundary = theGrid.deltayL(nLy - 1) / 2 + theGrid.deltayL(0) / 2;		//Periodicity
						deltazZboundary = theGrid.deltazL(nLz - 1) / 2 + theGrid.deltazL(0) / 2;		//Periodicity

						CExE_for_GridEdge[yBoundary_zBoundary_Ex](i) = 1;
						CExHy_for_GridEdge[yBoundary_zBoundary_Ex](i) = delt / (eps0*deltazZboundary);
						CExHz_for_GridEdge[yBoundary_zBoundary_Ex](i) = delt / (eps0*deltayYboundary);

					}//End if

					 //No x-interface values neded beneath structurre

				}//End if


			}//End function

			*/

}//End namespace FDTD