/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** Fire & Heat Transfer modules developed by:                         **
**   Yaqiang Jiang (y.jiang@ed.ac.uk)                                 **
**   Asif Usmani (asif.usmani@ed.ac.uk)                               **
**                                                                    **
** ****************************************************************** */

//
// Template written by Yaqiang Jiang (y.jiang@ed.ac.uk)
//
// AluminumEC9 written by Mohamed Anwar Orabi (maorabi@polyu.edu.hk)
//


//
//
// The properties of aluminum are only known until about 500 C, beyond which
// the Eurocode does not provide details. This is often reached in fire situations 
// and as such we had to go beyond the codified properties. In general, when
// tempearture reaches above 500, the properties are taken as if they were at 500.
//
//

#include <AluminumEC9.h>
#include <Matrix.h>
#include <OPS_Globals.h>
#include <cmath>

double AluminumEC9::epsilon = 1e-5;

AluminumEC9::AluminumEC9(int tag)
:HeatTransferMaterial(tag), trial_temp(0.0), 
 ini_temp(0.0), rho(2700.0), cp(0.0), enthalpy(0.0)
{
    if ( k == 0){
		k = new Matrix(3,3);
		if (k == 0) {
			opserr << "AluminumEC9::CarbonSteelEN93() - out of memory\n";
			exit(-1);
			}
		}
}


AluminumEC9::~AluminumEC9()
{
    if (k != 0)
		delete k;
}

int 
AluminumEC9::setTrialTemperature(double temp, int par)
{
    trial_temp = temp - 273.15;
    return 0;
}


const Matrix& 
AluminumEC9::getConductivity(void)
{
    if (trial_temp < 20.0) {
		// if ((20.0 - trial_temp) < epsilon) {
		//if (int(trial_temp + 0.5) == 20) {
			//double kc = 54.0 - 3.33 * 1e-2 * trial_temp;
			(*k)(0,0) = 191.4;
			(*k)(1,1) = 191.4;
			(*k)(2,2) = 191.4;
			//} else {
			//opserr << "CarbonSteelEN93::getConductivity() - trial temperature " 
			//	<< trial_temp << " out of bounds\n";
			//exit(-1);}
		} else if ((20.0 <= trial_temp) && (trial_temp < 500.0)) {
			double kc = 0.07* trial_temp + 190;
			(*k)(0,0) = kc;
			(*k)(1,1) = kc;
			(*k)(2,2) = kc;
		}else if ((500.0 <= trial_temp) && (trial_temp <= 1200.0)) {
			(*k)(0,0) = 225.0;
			(*k)(1,1) = 225.0;
			(*k)(2,2) = 225.0;
			} else {
				opserr << "AluminumEC9::getConductivity() - trial temperature "
					<< trial_temp << " out of bounds\n";
				exit(-1);
		 }

    return *k; // change
}


double  
AluminumEC9::getRho(void)
{
    return rho;
}


double 
AluminumEC9::getSpecificHeat(void)
{
    if (trial_temp < 20.0){
		cp = 911.2;
	} else if ((20.0 <= trial_temp) && (trial_temp < 500.0)) {
		cp = trial_temp * 0.41 + 903;
	} else {
		cp = 1108.0;
		}

    return cp;
}


double
AluminumEC9::getEnthalpy()
{
    // * The enthalpy data here is obtained by integrating the spesific heat over
    // * the temperature range [0,500] given in EN 1999-1-2.
	// * The analytical expresion is easily calcualted because the heat capacity is linear.
	// * because a limit of 500 is very small, enthalpy is returned assuming a max temp of 500
	if ((trial_temp >= 0.0) && (trial_temp <= 500))
	{
		enthalpy = 0.205 * trial_temp * trial_temp + 903 * trial_temp;
	}
	else if ((trial_temp > 500.0) && (trial_temp <= 1200.0))
	{
		enthalpy = 0.205 * 500.0 * 500.0 + 903.0 * 500.0;
	}
	else
	{
		opserr << "AluminumEC9::getEnthalpy(double temp) - nodal temperature "
			<< trial_temp << " out of bounds of [0, 1200]\n";
		exit(-1);
	}
			

    return enthalpy;	
}


double
AluminumEC9::getEnthalpy(double temp)
{
    double enthp;
	double nod_temp = temp - 273.15;
    
	// The temperature range is expanded to [0,1200] rather than the original one [20,1200] in Eurocode.
	// The reason is, for an analysis with initial temperature at 20, the solution could be lower than
	// 20 after initial iterations. Eventhough, the slope of H-T within the expanded temperature range is 
	// kept constant, the same as the heat capacity at T = 20;

	if ((nod_temp >= 0.0) && (nod_temp <= 500))
	{
		enthp = 0.205 * nod_temp * nod_temp + 903 * nod_temp;
	}
	else if ((nod_temp > 500.0) && (nod_temp <= 1200.0))
	{
		enthp = 0.205 * 500.0 * 500.0 + 903.0 * 500.0;
	}
	else
	{
		opserr << "AluminumEC9::getEnthalpy() - nodal temperature "
			<< nod_temp << " out of bounds of [0, 1200]\n";
		exit(-1);
	}
    return enthp;	
}


HeatTransferMaterial*
AluminumEC9::getCopy(void)
{
    AluminumEC9* theCopy = new AluminumEC9(this->getTag());
    theCopy->trial_temp = trial_temp;
    return theCopy;
}


void
AluminumEC9::update()
{
    return; 
}


int 
AluminumEC9::commitState(void)
{
    return 0;
}


int 
AluminumEC9::revertToLastCommit(void)
{
    return 0;
}


int 
AluminumEC9::revertToStart(void)
{
    return 0;
}


