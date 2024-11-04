/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
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
** ****************************************************************** */
                                                                        
//Modified by:  Anand Kumar IITJ---------02,2023// 
//              


#include <DPMsteelEC.h>
#include <Vector.h>
#include <Matrix.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>

#include <math.h>
#include <float.h>
#include <complex.h>

#include <elementAPI.h>
#include <OPS_Globals.h>


void *
OPS_DPMsteelEC()
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int    iData[1];
  double dData[7];
  int numData = 1;

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial DPMsteelEC tag" << endln;
    return 0;
  }

  numData = OPS_GetNumRemainingInputArgs();

  if (numData != 3 && numData != 7) {
    opserr << "Invalid #args, want: uniaxialMaterial DPMsteelEC " << iData[0] << " fy? E? b? <a1? a2? a3? a4?>>" << endln;
    return 0;
  }

  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid #args, want: uniaxialMaterial DPMsteelEC " << iData[0] << " fy? E? b? <a1? a2? a3? a4?>>" << endln;
    return 0;
  }

  if (numData == 3) {
    dData[3] = STEEL_01_DEFAULT_A1;
    dData[4] = STEEL_01_DEFAULT_A2;
    dData[5] = STEEL_01_DEFAULT_A3;
    dData[6] = STEEL_01_DEFAULT_A4;
  }

  // Parsing was successful, allocate the material
  theMaterial = new DPMsteelEC(iData[0], dData[0], dData[1], 
				   dData[2], dData[3], dData[4], 
				   dData[5], dData[6]);
  
  
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type DPMsteelEC Material\n";
    return 0;
  }

  return theMaterial;
}


DPMsteelEC::DPMsteelEC
(int tag, double FY, double E, double B,
    double A1, double A2, double A3, double A4) :
    UniaxialMaterial(tag, MAT_TAG_DPMsteelEC),
    fyT(FY), E0T(E), b(B), a1(A1), a2(A2), a3(A3), a4(A4)
{
    // Sets all history and state variables to initial values
    // History variables
    CminStrain = 0.0;
    CmaxStrain = 0.0;
    CshiftP = 1.0;
    CshiftN = 1.0;
    Cloading = 0;
    Ctemperature = 0;    //Added by Liming,2013
    Cmono = true;    //Added by Liming,2013

    TminStrain = 0.0;
    TmaxStrain = 0.0;
    TshiftP = 1.0;
    TshiftN = 1.0;
    Tloading = 0;
    Ttemperature = 0;    //Added by Liming,2013
    Tmono = true;    //Added by Liming,2013

    // State variables
    Cstrain = 0.0;
    Cstress = 0.0;
    //Ctangent = E0;

    Ctangent = E0T;///JZ, 07/10//
    Cmodulus = E0T; //added by Princeton


    Tstrain = 0.0;
    Tstress = 0.0;
    // Ttangent = E0;

    Ttangent = E0T;///JZ, 07/10//
    Tmodulus = E0T; //added by Princeton

    // AddingSensitivity:BEGIN /////////////////////////////////////
    parameterID = 0;
    SHVs = 0;
    // AddingSensitivity:END //////////////////////////////////////

    ThermalElongation = 0; //initialize //JZ, 07/10//
    E0 = E0T;//JZ, 07/10//
    fy = fyT;//JZ, 07/10//
    fp = 0.0;//JZ, 11/10//
    fu = 0.0; // Added by Anand Kumar IITJ  
    eplasP = 0.0; //Added by Anand Kumar IITJ : Previous Plastic Strain  
    eplas = 0.0;  //Added by Anand Kumar IITJ : Current Plastic Srain
    alpP = 0.0;   //Added by Anand Kumar IITJ : Previous alpha
    alp = 0.0;    //Added by Anand Kumar IITJ : Current alpha
}

DPMsteelEC::DPMsteelEC():UniaxialMaterial(0,MAT_TAG_DPMsteelEC),
 //fy(0.0), E0(0.0), b(0.0), a1(0.0), a2(0.0), a3(0.0), a4(0.0)
 fyT(0.0), E0T(0.0), b(0.0), a1(0.0), a2(0.0), a3(0.0), a4(0.0) //JZ, 07/10//
{
    // AddingSensitivity:BEGIN /////////////////////////////////////
    parameterID = 0;
    SHVs = 0;
    // AddingSensitivity:END //////////////////////////////////////

	  ThermalElongation = 0; //initialize //JZ, 07/10//
	  E0 = E0T;//JZ, 07/10//
	  fy = fyT;//JZ, 07/10//
	  fp = 0;//JZ, 11/10//
	  Ttemperature = 0;  //Added by Liming,2013
	  Ctemperature = 0;  //Added by Liming,2013
      Cmono = true;      //Added by Liming,2013
	  Tmono = true;      //Added by Liming,2013
}

DPMsteelEC::~DPMsteelEC()
{
    // AddingSensitivity:BEGIN /////////////////////////////////////
    if (SHVs != 0)
        delete SHVs;
    // AddingSensitivity:END //////////////////////////////////////
}

int DPMsteelEC::setTrialStrain(double strain, double FiberTemperature, double strainRate)
{

  Ttemperature  = FiberTemperature;    //Added by Liming,2013


   // Reset history variables to last converged state
   TminStrain = CminStrain;
   TmaxStrain = CmaxStrain;
   TshiftP = CshiftP;
   TshiftN = CshiftN;
   Tloading = Cloading;
   Tstrain = Cstrain;
   Tstress = Cstress;
   Ttangent = Ctangent;
   Tmono = Cmono;       //Added by Liming,2013
   // Determine change in strain from last converged state
   double dStrain = strain - Cstrain;
   
   if (fabs(dStrain) > DBL_EPSILON || (fabs(Ttemperature-Ctemperature)>1e-5)) {
     // Set trial strain
     Tstrain = strain;
     // Calculate the trial state given the trial strain 
     determineTrialState (dStrain);
	 
     Ctemperature = Ttemperature;    //Added by Liming,2013
	 // To enable 0 dStrain in the first iteration for new thermal action step
	 // The material stress may change due to the degraded modulus.---added by Liming, 2013, OCT//
     
     
   }

   return 0;
}

int DPMsteelEC::setTrial (double strain, double &stress, double &tangent, double strainRate)
{
   // Reset history variables to last converged state
   TminStrain = CminStrain;
   TmaxStrain = CmaxStrain;
   TshiftP = CshiftP;
   TshiftN = CshiftN;
   Tloading = Cloading;
   Tstrain = Cstrain;
   Tstress = Cstress;
   Ttangent = Ctangent;

   // Determine change in strain from last converged state
   double dStrain = strain - Cstrain;

   if (fabs(dStrain) > DBL_EPSILON) {
     // Set trial strain
     Tstrain = strain;
     
     // Calculate the trial state given the trial strain
     determineTrialState (dStrain);

   }

   stress = Tstress;
   tangent = Ttangent;

   return 0;
}

void DPMsteelEC::determineTrialState (double dStrain)
{
      double fyOneMinusB = fy * (1.0 - b);

      double Esh = b*E0;
      double epsy = fy/E0;
      
      double c1 = Esh*Tstrain;
      
      double c2 = TshiftN*fyOneMinusB;

      double c3 = TshiftP*fyOneMinusB;
	  
      double c = Cstress + E0*dStrain;
	  
	  //----if E0 changed due to temperature rise and the material 
	  //----has not come into unloading, the rectant stress would not be c
	  //--------------------------------------------------------------
      if (Tmono&&(E0!=E0T)) 
		   c = E0*Tstrain;
	 
	  //---------------
	  //Modified by Liming for considering stiffness degradation(E0).

      /**********************************************************
         removal of the following lines due to problems with
	 optimization may be required (e.g. on gnucc compiler
         with optimization turned on & -ffloat-store option not
         used) .. replace them with line that follows but which 
         now requires 2 function calls to achieve same result !!
      ************************************************************/

      double c1c3 = c1 + c3;

      if (c1c3 < c)
	Tstress = c1c3;
      else
	Tstress = c;

      double c1c2 = c1-c2;

      if (c1c2 > Tstress)
	Tstress = c1c2;

      /* ***********************************************************
      and replace them with:

      Tstress = fmax((c1-c2), fmin((c1+c3),c));
      **************************************************************/

      if (fabs(Tstress-c) < DBL_EPSILON)
	  Ttangent = E0;
      else
	  Ttangent = Esh;
      
      //Anand Kumar IITJ code starts:
      if ((E0 != E0T)) {
          double e = Tstrain;     
          //double epl = Cstrain - Cstress / (Cmodulus);   
          double ep = (fp / E0);
          double ey = 0.02;
          
          int sign = 0;
          if (e < 0) {
              sign = -1;
          }
          else {
              sign = 1;
          }          
          double cfun = (pow((fy - fp), 2)) / ((ey - ep) * E0 - 2 * (fy - fp));
          double bfun = sqrt(cfun * (ey - ep) * E0 + pow(cfun, 2));
          double afun = sqrt((ey - ep) * (ey - ep + cfun / E0));
  
          double sig_trial = E0 * (e - eplasP);
          double r, G;
          this->yield_funtion(E0, afun, bfun, cfun, alpP, fp, fy, fu, ey, r, G);
          if ((fabs(sig_trial) - r) < 0) {
              c = sig_trial;
              Ttangent = E0;
              eplas = eplasP;
              alp = alpP;
          }
          else {
              double sig_trial1 = E0 * (e - eplasP);
              double stress11;
              double E_tgt11;
              double eplas11;
              double alp11;
              this->Nr(sig_trial1, E0, afun, bfun, cfun, eplasP, fp, fy, fu, ey, alpP, stress11, E_tgt11, eplas11, alp11);
              c = stress11;
              if (E_tgt11 <= 10.0) {
                  Ttangent = 10.0;
              }
              else {
                  Ttangent = E_tgt11;
              }               
              eplas = eplas11;
              alp = alp11;
          }

          Tstress = c;          
          Tmodulus = E0;  
      }
      // Anand Kumar IITJ code ends.
  
      if (Tloading == 0 && dStrain != 0.0) {
	  if (dStrain > 0.0)
	    Tloading = 1;
	  else
	    Tloading = -1;
      }

      // Transition from loading to unloading, i.e. positive strain increment
      // to negative strain increment
      if (Tloading == 1 && dStrain < 0.0) {
	  Tloading = -1;
	  if (Tmono&& (fabs(Ttemperature - Ctemperature)<1e-5))
		  Tmono = false;
	  
	  if (Cstrain > TmaxStrain)
	    TmaxStrain = Cstrain;
	  TshiftN = 1 + a1*pow((TmaxStrain-TminStrain)/(2.0*a2*epsy),0.8);
      }

      // Transition from unloading to loading, i.e. negative strain increment
      // to positive strain increment
      if (Tloading == -1 && dStrain > 0.0) {
	  Tloading = 1;
	  //-------------------
	  if (Tmono&& (fabs(Ttemperature-Ctemperature)<1e-5))
		  Tmono = false;
	  //added by liming for identifying the material loading status
	  if (Cstrain < TminStrain)
	    TminStrain = Cstrain;
	  TshiftP = 1 + a3*pow((TmaxStrain-TminStrain)/(2.0*a4*epsy),0.8);
      }
   
}


double DPMsteelEC::getStrain ()
{
   return Tstrain;
}

double DPMsteelEC::getStress ()
{
   return Tstress;
}

double DPMsteelEC::getTangent ()
{
   return Ttangent; 
}

double 
DPMsteelEC::getThermalElongation(void) //***JZ
{
  return ThermalElongation;
}

//JZ 07/10 /////////////////////////////////////////////////////////////start
double 
DPMsteelEC::getElongTangent(double TempT, double &ET, double &Elong, double TempTmax) //PK add to include max temp
{
  //JZ updated, from rebar to C steel
  
  // EN 1992 pt 1-2-1. Class N hot rolled  reinforcing steel at elevated temperatures   
    if (TempT <= 80) {
        fy = fyT;
        E0 = E0T;

        //b=TempT*0.00325/80;

        fp = fyT;
        fu =  1.25 * fy; //Added by Anand Kumar
    }
    else if (TempT <= 180) {
        fy = fyT;
        E0 = E0T * (1 - (TempT - 80) * 0.1 / 100);

        //b=0.00325+(TempT - 80)*0.00325/100;

        fp = fyT * (1 - (TempT - 80) * (1 - 0.807) / 100);
        fu = 1.25 * fy; //Added by Anand Kumar
    }
    else if (TempT <= 280) {
        fy = fyT;
        E0 = E0T * (0.9 - (TempT - 180) * 0.1 / 100);

        //b=0.0065+(TempT - 180)*0.00325/100;

        fp = fyT * (0.807 - (TempT - 180) * (0.807 - 0.613) / 100);
        fu = 1.25 * fy; //Added by Anand Kumar
    }
    else if (TempT <= 380) {
        fy = fyT;
        E0 = E0T * (0.8 - (TempT - 280) * 0.1 / 100);

        //b=0.00975+(TempT - 280)*0.00355/100;

        fp = fyT * (0.613 - (TempT - 280) * (0.613 - 0.42) / 100);
        fu = (2 - 0.0025 * TempT) * fy; //Added by Anand Kumar
        //fu = 1.01 * fy;
    }
    else if (TempT <= 480) {
        fy = fyT * (1 - (TempT - 380) * 0.22 / 100);
        E0 = E0T * (0.7 - (TempT - 380) * 0.1 / 100);

        //b=0.0133+(TempT - 380)*0.0133/100;

        fp = fyT * (0.42 - (TempT - 380) * (0.42 - 0.36) / 100);
        fu = 1.01 * fy; //Added by Anand Kumar
    }
    else if (TempT <= 580) {
        fy = fyT * (0.78 - (TempT - 480) * 0.31 / 100);
        E0 = E0T * (0.6 - (TempT - 480) * 0.29 / 100);

        //b=0.0266+(TempT - 480)*0.0136/100;

        fp = fyT * (0.36 - (TempT - 480) * (0.36 - 0.18) / 100);
        fu = 1.01 * fy; //Added by Anand Kumar
    }
    else if (TempT <= 680) {
        fy = fyT * (0.47 - (TempT - 580) * 0.24 / 100);
        E0 = E0T * (0.31 - (TempT - 580) * 0.18 / 100);

        // b=0.0402-(TempT - 580)*0.0067/100;

        fp = fyT * (0.18 - (TempT - 580) * (0.18 - 0.075) / 100);
        fu = 1.01 * fy; //Added by Anand Kumar
    }
    else if (TempT <= 780) {
        fy = fyT * (0.23 - (TempT - 680) * 0.12 / 100);
        E0 = E0T * (0.13 - (TempT - 680) * 0.04 / 100);

        // b=0.0335-(TempT - 680)*0.0067/100;

        fp = fyT * (0.075 - (TempT - 680) * (0.075 - 0.05) / 100); 
        fu = 1.01 * fy; //Added by Anand Kumar
    }
    else if (TempT <= 880) {
        fy = fyT * (0.11 - (TempT - 780) * 0.05 / 100);
        E0 = E0T * (0.09 - (TempT - 780) * 0.0225 / 100);

        //  b=0.0268-(TempT - 780)*0.0067/100;

        fp = fyT * (0.05 - (TempT - 780) * (0.05 - 0.0375) / 100);
        fu = 1.01 * fy; //Added by Anand Kumar
    }
    else if (TempT <= 980) {
        fy = fyT * (0.06 - (TempT - 880) * 0.02 / 100);
        E0 = E0T * (0.0675 - (TempT - 880) * (0.0675 - 0.045) / 100);

        //  b=0.0201-(TempT - 880)*0.0067/100;

        fp = fyT * (0.0375 - (TempT - 880) * (0.0375 - 0.025) / 100);
        fu = 1.01 * fy; //Added by Anand Kumar
    }
    else if (TempT <= 1080) {
        fy = fyT * (0.04 - (TempT - 980) * 0.02 / 100);
        E0 = E0T * (0.045 - (TempT - 980) * (0.045 - 0.0225) / 100);

        // b=0.0134-(TempT - 980)*0.0067/100;

        fp = fyT * (0.025 - (TempT - 980) * (0.025 - 0.0125) / 100);
        fu = 1.01 * fy; //Added by Anand Kumar
    }
    else if (TempT <= 1180) {
        fy = fyT * (0.02 - (TempT - 1080) * 0.02 / 100);
        E0 = E0T * (0.0225 - (TempT - 1080) * 0.0225 / 100);

        //  b=0.0067-(TempT - 980)*0.0067/100;

        fp = fyT * (0.0125 - (TempT - 1080) * 0.0125 / 100);
        fu = 1.01 * fy; //Added by Anand Kumar
    }
    else {
        opserr << "the temperature is invalid23\n";
    }

  // caculation of thermal elongation of reinforcing steel. JZ
///*	
  if (TempT <= 1) {
		  ThermalElongation = TempT * 1.2164e-5;
	  }
  else if (TempT <= 730) {
      ThermalElongation = -2.416e-4 + 1.2e-5 *(TempT+20) + 0.4e-8 *(TempT+20)*(TempT+20);
  }
  else if (TempT <= 840) {
      ThermalElongation = 11e-3;
  }
  else if (TempT <= 1180) {
      ThermalElongation = -6.2e-3 + 2e-5*(TempT+20);
  }
  else {
	  opserr << "the temperature is invalid25\n";
  }

  //ThermalElongation = 0 ;   //debug  Liming
  ET = E0;   
  Elong = ThermalElongation;
  //TemperautreC = TempT;  //removed by liming

  //opserr << "\getelongation: " << ET << "\ temp:" << TemperautreC <<endln; //PK Check

  return 0;
}
//JZ 07/10 /////////////////////////////////////////////////////////////end 


int DPMsteelEC::commitState ()
{
   // History variables
   CminStrain = TminStrain;
   CmaxStrain = TmaxStrain;
   CshiftP = TshiftP;
   CshiftN = TshiftN;
   Cloading = Tloading;
   Cmono = Tmono;//added by liming,2013

   // State variables
   Cstrain = Tstrain;
   Cstress = Tstress;
   Ctangent = Ttangent;
   Cmodulus = Tmodulus; //added by Anand Kumar IITJ
   eplasP = eplas; //added by Anand Kumar IITJ
   alpP = alp; //added by Anand Kumar IITJ

   return 0;
}

int DPMsteelEC::revertToLastCommit ()
{
   // Reset trial history variables to last committed state
   TminStrain = CminStrain;
   TmaxStrain = CmaxStrain;
   TshiftP = CshiftP;
   TshiftN = CshiftN;
   Tloading = Cloading;
   Tmono = Cmono;
   
   // Reset trial state variables to last committed state
   Tstrain = Cstrain;
   Tstress = Cstress;
   Ttangent = Ctangent;
   Tmodulus = Cmodulus;  //added by Anand Kumar IITJ
   eplas = eplasP; //added by Anand Kumar IITJ
   alp = alpP; //added by Anand Kumar IITJ


   // AddingSensitivity:BEGIN /////////////////////////////////
   if (SHVs != 0)
       SHVs->Zero();
   // AddingSensitivity:END //////////////////////////////////

   return 0;
}

int DPMsteelEC::revertToStart ()
{
   // History variables
   CminStrain = 0.0;
   CmaxStrain = 0.0;
   CshiftP = 1.0;
   CshiftN = 1.0;
   Cloading = 0;

   TminStrain = 0.0;
   TmaxStrain = 0.0;
   TshiftP = 1.0;
   TshiftN = 1.0;
   Tloading = 0;

   // State variables
   Cstrain = 0.0;
   Cstress = 0.0;
   Ctangent = E0T;
   Cmodulus = E0;  //added by Anand Kumar IITJ

   Tstrain = 0.0;
   Tstress = 0.0;
   Ttangent = E0T;

   return 0;
}

UniaxialMaterial* DPMsteelEC::getCopy ()
{
   DPMsteelEC* theCopy = new DPMsteelEC(this->getTag(), fy, E0, b,
				  a1, a2, a3, a4);

   // Converged history variables
   theCopy->CminStrain = CminStrain;
   theCopy->CmaxStrain = CmaxStrain;
   theCopy->CshiftP = CshiftP;
   theCopy->CshiftN = CshiftN;
   theCopy->Cloading = Cloading;
   theCopy->Cmono = Cmono;   //added by liming,2013
   
   
   // Trial history variables
   theCopy->TminStrain = TminStrain;
   theCopy->TmaxStrain = TmaxStrain;
   theCopy->TshiftP = TshiftP;
   theCopy->TshiftN = TshiftN;
   theCopy->Tloading = Tloading;
   theCopy->Cmono = Tmono;   //added by liming,2013
   
   
   // Converged state variables
   theCopy->Cstrain = Cstrain;
   theCopy->Cstress = Cstress;
   theCopy->Ctangent = Ctangent;
   theCopy->Ctemperature = Ctemperature;   theCopy->Cmodulus = Cmodulus; //added by Anand Kumar IITJ  //added by liming,2013
   
   // Trial state variables
   theCopy->Tstrain = Tstrain;
   theCopy->Tstress = Tstress;
   theCopy->Ttangent = Ttangent;
   theCopy->Ctemperature = Ttemperature;   theCopy->Tmodulus = Tmodulus; //added by Anand Kumar IITJ  //added by liming,2013
   return theCopy;
}

int DPMsteelEC::sendSelf (int commitTag, Channel& theChannel)
{
   int res = 0;
   static Vector data(16);
   data(0) = this->getTag();

   // Material properties
   data(1) = fy;
   data(2) = E0;
   data(3) = b;
   data(4) = a1;
   data(5) = a2;
   data(6) = a3;
   data(7) = a4;

   // History variables from last converged state
   data(8) = CminStrain;
   data(9) = CmaxStrain;
   data(10) = CshiftP;
   data(11) = CshiftN;
   data(12) = Cloading;

   // State variables from last converged state
   data(13) = Cstrain;
   data(14) = Cstress;
   data(15) = Ctangent;

   // Data is only sent after convergence, so no trial variables
   // need to be sent through data vector

   res = theChannel.sendVector(this->getDbTag(), commitTag, data);
   if (res < 0) 
      opserr << "DPMsteelEC::sendSelf() - failed to send data\n";

   return res;
}

int DPMsteelEC::recvSelf (int commitTag, Channel& theChannel,
                                FEM_ObjectBroker& theBroker)
{
   int res = 0;
   static Vector data(16);
   res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  
   if (res < 0) {
      opserr << "DPMsteelEC::recvSelf() - failed to receive data\n";
      this->setTag(0);      
   }
   else {
      this->setTag(int(data(0)));

      // Material properties
      fy = data(1);
      E0 = data(2);
      b = data(3);
      a1 = data(4);
      a2 = data(5);
      a3 = data(6);
      a4 = data(7);

      // History variables from last converged state
      CminStrain = data(8);
      CmaxStrain = data(9);
      CshiftP = data(10);
      CshiftN = data(11);
      Cloading = int(data(12));

      // Copy converged history values into trial values since data is only
      // sent (received) after convergence
      TminStrain = CminStrain;
      TmaxStrain = CmaxStrain;
      TshiftP = CshiftP;
      TshiftN = CshiftN;
      Tloading = Cloading;

      // State variables from last converged state
      Cstrain = data(13);
      Cstress = data(14);
      Ctangent = data(15);      

      // Copy converged state values into trial values
      Tstrain = Cstrain;
      Tstress = Cstress;
      Ttangent = Ctangent;
      Tmodulus = Cmodulus; //added by Anand Kumar IITJ
   }
    
   return res;
}

void DPMsteelEC::Print (OPS_Stream& s, int flag)
{
   s << "DPMsteelEC tag: " << this->getTag() << endln;
   s << "  fy: " << fy << " ";
   s << "  E0: " << E0 << " ";
   s << "  b:  " << b << " ";
   s << "  a1: " << a1 << " ";
   s << "  a2: " << a2 << " ";
   s << "  a3: " << a3 << " ";
   s << "  a4: " << a4 << " ";
}


///// added by Anand Kumar IITJ
void
DPMsteelEC::yield_funtion(double E2, double a2, double b2, double c2, double epp, double fp2, double fy2, double fu2, double ey2, double &sigy, double &G)
{ 
    if (epp <= (0.04 - fu2 / E2)) {
        double q = (E2 * E2 * a2 * a2 - E2 * E2 * epp * epp + 2 * E2 * E2 * epp * ey2 - E2 * E2 * ey2 * ey2 + 2 * E2 * c2 * epp - 2 * E2 * c2 * ey2 - 2 * E2 * epp * fp2 + 2 * E2 * ey2 * fp2 + b2 * b2 - c2 * c2 + 2 * c2 * fp2 - fp2 * fp2);

        if (q <= 0) {
            q = 0.0;
        }

        sigy = fp2 + ((E2 * (b2 * b2 * ey2 - b2 * b2 * epp - E2 * a2 * a2 * c2 + E2 * a2 * a2 * fp2 + a2 * b2 * sqrt(q))) / (E2 * E2 * a2 * a2 + b2 * b2) - fp2) * (1 / (1 + exp(1e10 * (epp - 0.02 + fy2 / E2)))) + (-fp2 + (50 * (fu2 - fy2) * epp + 2 * fy2 - fu2) * (1 / (1 - 50 * (fu2 - fy2) / E2))) * (1 / (1 + exp(1e10 * -(epp - 0.02 + fy2 / E2)))) + ((epp - 0.04 + fu2 / E2) * -(50 * fu2 - 50 * fy2) / (fy2 / 4000 - fu2 / 4000 + 1)) * (1 / (1 + exp(1e10 * -(epp - 0.04 + fu2 / E2))));

        double ep_p = epp + 1e-9;

        double q1 = (E2 * E2 * a2 * a2 - E2 * E2 * ep_p * ep_p + 2 * E2 * E2 * ep_p * ey2 - E2 * E2 * ey2 * ey2 + 2 * E2 * c2 * ep_p - 2 * E2 * c2 * ey2 - 2 * E2 * ep_p * fp2 + 2 * E2 * ey2 * fp2 + b2 * b2 - c2 * c2 + 2 * c2 * fp2 - fp2 * fp2);
        if (q1 <= 0) {
            q1 = 0;
        }

        double sigy1 = fp2 + ((E2 * (b2 * b2 * ey2 - b2 * b2 * ep_p - E2 * a2 * a2 * c2 + E2 * a2 * a2 * fp2 + a2 * b2 * sqrt(q1))) / (E2 * E2 * a2 * a2 + b2 * b2) - fp2) * (1 / (1 + exp(1e10 * (ep_p - 0.02 + fy2 / E2)))) + (-fp2 + (50 * (fu2 - fy2) * ep_p + 2 * fy2 - fu2) * (1 / (1 - 50 * (fu2 - fy2) / E2))) * (1 / (1 + exp(1e10 * -(ep_p - 0.02 + fy2 / E2)))) + ((ep_p - 0.04 + fu2 / E2) * -(50 * fu2 - 50 * fy2) / (fy2 / 4000 - fu2 / 4000 + 1)) * (1 / (1 + exp(1e10 * -(ep_p - 0.04 + fu2 / E2))));
        G = (fabs(sigy1 - sigy) / 1e-9);
    }
    else {
        double fu_1 = (1.0 + b) * fy2; // 1+b gives the slope after (0.04, fu)
        double H_1 = (fu_1 - fy2) / 0.02;
        sigy = fu2 + H_1 * (epp - 0.04 + fu_1 / E2);
        double ep_p = epp + 1e-9;
        double sigy13 = fu2 + H_1 * (ep_p - 0.04 + fu_1 / E2);
        G = (fabs(sigy13 - sigy) / 1e-9);
    }
}

///// added by Anand Kumar IITJ
void
DPMsteelEC::Nr(double sg_trial1, double E1, double a1, double b1, double c1, double e_p, double fp1, double fy1, double fu1, double ey1, double alp_p, double &stress1, double &E_tgt, double& eplas1, double& alp1)
{
    int sign = 0;
    if (sg_trial1 < 0) {
        sign = -1;
    }
    else {
        sign = 1;
    }
    double d_lambda = 0.0;
    double d_d_lambda = 0.0;
    double r1 = 0.0;
    double G1 = 0.0;
    double r2 = 0.0;
    double G2 = 0.0;
    for (int i = 0; i < 25; i++) {

        this->yield_funtion(E1, a1, b1, c1, (alp_p + d_lambda), fp1, fy1, fu1, ey1, r1, G1);
        double sg_trial = fabs(sg_trial1);
        double rs = sg_trial - E1 * d_lambda - r1;
        d_d_lambda = rs / (G1 + E1);
        d_lambda = d_lambda + d_d_lambda;
        if (fabs(d_d_lambda) <= 1e-6) {
            break;
        }
    }
    this->yield_funtion(E1, a1, b1, c1, (e_p + d_lambda), fp1, fy1, fu1, ey1, r2, G2);
    E_tgt = E1 * G2 / (E1 + G2);
    stress1 = sg_trial1 - E1 * d_lambda * sign;
    eplas1 = e_p + d_lambda * sign;
    alp1 = alp_p + d_lambda;
}

// AddingSensitivity:BEGIN ///////////////////////////////////
int
DPMsteelEC::setParameter(const char** argv, int argc, Parameter& param)
{

    if (strcmp(argv[0], "sigmaY") == 0 || strcmp(argv[0], "fy") == 0)
        return param.addObject(1, this);

    if (strcmp(argv[0], "E") == 0)
        return param.addObject(2, this);

    if (strcmp(argv[0], "b") == 0)
        return param.addObject(3, this);

    if (strcmp(argv[0], "a1") == 0)
        return param.addObject(4, this);

    if (strcmp(argv[0], "a2") == 0)
        return param.addObject(5, this);

    if (strcmp(argv[0], "a3") == 0)
        return param.addObject(6, this);

    if (strcmp(argv[0], "a4") == 0)
        return param.addObject(7, this);

    return -1;
}



int
DPMsteelEC::updateParameter(int parameterID, Information& info)
{
    switch (parameterID) {
    case -1:
        return -1;
    case 1:
        this->fy = info.theDouble;
        break;
    case 2:
        this->E0 = info.theDouble;
        break;
    case 3:
        this->b = info.theDouble;
        break;
    case 4:
        this->a1 = info.theDouble;
        break;
    case 5:
        this->a2 = info.theDouble;
        break;
    case 6:
        this->a3 = info.theDouble;
        break;
    case 7:
        this->a4 = info.theDouble;
        break;
    default:
        return -1;
    }

    Ttangent = E0;          // Initial stiffness

    return 0;
}




int
DPMsteelEC::activateParameter(int passedParameterID)
{
    parameterID = passedParameterID;

    return 0;
}



double
DPMsteelEC::getStressSensitivity(int gradIndex, bool conditional)
{
    // Initialize return value
    double gradient = 0.0;


    // Pick up sensitivity history variables
    double CstrainSensitivity = 0.0;
    double CstressSensitivity = 0.0;
    if (SHVs != 0) {
        CstrainSensitivity = (*SHVs)(0, gradIndex);
        CstressSensitivity = (*SHVs)(1, gradIndex);
    }


    // Assign values to parameter derivatives (depending on what's random)
    double fySensitivity = 0.0;
    double E0Sensitivity = 0.0;
    double bSensitivity = 0.0;
    if (parameterID == 1) {
        fySensitivity = 1.0;
    }
    else if (parameterID == 2) {
        E0Sensitivity = 1.0;
    }
    else if (parameterID == 3) {
        bSensitivity = 1.0;
    }


    // Compute min and max stress
    double Tstress;
    double dStrain = Tstrain - Cstrain;
    double sigmaElastic = Cstress + E0 * dStrain;
    double fyOneMinusB = fy * (1.0 - b);
    double Esh = b * E0;
    double c1 = Esh * Tstrain;
    double c2 = TshiftN * fyOneMinusB;
    double c3 = TshiftP * fyOneMinusB;
    double sigmaMax = c1 + c3;
    double sigmaMin = c1 - c2;


    // Evaluate stress sensitivity 
    if ((sigmaMax < sigmaElastic) && (fabs(sigmaMax - sigmaElastic) > 1e-5)) {
        Tstress = sigmaMax;
        gradient = E0Sensitivity * b * Tstrain
            + E0 * bSensitivity * Tstrain
            + TshiftP * (fySensitivity * (1 - b) - fy * bSensitivity);
    }
    else {
        Tstress = sigmaElastic;
        gradient = CstressSensitivity
            + E0Sensitivity * (Tstrain - Cstrain)
            - E0 * CstrainSensitivity;
    }
    if (sigmaMin > Tstress) {
        gradient = E0Sensitivity * b * Tstrain
            + E0 * bSensitivity * Tstrain
            - TshiftN * (fySensitivity * (1 - b) - fy * bSensitivity);
    }

    return gradient;
}

double
DPMsteelEC::getInitialTangentSensitivity(int gradIndex)
{
    // For now, assume that this is only called for initial stiffness 
    if (parameterID == 2) {
        return 1.0;
    }
    else {
        return 0.0;
    }
}


int
DPMsteelEC::commitSensitivity(double TstrainSensitivity, int gradIndex, int numGrads)
{
    if (SHVs == 0) {
        SHVs = new Matrix(2, numGrads);
    }


    // Initialize unconditaional stress sensitivity
    double gradient = 0.0;


    // Pick up sensitivity history variables
    double CstrainSensitivity = 0.0;
    double CstressSensitivity = 0.0;
    if (SHVs != 0) {
        CstrainSensitivity = (*SHVs)(0, gradIndex);
        CstressSensitivity = (*SHVs)(1, gradIndex);
    }


    // Assign values to parameter derivatives (depending on what's random)
    double fySensitivity = 0.0;
    double E0Sensitivity = 0.0;
    double bSensitivity = 0.0;
    if (parameterID == 1) {
        fySensitivity = 1.0;
    }
    else if (parameterID == 2) {
        E0Sensitivity = 1.0;
    }
    else if (parameterID == 3) {
        bSensitivity = 1.0;
    }


    // Compute min and max stress
    double Tstress;
    double dStrain = Tstrain - Cstrain;
    double sigmaElastic = Cstress + E0 * dStrain;
    double fyOneMinusB = fy * (1.0 - b);
    double Esh = b * E0;
    double c1 = Esh * Tstrain;
    double c2 = TshiftN * fyOneMinusB;
    double c3 = TshiftP * fyOneMinusB;
    double sigmaMax = c1 + c3;
    double sigmaMin = c1 - c2;


    // Evaluate stress sensitivity ('gradient')
    if ((sigmaMax < sigmaElastic) && (fabs(sigmaMax - sigmaElastic) > 1e-5)) {
        Tstress = sigmaMax;
        gradient = E0Sensitivity * b * Tstrain
            + E0 * bSensitivity * Tstrain
            + E0 * b * TstrainSensitivity
            + TshiftP * (fySensitivity * (1 - b) - fy * bSensitivity);
    }
    else {
        Tstress = sigmaElastic;
        gradient = CstressSensitivity
            + E0Sensitivity * (Tstrain - Cstrain)
            + E0 * (TstrainSensitivity - CstrainSensitivity);
    }
    if (sigmaMin > Tstress) {
        gradient = E0Sensitivity * b * Tstrain
            + E0 * bSensitivity * Tstrain
            + E0 * b * TstrainSensitivity
            - TshiftN * (fySensitivity * (1 - b) - fy * bSensitivity);
    }


    // Commit history variables
    (*SHVs)(0, gradIndex) = TstrainSensitivity;
    (*SHVs)(1, gradIndex) = gradient;

    return 0;
}

// AddingSensitivity:END /////////////////////////////////////////////


//this function is no use, just for the definiation of pure virtual function.
int DPMsteelEC::setTrialStrain (double strain, double strainRate)
{
  opserr << "DPMsteelEC::setTrialStrain (double strain, double strainRate) - should never be called\n";
  return 0;
}


int 
DPMsteelEC::getVariable(const char *variable, Information &info)
{
  if (strcmp(variable,"ThermalElongation") == 0) {
    info.theDouble = ThermalElongation;    
    return 0;
  }else if (strcmp(variable,"ElongTangent") == 0) {
     Vector *theVector = info.theVector;
     if (theVector != 0) {
      double tempT, ET, Elong, TempTmax;
      tempT = (*theVector)(0);
	  ET = (*theVector)(1);
	  Elong = (*theVector)(2);
      TempTmax = (*theVector)(3);
      this->getElongTangent(tempT, ET, Elong, TempTmax);
	  (*theVector)(0) = tempT;
      (*theVector)(1) = ET;
      (*theVector)(2) = Elong;
	  (*theVector)(3) = TempTmax;
    }
    return 0;
	//The following code block is added only for recording temperature and Thermalelongation
  }else if (strcmp(variable,"TempAndElong") == 0) {
    Vector *theVector = info.theVector;
	if (theVector!= 0) {
		(*theVector)(0) = Ttemperature;
        (*theVector)(1) = ThermalElongation;
	}else{
		opserr<<"null Vector in DPMsteelEC"<<endln;
	}
	return 0;
 }
 //end of adding "TempAndElong"
   return -1;
}


