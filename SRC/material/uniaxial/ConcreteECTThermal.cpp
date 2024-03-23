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
                                                                       
// Added by Liming Jiang (UoE)
// Created: 06/13
//
// Description: This file contains the class definition for 
// ConcreteECTThermal. ConcreteECTThermal is modified from Concrete02Thermal
// ConcreteECTthermal is dedicated to provide a concrete material which 
// strictly satisfy Eurocode regarding the temperature dependent properties.

// Concrete02 is written by FMK in the year of 2006 and based on Concr2.f



#include <stdlib.h>
#include <ConcreteECTThermal.h>
#include <OPS_Globals.h>
#include <float.h>
#include <Channel.h>
#include <Information.h>

#include <elementAPI.h>
#include <OPS_Globals.h>

//Added by AK 2023
//double Temp_i = 0.0;

//int it = 0;


void *
OPS_ConcreteECTThermal(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int    iData[1];
  double dData[7];
  int numData = 1;

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial ConcreteECTThermal tag" << endln;
    return 0;
  }

  numData = OPS_GetNumRemainingInputArgs();

  if (numData != 7) {
    opserr << "Invalid #args, want: uniaxialMaterial ConcreteECTThermal " << iData[0] << "fpc? epsc0? fpcu? epscu? rat? ft? Ets?\n";
    return 0;
  }

  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid #args, want: uniaxialMaterial ConcreteECTThermal " << iData[0] << "fpc? epsc0? fpcu? epscu? rat? ft? Ets?\n";
    return 0;
  }


  // Parsing was successful, allocate the material
  theMaterial = new ConcreteECTThermal(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], dData[6]);
  
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type ConcreteECTThermal Material\n";
    return 0;
  }

  return theMaterial;
}


ConcreteECTThermal::ConcreteECTThermal(int tag, double _fc, double _epsc0, double _fcu,
				     double _epscu, double _rat, double _ft, double _Ets):
  UniaxialMaterial(tag, MAT_TAG_ConcreteECTThermal),
  //fc(_fc), epsc0(_epsc0), fcu(_fcu), epscu(_epscu), rat(_rat), ft(_ft), Ets(_Ets)
  fcT(_fc), epsc0T(_epsc0), fcuT(_fcu), epscuT(_epscu), rat(_rat), ftT(_ft), EtsT(_Ets) //JZ
{

//JZ 07/10 /////////////////////////////////////////////////////////////start
  fc = fcT;
  epsc0 = epsc0T;
  fcu = fcuT;
  epscu = epscuT;
  ft = ftT;
  Ets = EtsT;
//JZ 07/10 /////////////////////////////////////////////////////////////end 

  ecminP = 0.0;
  deptP = 0.0;

  eP = 2.0*fc/epsc0;  
  //eP = 1.5*fc/epsc0; //for the euro code, the 2.0 should be changed into 1.5
  epsP = 0.0;
  sigP = 0.0;
  eps = 0.0;
  sig = 0.0;
  e = 2.0*fc/epsc0;
  //e = 1.5*fc/epsc0;//for the euro code, the 2.0 should be changed into 1.5

  //if epsc0 is not 0.0025, then epsc0 = strainRatio*0.0025
  strainRatio = epsc0/0.0025;
  ThermalElongation = 0; //initialize 

  cooling=0; //PK add
  TempP = 0.0; //Pk add previous temp
  TmaxP = 0.0;
  Tmax = 0.0;
  epsLitsp = 0.0; // AK add transient strain
  Eps_lits = 0.0;
}

ConcreteECTThermal::ConcreteECTThermal(void):
  UniaxialMaterial(0, MAT_TAG_ConcreteECTThermal)
{
 
}

ConcreteECTThermal::~ConcreteECTThermal(void)
{
  // Does nothing
}


	

UniaxialMaterial*
ConcreteECTThermal::getCopy(void)
{
  ConcreteECTThermal *theCopy = new ConcreteECTThermal(this->getTag(), fc, epsc0, fcu, epscu, rat, ft, Ets);
  
  return theCopy;
}

double
ConcreteECTThermal::getInitialTangent(void)
{
  return 2.0*fc/epsc0;
}

int
ConcreteECTThermal::setTrialStrain(double trialStrain, double FiberTemperature, double strainRate)
{
   double 	ec0 = fc * 2. / epsc0;//?
   //double 	ec0 = fc * 1.5 / epsc0; //JZ. 27/07/10 ??
  
  // retrieve concrete hitory variables
  
  ecmin = ecminP;
  dept = deptP;
  
  // calculate current strain
  
  eps = trialStrain;
  double deps = eps - epsP;
  
  // if the current strain is less than the smallest previous strain 
  // call the monotonic envelope in compression and reset minimum strain 
  
  if (eps < ecmin) {
    this->Compr_Envlp(eps, sig, e);
    ecmin = eps;
  } else {;
    
    // else, if the current strain is between the minimum strain and ept 
    // (which corresponds to zero stress) the material is in the unloading- 
    // reloading branch and the stress remains between sigmin and sigmax 
    
    // calculate strain-stress coordinates of point R that determines 
    // the reloading slope according to Fig.2.11 in EERC Report 
    // (corresponding equations are 2.31 and 2.32 
    // the strain of point R is epsR and the stress is sigmR 
    
    double epsr = (fcu - rat * ec0 * epscu) / (ec0 * (1.0 - rat));
    double sigmr = ec0 * epsr;
    
    // calculate the previous minimum stress sigmm from the minimum 
    // previous strain ecmin and the monotonic envelope in compression 
    
    double sigmm;
    double dumy;
    this->Compr_Envlp(ecmin, sigmm, dumy);
    
    // calculate current reloading slope Er (Eq. 2.35 in EERC Report) 
    // calculate the intersection of the current reloading slope Er 
    // with the zero stress axis (variable ept) (Eq. 2.36 in EERC Report) 
    
    double er = (sigmm - sigmr) / (ecmin - epsr);
    double ept = ecmin - sigmm / er;
    
    if (eps <= ept) {
      double sigmin = sigmm + er * (eps - ecmin);
      double sigmax = er * .5f * (eps - ept);
      sig = sigP + ec0 * deps;
      e = ec0;
      if (sig <= sigmin) { 
	sig = sigmin;
	e = er;
      }      if (sig >= sigmax) {
	sig = sigmax;
	e = 0.5 * er;
      }
    } else {
      
      // else, if the current strain is between ept and epn 
      // (which corresponds to maximum remaining tensile strength) 
      // the response corresponds to the reloading branch in tension 
      // Since it is not saved, calculate the maximum remaining tensile 
      // strength sicn (Eq. 2.43 in EERC Report) 
      
      // calculate first the strain at the peak of the tensile stress-strain 
      // relation epn (Eq. 2.42 in EERC Report) 
      
      double epn = ept + dept;
      double sicn;
      if (eps <= epn) {
	this->Tens_Envlp(dept, sicn, e);
	if (dept != 0.0) {
	  e = sicn / dept;
	} else {
	  e = ec0;
	}
	sig = e * (eps - ept);
      } else {
	
	// else, if the current strain is larger than epn the response 
	// corresponds to the tensile envelope curve shifted by ept 
	
	double epstmp = eps - ept;
	this->Tens_Envlp(epstmp, sig, e);
	dept = eps - ept;
      }
    }
  }
//opserr<<"trialStrain: "<<eps << "  Stress: "<<sig<< "Modulus: "<<e<<endln;
  return 0;
}



double 
ConcreteECTThermal::getStrain(void)
{
  return eps;
}

double 
ConcreteECTThermal::getStress(void)
{
  return sig;
}

double 
ConcreteECTThermal::getTangent(void)
{
  return e;
}

double 
ConcreteECTThermal::getThermalElongation(void) //***JZ
{
  return ThermalElongation;
}

double 
ConcreteECTThermal::getElongTangent(double TempT, double& ET, double& Elong, double TempTmax) //PK add to include max temp
{
  //material properties with temperature
    if (TempT < 1) {
        Temp = 0.0;
    }
    else {
        Temp = TempT;
    }
  //Temp = TempT;  //make up the 20 degree which is minus in the class of thermalfield
  Tempmax = TempTmax; //PK add max temp for cooling
  // The datas are from EN 1992 part 1-2-1 
  // Tensile strength at elevated temperature
  
  bool Lits = true;  // AK add for trainsient strain
  double PhiT = 0.0;  // AK add for trainsient strain
  double PhiTP = 0.0; // AK add for trainsient strain
  double DelT = 0.0; // AK add for trainsient strain
     
  if (Temp <= 80) {
	  ft = ftT;
  }
  else if (Temp <= 580) {
	  ft = (1.0 - 1.0*(Temp -80)/500)*ftT;
	  Ets = (1.0 - 1.0*(Temp -80)/500)*fcT * 1.5 / epsc0T;
	  //Ets = (1.0 - 1.0*(Temp -80)/500)*EtsT;
  }
  else {
	  ft = 1.0e-10;
	  Ets = 1.0e-10;
	  //ft = 0;
	  //Ets = 0;
  }

  // compression strength, at elevated temperature
  //   strain at compression strength, at elevated temperature
  //   ultimate (crushing) strain, at elevated temperature
  if (Temp <= 0) {
	  fc = fcT;
	  epsc0 = -0.0025;
	  fcu = fcuT;
	  epscu = -0.02;
	  //Ets = EtsT;  jz what is there the statement?
  }
  else if (Temp <= 80) {
	  fc = fcT;
	  epsc0 = -(0.0025 + (0.003-0.0025)*(Temp - 0)/(80 - 0));
	  fcu = fcuT;
	  epscu = -(0.0200 + (0.0225-0.0200)*(Temp - 0)/(80 - 0));
  }
  else if (Temp <= 180) {
      fc = fcT*(1 - (Temp - 80)*0.05/100);
	  epsc0 = -(0.0030 + (0.003833-0.0030)*(Temp - 80)/100);
      fcu = fcuT*(1 - (Temp - 80)*0.05/100);
	  epscu = -(0.0225 + (0.0225-0.0200)*(Temp - 80)/100);
  }
  else if (Temp <= 280) {
      fc = fcT*(0.95 - (Temp - 180)*0.1/100);
	  epsc0 = -(0.003833 + (0.0050-0.003833)*(Temp - 180)/100);
	  fcu = fcuT*(0.95 - (Temp - 180)*0.1/100);
	  epscu = -(0.0250 + 0.0025*(Temp - 180)/100);
  }
  else if (Temp <= 380) {
      fc = fcT*(0.85 - (Temp - 280)*0.1/100);
	  epsc0 = -(0.0050 + (0.006333-0.0050)*(Temp - 280)/100);
	  fcu = fcuT*(0.85 - (Temp - 280)*0.1/100);
	  epscu = -(0.0275 + 0.0025*(Temp - 280)/100);
  }
  else if (Temp <= 480) {
      fc = fcT*(0.75 - (Temp - 380)*0.15/100);
	  epsc0 = -(0.006333 + (0.008667-0.006333)*(Temp - 380)/100);
	  fcu = fcuT*(0.75 - (Temp - 380)*0.15/100);
	  epscu = -(0.03 + 0.0025*(Temp - 380)/100);
  }
  else if (Temp <= 580) {
      fc = fcT*(0.60 - (Temp - 480)*0.15/100);
	  epsc0 = -(0.008667 + (0.012667-0.008667)*(Temp - 480)/100);
	  fcu = fcuT*(0.60 - (Temp - 480)*0.15/100);
	  epscu = -(0.0325 + 0.0025*(Temp - 480)/100);
  }
  else if (Temp <= 680) {
      fc = fcT*(0.45 - (Temp - 580)*0.15/100);
	  epsc0 = -(0.012667 + (0.013333 - 0.012667) * (Temp - 580) / 100);
	  fcu = fcuT*(0.45 - (Temp - 580)*0.15/100);
	  epscu = -(0.035 + 0.0025*(Temp - 580)/100);
  }
  else if (Temp <= 780) {
      fc = fcT*(0.30 - (Temp - 680)*0.15/100);
	  epsc0 = -(0.013333 + (0.014 - 0.013333) * (Temp - 680) / 100);
	  fcu = fcuT*(0.30 - (Temp - 680)*0.15/100);
	  epscu = -(0.0375 + 0.0025*(Temp - 680)/100);
  }
  else if (Temp <= 880) {
      fc = fcT*(0.15 - (Temp - 780)*0.07/100);
	  epsc0 = -(0.014 + (0.015 - 0.014) * (Temp - 680) / 100);
	  fcu = fcuT*(0.15 - (Temp - 780)*0.07/100);
	  epscu = -(0.04 + 0.0025*(Temp - 780)/100);
  }
  else if (Temp <= 980) {
      fc = fcT*(0.08 - (Temp - 880)*0.04/100);
	  epsc0 = -0.0150;
	  fcu = fcuT*(0.08 - (Temp - 880)*0.04/100);
	  epscu = -(0.0425 + 0.0025*(Temp - 880)/100);
  }
  else if (Temp <= 1080) {
      fc = fcT*(0.04 - (Temp - 980)*0.03/100);
	  epsc0 = -0.0150;
	  fcu = fcuT*(0.04 - (Temp - 980)*0.03/100);
	  epscu = -(0.045 + 0.0025*(Temp - 980)/100);
  }
  else  {
      opserr << "the temperature is invalid\n"; 
  }
//jz assign a miner to the valuables

	 // epsc0 = epsc0T*strainRatio;
	 // epscu = epscuT*strainRatio;

  // caculation of thermal elongation
	  if (Temp <= 1) {
		  ThermalElongation = (Temp - 0) * 9.213e-6;
	  }
  else if (Temp <= 680) {
      ThermalElongation = -1.8e-4 + 9e-6 *(Temp+20) + 2.3e-11 *(Temp+20)*(Temp+20)*(Temp+20);
  }
  else if (Temp <= 1180) {
      ThermalElongation = 14e-3;
  }
  else {
	  opserr << "the temperature is invalid\n";
  }

  if (Lits) {      
      if (sigP < 0) {
          //opserr << "under the LITs\n";
          if (Temp <= 0) {
              PhiT = 0.0 ;
          }
          else if (Temp <= 80) {
              PhiT = (0.0 + (0.001 - 0.0) * (Temp - 0) / (80 - 0));
          }
          else if (Temp <= 180) {
              PhiT = (0.001 + (0.0018 - 0.001) * (Temp - 80) / 100) ;
          }
          else if (Temp <= 280) {
              PhiT = (0.0018 + (0.0024 - 0.0018) * (Temp - 180) / 100);
          }
          else if (Temp <= 380) {
              PhiT = (0.0024 + (0.0049 - 0.0024) * (Temp - 280) / 100);
          }
          else if (Temp <= 480) {
              PhiT = (0.0049 + (0.0106 - 0.0049) * (Temp - 380) / 100);
          }
          else if (Temp <= 580) {
              PhiT = (0.0106 + (0.0274 - 0.0106) * (Temp - 480) / 100);
          }
          else if (Temp <= 680) {
              PhiT = (0.0274 + (0.0389 - 0.0274) * (Temp - 580) / 100);
          }
          else if (Temp <= 780) {
              PhiT = (0.0389 + (0.0733 - 0.0389) * (Temp - 680) / 100);
          }
          else if (Temp <= 880) {
              PhiT = 0.0733 ;
          }
          else if (Temp <= 980) {
              PhiT = 0.0733 ;
          }
          else if (Temp <= 1080) {
              PhiT = 0.0733 ;
          }

          if (TempP <= 0) {
              PhiTP = 0.0;
          }
          else if (TempP <= 80) {
              PhiTP = (0.0 + (0.001 - 0.0) * (TempP - 0) / (80 - 0));
          }
          else if (TempP <= 180) {
              PhiTP = (0.001 + (0.0018 - 0.001) * (TempP - 80) / 100);
          }
          else if (TempP <= 280) {
              PhiTP = (0.0018 + (0.0024 - 0.0018) * (TempP - 180) / 100);
          }
          else if (TempP <= 380) {
              PhiTP = (0.0024 + (0.0049 - 0.0024) * (TempP - 280) / 100);
          }
          else if (TempP <= 480) {
              PhiTP = (0.0049 + (0.0106 - 0.0049) * (TempP - 380) / 100);
          }
          else if (TempP <= 580) {
              PhiTP = (0.0106 + (0.0274 - 0.0106) * (TempP - 480) / 100);
          }
          else if (TempP <= 680) {
              PhiTP = (0.0274 + (0.0389 - 0.0274) * (TempP - 580) / 100);
          }
          else if (TempP <= 780) {
              PhiTP = (0.0389 + (0.0733 - 0.0389) * (TempP - 680) / 100);
          }
          else if (TempP <= 880) {
              PhiTP = 0.0733;
          }
          else if (TempP <= 980) {
              PhiTP = 0.0733;
          }
          else if (TempP <= 1080) {
              PhiTP = 0.0733;
          }
                     
          Eps_lits = PhiTP * (sig / fc) + (PhiT- PhiTP) * (sigP / fc);
          
          if (Eps_lits > epsLitsp) {
              epsLitsp = Eps_lits;
          }         
      }
      if (epsLitsp > 0)
          ThermalElongation = ThermalElongation - epsLitsp;
  }
  ET = 2*fc/epsc0; 
  Elong = ThermalElongation;
  
  DelT = Temp - TempP;
  if (DelT > 1.0) {
      Tmax = Temp;
  }
 
  ///PK COOLING PART FOR DESCENDING BRANCH OF A FIRE//// 
  /// If temperature is less that previous commited temp then we have cooling taking place
  if ((Temp - TempP) < -1.0) {
      Tempmax = TmaxP;
      //opserr << "cooling " << Temp << " " << TempP <<" " << TmaxP << endln;
      
      double kappa;
      double fcmax; //compr strength at max temp
      double fcumax; //ultimate compr strength at max temp
      double fcamb; //compr strength at cooled ambient temp
      double fcuamb; //ultimate compr strength at cooled ambient temp
      double epsc0max; //strain at compression strength for the max temp
      double epscumax; //ultimate strain at ultimate compression strength for the max temp
      if (TempP == Tempmax) {
          //opserr << "cooling,T,TP,Tmax " << Temp << " " << TempP << " " << Tempmax <<endln;
      }
      // PK Determine residual compressive strength of concrete heated to the max temp and then having cooled down to ambient
      // This will be the same for all the timesteps during the cooling phase
      // PK 1st step is to determine Kc,Tempmax according to table in 3.2.2 (EN1994-1-2:2005)

      if (Tempmax < 0) {
          opserr << "max temperature cannot be less than zero " << " " << Tempmax << endln;
      }
      else if (Tempmax <= 80) {
          kappa = 1;
          fcmax = fcT;
          fcumax = fcuT;
      }
      else if (Tempmax <= 180) {
          kappa = 1 - (Tempmax - 80) * 0.05 / 100;
          fcmax = fcT * (1 - (Tempmax - 80) * 0.05 / 100);
          fcumax = fcuT * (1 - (Tempmax - 80) * 0.05 / 100);
      }
      else if (Tempmax <= 280) {
          kappa = 0.95 - (Tempmax - 180) * 0.1 / 100;
          fcmax = fcT * (0.95 - (Tempmax - 180) * 0.1 / 100);
          fcumax = fcuT * (0.95 - (Tempmax - 180) * 0.1 / 100);
      }
      else if (Tempmax <= 380) {
          kappa = 0.85 - (Tempmax - 280) * 0.1 / 100;
          fcmax = fcT * (0.85 - (Tempmax - 280) * 0.1 / 100);
          fcumax = fcuT * (0.85 - (Tempmax - 280) * 0.1 / 100);
      }
      else if (Tempmax <= 480) {
          kappa = 0.75 - (Tempmax - 380) * 0.15 / 100;
          fcmax = fcT * (0.75 - (Tempmax - 380) * 0.15 / 100);
          fcumax = fcuT * (0.75 - (Tempmax - 380) * 0.15 / 100);
      }
      else if (Tempmax <= 580) {
          kappa = 0.60 - (Tempmax - 480) * 0.15 / 100;
          fcmax = fcT * (0.60 - (Tempmax - 480) * 0.15 / 100);
          fcumax = fcuT * (0.60 - (Tempmax - 480) * 0.15 / 100);
      }
      else if (Tempmax <= 680) {
          kappa = 0.45 - (Tempmax - 580) * 0.15 / 100;
          fcmax = fcT * (0.45 - (Tempmax - 580) * 0.15 / 100);
          fcumax = fcuT * (0.45 - (Tempmax - 580) * 0.15 / 100);
      }
      else if (Tempmax <= 780) {
          kappa = 0.30 - (Tempmax - 680) * 0.15 / 100;
          fcmax = fcT * (0.30 - (Tempmax - 680) * 0.15 / 100);
          fcumax = fcuT * (0.30 - (Tempmax - 680) * 0.15 / 100);
      }
      else if (Tempmax <= 880) {
          kappa = 0.15 - (Tempmax - 780) * 0.07 / 100;
          fcmax = fcT * (0.15 - (Tempmax - 780) * 0.07 / 100);
          fcumax = fcuT * (0.15 - (Tempmax - 780) * 0.07 / 100);
      }
      else if (Tempmax <= 980) {
          kappa = 0.08 - (Tempmax - 880) * 0.04 / 100;
          fcmax = fcT * (0.08 - (Tempmax - 880) * 0.04 / 100);
          fcumax = fcuT * (0.08 - (Tempmax - 880) * 0.04 / 100);
      }
      else if (Tempmax <= 1080) {
          kappa = 0.04 - (Tempmax - 980) * 0.03 / 100;
          fcmax = fcT * (0.04 - (Tempmax - 980) * 0.03 / 100);
          fcumax = fcuT * (0.04 - (Tempmax - 980) * 0.03 / 100);
      }
      else {
          opserr << "the temperature is invalid\n";
      }
      // PK 2nd step is to determine compressive strength at ambient after cooling as shown in ANNEX C (EN1994-1-2:2005)  
      if (Tempmax < 0) {
          opserr << "max temperature cannot be less than zero " << " " << Tempmax << endln;
      }
      else if (Tempmax <= 80) {
          fcamb = kappa * fcT;
          fcuamb = kappa * fcuT;
      }
      else if (Tempmax <= 280) {
          fcamb = (1 - (0.235 * (Tempmax - 80) / 200)) * fcT;
          fcuamb = (1 - (0.235 * (Tempmax - 80) / 200)) * fcuT;
      }
      else if (Tempmax <= 1080) {
          fcamb = 0.9 * kappa * fcT;
          fcuamb = 0.9 * kappa * fcuT;
      }
      else {
          opserr << "the temperature is invalid\n";
      }

      // Calculation of current compressive strength
      // linear interpolation between ambient and maximum compressive strength (after and before cooling)

      fc = fcmax - ((fcmax - fcamb) * (Tempmax - Temp) / Tempmax);
      fcu = fcumax - ((fcumax - fcuamb) * (Tempmax - Temp) / Tempmax);

      // Calculation of epsc0 for Tempmax and then keep it the same for all next time steps
      if (Tempmax < 0) {
          opserr << "max temperature cannot be less than zero " << " " << Tempmax << endln;
      }
      else if (Tempmax <= 80) {
          epsc0max = -(0.0025 + (0.004 - 0.0025) * (Tempmax - 0) / (80 - 0));
          epscumax = -(0.0200 + (0.0225 - 0.0200) * (Tempmax - 0) / (80 - 0));
      }
      else if (Tempmax <= 180) {
          epsc0max = -(0.0040 + (0.0055 - 0.0040) * (Tempmax - 80) / 100);
          epscumax = -(0.0225 + (0.0225 - 0.0200) * (Tempmax - 80) / 100);
      }
      else if (Tempmax <= 280) {
          epsc0max = -(0.0055 + (0.0070 - 0.0055) * (Tempmax - 180) / 100);
          epscumax = -(0.0250 + 0.0025 * (Tempmax - 180) / 100);
      }
      else if (Tempmax <= 380) {
          epsc0max = -(0.0070 + (0.0100 - 0.0070) * (Tempmax - 280) / 100);
          epscumax = -(0.0275 + 0.0025 * (Tempmax - 280) / 100);
      }
      else if (Tempmax <= 480) {
          epsc0max = -(0.0100 + (0.0150 - 0.0100) * (Tempmax - 380) / 100);
          epscumax = -(0.03 + 0.0025 * (Tempmax - 380) / 100);
      }
      else if (Tempmax <= 580) {
          epsc0max = -(0.0150 + (0.0250 - 0.0150) * (Tempmax - 480) / 100);
          epscumax = -(0.0325 + 0.0025 * (Tempmax - 480) / 100);
      }
      else if (Tempmax <= 680) {
          epsc0max = -0.0250;
          epscumax = -(0.035 + 0.0025 * (Tempmax - 580) / 100);
      }
      else if (Tempmax <= 780) {
          epsc0max = -0.0250;
          epscumax = -(0.0375 + 0.0025 * (Tempmax - 680) / 100);
      }
      else if (Tempmax <= 880) {
          epsc0max = -0.0250;
          epscumax = -(0.04 + 0.0025 * (Tempmax - 780) / 100);
      }
      else if (Tempmax <= 980) {
          epsc0max = -0.0250;
          epscumax = -(0.0425 + 0.0025 * (Tempmax - 880) / 100);
      }
      else if (Tempmax <= 1080) {
          epsc0max = -0.0250;
          epscumax = -(0.045 + 0.0025 * (Tempmax - 980) / 100);
      }
      else {
          opserr << "the temperature is invalid\n";
      }

      //make eps0 = eps0max

      epsc0 = epsc0max;

      // Calculating epscu
      epscu = epsc0 + ((epscumax - epsc0max) * fc / fcmax);

      ft = 0;

      double eres = 0.0;
      ThermalElong = 0.0;
      // caculation of thermal elongation
      if (Tempmax <= 1) {
          ThermalElong = (Tempmax - 0) * 9.213e-6;
      }
      else if (Tempmax <= 680) {
          ThermalElong = -1.8e-4 + 9e-6 * (Tempmax + 20) + 2.3e-11 * (Tempmax + 20) * (Tempmax + 20) * (Tempmax + 20);
      }
      else if (Tempmax <= 1180) {
          ThermalElong = 14e-3;
      }
      else {
          opserr << "the temperature is invalid\n";
      }

      // Make thermal elongation zero during the cooling phase
      
      if (Tempmax < 0) {
          opserr << "max temperature cannot be less than zero " << " " << Tempmax << endln;
      }
      else if (Tempmax <= 280) {
          //eres = (-0.00058 / 280) * Tempmax;
          eres = 0;       
      }
      else if (Tempmax <= 380) {
          //eres = -0.00058 + ((- 0.00029 + 0.00058) / 100)* (Tempmax - 280);
          eres = 0;        
      }
      else if (Tempmax <= 580) {
         // eres = -0.00029 + ((0.00171 + 0.00029) / 200) * (Tempmax - 380);
          eres = 0 + ((0.00171 + 0) / 200) * (Tempmax - 380);          
      }
      else if (Tempmax <= 780) {
          eres = 0.00171 + ((0.00329 - 0.00171) / 200) * (Tempmax - 580);          
      }
      else if (Tempmax <= 880) {
          eres = 0.00329 + ((0.005 - 0.00329) / 100) * (Tempmax - 780);
      }
      else {
          eres = 0.005;
         
      }
      
     // Elong = ThermalElong + ((eres - ThermalElong) / Tempmax) * (Tempmax - Temp) - epsLitsp;
      //if (Elong < 0) {
          //Elong = 0.0;
       //}

  }

  return 0;
}

int 
ConcreteECTThermal::commitState(void)
{
  ecminP = ecmin;
  deptP = dept;
  
  eP = e;
  sigP = sig;
  epsP = eps;
  
  TempP = Temp; //PK add set the previous temperature
  //epsLitsp = Eps_lits;
  TmaxP = Tmax;
  return 0;
}

int 
ConcreteECTThermal::revertToLastCommit(void)
{
  ecmin = ecminP;
  dept = deptP;
  
  e = eP;
  sig = sigP;
  eps = epsP;

  //Temp = TempP; //PK add set the previous temperature

  // NA ELENXW MIPWS EDW XANETAI TO TEMP LOGW MIN CONVERGENCE

  return 0;
}

int 
ConcreteECTThermal::revertToStart(void)
{
  ecminP = 0.0;
  deptP = 0.0;

  eP = 2.0*fc/epsc0;
  epsP = 0.0;
  sigP = 0.0;
  eps = 0.0;
  sig = 0.0;
  e = 2.0*fc/epsc0;

  return 0;
}

int 
ConcreteECTThermal::sendSelf(int commitTag, Channel &theChannel)
{
  static Vector data(13);
  data(0) =fc;    
  data(1) =epsc0; 
  data(2) =fcu;   
  data(3) =epscu; 
  data(4) =rat;   
  data(5) =ft;    
  data(6) =Ets;   
  data(7) =ecminP;
  data(8) =deptP; 
  data(9) =epsP;  
  data(10) =sigP; 
  data(11) =eP;   
  data(12) = this->getTag();

  if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "ConcreteECTThermal::sendSelf() - failed to sendSelf\n";
    return -1;
  }
  return 0;
}

int 
ConcreteECTThermal::recvSelf(int commitTag, Channel &theChannel, 
	     FEM_ObjectBroker &theBroker)
{

  static Vector data(13);

  if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "ConcreteECTThermal::recvSelf() - failed to recvSelf\n";
    return -1;
  }

  fc = data(0);
  epsc0 = data(1);
  fcu = data(2);
  epscu = data(3);
  rat = data(4);
  ft = data(5);
  Ets = data(6);
  ecminP = data(7);
  deptP = data(8);
  epsP = data(9);
  sigP = data(10);
  eP = data(11);
  this->setTag(data(12));

  e = eP;
  sig = sigP;
  eps = epsP;
  
  return 0;
}

void 
ConcreteECTThermal::Print(OPS_Stream &s, int flag)
{
  s << "ConcreteECTThermal:(strain, stress, tangent) " << eps << " " << sig << " " << e << endln;
}

void
ConcreteECTThermal::Tens_Envlp (double epsc, double &sigc, double &Ect)
{
/*-----------------------------------------------------------------------
! monotonic envelope of concrete in tension (positive envelope)
!
!   ft    = concrete tensile strength
!   Ec0   = initial tangent modulus of concrete 
!   Ets   = tension softening modulus
!   eps   = strain
!
!   returned variables
!    sigc  = stress corresponding to eps
!    Ect  = tangent concrete modulus
!-----------------------------------------------------------------------*/
  
  double Ec0  = 2.0*fc/epsc0;
  //double Ec0  = 1.5*fc/epsc0;

  double eps0 = ft/Ec0;
  double epsu = ft*(1.0/Ets+1.0/Ec0);
  if (epsc<=eps0) {
    sigc = epsc*Ec0;
    Ect  = Ec0;
  } else {
    if (epsc<=epsu) {
      Ect  = -0.1*Ets;
      sigc = ft-Ets*(epsc-eps0);
    } else {
      //      Ect  = 0.0
      Ect  = 1.0e-10;
      sigc = 1e-10;
    }
  }
  return;
}

  
void
ConcreteECTThermal::Compr_Envlp (double epsc, double &sigc, double &Ect) 
{
/*-----------------------------------------------------------------------
! monotonic envelope of concrete in compression (negative envelope)
!
!   fc    = concrete compressive strength
!   epsc0 = strain at concrete compressive strength
!   fcu   = stress at ultimate (crushing) strain 
!   epscu = ultimate (crushing) strain
!   Ec0   = initial concrete tangent modulus
!   epsc  = strain
!
!   returned variables
!   sigc  = current stress
!   Ect   = tangent concrete modulus
-----------------------------------------------------------------------*/

  double Ec0  = 2.0*fc/epsc0;
  //double Ec0  = 1.5*fc/epsc0;

  double ratLocal = epsc/epsc0;
  if (epsc>epsc0) {
    double ratSquare = ratLocal*ratLocal;
    double rSquare = (1 + ratSquare) * (1 + ratSquare);
    sigc = 2*ratLocal*fc/(1.0+ratSquare);
    Ect  = (2*fc/epsc0)*(1-ratSquare)/rSquare;
  } else {
    
    //   linear descending branch between epsc0 and epscu
    if (epsc>epscu) {
      sigc = (fcu-fc)*(epsc-epsc0)/(epscu-epsc0)+fc;
      Ect  = (fcu-fc)/(epscu-epsc0);
    } else {
	   
      // flat friction branch for strains larger than epscu
      
      sigc = fcu;
      Ect  = 1.0e-10;
      //       Ect  = 0.0
    }
  }

  return;
}

int
ConcreteECTThermal::getVariable(const char *varName, Information &theInfo)
{
  if (strcmp(varName,"ec") == 0) {
    theInfo.theDouble = epsc0;
    return 0;
  } else if (strcmp(varName,"ElongTangent") == 0) {
    Vector *theVector = theInfo.theVector;
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
  }
  else if(strcmp(varName,"ThermalElongation") == 0) {
	    theInfo.theDouble = ThermalElongation;
		return 0;
	}
  else if (strcmp(varName,"TempAndElong") == 0) {
    Vector *theVector = theInfo.theVector;
	if (theVector!= 0) {
		(*theVector)(0) = Temp;
        (*theVector)(1) = ThermalElongation;
	}else{
		opserr<<"null Vector in EC"<<endln;
	}

	return 0;
  }
  return -1;
}



//this function is no use, just for the definiation of pure virtual function.
int
ConcreteECTThermal::setTrialStrain(double strain, double strainRate)
{
  opserr << "ConcreteECTThermal::setTrialStrain(double strain, double strainRate) - should never be called\n";
  return -1;
}
