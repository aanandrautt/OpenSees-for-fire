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
                                                                       
// Added by Anand Kumar IITJ 
// Created: 02/2023


// Concrete02 is written by FMK in the year of 2006 and based on Concr2.f



#include <stdlib.h>
#include <DamagePlasticityConcreteECT.h>
#include <OPS_Globals.h>
#include <float.h>
#include <Channel.h>
#include <Information.h>
#include <MaterialResponse.h>   //Added by Anand Kumar [IITJ 2024] for damage recorder

#include <elementAPI.h>
#include <OPS_Globals.h>

void *
OPS_DamagePlasticityConcreteECT(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int    iData[1];
  double dData[7];
  int numData = 1;

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial DamagePlasticityConcreteECT tag" << endln;
    return 0;
  }

  numData = OPS_GetNumRemainingInputArgs();

  if (numData != 7) {
    opserr << "Invalid #args, want: uniaxialMaterial DamagePlasticityConcreteECT " << iData[0] << "fpc? epsc0? fpcu? epscu? rat? ft? Ets?\n";
    return 0;
  }

  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid #args, want: uniaxialMaterial DamagePlasticityConcreteECT " << iData[0] << "fpc? epsc0? fpcu? epscu? rat? ft? Ets?\n";
    return 0;
  }


  // Parsing was successful, allocate the material
  theMaterial = new DamagePlasticityConcreteECT(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], dData[6]);
  
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type DamagePlasticityConcreteECT Material\n";
    return 0;
  }

  return theMaterial;
}


DamagePlasticityConcreteECT::DamagePlasticityConcreteECT(int tag, double _fc, double _epsc0, double _fcu,
				     double _epscu, double _rat, double _ft, double _Ets):
  UniaxialMaterial(tag, MAT_TAG_DamagePlasticityConcreteECT),
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
  Tstrain = 0.0;

  cooling=0; //PK add
  TempP = 0.0; //Pk add previous temp
  TmaxP = 0.0;
  Tmax = 0.0;
  epsLitsp = 0.0; // AK add transient strain
  Eps_lits = 0.0;
  Eps_litsP = 0.0;

  e_res = 0.0; // AK added to keep track of residual strain.
  e_resP = 0.0;
  dc = 0.0;
  dcP = 0.0;
  dt = 0;
  dtP = 0.0;
  damageC = 0.0;
  damageT = 0.0;
  epti = 0.0;
  eptiP = 0.0;
}

DamagePlasticityConcreteECT::DamagePlasticityConcreteECT(void):
  UniaxialMaterial(0, MAT_TAG_DamagePlasticityConcreteECT)
{
 
}

DamagePlasticityConcreteECT::~DamagePlasticityConcreteECT(void)
{
  // Does nothing
}


	

UniaxialMaterial*
DamagePlasticityConcreteECT::getCopy(void)
{
  DamagePlasticityConcreteECT *theCopy = new DamagePlasticityConcreteECT(this->getTag(), fc, epsc0, fcu, epscu, rat, ft, Ets);
  
  return theCopy;
}

double
DamagePlasticityConcreteECT::getInitialTangent(void)
{
  return 2.0*fc/epsc0;
}

int
DamagePlasticityConcreteECT::setTrialStrain(double trialStrain, double FiberTemperature, double strainRate)
{
    if (fc == fcT && ft == ftT) {

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
        }
        else {
            ;

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
            }
            else {

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
                    }
                    else {
                        e = ec0;
                    }
                    sig = e * (eps - ept);
                }
                else {

                    // else, if the current strain is larger than epn the response 
                    // corresponds to the tensile envelope curve shifted by ept 

                    double epstmp = eps - ept;
                    this->Tens_Envlp(epstmp, sig, e);
                    dept = eps - ept;
                }
            }
        }
        //opserr<<"trialStrain: "<<eps << "  Stress: "<<sig<< "Modulus: "<<e<<endln;
    }

    else {       
        double 	ec0 = fc * 2. / epsc0;//?      

        ecmin = ecminP;
        dept = deptP;        

        // calculate current strain

        eps = trialStrain;
        double deps = eps - epsP;

        double epsr = (fcu - rat * ec0 * epscu) / (ec0 * (1.0 - rat));
        double sigmr = ec0 * epsr;
        
        double sigmm;
        double dumy;
        this->Compr_Envlp(ecmin, sigmm, dumy);

        double er = (sigmm - sigmr) / (ecmin - epsr);
        double ept = ecmin - sigmm / er;

        // if the current strain is less than the smallest previous strain 
        // call the monotonic envelope in compression and reset minimum strain 

        if (eps < ecmin) {           
            double sig1, e1;
            this->Compr_Envlp(eps, sig1, e1);          
                        
            double E_rr = (1 - dcP) * ec0;
            double e_elas = eps - e_resP;
            double sig_star = E_rr * e_elas;

            if (sig_star < sig1) {                
                sig = sig1;
                dc = (1 - er / ec0);
                e_res = eps - sigmm / er;
                ecmin = eps;
                e = e1;
            }
            else {
                sig = sig_star;
                dc = dcP;
                e_res = e_resP;  
                ecmin = eps;     
                e = E_rr;
            }
           damageC = dc;
        }
        else {   
            double sig3, e3;
            this->Compr_Envlp(eps, sig3, e3);

            double E_rr = (1 - dcP) * ec0;
            double e_elas = eps - e_resP;
            double sig_star = E_rr * e_elas;

            double sigmin = sigmm + er * (eps - ecmin);
            double sigmax = er * 0.5 * (eps - ept);
            double sig2 = sigP + ec0 * deps;
            double e2 = 0.0;

            if (sig2 <= sigmin) {
                sig2 = sigmin;
                e2 = ec0;
            }   if (sig2 >= sigmax) {
                sig2 = sigmax;
                e2 = 0.5 * ec0;
            }

            
            if (eps <= ept) {
                if (DelTP <= -0.1) {
                    if (sig_star > sig3) {                                             
                        sig = sig_star;
                        e_res = e_resP;
                        dc = dcP;
                        e = E_rr;
                    }
                    else {
                        sig = sig3;
                        ecmin = eps;
                        dc = (1 - er / ec0);
                        e_res = eps - sigmm / er;
                        e = e3;
                    }
                }
                else if (DelTP == 0) {
                    sig = sig2;
                    e_res = e_resP;
                    dc = dcP;
                    e = e2;
                }
                else {
                    if (sig_star < sig2) {
                        sig = sig2;
                        dc = dcP;
                        e_res = e_resP;
                        e = e2;
                    }
                    else {
                        sig = sig_star;
                        dc = dcP;
                        e_res = e_resP;
                        e = E_rr;
                    }
                }
                damageC = dc;
            }
            else {      

                double epn = ept + dept;
                double sicn;
           
                double eptu = ft * (1 / Ets + 1 / ec0);
                double eptd = fabs(eps - ept);
                
                epti = eptu / (1 + (1 - dtP) * (ec0 / Ets));
                

                if (eptd >= (ft / ec0)) {
                    dt = 1 - (Ets * fabs(eptu - eptd) / (ec0 * eptd));
                    if (dt > dtP) {
                        dtP = dt;
                    }
                   
                   damageT = dtP;
                }

                if (eps <= epn) {
                    this->Tens_Envlp(dept, sicn, e);
                    if (dept != 0.0) {
                        e = sicn / dept;                  
                    }
                    else {
                        e = ec0;                     
                    }
                    sig = e * (eps - ept);                  
                }
                else {
                    double epstmp = eps - ept;
                    this->Tens_Envlp(epstmp, sig, e);
                    dept = eps - ept;                 
                } 
            }
        }
    }
  return 0;
}



double 
DamagePlasticityConcreteECT::getStrain(void)
{
  return eps;
}

double 
DamagePlasticityConcreteECT::getStress(void)
{
  return sig;
}

double 
DamagePlasticityConcreteECT::getTangent(void)
{
  return e;
}

double
DamagePlasticityConcreteECT::getDamageC(void)
{
   return damageC;
}

double
DamagePlasticityConcreteECT::getDamageT(void)
{
    return damageT;
}

double 
DamagePlasticityConcreteECT::getThermalElongation(void) //***JZ
{
  return ThermalElongation;
}

double
DamagePlasticityConcreteECT::getTstrain(void) //***JZ
{
    return Tstrain;
}

double
DamagePlasticityConcreteECT::getLITS(void) //***JZ
{
    return epsLitsp;
}

double 
DamagePlasticityConcreteECT::getElongTangent(double TempT, double& ET, double& Elong, double TempTmax) //PK add to include max temp
{
  //material properties with temperature
    if (TempT < 0.1) {
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
  //double DelT = 0.0; // AK add for trainsient strain
       
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
	  fc = fcT * 0.99999;    // to make fc =! fcT to detect start of heating phase // AK
	  epsc0 = -(0.0025 + (0.003-0.0025)*(Temp - 0)/(80 - 0));
	  fcu = fcuT;
	  epscu = -(0.0200 + (0.0225-0.0200)*(Temp - 0)/(80 - 0));
  }
  else if (Temp <= 180) {
      fc = fcT * (1 - (Temp - 80) * 0.05 / 100);
      epsc0 = -(0.0030 + (0.003833 - 0.0030) * (Temp - 80) / 100);
      fcu = fcuT * (1 - (Temp - 80) * 0.05 / 100);
      epscu = -(0.0225 + (0.0225 - 0.0200) * (Temp - 80) / 100);
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
	  epsc0 = -(0.014 + (0.015 - 0.014) * (Temp - 780) / 100);
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
      // Added by Anand Kumar to handle temperature greater than 1100
      fc = 0.00001 * fcT;
      epsc0 = -0.0150;
      fcu = fcuT * 0.0001;
      epscu = -0.05;
  }
//jz assign a miner to the valuables

	 // epsc0 = epsc0T*strainRatio;
	 // epscu = epscuT*strainRatio;

  // caculation of thermal elongation
	  if (Temp <= 1) {
		  ThermalElongation = (Temp - 0) * 9.213e-6;
          Tstrain = (Temp - 0) * 9.213e-6;
	  }
  else if (Temp <= 680) {
      ThermalElongation = -1.8e-4 + 9e-6 *(Temp+20) + 2.3e-11 *(Temp+20)*(Temp+20)*(Temp+20);
      Tstrain = -1.8e-4 + 9e-6 * (Temp + 20) + 2.3e-11 * (Temp + 20) * (Temp + 20) * (Temp + 20);
  }
  else if (Temp <= 1380) {
      ThermalElongation = 14e-3;   //Temp <= 1180 is modified as 1380 to handle temperature greater than 1200
      Tstrain = 14e-3;
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
                     
          Eps_lits = PhiTP * (sig / fcT) + (PhiT- PhiTP) * (sigP / fcT);
          //Eps_lits = Eps_litsP + (PhiT - PhiTP) * (sigP / fcT);
          
          if (Eps_lits > epsLitsp) {
              epsLitsp = Eps_lits;
          }
      }
      if (epsLitsp > 0) {
          ThermalElongation = ThermalElongation - epsLitsp;
      }
  }
  ET = 2*fc/epsc0; 
  Elong = ThermalElongation;

  DelT = Temp - TempP;
  if (DelT > 0.091) {
      Tmax = Temp;
  }
 
  ///PK COOLING PART FOR DESCENDING BRANCH OF A FIRE//// 
  /// If temperature is less that previous commited temp then we have cooling taking place
  if ((Temp - TempP) < -0.091) {
      //cooling = 1.0;
      //opserr << "cooling " << cooling << endln;
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
          //opserr << "the temperature is invalid\n";
          // added by anand kumar to handle temperature greater than 1100
          kappa = 0.000001;
          fcmax = fcT * (0.000001);
          fcumax = fcuT * (0.000001);
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
          fcamb = (1 - (0.235 * (Tempmax - 80) / 200)) * fcT;  //0.235
          fcuamb = (1 - (0.235 * (Tempmax - 80) / 200)) * fcuT;
      }
      else if (Tempmax <= 1080) {
          fcamb = 0.9 * kappa * fcT;
          fcuamb = 0.9* kappa * fcuT;
      }
      else {
          //opserr << "the temperature is invalid\n";
          // added by anand kumar to handle temperature greater than 1100
          fcamb = 0.9 * kappa * fcT;
          fcuamb = 0.9 * kappa * fcuT;
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
          //epsc0max = -(0.0025 + (0.004 - 0.0025) * (Tempmax - 0) / (80 - 0));
          epsc0max = -(0.0025 + (0.003 - 0.0025) * (Tempmax - 0) / (80 - 0));
          epscumax = -(0.0200 + (0.0225 - 0.0200) * (Tempmax - 0) / (80 - 0));
      }
      else if (Tempmax <= 180) {
          //epsc0max = -(0.0040 + (0.0055 - 0.0040) * (Tempmax - 80) / 100);
          epsc0max = -(0.0030 + (0.003833 - 0.0030) * (Tempmax - 80) / 100);
          epscumax = -(0.0225 + (0.0225 - 0.0200) * (Tempmax - 80) / 100);
      }
      else if (Tempmax <= 280) {
          //epsc0max = -(0.0055 + (0.0070 - 0.0055) * (Tempmax - 180) / 100);
          epsc0max = -(0.003833 + (0.0050 - 0.003833) * (Tempmax - 180) / 100);
          epscumax = -(0.0250 + 0.0025 * (Tempmax - 180) / 100);
      }
      else if (Tempmax <= 380) {
          //epsc0max = -(0.0070 + (0.0100 - 0.0070) * (Tempmax - 280) / 100);
          epsc0max = -(0.0050 + (0.006333 - 0.0050) * (Tempmax - 280) / 100);
          epscumax = -(0.0275 + 0.0025 * (Tempmax - 280) / 100);
      }
      else if (Tempmax <= 480) {
          //epsc0max = -(0.0100 + (0.0150 - 0.0100) * (Tempmax - 380) / 100);
          epsc0max = -(0.006333 + (0.008667 - 0.006333) * (Tempmax - 380) / 100);
          epscumax = -(0.03 + 0.0025 * (Tempmax - 380) / 100);
      }
      else if (Tempmax <= 580) {
          //epsc0max = -(0.0150 + (0.0250 - 0.0150) * (Tempmax - 480) / 100);
          epsc0max = -(0.008667 + (0.012667 - 0.008667) * (Tempmax - 480) / 100);
          epscumax = -(0.0325 + 0.0025 * (Tempmax - 480) / 100);
      }
      else if (Tempmax <= 680) {
          //epsc0max = -0.0250;
          epsc0max = -(0.012667 + (0.013333 - 0.012667) * (Tempmax - 580) / 100);
          epscumax = -(0.035 + 0.0025 * (Tempmax - 580) / 100);
      }
      else if (Tempmax <= 780) {
          //epsc0max = -0.0250;
          epsc0max = -(0.013333 + (0.014 - 0.013333) * (Temp - 680) / 100);
          epscumax = -(0.0375 + 0.0025 * (Tempmax - 680) / 100);
      }
      else if (Tempmax <= 880) {
          //epsc0max = -0.0250;
          epsc0max = -(0.014 + (0.015 - 0.014) * (Temp - 780) / 100);
          epscumax = -(0.04 + 0.0025 * (Tempmax - 780) / 100);
      }
      else if (Tempmax <= 980) {
          //epsc0max = -0.0250;
          epsc0max = -0.0150;
          epscumax = -(0.0425 + 0.0025 * (Tempmax - 880) / 100);
      }
      else if (Tempmax <= 1080) {
          //epsc0max = -0.0250;
          epsc0max = -0.0150;
          epscumax = -(0.045 + 0.0025 * (Tempmax - 980) / 100);
      }
      else {
          //opserr << "the temperature is invalid\n";
          // added by anand kumar to handle temperature greater than 1100 
          //epsc0max = -0.0250;
          epsc0max = -0.0150;
          epscumax = -0.05;

      }

      //make eps0 = eps0max

      epsc0 = epsc0max;

      // Calculating epscu
      epscu = epsc0 + ((epscumax - epsc0max) * fc / fcmax);
      
      ft = 0;
  }
  return 0;
}

int 
DamagePlasticityConcreteECT::commitState(void)
{
  ecminP = ecmin;
  deptP = dept;
  
  eP = e;
  sigP = sig;
  epsP = eps;
  
  TempP = Temp; //PK add set the previous temperature
  //epsLitsp = Eps_lits;
  Eps_litsP = Eps_lits;

  e_resP = e_res; // Anand Kumar added to keep track of residual strain.
  dcP = dc;
  DelTP = DelT;
  //eptiP = epti;

  TmaxP = Tmax;
  return 0;
}

int 
DamagePlasticityConcreteECT::revertToLastCommit(void)
{
  ecmin = ecminP;
  dept = deptP;
  
  e = eP;
  sig = sigP;
  eps = epsP;
  Eps_lits = Eps_litsP;

  e_res = e_resP; // Anand Kumar added to keep track of residual strain.
  dc = dcP;
  //dt = dtP;
  DelT = DelTP;
  Temp = TempP; //PK add set the previous temperature

  // NA ELENXW MIPWS EDW XANETAI TO TEMP LOGW MIN CONVERGENCE

  return 0;
}

int 
DamagePlasticityConcreteECT::revertToStart(void)
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
DamagePlasticityConcreteECT::sendSelf(int commitTag, Channel &theChannel)
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
    opserr << "DamagePlasticityConcreteECT::sendSelf() - failed to sendSelf\n";
    return -1;
  }
  return 0;
}

int 
DamagePlasticityConcreteECT::recvSelf(int commitTag, Channel &theChannel, 
	     FEM_ObjectBroker &theBroker)
{

  static Vector data(13);

  if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "DamagePlasticityConcreteECT::recvSelf() - failed to recvSelf\n";
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
DamagePlasticityConcreteECT::Print(OPS_Stream &s, int flag)
{
  s << "DamagePlasticityConcreteECT:(strain, stress, tangent) " << eps << " " << sig << " " << e << endln;
}

void
DamagePlasticityConcreteECT::Tens_Envlp (double epsc, double &sigc, double &Ect)
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
DamagePlasticityConcreteECT::Compr_Envlp (double epsc, double &sigc, double &Ect) 
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

double
DamagePlasticityConcreteECT::newton_raphson_(double Err, double e_res, double initial_guess, double tolerance, int max_iterations) {
    double e = initial_guess;
    int iteration = 0;
    double el = 0.0;
    double tol = 0.0;

    while (iteration < max_iterations && tol == 0) {
        double r;
        double dr_de;
        this->Residual_(e, r, dr_de, Err, e_res);

        // Update del_e using the Newton-Raphson formula
        double del_e = r / dr_de;

        // Check for convergence based on tolerance
        if (fabs(r) < tolerance) {
            //std::cout << "Converged after " << iteration + 1 << " iterations." << std::endl;
            tol = 1;
            el = e - del_e;
        }
        e = e - del_e;
        iteration++;
    }
    return el;
}

void
DamagePlasticityConcreteECT::Residual_(double ee, double& r, double& dr_de, double Err, double e_res) {
    if (ee > epsc0) {
        double ratLocal = ee / epsc0;
        double ratSquare = ratLocal * ratLocal;
        double rSquare = (1 + ratSquare) * (1 + ratSquare);
        double sigc = 2 * ratLocal * fc / (1.0 + ratSquare);
        double Ect = (2 * fc / epsc0) * (1 - ratSquare) / rSquare;
        r = sigc - Err * (ee - e_res);
        dr_de = Ect - Err;
    }
    else {
        if (ee > epscu) {
            double sigc = (fcu - fc) * (ee - epsc0) / (epscu - epsc0) + fc;
            double Ect = (fcu - fc) / (epscu - epsc0);
            r = sigc - Err * (ee - e_res);
            dr_de = Ect - Err;
        }
        else {
            double sigc = fcu;
            double Ect = 1.0e-10;
            r = sigc - Err * (ee - e_res);
            dr_de = 0.0 - Err;
        }
    }
}

int
DamagePlasticityConcreteECT::getVariable(const char *varName, Information &theInfo)
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
DamagePlasticityConcreteECT::setTrialStrain(double strain, double strainRate)
{
  opserr << "DamagePlasticityConcreteECT::setTrialStrain(double strain, double strainRate) - should never be called\n";
  return -1;
}


/* Methods added by Anand Kumar IITJ to record Creep, Mechanical and Thermal strain: */

Response*
DamagePlasticityConcreteECT::setResponse(const char** argv, int argc,
    OPS_Stream& theOutput)
{
    Response* theResponse = 0;

    theOutput.tag("UniaxialMaterialOutput");
    theOutput.attr("matType", this->getClassType());
    theOutput.attr("matTag", this->getTag());

    // To record compression damage in concrete {Added by Anand Kumar IITJ}
    if (strcmp(argv[0], "damageC") == 0) {
        theOutput.tag("ResponseType", "dc");
        theResponse = new MaterialResponse(this, 1, this->getDamageC());
    }
    // To record tension damage in concrete {Added by Anand Kumar IITJ}
    else if (strcmp(argv[0], "damageT") == 0) {
        theOutput.tag("ResponseType", "dt");
        theResponse = new MaterialResponse(this, 2, this->getDamageT());
    }

    // Dummies: not required here
    else if (strcmp(argv[0], "strain") == 0) {
        theOutput.tag("ResponseType", "eps11");
        theResponse = new MaterialResponse(this, 3, this->getStrain());
    }

    // Dummies: not required here
    else if ((strcmp(argv[0], "stressStrain") == 0) ||
        (strcmp(argv[0], "stressANDstrain") == 0) ||
        (strcmp(argv[0], "stressAndStrain") == 0)) {
        theOutput.tag("ResponseType", "sig11");
        theOutput.tag("ResponseType", "eps11");
        theResponse = new MaterialResponse(this, 4, Vector(2));
    }
     // Dummies: not required here
    else if (strcmp(argv[0], "CreepStressStrainTangent") == 0) {
        theOutput.tag("ResponseType", "Tstrain");
        theOutput.tag("ResponseType", "LITS");
        theResponse = new MaterialResponse(this, 6, Vector(2));
    }
    // Dummies: not required here
    else if ((strcmp(argv[0], "stressStrainTangent") == 0) ||
        (strcmp(argv[0], "stressANDstrainANDtangent") == 0)) {
        theOutput.tag("ResponseType", "sig11");
        theOutput.tag("ResponseType", "eps11");
        theOutput.tag("ResponseType", "C11");
        theResponse = new MaterialResponse(this, 5, Vector(3));
    }

    // stress sensitivity for local sensitivity recorder purpose.  Quan 2009
    // limit:  no more than 10000 random variables/sensitivity parameters
    else if (strstr(argv[0], "stressSensitivity") != 0) {
        char* token = strtok((char*)argv[0], " ");
        if (token != NULL) token = strtok(NULL, " ");
        int gradient = atoi(token);
        theOutput.tag("ResponseType", "sigsens11");
        theResponse = new MaterialResponse(this, gradient + 10000, this->getStress());
    }
    // strain sensivitiy
    else if (strstr(argv[0], "strainSensitivity") != 0) {
        char* token = strtok((char*)argv[0], " ");
        if (token != NULL) token = strtok(NULL, " ");
        int gradient = atoi(token);
        theOutput.tag("ResponseType", "epssens11");
        theResponse = new MaterialResponse(this, gradient + 20000, this->getStrain());
    }


    theOutput.endTag();
    return theResponse;

}

int
DamagePlasticityConcreteECT::getResponse(int responseID, Information& matInfo)
{
    static Vector stressStrain(2);
    static Vector stressStrainTangent(3);
    static Vector CreepStressStrainTangent(2); //Added by AMK
    // each subclass must implement its own stuff   

    // added for sensitivity recorder. Quan 2009
    if ((responseID > 10000) && (responseID < 20000)) {
        matInfo.setDouble(this->getStressSensitivity(responseID - 10000, false));
        return 0;
    }
    else if (responseID > 20000) {
        matInfo.setDouble(this->getStrainSensitivity(responseID - 20000));
        return 0;
    }

    switch (responseID) {
    case 1:
        matInfo.setDouble(this->getDamageC());
        return 0;

    case 2:
        matInfo.setDouble(this->getDamageT());
        return 0;

    case 3:
        matInfo.setDouble(this->getStrain());
        return 0;

    case 4:
        stressStrain(0) = this->getStress();
        stressStrain(1) = this->getStrain();
        matInfo.setVector(stressStrain);
        return 0;

    case 5:
        stressStrainTangent(0) = this->getStress();
        stressStrainTangent(1) = this->getStrain();
        stressStrainTangent(2) = this->getTangent();
        matInfo.setVector(stressStrainTangent);
        return 0;

    case 6:
        //CreepStressStrainTangent(0) = this->getStress();
        //CreepStressStrainTangent(1) = this->getStrain();
        //CreepStressStrainTangent(2) = this->getTangent();
        CreepStressStrainTangent(0) = this->getTstrain();
        CreepStressStrainTangent(1) = this->getLITS();
        //CreepStressStrainTangent(5) = this->getShrink();
        matInfo.setVector(CreepStressStrainTangent);
        return 0;

    default:
        return -1;
    }
}