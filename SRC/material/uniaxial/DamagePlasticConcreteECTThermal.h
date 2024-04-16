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
// DamagePlasticConcreteECTThermal. DamagePlasticConcreteECTThermal is modified from Concrete02Thermal
// DamagePlasticConcreteECTThermal is dedicated to provide a concrete material which 
// strictly satisfy Eurocode regarding the temperature dependent properties.

// Concrete02 is written by FMK in the year of 2006 and based on Concr2.f




#ifndef DamagePlasticConcreteECTThermal_h
#define DamagePlasticConcreteECTThermal_h

#include <UniaxialMaterial.h>

class DamagePlasticConcreteECTThermal : public UniaxialMaterial
{
  public:
    DamagePlasticConcreteECTThermal(int tag, double _fc, double _epsc0, double _fcu,
	     double _epscu, double _rat, double _ft, double _Ets);

    DamagePlasticConcreteECTThermal(void);

    virtual ~DamagePlasticConcreteECTThermal();

    const char *getClassType(void) const {return "DamagePlasticConcreteECTThermal";};    
    double getInitialTangent(void);
    UniaxialMaterial *getCopy(void);


    int setTrialStrain(double strain, double rate);     //JZ this function is no use, just for the definiation of pure virtual function.
    int setTrialStrain(double strain, double FiberTemperature, double strainRate); //***JZ

    double newton_raphson_(double Err, double e_elas, double initial_guess, double tolerance, int max_iterations); // added by AK

    double getStrain(void);      
    double getStress(void);
    double getTangent(void);

    double getThermalElongation(void); //***JZ
    double getElongTangent(double, double&, double&, double);//***JZ //PK add to include max temp
    
    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);        
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    
    void Print(OPS_Stream &s, int flag =0);

    int getVariable(const char *variable, Information &);
    

    
 protected:
    
 private:
    void Tens_Envlp (double epsc, double &sigc, double &Ect);
    void Compr_Envlp (double epsc, double &sigc, double &Ect);
    //double Phi_T(double T); //AK add for transient strain

    void Residual_(double ee, double& sigc, double& Ect, double Err, double e_elas); // added by AK

    double Temp;  // concrete temp
    double steps;    //the amount of the steps.
    double strainRatio; //input strain over 0.0025(EU 1992)  
    double ThermalElongation; // eps(theata) = alpha * temperature
    double ThermalElong;
    double fcT;    
    double epsc0T; 
    double fcuT;   
    double epscuT;         
    double ftT;    
    double EtsT;  
    double cooling; //PK add
    double Tempmax;  // PK add max temp
    double TmaxP; //AK add for cooling
    double Tmax;  //AK add for cooling

    
    // matpar : Concrete FIXED PROPERTIES
    double fc;    // concrete compression strength           : mp(1)
    double epsc0; // strain at compression strength          : mp(2)
    double fcu;   // stress at ultimate (crushing) strain    : mp(3)
    double epscu; // ultimate (crushing) strain              : mp(4)       
    double rat;   // ratio between unloading slope at epscu and original slope : mp(5)
    double ft;    // concrete tensile strength               : mp(6)
    double Ets;   // tension stiffening slope                : mp(7)



    // hstvP : Concerete HISTORY VARIABLES last committed step
    double ecminP;  //  hstP(1)
    double deptP;   //  hstP(2)
    double epsP;  //  = strain at previous converged step
    double sigP;  //  = stress at previous converged step
    double eP;    //   stiffness modulus at last converged step;
    
    double e_resP; // AK added to keep track of residual strain.
    double dcP;
    
	double TempP; //PK add the previous temperature

    // hstv : Concerete HISTORY VARIABLES  current step
    double ecmin;  
    double dept;   
    double sig;   
    double e;     
    double eps;   
   
    
    double epsLitsp; // AK add for transient strain
    double Eps_lits; // AK add for transient strain
    double Eps_litsP;
    double e_res; // AK added to keep track of residual strain.
    double dc;
};


#endif

