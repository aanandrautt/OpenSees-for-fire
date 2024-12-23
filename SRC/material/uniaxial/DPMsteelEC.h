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
//Modified by:  Anand Kumar IITJ(j.zhang@ed.ac.uk)---------02,2023// 
//             
                                                                                                                                                
#ifndef DPMsteelEC_h
#define DPMsteelEC_h


#include <UniaxialMaterial.h>

// Default values for isotropic hardening parameters a1, a2, a3, and a4
#define STEEL_01_DEFAULT_A1        0.0
#define STEEL_01_DEFAULT_A2       55.0
#define STEEL_01_DEFAULT_A3        0.0
#define STEEL_01_DEFAULT_A4       55.0

class DPMsteelEC : public UniaxialMaterial
{
  public:
    DPMsteelEC(int tag, double fy, double E0, double b,
		   double a1 = STEEL_01_DEFAULT_A1, double a2 = STEEL_01_DEFAULT_A2,
		   double a3 = STEEL_01_DEFAULT_A3, double a4 = STEEL_01_DEFAULT_A4);
    DPMsteelEC();
    ~DPMsteelEC();

    const char *getClassType(void) const {return "DPMsteelEC";};


    double getThermalElongation(void); //***JZ
    double getElongTangent(double, double&, double&, double);//***JZ //PK add to include max temp

    int setTrialStrain(double strain, double strainRate =0); 
    int setTrialStrain(double strain, double FiberTemperature, double strainRate); //***JZ

    int setTrial (double strain, double &stress, double &tangent, double strainRate = 0.0);
    double getStrain(void);              
    double getStress(void);
    double getTangent(void);
    double getInitialTangent(void) {return E0;};

    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);        

    UniaxialMaterial *getCopy(void);
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    
    void Print(OPS_Stream &s, int flag =0);

    int getVariable(const char *variable, Information &);
    
    // AddingSensitivity:BEGIN //////////////////////////////////////////
    int setParameter(const char **argv, int argc, Parameter &param);
    int    updateParameter          (int parameterID, Information &info);
    int    activateParameter        (int parameterID);
    double getStressSensitivity     (int gradIndex, bool conditional);
    double getInitialTangentSensitivity(int gradIndex);
    int    commitSensitivity        (double strainGradient, int gradIndex, int numGrads);
    // AddingSensitivity:END ///////////////////////////////////////////
    
 protected:
    
 private:
    
    void yield_funtion(double E0, double afun, double bfun, double cfun, double ep, double fp, double fy, double fu, double ey, double& r, double& G); // Anand Kumar IITJ
    //double yield_funtions(double E0, double afun, double bfun, double cfun, double ep, double fp, double fy, double fu, double ey);
    void Nr(double sg_trial1, double E, double a, double b, double c, double ep, double fp, double fy, double fu, double ey, double alp_p, double &stress1, double &E_tgt, double& eplas1, double& alp1); // Anand Kumar IITJ

    //JZ 07/10 /////////////////////////////////////////////////////////////start
    double Temp;  // material temp  
    //double steps;    //the amount of the steps. 
    double ThermalElongation; // eps(theata) = alpha * temperature  
    double fyT;
    double E0T;
    double fp; // 11/10
    double fu; // Added by Anand Kumar IITJ
    double eplasP;
    double eplas;
    double alpP;
    double alp;
    //int position;
	
	
    //JZ 07/10 /////////////////////////////////////////////////////////////end 
    
    
    /*** Material Properties ***/
    double fy;  // Yield stress
    double E0;  // Initial stiffness
    double b;   // Hardening ratio (b = Esh/E0)
    double a1;
    double a2;
    double a3;
    double a4;  // a1 through a4 are coefficients for isotropic hardening
    
    /*** CONVERGED History Variables ***/
    double CminStrain;  // Minimum strain in compression
    double CmaxStrain;  // Maximum strain in tension
    double CshiftP;     // Shift in hysteresis loop for positive loading
    double CshiftN;     // Shift in hysteresis loop for negative loading
    int Cloading;       // Flag for loading/unloading
                        // 1 = loading (positive strain increment)
                        // -1 = unloading (negative strain increment)
                        // 0 initially
    bool Cmono;   //state of monotonic loading,added by liming,2013
	
    /*** CONVERGED State Variables ***/    
    double Cstrain;
    double Cstress;
    double Ctangent;    
    double Ctemperature;
    double Cmodulus; //added by Princeton
	
    /*** TRIAL History Variables ***/
    double TminStrain;
    double TmaxStrain;
    double TshiftP;
    double TshiftN;
    int Tloading;
    bool Tmono;   //state of monotonic loading,added by liming,2013
	
    /*** TRIAL State Variables ***/
    double Tstrain;
    double Tstress;
    double Ttangent; // Not really a state variable, but declared here
                     // for convenience
    double Ttemperature; //(Trial) Temperature,added by liming,2013
    double Tmodulus; //added by Princeton
    // Calculates the trial state variables based on the trial strain
    void determineTrialState (double dStrain);

    // Determines if a load reversal has occurred based on the trial strain
    void detectLoadReversal (double dStrain);

// AddingSensitivity:BEGIN //////////////////////////////////////////
    int parameterID;
	Matrix *SHVs;
// AddingSensitivity:END ///////////////////////////////////////////
};

#endif
