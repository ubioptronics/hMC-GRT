#include "trmc_2fls_opts.h"

#ifndef __TRMC_2FLS_GLOBALS_H__
#define __TRMC_2FLS_GLOBALS_H__

//***************************************************

#define SOURCE_FIBER_RADIUS 0.03
	// controls the (x,y) coordinates of photon launch in 
	// function initPhoton(...) 

//#define NORMALPHOTONINCIDENCE 
	// do not define this if photon must be launched corresponding to 
	// the Numerical Aperture of the fiber. See - initPhoton(..)

#define NUMERICAL_APERTURE 1.45
	//for full fiber-surface detection, set NA = 1

#define TERMINALPHOTONWT 1E-5
	//simulation terminates photon if weight is below this. 
#define SURVIVALCHANCE 50
	// defines one in "how many" photons survive the roulette 
	// is neccessary if you want to use the roulette routines. 
#define SPEEDOFLIGHT 3E10 	//cm/s
#define TIMERES 1E-12	  //set to appropriate abs. value of delta-t bins.
#define TIMEWINDOW 2500		//set number of delta-t bins 
#define BEERS_LAW_MODEL
#define USE_RAN2
	// if USE_RAN3 defined, uses ran3() to generate random #s, 
	// if not defined, (or is commented) uses ran2()
#define PRINT_MC_RUNNINGINFO
	// define to control creation of photon.number indicating 
	// current photon number in the doAnisotropicScatter() loop

//#define MC_TRACK_FL_PHOTONS_ORIGIN
	// define this if you want FL*origin.* to be created. 

#define MC_VERSION_STR "1.0"

//#define SOURCE_FIBER_IN     
	//define when source fiber is used as a source. Do not define if ZEMAX input is used
#define ZEMAX_INPUT          
	//define when ZEMAX array data is used as an input. Do not define if source fiber is used
#define NUMOFRAYS 999928	 //number of rays in ZEMAX input array data

/***************************************************/
#define NUMBER_DETECTOR_POSITIONS (sizeof(detectorradius)/sizeof(double))
	//the order of definitions do not matter, as gcc makes multiple passes to compile... 

double detectorradius[] = {0.03, 0.09, 0.15, 0.21, 0.27, 0.33}; 
long idum;	// seed for ran3. */
double MAX_ACCEPTANCE_ANGLE = 0.0;
double ***flphotoncurrent = NULL;
double diffuse_tra, diffuse_ref;
double photoncurrent[2][TIMEWINDOW][NUMBER_DETECTOR_POSITIONS];
double zemax_array[NUMOFRAYS][7]; 

char fname_zemaxin[100] = "4936_EX_1mmoff_0.dat";   //put name of zemax input array file  
char fname_zemaxout[100] = "Exited_photons.dat";    //select name of zemax output array file 

//***************************************************/
#endif
