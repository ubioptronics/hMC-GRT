#ifndef __TRMC_2FLS_BASE_H__
#define __TRMC_2FLS_BASE_H__


typedef struct 
{
	double x, y, z;
	double prevx, prevy, prevz;
	double fl_x, fl_y, fl_z, fl_timecreation;
	double cosalpha, cosbeta, cosgamma;
	double scatterlength;
	int tag;	//1 => excitation photon, 2=> fl. photon.
	double oldcosalpha, oldcosbeta, oldcosgamma;
	int isotropicremission;
	int photonalive; 
	double currentweight;
	double timeofflight;
	double decaytime;	
		// used to store the generated "_tau_" for this
		// photon. Set to 0 in initphoton..
}photon;

typedef struct
{
	int n_lambda;
	double **opt_trans_coeffs;
		// this is a n_lambdax4 matrix. The first index refers to the 
		// wavelength-tag, the seconnd access the elements. 
		// MUA = 0, MUS = 1, G_VAL = 2, REF_IND = 3;
	double **muafs, **taus, **flqys;
	// these must be a 2-d arrays. Each of the arrays is a FxF square matrix
	// where F = n_lambda-1; This matrix is only filled along the diagnol and 
	// up due to the nature of fluorescence reabsorption. The parser
	// assures that this will be done. Read ... 
	double thickness, z_start, z_end;
} layer;

typedef struct tissue_struct *TissuePtr;
typedef struct tissue_struct
{
	long numphotons; //number of photons
	double n_top, n_bottom; //refractive indicies above and below the tissue
	int n_lambda;	//total count of lambdas capable of propogation in this tissue
	double *lambdas;	// the values of the lamdas. *lambda <-> lambda_ex
	int numlayers;	// total number of layers
	layer *layerptr;	//array of pointers to the layers
	TissuePtr next_tissue; //make a linked list... 
} tissue;

//***************************************************/
void doAnisotropicScatter(long num_p, tissue const *tp);
void allocateStorageBins(const tissue *tp);
void initPhoton(photon *p,long num_p);
double getRandomNumber();
layer *getCurrentLayer(const tissue *tp, double p_z);
int getCurrentLayerNumber(const tissue *tp, double p_z);
void path(photon *p, const layer *lp, const tissue *tp);
void fluorescenceEvent(photon *p, layer *lp);
void scatter(photon *p, const layer *lp);
int score(const photon *pp, const tissue *tp);
void dumpPhotonIntensity(const tissue *tp);
double getLayeropt_coeff(const layer *lp, int r, int c);
void printPhoton(photon p);
double findFresnelReflectionCoef(double n1, double cosgamma, double n2);
double getRefractedAngle(double n1, double cosgamma, double n2);
float ran3(long *idum);
float ran2(long *idum);
void readZEMAXinput(char *fname);
void saveExitedPhotons(const photon *pp);

// functions from parsefile.c
int readLines(FILE *fp, char *input);
tissue *parseFile(int argc, char **argv);
int getNumberPhotons(const char *str, tissue *tp);
int getTopBottRefIndex(const char *str, tissue *tp);
int getn_lambdas(const char *str, tissue *tp);
int getlambdaVals(const char *str, tissue *tp);
int getnumLayers(const char *str, tissue *tp);
int getopt_Coeffs(const char *str, double *dp);
int getfluor_Coeffs(const char *str, layer *lp, int r, int c);
int getLayerThickness(const char *str, tissue *tp, int l_num);
void fprintTissue(FILE *fp, const tissue *tp);
void fprintLayer(FILE *fp, const layer *lp, const double *lambdas);
double getTotalmuaf(const layer *lp, int tag);
//***************************************************/
#endif
