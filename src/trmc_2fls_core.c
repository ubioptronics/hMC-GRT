// The multi-fluorophore version of the time-resolved fluorescence 
// monte carlo code originally written by Karthik Vishwanath. 
// The primary purpose was to model the presence of varying 
// concentrations of *physically different* fluorophores 
// in the same slab (layer) of the tissue.  
// 
// Interface function with ZEMAX ray tracing software was added by Paul Lee 


/***************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
//#include <sys/utsname.h>

#include "trmc_2fls_base.h"
#include "trmc_2fls_globals.h"
//***************************************************/
int main(int argc, char **argv)
{
	tissue *t_head, *tp;
	int numtissues = 0;
	time_t simstarttime, simendtime, tissuestarttime, tissueendtime;
	double timediff;
	int mins, hours, secs;
	
	
	t_head = tp = parseFile(argc, argv);
	simstarttime = time((time_t *)NULL);

	while (tp != NULL) {
		//srand((unsigned int)time((time_t *)NULL));
    		idum = -(int)time(NULL)%(1<<15);	//dummy var for ranN()
		//idum = -1000;
		tissuestarttime = time((time_t *)NULL);        
		fprintf(stdout, "Started tissue#%d: %s", numtissues, ctime(&tissuestarttime)); 
		doAnisotropicScatter(tp->numphotons, tp);
		//fprintTissue(stdout, tp); //print tissue details...
		tissueendtime = time((time_t *)NULL);        
		fprintf(stdout, "Ended tissue#%d: %s", numtissues, ctime(&tissueendtime)); 
		timediff = (double) (tissueendtime-tissuestarttime);
		mins = (int) (timediff/60.0)%60;
		hours = floorl(timediff/3600);
		secs = ((long)timediff)%60;
		fprintf(stdout, "Time taken to simulate tissue %d: %.2d:%.2d:%.2d (hh:mm:ss) \n", numtissues++, hours, mins, secs);
		fprintf(stdout, "--------------------------------------\n");
		tp = tp->next_tissue;
	}
	simendtime = time((time_t *)NULL);
	fprintf(stdout, "\n\nSimulations started:  %s", ctime(&simstarttime)); 
		//tell when it started.
	fprintf(stdout, "Simulation completed: %s", ctime(&simendtime)); 
		//tell when it ended

	timediff = (double) (simendtime - simstarttime);
	mins = (int) (timediff/60.0)%60;
	hours = floorl(timediff/3600);
	secs = ((long)timediff)%60;
	fprintf(stdout, "total time taken to complete %d tissues: %.2d:%.2d:%.2d (hh:mm:ss) \n", numtissues, hours, mins, secs);

	return(EXIT_SUCCESS);
} //int main()
//***************************************************/
void doAnisotropicScatter(long num_p, tissue const *tp)
{
	long nump;
	long exphotoncount = 0;
	long flphotoncount = 0;
	photon p1;
	layer *currentlayer;
	double maxdepth = 0.0;	//heights of all layers summed. 
	FILE *fp;
	
#ifdef PRINT_MC_RUNNINGINFO
	if ((fp = fopen("photon.number", "w")) == NULL)
	{
		perror("Cannot open photon.number in doAnisotropicScatter- quitting. ");
		exit(EXIT_FAILURE);
	}
#endif 
	maxdepth = ((tp->layerptr)[tp->numlayers-1]).z_end;
	fprintTissue(stdout, tp); //print tissue details...
	MAX_ACCEPTANCE_ANGLE = asin(NUMERICAL_APERTURE/tp->n_top);
	allocateStorageBins(tp);
		//allocate and clear photon storage bins...
		
#ifdef ZEMAX_INPUT
	memset(zemax_array, 0, sizeof(double)*tp->numphotons*6);
	readZEMAXinput(fname_zemaxin); 
#endif

	for (nump=0; nump<tp->numphotons; nump++) {	
#ifdef PRINT_MC_RUNNINGINFO
		rewind(fp);
		fprintf(fp, "%ld\n", nump);
#endif 
		fflush(NULL);
		initPhoton(&p1,nump);
		currentlayer = getCurrentLayer(tp, p1.z);
		while (p1.z < maxdepth && p1.z >= 0.0 && p1.photonalive) {
			path(&p1, currentlayer, tp); 
				// is what changes p1.x,y,z, scatterlength
				// updates weight, timeofflight
			if (tp->n_lambda > 1) // do fluorescence only if specified...
				fluorescenceEvent(&p1, currentlayer);
				// generate Fl. event.. 
			scatter(&p1, currentlayer); 
				// is what changes p1.cosalpha, beta gamma.
				// check is photon is healthy enough to propogate ;-)
			currentlayer = getCurrentLayer(tp, p1.z);
				//update layer 
			//printPhoton(p1);
			if (p1.currentweight < TERMINALPHOTONWT) {
				if (getRandomNumber() > 1/SURVIVALCHANCE)
					p1.photonalive = 0;
				else
					p1.currentweight *= SURVIVALCHANCE;
			}
		} // while (p1.z < maxdepth && p1.z >= 0.0 && p1.photonalive)

		if (p1.z <= 0.0 && p1.photonalive){ 	// exited from the top surface and is alive!
			if (p1.tag == 0)
				diffuse_ref += p1.currentweight;
			if  (score(&p1, tp))
				(p1.tag)? flphotoncount++: exphotoncount++;
			saveExitedPhotons(&p1);     // save the exited photons for ZEMAX output
		}

		if (p1.photonalive && !(p1.tag) && p1.z > maxdepth) 
			diffuse_tra += p1.currentweight;
	} // for (nump=0; nump<tp->numphotons; nump++) 	

#ifdef PRINT_MC_RUNNINGINFO
	fclose(fp);
#endif 

	fprintf(stdout, ">> # of EX_photons detected = %0.2E\n", (double)exphotoncount);
	fprintf(stdout, ">> # of FL_photons detected = %0.2E\n", (double)flphotoncount);

	dumpPhotonIntensity(tp);

	//printf("*************Done for tissue************* \n");
	
} //void doAnisotropicScatter(tissue t)
//***************************************************/
void allocateStorageBins(const tissue *tp)
{
	int i,j;

	memset(photoncurrent, 0, sizeof(double)*2*NUMBER_DETECTOR_POSITIONS*TIMEWINDOW);
	diffuse_tra = 0.0;
	diffuse_ref = 0.0;
	if (tp->n_lambda >1) {
		MALLOC(flphotoncurrent, sizeof(double **)*(tp->n_lambda-1));
		for (i=0; i< tp->n_lambda-1; i++) {
			MALLOC(flphotoncurrent[i], sizeof(double *)*TIMEWINDOW);
			for (j=0; j<TIMEWINDOW; j++) {
				MALLOC(flphotoncurrent[i][j], sizeof(double )*NUMBER_DETECTOR_POSITIONS);
				memset(flphotoncurrent[i][j], 0, sizeof(double)*NUMBER_DETECTOR_POSITIONS);
			} // for (j=0; j<TIMEWINDOW; j++) 
		} //for (i=0; i< tp->n_lambda-1; i++)
	} //if (tp->n_lambda >1) 

} //void allocateStorageBins(const tissue *tp)
//***************************************************/
void initPhoton(photon *p,long num_p)
{

	double temp1, temp2; 
#ifdef SOURCE_FIBER_INPUT
	
	temp1 = 2*M_PI*getRandomNumber();
	temp2 = SOURCE_FIBER_RADIUS*sqrt(getRandomNumber());
	p->x = cos(temp1)*temp2;
	p->y = sin(temp1)*temp2;
	p->z = 0.0f;
#endif

#ifdef NORMALPHOTONINCIDENCE
    p->cosgamma = 1.0f; 
	p->cosalpha = p->cosbeta = 0.0f;
#else
	p->cosgamma = cos(MAX_ACCEPTANCE_ANGLE*getRandomNumber() );
	temp1 = sqrt(1 - p->cosgamma*p->cosgamma); //sin(gamma)
	temp2 = 2*M_PI*getRandomNumber();
	p->cosalpha = temp1*cos(temp2);
	p->cosbeta = temp1*sin(temp2);
#endif

#ifdef ZEMAX_INPUT
	p->x = 0.1*zemax_array[num_p][0];
	p->y = 0.1*zemax_array[num_p][1];
	p->z = 0.0f;

	p->cosalpha = zemax_array[num_p][3];
    p->cosbeta = zemax_array[num_p][4];
    p->cosgamma = zemax_array[num_p][5]; 
#endif 

	p->tag = 0;		//start as excitation photon.
		// the tag can maximally be n_lambda-1...

//	p->scattercount_ex = 0; 
//	p->scattercount_fl = 0; 
//	
//	p->firstfl_event = 0;
	
	p->photonalive = 1;		//weight tracking
	p->currentweight = 1.0;
	p->timeofflight = 0.0;
	p->decaytime = 0.0;
	p->fl_x = 0.0;
	p->fl_y = 0.0;
	p->fl_z = 0.0;
	p->fl_timecreation = 0.0;
	p->isotropicremission = 0;

} //void initPhoton(photon *p)
//***************************************************/
layer *getCurrentLayer(const tissue *tp, double p_z)
{
	//depending on photon's current position (in p_z)
	// decide which layer the photon lies in. 
	
	int i;
	
	if (p_z < 0)
		return(NULL); 	// will return NULL if the photon 
						// has exited from the top

	for (i=0; i<tp->numlayers; i++)
		if(((tp->layerptr)[i]).z_end > p_z)	
			return(&(tp->layerptr)[i]);	// return currentlayer.
			// asumming layers are stacked along increasing z, 
			// the above conditon alone is sufficient to determine 
			// the currentlayer.
			
	return(NULL);	// will return NULL if the photon 
					//has exited from the bottom

} // layer *getCurrentLayer(const tissue *tp, double p_z)
//***************************************************/
int getCurrentLayerNumber(const tissue *tp, double p_z)
{
	// essentially a copy of getCurrentLayer(),
	// modified to return the layer#. 
	// return value is NON-ZERO for meaningful input 
	// (i.e. photon lies within the tissue medium.)
	// this function is used to determine which fluorophore (in 
	// which layer) gave rise to the fluorescent photon
	
	int i;
	
	if (p_z < 0)
		return(0); 	// will return NULL if the photon 
						// has exited from the top

	for (i=0; i<tp->numlayers; i++)
		if(((tp->layerptr)[i]).z_end > p_z)	
			return(i+1);	// return currentlayer.
			// asumming layers are stacked along increasing z, 
			// the above conditon alone is sufficient to determine 
			// the currentlayer.
			
	return(0);	// will return NULL if the photon 
					//has exited from the bottom

} // int getCurrentLayerNumber(const tissue *tp, double p_z)
//***************************************************/
void path(photon *p, const layer *lp, const tissue *tp)
{
	double current_mua, current_mus;
	double stepx, stepy, stepz;
	double xc, yc, zc, path_in_tissue;
	double z_boundary, next_mua, next_mus;
	double step1, step2;
	double current_refind, next_refind;
	layer *nextlp;

	current_mua = getLayeropt_coeff(lp, p->tag, MUA);
	current_mus = getLayeropt_coeff(lp, p->tag, MUS);
	current_refind = getLayeropt_coeff(lp, p->tag, REFIND);

	p->prevx = p->x;
	p->prevy = p->y;
	p->prevz = p->z;

#ifdef BEERS_LAW_MODEL
	p->scatterlength = -log(getRandomNumber())/(current_mus);
#else
	p->scatterlength = -log(getRandomNumber())/(current_mua+current_mus);
#endif

	stepx = p->x + p->scatterlength*p->cosalpha;	
		//move photon in x to "expected" position
	stepy = p->y + p->scatterlength*p->cosbeta;	
		//move photon in y to "expected" position
	stepz = p->z + p->scatterlength*p->cosgamma;	
		//move photon in z to "expected" position
	nextlp = getCurrentLayer(tp, stepz);	//check if layers have been crossed	

	if (nextlp == lp) { // photon is in same layer, 
		p->x = stepx;
		p->y = stepy;
		p->z = stepz;
		p->timeofflight += p->scatterlength/((double)(SPEEDOFLIGHT/current_refind));
		
#ifdef BEERS_LAW_MODEL
		p->currentweight *= exp(-p->scatterlength*current_mua);
#else
		p->currentweight *= current_mus/(current_mus+current_mua);
#endif
		return;	// the code does not *require* this... better safe than sorry.
	} //if (nextlp == lp)

	if (nextlp == NULL) {
		// lets detect photon exit first...
		//photon has EXITED the tissue (either at the top or bottom) and
		// total internal reflectance checking must be done...
	
		if (p->cosgamma > 0) {	// photon is exitting from the bottom 
			z_boundary = lp->z_end;	
			next_refind = tp->n_bottom;
		}
		else {	// photon is exitting from the top
			z_boundary = lp->z_start;	
			next_refind = tp->n_top;
		}

		if (getRandomNumber() <= findFresnelReflectionCoef(current_refind, p->cosgamma, next_refind)) {	//if photon is reflected
			p->x = stepx;
			p->y = stepy;
			p->z = 2*z_boundary - stepz;
			p->cosgamma *= -1;
			p->timeofflight += p->scatterlength/((double)(SPEEDOFLIGHT/current_refind));					
#ifdef BEERS_LAW_MODEL
			p->currentweight *= exp(-p->scatterlength*current_mua);
#else
			p->currentweight *= current_mus/(current_mus+current_mua);
#endif
			}
		else { 	
		// normal photon transmission -- will lead to exit from tissue...
		//note that code may have to be mod wrt calculation of 
		// timeofflight for **this** step, as it makes the photon escape
			p->x += p->scatterlength*p->cosalpha;
			p->y += p->scatterlength*p->cosbeta;
			p->z += p->scatterlength*p->cosgamma;
			p->timeofflight += p->scatterlength/((double)(SPEEDOFLIGHT/current_refind));
			// for BEER's law model, the photon weight attenuation must depend on the 
			// exact path length to the surface... 
			xc = (z_boundary-p->prevz)*(p->x-p->prevx)/(p->z-p->prevz);
			yc = (z_boundary-p->prevz)*(p->y-p->prevy)/(p->z-p->prevz);
			zc = (z_boundary-p->prevz);
			path_in_tissue = sqrt(xc*xc + yc*yc + zc*zc);
#ifdef BEERS_LAW_MODEL
			p->currentweight *= exp(-path_in_tissue*current_mua);
#else
			p->currentweight *= current_mus/(current_mus+current_mua);
#endif
		} //else // normal photon transmission -- will lead to exit from tissue
		return;	// the code does not *require* this... better safe than sorry.
	} //if (nextlp == NULL)
	else {	
		//i.e. nextlp != lp 
		// layers have been crossed... photon NOT in air 
		next_mua = getLayeropt_coeff(lp, p->tag, MUA);
		next_mus = getLayeropt_coeff(lp, p->tag, MUS);
		next_refind = getLayeropt_coeff(lp, p->tag, REFIND);

		if (next_refind != current_refind) {
			// layers have diff. ref_indicies ... lots of processing...
			if (p->cosgamma > 0)
				z_boundary = lp->z_end;	
			else
				z_boundary = lp->z_start;	

			if (getRandomNumber() <= findFresnelReflectionCoef(current_refind, p->cosgamma, next_refind)) {
				//photon is reflected
				p->x = stepx;	//symmetry on z-reflection will keep x unchanged
				p->y = stepy;	//symmetry on z-reflection will keep y unchanged
				p->z = 2*z_boundary - stepz;	//calculated as z_boundary + (z_boundary-stepz)
				p->cosgamma *= -1;	//reverse direction due to reflection
				p->timeofflight += p->scatterlength/((double)(SPEEDOFLIGHT/current_refind));
#ifdef BEERS_LAW_MODEL
				p->currentweight *= exp(-p->scatterlength*current_mua);
#else
				p->currentweight *= current_mus/(current_mus+current_mua);
#endif
				return;
			} // if (getRandomNumber() <= findFresnelReflectionCoef(current_refind, ...)
			else { 	// else it is transmitted	
				step1 = (p->scatterlength - fabs(z_boundary - p->z)/fabs(p->cosgamma));
				//foreshortened length in current layer
				step2 = (p->scatterlength - step1)*current_mus/next_mus;
				//adjusted (remaining) length in next layer
				stepx = p->x + step1*p->cosalpha;
				stepy = p->y + step1*p->cosbeta;
				stepz = p->z + step1*p->cosgamma;
				p->cosalpha *= current_refind/next_refind;
				p->cosbeta *= current_refind/next_refind;
				p->cosgamma = (p->cosgamma/fabs(p->cosgamma))*cos(getRefractedAngle(current_refind, p->cosgamma, next_refind));
				p->x = stepx + step2*p->cosalpha;
				p->y = stepy + step2*p->cosbeta;
				p->z = stepz + step2*p->cosgamma;
				p->timeofflight += step1/((double)(SPEEDOFLIGHT/current_refind)) + step2/((double)(SPEEDOFLIGHT/next_refind));
#ifdef BEERS_LAW_MODEL
                                p->currentweight *= exp( -(step1*current_mua + step2*next_mua));
#else
                                p->currentweight *= (next_mua/(next_mua+next_mus))* (current_mus/(current_mus+current_mua));
#endif
				return;
			}	//else (photon transmitted)
		} // if (next_refind != current_refind) 
		else {
			// layers have been crossed, and no refractive index mismatch. 
			// the calculations here are slightly different -- 
			// must attenuate photons' weight by actual path lengths over
			// *every* layer it traverses. Theres no need for adjusting step-sizes
			// per layer. Its quite conceivable here to have a photon cross two layer 
			// boundaries. The code does not perform the "correct" calculations for 
			// cases of multiple layer crossings, but for now, prints a diagnostic message 
			// so that we may atleast know how "frequent" such crossings are... 

			if (p->cosgamma > 0)
				z_boundary = lp->z_end;	
			else
				z_boundary = lp->z_start;	
			step1 = (p->scatterlength - fabs(z_boundary - p->z)/fabs(p->cosgamma));
				// foreshortened length in current layer
			step2 = (p->scatterlength - step1)*current_mus/next_mus;
				// adjusted (remaining) length in next layer
			p->timeofflight += step1/((double)(SPEEDOFLIGHT/current_refind)) + step2/((double)(SPEEDOFLIGHT/next_refind));
			
#ifdef BEERS_LAW_MODEL
			p->currentweight *= exp( -(step1*current_mua + step2*next_mua));
#else
			p->currentweight *= (next_mua/(next_mua+next_mus))* (current_mus/(current_mus+current_mua));
#endif
					
			// Move the photon to where its supposed to go...
			p->x += p->scatterlength*p->cosalpha;
			p->y += p->scatterlength*p->cosbeta;
			p->z += p->scatterlength*p->cosgamma;
			return;
		}// else [if (next_refind != current_refind) ]
	} // else [if (nextlp != lp)]

	printf("Can never reach here. Must abort by definition!\n");
	perror(".. inside void path(...)\n");
	exit(EXIT_FAILURE);

} //void path(photon *p, layer *lp)
//***************************************************/
void fluorescenceEvent(photon *p, layer *lp)
{
	double temp, eff_muaf, norm_fac;
	int newtag, oldtag, k;

	oldtag = p->tag+1; //the tag relates the the muafs/taus/flqys array by 1+tag...
	if (oldtag == lp->n_lambda) 	//the last wavelength ?
		return;		//can't have any more reabsorptions, so return
	
	eff_muaf = getTotalmuaf(lp, oldtag);
	temp = 1 - exp(-eff_muaf*(p->scatterlength));
	if (getRandomNumber() > temp)
		return; //no fluorescence, carry on with life..

	for (newtag=oldtag; newtag<lp->n_lambda-1; newtag++){	
		//start from tag+1 loop till n_lambda-1
		norm_fac = 0.0;
		for (k=newtag; k<lp->n_lambda; k++)	//sum the muafs starting from the newtag'th col
			norm_fac += lp->muafs[oldtag][k];
		temp = lp->muafs[oldtag][newtag]/norm_fac;
		if (getRandomNumber() < temp)
			break;
	} // for (newtag=1; newtag<lp->n_lambda-1; newtag++){
	
	p->tag = newtag;	//it is assured that newtag > old value of p->tag...
	// store the spatial origin of the fl-photon... 
	if (p->z > 0) {
		p->fl_x = p->x;
		p->fl_y = p->y;
		p->fl_z = p->z;
	}
	else {
		p->fl_x = p->prevx;
		p->fl_y = p->prevy;
		p->fl_z = p->prevz;
	}
	// store the temporal origin of the fl-photon... 
	p->fl_timecreation = p->timeofflight;
	p->isotropicremission = 1;
	//set this here. Will be detected and unset on the next run of scatter(...) 
	p->currentweight *= (lp->flqys[oldtag][newtag]);
	p->decaytime +=  -log(getRandomNumber())*(lp->taus[oldtag][newtag]);
	//p->decaytime += lp->taus[oldtag][newtag];

} // void fluoroscenceEvent(photon *p, layer *lp)
//***************************************************/
void scatter(photon *p, const layer *lp)
{
	double temp, current_g, cost, sint, cospsi, sinpsi, psi;

	current_g = getLayeropt_coeff(lp, p->tag, ANIS_G);

	p->oldcosalpha = p->cosalpha; //save for later calculations
	p->oldcosbeta = p->cosbeta;   //save for later calculations
	p->oldcosgamma = p->cosgamma;
		// we want the "old" value of cosgamma, because thats what 
		// made the photon exit the tissue. 

 	if (p->isotropicremission) {
		//A fl. event orients the photon randomly. Choose cosalpha, cosbeta
		//randomly. 
		p->isotropicremission = 0;
		temp = M_PI*getRandomNumber();  //theta, same thing..
		psi = 2*M_PI*getRandomNumber();
		p->cosalpha = sin(temp)*cos(psi);
		p->cosbeta = sin(temp)*sin(psi);
		p->cosgamma = cos(temp);
		return;
	}

	if (current_g == 0.0)
		cost = 2*getRandomNumber()- 1;
	else
	{
		temp = (1 - current_g*current_g)/(1 - current_g + (2*getRandomNumber()*current_g));
		cost = (1 + current_g*current_g - temp*temp)/(2*current_g);
	}
	
	sint = sqrt(1 - cost*cost);
	psi = 2*M_PI*getRandomNumber();
	cospsi = cos(psi);
	
	if (psi < M_PI)
		sinpsi = sqrt(1- cospsi*cospsi);
	else
		sinpsi = -sqrt(1- cospsi*cospsi);
	
	if (fabs(p->cosgamma) > 0.99999) {
		p->cosalpha = sint * cospsi;
		p->cosbeta = sint * sinpsi;
		if (p->cosgamma > 0)
			p->cosgamma = cost;
		else
			p->cosgamma = -cost;
	} //if (fabs(p->cosgamma) > 0.99999) {
	else {
		temp = sqrt(1 - p->cosgamma*p->cosgamma);
		p->cosalpha = sint*(p->oldcosalpha*p->cosgamma*cospsi - p->oldcosbeta*sinpsi)/temp + p->oldcosalpha*cost;
		p->cosbeta = sint*(p->oldcosbeta*p->cosgamma*cospsi + p->oldcosalpha*sinpsi)/temp + p->oldcosbeta*cost;
		p->cosgamma = p->cosgamma*cost - sint*cospsi*temp;
	}
} // void scatter(photon *p, layer *lp)
//***************************************************/
int score(const photon *pp, const tissue *tp)
{
	// return 0 if photon not detected
	// return 1 if the photon is detected
	double xintersect, yintersect, rho;
	int i, mytag; 
	long istot;
	double lambda;
	double collector_x, collector_y; // x,y coordinates of collection fiber
	double newx, newy;
	FILE *fp = NULL;
	char fname[100];

	lambda = -pp->z/(pp->prevz - pp->z);
	xintersect = pp->x + lambda*(pp->prevx - pp->x); 
		//x-coordinate of where photon path intersects the z=0 plane.
	yintersect = pp->y + lambda*(pp->prevy - pp->y);
		//y-coordinate of where photon path intersects the z=0 plane.
	collector_x = 0;
	collector_y = 0;

	newx = xintersect-collector_x; //shift origin to the collector's position
	newy = yintersect-collector_y; 
	rho = sqrt(newx*newx + newy*newy);
//	if (pp->tag)

	for (i=0; i<NUMBER_DETECTOR_POSITIONS; i++)	{
		// the value of pp->oldcosgamma here has the directional cosine of exit 
		// as calculated from the z-axis. Since, for now, all the fibers axes 
		// carry no tilt, we can calculate the sin of the exit angle pretty easily
		// from pp->oldcosgamma_100. In fact, pp->oldcosgamma MUST always be NEGATIVE
		// since according to our local coordinate system, increasing z is downward
		// into the tissue, and here we are only counting photons escaping the
		// surface (no transmittance, i.e.). Therefore the sin of the exit angle
		// is just -pp->oldcosgamma!! All we need to do therefore is to check
		// -pp->oldcosgamma is smaller than NA/top_layers_refractive_index to 
		// ensure detection by that fiber...
		if ( (rho <=detectorradius[i]) && (sqrt(1-pp->oldcosgamma*pp->oldcosgamma) <= NUMERICAL_APERTURE/tp->n_top) ) 	 {
			istot = (long)floor((pp->timeofflight + pp->decaytime)/TIMERES);
			if (istot > TIMEWINDOW-1)	//throw away these photons
				istot = TIMEWINDOW-1;
			mytag = (pp->tag) ?  1 : 0; //if pp->tag is > 0, then set mytag = 1;
			photoncurrent[mytag][istot][i] += pp->currentweight;
			if (mytag) { // if theres fluorescence, update flcurrent for each wavelength (tag) ...
				flphotoncurrent[pp->tag-1][istot][i] += pp->currentweight;
#ifdef MC_TRACK_FL_PHOTONS_ORIGIN
				sprintf(fname, "FL_origin.%d", i);
				if ((fp = fopen(fname, "a")) == NULL) {
					perror("Cannot open FL_origin.%d for writing, quitting");
					exit(EXIT_FAILURE);
				}
				// lets keep track of the spatial origin of fluorescence, as well as the 
				// the "instant" of fluorescence, and the "time bin" value in seconds, 
				// rounded to the nearest TIMERES point...
				fprintf(fp, "%1.3g\t%1.3g\t %1.3g\t %1.3g\n", pp->fl_x, pp->fl_y, pp->fl_z, pp->fl_timecreation+pp->decaytime);
				fclose(fp);
#endif
			} //if (mytag) 
		//this photon has been detected 
		//and its weight recorded by the ith radius every radius > 
		//detectorradius[i] is going to pick it up, which we don't need 
		return(1); 
		} // if ( (rho <=detectorradius[i]) && ...)
	} // for (i=0; i<NUMBER_DETECTOR_POSITIONS; i++)

	return(0);
} // int score(photon p, tissue tp)
//***************************************************/
void dumpPhotonIntensity(const tissue *tp)
{
	//just dump the intensities to file, keeping time signatures...
	//do it for each detector_position 
	
	int i, j, k;
	FILE *fp;
	char *TRdumpfnames, *timestamp, *temp_floriginname, *real_floriginname;
//	struct utsname my_utsname;
	time_t current_t;
	struct tm *current_tm;
	
	MALLOC(TRdumpfnames, sizeof(char)*300);
	MALLOC(timestamp, sizeof(char)*300);
	MALLOC(temp_floriginname, sizeof(char)*300);
	MALLOC(real_floriginname, sizeof(char)*300);

	/*
	if (uname(&my_utsname)) {
		sprintf(my_utsname.sysname, "unknown");
		sprintf(my_utsname.nodename, "unknown");
		sprintf(my_utsname.release, "unknown");
		sprintf(my_utsname.version, "unknown");
		sprintf(my_utsname.machine, "unknown");
	}
	*/

	time(&current_t);
	current_tm = localtime(&current_t);
	sprintf(timestamp, "__%.2d%.2d%.2d__%.2d-%.2d-%d.dump", current_tm->tm_hour, current_tm->tm_min, current_tm->tm_sec,  current_tm->tm_mday, current_tm->tm_mon+1, current_tm->tm_year+1900);
	//sprintf(timestamp, "%ld.dump", (long) current_t);

	if ((fp = fopen("TR_dump.log", "a")) == NULL) {
		perror("montecarlo, quitting.\n");
		exit(EXIT_FAILURE);
	}
	fprintf(fp, "/************************************************/\n");
	fprintf(fp, "*                   NEW RUN ...                 */\n");
	fprintf(fp, "/************************************************/\n");
	fprintf(fp, "Multi-flurophore MC code: version %s\n", MC_VERSION_STR);
//	fprintf(fp, "Executed on %s; CPU = %s \n \tOS = %s (%s)\n", my_utsname.nodename, my_utsname.machine, my_utsname.sysname, my_utsname.release);
	fprintf(fp, "created on %s \n", ctime(&current_t));
	fprintf(fp, "names = TR_ex%s [TR_fl%s] \n", timestamp, timestamp);
	//fprintf(fp, "numberofphotons = %#1.4E\n", (double)tp->numphotons);
	fprintTissue(fp, tp);
	fprintf(fp, "number of detectors %d: \n  [", (int)NUMBER_DETECTOR_POSITIONS);
	for(i=0; i<NUMBER_DETECTOR_POSITIONS; i++)
		fprintf(fp, " %#1.3f; ", detectorradius[i]);
	fprintf(fp, "]; (cm)\ntime resolution: %#1.3E (ns)\n", TIMERES/1E-9);
	fprintf(fp, "SOURCE_FIBER_RADIUS = %#1.3E (cm)\n", (double)SOURCE_FIBER_RADIUS);
	fprintf(fp, "Numerical aperture = %#1.3E \n", NUMERICAL_APERTURE);
	fprintf(fp, "\tmax acceptance angle in air = %#0.3g [deg]\n", (180/M_PI)*asin(NUMERICAL_APERTURE));
#ifndef NORMALPHOTONINCIDENCE
	fprintf(fp, "Photons launched as distributed within angle (for fiber in refractive index of n_top (= %0.3f) )= %#1.2E [deg]\n", tp->n_top, MAX_ACCEPTANCE_ANGLE*180/M_PI);
#else
	fprintf(fp,"Photons launched normally (NORMALPHOTONINCIDENCE defined)\n");
#endif
	fprintf(fp, "Vel of light, c = %#1.2E (cm/s) \n", SPEEDOFLIGHT);
#ifndef BEERS_LAW_MODEL
	fprintf(fp, "photon weight attenuation using SLJ_MC_METHOD\n");
#else 
	fprintf(fp, "photon weight attenuation using BEERS LAW\n");
#endif

#ifndef USE_RAN3
	fprintf(fp, "using ran2() as generator \n");
#else 
	fprintf(fp, "using ran3() as generator \n");
#endif

#ifdef MC_TRACK_FL_PHOTONS_ORIGIN
	fprintf(fp, "MC_TRACK_FL_PHOTONS_ORIGIN defined\n");
#endif

	fprintf(fp, "Diffuse reflectance/transmittance: ");
	fprintf(fp, "R_d (ex/fl) = %#1.4f; T_d = %#1.4f\n", diffuse_ref/tp->numphotons, diffuse_tra/tp->numphotons);

	fclose(fp);

	for (i=0; i< NUMBER_DETECTOR_POSITIONS; i++) {
		// Dump excitation photons first...
		sprintf(TRdumpfnames, "TR_ex%s.%d", timestamp, i);
		if ((fp = fopen(TRdumpfnames, "w"))==NULL) {
			printf("cannot open %s. fatal error. quitting.\n", TRdumpfnames);
			exit(EXIT_FAILURE);
		}
		for(j=0; j<TIMEWINDOW-1; j++)	//last bin is temp storage.. so ignore.
			if (photoncurrent[0][j][i] > 0.0)
				fprintf(fp,"%#0.4E \t %#0.6E\n", j*(TIMERES/1E-9), photoncurrent[0][j][i]);
		fclose(fp);
		
		//... fluorescence photons if n_lambda > 1
		if (tp->n_lambda > 1){
			sprintf(TRdumpfnames, "TR_fl%s.%d", timestamp, i);
			if ((fp = fopen(TRdumpfnames, "w"))==NULL) {
				printf("cannot open %s. fatal error. quitting.\n", TRdumpfnames);
				exit(EXIT_FAILURE);
			}
			for(j=0; j<TIMEWINDOW-1; j++)	//last bin is temp storage.. so ignore.
				if (photoncurrent[1][j][i] > 0.0) {
					fprintf(fp,"%#0.4E \t %#0.6E ", j*(TIMERES/1E-9), photoncurrent[1][j][i]);
					for (k=0; k<tp->n_lambda-1; k++)
						fprintf(fp, "\t %#0.4E ", flphotoncurrent[k][j][i]);
					fprintf(fp, "\n");
				} //if (photoncurrent[1][j][i] > 0.0)
					
			fclose(fp);

#ifdef MC_TRACK_FL_PHOTONS_ORIGIN
			// all the FL_origin.%d files have the origin of the fl-photons detected
			// at each detector position, %d. Rename these to the appropriate
			// FL_*.origin.%d
			sprintf(temp_floriginname, "FL_origin.%d", i);
			sprintf(real_floriginname, "FL%s.origin.%d", timestamp, i);
			if (rename(temp_floriginname, real_floriginname)) {
				fprintf(stderr, "old name: %s ----> New name: %s \n", (char *)temp_floriginname, (char *)real_floriginname);
				perror("Could not rename file. fl-photon tracking error (most likely, the old name does not exist) \n");
			} //if (rename(temp_floriginname, real_floriginname))
#endif
		} //if (tp->n_lambda > 1)
	} //for (i=0; i< NUMBER_DETECTOR_POSITIONS; i++)
	
	free(TRdumpfnames);
	free(timestamp);
	free(temp_floriginname);
	free(real_floriginname);
} // void dumpPhotonIntensity(const tissue *tp)
//***************************************************/
void printPhoton(photon p)
{
	printf("\nx=%3.3g y=%3.3g z=%3.3g \n", p.x, p.y, p.z);
	printf("prevx=%3.3g prevy=%3.3g prevz=%3.3g\n", p.prevx, p.prevy, p.prevz);
	printf("scatterlength = %3.3g\n",p.scatterlength);
	//printf("totexcitationlength = %3.3g, totfluorescencelength = %3.3g\n",p.totexcitationlength, p.totfluorescencelength);
	printf("cosalpha=%3.3g , cosbeta=%3.3g , cosgamma=%3.3g\n",p.cosalpha, p.cosbeta, p.cosgamma);
	//printf("tag = %d, EX-scattercount = %ld, FL-scattercount = %ld \n", p.tag, p.scattercount_ex, p.scattercount_fl);
	printf("weight = %E, time_of_flight = %E \n", p.currentweight, p.timeofflight);
} //void printPhoton(photon p)
//***************************************************/
double getRandomNumber()
{

#ifndef USE_RAN3
	return( (double)ran2(&idum) );
#else 
	return( (double)ran3(&idum) );
#endif
				
} //double getRandomNumber()
//***************************************************/
float ran3(long *idum)
{
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC 1.0E-9

  static int inext,inextp;
  static long ma[56];
  static int iff=0;
  long mj,mk;
  int i,ii,k;
  
  if (*idum < 0 || iff == 0) {
    iff=1;
    mj=MSEED-(*idum < 0 ? -*idum : *idum);
    mj %= MBIG;
    ma[55]=mj;
    mk=1;
    for (i=1;i<=54;i++) {
      ii=(21*i) % 55;
      ma[ii]=mk;
      mk=mj-mk;
      if (mk < MZ) mk += MBIG;
      mj=ma[ii];
    }
    for (k=1;k<=4;k++)
      for (i=1;i<=55;i++) {
	ma[i] -= ma[1+(i+30) % 55];
	if (ma[i] < MZ) ma[i] += MBIG;
      }
    inext=0;
    inextp=31;
    *idum=1;
  }
  if (++inext == 56) inext=1;
  if (++inextp == 56) inextp=1;
  mj=ma[inext]-ma[inextp];
  if (mj < MZ) mj += MBIG;
  ma[inext]=mj;
  return mj*FAC;

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC
} //float ran3(long *idum)
//***************************************************/
float ran2(long *idum) 
// long period (> 1e18) random number generator of l ecuyer with bays-durham shue 
// and added safeguards. returns a uniform random deviate between 0.0 and 1.0 
// (exclusive of the endpoint values). call with idum a negative integer to initialize; 
// thereafter, do not alter idum between successive deviates in a sequence. 
// rnmx should approximate the largest  oating value that is less than 1. 
//
//
// obtained from chapter-7-pdf file available online at:
// http://www.library.cornell.edu/nr/bookcpdf.html
// this code is of version 2.10 as mentioned in the preface at the same website
// and is revised 2002.
//
//	... and i do have a numerical recepies book purchased (so, don't believe am 
// 		breaking any (c) laws ;-)
//
{
#define im1 2147483563 
#define im2 2147483399 
#define am (1.0/im1)
#define imm1 (im1-1) 
#define ia1 40014 
#define ia2 40692
#define iq1 53668
#define iq2 52774 
#define ir1 12211 
#define ir2 3791 
#define ntab 32 
#define ndiv (1+imm1/ntab) 
#define eps 1.2e-7 
#define rnmx (1.0-eps)

 	int j; 
	long k; 
	static long idum2=123456789; 
	static long iy=0; 
	static long iv[ntab]; 
	float temp; 

	if (*idum <= 0) //initialize. 
	{ 
		if (-(*idum) < 1) *idum=1; // be sure to prevent idum = 0. 
			else *idum = -(*idum); 
		
		idum2=(*idum); 
		for (j=ntab+7;j>=0;j--)  //load the shue table (after 8 warm-ups). 
		{
			k=(*idum)/iq1; 
			*idum=ia1*(*idum-k*iq1)-k*ir1; 
			if (*idum < 0) *idum += im1; 
			if (j < ntab) iv[j] = *idum; 
		} 
		iy=iv[0]; 
	} 
	k=(*idum)/iq1; //start here when not initializing. 

	*idum=ia1*(*idum-k*iq1)-k*ir1; 	//compute idum=(ia1*idum) % im1 without 
	if (*idum < 0) *idum += im1; 	//overflows by schrage s method. 
	k=idum2/iq2; 
	idum2=ia2*(idum2-k*iq2)-k*ir2;	//compute idum2=(ia2*idum) % im2 likewise. 
	if (idum2 < 0) idum2 += im2; 
	j=iy/ndiv; 		//will be in the range 0..ntab-1. 
	iy=iv[j]-idum2; //here idum is shuffled, idum and idum2 are 
	iv[j] = *idum; 	// combined to generate output. 
	if (iy < 1) iy += imm1; 
	if ((temp=am*iy) > rnmx) 
	return rnmx; //because users don't expect endpoint values. 
	else return temp; 
#undef  im1
#undef  im2
#undef  am
#undef  imm1
#undef  ia1
#undef  ia2
#undef  iq1
#undef  iq2
#undef  ir1
#undef  ir2
#undef  ntab
#undef  ndiv
#undef  eps
#undef  rnmx
} //float ran2(long *idum) 
//***************************************************/
double findFresnelReflectionCoef(double n1, double cosgamma, double n2)
{	// photon moves from layer of refindex n1 -> n2. 
	// acos(fabs(cosgamma)) is the incident angle. 
	// returns the value of the internal (Fresnel) reflectance 
	// R(theta_i)
	// checks for total internal reflectance
	
	double theta_i, theta_r;
	double t1, t2;

	if (n1 < 1.0 || n2 < 1.0)
	{
		printf("findFresnelReflectionCoef(...): Refractive index lesser than unity! Must abort.\n");
		exit(2);
	}

	theta_i = acos(fabs(cosgamma));
	//printf("n1 = %3g, n2 = %3g, cosgamma = %3g\n", n1, n2, cos(theta_i));

	if ((n1 > n2)&&(theta_i >= asin(n2/n1)))
		return(1.0);
		// total internal reflection occurs

	theta_r = asin((n1/n2)*sin(theta_i));

	t1 = sin(theta_i - theta_r)/sin(theta_i + theta_r);
	t2 = tan(theta_i - theta_r)/tan(theta_i + theta_r);
	t1 *= t1;
	t2 *= t2;
	//printf("theta_i = %3g, theta_r = %3g; FresnelReflectionCoef = %3g\n", acos(fabs(cosgamma)), theta_r, 0.5*(t1+t2));
	return(0.5*(t1+t2));
}
//***************************************************/
double getRefractedAngle(double n1, double cosgamma, double n2)
{
	double theta_i;

	theta_i = acos(fabs(cosgamma));
	if ((n1 > n2)&&(theta_i >= asin(n2/n1)))
		// total internal reflection occurs
	{
		printf("\tn1 = %3g, n2 = %3g, cosgamma = %3g, theta_i = %3g \n", n1, n2, cosgamma, theta_i);
		printf("\nShould not be inside getRefractedAngle(...) with this condition matched. Aborting.\n");
		perror("Should not be inside getRefractedAngle(...) with this condition matched. Aborting.");
		exit(2);
	} //if ((n1 > n2)&&(theta_i >= asin(n2/n1)))
			
	return(asin((n1/n2)*sin(theta_i)));
} //double getRefractedAngle(double n1, double cosgamma, double n2)
//***************************************************/
void readZEMAXinput(char *fname)
{	//read ZEMAX input array for initial photon distribution 

	int n = 0;
 	FILE *fp;

    if ((fp = fopen(fname, "r"))==NULL) {
			printf("cannot open ZEMAX Input File. fatal error. quitting.\n");
			exit(EXIT_FAILURE);
		}

	while ((!feof(fp)))
	{
		fscanf(fp,"%lf%lf%lf%lf%lf%lf%lf\n", &zemax_array[n][0], &zemax_array[n][1], &zemax_array[n][2], &zemax_array[n][3], &zemax_array[n][4], &zemax_array[n][5], &zemax_array[n][6]);
		++n;
	}
	fclose(fp);
}

void saveExitedPhotons(const photon *pp)
{   //save geometrical information of the photons exiting the surface
	double xintersect, yintersect;
	double lambda;
	FILE *fp = NULL;
	int mytag;
	double wl;

	lambda = -pp->z/(pp->prevz - pp->z);
	xintersect = pp->x + lambda*(pp->prevx - pp->x); 
		//x-coordinate of where photon path intersects the z=0 plane.
	yintersect = pp->y + lambda*(pp->prevy - pp->y);

	mytag = pp->tag;

	if ((fp = fopen(fname_zemaxout, "a")) == NULL) {
		perror("Cannot open FL_origin.%d for writing, quitting. \n");
		exit(EXIT_FAILURE);
	}
	
	if (mytag==0) {
		wl = 0.36;
	}
	else if (mytag==1){
		wl = 0.42;
	}
	else if (mytag==2){
	  wl = 0.55;}
	else if (mytag==3){
	  wl = 0.45;}
	else {
	  wl = 0.53;}
	fprintf(fp, "%1.3g\t %1.3g\t %1.3g\t %1.3g\t %1.3g\t %1.3g\t %1.3g\t %1.3g\n", xintersect, yintersect, 0.0f, pp->oldcosalpha, pp->oldcosbeta, pp->oldcosgamma,pp->currentweight,wl);
	fclose(fp);
}
