// trmc_2fls_tc.c : time resolved monte carlo code capable of simulating 
// multiple fluorophores/layer -- the 'tissue creator' 
// Parses command line for the input model-tissue-file and creates 
// the tissues that will be used for by the simulation. 


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h> 

#include "trmc_2fls_base.h"
#include "trmc_2fls_opts.h"

//***************************************************/
tissue *parseFile(int argc, char **argv)
{
	char fname[300];	//input file name
	FILE *fp;
	int inpfiles, i, j, k;
	int linecount = 0,  skipped = 0;
	int rcfiletissuecount;
	char toparse[1024] = "";
	tissue *list_of_ts = NULL;
	tissue *c_tp, *tp = NULL;

	if (argc < 2)
	{
		//showUsage();
		fprintf(stdout, "Need at least one input file to simulate! Aborting ...\n");
		exit(EXIT_FAILURE);
	}
	
	for (inpfiles=1; inpfiles<argc; inpfiles++) {
		strcpy(fname, argv[inpfiles]);
		if ((fp = fopen(fname,"r")) == NULL) {
			fprintf(stdout, "Cannot open file %s - fatal error, aborting...  \n", fname);
			exit(EXIT_FAILURE);
		}
	linecount=0;
	rcfiletissuecount = 0;
		while (!feof(fp)) {
			linecount += (skipped = readLines(fp, toparse));
			if (!skipped) break; //then try and look for the next file

			MALLOC(tp, sizeof(tissue));
			rcfiletissuecount++; 
			tp->next_tissue = NULL;
			if (list_of_ts == NULL)
				c_tp = list_of_ts = tp;
			else {
				c_tp->next_tissue = tp; //save prev tissue
				c_tp = tp;
			}

			// get the #of photons, 
			if (getNumberPhotons(toparse, tp) < 0) {
				fprintf(stdout, "failed getting number of photons in \"%s\": %s: line %d", toparse, fname, linecount);
				exit(EXIT_FAILURE);
			}
			//get the refractive indicies
			linecount += (skipped = readLines(fp, toparse));
			if (getTopBottRefIndex(toparse, tp) < 0 ) {
				fprintf(stdout, "failed getting two refractive indicies from \"%s\": %s: line %d \nAborting. \n", toparse, fname, linecount);
				exit(EXIT_FAILURE);
			}
			// get the number of wavelengths 
			linecount += (skipped = readLines(fp, toparse));
			if (getn_lambdas(toparse, tp) < 0 ) {
				fprintf(stdout, "failed getting #of lambdas from \"%s\": %s: line %d \nAborting. \n", toparse, fname, linecount);
				exit(EXIT_FAILURE);
			}
			// allocate the lambdas array and fill it... 
			MALLOC(tp->lambdas, sizeof(double)*tp->n_lambda);
			linecount += (skipped = readLines(fp, toparse));
			if (getlambdaVals(toparse, tp) < 0 ) {
				fprintf(stdout, "failed getting lambdas from \"%s\": %s: line %d \nAborting. \n", toparse, fname, linecount);
				exit(EXIT_FAILURE);
			}
			// get the number of layers for this tissue 
			linecount += (skipped = readLines(fp, toparse));
			if (getnumLayers(toparse, tp) < 0 ) {
				fprintf(stdout, "failed getting # of layers from \"%s\": %s: line %d \nAborting. \n", toparse, fname, linecount);
				exit(EXIT_FAILURE);
			}
			// now its time to create the layers and all the gory, necessary, 
			// 2d arrays of optical/fluorophore coefficients ...
			MALLOC(tp->layerptr, sizeof(layer)*tp->numlayers);

			for (i=0; i<tp->numlayers; i++){
				MALLOC(tp->layerptr[i].opt_trans_coeffs, sizeof(double *)*tp->n_lambda);
				MALLOC(tp->layerptr[i].muafs, sizeof(double *)*tp->n_lambda);
				MALLOC(tp->layerptr[i].taus, sizeof(double *)*tp->n_lambda);
				MALLOC(tp->layerptr[i].flqys, sizeof(double *)*tp->n_lambda);
				tp->layerptr[i].n_lambda = tp->n_lambda;
			//one loop for the 2-d fluorophore arrays... 
				for(j=0; j<tp->n_lambda; j++) {
					MALLOC(tp->layerptr[i].opt_trans_coeffs[j], sizeof(double)*4);
					memset(tp->layerptr[i].opt_trans_coeffs[j], 0, sizeof(double)*4);
					MALLOC(tp->layerptr[i].muafs[j], sizeof(double)*tp->n_lambda);
					memset(tp->layerptr[i].muafs[j], 0, sizeof(double)*tp->n_lambda);
					MALLOC(tp->layerptr[i].taus[j], sizeof(double)*tp->n_lambda);
					memset(tp->layerptr[i].taus[j], 0, sizeof(double)*tp->n_lambda);
					MALLOC(tp->layerptr[i].flqys[j], sizeof(double)*tp->n_lambda);
					memset(tp->layerptr[i].flqys[j], 0, sizeof(double)*tp->n_lambda);
				} //for(j=0; j<tp->n_lambda; j++) 
			} //for (i=0; i<tp->numlayers; i++)

			//now the tissue is set up - have to fill in the values... 
			for (i=0; i<tp->numlayers; i++){
			// first the optical transport coeffs... 
				for(j=0; j<tp->n_lambda; j++) {
					linecount += (skipped = readLines(fp, toparse));
				//	fprintf(stdout,"%d: to parse: '%s', before getopt_Coeffs (n_lambda = %d) \n", j,toparse, tp->n_lambda);
					if (getopt_Coeffs(toparse, tp->layerptr[i].opt_trans_coeffs[j]) < 0 ) {
						fprintf(stdout, "failed getting opt-coeffs of layer %d, for lambda(%d) = %#1.3g \n", i, j, tp->lambdas[j]);
						fprintf(stdout, "from \"%s\": %s: line %d \nAborting. \n", toparse, fname, linecount);
						exit(EXIT_FAILURE);
					}
					//fprintf(stdout,"%d: to parse: '%s', after getopt_Coeffs (n_lambda = %d) \n", j,toparse, tp->n_lambda);
				} //for(j=0; j<tp->n_lambda; j++) 
			//next to the fluorophore coeffs... 
				for(j=1; j<tp->n_lambda; j++) 
				{ //start from index 1, zeroth index ignored since it corresponds to 
					// excitation lambda...
					k = j; // make fluorophore matrix diagnol in the top-half... 
					do {
					linecount += (skipped = readLines(fp, toparse));
						if (getfluor_Coeffs(toparse, &(tp->layerptr[i]), j, k) < 0 ) {
							fprintf(stdout, "failed getting fluorophore-coeffs of layer %d, for lambda(%d) = %#1.3g \n", i, j, tp->lambdas[j]);
							fprintf(stdout, "from \"%s\": %s: line %d \nAborting. \n", toparse, fname, linecount);
							exit(EXIT_FAILURE);
						}
					} while (++k < tp->n_lambda);
				} //for(j=0; j<tp->n_lambda; j++) 

			// lastly, the layer thicknesses... 
				linecount += (skipped = readLines(fp, toparse));
					if (getLayerThickness(toparse, tp, i) < 0 ){
						fprintf(stdout, "failed getting thicness of layer %d, for lambda(%d) = %#1.3g \n", i, j, tp->lambdas[j]);
						fprintf(stdout, "from \"%s\": %s: line %d \nAborting. \n", toparse, fname, linecount);
						exit(EXIT_FAILURE);
					}
					else {
						fprintf(MCDEBUG, "For layer %d: ", i);
						fprintf(MCDEBUG, " thickness = %#2.3g, z_start = %#2.3g, z_end = %#2.3g \n", tp->layerptr[i].thickness,  tp->layerptr[i].z_start, tp->layerptr[i].z_end);
					}

			} //for (i=0; i<tp->numlayers; i++)
		} //while (!feof(fp))
		fprintf(stdout, "Read a total of %d tissues, from %s \n", rcfiletissuecount, fname);
		fclose(fp);
	} // for (inpfiles=1; inpfiles<argc; inpfiles++) 
	
	return(list_of_ts);
} // tissue *parseFile(int argc, char **argv)
//***************************************************/
int readLines(FILE *fp, char *input)
{
	// reads fp for valid line of input. 
	// reads a file skipping all blank lines and lines starting with # 
	// returns the string that contains the line to be parsed... 

	int iscomment = 1;
	int linesread = 0;
	char *cp;
	static char rcinput[1024] = "";
	
	strcpy(input, "");
	strcpy(rcinput, "");
	while(iscomment && !feof(fp)) {
		if ((fgets(input, 1024, fp) == NULL) && !feof(fp)) {
			fprintf(MCDEBUG, "readLines(): can't read file, aborting. \n");
			exit(EXIT_FAILURE);
		}
		if (feof(fp) ) break; 
		linesread++;
		// find the \n and terminate str
		if ((cp=strchr(input, '\n')) != NULL) 
			*cp = '\0';	 
		if ((cp=strchr(input, '#')) != NULL) 
			*cp = '\0';	 //set the location of # to be the string terminator
		sscanf(input, "%s", rcinput); //get rid of initial whitespaces...
		if (strspn(rcinput, "0123456789zerosZEROS")) 
			iscomment = 0;
		
		//fprintf(MCDEBUG, "rcinput: %s \ninput = %s\n", rcinput, input);
	} //while(iscomment && !feof(fp)) 
		if (feof(fp))
			return(0);
		else
			return(linesread); //this is zero IFF feof was hit in this fgets...
} // int readLines(FILE *fp, char *input)
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/
int getNumberPhotons(const char *str, tissue *tp)
{
	// convert str to a double and return it. 
	// return a -ve number if theres a problem, 
	// this will work because there is no possible -ve input for
	// any of the tissue parameters. 
	
	double d;
	int p;

	p = sscanf(str, "%lg", &d);
	if (!p || p == EOF)
		return(-100);
	
	tp->numphotons = (long)d;

	fprintf(MCDEBUG, "getNumberPhotons() returning %ld \n", tp->numphotons);
	return(p);


}
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/
int getTopBottRefIndex(const char *str, tissue *tp)
{
	//gets the ref-indicies of the top/bottom layers
	// returns the total number of indicies read on success,
	// and sets the tp-> values, returns a -ve value on failure
	// and tp-> are undefined. 

	char l_str[1024], *cp;
	char delims[] = DELIMITERS;
	double d[2];
	int p, vals=0;


	strcpy(l_str, str);
	cp = strtok(l_str, delims);

	while(cp != NULL && vals <2) {
		p = sscanf(cp, "%lg", &d[vals]);
		if (p == EOF) return(-100);
		if (p) vals++;
		cp = strtok(NULL, delims);
	}

	if (cp) fprintf(MCDEBUG, "trailing characters in %s, ignoring (could be incorrect input?) \n", str);
	
	if (vals == 2) {
		tp->n_top = d[0];
		tp->n_bottom = d[1];

		fprintf(MCDEBUG, "getTopBottRefIndex() returning %#2.3g, %#2.3g \n", tp->n_top, tp->n_bottom);
		return(vals);
	}
	else
		return(-100);

}
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/
int getn_lambdas(const char *str, tissue *tp)
{
	// convert str to a double and return it. 
	// return a -ve number if theres a problem, 
	// this will work because there is no possible -ve input for
	// any of the tissue parameters. 
	
	double d;
	int p;

	p = sscanf(str, "%lg", &d);
	if (!p || p == EOF)
		return(-100);
	
	tp->n_lambda = (int)d;

	fprintf(MCDEBUG, "getn_lambdas returning %d \n", tp->n_lambda);
	return(p);
} //int getn_lambdas(const char *str, tissue *tp)
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/
int getlambdaVals(const char *str, tissue *tp)
{ 	//gets the list of lambdas as specified by tp->n_lambda
	// returns the total number read on success
	// returns a -ve number on failure and tp->* will be undefined. 

	char l_str[1024], *cp;
	char delims[] = DELIMITERS;
	int p, vals=0;


	strcpy(l_str, str);
	cp = strtok(l_str, delims);

	while(cp != NULL && vals < tp->n_lambda) {
		p = sscanf(cp, "%lg", tp->lambdas+vals);
		if (p == EOF) return(-100);
		if (p) vals++;
		cp = strtok(NULL, delims);
	}

	if (cp) fprintf(MCDEBUG, "trailing characters in %s, ignoring (could be incorrect input?) \n", str);
	
	if (vals == tp->n_lambda) {
		fprintf(MCDEBUG, "getlambdaVals() returning lambda list: [");
		for (vals=0; vals<tp->n_lambda; vals++)
			fprintf(MCDEBUG, " %#0.4g;", tp->lambdas[vals]);
		fprintf(MCDEBUG, "]\n");
		return(vals);
	}
	else
		return(-100);

} //int getlambdaVals(const char *str, tissue *tp)
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/
int getnumLayers(const char *str, tissue *tp)
{
	
	// convert str to int return it. 
	// return a -ve number if theres a problem, 
	// this will work because there is no possible -ve input for
	// any of the tissue parameters. 
	
	double d;
	int p;

	p = sscanf(str, "%lg", &d);
	if (!p || p == EOF)
		return(-100);
	
	tp->numlayers = (int)d;

	fprintf(MCDEBUG, "getnumLayers returning %d \n", tp->numlayers);
	return(p);
} // int getnumLayers(const char *str, tissue *tp)
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/
int getopt_Coeffs(const char *str, double *dp)
{
	//get 4 coeffs from the line into the dp array.
	// caller ensures dp is the right one.. etc.

	char l_str[1024], *cp;
	char delims[] = DELIMITERS;
	int p, vals=0;


	strcpy(l_str, str);
	cp = strtok(l_str, delims);

	while(cp != NULL && vals < 4) { //need 4 coeffs...
		p = sscanf(cp, "%lg", dp+vals);
		if (p == EOF) return(-100);
		if (p) vals++;
		cp = strtok(NULL, delims);
	}

	//if (cp) fprintf(MCDEBUG, "trailing characters in %s, ignoring (could be incorrect input?) \n", str);
	
	if (vals == 4) {
		return(vals);
	}
	else
		return(-100);
	
} //int getopt_Coeffs(const char *str, double *dp)
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/
int getfluor_Coeffs(const char *str, layer *lp, int r, int c)
{
	// each line in the input file specifying flurophore coeffs. 
	// has entries for the muafs, taus and flqys. 
	// Hence we need the row, r and col c, to set them right. 
	// caller takes care of looping the r/c indicies... 

	char l_str[1024], *cp;
	char delims[] = DELIMITERS;
	int p, vals=0;
	double d[3];

	strcpy(l_str, str);
	cp = strtok(l_str, delims);

	if (strstr(cp, "zeros") || strstr(cp, "ZEROS") )
	{
		//found the word zeros/ZEROS in the input string. 
		d[0] = d[1] = d[2] = 0.0;
		vals = 3; //we have all the vals... 
	}
	while(cp != NULL && vals < 3) { //need 3 coeffs...
		p = sscanf(cp, "%lg", d+vals);
		if (p == EOF) return(-100);
		if (p) vals++;
		cp = strtok(NULL, delims);
	}

	//if (cp) fprintf(MCDEBUG, "trailing characters in %s, ignoring (could be incorrect input?) \n", str);
	
	if (vals == 3) {
		lp->muafs[r][c] = d[0];
		lp->taus[r][c] = d[1];
		lp->flqys[r][c] = d[2];
		return(vals);
	}
	else
		return(-100);

} //int getfluor_Coeffs(const char *str, layer *lp, int r, int c)
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/
int getLayerThickness(const char *str, tissue *tp, int l_num)
{
	// get the value of thickness, and set z-start/z-end as well
	// l_num controls the calculation of z_start, z_end

	double d;
	int p;

	p = sscanf(str, "%lg", &d);
	if (!p || p == EOF)
		return(-100);
	
	tp->layerptr[l_num].thickness = d;

	if (!l_num)	 //i.e l_num == 0 => this is the first layer... 
		tp->layerptr[l_num].z_start= 0.0;
	else
		tp->layerptr[l_num].z_start = tp->layerptr[l_num-1].z_end;

	tp->layerptr[l_num].z_end = tp->layerptr[l_num].z_start + d;

	return(p);
	
} //int getLayerThickness(const char *str, layer *lp)
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/
void fprintTissue(FILE *fp, const tissue *tp)
{
	// print all info for the given tissue to fp 
	
	int i,j;

	//fprintf(fp, "\n----------New tissue def start:----------\n");
	fprintf(fp, "Number of photons: %#2.2E \n", (double)tp->numphotons);
	fprintf(fp, "n_top = %#0.3f; n_bottom = %0.3f \n", tp->n_top, tp->n_bottom);
	fprintf(fp, "supported wavelengths, n_lambda = %d \n\t[", tp->n_lambda);
	for (i=0; i<tp->n_lambda; i++)
		fprintf(fp, "%#3.1f; ", tp->lambdas[i]);
	fprintf(fp, "]\n");
	fprintf(fp, "tissue has %d layers: \n", tp->numlayers);
	for (i=0; i<tp->numlayers;i++) {
		fprintf(fp, "Layer %d: \n", i);
		fprintLayer(fp, &(tp->layerptr[i]), tp->lambdas);
		for (j=0; j<tp->n_lambda; j++)
			fprintf(fp, "lambda(%d): albedo = %#2.3f; optical depth = %#2.3f \n", j, (tp->layerptr[i]).opt_trans_coeffs[j][MUS]/((tp->layerptr[i]).opt_trans_coeffs[j][MUA] + (tp->layerptr[i]).opt_trans_coeffs[j][MUS]), (tp->layerptr[i]).thickness*((tp->layerptr[i]).opt_trans_coeffs[j][MUA] + (tp->layerptr[i]).opt_trans_coeffs[j][MUS]));
	}
	//fprintf(fp, "----------New tissue def end:----------\n");
} //void fprintTissue(FILE *fp, const tissue *tp)
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/
void fprintLayer(FILE *fp, const layer *lp, const double *lambdas)
{
	// print all info for the layer, lp 

	int i,j;
	
	fprintf(fp, "  transport coeffs: (mua/mus/g/ref_ind)\n");
	for (i=0; i<lp->n_lambda; i++)
		fprintf(fp, "    at %#3.1fnm [lambda(%d)]: %#1.3E/ %#1.3E/ %#1.3E/ %#1.3E \n", lambdas[i], i, lp->opt_trans_coeffs[i][0], lp->opt_trans_coeffs[i][1], lp->opt_trans_coeffs[i][2], lp->opt_trans_coeffs[i][3]);
	
	if (lp->n_lambda > 1) {
		fprintf(fp, "  Total muaf's stats: \n");
		for (i=1; i<lp->n_lambda; i++)
			fprintf(fp, "    abs of %#3.1fnm: %#2.3f 1/cm\n", lambdas[i-1], getTotalmuaf(lp, i));
		fprintf(fp, "  flurophore-coeffs: (muaf/tau/flqy)\n");
		for (i=1; i<lp->n_lambda; i++){
			for (j=i; j<lp->n_lambda; j++){
				fprintf(fp, "    abs->em  %#1.3E nm->%#1.3E nm: %#1.3E/ %#3.2E/ %#1.3E \n", lambdas[i-1], lambdas[j], lp->muafs[i][j], lp->taus[i][j], lp->flqys[i][j]);
			} //for (j=1; j<lp->n_lambda; j++)
		} //for (i=1; i<lp->n_lambda; i++)
	} //if (lp->n_lambda > 1) 
	
	fprintf(fp, "  z_start = %#3.3f cm; z_end = %#3.3f cm \n", lp->z_start, lp->z_end);
		
} //void fprintLayer(FILE *fp, const layer *lp, const double *lambdas)
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/
double getTotalmuaf(const layer *lp, int tag)
{
	//sums up the muafs[tag] array of lp and returns it. 
	// the matrix is gaurenteed to have zeros for non-available data.. 
	// tag = 1, for absorption of excitation, 
	// tag MUST be < n_lambda
	int i;
	double d=0.0; 

	for (i=1; i<lp->n_lambda; i++)
		d += lp->muafs[tag][i];
	
	return(d);

}
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/
double getLayeropt_coeff(const layer *lp, int r, int c)
{
	//return the optical coefficients from the opt_trans_coeffs member 
	// of lp. 

	if (r > lp->n_lambda-1 || c >3) {
		fprintf(stderr, "incorrect arguments passed to getLayeropt_coeff() \n");
		fprintf(stderr, "n_lambda = %d, r = %d, c = %d; aborting \n", lp->n_lambda, r, c);
		exit(EXIT_FAILURE);
	}

	return(lp->opt_trans_coeffs[r][c]);
} //double getLayeropt_coeff(const layer *lp, int r, int c)
//***************************************************/
