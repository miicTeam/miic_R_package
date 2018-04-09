#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <time.h>
#include <string>
#include <map>
#include <iostream>
#include <sstream>
#include <tuple>
#include <algorithm>
#include <unistd.h>
#include <Rcpp.h>


#include "utilities.h"
#include "modules_MI.h"
#include "computeInfo.h"


#define STEPMAX 20
#define EPS 0.01
#define INIT_EQUAL_WIDTH true
#define ALPHA_EFF_LVLS 1

using namespace Rcpp;
using namespace std;

// module to compute the conditional mutual information for mixed continuous and discrete variables

//////////////////////////////////////////////////////////////////////////////////////////////


//#define _MY_DEBUG_MInoU 1
//#define _MY_DEBUG_MI 1
//#define _MY_DEBUG_NEW_OPTFUN 1


///////////////////////////////////////////////////////////////////////////////////////////

//GENERAL VARIABLES NOTATION

//int n: number of sample
//int coarse : coarse graining step : minimum length of a bin : considering only n/coarse possible cuts

//double* looklog : lookup table for the logarithms of natural numbers up to n
//map<string,double> &looklbc : lookup table for the MDL complexity terms
//double* c2terms : precomputed vectors for the MDL complexity

///////////////////////////////////////////////////////////////////////////////////////////

// DYNAMIC PROGRAMMING OPTIMIZATION FUNCTION

//complicated (to handle) but flexible function
//optimize linear function depending on different variables

//////////////////////////////////////////////
//INPUT

//optimizing variable defined by <sortidx_var>

//<nbrV> number of terms in the optimization functions, consisting in variables and joint variables

//<factors> : ( <nbrV> x n )  matrix :
//all factors of the <nbrV> terms
//example <nbrV>=6
//x factors : factors[0]
//y factors : factors[1]
//u factors : factors[2]
//xu factors : factors[3]
//yu factors : factors[4]
//xyu factors : factors[5]

//<r> : (<nbrV> x 1) vector : number of levels of the <nbrV> terms

//<optfun> : ( <nbrV> + 1 x 1) vector : coefficient of the different terms in the optimization function f
//optimization function f=sum_k optfun[k]*H[k]
//f= optfun[nbrV+1]*H_optvar +optfun[0]*Hx + optfun[1]*Hy + optfun[2]*Hu + optfun[3]*Hxu + optfun[4]*Hyu + optfun[5]*Hxyu
//note optfun[nbrV+1] is reserved for optimized variable only

//USE EXAMPLE
//
//(1) I(x;yu)=Hx+Hyu-Hxyu :
//			optimize f on x : nbrV=2, terms to take into account Hx-Hxyu :
//				optfun[0]=0 //factors x= 0
//				optfun[6]=1 => optimization on x
//				optfun[4]= -1 => factors[4]	// factors xyu => factors yu
//				opfun=[0,0,0,0,-1,1]
//			optimize f on y : 2 terms to take into account Hyu-Hxyu :
//				optfun[0]=0 //factors x= 0
//				optfun[6]=0 => no terms on y
//				optfun[4]=1 => factors[4] //factors yu=> factors u
//				optfun[5]=-1 => factors[5] //factors xyu=> factors xu
//				opfun=[0,0,0,1,-1,0]
//
//(2) I(xu;y)=Hy+Hxu-Hxyu : opfun=[0,1,0,1,0,-1]
//
//(3) I(x;u)=Hx+Hu-Hxu : opfun=[1,0,1,-1,0,0]

//<sc> : terms in the MDL compelxity, that depends on the definition of the optimization function:
// example for I(x,y) optimized on x => sc=0.5*(r_y-1)

//OUTPUT
//int *memory_cut: output vector of optimal cuts :
//
inline __attribute__((always_inline))
double optfun_onerun_kmdl_coarse(int *sortidx_var, int *data, int nbrV, int **factors, int *r, double pxy,
								 double sc, double sc_levels1, double sc_levels2, int n, int nnr,int *cut,
								 int *r_opt, int maxbins, double* looklog, double** looklbc, double* lookH,
								 double** sc_look, int cplx) {

		// cout << "infunctionlevels: " << nnr << endl;

		int i,j,k,m;

		sc_levels1 = round(sc_levels1);
		sc_levels2 = round(sc_levels2);

		int coarse=ceil(1.0*nnr/maxbins);//step coarse graining
		if (coarse<1) coarse=1;
		//coarse= cbrt(n);
		int np=ceil(1.0*nnr/coarse); //number of possible cuts

		//temp variables to memorize optimization cuts
		int *memory_cuts_idx=(int *)calloc(np,sizeof(int));//indexes of the cuts (1..np)
		int *memory_cuts_pos=(int *)calloc(np,sizeof(int));//positions of the cuts (1..n)

		//variables for the computation of the complexity
		double sc2; // The cost for adding one bin.


		//dynamic programming optimize function and memorize of cuts
		double Imax;//I-kmdl
		int* nc =(int *)calloc(np,sizeof(int));
		int nctemp;//Hx-Hxy-LogCx

		//entropy in kj interval for the <nbrV>+1 terms
		//double* H_kj=(double *)calloc(nbrV+1,sizeof(double));
		double* H_kj=(double *)calloc(nbrV,sizeof(double)); //herve

		double* weights=(double *)calloc(4,sizeof(double)); //herve

		//function max value at each step
		double* I=(double *)calloc(np,sizeof(double)); // The optimal information value found at each idx.
		double* I_0k =(double *)calloc(np,sizeof(double)); // The information value for a unique bin fom idx 0 to k.
		double I_kj; // The information value for a bin fom idx k to j.
		double t;

		//number of points  in the intervals
		int nxj,nx;

		int xyu;

		int ir; // Iterator on non repeated values
		int nk,nj,nkj; // Number of values in intervals [0,k], [0,j], [k,j]
		int njforward,nkforward; // Indexes at current position of j and k


		//dynamic programming steps terms in the 0k interval

		//H_0k needs initialization at zero (problem: double H_0k[nbrV+1][np];)
		//double **H_0k=(double **)calloc(nbrV+1,sizeof(double*));// x y u xu yu xyu
		//for(m=0;m<nbrV+1;m++){
		double **H_0k=(double **)calloc(nbrV,sizeof(double*));// x y u xu yu xyu //herve
		for(m=0;m<nbrV;m++){
			H_0k[m]=(double *)calloc(np,sizeof(double));
		}

		int **nxyu=(int**)calloc(nbrV,sizeof(int*));// x y u xu yu xyu
		for(m=0;m<nbrV;m++){
			nxyu[m]=(int*)calloc(r[m],sizeof(int));
		}

		int **nxyu_k=(int**)calloc(nbrV,sizeof(int*));// x y u xu yu xyu
		for(m=0;m<nbrV;m++){
		 	nxyu_k[m]=(int*)calloc(r[m],sizeof(int));
		 	//cout <<  "m: " << m << "\tniv: " << r[m] << "\n";
		}

		//boolean vector, check_repet[i] is true if data[i]!=data[i+1]
		bool *check_repet = (bool*)calloc(n, sizeof(bool));
		for(i=0;i<(n-1);i++){
			check_repet[i] = (data[sortidx_var[i+1]]!=data[sortidx_var[i]]);
		}

		//matrix with number of observations of each level (of the non-opt variable) at idx j
		int **rx_0j=(int **)calloc(r[0],sizeof(int*));
		for(k=0;k<r[0];k++){
			rx_0j[k]=(int *)calloc(n,sizeof(int));
		}

		//matrix with number of unique levels of the non-opt variable from idx k to j on the opt variable.
		int **rx_kj_diff=(int **)calloc(n,sizeof(int*));
		for(k=0;k<n;k++){
			rx_kj_diff[k]=(int *)calloc(n,sizeof(int));
		}
		///////////////////////////////////////////////
		#if _MY_DEBUG_NEW_OPTFUN
			printf("\n-----------\n");
			printf("=> coarse=%d np=%d nbrV=%d\n",coarse,np,nbrV);

			for(m=0;m<nbrV;m++){
				printf("r[%d]=%d :\n",m,r[m]);

				for(i=0;i<n;i++){
					printf("%d ",factors[m][i]);
				}
				printf("\n");

			}
			printf("factors[0][sortidx_var[i]] :\n");
			for(i=0;i<n;i++){
				printf("%d ",factors[0][sortidx_var[i]]);
			}
			printf("\n");

		#endif

		/////////////////////////////////////////////////
		//j=0;
		#if _MY_DEBUG_NEW_OPTFUN
			printf("j=%d\n",0);fflush(stdout);
		#endif

		//////////// WEIGHTS ///////////////////////////////////// //herve
			if(nbrV==2){
			  weights[0]=-pxy;
			  weights[1]=pxy;
 			}else if(nbrV==4){
			  weights[0]=-2*pxy;
			  weights[1]=2*pxy-1;
			  weights[2]=2*pxy-1;
			  weights[3]=2*(1-pxy);
			}

		if(cplx==1){ //NML
		// Filling the rx_0j matrix
			rx_0j[factors[0][sortidx_var[0]]][0] = 1;
			for(j=1; j<n; j++){
				for(k=0; k<r[0]; k++){
					rx_0j[k][j] = rx_0j[k][j-1];
				}
				rx_0j[factors[0][sortidx_var[j]]][j] ++;
			}
			// Filling the rx_kj difference in levels from k to j matrix
			for(k=0;k<n;k++){
					for(j=k;j<n-1;j++){
					int nlevels = 0;
					for(int lvl=0;lvl<r[0];lvl++){
						nlevels += (rx_0j[lvl][j] - rx_0j[lvl][k]) > 0;
					}
					rx_kj_diff[k][j] = nlevels;
				}
			}
		}
		//computing statistics of the <nbrV> terms and the entropy

		njforward=0;//iterator on values
		ir=0;//iterator on not repeated values
		while(ir<coarse){

			for(m=0;m<nbrV;m++){

			  xyu=factors[m][sortidx_var[njforward]];

			  nxyu[m][xyu]++;

			  if(nxyu[m][xyu] != 1){
			    //H_0k[m][0] -= nxyu[m][xyu]*looklog[nxyu[m][xyu]] - (nxyu[m][xyu]-1)*looklog[nxyu[m][xyu]-1];
			    H_0k[m][0] += weights[m]*lookH[nxyu[m][xyu]-1]; //herve
			  }
			  if(m == 1){ //herve
			    if(cplx==0 && nxyu[m][xyu]==1) H_0k[m][0]  -= pxy * sc * looklog[n];  //MDL
			    //else if(cplx==1) H_0k[m][0] -= pxy * ( computeLogC(nxyu[m][xyu], fmin(sc_levels1, ALPHA_EFF_LVLS*rx_kj_diff[0][njforward+1]), sc_look) +
							//							logchoose(sc_levels1, fmin(sc_levels1, ALPHA_EFF_LVLS*rx_kj_diff[0][njforward+1]), looklog) -
							//				   		   computeLogC(nxyu[m][xyu]-1, fmin(sc_levels1, ALPHA_EFF_LVLS*rx_kj_diff[0][njforward]), sc_look) -
							//							logchoose(sc_levels1, fmin(sc_levels1, ALPHA_EFF_LVLS*rx_kj_diff[0][njforward]), looklog) );//NML
			    else if(cplx==1) H_0k[m][0] -= pxy * ( computeLogC(nxyu[m][xyu], sc_levels1, sc_look) - computeLogC(nxyu[m][xyu]-1, sc_levels1, sc_look) );//NML
			  }else if(m == 2){ //herve
			    if(cplx==0 && nxyu[m][xyu]==1) H_0k[m][0]  -= pxy * sc * looklog[n];  //MDL
			    else if(cplx==1) H_0k[m][0] -= pxy * ( computeLogC(nxyu[m][xyu], sc_levels2, sc_look) - computeLogC(nxyu[m][xyu]-1, sc_levels2, sc_look) );//NML
			  }else if(m == 3){ //herve
			    if(cplx==0 && nxyu[m][xyu]==1) H_0k[m][0]  -= (1-pxy) * sc * looklog[n];  //MDL
			    else if(cplx==1) H_0k[m][0] -= (1-pxy) * ( computeLogC(nxyu[m][xyu], sc_levels1, sc_look) - computeLogC(nxyu[m][xyu]-1, sc_levels1, sc_look) + computeLogC(nxyu[m][xyu], sc_levels2, sc_look) - computeLogC(nxyu[m][xyu]-1, sc_levels2, sc_look) );//NML
			  }
			}

			ir += int(check_repet[njforward]);
			njforward++;

		}

		#if _MY_DEBUG_NEW_OPTFUN
			printf("\n(j=%d njforward=%d   ir=%d )\n",0,njforward,ir);
		#endif



		#if _MY_DEBUG_NEW_OPTFUN
			printf("  ");
			for(m=0;m<nbrV;m++) {
				for(xyu=0;xyu<r[m];xyu++) printf("nxyu[%d][%d]=%d  ",m,xyu,nxyu[m][xyu]);fflush(stdout);
			}
			printf("\n");
		#endif

		nj=njforward;

		//if(optfun[nbrV] != 0) H_0k[nbrV][0]=-nj*looklog[nj]; //herve

		//for(m=0;m<nbrV+1;m++) if(optfun[m] != 0)	I_0k[0]+=optfun[m]*H_0k[m][0];
		//for(m=0;m<nbrV;m++) if(optfun[m] != 0)	I_0k[0]+=optfun[m]*H_0k[m][0];//herve
		for(m=0;m<nbrV;m++) I_0k[0] += H_0k[m][0];//herve

		//initizialitation
		memory_cuts_idx[0]=0;
		I[0]=I_0k[0];
		nc[0]=1;

		//if(cplx==0){ //MDL
		//	I[0] -= sc * looklog[n];
		//}
		//else if(cplx==1){ //NML
		//	I[0] -= computeLogC(nj, sc_levels, sc_look);
		//}

		#if _MY_DEBUG_NEW_OPTFUN
				printf("  ");
				for(m=0;m<nbrV+1;m++) printf("H_0k[%d]=%lf ",m,H_0k[m][0]);fflush(stdout);
				printf("\n  [0 -0] = %lf\n",I_0k[0]);fflush(stdout);
		#endif

		//moving j over the np possible cuts

		for(j=1;j<=np-1;j++){ //j=1...n-1

			#if _MY_DEBUG_NEW_OPTFUN
					printf("j=%d\n",j);fflush(stdout);
			#endif

			//COMPUTING STATISTICS AND FUNCTION FOR THE INTEVAL [0 j] -> I_0k

			//initialization with previous interval [0 j-1]
			//for(m=0;m<nbrV+1;m++) if(optfun[m] != 0) H_0k[m][j]=H_0k[m][j-1];
			for(m=0;m<nbrV;m++) H_0k[m][j]=H_0k[m][j-1]; //herve

			//computing statistics of the <nbrV> terms and the entropy

			// njforward iterator on values
			ir=0;//iterator on not repeated values
			while((ir<coarse )&& (njforward<n)){

				for(m=0;m<nbrV;m++){

				  xyu=factors[m][sortidx_var[njforward]];

				  nxyu[m][xyu]++;

				  if(nxyu[m][xyu] != 1){
				    //H_0k[m][0] -= nxyu[m][xyu]*looklog[nxyu[m][xyu]] - (nxyu[m][xyu]-1)*looklog[nxyu[m][xyu]-1];
				    H_0k[m][j] += weights[m]*lookH[nxyu[m][xyu]-1]; //herve
				  }

				  if(m == 1){ //herve
				    if(cplx==0 && nxyu[m][xyu]==1) H_0k[m][j]  -= pxy * sc * looklog[n];  //MDL
			    	//else if(cplx==1) H_0k[m][j] -= pxy * ( computeLogC(nxyu[m][xyu], fmin(sc_levels1, ALPHA_EFF_LVLS*rx_kj_diff[0][njforward+1]), sc_look) +
								//						   	logchoose(sc_levels1, fmin(sc_levels1, ALPHA_EFF_LVLS*rx_kj_diff[0][njforward+1]), looklog) -
								//				   		   computeLogC(nxyu[m][xyu]-1, fmin(sc_levels1, ALPHA_EFF_LVLS*rx_kj_diff[0][njforward]), sc_look) -
							 // 							   	logchoose(sc_levels1, fmin(sc_levels1, ALPHA_EFF_LVLS*rx_kj_diff[0][njforward]), looklog) );//NML
				    else if(cplx==1) H_0k[m][j] -= pxy * ( computeLogC(nxyu[m][xyu], sc_levels1, sc_look) - computeLogC(nxyu[m][xyu]-1, sc_levels1, sc_look) );//NML
				  }else if(m == 2){ //herve
				    if(cplx==0 && nxyu[m][xyu]==1) H_0k[m][j]  -= pxy * sc * looklog[n];  //MDL
				    else if(cplx==1) H_0k[m][j] -= pxy * ( computeLogC(nxyu[m][xyu], sc_levels2, sc_look) - computeLogC(nxyu[m][xyu]-1, sc_levels2, sc_look) );//NML
				  }else if(m == 3){ //herve
				    if(cplx==0 && nxyu[m][xyu]==1) H_0k[m][j]  -= (1-pxy) * sc * looklog[n];  //MDL
				    else if(cplx==1) H_0k[m][j] -= (1-pxy) * ( computeLogC(nxyu[m][xyu], sc_levels1, sc_look) - computeLogC(nxyu[m][xyu]-1, sc_levels1, sc_look) + computeLogC(nxyu[m][xyu], sc_levels2, sc_look) - computeLogC(nxyu[m][xyu]-1, sc_levels2, sc_look) );//NML
				  }

				}


				if(njforward+1 < n){//check no repetition
					ir += int(check_repet[njforward]);
				}
				njforward++;

			}

			//njforward: number of points between 0 and j
			nj=njforward;

			#if _MY_DEBUG_NEW_OPTFUN
				printf("(njforward=%d   ir=%d )\n",njforward,ir);fflush(stdout);

				printf("\n   ");

				for(m=0;m<nbrV;m++) {
					for(xyu=0;xyu<r[m];xyu++) printf("nxyu[%d][%d]=%d  ",m,xyu,nxyu[m][xyu]);fflush(stdout);
				}
				printf("\n");
			#endif

			//if(optfun[nbrV] != 0) H_0k[nbrV][j]=-nj*looklog[nj]; //herve

			I_0k[j]=0;
			for(m=0;m<nbrV;m++)	I_0k[j] += H_0k[m][j]; //herve
			I[j] = I_0k[j];
			if(j==(np-1)){
				//I[j] = -DBL_MAX;
			}


			Imax=-DBL_MAX;

			for(m=0;m<nbrV;m++) H_kj[m]=H_0k[m][j]; //herve

			#if _MY_DEBUG_NEW_OPTFUN
				printf("  ");
				printf("\n  [0 -%d] = %lf\n",j,I_0k[j]);fflush(stdout);
			#endif

			for(m=0;m<nbrV;m++) {
			  for(xyu=0;xyu<r[m];xyu++)
			    nxyu_k[m][xyu]=nxyu[m][xyu];
			}

			//moving k

			//k iterator on possible cuts
			nkforward=0;//iterator on values
			//it iterator on not repeated values

			for(k=0;k<=j-1;k++){//k=1...n-2 possible cuts

						ir=0;//iterator on not repeated values
						while( ir<coarse ){

						  for(m=0;m<nbrV;m++){

						    xyu=factors[m][sortidx_var[nkforward]];
						    nxyu_k[m][xyu]--;

						    H_kj[m] -= weights[m]*lookH[nxyu_k[m][xyu]]; //herve

						    if(m == 1){ //herve
						      if(cplx==0 && nxyu_k[m][xyu]==0) H_kj[m]  += pxy * sc * looklog[n];  //MDL
			    			  //else if(cplx==1) H_kj[m] -= pxy * ( computeLogC(nxyu_k[m][xyu], fmin(sc_levels1, ALPHA_EFF_LVLS*rx_kj_diff[nkforward+1][njforward-1]), sc_look) +
							  			//							logchoose(sc_levels1, fmin(sc_levels1, ALPHA_EFF_LVLS*rx_kj_diff[nkforward+1][njforward-1]), looklog) -
            //                                    				  computeLogC(nxyu_k[m][xyu]+1, fmin(sc_levels1, ALPHA_EFF_LVLS*rx_kj_diff[nkforward][njforward-1]), sc_look) -
							  			//							logchoose(sc_levels1, fmin(sc_levels1, ALPHA_EFF_LVLS*rx_kj_diff[nkforward][njforward-1]), looklog) );//NML
						      else if(cplx==1) H_kj[m] -= pxy * ( computeLogC(nxyu_k[m][xyu], sc_levels1, sc_look) - computeLogC(nxyu_k[m][xyu]+1, sc_levels1, sc_look) );//NML
						    }else if(m == 2){ //herve
						      if(cplx==0 && nxyu_k[m][xyu]==0) H_kj[m]  += pxy * sc * looklog[n];  //MDL
						      else if(cplx==1) H_kj[m] -= pxy * ( computeLogC(nxyu_k[m][xyu], sc_levels2, sc_look) - computeLogC(nxyu_k[m][xyu]+1, sc_levels2, sc_look) );//NML
						    }else if(m == 3){ //herve
						      if(cplx==0 && nxyu_k[m][xyu]==0) H_kj[m]  += (1-pxy) * sc * looklog[n];  //MDL
						      else if(cplx==1) H_kj[m] -= (1-pxy) * ( computeLogC(nxyu_k[m][xyu], sc_levels1, sc_look) - computeLogC(nxyu_k[m][xyu]+1, sc_levels1, sc_look) + computeLogC(nxyu_k[m][xyu], sc_levels2, sc_look) - computeLogC(nxyu_k[m][xyu]+1, sc_levels2, sc_look) );//NML
						    }


						  }

						  ir += int(check_repet[nkforward]);
						  nkforward++;
						}

						//nx=nxj-(k+1)*coarse;
						nkj=njforward-nkforward;
						//if(optfun[nbrV] != 0) H_kj[nbrV]=-nkj*looklog[nkj]; //herve


						#if _MY_DEBUG_NEW_OPTFUN
							printf("\nnxj=%d  \n   ",nkj);
							for(m=0;m<nbrV;m++) {
								for(xyu=0;xyu<r[m];xyu++) printf("nxyu_k[%d][%d]=%d ",m,xyu,nxyu_k[m][xyu]);fflush(stdout);
							}
							printf("\n");
						#endif

						I_kj=0;
						//for(m=0;m<nbrV+1;m++) if(optfun[m] != 0)	I_kj+=optfun[m]*H_kj[m];
						for(m=0;m<nbrV;m++)	I_kj += H_kj[m]; //herve


						#if _MY_DEBUG_NEW_OPTFUN


							printf("  ");
							for(m=0;m<nbrV+1;m++) printf("H_kj[%d]=%lf ",m,H_kj[m]);fflush(stdout);
							printf("\n   [%d -%d][%d - %d] = %lf (%lf+%lf-%lf)\n",0,k,k+1,j,I_0k[k]+I_kj-sc2,I_0k[k],I_kj,sc2);

						#endif

						nk=nkforward-1;//position of the actual possible cut

						//if( I[k] + I_kj - sc2 > I[j] ){
						//	t=I[k] + I_kj - sc2; //[0.. cuts.. k-1][k j]
						if( I[k] + I_kj > I[j] ){ //herve
							t=I[k] + I_kj; //[0.. cuts.. k-1][k j] //herve
							if (Imax<t){
								Imax=t;
								nctemp=nc[k]+1;
								I[j]=Imax; // optimized function for the interval [0 j]
								nc[j]=nctemp; // number of optimal cuts
								memory_cuts_idx[j] = k + 1; // index  of the (last) optimal cut
								memory_cuts_pos[j] = nk; // position  of the (last) optimal cut
							}
						}

						#if _MY_DEBUG_NEW
							printf("   f=%lf\n",Imax);fflush(stdout);
						#endif
			}


			#if _MY_DEBUG_NEW_OPTFUN
				printf("\n>>>j=%d: Imax=%lf cut[%d]=%d ncut=%d\n",j,Imax,j,memory_cuts_idx[j],nc[j]);fflush(stdout);

				printf("\n-----------\n");

			#endif
		}

	// free memory
	double return_res = I[np-1]/n;
	free(nc);
	free(I);
	free(I_0k);
	free(H_kj);
	free(weights);
	free(check_repet);

	if(cplx==1){
	for(k=0;k<r[0];k++) free(rx_0j[k]);
	free(rx_0j);
	for(i=0;i<n;i++) free(rx_kj_diff[i]);
	free(rx_kj_diff);
	}

	for(m=0;m<nbrV;m++) free(H_0k[m]);
	free(H_0k);

	for(m=0;m<nbrV;m++) free(nxyu_k[m]);
	free(nxyu_k);

	for(m=0;m<nbrV;m++) free(nxyu[m]);
	free(nxyu);

	// reconstruction of the optimal cuts from the memory cuts indexes and positions
	*r_opt=reconstruction_cut_coarse(memory_cuts_idx,memory_cuts_pos, np, n, cut);

	free(memory_cuts_idx);
	free(memory_cuts_pos);
	return return_res;

}



////////////////////////////////////////////////////////////////////////////////








//////////////////////////////////////////////////////////////////////////////////////////////
//compute I(x,y)

//dynamic programming for optimizing variables binning

//optimize on x I(x,y): Hx - Hxy - kmdl
//optimize on y I(x,y): Hy - Hxy - kmdl
//until convergence


int** compute_Ixy_alg1(int** data, int** sortidx, int* ptr_cnt, int* ptrVarIdx,  int* AllLevels,
						 int n, int maxbins, int **cut, int *r, double* c2terms, int init_bin,
						 double* looklog, double** looklbc, double* lookH, double** sc_look, int cplx)
{

	int j,l,ll;

	/////////////////////////////////////

	double* res=(double *)calloc(2,sizeof(double));//results res[0]->I,res[1]->I-k
	int** iterative_cuts = (int **)calloc(STEPMAX+1, sizeof(int*));
	for(int i=0; i<STEPMAX+1; i++){
		iterative_cuts[i] = (int *)calloc(maxbins*2, sizeof(int));
	}

	//////////////////////////////////
	//allocation factors  x y

	int **datafactors;
	datafactors=(int **)calloc((2),sizeof(int*));

	for(l=0;l<(2);l++){
		datafactors[l]=(int *)calloc(n,sizeof(int));
	}

	int* singlefactor = (int *)calloc(n, sizeof(int));

	int* r_temp=(int *)calloc(3,sizeof(int));

	int *ptr=(int *)calloc(2,sizeof(int));

	int *xy_factors = (int *)calloc(n,sizeof(int));

	int rxy=0;
	double Hxy;

	ptr[0]=0;
	ptr[1]=1;

	////////////////////////////////////

	//initialization of datafactors && sortidx


	for(l=0;l<2;l++){

		if(ptr_cnt[ptrVarIdx[l]]==1){
			update_datafactors(sortidx, ptrVarIdx[l], datafactors, l, n, cut);
		}
		else{
			for(j=0;j<=n-1;j++){
				datafactors[l][j]=data[ptrVarIdx[l]][j];
			}
		}
	}

	int np;
	// int np=ceil(1.0*n/coarse);//max possible cuts

	double sc;
	double sc_comb=0;//term complexity due to optimization on many models

 	double MInew;

	double* MI = (double *)calloc(STEPMAX,sizeof(double));
	double* MIk = (double *)calloc(STEPMAX,sizeof(double));
	double I_av,Ik_av;

	int stop,i,flag;

	double pxy; //herve
	int sc_levels1,sc_levels2; //herve
	//

	int **factors1=(int **)calloc(2,sizeof(int*));
	for(l=0;l<2;l++){
		factors1[l]=(int *)calloc(n,sizeof(int));
	}
	int *rt1=(int *)calloc(2,sizeof(int));

	//

		MI[0]=0;
		MIk[0]=-DBL_MAX;

		flag=0;

		// #if _MY_DEBUG_MInoU

		// 	jointfactors_u(datafactors,ptr,n, 2, r, xy_factors, &rxy,&Hxy);

		// 	r_temp[0]=r[0];
		// 	r_temp[1]=r[1];
		// 	r_temp[2]=rxy;

		// 	res=computeMI_knml(datafactors[0],datafactors[1],xy_factors,r_temp,n,c2terms, looklog);
		// 	res[1]-=sc_comb/n;

		// 	printf("%d: I_xz=%lf Ik_xz=%lf\n",stop,res[0],res[1]);
		// 	for(ll=0;ll<2;ll++) printf("r[%d]=%d ",ll,r[ll]);
		// 	printf("\n");fflush(stdout);

		// #endif


		double sc_levels_x(init_bin); // Number of levels of the first variable
		double sc_levels_y(init_bin); // Number of levels of the second variable
		int rx;
		int ry;
		for(stop=1;stop<STEPMAX;stop++)
		{

			///////////////////////////////////////////

			sc_comb=0;

			if(ptr_cnt[ptrVarIdx[0]]==1){

				//optimize on x
				//I(x;y)


				factors1[0]=datafactors[1];//y
				factors1[1]=singlefactor;// one single bin for x initially //herve

				rt1[0]=r[1];//xy
				rt1[1]=1; // one single bin for x initially //herve
				if (stop == 1){
					sc_levels_y = init_bin;
					rt1[0] = init_bin;
				}
				//else sc_levels_y = 0.5*(rt1[0]-1);
				else sc_levels_y = (sc_levels_y*(stop-1) + rt1[0])/stop; //Harmonic mean of previous levels.
				sc = 0.5*(sc_levels_y-1);

				#if _MY_DEBUG_MInoU
					printf("start optfun\n ");
					fflush(stdout);
				#endif

				rx = r[0];
				pxy = 1; //herve
				sc_levels1 = sc_levels_y;
				sc_levels2 = sc_levels_x;
				//sc_levels1 = r[1];
				//sc_levels2 = r[0];
				// Run optimization on X.
				//MInew=optfun_onerun_kmdl_coarse(sortidx[ptrVarIdx[0]], data[ptrVarIdx[0]],1, factors1, rt1, optfun1, sc, sc_levels_y, n, AllLevels[ptrVarIdx[0]], cut[0], &(r[0]), maxbins, looklog, looklbc, lookH, sc_look, cplx);
				MInew=optfun_onerun_kmdl_coarse(sortidx[ptrVarIdx[0]], data[ptrVarIdx[0]],2, factors1, rt1, pxy, sc, sc_levels1, sc_levels2, n, AllLevels[ptrVarIdx[0]], cut[0], &(r[0]), maxbins, looklog, looklbc, lookH, sc_look, cplx); // 2 factors //herve

				//update_datafactors(sortidx, ptrVarIdx[0], datafactors, 0, n, cut);

				if(AllLevels[ptrVarIdx[0]] > maxbins) np=maxbins;
				else np=AllLevels[ptrVarIdx[0]];

				//sc_comb += COEFF_COMB*logdbico(np,r[0]-1,looklbc);

			}

			////////////////////////////////////////////////

			if(ptr_cnt[ptrVarIdx[1]]==1){


				//opt y
				//I(x;y)


				factors1[0]=datafactors[0];//x
				factors1[1]=singlefactor;// one single bin for y initially //herve

				//rt1[0] = r[0];//xy
				rt1[0] = rx;
				rt1[1]=1; // one single bin for y initially //herve
				if (stop == 1){
					sc_levels_x = init_bin;
					rt1[0] = init_bin;
				}
				//else sc_levels_x = 0.5*(rt1[0]-1);
				else sc_levels_x = (sc_levels_x*(stop-1) + rt1[0])/stop; //Harmonic mean of previous levels.
				sc = 0.5*(sc_levels_x-1);

				ry = r[1];
				pxy = 1; //herve
				sc_levels1 = sc_levels_x;
				sc_levels2 = sc_levels_y;
				//sc_levels1 = rx;
				//sc_levels2 = ry;
				// Run optimization on Y.
				//MInew=optfun_onerun_kmdl_coarse(sortidx[ptrVarIdx[1]], data[ptrVarIdx[1]],1, factors1, rt1, optfun1, sc, sc_levels_x, n, AllLevels[ptrVarIdx[1]], cut[1], &(r[1]), maxbins, looklog, looklbc, lookH, sc_look, cplx);
				MInew=optfun_onerun_kmdl_coarse(sortidx[ptrVarIdx[1]], data[ptrVarIdx[1]],2, factors1, rt1, pxy, sc,
												sc_levels1, sc_levels2, n, AllLevels[ptrVarIdx[1]], cut[1], &(r[1]),
												maxbins, looklog, looklbc, lookH, sc_look, cplx); // 2 factors //herve



				update_datafactors(sortidx, ptrVarIdx[1], datafactors, 1, n, cut);
				update_datafactors(sortidx, ptrVarIdx[0], datafactors, 0, n, cut);

				//
				if(AllLevels[ptrVarIdx[1]] > maxbins) np=maxbins;
				else np=AllLevels[ptrVarIdx[1]];

				//sc_comb += COEFF_COMB*logdbico(np,r[1]-1, looklbc);

			}
			// Save cut points
			for(j=0; j<maxbins; j++){
				iterative_cuts[stop-1][j] = cut[0][j];
				iterative_cuts[stop-1][j+maxbins] = cut[1][j];
			}

			//////////////////////////////////////////

			jointfactors_u(datafactors,ptr,n, 2, r, xy_factors, &rxy,&Hxy);

			r_temp[0]=r[0];
			r_temp[1]=r[1];
			r_temp[2]=rxy;

			if (cplx==0) //MDL
				res=computeMI_kmdl(datafactors[0],datafactors[1],xy_factors,r_temp,n, looklog);
			else if(cplx == 1) //NML
				res=computeMI_knml(datafactors[0],datafactors[1],xy_factors,r_temp,n,c2terms, looklog);
			res[1]-=sc_comb/n;

			#if _MY_DEBUG_MInoU
			printf("%d: I_xz=%lf Ik_xz=%lf\n",stop,res[0],res[1]);
			for(ll=0;ll<2;ll++) printf("r[%d]=%d ",ll,r[ll]);
			printf("\n");fflush(stdout);
			#endif

			printf("%d: I_xz=%lf Ik_xz=%lf\n",stop,res[0],res[1]);

			for(i=stop-1;i>0;i--){
				if( (fabs(res[1]-MIk[i]) < EPS)){ // If no real improvement over last information  AND
					//(rx == r[0]) && (ry == r[1])) { // If the number of bins hasn't changed on both variables
					flag=1;
					Ik_av=MIk[i];
					I_av=MI[i];

					for(j=i+1;j<stop;j++){
						Ik_av+=MIk[j];
						I_av+=MI[j];
					}
					Ik_av/=(stop-i); //average over the periodic cycle
					I_av/=(stop-i);
					break;
				}
			}
			if(flag) break;
			MIk[stop]=res[1];
			MI[stop]=res[0];

		}//for
		iterative_cuts[stop][0]=-1; // mark where we stopped iterating
		iterative_cuts[stop][maxbins]=-1;

		if(flag){
		res[0]=I_av;
		res[1]=Ik_av;
		}

		#if _MY_DEBUG_MInoU
			printf("final: I_xy=%lf Ik_xy=%lf\n",res[0],res[1]);

			printf("0 : r=%d ",r[0]);
				for(i=0;i<r[0];i++){
				printf("%d ",cut[0][i]);
			}
			printf("\n");
			printf("1 : r=%d ",r[1]);
				for(i=0;i<r[1];i++){
				printf("%d ",cut[1][i]);
			}
			printf("\n");
		#endif

		// free memory

		free(res);
		free(xy_factors);
		free(r_temp);
		free(ptr);
		free(MI);
		free(MIk);

		//for(l=0;l<2;l++) free(factors1[l]);
		free(factors1);

		for(l=0;l<2;l++) free(datafactors[l]);
		free(datafactors);
		free(singlefactor);

		free(rt1);

	return iterative_cuts;
}


int** compute_Ixy_cond_u_alg1(int** data, int** sortidx, int* ptr_cnt, int* ptrVarIdx,  int* AllLevels, int nbrUi,
								int n, int maxbins, int **cut, int *r, double* c2terms, int init_nbin,
								double* looklog, double** looklbc, double* lookH, double** sc_look, int cplx)
{

	int j,l,ll;
	int STEPMAX1=10;

	/////////////////////////////////////

	double* res=(double *)calloc(2,sizeof(double));//results res[0]->I,res[1]->I-k
	int** iterative_cuts = (int **)calloc(STEPMAX+1, sizeof(int*));
	for(int i=0; i<STEPMAX+1; i++){
		iterative_cuts[i] = (int *)calloc(maxbins*2, sizeof(int));
	}

	//////////////////////////////////
	//allocation factors  x y

	int **datafactors;
	datafactors=(int **)calloc((nbrUi+2),sizeof(int*));

	for(l=0;l<(nbrUi+2);l++){
		datafactors[l]=(int *)calloc(n,sizeof(int));
	}

	int* r_temp=(int *)calloc(3,sizeof(int));

	////////////////////////////////////
	//initialization of datafactors && sortidx


	for(l=0;l<(nbrUi+2);l++){ //compute datafactors based on the positions of cut points in vector <cut>
		if(ptr_cnt[ptrVarIdx[l]]==1){
			update_datafactors(sortidx, ptrVarIdx[l], datafactors, l, n, cut);
		}else{
			for(j=0;j<=n-1;j++){//discrete case
				datafactors[l][j]=data[ptrVarIdx[l]][j];
			}
		}
	}

	int **uiyxfactors;//({ui},{uiy}),{uix},{uiyx})

	uiyxfactors=(int **)calloc(4,sizeof(int*));
	for(l=0;l<4;l++){
		uiyxfactors[l]=(int *)calloc(n,sizeof(int));
	}
	int *ruiyx;
	ruiyx=(int *)calloc(4,sizeof(int));

	double sc;
	double sc_comb=0;//term complexity due to optimization on many models

 	double MInew;

	double* MI = (double*)calloc(STEPMAX, sizeof(double));
	double* MIk = (double*)calloc(STEPMAX, sizeof(double));
	double I_av,Ik_av;

	int stop1,stop,i,flag;

	double pxy; //herve
	int sc_levels1,sc_levels2; //herve
	//

	int **factors1=(int **)calloc(2,sizeof(int*));
	for(l=0;l<2;l++){
		factors1[l]=(int *)calloc(n,sizeof(int));
	}
	int **factors3=(int **)calloc(4,sizeof(int*));
	for(l=0;l<4;l++){
		factors3[l]=(int *)calloc(n,sizeof(int));
	}

	int *rt1=(int *)calloc(2,sizeof(int));
	int *rt3=(int *)calloc(4,sizeof(int));
	//int *optfun1=(int *)calloc(2,sizeof(int)); //herve

	//

	MI[0]=0;
	MIk[0]=-DBL_MAX;


	// #if _MY_DEBUG_MInoU

	// 	jointfactors_u(datafactors,ptr,n, 2, r, xy_factors, &rxy,&Hxy);

	// 	r_temp[0]=r[0];
	// 	r_temp[1]=r[1];
	// 	r_temp[2]=rxy;

	// 	res=computeMI_knml(datafactors[0],datafactors[1],xy_factors,r_temp,n,c2terms, looklog);
	// 	res[1]-=sc_comb/n;

	// 	printf("%d: I_xz=%lf Ik_xz=%lf\n",stop,res[0],res[1]);
	// 	for(ll=0;ll<2;ll++) printf("r[%d]=%d ",ll,r[ll]);
	// 	printf("\n");fflush(stdout);

	// #endif

	cout << "starting loop on stop1\n";

	double sc_levels_x; // Number of levels of the first variable
	double sc_levels_y; // Number of levels of the second variable
	int rx = init_nbin;
	int ry = init_nbin;
	for(stop1=1;stop1<STEPMAX1;stop1++)
	  {

	      //////////////////////////////////////////////////////////////////////////////////////
	      // optimize  pxy * [ I(y;xu) + I(x;yu) ] + (1 - pxy) [ I(y;u) + I(x;u) ]  over u

	    for(stop=1;stop<STEPMAX;stop++)
		{

		  cout << "starting loop on Ui, stop=" << stop << "\n";
		  flag=0;
		  int nb_ui_cnt=0;

		  for(l=0;l<nbrUi;l++){

		    if(ptr_cnt[ptrVarIdx[l+2]]==1){//continuous

			//compute joint factors u xu yu xyu

			jointfactors_uiyx(datafactors, l+2,  n, nbrUi, r, uiyxfactors, ruiyx);

			//init variables for the optimization run

			//sc=0.5*(r[0]-1)*(ruiyx[0])+0.5*(r[1]-1)*(ruiyx[0]); //herve

			factors3[0]=uiyxfactors[0];//u
			factors3[1]=uiyxfactors[1];//yu
			factors3[2]=uiyxfactors[2];//xu
			factors3[3]=uiyxfactors[3];//xyu

			rt3[0]=ruiyx[0];//u
			rt3[1]=ruiyx[1];//yu
			rt3[2]=ruiyx[2];//xu
			rt3[3]=ruiyx[3];//xyu

			sc_levels_x=rx; //herve
			sc_levels_y=ry; //herve

			pxy = 0.5; //herve can be varied from 0 to 1
			sc_levels1 = sc_levels_x; //herve
			sc_levels2 = sc_levels_y; //herve

			sc = 0.25*(sc_levels_x-1) + 0.25*(sc_levels_y-1); //herve for MDL (not quite correct but good enough?)


			//optimization run on ptrVarIdx[l+2]

			//MInew=optfun_onerun_kmdl_coarse(sortidx[ptrVarIdx[l+2]],  data[ptrVarIdx[l+2]],	3, factors3, rt3, optfun3, sc, n, AllLevels[ptrVarIdx[l+2]], cut[l+2], &(r[l+2]), maxbins, looklog, looklbc, lookH);
			cout << "calling onerun on U\n";
			MInew=optfun_onerun_kmdl_coarse(sortidx[ptrVarIdx[l+2]], data[ptrVarIdx[l+2]], 4, factors3, rt3, pxy,
											sc, sc_levels1, sc_levels2, n, AllLevels[ptrVarIdx[l+2]], cut[l+2],
											&(r[l+2]), maxbins, looklog, looklbc, lookH, sc_look, cplx); // 4 factors //herve

			cout << "done onerun\n";
			//update datafactors of ptrVarIdx[0]

			update_datafactors(sortidx, ptrVarIdx[l+2], datafactors, l+2, n, cut);

			#if _MY_DEBUG_MI
				printf("opt u: MInew=%lf ",MInew);
				for(ll=0;ll<nbrUi+2;ll++) printf("r[%d]=%d ",ll,r[ll]);
				printf("\n");fflush(stdout);
			#endif

			nb_ui_cnt++;

		    }

		  }
		  //compute joint factors u

		  //jointfactors_u(datafactors,ptrUiIdx,n, nbrUi, r, ufactors, &ru,&Hu);
			//////////////////////////////////////////

			//jointfactors_u(datafactors,ptr,n, 2, r, xy_factors, &rxy,&Hxy);
			//compute joint factors u xu yu xyu
			jointfactors_uiyx(datafactors,-1,  n, nbrUi, r, uiyxfactors, ruiyx); //herve
			cout << "jointfactors uiyx\n";
			//uiyxfactors[0];//u
			//uiyxfactors[1];//yu
			//uiyxfactors[2];//xu
 			//uiyxfactors[3];//xyu
			//ruiyx -> 0:u,1:uy,2:ux,3:uyx

      		double* res = (double*) calloc(2, sizeof(double));
      		double* res_x = (double*) calloc(2, sizeof(double));
      		double* res_y = (double*) calloc(2, sizeof(double));
      		double* res_xu = (double*) calloc(2, sizeof(double));
      		double* res_yu = (double*) calloc(2, sizeof(double));

			if (cplx == 0){ //MDL
			  r_temp[2]=ruiyx[3];
			  r_temp[0]=r[0];
			  r_temp[1]=ruiyx[1];
			  res_x=computeMI_kmdl(datafactors[0],uiyxfactors[1],uiyxfactors[3],r_temp,n, looklog);
              res[0] = pxy*res_x[0];
              res[1] = pxy*res_x[1];

			  r_temp[0]=r[1];
			  r_temp[1]=ruiyx[2];
			  res_y = computeMI_kmdl(datafactors[1],uiyxfactors[2],uiyxfactors[3],r_temp,n, looklog);
			  res[0] += pxy*res_y[0];
			  res[1] += pxy*res_y[1];

			  r_temp[2]=ruiyx[2];
			  r_temp[0]=r[0];
			  r_temp[1]=ruiyx[0];
			  res_xu = computeMI_kmdl(datafactors[0],uiyxfactors[0],uiyxfactors[2],r_temp,n, looklog);
			  res[0] += (1-pxy)*res_xu[0];
			  res[1] += (1-pxy)*res_xu[1];

			  r_temp[2]=ruiyx[1];
			  r_temp[0]=r[1];
			  r_temp[1]=ruiyx[0];
			  res_yu = computeMI_kmdl(datafactors[1],uiyxfactors[0],uiyxfactors[1],r_temp,n, looklog);
			  res[0] += (1-pxy)*res_yu[0];
			  res[1] += (1-pxy)*res_yu[1];

			}else if(cplx == 1){ //NML
			  r_temp[2]=ruiyx[3];
			  r_temp[0]=r[0];
			  r_temp[1]=0;
			  r_temp[1]=ruiyx[1]; //vincent
			  res_x = computeMI_knml(datafactors[0],uiyxfactors[1],uiyxfactors[3],r_temp,n,c2terms, looklog);
			  res[0] += pxy*res_x[0];
			  res[1] += pxy*res_x[1];

			  r_temp[0]=r[1];
			  r_temp[1]=0;
			  r_temp[1]=ruiyx[2]; //vincent
			  res_y = computeMI_knml(datafactors[1],uiyxfactors[2],uiyxfactors[3],r_temp,n,c2terms, looklog);
			  res[0] += pxy*res_y[0];
			  res[1] += pxy*res_y[1];

			  r_temp[2]=ruiyx[2];
			  r_temp[0]=r[0];
			  r_temp[1]=ruiyx[0];
			  res_xu = computeMI_knml(datafactors[0],uiyxfactors[0],uiyxfactors[2],r_temp,n, c2terms, looklog);
			  res[0] += (1-pxy)*res_xu[0];
			  res[1] += (1-pxy)*res_xu[1];

			  r_temp[2]=ruiyx[1];
			  r_temp[0]=r[1];
			  r_temp[1]=ruiyx[0];
			  res_yu = computeMI_knml(datafactors[1],uiyxfactors[0],uiyxfactors[1],r_temp,n, c2terms, looklog);
			  res[0] += (1-pxy)*res_yu[0];
			  res[1] += (1-pxy)*res_yu[1];
			}
			res[1]-=sc_comb/n;

			//#if _MY_DEBUG_MInoU
			printf("%d: I_xz=%lf Ik_xz=%lf\n",stop,res[0],res[1]);
			for(ll=0;ll<(nbrUi+2);ll++) printf("r[%d]=%d ",ll,r[ll]);
			printf("\n");fflush(stdout);
			//#endif

			cout << "compute res\n";

			for(i=stop-1;i>0;i--){
				if( (fabs(res[1]-MIk[i]) < EPS)){ // If no real improvement over last information  AND
					//(rx == r[0]) && (ry == r[1])) { // If the number of bins hasn't changed on both variables
					flag=1;
					Ik_av=MIk[i];
					I_av=MI[i];

					for(j=i+1;j<stop;j++){
						Ik_av+=MIk[j];
						I_av+=MI[j];
					}
					Ik_av/=(stop-i); //average over the periodic cycle
					I_av/=(stop-i);
					break;
				}
			}
			cout << "testing res\n";
			if(flag) break;
			MIk[stop]=res[1];
			MI[stop]=res[0];

		}// end stop for u

        // Reset MIk and MI
		for(i=0;i<STEPMAX;i++){
	        MIk[i]=0;
	        MI[i]=0;
        }

		cout << "starting loop on X,Y\n";
	      for(stop=1;stop<STEPMAX;stop++) //optimize I(y;xu) over x and I(x;yu) over y
		{
			///////////////////////////////////////////

		    flag=0;
			sc_comb=0;


			//optimize I(y;xu) over x and I(x;yu) over y

			//compute joint factors u xu yu xyu
			cout << "jointfactors_uiyx\n";
			if(ptr_cnt[ptrVarIdx[0]]==1 || ptr_cnt[ptrVarIdx[1]]==1){
			  jointfactors_uiyx(datafactors,-1,  n, nbrUi, r, uiyxfactors, ruiyx);
			}

			if(ptr_cnt[ptrVarIdx[0]]==1){

				//optimize on x
				//I(y;xu)

 			        //compute joint factors u xu yu xyu
			        //jointfactors_uiyx(datafactors,0,  n, nbrUi, r, uiyxfactors, ruiyx);

				//init variables for the optimization run

				//sc=0.5*(r[1]-1)*(ruiyx[0]); // herve


				//factors1[0]=datafactors[1];//y
				//factors1[1]=singlefactor;// one single bin for x initially //herve
				factors1[0]=uiyxfactors[1];//yu
				factors1[1]=uiyxfactors[0];//u //herve

				//optfun1[0]=-1;//xy
				//optfun1[1]=1;//x
				rt1[0]=r[1];//y
				rt1[1]=1; // one single bin for x initially //herve
				if (stop == 1){
					sc_levels_y = init_nbin;
					rt1[0] = init_nbin;
				}
				//else sc_levels_y = 0.5*(rt1[0]-1);
				else sc_levels_y = (sc_levels_y*(stop-1) + rt1[0])/stop; //Harmonic mean of previous levels.
				sc = 0.5*(sc_levels_y-1);

				#if _MY_DEBUG_MInoU
					printf("start optfun\n ");
					fflush(stdout);
				#endif

				rx = r[0];
				pxy = 1; //herve
				sc_levels1 = sc_levels_y; //herve
				sc_levels2 = sc_levels_x; //herve
				// Run optimization on X.
				//MInew=optfun_onerun_kmdl_coarse(sortidx[ptrVarIdx[0]], data[ptrVarIdx[0]],1, factors1, rt1, optfun1, sc, sc_levels_y, n, AllLevels[ptrVarIdx[0]], cut[0], &(r[0]), maxbins, looklog, looklbc, lookH, sc_look, cplx);

			cout << "calling onerun for X\n";
			cout << sc_levels1 << "\t" << sc_levels2 << "\n";
				MInew=optfun_onerun_kmdl_coarse(sortidx[ptrVarIdx[0]], data[ptrVarIdx[0]],2, factors1, rt1, pxy, sc, sc_levels1, sc_levels2, n, AllLevels[ptrVarIdx[0]], cut[0], &(r[0]), maxbins, looklog, looklbc, lookH, sc_look, cplx); // 2 factors //herve
			cout << "done onerun for X\n";

				//update_datafactors(sortidx, ptrVarIdx[0], datafactors, 0, n, cut);

				//sc_comb += COEFF_COMB*logdbico(np,r[0]-1,looklbc);

			}

			////////////////////////////////////////////////

			if(ptr_cnt[ptrVarIdx[1]]==1){


				//opt y
				//I(x;yu)

 			        //compute joint factors u xu yu xyu
			        //jointfactors_uiyx(datafactors,1,  n, nbrUi, r, uiyxfactors, ruiyx);

				//init variables for the optimization run

				//sc=0.5*(r[0]-1)*(ruiyx[0]); //herve


				//factors1[0]=datafactors[0];//x
				//factors1[1]=singlefactor;// one single bin for y initially //herve
				factors1[0]=uiyxfactors[2];//xu
				factors1[1]=uiyxfactors[0];//u //herve

				//optfun1[0]=-1;//xy
				//optfun1[1]=1;//y
				//rt1[0] = r[0];//xy
				rt1[0] = rx;
				rt1[1]=1; // one single bin for y initially //herve
				if (stop == 1){
					sc_levels_x = init_nbin;
					rt1[0] = init_nbin;
				}
				//else sc_levels_x = 0.5*(rt1[0]-1);
				else sc_levels_x = (sc_levels_x*(stop-1) + rt1[0])/stop; //Harmonic mean of previous levels.
				sc = 0.5*(sc_levels_x-1);

				ry = r[1];
				pxy = 1; //herve
				sc_levels1 = sc_levels_x; //herve
				sc_levels2 = sc_levels_y; //herve
				// Run optimization on Y.
				//MInew=optfun_onerun_kmdl_coarse(sortidx[ptrVarIdx[1]], data[ptrVarIdx[1]],1, factors1, rt1, optfun1, sc, sc_levels_x, n, AllLevels[ptrVarIdx[1]], cut[1], &(r[1]), maxbins, looklog, looklbc, lookH, sc_look, cplx);
			cout << "calling onerun for Y\n";
				MInew=optfun_onerun_kmdl_coarse(sortidx[ptrVarIdx[1]], data[ptrVarIdx[1]],2, factors1, rt1, pxy, sc, sc_levels1, sc_levels2, n, AllLevels[ptrVarIdx[1]], cut[1], &(r[1]), maxbins, looklog, looklbc, lookH, sc_look, cplx); // 2 factors //herve



				update_datafactors(sortidx, ptrVarIdx[1], datafactors, 1, n, cut);
				update_datafactors(sortidx, ptrVarIdx[0], datafactors, 0, n, cut);

				//

				//sc_comb += COEFF_COMB*logdbico(np,r[1]-1, looklbc);

			}
			// Save cut points
			for(j=0; j<maxbins; j++){
				iterative_cuts[stop-1][j] = cut[0][j];
				iterative_cuts[stop-1][j+maxbins] = cut[1][j];
			}

			//////////////////////////////////////////

			//jointfactors_u(datafactors,ptr,n, 2, r, xy_factors, &rxy,&Hxy);
			//compute joint factors u xu yu xyu
			jointfactors_uiyx(datafactors,-1,  n, nbrUi, r, uiyxfactors, ruiyx); //herve
			//uiyxfactors[0];//u
			//uiyxfactors[1];//yu
			//uiyxfactors[2];//xu
 			//uiyxfactors[3];//xyu
			//ruiyx -> 0:u,1:uy,2:ux,3:uyx

      		double* res = (double*) calloc(2, sizeof(double));
      		double* res_y = (double*) calloc(2, sizeof(double));


			if (cplx == 0){ //MDL //cplx
			  r_temp[2]=ruiyx[3];
			  r_temp[0]=r[0];
			  r_temp[1]=ruiyx[1];
			  res=computeMI_kmdl(datafactors[0],uiyxfactors[1],uiyxfactors[3],r_temp,n, looklog);

			  r_temp[0]=r[1];
			  r_temp[1]=ruiyx[2];
			  res_y = computeMI_kmdl(datafactors[1],uiyxfactors[2],uiyxfactors[3],r_temp,n, looklog);
			  res[0] += res_y[0];
			  res[1] += res_y[1];

			}else if(cplx == 1){ //NML
			  r_temp[2]=ruiyx[3];
			  r_temp[0]=r[0];
			  r_temp[1]=0;
			  res=computeMI_knml(datafactors[0],uiyxfactors[1],uiyxfactors[3],r_temp,n,c2terms, looklog);

			  r_temp[0]=r[1];
			  r_temp[1]=0;
			  res_y = computeMI_knml(datafactors[1],uiyxfactors[2],uiyxfactors[3],r_temp,n,c2terms, looklog);
			  res[0] += res_y[0];
			  res[1] += res_y[1];
			}
			//res[1]-=sc_comb/n;

			#if _MY_DEBUG_MInoU
			printf("%d: I_xz=%lf Ik_xz=%lf\n",stop,res[0],res[1]);
			for(ll=0;ll<2;ll++) printf("r[%d]=%d ",ll,r[ll]);
			printf("\n");fflush(stdout);
			#endif

			printf("%d: I_xz=%lf Ik_xz=%lf\n",stop,res[0],res[1]);

			for(i=stop-1;i>0;i--){
				if( (fabs(res[1]-MIk[i]) < EPS) || // If no real improvement over last information  AND
					((rx == r[0]) && (ry == r[1])) ) { // If the number of bins hasn't changed on both variables
					flag=1;
					Ik_av=MIk[i];
					I_av=MI[i];

					for(j=i+1;j<stop;j++){
						Ik_av+=MIk[j];
						I_av+=MI[j];
					}
					Ik_av/=(stop-i); //average over the periodic cycle
					I_av/=(stop-i);
					break;
				}
			}
			if(flag) break;
			MIk[stop]=res[1];
			MI[stop]=res[0];

		}//for
		iterative_cuts[stop][0]=-1; // mark where we stopped iterating
		iterative_cuts[stop][maxbins]=-1;

		if(flag){
		res[0]=I_av;
		res[1]=Ik_av;
		}


		#if _MY_DEBUG_MInoU
			printf("final: I_xy=%lf Ik_xy=%lf\n",res[0],res[1]);

			printf("0 : r=%d ",r[0]);
				for(i=0;i<r[0];i++){
				printf("%d ",cut[0][i]);
			}
			printf("\n");
			printf("1 : r=%d ",r[1]);
				for(i=0;i<r[1];i++){
				printf("%d ",cut[1][i]);
			}
			printf("\n");
		#endif


		if(nbrUi==0) {
			break;
		}


	  }//end stop1

	// free memory

	free(r_temp);

	for(l=0;l<(2);l++) free(datafactors[l]);
	free(datafactors);

	free(factors1);
	free(rt1);
	//free(optfun1);

	return iterative_cuts;
}

//////////////////////////////////////////////////////////////////////////////////////////////
//compute I(x,y|u), conditional mutual information between two variables conditioned on a set u

//algorithm (1)
//
//initialization:
//optimize x optimization function=I(x,u)
//optimize y optimization function=I(y,u)
//optimize over the variables in u, optimization function= I(x,u)+I(y,u)
//iteration :
//optimize on x (y) optimization function=I(x,yu) ( I(y,xu) ) -> repeat until convergence


//

//int* ptrVarIdx // vector of length nbrUi+2 with indexes of the variables x,y,{u}
//int **cut // vectors with the positions of the cuts
//int *r // number of levels of each variables


int** compute_mi_cond_alg1(int** data, double** dataDouble, int** sortidx,  int* AllLevels, int* ptr_cnt, int* ptrVarIdx,
						   int nbrUi, int n, int maxbins, double* c2terms, int init_bin,
						   double* looklog, double** looklbc, double* lookH, double** sc_look, int cplx)
{

	int** iterative_cuts = (int **)calloc(STEPMAX, sizeof(int*));
	for(int i=0; i<STEPMAX; i++){
		iterative_cuts[i] = (int *)calloc(maxbins*2, sizeof(int));
	}

	int j,k,l,ll;

	int *r=(int *)calloc((nbrUi+2),sizeof(int));

	int **cut;
	cut=(int **)calloc(nbrUi+2,sizeof(int*));
	for(l=0;l<(nbrUi+2);l++){
		cut[l]=(int *)calloc(maxbins,sizeof(int));
	}
	int init_nbin=init_bin;

	int lbin=floor(n/init_nbin);
	if(lbin<1) {
		lbin=1;
		init_nbin=n;
	}

	double* wbin=(double *)calloc(nbrUi+2, sizeof(double)); //unit width of bin for both variables
	for(j=0; j<(nbrUi+2); j++){
		wbin[j] = (dataDouble[sortidx[j][n-1]][j] - dataDouble[sortidx[j][0]][j]) / init_bin;
	}


	// if(isNumberOfLevelsLessTwo(data, nbrUi, n, 2)){
	// 	res[0]= 0;
	// 	res[1]= 0;
	// 	res[2]= n;
	// 	return res;
	// }

	// updateContinuousVariables(data, AllLevels,ptr_cnt, nbrUi, n, 2);

	// ///////////////////////////////////////////////////////////////////////////////////
	// no conditioning, empty set of variables in u
	// calling funcion compute_Ixy_alg1
	// NO u
	if(nbrUi==0){
		//////////////////////////////////////
		// initialization cut and r

		for(l=0;l<2;l++){
			if(ptr_cnt[ptrVarIdx[l]]==1){

				if(!INIT_EQUAL_WIDTH){
					for(j=0;j<init_nbin-1;j++) {
						cut[l][j]=j*lbin+lbin-1;
					}
				}
				if(INIT_EQUAL_WIDTH){
					int dataDoubleit=0;
					double cut_value;
					double data_value;
					for(j=0;j<init_nbin-1;j++) {
						cut_value = dataDouble[sortidx[l][0]][l] + j*wbin[l]+wbin[l];
						data_value = dataDouble[sortidx[l][dataDoubleit]][l];
						while((data_value < cut_value) && (dataDoubleit < n-1)){
							data_value = dataDouble[sortidx[l][dataDoubleit]][l];
							dataDoubleit++;
						}
						cut[l][j] = dataDoubleit-1;
					}
				}
				cut[l][init_nbin-1]=n-1;
				r[l]=init_nbin;
			}else{
				r[l]=AllLevels[ptrVarIdx[l]];
			}
		}


		#if _MY_DEBUG_MI
				printf("	compute_Ixy_alg1.. \n ");
				printf("	cnt: %d %d\n",ptr_cnt[ptrVarIdx[0]],ptr_cnt[ptrVarIdx[1]]);
				fflush(stdout);
		#endif
		iterative_cuts = compute_Ixy_alg1(data, sortidx, ptr_cnt, ptrVarIdx, AllLevels, n,
										  maxbins, cut, r, c2terms, init_bin, looklog, looklbc, lookH, sc_look, cplx);


		// free memory

		#if _MY_DEBUG_MI
				printf("	free cut nbr=%d \n ",nbrUi);
				fflush(stdout);
		#endif

		#if _MY_DEBUG_MI
				printf("	return \n ");
				fflush(stdout);
		#endif

		free(r);
		for(l=0;l<2;l++){
			free(cut[l]);
		}
		free(cut);
		free(wbin);


		return iterative_cuts;
	}

	else{
		//////////////////////////////////////
		// initialization cut and r

		for(l=0;l<(nbrUi+2);l++){
			if(ptr_cnt[ptrVarIdx[l]]==1){

				if(!INIT_EQUAL_WIDTH){
					for(j=0;j<init_nbin-1;j++) {
						cut[l][j]=j*lbin+lbin-1;
					}
				}
				if(INIT_EQUAL_WIDTH){
					int dataDoubleit=0;
					double cut_value;
					double data_value;
					for(j=0;j<init_nbin-1;j++) {
						cut_value = dataDouble[sortidx[l][0]][l] + j*wbin[l]+wbin[l];
						data_value = dataDouble[sortidx[l][dataDoubleit]][l];
						while((data_value < cut_value) && (dataDoubleit < n-1)){
							data_value = dataDouble[sortidx[l][dataDoubleit]][l];
							dataDoubleit++;
						}
						cut[l][j] = dataDoubleit-1;
					}
				}
				cut[l][init_nbin-1]=n-1;
				r[l]=init_nbin;
			}else{
				r[l]=AllLevels[ptrVarIdx[l]];
			}
		}


		cout << " calling micond u\n";
		iterative_cuts = compute_Ixy_cond_u_alg1(data, sortidx, ptr_cnt, ptrVarIdx, AllLevels, nbrUi, n,
										  		 maxbins, cut, r, c2terms, init_bin, looklog, looklbc, lookH, sc_look, cplx);


		// free memory

		#if _MY_DEBUG_MI
				printf("	free cut nbr=%d \n ",nbrUi);
				fflush(stdout);
		#endif

		#if _MY_DEBUG_MI
				printf("	return \n ");
				fflush(stdout);
		#endif

		free(r);
		for(l=0;l<2;l++){
			free(cut[l]);
		}
		free(cut);
		free(wbin);


		return iterative_cuts;
	}

}




//####################################################################################################
//# Dealing with input variables

bool transformToFactorsContinuous(double** data, int** dataNumeric, int i, int n){

	std::multimap<double,int> myMap;

	//clean the dictionary since it is used column by column
	myMap.clear();
	// myMap["NA"] = -1;
	// myMap[""] = -1;

	vector <double> clmn;
	for(int j = 0; j < n; j++){
		double entry = data[j][i];
    clmn.push_back(entry);
	}

	sort(clmn.begin(), clmn.end());

	for(int j = 0; j < clmn.size(); j++){
			myMap.insert(pair<double,int>(clmn[j],j));
	}
	// for (std::map<double,int>::iterator it=myMap.begin(); it!=myMap.end(); ++it){
	// 	cout << "(" << it->first << "," << it->second << ")\n";
	// }


 	for(int j = 0; j < n; j++){
 		double entry = data[j][i];
    dataNumeric[j][i] = myMap.find(entry)->second;

    typedef std::multimap<double, int>::iterator iterator;
    std::pair<iterator, iterator> iterpair = myMap.equal_range(entry);
    iterator it = iterpair.first;
    for (; it != iterpair.second; ++it) {
        if (it->second == dataNumeric[j][i]) {
            myMap.erase(it);
            break;
        }
    }
	}
}

/*
 * Transforms the string into factors
 */
bool transformToFactors(double** data, int** dataNumeric, int n, int i){
	 // create a dictionary to store the factors of the strings
 	map<double,int> myMap;

	//clean the dictionary since it is used column by column
	myMap.clear();
	int factor = 0;

 	for(int j = 0; j < n; j++){
		map<double,int>::iterator it = myMap.find(data[j][i]);
		if ( it != myMap.end() ){
			dataNumeric[j][i] = it->second;
		}
		else {
			myMap[data[j][i]] = factor;
			dataNumeric[j][i] = factor;
			factor++;

		}
	}

}

bool transformToFactorsContinuousIdx(int** dataNumeric, int** dataNumericIdx, int n, int i){

	map<int,int> myMap;

	//clean the dictionary since it is used column by column
	myMap.clear();

	// vector <int> clmn;
	for(int j = 0; j < n; j++){
		int entry = dataNumeric[j][i];
		if (entry != -1) myMap[entry] = j;
	}

	int j = 0;
	for (std::map<int,int>::iterator it=myMap.begin(); it!=myMap.end(); ++it){
		dataNumericIdx[i][j] =  it->second;
		j++;
	}

}

extern "C" SEXP mydiscretizeMutual(SEXP RmyDist1, SEXP RmyDist2, SEXP RflatU, SEXP RnbrU, SEXP RmaxBins, SEXP Rinitbin, SEXP Rcplx){

	std::vector<double> myDist1Vec = Rcpp::as< vector <double> >(RmyDist1);
	std::vector<double> myDist2Vec = Rcpp::as< vector <double> >(RmyDist2);
	int maxbins = Rcpp::as<int> (RmaxBins);
	int init_bin = Rcpp::as<int> (Rinitbin);
	int cplx = Rcpp::as<int> (Rcplx);
	int nbrU = Rcpp::as<int> (RnbrU);
	int n = myDist1Vec.size();

	double* myDist1 = &myDist1Vec[0];
	double* myDist2 = &myDist2Vec[0];

	cout << "start...\n";
	std::vector<double> vectorflatU = Rcpp::as< vector <double> > (RflatU);
	cout << "vectorflatU\n";
	double** matrixU = new double*[nbrU];
	for(int l=0; l<nbrU; l++){
		matrixU[l] = new double[n];
	}
	cout << "declared matrixU\n";
	if(nbrU > 0){
		for(int l=0; l<nbrU; l++){
			matrixU[l] = new double[n];
			copy(vectorflatU.begin()+(l*n), vectorflatU.begin()+((l+1)*n), matrixU[l]);
		}
	}
	cout << "copied flatU\n";


	// data ##########################################################################################
	//create the data matrix for factors
	int i; int j; int l;
	double** data= new double*[n];
	for(i = 0; i < n; i++){
    	data[i] = new double[nbrU+2];
		data[i][0] = myDist1[i];
		data[i][1] = myDist2[i];
		for(l=0; l<nbrU; l++){
			data[i][l+2] = matrixU[l][i];
		}
	}
	cout << "filled data\n";

	int** dataNumeric = new int*[n];
	for(j = 0; j < n; j++){
		dataNumeric[j] = new int[nbrU+2];
	}
	for(i = 0; i < (nbrU+2); i++){
	    transformToFactorsContinuous(data, dataNumeric, i, n);
	}

	// sortidx #######################################################################################
	// for continuous non all gaussians
	//create the data matrix for factors indexes
	int** dataNumericIdx = new int*[nbrU+2];
	for(i = 0; i < (nbrU+2); i++){
	  dataNumericIdx[i] = new int[n];
	  for(j = 0; j < n; j++){
	    dataNumericIdx[i][j]=-1;
	  }
	}
	cout << "filled dataNumericIdx\n";

	for(i = 0; i < (nbrU+2); i++){
	  transformToFactorsContinuousIdx(dataNumeric, dataNumericIdx, n, i);
	  transformToFactors(data, dataNumeric, n, i);//update environment.dataNumeric taking into account repetition
	}
	cout << "double factors\n";

	// AllLevels #####################################################################################
	int* AllLevels = new int[nbrU+2];
	for(i = 0; i < (nbrU+2); i++){
		AllLevels[i] = n;
	}

	// ptr_cnt #######################################################################################
	int* ptr_cnt = new int[nbrU+2];
	for(i = 0; i < (nbrU+2); i++){
		ptr_cnt[i] = 1;
	}

	// nrbUi #########################################################################################
	int nbrUi = nbrU;

	// c2terms #######################################################################################
	double* c2terms = new double[n+1];
	for(int i = 0; i < n+1; i++){
		c2terms[i] = -1;
	}

	// lookup tables #################################################################################
	double* looklog = new double[n+2];
	looklog[0] = 0.0;
	for(int i = 1; i < n+2; i++){
	  looklog[i] = log(1.0*i);
	}

	double* lookH = new double[n+2];
	lookH[0] = 0.0;
	for(int i = 1; i < n+2; i++){
	  lookH[i] = i*looklog[i]-(i+1)*looklog[(i+1)];
	}
	double** looklbc = new double*[n+1];
	for(int k = 0; k < n+1; k++){
	  looklbc[k] = new double[maxbins+1];
	  for(int z = 0; z < maxbins+1; z++){
	    looklbc[k][z]=0;
	  }
	}


	// ###############################################################################################
	// ###############################################################################################
	// ###############################################################################################
	// ###############################################################################################
	// Declare tables_red ############################################################################

	int* posArray = new int[nbrU+2];
	for(i = 0; i < (nbrU+2); i++){
		posArray[i] = i;
	}

	int nbrRetValues = 3;

	int** dataNumeric_red ;//prograssive data rank with repetition for same values
	int** dataNumericIdx_red ;//index of sorted data
	int* AllLevels_red ;//number of levels
	int* cnt_red ;//bool continuous or not
	int* posArray_red ;//node references

	double* res_new;

	int samplesNotNA = 0;
	int* samplesToEvaluate = new int[n];
	int* samplesToEvaluateTamplate = new int[n];
	bool cnt = true;
	for(int i = 0; i < n; i++){
    samplesToEvaluate[i] = 1;
    if(i!=0)
      samplesToEvaluateTamplate[i] = samplesToEvaluateTamplate[i-1];
    else
      samplesToEvaluateTamplate[i] = 0;

    cnt = true;
		for(int j = 0; (j < nbrU+2) && (cnt); j++){
      if(dataNumeric[i][posArray[j]] == -1){
        cnt = false;
        break;
      }
    }
    if(cnt==false){
      samplesToEvaluate[i] = 0;
      samplesToEvaluateTamplate[i] ++ ;
    }
		else samplesNotNA++;
  }

  ////////////////////////////////////////////////////////////

  if(samplesNotNA <= 2){
    res_new = new double[3];
    res_new[0]=samplesNotNA;//N
    res_new[1]=0;//Ixyu
    res_new[2]=1;//cplx
  }

  else{

    ////////////////////////////////////////////////////////
    // allocate data reducted *_red
    // all *_red variables are passed to the optimization routine

    // not using memory space

    dataNumericIdx_red = new int*[nbrU+2];
    dataNumeric_red = new int*[nbrU+2];
    for(int j = 0; j < (nbrU+2); j++){
      dataNumericIdx_red[j] = new int[samplesNotNA];
      dataNumeric_red[j] = new int[samplesNotNA];
    }

    AllLevels_red = new int[nbrU+2];
    cnt_red = new int[nbrU+2];
    posArray_red = new int[nbrU+2];


    ////////////////////////////////////////////////////////
    // copy data to evaluate in new vectors to pass to optimization

    int k1,k2,si;

    int nnr;// effective number of not repeated values
    int prev_val;

    // cout << "lev bef:" << cnt_red[0] << " " << cnt_red[1] << endl;

    for(int j = 0; j < (nbrU+2); j++){
      posArray_red[j]=j;
      AllLevels_red[j]=AllLevels[posArray[j]];
      cnt_red[j]=ptr_cnt[posArray[j]];

      k1=0;
      k2=0;

      nnr=0;
      prev_val= -1;
      for(int i = 0; i < n; i++){
        if(	samplesToEvaluate[i] == 1){
          dataNumeric_red[j][k1]=dataNumeric[i][posArray[j]];
          k1++;
        }
        if(cnt_red[j] == 1){
          si=dataNumericIdx[posArray[j]][i];
          if(samplesToEvaluate[si] == 1){

            dataNumericIdx_red[j][k2] = si - samplesToEvaluateTamplate[si];
            k2++;

            if(dataNumeric[si][posArray[j]] != prev_val){//check whether is a different values or repeated
              nnr++;
              prev_val=dataNumeric[si][posArray[j]];
            }


          }
        }
      }

      if(cnt_red[j] == 1){
        if(isContinuousDiscrete(dataNumeric_red, samplesNotNA, j)){
          cnt_red[j] = 0;
        }
      }


      // update with the effective number of levels
      if(cnt_red[j] == 1)
        AllLevels_red[j]=nnr;
    }

    //check effective number of levels variables continuous enough to be treat as continuous?
    if(samplesNotNA != n){
      //updateNumberofLevelsAndContinuousVariables(dataNumeric_red, AllLevels_red, cnt_red, myNbrUi, samplesNotNA);
    }

			// cout << "lev:" << cnt_red[0] << " " << cnt_red[1] << endl;

  }

  // Declare the lookup table of the parametric complexity
  double** sc_look = new double*[maxbins+1];
  for(int K = 0; K < (maxbins+1); K++){
    sc_look[K] = new double[n+1];
    for(int i = 0; i < (n+1); i++){
          if(K==1) sc_look[K][i] = 0;
          else if (i==0) sc_look[K][i] = 0;
          else sc_look[K][i] = -1;
	}
  }
  for(int i=0; i<n; i++){
	  double d = computeLogC(i, 2, sc_look); // Initialize the c2 terms
  }


  int** iterative_cuts = compute_mi_cond_alg1(dataNumeric_red, data, dataNumericIdx_red, AllLevels_red,
  											  cnt_red, posArray_red, nbrUi, n, maxbins, c2terms, init_bin,
											  looklog, looklbc, lookH, sc_look, cplx);

  int* iterative_cutpoints1 = new int[STEPMAX*maxbins];
  int* iterative_cutpoints2 = new int[STEPMAX*maxbins];
  int niterations=0;
  int ncutpoints1;
  int ncutpoints2;
  for(int l=0; l<STEPMAX+1; l++){
	if((iterative_cuts[l][0]==-1) && (iterative_cuts[l][maxbins]==-1) ){
		niterations=l;
		break;
	}
	ncutpoints1=1;
	ncutpoints2=1;
	i=0;
	while(iterative_cuts[l][i] < iterative_cuts[l][i+1]){
		iterative_cutpoints1[maxbins*l+i] = iterative_cuts[l][i];
		i++;
	}
	for(int j=i; j<maxbins; j++){
		iterative_cutpoints1[maxbins*l+j] = -1;
	}
	i=0;
	while(iterative_cuts[l][i+maxbins] < iterative_cuts[l][i+maxbins+1]){
		iterative_cutpoints2[maxbins*l+i] = iterative_cuts[l][i+maxbins];
		i++;
	}
	for(int j=i; j<maxbins; j++){
		iterative_cutpoints2[maxbins*l+j] = -1;
	}
  }
  vector<double> cutpoints1(iterative_cutpoints1, (iterative_cutpoints1 + niterations*maxbins));
  vector<double> cutpoints2(iterative_cutpoints2, (iterative_cutpoints2 + niterations*maxbins));
  // structure the output
  List result = List::create(
    _["cutpoints1"] = cutpoints1,
    _["cutpoints2"] = cutpoints2
  ) ;

  for(int i=0; i<(STEPMAX+1); i++){
	  free(iterative_cuts[i]);
  }
  free(iterative_cuts);

  for(int i=0; i<(maxbins+1); i++){
	  delete[] sc_look[i];
  }
  delete[] sc_look;
  delete[] c2terms;

  for(int i=0; i<n; i++){
	  delete [] data[i];
	  delete [] dataNumeric[i];
  }
  delete[] data;
  delete[] dataNumeric;

  for(int i=0; i<2; i++){
	  delete [] dataNumericIdx[i];
  }
  delete[] dataNumericIdx;

  delete[] looklog;
  delete[] lookH;

  for(int i=0; i<n+1; i++){
	  delete[] looklbc[i];
  }
  delete[] looklbc;

   delete [] AllLevels_red;
   delete [] cnt_red;
   delete [] posArray_red;

   for(int j = 0; (j < 2); j++){
   	delete [] dataNumericIdx_red[j];
   	delete [] dataNumeric_red[j];
   }

   delete [] dataNumericIdx_red;
   delete [] dataNumeric_red;

   delete[] samplesToEvaluateTamplate;
   delete[] samplesToEvaluate;

  if(samplesNotNA <= 2) delete[] res_new;

  return result;
}
