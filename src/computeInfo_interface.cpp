

//#define _MY_DEBUG_ 1
//#define _MY_DEBUG_NEW 1
#include <math.h>
#include <string.h>

#include <stdio.h>
#include <stdlib.h>
#include "memory.h"
#include "computeInfo.h"

#include <numeric>
#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include "utilities.h"
#include <set>
#include <math.h>

using namespace std;

#define LARGE 1E+300

#define LARGE 1E+300

//structures for threads
struct ContainerT {
	int** sortedSample;
	int* ptrAllData;
	int* ptrAllLevels;
	int sampleSizeEff;
	int zpos;
	int nSample0;
	int ptrzi;
	int sampleSize;
	int nbrAllVar;
	int nbrUi;
	int k23;
	int** Nxuiz;
	int* Ny;
	int* Nxui;
	int* Nx;
	int* Nxyuiz;
	int* Nyuiz;
	int* Nuiz;
	int* Nz;
	int* ptrZiIdx;
	int modCplx;
	double* c2terms;
	struct ScoresOfZ** scoresAllZi;
	double* weights;
	bool testDistribution;
	double** freqs1;

	int** sample;
	int* ptrVarIdx;
};

struct ScoresOfZ {
	int N_xyuiz; 
	double info_xy_ui;
	double logC_xy_ui;

	double info_xy_uiz;
	double logC_xy_uiz;
	double Rzi;
};

struct FromToPlusContainer {
	int from; 
	int to;
	struct ContainerT* cont;
	MemorySpace* memoryThread;
	int* Opt_dBin;
	int* dBin;
	int* nbrLevCorrect;
	double* Pxyuiz;
};



// modCplx = myCplx = 0 --> MDL, myCplx = 1 --> NML
// If nbrZi== 0, return nSample[0]     & nSample[0]*I(xy|{ui})      & k_xy_ui
// If nbrZi > 0, return nSample[z_top] & nSample[z_top]*I(xy|{ui})  & k_xy_ui
//                      z_top          & nSample[z_top]*I(xy|{ui}z) & k_xy_uiz
//                      R_top          & nSample[z_top]*I(xyz|{ui}) & k_xyz_ui




double* getAllInfoNEW( int* ptrAllData, int* ptrAllLevels, int* ptrVarIdx, int nbrUi, int* ptrZiIdx, int nbrZi, int ziPos, int sampleSize, 
	int sampleSizeEff, int modCplx, int k23, double* c2terms, MemorySpace* memory, double* weights, double** freqs1, bool testDistribution)
{	

	int randomrescaling=1;
	float r,rr;

	int bin_max=100,MDL=0, TRUE=1, FALSE=0;
	int l,ok;
	int **sample,**sortedSample,**Opt_sortedSample; //[N+1][7]

	int iii;
	int nrow=sampleSize+1;
	int ncol=7;


	sample = (*memory).sample;

	int * sampleWithZ = new int[sampleSize];

	sortedSample = (*memory).sortedSample;

	Opt_sortedSample = (*memory).Opt_sortedSample;


	// sample       = imatrix(0,sampleSize,0,6);
	// sortedSample = imatrix(0,sampleSize,0,6);

	int nSample0,*nSample,*orderSample,*sampleKey;//[N+1 elements: 0 to N]

	//nSample = (int *)malloc(sampleSize+1* sizeof(int));
	nSample = (int *)malloc(nbrZi* sizeof(int));
	orderSample = (*memory).orderSample;
	sampleKey = (*memory).sampleKey;

	//nSample = ivector(0,sampleSize);
	// orderSample = ivector(0,sampleSize+1);
	// sampleKey   = ivector(0,sampleSize+1);

	int bin,PBin,Prui,increment,X,Y,Z;
	int ptrzi,zi,z;

	double NlogN,logN;

	int Lxyui,Lyui,Lui;
	// int  Nxyui;
	double Pxyui;

	int   Nyui,  Nui;
	int  Nxyuis,  Nyuis,  Nuis;

	int NNxyui,NNxyuiz,NNxyuizl,Ntot;  //for rescaling NML change 20160228

	//int NxNui,NyNui,NNxui,NNyui; //for rescaling MDL change 20160221


	//int *Nxyuiz;
	int *Nyuiz, *Nuiz, *Nz;  //[Z]
	//Nxyuiz = (int *)malloc(bin_max+1* sizeof(int));
	//Nyuiz = (int *)malloc(bin_max+1* sizeof(int));
	//Nuiz = (int *)malloc(bin_max+1 * sizeof(int));
	//Nz = (int *)malloc(bin_max+1 * sizeof(int));
	// Nxyuiz = (*memory).Nxyuiz;
	double* Pxyuiz = (*memory).Pxyuiz;

	int* bridge = (*memory).bridge;

	Nyuiz = (*memory).Nyuiz;
	Nuiz = (*memory).Nuiz;
	Nz = (*memory).Nz;
	int Nzs,Nuizs,Nyuizs,Nxyuizs,Nxuizs;

	//Nxyuiz = ivector(0,bin_max);
	//Nyuiz  = ivector(0,bin_max);
	//Nuiz   = ivector(0,bin_max); 
	//Nz     = ivector(0,bin_max);
	// int  Nxyuizl;
	double Pxyuizl;
	int Nyuizl, Nuizl, Nzl; //[Z]

	int *Ny;     
	//Ny = (int *)malloc(bin_max+1* sizeof(int));
	Ny = (*memory).Ny;
	//Ny     = ivector(0,bin_max);
	int  Nyj,Nys;                     //[Y]

	int *Nxui, *Nx;                   //[X]
	//Nxui = (int *)malloc(bin_max+1* sizeof(int));
	//Nx = (int *)malloc(bin_max+1* sizeof(int));
	Nxui = (*memory).Nxui;
	Nx = (*memory).Nx;
	//Nxui   = ivector(0,bin_max);
	//Nx     = ivector(0,bin_max);
	int  Nxuij, Nxj, Nxuis, Nxs;      //[X]

	int **Nxuiz;                      //[X][Z]

	nrow = bin_max+1;
	ncol = bin_max+1;
	//Nxuiz = (int **)malloc(nrow* sizeof(int*));
	//for(iii = 0; iii < nrow; iii++) 
	//	Nxuiz[iii] = (int *)malloc(ncol* sizeof(int));
	Nxuiz = (*memory).Nxuiz;
	//Nxuiz=imatrix(0,bin_max,0,bin_max);
	int Nxuizjl;                     //[X][Z]

	double info_xui_y,info_yui_x,info_ui_y,info_ui_x;
	double logC_xui_y,logC_yui_x,logC_ui_y,logC_ui_x;

	double info_xuiz_y,info_yuiz_x,info_uiz_y,info_uiz_x;
	double logC_xuiz_y,logC_yuiz_x,logC_uiz_y,logC_uiz_x;

	double info_xui_z,info_yui_z,info_ui_z;
	double logC_xui_z,logC_yui_z,logC_ui_z;


	double info3xy_ui,info2xy_ui,info3xy_uiz,info2xy_uiz;
	double info3xz_ui,info2xz_ui,info3yz_ui,info2yz_ui;
	double logC3xy_ui,logC2xy_ui,logC3xy_uiz,logC2xy_uiz;
	double logC3xz_ui,logC2xz_ui,logC3yz_ui,logC2yz_ui;


	double info_xy_ui,info_xy_uiz;
	double info_xz_ui,info_yz_ui;
	double logC_xy_ui,logC_xy_uiz;
	double logC_xz_ui,logC_yz_ui;
	double testinfo_xy_ui,testinfo_xy_uiz;
	double testinfo_xz_ui,testinfo_yz_ui;
	double yz,xz,xyz,first,second,dpi,Rzi,testz;

	//int N_xy_ui=-1,N_xyuiz=-1;
	int N_xy_ui=0,N_xyuiz=0;
	//double LARGE DBL_MAX;
	double min_info_logC=LARGE, max_info_logC= - LARGE;
	//double min_info_logC = 100000.0, max_info_logC = -100000.0;

	int countmin,kz0;

	int z_top,N_xyuiz_top;
	double R_top;
	double NIxy_ui=-1.0,NIxy_ui_top,NIxy_uiz_top,NIxyz_ui_top;
	double k_xy_ui=-1.0,k_xy_ui_top,k_xy_uiz_top,k_xyz_ui_top;
  
  /////////////////////////

	int i, j, k;	// for loops

	int nbrAllVar;
	// int nbrSlot, maxNbrCombn, idxCount;
	// Output pointer
	int nbrRetValues = 3;

	// Local pointers
	//int *cpVarIdx = NULL;



// previous code:
// If no zi, return nSample0*I(..|{ui})[0] & nml[0] & nSample0
// If nbrZi > 0, return {nSample[z1]*I(..|{ui})[z1], NML(..|{ui})[z1], nSample[z1],nSample[z2]*I(..|{ui})[z2], NML(..|{ui})[z2], nSample[z2], ... } 
	//	if( nbrZi > 0 )	{ nbrRetValues = nbrRetValues * nbrZi; }	
// new code:
// If no zi, return nSample0 & nSample0*I(xy|{ui}) & k_xy_ui 
// If nbrZi > 0, return nSample[z_top] & nSample[z_top]*I(xy|{ui}) & k_xy_ui 
//                      z_top & nSample[z_top]*I(xy|{ui}z) & k_xy_uiz 
//                      R_top & nSample[z_top]*I(xyz|{ui}) & k_xyz_ui 
	if( nbrZi > 0 )	{ nbrRetValues = 9; }	

	double *ptrRetValues = (double*)malloc(nbrRetValues* sizeof(double));	

	ptrRetValues[0] = -1; 
	ptrRetValues[1] = -1;
	ptrRetValues[2] = -1;

	if( nbrZi > 0 ){
		ptrRetValues[3] = -1; 
		ptrRetValues[4] = -1;
		ptrRetValues[5] = -1;
		ptrRetValues[6] = -1; 
		ptrRetValues[7] = -1;
		ptrRetValues[8] = -1;
	}	

	// Define the total number of variables (xi, xj, {ui}, {zi})
	nbrAllVar = (nbrUi+2);
	//if( ziPos == nbrAllVar ) { nbrAllVar++; }; removed from original code

	int **dBin,**Opt_dBin; //[1][Nb_ui+2]

	nrow = 0 + 1;
	ncol = nbrAllVar +1;
	//dBin = (int **)malloc(nrow* sizeof(int*));
	//for(iii = 0; iii < nrow; iii++) 
	//	dBin[iii] = (int *)malloc(ncol* sizeof(int));
	dBin = (int **)malloc(nrow* sizeof(int*));
	for(iii = 0; iii < nrow; iii++) 
	  dBin[iii] = (int *)malloc(ncol* sizeof(int));

	int* nbrLevCorrect = (int *)malloc(ncol * sizeof(int));
		
	//Opt_dBin = (int **)malloc(nrow * sizeof(int*));
	//for(iii = 0; iii < nrow; iii++) 
	//	Opt_dBin[iii] = (int *)malloc(ncol * sizeof(int));
	Opt_dBin = (int **)malloc(nrow* sizeof(int*));
	for(iii = 0; iii < nrow; iii++) 
	  Opt_dBin[iii] = (int *)malloc(ncol* sizeof(int));
	//dBin     = imatrix(0,0,0,nbrAllVar);
	///Opt_dBin = imatrix(0,0,0,nbrAllVar);

	#if _MY_DEBUG_NEW
		printf("\n# ---- Printed from C [getAllInfo]: ----\n");
		
		printf("\n# ----> DATA\n");
		// i --> sample number; j --> variable number
		printf( "X0=%d Y0=%d\t", ptrAllData[0+ptrVarIdx[0]*sampleSize], ptrAllData[0+ptrVarIdx[1]*sampleSize]);

		for( i=0; i < sampleSize; i++ )
		  {
		    if(ptrAllData[i+ptrVarIdx[0]*sampleSize]>-1 && ptrAllData[i+ptrVarIdx[1]*sampleSize]>-1)  printf( "X=%d Y=%d\t", ptrAllData[i+ptrVarIdx[0]*sampleSize], ptrAllData[i+ptrVarIdx[1]*sampleSize]);
		    i = sampleSize;
		  }

		for( i=0; i < 6; i++ )
		{
			// Display data for xi, xj, {ui}
			printf("# "); 
			for( j=0; j < (2+nbrUi); j++ ) { 
				printf( "%d\t", ptrAllData[i+ptrVarIdx[j]*sampleSize] );
			} 
			// Display data for zi
			for( k=0; k < nbrZi; k++ ) { 
				printf( "%d\t", ptrAllData[i+ptrZiIdx[k]*sampleSize] );
			} 
			printf("\n");
		}
		printf("# ../..\n"); 

		printf("\n# ----> LEVELS\n");
		printf("# ");	for( i=0; i < (2+nbrUi); i++ ) { printf( "%d\t", ptrAllLevels[ptrVarIdx[i]] ); } 
		for( i=0; i < nbrZi; i++ ) { printf( "%d\t", ptrAllLevels[ptrZiIdx[i]] ); }
		printf("\n");

		printf("\n# ----> VAR IDX (xi, xj, {ui})\n");
		printf("# ");	for( i=0; i < (2+nbrUi); i++ ) { printf( "%d\t", ptrVarIdx[i] ); } 
		printf("\n");

		printf("\n# ----> NBR UI\n");
		printf("# %d\n", nbrUi);

		printf("\n# ----> ZI IDX (nbrZi=%d)\n",nbrZi);
		printf("# ");	for( i=0; i < nbrZi; i++ ) { printf( "%d\t", ptrZiIdx[i] ); } 
		printf("\n");

		printf("\n# ----> ZI POS\n");
		printf("# %d\n", ziPos);

		printf("\n# ----> SAMPLE SIZE\n");
		printf("# %d\n", sampleSize);

		printf("\n# ----> EFFECTIVE SAMPLE SIZE\n");
		printf("# %d\n", sampleSizeEff);

		printf("\n# ----> NBR ALLVAR\n");
		printf("# %d\n\n", nbrAllVar);

		printf("\n# ----> Check for MODEL COMPLEXITY\n");
		printf("# %d\n\n", modCplx);
	#endif // _MY_DEBUG_NEW



	// find samples without NA in x,y,ui and store their id in sample[k][0]
	for( i = 0, k = 0; i < sampleSize; i++)
	{
		ok=TRUE;
		for( j = 0; j < nbrAllVar; j++ )
		{
			// check that X,Y,Uj do not contain NA
			if( ptrAllData[ i + ptrVarIdx[j]*sampleSize ] == -1 )
			{
				ok=FALSE; 
				j=nbrAllVar;
			} 
		}  



		if(nbrZi==1) // if only one variable zi (to estimate NI(xy|ui) and 
		{          // NI(xy|uiz) on the exact same samples in case of NA)
			ptrzi=ptrZiIdx[0];
			if( ptrAllData[ i + ptrzi*sampleSize ] == -1 )
			{
				ok=FALSE; 		
			}
		}

		if(ok==TRUE)
		{
			sample[k][0] = i; // sample number
			k++;
		}
	}
    nSample0=k;



    //update dbin
    if(nSample0 < sampleSize){

    	for( j = 0; j < nbrAllVar; j++ ){
    		std::set<int> s;

		   	for( k = 0; k < nSample0; k++ )
			{
				i = sample[k][0];
				s.insert(ptrAllData[ i + ptrVarIdx[j]*sampleSize ]);
			}
			nbrLevCorrect[j] = s.size();
		}
	} else {
		for( j = 0; j < nbrAllVar; j++ )
			nbrLevCorrect[j]=ptrAllLevels[ptrVarIdx[j]];
	}


	#if _MY_DEBUG_NEW
		printf("\n# =====> test before initialisation of bin numbers for continuous variables \n");
		printf("# ==> sampleSize=%d    nSample0=%d  \n",sampleSize,nSample0);
	#endif // _MY_DEBUG_NEW



	if(nSample0 > 0){

		// initialisation of bin numbers for continuous variables in x,y, ui
		for(j=0; j < nbrAllVar; j++ )
		{

			dBin[0][j]=ptrAllLevels[ptrVarIdx[j]];

			Opt_dBin[0][j] = dBin[0][j];
		}

		#if _MY_DEBUG_NEW
			printf("\n# =====> test before optimisation on xy_ui \n");
			for(j=0; j < nbrAllVar; j++ ){
				printf("# dBin[0][%d] = %d \n",j,dBin[0][j]);
			}
			/*    for(j=0; j < nSample0; j++ )
			  {
				printf("# sample[%d][0] = %d \n",j,sample[j][0]);
			  }*/

		#endif // _MY_DEBUG_NEW

		countmin=0; // change 20160216

		//compute Lxyui, Lyui, Lui indices for later counting purpose
		for( k = 0; k < nSample0; k++ )
		{ 
			i=sample[k][0]; // sample number

			bin=ptrAllData[ i + ptrVarIdx[0]*sampleSize ];


			sample[k][1] = bin; // Lxyui initialisation
			sample[k][4] = bin; // X

			bin=ptrAllData[ i + ptrVarIdx[1]*sampleSize ];

			sample[k][5] = bin; // Y
			PBin=dBin[0][0];		
			increment = bin*PBin;
			sample[k][1] += increment; // Lxyui
			sample[k][2] = increment; // Lyui initialisation
			sample[k][3] = 0; // Lui initialisation
	  
			for( j = 2; j < nbrAllVar; j++ )
			{
				bin=ptrAllData[ i + ptrVarIdx[j]*sampleSize ];


				PBin *= dBin[0][j-1];
				increment = bin*PBin;
				sample[k][1] += increment; // Lxyui
				sample[k][2] += increment; // Lyui
				sample[k][3] += increment; // Lui
			} 	  
		}


		//k++;
		bin = PBin*dBin[0][nbrAllVar-1];

		sample[k][0] = nSample0; // extra sample id (useful for termination of counts)
		sample[k][1] = bin; // max Lxyui (useful for termination of counts)
		sample[k][2] = bin; // max Lyui  (useful for termination of counts)
		sample[k][3] = bin; // max Lui   (useful for termination of counts)
		/*	
		sample[k][4] = 0; // X
		sample[k][5] = 0; // Y
		sample[k][6] = 0; // Z 
		*/


//////////////////////////////////////////////////
		#if _MY_DEBUG_NEW
			printf("\n# =====> test before sorting \n");
		#endif // _MY_DEBUG_NEW

		/////////sort sample in increasing Lxyui stored in sample[k][1]
		//for( k = 0; k <= nSample0; k++ ){ 
		for( k = 1; k <=nSample0+1; k++ ){ 
			orderSample[k] = k-1;
			sampleKey[k]=sample[k-1][1]; // will sort in increasing Lxyui
		}

		#if _MY_DEBUG_NEW
			printf("\n# =====> test before sort2int \n");
		#endif // _MY_DEBUG_NEW

		sort2arrays(nSample0+1,sampleKey,orderSample, bridge);

		#if _MY_DEBUG_NEW
			printf("\n# =====> test after sort2int \n");
			for( k = 0; k <= nSample0; k++ ){ 
			  // if(sample[k][0]<0) 
			  //if(k>(nSample0-10))  printf("@ nSample0=%d   sampleKey[k=%d]=%d  orderSample[k=%d]=%d  sample[i=%d][1]=%d  sample[i=%d][0]=%d \n",nSample0,k,sampleKey[k],k,orderSample[k],k,sample[k][1],k,sample[k][0]);
			}
		#endif // _MY_DEBUG_NEW


		//for( k = 0; k <= nSample0; k++ ){ 
		for( k = 1; k <= nSample0+1; k++ ){ 

			#if _MY_DEBUG_NEW
				  //	  printf("\n k=%d ",k);
			#endif // _MY_DEBUG_NEW

				  i=orderSample[k];

			#if _MY_DEBUG_NEW
				  if(i<0 || i>sampleSize) printf("\n# @@@@@@@@@@@@@@@@@@@@> Probs i=%d sampleSize=%d \n",i,sampleSize);
				  if(k>(nSample0-10))  printf("\n @ ordS[k=%d]=%d ",k,orderSample[k]);
			#endif // _MY_DEBUG_NEW

			for( j=0; j< 6; j++ ){

			    sortedSample[k-1][j]=sample[i][j];
				#if _MY_DEBUG_NEW
					    if(k>(nSample0-10))  printf("sdS[k=%d][j=%d]=%d ",k-1,j,sortedSample[k-1][j]);
				#endif // _MY_DEBUG_NEW
			}
		}
//////////////////////////////////////////////////

		#if _MY_DEBUG_NEW
			printf("\n# =====> test before initialization of counts and mutual infos & logCs \n");
		#endif // _MY_DEBUG_NEW


		/////////initialization of counts and mutual infos & logCs 
		Lxyui = sortedSample[0][1]; // min Lxyui
		Lyui  = sortedSample[0][2]; // min Lyui
		Lui   = sortedSample[0][3]; // min Lui

		//Nxyui = 1;
		Pxyui = weights[sortedSample[0][0]];
		NNxyui= 0;
		Nyui  = 0;
		Nui   = 0;

		for( k = 0; k < dBin[0][0]; k++ ){ 
		  Nxui[k]=0;
		  Nx[k]=0;
		}

		for( k = 0; k < dBin[0][1]; k++ ) Ny[k]=0;
		X=sortedSample[0][4];
		Y=sortedSample[0][5];

		info_xui_y = 0.0;
		info_yui_x = 0.0;
		info_ui_y  = 0.0;
		info_ui_x  = 0.0;
		logC_xui_y = 0.0;
		logC_yui_x = 0.0;
		logC_ui_y  = 0.0;
		logC_ui_x  = 0.0;

		Nxyuis=0; // 6 test variables for counts, should all sum to nSample0
		Nyuis=0;
		Nxuis=0;
		Nuis=0;
		Nxs=0;
		Nys=0;
		Ntot=0;

		/////////make the counts and compute mutual infos & logCs 
		for( k = 1; k <= nSample0; k++ ){ 


			#if _MY_DEBUG_NEW
				if(k>(nSample0-10) || Nxyui<0 || Nyui<0 || Nui<0 || X<0 || Y<0 || X>=dBin[0][0] || Y>=dBin[0][1]){
				  printf("before counts nSam0=%d  X=%d(%d) Y=%d(%d) k=%d sortS[%d][1]=%d Lxyui=%d sortS[%d][2]=%d Lyui=%d sortS[%d][3]=%d Lui=%d Nxyui=%d Nyui=%d Nui=%d  \n",nSample0,X,dBin[0][0],Y,dBin[0][1],k,k,sortedSample[k][1],Lxyui,k,sortedSample[k][2],Lyui,k,sortedSample[k][3],Lui,Nxyui,Nyui,Nui);
				}
			#endif // _MY_DEBUG_NEW

			if(sortedSample[k][1] > Lxyui){

			        NNxyui=0;
			        if(sampleSizeEff!=sampleSize){
				  if(randomrescaling){
				    NNxyui=(int)floor(Pxyui);
				    r=Pxyui - NNxyui;
				    rr=(double)rand() / (double)RAND_MAX ;
				    // printf("%f\n",rr); 
				    if (r > rr) NNxyui++;
				    //printf("# =====> rando  r=%g rr=%g NNxyui=%d \n",r,rr,NNxyui);
				  }
				  // else{
				  //   NNxyui=(int)floor(0.5+((double)Nxyui*sampleSizeEff)/sampleSize);
				  // }
				}else{
				  //NNxyui=Nxyui; 
				  NNxyui=(int)Pxyui; 
				}
				
				Nui     += NNxyui; 
				Nyui    += NNxyui; 
				Nxui[X] += NNxyui; 
				
				Nx[X] += NNxyui;	//4
				Ny[Y] += NNxyui;	//4
				Ntot +=  NNxyui;

				if( NNxyui > 0 ){
				  NlogN = NNxyui * log(NNxyui);
				  info_xui_y += NlogN;
				  info_yui_x += NlogN;
				}
				Lxyui = sortedSample[k][1];
				Nxyuis += NNxyui;
				Pxyui = weights[sortedSample[k][0]]; //weights[k];	
				//Nxyui = 1;

				if(k<nSample0) X=sortedSample[k][4];
				if(k<nSample0) Y=sortedSample[k][5];
		


				if(sortedSample[k][2] > Lyui){

				        if( Nyui > 0 ){
					  NlogN = Nyui * log(Nyui);
					  info_yui_x -= NlogN;
					  info_ui_y  += NlogN;
					  
					  if(modCplx != MDL){
					    if(Nyui>Ntot) printf("# ==$$$===> Ntot=%d Nyui=%d \n",Ntot,Nyui);
					    logC_yui_x += compute_LogC_C2(Nyui,dBin[0][0],c2terms);
					  }
					}
					Lyui = sortedSample[k][2];
					Nyuis += Nyui;
					Nyui = 0;

					if(sortedSample[k][3] > Lui){ 
						for( j = 0; j < dBin[0][0]; j++ ){ 
							Nxuij=Nxui[j];
							if( Nxuij > 0 ){
								NlogN = Nxuij * log(Nxuij);
								info_xui_y -= NlogN;
								info_ui_x  += NlogN;
								if(modCplx != MDL){
								    if(Nxuij>Ntot) printf("# ==$$$===> Ntot=%d Nxuij=%d \n",Ntot,Nxuij);
								    logC_xui_y += compute_LogC_C2(Nxuij,dBin[0][1],c2terms);
								}
				   
								Nxuis += Nxuij;
								Nxui[j]=0;
							}
						}
						
						if( Nui > 0 ){
						  NlogN = Nui * log(Nui);
						  info_ui_y  -= NlogN;
						  info_ui_x  -= NlogN;
						  if(modCplx != MDL){
						    if(Nui>Ntot) printf("# ==$$$===> Ntot=%d Nui=%d \n",Ntot,Nui);
						    logC_ui_x += compute_LogC_C2(Nui,dBin[0][0],c2terms);
						    logC_ui_y += compute_LogC_C2(Nui,dBin[0][1],c2terms);
						  }
						}
						Lui = sortedSample[k][3];
						Nuis += Nui;
						Nui = 0;

					}
				}

			} else {
			  //Nxyui++;
			  Pxyui += weights[sortedSample[k][0]]; //weights[k];	
			}

		}


		// increment for info for Nx[X] and Ny[Y] contributions	
		for( j = 0; j < dBin[0][0]; j++ ){ 
			Nxj=Nx[j];
			if( Nxj > 0 ){
			      //NlogN = Nxj * log(Nxj/(1.0*nSample0));
				NlogN = Nxj * log(Nxj/(1.0*Ntot));
				info_yui_x -= NlogN;
				info_ui_x  -= NlogN;
				Nxs += Nxj;
			}
		}
		for( j = 0; j < dBin[0][1]; j++ ){ 
			Nyj=Ny[j];
			if( Nyj > 0 ){
			      //NlogN = Nyj * log(Nyj/(1.0*nSample0));
				NlogN = Nyj * log(Nyj/(1.0*Ntot));
				info_xui_y -= NlogN;
				info_ui_y  -= NlogN;
				Nys += Nyj;
			}
		}
    

		/////////check maximum mutual infos - cplx terms
		
		//info3xy_ui = 0.5*(info_xui_y + info_yui_x)*sampleSizeEff/(1.0*sampleSize) ;
		//info2xy_ui = 0.5*(info_ui_y  + info_ui_x)*sampleSizeEff/(1.0*sampleSize)  ;
		info3xy_ui = 0.5*(info_xui_y + info_yui_x);
		info2xy_ui = 0.5*(info_ui_y  + info_ui_x);


		#if _MY_DEBUG_NEW
		    printf("# =====> test before check maximum mutual infos - cplx terms nSample0=%d Nxyuis=%d   Nyuis=%d   Nuis=%d  Nxuis=%d  Nxs=%d  Nys=%d  Ntot=%d  info_ui_y=%g info_ui_x=%g info_ui_y+info_ui_x =%g  \n",nSample0,Nxyuis,Nyuis,Nuis,Nxuis,Nxs,Nys,Ntot,info_ui_y,info_ui_x,info_ui_y+info_ui_x);
			printf("# =====> test before check maximum mutual infos - cplx terms sampleSizeEff=%d sampleSize=%d info3xy_ui=%g   info2xy_ui=%g \n",sampleSizeEff,sampleSize,info3xy_ui,info2xy_ui);
			printf("# =====> test before check maximum mutual infos - cplx terms info_xui_y=%g info_yui_x=%g info_ui_y=%g info_ui_x=%g \n",info_xui_y,info_yui_x,info_ui_y,info_ui_x);
			printf("# =====> test before check maximum mutual infos - cplx terms logC_xui_y=%g logC_yui_x=%g logC_ui_y=%g logC_ui_x=%g \n",logC_xui_y,logC_yui_x,logC_ui_y,logC_ui_x);
		#endif // _MY_DEBUG_NEW


		if(modCplx == MDL) {
			Prui=1;
			//NN=nSample0;
			//NN = (int)floor(0.5 + ((double)nSample0*sampleSizeEff)/sampleSize);
			//logN=log(NN);
			logN=log(Ntot);
			for(j=2; j < nbrAllVar; j++ ) Prui *= dBin[0][j];
			logC_xui_y= 0.5*(dBin[0][1] - 1)*(dBin[0][0]*Prui - 1)*logN;
			logC_yui_x= 0.5*(dBin[0][0] - 1)*(dBin[0][1]*Prui - 1)*logN;
			logC_ui_y = 0.5*(dBin[0][1] - 1)*(Prui - 1)*logN;
			logC_ui_x = 0.5*(dBin[0][0] - 1)*(Prui - 1)*logN;
		}

		logC3xy_ui = 0.5*(logC_xui_y + logC_yui_x) ;
		logC2xy_ui = 0.5*(logC_ui_y  + logC_ui_x)  ;
	 
		if(nbrUi==0) testinfo_xy_ui =info3xy_ui - info2xy_ui - logC3xy_ui + logC2xy_ui; //change 20160221
		else         testinfo_xy_ui =info3xy_ui - logC3xy_ui + info2xy_ui - logC2xy_ui;

		#if _MY_DEBUG_NEW
			printf("\n# =====> test before change bin numbers 1,  info3xy_ui - logC3xy_ui=%g  info2xy_ui - logC2xy_ui=%g testinfo_xy_ui=%g info3xy_ui=%g  logC3xy_ui=%g  info2xy_ui=%g logC2xy_ui=%g max_info_logC=%g min_info_logC=%g \n",(info3xy_ui - logC3xy_ui),(info2xy_ui - logC2xy_ui),testinfo_xy_ui,info3xy_ui,logC3xy_ui,info2xy_ui,logC2xy_ui,max_info_logC,min_info_logC);
		#endif // _MY_DEBUG_NEW

		if(max_info_logC < testinfo_xy_ui){
			//N_xy_ui =nSample0;
   			        //N_xy_ui = (int)floor(0.5 + ((double)nSample0*sampleSizeEff)/sampleSize);
			N_xy_ui = Ntot;
			NIxy_ui = info3xy_ui - info2xy_ui; // info to be returned if no z
			k_xy_ui = logC3xy_ui - logC2xy_ui; // cplx to be returned if no z

			max_info_logC = testinfo_xy_ui;
			min_info_logC = max_info_logC;
			/////////////////////////////////////
			for(j=0; j < nbrAllVar; j++ ) Opt_dBin[0][j] = dBin[0][j];
			if(nbrZi>0) {
			  for( k = 0; k <= nSample0; k++ ){ 
			    for( j=0; j< 6; j++ ){
			      Opt_sortedSample[k][j]=sortedSample[k][j];
			    }
			  }
			}
			////////////////////////////////////

		} else if(min_info_logC > testinfo_xy_ui){
			countmin++;
			min_info_logC = testinfo_xy_ui;
		} else {
			countmin=0;
			min_info_logC = testinfo_xy_ui;
		}

		#if _MY_DEBUG_NEW
			printf("\n# =====> test before change bin numbers 2,  info3xy_ui - logC3xy_ui=%g  info2xy_ui - logC2xy_ui=%g testinfo_xy_ui=%g info3xy_ui=%g  logC3xy_ui=%g  info2xy_ui=%g logC2xy_ui=%g max_info_logC=%g min_info_logC=%g \n",(info3xy_ui - logC3xy_ui),(info2xy_ui - logC2xy_ui),testinfo_xy_ui,info3xy_ui,logC3xy_ui,info2xy_ui,logC2xy_ui,max_info_logC,min_info_logC);
		#endif // _MY_DEBUG_NEW

	
		// for(p = 0; p < nbrAllVar; p++)
		//printf("vars: nbrAllVar=%d nbrZi=%d max_info_logC=%g min_info_logC=%g", nbrAllVar, nbrZi , max_info_logC,min_info_logC); 
		//printf("\n"); 

		// 	for( j=0; j < nbrAllVar; j++ ){
		//   printf(" Opt_dBin[0][j=%d]=%d",j,Opt_dBin[0][j]);
		// }
		// printf("\n");

	#if _MY_DEBUG_NEW
	    printf("\n# =====> test before z nbrZi=%d \n",nbrZi);
	#endif // _MY_DEBUG_NEW

 
///////////////////////////////////////////////////////////////////////
    if(nbrZi==0){

		ptrRetValues[0] = N_xy_ui; 
		ptrRetValues[1] = NIxy_ui;
		ptrRetValues[2] = k_xy_ui;

		// cout <<  NIxy_ui << "\t" << nbrLevCorrect[0] << ", " << nbrLevCorrect[1] << ", " << N_xy_ui << endl;

		if(nbrLevCorrect[0] == 1 || nbrLevCorrect[1] == 1 || nbrLevCorrect[0] == N_xy_ui || nbrLevCorrect[1] == N_xy_ui ){
			ptrRetValues[0] = N_xy_ui; 
			ptrRetValues[1] = 0;
			ptrRetValues[2] = k_xy_ui;
		}

/////////pick next z and compute score ///////////////////////////////
    } else { //(nbrZi>0)

	      /// re-establish optimum x,y,ui bins and Lxyui, Lyui, Lui
		// printf("££££££££££££££££££ re-establish optimum x,y,ui bins:");

		/////////////////////////////////////
	 	for( j=0; j < nbrAllVar; j++ ){
			dBin[0][j] = Opt_dBin[0][j];
			// printf(" dBin[0][j=%d]=%d",j,dBin[0][j]);
		}
		// printf("\n");
		for( k = 0; k <= nSample0; k++ ){ 
		  for( j=0; j< 6; j++ ){
		    sortedSample[k][j]=Opt_sortedSample[k][j];
		  }
		}
		////////////////////////////////////

		/// find optimum zi /////////////////////////////////////
		z=nbrAllVar;

		//N_xyuiz_top = -1; 
		N_xyuiz_top = 0; 
		NIxy_ui_top  = -1;
		k_xy_ui_top  = -1;
		
		z_top = -1;
		NIxy_uiz_top = -1;
		k_xy_uiz_top = -1;
		
		R_top=-LARGE;
		NIxyz_ui_top = -1;
		k_xyz_ui_top = -1;

		for(zi=0; zi < nbrZi ; zi++ ){

			// initialisation of bin numbers for continuous variable zi
			ptrzi=ptrZiIdx[zi];

			dBin[0][z]=ptrAllLevels[ptrzi];
			nSample[zi]=0;

			for( k = 0; k < nSample0; k++ ){ 
				i=sortedSample[k][0];
				// find the first sample for which zi does not contain NA
				if( ptrAllData[ i + ptrzi*sampleSize ] > -1 )
				{
					kz0=k;
					nSample[zi]=1;
					k=nSample0;
				}
			}

			for( i = 0, k = 0; i < sampleSize; i++)
			{
				// i = sample[k][0];
				ok=TRUE;
				for( j = 0; j < nbrAllVar; j++ )
				{
					// check that X,Y,Uj do not contain NA
					if( ptrAllData[ i + ptrVarIdx[j]*sampleSize ] == -1 )
					{
						ok=FALSE; 
						j=nbrAllVar;
					} 
				}  

				ptrzi=ptrZiIdx[zi];
				if( ptrAllData[ i + ptrzi*sampleSize ] == -1 )
				{
					ok=FALSE; 		
				}
				

				if(ok==TRUE)
				{
					sampleWithZ[k] = i; // sample number
					k++;
				}
			}

			int nSampleZ = k;

			bool isGoodCandidate = true;
			//update dbin
			// cout << "HERE\n\n\n\n";

		    if(nSampleZ < sampleSize){

		    	for( j = 0; j < 2 && isGoodCandidate; j++ ){
		    		std::set<int> s;
			    	// double* propObserved = (double *)calloc(ptrAllLevels[j], sizeof(int));



				   	for( k = 0; k < nSampleZ; k++ )
					{
						i = sampleWithZ[k];
						// propObserved[ptrAllData[ i + ptrVarIdx[j]*sampleSize ]]++;
						s.insert(ptrAllData[ i + ptrVarIdx[j]*sampleSize ]);
					}
					nbrLevCorrect[j] = s.size();

					// chi2 test for distribution
					// if(!testDistribution(propObserved, proportions[ptrVarIdx[j]], ptrAllLevels[j])){
					// 	// the sample is not representative of the population, hence we do not want this z as possible z
					// 	isGoodCandidate = false;
					// 	break;
					// }
					// free(propObserved);

					// KL distance 
				}

				if(testDistribution){
					int** counts2 = (int **)calloc(ptrAllLevels[ptrVarIdx[0]], sizeof(int*));
			    	for( j = 0; j < ptrAllLevels[ptrVarIdx[0]]; j++ )
			    		counts2[j] = (int *)calloc(ptrAllLevels[ptrVarIdx[1]], sizeof(int));

			    	// fill table
					for(int k = 0; k < nSampleZ; k++){
						i = sampleWithZ[k];
						counts2[ ptrAllData[i + ptrVarIdx[0]*sampleSize]][ptrAllData[ i + ptrVarIdx[1]*sampleSize]]++;
					}

					// double cplxMdl = (nbrUi + 3 ) * 5;
					double cplxMdl = log(nSampleZ);
					// double cplxMdl = ptrAllLevels[ptrVarIdx[0]] * ptrAllLevels[ptrVarIdx[1]] * log(nSampleZ) ; //*(ptrAllLevels[ptrzi] - 1) ;
					// for(int counterU = 0 ; counterU < nbrUi; counterU++)
					// 	cplxMdl *= ptrAllLevels[ptrVarIdx[counterU+2]];

					double kldiv = exp(-nSampleZ * kl(counts2, freqs1, ptrAllLevels[ptrVarIdx[0]], ptrAllLevels[ptrVarIdx[1]]) + cplxMdl);
					
					
					// cout << kldiv << endl;

					if(kldiv < 1)
						isGoodCandidate = false;

					for( j = 0; j < ptrAllLevels[ptrVarIdx[0]]; j++ )
						delete [] counts2[j];
					delete [] counts2;
				}
			} 

			if(isGoodCandidate){



				#if _MY_DEBUG_NEW
					  printf("\n# =====> test before z 2 nSample[zi=%d]=%d nSample0=%d \n",zi,nSample[zi],nSample0);
				#endif // _MY_DEBUG_NEW

				if(nSample[zi] == 1){

					max_info_logC = - LARGE;
					min_info_logC = LARGE;

					countmin=0; // change 20160216

				    nSample[zi]=0; 
					Lxyui = sortedSample[kz0][1]; // min Lxyui
					Lyui  = sortedSample[kz0][2]; // min Lyui
					Lui   = sortedSample[kz0][3]; // min Lui

					NNxyuiz = 0;
					NNxyuizl= 0;
					// Nxyui = 0;
					Pxyui = 0;
					Nyui  = 0;
					Nui   = 0;
					for( k = 0; k < dBin[0][0]; k++ ){ 
						Nxui[k]=0;
						Nx[k]=0;
					
						for( l = 0; l < dBin[0][z]; l++ ){
							Nxuiz[k][l]=0;
						}
					}

					for( k = 0; k < dBin[0][1]; k++ ){
						Ny[k]=0;
					}
					for( l = 0; l < dBin[0][z]; l++ ){
						//Nxyuiz[l]=0;
						Pxyuiz[l]=0;
						Nyuiz[l]=0;
						Nuiz[l]=0;
						Nz[l]=0;
					}

					X=sortedSample[kz0][4];
					Y=sortedSample[kz0][5];
					i=sortedSample[kz0][0];
					Z=ptrAllData[ i + ptrzi*sampleSize ];  // first Z

					//Nxyuiz[Z]=1;
					Pxyuiz[Z]=weights[sortedSample[kz0][0]];


					sortedSample[nSample0][0]=i; // to terminate loop properly below

					info_xui_y = 0.0;
					info_yui_x = 0.0;
					info_ui_y  = 0.0;
					info_ui_x  = 0.0;
					logC_xui_y = 0.0;
					logC_yui_x = 0.0;
					logC_ui_y  = 0.0;
					logC_ui_x  = 0.0;

					info_xuiz_y = 0.0;
					info_yuiz_x = 0.0;
					info_uiz_y  = 0.0;
					info_uiz_x  = 0.0;
					logC_xuiz_y = 0.0;
					logC_yuiz_x = 0.0;
					logC_uiz_y  = 0.0;
					logC_uiz_x  = 0.0;

					info_xui_z = 0.0;
					info_yui_z = 0.0;
					info_ui_z  = 0.0;
					logC_xui_z = 0.0;
					logC_yui_z = 0.0;
					logC_ui_z  = 0.0;

					info_xy_ui  = 0.0;
					info_xy_uiz = 0.0;
					info_xz_ui  = 0.0;
					info_yz_ui  = 0.0;
					logC_xy_ui  = 0.0;
					logC_xy_uiz = 0.0;
					logC_xz_ui  = 0.0;
					logC_yz_ui  = 0.0;

					Nxyuis=0; // 11 test variables for counts, should all sum to nSample0
					Nyuis=0;
					Nxuis=0;
					Nuis=0;
					Nxs=0;
					Nys=0;

					Nzs=0;
					Nuizs=0;
					Nyuizs=0;
					Nxyuizs=0;
					Nxuizs=0;

					#if _MY_DEBUG_NEW
			    	    //printf("\n# =====> test before z 4 nSample0=%d sortedSample[kz0=%d][0]=%d sortedSample[kz0+1=%d][0]=%d sortedSample[kz0=%d][1]=%d sortedSample[kz0=%d][2]=%d sortedSample[kz0=%d][3]=%d sortedSample[kz0=%d][4]=%d sortedSample[kz0=%d][5]=%d sortedSample[kz0=%d][6]=%d \n",nSample0,kz0,sortedSample[kz0][0],kz0+1,sortedSample[kz0+1][0],kz0,sortedSample[kz0][1],kz0,sortedSample[kz0][2],kz0,sortedSample[kz0][3],kz0,sortedSample[kz0][4],kz0,sortedSample[kz0][5],kz0,sortedSample[kz0][6]);
					#endif // _MY_DEBUG_NEW


					/////////make the counts and compute mutual infos & logCs 
					for( k = kz0+1; k <= nSample0; k++ ){ // change 20160220

						#if _MY_DEBUG_NEW
							//printf("z nSample0=%d  k=%d ptrAllData[ %d + ptrVarIdx[0]*sampleSize]=%d sortedSample[%d][0]=%d  sortedSample[%d][1]=%d Lxyui=%d sortedSample[%d][2]=%d Lyui=%d sortedSample[%d][3]=%d Lui=%d Nxyui=%d Nyui=%d Nui=%d  \n",nSample0,k,sortedSample[k][0],ptrAllData[ sortedSample[k][0] + ptrVarIdx[0]*sampleSize],k,sortedSample[k][0],k,sortedSample[k][1],Lxyui,k,sortedSample[k][2],Lyui,k,sortedSample[k][3],Lui,Nxyui,Nyui,Nui);
						#endif // _MY_DEBUG_NEW

						i=sortedSample[k][0];

						// check whether zi does not contain NA
						if( ptrAllData[ i + ptrzi*sampleSize ] > -1 )
						{

							Z=ptrAllData[ i + ptrzi*sampleSize ];  // Z

							if(sortedSample[k][1] > Lxyui){

							        NNxyuiz=0;
								for( l = 0; l < dBin[0][z]; l++ ){
									//Nxyuizl=Nxyuiz[l];//
									Pxyuizl=Pxyuiz[l];//
									if( Pxyuizl > 0 ){
									    if(sampleSizeEff!=sampleSize){
									        if(randomrescaling){
										    NNxyuizl=(int)floor(Pxyuizl);
										    r=Pxyuizl - NNxyuiz;
										    rr=(double)rand() / (double)RAND_MAX ;
										    // printf("%f\n",rr); 
										    if (r > rr) NNxyuizl++;
										    //printf("# =====> rando  r=%g rr=%g NNxyui=%d \n",r,rr,NNxyui);
										}
										// else{
										//     NNxyuizl=(int)floor(0.5+((double)Nxyuizl*sampleSizeEff)/sampleSize);
										// }
									    }else{
									        NNxyuizl=(int)Pxyuizl;
									    }

									    if( NNxyuizl > 0 ){
									        NlogN = NNxyuizl * log(NNxyuizl);
										info_xuiz_y += NlogN;//
										info_yuiz_x += NlogN;//
										NNxyuiz += NNxyuizl;

										Nz[l]   += NNxyuizl;	//4
										Nuiz[l] += NNxyuizl;
										Nyuiz[l]+= NNxyuizl;
										Nxuiz[X][l] += NNxyuizl;
									    }

									    //Nxyuiz[l]=0;//
									    Pxyuiz[l]=0;//
									}
								}
								//Nxyuiz[Z]=1;  //1  change
								// printf("%f ",weights[sortedSample[k][0]]);

								Pxyuiz[Z]=weights[sortedSample[k][0]];

								if(NNxyuiz > 0){
								    NlogN = NNxyuiz * log(NNxyuiz);
								    info_xui_y += NlogN;//
								    info_yui_x += NlogN;//

								    nSample[zi] += NNxyuiz;
							
								    Nxyuizs += NNxyuiz;

								    Nui     += NNxyuiz; 
								    Nyui    += NNxyuiz; 
								    Nxui[X] += NNxyuiz; 
									  
								    Nx[X]   += NNxyuiz;	//4
								    Ny[Y]   += NNxyuiz;	//4
								    Ntot    += NNxyuiz;
									  
								    Nxyuis  += NNxyuiz;
								}


								Lxyui = sortedSample[k][1];

								if(k<nSample0) X=sortedSample[k][4];
								if(k<nSample0) Y=sortedSample[k][5];


									if(sortedSample[k][2] > Lyui){
								    if( Nyui > 0 ){
								        NlogN = Nyui * log(Nyui);
									info_yui_x -= NlogN;//
									info_yui_z -= NlogN;//
									info_ui_y  += NlogN;//
									if(modCplx != MDL){
									    logC_yui_x += compute_LogC_C2(Nyui,dBin[0][0],c2terms);
									    logC_yui_z += compute_LogC_C2(Nyui,dBin[0][z],c2terms);
									}

									for( l = 0; l < dBin[0][z]; l++ ){
									    Nyuizl=Nyuiz[l];
									    if( Nyuizl > 0 ){
									        NlogN = Nyuizl * log(Nyuizl);
										info_yui_z += NlogN;//
										info_uiz_y += NlogN;//
										info_yuiz_x -= NlogN;//
										if(modCplx != MDL){
										    logC_yuiz_x += compute_LogC_C2(Nyuizl,dBin[0][0],c2terms);
										}
										Nyuizs += Nyuizl;
										Nyuiz[l]=0;//
									    }
									}
									//Nyuiz[Z]=1;  //2
								      
									Nyuis += Nyui;
									Nyui = 0; //2
								    }
								    Lyui = sortedSample[k][2];

								    if(sortedSample[k][3] > Lui){
								         if( Nui > 0 ){
									        NlogN = Nui * log(Nui);
										info_ui_x  -= NlogN;//
										info_ui_y  -= NlogN;//
										info_ui_z  -= NlogN;//
										if(modCplx != MDL){
											logC_ui_x += compute_LogC_C2(Nui,dBin[0][0], c2terms);//
											logC_ui_y += compute_LogC_C2(Nui,dBin[0][1], c2terms);//
											logC_ui_z += compute_LogC_C2(Nui,dBin[0][z], c2terms);//
										}
										Nuis += Nui;
										Nui = 0;   //3

										for( l = 0; l < dBin[0][z]; l++ ){
											Nuizl=Nuiz[l];
											if( Nuizl > 0 ){
												NlogN = Nuizl * log(Nuizl);
												info_ui_z  += NlogN;//
												info_uiz_x -= NlogN;//
												info_uiz_y -= NlogN;//
												if(modCplx != MDL){
												  logC_uiz_x += compute_LogC_C2(Nuizl,dBin[0][0], c2terms);//
												  logC_uiz_y += compute_LogC_C2(Nuizl,dBin[0][1], c2terms);//
												}
												Nuizs += Nuizl;
												Nuiz[l]=0;//
											}
										}
										//Nuiz[Z]=1;   //3


										for( j = 0; j < dBin[0][0]; j++ ){ 
											Nxuij=Nxui[j];
											if( Nxuij > 0 ){
												NlogN = Nxuij * log(Nxuij);
												info_xui_y -= NlogN;//
												info_xui_z -= NlogN;//
												info_ui_x  += NlogN;//
												if(modCplx != MDL){
												  logC_xui_y += compute_LogC_C2(Nxuij,dBin[0][1], c2terms);//
												  logC_xui_z += compute_LogC_C2(Nxuij,dBin[0][z], c2terms);//
												}
												Nxuis += Nxuij;
												Nxui[j]=0;//

												for( l = 0; l < dBin[0][z]; l++ ){
													Nxuizjl=Nxuiz[j][l];
													if( Nxuizjl > 0 ){
														NlogN = Nxuizjl * log(Nxuizjl);
														info_xui_z += NlogN;//
														info_uiz_x += NlogN;//
														info_xuiz_y -= NlogN;//
														if(modCplx != MDL){
														  logC_xuiz_y += compute_LogC_C2(Nxuizjl,dBin[0][1], c2terms);//
														}
														Nxuizs += Nxuizjl;
														Nxuiz[j][l]=0;//
												  	}
												}
												//Nxuiz[X][Z]=1;// change

											}
			       							}

								         }
								         Lui = sortedSample[k][3];

								    }

								}

							} else {
								// Nxyuiz[Z]++;	        //1
								Pxyuiz[Z]+=weights[sortedSample[k][0]];
							}
						}
					}


			    	// increment info with Nx[X], Ny[Y] and Nz[Z] contributions	
					for( j = 0; j < dBin[0][0]; j++ ){ 
						Nxj=Nx[j];
						if( Nxj > 0 ){
							NlogN = Nxj * log(Nxj/(1.0*nSample[zi]));
							info_yui_x  -= NlogN;//
							info_ui_x   -= NlogN;//
							info_uiz_x  -= NlogN;//
							info_yuiz_x -= NlogN;//
							Nxs += Nxj;
							Nx[j]=0;
						}
					}
					for( j = 0; j < dBin[0][1]; j++ ){ 
						Nyj=Ny[j];
						if( Nyj > 0 ){
							NlogN = Nyj * log(Nyj/(1.0*nSample[zi]));
							info_xui_y  -= NlogN;//
							info_ui_y   -= NlogN;//
							info_uiz_y  -= NlogN;//
							info_xuiz_y -= NlogN;//
							Nys += Nyj;
							Ny[j]=0;
						}
					}
					for( l = 0; l < dBin[0][z]; l++ ){ 
						Nzl=Nz[l];
							if( Nzl > 0 ){
							NlogN = Nzl * log(Nzl/(1.0*nSample[zi]));
							info_xui_z  -= NlogN;//
							info_yui_z  -= NlogN;//
							info_ui_z   -= NlogN;//
							Nzs += Nzl;
							Nz[l]=0;
						}
					}
			  
					//////////////////////////////////////
					#if _MY_DEBUG_NEW
						printf("# =====> Z test before check maximum mutual infos - cplx terms nSample0=%d nSample[zi=%d]=%d\n Nxyuis=%d   Nyuis=%d   Nuis=%d  Nxuis=%d  Nxs=%d  Nys=%d  Nzs=%d Nuizs=%d Nyuizs=%d Nxyuizs=%d Nxuizs=%d \n info_ui_y=%g info_ui_x=%g info_ui_y+info_ui_x =%g  \n",nSample0,zi,nSample[zi],Nxyuis,Nyuis,Nuis,Nxuis,Nxs,Nys,Nzs,Nuizs,Nyuizs,Nxyuizs,Nxuizs,info_ui_y,info_ui_x,info_ui_y+info_ui_x);
						printf("# =====> Z test before check maximum mutual infos - cplx terms sampleSizeEff=%d sampleSize=%d info3xy_ui=%g   info2xy_ui=%g \n",sampleSizeEff,sampleSize,info3xy_ui,info2xy_ui);
						printf("# =====> Z test before check maximum mutual infos - cplx terms info_xui_y=%g info_yui_x=%g info_ui_y=%g info_ui_x=%g \n",info_xui_y,info_yui_x,info_ui_y,info_ui_x);
						printf("# =====> Z test before check maximum mutual infos - cplx terms logC_xui_y=%g logC_yui_x=%g logC_ui_y=%g logC_ui_x=%g \n",logC_xui_y,logC_yui_x,logC_ui_y,logC_ui_x);
					#endif // _MY_DEBUG_NEW


					#if _MY_DEBUG_NEW

					if (info_xui_y < -0.000001 || info_yui_x < -0.000001 || info_ui_y  < -0.000001 || info_ui_x  < -0.000001 || info_xuiz_y < -0.000001 || info_yuiz_x < -0.000001 || info_uiz_y  < -0.000001 || info_uiz_x  < -0.000001 || info_xui_z < -0.000001 || info_yui_z < -0.000001 || info_ui_z  < -0.000001) printf("\n# ===@@@@@@@===> 2 Probl zi=%d, info_xui_y=%g  info_yui_x=%g  info_ui_y =%g  info_ui_x =%g  info_xuiz_y=%g  info_yuiz_x=%g  info_uiz_y =%g  info_uiz_x =%g  info_xui_z=%g  info_yui_z=%g  info_ui_z=%g \n",zi,info_xui_y,info_yui_x,info_ui_y ,info_ui_x ,info_xuiz_y,info_yuiz_x,info_uiz_y ,info_uiz_x ,info_xui_z,info_yui_z,info_ui_z);
					
					#endif // _MY_DEBUG_NEW


					/////////check maximum mutual infos - cplx terms	
					if(modCplx == MDL) {
						Prui=1;
						//NN = (int)floor(0.5 + ((double)nSample[zi]*sampleSizeEff)/sampleSize);
						//logN=log(NN);
						logN=log(nSample[zi]);

						for(j=2; j < nbrAllVar; j++ ) 
							Prui *= dBin[0][j];
					}

			            //// NI(xy|ui)
					//info3xy_ui = 0.5*(info_xui_y + info_yui_x)*sampleSizeEff/(1.0*sampleSize) ;
					//info2xy_ui = 0.5*(info_ui_y  + info_ui_x)*sampleSizeEff/(1.0*sampleSize)  ;
					info3xy_ui = 0.5*(info_xui_y + info_yui_x);
					info2xy_ui = 0.5*(info_ui_y  + info_ui_x);

					if(modCplx == MDL) {
						logC_xui_y= 0.5*(dBin[0][1] - 1)*(dBin[0][0]*Prui - 1)*logN;
						logC_yui_x= 0.5*(dBin[0][0] - 1)*(dBin[0][1]*Prui - 1)*logN;
						logC_ui_y = 0.5*(dBin[0][1] - 1)*(Prui - 1)*logN;
						logC_ui_x = 0.5*(dBin[0][0] - 1)*(Prui - 1)*logN;
					}
					logC3xy_ui = 0.5*(logC_xui_y + logC_yui_x) ;
					logC2xy_ui = 0.5*(logC_ui_y  + logC_ui_x)  ;

					//testinfo_xy_ui =info3xy_ui - logC3xy_ui + info2xy_ui - logC2xy_ui;
					if(nbrUi==0) testinfo_xy_ui =info3xy_ui - info2xy_ui - logC3xy_ui + logC2xy_ui; //change 20160221
					else         testinfo_xy_ui =info3xy_ui - logC3xy_ui + info2xy_ui - logC2xy_ui;

					//// NI(yz|ui)
					//info3yz_ui = 0.5*(info_uiz_y + info_yui_z)*sampleSizeEff/(1.0*sampleSize) ;
					//info2yz_ui = 0.5*(info_ui_y  + info_ui_z)*sampleSizeEff/(1.0*sampleSize)  ;
					info3yz_ui = 0.5*(info_uiz_y + info_yui_z);
					info2yz_ui = 0.5*(info_ui_y  + info_ui_z);

					if(modCplx == MDL) {
						logC_uiz_y= 0.5*(dBin[0][1] - 1)*(dBin[0][z]*Prui - 1)*logN;
						logC_yui_z= 0.5*(dBin[0][z] - 1)*(dBin[0][1]*Prui - 1)*logN;
						logC_ui_z = 0.5*(dBin[0][z] - 1)*(Prui - 1)*logN;
					}
					logC3yz_ui = 0.5*(logC_uiz_y + logC_yui_z) ;
					logC2yz_ui = 0.5*(logC_ui_y  + logC_ui_z)  ;

					//testinfo_yz_ui =info3yz_ui - logC3yz_ui + info2yz_ui - logC2yz_ui;
					if(nbrUi==0) testinfo_yz_ui =info3yz_ui - info2yz_ui - logC3yz_ui + logC2yz_ui; //change 20160221
					else         testinfo_yz_ui =info3yz_ui - logC3yz_ui + info2yz_ui - logC2yz_ui;

					    //// NI(xz|ui)
					//info3xz_ui = 0.5*(info_xui_z + info_uiz_x)*sampleSizeEff/(1.0*sampleSize) ;
					//info2xz_ui = 0.5*(info_ui_z  + info_ui_x)*sampleSizeEff/(1.0*sampleSize)  ;
					info3xz_ui = 0.5*(info_xui_z + info_uiz_x);
					info2xz_ui = 0.5*(info_ui_z  + info_ui_x);

					if(modCplx == MDL) {
						logC_uiz_x= 0.5*(dBin[0][0] - 1)*(dBin[0][z]*Prui - 1)*logN;
						logC_xui_z= 0.5*(dBin[0][z] - 1)*(dBin[0][0]*Prui - 1)*logN;
					}
					logC3xz_ui = 0.5*(logC_uiz_x + logC_xui_z) ;
					logC2xz_ui = 0.5*(logC_ui_x  + logC_ui_z)  ;

					//testinfo_xz_ui =info3xz_ui - logC3xz_ui + info2xz_ui - logC2xz_ui;
					if(nbrUi==0) testinfo_xz_ui =info3xz_ui - info2xz_ui - logC3xz_ui + logC2xz_ui; //change 20160221
					else         testinfo_xz_ui =info3xz_ui - logC3xz_ui + info2xz_ui - logC2xz_ui;

					    //// NI(xy|uiz)
					//info3xy_uiz = 0.5*(info_xuiz_y + info_yuiz_x)*sampleSizeEff/(1.0*sampleSize) ;
					//info2xy_uiz = 0.5*(info_uiz_y  + info_uiz_x)*sampleSizeEff/(1.0*sampleSize)  ;
					info3xy_uiz = 0.5*(info_xuiz_y + info_yuiz_x);
					info2xy_uiz = 0.5*(info_uiz_y  + info_uiz_x);

					if(modCplx == MDL) {
						Prui *= dBin[0][z];
						logC_xuiz_y= 0.5*(dBin[0][1] - 1)*(dBin[0][0]*Prui - 1)*logN;
						logC_yuiz_x= 0.5*(dBin[0][0] - 1)*(dBin[0][1]*Prui - 1)*logN;
						logC_uiz_y = 0.5*(dBin[0][1] - 1)*(Prui - 1)*logN;
						logC_uiz_x = 0.5*(dBin[0][0] - 1)*(Prui - 1)*logN;
					}
					logC3xy_uiz = 0.5*(logC_xuiz_y + logC_yuiz_x) ;
					logC2xy_uiz = 0.5*(logC_uiz_y  + logC_uiz_x)  ;

					testinfo_xy_uiz =info3xy_uiz -logC3xy_uiz +info2xy_uiz -logC2xy_uiz;

					//// test & store max

					testz=testinfo_xy_ui+testinfo_yz_ui+testinfo_xz_ui+testinfo_xy_uiz; 
					if(max_info_logC < testz){

					  //printf("max countmin %d bins_z %d zi=%d  testz=%g  min_info_logC=%g  max_info_logC=%g \n",countmin,dBin[0][z],zi,testz,min_info_logC,max_info_logC);
					  //for(j=0;j<=z;j++) printf(" dBin[0][j=%d]=%d",j,dBin[0][j]);
					  //printf("\n");

					        //N_xyuiz = (int)floor(0.5 + ((double)nSample[zi]*sampleSizeEff)/sampleSize);
						N_xyuiz = nSample[zi];
						info_xy_ui = info3xy_ui - info2xy_ui; // info NI(xy|ui)
						logC_xy_ui = logC3xy_ui - logC2xy_ui; // cplx k_(xy|ui)
						info_yz_ui = info3yz_ui - info2yz_ui; // info NI(yz|ui)
						logC_yz_ui = logC3yz_ui - logC2yz_ui; // cplx k_(yz|ui)
						info_xz_ui = info3xz_ui - info2xz_ui; // info NI(xz|ui)
						logC_xz_ui = logC3xz_ui - logC2xz_ui; // cplx k_(xz|ui)
						info_xy_uiz =info3xy_uiz -info2xy_uiz; // info NI(xy|uiz)
						logC_xy_uiz =logC3xy_uiz -logC2xy_uiz; // cplx k_(xy|uiz)

						max_info_logC = testz;

						//for(j=0; j < nbrAllVar; j++ ) Opt_dBin[0][j] = dBin[0][j]; change 20160216
						Opt_dBin[0][z] = dBin[0][z];
						min_info_logC = max_info_logC;
					} else if(min_info_logC > testz){
						countmin++;
						min_info_logC = testz;
					} else {
						countmin=0;
						min_info_logC = testz;
					}


					#if _MY_DEBUG_NEW
					
						    printf("\n# =====> test before change zi=%d bin numbers, testz=%g testinfo_xy_ui=%g info3xy_ui=%g  logC3xy_ui=%g  info2xy_ui=%g logC2xy_ui=%g max_info_logC=%g min_info_logC=%g \n",zi,testz,testinfo_xy_ui,info3xy_ui,logC3xy_ui,info2xy_ui,logC2xy_ui,max_info_logC,min_info_logC);
						    printf("\n# =====> test before change zi=%d bin numbers, testz=%g testinfo_xy_uiz=%g info3xy_uiz=%g  logC3xy_uiz=%g  info2xy_uiz=%g logC2xy_uiz=%g max_info_logC=%g min_info_logC=%g \n",zi,testz,testinfo_xy_uiz,info3xy_uiz,logC3xy_uiz,info2xy_uiz,logC2xy_uiz,max_info_logC,min_info_logC);
						    printf("\n# =====> test before change zi=%d bin numbers, testz=%g testinfo_xz_ui=%g info3xz_ui=%g  logC3xz_ui=%g  info2xz_ui=%g logC2xz_ui=%g max_info_logC=%g min_info_logC=%g \n",zi,testz,testinfo_xz_ui,info3xz_ui,logC3xz_ui,info2xz_ui,logC2xz_ui,max_info_logC,min_info_logC);
						    printf("\n# =====> test before change zi=%d bin numbers, testz=%g testinfo_yz_ui=%g info3yz_ui=%g  logC3yz_ui=%g  info2yz_ui=%g logC2yz_ui=%g max_info_logC=%g min_info_logC=%g \n",zi,testz,testinfo_yz_ui,info3yz_ui,logC3yz_ui,info2yz_ui,logC2yz_ui,max_info_logC,min_info_logC);
					
					#endif // _MY_DEBUG_NEW

		
					// int p;
					// for(p = 0; p < nbrAllVar; p++)
					// 	dBin[0][p] = 6;
				
					///////// compute score and store z with max score
					xz  = info_xz_ui - info_xy_ui;
					yz  = info_yz_ui - info_xy_ui;
					xyz = info_xy_ui - info_xy_uiz;

					// cout << "xz: " << xz << endl;
					// cout << "yz: " << yz << endl;
					// cout << "xyz: " << xyz << endl;

					// cout << "nbrLevCorrect[0]: " << nbrLevCorrect[0] << endl;
					// cout << "nbrLevCorrect[1] " << nbrLevCorrect[1] << endl;



					if(k23==TRUE && nbrLevCorrect[0] > 1 && nbrLevCorrect[1] > 1){
						xz  -= logC_xz_ui - logC_xy_ui;
						yz  -= logC_yz_ui - logC_xy_ui;
						xyz -= logC_xy_ui - logC_xy_uiz;
					}
					if(xz < yz){
						first = xz;
						second = yz;
					} else {
						first = yz;
						second = xz;
					}
					dpi = first - log1p(exp(first-second));

					if(xyz < dpi){
						Rzi = xyz;
					} else {
						Rzi = dpi; // final score for each zi ;)
					}

					// cout << endl << ptrzi << " Rzi: " << Rzi << endl << endl;


					if(Rzi > R_top){ //
						N_xyuiz_top = N_xyuiz; 
						NIxy_ui_top  = info_xy_ui;
						k_xy_ui_top  = logC_xy_ui;

						z_top = zi;
						NIxy_uiz_top = info_xy_uiz;
						k_xy_uiz_top = logC_xy_uiz;

						R_top = Rzi;
						NIxyz_ui_top = info_xy_ui - info_xy_uiz;
						//k_xyz_ui_top = logC_xy_ui - logC_xy_uiz;
						k_xyz_ui_top = - logC_xy_ui + logC_xy_uiz; // to fit eq(20) and eq(22) in BMC Bioinfo 2016
					}
				}
			}
	      
		}// end loop on all z


		ptrRetValues[0] = N_xyuiz_top; 
		ptrRetValues[1] = NIxy_ui_top;
		ptrRetValues[2] = k_xy_ui_top;

		ptrRetValues[3] = z_top;
		ptrRetValues[4] = NIxy_uiz_top;
		ptrRetValues[5] = k_xy_uiz_top;

		ptrRetValues[6] = R_top;
		ptrRetValues[7] = NIxyz_ui_top;
		ptrRetValues[8] = k_xyz_ui_top;

		// ptrRetValues old assignment printf("# Info = %g\t, N*Info = %g\t, NML = %g\t, nSample = %g\n\n", (ptrRetValues[0]/ptrRetValues[2]), ptrRetValues[0], ptrRetValues[1], ptrRetValues[2] );

	}
	}
    
	#if _MY_DEBUG_NEW
		printf("\n# =====> test after z \n");
	#endif // _MY_DEBUG_NEW
	// free memory


	for(i=0; i<1;i++)
		free(dBin[i]);
	free(dBin);

	for(i=0; i<1;i++)
		free(Opt_dBin[i]);
	free(Opt_dBin);

	free(nbrLevCorrect);

	delete [] sampleWithZ;


	#if _MY_DEBUG_NEW
		printf("\n# =====> end getAllInfoNEW \n");
		if(nbrZi==0){
		  printf("# N=ptrRetValues[%d]=%g ",0,ptrRetValues[0]);
		  // for(i =1; i <3; i++) printf("# ptrRetValues[%d]=%g ",i,ptrRetValues[i]);
		}
		if(nbrZi > 0){
		  printf("# N=ptrRetValues[%d]=%g ",0,ptrRetValues[0]);
		  for(i =1; i <3; i++) printf("# ptrRetValues[%d]=%g ",i,ptrRetValues[i]);
		  printf("# z=ptrRetValues[%d]=%g ",3,ptrRetValues[3]);
		  for(i =4; i <9; i++) printf(" ptrRetValues[%d]=%g ",i,ptrRetValues[i]);
		}
	#endif // _MY_DEBUG_NEW

	return( ptrRetValues ); 

}



void* evaluateZ(void* container) {
	struct FromToPlusContainer* cont1 = (struct FromToPlusContainer*)container;
	MemorySpace* memoryThread = cont1->memoryThread;
	
	struct ContainerT* cont = cont1->cont;

	int** sortedSample = memoryThread->sortedSample;
	int* ptrAllData = cont->ptrAllData;
	int* dBin = cont1->dBin;

	int* ptrAllLevels = cont->ptrAllLevels;
	int sampleSizeEff = cont->sampleSizeEff;
	int nSample0 = cont->nSample0;
	int ptrzi = cont->ptrzi;
	int sampleSize = cont->sampleSize;

	int* ptrVarIdx = cont->ptrVarIdx;
	bool testDistribution = cont->testDistribution;
	double** freqs1 = cont->freqs1;


	int nbrAllVar = cont->nbrAllVar;
	int nbrUi = cont->nbrUi;
	int k23 = cont->k23;
	int* Opt_dBin = cont1->Opt_dBin;
	int* nbrLevCorrect = cont1->nbrLevCorrect;

	int** Nxuiz = memoryThread->Nxuiz;
	int* Ny = memoryThread->Ny;
	int* Nxui = memoryThread->Nxui;
	int* Nx = memoryThread->Nx;
	int* Nyuiz = memoryThread->Nyuiz;
	int* Nuiz = memoryThread->Nuiz;
	int* Nz = memoryThread->Nz;


	int modCplx = cont->modCplx;
	int* ptrZiIdx = cont->ptrZiIdx;
	double* c2terms = cont->c2terms;
	struct ScoresOfZ** scoresAllZi = cont->scoresAllZi;
	double* weights = cont->weights;

	double* Pxyuiz = new double[200];

	int zpos = cont1->from;
	int to = cont1->to;

	int randomrescaling=1;
	float r,rr;
	int MDL=0,FALSE=0, TRUE=1;

	int z;
	int k = 0, l = 0 ;
	double max_info_logC;
	double min_info_logC;

	double Pxyuizl;

	int Lxyui,Lyui,Lui;
	int Nyui,  Nui;
	int Nxyuis,  Nyuis,  Nuis;

	int NNxyuiz,NNxyuizl,Ntot;  //for rescaling NML change 20160228

	int Prui,X,Y,Z,i;

	double NlogN,logN;

	double info_xui_y,info_yui_x,info_ui_y,info_ui_x;
	double logC_xui_y,logC_yui_x,logC_ui_y,logC_ui_x;

	double info_xuiz_y,info_yuiz_x,info_uiz_y,info_uiz_x;
	double logC_xuiz_y,logC_yuiz_x,logC_uiz_y,logC_uiz_x;

	double info_xui_z,info_yui_z,info_ui_z;
	double logC_xui_z,logC_yui_z,logC_ui_z;

	int  Nxuij, Nxj, Nxuis, Nxs;      //[X]
	int  Nyj,Nys;                     //[Y]
	int Nzs,Nuizs,Nyuizs,Nxyuizs,Nxuizs;

	int  Nyuizl, Nuizl, Nzl; //[Z]
	int Nxuizjl;
	int j;

	double info3xy_ui,info2xy_ui,info3xy_uiz,info2xy_uiz;
	double info3xz_ui,info2xz_ui,info3yz_ui,info2yz_ui;
	double logC3xy_ui,logC2xy_ui,logC3xy_uiz,logC2xy_uiz;
	double logC3xz_ui,logC2xz_ui,logC3yz_ui,logC2yz_ui;


	int N_xyuiz=0;

	double info_xy_ui,info_xy_uiz;
	double info_xz_ui,info_yz_ui;
	double logC_xy_ui,logC_xy_uiz;
	double logC_xz_ui,logC_yz_ui;
	double testinfo_xy_ui,testinfo_xy_uiz;
	double testinfo_xz_ui,testinfo_yz_ui;
	double yz,xz,xyz,first,second,dpi,Rzi,testz;

	int * sampleWithZ = new int[sampleSize];

	int countmin,kz0;

	z=nbrAllVar;
	bool ok;
	while(zpos < to){
		scoresAllZi[zpos] = new ScoresOfZ();

		ptrzi=ptrZiIdx[zpos];

		dBin[z]=ptrAllLevels[ptrzi];



		for( i = 0, k = 0; i < sampleSize; i++)
		{
			// i = sample[k][0];
			ok=TRUE;
			for( j = 0; j < nbrAllVar; j++ )
			{
				// check that X,Y,Uj do not contain NA
				if( ptrAllData[ i + ptrVarIdx[j]*sampleSize ] == -1 )
				{
					ok=FALSE; 
					j=nbrAllVar;
				} 
			}  

			if( ptrAllData[ i + ptrzi*sampleSize ] == -1 )
			{
				ok=FALSE; 		
			}
			

			if(ok==TRUE)
			{
				sampleWithZ[k] = i; // sample number
				k++;
			}
		}

		int nSampleZ = k;

		//update dbin
		bool isGoodCandidate = true;
	    if(nSampleZ < nSample0){
	    	for( j = 0; j < nbrAllVar; j++ ){
	    		std::set<int> s;

			   	for( k = 0; k < nSampleZ; k++ )
				{
					i = sampleWithZ[k];
					s.insert(ptrAllData[ i + ptrVarIdx[j]*sampleSize ]);
				}
				nbrLevCorrect[j] = s.size();
			}


			if(testDistribution){
					int** counts2 = (int **)calloc(ptrAllLevels[ptrVarIdx[0]], sizeof(int*));
			    	for( j = 0; j < ptrAllLevels[ptrVarIdx[0]]; j++ )
			    		counts2[j] = (int *)calloc(ptrAllLevels[ptrVarIdx[1]], sizeof(int));

			    	// fill table
					for(int k = 0; k < nSampleZ; k++){
						i = sampleWithZ[k];
						counts2[ ptrAllData[ i + ptrVarIdx[0]*sampleSize ]][ ptrAllData[ i + ptrVarIdx[1]*sampleSize ]]++;
					}

					// double cplxMdl = (nbrUi + 3 ) * 5;
					double cplxMdl = log(nSampleZ);
					// double cplxMdl = ptrAllLevels[ptrVarIdx[0]] * ptrAllLevels[ptrVarIdx[1]] * log(nSampleZ) ; //*(ptrAllLevels[ptrzi] - 1) ;
					// for(int counterU = 0 ; counterU < nbrUi; counterU++)
					// 	cplxMdl *= ptrAllLevels[ptrVarIdx[counterU+2]];

					double kldiv = exp(-nSampleZ*kl(counts2, freqs1, ptrAllLevels[ptrVarIdx[0]], ptrAllLevels[ptrVarIdx[1]]) + cplxMdl);
					
					
					// cout << kldiv << endl;

					if(kldiv < 1)
						isGoodCandidate = false;

					for( j = 0; j < ptrAllLevels[ptrVarIdx[0]]; j++ )
						delete [] counts2[j];
					delete [] counts2;
				}








		}







		if(isGoodCandidate){
			for( k = 0; k < nSample0; k++ ){ 
				i=sortedSample[k][0];
				// find the first sample for which zi does not contain NA
				if( ptrAllData[ i + ptrzi*sampleSize ] > -1 )
				{
					kz0=k;
					nSampleZ=1;
					k=nSample0;
				}
			}

			if(nSampleZ == 1){

				max_info_logC = -LARGE;
				min_info_logC = LARGE;

				countmin=0;

			    nSampleZ=0; 
				Lxyui = sortedSample[kz0][1]; // min Lxyui
				Lyui  = sortedSample[kz0][2]; // min Lyui
				Lui   = sortedSample[kz0][3]; // min Lui

				NNxyuiz = 0;
				NNxyuizl= 0;

				Nyui  = 0;
				Nui   = 0;
				for( k = 0; k < dBin[0]; k++ ){ 
					Nxui[k]=0;
					Nx[k]=0;
				
					for( l = 0; l < dBin[z]; l++ ){
						Nxuiz[k][l]=0;
					}
				}


				for( k = 0; k < dBin[1]; k++ ){
					Ny[k]=0;
				}
				for( l = 0; l < dBin[z]; l++ ){
					Pxyuiz[l]=0;
					Nyuiz[l]=0;
					Nuiz[l]=0;
					Nz[l]=0;
				}

				X=sortedSample[kz0][4];
				Y=sortedSample[kz0][5];
				i=sortedSample[kz0][0];
				Z=ptrAllData[ i + ptrzi*sampleSize ];  // first Z

				Pxyuiz[Z]=weights[sortedSample[kz0][0]]; //segfault

				sortedSample[nSample0][0]=i; // to terminate loop properly below


				info_xui_y = 0.0;
				info_yui_x = 0.0;
				info_ui_y  = 0.0;
				info_ui_x  = 0.0;
				logC_xui_y = 0.0;
				logC_yui_x = 0.0;
				logC_ui_y  = 0.0;
				logC_ui_x  = 0.0;

				info_xuiz_y = 0.0;
				info_yuiz_x = 0.0;
				info_uiz_y  = 0.0;
				info_uiz_x  = 0.0;
				logC_xuiz_y = 0.0;
				logC_yuiz_x = 0.0;
				logC_uiz_y  = 0.0;
				logC_uiz_x  = 0.0;

				info_xui_z = 0.0;
				info_yui_z = 0.0;
				info_ui_z  = 0.0;
				logC_xui_z = 0.0;
				logC_yui_z = 0.0;
				logC_ui_z  = 0.0;

				info_xy_ui  = 0.0;
				info_xy_uiz = 0.0;
				info_xz_ui  = 0.0;
				info_yz_ui  = 0.0;
				logC_xy_ui  = 0.0;
				logC_xy_uiz = 0.0;
				logC_xz_ui  = 0.0;
				logC_yz_ui  = 0.0;

				Nxyuis=0; // 11 test variables for counts, should all sum to nSample0
				Nyuis=0;
				Nxuis=0;
				Nuis=0;
				Nxs=0;
				Nys=0;

				Nzs=0;
				Nuizs=0;
				Nyuizs=0;
				Nxyuizs=0;
				Nxuizs=0;

				/////////make the counts and compute mutual infos & logCs 
				for( k = kz0+1; k <= nSample0; k++ ){ // change 20160220

					i=sortedSample[k][0];

					// check whether zi does not contain NA
					if( ptrAllData[ i + ptrzi*sampleSize ] > -1 )
					{

						Z=ptrAllData[ i + ptrzi*sampleSize ];  // Z

						if(sortedSample[k][1] > Lxyui){

						        NNxyuiz=0;
							for( l = 0; l < dBin[z]; l++ ){
								//Nxyuizl=Nxyuiz[l];//
								Pxyuizl=Pxyuiz[l];//
								if( Pxyuizl > 0 ){
								    if(sampleSizeEff!=sampleSize){
								        if(randomrescaling){
									    NNxyuizl=(int)floor(Pxyuizl);
									    r=Pxyuizl - NNxyuiz;
									    rr=(double)rand() / (double)RAND_MAX ;
									    // printf("%f\n",rr); 
									    if (r > rr) NNxyuizl++;
									    //printf("# =====> rando  r=%g rr=%g NNxyui=%d \n",r,rr,NNxyui);
									}
									// else{
									//     NNxyuizl=(int)floor(0.5+((double)Nxyuizl*sampleSizeEff)/sampleSize);
									// }
								    }else{
								        NNxyuizl=(int)Pxyuizl;
								    }

								    if( NNxyuizl > 0 ){
								        NlogN = NNxyuizl * log(NNxyuizl);
									info_xuiz_y += NlogN;//
									info_yuiz_x += NlogN;//
									NNxyuiz += NNxyuizl;

									Nz[l]   += NNxyuizl;	//4
									Nuiz[l] += NNxyuizl;
									Nyuiz[l]+= NNxyuizl;
									Nxuiz[X][l] += NNxyuizl;
								    }

								    //Nxyuiz[l]=0;//
								    Pxyuiz[l]=0;//
								}
							}

							Pxyuiz[Z]=weights[sortedSample[k][0]];

							if(NNxyuiz > 0){


							    NlogN = NNxyuiz * log(NNxyuiz);
							    info_xui_y += NlogN;//
							    info_yui_x += NlogN;//

							    nSampleZ += NNxyuiz;
						
							    Nxyuizs += NNxyuiz;

							    Nui     += NNxyuiz; 
							    Nyui    += NNxyuiz; 
							    Nxui[X] += NNxyuiz; 
								  
							    Nx[X]   += NNxyuiz;	//4
							    Ny[Y]   += NNxyuiz;	//4
							    Ntot    += NNxyuiz;
								  
							    Nxyuis  += NNxyuiz;
							}


							Lxyui = sortedSample[k][1];

							if(k<nSample0) X=sortedSample[k][4];
							if(k<nSample0) Y=sortedSample[k][5];


								if(sortedSample[k][2] > Lyui){
							    if( Nyui > 0 ){
							        NlogN = Nyui * log(Nyui);
								info_yui_x -= NlogN;//
								info_yui_z -= NlogN;//
								info_ui_y  += NlogN;//
								if(modCplx != MDL){
								    logC_yui_x += compute_LogC_C2(Nyui,dBin[0],c2terms);
								    logC_yui_z += compute_LogC_C2(Nyui,dBin[z],c2terms);
								}

								for( l = 0; l < dBin[z]; l++ ){
								    Nyuizl=Nyuiz[l];
								    if( Nyuizl > 0 ){
								        NlogN = Nyuizl * log(Nyuizl);
									info_yui_z += NlogN;//
									info_uiz_y += NlogN;//
									info_yuiz_x -= NlogN;//
									if(modCplx != MDL){
									    logC_yuiz_x += compute_LogC_C2(Nyuizl,dBin[0],c2terms);
									}
									Nyuizs += Nyuizl;
									Nyuiz[l]=0;//
								    }
								}
								//Nyuiz[Z]=1;  //2
							      
								Nyuis += Nyui;
								Nyui = 0; //2
							    }
							    Lyui = sortedSample[k][2];

							    if(sortedSample[k][3] > Lui){
							         if( Nui > 0 ){
								        NlogN = Nui * log(Nui);
									info_ui_x  -= NlogN;//
									info_ui_y  -= NlogN;//
									info_ui_z  -= NlogN;//
									if(modCplx != MDL){
										logC_ui_x += compute_LogC_C2(Nui,dBin[0], c2terms);//
										logC_ui_y += compute_LogC_C2(Nui,dBin[1], c2terms);//
										logC_ui_z += compute_LogC_C2(Nui,dBin[z], c2terms);//
									}
									Nuis += Nui;
									Nui = 0;   //3

									for( l = 0; l < dBin[z]; l++ ){
										Nuizl=Nuiz[l];
										if( Nuizl > 0 ){
											NlogN = Nuizl * log(Nuizl);
											info_ui_z  += NlogN;//
											info_uiz_x -= NlogN;//
											info_uiz_y -= NlogN;//
											if(modCplx != MDL){
											  logC_uiz_x += compute_LogC_C2(Nuizl,dBin[0], c2terms);//
											  logC_uiz_y += compute_LogC_C2(Nuizl,dBin[1], c2terms);//
											}
											Nuizs += Nuizl;
											Nuiz[l]=0;//
										}
									}
									//Nuiz[Z]=1;   //3

									
									for( j = 0; j < dBin[0]; j++ ){ 
										Nxuij=Nxui[j];
										if( Nxuij > 0 ){
											NlogN = Nxuij * log(Nxuij);
											info_xui_y -= NlogN;//
											info_xui_z -= NlogN;//
											info_ui_x  += NlogN;//
											if(modCplx != MDL){
											  logC_xui_y += compute_LogC_C2(Nxuij,dBin[1], c2terms);//
											  logC_xui_z += compute_LogC_C2(Nxuij,dBin[z], c2terms);//
											}
											Nxuis += Nxuij;
											Nxui[j]=0;//

											for( l = 0; l < dBin[z]; l++ ){
												Nxuizjl=Nxuiz[j][l];
												if( Nxuizjl > 0 ){
													NlogN = Nxuizjl * log(Nxuizjl);
													info_xui_z += NlogN;//
													info_uiz_x += NlogN;//
													info_xuiz_y -= NlogN;//
													if(modCplx != MDL){
													  logC_xuiz_y += compute_LogC_C2(Nxuizjl,dBin[1], c2terms);//
													}
													Nxuizs += Nxuizjl;
													Nxuiz[j][l]=0;//
											  	}
											}
											//Nxuiz[X][Z]=1;// change

										}
		       							}

							         }
							         Lui = sortedSample[k][3];

							    }

							}

						} else {
							Pxyuiz[Z]+=weights[sortedSample[k][0]];
						}
					}
				}


		    	// increment info with Nx[X], Ny[Y] and Nz[Z] contributions	
				for( j = 0; j < dBin[0]; j++ ){ 
					Nxj=Nx[j];
					if( Nxj > 0 ){
						NlogN = Nxj * log(Nxj/(1.0*nSampleZ));
						info_yui_x  -= NlogN;//
						info_ui_x   -= NlogN;//
						info_uiz_x  -= NlogN;//
						info_yuiz_x -= NlogN;//
						Nxs += Nxj;
						Nx[j]=0;
					}
				}
				for( j = 0; j < dBin[1]; j++ ){ 
					Nyj=Ny[j];
					if( Nyj > 0 ){
						NlogN = Nyj * log(Nyj/(1.0*nSampleZ));
						info_xui_y  -= NlogN;//
						info_ui_y   -= NlogN;//
						info_uiz_y  -= NlogN;//
						info_xuiz_y -= NlogN;//
						Nys += Nyj;
						Ny[j]=0;
					}
				}
				for( l = 0; l < dBin[z]; l++ ){ 
					Nzl=Nz[l];
						if( Nzl > 0 ){
						NlogN = Nzl * log(Nzl/(1.0*nSampleZ));
						info_xui_z  -= NlogN;//
						info_yui_z  -= NlogN;//
						info_ui_z   -= NlogN;//
						Nzs += Nzl;
						Nz[l]=0;
					}
				}

				/////////check maximum mutual infos - cplx terms	
				if(modCplx == MDL) {
					Prui=1;
					//NN = (int)floor(0.5 + ((double)nSampleZ*sampleSizeEff)/sampleSize);
					//logN=log(NN);
					logN=log(nSampleZ);

					for(j=2; j < nbrAllVar; j++ ) 
						Prui *= dBin[j];
				}

		            //// NI(xy|ui)
				//info3xy_ui = 0.5*(info_xui_y + info_yui_x)*sampleSizeEff/(1.0*sampleSize) ;
				//info2xy_ui = 0.5*(info_ui_y  + info_ui_x)*sampleSizeEff/(1.0*sampleSize)  ;
				info3xy_ui = 0.5*(info_xui_y + info_yui_x);
				info2xy_ui = 0.5*(info_ui_y  + info_ui_x);

				if(modCplx == MDL) {
					logC_xui_y= 0.5*(dBin[1] - 1)*(dBin[0]*Prui - 1)*logN;
					logC_yui_x= 0.5*(dBin[0] - 1)*(dBin[1]*Prui - 1)*logN;
					logC_ui_y = 0.5*(dBin[1] - 1)*(Prui - 1)*logN;
					logC_ui_x = 0.5*(dBin[0] - 1)*(Prui - 1)*logN;
				}
				logC3xy_ui = 0.5*(logC_xui_y + logC_yui_x) ;
				logC2xy_ui = 0.5*(logC_ui_y  + logC_ui_x)  ;

				//testinfo_xy_ui =info3xy_ui - logC3xy_ui + info2xy_ui - logC2xy_ui;
				if(nbrUi==0) testinfo_xy_ui =info3xy_ui - info2xy_ui - logC3xy_ui + logC2xy_ui; //change 20160221
				else         testinfo_xy_ui =info3xy_ui - logC3xy_ui + info2xy_ui - logC2xy_ui;

				//// NI(yz|ui)
				//info3yz_ui = 0.5*(info_uiz_y + info_yui_z)*sampleSizeEff/(1.0*sampleSize) ;
				//info2yz_ui = 0.5*(info_ui_y  + info_ui_z)*sampleSizeEff/(1.0*sampleSize)  ;
				info3yz_ui = 0.5*(info_uiz_y + info_yui_z);
				info2yz_ui = 0.5*(info_ui_y  + info_ui_z);

				if(modCplx == MDL) {
					logC_uiz_y= 0.5*(dBin[1] - 1)*(dBin[z]*Prui - 1)*logN;
					logC_yui_z= 0.5*(dBin[z] - 1)*(dBin[1]*Prui - 1)*logN;
					logC_ui_z = 0.5*(dBin[z] - 1)*(Prui - 1)*logN;
				}
				logC3yz_ui = 0.5*(logC_uiz_y + logC_yui_z) ;
				logC2yz_ui = 0.5*(logC_ui_y  + logC_ui_z)  ;

				//testinfo_yz_ui =info3yz_ui - logC3yz_ui + info2yz_ui - logC2yz_ui;
				if(nbrUi==0) testinfo_yz_ui =info3yz_ui - info2yz_ui - logC3yz_ui + logC2yz_ui; //change 20160221
				else         testinfo_yz_ui =info3yz_ui - logC3yz_ui + info2yz_ui - logC2yz_ui;

				    //// NI(xz|ui)
				//info3xz_ui = 0.5*(info_xui_z + info_uiz_x)*sampleSizeEff/(1.0*sampleSize) ;
				//info2xz_ui = 0.5*(info_ui_z  + info_ui_x)*sampleSizeEff/(1.0*sampleSize)  ;
				info3xz_ui = 0.5*(info_xui_z + info_uiz_x);
				info2xz_ui = 0.5*(info_ui_z  + info_ui_x);

				if(modCplx == MDL) {
					logC_uiz_x= 0.5*(dBin[0] - 1)*(dBin[z]*Prui - 1)*logN;
					logC_xui_z= 0.5*(dBin[z] - 1)*(dBin[0]*Prui - 1)*logN;
				}
				logC3xz_ui = 0.5*(logC_uiz_x + logC_xui_z) ;
				logC2xz_ui = 0.5*(logC_ui_x  + logC_ui_z)  ;

				//testinfo_xz_ui =info3xz_ui - logC3xz_ui + info2xz_ui - logC2xz_ui;
				if(nbrUi==0) testinfo_xz_ui =info3xz_ui - info2xz_ui - logC3xz_ui + logC2xz_ui; //change 20160221
				else testinfo_xz_ui =info3xz_ui - logC3xz_ui + info2xz_ui - logC2xz_ui;

				    //// NI(xy|uiz)
				//info3xy_uiz = 0.5*(info_xuiz_y + info_yuiz_x)*sampleSizeEff/(1.0*sampleSize) ;
				//info2xy_uiz = 0.5*(info_uiz_y  + info_uiz_x)*sampleSizeEff/(1.0*sampleSize)  ;
				info3xy_uiz = 0.5*(info_xuiz_y + info_yuiz_x);
				info2xy_uiz = 0.5*(info_uiz_y  + info_uiz_x);

				if(modCplx == MDL) {
					Prui *= dBin[z];
					logC_xuiz_y= 0.5*(dBin[1] - 1)*(dBin[0]*Prui - 1)*logN;
					logC_yuiz_x= 0.5*(dBin[0] - 1)*(dBin[1]*Prui - 1)*logN;
					logC_uiz_y = 0.5*(dBin[1] - 1)*(Prui - 1)*logN;
					logC_uiz_x = 0.5*(dBin[0] - 1)*(Prui - 1)*logN;
				}
				logC3xy_uiz = 0.5*(logC_xuiz_y + logC_yuiz_x) ;
				logC2xy_uiz = 0.5*(logC_uiz_y  + logC_uiz_x)  ;

				testinfo_xy_uiz =info3xy_uiz -logC3xy_uiz +info2xy_uiz -logC2xy_uiz;

				//// test & store max

				testz=testinfo_xy_ui+testinfo_yz_ui+testinfo_xz_ui+testinfo_xy_uiz; 
				if(max_info_logC < testz){

					N_xyuiz = nSampleZ;
					info_xy_ui = info3xy_ui - info2xy_ui; // info NI(xy|ui)
					logC_xy_ui = logC3xy_ui - logC2xy_ui; // cplx k_(xy|ui)
					info_yz_ui = info3yz_ui - info2yz_ui; // info NI(yz|ui)
					logC_yz_ui = logC3yz_ui - logC2yz_ui; // cplx k_(yz|ui)
					info_xz_ui = info3xz_ui - info2xz_ui; // info NI(xz|ui)
					logC_xz_ui = logC3xz_ui - logC2xz_ui; // cplx k_(xz|ui)
					info_xy_uiz =info3xy_uiz -info2xy_uiz; // info NI(xy|uiz)
					logC_xy_uiz =logC3xy_uiz -logC2xy_uiz; // cplx k_(xy|uiz)

					max_info_logC = testz;

					//for(j=0; j < nbrAllVar; j++ ) Opt_dBin[j] = dBin[j]; change 20160216
					Opt_dBin[z] = dBin[z];
					min_info_logC = max_info_logC;
				} else if(min_info_logC > testz){
					countmin++;
					min_info_logC = testz;
				} else {
					countmin=0;
					min_info_logC = testz;
				}
			
				///////// compute score and store z with max score
				xz  = info_xz_ui - info_xy_ui;
				yz  = info_yz_ui - info_xy_ui;
				xyz = info_xy_ui - info_xy_uiz;

				if(k23==TRUE && nbrLevCorrect[0] > 1 && nbrLevCorrect[1] > 1){
					xz  -= logC_xz_ui - logC_xy_ui;
					yz  -= logC_yz_ui - logC_xy_ui;
					xyz -= logC_xy_ui - logC_xy_uiz;
				}

				if(xz < yz){
					first = xz;
					second = yz;
				} else {
					first = yz;
					second = xz;
				}
				dpi = first - log1p(exp(first-second));

				if(xyz < dpi){
					Rzi = xyz;
				} else {
					Rzi = dpi; // final score for each zi ;)
				}
				
				scoresAllZi[zpos]->N_xyuiz = N_xyuiz; 
				scoresAllZi[zpos]->info_xy_ui  = info_xy_ui;
				scoresAllZi[zpos]->logC_xy_ui  = logC_xy_ui;

				scoresAllZi[zpos]->info_xy_uiz = info_xy_uiz;
				scoresAllZi[zpos]->logC_xy_uiz = logC_xy_uiz;

				scoresAllZi[zpos]->Rzi = Rzi;	
			
			} else {
				// not a good candidate
				scoresAllZi[zpos]->N_xyuiz = 0; 
				scoresAllZi[zpos]->info_xy_ui  = 0;
				scoresAllZi[zpos]->logC_xy_ui  = 0;

				scoresAllZi[zpos]->info_xy_uiz = 0;
				scoresAllZi[zpos]->logC_xy_uiz = 0;

				scoresAllZi[zpos]->Rzi = 0;	
			}

		}
		zpos += 1;
	}
	delete [] sampleWithZ;
	delete [] Pxyuiz;
	return NULL;
}



// modCplx = myCplx = 0 --> MDL, myCplx = 1 --> NML
// If nbrZi== 0, return nSample[0]     & nSample[0]*I(xy|{ui})      & k_xy_ui 
// If nbrZi > 0, return nSample[z_top] & nSample[z_top]*I(xy|{ui})  & k_xy_ui 
//                      z_top          & nSample[z_top]*I(xy|{ui}z) & k_xy_uiz 
//                      R_top          & nSample[z_top]*I(xyz|{ui}) & k_xyz_ui 
double* getAllInfoNEWThreads( int* ptrAllData, int* ptrAllLevels, int* ptrVarIdx, int nbrUi, int* ptrZiIdx, int nbrZi, int ziPos, int sampleSize, 
	int sampleSizeEff, int modCplx, int k23, double* c2terms, MemorySpace* memory, int nthreadsMax, MemorySpace* memoryThreads, double* weights, 
	double** freqs1, bool testDistribution)
{	
	int randomrescaling=1;
	float r,rr;

	int bin_max=100,MDL=0, TRUE=1, FALSE=0;
	int ok;
	int **sample,**sortedSample,**Opt_sortedSample; //[N+1][7]

	int ncol=7;

	sample = (*memory).sample;

	sortedSample = (*memory).sortedSample;

	Opt_sortedSample = (*memory).Opt_sortedSample;

	int nSample0,*orderSample,*sampleKey;//[N+1 elements: 0 to N]

	orderSample = (*memory).orderSample;
	sampleKey = (*memory).sampleKey;

	int bin,PBin,Prui,increment,X,Y;
	int ptrzi,zi;

	double NlogN,logN;

	int Lxyui,Lyui,Lui;
	double Pxyui;
	int  Nyui,  Nui;
	int  Nxyuis,  Nyuis,  Nuis;

	int NNxyui,Ntot;  //for rescaling NML change 20160228

	int *Nxyuiz, *Nyuiz, *Nuiz, *Nz, *bridge;  //[Z]

	Nxyuiz = (*memory).Nxyuiz;
	Nyuiz = (*memory).Nyuiz;
	Nuiz = (*memory).Nuiz;
	Nz = (*memory).Nz;
	bridge = (*memory).bridge;

	int *Ny;     
	Ny = (*memory).Ny;
	int  Nyj,Nys;                     //[Y]

	int *Nxui, *Nx;                   //[X]

	Nxui = (*memory).Nxui;
	Nx = (*memory).Nx;

	int  Nxuij, Nxj, Nxuis, Nxs;      //[X]

	int **Nxuiz;                      //[X][Z]

	ncol = bin_max+1;

	Nxuiz = (*memory).Nxuiz;

	double info_xui_y,info_yui_x,info_ui_y,info_ui_x;
	double logC_xui_y,logC_yui_x,logC_ui_y,logC_ui_x;

	double info3xy_ui,info2xy_ui;
	double logC3xy_ui,logC2xy_ui;


	double testinfo_xy_ui;

	//int N_xy_ui=-1,N_xyuiz=-1;
	int N_xy_ui=0;
	//double LARGE DBL_MAX;
	double min_info_logC=LARGE, max_info_logC= - LARGE;
	//double min_info_logC = 100000.0, max_info_logC = -100000.0;

	int countmin;

	int z_top,N_xyuiz_top;
	double R_top;
	double NIxy_ui=-1.0,NIxy_ui_top,NIxy_uiz_top,NIxyz_ui_top;
	double k_xy_ui=-1.0,k_xy_ui_top,k_xy_uiz_top,k_xyz_ui_top;
  
  /////////////////////////

	int i, j, k;	// for loops

	int nbrAllVar;

	int nbrRetValues = 3;

	// If no zi, return nSample0 & nSample0*I(xy|{ui}) & k_xy_ui 
	// If nbrZi > 0, return nSample[z_top] & nSample[z_top]*I(xy|{ui}) & k_xy_ui 
	//                      z_top & nSample[z_top]*I(xy|{ui}z) & k_xy_uiz 
	//                      R_top & nSample[z_top]*I(xyz|{ui}) & k_xyz_ui 
	if( nbrZi > 0 )	{ nbrRetValues = 9; }	

	double *ptrRetValues = (double*)malloc(nbrRetValues* sizeof(double));	

	ptrRetValues[0] = -1; 
	ptrRetValues[1] = -1;
	ptrRetValues[2] = -1;

	if( nbrZi > 0 ){
		ptrRetValues[3] = -1; 
		ptrRetValues[4] = -1;
		ptrRetValues[5] = -1;
		ptrRetValues[6] = -1; 
		ptrRetValues[7] = -1;
		ptrRetValues[8] = -1;
	}	

	// Define the total number of variables (xi, xj, {ui}, {zi})
	nbrAllVar = (nbrUi+2);

	int *dBin,*Opt_dBin; //[1][Nb_ui+2]

	ncol = nbrAllVar +1;
	
    dBin = (int *)malloc(ncol* sizeof(int));
	Opt_dBin = (int *)malloc(ncol* sizeof(int));
	int* nbrLevCorrect = (int *)malloc(ncol * sizeof(int));

	// find samples without NA in x,y,ui and store their id in sample[k][0]
	for( i = 0, k = 0; i < sampleSize; i++)
	{
		ok=TRUE;
		for( j = 0; j < nbrAllVar; j++ )
		{
			// check that X,Y,Uj do not contain NA
			if( ptrAllData[ i + ptrVarIdx[j]*sampleSize ] == -1 )
			{
				ok=FALSE; 
				j=nbrAllVar;
			} 
		}  

		if(nbrZi==1) // if only one variable zi (to estimate NI(xy|ui) and 
		{          // NI(xy|uiz) on the exact same samples in case of NA)
			ptrzi=ptrZiIdx[0];
			if( ptrAllData[ i + ptrzi*sampleSize ] == -1 )
			{
				ok=FALSE; 		
			}
		}

		if(ok==TRUE)
		{
			sample[k][0] = i; // sample number
			k++;
		}
	}
    nSample0=k;

    //update dbin
    if(nSample0 < sampleSize){
    	for( j = 0; j < nbrAllVar; j++ ){
    		std::set<int> s;

		   	for( k = 0; k < nSample0; k++ )
			{
				i = sample[k][0];
				s.insert(ptrAllData[ i + ptrVarIdx[j]*sampleSize ]);
			}
			nbrLevCorrect[j] = s.size();
		}
	} else {
		for( j = 0; j < nbrAllVar; j++ )
			nbrLevCorrect[j]=ptrAllLevels[ptrVarIdx[j]];
	}

	if(nSample0 > 0){

		for(j=0; j < nbrAllVar; j++ )
		{

			dBin[j]=ptrAllLevels[ptrVarIdx[j]];

			Opt_dBin[j] = dBin[j];
		}

		//compute Lxyui, Lyui, Lui indices for later counting purpose
		for( k = 0; k < nSample0; k++ )
		{ 
			i=sample[k][0]; // sample number

			bin=ptrAllData[ i + ptrVarIdx[0]*sampleSize ];


			sample[k][1] = bin; // Lxyui initialisation
			sample[k][4] = bin; // X

			bin=ptrAllData[ i + ptrVarIdx[1]*sampleSize ];


			sample[k][5] = bin; // Y
			PBin=dBin[0];		
			increment = bin*PBin;
			sample[k][1] += increment; // Lxyui
			sample[k][2] = increment; // Lyui initialisation
			sample[k][3] = 0; // Lui initialisation
	  
			for( j = 2; j < nbrAllVar; j++ )
			{
				bin=ptrAllData[ i + ptrVarIdx[j]*sampleSize ];

				PBin *= dBin[j-1];
				increment = bin*PBin;
				sample[k][1] += increment; // Lxyui
				sample[k][2] += increment; // Lyui
				sample[k][3] += increment; // Lui
			} 	  
		}

		//k++;
		bin = PBin*dBin[nbrAllVar-1];

		sample[k][0] = nSample0; // extra sample id (useful for termination of counts)
		sample[k][1] = bin; // max Lxyui (useful for termination of counts)
		sample[k][2] = bin; // max Lyui  (useful for termination of counts)
		sample[k][3] = bin; // max Lui   (useful for termination of counts)


		/////////sort sample in increasing Lxyui stored in sample[k][1]
		//for( k = 0; k <= nSample0; k++ ){ 
		for( k = 1; k <=nSample0+1; k++ ){ 
			orderSample[k] = k-1;
			sampleKey[k]=sample[k-1][1]; // will sort in increasing Lxyui
		}

		sort2arrays(nSample0+1,sampleKey,orderSample, bridge);

		for( k = 1; k <= nSample0+1; k++ ){ 

			i=orderSample[k];

			for( j=0; j< 6; j++ ){

			    sortedSample[k-1][j]=sample[i][j];
			}
		}

		/////////initialization of counts and mutual infos & logCs 
		Lxyui = sortedSample[0][1]; // min Lxyui
		Lyui  = sortedSample[0][2]; // min Lyui
		Lui   = sortedSample[0][3]; // min Lui

		Pxyui = weights[sortedSample[0][0]];
		NNxyui= 0;
		Nyui  = 0;
		Nui   = 0;

		for( k = 0; k < dBin[0]; k++ ){ 
		  Nxui[k]=0;
		  Nx[k]=0;
		}

		for( k = 0; k < dBin[1]; k++ ) Ny[k]=0;
		X=sortedSample[0][4];
		Y=sortedSample[0][5];

		info_xui_y = 0.0;
		info_yui_x = 0.0;
		info_ui_y  = 0.0;
		info_ui_x  = 0.0;
		logC_xui_y = 0.0;
		logC_yui_x = 0.0;
		logC_ui_y  = 0.0;
		logC_ui_x  = 0.0;

		Nxyuis=0; // 6 test variables for counts, should all sum to nSample0
		Nyuis=0;
		Nxuis=0;
		Nuis=0;
		Nxs=0;
		Nys=0;
		Ntot=0;

		/////////make the counts and compute mutual infos & logCs 
		for( k = 1; k <= nSample0; k++ ){ 

			if(sortedSample[k][1] > Lxyui){

			    NNxyui=0;
			    if(sampleSizeEff!=sampleSize){
				  if(randomrescaling){
				    NNxyui=(int)floor(Pxyui);
				    r=Pxyui - NNxyui;
				    rr=(double)rand() / (double)RAND_MAX ;
				    // printf("%f\n",rr); 
				    if (r > rr) NNxyui++;
				    //printf("# =====> rando  r=%g rr=%g NNxyui=%d \n",r,rr,NNxyui);
				  }
				}else{
				  //NNxyui=Nxyui; 
				  NNxyui=(int)Pxyui; 
				}
				
				Nui     += NNxyui; 
				Nyui    += NNxyui; 
				Nxui[X] += NNxyui; 
				
				Nx[X] += NNxyui;	//4
				Ny[Y] += NNxyui;	//4
				Ntot +=  NNxyui;

				if( NNxyui > 0 ){
				  NlogN = NNxyui * log(NNxyui);
				  info_xui_y += NlogN;
				  info_yui_x += NlogN;
				}
				Lxyui = sortedSample[k][1];
				Nxyuis += NNxyui;
				Pxyui = weights[sortedSample[k][0]]; //weights[k];	

				if(k<nSample0) X=sortedSample[k][4];
				if(k<nSample0) Y=sortedSample[k][5];
		

				if(sortedSample[k][2] > Lyui){

				        if( Nyui > 0 ){
					  NlogN = Nyui * log(Nyui);
					  info_yui_x -= NlogN;
					  info_ui_y  += NlogN;
					  
					  if(modCplx != MDL){
					    logC_yui_x += compute_LogC_C2(Nyui,dBin[0],c2terms);
					  }
					}
					Lyui = sortedSample[k][2];
					Nyuis += Nyui;
					Nyui = 0;

					if(sortedSample[k][3] > Lui){ 
						for( j = 0; j < dBin[0]; j++ ){ 
							Nxuij=Nxui[j];
							if( Nxuij > 0 ){
								NlogN = Nxuij * log(Nxuij);
								info_xui_y -= NlogN;
								info_ui_x  += NlogN;
								if(modCplx != MDL){
								    logC_xui_y += compute_LogC_C2(Nxuij,dBin[1],c2terms);
								}
				   
								Nxuis += Nxuij;
								Nxui[j]=0;
							}
						}
						
						if( Nui > 0 ){
						  NlogN = Nui * log(Nui);
						  info_ui_y  -= NlogN;
						  info_ui_x  -= NlogN;
						  if(modCplx != MDL){
						    logC_ui_x += compute_LogC_C2(Nui,dBin[0],c2terms);
						    logC_ui_y += compute_LogC_C2(Nui,dBin[1],c2terms);
						  }
						}
						Lui = sortedSample[k][3];
						Nuis += Nui;
						Nui = 0;

					}
				}

			} else {
			  Pxyui += weights[sortedSample[k][0]]; //weights[k];	
			}
		}

		// increment for info for Nx[X] and Ny[Y] contributions	
		for( j = 0; j < dBin[0]; j++ ){ 
			Nxj=Nx[j];
			if( Nxj > 0 ){
			      //NlogN = Nxj * log(Nxj/(1.0*nSample0));
				NlogN = Nxj * log(Nxj/(1.0*Ntot));
				info_yui_x -= NlogN;
				info_ui_x  -= NlogN;
				Nxs += Nxj;
			}
		}
		for( j = 0; j < dBin[1]; j++ ){ 
			Nyj=Ny[j];
			if( Nyj > 0 ){
			      //NlogN = Nyj * log(Nyj/(1.0*nSample0));
				NlogN = Nyj * log(Nyj/(1.0*Ntot));
				info_xui_y -= NlogN;
				info_ui_y  -= NlogN;
				Nys += Nyj;
			}
		}
    
		/////////check maximum mutual infos - cplx terms
		info3xy_ui = 0.5*(info_xui_y + info_yui_x);
		info2xy_ui = 0.5*(info_ui_y  + info_ui_x);

		if(modCplx == MDL) {
			Prui=1;

			logN=log(Ntot);
			for(j=2; j < nbrAllVar; j++ ) Prui *= dBin[j];
			logC_xui_y= 0.5*(dBin[1] - 1)*(dBin[0]*Prui - 1)*logN;
			logC_yui_x= 0.5*(dBin[0] - 1)*(dBin[1]*Prui - 1)*logN;
			logC_ui_y = 0.5*(dBin[1] - 1)*(Prui - 1)*logN;
			logC_ui_x = 0.5*(dBin[0] - 1)*(Prui - 1)*logN;
		}

		logC3xy_ui = 0.5*(logC_xui_y + logC_yui_x) ;
		logC2xy_ui = 0.5*(logC_ui_y  + logC_ui_x)  ;
	 
		if(nbrUi==0) testinfo_xy_ui =info3xy_ui - info2xy_ui - logC3xy_ui + logC2xy_ui; //change 20160221
		else         testinfo_xy_ui =info3xy_ui - logC3xy_ui + info2xy_ui - logC2xy_ui;

		if(max_info_logC < testinfo_xy_ui){

			N_xy_ui = Ntot;
			NIxy_ui = info3xy_ui - info2xy_ui; // info to be returned if no z
			k_xy_ui = logC3xy_ui - logC2xy_ui; // cplx to be returned if no z

			max_info_logC = testinfo_xy_ui;
			min_info_logC = max_info_logC;
			/////////////////////////////////////
			for(j=0; j < nbrAllVar; j++ ) Opt_dBin[j] = dBin[j];
			if(nbrZi>0) {
			  for( k = 0; k <= nSample0; k++ ){ 
			    for( j=0; j< 6; j++ ){
			      Opt_sortedSample[k][j]=sortedSample[k][j];
			    }
			  }
			}
			////////////////////////////////////

		} else if(min_info_logC > testinfo_xy_ui){
			countmin++;
			min_info_logC = testinfo_xy_ui;
		} else {
			countmin=0;
			min_info_logC = testinfo_xy_ui;
		}

	///////////////////////////////////////////////////////////////////////
    if(nbrZi==0){

		ptrRetValues[0] = N_xy_ui; 
		ptrRetValues[1] = NIxy_ui;
		ptrRetValues[2] = k_xy_ui;

	// pick next z and compute score ///////////////////////////////
    } else { //(nbrZi>0)


		/////////////////////////////////////
	 	for( j=0; j < nbrAllVar; j++ ){
			dBin[j] = Opt_dBin[j];
		}

		for( k = 0; k <= nSample0; k++ ){ 
		  for( j=0; j< 6; j++ ){
		    sortedSample[k][j]=Opt_sortedSample[k][j];
		  }
		}
		////////////////////////////////////

		/// find optimum zi /////////////////////////////////////

		//N_xyuiz_top = -1; 
		N_xyuiz_top = 0; 
		NIxy_ui_top  = -1;
		k_xy_ui_top  = -1;
		
		z_top = -1;
		NIxy_uiz_top = -1;
		k_xy_uiz_top = -1;
		
		R_top= -LARGE;
		NIxyz_ui_top = -1;
		k_xyz_ui_top = -1;

		struct ScoresOfZ** scoresAllZi =new ScoresOfZ*[nbrZi];
		struct ContainerT cont;

		cont.sortedSample = sortedSample;
		cont.ptrAllData = ptrAllData;
		// cont.dBin = dBin;
		cont.ptrAllLevels = ptrAllLevels;
		cont.sampleSizeEff = sampleSizeEff;
		cont.nSample0 = nSample0;
		cont.sampleSize = sampleSize;
		cont.nbrAllVar = nbrAllVar;
		cont.nbrUi = nbrUi;
		cont.k23 = k23;
		// cont.Opt_dBin = Opt_dBin;
		cont.Nxuiz = Nxuiz;
		cont.Ny = Ny;
		cont.Nxui = Nxui;
		cont.Nx = Nx;
		cont.Nxyuiz = Nxyuiz;
		cont.Nyuiz = Nyuiz;
		cont.Nuiz = Nuiz;
		cont.Nz = Nz;
		cont.modCplx = modCplx;
		cont.c2terms = c2terms;
		cont.scoresAllZi = scoresAllZi;
		cont.weights = weights;
		cont.ptrZiIdx = ptrZiIdx;
		cont.weights = weights;
		cont.testDistribution = testDistribution;
		cont.freqs1 = freqs1;
		// cont.bridge=bridge;

		cont.sample = sample;
		cont.ptrVarIdx = ptrVarIdx;

		if(nthreadsMax == 0)
			nthreadsMax = 1;

		int step;
		step = nbrZi/(nthreadsMax);
		if(nbrZi%(nthreadsMax) != 0)
			step++;


		if(nbrZi%(step) == 0)
			nthreadsMax=nbrZi/(step);
		else 
			nthreadsMax=(nbrZi/(step) + 1);


		int iterator = 0;

		// if(environment.numSearchMore%(nthreadsMax-1) == 0)
		// 	nthreadsMax--;

		pthread_t* pt;
		pt = (pthread_t*)(malloc(nthreadsMax * sizeof(pthread_t)));
		struct FromToPlusContainer* fromTo = new FromToPlusContainer[nthreadsMax];

		int i = 0;

		while(i < nbrZi){

			if(i + step > nbrZi)
				step = nbrZi - i;
			
			fromTo[iterator].dBin = (int *)malloc(ncol* sizeof(int));
			fromTo[iterator].Opt_dBin = (int *)malloc(ncol* sizeof(int));
			fromTo[iterator].nbrLevCorrect = nbrLevCorrect;
			int kk,pp;
			for(kk = 0 ; kk < nbrAllVar + 1; kk++){
				fromTo[iterator].nbrLevCorrect[kk] = nbrLevCorrect[kk];
				fromTo[iterator].dBin[kk] = dBin[kk];
				fromTo[iterator].Opt_dBin[kk] = Opt_dBin[kk];
			}

			for(kk = 0 ; kk < (*memory).maxlevel + 1; kk++){

				// memoryThreads[iterator].Nxyuiz[kk] = Nxyuiz[kk];
				memoryThreads[iterator].Nyuiz[kk] = Nyuiz[kk];
				memoryThreads[iterator].Nuiz[kk] = Nuiz[kk];
				memoryThreads[iterator].Nz[kk] = Nz[kk];
				memoryThreads[iterator].Ny[kk] = Ny[kk];
				memoryThreads[iterator].Nxui[kk] = Nxui[kk];
				memoryThreads[iterator].Nx[kk] = Nx[kk];

			}

			for(kk = 0 ; kk < sampleSize + 1; kk++){
				for(pp = 0 ; pp < 7; pp++){
					memoryThreads[iterator].sortedSample[kk][pp] = sortedSample[kk][pp];
				}
			}

			for(kk = 0 ; kk < (*memory).maxlevel + 1; kk++){
				for(pp = 0 ; pp < (*memory).maxlevel + 1; pp++){
					memoryThreads[iterator].Nxuiz[kk][pp] = Nxuiz[kk][pp];
				}
			}

			fromTo[iterator].from = i;
			fromTo[iterator].to = i + step;
			fromTo[iterator].cont= &cont;
			fromTo[iterator].memoryThread= &memoryThreads[iterator];

			//// Search for new contributing node and its rank
			if (pthread_create(&pt[iterator], NULL, evaluateZ, (void *)&fromTo[iterator]) < 0)
			{ 
	  			exit (1);
	  		}

			iterator++;
			i += step;	
		}

		int pos;
		for(pos = 0; pos < nthreadsMax; pos++)
			pthread_join( pt[pos], NULL);

		

		for(zi=0; zi < nbrZi ; zi++){
			// set values of best z
			if(scoresAllZi[zi]->Rzi > R_top){ //
					N_xyuiz_top = scoresAllZi[zi]->N_xyuiz; 
					NIxy_ui_top  = scoresAllZi[zi]->info_xy_ui;
					k_xy_ui_top  = scoresAllZi[zi]->logC_xy_ui;

					z_top = zi;
					NIxy_uiz_top = scoresAllZi[zi]->info_xy_uiz;
					k_xy_uiz_top = scoresAllZi[zi]->logC_xy_uiz;

					R_top = scoresAllZi[zi]->Rzi;
					NIxyz_ui_top = scoresAllZi[zi]->info_xy_ui - scoresAllZi[zi]->info_xy_uiz;
					//k_xyz_ui_top = logC_xy_ui - logC_xy_uiz;
					k_xyz_ui_top = - scoresAllZi[zi]->logC_xy_ui + scoresAllZi[zi]->logC_xy_uiz; // to fit eq(20) and eq(22) in BMC Bioinfo 2016
				}
	      
		}

		ptrRetValues[0] = N_xyuiz_top; 
		ptrRetValues[1] = NIxy_ui_top;
		ptrRetValues[2] = k_xy_ui_top;

		ptrRetValues[3] = z_top;
		ptrRetValues[4] = NIxy_uiz_top;
		ptrRetValues[5] = k_xy_uiz_top;

		ptrRetValues[6] = R_top;
		ptrRetValues[7] = NIxyz_ui_top;
		ptrRetValues[8] = k_xyz_ui_top;

		iterator -=1;
		while(iterator >=0 ){
			free(fromTo[iterator].dBin);
			free(fromTo[iterator].Opt_dBin);
			delete [] scoresAllZi[iterator];
			iterator--;
		}
		delete [] fromTo;
		free(pt);
		delete [] scoresAllZi;
		}
	
	}


	free(dBin);

	free(Opt_dBin);
	free(nbrLevCorrect);

	return( ptrRetValues ); 

}