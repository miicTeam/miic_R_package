

//#define _MY_DEBUG_ 1
//#define _MY_DEBUG_NEW 1
#include <cmath>
#include <string.h>

#include <stdio.h>
#include <stdlib.h>
#include "memory.h"
#include <R.h>

#include "computeInfo.h"

#include <numeric>
#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include "utilities.h"
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
	struct ContainerMemory* memoryThread;
	int* Opt_dBin;
	int* dBin;
};



// modCplx = myCplx = 0 --> MDL, myCplx = 1 --> NML
// If nbrZi== 0, return nSample[0]     & nSample[0]*I(xy|{ui})      & k_xy_ui
// If nbrZi > 0, return nSample[z_top] & nSample[z_top]*I(xy|{ui})  & k_xy_ui
//                      z_top          & nSample[z_top]*I(xy|{ui}z) & k_xy_uiz
//                      R_top          & nSample[z_top]*I(xyz|{ui}) & k_xyz_ui
double* getAllInfoNEW( int* ptrAllData, int* ptrAllLevels, int* ptrVarIdx, int nbrUi, int* ptrZiIdx, int nbrZi, int ziPos, int sampleSize,
	int sampleSizeEff, int modCplx, int k23, double* c2terms, MemorySpace* memory)
{
	int randomrescaling=1;
	float r,rr;

	int bin_max=100,MDL=0, NML=1, TRUE=1, FALSE=0;
	int l,ok;
	int **sample,**sortedSample,**Opt_sortedSample; //[N+1][7]

	if(sampleSize < bin_max)
		bin_max= sampleSize;

	int iii;
	int nrow=sampleSize+1;
	int ncol=7;


	sample = (*memory).sample;

	sortedSample = (*memory).sortedSample;

	Opt_sortedSample = (*memory).Opt_sortedSample;


	int nSample0,*nSample,*orderSample,*sampleKey;//[N+1 elements: 0 to N]

	nSample = (int *)malloc(nbrZi* sizeof(int));
	orderSample = (*memory).orderSample;
	sampleKey = (*memory).sampleKey;


	int* bridge = (*memory).bridge;

	int bin,PBin,Prui,increment,NN,X,Y,Z;
	int ptrzi,zi,z;

	double NlogN,logN;

	int Lxyui,Lyui,Lui;
	int  Nxyui,  Nyui,  Nui;
	int  Nxyuis,  Nyuis,  Nuis;

	int NNxyui,NNxyuiz,NNxyuizl,Ntot;  //for rescaling NML change 20160228


	int *Nxyuiz, *Nyuiz, *Nuiz, *Nz;  //[Z]

	Nxyuiz = (*memory).Nxyuiz;
	Nyuiz = (*memory).Nyuiz;
	Nuiz = (*memory).Nuiz;
	Nz = (*memory).Nz;
	int Nzs,Nuizs,Nyuizs,Nxyuizs,Nxuizs;

	int  Nxyuizl, Nyuizl, Nuizl, Nzl; //[Z]

	int *Ny;
	Ny = (*memory).Ny;
	int  Nyj,Nys;                     //[Y]

	int *Nxui, *Nx;                   //[X]

	Nxui = (*memory).Nxui;
	Nx = (*memory).Nx;

	int  Nxuij, Nxj, Nxuis, Nxs;      //[X]

	int **Nxuiz;                      //[X][Z]

	nrow = bin_max+1;
	ncol = bin_max+1;

	Nxuiz = (*memory).Nxuiz;
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

	int N_xy_ui=0,N_xyuiz=0;
	double min_info_logC=LARGE, max_info_logC= - LARGE;

	int countmin,kz0;

	int z_top,N_xyuiz_top;
	double R_top;
	double NIxy_ui=-1.0,NIxy_ui_top,NIxy_uiz_top,NIxyz_ui_top;
	double k_xy_ui=-1.0,k_xy_ui_top,k_xy_uiz_top,k_xyz_ui_top;


	int i, j, k;	// for loops

	int nbrAllVar;
	// Output pointer
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

	nbrAllVar = (nbrUi+2);

	int **dBin,**Opt_dBin; //[1][Nb_ui+2]

	nrow = 0 + 1;
	ncol = nbrAllVar +1;

	dBin = (int **)malloc(nrow* sizeof(int*));
	for(iii = 0; iii < nrow; iii++)
	  dBin[iii] = (int *)malloc(ncol* sizeof(int));


	Opt_dBin = (int **)malloc(nrow* sizeof(int*));
	for(iii = 0; iii < nrow; iii++)
	  Opt_dBin[iii] = (int *)malloc(ncol* sizeof(int));

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



	if(nSample0 > 0){

		// initialisation of bin numbers for continuous variables in x,y, ui
		for(j=0; j < nbrAllVar; j++ )
		{

			dBin[0][j]=ptrAllLevels[ptrVarIdx[j]];

			Opt_dBin[0][j] = dBin[0][j];
		}


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

		sample[k][4] = 0; // X
		sample[k][5] = 0; // Y
		sample[k][6] = 0; // Z



		//////////////////////////////////////////////////

		/////////sort sample in increasing Lxyui stored in sample[k][1]
		//for( k = 0; k <= nSample0; k++ ){
		for( k = 1; k <=nSample0 +1; k++ ){
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


		//////////////////////////////////////////////////


		/////////initialization of counts and mutual infos & logCs
		Lxyui = sortedSample[0][1]; // min Lxyui
		Lyui  = sortedSample[0][2]; // min Lyui
		Lui   = sortedSample[0][3]; // min Lui

		Nxyui = 1;
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

			if(sortedSample[k][1] > Lxyui){

			        NNxyui=0;
			        if(sampleSizeEff!=sampleSize){
				  if(randomrescaling){
			  		GetRNGstate();
				    NNxyui=(int)floor(((double)Nxyui*sampleSizeEff)/sampleSize);
				    r=((double)Nxyui*sampleSizeEff)/sampleSize - NNxyui;
				    rr=(double)unif_rand() ;
				    PutRNGstate();
				    if (r > rr) NNxyui++;
				  }else{
				    NNxyui=(int)floor(0.5+((double)Nxyui*sampleSizeEff)/sampleSize);
				  }
				}else{
				  NNxyui=Nxyui;
				}

				Nui     += NNxyui;
				Nyui    += NNxyui;
				Nxui[X] += NNxyui;

				Nx[X] += NNxyui;	//4
				Ny[Y] += NNxyui;	//4
				Ntot +=  NNxyui;

				if( NNxyui > 0 ){
				  NlogN = NNxyui * std::log(static_cast<double>(NNxyui));//log(NNxyui);
				  info_xui_y += NlogN;
				  info_yui_x += NlogN;
				}
				Lxyui = sortedSample[k][1];
				Nxyuis += NNxyui;
				Nxyui = 1;

				if(k<nSample0) X=sortedSample[k][4];
				if(k<nSample0) Y=sortedSample[k][5];



				if(sortedSample[k][2] > Lyui){

				        if( Nyui > 0 ){
					  NlogN = Nyui * std::log(static_cast<double>(Nyui));//log(Nyui);
					  info_yui_x -= NlogN;
					  info_ui_y  += NlogN;

					  if(modCplx != MDL){
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
								NlogN = Nxuij * std::log(static_cast<double>(Nxuij));//log(Nxuij);
								info_xui_y -= NlogN;
								info_ui_x  += NlogN;
								if(modCplx != MDL){
								    logC_xui_y += compute_LogC_C2(Nxuij,dBin[0][1],c2terms);
								}

								Nxuis += Nxuij;
								Nxui[j]=0;
							}
						}

						if( Nui > 0 ){
						  NlogN = Nui * std::log(static_cast<double>(Nui));//log(Nui);
						  info_ui_y  -= NlogN;
						  info_ui_x  -= NlogN;
						  if(modCplx != MDL){
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
			  Nxyui++;
			}

		}

		// increment for info for Nx[X] and Ny[Y] contributions
		for( j = 0; j < dBin[0][0]; j++ ){
			Nxj=Nx[j];
			if( Nxj > 0 ){
			      //NlogN = Nxj * log(Nxj/(1.0*nSample0));
				NlogN = Nxj * std::log(static_cast<double>(Nxj/(1.0*Ntot)));//log(Nxj/(1.0*Ntot));
				info_yui_x -= NlogN;
				info_ui_x  -= NlogN;
				Nxs += Nxj;
			}
		}
		for( j = 0; j < dBin[0][1]; j++ ){
			Nyj=Ny[j];
			if( Nyj > 0 ){
			      //NlogN = Nyj * log(Nyj/(1.0*nSample0));
				NlogN = Nyj * std::log(static_cast<double>(Nyj/(1.0*Ntot)));//log(Nyj/(1.0*Ntot));
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


		if(modCplx == MDL) {
			Prui=1;
			//NN=nSample0;
			//NN = (int)floor(0.5 + ((double)nSample0*sampleSizeEff)/sampleSize);
			//logN=log(NN);
			logN=std::log(static_cast<double>(Ntot));//log(Ntot);
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
		

///////////////////////////////////////////////////////////////////////
    if(nbrZi==0){

	ptrRetValues[0] = N_xy_ui;
	ptrRetValues[1] = NIxy_ui;
	ptrRetValues[2] = k_xy_ui;

/////////pick next z and compute score ///////////////////////////////
    } else { //(nbrZi>0)

		/////////////////////////////////////
	 	for( j=0; j < nbrAllVar; j++ ){
			dBin[0][j] = Opt_dBin[0][j];
		}
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

			if(nSample[zi] == 1){

				max_info_logC = - LARGE;
				min_info_logC = LARGE;

				countmin=0; // change 20160216
					/////////initialization of counts and mutual infos & logCs


			    nSample[zi]=0;
				Lxyui = sortedSample[kz0][1]; // min Lxyui
				Lyui  = sortedSample[kz0][2]; // min Lyui
				Lui   = sortedSample[kz0][3]; // min Lui

				NNxyuiz = 0;
				NNxyuizl= 0;
				Nxyui = 0;
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
					Nxyuiz[l]=0;
					Nyuiz[l]=0;
					Nuiz[l]=0;
					Nz[l]=0;
				}

				X=sortedSample[kz0][4];
				Y=sortedSample[kz0][5];
				i=sortedSample[kz0][0];


				Z=ptrAllData[ i + ptrzi*sampleSize ];  // first Z

				Nxyuiz[Z]=1;


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
				for( k = kz0+1; k <= nSample0; k++ ){

					i=sortedSample[k][0];

					// check whether zi does not contain NA
					if( ptrAllData[ i + ptrzi*sampleSize ] > -1 )
					{

						Z=ptrAllData[ i + ptrzi*sampleSize ];  // Z

						if(sortedSample[k][1] > Lxyui){

						        NNxyuiz=0;
							for( l = 0; l < dBin[0][z]; l++ ){
								Nxyuizl=Nxyuiz[l];//
								if( Nxyuizl > 0 ){
								    if(sampleSizeEff!=sampleSize){
								        if(randomrescaling){
							        	GetRNGstate();
									    NNxyuizl=(int)floor(((double)Nxyuizl*sampleSizeEff)/sampleSize);
									    r=((double)Nxyuizl*sampleSizeEff)/sampleSize - NNxyuiz;
									    rr=(double)unif_rand() ;
									    PutRNGstate();
									    if (r > rr) NNxyuizl++;
									}else{
									    NNxyuizl=(int)floor(0.5+((double)Nxyuizl*sampleSizeEff)/sampleSize);
									}
								    }else{
								        NNxyuizl=Nxyuizl;
								    }

								    if( NNxyuizl > 0 ){
								        NlogN = NNxyuizl * std::log(static_cast<double>(NNxyuizl));//log(NNxyuizl);

									info_xuiz_y += NlogN;//
									info_yuiz_x += NlogN;//
									NNxyuiz += NNxyuizl;




									Nz[l]   += NNxyuizl;	//4
									Nuiz[l] += NNxyuizl;
									Nyuiz[l]+= NNxyuizl;
									Nxuiz[X][l] += NNxyuizl;

								    }

								    Nxyuiz[l]=0;//
								}
							}
							Nxyuiz[Z]=1;  //1  change


							if(NNxyuiz > 0){
							    NlogN = NNxyuiz * std::log(static_cast<double>(NNxyuiz));//log(NNxyuiz);
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
							        NlogN = Nyui * std::log(static_cast<double>(Nyui));//log(Nyui);
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
								        NlogN = Nyuizl * std::log(static_cast<double>(Nyuizl));//log(Nyuizl);
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
								        NlogN = Nui * std::log(static_cast<double>(Nui));//log(Nui);
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


											NlogN = Nuizl * std::log(static_cast<double>(Nuizl));//log(Nuizl);
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

									for( j = 0; j < dBin[0][0]; j++ ){
										Nxuij=Nxui[j];
										if( Nxuij > 0 ){
											NlogN = Nxuij * std::log(static_cast<double>(Nxuij));//log(Nxuij);
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
													NlogN = Nxuizjl * std::log(static_cast<double>(Nxuizjl));//log(Nxuizjl);
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
										}
		       							}

							         }
							         Lui = sortedSample[k][3];
							    }

							}

						} else {
							Nxyuiz[Z]++;	        //1
						}
					}
				}

		    	// increment info with Nx[X], Ny[Y] and Nz[Z] contributions
				for( j = 0; j < dBin[0][0]; j++ ){
					Nxj=Nx[j];
					if( Nxj > 0 ){
						NlogN = Nxj * std::log(static_cast<double>(Nxj/(1.0*nSample[zi])));//log(Nxj/(1.0*nSample[zi]));
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
						NlogN = Nyj * std::log(static_cast<double>(Nyj/(1.0*nSample[zi])));//log(Nyj/(1.0*nSample[zi]));
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
						NlogN = Nzl * std::log(static_cast<double>(Nzl/(1.0*nSample[zi])));//log(Nzl/(1.0*nSample[zi]));
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

					logN=std::log(static_cast<double>(nSample[zi]));//log(nSample[zi]);

					for(j=2; j < nbrAllVar; j++ )
						Prui *= dBin[0][j];
				}

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


				///////// compute score and store z with max score
				xz  = info_xz_ui - info_xy_ui;
				yz  = info_yz_ui - info_xy_ui;
				xyz = info_xy_ui - info_xy_uiz;
				if(k23==TRUE){
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

	}
	}



	for(i=0; i<1;i++)
		free(dBin[i]);
	free(dBin);

	for(i=0; i<1;i++)
		free(Opt_dBin[i]);
	free(Opt_dBin);

	free(nSample);


	return( ptrRetValues );

}



void* evaluateZ(void* container) {
	struct FromToPlusContainer* cont1 = (struct FromToPlusContainer*)container;
	struct ContainerMemory* memoryThread = cont1->memoryThread;
	
	struct ContainerT* cont = cont1->cont;

	int** sortedSample = memoryThread->sortedSample;
	int* ptrAllData = cont->ptrAllData;
	int* dBin = cont1->dBin;

	int* ptrAllLevels = cont->ptrAllLevels;
	int sampleSizeEff = cont->sampleSizeEff;
	int nSample0 = cont->nSample0;
	int ptrzi = cont->ptrzi;
	int sampleSize = cont->sampleSize;

	int nbrAllVar = cont->nbrAllVar;
	int nbrUi = cont->nbrUi;
	int k23 = cont->k23;
	int* Opt_dBin = cont1->Opt_dBin;

	int** Nxuiz = memoryThread->Nxuiz;
	int* Ny = memoryThread->Ny;
	int* Nxui = memoryThread->Nxui;
	int* Nx = memoryThread->Nx;
	int* Nxyuiz = memoryThread->Nxyuiz;
	int* Nyuiz = memoryThread->Nyuiz;
	int* Nuiz = memoryThread->Nuiz;
	int* Nz = memoryThread->Nz;

	int modCplx = cont->modCplx;
	int* ptrZiIdx = cont->ptrZiIdx;
	double* c2terms = cont->c2terms;
	struct ScoresOfZ** scoresAllZi = cont->scoresAllZi;

	
	int zpos = cont1->from;
	int to = cont1->to;

	int randomrescaling=1;
	float r,rr;
	int TRUE=1;

	int z;
	int k = 0, l = 0 ;
	double max_info_logC;
	double min_info_logC;

	int nSampleZ=0;

	int Lxyui,Lyui,Lui;
	int  Nyui,  Nui;
	int  Nxyuis,  Nyuis,  Nuis;

	int NNxyuiz,NNxyuizl,Ntot;  //for rescaling NML change 20160228

	int Prui,X,Y,Z,i,MDL=0;

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

	int  Nxyuizl, Nyuizl, Nuizl, Nzl; //[Z]
	int Nxuizjl;
	int j;

	double info3xy_ui,info2xy_ui,info3xy_uiz,info2xy_uiz;
	double info3xz_ui,info2xz_ui,info3yz_ui,info2yz_ui;
	double logC3xy_ui,logC2xy_ui,logC3xy_uiz,logC2xy_uiz;
	double logC3xz_ui,logC2xz_ui,logC3yz_ui,logC2yz_ui;


	int N_xy_ui=0,N_xyuiz=0;

	double info_xy_ui,info_xy_uiz;
	double info_xz_ui,info_yz_ui;
	double logC_xy_ui,logC_xy_uiz;
	double logC_xz_ui,logC_yz_ui;
	double testinfo_xy_ui,testinfo_xy_uiz;
	double testinfo_xz_ui,testinfo_yz_ui;
	double yz,xz,xyz,first,second,dpi,Rzi,testz;

	int countmin,stop,kz0;

	z=nbrAllVar;
	while(zpos < to){
		scoresAllZi[zpos] = new ScoresOfZ();

		ptrzi=ptrZiIdx[zpos];

		dBin[z]=ptrAllLevels[ptrzi];
		
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
				Nxyuiz[l]=0;
				Nyuiz[l]=0;
				Nuiz[l]=0;
				Nz[l]=0;
			}

			X=sortedSample[kz0][4];
			Y=sortedSample[kz0][5];
			i=sortedSample[kz0][0];
			Z=ptrAllData[ i + ptrzi*sampleSize ];  // first Z

			Nxyuiz[Z]=1;

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

							Nxyuizl=Nxyuiz[l];//
							if( Nxyuizl > 0 ){

							    if(sampleSizeEff!=sampleSize){
							        if(randomrescaling){
							        	GetRNGstate();
									    NNxyuizl=(int)floor(((double)Nxyuizl*sampleSizeEff)/sampleSize);
									    r=((double)Nxyuizl*sampleSizeEff)/sampleSize - NNxyuiz;
									    rr=(double)unif_rand() ;
									    PutRNGstate();
								    if (r > rr) NNxyuizl++;
								}else{
								    NNxyuizl=(int)floor(0.5+((double)Nxyuizl*sampleSizeEff)/sampleSize);
								}
							    }else{
							        NNxyuizl=Nxyuizl; 
							    }

							    if( NNxyuizl > 0 ){
							        NlogN = NNxyuizl * std::log(static_cast<double>(NNxyuizl));//log(NNxyuizl);
								info_xuiz_y += NlogN;//
								info_yuiz_x += NlogN;//
								NNxyuiz += NNxyuizl;

								Nz[l]   += NNxyuizl;	//4
								Nuiz[l] += NNxyuizl;
								Nyuiz[l]+= NNxyuizl;
								Nxuiz[X][l] += NNxyuizl;
							    }

							    Nxyuiz[l]=0;//
							}
						}
						Nxyuiz[Z]=1;  //1  change

						if(NNxyuiz > 0){


						    NlogN = NNxyuiz * std::log(static_cast<double>(NNxyuiz));//log(NNxyuiz);
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
						        NlogN = Nyui * std::log(static_cast<double>(Nyui));//log(Nyui);
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
							        NlogN = Nyuizl * std::log(static_cast<double>(Nyuizl));//log(Nyuizl);
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
							        NlogN = Nui * std::log(static_cast<double>(Nui));//log(Nui);
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
										NlogN = Nuizl * std::log(static_cast<double>(Nuizl));//log(Nuizl);
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
										NlogN = Nxuij * std::log(static_cast<double>(Nxuij));//log(Nxuij);
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
												NlogN = Nxuizjl * std::log(static_cast<double>(Nxuizjl));//log(Nxuizjl);
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
						Nxyuiz[Z]++;	        //1
					}
				}
			}


	    	// increment info with Nx[X], Ny[Y] and Nz[Z] contributions	
			for( j = 0; j < dBin[0]; j++ ){ 
				Nxj=Nx[j];
				if( Nxj > 0 ){
					NlogN = Nxj * std::log(static_cast<double>(nSampleZ));//log(Nxj/(1.0*nSampleZ));
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
					NlogN = Nyj * std::log(static_cast<double>(nSampleZ));//log(Nyj/(1.0*nSampleZ));
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
					NlogN = Nzl * std::log(static_cast<double>(nSampleZ));//log(Nzl/(1.0*nSampleZ));
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
				logN=std::log(static_cast<double>(nSampleZ));//log(nSampleZ);

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
			if(k23==TRUE){
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

		}
		zpos += 1;
	}
	void* dummyPointer = NULL;
	return dummyPointer;
}



// modCplx = myCplx = 0 --> MDL, myCplx = 1 --> NML
// If nbrZi== 0, return nSample[0]     & nSample[0]*I(xy|{ui})      & k_xy_ui 
// If nbrZi > 0, return nSample[z_top] & nSample[z_top]*I(xy|{ui})  & k_xy_ui 
//                      z_top          & nSample[z_top]*I(xy|{ui}z) & k_xy_uiz 
//                      R_top          & nSample[z_top]*I(xyz|{ui}) & k_xyz_ui 
double* getAllInfoNEWThreads( int* ptrAllData, int* ptrAllLevels, int* ptrVarIdx, int nbrUi, int* ptrZiIdx, int nbrZi, int ziPos, int sampleSize, 
	int sampleSizeEff, int modCplx, int k23, double* c2terms, MemorySpace* memory, int nthreadsMax, struct ContainerMemory* memoryThreads)
{	
	int randomrescaling=1;
	float r,rr;

	int bin_max=100,MDL=0, NML=1, TRUE=1, FALSE=0;
	int l,ok;
	int **sample,**sortedSample,**Opt_sortedSample; //[N+1][7]

	int nrow=sampleSize+1;
	int ncol=7;

	sample = (*memory).sample;

	sortedSample = (*memory).sortedSample;

	Opt_sortedSample = (*memory).Opt_sortedSample;

	int nSample0,*nSample,*orderSample,*sampleKey;//[N+1 elements: 0 to N]

	nSample = (int *)malloc(nbrZi* sizeof(int));
	orderSample = (*memory).orderSample;
	sampleKey = (*memory).sampleKey;

	int bin,PBin,Prui,increment,X,Y,Z;
	int ptrzi,zi,z;

	double NlogN,logN;

	int Lxyuiz,Lxyui,Lyui,Lui;
	int  Nxyui,  Nyui,  Nui;
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

	nrow = bin_max+1;
	ncol = bin_max+1;

	Nxuiz = (*memory).Nxuiz;
	int Nxuizjl;                     //[X][Z]

	double info_xui_y,info_yui_x,info_ui_y,info_ui_x;
	double logC_xui_y,logC_yui_x,logC_ui_y,logC_ui_x;

	double test_xy_ui;

	double info3xy_ui,info2xy_ui;
	double info3xz_ui,info2xz_ui;
	double logC3xy_ui,logC2xy_ui;
	double logC3xz_ui,logC2xz_ui;


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

	int countmin,stop;

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

	nrow = 0 + 1;
	ncol = nbrAllVar +1;
	
    dBin = (int *)malloc(ncol* sizeof(int));
	Opt_dBin = (int *)malloc(ncol* sizeof(int));

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

		Nxyui = 1;
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
				    NNxyui=(int)floor(((double)Nxyui*sampleSizeEff)/sampleSize);
				    r=((double)Nxyui*sampleSizeEff)/sampleSize - NNxyui;
				    rr=(double)unif_rand() ;
				    PutRNGstate();
				    if (r > rr) NNxyui++;
				  }else{
				    NNxyui=(int)floor(0.5+((double)Nxyui*sampleSizeEff)/sampleSize);
				  }
				}else{
				  NNxyui=Nxyui; 
				}
				
				Nui     += NNxyui; 
				Nyui    += NNxyui; 
				Nxui[X] += NNxyui; 
				
				Nx[X] += NNxyui;	//4
				Ny[Y] += NNxyui;	//4
				Ntot +=  NNxyui;

				if( NNxyui > 0 ){
				  NlogN = NNxyui * std::log(static_cast<double>(NNxyui));//log(NNxyui);
				  info_xui_y += NlogN;
				  info_yui_x += NlogN;
				}
				Lxyui = sortedSample[k][1];
				Nxyuis += NNxyui;
				Nxyui = 1;

				if(k<nSample0) X=sortedSample[k][4];
				if(k<nSample0) Y=sortedSample[k][5];
		

				if(sortedSample[k][2] > Lyui){

				        if( Nyui > 0 ){
					  NlogN = Nyui * std::log(static_cast<double>(Nyui));//log(Nyui);
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
								NlogN = Nxuij * std::log(static_cast<double>(Nxuij));//log(Nxuij);
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
						  NlogN = Nui * std::log(static_cast<double>(Nui));//(Nui);
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
			  Nxyui++;
			}
		}

		// increment for info for Nx[X] and Ny[Y] contributions	
		for( j = 0; j < dBin[0]; j++ ){ 
			Nxj=Nx[j];
			if( Nxj > 0 ){
			      //NlogN = Nxj * log(Nxj/(1.0*nSample0));
				NlogN = Nxj * std::log(static_cast<double>(Nxj/(1.0*Ntot)));//log(Nxj/(1.0*Ntot));
				info_yui_x -= NlogN;
				info_ui_x  -= NlogN;
				Nxs += Nxj;
			}
		}
		for( j = 0; j < dBin[1]; j++ ){ 
			Nyj=Ny[j];
			if( Nyj > 0 ){
			      //NlogN = Nyj * log(Nyj/(1.0*nSample0));
				NlogN = Nyj * std::log(static_cast<double>(Nyj/(1.0*Ntot)));//log(Nyj/(1.0*Ntot));
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

			logN=std::log(static_cast<double>(Ntot));//log(Ntot);
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

		int p;
 
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

		z=nbrAllVar;

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
		cont.ptrZiIdx = ptrZiIdx;
		// cont.bridge=bridge;

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
			int kk,pp;
			for(kk = 0 ; kk < nbrAllVar + 1; kk++){
				fromTo[iterator].dBin[kk] = dBin[kk];
				fromTo[iterator].Opt_dBin[kk] = Opt_dBin[kk];
			}

			for(kk = 0 ; kk < (*memory).maxlevel + 1; kk++){

				memoryThreads[iterator].Nxyuiz[kk] = Nxyuiz[kk];
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
	  			return NULL;
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

	return( ptrRetValues ); 

}
