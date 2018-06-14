#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <time.h>
#include <algorithm>
#include <iostream>

#include "utilities.h"

#include "computeInfo.h"

using namespace std;
//functions utilities for dynamic programming optimization


///////////////////////////////////////////////////////////////////////////////////////////

//#define _MY_DEBUG_NEW 1
//#define _MY_DEBUG_NEW_UI_2 1
//#define _MY_DEBUG_NEW_UI 1
//#define _MY_DEBUG_NEW_OPTFUN 1
//#define DEBUG_JOINT 1
//#define _MY_DEBUG_NEW_UI_NML 1
//#define _MY_DEBUG_NEW_3p 1
//#define _MY_DEBUG_NEW_3p_f 1
//#define DEBUG 1



/////////////////////////////////////////////////////////////////////////////
//INPUT:
//memory_cuts: vector length n with recorded recursvely best cuts, 
//possible values 0...n :
//0->stops (one bins); -k -> stops (two bin ([0 k-1][k ..];k->continute [.. k-1][k ...])
// 
//OUTPUT:
//r : number cuts
//cut: vector with cuts point-> [0 cut[0]][cut[0]+1 cut[1]]...[cut[r-2] cut[r-1]]

int reconstruction_cut_coarse(int *memory_cuts, int *memory_cuts2, int np, int n,int *cut)
{
		int ncuts=0;
		int l,s;

		if(memory_cuts[np-1]==0){

			//fprintf(stderr,"\n one bin!");// bin [0 n]
			//exit(1);
		}

		l=memory_cuts[np-1];

		while(l>0){
			ncuts++;
			l=memory_cuts[l-1];
		}
		if(l<0) ncuts++;

		cut[ncuts]=n-1;//conventional last cut 

		l=ncuts-1;
		s=memory_cuts[np-1];
		cut[l]=memory_cuts2[np-1];//truly last cut
		l--;
		while(s>0 && l>=0){
			cut[l]=memory_cuts2[s-1];
			s=memory_cuts[s-1];
			l--;
		}
		
		return ncuts+1; //number of levels (r)
		
}

/////////////////////////////////////////////////////////////////
//update datafactors of a variable from the cut positions vector <cut>
//INPUT:
//d: index of variable in datafactors
//varidx: index of variable in sortidx

void update_datafactors(int **sortidx, int varidx, int **datafactors,int d, int n, int **cut)
{
	int j,sjj,uu=0;

	for(j=0;j<=n-1;j++){
			sjj=sortidx[varidx][j];
			if(j>cut[d][uu]) uu++;
			datafactors[d][sjj]=uu;
	}
	return;
}



// /////////////////////////////////////////////////////////////////
// //update datafactors d

// void update_datafactors(int **sortidx, int **data, int varidx, int **datafactors,int d, int n, int **cut)
// {
// 	int j,sjj=sortidx[varidx][0];uu=0,rj=0;
// 	datafactors[d][sjj]=0;
// 	for(j=1;j<=n-1;j++){
// 		if(data[varidx][sjj]!=data[varidx][sortidx[varidx][j]]) rj++;;//account for repeated values
// 		sjj=sortidx[varidx][j];
// 		if(rj>cut[d][uu]) uu++;
// 		datafactors[d][sjj]=uu;
// 	}
	
// }


/////////////////////////////////////////////////////////////////////////////////////////

//INPUT:
//datarank, datafactors
//take 0 1 as x y and 2..<Mui>+1 for {u} indexes
//ignore index <dui> in the u

//OUTPUT
//return joint datafactors ui and uiy uix uixy, with number of levels ruiyx[0,1,2,3]
//ruiyx -> 0:u,1:uy,2:ux,3:uyx

void jointfactors_uiyx(int **datafactors,int dui, int n, int Mui, int *r, int **uiyxfactors,int *ruiyx)
{

		int df,Pbin_ui;
		int *datauix=(int*)calloc(n,sizeof(int));
		int *datauiy=(int*)calloc(n,sizeof(int));
		int *datauiyx=(int*)calloc(n,sizeof(int));
		int *dataui=(int*)calloc(n,sizeof(int));

		int j,jj,l;
		int uu=0;
		bool tooManyLevels = false;

		/////////////////
		//from single datafactors to joint datafactors ui; (ui,y);(ui,x);(ui,y,x)
		int nbrLevelsJoint = 1;
		for(jj=0;jj<Mui+2 && !tooManyLevels;jj++){
			nbrLevelsJoint *= r[jj];
			if(nbrLevelsJoint>8*n){
				tooManyLevels = true;
				// cout << "tooManyLevels: " << nbrLevelsJoint << endl;
			}
		}

		// if(!tooManyLevels)
		// 	cout << "NOT " << nbrLevelsJoint << endl;
		
		int *vecZeroOnesui;
		int *vecZeroOnesuiy;
		int *vecZeroOnesuix;
		int *vecZeroOnesuiyx;
		//
		int * orderSample_ux;
		int * orderSample_uyx;
		//declaration
		if(!tooManyLevels){
			vecZeroOnesui=(int*)calloc(nbrLevelsJoint,sizeof(int));
			vecZeroOnesuiy=(int*)calloc(nbrLevelsJoint,sizeof(int));
			vecZeroOnesuix=(int*)calloc(nbrLevelsJoint,sizeof(int));
			vecZeroOnesuiyx=(int*)calloc(nbrLevelsJoint,sizeof(int));
		}
		else{
			orderSample_ux = new int[n];
			orderSample_uyx = new int[n];
		}

		for(jj=0;jj<=n-1;jj++){

			datauiyx[jj]=datafactors[0][jj];//yx
			datauiyx[jj]+=datafactors[1][jj]*r[0];
			datauix[jj]=datafactors[0][jj];//x
			datauiy[jj]=datafactors[1][jj];//y
			dataui[jj]=0;

			Pbin_ui=1;
			for (l=Mui+1;l>=2;l--){
				if(l!=dui){
					df=datafactors[l][jj]*Pbin_ui;
					datauiyx[jj]+=df*r[1]*r[0];
					datauix[jj]+=df*r[0];
					datauiy[jj]+=df*r[1];


					dataui[jj]+=df;
					Pbin_ui*=r[l];
				}
			}
			//init

			if(!tooManyLevels){
				vecZeroOnesui[dataui[jj]] = 1;
				vecZeroOnesuiy[datauiy[jj]] = 1;
				vecZeroOnesuix[datauix[jj]] = 1;
				vecZeroOnesuiyx[datauiyx[jj]] = 1;		
			}
			else{
				orderSample_ux[jj] = jj;
				orderSample_uyx[jj] = jj;
			}
		}

		if(!tooManyLevels){

			// create the vector for storing the position
			int pos = 0;
			int pos1 = 0;
			int pos2 = 0;
			int pos3 = 0;
			for(jj=0;jj<=nbrLevelsJoint;jj++){
				
				if(vecZeroOnesui[jj] == 1){
					vecZeroOnesui[jj] = pos;
					pos += 1;
				}

				if(vecZeroOnesuiy[jj] == 1){
					vecZeroOnesuiy[jj] = pos1;
					pos1 += 1;
				}

				if(vecZeroOnesuix[jj] == 1){
					vecZeroOnesuix[jj] = pos2;
					pos2 += 1;
				}

				if(vecZeroOnesuiyx[jj] == 1){
					vecZeroOnesuiyx[jj] = pos3;
					pos3 += 1;
				}

			}

			ruiyx[0]=pos;//ui
			ruiyx[1]=pos1;//uiy
			ruiyx[2]=pos2;//uix
			ruiyx[3]=pos3;//uiyx


			for(j=0;j<=n-1;j++){
				uiyxfactors[0][j]=vecZeroOnesui[dataui[j]];//ui
				uiyxfactors[1][j]=vecZeroOnesuiy[datauiy[j]];;//uiy
				uiyxfactors[2][j]=vecZeroOnesuix[datauix[j]];;//uix
				uiyxfactors[3][j]=vecZeroOnesuiyx[datauiyx[j]];;//uiyx

			}
		} else {

			
			sort2arraysConfidence(n,datauix,orderSample_ux);
			sort2arraysConfidence(n,datauiyx,orderSample_uyx);

			//joint datafactors without gaps
			//compute term constant H(y,Ui) and H(Ui)
			
			int ix,iyx,iix,iiyx;
			
			ruiyx[0]=0;//ui
			ruiyx[1]=0;//uiy
			ruiyx[2]=0;//uix
			ruiyx[3]=0;//uiyx

			ix=orderSample_ux[0];
			iyx=orderSample_uyx[0];
			uiyxfactors[0][ix]=0;//ui
			uiyxfactors[1][iyx]=0;//uiy
			uiyxfactors[2][ix]=0;//uix
			uiyxfactors[3][iyx]=0;//uiyx

			for(j=0;j < n-1;j++){
				iix=ix;
				iiyx=iyx;
				ix=orderSample_ux[j+1];
				iyx=orderSample_uyx[j+1];

				if(dataui[ix]>dataui[iix]){
					ruiyx[0]++;
				}
				uiyxfactors[0][ix]=ruiyx[0];//ui

				if(datauiy[iyx]>datauiy[iiyx]){
					ruiyx[1]++;
				}
				uiyxfactors[1][iyx]=ruiyx[1];//uiy

				if(datauix[ix]>datauix[iix]){
					ruiyx[2]++;
				}
				uiyxfactors[2][ix]=ruiyx[2];//uix

				if(datauiyx[iyx]>datauiyx[iiyx]){
					ruiyx[3]++;
				}
				uiyxfactors[3][iyx]=ruiyx[3];//uiyx
			}	


			//number joint levels 
			ruiyx[0]++;//ui
			ruiyx[1]++;//uiy
			ruiyx[2]++;//uix
			ruiyx[3]++;//uiyx
			delete [] orderSample_ux;
			delete [] orderSample_uyx;

		}



		#if _MY_DEBUG_NEW_UI
			printf("j,dataui[j],datauiy[j]\n");
			for(j=0;j<=n-1;j++){
				printf("%d-> %d %d\n",j,dataui[j],datauiy[j]);
			}
		#endif

		#if _MY_DEBUG_NEW_UI_2
			printf("j,uiyxfactors[0][j],uiyxfactors[1][j],uiyxfactors[2][j],uiyxfactors[3][j]\n");
			for(j=0;j<=n-1;j++){
				printf("%d-> %d %d %d %d\n",j,uiyxfactors[0][j],uiyxfactors[1][j],uiyxfactors[2][j],uiyxfactors[3][j]);
			}
		#endif


		free(datauix);
		free(datauiy);
		free(datauiyx);
		free(dataui);
		if(!tooManyLevels){
			free(vecZeroOnesui);
			free(vecZeroOnesuiy);
			free(vecZeroOnesuix);
			free(vecZeroOnesuiyx);
		}

		return;
}
////////////////////////////////////////////////////////////////////////////////////////////////

















/////////////////////////////////////////////////////////////////////////////////////////

//INPUT:
//datarank, datafactors, cut
// ptrVar Index_x, Index_y, Index_u1,..Index_uMui

//OUTPUT
//return joint datafactors ui and uiy uix uixy, with number of levels ruiyx[0,1,2,3]

void jointfactors_uyx(int **datafactors, int* ptrVar, int n, int Mui, int *r, int **uyxfactors,int *ruyx)
{


		int df,Pbin_ui;
		int *datauix=(int*)calloc(n,sizeof(int));
		int *datauiy=(int*)calloc(n,sizeof(int));
		int *datauiyx=(int*)calloc(n,sizeof(int));
		int *dataui=(int*)calloc(n,sizeof(int));

		int j,jj,l;
		int uu=0;
		bool tooManyLevels = false;

		/////////////////
		//from single datafactors to joint datafactors ui; (ui,y);(ui,x);(ui,y,x)
		int nbrLevelsJoint = 1;
		for(jj=0; jj < Mui+ 2 && !tooManyLevels;jj++){
			nbrLevelsJoint *= r[ptrVar[jj]];
			if(nbrLevelsJoint>8*n)
				tooManyLevels = true;
		}

		//decl
		int *vecZeroOnesui;
		int *vecZeroOnesuiy;
		int *vecZeroOnesuix;
		int *vecZeroOnesuiyx;
		//
		int * orderSample_ux;
		int * orderSample_uyx;

		if(!tooManyLevels){
			vecZeroOnesui=(int*)calloc(nbrLevelsJoint,sizeof(int));
			vecZeroOnesuiy=(int*)calloc(nbrLevelsJoint,sizeof(int));
			vecZeroOnesuix=(int*)calloc(nbrLevelsJoint,sizeof(int));
			vecZeroOnesuiyx=(int*)calloc(nbrLevelsJoint,sizeof(int));
		}
		else{
			orderSample_ux = new int[n];
			orderSample_uyx = new int[n];
		}

		for(jj=0;jj<=n-1;jj++){

			datauiyx[jj]=datafactors[ptrVar[0]][jj];//yx
			datauiyx[jj]+=datafactors[ptrVar[1]][jj]*r[ptrVar[0]];
			datauix[jj]=datafactors[ptrVar[0]][jj];//x
			datauiy[jj]=datafactors[ptrVar[1]][jj];//y
			dataui[jj]=0;

			Pbin_ui=1;
			for (l=Mui+1;l>=2;l--){
					df=datafactors[ptrVar[l]][jj]*Pbin_ui;
					datauiyx[jj]+=df*r[ptrVar[1]]*r[ptrVar[0]];
					datauix[jj]+=df*r[ptrVar[0]];
					datauiy[jj]+=df*r[ptrVar[1]];

					dataui[jj]+=df;
					Pbin_ui*=r[ptrVar[l]];
			}
			//init

			if(!tooManyLevels){
				vecZeroOnesui[dataui[jj]] = 1;
				vecZeroOnesuiy[datauiy[jj]] = 1;
				vecZeroOnesuix[datauix[jj]] = 1;
				vecZeroOnesuiyx[datauiyx[jj]] = 1;		
			}
			else{
				orderSample_ux[jj] = jj;
				orderSample_uyx[jj] = jj;
			}
		
		}

		if(!tooManyLevels){
			// create the vector for storing the position
			int pos = 0;
			int pos1 = 0;
			int pos2 = 0;
			int pos3 = 0;
			for(jj=0;jj<=nbrLevelsJoint;jj++){
				
				if(vecZeroOnesui[jj] == 1){
					vecZeroOnesui[jj] = pos;
					pos += 1;
				}

				if(vecZeroOnesuiy[jj] == 1){
					vecZeroOnesuiy[jj] = pos1;
					pos1 += 1;
				}

				if(vecZeroOnesuix[jj] == 1){
					vecZeroOnesuix[jj] = pos2;
					pos2 += 1;
				}

				if(vecZeroOnesuiyx[jj] == 1){
					vecZeroOnesuiyx[jj] = pos3;
					pos3 += 1;
				}

			}


			ruyx[0]=pos;//ui
			ruyx[1]=pos1;//uiy
			ruyx[2]=pos2;//uix
			ruyx[3]=pos3;//uiyx


			for(j=0;j<=n-1;j++){
				uyxfactors[0][j]=vecZeroOnesui[dataui[j]];//ui
				uyxfactors[1][j]=vecZeroOnesuiy[datauiy[j]];;//uiy
				uyxfactors[2][j]=vecZeroOnesuix[datauix[j]];;//uix
				uyxfactors[3][j]=vecZeroOnesuiyx[datauiyx[j]];;//uiyx
			}
		} else {
			
			sort2arraysConfidence(n,datauix,orderSample_ux);
			sort2arraysConfidence(n,datauiyx,orderSample_uyx);

			//joint datafactors without gaps
			//compute term constant H(y,Ui) and H(Ui)
			
			int ix,iyx,iix,iiyx;
			
			ruyx[0]=0;//ui
			ruyx[1]=0;//uiy
			ruyx[2]=0;//uix
			ruyx[3]=0;//uiyx

			ix=orderSample_ux[0];
			iyx=orderSample_uyx[0];
			uyxfactors[0][ix]=0;//ui
			uyxfactors[1][iyx]=0;//uiy
			uyxfactors[2][ix]=0;//uix
			uyxfactors[3][iyx]=0;//uiyx

			for(j=0;j<n-1;j++){
				iix=ix;
				iiyx=iyx;
				ix=orderSample_ux[j+1];
				iyx=orderSample_uyx[j+1];

				if(dataui[ix]>dataui[iix]){
					ruyx[0]++;
				}
				uyxfactors[0][ix]=ruyx[0];//ui

				if(datauiy[iyx]>datauiy[iiyx]){
					ruyx[1]++;
				}
				uyxfactors[1][iyx]=ruyx[1];//uiy

				if(datauix[ix]>datauix[iix]){
					ruyx[2]++;
				}
				uyxfactors[2][ix]=ruyx[2];//uix

				if(datauiyx[iyx]>datauiyx[iiyx]){
					ruyx[3]++;
				}
				uyxfactors[3][iyx]=ruyx[3];//uiyx
			}	


			//number joint levels 
			ruyx[0]++;//ui
			ruyx[1]++;//uiy
			ruyx[2]++;//uix
			ruyx[3]++;//uiyx

			delete [] orderSample_ux;
			delete [] orderSample_uyx;

		}

		#if _MY_DEBUG_NEW_UI
			printf("j,dataui[j],datauiy[j]\n");
			for(j=0;j<=n-1;j++){
				printf("%d-> %d %d\n",j,dataui[j],datauiy[j]);
			}
		#endif

	
		#if _MY_DEBUG_NEW_UI_2
			printf("j,uyxfactors[0][j],uyxfactors[1][j],uyxfactors[2][j],uyxfactors[3][j]\n");
			for(j=0;j<=n-1;j++){
				printf("%d-> %d %d %d %d\n",j,uyxfactors[0][j],uyxfactors[1][j],uyxfactors[2][j],uyxfactors[3][j]);
			}
		#endif

		// free memory

		free(datauix);
		free(datauiy);
		free(datauiyx);
		free(dataui);
		if(!tooManyLevels){
			free(vecZeroOnesui);
			free(vecZeroOnesuiy);
			free(vecZeroOnesuix);
			free(vecZeroOnesuiyx);
		}

		return;
}




//////////////////////////////////////////////////////////////////

//INPUT:
//datarank, datafactors, cut

//OUTPUT
//return joint datafactors ui , with number of levels rui
//entropy term Hui
void jointfactors_u(int **datafactors,int *ptrIdx,int n, int Mui, int *r, int *ufactors,int *ru,double *Hu)
{

		//from cuts to joint datafactors 

		int jj;
		int i,j,l;

		if(Mui==1){
			for(jj=0;jj<=n-1;jj++){
				ufactors[jj]=datafactors[ptrIdx[0]][jj];
			}
			*ru=r[ptrIdx[0]];
			*Hu=0; //obsolete
			return;
		}
		/////////////////
		//update joint datafactors (with gaps) ui and (ui,y)

		int df,Pbin_ui;
		int *datau=(int*)calloc(n,sizeof(int));

		bool tooManyLevels = false;

		int nbrLevelsJoint = 1;
		for(l=0;l<Mui && !tooManyLevels;l++){
			nbrLevelsJoint *= r[ptrIdx[l]];
			if(nbrLevelsJoint>8*n)
				tooManyLevels = true;
		}

		//decl
		int *vecZeroOnesui;
		int * orderSample_u;
		if(!tooManyLevels){
			vecZeroOnesui=(int*)calloc(nbrLevelsJoint,sizeof(int));
		}
		else{
			orderSample_u= new int[n]; 
		}	

		for(jj=0;jj<=n-1;jj++){

			datau[jj]=0;

			Pbin_ui=1;
			for (l=Mui-1;l>=0;l--){
					df=datafactors[ptrIdx[l]][jj]*Pbin_ui;
					datau[jj]+=df;
					Pbin_ui*=r[ptrIdx[l]];
			}

			//init
			if(!tooManyLevels){
				vecZeroOnesui[datau[jj]] = 1;
			}
			else{
				orderSample_u[jj] = jj;
			}
		}

		if(!tooManyLevels){
			// create the vector for storing the position
			int pos = 0;
			
			for(jj=0;jj<=nbrLevelsJoint;jj++){
				
				if(vecZeroOnesui[jj] == 1){
					vecZeroOnesui[jj] = pos;
					pos += 1;
				}

			}

			*ru=pos;//ui

			for(j=0;j<=n-1;j++){
				ufactors[j]=vecZeroOnesui[datau[j]];//ui
			}

		 } else {
		
			sort2arraysConfidence(n,datau,orderSample_u);

			int ix, iix;
			ix=orderSample_u[0];
			*ru=0;
			ufactors[ix]=0;//ui



			for(j=0;j < n-1;j++){
				iix=ix;
				ix=orderSample_u[j+1];

				if(datau[ix]>datau[iix]){
					*ru=*ru+1;
				}
				ufactors[ix]=*ru;//ui

				// cout << ufactors[j] << " ";
			}	


			*ru=*ru+1;

			// cout << "\n" << *ru << "\n";
			delete [] orderSample_u;
		}

		#if DEBUG_JOINT
			printf("\nj -> datau[j]\n");
			for(j=0;j<=n-1;j++){
				printf("%d -> %d \n",j,datau[j]);
			}
			fflush(stdout);
		#endif


		#if DEBUG_JOINT
			printf("j,orderSample[j+1],sampleKey[j+1],datau[orderSample[j+1]]\n");
			for(j=0;j<=n-1;j++){
				printf("%d: %d %d, %d \n",j,orderSample[j+1],sampleKey[j+1],datau[orderSample[j+1]]);
			}
			fflush(stdout);
		#endif
	
		//joint datafactors without gaps
		//compute term constant H(y,Ui) and H(Ui)
	
		*Hu=0;//useless
		
		#if DEBUG
			printf("Hu=%lf\n",*Hu);
		#endif

		//number joint levels ui and uiy
		// *ru=(*ru)+1;

		#if DEBUG_JOINT
			printf("j,uiyfactors[0][i] (ru=%d)\n",*ru);
			for(j=0;j<=n-1;j++){
				printf("%d-> %d\n",j,ufactors[j]);
			}
			fflush(stdout);
		#endif

		free(datau);
		if(!tooManyLevels){
			free(vecZeroOnesui);
		}

		return;
}


//////////////////////////////////////////////////
// !! included in Info_cnt.cpp !!


// ///////////////////////////////////////////////
// //complicated (to handle) but flexible function 
// //optimize linear function depending on different variables

// //////////////////////////////////////////////
// //INPUT

// //optimizing on <sortidx_var>

// //<nbrV> number of factors variable   

// //<factors> : ( <nbrV> x n )  matrix :
// //example <nbrV>=6
// //x factors : factors[0] 
// //y factors : factors[1]
// //u factors : factors[2]
// //xu factors : factors[3]
// //yu factors : factors[4]
// //xyu factors : factors[5]

// //<r> : <nbrV> component vectors: levels in the <nbrV> components of factors

// //<optfun> : <nbrV> + 1  component vectors:
// //optimized function f=sum_k optfun[k]*H[k] 
// //f= optfun[nbrV+1]*H_optvar +optfun[0]*Hx + optfun[1]*Hy + optfun[2]*Hu + optfun[3]*Hxu + optfun[4]*Hyu + optfun[5]*Hxyu
// //optfun[nbrV+1]=reserved for optimized variable only

// //USE EXAMPLE

// //I(x;yu)=Hx+Hyu-Hxyu :
// //			optimize f on x : Hx-Hxyu : 
// //				factors[0]=0 //factors x= 0
// //				factors[3]=factors[2] //factors xu= factors u
// //				factors[5]=factors[4]	//factors xyu= factors yu
// //				opfun=[1,0,0,0,0,-1]
// //			optimize f on x : Hx-Hxyu : 
// //				factors[0]=0 //factors x= 0
// //				factors[3]=factors[2] //factors xu= factors u
// //				factors[5]=factors[4]	//factors xyu= factors yu
// //				opfun=[1,0,0,0,0,-1]

// //I(xu;y)=Hy+Hxu-Hxyu : opfun=[0,1,0,1,0,-1]


// //I(x;u)=Hx+Hu-Hxu : opfun=[1,0,1,-1,0,0]




// double optfun_onerun_kmdl_coarse(int *sortidx_var,int nbrV, int **factors, int *r, int *optfun, double sc,int n, int *memory_cut, int coarse, double* looklog, map<string,double> &looklbc)
// {


// 		int i,j,k,m;

// 		int np=ceil(1.0*n/coarse);

// 		double k_sc=0;
// 		double sc2,scr,sctemp;

	
// 		//dynamic programming optimize function and memorize of cuts
// 		int c,ct;
// 		double fmax;//I-kmdl
// 		int nc[np],nctemp;//Hx-Hxy-LogCx

// 		// double *H_kj=calloc(nbrV+1,sizeof(double));// x y u xu yu xyu
// 		double H_kj[nbrV+1];//x y u xu yu xyu

// 		// double *I=calloc(np,sizeof(double));
// 		// double *I_0k=calloc(np,sizeof(double));
// 		double I[np];
// 		I[0]=0;
// 		double I_0k[np];
// 		double I_kj,t;
		
// 		int nxj,nx;

// 		int xyu;



// 		//dynamic programming steps terms

// 		//H_0k needs initialization at zero (problem: double H_0k[nbrV+1][np];)
// 		double **H_0k=calloc(nbrV+1,sizeof(double*));// x y u xu yu xyu
// 		for(m=0;m<nbrV+1;m++) H_0k[m]=calloc(np,sizeof(double));

// 		int **nxyu=calloc(nbrV,sizeof(int*));// x y u xu yu xyu
// 		for(m=0;m<nbrV;m++) nxyu[m]=calloc(r[m],sizeof(int));

// 		int **nxyu_k=calloc(nbrV,sizeof(int*));// x y u xu yu xyu
// 		for(m=0;m<nbrV;m++) nxyu_k[m]=calloc(r[m],sizeof(int));	
		

// 		///////////////////////////////////////////////
// 		#if _MY_DEBUG_NEW_OPTFUN
// 			for(m=0;m<nbrV;m++){
// 				printf("r[%d]=%d :\n",m,r[m]);

// 				for(i=0;i<n;i++){
// 					printf("%d ",factors[m][i]);
// 				}
// 				printf("\n");

// 			}
// 		#endif


// 		/////////////////////////////////////////////////
// 		//j=0;
// 		#if _MY_DEBUG_NEW_OPTFUN
// 			printf("j=%d\n",0);fflush(stdout);
// 		#endif


// 		for(i=0;i<coarse;i++){

// 			for(m=0;m<nbrV;m++){

// 				if(optfun[m] != 0){ //compute only necessary terms

// 					xyu=factors[m][sortidx_var[i]];

// 					nxyu[m][xyu]++;

// 					if(nxyu[m][xyu] != 1) H_0k[m][0]-=nxyu[m][xyu]*looklog[nxyu[m][xyu]]-(nxyu[m][xyu]-1)*looklog[nxyu[m][xyu]-1];

// 				}

// 			}
// 		}

// 		if(optfun[nbrV] != 0) H_0k[nbrV][0]=-coarse*looklog[coarse];

// 		for(m=0;m<nbrV+1;m++) if(optfun[m] != 0)	I[0]+=optfun[m]*H_0k[m][0];
		
// 		memory_cut[0]=0;	
// 		I_0k[0]=I[0];
// 		nc[0]=1;

// 		#if _MY_DEBUG_NEW_OPTFUN
// 				printf("  ");
// 				for(m=0;m<nbrV+1;m++) printf("H_0k[%d]=%lf ",m,H_0k[m][0]);fflush(stdout);
// 				printf("\n  [0 -0] = %lf\n",I_0k[0]);fflush(stdout);
// 		#endif

		

// 		for(j=1;j<=np-1;j++){ //j=1...n-1
			
// 			#if _MY_DEBUG_NEW_OPTFUN
// 					printf("j=%d\n",j);fflush(stdout);
// 			#endif 

// 			for(m=0;m<nbrV+1;m++) if(optfun[m] != 0) H_0k[m][j]=H_0k[m][j-1];	
		
// 			ct=0;	
// 			for(i=(j*coarse);(i<(j+1)*coarse)&&(i<n);i++){
			

// 				for(m=0;m<nbrV;m++){

// 					if(optfun[m] != 0){ //compute only necessary terms

// 						xyu=factors[m][sortidx_var[i]];

// 						nxyu[m][xyu]++;

// 						if(nxyu[m][xyu] != 1) H_0k[m][j]-=nxyu[m][xyu]*looklog[nxyu[m][xyu]]-(nxyu[m][xyu]-1)*looklog[nxyu[m][xyu]-1];

// 					}

// 				}
	
// 				ct++;
// 			}
// 			#if _MY_DEBUG_NEW_OPTFUN
// 				printf("  ");
// 				for(m=0;m<nbrV;m++) {
// 					for(xyu=0;xyu<r[m];xyu++) printf("nxyu[%d][%d]=%d ",m,xyu,nxyu[m][xyu]);fflush(stdout);
// 				}
// 				printf("\n");
// 			#endif


// 			nxj=j*coarse+ct;
// 			if(optfun[nbrV] != 0) H_0k[nbrV][j]=-nxj*looklog[nxj];

// 			I_0k[j]=0;
// 			for(m=0;m<nbrV+1;m++) if(optfun[m] != 0)	I_0k[j]+=optfun[m]*H_0k[m][j];

// 			//complexity terms (local version)

// 			k_sc=sc*looklog[(nxj+1)];

// 			// 2 bins combinatorial term

// 			string=strcat(to_string(nxj),"-",to_string(1));
// 			map<string,double>::iterator it=looklbc.find(string);
// 			if( it != looklbc.end()){
// 				sc2=2*k_sc+ it->second;
// 			}
// 			else{
// 				double v=log(dbico(nxj,1));
// 				looklbc[string]=v;
// 				sc2=2*k_sc+ v;
// 			}

// 			// if(looklbc[nxj][1]==-1) looklbc[nxj][1]=log(dbico(nxj,1));
// 			// sc2=2*k_sc+looklbc[nxj][1];

// 			// sc2=2*k_sc+log(dbico(nxj,1));
			
// 			//init

// 			fmax=-DBL_MAX;

// 			for(m=0;m<nbrV+1;m++) if(optfun[m] != 0) H_kj[m]=H_0k[m][j];	

// 			#if _MY_DEBUG_NEW_OPTFUN
// 				printf("  ");
// 				//for(m=0;m<nbrV+1;m++) printf("H_kj[%d]=%lf ",m,H_kj[m]);fflush(stdout);
// 				printf("\n  [0 -%d] = %lf\n",j,I_0k[j]);fflush(stdout);
// 			#endif
					
// 			for(m=0;m<nbrV;m++) if(optfun[m] != 0) {
// 				for(xyu=0;xyu<r[m];xyu++) nxyu_k[m][xyu]=nxyu[m][xyu];
// 			}
		
// 			for(k=0;k<=j-1;k++){//k=1...n-2 possible cuts


// 					#if _MY_DEBUG_NEW_OPTFUN
// 						printf("k=%d  ",k);fflush(stdout);
// 					#endif

// 						for(i=0;i<coarse;i++){

// 								for(m=0;m<nbrV;m++){

// 									if(optfun[m] != 0){ //compute only necessary terms

// 										xyu=factors[m][sortidx_var[(k*coarse)+i]];
// 										nxyu_k[m][xyu]--;

// 										if(nxyu_k[m][xyu] != 0) H_kj[m]-=nxyu_k[m][xyu]*looklog[nxyu_k[m][xyu]]-(nxyu_k[m][xyu]+1)*looklog[(nxyu_k[m][xyu]+1)];
// 									}

// 								}

// 						}
						
// 						#if _MY_DEBUG_NEW_OPTFUN
// 							printf("\n  ");
// 							for(m=0;m<nbrV;m++) {
// 							for(xyu=0;xyu<r[m];xyu++) printf("nxyu_k[%d][%d]=%d ",m,xyu,nxyu_k[m][xyu]);fflush(stdout);
// 							}
// 							printf("\n");
// 						#endif

// 						nx=nxj-(k+1)*coarse;
// 						if(optfun[nbrV] != 0) H_kj[nbrV]=-nx*looklog[nx];

// 						I_kj=0;
// 						for(m=0;m<nbrV+1;m++) if(optfun[m] != 0)	I_kj+=optfun[m]*H_kj[m];

// 						//recursive number of bins combinatorial term
// 						string=strcat(to_string(nxj),"-",to_string(nc[k]));
// 						map<string,double>::iterator it=looklbc.find(string);
// 						if( it != looklbc.end()){
// 							scr=(nc[k]+1)*k_sc+ it->second;
// 						}
// 						else{
// 							double v=log(dbico(nxj,nc[k]));
// 							looklbc[string]=v;
// 							scr=(nc[k]+1)*k_sc+ v;
// 						}

// 						// if(looklbc[nxj][nc[k]]==-1) looklbc[nxj][nc[k]]=log(dbico(nxj,nc[k]));
// 						// scr=(nc[k]+1)*k_sc+looklbc[nxj][nc[k]];
						
// 						// scr=(nc[k]+1)*k_sc+log(dbico(nxj,(nc[k])));


// 						#if _MY_DEBUG_NEW_OPTFUN


// 							printf("  ");
// 							for(m=0;m<nbrV+1;m++) printf("H_kj[%d]=%lf ",m,H_kj[m]);fflush(stdout);
// 							printf("\n   [%d -%d][%d - %d] = %lf (%lf+%lf-%lf)\n",0,k,k+1,j,I_0k[k]+I_kj-sc2,I_0k[k],I_kj,sc2);
// 							printf("   [0 -?-%d][%d - %d] = %lf (%lf+%lf-%lf)\n",k,k+1,j,I[k]+I_kj-scr,I[k],I_kj,scr);
					
// 						#endif

// 						if( I_0k[k] - sc2 > I[k] - scr ){ //one cut in k ore more?
// 								t=I_0k[k] + I_kj-sc2;
// 								if (fmax<t){
// 									c=-k-1;// convention to refers to the two bins [0 k] [k+1 j]
// 									fmax=t;
// 									sctemp=sc2;
// 									nctemp=2;
// 								}
// 						}else{//more cuts
// 							t=I[k] + I_kj-scr;//[0.. cuts.. k-1][k j]
// 							if (fmax<t){
// 								c=k+1;
// 								fmax=t;
// 								sctemp=scr;
// 								nctemp=nc[k]+1;
// 							}
// 						}

// 						#if _MY_DEBUG_NEW
// 							printf("   f=%lf\n",fmax);fflush(stdout);
// 						#endif
// 			}

// 			I[j]=fmax+sctemp; 
// 			nc[j]=nctemp;
// 			memory_cut[j]=c;
// 			#if _MY_DEBUG_NEW_OPTFUN
// 				printf("\n>>>j=%d: fmax=%lf cut[%d]=%d ncut=%d\n",j,fmax,j,memory_cut[j],nc[j]);fflush(stdout);
// 			#endif
// 		}

// 	// free memory

// 	for(m=0;m<nbrV+1;m++) free(H_0k[m]);
// 	free(H_0k);

// 	for(m=0;m<nbrV;m++) free(nxyu[m]);
// 	free(nxyu);

// 	for(m=0;m<nbrV;m++) free(nxyu_k[m]);
// 	free(nxyu_k);	

// 	// 

// 	return I[np-1]/n;

// }




///////////////////////////////////////////////

//ruyx -> 0:u,1:uy,2:ux,3:uyx

double* computeMIcond_knml(int **uiyxfactors, int *ruiyx, int *r,int n, double* c2terms, double* looklog){

	double *I=(double*)calloc(2,sizeof(double));

	int j,u,ux,uy,uyx;

	double Hux=0,Huy=0,Huyx=0,Hu=0,SC=0;

	int *nui=(int*)calloc(ruiyx[0],sizeof(int));
	int *nuiy=(int*)calloc(ruiyx[1],sizeof(int));
	int *nuix=(int*)calloc(ruiyx[2],sizeof(int));
	int *nuiyx=(int*)calloc(ruiyx[3],sizeof(int));


	for(j=0;j<n;j++){
		nui[uiyxfactors[0][j]]++;
		nuiy[uiyxfactors[1][j]]++;
		nuix[uiyxfactors[2][j]]++;
		nuiyx[uiyxfactors[3][j]]++;
	}

	for(u=0;u<ruiyx[0];u++){
		if(nui[u]>0) Hu-=nui[u]*looklog[nui[u]];
		SC-=0.5*computeLogC(nui[u],r[0],c2terms);
		SC-=0.5*computeLogC(nui[u],r[1],c2terms);
	}

	for(uy=0;uy<ruiyx[1];uy++){
		if(nuiy[uy]>0) Huy-=nuiy[uy]*looklog[nuiy[uy]];
		SC+=0.5*computeLogC(nuiy[uy],r[0],c2terms);
	}

	for(ux=0;ux<ruiyx[2];ux++){
		if(nuix[ux]>0) Hux-=nuix[ux]*looklog[nuix[ux]];
		SC+=0.5*computeLogC(nuix[ux],r[1],c2terms);
	}

	for(uyx=0;uyx<ruiyx[3];uyx++){
		if(nuiyx[uyx]>0) Huyx-=nuiyx[uyx]*looklog[nuiyx[uyx]];
	}

	I[0]=(Hux+Huy-Hu-Huyx)/n;

	I[1]=I[0]-SC/n;


	free(nui);
	free(nuiy);
	free(nuix);
	free(nuiyx);

	return I;

}


///////////////////////////////////////////////

//rux -> 0:x,1;u,2:ux

double* computeMI_knml(int* xfactors,int* ufactors,int* uxfactors,int* rux,int n, double* c2terms,double* looklog, int flag){

	double *I=(double*)calloc(2,sizeof(double));

	int j,x,u,ux;

	double Hux=0,Hu=0,Hx=0,SC=0;

	int *nx=(int*)calloc(rux[0],sizeof(int));
	int *nu=(int*)calloc(rux[1],sizeof(int));
	int *nux=(int*)calloc(rux[2],sizeof(int));

	for(j=0;j<n;j++){
		nx[xfactors[j]]++;
		nu[ufactors[j]]++;
		nux[uxfactors[j]]++;
	}

	for(x=0;x<rux[0];x++){
		if(nx[x]>0) Hx-=nx[x]*looklog[nx[x]];
		if(flag==0 || flag==1) SC+=computeLogC(nx[x],rux[1],c2terms);
	}
	for(u=0;u<rux[1];u++){
		if(nu[u]>0) Hu-=nu[u]*looklog[nu[u]];
		if(flag==0 || flag==2) SC+=computeLogC(nu[u],rux[0],c2terms);
	}

	for(ux=0;ux<rux[2];ux++){
		if(nux[ux]>0) Hux-=nux[ux]*looklog[nux[ux]];
	}
	
	if(flag==0 || flag==1) SC-=computeLogC(n,rux[0],c2terms);
	if(flag==0 || flag==2) SC-=computeLogC(n,rux[1],c2terms);

	I[0]=looklog[n]+(Hu+Hx-Hux)/n;

	if(flag==0) I[1]=I[0]-0.5*SC/n;
	else I[1]=I[0]-SC/n;

	free(nx);
	free(nu);
	free(nux);

	return I;

}


////////////////////////////////////////////////////////////////////////////////////////////


double* computeMIcond_kmdl(int **uiyxfactors, int *ruiyx, int *r,int n, double* looklog)
{

	double *I=(double*)calloc(2,sizeof(double));

	int j,l,u,ux,uy,uyx;

	double Hux=0,Huy=0,Huyx=0,Hu=0,SC=0;

	int *nui=(int*)calloc(ruiyx[0],sizeof(int));
	int *nuiy=(int*)calloc(ruiyx[1],sizeof(int));
	int *nuix=(int*)calloc(ruiyx[2],sizeof(int));
	int *nuiyx=(int*)calloc(ruiyx[3],sizeof(int));


	for(j=0;j<n;j++){
		nui[uiyxfactors[0][j]]++;
		nuiy[uiyxfactors[1][j]]++;
		nuix[uiyxfactors[2][j]]++;
		nuiyx[uiyxfactors[3][j]]++;
	}

	for(u=0;u<ruiyx[0];u++){
		if(nui[u]>0) Hu-=nui[u]*looklog[nui[u]];
	}
	for(uy=0;uy<ruiyx[1];uy++){
		if(nuiy[uy]>0) Huy-=nuiy[uy]*looklog[nuiy[uy]];
	}
	for(ux=0;ux<ruiyx[2];ux++){
		if(nuix[ux]>0) Hux-=nuix[ux]*looklog[nuix[ux]];
	}
	for(uyx=0;uyx<ruiyx[3];uyx++){
		if(nuiyx[uyx]>0) Huyx-=nuiyx[uyx]*looklog[nuiyx[uyx]];
	}

	SC=0.5*(r[0]-1)*(r[1]-1)*looklog[n];
	// SC*=ruiyx[0];

	I[0]=(Hux+Huy-Hu-Huyx)/n;

	I[1]=I[0]-SC/n;

	free(nui);
	free(nuiy);
	free(nuix);
	free(nuiyx);

	return I;

}

///////////////////////////////////////////////

//rux -> 0:x,1;u,2:ux

double* computeMI_kmdl(int* xfactors,int* ufactors,int* uxfactors,int* rux,int n,double* looklog, int flag){

	double *I=(double*)calloc(2,sizeof(double));

	int j,x,u,ux;

	double Hux=0,Hu=0,Hx=0,SC=0;

	int *nx=(int*)calloc(rux[0],sizeof(int));
	int *nu=(int*)calloc(rux[1],sizeof(int));
	int *nux=(int*)calloc(rux[2],sizeof(int));

	for(j=0;j<n;j++){
		nx[xfactors[j]]++;
		nu[ufactors[j]]++;
		nux[uxfactors[j]]++;
	}

	for(x=0;x<rux[0];x++){
		if(nx[x]>0) Hx-=nx[x]*looklog[nx[x]];
	}
	for(u=0;u<rux[1];u++){
		if(nu[u]>0) Hu-=nu[u]*looklog[nu[u]];
	}

	for(ux=0;ux<rux[2];ux++){
		if(nux[ux]>0) Hux-=nux[ux]*looklog[nux[ux]];
	}
	
	SC=0.5*looklog[n];
	if(flag==0 || flag==1) SC*=(rux[0]-1);
	if(flag==0 || flag==1) SC*=(rux[1]-1);

	I[0]=looklog[n]+(Hu+Hx-Hux)/n;

	I[1]=I[0]-SC/n;

	free(nx);
	free(nu);
	free(nux);

	return I;

}