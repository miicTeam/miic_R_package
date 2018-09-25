#include <stdio.h>
#include <stdlib.h>	
#include <string.h>	
#include <math.h>
#include <stdint.h>
#include <cmath>
#include <iostream>
#include "computeInfo.h"

# define M_PI		3.14159265358979323846	/* pi */
# define N_COL_NML 1000
using namespace std;

unsigned long binomialCoeff(int n, int k)
{
    int res = 1;
 
    // Since C(n, k) = C(n, n-k)
    if ( k > n - k )
        k = n - k;
 
    // Calculate value of [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]
    for (int i = 0; i < k; ++i)
    {
        res *= (n - i);
        res /= (i + 1);
    }
    return res;
}

int getDVect( int* myPtrAllLevels, int* myPtrVarIdx, int myNbrAllVar, int* ioPtrDVectVar )
{
	int errCode = 0;
	int i;
	// Compute the di for the variables
	// dmax = 1, dmax-1=rmax, ../.., d2 = rmax*..*r3, d1 = rmax*..*r2, d0 = rmax*..*r2*r1
   	ioPtrDVectVar[(myNbrAllVar-1)] = 1;
	for( i = (myNbrAllVar-2); i >= 0; i-- )
	{
		ioPtrDVectVar[i] = ioPtrDVectVar[i+1]*myPtrAllLevels[myPtrVarIdx[i+1]];
	}

	return errCode;	
}

// myCplx = 0 --> MDL, myCplx = 1 --> NML
int computeInfoAndCplx( int *myNxyui, int *myDVect, int *myLevels, int *myVarIdx, int mySSize, int mySSizeEff, int myNbrVar, double* infoAndCplx, int myCplxType, double* c2terms )
{
	int errCode = 0;

	int i, j;
	int Lui;		            // Partial index to read the count (ie, only with ui)
	int Lxui;		            // Partial index to read the count (ie, only with x and ui)
	int Lyui;		            // Partial index to read the count (ie, only with y and ui)
	int LLxyui;		            // Full index to read the count (ie, X, Y and {ui})

	int Nui;		            // Sum over X and Y for a {ui}
	int Nxui;		            // Sum over Y of Nxyui
	int Nyui;		            // Sum over X of Nxyui
	int NN = 0;		            // Sum over {ui} of Nui
	
	int tmpRescaleN;            // temp variable for rescaling the 'n' considering the effective sample size

	// Compute the first part
	double NI = 0.0;			// NI = NIxy - NIx - NIy	
	double NML = 0.0, logCNML = 0.0, MDL = 0.0;

	int stopUiStatesLoop = 0;


	//  Define a vector to set the combinaison of the {ui}
	// (var := X, Y, U1, U2, etc...)
	int *uiVal = NULL;
	if( myNbrVar > 2 )
	{
		uiVal = (int*)malloc((myNbrVar-2)*sizeof(int));
		memset( uiVal, 0, (myNbrVar-2)*sizeof(int) );


	} 
	// Loop over {ui} combinations
	do
	{
		// --------------- Partial index with ui ---------------
		// Compute a partial index that takes into account only the combination of {ui}
		Lui=0;
		if( uiVal != NULL )
		{
    		for( i = (myNbrVar-1); i >=2; i-- ) 
			{ Lui = Lui + myDVect[i]*uiVal[i-2]; }
		}
		// ---------------  Nui ---------------
		// Compute Nui (ie, sum over X and Y cases of Nxyui)
		Nui=0;
		for( i = 0; i < myLevels[myVarIdx[0]]; i++ )
		{

			for( j = 0; j < myLevels[myVarIdx[1]]; j++ )
			{ 

				LLxyui = Lui + i*myDVect[0] + j*myDVect[1];
				Nui = Nui + myNxyui[LLxyui];
      		}
    	}

		// --------------- NN ---------------
		// Add the current Nui to the sum over X, Y and {ui} of Nxyui
		NN = NN + Nui;


		// --------------- NI ---------------
		// Update the logCNML coefficient if needed (1/C(Nui, ry))
		//if( myCplx == 1 ) { logCNML = - computeLogC( Nui, myLevels[myVarIdx[1]] ); }
		if( myCplxType == 1 )
		{
		    tmpRescaleN = (int)floor( 0.5 + ((double)Nui*mySSizeEff)/mySSize ); 
		    logCNML = - computeLogC( tmpRescaleN, myLevels[myVarIdx[1]], c2terms);
		}
		// Update the NI (sum over ui and over x of Nxui)
		for( i = 0; i < myLevels[myVarIdx[0]]; i++ )
		{
			Nxui = 0;
			for( j = 0; j < myLevels[myVarIdx[1]]; j++ )
			{
				Lxui = Lui + i*myDVect[0] + j*myDVect[1];
				Nxui = Nxui + myNxyui[Lxui]; 
				//printf("# --!!--> Nxui = %d\n", Nxui);
			}
			// Add to the NI
			if( Nxui > 0 )
			{ 
				NI = NI - Nxui*std::log(static_cast<double>(Nxui/(1.0*Nui)));//log(Nxui/(1.0*Nui)); 
				//if( myCplx == 1 ) { logCNML = logCNML + computeLogC( Nxui, myLevels[myVarIdx[1]] ); }
				if( myCplxType == 1 )
				{ 
				    tmpRescaleN = (int)floor( 0.5 + ((double)Nxui*mySSizeEff)/mySSize ); 
				    logCNML = logCNML + computeLogC( tmpRescaleN, myLevels[myVarIdx[1]], c2terms);
				}
			}

		}
		if( Nui > 0 && myCplxType == 1 ){ NML = NML + 0.5*logCNML; }

        // Update the logCNML coefficient if needed (1/C(Nui, rx))
		//if( myCplx == 1 ) { logCNML = - computeLogC( Nui, myLevels[myVarIdx[0]] ); }
		if( myCplxType == 1 )
		{ 
		    tmpRescaleN = (int)floor( 0.5 + ((double)Nui*mySSizeEff)/mySSize );
		    logCNML = - computeLogC( tmpRescaleN, myLevels[myVarIdx[0]], c2terms);
		}
		// Increment the NIy (sum over ui and over y of Nxyui)
		for( j = 0; j < myLevels[myVarIdx[1]]; j++ )
		{
			Nyui = 0;
			for( i = 0; i < myLevels[myVarIdx[0]]; i++ )
			{
				Lyui = Lui + i*myDVect[0] + j*myDVect[1];
				Nyui = Nyui + myNxyui[Lyui];
				//printf("# --!!--> Nyui = %d\n", Nyui);
			}
			// Add to the NI
			if( Nyui > 0 )
			{ 
				NI = NI - Nyui*std::log(static_cast<double>(Nyui/((double)Nui)));//log(Nyui/((double)Nui));
				//if( myCplx == 1 ) { logCNML = logCNML + computeLogC( Nyui, myLevels[myVarIdx[0]] ); }
				if( myCplxType == 1 )
				{
				    tmpRescaleN = (int)floor( 0.5 + ((double)Nyui*mySSizeEff)/mySSize ); 
				    logCNML = logCNML + computeLogC( tmpRescaleN, myLevels[myVarIdx[0]], c2terms);
				}
			}

		}
		if( Nui > 0 && myCplxType == 1){ NML = NML + 0.5*logCNML; }

		// Update the NI (sum over ui, over x and over y of Nxyui)
		for( i = 0; i < myLevels[myVarIdx[0]]; i++ )
		{
			for( j = 0; j < myLevels[myVarIdx[1]]; j++ )
			{
				LLxyui = Lui + i*myDVect[0] + j*myDVect[1];
				if( myNxyui[LLxyui] > 0 ){ NI = NI + myNxyui[LLxyui]*std::log(static_cast<double>(myNxyui[LLxyui]/((double)Nui)));}//log(myNxyui[LLxyui]/((double)Nui)); }
			}
		}	

		stopUiStatesLoop = 1;
		if( uiVal != NULL )
		{
			// Compute the next ui combinaison
			for( i = 2; i < myNbrVar; i++ )
			{
				if( uiVal[i-2] < (myLevels[myVarIdx[i]]-1) )
				{
					uiVal[i-2]++;
					stopUiStatesLoop = 0;
					i = myNbrVar;

				} else{
					uiVal[i-2] = 0;				
				}
			}
		} 		
	} while( stopUiStatesLoop == 0 );

    // Rescale the NI Value and the NN
	NI = ((double)NI*mySSizeEff)/mySSize;
	NN = (int)floor( 0.5 + (((double)NN*mySSizeEff)/mySSize));

	if( myCplxType == 0 ) 
	{ 
		// Compute the MDL
		MDL = 0.5*( myLevels[myVarIdx[0]] - 1 )*( myLevels[myVarIdx[1]] - 1 )*std::log(static_cast<double>(NN));//log( NN );
		for( i = 2; i < myNbrVar; i++ )
		{
			MDL = MDL*myLevels[myVarIdx[i]];
		}
	}

	// Set the io array
	infoAndCplx[0] = NI;
	infoAndCplx[1] = -1;
	if( myCplxType == 0) { infoAndCplx[1] = MDL; } else if( myCplxType == 1 ) { infoAndCplx[1] = NML; }

	// --- FREE FREE FREE ---
	free(uiVal);
	uiVal = NULL;

	return errCode;
}



double compute_rescaledLogC( int sampleSize, int sampleSizeEff, int N, int r , double* c2terms )
{
        //double C[r+1];
	double C2, logC, D;
	int rr,h,NN;

	NN = (int)floor(0.5 + ((double)N*sampleSizeEff)/sampleSize);

	if(NN<=10)
	{
		C2=0;
		if(c2terms[NN] == -1){
			for( h = 0; h <= NN; h++ )
			{
			  C2 = C2 + binomialCoeff( NN, h )*pow( (h/((double)NN)), h )*pow( ((NN-h)/((double)NN ) ), ( NN-h ) );
	   		}
	   		c2terms[NN] = C2;
   		} else {
	   		C2 = c2terms[NN];
	   	}

	} else {
    	if(c2terms[NN] == -1){
       		C2 = pow( (NN * M_PI_2),0.5)*exp( pow(( 8/( 9*NN*M_PI ) ),0.5) + ( 3*M_PI-16 )/( 36*NN*M_PI ) );
       		c2terms[NN] = C2;
	   	} else {
	   		C2 = c2terms[NN];
	   	}
	}

	D=C2;
	logC=std::log(static_cast<double>(D));//log(D);

	if( r > 2 )
	{
    	for( rr = 3; rr <= r; rr++)
		{
		  D=1+(NN/(1.0*(rr-2)*D));
		  logC=logC+std::log(static_cast<double>(D));//log(D);
		}
  	}

	return logC;
}

double compute_LogC_C2( int N, int r, double* c2terms)
{
	double C2, logC, D;
	int rr,h,NN;

	NN=N;

        if( c2terms[NN] > 0){

 	  C2 = c2terms[NN];

	} else {

		if(NN<=10)
		{
			C2=0;
			for( h = 0; h <= NN; h++ )
				C2 = C2 + binomialCoeff( NN, h )*pow( (h/((double)NN)), h )*pow( ((NN-h)/((double)NN ) ), ( NN-h ) );

		} else {
			C2 = pow(( NN * M_PI_2),0.5)*exp( pow(( 8/( 9*NN*M_PI ) ),0.5) + ( 3*M_PI-16 )/( 36*NN*M_PI ) );
		}
	  c2terms[NN] = C2;

	}

	D=C2;
	logC=std::log(static_cast<double>(D));//log(D);
	
	if( r > 2 ){
		for( rr = 3; rr <= r; rr++){
			D=1+(NN/(1.0*(rr-2)*D));
			logC=logC+std::log(static_cast<double>(D));//log(D);
		}
	}
	  
	return logC;
}

double computeLogC( int N, int r, double* c2terms)
{
	double C2, logC, D;
	int rr,h;

	if(N<=10)
	{
		if(c2terms[N] == -1){
			C2=0;
			for( h = 0; h <= N; h++ )
			{
				C2 = C2 + binomialCoeff( N, h )*pow( ( h/( (double)N ) ), h )*pow( ( ( N-h )/( (double)N ) ), ( N-h ) );
	   		}
	   		c2terms[N] = C2;
	   	} else {
	   		C2 = c2terms[N];
	   	}
	} else {
    	if(c2terms[N] == -1){
       		C2 = pow(( N * M_PI_2), 0.5)*exp( pow(( 8/( 9*N*M_PI ) ), 0.5) + ( 3*M_PI-16 )/( 36*N*M_PI ) );
       		c2terms[N] = C2;
	   	} else {
	   		C2 = c2terms[N];
	   	}
	}

	D=C2;
	logC=std::log(static_cast<double>(D));//log(D);

	if( r > 2 )
	{
    	for( rr = 3; rr <= r; rr++)
		{
		  D=1+(N/(1.0*(rr-2)*D));
		  logC=logC+std::log(static_cast<double>(D));//log(D);
		}
  	}

	return logC;
}


double computeLogC( int N, int r, double** cterms)
{
	if(r<1000){
		if(cterms[r][N] != -1){
		return(cterms[r][N]);
		}
	}

    double* c2terms = cterms[2];
	double C2, logC, D;
	int rr,h;

	if(N<=10)
	{
		if(c2terms[N] == -1){
			C2=0;
			for( h = 0; h <= N; h++ )
			{
				C2 = C2 + binomialCoeff( N, h )*pow( ( h/( (double)N ) ), h )*pow( ( ( N-h )/( (double)N ) ), ( N-h ) );
	   		}
            C2 = log(C2);
	   	} else {
	   		C2 = c2terms[N];
	   	}
	} else {
    	if(c2terms[N] == -1){
       		C2 = sqrt( N * M_PI_2)*exp( sqrt( 8/( 9*N*M_PI ) ) + ( 3*M_PI-16 )/( 36*N*M_PI ) );
            C2 = log(C2);
	   	} else {
	   		C2 = c2terms[N];
	   	}
	}

	logC=C2;
	D=exp(C2);

	if( r > 2 )
	{
    	for( rr = 3; rr <= r; rr++)
		{
		  D=1+(N/(1.0*(rr-2)*D));
		  logC=logC+log(D);
		}
  	}

	if(r<1000){
		cterms[r][N]=logC;
	}
	return logC;
}

double compute_belief_factor(int k, int nx, int n_samples, double** belief_factor){

	double factor;
	if(k >= N_COL_NML){
		if(k==0 && nx!=0) factor = 1.0;
		else if(nx==0) factor = 0.0;
		else factor = (0.5 + 0.8 * pow((0.5*(1.0/k) * (1.0*n_samples/nx)), 2.0/3) );
		//else factor = (0.5 + 0.8 * (0.5*(1.0/k) * (1.0*n_samples/nx)));
	}
	else{
		factor = belief_factor[k][nx];
		if( (factor == 0.0) && (nx != 0) ){
			if(k==0) factor = 1.0;
			else factor = (0.5 + 0.8 * pow((0.5*(1.0/k) * (1.0*n_samples/nx)), 2.0/3) );
			//else factor = (0.5 + 0.8 * (0.5*(1.0/k) * (1.0*n_samples/nx)));
			belief_factor[k][nx] = factor;
		}
	}
	return factor;
}