#include <stdio.h>
#include <stdlib.h>	/*malloc*/
#include <string.h>	/* memset*/
#include <math.h>
#include <stdint.h>
#include <iostream>

#include "computeInfo.h"

//# define M_PI		3.14159265358979323846	/* pi */
# define N_COL_NML 1000

unsigned long binomialCoeff(int n, int k)
{
    unsigned long res = 1;
 
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
       		C2 = sqrt( NN * M_PI_2)*exp( sqrt( 8/( 9*NN*M_PI ) ) + ( 3*M_PI-16 )/( 36*NN*M_PI ) );
       		c2terms[NN] = C2;
	   	} else {
	   		C2 = c2terms[NN];
	   	}
	}

	D=C2;
	logC=log(D);

	if( r > 2 )
	{
    	for( rr = 3; rr <= r; rr++)
		{
		  D=1+(NN/(1.0*(rr-2)*D));
		  logC=logC+log(D);
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
			C2 = sqrt( NN * M_PI_2)*exp( sqrt( 8/( 9*NN*M_PI ) ) + ( 3*M_PI-16 )/( 36*NN*M_PI ) );
		}
	  c2terms[NN] = C2;

	}

	D=C2;
	logC=log(D);
	
	if( r > 2 ){
		for( rr = 3; rr <= r; rr++){
			D=1+(NN/(1.0*(rr-2)*D));
			logC=logC+log(D);
		}
	}
	  
	return logC;
}

double lnfactorial(int n, double* looklog) {
	int y;
	double z;
	if (n == 0){
		z=1;
	}
	else {
		z = 0;
		for (y=2; y<=n; y++) z += looklog[y];
	}
	return z;
}

double logchoose(int n, int k, double* looklog){
	if(n==k || k==0){
		return 0;
	}
	double res = lnfactorial(n, looklog) - lnfactorial(k, looklog) - lnfactorial(n-k, looklog);
	return(res);
}

double logchoose(int n, int k, double* looklog, double** lookchoose){

	if(k<N_COL_NML){
		double val = lookchoose[k][n];
		if(val >= 0){
			return(val);
		}
	}

	if(n==k || k==0){
		return 0;
	}
	double res = lnfactorial(n, looklog) - lnfactorial(k, looklog) - lnfactorial(n-k, looklog);

	if(k<N_COL_NML){
		lookchoose[k][n]=res;
	}
	return(res);
}

double computeLogC(int N, int r, double* looklog, double* c2terms)
{
	double C2, logC, D;
	int rr,h;

	if(N<=1000)
	{
		if(c2terms[N] == -1){
			C2=0;
			for( h = 0; h <= N; h++ )
			{
				//C2 = C2 + binomialCoeff( N, h )*pow( ( h/( (double)N ) ), h )*pow( ( ( N-h )/( (double)N ) ), ( N-h ) );
				C2 = C2 + exp( logchoose(N, h, looklog) + log(pow( ( h/( (double)N ) ), h )) + log(pow( ( ( N-h )/( (double)N ) ), ( N-h ) )) );
	   		}
	   		c2terms[N] = C2;
	   	} else {
	   		C2 = c2terms[N];
	   	}
	} else {
    	if(c2terms[N] == -1){
       		C2 = sqrt( N * M_PI_2)*exp( sqrt( 8/( 9*N*M_PI ) ) + ( 3*M_PI-16 )/( 36*N*M_PI ) );
       		c2terms[N] = C2;
	   	} else {
	   		C2 = c2terms[N];
	   	}
	}

	D=C2;
	logC=log(D);

	if( r > 2 )
	{
    	for( rr = 3; rr <= r; rr++)
		{
		  D=1+(N/(1.0*(rr-2)*D));
		  logC=logC+log(D);
		}
  	}

	return logC;
}

double computeLogC( int N, int r, double* looklog, double** cterms)
{
	if(r<N_COL_NML){
		double val = cterms[r][N];
		if(val >= 0){
			return(val);
		}
	}

	double* sc_look = cterms[2];
	double C2, logC, D;
	int rr,h;

	if(N<=1000)
	{
		if(sc_look[N] == -1){
			C2=0;
			for( h = 0; h <= N; h++ )
			{
				//C2 = C2 + binomialCoeff( N, h )*pow( ( h/( (double)N ) ), h )*pow( ( ( N-h )/( (double)N ) ), ( N-h ) );
				C2 = C2 + exp( logchoose(N, h, looklog) + log(pow( ( h/( (double)N ) ), h )) + log(pow( ( ( N-h )/( (double)N ) ), ( N-h ) )) );
	   		}
            C2 = log(C2);
	   	} else {
	   		C2 = sc_look[N];
	   	}
	} else {
    	if(sc_look[N] == -1){
       		C2 = sqrt( N * M_PI_2)*exp( sqrt( 8/( 9*N*M_PI ) ) + ( 3*M_PI-16 )/( 36*N*M_PI ) );
            C2 = log(C2);
	   	} else {
	   		C2 = sc_look[N];
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

	if(r<N_COL_NML){
		cterms[r][N]=logC;
	}
	return logC;
}

double computeLogRC(int N, int reff, double* looklog, double** cterms){

	double RC = 0;
	for(int k=0; k<reff; k++){
		RC += pow(-1, k) * binomialCoeff(reff, k) * exp(computeLogC(N,  reff-k, looklog,  cterms));
	}
	return(log(RC));
}

double computeLogRCr(int N, int reff, double* looklog, double** cterms){

	double RC = 0;
	for(int k=0; k<reff; k++){
		RC += pow(-1, k) * exp(computeLogC(N, reff-k, looklog, cterms) - computeLogC(N, reff, looklog, cterms) + logchoose(reff,  k, looklog));
	}
	return(log(RC));
}

double computeLogHDC(int N, int r, int reff, double* looklog, double** cterms){

	//return(logchoose(r, reff, looklog) + computeLogRC(N, reff, looklog, cterms));
	if(N==0) return(0);
	if(N==1) return(log(r));
	//return(logchoose(r, reff, looklog) + computeLogRCr(N, reff, looklog, cterms) + computeLogC(N, reff, looklog, cterms));
	return(logchoose(r, reff, looklog) + computeLogC(N, reff, looklog, cterms));
}