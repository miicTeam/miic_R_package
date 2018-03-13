#include "probaOrientation.h"
#include <stdlib.h>

double* getOrientTplLVDegPropag( int nbrTpl, int* ptrAllTpl, double* ptrAllI3, int LV , int isDeg , int isPropag, int halfVStructures)
{
   	int i, j, k, errCode = 0;	// for loops


	int nbrRetProbaValues = -1;                             // Nbr proba to return
	double *ptrRetProbValues;						        // To return ProbArrowhead
                                                            // >> ProbArrowhead[1][i] 1<-->2 ProbArrowhead[2][i]
	                                                        // >> ProbArrowhead[3][i] 2<-->3 ProbArrowhead[4][i]
	                                                        //
	                                                        // ProbArrowhead > 0.5 corresponds to a arrow head  (>)
	                                                        // ProbArrowhead < 0.5 corresponds to a arrow tail  (-)
	nbrRetProbaValues = (4*nbrTpl);

	ptrRetProbValues = new double[nbrRetProbaValues];//malloc(nbrRetProbaValues * sizeof *ptrRetProbValues);
	
    // Initialise the arrowhead probabilities to 0.5
    for( i = 0; i < nbrTpl; i++ )
	{ for( j = 0; j < 4; j++ ) { ptrRetProbValues[i+j*nbrTpl] = 0.5; } }


    //////////////////////////////////////////////////////////////
    //  iteratively converge towards partially oriented graphs 
    //  including possible latent variables and Propagation/Non-Propagation rules.
	OrientTpl_LV_Deg_Propag( nbrTpl, ptrAllTpl, ptrAllI3, ptrRetProbValues, LV, isDeg, isPropag, halfVStructures);
    
	return ptrRetProbValues;
}

