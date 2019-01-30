#include "probaOrientation.h"
#include <stdlib.h>
#include <stdio.h>

#define _MY_DEBUG_ 0

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

#if _MY_DEBUG_
	// Print inputs
	printf("\n# ---- Printed from C [getAllInfo]: ----\n");
	printf("\n# Input args:\n# --------");
	printf("# --> nbrTpl: %d\n", nbrTpl);
	printf("# --> isDeg: %d\n", isDeg);
	printf("# --> isPropag: %d\n", isPropag);
	printf("# --> LV: %d\n", LV);
	printf("# --> All I3:\n");
	for( i = 0; i < nbrTpl; i++ ){ printf("# I3[%d] = %g\n", i, ptrAllI3[i]); }
	
	printf("# --> All Tpl:\n");
	for( i = 0; i < nbrTpl; i++ )
	{ 
	    printf("# Tpl[1] = %d -- Tpl[2] = %d -- Tpl[3] = %d\n", ptrAllTpl[i+0*nbrTpl], ptrAllTpl[i+1*nbrTpl], ptrAllTpl[i+2*nbrTpl]);
	}
	
	printf("# --> All Prob:\n");
    for( i = 0; i < nbrTpl; i++ )
	{
	    printf("# Proba[1] = %g <--> Proba[2] = %g -- Proba[3] = %g <--> Proba[4] = %g\n", ptrRetProbValues[i+0*nbrTpl], ptrRetProbValues[i+1*nbrTpl], ptrRetProbValues[i+2*nbrTpl], ptrRetProbValues[i+3*nbrTpl]);
	}
	printf("# --------\n");
#endif // _MY_DEBUG_

    //////////////////////////////////////////////////////////////
    //  iteratively converge towards partially oriented graphs 
    //  including possible latent variables and Propagation/Non-Propagation rules.
	OrientTpl_LV_Deg_Propag( nbrTpl, ptrAllTpl, ptrAllI3, ptrRetProbValues, LV, isDeg, isPropag, halfVStructures);
    //////////////////////////////////////////////////////////////
    
#if _MY_DEBUG_
    printf("\n# Orientation  (P>0.5: arrow-head; P<0.5: arrow-tail) :\n");
    for( i = 0; i < nbrTpl; i++ )
    { 
        printf("! Tpl %i (%g-%g) %i (%g-%g) %i -- I3=%g\n", ptrAllTpl[0*nbrTpl+i], ptrRetProbValues[0*nbrTpl+i], ptrRetProbValues[1*nbrTpl+i], ptrAllTpl[1*nbrTpl+i], ptrRetProbValues[2*nbrTpl+i], ptrRetProbValues[3*nbrTpl+i], ptrAllTpl[2*nbrTpl+i],ptrAllI3[i] );
    }
#endif // _MY_DEBUG_

	return ptrRetProbValues;
}

