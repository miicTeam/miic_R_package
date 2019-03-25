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

#include "utilities.h" // isNumberOfLevelsLessTwo
#include "modules_MI.h"
#include "computeInfo.h"

#define STEPMAX 50
#define STEPMAXOUTER 50
#define EPS 1e-5
#define FLAG_CPLX 1

using namespace std;

// module to compute the conditional mutual information for mixed continuous and discrete variables 

//////////////////////////////////////////////////////////////////////////////////////////////


//#define _MY_DEBUG_MInoU 1
//#define _MY_DEBUG_MI 1
//#define _MY_DEBUG_NEW_OPTFUN 1
//#define _MY_DEBUG_NEW 1

/////////////////////////////////////////////////////////////////////////////

// template <typename T>
// string ToString(T val)
// {
//     stringstream stream;
//     stream << val;
//     return stream.str();
// }

///////////////////////////////////////////////////////////////////////////////////////////

void reset_u_cutpoints(int** cut, int nbrUi, int* ptr_cnt, int* ptrVarIdx, int init_nbin, int maxbins,
                       int lbin, int* r, int* AllLevels, int n){

    for(int l=2;l<(nbrUi+2);l++){
        if(ptr_cnt[ptrVarIdx[l]]==1){
            for(int j=0;j<init_nbin-1;j++) {
                cut[l][j]=j*lbin+lbin-1;
            }
            cut[l][init_nbin-1]=n-1;
            for(int j=init_nbin;j<maxbins;j++) {
                cut[l][j]=0;
            }
            r[l]=init_nbin;
        }else{
            r[l]=AllLevels[ptrVarIdx[l]];
        }
    }
}

//GENERAL VARIABLES NOTATION

//int n: number of sample
//int coarse : coarse graining step : minimum length of a bin : considering only n/coarse possible cuts 

//double* looklog : lookup table for the logarithms of natural numbers up to n
//double* c2terms : precomputed vectors for the MDL complexity

///////////////////////////////////////////////////////////////////////////////////////////


void optfun_onerun_kmdl_coarse(int *sortidx_var, int *data, int nbrV, int **factors, int *r, double sc,
                                      int sc_levels1, int previous_levels, int n, int nnr, int *cut, int *r_opt, 
                                      int maxbins, double* looklog, double* lookH, double** cterms, int cplx, 
                                      double* sample_weights, bool flag_efN)
{
/*  DYNAMIC PROGRAMMING OPTIMIZATION FUNCTION
 * 
 *  Takes as input the factors and joint factors and optimizes a variable by maximizing a given information.
 *  nbrV determines the number of factors to consider.
 * 
 * ////////////////////////////////////////////
 * INPUT
 * 
 * The optimizing variable is defined by <sortidx_var> and <data>.
 * <sortidx> Contains the rank indices (ordering) of the optimized variable, <data> contains the ranks.
 * 
 * <nbrV> is the number of terms in the optimization functions, consisting in variables and joint variables,
 * which (discretized) values are stored in <factors> : ( <nbrV> x n )  matrix :
 * Example <nbrV>=2 : Optimizing X on I(X;Y)
 *      xy factors : factors[0] 
 *      x factors : factors[1]
 * 
 * <r> is the number of unique values of each factor.
 * 
 * <sc> is the stochastic complexity used for the simple MDL cost.
 * Its value is 0.5 * number_of_other_levels where number of other levels is the product of number of unique
 * observed levels of all variables except the variable being optimized.
 * 
 * <sc_levels> and <previous_levels> are respectively the number of levels of the other variable in the mutual 
 * information term and the number of levels of the otpimized variable in the previous iteration (or initial number
 * for the first iteration).
 * 
 * <n> is the number of samples in total.
 * 
 * <nnr> is the number of non repeated samples.
 * 
 * <cut> is the vector containing cutpoints for discretizing the optimized variables. Modified at the end wih the
 * optimal solution with a call to reconstruction_cut_coarse().
 * 
 * <r_opt> is a pointer to the number of bins of the optimized variable, also updated by the call to 
 * reconstruction_cut_coarse().
 * 
 * <maxbins> is the maximum number of bins allowed. It controls the resolution of possible cutpoints : 
 * For example, 50 maximum cutpoints on a 1000 samples vector results in possible cutpoints every 1000/50 = 20 positions.
 * This speeds up significantly the dynamic programming at the cost of finding approximated solutions.
 * 
 * <looklog, lookH, cterms> are lookup tables for computing respectively log, entropy and stochastic complexity (as in 
 * Kontkanen & Myllymäki, 2005)
 * 
 * <cplx> is the choice of complexity : 0 for simple MDL (product of all observed values, see <sc>) and 1 for NML with
 * stochastic complexity and a combinatorial term with the previous number of levels.
 * 
 * <sample_weights> are unique weights for each sample to account for effective N != N, if the <flag_efN> is set to true.
 * 
 * 
 * Use-cases : 
 *  Optimizing Ik(X;Y) by choosing discretization on X:
 *      <sortidx_var>, <data> -> Contains the ranks of X
 *      <nbrV> -> 2
 *      <factors>
 *           [0] : The discretized Y
 *           [1] : A vector of one repeated level (starting from X with one single bin)
 *      <sc_levels1>        -> the number of levels of discretized Y
 *      <previous_levels>   -> the number of levels of discretizedd X in the previous iteration (for cost correction)
 *      <cut>               -> the vector that will contain the cutpoints for X
 * 
 */

    int i,j,k,m;

    int coarse=ceil(1.0*nnr/maxbins);//step coarse graining
    if (coarse<1) coarse=1;
    int np=ceil(1.0*nnr/coarse); //number of possible cuts 

    //temp variables to memorize optimization cuts
    int memory_cuts_idx[np];//indexes of the cuts (1..np)
    int memory_cuts_pos[np];//positions of the cuts (1..n)

    //dynamic programming optimize function and memorize of cuts
    double Imax;//I-kmdl

    double weights[4];

    //function max value at each step
    //double* Ik=(double *)calloc(np,sizeof(double));
    double Ik[np]; // The optimal information value found at each idx.
    double Ik_kj; // The information value for a bin fom idx k to j.
    double t;
        
    int xyu;

    int ir; // Iterator on non repeated values
    int nk,nj; // Number of values in intervals [0,k], [0,j], [k,j]
    int njforward,nkforward; // Indexes at current position of j and k


    //entropy in kj interval for the <nbrV>+1 terms
    double Hk_kj[nbrV];


    int **nxyu=(int**)calloc(nbrV,sizeof(int*));// x y u xu yu xyu
    for(m=0;m<nbrV;m++){
        nxyu[m]=(int*)calloc(r[m],sizeof(int));
    }

    int **nxyu_k=(int**)calloc(nbrV,sizeof(int*));// x y u xu yu xyu
    for(m=0;m<nbrV;m++){
        nxyu_k[m]=(int*)calloc(r[m],sizeof(int));
    }
    int nxyu_ef;
        
    //boolean vector, check_repet[i] is true if data[i]!=data[i+1]
    bool check_repet[n];
    for(i=0;i<(n-1);i++){
        check_repet[i] = (data[sortidx_var[i+1]]!=data[sortidx_var[i]]);
    }

    ////matrix with number of observations of each level (of the non-opt variable) at idx j
    //int ***rx_0j=(int ***)calloc(sc_levels1,sizeof(int**));
    //for(k=0;k<sc_levels1;k++){
    //    rx_0j[k]=(int **)calloc(n,sizeof(int*));
    //    for(i=0;i<n;i++){
    //        rx_0j[k][i]=(int *)calloc(r[1],sizeof(int));
    //    }
    //}

    ////matrix with number of unique levels of the non-opt variable from idx k to j on the opt variable.
    //int ***rx_kj_diff=(int ***)calloc(n,sizeof(int**));
    //for(k=0;k<n;k++){
    //    rx_kj_diff[k]=(int **)calloc(n,sizeof(int*));
    //    for(i=0;i<n;i++){
    //        rx_kj_diff[k][i]=(int *)calloc(r[1],sizeof(int));
    //    }
    //}

    double ef_nj = 0;
    double efN_factor;
    double ef_nk;

    //////////// WEIGHTS /////////////////////////////////////
    if(nbrV==2){
        weights[0]=-1;
        weights[1]=1;
    }

    //// Filling the rx_0j matrix
    //for(xyu=0; xyu<r[1]; xyu++){
    //    rx_0j[factors[2][sortidx_var[0]]][0][factors[1][sortidx_var[0]]] = 1;
    //}
    //for(j=1; j<n; j++){
    //    for(k=0; k<sc_levels1; k++){
    //        for(xyu=0; xyu<r[1]; xyu++){
    //            rx_0j[k][j][xyu] = rx_0j[k][j-1][xyu];
    //        }
    //    }
    //    rx_0j[factors[2][sortidx_var[j]]][j][factors[1][sortidx_var[j]]] ++;
    //}
    //// Filling the rx_kj difference in levels from k to j matrix
    //int nlevels[r[1]];
    //for(k=0;k<n;k++){
    //    for(j=k;j<n;j++){
    //        for(xyu=0; xyu<r[1]; xyu++){
    //            nlevels[xyu] = 0;
    //            for(int lvl=0;lvl<sc_levels1;lvl++){
    //                    nlevels[xyu] += (rx_0j[lvl][j][xyu] - rx_0j[lvl][k][xyu]) > 0;
    //                }
    //            rx_kj_diff[k][j][xyu] = nlevels[xyu];
    //        }
    //    }
    //}

    //computing statistics of the <nbrV> terms and the entropy

    njforward=0;//iterator on values

    //moving j over the np possible cuts 
    for(j=0;j<np;j++){ //j=1...n-1
        //COMPUTING STATISTICS AND FUNCTION FOR THE INTEVAL [0 j] -> I_0k

        //computing statistics of the <nbrV> terms and the entropy

        ir=0;//iterator on not repeated values
        // njforward iterator on values
        while((ir<coarse) && (njforward<n)){

            for(m=0;m<nbrV;m++){

                xyu=factors[m][sortidx_var[njforward]];
                nxyu[m][xyu]++;
            }

           if(flag_efN) ef_nj += sample_weights[njforward];
            if(njforward+1 < n){//check no repetition
                ir += int(check_repet[njforward]);
            }
            njforward++;
        }

        if(flag_efN) efN_factor = ef_nj / njforward;

        Ik[j]=0;
        for(m=0;m<nbrV;m++){
            Hk_kj[m]=0;
            for(xyu=0; xyu<r[m]; xyu++){

                //nxyu_ef = nxyu[m][xyu];
                nxyu_ef = flag_efN ? int(efN_factor * nxyu[m][xyu] + 0.5) : nxyu[m][xyu];
                Hk_kj[m] -= weights[m] * lookH[nxyu_ef];

                if(m == 1){ //herve
                    if(cplx==0 && nxyu[m][xyu]>0) Hk_kj[m] -= sc * looklog[n];
                    else if(cplx==1){
                        Hk_kj[m] -= computeLogC(nxyu_ef,  sc_levels1, looklog,  cterms);
                        //if(sc_levels1>1) cout << njforward << "   " << nxyu_ef << " " << sc_levels1 << " " << rx_kj_diff[0][njforward-1][xyu] << ", " << 
                        //                         computeLogC(nxyu_ef, sc_levels1, looklog, cterms) << " " << 
                        //                         computeLogHDC(nxyu_ef, sc_levels1, rx_kj_diff[0][njforward-1][xyu], looklog, cterms) << endl;
                        //Hk_kj[m] -= computeLogHDC(nxyu_ef, sc_levels1, rx_kj_diff[0][njforward-1][xyu], looklog, cterms);
                    }
                }
            }
            Ik[j] += Hk_kj[m];
        }

        //njforward: number of points between 0 and j 
        nj=njforward;


        Imax=-DBL_MAX;

        for(m=0;m<nbrV;m++) {
            for(xyu=0;xyu<r[m];xyu++) nxyu_k[m][xyu]=nxyu[m][xyu];
        }
        ef_nk = ef_nj;
        

        //moving k 

        //k iterator on possible cuts
        nkforward=0;//iterator on values
        //it iterator on not repeated vales

        //Before trying to create bins : solution is one single bin from 0 to j.
        memory_cuts_idx[j] = 0;
        memory_cuts_pos[j] = 0;

        for(k=0;k<j;k++){//k=1...n-2 possible cuts

            ir=0;//iterator on not repeated values
            while( ir<coarse ){

                for(m=0;m<nbrV;m++){

                    xyu=factors[m][sortidx_var[nkforward]];
                    nxyu_k[m][xyu]--;
                }

                if(flag_efN) ef_nk -= sample_weights[nkforward];
                ir += int(check_repet[nkforward]);
                nkforward++;
            }

            if(flag_efN) efN_factor = ef_nk / (njforward-nkforward);

            Ik_kj=0;
            for(m=0;m<nbrV;m++){
                Hk_kj[m]=0;
                for(xyu=0; xyu<r[m]; xyu++){

                    nxyu_ef = flag_efN ? int(efN_factor * nxyu_k[m][xyu] + 0.5) : nxyu_k[m][xyu];
                    Hk_kj[m] -= weights[m] * lookH[nxyu_ef];//j_efN * nxyu_k[m][xyu]*looklog[int(j_efN * nxyu_k[m][xyu] + 0.5)];

                    if(m == 1){ //herve
                        if(cplx==0 && nxyu_k[m][xyu]>0) Hk_kj[m] -= sc * looklog[n];
                        else if(cplx==1){
                            Hk_kj[m] -= computeLogC(nxyu_ef,  sc_levels1, looklog,  cterms);
                            //if(sc_levels1>1) cout << "\t" << nkforward << "   " << nxyu_ef << " " << sc_levels1 << " " << rx_kj_diff[nkforward-1][njforward-1] << ", " << 
                            //                        computeLogC(nxyu_ef, sc_levels1, looklog, cterms) << " " << 
                            //                        computeLogHDC(nxyu_ef, sc_levels1, rx_kj_diff[0][njforward-1][xyu], looklog, cterms) << endl;
                            //Hk_kj[m] -= computeLogHDC(nxyu_ef, sc_levels1, rx_kj_diff[nkforward-1][njforward-1][xyu], looklog, cterms);
                        }
                    }
                }
                Ik_kj += Hk_kj[m];
            }

            if(cplx==1){
                Ik_kj -= log((1.0*np-1)/(previous_levels-1) - 1) + 1; // Combinatorial approximation
            }

            nk=nkforward-1;//position of the actual possible cut

            if((Ik[k] + Ik_kj) > Ik[j]) {
                t=Ik[k] + Ik_kj; //[0.. cuts.. k-1][k j] //herve
                if (Imax<t){
                    Imax=t;
                    Ik[j]=Imax; // optimized function for the interval [0 j]
                    if(memory_cuts_idx[k+1] == 0){
                        memory_cuts_idx[j] = -k - 1; // index  of the (last) optimal cut
                    }
                    else{
                        memory_cuts_idx[j] = k + 1; // index  of the (last) optimal cut
                    }
                    memory_cuts_pos[j] = nk; // position  of the (last) optimal cut
                }
            }
        }


    }

    for(m=0;m<nbrV;m++) free(nxyu_k[m]);
    free(nxyu_k);

    for(m=0;m<nbrV;m++) free(nxyu[m]);
    free(nxyu);

    //for(k=0;k<sc_levels1;k++){
    //    for(i=0;i<n;i++){
    //        free(rx_0j[k][i]);
    //    }
    //    free(rx_0j[k]);
    //};
    //free(rx_0j);

    //for(k=0;k<n;k++){
    //    for(i=0;i<n;i++){
    //        free(rx_kj_diff[k][i]);
    //    }
    //    free(rx_kj_diff[k]);
    //}
    //free(rx_kj_diff);

    // reconstruction of the optimal cuts from the memory cuts indexes and positions
    *r_opt=reconstruction_cut_coarse(memory_cuts_idx, memory_cuts_pos, np, n, cut);
}



////////////////////////////////////////////////////////////////////////////////








//////////////////////////////////////////////////////////////////////////////////////////////
//compute I(x,y)

//dynamic programming for optimizing variables binning 

//optimize on x I(x,y): Hx - Hxy - kmdl
//optimize on y I(x,y): Hy - Hxy - kmdl 
//until convergence

double* compute_Ixy_alg1(int** data, int** sortidx, int* ptr_cnt, int* ptrVarIdx,  int* AllLevels, 
                         int n, int maxbins, int initbins, int **cut, int *r, double* c2terms, double** cterms, 
                         double* looklog, double* lookH, int cplx, double* sample_weights, bool flag_effN)
{

    int j,l,ll;

    /////////////////////////////////////

    double* res_temp=(double *)calloc(2,sizeof(double));//res_tempults res_temp[0]->I,res_temp[1]->I-k

    //////////////////////////////////
    //allocation factors  x y

    int **datafactors;
    datafactors=(int **)calloc((2),sizeof(int*));

    for(l=0;l<(2);l++){
        datafactors[l]=(int *)calloc(n,sizeof(int));
    }

    int* r_temp=(int *)calloc(3,sizeof(int));

    int *ptr=(int *)calloc(2,sizeof(int));

    int *xy_factors = (int *)calloc(n,sizeof(int));

    int rxy=0;
    
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

    jointfactors_u(datafactors,ptr,n, 2, r, xy_factors, &rxy);
    r_temp[0]=r[0];
    r_temp[1]=r[1];
    r_temp[2]=rxy;
    if(cplx == 1)
        res_temp=computeMI_knml(datafactors[0],datafactors[1],xy_factors,r_temp,n,c2terms, looklog, 0);
    else
        res_temp=computeMI_kmdl(datafactors[0],datafactors[1],xy_factors,r_temp,n, looklog, 0);


    /////////////////////////////////////
    if(ptr_cnt[ptrVarIdx[0]]==0 && ptr_cnt[ptrVarIdx[1]]==0){  //all discrete

        // free memory
        free(xy_factors);
        free(r_temp);
        free(ptr);

        for(l=0;l<(2);l++) free(datafactors[l]);
        free(datafactors);


        return res_temp;
    }

    // Find the best initial conditions with the same number of bins (equalfreq) on all continuous variables.
    double max_res=res_temp[1];
    int max_initbins=initbins;
    for(int new_initbins=2; (new_initbins<initbins) && (new_initbins<20); new_initbins++){

        int lbin=floor(n/new_initbins);
        if(lbin<1) {
            lbin=1;
            new_initbins=n;
        }

        // Reinitialization cut and r
        for(l=0;l<2;l++){
        if(ptr_cnt[ptrVarIdx[l]]==1){

            for(j=0;j<new_initbins-1;j++) {
                cut[l][j]=j*lbin+lbin-1;
            }
            cut[l][new_initbins-1]=n-1;
            r[l]=new_initbins;
        }else{
            r[l]=AllLevels[ptrVarIdx[l]];
        }
        }

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
        
        jointfactors_u(datafactors,ptr,n, 2, r, xy_factors, &rxy);
        r_temp[0]=r[0];
        r_temp[1]=r[1];
        r_temp[2]=rxy;
        if(cplx == 1)
            res_temp=computeMI_knml(datafactors[0],datafactors[1],xy_factors,r_temp,n,c2terms, looklog, 0);
        else
            res_temp=computeMI_kmdl(datafactors[0],datafactors[1],xy_factors,r_temp,n, looklog, 0);

        if(res_temp[1] > max_res){
            max_initbins = new_initbins;
            max_res = res_temp[1];
        }
        free(res_temp);
    }

    int lbin=floor(n/max_initbins);
    if(lbin<1) {
        lbin=1;
        max_initbins=n;
    }

    // Reinitialization cut and r
    for(l=0;l<2;l++){
    if(ptr_cnt[ptrVarIdx[l]]==1){

        for(j=0;j<max_initbins-1;j++) {
            cut[l][j]=j*lbin+lbin-1;
        }
        cut[l][max_initbins-1]=n-1;
        r[l]=max_initbins;
    }else{
        r[l]=AllLevels[ptrVarIdx[l]];
    }
    }

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


    // Run dynamic optimization with the best initial conditions.
    double sc;

    double* MI = (double *)calloc(STEPMAX,sizeof(double));
    double* MIk = (double *)calloc(STEPMAX,sizeof(double));
    double I_av,Ik_av;
     
    int stop,i,flag;

    // 

    int** factors1=(int **)calloc(3,sizeof(int*));
    int* singlefactor = (int *)calloc(n, sizeof(int));
    int* rt1=(int *)calloc(2,sizeof(int));
    int sc_levels1, sc_levels2;

    // 
    MI[0]=0;
    MIk[0]=-DBL_MAX;

    flag=0;
    int sc_levels_x = r[0]; // Number of levels of the first variable
    int sc_levels_y = r[1]; // Number of levels of the second variable
    int rx = r[0];
    int ry = r[1];
    int np; //number of possible cuts
    for(stop=1;stop<STEPMAX;stop++)
    {

        ///////////////////////////////////////////
        if(ptr_cnt[ptrVarIdx[0]]==1){

            //optimize on x
            //I(x;y)

            factors1[0]=datafactors[1]; //y
            factors1[1]=singlefactor;//One single bin at start
            factors1[2]=datafactors[1];//y

            rt1[0]=ry; //y
            rt1[1]=1; //One single bin at start

            //sc_levels_y = (sc_levels_y*(stop-1) + ry)/stop; //Harmonic mean of previous levels.
            //sc_levels_x = rx;
            sc_levels_y = ry;

            sc_levels1 = sc_levels_y;
            sc_levels2 = sc_levels_x;
            sc = 0.5*(sc_levels_y-1);

            rx = r[0];
            // Optimization run on X.
            optfun_onerun_kmdl_coarse(sortidx[ptrVarIdx[0]], data[ptrVarIdx[0]], 2, factors1, rt1,
                                      sc, sc_levels1, sc_levels2, n, AllLevels[ptrVarIdx[0]], cut[0], &(r[0]),
                                      maxbins, looklog, lookH, cterms, cplx, sample_weights, flag_effN); // 2 factors
        }

        ////////////////////////////////////////////////
        if(ptr_cnt[ptrVarIdx[1]]==1){
    
            //opt y
            //I(x;y)

            factors1[0]=datafactors[0];//x before its optimization
            factors1[1]=singlefactor;//One single bin at start
            factors1[2]=datafactors[0];//x

            rt1[0]=rx;//x before its optimization
            rt1[1]=1; //One single bin at start

            //sc_levels_x = (sc_levels_x*(stop-1) + rx)/stop;
            sc_levels_x = rx;
            //sc_levels_y = ry;

            sc_levels1 = sc_levels_x;
            sc_levels2 = sc_levels_y;
            sc = 0.5*(sc_levels_x-1);

            ry = r[1];
            // Optimization run on Y.
            optfun_onerun_kmdl_coarse(sortidx[ptrVarIdx[1]], data[ptrVarIdx[1]],2, factors1, rt1, sc,
                                      sc_levels1, sc_levels2, n, AllLevels[ptrVarIdx[1]], cut[1], &(r[1]),
                                      maxbins, looklog, lookH, cterms, cplx, sample_weights, flag_effN); // 2 factors
        }

        //////////////////////////////////////////
        // update both datafactors
        if(ptr_cnt[ptrVarIdx[0]]==1){
            update_datafactors(sortidx, ptrVarIdx[0], datafactors, 0, n, cut);
            rx = r[0];
        }
        if(ptr_cnt[ptrVarIdx[1]]==1){
            update_datafactors(sortidx, ptrVarIdx[1], datafactors, 1, n, cut);    
            ry = r[1];
        }

        jointfactors_u(datafactors,ptr,n, 2, r, xy_factors, &rxy);
        r_temp[0]=r[0];
        r_temp[1]=r[1];
        r_temp[2]=rxy;
        if(cplx == 1)
            res_temp=computeMI_knml(datafactors[0],datafactors[1],xy_factors,r_temp,n,c2terms,looklog, 0);
        else
            res_temp=computeMI_kmdl(datafactors[0],datafactors[1],xy_factors,r_temp,n,looklog, 0);
        //Adding combinatorial term
        if(ptr_cnt[ptrVarIdx[0]] == 1 && rx>1){
            np = min(maxbins, AllLevels[ptrVarIdx[0]]);
            res_temp[1] -= 0.5*(rx-1)*(log((1.0*np-1) / (rx-1) - 1) + 1)/n;
        }
        if(ptr_cnt[ptrVarIdx[1]] == 1 && ry>1){
            np = min(maxbins, AllLevels[ptrVarIdx[1]]);
            res_temp[1] -= 0.5*(ry-1)*(log((1.0*np-1) / (ry-1) - 1) + 1)/n;
        }

        for(i=stop-1;i>0;i--){
            if( fabs(res_temp[1]-MIk[i]) < EPS ){
                flag=1;
                Ik_av=MIk[i];
                I_av=MI[i];

                for(j=i+1;j<stop;j++){
                    Ik_av+=MIk[j];
                    I_av+=MI[j];
                }
                Ik_av/=(stop-i);//average over the periodic cycle
                I_av/=(stop-i);
                break;
            }
        }
        MIk[stop]=res_temp[1];
        MI[stop]=res_temp[0];
        free(res_temp);

        if( flag || (ptr_cnt[ptrVarIdx[0]]==0) || (ptr_cnt[ptrVarIdx[1]]==0) ){
            break;
        }
    }//for

    #if _MY_DEBUG_MInoU    
        printf("final: I_xy=%lf Ik_xy=%lf\n",res_temp[0],res_temp[1]);

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

    double* return_res = (double*) calloc(2, sizeof(double));
    if(flag){
        return_res[0] = I_av;
        return_res[1] = Ik_av;
    }
    else{
        return_res[0] = MI[stop];
        return_res[1] = MIk[stop];
    }

    // I and Ik can always be 0 by choosing 1 bin on either X or Y.
    if( (return_res[1] < 0) && ((ptr_cnt[ptrVarIdx[0]]==1) && (ptr_cnt[ptrVarIdx[1]]==1)) ){
        return_res[0] = 0;
        return_res[1] = 0;
    }

    // free memory
    free(xy_factors);
    free(r_temp);
    free(ptr);

    for(l=0;l<(2);l++) free(datafactors[l]);
    free(datafactors);

    free(factors1);
    free(rt1);
    free(MI);
    free(MIk);
    free(singlefactor);

    return return_res;
}




double* compute_Ixy_cond_u_new_alg1(int** data, int** sortidx, int* ptr_cnt, int* ptrVarIdx, int* AllLevels, 
                                    int nbrUi, int n, int maxbins, int initbins, int **cut, int *r, double* c2terms,
                                    double**cterms, int lbin, double* looklog, double* lookH, int cplx, 
                                    double* sample_weights, bool flag_effN)
{

    int j,l,ll;
    int STEPMAX1=50;

    /////////////////////////////////////

    double* res_temp=(double *)calloc(2,sizeof(double));//res_tempults res_temp[0]->I,res_temp[1]->I-k

    //////////////////////////////////
    //allocation factors  x y

    int **datafactors=(int **)calloc((nbrUi+2),sizeof(int*));
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

    double* MI = (double*)calloc(STEPMAX, sizeof(double));
    double* MIk = (double*)calloc(STEPMAX, sizeof(double));
    double I_av,Ik_av;
    double* MI1 = (double*)calloc(STEPMAX1, sizeof(double));
    double* MIk1 = (double*)calloc(STEPMAX1, sizeof(double));
    double I_av1,Ik_av1;
    MI[0]=0;
    MIk[0]=-DBL_MAX;

    // Loop variables
    int stop1,stop,i,flag,flag1;
    // Mutual informations for computing I(x;y|u).
    double I_x_yu;
    double I_y_xu;
    double I_x_u;
    double I_y_u;
    double cond_I;
    double Ik_x_yu;
    double Ik_y_xu;
    double Ik_x_u;
    double Ik_y_u;
    double cond_Ik;

    // Complexity factor (# of levels) passed to each optimization
    int sc_levels1, sc_levels2;


    int **factors1=(int **)calloc(3,sizeof(int*));

    // Find the best initial conditions with the same number of bins (equalfreq) on all continuous variables.
    double max_res=res_temp[1];
    int max_initbins=initbins;
    for(int new_initbins=2; (new_initbins<initbins) && (new_initbins<20); new_initbins++){

        //////////////////////////////////////
        // initialization cut and r
        for(l=0;l<(nbrUi+2);l++){
            if(ptr_cnt[ptrVarIdx[l]]==1){

                for(j=0;j<new_initbins-1;j++) {
                    cut[l][j]=j*lbin+lbin-1;
                }
                cut[l][new_initbins-1]=n-1;
                r[l]=new_initbins;
            }else{
                r[l]=AllLevels[ptrVarIdx[l]];
            }
        }

        for(l=0;l<(nbrUi+2);l++){ //compute datafactors based on the positions of cut points in vector <cut>
            if(ptr_cnt[ptrVarIdx[l]]==1){
                update_datafactors(sortidx, ptrVarIdx[l], datafactors, l, n, cut);
            }else{
                for(j=0;j<=n-1;j++){//discrete case
                    datafactors[l][j]=data[ptrVarIdx[l]][j];
                }
            }
        }

        jointfactors_uiyx(datafactors, -1, n, nbrUi, r, uiyxfactors, ruiyx);
        r_temp[0]=r[1];
        r_temp[1]=ruiyx[2];
        r_temp[2]=ruiyx[3];
        if(cplx == 1)
            res_temp=computeMI_knml(datafactors[1],uiyxfactors[2], uiyxfactors[3],r_temp,n,c2terms, looklog, FLAG_CPLX);
        else
            res_temp=computeMI_kmdl(datafactors[1],uiyxfactors[2], uiyxfactors[3],r_temp,n, looklog, 0);
        I_y_xu = res_temp[0]; // Before optimization on X.
        Ik_y_xu = res_temp[1];
        free(res_temp);

        r_temp[0]=r[0];
        r_temp[1]=ruiyx[1];
        r_temp[2]=ruiyx[3];
        if(cplx == 1)
            res_temp=computeMI_knml(datafactors[0],uiyxfactors[1], uiyxfactors[3],r_temp,n,c2terms, looklog, FLAG_CPLX);
        else
            res_temp=computeMI_kmdl(datafactors[0],uiyxfactors[1], uiyxfactors[3],r_temp,n, looklog, 0);
        I_x_yu = res_temp[0]; //Before updating Y (and X).
        Ik_x_yu = res_temp[1];
        free(res_temp);


        if( (Ik_y_xu + Ik_x_yu) > max_res){
            max_initbins = new_initbins;
            max_res = (Ik_y_xu + Ik_x_yu);
        }
    }

    lbin=floor(n/max_initbins);
    if(lbin<1) {
        lbin=1;
        max_initbins=n;
    }

    //////////////////////////////////////
    // initialization cut and r
    for(l=0;l<(nbrUi+2);l++){
        if(ptr_cnt[ptrVarIdx[l]]==1){

            for(j=0;j<max_initbins-1;j++) {
                cut[l][j]=j*lbin+lbin-1;
            }
            cut[l][max_initbins-1]=n-1;
            r[l]=max_initbins;
        }else{
            r[l]=AllLevels[ptrVarIdx[l]];
        }
    }

    for(l=0;l<(nbrUi+2);l++){ //compute datafactors based on the positions of cut points in vector <cut>
        if(ptr_cnt[ptrVarIdx[l]]==1){
            update_datafactors(sortidx, ptrVarIdx[l], datafactors, l, n, cut);
        }else{
            for(j=0;j<=n-1;j++){//discrete case
                datafactors[l][j]=data[ptrVarIdx[l]][j];
            }
        }
    }

    // Run optimization with best initial equal freq.
    int *rt1=(int *)calloc(2,sizeof(int));
    int *r_old = (int *)calloc((nbrUi+2),sizeof(int));
    for(int ll=0; ll<(nbrUi+2); ll++){
        r_old[ll] = r[ll];
    }
    int U_counter;


    int sc_levels_x; // Number of levels of the first variable
    int sc_levels_y; // Number of levels of the second variable
    flag1=0;
    int max_U_counter = 3;
    int np; //number of possible cuts for combinatorial term
    for(stop1=1;stop1<STEPMAX1;stop1++)
    {

        ///////////////////////////////////////////
        //optimize I(y;xu) over x and u
        U_counter=0;
        while(U_counter < max_U_counter) {
            for(l=0;l<nbrUi;l++){

                if(ptr_cnt[ptrVarIdx[l+2]]==1){
                    //opt u
                    //I(y;xu)
                    jointfactors_uiyx(datafactors, l+2,  n, nbrUi, r_old, uiyxfactors, ruiyx);

                    //init variables for the optimization run
                    factors1[0]=uiyxfactors[3];//xyu
                    factors1[1]=uiyxfactors[2];//xu
                    factors1[2]=datafactors[1];//y

                    rt1[0]=ruiyx[3];//xyu
                    rt1[1]=ruiyx[2];//xu
                    sc_levels_x = r_old[0];
                    sc_levels_y = r_old[1];

                    sc = 0.5*(sc_levels_y-1)*ruiyx[2];

                    sc_levels1 = sc_levels_y;
                    sc_levels2 = r_old[l+2]; //old nlevels for combinatorial term

                    // Run optimization on U.
                    optfun_onerun_kmdl_coarse(sortidx[ptrVarIdx[l+2]], data[ptrVarIdx[l+2]], 2, factors1, rt1,
                                              sc, sc_levels1, sc_levels2, n, AllLevels[ptrVarIdx[l+2]], cut[l+2],
                                              &(r[l+2]), maxbins, looklog, lookH, cterms, cplx, sample_weights, flag_effN); // 2 factors

                    //update_datafactors(sortidx, ptrVarIdx[l+2], datafactors, l+2, n, cut); //moved outside of U loop
                }
            }//for all Uis
            for(int ll=0; ll<nbrUi; ll++){
                if(ptr_cnt[ptrVarIdx[ll+2]]==1) update_datafactors(sortidx, ptrVarIdx[ll+2], datafactors, ll+2, n, cut);
                r_old[ll+2] = r[ll+2]; 
            }
            U_counter++;
            if(nbrUi==1) U_counter = max_U_counter;
        }//U_counter loop

        jointfactors_uiyx(datafactors, -1, n, nbrUi, r_old, uiyxfactors, ruiyx);
        r_temp[0]=r_old[1];
        r_temp[1]=ruiyx[2];
        r_temp[2]=ruiyx[3];
        if(cplx == 1)
            res_temp=computeMI_knml(datafactors[1],uiyxfactors[2], uiyxfactors[3],r_temp,n,c2terms, looklog, FLAG_CPLX);
        else
            res_temp=computeMI_kmdl(datafactors[1],uiyxfactors[2], uiyxfactors[3],r_temp,n, looklog, 0);
        I_y_xu = res_temp[0]; // Before optimization on X.
        Ik_y_xu = res_temp[1];
        free(res_temp);
        //// Adding combinatorial term
        for(ll=0;(ll<nbrUi); ll++){
            np = min(AllLevels[ptrVarIdx[ll+2]], maxbins);
            if((ptr_cnt[ptrVarIdx[ll+2]]==1) && (r_old[ll+2]>1))
                Ik_y_xu -= (r_old[ll+2]-1)*(log((1.0*np-1) / (r_old[ll+2]-1) - 1) + 1)/n;
        }
        if((ptr_cnt[ptrVarIdx[0]]==1) && (r_old[0]>1)){
            np = min(AllLevels[ptrVarIdx[0]], maxbins);
            Ik_y_xu -= (r_old[0]-1)*(log((1.0*np-1) / (r_old[0]-1) - 1) + 1)/n;
        }



        if(ptr_cnt[ptrVarIdx[0]]==1){
            //opt x
            //I(y;xu)

            //compute joint factors u yu xu xyu
            jointfactors_uiyx(datafactors, -1, n, nbrUi, r_old, uiyxfactors, ruiyx);
            //init variables for the optimization run
            factors1[0]=uiyxfactors[1];//uy
            factors1[1]=uiyxfactors[0];//u
            factors1[2]=datafactors[1];//y

            rt1[0]=ruiyx[1];//uy
            rt1[1]=ruiyx[0];//u
            sc_levels_x = r_old[0];
            sc_levels_y = r_old[1];
            sc = 0.5*(sc_levels_y-1)*ruiyx[0];

            #if _MY_DEBUG_MInoU
                printf("start optfun\n ");
                fflush(stdout);
            #endif

            sc_levels1 = sc_levels_y; //herve
            sc_levels2 = sc_levels_x; //herve

            // Run optimization on X.
            optfun_onerun_kmdl_coarse(sortidx[ptrVarIdx[0]], data[ptrVarIdx[0]], 2, factors1, rt1, 
                                      sc, sc_levels1, sc_levels2, n, AllLevels[ptrVarIdx[0]], cut[0], 
                                      &(r[0]), maxbins, looklog, lookH, cterms, cplx, sample_weights, flag_effN); // 2 factors

            //update_datafactors(sortidx, ptrVarIdx[0], datafactors, 0, n, cut); //moved to after Y opt
        }

        // Reset cutpoints on U
        reset_u_cutpoints(cut, nbrUi, ptr_cnt, ptrVarIdx, initbins, maxbins, lbin, r, AllLevels, n);
        for(l=0; l<nbrUi; l++){ 
            if(ptr_cnt[ptrVarIdx[l+2]]==1) update_datafactors(sortidx, ptrVarIdx[l+2], datafactors, l+2, n, cut);
            r_old[l+2] = r[l+2];//r[l+2] is set to init_nbins during reset_u_cutpoints
        }

        ///////////////////////////////////////////
        //optimize I(x;yu) over y and u
        U_counter=0;
        while(U_counter < max_U_counter) {
            for(l=0;l<nbrUi;l++){

                if(ptr_cnt[ptrVarIdx[l+2]]==1){
                    //opt u
                    //I(x;yu)
                    jointfactors_uiyx(datafactors, l+2,  n, nbrUi, r_old, uiyxfactors, ruiyx);

                    //init variables for the optimization run
                    factors1[0]=uiyxfactors[3];//xyu
                    factors1[1]=uiyxfactors[1];//yu
                    factors1[2]=datafactors[0];//x

                    rt1[0]=ruiyx[3]; //xyu
                    rt1[1]=ruiyx[1]; //yu
                    sc_levels_y = r_old[1];
                    sc_levels_x = r_old[0];

                    sc_levels1 = sc_levels_x; //herve
                    sc_levels2 = r_old[l+2]; //herve
                    sc = 0.5*(sc_levels_x-1)*ruiyx[1];

                    // Run optimization on U.
                    optfun_onerun_kmdl_coarse(sortidx[ptrVarIdx[l+2]], data[ptrVarIdx[l+2]], 2, factors1, rt1,
                                              sc, sc_levels1, sc_levels2, n, AllLevels[ptrVarIdx[l+2]], cut[l+2],
                                              &(r[l+2]), maxbins, looklog, lookH, cterms, cplx, sample_weights, flag_effN); // 2 factors

                    //update_datafactors(sortidx, ptrVarIdx[l+2], datafactors, l+2, n, cut); //moved to outside U loop
                }
            }//for all Uis
            for(int ll=0; ll<nbrUi; ll++){
                if(ptr_cnt[ptrVarIdx[ll+2]]==1) update_datafactors(sortidx, ptrVarIdx[ll+2], datafactors, ll+2, n, cut);
                r_old[ll+2] = r[ll+2];
            }
            U_counter++;
            if(nbrUi==1) U_counter = max_U_counter;
        }//U_counter loop

        jointfactors_uiyx(datafactors, -1, n, nbrUi, r_old, uiyxfactors, ruiyx);
        r_temp[0]=r_old[0];
        r_temp[1]=ruiyx[1];
        r_temp[2]=ruiyx[3];
        if(cplx == 1)
            res_temp=computeMI_knml(datafactors[0],uiyxfactors[1], uiyxfactors[3],r_temp,n,c2terms, looklog, FLAG_CPLX);
        else
            res_temp=computeMI_kmdl(datafactors[0],uiyxfactors[1], uiyxfactors[3],r_temp,n, looklog, 0);
        I_x_yu = res_temp[0]; //Before updating Y (and X).
        Ik_x_yu = res_temp[1];
        free(res_temp);
        //// Adding combinatorial term
        for(ll=0;(ll<nbrUi); ll++){
            np = min(AllLevels[ptrVarIdx[ll+2]], maxbins);
            if((ptr_cnt[ptrVarIdx[ll+2]]==1) && (r_old[ll+2]>1))
                Ik_x_yu -= (r_old[ll+2]-1)*(log((1.0*np-1) / (r_old[ll+2]-1) - 1) + 1)/n;
        }
        if((ptr_cnt[ptrVarIdx[1]]==1) && (r_old[1]>1)) {
            np = min(AllLevels[ptrVarIdx[1]], maxbins);
            Ik_x_yu -= (r_old[1]-1)*(log((1.0*np-1) / (r_old[1]-1) - 1) + 1)/n;
        }


        if(ptr_cnt[ptrVarIdx[1]]==1){
            //optimize on y
            //I(x;yu)
            jointfactors_uiyx(datafactors, -1, n, nbrUi, r_old, uiyxfactors, ruiyx);

            //init variables for the optimization run
            factors1[0]=uiyxfactors[2];//ux
            factors1[1]=uiyxfactors[0];//u
            factors1[2]=datafactors[0];//x

            rt1[0]=ruiyx[2]; //ux
            rt1[1]=ruiyx[0]; //u

            sc_levels_x = r_old[0];
            sc_levels_y = r_old[1];

            sc_levels1 = sc_levels_x; //herve
            sc_levels2 = sc_levels_y; //herve
            sc = 0.5*(sc_levels_x-1)*ruiyx[0];

            // Run optimization on Y.
            optfun_onerun_kmdl_coarse(sortidx[ptrVarIdx[1]], data[ptrVarIdx[1]], 2, factors1, rt1, 
                                      sc, sc_levels1, sc_levels2, n, AllLevels[ptrVarIdx[1]], cut[1], 
                                      &(r[1]), maxbins, looklog, lookH, cterms, cplx, sample_weights, flag_effN); // 2 factors

            //update_datafactors(sortidx, ptrVarIdx[1], datafactors, 1, n, cut); //moved to end of loop1
        }

        // Reset cutpoints on U
        reset_u_cutpoints(cut, nbrUi, ptr_cnt, ptrVarIdx, initbins, maxbins, lbin, r, AllLevels, n);
        for(l=0; l<nbrUi; l++){
            if(ptr_cnt[ptrVarIdx[l+2]]==1) update_datafactors(sortidx, ptrVarIdx[l+2], datafactors, l+2, n, cut);
            r_old[l+2] = r[l+2];
        }

        ///////////////////////////////////////////
        //optimize I(x;u) over u
        U_counter=0;
        while(U_counter < max_U_counter) {
            for(l=0;l<nbrUi;l++){

                if(ptr_cnt[ptrVarIdx[l+2]]==1){
                    jointfactors_uiyx(datafactors, l+2,  n, nbrUi, r_old, uiyxfactors, ruiyx);

                    //init variables for the optimization run
                    factors1[0]=uiyxfactors[2];//xu
                    factors1[1]=uiyxfactors[0];//u
                    factors1[2]=datafactors[0];//x

                    rt1[0]=ruiyx[2];//xu
                    rt1[1]=ruiyx[0];//u

                    sc_levels1 = r_old[0];//x
                    sc_levels2 = r_old[l+2];//u
                    sc = 0.5*(sc_levels_x-1)*ruiyx[0];

                    //optimization run on ptrVarIdx[l+2]
                    optfun_onerun_kmdl_coarse(sortidx[ptrVarIdx[l+2]], data[ptrVarIdx[l+2]], 2, factors1, rt1,
                                              sc, sc_levels1, sc_levels2, n, AllLevels[ptrVarIdx[l+2]], cut[l+2],
                                              &(r[l+2]), maxbins, looklog, lookH, cterms, cplx, sample_weights, flag_effN); // 2 factors //herve

                    //update_datafactors(sortidx, ptrVarIdx[l+2], datafactors, l+2, n, cut);
                    }
            }//for all Uis
            for(int ll=0; ll<nbrUi; ll++){
                if(ptr_cnt[ptrVarIdx[ll+2]]==1) update_datafactors(sortidx, ptrVarIdx[ll+2], datafactors, ll+2, n, cut);
                r_old[ll+2] = r[ll+2];
            }
            U_counter++;
            if(nbrUi==1) U_counter = max_U_counter;
        }//U_counter loop

        jointfactors_uiyx(datafactors, -1, n, nbrUi, r_old, uiyxfactors, ruiyx);
        r_temp[0]=r_old[0];
        r_temp[1]=ruiyx[0];
        r_temp[2]=ruiyx[2];
        if(cplx == 1)
            res_temp=computeMI_knml(datafactors[0],uiyxfactors[0], uiyxfactors[2],r_temp,n,c2terms, looklog, FLAG_CPLX);
        else
            res_temp=computeMI_kmdl(datafactors[0],uiyxfactors[0], uiyxfactors[2],r_temp,n, looklog, 0);
        I_x_u = res_temp[0]; //After optimization on U.
        Ik_x_u = res_temp[1];
        free(res_temp);
        // Adding combinatorial term
        for(ll=0;(ll<nbrUi); ll++){
            np = min(AllLevels[ptrVarIdx[ll+2]], maxbins);
            if((ptr_cnt[ptrVarIdx[ll+2]]==1) && (r_old[ll+2]>1))
                Ik_x_u -= (r_old[ll+2]-1)*(log((1.0*np-1) / (r_old[ll+2]-1) - 1) + 1)/n;
        }


        // Reset cutpoints on U
        reset_u_cutpoints(cut, nbrUi, ptr_cnt, ptrVarIdx, initbins, maxbins, lbin, r, AllLevels, n);
        for(l=0; l<nbrUi; l++){ 
            if(ptr_cnt[ptrVarIdx[l+2]]==1) update_datafactors(sortidx, ptrVarIdx[l+2], datafactors, l+2, n, cut);
            r_old[l+2] = r[l+2];
        }


        ///////////////////////////////////////////
        //optimize I(y;u) over u
        U_counter=0;
        while(U_counter < max_U_counter) {
            for(l=0;l<nbrUi;l++){

                if(ptr_cnt[ptrVarIdx[l+2]]==1){
                    jointfactors_uiyx(datafactors, l+2,  n, nbrUi, r_old, uiyxfactors, ruiyx);

                    factors1[0]=uiyxfactors[1];//yu
                    factors1[1]=uiyxfactors[0];//u
                    factors1[2]=datafactors[1];//y

                    rt1[0]=ruiyx[1];//yu
                    rt1[1]=ruiyx[0];//u

                    sc_levels1 = r_old[1];//y
                    sc_levels2 = r_old[l+2];//u
                    sc = 0.5*(sc_levels1-1)*ruiyx[0];

                    //optimization run on ptrVarIdx[l+2]
                    optfun_onerun_kmdl_coarse(sortidx[ptrVarIdx[l+2]], data[ptrVarIdx[l+2]], 2, factors1, rt1,
                                              sc, sc_levels1, sc_levels2, n, AllLevels[ptrVarIdx[l+2]], cut[l+2],
                                              &(r[l+2]), maxbins, looklog, lookH, cterms, cplx, sample_weights, flag_effN); //2 factors

                    //update_datafactors(sortidx, ptrVarIdx[l+2], datafactors, l+2, n, cut);
                }
            }//for all Uis
            for(int ll=0; ll<nbrUi; ll++){
                if(ptr_cnt[ptrVarIdx[ll+2]]==1) update_datafactors(sortidx, ptrVarIdx[ll+2], datafactors, ll+2, n, cut);
                r_old[ll+2] = r[ll+2];
            }
            U_counter++;
            if(nbrUi==1) U_counter = max_U_counter;
        }//U_counter loop

        jointfactors_uiyx(datafactors, -1, n, nbrUi, r_old, uiyxfactors, ruiyx);
        r_temp[0]=r_old[1];
        r_temp[1]=ruiyx[0];
        r_temp[2]=ruiyx[1];
        if(cplx == 1)
            res_temp=computeMI_knml(datafactors[1],uiyxfactors[0], uiyxfactors[1],r_temp,n,c2terms, looklog, FLAG_CPLX);
        else
            res_temp=computeMI_kmdl(datafactors[1],uiyxfactors[0], uiyxfactors[1],r_temp,n, looklog, 0);
        I_y_u = res_temp[0]; //After optimization on U.
        Ik_y_u = res_temp[1];
        free(res_temp);
        // Adding combinatorial term
        for(ll=0;(ll<nbrUi); ll++){
            np = min(AllLevels[ptrVarIdx[ll+2]], maxbins);
            if((ptr_cnt[ptrVarIdx[ll+2]]==1) && (r_old[ll+2]>1))
                Ik_y_u -= (r_old[ll+2]-1)*(log((1.0*np-1) / (r_old[ll+2]-1) - 1) + 1)/n;
        }


        // Reset cutpoints on U
        reset_u_cutpoints(cut, nbrUi, ptr_cnt, ptrVarIdx, initbins, maxbins, lbin, r, AllLevels, n);
        for(l=0; l<nbrUi; l++){ 
            if(ptr_cnt[ptrVarIdx[l+2]]==1) update_datafactors(sortidx, ptrVarIdx[l+2], datafactors, l+2, n, cut);
            r_old[l+2] = r[l+2];
        }

        //Update X and Y
        if(ptr_cnt[ptrVarIdx[0]] == 1) {
            update_datafactors(sortidx, ptrVarIdx[0], datafactors, 0, n, cut);
            r_old[0] = r[0];
        }
        if(ptr_cnt[ptrVarIdx[1]] == 1) {
            update_datafactors(sortidx, ptrVarIdx[1], datafactors, 1, n, cut);
            r_old[1] = r[1];
        }

        //#################
        // Compute I(X;Y|U)
        cond_I  = 0.5* (I_x_yu - I_x_u + I_y_xu - I_y_u);
        cond_Ik = 0.5* (Ik_x_yu - Ik_x_u + Ik_y_xu - Ik_y_u);
        //printf("0.5*(%.2f - %.2f + %.2f - %.2f) = %.2f\n", I_x_yu, I_x_u, I_y_xu, I_y_u, res_temp[0]);

        // Test stop condition on stop1
        for(i=stop1-1;i>0;i--){
            if( fabs(cond_Ik-MIk1[i]) < EPS ) { // If no real improvement over last information
                flag1=1;
                Ik_av1=MIk1[i];
                I_av1=MI1[i];

                for(j=i+1;j<stop1;j++){
                    Ik_av1+=MIk1[j];
                    I_av1+=MI1[j];
                }
                Ik_av1/=(stop1-i); //average over the periodic cycle
                I_av1/=(stop1-i);
                break;
            }
        }
        if(flag1) break;
        MIk1[stop1] = cond_Ik;
        MI1[stop1]  = cond_I;

    }//end stop1
    double* return_res = (double*)calloc(2,sizeof(double));
    if(flag1){
        return_res[0] = I_av1;
        return_res[1] = Ik_av1;
    }
    else{
        return_res[0] = cond_I;
        return_res[1] = cond_Ik;
    }

    // I and Ik can always be 0 by choosing 1 bin on either X or Y.
    if( (return_res[1] < 0) && ((ptr_cnt[ptrVarIdx[0]]==1) && (ptr_cnt[ptrVarIdx[1]]==1)) ){
        return_res[0] = 0;
        return_res[1] = 0;
    }

    // free memory
    free(r_temp);

    for(l=0;l<(nbrUi+2);l++) free(datafactors[l]);
    free(datafactors);
    for(l=0;l<4;l++) free(uiyxfactors[l]);
    free(uiyxfactors);
    free(factors1);
    free(rt1);
    free(r_old);

    free(ruiyx);
    free(MI);
    free(MIk);
    free(MI1);
    free(MIk1);

    return return_res;
}


double* compute_mi_cond_alg1(int** data, int** sortidx,  int* AllLevels, int* ptr_cnt, int* ptrVarIdx, int nbrUi,
                             int n, int maxbins, int initbins, double* c2terms, double** cterms, double* looklog,
                             double* lookH, int cplx, double* sample_weights, bool flag_effN)
{
    double* res= new double[3]();//results res[0]->I,res[1]->I-k
    double* res_temp;//=(double *)calloc(2,sizeof(double));//results res[0]->I,res[1]->I-k
    
    int j,k,l,ll;
    int np;    // int np=ceil(1.0*n/coarse);

    int *r=(int *)calloc((nbrUi+2),sizeof(int));

    int **cut;
    cut=(int **)calloc(nbrUi+2,sizeof(int*));
    for(l=0;l<(nbrUi+2);l++){
        cut[l]=(int *)calloc(maxbins,sizeof(int));
    }

    int lbin=floor(n/initbins);
    if(lbin<1) {
        lbin=1;
        initbins=n;
    }

    // ///////////////////////////////////////////////////////////////////////////////////
    // no conditioning, empty set of variables in u 
    // calling funcion compute_Ixy_alg1
    // NO u 
    if(nbrUi==0){

        //////////////////////////////////////
        // initialization cut and r
        for(l=0;l<(nbrUi+2);l++){
            if(ptr_cnt[ptrVarIdx[l]]==1){

                for(j=0;j<initbins-1;j++) {
                    cut[l][j]=j*lbin+lbin-1;
                }
                cut[l][initbins-1]=n-1;
                r[l]=initbins;
            }else{
                r[l]=AllLevels[ptrVarIdx[l]];
            }
        }

        res_temp=compute_Ixy_alg1(data, sortidx, ptr_cnt, ptrVarIdx, AllLevels, n, maxbins, initbins, 
                                  cut, r, c2terms, cterms, looklog, lookH, cplx, sample_weights, flag_effN);

        res[0]=n;
        res[1]=res_temp[0];
        res[2]=res_temp[0]-res_temp[1];

        // free memory
        free(res_temp);
        free(r);
        for(l=0;l<(nbrUi+2);l++){
            free(cut[l]);
        }
        free(cut);

        return res;
    }

    else { // with U

        //////////////////////////////////////
        // initialization cut and r
        for(l=0;l<(nbrUi+2);l++){
            if(ptr_cnt[ptrVarIdx[l]]==1){

                for(j=0;j<initbins-1;j++) {
                    cut[l][j]=j*lbin+lbin-1;
                }
                cut[l][initbins-1]=n-1;
                r[l]=initbins;
            }else{
                r[l]=AllLevels[ptrVarIdx[l]];
            }
        }

        res_temp=compute_Ixy_cond_u_new_alg1(data, sortidx, ptr_cnt, ptrVarIdx, AllLevels, nbrUi, n, 
                                             maxbins, initbins, cut, r, c2terms, cterms,
                                             lbin, looklog, lookH, cplx, sample_weights, flag_effN);

        res[0]=n;
        res[1]=res_temp[0];
        res[2]=res_temp[0]-res_temp[1];

        // free memory
        free(res_temp);
        free(r);
        for(l=0;l<(nbrUi+2);l++){
            free(cut[l]);
        }
        free(cut);

        return res;
    }
}


// //////////////////////////////////////////////////////////////////////////////////////////77
// compute Rscore and three point mutual information I(x;y;z | u)
//input x y z u-> compute Rscore and Ixyz
//
//
// Returns:
//res[0]=Rscore
//res[1]=N*Ixyz
//res[2]=N*kxyz


double* compute_Rscore_Ixyz_new_alg5(int** data, int** sortidx, int* AllLevels, int* ptr_cnt, int* ptrVarIdx,
                                     int nbrUi, int ptrZiIdx, int n, int maxbins, int initbins, double* c2terms, 
                                     double** cterms, double* looklog, double* lookH, int cplx, double* sample_weights,
                                     bool flag_effN)
{

    double Rscore;
    double nv,dpi,first,second,xz,yz;

    int j,k,l,ll;

    double I_xy_u,I_xz_u,I_yz_u,I_xy_zu;
    double Ik_xy_u,Ik_xz_u,Ik_yz_u,Ik_xy_zu;
    double I_xyz_u,Ik_xyz_u;

    double* res_temp=(double *)calloc(2,sizeof(double));//results res[0]->I,res[1]->I-k
    double* res = new double[3]();//results res[0]->I,res[1]->I-k

    int* ptrVarIdx_t = (int *) calloc((nbrUi +2),sizeof(int));

    int* r=(int *)calloc((nbrUi+3),sizeof(int));
    int** cut;
    cut=(int **)calloc(nbrUi+3,sizeof(int*));
    for(l=0;l<(nbrUi+3);l++){
        cut[l]=(int *)calloc(maxbins,sizeof(int));
    }

    int lbin=floor(n/initbins);
    if(lbin<1) {
        lbin=1;
        initbins=n;
    }


    //initialitize cuts vectors
    for(l=0;l<(nbrUi+2);l++){
        if(ptr_cnt[ptrVarIdx[l]]==1){
            for(j=0;j<initbins-1;j++) {
                cut[l][j]=j*lbin+lbin-1;
            }
            cut[l][initbins-1]=n-1;
            r[l]=initbins;
        }else{
            r[l]=AllLevels[ptrVarIdx[l]];
        }
    }
    //z
    l=nbrUi+2;
    if(ptr_cnt[ptrZiIdx]==1){
        for(j=0;j<initbins-1;j++) {
            cut[l][j]=j*lbin+lbin-1;
        }
        cut[l][initbins-1]=n-1;
        r[l]=initbins;
    }else{
        r[l]=AllLevels[ptrZiIdx];
    }

    int **datafactors=(int **)calloc((nbrUi+3),sizeof(int*));
    for(l=0;l<(nbrUi+3);l++){
        datafactors[l]=(int *)calloc(n,sizeof(int));
    }
    int *ptr_u2_t=(int *)calloc(nbrUi+2,sizeof(int));
    int **factors4_t=(int **)calloc(4,sizeof(int*));
    for(l=0;l<4;l++){
        factors4_t[l]=(int *)calloc(n,sizeof(int));
    }
    int *r4_t=(int *)calloc(4,sizeof(int));
    int* r_temp=(int *)calloc(2,sizeof(int));

    // if opt
    int **cut_t = (int **)calloc(nbrUi+2,sizeof(int*));
    int *r_t = (int *)calloc((nbrUi+2),sizeof(int));
    for(l=0;l<(nbrUi+2);l++){
        cut_t[l]=(int *)calloc(maxbins,sizeof(int));
    }

    //if(flag_opt==0){
    //    // For now, compute all MIs with equal freq cutpoints

    //    // Fill datafactors
    //    for(l=0;l<(nbrUi+2);l++){
    //        if(ptr_cnt[ptrVarIdx[l]]==1){
    //            update_datafactors(sortidx, ptrVarIdx[l], datafactors, l, n, cut);    
    //        }else{
    //            for(j=0;j<=n-1;j++){
    //                datafactors[l][j]=data[ptrVarIdx[l]][j];
    //            }
    //        }
    //    }
    //    //z
    //    l=nbrUi+2;
    //    if(ptr_cnt[ptrZiIdx]==1){
    //        update_datafactors(sortidx, ptrZiIdx, datafactors, l, n, cut);    
    //    }else{
    //        for(j=0;j<=n-1;j++){
    //            datafactors[l][j]=data[ptrZiIdx][j];
    //        }
    //    }

    //    //I(x;y|u,z)
    //    jointfactors_uiyx(datafactors,-1,  n, nbrUi+1, r, factors4_t, r4_t);
    //    //compute MI and knml complexity 
    //    if(cplx == 1)
    //        res_temp=computeMIcond_knml(factors4_t,r4_t,r,n,c2terms, looklog);
    //    else
    //        res_temp=computeMIcond_kmdl(factors4_t,r4_t,r,n, looklog);
    //    I_xy_zu=res_temp[0];
    //    Ik_xy_zu=res_temp[1];
    //    free(res_temp);


    //    //I(x;y|u)
    //    ptr_u2_t[0]=0;
    //    ptr_u2_t[1]=1;
    //    for(ll=0;ll<nbrUi;ll++) ptr_u2_t[ll+2]=ll+2;

    //    jointfactors_uyx(datafactors,ptr_u2_t, n, nbrUi, r,factors4_t, r4_t);

    //    r_temp[0]=r[0];
    //    r_temp[1]=r[1];
    //    if(cplx == 1)
    //        res_temp=computeMIcond_knml(factors4_t,r4_t,r_temp,n,c2terms, looklog);
    //    else
    //        res_temp=computeMIcond_kmdl(factors4_t,r4_t,r_temp,n, looklog);
    //    I_xy_u=res_temp[0];
    //    Ik_xy_u=res_temp[1];
    //    free(res_temp);


    //    //I(z;x|u)
    //    ptr_u2_t[0]=nbrUi+2;
    //    ptr_u2_t[1]=0;
    //    for(ll=0;ll<nbrUi;ll++) ptr_u2_t[ll+2]=ll+2;

    //    jointfactors_uyx(datafactors,ptr_u2_t, n, nbrUi, r, factors4_t, r4_t);

    //    r_temp[0]=r[nbrUi+2];
    //    r_temp[1]=r[0];
    //    if(cplx == 1)
    //        res_temp=computeMIcond_knml(factors4_t,r4_t,r_temp,n,c2terms, looklog);
    //    else
    //        res_temp=computeMIcond_kmdl(factors4_t,r4_t,r_temp,n, looklog);
    //    I_xz_u=res_temp[0];
    //    Ik_xz_u=res_temp[1];
    //    free(res_temp);


    //    //I(z;y|u)
    //    ptr_u2_t[1]=1;
    //    
    //    jointfactors_uyx(datafactors,ptr_u2_t, n, nbrUi, r, factors4_t, r4_t);

    //    r_temp[0]=r[nbrUi+2];
    //    r_temp[1]=r[1];
    //    if(cplx == 1)
    //        res_temp=computeMIcond_knml(factors4_t,r4_t,r_temp,n,c2terms, looklog);
    //    else
    //        res_temp=computeMIcond_kmdl(factors4_t,r4_t,r_temp,n, looklog);
    //    I_yz_u=res_temp[0];
    //    Ik_yz_u=res_temp[1];
    //    free(res_temp);

    //}

    {
        // Optimize variables for each MI estimation for the R score

        // ######################
        //I(x,y|u,z)
        res_temp=compute_Ixy_cond_u_new_alg1(data, sortidx, ptr_cnt, ptrVarIdx, AllLevels, nbrUi+1, n, maxbins, //
                                            initbins, cut, r, c2terms, cterms, lbin, looklog, lookH, cplx, 
                                            sample_weights, flag_effN);
        I_xy_zu=res_temp[0];
        Ik_xy_zu=res_temp[1];
        free(res_temp);
        

        // ######################
        //I(x,y|u)

        for(l=0;l<(nbrUi+2);l++){
            if(ptr_cnt[ptrVarIdx[l]]==1){
                for(j=0;j<initbins-1;j++) {
                    cut_t[l][j]=j*lbin+lbin-1;
                }
                cut_t[l][initbins-1]=n-1;
                for(int j=initbins;j<maxbins;j++) {
                    cut_t[l][j]=0;
                }
                r_t[l]=initbins;
            }else{
                r_t[l]=AllLevels[ptrVarIdx[l]];
            }
        }
        // Do opt run on I(X;Y|U)
        if(nbrUi>0){
            res_temp=compute_Ixy_cond_u_new_alg1(data, sortidx, ptr_cnt, ptrVarIdx, AllLevels, nbrUi, n, maxbins,
                                                initbins, cut_t, r_t, c2terms, cterms, lbin, looklog, lookH, cplx, 
                                                sample_weights, flag_effN);
        }
        else{
            res_temp=compute_Ixy_alg1(data, sortidx, ptr_cnt, ptrVarIdx, AllLevels, n, maxbins, initbins, cut_t, 
                                      r_t, c2terms, cterms, looklog, lookH, cplx, sample_weights, flag_effN);
        }
        I_xy_u=res_temp[0];
        Ik_xy_u=res_temp[1];
        free(res_temp);


        //########################
        //I(z,x|u)
        ptrVarIdx_t[0] = ptrVarIdx[0]; //X
        ptrVarIdx_t[1] = ptrVarIdx[nbrUi+2]; //Z
        for(ll=0; ll<nbrUi; ll++) ptrVarIdx_t[ll+2]=ll+2;
        // Reset cut
        for(l=0;l<(nbrUi+2);l++){
            if(ptr_cnt[ptrVarIdx_t[l]]==1){
                for(j=0;j<initbins-1;j++) {
                    cut_t[l][j]=j*lbin+lbin-1;
                }
                cut_t[l][initbins-1]=n-1;
                for(int j=initbins;j<maxbins;j++) {
                    cut_t[l][j]=0;
                }
                r_t[l]=initbins;
            }else{
                r_t[l]=AllLevels[ptrVarIdx_t[l]];
            }
        }
        // Do opt run on I(X;Z|U)
        if(nbrUi>0){
            res_temp=compute_Ixy_cond_u_new_alg1(data, sortidx, ptr_cnt, ptrVarIdx_t, AllLevels, nbrUi, n, maxbins,
                                                initbins, cut_t, r_t, c2terms, cterms, lbin, looklog, lookH, cplx, 
                                                sample_weights, flag_effN);
        }
        else{
            res_temp=compute_Ixy_alg1(data, sortidx, ptr_cnt, ptrVarIdx_t, AllLevels, n, maxbins, initbins, cut_t, 
                                      r_t, c2terms, cterms, looklog, lookH, cplx, sample_weights, flag_effN);
        }
        I_xz_u=res_temp[0];
        Ik_xz_u=res_temp[1];
        free(res_temp);


        //########################
        //I(z,y|u)
        ptrVarIdx_t[0] = ptrVarIdx[1]; //Y
        ptrVarIdx_t[1] = ptrVarIdx[nbrUi+2]; //Z
        for(ll=0; ll<nbrUi; ll++) ptrVarIdx_t[ll+2]=ll+2;
        // Reset cut
        for(l=0;l<(nbrUi+2);l++){
            if(ptr_cnt[ptrVarIdx_t[l]]==1){
                for(j=0;j<initbins-1;j++) {
                    cut_t[l][j]=j*lbin+lbin-1;
                }
                cut_t[l][initbins-1]=n-1;
                for(int j=initbins;j<maxbins;j++) {
                    cut_t[l][j]=0;
                }
                r_t[l]=initbins;
            }else{
                r_t[l]=AllLevels[ptrVarIdx_t[l]];
            }
        }
        // Do opt run on I(Y;Z|U)
        if(nbrUi>0){
            res_temp=compute_Ixy_cond_u_new_alg1(data, sortidx, ptr_cnt, ptrVarIdx_t, AllLevels, nbrUi, n, maxbins, 
                                                 initbins, cut_t, r_t, c2terms, cterms, lbin, looklog, lookH, cplx, 
                                                 sample_weights, flag_effN);
        }
        else{
            res_temp=compute_Ixy_alg1(data, sortidx, ptr_cnt, ptrVarIdx_t, AllLevels, n, maxbins, initbins, cut_t, 
                                      r_t, c2terms, cterms, looklog, lookH, cplx, sample_weights, flag_effN);
        }
        I_yz_u=res_temp[0];
        Ik_yz_u=res_temp[1];
        free(res_temp);
    }

    //compute conditional three point mutual information
     I_xyz_u  = I_xy_u  - I_xy_zu;
     Ik_xyz_u = Ik_xy_u - Ik_xy_zu;

    ///////////////////////////////////////

    //compute Rscore

    //compute probability of 
    //not v-structure: nv 
    //not dpi inequality : dpi
            
    // nv
    nv=n*Ik_xyz_u;
    //if(I_xyz_u == 0){
    //    nv = n*-0.00001;
    //}

    xz=n*(Ik_xz_u-Ik_xy_u);
    yz=n*(Ik_yz_u-Ik_xy_u);

    if(xz < yz){
        first = xz;
        second = yz;
    } else {
        first = yz;
        second = xz;
    }
    
    // dpi
    dpi = first - log1p(exp(first-second));

    if(dpi<nv){
        //Pdpi>Pnv => Rscore=Pdpi
        Rscore=dpi;
    }
    else{
        //Pdpi<Pnv => Rscore=Pnv
        Rscore=nv;
    }

    res[0]=Rscore;
    res[1]=n*I_xyz_u;
    res[2]=n*I_xyz_u-nv;


    //free memory

    for(l=0; l<(nbrUi+3); l++)
        free(cut[l]);
    for(l=0; l<(nbrUi+2); l++)
        free(cut_t[l]);
    free(cut);
    free(cut_t);


    free(r);
    free(ptrVarIdx_t);

    for(l=0; l<(nbrUi+3); l++) free(datafactors[l]);
    free(datafactors);
    free(ptr_u2_t);
    for(l=0; l<4; l++) free(factors4_t[l]);
    free(factors4_t);
    free(r4_t);
    free(r_temp);
    free(r_t);

    return res;
}