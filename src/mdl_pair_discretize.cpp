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
#include <Info_cnt.h>
#include <Rcpp.h>


#include "utilities.h"
#include "modules_MI.h"
#include "computeInfo.h"


#define STEPMAX 50
#define EPS 1e-5
#define INIT_EQUAL_WIDTH false
#define ALPHA_EFF_LVLS 1
#define COEFF_COMB 0

using namespace Rcpp;
using namespace std;

// module to compute the conditional mutual information for mixed continuous and discrete variables

//////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////
//compute I(x,y)

//dynamic programming for optimizing variables binning

//optimize on x I(x,y): Hx - Hxy - kmdl
//optimize on y I(x,y): Hy - Hxy - kmdl
//until convergence


int** compute_Ixy_alg1(int** data, int** sortidx, int* ptr_cnt, int* ptrVarIdx,  int* AllLevels,
                       int n, int maxbins, int **cut, int *r, double* c2terms, int initbins,
                       double* looklog, double** looklbc, double* lookH, double** sc_look, int cplx,
                       double* sample_weights, bool flag_effN)
{

    int j,l,ll;

    /////////////////////////////////////

    double* res_temp=(double *)calloc(2,sizeof(double));//results res[0]->I,res[1]->I-k
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



    // Find the best initial conditions with the same number of bins (equalfreq) on all continuous variables.
    double max_res=res_temp[1];
    int max_initbins=initbins;
    for(int new_initbins=2; (new_initbins<initbins) && (initbins<20); new_initbins++){

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

        jointfactors_u(datafactors, ptr, n, 2, r, xy_factors, &rxy);
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
    int* rt1=(int *)calloc(1,sizeof(int));
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
                                      maxbins, looklog, lookH, sc_look, cplx, sample_weights, flag_effN); // 2 factors
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
                                      maxbins, looklog, lookH, sc_look, cplx, sample_weights, flag_effN); // 2 factors
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

        // Save cut points
        for(j=0; j<maxbins; j++){
            iterative_cuts[stop-1][j] = cut[0][j];
            iterative_cuts[stop-1][j+maxbins] = cut[1][j];
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

    iterative_cuts[stop][0]=-1; // mark where we stopped iterating
    iterative_cuts[stop][1]=100000*return_res[1]; // Pass Ik[X;Y]
    iterative_cuts[stop][2]=100000*return_res[0]; // Pass I[X;Y]
    iterative_cuts[stop][3]=100000*max_res; // Pass max res before optimization with equal freq
    iterative_cuts[stop][maxbins]=-1;

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


    for(l=0;l<2;l++) free(datafactors[l]);
    free(datafactors);
    free(singlefactor);

    free(factors1);

    free(rt1);

    return iterative_cuts;
}


int** compute_Ixy_cond_u_new_alg1(int** data, int** sortidx, int* ptr_cnt, int* ptrVarIdx,  int* AllLevels, int nbrUi,
                                  int n, int maxbins, int **cut, int *r, double* c2terms, int initbins, int lbin,
                                  double* looklog, double** looklbc, double* lookH, double** sc_look, int cplx,
                                  double* sample_weights, bool flag_effN)
{

    int j,l,ll;
    int STEPMAX1=50;

    /////////////////////////////////////

    double* res_temp=(double *)calloc(2,sizeof(double));//res_tempults res_temp[0]->I,res_temp[1]->I-k

    double* res=(double *)calloc(2,sizeof(double));//results res[0]->I,res[1]->I-k
    int** iterative_cuts = (int **)calloc(STEPMAX+1, sizeof(int*));
    for(int i=0; i<STEPMAX+1; i++){
        iterative_cuts[i] = (int *)calloc(maxbins*(nbrUi+2), sizeof(int));
    }

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
    for(int new_initbins=2; (new_initbins<initbins) && (initbins<20); new_initbins++){

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
            res_temp=computeMI_knml(datafactors[1],uiyxfactors[2], uiyxfactors[3],r_temp,n,c2terms, looklog, 1);
        else
            res_temp=computeMI_kmdl(datafactors[1],uiyxfactors[2], uiyxfactors[3],r_temp,n, looklog, 0);
        I_y_xu = res_temp[0]; // Before optimization on X.
        Ik_y_xu = res_temp[1];
        free(res_temp);

        r_temp[0]=r[0];
        r_temp[1]=ruiyx[1];
        r_temp[2]=ruiyx[3];
        if(cplx == 1)
            res_temp=computeMI_knml(datafactors[0],uiyxfactors[1], uiyxfactors[3],r_temp,n,c2terms, looklog, 1);
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
                                              &(r[l+2]), maxbins, looklog, lookH, sc_look, cplx, sample_weights, flag_effN); // 2 factors

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
            res_temp=computeMI_knml(datafactors[1],uiyxfactors[2], uiyxfactors[3],r_temp,n,c2terms, looklog, 1);
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
                                      &(r[0]), maxbins, looklog, lookH, sc_look, cplx, sample_weights, flag_effN); // 2 factors

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
                                              &(r[l+2]), maxbins, looklog, lookH, sc_look, cplx, sample_weights, flag_effN); // 2 factors

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
            res_temp=computeMI_knml(datafactors[0],uiyxfactors[1], uiyxfactors[3],r_temp,n,c2terms, looklog, 1);
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
                                      &(r[1]), maxbins, looklog, lookH, sc_look, cplx, sample_weights, flag_effN); // 2 factors

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
                                              &(r[l+2]), maxbins, looklog, lookH, sc_look, cplx, sample_weights, flag_effN); // 2 factors //herve

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
            res_temp=computeMI_knml(datafactors[0],uiyxfactors[0], uiyxfactors[2],r_temp,n,c2terms, looklog, 1);
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
                                              &(r[l+2]), maxbins, looklog, lookH, sc_look, cplx, sample_weights, flag_effN); //2 factors

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
            res_temp=computeMI_knml(datafactors[1],uiyxfactors[0], uiyxfactors[1],r_temp,n,c2terms, looklog, 1);
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
        cond_I  = 0.5* (I_x_yu  - I_x_u  + I_y_xu  - I_y_u);
        cond_Ik = 0.5* (Ik_x_yu - Ik_x_u + Ik_y_xu - Ik_y_u);
        //printf("0.5*(%.2f - %.2f + %.2f - %.2f) = %.2f\n", I_x_yu, I_x_u, I_y_xu, I_y_u, cond_I);
                // Save cut points
        for(j=0; j<maxbins; j++){
            for(l=0; l<(nbrUi+2); l++){
                iterative_cuts[stop1-1][j+l*maxbins] = cut[l][j];
            }
        }

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


    for(l=0; l<(nbrUi+2); l++){
        iterative_cuts[stop1][l*maxbins]=-1; // mark where we stopped iterating
        iterative_cuts[stop1][l*maxbins+1]=100000*return_res[1]; // pass Ik[X;Y|U]
        iterative_cuts[stop1][l*maxbins+2]=100000*return_res[0]; // pass I[X;Y|U]
    }

    // free memory
    free(r_temp);

    for(l=0;l<(nbrUi+2);l++) free(datafactors[l]);
    free(datafactors);
    for(l=0;l<4;l++) free(uiyxfactors[l]);
    free(uiyxfactors);
    free(factors1);
    free(rt1);

    free(ruiyx);
    free(MI);
    free(MIk);
    free(MI1);
    free(MIk1);

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
                           int nbrUi, int n, int maxbins, double* c2terms, int init_bin, double* looklog, double** looklbc,
                           double* lookH, double** sc_look, int cplx, double* sample_weights, bool flag_effN)
{

    //int** iterative_cuts = (int **)calloc(STEPMAX, sizeof(int*));
    //for(int i=0; i<STEPMAX; i++){
    //    iterative_cuts[i] = (int *)calloc(maxbins*2, sizeof(int));
    //}

    int j,k,l,ll;

    int *r=(int *)calloc((nbrUi+2),sizeof(int));

    int **cut;
    cut=(int **)calloc(nbrUi+2,sizeof(int*));
    for(l=0;l<(nbrUi+2);l++){
        cut[l]=(int *)calloc(maxbins,sizeof(int));
    }

    int lbin=floor(n/init_bin);
    if(lbin<1) {
        lbin=1;
        init_bin=n;
    }

    double* wbin=(double *)calloc(nbrUi+2, sizeof(double)); //unit width of bin for both variables
    for(j=0; j<(nbrUi+2); j++){
        if(ptr_cnt[j]==1) wbin[j] = (dataDouble[sortidx[j][n-1]][j] - dataDouble[sortidx[j][0]][j]) / init_bin;
    }


    // if(isNumberOfLevelsLessTwo(data, nbrUi, n, 2)){
    //     res[0]= 0;
    //     res[1]= 0;
    //     res[2]= n;
    //     return res;
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
                    for(j=0;j<init_bin-1;j++) {
                        cut[l][j]=j*lbin+lbin-1;
                    }
                }
                if(INIT_EQUAL_WIDTH){
                    int dataDoubleit=0;
                    double cut_value;
                    double data_value;
                    for(j=0;j<init_bin-1;j++) {
                        cut_value = dataDouble[sortidx[l][0]][l] + j*wbin[l]+wbin[l];
                        data_value = dataDouble[sortidx[l][dataDoubleit]][l];
                        while((data_value < cut_value) && (dataDoubleit < n-1)){
                            data_value = dataDouble[sortidx[l][dataDoubleit]][l];
                            dataDoubleit++;
                        }
                        cut[l][j] = dataDoubleit-1;
                    }
                }
                cut[l][init_bin-1]=n-1;
                r[l]=init_bin;
            }else{
                r[l]=AllLevels[ptrVarIdx[l]];
            }
        }


        #if _MY_DEBUG_MI
                printf("    compute_Ixy_alg1.. \n ");
                printf("    cnt: %d %d\n",ptr_cnt[ptrVarIdx[0]],ptr_cnt[ptrVarIdx[1]]);
                fflush(stdout);
        #endif
        int **iterative_cuts = compute_Ixy_alg1(data, sortidx, ptr_cnt, ptrVarIdx, AllLevels, n, maxbins, cut, r,
                                                c2terms, init_bin, looklog, looklbc, lookH, sc_look, cplx,
                                                sample_weights, flag_effN);


        // free memory

        #if _MY_DEBUG_MI
                printf("    free cut nbr=%d \n ",nbrUi);
                fflush(stdout);
        #endif

        #if _MY_DEBUG_MI
                printf("    return \n ");
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
                    for(j=0;j<init_bin-1;j++) {
                        cut[l][j]=j*lbin+lbin-1;
                    }
                }
                if(INIT_EQUAL_WIDTH){
                    int dataDoubleit=0;
                    double cut_value;
                    double data_value;
                    for(j=0;j<init_bin-1;j++) {
                        cut_value = dataDouble[sortidx[l][0]][l] + j*wbin[l]+wbin[l];
                        data_value = dataDouble[sortidx[l][dataDoubleit]][l];
                        while((data_value < cut_value) && (dataDoubleit < n-1)){
                            data_value = dataDouble[sortidx[l][dataDoubleit]][l];
                            dataDoubleit++;
                        }
                        cut[l][j] = dataDoubleit-1;
                    }
                }
                cut[l][init_bin-1]=n-1;
                r[l]=init_bin;
            }else{
                r[l]=AllLevels[ptrVarIdx[l]];
            }
        }


        int **iterative_cuts = compute_Ixy_cond_u_new_alg1(data, sortidx, ptr_cnt, ptrVarIdx, AllLevels, nbrUi, n, maxbins,
                                                           cut, r, c2terms, init_bin, lbin,looklog, looklbc, lookH, sc_look,
                                                           cplx, sample_weights, flag_effN);
        //int **iterative_cuts = old_compute_mi_cond_alg1(data, sortidx, AllLevels, ptr_cnt, ptrVarIdx, nbrUi, n,
        //                                                maxbins, c2terms, init_bin, looklog, looklbc, lookH, cplx);


        // free memory

        #if _MY_DEBUG_MI
                printf("    free cut nbr=%d \n ",nbrUi);
                fflush(stdout);
        #endif

        #if _MY_DEBUG_MI
                printf("    return \n ");
                fflush(stdout);
        #endif

        free(r);
        for(l=0;l<(nbrUi+2);l++){
            free(cut[l]);
        }
        free(cut);
        free(wbin);


        return iterative_cuts;
    }

}




//####################################################################################################
//# Dealing with input variables

void transformToFactorsContinuous(double** data, int** dataNumeric, int i, int n){

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
    //     cout << "(" << it->first << "," << it->second << ")\n";
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
void transformToFactors(double** data, int** dataNumeric, int n, int i){
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

void transformToFactorsContinuousIdx(int** dataNumeric, int** dataNumericIdx, int n, int i){

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

extern "C" SEXP mydiscretizeMutual(SEXP RmyDist1, SEXP RmyDist2, SEXP RflatU, SEXP RnbrU, SEXP RmaxBins,
                                   SEXP Rinitbin, SEXP Rcplx, SEXP Rcnt_vec, SEXP Rnlevels, SEXP ReffN){

    std::vector<double> myDist1Vec = Rcpp::as< vector <double> >(RmyDist1);
    std::vector<double> myDist2Vec = Rcpp::as< vector <double> >(RmyDist2);
    std::vector<double> cnt_vec = Rcpp::as< vector <double> >(Rcnt_vec);
    std::vector<double> nlevels = Rcpp::as< vector <double> >(Rnlevels);
    int maxbins = Rcpp::as<int> (RmaxBins);
    int init_bin = Rcpp::as<int> (Rinitbin);
    int cplx = Rcpp::as<int> (Rcplx);
    int nbrU = Rcpp::as<int> (RnbrU);
    int n = myDist1Vec.size();
    int effN = Rcpp::as<int> (ReffN);

    double* myDist1 = &myDist1Vec[0];
    double* myDist2 = &myDist2Vec[0];

    std::vector<double> vectorflatU = Rcpp::as< vector <double> > (RflatU);
    double** matrixU = new double*[nbrU];
    for(int l=0; l<nbrU; l++){
        matrixU[l] = new double[n];
    }
    if(nbrU > 0){
        for(int l=0; l<nbrU; l++){
            matrixU[l] = new double[n];
            copy(vectorflatU.begin()+(l*n), vectorflatU.begin()+((l+1)*n), matrixU[l]);
        }
    }


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

    int** dataNumeric = new int*[n];
    for(j = 0; j < n; j++){
        dataNumeric[j] = new int[nbrU+2];
    }
    for(i = 0; i < (nbrU+2); i++){
        if(cnt_vec[i]==1) transformToFactorsContinuous(data, dataNumeric, i, n); //update environment.dataNumeric not taking into account repetition
        else transformToFactors(data, dataNumeric, n, i);
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

    for(i = 0; i < (nbrU+2); i++){
        if(cnt_vec[i]==1){
            transformToFactorsContinuousIdx(dataNumeric, dataNumericIdx, n, i);
            transformToFactors(data, dataNumeric, n, i);//update environment.dataNumeric taking into account repetition
        }
    }

    // AllLevels #####################################################################################
    int* AllLevels = new int[nbrU+2];
    for(i = 0; i < (nbrU+2); i++){
        AllLevels[i] = nlevels[i];
        //AllLevels[i] = n;
    }

    // ptr_cnt #######################################################################################
    int* ptr_cnt = new int[nbrU+2];
    for(i = 0; i < (nbrU+2); i++){
        ptr_cnt[i] = cnt_vec[i];
        //ptr_cnt[i] = 1;
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
      lookH[i] = i*looklog[i];//-(i+1)*looklog[(i+1)];
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

    int** dataNumeric_red ;//progressive data rank with repetition for same values
    int** dataNumericIdx_red ;//index of sorted data
    int* AllLevels_red ;//number of levels
    int* cnt_red ;//bool continuous or not
    int* posArray_red ;//node references

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
        if(    samplesToEvaluate[i] == 1){
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

      // update with the effective number of levels
      if(cnt_red[j] == 1)
        AllLevels_red[j]=nnr;
    }

    //check effective number of levels variables continuous enough to be treat as continuous?
    if(samplesNotNA != n){
      //updateNumberofLevelsAndContinuousVariables(dataNumeric_red, AllLevels_red, cnt_red, myNbrUi, samplesNotNA);
    }

            // cout << "lev:" << cnt_red[0] << " " << cnt_red[1] << endl;


  for(int i=0; i<n; i++){
      delete [] dataNumeric[i];
  }
  delete[] dataNumeric;

  for(int i=0; i<(nbrU+2); i++){
      delete [] dataNumericIdx[i];
  }
  delete[] dataNumericIdx;

  delete[] samplesToEvaluateTamplate;
  delete[] samplesToEvaluate;

  delete[] AllLevels;
  delete[] ptr_cnt;
  delete[] posArray;



  // Declare the lookup table of the parametric complexity
  int ncol = max(maxbins, init_bin)+1;
  for(int j=0; j<(nbrU+2); j++){
      if(cnt_red[j]==0) ncol = max(ncol, (AllLevels_red[j]+1));
  }
  ncol = 1000; // combinations of Us can exceed ncol
  double** sc_look = new double*[ncol];
  for(int K = 0; K < (ncol); K++){
    sc_look[K] = new double[n+1];
    for(int i = 0; i < (n+1); i++){
          if(K==1) sc_look[K][i] = 0;
          else if (i==0) sc_look[K][i] = 0;
          else sc_look[K][i] = -1;
    }
  }
  for(int i=0; i<n; i++){
      double d = computeLogC(i,  2, looklog,  sc_look); // Initialize the c2 terms
  }

  double* sample_weights;
  if(effN != n){
      sample_weights = new double[n];
      for(int i=0; i<n; i++) sample_weights[i] = double(effN)/n;
  }
  else{
      sample_weights = NULL;
  }

  int** iterative_cuts = compute_mi_cond_alg1(dataNumeric_red, data, dataNumericIdx_red, AllLevels_red,
                                              cnt_red, posArray_red, nbrUi, n, maxbins, c2terms, init_bin,
                                              looklog, looklbc, lookH, sc_look, cplx, sample_weights, effN != n);

  //cout << "==================================" << endl;
  //cout << computeLogHDC(10, 9, 8, looklog, sc_look) << endl;
  //cout << computeLogHDC(1, 9, 4, looklog, sc_look) << endl;

  int niterations=0;
  double* res = new double[2];
  double max_res_ef;
  int** iterative_cutpoints = new int*[STEPMAX*maxbins];
  for(int i = 0; i < (STEPMAX*maxbins); i++){
      iterative_cutpoints[i] = new int[nbrUi+2];
  }
  for(int l=0; l<STEPMAX+1; l++){
    if(iterative_cuts[l][0]==-1){
        niterations=l;
        res[1] = iterative_cuts[l][1]/100000.0;
        res[0] = iterative_cuts[l][2]/100000.0;
        max_res_ef = iterative_cuts[l][3]/100000.0;
        break;
    }
    for(int k=0; k<(nbrUi+2); k++){
        i=0;
        while(iterative_cuts[l][i+maxbins*k] < iterative_cuts[l][i+maxbins*k+1]){
            iterative_cutpoints[maxbins*l+i][k] = iterative_cuts[l][i+maxbins*k];
            i++;
        }
        for(int j=i; j<maxbins; j++){
            iterative_cutpoints[maxbins*l+j][k] = -1;
        }
    }
  }

  NumericMatrix cutpoints(niterations*maxbins, nbrUi+2);
  for(i=0; i<cutpoints.nrow(); i++){
      for(j=0; j<(nbrUi+2); j++){
          cutpoints[i+j*cutpoints.nrow()] = iterative_cutpoints[i][j];//
      }
  }
  // structure the output
  List result = List::create(
    _["cutpointsmatrix"] = cutpoints,
    _["info"] = res[0],
    _["infok"] = res[1],
    _["efinfo"] = max_res_ef
  );

  for(int i=0; i<(STEPMAX+1); i++){
      free(iterative_cuts[i]);
  }
  free(iterative_cuts);

  for(int i=0; i<ncol; i++){
      delete[] sc_look[i];
  }
  delete[] sc_look;
  delete[] c2terms;

  for(int i=0; i<n; i++){
      delete [] data[i];
  }
  delete[] data;

  delete[] looklog;
  delete[] lookH;

  for(int l=0; l<nbrU; l++){
      delete[] matrixU[l];
  }
  delete[] matrixU;

  for(int i=0; i<n+1; i++){
      delete[] looklbc[i];
  }
  delete[] looklbc;
  delete[] sample_weights;

  for(int i=0; i<STEPMAX*maxbins; i++){
      delete[] iterative_cutpoints[i];
  }
  delete[] iterative_cutpoints;

   delete [] AllLevels_red;
   delete [] cnt_red;
   delete [] posArray_red;

   for(int j = 0; j<(nbrU+2); j++){
       delete [] dataNumericIdx_red[j];
       delete [] dataNumeric_red[j];
   }

   delete [] dataNumericIdx_red;
   delete [] dataNumeric_red;

  return result;
}
