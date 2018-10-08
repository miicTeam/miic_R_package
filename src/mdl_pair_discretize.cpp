
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


#define STEPMAX 50
#define EPS 1e-5
#define INIT_EQUAL_WIDTH false
#define ALPHA_EFF_LVLS 1
#define COEFF_COMB 0

using namespace Rcpp;
using namespace std;

// module to compute the conditional mutual information for mixed continuous and discrete variables

//////////////////////////////////////////////////////////////////////////////////////////////

void reset_u_cutpoints(int** cut, int nbrUi, int* ptr_cnt, int* ptrVarIdx, int init_bin, int maxbins,
                       int lbin, int* r, int* AllLevels, int n){
    for(int l=2;l<(nbrUi+2);l++){
        if(ptr_cnt[ptrVarIdx[l]]==1){
            for(int j=0;j<init_bin-1;j++) {
                cut[l][j]=j*lbin+lbin-1;
            }
            cut[l][init_bin-1]=n-1;
            for(int j=init_bin;j<maxbins;j++) {
                cut[l][j]=0;
            }
            r[l]=init_bin;
        }else{
            r[l]=AllLevels[ptrVarIdx[l]];
        }
    }
}

void reset_cutpoints(int** cut, int nbrUi, int* ptr_cnt, int* ptrVarIdx, int init_bin, int maxbins,
                     int lbin, int* r, int* AllLevels, int n){
    for(int l=0;l<(nbrUi+2);l++){
        if(ptr_cnt[ptrVarIdx[l]]==1){
            for(int j=0;j<init_bin-1;j++) {
                cut[l][j]=j*lbin+lbin-1;
            }
            cut[l][init_bin-1]=n-1;
            for(int j=init_bin;j<maxbins;j++) {
                cut[l][j]=0;
            }
            r[l]=init_bin;
        }else{
            r[l]=AllLevels[ptrVarIdx[l]];
        }
    }
}

inline __attribute__((always_inline))
double logdbico(int n,    int k, double** looklbc, double* looklog){


    // string str=ToString(n) + "-" + ToString(k);
    // tuple<int,int> tpl (n,k);
    // map<tuple<int,int>,double>::iterator it=looklbc.find(tpl);
    if( looklbc[n][k] != 0){
        return looklbc[n][k];
    }
    else{
        double v=logchoose(n,k, looklog);
        looklbc[n][k]=v;
        return v;
    }
}

//#define _MY_DEBUG_MInoU 1
//#define _MY_DEBUG_MI 1
//#define _MY_DEBUG_NEW_OPTFUN 1


///////////////////////////////////////////////////////////////////////////////////////////

inline __attribute__((always_inline))
double old_optfun_onerun_kmdl_coarse(int *sortidx_var, int *data, int nbrV, int **factors, int *r, int *optfun, double sc,
                                      int n, int nnr,int *cut, int *r_opt, int maxbins, double* looklog, double** looklbc, double* lookH)
{

        // cout << "infunctionlevels: " << nnr << endl;

        int i,j,k,m;

        int coarse=ceil(1.0*nnr/maxbins);//step coarse graining
        if (coarse<1) coarse=1;
        int np=ceil(1.0*nnr/coarse); //number of possible cuts

        //temp variables to memorize optimization cuts
        int *memory_cuts_idx=(int *)calloc(np,sizeof(int));//indexes of the cuts (1..np)
        int *memory_cuts_pos=(int *)calloc(np,sizeof(int));//positions of the cuts (1..n)
        int c,c2;

        //variables for the computation of the complexity
        double k_sc=0;
        double sc2,scr,sctemp;


        //dynamic programming optimize function and memorize of cuts
        double fmax;//I-kmdl
        int* nc =(int *)calloc(np,sizeof(int));
        int nctemp;//Hx-Hxy-LogCx

        //entropy in kj interval for the <nbrV>+1 terms
        double* H_kj=(double *)calloc(nbrV+1,sizeof(double));

        //function max value at each step
        double* I=(double *)calloc(np,sizeof(double));
        I[0]=0;
        double* I_0k =(double *)calloc(np,sizeof(double));
        double I_kj,t;

        //number of points in the intervals
        int nxj,nx;

        int xyu;

        //
        int ir;
        int nk,nj,nkj;
        int njforward,nkforward;


        //dynamic programming steps terms in the 0k interval

        //H_0k needs initialization at zero (problem: double H_0k[nbrV+1][np];)
        double **H_0k=(double **)calloc(nbrV+1,sizeof(double*));// x y u xu yu xyu
        for(m=0;m<nbrV+1;m++){
            H_0k[m]=(double *)calloc(np,sizeof(double));
        }

        int **nxyu=(int**)calloc(nbrV,sizeof(int*));// x y u xu yu xyu
        for(m=0;m<nbrV;m++){
            nxyu[m]=(int*)calloc(r[m],sizeof(int));
        }

        int **nxyu_k=(int**)calloc(nbrV,sizeof(int*));// x y u xu yu xyu
        for(m=0;m<nbrV;m++){
             nxyu_k[m]=(int*)calloc(r[m],sizeof(int));
             // cout << "\tniv: " << r[m] << "\t";
        }

        bool *check_repet = (bool*)calloc(n, sizeof(bool));
        for(m=0;m<(n-1);m++){
            check_repet[m] = (data[sortidx_var[m+1]]!=data[sortidx_var[m]]);
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

        //computing statistics of the <nbrV> terms and the entropy

        njforward=0;//iterator on values
        ir=0;//iterator on not repeated values
        while(ir<coarse){

            for(m=0;m<nbrV;m++){

                if(optfun[m] != 0){ //compute only necessary terms

                    xyu=factors[m][sortidx_var[njforward]];

                    nxyu[m][xyu]++;

                    if(nxyu[m][xyu] != 1){
                        H_0k[m][0] -= nxyu[m][xyu]*looklog[nxyu[m][xyu]] - (nxyu[m][xyu]-1)*looklog[nxyu[m][xyu]-1];
                    }
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

        if(optfun[nbrV] != 0) H_0k[nbrV][0]=-nj*looklog[nj];

        for(m=0;m<nbrV+1;m++) if(optfun[m] != 0)    I[0]+=optfun[m]*H_0k[m][0];

        //initizialitation
        memory_cuts_idx[0]=0;
        I_0k[0]=I[0];
        nc[0]=1;

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
            for(m=0;m<nbrV+1;m++) if(optfun[m] != 0) H_0k[m][j]=H_0k[m][j-1];

            //computing statistics of the <nbrV> terms and the entropy

            // njforward iterator on values
            ir=0;//iterator on not repeated values
            while((ir<coarse )&& (njforward<n)){

                for(m=0;m<nbrV;m++){

                    if(optfun[m] != 0){ //compute only necessary terms

                        xyu=factors[m][sortidx_var[njforward]];

                        nxyu[m][xyu]++;

                        if(nxyu[m][xyu] != 1) H_0k[m][j]-=nxyu[m][xyu]*looklog[nxyu[m][xyu]]-(nxyu[m][xyu]-1)*looklog[nxyu[m][xyu]-1];

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

            if(optfun[nbrV] != 0) H_0k[nbrV][j]=-nj*looklog[nj];

            I_0k[j]=0;
            for(m=0;m<nbrV+1;m++) if(optfun[m] != 0)    I_0k[j]+=optfun[m]*H_0k[m][j];

            //complexity terms (local version)

            k_sc=sc*looklog[(nj+1)];

            //computing complexity term for the [0 j] interval
            // for the 2 bins combinatorial term : [0 k][k+1 j] mode

            // string str=ToString(nj) + "-" + ToString(1);
            // map<string,double>::iterator it=looklbc.find(str);

            // tuple<int,int> tpl (nj,1);
            // map<tuple<int,int>,double>::iterator it=looklbc.find(tpl);



            // if( looklbc[tpl][nj][1] != 0){
            //     sc2=2*k_sc+ looklbc[tpl][nj][1];
            // }
            // else{
            //     double v=log(dbico(nj,1));
            //     looklbc[tpl][nj][1]=v;
            //     sc2=2*k_sc+ v;
            // }

            // formule 7 part, here combinat
            //cout << "k_sc : " << k_sc << "\t";
            //cout << "dbico(" << nj << ", 1) : " << logdbico(nj, 1, looklbc) << "\t";
            sc2=2*k_sc+ COEFF_COMB*logdbico(nj,1,looklbc, looklog);
            //cout << "sc2 : " << sc2 << "\n";



            // if(looklbc[nxj][1]==-1) looklbc[nxj][1]=log(dbico(nxj,1));
            // sc2=2*k_sc+looklbc[nxj][1];
            // sc2=2*k_sc+log(dbico(nxj,1));

            //init

            fmax=-DBL_MAX;

            for(m=0;m<nbrV+1;m++) if(optfun[m] != 0) H_kj[m]=H_0k[m][j];

            #if _MY_DEBUG_NEW_OPTFUN
                printf("  ");
                printf("\n  [0 -%d] = %lf\n",j,I_0k[j]);fflush(stdout);
            #endif

            for(m=0;m<nbrV;m++) {
                if(optfun[m] != 0) {
                    for(xyu=0;xyu<r[m];xyu++)
                        nxyu_k[m][xyu]=nxyu[m][xyu];
                }
            }

            //moving k

            //k iterator on possible cuts
            nkforward=0;//iterator on values
            //it iterator on not repeated values

            for(k=0;k<=j-1;k++){//k=1...n-2 possible cuts

                        ir=0;//iterator on not repeated values
                        while( ir<coarse ){

                                for(m=0;m<nbrV;m++){

                                    if(optfun[m] != 0){ //compute only necessary terms

                                        xyu=factors[m][sortidx_var[nkforward]];
                                        nxyu_k[m][xyu]--;

                                        H_kj[m]-= lookH[nxyu_k[m][xyu]];
                                    }

                                }

                                ir += int(check_repet[nkforward]);
                                nkforward++;
                        }

                        //nx=nxj-(k+1)*coarse;
                        nkj=njforward-nkforward;
                        if(optfun[nbrV] != 0) H_kj[nbrV]=-nkj*looklog[nkj];

                        //if(optfun[nbrV] != 0) H_kj[nbrV]=-lookH[nkj];

                        #if _MY_DEBUG_NEW_OPTFUN
                            printf("\nnxj=%d  \n   ",nkj);
                            for(m=0;m<nbrV;m++) {
                                for(xyu=0;xyu<r[m];xyu++) printf("nxyu_k[%d][%d]=%d ",m,xyu,nxyu_k[m][xyu]);fflush(stdout);
                            }
                            printf("\n");
                        #endif

                        I_kj=0;
                        for(m=0;m<nbrV+1;m++) if(optfun[m] != 0)    I_kj+=optfun[m]*H_kj[m];

                        //computing complexity term for the [0..k][k+1 j] mode

                        //recursive number of bins combinatorial term

                        // string str=ToString(nj) + "-" + ToString(nc[k]);
                        // map<string,double>::iterator it=looklbc.find(str);

                        // tuple<int,int> tpl (nj,nc[k]);
                        // map<tuple<int,int>,double>::iterator it=looklbc.find(tpl);


                        // if( looklbc[tpl][nj][nc[k]] != 0){
                        //     scr=(nc[k]+1)*k_sc+ looklbc[tpl][nj][nc[k]];
                        // }
                        // else{
                        //     double v=log(dbico(nj,nc[k]));
                        //     looklbc[tpl][nj][nc[k]]=v;
                        //     scr=(nc[k]+1)*k_sc+ v;
                        // }

                        //cout << "optfun\tn : " << nj << "\tk : " << nc[k] << "\tidx : " << k << "\n";
                        scr=(nc[k]+1)*k_sc+ COEFF_COMB*logdbico(nj,nc[k],looklbc, looklog);

                        // if(looklbc[nxj][nc[k]]==-1) looklbc[nxj][nc[k]]=log(dbico(nxj,nc[k]));
                        // scr=(nc[k]+1)*k_sc+looklbc[nxj][nc[k]];

                        // scr=(nc[k]+1)*k_sc+log(dbico(nxj,(nc[k])));


                        #if _MY_DEBUG_NEW_OPTFUN


                            printf("  ");
                            for(m=0;m<nbrV+1;m++) printf("H_kj[%d]=%lf ",m,H_kj[m]);fflush(stdout);
                            printf("\n   [%d -%d][%d - %d] = %lf (%lf+%lf-%lf)\n",0,k,k+1,j,I_0k[k]+I_kj-sc2,I_0k[k],I_kj,sc2);
                            printf("   [0 -?-%d][%d - %d] = %lf (%lf+%lf-%lf)\n",k,k+1,j,I[k]+I_kj-scr,I[k],I_kj,scr);

                        #endif

                        //which mode is more convinient?
                        // [0...k][k+1 j] mode or [0 k][k+1 j] mode?
                        // in other terms more cuts or only one cut?

                        nk=nkforward-1;//position of the actual possible cut

                        if( I_0k[k] - sc2 > I[k] - scr ){ //one cut in k ore more?
                                t=I_0k[k] + I_kj-sc2;
                                if (fmax<t){
                                    c=-k-1;// convention to refers to the two bins [0 k] [k+1 j]
                                    c2=nk;
                                    fmax=t;
                                    sctemp=sc2;
                                    nctemp=2;
                                }
                        }else{//more cuts
                            t=I[k] + I_kj-scr;//[0.. cuts.. k-1][k j]
                            if (fmax<t){
                                c=k+1;
                                c2=nk;
                                fmax=t;
                                sctemp=scr;
                                nctemp=nc[k]+1;
                            }
                        }

                        #if _MY_DEBUG_NEW
                            printf("   f=%lf\n",fmax);fflush(stdout);
                        #endif
            }

            // optimized function for the interval [0 j]
            I[j]=fmax+sctemp;
            // number of optimal cuts
            nc[j]=nctemp;
            //cout << "fmax : " << fmax << "\tnctemp : " << nctemp << "\tnc[j] : " << nc[j] << "\tt : " << t <<"\n" ;
            // index  of the (last) optimal cut
            memory_cuts_idx[j]=c;
            // position  of the (last) optimal cut
            memory_cuts_pos[j]=c2;

            #if _MY_DEBUG_NEW_OPTFUN
                printf("\n>>>j=%d: fmax=%lf cut[%d]=%d ncut=%d\n",j,fmax,j,memory_cuts_idx[j],nc[j]);fflush(stdout);

                printf("\n-----------\n");

            #endif
        }

    // free memory
    free(nc);
    free(I);
    free(I_0k);
    free(H_kj);
    free(check_repet);

    for(m=0;m<nbrV+1;m++) free(H_0k[m]);
    free(H_0k);

    for(m=0;m<nbrV;m++) free(nxyu_k[m]);
    free(nxyu_k);

    for(m=0;m<nbrV;m++) free(nxyu[m]);
    free(nxyu);

    // reconstruction of the optimal cuts from the memory cuts indexes and positions
    *r_opt=reconstruction_cut_coarse(memory_cuts_idx,memory_cuts_pos, np, n, cut);

    free(memory_cuts_idx);
    free(memory_cuts_pos);
    return I[np-1]/n;

}


int** old_compute_mi_cond_alg1(int** data, int** sortidx,  int* AllLevels, int* ptr_cnt, int* ptrVarIdx, int nbrUi,
                                 int n, int maxbins, double* c2terms, int initbins,
                                 double* looklog, double** looklbc, double* lookH, int cplx)
{


    int** iterative_cuts = (int **)calloc(STEPMAX+1, sizeof(int*));
    for(int i=0; i<STEPMAX+1; i++){
        iterative_cuts[i] = (int *)calloc(maxbins*(nbrUi+2), sizeof(int));
    }

    double* res_temp=(double *)calloc(2,sizeof(double));//results res[0]->I,res[1]->I-k
    double* res=(double *)calloc(3,sizeof(double));//results res[0]->I,res[1]->I-k

    int j,k,l,ll;
    int np;    // int np=ceil(1.0*n/coarse);

    int *r=(int *)calloc((nbrUi+2),sizeof(int));

    int **cut;
    cut=(int **)calloc(nbrUi+2,sizeof(int*));
    for(l=0;l<(nbrUi+2);l++){
        cut[l]=(int *)calloc(maxbins,sizeof(int));
    }
    int init_bin=initbins;

    int lbin=floor(n/init_bin);
    if(lbin<1) {
        lbin=1;
        init_bin=n;
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

        cout << "OLD NBRU==0 : NOT IMPLEMENTED !\n";
        return iterative_cuts;
    }

    //////////////////////////////////////////////////////////////////



    double sc;

    double sc_comb;//term complexity due to the optimization on many models
    //double  sc_comb_u;
    double Hu;


     double MInew;

     // double logn=log(1.0*n);

    double *MI = new double[STEPMAX];
    double *MIk = new double[STEPMAX];
    double I_av,Ik_av;


    int stop,i,flag;

    // allocation

    // int *memory_cut=(int *)calloc(np,sizeof(int));

    ////////////////////////////////////
    //initialitize cuts vectors


    for(l=0;l<(nbrUi+2);l++){
            if(ptr_cnt[ptrVarIdx[l]]==1){
                for(j=0;j<init_bin-1;j++) {
                    cut[l][j]=j*lbin+lbin-1;
                }
                cut[l][init_bin-1]=n-1;
                r[l]=init_bin;
            }else{
                r[l]=AllLevels[ptrVarIdx[l]];
            }
    }


    //initialization of datafactors
    //note: data are sorted based on sortidx[ptrVarIdx[l]] vectors (l=0,1,..nbrUi+2)


    int **datafactors;
    datafactors=(int **)calloc((nbrUi+2),sizeof(int*));

    for(l=0;l<(nbrUi+2);l++){
        datafactors[l]=(int *)calloc(n,sizeof(int));
    }


    int uu;
    int sjj;

    for(l=0;l<(nbrUi+2);l++){ //compute datafactors based on the positions of cut points in vector <cut>
        uu=0;
        if(ptr_cnt[ptrVarIdx[l]]==1){

            update_datafactors(sortidx, ptrVarIdx[l], datafactors, l, n, cut);

            // for(j=0;j<=n-1;j++){
            //     sjj=sortidx[ptrVarIdx[l]][j];
            //     if(j>cut[l][uu]) uu++;
            //     datafactors[l][sjj]=uu;
            // }
        }else{
            for(j=0;j<=n-1;j++){//discrete case
                datafactors[l][j]=data[ptrVarIdx[l]][j];
            }
        }
    }


    //

    int **uiyxfactors;//({ui},{uiy}),{uix},{uiyx})

    uiyxfactors=(int **)calloc(4,sizeof(int*));
    for(l=0;l<4;l++){
        uiyxfactors[l]=(int *)calloc(n,sizeof(int));
    }
    int *ruiyx;
    ruiyx=(int *)calloc(4,sizeof(int));

    //

    int flag_d=0;
    for(l=0;l<(nbrUi+2);l++)
        flag_d += ptr_cnt[ptrVarIdx[l]];

    if(flag_d == 0){//all discrete

        //compute joint factors u xu yu xyu

        jointfactors_uiyx(datafactors,-1,  n, nbrUi, r, uiyxfactors, ruiyx);

        //compute MI and knml complexity

        if(cplx == 1)
            res_temp=computeMIcond_knml(uiyxfactors,ruiyx,r,n,c2terms, looklog);
        else
            res_temp=computeMIcond_kmdl(uiyxfactors,ruiyx,r,n, looklog);

        res[0]=n;
        res[1]=res_temp[0];
        res[2]=res_temp[0]-res_temp[1];


        // free memory
        for(l=0;l<(nbrUi+2);l++)
            free(cut[l]);
        free(cut);

        for(l=0;l<4;l++) free(uiyxfactors[l]);
        free(uiyxfactors);

        // free(memory_cut);

        for(l=0;l<(nbrUi+2);l++) free(datafactors[l]);
        free(datafactors);

        free(res_temp);
        free(r);
        free(ruiyx);

        return iterative_cuts;

    }


    // #if _MY_DEBUG_MI
    //             printf("datafactors 0:\n ");
    //             for(j=0;j<n;j++) printf("%d\n",datafactors[0][j]);
    //             printf("\n");fflush(stdout);
    // #endif


    ///////////////////////////////////

    int **factors1=(int **)calloc(1,sizeof(int*));
    int *rt1=(int *)calloc(1,sizeof(int));
    int *old_optfun1=(int *)calloc(2,sizeof(int));

    int ru=0;

    //allocation factors


    int **factors3=(int **)calloc(3,sizeof(int*));
    int *rt3=(int *)calloc(3,sizeof(int));
    int *old_optfun3=(int *)calloc(4,sizeof(int));

    //

    int *ufactors;
    ufactors=(int *)calloc(n,sizeof(int));


    //

    //int **uxfactors;
    //uxfactors=(int **)calloc(2,sizeof(int*));

    //uxfactors[0]=(int *)calloc(n,sizeof(int));
    //uxfactors[1]=(int *)calloc(n,sizeof(int));

    //int *rux=(int *)calloc(2,sizeof(int));


    //

    int *ptrUiIdx=(int *)calloc(nbrUi,sizeof(int));
    for(l=0;l<nbrUi;l++) ptrUiIdx[l]=l+2;

    int *ptr_u1_t=(int *)calloc(nbrUi+1,sizeof(int));


    // #if _MY_DEBUG_MI
    //             printf("allocation done:\n ");
    //             for(ll=0;ll<nbrUi+2;ll++) printf("r[%d]=%d ",ll,r[ll]);
    //             printf("\n");fflush(stdout);
    // #endif


    /////////////////////////////////////////////////////////////////////////
    //compute I(x,y | u)

    //compute joint factors u

    jointfactors_u(datafactors,ptrUiIdx,n, nbrUi, r, ufactors, &ru,&Hu);


    ////////////////////////////////////////

    if(ptr_cnt[ptrVarIdx[0]]==1){//continuous

        //optimize on x
        //I(x;u)

        //init variables for the optimization run

        sc=0.5*(ru-1);

        factors1[0]=ufactors;//u
        old_optfun1[0]=-1;//u
        old_optfun1[1]=1;//x
        rt1[0]=ru;

        //optimization run on ptrVarIdx[0]


        MInew=old_optfun_onerun_kmdl_coarse(sortidx[ptrVarIdx[0]], data[ptrVarIdx[0]], 1, factors1, rt1, old_optfun1, sc,n,
                                        AllLevels[ptrVarIdx[0]], cut[0], &(r[0]), maxbins, looklog, looklbc, lookH);

        //update datafactors of ptrVarIdx[0]

        update_datafactors(sortidx,    ptrVarIdx[0], datafactors, 0, n, cut);

        #if _MY_DEBUG_MI
                printf("opt x: MInew=%lf ",MInew);
                for(ll=0;ll<nbrUi+2;ll++) printf("r[%d]=%d ",ll,r[ll]);
                printf("\n");fflush(stdout);
        #endif

    }

    ////////////////////////////////////

    if(ptr_cnt[ptrVarIdx[1]]==1){//continuous

        //optimize on y
        //I(y;u)


        //init variables for the optimization run

        sc=0.5*(ru-1);

        factors1[0]=ufactors;//u
        old_optfun1[0]=-1;//u
        old_optfun1[1]=1;//y
        rt1[0]=ru;

        //optimization run on ptrVarIdx[1]

        MInew=old_optfun_onerun_kmdl_coarse(sortidx[ptrVarIdx[1]], data[ptrVarIdx[1]],    1, factors1, rt1, old_optfun1, sc,n,
                                        AllLevels[ptrVarIdx[1]], cut[1], &(r[1]),maxbins, looklog, looklbc, lookH);

        //update datafactors of ptrVarIdx[1]

        update_datafactors(sortidx, ptrVarIdx[1], datafactors, 1, n, cut);

        #if _MY_DEBUG_MI
                printf("opt y: MInew=%lf ",MInew);
                for(ll=0;ll<nbrUi+2;ll++) printf("r[%d]=%d ",ll,r[ll]);
                printf("\n");fflush(stdout);
        #endif

    }

    ////////////////////////////////////

    //optimization function:
    //I(x;u)+I(y;u)

    //optimization over variables in u
    //f=2Hu-Hxu-Hyu

    int nb_ui_cnt=0;

    for(l=0;l<nbrUi;l++){

        if(ptr_cnt[ptrVarIdx[l+2]]==1){//continuous

            //compute joint factors u xu yu xyu

            jointfactors_uiyx(datafactors,l+2,  n, nbrUi, r, uiyxfactors, ruiyx);

            //init variables for the optimization run

            sc=0.5*(r[0]-1)*(ruiyx[0])+0.5*(r[1]-1)*(ruiyx[0]);

            factors3[0]=uiyxfactors[0];//u
            factors3[1]=uiyxfactors[1];//yu
            factors3[2]=uiyxfactors[2];//xu

            old_optfun3[0]=2;//u
            old_optfun3[1]=-1;//yu
            old_optfun3[2]=-1;//xu
            old_optfun3[3]=0;//ui

            rt3[0]=ruiyx[0];//u
            rt3[1]=ruiyx[1];//yu
            rt3[2]=ruiyx[2];//xu

            //optimization run on ptrVarIdx[l+2]

            MInew=old_optfun_onerun_kmdl_coarse(sortidx[ptrVarIdx[l+2]],  data[ptrVarIdx[l+2]],    3, factors3, rt3, old_optfun3, sc, n,
                                            AllLevels[ptrVarIdx[l+2]], cut[l+2], &(r[l+2]), maxbins, looklog, looklbc, lookH);

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

    jointfactors_u(datafactors,ptrUiIdx,n, nbrUi, r, ufactors, &ru,&Hu);

    //compute complexity term from the u set
    if(n > maxbins) np=maxbins;
    else np=n;
    // np=n;


    // if(nb_ui_cnt!=0 && np>ru-1){

    //         // string str=ToString(np) + "-" + ToString(ru-1);
    //         // map<string,double>::iterator it=looklbc.find(str);

    //         tuple<int,int> tpl (np,ru-1);
    //         map<tuple<int,int>,double>::iterator it=looklbc.find(tpl);


    //         if( it != looklbc.end()){
    //             sc_comb_u = it->second;
    //         }
    //         else{
    //             double v=log(dbico(np,ru-1));
    //             looklbc[tpl]=v;
    //             sc_comb_u = v;
    //         }


    //     // if(looklbc[np][ru-1]==-1) looklbc[np][ru-1]=log(dbico(np,ru-1));
    //     // sc_comb_u=looklbc[np][ru-1];
    // }
    // else sc_comb_u=0;

    /////////////////////////////////////////

    for(ll=0;ll<nbrUi;ll++) ptr_u1_t[ll+1]=ll+2;

     MI[0]=0;
     MIk[0]=-DBL_MAX;

    flag=0;

    //ITERATION OF OPTIMIZATION PROCESS OVER THE VARIABLES (maximum iterations <STEPMAX>)

    for(stop=1;stop<STEPMAX;stop++)
    {

        // sc_comb=sc_comb_u;
        sc_comb=0;

    //////////////////

        if(ptr_cnt[ptrVarIdx[0]]==1){//continuous

            //optimize on x
            //I(x;yu)

            ptr_u1_t[0]=1;

            //compute joint factors u

            jointfactors_u(datafactors,ptr_u1_t,n, nbrUi+1, r, ufactors, &ru,&Hu);

            //init variables for the optimization run

            sc=0.5*(ru-1);

            factors1[0]=ufactors;//uy
            old_optfun1[0]=-1;//uy
            old_optfun1[1]=1;//x
            rt1[0]=ru;//uy

            //optimization run on ptrVarIdx[0]

            MInew=old_optfun_onerun_kmdl_coarse(sortidx[ptrVarIdx[0]], data[ptrVarIdx[0]], 1, factors1, rt1, old_optfun1, sc,n,
                                            AllLevels[ptrVarIdx[0]], cut[0], &(r[0]), maxbins, looklog, looklbc, lookH);

            //update datafactors of ptrVarIdx[0]

            update_datafactors(sortidx, ptrVarIdx[0], datafactors, 0, n, cut);

            // complexity term

            if(AllLevels[ptrVarIdx[0]] > maxbins) np=maxbins;
            else np=AllLevels[ptrVarIdx[0]];
            // np=n;

            //

            // string str=ToString(np) + "-" + ToString(r[0]-1);
            // map<string,double>::iterator it=looklbc.find(str);

            // tuple<int,int> tpl (np,r[0]-1);
            // map<tuple<int,int>,double>::iterator it=looklbc.find(tpl);

            // if( looklbc[tpl][np][r[0]-1]){
            //     sc_comb += looklbc[tpl][np][r[0]-1];
            // }
            // else{
            //     double v=log(dbico(np,r[0]-1));
            //     looklbc[tpl][np][r[0]-1]=v;
            //     sc_comb += v;
            // }
            sc_comb += COEFF_COMB*logdbico(np,r[0]-1, looklbc, looklog);

            // if(looklbc[np][r[0]-1]==-1) looklbc[np][r[0]-1]=log(dbico(np,r[0]-1));
            // sc_comb+=looklbc[np][r[0]-1];


            #if _MY_DEBUG_MI
                    printf("%d\n",stop);fflush(stdout);
                    printf("MInew=%lf ",MInew);
                    for(ll=0;ll<nbrUi+2;ll++) printf("r[%d]=%d ",ll,r[ll]);
                    printf("\n");fflush(stdout);
            #endif

        }
    //////////////////////////////////////////////

        if(ptr_cnt[ptrVarIdx[1]]==1){//continous


            //optimize on y
            //I(y;xu)

            ptr_u1_t[0]=0;

            //compute joint factors u

            jointfactors_u(datafactors,ptr_u1_t,n, nbrUi+1, r, ufactors, &ru,&Hu);

            //init variables for the optimization run

            sc=0.5*(ru-1);

            factors1[0]=ufactors;//ux
            old_optfun1[0]=-1;//ux
            old_optfun1[1]=1;//x
            rt1[0]=ru;//ux

            //optimization run on ptrVarIdx[1]

            MInew=old_optfun_onerun_kmdl_coarse(sortidx[ptrVarIdx[1]], data[ptrVarIdx[1]], 1, factors1, rt1, old_optfun1, sc,n,
                                            AllLevels[ptrVarIdx[1]], cut[1], &(r[1]), maxbins, looklog, looklbc, lookH);

            //update datafactors of ptrVarIdx[1]

            update_datafactors(sortidx, ptrVarIdx[1], datafactors, 1, n, cut);

            #if _MY_DEBUG_MI
                    printf("MInew=%lf ",MInew);
                    for(ll=0;ll<nbrUi+2;ll++) printf("r[%d]=%d ",ll,r[ll]);
                    printf("\n");fflush(stdout);
            #endif

            // complexity term
            if(AllLevels[ptrVarIdx[1]] > maxbins) np=maxbins;
            else np=AllLevels[ptrVarIdx[1]];
            // np=n;

            //np=ceil(1.0*AllLevels[ptrVarIdx[1]]/coarse);

            //

            // string str=ToString(np) + "-" + ToString(r[1]-1);
            // map<string,double>::iterator it=looklbc.find(str);

            // tuple<int,int> tpl (np,r[1]-1);
            // map<tuple<int,int>,double>::iterator it=looklbc.find(tpl);

            // if( looklbc[tpl][np][r[1]-1] != 0){
            //     sc_comb += looklbc[tpl][np][r[1]-1];
            // }
            // else{
            //     double v=log(dbico(np,r[1]-1));
            //     looklbc[tpl][np][r[0]-1]=v;
            //     sc_comb += v;
            // }

            sc_comb += COEFF_COMB*logdbico(np,r[1]-1, looklbc, looklog);


            // if(looklbc[np][r[1]-1]==-1) looklbc[np][r[1]-1]=log(dbico(np,r[1]-1));
            // sc_comb+=looklbc[np][r[1]-1];


        }

        ///////////////////////////////////////////

        //compute joint factors u xu yu xyu

        jointfactors_uiyx(datafactors,-1,  n, nbrUi, r, uiyxfactors, ruiyx);

        // Save cut points
        for(j=0; j<maxbins; j++){
            for(l=0; l<(nbrUi+2); l++){
                iterative_cuts[stop-1][j+l*maxbins] = cut[l][j];
            }
        }

        if(cplx == 1)
            res_temp=computeMIcond_knml(uiyxfactors,ruiyx,r,n,c2terms, looklog);
        else
            res_temp=computeMIcond_kmdl(uiyxfactors,ruiyx,r,n, looklog);

        //combinatorial term
        // res_temp[1]-=(log(dbico(n,r[0]-1))+log(dbico(n,r[1]-1)))/n;//-log(dbico(np,ru-1)))/n;
        // res_temp[1]-=(log(dbico(n,r[0]-1))+log(dbico(n,r[1]-1))-log(dbico(n,ruiyx[0]-1)))/n;
        // res_temp[1]-=(log(dbico(n,r[0]-1))+log(dbico(n,r[1]-1)))/n;
        res_temp[1]-=sc_comb/n;
        // res_temp[1]-=log(n)/n;

        #if _MY_DEBUG_MI
            printf("\nI=%lf Ik=%lf\n",res_temp[0],res_temp[1]);
            printf("\n");fflush(stdout);
        #endif

        //iteration reach fixed point or periodic cycle?

        for(i=stop-1;i>0;i--){
                if(fabs(res_temp[1]-MIk[i]) < EPS ){
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
            if(flag) break;//stop iteration
            MIk[stop]=res_temp[1];
            MI[stop]=res_temp[0];
    }//for

    if(flag){
        res[1]=I_av;
        res[2]=I_av-Ik_av;
    }
    else{
        res[1]=res_temp[0];
        res[2]=res_temp[0]-res_temp[1];
    }
    res[0]=n;

    if(!flag) stop--; //If stopped because of loop condition and not flag, stop1 must be decremented.
    for(l=0; l<(nbrUi+2); l++){
        iterative_cuts[stop][l*maxbins]=-1; // mark where we stopped iterating
        iterative_cuts[stop][l*maxbins+1]=100000*(res[1]-res[2]); // pass Ik[X;Y|U]
        iterative_cuts[stop][l*maxbins+2]=100000*res[1]; // pass I[X;Y|U]
    }




    #if _MY_DEBUG_MI
        printf("\nI=%lf Ik=%lf ",res_temp[0],res_temp[1]);
        for(ll=0;ll<nbrUi+2;ll++) printf("r[%d]=%d ",ll,r[ll]);
            printf("\n");fflush(stdout);
    #endif

    // free memory


    for(l=0;l<(nbrUi+2);l++)
        free(cut[l]);
    free(cut);

    free(ufactors);
    //free(uxfactors[0]);
    //free(uxfactors[1]);
    //free(uxfactors);


    // free(memory_cut);

    for(l=0;l<(nbrUi+2);l++) free(datafactors[l]);
    free(datafactors);

    free(factors1);
    free(factors3);
    free(res_temp);
    free(r);
    free(ruiyx);
    free(rt1);
    free(old_optfun1);
    free(rt3);
    free(old_optfun3);
    //free(rux);
    free(ptrUiIdx);
    free(ptr_u1_t);

    for(l=0;l<4;l++) free(uiyxfactors[l]);
    free(uiyxfactors);

    delete [] MI;
    delete [] MIk;

    return iterative_cuts;
}

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
//            optimize f on x : nbrV=2, terms to take into account Hx-Hxyu :
//                optfun[0]=0 //factors x= 0
//                optfun[6]=1 => optimization on x
//                optfun[4]= -1 => factors[4]    // factors xyu => factors yu
//                opfun=[0,0,0,0,-1,1]
//            optimize f on y : 2 terms to take into account Hyu-Hxyu :
//                optfun[0]=0 //factors x= 0
//                optfun[6]=0 => no terms on y
//                optfun[4]=1 => factors[4] //factors yu=> factors u
//                optfun[5]=-1 => factors[5] //factors xyu=> factors xu
//                opfun=[0,0,0,1,-1,0]
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
double* optfun_onerun_kmdl_coarse(int *sortidx_var, int *data, int nbrV, int **factors, int *r,
                                  double sc, int sc_levels1, int sc_levels2, int n, int nnr,int *cut,
                                  int *r_opt, int maxbins, double* looklog, double** looklbc, double* lookH,
                                  double** cterms, int cplx, int flag_allow_unique_bin=1) {

    int i,j,k,m;
    double pxy = 1.0;

    int coarse=ceil(1.0*nnr/maxbins);//step coarse graining
    //if (coarse<cbrt(n)) coarse=cbrt(n);
    if (coarse<1) coarse=1;
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

    double* weights=(double *)calloc(4,sizeof(double));

    //function max value at each step
    double* Ik=(double *)calloc(np,sizeof(double)); // The optimal information value found at each idx.
    double* Ik_0k =(double *)calloc(np,sizeof(double)); // The information value for a unique bin fom idx 0 to k.
    double Ik_kj; // The information value for a bin fom idx k to j.
    double t;

    //number of points in the intervals
    int nxj,nx;

    int xyu;

    int ir; // Iterator on non repeated values
    int nk,nj,nkj; // Number of values in intervals [0,k], [0,j], [k,j]
    int njforward,nkforward; // Indexes at current position of j and k
    int n_without_u;


    //entropy in kj interval for the <nbrV>+1 terms
    double* Hk_kj=(double *)calloc(nbrV,sizeof(double));

    //dynamic programming steps terms in the 0k interval
    double **Hk_0k=(double **)calloc(nbrV,sizeof(double*));// x y u xu yu xyu //herve
    for(m=0;m<nbrV;m++){
        Hk_0k[m]=(double *)calloc(np,sizeof(double));
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
        weights[0]=-1;
        weights[1]=1;
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
                Hk_0k[m][0] += weights[m]*lookH[nxyu[m][xyu]-1];
            }

            if(m == 1){ //herve
                if(cplx==0 && nxyu[m][xyu]==1) Hk_0k[m][0]  -= sc * looklog[n];  //MDL
                else if(cplx==1){
					Hk_0k[m][0] -= pxy * ( computeLogC(nxyu[m][xyu]  , sc_levels1, cterms) -
                                           computeLogC(nxyu[m][xyu]-1, sc_levels1, cterms) );//NML
                }
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
            for(xyu=0;xyu<r[m];xyu++){
                printf("nxyu[%d][%d]=",m,xyu);fflush(stdout);
                printf("%d  ",nxyu[m][xyu]);fflush(stdout);
            }
        }
        printf("\n");
    #endif

    nj=njforward;

    for(m=0;m<nbrV;m++) Ik_0k[0] += Hk_0k[m][0];//herve

    //initizialitation
    memory_cuts_idx[0]=0;
    Ik[0]=Ik_0k[0];
    nc[0]=1;


    #if _MY_DEBUG_NEW_OPTFUN
            printf("  ");
            for(m=0;m<nbrV;m++) printf("H_0k[%d]=%lf ",m,H_0k[m][0]);fflush(stdout);
            printf("\n  [0 -0] = %lf\n",I_0k[0]);fflush(stdout);
    #endif

    //moving j over the np possible cuts

    for(j=1;j<=np-1;j++){ //j=1...n-1

        #if _MY_DEBUG_NEW_OPTFUN
                printf("j=%d\n",j);fflush(stdout);
        #endif

        //COMPUTING STATISTICS AND FUNCTION FOR THE INTEVAL [0 j] -> I_0k

        //initialization with previous interval [0 j-1]
        for(m=0;m<nbrV;m++) Hk_0k[m][j]=Hk_0k[m][j-1]; //herve

        //computing statistics of the <nbrV> terms and the entropy

        // njforward iterator on values
        ir=0;//iterator on not repeated values
        while((ir<coarse )&& (njforward<n)){

            for(m=0;m<nbrV;m++){

                xyu=factors[m][sortidx_var[njforward]];
                nxyu[m][xyu]++;

                if(nxyu[m][xyu] != 1){
                    Hk_0k[m][j] += weights[m]*lookH[nxyu[m][xyu]-1];
                }

                if(m == 1){ //herve
                    if(cplx==0 && nxyu[m][xyu]==1) Hk_0k[m][j]  -= sc * looklog[n];  //MDL
                    else if(cplx==1){
						Hk_0k[m][j] -= pxy * ( computeLogC(nxyu[m][xyu]  , sc_levels1, cterms) -
                                               computeLogC(nxyu[m][xyu]-1, sc_levels1, cterms) );//NML
                    }
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
                for(xyu=0;xyu<r[m];xyu++) printf("nxyu[%d][%d]=",m,xyu);fflush(stdout);
                for(xyu=0;xyu<r[m];xyu++) printf("%d  ",nxyu[m][xyu]);fflush(stdout);
            }
            printf("\n");
        #endif

        Ik_0k[j]=0;
        for(m=0;m<nbrV;m++) Ik_0k[j] += Hk_0k[m][j]; //herve
        Ik[j] = Ik_0k[j];

        Imax=-DBL_MAX;

        for(m=0;m<nbrV;m++) Hk_kj[m]=Hk_0k[m][j];

        #if _MY_DEBUG_NEW_OPTFUN
            printf("  ");
            printf("\n  [0 -%d] = %lf\n",j,I_0k[j]);fflush(stdout);
        #endif

        for(m=0;m<nbrV;m++) {
            for(xyu=0;xyu<r[m];xyu++) nxyu_k[m][xyu]=nxyu[m][xyu];
        }

        //moving k

        //k iterator on possible cuts
        nkforward=0;//iterator on values

        for(k=0;k<=j-1;k++){//k=1...n-2 possible cuts

            ir=0;//iterator on not repeated values
            while( ir<coarse ){

            	n_without_u=njforward-nkforward-1;
                for(m=0;m<nbrV;m++){

                    xyu=factors[m][sortidx_var[nkforward]];
                    nxyu_k[m][xyu]--;

                    Hk_kj[m] -= weights[m]*lookH[nxyu_k[m][xyu]];

                    if(m == 1){ //herve
                        if(cplx==0 && nxyu_k[m][xyu]==0) Hk_kj[m]  += sc * looklog[n];  //MDL
                        else if(cplx==1){
							Hk_kj[m] -= pxy * ( computeLogC(nxyu_k[m][xyu]  , sc_levels1, cterms) -
                                                computeLogC(nxyu_k[m][xyu]+1, sc_levels1, cterms) );//NML
                        }
                    }
                }

                ir += int(check_repet[nkforward]);
                nkforward++;
            }

            nkj=njforward-nkforward;

            Ik_kj=0;
            for(m=0;m<nbrV;m++) Ik_kj += Hk_kj[m];
			if(cplx==1){
                Ik_kj -= log((1.0*np-1)/(sc_levels2-1) - 1) + 1;
            }

            nk=nkforward-1;//position of the actual possible cut

            if( (Ik[k] + Ik_kj) > Ik[j] ) {
                t=Ik[k] + Ik_kj; //[0.. cuts.. k-1][k j] //herve
                if (Imax<t){
                    Imax=t;
                    nctemp=nc[k]+1;
                    Ik[j]=Imax; // optimized function for the interval [0 j]
                    nc[j]=nctemp; // number of optimal cuts
                    memory_cuts_idx[j] = k + 1; // index  of the (last) optimal cut
                    memory_cuts_pos[j] = nk; // position  of the (last) optimal cut
                }
            }
        }


    }

    // free memory
    double* return_res = new double[2];
    return_res[0] = 0;
    return_res[1] = Ik[np-1]/n;
    free(nc);
    free(Ik);
    free(Ik_0k);
    free(Hk_kj);
    free(weights);
    free(check_repet);

    for(m=0;m<nbrV;m++) free(Hk_0k[m]);
    free(Hk_0k);

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


//////////////////////////////////////////////////////////////////////////////////////////////
//compute I(x,y)

//dynamic programming for optimizing variables binning

//optimize on x I(x,y): Hx - Hxy - kmdl
//optimize on y I(x,y): Hy - Hxy - kmdl
//until convergence


int** compute_Ixy_alg1(int** data, int** sortidx, int* ptr_cnt, int* ptrVarIdx,  int* AllLevels,
                       int n, int maxbins, int **cut, int *r, double* c2terms, int init_bin,
                       double* looklog, double** looklbc, double* lookH, double** cterms, int cplx)
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


    double sc;
    double sc_comb=0;//term complexity due to optimization on many models

     double* MInew;

    double* MI = (double *)calloc(STEPMAX,sizeof(double));
    double* MIk = (double *)calloc(STEPMAX,sizeof(double));
    double I_av,Ik_av;

    int stop,i,flag;

    int sc_levels1,sc_levels2; //herve
    //

    int **factors1=(int **)calloc(2,sizeof(int*));
    int* singlefactor = (int *)calloc(n, sizeof(int));
    int *rt1=(int *)calloc(2,sizeof(int));


    MI[0]=0;
    MIk[0]=-DBL_MAX;

    flag=0;

    int sc_levels_x = r[0]; // Number of levels of the first variable
    int sc_levels_y = r[1]; // Number of levels of the second variable
    int rx = r[0];
    int ry = r[1];
    int* r_t = (int*) calloc(2, sizeof(int));
    int np;
    for(stop=1;stop<STEPMAX;stop++)
    {

        ///////////////////////////////////////////
        if(ptr_cnt[ptrVarIdx[0]]==1){

            //optimize on x
            //I(x;y)

            factors1[0]=datafactors[1]; //y
            factors1[1]=singlefactor; //One single bin at start

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
            MInew=optfun_onerun_kmdl_coarse(sortidx[ptrVarIdx[0]], data[ptrVarIdx[0]], 2, factors1, rt1, sc,
                                            sc_levels1, sc_levels2, n, AllLevels[ptrVarIdx[0]], cut[0], &(r[0]),
                                            maxbins, looklog, looklbc, lookH, cterms, cplx);
        }

        ////////////////////////////////////////////////
        if(ptr_cnt[ptrVarIdx[1]]==1){

            //opt y
            //I(x;y)
            factors1[0]=datafactors[0];//x before its optimization
            factors1[1]=singlefactor; //One single bin at start

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
            MInew=optfun_onerun_kmdl_coarse(sortidx[ptrVarIdx[1]], data[ptrVarIdx[1]],2, factors1, rt1, sc,
                                            sc_levels1, sc_levels2, n, AllLevels[ptrVarIdx[1]], cut[1], &(r[1]),
                                            maxbins, looklog, looklbc, lookH, cterms, cplx);
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
        jointfactors_u(datafactors,ptr,n, 2, r, xy_factors, &rxy, &Hxy);
        r_temp[0]=rx;
        r_temp[1]=ry;
        r_temp[2]=rxy;
        if(cplx == 1)
            res=computeMI_knml(datafactors[0],datafactors[1],xy_factors,r_temp,n,c2terms,looklog, 0);
        else
            res=computeMI_kmdl(datafactors[0],datafactors[1],xy_factors,r_temp,n,looklog, 0);
		//Adding combinatorial term
        if(ptr_cnt[ptrVarIdx[0]] == 1 && rx>1){
            np = min(maxbins, AllLevels[ptrVarIdx[0]]);
            res[1] -= 0.5*(rx-1)*(log((1.0*np-1) / (rx-1) - 1) + 1)/n;
        }
        if(ptr_cnt[ptrVarIdx[1]] == 1 && ry>1){
            np = min(maxbins, AllLevels[ptrVarIdx[1]]);
            res[1] -= 0.5*(ry-1)*(log((1.0*np-1) / (ry-1) - 1) + 1)/n;
        }


        for(i=stop-1;i>0;i--){
            if( fabs(res[1]-MIk[i]) < EPS ){
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
        if(flag){
            break;
        }
        MIk[stop]=res[1];
        MI[stop]=res[0];

    }//for
    if(!flag) stop--; //If stopped because of loop condition and not flag, stop1 must be decremented.
    iterative_cuts[stop][0]=-1; // mark where we stopped iterating
    iterative_cuts[stop][1]=100000*res[1]; // Pass Ik[X;Y]
    iterative_cuts[stop][2]=100000*res[0]; // Pass I[X;Y]
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


    for(l=0;l<2;l++) free(datafactors[l]);
    free(datafactors);
    free(singlefactor);

    free(factors1);

    free(rt1);

    return iterative_cuts;
}


int** compute_Ixy_cond_u_new_alg1(int** data, int** sortidx, int* ptr_cnt, int* ptrVarIdx,  int* AllLevels, int nbrUi,
                                  int n, int maxbins, int **cut, int *r, double* c2terms, int init_bin, int lbin,
                                  double* looklog, double** looklbc, double* lookH, double** sc_look, int cplx)
{

    int j,l,ll;
    int STEPMAX1=50;

    /////////////////////////////////////

    double* res=(double *)calloc(2,sizeof(double));//results res[0]->I,res[1]->I-k
    int** iterative_cuts = (int **)calloc(STEPMAX+1, sizeof(int*));
    for(int i=0; i<STEPMAX+1; i++){
        iterative_cuts[i] = (int *)calloc(maxbins*(nbrUi+2), sizeof(int));
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

     double* MInew;

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
    double Ik_x_yu;
    double Ik_y_xu;
    double Ik_x_u;
    double Ik_y_u;

    // Complexity factor (# of levels) passed to each optimization
    int sc_levels1, sc_levels2;


    int **factors1=(int **)calloc(2,sizeof(int*));

    int *rt1=(int *)calloc(2,sizeof(int));
    int *r_old = (int *)calloc((nbrUi+2),sizeof(int));
    for(int ll=0; ll<(nbrUi+2); ll++){
        r_old[ll] = r[ll];
    }
    int U_counter;


    int sc_levels_x; // Number of levels of the first variable
    int sc_levels_y; // Number of levels of the second variable
    flag1=0;
    int flag_allow_unique_bin = 1;
    int max_U_counter = 3;
    int np;
    for(stop1=1;stop1<STEPMAX1;stop1++)
    {

        ///////////////////////////////////////////
        //optimize I(y;xu) over x and u
        U_counter = 0;
        while(U_counter < max_U_counter){
            for(l=0;l<nbrUi;l++){

                if(ptr_cnt[ptrVarIdx[l+2]]==1){
                    //opt u
                    //I(y;xu)
                    jointfactors_uiyx(datafactors, l+2,  n, nbrUi, r_old, uiyxfactors, ruiyx);

                    //init variables for the optimization run
                    factors1[0]=uiyxfactors[3];//xyu
                    factors1[1]=uiyxfactors[2];//xu

                    rt1[0]=ruiyx[3];//xyu
                    rt1[1]=ruiyx[2];//xu
                    sc_levels_x = r_old[0];
                    sc_levels_y = r_old[1];

                    sc = 0.5*(sc_levels_y-1);

                    sc_levels1 = sc_levels_y;
                    sc_levels2 = r_old[l+2];
                    //sc_levels2 = init_bin;

                    // Run optimization on U.
                    MInew=optfun_onerun_kmdl_coarse(sortidx[ptrVarIdx[l+2]], data[ptrVarIdx[l+2]], 2, factors1, rt1,
                                                    sc, sc_levels1, sc_levels2, n, AllLevels[ptrVarIdx[l+2]], cut[l+2],
                                                    &(r[l+2]), maxbins, looklog, looklbc, lookH, sc_look,
                                                    cplx, flag_allow_unique_bin); // 2 factors

                    //update_datafactors(sortidx, ptrVarIdx[l+2], datafactors, l+2, n, cut); //move to outside U loop
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
        r_temp[0]=r_old[1]; //Y
        r_temp[1]=ruiyx[2]; //UX
        r_temp[2]=ruiyx[3]; //UXY
        if(cplx == 1)
            res=computeMI_knml(datafactors[1],uiyxfactors[2], uiyxfactors[3],r_temp,n,c2terms, looklog, 1);
        else
            res=computeMI_kmdl(datafactors[1],uiyxfactors[2], uiyxfactors[3],r_temp,n, looklog, 1);
        I_y_xu = res[0]; // Before optimization on X.
        Ik_y_xu = res[1];
        // Adding combinatorial term
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

            rt1[0]=ruiyx[1];//uy
            rt1[1]=ruiyx[0];//u
            sc_levels_x = r_old[0];
            sc_levels_y = r_old[1];
            sc = 0.5*(sc_levels_y-1);

            #if _MY_DEBUG_MInoU
                printf("start optfun\n ");
                fflush(stdout);
            #endif

            sc_levels1 = sc_levels_y; //herve
            sc_levels2 = sc_levels_x; //herve
            //sc_levels2 = init_bin;
            // Run optimization on X.
            res=optfun_onerun_kmdl_coarse(sortidx[ptrVarIdx[0]], data[ptrVarIdx[0]], 2, factors1, rt1, sc,
                                          sc_levels1, sc_levels2, n, AllLevels[ptrVarIdx[0]], cut[0], &(r[0]),
                                          maxbins, looklog, looklbc, lookH, sc_look, cplx,
                                          flag_allow_unique_bin); // 2 factors

            //update_datafactors(sortidx, ptrVarIdx[0], datafactors, 0, n, cut); //moved to after Y opt
        }

        // Reset cutpoints on U
        reset_u_cutpoints(cut, nbrUi, ptr_cnt, ptrVarIdx, init_bin, maxbins, lbin, r, AllLevels, n);
        for(l=0; l<nbrUi; l++){
            if(ptr_cnt[ptrVarIdx[l+2]]==1) update_datafactors(sortidx, ptrVarIdx[l+2], datafactors, l+2, n, cut);
            r_old[l+2] = r[l+2];//r[l+2] is set to init_bin during reset_u_cutpoints
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

                    rt1[0]=ruiyx[3]; //xyu
                    rt1[1]=ruiyx[1]; //yu
                    sc_levels_y = r_old[1];
                    sc_levels_x = r_old[0];

                    sc_levels1 = sc_levels_x; //herve
                    sc_levels2 = r_old[l+2]; //herve
                    //sc_levels2 = init_bin;
                    sc = 0.5*(sc_levels_x-1)*ruiyx[1];
                    // Run optimization on U.
                    MInew=optfun_onerun_kmdl_coarse(sortidx[ptrVarIdx[l+2]], data[ptrVarIdx[l+2]], 2, factors1, rt1, sc,
                                                    sc_levels1, sc_levels2, n, AllLevels[ptrVarIdx[l+2]], cut[l+2], &(r[l+2]),
                                                    maxbins, looklog, looklbc, lookH, sc_look, cplx,
                                                    flag_allow_unique_bin); // 2 factors

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
        r_temp[1]=ruiyx[1];
        r_temp[2]=ruiyx[3];
        if(cplx == 1)
            res=computeMI_knml(datafactors[0],uiyxfactors[1], uiyxfactors[3],r_temp,n,c2terms, looklog, 1);
        else
            res=computeMI_kmdl(datafactors[0],uiyxfactors[1], uiyxfactors[3],r_temp,n, looklog, 1);
        I_x_yu = res[0]; //Before updating X and Y.
        Ik_x_yu = res[1];
        // Adding combinatorial term
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

            rt1[0]=ruiyx[2]; //ux
            rt1[1]=ruiyx[0]; //u
            sc_levels_y = r_old[1];
            sc_levels_x = r_old[0];
            sc = 0.5*(sc_levels_x-1);

            #if _MY_DEBUG_MInoU
                printf("start optfun\n ");
                fflush(stdout);
            #endif

            sc_levels1 = sc_levels_x; //herve
            sc_levels2 = sc_levels_y; //herve
            //sc_levels2 = init_bin;
            // Run optimization on Y.
            res=optfun_onerun_kmdl_coarse(sortidx[ptrVarIdx[1]], data[ptrVarIdx[1]], 2, factors1, rt1, sc,
                                          sc_levels1, sc_levels2, n, AllLevels[ptrVarIdx[1]], cut[1], &(r[1]),
                                          maxbins, looklog, looklbc, lookH, sc_look, cplx,
                                          flag_allow_unique_bin); // 2 factors

            //update_datafactors(sortidx, ptrVarIdx[1], datafactors, 1, n, cut); //moved to end of loop1
        }

        // Reset cutpoints on U
        reset_u_cutpoints(cut, nbrUi, ptr_cnt, ptrVarIdx, init_bin, maxbins, lbin, r, AllLevels, n);
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

                    rt1[0]=ruiyx[2];//xu
                    rt1[1]=ruiyx[0];//u

                    sc_levels1 = r_old[0];//x
                    sc_levels2 = r_old[l+2];//u
                    //sc_levels2 = init_bin;

                    sc = 0.5*(sc_levels1-1);

                    //optimization run on ptrVarIdx[l+2]
                    res=optfun_onerun_kmdl_coarse(sortidx[ptrVarIdx[l+2]], data[ptrVarIdx[l+2]], 2, factors1, rt1,
                                                  sc, sc_levels1, sc_levels2, n, AllLevels[ptrVarIdx[l+2]], cut[l+2],
                                                  &(r[l+2]), maxbins, looklog, looklbc, lookH, sc_look,
                                                  cplx, flag_allow_unique_bin);

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
            res=computeMI_knml(datafactors[0],uiyxfactors[0], uiyxfactors[2],r_temp,n,c2terms, looklog, 1);
        else
            res=computeMI_kmdl(datafactors[0],uiyxfactors[0], uiyxfactors[2],r_temp,n, looklog, 1);
        I_x_u = res[0]; //After optimization on U.
        Ik_x_u = res[1];
        // Adding combinatorial term
        for(ll=0;(ll<nbrUi); ll++){
            np = min(AllLevels[ptrVarIdx[ll+2]], maxbins);
            if((ptr_cnt[ptrVarIdx[ll+2]]==1) && (r_old[ll+2]>1))
                Ik_x_u -= (r_old[ll+2]-1)*(log((1.0*np-1) / (r_old[ll+2]-1) - 1) + 1)/n;
        }

        // Reset cutpoints on U
        reset_u_cutpoints(cut, nbrUi, ptr_cnt, ptrVarIdx, init_bin, maxbins, lbin, r, AllLevels, n);
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

                    //init variables for the optimization run
                    factors1[0]=uiyxfactors[1];//yu
                    factors1[1]=uiyxfactors[0];//u

                    rt1[0]=ruiyx[1];//yu
                    rt1[1]=ruiyx[0];//u

                    sc_levels1 = r_old[1];//y
                    sc_levels2 = r_old[l+2];//u
                    //sc_levels2 = init_bin;

                    sc = 0.5*(sc_levels1-1);

                    //optimization run on ptrVarIdx[l+2]
                    res=optfun_onerun_kmdl_coarse(sortidx[ptrVarIdx[l+2]], data[ptrVarIdx[l+2]], 2, factors1, rt1, sc,
                                                  sc_levels1, sc_levels2, n, AllLevels[ptrVarIdx[l+2]], cut[l+2], &(r[l+2]),
                                                  maxbins, looklog, looklbc, lookH, sc_look, cplx,
                                                  flag_allow_unique_bin);

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
            res=computeMI_knml(datafactors[1],uiyxfactors[0], uiyxfactors[1],r_temp,n,c2terms, looklog, 1);
        else
            res=computeMI_kmdl(datafactors[1],uiyxfactors[0], uiyxfactors[1],r_temp,n, looklog);
        I_y_u = res[0]; //After optimization on U.
        Ik_y_u = res[1];
        // Adding combinatorial term
        for(ll=0;(ll<nbrUi); ll++){
            np = min(AllLevels[ptrVarIdx[ll+2]], maxbins);
            if((ptr_cnt[ptrVarIdx[ll+2]]==1) && (r_old[ll+2]>1))
                Ik_y_u -= (r_old[ll+2]-1)*(log((1.0*np-1) / (r_old[ll+2]-1) - 1) + 1)/n;
        }

        // Reset cutpoints on U
        reset_u_cutpoints(cut, nbrUi, ptr_cnt, ptrVarIdx, init_bin, maxbins, lbin, r, AllLevels, n);
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
        res[0] = 0.5* (I_x_yu - I_x_u + I_y_xu - I_y_u);
        res[1] = 0.5* (Ik_x_yu - Ik_x_u + Ik_y_xu - Ik_y_u);
        //printf("I = 0.5*(%.2f - %.2f + %.2f - %.2f) = %.2f\n", I_x_yu, I_x_u, I_y_xu, I_y_u, res[0]);
        //printf("Ik = 0.5*(%.2f - %.2f + %.2f - %.2f) = %.2f\n", Ik_x_yu, Ik_x_u, Ik_y_xu, Ik_y_u, res[1]);

        // Save cut points
        for(j=0; j<maxbins; j++){
            for(l=0; l<(nbrUi+2); l++){
                iterative_cuts[stop1-1][j+l*maxbins] = cut[l][j];
            }
        }

        // Test stop condition on stop1
        for(i=stop1-1;i>0;i--){
            if( fabs(res[1]-MIk1[i]) < EPS ) { // If no real improvement over last information
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
        //flag1=true;
        if(flag1) break;
        MIk1[stop1]=res[1];
        MI1[stop1]=res[0];

    }//end stop1
    if(flag1){
        res[0] = I_av1;
        res[1] = Ik_av1;
    }

    if(!flag1) stop1--; //If stopped because of loop condition and not flag, stop1 must be decremented.
    for(l=0; l<(nbrUi+2); l++){
        iterative_cuts[stop1][l*maxbins]=-1; // mark where we stopped iterating
        iterative_cuts[stop1][l*maxbins+1]=100000*res[1]; // pass Ik[X;Y|U]
        iterative_cuts[stop1][l*maxbins+2]=100000*res[0]; // pass I[X;Y|U]
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
                           double* lookH, double** sc_look, int cplx)
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
                                                c2terms, init_bin, looklog, looklbc, lookH, sc_look, cplx);


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


        int **iterative_cuts = compute_Ixy_cond_u_new_alg1(data, sortidx, ptr_cnt, ptrVarIdx, AllLevels, nbrUi, n, maxbins, cut, r, c2terms,
                                                           init_bin, lbin,looklog, looklbc, lookH, sc_look, cplx);
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
                                   SEXP Rinitbin, SEXP Rcplx, SEXP Rcnt_vec, SEXP Rnlevels){

    std::vector<double> myDist1Vec = Rcpp::as< vector <double> >(RmyDist1);
    std::vector<double> myDist2Vec = Rcpp::as< vector <double> >(RmyDist2);
    std::vector<double> cnt_vec = Rcpp::as< vector <double> >(Rcnt_vec);
    std::vector<double> nlevels = Rcpp::as< vector <double> >(Rnlevels);
    int maxbins = Rcpp::as<int> (RmaxBins);
    int init_bin = Rcpp::as<int> (Rinitbin);
    int cplx = Rcpp::as<int> (Rcplx);
    int nbrU = Rcpp::as<int> (RnbrU);
    int n = myDist1Vec.size();

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
      double d = computeLogC(i, 2, sc_look); // Initialize the c2 terms
  }


  int** iterative_cuts = compute_mi_cond_alg1(dataNumeric_red, data, dataNumericIdx_red, AllLevels_red,
                                              cnt_red, posArray_red, nbrUi, n, maxbins, c2terms, init_bin,
                                              looklog, looklbc, lookH, sc_look, cplx);

  int niterations=0;
  double* res = new double[2];
  int** iterative_cutpoints = new int*[STEPMAX*maxbins];
  for(int i = 0; i < (STEPMAX*maxbins); i++){
      iterative_cutpoints[i] = new int[nbrUi+2];
  }
  for(int l=0; l<STEPMAX+1; l++){
    if(iterative_cuts[l][0]==-1){
        niterations=l;
        res[1] = iterative_cuts[l][1]/100000.0;
        res[0] = iterative_cuts[l][2]/100000.0;
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
    _["infok"] = res[1]
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
