#include <cmath>
#include <iostream>
#include <algorithm>
#include <limits>
#include <cfloat>
#include <vector>
#include <Rcpp.h>

#include <utilities.h>

using namespace Rcpp;
using namespace std;


/*================================================================================================*
 *                   Copyright (C) 2011 Paul Mineiro                                              *
 * https://github.com/etheory/fastapprox/blob/master/fastapprox                                   */
static inline float __attribute__((always_inline)) fastlog2 (float x)
{
  union { float f; uint32_t i; } vx = { x };
  union { uint32_t i; float f; } mx = { (vx.i & 0x007FFFFF) | 0x3f000000 };
  float y = vx.i;
  y *= 1.1920928955078125e-7f;

  return y - 124.22551499f
           - 1.498030302f * mx.f 
           - 1.72587999f / (0.3520887068f + mx.f);
}
static inline float __attribute__((always_inline)) fastlog (float x)
{
  return 0.69314718f * fastlog2 (x);
}
/*================================================================================================*/



double* compute_B(int K, int e, double* candidate_cut_points, double* myDist, int n, double epsilon,
                  double xmin, double* B_previous_K, double** sc_look, int* count_for_e, int E){

  double* res = new double[2];

  int ne = count_for_e[e-1];
  double B;

  if(K==1){
    B = -ne * (log(epsilon * ne / ((candidate_cut_points[e-1] - (xmin - epsilon/2))*n)));
    res[0] = 0;
    res[1] = B;
    return(res);
  } 
  else {
    double R_ne_K = compute_parametric_complexity(ne, K, sc_look);
    double min = DBL_MAX;
    int minidx;
    int nep;
    // Loop on ep
    for(int ep=K-1; ep<e-1; ep++){
      nep = count_for_e[ep-1];
      double R_nep_prevK = compute_parametric_complexity(nep, K-1, sc_look);
      double B = B_previous_K[ep-1] - (ne-nep) *
              fastlog((epsilon*(ne - nep))/((candidate_cut_points[e-1] - candidate_cut_points[ep-1])*n)) +
              fastlog((R_ne_K / R_nep_prevK) * ((E-K+2)/(K-1)));

      if(B < min){
        min = B;
        minidx = ep;
      }
    }
    int best_ep = minidx;
    double B_current_K = min;
    res[0] = best_ep;
    res[1] = B_current_K;
    return(res);
  }
} 


extern "C" SEXP mydiscretizeMDL(SEXP RmyDist, SEXP RmaxBins){

  std::vector<double> myDistVec = Rcpp::as< vector <double> >(RmyDist);
  int maxBins = Rcpp::as<int> (RmaxBins);
  double epsilon = 0.001;
  int n = myDistVec.size();
  sort(myDistVec.begin(), myDistVec.end());

  double* myDist = &myDistVec[0];

  double* middle_points = new double[n-1];
  int n_unique_points(0);
  double middle_point;
  
  for(int i=1; i<n-1; i++){
    if (myDist[i] != myDist[i+1]){
      middle_points[n_unique_points] = (myDist[i]+myDist[i+1])/2;
      n_unique_points++;  
    }
  }
  int n_candidate_cut_points = n_unique_points+1;
  double* candidate_cut_points = new double[n_candidate_cut_points];
  
  for(int i=0; i<n_unique_points; i++){
    candidate_cut_points[i] = middle_points[i];
  }
  
  int* count_for_e = new int[n_candidate_cut_points];
  for(int e=0; e<n_unique_points; e++){
    int count=0;
    for (int i=0; i<n; i++){
      if(myDist[i]<candidate_cut_points[e]){
        count++;
      } else break;
    }
    count_for_e[e] = count;
  }
  
  double xmin = myDist[0] - epsilon;
  double xmax = myDist[n-1] + epsilon;

  double** sc_look = new double*[n];
  for(int i = 0; i < n; i++){
    sc_look[i] = new double[maxBins];
    for(int k = 0; k < maxBins; k++){
      sc_look[i][k]=0.;
    }
  }


  double** dyntable = new double*[n_candidate_cut_points];
  int** dyntable_trace = new int*[n_candidate_cut_points];
  for(int i = 0; i < n_candidate_cut_points; i++){
    dyntable[i] = new double[maxBins];
    dyntable_trace[i] = new int[maxBins];
    for(int k = 0; k < maxBins; k++){
      dyntable[i][k]=0.;
      dyntable_trace[i][k]=0;
    }
  }


  // Loop on K
  for(int K=1; K<=maxBins; K++){
    double* previous_column = new double[n_candidate_cut_points];
    for(int i=0; i<n_candidate_cut_points; i++){ // Retrieve previous column
      if(K==1) previous_column[i]=0;
      else previous_column[i] = dyntable[i][K-2];
    }
    // Loop on e
    for(int e=K; e<n_candidate_cut_points; e++){
      double* res = compute_B(K, e, candidate_cut_points, myDist, n, epsilon, xmin, 
                              previous_column, sc_look, count_for_e, n_candidate_cut_points);
      dyntable[e-1][K-1] = res[1];
      dyntable_trace[e-1][K-1] = res[0];
    }
  }
  int best_K = 0;
  double best_score = DBL_MAX;
  for(int i=0; i<maxBins; i++){
    if(dyntable[n_candidate_cut_points-2][i] < best_score){
      best_score = dyntable[n_candidate_cut_points-2][i];
      best_K = i+1;
    }
  }
  double* cut_points = new double[best_K+2];
  int K = best_K;
  int curr_idx = dyntable_trace[n_candidate_cut_points-2][K-1];
  cut_points[K] = candidate_cut_points[curr_idx-1];
  int prev_idx = curr_idx;
  K--;
  for(K; K>0; K--){
    curr_idx = dyntable_trace[prev_idx][K];
    cut_points[K] = candidate_cut_points[curr_idx-1];
    prev_idx = curr_idx;
  }
  cut_points[0] = xmin;
  cut_points[best_K+1] = xmax;

  //cout << "cut points : [" ;
  //for(int i=0; i<best_K+2; i++){
  //    cout << cut_points[i] << "  ";
  //}
  //cout << "]\n";

  vector<double> values(cut_points, cut_points + best_K+2);

  // structure the output
  List result = List::create(
    _["cutpoints"] = values
  ) ;

  return result;
}
