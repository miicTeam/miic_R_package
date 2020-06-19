#include <Rcpp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>

#define FIRSTLINE 0    // keep at 0
#define FIRSTCOLUMN 0  // 1=first column corresponds to sample name/number
using namespace std;
using namespace Rcpp;
using uint = unsigned int;

// Transform data to factors
void transformToFactors(uint numSamples,
    std::vector<std::vector<std::string> > data, double** dataNumeric, int i) {
  // create a dictionary to store the factors of the strings
  map<string, int> myMap;
  // clean the dictionary since it is used column by column
  myMap.clear();
  myMap["NA"] = -1;
  myMap[""] = -1;
  int factor = 0;

  for (uint j = 0; j < numSamples; j++) {
    map<string, int>::iterator it = myMap.find(data[j][i]);
    if (it != myMap.end()) {
      dataNumeric[j][i] = it->second;
    } else {
      myMap[data[j][i]] = factor;
      dataNumeric[j][i] = factor;
      factor++;
    }
  }
}

// Read input and fill all structures
double** reading_input(int& row_num, int col_num,
    std::vector<std::string> vectorData, vector<vector<string> >& state) {
  vector<string> vec;
  // convert input data
  for (int i = 0; i < (int)vectorData.size(); i++) {
    if (i >= col_num) {
      if (i % col_num == 0) {
        if (i != col_num) {
          state.push_back(vec);
          vec.clear();
        }
      }
      vec.push_back(vectorData[i]);
    }
  }
  state.push_back(vec);

  double** dataNumeric = new double*[row_num];
  for (int pos = 0; pos < row_num; pos++) {
    dataNumeric[pos] = new double[col_num];
  }

  for (int i = 0; i < col_num; i++) {
    transformToFactors(row_num, state, dataNumeric, i);
  }

  return dataNumeric;
}

// Evaluate the effective number of samples in the dataset
int crosscorrelation(const int row_num, const int col_num,
    double** configuration, vector<double>& correlationV, int c_max) {
  double sumk, sumi, norm;
  double* Corr = new double[c_max];
  double* meanconf = new double[col_num];

  int i, k, c, neff;
  int start = 0, shift = 0;

  for (k = FIRSTCOLUMN; k < col_num; k++) {
    meanconf[k] = 0.0;
    for (i = shift + FIRSTLINE; i < row_num; i++) {
      meanconf[k] = meanconf[k] + configuration[i][k];
    }
    meanconf[k] = meanconf[k] * 1.0 / (row_num - shift - FIRSTLINE);
  }

  // c=0 : a configuration is equal to itself
  for (c = start; c < c_max; c++) {
    sumi = 0;
    for (i = shift + FIRSTLINE; i < row_num - c; i++) {
      sumk = 0;
      for (k = FIRSTCOLUMN; k < col_num; k++) {
        sumk = sumk + (configuration[i][k] - meanconf[k]) *
                          (configuration[i + c][k] - meanconf[k]);
      }
      sumk = sumk / (col_num - FIRSTCOLUMN);
      sumi = sumi + sumk;
    }
    sumi = sumi / (row_num - c - shift - FIRSTLINE);
    if (c == start) norm = sumi;
    Corr[c] = sumi / norm;

    correlationV.push_back(Corr[c]);
  }

  neff = floor(0.5 + row_num * (1 - Corr[1]) / (1 + Corr[1]));

  delete[] Corr;
  delete[] meanconf;

  return neff;
}

// Evaluate the effective number of samples in the dataset
// [[Rcpp::export]]
List evaluateEffn(
    SEXP inputDataR, SEXP variable_numR, SEXP sample_numR) {
  vector<double> correlationV;
  std::vector<std::string> vectorData;
  std::vector<std::vector<std::string> > state;

  vectorData = Rcpp::as<vector<string> >(inputDataR);

  int variable_num = Rcpp::as<int>(variable_numR);
  int sample_num = Rcpp::as<int>(sample_numR);

  int Neff, c_max = 30;

  double** dataNumeric =
      reading_input(sample_num, variable_num, vectorData, state);
  //  Compute crosscorrelation C(c) and identify R:
  Neff = crosscorrelation(
      sample_num, variable_num, dataNumeric, correlationV, c_max);

  if (Neff > sample_num) Neff = sample_num;

  List result = List::create(
      _["correlation"] = correlationV,
      _["neff"]        = Neff
  );

  for (int pos = 0; pos < sample_num; pos++) delete[] dataNumeric[pos];
  delete[] dataNumeric;

  return result;
}
