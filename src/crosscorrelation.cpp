#include <Rcpp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>

using namespace Rcpp;
using std::string;
using std::vector;

// Transform data to factors
void transformToFactors(const vector<vector<string>>& data,
    vector<vector<double>>& data_numeric, int i) {
  // create a dictionary to store the factors of the strings
  std::map<string, int> myMap;
  // clean the dictionary since it is used column by column
  myMap.clear();
  myMap["NA"] = -1;
  myMap[""] = -1;
  int factor = 0;

  for (size_t j = 0; j < data.size(); j++) {
    auto it = myMap.find(data[j][i]);
    if (it != myMap.end()) {
      data_numeric[j][i] = it->second;
    } else {
      myMap[data[j][i]] = factor;
      data_numeric[j][i] = factor;
      factor++;
    }
  }
}

// Evaluate the effective number of samples in the dataset
int crosscorrelation(const vector<vector<double>>& data_numeric,
    vector<double>& correlation, int c_max) {
  int n_rows = data_numeric.size();
  int n_cols = data_numeric[0].size();
  vector<double> corr(c_max);
  vector<double> mean_conf(n_cols);

  int shift = 0;
  for (int k = 0; k < n_cols; k++) {
    mean_conf[k] = 0.0;
    for (int i = shift; i < n_rows; i++) {
      mean_conf[k] = mean_conf[k] + data_numeric[i][k];
    }
    mean_conf[k] = mean_conf[k] * 1.0 / (n_rows - shift);
  }

  double sumk, sumi, norm;
  // c=0 : a data_numeric is equal to itself
  for (int c = 0; c < c_max; c++) {
    sumi = 0;
    for (int i = shift; i < n_rows - c; i++) {
      sumk = 0;
      for (int k = 0; k < n_cols; k++) {
        sumk += (data_numeric[i][k] - mean_conf[k]) *
                (data_numeric[i + c][k] - mean_conf[k]);
      }
      sumk /= n_cols;
      sumi += sumk;
    }
    sumi /= (n_rows - c - shift);
    if (c == 0) norm = sumi;
    corr[c] = sumi / norm;

    correlation.push_back(corr[c]);
  }

  int n_eff = floor(0.5 + n_rows * (1 - corr[1]) / (1 + corr[1]));

  return n_eff;
}

// Evaluate the effective number of samples in the dataset
// [[Rcpp::export]]
List evaluateEffn(DataFrame input_data) {
  auto data = as<vector<vector<string>>>(input_data);
  int n_samples = data.size();
  int n_nodes = data[0].size();
  vector<vector<double>> data_numeric(n_samples, vector<double>(n_nodes));

  for (int i = 0; i < n_nodes; i++) {
    transformToFactors(data, data_numeric, i);
  }

  vector<double> correlation;
  int c_max = 30;
  int n_eff = crosscorrelation(data_numeric, correlation, c_max);

  if (n_eff > n_samples) n_eff = n_samples;

  List result = List::create(
      _["correlation"] = correlation,
      _["n_eff"]       = n_eff
  );
  return result;
}
