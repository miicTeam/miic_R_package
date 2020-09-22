//*****************************************************************************
// Filename   : tmiic_lag_data.cpp                   Namespace: none                                   
//
// Author     : Franck SIMON                     Creation date: 16 dec 2020 
//
// Description: Lag input data for temporal mode of miic (tmiic)
//
// Changes history:
//*****************************************************************************

//=============================================================================
// INCLUDE SECTION
//=============================================================================
#include <vector>
#include <Rcpp.h>

using std::vector;
using std::string;

using Rcpp::_;
using Rcpp::as;
using Rcpp::is;
using Rcpp::Date;
using Rcpp::Datetime;
using Rcpp::String;
using Rcpp::List;
using Rcpp::DataFrame;
using Rcpp::DateVector;
using Rcpp::DatetimeVector;
using Rcpp::LogicalVector;
using Rcpp::IntegerVector;
using Rcpp::NumericVector;
using Rcpp::ComplexVector;
using Rcpp::CharacterVector;
using Rcpp::warning;
using Rcpp::stop;

//=============================================================================
// CONSTANTS
//=============================================================================
enum R_COL_TYPES {COL_DATE, COL_DATETIME, COL_LOGICAL, COL_INTEGER,
                  COL_NUMERIC, COL_COMPLEX, COL_CHARACTER};

//=============================================================================
// FUNCTIONS
//=============================================================================
// getTypeOfDfCols
//-----------------------------------------------------------------------------
// Description: get the columns types of a dataframe 
//
// Params: 
// - DataFrame: a dataframe
//
// Return: 
// - vector<int>: the type of the columns
//-----------------------------------------------------------------------------
vector<int> getTypeOfDfCols (DataFrame df) {
  vector<int> list_col_types;
  for (int col_idx = 0; col_idx < df.size(); ++col_idx) {
    auto col_values = df[col_idx];
    if ( is<DateVector> (col_values) )
      list_col_types.push_back (COL_DATE);
    else if ( is<DatetimeVector> (col_values) )
      list_col_types.push_back (COL_DATETIME);
    else if ( is<LogicalVector> (col_values) )
      list_col_types.push_back (COL_LOGICAL);
    else if ( is<IntegerVector> (col_values) )
      list_col_types.push_back (COL_INTEGER);
    else if ( is<NumericVector> (col_values) )
      list_col_types.push_back (COL_NUMERIC);
    else if ( is<ComplexVector> (col_values) )
      list_col_types.push_back (COL_COMPLEX);
    else if ( is<CharacterVector> (col_values) )
      list_col_types.push_back (COL_CHARACTER);
    else {
      vector<string> list_cols = as< vector<string> >( df.names() );
      stop ("Error: Type unknown for column " + list_cols[col_idx]); 
      }
    }
  return (list_col_types);
  }

//-----------------------------------------------------------------------------
// cloneDf
//-----------------------------------------------------------------------------
// Description: Clone a DataFrame
//
// Params: 
// - DataFrame: the dataframe to be cloned
// - vector<int> : List of the dataframe columns types
//
// Return: 
// - DataFrame: the clone of the input dataframe
//-----------------------------------------------------------------------------
DataFrame cloneDf (DataFrame df, vector<int> list_cols_types) {
  vector<string> list_cols = as< vector<string> >( df.names() );
  DataFrame df_ret;
  for (int col_idx = 0; col_idx < df.size(); ++col_idx) {
    int col_type = list_cols_types[col_idx]; 
    if (col_type == COL_LOGICAL)
      df_ret.push_back (clone ( LogicalVector (df[col_idx]) ), list_cols[col_idx]);
    else if (col_type == COL_INTEGER)
      df_ret.push_back (clone ( IntegerVector (df[col_idx]) ), list_cols[col_idx]);
    else if (col_type == COL_COMPLEX)
      df_ret.push_back (clone ( ComplexVector (df[col_idx]) ), list_cols[col_idx]);
    else if (col_type == COL_CHARACTER)
      df_ret.push_back (clone ( CharacterVector (df[col_idx]) ), list_cols[col_idx]);
    else // COL_DATE || COL_DATETIME || COL_NUMERIC
      df_ret.push_back (clone ( NumericVector (df[col_idx]) ), list_cols[col_idx]);
    }
  return (df_ret);
  }

//-----------------------------------------------------------------------------
// dropTimestepsCol
//-----------------------------------------------------------------------------
// Description: Drop timesteps from a dataframe column
//
// For each timeseries i, deleted timesteps will be the rows
// having list_drops_start[i] <= row index < list_drops_end[i].
// As the remove function of Rcpp is slow, the drops are performed by copying
// the rows that will remain into a new column
//
// Params: 
// - Dataframe: the dataframe where the modified column will be added
// - string: the column name
// - template class: the column values
// - long: the number of rows after the drops
// - vector<long>: the list of row indexes where the drops start
// - vector<long>: the list of row indexes where the drops end
//-----------------------------------------------------------------------------
template <class T>
void dropTimestepsCol (DataFrame &df, string col_name, T col_values, 
                       long n_rows_after_drop, vector<long> list_drops_start, 
                       vector<long> list_drops_end) {
  T col_copy (n_rows_after_drop);
  long break_idx = 0;
  long drop_start = list_drops_start[break_idx];
  long drop_end = list_drops_end[break_idx];
  long row_copy_idx = -1;
  for (long row_idx = 0; row_idx < col_values.size(); ++row_idx) {
    if (row_idx >= drop_start) {
      if (row_idx < drop_end)
        continue;
      ++break_idx;
      drop_start = list_drops_start[break_idx];
      drop_end = list_drops_end[break_idx];
      if ( (row_idx >= drop_start) && (row_idx < drop_end) )
        continue;
      }
    col_copy[++row_copy_idx] = col_values[row_idx];
    }
  df.push_back (col_copy, col_name);
  }

//-----------------------------------------------------------------------------
// dropTimestepsDf
//-----------------------------------------------------------------------------
// Description: Drop the x first/last timesteps of all timeseries in a dataframe
//
// Params: 
// - DataFrame   : A dataframe
// - vector<int> : List of the dataframe columns types
// - long        : Number of timesteps to delete. 
//                 If > 0, the x first timesteps of each timeseries are dropped. 
//                 When < 0, the x last ones are deleted
// - vector<long>: List of row indexes indicating breaks between the timeseries
//                 This list is updated by dropTimestepsDf to still match
//                 the timeseries breaks after x timesteps have been removed
//
// Return: 
// - DataFrame   : The input dataframe with the x first or last timesteps 
//                 removed from each timeseries
//-----------------------------------------------------------------------------
DataFrame dropTimestepsDf (DataFrame df, vector<int> list_cols_types, 
                           long n_drop, vector<long> &list_breaks) {
  if (n_drop == 0)
    return (df);
  //
  // Compute the row index ranges to be dropped for each timeseries
  //
  vector<long> list_drops_start;
  vector<long> list_drops_end;
  if (n_drop > 0) {
    //
    // When drop is > 0 (first timesteps dropped), the starts are equal to the 
    // breaks and the ends are equal to the breaks + n_drop unless 
    // breaks + n_drop is greater than the next break. In such case, 
    // the end is set to the next break.
    // Note that the list_drops_start and list_drops_end sizes will be the  
    // number of timesteps + 1 as we need one extra to process the timesteps 
    // not dropped of the last timeseries
    //
    list_drops_start = list_breaks;
    long size_drops_start = list_drops_start.size();
    for (long i = 0; i < size_drops_start; ++i) {
      long drop_end = list_drops_start[i] + n_drop;
      if (i == size_drops_start - 1) 
        list_drops_end.push_back (drop_end);
      else {
        if (drop_end > list_breaks[i + 1])
          list_drops_end.push_back (list_breaks[i + 1]);
        else
          list_drops_end.push_back (drop_end);
        }
      }
    }    
  else {
    //
    // When drop is < 0  (last timesteps dropped), the ends are equal to the 
    // breaks and the starts are equal to the breaks + n_drop unless 
    // breaks + n_drop is smaller than the previous break. In such case,
    // the start is set to the previous break.
    //
    list_drops_end = list_breaks;
    list_drops_end.erase ( list_drops_end.begin() );
    
    long size_drops_end = list_drops_end.size();
    for (long i = 0; i < size_drops_end; ++i) {
      long drop_start = list_drops_end[i] + n_drop;
      if (drop_start < list_breaks[i])
        list_drops_start.push_back (list_breaks[i]);
      else
        list_drops_start.push_back (drop_start);
      }
    }
  //
  // Compute the future breaks after the drops and the number of remaining rows
  //
  long n_drop_abs = abs (n_drop);
  long n_rows_after_drop = 0;
  vector<long> list_breaks_upd = {0};
  for (long break_idx = 0; break_idx < list_breaks.size() - 1; ++break_idx) {
    long new_ts_size = list_breaks[break_idx+1] - list_breaks[break_idx] - n_drop_abs;
    if (new_ts_size > 0) {
      n_rows_after_drop += new_ts_size;
      list_breaks_upd.push_back (n_rows_after_drop);
      }
    }
  list_breaks = list_breaks_upd;
  //
  // Apply the drops to each column of the dataframe
  //
  DataFrame df_ret;
  vector<string> list_cols = as<vector<string> >( df.names() );
  for (int col_idx = 0; col_idx < list_cols.size(); ++col_idx) {
    string col_name = list_cols[col_idx];
    int col_type = list_cols_types[col_idx];
    
    if (col_type == COL_LOGICAL)
      dropTimestepsCol <LogicalVector> (df_ret, col_name, LogicalVector (df[col_idx]), 
                              n_rows_after_drop, list_drops_start, list_drops_end);
    else if (col_type == COL_INTEGER)
      dropTimestepsCol <IntegerVector> (df_ret, col_name, IntegerVector (df[col_idx]), 
                              n_rows_after_drop, list_drops_start, list_drops_end);
    else if (col_type == COL_COMPLEX)
      dropTimestepsCol <ComplexVector> (df_ret, col_name, ComplexVector (df[col_idx]), 
                              n_rows_after_drop, list_drops_start, list_drops_end);
    else if (col_type == COL_CHARACTER)
      dropTimestepsCol <CharacterVector> (df_ret, col_name, CharacterVector (df[col_idx]), 
                              n_rows_after_drop, list_drops_start, list_drops_end);
    else // COL_DATE || COL_DATETIME || COL_NUMERIC
      dropTimestepsCol <NumericVector> (df_ret, col_name, NumericVector (df[col_idx]), 
                              n_rows_after_drop, list_drops_start, list_drops_end);
    }
  return (df_ret);
  }

//-----------------------------------------------------------------------------
// movavg
//-----------------------------------------------------------------------------
// Description: Apply moving average operation(s) on a dataframe
//
// Params: 
// - DataFrame   : A dataframe
// - vector<int> : List of dataframe columns types. This list is updated when 
//                 a moving average is applied on an integer or logical column
// - vector<int> : List of moving averages to apply to each column.
//                 Values below 2 are ignored
// - vector<long>: List of row indexes indicating breaks between the timeseries
//                 This list is updated to still match the timeseries breaks
//                 if movavg drops the rows where NAs have been introduced
// - bool        : When true, keep rows with NAs introduced by moving averages
//
// Return: 
// - DataFrame   : The dataframe with moving averages applied
//-----------------------------------------------------------------------------
DataFrame movavg (DataFrame df, vector<int> &list_cols_types, 
                  vector<int> list_movavgs, vector<long> &list_breaks, 
                  bool keep_max_data) {
  //
  // Because moving average modifies the dataframe, take a copy first
  //
  DataFrame df_ret = cloneDf (df, list_cols_types);
  //
  // Apply moving average when the window is > 1 for the column  
  // (but not on the first column which is the timesteps information)
  //
  vector<string> list_cols = as<vector<string> >( df.names() );
  for (int col_idx = 1; col_idx < df.size(); ++col_idx) {
    int movavg = list_movavgs[col_idx-1];
    if (movavg <= 1)
      continue;

    string col_name = list_cols[col_idx];
    auto col_values = df[col_idx];
    int col_type = list_cols_types[col_idx];
    //
    // Check if the type of column allows a moving average
    //
    if (  (col_type == COL_DATE) || (col_type == COL_DATETIME) 
       || (col_type == COL_COMPLEX) || (col_type == COL_CHARACTER) ) {
      warning ("Warning: Moving average request on column " 
               + col_name + " has been ignored");
      list_movavgs[col_idx-1] = -1;
      continue;
      }
    if (col_type == COL_LOGICAL)
      warning ("Warning: Moving average applied on logical column "
               + col_name);
    //
    // Center the moving average window
    //
    int movavg_radius = (int) ( (movavg-1) / 2 );
    int window_start = -movavg_radius;
    int window_end = movavg_radius;
    if (movavg % 2 == 0)
      --window_start;
    //
    // Apply movavg to each timeseries over this column
    // 
    NumericVector col_values_in = col_values;
    NumericVector col_values_ret = df_ret[col_idx];
    bool warning_displayed = false;
    for (long ts_idx = 0; ts_idx < list_breaks.size() - 1; ++ts_idx) {
      //
      // Identify the timesteps range on wich a moving average can be computed
      //
      long first_step = list_breaks[ts_idx] - window_start;
      if (first_step >= list_breaks[ts_idx + 1])
        first_step = list_breaks[ts_idx + 1] - 1;
      long last_step = list_breaks[ts_idx + 1] - 1 - window_end;
      if (last_step < list_breaks[ts_idx])
        last_step = list_breaks[ts_idx];
      
      if ( (first_step > last_step) && (!warning_displayed) ) {
        warning ("Warning: Not enough timesteps to perform a moving average of " 
                 + std::to_string (movavg) + " on variable " + col_name
                 + " (i.e.: timeseries " + std::to_string(ts_idx+1) + ")" );
        warning_displayed = true;
        }
      //
      // Put NAs outside of the range and compute moving averages inside
      //
      long row_idx = list_breaks[ts_idx + 1] - 1;
      for (; row_idx > last_step; --row_idx)
        col_values_ret[row_idx] = NA_REAL;
      for (; row_idx >= first_step; --row_idx) {
        double cumsum = 0;
        for (int i = window_start; i <= window_end; ++i)
          cumsum += col_values_in[row_idx + i];
        col_values_ret[row_idx] = cumsum / movavg;
        }
      for (; row_idx >= list_breaks[ts_idx]; --row_idx)
        col_values_ret[row_idx] = NA_REAL;
      }
    df_ret[col_idx] = col_values_ret;
    list_cols_types[col_idx] = COL_NUMERIC;
    }
  //
  // Remove row(s) where NAs have been introduced because of moving average
  //
  if (!keep_max_data) {
    int max_movavg = *max_element (list_movavgs.begin(), list_movavgs.end());
    int max_movavg_radius = (int) ( (max_movavg-1) / 2 );
    int window_start = -max_movavg_radius;
    int window_end = max_movavg_radius;
    if (max_movavg % 2 == 0)
      --window_start;
    
    df_ret = dropTimestepsDf (df_ret, list_cols_types, -window_start, list_breaks);
    df_ret = dropTimestepsDf (df_ret, list_cols_types, -window_end, list_breaks);
    }
  
  return (df_ret);
  }

//-----------------------------------------------------------------------------
// lagCol
//-----------------------------------------------------------------------------
// Description: lag a dataframe column
//
// As the insert and remove functions of Rcpp are slow, the lags are performed
// by shifting and copying the rows into a new column
//
// Params: 
// - Dataframe: the dataframe where the modified column will be added
// - string: the column name
// - template class: the column values
// - int: the lag to apply
// - vector<long>: the list of timeseries breaks
// - template class: the NA value to use to fill lags 
//-----------------------------------------------------------------------------
template <class T, class N>
void lagCol (DataFrame &df, string col_name, T col_values, int lag, 
             vector<long> list_breaks, N na_value) {
  T col_lagged (col_values.size());
  long col_lagged_idx = -1;
  for (long ts_idx = 0; ts_idx < list_breaks.size() - 1; ++ts_idx) {
    long start_of_ts = list_breaks[ts_idx]; 
    long start_of_next_ts = list_breaks[ts_idx + 1]; 
    long size_of_timeseries = start_of_next_ts - start_of_ts;
    if (lag > size_of_timeseries)
      lag = size_of_timeseries;
    
    long i = 0;
    for (; i < lag; ++i)
      col_lagged[++col_lagged_idx] = na_value;
    for (; i < size_of_timeseries; ++i)
      col_lagged[++col_lagged_idx] = col_values[start_of_ts++];
    }
  df.push_back (col_lagged, col_name);
  }

//-----------------------------------------------------------------------------
// lagData
//-----------------------------------------------------------------------------
// Description: lag the input data for miic in temporal mode
//
// Params: 
// - List: list with an "input_data" entry containing the input dataframe
// - List: list of variable arguments. lagData requires these arguments:
//   * tau          : list of tau, the number of history layers per variable
//                    -1 or 0 if the variable is not lagged
//   * delta_tau    : list of delta_tau, the number of timesteps between each 
//                    history layer. Ignored if the variable is not lagged
//   * movavg       : list of moving averages to apply before lagging data.
//                    Values below 2 are ignored
//   * keep_max_data: if true, keep rows with NA introduced during lagging
//                    or moving average process, otherwise delete these rows 
//
// Return: 
// - List: A list with an "output_data" entry containing the lagged dataframe
//-----------------------------------------------------------------------------
// [[Rcpp::export]]
List lagData (List input_data, List arg_list) {
  vector<int> list_taus = as<vector<int>>(arg_list["tau"]);
  vector<int> list_delta_taus = as<vector<int>>(arg_list["delta_tau"]);
  for (int tau_idx = 0; tau_idx < list_taus.size(); ++tau_idx)
    if (list_taus[tau_idx] < 1)
      list_delta_taus[tau_idx] = 0;
  vector<int> list_movavgs = as<vector<int>>(arg_list["movavg"]);
  bool keep_max_data =  as<bool>(arg_list["keep_max_data"]);
  
  int tau_max = *max_element (list_taus.begin(), list_taus.end());
  
  DataFrame df_in = as<DataFrame> (input_data["input_data"]);
  vector<string> list_cols = as<vector<string> >(df_in.names());
  //
  // As we can not rely durably on the Rccp::is<T>() function to know the
  // type of the columns (i.e.: CharacterVector turns to IntegerVector when
  // stringsAsFactors option is TRUE in R), we will memorize these types 
  // and refers to the initial types when a type check is needed
  //
  vector<int> list_cols_types = getTypeOfDfCols(df_in);
  //
  // Identify the different timeseries, the first and last rows are 
  // considered as breaks, so for n timeseries, we have n+1 breaks
  //
  IntegerVector timesteps = df_in[0];
  vector<long> list_breaks = {0}; 
  long n_timesteps = timesteps.size();
  if (n_timesteps == 0)
      stop ("Error: Input data is empty");
  long prev_ts = timesteps[0] - 1;
  bool warning_displayed = false;
  for (long i = 0; i < n_timesteps; ++i) {
    if (timesteps[i] < prev_ts)
      list_breaks.push_back (i);
    else if ( (prev_ts + 1 != timesteps[i]) && (!warning_displayed) ) {
      warning ("Warning: Timesteps increment differs from 1 (i.e.: timeseries "
               + std::to_string (list_breaks.size()) + ", row " + std::to_string(i+1) + ")" );
      warning_displayed = true;
      }
    prev_ts = timesteps[i]; 
    }
  list_breaks.push_back (n_timesteps);
  //
  // Apply moving average
  //
  if ( *max_element (list_movavgs.begin(), list_movavgs.end()) > 1)
    df_in = movavg (df_in, list_cols_types, list_movavgs, list_breaks, keep_max_data);
  //
  // Check if lag is possible
  //
  vector<int> list_max_lags;
  for (int i = 0; i < list_taus.size(); ++i)
    list_max_lags.push_back (list_taus[i] * list_delta_taus[i]);
  for (int max_lags_idx = 0; max_lags_idx < list_max_lags.size(); ++max_lags_idx) {
    int max_lag = list_max_lags[max_lags_idx];
    for (long ts_idx = 0; ts_idx < list_breaks.size() - 1; ++ts_idx) {
      long n_steps = list_breaks[ts_idx + 1] - list_breaks[ts_idx];
      if (max_lag >= n_steps) {
        warning ("Warning: Not enough timesteps to perform a lag of " + std::to_string (max_lag) 
                 + " steps on variable " + list_cols[max_lags_idx+1] 
                 + " (i.e.: timeseries " + std::to_string(ts_idx+1) + ")" );
        break;          
        }
      }
    }
  //
  // Lagging part, first step: for not lagged and lag0 vars, copy data as it
  //
  list_cols.erase ( list_cols.begin() ); // remove timesteps col
  list_cols_types.erase ( list_cols_types.begin() );
  int n_cols = list_cols.size();
  
  DataFrame df_work;
  for (int col_idx = 0; col_idx < n_cols; ++col_idx) {
    string col_name = list_cols[col_idx];
    if (list_taus[col_idx] > 0)
      col_name = col_name + "_lag0";
    auto col_values = df_in[col_idx+1];
    df_work.push_back (col_values, col_name);
    }
  //    
  // Second step, had extra columns for lagged vars
  //   
  for (int tau_idx = 1; tau_idx <= tau_max; ++tau_idx) {
    for (int col_idx = 0; col_idx < n_cols; ++col_idx) {
      string col_name = list_cols[col_idx];
      if (tau_idx <= list_taus[col_idx]) {
        int col_not_lagged_type = list_cols_types[col_idx];
        list_cols_types.push_back (col_not_lagged_type);
        
        int lag = tau_idx * list_delta_taus[col_idx];
        col_name = col_name + "_lag" + std::to_string (lag);
        auto col_values = df_in[col_idx+1];
        
        if (col_not_lagged_type == COL_LOGICAL)
          lagCol <LogicalVector, int> (df_work, col_name, 
            LogicalVector(col_values), lag, list_breaks, NA_LOGICAL);
        else if (col_not_lagged_type == COL_INTEGER)
          lagCol <IntegerVector, int> (df_work, col_name, 
            IntegerVector(col_values), lag, list_breaks, NA_INTEGER);
        else if (col_not_lagged_type == COL_COMPLEX) {
          Rcomplex NA_COMPLEX;
          NA_COMPLEX.r = NA_REAL;
          NA_COMPLEX.i = NA_REAL;
          lagCol <ComplexVector, Rcomplex> (df_work, col_name, 
            ComplexVector(col_values), lag, list_breaks, NA_COMPLEX);
          }
        else if (col_not_lagged_type == COL_CHARACTER)
          lagCol <CharacterVector, String> (df_work, col_name, 
            CharacterVector(col_values), lag, list_breaks, NA_STRING);
        else // COL_DATE || COL_DATETIME || COL_NUMERIC
          lagCol <NumericVector, double> (df_work, col_name, 
            NumericVector(col_values), lag, list_breaks, NA_REAL);
        }
      }
    }
  //
  // If we drop rows having NA introduced by the lagging, we look for the 
  // maximum lag applied and drop the first "max_lag" rows of each column
  //
  if (!keep_max_data) {
    int max_all_lags = *max_element (list_max_lags.begin(), list_max_lags.end());
    df_work = dropTimestepsDf (df_work, list_cols_types, max_all_lags, list_breaks);
    }
  //
  // Cast back dates and datetimes transformed into numeric
  //
  DataFrame df_out;
  vector<string> list_cols_lagged = as<vector<string> >( df_work.names() );
  for (int col_idx = 0; col_idx < list_cols_lagged.size(); ++col_idx) {
    auto col_values = df_work[col_idx];
    string col_name = list_cols_lagged[col_idx];
    
    if (list_cols_types[col_idx] == COL_DATE) {
      DateVector col_values_typed = col_values;
      std::vector<Date> date_values = col_values_typed.getDates();
      df_out.push_back (date_values, col_name);
      }
    else if (list_cols_types[col_idx] == COL_DATETIME) {
      DatetimeVector col_values_typed = col_values;
      std::vector<Datetime> datetime_values = col_values_typed.getDatetimes();
      df_out.push_back (datetime_values, col_name);
      }
    else
      df_out.push_back (col_values, col_name);
    }
  
  List result = List::create(_["output_data"] = df_out);
  return result;
  }
