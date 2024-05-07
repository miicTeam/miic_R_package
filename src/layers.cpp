//******************************************************************************
// Filename   : layer.cpp                             Creation date: 4 may 2024
//
// Description: Layer class implementation
//
// Author     : Franck SIMON
//******************************************************************************

//==============================================================================
// INCLUDE
//==============================================================================
#include <algorithm>
#include <cctype>
#include <string>
#include <sstream>
#include <iostream>

#include "layers.h"

//==============================================================================
// CONSTRUCTORS
//==============================================================================
// Basic constructor
//------------------------------------------------------------------------------
Layer::Layer (int id,
              std::vector<int> contributors,
              std::vector<bool> preoriented)
   {
   this->m_id = id;
   this->m_contributors = contributors;
   this->m_preoriented = preoriented;
   }

//==============================================================================
// PUBLIC FUNCTIONS
//==============================================================================
// Test if a contributor Z is valid for an edge X, Y
//------------------------------------------------------------------------------
bool Layer::is_valid_contributor (Layer &layer_X, Layer &layer_Y, Layer &layer_Z)
  {
  auto itr = std::find (layer_X.m_contributors.begin(),
                        layer_X.m_contributors.end(),
                        layer_Z.m_id);
  if ( itr != layer_X.m_contributors.end() )
    return (true);

  itr = std::find (layer_Y.m_contributors.begin(),
                   layer_Y.m_contributors.end(),
                   layer_Z.m_id);
  if ( itr != layer_Y.m_contributors.end() )
    return (true);

  return (false);
  }

//------------------------------------------------------------------------------
// Indicates if an edge from the layer id must be pre-oriented toward the
// instance layer
//------------------------------------------------------------------------------
bool Layer::is_preoriented (int id)
  {
  auto itr = std::find (this->m_contributors.begin(), this->m_contributors.end(), id);
  if ( itr == this->m_contributors.end() )
    {
    // Rcpp::Rcout << "Wrong call to  Layer::is_preoriented with id " << id
    //             << " for Layer:\n";
    // this->print ();
    return (false);
    }
  auto i = std::distance (this->m_contributors.begin(), itr);
  return (this->m_preoriented[i]);
  }

//------------------------------------------------------------------------------
// Print the instance attributes
//------------------------------------------------------------------------------
void Layer::print ()
  {
  Rcpp::Rcout << "id: " << this->m_id
              << ", contributors:";
  if (this->m_contributors.size() <= 0)
    Rcpp::Rcout << " none";
  else
    {
    for (std::size_t i = 0; i < this->m_contributors.size(); ++i)
      {
      if (this->m_preoriented[i])
        Rcpp::Rcout << " { " << this->m_contributors[i] << " pre-oriented }";
      else
        Rcpp::Rcout << " { " << this->m_contributors[i] << " not pre-oriented }";
      }
    }
  Rcpp::Rcout << "\n";
  }

//==============================================================================
// PUBLIC STATIC FUNCTIONS
//==============================================================================
// Initialization from R layer data frame
//
// Params:
// - a R dataframe with specific columns:
//   * connected: boolean
//   * contributors: a ';' separated string containing the list of layers id
//     of nodes that can be contributors
//   * preoriented: a ';' separated string containing the same number of items
//     the contributors. Indicates for each layer if edges are pre-oriented
//     from the contributor layer toward the layer considered
//   The rows indices in the data frame are considered as layers IDs
// Returns:
// - vector<Layer> : the list of Layer object
//   (with positions in layers vector matching to layer IDs)
//
// NB: All the id are decreased by 1 to match the C++ convention starting from 0
//------------------------------------------------------------------------------
std::vector<Layer> Layer::initFromR (const Rcpp::DataFrame &df)
  {
  Rcpp::CharacterVector contributors_str =  df["contributors"];
  Rcpp::CharacterVector preoriented_str = df["preoriented"];

  std::vector<Layer> layers;
  for (auto i = 0; i < df.nrows(); ++i)
    {
    std::vector<int> contributors;
    if ( ! Rcpp::CharacterVector::is_na (contributors_str[i]) )
      {
      std::stringstream ss_to_split ( Rcpp::as<std::string> (contributors_str[i]) );
      std::string item_str;
      while (getline (ss_to_split, item_str, ';'))
        contributors.push_back ( std::stoi( item_str ) - 1 );
      }

    std::vector<bool> preoriented;
    if ( ! Rcpp::CharacterVector::is_na (preoriented_str[i]) )
      {
      std::stringstream ss_to_split ( Rcpp::as<std::string> (preoriented_str[i]) );
      std::string item_str;
      bool item_bool;
      while (getline (ss_to_split, item_str, ';'))
        {
        // Transform TRUE/FALSE in the string coming from R as true/false
        // so boolalpha will convert correctly into C++ bool
        std::transform (item_str.begin(), item_str.end(), item_str.begin(),
          [](unsigned char c) { return std::tolower(c); });
        // for(auto& c : item_str)
        //   c = tolower(c);
        std::stringstream (item_str) >> std::boolalpha >> item_bool;
        preoriented.push_back (item_bool);
        }
      }

    layers.push_back (Layer (i, contributors, preoriented) );
    }
  return (layers);
  }

//------------------------------------------------------------------------------
// Evaluate if (nodes from) layers are connected
//------------------------------------------------------------------------------
bool Layer::are_connected (Layer &layer1, Layer &layer2)
  {
  auto itr = std::find (layer1.m_contributors.begin(),
                        layer1.m_contributors.end(),
                        layer2.m_id);
  if ( itr != layer1.m_contributors.end() )
    return (true);

  itr = std::find (layer2.m_contributors.begin(),
                   layer2.m_contributors.end(),
                   layer1.m_id);
  if ( itr != layer2.m_contributors.end() )
    return (true);

  return (false);
  }

//------------------------------------------------------------------------------
// Evaluate if layers are pre-oirented
// Return: an int
// - 0 : no pre-orientation
// - >0: pre-oriented from layer1 -> layer2
// - <0: pre-oriented from layer2 -> layer1
//------------------------------------------------------------------------------
int Layer::get_preorientation (Layer &layer1, Layer &layer2)
  {
  if (layer1.is_preoriented (layer2.m_id))
    return (-2);
  if (layer2.is_preoriented (layer1.m_id))
    return (2);
  return (0);
  }
