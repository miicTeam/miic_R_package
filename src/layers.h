//******************************************************************************
// Filename   : layer.h                               Creation date: 4 may 2024
//
// Description: Layer class include
//
// Author     : Franck SIMON
//******************************************************************************
#ifndef MIIC_LAYERS
#define MIIC_LAYERS

//==============================================================================
// INCLUDE
//==============================================================================
#include <vector>
#include <Rcpp.h>

//==============================================================================
// CLASS
//==============================================================================
class Layer
  {
  public:
    //--------------------------------------------------------------------------
    // Basic constructor
    //--------------------------------------------------------------------------
    Layer (int id,
           std::vector<int> contributors,
           std::vector<bool> preoriented);

    //--------------------------------------------------------------------------
    // Public member functions
    //--------------------------------------------------------------------------
    // Indicates if nodes from layer id induce a pre-orientation toward the
    // instance layer
    //
    bool is_preoriented (int id);
    //
    // Print the instance
    //
    void print();

    //--------------------------------------------------------------------------
    // Public static functions
    //--------------------------------------------------------------------------
    // Initialize a vector of layers from the R data frame
    //
    static std::vector<Layer> initFromR (const Rcpp::DataFrame &);
    //
    // Evaluate if layers are connected
    //
    static bool are_connected (Layer &layer1, Layer &layer2);
    //
    // Test if a layer Z is valid to look for contributors in regard of an edge
    // having its nodes in layer X and Y
    //
    static bool is_valid_contributor (Layer &layer_X, Layer &layer_Y, Layer &layer_Z);
    //
    // Evaluate if layers are pre-oirented
    //
    static int get_preorientation (Layer &layer1, Layer &layer2);

  private:
    //--------------------------------------------------------------------------
    // Private attributes
    //--------------------------------------------------------------------------
    // Id of the layer
    //
    int m_id;
    //
    // List of layers from which nodes can be contributors
    //
    std::vector <int> m_contributors;
    //
    // Same order as m_contributors, indicates if nodes from each contributor
    // layer induce a pre-orientation toward the instance layer
    //
    std::vector <bool> m_preoriented;
  };

#endif
