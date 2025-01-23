#' Style writing function for the miic network
#' @description This function writes the  miic style for a correct visualization using the cytoscape tool (http://www.cytoscape.org/).
#' @details The style is written in the xml file format.
#' @param file [a string] The file path of the output file (containing the file name without extension).
#' @export
#' @useDynLib miic

miic.write.style.cytoscape <- function(file) {
  if (missing(file)) {
    cat("The file path is necessary")
  } else {
    str <- "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>
    <vizmap documentVersion=\"3.0\" id=\"VizMap-2017_02_27-15_04\">
      <visualStyle name=\"miic_style\">
          <network>
              <visualProperty name=\"NETWORK_BACKGROUND_PAINT\" default=\"#FFFFFF\"/>
              <visualProperty name=\"NETWORK_EDGE_SELECTION\" default=\"true\"/>
              <visualProperty name=\"NETWORK_CENTER_X_LOCATION\" default=\"0.0\"/>
              <visualProperty name=\"NETWORK_CENTER_Y_LOCATION\" default=\"0.0\"/>
              <visualProperty name=\"NETWORK_DEPTH\" default=\"0.0\"/>
              <visualProperty name=\"NETWORK_HEIGHT\" default=\"400.0\"/>
              <visualProperty name=\"NETWORK_TITLE\" default=\"\"/>
              <visualProperty name=\"NETWORK_WIDTH\" default=\"550.0\"/>
              <visualProperty name=\"NETWORK_SCALE_FACTOR\" default=\"1.0\"/>
              <visualProperty name=\"NETWORK_CENTER_Z_LOCATION\" default=\"0.0\"/>
              <visualProperty name=\"NETWORK_NODE_SELECTION\" default=\"true\"/>
          </network>
          <node>
              <dependency name=\"nodeCustomGraphicsSizeSync\" value=\"true\"/>
              <dependency name=\"nodeSizeLocked\" value=\"true\"/>
              <visualProperty name=\"NODE_CUSTOMGRAPHICS_1\" default=\"org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],\"/>
              <visualProperty name=\"NODE_CUSTOMGRAPHICS_POSITION_3\" default=\"C,C,c,0.00,0.00\"/>
              <visualProperty name=\"NODE_TOOLTIP\" default=\"\"/>
              <visualProperty name=\"NODE_CUSTOMGRAPHICS_2\" default=\"org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],\"/>
              <visualProperty name=\"COMPOUND_NODE_PADDING\" default=\"10.0\"/>
              <visualProperty name=\"NODE_CUSTOMGRAPHICS_7\" default=\"org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],\"/>
              <visualProperty name=\"NODE_PAINT\" default=\"#787878\"/>
              <visualProperty name=\"NODE_LABEL_FONT_FACE\" default=\"Ubuntu,plain,12\"/>
              <visualProperty name=\"NODE_LABEL_FONT_SIZE\" default=\"14\"/>
              <visualProperty name=\"NODE_SIZE\" default=\"40.0\"/>
              <visualProperty name=\"NODE_CUSTOMPAINT_7\" default=\"DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_7, name=Node Custom Paint 7)\"/>
              <visualProperty name=\"NODE_CUSTOMGRAPHICS_9\" default=\"org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],\"/>
              <visualProperty name=\"NODE_LABEL_TRANSPARENCY\" default=\"255\"/>
              <visualProperty name=\"NODE_TRANSPARENCY\" default=\"180\"/>
              <visualProperty name=\"NODE_CUSTOMPAINT_9\" default=\"DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_9, name=Node Custom Paint 9)\"/>
              <visualProperty name=\"NODE_CUSTOMGRAPHICS_8\" default=\"org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],\"/>
              <visualProperty name=\"NODE_BORDER_PAINT\" default=\"#666666\"/>
              <visualProperty name=\"NODE_CUSTOMGRAPHICS_5\" default=\"org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],\"/>
              <visualProperty name=\"NODE_CUSTOMGRAPHICS_POSITION_5\" default=\"C,C,c,0.00,0.00\"/>
              <visualProperty name=\"NODE_LABEL_COLOR\" default=\"#000000\"/>
              <visualProperty name=\"NODE_FILL_COLOR\" default=\"#f7f7f7\"/>
              <visualProperty name=\"NODE_CUSTOMGRAPHICS_POSITION_1\" default=\"C,C,c,0.00,0.00\"/>
              <visualProperty name=\"NODE_CUSTOMGRAPHICS_POSITION_2\" default=\"C,C,c,0.00,0.00\"/>
              <visualProperty name=\"COMPOUND_NODE_SHAPE\" default=\"ROUND_RECTANGLE\"/>
              <visualProperty name=\"NODE_Z_LOCATION\" default=\"0.0\"/>
              <visualProperty name=\"NODE_BORDER_STROKE\" default=\"SOLID\"/>
              <visualProperty name=\"NODE_LABEL_WIDTH\" default=\"200.0\"/>
              <visualProperty name=\"NODE_LABEL\" default=\"\">
                  <passthroughMapping attributeType=\"string\" attributeName=\"name\"/>
              </visualProperty>
              <visualProperty name=\"NODE_CUSTOMGRAPHICS_SIZE_2\" default=\"0.0\"/>
              <visualProperty name=\"NODE_BORDER_WIDTH\" default=\"2.0\"/>
              <visualProperty name=\"NODE_SELECTED_PAINT\" default=\"#FFFF00\"/>
              <visualProperty name=\"NODE_CUSTOMGRAPHICS_SIZE_9\" default=\"0.0\"/>
              <visualProperty name=\"NODE_Y_LOCATION\" default=\"0.0\"/>
              <visualProperty name=\"NODE_CUSTOMGRAPHICS_POSITION_6\" default=\"C,C,c,0.00,0.00\"/>
              <visualProperty name=\"NODE_CUSTOMGRAPHICS_4\" default=\"org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],\"/>
              <visualProperty name=\"NODE_DEPTH\" default=\"0.0\"/>
              <visualProperty name=\"NODE_CUSTOMGRAPHICS_SIZE_4\" default=\"0.0\"/>
              <visualProperty name=\"NODE_SHAPE\" default=\"ELLIPSE\"/>
              <visualProperty name=\"NODE_CUSTOMPAINT_6\" default=\"DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_6, name=Node Custom Paint 6)\"/>
              <visualProperty name=\"NODE_CUSTOMPAINT_1\" default=\"DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_1, name=Node Custom Paint 1)\"/>
              <visualProperty name=\"NODE_CUSTOMGRAPHICS_SIZE_8\" default=\"0.0\"/>
              <visualProperty name=\"NODE_CUSTOMPAINT_2\" default=\"DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_2, name=Node Custom Paint 2)\"/>
              <visualProperty name=\"NODE_CUSTOMGRAPHICS_SIZE_3\" default=\"0.0\"/>
              <visualProperty name=\"NODE_CUSTOMGRAPHICS_SIZE_6\" default=\"0.0\"/>
              <visualProperty name=\"NODE_CUSTOMPAINT_3\" default=\"DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_3, name=Node Custom Paint 3)\"/>
              <visualProperty name=\"NODE_CUSTOMGRAPHICS_SIZE_7\" default=\"0.0\"/>
              <visualProperty name=\"NODE_NESTED_NETWORK_IMAGE_VISIBLE\" default=\"true\"/>
              <visualProperty name=\"NODE_CUSTOMGRAPHICS_6\" default=\"org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],\"/>
              <visualProperty name=\"NODE_CUSTOMPAINT_4\" default=\"DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_4, name=Node Custom Paint 4)\"/>
              <visualProperty name=\"NODE_VISIBLE\" default=\"true\"/>
              <visualProperty name=\"NODE_CUSTOMGRAPHICS_3\" default=\"org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],\"/>
              <visualProperty name=\"NODE_CUSTOMPAINT_8\" default=\"DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_8, name=Node Custom Paint 8)\"/>
              <visualProperty name=\"NODE_CUSTOMGRAPHICS_POSITION_9\" default=\"C,C,c,0.00,0.00\"/>
              <visualProperty name=\"NODE_CUSTOMGRAPHICS_SIZE_5\" default=\"0.0\"/>
              <visualProperty name=\"NODE_BORDER_TRANSPARENCY\" default=\"255\"/>
              <visualProperty name=\"NODE_HEIGHT\" default=\"40.0\"/>
              <visualProperty name=\"NODE_WIDTH\" default=\"60.0\"/>
              <visualProperty name=\"NODE_CUSTOMGRAPHICS_POSITION_4\" default=\"C,C,c,0.00,0.00\"/>
              <visualProperty name=\"NODE_LABEL_POSITION\" default=\"C,C,c,0.00,0.00\"/>
              <visualProperty name=\"NODE_X_LOCATION\" default=\"0.0\"/>
              <visualProperty name=\"NODE_CUSTOMPAINT_5\" default=\"DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_5, name=Node Custom Paint 5)\"/>
              <visualProperty name=\"NODE_CUSTOMGRAPHICS_SIZE_1\" default=\"0.0\"/>
              <visualProperty name=\"NODE_CUSTOMGRAPHICS_POSITION_7\" default=\"C,C,c,0.00,0.00\"/>
              <visualProperty name=\"NODE_CUSTOMGRAPHICS_POSITION_8\" default=\"C,C,c,0.00,0.00\"/>
              <visualProperty name=\"NODE_SELECTED\" default=\"false\"/>
          </node>
          <edge>
              <dependency name=\"arrowColorMatchesEdge\" value=\"false\"/>
              <visualProperty name=\"EDGE_SELECTED\" default=\"false\"/>
              <visualProperty name=\"EDGE_LABEL_FONT_FACE\" default=\"SansSerif.plain,plain,10\"/>
              <visualProperty name=\"EDGE_STROKE_UNSELECTED_PAINT\" default=\"#808080\">
                  <discreteMapping attributeType=\"string\" attributeName=\"sign\">
                      <discreteMappingEntry value=\"#ff3300\" attributeValue=\"+\"/>
                      <discreteMappingEntry value=\"#9999FF\" attributeValue=\"-\"/>
                  </discreteMapping>
              </visualProperty>
              <visualProperty name=\"EDGE_WIDTH\" default=\"10.0\">
                  <continuousMapping attributeType=\"float\" attributeName=\"weight\">
                      <continuousMappingPoint lesserValue=\"2.0\" greaterValue=\"2.0\" equalValue=\"2.0\" attrValue=\"0.0\"/>
                      <continuousMappingPoint lesserValue=\"12.0\" greaterValue=\"12.0\" equalValue=\"12.0\" attrValue=\"20.0\"/>
                  </continuousMapping>
              </visualProperty>
              <visualProperty name=\"EDGE_TARGET_ARROW_SHAPE\" default=\"NONE\">
                  <discreteMapping attributeType=\"string\" attributeName=\"targetArrowShape\">
                      <discreteMappingEntry value=\"T\" attributeValue=\"T\"/>
                      <discreteMappingEntry value=\"ARROW\" attributeValue=\"arrow\"/>
                      <discreteMappingEntry value=\"NONE\" attributeValue=\"none\"/>
                  </discreteMapping>
              </visualProperty>
              <visualProperty name=\"EDGE_LINE_TYPE\" default=\"SOLID\">
                  <discreteMapping attributeType=\"integer\" attributeName=\"edgeType\">
                      <discreteMappingEntry value=\"SOLID\" attributeValue=\"1\"/>
                      <discreteMappingEntry value=\"SOLID\" attributeValue=\"2\"/>
                      <discreteMappingEntry value=\"BACKWARD_SLASH\" attributeValue=\"6\"/>
                  </discreteMapping>
              </visualProperty>
              <visualProperty name=\"EDGE_CURVED\" default=\"true\"/>
              <visualProperty name=\"EDGE_SOURCE_ARROW_SHAPE\" default=\"NONE\">
                  <discreteMapping attributeType=\"string\" attributeName=\"sourceArrowShape\">
                      <discreteMappingEntry value=\"T\" attributeValue=\"T\"/>
                      <discreteMappingEntry value=\"ARROW\" attributeValue=\"arrow\"/>
                      <discreteMappingEntry value=\"NONE\" attributeValue=\"none\"/>
                  </discreteMapping>
              </visualProperty>
              <visualProperty name=\"EDGE_LABEL\" default=\"\"/>
              <visualProperty name=\"EDGE_SOURCE_ARROW_SELECTED_PAINT\" default=\"#FFFF00\"/>
              <visualProperty name=\"EDGE_TRANSPARENCY\" default=\"230\"/>
              <visualProperty name=\"EDGE_VISIBLE\" default=\"true\"/>
              <visualProperty name=\"EDGE_LABEL_COLOR\" default=\"#000000\"/>
              <visualProperty name=\"EDGE_LABEL_WIDTH\" default=\"200.0\"/>
              <visualProperty name=\"EDGE_UNSELECTED_PAINT\">
                  <discreteMapping attributeType=\"string\" attributeName=\"sign\">
                      <discreteMappingEntry value=\"#919191\" attributeValue=\"+\"/>
                      <discreteMappingEntry value=\"#2B9FFF\" attributeValue=\"-\"/>
                  </discreteMapping>
              </visualProperty>
              <visualProperty name=\"EDGE_SOURCE_ARROW_UNSELECTED_PAINT\" default=\"#000000\"/>
              <visualProperty name=\"EDGE_TARGET_ARROW_UNSELECTED_PAINT\" default=\"#000000\"/>
              <visualProperty name=\"EDGE_LABEL_TRANSPARENCY\" default=\"255\"/>
              <visualProperty name=\"EDGE_LABEL_FONT_SIZE\" default=\"10\"/>
              <visualProperty name=\"EDGE_BEND\"/>
              <visualProperty name=\"EDGE_TARGET_ARROW_SELECTED_PAINT\" default=\"#FFFF00\"/>
              <visualProperty name=\"EDGE_TOOLTIP\" default=\"\"/>
              <visualProperty name=\"EDGE_STROKE_SELECTED_PAINT\" default=\"#FF0000\"/>
          </edge>
      </visualStyle>
  </vizmap>"

    write(x = str, file = paste(file, ".xml", sep = ""))
  }
}
