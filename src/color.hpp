// Functions to generate color strings that can be used with Graphviz and HTML.


#ifndef CZI_EXPRESSION_MATRIX2_COLOR_HPP
#define CZI_EXPRESSION_MATRIX2_COLOR_HPP

#include "string.hpp"

namespace ChanZuckerberg {
    namespace ExpressionMatrix2 {

        // Return a color string given the red, green, blue components in (0,1).
        string color(double red, double green, double blue);

        // Return a color string for grey in (0,1).
        string grey(double grey);

        // Return a color string going from red (x=0) to yellow (x=0.5) to green (x=1).
        string spectralColor(double x);

        // Return the color string corresponding to a given index
        // using the Brewer 12-class qualitative (Set) color scale.
        // See colorbrewer2.org.
        string brewerSetColor(size_t index);

        // I created this by hand. It uses brighter colors than brewerSetColor,
        // but it does not use red and black which are used for special purposes.
        // It also has 12 colors.
        string colorPalette1(size_t index);
    }
}

#endif

