
#include "color.hpp"
#include "CZI_ASSERT.hpp"
using namespace ChanZuckerberg::ExpressionMatrix2;

#include "boost_array.hpp"

#include "algorithm.hpp"
#include "iostream.hpp"
#include "sstream.hpp"

// Return a color string given the red,m green, blue components in (0,1).
string ChanZuckerberg::ExpressionMatrix2::color(double red, double green, double blue)
{
    // Clip to (0,1).
    red = max(0., min(red, 1.));
    green = max(0., min(green, 1.));
    blue = max(0., min(blue, 1.));

    // Convbert to integer.
    const int r = int(red * 255);
    const int g = int(green * 255);
    const int b = int(blue * 255);

    // Format them into a Graphviz color.
    ostringstream s;
    s << "#";
    s.fill('0');
    s << hex;
    s.width(2);
    s << r;
    s.width(2);
    s << g;
    s.width(2);
    s << b;

    // Return the corresponding string.
    return s.str();
}



// Return a color string for grey in (0,1).
string ChanZuckerberg::ExpressionMatrix2::grey(double g)
{
    return color(g, g, g);
}



// Return a color string going from red (x=0) to yellow (x=0.5) to green (x=1).
string ChanZuckerberg::ExpressionMatrix2::spectralColor(double x)
{
    double red, green, blue;

    // Always set blue to 0.
    blue = 0;

    // Set red and green based on the value of x.
    if(x <= 0.) {
        // x is zero or less: red.
        red = 1.;
        green = 0.;
    } else if(x<=0.5) {
        // x is between 0 and 0.5.
        red = 1.;
        green = 2. * x;

    } else if(x < 1.) {
        red = 2. * (1.-x);
        green = 1.;
    } else {
        // x is 1 or more: green.
        red = 0.;
        green = 1.;
    }

    return color(red, green, blue);
}



// Return the color string corresponding to a given index
// using the Brewer 12-class qualitative (Set) color scale.
// See colorbrewer2.org.
string ChanZuckerberg::ExpressionMatrix2::brewerSetColor(size_t index)
{

    static const array<string, 12> table =
    {
        "#8dd3c7",
        "#ffffb3",
        "#bebada",
        "#fb8072",
        "#80b1d3",
        "#fdb462",
        "#b3de69",
        "#fccde5",
        "#d9d9d9",
        "#bc80bd",
        "#ccebc5",
        "#ffed6f"
    };
    CZI_ASSERT(index < table.size());

    return table[index];
}
