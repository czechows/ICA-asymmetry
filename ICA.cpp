#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::depends(dlib)]]
#include<dlib/optimization.h>

using namespace std;
using namespace dlib;
using namespace Rcpp;

typedef matrix<double,0,1> column_vector;

double my_func (const column_vector& m) // to be substituted with the assymetry
{
    const double x = m(0); 
    const double y = m(1);

    return pow(x,2) + pow(y,2);
}


const column_vector my_func_derivative (const column_vector& m) // to be substituted with the assymetry derivative
/*!
    ensures
        - returns the gradient vector for the rosen function
!*/
{
    const double x = m(0);
    const double y = m(1);

    // make us a column vector of length 2
    column_vector res(2);

    // now compute the gradient vector
    res(0) = 2*x; // derivative of rosen() with respect to x
    res(1) = 2*y;              // derivative of rosen() with respect to y
    return res;
}

void try_algorithm()
{
    column_vector starting_point(2);
    starting_point = 4, 8;

    // The first example below finds the minimum of the rosen() function and uses the
    // analytical derivative computed by rosen_derivative().  Since it is very easy to
    // make a mistake while coding a function like rosen_derivative() it is a good idea
    // to compare your derivative function against a numerical approximation and see if
    // the results are similar.  If they are very different then you probably made a 
    // mistake.  So the first thing we do is compare the results at a test point: 
 //   cout << "Difference between analytic derivative and numerical approximation of derivative: " 
 //         << length(derivative(rosen)(starting_point) - rosen_derivative(starting_point)) << endl;


    // Now we use the find_min() function to find the minimum point.  The first argument
    // to this routine is the search strategy we want to use.  The second argument is the 
    // stopping strategy.  Below I'm using the objective_delta_stop_strategy which just 
    // says that the search should stop when the change in the function being optimized 
    // is small enough.

    // The other arguments to find_min() are the function to be minimized, its derivative, 
    // then the starting point, and the last is an acceptable minimum value of the function() 
    // function.  That is, if the algorithm finds any inputs to function() that gives an output 
    // value <= -1 then it will stop immediately.  Usually you supply a number smaller than 
    // the actual global minimum.  So since the smallest output of the function() function is 0 
    // we just put -1 here which effectively causes this last argument to be disregarded.

    find_min(bfgs_search_strategy(),  // Use BFGS search algorithm
             objective_delta_stop_strategy(1e-7), // Stop when the change in function() is less than 1e-7
             my_func, my_func_derivative, starting_point, -1);

    // Once the function ends the starting_point vector will contain the optimum point 
    // of (1,1).
    cout << "Function minimum:\n" << starting_point << endl;
}

// [[Rcpp::export]]
void ICA() {
  try_algorithm();
}
