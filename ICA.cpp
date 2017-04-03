#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::depends(dlib)]]
#include<dlib/optimization.h>

using namespace std;
using namespace dlib;
using namespace Rcpp;

typedef matrix<double,0,1> column_vector;

class ICAC{

private:  
  matrix<double> X;
  column_vector m;
  matrix<double> W;
  column_vector sigma;
  column_vector tau;
  double n;
  double d;
  matrix<double> Omega;
  column_vector sj1;
  column_vector sj2;


public:
  ICAC()
  {
  }
  ICAC( const matrix<double>& X, const column_vector& m, const matrix<double>& W, const column_vector& sigma, const column_vector& tau ):
    X(X),
    m(m),
    W(W),
    sigma(sigma),
    tau(tau),
    n( X.nr() ),
    d( X.nc() ),
    Omega( d, d )
  {
  }


  void sj_compute()
  {
    double result = 0.;

    for( int i = 0; i < n; i++ )
    {
      result += 0;
    }
  }



  static double my_func ( const column_vector& mW ) // to be substituted with the assymetry
  {

      const double x = mW(0); 
      const double y = mW(1);

      return pow(x,2) + pow(y,2);
  }





  static const column_vector my_func_derivative (const column_vector& mW ) // to be substituted with the assymetry derivative
  /*!
      ensures
          - returns the gradient vector for the rosen function
  !*/
  {
      const double x = mW(0);
      const double y = mW(1);

      // make us a column vector of length 2
      column_vector res(2);

      // now compute the gradient vector
      res(0) = 2*x; // derivative of rosen() with respect to x
      res(1) = 2*y;              // derivative of rosen() with respect to y
      return res;
  }


  pair<double, column_vector> try_algorithm( column_vector starting_point )
  {
      // The other arguments to find_min() are the function to be minimized, its derivative, 
      // then the starting point, and the last is an acceptable minimum value of the function() 
      // function.  That is, if the algorithm finds any inputs to function() that gives an output 
      // value <= -1 then it will stop immediately.  Usually you supply a number smaller than 
      // the actual global minimum.  So since the smallest output of the function() function is 0 
      // we just put -1 here which effectively causes this last argument to be disregarded.

      double result = find_min(bfgs_search_strategy(),  // Use BFGS search algorithm
               objective_delta_stop_strategy(1e-7), // Stop when the change in function() is less than 1e-7
               my_func, my_func_derivative, starting_point, -1);

      return make_pair( result, starting_point );

      // Once the function ends the starting_point vector will contain the optimum point 
  }
};


// [[Rcpp::export]]
double ICA() {
  column_vector startp(2);
  startp = 4, 8;

  ICAC my_ica;
  double min;
  column_vector argmin;

  matrix<double> myMatrix(2,2);
  myMatrix = 1,2,3,4;
  cout << "MY MATRIX:\n" << myMatrix << "\n";
  matrix<double> myMatrix2 = reshape_to_column_vector( myMatrix );
  cout << "MY MATRIX:\n" << myMatrix2 << "\n";
  matrix<double> myMatrix3 = join_cols( myMatrix2, startp );
  cout << "MY MATRIX:\n" << myMatrix3 << "\n";

  
  pair<double, column_vector> ica_args = my_ica.try_algorithm( startp );

  cout << "Function attains minimum:\n" << ica_args.second <<  ", which equals to: " << ica_args.first << endl;

  return ica_args.first;
}
