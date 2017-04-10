#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::depends(dlib)]]
#include<dlib/optimization.h>

using namespace std;
using namespace dlib;
using namespace Rcpp;


typedef matrix<double,0,1> column_vector;

// TODO: get rid of this struct?
struct ICAC{

  private:  
    column_vector m;
    matrix<double> W;

  public:
    static double d;
    column_vector mW;
    static matrix<double> X;
    ICAC( const column_vector& m, const matrix<double>& W ):
      m(m),
      W(W)
  {
    mW = reshape_to_column_vector( join_cols( trans(m), W ) );
  }

};


double ICAC::d = 0.; 
matrix<double> ICAC::X;


double lnl ( const column_vector& _mW ) // to be substituted with the assymetry
{
  double d = ICAC::d;

  matrix<double> mW = reshape( _mW, d+1, d );
  matrix<double> W = rowm(mW, range(1,d));
  matrix<double> m = rowm(mW, 0);

  int n = ICAC::X.nr();

  column_vector s1(d);
  column_vector s2(d);
  column_vector g(d);

  for( int j=0; j < d; j++ )
  {
    s1(j)=0;
    s2(j)=0;

    for( int i=0; i < n; i++ )
    {
      double val = trans( colm(W,j) )*trans( rowm(ICAC::X,i) - m );

      if( val <= 0. )
        s1(j) += val*val;
      else
        s2(j) += val*val;

    }

    g(j) = pow(s1(j), 1./3.) + pow(s2(j), 1./3.);
  }

  double result = 1./pow( abs(det(W)), 2./3. );

  for( int j=0; j<d; j++ )
  {
    result = result*g(j);
  }

  return log(result);
}





column_vector grad_lnl (const column_vector& _mW ) // to be substituted with the assymetry derivative
{
  double d = ICAC::d;

  column_vector result = _mW;

  matrix<double> mW = reshape( _mW, d+1, d );
  matrix<double> W = rowm(mW, range(1,d));
  matrix<double> m = rowm(mW, 0);

  matrix<double> inv_tran_W = trans( inv( W ) );

  int n = ICAC::X.nr();

  column_vector s1(d);
  column_vector s2(d);

  column_vector grad_s1(d);
  column_vector grad_s2(d);

  matrix<double> der_s1 = zeros_matrix<double>(d,d);
  matrix<double> der_s2 = zeros_matrix<double>(d,d);

  for( int j=0; j < d; j++ )
  {
    s1(j)=0;
    s2(j)=0;
    grad_s1(j)=0;
    grad_s2(j)=0;

    for( int i=0; i < n; i++ )
    {
      double val = trans( colm(W,j) )*trans( rowm(ICAC::X,i) - m );

      if( val <= 0. )
      {
        s1(j) += val*val;
        grad_s1(j) += 2*val;

        for( int k=0; k<d; k++ )
        {
          der_s1(j,k) += 2*val*( ICAC::X(i,k) - m(k) );
        }

      }
      else
      {
        s2(j) += val*val;
        grad_s2(j) += 2*val;


        for( int k=0; k<d; k++ )
        {
          der_s2(j,k) += 2*val*( ICAC::X(i,k) - m(k) );
        }
      }
    }
  }

  for( int k=0; k < d; k++ ) // TODO: rewrite as matrix multiplication (vectorize)
  {
    result(k) = 0.;

    for( int j=0; j < d; j++ )
    {
      double factor =  -1. /( pow( s1(j), 1./3. ) + pow( s2(j), 1./3. ) );
      double factor1 = 1. /( 3. * pow( s1(j), 2./3. ) );
      double factor2 = 1. /( 3. * pow( s2(j), 2./3. ) );

      result(k) += factor*( factor1*grad_s1(j)*W(k,j) + factor2*grad_s2(j)*W(k,j) );
    }
  }


  for( int p=0; p < d; p++ ) // TODO: rewrite as matrix multiplication (vectorize)
  {
    for( int k=0; k < d; k++ )
    {
      double factor =  1. /( pow( s1(p), 1./3. ) + pow( s2(p), 1./3. ) );
      double factor1 = 1. /( 3. * pow( s1(p), 2./3. ) );
      double factor2 = 1. /( 3. * pow( s2(p), 2./3. ) );

      result( d + k*d + p ) = factor*( factor1*der_s1(p,k) + factor2*der_s2(p,k) ) - (2./3.)*inv_tran_W( k, p ); 
    }
  }

  return result;
}







// [[Rcpp::export]]
RcppExport SEXP ICA( NumericMatrix XX, NumericVector mm, NumericMatrix WW ) {

  std::vector< double > _X = as< std::vector<double> >(XX);
  ICAC::X = reshape( mat( _X ), XX.nrow(), XX.ncol() );

  std::vector< double > _W = as< std::vector<double> >(WW);
  matrix<double> W = reshape( mat( _W ), WW.nrow(), WW.ncol() );

  std::vector<double> _m = as< std::vector<double> >(mm);
  column_vector m( _m.size() );
  m = mat( _m );

  ICAC::d = W.nr();
  ICAC my_ica( m ,W );

  // TODO: erase this
  cout << "lnl = " << lnl( my_ica.mW ) << "\n";
  cout << "der lnl = \n" << derivative(lnl, 1e-5)( my_ica.mW ) << "\n";
  cout << "grad_lnl = \n" << grad_lnl( my_ica.mW ) << "\n";

  double result = 0.;
  // TODO: this has to work
  /*
  double result = find_min(bfgs_search_strategy(),  // Use BFGS search algorithm
      objective_delta_stop_strategy(1e-7), // Stop when the change in function() is less than 1e-7
      lnl, grad_lnl, my_ica.mW, -1);

  cout << "minimum = " << result << "\n";
  */

  return wrap(result);
}
