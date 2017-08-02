#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::plugins(cpp11)]]
//[[Rcpp::depends(dlib)]]
#include<dlib/optimization.h>
//[[Rcpp::depends(BH)]]
#include <boost/math/special_functions/gamma.hpp>
//[[Rcpp::depends(BH)]]
#include <boost/math/special_functions/digamma.hpp>


using namespace std;
using namespace dlib;
using namespace Rcpp;

typedef matrix<double,0,1> column_vector;

struct ICAC{

  private:  
    column_vector m;
    matrix<double> W;

  public:
    static int d;
    static int D;
    static bool generalized;
    column_vector args;
    static matrix<double> X;
    ICAC( const column_vector& m, const matrix<double>& W, const double c = 1. ):
      m(m),
      W(W)
  {
    matrix<double> cMat(1,1);
    cMat(0,0) = c;
    args = reshape_to_column_vector( join_cols( trans(m), W ) );

    if( generalized )
    {
      args = reshape_to_column_vector( join_cols( args, cMat ) );
    }
  }
};

int ICAC::d = 0.; 
int ICAC::D = 0.;
bool ICAC::generalized = 0;
matrix<double> ICAC::X;

double lnl ( const column_vector& args ) // to be substituted with the asymmetry
{
  int d = ICAC::d;
  int D = ICAC::D;
  bool generalized = ICAC::generalized;

  matrix<double> mW; 
  if( generalized )
    mW = subm( args, range(0, args.nr() - 2), range(0,0) );
  else
    mW = args;

  mW = reshape( mW, D+1, D );
  matrix<double> W = rowm(mW, range(1,D));
  matrix<double> m = rowm(mW, 0);

  double c;
  if( ICAC::generalized )
    c= args( args.nr() - 1 );
  else 
    c=2.;

  int n = ICAC::X.nr();

  column_vector s1(D);
  column_vector s2(D);
  column_vector g(d);

  for( int j=0; j < D; j++ )
  {
    s1(j)=0;
    s2(j)=0;

    for( int i=0; i < n; i++ )
    {
      double val = trans( colm(W,j) )*trans( rowm(ICAC::X,i) - m );

      if( val <= 0. )
        s1(j) += pow( abs(val), c );
      else
        s2(j) += pow( abs(val), c );

    }

    if( j<d )
    {
      g(j) = pow(s1(j), 1./(c+1.) ) + pow(s2(j), 1./(c+1.) );
    }
  }

  double result = 1./pow( abs(det(W)), c/(c+1.) );

  for( int j=0; j<d; j++ )
  {
    result = result*g(j);
  }

  for( int j=d; j<D; j++ )
  {
    result = result*pow( (s1(j) + s2(j)), 1./(c+1.) );
  }

  return log(result);
}

// only for the generalized distribution (with c), the other methods are equivalent with lnl
double lnL( const column_vector& _mWC )
{
  double result = 0.;
  double c = _mWC( _mWC.nr() - 1 );
  int n = ICAC::X.nr();
  int d = ICAC::d;
  double e = std::exp(1);
  double kappa = pow( (boost::math::tgamma( 1./c ) )/c, -c );
 
  result += d;
  result = result*log( kappa*n/(c*e) ); 
  column_vector _mW = subm( _mWC, range(0, _mWC.nr() - 2), range(0,0) );

  result = result - ( c+1 )*lnl( _mWC );

  return result;
}

column_vector grad_lnl (const column_vector& args ) 
{
  int d = ICAC::d;
  int D = ICAC::D;
  bool generalized = ICAC::generalized;

  column_vector result = args;

  matrix<double> mW; 
  if( generalized )
    mW = subm( args, range(0, args.nr() - 2), range(0,0) );
  else
    mW = args;

  mW = reshape( mW, D+1, D );
  matrix<double> W = rowm(mW, range(1,D));
  matrix<double> m = rowm(mW, 0);

  double c;
  if( ICAC::generalized )
    c= args( args.nr() - 1 );
  else 
    c=2.;

  matrix<double> inv_tran_W = trans( inv( W ) );

  int n = ICAC::X.nr();

  column_vector s1(D);
  column_vector s2(D);

  column_vector grad_s1(D);
  column_vector grad_s2(D);
 
  column_vector der_s1_c(D);
  column_vector der_s2_c(D);
  
  matrix<double> der_s1 = zeros_matrix<double>(D,D);
  matrix<double> der_s2 = zeros_matrix<double>(D,D);

  for( int j=0; j < D; j++ )
  {
    s1(j)=0;
    s2(j)=0;
    grad_s1(j)=0;
    grad_s2(j)=0;
    der_s1_c(j)=0;
    der_s2_c(j)=0;

    for( int i=0; i < n; i++ )
    {
      double val = trans( colm(W,j) )*trans( rowm(ICAC::X,i) - m );

      if( val <= 0. )
      {
        s1(j) += pow(abs(val), c);
        grad_s1(j) += c*pow( abs(val), c-1. );
        der_s1_c(j) += pow( abs(val), c )*log( abs( val ) );

        for( int k=0; k<D; k++ )
        {
          der_s1(j,k) += c*pow( abs(val), c-1. )*( ICAC::X(i,k) - m(k) );
        }

      }
      else
      {
        s2(j) += pow(abs(val), c);
        grad_s2(j) += c*pow( abs(val), c-1. );
        der_s2_c(j) += pow( abs(val), c )*log( abs( val ) );


        for( int k=0; k<D; k++ )
        {
          der_s2(j,k) += c*pow( abs(val), c-1. )*( ICAC::X(i,k) - m(k) );
        }
      }
    }
  }

  for( int k=0; k < D; k++ ) 
  {
    result(k) = 0.;

    for( int j=0; j < d; j++ )
    {
      double factor =  -1. /( pow( s1(j), 1./(c+1.) ) + pow( s2(j), 1./(c+1.) ) );
      
      double factor1 = 1. /( (c+1) * pow( s1(j), c/(c+1.) ) );
      if( s1(j) == 0. )
        factor1 = 0.;

      double factor2 = 1. /( (c+1) * pow( s2(j), c/(c+1.) ) );
      if( s2(j) == 0. )
        factor2 = 0.;

      result(k) += factor*( factor1*grad_s1(j)*W(k,j) + factor2*grad_s2(j)*W(k,j) ); 
    }

    for( int j=d; j < D; j++ ) // TODO: add c once formulas become available
    {
      double factor = -1. /( 3.*( s1(j) + s2(j) ) );
      if( s1(j) + s2(j) == 0. )
        factor = 0.;

      result(k) += factor*( grad_s1(j)*W(k,j) + grad_s2(j)*W(k,j) ); 
    }
  }

  for( int p=0; p < D; p++ ) 
  {
    for( int k=0; k < D; k++ )
    {
      double factor =  1. /( pow( s1(k), 1./( c + 1. ) ) + pow( s2(k), 1./( c + 1. ) ) );

      double factor1 = 1. /( (c+1) * pow( s1(k), c/( c + 1. ) ) );
      if( s1(k) == 0. )
        factor1 = 0.;


      double factor2 = 1. /( (c+1) * pow( s2(k), c/( c + 1. ) ) );  // TODO: redo indices to match the paper 
      if( s2(k) == 0. )
        factor2 = 0.;

      result( D + p*D + k ) = factor*( -factor1*der_s1(k,p) + factor2*der_s2(k,p) ) - ( c/(c+1.) )*inv_tran_W( p, k ); 

      if( d<D ) // noise terms, TODO: check with Przemek whether this is the correct formula
      {
        result( D + p*D + k ) +=  ( der_s1(k,p) + der_s2(k,p) ) /( 3.*( s1(p) + s2(p) ) ); 
      }

    }
  }

  if( ICAC::generalized )
  {
    result( result.nr() - 1 ) = 0.; 

    for( int j=0; j<D; j++ ) 
    {
      double factor =  1. /( pow( s1(j), 1./( c + 1. ) ) + pow( s2(j), 1./( c + 1. ) ) );
      double sum = 0.;

      if( !( s1(j)==0. ) )
      {
        sum += ( 1./(c+1.) )*der_s1_c(j)*pow( s1(j), -c/(c+1.) );
        sum += -pow( s1(j), 1./(c+1.) ) * log( abs(s1(j)) ) / pow( c+1., 2. );
      }

      if( !( s2(j)==0. ) )
      {
        sum += (1./(c+1.)) * pow( s2(j), -c/(c+1.) ) * der_s2_c(j);
        sum += -pow( s2(j), 1./(c+1.) ) * log( abs(s2(j)) ) / pow( c+1., 2. );
      }

      result( result.nr() - 1 ) += factor*sum;
    }
  }

  return result;
}

column_vector grad_lnL (const column_vector& args ) 
{
  int D = ICAC::D;
  int d = ICAC::d;
  int n = ICAC::X.nr();
  double c = args( args.nr() - 1 );
  double e = std::exp(1);

  //cout << "grad = " << grad_lnl(args) << "\n";
  column_vector result = -(c+1)*grad_lnl( args );

  result( result.nr() - 1 ) += lnl( args )/c;
  result( result.nr() - 1 ) += (log(c*e/n) - 1. + c + boost::math::digamma(1./c) )*d/c;


  return result;

}

// [[Rcpp::export(ICA)]]
RcppExport SEXP ICA( const NumericMatrix& XX, NumericVector& mm, NumericMatrix& WW, double& c, int gauss_noise = 0, bool generalized = 0 ) {
  std::vector< double > _X = as< std::vector<double> >( transpose( XX ) );
  ICAC::X = reshape( mat( _X ), XX.nrow(), XX.ncol() );

  std::vector< double > _W = as< std::vector<double> >( transpose( WW ) );
  matrix<double> W = reshape( mat( _W ), WW.nrow(), WW.ncol() );

  std::vector<double> _m = as< std::vector<double> >(mm);
  column_vector m( _m.size() );
  m = mat( _m );
 
  ICAC::D = W.nr();
  ICAC::d = W.nr() - gauss_noise;
  ICAC::generalized = generalized;

  if( generalized && (gauss_noise > 0.) )
  {
    cout << "Denoising not supported for the generalized model, setting D=d.\n";
    ICAC::d = ICAC::D;
  }

  ICAC my_ica( m, W, c );

  double result = 1.; 

  // GRADIENT CHECKING
  column_vector temp_args = my_ica.args;
/*
  if( ICAC::generalized )
  {
      result = find_min(bfgs_search_strategy(),  // Use BFGS search algorithm
      objective_delta_stop_strategy(1e-7), // Stop when the change in function() is less than 1e-7
      lnL, grad_lnL, my_ica.args, -30.);
  }
  else
  {
      result = find_min(bfgs_search_strategy(),  // Use BFGS search algorithm
      objective_delta_stop_strategy(1e-7), // Stop when the change in function() is less than 1e-7
      lnl, grad_lnl, my_ica.args, -30.);

  }
*/
  //cout << "Solution : " << my_ica.args << " with minimum " << result << "\n \n";

  double approx_result = 1.;

  if( ICAC::generalized )
  {
      approx_result = find_min_using_approximate_derivatives(bfgs_search_strategy(),
                                             objective_delta_stop_strategy(1e-7),
                                             lnL, temp_args, -30.);
  }
  else
  {
      approx_result = find_min_using_approximate_derivatives(bfgs_search_strategy(),
                                             objective_delta_stop_strategy(1e-7),
                                             lnl, temp_args, -30.);
  }


  cout << "Solution : " << temp_args << " with minimum " << approx_result << "\n \n";


 
  if( generalized )
    c = temp_args( temp_args.nr() - 1. );

  matrix<double> args = reshape( temp_args, ICAC::D+1, ICAC::D );

  for( int i = 0; i < ICAC::D; i++ )
  {
    mm( i ) = args( 0, i );

    for( int j = 0; j < ICAC::D; j++ )
    {
      WW(j,i) = args( j+1, i );
    }
  }

  return wrap(result);
}
