#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::plugins(cpp11)]]
//[[Rcpp::depends(dlib)]]
#include<dlib/optimization.h>

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
    column_vector mW;
    static matrix<double> X;
    ICAC( const column_vector& m, const matrix<double>& W ):
      m(m),
      W(W)
  {
    mW = reshape_to_column_vector( join_cols( trans(m), W ) );
  }

};

int ICAC::d = 0.; 
int ICAC::D = 0.;
matrix<double> ICAC::X;

double lnl ( const column_vector& _mW ) // to be substituted with the asymmetry
{
  int d = ICAC::d;
  int D = ICAC::D;

  matrix<double> mW = reshape( _mW, D+1, D );
  matrix<double> W = rowm(mW, range(1,D));
  matrix<double> m = rowm(mW, 0);

  int n = ICAC::X.nr();

  column_vector s1(D);
  column_vector s2(D);
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

  for( int j=d; j<D; j++ )
  {
    result = result*pow( s1(j) + s2(j), 1./3. );
  }

  return log(result);
}

const column_vector grad_lnl (const column_vector& _mW ) 
{
  int d = ICAC::d;
  int D = ICAC::D;

  column_vector result = _mW;

  matrix<double> mW = reshape( _mW, D+1, D );
  matrix<double> W = rowm(mW, range(1,D));
  matrix<double> m = rowm(mW, 0);

  matrix<double> inv_tran_W = trans( inv( W ) );

  int n = ICAC::X.nr();

  column_vector s1(D);
  column_vector s2(d);

  column_vector grad_s1(D);
  column_vector grad_s2(D);

  matrix<double> der_s1 = zeros_matrix<double>(D,D);
  matrix<double> der_s2 = zeros_matrix<double>(D,D);

  for( int j=0; j < D; j++ )
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

        for( int k=0; k<D; k++ )
        {
          der_s1(j,k) += 2*val*( ICAC::X(i,k) - m(k) );
        }

      }
      else
      {
        s2(j) += val*val;
        grad_s2(j) += 2*val;


        for( int k=0; k<D; k++ )
        {
          der_s2(j,k) += 2*val*( ICAC::X(i,k) - m(k) );
        }
      }
    }
  }

  for( int k=0; k < D; k++ ) 
  {
    result(k) = 0.;

    for( int j=0; j < d; j++ )
    {
      double factor =  -1. /( pow( s1(j), 1./3. ) + pow( s2(j), 1./3. ) );
      
      double factor1 = 1. /( 3. * pow( s1(j), 2./3. ) );
      if( s1(j) == 0. )
        factor1 = 0.;

      double factor2 = 1. /( 3. * pow( s2(j), 2./3. ) );
      if( s2(j) == 0. )
        factor2 = 0.;

      result(k) += factor*( factor1*grad_s1(j)*W(k,j) + factor2*grad_s2(j)*W(k,j) ); 
    }

    for( int j=d; j < D; j++ )
    {
      double factor = -1. /( 3.*( s1(j) + s2(j) ) );

      result(k) += factor*( grad_s1(j)*W(k,j) + grad_s2(j)*W(k,j) ); 
    }
  }

  for( int p=0; p < d; p++ ) 
  {
    for( int k=0; k < d; k++ )
    {
      double factor =  1. /( pow( s1(k), 1./3. ) + pow( s2(k), 1./3. ) );

      double factor1 = 1. /( 3. * pow( s1(k), 2./3. ) );
      if( s1(k) == 0. )
        factor1 = 0.;


      double factor2 = 1. /( 3. * pow( s2(k), 2./3. ) );  // TODO: redo indices to match the paper 
      if( s2(k) == 0. )
        factor2 = 0.;

      result( d + p*d + k ) = factor*( factor1*der_s1(k,p) + factor2*der_s2(k,p) ) - (2./3.)*inv_tran_W( p, k ); 
      result( d + p*d + k ) +=  ( der_s1(k,p) + der_s2(k,p) ) /( 3.*( s1(p) + s2(p) ) ); 

    }
  }

  return result;
}

// [[Rcpp::export(ICA)]]
RcppExport SEXP ICA( const NumericMatrix& XX, NumericVector& mm, NumericMatrix& WW, int gauss_noise = 0. ) {
  std::vector< double > _X = as< std::vector<double> >( transpose( XX ) );
  ICAC::X = reshape( mat( _X ), XX.nrow(), XX.ncol() );

  std::vector< double > _W = as< std::vector<double> >( transpose( WW ) );
  matrix<double> W = reshape( mat( _W ), WW.nrow(), WW.ncol() );

  std::vector<double> _m = as< std::vector<double> >(mm);
  column_vector m( _m.size() );
  m = mat( _m );
 
  ICAC::D = W.nr();
  ICAC::d = W.nr() - gauss_noise;

  ICAC my_ica( m, W );

  column_vector test = my_ica.mW;

  double result = find_min(bfgs_search_strategy(),  // Use BFGS search algorithm
      objective_delta_stop_strategy(1e-7), // Stop when the change in function() is less than 1e-7
      lnl, grad_lnl, my_ica.mW, -1);

  matrix<double> mW = reshape( my_ica.mW, ICAC::D+1, ICAC::D );

  for( int i = 0; i < ICAC::D; i++ )
  {
    mm( i ) = mW( 0, i );

    for( int j = 0; j < ICAC::D; j++ )
    {
      WW(j,i) = mW( j+1, i );
    }
  }

  return wrap(result);
}
