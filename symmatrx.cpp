
#include <math.h>

using namespace std;

#include "symmatrx.hpp"

SymMatrix::SymMatrix( int newSize ) :
   size(0),
   nelement(0),
   maxsize(0),
   scratch(0),
   value(0)
   {
   assert( newSize > 0 );
   Zero( newSize );
   }

SymMatrix::~SymMatrix() {
   if( scratch ) delete []  scratch;
   }

void SymMatrix::SetSize( int newSize ) {
   size = newSize;
   nelement = (size*(size+1))/2;
   if( newSize <= maxsize ) return;
   if( value ) delete [] scratch;
   maxsize = newSize;
   scratch = new double[nelement+size];
   value = scratch + size;
   }

void SymMatrix::Zero( int newSize ) {
   if( newSize != 0 ) SetSize(newSize);
   for( int i = 0; i < nelement; i++ ) value[i] = 0.0;
   }

void SymMatrix::ScaleBy( double factor ) {
   for( int i = 0; i < nelement; i++ ) value[i] *= factor;
   }

////////////////////////////////////////////////////////////////////
// Not the best inversion - scrounged from elsewhere (C code..)

inline
double *ROW(double *value, int i) {
    return  value+(i*(i+1))/2;
    }

inline
double &Lij(double *value, int i, int j ) {
    return ROW(value,i)[j];
    }

static int chol_dec( double *value, int np ) {
   double *r1;
   double *r2;
   double sum;
   int i,j,k,kmin;
   static double small = 1.0e-10;
   static double abssmall = 1.0e-30;
   for (i = 0; i<np; i++ ) {

      /* For each row being decomposed, skip over all initial zero's in
	 the row */

      r1 = ROW(value,i);
      for ( j=0; j<i; j++ ) {
	 if( *r1++ != 0.0 ) break;
	 }
      kmin = j;
      for ( ; j<= i; j++ ) {
	r1=ROW(value,i)+kmin; r2=ROW(value,j)+kmin; sum=0.0;
	for( k=kmin; k<j; k++ ) sum -= *r1++ * *r2++;
	sum += *r1;
	if ( i==j ) {
	    if( sum < *r1 * small || sum < abssmall ) return ++i;
	    *r1 = sqrt( sum );
	    }
	else
	    *r1 = sum/ *r2;
	}
      }
   return 0;
   }

static void chol_inv( double *value, double *scratch, int np ) {
   double *r0;
   double sum;
   int i,j,k;
   for (i=np; i--; )  {
       r0 = &Lij( value, np-1, i );
       for (j=np; j-- > i; ) { scratch[j] = *r0; r0=r0-j; }
       for (j=np; j-- > i; ) {
	  sum = (i==j) ? 1.0/scratch[i] : 0.0;
	  for (k=i; ++k <= j; ) sum -= scratch[k] * Lij(value,j,k);
	  for (k=j; ++k < np; ) sum -= scratch[k] * Lij(value,k,j);
	  Lij(value,j,i) = sum/scratch[i];
	  }
       }
   }
//////////////////////////////////////////////////////////////////////

int SymMatrix::Invert() {
   badrow = chol_dec( value, size );
   if( !badrow ) chol_inv( value, scratch, size );
   return badrow ? 0 : 1;
   }

