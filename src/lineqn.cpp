// Simplistic linear equations module...

// NOTE: Currently the Normal equations are not inverted at the
// conclusion of the Solve routine.  This speeds it up, but be
// prepared for garbage if you do anything else with the matrix
// elements...

using namespace std;

#include "lineqn.hpp"


void leVector::Zero( int newSize ) {
   assert( newSize >= 0 );
   if( newSize && newSize > maxsize ) {
      if( vec ) delete [] vec;
      vec = new double[newSize];
      size = maxsize = newSize;
      }
   for( int i = 0; i < size; i++ ) vec[i] = 0;
   }

//////////////////////////////////////////////////////////////////////

LinearEquations::LinearEquations( BLT_Def &def ) :
   N(def),
   b(def.Size()),
   param(def.Size()),
   scratch(def.Size()),
   ssr(0.0),
   nobs(0),
   status( leSumming ) {
   N.Zero();
   }

LinearEquations::LinearEquations( int size ) :
   N(size),
   b(size),
   param(size),
   scratch(size),
   ssr(0.0),
   nobs(0),
   status( leSumming ) {
   N.Zero();
   b.Zero();
   }

int LinearEquations::Solve() {
   if( Solved() ) return 1;
   if( !Summing() ) return 0;
   if( ! N.FormCholeskiDecomposition() ) {
      status = leGarbage;
      badrow = N.BadRow();
      return 0;
      }
   N.Solve( b.vec, param.vec );
   // N.Invert();
   status = leSolved;
   return 1;
   }


int LinearEquations::Invert() {
   if( Inverted() ) return 1;
   if( ! Solved() ) Solve();
   if( ! Solved() ) return 0;
   if( ! N.Invert() ) {
       status = leGarbage;
       return 0;
       }
   status = leInverted;
   return 1;
   }


