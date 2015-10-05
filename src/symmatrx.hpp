#ifndef SYMMATRX_HPP
#define SYMMATRX_HPP
/***************************************************************************
 Copyright 2015 Crown copyright (c)
 Land Information New Zealand and the New Zealand Government.
 All rights reserved

 This program is released under the terms of the new BSD license. See 
 the LICENSE file for more information.
****************************************************************************/

#include <assert.h>

class SymMatrix {
     double *scratch;
     double *value;
     int size;
     int nelement;
     int maxsize;
     int badrow;
     void SetSize( int size );
     SymMatrix( const SymMatrix & );
     void operator = ( SymMatrix & );
  public:
     SymMatrix( int size );
     ~SymMatrix();
     void Zero( int size = 0 );
     int Invert();  // Zero = error;
     int BadRow(){ return badrow; } // Row causing failure...
     void ScaleBy( double value );
     int Size(){ return size; }
     double &operator()( int i, int j );
     };

inline
double &SymMatrix::operator()( int i, int j ) {
     assert( i > 0 && i <= size );
     assert( j > 0 && j <= size );
     int element;
     if( i < j ) {
        element = (j*(j-1))/2 + i - 1 ;
        }
     else {
        element = (i*(i-1))/2 + j - 1;
        }
     return value[element];
     }

#endif

