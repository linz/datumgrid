#ifndef BLTMATRX_HPP
#define BLTMATRX_HPP
/***************************************************************************
 Copyright 2015 Crown copyright (c)
 Land Information New Zealand and the New Zealand Government.
 All rights reserved

 This program is released under the terms of the new BSD license. See 
 the LICENSE file for more information.
****************************************************************************/

#include <assert.h>

// Header file for banded lower triangular matrix structure.
// Notes: All matrix operations are 1 based, not 0 based.

// These routines are very strongly focussed upon solving least squares
// equations.  In fact that's all they are good for!

class ProgressMeter;

class BLT_Def {
  friend class BLT_Matrix;
     int nRows;
     int *rows;
   public:
     BLT_Def( int nrows );
     ~BLT_Def();
     void SetNonZero( int row, int col );
     int Size(){ return nRows; }
     };

class BLT_Matrix {
    public:
       // Constructors for 1) full matrix
       //                  2) bandwidth limited matrix
       //                  3) matrix with specified first column non
       //                     zero for each row
       BLT_Matrix( int rows );
       BLT_Matrix( int rows, int bandwidth );
       BLT_Matrix( BLT_Def &def );
       ~BLT_Matrix();

       double & operator()( int i, int j );

       void Zero();
       int Size(){ return nRows; }

       // Return 1 = OK, 0 = failure.  BadRow returns number of row in
       // which singularity is detected.
       int FormCholeskiDecomposition();
       // Note: Inverse doesn't include non-zero elements outside
       // original bandwidth - use with caution. It cannot be used
       // for solving equations.  This must be done from the Choleski
       // decomposition.
       int Invert();
       // Note: Solve can only be called when the Choleski decomposition
       // has been formed.
       int Solve( double *b, double *s );
       // Multiplies a vector b by the inverse choleski matrix...
       int CholInvMult( double *b, double *s );
       int BadRow(){ return badrow; }
       long NonZeroCount();  // Total number of potentially n.z. elements in matrix
       static void SetProgressMeter( ProgressMeter *newMeter ){ meter = newMeter; }

    private:

       class BLT_Row {
         friend class BLT_Matrix;
         int first;
         int last;
         double *values;

         BLT_Row();
         ~BLT_Row();
         void SetSize( int first, int last );
         void Zero();
         double & operator()( int i );
         // double operator *( BLT_Row &r2 );
         // double operator *( double *vector );
         };

       int nRows;
       int badrow;

       BLT_Row *rows;
       BLT_Row &Row( int i );
       static double small;
       static double absSmall;
       static ProgressMeter *meter;

       enum { Entering, Choleski, Inverse, Garbage } status;
    };

inline
double & BLT_Matrix::BLT_Row::operator ()( int i ) {
     assert( i >= first && i <= last );
     return values[i-first];
     }

inline
BLT_Matrix::BLT_Row &BLT_Matrix::Row( int i ) {
    assert( i > 0 && i <= nRows );
    return rows[i-1];
    }

inline
double & BLT_Matrix::operator ()( int i, int j ) {
     assert( i > 0 && i <= nRows );
     assert( j > 0 && j <= nRows );
     return ( i < j ) ? Row(j)(i) : Row(i)(j);
     }

#endif
