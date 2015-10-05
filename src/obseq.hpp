#ifndef OBSEQ_HPP
#define OBSEQ_HPP
/***************************************************************************
 Copyright 2015 Crown copyright (c)
 Land Information New Zealand and the New Zealand Government.
 All rights reserved

 This program is released under the terms of the new BSD license. See 
 the LICENSE file for more information.
****************************************************************************/

#include <assert.h>
#include <iostream>

// This implementation only supports a diagonal weight matrix - weights
// are associated with each row.

#ifndef LINEQN_HPP
#include "lineqn.hpp"
#endif

#ifndef SYMMATRX_HPP
#include "symmatrx.hpp"
#endif

class ObsRow {
   friend class Obseqn;
   friend ostream & operator << ( ostream &os, ObsRow &obs );
 public:
   ObsRow( int nCols = 5, int nImplicit = 0, int checkRowNo = 1 );
   ~ObsRow();
   double &Value(){ return value;}
   double &Weight(){ return weight;}
   double &Param( int i );
   double &Implicit( int i );
   void Zero( int nImplicit = 0, int checkRowNo = 1 );
   void SumInto( LinearEquations &le );
 private:
   double value;
   double weight;
   int *cols;
   double *param;
   int ncols;
   int nextcol;
   double *implicit;
   int nimp;
   int check;
   void SetSize( int size );
   void SetImplicit( int nImplicit );
   void AddWeightedRow( double wgt, ObsRow &row );
   void SumWeighted( double wgt, LinearEquations &le );
   void SumWeightedProduct( double wgt, ObsRow &prd, LinearEquations &le );

   // Disable copy and assignment

   ObsRow( const ObsRow & );
   ObsRow & operator = ( ObsRow & );
   };


class Obseqn {
  friend class LinearEquations;
  friend ostream & operator << ( ostream &os, Obseqn &obs );
  public:
    Obseqn( int nRows = 1, int nImplicit = 0, int checkRowNo = 1 );
    ~Obseqn();
    void Zero( int nRows = -1, int nImplicit = -1, int checkRowNo = -1 );

    // NOTE: These functions are really designed for writing value into
    // arrays.  The function A MUST NOT be used for reading the value if
    // checkRowNo is set to zero.  Even when checkRowNo is set the function
    // will create a column if it is not already defined, even when reading.

    // I haven't had time to do this properly yet!!

    double &A(int i, int j);
    double &f(int i, int j);
    double &y(int i);
    double &w(int i);
    int nrow(){ return nrows; }

    // Functions return 1 on success, 0 on failure.  Only SumInto can
    // fail, and then only if it is trying to sum equations with
    // implicit parameters.
    
    int SumInto( LinearEquations &le );
    int WeightedSumInto( LinearEquations &le, SymMatrix &weight );

    // NOTE: This routine currently does not return the covariance matrix
    // of the parameters.

    int CalcValue( LinearEquations &le, leVector *prm = 0, SymMatrix *cvr = 0 );
    int CalcImplicitParams( LinearEquations &le, leVector &prm );

  private:

    ObsRow *row;
    int nrows;
    int maxrows;
    int nimplicit;
    int check;

    SymMatrix *impWgt;
    leVector scratch;
    Obseqn *impObs;

    void SetSize( int nRows );
    ObsRow &Row( int i );
    void FormImplicitObs();

    // Disable copy and assign

    Obseqn( const Obseqn & );
    Obseqn & operator = (Obseqn &);
    };

inline
double &ObsRow::Implicit( int i ) {
    assert( i > 0 && i <= nimp );
    return implicit[i-1];
    }

inline
ObsRow &Obseqn::Row( int i ) {
    assert( i > 0 && i <= nrows );
    return row[i-1];
    }

inline
double & Obseqn::A( int i, int j ) { return Row(i).Param(j); }

inline
double & Obseqn::f( int i, int j ) { return Row(i).Implicit(j); }

inline
double & Obseqn::y( int i ) { return Row(i).Value(); }

inline
double & Obseqn::w( int i ) { return Row(i).Weight(); }

#endif
