#ifndef LINEQN_HPP
#define LINEQN_HPP
/***************************************************************************
 Copyright 2015 Crown copyright (c)
 Land Information New Zealand and the New Zealand Government.
 All rights reserved

 This program is released under the terms of the new BSD license. See 
 the LICENSE file for more information.
****************************************************************************/

// ******** WARNING **********
// In this implementation the linear equations are currently not inverted
// after the LinearEquations::Solve()....

#ifndef BLTMATRX_HPP
#include "bltmatrx.hpp"
#endif

class leVector {
  friend class LinearEquations;   // Yucky - but I'm in a hurry.
     double *vec;
     int size;
     int maxsize;
  public:
     leVector( int size );
     ~leVector();
     double &operator ()( int i );
     double *Vector(){ return vec; }
     void Zero( int size = 0 );
     };

class LinearEquations {
  friend class Obseqn;
  friend class ObsRow;     // Even yuckier...
  public:
    LinearEquations( BLT_Def &def );
    LinearEquations( int size );
    int Solve();
    int Invert();
    int BadRow(){ return badrow; }
    int NParam(){ return N.Size(); }
    double Param( int i ){ return param(i); }
    int Summing(){ return status == leSumming ? 1 : 0; }
    int Solved(){ return status == leSolved || status == leInverted ? 1 : 0; }
    int Inverted() { return status == leInverted; }
  private:
    BLT_Matrix N;
    leVector b;
    leVector param;
    leVector scratch;
    double ssr;
    int nobs;
    int nimplicit;
    int badrow;
    enum { leSumming, leSolved, leInverted, leGarbage } status;
    };

inline
leVector::leVector( int size ) : size(0), maxsize(0), vec(0) {
    assert( size > 0 );
    Zero( size );
    }

inline
leVector::~leVector() { delete [] vec; }

inline
double &leVector:: operator () ( int i ) {
    assert( i > 0 && i <= size );
    return vec[i-1];
    }

#endif
