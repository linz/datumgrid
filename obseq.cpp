// Simplistic code for managing a sparse observation equation, and
// summing into a banded least squares matrix.

// Note: covariance of output results currently not implemented.

#include <assert.h>

using namespace std;

#include "obseq.hpp"
#include "symmatrx.hpp"

ObsRow::ObsRow( int nCols, int nImplicit, int checkRowNo ) :
   value(0.0),
   weight(1.0),
   cols(0),
   param(0),
   ncols(0),
   nextcol(0),
   implicit(0),
   nimp(0),
   check(checkRowNo) {
   SetImplicit(nImplicit);
   SetSize(nCols);
   };

ObsRow::~ObsRow() {
   if( cols ) delete [] cols;
   if( param ) delete [] param;
   if( implicit ) delete [] implicit;
   }

void ObsRow::SetSize( int size ) {
   if( size <= 0 ) size = 8; // Very arbitrary....
   int * newCols = new int[size];
   double * newParam = new double[size];
   if( nextcol ) {
      for( int i = 0; i < nextcol; i++ ) {
         newCols[i] = cols[i];
         newParam[i] = param[i];
         }
      }
   if( cols ) delete [] cols;
   if( param ) delete [] param;
   cols = newCols;
   param = newParam;
   ncols = size;
   }

void ObsRow::Zero( int nImplicit, int checkRowNo ) {
   nextcol = 0;
   value = 0.0;
   weight = 1.0;
   SetImplicit( nImplicit );
   check = checkRowNo;
   }

void ObsRow::SetImplicit( int nImplicit ) {
   if( nimp < nImplicit ) {
      if( implicit ) delete [] implicit;
      implicit = new double[nImplicit];
      nimp = nImplicit;
      }
   for( int i = 0; i < nImplicit; i++ ) implicit[i] = 0.0;
   }

double &ObsRow::Param( int col ) {
   if( check ) {
      for( int i = 0; i < nextcol; i++ ) {
         if( cols[i] == col ) return param[i];
         }
      }
   if( nextcol >= ncols ) SetSize( ncols+ncols );
   cols[nextcol] = col;
   param[nextcol] = 0.0;
   return param[nextcol++];
   }


void ObsRow::SumInto( LinearEquations &le ) {
   SumWeighted( weight, le );
   }

void ObsRow::AddWeightedRow( double wgt, ObsRow &row ) {
   assert( !implicit ); // This isn't supposed to be used for anything much
   for( int i = 0; i < row.nextcol; i++ ) {
      Param(row.cols[i]) += wgt * row.param[i];
      }
   Value() += wgt * row.Value();
   }


void ObsRow::SumWeighted( double wgt, LinearEquations &le ) {
   for( int i = 0; i < nextcol; i++ ) {
      int row = cols[i];
      double pw = param[i] * wgt;
      le.b(row) += pw * value;
      for( int j = 0; j <= i; j++ ) {
         le.N(row,cols[j]) += pw * param[j];
         }
      }
   le.ssr += value*value*wgt;
   le.nobs++;
   }

// This adds the products r1 . w . row2 and row2 . w . row1 into the
// linear equations.  It  is only for use by Obseqn.

void ObsRow::SumWeightedProduct( double wgt, ObsRow &prd, LinearEquations &le ) {
   for( int i = 0; i < nextcol; i++ ) {
       int row = cols[i];
       double pw = param[i] * wgt;
       le.b(row) += pw * prd.value;
       for( int j = 0; j < prd.nextcol; j++ ) {
           double v = pw * prd.param[j];
           if(  prd.cols[j] == row ) v += v;
           le.N(row,prd.cols[j]) += v;
            }
        }
   double vw = value * wgt;
   for( int j = 0; j < prd.nextcol; j++ ) {
       le.b(prd.cols[j]) += vw * prd.param[j];
       }
   le.ssr += vw * prd.value;
   }


ostream & operator << ( ostream &os, ObsRow &obs ) {
   os << "V: " << obs.value << "  W: " << obs.weight;
   for( int i = 0; i < obs.nextcol; i++ ) {
      os << "  p" << obs.cols[i] << ": " << obs.param[i];
      }
   os << "\n";
   return os;
   }


//////////////////////////////////////////////////////////////////////

Obseqn::Obseqn( int nRows, int nImplicit, int checkRowNo ) :
    row(0),
    nrows(0),
    maxrows(0),
    nimplicit( nImplicit ),
    impObs(0),
    impWgt(0),
    scratch(1),
    check( checkRowNo ) {
    assert( nRows > 0 );
    SetSize( nRows );

    // Update rows if default constructor will not be valid
    if( nimplicit || !check ) Zero();
    }

Obseqn::~Obseqn() {
    if( row ) delete [] row;
    if( impObs ) delete impObs;
    if( impWgt ) delete impWgt;
    }

void Obseqn::Zero( int nRows, int nImplicit, int checkRowNo ) {
    if( nRows > maxrows ) SetSize( nRows );
    if( nRows != -1 ) nrows = nRows;
    if( nImplicit != -1 ) nimplicit = nImplicit;
    if( checkRowNo != -1 ) check = checkRowNo;
    for( int i = 0; i < nrows; i++ ) row[i].Zero( nimplicit, check );
    }

void Obseqn::SetSize( int nRows ) {
    nrows = nRows;
    if( nRows <= maxrows ) return;
    if( row ) delete [] row;
    row = new ObsRow[nRows];
    maxrows = nRows;
    }

void Obseqn::FormImplicitObs() {

    if( !impObs ) {
       impObs = new Obseqn( nimplicit );
       impWgt = new SymMatrix( nimplicit );
       }
    else {
       impObs->Zero( nimplicit, 0 );
       impWgt->Zero( nimplicit );
       }

    Obseqn &obs = *impObs;
    SymMatrix &wgt = *impWgt;

    // Form the rows used to eliminate the implicit parameters

    for( int i = 1; i <= nrows; i++ ) {
       for( int j = 1; j <= nimplicit; j++ ) {
          double wa = Row(i).Weight() * Row(i).Implicit(j);
          obs.Row(j).AddWeightedRow( wa, Row(i) );
          for( int k = 1; k <= j; k++ ) {
             wgt(j,k) += wa * Row(i).Implicit(k);
             }
          }
       }
    }

int Obseqn::SumInto( LinearEquations &le ) {

    assert( le.Summing() );

    // Handle implicit parameters
    if( nimplicit ) {

       // Set up the arrays needed for the implicit equations

       FormImplicitObs();

       // Form the weight matrix used for the elimination...

       if( !impWgt->Invert()) return 0;
       impWgt->ScaleBy(-1.0);

       impObs->WeightedSumInto( le, *impWgt );

       // Remove the dummy observations from the obs count..
       le.nobs -= nrows;
       le.nimplicit++;
       }

    // Basic summation...
    for( int i = 1; i <= nrows; i++ ) Row(i).SumInto( le );
    return 1;
    }


int Obseqn::WeightedSumInto( LinearEquations &le, SymMatrix &weight ) {
    assert( nimplicit == 0 );
    assert( weight.Size() == nrows );
    assert( le.Summing() );
    for( int i = 1; i <= nrows; i++ ) {
       for( int j = 1; j < i; j++ )
          Row(i).SumWeightedProduct( weight(i,j), Row(j), le );
       Row(i).SumWeighted( weight(i,i), le );
       }
    le.nobs += nrows;
    return 1;
    }


int Obseqn::CalcValue( LinearEquations &le, leVector *prm, SymMatrix *sym ) {
    if( !le.Solved() ) return 0;
    // Haven't supported this yet..
    assert( !nimplicit || !sym );
    leVector & impParam = prm ? *prm : scratch;
    if( nimplicit && !CalcImplicitParams( le, impParam )) return 0;
    for( int i = 1; i <= nrows; i++ ) {
       ObsRow &row = Row(i);
       double value = 0;
       for( int j = 0; j < row.nextcol; j++ ) {
           value += row.param[j]*le.Param( row.cols[j] );
           }
       if( nimplicit ) {
           for( int j = 1; j <= nimplicit; j++ ) {
               value += row.Implicit(j) * impParam(j);
               }
           }
       (*prm)(i) = value;
       }
    //
    if( !sym ) return 1;
    if( nimplicit ) return 0; // In case assertion is not being tested.
    if( !le.Invert() ) return 0;

    sym->Zero(nrows);
    for( int i = 1; i <= nrows; i++ ) {
       ObsRow &oi = Row(i);
       for( int j = 1; j <= i; j++ ) {
          ObsRow &oj = Row(j);
          double cvr = 0;
          for( int ic = 0; ic < oi.nextcol; ic++ ) {
              for( int jc = 0; jc < oj.nextcol; jc++ ) {
                 cvr += oi.param[ic]*oj.param[jc]*le.N(oi.cols[ic],oi.cols[jc]);
                 }
              }
          (*sym)(i,j) = cvr;
          }
       }

    return 1;
    }


#if 0
int Obseqn::CalcResiduals( LinearEquations &le, leVector *prm ) {
    if( !le.Solved() ) return 0;
    leVector & impParam = prm ? *prm : scratch;
    if( nimplicit && !CalcImplicitParams( le, impParam )) return 0;
    for( int i = 1; i <= nrows; i++ ) {
       ObsRow &row = Row(i);
       for( int j = 0; j < row.nextcol; j++ ) {
           row.Value() -= row.param[j]*le.Param( row.cols[j] );
           }
       if( nimplicit ) {
           for( int j = 1; j <= nimplicit; j++ ) {
               row.Value() -= row.Implicit(j) * impParam(j);
               }
           }
       }
    return 1;
    }
#endif


int Obseqn::CalcImplicitParams( LinearEquations &le, leVector &prm ) {
    assert(0);
    return 0;
#if 0
    if( !nimplicit ) return 0;
    if( !le.Solved() ) return 0;
    FormImplicitObs();
    prm.Zero(nimplicit);
    if( !impWgt->Invert() ) return 0;
    impObs->CalcResiduals( le );

    // Set up the parameters...
    for( int i = 1; i <= nimplicit; i++ ) {
       double vi = impObs->Row(i).Value();
       SymMatrix &invWgt = *impWgt;
       for( int j = 1; j <= nimplicit; j++ ) {
           prm(j) += invWgt(j,i) * vi;
           }
       }

    // May want to get the covariance at a later date.
    // NOTE: To get the covariance may require that the normal matrix in
    // lineqn.cpp be inverted in Solve()....

    return 1;
#endif
    }

ostream & operator << ( ostream &os, Obseqn &obs ) {
   for( int i=0; i < obs.nrows; i++ )
   {
       os << (i ? "       " : "Obseq: ");
       os << i+1 << " " << obs.Row(i+1);
   }
   return os;
   }
