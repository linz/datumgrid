#ifndef GRID_HPP
#define GRID_HPP

#include <iostream>
#include <cstring>
#include "contrlpt.hpp"

using namespace std;

class DoubleRange {
  public:
    DoubleRange() : min(0.0), max(0.0), npt(0) {}
    void expandToInclude( double v );
    double min, max;
    long npt;

    };

inline void DoubleRange::expandToInclude( double v ) {
    if( npt++ == 0 ) {
        min = max = v;
        }
    else if( v < min ) {
        min = v;
        }
    else if( v > max ) {
        max = v;
        }
    }

// Parameters required to set up a grid...

class GridParams {
    public:
       double spacing;
       double maxPointProximity;
       bool zeroOutsideProximity;
       double distortionError;
       int  pointInfluenceRange;
       int ndpCoord;
       int ndpValue;
       GridParams() :
          spacing(50000.0),
          maxPointProximity(100000.0),
          zeroOutsideProximity(false),
          distortionError(1.0),
          pointInfluenceRange(1),
          ndpCoord(0),
          ndpValue(4)
          {}
    };


class Grid;

struct GridPoint {
  GridPoint(){ dxy[0]=dxy[1]=0.0; sr=0.0; inrange=false; paramno=-1;}
  double dxy[2];
  double sr;  // Standardised residual of distortion of the adjacent cell
  long  paramno;
  bool inrange;
  };


// GridRow - used internally by Grid - should be nested class?

class GridRow {
    friend class Grid;
  public:
  private:
    GridRow() : cmin(-1), cmax(-2), paramno(-1), data(0) { ;}
    void expandRange( long nmin, long nmax );
    void crdRange( DoubleRange &range, int crd );
    long setParamNo( long paramno, bool inRangeOnly );
    void writeSurferRow( ostream &os, long nrow, int crd );

    void allocate();
    long cmin, cmax;
    long paramno;
    GridPoint *data;
  };



class Grid {
  public:
    Grid( GridParams &param, ControlPointList &pts );
    ~Grid();
    long nrows(){ return ngrd[1]; }
    long ncols(){ return ngrd[0]; }
    void convert( double xy[2], long cr[2] );
    void convert( long  c, long r , double xy[2] );
    void convert( long  cr[2], double xy[2] ){ convert( cr[0], cr[1], xy ); }
    void gridCoords( double xy[2], double gxy[2] );
    const double *getSpacing(){ return spacing; }
    char isValidPoint( long c, long r );
    GridPoint & operator() (long c, long r );
    long paramNo( long c, long r );
    long paramCount(){ return paramcount; }
    ostream &dumpSpecTo( ostream &os, bool full=false );
    int writeSurferFiles( string &rootName );
  private:
    void writeSurferFile( ostream &os, int crd );
    double xy0[2], spacing[2];
    long ngrd[2];
    long paramcount;
    GridPoint dummy;
    GridRow *rows;
    };

#endif
