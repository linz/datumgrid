#ifndef GRID_HPP
#define GRID_HPP
/***************************************************************************
 Copyright 2015 Crown copyright (c)
 Land Information New Zealand and the New Zealand Government.
 All rights reserved

 This program is released under the terms of the new BSD license. See 
 the LICENSE file for more information.
****************************************************************************/

#include <iostream>
#include <cstring>
#include <string>
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
        enum GridBoundaryOption { grdFit, grdZero, grdIgnore };
       double xSpacing;
       double ySpacing;
       double xmin;
       double ymin;
       double xmax;
       double ymax;
       double xoffset;
       double yoffset;
       int ngridx;
       int ngridy;
       double xScale;
       double yScale;
       double maxPointProximity;
       GridBoundaryOption boundaryOption;
       double heightZero;
       double distortionError;
       double shearWeight;
       double scaleWeight;
       double nonConstantWeight;
       double nonLinearWeight;
       double controlNodeTolerance;
       int  pointInfluenceRange;
       int ndpCoord;
       int ndpValue;
       string xcolname;
       string ycolname;
       string dxcolname;
       string dycolname;
       string xerrname;
       string yerrname;
       string xycorrname;
       bool heightGrid;
       bool printGridParams;
       bool fillGrid;
       bool calcGridCovar;
       bool calcStdRes;
       bool fixedGrid;
       bool fixControlNodes;
       bool controlNodesOnly;
       bool pointsHaveIds;
       GridParams() :
          xSpacing(50000.0),
          ySpacing(50000.0),
          xmin(0.0),
          ymin(0.0),
          xmax(0.0),
          ymax(0.0),
          xoffset(0.0),
          yoffset(0.0),
          ngridx(0),
          ngridy(0),
          xScale(1.0),
          yScale(1.0),
          maxPointProximity(100000.0),
          controlNodeTolerance(1.0),
          boundaryOption(GridParams::grdFit),
          heightZero(0.0),
          distortionError(1.0),
          shearWeight(1.0),
          scaleWeight(1.0),
          nonConstantWeight(1.0),
          nonLinearWeight(1.0),
          pointInfluenceRange(1),
          ndpCoord(0),
          ndpValue(4),
          xcolname("x"),
          ycolname("y"),
          dxcolname("dx"),
          dycolname("dy"),
          xerrname("stdx"),
          yerrname("stdy"),
          xycorrname("corrxy"),
          heightGrid(false),
          printGridParams(false),
          fillGrid(false),
          calcGridCovar(false),
          calcStdRes(true),
          fixedGrid(false),
          fixControlNodes(false),
          controlNodesOnly(false),
          pointsHaveIds(true)
          {}
    };


class Grid;

struct GridPoint {
  GridPoint(){ dxy[0]=dxy[1]=0.0; sr=0.0; inrange=false; fixed=false; paramno=-1;}
  double dxy[2];
  double cvr[3];
  double sr;  // Standardised residual of distortion of the adjacent cell
  long  paramno;
  bool fixed;
  bool inrange;
  };


// GridRow - used internally by Grid - should be nested class?

class GridRow {
    friend class Grid;
  public:
  private:
    GridRow() : cmin(-1), cmax(-2), paramno(-1), data(0) { ;}
    ~GridRow();
    void expandRange( long nmin, long nmax );
    void crdRange( DoubleRange &range, int crd );
    long setParamNo( long paramno, GridParams::GridBoundaryOption boundaryOption, bool heightGrid );
    long nFixed();
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
    void convert( const double xy[2], long cr[2] );
    void convert( const long  c, const long r , double xy[2] );
    void convert( const long  cr[2], double xy[2] ){ convert( cr[0], cr[1], xy ); }
    void gridCoords( double xy[2], double gxy[2] );
    const double *getSpacing(){ return spacing; }
    const double *getScale(){ return scale; }
    bool isHeightGrid(){ return heightGrid; }
    char isValidPoint( long c, long r );
    char isValidPoint( long cr[2] ){ return isValidPoint( cr[0], cr[1] ); }
    GridPoint & operator() (long c, long r );
    long paramNo( long c, long r );
    long paramCount(){ return paramcount; }
    long nFixed(){ return nfixed; }
    ostream &dumpSpecTo( ostream &os, bool full=false );
    int writeSurferFiles( string &rootName );
  private:
    void writeSurferFile( ostream &os, int crd );
    double xy0[2], spacing[2], scale[2];
    long ngrd[2];
    long paramcount;
    long nfixed;
    bool heightGrid;
    GridPoint dummy;
    GridRow *rows;
    };

#endif
