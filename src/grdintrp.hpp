#ifndef GRDINTRP_HPP
#define GRDINTRP_HPP
/***************************************************************************
 Copyright 2015 Crown copyright (c)
 Land Information New Zealand and the New Zealand Government.
 All rights reserved

 This program is released under the terms of the new BSD license. See 
 the LICENSE file for more information.
****************************************************************************/

#ifndef GRID_HPP
#include "grid.hpp"
#endif

// Base class for a linear interpolator on grid nodes (linear in grid node
// values, not in position.


class GridInterpolationPoint {
  public:
    long col, row;
    double dx[2], dy[2];  // Influence of grid x, y offsets on dx, and on dy
    };

class GridInterpolator {
  public:
    virtual void setupInterpolationPoint( double xy[2], Grid &grid ) = 0;
    int calcOffset( double xy[2], Grid &grid, double dxy[2] );
    const int pointInfluenceRange(){ return range; }
    const int nInterpolationPoints(){ return npoints; }
    GridInterpolationPoint *operator[] (int i) { return & pts[i]; }
  protected:
    GridInterpolator( int npoints, int range );
    virtual ~GridInterpolator();
    GridInterpolationPoint *pts;
  private:
    int npoints;
    int range;
    };


class BilinearInterpolator : public GridInterpolator {
  public:
     BilinearInterpolator();
     virtual void setupInterpolationPoint( double xy[2], Grid &grid );
  };

#endif
