/***************************************************************************
 Copyright 2015 Crown copyright (c)
 Land Information New Zealand and the New Zealand Government.
 All rights reserved

 This program is released under the terms of the new BSD license. See 
 the LICENSE file for more information.
****************************************************************************/
#include <math.h>

using namespace std;

#include "grdintrp.hpp"

GridInterpolator::GridInterpolator(  int npoints, int range ) :
   npoints(npoints), range(range) {
   pts = new GridInterpolationPoint[npoints];
   }

GridInterpolator::~GridInterpolator() {
   delete [] pts;
   }

int GridInterpolator::calcOffset( double xy[2], Grid &grid, double dxy[2] ) {
   setupInterpolationPoint( xy, grid );
   int ok = 1;
   GridInterpolationPoint *pt = pts;
   double dx = 0, dy = 0;
   for( int i = 0; i < npoints; i++, pt++ ) {
      if( ! grid.isValidPoint( pt->col, pt->row ) ) { ok = 0; break; }
      double *gdxy = grid(pt->col, pt->row).dxy;
      dx += gdxy[0] * pt->dx[0] + gdxy[1] * pt->dx[1];
      dy += gdxy[0] * pt->dy[0] + gdxy[1] * pt->dy[1];
      }
   if( ok ) {
      dxy[0] = dx; dxy[1] = dy;
      }
   else {
      dxy[0] = dxy[1] = 0.0;
      }
   return ok;
   }


BilinearInterpolator::BilinearInterpolator() :
   GridInterpolator( 4, 1 ) {}

void BilinearInterpolator::setupInterpolationPoint( double xy[2], Grid &grid ) {
   long cr[2];
   double gxy[2];
   int i, j;
   grid.gridCoords( xy, gxy );
   GridInterpolationPoint *p = pts;
   for( i = 0; i < 2; i++ ) { cr[i] = (long) floor(gxy[i]); gxy[i] -= cr[i];}
   for( i = 0; i < 2; i++ ) {
     gxy[0] = 1 - gxy[0];
     for( j = 0; j < 2; j++ ) {
        gxy[1] = 1 - gxy[1];
        p->col = cr[0] + i;
        p->row = cr[1] + j;
        p->dx[0] = p->dy[1] = gxy[0] * gxy[1];
        p->dx[1] = p->dy[0] = 0.0;
        p++;
        }
     }
   }
