#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

#include "fixfmt.hpp"
#include "grid.hpp"

static string missingValue("1.70141e+038");

void GridRow::expandRange( long nmin, long nmax ) {
    if( cmax < cmin ) {
       cmax = cmin = nmin;
       }
    else if( nmin < cmin ) cmin = nmin;
    else if( nmin > cmax ) cmax = nmin;
    if( nmax < cmin ) cmin = nmax;
    else if( nmax > cmax ) cmax = nmax;
    }

void GridRow::allocate() {
    if( cmax >= cmin ) {
       data = new GridPoint[cmax-cmin+1];
       }
    }

void GridRow::crdRange( DoubleRange &range, int crd ) {
    for( long i = 0; i <= cmax-cmin; i++ ) {
       range.expandToInclude( data[i].dxy[crd] );
       }
    }

long GridRow::setParamNo( long paramno, bool inRangeOnly )
{
    this->paramno = paramno;
    if( data )
    {
        for( int i=cmin; i <= cmax; i++ )
        {
            if( inRangeOnly && ! data[i].inrange )
            {
                data[i].paramno=0;
            }
            else
            {
                data[i].paramno=paramno;
                paramno += 2;
            }
        }
    }
    return paramno;
}

void GridRow::writeSurferRow( ostream &os, long nrow, int crd ) {
    for( long i = 0; i < nrow; i++ ) {
        if( i ) {
            os << (i % 8 ? ' ' : '\n');
            }
        if( i >= cmin && i <= cmax ) {
            os << FixedFormat(3) << data[i-cmin].dxy[crd];
            }
        else {
            os << missingValue;
            }
        }
    os << '\n';
    }


Grid::Grid( GridParams &param, ControlPointList &pts ) {
    rows = 0;

    spacing[0] = spacing[1] = param.spacing;

    // Work out the number of grid cells required either side of each control
    // point.

    int ptInf = (param.pointInfluenceRange+1)/2;
    double proximity = ptInf*param.spacing*1.1;
    if( proximity <  param.maxPointProximity ) proximity = param.maxPointProximity;

    // Determine the extents covererd by the control points..

    double xmin, ymin, xmax, ymax;
    int i;

    for( i = 0; i < pts.size(); i++ ) {
       const double *xy = pts[i]->coord();
       if( i == 0 ) {
          xmin = xmax = xy[0];
          ymin = ymax = xy[1];
          }
       else {
          if( xy[0] < xmin ) xmin = xy[0]; else if( xy[0] > xmax ) xmax = xy[0];
          if( xy[1] < ymin ) ymin = xy[1]; else if( xy[1] > ymax ) ymax = xy[1];
          }
       }

    // Work out the grid location and extents...  Grid points are organised in
    // the same way as words on a page.

    xmin -= proximity; xmax += proximity;
    ymin -= proximity; ymax += proximity;

    double x0, y0;
    long ngx, ngy;

    x0 = spacing[0] * floor(xmin/spacing[0]);
    y0 = spacing[1] * ceil(ymax/spacing[1]);

    ngx = (long) ceil( (xmax - x0)/spacing[0]) + 1;
    ngy = (long) ceil( (y0 - ymin)/spacing[1]) + 1;

    // Allocate the grid rows..

    rows = new GridRow[ngy];

    // Define the range of require cells for each row.

    for( i = 0; i < pts.size(); i++ ) {
       const double *xy = pts[i]->coord();
       long rmin, rmax, cmin, cmax;
       rmin = (long) floor( (y0 - xy[1] - proximity)/spacing[1] );
       rmax = (long)  ceil( (y0 - xy[1] + proximity)/spacing[1] );
       cmin = (long) floor( (xy[0] - x0 - proximity)/spacing[0] );
       cmax = (long)  ceil( (xy[0] - x0 + proximity)/spacing[0] );

       for( ; rmin <= rmax; rmin++ ) rows[rmin].expandRange( cmin, cmax );
       }

    // Now allocate space for the grid row and count the parameters in
    // the adjustment.

    paramcount = 0;
    long r;
    for( r = 0; r < ngy; r++ ) {
       GridRow &row = rows[r];
       row.allocate();
       paramcount=row.setParamNo(paramcount,param.zeroOutsideProximity);
       }

    //  Set up the parameters of the grid..

    spacing[1] = -spacing[1];
    xy0[0] = x0;
    xy0[1] = y0;
    ngrd[0] = ngx;
    ngrd[1] = ngy;

    // Define which points are within range of a control point

    for( i = 0; i < pts.size(); i++ ) {
       if( pts[i]->isRejected() ) continue;

       const double *xy = pts[i]->coord();
       long rmin, rmax, cmin, cmax;
       rmin = (long) floor( (y0 - xy[1] - proximity)/(-spacing[1]) );
       rmax = (long)  ceil( (y0 - xy[1] + proximity)/(-spacing[1]) );
       cmin = (long) floor( (xy[0] - x0 - proximity)/spacing[0] );
       cmax = (long)  ceil( (xy[0] - x0 + proximity)/spacing[0] );

       double p2=proximity*proximity;
       for( ; rmin <= rmax; rmin++ ) 
           for( long c=cmin; c <= cmax; c++ )
           {
               GridPoint &gp=(*this)(c,rmin);
               if( gp.inrange ) continue;
               double gxy[2];
               convert(c,rmin,gxy);
               double distance=(gxy[0]-xy[0])*(gxy[0]-xy[0])+
                   (gxy[1]-xy[1])*(gxy[1]-xy[1]);
               if( distance < p2 )
               { 
                   gp.inrange=true;
               }
           }
       }
    }


Grid::~Grid() {
    delete [] rows;
    }


void Grid::convert( double xy[2], long cr[2] ) {
    cr[0] = (long) floor((xy[0] - xy0[0])/spacing[0]);
    cr[1] = (long) floor((xy[1] - xy0[1])/spacing[1]);
    }


void Grid::convert( long c, long r, double xy[2] ){
    xy[0] = xy0[0] + c * spacing[0];
    xy[1] = xy0[1] + r * spacing[1];
    }


void Grid::gridCoords( double xy[2], double gxy[2] ) {
    gxy[0] = (xy[0] - xy0[0])/spacing[0];
    gxy[1] = (xy[1] - xy0[1])/spacing[1];
    }


char Grid::isValidPoint( long c, long r ) {
    if( r < 0 || r >= ngrd[1] ) return 0;
    if( c < rows[r].cmin || c > rows[r].cmax ) return 0;
    return 1;
    }

GridPoint & Grid::operator () ( long c, long r ) {
    if( isValidPoint(c,r) ) {
       return rows[r].data[c-rows[r].cmin];
       }
    else {
       return dummy;
       }
    }

long Grid::paramNo( long c, long r ) {
    if( isValidPoint(c,r) ) {
       (*this)(c,r).paramno;
       }
    else {
       return -1;
       }
    }

ostream & Grid::dumpSpecTo( ostream &os, bool full ) {
    os << "Grid definition\n";
    os << "Size: " << nrows() << " rows by " << ncols() << " cols\n";
    os << "Top left point " << xy0[0] << ", " << xy0[1] << "\n";
    os << "Spacing " << spacing[0] << ", " << spacing[1] << "\n";
    os << "Parameters " << paramcount << "\n";
    os << "Row definitions: \n";
    if( full )
    {
        for( long i = 0; i < nrows(); i++ ) {
            os << " row " << i << ": " << rows[i].cmin << " " << rows[i].cmax
                 << " " << rows[i].paramno << "\n";
            }
    }
    return os;
    }


void Grid::writeSurferFile( ostream &os, int crd ) {
   os << "DSAA\n";
   os << ngrd[0] << " " << ngrd[1] << "\n";
   os << FixedFormat(3) << xy0[0]  << " " << (xy0[0] + (ngrd[0]-1)*spacing[0]) << '\n';
   os << FixedFormat(3) << (xy0[1] + (ngrd[1]-1)*spacing[1]) << " " << xy0[1] << '\n';
   DoubleRange zrng;
   long i;
   for( i = 0; i < ngrd[1]; i++ ) {
      rows[i].crdRange( zrng, crd );
      }
   os << FixedFormat(3) << zrng.min << " " << zrng.max << "\n";
   for( i = ngrd[1]; i-- > 0; ) {
      rows[i].writeSurferRow( os, ngrd[0], crd );
      }
   }

int Grid::writeSurferFiles( string &rootName ) {
    string fileName;
    {
      fileName = rootName + "_E.grd";
      ofstream of(fileName.c_str());
      if( of.bad() ) return 0;
      writeSurferFile( of, 0 );
      of.close();
      }
    {
      fileName = rootName + "_N.grd";
      ofstream of(fileName.c_str());
      if( of.bad() ) return 0;
      writeSurferFile( of, 1 );
      of.close();
      }
    return 0;
    }
