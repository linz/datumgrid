// Routine to

#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <cmath>

using namespace std;

#include "fixfmt.hpp"

const double dummyValue = 1.0e10;
const string missingValue("1.70141e+038");

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

class GridFileParam {
  public:
    GridFileParam( string name ) : fname(name), fs(name.c_str()){};
    string fname;
    long loc;
    ofstream fs;
    DoubleRange rng;
  };


static double distortionParam[4][8];
static double affineParam[2][8];

void setupDistortionParam( double dx, double dy ) {
   double dp[4][8] = {
       { -1, -1, 1, -1, -1, 1, 1, 1 },  // gamma1
       { 1, -1, 1, 1, -1, -1, -1, 1 },  // gamma2
       { 0, -1, 0, 1, 0, 1, 0, -1 },    // d2vdxdy
       { -1, 0, 1, 0, 1, 0, -1, 0 }     // d2udxdy
       };
   double ap[2][8] = {
       { -1, 1, 1, 1, -1, -1, 1, -1 },  // scale
       { -1, -1, -1, 1, 1, -1, 1, 1 }   // rotation
       };

   const double grdspc[2] = { dx, dy };

   double mult = 1000000.0/(sqrt(0.5*fabs(grdspc[0]*grdspc[1])));

   for( int i = 0; i < 4; i++ ) {
      double mult1 = 0.0;
      for( int j = 0; j < 8; j++ ) mult1 += dp[i][j]*dp[i][j];
      mult1 = mult / sqrt(mult1);
      for( int j = 0; j < 8; j++ ) distortionParam[i][j] = dp[i][j] * mult1;
      }

   mult = 250000.0/(sqrt(fabs(grdspc[0]*grdspc[1])));
   for( int i = 0; i < 2; i++ ) {
      for( int j = 0; j < 8; j++ ) affineParam[i][j] = ap[i][j] * mult;
      }
   }



void calcDeformation( double *r0, double *r1, double *result ) {
    double distortion[4] = { 0, 0, 0, 0 };
    double &dilatation = result[0];
    double &rotation = result[1];
    int np = 0;
    dilatation = rotation = 0.0;
    for( int ir = 0; ir < 2; ir ++ ) {
       double *r = (ir==0) ? r0 : r1;
       for( int ic = 0; ic < 4; ic++, np++ ) {
          for( int d = 0; d < 4; d++ ) {
             distortion[d] += distortionParam[d][np]*r[ic];
             }
          dilatation += affineParam[0][np]*r[ic];
          rotation   += affineParam[1][np]*r[ic];
          }
       }
    double cumdist = 0.0;
    for( int d = 0; d < 4; d++ ) cumdist += distortion[d]*distortion[d];
    result[3] = sqrt(cumdist);
    result[2] = hypot(distortion[0],distortion[1]);
    }

int main( int argc, char *argv[] ) {
    if( argc != 2 ) {
       cout << "Syntax: griddef root_name" << endl;
       return 0;
       }

    string rootname(argv[1]);

    string filename = rootname + "_E.grd";
    ifstream fe(filename.c_str());
    filename = rootname + "_N.grd";
    ifstream fn(filename.c_str());
    if( fe.bad() || fn.bad() ) {
       cout << "Cannot open grid files " << rootname << "_E.grd and "
            << rootname << "_N.grd" << endl;
       return 0;
       }

    string dsaa;
    long nx, ny;
    double xmin, xmax, ymin, ymax, emin, emax, nmin, nmax;
    fe >> dsaa >> nx >> ny >> xmin >> xmax >> ymin >> ymax >> emin >> emax;
    if( fe.bad() || dsaa != "DSAA" ) {
        cout << "Invalid east grid file header in " << rootname << "_E.grd" << endl;
        return 0;
        }
    long nx1, ny1;
    fn >> dsaa >> nx1 >> ny1 >> xmin >> xmax >> ymin >> ymax >> nmin >> nmax;

    double dx = (xmax - xmin)/(nx-1);
    double dy = (ymax - ymin)/(ny-1);
    xmin += dx/2; xmax -= dx/2;
    ymin += dy/2; ymax -= dy/2;

    setupDistortionParam( dx, dy );

    double *row0 = new double[nx*2];
    double *row1 = new double[nx*2];

    GridFileParam gdil( rootname + "_dil.grd" );
    GridFileParam grot( rootname + "_rot.grd" );
    GridFileParam gshr( rootname + "_shr.grd" );
    GridFileParam gdst( rootname + "_dst.grd" );

    const int ngf = 4;
    GridFileParam *gf[ngf] = { &gdil, &grot, &gshr, &gdst  };
    double gval[ngf];

    int f;
    string spacer(' ',50);
    for( f = 0; f < ngf; f++ ) {
       GridFileParam *g = gf[f];
       if( g->fs.bad() ) {
           cout << "Cannot create output file " << g->fname << endl;
           return 0;
           }
       g->fs << "DSAA\n"
             << (nx-1) << " " << (ny-1) << "\n"
             << FixedFormat(3)
             << xmin << " " << xmax << "\n"
             << ymin << " " << ymax << "\n";
       g->loc = g->fs.tellp();
       g->fs << spacer << "\n";
       }
    long r, c, c0;
    for( r = 0; r < ny; r++ ) {
       double *rtmp;
       rtmp = row0; row0 = row1; row1 = rtmp;
       for( c = 0, rtmp = row0; c < nx; c++ ) {
          fe >> (*rtmp++); fn >> (*rtmp++);
          }
       if( fe.bad() || fn.bad() ) {
          cout << "Error reading grids at row " << (r+1) <<endl;
          return 0;
          }
       if( r == 0 ) continue;
       for( c = 0; c < nx-1; c++ ) {
           c0 = c * 2;
           if( c ) {
              for( f = 0; f < ngf; f++ ) gf[f]->fs << (c % 8 ? ' ' : '\n');
              }
           if( row0[c0] > dummyValue || row0[c0+1] > dummyValue ||
               row0[c0+2] > dummyValue || row0[c0+3] > dummyValue ||
               row1[c0] > dummyValue || row1[c0+1] > dummyValue ||
               row1[c0+2] > dummyValue || row1[c0+3] > dummyValue ) {
              for( f = 0; f < ngf; f++ ) { gf[f]->fs << missingValue; }
               continue;
               }
           calcDeformation( row0+c0, row1+c0, gval );
           for( f = 0; f < ngf; f++ ) {
              gf[f]->fs << gval[f];
              gf[f]->rng.expandToInclude( gval[f] );
              }
           }
       for( f = 0; f < ngf; f++ ) gf[f]->fs << "\n";
       }
    for( f = 0; f < ngf; f++ ) {
       gf[f]->fs.seekp( gf[f]->loc );
       gf[f]->fs << gf[f]->rng.min << ' ' << gf[f]->rng.max;
       }
    return 0;
    }




