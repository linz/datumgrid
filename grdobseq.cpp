// Module to set up the observation equations for calculating offsets at
// grid points...

#include <iostream>
#include <math.h>

using namespace std;

#include "grdobseq.hpp"
#include "obseq.hpp"
#include "progress.hpp"
#include "fixfmt.hpp"

// Grabbed from elsewhere - could be much more efficient!

static double vector_standardised_residual( int size, double *vec, double *cvr, int &rank ) {

   double *r1, *r2, *ri, *rj;
   double sum, prod;
   int i,j,k;
   static double small = 1.0e-10;

   /* Cholesky decomposition */

   rank = size;
   for (i = 0, ri=cvr; i<size; i++, ri+=i ) {
      for ( j=0, rj = cvr; j<= i; j++, rj+=j ) {
        r1=ri; r2=rj; sum=0.0;
        for( k=0; k<j; k++ ) sum -= *r1++ * *r2++;
        sum += *r1;
        if ( i==j ) {
            if( sum < small ) { (rank)--; *r1 = 1.0; } else *r1 = sqrt( sum );
            }
        else
            *r1 = sum/ *r2;
        }
      }

   /* Form L"v, and its dot product with itself (this is the number we want) */

   prod = 0.0;

   for ( i=0, r1=cvr; i<size; i++ ) {
       for (k = 0, sum=0.0; k<i; k++ ) sum -= *r1++ * vec[k];
       vec[i] = (vec[i]+sum) / *r1++;
       prod += vec[i] * vec[i];
       }

   /* Convert to standardised residual by dividing by its rank and
      square rooting */

   if( rank ) {
      prod /= rank;
      if( prod > 0 ) prod = sqrt(prod); else prod = 0.0;
      }

   else {
      prod = 1.0;
      }

   return prod;
   }

// Setup up the bandwidth definition for each matrix element based upon the
// first and last parameters that can be spanned by the observation equations.
// (This is expressed as the range of cells influenced by a point - generally
// will depend upon the interpolation method)

void setupBandwidthDefinition( Grid &grd, BLT_Def &bltdef, int ptInfluenceRange ) {
   int pir = ptInfluenceRange;
   if( pir < 1 ) pir = 1;

   for( long i = 0; i < grd.nrows(); i++ ) {
      for( long j = 0; j < grd.ncols(); j++ ) {
          long prm1 = grd.paramNo(j,i);
          if( prm1 <= 0 ) continue;
          long prm2 = -1;
          for( long ic = i - pir; ic <= i && prm2 <= 0; ic++ ) {
            for( long jr = j - pir; jr <= j; jr++ ) {
               prm2 = grd.paramNo( jr, ic );
               if( prm2 > 0 ) break;
               }
            }
          bltdef.SetNonZero( prm1, prm2 );
          bltdef.SetNonZero( prm1 + 1, prm2 );
         }
      }
   }

// Setup the distortion parameters which are those that we want to minimize
// in our solution.  These are set up for the grid element defined by the
// four grid nodes (c,r), (c+1,r), (c,r+1), (c+1,r+1).

static double distortionParam[4][8];
static double affineParam[2][8];

void setupDistortionParam( Grid &grd, double dstError ) {
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

   const double *grdspc = grd.getSpacing();
   double mult = 1000000.0/(dstError*sqrt(0.5*fabs(grdspc[0]*grdspc[1])));

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


int DistortionObseqn( Grid &grd, long c, long r, Obseqn &oe ) {
   long prmNo[4];
   prmNo[0] = grd.paramNo( c, r ); if( prmNo[0] <= 0 ) return 0;
   prmNo[1] = grd.paramNo( c+1, r ); if( prmNo[1] <= 0 ) return 0;
   prmNo[2] = grd.paramNo( c, r+1 ); if( prmNo[2] <= 0 ) return 0;
   prmNo[3] = grd.paramNo( c+1, r+1 ); if( prmNo[3] <= 0 ) return 0;
   double xy0[2];
   double xy1[2];
   grd.convert(c,r,xy0);
   grd.convert(c+1,r+1,xy1);
   // cout << "distortion " << (xy0[0]+xy1[0])/2.0 << " " << (xy0[1]+xy1[1])/2.0 
   //    << " " << prmNo[0]
   //    << " " << prmNo[1]
   //    << " " << prmNo[2]
   //    << " " << prmNo[3]
   //    << "\n";
   oe.Zero( 4, 0, 0 );
   for( int d = 0; d < 4; d++ ) {
      double *dp = & distortionParam[d][0];
      for( int nod = 0; nod < 4; nod++ ) {
         oe.A(d+1,prmNo[nod]) = *dp++;
         oe.A(d+1,prmNo[nod]+1) = *dp++;
         }
      }
   return 1;
   }


long sumDistortionConstraints( Grid &grd, LinearEquations &le, ProgressMeter &pm ) {
   long nr = grd.nrows() - 1;
   long nc = grd.ncols() - 1;
   pm.Start("Applying distortion constraints", nr );
   Obseqn oe(4);
   long nsum=0;
   for( long r = 0; r < nr; r++ ) {
      pm.Update(r);
      for( long c = 0; c < nc; c++ ) {
         if( DistortionObseqn( grd, c, r, oe )) 
         {
             oe.SumInto( le );
             nsum++;
         }
         }
      }
   pm.Finish();
   return nsum;
   }


void calcGridStandardisedResiduals( Grid &grd, LinearEquations &le, ProgressMeter &pm ) {
   long nr = grd.nrows() - 1;
   long nc = grd.ncols() - 1;
   pm.Start("Calculating distortion residuals", nr );
   Obseqn oe(4);
   for( long r = 0; r < nr; r++ ) {
      pm.Update(r);
      for( long c = 0; c < nc; c++ ) {
         double sr = 0.0;
         if( DistortionObseqn( grd, c, r, oe ) ) {
            leVector calc(4);
            SymMatrix cvr(4);
            if( oe.CalcValue(le,&calc,&cvr ) ) {
               double tmpvec[4];
               double tmpcvr[10];
               int i, j, k;
               k = 0;
               for( i = 0; i < 4; i++ ) {
                  tmpvec[i] = calc(i+1);
                  for( j = 0; j <= i; j++ ) {
                     tmpcvr[k] = cvr(i+1,j+1);
                     if( i == j ) tmpcvr[k]++;
                     k++;
                     }
                  }
               int rank;
               sr = vector_standardised_residual( 4, tmpvec, tmpcvr, rank );
               }
            grd(c,r).sr = sr;
            }
         }
      }
   pm.Finish();
   }


void writeGridDistortion( Grid &grd, ostream &os ) {
   long nr = grd.nrows() - 1;
   long nc = grd.ncols() - 1;
   double xc = grd.getSpacing()[0]/2.0;
   double yc = grd.getSpacing()[1]/2.0;
   setupDistortionParam( grd, 1.0 );
   os << "x,y,scale,rotation,distortion,shear,rlaxis,stdres\n";
   for( long r = 0; r < nr; r++ ) {
      for( long c = 0; c < nc; c++ ) {
         if( grd.paramNo( c, r ) <= 0 ||
             grd.paramNo( c+1, r ) <= 0 ||
             grd.paramNo( c, r+1 ) <= 0 ||
             grd.paramNo( c+1, r+1 ) <= 0 ) continue;
         double distortion[4] = { 0, 0, 0, 0 };
         double dilatation = 0;
         double rotation = 0;
         int np = 0;
         for( int ir = 0; ir < 2; ir ++ ) {
            for( int ic = 0; ic < 2; ic++ ) {
               GridPoint &gp = grd(c+ic, r+ir );
               int npx = np++;
               int npy = np++;
               for( int d = 0; d < 4; d++ ) {
                  distortion[d] += distortionParam[d][npx]*gp.dxy[0] +
                                   distortionParam[d][npy]*gp.dxy[1];
                  }
               dilatation += affineParam[0][npx]*gp.dxy[0] +
                             affineParam[0][npy]*gp.dxy[1];
               rotation   += affineParam[1][npx]*gp.dxy[0] +
                             affineParam[1][npy]*gp.dxy[1];
               }
            }
         double cumdist = 0.0;
         for( int d = 0; d < 4; d++ ) cumdist += distortion[d]*distortion[d];
         cumdist = sqrt(cumdist);
         double shear = hypot(distortion[0],distortion[1]);
         double rlaxis = atan2(distortion[1],distortion[0])*(90/M_PI) - 45;
         long cr[2] = {c,r};
         double xy[2];
         grd.convert( cr, xy );
         os << FixedFormat(0) << (xy[0]+xc) << "," << (xy[1]+yc) << ","
            << FixedFormat(2) << dilatation << "," << rotation << ","
            << cumdist << "," << shear << "," << rlaxis << ","
            << grd(c,r).sr << "\n";
         }

      }
   }


void ControlPointObseqn( Grid &grd, ControlPoint &cp, GridInterpolator &gi,
    Obseqn &oe ) {
    oe.Zero( 2, 0, 0 );
    gi.setupInterpolationPoint( cp.coord(), grd );
    cout << "Control point:";
    for( int ip = 0; ip < gi.nInterpolationPoints(); ip++ ) {
       GridInterpolationPoint &gip = *gi[ip];
       long prmNo = grd.paramNo( gip.col, gip.row );
       if( prmNo <= 0 ) continue;
       cout << " " << prmNo;
       oe.A(1, prmNo ) = gip.dx[0];
       oe.A(2, prmNo ) = gip.dy[0];
       prmNo++;
       oe.A(1, prmNo ) = gip.dx[1];
       oe.A(2, prmNo ) = gip.dy[1];
       }
    oe.y(1) = cp.offset()[0];
    oe.y(2) = cp.offset()[1];
    double wgt = 1.0/cp.getError();
    wgt *= wgt;
    oe.w(1) = wgt;
    oe.w(2) = wgt;
    cout <<  " weight " << wgt << "\n";
    }


long sumControlPoints( Grid &grd, ControlPointList &pts, GridInterpolator &gi,
    LinearEquations &le, ProgressMeter &pm ) {
    pm.Start("Adding control points",pts.size() );
    Obseqn oe(2);
    int nsum=0;
    for( int i = 0; i < pts.size(); i++ ) {
       pm.Update(i);
       ControlPoint &cp = *pts[i];
       if( cp.isRejected() ) continue;
       ControlPointObseqn( grd, cp, gi, oe );
       oe.SumInto( le );
       nsum++;
       }
    pm.Finish();
    return nsum;
    }



void calcControlPointResidual( Grid &grd, ControlPoint &cp, GridInterpolator &gi,
    LinearEquations &le ) {
    Obseqn oe(2,0,0);
    ControlPointObseqn( grd, cp, gi, oe );
    leVector calc(2);
    SymMatrix cvr(2);
    if( ! oe.CalcValue( le, &calc, &cvr )) { return; }
    cp.calcOffset()[0] = calc(1);
    cp.calcOffset()[1] = calc(2);
    double tmpvec[2], tmpcvr[3];
    tmpvec[0] = cp.calcOffset()[0] - cp.offset()[0];
    tmpvec[1] = cp.calcOffset()[1] - cp.offset()[1];
    cp.distanceResidual() = hypot( tmpvec[0], tmpvec[1] );
    double sign =  cp.isRejected() ? 1 : -1;
    double wgt = cp.getError();
    wgt *= wgt;
    tmpcvr[0] = wgt + sign * cvr(1,1);
    tmpcvr[1] = sign * cvr(1,2);
    tmpcvr[2] = wgt + sign * cvr(2,2);
    int rank;

    cp.stdResidual() = vector_standardised_residual( 2, tmpvec, tmpcvr, rank );
    }



void calcControlPointListResiduals( Grid &grd, ControlPointList &pts, GridInterpolator &gi,
    LinearEquations &le, ProgressMeter &pm ) {
    pm.Start("Calculating control point residuals",pts.size() );
    for( int i = 0; i < pts.size(); i++ ) {
       pm.Update(i);
       ControlPoint &cp = *pts[i];
       calcControlPointResidual( grd, cp, gi, le );
       }
    pm.Finish();
    }


int CalculateGridModel( Grid &grd, ControlPointList &pts,
                        GridInterpolator &gi, double dstError ) {
    AsciiBarMeter pm;

    // Set up the linear equations

    cout << "Setting up the linear equations" << endl;

    BLT_Def *blt = new BLT_Def( grd.paramCount() );
    setupBandwidthDefinition( grd, *blt, gi.pointInfluenceRange() );
    LinearEquations le( *blt );
    delete blt;

    //  Sum the equations

    long nsum=sumControlPoints( grd, pts, gi, le, pm );
    cout << "Summed " << nsum << " control points" << endl;
    if( nsum == 0 )
    {
        cout << "No data - failed!" << endl;
        return 0;
    }
    
    // Apply the distortion constraints

    setupDistortionParam( grd, dstError );
    nsum=sumDistortionConstraints( grd, le, pm );
    cout << "Summed " << nsum << " distortion constraints" << endl;

    // Solve the equations

    cout << "Solving the linear equations" << endl;

    BLT_Matrix::SetProgressMeter( &pm );
    int status = le.Solve();
    BLT_Matrix::SetProgressMeter( 0 );

    if( ! status ) {
       cout << "Solution failed at row " << le.BadRow() << endl;
       return 0;
       }

    // Copy the solution to the grid

    pm.Start("Applying the solution to the grid",grd.nrows() );
    for( long r = 0; r < grd.nrows(); r++ ) {
        pm.Update(r);
        for( long c = 0; c < grd.ncols(); c++ ) {
            long prmNo = grd.paramNo( c, r );
            if( prmNo <= 0 ) continue;
            grd(c,r).dxy[0] = le.Param(prmNo);
            grd(c,r).dxy[1] = le.Param(prmNo+1);
            }
        }
    pm.Finish();

    // Invert the covariance matrix

    BLT_Matrix::SetProgressMeter( &pm );
    le.Invert();
    BLT_Matrix::SetProgressMeter( 0 );


    // Update the control points with calculated values and residuals

    calcControlPointListResiduals( grd, pts, gi, le, pm );
    calcGridStandardisedResiduals( grd, le, pm );

    return 1;
    }
