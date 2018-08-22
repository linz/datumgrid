/***************************************************************************
 Copyright 2015 Crown copyright (c)
 Land Information New Zealand and the New Zealand Government.
 All rights reserved

 This program is released under the terms of the new BSD license. See 
 the LICENSE file for more information.
****************************************************************************/

// Module to set up the observation equations for calculating offsets at
// grid points...

#include <iostream>
#define _USE_MATH_DEFINES
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
            long jmax=ic == i ? j : j+pir;
            for( long jr = j - pir; jr <= jmax; jr++ ) {
               prm2 = grd.paramNo( jr, ic );
               if( prm2 > 0 ) break;
               }
            }
          bltdef.SetNonZero( prm1, prm2 );
          if( ! grd.isHeightGrid() ) bltdef.SetNonZero( prm1 + 1, prm2 );
         }
      }
   }

// Setup the distortion parameters which are those that we want to minimize
// in our solution.  These are set up for the grid element defined by the
// four grid nodes (c,r), (c+1,r), (c,r+1), (c+1,r+1).

static double distortionParam[5][8];
static double affineParam[2][8];
static int nDistortionParam;

void setupDistortionParam( Grid &grd, GridParams &param ) {
   // 1000000.0 converts to ppm
   double dstError=param.distortionError;
   double ddx=1000000.0/(grd.getScale()[0]*grd.getSpacing()[0]);
   double ddy=1000000.0/(grd.getScale()[1]*grd.getSpacing()[1]);
   double ddxy=sqrt(ddx*ddx+ddy*ddy);
   double dudx[8]={-1,0,-1,0,1,0,1,0};
   double dudy[8]={-1,0,1,0,-1,0,1,0};
   double dvdx[8]={0,-1,0,-1,0,1,0,1};
   double dvdy[8]={0,-1,0,1,0,-1,0,1};
   double d2udxdy[8]={-1,0,1,0,1,0,-1,0};
   double d2vdxdy[8]={0,-1,0,1,0,1,0,-1};
   int nd=0;
   for( int i=0; i<8; i++ )
   {
       nd=0;
       dudx[i] *= 0.5*ddx;
       dvdx[i] *= 0.5*ddx;
       dudy[i] *= 0.5*ddy;
       dvdy[i] *= 0.5*ddy;
       if( param.shearWeight > 0.0 )
       {
           double obsweight=sqrt(param.shearWeight)/(dstError*2.0);
           distortionParam[nd][i]=(dudx[i]-dvdy[i])*obsweight ; // Shear1
           distortionParam[nd+1][i]=(dudy[i]+dvdx[i])*obsweight ; // Shear2
           nd += 2;
       }
       if( param.nonLinearWeight > 0.0 )
       {
           double obsweight=sqrt(param.nonLinearWeight)/(dstError*4.0);
           distortionParam[nd][i]=d2udxdy[i]*ddxy*obsweight; // Non linear 1
           distortionParam[nd+1][i]=d2vdxdy[i]*ddxy*obsweight; // Non linear 1
           nd += 2;
       }
       if( param.scaleWeight > 0.0 )
       {
           double obsweight=sqrt(param.scaleWeight)/(dstError*2.0);
           distortionParam[nd][i]=(dudx[i]+dvdy[i])*obsweight; // Linear scale
           nd++;
       }
       affineParam[0][i]=(dudx[i]+dvdy[i])/2.0; // Linear scale
       affineParam[1][i]=(dudy[i]-dvdx[i])/2.0; // Rotation
   }
   nDistortionParam=nd;
   }

void setupHeightDistortionParam( Grid &grd, GridParams &param ) {
   // 1000000.0 converts to ppm
   double dstError=param.distortionError;
   double ddx=1000000.0/(grd.getScale()[0]*grd.getSpacing()[0]);
   double ddy=1000000.0/(grd.getScale()[1]*grd.getSpacing()[1]);
   double ddxy=sqrt(ddx*ddx+ddy*ddy);
   double c1[8]={-0.75,0,0.25,0,0.25,0,0.25,0};
   double c2[8]={0.25,0,-0.75,0,0.25,0,0.25,0};
   double c3[8]={0.25,0,0.25,0,-0.75,0,0.25,0};
   double c4[8]={0.25,0,0.25,0,0.25,0,-0.75,0};
   double t1[8]={-1,0,1,0,-1,0,1,0};
   double cobsweight=sqrt(param.nonConstantWeight)*ddxy/(dstError*2.0);
   double tobsweight=sqrt(param.nonLinearWeight)*ddxy/(dstError*2.0);
   int nd=0;
   for( int i=0; i<8; i++ )
   {
       nd=0;
       if( param.nonConstantWeight > 0.0 )
       {
           distortionParam[nd][i]=c1[i]*cobsweight ; 
           distortionParam[nd+1][i]=c2[i]*cobsweight ; 
           distortionParam[nd+2][i]=c3[i]*cobsweight ; 
           distortionParam[nd+3][i]=c4[i]*cobsweight ; 
           nd += 4;
       }
       if( param.nonLinearWeight > 0.0 )
       {
           distortionParam[nd][i]=t1[i]*tobsweight; // Non linear 1
           nd += 1;
       }
   }
   nDistortionParam=nd;
   }


int DistortionObseqn( Grid &grd, long c, long r, Obseqn &oe ) {
   GridPoint *gp[4]={&grd(c,r),&grd(c+1,r),&grd(c,r+1),&grd(c+1,r+1)};
   long prmNo[4];
   for( int i=0; i<4; i++ )
   {
       prmNo[i]=gp[i]->paramno; 
       if(prmNo[i] < 0 ) return 0;
   }
#ifdef DEBUG_GRDOBSEQ
   double xy0[2];
   double xy1[2];
   grd.convert(c,r,xy0);
   grd.convert(c+1,r+1,xy1);
   cout << "distortion " << (xy0[0]+xy1[0])/2.0 << " " << (xy0[1]+xy1[1])/2.0 
       << " " << prmNo[0]
       << " " << prmNo[1]
       << " " << prmNo[2]
       << " " << prmNo[3]
       << "\n";
#endif
   int valid=0;
   bool heightGrid=grd.isHeightGrid();
   oe.Zero( nDistortionParam, 0, 0 );
   for( int d = 0; d < nDistortionParam; d++ ) {
      double *dp = & distortionParam[d][0];
      for( int nod = 0; nod < 4; nod++ ) {
         bool sumprm=prmNo[nod] != 0;
         if( sumprm ) valid=1;
         double *dxy=gp[nod]->dxy;
         if( sumprm ) oe.A(d+1,prmNo[nod]) = *dp;
         oe.y(d+1) -= (*dp)*dxy[0];
         dp++;
         if( ! heightGrid ) 
         {
             oe.y(d+1) -= (*dp)*dxy[1];
             if( sumprm ) oe.A(d+1,prmNo[nod]+1) = *dp;
         }
         dp++;
         }
      }
   return valid;
   }


long sumDistortionConstraints( Grid &grd, LinearEquations &le, ProgressMeter &pm ) {
   long nr = grd.nrows() - 1;
   long nc = grd.ncols() - 1;
   pm.Start("Applying distortion constraints", nr );
   Obseqn oe(4);
   long nsum=0;
   for( long r = 0; r < nr; r++ ) {
#ifndef DEBUG_GRDOBSEQ
      pm.Update(r);
#endif
      for( long c = 0; c < nc; c++ ) {
         if( DistortionObseqn( grd, c, r, oe )) 
         {
#ifdef DEBUG_GRDOBSEQ
             cout << oe;
#endif
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


void writeGridDistortion( Grid &grd, GridParams &param , ostream &os ) {
   if( grd.isHeightGrid() ) return;
   long nr = grd.nrows() - 1;
   long nc = grd.ncols() - 1;
   double xc = grd.getSpacing()[0]/2.0;
   double yc = grd.getSpacing()[1]/2.0;
   double olderror=param.distortionError;
   param.distortionError=1.0;
   setupDistortionParam( grd, param );
   param.distortionError=olderror;
   os << param.xcolname << "," << param.ycolname 
       << ",scale,rotation,distortion,shear,rlaxis,stdres\n";
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
         double xy[2];
         grd.convert( c, r, xy );
         os << FixedFormat(param.ndpCoord) << (xy[0]+xc) << "," << (xy[1]+yc) << ","
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
#ifdef DEBUG_GRDOBSEQ
    cout << "cpt " << cp.getId() << ":";
#endif
    for( int ip = 0; ip < gi.nInterpolationPoints(); ip++ ) {
       GridInterpolationPoint &gip = *gi[ip];
       long prmNo = grd.paramNo( gip.col, gip.row );
       if( prmNo <= 0 ) continue;
#ifdef DEBUG_GRDOBSEQ
       cout << " " << prmNo 
           << " " << gip.dx[0] 
           << " " << gip.dy[0] 
           << " " << gip.dx[1] 
           << " " << gip.dy[1];
#endif
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
#ifdef DEBUG_GRDOBSEQ
     cout <<  " weight " << wgt << "\n";
#endif
    }

void HeightControlPointObseqn( Grid &grd, ControlPoint &cp, GridInterpolator &gi,
    Obseqn &oe ) {
    oe.Zero( 1, 0, 0 );
    gi.setupInterpolationPoint( cp.coord(), grd );
#ifdef DEBUG_GRDOBSEQ
    cout << "cpt " << cp.getId() << ":";
#endif
    for( int ip = 0; ip < gi.nInterpolationPoints(); ip++ ) {
       GridInterpolationPoint &gip = *gi[ip];
       long prmNo = grd.paramNo( gip.col, gip.row );
       if( prmNo <= 0 ) continue;
#ifdef DEBUG_GRDOBSEQ
       cout << " " << prmNo 
           << " " << gip.dx[0];
#endif
       oe.A(1, prmNo ) = gip.dx[0];
       }
    oe.y(1) = cp.offset()[0];
    double wgt = 1.0/cp.getError();
    wgt *= wgt;
    oe.w(1) = wgt;
#ifdef DEBUG_GRDOBSEQ
     cout <<  " weight " << wgt << "\n";
#endif
    }


long sumControlPoints( Grid &grd, ControlPointList &pts, GridInterpolator &gi,
    LinearEquations &le, ProgressMeter &pm ) {
    pm.Start("Adding control points",pts.size() );
    Obseqn oe(2);
    int nsum=0;
    bool heightGrid = grd.isHeightGrid();
    for( int i = 0; i < pts.size(); i++ ) {
#ifndef DEBUG_GRDOBSEQ
       pm.Update(i);
#endif
       ControlPoint &cp = *pts[i];
       if( cp.isRejected() ) continue;
       if( cp.isNode() ) continue;
       if( heightGrid ) HeightControlPointObseqn( grd, cp, gi, oe );
       else ControlPointObseqn( grd, cp, gi, oe );
#ifdef DEBUG_GRDOBSEQ
       cout << oe;
#endif
       oe.SumInto( le );
       nsum++;
       }
    pm.Finish();
    return nsum;
    }


void calcControlPointResidual( Grid &grd, ControlPoint &cp, GridInterpolator &gi,
    LinearEquations &le, bool calcStdRes ) {
    Obseqn oe(2,0,0);
    ControlPointObseqn( grd, cp, gi, oe );
    leVector calc(2);
    SymMatrix cvr(2);
    SymMatrix *ptrCvr = calcStdRes ? &cvr : 0;
    if( ! oe.CalcValue( le, &calc, ptrCvr )) { return; }
    cp.calcOffset()[0] = calc(1);
    cp.calcOffset()[1] = calc(2);
    double tmpvec[2], tmpcvr[3];
    tmpvec[0] = cp.calcOffset()[0] - cp.offset()[0];
    tmpvec[1] = cp.calcOffset()[1] - cp.offset()[1];
    cp.distanceResidual() = hypot( tmpvec[0], tmpvec[1] );
    if( calcStdRes )
    {
        double sign =  cp.isRejected() ? 1 : -1;
        double wgt = cp.getError();
        wgt *= wgt;
        tmpcvr[0] = wgt + sign * cvr(1,1);
        tmpcvr[1] = sign * cvr(1,2);
        tmpcvr[2] = wgt + sign * cvr(2,2);
        int rank;
        cp.stdResidual() = vector_standardised_residual( 2, tmpvec, tmpcvr, rank );
        double proberr=0.0;
        if( cp.stdResidual() > 0.0 )
        {
            // ratio=obserr/residualerror
            double ratio=cp.getError()/(cp.distanceResidual()/cp.stdResidual());
            proberr=cp.distanceResidual()*(ratio*ratio);
        }
        cp.probError()=proberr;
    }
    else
    {
        cp.stdResidual()=0.0;
    }
    }

void calcHeightControlPointResidual( Grid &grd, ControlPoint &cp, GridInterpolator &gi,
    LinearEquations &le, bool calcStdRes ) {
    Obseqn oe(1,0,0);
    HeightControlPointObseqn( grd, cp, gi, oe );
    leVector calc(1);
    SymMatrix cvr(1);
    SymMatrix *ptrCvr = calcStdRes ? &cvr : 0;
    if( ! oe.CalcValue( le, &calc, ptrCvr )) { return; }
    cp.calcOffset()[0] = calc(1);
    double tmpvec[1], tmpcvr[1];
    tmpvec[0] = cp.calcOffset()[0] - cp.offset()[0];
    cp.distanceResidual() = fabs( tmpvec[0] );
    if( calcStdRes )
    {
        double sign =  cp.isRejected() ? 1 : -1;
        double wgt = cp.getError();
        wgt *= wgt;
        wgt = wgt + sign * cvr(1,1);
        if( wgt > 0.0 ) wgt=sqrt(wgt);
        cp.stdResidual() = wgt > 0 ? cp.distanceResidual()/wgt : 1.0;
        cp.probError()=wgt > 0 ? cp.distanceResidual()*(cp.getError()/wgt) : 0.0;
    }
    else
    {
        cp.stdResidual()=0.0;
    }
    }



void calcControlPointListResiduals( Grid &grd, ControlPointList &pts, GridInterpolator &gi,
    LinearEquations &le, bool calcStdRes, ProgressMeter &pm ) {
    pm.Start("Calculating control point residuals",pts.size() );
    bool heightGrid = grd.isHeightGrid();
    for( int i = 0; i < pts.size(); i++ ) {
       pm.Update(i);
       ControlPoint &cp = *pts[i];
       if( cp.isNode() ) continue;
       if( heightGrid ) calcHeightControlPointResidual( grd, cp, gi, le, calcStdRes );
       else calcControlPointResidual( grd, cp, gi, le, calcStdRes );
       }
    pm.Finish();
    }


int CalculateGridModel( Grid &grd, ControlPointList &pts,
                        GridInterpolator &gi, GridParams &param ) {
    AsciiBarMeter pm;
    bool heightGrid=grd.isHeightGrid();

    // Set up the linear equations

    cout << "Setting up the linear equations" << endl;

    BLT_Def *blt = new BLT_Def( grd.paramCount() );
    setupBandwidthDefinition( grd, *blt, gi.pointInfluenceRange() );
    LinearEquations le( *blt );
    delete blt;

    //  Sum the equations

    long nsum=sumControlPoints( grd, pts, gi, le, pm );
    cout << "Summed " << nsum << " control points" << endl;
    if( nsum == 0 && grd.nFixed() == 0)
    {
        cout << "No data - failed!" << endl;
        return 0;
    }
    
    // Apply the distortion constraints

    if( heightGrid) setupHeightDistortionParam( grd, param );
    else setupDistortionParam( grd, param );
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

    {
    bool calcCvr = param.calcGridCovar;
    Obseqn oe(2);
    int ngprm=heightGrid ? 1 : 2;
    leVector dxy(ngprm);
    SymMatrix cvr(ngprm);
    SymMatrix *ptrcvr = calcCvr ? &cvr : 0;

    pm.Start("Applying the solution to the grid",grd.nrows() );
    for( long r = 0; r < grd.nrows(); r++ ) {
        pm.Update(r);
        for( long c = 0; c < grd.ncols(); c++ ) {
            long prmNo = grd.paramNo( c, r );
            if( prmNo <= 0 ) continue;
            oe.Zero( ngprm, 0, 0 );
            oe.A(1, prmNo ) = 1.0;
            if( ! heightGrid ) oe.A(2, prmNo+1 ) = 1.0;
            if( ! oe.CalcValue( le, &dxy, ptrcvr ) )
            {
            }
#ifdef DEBUG_GRDOBSEQ
            cout << "Applying offset " << dxy(1) << " " << dxy(2) << "\n";
#endif
            grd(c,r).dxy[0] += dxy(1);
            if( ! heightGrid ) grd(c,r).dxy[1] += dxy(2);
            if( calcCvr )
            {
                grd(c,r).cvr[0] = cvr(1,1);
                if( ! heightGrid )
                {
                    grd(c,r).cvr[1] = cvr(1,2);
                    grd(c,r).cvr[2] = cvr(2,2);
                }
            }
            }
        }
    pm.Finish();
    }

    // Invert the covariance matrix

    // BLT_Matrix::SetProgressMeter( &pm );
    // le.Invert();
    // BLT_Matrix::SetProgressMeter( 0 );


    // Update the control points with calculated values and residuals

    calcControlPointListResiduals( grd, pts, gi, le, param.calcStdRes, pm );
    // calcGridStandardisedResiduals( grd, le, pm );

    return 1;
    }
