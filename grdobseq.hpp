#ifndef GRDOBSEQ_HPP
#define GRDOBSEQ_HPP

#include <iostream>

#ifndef GRID_HPP
#include "grid.hpp"
#endif

#ifndef LINEQN_HPP
#include "lineqn.hpp"
#endif

#ifndef PROGRESS_HPP
#include "progress.hpp"
#endif

#ifndef GRDINTRP_HPP
#include "grdintrp.hpp"
#endif

#ifndef CONTRLPT_HPP
#include "contrlpt.hpp"
#endif

void setupBandwidthDefinition( Grid &grd, BLT_Def &blt, int ptInfluenceRange );
void setupDistortionParam( Grid &grd, double dstWeight );
long sumDistortionConstraints( Grid &grd, LinearEquations &le, ProgressMeter &pm );
long sumControlPoints( Grid &grd, ControlPointList &pts, GridInterpolator &gi,
    LinearEquations &le, ProgressMeter &pm );
int CalculateGridModel( Grid &grd, ControlPointList &pts,
                        GridInterpolator &gi, double dstError );
void writeGridDistortion( Grid &grd, int coord_precision, ostream &os );

#endif
