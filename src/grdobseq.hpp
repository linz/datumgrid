#ifndef GRDOBSEQ_HPP
#define GRDOBSEQ_HPP
/***************************************************************************
 Copyright 2015 Crown copyright (c)
 Land Information New Zealand and the New Zealand Government.
 All rights reserved

 This program is released under the terms of the new BSD license. See 
 the LICENSE file for more information.
****************************************************************************/

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
void setupDistortionParam( Grid &grd, GridParams &param );
long sumDistortionConstraints( Grid &grd, LinearEquations &le, ProgressMeter &pm );
long sumControlPoints( Grid &grd, ControlPointList &pts, GridInterpolator &gi,
    LinearEquations &le, ProgressMeter &pm );
int CalculateGridModel( Grid &grd, ControlPointList &pts,
                        GridInterpolator &gi, GridParams &param );
void writeGridDistortion( Grid &grd, GridParams &param, ostream &os );

#endif
