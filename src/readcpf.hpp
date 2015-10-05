#ifndef READCPF_HPP
#define READCPF_HPP

/***************************************************************************
 Copyright 2015 Crown copyright (c)
 Land Information New Zealand and the New Zealand Government.
 All rights reserved

 This program is released under the terms of the new BSD license. See 
 the LICENSE file for more information.
****************************************************************************/
#include "contrlpt.hpp"

int readControlPointFile( string filename, ControlPointList &list, bool heightPoints=false );

#endif
