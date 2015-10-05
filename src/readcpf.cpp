
/***************************************************************************
 Copyright 2015 Crown copyright (c)
 Land Information New Zealand and the New Zealand Government.
 All rights reserved

 This program is released under the terms of the new BSD license. See 
 the LICENSE file for more information.
****************************************************************************/
// Module to read a control point file


#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

using namespace std;

#include "readcpf.hpp"

const int MAXREC = 100;

int readControlPointFile( string filename, ControlPointList &list, bool heightPoints ){
    ifstream ifs( filename.c_str() );
    if( ! ifs.good() ) {
       cerr << "Cannot open control point file " << filename << endl;
       return 0;
       }

    string id;
    string buffer;
    double x, y, dx, dy;
    double error;
    string ptclass;
    int nrec = 0;
    dy=0.0;

    // Get a record into an input string stream...
    while( getline(ifs,buffer) ){
       nrec++;
       replace(buffer.begin(),buffer.end(),',',' ');
       istringstream record( buffer );
       record >> id;

       // Skip blank strings and comments
       if( ! record.good() || id[0] == '!' || id[0] == '#' ) continue;

       record >> x >> y >> dx;
       if( ! heightPoints ) record >> dy;
       if( ! record.fail() ) {
          record  >> ptclass;
          if( record.fail() ) ptclass="default";
          ControlPoint *cpt=new ControlPoint( id, x, y, dx, dy, ptclass );
          record >> error;
          if( ! record.fail() ) cpt->setError( error );
          list.add( cpt );
          }
       else if( nrec > 1 ) {
          cerr << "Invalid data in control point file " << filename
               << " at record " << nrec << "\n";
          }

       }
    return list.size();
    }
