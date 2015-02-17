// Module to read a control point file


#include <fstream>
#include <sstream>
#include <iostream>

using namespace std;

#include "readcpf.hpp"

const int MAXREC = 100;

int readControlPointFile( string filename, ControlPointList &list ){
    ifstream ifs( filename.c_str() );
    if( ! ifs.good() ) {
       cerr << "Cannot open control point file " << filename << endl;
       return 0;
       }


    char buf[MAXREC];
    string id;
    double x, y, dx, dy;
    double error;
    string ptclass;
    int nrec = 0;

    while( ifs.good() ){
       // Get a record into an input string stream...
       ifs.getline( buf, MAXREC );
       nrec++;
       int nchar=ifs.gcount()-1;
       if( nchar <= 0 ) continue;

       istringstream record( buf );

       record >> id;

       // Skip blank strings and comments
       if( ! record.good() || id[0] == '!' || id[0] == '#' ) continue;

       record >> x >> y >> dx >> dy >> ptclass;
       if( ! record.fail() ) {
          ControlPoint *cpt=new ControlPoint( id, x, y, dx, dy, ptclass );
          record >> error;
          if( ! record.fail() ) cpt->setError( error );
          list.add( cpt );
          }
       else {
          cerr << "Invalid data in control point file " << filename
               << " at record " << nrec << "\n";
          }

       }
    return list.size();
    }
