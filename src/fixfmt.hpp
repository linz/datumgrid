#ifndef _FIXFMT_HPP
#define _FIXFMT_HPP

/***************************************************************************
 Copyright 2015 Crown copyright (c)
 Land Information New Zealand and the New Zealand Government.
 All rights reserved

 This program is released under the terms of the new BSD license. See 
 the LICENSE file for more information.
****************************************************************************/

////////////////////////////////////////////////////////////////////
// class FixedFormat:
//
// Instances of FixedFormat act as manipulators for fixed format
// output.  Each instance has an associated width and precision
// (the width may be 0, in which case any pending setw() width will
// be used.
//

class FixedFormat {
   int fpWidth;
   int fpPrecision;
 public:
   FixedFormat( int w, int p) : fpWidth(w), fpPrecision(p){};
   FixedFormat( int p ) : fpWidth(0), fpPrecision( p ){};
   void width( int w ){ fpWidth = w; };
   void precision( int p ) { fpPrecision = p; }
   friend ostream & operator <<( ostream &os, const FixedFormat &fp );
   };

inline
ostream & operator << ( ostream &os, const FixedFormat &fp ) {
   os.setf( ios::showpoint | ios::fixed | ios::right );
   os.precision( fp.fpPrecision );
   if( fp.fpWidth ) os.width( fp.fpWidth );
   return os;
   }

#endif
