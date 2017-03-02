
/***************************************************************************
 Copyright 2015 Crown copyright (c)
 Land Information New Zealand and the New Zealand Government.
 All rights reserved

 This program is released under the terms of the new BSD license. See 
 the LICENSE file for more information.
****************************************************************************/


#include <string>
#include <math.h>

using namespace std;

#include "contrlpt.hpp"
#include "smplarry.hpp"

/////////////////////////////////////////////////////////////////////////
// class ControlPointClass

double ControlPointClass::defaultError = 1.0;

SimplePointerArray<ControlPointClass> ControlPointClass::list( 8, 1);


ControlPointClass::ControlPointClass( const string &ptclass ) :
    name(ptclass), error(0.0), rejected(0),
    usedCount(0), usedSSR(0.0),
    unusedCount(0), unusedSSR(0.0) {};

ControlPointClass &ControlPointClass::named( const string &ptclass ) {
   ControlPointClass *cls;
   cls = findClassNamed( ptclass );
   if( ! cls ) {
      cls = new ControlPointClass( ptclass );
      list.add( cls );
      }
   return * cls;
   }

ControlPointClass *ControlPointClass::findClassNamed( const string &ptclass ) {
   // int cs = string::set_case_sensitive(0);
   ControlPointClass *cpc = 0;
   for( int i = 0; i < list.size(); i++ ) {
      cpc = list[i];
      if( cpc->name == ptclass ) break;
      cpc = 0;
      }
   // string::set_case_sensitive( cs );
   return cpc;
   }

void ControlPointClass::setError( double wgt ) {
   error = fabs(wgt);
   }

void ControlPointClass::setDefaultError( double wgt ) {
   defaultError = fabs(wgt);
   }

double ControlPointClass::getError() {
   return error > 0.0 ? error : defaultError;
   }

void ControlPointClass::clearSumStdRes() {
   usedCount = unusedCount = 0;
   usedSSR = unusedSSR = 0;
   }

void ControlPointClass::addStdRes( double sr, char used ) {
   if( used ) {
      usedCount++; usedSSR += sr*sr;
      }
   else {
      unusedCount++; unusedSSR += sr*sr;
      }
   }

double ControlPointClass::RMS_StdRes( char used ) {
   double result = 1.0;
   if( used ) {
      if( usedCount ) result = sqrt( usedSSR / usedCount );
      }
   else {
      if( unusedCount ) result = sqrt( unusedSSR / unusedCount );
      }
   return result;
   }

/////////////////////////////////////////////////////////////////////////
// class ControlPoint

ControlPoint::ControlPoint( const string &id, double x, double y, double dx, double dy, string &ptclass ) :
   id(id), ptclass( ControlPointClass::named(ptclass)) {
   xy[0] = x; xy[1] = y;
   dxy[0] = dx; dxy[1] = dy;
   cxy[0] = cxy[1] = distres = stdres = 0.0;
   error = 0.0;
   rejected = 0;
   }

void ControlPoint::setError( double wgt ) {
   error = fabs( wgt );
   }

double ControlPoint::getError() {
   return error > 0.0 ? error : ptclass.getError();
   }

bool ControlPoint::isRejected() {
   return rejected || ptclass.isRejected() || unused;
   }


/////////////////////////////////////////////////////////////////////////
// class ControlPointList

ControlPoint *ControlPointList::operator [] ( const string &id ) {
    // int cs = string::set_case_sensitive( 0 );
    ControlPoint *cpt;
    for( int i = 0; i < size(); i++ ) {
       cpt = itemAt(i);
       if( cpt->getId() == id ) break;
       cpt = 0;
       }
    // string::set_case_sensitive( cs );
    return cpt;
    }

bool ControlPointList::isUsed()
{
    for( int i = 0; i < size(); i++ ) {
        if( ! itemAt(i)->isUsed() ) return true;
    }
    return false;
}
