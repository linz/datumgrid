#ifndef CONTRLPT_HPP
#define CONTRLPT_HPP
/***************************************************************************
 Copyright 2015 Crown copyright (c)
 Land Information New Zealand and the New Zealand Government.
 All rights reserved

 This program is released under the terms of the new BSD license. See 
 the LICENSE file for more information.
****************************************************************************/

#include <cstring>
#include "smplarry.hpp"

class ControlPointClass {
   public:
     static ControlPointClass &named( const string &ptclass );
     void setError( double wgt );
     double getError();
     void setRejected( char rejected = 1 ){ this->rejected = rejected;}
     char isRejected(){ return rejected; }
     static void setDefaultError( double wgt );
     const string &getName(){ return name; }
     void clearSumStdRes();
     void addStdRes( double sr, char used = 1 );
     long stdResCount( char used = 1 ){ return used ? usedCount : unusedCount; }
     double RMS_StdRes( char used = 1 );
     static int count() { return list.size(); }
     static ControlPointClass *classNumber( int i ){ return list[i]; }
   private:
     ControlPointClass( const string &ptclass );
     static ControlPointClass *findClassNamed( const string &ptclass );
     string name;
     double error;
     char rejected;
     long usedCount;
     double usedSSR;
     long unusedCount;
     double unusedSSR;
     static double defaultError;
     static SimplePointerArray<ControlPointClass> list;
   };


class ControlPoint {
   public:
      ControlPoint( const string &id, double x, double y, double dx, double dy, string &ptclass );
      void setError( double wgt );
      double getError();
      void setRejected( char rejected = 1 ){ this->rejected = rejected; }
      char isRejected();
      ControlPointClass &getClass(){ return ptclass; }
      const string &getId(){ return id; }
      double *coord(){ return xy;}
      double *offset(){ return dxy;}
      double *calcOffset(){ return cxy; }
      double &distanceResidual() { return distres; }
      double &stdResidual() { return stdres; }

   private:
      string id;
      double xy[2];
      double dxy[2];
      double cxy[2];  /* Calculated value */
      double distres;
      double stdres;
      double error;
      char rejected;
      ControlPointClass &ptclass;
   };


class ControlPointList : public SimplePointerArray<ControlPoint> {
   public:
      ControlPointList() : SimplePointerArray<ControlPoint>( 32, 1 ) {;}
      ControlPoint * operator[] ( const string &id );
      ControlPoint * operator[] ( int i ){ return SimplePointerArray<ControlPoint>::itemAt(i); }
   };


#endif
