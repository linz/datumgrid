/***************************************************************************
 Copyright 2015 Crown copyright (c)
 Land Information New Zealand and the New Zealand Government.
 All rights reserved

 This program is released under the terms of the new BSD license. See 
 the LICENSE file for more information.
****************************************************************************/

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cmath>
#include <iomanip>

using namespace std;

#include "contrlpt.hpp"
#include "grdintrp.hpp"
#include "grdobseq.hpp"
#include "readcpf.hpp"
#include "grid.hpp"
#include "fixfmt.hpp"

const int MAXREC = 256;

int readPositiveNumber( istream &is, double &val, string &errmess, bool optional=false, bool zeroOk=false ) {
   is >> val;
   if( is.fail() ) {
      if( ! optional ) errmess = "Missing value";
      return 0;
      }
   else if ( val < 0.0 || (val == 0.0 && ! zeroOk) ) {
      errmess = "Value must be positive";
      return 0;
      }
   else {
      return 1;
      }
   }

int readNumber( istream &is, double &val, string &errmess, bool optional=false ) {
   is >> val;
   if( is.fail() ) {
      if( ! optional ) errmess = "Missing value";
      return 0;
      }
   else {
      return 1;
      }
   }

int readErrorOrReject( istream &is, double &val, char &reject, string &errmess ) {
   string option;
   is >> option;
   // string::set_case_sensitive( 0 );
   reject = 0;
   if( option == "reject" ) {
      reject = 1;
      return 1;
      }
   else if( option == "error" ) {
      return readPositiveNumber( is, val, errmess );
      }
   else {
      errmess = "Invalid option ";
      errmess += option;
      return 0;
      }
   }

int readBool( istream &is, bool &val, string &errmess ) {
   string option;
   is >> option;
   // string::set_case_sensitive( 0 );
   bool result=true;
   if( option == "false" || option == "no" ) result=false;
   val=result;
   return 1;
}


int readCommandFile( char *filename, GridParams &param, ControlPointList &pts ) {
   ifstream ifs(filename);
   if( ! ifs.good() ) {
       cerr << "Cannot open command file " << filename << endl;
       return 0;
       }

   // string::set_case_sensitive(0);

   char buf[MAXREC];
   string command, error;
   int ok = 1;
   int nrec = 0;

   while( ifs.good() ) {
       ifs.getline( buf, MAXREC );
       nrec++;
       int count=ifs.gcount();
       if( count < 2 ) continue;
       istringstream record( buf );
       // cout << "\"" << buf << "\" " << count << "\n";
       record >> command;
       if( record.fail() || command[0] == '!' ) continue;
       if( command[0] == '#' ) continue;

       error = "";
       if( command == "data_file" ) {
           string df;
           record >> df;
           if( record.fail() )
           {
              error = "Data file name missing";
              }
           else {
              if( ! readControlPointFile( df, pts, param.heightGrid ) ) {
                 error = "No control points in control point file\n";
                 }
              }
           }

       else if ( command == "height_grid" ) {
           if( pts.size() > 0 )
           {
               error = "height_grid option must be before data_file";
           }
           readBool(record,param.heightGrid,error);
           if( param.heightGrid ) param.dxcolname="dh";
           if( param.heightGrid ) param.xerrname="stdh";
           }
       else if ( command == "grid_spacing" ) {
           readPositiveNumber( record, param.xSpacing, error );
           if( ! readPositiveNumber( record, param.ySpacing, error, true ))
           {
               param.ySpacing=param.xSpacing;
           }
           }

       else if ( command == "coordinate_to_metres" ) {
           readPositiveNumber( record, param.xScale, error );
           if( ! readPositiveNumber( record, param.yScale, error, true ))
           {
               param.yScale=param.xScale;
           }
           }

       else if ( command ==  "required_point_proximity" ) {
           readPositiveNumber( record, param.maxPointProximity, error );
           }

       else if ( command == "beyond_proximity" ) {
           string option;
           record >> option;
           if( option == "fit" )
           {
               param.boundaryOption=GridParams::grdFit;
           }
           else if( option == "zero" )
           {
               param.boundaryOption=GridParams::grdZero;
           }
           else if( option == "ignore" )
           {
               param.boundaryOption=GridParams::grdIgnore;
           }
           else
           {
               error = "beyond_proximity option must be one of \"fit\", \"zero\", or \"ignore\"";
           }
           }
       else if ( command == "height_zero_value" ) {
           readNumber( record, param.heightZero, error );
           }
       else if ( command == "distortion_error" ) {
           readPositiveNumber( record, param.distortionError, error );
           }

       else if ( command == "shear_weight" ) {
           readPositiveNumber( record, param.shearWeight, error, false, true  );
           }
       else if ( command == "scale_weight" ) {
           readPositiveNumber( record, param.scaleWeight, error, false, true  );
           }
       else if ( command == "constant_weight" ) {
           readPositiveNumber( record, param.nonConstantWeight, error, false, true  );
           }
       else if ( command == "non_linear_weight" ) {
           readPositiveNumber( record, param.nonLinearWeight, error, false, true  );
           }

       else if ( command == "default_point_error" ) {
           double dfltError;
           if( readPositiveNumber( record, dfltError, error ) ) {
              ControlPointClass::setDefaultError( dfltError );
              }
           }

       else if ( command == "coordinate_precision" ) {
           record >> param.ndpCoord;
           if( record.fail() || param.ndpCoord < 0 || param.ndpCoord > 10 )
           {
               error="Invalid coordinate precision, must be between 0 and 10";
           }
           }
       else if ( command == "value_precision" ) {
           record >> param.ndpValue;
           if( record.fail() || param.ndpValue < 0 || param.ndpValue > 10 )
           {
               error="Invalid value precision, must be between 0 and 10";
           }
           }
       else if ( command == "point" ) {
           string id;
           record >> id;
           ControlPoint *cpt = pts[id];
           if( !cpt ) {
               error = "Invalid point id ";
               error += id;
               }
           else {
               double errVal;
               char reject;
               if( readErrorOrReject( record, errVal, reject, error ) ) {
                   if( reject ) {
                      cpt->setRejected();
                      }
                   else {
                      cpt->setError( errVal );
                      }
                    }

                }
           }

       else if ( command == "class" ) {
           string clsName;
           record >> clsName;
           if( clsName == "" ) {
               error = "Missing class name";
               }
           else {
               ControlPointClass &cpc = ControlPointClass::named(clsName);
               double errVal;
               char reject;
               if( readErrorOrReject( record, errVal, reject, error )) {
                   if( reject ) {
                      cpc.setRejected();
                      }
                   else {
                      cpc.setError( errVal );
                      }
                    }
               }

           }

       else if ( command == "columns" ) {
           record >> param.xcolname >> param.ycolname >> param.dxcolname;
           if( ! param.heightGrid ) record >> param.dycolname;
           if( record.fail() )
           {
               error="Missing column names";
           }
           }
       else if ( command == "uncertainty_columns" ) {
           record >> param.xerrname;
           if( ! param.heightGrid ) record >> param.yerrname >> param.xycorrname;
           if( record.fail() )
           {
               error="Missing column names";
           }
           }

       else if( command == "print_grid_params" )
       {
           readBool(record,param.printGridParams,error);
       }
       else if( command == "fill_grid" )
       {
           readBool(record,param.fillGrid,error);
       }
       else if( command == "calculate_grid_uncertainty" )
       {
           readBool(record,param.calcGridCovar,error);
       }
       else if( command == "calculate_control_point_stdres" )
       {
           readBool(record,param.calcStdRes,error);
       }
       else {
           error = "Invalid command";
           }
       if( error != "" ) {
           cerr << "\nError in command file " << filename << " record " << nrec << "\n";
           cerr << "Command " << command << ": " << error << "\n";
           ok = 0;
           }
       }
   return ok;
   }


int main( int argc, char *argv[] ) {

   ControlPointList points;
   GridParams param;

   if( argc < 3 ) {
      cout << "Syntax: datumgrid command_file_name output_file\n\n";
      cout << "The output file name is the root name for three output files\n";
      cout << "   xxxxx_cpt.csv  xxxxx_grd.csv  xxxxx_def.csv\n";
      return 0;
      }

   if( ! readCommandFile( argv[1], param, points ) ) {
      return 0;
      }


   bool heightGrid = param.heightGrid;
   string rootfilename = argv[2];
   string logfilename = rootfilename + "_log.txt";
   ofstream logfile( logfilename.c_str() );
   logfile << "Config file: " << argv[0] << endl;

// Set up an interpolator

   BilinearInterpolator gi;

// Try to create a grid...

   param.pointInfluenceRange = gi.pointInfluenceRange();
   Grid grid( param, points );
   cout << "Grid specification" << endl;
   grid.dumpSpecTo(cout);
   logfile << "Grid specification" << endl;
   grid.dumpSpecTo(logfile);

   // grid.dumpSpecTo( cout );

   int success=CalculateGridModel( grid, points, gi, param );
   if( ! success ) {
      cout << "Failed to calculate grid model" << endl;
      logfile << "Failed to calculate grid model" << endl;
      }

   // Write grid file... */
   {
      string outputfile;
      outputfile = rootfilename + "_grd.csv";
      ofstream grdfile( outputfile.c_str() );
      cout << "Writing grid to " << outputfile << endl;
      // grdfile << "x,y,dx,dy,calc,paramno\n";
      grdfile << param.xcolname << ","
              << param.ycolname;
     
      if( heightGrid ) 
      {
          grdfile << "," << param.dxcolname;
          if( param.calcGridCovar )
          {
              grdfile << "," << param.xerrname;
          }
      }
      else
      {
          grdfile << "," << param.dxcolname << "," << param.dycolname;
          if( param.calcGridCovar )
          {
              grdfile << "," << param.xerrname
                      << "," << param.yerrname
                      << "," << param.xycorrname;
          }
      }
      if( param.printGridParams ) grdfile << ",c,r,mode,paramno";
      grdfile << "\n";
      double zeroOffset[2]={0,0};
      for( long r = 0; r < grid.nrows(); r++ ) for (long c = 0; c < grid.ncols(); c++ ) {
         double xy[2];
         double *covar;
         grid.convert( c,r, xy );
         double *offset=zeroOffset;
         int paramno=0;
         string mode="fill";
         if( grid.isValidPoint(c,r) ) {
            GridPoint &gp=grid(c,r);
            offset = gp.dxy;
            covar = gp.cvr;
            paramno=gp.paramno;
            if(paramno < 0 ) mode="ignore"; 
            else if(paramno <= 0 ) mode="zero"; 
            else mode="calc";
            }
         else if( ! param.fillGrid ) continue;
         grdfile << FixedFormat(param.ndpCoord) << xy[0] << "," << xy[1] << ","
                    << FixedFormat(param.ndpValue) << offset[0];
         if( ! heightGrid ) grdfile << "," << offset[1];
         if( param.calcGridCovar )
         {
             grdfile << "," << sqrt(covar[0]);
             if( ! heightGrid )
             {
                 double xycorr = sqrt(covar[0]*covar[2]);
                 if( xycorr > 0.0 ) xycorr=covar[1]/xycorr;
                 grdfile << "," << sqrt(covar[2]) << "," << xycorr;
             }
         }
         if( param.printGridParams ) grdfile << "," << c << "," << r << "," << mode << "," << paramno;
         grdfile << endl;
         }
      logfile << "Grid written to " << outputfile << endl;
      }

   if( ! success ) return 0;

   // Write control point file.  Also sum standardised residuals for classes
   {
      for( int i = 0; i < ControlPointClass::count(); i++ ) {
         ControlPointClass::classNumber(i)->clearSumStdRes();
         }
      }
   {
      string outputfile;
      outputfile = rootfilename + "_cpt.csv";
      cout << "Writing control points to " << outputfile << endl;
      ofstream cptfile( outputfile.c_str() );
      cptfile << "id,";
      cptfile << param.xcolname << ",";
      cptfile << param.ycolname << ",";
      cptfile << param.dxcolname << ",";
      if( ! heightGrid ) cptfile << param.dycolname << ",";
      cptfile << "calc" << param.dxcolname << ",";
      if( ! heightGrid ) cptfile << "calc" << param.dycolname << ",";
      cptfile << "res" << param.dxcolname << ",";
      if( ! heightGrid ) cptfile << "res" << param.dycolname << "," << "residual,";
      cptfile << "stdres,class,error,used\n";

      for( long i = 0; i < points.size(); i++ ) {
         ControlPoint &cpt = * points[i];
         cptfile << "\"" << cpt.getId() << "\",";
         cptfile << FixedFormat(param.ndpCoord) << cpt.coord()[0] << "," << cpt.coord()[1] << ",";
         cptfile << FixedFormat(param.ndpValue);
         cptfile<< cpt.offset()[0] << ",";
         if( ! heightGrid ) cptfile << cpt.offset()[1] << ",";
         cptfile << cpt.calcOffset()[0] << ",";
         if( ! heightGrid ) cptfile << cpt.calcOffset()[1] << ",";
         cptfile << (cpt.offset()[0]-cpt.calcOffset()[0]) << ",";
         if( ! heightGrid ) cptfile << (cpt.offset()[1]-cpt.calcOffset()[1]) << "," << cpt.distanceResidual() << ",";
         cptfile << cpt.stdResidual() << ",\"";
         cptfile << cpt.getClass().getName() << "\",";
         cptfile << cpt.getError() << "," ;
         cptfile << (cpt.isRejected() ? 0 : 1) <<  endl;
         cpt.getClass().addStdRes( cpt.stdResidual(), cpt.isRejected() ? 0 : 1 );
         }
      }

   if( 0 ){
      logfile << "Summary of residuals by class\n";
      for( int i = 0; i < ControlPointClass::count(); i++ ) {
         ControlPointClass &cpc = *ControlPointClass::classNumber(i);
         logfile << setw(10)
                 << cpc.getName() << "    used  "
                 << setw(5) << cpc.stdResCount(1) << "  "
                 << FixedFormat(7,2) << cpc.RMS_StdRes(1)
                 << "     unused  "
                 << setw(5) << cpc.stdResCount(0) << "  "
                 << FixedFormat(7,2) << cpc.RMS_StdRes(0)
                 << "\n";
         }
      }

   // Write deformation file... */
   if( ! heightGrid ) {
      string outputfile;
      outputfile = rootfilename + "_def.csv";
      ofstream deffile( outputfile.c_str() );
      cout << "Writing deformation to " << outputfile << endl;
      writeGridDistortion( grid, param, deffile );
      }

   /*
   // Write Surfer format grid files
   cout << "Writing surfer format grid files" << endl;
   grid.writeSurferFiles( rootfilename );
   */
   return 0;
   }
