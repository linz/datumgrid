#include <iostream>
#include <fstream>
#include <strstream>
#include <cstring>
#include <iomanip>

using namespace std;

#include "contrlpt.hpp"
#include "grdintrp.hpp"
#include "grdobseq.hpp"
#include "readcpf.hpp"
#include "grid.hpp"
#include "fixfmt.hpp"

const int MAXREC = 256;

int readPositiveNumber( istream &is, double &val, string &errmess ) {
   is >> val;
   if( is.fail() ) {
      errmess = "Missing value";
      return 0;
      }
   else if ( val <= 0.0 ) {
      errmess = "Value must be positive";
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
       istrstream record( buf, count-1 );
       // cout << "\"" << buf << "\" " << count << "\n";
       record >> command;
       if( !record.good() || command[0] == '!' ) continue;

       error = "";
       if( command == "data_file" ) {
           string df;
           record >> df;
           if( record.fail() )
           {
              error = "Data file name missing";
              }
           else {
              if( ! readControlPointFile( df, pts ) ) {
                 error = "No control points in control point file\n";
                 }
              }
           }

       else if ( command == "grid_spacing" ) {
           readPositiveNumber( record, param.spacing, error );
           }

       else if ( command ==  "required_point_proximity" ) {
           readPositiveNumber( record, param.maxPointProximity, error );
           }

       else if ( command == "zero_outside_proximity" ) {
           param.zeroOutsideProximity=true;
           }


       else if ( command == "distortion_error" ) {
           readPositiveNumber( record, param.distortionError, error );
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
      cout << "Syntax: datumgrd command_file_name output_file\n\n";
      cout << "The output file name is the root name for three output files\n";
      cout << "   xxxxx_cpt.csv  xxxxx_grd.csv  xxxxx_def.csv\n";
      return 0;
      }

   if( ! readCommandFile( argv[1], param, points ) ) {
      return 0;
      }

   string rootfilename = argv[2];

// Set up an interpolator

   BilinearInterpolator gi;

// Try to create a grid...

   param.pointInfluenceRange = gi.pointInfluenceRange();
   Grid grid( param, points );
   cout << "Grid specification" << endl;
   grid.dumpSpecTo(cout);

   // grid.dumpSpecTo( cout );

   int success=CalculateGridModel( grid, points, gi, param.distortionError  );
   if( ! success ) {
      cout << "Failed to calculate grid model" << endl;
      }

   string logfilename = rootfilename + "_log.txt";
   ofstream logfile( logfilename.c_str() );

   // Write grid file... */
   {
      string outputfile;
      outputfile = rootfilename + "_grd.csv";
      ofstream grdfile( outputfile.c_str() );
      cout << "Writing grid to " << outputfile << endl;
      grdfile << "x,y,dx,dy,calc\n";
      for( long r = 0; r < grid.nrows(); r++ ) for (long c = 0; c < grid.ncols(); c++ ) {
         if( grid.isValidPoint(c,r) ) {
            long cr[2] = {c,r};
            double xy[2];
            grid.convert( cr, xy );
            double *offset = grid(c,r).dxy;
            int calc = grid(c,r).inrange ? 1 : 0;
            grdfile << FixedFormat(param.ndpCoord) << xy[0] << "," << xy[1] << ","
                    << FixedFormat(param.ndpValue) << offset[0] << "," << offset[1]
                    << "," << calc
                    << endl;
            }
         }
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
      cptfile << "id,x,y,dx,dy,calcdx,calcdy,residual,stdres,class,error,used\n";
      for( long i = 0; i < points.size(); i++ ) {
         ControlPoint &cpt = * points[i];
         cptfile << "\"" << cpt.getId() << "\","
                 << FixedFormat(param.ndpCoord) << cpt.coord()[0] << "," << cpt.coord()[1] << ","
                 << FixedFormat(param.ndpValue) << cpt.offset()[0] << "," << cpt.offset()[1] << ","
                 << cpt.calcOffset()[0] << "," << cpt.calcOffset()[1] << ","
                 << cpt.distanceResidual() << "," << cpt.stdResidual() << ",\""
                 << cpt.getClass().getName() << "\","
                 << cpt.getError() << "," 
                 << (cpt.isRejected() ? 0 : 1) <<  endl;
         cpt.getClass().addStdRes( cpt.stdResidual(), cpt.isRejected() ? 0 : 1 );
         }
      }

   {
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
   {
      string outputfile;
      outputfile = rootfilename + "_def.csv";
      ofstream deffile( outputfile.c_str() );
      cout << "Writing deformation to " << outputfile << endl;
      writeGridDistortion( grid, deffile );
      }

   // Write Surfer format grid files
   cout << "Writing surfer format grid files" << endl;
   grid.writeSurferFiles( rootfilename );
   return 0;
   }
