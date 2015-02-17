
#include <stdlib.h>
#include <string.h>

using namespace std;

#include "smplarry.hpp"

typedef void *pvoid;

VoidPointerArray::VoidPointerArray( int size ) {
   if( size <= 0 ) size = 16;
   data = new pvoid[size];
   maxdata = size;
   ndata = 0;
   }

VoidPointerArray::~VoidPointerArray() {
   delete [] data;
   }

void VoidPointerArray::setSize( int size ) {
   pvoid *newData;
   if( size > ndata ) {
      newData = new pvoid[size];
      if( ndata > 0 ) memcpy( newData, data, ndata * sizeof(pvoid) );
      delete [] data;
      maxdata = size;
      data = newData;
      }
   }

int VoidPointerArray::add( void * item ) {
   if( ndata >= maxdata ) setSize( maxdata * 2 );
   data[ndata] = item;
   ndata++;
   return ndata;
   }

void * VoidPointerArray::itemAt( int i ) {
   if( i < 0 || i >= ndata ) return 0;
   return data[i];
   }

void *VoidPointerArray::remove( void *item ) {
   for( int i = 1; i <= ndata; i++ ) {
      if( itemAt(i) == item ) return remove( i );
      }
   return 0;
   }

void *VoidPointerArray::remove( int i ) {
   if( i < 0 || i >= ndata ) return 0;
   void *result = data[i];
   for( ; i < ndata; i++ ) { data[i] = data[i+1];  }
   ndata--;
   return result;
   }
