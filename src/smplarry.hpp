#ifndef SMPLARRY_HPP
#define SMPLARRY_HPP
/***************************************************************************
 Copyright 2015 Crown copyright (c)
 Land Information New Zealand and the New Zealand Government.
 All rights reserved

 This program is released under the terms of the new BSD license. See 
 the LICENSE file for more information.
****************************************************************************/

// A very simple array class.  Allocates and expands arrays. Uses pointers


class VoidPointerArray {
    void **data;
    int ndata;
    int maxdata;
    void setSize( int size );
  public:
    VoidPointerArray( int size = 16 );
    ~VoidPointerArray();
    int size(){ return ndata; }
    int add( void *item );
    void *remove( void *item );
    void *remove( int i );
    void *itemAt( int i );
    };


template <class T>
class SimplePointerArray  : private VoidPointerArray {
    char ownsdata;
    void SetSize( int size );
  public:
    SimplePointerArray( int size = 16, char ownsdata = 0 ) :
        VoidPointerArray(size), ownsdata(ownsdata) {}
    ~SimplePointerArray();
    using VoidPointerArray::size;     // Want to be able to see this
    int add( T *item ){ return VoidPointerArray::add( (void *) item ); }
    T * remove( T *item ){ return (T *) VoidPointerArray::remove( (void *) item);}
    T * remove( int i ){ return (T *) VoidPointerArray::remove( i );}
    // int deleteItem( T *item );
    // int deleteItem( int i );
    T *itemAt( int i ){ return (T *) VoidPointerArray::itemAt(i); }
    T *operator []( int i ){ return (T *) VoidPointerArray::itemAt(i);}
    };

template <class T>
SimplePointerArray<T>::~SimplePointerArray() {
   if( ownsdata ) {
      for( int i = 1; i <= size(); i++ ) { delete (T *) itemAt(i); }
      }
   }


#endif
