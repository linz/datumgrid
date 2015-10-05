#ifndef _PROGRESS_HPP
#define _PROGRESS_HPP
/***************************************************************************
 Copyright 2015 Crown copyright (c)
 Land Information New Zealand and the New Zealand Government.
 All rights reserved

 This program is released under the terms of the new BSD license. See 
 the LICENSE file for more information.
****************************************************************************/

class ProgressMeter {
  private:
    int CalcDisplayValue( long newValue );
  protected:
    long size;
    long curvalue;
    int resolution;
    int display;
    int nesting;
    int showing;
    int checkAbort;
    virtual void Show( const char *status ) = 0;
    virtual void Redisplay() = 0;
    virtual void Hide() = 0;
    virtual int Abort(){ return 0;}
  public:
    ProgressMeter() :
        size(0),
        nesting(0),
        resolution(100),
        showing(0),
        checkAbort(0) {}
    virtual ~ProgressMeter(){}
    void Start(const char *status, long range, long value = 0);
    int Update(long value);
    void Finish();
    };


class AsciiBarMeter : public ProgressMeter {
    int curPos;
    char *prefix;
  protected:
    virtual void Show( const char *status );
    virtual void Redisplay();
    virtual void Hide();
public:
    AsciiBarMeter( const char *prefix = 0, int size = 50 );
    ~AsciiBarMeter();
    };

#endif
