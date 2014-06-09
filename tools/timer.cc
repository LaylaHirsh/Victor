#define VERSIONSTRING "$Date: 2008-05-09 13:58:51 $ $Revision: 1.4 $" 
//  Timer class for benchmarking and more.
//
//  @version 0.1

#include "timer.h"

int Timer::seconds() const
{
  if (endTime != 0) {
    return endTime-startTime;
  } 
  else {
    time_t now = time(0);
    return now-startTime;
  }
}

int Timer::minutes() const
{
  if (endTime != 0) {
    return (endTime-startTime)/60;
  } 
  else {
    time_t now =time(0);
    return (now-startTime)/60;
  }
}

int Timer::hours() const
{
  if (endTime != 0) {
    return (endTime-startTime)/3600;
  } 
  else {
    time_t now = time(0);
    return (now-startTime)/3600;
  }
}




