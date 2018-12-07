#include <iostream>
#include <sstream>
#include <iomanip>

#include <string>
#include <vector>

#include <sys/time.h>
#include <time.h>
#include <signal.h>
#include <unistd.h>
#include <errno.h>

#include <cmath>

#include <Singleton.h>
#include <ErrorLogger.h>

#include "Timer.h"

long Timer::fcountReal = 0;
long Timer::fcountVirt = 0;
long Timer::fcountProf = 0;
struct sigaction Timer::fsactReal;
struct sigaction Timer::fsactVirt;
struct sigaction Timer::fsactProf;

const int Timer::ftimeStep = 1000000; //To make sure the signal handling isn't needed

void Timer::initialize() {
   sigemptyset( &fsactReal.sa_mask );
   sigemptyset( &fsactVirt.sa_mask );
   sigemptyset( &fsactProf.sa_mask );
   fsactReal.sa_flags = 0;
   fsactVirt.sa_flags = 0;
   fsactProf.sa_flags = 0;
   fsactReal.sa_handler = &Timer::catchReal;
   fsactVirt.sa_handler = &Timer::catchVirt;
   fsactProf.sa_handler = &Timer::catchProf;
   sigaction( SIGALRM, &fsactReal, NULL );
   sigaction( SIGVTALRM, &fsactVirt, NULL );
   sigaction( SIGPROF, &fsactProf, NULL );

   itimerval dval;
   dval.it_interval.tv_sec=ftimeStep;
   dval.it_interval.tv_usec=0;
   dval.it_value.tv_sec=ftimeStep;
   dval.it_value.tv_usec=0;
   setitimer( ITIMER_REAL, &dval, 0 );
   setitimer( ITIMER_VIRTUAL, &dval, 0 );
   setitimer( ITIMER_PROF, &dval, 0 );
}

size_t Timer::start() {
   ftStartTimes.push_back(getCurrentTime());
   return ftStartTimes.size()-1;
}

void Timer::stop(size_t timer, const std::string functionName, const std::string fileName, const int lineNumber) {
   timevals end = getCurrentTime();
   timevals duration; 
   duration.realTime = end.realTime - ftStartTimes[timer].realTime;
   duration.virtTime = end.virtTime - ftStartTimes[timer].virtTime;
   duration.profTime = end.profTime - ftStartTimes[timer].profTime;

   std::ostringstream omsg;
   omsg << "Time elapsed: "
        << prettyPrint(duration.realTime) <<" (Real);  "
        << prettyPrint(duration.virtTime) <<" (Process); "
        << prettyPrint(duration.profTime-duration.virtTime) <<" (System)";

   utl::ErrorLogger::GetInstance().Log(utl::ErrorLogger::eInfo, functionName, fileName, lineNumber, omsg.str());
}

std::string Timer::prettyPrint( timeval t ) {
   std::ostringstream ostr;
   bool hour=false;

   if (t.tv_sec > 3600) {
      ostr << t.tv_sec/3600 << "h:";
      t.tv_sec -= t.tv_sec/3600 * 3600;
      hour=true;
   }

   if (t.tv_sec > 60 || hour) {
      ostr << std::setw(2) << t.tv_sec/60 << "m:";
      t.tv_sec -= t.tv_sec/60 * 60;
   }

   ostr << std::setw(6) << std::setprecision(3) << float(t.tv_sec) + t.tv_usec/1e6 <<"s";

   return ostr.str();
}

timevals Timer::getCurrentTime() {
   itimerval ittv;
   timevals tvs;

   // getitimer returns the interval and the time left, hence the
   // subtraction
   getitimer( ITIMER_REAL, &ittv );
   tvs.realTime = ittv.it_interval-ittv.it_value;
   tvs.realTime.tv_sec += ftimeStep*fcountReal;
   
   getitimer( ITIMER_VIRTUAL, &ittv );
   tvs.virtTime = ittv.it_interval-ittv.it_value;
   tvs.virtTime.tv_sec += ftimeStep*fcountVirt;
   
   getitimer( ITIMER_PROF, &ittv );
   tvs.profTime = ittv.it_interval-ittv.it_value;
   tvs.profTime.tv_sec += ftimeStep*fcountProf;

   return tvs;
}

timeval& operator += (timeval& a, const timeval& b) {
   a.tv_sec += b.tv_sec;
   a.tv_usec += b.tv_usec;
   if (a.tv_usec > 999999) {
      ++a.tv_sec;
      a.tv_usec -= 1000000;
   }
}

timeval operator + (const timeval& a, const timeval& b) {
   timeval x(a);
   x += b;
   return x;
}

timeval& operator -= (timeval& a, const timeval& b) {
   a.tv_sec -= b.tv_sec;
   a.tv_usec -= b.tv_usec;
   if (a.tv_usec < 0) {
      --a.tv_sec;
      a.tv_usec += 1000000;
   }
}

timeval operator - (const timeval& a, const timeval& b) {
   timeval x(a);
   x -= b;
   return x;
}
