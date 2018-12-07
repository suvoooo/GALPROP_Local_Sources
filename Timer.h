
#ifndef _Timer_h_
#define _Timer_h_

#include <string>
#include <vector>

#include <sys/time.h>
#include <time.h>
#include <signal.h>

#include <Singleton.h>


struct timevals {
   timeval realTime;
   timeval virtTime;
   timeval profTime;
};

timeval& operator += (timeval& a, const timeval& b);
timeval& operator -= (timeval& a, const timeval& b);
timeval operator + (const timeval& a, const timeval& b);
timeval operator - (const timeval& a, const timeval& b);

class Timer : public utl::LeakingSingleton<Timer> {
   public:
      size_t start();
      void stop(const size_t timer, 
		const std::string functionName,
		const std::string fileName,
		const int lineNumber);
      static void initialize();
      static void catchReal(int sig) { ++fcountReal; }
      static void catchVirt(int sig) { ++fcountVirt; }
      static void catchProf(int sig) { ++fcountProf; }
   private:
      Timer () {};
      ~Timer () {};

      timevals getCurrentTime();
      std::string prettyPrint(timeval t);

      static const int ftimeStep;
      static long fcountReal, fcountVirt, fcountProf;
      static struct sigaction fsactReal, fsactVirt, fsactProf;
      std::vector<timevals> ftStartTimes;

      friend class utl::LeakingSingleton<Timer>;
};

// A macro to time a subroutine (discarding return value).  Function should
// be both function name and the variables, e.g. TIME_SUBROUTINE(min(a,b));
#define TIME_SUBROUTINE(function, ...) \
do {const size_t index = Timer::GetInstance().start(); \
function(__VA_ARGS__); \
Timer::GetInstance().stop(index, #function, __FILE__, __LINE__);} while(0)

// A macro to time a function returning a value.  Assignment is everything left
// of the equal sign.  Note that the lines are contained in a block, so any
// variable declarations in assignment are lost.  Example:
// double m; TIME_FUNCTION(m, min(a,b));
#define TIME_FUNCTION(assignment, function, ...) \
do {const size_t index = Timer::GetInstance().start(); \
assignment = function(__VA_ARGS__); \
Timer::GetInstance().stop(index, #function, __FILE__, __LINE__);} while (0)


#endif
