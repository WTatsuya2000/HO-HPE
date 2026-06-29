#include<string>
#include<sys/time.h>

class StopWatch
{
private:
  struct timeval tv_origin;
  double elapsed;
  bool is_stopping;

  double GetTime() const
  {
    struct timeval tv_tmp;
    gettimeofday( &tv_tmp , 0 );
    return (double)( tv_tmp.tv_sec - tv_origin.tv_sec )
      + 1.e-6 * (double)( tv_tmp.tv_usec - tv_origin.tv_usec );
  }

public:
  StopWatch(){ Reset(); }

  void Reset()
  {
    elapsed = 0.;
    gettimeofday( &tv_origin , 0 );
    is_stopping = false;
  }
  
  double operator()() const
  {
    if( is_stopping ) return elapsed;
    else              return GetTime() + elapsed;
  }

  string operator()( const char* head ) const
  {
    std::ostringstream message;
    if( is_stopping ) message << head << elapsed << "\n";
    else              message << head << GetTime() + elapsed << "\n";
    return message.str();
  }

  double Stop()
  {
    elapsed += GetTime();
    is_stopping = true;
    return elapsed;
  }

  double Start()
  {
    if( !is_stopping ){ Reset();  return 0.;  }
    is_stopping = false;
    gettimeofday( &tv_origin , 0 );
    return elapsed;
  }
} SWatch;

