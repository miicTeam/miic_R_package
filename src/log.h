#ifndef _LOG_H_
#define _LOG_H_

#include <fstream>
using namespace std;

class Log {
  public:
    Log(char* filename);
    ~Log();
    void write(char* logline);
  private:
    ofstream m_stream;
};

#endif
