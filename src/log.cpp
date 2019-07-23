#include "log.h"

Log::Log(char* filename) {
  m_stream.open(filename);
}

void Log::write(char* logline){
  m_stream << logline << endl;
}

Log::~Log(){
  m_stream.close();
}
