#include <sstream>
#include "c_err.h"

double
read_time(const std::string & str){
  std::istringstream ss(str);
  double v=0;
  std::string u;
  ss >> v >> u;
  if (u=="s") return v;
  if (u=="m") return 60*v;
  if (u=="h") return 3600*v;
  if (u=="d") return 24*3600*v;
  if (u=="ms") return 1e-3*v;
  if (u=="us") return 1e-6*v;
  if (u=="ns") return 1e-9*v;
  throw Err() << "Unknown time unit: " << str;
}

double
read_temp(const std::string & str){
  std::istringstream ss(str);
  double v=0;
  std::string u;
  ss >> v >> u;
  if (u=="C") return v + 273.15;
  if (u=="K") return v;
  if (u=="mK") return 1e-3*v;
  if (u=="uK") return 1e-6*v;
  if (u=="nK") return 1e-9*v;
  throw Err() << "Unknown temperature unit: " << str;
}

double
read_heat_cap(const std::string & str){
  std::istringstream ss(str);
  double v=0;
  std::string u;
  ss >> v >> u;
  if (u=="J/K")  return v;
  if (u=="mJ/K") return v*1e-3;
  if (u=="uJ/K") return v*1e-6;
  if (u=="nJ/K") return v*1e-9;
  if (u=="pJ/K") return v*1e-12;
  if (u=="fJ/K") return v*1e-15;
  throw Err() << "Unknown heat capacity unit: " << str;
}

double
read_power(const std::string & str){
  std::istringstream ss(str);
  double v=0;
  std::string u;
  ss >> v >> u;
  if (u=="W")  return v;
  if (u=="mW") return v*1e-3;
  if (u=="uW") return v*1e-6;
  if (u=="nW") return v*1e-9;
  if (u=="pW") return v*1e-12;
  if (u=="fW") return v*1e-15;
  throw Err() << "Unknown power unit: " << str;
}

double
read_heat_cond(const std::string & str){
  std::istringstream ss(str);
  double v=0;
  std::string u;
  ss >> v >> u;
  if (u=="W/K")  return v;
  if (u=="mW/K") return v*1e-3;
  if (u=="uW/K") return v*1e-6;
  if (u=="nW/K") return v*1e-9;
  if (u=="pW/K") return v*1e-12;
  if (u=="fW/K") return v*1e-15;
  throw Err() << "Unknown heat conductivity unit: " << str;
}

double
read_magn_field(const std::string & str){
  std::istringstream ss(str);
  double v=0;
  std::string u;
  ss >> v >> u;
  if (u=="T")  return v;
  if (u=="mT") return v*1e-3;
  if (u=="uT") return v*1e-6;
  if (u=="nT") return v*1e-9;
  if (u=="G")   return v*1e-4;
  if (u=="kG")  return v*1e-1;
  if (u=="Oe")  return v*1e-4;
  if (u=="kOe") return v*1e-1;
  throw Err() << "Unknown magnetic field unit: " << str;
}

double
read_magn_field_rate(const std::string & str){
  std::istringstream ss(str);
  double v=0;
  std::string u;
  ss >> v >> u;
  if (u=="T/s")   return v;
  if (u=="mT/s")  return v*1e-3;
  if (u=="uT/s")  return v*1e-6;
  if (u=="nT/s")  return v*1e-9;
  if (u=="G/s")   return v*1e-4;
  if (u=="kG/s")  return v*1e-1;
  if (u=="Oe/s")  return v*1e-4;
  if (u=="kOe/s") return v*1e-1;
  if (u=="T/m")   return v/60.0;
  if (u=="mT/m")  return v*1e-3/60.0;
  if (u=="uT/m")  return v*1e-6/60.0;
  if (u=="nT/m")  return v*1e-9/60.0;
  if (u=="G/m")   return v*1e-4/60.0;
  if (u=="kG/m")  return v*1e-1/60.0;
  if (u=="Oe/m")  return v*1e-4/60.0;
  if (u=="kOe/m") return v*1e-1/60.0;
  throw Err() << "Unknown magnetic field rate unit: " << str;
}
