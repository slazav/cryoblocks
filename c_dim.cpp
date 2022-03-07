#include <sstream>
#include "c_err.h"

double
read_dimensionless(const std::string & str){
  std::istringstream ss(str);
  double v=0;
  ss >> v;
  if (!ss || !ss.eof()) throw Err() << "Dimensionless value expected: " << str;
  return v;
}

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
read_length(const std::string & str){
  std::istringstream ss(str);
  double v=0;
  std::string u;
  ss >> v >> u;
  if (u=="m") return v;
  if (u=="cm") return 1e-2*v;
  if (u=="mm") return 1e-3*v;
  if (u=="um") return 1e-6*v;
  if (u=="nm") return 1e-9*v;
  throw Err() << "Unknown length unit: " << str;
}

double
read_area(const std::string & str){
  std::istringstream ss(str);
  double v=0;
  std::string u;
  ss >> v >> u;
  if (u=="m^2") return v;
  if (u=="cm^2") return 1e-4*v;
  if (u=="mm^2") return 1e-6*v;
  if (u=="um^2") return 1e-12*v;
  if (u=="nm^2") return 1e-18*v;
  throw Err() << "Unknown area unit: " << str;
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

double
read_gyro(const std::string & str){
  std::istringstream ss(str);
  double v=0;
  std::string u;
  ss >> v >> u;
  if (u=="rad/s/T") return v;
  if (u=="rad/s/G")  return v*1e4;
  throw Err() << "Unknown unit for gyromagnetic ratio: " << str;
}

double
read_resistance(const std::string & str){
  std::istringstream ss(str);
  double v=0;
  std::string u;
  ss >> v >> u;
  if (u=="Ohm")  return v;
  if (u=="kOhm") return v*1e3;
  if (u=="MOhm") return v*1e6;
  if (u=="GOhm") return v*1e9;
  if (u=="mOhm") return v*1e-3;
  if (u=="uOhm") return v*1e-6;
  if (u=="nOhm") return v*1e-9;
  throw Err() << "Unknown unit for gyromagnetic ratio: " << str;
}

double read_kappa(const std::string & str){
  std::istringstream ss(str);
  double v=0;
  std::string u;
  ss >> v >> u;
  if (u=="K*s")  return v;
  throw Err() << "Unknown unit for Korringa constant: " << str;
}

double
read_mass(const std::string & str){
  std::istringstream ss(str);
  double v=0;
  std::string u;
  ss >> v >> u;
  if (u=="kg")  return v;
  if (u=="g")   return v*1e-3;
  if (u=="mg") return v*1e-6;
  throw Err() << "Unknown unit for mass: " << str;
}

double
read_pressure(const std::string & str){
  std::istringstream ss(str);
  double v=0;
  std::string u;
  ss >> v >> u;
  if (u=="Pa")  return v;
  if (u=="kPa")   return v*1e3;
  if (u=="MPa")   return v*1e6;
  if (u=="GPa")   return v*1e9;
  if (u=="bar")   return v*1e5;
  if (u=="mbar")   return v*1e2;
  if (u=="psi")   return v*6894.76;
  if (u=="atm")   return v*101325;
  if (u=="torr")  return v*133.322;
  throw Err() << "Unknown unit for pressure: " << str;
}

