#include <string>
#include <sstream>
#include "inc/err.h"

double
read_value(const std::string & str, const std::string & unit){
  std::istringstream ss(str);
  double v=0;
  std::string u;
  ss >> v >> u;

  if (u == unit) return v;

  // dimensionless
  if (unit == ""){
    throw Err() << "Dimensionless value expected: " << str;
  }

  // time
  else if (unit == "s"){
    if (u=="m") return 60*v;
    if (u=="h") return 3600*v;
    if (u=="d") return 24*3600*v;
    if (u=="ms") return 1e-3*v;
    if (u=="us") return 1e-6*v;
    if (u=="ns") return 1e-9*v;
  }

  // temperature
  else if (unit == "K"){
    if (u=="C") return v + 273.15;
    if (u=="mK") return 1e-3*v;
    if (u=="uK") return 1e-6*v;
    if (u=="nK") return 1e-9*v;
  }

  // length
  else if (unit == "m"){
    if (u=="cm") return 1e-2*v;
    if (u=="mm") return 1e-3*v;
    if (u=="um") return 1e-6*v;
    if (u=="nm") return 1e-9*v;
  }

  //area
  else if (unit == "m^2"){
    if (u=="cm^2") return 1e-4*v;
    if (u=="mm^2") return 1e-6*v;
    if (u=="um^2") return 1e-12*v;
    if (u=="nm^2") return 1e-18*v;
  }

  // volume
  else if (unit == "m^3"){
    if (u=="l")   return 1e-3*v;
    if (u=="cm^3") return 1e-6*v;
    if (u=="mm^3") return 1e-9*v;
    if (u=="um^3") return 1e-18*v;
    if (u=="nm^3") return 1e-27*v;
  }

  // mass
  else if (unit == "kg"){
    if (u=="g")   return v*1e-3;
    if (u=="mg") return v*1e-6;
  }

  // pressure
  else if (unit == "Pa"){
    if (u=="kPa")   return v*1e3;
    if (u=="MPa")   return v*1e6;
    if (u=="GPa")   return v*1e9;
    if (u=="bar")   return v*1e5;
    if (u=="mbar")   return v*1e2;
    if (u=="psi")   return v*6894.76;
    if (u=="atm")   return v*101325;
    if (u=="torr")  return v*133.322;
  }

  // energy
  else if (unit == "J"){
    if (u=="mJ") return v*1e-3;
    if (u=="uJ") return v*1e-6;
    if (u=="nJ") return v*1e-9;
    if (u=="pJ") return v*1e-12;
    if (u=="fJ") return v*1e-15;
  }

  // power
  else if (unit == "W"){
    if (u=="J/s")  return v;
    if (u=="mJ/s") return v*1e-3;
    if (u=="mW") return v*1e-3;
    if (u=="uW") return v*1e-6;
    if (u=="nW") return v*1e-9;
    if (u=="pW") return v*1e-12;
    if (u=="fW") return v*1e-15;
  }

  // heat capacity
  else if (unit == "J/K"){
    if (u=="mJ/K") return v*1e-3;
    if (u=="uJ/K") return v*1e-6;
    if (u=="nJ/K") return v*1e-9;
    if (u=="pJ/K") return v*1e-12;
    if (u=="fJ/K") return v*1e-15;
  }

  // thermal conductivity
  else if (unit == "W/K"){
    if (u=="mW/K") return v*1e-3;
    if (u=="uW/K") return v*1e-6;
    if (u=="nW/K") return v*1e-9;
    if (u=="pW/K") return v*1e-12;
    if (u=="fW/K") return v*1e-15;
  }

  // magnetic field
  else if (unit == "T"){
    if (u=="mT") return v*1e-3;
    if (u=="uT") return v*1e-6;
    if (u=="nT") return v*1e-9;
    if (u=="G")   return v*1e-4;
    if (u=="kG")  return v*1e-1;
    if (u=="Oe")  return v*1e-4;
    if (u=="kOe") return v*1e-1;
  }

  // magnetic field rate
  else if (unit == "T/s"){
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
    if (u=="T/h")   return v/3600.0;
    if (u=="mT/h")  return v*1e-3/3600.0;
    if (u=="uT/h")  return v*1e-6/3600.0;
    if (u=="nT/h")  return v*1e-9/3600.0;
    if (u=="G/h")   return v*1e-4/3600.0;
    if (u=="kG/h")  return v*1e-1/3600.0;
    if (u=="Oe/h")  return v*1e-4/3600.0;
    if (u=="kOe/h") return v*1e-1/3600.0;
  }

  // gyromagnetic ratio
  else if (unit == "rad/s/T"){
    if (u=="rad/s/G")  return v*1e4;
  }

  // resistance
  else if (unit == "Ohm"){
    if (u=="kOhm") return v*1e3;
    if (u=="MOhm") return v*1e6;
    if (u=="GOhm") return v*1e9;
    if (u=="mOhm") return v*1e-3;
    if (u=="uOhm") return v*1e-6;
    if (u=="nOhm") return v*1e-9;
  }

  // capacitance
  else if (unit == "F"){
    if (u=="mF") return v*1e-3;
    if (u=="uF") return v*1e-6;
    if (u=="nF") return v*1e-9;
    if (u=="pF") return v*1e-12;
    if (u=="fF") return v*1e-15;
  }

  // inductance
  else if (unit == "H"){
    if (u=="mH") return v*1e-3;
    if (u=="uH") return v*1e-6;
    if (u=="nH") return v*1e-9;
    if (u=="pH") return v*1e-12;
    if (u=="fH") return v*1e-15;
  }

  // Korringa constant
  else if (unit == "K*s"){
    if (u=="mK*s") return v*1e-3;
  }


  // molar flow
  else if (unit == "mol/s"){
    if (u=="mol/s") return v;
    if (u=="mmol/s") return v*1e-3;
    if (u=="umol/s") return v*1e-6;
  }

  throw Err() << "Can't convert to " << unit << ": " << str;
}
