#ifndef C_READ_DIM_H
#define C_READ_DIM_H

#include <string>

// Read parameters with dimensions

double read_dimensionless(const std::string & str);
double read_time(const std::string & str);
double read_temp(const std::string & str);
double read_length(const std::string & str);
double read_area(const std::string & str);
double read_volume(const std::string & str);
double read_heat_cap(const std::string & str);
double read_power(const std::string & str);
double read_heat_cond(const std::string & str);
double read_magn_field(const std::string & str);
double read_magn_field_rate(const std::string & str);
double read_gyro(const std::string & str);
double read_resistance(const std::string & str);
double read_kappa(const std::string & str); // Karringa constant K*s
double read_mass(const std::string & str);
double read_pressure(const std::string & str);
double read_circ(const std::string & str);

#endif
