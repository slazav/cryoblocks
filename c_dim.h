#ifndef C_READ_DIM_H
#define C_READ_DIM_H

// Read parameters with dimensions

double read_dimensionless(const std::string & str);
double read_time(const std::string & str);
double read_temp(const std::string & str);
double read_length(const std::string & str);
double read_area(const std::string & str);
double read_heat_cap(const std::string & str);
double read_power(const std::string & str);
double read_heat_cond(const std::string & str);
double read_magn_field(const std::string & str);
double read_magn_field_rate(const std::string & str);
double read_gyro(const std::string & str);
double read_resistance(const std::string & str);

#endif
