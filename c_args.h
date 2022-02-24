#ifndef C_ARGS_H
#define C_ARGS_H

#include <vector>
#include <string>
#include <map>

/******************************************************************/
// Options (key-value pairs)
class Opt: public std::map<std::string, std::string> {
};

std::ostream & operator<< (std::ostream & s, const std::vector<std::string> & v);

Opt
get_position_args(const std::vector<std::string>::const_iterator & b,
                  const std::vector<std::string>::const_iterator & e,
                  const std::vector<std::string> & arg_list);

Opt
get_key_val_args(const std::vector<std::string>::const_iterator & b,
                 const std::vector<std::string>::const_iterator & e,
                 const std::vector<std::string> & arg_list);

std::string
get_key_val(const std::vector<std::string>::const_iterator & b,
            const std::vector<std::string>::const_iterator & e,
            const std::string & key, const std::string & def = std::string());

#endif
