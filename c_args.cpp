#include <iostream>
#include <sstream>
#include "c_err.h"
#include "c_args.h"

std::ostream &
operator<< (std::ostream & s, const std::vector<std::string> & v){
  for (int i=0; i<v.size(); i++) s << (i>0? ", ":"") << v[i];
  return s;
}

Opt
get_position_args(const std::vector<std::string>::const_iterator & b,
                  const std::vector<std::string>::const_iterator & e,
                  const std::vector<std::string> & arg_list){

  // Parse position arguments
  Opt opts;
  auto arg = b;
  for (auto const & key:arg_list){
    if (arg==e)
      throw Err() << "Position argument is missing (expected: " << arg_list << ")";
    opts[key] = *arg;
    ++arg;
  }
  return opts;
}

Opt
get_key_val_args(const std::vector<std::string>::const_iterator & b,
                 const std::vector<std::string>::const_iterator & e,
                 const std::vector<std::string> & arg_list){

  // Extract default values of key-value arguments
  Opt opts;
  for (auto const & a:arg_list){
    // parse key=value
    auto n = a.find('=', 0);
    if (n == std::string::npos)
      throw Err() << "Can't find key=value pair in arg_list: " << a;
    auto key = a.substr(0,n);
    auto val = a.substr(n+1);
    opts[key] = val;
  }

  // Parse key-value arguments
  for (auto arg=b; arg!=e; ++arg){
    // parse key=value
    auto n = arg->find('=', 0);
    if (n == std::string::npos)
      throw Err() << "Can't find key=value pair: " << *arg;
    auto key = arg->substr(0,n);
    auto val = arg->substr(n+1);

    // replace existing value in opts
    if (opts.count(key)==0)
      throw Err() << "Unknown parameter: " << key << " (known parameters with defaults: " << arg_list << ")";
    opts[key] = val;
  }
  return opts;
}

// get a single argument
std::string
get_key_val(const std::vector<std::string>::const_iterator & b,
            const std::vector<std::string>::const_iterator & e,
            const std::string & key, const std::string & def){

  std::string ret = def;
  for (auto arg=b; arg!=e; ++arg){
    // parse key=value
    auto n = arg->find('=', 0);
    if (n == std::string::npos) continue;
    if (key != arg->substr(0,n)) continue;
    ret = arg->substr(n+1);
  }
  return ret;
}

