#ifndef C_ERR_H
#define C_ERR_H

// Error class for exceptions
class Err {
  std::ostringstream s;
  public:
    Err(){}
    Err(const Err & o) { s << o.s.str(); }
    template <typename T>
    Err & operator<<(const T & o){ s << o; return *this; }
    std::string str()  const { return s.str(); }
};

#endif
