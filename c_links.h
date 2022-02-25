#ifndef C_LINKS_H
#define C_LINKS_H

#include <memory>

#include "c_dim.h"

/********************************************************************/
// Base class for a thermal link
class LinkBase {
  public:
    // get heat flow (in W) from block1 to block2
    virtual double get_qdot(const double T1, const double T2) const = 0;
};

/********************************************************************/
// Constant energy flow which does not depend on temperatures.
// Not a physical process, but good for tests.
class LinkConst: public LinkBase {
    double Q;
  public:
    LinkConst(const double Q): Q(Q){;}

    double get_qdot(const double T1, const double T2) const override { return Q;}
};

/********************************************************************/
/********************************************************************/

// Create a link using vector<string> parameters:

std::shared_ptr<LinkBase> create_link(
  const std::vector<std::string>::const_iterator & b,
  const std::vector<std::string>::const_iterator & e){

  // extract type
  auto type = get_key_val(b,e, "type");
  if (type == "")   throw Err() << "link type is not set";

  if (type=="const"){
    auto opts = get_key_val_args(b,e, {"type=", "Qdot=0W"});
    auto Q = read_power(opts["Qdot"]);
    return std::shared_ptr<LinkBase>(new LinkConst(Q));
  }
  throw Err() << "unknown link type: " << type;
}


#endif
