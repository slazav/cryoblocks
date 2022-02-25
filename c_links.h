#ifndef C_LINKS_H
#define C_LINKS_H

#include <memory>

#include "c_dim.h"

/********************************************************************/
// Base class for a thermal link
class LinkBase {
    std::string block1, block2;
  public:

    LinkBase(const std::string & bl1, const std::string & bl2): block1(bl1), block2(bl2){}

    // get name of the first block
    const std::string & get_block1() const {return block1;}

    // get name of the second block
    const std::string & get_block2() const {return block2;}

    // get heat flow (in W) from block1 to block2
    virtual double get_qdot(const double T1, const double T2) const = 0;
};

/********************************************************************/
// Constant energy flow which does not depend on temperatures.
// Not a physical process, but good for tests.
class LinkConst: public LinkBase {
    double Q;
  public:
    LinkConst(const std::string & bl1, const std::string & bl2, const double Q):
      LinkBase(bl1,bl2), Q(Q){;}

    double get_qdot(const double T1, const double T2) const override { return Q;}
};

/********************************************************************/
/********************************************************************/

// Create a link using vector<string> parameters:

std::shared_ptr<LinkBase> create_link(
  const std::vector<std::string>::const_iterator & b,
  const std::vector<std::string>::const_iterator & e,
  const std::string & bl1, const std::string & bl2){

  // extract type
  auto type = get_key_val(b,e, "type");
  if (type == "")   throw Err() << "link type is not set";

  if (type=="const"){
    auto opts = get_key_val_args(b,e, {"type=", "Qdot=0W"});
    auto Q = read_power(opts["Qdot"]);
    return std::shared_ptr<LinkBase>(new LinkConst(bl1, bl2, Q));
  }
  throw Err() << "unknown link type: " << type;
}


#endif
