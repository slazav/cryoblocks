#ifndef C_LINKS_H
#define C_LINKS_H

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


#endif
