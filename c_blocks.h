#ifndef C_BLOCKS_H
#define C_BLOCKS_H

// Blocks

/********************************************************************/
// Base class for a thermal block
class BlockBase {
  // get heat capacity as a function of temperature and magnetic field
  public: virtual double get_heat_cap(const double T, const double H) const = 0;
};

/********************************************************************/
// Simple block with a constant heat capacity. This can be also used for
// a thermal bath (C=INFINITY)
class BlockSimple: public BlockBase {
    double C;
  public:
    BlockSimple(const double C): C(C){}
    double get_heat_cap(const double T, const double H) const override { return C;}
};

#endif
