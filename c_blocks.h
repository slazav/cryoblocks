#ifndef C_BLOCKS_H
#define C_BLOCKS_H

#include <memory>

#include "c_dim.h"

// Blocks
/********************************************************************/
/********************************************************************/
// Base class for a thermal block.
class BlockBase {
  public:

  // Calculate temperature change.
  // dQ = T dS = T(partial dS/dT) dT + T(partial dS/dB) dB
  virtual double get_dt(
    const double dQ, const double T, const double B, const double dB) const = 0;

  // is the block has zero heat capacity?
  virtual bool is_zero_c() const {return false;}
};

/********************************************************************/
// Simple block with a constant heat capacity. This can be also used for
// a thermal bath (C=INFINITY)
class BlockSimple: public BlockBase {
    double C;

  public:

    BlockSimple(const double C): C(C){}
    double get_dt(const double dQ, const double T, const double B, const double dB) const override {
       return dQ/C;
    }

    // simple block can have zero heat capacity
    bool is_zero_c() const {return C==0;}
};


// A block with paramegnetic heat capacity
class BlockParamagn: public BlockBase {

  double Bint; // internal field [T]
  double gyro; // gyromagnetic ratio [rad/s/T]
  double spin; // 1/2, 3/2, etc.
  double nmole; // number of moles

  double hbar = 1.05457181710e-34; // [J s]
  double NA = 6.02214086e+23;      // [1/mol] Avogadro's Constant
  double kB = 1.38064852e-23;      // [m^2 kg s-2 K-1] [J/K] Boltzmann constant
  double R = NA*kB;                // R-constant, [J/mol/K]

  public:
    BlockParamagn(const double Bint, const double gyro, const double spin, const double nmole):
      Bint(Bint), gyro(gyro), spin(spin), nmole(nmole){}

    double get_dt(const double dQ, const double T, const double B, const double dB) const override{
      double x = gyro*hbar*sqrt(B*B + Bint*Bint)/kB/T/2.0;
      double y = (2.0*spin + 1.0)*x;
      double C = R*nmole *(pow(x/sinh(x),2) - pow(y/sinh(y),2));
      // dT = dQ/C - D/C dB
      // We want to find D/C = (dS/dB)/(dS/dT). S depends only on A = T/sqrt(B*B + Bint*Bint)
      // Then D/C = (dA/dB) / (dA/dT) = - T*B/(B*B + Bint*Bint)
      double DC = T*B/(B*B + Bint*Bint);
      return dQ/C + DC*dB;
    }

};

/********************************************************************/
/********************************************************************/

// Create a block using vector<string> parameters:

std::shared_ptr<BlockBase> create_block(
  const std::vector<std::string>::const_iterator & b,
  const std::vector<std::string>::const_iterator & e){

  // extract type parameter
  auto type = get_key_val(b,e, "type");
  if (type == "") throw Err() << "block type is not set";

  // thermal bath, infinite heat capacity
  if (type == "bath") {
    auto opts = get_key_val_args(b,e, {"type="});
    return std::shared_ptr<BlockBase>(new BlockSimple(INFINITY));
  }

  // a block with zero heat capacity
  if (type == "zero-c") {
    auto opts = get_key_val_args(b,e, {"type="});
    return std::shared_ptr<BlockBase>(new BlockSimple(0));
  }

  // simple block with a constant heat capacity
  if (type == "simple") {
    auto opts = get_key_val_args(b,e, {"type=", "C=1J/K"});
    return std::shared_ptr<BlockBase>(new BlockSimple(read_heat_cap(opts["C"])));
  }

  // A block with paramagnetic heat capacity
  if (type == "paramagnet") {
    auto opts = get_key_val_args(b,e, {"type=", "Bint=", "gyro=", "spin=", "material=", "nmol="});

    double Bint=0, gyro=0, spin=0, nmol=0;

    if (opts["material"] == "copper"){
      Bint = 0.36e-3;    // [T], dipolar feld in copper
      gyro = 71.118e6;   // [rad/s/T] gyromagnetic ratio of copper
      spin = 1.5;        // spin 3/2
    }
    if (opts["material"] == "he3"){
      Bint = 720e-3;    // [T], dipolar feld in solid helium-3
      gyro = 203.789e6; // [rad/s/T] gyromagnetic ratio, helium-3
      spin = 0.5;       // spin 1/2, helium-3
    }
    if (opts["Bint"] != "") Bint = read_magn_field(opts["Bint"]);
    if (opts["gyro"] != "") gyro = read_gyro(opts["gyro"]);
    if (opts["spin"] != "") spin = read_dimensionless(opts["spin"]);
    if (opts["nmol"] != "") nmol = read_dimensionless(opts["nmol"]);
    if (Bint <= 0) throw Err() << "A positive value expected: Bint";
    if (gyro <= 0) throw Err() << "A positive value expected: gyro";
    if (spin <= 0) throw Err() << "A positive value expected: spin";
    if (nmol <= 0) throw Err() << "A positive value expected: nmol";
    return std::shared_ptr<BlockBase>(new BlockParamagn(Bint, gyro, spin, nmol));
  }
  throw Err() << "unknown block type: " << type;
}
#endif
