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

  typedef std::vector<std::string>::const_iterator str_cit;

  // Calculate temperature change.
  // dQ = T dS = T(partial dS/dT) dT + T(partial dS/dB) dB
  virtual double get_dt(
    const double dQ, const double T, const double B, const double dB) const = 0;

  // is the block has zero heat capacity?
  virtual bool is_zero_c() const {return false;}
};

/********************************************************************/
// Simple block with a constant heat capacity. This can be also used for
// a thermal bath (C=INFINITY) or zero-c blocks
class BlockSimple: public BlockBase {
    double C;

  public:

    static std::shared_ptr<BlockBase> create(const str_cit & b, const str_cit & e){
      auto opts = get_key_val_args(b,e, {"type=", "C=1J/K"});
      return std::shared_ptr<BlockBase>(new BlockSimple(read_heat_cap(opts["C"])));
    }

    static std::shared_ptr<BlockBase> create_bath(const str_cit & b, const str_cit & e){
      auto opts = get_key_val_args(b,e, {"type="});
      return std::shared_ptr<BlockBase>(new BlockSimple(INFINITY));
    }

    static std::shared_ptr<BlockBase> create_zero(const str_cit & b, const str_cit & e){
      auto opts = get_key_val_args(b,e, {"type="});
      return std::shared_ptr<BlockBase>(new BlockSimple(0));
    }

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
  double moles; // number of moles

  double hbar = 1.05457181710e-34; // [J s]
  double NA = 6.02214086e+23;      // [1/mol] Avogadro's Constant
  double kB = 1.38064852e-23;      // [m^2 kg s-2 K-1] [J/K] Boltzmann constant
  double R = NA*kB;                // R-constant, [J/mol/K]

  public:

    static std::shared_ptr<BlockBase> create(const str_cit & b, const str_cit & e){
      auto opts = get_key_val_args(b,e, {"type=", "Bint=", "gyro=", "spin=", "material=", "moles=", "mass="});
      double Bint=0, gyro=0, spin=0, moles=0;

      if (opts["material"] == "copper"){
        Bint = 0.36e-3;    // [T], dipolar feld in copper
        gyro = 71.118e6;   // [rad/s/T] gyromagnetic ratio of copper
        spin = 1.5;        // spin 3/2
        if (opts["mass"]!="") moles = read_mass(opts["mass"])/63.546e-3;
      }
      else if (opts["material"] == "he3"){
        Bint = 720e-3;    // [T], dipolar feld in solid helium-3
        gyro = 203.789e6; // [rad/s/T] gyromagnetic ratio, helium-3
        spin = 0.5;       // spin 1/2, helium-3
        if (opts["mass"]!="") moles = read_mass(opts["mass"])/3.016029e-3;
      }
      else {
        if (opts["mass"]!="") throw Err() << "mass parameter can be used only together with material";
      }
      if (opts["Bint"]  != "") Bint  = read_magn_field(opts["Bint"]);
      if (opts["gyro"]  != "") gyro  = read_gyro(opts["gyro"]);
      if (opts["spin"]  != "") spin  = read_dimensionless(opts["spin"]);
      if (opts["moles"] != "") moles = read_dimensionless(opts["moles"]);
      if (Bint  <= 0) throw Err() << "A positive value expected: Bint";
      if (gyro  <= 0) throw Err() << "A positive value expected: gyro";
      if (spin  <= 0) throw Err() << "A positive value expected: spin";
      if (moles <= 0) throw Err() << "A positive value expected: moles";
      return std::shared_ptr<BlockBase>(new BlockParamagn(Bint, gyro, spin, moles));
    }

    BlockParamagn(const double Bint, const double gyro, const double spin, const double moles):
      Bint(Bint), gyro(gyro), spin(spin), moles(moles){}

    double get_dt(const double dQ, const double T, const double B, const double dB) const override{
      double x = gyro*hbar*sqrt(B*B + Bint*Bint)/kB/T/2.0;
      double y = (2.0*spin + 1.0)*x;
      double C = R*moles *(pow(x/sinh(x),2) - pow(y/sinh(y),2));
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
  if (type == "bath") { return BlockSimple::create_bath(b,e);}
  if (type == "zero-c") { return BlockSimple::create_zero(b,e);}
  if (type == "simple") { return BlockSimple::create(b,e);}
  if (type == "paramagnet") { return BlockParamagn::create(b,e); }

  throw Err() << "unknown block type: " << type;
}
#endif
