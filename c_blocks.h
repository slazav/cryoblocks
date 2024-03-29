#ifndef C_BLOCKS_H
#define C_BLOCKS_H

#include <memory>
#include <cmath>
#include <algorithm>

#include "c_dim.h"

// Blocks
/********************************************************************/
/********************************************************************/
// Base class for a thermal block.
class BlockBase {
  public:

  typedef std::vector<std::string>::const_iterator str_cit;

  // Get heat capacity, C = T*(partial dS/dT)
  // as a function of temperature and magnetic field
  virtual double get_C(const double T, const double B) const = 0;

  // Get demag cooling power, D = T*(partial dS/dB) / C = (partial dS/dB)/(partial dS/dT)
  // as a function of temperature and magnetic field
  virtual double get_D(const double T, const double B) const = 0;

  // Calculate temperature change.
  // dQ = T dS = T(partial dS/dT) dT + T(partial dS/dB) dB = C dT + C*D dB
  // dT = dQ/C - D*dB
  double get_dt(const double dQ, const double T, const double B, const double dB) const{
    return dQ/get_C(T,B) - get_D(T,B)*dB;
  }

  // is the block has zero heat capacity?
  virtual bool is_zero_c() const {return false;}
};

std::shared_ptr<BlockBase> create_block(
  const std::vector<std::string>::const_iterator & b,
  const std::vector<std::string>::const_iterator & e);

/********************************************************************/
// Compound block.
class BlockComp: public BlockBase {
  std::vector< std::shared_ptr<BlockBase> > blocks;

  public:

  // Constructor.
  // Arguments for different blocks are separated with "+"
  BlockComp(const str_cit & b, const str_cit & e){
    auto b1=std::find(b, e, ":")+1;
    while (1) {
      auto e1 = std::find(b1, e, "+");
      auto bl = create_block(b1,e1);
      if (!bl->is_zero_c()) blocks.push_back(bl);
      if (e1==e) break;
      b1 = e1 + 1;
    }
  }

  // Heat capacity is a sum of blocks' heat capacities
  double get_C(const double T, const double B) const override {
    double C = 0.0;
    for (const auto bl: blocks) C += bl->get_C(T,B);
    return C;
  }

  // C*D is additive!
  virtual double get_D(const double T, const double B) const override {
    double C  = 0.0;
    double CD = 0.0;
    for (const auto bl: blocks){
      C += bl->get_C(T,B);
      CD += C*bl->get_D(T,B);
    }
    if (C==0.0) return 0.0;
    return CD/C;
  }

  bool is_zero_c() const override {return blocks.size()==0;}
};

/********************************************************************/
// Simple block with a constant or power function heat capacity.
// Can have zero heat capacity.
class BlockSimple: public BlockBase {
    double factor=0;
    double power=0;

  public:

    BlockSimple(const str_cit & b, const str_cit & e){
      auto opts = get_key_val_args(b,e, {"type=", "C=", "factor=", "power="});
      if (opts["C"]     != "") {factor = read_value(opts["C"], "J/K"); power=0;}
      if (opts["factor"] != "") factor = read_value(opts["factor"], "");
      if (opts["power"]  != "") power = read_value(opts["power"],  "");
    }

    double get_C(const double T, const double B) const override {
      return factor*pow(T,power);
    }

    double get_D(const double T, const double B) const override {
      return 0.0;
    }

    // simple block can have zero heat capacity
    bool is_zero_c() const override {return factor==0;}
};

/********************************************************************/
//  Thermal bath, block with infinite heat capacity.
class BlockBath: public BlockBase {
  public:

    BlockBath(const str_cit & b, const str_cit & e){
      auto opts = get_key_val_args(b,e, {"type="}); }

    double get_C(const double T, const double B) const override {
      return INFINITY;
    }

    double get_D(const double T, const double B) const override {
      return 0.0;
    }

};

/********************************************************************/
//  block with zero heat capacity.
class BlockZero: public BlockBase {
  public:

    BlockZero(const str_cit & b, const str_cit & e){
      auto opts = get_key_val_args(b,e, {"type="}); }

    double get_C(const double T, const double B) const override {
      return 0.0;
    }

    double get_D(const double T, const double B) const override {
      return 0.0;
    }

    bool is_zero_c() const override {return true;}
};

#ifdef HE3
extern "C" {
#include <he3.h>
}

/********************************************************************/
// Paramagnetic material
class BlockParamagn: public BlockBase {

  double Bint=0; // internal field [T]
  double gyro=0; // gyromagnetic ratio [rad/s/T]
  double spin=0; // 1/2, 3/2, etc.
  double moles=0; // number of moles
  double R = 8.314472; // R-constant, [J/mol/K]
  double Bf = 1; // factor for external magnetic field

  public:

    BlockParamagn(const str_cit & b, const str_cit & e){
      auto opts = get_key_val_args(b,e, {"type=", "Bint=", "gyro=", "spin=", "material=", "moles=", "mass=", "Bfactor="});

      if (opts["material"] == "copper"){
        Bint = 0.36e-3;    // [T], dipolar feld in copper
        gyro = 71.118e6;   // [rad/s/T] gyromagnetic ratio of copper
        spin = 1.5;        // spin 3/2
        if (opts["mass"]!="") moles = read_value(opts["mass"], "kg")/63.546e-3;
      }
      else if (opts["material"] == "he3"){
        Bint = 720e-3;    // [T], dipolar feld in solid helium-3
        gyro = 203.789e6; // [rad/s/T] gyromagnetic ratio, helium-3
        spin = 0.5;       // spin 1/2, helium-3
        if (opts["mass"]!="") moles = read_value(opts["mass"], "kg")/3.016029e-3;
      }
      else {
        if (opts["mass"]!="") throw Err() << "mass parameter can be used only together with material";
      }
      if (opts["Bint"]  != "") Bint  = read_value(opts["Bint"], "T");

      if (opts["Bfactor"]  != "") Bf  = read_value(opts["Bfactor"], "");

      if (opts["gyro"]  != "") gyro  = read_value(opts["gyro"], "rad/s/T");
      if (gyro  <= 0) throw Err() << "A positive value expected: gyro";

      if (opts["spin"]  != "") spin  = read_value(opts["spin"], "");
      if (spin  <= 0) throw Err() << "A positive value expected: spin";

      if (opts["moles"] != "") moles = read_value(opts["moles"], "");
      if (moles <= 0) throw Err() << "A positive value expected: moles";
    }

    double get_C(const double T, const double B) const override {
      double t(T), b(B*Bf), bi(Bint), g(gyro), s(spin);
      return moles*R*magn_par_c_(&t,&b,&bi,&g,&s);
    }

    double get_D(const double T, const double B) const override {
      double t(T), b(B*Bf), bi(Bint), g(gyro), s(spin);
      return magn_par_d_(&t,&b,&bi,&g,&s);
    }

};

/********************************************************************/
// Curie-Weiss magnet with spin 1/2
class BlockCurieWeiss: public BlockBase {

  double Tc=0;   // critical temp [K]
  double gyro=0; // gyromagnetic ratio [rad/s/T]
  double moles=0; // number of moles
  double R = 8.314472; // R-constant, [J/mol/K]
  double Bf = 1;

  public:

    BlockCurieWeiss(const str_cit & b, const str_cit & e){
      auto opts = get_key_val_args(b,e, {"type=", "Tc=", "gyro=", "material=", "moles=", "mass=", "Bfactor="});

      if (opts["Tc"]  != "") Tc = read_value(opts["Tc"], "K");
      if (Tc <= 0) throw Err() << "A positive value expected: Tc";

      if (opts["gyro"]  != "") gyro  = read_value(opts["gyro"], "rad/s/T");
      if (gyro  <= 0) throw Err() << "A positive value expected: gyro";

      if (opts["moles"] != "") moles = read_value(opts["moles"], "");
      if (moles <= 0) throw Err() << "A positive value expected: moles";

      if (opts["Bfactor"]  != "") Bf  = read_value(opts["Bfactor"], "");
    }

    double get_C(const double T, const double B) const override {
      double t(T), b(B*Bf), tc(Tc), g(gyro);
      return moles*R*magn_cw_c_(&t,&b,&tc,&g);
    }

    double get_D(const double T, const double B) const override {
      double t(T), b(B*Bf), tc(Tc), g(gyro);
      return magn_cw_d_(&t,&b,&tc,&g);
    }

};


/********************************************************************/
// Liquid He3
// No A-phase support!
class BlockLHe3: public BlockBase {

  double P=0;  // pressure [bar]
  double Tc;   // Tc, K
  double moles=0;
  double R = 8.314472; // R-constant, [J/mol/K]

  public:

    BlockLHe3(const str_cit & b, const str_cit & e){
      auto opts = get_key_val_args(b,e, {"type=", "P=", "moles=", "mass=", "volume="});

      if (opts["P"]       != "") P     = read_value(opts["P"], "Pa")/1e5; // Pa -> bar
      if (opts["mass"]    != "") moles = read_value(opts["mass"], "kg") / 3.016029e-3;
      if (opts["volume"]  != "") moles = read_value(opts["volume"], "m^3")*1e6 / he3_vm_(&P);

      if (opts["moles"]   != "") moles = read_value(opts["moles"], "");
      if (moles <= 0) throw Err() << "A positive value expected: moles";
      Tc = he3_tc_(&P)*1e-3;
    }

    // no support for A-phase!
    double get_C(const double T, const double B) const override {
      double p(P), t(T), ttc(T/Tc);
      if (T > Tc){
         double Vm = he3_vm_(&p);
         return he3_cv_n_(&t,&Vm)*R*moles;
      }
      else {
        double ttc = T/Tc;
        return he3_c_b_(&ttc,&p)*R*moles;
      }
    }

    double get_D(const double T, const double B) const override {
      return 0.0;
    }

};

#endif


/********************************************************************/
/********************************************************************/

// Create a block using vector<string> parameters:

std::shared_ptr<BlockBase> create_block(
  const std::vector<std::string>::const_iterator & b,
  const std::vector<std::string>::const_iterator & e){

  // extract type parameter
  auto type = get_key_val(b,e, "type");
  if (type == "") throw Err() << "block type is not set";

  if (type == "bath") {
    return std::shared_ptr<BlockBase>(new BlockBath(b,e));}
  if (type == "zero-c") {
    return std::shared_ptr<BlockBase>(new BlockZero(b,e));}
  if (type == "simple") {
    return std::shared_ptr<BlockBase>(new BlockSimple(b,e));}
  if (type == "compound") {
    return std::shared_ptr<BlockBase>(new BlockComp(b,e));}
#ifdef HE3
  if (type == "paramagnet") {
    return std::shared_ptr<BlockBase>(new BlockParamagn(b,e)); }
  if (type == "curie-weiss") {
    return std::shared_ptr<BlockBase>(new BlockCurieWeiss(b,e)); }
  if (type == "liquid_he3") {
    return std::shared_ptr<BlockBase>(new BlockLHe3(b,e)); }
#endif

  throw Err() << "unknown block type: " << type;
}

std::shared_ptr<BlockBase> create_block(const std::vector<std::string> & v) {
  return create_block(v.begin(), v.end());}

#endif
