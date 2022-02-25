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
// Bar made of material with heat conductivity K = Ka*T^Kb,
// and Area-to-length ratio SL [m].
class LinkSimpleBar: public LinkBase {
    double Ka, Kb, SL;
  public:
    LinkSimpleBar(const double Ka, const double Kb, const double SL): Ka(Ka), Kb(Kb), SL(SL){;}

    double get_qdot(const double T1, const double T2) const override {
      // Q = SL * \int_T1^T2 K(T) dT
      return SL*Ka/(Kb+1) * (pow(T1, Kb+1) - pow(T2, Kb+1));
    }
};

/********************************************************************/
// Bar made of metal. Total resistance R.
class LinkMetalBar: public LinkBase {
    double R;
    double l = 2.44e-8;
  public:
    LinkMetalBar(const double R): R(R){;}

    double get_qdot(const double T1, const double T2) const override {
      // K = T*lambda/R
      return (T1*T1-T2*T2)/2 * l/R;
    }
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

  /*****************/
  // constant heat flow
  if (type=="const"){
    auto opts = get_key_val_args(b,e, {"type=", "Qdot=0W"});
    auto Q = read_power(opts["Qdot"]);
    return std::shared_ptr<LinkBase>(new LinkConst(Q));
  }

  /*****************/
  // A bar made of some material with Ka*T^Kb heat conductivity [W/m/K].
  // Length L, cross-section area S.
  if (type=="simple_bar"){
    auto opts = get_key_val_args(b,e, {"type=", "Ka=", "Kb=", "material=", "L=", "S="});

    double Ka=0, Kb=0, L=0, S=0;

    if      (opts["material"] == "torlon4203"){ Ka = 6.13e-3;  Kb=2.18; } // 0.1 - 5K
    else if (opts["material"] == "GRP")       { Ka = 2.42e-3;  Kb=1.77; } // 0.2 - 5K
    else if (opts["material"] == "nylon")     { Ka = 2.42e-3;  Kb=1.77; } // 0.2 - 5K
    else if (opts["material"] == "G10-CR")    { Ka = 2.523e-2; Kb=2.21; }
    else if (opts["material"] == "macor")     { Ka = 1.127e-3; Kb=2.82; }
    else if (opts["material"] == "stycast1266")    { Ka = 4.9e-2; Kb=1.98; } // 0.045-0.45K
    else if (opts["material"] == "stycast2850ft")  { Ka = 9.2e-3; Kb=2.65; }
    else if (opts["material"] == "araldite_ct200") { Ka = 2.4e-2; Kb=1.74; }
    else if (opts["material"] == "CuNi")           { Ka = 6.5e-2; Kb=1.1; }
    else if (opts["material"] == "Manganin")       { Ka = 5.976e-2; Kb=1.02; }

    if (opts["Ka"] != "") Ka = read_dimensionless(opts["Ka"]);
    if (opts["Kb"] != "") Ka = read_dimensionless(opts["Kb"]);
    if (opts["S"] != "") S = read_area(opts["S"]);
    if (opts["L"] != "") L = read_length(opts["L"]);
    if (Ka <= 0) throw Err() << "A positive value expected: Ka";
    if (L  <= 0) throw Err() << "A positive value expected: L";
    if (S  <= 0) throw Err() << "A positive value expected: S";
    return std::shared_ptr<LinkBase>(new LinkSimpleBar(Ka, Kb, S/L));
  }

  /*****************/
  // A bar made of metal. Total resistance R, Wiedemann-Franz low is used.
  if (type=="metal_bar"){
    auto opts = get_key_val_args(b,e, {"type=", "R="});
    double R = read_resistance(opts["R"]);
    return std::shared_ptr<LinkBase>(new LinkMetalBar(R));
  }

  throw Err() << "unknown link type: " << type;
}


#endif
