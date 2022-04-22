#ifndef C_LINKS_H
#define C_LINKS_H

#include <memory>

#include "c_dim.h"

/********************************************************************/
// Base class for a thermal link
class LinkBase {
  public:
    typedef std::vector<std::string>::const_iterator str_cit;

    // get heat flow (in W) from block1 to block2
    virtual double get_qdot(const double T1, const double T2, const double B) const = 0;

};

/********************************************************************/
// Constant energy flow which does not depend on temperatures.
// Not a physical process, but good for tests.
class LinkConst: public LinkBase {
    double Q = 0;
  public:

    /*****************/
    LinkConst(const str_cit & b, const str_cit & e) {
      auto opts = get_key_val_args(b,e, {"type=", "Qdot="});
      if (opts["Qdot"]!="") Q = read_value(opts["Qdot"], "W");
    }

    /*****************/
    double get_qdot(const double T1, const double T2, const double B) const override { return Q;}

};

/********************************************************************/
// Energy flow proportional to temperature difference, Qdot=K*(T1-T2).
class LinkSimple: public LinkBase {
    double K = 0;
  public:

    /*****************/
    LinkSimple(const str_cit & b, const str_cit & e) {
      auto opts = get_key_val_args(b,e, {"type=", "K="});
      if (opts["K"]!="") K = read_value(opts["K"], "W/K");
    }

    /*****************/
    double get_qdot(const double T1, const double T2, const double B) const override { return K*(T1-T2);}
};

/********************************************************************/
// Bar made of material with heat conductivity K = Ka*T^Kb,
// and area-to-length ratio SL [m].
class LinkSimpleBar: public LinkBase {
    double Ka=0, Kb=0, SL=0;
  public:

    /*****************/
    LinkSimpleBar(const str_cit & b, const str_cit & e) {
      auto opts = get_key_val_args(b,e, {"type=", "Ka=", "Kb=", "material=", "L=", "S="});

      double L=0, S=0;

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

      if (opts["Ka"] != "") Ka = read_value(opts["Ka"], "");
      if (Ka <= 0) throw Err() << "A positive value expected: Ka";

      if (opts["Kb"] != "") Kb = read_value(opts["Kb"], "");

      if (opts["S"] != "") S = read_value(opts["S"], "m^2");
      if (S  <= 0) throw Err() << "A positive value expected: S";

      if (opts["L"] != "") L = read_value(opts["L"], "m");
      if (L  <= 0) throw Err() << "A positive value expected: L";

      SL = S/L;
    }

    /*****************/
    double get_qdot(const double T1, const double T2, const double B) const override {
      // Q = SL * \int_T1^T2 K(T) dT
      return SL*Ka/(Kb+1) * (pow(T1, Kb+1) - pow(T2, Kb+1));
    }
};

/********************************************************************/
// Bar made of metal. Total resistance R.
class LinkMetalBar: public LinkBase {
    double R=0;
    double l = 2.44e-8;
  public:

    /*****************/
    LinkMetalBar(const str_cit & b, const str_cit & e) {
      auto opts = get_key_val_args(b,e, {"type=", "R="});
      if (opts["R"] != "") R = read_value(opts["R"], "Ohm");
      if (R <= 0) throw Err() << "A positive value expected: R";
    }

    /*****************/
    double get_qdot(const double T1, const double T2, const double B) const override {
      // K = T*lambda/R
      return (T1*T1-T2*T2)/2 * l/R;
    }
};


/********************************************************************/
// Spin-lattice relaxation, Korringa law.
// block1 should be spin system, block2 - electrons
class LinkKorringa: public LinkBase {
  double Bint=0; // internal field [T]
  double gyro=0; // gyromagnetic ratio [rad/s/T]
  double spin=0; // 1/2, 3/2, etc.
  double kappa0=0; // high-field value of Korringa constant
  double alpha=1;  // alpha parameter in field dependence of kappa
  double moles=0;  // number of moles

  double hbar = 1.05457181710e-34; // [J s]
  double NA = 6.02214086e+23;      // [1/mol] Avogadro's Constant
  double kB = 1.38064852e-23;      // [m^2 kg s-2 K-1] [J/K] Boltzmann constant
  double R = NA*kB;                // R-constant, [J/mol/K]
  double mu0 = 4*M_PI*1e-7;        // [V s/A/m] Vacuum permeability
  double muN = 5.05e-27;           // [A m^2] nuclear magneton

  public:

    /*****************/
    LinkKorringa(const str_cit & b, const str_cit & e) {
      auto opts = get_key_val_args(b,e,
        {"type=", "Bint=", "gyro=", "spin=", "kappa0=", "alpha=", "moles=", "material=", "mass="});

      if (opts["material"] == "copper"){
        Bint = 0.36e-3;    // [T], dipolar feld in copper
        gyro = 71.118e6;   // [rad/s/T] gyromagnetic ratio of copper
        spin = 1.5;        // spin 3/2
        kappa0 = 1.2;
        alpha  = 2.6;
        if (opts["mass"]!="") moles = read_value(opts["mass"], "kg")/63.546e-3;
      }
      else {
        if (opts["mass"]!="") throw Err() << "mass parameter can be used only together with material";
      }

      if (opts["Bint"]   != "") Bint = read_value(opts["Bint"], "T");

      if (opts["gyro"]   != "") gyro = read_value(opts["gyro"], "rad/s/T");
      if (gyro   <= 0) throw Err() << "A positive value expected: gyro";

      if (opts["spin"]   != "") spin = read_value(opts["spin"], "");
      if (spin   <= 0) throw Err() << "A positive value expected: spin";

      if (opts["kappa0"] != "") kappa0 = read_value(opts["kappa0"], "K*s");
      if (kappa0 <= 0) throw Err() << "A positive value expected: kappa0";

      if (opts["alpha"]  != "") alpha  = read_value(opts["alpha"], "");

      if (opts["moles"]  != "") moles = read_value(opts["moles"], "");
      if (moles  <= 0) throw Err() << "A positive value expected: moles";
    }

    /*****************/
    double get_qdot(const double T1, const double T2, const double B) const override {
      // heat capacity of copper nuclei:
      double x = gyro*hbar*sqrt(B*B + Bint*Bint)/kB/T1/2.0;
      double y = (2.0*spin + 1.0)*x;
      double C = R*moles *(pow(x/sinh(x),2) - pow(y/sinh(y),2));

      // See Pobell book f.10.10 and below
      double kappa = kappa0 * (B*B + Bint*Bint)/(B*B + alpha*Bint*Bint);
      double z = (mu0*muN*B)/(2*kB);
      double tau = kappa/z * tanh(z/T2);
      return C * T1/T2*(T1-T2) / tau;
    }

};

/********************************************************************/
// Electron-phonon coupling  Qdot = C*(Tph^5-Tel^5)
// C = 2e3 * 7.11 cm^3/mole * Nmoles  for copper
// See Pobell book, f.10.9
class LinkElPh: public LinkBase {
  double CM;

  public:

    /*****************/
    LinkElPh(const str_cit & b, const str_cit & e) {
      auto opts = get_key_val_args(b,e,
        {"type=", "C=", "moles=", "material=", "mass="});

      double C=0, moles=0;

      if (opts["material"] == "copper"){
        C = 2e3 * 7.11;
        if (opts["mass"]!="") moles = read_value(opts["mass"], "kg")/63.546e-3;
      }
      else {
        if (opts["mass"]!="") throw Err() << "mass parameter can be used only together with material";
      }

      if (opts["C"]      != "") C = read_value(opts["C"], "");
      if (C      <= 0) throw Err() << "A positive value expected: C";

      if (opts["moles"]  != "") moles = read_value(opts["moles"], "");
      if (moles  <= 0) throw Err() << "A positive value expected: moles";

      CM=C*moles;
    }

    /*****************/
    double get_qdot(const double T1, const double T2, const double B) const override {
      return (pow(T1,5)-pow(T2,5))*CM;
    }

};

/********************************************************************/
// Kapitza resistance between He3 and solid
// power=1:  R = 900/T/area [K^2/m^2/W]
// power=2:  R = 41/T^2/area -- Dilution fridge heat exchanger with fine silver powder, see Y.Takano-1994
// power=3:  R = 0.1/T^3/area -- Dilution fridge heat exchanger with coarse silver powder, see Y.Takano-1994

// Q = dT / R = area * T^n*dT / C ~= area (T1^(n+1)-T2^(n+1))/(n+1) /900

class LinkKapRes: public LinkBase {
  double area = 0;
  int power = 1;
  double C = 0;

  public:



    /*****************/
    LinkKapRes(const str_cit & b, const str_cit & e) {
      auto opts = get_key_val_args(b,e, {"type=", "area=", "power=", "C="});

      if (opts["area"]  != "") area = read_value(opts["area"], "m^2");
      if (opts["power"] != "") power = read_value(opts["power"], "");
      switch (power) {
        case 1: C=900.0; break;
        case 2: C=41.0; break;
        case 3: C=0.1; break;
        default: throw Err() << "unknown power setting for Kapitza resistance: " << power;
      }
      if (opts["C"] != "") C = read_value(opts["C"], "");
      if (C    <= 0) throw Err() << "A positive value expected: C";
      if (area <= 0) throw Err() << "A positive value expected: area";
    }

    /*****************/
    double get_qdot(const double T1, const double T2, const double B) const override {
      return area/((power+1)*C) * (pow(T1,power+1)-pow(T2,power+1));
    }
};

// Cooling power of a dilution refrigerator.
// block1 should be mixing chamber, block2 - thermal bath at any temperature.
class LinkDil: public LinkBase {
  double ndot = 0;
  public:
    /*****************/
    LinkDil(const str_cit & b, const str_cit & e) {
      auto opts = get_key_val_args(b,e, {"type=", "ndot="});
      if (opts["ndot"]  != "") ndot = read_value(opts["ndot"], "mol/s");
      if (ndot <= 0) throw Err() << "A positive value expected: ndot";
    }
    /*****************/
    double get_qdot(const double T1, const double T2, const double B) const override {
      return 84*T1*T1*ndot; }
};

// Heat transfer by circulation in a dilution refrigerator.
// Phase parameter is "C" or "D".
// T2 is not used in the calculation.
class LinkCirc: public LinkBase {
  double ndot = 0;
  int    phase = 0;
  public:
    /*****************/
    LinkCirc(const str_cit & b, const str_cit & e) {
      auto opts = get_key_val_args(b,e, {"type=", "ndot=", "phase="});
      if (opts["ndot"]  != "") ndot = read_value(opts["ndot"], "mol/s");
      if (opts["phase"] == "C")  phase=1;
      if (opts["phase"] == "D")  phase=2;
      if (ndot <= 0) throw Err() << "A positive value expected: ndot";
      if (phase<1 || phase>2) throw Err() << "Value is missing: phase (expected C or D)";
    }
    /*****************/
    double get_qdot(const double T1, const double T2, const double B) const override {
      if (phase == 1) return  23*T1*T1*ndot;
      if (phase == 2) return 107*T1*T1*ndot;
      throw Err() << "Wrong value of phase parameter: " << phase;
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
  if (type=="const")       return std::shared_ptr<LinkBase>(new LinkConst(b,e));
  if (type=="simple")      return std::shared_ptr<LinkBase>(new LinkSimple(b,e));
  if (type=="simple_bar")  return std::shared_ptr<LinkBase>(new LinkSimpleBar(b,e));
  if (type=="metal_bar")   return std::shared_ptr<LinkBase>(new LinkMetalBar(b,e));
  if (type=="korringa")    return std::shared_ptr<LinkBase>(new LinkKorringa(b,e));
  if (type=="el_ph")       return std::shared_ptr<LinkBase>(new LinkElPh(b,e));
  if (type=="kap_res_he3") return std::shared_ptr<LinkBase>(new LinkKapRes(b,e));
  if (type=="dilution_cooling") return std::shared_ptr<LinkBase>(new LinkDil(b,e));
  if (type=="dilution_circ")    return std::shared_ptr<LinkBase>(new LinkCirc(b,e));

  throw Err() << "unknown link type: " << type;
}

std::shared_ptr<LinkBase> create_link(const std::vector<std::string> & v) {
  return create_link(v.begin(), v.end());}


#endif
