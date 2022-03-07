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
    double Q;
  public:
    LinkConst(const double Q): Q(Q){;}


    static std::shared_ptr<LinkBase> create (const str_cit & b, const str_cit & e) {
      auto opts = get_key_val_args(b,e, {"type=", "Qdot=0W"});
      auto Q = read_power(opts["Qdot"]);
      return std::shared_ptr<LinkBase>(new LinkConst(Q));
    }

    double get_qdot(const double T1, const double T2, const double B) const override { return Q;}

};

/********************************************************************/
// Bar made of material with heat conductivity K = Ka*T^Kb,
// and Area-to-length ratio SL [m].
class LinkSimpleBar: public LinkBase {
    double Ka, Kb, SL;
  public:
    LinkSimpleBar(const double Ka, const double Kb, const double SL): Ka(Ka), Kb(Kb), SL(SL){;}

    double get_qdot(const double T1, const double T2, const double B) const override {
      // Q = SL * \int_T1^T2 K(T) dT
      return SL*Ka/(Kb+1) * (pow(T1, Kb+1) - pow(T2, Kb+1));
    }

    static std::shared_ptr<LinkBase> create (const str_cit & b, const str_cit & e) {
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
};

/********************************************************************/
// Bar made of metal. Total resistance R.
class LinkMetalBar: public LinkBase {
    double R;
    double l = 2.44e-8;
  public:
    LinkMetalBar(const double R): R(R){;}

    double get_qdot(const double T1, const double T2, const double B) const override {
      // K = T*lambda/R
      return (T1*T1-T2*T2)/2 * l/R;
    }

    static std::shared_ptr<LinkBase> create (const str_cit & b, const str_cit & e) {
      auto opts = get_key_val_args(b,e, {"type=", "R="});
      double R = read_resistance(opts["R"]);
      return std::shared_ptr<LinkBase>(new LinkMetalBar(R));
    }
};


/********************************************************************/
// Spin-lattice relaxation, Korringa law.
// block1 should be spin system, block2 - electrons
class LinkKorringa: public LinkBase {
  double Bint; // internal field [T]
  double gyro; // gyromagnetic ratio [rad/s/T]
  double spin; // 1/2, 3/2, etc.
  double kappa0; // high-field value of Korringa constant
  double alpha;  // alpha parameter in field dependence of kappa
  double nmole;  // number of moles

  double hbar = 1.05457181710e-34; // [J s]
  double NA = 6.02214086e+23;      // [1/mol] Avogadro's Constant
  double kB = 1.38064852e-23;      // [m^2 kg s-2 K-1] [J/K] Boltzmann constant
  double R = NA*kB;                // R-constant, [J/mol/K]
  double mu0 = 4*M_PI*1e-7;        // [V s/A/m] Vacuum permeability
  double muN = 5.05e-27;           // [A m^2] nuclear magneton

  public:
    LinkKorringa(const double Bint, const double gyro,
                 const double spin, const double kappa0,
                 const double alpha, const double nmole):
                   Bint(Bint), gyro(gyro), spin(spin), kappa0(kappa0),
                   alpha(alpha), nmole(nmole){;}

    double get_qdot(const double T1, const double T2, const double B) const override {
      // heat capacity of copper nuclei:
      double x = gyro*hbar*sqrt(B*B + Bint*Bint)/kB/T1/2.0;
      double y = (2.0*spin + 1.0)*x;
      double C = R*nmole *(pow(x/sinh(x),2) - pow(y/sinh(y),2));

      // See Pobell book f.10.10 and below
      double kappa = kappa0 * (B*B + Bint*Bint)/(B*B + alpha*Bint*Bint);
      double z = (mu0*muN*B)/(2*kB);
      double tau = kappa/z * tanh(z/T2);
      return C * T1/T2*(T1-T2) / tau;
    }

    static std::shared_ptr<LinkBase> create (const str_cit & b, const str_cit & e) {
      auto opts = get_key_val_args(b,e,
        {"type=", "Bint=", "gyro=", "spin=", "kappa0=", "alpha=", "moles=", "material="});

      double Bint=0, gyro=0, spin=0, kappa0=0, alpha=0, nmol=0;

      if (opts["material"] == "copper"){
        Bint = 0.36e-3;    // [T], dipolar feld in copper
        gyro = 71.118e6;   // [rad/s/T] gyromagnetic ratio of copper
        spin = 1.5;        // spin 3/2
        kappa0 = 1.2;
        alpha  = 2.6;
      }

      if (opts["Bint"]   != "") Bint = read_magn_field(opts["Bint"]);
      if (opts["gyro"]   != "") gyro = read_gyro(opts["gyro"]);
      if (opts["spin"]   != "") spin = read_dimensionless(opts["spin"]);
      if (opts["kappa0"] != "") kappa0 = read_kappa(opts["kappa0"]);
      if (opts["alpha"]  != "") alpha  = read_dimensionless(opts["alpha"]);
      if (opts["moles"]  != "") nmol = read_dimensionless(opts["moles"]);
      if (gyro   <= 0) throw Err() << "A positive value expected: gyro";
      if (spin   <= 0) throw Err() << "A positive value expected: spin";
      if (kappa0 <= 0) throw Err() << "A positive value expected: kappa0";
      if (nmol   <= 0) throw Err() << "A positive value expected: moles";
      return std::shared_ptr<LinkBase>(new LinkKorringa(Bint, gyro, spin, kappa0, alpha, nmol));
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
  if (type=="const")      return LinkConst::create(b,e);
  if (type=="simple_bar") return LinkSimpleBar::create(b,e);
  if (type=="metal_bar")  return LinkMetalBar::create(b,e);
  if (type=="korringa")   return LinkKorringa::create(b,e);

  throw Err() << "unknown link type: " << type;
}


#endif
