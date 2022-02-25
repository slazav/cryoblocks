#ifndef C_BLOCKS_H
#define C_BLOCKS_H

// Blocks

/********************************************************************/
// Base class for a thermal block.
class BlockBase {
  // Calculate temperature change.
  // dQ = T dS = T(partial dS/dT) dT + T(partial dS/dB) dB = C dT + D dB
  public: virtual double get_dt(
    const double dQ, const double T, const double B, const double dB) const = 0;
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

#endif
