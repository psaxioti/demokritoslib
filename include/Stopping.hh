#ifndef STOPPING_H
#define STOPPING_H 1

#include "Element.hh"
#include "Isotope.hh"
#include "Target.hh"

#include <map>
#include <string>
#include <vector>

#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>

// Everything in keV and 10^18 atoms/cm2

class Stopping {
public:
   Stopping();
   ~Stopping();

   double GetStopping(double En, Isotope *BeamIso, Element *TargetEl);
   double GetStopping(double En, Isotope *BeamIso, Target *target, int Layer = 0);
   void SetStopping(std::string &Stopping);

   double GetStraggling(double En, double InitialEnSpread, double LayerDE, Isotope *Beam, Element *TargetEl);
   double GetStraggling(double En, double InitialEnSpread, double LayerDE, Isotope *Beam, Target *target, int Layer = 0);

private:
   double ZBL_Stopping(double En, Isotope *BeamIso, Element *TargetEl);
   double ZBL_ElectronicStopping(double En, Isotope *BeamIso, Element *TargetEl);
   double ZBL_NuclearStopping(double En, Isotope *BeamIso, Element *TargetEl);
   double SRIM_Stopping(double En, Isotope *BeamIso, Element *TargetEl);

   double GetStraggling2(double En, double InitialEnSpread, double LayerDE, Isotope *Beam, Element *TargetEl);
   double GetNonStatisticalStraggling2(double IncomingStopping, double OutgoingEn, double InitialEnSpread, Isotope *Beam, Element *TargetEl);
   double GetBohrElectronicStraggling2(double LayerThickness, Isotope *Beam, Element *TargetEl);
   double GetBohrNuclearStraggling2(double LayerThickness, Isotope *Beam, Element *TargetEl);
   double GetChuFactor(double En, Isotope *Beam, Element *TargetEl);

   void ReadZBL_Coeffs();
   void ReadSRIM_Coeffs();

   std::string StoppingName;
   std::map<Element *, double[12]> ZBL_coeffs;
   std::map<int, std::map<int, std::vector<double>>> SRIM_coeffs;
   std::vector<double> SRIM_Energies;
   int SRIM_ZBeam_min;
   int SRIM_ZBeam_max;
   int SRIM_ZTarget_min;
   int SRIM_ZTarget_max;
   int SRIM_Previous_ZBeam, SRIM_Previous_ZTarget;
   gsl_interp_accel *acc = nullptr;
   gsl_spline *spline = nullptr;
};

#endif