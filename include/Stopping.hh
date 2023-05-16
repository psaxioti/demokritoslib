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

class Stopping {
private:
   /**
    * Stopping class constructor.
    *
    * It will read and hold the values needed for the calculation of stopping powers.
    */
   Stopping();
   static Stopping *StoppingInstance;

public:
   ~Stopping();

   /**
    * This is the static method that controls the access to the singleton
    * instance. On the first run, it creates a singleton object and places it
    * into the static field. On subsequent runs, it returns the client existing
    * object stored in the static field.
    */
   static Stopping *GetInstance();

   /// @brief Get stopping power
   /// @param En Beam energy for which the stopping power will be calculated
   /// @param BeamIso Isotope of beam
   /// @param target Complete target that the beam navigates
   /// @param Layer Layer of target the beam navigates
   ///   If not provided it is assumed it travels through the first layer of the target
   /// @return The complete stopping power of the target on the beam for the specified energy
   double GetStopping(double En, Isotope *BeamIso, Target *target, int Layer = 0);

   /// @brief Set the desired stopping calculation
   /// @param Stopping Can be either ZBL or SRIM.
   ///   Default is ZBL.
   ///   If parameters for the specific beam-target compination are missing the ZBL will be used.
   void SetStopping(const std::string &Stopping);

   /// @brief Calulate energy spread after traveling a target layer that absorbes energy LayerDE
   /// @param En Initial mean energy of beam in keV
   /// @param InitialEnSpread Initial energy spread of beab in keV
   /// @param LayerDE Energy lost in the layer traveled in keV
   /// @param Beam Beam isotope
   /// @param target Complete target that the beam navigates
   /// @param Layer Layer of target the beam navigates
   ///   If not provided it is assumed it travels through the first layer of the target
   /// @return The energy spread of the beam in keV after traveling in the target and losing LayerDE energy
   double GetStraggling(double En, double InitialEnSpread, double LayerDE, Isotope *Beam, Target *target, int Layer = 0);

private:
   double GetStraggling(double En, double InitialEnSpread, double LayerDE, Isotope *Beam, Element *TargetEl);
   /// @brief Get stopping power
   /// @param En Beam energy for which the stopping power will be calculated
   /// @param BeamIso Isotope of beam
   /// @param TargetEl Target element for which the stopping power is calculated
   /// @return The stopping power in keV per 10^18 atoms/cm2 of the beam isotope in a specific element
   double GetStopping(double En, Isotope *BeamIso, Element *TargetEl);
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
