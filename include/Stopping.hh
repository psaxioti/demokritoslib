#ifndef STOPPING_H
#define STOPPING_H 1

#include "Element.hh"
#include "Isotope.hh"
#include "Target.hh"

#include <map>
#include <string>
#include <vector>

#include <gsl/gsl_integration.h>
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
   ///   Default is SRIM.
   ///   If parameters for the specific beam-target compination are missing the ZBL will be used.
   void SetStopping(const std::string &Stopping);

   /// @brief Calulate energy spread after traveling a target layer that absorbes energy LayerDE
   /// @param InitialEnergy Initial mean energy in keV
   /// @param InitialEnergySpread Initial energy spread in keV
   /// @param LayerDE Energy lost in the layer traveled in keV
   /// @param Beam Isotope of beam
   /// @param target Complete target that the beam navigates
   /// @param Layer Layer of target the beam navigates
   ///   If not provided it is assumed it travels through the first layer of the target
   /// @return The energy spread of the beam in keV after traveling in the target and losing LayerDE energy
   double GetDEStraggling(double InitialEnergy, double InitialEnergySpread, double LayerDE, Isotope *Beam, Target *target, int Layer = 0);

   /// @brief Calculate range of beam in a target layer
   /// @param En Incoming energy in keV
   /// @param Beam Isotope of beam
   /// @param target Complete target that the beam navigates
   /// @param Layer Layer of target the beam navigates
   /// @return Range in 10^18 atoms/cm2
   double GetRange(double En, Isotope *Beam, Target *target, int Layer = 0);

   /// @brief Get thickness of target traversed
   /// @param Eout Energy of beam at the end of the traveled thickness
   /// @param Ein Incoming energy of beam
   /// @param Beam Isotope of beam
   /// @param target Complete target that the beam navigates
   /// @param Layer Layer of target the beam navigates
   /// @return Thickness of target layer traversed in 10^18 atoms/cm2
   double GetTraveledThickness(double Eout, double Ein, Isotope *Beam, Target *target, int Layer = 0);

   /// @brief Get outgoing energy after traversing a target layer
   /// @param En Incoming energy of beam in keV
   /// @param Beam Isotope of beam
   /// @param target Complete target that the beam navigates
   /// @param Layer Layer of target the beam navigates
   /// @return Outgoing energy from the target layer in keV
   double GetEout(double En, Isotope *Beam, Target *target, int Layer = 0);

   /// @brief Calulate energy spread after traveling a target layer
   /// @param Ein Initial mean energy in keV
   /// @param InitialEnergySpread Initial energy spread at Ein in keV
   /// @param Eout Energy in keV at which the straggling will be calculated
   /// @param Beam Isotope of beam
   /// @param target Complete target that the beam navigates
   /// @param Layer Layer of target the beam navigates
   /// @return The energy spread of the beam in keV after traveling in the target and losing energy to Eout from Ein
   double GetStraggling(double Ein, double InitialEnergySpread, double Eout, Isotope *Beam, Target *target, int Layer = 0);

   // To be removed!!!!!!!!!!!!!!!!!!!
   double CalculateStraggling(double Ein, double s0, double Eout, Isotope *Beam, Target *target, int Layer = 0);

private:
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

   double GetStatisticalStraggling2(double En, double LayerDE, Isotope *Beam, Element *TargetEl);
   double GetNonStatisticalStraggling2(double En, double InitialEnSpread, double LayerDE, Isotope *Beam, Target *target, int Layer = 0);
   double GetBohrElectronicStraggling2(double LayerThickness, Isotope *Beam, Element *TargetEl);
   double GetBohrNuclearStraggling2(double LayerThickness, Isotope *Beam, Element *TargetEl);
   double GetChuFactor(double En, Isotope *Beam, Element *TargetEl);

   void ReadZBL_Coeffs();
   void ReadSRIM_Coeffs();
   void ReadCHU_coeffs();

   std::string StoppingName;
   std::map<Element *, double[12]> ZBL_coeffs;
   std::map<int, std::map<int, std::vector<double>>> SRIM_coeffs;
   std::vector<double> SRIM_Energies;
   std::map<int, std::map<double, double>> CHU_coeffs;
   std::vector<double> CHU_EnergyPerMass;
   int CHU_Previous_ZTarget = 0;
   int SRIM_ZBeam_min;
   int SRIM_ZBeam_max;
   int SRIM_ZTarget_min;
   int SRIM_ZTarget_max;
   int SRIM_Previous_ZBeam, SRIM_Previous_ZTarget;
   gsl_interp_accel *srim_acc = nullptr;
   gsl_spline *srim_spline = nullptr;
   gsl_interp_accel *chu_acc = nullptr;
   gsl_spline *chu_spline = nullptr;

   gsl_integration_workspace *EnergyLossIntegrationWorkspace = nullptr;
   int IntegrationWorkspaceSize = 10000;
   double IntegrationepsabsDefault = 0.;
   double IntegrationepsrelDefault = 1e-3;
   double IntegrationEpsRelMultFactor = 3.;
   int QAGIntegrationKey = 1;

   double dxde(double en);
   friend double dxdeWrapper(double, void *);
   Isotope *ThicknessCalculationBeam;
   Target *ThicknessCalculationTarget;
   int ThicknessCalculationLayer;
   double ThicknessCalculationEin;

   double thicknessRoot(double en);
   friend double thicknessRootWrapper(double, void *);

   double IntegratedStraggling(double Ein, double s0, double Eout, Isotope *Beam, Element *el);
   friend double stragglingWrapper(double, void *);
   double straggling(double);
   Isotope *StragglingCalculationBeam;
   Element *StragglingCalculationElement;
};

#endif
