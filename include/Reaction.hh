#ifndef REACTION_H
#define REACTION_H 1

#include <filesystem>
#include <string>
#include <vector>

#include "gsl/gsl_errno.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>

class Isotope;
class Target;
class Stopping;

class Reaction {
public:
   /// @brief Reaction class constructor.
   /// @param Z Atomic number of target.
   /// @param A Mass number of target.
   /// @param gEnergy Energy of gamma detected.
   /// @param gTheta Theta angle of emitted gamma.
   /// @param reaction Reaction name.
   /// @param Filename Name of r33 file.
   /// @param source Reference.
   Reaction(int Z, int A, float gEnergy, float gTheta, std::string reaction, std::string Filename, std::string source);
   /// @brief Reaction class constructor.
   /// @param FileName Name of r33 file.
   explicit Reaction(std::filesystem::path FileName);
   ~Reaction();

   /// Get Atomic number of reaction's target
   int GetAtomicNumber();
   /// Get Mass number of reaction's target
   int GetMassNumber();
   /// Get name of reaction
   std::string GetReactionName();
   /// Get theta angle of emmited gamma
   double GetTheta();
   /// Get name of r33 file
   std::string GetFileName();
   /// Get energy of emmited gamma
   double GetEgamma();
   /// Get minimum beam energy of the cross sections
   double GetEmin();
   /// Get maximum beam energy of the cross sections
   double GetEmax();
   /// Get number of cross sections' points
   int GetNumberOfPoints();
   /// Get reference of the cross sections
   std::string GetSource();

   /// Get cross section for energy En (keV)
   // Linear interpolation of the cross section points, 0. outside the energy range
   double GetCrossSection(double En);
   /// Get cross section for energy En (keV) convoluted with a gaussian of sigma
   // The convolution is done by arithmetic sum at the moment!!!
   double GetStragglingCrossSection(double En, double sigma);
   /**
    * Get yield.
    * Calculated as the simple calculation of Factors*CrossSection*dE/Stopping
    *
    * @param En Energy at which the cross section will be calculated.
    * @param dE Energy step of the calculation.
    * @param AtomicPerCent Atomic percentage of the element in the target (0 - 100).
    * @param AtomicAbundunce Isotopic abundunce of the isotope in the target (0 - 100).
    * @param Stopping Stopping power of the beam in the target in keV per 10^18 atoms/cm2.
    *    It is advised to give it for the same En energy
    */
   double GetYield(double En, double dE, double AtomicPerCent, double AtomicAbundunce, double Stopping);
   /**
    * Set variables needed for yield calculation through integration.
    *
    * @param AtomicPerCent Atomic percentage of the element in the target (0 - 100).
    * @param AtomicAbundunce Isotopic abundunce of the isotope in the target (0 - 100).
    * @param Stop Instance of Stopping class for calculating the stopping power.
    * @param beamiso Instance of Isotope class with beam isotope.
    * @param target Instance of Target class with the specific target.
    * @param layer Layer of the target that the particle navigates.
    *    [optional] If not provided it will be calculated for first layer of the target.
    */
   void SetIntegratedYield(double AtomicPerCent, double AtomicAbundunce, Stopping *Stop, Isotope *beamiso, Target *target, int layer = 0);
   /**
    * Get yield.
    * Calculated as the integral of Factors*CrossSection/Stopping between Emin and Emax
    *
    * @param Emin Lower limit of the integration.
    * @param Emax Upper limit of the integration.
    * @param error Returns the absolute estimated error of the integration.
    * @param epsrel Returns the epsrel used for the integration.
    */
   double GetIntegratedYield(double Emin, double Emax, double &error, double &epsrel);
   /**
    * Set gsl Integration parameters. The algorithm used is QAG with key 5
    *
    * @param WorkSpaceSize Sets the size of the workspace. To the same value the limit is also set.
    *    Default value is 10000
    * @param epsabs Sets epsabs value. Default is 0.
    * @param epsrel Sets epsrel initial value. Default is 1e-3.
    * @param epsrelMult Sets multiplication factor to change epsrel if no convergeance is reached.
    */
   void SetIntegrationParameters(double WorkSpaceSize, double epsabs, double epsrel, double epsrelMult);

private:
   void ReadR33File();

   int AtomicNumber;
   int MassNumber = 0;
   std::string ReactionName;
   double Theta;
   std::string FileName;
   double Egamma;
   std::string Source;

   double EnergyFactor;
   double CrossSectionFactor;

   std::vector<double> CrossSectionEnergy;
   std::vector<double> CrossSection;

   void InitializeInterpolation();

   gsl_interp_accel *acc = nullptr;
   gsl_spline *spline = nullptr;

   double YieldFunction(double x);
   friend double ReactionYieldWrap(double, void *);
   gsl_integration_workspace *gsl_workspace = nullptr;
   Stopping *YieldStopping = nullptr;
   Isotope *BeamIsotope = nullptr;
   Target *YieldTarget = nullptr;
   int YieldTargetLayer = 0;
   double YieldMultFactor;
   int IntegrationWorkspaceSize = 10000;
   double IntegrationepsabsDefault = 0.;
   double IntegrationepsrelDefault = 1e-3;
   double IntegrationEpsRelMultFactor = 3.;
};

#endif