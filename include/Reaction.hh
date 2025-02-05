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
class Element;
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
   //   Reaction(int Z, int A, float gEnergy, float gTheta, std::string reaction, std::string Filename, std::string source);

   /// @brief Reaction class constructor.
   /// @param FileName Name of r33 file.
   explicit Reaction(std::filesystem::path FileName);
   ~Reaction();

   /// Get Atomic number of reaction's target
   int GetTargetAtomicNumber();
   /// Get Mass number of reaction's target
   int GetTargetMassNumber();
   /// Get reaction's target isotope
   Isotope *GetTargetIsotope();
   /// Get reaction's target element
   Element *GetTargetElement();
   /// Get Atomic number of reaction's beam
   int GetBeamAtomicNumber();
   /// Get Mass number of reaction's beam
   int GetBeamMassNumber();
   /// Get reaction's beam isotope
   Isotope *GetBeamIsotope();
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
   /// Get minimum energy step of cross section
   double GetEnergyMinStep();
   /// Get number of cross sections' points
   int GetNumberOfPoints();
   /// Get reference of the cross sections
   std::string GetSource();

   /// Get cross section for energy En (keV)
   // Linear interpolation of the cross section points, 0. outside the energy range
   double GetCrossSection(double En);

   /// @brief Get cross section convoluted with a gaussian
   /// @param En Central energy of the cross section
   /// @param sigma Sigma of the gaussian
   /// @return Convolution of the cross section at the specified energy with the given sigma
   double GetStragglingCrossSection(double En, double sigma);

   /// @brief Get cross section convoluted with a gaussian
   /// @param En Central energy of the cross section
   /// @param InitialBeamEnergy Initial Beam energy. Has to be greater than En
   /// @param InitialBeamEnegySpread Initial energy spread of the beam energy
   /// @param target Target in which the beam travels
   /// @return Convolution of the cross section at the specified energy with sigma calculated from the straggling of the beam in the target
   double GetStragglingCrossSection(double En, double InitialBeamEnergy, double InitialBeamEnegySpread, Target *target);

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

   /// @brief Get yield.
   /// Calculated as the integral of Factors*CrossSection/Stopping between Emin and Emax.
   /// @param Emin Lower limit of the integration. If Emin is lower than the stopping power in the target, then the yield in the Layer will be calculated.
   /// @param Emax Upper limit of the integration. Considered as beam energy.
   /// @param error Returns the absolute estimated error of the integration.
   /// @param epsrel Returns the epsrel used for the integration.
   /// @param target Target in which the beam travels
   /// @return Yield
   double GetYield(double Emin, double Emax, double &error, double &epsrel, Target *target);

   /// @brief Get yield.
   /// Calculated as the integral of Factors*CrossSection/Stopping between Emin and InitialBeamEnergy.
   /// Straggling is also considered
   /// @param Emin Lower limit of the integration. If Emin is lower than the stopping power in the target, then the yield in the Layer will be calculated.
   /// @param InitialBeamEnergy Initial Beam energy. Has to be greater than Emin
   /// @param InitialBeamEnegySpread Initial energy spread of the beam energy.
   ///     If negative value is provided, then no straggling will be taken into account!!!
   /// @param error Returns the absolute estimated error of the integration.
   /// @param epsrel Returns the epsrel used for the integration.
   /// @param target Target in which the beam travels
   /// @return Yield
   double GetYield(double Emin, double InitialBeamEnergy, double InitialBeamEnegySpread, double &error, double &epsrel, Target *target);

   /// @brief Set gsl Integration parameters. The algorithm used is QAG with default key 5
   /// @param WorkSpaceSize Sets the size of the workspace. To the same value the limit is also set. Default value is 10000
   /// @param epsabs Sets epsabs value. Default is 0.
   /// @param epsrel Sets epsrel initial value. Default is 1e-6.
   /// @param epsrelMult Sets multiplication factor to change epsrel if no convergeance is reached.
   /// @param QAG_Key Sets the QAG algorithm key. Default is 5.
   void SetIntegrationParameters(double WorkSpaceSize = 10000, double epsabs = 0, double epsrel = 1.e-6, double epsrelMult = 3, int QAG_Key = 5);

   // To be removed!!!!!!!!!
   double GetYield(double En, double dE, double AtomicPerCent, double AtomicAbundunce, double Stopping, double enSigma);

private:
   void ReadR33File(bool ReadCSs = false);

   Isotope *ReactionBeamIsotope;
   Isotope *ReactionTargetIsotope;
   Element *ReactionTargetElement;

   std::string FileName;
   std::string ReactionName;
   std::string Source;
   double Theta;
   double Egamma;

   double ReactionMinEnergyStep = 10000.;
   std::vector<double> CrossSectionEnergy;
   std::vector<double> CrossSection;

   // Functions/Parameters relative to cross section interpolation
   void InitializeInterpolation();
   gsl_interp_accel *CS_InterpolationAccellerator = nullptr;
   gsl_spline *CS_InterpolationSpline = nullptr;

   // General integration Parameters
   int IntegrationWorkspaceSize = 50000;
   double IntegrationepsabsDefault = 0.;
   double IntegrationepsrelDefault = 1e-6;
   double IntegrationEpsRelMultFactor = 3.;
   int QAGIntegrationKey = 2;

   void SetYieldCalculationParameters(Target *target, int layer, double ein, double sigma);
   double YieldCalculation(double);
   friend double YieldCalculationWrap(double, void *);
   gsl_integration_workspace *YieldCalculationWorkspace = nullptr;
   Target *YieldCalculationTarget;
   int YieldCalculationLayer;
   double YieldCalculationMultFactor;
   double YieldCalculationInitialSigma;
   double YieldCalculationEin;
   double GetStragglingCrossSectionOverStopping(double En, double sigma);
   double StragglingCrossSectionOverStopping(double E);
   friend double StragglingCrossSectionOverStoppingWrap(double z, void *user_data);
   double YieldCalculationCurrentE;
   double YieldCalculationCurrentSigma;
   gsl_integration_workspace *StragglingCrossSectionOverStoppingWorkspace = nullptr;

   void SetIntegratedYield(Target *target, int layer = 0, double sigma = 0);
   double YieldFunction(double x);
   friend double ReactionYieldWrap(double, void *);
   gsl_integration_workspace *YieldIntegrationWorkspace = nullptr;
   Target *YieldTarget = nullptr;
   int YieldTargetLayer = 0;
   double YieldMultFactor;
   double YieldEnergySigma = 0;

   double StragglingYieldFunction(double x);
   friend double StragglingReactionYieldWrap(double, void *);
   void SetStragglingParameters(double En, double sigma);
   void SetSIntegratedYield(Target *target, int layer, double ein, double sigma);
   double StragglingCrossSection(double E);
   friend double StragglingCrossSectionWrap(double z, void *user_data);
   gsl_integration_workspace *StragglingIntegrationWorkspace = nullptr;
   double StragglingEin;
   double StragglingEn;
   double StragglingSigma;
};

#endif