#ifndef REACTION_H
#define REACTION_H 1

#include <filesystem>
#include <string>
#include <vector>

#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>

class Reaction {
public:
   Reaction(int, int, float, float, int, float, float, std::string, std::string, std::string);
   explicit Reaction(std::filesystem::path FileName);
   ~Reaction();

   int GetAtomicNumber();
   int GetMassNumber();
   std::string GetReactionName();
   double GetTheta();
   std::string GetFileName();
   double GetEgamma();
   double GetEmin();
   double GetEmax();
   int GetNumberOfPoints();
   std::string GetSource();

   double GetCrossSection(double En);
   double GetStragglingCrossSection(double En, double sigma);

private:
   void ReadR33File();

   int AtomicNumber;
   int MassNumber = 0;
   std::string ReactionName;
   double Theta;
   std::string FileName;
   double Egamma;
   double Emin, Emax;
   int NumberOfPoints;
   std::string Source;

   double EnergyFactor;
   double CrossSectionFactor;

   std::vector<double> CrossSectionEnergy;
   std::vector<double> CrossSection;

   void InitializeInterpolation();

   gsl_interp_accel *acc = nullptr;
   gsl_spline *spline = nullptr;
};

#endif