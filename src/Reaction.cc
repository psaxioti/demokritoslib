#include "Reaction.hh"

#include <algorithm>
#include <fstream>
#include <map>

Reaction::Reaction(int Z, int A, float Energy, float th, int No, float emin, float emax, std::string react, std::string file, std::string sour)
    : AtomicNumber(Z), MassNumber(A), Egamma(Energy), Theta(th), FileName(file),
      Emin(emin), Emax(emax), NumberOfPoints(No), ReactionName(react), Source(sour) {
   ReadR33File();
}

Reaction::Reaction(std::filesystem::path file) {
   FileName = std::string(file.parent_path()) + "/" + std::string(file.filename());
   ReadR33File();
}

Reaction::~Reaction() {
   if (spline)
      gsl_spline_free(spline);
   if (acc)
      gsl_interp_accel_free(acc);
}

int Reaction::GetAtomicNumber() {
   return AtomicNumber;
}

int Reaction::GetMassNumber() {
   return MassNumber;
}

std::string Reaction::GetReactionName() {
   return ReactionName;
}

double Reaction::GetTheta() {
   return Theta;
}

std::string Reaction::GetFileName() {
   return FileName;
}

double Reaction::GetEgamma() {
   return Egamma;
}

double Reaction::GetEmin() {
   return Emin;
}

double Reaction::GetEmax() {
   return Emax;
}

int Reaction::GetNumberOfPoints() {
   return NumberOfPoints;
}

std::string Reaction::GetSource() {
   return Source;
}

void Reaction::ReadR33File() {
   if (!ReactionName.empty() && CrossSection.size() > 0)
      return;
   std::ifstream file;
   file.open(FileName);

   bool ReadingData = false;
   while (!file.eof()) {
      std::string line;
      std::getline(file, line, '\n');

      if (line.substr(0, 8) == "EndData:") {
         if (CrossSectionEnergy[0] > CrossSectionEnergy[1]) {
            std::reverse(CrossSectionEnergy.begin(), CrossSectionEnergy.end());
            std::reverse(CrossSection.begin(), CrossSection.end());
         }
         break;
      } else if (ReadingData) {
         for (int i = 0; i < 3; ++i) {
            auto NextDelim = line.find_first_of(" ,\t");
            if (i == 0)
               CrossSectionEnergy.push_back(stod(line.substr(0, NextDelim)) * EnergyFactor);
            if (i == 2) {
               CrossSection.push_back(stod(line.substr(0, NextDelim)) / CrossSectionFactor);
            }
            line.erase(0, NextDelim);
         }
      } else if (line.substr(0, 5) == "Zeds:") {
         line.erase(0, line.find(',') + 1);
         AtomicNumber = std::stoi(line.substr(0, line.find(',')));
      } else if (line.substr(0, 10) == "Enfactors:") {
         line.erase(0, line.find(',') + 1);
         if (line.empty())
            EnergyFactor = 1.;
         else
            EnergyFactor = stod(line.substr(0, line.find(',')));
      } else if (line.substr(0, 7) == "Egamma:") {
         line.erase(0, 7);
         Egamma = stod(line);
      } else if (line.substr(0, 6) == "Theta:") {
         line.erase(0, 6);
         Theta = stod(line);
      } else if (line.substr(0, 9) == "Reaction:") {
         ReactionName = line.substr(10);
         auto first = ReactionName.find_first_not_of(" \t");
         auto last = ReactionName.find_last_not_of(" \t") + 1;
         ReactionName = ReactionName.substr(first, last - first);
      } else if (line.substr(0, 7) == "Source:") {
         Source = line.substr(8);
         auto first = Source.find_first_not_of(" \t");
         auto last = Source.find_last_not_of(" \t") + 1;
         Source = Source.substr(first, last - first);
      } else if (line.substr(0, 6) == "Units:") {
         line.erase(0, 6);
         if (line.find("mb"))
            CrossSectionFactor = 1.;
         else if (line.find("tot"))
            CrossSectionFactor = 4. * M_PI;
         else
            CrossSectionFactor = 0.;
      } else if (line.substr(0, 5) == "Data:" && MassNumber > 0)
         ReadingData = true;
      else if (line.substr(0, 5) == "Data:" && MassNumber == 0)
         break;
   }
   file.close();

   MassNumber = stoi(ReactionName.substr(0, ReactionName.find_first_not_of("0123456789")));
}

void Reaction::InitializeInterpolation() {
   double *x = &CrossSectionEnergy[0];
   double *y = &CrossSection[0];

   if (spline)
      gsl_spline_free(spline);
   if (acc)
      gsl_interp_accel_free(acc);
   acc = gsl_interp_accel_alloc();
   spline = gsl_spline_alloc(gsl_interp_linear, CrossSection.size());
   //   spline = gsl_spline_alloc(gsl_interp_cspline, CrossSection.size());
   //   spline = gsl_spline_alloc(gsl_interp_akima, CrossSection.size());
   //   spline = gsl_spline_alloc(gsl_interp_steffen, CrossSection.size());

   gsl_spline_init(spline, x, y, CrossSection.size());
}

double Reaction::GetCrossSection(double En) {
   if (CrossSection.size() == 0)
      ReadR33File();

   if (!acc)
      InitializeInterpolation();

   return gsl_spline_eval(spline, En, acc);
}

double Reaction::GetStragglingCrossSection(double En, double sigma) {
   double Esum_gauss, cross_sum;
   int dE_strag;

   cross_sum = 0;
   Esum_gauss = En - 3 * sigma;
   dE_strag = 1; ////// Straggling integration step
   while (Esum_gauss < En + 3 * sigma) {
      cross_sum += GetCrossSection(Esum_gauss) * (1 / (sigma * sqrt(2 * M_PI))) * exp(-pow((Esum_gauss - En), 2) / (2 * pow(sigma, 2)));
      Esum_gauss += dE_strag;
   }
   return cross_sum;
}
