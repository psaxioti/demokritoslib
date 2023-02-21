#include "Reaction.hh"

#include <algorithm>
#include <fstream>
#include <map>

#include "Isotope.hh"
#include "Stopping.hh"
#include "Target.hh"

// const double AVOGADRO = 6.02214141070409084099072e23;
const double AVOGADRO = 6.02214141070409084099072e12;

Reaction::Reaction(int Z, int A, float gEnergy, float gTheta, std::string reaction, std::string Filename, std::string source)
    : AtomicNumber(Z), MassNumber(A), Egamma(gEnergy), Theta(gTheta), FileName(Filename),
      ReactionName(reaction), Source(source) {
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
   if (gsl_workspace)
      gsl_integration_workspace_free(gsl_workspace);
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
   return CrossSectionEnergy[0];
}

double Reaction::GetEmax() {
   return CrossSectionEnergy[GetNumberOfPoints() - 1];
}

int Reaction::GetNumberOfPoints() {
   return CrossSectionEnergy.size();
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
         line.erase(0, line.find_first_not_of(' '));
         for (int i = 0; i < 3; ++i) {
            auto NextDelim = line.find_first_of(" ,\t");
            double NextNumb = stod(line.substr(0, NextDelim + 1));
            if (NextNumb < 0. && i == 0)
               continue;
            if (i == 0)
               CrossSectionEnergy.push_back(stod(line.substr(0, NextDelim + 1)) * EnergyFactor);
            if (i == 2) {
               CrossSection.push_back(stod(line.substr(0, NextDelim + 1)) / CrossSectionFactor);
            }
            line.erase(0, NextDelim + 1);
         }
         int kk = CrossSectionEnergy.size();
         int kk1 = CrossSection.size();
      } else if (line.substr(0, 5) == "Zeds:") {
         line.erase(0, line.find(',') + 1);
         AtomicNumber = std::stoi(line.substr(0, line.find(',')));
      } else if (line.substr(0, 10) == "Enfactors:") {
         line.erase(0, line.find(':') + 1);
         if (line.empty())
            EnergyFactor = 1.;
         else
            EnergyFactor = stod(line.substr(0, line.find(',') + 1));
      } else if (line.substr(0, 7) == "Egamma:") {
         line.erase(0, 8);
         Egamma = stod(line);
      } else if (line.substr(0, 6) == "Theta:") {
         line.erase(0, 7);
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
         line.erase(0, 7);
         if (line.find("mb") != std::string::npos)
            CrossSectionFactor = 1.;
         else if (line.find("tot") != std::string::npos)
            CrossSectionFactor = 4. * M_PI;
         else
            CrossSectionFactor = 0.;
         //      } else if (line.substr(0, 5) == "Data:" && MassNumber > 0)
      } else if (line.substr(0, 5) == "Data:")
         ReadingData = true;
      //      else if (line.substr(0, 5) == "Data:" && MassNumber == 0)
      //         break;
   }
   file.close();

   MassNumber = stoi(ReactionName.substr(0, ReactionName.find_first_not_of("0123456789")));
   InitializeInterpolation();
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

   if (En < GetEmin() || En > GetEmax())
      return 0.;

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

double Reaction::GetYield(double En, double dE, double AtomicPerCent, double AtomicAbundunce, double Stopping) {
   if (En < GetEmin() || En > GetEmax())
      return 0.;

   double yield = GetCrossSection(En) * dE * 1e-27 * 1e15 * 6.24E12 * 1000 * (AtomicPerCent / 100.) * (AtomicAbundunce / 100);
   return yield / Stopping;
}

double ReactionYieldWrap(double z, void *user_data) {
   Reaction *this_ptr = (Reaction *)user_data;
   return this_ptr->YieldFunction(z);
}

double Reaction::YieldFunction(double En) {
   double f = YieldMultFactor * GetCrossSection(En) / YieldStopping->GetStopping(En, BeamIsotope, YieldTarget, YieldTargetLayer);
   return f;
}

void Reaction::SetIntegrationParameters(double WorkSpaceSize, double epsabs, double epsrel, double epsrelMult) {
   if (IntegrationWorkspaceSize != WorkSpaceSize || !gsl_workspace) {
      IntegrationWorkspaceSize = WorkSpaceSize;
      gsl_integration_workspace_free(gsl_workspace);
      gsl_workspace = gsl_integration_workspace_alloc(IntegrationWorkspaceSize);
   }
   IntegrationepsabsDefault = epsabs;
   IntegrationepsrelDefault = epsrel;
   IntegrationEpsRelMultFactor = epsrelMult;
   return;
}

void Reaction::SetIntegratedYield(double AtomicPerCent, double AtomicAbundunce, Stopping *Stop, Isotope *beamiso, Target *target, int layer) {
   YieldStopping = Stop;
   BeamIsotope = beamiso;
   YieldTarget = target;
   YieldTargetLayer = layer;
   YieldMultFactor = AVOGADRO * AtomicPerCent * AtomicAbundunce / 10000.;
   YieldMultFactor *= (1e-27 * 1e15 * 1000);

   SetIntegrationParameters(IntegrationWorkspaceSize, IntegrationepsabsDefault, IntegrationepsrelDefault, IntegrationEpsRelMultFactor);
   return;
}

double Reaction::GetIntegratedYield(double Emin, double Emax, double &error, double &epsrel) {
   if (Emax < GetEmin() ||
       Emax < Emin)
      return 0.;

   Emin = std::max(Emin, GetEmin());
   Emax = std::min(Emax, GetEmax());

   gsl_function F;
   F.function = &ReactionYieldWrap;
   F.params = (void *)this;

   double epsabs = IntegrationepsabsDefault;
   epsrel = IntegrationepsrelDefault / IntegrationEpsRelMultFactor;
   double limit = IntegrationWorkspaceSize;
   int QAGIntegrationKey = 5;

   int status = 1;
   gsl_set_error_handler_off();
   double result;
   while (status) {
      epsrel *= IntegrationEpsRelMultFactor;
      //      status = gsl_integration_qags(&F, CrossSectionEnergy[0], En,
      //                                    epsabs, epsrel, limit, gsl_workspace, &result, &error);
      status = gsl_integration_qag(&F, Emin, Emax,
                                   epsabs, epsrel, limit, QAGIntegrationKey, gsl_workspace, &result, &error);
   }
   return result;
}
