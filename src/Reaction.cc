#include "Reaction.hh"

#include <algorithm>
#include <fstream>
#include <map>

#include "MassFunctions.hh"
#include "Stopping.hh"
#include "Target.hh"

// const double AVOGADRO = 6.02214141070409084099072e23;
const double AVOGADRO = 6.02214141070409084099072e12;

// Reaction::Reaction(int Z, int A, float gEnergy, float gTheta, std::string reaction, std::string Filename, std::string source)
//     : AtomicNumber(Z), MassNumber(A), Egamma(gEnergy), Theta(gTheta), FileName(Filename),
//       ReactionName(reaction), Source(source) {
//    ReadR33File();
// }

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
   return ReactionTargetIsotope->GetZ();
}

int Reaction::GetMassNumber() {
   return ReactionTargetIsotope->GetA();
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
   if (CrossSection.size() == 0)
      ReadR33File(true);
   return CrossSectionEnergy[0];
}

double Reaction::GetEmax() {
   if (CrossSection.size() == 0)
      ReadR33File(true);
   return CrossSectionEnergy[GetNumberOfPoints() - 1];
}

double Reaction::GetEnergyMinStep() {
   if (CrossSection.size() == 0)
      ReadR33File(true);
   return ReactionMinEnergyStep;
}

int Reaction::GetNumberOfPoints() {
   if (CrossSection.size() == 0)
      ReadR33File(true);
   return CrossSectionEnergy.size();
}

std::string Reaction::GetSource() {
   return Source;
}

void Reaction::ReadR33File(bool ReadCSs) {
   if (!ReactionName.empty() && CrossSection.size() > 0)
      return;
   std::ifstream file;
   file.open(FileName);

   bool ReadingData = false;
   static double EnergyFactor;
   static double CrossSectionFactor;

   int BeamZ, TargetZ;
   int BeamA, TargetA;

   while (!file.eof()) {
      std::string line;
      std::getline(file, line, '\n');

      auto div_char = line.find(':');
      std::string read_val = line.substr(0, div_char);

      if (read_val == "EndData") {
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
               CrossSectionEnergy.push_back(NextNumb * EnergyFactor);
            if (i == 2) {
               CrossSection.push_back(NextNumb / CrossSectionFactor);
               if (CrossSectionEnergy.size() > 1) {
                  double newmin = std::abs(CrossSectionEnergy[CrossSectionEnergy.size() - 1] - CrossSectionEnergy[CrossSectionEnergy.size() - 2]);
                  ReactionMinEnergyStep = std::min(ReactionMinEnergyStep, newmin);
               }
            }
            line.erase(0, NextDelim + 1);
         }
      } else if (read_val == "Zeds") {
         line.erase(0, div_char + 1);
         for (int i = 0; i < 2; ++i) {
            auto NextDelim = line.find(',');
            double NextNumb = stoi(line.substr(0, NextDelim + 1));
            if (i == 0)
               BeamZ = NextNumb;
            else if (i == 1)
               TargetZ = NextNumb;
            line.erase(0, NextDelim + 1);
         }
      } else if (read_val == "Enfactors") {
         line.erase(0, div_char + 1);
         line.erase(0, line.find_first_not_of(' '));
         if (line.empty())
            EnergyFactor = 1.;
         else
            EnergyFactor = stod(line.substr(0, line.find(',') + 1));
      } else if (read_val == "Egamma") {
         line.erase(0, div_char + 1);
         Egamma = stod(line);
      } else if (read_val == "Theta") {
         line.erase(0, div_char + 1);
         Theta = stod(line);
      } else if (read_val == "Reaction") {
         ReactionName = line.erase(0, div_char + 1);
         auto first = ReactionName.find_first_not_of(" \t");
         auto last = ReactionName.find_last_not_of(" \t") + 1;
         ReactionName = ReactionName.substr(first, last - first);
      } else if (read_val == "Source") {
         Source = line.erase(0, div_char + 1);
         auto first = Source.find_first_not_of(" \t");
         auto last = Source.find_last_not_of(" \t") + 1;
         Source = Source.substr(first, last - first);
      } else if (read_val == "Units") {
         line.erase(0, div_char + 1);
         if (line.find("mb") != std::string::npos)
            CrossSectionFactor = 1.;
         else if (line.find("tot") != std::string::npos)
            CrossSectionFactor = 4. * M_PI;
         else
            CrossSectionFactor = 0.;
      } else if (read_val == "Data" && ReadCSs)
         ReadingData = true;
      else if (read_val == "Data")
         break;
   }
   file.close();

   TargetA = stoi(ReactionName.substr(0, ReactionName.find_first_not_of("0123456789")));
   auto firstParam = ReactionName.find('(');
   auto firstComma = ReactionName.find(',');
   std::string beam = ReactionName.substr(firstParam + 1, firstComma - firstParam - 1);
   if (beam == "p")
      BeamA = 1;
   else if (beam == "d")
      BeamA = 2;

   ReactionBeamIsotope = GetIsotope(BeamZ, BeamA);
   ReactionTargetElement = GetElement(TargetZ);
   ReactionTargetIsotope = GetIsotope(TargetZ, TargetA);
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
      ReadR33File(true);

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
   double f = YieldMultFactor * GetCrossSection(En) / Stopping::GetInstance()->GetStopping(En, ReactionBeamIsotope, YieldTarget, YieldTargetLayer);
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

void Reaction::SetIntegratedYield(Target *target, int layer) {
   YieldTarget = target;
   YieldTargetLayer = layer;
   YieldMultFactor = AVOGADRO * YieldTarget->GetElementAtomicPercentInLayer(ReactionTargetElement, layer) * ReactionTargetIsotope->GetAbundance() / 10000.;
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
