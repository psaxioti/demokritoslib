#include "Reaction.hh"

#include <algorithm>
#include <fstream>
#include <map>

#include "MassFunctions.hh"
#include "Stopping.hh"
#include "Target.hh"

#include "gsl/gsl_const_mksa.h"

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
   if (CS_InterpolationSpline)
      gsl_spline_free(CS_InterpolationSpline);
   if (CS_InterpolationAccellerator)
      gsl_interp_accel_free(CS_InterpolationAccellerator);
   if (YieldIntegrationWorkspace)
      gsl_integration_workspace_free(YieldIntegrationWorkspace);
   if (StragglingIntegrationWorkspace)
      gsl_integration_workspace_free(StragglingIntegrationWorkspace);
}

int Reaction::GetTargetAtomicNumber() {
   return ReactionTargetIsotope->GetZ();
}

int Reaction::GetTargetMassNumber() {
   return ReactionTargetIsotope->GetA();
}

Isotope *Reaction::GetTargetIsotope() {
   return ReactionTargetIsotope;
}

Element *Reaction::GetTargetElement() {
   return GetElement(ReactionTargetIsotope->GetZ());
}

int Reaction::GetBeamAtomicNumber() {
   return ReactionBeamIsotope->GetZ();
}

int Reaction::GetBeamMassNumber() {
   return ReactionBeamIsotope->GetA();
}

Isotope *Reaction::GetBeamIsotope() {
   return ReactionBeamIsotope;
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
         for (int i = 0; i < 3; ++i) {
            line.erase(0, line.find_first_not_of(" ,\t"));
            auto NextDelim = line.find_first_of(" ,\t");
            if (NextDelim == std::string::npos)
               NextDelim = line.length();
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
               break;
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
   if (CS_InterpolationSpline)
      gsl_spline_free(CS_InterpolationSpline);
   if (CS_InterpolationAccellerator)
      gsl_interp_accel_free(CS_InterpolationAccellerator);
   CS_InterpolationAccellerator = gsl_interp_accel_alloc();
   CS_InterpolationSpline = gsl_spline_alloc(gsl_interp_linear, CrossSection.size());
   //   CS_InterpolationSpline = gsl_spline_alloc(gsl_interp_cspline, CrossSection.size());
   //   CS_InterpolationSpline = gsl_spline_alloc(gsl_interp_akima, CrossSection.size());
   //   CS_InterpolationSpline = gsl_spline_alloc(gsl_interp_steffen, CrossSection.size());

   gsl_spline_init(CS_InterpolationSpline, x, y, CrossSection.size());
}

double Reaction::GetCrossSection(double En) {
   if (CrossSection.size() == 0)
      ReadR33File(true);

   if (!CS_InterpolationAccellerator)
      InitializeInterpolation();

   if (En < GetEmin() || En > GetEmax())
      return 0.;

   return gsl_spline_eval(CS_InterpolationSpline, En, CS_InterpolationAccellerator);
}

double StragglingCrossSectionWrap(double z, void *user_data) {
   Reaction *this_ptr = (Reaction *)user_data;
   return this_ptr->StragglingCrossSection(z);
}

double Reaction::StragglingCrossSection(double E) {
   double f = (1 / (StragglingSigma * sqrt(2 * M_PI))) * exp(-pow((E - StragglingEn), 2) / (2 * pow(StragglingSigma, 2)));
   f *= GetCrossSection(E);
   //   f /= Stopping::GetInstance()->GetStopping(E, ReactionBeamIsotope, YieldTarget, YieldTargetLayer);
   return f;
}

double Reaction::GetStragglingCrossSection(double En, double InitialBeamEnergy, double InitialBeamEnegySpread, Target *target) {
   //   double sigma = Stopping::GetInstance()->GetStraggling(InitialBeamEnergy, InitialBeamEnegySpread, InitialBeamEnergy - En, ReactionBeamIsotope, target);
   double sigma = Stopping::GetInstance()->GetStraggling(InitialBeamEnergy, InitialBeamEnegySpread, En, ReactionBeamIsotope, target);
   return GetStragglingCrossSection(En, sigma);
}

double Reaction::GetStragglingCrossSection(double En, double sigma) {
   if (sigma == 0.)
      return GetCrossSection(En);

   gsl_function F;
   SetStragglingParameters(En, sigma);
   F.function = &StragglingCrossSectionWrap;
   F.params = (double *)this;

   double result;
   double error;
   gsl_integration_qag(&F, En - (100 * sigma), En + (100 * sigma),
                       IntegrationepsabsDefault, IntegrationepsrelDefault, IntegrationWorkspaceSize, QAGIntegrationKey,
                       StragglingIntegrationWorkspace, &result, &error);
   //   gsl_integration_qagi(&F, IntegrationepsabsDefault, IntegrationepsrelDefault, IntegrationWorkspaceSize, StragglingIntegrationWorkspace, &result, &error);
   return result;
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

double StragglingReactionYieldWrap(double z, void *user_data) {
   Reaction *this_ptr = (Reaction *)user_data;
   return this_ptr->StragglingYieldFunction(z);
}

double Reaction::StragglingYieldFunction(double En) {
   double sigma = Stopping::GetInstance()->GetStraggling(StragglingEin, YieldEnergySigma, En, ReactionBeamIsotope, YieldTarget, YieldTargetLayer);
   double f = YieldMultFactor * GetStragglingCrossSection(En, sigma) / Stopping::GetInstance()->GetStopping(En, ReactionBeamIsotope, YieldTarget, YieldTargetLayer);
   //   double f = YieldMultFactor * GetStragglingCrossSection(En, sigma);
   return f;
}

void Reaction::SetIntegrationParameters(double WorkSpaceSize, double epsabs, double epsrel, double epsrelMult, int QAGkey) {
   if (IntegrationWorkspaceSize != WorkSpaceSize || !YieldIntegrationWorkspace) {
      IntegrationWorkspaceSize = WorkSpaceSize;
      gsl_integration_workspace_free(YieldIntegrationWorkspace);
      YieldIntegrationWorkspace = gsl_integration_workspace_alloc(IntegrationWorkspaceSize);
      gsl_integration_workspace_free(StragglingIntegrationWorkspace);
      StragglingIntegrationWorkspace = gsl_integration_workspace_alloc(IntegrationWorkspaceSize);
      gsl_integration_workspace_free(YieldCalculationWorkspace);
      YieldCalculationWorkspace = gsl_integration_workspace_alloc(IntegrationWorkspaceSize);
      gsl_integration_workspace_free(StragglingCrossSectionOverStoppingWorkspace);
      StragglingCrossSectionOverStoppingWorkspace = gsl_integration_workspace_alloc(IntegrationWorkspaceSize);
   }
   IntegrationepsabsDefault = epsabs;
   IntegrationepsrelDefault = epsrel;
   IntegrationEpsRelMultFactor = epsrelMult;
   QAGIntegrationKey = QAGkey;
   return;
}

void Reaction::SetIntegratedYield(Target *target, int layer, double sigma) {
   YieldTarget = target;
   YieldTargetLayer = layer;
   YieldMultFactor = YieldTarget->GetElementAtomicPercentInLayer(ReactionTargetElement, layer) * ReactionTargetIsotope->GetAbundance();
   YieldMultFactor /= (1.E19 * GSL_CONST_MKSA_ELECTRON_CHARGE);
   //   YieldMultFactor *= .624;
   YieldEnergySigma = sigma;

   SetIntegrationParameters(IntegrationWorkspaceSize, IntegrationepsabsDefault, IntegrationepsrelDefault, IntegrationEpsRelMultFactor, QAGIntegrationKey);
   return;
}
void Reaction::SetSIntegratedYield(Target *target, int layer, double ein, double sigma) {
   YieldTarget = target;
   YieldTargetLayer = layer;
   YieldMultFactor = YieldTarget->GetElementAtomicPercentInLayer(ReactionTargetElement, layer) * ReactionTargetIsotope->GetAbundance();
   YieldMultFactor /= (1.E19 * GSL_CONST_MKSA_ELECTRON_CHARGE);
   //   YieldMultFactor *= .624;
   YieldEnergySigma = sigma;
   StragglingEin = ein;

   SetIntegrationParameters(IntegrationWorkspaceSize, IntegrationepsabsDefault, IntegrationepsrelDefault, IntegrationEpsRelMultFactor, QAGIntegrationKey);
   return;
}

void Reaction::SetStragglingParameters(double En, double sigma) {
   StragglingEn = En;
   StragglingSigma = sigma;
   SetIntegrationParameters(IntegrationWorkspaceSize, IntegrationepsabsDefault, IntegrationepsrelDefault, IntegrationEpsRelMultFactor, QAGIntegrationKey);
}

double Reaction::GetYield(double Emin, double Emax, double &error, double &epsrel, Target *target) {
   if (Emax < GetEmin() ||
       Emax < Emin)
      return 0.;

   Emin = std::max(Emin, GetEmin());
   Emax = std::min(Emax, GetEmax());

   SetIntegratedYield(target);

   gsl_function F;
   F.function = &ReactionYieldWrap;
   F.params = (void *)this;

   epsrel = IntegrationepsrelDefault / IntegrationEpsRelMultFactor;

   int status = 1;
   gsl_set_error_handler_off();
   double result;
   while (status) {
      epsrel *= IntegrationEpsRelMultFactor;
      //      status = gsl_integration_qags(&F, CrossSectionEnergy[0], En,
      //                                    epsabs, epsrel, limit, gsl_workspace, &result, &error);
      status = gsl_integration_qag(&F, Emin, Emax,
                                   IntegrationepsabsDefault, epsrel, IntegrationWorkspaceSize, QAGIntegrationKey,
                                   YieldIntegrationWorkspace, &result, &error);
   }
   return result;
}

// double Reaction::GetYield(double Emin, double Emax, double InitialBeamEnegySpread, double &error, double &epsrel, Target *target) {
//    if (Emax < GetEmin() ||
//        Emax < Emin)
//       return 0.;

//    Emin = std::max(Emin, GetEmin());
//    Emax = std::min(Emax, GetEmax());

//    SetSIntegratedYield(target, 0, Emax, InitialBeamEnegySpread);

//    gsl_function F;
//    F.function = &StragglingReactionYieldWrap;
//    F.params = (void *)this;

//    epsrel = IntegrationepsrelDefault / IntegrationEpsRelMultFactor;

//    int status = 1;
//    gsl_set_error_handler_off();
//    double result;
//    while (status) {
//       epsrel *= IntegrationEpsRelMultFactor;
//       //      status = gsl_integration_qags(&F, CrossSectionEnergy[0], En,
//       //                                    epsabs, epsrel, limit, gsl_workspace, &result, &error);
//       status = gsl_integration_qag(&F, Emin, Emax,
//                                    IntegrationepsabsDefault, epsrel, IntegrationWorkspaceSize, QAGIntegrationKey,
//                                    YieldIntegrationWorkspace, &result, &error);
//       //      status = gsl_integration_qags(&F, Emin, Emax,
//       //                                    IntegrationepsabsDefault, epsrel, IntegrationWorkspaceSize,
//       //                                    YieldIntegrationWorkspace, &result, &error);
//       //   auto work = gsl_integration_romberg_alloc(20);
//       //   size_t nevals = 0;
//       //   status = gsl_integration_romberg(&F, Emin, Emax,
//       //                                    IntegrationepsabsDefault, epsrel,
//       //                                    &result, &nevals, work);
//       //   gsl_integration_romberg_free(work);
//    }
//    return result;
// }

void Reaction::SetYieldCalculationParameters(Target *target, int layer, double ein, double sigma) {
   YieldCalculationTarget = target;
   YieldCalculationLayer = layer;
   YieldCalculationMultFactor = target->GetElementAtomicPercentInLayer(ReactionTargetElement, layer) * ReactionTargetIsotope->GetAbundance();
   YieldCalculationMultFactor /= (1.E19 * GSL_CONST_MKSA_ELECTRON_CHARGE);
   YieldCalculationInitialSigma = sigma;
   YieldCalculationEin = ein;

   SetIntegrationParameters(IntegrationWorkspaceSize, IntegrationepsabsDefault, IntegrationepsrelDefault, IntegrationEpsRelMultFactor, QAGIntegrationKey);
}

double Reaction::YieldCalculation(double En) {
   double CSoverStopping;
   if (YieldCalculationInitialSigma < 0)
      CSoverStopping = GetCrossSection(En) /
                       Stopping::GetInstance()->GetStopping(
                           En, ReactionBeamIsotope, YieldCalculationTarget, YieldCalculationLayer);
   else {
      double sigma = Stopping::GetInstance()->GetStraggling(
          YieldCalculationEin, YieldCalculationInitialSigma, En,
          ReactionBeamIsotope, YieldCalculationTarget, YieldCalculationLayer);
      // CSoverStopping = GetStragglingCrossSection(En, sigma) /
      //                  Stopping::GetInstance()->GetStopping(
      //                      En, ReactionBeamIsotope, YieldCalculationTarget, YieldCalculationLayer);
      CSoverStopping = GetStragglingCrossSectionOverStopping(En, sigma);
   }
   return YieldCalculationMultFactor * CSoverStopping;
}

double YieldCalculationWrap(double En, void *user_data) {
   Reaction *this_ptr = (Reaction *)user_data;
   return this_ptr->YieldCalculation(En);
}

double StragglingCrossSectionOverStoppingWrap(double z, void *user_data) {
   Reaction *this_ptr = (Reaction *)user_data;
   return this_ptr->StragglingCrossSectionOverStopping(z);
}

double Reaction::StragglingCrossSectionOverStopping(double En) {
   double f = (1 / (YieldCalculationCurrentSigma * sqrt(2 * M_PI))) * exp(-pow((En - YieldCalculationCurrentE), 2) / (2 * pow(YieldCalculationCurrentSigma, 2)));
   f *= GetCrossSection(En);
   f /= Stopping::GetInstance()->GetStopping(En, ReactionBeamIsotope, YieldCalculationTarget, YieldCalculationLayer);
   return f;
}

double Reaction::GetStragglingCrossSectionOverStopping(double En, double sigma) {
   if (sigma == 0.)
      return GetCrossSection(En);

   gsl_function F;

   YieldCalculationCurrentE = En;
   YieldCalculationCurrentSigma = sigma;
   SetIntegrationParameters(IntegrationWorkspaceSize, IntegrationepsabsDefault, IntegrationepsrelDefault, IntegrationEpsRelMultFactor, QAGIntegrationKey);

   F.function = &StragglingCrossSectionOverStoppingWrap;
   F.params = (double *)this;

   double result;
   double error;
   // gsl_integration_qag(&F, En - (100 * sigma), En + (100 * sigma),
   //                     IntegrationepsabsDefault, IntegrationepsrelDefault, IntegrationWorkspaceSize, QAGIntegrationKey,
   //                     StragglingCrossSectionOverStoppingWorkspace, &result, &error);
   gsl_integration_qag(&F, En - (10 * sigma), En + (10 * sigma),
                       0., 1.e-9, 10000, 1,
                       StragglingCrossSectionOverStoppingWorkspace, &result, &error);

   return result;

   //   return result / Stopping::GetInstance()->GetStopping(En, ReactionBeamIsotope, YieldCalculationTarget, YieldCalculationLayer);
}

double Reaction::GetYield(double Emin, double Emax, double InitialBeamEnegySpread, double &error, double &epsrel, Target *target) {
   if (Emax < GetEmin() ||
       Emax < Emin ||
       (Emin == Emax && InitialBeamEnegySpread <= 0.))
      return 0.;

   Emin = std::max(Emin, GetEmin());
   Emax = std::min(Emax, GetEmax());

   double yield = 0.;

   for (int Layer = 0; Layer < target->GetNumberOfLayers(); ++Layer) {
      double LayerEout = Stopping::GetInstance()->GetEout(Emax, ReactionBeamIsotope, target, Layer);
      double LayerEmin = std::max(Emin, LayerEout);
      double LayerEmax = Emax;
      double LayerSigma = InitialBeamEnegySpread;
      SetYieldCalculationParameters(target, Layer, LayerEmax, LayerSigma);

      if (LayerEmin == LayerEmax) {
         return YieldCalculationMultFactor * GetStragglingCrossSectionOverStopping(LayerEmin, LayerSigma);
      }

      gsl_function F;
      F.function = &YieldCalculationWrap;
      F.params = (void *)this;

      epsrel = IntegrationepsrelDefault / IntegrationEpsRelMultFactor;

      int status = 1;
      gsl_set_error_handler_off();
      double result;
      while (status) {
         epsrel *= IntegrationEpsRelMultFactor;
         status = gsl_integration_qag(&F, Emin, Emax,
                                      IntegrationepsabsDefault, epsrel, IntegrationWorkspaceSize, QAGIntegrationKey,
                                      YieldIntegrationWorkspace, &result, &error);
      }
      yield += result;
   }
   return yield;
}

double Reaction::GetYield(double En, double dE, double AtomicPerCent, double AtomicAbundunce, double Stopping, double enSigma) {
   if (En < GetEmin() || En > GetEmax())
      return 0.;

   double yield = GetStragglingCrossSection(En, enSigma) * dE * 1e-27 * 1e15 * 6.24E12 * 1000 * (AtomicPerCent / 100.) * (AtomicAbundunce / 100);
   return yield / Stopping;
}
