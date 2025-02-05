#include "Stopping.hh"

#include <algorithm>
#include <cmath>
#include <istream>
#include <streambuf>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

#include "MassFunctions.hh"

extern char zbla_start[] asm("_binary_data_scoef_95a_start");
extern char zbla_end[] asm("_binary_data_scoef_95a_end");
extern char zblb_start[] asm("_binary_data_scoef_95b_start");
extern char zblb_end[] asm("_binary_data_scoef_95b_end");
extern char srim_start[] asm("_binary_data_srim_txt_start");
extern char srim_end[] asm("_binary_data_srim_txt_end");
extern char chu_start[] asm("_binary_data_chu_dat_start");
extern char chu_end[] asm("_binary_data_chu_dat_end");

struct membuf : std::streambuf {
   membuf(char *begin, char *end) {
      this->setg(begin, begin, end);
   }
};

Stopping *Stopping::StoppingInstance = nullptr;

Stopping *Stopping::GetInstance() {
   if (!StoppingInstance)
      StoppingInstance = new Stopping();
   return StoppingInstance;
}

Stopping::Stopping() {
   StoppingName = "SRIM";
   ReadZBL_Coeffs();
   ReadSRIM_Coeffs();
   ReadCHU_coeffs();
}

Stopping::~Stopping() {
   if (srim_spline)
      gsl_spline_free(srim_spline);
   if (srim_acc)
      gsl_interp_accel_free(srim_acc);
   if (chu_spline)
      gsl_spline_free(chu_spline);
   if (chu_acc)
      gsl_interp_accel_free(chu_acc);
}

double Stopping::GetStopping(double En, Isotope *BeamIso, Element *TargetEl) {
   if (StoppingName == "ZBL")
      return ZBL_Stopping(En, BeamIso, TargetEl);
   else
      return SRIM_Stopping(En, BeamIso, TargetEl);
}

double Stopping::GetStopping(double En, Isotope *BeamIso, Target *target, int Layer) {
   auto LayerElements = target->GetElementsInLayer(Layer);
   double StoppingSum = 0.;
   for (size_t i = 0; i < LayerElements.size(); ++i)
      StoppingSum += target->GetElementAtomicPercentInLayer(LayerElements[i], Layer) * GetStopping(En, BeamIso, LayerElements[i]);
   return (StoppingSum / target->GetLayerAtomicPerCent(Layer));
}

void Stopping::SetStopping(const std::string &Stopping) {
   if (Stopping == "ZBL" || Stopping == "SRIM")
      StoppingName = Stopping;
}

double Stopping::GetDEStraggling(double En, double InitialEnSpread, double LayerDE, Isotope *Beam, Target *target, int Layer) {
   double straggling2 = GetNonStatisticalStraggling2(En, InitialEnSpread, LayerDE, Beam, target, Layer);
   auto LayerElements = target->GetElementsInLayer(Layer);
   for (size_t i = 0; i < LayerElements.size(); ++i) {
      double nonstatisticalstraggling2 = GetStatisticalStraggling2(En, LayerDE, Beam, LayerElements[i]);
      straggling2 += target->GetElementAtomicPercentInLayer(LayerElements[i], Layer) * nonstatisticalstraggling2 / target->GetLayerAtomicPerCent(Layer);
   }
   return std::sqrt(straggling2);
}

double Stopping::ZBL_Stopping(double En, Isotope *BeamIso, Element *TargetEl) {
   //   return ZBL_ElectronicStopping(En, BeamIso, TargetEl);
   return (ZBL_ElectronicStopping(En, BeamIso, TargetEl) + ZBL_NuclearStopping(En, BeamIso, TargetEl));
}

double Stopping::ZBL_ElectronicStopping(double En, Isotope *BeamIso, Element *TargetEl) {
   double EnPerMass = En / BeamIso->GetMass();
   double stopping = 0.;
   // Have to implement something better
   // Just added it here very roughly
   double effectivecharge = BeamIso->GetZ() * BeamIso->GetZ();
   if (EnPerMass <= 10.) {
      double y = (TargetEl->GetZ() > 6) ? 0.45 : 0.35;
      stopping = ZBL_ElectronicStopping(10. * BeamIso->GetMass(), BeamIso, TargetEl) * std::pow(EnPerMass / 10., y);
   } else if (EnPerMass < 10000.) {
      double Slow = (ZBL_coeffs[TargetEl][0] * pow(EnPerMass, ZBL_coeffs[TargetEl][1])) + (ZBL_coeffs[TargetEl][2] * pow(EnPerMass, ZBL_coeffs[TargetEl][3]));
      double Shigh = ZBL_coeffs[TargetEl][4] * log((ZBL_coeffs[TargetEl][6] / EnPerMass) + (ZBL_coeffs[TargetEl][7] * EnPerMass)) / pow(EnPerMass, ZBL_coeffs[TargetEl][5]);
      stopping = Slow * Shigh / (Slow + Shigh);
   } else {
      double x = std::log(EnPerMass) / EnPerMass;
      stopping = ZBL_coeffs[TargetEl][8] + (ZBL_coeffs[TargetEl][9] * x) + (ZBL_coeffs[TargetEl][10] * std::pow(x, 2)) + (ZBL_coeffs[TargetEl][11] / x);
   }
   return effectivecharge * stopping;
}

double Stopping::ZBL_NuclearStopping(double En, Isotope *BeamIso, Element *TargetEl) {
   double epsilon = BeamIso->GetZ() * TargetEl->GetZ() * (BeamIso->GetMass() + TargetEl->GetMass());
   epsilon *= (std::pow(BeamIso->GetZ(), .23) + std::pow(TargetEl->GetZ(), .23));
   epsilon = (32.53 * TargetEl->GetMass() * En) / epsilon;
   double reducedNuclearStopping;
   if (epsilon <= 30.)
      reducedNuclearStopping = std::log1p(1.1383 * epsilon) / (2 * (epsilon + .01321 * std::pow(epsilon, .21226) + (.19593 * std::pow(epsilon, .5))));
   else
      reducedNuclearStopping = std::log(epsilon) / (2 * epsilon);
   double stopping = 8.462 * BeamIso->GetZ() * TargetEl->GetZ() * BeamIso->GetMass() * reducedNuclearStopping;
   stopping /= ((BeamIso->GetMass() + TargetEl->GetMass()) * (std::pow(BeamIso->GetZ(), .23) + std::pow(TargetEl->GetZ(), .23)));
   return stopping;
}

double Stopping::SRIM_Stopping(double En, Isotope *BeamIso, Element *TargetEl) {
   if (BeamIso->GetZ() < SRIM_ZBeam_min || BeamIso->GetZ() > SRIM_ZBeam_max ||
       TargetEl->GetZ() < SRIM_ZTarget_min || TargetEl->GetZ() > SRIM_ZTarget_max)
      return ZBL_Stopping(En, BeamIso, TargetEl);
   if (SRIM_Previous_ZBeam != BeamIso->GetZ() || SRIM_Previous_ZTarget != TargetEl->GetZ()) {
      double *x = &SRIM_Energies[0];
      double *y = &SRIM_coeffs[BeamIso->GetZ()][TargetEl->GetZ()][0];
      gsl_interp_accel_reset(srim_acc);
      gsl_spline_init(srim_spline, x, y, SRIM_Energies.size());
      SRIM_Previous_ZBeam = BeamIso->GetZ();
      SRIM_Previous_ZTarget = TargetEl->GetZ();
   }
   if (En / BeamIso->GetMass() < SRIM_Energies[0])
      return SRIM_coeffs[BeamIso->GetZ()][TargetEl->GetZ()][0];
   return gsl_spline_eval(srim_spline, En / BeamIso->GetMass(), srim_acc);
}

void Stopping::ReadZBL_Coeffs() {
   membuf abuf(zbla_start, zbla_end);
   std::istream afile(&abuf);

   membuf bbuf(zblb_start, zblb_end);
   std::istream bfile(&bbuf);

   std::string line, dummy;
   int Z;
   double A1, A2, A3, A4, A5, A6, A7, A8;
   double A9, A10, A11, A12;
   double mass;

   getline(afile, line);
   getline(afile, line);
   getline(bfile, line);
   getline(bfile, line);
   while (afile >> Z >> dummy >> dummy >> mass >> dummy >> dummy >> dummy >> dummy >> A1 >> A2 >> A3 >> A4 >> A5 >> A6 >> A7 >> A8) {
      bfile >> A9 >> A10 >> A11 >> A12 >>
          dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >>
          dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >>
          dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >>
          dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy;
      if (Z < 1)
         continue;
      Element *el = GetElement(Z);
      el->SetMass(mass);
      ZBL_coeffs[el][0] = A1;
      ZBL_coeffs[el][1] = A2;
      ZBL_coeffs[el][2] = A3;
      ZBL_coeffs[el][3] = A4;
      ZBL_coeffs[el][4] = A5;
      ZBL_coeffs[el][5] = A6;
      ZBL_coeffs[el][6] = A7;
      ZBL_coeffs[el][7] = A8;
      ZBL_coeffs[el][8] = A9;
      ZBL_coeffs[el][9] = A10;
      ZBL_coeffs[el][10] = A11;
      ZBL_coeffs[el][11] = A12;
   }
}

void Stopping::ReadSRIM_Coeffs() {
   membuf abuf(srim_start, srim_end);
   std::istream afile(&abuf);

   double x_min = 1, x_max = 1;
   int x_points = 0;

   while (!afile.eof()) {
      std::string line;
      std::getline(afile, line);

      auto equal_char = line.find('=');
      std::string read_val = line.substr(0, equal_char);

      if (read_val == "z1-min")
         SRIM_ZBeam_min = std::stod(line.substr(equal_char + 1));
      else if (read_val == "z1-max")
         SRIM_ZBeam_max = std::stod(line.substr(equal_char + 1));
      else if (read_val == "z2-min")
         SRIM_ZTarget_min = std::stod(line.substr(equal_char + 1));
      else if (read_val == "z2-max")
         SRIM_ZTarget_max = std::stod(line.substr(equal_char + 1));
      else if (read_val == "x-min")
         x_min = std::stod(line.substr(equal_char + 1));
      else if (read_val == "x-max")
         x_max = std::stod(line.substr(equal_char + 1));
      else if (read_val == "x-points")
         x_points = std::stoi(line.substr(equal_char + 1));
      else if (line.substr(0, 9) == "#STOPPING") {
         int z1_current, z2_current;
         line.erase(0, line.find('=') + 1);
         z1_current = std::stoi(line.substr(0, line.find(' ')));
         line.erase(0, line.find('=') + 1);
         z2_current = std::stoi(line);
         for (size_t i = 0; i < x_points; ++i) {
            std::getline(afile, line);
            SRIM_coeffs[z1_current][z2_current].push_back(std::stod(line));
         }
      }
   }

   double logStep = (log10(x_max) - log10(x_min)) / (x_points - 1);
   double logEner = log10(x_min);
   SRIM_Energies.clear();
   while (SRIM_Energies.size() < x_points) {
      SRIM_Energies.push_back(std::pow(10., logEner));
      logEner += logStep;
   }
   srim_acc = gsl_interp_accel_alloc();
   srim_spline = gsl_spline_alloc(gsl_interp_linear, x_points);
   //   srim_spline = gsl_spline_alloc(gsl_interp_steffen, x_points);
}

void Stopping::ReadCHU_coeffs() {
   membuf abuf(chu_start, chu_end);
   std::istream file(&abuf);

   int Z;
   double A1, A2, A3, A4, A5, A6, A7, A8;
   double A9, A10, A11, A12, A13, A14;
   int nb_points = 0;

   file >> A1 >> A2 >> A3 >> A4 >> A5 >> A6 >> A7 >> A8 >> A9 >> A10 >> A11 >> A12 >> A13 >> A14;
   CHU_EnergyPerMass.push_back(A1);
   CHU_EnergyPerMass.push_back(A2);
   CHU_EnergyPerMass.push_back(A3);
   CHU_EnergyPerMass.push_back(A4);
   CHU_EnergyPerMass.push_back(A5);
   CHU_EnergyPerMass.push_back(A6);
   CHU_EnergyPerMass.push_back(A7);
   CHU_EnergyPerMass.push_back(A8);
   CHU_EnergyPerMass.push_back(A9);
   CHU_EnergyPerMass.push_back(A10);
   CHU_EnergyPerMass.push_back(A11);
   CHU_EnergyPerMass.push_back(A12);
   CHU_EnergyPerMass.push_back(A13);
   CHU_EnergyPerMass.push_back(A14);
   if (CHU_EnergyPerMass[0] > CHU_EnergyPerMass[1])
      std::reverse(CHU_EnergyPerMass.begin(), CHU_EnergyPerMass.end());

   while (file >> Z >> CHU_coeffs[Z][A1] >> CHU_coeffs[Z][A2] >>
          CHU_coeffs[Z][A3] >> CHU_coeffs[Z][A4] >>
          CHU_coeffs[Z][A5] >> CHU_coeffs[Z][A6] >>
          CHU_coeffs[Z][A7] >> CHU_coeffs[Z][A8] >>
          CHU_coeffs[Z][A9] >> CHU_coeffs[Z][A10] >>
          CHU_coeffs[Z][A11] >> CHU_coeffs[Z][A12] >>
          CHU_coeffs[Z][A13] >> CHU_coeffs[Z][A14])
      nb_points++;

   chu_acc = gsl_interp_accel_alloc();
   chu_spline = gsl_spline_alloc(gsl_interp_linear, CHU_EnergyPerMass.size());
}

double Stopping::GetStatisticalStraggling2(double En, double LayerDE, Isotope *Beam, Element *TargetEl) {
   double Initial_Stopping = GetStopping(En, Beam, TargetEl);
   double s_el2 = GetChuFactor(En, Beam, TargetEl) * GetBohrElectronicStraggling2(LayerDE / Initial_Stopping, Beam, TargetEl);
   return (s_el2);
}

double Stopping::GetNonStatisticalStraggling2(double En, double InitialEnSpread, double LayerDE, Isotope *Beam, Target *target, int Layer) {
   double Initial_Stopping = GetStopping(En, Beam, target, Layer);
   double Final_Stopping = GetStopping(En - LayerDE, Beam, target, Layer);
   double Straggling = Final_Stopping * InitialEnSpread / Initial_Stopping;
   return std::pow(Straggling, 2);
}

double Stopping::GetBohrElectronicStraggling2(double LayerThickness, Isotope *Beam, Element *TargetEl) {
   return (0.26 * std::pow(Beam->GetZ(), 2) * TargetEl->GetZ() * LayerThickness);
}

double Stopping::GetBohrNuclearStraggling2(double LayerThickness, Isotope *Beam, Element *TargetEl) {
   double straggling = GetBohrElectronicStraggling2(LayerThickness, Beam, TargetEl);
   straggling *= (TargetEl->GetZ() * std::pow(Beam->GetA() / (Beam->GetMass() + TargetEl->GetMass()), 2));
   return straggling;
}

double Stopping::GetChuFactor(double En, Isotope *Beam, Element *TargetEl) {
   if (CHU_Previous_ZTarget != TargetEl->GetZ()) {
      double *x = &CHU_EnergyPerMass[0];
      std::vector<double> coeffs;
      coeffs.clear();
      for (const auto &en : CHU_EnergyPerMass)
         coeffs.push_back(CHU_coeffs[TargetEl->GetZ()][en]);
      double *y = &coeffs[0];
      gsl_interp_accel_reset(chu_acc);
      gsl_spline_init(chu_spline, x, y, CHU_EnergyPerMass.size());
      CHU_Previous_ZTarget = TargetEl->GetZ();
   }
   if (En / Beam->GetMass() < CHU_EnergyPerMass[0])
      return 1.;
   return gsl_spline_eval(chu_spline, En / Beam->GetMass(), chu_acc);
}

double dxdeWrapper(double z, void *user_data) {
   Stopping *this_ptr = (Stopping *)user_data;
   return this_ptr->dxde(z);
}

double Stopping::dxde(double en) {
   auto res = GetStopping(en, ThicknessCalculationBeam, ThicknessCalculationTarget, ThicknessCalculationLayer);
   return (1 / res);
}

double Stopping::GetRange(double En, Isotope *Beam, Target *target, int Layer) {
   return GetTraveledThickness(0., En, Beam, target, Layer);
}

double Stopping::GetTraveledThickness(double Eout, double Ein, Isotope *Beam, Target *target, int Layer) {
   ThicknessCalculationBeam = Beam;
   ThicknessCalculationTarget = target;
   ThicknessCalculationLayer = Layer;
   ThicknessCalculationEin = Ein;

   auto dxdeIntegrationWorkspace = gsl_integration_workspace_alloc(IntegrationWorkspaceSize);

   gsl_function F;
   F.function = &dxdeWrapper;
   F.params = (double *)this;

   double result;
   double error;
   gsl_integration_qag(&F, Eout, Ein,
                       IntegrationepsabsDefault, IntegrationepsrelDefault, IntegrationWorkspaceSize, QAGIntegrationKey,
                       dxdeIntegrationWorkspace, &result, &error);
   gsl_integration_workspace_free(dxdeIntegrationWorkspace);
   return result;
}

double thicknessRootWrapper(double z, void *user_data) {
   Stopping *this_ptr = (Stopping *)user_data;
   return this_ptr->thicknessRoot(z);
}

double Stopping::thicknessRoot(double en) {
   auto result = GetTraveledThickness(en, ThicknessCalculationEin, ThicknessCalculationBeam, ThicknessCalculationTarget, ThicknessCalculationLayer);
   result -= ThicknessCalculationTarget->GetLayerThickness(ThicknessCalculationLayer);
   return result;
}

double Stopping::GetEout(double En, Isotope *Beam, Target *target, int Layer) {
   if (GetRange(En, Beam, target, Layer) <= target->GetLayerThickness(Layer))
      return 0.;

   int status;
   int iter = 0;
   int max_iter = 100;
   double result = 0;

   const gsl_root_fsolver_type *RootFinderType;
   gsl_root_fsolver *RootFinder;
   gsl_function F;
   F.function = &thicknessRootWrapper;
   F.params = (double *)this;
   RootFinderType = gsl_root_fsolver_brent;
   //   RootFinderType = gsl_root_fsolver_bisection;
   RootFinder = gsl_root_fsolver_alloc(RootFinderType);
   gsl_root_fsolver_set(RootFinder, &F, 0, En);
   do {
      iter++;
      gsl_root_fsolver_iterate(RootFinder);
      result = gsl_root_fsolver_root(RootFinder);
      auto x_lo = gsl_root_fsolver_x_lower(RootFinder);
      auto x_hi = gsl_root_fsolver_x_upper(RootFinder);
      status = gsl_root_test_interval(x_lo, x_hi, 0, 1.E-9);
      //   status = gsl_root_test_residual(thick(result), 1.E-9);
   } while (status == GSL_CONTINUE && iter < max_iter);

   gsl_root_fsolver_free(RootFinder);

   return result;
}

double stragglingWrapper(double z, void *user_data) {
   Stopping *this_ptr = (Stopping *)user_data;
   return this_ptr->straggling(z);
}

double Stopping::straggling(double e) {
   return (GetChuFactor(e, StragglingCalculationBeam, StragglingCalculationElement) /
           std::pow(GetStopping(e, StragglingCalculationBeam, StragglingCalculationElement), 3));
}

double Stopping::IntegratedStraggling(double Ein, double s0, double Eout, Isotope *Beam, Element *el) {
   auto stragIntegrationWorkspace = gsl_integration_workspace_alloc(IntegrationWorkspaceSize);

   StragglingCalculationBeam = Beam;
   StragglingCalculationElement = el;

   gsl_function F;
   F.function = &stragglingWrapper;
   F.params = (double *)this;

   double result;
   double error;
   gsl_integration_qag(&F, Eout, Ein,
                       IntegrationepsabsDefault, IntegrationepsrelDefault, IntegrationWorkspaceSize, QAGIntegrationKey,
                       stragIntegrationWorkspace, &result, &error);

   result *= std::pow(GetStopping(Eout, Beam, el), 2) *
             0.26 * std::pow(Beam->GetZ(), 2) * el->GetZ();

   gsl_integration_workspace_free(stragIntegrationWorkspace);

   return result;
}

double Stopping::GetStraggling(double Ein, double s0, double Eout, Isotope *Beam, Target *target, int Layer) {
   double s2 = std::pow(GetStopping(Eout, Beam, target, Layer) * s0 / GetStopping(Ein, Beam, target, Layer), 2.);

   auto elem = target->GetElementsInLayer(Layer);
   for (const auto &el : elem)
      s2 += target->GetElementAtomicPercentInLayer(el, Layer) * IntegratedStraggling(Ein, s0, Eout, Beam, el) /
            target->GetLayerAtomicPerCent(Layer);

   return std::sqrt(s2);
}

double Stopping::CalculateStraggling(double Ein, double s0, double Eout, Isotope *Beam, Target *target, int Layer) {
   double curE = Ein;
   double s = s0;
   double de = .05;

   while (curE >= Eout + de) {
      s = GetDEStraggling(curE, s, de, Beam, target, Layer);
      if (curE - de < Eout + de && curE > Eout + de)
         de = curE - Eout - de;
      curE -= de;
   }
   return s;
}
