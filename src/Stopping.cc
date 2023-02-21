#include "Stopping.hh"

#include <cmath>
#include <istream>
#include <streambuf>

#include "MassFunctions.hh"

extern char zbla_start[] asm("_binary_data_scoef_95a_start");
extern char zbla_end[] asm("_binary_data_scoef_95a_end");
extern char zblb_start[] asm("_binary_data_scoef_95b_start");
extern char zblb_end[] asm("_binary_data_scoef_95b_end");
extern char srim_start[] asm("_binary_data_srim_txt_start");
extern char srim_end[] asm("_binary_data_srim_txt_end");

struct membuf : std::streambuf {
   membuf(char *begin, char *end) {
      this->setg(begin, begin, end);
   }
};

Stopping::Stopping() {
   StoppingName = "ZBL";
   ReadZBL_Coeffs();
   ReadSRIM_Coeffs();
}

Stopping::~Stopping() {
   if (spline)
      gsl_spline_free(spline);
   if (acc)
      gsl_interp_accel_free(acc);
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

void Stopping::SetStopping(std::string &Stopping) {
   if (Stopping == "ZBL" || Stopping == "SRIM")
      StoppingName = Stopping;
}

double Stopping::GetStraggling(double En, double InitialEnSpread, double LayerDE, Isotope *Beam, Element *TargetEl) {
   return std::sqrt(GetStraggling2(En, InitialEnSpread, LayerDE, Beam, TargetEl));
}

double Stopping::GetStraggling(double En, double InitialEnSpread, double LayerDE, Isotope *Beam, Target *target, int Layer) {
   double straggling2 = 0.;
   auto LayerElements = target->GetElementsInLayer(Layer);
   for (size_t i = 0; i < LayerElements.size(); ++i)
      straggling2 += (target->GetElementAtomicPercentInLayer(LayerElements[i], Layer) * GetStraggling2(En, InitialEnSpread, LayerDE, Beam, LayerElements[i]));
   return std::sqrt(straggling2 / target->GetLayerAtomicPerCent(Layer));
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
   double effectivecharge = BeamIso->GetZ();
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
      gsl_interp_accel_reset(acc);
      gsl_spline_init(spline, x, y, SRIM_Energies.size());
   }
   return gsl_spline_eval(spline, En / BeamIso->GetMass(), acc);
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

   double x_min, x_max;
   int x_points;

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
   for (size_t i = 0; i < x_points; ++i) {
      double kk = std::log10(x_min) + i * (log10(x_max) - log10(x_min)) / (x_points - 1);
      SRIM_Energies.push_back(std::pow(10., kk));
   }
   acc = gsl_interp_accel_alloc();
   spline = gsl_spline_alloc(gsl_interp_linear, x_points);
}

double Stopping::GetStraggling2(double En, double InitialEnSpread, double LayerDE, Isotope *Beam, Element *TargetEl) {
   double Initial_Stopping = GetStopping(En, Beam, TargetEl);
   double s_non2 = GetNonStatisticalStraggling2(Initial_Stopping, En - LayerDE, InitialEnSpread, Beam, TargetEl);
   double s_el2 = GetChuFactor(En, Beam, TargetEl) * GetBohrElectronicStraggling2(LayerDE / Initial_Stopping, Beam, TargetEl);
   return (s_non2 + s_el2);
}

double Stopping::GetNonStatisticalStraggling2(double IncomingStopping, double OutgoingEn, double InitialEnSpread, Isotope *Beam, Element *TargetEl) {
   double Straggling = IncomingStopping * InitialEnSpread / GetStopping(OutgoingEn, Beam, TargetEl);
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
   return 1.;
}
