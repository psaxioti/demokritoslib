#include "Isotope.hh"

static double amu_to_kev = 931494.061;

Isotope::Isotope(int Z, int N, std::string Symb, double mass, double Abund, double massexc) : AtomicNumber(Z), NeutronNumber(N), Mass(mass), Abundance(Abund), NaturalAbundance(Abund) {
   Symbol = std::to_string(AtomicNumber + NeutronNumber) + Symb;
   if (massexc == 0.)
      MassExcess = (Mass - AtomicNumber - NeutronNumber) * amu_to_kev;
   else
      MassExcess = massexc;
}

Isotope::~Isotope() {}

double Isotope::GetAbundance() {
   return Abundance;
}

double Isotope::GetNaturalAbundance() {
   return NaturalAbundance;
}

void Isotope::SetAbundance(double NewAbund) {
   if (NewAbund < 0. || NewAbund > 100.)
      return;
   else
      Abundance = NewAbund;
   return;
}

double Isotope::GetMass() {
   return Mass;
}

double Isotope::GetMassExcess() {
   return MassExcess;
}

std::string Isotope::GetSymbol() {
   return Symbol;
}

int Isotope::GetZ() {
   return AtomicNumber;
}

int Isotope::GetN() {
   return NeutronNumber;
}

int Isotope::GetA() {
   return AtomicNumber + NeutronNumber;
}
