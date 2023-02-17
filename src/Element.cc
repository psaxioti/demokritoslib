#include "Element.hh"

#include <algorithm>

Element::Element(int Z, std::string Symbol) : AtomicNumber(Z), Symbol(Symbol) {
}

Element::~Element() {
}

void Element::AddIsotope(Isotope *Iso) {
   if (std::find(Isotopes.begin(), Isotopes.end(), Iso) == Isotopes.end())
      Isotopes.push_back(Iso);
   return;
}

std::vector<Isotope *> Element::GetElementIsotopes() {
   return Isotopes;
}

int Element::GetZ() {
   return AtomicNumber;
}

std::string Element::GetSymbol() {
   return Symbol;
}

void Element::SetMass(double Mass) {
   ElementMass = Mass;
}

double Element::GetMass() {
   return ElementMass;
}
