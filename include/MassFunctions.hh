#ifndef MASSFUNCTIONS_H
#define MASSFUNCTIONS_H 1

#include <vector>

#include "Element.hh"
#include "Isotope.hh"

/// Get the vector with the elements
std::vector<Element *> GetElements();
/// Get the vector with the isotopes
std::vector<Isotope *> GetIsotopes();

/// Get the element by symbol (eg Au)
Element *GetElement(const std::string &Symbol);
/// Get the element by atomic number
Element *GetElement(const int &Z);

/// Get the isotope by symbol (eg 197Au)
Isotope *GetIsotope(const std::string &Symbol);
/// Get the isotope by atomic and mass number
Isotope *GetIsotope(const int &Z, const int &A);

#endif
