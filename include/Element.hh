#ifndef ELEMENT_H
#define ELEMENT_H 1

#include <string>
#include <vector>

#include "Isotope.hh"

class Element {
public:
   /**
    * Element class constructor.
    *
    * @param Z Atomic number of element.
    * @param Symbol Symbol of the element in the form of eg. Au.
    */
   Element(int Z, std::string Symbol);
   ~Element();

   /// Add Isotope to element
   void AddIsotope(Isotope *Iso);
   /// Get Isotopes of Element
   std::vector<Isotope *> GetElementIsotopes();

   /// Get Element atomic number
   int GetZ();
   /// Get Element symbol
   std::string GetSymbol();

private:
   int AtomicNumber;
   double ElementMass = 0.;
   std::string Symbol;
   std::vector<Isotope *> Isotopes;

   friend class Stopping;
   void SetMass(double Mass);
   double GetMass();
};

#endif
