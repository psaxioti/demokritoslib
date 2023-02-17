#ifndef ISOTOPE_H
#define ISOTOPE_H 1

#include <string>

class Isotope {
public:
   /**
    * Isoope class constructor.
    *
    * @param Z Atomic number of isotope.
    * @param N Neutron numbers of isotope.
    * @param Symb Symbol of the isotope in the form of eg. 16O.
    * @param mass Mass of isotope in amu.
    * @param Abund Natural abundance of isotope in per cent (0.-100.).
    *    [optional] If not provided it will be set to 0.
    * @param massexc Mass excess of isotope in keV.
    *    [optional] If not provided it will be calculated from the mass provided.
    */
   Isotope(int Z, int N, std::string Symb, double mass, double Abund = 0., double massexc = 0.);
   ~Isotope();
   /// Get Isotope per cent abundance set (0 - 100) if different than the natural
   double GetAbundance();
   /// Get Isotope per cent natural abundance (0 - 100)
   double GetNaturalAbundance();
   /// Set Isotope per cent abundance (0 - 100) if different than the natural
   void SetAbundance(double NewAbund);
   /// Get Isotope mass in amu
   double GetMass();
   /// Get Isotope mass excess in keV
   double GetMassExcess();
   /// Get Isotope symbol (eg. 16O)
   std::string GetSymbol();
   /// Get Isotope atomic number
   int GetZ();
   /// Get Isotope neutron number
   int GetN();
   /// Get Isotope mass number
   int GetA();

private:
   double Mass, MassExcess;
   double Abundance, NaturalAbundance;
   int AtomicNumber, NeutronNumber;
   std::string Symbol;
};

#endif
