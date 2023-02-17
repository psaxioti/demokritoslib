#include "MassFunctions.hh"

#include <algorithm>
#include <istream>
#include <map>
#include <streambuf>

/// Vector containing the elements read
static std::vector<Element *> Elements;
/// Vector containing the isotopes read
static std::vector<Isotope *> Isotopes;

extern char masses_start[] asm("_binary_data_masses_dat_start");
extern char masses_end[] asm("_binary_data_masses_dat_end");

extern char abundance_start[] asm("_binary_data_abundance_dat_start");
extern char abundance_end[] asm("_binary_data_abundance_dat_end");

struct membuf : std::streambuf {
   membuf(char *begin, char *end) {
      this->setg(begin, begin, end);
   }
};

void ReadMassFiles() {
   if (Elements.size() > 0)
      return;

   membuf abuf(abundance_start, abundance_end);
   std::istream afile(&abuf);
   std::map<int, std::map<int, double>> abundances;
   while (afile.good()) {
      int z, n, a;
      double abund;
      while (afile.peek() == '#')
         afile.ignore(1000, '\n');
      std::string dum;
      afile >> z >> a >> dum >> abund >> dum;
      n = a - z;
      abundances[z][n] = abund;
   }

   membuf sbuf(masses_start, masses_end);
   std::istream file(&sbuf);

   bool ReadingIntro = false;
   while (file.good()) {
      char next_char = file.peek();
      if (next_char == '1' || ReadingIntro) {
         file.ignore(1000, '\n');
         if (next_char == '1' && !ReadingIntro)
            ReadingIntro = true;
         else if (next_char == '1' && ReadingIntro) {
            ReadingIntro = false;
            file.ignore(1000, '\n');
         }
      } else {
         std::string line;
         std::getline(file, line);
         if (std::count(line.begin(), line.end(), ' ') == line.length())
            continue;
         int n = stoi(line.substr(4, 5));
         int z = stoi(line.substr(9, 5));
         std::string sym = line.substr(20, 3);
         auto first = sym.find_first_not_of(" \t");
         auto last = sym.find_last_not_of(" \t") + 1;
         sym = sym.substr(first, last - first);
         std::string dum = line.substr(29, 14);
         dum = dum.substr(0, dum.find('#'));
         double massex = stod(dum);
         int mass_amu = stoi(line.substr(106, 3));
         dum = line.substr(110, 13);
         dum = dum.substr(0, dum.find('#'));
         double mass_amu_sub = stod(dum);
         double mass = (mass_amu_sub / 1000000.) + mass_amu;
         if (z == 0)
            continue;

         Element *readingel = nullptr;
         if (Elements.size() > 0)
            readingel = GetElement(sym);
         if (!readingel) {
            readingel = new Element(z, sym);
            Elements.push_back(readingel);
         }

         Isotope *readingiso = new Isotope(z, n, sym, mass, abundances[z][n], massex);
         Isotopes.push_back(readingiso);
         readingel->AddIsotope(readingiso);
      }
   }
   return;
}

std::vector<Element *> GetElements() {
   ReadMassFiles();
   return Elements;
}

std::vector<Isotope *> GetIsotopes() {
   ReadMassFiles();
   return Isotopes;
}

Element *GetElement(const std::string &Symbol) {
   ReadMassFiles();
   for (size_t i = 0; i < Elements.size(); ++i)
      if (Elements[i]->GetSymbol() == Symbol)
         return Elements[i];
   return nullptr;
}

Element *GetElement(const int &Z) {
   ReadMassFiles();
   for (size_t i = 0; i < Elements.size(); ++i)
      if (Elements[i]->GetZ() == Z)
         return Elements[i];
   return nullptr;
}

Isotope *GetIsotope(const std::string &Symbol) {
   ReadMassFiles();
   for (size_t i = 0; i < Isotopes.size(); ++i)
      if (Isotopes[i]->GetSymbol() == Symbol)
         return Isotopes[i];
   return nullptr;
}

Isotope *GetIsotope(const int &Z, const int &A) {
   ReadMassFiles();
   for (size_t i = 0; i < Isotopes.size(); ++i)
      if (Isotopes[i]->GetZ() == Z && Isotopes[i]->GetA() == A)
         return Isotopes[i];
   return nullptr;
}
