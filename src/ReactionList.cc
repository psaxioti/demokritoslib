#include "ReactionList.hh"

#include <algorithm>
#include <filesystem>
#include <fstream>

#include <limits.h>
#include <unistd.h>

#include "MassFunctions.hh"
#include "Reaction.hh"

ReactionList::ReactionList(std::string R33folder) {
   if (R33folder.empty()) {
      char result[PATH_MAX];
      ssize_t count = readlink("/proc/self/exe", result, PATH_MAX);
      R33Folder = std::string(result, (count > 0) ? count : 0);
      if (R33Folder == "/usr/local/bin")
         R33Folder = std::string(std::getenv("HOME")) + "/.R33Data";
      else
         R33Folder += "/R33Data";
   } else
      R33Folder = R33folder;
   std::string R33ReactionsFile = R33Folder + "/.reactions.dat";
   if (std::filesystem::exists(std::filesystem::path(R33ReactionsFile)))
      ReadReactionList();
   else
      CreateReactionList();
}

void ReactionList::ReadReactionList() {
   int No;
   float E, th, Emin, Emax;
   std::string BeamSymbol, TargetSymbol;
   std::string source,
       name, file_name;
   std::ifstream ReactionFile;

   Reactions.clear();

   ReactionFile.open(R33Folder + "/.reactions.dat");
   while (ReactionFile >> BeamSymbol >> TargetSymbol >> E >> th >> No >> Emin >> Emax >> name >> file_name) {
      std::getline(ReactionFile, source, '\n');

      Reactions.push_back(new Reaction(file_name));
   }
   ReactionFile.close();
   return;
}

void ReactionList::CreateReactionList() {
   std::ofstream ReactionFile;
   ReactionFile.open(R33Folder + "/.reactions.dat");

   Reactions.clear();

   for (const auto &Entry : std::filesystem::directory_iterator(R33Folder)) {
      if (Entry.is_regular_file()) {
         std::filesystem::path File = Entry.path();
         if (File.extension() == ".R33" || File.extension() == ".r33") {
            Reactions.push_back(new Reaction(File));
            ReactionFile << Reactions[Reactions.size() - 1]->GetBeamIsotope()->GetSymbol() << " "
                         << Reactions[Reactions.size() - 1]->GetTargetIsotope()->GetSymbol() << " "
                         << Reactions[Reactions.size() - 1]->GetEgamma() << " "
                         << Reactions[Reactions.size() - 1]->GetTheta() << " "
                         << Reactions[Reactions.size() - 1]->GetNumberOfPoints() << " "
                         << Reactions[Reactions.size() - 1]->GetEmin() << " "
                         << Reactions[Reactions.size() - 1]->GetEmax() << " "
                         << Reactions[Reactions.size() - 1]->GetReactionName() << " "
                         << Reactions[Reactions.size() - 1]->GetFileName() << " "
                         << Reactions[Reactions.size() - 1]->GetSource() << std::endl;
         }
      }
   }
   ReactionFile.close();
}

const std::vector<Reaction *> ReactionList::GetReactionList() {
   return Reactions;
}

const std::vector<Reaction *> ReactionList::GetFilteredReactionList(Isotope *BeamIsotope, Isotope *TargetIsotope) {
   std::vector<Reaction *> FilteredList;
   FilteredList.clear();
   std::copy_if(Reactions.begin(), Reactions.end(), std::back_inserter(FilteredList),
                [BeamIsotope, TargetIsotope](Reaction *reaction) { return (reaction->GetBeamIsotope() == BeamIsotope && reaction->GetTargetIsotope() == TargetIsotope); });
   //   for (const auto &reaction : Reactions)
   //      if (reaction->GetAtomicNumber() == TargetZ && reaction->GetMassNumber() == TargetA)
   //         FilteredList.push_back(reaction);
   return FilteredList;
}

const std::vector<Reaction *> ReactionList::GetFilteredReactionList(Isotope *BeamIsotope, Element *TargetElement) {
   std::vector<Reaction *> FilteredList;
   FilteredList.clear();
   std::copy_if(Reactions.begin(), Reactions.end(), std::back_inserter(FilteredList),
                [BeamIsotope, TargetElement](Reaction *reaction) { return (reaction->GetBeamIsotope() == BeamIsotope && reaction->GetTargetElement() == TargetElement); });
   //   for (const auto &reaction : Reactions)
   //      if (reaction->GetAtomicNumber() == TargetZ && reaction->GetMassNumber() == TargetA)
   //         FilteredList.push_back(reaction);
   return FilteredList;
}

std::string ReactionList::GetR33Folder() {
   return R33Folder;
}

void ReactionList::SetR33Folder(std::string NewR33Folder) {
   R33Folder = NewR33Folder;
   std::string R33ReactionsFile = R33Folder + "/.reactions.dat";
   if (std::filesystem::exists(std::filesystem::path(R33ReactionsFile)))
      ReadReactionList();
   else
      CreateReactionList();
}

void ReactionList::RecreateList() {
   CreateReactionList();
}
