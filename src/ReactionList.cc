#include "ReactionList.hh"

#include <filesystem>
#include <fstream>

#include <limits.h>
#include <unistd.h>

#include "Reaction.hh"

ReactionList::ReactionList(std::string R33folder) {
   if (R33folder.empty()) {
      char result[PATH_MAX];
      ssize_t count = readlink("/proc/self/exe", result, PATH_MAX);
      R33Folder = std::string(result, (count > 0) ? count : 0);
      if (R33Folder == "/src/local/bin")
         R33Folder == std::string(std::getenv("HOME")) + "/.R33Data";
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
   int Z, A, No;
   float E, th, Emin, Emax;
   std::string source, name, file_name;
   std::ifstream ReactionFile;

   ReactionFile.open(R33Folder + "/.reactions.dat");
   while (ReactionFile >> Z >> A >> E >> th >> No >> Emin >> Emax >> name >> file_name) {
      std::getline(ReactionFile, source, '\n');

      Reactions.push_back(new Reaction(file_name));
   }
   ReactionFile.close();
   return;
}

void ReactionList::CreateReactionList() {
   std::ofstream ReactionFile;
   ReactionFile.open(R33Folder + "/.reactions.dat");

   for (const auto &Entry : std::filesystem::directory_iterator(R33Folder)) {
      if (Entry.is_regular_file()) {
         std::filesystem::path File = Entry.path();
         if (File.extension() == "R33" || File.extension() == "r33") {
            Reactions.push_back(new Reaction(File));
            ReactionFile << Reactions[Reactions.size() - 1]->GetAtomicNumber() << " "
                         << Reactions[Reactions.size() - 1]->GetMassNumber() << " "
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