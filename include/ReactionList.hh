#ifndef REACTIONLIST_H
#define REACTIONLIST_H 1

#include "Reaction.hh"

#include <vector>

class ReactionList {

public:
   /// @brief ReactionList class constructor.
   ///   It will read or create the .reactions.dat file in R33Folder and hold the available r33 files
   /// @param R33Folder Folder in which to look for r33 files.
   ///   If none provided it will look in local R33Data or ${HOME}/.R33Data folders
   ReactionList(std::string R33Folder = "");

   /// Returns the folder were the r33 files are located.
   std::string GetR33Folder();
   /// Sets the folder with the r33 files.
   void SetR33Folder(std::string NewR33Folder);
   /// @brief Force recreation of reactions list DB file
   void RecreateList();

   /// @brief Get the complete reaction list of files in set folder
   /// @return Vector of reactions
   const std::vector<Reaction *> GetReactionList();
   /// @brief Get a subset of reactions for specific beam isotope and target isotope
   /// @param BeamIsotope isotope of the beam
   /// @param TargetIsotope isotope of target
   /// @return Vector of reactions for the specified reaction
   const std::vector<Reaction *> GetFilteredReactionList(Isotope *BeamIsotope, Isotope *TargetIsotope);
   /// @brief Get a subset of reactions for specific beam isotope and target element
   /// @param BeamIsotope isotope of the beam
   /// @param TargetElement element of target
   /// @return Vector of reactions for the specified reaction
   const std::vector<Reaction *> GetFilteredReactionList(Isotope *BeamIsotope, Element *TargetElement);

private:
   std::string R33Folder;
   std::vector<Reaction *> Reactions;

   void CreateReactionList();
   void ReadReactionList();
};

#endif