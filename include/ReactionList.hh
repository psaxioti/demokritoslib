#ifndef REACTIONLIST_H
#define REACTIONLIST_H 1

#include "Reaction.hh"

#include <vector>

class ReactionList {

public:
   ReactionList(std::string R33Folder = "");

   std::string GetR33Folder();
   void SetR33Folder(std::string NewR33Folder);

private:
   std::string R33Folder;
   std::vector<Reaction *> Reactions;

   void CreateReactionList();
   void ReadReactionList();
};

#endif