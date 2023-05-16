#ifndef FITSTATE_H
#define FITSTATE_H 1

#include <string>

enum class FitState { Free,
                      Fixed,
                      Fit };

static std::string GetFitStateString(FitState state) {
   switch (state) {
      case FitState::Free:
         return "Free";
      case FitState::Fixed:
         return "Fixed";
      case FitState::Fit:
         return "Fit";
      default:
         return "Free";
   }
};

static FitState GetFitStateFromString(std::string state) {
   if (state == "Free")
      return FitState::Fixed;
   else if (state == "Fixed")
      return FitState::Fixed;
   else if (state == "Fit")
      return FitState::Fit;
   return FitState::Fixed;
};

#endif