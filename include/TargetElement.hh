#ifndef TARGET_ELEMENT_H
#define TARGET_ELEMENT_H 1

#include "Element.hh"

class TargetElement {
   friend class TargetLayer;

public:
   enum class FitState { Free,
                         Fixed,
                         Fit };

private:
   TargetElement(Element *el = nullptr, double atomic_percent = 100.);
   ~TargetElement();

   void SetFitState(FitState NewFitState);
   FitState GetFitState();

   void SetAtomicPercent(double atomic_percent);
   double GetAtomicPercent();

   void SetElement(Element *el);
   Element *GetElement();

private:
   Element *Target_Element;
   double AtomicPercent;
   FitState TargetElementFitState;
};

#endif