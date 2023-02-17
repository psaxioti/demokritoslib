#ifndef TARGET_ELEMENT_H
#define TARGET_ELEMENT_H 1

#include "FitStateEnum.hh"

class Element;

class TargetElement {
   friend class TargetLayer;

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