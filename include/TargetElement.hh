#ifndef TARGET_ELEMENT_H
#define TARGET_ELEMENT_H 1

#include "FitStateEnum.hh"

#include <map>

class Element;

class TargetElement {
   friend class TargetLayer;
   friend double GetLinkedByAtomicPercent(double val, const std::map<TargetElement *, double>::value_type &map);

private:
   TargetElement(Element *el = nullptr, double atomic_percent = 100.);
   ~TargetElement();

   void SetFitState(FitState NewFitState);
   FitState GetFitState();

   void SetAtomicPercent(double atomic_percent);
   double GetAtomicPercent();

   void SetElement(Element *el);
   Element *GetElement();

   void LinkElement(TargetElement *, double Factor = 1.);
   void UnlinkElement(TargetElement *);
   void SetLinkedByElement(TargetElement *, double Factor = 1.);
   void RemoveLinkedByElement(TargetElement *);
   std::map<Element *, double> GetLinkedElementsAndFactors();

private:
   Element *Target_Element;
   double AtomicPercent;
   FitState TargetElementFitState;

   bool IsLinkedByElements;
   bool LinksToElements;
   std::map<TargetElement *, double> ElementsLinkedBy;
   std::map<TargetElement *, double> ElementsLinkedTo;
};

#endif