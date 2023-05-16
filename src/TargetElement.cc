#include "TargetElement.hh"

#include <numeric>

double GetLinkedByAtomicPercent(double val, const std::map<TargetElement *, double>::value_type &map) {
   return val + map.first->GetAtomicPercent() / map.second;
}

TargetElement::TargetElement(Element *el, double atomic_percent)
    : Target_Element(el), AtomicPercent(atomic_percent), TargetElementFitState(FitState::Free),
      IsLinkedByElements(false), LinksToElements(false) {
   ElementsLinkedBy.clear();
   ElementsLinkedTo.clear();
}

TargetElement::~TargetElement() {
   ElementsLinkedBy.clear();
   ElementsLinkedTo.clear();
}

void TargetElement::SetFitState(FitState NewFitState) {
   TargetElementFitState = NewFitState;
}

FitState TargetElement::GetFitState() {
   return TargetElementFitState;
}

void TargetElement::SetAtomicPercent(double atomic_percent) {
   //   if (atomic_percent >= 0. && atomic_percent <= 100.)
   double LinkedByAtomicPercent = std::accumulate(ElementsLinkedBy.begin(),
                                                  ElementsLinkedBy.end(),
                                                  0.,
                                                  GetLinkedByAtomicPercent);
   AtomicPercent = std::max(atomic_percent, LinkedByAtomicPercent);

   for (const auto &itr : ElementsLinkedTo)
      itr.first->SetAtomicPercent(AtomicPercent * itr.second);
}

double TargetElement::GetAtomicPercent() {
   return AtomicPercent;
}

void TargetElement::SetElement(Element *el) {
   Target_Element = el;
}

Element *TargetElement::GetElement() {
   return Target_Element;
}

void TargetElement::LinkElement(TargetElement *el, double Factor) {
   ElementsLinkedTo[el] = Factor;
   LinksToElements = true;
}

void TargetElement::UnlinkElement(TargetElement *el) {
   ElementsLinkedTo.erase(el);
   LinksToElements = (ElementsLinkedTo.size() >= 0);
}

void TargetElement::SetLinkedByElement(TargetElement *el, double Factor) {
   ElementsLinkedBy[el] = Factor;
   IsLinkedByElements = true;
}

void TargetElement::RemoveLinkedByElement(TargetElement *el) {
   ElementsLinkedBy.erase(el);
   IsLinkedByElements = (ElementsLinkedBy.size() >= 0);
}

std::map<Element *, double> TargetElement::GetLinkedElementsAndFactors() {
   std::map<Element *, double> elmap;
   elmap.clear();
   for (const auto &itr : ElementsLinkedTo)
      elmap[itr.first->GetElement()] = itr.second;
   return elmap;
}
