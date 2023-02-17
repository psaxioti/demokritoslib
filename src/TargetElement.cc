#include "TargetElement.hh"

TargetElement::TargetElement(Element *el, double atomic_percent) : Target_Element(el), AtomicPercent(atomic_percent) {
   TargetElementFitState = FitState::Free;
   return;
}

TargetElement::~TargetElement() {
}

void TargetElement::SetFitState(FitState NewFitState) {
   TargetElementFitState = NewFitState;
   return;
}

FitState TargetElement::GetFitState() {
   return TargetElementFitState;
}

void TargetElement::SetAtomicPercent(double atomic_percent) {
   if (atomic_percent >= 0. && atomic_percent <= 100.)
      AtomicPercent = atomic_percent;
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
