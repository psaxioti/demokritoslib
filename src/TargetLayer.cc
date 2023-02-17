#include "TargetLayer.hh"

#include "TargetElement.hh"

TargetLayer::TargetLayer(double thick) : Thickness(thick) {
}

TargetLayer::~TargetLayer() {
}

void TargetLayer::SetThickness(double thick) {
   Thickness = thick;
}

double TargetLayer::GetThickness() {
   return Thickness;
}

std::vector<Element *> TargetLayer::GetLayerElements() {
   std::vector<Element *> elems;
   for (std::size_t i = 0; i < TargetElements.size(); ++i)
      elems.push_back(TargetElements[i]->GetElement());
   return elems;
}

void TargetLayer::AddElement(Element *el, double atomic_percent) {
   if (ElementInLayer(el) < 0)
      TargetElements.push_back(new TargetElement(el, atomic_percent));
   CalculateAtomicPercent();
}

void TargetLayer::RemoveElement(Element *el) {
   int i = ElementInLayer(el);
   if (i >= 0)
      TargetElements.erase(TargetElements.begin() + i);
}

int TargetLayer::ElementInLayer(Element *el) {
   if (TargetElements.size() > 0)
      for (std::size_t i = 0; i < TargetElements.size(); ++i)
         if (TargetElements[i]->GetElement() == el)
            return i;
   return -1;
}

void TargetLayer::SetElementAtomicPercent(Element *el, double atomic_percent) {
   int i = ElementInLayer(el);
   if (i >= 0)
      TargetElements[i]->SetAtomicPercent(atomic_percent);
   CalculateAtomicPercent();
}

double TargetLayer::GetElementAtomicPercent(Element *el) {
   int i = ElementInLayer(el);
   if (i >= 0)
      return TargetElements[i]->GetAtomicPercent();
   return -1;
}

void TargetLayer::SetElementFitState(Element *el, FitState fitState) {
   int i = ElementInLayer(el);
   if (i >= 0)
      TargetElements[i]->SetFitState(fitState);
}

FitState TargetLayer::GetElementFitState(Element *el) {
   int i = ElementInLayer(el);
   if (i >= 0)
      return TargetElements[i]->GetFitState();
   return FitState::Free;
}

double TargetLayer::GetLayerAtomicPerCent() {
   return LayerAtomicPerCent;
}

void TargetLayer::CalculateAtomicPercent() {
   LayerAtomicPerCent = 0.;
   for (std::size_t i = 0; i < TargetElements.size(); ++i)
      LayerAtomicPerCent += TargetElements[i]->GetAtomicPercent();
}
