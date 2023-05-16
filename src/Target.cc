#include "Target.hh"

#include "TargetLayer.hh"

Target::Target() {
   TargetLayers.clear();
}

Target::~Target() {
}

void Target::AddLayer(double Thickness) {
   TargetLayers.push_back(new TargetLayer(Thickness));
}

void Target::InsertLayer(int PreviousLayer, double Thickness) {
   if (PreviousLayer++ <= TargetLayers.size())
      TargetLayers.insert(TargetLayers.begin() + PreviousLayer, new TargetLayer(Thickness));
}

void Target::RemoveLayer(int Layer) {
   if (Layer < TargetLayers.size())
      TargetLayers.erase(TargetLayers.begin() + Layer);
}

void Target::Clear() {
   for (int i = 0; i < TargetLayers.size(); ++i)
      RemoveLayer(i);
}

void Target::SetLayerThickness(double Thickness, int Layer) {
   if (Layer < TargetLayers.size())
      TargetLayers[Layer]->SetThickness(Thickness);
}

double Target::GetLayerThickness(int Layer) {
   if (Layer < TargetLayers.size())
      return TargetLayers[Layer]->GetThickness();
   return 0.;
}

void Target::AddElementToLayer(Element *el, int Layer) {
   if (TargetLayers.size() == 0)
      AddLayer();
   TargetLayers[Layer]->AddElement(el);
}

void Target::ChangeElementInLayer(int ElementIndex, Element *el, int Layer) {
   if (Layer < TargetLayers.size())
      TargetLayers[Layer]->SetElement(ElementIndex, el);
}

void Target::RemoveElementFromLayer(Element *el, int Layer) {
   if (Layer < TargetLayers.size())
      TargetLayers[Layer]->RemoveElement(el);
}

bool Target::ElementExistsInLayer(Element *el, int Layer) {
   if (Layer < TargetLayers.size())
      return (TargetLayers[Layer]->ElementInLayer(el) >= 0);
   else
      return false;
}

std::vector<Element *> Target::GetElementsInLayer(int Layer) {
   if (Layer < TargetLayers.size())
      return TargetLayers[Layer]->GetLayerElements();
   std::vector<Element *> elvec;
   return elvec;
}

void Target::SetElementAtomicPercentInLayer(Element *el, double AtomicPercent, int Layer) {
   if (Layer < TargetLayers.size())
      TargetLayers[Layer]->SetElementAtomicPercent(el, AtomicPercent);
}

double Target::GetElementAtomicPercentInLayer(Element *el, int Layer) {
   if (Layer < TargetLayers.size())
      return TargetLayers[Layer]->GetElementAtomicPercent(el);
   return 0.;
}

void Target::SetElementFitStateInLayer(Element *el, FitState FitState, int Layer) {
   if (Layer < TargetLayers.size())
      TargetLayers[Layer]->SetElementFitState(el, FitState);
}

FitState Target::GetElementFitStateInLayer(Element *el, int Layer) {
   if (Layer < TargetLayers.size())
      return TargetLayers[Layer]->GetElementFitState(el);
   return FitState::Free;
}

double Target::GetLayerAtomicPerCent(int Layer) {
   if (Layer < TargetLayers.size())
      return TargetLayers[Layer]->GetLayerAtomicPerCent();
   return 0.;
}

void Target::LinkElementsInLayer(Element *element, Element *LinkedElement, float linkfactor, int Layer) {
   if (Layer < TargetLayers.size())
      TargetLayers[Layer]->LinkElements(element, LinkedElement, linkfactor);
}

void Target::SetLinkFactorInLayer(Element *element, Element *LinkedElement, float linkfactor, int Layer) {
   LinkElementsInLayer(element, LinkedElement, linkfactor, Layer);
}

void Target::UnlinkElementsInLayer(Element *element, Element *LinkedElement, int Layer) {
   if (Layer < TargetLayers.size())
      TargetLayers[Layer]->UnlinkElements(element, LinkedElement);
}

std::map<Element *, std::map<Element *, double>> Target::GetLinkedElementsAndFactorsInLayer(int Layer) {
   if (Layer < TargetLayers.size())
      TargetLayers[Layer]->GetLinkedElementsAndFactors();
   std::map<Element *, std::map<Element *, double>> elmap;
   return elmap;
}
