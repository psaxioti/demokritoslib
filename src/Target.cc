#include "Target.hh"

Target::Target() {
   TargetLayers.clear();
}

Target::~Target() {
}

void Target::AddLayer(double Thickness) {
   TargetLayers.push_back(new TargetLayer(Thickness));
}

void Target::InsertLayer(int PreviousLayer, double Thickness) {
   TargetLayers.insert(TargetLayers.begin() + PreviousLayer + 1, new TargetLayer(Thickness));
}

void Target::RemoveLayer(int Layer) {
   TargetLayers.erase(TargetLayers.begin() + Layer);
}

void Target::SetLayerThickness(double Thickness, int Layer) {
   TargetLayers[Layer]->SetThickness(Thickness);
}

double Target::GetLayerThickness(int Layer) {
   return TargetLayers[Layer]->GetThickness();
}

void Target::AddElementToLayer(Element *el, int Layer) {
   if (TargetLayers.size() == 0)
      AddLayer();
   TargetLayers[Layer]->AddElement(el);
}

void Target::RemoveElementFromLayer(Element *el, int Layer) {
   TargetLayers[Layer]->RemoveElement(el);
}

bool Target::ElementExistsInLayer(Element *el, int Layer) {
   if (Layer < TargetLayers.size())
      return (TargetLayers[Layer]->ElementInLayer(el) >= 0);
   else
      return false;
}

std::vector<Element *> Target::GetElementsInLayer(int Layer) {
   return TargetLayers[Layer]->GetLayerElements();
}

void Target::SetElementAtomicPercentInLayer(Element *el, double AtomicPercent, int Layer) {
   TargetLayers[Layer]->SetElementAtomicPercent(el, AtomicPercent);
}

double Target::GetElementAtomicPercentInLayer(Element *el, int Layer) {
   return TargetLayers[Layer]->GetElementAtomicPercent(el);
}

void Target::SetElementFitStateInLayer(Element *el, TargetElement::FitState FitState, int Layer) {
   TargetLayers[Layer]->SetElementFitState(el, FitState);
}

TargetElement::FitState Target::GetElementFitStateInLayer(Element *el, int Layer) {
   return TargetLayers[Layer]->GetElementFitState(el);
}

double Target::GetLayerAtomicPerCent(int Layer) {
   if (Layer < TargetLayers.size())
      return TargetLayers[Layer]->GetLayerAtomicPerCent();
   else
      return 0.;
}
