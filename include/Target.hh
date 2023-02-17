#ifndef TARGET_H
#define TARGET_H 1

#include "TargetLayer.hh"

class Target {

public:
   Target();
   ~Target();

   void AddLayer(double Thickness = 100);
   void InsertLayer(int PreviousLayer, double Thickness);
   void RemoveLayer(int Layer = 0);

   void SetLayerThickness(double Thickness, int Layer = 0);
   double GetLayerThickness(int Layer = 0);

   void AddElementToLayer(Element *el, int Layer = 0);
   void RemoveElementFromLayer(Element *el, int Layer = 0);
   bool ElementExistsInLayer(Element *el, int Layer = 0);
   std::vector<Element *> GetElementsInLayer(int Layer = 0);

   void SetElementAtomicPercentInLayer(Element *el, double AtomicPercent, int Layer = 0);
   double GetElementAtomicPercentInLayer(Element *el, int Layer = 0);
   double GetLayerAtomicPerCent(int Layer = 0);

   void SetElementFitStateInLayer(Element *el, TargetElement::FitState FitState, int Layer = 0);
   TargetElement::FitState GetElementFitStateInLayer(Element *el, int Layer = 0);

private:
   std::vector<TargetLayer *> TargetLayers;
};

#endif