#ifndef TARGET_LAYER_H
#define TARGET_LAYER_H 1

#include <map>
#include <vector>

#include "FitStateEnum.hh"

class TargetElement;
class Element;

class TargetLayer {
   friend class Target;

private:
   TargetLayer(double thick = 100.);
   ~TargetLayer();

   void SetThickness(double thick);
   double GetThickness();

   std::vector<Element *> GetLayerElements();
   void AddElement(Element *el, double atomic_percent = 100.);
   void SetElement(int Index, Element *el);
   void RemoveElement(Element *el);

   void SetElementAtomicPercent(Element *el, double atomic_percent);
   double GetElementAtomicPercent(Element *el);

   void SetElementFitState(Element *el, FitState fitState);
   FitState GetElementFitState(Element *el);

   double GetLayerAtomicPerCent();

   void LinkElements(Element *element, Element *LinkedElement, float linkfactor);
   void UnlinkElements(Element *element, Element *LinkedElement);
   std::map<Element *, std::map<Element *, double>> GetLinkedElementsAndFactors();

private:
   int ElementInLayer(Element *el);
   void CalculateAtomicPercent();

   std::vector<TargetElement *> TargetElements;
   double Thickness;
   double LayerAtomicPerCent;
};

#endif