#ifndef TARGET_H
#define TARGET_H 1

#include <vector>

#include "FitStateEnum.hh"

class TargetLayer;
class TargetElement;
class Element;

class Target {

public:
   /**
    * Target class constructor.
    *
    * It will hold the information of the target (layers, elements, ...)
    */
   Target();
   ~Target();

   /// @brief Add a layer to the target
   /// @param Thickness Thicknes of the layer. Default is 100
   void AddLayer(double Thickness = 100);
   /// @brief Insert layer
   /// @param PreviousLayer Layer after which the new layer will be created
   /// @param Thickness Thickness of the layer
   void InsertLayer(int PreviousLayer, double Thickness);
   /// @brief Remove layer from target
   /// @param Layer Layer number to remove. If none given the first layer will be removed
   void RemoveLayer(int Layer = 0);

   /// @brief Set thickness of layer
   /// @param Thickness Thickness of the layer
   /// @param Layer Layer number for which the thickness is set
   void SetLayerThickness(double Thickness, int Layer = 0);
   /// @brief Get thickness of layer
   /// @param Layer Layer number for which the thickness is returned. Default is first layer of target.
   /// @return Thickness of layer
   double GetLayerThickness(int Layer = 0);

   /// @brief Adds an element to a target layer
   /// @param el Element to add
   /// @param Layer Layer at which the element will be added. Default is first layer of target
   void AddElementToLayer(Element *el, int Layer = 0);
   /// @brief Remove element from a layer
   /// @param el Element to remove
   /// @param Layer Layer from which the element will be removed. Default is first layer of target
   void RemoveElementFromLayer(Element *el, int Layer = 0);
   /// @brief Check if elemnt exists in a target layer
   /// @param el Element to check if already exists
   /// @param Layer Layer in which the element will be searched. Default is first layer of the target
   /// @return True if element already exists in layer, else false
   bool ElementExistsInLayer(Element *el, int Layer = 0);
   /// @brief Get elements in a target layer
   /// @param Layer Layer to get elements. Default is first layer of the target
   /// @return A vector with all the elements in the target layer
   std::vector<Element *> GetElementsInLayer(int Layer = 0);

   /// @brief Set atomic percent of element in layer
   /// @param el Element for which the atomic percentage is set
   /// @param AtomicPercent Atomic percentage to set (0 - 100)
   /// @param Layer Layer in which the atomic percentage of the element will be set. Default is first layer of target.
   void SetElementAtomicPercentInLayer(Element *el, double AtomicPercent, int Layer = 0);
   /// @brief Get atomic percent of element in layer
   /// @param el Element for which the atomic percentage will be returned
   /// @param Layer Layer in which the atomic percentage of the element will be returned. Default is first layer of target.
   /// @return The atomic percentage in the target layer for requested element (0 - 100)
   double GetElementAtomicPercentInLayer(Element *el, int Layer = 0);
   /// @brief Get the sum of atomic percentages in target layer
   /// @param Layer Layer in which the atomic percentage will be returned. Default is first layer of target
   /// @return Sum of the atomic percentages of all the elements in the target layer
   double GetLayerAtomicPerCent(int Layer = 0);

   void SetElementFitStateInLayer(Element *el, FitState FitState, int Layer = 0);
   FitState GetElementFitStateInLayer(Element *el, int Layer = 0);

private:
   std::vector<TargetLayer *> TargetLayers;
};

#endif