//
// Created by fogoz on 06/08/2025.
//

#ifndef CURVEVISUALISER_CURVEFACTORY_H
#define CURVEVISUALISER_CURVEFACTORY_H
#include <memory>

#include "BaseCurve.h"
#include <vector>


namespace CurveFactory {
    std::shared_ptr<BaseCurve> getBaseCurve(std::vector<double>& parameters);
}


#endif //CURVEVISUALISER_CURVEFACTORY_H