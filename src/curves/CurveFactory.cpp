//
// Created by fogoz on 06/08/2025.
//

#include <curves/CurveFactory.h>

#include "curves/ArcCurve.h"
#include "curves/BezierCurve.h"
#include "curves/ClothoidCurve.h"
#include "curves/CurveList.h"
#include "curves/HermiteSplineCurve.h"

std::shared_ptr<BaseCurve> CurveFactory::getBaseCurve(std::vector<double>& parameters) {
    auto type = static_cast<CurveFactory::CurveTypes>(parameters[0]);
    parameters.erase(parameters.begin());
    switch (type) {
        case BASE:
            return nullptr;
            break;
        case ARC:
            return std::make_shared<ArcCurve>(parameters);
            break;
        case BEZIER:
            return std::make_shared<BezierCurve>(parameters);
            break;
        case CLOTHOID:
            return std::make_shared<ClothoidCurve>(parameters);
            break;
        case LIST:
            return std::make_shared<CurveList>(parameters);
            break;
        case SPLINE: {
            int size = parameters[0];
            parameters.erase(parameters.begin());
            if (size == 5) {
                return std::make_shared<QuinticHermiteSpline>(parameters);
            }else {
                return nullptr;
            }
        }
    }
}
