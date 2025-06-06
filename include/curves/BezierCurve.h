//
// Created by fogoz on 10/05/2025.
//

#ifndef BEZIERCURVE_H
#define BEZIERCURVE_H
#include <memory>

#include "BaseCurve.h"


class BezierCurve : public BaseCurve{
    Position pos1;
    Position pos2;
    Position pos3;
    Position pos4;
    bool reversed;
public:
    BezierCurve(Position pos1, Position pos2, Position pos3, Position pos4, bool reversed=false);

    Position getPosition(double value, double h) override;

    Position getDerivative(double value) override;

    Position getSecondDerivative(double value);

    Position getThirdDerivative(double value);

#ifdef ENABLE_CURVATURE_POS
    static std::shared_ptr<BezierCurve> getBezierCurve(Position pos1, Position pos2, bool reversed=false);

    static std::shared_ptr<BezierCurve> getBezierCurveGridSearch(Position pos1, Position pos2, bool reversed=false);
#endif
};


#endif //BEZIERCURVE_H
