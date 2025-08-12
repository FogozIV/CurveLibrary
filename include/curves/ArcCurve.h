//
// Created by fogoz on 10/05/2025.
//

#ifndef ARCCURVE_H
#define ARCCURVE_H

#include <optional>
#include <optional>
#include <tuple>
#include <tuple>

#include "BaseCurve.h"


class ArcCurve : public BaseCurve{
    Position center;
    double radius;
    Angle angleStart;
    Angle angleEnd;
public:
    ArcCurve(Position center, double radius, Angle angleStart, Angle angleEnd);

    ArcCurve(Position begin, double radius, Angle angleEnd);

    ArcCurve(std::vector<double>& parameters);

    Position getPosition(double value, double h) override;

    Position getDerivative(double value) override;

    std::vector<double> getParameters() override;
};

std::optional<ArcCurve> getArcCurve(Position begin, Position end);

#endif //ARCCURVE_H
