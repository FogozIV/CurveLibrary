//
// Created by fogoz on 10/05/2025.
//

#ifndef CLOTHOID_H
#define CLOTHOID_H
#include <memory>
#include <optional>
#include <tuple>

#include "BaseCurve.h"


class ClothoidCurve : public BaseCurve {
    Position start;
    double startCurvature;
    double curvatureRate;
    double length;
    std::optional<std::tuple<double, Position, double>> previousPos= std::nullopt;
    void clothoidODE(double t, std::vector<double> &y, std::vector<double> &dydt) const;
public:
    ClothoidCurve(Position start, double startCurvature, double curvatureRate, double length);

    ClothoidCurve(Position start, double initialCurvature, std::array<double, 2> cl);

    ClothoidCurve(Position start, std::array<double, 3> kcl);

    static std::shared_ptr<ClothoidCurve> getClothoidCurve(Position start, Angle endAngle, double initialCurvature, double length);

    static std::shared_ptr<ClothoidCurve> getClothoidCurveDelta(Position start, Angle deltaAngle, double initialCurvature,
                                                                double length);
    static std::optional<std::shared_ptr<ClothoidCurve>> getClothoidCurve(Position start, Position end, double initialCurvature);

    static double getCurvatureRateAngle(Angle delta, double length, double initialCurvature);

    Position localToGlobal(const Position &local) const;

    Position getPosition(double value, double h) override;


    Position getDerivative(double value) override;

    double getLength(double ti, double t_end, double h) override;

    double getValueForLength(double ti, double length, double h) override;

};
std::optional<std::array<double, 3>> buildClothoid(Position begin, Position end);

std::optional<std::array<double, 2>> buildClothoid(Position begin, Position end, double kappa0);


#endif //CLOTHOID_H
