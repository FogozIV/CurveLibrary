//
// Created by fogoz on 10/05/2025.
//

#ifndef CLOTHOID_H
#define CLOTHOID_H
#include <memory>
#include <optional>
#include <tuple>

#include "BaseCurve.h"
#include "utils/G2Solve3Arc.h"
#include "utils/Matrix.h"

struct ClothoidSegmentOld {
    Matrix<2,1> xy;
    double length;
    double theta;
    double kappa;
    double dkappa;
};

class ClothoidCurve : public BaseCurve {
    Position start;
    double startCurvature;
    double curvatureRate;
    double length;
    std::optional<std::tuple<double, Position, double>> previousPos= std::nullopt;
    void clothoidODE(double t, std::vector<double> &y, std::vector<double> &dydt) const;
    ClothoidCurveV2 curve;
public:
    ClothoidCurve(Position start, double startCurvature, double curvatureRate, double length);

    ClothoidCurve(Position start, double initialCurvature, std::array<double, 2> cl);

    ClothoidCurve(Position start, std::array<double, 3> kcl);

    ClothoidCurve(const ClothoidSegmentOld& segment);

    ClothoidCurve(const ClothoidSegment& segment);

    ClothoidCurve(ClothoidCurveV2 curve);

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

    [[nodiscard]] Position start1() const {
        return start;
    }

    [[nodiscard]] double start_curvature() const {
        return startCurvature;
    }

    [[nodiscard]] double curvature_rate() const {
        return curvatureRate;
    }

    [[nodiscard]] std::optional<std::tuple<double, Position, double>> previous_pos() const {
        return previousPos;
    }
};

double integralXk(double a, double b, double c, double k=0, double h=0.0001);

double integralYk(double a, double b, double c, double k=0, double h=0.0001);
std::optional<std::array<double, 3>> buildClothoid(Position begin, Position end);


#endif //CLOTHOID_H
