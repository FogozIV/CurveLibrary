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

class ClothoidCurve : public BaseCurve {
    ClothoidCurveV2 curve;
public:
    ClothoidCurve(Position start, double startCurvature, double curvatureRate, double length);

    ClothoidCurve(ClothoidCurveV2 curve);

    static std::shared_ptr<ClothoidCurve> getClothoidCurve(Position start, Angle endAngle, double initialCurvature, double length);

    static std::shared_ptr<ClothoidCurve> getClothoidCurveDelta(Position start, Angle deltaAngle, double initialCurvature,
                                                                double length);

    static std::optional<std::shared_ptr<ClothoidCurve>> getClothoidCurve(Position start, Position end);

    static double getCurvatureRateAngle(Angle delta, double length, double initialCurvature);

    Position getPosition(double value, double h) override;

    Position getDerivative(double value) override;

    double getLength(double ti, double t_end, double h) override;

    double getValueForLength(double ti, double length, double h) override;

    [[nodiscard]] Position start1() const {
        return {curve.m_CD.m_x0, curve.m_CD.m_y0, Angle::fromRadians(curve.m_CD.m_theta0), curve.m_CD.m_kappa0};
    }

    [[nodiscard]] double start_curvature() const {
        return curve.m_CD.m_kappa0;
    }

    [[nodiscard]] double curvature_rate() const {
        return curve.m_CD.m_dk;
    }
};


#endif //CLOTHOID_H
