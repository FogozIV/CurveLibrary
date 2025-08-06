//
// Created by fogoz on 02/08/2025.
//

#ifndef HERMITESPLINECURVE_H
#define HERMITESPLINECURVE_H
#include <memory>

#include "BaseCurve.h"

class HermiteCondition {
    double L;

public:
    const Position p, dp, ddp;

    HermiteCondition(Position pos, double L) : L(L), p(pos), dp(L*p.getSinCosAngle()), ddp(p.getCurvature() * pow(L, 2) * p.getNormalVector()) {
    }
};


template<size_t Order>
class HermiteSplineCurve : public BaseCurve {
public:
    std::array<Position, Order + 1> pos_coeff;

    HermiteSplineCurve(std::array<Position, Order + 1> pos_coeff) : BaseCurve(0, 1), pos_coeff(pos_coeff) {
    }

    static std::shared_ptr<HermiteSplineCurve<Order>> getSplineFromStartEnd(Position start, Position end, double L) {
        assert(Order >= 5);
        auto s = HermiteCondition(start, L);
        auto e = HermiteCondition(end, L);
        std::array<Position, Order + 1> pos_coeff;
        pos_coeff[0] = s.p;
        pos_coeff[1] = s.dp;
        pos_coeff[2] = s.ddp / 2;
        pos_coeff[3] = 10 * (end - start) - 6 * s.dp - 4 * e.dp - 1.5 * s.ddp + 0.5 * e.ddp;
        pos_coeff[4] = -15 * (end - start) + 8 * s.dp + 7 * e.dp + 1.5 * s.ddp - e.ddp;
        pos_coeff[5] = 6 * (end - start) - 3 * s.dp - 3 * e.dp - 0.5 * s.ddp + 0.5 * e.ddp;
        return std::make_shared<HermiteSplineCurve<5> >(pos_coeff);
    }

    template<size_t derivative_count>
    std::array<Position, derivative_count + 1> eval(double u) const {
        std::array<Position, derivative_count + 1> D = {};

        D[0] = pos_coeff[Order]; // start with top coefficient

        for (int i = Order - 1; i >= 0; --i) {
            for (int k = derivative_count; k > 0; --k) {
                D[k] = D[k] * u + k * D[k - 1];
            }
            D[0] = D[0] * u + pos_coeff[i];
        }
        return D;
    }

    Position getPosition(double value, double h) override {
        auto data = eval<2>(value);
        const Position& derivative = data[1];
        const Position& secondDerivative = data[2];
        return {data[0].getX(), data[0].getY(), Angle::fromRadians(atan2(derivative.getY(), derivative.getX())), (derivative.getX() * secondDerivative.getY() - derivative.getY() * secondDerivative.getX() ) / pow((pow(derivative.getX(), 2) + pow(derivative.getY(), 2)), 3.0/2.0)};
    }

    Position getDerivative(double value) override {
        auto data = eval<3>(value);
        const Position& derivative = data[1];
        const Position& secondDerivative = data[2];
        const Position& thirdDerivative = data[3];

        double dx = derivative.getX();
        double dy = derivative.getY();
        double ddx = secondDerivative.getX();
        double ddy = secondDerivative.getY();
        double dddx = thirdDerivative.getX();
        double dddy = thirdDerivative.getY();

        // θ' = (x' y'' - y' x'') / (x'^2 + y'^2)
        double headingDerivative = (dx * ddy - dy * ddx) / (dx * dx + dy * dy);
        double curvature_derivative = ((ddx * dddy - ddy * dddx) * (dx * dx + dy * dy) -
                                       3 * (dx * ddy - dy * ddx) * (dx * ddx + dy * ddy)) /
                                      pow(dx * dx + dy * dy, 3);
        return {dx, dy, Angle::fromRadians(headingDerivative), curvature_derivative};
    }
};

class QuaticHermiteSpline : public HermiteSplineCurve<5> {
public:
    explicit QuaticHermiteSpline(const std::array<Position, 6> &pos_coeff)
        : HermiteSplineCurve<5>(pos_coeff) {
    }
};

#endif //HERMITESPLINECURVE_H
