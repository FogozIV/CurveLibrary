//
// Created by fogoz on 10/05/2025.
//

#include "curves/BezierCurve.h"

#include <optional>
#include <optional>
#include <tuple>
#include <tuple>
#include <utility>

BezierCurve::BezierCurve(Position pos1, Position pos2, Position pos3, Position pos4, bool reversed): BaseCurve(0, 1.0), reversed(reversed) {
    this->pos1 = std::move(pos1);
    this->pos2 = std::move(pos2);
    this->pos3 = std::move(pos3);
    this->pos4 = std::move(pos4);
}

Position BezierCurve::getPosition(double value, double h) {
    if (reversed) {
        value = 1 - value;
    }
    Position result =  pos1 * pow(1-value, 3) + pos2 * (3*pow(1-value, 2)*value) + pos3 * (3*(1-value)*pow(value, 2)) + pos4 * pow(value, 3);
    Position derivative = getDerivative(value);
    Position secondDerivative = getSecondDerivative(value);
    return {result.getX(), result.getY(), Angle::fromRadians(atan2(derivative.getY(), derivative.getX())), (derivative.getX() * secondDerivative.getY() - derivative.getY() * secondDerivative.getX() ) / pow((pow(derivative.getX(), 2) + pow(derivative.getY(), 2)), 3.0/2.0)};

}

Position BezierCurve::getDerivative(double value) {
    if (reversed) {
        value = 1 - value;
    }
    // First derivative (velocity)
    Position derivative = (pos2 - pos1) * (3 * pow(1 - value, 2)) +
                          (pos3 - pos2) * (6 * (1 - value) * value) +
                          (pos4 - pos3) * (3 * pow(value, 2));

    // Second and third derivatives
    Position secondDerivative = getSecondDerivative(value);
    Position thirdDerivative = getThirdDerivative(value);

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

Position BezierCurve::getSecondDerivative(double value) {
    if (reversed) {
        value = 1 - value;
    }
    return 6 * (1-value) * (pos3 - 2*pos2 + pos1) + 6 * value * (pos4 - 2*pos3 + pos2);
}

Position BezierCurve::getThirdDerivative(double value) {
    return 6 * (pos4- 3* pos3 + 3*pos2 - pos1);
}

#ifdef ENABLE_CURVATURE_POS
std::shared_ptr<BezierCurve> BezierCurve::getBezierCurve(Position start, Position end, bool reversed) {
    int max_iters = 1000;
    Position pos1, pos2, pos3, pos4;
    pos1 = {start.getX(), start.getY()};
    pos4 = {end.getX(), end.getY()};
    double dx1 = cos(start.getAngle().toRadians());
    double dy1 = sin(start.getAngle().toRadians());
    double dx4 = cos(end.getAngle().toRadians());
    double dy4 = sin(end.getAngle().toRadians());

    double L1 = 10;
    double L4 = 10;
    for (int iter = 0; iter < max_iters; ++iter) {
        pos2 = {pos1.getX() + L1 * dx1, pos1.getY() + L1 * dy1};
        pos3 = {pos4.getX() - L4 * dx4, pos4.getY() - L4 * dy4};
        BezierCurve curve(pos1, pos2, pos3, pos4);
        double k1_est = curve.getPosition(0, 0.01).getCurvature();
        double k4_est = curve.getPosition(1, 0.01).getCurvature();
#ifndef ARDUINO
        std::cout << "k1_est: " << k1_est  << "k4_est: " << k4_est<< std::endl;
#endif
        double e_1 = k1_est - start.getCurvature();
        double e_4 = k4_est - end.getCurvature();
        bool curvature_ok = true;
        for (double t = 0.0f; t <= 1.0f; t += 0.05f) {
            double k = curve.getPosition(t, 0.01).getCurvature();
            if (fabs(k) > 0.004) {
                curvature_ok = false;
                break;
            }
        }

        if (!curvature_ok) {
            // Flatten the curve:
            L1 *= 1.05f;
            L4 *= 1.05f;
            continue;
        }
        if (fabs(e_1) < 1e-6 && fabs(e_4) < 1e-6) {
            break;
        }

        double delta = 0.01f;

        double L1p = L1 + delta;
        Position pos2p = {pos1.getX() + L1p * dx1, pos1.getY() + L1p * dy1};
        double k1p = BezierCurve(pos1, pos2p, pos3, pos4).getPosition(0, 0.01).getCurvature();
        double dk1dl = k1p/delta;

        double L4p = L4 - delta;
        Position pos3p = {pos4.getX() - L4p * dx4, pos4.getY() - L4p * dy4};
        double k4p = BezierCurve(pos1, pos2, pos3p, pos4).getPosition(1, 0.01).getCurvature();
        double dk4dl = k4p/delta;

        L1-= 0.00001 * e_1 * dk1dl;
        L4-= 0.00001 * e_4 * dk4dl;

        if (L1 <= 0) {
            L1 *= -1;
        }
        if (L4 <= 0) {
            L4 *= -1;
        }


    }
#ifndef ARDUINO
    std::cout << pos1 << " " << pos2 << " " << pos3 << " " << pos4 << std::endl;
#endif

    return std::make_shared<BezierCurve>(pos1, pos2, pos3, pos4, reversed);
}

std::shared_ptr<BezierCurve> BezierCurve::getBezierCurveGridSearch(Position start, Position end, bool reversed) {
    const double min_L = 20;
    const double max_L = 1000;
    const double step = 2;
    const int samples = 10;

    double best_err_L2 = min_L;
    double best_err_L3 = min_L;
    double best_err = 1e6;
    Position best_pos1, best_pos2, best_pos3, best_pos4;
    Position pos1, pos2, pos3, pos4;
    pos1 = {start.getX(), start.getY()};
    pos4 = {end.getX(), end.getY()};
    double dx1 = cos(start.getAngle().toRadians());
    double dy1 = sin(start.getAngle().toRadians());
    double dx4 = cos(end.getAngle().toRadians());
    double dy4 = sin(end.getAngle().toRadians());

    for (double L0 = min_L; L0 <= max_L; L0 += step) {
        for (double L3 = min_L; L3 <= max_L; L3 += step) {
            // Build control points
            pos2 = { pos1.getX() + L0 * dx1, pos1.getY() + L0 * dy1 };
            pos3 = { pos4.getX() - L3 * dx4, pos4.getY() - L3 * dy4 };
            BezierCurve curve(pos1, pos2, pos3, pos4);
            // Check curvature
            bool too_sharp = false;
            for (int i = 1; i < samples; ++i) {
                double t = i / (double)samples;
                double kappa = curve.getPosition(t, 0.01).getCurvature();
                if (kappa > 0.04) {
                    too_sharp = true;
                    break;
                }
            }
            if (too_sharp) continue;

            // End pose error (position + heading)
            Position end_pos = curve.getPosition(1, 0.01);
            double err = (end_pos - end).normCompleteRad(1, 100, 10000);

            if (err < best_err) {
                best_err = err;
                best_pos1 = curve.getPosition(0, 0.01);
                best_pos4 = curve.getPosition(1, 0.01);
                best_pos2 = pos2;
                best_pos3 = pos3;
                best_err_L2 = L0;
                best_err_L3 = L3;
                // Save best so far — P0..P3 already updated
            }
        }
    }
#ifndef ARDUINO
    std::cout << best_pos1 << " " << pos2 << " " << pos3 << " " << best_pos4 << std::endl;
#endif
    return std::make_shared<BezierCurve>(best_pos1, best_pos2, best_pos3, pos4, reversed);
}


#endif
