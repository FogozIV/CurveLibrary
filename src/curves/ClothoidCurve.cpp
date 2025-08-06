//
// Created by fogoz on 10/05/2025.
//

#include "curves/ClothoidCurve.h"

#include <tuple>

#include "utils/G2Solve3Arc.h"

ClothoidCurve::ClothoidCurve(Position start, double startCurvature, double curvatureRate,
                             double length) : BaseCurve(0, length) {
    this->curveType = CurveFactory::CLOTHOID;
    this->curve.m_CD.m_x0 = start.getX();
    this->curve.m_CD.m_y0 = start.getY();
    this->curve.m_CD.m_theta0 = start.getAngle().toRadians();
    this->curve.m_CD.m_kappa0 = startCurvature;
    this->curve.m_CD.m_dk = curvatureRate;
    this->curve.m_L = length;
}

ClothoidCurve::ClothoidCurve(ClothoidCurveV2 curve) : BaseCurve(0, curve.length()) {
    this->curveType = CurveFactory::CLOTHOID;
    this->curve = curve;
}

ClothoidCurve::ClothoidCurve(std::vector<double>& parameters) : BaseCurve(0, parameters[5]) {
    this->curveType = CurveFactory::CLOTHOID;
    this->curve.m_CD.m_x0 = parameters[0];
    this->curve.m_CD.m_y0 = parameters[1];
    this->curve.m_CD.m_theta0 = parameters[2];
    this->curve.m_CD.m_kappa0 = parameters[3];
    this->curve.m_CD.m_dk = parameters[4];
    parameters.erase(parameters.begin(), parameters.begin() + 6);
}

std::shared_ptr<ClothoidCurve> ClothoidCurve::getClothoidCurve(Position start, Angle endAngle, double initialCurvature,
                                                               double length) {
    Angle delta = endAngle - start.getAngle();
    double curvatureRate = delta.toRadians() - initialCurvature * length;
    curvatureRate /= pow(length, 2);
    curvatureRate *= 2;
    return std::make_shared<ClothoidCurve>(start, initialCurvature, curvatureRate, length);
}

std::shared_ptr<ClothoidCurve> ClothoidCurve::getClothoidCurveDelta(Position start, Angle deltaAngle,
                                                                    double initialCurvature, double length) {
    return getClothoidCurve(start, start.getAngle() + deltaAngle, initialCurvature, length);
}

std::shared_ptr<ClothoidCurve> ClothoidCurve::getClothoidCurveG0(Position start, double x, double y) {
    auto cc = ClothoidCurveV2();
    cc.build_G0Kappa(start.getX(), start.getY(), start.getAngle().toRadians(), start.getCurvature(), x, y);
    return std::make_shared<ClothoidCurve>(cc);
}

double ClothoidCurve::getCurvatureRateAngle(Angle delta, double length, double initialCurvature) {
    return (delta.toRadians() - initialCurvature * length) / pow(length, 2) * 2;
}


Position ClothoidCurve::getPosition(double value, double h) {
    double theta, kappa, x0, y0;
    curve.m_CD.evaluate(value, theta, kappa, x0, y0);
    return {x0, y0, Angle::fromRadians(theta), kappa};
}

Position ClothoidCurve::getDerivative(double value) {
    double angle = curve.m_CD.m_theta0 + start_curvature() * value + 0.5 * curvature_rate() * value * value;
    double current_curvature = start_curvature() + curvature_rate() * value;
    return {cos(angle), sin(angle), Angle::fromRadians(current_curvature), curvature_rate()};
}

double ClothoidCurve::getLength(double ti, double t_end, double h) {
    return t_end - ti;
}

double ClothoidCurve::getValueForLength(double ti, double length, double h) {
    return std::min(ti + length, this->curve.m_L);
}

std::vector<double> ClothoidCurve::getParameters() {
    return {curve.m_CD.m_x0, curve.m_CD.m_y0, curve.m_CD.m_theta0, curve.m_CD.m_kappa0, curve.m_CD.m_dk, curve.m_L};
}
