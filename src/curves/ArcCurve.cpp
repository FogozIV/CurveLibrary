//
// Created by fogoz on 10/05/2025.
//

#include "curves/ArcCurve.h"

#include <optional>
#include <optional>
#include <tuple>
#include <tuple>


ArcCurve::ArcCurve(Position center, double radius, Angle angleStart, Angle angleEnd): BaseCurve(0, 1), center(center), radius(radius), angleStart(angleStart), angleEnd(angleEnd) {
}

Position ArcCurve::getPosition(double value, double h) {
    double angle = angleStart.toRadians() + value * (angleEnd.toRadians() - angleStart.toRadians());
    double x = center.getX() + radius * cos(angle);
    double y = center.getY() + radius * sin(angle);
    double diff = (angleEnd - angleStart).warpAngle().toDegrees();
    bool isCCW = diff > 0;
    double orientation = angle + (isCCW ? M_PI/2 : -M_PI/2);
#ifdef ENABLE_CURVATURE_POS
    return Position(x, y, Angle::fromRadians(orientation), 1/radius);
#else
    return {x, y, Angle::fromRadians(orientation)};//Angle not used
#endif
}

Position ArcCurve::getDerivative(double value) {
    double start = angleStart.toRadians();
    double end = angleEnd.toRadians();
    double angle = start + value * (end - start);
    double dAngle = end - start;

    double dx = -radius * sin(angle) * dAngle;
    double dy =  radius * cos(angle) * dAngle;
    double heading = atan2(dy, dx);

    // Heading derivative θ' = curvature * speed
    double speed = std::sqrt(dx*dx + dy*dy); // ||v||
    double curvature = (radius == 0.0) ? 0.0 : 1.0 / radius;
    double headingDerivative = curvature * speed;

#ifdef ENABLE_CURVATURE_POS
    return Position(dx, dy, Angle::fromRadians(headingDerivative), 0.0);
#else
    return {dx, dy, Angle::fromRadians(headingDerivative)};
#endif
}

std::optional<ArcCurve> getArcCurve(Position begin, Position end, bool backward) {
    auto o_center = intersectPerpendicularLine(begin, end);
    if(!o_center.has_value()){
        return std::nullopt;
    }
    auto center = o_center.value();
    auto d0 = begin - center;
    auto d1 = end - center;

    auto angleStart = d0.getVectorAngle() ;
    auto angleEnd = d1.getVectorAngle() ;

    auto delta = angleEnd - angleStart;
    delta.warpAngle();

    if (backward) {
        if (delta.toDegrees() < 0) {
            delta += AngleConstants::FULL_TURN;
        }
    }else {
        if (delta.toDegrees() > 0) {
            delta -= AngleConstants::FULL_TURN;
        }
    }
    return ArcCurve(center, d0.getDistance(), angleStart, angleStart+delta);
}
