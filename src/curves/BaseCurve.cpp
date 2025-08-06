//
// Created by fogoz on 10/05/2025.
//

#include <curves/BaseCurve.h>
#ifndef ARDUINO
#include <iostream>
#endif


BaseCurve::BaseCurve(double minValue, double maxValue) : minValue(minValue), maxValue(maxValue){
    this->curveType = CurveFactory::BASE;
}

double BaseCurve::getMinValue() const {
    return minValue;
}

double BaseCurve::getMaxValue() const {
    return maxValue;
}

double BaseCurve::getLength(const double h) {
    return getLength(minValue, maxValue, h);
}

double BaseCurve::getLength(double ti, double t_end, double h) {
    if (std::isnan(ti)) {
        ti = minValue;
    }
    std::vector<double> y0 = {0.0};
    auto lambda = [this](double t, std::vector<double> &y, std::vector<double> &dydt) {
        dydt[0] = this->getDerivative(t).getDistance();
    };
    runge_kutta_4(ti, y0, t_end, h, lambda);
    return y0[0];
}

double BaseCurve::getValueForLength(double ti, double length, double h) {
    if (std::isnan(ti))
        ti = minValue;
    std::vector<double> y0 = {0.0};
    return runge_kutta_4_maximized(ti, y0, maxValue, h, [this](double t, std::vector<double> &y, std::vector<double> &dydt) {
        dydt[0] = this->getDerivative(t).getDistance();
    }, {length});
}

double BaseCurve::findNearest(Position pos, double h) {
    double dist = (getPosition(maxValue) - pos).getDistance();
    double index = maxValue;
    for (double i = minValue; i < maxValue; i += h) {
        double d = (getPosition(i) - pos).getDistance();
        if (d < dist) {
            dist = d;
            index = i;
        }
    }
    return index;
}

double BaseCurve::findNearest(Position pos, double guess, double h, int maxIter, double t_min, double t_max) {
    double bestT = guess;
    double bestDist = (getPosition(bestT) - pos).getDistance();

    for (int i = 0; i < maxIter; ++i) {
        double t_left = bestT - h;
        double t_right = bestT + h;

        // Clamp within bounds
        t_left = std::max(std::max(t_left, minValue), t_min);
        t_right = std::max(std::min(t_right, maxValue), t_max);

        double d_left = (getPosition(t_left) - pos).getDistance();
        double d_right = (getPosition(t_right) - pos).getDistance();

        if (d_left < bestDist) {
            bestT = t_left;
            bestDist = d_left;
        } else if (d_right < bestDist) {
            bestT = t_right;
            bestDist = d_right;
        } else {
            // Neither side is better: we are at local min
            break;
        }
    }

    return bestT;
}

bool BaseCurve::isBackward() {
    return backward;
}
