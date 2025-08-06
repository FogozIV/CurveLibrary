//
// Created by fogoz on 10/05/2025.
//

#ifndef BASECURVE_H
#define BASECURVE_H
#include <cstdint>

#include "utils/Integrator.h"
#include "utils/Position.h"

namespace CurveFactory {
    enum CurveTypes {
        BASE, ARC, BEZIER, CLOTHOID, LIST, SPLINE
    };
}


class BaseCurve {
protected:
    double minValue;
    double maxValue;
    bool backward = false;
    CurveFactory::CurveTypes curveType;

public:
    virtual ~BaseCurve() = default;

    BaseCurve(double minValue, double maxValue);

    virtual Position getPosition(double value, double h=0.01) = 0; //not const in case we buffer like in the clothoid

    virtual Position getDerivative(double value) = 0; //not const in case we buffer like in the clothoid

    virtual double getMinValue() const;

    virtual double getMaxValue() const;

    virtual double getLength(const double h=0.01);

    virtual double getLength(double ti, double t_end, double h = 0.01);

    virtual double getValueForLength(double ti, double length, double h);

    double findNearest(Position pos, double h);

    double findNearest(Position pos, double guess, double h, int maxIter, double t_min, double t_max);

    virtual std::vector<double> getParameters() = 0;

    virtual std::vector<double> getFullCurve() {
        auto params = this->getParameters();
        params.insert(params.begin(), static_cast<double>(static_cast<uint64_t>(this->curveType)));
        return params;
    }

    virtual bool isBackward();

    virtual Position getLastPosition() {
        return getPosition(getMaxValue());
    }

};




#endif //BASECURVE_H
