//
// Created by fogoz on 17/05/2025.
//

#ifndef CLOTHOIDDATA_H
#define CLOTHOIDDATA_H
#include "Fresnel.h"
#include "utils/Position.h"


struct ClothoidData {
    double m_x0 = 0;
    double m_y0 = 0;
    double m_theta0 = 0;
    double m_kappa0 = 0;
    double m_dk = 0;

public:
    int build_G1(Position pos0, Position pos1, double const tol, double &L);

    int build_G1(double const _x0, double const _y0, double const _theta0, double const x1, double const y1,
                 double const theta1, double const tol, double &L, bool compute_deriv = false, double L_D[2] = nullptr,
                 double k_D[2] = nullptr, double dk_D[2] = nullptr);

    double kappa_begin() {
        return m_kappa0;
    }

    double kappa(double s) {
        return m_kappa0 + s * m_dk;
    }

    double theta(double s) {
        return m_theta0 + s * (m_kappa0 + 0.5 * s * m_dk);
    }

    void origin_at(double const s_origin) {
        double C, S;
        double const sdk = s_origin * m_dk;
        GeneralizedFresnelCS(sdk * s_origin, m_kappa0 * s_origin, m_theta0, C, S);
        m_x0 += s_origin * C;
        m_y0 += s_origin * S;
        m_theta0 += s_origin * (m_kappa0 + 0.5 * sdk);
        m_kappa0 += sdk;
    }

    void reverse(double);

    void eval(double const s, ClothoidData &C) const;

    void evaluate(double const s, double &theta, double &kappa, double &x, double &y) const;

    int build_G0Kappa(double const _x0, double const _y0, double const _theta0, double const _kappa0, double const x1,
                      double const y1, double const tol, double &L, double &sigma);
};

struct ClothoidCurveV2 {
    ClothoidData m_CD;
    double m_L = 0;

public:
    int build_G1(double const x0, double const y0, double const theta0, double const x1, double const y1,
                 double const theta1, double const tol = 1e-12);

    int build_G0Kappa(double const _x0, double const _y0, double const _theta0, double const _kappa0, double const x1,
                      double const y1, double const tol=1e-7);

    double kappa_begin();

    double kappa_end();

    double dkappa();

    double length();

    double theta(double s);

    void origin_at(double const s_origin);

    void changeCurvilinearOrigin(double x, double newL);

    void reverse();
};

#ifndef ARDUINO
#include <iostream>

inline std::ostream &operator<<(std::ostream &stream, const ClothoidCurveV2 &curve) {
    stream << "x0 : " << curve.m_CD.m_x0 << " y0 : " << curve.m_CD.m_y0 << " theta0 : " << curve.m_CD.m_theta0 <<
            " kappa0 : " << curve.m_CD.m_kappa0 << " dk : " << curve.m_CD.m_dk << " L : " << curve.m_L;
    return stream;
}


#endif


#endif //CLOTHOIDDATA_H
