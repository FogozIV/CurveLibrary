//
// Created by fogoz on 16/05/2025.
//

#ifndef G2SOLVE3ARC_H
#define G2SOLVE3ARC_H

#include <cmath>

#include "ClothoidData.h"
#include "utils/Position.h"

class G2Solve3Arc {
    ClothoidCurveV2 m_segment0, m_segment1, m_segmentM;
    double m_tolerance = 1e-10;
    int m_max_iter = 100;

    double m_x0 = 0;
    double m_y0 = 0;
    double m_theta0 = 0;
    double m_kappa0 = 0;
    double m_dkappa0 = 0;
    double m_x1 = 0;
    double m_y1 = 0;
    double m_theta1 = 0;
    double m_kappa1 = 0;
    double m_dkappa1 = 0;

    double m_phi = 0;
    double m_Lscale = 0;
    double m_th0 = 0;
    double m_th1 = 0;
    double m_s0 = 0;
    double m_s1 = 0;

    double m_K0 = 0, m_K1 = 0, m_c0 = 0, m_c1 = 0, m_c2 = 0, m_c3 = 0, m_c4 = 0;
    double m_c5 = 0, m_c6 = 0, m_c7 = 0, m_c8 = 0, m_c9 = 0, m_c10 = 0;
    double m_c11 = 0, m_c12 = 0, m_c13 = 0, m_c14 = 0;


    void evalFJ(double const vars[2], double F[2], double J[2][2]);

    void evalF(double const vars[2], double F[2]);

    void buildSolution(double sM, double thM);

    int solve(double sM_guess, double thM_guess);

public:

    int build(Position start, Position end);

    int build(double x0, double y0, double theta0, double kappa0, double x1, double y1, double theta1, double kappa1,
              double Dmax = 0, double dmax = 0);

    int build_fixed_length(double s0, double x0, double y0, double theta0, double kappa0, double s1, double x1, double y1,
              double theta1, double kappa1);

    [[nodiscard]] ClothoidCurveV2 getSegment0() const {
        return m_segment0;
    }

    [[nodiscard]] ClothoidCurveV2 getSegment1() const {
        return m_segment1;
    }

    [[nodiscard]] ClothoidCurveV2 getSegmentMiddle() const {
        return m_segmentM;
    }

    void reverse();
};


#endif //G2SOLVE3ARC_H
