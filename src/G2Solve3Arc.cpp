//
// Created by fogoz on 16/05/2025.
//

#include "../include/utils/G2Solve3Arc.h"

#include "curves/ClothoidCurve.h"
#include "utils/Angle.h"
#include "utils/ClothoidData.h"
#include "utils/Fresnel.h"
#include "utils/Solver2x2.h"


void G2Solve3Arc::evalFJ(double const vars[2], double F[2], double J[2][2]) {
    double const sM{vars[0]};
    double const thM{vars[1]};

    double const dsM{1.0 / (m_c13 + (m_c14 + sM) * sM)};
    double const dsMsM{dsM * sM};
    double const dK0{dsM * (m_c0 * thM + sM * (m_c1 * thM + m_c2 - sM * m_K0) + m_c3)};
    double const dK1{dsM * (m_c0 * thM + sM * (m_c4 * thM + m_c5 + sM * m_K1) + m_c6)};
    double const dKM{dsMsM * (thM * (m_c7 - 2 * sM) + m_c8 * sM + m_c9)};
    double const KM{dsMsM * (m_c10 * thM + m_c11 * sM + m_c12)};

    double X0[3], Y0[3],
            X1[3], Y1[3],
            XMp[3], YMp[3],
            XMm[3], YMm[3];
    GeneralizedFresnelCS(3, dK0, m_K0, m_th0, X0, Y0);
    GeneralizedFresnelCS(3, dK1, -m_K1, m_th1, X1, Y1);
    GeneralizedFresnelCS(3, dKM, KM, thM, XMp, YMp);
    GeneralizedFresnelCS(3, dKM, -KM, thM, XMm, YMm);

    // in the standard problem dx = 2, dy = 0
    double const t0{XMp[0] + XMm[0]};
    double const t1{YMp[0] + YMm[0]};
    F[0] = m_s0 * X0[0] + m_s1 * X1[0] + sM * t0 - 2;
    F[1] = m_s0 * Y0[0] + m_s1 * Y1[0] + sM * t1 - 0;

    // calcolo J(F)
    double const dsM2{dsM * dsM};
    double const g0{-(2 * sM + m_c14) * dsM2};
    double const g1{(m_c13 - sM * sM) * dsM2};
    double const g2{sM * (sM * m_c14 + 2 * m_c13) * dsM2};

    double const dK0_sM{(m_c0 * thM + m_c3) * g0 + (m_c1 * thM + m_c2) * g1 - m_K0 * g2};
    double const dK1_sM{(m_c0 * thM + m_c6) * g0 + (m_c4 * thM + m_c5) * g1 + m_K1 * g2};
    double const dKM_sM{(m_c7 * thM + m_c9) * g1 + (m_c8 - 2 * thM) * g2};
    double const KM_sM{(m_c10 * thM + m_c12) * g1 + m_c11 * g2};

    double const dK0_thM{(m_c0 + m_c1 * sM) * dsM};
    double const dK1_thM{(m_c0 + m_c4 * sM) * dsM};
    double const dKM_thM{(m_c7 - 2 * sM) * dsMsM};
    double const KM_thM{m_c10 * dsMsM};

    // coeff fresnel per f_j per lo jacobiano
    double const f0{-0.5 * m_s0 * Y0[2]};
    double const f1{-0.5 * m_s1 * Y1[2]};
    double const f2{-0.5 * sM * (YMm[2] + YMp[2])};
    double const f3{sM * (YMm[1] - YMp[1])};
    double const f4{0.5 * m_s0 * X0[2]};
    double const f5{0.5 * m_s1 * X1[2]};
    double const f6{0.5 * sM * (XMm[2] + XMp[2])};
    double const f7{sM * (XMp[1] - XMm[1])};

    J[0][0] = f0 * dK0_sM + f1 * dK1_sM + f2 * dKM_sM + f3 * KM_sM + t0;
    J[0][1] = f0 * dK0_thM + f1 * dK1_thM + f2 * dKM_thM + f3 * KM_thM - sM * t1;
    J[1][0] = f4 * dK0_sM + f5 * dK1_sM + f6 * dKM_sM + f7 * KM_sM + t1;
    J[1][1] = f4 * dK0_thM + f5 * dK1_thM + f6 * dKM_thM + f7 * KM_thM + sM * t0;
}

void G2Solve3Arc::evalF(double const vars[2], double F[2]) {
    double const sM{vars[0]};
    double const thM{vars[1]};

    double const dsM{1.0 / (m_c13 + (m_c14 + sM) * sM)};
    double const dK0{dsM * (m_c0 * thM + sM * (m_c1 * thM - m_K0 * sM + m_c2) + m_c3)};
    double const dK1{dsM * (m_c0 * thM + sM * (m_c4 * thM + m_K1 * sM + m_c5) + m_c6)};
    double const dKM{dsM * sM * (thM * (m_c7 - 2 * sM) + m_c8 * sM + m_c9)};
    double const KM{dsM * sM * (m_c10 * thM + m_c11 * sM + m_c12)};

    double X0, Y0, X1, Y1, XMp, YMp, XMm, YMm;
    GeneralizedFresnelCS(dK0, m_K0, m_th0, X0, Y0);
    GeneralizedFresnelCS(dK1, -m_K1, m_th1, X1, Y1);
    GeneralizedFresnelCS(dKM, KM, thM, XMp, YMp);
    GeneralizedFresnelCS(dKM, -KM, thM, XMm, YMm);

    // in the standard problem dx = 2, dy = 0
    F[0] = m_s0 * X0 + m_s1 * X1 + sM * (XMm + XMp) - 2;
    F[1] = m_s0 * Y0 + m_s1 * Y1 + sM * (YMm + YMp) - 0;
}

void G2Solve3Arc::buildSolution(double sM, double thM) {
    double const dsM { 1.0 / (m_c13+(m_c14+sM)*sM) };
    double       dK0 { dsM*(m_c0*thM + sM*(m_c1*thM - m_K0*sM + m_c2) + m_c3) };
    double       dK1 { dsM*(m_c0*thM + sM*(m_c4*thM + m_K1*sM + m_c5) + m_c6) };
    double       dKM { dsM*sM*(m_c7*thM + sM*(m_c8 - 2*thM) + m_c9) };
    double       KM  { dsM*sM*(m_c10*thM + m_c11*sM + m_c12) };

    double xa, ya, xmL, ymL;
    GeneralizedFresnelCS( dK0,  m_K0, m_th0, xa,  ya  );
    GeneralizedFresnelCS( dKM,   -KM,   thM, xmL, ymL );

    double const xM { m_s0 * xa + sM * xmL - 1 };
    double const yM { m_s0 * ya + sM * ymL };

    // rovescia trasformazione standard
    double const L0{ m_s0/m_Lscale };
    double const L1{ m_s1/m_Lscale };
    double const LM{ sM/m_Lscale   };

    dK0 *= pow(m_Lscale/m_s0, 2);
    dK1 *= pow(m_Lscale/m_s1, 2);
    dKM *= pow(m_Lscale/sM, 2);
    KM  *= m_Lscale/sM;
    //th0 = theta0 - phi;
    //th1 = theta1 - phi;
    m_segment0 = { m_x0, m_y0, m_phi+m_th0, m_kappa0, dK0, L0 };
    m_segment1 = { m_x1, m_y1, m_phi+m_th1, m_kappa1, dK1, L1 };
    m_segment1.changeCurvilinearOrigin(-L1, L1);

    // la trasformazione inversa da [-1,1] a (x0,y0)-(x1,y1)
    // g(x,y) = RotInv(phi)*(1/lambda*[X;Y] - [xbar;ybar]) = [x;y]

    double const C  { cos(m_phi) };
    double const S  { sin(m_phi) };
    double const dx { (xM + 1) / m_Lscale };
    double const dy { yM / m_Lscale };
    m_segmentM = {
      m_x0 + C * dx - S * dy,
      m_y0 + C * dy + S * dx,
      thM + m_phi, KM, dKM, 2*LM
    };
    m_segmentM.changeCurvilinearOrigin(-LM, 2*LM);
}

int G2Solve3Arc::solve(double sM_guess, double thM_guess) {
    double X[2];
    X[0] = sM_guess;
    X[1] = thM_guess;


    int iter{0};
    bool converged{false};
    Solve2x2 solver;
    do {
        double J[2][2];
        double d[2];
        double F[2];
        evalFJ(X, F, J);
        double const lenF{hypot(F[0], F[1])};
        converged = lenF < m_tolerance;
        if (converged || !solver.factorize(J)) break;
        solver.solve(F, d);
#if 1
        // use undamped Newton
        X[0] -= d[0];
        X[1] -= d[1];
#else
            double FF[2], dd[2], XX[2];
            // Affine invariant Newton solver
            double nd = hypot( d[0], d[1] );
            bool step_found = false;
            double tau = 2;
            do {
                tau  /= 2;
                XX[0] = X[0]-tau*d[0];
                XX[1] = X[1]-tau*d[1];
                evalF(XX, FF);
                solver.solve(FF, dd);
                step_found = hypot( dd[0], dd[1] ) <= (1-tau/2)*nd + 1e-6;
                //&& XX[0] > 0; // && XX[0] > X[0]/4 && XX[0] < 4*X[0];
                //&& XX[1] > thmin && XX[1] < thmax;
            } while ( tau > 1e-6 && !step_found );
            if ( !step_found ) break;
            X[0] = XX[0];
            X[1] = XX[1];
#endif
    } while (++iter < m_max_iter);

    // re-check solution
    if (converged)
        converged = FP_INFINITE != fpclassify(X[0]) &&
                    FP_NAN != fpclassify(X[0]) &&
                    FP_INFINITE != fpclassify(X[1]) &&
                    FP_NAN != fpclassify(X[1]);
    if ( converged ) buildSolution(X[0], X[1]);
    return converged ? iter : -1;
}


int G2Solve3Arc::build(Position start, Position end) {
    return build(start.getX(), start.getY(), start.getAngle().toRadians(), start.getCurvature(),
                 end.getX(), end.getY(), end.getAngle().toRadians(), end.getCurvature());
}


int G2Solve3Arc::build(double x0, double y0, double theta0, double kappa0, double x1, double y1, double theta1,
                       double kappa1, double Dmax, double dmax) {
    m_x0 = x0;
    m_y0 = y0;
    m_theta0 = theta0;
    m_kappa0 = kappa0;
    m_x1 = x1;
    m_y1 = y1;
    m_theta1 = theta1;
    m_kappa1 = kappa1;

    // transform to reference frame
    double const dx{m_x1 - m_x0};
    double const dy{m_y1 - m_y0};
    m_phi = atan2(dy, dx);
    m_Lscale = 2 / hypot(dx, dy);

    m_th0 = m_theta0 - m_phi;
    m_th1 = m_theta1 - m_phi;

    // put in range
    m_th0 = WARP_ANGLE(m_th0);
    m_th1 = WARP_ANGLE(m_th1);

    m_K0 = (m_kappa0 / m_Lscale); // k0
    m_K1 = (m_kappa1 / m_Lscale); // k1

    if (Dmax <= 0) Dmax = M_PI;
    if (dmax <= 0) dmax = M_PI / 8;

    if (Dmax > 2 * M_PI) Dmax = M_PI * 2;
    if (dmax > M_PI / 4) dmax = M_PI / 4;

    // compute guess G1

    ClothoidCurveV2 SG;
    SG.build_G1( -1, 0, m_th0, 1, 0, m_th1 );

    double const kA { SG.kappa_begin() };
    double const kB { SG.kappa_end() };
    double const dk { abs(SG.dkappa()) };
    double const L3 { SG.length()/3 };

    double tmp{0.5 * abs(m_K0 - kA) / dmax};
    m_s0 = L3;
    if (tmp * m_s0 > 1) m_s0 = 1 / tmp;
    tmp = (abs(m_K0 + kA) + m_s0 * dk) / (2 * Dmax);
    if (tmp * m_s0 > 1) m_s0 = 1 / tmp;

    tmp = 0.5 * abs(m_K1 - kB) / dmax;
    m_s1 = L3;
    if (tmp * m_s1 > 1) m_s1 = 1 / tmp;
    tmp = (abs(m_K1 + kB) + m_s1 * dk) / (2 * Dmax);
    if (tmp * m_s1 > 1) m_s1 = 1 / tmp;

    double const dth{abs(m_th0 - m_th1) / (2 * M_PI)};
    double const scale{pow(cos(pow(dth, 4) * M_PI_2), 3)};
    m_s0 *= scale;
    m_s1 *= scale;

    double const L{(3 * L3 - m_s0 - m_s1) / 2};
    double const thM{SG.theta(m_s0+L)};
    m_th0 = SG.theta(0);
    m_th1 = SG.theta(SG.length());

    // setup

    m_K0 *= m_s0;
    m_K1 *= m_s1;

    double const t0{2 * m_th0 + m_K0};
    double const t1{2 * m_th1 - m_K1};

    m_c0 = m_s0 * m_s1;
    m_c1 = 2 * m_s0;
    m_c2 = 0.25 * ((m_K1 - 6 * (m_K0 + m_th0) - 2 * m_th1) * m_s0 - 3 * m_K0 * m_s1);
    m_c3 = -m_c0 * (m_K0 + m_th0);
    m_c4 = 2 * m_s1;
    m_c5 = 0.25 * ((6 * (m_K1 - m_th1) - m_K0 - 2 * m_th0) * m_s1 + 3 * m_K1 * m_s0);
    m_c6 = m_c0 * (m_K1 - m_th1);
    m_c7 = -0.5 * (m_s0 + m_s1);
    m_c8 = m_th0 + m_th1 + 0.5 * (m_K0 - m_K1);
    m_c9 = 0.25 * (t1 * m_s0 + t0 * m_s1);
    m_c10 = 0.5 * (m_s1 - m_s0);
    m_c11 = 0.5 * (m_th1 - m_th0) - 0.25 * (m_K0 + m_K1);
    m_c12 = 0.25 * (t1 * m_s0 - t0 * m_s1);
    m_c13 = 0.5 * m_s0 * m_s1;
    m_c14 = 0.75 * (m_s0 + m_s1);
    return solve(L, thM);
}

int G2Solve3Arc::build_fixed_length(double s0, double x0, double y0, double theta0, double kappa0, double s1, double x1,
                                    double y1, double theta1, double kappa1) {
    m_x0 = x0;
    m_y0 = y0;
    m_theta0 = theta0;
    m_kappa0 = kappa0;
    m_x1 = x1;
    m_y1 = y1;
    m_theta1 = theta1;
    m_kappa1 = kappa1;

    // transform to reference frame
    double const dx{m_x1 - m_x0};
    double const dy{m_y1 - m_y0};
    m_phi = atan2(dy, dx);
    m_Lscale = 2 / hypot(dx, dy);

    m_th0 = m_theta0 - m_phi;
    m_th1 = m_theta1 - m_phi;

    // put in range
    m_th0 = WARP_ANGLE(m_th0);
    m_th1 = WARP_ANGLE(m_th1);

    m_K0 = (m_kappa0 / m_Lscale); // k0
    m_K1 = (m_kappa1 / m_Lscale); // k1

    // compute guess G1
    ClothoidCurveV2 SG;
    SG.build_G1( -1, 0, m_th0, 1, 0, m_th1 );

    m_s0 = s0 * m_Lscale;
    m_s1 = s1 * m_Lscale;

    double const L{(SG.length() - m_s0 - m_s1) / 2};
    double const thM{SG.theta(m_s0 + L)};
    m_th0 = SG.theta(0);
    m_th1 = SG.theta(SG.length());

    // setup

    m_K0 *= m_s0;
    m_K1 *= m_s1;

    double const t0{2 * m_th0 + m_K0};
    double const t1{2 * m_th1 - m_K1};

    m_c0 = m_s0 * m_s1;
    m_c1 = 2 * m_s0;
    m_c2 = 0.25 * ((m_K1 - 6 * (m_K0 + m_th0) - 2 * m_th1) * m_s0 - 3 * m_K0 * m_s1);
    m_c3 = -m_c0 * (m_K0 + m_th0);
    m_c4 = 2 * m_s1;
    m_c5 = 0.25 * ((6 * (m_K1 - m_th1) - m_K0 - 2 * m_th0) * m_s1 + 3 * m_K1 * m_s0);
    m_c6 = m_c0 * (m_K1 - m_th1);
    m_c7 = -0.5 * (m_s0 + m_s1);
    m_c8 = m_th0 + m_th1 + 0.5 * (m_K0 - m_K1);
    m_c9 = 0.25 * (t1 * m_s0 + t0 * m_s1);
    m_c10 = 0.5 * (m_s1 - m_s0);
    m_c11 = 0.5 * (m_th1 - m_th0) - 0.25 * (m_K0 + m_K1);
    m_c12 = 0.25 * (t1 * m_s0 - t0 * m_s1);
    m_c13 = 0.5 * m_s0 * m_s1;
    m_c14 = 0.75 * (m_s0 + m_s1);

    return solve(L, thM);
}

std::shared_ptr<ClothoidCurve> G2Solve3Arc::getSegment0Curve() const {
    return std::make_shared<ClothoidCurve>(m_segment0);
}

std::shared_ptr<ClothoidCurve> G2Solve3Arc::getSegment1Curve() const {
    return std::make_shared<ClothoidCurve>(m_segment1);
}

std::shared_ptr<ClothoidCurve> G2Solve3Arc::getSegmentMiddleCurve() const {
    return std::make_shared<ClothoidCurve>(m_segmentM);
}

void G2Solve3Arc::reverse() {
    std::swap(m_segment0, m_segment1);
    m_segmentM.reverse();
    m_segment0.reverse();
    m_segment1.reverse();
}
