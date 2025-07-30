//
// Created by fogoz on 17/05/2025.
//

#include <utils/ClothoidData.h>

#include <cmath>

#include "utils/Angle.h"
#include "utils/Fresnel.h"


int ClothoidData::build_G1(double const _x0, double const _y0, double const _theta0, double const x1, double const y1,
                           double const theta1, double const tol, double &L, bool const compute_deriv, double L_D[2],
                           double k_D[2], double dk_D[2]) {
    static constexpr double CF[]{
        2.989696028701907, 0.716228953608281, -0.458969738821509, -0.502821153340377, 0.261062141752652,
        -0.045854475238709
    };

    m_x0 = _x0;
    m_y0 = _y0;
    m_theta0 = _theta0;

    // traslazione in (0,0)
    double const dx = x1 - m_x0;
    double const dy = y1 - m_y0;
    double const r = hypot(dx, dy);
    double const phi = atan2(dy, dx);
    double phi0 = m_theta0 - phi;
    double phi1 = theta1 - phi;

    phi0 = WARP_ANGLE(phi0);
    phi1 = WARP_ANGLE(phi1);

    double delta = phi1 - phi0;

    // punto iniziale
    double X{phi0 * M_1_PI};
    double Y{phi1 * M_1_PI};
    double const xy{X * Y};
    Y *= Y;
    X *= X;
    double A{
        (phi0 + phi1) * (CF[0] + xy * (CF[1] + xy * CF[2]) + (CF[3] + xy * CF[4]) * (X + Y) + CF[5] * (X * X + Y * Y))
    };
    // newton
    double g{0}, intC[3], intS[3];
    int niter{0};
    do {
        GeneralizedFresnelCS(3, 2 * A, delta - A, phi0, intC, intS);
        g = intS[0];
        double const dg{intC[2] - intC[1]};
        A -= g / dg;
    } while (++niter <= 100 && fabs(g) > tol);
    ASSERT(fabs(g) <= tol, -1);
    GeneralizedFresnelCS(2 * A, delta - A, phi0, intC[0], intS[0]);
    L = r / intC[0];

    ASSERT(L > 0, -1);
    m_kappa0 = (delta - A) / L;
    m_dk = 2 * A / L / L;

    if (compute_deriv) {
        double const alpha{intC[0] * intC[1] + intS[0] * intS[1]};
        double const beta{intC[0] * intC[2] + intS[0] * intS[2]};
        double const gamma{intC[0] * intC[0] + intS[0] * intS[0]};
        double const tx{intC[1] - intC[2]};
        double const ty{intS[1] - intS[2]};
        double const txy{L * (intC[1] * intS[2] - intC[2] * intS[1])};
        double const omega{L * (intS[0] * tx - intC[0] * ty) - txy};

        delta = intC[0] * tx + intS[0] * ty;

        L_D[0] = omega / delta;
        L_D[1] = txy / delta;

        delta *= L;
        k_D[0] = (beta - gamma - m_kappa0 * omega) / delta;
        k_D[1] = -(beta + m_kappa0 * txy) / delta;

        delta *= L / 2;
        dk_D[0] = (gamma - alpha - m_dk * omega * L) / delta;
        dk_D[1] = (alpha - m_dk * txy * L) / delta;
    }

    return niter;
}

void ClothoidData::reverse(double L) {
    double C, S;
    GeneralizedFresnelCS(m_dk * L * L, m_kappa0 * L, m_theta0, C, S);
    m_x0 += L * C;
    m_y0 += L * S;
    m_theta0 += L * (m_kappa0 + 0.5 * L * m_dk);
    m_kappa0 += L * m_dk;
    m_theta0 += M_PI;
    m_theta0 = WARP_ANGLE(m_theta0);
    m_kappa0 = -m_kappa0;
}

void ClothoidData::eval(double const s, ClothoidData &C) const {
    this->evaluate(s, C.m_theta0, C.m_kappa0, C.m_x0, C.m_y0);
    C.m_dk = m_dk;
}

void ClothoidData::evaluate(double const s, double &theta, double &kappa, double &x, double &y) const {
    double C, S;
    double const sdk = s * m_dk;
    GeneralizedFresnelCS(sdk * s, m_kappa0 * s, m_theta0, C, S);
    x = m_x0 + s * C;
    y = m_y0 + s * S;
    theta = m_theta0 + s * (m_kappa0 + 0.5 * sdk);
    kappa = m_kappa0 + sdk;
}


int ClothoidCurveV2::build_G1(double const x0, double const y0, double const theta0, double const x1, double const y1,
                              double const theta1, double const tol) {
    return m_CD.build_G1(x0, y0, theta0, x1, y1, theta1, tol, m_L);
}

int ClothoidCurveV2::build_G0Kappa(
    double const _x0, double const _y0,
    double const _theta0, double const _kappa0,
    double const x1, double const y1,
    double const tol) {
    return m_CD.build_G0Kappa(_x0, _y0, _theta0, _kappa0, x1, y1, tol, m_L, m_CD.m_dk);
}

double ClothoidCurveV2::kappa_begin() {
    return m_CD.kappa_begin();
}

double ClothoidCurveV2::kappa_end() {
    return m_CD.kappa(m_L);
}

double ClothoidCurveV2::dkappa() {
    return m_CD.m_dk;
}

double ClothoidCurveV2::length() {
    return m_L;
}

double ClothoidCurveV2::theta(double s) {
    return m_CD.theta(s);
}

void ClothoidCurveV2::origin_at(double const s_origin) {
    m_CD.origin_at(s_origin);
}

void ClothoidCurveV2::changeCurvilinearOrigin(double x, double new_L) {
    m_CD.origin_at(x);
    m_L = new_L;
}

void ClothoidCurveV2::reverse() {
    m_CD.reverse(m_L);
}

int ClothoidData::build_G0Kappa(
    double const _x0, double const _y0,
    double const _theta0, double const _kappa0,
    double const x1, double const y1,
    double const tol,
    double &L, double &sigma)
{
    // ---------------------------------------------------------
    // 1. Normalize coordinates (translate + rotate)
    // ---------------------------------------------------------
    m_x0     = _x0;
    m_y0     = _y0;
    m_theta0 = _theta0;
    m_kappa0 = _kappa0;

    double dx  = x1 - m_x0;
    double dy  = y1 - m_y0;
    double r   = hypot(dx, dy);
    double phi = atan2(dy, dx);

    // rotate starting heading into normalized frame
    double phi0 = WARP_ANGLE(m_theta0 - phi);

    // ---------------------------------------------------------
    // 2. Initial guess for Newton
    // ---------------------------------------------------------
    L     = r;     // guess: straight line distance
    sigma = 0.0;   // guess: zero curvature rate

    // ---------------------------------------------------------
    // 3. Newton iteration on (L, sigma)
    // ---------------------------------------------------------
    int niter = 0;
    for (; niter < 50; ++niter) {
        // safety clamp (don’t allow negative lengths)
        if (L < 1e-9) L = 1e-9;

        // -----------------------------------------------------
        // 3a. Setup normalized Fresnel parameters
        // -----------------------------------------------------
        double a    = sigma * 0.5;  // because θ(s) = θ0 + κ0*s + 0.5 σ s^2
        double A    = a * L * L;    // helpful shorthand
        double delta = m_kappa0 * L;     // linear curvature contribution

        // Evaluate generalized Fresnel up to t^2 terms (nk=3)
        double intC[3], intS[3];
        GeneralizedFresnelCS(3, 2*A, delta, phi0, intC, intS);

        // -----------------------------------------------------
        // 3b. Compute residuals: end point position error
        // -----------------------------------------------------
        // intC[0], intS[0] give normalized displacement (scaled by L)
        double X = L * intC[0];
        double Y = L * intS[0];

        double Fx = X - r;
        double Fy = Y;

        if (sqrt(Fx*Fx + Fy*Fy) < tol)
            break;  // converged

        // -----------------------------------------------------
        // 3c. Partial derivatives for Newton Jacobian
        // -----------------------------------------------------

        // ∂(x,y)/∂L = endpoint tangent direction
        double theta_end = phi0 + m_kappa0*L + 0.5*sigma*L*L;
        double dFx_dL = cos(theta_end);
        double dFy_dL = sin(theta_end);

        // ∂(x,y)/∂σ comes from Fresnel t² terms
        // derivative integrals: -1/2 s² sin(θ) for x, +1/2 s² cos(θ) for y
        // intC[2], intS[2] already contain ∫ t² cos(), ∫ t² sin()
        double dFx_dSigma = 0.5 * L*L*L * (-intS[2]); // scale factor L^3/2
        double dFy_dSigma = 0.5 * L*L*L * ( intC[2]);

        // -----------------------------------------------------
        // 3d. Newton update: solve 2x2 system
        // -----------------------------------------------------
        double J00 = dFx_dL;     double J01 = dFx_dSigma;
        double J10 = dFy_dL;     double J11 = dFy_dSigma;

        double det = J00*J11 - J01*J10;
        if (fabs(det) < 1e-14) break; // singular Jacobian

        double dL     = (-Fx * J11 + Fy * J01) / det;
        double dSigma = (-J00 * Fy + J10 * Fx) / det;

        // Damping for stability
        L     += 0.7 * dL;
        sigma += 0.7 * dSigma;
    }

    // ---------------------------------------------------------
    // 4. Store final curvature data
    // ---------------------------------------------------------
    m_kappa0 = _kappa0;
    m_dk     = sigma; // curvature rate

    return niter;
}

