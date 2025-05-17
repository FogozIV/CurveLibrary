//
// Created by fogoz on 17/05/2025.
//

#include "../include/utils/ClothoidData.h"

#include <cmath>

#include "utils/Angle.h"
#include "utils/Fresnel.h"


int ClothoidData::build_G1(
    double const _x0,
    double const _y0,
    double const _theta0,
    double const x1,
    double const y1,
    double const theta1,
    double const tol,
    double &     L,
    bool      const compute_deriv,
    double       L_D[2],
    double       k_D[2],
    double       dk_D[2]
  ) {
    static constexpr double CF[]{
      2.989696028701907,   0.716228953608281,
      -0.458969738821509, -0.502821153340377,
      0.261062141752652,  -0.045854475238709
    };

    m_x0     = _x0;
    m_y0     = _y0;
    m_theta0 = _theta0;

    // traslazione in (0,0)
    double const dx   = x1 - m_x0;
    double const dy   = y1 - m_y0;
    double const r    = hypot( dx, dy );
    double const phi  = atan2( dy, dx );
    double       phi0 = m_theta0 - phi;
    double       phi1 = theta1 - phi;

    phi0 -= 2*M_PI*round(phi0/(2*M_PI));
    phi1 -= 2*M_PI*round(phi1/(2*M_PI));

    WARP_ANGLE(phi0);
    WARP_ANGLE(phi1);

    double delta = phi1 - phi0;

    // punto iniziale
    double       X  { phi0*M_1_PI };
    double       Y  { phi1*M_1_PI };
    double const xy { X*Y };
    Y *= Y; X *= X;
    double A{ (phi0+phi1) * ( CF[0] + xy*(CF[1] + xy*CF[2]) +
                               ( CF[3]+xy*CF[4])*(X+Y) + CF[5]*(X*X+Y*Y) ) };
    // newton
    double g{0}, intC[3], intS[3];
    int   niter{0};
    do {
      GeneralizedFresnelCS( 3, 2*A, delta-A, phi0, intC, intS );
      g   = intS[0];
      double const dg{ intC[2] - intC[1] };
      A  -= g / dg;
    } while ( ++niter <= 10 && abs(g) > tol );

    assert(abs(g) <= tol);
    GeneralizedFresnelCS( 2*A, delta-A, phi0, intC[0], intS[0] );
    L = r/intC[0];

    assert( L > 0);
    m_kappa0 = (delta-A)/L;
    m_dk     = 2*A/L/L;

    if ( compute_deriv ) {

      double const alpha { intC[0]*intC[1] + intS[0]*intS[1] };
      double const beta  { intC[0]*intC[2] + intS[0]*intS[2] };
      double const gamma { intC[0]*intC[0] + intS[0]*intS[0] };
      double const tx    { intC[1]-intC[2] };
      double const ty    { intS[1]-intS[2] };
      double const txy   { L*(intC[1]*intS[2]-intC[2]*intS[1]) };
      double const omega { L*(intS[0]*tx-intC[0]*ty) - txy };

      delta = intC[0]*tx + intS[0]*ty;

      L_D[0] = omega/delta;
      L_D[1] = txy/delta;

      delta *= L;
      k_D[0] = (beta-gamma-m_kappa0*omega)/delta;
      k_D[1] = -(beta+m_kappa0*txy)/delta;

      delta  *= L/2;
      dk_D[0] = (gamma-alpha-m_dk*omega*L)/delta;
      dk_D[1] = (alpha-m_dk*txy*L)/delta;
    }

    return niter;
  }

void ClothoidData::reverse(double L) {
  double C, S;
  GeneralizedFresnelCS( m_dk*L*L, m_kappa0*L, m_theta0, C, S );
  m_x0     += L*C;
  m_y0     += L*S;
  m_theta0 += L*(m_kappa0+0.5*L*m_dk);
  m_kappa0 += L*m_dk;
  m_theta0 += M_PI;
  WARP_ANGLE(m_theta0);
  m_kappa0  = -m_kappa0;
}

void ClothoidData::eval(double const s, ClothoidData &C) const {
  this->evaluate( s, C.m_theta0, C.m_kappa0, C.m_x0, C.m_y0 );
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
