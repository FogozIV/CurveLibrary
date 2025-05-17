//
// Created by fogoz on 16/05/2025.
//

#ifndef FRESNEL_H
#define FRESNEL_H
#include <cassert>
#include <cmath>

#define A_THRESHOLD 0.1
#define A_SERIE_SIZE 3
static constexpr double m_1_sqrt_pi{0.564189583547756286948079451561};

static constexpr double fn[] = {
    0.49999988085884732562,
    1.3511177791210715095,
    1.3175407836168659241,
    1.1861149300293854992,
    0.7709627298888346769,
    0.4173874338787963957,
    0.19044202705272903923,
    0.06655998896627697537,
    0.022789258616785717418,
    0.0040116689358507943804,
    0.0012192036851249883877
  };

static constexpr double fd[] = {
    1.0,
    2.7022305772400260215,
    4.2059268151438492767,
    4.5221882840107715516,
    3.7240352281630359588,
    2.4589286254678152943,
    1.3125491629443702962,
    0.5997685720120932908,
    0.20907680750378849485,
    0.07159621634657901433,
    0.012602969513793714191,
    0.0038302423512931250065
  };

static constexpr double gn[] = {
    0.50000014392706344801,
    0.032346434925349128728,
    0.17619325157863254363,
    0.038606273170706486252,
    0.023693692309257725361,
    0.007092018516845033662,
    0.0012492123212412087428,
    0.00044023040894778468486,
    -8.80266827476172521e-6,
    -1.4033554916580018648e-8,
    2.3509221782155474353e-10
  };

static constexpr double gd[] = {
    1.0,
    2.0646987497019598937,
    2.9109311766948031235,
    2.6561936751333032911,
    2.0195563983177268073,
    1.1167891129189363902,
    0.57267874755973172715,
    0.19408481169593070798,
    0.07634808341431248904,
    0.011573247407207865977,
    0.0044099273693067311209,
    -0.00009070958410429993314
  };
  //!
  //! **Purpose:**
  //!
  //! Compute Fresnel integrals C(x) and S(x)
  //!
  //! \f[
  //!   S(x) = \int_0^x \sin t^2 \,\mathrm{d} t, \qquad
  //!   C(x) = \int_0^x \cos t^2 \,\mathrm{d} t
  //! \f]
  //!
  //! **Example:**
  //!
  //! | \f$ x \f$ | \f$ C(x) \f$  |  \f$ S(x) \f$ |
  //! | :-------: | :-----------: | :-----------: |
  //! | 0.0       | 0.00000000    |  0.00000000   |
  //! | 0.5       | 0.49234423    |  0.06473243   |
  //! | 1.0       | 0.77989340    |  0.43825915   |
  //! | 1.5       | 0.44526118    |  0.69750496   |
  //! | 2.0       | 0.48825341    |  0.34341568   |
  //! | 2.5       | 0.45741301    |  0.61918176   |
  //!
  //! **Adapted from:**
  //!
  //! - *William J. Thompson*, Atlas for computing mathematical functions :
  //!   an illustrated guide for practitioners, with programs in C and Mathematica,
  //!   Wiley, 1997.
  //!
  //! **Author:**
  //!
  //! - *Venkata Sivakanth Telasula*,
  //!   email: sivakanth.telasula@gmail.com,
  //!   date: August 11, 2005
  //!
  //! \param[in]  y Argument of \f$C(y)\f$ and \f$S(y)\f$
  //! \param[out] C \f$C(x)\f$
  //! \param[out] S \f$S(x)\f$
  //!
static constexpr void FresnelCS( double y, double & C, double & S ) {

    constexpr double eps { 1E-15 };
    double const x{ y > 0 ? y : -y };

    if ( x < 1.0 ) {
      double term;

      double const s { M_PI_2*(x*x) };
      double const t { -s*s };

      // Cosine integral series
      double twofn   { 0.0 };
      double fact    { 1.0 };
      double denterm { 1.0 };
      double numterm { 1.0 };
      double sum     { 1.0 };
      do {
        twofn   += 2.0;
        fact    *= twofn*(twofn-1.0);
        denterm += 4.0;
        numterm *= t;
        term     = numterm/(fact*denterm);
        sum     += term;
      } while ( abs(term) > eps*abs(sum) );

      C = x*sum;

      // Sine integral series
      twofn   = 1.0;
      fact    = 1.0;
      denterm = 3.0;
      numterm = 1.0;
      sum     = 1.0/3.0;
      do {
        twofn   += 2.0;
        fact    *= twofn*(twofn-1.0);
        denterm += 4.0;
        numterm *= t;
        term     = numterm/(fact*denterm);
        sum     += term;
      } while ( abs(term) > eps*abs(sum) );

      S = M_PI_2*sum*(x*x*x);

    } else if ( x < 6.0 ) {

      // Rational approximation for f
      double sumn{ 0.0 };
      double sumd{ fd[11] };
      for ( int k=10; k >= 0; --k ) {
        sumn = fn[k] + x*sumn;
        sumd = fd[k] + x*sumd;
      }
      double const f{ sumn/sumd };

      // Rational approximation for g
      sumn = 0.0;
      sumd = gd[11];
      for ( int k=10; k >= 0; --k ) {
        sumn = gn[k] + x*sumn;
        sumd = gd[k] + x*sumd;
      }
      double const g    { sumn/sumd };
      double const U    { M_PI_2*(x*x) };
      double const SinU { sin(U) };
      double const CosU { cos(U) };
      C = 0.5 + f*SinU - g*CosU;
      S = 0.5 - f*CosU - g*SinU;

    } else {

      double absterm;

      // x >= 6; asymptotic expansions for  f  and  g

      double const s { M_PI*x*x };
      double const t { -1/(s*s) };

      // Expansion for f
      double       numterm {-1.0 };
      double       term    { 1.0 };
      double       sum     { 1.0 };
      double       oldterm { 1.0 };
      double const eps10   { 0.1 * eps };

      do {
        numterm += 4.0;
        term    *= numterm*(numterm-2.0)*t;
        sum     += term;
        absterm  = abs(term);
        assert(oldterm >= absterm);
        oldterm  = absterm;
      } while ( absterm > eps10 * abs(sum) );

      double const f{ sum / (M_PI*x) };

      //  Expansion for  g
      numterm = -1.0;
      term    =  1.0;
      sum     =  1.0;
      oldterm =  1.0;

      do {
        numterm += 4.0;
        term    *= numterm*(numterm+2.0)*t;
        sum     += term;
        absterm  = abs(term);
        assert(oldterm >= absterm);
        oldterm  = absterm;
      } while ( absterm > eps10 * abs(sum) );

      double       g    { M_PI*x }; g = sum/(g*g*x);
      double const U    { M_PI_2*(x*x) };
      double const SinU { sin(U) };
      double const CosU { cos(U) };
      C = 0.5 + f*SinU - g*CosU;
      S = 0.5 - f*CosU - g*SinU;

    }
    if ( y < 0 ) { C = -C; S = -S; }
  }

static constexpr void FresnelCS(
  int   const nk,
  double const t,
  double       C[],
  double       S[]
) {
    FresnelCS(t,C[0],S[0]);
    if ( nk > 1 ) {
        double const tt { M_PI_2*(t*t) };
        double const ss { sin(tt) };
        double const cc { cos(tt) };
        C[1] = ss*M_1_PI;
        S[1] = (1-cc)*M_1_PI;
        if ( nk > 2 ) {
            C[2] = (t*ss-S[0])*M_1_PI;
            S[2] = (C[0]-t*cc)*M_1_PI;
        }
    }
}
static constexpr void evalXYaLarge(
    double const a,
    double const b,
    double &     X,
    double &     Y
  ) {
    double const s    = a > 0 ? +1 : -1;
    double const absa = abs(a);
    double const z    = m_1_sqrt_pi*sqrt(absa);
    double const ell  = s*b*m_1_sqrt_pi/sqrt(absa);
    double const g    = -0.5*s*(b*b)/absa;
    double const cg   = cos(g)/z;
    double const sg   = sin(g)/z;

    double Cl, Sl, Cz, Sz;
    FresnelCS( ell,   Cl, Sl );
    FresnelCS( ell+z, Cz, Sz );

    double const dC0{ Cz - Cl };
    double const dS0{ Sz - Sl };

    X = cg * dC0 - s * sg * dS0;
    Y = sg * dC0 + s * cg * dS0;
  }

  // -------------------------------------------------------------------------
  // nk max 3
  static constexpr void evalXYaLarge(
    int   const nk,
    double const a,
    double const b,
    double       X[],
    double       Y[]
  ) {

    assert(nk < 4 && nk > 0);

    double const s    { static_cast<double>(a > 0 ? +1 : -1) };
    double const absa { abs(a) };
    double const z    { m_1_sqrt_pi*sqrt(absa) };
    double const ell  { s*b*m_1_sqrt_pi/sqrt(absa) };
    double const g    { -0.5*s*(b*b)/absa };
    double       cg   { cos(g)/z };
    double       sg   { sin(g)/z };

    double Cl[3], Sl[3], Cz[3], Sz[3];

    FresnelCS( nk, ell,   Cl, Sl );
    FresnelCS( nk, ell+z, Cz, Sz );

    double const dC0 { Cz[0] - Cl[0] };
    double const dS0 { Sz[0] - Sl[0] };
    X[0] = cg * dC0 - s * sg * dS0;
    Y[0] = sg * dC0 + s * cg * dS0;
    if ( nk > 1 ) {
      cg /= z;
      sg /= z;
      double const dC1 { Cz[1] - Cl[1] };
      double const dS1 { Sz[1] - Sl[1] };
      double       DC  { dC1-ell*dC0   };
      double       DS  { dS1-ell*dS0   };
      X[1] = cg * DC - s * sg * DS;
      Y[1] = sg * DC + s * cg * DS;
      if ( nk > 2 ) {
        double const dC2{ Cz[2] - Cl[2] };
        double const dS2{ Sz[2] - Sl[2] };
        DC   = dC2+ell*(ell*dC0-2*dC1);
        DS   = dS2+ell*(ell*dS0-2*dS1);
        cg   = cg/z;
        sg   = sg/z;
        X[2] = cg * DC - s * sg * DS;
        Y[2] = sg * DC + s * cg * DS;
      }
    }
  }

static double LommelReduced( double const mu, double const nu, double const b ) {
    double tmp{ 1/((mu+nu+1)*(mu-nu+1)) };
    double res{ tmp };
    for ( int n = 1; n <= 100; ++n ) {
        tmp *= (-b/(2*n+mu-nu+1)) * (b/(2*n+mu+nu+1));
        res += tmp;
        if ( abs(tmp) < abs(res) * 1e-50 ) break;
    }
    return res;
}

static constexpr void evalXYazero(
    int const nk,
    double const b,
    double X[],
    double Y[]
) {
    double const sb{sin(b)};
    double const cb{cos(b)};
    double const b2{b * b};
    if (abs(b) < 1e-3) {
        X[0] = 1 - (b2 / 6) * (1 - (b2 / 20) * (1 - (b2 / 42)));
        Y[0] = (b / 2) * (1 - (b2 / 12) * (1 - (b2 / 30)));
    } else {
        X[0] = sb / b;
        Y[0] = (1 - cb) / b;
    }
    // use recurrence in the stable part
    int m{static_cast<int>(floor(2 * b))};
    if (m >= nk) m = nk - 1;
    if (m < 1) m = 1;
    for (int k{1}; k < m; ++k) {
        X[k] = (sb - k * Y[k - 1]) / b;
        Y[k] = (k * X[k - 1] - cb) / b;
    }
    //  use Lommel for the unstable part
    if (m < nk) {
        double const A{b * sb};
        double const D{sb - b * cb};
        double const B{b * D};
        double const C{-b2 * sb};
        double rLa{LommelReduced(m + 0.5, 1.5, b)};
        double rLd{LommelReduced(m + 0.5, 0.5, b)};
        for (int k{m}; k < nk; ++k) {
            double const rLb{LommelReduced(k + 1.5, 0.5, b)};
            double const rLc{LommelReduced(k + 1.5, 1.5, b)};
            X[k] = (k * A * rLa + B * rLb + cb) / (1 + k);
            Y[k] = (C * rLc + sb) / (2 + k) + D * rLd;
            rLa = rLc;
            rLd = rLb;
        }
    }
}

// -------------------------------------------------------------------------
// -------------------------------------------------------------------------

static constexpr void evalXYaSmall(
    double const a,
    double const b,
    int const p,
    double &X,
    double &Y
) {
    assert(p < 11 && p > 0);

    double X0[43], Y0[43];

    int const nkk{4 * p + 3}; // max 43
    evalXYazero(nkk, b, X0, Y0);

    X = X0[0] - (a / 2) * Y0[2];
    Y = Y0[0] + (a / 2) * X0[2];

    double t{1};
    double const aa{-a * a / 4}; // controllare!
    for (int n{1}; n <= p; ++n) {
        t *= aa / (2 * n * (2 * n - 1));
        double const bf{a / (4 * n + 2)};
        int const jj{4 * n};
        X += t * (X0[jj] - bf * Y0[jj + 2]);
        Y += t * (Y0[jj] + bf * X0[jj + 2]);
    }
}

// -------------------------------------------------------------------------

static constexpr void evalXYaSmall(
    int const nk,
    double const a,
    double const b,
    int const p,
    double X[],
    double Y[]
) {
    int nkk{nk + 4 * p + 2}; // max 45
    double X0[45], Y0[45];

    assert(nkk < 46);

    evalXYazero(nkk, b, X0, Y0);

    for (int j = 0; j < nk; ++j) {
        X[j] = X0[j] - (a / 2) * Y0[j + 2];
        Y[j] = Y0[j] + (a / 2) * X0[j + 2];
    }

    double t{1};
    double const aa{-a * a / 4}; // controllare!
    for (int n{1}; n <= p; ++n) {
        t *= aa / (2 * n * (2 * n - 1));
        double const bf{a / (4 * n + 2)};
        for (int j{0}; j < nk; ++j) {
            int const jj{4 * n + j};
            X[j] += t * (X0[jj] - bf * Y0[jj + 2]);
            Y[j] += t * (Y0[jj] + bf * X0[jj + 2]);
        }
    }
}

static constexpr void GeneralizedFresnelCS(
    int const nk,
    double const a,
    double const b,
    double const c,
    double intC[],
    double intS[]
) {
    assert(nk > 0 && nk < 4);

    if (abs(a) < A_THRESHOLD) evalXYaSmall(nk, a, b, A_SERIE_SIZE, intC, intS);
    else evalXYaLarge(nk, a, b, intC, intS);

    double const cosc{cos(c)};
    double const sinc{sin(c)};

    for (int k{0}; k < nk; ++k) {
        double const xx{intC[k]};
        double const yy{intS[k]};
        intC[k] = xx * cosc - yy * sinc;
        intS[k] = xx * sinc + yy * cosc;
    }
}

static constexpr void GeneralizedFresnelCS(
    double const a,
    double const b,
    double const c,
    double &intC,
    double &intS
) {
    double xx, yy;
    if (abs(a) < A_THRESHOLD) evalXYaSmall(a, b, A_SERIE_SIZE, xx, yy);
    else evalXYaLarge(a, b, xx, yy);

    double const cosc{cos(c)};
    double const sinc{sin(c)};

    intC = xx * cosc - yy * sinc;
    intS = xx * sinc + yy * cosc;
}
#endif //FRESNEL_H
