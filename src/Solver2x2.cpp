//
// Created by fogoz on 16/05/2025.
//

#include "../include/utils/Solver2x2.h"

#include <cmath>
#include <math.h>

bool
  Solve2x2::factorize( double A[2][2] ) {
    // full pivoting
    double Amax = fabs(A[0][0]);
    double tmp  = fabs(A[0][1]);
    int   ij{0};
    if ( tmp > Amax ) { ij = 1; Amax = tmp; }
    tmp = fabs(A[1][0]);
    if ( tmp > Amax ) { ij = 2; Amax = tmp; }
    tmp = fabs(A[1][1]);
    if ( tmp > Amax ) { ij = 3; Amax = tmp; }
    if ( fabs(Amax) < 1e-8) return false;
    if ( (ij&0x01) == 0x01 ) { j[0] = 1; j[1] = 0; }
    else                     { j[0] = 0; j[1] = 1; }
    if ( (ij&0x02) == 0x02 ) { i[0] = 1; i[1] = 0; }
    else                     { i[0] = 0; i[1] = 1; }
    // apply factorization
    LU[0][0] = A[i[0]][j[0]];
    LU[0][1] = A[i[0]][j[1]];
    LU[1][0] = A[i[1]][j[0]];
    LU[1][1] = A[i[1]][j[1]];

    LU[1][0] /= LU[0][0];
    LU[1][1] -= LU[1][0]*LU[0][1];
    // check for singularity
    singular = fabs( LU[1][1] ) < epsi;
    return true;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  //!
  //! Solve the linear system \f$ Ax=b \f$ with
  //! \f$ A \f$ stored amd facted with a previous call
  //! of method ``factorize``.
  //!
  //! \param[in]  b the r.h.s. of \f$ Ax=b \f$
  //! \param[out] x the solution of \f$ Ax=b \f$
  //! \return `true` if solution is found
  //!
  bool
  Solve2x2::solve( double const b[2], double x[2] ) const {
    if ( singular ) {
      // L^+ Pb
      double tmp = (b[i[0]] + LU[1][0]*b[i[1]]) /
                      ( (1+pow(LU[1][0],2) ) * ( pow(LU[0][0],2)+pow(LU[0][1],2) ) );
      x[j[0]] = tmp*LU[0][0];
      x[j[1]] = tmp*LU[0][1];
      // check consistency
      tmp = (LU[0][0]*x[j[0]]+LU[0][1]*x[j[1]]);
      return hypot( b[i[0]]-tmp, b[i[1]]+tmp*LU[1][0] ) < hypot(b[0],b[1])*epsi;
    }
    // non singular
    // L^(-1) Pb
    x[j[0]] = b[i[0]];
    x[j[1]] = b[i[1]]-LU[1][0]*x[j[0]];
    // U^(-1) x
    x[j[1]] /= LU[1][1];
    x[j[0]]  = (x[j[0]]-LU[0][1]*x[j[1]])/LU[0][0];
    return FP_INFINITE != fpclassify(x[0]) &&
           FP_NAN      != fpclassify(x[0]) &&
           FP_INFINITE != fpclassify(x[1]) &&
           FP_NAN      != fpclassify(x[1]);
  }