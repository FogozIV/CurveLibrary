//
// Created by fogoz on 16/05/2025.
//

#ifndef SOLVER2X2_H
#define SOLVER2X2_H



class Solve2x2 {
    int i[2], j[2];
    double LU[2][2];
    double epsi;
    bool      singular;

public:

    Solve2x2() : epsi(1e-10) {}
    bool factorize( double A[2][2] );
    bool solve( double const b[2], double x[2] ) const;
};


#endif //SOLVER2X2_H
