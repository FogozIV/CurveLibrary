//
// Created by fogoz on 21/05/2025.
//

#ifndef SOLVER3X3_H
#define SOLVER3X3_H
#include <cmath>
#include <cstring>
// Simple 3x3 solver using Gaussian elimination
inline bool solve3x3(const double A[3][3], const double b[3], double x[3]) {
    // Make a copy to preserve original matrix
    double mat[3][3];
    double vec[3];
    memcpy(mat, A, sizeof(mat));
    memcpy(vec, b, sizeof(vec));
    
    // Forward elimination
    for (int col = 0; col < 2; ++col) {
        // Find pivot row
        int pivot = col;
        for (int row = col+1; row < 3; ++row) {
            if (fabs(mat[row][col]) > fabs(mat[pivot][col])) {
                pivot = row;
            }
        }
        
        // Swap rows if needed
        if (pivot != col) {
            for (int i = 0; i < 3; ++i) {
                std::swap(mat[col][i], mat[pivot][i]);
            }
            std::swap(vec[col], vec[pivot]);
        }
        
        // Eliminate column
        for (int row = col+1; row < 3; ++row) {
            double factor = mat[row][col] / mat[col][col];
            for (int i = col; i < 3; ++i) {
                mat[row][i] -= factor * mat[col][i];
            }
            vec[row] -= factor * vec[col];
        }
    }
    
    // Back substitution
    for (int row = 2; row >= 0; --row) {
        double sum = vec[row];
        for (int col = row+1; col < 3; ++col) {
            sum -= mat[row][col] * x[col];
        }
        
        if (fabs(mat[row][row]) < 1e-6) return false; // Singular matrix
        x[row] = sum / mat[row][row];
    }
    
    return true;
}
#endif //SOLVER3X3_H
