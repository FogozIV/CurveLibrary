//
// Created by fogoz on 10/05/2025.
//

#include "curves/ClothoidCurve.h"

#include <tuple>
#include "utils/Solver.h"

ClothoidCurve::ClothoidCurve(Position start, double startCurvature, double curvatureRate,
                             double length) : BaseCurve(0, length), start(start), startCurvature(startCurvature),
                                              curvatureRate(curvatureRate), length(length) {
}

ClothoidCurve::ClothoidCurve(Position start, double initialCurvature, std::array<double, 2> cl) : ClothoidCurve(
    start, initialCurvature, cl[0], cl[1]) {
}

ClothoidCurve::ClothoidCurve(Position start, std::array<double, 3> kcl) : ClothoidCurve(start, kcl[0], kcl[1], kcl[2]) {
}

std::shared_ptr<ClothoidCurve> ClothoidCurve::getClothoidCurve(Position start, Angle endAngle, double initialCurvature,
                                                               double length) {
    Angle delta = endAngle - start.getAngle();
    double curvatureRate = delta.toRadians() - initialCurvature * length;
    curvatureRate /= pow(length, 2);
    curvatureRate *= 2;
    return std::make_shared<ClothoidCurve>(start, initialCurvature, curvatureRate, length);
}

std::shared_ptr<ClothoidCurve> ClothoidCurve::getClothoidCurveDelta(Position start, Angle deltaAngle,
                                                                    double initialCurvature, double length) {
    return getClothoidCurve(start, start.getAngle() + deltaAngle, initialCurvature, length);
}

std::optional<std::shared_ptr<ClothoidCurve> > ClothoidCurve::getClothoidCurve(
    Position start, Position end, double initialCurvature) {
    return nullptr;
}

double ClothoidCurve::getCurvatureRateAngle(Angle delta, double length, double initialCurvature) {
    return (delta.toRadians() - initialCurvature * length) / pow(length, 2) * 2;
}

Position ClothoidCurve::localToGlobal(const Position &local) const {
    double cos_angle = cos(start.getAngle().toRadians());
    double sin_angle = sin(start.getAngle().toRadians());

    // Rotate and translate position
    double global_x = start.getX() + local.getX() * cos_angle - local.getY() * sin_angle;
    double global_y = start.getY() + local.getX() * sin_angle + local.getY() * cos_angle;

    Angle global_angle = local.getAngle() + start.getAngle();

    return Position(global_x, global_y, global_angle);
}

Position ClothoidCurve::getPosition(double value, double h) {
    if (previousPos.has_value() && std::get<0>(previousPos.value()) == value) {
        return localToGlobal(std::get<1>(previousPos.value()));
    }
    double t_i = 0;
    std::vector<double> y = {0.0, 0.0, 0.0};
    if (previousPos.has_value()) {
        t_i = std::get<2>(previousPos.value());
        if (t_i >= value) {
            t_i = 0;
        } else {
            auto previous_pos = std::get<1>(previousPos.value());
            y[0] = previous_pos.getX();
            y[1] = previous_pos.getY();
            y[2] = previous_pos.getAngle().toRadians();
        }
    }
    runge_kutta_4(t_i, y, value, h, [this](double t, std::vector<double> &y, std::vector<double> &dydt) {
        this->clothoidODE(t, y, dydt);
    });

    auto local = Position(y[0], y[1], Angle::fromRadians(y[2]));
    previousPos = std::make_tuple(value, local, t_i);
    /*
    double startAngle = start.getAngle().toRadians();
    double cosAngle = cos(startAngle);
    double sinAngle = sin(startAngle);

    double rotatedX = y[0] * cosAngle - y[1] * sinAngle;
    double rotatedY = y[0] * sinAngle + y[1] * cosAngle;
    auto pos = Position(start.getX() + rotatedX, start.getY() + rotatedY, Angle::fromRadians(y[2]));
    */
    auto pos = localToGlobal(local);
    return pos;
}

void ClothoidCurve::clothoidODE(double t, std::vector<double> &y, std::vector<double> &dydt) const {
    double current_curvature = startCurvature + curvatureRate * t;

    dydt[0] = cos(y[2]); // dx/dt
    dydt[1] = sin(y[2]); // dy/dt
    dydt[2] = current_curvature; // dtheta/dt
}

Position ClothoidCurve::getDerivative(double value) {
    std::vector<double> y(3), dydt(3);
    Position pos = getPosition(value, 0.01); // small step
    double angle = pos.getAngle().toRadians();

    double current_curvature = startCurvature + curvatureRate * value;

    return {cos(angle), sin(angle), Angle::fromRadians(current_curvature)};
}

double ClothoidCurve::getLength(double ti, double t_end, double h) {
    return t_end - ti;
}

double ClothoidCurve::getValueForLength(double ti, double length, double h) {
    return std::min(ti + length, this->length);
}

double integralXk(double a, double b, double c, double k=0, double h=0.001) {
    double t0 = 0.0;
    std::vector<double> y = {0};
    runge_kutta_4(t0, y, 1, h, [a,b,c, k](double t, std::vector<double> &y, std::vector<double> &dydt) {
        dydt[0] = pow(t, k) * cos(a/2*pow(t,2) + b*t + c);
    });
    return y[0];
}
double integralYk(double a, double b, double c, double k=0, double h=0.001) {
    double t0 = 0.0;
    std::vector<double> y = {0};
    runge_kutta_4(t0, y, 1, h, [a,b,c,k](double t, std::vector<double> &y, std::vector<double> &dydt) {
        dydt[0] = pow(t, k) * sin(a/2*pow(t,2) + b*t + c);
    });
    return y[0];
}

std::optional<std::array<double, 3>> buildClothoid(Position begin, Position end) {
    double kappa = 0.0;
    double dkappa = 0.0;
    double L = 0.0;
    // Step 1: Convert to relative coordinates
    auto dPos = end - begin;
    double dx = dPos.getX();
    double dy = dPos.getY();
    auto theta0 = begin.getAngle();
    auto theta1 = end.getAngle();
    double r = hypot(dx, dy);
    Angle phi = Angle::fromRadians(atan2(dy, dx));

    // Step 2: Normalize angles relative to phi
    auto phi0 = (theta0 - phi).warpAngle().toRadians();
    auto phi1 = (theta1 - phi).warpAngle().toRadians();
    auto delta = (phi1 - phi0);

    // Step 3: Define the function g(A) = Y(2A, delta - A, phi0)
    auto g = [phi0, delta](double A) {
        auto res = delta - A;
        return integralYk((2*A), res, phi0);
    };

    // Step 4: Initial guess for A
    auto A = (3 * (phi0 + phi1)); // Simple initial guess

    // Step 5: Newton-Raphson to find root of g(A)
    double tolerance = 1e-10;
    int max_iter = 100;
    double g_val = g(A);

    for (int i = 0; i < max_iter && fabs(g_val) > tolerance; i++) {
        // Finite difference approximation of g'(A)
        double h = 1e-8;
        double g_prime = (g(A + h) - g(A - h)) / (2 * h);

        if (fabs(g_prime) < 1e-15) break; // Avoid division by zero

        A = A - g_val / g_prime;
        g_val = g(A);
    }

    // Step 6: Compute final parameters
    auto X_val = integralXk((2*A), delta - A, phi0);
    //auto Y_val = integralYk((2*A), delta - A, phi0, 0.01);

    if (fabs(X_val) < 1e-15) return std::nullopt; // No solution

    L = r / X_val;
    kappa = ((delta - A) / L);
    dkappa = (2 * A / (L * L));

    return std::array<double, 3>({kappa, dkappa, L});
}

std::optional<std::array<double, 2>> buildClothoid(Position begin, Position end, double kappa0) {
    // Step 1: Convert to relative coordinates
    auto dPos = end - begin;
    double dx = dPos.getX();
    double dy = dPos.getY();
    auto theta0 = begin.getAngle();
    auto theta1 = end.getAngle();
    double r = hypot(dx, dy);
    Angle phi = Angle::fromRadians(atan2(dy, dx));

    // Step 2: Normalize angles relative to phi
    auto phi0 = (theta0 - phi).warpAngle().toRadians();
    auto phi1 = (theta1 - phi).warpAngle().toRadians();
    auto delta = (phi1 - phi0);


    auto computeError = [=](double L)->std::pair<double, double> {
        double dkappa = 2 * (delta - kappa0 * L) / (L * L);

        double X = integralXk(dkappa * L * L, kappa0 * L, phi0);
        double Y = integralYk(dkappa * L * L, kappa0 * L, phi0);

        double x_err = X * L - r;
        double y_err = Y * L;
        return {hypot(x_err, y_err), dkappa};
    };

    // Step 4: Initial guess for A
    auto L_min = 1e-6; // Simple initial guess
    auto L_max = 2300;
    double L = r;

    // Step 5: Hybrid root finding (combination of bisection and Newton)
    double tolerance = 1e-10;
    int max_iter = 1000;
    double error = 0;
    double dkappa = 0;

    for (int i = 0; i < max_iter; ++i) {
        auto [current_error, current_dkappa] = computeError(L);
        std::cout << "L: " << L << " Error: " << current_error << std::endl;

        if (current_error < tolerance) {
            return std::array<double, 2>{current_dkappa, L};
        }

        // Try Newton step
        double h = std::max(L * 1e-4, 1e-8);
        auto [error_plus, _] = computeError(L + h);
        auto [error_minus, __] = computeError(L - h);

        double derivative = (error_plus - error_minus) / (2 * h);

        if (std::abs(derivative) > 1e-15) {
            double L_new = L - 0.5 * current_error / derivative;
            // Keep within bounds
            if (L_new > L_min && L_new < L_max) {
                L = L_new;
                continue;
            }
        }

        // Fall back to bisection if Newton fails
        auto err_low = computeError(L_min);
        auto err_high = computeError(L_max);

        if (err_low.first * err_high.first > 0) {
            // No solution in interval
            break;
        }

        L = (L_min + L_max) / 2;
        auto [mid_err, mid_dkappa] = computeError(L);

        if (mid_err < tolerance) {
            return std::array<double, 2>{mid_dkappa, L};
        }

        if (mid_err * err_low.first < 0) {
            L_max = L;
        } else {
            L_min = L;
        }
    }

    return std::nullopt;
}
