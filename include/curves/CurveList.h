//
// Created by fogoz on 17/05/2025.
//

#ifndef CURVELIST_H
#define CURVELIST_H
#include <vector>
#include <memory>
#include <functional>

class CurveList : public BaseCurve {
private:
    std::vector<std::shared_ptr<BaseCurve>> curves;
    std::vector<double> segment_arc_lengths;  // Cumulative arc lengths
    double total_arc_length = 0.0;

public:
    // Initialize with default min/max (updated in addCurve)
    CurveList() : BaseCurve(0.0, 0.0) {}

    // Add a curve and update total arc length
    void addCurve(const std::shared_ptr<BaseCurve>& curve) {
        if (!curve) throw std::invalid_argument("Null curve provided");
        curves.push_back(curve);

        // Update min/max parametric values (summed for the list)
        minValue = 0.0;  // Always starts at 0 for the list
        maxValue = total_arc_length += curve->getLength(
            curve->getMinValue(),
            curve->getMaxValue(),
            0.001
        );
        segment_arc_lengths.push_back(total_arc_length);
    }

    // --- Required BaseCurve overrides ---
    Position getPosition(double value, double h = 0.01) override {
        auto [curve, local_arc_length] = _getCurveAtArcLength(value);
        double parametric_value = curve->getValueForLength(
            curve->getMinValue(),
            local_arc_length,
            h
        );
        return curve->getPosition(parametric_value, h);
    }

    Position getDerivative(double value) override {
        auto [curve, local_arc_length] = _getCurveAtArcLength(value);
        double parametric_value = curve->getValueForLength(
            curve->getMinValue(),
            local_arc_length,
            0.001
        );
        return curve->getDerivative(parametric_value);
    }

    // --- Optional overrides for better accuracy ---
    double getLength(double ti, double t_end, double h = 0.01) override {
        ti = std::clamp(ti, minValue, maxValue);
        t_end = std::clamp(t_end, minValue, maxValue);
        return t_end - ti;  // Arc length = parametric value for the list
    }

    double getValueForLength(double ti, double length, double h = 0.01) override {
        ti = std::clamp(ti, minValue, maxValue);
        return std::clamp(ti + length, minValue, maxValue);
    }

    // --- Helper methods ---
    std::pair<std::shared_ptr<BaseCurve>, double> _getCurveAtArcLength(double arc_length) {
        if (curves.empty()) throw std::runtime_error("CurveList is empty");
        arc_length = std::clamp(arc_length, minValue, maxValue);

        // Binary search for O(log N) lookup
        auto it = std::upper_bound(segment_arc_lengths.begin(), segment_arc_lengths.end(), arc_length);
        size_t segment_idx = std::distance(segment_arc_lengths.begin(), it);
        if (segment_idx >= curves.size()) segment_idx = curves.size() - 1;

        double segment_start = (segment_idx == 0) ? 0.0 : segment_arc_lengths[segment_idx - 1];
        return {curves[segment_idx], arc_length - segment_start};
    }
};
#endif //CURVELIST_H
