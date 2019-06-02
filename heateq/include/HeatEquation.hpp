#pragma once

#include <cmath>
#include <vector>
#include <cstdint>
#include <functional>

class HeatEquation {
public:
    HeatEquation(
        std::function<double(double, double)> const & iHeatSources,
        std::function<double(double, double)> const & iExactSolution,
        double iDiffusivityCoefficient,
        std::size_t iSpatialBreaksNum,
        std::size_t iTimeLayersNum
    );

    double GetError( std::function<double(double, double)> const & iExactSolution ) const noexcept;

public:
    HeatEquation() = delete;
    ~HeatEquation() = default;

private:
    std::vector<double> ModifiedThomasAlgorithm( std::vector<double> const & iFreePart ) const noexcept;

private:
    double static constexpr m_leftBound      { 0.0 };
    double static constexpr m_rightBound     { 1.0 };
    double static constexpr m_timeStartPoint { 0.0 };
    double static constexpr m_timeEndPoint   { 1.0 };

    std::size_t const m_spatialBreaksNum;
    std::size_t const m_timeLayersNum;
    double      const m_diffusivityCoefficient;

    double const m_spaceStep;
    double const m_timeStep;
    double const m_matrixAboveDiagonalElement;
    double const m_matrixMainDiagonalElement;

    std::vector<std::vector<double>> m_result;
};
