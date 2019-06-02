#include <HeatEquation.hpp>

HeatEquation::HeatEquation(
   std::function<double(double, double)> const & iHeatSources,
   std::function<double(double, double)> const & iExactSolution,
   double const iDiffusivityCoefficient,
   std::size_t const iSpatialBreaksNum,
   std::size_t const iTimeLayersNum
)
    : m_spatialBreaksNum{ iSpatialBreaksNum + 1 }
    , m_timeLayersNum{ iTimeLayersNum + 1 }
    , m_diffusivityCoefficient{ iDiffusivityCoefficient }
    , m_spaceStep{ ( m_rightBound - m_leftBound ) / m_spatialBreaksNum }
    , m_timeStep{ ( m_timeEndPoint - m_timeStartPoint ) / m_timeLayersNum }
    , m_matrixAboveDiagonalElement{ ( -m_diffusivityCoefficient * m_timeStep ) / std::pow( m_spaceStep, 2 ) }
    , m_matrixMainDiagonalElement{ 1.0 + ( ( 2.0 * m_diffusivityCoefficient * m_timeStep ) / std::pow( m_spaceStep, 2 ) ) }
    , m_result( m_timeLayersNum )
{
    std::vector<double> taFreePart;
    taFreePart.reserve( m_spatialBreaksNum - 2 );

    std::vector<double> & firstTimeLayer = m_result.front();
    firstTimeLayer.resize( m_spatialBreaksNum );

    for( std::size_t index { 0 }; index < m_spatialBreaksNum; ++index )
        firstTimeLayer[ index ] = iExactSolution( index * m_spaceStep, m_timeStartPoint );

    double const courantNum { ( m_diffusivityCoefficient * m_timeStep ) / std::pow( m_spaceStep, 2 ) };
    for( std::size_t timeLayerIndex { 1 }; timeLayerIndex < m_result.size(); ++timeLayerIndex )
    {
        for( std::size_t spaceIndex { 1 }; spaceIndex < m_spatialBreaksNum - 1; ++spaceIndex )
            taFreePart.emplace_back( m_result[ timeLayerIndex - 1 ][ spaceIndex ] + m_timeStep
                * iHeatSources( m_leftBound + spaceIndex * m_spaceStep, m_timeStartPoint + timeLayerIndex * m_timeStep )
            );

        double const leftBound  { iExactSolution( m_leftBound, m_timeStartPoint + timeLayerIndex * m_timeStep ) };
        double const rightBound { iExactSolution( m_rightBound, m_timeStartPoint + timeLayerIndex * m_timeStep ) };
        taFreePart.front() += courantNum * leftBound;
        taFreePart.back()  += courantNum * rightBound;

        auto mid = this->ModifiedThomasAlgorithm( taFreePart );
        taFreePart.clear();

        std::vector<double> & timeLayer = m_result[ timeLayerIndex ];
        timeLayer.reserve( m_spatialBreaksNum );

        timeLayer.emplace_back( leftBound );
        std::move( mid.begin(), mid.end(), std::back_inserter( timeLayer ) );
        timeLayer.emplace_back( rightBound );
    }
}

std::vector<double> HeatEquation::ModifiedThomasAlgorithm( std::vector<double> const & iFreePart ) const noexcept
{
    std::size_t const n = iFreePart.size();
    std::vector<double> alpha( n - 1 ), beta( n - 1 ), result( n );

    alpha.back() = - m_matrixAboveDiagonalElement / m_matrixMainDiagonalElement;
    beta.back() = iFreePart.back() / m_matrixMainDiagonalElement;

    for( auto index { n - 2 }; index > 0; --index )
    {
        double commonFactor = 1.0 / ( m_matrixMainDiagonalElement + m_matrixAboveDiagonalElement * alpha[ index ] );
        beta[ index - 1 ] = ( iFreePart[ index ] - beta[ index ] * m_matrixAboveDiagonalElement ) * commonFactor;
        alpha[ index - 1 ] = -m_matrixAboveDiagonalElement * commonFactor;
    }

    result.front() = ( iFreePart.front() - m_matrixMainDiagonalElement * beta.front() )
        / ( m_matrixMainDiagonalElement + m_matrixAboveDiagonalElement * alpha.front() );

    for( std::size_t index { 0 }; index < n - 1; ++index )
        result[ index + 1 ] = alpha[ index ] * result[ index ] + beta[ index ];
    return std::move( result );
}


double HeatEquation::GetError( std::function<double(double, double)> const & iExactSolution ) const noexcept
{
    double maxError { 0.0 };
    for ( std::size_t time_iter{ 1 }; time_iter < m_timeLayersNum; ++time_iter )
        for ( std::size_t space_iter{ 1 }; space_iter < m_spatialBreaksNum - 1; ++space_iter )
        {
            double exactSolution
                = iExactSolution( m_leftBound + space_iter * m_spaceStep, m_timeStartPoint + time_iter * m_timeStep );
            double approximateSolution  = m_result[time_iter][space_iter];
            double absDifference  = std::abs( exactSolution - approximateSolution );
            maxError  = ( absDifference > maxError ) ? absDifference : maxError;
        }
    return maxError * m_spaceStep * m_timeStep;
}
