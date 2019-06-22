#include <HeatEquation.hpp>

#include <fstream>
#include <atomic>
#include <mutex>

#include <heq.pb.h>

HeatEquation::HeatEquation(
   std::function<double(double, double)> iHeatSources,
   std::function<double(double, double)> iExactSolution,
   double const iDiffusivityCoefficient,
   std::size_t const iSpatialBreaksNum,
   std::size_t const iTimeLayersNum,
   std::size_t const iThreadCount
)
    : m_source{ std::move( iHeatSources ) }
    , m_exact{ std::move( iExactSolution ) }
    , m_spatialBreaksNum{ iSpatialBreaksNum + 1 }
    , m_timeLayersNum{ iTimeLayersNum + 1 }
    , m_diffusivityCoefficient{ iDiffusivityCoefficient }
    , m_spaceStep{ ( m_rightBound - m_leftBound ) / m_spatialBreaksNum }
    , m_timeStep{ ( m_timeEndPoint - m_timeStartPoint ) / m_timeLayersNum }
    , m_matrixAboveDiagonalElement{ ( -m_diffusivityCoefficient * m_timeStep ) / std::pow( m_spaceStep, 2 ) }
    , m_matrixMainDiagonalElement{ 1.0 + ( ( 2.0 * m_diffusivityCoefficient * m_timeStep ) / std::pow( m_spaceStep, 2 ) ) }
    , m_result( m_timeLayersNum )
    , m_taskArena( iThreadCount, 1 )
{
    std::vector<double> taFreePart;
    taFreePart.reserve( m_spatialBreaksNum - 2 );

    std::vector<double> & firstTimeLayer = m_result.front();
    firstTimeLayer.resize( m_spatialBreaksNum );

    for( std::size_t index { 0 }; index < m_spatialBreaksNum; ++index )
        firstTimeLayer[ index ] = m_exact( index * m_spaceStep, m_timeStartPoint );

    double const courantNum { ( m_diffusivityCoefficient * m_timeStep ) / std::pow( m_spaceStep, 2 ) };
    for( std::size_t timeLayerIndex { 1 }; timeLayerIndex < m_result.size(); ++timeLayerIndex )
    {
        taFreePart.resize( m_spatialBreaksNum - 2 );
        m_taskArena.execute( [&]{
            tbb::parallel_for( std::size_t{ 0 }, taFreePart.size(), [&]( std::size_t const index) {
                taFreePart[ index ] = m_result[ timeLayerIndex - 1 ][ index + 1 ] + m_timeStep *
                        m_source( m_leftBound + ( index + 1 ) * m_spaceStep, m_timeStartPoint + timeLayerIndex * m_timeStep );
            } );
        });

        double const leftBound  { m_exact( m_leftBound, m_timeStartPoint + timeLayerIndex * m_timeStep ) };
        double const rightBound { m_exact( m_rightBound, m_timeStartPoint + timeLayerIndex * m_timeStep ) };
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

std::vector<double> HeatEquation::ModifiedThomasAlgorithm( std::vector<double> const & iFreePart ) noexcept
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
    return result;
}


double HeatEquation::GetError( std::function<double(double, double)> const & iExactSolution ) noexcept
{
    std::atomic< double > error{ 0.0 };
    m_taskArena.execute( [&]{
        tbb::parallel_for( std::size_t{ 1 }, m_timeLayersNum, [&]( std::size_t const timeIndex ) {
            double timeLayerError { 0.0 };
            for( std::size_t spaceIndex{ 1 }; spaceIndex < m_spatialBreaksNum - 1; ++spaceIndex )
            {
                double const exactSolution =
                    iExactSolution( m_leftBound + spaceIndex * m_spaceStep, m_timeStartPoint + timeIndex * m_timeStep );
                double const approximateSolution = m_result[ timeIndex ][ spaceIndex ];
                double absDifference  = std::abs( exactSolution - approximateSolution );
                timeLayerError  = ( absDifference > timeLayerError ) ? absDifference : timeLayerError;
            }

            if( timeLayerError > error )
                error = timeLayerError;
        } );
    } );

    return error * m_spaceStep * m_timeStep;
}

void HeatEquation::WriteResult( fs::path const & iPath ) noexcept
{
    GOOGLE_PROTOBUF_VERIFY_VERSION;
    std::mutex mutex;

    HEQSerialize::Mesh mesh;
    mesh.set_type( HEQSerialize::Mesh_Type::Mesh_Type_APPROXIMATION );
    mesh.set_spatial_splits_num( m_spatialBreaksNum );
    mesh.set_time_layers_num( m_timeLayersNum );

    m_taskArena.execute( [&]{
        tbb::parallel_for( std::size_t{ 0 }, m_result.size(), [&]( std::size_t const index ) {
            HEQSerialize::TimeLayer timeLayer;
            for( double const value : m_result[ index ] )
                timeLayer.add_value( value );

            std::unique_lock lock{ mutex };
            (*mesh.mutable_mesh())[ index ] = timeLayer;
        } );
    } );

    std::ofstream ofs{ iPath / "approximation.pb", std::ios::binary | std::ios::trunc };
    mesh.SerializeToOstream( &ofs );
    ofs.close();
}