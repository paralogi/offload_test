#include <iostream>
#include <iomanip>
#include <type_traits>
#include <vector>
#include <chrono>
#include <cmath>

#if _OPENMP
#include <omp.h>
#else
// ---
void omp_set_num_threads( int ) {}
int omp_get_num_threads() { return 1; }
int omp_get_max_threads() { return 1; }
int omp_get_thread_num() { return 0; }
int omp_get_num_procs() { return 1; }
// ---
void omp_set_default_device( int ) {}
int omp_get_default_device() { return 0; }
int omp_get_num_devices() { return 0; }
int omp_get_num_teams() { return 1; }
int omp_get_team_num() { return 0; }
// ---
int omp_is_initial_device() { return 0; }
int omp_get_initial_device() { return 0; }
int omp_get_max_task_priority() { return 0; }
#endif

#pragma omp declare target

constexpr double radian = M_PI / 180.;
constexpr double degree = 1 / radian;

template< int O >
using Orientation = std::integral_constant< int, O >;
using Ox = Orientation< 0 >;
using Oy = Orientation< 1 >;
using Oz = Orientation< 2 >;

template< typename Real >
inline auto eval( Real val ) {
    return val;
}

template< typename Real >
inline auto select( bool if_cond, Real then_val, Real else_val ) {
    return if_cond ? then_val : else_val;
}

#pragma omp end declare target

/**
 *
 */
struct Config {
    double origin_logitude;
    double origin_latitude;
    double origin_altitude;
    double scale_factor;

    static Config create() {
        Config config;
        config.origin_logitude = 0;
        config.origin_latitude = 0;
        config.origin_altitude = 0;
        config.scale_factor = 0.75;
        return config;
    }
};

#pragma omp declare target

/**
 *
 */
template< typename Scalar >
struct Point3 {
    Scalar storage[ 3 ];

    static Point3< Scalar > create( Scalar x, Scalar y, Scalar z ) {
        Point3< Scalar > point;
        point.storage[ 0 ] = x;
        point.storage[ 1 ] = y;
        point.storage[ 2 ] = z;
        return point;
    }

    auto& x() {
        return storage[ 0 ];
    }

    auto& y() {
        return storage[ 1 ];
    }

    auto& z() {
        return storage[ 2 ];
    }

    template< int O >
    auto& get( Orientation< O > o ) {
        return storage[ o ];
    }

    template< int O1, int O2 >
    void swap( Orientation< O1 > o1, Orientation< O2 > o2 ) {
        std::swap( get( o1 ), get( o2 ) );
    }

    template< int O1, int O2, typename Real >
    void rotate( Orientation< O1 > o1, Orientation< O2 > o2, Real cos, Real sin ) {
        auto v1 = get( o1 );
        auto v2 = get( o2 );
        get( o1 ) = cos * v1 + sin * v2;
        get( o2 ) = -sin * v1 + cos * v2;
    }
};

/**
 *
 */
struct Georef {
    Point3< double > origin;
    double sin_lon0;
    double cos_lon0;
    double sin_lat0;
    double cos_lat0;
    double scale_factor;
    double major_radius;
    double minor_radius;
    double polar_radius;
    double normal_radius;
    double flattening;
    double eccentricity1;
    double eccentricity2;

    static Georef create( Config& config ) {
        double longitude = config.origin_logitude * radian;
        double latitude = config.origin_latitude * radian;
        double altitude = config.origin_altitude;
        double flattening = 1. / 298.257223563;
        double major_radius = 6378137. / config.scale_factor;
        double minor_radius = major_radius * ( 1. - flattening );
        double polar_radius = major_radius / ( 1. - flattening );
        double eccentricity1 = flattening * ( 2. - flattening );
        double eccentricity2 = eccentricity1 / ( 1. - eccentricity1 );
        double sin_lon0 = std::sin( longitude );
        double cos_lon0 = std::cos( longitude );
        double sin_lat0 = std::sin( latitude );
        double cos_lat0 = std::cos( latitude );
        double normal_radius = polar_radius / std::sqrt( 1. + eccentricity2 * cos_lat0 * cos_lat0 );

        Georef georef;
        georef.origin = Point3< double >::create( longitude, latitude, altitude );
        georef.sin_lon0 = sin_lon0;
        georef.cos_lon0 = cos_lon0;
        georef.sin_lat0 = sin_lat0;
        georef.cos_lat0 = cos_lat0;
        georef.scale_factor = config.scale_factor;
        georef.major_radius = major_radius;
        georef.minor_radius = minor_radius;
        georef.polar_radius = polar_radius;
        georef.normal_radius = normal_radius;
        georef.flattening = flattening;
        georef.eccentricity1 = eccentricity1;
        georef.eccentricity2 = eccentricity2;
        return georef;
    }

    /**
     * projection epsg:4326 to epsg:4978
     */
    template< typename Object >
    Georef& geod2ecef( Object& object ) {
        auto longitude = eval( ( origin.x() + ( object.x() - origin.x() ) * scale_factor ) * radian );
        auto latitude = eval( ( origin.y() + ( object.y() - origin.y() ) * scale_factor ) * radian );
        auto cos_lat = eval( cos( latitude ) );
        auto sin_lat = sin( latitude );
        auto normal = eval( polar_radius / sqrt( 1. + eccentricity2 * cos_lat * cos_lat ) );
        auto hplane = eval( ( normal + object.z() ) * cos_lat );
        object.x() = hplane * cos( longitude );
        object.y() = hplane * sin( longitude );
        object.z() = ( object.z() + normal * ( 1 - eccentricity1 ) ) * sin_lat;
        return *this;
    }

    /**
     * projection epsg:4978 to epsg:4326
     */
    template< typename Object >
    Georef& ecef2geod( Object& object ) {
        auto hplane = eval( hypot( object.x(), object.y() ) );
        auto tangent0 = object.z() / hplane * ( 1. + ( eccentricity2 * minor_radius ) / hypot( hplane, object.z() ) );
        auto latitude0 = eval( atan( tangent0 * ( 1. - flattening ) ) );
        auto tangent1 = ( ( object.z() + ( eccentricity2 * minor_radius ) * pow( sin( latitude0 ), 3 ) ) /
                ( hplane - ( eccentricity1 * major_radius ) * pow( cos( latitude0 ), 3 ) ) );
        auto latitude1 = eval( atan( tangent1 * ( 1. - flattening ) ) );
        auto tangent = eval( ( object.z() + ( eccentricity2 * minor_radius ) * pow( sin( latitude1 ), 3 ) ) /
                ( hplane - ( eccentricity1 * major_radius ) * pow( cos( latitude1 ), 3 ) ) );
        auto longitude = atan2( object.y(), object.x() );
        auto latitude = eval( atan( tangent ) );
        auto cos_lat = eval( select( latitude == M_PI_2, 1., cos( latitude ) ) );
        auto sin_lat = select( latitude == 0, 1., sin( latitude ) );
        auto normal = polar_radius / sqrt( 1. + eccentricity2 * cos_lat * cos_lat );
        object.x() = origin.x() + ( longitude * degree - origin.x() ) / scale_factor;
        object.y() = origin.y() + ( latitude * degree - origin.y() ) / scale_factor;
        object.z() = select( abs( tangent ) <= 1, hplane / cos_lat - normal, object.z() / sin_lat - normal * ( 1. - eccentricity1 ) );
        return *this;
    }

    /**
     * projection epsg:4978 to epsg:5819
     */
    template< typename Object >
    Georef& ecef2topo( Object& object ) {
        object.z() += eccentricity1 * normal_radius * sin_lat0;
        object.rotate( Ox{}, Oy{}, cos_lon0, sin_lon0 );
        object.rotate( Oz{}, Ox{}, sin_lat0, cos_lat0 );
        object.swap( Ox{}, Oy{} );
        object.y() = -object.y();
        object.z() -= normal_radius + origin.z();
        return *this;
    }

    /**
     * projection epsg:5819 to epsg:4978
     */
    template< typename Object >
    Georef& topo2ecef( Object& object ) {
        object.z() += normal_radius + origin.z();
        object.y() = -object.y();
        object.swap( Ox{}, Oy{} );
        object.rotate( Oz{}, Ox{}, sin_lat0, -cos_lat0 );
        object.rotate( Ox{}, Oy{}, cos_lon0, -sin_lon0 );
        object.z() -= eccentricity1 * normal_radius * sin_lat0;
        return *this;
    }
};

#pragma omp end declare target

/**
 *
 */
int main( int argc, char** argv ) {
    auto config = Config::create();
    auto georef = Georef::create( config );
    std::cout << std::fixed << std::setprecision( 3 );

    {
        int num_teams, num_threads;

        #pragma omp target teams map(from: num_teams, num_threads)
        {
            num_teams = omp_get_num_teams();
            #pragma omp parallel
            #pragma omp single
            num_threads = omp_get_num_threads();
        }

        std::cout << "default device:" << omp_get_default_device() << "\n";
        std::cout << "initial device:" << omp_get_initial_device() << "\n";
        std::cout << "device count:" << omp_get_num_devices() << "\n";
        std::cout << "teams count: " << num_teams << "\n";
        std::cout << "threads count: " << num_threads << "\n";
    }

    uint32_t count = 1'000'000;
    std::vector< double > storage( count * 3 );

    {
        for ( uint32_t i = 0; i < count; i++ )
            storage[ i ] = -180 + 360 * double( i ) / ( count - 1 );

        for ( uint32_t i = 0; i < count; i++ )
            storage[ i + count ] = -90 + 180 * double( i ) / ( count - 1 );

        std::cout << "in:\n";
        for ( uint32_t i = 0; i < count; i += ( count - 1 ) / 2 )
            std::cout << storage[ i ] << " " << storage[ i + count ] << " " << storage[ i + count * 2 ] << "\n";
    }

    {
        auto data = storage.data();
        auto size = storage.size();
        auto timer = std::chrono::steady_clock::now();

        #pragma omp target map(tofrom: data[:size]), map(to: georef)
        #pragma omp teams distribute parallel for firstprivate(georef)
        for ( uint32_t i = 0; i < count; i++ ) {
            uint32_t xi = i;
            uint32_t yi = i + count;
            uint32_t zi = i + count * 2;
            auto point = Point3< double >::create( data[ xi ], data[ yi ], data[ zi ] );
            georef.geod2ecef( point ).ecef2topo( point ).topo2ecef( point ).ecef2geod( point );
            data[ xi ] = point.x();
            data[ yi ] = point.y();
            data[ zi ] = point.z();
        }

        std::cout << std::chrono::duration_cast< std::chrono::duration< double > >(
                std::chrono::steady_clock::now() - timer ).count() << "s\n";
    }

    {
        std::cout << "out:\n";
        for ( uint32_t i = 0; i < count; i += ( count - 1 ) / 2 )
            std::cout << storage[ i ] << " " << storage[ i + count ] << " " << storage[ i + count * 2 ] << "\n";
    }

    return 0;
}
