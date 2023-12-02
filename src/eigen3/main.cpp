#include <cmath>
#include <chrono>
#include <type_traits>
#include <iostream>
#include <iomanip>

#include <Eigen/Core>

constexpr double radian = M_PI / 180.;
constexpr double degree = 1 / radian;

template< int O >
using Orientation = std::integral_constant< int, O >;
using Ox = Orientation< 0 >;
using Oy = Orientation< 1 >;
using Oz = Orientation< 2 >;

using Eigen::cos;
using Eigen::sin;
using Eigen::atan;
using Eigen::pow;
using Eigen::sqrt;
using Eigen::abs;

template< typename Real, typename = typename std::enable_if< std::is_scalar< Real >::value >::type >
inline auto eval( Real val ) {
    return val;
}

template< typename Derived >
inline auto eval( const Eigen::DenseBase< Derived >& val ) {
    return val.eval();
}

template< typename Real >
inline auto select( bool if_cond, Real then_val, Real else_val ) {
    return if_cond ? then_val : else_val;
}

template< typename Derived, typename Then, typename Else >
inline auto select( const Eigen::DenseBase< Derived >& if_cond, const Then& then_val, const Else& else_val ) {
    return if_cond.select( then_val, else_val );
}

template< typename Derived >
inline auto hypot( const Eigen::ArrayBase< Derived >& x, const Eigen::ArrayBase< Derived >& y ) {
    return sqrt( x * x + y * y );
}

template< typename Derived, typename Other >
inline auto hypot( const Eigen::ArrayBase< Derived >& x, const Other& y ) {
    return sqrt( x * x + y * y );
}

template< typename Derived >
inline auto atan2( const Eigen::ArrayBase< Derived >& y, const Eigen::ArrayBase< Derived >& x ) {
    auto phi = atan( y / x );
    return select( x < 0, select( y < 0, phi - M_PI, phi + M_PI ), phi );
}

template< typename Scalar >
struct Array {};

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
template< typename Scalar >
struct Point3< Array< Scalar > > {
    Eigen::Array< Scalar, Eigen::Dynamic, 3 > storage;

    static Point3< Array< Scalar > > create( Eigen::Index size ) {
        Point3< Array< Scalar > > point;
        point.storage.resize( size, 3 );
        return point;
    }

    auto x() {
        return storage.col( 0 );
    }

    auto y() {
        return storage.col( 1 );
    }

    auto z() {
        return storage.col( 2 );
    }

    template< int O >
    auto get( Orientation< O > o ) {
        return storage.col( o );
    }

    template< int O1, int O2 >
    void swap( Orientation< O1 > o1, Orientation< O2 > o2 ) {
        get( o1 ).swap( get( o2 ) );
    }

    template< int O1, int O2, typename Real >
    void rotate( Orientation< O1 > o1, Orientation< O2 > o2, Real cos, Real sin ) {
        auto v1 = get( o1 ).eval();
        auto v2 = get( o2 );
        get( o1 ) = cos * v1 + sin * v2;
        get( o2 ) = -sin * v1 + cos * v2;
    }
};

/**
 *
 */
struct Config {
    Point3< double > origin;
    double scale_factor;

    static Config create() {
        Config config;
        config.origin = Point3< double >::create( 0, 0, 0 );
        config.scale_factor = 0.75;
        return config;
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
        double longitude = config.origin.x() * radian;
        double latitude = config.origin.y() * radian;
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
        georef.origin = config.origin;
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

/**
 *
 */
int main( int argc, char** argv ) {
    auto config = Config::create();
    auto georef = Georef::create( config );
    std::cout << std::fixed << std::setprecision( 3 );

    {
        auto point = Point3< double >::create( 90, 45, 0 );
        georef.geod2ecef( point ).ecef2topo( point );
        std::cout << point.x() << " " << point.y() << " " << point.z() << "\n";
        georef.topo2ecef( point ).ecef2geod( point );
        std::cout << point.x() << " " << point.y() << " " << point.z() << "\n\n";
    }

    {
        Eigen::Index size = 1'000'000;
        auto points = Point3< Array< double > >::create( size );
        points.x().setLinSpaced( size, -180, 180 );
        points.y().setLinSpaced( size, -90, 90 );

        auto timer = std::chrono::steady_clock::now();
        georef.geod2ecef( points ).ecef2topo( points );
        georef.topo2ecef( points ).ecef2geod( points );
        std::cout << std::chrono::duration_cast< std::chrono::duration< double > >(
                std::chrono::steady_clock::now() - timer ).count() << "s\n";
    }

    return 0;
}
