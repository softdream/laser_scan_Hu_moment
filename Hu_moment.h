#ifndef __HU_MOMENT_H
#define __HU_MOMENT_H

#include "dataContainer.h"

#include <cmath>

namespace moment
{

template<typename T>
class HuMoment
{
public:
	using DataType = T;
	using ScanContainerType = typename sensor::DataContainer<T>;
	using PointType = typename Eigen::Matrix<T, 2, 1>;	


	HuMoment()
	{

	}

	~HuMoment()
	{

	}	

        const Eigen::Matrix<DataType, 7, 1> makeHuMoment( const ScanContainerType &scan )
        {
                getCentroid( scan );
                getCentralMomentByScan( scan );
                getNormalizedCentralMoment();

                DataType I1 = y20 + y02;
                DataType I2 = std::pow( y20 - y02, 2 ) + 4 * y11 * y11;
                DataType I3 = std::pow( y30 - 3 * y12, 2 ) + std::pow( 3 * y21 - y03, 2 );
                DataType I4 = std::pow( y30 - y12, 2 ) + std::pow( y21 + y03, 2 );
                DataType I5 = ( y30 - 3 * y21 ) * ( y30 + y12 ) * ( std::pow( y30 + y12, 2 ) - 3 * std::pow( y21 + y03, 2 ) )
                              + ( 3 * y21 - y03 ) * ( y21 + y03 ) * ( 3 * std::pow( y30 + y12, 2 ) - std::pow( y21 + y03, 2 ) );

                DataType I6 = ( y20 - y02 ) * ( std::pow( y30 + y12, 2 ) - std::pow( y21 + y03, 2 ) ) + 4 * y11 * ( y30 + y12 ) * ( y21 + y03 );
                DataType I7 = ( 3 * y21 - y03 ) * ( y30 + y12 ) * ( std::pow( y30 + y12, 2 ) - 3 * std::pow( y21 + y03, 2 ) ) - ( y30 - 3 * y12 ) * ( y21 + y03 ) * ( 3 * std::pow( y30 + y12, 2 ) - std::pow( y21 + y03, 2 ) );

                Eigen::Matrix<DataType, 7, 1> invarint_moment;
                invarint_moment << I1, I2, I3, I4, I5, I6, I7;

                return invarint_moment;
        }

private:
	void getCentroid( const ScanContainerType &scan )	
	{
		DataType m00 = 0;
		DataType m10 = 0;
		DataType m01 = 0;
		for( int i = 0; i < scan.getSize(); i ++ ){
			PointType point = scan[i];
			Eigen::Vector2i point_in_map = observedPointWorld2Map( point );			

			m00 += 1;	
			m10 += point_in_map[0];
			m01 += point_in_map[1];
		}
		
		x0 = m10 / m00;
		y0 = m01 / m00;
	}

	void getCentralMomentByScan( const ScanContainerType &scan )
	{

		for( int i = 0; i < scan.getSize(); i ++ ){
                        PointType point = scan[i];
			Eigen::Vector2i point_in_map = observedPointWorld2Map( point );			

			//mu20 += std::pow( point_in_map[0] - x0, 2 ) * std::pow( point_in_map[1] - y0, 0 );
			mu20 += std::pow( point_in_map[0] - x0, 2 );	
			//mu02 += std::pow( point_in_map[0] - x0, 0 ) * std::pow( point_in_map[1] - y0, 2 );
			mu02 += std::pow( point_in_map[1] - y0, 2 );
			//mu30 += std::pow( point_in_map[0] - x0, 3 ) * std::pow( point_in_map[1] - y0, 0 );
			mu30 += std::pow( point_in_map[0] - x0, 3 );
			//mu12 += std::pow( point_in_map[0] - x0, 1 ) * std::pow( point_in_map[1] - y0, 2 );
			mu12 += ( point_in_map[0] - x0) * std::pow( point_in_map[1] - y0, 2 );
			//mu21 += std::pow( point_in_map[0] - x0, 2 ) * std::pow( point_in_map[1] - y0, 1 );
			mu21 += std::pow( point_in_map[0] - x0, 2 ) * ( point_in_map[1] - y0 );
			//mu03 += std::pow( point_in_map[0] - x0, 0 ) * std::pow( point_in_map[1] - y0, 3 );
			mu03 += std::pow( point_in_map[1] - y0, 3 );
			//mu11 += std::pow( point_in_map[0] - x0, 1 ) * std::pow( point_in_map[1] - y0, 1 );
			mu11 += ( point_in_map[0] - x0 ) * ( point_in_map[1] - y0 );

			mu00 += 1;
		}

	}
	
	void getNormalizedCentralMoment(  )
	{
		//DataType y_p_q = mu_p_q / std::pow( m00, ( p + q + 2 ) * 0.5 );
		y20 = mu20 / std::pow( mu00, 2 );
		y02 = mu02 / std::pow( mu00, 2 );	
		y30 = mu30 / std::pow( mu30, 2.5 );
		y12 = mu12 / std::pow( mu12, 2.5 );
		y21 = mu21 / std::pow( mu21, 2.5 );	
		y03 = mu03 / std::pow( mu03, 2.5 );
	}
	

private:
	const Eigen::Vector2i observedPointWorld2Map( const PointType &point_in_world )
	{
		PointType tmp( point_in_world * map_scale_ );
		Eigen::Vector2i tmp_int( std::round( tmp[0] ), std::round( tmp[1] ) );	
		
		return ( tmp_int + map_center_ );
	}

private:
	DataType laser_max_dist_ = 10.0;
	DataType map_scale_ = 10.0;
	Eigen::Vector2i map_center_ = Eigen::Vector2i( 0, 0 );

	DataType x0 = 0;
	DataType y0 = 0;
	
	DataType mu00 = 0;
	
	DataType mu20 = 0;
	DataType mu02 = 0;
	DataType mu11 = 0;	
	DataType mu30 = 0;
	DataType mu12 = 0;
	DataType mu21 = 0;
	DataType mu03 = 0;

	DataType y20 = 0;
	DataType y02 = 0;
	DataType y11 = 0;
	DataType y30 = 0;
	DataType y12 = 0;
	DataType y21 = 0;
	DataType y03 = 0;
};

}

#endif
