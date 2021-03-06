#include "laserSimulation.h"

#include "Hu_moment.h"

void laserData2Container( const sensor::LaserScan &scan, sensor::ScanContainer &container )
{
        size_t size = 1440;

        float angle = -3.14159f;
        container.clear();

        for( int i = 0; i < size; i ++ ){
                float dist = scan.ranges[ i ];

                if( dist >= 0.0099999998f && dist <= 14.0000000000f ){
                        //dist *= scaleToMap;
                        container.addData( Eigen::Vector2f( cos(angle) * dist, sin(angle) * dist ) );
                }

                angle += 0.0043633231f;
        }

        std::cout<<"Scan Container Size: "<<container.getSize()<<std::endl;
}


int main()
{
	simulation::Simulation simulation;
        simulation.openSimulationFile( "frame1.txt" );

	sensor::LaserScan scan;
        sensor::ScanContainer scan_container;

        simulation.readAFrameData( scan );
        laserData2Container( scan, scan_container );

	return 0;
}
