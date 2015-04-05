
#include "scene.h"
#include "pcg_solver.h"

namespace
{
	//////////////////////////////////////////////
	// Add any custom classes or functions here! //
	//////////////////////////////////////////////
    
}

// advection
void Fluid2::fluidAdvection( const float dt )
{
    // ink advection
    {
    }

    // velocity advection
    {
    }
}

// emission
void Fluid2::fluidEmission()
{
	if( Scene::testcase >= Scene::SMOKE )
	{
        // emit source ink
        {
        }
        // emit source velocity
        {
        }
	}
}

// volume forces
void Fluid2::fluidVolumeForces( const float dt )
{
	if( Scene::testcase >= Scene::SMOKE )
	{
        // gravity
	}
}

// viscosity
void Fluid2::fluidViscosity( const float dt )
{
	if( Scene::testcase >= Scene::SMOKE )
	{
        // viscosity
	}
}

// pressure
void Fluid2::fluidPressureProjection( const float dt )
{
	if( Scene::testcase >= Scene::SMOKE )
	{
        // pressure
	}
}
