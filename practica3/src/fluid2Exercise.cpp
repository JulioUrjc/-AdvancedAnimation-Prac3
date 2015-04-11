
#include "scene.h"
#include "pcg_solver.h"

namespace{
	//////////////////////////////////////////////
	// Add any custom classes or functions here! //
	//////////////////////////////////////////////
	
	inline int clamp(int a, int b, int c){	return a<b ? b:(a>c ? c:a); }

	float interpolation(Vec2 idPos, Array2<float>& inkAux){
		// Calculate min and max position
		Vec2 posMin(floorf(idPos.x), floorf(idPos.y));
		Vec2 posMax(ceilf(idPos.x), ceilf(idPos.y));
		Vec2 mid(idPos-posMin);

		Index2 min((int)posMin.x, (int)posMin.y);
		Index2 max((int)posMax.x, (int)posMax.y);
		// Calculate index cells
		Index2 s= inkAux.getSize();
		Index2 id1(clamp(min.x, 0, s.x-1), clamp(min.y, 0, s.y-1));
		Index2 id2(clamp(max.x, 0, s.x-1), clamp(min.y, 0, s.y-1));
		Index2 id3(clamp(min.x, 0, s.x-1), clamp(max.y, 0, s.y-1));
		Index2 id4(clamp(max.x, 0, s.x-1), clamp(max.y, 0, s.y-1));

		// Bilinear interpolation
		// alfa= x-i	beta= y-j
		// U_alfa_j=  (1-alfa)*u_ij + alfa*u_i + 1j
		// U_alfa_j+1= (1-alfa)*u_ij + 1 + alfa*u_i + 1j + 1
		// I= (1-beta)u_alfa_j + beta*u_alfa_j + 1
		return (inkAux[id1]*(1.0f-mid.x)+inkAux[id2]*mid.x)* (1.0f-mid.y) + (inkAux[id3]*(1.0f-mid.x)+inkAux[id4]* mid.x)* mid.y;
	}	
}

/***** Advection *****/
void Fluid2::fluidAdvection( const float dt ){
    // ink advection
		Array2<float> inkAux(ink);																// Copy for the interpolation
		Index2 size = ink.getSize();															// Ink size

		for(int i= 0; i<size.x; ++i){
			for(int j= 0; j<size.y; ++j){
				Index2 id(i, j);																// Global index

				Vec2 position(grid.getCellPos(id));												// Global coordinates
				Vec2 velocity((velocityX[id] + velocityX[Index2(i+1, j)])*0.5f, (velocityY[id] + velocityY[Index2(i, j+1)])*0.5f); // Advencion Method
				Vec2 newPosition(position - dt*velocity);										// Integrate position

				ink[id] = interpolation(grid.getCellIndex(newPosition), inkAux);				// Bilinear interpolation
			}
		}

    // velocity advection
			Array2<float> xCopy(velocityX);					 									// Copy X for the interpolation
			Array2<float> yCopy(velocityY);					 									// Copy Y for the interpolation
			
			Index2 sizeX = velocityX.getSize();				 									// X size
			Index2 sizeY = velocityY.getSize();				 									// Y size

			for(int i= 0; i<sizeX.x; ++i){
				for(int j= 0; j<sizeX.y; ++j){
					Index2 id(i, j);
					// Calculate index of the previus velocities in X to compute the velocityX
					Index2 idX1(clamp(i-1, 0, sizeY.x-1), clamp(j  , 0, sizeY.y-1));
					Index2 idX2(clamp(i  , 0, sizeY.x-1), clamp(j  , 0, sizeY.y-1));
					Index2 idX3(clamp(i-1, 0, sizeY.x-1), clamp(j+1, 0, sizeY.y-1));
					Index2 idX4(clamp(i  , 0, sizeY.x-1), clamp(j+1, 0, sizeY.y-1));

					Vec2 position(grid.getFaceXPos(id));												 // Calculate global coordinates VelocityX
					Vec2 velocity(xCopy[id], (yCopy[idX1]+yCopy[idX2]+yCopy[idX3]+yCopy[idX4])*0.25f);   // Advencion Method
					Vec2 newXPosition(position - dt*velocity);     										 // Integrate X position

					velocityX[id] = interpolation(grid.getFaceIndex(newXPosition, 0), xCopy);			 // Bilinear interpolation
				}
			}

			for(int i= 0; i<sizeY.x; ++i){
				for(int j= 0; j<sizeY.y; ++j){
					Index2 id(i, j);
					// Calculate index of the previus velocities in Y to compute the velocityY
					Index2 idY1(clamp(i  , 0, sizeX.x-1), clamp(j-1, 0, sizeX.y-1));
					Index2 idY2(clamp(i  , 0, sizeX.x-1), clamp(j  , 0, sizeX.y-1));
					Index2 idY3(clamp(i+1, 0, sizeX.x-1), clamp(j-1, 0, sizeX.y-1));
					Index2 idY4(clamp(i+1, 0, sizeX.x-1), clamp(j  , 0, sizeX.y-1));

					Vec2 position(grid.getFaceYPos(id));											   // Calculate global coordinates VelocityY
					Vec2 velocity((xCopy[idY1]+xCopy[idY2]+xCopy[idY3]+xCopy[idY4])*0.25f, yCopy[id]); // Advencion Method
					Vec2 newYPosition(position-dt*velocity);     									   // Integrate Y position

					velocityY[id] = interpolation(grid.getFaceIndex(newYPosition, 1), yCopy);		   // Bilinear interpolation
				}
			}
}

/***** Emission *****/
void Fluid2::fluidEmission(){
	if( Scene::testcase >= Scene::SMOKE ){
        // emit source ink
        {
        }
        // emit source velocity
        {
        }
	}
}

// volume forces
void Fluid2::fluidVolumeForces( const float dt ){
	if( Scene::testcase >= Scene::SMOKE ){
        // gravity
	}
}

// viscosity
void Fluid2::fluidViscosity( const float dt ){
	if( Scene::testcase >= Scene::SMOKE ){
        // viscosity
	}
}

// pressure
void Fluid2::fluidPressureProjection( const float dt ){
	if( Scene::testcase >= Scene::SMOKE ){
        // pressure
	}
}