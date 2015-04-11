
#include "scene.h"
#include "pcg_solver.h"

namespace
{
	//////////////////////////////////////////////
	// Add any custom classes or functions here! //
	//////////////////////////////////////////////

	inline int clamp(int x, int a, int b){
		return (x<a) ? a : ((x>b) ? b : x);
	}

	inline float bilerp(Vec2 pos, Array2<float> field, Index2 size){
		const Vec2 posMin(floorf(pos.x), floorf(pos.y));
		const Vec2 posMax(ceilf(pos.x), ceilf(pos.y));
		const Vec2 t(pos - posMin);

		const Index2 id1(clamp((int)posMin.x, 0, size.x - 1), clamp((int)posMin.y, 0, size.y - 1));
		const Index2 id2(clamp((int)posMax.x, 0, size.x - 1), clamp((int)posMin.y, 0, size.y - 1));
		const Index2 id3(clamp((int)posMin.x, 0, size.x - 1), clamp((int)posMax.y, 0, size.y - 1));
		const Index2 id4(clamp((int)posMax.x, 0, size.x - 1), clamp((int)posMax.y, 0, size.y - 1));

		///// Bilineal interpolation - General case
		/// alfa = x - i
		/// beta = y - j
		/// U_alfa_j =  (1 - alfa) * u_ij + alfa * u_i+1j
		/// U_alfa_j+1 = (1 - alfa) * u_ij+1 + alfa * u_i+1j+1
		/// I = (1 - beta)u_alfa_j + beta * u_alfa_j+1

		const float y1 = (1.0f - t.x) * field[id1] + t.x * field[id2];
		const float y2 = (1.0f - t.x) * field[id3] + t.x * field[id4];
		const float value = (1 - t.y) * y1 + t.y * y2;
		return value;

	}

    
}

// advection
void Fluid2::fluidAdvection( const float dt )
{
    // ink advection
    {
		Array2<float> inkAux(ink);
		for (unsigned int i = 0; i < ink.getSize().x; ++i){
			for (unsigned int j = 0; j < ink.getSize().y; ++j){
				const Index2 id(i, j);

				/// Calculate global coordinates
				const Vec2 pos(grid.getCellPos(id));
				/// Calculate index of the velocitis
				const Index2 idvXaux(clamp(i + 1, 0, velocityX.getSize().x - 1), clamp(j, 0, velocityX.getSize().y - 1));
				const Index2 idvYaux(clamp(i, 0, velocityY.getSize().x - 1), clamp(j+1, 0, velocityY.getSize().y - 1));
				/// Calculate advencion method
				const Vec2 vel((velocityX[id] + velocityX[idvXaux]) * 0.5f, (velocityY[id] + velocityY[idvYaux]) * 0.5f);
				const Vec2 endpos(pos - dt * vel);
				/// Bilineal interpolation with index
				ink[id] = bilerp(grid.getCellIndex(endpos),inkAux,inkAux.getSize());;
			}
		}
    }
	
    // velocity advection
    {
		Array2< float > vXaux(velocityX);
		Array2< float > vYaux(velocityY);
		const Index2 sizeVx = velocityX.getSize();
		const Index2 sizeVy = velocityY.getSize();
		for (unsigned int i = 0; i < sizeVx.x; ++i){
			for (unsigned int j = 0; j < sizeVx.y; ++j){
				const Index2 id(i, j);
				/// Calculate global coordinates VelocityX
				const Vec2 pos(grid.getFaceXPos(id));
				/// Calculate index of the previus velocities in X to compute the velocityX
				const Index2 id1(clamp(i-1, 0, sizeVy.x-1), clamp(j,   0, sizeVy.y-1));
				const Index2 id2(clamp(i,   0, sizeVy.x-1), clamp(j,   0, sizeVy.y-1));
				const Index2 id3(clamp(i-1, 0, sizeVy.x-1), clamp(j+1, 0, sizeVy.y-1));
				const Index2 id4(clamp(i,   0, sizeVy.x-1), clamp(j+1, 0, sizeVy.y-1));
				/// Compute advection methond
				const Vec2 vel(vXaux[id], (vYaux[id1] + vYaux[id2] + vYaux[id3] + vYaux[id4])*0.25f);
				const Vec2 endpos(pos - dt*vel);
				/// Bilineal interpolation with index, face index in axis 0, horizontal
				velocityX[id] = bilerp(grid.getFaceIndex(pos, 0),vXaux,vXaux.getSize());
			}
		}
		for (unsigned int i = 0; i < sizeVy.x; ++i){
			for (unsigned int j = 0; j < sizeVy.y; ++j){
				const Index2 id(i, j);
				/// Calculate global coordinates VelocityY
				const Vec2 pos(grid.getFaceYPos(id));
				/// Calculate index of the previus velocities in X to compute the velocityX
				const Index2 id1(clamp(i,   0, sizeVx.x - 1), clamp(j-1, 0, sizeVx.y - 1));
				const Index2 id2(clamp(i,   0, sizeVx.x - 1), clamp(j-1, 0, sizeVx.y - 1));
				const Index2 id3(clamp(i+1, 0, sizeVx.x - 1), clamp(j,   0, sizeVx.y - 1));
				const Index2 id4(clamp(i+1, 0, sizeVx.x - 1), clamp(j,   0, sizeVx.y - 1));
				/// Compute advection methond
				const Vec2 vel(vXaux[id], (vYaux[id1] + vYaux[id2] + vYaux[id3] + vYaux[id4])*0.25f);
				const Vec2 endpos(pos - dt*vel);
				/// Bilineal interpolation with index, face index in axis 0, horizontal
				velocityY[id] = bilerp(grid.getFaceIndex(pos, 1), vXaux, vXaux.getSize());
			}
		}
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
