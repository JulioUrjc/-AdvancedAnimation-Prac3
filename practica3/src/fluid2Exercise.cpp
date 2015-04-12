
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
		Index2 sizeInk = ink.getSize();															// Ink size

		for (int i= 0; i<sizeInk.x; ++i){
			for (int j= 0; j<sizeInk.y; ++j){
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
			
			Index2 sizeX= velocityX.getSize();				 									// X size
			Index2 sizeY= velocityY.getSize();				 									// Y size

			for(int i= 0; i<sizeX.x; ++i){
				for(int j= 0; j<sizeX.y; ++j){
					Index2 id(i, j);																	 // Global index
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
					Index2 id(i, j);																	 // Global index
					// Calculate index of the previus velocities in Y to compute the velocityY
					Index2 idY1(clamp(i  , 0, sizeX.x-1), clamp(j-1, 0, sizeX.y-1));
					Index2 idY2(clamp(i  , 0, sizeX.x-1), clamp(j  , 0, sizeX.y-1));
					Index2 idY3(clamp(i+1, 0, sizeX.x-1), clamp(j-1, 0, sizeX.y-1));
					Index2 idY4(clamp(i+1, 0, sizeX.x-1), clamp(j  , 0, sizeX.y-1));

					Vec2 position(grid.getFaceYPos(id));												// Calculate global coordinates VelocityY
					Vec2 velocity((xCopy[idY1]+xCopy[idY2]+xCopy[idY3]+xCopy[idY4])*0.25f, yCopy[id]);	// Advencion Method
					Vec2 newYPosition(position-dt*velocity);     										// Integrate Y position

					velocityY[id] = interpolation(grid.getFaceIndex(newYPosition, 1), yCopy);			// Bilinear interpolation
				}
			}
}

/***** Emission *****/
void Fluid2::fluidEmission(){
	if( Scene::testcase >= Scene::SMOKE ){
		Bbox2 source(-0.18f, -1.9f, 0.18f, -1.7f);
		Vec2 ismin = grid.getCellIndex(source.minPosition);
		Vec2 ismax = grid.getCellIndex(source.maxPosition);
		Index2 cornermin((int)floor(ismin.x), (int)floor(ismin.y));
		Index2 cornermax((int)ceil(ismax.x), (int)ceil(ismax.y));
		// emit source ink
		for(int i= cornermin.x; i<= cornermax.x; ++i)
		for(int j= cornermin.y; j<= cornermax.y; ++j)
			ink[Index2(i, j)]= 1.0f;
		// emit source velocity
		for(int i= cornermin.x; i<= cornermax.x+1; ++i)
		for(int j= cornermin.y; j<= cornermax.y; ++j)
			velocityX[Index2(i, j)]= 0.0f;

		for(int i= cornermin.x; i<= cornermax.x; ++i)
		for(int j= cornermin.y; j<= cornermax.y+1; ++j)
			velocityY[Index2(i, j)]= 10.0f;
	}
}

/***** Volume Forces *****/
void Fluid2::fluidVolumeForces( const float dt ){
	if( Scene::testcase >= Scene::SMOKE ){
		float incrementGravity = dt*Scene::kGravity;			// Gravity Increment
		Index2 sizeY = velocityY.getSize();						// Y size

		for (int i = 0; i < sizeY.x; ++i){
			for (int j = 0; j < sizeY.y; ++j){
				Index2 id(i, j);
				velocityY[id] += incrementGravity;				// Add increment
			}
		}
	}
}


/***** Viscosity *****/
void Fluid2::fluidViscosity( const float dt ){
	if( Scene::testcase >= Scene::SMOKE ){
		Array2<float> xCopy(velocityX);														// Copy X for the interpolation
		Array2<float> yCopy(velocityY);					 									// Copy Y for the interpolation

		Index2 sizeX= velocityX.getSize();				 									// X size
		Index2 sizeY= velocityY.getSize();				 									// Y size

		float dx= 1.0f/pow(grid.getCellDx().x, 2);
		float dy= 1.0f/pow(grid.getCellDx().y, 2);
		float cte= dt*Scene::kViscosity / Scene::kDensity;

		for(int i= 0; i<sizeX.x; ++i){
			for(int j= 0; j<sizeX.y; ++j){
				Index2 id(i, j);
				// Calculate index of velocities
				Index2 idX1(clamp(i+1, 0, sizeX.x-1), j);
				Index2 idX2(clamp(i-1, 0, sizeX.x-1), j);
				Index2 idX3(i, clamp(j+1, 0, sizeX.y-1));
				Index2 idX4(i, clamp(j-1, 0, sizeX.y-1));
				// Calculate X velocity component
				velocityX[id] += cte*((xCopy[idX1]-2.0f*xCopy[id]+xCopy[idX2])*dx+(xCopy[idX3]-2.0f*xCopy[id]+xCopy[idX4])*dy);
			}
		}

		for(int i= 0; i<sizeY.x; ++i){
			for(int j= 0; j<sizeY.y; ++j){
				Index2 id(i, j);
				// Calculate index of velocities
				Index2 idY1(clamp(i+1, 0, sizeY.x-1), j);
				Index2 idY2(clamp(i-1, 0, sizeY.x-1), j);
				Index2 idY3(i, clamp(j+1, 0, sizeY.y-1));
				Index2 idY4(i, clamp(j-1, 0, sizeY.y-1));
				// Calculate Y velocity component
				velocityY[id] += cte*((yCopy[idY1]-2.0f*yCopy[id]+yCopy[idY2])*dx+(yCopy[idY3]-2.0f*yCopy[id]+yCopy[idY4])*dy);
			}
		}
	}
}

/***** Pressure *****/
void Fluid2::fluidPressureProjection( const float dt ){
	if( Scene::testcase >= Scene::SMOKE ){

		float invDx   = 1.0f/grid.getCellDx().x;
		float invDy   = 1.0f/grid.getCellDx().y;
		float invsqrDx= 1.0f/pow(grid.getCellDx().x, 2);
		float invsqrDy= 1.0f/pow(grid.getCellDx().y, 2);

		Index2 pSize= pressure.getSize();
		Index2 sizeX= velocityX.getSize();
		Index2 sizeY= velocityY.getSize();

		// Boundaries
		for(int j= 0; j<sizeX.y; ++j){
			velocityX[Index2(0, j)]= 0.0f;
			velocityX[Index2(sizeX.x-1, j)]= 0.0f;
		}
		for(int i= 0; i<sizeY.x; ++i){
			velocityY[Index2(i, 0)] = 0.0f;
			// velocityY[ Index2(i, sizeY.y-1)]= 0.0f;				// Top solid
		}

		// rhs
		float cte = Scene::kDensity / dt;
		std::vector<double> rhs(pSize.x*pSize.y);
		for(int i= 0; i< pSize.x; ++i){
			for(int j= 0; j< pSize.y; ++j){
				Index2 id(i, j);
				rhs[pressure.getLinearIndex(i, j)]= -cte*((velocityX[Index2(i+1, j)]-velocityX[id])*invDx+(velocityY[Index2(i, j+1)]-velocityY[id])*invDy);
			}
		}

		// Matrix A
		SparseMatrix<double> A(pSize.x*pSize.y, 5);
		for (int i= 0; i< pSize.x; ++i){
			for (int j= 0; j< pSize.y; ++j){
				int id = pressure.getLinearIndex(i, j);
				if (i> 0){
					int id2 = pressure.getLinearIndex(i-1, j);
					A.add_to_element(id, id, invsqrDx);
					A.add_to_element(id, id2, -1.0*invsqrDx);
				}
				if (i< pSize.x-1) {
					int id2 = pressure.getLinearIndex(i+1, j);
					A.add_to_element(id, id, invsqrDx);
					A.add_to_element(id, id2, -1.0*invsqrDx);
				}
				if (j> 0) {
					int id2 = pressure.getLinearIndex(i, j-1);
					A.add_to_element(id, id, invsqrDy);
					A.add_to_element(id, id2, -1.0*invsqrDy);
				}
				// Top air
				A.add_to_element(id, id, invsqrDy);
				if (j< pSize.y-1){
					int id2 = pressure.getLinearIndex(i, j+1);
					A.add_to_element(id, id2, -1.0*invsqrDy);
				}
				
				//if( j< pSize.y-1 ) {												// Top solid
				//    int id2= pressure.getLinearIndex( i, j+1);
				//    A.add_to_element(id, id, invsqrDy);
				//	A.add_to_element(id, id2, -1.0*invsqrDy);
				//}
			}
		}

		// pcg solver
		PCGSolver<double> solver;
		solver.set_solver_parameters(1e-6, 10000);

		double residual_out;
		int iterations_out;
		std::vector<double> p(pSize.x*pSize.y);
		solver.solve(A, rhs, p, residual_out, iterations_out);
		//std::cout << "Pressure system result: res=" << residual_out << ", iter=" << iterations_out << std::endl;
		
		// set pressure
		for(int i = 0, n= pSize.x*pSize.y; i< n; ++i)
			pressure[i] = (float)p[i];

		// apply pressure gradient
		float K= dt/Scene::kDensity;

		for(int i = 1; i< sizeX.x-1; ++i)
		for(int j = 0; j< sizeX.y; ++j){
			Index2 id(i, j);
			float gradpressure= (pressure[id]-pressure[Index2(i-1, j)])*invDx;
			velocityX[id] -= K*gradpressure;
		}
		for(int i = 0; i< sizeY.x; ++i)
		for(int j = 1; j< sizeY.y-1; ++j){
			Index2 id(i, j);
			float gradpressure= (pressure[id]-pressure[Index2(i, j-1)])*invDy;
			velocityY[id] -= K*gradpressure;
		}
		// apply pressure gradient: Top air
		for (int i= 0; i< sizeY.x; ++i){
			Index2 id(i, sizeY.y - 1);
			float gradpressure= (0.0f-pressure[Index2(i, pSize.y-1)])*invDy;
			velocityY[id] -= K*gradpressure;
		}
	}
}