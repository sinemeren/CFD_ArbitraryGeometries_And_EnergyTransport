#include "uvp.hpp"
#include "helper.hpp"
#include "datastructures.hpp"
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cfloat>
#include "boundary_val.hpp"

// Determines the value of F and G
void calculate_fg(
        double Re,
        double GX,
        double GY,
        double alpha,
        double dt,
        double dx,
        double dy,
        int imax,
        int jmax,
        Grid& grid,
        matrix<double> &F,
        matrix<double> &G,
        double beta,
        double heatTransfer, const matrix<unsigned int> &flags)
{

    matrix<double> U,V, T;
    grid.velocity(U,velocity_type::U);
    grid.velocity(V,velocity_type::V);
    grid.temperature(T);

    /* ---------------------------------------set boundary values for F and G-------------------------------------------------------*/

    // -------------------------------------------------- Left boundary --------------------------------------------------------------
	for(size_t j = 1; j <= jmax; j++){

		if( B_N(flags.at(0).at(j)) ){ G.at(0).at(j) = V.at(0).at(j); }

		if( B_S(flags.at(0).at(j)) ){ G.at(0).at(j-1) = V.at(0).at(j-1); }

		if( B_W(flags.at(0).at(j)) ){ throw std::runtime_error(std::string("impossible to have B_W at left boundary in FG"));}

		if( B_E(flags.at(0).at(j)) ){ F.at(0).at(j) = U.at(0).at(j); }

		if( B_NE(flags.at(0).at(j)) ){ F.at(0).at(j) = U.at(0).at(j); G.at(0).at(j) = V.at(0).at(j); }

		if( B_NW(flags.at(0).at(j)) ){throw std::runtime_error(std::string("impossible to have B_NW at left boundary in FG")); }

		if( B_SE(flags.at(0).at(j)) ){ F.at(0).at(j) = U.at(0).at(j); G.at(0).at(j-1) = V.at(0).at(j-1); }

		if( B_SW(flags.at(0).at(j))){throw std::runtime_error(std::string("impossible to have B_SW at left boundary in FG")); }

		// Inflow - only x direction
		if((flags.at(0).at(j) & flag_10bit::INFLOW) == flag_10bit::INFLOW) { F.at(0).at(j) = U.at(0).at(j);}

		// Outflow - only x direction
		if((flags.at(0).at(j) & flag_10bit::OUTFLOW) == flag_10bit::OUTFLOW) { throw std::runtime_error(std::string("impossible to have outflow at left boundary in FG")); }

	}

	// ---------------------------------------------------- Right boundary -------------------------------------------------------------
	for(size_t j = 1; j <= jmax; j++){

		if( B_N(flags.at(imax+1).at(j)) ){ G.at(imax+1).at(j) = V.at(imax+1).at(j); }

		if( B_S(flags.at(imax+1).at(j)) ){ G.at(imax+1).at(j-1) = V.at(imax+1).at(j-1); }

		if( B_W(flags.at(imax+1).at(j)) ){ F.at(imax).at(j) = U.at(imax).at(j); }

		if( B_E(flags.at(imax+1).at(j)) ){ throw std::runtime_error(std::string("impossible to have B_E at right boundary in FG")); }

		if( B_NE(flags.at(imax+1).at(j)) ){ throw std::runtime_error(std::string("impossible to have B_NE at right boundary in FG")); }

		if( B_NW(flags.at(imax+1).at(j)) ){ F.at(imax).at(j) = U.at(imax).at(j); G.at(imax+1).at(j) = V.at(imax+1).at(j); }

		if( B_SE(flags.at(imax+1).at(j)) ){ throw std::runtime_error(std::string("impossible to have B_SE at right boundary in FG"));  }

		if( B_SW(flags.at(imax+1).at(j))){ F.at(imax).at(j) = U.at(imax).at(j); G.at(imax+1).at(j-1) = V.at(imax+1).at(j-1); }

		// Inflow - only x direction
		if((flags.at(imax+1).at(j) & flag_10bit::INFLOW) == flag_10bit::INFLOW) { throw std::runtime_error(std::string("impossible to have inflow at right boundary in FG")); }

		// Outflow - only x direction
		if((flags.at(imax+1).at(j) & flag_10bit::OUTFLOW) == flag_10bit::OUTFLOW) { F.at(imax+1).at(j) = U.at(imax).at(j);}

	}

	// ------------------------------------------------------- Top boundary ------------------------------------------------------
	for( size_t i = 1; i <= imax; i++){

		if( B_N(flags.at(i).at(jmax+1)) ){ throw std::runtime_error(std::string("impossible to have B_N at top boundary in FG")); }

		if( B_S(flags.at(i).at(jmax+1)) ){ G.at(i).at(jmax) = V.at(i).at(jmax); }

		if( B_W(flags.at(i).at(jmax+1)) ){ F.at(i-1).at(jmax+1) = U.at(i-1).at(jmax+1); }

		if( B_E(flags.at(i).at(jmax+1)) ){ F.at(i).at(jmax+1) = U.at(i).at(jmax+1); }

		if( B_NE(flags.at(i).at(jmax+1)) ){ throw std::runtime_error(std::string("impossible to have B_NE at top boundary in FG")); }

		if( B_NW(flags.at(i).at(jmax+1)) ){ throw std::runtime_error(std::string("impossible to have B_NW at top boundary in FG")); }

		if( B_SE(flags.at(i).at(jmax+1)) ){ F.at(i).at(jmax+1) = U.at(i).at(jmax+1); G.at(i).at(jmax) = V.at(i).at(jmax); }

		if( B_SW(flags.at(i).at(jmax+1))){ F.at(i-1).at(jmax+1) = U.at(i-1).at(jmax+1); G.at(i).at(jmax) = V.at(i).at(jmax); }

		// Inflow - only x direction
		if((flags.at(i).at(jmax+1) & flag_10bit::INFLOW) == flag_10bit::INFLOW) { F.at(i).at(jmax+1) = U.at(i).at(jmax+1);}

		// Outflow - only x direction
		if((flags.at(i).at(jmax+1) & flag_10bit::OUTFLOW) == flag_10bit::OUTFLOW) { F.at(i).at(jmax+1) = U.at(i-1).at(jmax+1);}

	}

	// ---------------------------------------------------- Bottom boundary ------------------------------------------------------------
	for( size_t i = 1; i <= imax; i++){

		if( B_N(flags.at(i).at(0)) ){ G.at(i).at(0) = V.at(i).at(0); }

		if( B_S(flags.at(i).at(0)) ){ throw std::runtime_error(std::string("impossible to have B_S at bottom boundary in FG")); }

		if( B_W(flags.at(i).at(0)) ){ F.at(i-1).at(0) = U.at(i-1).at(0); }

		if( B_E(flags.at(i).at(0)) ){ F.at(i).at(0) = U.at(i).at(0); }

		if( B_NE(flags.at(i).at(0)) ){ F.at(i).at(0) = U.at(i).at(0); G.at(i).at(0) = V.at(i).at(0); }

		if( B_NW(flags.at(i).at(0)) ){ F.at(i-1).at(0) = U.at(i-1).at(0); G.at(i).at(0) = V.at(i).at(0); }

		if( B_SE(flags.at(i).at(0)) ){ throw std::runtime_error(std::string("impossible to have B_SE at bottom boundary in FG")); }

		if( B_SW(flags.at(i).at(0))){ throw std::runtime_error(std::string("impossible to have B_WW at bottom boundary in FG")); }

		// Inflow - only x direction
		if((flags.at(i).at(0) & flag_10bit::INFLOW) == flag_10bit::INFLOW) { F.at(i).at(0) = U.at(i).at(0);}

		// Outflow - only x direction
		if((flags.at(i).at(0) & flag_10bit::OUTFLOW) == flag_10bit::OUTFLOW) { F.at(i).at(0) = U.at(i-1).at(0);}

	}

	// --------------------------------------------------- Top-left corner --------------------------------------------------------
	if( B_N(flags.at(0).at(jmax+1)) ){ throw std::runtime_error(std::string("impossible to have B_N at top-left boundary in FG")); }

	if( B_S(flags.at(0).at(jmax+1)) ){ G.at(0).at(jmax) = V.at(0).at(jmax); }

	if( B_W(flags.at(0).at(jmax+1)) ){ throw std::runtime_error(std::string("impossible to have B_W at top-left boundary in FG")); }

	if( B_E(flags.at(0).at(jmax+1)) ){ F.at(0).at(jmax+1) = U.at(0).at(jmax+1); }

	if( B_NE(flags.at(0).at(jmax+1)) ){ throw std::runtime_error(std::string("impossible to have B_NE at top-left boundary in FG")); }

	if( B_NW(flags.at(0).at(jmax+1)) ){ throw std::runtime_error(std::string("impossible to have B_NW at top-left boundary in FG")); }

	if( B_SE(flags.at(0).at(jmax+1)) ){ F.at(0).at(jmax+1) = U.at(0).at(jmax+1); G.at(0).at(jmax) = V.at(0).at(jmax); }

	if( B_SW(flags.at(0).at(jmax+1))){ throw std::runtime_error(std::string("impossible to have B_SW at top-left boundary in FG")); }

	// Inflow - only x direction
	if((flags.at(0).at(jmax+1) & flag_10bit::INFLOW) == flag_10bit::INFLOW) { F.at(0).at(jmax+1) = U.at(0).at(jmax+1);}

	// Outflow - only x direction
	if((flags.at(0).at(jmax+1) & flag_10bit::OUTFLOW) == flag_10bit::OUTFLOW) { throw std::runtime_error(std::string("impossible to have outflow at top-left boundary in FG"));}


	// ---------------------------------------------------- Top-right corner -----------------------------------------------------------
	if( B_N(flags.at(imax+1).at(jmax+1)) ){ throw std::runtime_error(std::string("impossible to have B_N at top-right boundary in FG")); }

	if( B_S(flags.at(imax+1).at(jmax+1)) ){ G.at(imax+1).at(jmax) = V.at(imax+1).at(jmax); }

	if( B_W(flags.at(imax+1).at(jmax+1)) ){ F.at(imax).at(jmax+1) = U.at(imax).at(jmax+1); }

	if( B_E(flags.at(imax+1).at(jmax+1)) ){ throw std::runtime_error(std::string("impossible to have B_W at top-right boundary in FG")); }

	if( B_NE(flags.at(imax+1).at(jmax+1)) ){ throw std::runtime_error(std::string("impossible to have B_N at top-right boundary in FG")); }

	if( B_NW(flags.at(imax+1).at(jmax+1)) ){ throw std::runtime_error(std::string("impossible to have B_NW at top-right boundary in FG")); }

	if( B_SE(flags.at(imax+1).at(jmax+1)) ){ throw std::runtime_error(std::string("impossible to have B_SE at top-right boundary in FG")); }

	if( B_SW(flags.at(imax+1).at(jmax+1))){ F.at(imax).at(jmax+1) = U.at(imax).at(jmax+1); G.at(imax+1).at(jmax) = V.at(imax+1).at(jmax); }

	// Inflow - only x direction
	if((flags.at(imax+1).at(jmax+1) & flag_10bit::INFLOW) == flag_10bit::INFLOW) { throw std::runtime_error(std::string("impossible to have inflow at top-right boundary in FG"));}

	// Outflow - only x direction
	if((flags.at(imax+1).at(jmax+1) & flag_10bit::OUTFLOW) == flag_10bit::OUTFLOW) { F.at(imax+1).at(jmax+1) = U.at(imax).at(jmax+1);}


	// ----------------------------------------------------- Bottom-left corner -----------------------------------------------------------

	if( B_N(flags.at(0).at(0)) ){ G.at(0).at(0) = V.at(0).at(0); }

	if( B_S(flags.at(0).at(0)) ){ throw std::runtime_error(std::string("impossible to have B_S at bottom-left boundary in FG")); }

	if( B_W(flags.at(0).at(0)) ){ throw std::runtime_error(std::string("impossible to have B_W at bottom-left boundary in FG")); }

	if( B_E(flags.at(0).at(0)) ){ F.at(0).at(0) = U.at(0).at(0); }

	if( B_NE(flags.at(0).at(0)) ){ F.at(0).at(0) = U.at(0).at(0); G.at(0).at(0) = V.at(0).at(0); }

	if( B_NW(flags.at(0).at(0)) ){ throw std::runtime_error(std::string("impossible to have B_NW at bottom-left boundary in FG")); }

	if( B_SE(flags.at(0).at(0)) ){ throw std::runtime_error(std::string("impossible to have B_SE at bottom-left boundary in FG")); }

	if( B_SW(flags.at(0).at(0))){ throw std::runtime_error(std::string("impossible to have B_SW at bottom-left boundary in FG")); }

	// Inflow - only x direction
	if((flags.at(0).at(0) & flag_10bit::INFLOW) == flag_10bit::INFLOW) { F.at(0).at(0) = U.at(0).at(0);}

	// Outflow - only x direction
	if((flags.at(0).at(0) & flag_10bit::OUTFLOW) == flag_10bit::OUTFLOW) { throw std::runtime_error(std::string("impossible to have outflow at bottom-left boundary in FG"));}

	// ------------------------------------------------------ Bottom-right corner ------------------------------------------------------------

	if( B_N(flags.at(imax+1).at(0)) ){ G.at(imax+1).at(0) = V.at(imax+1).at(0); }

	if( B_S(flags.at(imax+1).at(0)) ){ throw std::runtime_error(std::string("impossible to have B_S at bottom-right boundary in FG")); }

	if( B_W(flags.at(imax+1).at(0)) ){ F.at(imax).at(0) = U.at(imax).at(0); }

	if( B_E(flags.at(imax+1).at(0)) ){ throw std::runtime_error(std::string("impossible to have B_E at bottom-right boundary in FG")); }

	if( B_NE(flags.at(imax+1).at(0)) ){ throw std::runtime_error(std::string("impossible to have B_NE at bottom-right boundary in FG")); }

	if( B_NW(flags.at(imax+1).at(0)) ){ F.at(imax).at(0) = U.at(imax).at(0); G.at(imax+1).at(0) = V.at(imax+1).at(0); }

	if( B_SE(flags.at(imax+1).at(0)) ){ throw std::runtime_error(std::string("impossible to have B_SE at bottom-right boundary in FG")); }

	if( B_SW(flags.at(imax+1).at(0))){ throw std::runtime_error(std::string("impossible to have B_SW at bottom-right boundary in FG")); }

	// Inflow - only x direction
	if((flags.at(imax+1).at(0) & flag_10bit::INFLOW) == flag_10bit::INFLOW) { throw std::runtime_error(std::string("impossible to have inflow at bottom-right boundary in FG"));}

	// Outflow - only x direction
	if((flags.at(imax+1).at(0) & flag_10bit::OUTFLOW) == flag_10bit::OUTFLOW) { F.at(imax+1).at(0) = U.at(imax).at(0);}


    // ----------------------------------------- Inner cells -------------------------------------------------------------------
    for( size_t i = 1; i <= imax; i++){

        for(size_t j = 1; j <= jmax; j++){


            if( B_N(flags.at(i).at(j)) ){ G.at(i).at(j) = V.at(i).at(j); }

            if( B_S(flags.at(i).at(j)) ){ G.at(i).at(j-1) = V.at(i).at(j-1); }

            if( B_W(flags.at(i).at(j)) ){ F.at(i-1).at(j) = U.at(i-1).at(j); }

            if( B_E(flags.at(i).at(j)) ){ F.at(i).at(j) = U.at(i).at(j); }

            if( B_NE(flags.at(i).at(j)) ){ F.at(i).at(j) = U.at(i).at(j); G.at(i).at(j) = V.at(i).at(j); }

            if( B_NW(flags.at(i).at(j)) ){ F.at(i-1).at(j) = U.at(i-1).at(j); G.at(i).at(j) = V.at(i).at(j); }

            if( B_SE(flags.at(i).at(j)) ){ F.at(i).at(j) = U.at(i).at(j); G.at(i).at(j-1) = V.at(i).at(j-1); }

            if( B_SW(flags.at(i).at(j))){ F.at(i-1).at(j) = U.at(i-1).at(j); G.at(i).at(j-1) = V.at(i).at(j-1); }

            // Inflow - only x direction
            if((flags.at(i).at(j) & flag_10bit::INFLOW) == flag_10bit::INFLOW) { F.at(i).at(j) = U.at(i).at(j);}

            // Outflow - only x direction
            if((flags.at(i).at(j) & flag_10bit::OUTFLOW) == flag_10bit::OUTFLOW) { F.at(i).at(j) = U.at(i-1).at(j);}

        }
    }



    int lmax,kmax;
    if (heatTransfer == 1.0){
    	lmax = imax - 1;
    	kmax = jmax - 1;

    }
    else{
    	lmax = imax;
    	kmax = jmax;
    }

    // ------------------------------------------------------ F and G calculations -------------------------------------------------------------
    for( size_t i = 1; i <= lmax; i++){

        for( size_t j = 1; j <= jmax; j++){

            if(((flags.at(i).at(j) & flag_10bit::CELL_TYPE) == flag_10bit::CELL_TYPE)){

                F.at(i).at(j) = U.at(i).at(j) + dt *(

                        (1.0 / Re) *( ( ( U.at(i-1).at(j) - 2.0 * U.at(i).at(j) + U.at(i+1).at(j) ) / ( dx * dx ) ) +  ( ( U.at(i).at(j-1) - 2.0 * U.at(i).at(j) + U.at(i).at(j+1) ) / ( dy * dy ) ) )

                        // d(u^2)/d(x)
                        - (0.25) * ( 1.0 / dx ) * (  pow( ( U.at(i).at(j) + U.at(i+1).at(j) ), 2.0 )  -   pow( ( U.at(i-1).at(j) + U.at(i).at(j) ), 2.0 )
                                                     + alpha * ( ( fabs( U.at(i).at(j) + U.at(i+1).at(j) )  * ( U.at(i).at(j) - U.at(i+1).at(j)) ) -   (  fabs( U.at(i-1).at(j) + U.at(i).at(j)) * ( U.at(i-1).at(j) - U.at(i).at(j) ) ) )
                        )

                        // d(uv)/d(y)
                        - (0.25) * ( 1.0 / dy ) * ( ( V.at(i).at(j) + V.at(i+1).at(j) ) * ( U.at(i).at(j) + U.at(i).at(j+1) ) - ( V.at(i).at(j-1) + V.at(i+1).at(j-1)) * ( U.at(i).at(j-1) + U.at(i).at(j) )
                                                    + alpha * ( fabs(V.at(i).at(j) + V.at(i+1).at(j)) * ( U.at(i).at(j) - U.at(i).at(j+1) ) - fabs( V.at(i).at(j-1) + V.at(i+1).at(j-1)) * ( U.at(i).at(j-1) - U.at(i).at(j) )   )
                        )

                        + GX)

                                -  heatTransfer * GX * beta * dt *0.5 * ( T.at(i).at(j) + T.at(i+1).at(j));
            }
        }
    }



    for(size_t i = 1; i <= imax; i++){

        for( size_t j = 1; j <= kmax; j++){

            if(((flags.at(i).at(j) & flag_10bit::CELL_TYPE) == flag_10bit::CELL_TYPE))
            {
                G.at(i).at(j) = V.at(i).at(j) + dt * (
                        (1.0 / Re) * (((V.at(i - 1).at(j) - 2.0 * V.at(i).at(j) + V.at(i + 1).at(j)) / (dx * dx)) +
                                      ((V.at(i).at(j - 1) - 2.0 * V.at(i).at(j) + V.at(i).at(j + 1)) / (dy * dy)))

                        // d(v^2) / dy
                        - (0.25) * (1.0 / dy) *
                          (pow((V.at(i).at(j) + V.at(i).at(j + 1)), 2.0) - pow((V.at(i).at(j - 1) + V.at(i).at(j)), 2.0)
                           + alpha * ((fabs(V.at(i).at(j) + V.at(i).at(j + 1)) * (V.at(i).at(j) - V.at(i).at(j + 1))) -
                                      (fabs(V.at(i).at(j - 1) + V.at(i).at(j)) * (V.at(i).at(j - 1) - V.at(i).at(j))))
                          )

                        // d(uv)/dx
                        - (0.25) * (1.0 / dx) *
                          ((U.at(i).at(j) + U.at(i).at(j + 1)) * (V.at(i).at(j) + V.at(i + 1).at(j)) -
                           (U.at(i - 1).at(j) + U.at(i - 1).at(j + 1)) * (V.at(i - 1).at(j) + V.at(i).at(j))
                           + alpha * (fabs(U.at(i).at(j) + U.at(i).at(j + 1)) * (V.at(i).at(j) - V.at(i + 1).at(j)) -
                                      fabs(U.at(i - 1).at(j) + U.at(i - 1).at(j + 1)) *
                                      (V.at(i - 1).at(j) - V.at(i).at(j)))
                          )

                        + GY)
                                - heatTransfer * GY * beta * dt * 0.5 * (T.at(i).at(j) + T.at(i).at(j + 1));

            }
        }
    }

}

// This operation computes the right hand side of the pressure poisson equation.
void calculate_rs(
        double dt,
        double dx,
        double dy,
        int imax,
        int jmax,
        matrix<double> &F,
        matrix<double> &G,
        matrix<double> &RS, const matrix<unsigned int> &flags)
{


    for(size_t i = 1; i <= imax ; i++){

        for( size_t j = 1; j <= jmax ; j++)
        {
            if(((flags.at(i).at(j) & flag_10bit::CELL_TYPE) == flag_10bit::CELL_TYPE))
            {
                RS.at(i).at(j) = ( 1.0 / dt ) * ( ( F.at(i).at(j) - F.at(i-1).at(j) ) / dx  + (G.at(i).at(j) - G.at(i).at(j-1)) / dy );
            }
        }
    }
}

// Determines the maximal time step size
void calculate_dt(double Re, double Pr,  double tau, double *dt, double dx, double dy, int imax, int jmax, Grid &grid, double heatTransfer)
{
    matrix<double> U,V;
    grid.velocity(U,velocity_type::U);
    grid.velocity(V,velocity_type::V);

    double u_max = maxElementOfMatrix(U);
    double v_max = maxElementOfMatrix(V);

    double term = ( 1.0 / ((1.0 / (dx * dx) ) + ( 1.0 / (dy * dy) ) ) ) ;
    double term1 = (Re / 2.0) * term;
    double term2 = dx / fabs(u_max);
    double term3 = dy / fabs(v_max);


    double min1 = std::min(term1,term2);
    double min2 = std::min(min1,term3);

    if( heatTransfer == 1.0){

        min2 = std::min(min2, ( 0.5 * Re * Pr * term ) );
    }

    //*dt = tau * min2;
    if( (tau > 0) && (tau < 1))
    {
        *dt = tau * min2;
    }

}

void calculate_uv(
        double dt,
        double dx,
        double dy,
        int imax,
        int jmax,
        Grid& grid,
        matrix<double> &F,
        matrix<double> &G,
		double heatTransfer,
        const matrix<unsigned int> &flags)
{
    matrix<double> U,V,P;
    grid.velocity(U,velocity_type::U);
    grid.velocity(V,velocity_type::V);
    grid.pressure(P);

    int lmax,kmax;
	if (heatTransfer == 1.0){
		lmax = imax - 1;
		kmax = jmax - 1;

	}
	else{
		lmax = imax;
		kmax = jmax;
	}

    for( size_t i = 1; i <= lmax; i++){

        for( size_t j = 1; j <= jmax ; j++){

            if(((flags.at(i).at(j) & flag_10bit::CELL_TYPE) == flag_10bit::CELL_TYPE) ){
                U.at(i).at(j) = F.at(i).at(j) - (dt / dx) * ( P.at(i+1).at(j) - P.at(i).at(j));

            }
        }
    }

    for( size_t i = 1; i <= imax; i++){

        for( size_t j = 1; j <= kmax; j++) {

            if(((flags.at(i).at(j) & flag_10bit::CELL_TYPE) == flag_10bit::CELL_TYPE) ){

                V.at(i).at(j) = G.at(i).at(j) - (dt / dy) * (P.at(i).at(j+1) - P.at(i).at(j));

            }
        }
    }

    grid.set_velocity(U,velocity_type::U);
    grid.set_velocity(V,velocity_type::V);

}


void calculate_temp(double dt, double dx, double dy, int imax, int jmax, double Re, double Pr, double alpha, Grid &grid, const matrix<unsigned int> &flags){


    matrix<double> T, T_n, U, V;


    grid.temperature(T);
    grid.temperature(T_n);
    grid.velocity(U, velocity_type::U);
    grid.velocity(V, velocity_type::V);

    double diffusion_x, diffusion_y, convection1, convection2;

    for(size_t i = 1; i <= imax; i++){

        for( size_t j = 1; j <= jmax; j++){


            if(((flags.at(i).at(j) & flag_10bit::CELL_TYPE) == flag_10bit::CELL_TYPE) ){

                // diffusion terms
                diffusion_x = (T.at(i+1).at(j) - 2.0 * T.at(i).at(j) + T.at(i-1).at(j)) / (dx * dx);
                diffusion_y = (T.at(i).at(j+1) - 2.0 * T.at(i).at(j) + T.at(i).at(j-1)) / (dy * dy);

                // convection terms
                convection1 = (1.0 / dx) * ((U.at(i).at(j) * ( 0.5 *  (T.at(i).at(j) + T.at(i+1).at(j))  )) - (U.at(i-1).at(j) * ( 0.5 *  (T.at(i-1).at(j) + T.at(i).at(j))  ) ))
                              + (alpha / dx) * ( fabs(U.at(i).at(j)) * ( 0.5 *  (T.at(i).at(j)- T.at(i+1).at(j))  ) - fabs(U.at(i-1).at(j)) * ( 0.5 *  (T.at(i-1).at(j) - T.at(i).at(j))  ));

                convection2 = (1.0 / dy) * ((V.at(i).at(j) * ( 0.5 *  (T.at(i).at(j) + T.at(i).at(j+1))  )) - (V.at(i).at(j-1) * ( 0.5 *  (T.at(i).at(j-1) + T.at(i).at(j))  ) ))
                              + (alpha/ dy) * ( fabs(V.at(i).at(j)) * ( 0.5 *  (T.at(i).at(j) - T.at(i).at(j+1))  ) - fabs(V.at(i).at(j-1)) * ( 0.5 *  (T.at(i).at(j-1) - T.at(i).at(j))  ));

                T_n.at(i).at(j) = T.at(i).at(j) + dt * ( (1.0 / (Re * Pr)) * (diffusion_x + diffusion_y) - (convection1 + convection2) );

            }

        }
    }

    grid.set_temperature(T_n);

}
