#include "boundary_val.hpp"
#include "datastructures.hpp"
#include "grid.hpp"

void boundaryvaluesTemperature(int imax, int jmax, Grid& grid, const matrix<unsigned int> &flag, double T_h, double T_c){

    matrix<double> T;
    grid.temperature(T);

    // Top boundary
    for (int i = 1; i <= imax; ++i) {

		switch(flag.at(i).at(jmax+1) & (flag_10bit::NO_SLIP | flag_10bit::FREE_SLIP)){
		case flag_10bit::NO_SLIP: // Obstacle conditions are adiabatic
		case flag_10bit::FREE_SLIP:

		if (B_E(flag[i][jmax+1]) == 1) {T[i][jmax+1] = T[i + 1][jmax+1];}

		if (B_W(flag[i][jmax+1]) == 1) {T[i][jmax+1] = T[i - 1][jmax+1];}

		if (B_N(flag[i][jmax+1]) == 1) {throw std::runtime_error(std::string("impossible to have B_N at top boundary in temp"));}

		if (B_S(flag[i][jmax+1]) == 1) {T[i][jmax+1] = T[i][jmax];}

		if (B_NE(flag[i][jmax+1]) == 1) {throw std::runtime_error(std::string("impossible to have B_NE at top boundary in temp"));}

		if (B_NW(flag[i][jmax+1]) == 1) {throw std::runtime_error(std::string("impossible to have B_NW at top boundary in temp"));}

		if (B_SE(flag[i][jmax+1]) == 1) {T[i][jmax+1] = (T[i][jmax] + T[i + 1][jmax+1]) / 2;}

		if (B_SW(flag[i][jmax+1]) == 1) {T[i][jmax+1] = (T[i][jmax] + T[i - 1][jmax+1]) / 2;}


		if ((flag.at(i).at(jmax+1) & flag_10bit::DIRICHLET) == flag_10bit::DIRICHLET){
			T[i][jmax+1] = 2 * T_c - T[i][jmax];
		}

		break;

		default:
			break;

		}
	}

    // Bottom boundary
    for (int i = 1; i <= imax; ++i) {

		switch(flag.at(i).at(0) & (flag_10bit::NO_SLIP | flag_10bit::FREE_SLIP)){
		case flag_10bit::NO_SLIP: // Obstacle conditions are adiabatic
		case flag_10bit::FREE_SLIP:

		if (B_E(flag[i][0]) == 1) {T[i][0] = T[i + 1][0];}

		if (B_W(flag[i][0]) == 1) {T[i][0] = T[i - 1][0];}

		if (B_N(flag[i][0]) == 1) {T[i][0] = T[i][1];}

		if (B_S(flag[i][0]) == 1) {throw std::runtime_error(std::string("impossible to have B_S at bottom boundary in temp"));}

		if (B_NE(flag[i][0]) == 1) {T[i][0] = (T[i][1] + T[i + 1][0]) / 2;}

		if (B_NW(flag[i][0]) == 1) {T[i][0] = (T[i][1] + T[i - 1][0]) / 2;}

		if (B_SE(flag[i][0]) == 1) {throw std::runtime_error(std::string("impossible to have B_SE at bottom boundary in temp"));}

		if (B_SW(flag[i][0]) == 1) {throw std::runtime_error(std::string("impossible to have B_SW at bottom boundary in temp"));}

		if ((flag.at(i).at(0) & flag_10bit::DIRICHLET) == flag_10bit::DIRICHLET){
			T[i][0] = 2 * T_h - T[i][1];
		}

		break;

		default:
			break;

			}
		}

    // Left boundary
	for (int j = 1; j <= jmax; ++j) {

		switch(flag.at(0).at(j) & (flag_10bit::NO_SLIP | flag_10bit::FREE_SLIP)){
		case flag_10bit::NO_SLIP: // Obstacle conditions are adiabatic
		case flag_10bit::FREE_SLIP:

			if (B_E(flag[0][j]) == 1) {T[0][j] = T[1][j];}

			if (B_W(flag[0][j]) == 1) {throw std::runtime_error(std::string("impossible to have B_W at left boundary in temp"));}

			if (B_N(flag[0][j]) == 1) {T[0][j] = T[0][j + 1];}

			if (B_S(flag[0][j]) == 1) {T[0][j] = T[0][j - 1];}

			if (B_NE(flag[0][j]) == 1) {T[0][j] = (T[0][j + 1] + T[1][j]) / 2;}

			if (B_NW(flag[0][j]) == 1) {throw std::runtime_error(std::string("impossible to have B_NW at left boundary in temp"));}

			if (B_SE(flag[0][j]) == 1) {T[0][j] = (T[0][j - 1] + T[1][j]) / 2;}

			if (B_SW(flag[0][j]) == 1) {throw std::runtime_error(std::string("impossible to have B_SW at left boundary in temp"));}


			if ((flag.at(0).at(j) & flag_10bit::DIRICHLET) == flag_10bit::DIRICHLET){
					T[0][j] = 2 * T_h - T[1][j];
			}

		break;

		default:
			break;

			}
		}

	// Right Boundary

	for (int j = 1; j <= jmax; ++j) {

			switch(flag.at(imax+1).at(j) & (flag_10bit::NO_SLIP | flag_10bit::FREE_SLIP)){
			case flag_10bit::NO_SLIP: // Obstacle conditions are adiabatic
			case flag_10bit::FREE_SLIP:


			if (B_E(flag[imax+1][j]) == 1) {throw std::runtime_error(std::string("impossible to have B_E at right boundary in temp"));}

			if (B_W(flag[imax+1][j]) == 1) {T[imax+1][j] = T[imax][j];}

			if (B_N(flag[imax+1][j]) == 1) {T[imax+1][j] = T[imax+1][j + 1];}

			if (B_S(flag[imax+1][j]) == 1) {T[imax+1][j] = T[imax+1][j - 1];}

			if (B_NE(flag[imax+1][j]) == 1) {throw std::runtime_error(std::string("impossible to have B_E at right boundary in temp"));}

			if (B_NW(flag[imax+1][j]) == 1) {T[imax+1][j] = (T[imax+1][j + 1] + T[imax][j]) / 2;}

			if (B_SE(flag[imax+1][j]) == 1) {throw std::runtime_error(std::string("impossible to have B_E at right boundary in temp"));}

			if (B_SW(flag[imax+1][j]) == 1) {T[imax+1][j] = (T[imax+1][j - 1] + T[imax][j]) / 2;}


			if ((flag.at(imax+1).at(j) & flag_10bit::DIRICHLET) == flag_10bit::DIRICHLET){
				T[imax+1][j] = 2 * T_c - T[imax][j];
			}

			break;

			default:
				break;

		}
	}

	// Inner cell
	for (int i = 1; i <= imax; ++i) {

	        for (int j = 1; j <= jmax; ++j) {

	        	switch(flag.at(i).at(j) & (flag_10bit::NO_SLIP | flag_10bit::FREE_SLIP)){
	        	case flag_10bit::NO_SLIP: //No slip conditions
	        	case flag_10bit::FREE_SLIP:

	            if (B_E(flag[i][j]) == 1) {T[i][j] = T[i + 1][j];}

	            if (B_W(flag[i][j]) == 1) {T[i][j] = T[i - 1][j];}

	            if (B_N(flag[i][j]) == 1) {T[i][j] = T[i][j + 1];}

	            if (B_S(flag[i][j]) == 1) {T[i][j] = T[i][j - 1];}

	            if (B_NE(flag[i][j]) == 1) {T[i][j] = (T[i][j + 1] + T[i + 1][j]) / 2;}

	            if (B_NW(flag[i][j]) == 1) {T[i][j] = (T[i][j + 1] + T[i - 1][j]) / 2;}

	            if (B_SE(flag[i][j]) == 1) {T[i][j] = (T[i][j - 1] + T[i + 1][j]) / 2;}

	            if (B_SW(flag[i][j]) == 1) {T[i][j] = (T[i][j - 1] + T[i - 1][j]) / 2;}

	            break;

	        	default:
	        		break;
	        	}
	        }
	    }

    grid.set_temperature(T);

}



// ------------------------------------------- check conditions for obstacle cells -----------------------------------------------------

int B_E(unsigned int flag){
    return (((flag & flag_10bit::E) == flag_10bit::E) && !(((flag & flag_10bit::N) == flag_10bit::N) || ((flag & flag_10bit::S) == flag_10bit::S)));
}

int B_W(unsigned int flag){
    return (((flag & flag_10bit::W) == flag_10bit::W) && !(((flag & flag_10bit::N) == flag_10bit::N) || ((flag & flag_10bit::S) == flag_10bit::S)));
}

int B_N(unsigned int flag){
    return (((flag & flag_10bit::N) == flag_10bit::N) && !(((flag & flag_10bit::E) == flag_10bit::E) || ((flag & flag_10bit::W) == flag_10bit::W)));
}

int B_S(unsigned int flag){
    return (((flag & flag_10bit::S) == flag_10bit::S) && !(((flag & flag_10bit::E) == flag_10bit::E) || ((flag & flag_10bit::S) == flag_10bit::W)));
}

int B_NE(unsigned int flag){
    return (((flag & flag_10bit::N) == flag_10bit::N) &&  ((flag & flag_10bit::E) == flag_10bit::E));
}

int B_NW(unsigned int flag){
    return (((flag & flag_10bit::N) == flag_10bit::N) &&  ((flag & flag_10bit::W) == flag_10bit::W));
}

int B_SE(unsigned int flag){
    return (((flag & flag_10bit::S) == flag_10bit::S) &&  ((flag & flag_10bit::E) == flag_10bit::E));
}

int B_SW(unsigned int flag){
    return (((flag & flag_10bit::S) == flag_10bit::S) &&  ((flag & flag_10bit::W) == flag_10bit::W));
}


void boundaryvalues_uv(int imax,int jmax, Grid& grid,const matrix<unsigned int> &flag){

    matrix<double> U,V;
    grid.velocity(U,velocity_type::U);
    grid.velocity(V,velocity_type::V);



    // ---------------------------------------------------- left boundary --------------------------------------------------------------
    for(int j = 1; j<=jmax; ++j){

        switch(flag.at(0).at(j) & (flag_10bit::NO_SLIP | flag_10bit::FREE_SLIP | flag_10bit::OUTFLOW | flag_10bit::INFLOW )){
            case flag_10bit::NO_SLIP: //No slip conditions
                if (B_N(flag.at(0).at(j)) == 1){
                    V.at(0).at(j) = 0;
                    U.at(0).at(j) = -U.at(0).at(j+1);
                }

                if (B_W(flag.at(0).at(j)) == 1){
                    throw std::runtime_error(std::string("impossible to have B_W at left boundary "));
                }

                if (B_S(flag.at(0).at(j)) == 1){
                    V.at(0).at(j-1) = 0;
                    U.at(0).at(j) = -U.at(0).at(j-1);
                }

                if (B_E(flag.at(0).at(j)) == 1){
                    U.at(0).at(j) = 0;
                    V.at(0).at(j-1) = -V.at(1).at(j-1);
                    V.at(0).at(j) = -V.at(1).at(j);
                }

                if (B_NE(flag.at(0).at(j)) == 1){
                    U.at(0).at(j) = 0;
                    V.at(0).at(j) = 0;
                    V.at(0).at(j-1) = -V.at(1).at(j-1);
                }

                if (B_NW(flag.at(0).at(j)) == 1){
                    throw std::runtime_error(std::string("impossible to have B_NW at left boundary "));
                }

                if (B_SE(flag.at(0).at(j)) == 1){
                    U.at(0).at(j)=0;
                    V.at(0).at(j-1)=0;
                    V.at(0).at(j) = -V.at(1).at(j);
                }

                if (B_SW(flag.at(0).at(j)) == 1){
                    throw std::runtime_error(std::string("impossible to have B_NW at left boundary "));
                }
                break;
            case flag_10bit::FREE_SLIP://Free Slip conditions
                if (B_E(flag.at(0).at(j)) == 1){
                    U.at(0).at(j) = 0;
                    V.at(0).at(j) = V.at(1).at(j);
                    V.at(0).at(j-1) = V.at(1).at(j-1);
                }

                if (B_W(flag.at(0).at(j)) == 1) {
                    throw std::runtime_error(std::string("impossible to have B_W at left boundary "));
                }

                if (B_N(flag.at(0).at(j)) == 1) {
                    V.at(0).at(j) = 0;
                    U.at(0).at(j) = U.at(0).at(j+1);
                }

                if (B_S(flag.at(0).at(j)) == 1){
                    V.at(0).at(j-1) = 0;
                    U.at(0).at(j) = U.at(0).at(j-1);
                }

                if (B_NE(flag.at(0).at(j)) == 1){
                    U.at(0).at(j) = 0;
                    V.at(0).at(j) = 0;
                    V.at(0).at(j-1) = V.at(1).at(j-1);
                }

                if (B_NW(flag.at(0).at(j)) == 1){
                    throw std::runtime_error(std::string("impossible to have B_NW at left boundary "));
                }

                if (B_SE(flag.at(0).at(j)) == 1){
                    U.at(0).at(j)=0;
                    V.at(0).at(j-1)=0;
                    V.at(0).at(j) = V.at(1).at(j);
                }

                if (B_SW(flag.at(0).at(j)) == 1){
                    throw std::runtime_error(std::string("impossible to have B_SW at left boundary "));
                }
                break;

            case flag_10bit::OUTFLOW: // assuming only x-direction
                throw std::runtime_error(std::string("impossible to have outflow at left boundary (only x direction assumed"));
                break;

            case flag_10bit::INFLOW: // assuming only x-direction

                U.at(0).at(j)=1;
                V.at(0).at(j) = 0;
                V.at(0).at(j-1) = 0;
                break;

            default:
                break;

        }

    }

    // ------------------------------------------------------- Right Boundary ----------------------------------------------------------------
    for(int j = 1; j<=jmax; ++j){

        switch(flag.at(imax+1).at(j) & (flag_10bit::NO_SLIP | flag_10bit::FREE_SLIP | flag_10bit::OUTFLOW | flag_10bit::INFLOW)){
            case flag_10bit::NO_SLIP: //No slip conditions
                if (B_N(flag.at(imax+1).at(j)) == 1){
                    V.at(imax+1).at(j) = 0;
                    U.at(imax).at(j) = -U.at(imax).at(j+1);
                    U.at(imax+1).at(j) = -U.at(imax+1).at(j+1);
                }

                if (B_W(flag.at(imax+1).at(j)) == 1){
                    U.at(imax).at(j) = 0;
                    V.at(imax+1).at(j-1) = -V.at(imax).at(j-1);
                    V.at(imax+1).at(j) = -V.at(imax).at(j);
                }

                if (B_S(flag.at(imax+1).at(j)) == 1){
                    V.at(imax+1).at(j-1) = 0;
                    U.at(imax).at(j) = -U.at(imax).at(j-1);
                    U.at(imax+1).at(j) = -U.at(imax+1).at(j-1);
                }

                if (B_E(flag.at(imax+1).at(j)) == 1){
                    throw std::runtime_error(std::string("impossible to have B_E at right boundary "));
                }

                if (B_NE(flag.at(imax+1).at(j)) == 1){
                    throw std::runtime_error(std::string("impossible to have B_NE at right boundary "));
                }

                if (B_NW(flag.at(imax+1).at(j)) == 1){
                    U.at(imax).at(j) = 0;
                    U.at(imax+1).at(j) = - U.at(imax+1).at(j+1);
                    V.at(imax+1).at(j) = 0;
                    V.at(imax+1).at(j-1) = -V.at(imax).at(j-1);
                }

                if (B_SE(flag.at(imax+1).at(j)) == 1){
                    throw std::runtime_error(std::string("impossible to have B_SE at right boundary "));
                }

                if (B_SW(flag.at(imax+1).at(j)) == 1){
                    U.at(imax).at(j) = 0;
                    U.at(imax+1).at(j) = -U.at(imax+1).at(j-1);
                    V.at(imax+1).at(j-1) = 0;
                    V.at(imax+1).at(j) = -V.at(imax).at(j);
                }
                break;

            case flag_10bit::FREE_SLIP://Free Slip conditions
                if (B_E(flag.at(imax+1).at(j)) == 1){
                    throw std::runtime_error(std::string("impossible to have B_E at right boundary "));
                }

                if (B_W(flag.at(imax+1).at(j)) == 1) {
                    U.at(imax).at(j) = 0;
                    V.at(imax+1).at(j-1) = V.at(imax).at(j-1);
                    V.at(imax+1).at(j) = V.at(imax).at(j);
                }

                if (B_N(flag.at(imax+1).at(j)) == 1) {
                    V.at(imax+1).at(j) = 0;
                    U.at(imax).at(j) = U.at(imax).at(j+1);
                    U.at(imax+1).at(j) = U.at(imax+1).at(j+1);
                }

                if (B_S(flag.at(imax+1).at(j)) == 1){
                    V.at(imax+1).at(j-1) = 0;
                    U.at(imax).at(j) = U.at(imax).at(j-1);
                    U.at(imax+1).at(j) = U.at(imax+1).at(j-1);
                }

                if (B_NE(flag.at(imax+1).at(j)) == 1){
                    throw std::runtime_error(std::string("impossible to have B_NE at right boundary "));
                }

                if (B_NW(flag.at(imax+1).at(j)) == 1){
                    U.at(imax).at(j) = 0;
                    U.at(imax+1).at(j) = U.at(imax+1).at(j+1);
                    V.at(imax+1).at(j) = 0;
                    V.at(imax+1).at(j-1) = V.at(imax).at(j-1);
                }

                if (B_SE(flag.at(imax+1).at(j)) == 1){
                    throw std::runtime_error(std::string("impossible to have B_SE at right boundary "));
                }

                if (B_SW(flag.at(imax+1).at(j)) == 1){
                    U.at(imax).at(j) = 0;
                    U.at(imax+1).at(j) = U.at(imax+1).at(j-1);
                    V.at(imax+1).at(j-1) = 0;
                    V.at(imax+1).at(j) = V.at(imax).at(j);
                }
                break;

            case flag_10bit::OUTFLOW: // assuming only x-direction
                U.at(imax+1).at(j) = U.at(imax).at(j);
                V.at(imax+1).at(j) = V.at(imax).at(j);
                V.at(imax+1).at(j-1) = V.at(imax).at(j-1);
                break;

            case flag_10bit::INFLOW: // assuming only x-direction
                U.at(imax+1).at(j)=1;
                V.at(imax+1).at(j) = 0;
                V.at(imax+1).at(j-1) = 0;
                break;

            default:
                break;
        }
    }

    //  ---------------------------------------------------- Top Boundary ----------------------------------------------------------------
    for(int i = 1; i<=imax; ++i){
        switch(flag.at(i).at(jmax+1) & (flag_10bit::NO_SLIP | flag_10bit::FREE_SLIP | flag_10bit::OUTFLOW | flag_10bit::INFLOW)){
            case flag_10bit::NO_SLIP: //No slip conditions
                if (B_N(flag.at(i).at(jmax+1)) == 1){
                    throw std::runtime_error(std::string("impossible to have B_N at top boundary "));
                }

                if (B_W(flag.at(i).at(jmax+1)) == 1){
                    U.at(i-1).at(jmax+1) = 0;
                    V.at(i).at(jmax) = -V.at(i-1).at(jmax);
                    V.at(i).at(jmax+1) = -V.at(i-1).at(jmax+1);
                }

                if (B_S(flag.at(i).at(jmax+1)) == 1){
                    V.at(i).at(jmax) = 0;
                    U.at(i-1).at(jmax+1) = -U.at(i-1).at(jmax);
                    U.at(i).at(jmax+1) = -U.at(i).at(jmax);
                }

                if (B_E(flag.at(i).at(jmax+1)) == 1){
                    U.at(i).at(jmax+1) = 0;
                    V.at(i).at(jmax) = -V.at(i+1).at(jmax);
                    V.at(i).at(jmax+1) = -V.at(i+1).at(jmax+1);
                }

                if (B_NE(flag.at(i).at(jmax+1)) == 1){
                    throw std::runtime_error(std::string("impossible to have B_NE at top boundary "));
                }

                if (B_NW(flag.at(i).at(jmax+1)) == 1){
                    throw std::runtime_error(std::string("impossible to have B_NW at top boundary "));
                }

                if (B_SE(flag.at(i).at(jmax+1)) == 1){
                    U.at(i).at(jmax+1)=0;
                    U.at(i-1).at(jmax+1) = -U.at(i-1).at(jmax);
                    V.at(i).at(jmax)=0;
                    V.at(i).at(jmax+1) = -V.at(i+1).at(jmax+1);
                }

                if (B_SW(flag.at(i).at(jmax+1)) == 1){
                    U.at(i-1).at(jmax+1) = 0;
                    U.at(i).at(jmax+1) = -U.at(i).at(jmax);
                    V.at(i).at(jmax) = 0;
                    V.at(i).at(jmax+1) = -V.at(i-1).at(jmax+1);
                }
                break;

            case flag_10bit::FREE_SLIP://Free Slip conditions
                if (B_E(flag.at(i).at(jmax+1)) == 1){
                    U.at(i).at(jmax+1) = 0;
                    V.at(i).at(jmax+1) = V.at(i+1).at(jmax+1);
                    V.at(i).at(jmax) = V.at(i+1).at(jmax);
                }

                if (B_W(flag.at(i).at(jmax+1)) == 1) {
                    U.at(i-1).at(jmax+1) = 0;
                    V.at(i).at(jmax) = V.at(i-1).at(jmax);
                    V.at(i).at(jmax+1) = V.at(i-1).at(jmax+1);
                }

                if (B_N(flag.at(i).at(jmax+1)) == 1) {
                    throw std::runtime_error(std::string("impossible to have B_N at top boundary "));
                }

                if (B_S(flag.at(i).at(jmax+1)) == 1){
                    V.at(i).at(jmax) = 0;
                    U.at(i-1).at(jmax+1) = U.at(i-1).at(jmax);
                    U.at(i).at(jmax+1) = U.at(i).at(jmax);
                }

                if (B_NE(flag.at(i).at(jmax+1)) == 1){
                    throw std::runtime_error(std::string("impossible to have B_NE at top boundary "));
                }

                if (B_NW(flag.at(i).at(jmax+1)) == 1){
                    throw std::runtime_error(std::string("impossible to have B_NW at top boundary "));
                }

                if (B_SE(flag.at(i).at(jmax+1)) == 1){
                    U.at(i).at(jmax+1)=0;
                    U.at(i-1).at(jmax+1) = U.at(i-1).at(jmax);
                    V.at(i).at(jmax)=0;
                    V.at(i).at(jmax+1) = V.at(i+1).at(jmax+1);
                }

                if (B_SW(flag.at(i).at(jmax+1)) == 1){
                    U.at(i-1).at(jmax+1) = 0;
                    U.at(i).at(jmax+1) = U.at(i).at(jmax);
                    V.at(i).at(jmax) = 0;
                    V.at(i).at(jmax+1) = V.at(i-1).at(jmax+1);
                }
                break;

            case flag_10bit::OUTFLOW: // assuming only x-direction
                U.at(i).at(jmax+1) = U.at(i-1).at(jmax+1);
                V.at(i).at(jmax+1) = V.at(i-1).at(jmax+1);
                V.at(i).at(jmax) = V.at(i-1).at(jmax);
                break;

            case flag_10bit::INFLOW: // assuming only x-direction
                U.at(i).at(jmax+1)=1;
                V.at(i).at(jmax+1) = 0;
                V.at(i).at(jmax) = 0;
                break;

            default:
                break;
        }

    }

    // ------------------------------------------------- Bottom Boundary -----------------------------------------------------------------
    for(int i = 1; i<=imax; ++i){
        switch(flag.at(i).at(0) & (flag_10bit::NO_SLIP | flag_10bit::FREE_SLIP | flag_10bit::OUTFLOW | flag_10bit::INFLOW)){
            case flag_10bit::NO_SLIP: //No slip conditions
                if (B_N(flag.at(i).at(0)) == 1){
                    V.at(i).at(0) = 0;
                    U.at(i-1).at(0) = -U.at(i-1).at(1);
                    U.at(i).at(0) = -U.at(i).at(1);
                }

                if (B_W(flag.at(i).at(0)) == 1){
                    U.at(i-1).at(0) = 0;
                    V.at(i).at(0) = -V.at(i-1).at(0);
                }

                if (B_S(flag.at(i).at(0)) == 1){
                    throw std::runtime_error(std::string("impossible to have B_S at bottom boundary "));
                }

                if (B_E(flag.at(i).at(0)) == 1){
                    U.at(i).at(0) = 0;
                    V.at(i).at(0) = -V.at(i+1).at(0);
                }

                if (B_NE(flag.at(i).at(0)) == 1){
                    U.at(i).at(0) = 0;
                    U.at(i-1).at(0) = -U.at(i-1).at(1);
                    V.at(i).at(0) = 0;
                }

                if (B_NW(flag.at(i).at(0)) == 1){
                    U.at(i-1).at(0) = 0;
                    U.at(i).at(0) = - U.at(i).at(1);
                    V.at(i).at(0) = 0;
                }

                if (B_SE(flag.at(i).at(0)) == 1){
                    throw std::runtime_error(std::string("impossible to have B_SE at bottom boundary "));
                }

                if (B_SW(flag.at(i).at(0)) == 1){
                    throw std::runtime_error(std::string("impossible to have B_SW at bottom boundary "));
                }
                break;

            case flag_10bit::FREE_SLIP://Free Slip conditions
                if (B_E(flag.at(i).at(0)) == 1){
                    U.at(i).at(0) = 0;
                    V.at(i).at(0) = V.at(i+1).at(0);
                }

                if (B_W(flag.at(i).at(0)) == 1) {
                    U.at(i-1).at(0) = 0;
                    V.at(i).at(0) = V.at(i-1).at(0);
                }

                if (B_N(flag.at(i).at(0)) == 1) {
                    V.at(i).at(0) = 0;
                    U.at(i-1).at(0) = U.at(i-1).at(1);
                    U.at(i).at(0) = U.at(i).at(1);
                }

                if (B_S(flag.at(i).at(0)) == 1){
                    throw std::runtime_error(std::string("impossible to have B_S at top boundary "));
                }

                if (B_NE(flag.at(i).at(0)) == 1){
                    U.at(i).at(0) = 0;
                    U.at(i-1).at(0) = U.at(i-1).at(1);
                    V.at(i).at(0) = 0;
                }

                if (B_NW(flag.at(i).at(0)) == 1){
                    U.at(i-1).at(0) = 0;
                    U.at(i).at(0) = U.at(i).at(1);
                    V.at(i).at(0) = 0;
                }

                if (B_SE(flag.at(i).at(0)) == 1){
                    throw std::runtime_error(std::string("impossible to have B_SE at top boundary "));
                }

                if (B_SW(flag.at(i).at(0)) == 1){
                    throw std::runtime_error(std::string("impossible to have B_SW at top boundary "));
                }
                break;

            case flag_10bit::OUTFLOW: // assuming only x-direction
                U.at(i).at(0) = U.at(i-1).at(0);
                V.at(i).at(0) = V.at(i-1).at(0);
                break;

            case flag_10bit::INFLOW: // assuming only x-direction
                U.at(i).at(0)=1;
                V.at(i).at(0) = 0;
                break;

            default:
                break;
        }
    }

    // ---------------------------------------------- Top-Left Corner ----------------------------------------------------------------------
    switch(flag.at(0).at(jmax+1) & (flag_10bit::NO_SLIP | flag_10bit::FREE_SLIP | flag_10bit::OUTFLOW | flag_10bit::INFLOW)){
        case flag_10bit::NO_SLIP: //No slip conditions
            if (B_N(flag.at(0).at(jmax+1)) == 1){
                throw std::runtime_error(std::string("impossible to have B_N at top-left boundary "));
            }

            if (B_W(flag.at(0).at(jmax+1)) == 1){
                throw std::runtime_error(std::string("impossible to have B_W at top-left boundary "));
            }

            if (B_S(flag.at(0).at(jmax+1)) == 1){
                V.at(0).at(jmax) = 0;
                U.at(0).at(jmax+1) = -U.at(0).at(jmax);
            }

            if (B_E(flag.at(0).at(jmax+1)) == 1){
                U.at(0).at(jmax+1) = 0;
                V.at(0).at(jmax) = -V.at(1).at(jmax);
                V.at(0).at(jmax+1) = -V.at(1).at(jmax+1);
            }

            if (B_NE(flag.at(0).at(jmax+1)) == 1){
                throw std::runtime_error(std::string("impossible to have B_NE at top-left boundary "));
            }

            if (B_NW(flag.at(0).at(jmax+1)) == 1){
                throw std::runtime_error(std::string("impossible to have B_NW at top-left boundary "));
            }

            if (B_SE(flag.at(0).at(jmax+1)) == 1){
                U.at(0).at(jmax+1)=0;
                V.at(0).at(jmax)=0;
                V.at(0).at(jmax+1) = -V.at(1).at(jmax+1);
            }

            if (B_SW(flag.at(0).at(jmax+1)) == 1){
                throw std::runtime_error(std::string("impossible to have B_SW at top-left boundary "));
            }
            break;

        case flag_10bit::FREE_SLIP://Free Slip conditions
            if (B_E(flag.at(0).at(jmax+1)) == 1){
                U.at(0).at(jmax+1) = 0;
                V.at(0).at(jmax+1) = V.at(1).at(jmax+1);
                V.at(0).at(jmax) = V.at(1).at(jmax);
            }

            if (B_W(flag.at(0).at(jmax+1)) == 1) {
                throw std::runtime_error(std::string("impossible to have B_W at top-left boundary "));
            }

            if (B_N(flag.at(0).at(jmax+1)) == 1) {
                throw std::runtime_error(std::string("impossible to have B_N at top-left boundary "));
            }

            if (B_S(flag.at(0).at(jmax+1)) == 1){
                V.at(0).at(jmax) = 0;
                U.at(0).at(jmax+1) = U.at(0).at(jmax);
            }

            if (B_NE(flag.at(0).at(jmax+1)) == 1){
                throw std::runtime_error(std::string("impossible to have B_NE at top-left boundary "));
            }

            if (B_NW(flag.at(0).at(jmax+1)) == 1){
                throw std::runtime_error(std::string("impossible to have B_NW at top-left boundary "));
            }

            if (B_SE(flag.at(0).at(jmax+1)) == 1){
                U.at(0).at(jmax+1)=0;
                V.at(0).at(jmax)=0;
                V.at(0).at(jmax+1) = V.at(1).at(jmax+1);
            }

            if (B_SW(flag.at(0).at(jmax+1)) == 1){
                throw std::runtime_error(std::string("impossible to have B_SW at top-left boundary "));
            }
            break;

        case flag_10bit::OUTFLOW: // assuming only x-direction
            throw std::runtime_error(std::string("impossible to have outflow at top-left boundary (only x direction assumed"));
            break;

        case flag_10bit::INFLOW: // assuming only x-direction
            U.at(0).at(jmax+1)=1;
            V.at(0).at(jmax+1) = 0;
            V.at(0).at(jmax) = 0;
            break;

        default:
            break;
    }

    // ------------------------------------------------- Bottom-left Corner ------------------------------------------------------------------
    switch(flag.at(0).at(0) & (flag_10bit::NO_SLIP | flag_10bit::FREE_SLIP | flag_10bit::OUTFLOW | flag_10bit::INFLOW)){
        case flag_10bit::NO_SLIP: //No slip conditions
            if (B_N(flag.at(0).at(0)) == 1){
                V.at(0).at(0) = 0;
                U.at(0).at(0) = -U.at(0).at(1);
            }

            if (B_W(flag.at(0).at(0)) == 1){
                throw std::runtime_error(std::string("impossible to have B_W at bottom-left boundary"));
            }

            if (B_S(flag.at(0).at(0)) == 1){
                throw std::runtime_error(std::string("impossible to have B_S at bottom-left boundary"));
            }

            if (B_E(flag.at(0).at(0)) == 1){
                U.at(0).at(0) = 0;
                V.at(0).at(0) = -V.at(1).at(0);
            }

            if (B_NE(flag.at(0).at(0)) == 1){
                U.at(0).at(0) = 0;
                V.at(0).at(0) = 0;
            }

            if (B_NW(flag.at(0).at(0)) == 1){
                throw std::runtime_error(std::string("impossible to have B_NW at bottom-left boundary"));
            }

            if (B_SE(flag.at(0).at(0)) == 1){
                throw std::runtime_error(std::string("impossible to have B_SE at bottom-left boundary"));
            }

            if (B_SW(flag.at(0).at(0)) == 1){
                throw std::runtime_error(std::string("impossible to have B_SW at bottom-left boundary"));
            }
            break;

        case flag_10bit::FREE_SLIP://Free Slip conditions
            if (B_E(flag.at(0).at(0)) == 1){
                U.at(0).at(0) = 0;
                V.at(0).at(0) = V.at(1).at(0);
            }

            if (B_W(flag.at(0).at(0)) == 1) {
                throw std::runtime_error(std::string("impossible to have B_W at bottom-left boundary"));
            }

            if (B_N(flag.at(0).at(0)) == 1) {
                V.at(0).at(0) = 0;
                U.at(0).at(0) = U.at(0).at(1);
            }

            if (B_S(flag.at(0).at(0)) == 1){
                throw std::runtime_error(std::string("impossible to have B_S at bottom-left boundary"));
            }

            if (B_NE(flag.at(0).at(0)) == 1){
                U.at(0).at(0) = 0;
                V.at(0).at(0) = 0;
            }

            if (B_NW(flag.at(0).at(0)) == 1){
                throw std::runtime_error(std::string("impossible to have B_NW at bottom-left boundary"));
            }

            if (B_SE(flag.at(0).at(0)) == 1){
                throw std::runtime_error(std::string("impossible to have B_SE at bottom-left boundary"));
            }

            if (B_SW(flag.at(0).at(0)) == 1){
                throw std::runtime_error(std::string("impossible to have B_SW at bottom-left boundary"));
            }
            break;

        case flag_10bit::OUTFLOW: // assuming only x-direction
            throw std::runtime_error(std::string("impossible to have outflow at bottom-left boundary (only x direction)"));
            break;

        case flag_10bit::INFLOW: // assuming only x-direction
            U.at(0).at(0)= 1;
            V.at(0).at(0) = 0;
            break;

        default:
            break;
    }

    // ----------------------------------------------------- Top-right corner -----------------------------------------------------------
    switch(flag.at(imax+1).at(jmax+1) & (flag_10bit::NO_SLIP | flag_10bit::FREE_SLIP | flag_10bit::OUTFLOW | flag_10bit::INFLOW)){
        case flag_10bit::NO_SLIP: //No slip conditions
            if (B_N(flag.at(imax+1).at(jmax+1)) == 1){
                throw std::runtime_error(std::string("impossible to have B_N at top-right boundary"));
            }

            if (B_W(flag.at(imax+1).at(jmax+1)) == 1){
                U.at(imax).at(jmax+1) = 0;
                V.at(imax+1).at(jmax) = -V.at(imax).at(jmax);
                V.at(imax+1).at(jmax+1) = -V.at(imax).at(jmax+1);
            }

            if (B_S(flag.at(imax+1).at(jmax+1)) == 1){
                V.at(imax+1).at(jmax) = 0;
                U.at(imax).at(jmax+1) = -U.at(imax).at(jmax);
                U.at(imax+1).at(jmax+1) = -U.at(imax+1).at(jmax);
            }

            if (B_E(flag.at(imax+1).at(jmax+1)) == 1){
                throw std::runtime_error(std::string("impossible to have B_E at top-right boundary"));
            }

            if (B_NE(flag.at(imax+1).at(jmax+1)) == 1){
                throw std::runtime_error(std::string("impossible to have B_NE at top-right boundary"));
            }

            if (B_NW(flag.at(imax+1).at(jmax+1)) == 1){
                throw std::runtime_error(std::string("impossible to have B_NW at top-right boundary"));
            }

            if (B_SE(flag.at(imax+1).at(jmax+1)) == 1){
                throw std::runtime_error(std::string("impossible to have B_SE at top-right boundary"));
            }

            if (B_SW(flag.at(imax+1).at(jmax+1)) == 1){
                U.at(imax).at(jmax+1) = 0;
                U.at(imax+1).at(jmax+1) = -U.at(imax+1).at(jmax);
                V.at(imax+1).at(jmax) = 0;
                V.at(imax+1).at(jmax+1) = -V.at(imax).at(jmax+1);
            }
            break;

        case flag_10bit::FREE_SLIP://Free Slip conditions
            if (B_E(flag.at(imax+1).at(jmax+1)) == 1){
                throw std::runtime_error(std::string("impossible to have B_E at top-right boundary"));
            }

            if (B_W(flag.at(imax+1).at(jmax+1)) == 1) {
                U.at(imax).at(jmax+1) = 0;
                V.at(imax+1).at(jmax) = V.at(imax).at(jmax);
                V.at(imax+1).at(jmax+1) = V.at(imax).at(jmax+1);
            }

            if (B_N(flag.at(imax+1).at(jmax+1)) == 1) {
                throw std::runtime_error(std::string("impossible to have B_N at top-right boundary"));
            }

            if (B_S(flag.at(imax+1).at(jmax+1)) == 1){
                V.at(imax+1).at(jmax) = 0;
                U.at(imax).at(jmax+1) = U.at(imax).at(jmax);
                U.at(imax+1).at(jmax+1) = U.at(imax+1).at(jmax);
            }

            if (B_NE(flag.at(imax+1).at(jmax+1)) == 1){
                throw std::runtime_error(std::string("impossible to have B_NE at top-right boundary"));
            }

            if (B_NW(flag.at(imax+1).at(jmax+1)) == 1){
                throw std::runtime_error(std::string("impossible to have B_NW at top-right boundary"));
            }

            if (B_SE(flag.at(imax+1).at(jmax+1)) == 1){
                throw std::runtime_error(std::string("impossible to have B_SE at top-right boundary"));
            }

            if (B_SW(flag.at(imax+1).at(jmax+1)) == 1){
                U.at(imax).at(jmax+1) = 0;
                U.at(imax+1).at(jmax+1) = U.at(imax+1).at(jmax);
                V.at(imax+1).at(jmax) = 0;
                V.at(imax+1).at(jmax+1) = V.at(imax).at(jmax+1);
            }
            break;

        case flag_10bit::OUTFLOW: // assuming only x-direction
            U.at(imax+1).at(jmax+1) = U.at(imax).at(jmax+1);
            V.at(imax+1).at(jmax+1) = V.at(imax).at(jmax+1);
            V.at(imax+1).at(jmax) = V.at(imax).at(jmax);
            break;

        case flag_10bit::INFLOW: // assuming only x-direction
            U.at(imax+1).at(jmax+1)=1;
            V.at(imax+1).at(jmax+1) = 0;
            V.at(imax+1).at(jmax) = 0;
            break;

        default:
            break;
    }

    // ------------------------------------------------------ Bottom-right corner -----------------------------------------------------------
    switch(flag.at(imax+1).at(0) & (flag_10bit::NO_SLIP | flag_10bit::FREE_SLIP | flag_10bit::OUTFLOW | flag_10bit::INFLOW)){
        case flag_10bit::NO_SLIP: //No slip conditions
            if (B_N(flag.at(imax+1).at(0)) == 1){
                V.at(imax+1).at(0) = 0;
                U.at(imax).at(0) = -U.at(imax).at(1);
                U.at(imax+1).at(0) = -U.at(imax+1).at(1);
            }

            if (B_W(flag.at(imax+1).at(0)) == 1){
                U.at(imax).at(0) = 0;
                V.at(imax+1).at(0) = -V.at(imax).at(0);
            }

            if (B_S(flag.at(imax+1).at(0)) == 1){
                throw std::runtime_error(std::string("impossible to have B_S at bottom-right boundary"));
            }

            if (B_E(flag.at(imax+1).at(0)) == 1){
                throw std::runtime_error(std::string("impossible to have B_E at bottom-right boundary"));
            }

            if (B_NE(flag.at(imax+1).at(0)) == 1){
                throw std::runtime_error(std::string("impossible to have B_NE at bottom-right boundary"));
            }

            if (B_NW(flag.at(imax+1).at(0)) == 1){
                U.at(imax).at(0) = 0;
                U.at(imax+1).at(0) = - U.at(imax+1).at(1);
                V.at(imax+1).at(0) = 0;
            }

            if (B_SE(flag.at(imax+1).at(0)) == 1){
                throw std::runtime_error(std::string("impossible to have B_SE at bottom-right boundary"));
            }

            if (B_SW(flag.at(imax+1).at(0)) == 1){
                throw std::runtime_error(std::string("impossible to have B_SW at bottom-right boundary"));
            }
            break;

        case flag_10bit::FREE_SLIP://Free Slip conditions
            if (B_E(flag.at(imax+1).at(0)) == 1){
                throw std::runtime_error(std::string("impossible to have B_E at bottom-right boundary"));
            }

            if (B_W(flag.at(imax+1).at(0)) == 1) {
                U.at(imax).at(0) = 0;
                V.at(imax+1).at(0) = V.at(imax).at(0);
            }

            if (B_N(flag.at(imax+1).at(0)) == 1) {
                V.at(imax+1).at(0) = 0;
                U.at(imax).at(0) = U.at(imax).at(1);
                U.at(imax+1).at(0) = U.at(imax+1).at(1);
            }

            if (B_S(flag.at(imax+1).at(0)) == 1){
                throw std::runtime_error(std::string("impossible to have B_S at bottom-right boundary"));
            }

            if (B_NE(flag.at(imax+1).at(0)) == 1){
                throw std::runtime_error(std::string("impossible to have B_NE at bottom-right boundary"));
            }

            if (B_NW(flag.at(imax+1).at(0)) == 1){
                U.at(imax).at(0) = 0;
                U.at(imax+1).at(0) = U.at(imax+1).at(1);
                V.at(imax+1).at(0) = 0;
            }

            if (B_SE(flag.at(imax+1).at(0)) == 1){
                throw std::runtime_error(std::string("impossible to have B_SE at bottom-right boundary"));
            }

            if (B_SW(flag.at(imax+1).at(0)) == 1){
                throw std::runtime_error(std::string("impossible to have B_SW at bottom-right boundary"));
            }
            break;

        case flag_10bit::OUTFLOW: // assuming only x-direction
            U.at(imax+1).at(0) = U.at(imax).at(0);
            V.at(imax+1).at(0) = V.at(imax).at(0);
            break;

        case flag_10bit::INFLOW: // assuming only x-direction
            U.at(imax+1).at(0)=1;
            V.at(imax+1).at(0) = 0;
            break;

        default:
            break;
    }

    // ---------------------------------------------------------- Inner Cells --------------------------------------------------------------------
    for(int i = 1; i<=imax; ++i){

        for(int j = 1; j<=jmax; ++j){

            switch(flag.at(i).at(j) & (flag_10bit::NO_SLIP | flag_10bit::FREE_SLIP | flag_10bit::OUTFLOW | flag_10bit::INFLOW )){
                case flag_10bit::NO_SLIP: //No slip conditions
                    if (B_N(flag.at(i).at(j)) == 1){
                        V.at(i).at(j) = 0;
                        U.at(i-1).at(j) = -U.at(i-1).at(j+1);
                        U.at(i).at(j) = -U.at(i).at(j+1);
                    }

                    if (B_W(flag.at(i).at(j)) == 1){
                        U.at(i-1).at(j) = 0;
                        V.at(i).at(j-1) = -V.at(i-1).at(j-1);
                        V.at(i).at(j) = -V.at(i-1).at(j);
                    }

                    if (B_S(flag.at(i).at(j)) == 1){
                        V.at(i).at(j-1) = 0;
                        U.at(i-1).at(j) = -U.at(i-1).at(j-1);
                        U.at(i).at(j) = -U.at(i).at(j-1);
                    }

                    if (B_E(flag.at(i).at(j)) == 1){
                        U.at(i).at(j) = 0;
                        V.at(i).at(j-1) = -V.at(i+1).at(j-1);
                        V.at(i).at(j) = -V.at(i+1).at(j);
                    }

                    if (B_NE(flag.at(i).at(j)) == 1){
                        U.at(i).at(j) = 0;
                        U.at(i-1).at(j) = -U.at(i-1).at(j+1);
                        V.at(i).at(j) = 0;
                        V.at(i).at(j-1) = -V.at(i+1).at(j-1);
                    }

                    if (B_NW(flag.at(i).at(j)) == 1){
                        U.at(i-1).at(j) = 0;
                        U.at(i).at(j) = - U.at(i).at(j+1);
                        V.at(i).at(j) = 0;
                        V.at(i).at(j-1) = -V.at(i-1).at(j-1);
                    }

                    if (B_SE(flag.at(i).at(j)) == 1){
                        U.at(i).at(j)=0;
                        U.at(i-1).at(j) = -U.at(i-1).at(j-1);
                        V.at(i).at(j-1)=0;
                        V.at(i).at(j) = -V.at(i+1).at(j);
                    }

                    if (B_SW(flag.at(i).at(j)) == 1){
                        U.at(i-1).at(j) = 0;
                        U.at(i).at(j) = -U.at(i).at(j-1);
                        V.at(i).at(j-1) = 0;
                        V.at(i).at(j) = -V.at(i-1).at(j);
                    }
                    break;

                case flag_10bit::FREE_SLIP://Free Slip conditions
                    if (B_E(flag.at(i).at(j)) == 1){
                        U.at(i).at(j) = 0;
                        V.at(i).at(j) = V.at(i+1).at(j);
                        V.at(i).at(j-1) = V.at(i+1).at(j-1);
                    }

                    if (B_W(flag.at(i).at(j)) == 1) {
                        U.at(i-1).at(j) = 0;
                        V.at(i).at(j-1) = V.at(i-1).at(j-1);
                        V.at(i).at(j) = V.at(i-1).at(j);
                    }

                    if (B_N(flag.at(i).at(j)) == 1) {
                        V.at(i).at(j) = 0;
                        U.at(i-1).at(j) = U.at(i-1).at(j+1);
                        U.at(i).at(j) = U.at(i).at(j+1);
                    }

                    if (B_S(flag.at(i).at(j)) == 1){
                        V.at(i).at(j-1) = 0;
                        U.at(i-1).at(j) = U.at(i-1).at(j-1);
                        U.at(i).at(j) = U.at(i).at(j-1);
                    }

                    if (B_NE(flag.at(i).at(j)) == 1){
                        U.at(i).at(j) = 0;
                        U.at(i-1).at(j) = U.at(i-1).at(j+1);
                        V.at(i).at(j) = 0;
                        V.at(i).at(j-1) = V.at(i+1).at(j-1);
                    }

                    if (B_NW(flag.at(i).at(j)) == 1){
                        U.at(i-1).at(j) = 0;
                        U.at(i).at(j) = U.at(i).at(j+1);
                        V.at(i).at(j) = 0;
                        V.at(i).at(j-1) = V.at(i-1).at(j-1);
                    }

                    if (B_SE(flag.at(i).at(j)) == 1){
                        U.at(i).at(j)=0;
                        U.at(i-1).at(j) = U.at(i-1).at(j-1);
                        V.at(i).at(j-1)=0;
                        V.at(i).at(j) = V.at(i+1).at(j);
                    }

                    if (B_SW(flag.at(i).at(j)) == 1){
                        U.at(i-1).at(j) = 0;
                        U.at(i).at(j) = U.at(i).at(j-1);
                        V.at(i).at(j-1) = 0;
                        V.at(i).at(j) = V.at(i-1).at(j);
                    }
                    break;

                case flag_10bit::OUTFLOW: // assuming only x-direction
                    U.at(i).at(j) = U.at(i-1).at(j);
                    V.at(i).at(j) = V.at(i-1).at(j);
                    V.at(i).at(j-1) = V.at(i-1).at(j-1);
                    break;

                case flag_10bit::INFLOW: // assuming only x-direction
                    U.at(i).at(j)=1;
                    V.at(i).at(j) = 0;
                    V.at(i).at(j-1) = 0;
                    break;

                default:
                    break;
            }
        }
    }

    grid.set_velocity(U,velocity_type::U);
    grid.set_velocity(V,velocity_type::V);

}
