#include "sor.hpp"
#include <cmath>
#include "boundary_val.hpp"
void sor(
        double omg,
        double dx,
        double dy,
        int    imax,
        int    jmax,
        Grid& grid,
        matrix<double> &RS,
        double *res,
        const matrix<unsigned int> flag
) {
   // int i, j;
    double rloc;
    double coeff = omg / (2.0 * (1.0 / (dx * dx) + 1.0 / (dy * dy)));
    matrix<double> P;
    grid.pressure(P);

    /* SOR iteration */
    for(int i = 1; i <= imax; i++) {

        for(int j = 1; j<=jmax; j++) {

            if(((flag.at(i).at(j) & flag_10bit::CELL_TYPE) == flag_10bit::CELL_TYPE)){
                P.at(i).at(j) = (1.0-omg)*P.at(i).at(j)
                                + coeff*(( P.at(i+1).at(j)+P.at(i-1).at(j))/(dx*dx) + ( P.at(i).at(j+1)+P.at(i).at(j-1))/(dy*dy) - RS.at(i).at(j));
            }
        }
    }

    /* compute the residual */
    rloc = 0;
    int num_fluid_elem=0;
    for(int i = 1; i <= imax; i++) {
        for(int j = 1; j <= jmax; j++) {

            if(((flag.at(i).at(j) & flag_10bit::CELL_TYPE) == flag_10bit::CELL_TYPE)){
                num_fluid_elem++;
                rloc += ( (P.at(i+1).at(j)-2.0*P.at(i).at(j)+P.at(i-1).at(j))/(dx*dx) + ( P.at(i).at(j+1)-2.0*P.at(i).at(j)+P.at(i).at(j-1))/(dy*dy) - RS.at(i).at(j))*
                        ( (P.at(i+1).at(j)-2.0*P.at(i).at(j)+P.at(i-1).at(j))/(dx*dx) + ( P.at(i).at(j+1)-2.0*P.at(i).at(j)+P.at(i).at(j-1))/(dy*dy) - RS.at(i).at(j));
            }
        }
    }

    //rloc = rloc/(imax*jmax);
    rloc = rloc/(num_fluid_elem);
    rloc = sqrt(rloc);

    /* set residual */
    *res = rloc;
    num_fluid_elem = 0;


    // --------------------------------------------- BOUNDARY CONDITIONS ---------------------------------------------

    // Top boundary
    for (int i = 1; i <= imax; ++i) {

    	switch(flag.at(i).at(jmax+1) & (flag_10bit::NO_SLIP | flag_10bit::FREE_SLIP | flag_10bit::OUTFLOW | flag_10bit::INFLOW)){
		case flag_10bit::NO_SLIP: //No slip conditions
		case flag_10bit::FREE_SLIP:
		case flag_10bit::OUTFLOW:
		case flag_10bit::INFLOW:

		if (B_E(flag[i][jmax+1]) == 1) {P[i][jmax+1] = P[i + 1][jmax+1];}

		if (B_W(flag[i][jmax+1]) == 1) {P[i][jmax+1] = P[i - 1][jmax+1];}

		if (B_N(flag[i][jmax+1]) == 1) {throw std::runtime_error(std::string("impossible to have B_N at top boundary in sor"));}

		if (B_S(flag[i][jmax+1]) == 1) {P[i][jmax+1] = P[i][jmax];}

		if (B_NE(flag[i][jmax+1]) == 1) {throw std::runtime_error(std::string("impossible to have B_NE at top boundary in sor"));}

		if (B_NW(flag[i][jmax+1]) == 1) {throw std::runtime_error(std::string("impossible to have B_NW at top boundary in sor"));}

		if (B_SE(flag[i][jmax+1]) == 1) {P[i][jmax+1] = (P[i][jmax] + P[i + 1][jmax+1]) / 2;}

		if (B_SW(flag[i][jmax+1]) == 1) {P[i][jmax+1] = (P[i][jmax] + P[i - 1][jmax+1]) / 2;}

		/*
		if (flag[i][j] & (1 << 3)) {P[i][j] = 0;}

		if (flag[i][j] & (1 << 4)) {P[i][j] = P[i + 1][j];}
		*/

		break;

		default:
			break;

    	}
    }

    // Bottom boundary
    for (int i = 1; i <= imax; ++i) {

		switch(flag.at(i).at(0) & (flag_10bit::NO_SLIP | flag_10bit::FREE_SLIP | flag_10bit::OUTFLOW | flag_10bit::INFLOW)){
		case flag_10bit::NO_SLIP: //No slip conditions
		case flag_10bit::FREE_SLIP:
		case flag_10bit::OUTFLOW:
		case flag_10bit::INFLOW:

		if (B_E(flag[i][0]) == 1) {P[i][0] = P[i + 1][0];}

		if (B_W(flag[i][0]) == 1) {P[i][0] = P[i - 1][0];}

		if (B_N(flag[i][0]) == 1) {P[i][0] = P[i][1];}

		if (B_S(flag[i][0]) == 1) {throw std::runtime_error(std::string("impossible to have B_S at bottom boundary in sor"));}

		if (B_NE(flag[i][0]) == 1) {P[i][0] = (P[i][1] + P[i + 1][0]) / 2;}

		if (B_NW(flag[i][0]) == 1) {P[i][0] = (P[i][1] + P[i - 1][0]) / 2;}

		if (B_SE(flag[i][0]) == 1) {throw std::runtime_error(std::string("impossible to have B_SE at bottom boundary in sor"));}

		if (B_SW(flag[i][0]) == 1) {throw std::runtime_error(std::string("impossible to have B_SW at bottom boundary in sor"));}

		/*
		if (flag[i][j] & (1 << 3)) P[i][j] = 0;

		if (flag[i][j] & (1 << 4)) P[i][j] = P[i + 1][j];
		*/

		break;

		default:
			break;


		}
    }

    // Left boundary
    for (int j = 1; j <= jmax; ++j) {

    	switch(flag.at(0).at(j) & (flag_10bit::NO_SLIP | flag_10bit::FREE_SLIP | flag_10bit::OUTFLOW | flag_10bit::INFLOW)){
				case flag_10bit::NO_SLIP: //No slip conditions
				case flag_10bit::FREE_SLIP:
				case flag_10bit::OUTFLOW:
				case flag_10bit::INFLOW:

				if (B_E(flag[0][j]) == 1) {P[0][j] = P[1][j];}

				if (B_W(flag[0][j]) == 1) {throw std::runtime_error(std::string("impossible to have B_W at left boundary in sor"));}

				if (B_N(flag[0][j]) == 1) {P[0][j] = P[0][j + 1];}

				if (B_S(flag[0][j]) == 1) {P[0][j] = P[0][j - 1];}

				if (B_NE(flag[0][j]) == 1) {P[0][j] = (P[0][j + 1] + P[1][j]) / 2;}

				if (B_NW(flag[0][j]) == 1) {throw std::runtime_error(std::string("impossible to have B_NW at left boundary in sor"));}

				if (B_SE(flag[0][j]) == 1) {P[0][j] = (P[0][j - 1] + P[1][j]) / 2;}

				if (B_SW(flag[0][j]) == 1) {throw std::runtime_error(std::string("impossible to have B_SW at left boundary in sor"));}

				/*
				if (flag[i][j] & (1 << 3)) P[i][j] = 0;

				if (flag[i][j] & (1 << 4)) P[i][j] = P[i + 1][j];
				*/

				break;

				default:
					break;



    	}
    }

    // Right Boundary

    for (int j = 1; j <= jmax; ++j) {

        	switch(flag.at(imax+1).at(j) & (flag_10bit::NO_SLIP | flag_10bit::FREE_SLIP | flag_10bit::OUTFLOW | flag_10bit::INFLOW)){
			case flag_10bit::NO_SLIP: //No slip conditions
			case flag_10bit::FREE_SLIP:
			case flag_10bit::OUTFLOW:
			case flag_10bit::INFLOW:


			if (B_E(flag[imax+1][j]) == 1) {throw std::runtime_error(std::string("impossible to have B_E at right boundary in sor"));}

			if (B_W(flag[imax+1][j]) == 1) {P[imax+1][j] = P[imax][j];}

			if (B_N(flag[imax+1][j]) == 1) {P[imax+1][j] = P[imax+1][j + 1];}

			if (B_S(flag[imax+1][j]) == 1) {P[imax+1][j] = P[imax+1][j - 1];}

			if (B_NE(flag[imax+1][j]) == 1) {throw std::runtime_error(std::string("impossible to have B_E at right boundary in sor"));}

			if (B_NW(flag[imax+1][j]) == 1) {P[imax+1][j] = (P[imax+1][j + 1] + P[imax][j]) / 2;}

			if (B_SE(flag[imax+1][j]) == 1) {throw std::runtime_error(std::string("impossible to have B_E at right boundary in sor"));}

			if (B_SW(flag[imax+1][j]) == 1) {P[imax+1][j] = (P[imax+1][j - 1] + P[imax][j]) / 2;}

			/*
			if (flag[i][j] & (1 << 3)) P[i][j] = 0;

			if (flag[i][j] & (1 << 4)) P[i][j] = P[i + 1][j];
			*/

			break;

			default:
				break;
        }
    }

    // Top - Left corner
    switch(flag.at(0).at(jmax+1) & (flag_10bit::NO_SLIP | flag_10bit::FREE_SLIP | flag_10bit::OUTFLOW | flag_10bit::INFLOW)){
	case flag_10bit::NO_SLIP: //No slip conditions
	case flag_10bit::FREE_SLIP:
	case flag_10bit::OUTFLOW:
	case flag_10bit::INFLOW:

	if (B_E(flag[0][jmax+1]) == 1) {P[0][jmax+1] = P[1][jmax+1];}

	if (B_W(flag[0][jmax+1]) == 1) {throw std::runtime_error(std::string("impossible to have B_W at top-left boundary in sor"));}

	if (B_N(flag[0][jmax+1]) == 1) {throw std::runtime_error(std::string("impossible to have B_N at top-left boundary in sor"));}

	if (B_S(flag[0][jmax+1]) == 1) {P[0][jmax+1] = P[0][jmax];}

	if (B_NE(flag[0][jmax+1]) == 1) {throw std::runtime_error(std::string("impossible to have B_NE at top-left boundary in sor"));}

	if (B_NW(flag[0][jmax+1]) == 1) {throw std::runtime_error(std::string("impossible to have B_NW at top-left boundary in sor"));}

	if (B_SE(flag[0][jmax+1]) == 1) {P[0][jmax+1] = (P[0][jmax] + P[1][jmax+1]) / 2;}

	if (B_SW(flag[0][jmax+1]) == 1) {throw std::runtime_error(std::string("impossible to have B_SW at top-left boundary in sor"));}

	/*
	if (flag[i][j] & (1 << 3)) P[i][j] = 0;

	if (flag[i][j] & (1 << 4)) P[i][j] = P[i + 1][j];
	*/

	break;

	default:
		break;
	}

    // Top-right boundary
    switch(flag.at(imax+1).at(jmax+1) & (flag_10bit::NO_SLIP | flag_10bit::FREE_SLIP | flag_10bit::OUTFLOW | flag_10bit::INFLOW)){
	case flag_10bit::NO_SLIP: //No slip conditions
	case flag_10bit::FREE_SLIP:
	case flag_10bit::OUTFLOW:
	case flag_10bit::INFLOW:

	if (B_E(flag[imax+1][jmax+1]) == 1) {throw std::runtime_error(std::string("impossible to have B_E at top-right boundary in sor"));}

	if (B_W(flag[imax+1][jmax+1]) == 1) {P[imax+1][jmax+1] = P[imax][jmax+1];}

	if (B_N(flag[imax+1][jmax+1]) == 1) {throw std::runtime_error(std::string("impossible to have B_N at top-right boundary in sor"));}

	if (B_S(flag[imax+1][jmax+1]) == 1) {P[imax+1][jmax+1] = P[imax+1][jmax];}

	if (B_NE(flag[imax+1][jmax+1]) == 1) {throw std::runtime_error(std::string("impossible to have B_NE at top-right boundary in sor"));}

	if (B_NW(flag[imax+1][jmax+1]) == 1) {throw std::runtime_error(std::string("impossible to have B_NW at top-right boundary in sor"));}

	if (B_SE(flag[imax+1][jmax+1]) == 1) {throw std::runtime_error(std::string("impossible to have B_SE at top-right boundary in sor"));}

	if (B_SW(flag[imax+1][jmax+1]) == 1) {P[imax+1][jmax+1] = (P[imax+1][jmax] + P[imax][jmax+1]) / 2;}

	/*
	if (flag[i][j] & (1 << 3)) P[i][j] = 0;

	if (flag[i][j] & (1 << 4)) P[i][j] = P[i + 1][j];
	*/

	break;

	default:
		break;
	}

    // Bottom - left corner
    switch(flag.at(0).at(0) & (flag_10bit::NO_SLIP | flag_10bit::FREE_SLIP | flag_10bit::OUTFLOW | flag_10bit::INFLOW)){
	case flag_10bit::NO_SLIP: //No slip conditions
	case flag_10bit::FREE_SLIP:
	case flag_10bit::OUTFLOW:
	case flag_10bit::INFLOW:

	if (B_E(flag[0][0]) == 1) {P[0][0] = P[1][0];}

	if (B_W(flag[0][0]) == 1) {throw std::runtime_error(std::string("impossible to have B_W at bottom-left boundary in sor"));}

	if (B_N(flag[0][0]) == 1) {P[0][0] = P[0][1];}

	if (B_S(flag[0][0]) == 1) {throw std::runtime_error(std::string("impossible to have B_S at bottom-left boundary in sor"));}

	if (B_NE(flag[0][0]) == 1) {P[0][0] = (P[0][1] + P[1][0]) / 2;}

	if (B_NW(flag[0][0]) == 1) {throw std::runtime_error(std::string("impossible to have B_NW at bottom-left boundary in sor"));}

	if (B_SE(flag[0][0]) == 1) {throw std::runtime_error(std::string("impossible to have B_SE at bottom-left boundary in sor"));}

	if (B_SW(flag[0][0]) == 1) {throw std::runtime_error(std::string("impossible to have B_SW at bottom-left boundary in sor"));}

	/*
	if (flag[i][j] & (1 << 3)) P[i][j] = 0;

	if (flag[i][j] & (1 << 4)) P[i][j] = P[i + 1][j];
	*/

	break;

	default:
		break;
	}

    // Bottom-right corner
    switch(flag.at(imax+1).at(0) & (flag_10bit::NO_SLIP | flag_10bit::FREE_SLIP | flag_10bit::OUTFLOW | flag_10bit::INFLOW)){
	case flag_10bit::NO_SLIP: //No slip conditions
	case flag_10bit::FREE_SLIP:
	case flag_10bit::OUTFLOW:
	case flag_10bit::INFLOW:

	if (B_E(flag[imax+1][0]) == 1) {throw std::runtime_error(std::string("impossible to have B_E at bottom-right boundary in sor"));}

	if (B_W(flag[imax+1][0]) == 1) {P[imax+1][0] = P[imax][0];}

	if (B_N(flag[imax+1][0]) == 1) {P[imax+1][0] = P[imax+1][1];}

	if (B_S(flag[imax+1][0]) == 1) {throw std::runtime_error(std::string("impossible to have B_S at bottom-right boundary in sor"));}

	if (B_NE(flag[imax+1][0]) == 1) {throw std::runtime_error(std::string("impossible to have B_NE at bottom-right boundary in sor"));}

	if (B_NW(flag[imax+1][0]) == 1) {P[imax+1][0] = (P[imax+1][1] + P[imax][0]) / 2;}

	if (B_SE(flag[imax+1][0]) == 1) {throw std::runtime_error(std::string("impossible to have B_SE at bottom-right boundary in sor"));}

	if (B_SW(flag[imax+1][0]) == 1) {throw std::runtime_error(std::string("impossible to have B_SW at bottom-right boundary in sor"));}

	/*
	if (flag[i][j] & (1 << 3)) P[i][j] = 0;

	if (flag[i][j] & (1 << 4)) P[i][j] = P[i + 1][j];
	*/

	break;

	default:
		break;
	}

    // inner cells
    for (int i = 1; i <= imax; ++i) {

        for (int j = 1; j <= jmax; ++j) {

        	switch(flag.at(i).at(j) & (flag_10bit::NO_SLIP | flag_10bit::FREE_SLIP | flag_10bit::OUTFLOW | flag_10bit::INFLOW)){
        	case flag_10bit::NO_SLIP: //No slip conditions
        	case flag_10bit::FREE_SLIP:
        	case flag_10bit::OUTFLOW:
        	case flag_10bit::INFLOW:

            if (B_E(flag[i][j]) == 1) {P[i][j] = P[i + 1][j];}

            if (B_W(flag[i][j]) == 1) {P[i][j] = P[i - 1][j];}

            if (B_N(flag[i][j]) == 1) {P[i][j] = P[i][j + 1];}

            if (B_S(flag[i][j]) == 1) {P[i][j] = P[i][j - 1];}

            if (B_NE(flag[i][j]) == 1) {P[i][j] = (P[i][j + 1] + P[i + 1][j]) / 2;}

            if (B_NW(flag[i][j]) == 1) {P[i][j] = (P[i][j + 1] + P[i - 1][j]) / 2;}

            if (B_SE(flag[i][j]) == 1) {P[i][j] = (P[i][j - 1] + P[i + 1][j]) / 2;}

            if (B_SW(flag[i][j]) == 1) {P[i][j] = (P[i][j - 1] + P[i - 1][j]) / 2;}

            /*
            if (flag[i][j] & (1 << 3)) P[i][j] = 0;

            if (flag[i][j] & (1 << 4)) P[i][j] = P[i + 1][j];
            */

            break;

        	default:
        		break;
        	}
        }
    }

    grid.set_pressure(P);

}

