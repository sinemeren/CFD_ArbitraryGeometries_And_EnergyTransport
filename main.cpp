#include "helper.hpp"
#include "visual.hpp"
#include "init.hpp"
#include "sor.hpp"
#include <cstdio>
#include <iostream>
#include "uvp.hpp"
#include "boundary_val.hpp"

/**
 * The main operation reads the configuration file, initializes the scenario and
 * contains the main loop. So here are the individual steps of the algorithm:
 *
 * - DONE - read the program configuration file using read_parameters()
 * - DONE - set up the matrices (arrays) needed. Use the predefined matrix<typename> type and give initial values in the constructor.
 * - DONE - perform the main loop
 * - at the end: destroy any memory allocated and print some useful statistics
 *
 * The layout of the grid is decribed by the first figure below, the enumeration
 * of the whole grid is given by the second figure. All the unknowns corresond
 * to a two-dimensional degree of freedom layout, so they are not stored in
 * arrays, but in a matrix.
 *
 * @image html grid.jpg
 *
 * @image html whole-grid.jpg
 *
 * Within the main loop, the following steps are required (for some of the
 * operations, a definition is defined already within uvp.h):
 *
 * - DONE - calculate_dt() Determine the maximal time step size.
 * - DONE - boundaryvalues() Set the boundary values for the next time step.
 * - DONE - calculate_fg() Determine the values of F and G (diffusion and confection).
 *   This is the right hand side of the pressure equation and used later on for
 *   the time step transition.
 * - DONE - calculate_rs()
 * - DONE - Iterate the pressure poisson equation until the residual becomes smaller
 *   than eps or the maximal number of iterations is performed. Within the
 *   iteration loop the operation sor() is used.
 * - DONE - calculate_uv() Calculate the velocity at the next time step.
 */

void printProgressBar( double t, double t_end ){

    std::string bar;
    int percent = (t/t_end) * 100;

    for(int i = 0; i < 50; i++){
        if( i < (percent/2)){
            bar.replace(i,1,"=");
        }else if( i == (percent/2)){
            bar.replace(i,1,">");
        }else{
            bar.replace(i,1," ");
        }
    }


    std::cout<< "\r" "[" << bar << "] ";
    std::cout.width( 3 );
    std::cout<< percent << "%     \r" << std::flush;
}

int main(int argc, char* argv[]){



    // decleration of Variables
    double Re, UI, VI, PI, GX, GY, t_end, xlength, ylength, dt, dx, dy, alpha, omg, tau, eps, dt_value, TI, Pr, beta, heatTransfer, T_h, T_c, deltaP;
    int imax, jmax, itermax;
    std::string program;
    std::string geometry;


    std::string szFileName;

    if(argc < 2){

        throw std::runtime_error("please enter the example number");
    }
    else if(argc >3){

        throw std::runtime_error("too many input arguments, exampl usage: ./sim a");
    }
    else{

        if(*argv[1] == 'a'){

            szFileName = "PlaneShearFlow.dat";
        }
        else if(*argv[1] == 'b'){

            szFileName = "TheKarmanVortexStreeet.dat";
        }
        else if(*argv[1] == 'c'){

            szFileName = "flowOverStep.dat";
        }
        else if(*argv[1] == 'd'){

            szFileName = "NaturalConvection.dat";

        }
        else if(*argv[1] == 'g'){
            
            szFileName = "NaturalConvection2.dat";

        }
        else if(*argv[1] == 'e'){

            szFileName = "FluidTrap.dat";

        }
        else if(*argv[1] == 'f'){

            szFileName = "RayleighBenardConvection.dat";

        }
        else{
            throw std::runtime_error("examples are until f");
        }

    }

    std::cout << szFileName << " is running" << std::endl;

    // reading parameters
    read_parameters(szFileName, &Re, &UI, &VI, &PI, &GX, & GY, &t_end, &xlength, &ylength,
                    &dt, & dx, &dy, &imax, &jmax, &alpha, &omg, &tau, &itermax, &eps, &dt_value, &program, &geometry, &TI, &Pr, &beta, &heatTransfer, &T_h, &T_c, &deltaP);

    double val = deltaP * Re / dx;

    // setting up the matrices and assigning initial values
    matrix<double> F(imax+2, std::vector<double>(jmax+2, 0.0));
    matrix<double> G(imax+2, std::vector<double>(jmax+2, 0.0));
    matrix<double> RS(imax+2, std::vector<double>(jmax+2, 0.0));
    matrix<unsigned int> Flags(imax+2, std::vector<unsigned int>(jmax+2, 0));

    // setting up the initial grid
    Grid grid(imax, jmax, 1 , TI, PI, UI, VI);

    //initialize cell informations according to pgm file
    initCellInfo(geometry,imax,jmax, Flags);

    // initializing t and iteration variable for writing VTK file
    double t = 0.0;
    int iter_visual = 0;

    // performing iterations
    while ( t <= t_end){

        // calculation of maximal step size
        calculate_dt(Re, Pr, tau, &dt, dx, dy, imax, jmax, grid, heatTransfer);

        // setting the boundary values for temperature
        if(heatTransfer == 1.0){boundaryvaluesTemperature(imax, jmax, grid, Flags, T_h, T_c);}

        // setting the boundary values for velocity
        boundaryvalues_uv(imax, jmax, grid,Flags);

        // compute new temperature
        if(heatTransfer == 1.0){calculate_temp(dt, dx, dy, imax, jmax, Re, Pr, alpha, grid, Flags);}

        // calculation of F and G matrices
        calculate_fg(Re,GX, GY, alpha, dt, dx, dy, imax, jmax, grid, F, G, beta, heatTransfer, Flags);

        // calculation of right hand side
        calculate_rs(dt, dx, dy, imax, jmax, F, G, RS, Flags);

        // performing pressure poisson iteration
        int iter = 0;

        double res = 5.0;

        while ((iter <= itermax) && ( eps < res )){

        	if (iter == itermax){
            std::cout << " \r maximum iteration is reached! " << std::flush;
        	}

            sor(omg, dx, dy, imax, jmax, grid, RS, &res, Flags);
            iter = iter + 1;


        }

        // calculation of velocities
        calculate_uv(dt, dx, dy, imax, jmax, grid, F, G, heatTransfer, Flags);


        //update for time
        t = t + dt;

        // writing VTK file
        bool IsDoubleEqual = fabs(t - (iter_visual * dt_value)) < 0.0001;
        if( ( t!= 0.0 ) && ( ( t  > (iter_visual * dt_value) ) || (IsDoubleEqual) )){

            //VTKHelper::printVTKFile(grid, dx, dy, "t=" + std::to_string(t), "vtk", iter_visual);
            VTKHelper::printVTKFile(grid, dx, dy, "test", "vtk", iter_visual);

            iter_visual = iter_visual + 1;

        }

        // call for progress bar
        printProgressBar( t, t_end);

    }

    return -1;
}
