#include "helper.hpp"
#include "init.hpp"
#include <fstream>
#include <iostream>
int read_parameters( std::string szFileName,       /* name of the file */
                     double* Re,                /* reynolds number   */
                     double* UI,                /* velocity x-direction */
                     double* VI,                /* velocity y-direction */
                     double* PI,                /* pressure */
                     double* GX,                /* gravitation x-direction */
                     double* GY,                /* gravitation y-direction */
                     double* t_end,             /* end time */
                     double* xlength,           /* length of the domain x-dir.*/
                     double* ylength,           /* length of the domain y-dir.*/
                     double* dt,                /* time step */
                     double* dx,                /* length of a cell x-dir. */
                     double* dy,                /* length of a cell y-dir. */
                     int  *imax,                /* number of cells x-direction*/
                     int  *jmax,                /* number of cells y-direction*/
                     double *alpha,             /* uppwind differencing factor*/
                     double *omg,               /* relaxation factor */
                     double *tau,               /* safety factor for time step*/
                     int  *itermax,             /* max. number of iterations  */
        /* for pressure per time step */
                     double *eps,               /* accuracy bound for pressure*/
                     double *dt_value,           /* time for output */
                     std::string*  program,
                     std::string* geometry,
                     double* TI,                 /* initial temperature  */
                     double* Pr,                 /* prandtl number   */
                     double* beta,                /* the coefficient of thermal expansion */
                     double* heatTransfer,            /* Flag to enable heat transfer */
                     double* T_h,
                     double* T_c,
                     double* deltaP
)
{

    std::ifstream file(szFileName);
    if (!file.is_open()) return -1;
    std::string var;
    while (!file.eof() && file.good()) {
        file >> var;
        if (var[0] == '#') {     /* ignore comment line*/
            file.ignore(MAX_LINE_LENGTH, '\n');
        }
        else {
            if (var == "xlength")  file >> *xlength;
            if (var == "ylength")  file >> *ylength;
            if (var == "Re")       file >> *Re;
            if (var == "t_end")    file >> *t_end;
            if (var == "dt")       file >> *dt;
            if (var == "omg")      file >> *omg;
            if (var == "eps")      file >> *eps;
            if (var == "tau")      file >> *tau;
            if (var == "alpha")    file >> *alpha;
            if (var == "dt_value") file >> *dt_value;
            if (var == "UI")       file >> *UI;
            if (var == "VI")       file >> *VI;
            if (var == "GX")       file >> *GX;
            if (var == "GY")       file >> *GY;
            if (var == "PI")       file >> *PI;
            if (var == "itermax")  file >> *itermax;
            if (var == "imax")     file >> *imax;
            if (var == "jmax")     file >> *jmax;
            if( var == "program")  file >> *program;
            if( var == "geometry") file >> *geometry;
            if( var == "TI")       file >> *TI;
            if( var == "Pr")       file >> *Pr;
            if( var == "beta")     file >> *beta;
            if( var == "heatTransfer") file >> *heatTransfer;
            if( var == "T_h")      file >> *T_h;
            if( var == "T_c")      file >> *T_c;
            if( var == "deltaP")   file >> *deltaP;

        }
    }


    *dx = *xlength / (double)(*imax);
    *dy = *ylength / (double)(*jmax);


    if (!file.good() && !file.eof()) return -1;


    return 1;
}

void initCellInfo(const std::string &pgmFile, const int &imax, const int &jmax, matrix<unsigned int> &flag){

    // pgm reading
    const char *geometry_c = pgmFile.c_str();
    matrix<int> geo_pgm = read_pgm(geometry_c);

    // check if any forbidden cell exists
    check_forbidcell(imax, jmax, geo_pgm);

    // setting flags for cells
    for( size_t i = 0; i <= imax+1; i++){

        for(size_t j = 0; j <= jmax+1; j++){

            // fluid
            if(geo_pgm.at(i).at(j) == 4){flag.at(i).at(j) = flag_10bit::CELL_TYPE;}

            // no slip
            if(geo_pgm.at(i).at(j) == 0){flag.at(i).at(j) = flag_10bit::NO_SLIP;}

            // free-slip
            if(geo_pgm.at(i).at(j) == 1){flag.at(i).at(j) = flag_10bit::FREE_SLIP;}

            // out-flow
            if(geo_pgm.at(i).at(j) == 2){flag.at(i).at(j) = flag_10bit::OUTFLOW;}

            // in-flow
            if(geo_pgm.at(i).at(j) == 3){flag.at(i).at(j) = flag_10bit::INFLOW;}

            // dirichlet
            if(geo_pgm.at(i).at(j) == 5){flag.at(i).at(j) = flag_10bit::DIRICHLET; flag.at(i).at(j) |= flag_10bit::NO_SLIP;}



            // set boundaries for obstacle
            if( !IsFluid(geo_pgm.at(i).at(j))){

                // North Boundary B_N
                if (j<jmax+1){
                    if (geo_pgm.at(i).at(j+1) == 4){
                        flag.at(i).at(j) |= flag_10bit::N;
                    }
                }

                // South Boundary B_S
                if (j>0){
                    if(geo_pgm.at(i).at(j-1) == 4){
                        flag.at(i).at(j) |= flag_10bit::S;
                    }
                }

                // West Boundary B_W
                if (i>0){
                    if(geo_pgm.at(i-1).at(j) == 4){
                        flag.at(i).at(j) |= flag_10bit::W;
                    }
                }

                // East Boundary B_E
                if (i<imax+1){
                    if(geo_pgm.at(i+1).at(j) == 4){
                        flag.at(i).at(j) |= flag_10bit::E;
                    }
                }
            }

        }
    }

    std::ofstream file;
    std::string fileName = "flags.txt";
    file.open(fileName);

    for( size_t i = 0; i <= imax+1; i++) {

        for(size_t j = 0; j <= jmax+1; j++) {

            file << flag.at(i).at(j) << " " ;
        }
        file << std::endl;
    }
    file.close();
}

// check if the cell is fluid or not
bool IsFluid(unsigned int val){
    bool isfluid;
    isfluid = ((val == 4) || (val == 3) || (val == 2)) ? true:false;
    return isfluid;
}

// fluid on opposite sides Left and right
bool check_EW(int i, int j, matrix<int> geo_pgm)
{
    if( (geo_pgm.at(i-1).at(j)==4)&&(geo_pgm.at(i+1).at(j)==4) ) {return true;}
    else {return false;}
}

// fluid on opposite sides top and bottom
bool check_NS(int i, int j, matrix<int> geo_pgm)
{
    if( (geo_pgm.at(i).at(j+1)==4)&&(geo_pgm.at(i).at(j-1)==4) ) {return true;}
    else {return false;}
}

//Avoids any forbidden configuration
void check_forbidcell(int imax, int jmax, matrix<int> geo_pgm)
{
    for(int i=1; i<=imax; i++)
    {
        for(int j=1; j<=jmax; j++)
        {
            // inner cells
            if(geo_pgm.at(i).at(j)!=4)
            {
                if(check_EW(i,j,geo_pgm)||check_NS(i,j,geo_pgm))
                {
                    // TODO : throw exception
                    throw std::runtime_error(std::string("forbidden cell detected, boundary cell has more than two neighboring cell "));
                }
            }

            // left boundary
            if(geo_pgm.at(0).at(j)!=4)
            {
                if(check_NS(0,j,geo_pgm))
                {
                    // TODO : throw exception
                    throw std::runtime_error(std::string("forbidden cell detected, boundary cell has more than two neighboring cell "));

                }
            }

            // right boundary
            if(geo_pgm.at(imax+1).at(j)!=4)
            {
                if(check_NS(imax+1,j,geo_pgm))
                {
                    // TODO : throw exception
                    throw std::runtime_error(std::string("forbidden cell detected, boundary cell has more than two neighboring cell "));

                }
            }
        }
        // top boundary
        if(geo_pgm.at(i).at(jmax+1)!=4)
        {
            if(check_EW(i,jmax+1,geo_pgm))
            {
                // TODO : throw exception
                throw std::runtime_error(std::string("forbidden cell detected, boundary cell has more than two neighboring cell "));

            }
        }

        //bottom boundary
        if(geo_pgm.at(i).at(0)!=4)
        {
            if(check_EW(i,0,geo_pgm))
            {
                // TODO : throw exception
                throw std::runtime_error(std::string("forbidden cell detected, boundary cell has more than two neighboring cell "));

            }
        }
    }

}

