#ifndef __RANDWERTE_HPP__
#define __RANDWERTE_HPP__

#include "cstring"
#include "helper.hpp"
#include "datastructures.hpp"
#include "grid.hpp"

/**
 * The boundary values of the problem are set.
 */
void boundaryvaluesTemperature(int imax, int jmax, Grid& grid, const matrix<unsigned int> &flag, double T_h, double T_c);

void boundaryvalues_uv(int imax,int jmax, Grid& grid,const matrix<unsigned int> &flag);

int B_E(unsigned int flag);

int B_W(unsigned int flag);

int B_N(unsigned int flag);

int B_S(unsigned int flag);

int B_NE(unsigned int flag);

int B_NW(unsigned int flag);

int B_SE(unsigned int flag);

int B_SW(unsigned int flag);


#endif
