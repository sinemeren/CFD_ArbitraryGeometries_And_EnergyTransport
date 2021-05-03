#include "vector"

#ifndef DATA_STRUCTURES_H
#define DATA_STRUCTURES_H

template<typename T>
using matrix = std::vector<std::vector<T>>;


enum flag_10bit{

    CELL_TYPE = 1,  // 1<<0
    NO_SLIP = 2,	// 1<<1
    FREE_SLIP = 4,	// 1<<2
    OUTFLOW = 8,	// 1<<3
    INFLOW = 16,	// 1<<4
    N = 32,			// 1<<5
    S = 64,			// 1<<6
    W = 128,		// 1<<7
    E = 256,		// 1<<8
    DIRICHLET = 512  // 1<<9
};


#endif //DATA_STRUCTURES_H
