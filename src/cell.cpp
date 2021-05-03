#include "cell.hpp"

Cell::Cell() {};

Cell::Cell(double& TI,double& PI, double& UI, double& VI):
_pressure(PI),_temperature(TI) {
    set_velocity(UI, velocity_type::U);
    set_velocity(VI, velocity_type::V);
};

// Pressure Get and Set
double& Cell::pressure() {
    return _pressure;
}

void Cell::set_pressure(double& value) {
    _pressure = value;
}


// Temperature Get and Set
double& Cell::temperature() {
    return _temperature;
}

void Cell::set_temperature(double& value) {
    _temperature = value;
}


// Velocity Get and Set
double& Cell::velocity(velocity_type type) {
    return _velocity[static_cast<int>(type)];
};

void Cell::set_velocity(double& value, velocity_type type) {
    _velocity[static_cast<int>(type)] = value;
};

// borders Get and Set
bool& Cell::border(border_position position) {
    return _border[static_cast<int>(position)];
};

// Set border
void Cell::set_border(border_position position) {
    _border[static_cast<int>(position)] = true;
};
