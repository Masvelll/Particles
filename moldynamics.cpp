#include <pybind11/pybind11.h>
#include <iostream>

namespace py = pybind11;

struct Mol{
    Mol(const double &x, const double &y, const double &z) : x(x), y(y), z(z){
        vx = 0; vy = 0; vz = 0;
        ax = 0; ay = 0; az = 0;
        ax_prev = 0; ay_prev = 0; az_prev = 0;
    }

    double x, y, z;
    double vx, vy, vz;
    double ax, ay, az;
    double ax_prev, ay_prev, az_prev;
};

void show_pos(Mol mol){
    std::cout << "x: " << mol.x
              << " y: " << mol.y
              << " z: " << mol.z << "\n";
};

PYBIND11_MODULE(moldynamics, m) {
    m.doc() = "pybind plugin for molecular dynamics";

    py::class_<Mol>(m, "Mol")
        .def(py::init<const double &, const double &, const double &>());
    
    m.def("show_pos", &show_pos, "A function that shows position of a molecule");
        
}


// to make a file
// c++ -O3 -Wall -shared -std=c++11 -fPIC $(python3-config --includes) -I extern/pybind11/include moldynamics.cpp -o moldynamics$(python3-config --extension-suffix)