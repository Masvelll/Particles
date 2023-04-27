#include <pybind11/pybind11.h>
#include <iostream>
#include <vector>
#include <pybind11/stl.h>

namespace py = pybind11;

struct Mol{
    Mol(const double &x, const double &y, const double &z,
        const double &mass) : x(x), y(y), z(z), mass(mass){
        vx = 0; vy = 0; vz = 0;
        ax = 0; ay = 0; az = 0;
        ax_prev = 0; ay_prev = 0; az_prev = 0;
    }
    double mass;
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

class Box{
public:
    Box(const size_t &size, std::vector<Mol> &mols) : size(size), mols(mols), full_amount(0) {}
    Box(const size_t &size) : size(size), full_amount(0) { };

    void GenerateMols(const int64_t &amount, const double &mass){
        full_amount += amount;
        int64_t axis_step = pow(size, 1.0/3);
        for (double x = 0.0; x < size; x+= axis_step) {
            for (double y = 0.0; y < size; y+= axis_step) {
                for (double z = 0.0; z < size; z+= axis_step) {
                    Mol m = Mol(x, y, z, mass);
                    mols.push_back(m);
                };
            };
        };
    }

    std::vector<std::vector<double>> GetMols() {
        std::vector<std::vector<double>> poses;

        for (auto i = 0; i < full_amount; i++){
            Mol m = mols[i];   
            std::vector<double> pose;
            pose.push_back(m.x);
            pose.push_back(m.y);
            pose.push_back(m.z);
            
            poses.push_back(pose);
            }
            
        

        return poses;
    }

private:
    size_t size;
    int64_t full_amount;
    std::vector<Mol> mols;

      
};



/*void calc_forces(Mol &p1, Mol &p2){
    double rx = p1.x - p2.x;
    double ry = p1.y - p2.y;
    double rz = p1.z - p2.z;

    if (abs(rx) > size / 2) {rx -= size * (round(rx / size));}
    if (abs(ry) > size / 2) {ry -= size * (round(ry / size));}
    if (abs(rz) > size / 2) {rz -= size * (round(rz / size));}

    double r2 = rx * rx + ry * ry + rz * rz;
    double force_r = - 24 * epsilon * ((pow(sigma, 6) / pow(r2, 3) - 2 * (pow(sigma, 12) / pow(r2, 6))));
    
    double fx = force_r * rx / r2;
    double fy = force_r * ry / r2;
    double fz = force_r * rz / r2;

    p1.ax = fx / p1.mass; p2.ax = - fx / p2.mass;
    p1.ay = fy / p1.mass; p2.ay = - fy / p2.mass;
    p1.az = fz / p1.mass; p2.az = - fz / p2.mass;
};*/



PYBIND11_MODULE(moldynamics, m) {
    m.doc() = "pybind plugin for molecular dynamics";

    m.attr("epsilon") = 1;
    m.attr("sigma") = 1;

    py::class_<Mol>(m, "Mol")
        .def(py::init<const double &, const double &, const double &, const double &>());
    
    py::class_<Box>(m, "Box")
        .def(py::init<const size_t &, std::vector<Mol> &>())
        .def(py::init<const size_t &>())
        .def("GetMols", &Box::GetMols)
        .def("GenerateMols", &Box::GenerateMols);

    m.def("show_pos", &show_pos, "A function that shows position of a molecule");
        
}


// to make a file
// c++ -O3 -Wall -shared -std=c++11 -fPIC $(python3-config --includes) -I extern/pybind11/include moldynamics.cpp -o moldynamics$(python3-config --extension-suffix)