#include <pybind11/pybind11.h>
#include <iostream>
#include <vector>
#include <pybind11/stl.h>

namespace py = pybind11;

struct Mol{
    Mol(const double &x, const double &y, const double &z,
        const double &mass) : mass(mass), x(x), y(y), z(z) {
        vx = 0.0; vy = 0.0; vz = 0.0;
        ax = 0.0; ay = 0.0; az = 0.0;
        ax_prev = 0.0; ay_prev = 0.0; az_prev = 0.0;
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
    //Box(const size_t &size, std::vector<Mol> &mols) : size(size), mols(mols), full_amount(0) {}
    Box(const size_t &size) : size(size), full_amount(0), kinetic_energy(0), potential_energy(0), dt(0.01), t(0) {};

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

    void calc_forces(Mol &p1, Mol &p2){
        double epsilon = 1.0;
        double sigma = 1.0;

        double rx = p1.x - p2.x;
        double ry = p1.y - p2.y;
        double rz = p1.z - p2.z;

        if (abs(rx) > size / 2) {rx -= size * (round(rx / size));}
        if (abs(ry) > size / 2) {ry -= size * (round(ry / size));}
        if (abs(rz) > size / 2) {rz -= size * (round(rz / size));}

        double r2 = rx * rx + ry * ry + rz * rz;
        double force_r = - 24 * epsilon * ((pow(sigma, 6) / pow(r2, 3) - 2 * (pow(sigma, 12) / pow(r2, 6))));
        
        //std::cout << "r2 = " <<  r2 << "\n";

        double fx = force_r * rx / r2;
        double fy = force_r * ry / r2;
        double fz = force_r * rz / r2;

        p1.ax = fx / p1.mass; p2.ax = - fx / p2.mass;
        p1.ay = fy / p1.mass; p2.ay = - fy / p2.mass;
        p1.az = fz / p1.mass; p2.az = - fz / p2.mass;

        //std::cout << p1.az << " ";
    };

    double calc_distance(const Mol &p1, const Mol &p2){
        double rx = p1.x - p2.x;
        double ry = p1.y - p2.y;
        double rz = p1.z - p2.z;

        if (abs(rx) > size / 2) {rx -= size * (round(rx / size));}
        if (abs(ry) > size / 2) {ry -= size * (round(ry / size));}
        if (abs(rz) > size / 2) {rz -= size * (round(rz / size));}

        double r2 = rx * rx + ry * ry + rz * rz;
        return r2;
    };

    void calc_energy(const std::vector<Mol> &mols){
        double epsilon = 1.0;
        double sigma = 1.0;
        double kinetic = 0.0, potential = 0.0;
        for (int16_t i = 0; i < full_amount - 1; i++){
            for (int16_t j = i + 1; j < full_amount; j++){
                double r2 = calc_distance(mols[i], mols[j]);
                //std::cout << i << " " << j << " r2: " << r2 << "\n";
                potential += 4 * epsilon * (pow(sigma, 12) / pow(r2, 6) - pow(sigma, 6) / pow(r2, 3));
            };
            kinetic += mols[0].mass / 2 * (pow(mols[i].vx, 2) + 
                                pow(mols[i].vy, 2) +
                                pow(mols[i].vz, 2));
        };
        kinetic += mols[0].mass / 2 * (pow(mols[full_amount - 1].vx, 2) + 
                            pow(mols[full_amount - 1].vy, 2) +
                            pow(mols[full_amount - 1].vz, 2));
        std::cout << "Potential: " << potential << "\n";
        std::cout << "Kinetic: " << kinetic << "\n";
        kinetic_energy = kinetic;
        potential_energy = potential;
    };

    void update_force(std::vector<Mol> &mols, int16_t i){
        for (int16_t j = 0; i < full_amount; i++){
            if (i != j){
                calc_forces(mols[i], mols[j]);

            };
            
        };

    };

    void apply_boundary(Mol &p){
        if (p.x > size){ p.x -= size; }
        if (p.x < 0) {p.x += size;}
        if (p.y > size){ p.y -= size; }
        if (p.y < 0) {p.y += size;}
        if (p.z > size){ p.z -= size; }
        if (p.z < 0) {p.z += size;}
    };

    void update_position(std::vector<Mol> &mols){
        for (int16_t i = 0; i < full_amount; i++){
            mols[i].x = mols[i].x + mols[i].vx * dt 
            + 1 / 2 * mols[i].ax_prev * dt * dt;
            mols[i].y = mols[i].y + mols[i].vy * dt 
            + 1 / 2 * mols[i].ay_prev * dt * dt;
            mols[i].z = mols[i].z + mols[i].vz * dt 
            + 1 / 2 * mols[i].az_prev * dt * dt;

            std::cout << "Before az: " << mols[i].az_prev << ' ' << mols[i].az << "\n";
            std::cout << "z: " << mols[i].z << " vz: " << mols[i].vz << "\n";
            mols[i].ax_prev = mols[i].ax;
            mols[i].ay_prev = mols[i].ay;
            mols[i].az_prev = mols[i].az;

            apply_boundary(mols[i]);

        };

        for (int16_t i = 0; i < full_amount; i++){
            update_force(mols, i);
            //std::cout << "After: " << 1.0 / 2.0 * (mols[i].ax_prev + mols[i].ax) << "\n";
            mols[i].vx += 1.0 / 2.0 * (mols[i].ax_prev + mols[i].ax) * dt;
            mols[i].vy += 1.0 / 2.0 * (mols[i].ax_prev + mols[i].ax) * dt;
            mols[i].vz += 1.0 / 2.0 * (mols[i].az_prev + mols[i].az) * dt;
        };
        
    };

    void update(){
        t += dt;
        calc_energy(mols);
        update_position(mols);
    };


private:
    size_t size;
    int64_t full_amount;
    double kinetic_energy, potential_energy, dt, t;
    std::vector<Mol> mols;
    

      
};



PYBIND11_MODULE(moldynamics, m) {
    m.doc() = "pybind plugin for molecular dynamics";

    m.attr("epsilon") = 1;
    m.attr("sigma") = 1;

    py::class_<Mol>(m, "Mol")
        .def(py::init<double &, double &, double &, const double &>());
    
    py::class_<Box>(m, "Box")
        //.def(py::init<const size_t &, std::vector<Mol> &>())
        .def(py::init<const size_t &>())
        .def("GetMols", &Box::GetMols)
        .def("GenerateMols", &Box::GenerateMols)
        .def("update", &Box::update);


    m.def("show_pos", &show_pos, "A function that shows position of a molecule");
        
}


// to make a file
// c++ -O3 -Wall -shared -std=c++11 -fPIC $(python3-config --includes) -I extern/pybind11/include moldynamics.cpp -o moldynamics$(python3-config --extension-suffix)