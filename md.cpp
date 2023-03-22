#include <iostream>
#include <vector>
#include <cmath>

const double epsilon = 1.0;
const double sigma = 1.0;
const double dt = 0.01;
const double mass = 1.0;
const double size = 10;
const int num_of_particles = 10;

struct Pos {
    double x, y, z;
};

struct Force {
    double fx, fy, fz;
};

void show(Pos const &p){
    std::cout << "x = " << p.x << " y = " << p.y << " z = " << p.z << "\n";
};

void show(Force const &p){
    std::cout << "fx = " << p.fx << " fy = " << p.fy << " fz = " << p.fz << "\n";
};

Force calc_forces(Pos const &p1, Pos const &p2){
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
    Force f;
    f.fx = fx;
    f.fy = fy;
    f.fz = fz;
    return f;
}

int main() {

    std::vector<Pos> positions;
    std::vector<Force> forces;
    for (double i; i < num_of_particles; i++) {
        Pos p;
        p.x = i; p.y = i; p.z = i;
        positions.push_back(p);
    }
    
    for (int8_t i; i < num_of_particles; i++){
        show(positions[i]);
    }

    Force f = calc_forces(positions[0], positions[1]);
    show(f);

    

    
    return 0;
}