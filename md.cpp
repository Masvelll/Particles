#include <iostream>
#include <vector>
#include <cmath>

const double epsilon = 1.0;
const double sigma = 1.0;
const double dt = 0.001;
const double mass = 1.0;
const double size = 5;
const int num_of_particles = 10;
double t = 0.0;
double kinetic_energy, potential_energy = 0;

struct Particle {
    double x, y, z;
    double vx = 0.0, vy = 0.0, vz = 0.0;
    double ax = 0.0, ay = 0.0, az = 0.0;
    double ax_prev = 0.0, ay_prev = 0.0, az_prev = 0.0;
};

void show(Particle const &p){
    std::cout << "x = " << p.x << " y = " << p.y << " z = " << p.z << "\n";
    std::cout << "vx = " << p.vx << " vy = " << p.vy << " vz = " << p.vz << "\n";
    std::cout << "ax = " << p.ax << " ay = " << p.ay << " az = " << p.az << "\n";
    
};



void calc_forces(Particle &p1, Particle &p2){
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

    p1.ax = fx / mass; p2.ax = - fx / mass;
    p1.ay = fy / mass; p2.ay = - fy / mass;
    p1.az = fz / mass; p2.az = - fz / mass;
}

double calc_distance(const Particle &p1, const Particle &p2){
    double rx = p1.x - p2.x;
    double ry = p1.y - p2.y;
    double rz = p1.z - p2.z;

    if (abs(rx) > size / 2) {rx -= size * (round(rx / size));}
    if (abs(ry) > size / 2) {ry -= size * (round(ry / size));}
    if (abs(rz) > size / 2) {rz -= size * (round(rz / size));}

    double r2 = rx * rx + ry * ry + rz * rz;
    return r2;
};

void calc_energy(const std::vector<Particle> &particles, double &kinetic_energy, double &potential_energy){
    double kinetic = 0, potential = 0;
    for (int16_t i = 0; i < num_of_particles - 1; i++){
        for (int16_t j = i + 1; j < num_of_particles; j++){
            double r2 = calc_distance(particles[i], particles[j]);
            //std::cout << i << " " << j << " r2: " << r2 << "\n";
            potential += 4 * epsilon * (pow(sigma, 12) / pow(r2, 6) - pow(sigma, 6) / pow(r2, 3));
        };
        kinetic += mass / 2 * (pow(particles[i].vx, 2) + 
                               pow(particles[i].vy, 2) +
                               pow(particles[i].vz, 2));
    };
    kinetic += mass / 2 * (pow(particles[num_of_particles - 1].vx, 2) + 
                           pow(particles[num_of_particles - 1].vy, 2) +
                           pow(particles[num_of_particles - 1].vz, 2));
    std::cout << "Potential: " << potential << "\n";
    kinetic_energy = kinetic;
    potential_energy = potential;
};

void update_force(std::vector<Particle> &particles, int16_t i){
    for (int16_t j = 0; i < num_of_particles; i++){
        if (i != j){
            calc_forces(particles[i], particles[j]);

        };
        
    };

};

void apply_boundary(Particle &p){
    if (p.x > size){ p.x -= size; }
    if (p.x < 0) {p.x += size;}
    if (p.y > size){ p.y -= size; }
    if (p.y < 0) {p.y += size;}
    if (p.z > size){ p.z -= size; }
    if (p.z < 0) {p.z += size;}
};

void update_position(std::vector<Particle> &particles){
    for (int16_t i = 0; i < num_of_particles; i++){
        particles[i].x = particles[i].x + particles[i].vx * dt 
        + 1 / 2 * particles[i].ax_prev * dt * dt;
        particles[i].y = particles[i].y + particles[i].vy * dt 
        + 1 / 2 * particles[i].ay_prev * dt * dt;
        particles[i].z = particles[i].z + particles[i].vz * dt 
        + 1 / 2 * particles[i].az_prev * dt * dt;

        //std::cout << "Before: " << particles[i].ax_prev << ' ' << particles[i].ax << "\n";
        particles[i].ax_prev = particles[i].ax;
        particles[i].ay_prev = particles[i].ay;
        particles[i].az_prev = particles[i].az;

        apply_boundary(particles[i]);

    };

    for (int16_t i = 0; i < num_of_particles; i++){
        update_force(particles, i);
        //std::cout << "After: " << 1.0 / 2.0 * (particles[i].ax_prev + particles[i].ax) << "\n";
        particles[i].vx += 1.0 / 2.0 * (particles[i].ax_prev + particles[i].ax) * dt;
        particles[i].vy += 1.0 / 2.0 * (particles[i].ax_prev + particles[i].ax) * dt;
        particles[i].vz += 1.0 / 2.0 * (particles[i].az_prev + particles[i].az) * dt;
    };
    
};

void update(std::vector<Particle> &particles, int16_t n){
    for (int16_t i = 0; i < n; i++){
        t += dt;
        calc_energy(particles, kinetic_energy, potential_energy);
        update_position(particles);

    };
};


int main() {
    


    std::vector<Particle> particles;
    for (double i = 0.4; i < num_of_particles; i++) {
        Particle p;
        p.x = i / 2; p.y = i / 2; p.z = i / 2;
        p.vx = 0.0; p.vy = 0.0; p.vz = 0.0;
        p.ax = 0.0; p.ay = 0.0; p.az = 0.0;
        p.ax_prev = 0.0; p.ay_prev = 0.0; p.az_prev = 0.0;
        particles.push_back(p);
    };
    std::cout<<"\n";
    

    for (int16_t i = 0; i < num_of_particles; i++){
        std::cout << "Kinetic: " << kinetic_energy << " ";
        std::cout << "Potential: " << potential_energy << "\n";
        show(particles[i]);
        update(particles, 1);
    };

    

    
    return 0;
}