/*
N-body Simulation Problem
Author: Abhiroop Ghosh
Course: CMSE 822 Parallel Computing Fall 2020
*/
#include <iostream>
#include <cmath>
#include <string>
#include <cstdlib>

const double dt = 1e-2;  // Timestep used for simulation
const double epsi = 1e-12; // Softening factor to prevent division by 0 errors
const double G = 6.674e-11; // Gravitational Constant  m^3 kg^-1 s^-2

// Each particle is represented by a Particle struct defining its position and velocity
struct Particle 
{ 
    double mass = 0.;  // Mass of particle
    double position_x = 0., position_y = 0., position_z = 0.;  // 3D coordinates of the particle
    double force_x = 0., force_y = 0., force_z = 0.;  // velocity of each particle
    double velocity_x = 0., velocity_y = 0., velocity_z = 0.;  // velocity of each particle
};

// Calculate the force and velocity of each particle at a certain point in time
void calculate_forces_and_velocity(Particle particle_array[], int n_particles)
{
    // For every pair of particles (i,j) calculate the forces f(i,j) and velocity v(i,j)
    for (int i = 0; i < n_particles; i++)
    {
        Particle* particle_i = &particle_array[i];
        particle_i->force_x = 0.;
        particle_i->force_y = 0.;
        particle_i->force_z = 0.;
        /*
        particle_i->velocity_x = 0.;
        particle_i->velocity_y = 0.;
        particle_i->velocity_z = 0.;*/
        for (int j = 0; j < n_particles; j++)
        {
            if (i == j)  // Same particle
                continue;
            Particle* particle_j = &particle_array[j];

            // Calculate vector r(i,j) = r(j) - r(i) between particle i and j
            double r_ij_x = particle_j->position_x - particle_i->position_x;
            double r_ij_y = particle_j->position_y - particle_i->position_y;
            double r_ij_z = particle_j->position_z - particle_i->position_z;

            // Calculate absolute distance between particles i and j
            double r_l2_norm = std::sqrt(std::pow(r_ij_x, 2) + std::pow(r_ij_y, 2) + std::pow(r_ij_z, 2));

            // Calculate force f(i,j)
            double force_magnitude = G * (particle_i->mass * particle_j->mass) / (pow(r_l2_norm, 2) + pow(epsi, 2));
            //printf("i=%d, j = %d, force=%f, r=(%f,%f,%f), r_norm=%f, mass_i = %f, mass_j = %f, G = %f\n", i, j, force_magnitude, r_ij_x, r_ij_y, r_ij_z, r_l2_norm, particle_i->mass, particle_j->mass, G);
            particle_i->force_x += force_magnitude * r_ij_x / r_l2_norm;
            particle_i->force_y += force_magnitude * r_ij_y / r_l2_norm;
            particle_i->force_z += force_magnitude * r_ij_z / r_l2_norm;
        }
        
        // Update velocity v(i)
        particle_i->velocity_x += particle_i->force_x * dt / particle_i->mass;
        particle_i->velocity_y += particle_i->force_y * dt / particle_i->mass;
        particle_i->velocity_z += particle_i->force_z * dt / particle_i->mass;
    }
}

/*
Update particle positions based on their velocities. For a time t,
calculate_forces_and_velocity() should be called first to calculate velocities.
*/
void update_positions(Particle particle_array[], int n_particles)
{
    for (int i = 0 ; i < n_particles; i++)
    {
        particle_array[i].position_x += particle_array[i].velocity_x * dt;
        particle_array[i].position_y += particle_array[i].velocity_y * dt;
        particle_array[i].position_z += particle_array[i].velocity_z * dt;
    }
}

void print_particle_system(Particle particle_arr[], int n_particles)
{
    for (int i = 0; i < n_particles; i++)
    {
        printf("Particle %d, Mass = %f, pos = (%f, %f, %f)", i, particle_arr[i].mass, particle_arr[i].position_x, particle_arr[i].position_y, particle_arr[i].position_z);
        printf(", v = (%f, %f, %f)", particle_arr[i].velocity_x, particle_arr[i].velocity_y, particle_arr[i].velocity_z);
        printf(", , F = (%f, %f, %f)\n", particle_arr[i].force_x, particle_arr[i].force_y, particle_arr[i].force_z);
        //std::cout << "Particle " << i << ", Mass = " << particle_arr[i].mass << ", (x, y, z) = (" << particle_arr[i].position_x 
    }
}

int main(int argc, char *argv[])
{
    int max_iter = 1; bool debug = true;
    if (argc == 2)
    {
        max_iter = std::atoi(argv[1]);
    }
    const int N = 6;
    Particle p[N];
    //p[0].mass = 2e24; p[0].position_x = 1e3; p[0].position_y = 3e3; p[0].position_z = 0;
    //p[1].mass = 4e22; p[1].position_x = 4e3; p[1].position_y = 2e3; p[1].position_z = 0;
    //p[2].mass = 3e25; p[2].position_x = 5e3; p[2].position_y = 5e3; p[2].position_z = 0;
    //p[3].mass = 2e18; p[3].position_x = 7e3; p[3].position_y = 4e3; p[3].position_z = 0;
    /*
    p[0].mass = 2; p[0].position_x = 1; p[0].position_y = 3; p[0].position_z = 0;
    p[1].mass = 4; p[1].position_x = 4; p[1].position_y = 2; p[1].position_z = 0;
    p[2].mass = 3; p[2].position_x = 5; p[2].position_y = 5; p[2].position_z = 0;
    p[3].mass = 2; p[3].position_x = 7; p[3].position_y = 4; p[3].position_z = 0;*/

    p[0].mass = 5e23; p[0].position_x = 3e4; p[0].position_y = 1e3; p[0].position_z = 0;
    p[1].mass = 4e21; p[1].position_x = 7e4; p[1].position_y = 7e4; p[1].position_z = 0;
    p[2].mass = 3e23; p[2].position_x = 8e4; p[2].position_y = 9e4; p[2].position_z = 0;
    p[3].mass = 7e19; p[3].position_x = 6.8e4; p[3].position_y = 4e4; p[3].position_z = 0;
    p[4].mass = 8e22; p[4].position_x = 6.35e4; p[4].position_y = 3e4; p[4].position_z = 0;
    p[5].mass = 9e21; p[5].position_x = 9e4; p[5].position_y = 3e4; p[5].position_z = 0;

    std::cout << "\nParticle initial state (t = 0)" << std::endl;
    print_particle_system(p, N);

    for (int i = 0; i < max_iter; i++)
    {
        calculate_forces_and_velocity(p, N);
        update_positions(p, N);
        if (debug && i % 10 == 0)
        {
            std::cout << "\nParticle state, i = " << i << std::endl;
            print_particle_system(p, N);
        }
    }
    
    std::cout << "\nParticle final state" << std::endl;
    print_particle_system(p, N);
}