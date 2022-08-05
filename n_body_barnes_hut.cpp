/*
N-body Simulation Problem
Author: Abhiroop Ghosh
Course: CMSE 822 Parallel Computing Fall 2020
*/
#include <iostream>
#include <cmath>
#include <string>
#include <cstdlib>
#include <mpi.h>
#include "quadtree.cpp"
#include <chrono> 

using namespace std::chrono; 

const double dt = 1e-2;  // Timestep used for simulation
const double epsi = 1e-12; // Softening factor to prevent division by 0 errors
const double G = 6.674e-11; // Gravitational Constant  m^3 kg^-1 s^-2
const double thetaTol = 0.5; // Opening angle tolerance
bool debug = false;

//void compute_force(Node* node, Particle* particle, double force[2], int* n_computations)
void compute_force(Node* node, double mass, double coord_x, double coord_y, int particle_id, double force[2], int* n_computations)
{
    double node_len = std::abs(node->bottom_right_coord[0] - node->top_left_coord[0]);
    double dist = std::sqrt(std::pow(coord_x - node->centre_of_mass[0], 2) + std::pow(coord_y - node->centre_of_mass[1], 2));
    double theta = node_len / dist;
    //force[0] = 0.; force[1] = 0.;

    //SIMD??
    // check if node is a leaf node. If node not a leaf node, use centre of mass to calc force if theta < tol
    if (((node->has_children == false) && (node->has_particle == true)) || ((node->has_children == true) && (theta < thetaTol)))
    {
        double r[] = {node->centre_of_mass[0] - coord_x, node->centre_of_mass[1] - coord_y};
        // If distance between particles is zero
        if ((std::abs(r[0]) <= 1e-5) && (std::abs(r[1]) <= 1e-5))
            return;
        double r_squared = std::pow(r[0], 2) + std::pow(r[1], 2);
        double r_l2_norm = std::sqrt(r_squared);
        double force_magnitude = G * (node->total_particle_mass * mass) / (r_squared + std::pow(epsi, 2));
        force[0] += force_magnitude * r[0] / r_l2_norm;
        force[1] += force_magnitude * r[1] / r_l2_norm;
        *n_computations += 1;

        // for debug
        if ((debug == true) && (node->has_children == true) && (theta < thetaTol))
        {
            std::cout << node->total_particle_mass << std::endl;
        }
        return;
    }
    // SIMD??
    else if (node->has_children == true)
    {
        compute_force(node->child1, mass, coord_x, coord_y, particle_id, force, n_computations);
        compute_force(node->child2, mass, coord_x, coord_y, particle_id, force, n_computations);
        compute_force(node->child3, mass, coord_x, coord_y, particle_id, force, n_computations);
        compute_force(node->child4, mass, coord_x, coord_y, particle_id, force, n_computations);

    }
    else
    {
        return;
    }
    
}

#pragma omp declare simd
void calculate_velocity(double mass, double* velocity_x, double* velocity_y, double force[])
{
    // Update velocity
    *velocity_x += force[0] * dt / mass;
    *velocity_y += force[1] * dt / mass;
}

#pragma omp declare simd
void update_positions(double velocity_x, double velocity_y, double* coord_x, double* coord_y)
{
    // Update positions
    *coord_x += velocity_x * dt;
    *coord_y += velocity_y * dt;
}

int main(int argc, char *argv[])
{
    int max_iter = 1000; bool debug = false;
    if (argc == 2)
    {
        max_iter = std::atoi(argv[1]);
    }
    QuadTree tree;
    
    int n_particles = 6;
    double mass_arr[] = {5e23, 4e21, 3e23, 7e19, 8e22, 9e21};
    double coord_x_arr[] = {3e4, 7e6, 8e9, 6.8e1, 6.35e8, 9e9};
    
    double coord_y_arr[] = {1e3, 7e4, 9e3, 4e9, 3e8, 3e9};
    double velocity_x_arr[] = {0, 0, 0, 0, 0, 0};
    double velocity_y_arr[] = {0, 0, 0, 0, 0, 0};
    double force_x_arr[] = {0, 0, 0, 0, 0, 0};
    double force_y_arr[] = {0, 0, 0, 0, 0, 0};
    int particle_id_arr[] = {0, 1, 2, 3, 4, 5};
    int n_computations_arr[] = {0, 0, 0, 0, 0, 0};
    
    for (int i = 0; i < n_particles; i++)
    {
        tree.insert_particle(mass_arr[i], coord_x_arr[i], coord_y_arr[i], particle_id_arr[i]);
    }

    int n_global = n_particles;

    auto start = high_resolution_clock::now(); 
    // Initialize MPI and get current rank and no. of processors
    MPI_Init(&argc, &argv);
    int rank, n_processors;
    MPI_Comm_size(MPI_COMM_WORLD, &n_processors);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Evenly distribute total number of particles among processors
    int n_local = (int) n_global / n_processors;  
    int extra_particles = n_global % n_processors;
    /* For good load balancing, each processor takes 1 extra particle
    till all such particles are assigned.*/
    if (rank < extra_particles) n_local++;

    MPI_Request tree_request, mass_arr_request, coord_x_arr_request, coord_y_arr_request, velocity_x_arr_request, velocity_y_arr_request;
    MPI_Request tree_request_recv, mass_arr_request_recv, coord_x_arr_request_recv, coord_y_arr_request_recv, velocity_x_arr_request_recv, velocity_y_arr_request_recv;
    // Perform simulation over multiple timsteps
    for (int t = 0; t < max_iter; t++)
    {
        std::cout << "Rank = " << rank << ", Step = " << t+1 << ", Time = " << (t + 1)*dt << std::endl;

        // Create local array for particles
        double m_local[n_local];
        double c_x_local[n_local], c_y_local[n_local];
        double v_x_local[n_local], v_y_local[n_local];
        double f_x_local[n_local], f_y_local[n_local];
        int pid_local[n_local], n_comp_local[n_local];
        
        // Sender process
        if (rank == 0)
        {
            // Divide array into chunks and send to each processor
            for (int r = 1; r < n_processors; r++)
            { 
                int elements_per_process = (int) n_global / n_processors;
                int index = (elements_per_process + 1) * std::min(r, extra_particles) + elements_per_process * std::max(0, r - extra_particles);
                if (r < extra_particles)
                {
                    elements_per_process++;
                }
                //MPI_Send(&elements_per_process, 1, MPI_INT, r, 0, MPI_COMM_WORLD); 
                MPI_Send(&mass_arr[index], elements_per_process, MPI_DOUBLE, r, 0, MPI_COMM_WORLD); 
                MPI_Send(&coord_x_arr[index], elements_per_process, MPI_DOUBLE, r, 0, MPI_COMM_WORLD); 
                MPI_Send(&coord_y_arr[index], elements_per_process, MPI_DOUBLE, r, 0, MPI_COMM_WORLD); 
                MPI_Send(&velocity_x_arr[index], elements_per_process, MPI_DOUBLE, r, 0, MPI_COMM_WORLD); 
                MPI_Send(&velocity_y_arr[index], elements_per_process, MPI_DOUBLE, r, 0, MPI_COMM_WORLD); 
                MPI_Send(&particle_id_arr[index], elements_per_process, MPI_INT, r, 0, MPI_COMM_WORLD); 
                MPI_Send(&n_computations_arr[index], elements_per_process, MPI_INT, r, 0, MPI_COMM_WORLD); 
            }
            for(int i = 0; i < n_local; i++)
            {
                m_local[i] = mass_arr[i];
                c_x_local[i] = coord_x_arr[i];
                c_y_local[i] = coord_y_arr[i];
                v_x_local[i] = velocity_x_arr[i];
                v_y_local[i] = velocity_y_arr[i];
                f_x_local[i] = force_x_arr[i];
                f_y_local[i] = force_y_arr[i];
                pid_local[i] = particle_id_arr[i];
                n_comp_local[i] = n_computations_arr[i];
            }
        }
        // Receive array chunks sent by rank 0
        else
        {
            int a = 0;
            //MPI_Recv(&a, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
            MPI_Recv(&m_local[0], n_local, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
            MPI_Recv(&c_x_local[0], n_local, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
            MPI_Recv(&c_y_local[0], n_local, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
            MPI_Recv(&v_x_local[0], n_local, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
            MPI_Recv(&v_y_local[0], n_local, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
            MPI_Recv(&pid_local[0], n_local, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
            MPI_Recv(&n_computations_arr[0], n_local, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
        }

        // Perform calculation
        #pragma omp simd
        for (int i = 0; i < n_local; i++)
        {
            int n_computations = 0;
            // Force calculation
            double force[] = {0., 0.};
            compute_force(tree.root, m_local[i], c_x_local[i], c_y_local[i], pid_local[i], force, &n_computations);
            calculate_velocity(m_local[i], &v_x_local[i], &v_y_local[i], force);
            update_positions(v_x_local[i], v_y_local[i], &c_x_local[i],&c_y_local[i]);
            f_x_local[i] = force[0];
            f_y_local[i] = force[1];
            n_comp_local[i] = n_computations;
            //std::cout << "Rank = " << rank << "," << m_local[i] << "," << f_x_local[i] << "," << f_y_local[i] << std::endl;
        }
        
        // Gather results from all processors to rank 0
        if (rank == 0)
        {
            #pragma omp simd
            for (int r = 1; r < n_processors; r++)
            { 
                int elements_per_process = (int) n_global / n_processors;
                int index = (elements_per_process + 1) * std::min(r, extra_particles) + elements_per_process * std::max(0, r - extra_particles);
                if (r < extra_particles)
                {
                    elements_per_process++;
                }
                MPI_Recv(&mass_arr[index], elements_per_process, MPI_DOUBLE, r, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
                MPI_Recv(&coord_x_arr[index], elements_per_process, MPI_DOUBLE, r, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
                MPI_Recv(&coord_y_arr[index], elements_per_process, MPI_DOUBLE, r, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
                MPI_Recv(&velocity_x_arr[index], elements_per_process, MPI_DOUBLE, r, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
                MPI_Recv(&velocity_y_arr[index], elements_per_process, MPI_DOUBLE, r, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
                MPI_Recv(&force_x_arr[index], elements_per_process, MPI_DOUBLE, r, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
                MPI_Recv(&force_y_arr[index], elements_per_process, MPI_DOUBLE, r, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
                MPI_Recv(&particle_id_arr[index], elements_per_process, MPI_INT, r, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
                MPI_Recv(&n_computations_arr[index], elements_per_process, MPI_INT, r, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                //printf("index = %d, Mass recv = %f\n", index, mass_arr[index]);
            }
        }
        // Send results to rank 0
        else
        {
            //printf("Mass send = %f\n", m_local[0]);
            MPI_Send(&m_local[0], n_local, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD); 
            MPI_Send(&c_x_local[0], n_local, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD); 
            MPI_Send(&c_y_local[0], n_local, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD); 
            MPI_Send(&v_x_local[0], n_local, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD); 
            MPI_Send(&v_y_local[0], n_local, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD); 
            MPI_Send(&f_x_local[0], n_local, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD); 
            MPI_Send(&f_y_local[0], n_local, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD); 
            MPI_Send(&pid_local[0], n_local, MPI_INT, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&n_comp_local[0], n_local, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
        
        if (rank == 0)
        {
            #pragma omp simd
            for (int i = 0; i < n_local; i++)
            {
                mass_arr[i] = m_local[i];
                coord_x_arr[i] = c_x_local[i];
                coord_y_arr[i] = c_y_local[i];
                force_x_arr[i] = f_x_local[i];
                force_y_arr[i] = f_y_local[i];
                velocity_x_arr[i] = v_x_local[i];
                velocity_y_arr[i] = v_y_local[i];
                particle_id_arr[i] = pid_local[i];
                n_computations_arr[i] = n_comp_local[i];
            }
        }
    }

    // call Finalize
    MPI_Finalize();
    auto stop = high_resolution_clock::now(); 
  
    // Get duration
    auto duration = duration_cast<milliseconds>(stop - start); 

    if (rank == 0)
    {
        #pragma omp simd
        for (int i = 0; i < n_global; i++)
        {
            std::cout << "i = " << particle_id_arr[i] << ", Mass = " << mass_arr[i] << ", Force = " << force_x_arr[i] << ", " << force_y_arr[i];
            std::cout << ", velocity = " << velocity_x_arr[i] << "," << velocity_y_arr[i];
            std::cout << ", new pos = " << coord_x_arr[i] << "," << coord_y_arr[i] << ", no. of comps = " << n_computations_arr[i] << std::endl;
        }
        std::cout << "Total time taken = " << duration.count() << std::endl;
    }
    /*
    Node* curr_node = tree.root;
    std::cout << "Root, No. of particles = " << curr_node->n_particles << "," << "Total mass = " << curr_node->total_particle_mass;
    std::cout << ", centre of mass = " << curr_node->centre_of_mass[0] << "," << curr_node->centre_of_mass[1] << std::endl;
    
    curr_node = tree.root->child2;
    std::cout << "Root,child2, No. of particles = " << curr_node->n_particles << "," << "Total mass = " << curr_node->total_particle_mass;
    std::cout << ", centre of mass = " << curr_node->centre_of_mass[0] << "," << curr_node->centre_of_mass[1] << std::endl;
    
    curr_node = tree.root->child2->child2;
    std::cout << "Root,child2,child2, No. of particles = " << curr_node->n_particles << "," << "Total mass = " << curr_node->total_particle_mass;
    std::cout << ", centre of mass = " << curr_node->centre_of_mass[0] << "," << curr_node->centre_of_mass[1] << std::endl;
    
    curr_node = tree.root->child2->child3;
    std::cout << "Root,child2,child3, No. of particles = " << curr_node->n_particles << "," << "Total mass = " << curr_node->total_particle_mass;
    std::cout << ", centre of mass = " << curr_node->centre_of_mass[0] << "," << curr_node->centre_of_mass[1] << std::endl;
    
    curr_node = tree.root->child3;
    std::cout << "Root,child3, No. of particles = " << curr_node->n_particles << "," << "Total mass = " << curr_node->total_particle_mass;
    std::cout << ", centre of mass = " << curr_node->centre_of_mass[0] << "," << curr_node->centre_of_mass[1] << std::endl;
    
    curr_node = tree.root->child4;
    std::cout << "Root,child4, No. of particles = " << curr_node->n_particles << "," << "Total mass = " << curr_node->total_particle_mass;
    std::cout << ", centre of mass = " << curr_node->centre_of_mass[0] << "," << curr_node->centre_of_mass[1] << std::endl;
    
    curr_node = tree.root->child4->child1;
    std::cout << "Root,child4,child1, No. of particles = " << curr_node->n_particles << "," << "Total mass = " << curr_node->total_particle_mass;
    std::cout << ", centre of mass = " << curr_node->centre_of_mass[0] << "," << curr_node->centre_of_mass[1] << std::endl;

    curr_node = tree.root->child4->child1->child2;
    std::cout << "Root,child4,child1,child2, No. of particles = " << curr_node->n_particles << "," << "Total mass = " << curr_node->total_particle_mass;
    std::cout << ", centre of mass = " << curr_node->centre_of_mass[0] << "," << curr_node->centre_of_mass[1] << std::endl;

    curr_node = tree.root->child4->child1->child4;
    std::cout << "Root,child4,child1,child4, No. of particles = " << curr_node->n_particles << "," << "Total mass = " << curr_node->total_particle_mass;
    std::cout << ", centre of mass = " << curr_node->centre_of_mass[0] << "," << curr_node->centre_of_mass[1] << std::endl;

    curr_node = tree.root->child4->child2;
    std::cout << "Root,child4,child2, No. of particles = " << curr_node->n_particles << "," << "Total mass = " << curr_node->total_particle_mass;
    std::cout << ", centre of mass = " << curr_node->centre_of_mass[0] << "," << curr_node->centre_of_mass[1] << std::endl;

    std::cout << std::endl;
    std::cout << std::endl;*/
}