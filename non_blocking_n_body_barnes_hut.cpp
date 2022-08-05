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
#include <string>
/*
#include "H5Cpp.h"
using namespace H5;
using namespace std; */

//#include "./hdfql-2.2.0/include/HDFql.hpp"
#include "HDFql.hpp"

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

//#pragma omp declare simd
void calculate_velocity(double mass, double* velocity_x, double* velocity_y, double force[])
{
    // Update velocity
    *velocity_x += force[0] * dt / mass;
    *velocity_y += force[1] * dt / mass;
}

//#pragma omp declare simd
void update_positions(double velocity_x, double velocity_y, double* coord_x, double* coord_y)
{
    // Update positions
    *coord_x += velocity_x * dt;
    *coord_y += velocity_y * dt;
}
/*
void create_particles(int n_particles, double mass_arr[], double coord_x_arr[], double coord_y_arr[],
                      double velocity_x_arr[], double velocity_y_arr[], double force_x_arr[], double force_y_arr[],
                      int particle_id_arr[], int n_computations_arr[])
{
    srand(567);  // Seed random number generator
    double mass_ll = 1e20, mass_ul = 1e30, coord_x_ll = 0, coord_x_ul = 1e100, coord_y_ll = 0, coord_y_ul = 1e100;
    for (int i = 0; i < n_particles; i++)
    {
        mass_arr[i] = mass_ll + rand() / (double)RAND_MAX * (mass_ul - mass_ll);
        coord_x_arr[i] = coord_x_ll + rand() / (double)RAND_MAX * (coord_x_ul - coord_x_ll);
        coord_y_arr[i] = coord_y_ll + rand() / (double)RAND_MAX * (coord_y_ul - coord_y_ll);
        std::cout << rand() << "," << RAND_MAX << std::endl;
        particle_id_arr[i] = i;

        velocity_x_arr[i] = 0.;
        velocity_y_arr[i] = 0.;
        force_x_arr[i] = 0.;
        force_y_arr[i] = 0.;
        n_computations_arr[i] = i;
    }
}*/

int main(int argc, char *argv[])
{
    int max_iter = 1000; bool debug = false;
    if (argc >= 2)
    {
        max_iter = std::atoi(argv[1]);
    }
    int n_particles = 1000;
    if (argc >= 3) n_particles = std::atoi(argv[2]);

    QuadTree tree;
    
    /*
    int n_particles = 6;
    double mass_arr[] = {5e23, 4e21, 3e23, 7e19, 8e22, 9e21};
    double coord_x_arr[] = {3e4, 7e6, 8e9, 6.8e1, 6.35e8, 9e9};
    
    double coord_y_arr[] = {1e3, 7e4, 9e3, 4e9, 3e8, 3e9};
    double velocity_x_arr[] = {0, 0, 0, 0, 0, 0};
    double velocity_y_arr[] = {0, 0, 0, 0, 0, 0};
    double force_x_arr[] = {0, 0, 0, 0, 0, 0};
    double force_y_arr[] = {0, 0, 0, 0, 0, 0};
    int particle_id_arr[] = {0, 1, 2, 3, 4, 5};
    int n_computations_arr[] = {0, 0, 0, 0, 0, 0};*/
    
    double mass_arr[n_particles], coord_x_arr[n_particles], coord_y_arr[n_particles], velocity_x_arr[n_particles], velocity_y_arr[n_particles], force_x_arr[n_particles], force_y_arr[n_particles];
    int particle_id_arr[n_particles], n_computations_arr[n_particles];

    //create_particles(n_particles, mass_arr, coord_x_arr, coord_y_arr, velocity_x_arr, velocity_y_arr, force_x_arr, force_y_arr, particle_id_arr, n_computations_arr);
    double mass_ll = 1e10, mass_ul = 1e30, coord_x_ll = 4e10, coord_x_ul = 8e10, coord_y_ll = 4e10, coord_y_ul = 8e10;
    //double mass_ll = 1e20, mass_ul = 1e30, coord_x_ll = 5e10, coord_x_ul = 5e11, coord_y_ll = 5e10, coord_y_ul = 5e11;
    srand(567);  // Seed random number generator
    #pragma omp parallel for
    for (int i = 0; i < n_particles; i++)
    {
        mass_arr[i] = mass_ll + rand() / (double)RAND_MAX * (mass_ul - mass_ll);
        coord_x_arr[i] = coord_x_ll + rand() / (double)RAND_MAX * (coord_x_ul - coord_x_ll);
        coord_y_arr[i] = coord_y_ll + rand() / (double)RAND_MAX * (coord_y_ul - coord_y_ll);
        //std::cout << rand() << "," << RAND_MAX << std::endl;
        particle_id_arr[i] = i;

        velocity_x_arr[i] = 0.;
        velocity_y_arr[i] = 0.;
        force_x_arr[i] = 0.;
        force_y_arr[i] = 0.;
        n_computations_arr[i] = i;
    }

    //printf("Mass 0 = %f", velocity_x_arr[0]);
    for (int i = 0; i < n_particles; i++)
    {
        tree.insert_particle(mass_arr[i], coord_x_arr[i], coord_y_arr[i], particle_id_arr[i]);
    }

    int n_global = n_particles;

    auto start = high_resolution_clock::now(); 
    // Initialize MPI and get current rank and no. of processors
    //MPI_Init(&argc, &argv);
    int rank, n_processors;
    MPI_Comm_size(MPI_COMM_WORLD, &n_processors);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Evenly distribute total number of particles among processors
    int n_local = (int) n_global / n_processors;  
    int extra_particles = n_global % n_processors;
    /* For good load balancing, each processor takes 1 extra particle
    till all such particles are assigned.*/
    if (rank < extra_particles) n_local++;

    MPI_Request m_req_send, cx_req_send, cy_req_send, vx_req_send, vy_req_send, fx_req_send, fy_req_send, pid_req_send, n_comp_req_send; // Send from rank 0
    MPI_Request m_req_recv, cx_req_recv, cy_req_recv, vx_req_recv, vy_req_recv, fx_req_recv, fy_req_recv, pid_req_recv, n_comp_req_recv; // Receive from rank 0
    MPI_Request m_req_send_res, cx_req_send_res, cy_req_send_res, vx_req_send_res, vy_req_send_res, fx_req_send_res, fy_req_send_res, pid_req_send_res, n_comp_req_send_res; // Send results to rank 0
    MPI_Request m_req_recv_res, cx_req_recv_res, cy_req_recv_res, vx_req_recv_res, vy_req_recv_res, fx_req_recv_res, fy_req_recv_res, pid_req_recv_res, n_comp_req_recv_res; // Receive results in rank 0
    
    MPI_Status m_status_send, cx_status_send, cy_status_send, vx_status_send, vy_status_send, fx_status_send, fy_status_send, pid_status_send, n_comp_status_send; // Send from rank 0
    MPI_Status m_status_recv, cx_status_recv, cy_status_recv, vx_status_recv, vy_status_recv, fx_status_recv, fy_status_recv, pid_status_recv, n_comp_status_recv; // Receive from rank 0
    MPI_Status m_status_send_res, cx_status_send_res, cy_status_send_res, vx_status_send_res, vy_status_send_res, fx_status_send_res, fy_status_send_res, pid_status_send_res, n_comp_status_send_res; // Send results to rank 0
    MPI_Status m_status_recv_res, cx_status_recv_res, cy_status_recv_res, vx_status_recv_res, vy_status_recv_res, fx_status_recv_res, fy_status_recv_res, pid_status_recv_res, n_comp_status_recv_res; // Receive results in rank 0
    
    // Perform simulation over multiple timsteps
    for (int t = 0; t < max_iter; t++)
    {
        std::cout << "Rank = " << rank << ", Step = " << t+1 << ", Time = " << (t + 1)*dt << std::endl;
/*
        QuadTree tree;
        for (int i = 0; i < n_particles; i++)
        {
            tree.insert_particle(mass_arr[i], coord_x_arr[i], coord_y_arr[i], particle_id_arr[i]);
        }*/
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
            //#pragma omp parallel for
            for (int r = 1; r < n_processors; r++)
            { 
                int elements_per_process = (int) n_global / n_processors;
                int index = (elements_per_process + 1) * std::min(r, extra_particles) + elements_per_process * std::max(0, r - extra_particles);
                if (r < extra_particles)
                {
                    elements_per_process++;
                }
                //MPI_Send(&elements_per_process, 1, MPI_INT, r, 0, MPI_COMM_WORLD); 
                MPI_Isend(&mass_arr[index], elements_per_process, MPI_DOUBLE, r, 0, MPI_COMM_WORLD, &m_req_send); 
                MPI_Isend(&coord_x_arr[index], elements_per_process, MPI_DOUBLE, r, 0, MPI_COMM_WORLD, &cx_req_send); 
                MPI_Isend(&coord_y_arr[index], elements_per_process, MPI_DOUBLE, r, 0, MPI_COMM_WORLD, &cy_req_send); 
                MPI_Isend(&velocity_x_arr[index], elements_per_process, MPI_DOUBLE, r, 0, MPI_COMM_WORLD, &vx_req_send); 
                MPI_Isend(&velocity_y_arr[index], elements_per_process, MPI_DOUBLE, r, 0, MPI_COMM_WORLD, &vy_req_send); 
                MPI_Isend(&particle_id_arr[index], elements_per_process, MPI_INT, r, 0, MPI_COMM_WORLD, &pid_req_send); 
                MPI_Isend(&n_computations_arr[index], elements_per_process, MPI_INT, r, 0, MPI_COMM_WORLD, &n_comp_req_send); 
            }
            #pragma omp parallel for
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
            MPI_Irecv(&m_local[0], n_local, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &m_req_recv); 
            MPI_Irecv(&c_x_local[0], n_local, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &cx_req_recv); 
            MPI_Irecv(&c_y_local[0], n_local, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &cy_req_recv); 
            MPI_Irecv(&v_x_local[0], n_local, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &vx_req_recv); 
            MPI_Irecv(&v_y_local[0], n_local, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &vy_req_recv); 
            MPI_Irecv(&pid_local[0], n_local, MPI_INT, 0, 0, MPI_COMM_WORLD, &pid_req_recv); 
            MPI_Irecv(&n_computations_arr[0], n_local, MPI_INT, 0, 0, MPI_COMM_WORLD, &n_comp_req_recv); 
            
            MPI_Request req_recv_arr[] = {m_req_recv, cx_req_recv, cy_req_recv, vx_req_recv, vy_req_recv, pid_req_recv, n_comp_req_recv};
            MPI_Status status_recv_arr[] = {m_status_recv, cx_status_recv, cy_status_recv, vx_status_recv, vy_status_recv, pid_status_recv, n_comp_status_recv};
            MPI_Waitall(7, req_recv_arr, status_recv_arr);
        }

        // Perform calculation
        #pragma omp parallel for
        for (int i = 0; i < n_local; i++)
        {
            int n_computations = 0;
            // Force calculation
            double force[] = {0., 0.};
            compute_force(tree.root, m_local[i], c_x_local[i], c_y_local[i], pid_local[i], force, &n_computations);
            calculate_velocity(m_local[i], &v_x_local[i], &v_y_local[i], force);
            f_x_local[i] = force[0];
            f_y_local[i] = force[1];
            n_comp_local[i] = n_computations;
            //std::cout << "Rank = " << rank << "," << m_local[i] << "," << f_x_local[i] << "," << f_y_local[i] << std::endl;
        }
        #pragma omp parallel for
        for (int i = 0; i < n_local; i++)
        {
            update_positions(v_x_local[i], v_y_local[i], &c_x_local[i],&c_y_local[i]);
        }
        // Gather results from all processors to rank 0
        if (rank == 0)
        {
            for (int r = 1; r < n_processors; r++)
            { 
                int elements_per_process = (int) n_global / n_processors;
                int index = (elements_per_process + 1) * std::min(r, extra_particles) + elements_per_process * std::max(0, r - extra_particles);
                if (r < extra_particles)
                {
                    elements_per_process++;
                }
                MPI_Irecv(&mass_arr[index], elements_per_process, MPI_DOUBLE, r, 0, MPI_COMM_WORLD, &m_req_recv_res); 
                MPI_Irecv(&coord_x_arr[index], elements_per_process, MPI_DOUBLE, r, 0, MPI_COMM_WORLD, &cx_req_recv_res); 
                MPI_Irecv(&coord_y_arr[index], elements_per_process, MPI_DOUBLE, r, 0, MPI_COMM_WORLD, &cy_req_recv_res); 
                MPI_Irecv(&velocity_x_arr[index], elements_per_process, MPI_DOUBLE, r, 0, MPI_COMM_WORLD, &vx_req_recv_res); 
                MPI_Irecv(&velocity_y_arr[index], elements_per_process, MPI_DOUBLE, r, 0, MPI_COMM_WORLD, &vy_req_recv_res); 
                MPI_Irecv(&force_x_arr[index], elements_per_process, MPI_DOUBLE, r, 0, MPI_COMM_WORLD, &fx_req_recv_res); 
                MPI_Irecv(&force_y_arr[index], elements_per_process, MPI_DOUBLE, r, 0, MPI_COMM_WORLD, &fy_req_recv_res); 
                MPI_Irecv(&particle_id_arr[index], elements_per_process, MPI_INT, r, 0, MPI_COMM_WORLD, &pid_req_recv_res);
                MPI_Irecv(&n_computations_arr[index], elements_per_process, MPI_INT, r, 0, MPI_COMM_WORLD, &n_comp_req_recv_res);
                //printf("index = %d, Mass recv = %f\n", index, mass_arr[index]);

                MPI_Request req_recv_res_arr[] = {m_req_recv_res, cx_req_recv_res, cy_req_recv_res, vx_req_recv_res, vy_req_recv_res, fx_req_recv_res, fy_req_recv_res, pid_req_recv_res, n_comp_req_recv_res};
                MPI_Status status_recv_res_arr[] = {m_status_recv_res, cx_status_recv_res, cy_status_recv_res, vx_status_recv_res, vy_status_recv_res, fx_status_recv_res, fy_status_recv_res, pid_status_recv_res, n_comp_status_recv_res};
                
                // Perform calculation while waiting on receive
                if (r == 1)
                {
                    #pragma omp parallel for
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
                MPI_Waitall(9, req_recv_res_arr, status_recv_res_arr);
            }
            /*
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
                    }*/
            std::cout << " t = " << t << ", rank = " << rank << std::endl;

            // Write results to HDF5 file
            if (t == 0)
            {
                // create an HDF5 file
                HDFql::execute("CREATE TRUNCATE FILE results.h5");
            }          
            
            if ((t + 1) % 100 == 0)
            {
                // Open HDF5 file
                HDFql::execute("USEFILE results.h5");
                
                // Register the result arrays to be recorded
                int number_mass = HDFql::variableRegister(mass_arr);
                int number_coord_x = HDFql::variableRegister(coord_x_arr);
                int number_coord_y = HDFql::variableRegister(coord_y_arr);
                int number_v_x = HDFql::variableRegister(velocity_x_arr);
                int number_v_y = HDFql::variableRegister(velocity_y_arr);
                int number_f_x = HDFql::variableRegister(force_x_arr);
                int number_f_y = HDFql::variableRegister(force_y_arr);
                int number_pid = HDFql::variableRegister(particle_id_arr);
                int number_ncomp = HDFql::variableRegister(n_computations_arr);

                // Create group for current timestep and corresponding dataset
                std::string dataset_name = "Step" + std::to_string(t + 1) + "/mass";
                std::string dataset_create_query = "CREATE DATASET " + dataset_name + " AS FLOAT(" + std::to_string(n_global) + ")";
                HDFql::execute(dataset_create_query.c_str());
                std::string dataset_query = "SELECT FROM " + dataset_name;
                HDFql::execute(dataset_query.c_str());
                std::string ins = "INSERT INTO " + dataset_name + " VALUES FROM MEMORY " + std::to_string(number_mass); HDFql::execute(ins.c_str());
                /*
                std::cout << dataset_create_query << std::endl;
                std::cout << dataset_query << std::endl;
                std::cout << ins << std::endl;*/
                
                dataset_name = "Step" + std::to_string(t + 1) + "/cx";
                dataset_create_query = "CREATE DATASET " + dataset_name + " AS FLOAT(" + std::to_string(n_global) + ")";
                HDFql::execute(dataset_create_query.c_str());
                dataset_query = "SELECT FROM " + dataset_name;
                HDFql::execute(dataset_query.c_str());
                ins = "INSERT INTO " + dataset_name + " VALUES FROM MEMORY " + std::to_string(number_coord_x); HDFql::execute(ins.c_str());
                
                dataset_name = "Step" + std::to_string(t + 1) + "/cy";
                dataset_create_query = "CREATE DATASET " + dataset_name + " AS FLOAT(" + std::to_string(n_global) + ")";
                HDFql::execute(dataset_create_query.c_str());
                dataset_query = "SELECT FROM " + dataset_name;
                HDFql::execute(dataset_query.c_str());
                ins = "INSERT INTO " + dataset_name + " VALUES FROM MEMORY " + std::to_string(number_coord_y); HDFql::execute(ins.c_str());
                
                dataset_name = "Step" + std::to_string(t + 1) + "/vx";
                dataset_create_query = "CREATE DATASET " + dataset_name + " AS FLOAT(" + std::to_string(n_global) + ")";
                HDFql::execute(dataset_create_query.c_str());
                dataset_query = "SELECT FROM " + dataset_name;
                HDFql::execute(dataset_query.c_str());
                ins = "INSERT INTO " + dataset_name + " VALUES FROM MEMORY " + std::to_string(number_v_x); HDFql::execute(ins.c_str());
                
                dataset_name = "Step" + std::to_string(t + 1) + "/vy";
                dataset_create_query = "CREATE DATASET " + dataset_name + " AS FLOAT(" + std::to_string(n_global) + ")";
                HDFql::execute(dataset_create_query.c_str());
                dataset_query = "SELECT FROM " + dataset_name;
                HDFql::execute(dataset_query.c_str());
                ins = "INSERT INTO " + dataset_name + " VALUES FROM MEMORY " + std::to_string(number_v_y); HDFql::execute(ins.c_str());
                
                dataset_name = "Step" + std::to_string(t + 1) + "/fx";
                dataset_create_query = "CREATE DATASET " + dataset_name + " AS FLOAT(" + std::to_string(n_global) + ")";
                HDFql::execute(dataset_create_query.c_str());
                dataset_query = "SELECT FROM " + dataset_name;
                HDFql::execute(dataset_query.c_str());
                ins = "INSERT INTO " + dataset_name + " VALUES FROM MEMORY " + std::to_string(number_f_x); HDFql::execute(ins.c_str());
                
                dataset_name = "Step" + std::to_string(t + 1) + "/fy";
                dataset_create_query = "CREATE DATASET " + dataset_name + " AS FLOAT(" + std::to_string(n_global) + ")";
                HDFql::execute(dataset_create_query.c_str());
                dataset_query = "SELECT FROM " + dataset_name;
                HDFql::execute(dataset_query.c_str());
                ins = "INSERT INTO " + dataset_name + " VALUES FROM MEMORY " + std::to_string(number_f_y); HDFql::execute(ins.c_str());
                
                dataset_name = "Step" + std::to_string(t + 1) + "/pid";
                dataset_create_query = "CREATE DATASET " + dataset_name + " AS FLOAT(" + std::to_string(n_global) + ")";
                HDFql::execute(dataset_create_query.c_str());
                dataset_query = "SELECT FROM " + dataset_name;
                HDFql::execute(dataset_query.c_str());
                ins = "INSERT INTO " + dataset_name + " VALUES FROM MEMORY " + std::to_string(number_pid); HDFql::execute(ins.c_str());
                
                dataset_name = "Step" + std::to_string(t + 1) + "/ncomp";
                dataset_create_query = "CREATE DATASET " + dataset_name + " AS FLOAT(" + std::to_string(n_global) + ")";
                HDFql::execute(dataset_create_query.c_str());
                dataset_query = "SELECT FROM " + dataset_name;
                HDFql::execute(dataset_query.c_str());
                ins = "INSERT INTO " + dataset_name + " VALUES FROM MEMORY " + std::to_string(number_ncomp); HDFql::execute(ins.c_str());

                HDFql::variableUnregister(mass_arr);
            }
        }
        // Send results to rank 0
        else
        {
            //printf("Mass send = %f\n", m_local[0]);
            MPI_Isend(&m_local[0], n_local, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &m_req_send_res); 
            MPI_Isend(&c_x_local[0], n_local, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &cx_req_send_res); 
            MPI_Isend(&c_y_local[0], n_local, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &cy_req_send_res); 
            MPI_Isend(&v_x_local[0], n_local, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &vx_req_send_res); 
            MPI_Isend(&v_y_local[0], n_local, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &vy_req_send_res); 
            MPI_Isend(&f_x_local[0], n_local, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &fx_req_send_res); 
            MPI_Isend(&f_y_local[0], n_local, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &fy_req_send_res); 
            MPI_Isend(&pid_local[0], n_local, MPI_INT, 0, 0, MPI_COMM_WORLD, &pid_req_send_res);
            MPI_Isend(&n_comp_local[0], n_local, MPI_INT, 0, 0, MPI_COMM_WORLD, &n_comp_req_send_res);
        }
    }
    // call Finalize
    //MPI_Finalize();
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
    printf("Using %d particles\n", n_particles);
    printf("Using %d timesteps\n", max_iter);
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