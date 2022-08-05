#include<iostream>

struct Particle
{
    double mass = 0.;
    double coord[2] = {0., 0.};
    int particle_id = -1;
    double force[2] = {0., 0.};  // force on each particle
    double velocity[2] = {0., 0.};  // velocity of each particle
};
struct Node 
{ 
    // Store the bounding box of the node
    double top_left_coord[2] = {0, 1e100};
    //double top_right_coord[2] = {0., 0.};
    //double bottom_left_coord[2] = {0., 0.};
    double bottom_right_coord[2] = {1e100, 0};
    double centre_of_mass[2] = {0., 0.};

    // Pointers to children
    bool has_children = false;
    struct Node* child1;  // Top left
    struct Node* child2;  // Top right
    struct Node* child3;  // Bottom left
    struct Node* child4;  // Bottom right

    // Total mass and number of particles present in the bounding box
    double total_mass = 0.;
    int n_particles = 0;

    /* Particle associated with a node. Has values only if 
    no other particle present ithin bounding box*/
    bool has_particle = false;
    Particle particle;
    double total_particle_mass = 0.;
    int particle_id = -1;
    /*
    Node()
    {
        this->n_particles = 0;
        this->total_mass = 0.;
        this->total_particle_mass = 0.;
        this->particle_id = -1;
        this->has_particle = false;
        this->has_children = false;
    }*/
};

class QuadTree
{
    public:
        Node* root = new Node;
        int calc_centre(Node* node);
        bool check_coord_in_box(Node* node, double coord[]);
        void divide_box(Node* node);
        //int insert_particle(Particle* particle);
        int insert_particle(double mass_arr, double coord_x_arr, double coord_y_arr, int particle_id_arr);
        void print_tree();
        //void compute_force(Particle* Particle);
        QuadTree()
        {/*
            this->root->top_left_coord[0] = -1e5;
            this->root->top_left_coord[1] = 1e5;
            this->root->bottom_right_coord[0] = 1e5;
            this->root->bottom_right_coord[1] = -1e5;*/
        }
};

bool QuadTree::check_coord_in_box(Node* node, double coord[])
{
    if ((coord[0] >= node->top_left_coord[0]) && (coord[0] <= node->bottom_right_coord[0]) && (coord[1] >= node->bottom_right_coord[1]) && (coord[1] <= node->top_left_coord[1]))
    {
        return true;
    }
    return false;
}

int QuadTree::calc_centre(Node* node)
{
    node->centre_of_mass[0] = (node->top_left_coord[0] + node->bottom_right_coord[0]) / 2;
    node->centre_of_mass[1] = (node->top_left_coord[1] + node->bottom_right_coord[1]) / 2;

    return 0;
}

void QuadTree::divide_box(Node* node)
{
    node->child1 = new Node;
    node->child2 = new Node;
    node->child3 = new Node;
    node->child4 = new Node;
    node->has_children = true;
    // Top left quadrant
    node->child1->top_left_coord[0] = node->top_left_coord[0];
    node->child1->top_left_coord[1] = node->top_left_coord[1];
    node->child1->bottom_right_coord[0] = (node->top_left_coord[0] + node->bottom_right_coord[0]) / 2;
    node->child1->bottom_right_coord[1] = (node->top_left_coord[1] + node->bottom_right_coord[1]) / 2;

    // Top right quadrant
    node->child2->top_left_coord[0] = node->child1->bottom_right_coord[0];
    node->child2->top_left_coord[1] = node->top_left_coord[1];
    node->child2->bottom_right_coord[0] = node->bottom_right_coord[0];
    node->child2->bottom_right_coord[1] = node->child1->bottom_right_coord[1];
    
    // Bottom left quadrant
    node->child3->top_left_coord[0] = node->top_left_coord[0];
    node->child3->top_left_coord[1] = node->child1->bottom_right_coord[1];
    node->child3->bottom_right_coord[0] = node->child1->bottom_right_coord[0];
    node->child3->bottom_right_coord[1] = node->bottom_right_coord[1];
    
    // Bottom right quadrant
    node->child4->top_left_coord[0] = node->child1->bottom_right_coord[0];
    node->child4->top_left_coord[1] = node->child1->bottom_right_coord[1];
    node->child4->bottom_right_coord[0] = node->bottom_right_coord[0];
    node->child4->bottom_right_coord[1] = node->bottom_right_coord[1];

    if (node->has_particle == true)
    {
        Node* node_child;
        if (this->check_coord_in_box(node->child1, node->particle.coord) == true)
        {
            node_child = node->child1;
        }
        else if (this->check_coord_in_box(node->child2, node->particle.coord) == true)
        {
            node_child = node->child2;
        }
        else if (this->check_coord_in_box(node->child3, node->particle.coord) == true)
        {
            node_child = node->child3;
        }
        else if  (this->check_coord_in_box(node->child4, node->particle.coord) == true)
        {
            node_child = node->child4;
        }
        node_child->particle = node->particle;
        node_child->n_particles = 1;
        node_child->total_particle_mass = node_child->particle.mass;
        node_child->has_particle = true;

        //node->centre_of_mass[0] = node->centre_of_mass[0] * node->total_particle_mass + node_child->particle->coord[0] * node_child->particle->mass;
        //node->centre_of_mass[1] = node->centre_of_mass[1] * node->total_particle_mass + node_child->particle->coord[1] * node_child->particle->mass;
        node_child->centre_of_mass[0] = node->particle.coord[0];
        node_child->centre_of_mass[1] = node->particle.coord[1];

        // node->particle = NULL;
        node->particle = Particle();
        node->has_particle = false;
    }
    
}

//int QuadTree::insert_particle(Particle* particle)
int QuadTree::insert_particle(double mass, double coord_x, double coord_y, int particle_id)
{
    Particle particle = Particle{mass, {coord_x, coord_y}, particle_id};
    Node* current_node = this->root;
    while (1)
    {
        // If coordinate of particle in bounding box of current node
        if (this->check_coord_in_box(current_node, particle.coord) == true)
        {
            // If current node has no children
            if (current_node->has_children == false)
            {
                /* If current node has no assigned particle, assign one.
                Else, create child nodes and move particles to child.*/
                if (current_node->has_particle == false)
                {
                    current_node->particle = particle;
                    current_node->n_particles = 1;
                    current_node->total_particle_mass = mass;
                    current_node->has_particle = true;
                    current_node->centre_of_mass[0] = coord_x;
                    current_node->centre_of_mass[1] = coord_y;
                    break;
                }
                else
                {
                    this->divide_box(current_node);  // Divided into 4 parts and move previously assigned particle to correct leaf node
                }
                
            }
            else
            {
                current_node->centre_of_mass[0] = current_node->centre_of_mass[0] * current_node->total_particle_mass + coord_x * mass;
                current_node->centre_of_mass[1] = current_node->centre_of_mass[1] * current_node->total_particle_mass + coord_y * mass;

                current_node->n_particles += 1;
                current_node->total_particle_mass += mass;

                current_node->centre_of_mass[0] /= current_node->total_particle_mass;
                current_node->centre_of_mass[1] /= current_node->total_particle_mass;
                
                double coord[] = {coord_x, coord_y};
                // Check which box current particle should go to
                if (this->check_coord_in_box(current_node->child1, coord) == true)
                {
                    current_node = current_node->child1;
                }
                else if (this->check_coord_in_box(current_node->child2, coord) == true)
                {
                    current_node = current_node->child2;
                }
                else if (this->check_coord_in_box(current_node->child3, coord) == true)
                {
                    current_node = current_node->child3;
                }
                else if (this->check_coord_in_box(current_node->child4, coord) == true)
                {
                    current_node = current_node->child4;
                }
            }
        }
        else
        {
            std::cout << "Warning: particle outside bound box, mass = " << particle.mass << ", coord={" << coord_x << "," << coord_y << "}" << std::endl;
            break;
        }
        
    }
    return 0;
    
}

void QuadTree::print_tree()
{

}

//int main(int argc, char *argv[])
QuadTree demo()
{
    QuadTree tree;
    Particle p_arr[] = {Particle{5e23, {3e4, 1e3}, 0},
                            Particle{4e21, {7e4, 7e4}, 1},
                            Particle{3e23, {8e4, 9e4}, 2},
                            Particle{7e19, {6.8e4, 4e4}, 3},
                            Particle{8e22, {6.35e4, 3e4}, 4},
                            Particle{9e21, {9e4, 3e4}, 5}};
    
    for (int i = 0; i < 6; i++)
    {
        //tree.insert_particle(&p_arr[i]);
    }
    /*
    Particle p1 = {5e23, {3e4, 1e3}, 0};
    Particle p2 = {4e21, {7e4, 7e4}, 1};
    Particle p3 = {3e23, {8e4, 9e4}, 2};
    Particle p4 = {7e19, {6.8e4, 4e4}, 3};
    Particle p5 = {8e22, {6.35e4, 3e4}, 4};
    Particle p6 = {9e21, {9e4, 3e4}, 5};
    tree.insert_particle(&p1);
    tree.insert_particle(&p2);
    tree.insert_particle(&p3);
    tree.insert_particle(&p4);
    tree.insert_particle(&p5);
    tree.insert_particle(&p6);*/
    
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

    return tree;
}