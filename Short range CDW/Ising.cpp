#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <vector>

/*
   Program performs the Monte Carlo simulation of the Ising type hamiltonian:
   H=J_1 \sum_{\langle i j\rangle} s_i s_j+J_2 \sum_{\langle i j\rangle} \sigma_i s_j,
   where the spin variable s_i (B-site) takes the values {+1,-1} and sigma_i (A-site) takes values {1,0}.
   It is alse equiped with some analysis methods and basic I/O.

*/

const int size = 10;              // supercell size
const int N = size * size * size; // number of cells in the supercell
const int sites = 2;              // number of sites -> We have two Ising variables, therefore two sites.

/* Mersenne twister random generators */
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_int_distribution<> distrib(0, 1);
std::uniform_int_distribution<> pick(0, size - 1);
std::uniform_real_distribution<> dis(0.0, 1.0);

class Model
{   bool l[size][size][size][sites]; // 4D boolean array to store Ising variables
    float t;                         // temperature in units of 1/k_B
    float x;                         // A-site composition ratio (nominal)
    float x_real;                    // achieved A-site composition ratio
    float y;                         // B-site composition ratio (determined by charge neutrality)
    float y_real;                    // achieved B-site composition ratio
    float J1 = 0;                    // first Ising coupling
    float J2 = 0;                    // second Ising coupling
    bool glauber_test = 0;           // by default Metropolis test is used; Glauber better at higher temperatures

private:
    // Methods used by MC simulation

    /* Basic MC move: Swap B-site spins s_ijk -> s_klm */
    void swap_spins(int, int, int, int, int, int);
    /* Calculates the sum of first neighbor spin variables implementing periodic boundary conditions.
     Used for energy change calculations.*/
    float neighbor(int, int, int, int);
    /*Interaction energy of a given spin s_ijk*/
    float interaction(int, int, int);
    /*Accept or reject a proposed move.*/
    bool test(float);
    /*One MC step: Propose a move, calculate the energy change and reject/accept on the basis of a Metropolis/Glauber
     * test*/
    int update(int, int, int, int, int, int);

public:
    /*Class constructor: takes in temperature, doping and Ising couplings and
      initializes the Model. */
    Model(float, float, float, float);

    // various getters and setters
    /*Sets the option to use Glauber test instead of Metropolis test during MC sampling. */
    void set_glauber(bool);
    /*Set Ising couplings.*/
    void set_couplings(float, float);
    /*Set temperature.*/
    void set_temp(float);
    /*Set A-site composition (B-site composition is automatically fixed by this)*/
    void set_x(float);
    /*Get A-site composition*/
    float get_x();
    /*Get B-site composition*/
    float get_y();
    /*Convert booleans to spin variables*/
    int get_s(int, int, int);
    int get_sigma(int, int, int);

    // sampling methods

    void sampling_cycle(int);
    /*Calculates the staggered magnetization for a particular microstate*/
    float magnetization();
    /*Measures the order parameter using MC importance sampling.*/
    void order_parameter(int, int, int);
    /*Internal storage of the measured order_parameter*/
    float order = 0;
    // translation methods

    /*Returns the atom name associated with a given spin variable */
    std::string get_atom(int, int, int, int);
    /*Outputs the atom list to a file. This is used to create an initial structure in javelin.*/
    void output_structure(std::string, int);
};

/* Measures the doping dependence of the order parameter. */
void sweep_x(float, float, float, float, float, int);
/* Generates a specified number of thermalized configurations */
void get_sample_configurations(float, float, float, float, int, std::string);

int main()
{   float temperature, J1, J2, from, to, points, doping, samples;
    std::string output_path;
    std::string option = "";

    std::cout << "Enter temperature: ";
    std::cin >> temperature;

    std::cout << "Enter J1: ";
    std::cin >> J1;

    std::cout << "Enter J2: ";
    std::cin >> J2;

    std::cout << "Enter option (sweep_doping/get_samples): ";
    std::cin >> option;

    if (option == "sweep_doping")
    {

        std::cout << "Sweep doping from: ";
        std::cin >> from;

        std::cout << "Sweep doping to: ";
        std::cin >> to;

        std::cout << "How many data points? ";
        std::cin >> points;

        sweep_x(temperature, J1, J2, from, to, points);

        std::cout << "Doping sweep Done!" << std::endl;
    }
    if (option == "get_samples")
    {   std::cout << "Enter doping: ";
        std::cin >> doping;

        std::cout << "Enter desired number of samples: ";
        std::cin >> samples;

        std::cout << "Path to output directory: ";
        std::cin >> output_path;

        get_sample_configurations(temperature, doping, J1, J2, samples, output_path);
        std::cout << "Samples generated!" << std::endl;
    }

    return 0;
}
/*******************************************************************************************************************************************/
void get_sample_configurations(float temperature, float doping, float J1, float J2, int samples, std::string output_path)
{   // simulation parameters
    int thermalization_time = 200; // measured in MC steps per site
    int samples_taken = 30;        // number of time a value of a physical quantity is measured for averaging
    int sampling_time = 10;        // time[MC steps per site] to wait between measurements

    std::vector<Model> storage;
    for (int n = 0; n < samples; n++)
    {   storage.push_back(Model(temperature, doping, J1, J2));
    }

    // thermalizing samples in parallel
    #pragma omp parallel for
    for (int n = 0; n < samples; n++)
    {   storage[n].order_parameter(thermalization_time, samples_taken, sampling_time);
    }
    float op = 0;

    for (int n = 0; n < samples; n++)
    {   op += storage[n].order;
    }
    op = op / float(samples);

    // writing metadata to a file

    std::ofstream myfile(output_path + "_metadata.txt");

    // output simulation metadata
    myfile << "Temperature: " << temperature << std::endl;
    myfile << "J1: " << J1 << std::endl;
    myfile << "J2: " << J2 << std::endl;
    myfile << "Doping: " << doping << std::endl;
    myfile << "Number of samples: " << samples << std::endl;
    myfile << "Thermalization time [MC cycles]: " << thermalization_time << std::endl;
    myfile << "Samples taken for averaging: " << samples_taken << std::endl;
    myfile << "Sampling time [MC cycles]: " << sampling_time << std::endl;
    myfile << "Supercell size: " << size << std::endl;
    myfile << "Average order parameter: " << op << std::endl;
    myfile.close();
    // writing atom lists to .txt files for further analysis and simulations
    #pragma omp parallel for
    for (int n = 0; n < samples; n++)
    {   storage[n].output_structure(output_path +"_" + std::to_string(n) + ".csv", 5);
    }
}
/*******************************************************************************************************************************************/
void sweep_x(float temperature, float J1, float J2, float from, float to, int points)
{

    // simulation parameters
    int thermalization_time = 200; // measured in MC steps per site
    int samples_taken = 100;       // number of time a value of a physical quantity is measured for averaging
    int sampling_time = 20;        // time[MC steps per site] to wait between measurements

    float step = (to - from) / float(points);
    // creating an array of Models (each corresponding to different doping)
    std::vector<Model> storage;
    for (int n = 0; n < points; n++)
    {   storage.push_back(Model(temperature, from, J1, J2));
        from = from + step;
    }

// calculating the order parameter for different dopings in parallel
    #pragma omp parallel for
    for (int n = 0; n < points; n++)
    {   storage[n].order_parameter(thermalization_time, samples_taken, sampling_time);
    }
    // opening output files
    std::ofstream myfile("order_parameter.txt");
    std::ofstream myfile1("doping_list.txt");
    std::ofstream myfile2("metadata.txt");
    // output simulation metadata
    myfile2 << "Temperature: " << temperature << std::endl;
    myfile2 << "J1: " << J1 << std::endl;
    myfile2 << "J2: " << J2 << std::endl;
    myfile2 << "Doping swept from: " << from << std::endl;
    myfile2 << "Doping swept to: " << to << std::endl;
    myfile2 << "Number of data points: " << points << std::endl;
    myfile2 << "Thermalization time [MC cycles]: " << thermalization_time << std::endl;
    myfile2 << "Samples taken for averaging: " << samples_taken << std::endl;
    myfile2 << "Sampling time [MC cycles]: " << sampling_time << std::endl;
    myfile2 << "Supercell size: " << size << std::endl;

    // outputing data
    for (int n = 0; n < points; n++)
    {

        myfile1 << storage[n].get_x() << std::endl;
        myfile << storage[n].order << std::endl;
    }

    myfile.close();
    myfile1.close();
    myfile2.close();
}
/***************************************************************************************************************************************/

Model::Model(float temperature, float doping, float coupling1, float coupling2)
{   // initializing parameters
    this->t = temperature;
    this->x = doping;
    this->y = doping; // dictated by charge neutrality
    this->J1 = coupling1;
    this->J2 = coupling2;

    // helper variables
    int count1 = 0;
    int count2 = 0;
    float r;

    // populating the Model

    for (int i = 0; i < size; i++)
    {   for (int j = 0; j < size; j++)
        {   for (int k = 0; k < size; k++)
            {

                if (pow(-1, i + j + k) == 1)
                {   l[i][j][k][0] = 1;
                }
                else
                {   l[i][j][k][0] = 0;
                }

                r = dis(gen);
                if (r < this->x)
                {   l[i][j][k][1] = 1;
                    count1++;
                }
                else
                {   l[i][j][k][1] = 0;
                }
            }
        }
    }

    for (int i = 0; i < size; i++)
    {   for (int j = 0; j < size; j++)
        {   for (int k = 0; k < size; k++)
            {   if (l[i][j][k][0] == 0)
                {   r = dis(gen);
                    if (r < this->y)
                    {   l[i][j][k][0] = 1;
                    }
                }
                count2 = count2 + int(l[i][j][k][0]);
            }
        }
    }

    std::cout << "Model created!" << std::endl << "Temperature:  " << this->t << " 1/k_B" << std::endl;

    this->x_real = float(count1) / float(N);
    this->y_real = float(count2) / float(N);

    std::cout << "(Real) A-site composition ratio: " << float(count1) / float(N) << std::endl;
    std::cout << "(Real) B-site composition ratio: " << float(count2) / float(N) << std::endl;
}

void Model::set_glauber(bool a)
{   this->glauber_test = a;
}

void Model::set_temp(float T)
{   this->t = T;
}

void Model::set_couplings(float coupling1, float coupling2)
{   this->J1 = coupling1;
    this->J2 = coupling2;
}

void Model::set_x(float a)
{   this->x = a;
    int count1 = 0;
    int count2 = 0;

    float r = 0;
    // populating the Model
    for (int i = 0; i < size; i++)
    {   for (int j = 0; j < size; j++)
        {   for (int k = 0; k < size; k++)
            {

                if (pow(-1, i + j + k) == 1)
                {   l[i][j][k][0] = 1;
                }
                else
                {   l[i][j][k][0] = 0;
                }

                r = dis(gen);
                if (r < this->x)
                {   l[i][j][k][1] = 1;
                    count1++;
                }
                else
                {   l[i][j][k][1] = 0;
                }
            }
        }
    }

    for (int i = 0; i < size; i++)
    {   for (int j = 0; j < size; j++)
        {   for (int k = 0; k < size; k++)
            {   if (l[i][j][k][0] == 0)
                {   r = dis(gen);
                    if (r < (this->x))
                    {   l[i][j][k][0] = 1;
                    }
                }
                count2 = count2 + int(l[i][j][k][0]);
            }
        }
    }

    this->x_real = float(count1) / float(N);
    this->y_real = float(count2) / float(N);
}

float Model::get_x()
{   return this->x;
}

float Model::get_y()
{   return this->y_real;
}

int Model::get_s(int i, int j, int k)
{   if (l[i][j][k][0] == 1)
    {   return 1;
    }
    else
    {   return -1;
    }
}

int Model::get_sigma(int i, int j, int k)
{   if (l[i][j][k][1] == 1)
    {   return 1;
    }
    else
    {   return 0;
    }
}

// MONTE CARLO STEPS

void Model::swap_spins(int i, int j, int k, int u, int w, int v)
{   // std::swap(this->l[i][j][k][0], this->l[u][w][v][0]);
    bool a;
    bool b;
    a = this->l[i][j][k][0];
    b = this->l[u][w][v][0];

    this->l[i][j][k][0] = b;
    this->l[u][w][v][0] = a;
}

float Model::neighbor(int i, int j, int k, int site)
{   float sum = 0;
    // periodic boundary conditions
    int i_p = i == 0 ? size - 1 : i - 1;
    int i_n = i == size - 1 ? 0 : i + 1;
    int j_p = j == 0 ? size - 1 : j - 1;
    int j_n = j == size - 1 ? 0 : j + 1;
    int k_p = k == 0 ? size - 1 : k - 1;
    int k_n = k == size - 1 ? 0 : k + 1;

    if (site == 0)
    {   sum += this->get_s(i_p, j, k);
        sum += this->get_s(i_n, j, k);
        sum += this->get_s(i, j_p, k);
        sum += this->get_s(i, j_n, k);
        sum += this->get_s(i, j, k_p);
        sum += this->get_s(i, j, k_n);
        return sum;
    }
    else
    {   sum += this->get_sigma(i, j, k);
        sum += this->get_sigma(i_n, j, k);
        sum += this->get_sigma(i_n, j_n, k);
        sum += this->get_sigma(i, j_n, k);

        sum += this->get_sigma(i, j, k_n);
        sum += this->get_sigma(i_n, j, k_n);
        sum += this->get_sigma(i_n, j_n, k_n);
        sum += this->get_sigma(i, j_n, k_n);
        return sum;
    }
}

float Model::interaction(int i, int j, int k)
{   return this->get_s(i, j, k) * (this->J1 * this->neighbor(i, j, k, 0) + this->J2 * this->neighbor(i, j, k, 1));
}

bool Model::test(float E)
{

    float p = 0;
    float r = 0;

    if (this->glauber_test == 0)
    {

        if (E > 0)
        {   p = exp(-E / t);
            r = dis(gen);
            if (r < p)
            {   return 1;
            }
            else
            {   return 0;
            }
        }
        else
        {   return 1;
        }
    }
    else
    {   float p = 0;
        float r = 0;
        if (E > 0)
        {   p = exp(-E / t);
            p = p / (1 + p);
            r = dis(gen);
            if (r < p)
            {   return 1;
            }
            else
            {   return 0;
            }
        }
        else
        {   return 1;
        }
    }
}

int Model::update(int i, int j, int k, int u, int w, int v)
{   float E_initial, E_final, delta;
    E_initial = this->interaction(i, j, k) + this->interaction(u, v, w);
    this->swap_spins(i, j, k, u, v, w);
    E_final = this->interaction(i, j, k) + this->interaction(u, v, w);
    delta = E_final - E_initial;
    if (this->test(delta))
    {   if (delta < 0)
        {   return 1; // accepted good move
        }
        else if (delta > 0)
        {   return 2; // accepted bad move
        }
        else if (delta == 0)
        {   return 3; // accepted neutral
        }
    }
    this->swap_spins(i, j, k, u, v, w);
    return 0; // rejected move
}

void Model::sampling_cycle(int steps_per_site)
{   int i, j, k, u, v, w;
    int p = 0;

    int good = 0;
    int bad = 0;
    int neutral = 0;

    int success;

    do
    {   i = pick(gen);
        j = pick(gen);
        k = pick(gen);
        u = pick(gen);
        v = pick(gen);
        w = pick(gen);
        if (this->get_s(i, j, k) == this->get_s(u, v, w))
        {   continue;
        }
        else
        {   success = this->update(i, j, k, u, v, w);
            if (success == 1)
            {   good++;
            }
            else if (success == 2)
            {   bad++;
            }
            else if (success == 3)
            {   neutral++;
            }
        }
        p++;
    }
    while (p < (steps_per_site * N));

    // std::cout << "Good: " << good << "     Bad: " << bad << "    Neutral: " << neutral << std::endl;
    // std::cout << "Magnetization: " << this->magnetization() << std::endl;
}

// PHYSICAL QUANTITIES

float Model::magnetization()
{   float order = 0;

    for (int i = 0; i < size; i++)
    {   for (int j = 0; j < size; j++)
        {   for (int k = 0; k < size; k++)
            {   order += pow(-1, i + j + k) * this->get_s(i, j, k);
            }
        }
    }

    order = order / float(N);
    return order;
}

void Model::order_parameter(int thermalization_time, int samples_taken, int sampling_time)
{   float op = 0;
    this->sampling_cycle(thermalization_time);
    std::cout << "Temperature: " << this->t << "---> Thermalization complete" << std::endl;
    std::cout << "Doping: " << this->x_real << std::endl;

    for (int i = 0; i < samples_taken; i++)
    {   op += abs(this->magnetization());
        this->sampling_cycle(sampling_time);
    }

    // averaging
    op = op / float(samples_taken);
    this->order = op;
}

std::string Model::get_atom(int i, int j, int k, int n)
{   switch (n)
    {   case 4:
            if (this->get_sigma(i, j, k) == 0)
            {   return "Ba";
            }
            else
            {   return "K";
            }
            break;
        case 0:
            if (this->get_s(i, j, k) == -1)
            {   return "Hg";
            }
            else
            {   return "Pt";
            }
            break;
        case 1:
            return "O";
            break;
        case 2:
            return "O";
            break;
        case 3:
            return "O";
            break;
    }

    return "None";
}

void Model::output_structure(std::string filename, int atoms_per_cell)
{

    std::map<std::string, int> number;
    number["Ba"] = 56;
    number["K"] = 19;
    number["Pt"] = 78;
    number["Hg"] = 80;
    number["O"] = 8;

    std::map<int, std::string> pos;
    pos[0] = "0.0,0.0,0.0";
    pos[1] = "0.5,0.0,0.0";
    pos[2] = "0.0,0.5,0.0";
    pos[3] = "0.0,0.0,0.5";
    pos[4] = "0.5,0.5,0.5";

    std::ofstream myfile(filename);
    myfile << "i,j,k,site,Z,symbol,x,y,z" << std::endl;
    if (myfile.is_open())
    {   for (int i = 0; i < size; i++)
        {   for (int j = 0; j < size; j++)
            {   for (int k = 0; k < size; k++)
                {

                    for (int n = 0; n < atoms_per_cell; n++)
                    {   myfile << i << "," << j << "," << k << "," << n << "," << number[this->get_atom(i, j, k, n)]
                               << "," << this->get_atom(i, j, k, n) << "," << pos[n] << std::endl;
                    }
                }
            }
        }
    }
    else
    {   std::cout << "Unable to open file!" << std::endl;
    }

    myfile.close();
}

/****************************************************************************************************************************************/
