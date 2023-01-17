#include <cstdlib>
#include <iostream>
#include <map>
#include <math.h>
#include <memory>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <random>
#include <string>
#include <vector>

namespace py = pybind11;
using std::cout;
using std::endl;
using std::size_t;

/*Mersenne Twister random number generator */
std::random_device rd;
std::mt19937 gen(rd());

/*indexing pairs by a unique integer*/
int pairing(int a, int b)
{
    return (a + b) * (a + b + 1) + 2 * b;
}

/*A class to store and perform simple operations on 3D vectors*/
class Vector
{
    double v[3];

  public:
    /*Constructors*/
    Vector();
    Vector(double, double, double);

    void set(double, double, double);
    void set_comp(int, double);
    double get(int);
    double norm();
    double norm_squared();
    void print();
    std::vector<double> to_vector();
    void mul(double);
    void show();

    Vector operator+(const Vector &b)
    {
        Vector v;
        for (int i = 0; i < 3; i++)
        {
            v.v[i] = this->v[i] + b.v[i];
        }
        return v;
    }

    Vector operator*(const double a)
    {
        Vector v;
        for (int i = 0; i < 3; i++)
        {
            v.v[i] = this->v[i] * a;
        }
        return v;
    }

    Vector operator-(const Vector &b)
    {
        Vector v;
        for (int i = 0; i < 3; i++)
        {
            v.v[i] = this->v[i] - b.v[i];
        }
        return v;
    }

    double operator*(const Vector &b)
    {
        double dot = 0;
        for (int i = 0; i < 3; i++)
        {
            dot += this->v[i] * b.v[i];
        }
        return dot;
    }

    void operator=(const Vector &b)
    {
        for (int i = 0; i < 3; i++)
        {
            this->v[i] = b.v[i];
        }
    }
};

// implementation of methods
Vector::Vector()
{
    this->v[0] = 0;
    this->v[1] = 0;
    this->v[2] = 0;
}

Vector::Vector(double i, double j, double k)
{
    this->v[0] = i;
    this->v[1] = j;
    this->v[2] = k;
}

void Vector::mul(double factor)
{
    this->v[0] = this->v[0] * factor;
    this->v[1] = this->v[1] * factor;
    this->v[2] = this->v[2] * factor;
}

void Vector::set(double i, double j, double k)
{
    this->v[0] = i;
    this->v[1] = j;
    this->v[2] = k;
}

void Vector::set_comp(int i, double value)
{
    this->v[i] = value;
}

double Vector::get(int i)
{
    return this->v[i];
}

double Vector::norm()
{
    double n = 0;
    for (int i = 0; i < 3; i++)
    {
        n += this->v[i] * this->v[i];
    }
    return sqrt(n);
}

double Vector::norm_squared()
{
    double n = 0;
    for (int i = 0; i < 3; i++)
    {
        n += this->v[i] * this->v[i];
    }
    return n;
}

void Vector::show()
{
    py::print(py::make_tuple(this->v[0], this->v[1], this->v[2]));
}

std::vector<double> Vector::to_vector()
{
    std::vector<double> temp;
    temp.push_back(v[0]);
    temp.push_back(v[1]);
    temp.push_back(v[2]);
    return temp;
}

/*----------------------------------------------------------------------------------------------------------------------*/

/* A collection of different empirical potentials and distance dependent
   constraints that can potentially be used as penalty functions during
   constrained optimization of atomic valences as calculated from Brown's
   rules. All instances are inherited from the abstract class "Interaction" and
   override its default evaluate method. */
class Interaction
{
  public:
    std::vector<double> p;
    Interaction(py::array_t<double> par)
    {
        auto a = par.mutable_unchecked<1>();
        for (py::ssize_t i = 0; i < a.shape(0); i++)
        {
            this->p.push_back(a(i));
        }
    }

    virtual ~Interaction() = default;
    virtual double evaluate(double r) = 0;
};

/*Derived classes which override the evaluate method */

class SpringInteraction : public Interaction
{
  public:
    using Interaction::Interaction;

    double evaluate(double r) override
    {
        r = r - this->p[1];
        return this->p[0] * r * r;
    }
};

class MorseInteraction : public Interaction
{
  public:
    using Interaction::Interaction;

    double evaluate(double r) override
    {
        double temp;
        temp = (1 - exp((r - p[2]) / p[1]));
        return p[0] * temp * temp;
    }
};

class RepulsiveInteraction : public Interaction
{
  public:
    using Interaction::Interaction;
    double evaluate(double r) override
    {
        return pow((p[0] / r), p[1]);
    }
};

class LennardJonesInteraction : public Interaction
{
  public:
    using Interaction::Interaction;
    ~LennardJonesInteraction()
    {
    }

    double evaluate(double r) override
    {
        return p[0] * (pow((p[1] / r), 12) - pow((p[1] / r), 6));
    }
};

class YukawaInteraction : public Interaction
{
  public:
    using Interaction::Interaction;

    double evaluate(double r) override
    {
        return p[0] * p[1] * exp(-r / p[2]) / r;
    }
};

class TruncatedCoulombInteraction : public Interaction
{
  public:
    using Interaction::Interaction;

    double evaluate(double r) override
    {
        return p[0] * p[1] * erfc(r / p[2]) / r;
    }
};

class VolumeExclusion : public Interaction
{
    using Interaction::Interaction;

    double evaluate(double r) override
    {
        if (r > (p[0] + p[1]))
            return 0;
        else
            return 1e10;
    }
};

class SoftVolumeExclusion : public Interaction
{
    using Interaction::Interaction;

    double evaluate(double r) override
    {
        return p[0] * pow(p[1] / r, 12);
    }
};

/*-----------------------------------------------------------------------------------------------------------------------*/
class structure
{
  private:
    /*Primitive lattice vectors*/
    Vector a1, a2, a3;

    /*Average atom positions in the cell*/
    // std::map<int, Vector> basis;

    /*Mapping between sites and their bonds*/
    std::map<int, std::vector<std::vector<int>>> bonds;
    /*Bond valence sum parameters -> indexed by cantor pairing*/
    std::map<int, double> R_0;
    std::map<int, double> b;
    /*Nominal valence of atoms*/
    std::map<int, double> nominal;
    std::map<int, double> nominal_vector_valence;
    std::map<int, double> weight;
    std::map<int, double> vector_weight;
    /*Methods used only inside the class*/
    void load_atoms(py::array_t<int>);
    void load_positions(py::array_t<double>);

  public:
    /*Array of atom positions with respect to the cell origin.
     Given in fractional coordinates*/
    std::vector<Vector> positions;
    /*Array of atomic numbers*/
    std::vector<int> atoms;
    /*Number of unit cells in each direction */
    int shape;
    /*Number of sites per unit cell*/
    int sites;

    /*constructor*/
    structure(int, int, py::array_t<double>);

    /*Load atomic numbers and positions*/
    void load(py::array_t<int>, py::array_t<double>);

    /*Getters and setters -> They use indexing to treat 1D array as a 3D array of
      unit cells. Automatic implementation of periodic boundary conditions*/
    int index(int, int, int, int);
    int get_atom(int, int, int, int);
    void set_position(double, double, double);
    Vector get_position(int, int, int, int);
    Vector get_reduced_position(int, int, int, int);
    Vector get_real_position(int, int, int, int);
    Vector get_scaled_position(int, int, int, int);

    void set_bonds(int, py::array_t<int>);
    void set_parameters(int, int, double, double);
    void set_nominals_and_weights(py::array_t<int>, py::array_t<double>, py::array_t<double>, py::array_t<double>,
                                  py::array_t<double>);

    /*Return atoms and positions to python*/
    std::vector<int> get_atoms();
    std::vector<std::vector<double>> get_positions();
    std::vector<std::vector<double>> get_cell();
    std::vector<int> get_shape();
    /*Perturb the positions by a small amount drawn randomly from the normal
     * distribution*/
    void agitate(int, double, double, double);

    /*Get properties*/
    double compute_valence(int, int, int, int);
    double compute_valence_vector(int, int, int, int);
    std::vector<double> valence_list(int, int);
    double average_valence(int);
    double valence_deviation(int, int, int, int);
    double valence_vector_deviation(int, int, int, int);
    double local_instability_index(int, int, int, int);
    double local_valence_vector_change(int, int, int, int);
    double global_instability_index();
    std::vector<double> rms_displacement(int);
    std::vector<double> atomic_displacement_parameters(int);
    double average_distance(int, int, py::array_t<int>);
    std::vector<double> distance_list(int, int, py::array_t<int>);
    double local_dipole_moment(py::array_t<int>);
    double CS_parameter(py::array_t<int>);
};

/*Helper function for map initializaton*/
std::vector<std::vector<int>> return_empty_vector()
{
    std::vector<std::vector<int>> v;
    return v;
}

/*Method implementation */

structure::structure(int ncells, int sites_per_cell, py::array_t<double> primitive_latt_vectors)
{
    this->shape = ncells;
    this->sites = sites_per_cell;

    auto r = primitive_latt_vectors.mutable_unchecked<2>();
    for (py::ssize_t i = 0; i < 3; i++)
    {
        for (py::ssize_t j = 0; j < 3; j++)
        {
            if (r(i, j) < pow(10, -5))
            {
                r(i, j) = 0;
            }
        }
    }
    this->a1.set(r(0, 0), r(0, 1), r(0, 2));
    this->a2.set(r(1, 0), r(1, 1), r(1, 2));
    this->a3.set(r(2, 0), r(2, 1), r(2, 2));
}

void structure::load_atoms(py::array_t<int> x)
{
    auto r = x.mutable_unchecked<1>();
    for (py::ssize_t i = 0; i < r.shape(0); i++)
    {
        this->atoms.push_back(r(i));
    }
}

void structure::load_positions(py::array_t<double> x)
{
    auto r = x.mutable_unchecked<2>();
    Vector temp;
    for (py::ssize_t i = 0; i < r.shape(0); i++)
    {
        temp.set(r(i, 0), r(i, 1), r(i, 2));
        this->positions.push_back(temp);
    }
}

void structure::load(py::array_t<int> x, py::array_t<double> y)
{
    this->atoms.clear();
    this->positions.clear();
    this->load_atoms(x);
    this->load_positions(y);

    if (int(this->atoms.size()) != this->sites * this->shape * this->shape * this->shape)
    {
        throw "Size of the given array does not match the prescribed size of "
              "the "
              "structure!";
    }
}

/*3D indexing of a 1D array*/
int structure::index(int i, int j, int k, int s)
{
    return i * this->shape * this->shape * this->sites + j * this->shape * this->sites + k * this->sites + s;
}

/*Get atomic number of an atom at a specific site*/
int structure::get_atom(int i, int j, int k, int s)
{
    int I, J, K;
    I = i < 0 ? i + this->shape : i;
    i = I;
    I = i > (this->shape - 1) ? i - this->shape : i;
    J = j < 0 ? j + this->shape : j;
    j = J;
    J = j > (this->shape - 1) ? j - this->shape : j;
    K = k < 0 ? k + this->shape : k;
    k = K;
    K = k > (this->shape - 1) ? k - this->shape : k;

    return this->atoms[this->index(I, J, K, s)];
}

/*Get position in lattice units*/
Vector structure::get_position(int i, int j, int k, int s)
{
    int I, J, K;

    I = i < 0 ? i + this->shape : i;
    i = I;
    I = i > (this->shape - 1) ? i - this->shape : i;
    J = j < 0 ? j + this->shape : j;
    j = J;
    J = j > (this->shape - 1) ? j - this->shape : j;
    K = k < 0 ? k + this->shape : k;
    k = K;
    K = k > (this->shape - 1) ? k - this->shape : k;

    return this->positions[this->index(I, J, K, s)];
}

/*Get atom position in Cartesian coordinates*/
Vector structure::get_real_position(int i, int j, int k, int s)
{
    Vector scaled, real_pos;
    scaled = this->get_position(i, j, k, s);
    real_pos = this->a1 * scaled.get(0) + this->a2 * scaled.get(1) + this->a3 * scaled.get(2);
    real_pos = real_pos + this->a1 * double(i) + this->a2 * double(j) + this->a3 * double(k);
    return real_pos;
}

Vector structure::get_scaled_position(int i, int j, int k, int s)
{
    Vector scaled;
    scaled.set(i, j, k);
    return scaled + this->get_position(i, j, k, s);
}

Vector structure::get_reduced_position(int i, int j, int k, int s)
{
    Vector v = this->get_position(i, j, k, s);
    return this->a1 * v.get(0) + this->a2 * v.get(1) + this->a3 * v.get(2);
}

/* Methods used to retrieve structure in python after the simulation is done*/
std::vector<int> structure::get_atoms()
{
    return this->atoms;
}
std::vector<int> structure::get_shape()
{
    std::vector<int> shape;
    shape.push_back(this->shape);
    shape.push_back(this->shape);
    shape.push_back(this->shape);
    shape.push_back(this->sites);
    return shape;
}

std::vector<std::vector<double>> structure::get_cell()
{
    std::vector<std::vector<double>> temp;
    temp.push_back(a1.to_vector());
    temp.push_back(a2.to_vector());
    temp.push_back(a3.to_vector());
    return temp;
}

std::vector<std::vector<double>> structure::get_positions()
{
    std::vector<std::vector<double>> temp;
    for (Vector &pos : this->positions)
    {
        temp.push_back({pos.get(0), pos.get(1), pos.get(2)});
    }
    return temp;
}

void structure::set_bonds(int site, py::array_t<int> x)
{
    auto r = x.mutable_unchecked<2>();
    std::vector<int> temp;
    this->bonds[site].clear();
    for (py::ssize_t i = 0; i < r.shape(0); i++)
    {
        temp.clear();
        for (py::ssize_t j = 0; j < r.shape(1); j++)
        {
            temp.push_back(r(i, j));
        }
        if (temp.size() != 5)
        {
            throw "Invalid neighbor!";
        }
        this->bonds[site].push_back(temp);
    }
}

void structure::set_nominals_and_weights(py::array_t<int> atoms, py::array_t<double> nominals,
                                         py::array_t<double> vector_nominals, py::array_t<double> weights,
                                         py::array_t<double> vector_weights)
{
    auto a = atoms.mutable_unchecked<1>();
    auto n = nominals.mutable_unchecked<1>();
    auto vn = vector_nominals.mutable_unchecked<1>();
    auto w = weights.mutable_unchecked<1>();
    auto vw = weights.mutable_unchecked<1>();

    for (unsigned long int i = 0; i < sizeof(a); i++)
    {
        this->nominal[a(i)] = n(i);
        this->nominal_vector_valence[a(i)] = vn(i);
        this->weight[a(i)] = w(i);
        this->vector_weight[a(i)] = vw(i);
    }
}

void structure::set_parameters(int atom1, int atom2, double R, double B)
{
    this->R_0[pairing(atom1, atom2)] = R;
    this->R_0[pairing(atom2, atom1)] = R;

    this->b[pairing(atom1, atom2)] = B;
    this->b[pairing(atom2, atom1)] = B;
}

/*Displace all atoms at a specfic site by a random amount drawn from 3D
 * Gaussian distribution*/
void structure::agitate(int site, double s11, double s22, double s33)
{
    if (site < 0 || site >= this->sites)
    {
        py::print("Invalid site!");
        return;
    }

    std::normal_distribution<double> gauss1(0, s11);
    std::normal_distribution<double> gauss2(0, s22);
    std::normal_distribution<double> gauss3(0, s33);

    Vector disp;
    for (int i = 0; i < this->shape; i++)
    {
        for (int j = 0; j < this->shape; j++)
        {
            for (int k = 0; k < this->shape; k++)
            {
                disp.set(gauss1(gen), gauss2(gen), gauss3(gen));
                this->positions[this->index(i, j, k, site)] = this->positions[this->index(i, j, k, site)] + disp;
            }
        }
    }
}

double structure::compute_valence(int i, int j, int k, int s)
{
    int atom1, atom2;
    Vector pos1, pos2;
    atom1 = this->get_atom(i, j, k, s);
    pos1 = this->get_real_position(i, j, k, s);
    double valence = 0;
    int pair;

    for (const std::vector<int> &v : this->bonds[s])
    {
        pos2 = this->get_real_position(i + v[2], j + v[3], k + v[4], v[1]);
        atom2 = this->get_atom(i + v[2], j + v[3], k + v[4], v[1]);
        pair = pairing(atom1, atom2);
        valence += exp((this->R_0[pair] - (pos2 - pos1).norm()) / this->b[pair]);
    }

    return valence;
}

double structure::compute_valence_vector(int i, int j, int k, int s)
{
    int atom1, atom2, pair;
    double r;
    Vector pos1, pos2, bond_vector, valence_vector;
    atom1 = this->get_atom(i, j, k, s);
    pos1 = this->get_real_position(i, j, k, s);
    valence_vector.set(0, 0, 0);

    for (const std::vector<int> &v : this->bonds[s])
    {
        pos2 = this->get_real_position(i + v[2], j + v[3], k + v[4], v[1]);
        atom2 = this->get_atom(i + v[2], j + v[3], k + v[4], v[1]);
        pair = pairing(atom1, atom2);
        r = (pos2 - pos1).norm();
        bond_vector = pos2 - pos1;
        bond_vector.mul(1 / r);
        bond_vector.mul(exp((this->R_0[pair] - r) / this->b[pair]));
        valence_vector = valence_vector + bond_vector;
    }

    return valence_vector.norm();
}

double structure::valence_deviation(int i, int j, int k, int s)
{
    int atom = this->get_atom(i, j, k, s);
    double difference = this->compute_valence(i, j, k, s) - this->nominal[atom];
    return this->weight[atom] * difference * difference;
}

double structure::valence_vector_deviation(int i, int j, int k, int s)
{
    int atom = this->get_atom(i, j, k, s);
    double difference = this->compute_valence_vector(i, j, k, s) - this->nominal_vector_valence[atom];
    return this->vector_weight[atom] * difference * difference;
}

double structure::local_instability_index(int i, int j, int k, int s)
{
    double lii = this->valence_deviation(i, j, k, s);
    for (std::vector<int> n : this->bonds[s])
    {
        lii += this->valence_deviation(i + n[2], j + n[3], k + n[4], n[1]);
    }
    return lii;
}

double structure::local_valence_vector_change(int i, int j, int k, int s)
{
    double lvvc = this->valence_vector_deviation(i, j, k, s);
    for (std::vector<int> n : this->bonds[s])
    {
        lvvc += this->valence_vector_deviation(i + n[2], j + n[3], k + n[4], n[1]);
    }
    return lvvc;
}

/*Data analysis methods*/
std::vector<double> structure::valence_list(int site, int atom)
{
    std::vector<double> d;
    for (int i = 0; i < this->shape; i++)
    {
        for (int j = 0; j < this->shape; j++)
        {
            for (int k = 0; k < this->shape; k++)
            {
                if (this->get_atom(i, j, k, site) == atom)
                {
                    d.push_back(this->compute_valence(i, j, k, site));
                }
            }
        }
    }
    return d;
}

double structure::average_valence(int atom)
{
    double temp = 0;
    int count = 0;
    for (int i = 0; i < this->shape; i++)
    {
        for (int j = 0; j < this->shape; j++)
        {
            for (int k = 0; k < this->shape; k++)
            {
                for (int s = 0; s < this->sites; s++)
                {
                    if (this->get_atom(i, j, k, s) == atom)
                    {
                        temp += this->compute_valence(i, j, k, s);
                        count++;
                    }
                }
            }
        }
    }
    return temp / double(count);
}

double structure::global_instability_index()
{
    double temp = 0;
    int atom;
    double diff;
    for (int i = 0; i < this->shape; i++)
    {
        for (int j = 0; j < this->shape; j++)
        {
            for (int k = 0; k < this->shape; k++)
            {
                for (int s = 0; s < this->sites; s++)
                {
                    atom = this->get_atom(i, j, k, s);
                    diff = (this->compute_valence(i, j, k, s) - this->nominal[atom]);
                    temp += this->weight[atom] * diff * diff;
                }
            }
        }
    }
    return temp;
}

std::vector<double> structure::rms_displacement(int s)
{
    std::vector<double> std;
    int count = 0;
    Vector av;
    double sigma_x = 0;
    double sigma_y = 0;
    double sigma_z = 0;
    double av_x = 0;
    double av_y = 0;
    double av_z = 0;
    for (int i = 0; i < this->shape; i++)
    {
        for (int j = 0; j < this->shape; j++)
        {
            for (int k = 0; k < this->shape; k++)
            {
                av_x += this->get_position(i, j, k, s).get(0);
                av_y += this->get_position(i, j, k, s).get(1);
                av_z += this->get_position(i, j, k, s).get(2);

                sigma_x += pow(this->get_position(i, j, k, s).get(0), 2);
                sigma_y += pow(this->get_position(i, j, k, s).get(1), 2);
                sigma_z += pow(this->get_position(i, j, k, s).get(2), 2);
                count++;
            }
        }
    }
    av_x = av_x / double(count);
    av_y = av_y / double(count);
    av_z = av_z / double(count);

    sigma_x = std::sqrt(sigma_x / double(count) - av_x * av_x);
    sigma_y = std::sqrt(sigma_y / double(count) - av_y * av_y);
    sigma_z = std::sqrt(sigma_z / double(count) - av_z * av_z);

    std.push_back(sigma_x);
    std.push_back(sigma_y);
    std.push_back(sigma_z);
    return std;
}

/*Anisotropic B parameters*/
std::vector<double> structure::atomic_displacement_parameters(int s)
{
    std::vector<double> std;
    int count = 0;
    Vector av;
    double sigma_x = 0;
    double sigma_y = 0;
    double sigma_z = 0;
    double av_x = 0;
    double av_y = 0;
    double av_z = 0;
    for (int i = 0; i < this->shape; i++)
    {
        for (int j = 0; j < this->shape; j++)
        {
            for (int k = 0; k < this->shape; k++)
            {
                av_x += this->get_reduced_position(i, j, k, s).get(0);
                av_y += this->get_reduced_position(i, j, k, s).get(1);
                av_z += this->get_reduced_position(i, j, k, s).get(2);

                sigma_x += pow(this->get_reduced_position(i, j, k, s).get(0), 2);
                sigma_y += pow(this->get_reduced_position(i, j, k, s).get(1), 2);
                sigma_z += pow(this->get_reduced_position(i, j, k, s).get(2), 2);
                count++;
            }
        }
    }
    av_x = av_x / double(count);
    av_y = av_y / double(count);
    av_z = av_z / double(count);

    /*Anisotropic U parameters */
    sigma_x = sigma_x / double(count) - av_x * av_x;
    sigma_y = sigma_y / double(count) - av_y * av_y;
    sigma_z = sigma_z / double(count) - av_z * av_z;

    std.push_back(8 * M_PI * M_PI * sigma_x);
    std.push_back(8 * M_PI * M_PI * sigma_y);
    std.push_back(8 * M_PI * M_PI * sigma_z);
    return std;
}

std::vector<double> structure::distance_list(int atom1, int atom2, py::array_t<int> neighbors)
{
    std::vector<double> distances;
    auto r = neighbors.mutable_unchecked<2>();
    std::vector<int> temp;
    int a1, a2;
    Vector pos1, pos2;

    for (py::ssize_t i = 0; i < r.shape(0); i++)
    {
        temp.clear();
        for (py::ssize_t j = 0; j < r.shape(1); j++)
        {
            temp.push_back(r(i, j));
        }
        if (temp.size() != 5)
        {
            throw "Invalid neighbor!";
        }

        for (int l = 0; l < this->shape; l++)
        {
            for (int m = 0; m < this->shape; m++)
            {
                for (int n = 0; n < this->shape; n++)
                {
                    a1 = this->get_atom(l, m, n, temp[0]);
                    pos1 = this->get_real_position(l, m, n, temp[0]);
                    a2 = this->get_atom(l + temp[2], m + temp[3], n + temp[4], temp[1]);
                    pos2 = this->get_real_position(l + temp[2], m + temp[3], n + temp[4], temp[1]);

                    if (a1 == atom1 && a2 == atom2)
                    {
                        distances.push_back((pos2 - pos1).norm());
                    }
                    if (a2 == atom1 && a1 == atom2)
                    {
                        distances.push_back((pos2 - pos1).norm());
                    }
                }
            }
        }
    }
    return distances;
}

double structure::average_distance(int atom1, int atom2, py::array_t<int> neighbors)
{
    std::vector<double> distances;
    distances = this->distance_list(atom1, atom2, neighbors);
    double temp = 0;
    int count = 0;
    for (auto dist : distances)
    {
        temp += dist;
        count++;
    }
    return temp / double(count);
}

double structure::local_dipole_moment(py::array_t<int> neig)
{
    auto r = neig.mutable_unchecked<2>();
    std::vector<std::vector<int>> neighbors;
    Vector pos1, pos2, dipole;
    std::vector<int> temp;
    double result = 0;
    int count = 0;

    for (py::ssize_t i = 0; i < r.shape(0); i++)
    {
        temp.clear();
        for (py::ssize_t j = 0; j < r.shape(1); j++)
        {
            temp.push_back(r(i, j));
        }
        if (temp.size() != 5)
        {
            throw "Invalid neighbor!";
        }
        neighbors.push_back(temp);
    }
    for (int i = 0; i < this->shape; i++)
    {
        for (int j = 0; j < this->shape; j++)
        {
            for (int k = 0; k < this->shape; k++)
            {
                dipole.set(0, 0, 0);
                pos1 = this->get_scaled_position(i, j, k, neighbors[0][0]);
                for (auto &n : neighbors)
                {
                    pos2 = this->get_scaled_position(i + n[2], j + n[3], k + n[4], n[1]);
                    dipole = dipole + (pos2 - pos1);
                }
                result += dipole.norm();
                count++;
            }
        }
    }

    return result / double(count);
}

double structure::CS_parameter(py::array_t<int> neig)
{
    auto r = neig.mutable_unchecked<2>();
    std::vector<std::vector<int>> neighbors;
    Vector pos1, pos2, pos3;
    std::vector<int> temp;
    double CS;
    double result = 0;
    int count = 0;

    for (py::ssize_t i = 0; i < r.shape(0); i++)
    {
        temp.clear();
        for (py::ssize_t j = 0; j < r.shape(1); j++)
        {
            temp.push_back(r(i, j));
        }
        if (temp.size() != 5)
        {
            throw "Invalid neighbor!";
        }
        neighbors.push_back(temp);
    }
    for (int i = 0; i < this->shape; i++)
    {
        for (int j = 0; j < this->shape; j++)
        {
            for (int k = 0; k < this->shape; k++)
            {
                CS = 0;
                for (size_t l = 0; l < neighbors.size(); l = l + 2)
                {
                    pos1 = this->get_scaled_position(i, j, k, neighbors[l][0]);
                    pos2 = this->get_scaled_position(i + neighbors[l][2], j + neighbors[l][3], k + neighbors[l][4],
                                                     neighbors[l][1]);
                    pos3 = this->get_scaled_position(i + neighbors[l + 1][2], j + neighbors[l + 1][3],
                                                     k + neighbors[l + 1][4], neighbors[l + 1][1]);

                    CS += ((pos2 - pos1) + (pos3 - pos1)).norm_squared();
                }
                result += CS;
                count++;
            }
        }
    }

    return result / double(count);
}

/*-----------------------------------------------------------------------------------------------------------------------*/
/*helper function*/
std::normal_distribution<double> gaussian(double i, double j)
{
    std::normal_distribution<double> g(i, j);
    return g;
}
/* A class used to perform the elementary steps of the Monte Carlo algorithm.
   For now, the only implemented generator is the one that proposes
   displacements drawn from the normal distribution. In the future, plan is to
   make this class abstract and use polymorphism to implement several other
   useful generators.
*/

class Generator
{
  private:
    int i, j, k;
    Vector temp, storage;
    std::normal_distribution<double> gauss1;
    std::normal_distribution<double> gauss2;
    std::normal_distribution<double> gauss3;

  public:
    int site;
    Generator(int, py::array_t<double>, py::array_t<double>);
    void reset_rng();
    void modify(structure &, int, int, int);
    void undo_last(structure &);
};

Generator::Generator(int s, py::array_t<double> mean, py::array_t<double> std_dev)
{
    this->site = s;
    auto m = mean.mutable_unchecked<1>();
    auto std = std_dev.mutable_unchecked<1>();

    this->gauss1 = gaussian(m(0), std(0));
    this->gauss2 = gaussian(m(1), std(1));
    this->gauss3 = gaussian(m(2), std(2));
}

void Generator::modify(structure &s, int x, int y, int z)
{
    this->i = x;
    this->j = y;
    this->k = z;
    this->temp.set(this->gauss1(gen), this->gauss2(gen), this->gauss3(gen));
    this->storage = s.positions[s.index(i, j, k, site)];
    s.positions[s.index(i, j, k, site)] = this->temp;
}

void Generator::undo_last(structure &s)
{
    s.positions[s.index(i, j, k, site)] = this->storage;
}

void Generator::reset_rng()
{
    this->gauss1.reset();
    this->gauss2.reset();
    this->gauss3.reset();
}

/*-------------------------------------------------------------------------------------------------------------------------*/
/*A class that stores neighbors for a particular site and computes interaction
energies. It is mainly used to include
non-bonded interactions and to constrain the valence deviation optimization*/
class Calculator
{
  private:
    /*Neighbors*/
    std::vector<std::vector<int>> neighbors;
    /*Pair interactions which contribute -> indexed by atomic pairs using cantor
     * map */
    std::map<int, Interaction *> Interactions;

  public:
    /*Site for which to compute interactions*/
    int site;
    /*Constructor*/
    Calculator(int st, py::array_t<int> neig)
    {
        this->site = st;
        auto r = neig.mutable_unchecked<2>();
        std::vector<int> temp;
        for (py::ssize_t i = 0; i < r.shape(0); i++)
        {
            temp.clear();
            for (py::ssize_t j = 0; j < r.shape(1); j++)
            {
                temp.push_back(r(i, j));
            }
            if (temp.size() != 5)
            {
                throw "Invalid neighbor!";
            }
            this->neighbors.push_back(temp);
        }
    }

    void add_interaction(int atom1, int atom2, std::string interaction_type, py::array_t<double> p)
    {
        if (interaction_type == "Spring")
        {
            this->Interactions.insert(std::pair<int, Interaction *>(pairing(atom1, atom2), new SpringInteraction(p)));
            this->Interactions.insert(std::pair<int, Interaction *>(pairing(atom2, atom1), new SpringInteraction(p)));
        }
        else if (interaction_type == "LennardJones")
        {
            this->Interactions.insert(
                std::pair<int, Interaction *>(pairing(atom1, atom2), new LennardJonesInteraction(p)));
            this->Interactions.insert(
                std::pair<int, Interaction *>(pairing(atom2, atom1), new LennardJonesInteraction(p)));
        }
        else if (interaction_type == "Yukawa")
        {
            this->Interactions.insert(std::pair<int, Interaction *>(pairing(atom1, atom2), new YukawaInteraction(p)));
            this->Interactions.insert(std::pair<int, Interaction *>(pairing(atom2, atom1), new YukawaInteraction(p)));
        }
        else if (interaction_type == "TruncatedCoulomb")
        {
            this->Interactions.insert(
                std::pair<int, Interaction *>(pairing(atom1, atom2), new TruncatedCoulombInteraction(p)));
            this->Interactions.insert(
                std::pair<int, Interaction *>(pairing(atom2, atom1), new TruncatedCoulombInteraction(p)));
        }
        else if (interaction_type == "Repulsive")
        {
            this->Interactions.insert(
                std::pair<int, Interaction *>(pairing(atom1, atom2), new RepulsiveInteraction(p)));
            this->Interactions.insert(
                std::pair<int, Interaction *>(pairing(atom2, atom1), new RepulsiveInteraction(p)));
        }
        else if (interaction_type == "Morse")
        {
            this->Interactions.insert(std::pair<int, Interaction *>(pairing(atom1, atom2), new MorseInteraction(p)));
            this->Interactions.insert(std::pair<int, Interaction *>(pairing(atom2, atom1), new MorseInteraction(p)));
        }
        else if (interaction_type == "VolumeExclusion")
        {
            this->Interactions.insert(std::pair<int, Interaction *>(pairing(atom1, atom2), new VolumeExclusion(p)));
            this->Interactions.insert(std::pair<int, Interaction *>(pairing(atom2, atom1), new VolumeExclusion(p)));
        }
        else if (interaction_type == "SoftVolumeExclusion")
        {
            this->Interactions.insert(std::pair<int, Interaction *>(pairing(atom1, atom2), new SoftVolumeExclusion(p)));
            this->Interactions.insert(std::pair<int, Interaction *>(pairing(atom2, atom1), new SoftVolumeExclusion(p)));
        }
        else
        {
            throw "Unknown interaction type!";
        }

        return;
    }

    double execute(structure &s, int i, int j, int k)
    {
        double result = 0;

        /*Getting the species and position of the modified atom*/
        int atom1 = s.get_atom(i, j, k, this->site);
        Vector pos1 = s.get_real_position(i, j, k, this->site);

        int atom2;
        Vector pos2;

        int pair_index;

        /*Iterating over neighbors*/
        for (std::vector<int> &n : this->neighbors)
        {
            /*Getting the species and position of the neighbor*/
            atom2 = s.get_atom(i + n[2], j + n[3], k + n[4], n[1]);

            pos2 = s.get_real_position(i + n[2], j + n[3], k + n[4], n[1]);

            /*Calculate contribution to interaction energy*/
            pair_index = pairing(atom1, atom2);
            if (this->Interactions.find(pair_index) != this->Interactions.end())
            {
                result += this->Interactions[pair_index]->evaluate(((pos2 - pos1).norm()));
            }
            else
            {
                py::print("Unspecified interaction encountered!");
                throw("Error");
            }
        }

        return result;
    }

    ~Calculator()
    {
        for (std::pair<int, Interaction *> p : this->Interactions)
        {
            delete p.second;
        }
    }
};

/*-------------------------------------------------------------------------------------------------------------------------*/

/* Monte Carlo engine that implements the Metropolis-Hastings algorithm to
   optimize the global istability index subject to given hard and soft
   constraints given by the penalty function.*/
class simulation
{
  private:
    std::vector<Generator *> generators;
    std::vector<Calculator *> calculators;

  public:
    bool optimize_valence = true;
    bool optimize_vector_valence = true;
    simulation()
    {
        return;
    }
    void add_generator(Generator &m)
    {
        for (auto mod : this->generators)
        {
            if (mod == std::addressof(m))
            {
                py::print("Duplication detected! Modifier not added.");
                return;
            }
        }
        this->generators.push_back(&m);
    }
    void add_calculator(Calculator &c)
    {
        for (auto cal : this->calculators)
        {
            if (cal == std::addressof(c))
            {
                py::print("Duplication detected! Calculator not added.");
                return;
            }
        }
        this->calculators.push_back(&c);
    }

    void run_cycle(structure &, double, double, double);
};

void simulation::run_cycle(structure &s, double T1, double T2, double T3)
{
    /*Iterations per cycle = number of unit cells in a crystal */
    int iterations = int(pow(s.shape, 3));
    /*Random number generators*/
    std::uniform_int_distribution<int> pick(0, s.shape - 1);
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    int i, j, k, st;

    int good = 0;
    int neutral = 0;
    int bad = 0;

    double bv_cost1, bv_cost2, bvv_cost1, bvv_cost2, constraint1, constraint2, change;
    double p, boltzmann;
    /*Looping over iterations (number of proposed MC moves per cycle).
    By default equal to number of unit cells in the simulated crystal */
    for (int n = 0; n < iterations; n++)
    {
        /*Looping over modifiers*/
        for (auto &mod : this->generators)
        {
            /*Pick the unit cell randomly*/
            i = pick(gen);
            j = pick(gen);
            k = pick(gen);
            st = mod->site;

            /*Calculate the energy before the MC move*/
            bv_cost1 = 0;
            bvv_cost1 = 0;
            constraint1 = 0;
            /*Iterating over calculators*/
            for (auto &calc : this->calculators)
            {
                /*Compare calculator site to modified site*/
                if (calc->site == st)
                {
                    constraint1 += calc->execute(s, i, j, k);
                }
            }
            if (optimize_valence == true)
            {
                bv_cost1 = s.local_instability_index(i, j, k, st);
            }
            if (optimize_vector_valence == true)
            {
                bvv_cost1 = s.local_valence_vector_change(i, j, k, st);
            }
            /*Modify the structure s*/
            mod->modify(s, i, j, k);

            /*Calculate the energy after the proposed move*/
            bv_cost2 = 0;
            bvv_cost2 = 0;
            constraint2 = 0;
            for (auto &calc : this->calculators)
            {
                /*Compare calculator site to modified site*/
                if (calc->site == st)
                {
                    constraint2 += calc->execute(s, i, j, k);
                }
            }
            if (optimize_valence == true)
            {
                bv_cost2 = s.local_instability_index(i, j, k, st);
            }
            if (optimize_vector_valence == true)
            {
                bvv_cost2 = s.local_valence_vector_change(i, j, k, st);
            }
            /*Calculate the change*/

            change = (bv_cost2 - bv_cost1) / T1 + (bvv_cost2 - bvv_cost1) / T2 + (constraint2 - constraint1) / T3;
            /*Use Metropolis/Glauber testing to decide whether to accept the
             * proposed move*/
            if (change > 0)
            {
                boltzmann = exp(-change);
                p = dis(gen);

                if (p < boltzmann)
                {
                    bad++;
                }
                else
                {
                    mod->undo_last(s);
                }
            }
            else
            {
                if (change == 0)
                {
                    neutral++;
                }
                else
                {
                    good++;
                }
            }
        }
    }

    py::print("Good: ", good, "    Neutral: ", neutral, "    Bad: ", bad, " out of ",
              this->generators.size() * s.shape * s.shape * s.shape);

    /*Reseting rng's*/
    dis.reset();
    pick.reset();
    for (auto mod : this->generators)
    {
        mod->reset_rng();
    }
}

/*------------------------------------------------------------------------------------------------------------------------*/

PYBIND11_MODULE(BvsOptimizer, handle)
{
    handle.doc() = "A simple Pybind11 wrapped C++ module for atomistic Monte Carlo simuluations"
                   "of local structure of disordered materials. It implements a "
                   "few simple interatomic interaction potentials and "
                   "provides the possibility of Monte Carlo "
                   "optimization of atomic oxidation states using bond valence sum"
                   "method from crystal chemistry, thereby providing a simple way of"
                   "simulating the local bonding geometry of atoms within lattices. ",
    /*Exposed classes and their selected methods intended to be called from
     * Python.*/
        py::class_<structure>(handle, "Crystal",
                              "A class for storage, manipulation and anaysis of a "
                              "disordered structure.")
            .def(py::init<int, int, py::array_t<double>>(), "Class constructor.", py::arg("ncells"),
                 py::arg("sites_per_cell"), py::arg("unitcell"))

            .def("load", &structure::load,
                 "Loads the atomic numbers of atoms and theirpositions within the "
                 "unit "
                 "cell (lattice units).",
                 py::arg("atoms"), py::arg("positions"))

            .def("get_atoms", &structure::get_atoms,
                 "Returns a list of atomic numbers. Automatic casting to numpy "
                 "array. "
                 "Used for conversion to javelin.Structure object.")

            .def("get_positions", &structure::get_positions,
                 "Returns the atomic positions within the cell in lattice units. "
                 "Used "
                 "for conversion to javelin.Structure object.")

            .def("get_cell", &structure::get_cell, "Returns the 2D array containing primitive lattice vectors.")

            .def("get_shape", &structure::get_shape,
                 "Returns the shape of the supercell in the form [Nx, Ny, Nz, "
                 "sites_per_cell].")

            .def("get_scaled_position", &structure::get_scaled_position)

            .def("agitate", &structure::agitate,
                 "Displaces atoms at a given site by a random vector drawn from a "
                 "normal distribution of specified st. deviation.",
                 py::arg("site"), py::arg("sigma_x"), py::arg("sigma_y"), py::arg("sigma_z"))

            .def("set_bonds", &structure::set_bonds,
                 "Specifies bonded neighbors for a given site in the form of "
                 "javelin.Neighborlist.",
                 py::arg("site"), py::arg("neighborlist"))

            .def("set_parameters", &structure::set_parameters,
                 "Sets R and b values used in empirical BV calculations for a "
                 "given pair of atomic species",
                 py::arg("atom1"), py::arg("atom2"), py::arg("R"), py::arg("b"))

            .def("set_nominals_and_weights", &structure::set_nominals_and_weights,
                 "Sets the ideal valences/vector valences and weights for a list of "
                 "atomic species",
                 py::arg("atom_list"), py::arg("ideal_valences"), py::arg("ideal_vector_valences"), py::arg("weights"),
                 py::arg("vector_weights"))

            .def("get_valence_list", &structure::valence_list,
                 "Returns a list of oxidation states of a given atomic species at a "
                 "particular site.",
                 py::arg("site"), py::arg("atom"))

            .def("get_average_valence", &structure::average_valence,
                 "Returns average oxidation state/valence of a given atom", py::arg("atom"))

            .def("global_instability_index", &structure::global_instability_index,
                 "Returns a so called global instability index which represents a "
                 "measure of"
                 "deviation of empirical valences of atoms in a particular "
                 "structure/bonding environment "
                 "from their ideal/nominal values")

            .def("get_distance_list", &structure::distance_list,
                 "Returns a list of distances between neighboring atoms of a given "
                 "species.",
                 py::arg("atom1"), py::arg("atom2"), py::arg("neighborlist"))

            .def("average_distance", &structure::average_distance,
                 "Returns an average distance between neighboring atoms of a given "
                 "kind.",
                 py::arg("atom1"), py::arg("atom2"), py::arg("neighborlist"))

            .def("rms_displacement", &structure::rms_displacement,
                 "Returns the standard deviation of atomic displacements around a "
                 "given site.",
                 py::arg("site"))

            .def("local_polarization", &structure::local_dipole_moment,
                 "Returns the avereage magnitude of local dipole moment of a given "
                 "coordination complex.",
                 py::arg("neighborlist"))

            .def("CS_parameter", &structure::CS_parameter,
                 "Returns the centrosymmetry parameter for a give site coordination.", py::arg("neighborlist"))

            .def("atomic_displacement_parameters", &structure::atomic_displacement_parameters,
                 "Returns the anisotropic atomic displacement parameters "
                 "(B-factors) "
                 "for a given site",
                 py::arg("site"));

    py::class_<Generator>(handle, "Generator",
                          "A class that performs the structure modification in "
                          "each step of MC sampling.")
        .def(py::init<int, py::array_t<double>, py::array_t<double>>(), "Class constructor.", py::arg("site"),
             py::arg("average_position"), py::arg("std. deviation"));

    py::class_<Calculator>(handle, "Calculator",
                           "A class used to perform energy calculations during the simulation. It "
                           "needs to be created for each modified site and added to the "
                           "simulation "
                           "engine. It can also be used to manage constraints which are "
                           "implemented "
                           "as pseudo-interactions.")
        .def(py::init<int, py::array_t<int>>(), "Class constructor.", py::arg("site"), py::arg("interacting_neighbors"))
        .def("add_interaction", &Calculator::add_interaction,
             "Adds a supported interaction energy for a pair of atomic species.", py::arg("atom1"), py::arg("atom2"),
             py::arg("interaction_type"), py::arg("interaction_parameters"));

    py::class_<simulation>(handle, "Simulation", "Simulation engine class.")
        .def(py::init<>(), "Constructor")
        .def_readwrite("optimize_valence", &simulation::optimize_valence)
        .def_readwrite("optimize_vector_valence", &simulation::optimize_vector_valence)
        .def("add_generator", &simulation::add_generator,
             "Add a move generator to each site intended to be modified during "
             "the simulation",
             py::arg("generator"))
        .def("add_calculator", &simulation::add_calculator, "Add a calculator for each modified site.",
             py::arg("calculator"))
        .def("run_cycle", &simulation::run_cycle,
             "Main method used to start the MC cycle. It receives (by "
             "reference) "
             "the initial structure and modifies it in place. Also, it uses 3 "
             "temperatures to independently control the acceptance ratio for "
             "valence optimization, valence vector optimization and constraints "
             "(interactions)",
             py::arg("structure"), py::arg("BV_temperature"), py::arg("BVV_temperature"),
             py::arg("constraint_temperature"));
}
