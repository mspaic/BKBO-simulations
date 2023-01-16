cimport cython
import numpy as np 
cimport numpy as cnp 
from libc.math cimport exp, sqrt, pow, INFINITY, NAN
from javelin.energies cimport Energy



cdef class HarmonicInteractions(Energy):

	cdef readonly double [:] lattice_scale
	cdef readonly double [:,:] average_position
	cdef readonly double [:,:,:,:,:,:,:] FC_tensor

	cdef readonly Py_ssize_t size_x, size_y, size_z


	def __init__(self, double [:,:,:,:,:,:,:] fc,
				 double [:,:] average_pos,
				 double [:] scale,
				 Py_ssize_t [:] size ):

		self.FC_tensor = fc
		self.average_position = average_pos
		self.lattice_scale = scale
		self.size_x = size[0]
		self.size_y = size[1]
		self.size_z = size[2]




	cpdef double evaluate(self,
						  int a1, double x1, double y1, double z1,
						  int a2, double x2, double y2, double z2,
						  Py_ssize_t site1, Py_ssize_t site2,
						  Py_ssize_t target_x, Py_ssize_t target_y, Py_ssize_t target_z) except *:

		cdef Py_ssize_t i,j
		cdef double u_1[3], u_2[3]
		cdef double result
		result = 0

		u_1[:] = x1, y1, z1
		u_2[:] = x2, y2, z2

		for i in range(3):
			u_1[i] = (u_1[i]  - self.average_position[site1][i])*self.lattice_scale[i]
			u_2[i] = (u_2[i]  - self.average_position[site2][i])*self.lattice_scale[i]


		if (target_x > (self.size_x - 1)) or (target_x <  (-self.size_x + 1)) :
			return 0
		if (target_y > (self.size_y - 1)) or (target_y <  (-self.size_y + 1)) :
			return 0
		if (target_z > (self.size_z - 1)) or (target_z <  (-self.size_z + 1)) :
			return 0

		#calculate quadratic form \sum_ij FC(ij) u1(i) u2(j)
		for i in range(3):
			for j in range(3):
				result += (self.FC_tensor[site1][site2][target_x][target_y][target_z][i][j]*u_1[i]*u_2[j])/2
		return result



cdef class YukawaEnergy(Energy):

	cdef readonly double length, g
	cdef readonly int atom_type1, atom_type2
	cdef unsigned char check_atom_types

	def __init__(self, double length, double g, int atom_type1=-1, int atom_type2=-1):
		self.length = length
		self.g = g
		self.atom_type1 = atom_type1
		self.atom_type2 = atom_type2
		self.check_atom_types = self.atom_type1 != -1 or self.atom_type2 != -1
	def __str__(self):
		return "{}(range={},g={},atoms={})".format(self.__class__.__name__,
													 self.length,
													 self.g,
													 '{}-{}'.format(self.atom_type1,self.atom_type2)
													 if self.check_atom_types else 'all')
	cdef unsigned char valid_atoms(self, int atom1, int atom2):
		if atom1 == self.atom_type1 and atom2 == self.atom_type2:
			return True
		elif atom1 == self.atom_type2 and atom2 == self.atom_type1:
			return True
		else:
			return False
	cpdef double evaluate(self,
						  int a1, double x1, double y1, double z1,
						  int a2, double x2, double y2, double z2,
						  Py_ssize_t site1, Py_ssize_t site2,
						  Py_ssize_t target_x, Py_ssize_t target_y, Py_ssize_t target_z) except *:
		if self.check_atom_types:
			if not self.valid_atoms(a1, a2):
				return 0
		cdef double diff = distance(x1, y1, z1, x2+target_x, y2+target_y, z2+target_z)
		return self.g*exp(-diff/self.length)/diff


cdef class MorseEnergy(Energy):

	cdef readonly double D,b,desired

	cdef readonly int atom_type1, atom_type2
	cdef unsigned char check_atom_types
	def __init__(self, double D, double b, double desired, int atom_type1=-1, int atom_type2=-1):
		self.D = D
		self.b = b
		self.desired = desired
		self.atom_type1 = atom_type1
		self.atom_type2 = atom_type2
		self.check_atom_types = self.atom_type1 != -1 or self.atom_type2 != -1
	def __str__(self):
		return "{}(D={},b={}, desired ={}, atoms={})".format(self.__class__.__name__,
													 self.D,
													 self.b,
													 self.desired,
													 '{}-{}'.format(self.atom_type1,self.atom_type2)
													 if self.check_atom_types else 'all')
	cdef unsigned char valid_atoms(self, int atom1, int atom2):
		if atom1 == self.atom_type1 and atom2 == self.atom_type2:
			return True
		elif atom1 == self.atom_type2 and atom2 == self.atom_type1:
			return True
		else:
			return False
	cpdef double evaluate(self,
						  int a1, double x1, double y1, double z1,
						  int a2, double x2, double y2, double z2,
						  Py_ssize_t site1, Py_ssize_t site2,
						  Py_ssize_t target_x, Py_ssize_t target_y, Py_ssize_t target_z) except *:
		if self.check_atom_types:
			if not self.valid_atoms(a1, a2):
				return 0
		cdef double diff = distance(x1, y1, z1, x2+target_x, y2+target_y, z2+target_z) - self.desired
		return self.D*pow((1-exp(-diff/self.b)),2)





cdef class VolumeExclusionEnergy(Energy):

	cdef readonly r1, r2

	cdef readonly int atom_type1, atom_type2
	cdef unsigned char check_atom_types
	def __init__(self, double r1, double r2, int atom_type1=-1, int atom_type2=-1):
		self.r1 = r1
		self.r2 = r2
		self.atom_type1 = atom_type1
		self.atom_type2 = atom_type2
		self.check_atom_types = self.atom_type1 != -1 or self.atom_type2 != -1
	def __str__(self):
		return "{}(r1={},r2={}, atoms={})".format(self.__class__.__name__,
													 self.r1,
													 self.r2,
													 '{}-{}'.format(self.atom_type1,self.atom_type2)
													 if self.check_atom_types else 'all')
	cdef unsigned char valid_atoms(self, int atom1, int atom2):
		if atom1 == self.atom_type1 and atom2 == self.atom_type2:
			return True
		elif atom1 == self.atom_type2 and atom2 == self.atom_type1:
			return True
		else:
			return False
	cpdef double evaluate(self,
						  int a1, double x1, double y1, double z1,
						  int a2, double x2, double y2, double z2,
						  Py_ssize_t site1, Py_ssize_t site2,
						  Py_ssize_t target_x, Py_ssize_t target_y, Py_ssize_t target_z) except *:
		if self.check_atom_types:
			if not self.valid_atoms(a1, a2):
				return 0
		cdef double diff = distance(x1, y1, z1, x2+target_x, y2+target_y, z2+target_z)
		if diff>(self.r1 + self.r2):
			return 0
		else:
			return 1e10


cdef class BornMayerRepulsiveEnergy(Energy):

	cdef readonly A,l

	cdef readonly int atom_type1, atom_type2
	cdef unsigned char check_atom_types
	def __init__(self, double A, double l, int atom_type1=-1, int atom_type2=-1):
		self.A = A
		self.l = l
		self.atom_type1 = atom_type1
		self.atom_type2 = atom_type2
		self.check_atom_types = self.atom_type1 != -1 or self.atom_type2 != -1
	def __str__(self):
		return "{}(A={},l={},atoms={})".format(self.__class__.__name__,
													 self.A,
													 self.l,
													 '{}-{}'.format(self.atom_type1,self.atom_type2)
													 if self.check_atom_types else 'all')
	cdef unsigned char valid_atoms(self, int atom1, int atom2):
		if atom1 == self.atom_type1 and atom2 == self.atom_type2:
			return True
		elif atom1 == self.atom_type2 and atom2 == self.atom_type1:
			return True
		else:
			return False
	cpdef double evaluate(self,
						  int a1, double x1, double y1, double z1,
						  int a2, double x2, double y2, double z2,
						  Py_ssize_t site1, Py_ssize_t site2,
						  Py_ssize_t target_x, Py_ssize_t target_y, Py_ssize_t target_z) except *:
		if self.check_atom_types:
			if not self.valid_atoms(a1, a2):
				return 0
		cdef double diff = distance(x1, y1, z1, x2+target_x, y2+target_y, z2+target_z)
		return self.A*exp(-diff/self.l)



#calculates distance between 2 points
cdef double distance(double x1, double y1, double z1,
					 double x2, double y2, double z2):
	cdef double dX = x2 - x1
	cdef double dY = y2 - y1
	cdef double dZ = z2 - z1
	return sqrt( dX*dX + dY*dY + dZ*dZ )

"""
cdef double distance_sq(double x1, double y1, double z1,
					 double x2, double y2, double z2):
	cdef double dX = x2 - x1
	cdef double dY = y2 - y1
	cdef double dZ = z2 - z1
	return dX*dX + dY*dY + dZ*dZ

"""

