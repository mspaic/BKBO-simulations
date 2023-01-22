
import sys
from javelin.structure import Structure
import numpy as np
from os import remove
from os import path
import subprocess as sp
import signal


class VerboseCalledProcessError(sp.CalledProcessError):
    """
       Class used to redirect terminal error messages to
       Python console when an external program is called
       using bash function.
    """

    def __str__(self):
        if self.returncode and self.returncode < 0:
            try:
                msg = "Command '%s' died with %r." % (
                    self.cmd, signal.Signals(-self.returncode))
            except ValueError:
                msg = "Command '%s' died with unknown signal %d." % (
                    self.cmd, -self.returncode)
        else:
            msg = "Command '%s' returned non-zero exit status %d." % (
                self.cmd, self.returncode)

        return f'{msg}\n' \
               f'Stdout:\n' \
               f'{self.output}\n' \
               f'Stderr:\n' \
               f'{self.stderr}'


def bash(cmd, print_stdout=True, print_stderr=True):
    proc = sp.Popen(cmd, stderr=sp.PIPE, stdout=sp.PIPE, shell=True, universal_newlines=True,
                    executable='/bin/bash')
    """
        Function used for running external programs and 
        and redirecting their console output.
    """
    all_stdout = []
    all_stderr = []
    while proc.poll() is None:
        for stdout_line in proc.stdout:
            if stdout_line != '':
                if print_stdout:
                    print(stdout_line, end='')
                all_stdout.append(stdout_line)
        for stderr_line in proc.stderr:
            if stderr_line != '':
                if print_stderr:
                    print(stderr_line, end='', file=sys.stderr)
                all_stderr.append(stderr_line)

    stdout_text = ''.join(all_stdout)
    stderr_text = ''.join(all_stderr)
    if proc.wait() != 0:
        raise VerboseCalledProcessError(proc.returncode, cmd, stdout_text, stderr_text)


def clean(directory, name, n_samples, extension):
    """
    Deletes a batch of files whose path is of the form
    "directory/name__.extension" where __ is a placeholder
    for a number in range from 0 to n_samples.
    """
    for i in range(n_samples):
        file = directory+name+str(i)+extension
        if (path.isfile(file)):
            remove(file)


def read_structure(path_to_file=""):
    """Reads javelin.structure object from a .txt file """
    sep = " "
    if (path.isfile(path_to_file)):
        g = open(path_to_file, mode="r")
    else:
        print("File not found!")
        return 0
    cell = g.readline()
    supercell = []
    for i in cell.split(sep):
        if i != "" and i != "SUPERCELL:" and i != "\n":
            supercell.append(i)
    supercell = np.array(supercell, dtype=int)
    cell = g.readline()
    unitcell = []
    for i in cell.split(sep):
        if i != "" and i != "CELL:" and i != "\n":
            unitcell.append(i)
    unitcell = np.array(unitcell, dtype=float)
    g.readline()
    lines = g.read().splitlines()
    sym = []
    pos = []
    temp = []
    for line in lines:
        for i in line.split(sep):
            if i != "":
                temp.append(i)
        pos.append((temp[0], temp[1], temp[2]))
        sym.append(temp[3])
        temp = []
    g.close()
    structure = Structure(symbols=sym, positions=pos, unitcell=unitcell)
    structure.reindex(supercell)
    return structure


def write_structure(structure=None, path_to_file="structure.txt"):
    """Writes javelin.structure object to a .txt file """

    sep = "    "
    unit = structure.unitcell
    shape = structure.atoms.index.to_numpy().max() + np.array(1)
    f = open(path_to_file, "w")
    f.write("SUPERCELL: %4d %4d %4d %4d\n" % (shape[0], shape[1], shape[2], shape[3]))
    f.write("CELL: %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n" % unit.cell)
    f.write("%15s %15s %15s %s\n" % ('X [l.u.]', 'Y [l.u.].', 'Z [l.u.]', 'EL'))
    f.close()
    df = structure.atoms.copy()
    cols = ['x', 'y', 'z', 'symbol']
    df = df[cols]
    f = open(path_to_file, mode="a")
    np.savetxt(f, df.values, delimiter=sep, fmt="%15.10f %15.10f %15.10f %s")
    f.close()


class Scatty:
    """
    Class used for setting up a scatty calculation and dealing with scatty input and output files.
    Attributes mostly correspond to avaliable scatty config options (with few additions).
    """

    def __init__(self):
        self.name = "HKL"
        """
        A name for the calculation (e.g., the scattering plane to be calculated). This will be included in the filename of
        the output files.
        """
        self.title = "Structure"
        """
        A name given to input structure files ('title_atoms_0,1,2...txt')
        """

        self.path = ""
        """
        A path to desired input/output directory.
        """

        self.replacement_dict = {}
        """
        A dictionary to replace one atomic species (key) with another (value) when writing input files. 
        """

        self.basis = []
        """
        A list of average site positions within the unit cell (fractional coordinates). 
        """

        self.centre = (0.0, 0.0, 0.0)
        """
        The centre of the calculated scattering pattern in reciprocal space, in reciprocal-lattice (hkl) units.
        The format is three real numbers, giving the components of the origin shift from Q = 0. The default value is 0 0 0.
        """

        self.x_axis = (1, 0, 0, 100)
        """
        The reciprocal-space vector used as the x-axis of the calculated scattering pattern. The format is 3 real
        numbers followed by 1 integer. The three real numbers give the components of the reciprocal-space vector relative
        to the origin, in reciprocal-lattice (hkl) units. The integer p specifies the number of points along the x direction.
        Scatty will calculate 2p + 1 points along this direction, extending from ORIGIN – X_AXIS to ORIGIN + X_AXIS
        """
        self.y_axis = (0, 1, 0, 100)
        """
        The reciprocal-space vector used as the y-axis of the calculated scattering pattern (format as above), and the
        corresponding number of points.
        """

        self.z_axis = (0, 0, 0, 1)
        """
        The reciprocal-space vector used as the z-axis of the calculated scattering pattern (format as above), and the
        corresponding number of points.
        """

        self.radiation = "X"
        """
        The type of incident radiation for the scattering simulation, indicated by a single character: either X for
        X-ray scattering, N for neutron scattering, or E for electron scattering.
        """

        self.expansion_error = None
        """
        One real number x, giving the maximum acceptable error in the Taylor
        expansion of exp(iG · u), where G is a wavevector and u is an atomic displacement (see “Ultrafast calculation of
        diffuse scattering from atomistic models”). This should be a number much smaller than unity. If the atoms files
        include displacive disorder, EXPANSION_MAX_ERROR and/or EXPANSION_ORDER must be given. The way in which the
        maximum-error constraint |eps| ≤ x is implemented depends on whether EXPANSION_ORDER n is also given:
        • If EXPANSION_ORDER n is not given, then the order n of the Taylor expansion is determined at each |G| by the
          requirement that |eps| ≤ x;
        • If EXPANSION_ORDER n is given, then where the n-th order Taylor expansion of exp(iG · u) could yield an error
          |eps| > x, Scatty will replace the Taylor expansion with direct calculation of exp(iG · u), to ensure that |eps| ≤ x.
        In most cases, the calculation will run fastest if you specify EXPANSION_MAX_ERROR and do not specify EXPANSION_ORDER.
        The possible exception is for supercells that contain large displacements, because a very large number of terms in the
        Taylor expansion can then be required to obtain reasonable accuracy (e.g., for displacements of ~0.5 Å and typical
        |G max | ~ 15 Å -1 , a maximum error |eps| < 0.05 would require an expansion to order 20). In such cases, it may be preferable
        to specify EXPANSION_ORDER as well as EXPANSION_MAX_ERROR.
        """

        self.expansion_order = None
        """
        One integer n, specifying the number of terms included in the Taylor expan-
        sion of exp(iG · u), where G is a wavevector and u is an atomic displacement (see “Ultrafast calculation of dif-
        fuse scattering from atomistic models”). If the atoms files include displacive disorder, EXPANSION_ORDER and/or
        EXPANSION_MAX_ERROR must be given. If EXPANSION_MAX_ERROR is given, n is the maximum number of terms
        included in the Taylor expansion (see EXPANSION_MAX_ERROR above).
        """

        self.remove_bragg = False
        """
        If this keyword is given, Scatty will remove nuclear Bragg peaks from the scattering
        calculation. If given, it should be followed by one of the letters P, I, F, R, A, B, C, or H to specify the centring of the
        crystallographic unit cell. This keyword has no effect for magnetic scattering.
        """

        self.cell_centering = None
        """
        Type of cell centering. Need to specify if Bragg peaks are removed. 
        """

        self.window_type = "lanczos"
        """
        Type of window to be used for resampling. Allowed options are 'lanczos' and 'cosine'. 
        """

        self.window = 3
        """
        The value of m in the Lanczos interpolation formula, which is an integer ≥ 2 (see “Ultrafast
        calculation of diffuse scattering from atomistic models”). The default value is 3. If WINDOW 0 is given, nearest-
        neighbour interpolation is used instead of Lanczos resampling.
        """
        self.cutoff = 2
        """
        The value of m' in the Lanczos interpolation formula (see “Ultrafast calculation of diffuse scattering
        from atomistic models”). The default value is 2.
        """
        self.sum_type = "PARALLEL"
        """
        A character string giving the summation type in real space. The options are PARALLEL (sum over atom
        pairs within a parallelepiped) or SPHERE (sum over atom pairs within a sphere). The default value is PARALLEL.
        """
        self.symmetry = None
        """
        If this keyword is given, Scatty will symmetrise the calculated scattering patterns. For
        systems showing only short-range correlations, the symmetry of the scattering pattern should be the same as the
        Laue symmetry of the crystal. Hence, the options that can follow the SYMMETRY keyword are the 11 Laue classes:
        • m-3m
        • m-3
        • 6|mmm
        • 6|m
        • -3m (hexagonal axes are assumed)
        • -3 (hexagonal axes are assumed)
        • 4|mmm
        • 4|m
        • mmm
        • 2|m
        • -1
        Due to the finite size of the spin configurations, atomistic simulations do not reproduce the Laue symmetry exactly, so it
        is often possible to obtain a large improvement in apparent statistics by symmetrising. Of course, one should check before
        using this option that the scattering pattern actually possesses the expected symmetry!
        """

        self.magnetic_only = False
        """
        If this keyword is given, nuclear scattering will be excluded from the calculation, and only
        magnetic scattering will be calculated. This is useful to compare calculations with data measured using xyz neutron
        polarisation analysis, which can allow magnetic and nuclear diffuse scattering to be separated.
        """

        self.temp_subtract = False
        """
        If this keyword is given, a ideal paramagnetic background is subtracted from the
        calculated magnetic scattering intensity. This keyword should be given if the calculated magnetic scattering pattern
        is being compared to experimental magnetic diffuse-scattering data obtained by low–high temperature subtraction.
        Please note that the magnetic scattering pattern will contain both negative and positive intensities in this case.
        """

        self.num_threads = 0
        """
        Maximum number of threads to be used by the calculation. If set to zero then number of threads is automatically 
        assigned by the OS. 
        """

        self.ppm_output = False
        """
        Writes an image file in .ppm format, a simple image format that stores pixel information
        in plain text. It is possible to convert .ppm images into standard formats (e.g., .gif or .png) using an image-editing
        program such as ImageMagick (www.imagemagick.org). This option is only available for two-dimensional plots.
        """

        self.ppm_range = ()
        """
        If the option to write a .ppm image is selected, this keyword may be given to specify the
        minimum and maximum intensity values shown on the image (two real numbers, minimum followed by maximum).
        By default, the minimum and maximum intensities in the calculated plane are used.
        """

        self.ppm_colourmap = ""
        """
        (optional) If the option to write a .ppm image is selected, this keyword may be given to specify
        the colourmap. The options are:
            • default: By default, Scatty uses the “cool-to-warm” colourmap used in the ParaView program. This map has
            been designed for scientific visualisation in order to be easily interpreted and aesthetically pleasing [3].
            • heat: Black-red-yellow-white colourmap.
            • jet: The “rainbow” colourmap used in Matlab.
            • grey1: Greyscale with a black background.
            • grey2: Greyscale with a white background.
            • viridis
        
        """

        self.sc_bragg_output = False
        """
        Output calculated intensity as a .vtk file which can be read in free software for three-dimensional data visualisation, such as
        ParaView and MayaVi. Axes are labelled in reciprocal-lattice
        units.
        """

    def write_config(self):
        """
        Generate "scatty_config.txt" file based on the given settings.
        """
        if self.remove_bragg == True and self.cell_centering == None:
            raise ValueError(
                "Cell centering information required for removal of Bragg peaks!")
        g = open(self.path+"scatty_config.txt", 'w')
        g.write("NAME  " + self.name + "\n")
        g.write("RADIATION  " + self.radiation+"\n")

        g.write("CENTRE %7.4f %7.4f %7.4f\n" % self.centre)
        g.write("X_AXIS %7.4f %7.4f %7.4f %6d\n" % self.x_axis)
        g.write("Y_AXIS %7.4f %7.4f %7.4f %6d\n" % self.y_axis)
        g.write("Z_AXIS %7.4f %7.4f %7.4f %6d\n" % self.z_axis)
        g.write("WINDOW_TYPE " + self.window_type + "\n")
        g.write("WINDOW %4d\n" % self.window)
        g.write("CUTOFF %4d\n" % self.cutoff)
        g.write("SUM_TYPE " + self.sum_type + "\n")
        if self.symmetry != None:
            g.write("SYMMETRY " + self.symmetry+"\n")
        if self.magnetic_only == True:
            g.write("MAG_ONLY\n")
        if self.temp_subtract == True:
            g.write("TEMP_SUBTRACT\n")
        if self.expansion_error != None:
            g.write("EXPANSION_MAX_ERROR %8.6f\n" % self.expansion_error)
        if self.expansion_order != None:
            g.write("EXPANSION_ORDER %4d\n" % self.expansion_order)
        if self.remove_bragg == True:
            g.write("REMOVE_BRAGG  " + self.cell_centering + "\n")
        if self.num_threads != 0:
            g.write("MAX_NUM_THREADS %4d\n" % self.num_threads)
        if self.ppm_output:
            g.write("PPM_OUTPUT")
            if self.ppm_colourmap != "":
                g.write("PPM_COLOURMAP " + str(self.ppm_colourmap)+"\n")
            if len(self.ppm_range) != 0:
                g.write("PPM_RANGE %8.4f %8.4f\n" % self.ppm_range)
        if self.sc_bragg_output:
            g.write("SUPERCELL_BRAGG_OUTPUT\n")
        g.close()

    def write_atom_files(self, structures=None):
        """
        Convert a list of javelin.structure objects to input .txt files for scatty.
        """

        if len(structures) == 0:
            raise ValueError("List of structures is empty!")
        count = 0
        sep = " "
        for stru in structures:
            unit = stru.unitcell
            count = count + 1
            shape = stru.atoms.index.to_numpy().max() + np.array(1)
            if shape[3] != len(self.basis):
                raise ValueError("Number of sites does not match the given structure!")
            if count < 10:
                filename = self.title+"_atoms_0"+str(count)+".txt"
            else:
                filename = self.title+"_atoms_"+str(count)+".txt"

            f = open(self.path+filename, "w")
            f.write("TITLE"+sep + self.title+"\n")
            f.write("CELL %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f\n" % unit.cell)
            for i in range(len(self.basis)):
                f.write("SITE %8.5f %8.5f %8.5f\n" % self.basis[i])
            f.write("BOX %8d %8d %8d\n" % (shape[0], shape[1], shape[2]))
            for i in range(len(self.basis)):
                f.write("OCC ")
                for key in stru.get_average_structure()[i]:
                    if key in self.replacement_dict:
                        f.write("%6s %8.5f" % (
                            self.replacement_dict[key], stru.get_average_structure()[i][key]["occ"]))
                    else:
                        f.write("%6s %8.5f" % (key, stru.get_average_structure()[i][key]["occ"]))
                f.write("\n")
            f.close()

            df = stru.atoms.copy()
            df.pop("Z")
            df["designator"] = ["ATOM"]*df.shape[0]
            df = df.reset_index()
            cols = ['designator', 'site', 'i', 'j', 'k', 'x', 'y', 'z', 'symbol']
            df = df[cols]
            df.loc[df['designator'] == "ATOM", ['x', 'y', 'z']
                   ] -= np.array(self.basis*int(df.shape[0]/len(self.basis)))
            df.loc[df['designator'] == "ATOM", ['site']] += 1
            for a in stru.get_atom_symbols():
                if a in self.replacement_dict:
                    df.loc[df['symbol'] == a, ['symbol']] = self.replacement_dict[a]

            f = open(self.path+filename, "a")
            np.savetxt(f, df.values, delimiter=sep,
                       fmt="%4s %10d %10d %10d %10d %15.10f %15.10f %15.10f %8s")
            f.close()

    def write_input(self, structures=[]):
        """
        Write input atom files and scatty_config file.
        """
        print("Writing input files to directory "+self.path+"...")
        if (len(structures) != 0):
            self.write_config()
            self.write_atom_files(structures)

    def run(self):
        sys.path.append("scatty/")
        from Scatty import runner
        runner(self.name)

    def read_numpy(self, filename=None):
        """
        Read calculated intensity and returns 4 numpy arrays containing h,k,l axes and intensity.
        """
        filename = self.path + self.title+"_"+self.name+"_sc_list.txt"

        if (path.isfile(filename)):
            a, b, c, d = np.loadtxt(
                r""+filename, usecols=(0, 1, 2, 3), unpack=True)
            H = np.unique(a)
            K = np.unique(b)
            L = np.unique(c)
            Intensity = d.reshape(H.size, K.size, L.size)
            return H, K, L, Intensity
        else:
            print("Intensity file not found!")
            return 0

    def read_xarray(self, axes_names=["H", "K", "L"]):
        """
        Read calculated intensity and returns xarray.Datarray object that packages both h,k,l axes and intensity into
        a single object.
        """
        from xarray import DataArray

        filename = self.path + self.title+"_"+self.name+"_sc_list.txt"

        if (path.isfile(filename)):
            a, b, c, d = np.loadtxt(
                r""+filename, usecols=(0, 1, 2, 3), unpack=True)
            h = np.unique(a)
            k = np.unique(b)
            l = np.unique(c)
            scatt = d.reshape(h.size, k.size, l.size)
            I = DataArray(scatt, dims=[axes_names[0], axes_names[1], axes_names[2]],
                          coords={axes_names[0]: h, axes_names[1]: k, axes_names[2]: l})
            return I
        else:
            print("Intensity file not found!")
            return 0

    def clean_input(self):
        """
        Delete input files.
        """
        if path.isfile(self.path+"scatty_config.txt"):
            remove(self.path+"scatty_config.txt")
        for count in range(300):
            if count < 10:
                filename = self.title+"_atoms_0"+str(count)+".txt"
            else:
                filename = self.title+"_atoms_"+str(count)+".txt"
            if path.isfile(self.path+filename):
                remove(self.path+filename)

    def clean_output(self):
        """
        Delete scatty output files.
        """
        if path.isfile(self.path+self.title+"_"+self.name+"_sc_list.txt"):
            remove(self.path+self.title+"_"+self.name+"_sc_list.txt")
        if path.isfile(self.path+self.title+"_"+self.name+"_sc.txt"):
            remove(self.path+self.title+"_"+self.name+"_sc.txt")
        if path.isfile(self.path+self.title+"_"+self.name+"_sc.vtk"):
            remove(self.path+self.title+"_"+self.name+"_sc.vtk")
        if path.isfile(self.path+self.title+"_"+self.name+"_sc.ppm"):
            remove(self.path+self.title+"_"+self.name+"_sc.ppm")
        if path.isfile(self.path+self.title+"_"+self.name+"_sc_ppm_scale.txt"):
            remove(self.path+self.title+"_"+self.name+"_sc_ppm_scale.txt")
        if path.isfile(self.path+self.title+"_"+self.name+"_sc_colourbar.ppm"):
            remove(self.path+self.title+"_"+self.name+"_sc_colourbar.ppm")
        if path.isfile(self.path+self.title+"_"+self.name+"_scatty_info.txt"):
            remove(self.path+self.title+"_"+self.name+"_scatty_info.txt")
