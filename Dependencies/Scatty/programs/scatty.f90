PROGRAM scatty

   USE singleton_dble
   USE omp_lib

   IMPLICIT NONE

   ! Calculates diffuse I(Q)

   INTEGER :: i, j, a, e, m, n, p, num_cells(3), x, y, z, multiplier, sum_type, iterm
   INTEGER :: ibox, num_boxes, num_atoms_in_cell, i0(3), i1(3), i_max(3), i_min(3), isymm
   INTEGER :: iq_bin, num_q_bins, itaylor, jtaylor, max_taylor, n_taylor, max_num_terms, from, to
   INTEGER :: laue_class, supercell(4), colourmap, nq(3), qx, qy, qz, xtal_radiation, hkl_fft(3), max_num_elements
   INTEGER :: h, k, l, max_bragg(3), min_bragg(3), window_type, num_symm_ops(12), ihkl(3), tot_num_cells
   INTEGER :: count_taylor(3), centring, num_elements_in_cell, max_num_threads, it_count, progress, num_it
   INTEGER, PARAMETER :: nq_bin = 5
   INTEGER(8) :: num_hkl, num_q, iq, idft
   INTEGER, ALLOCATABLE :: c_r(:, :, :, :, :), num_elements(:), nests(:), prod(:, :), num_products(:), ielement_in_cell(:)
   INTEGER, ALLOCATABLE :: item_index(:), count_rms(:), itaylor_from_iterm(:), n_taylor_bin(:), num_terms(:)
   INTEGER(8), ALLOCATABLE :: max_num_dft(:, :), num_hkl_bin(:), num_dft(:, :, :)
   REAL(8), PARAMETER :: eps = 10D0*EPSILON(0.0), mag_prefactor = 0.07265289D0, mem_warning = 1D0
   REAL(8) :: s_i(3), hkl(3), r_max, q_sq, iq_sq, dist_from_bragg(3), tpi, norm, num_atoms, max_err, expansion_max_error
   REAL(8) :: tot_num_exp, cell_params(3), pi, r_cut, ss, cell_angles(3), temp_subtract_term, u_max, u_mod, u_rms_max
   REAL(8) :: q(3), r_orth(3, 3), q_orth(3, 3), r_metric(3, 3), q_metric(3, 3), vol, mr_cut, sphere_norm, q_mod, q_max
   REAL(8) :: axis(3, 3), axis_orth(3, 3), axis_length(3), origin(3), origin_orth(3), centre(3), centre_orth(3), spacing_orth(3)
   REAL(8) :: weight, qr, dhkl(3), mem_unit, mem_footprint, hkl_scale(3), ex, min_u_cutoff, clip_dist, q_sq_to_ss, separation(3, 3)
   REAL(8) :: scale, flat_bgr, linear_bgr, chi_sq, hkl2xyz(3, 3), axis2(3, 3), zero, t0, t1, t2, xtal_i_sq_sum
   REAL :: max, min
   REAL(8), ALLOCATABLE :: intensity_bragg(:, :, :), intensity_laue(:, :, :), qvec_from_iq(:, :)
   REAL(8), ALLOCATABLE :: u_r(:, :, :, :, :, :), spins(:, :, :, :, :, :), a_r(:, :, :), intensity(:)
   REAL(8), ALLOCATABLE :: xff(:, :), nsl(:), eff(:, :), mag_ff_hkl(:), ff0(:, :), ff2(:, :), q_mod_all(:), u_iso(:)
   REAL(8), ALLOCATABLE :: u_product(:, :, :), u_product_sgl(:, :, :), c2(:), u_rms(:), conc(:), conc_recalc(:)
   REAL(8), ALLOCATABLE :: cell_coords(:, :), box(:, :), window(:, :), hkl_product(:), u_r_scale(:, :, :, :), tot_moment_sq(:)
   REAL(8), ALLOCATABLE :: r_dft_list(:, :, :, :, :), u_r_dft_list(:, :, :, :, :), b_r_dft_list(:, :, :, :), b_r(:, :, :)
   REAL(8), ALLOCATABLE :: u_min_bin(:), q_max_bin(:), q_dot_u(:), intensity_expt(:), errors_sq(:)
   COMPLEX(8), ALLOCATABLE :: au_q(:, :, :, :, :, :), a_q(:, :, :, :, :), s_q(:, :, :, :, :, :)
   COMPLEX(8) :: mag_amplitude(3), ft, amplitude, bragg_sf0, bragg_sf, sf2
   COMPLEX(8), ALLOCATABLE :: pre_taylor(:), pre_taylor_dft(:), sf(:), struct_fac(:), mag_struct_fac(:)
   COMPLEX(8), ALLOCATABLE :: u_hkl_product(:), u_hkl_product_bragg(:), xff_hkl(:)
   CHARACTER(128) :: title, laue_class_str, colourmap_str, filename_stem, out_str, name, sum_str, window_str, elem
   CHARACTER(128) :: read_file, atom_file, spin_file, data_file
   CHARACTER(85) :: dash
   CHARACTER(999) :: str
   CHARACTER(1) :: radiation_str, centring_str, dummy
   CHARACTER(4) :: itc
   CHARACTER(4), ALLOCATABLE :: element_str(:, :)
   LOGICAL :: output_vtk, output_ppm, output_ascii, output_supercell_bragg, temp_subtract, mag_only, fit_data
   LOGICAL :: file_exists, spin_file_exists, atom_file_exists, use_expansion_max_error, found_n_taylor, clip_displacements
   LOGICAL :: use_disp, use_occ, use_mag, bragg, remove_bragg, gamma_point, found_element, use_direct_ft, use_expansion_order
   LOGICAL :: refine_flat_bgr, refine_linear_bgr, refine_scale, suppress_errors
   LOGICAL, ALLOCATABLE :: keep_hkl(:, :, :), keep_hkl_symm(:, :, :), valid(:)

   pi = ACOS(-1D0)
   tpi = 2D0*pi
   q_sq_to_ss = 1D0/(16D0*pi*pi)
   dash = REPEAT("=", 80)
   zero = 0D0
   num_symm_ops = [2, 4, 8, 8, 16, 6, 12, 12, 12, 24, 24, 48]
   !Read input from the command line
   i = iargc()
   IF (i < 1) THEN
      WRITE (*, *) 'Please start the program by writing ./scatty [input file name stem]'
      WRITE (*, *) ''
      STOP
   ELSE
      CALL getarg(1, filename_stem)
   END IF

   WRITE (*, 56) dash
   WRITE (*, *) 'Welcome to SCATTY!'
   WRITE (*, *) 'Diffuse scattering calculator'
   WRITE (*, *) '(version dated December 2022)'
   WRITE (*, 56) dash

   CALL read_spindiff_config

   ! read title, lattice parameters, number of cell_coords in unit cell and their coordinates,
   ! size of supercell from lowest-numbered config.
   num_boxes = 0
   DO k = 1, 9999
      atom_file = TRIM(filename_stem)//'_atoms_'//TRIM(itc(k))//'.txt'
      spin_file = TRIM(filename_stem)//'_spins_'//TRIM(itc(k))//'.txt'
      INQUIRE (file=spin_file, exist=spin_file_exists)
      INQUIRE (file=atom_file, exist=atom_file_exists)
      IF (spin_file_exists .AND. (.NOT. atom_file_exists)) THEN
         read_file = spin_file
         file_exists = .TRUE.
      ELSEIF (atom_file_exists .AND. (.NOT. spin_file_exists)) THEN
         read_file = atom_file
         file_exists = .TRUE.
      ELSEIF (atom_file_exists .AND. spin_file_exists) THEN
         WRITE (*, *) 'Input error: ATOMS and SPINS files cannot both be present.'
         STOP
      ELSE
         file_exists = .FALSE.
      END IF
      IF (file_exists) THEN
         num_boxes = num_boxes + 1
         IF (num_boxes == 1) THEN
            WRITE (*, *) 'Reading average-structure properties from lowest-numbered supercell...'
            WRITE (*, *) 'Entries read are TITLE, CELL, BOX, SITE, OCC, C2,'
            WRITE (*, *) 'FORM_FACTOR_J0, and FORM_FACTOR_J2.'
            CALL read_config_short
         END IF
      END IF
   END DO
   IF (num_boxes == 0) THEN
      WRITE (*, *) 'Error: No supercells found.'
      STOP
   END IF
   WRITE (*, '(i3,a)') num_boxes, ' supercells found.'

   CALL calc_metric_tensor(cell_params, cell_angles, r_orth, q_orth, r_metric, q_metric, vol)
   m = 0
   sphere_norm = 0D0
   IF (n .NE. 0) THEN
      n = NINT(n*multiplier/2D0)
      m = 2*n
      r_max = 1.0/multiplier
      r_cut = r_max - 1D0/m ! r_max=0.5 since box is normalised to side 1
      CALL calc_max_radius(r_orth, num_cells, mr_cut)
      mr_cut = mr_cut*(r_max - 1D0/m)
      sphere_norm = 4D0*pi*mr_cut**3D0/(tot_num_cells*vol)
      IF ((r_cut - 1D0/m) < -eps) STOP 'Config error: Real-space minimum cannot be less than zero.'
   END IF
   centre_orth = MATMUL(q_orth, centre(:))
   DO i = 1, 3
      axis_orth(i, :) = MATMUL(q_orth, axis(i, :))
      axis_length(i) = DOT_PRODUCT(axis_orth(i, :), axis_orth(i, :))
      IF (ABS(axis_length(i)) > eps) THEN
         axis_length(i) = SQRT(axis_length(i))
      ELSE
         axis_length(i) = 0D0
      END IF
      DO j = 1, i - 1
         IF (ABS(DOT_PRODUCT(axis_orth(i, :), axis_orth(j, :))) > eps) THEN
            WRITE (*, '(a, f8.4)') 'WARNING: Reciprocal lattice vectors do not seem to be orthogonal; dot product = '&
            &, DOT_PRODUCT(axis_orth(i, :), axis_orth(j, :))
            WRITE (*, *) ''
         END IF
      END DO
      IF (ABS(axis_length(i)) < eps) nq(i) = 0
   END DO

   data_file = 'scatty_data_01.txt'
   INQUIRE (file=data_file, exist=file_exists)
   IF (file_exists) THEN
      fit_data = .TRUE.
   ELSE
      fit_data = .FALSE.
   END IF
   IF (max_num_threads > 0) CALL OMP_SET_NUM_THREADS(max_num_threads)

   OPEN (7, file=TRIM(ADJUSTL(filename_stem))//'_'//TRIM(ADJUSTL(name))//'_scatty_info.txt', status='replace')
   DO i = 6, 7
      WRITE (i, 56) dash
      WRITE (i, *) 'Calculation name: ', TRIM(ADJUSTL(name))
      WRITE (i, '(a, i3)') ' Number of OpenMP threads: ', omp_get_max_threads()
      WRITE (i, *) 'Window type: ', TRIM(window_str)
      WRITE (i, '(a16,i3)') 'Window length: ', n
      WRITE (i, '(a16,i3)') 'Window cutoff: ', multiplier
      IF (n > 0) THEN
         WRITE (i, *) 'Sum type: ', TRIM(sum_str)
         WRITE (i, '(a29,f6.4,a5,f6.4)') ' Real-space cutoff distance: ', r_cut, ' +/- ', 1D0/m
      END IF

      WRITE (i, *) 'Laue class: ', TRIM(laue_class_str)
      WRITE (i, *) 'Radiation: ', TRIM(radiation_str)
      WRITE (i, *) 'Remove Bragg peaks? ', remove_bragg, ' centring ', centring_str
      WRITE (i, *) 'Calculate magnetic scattering only? ', mag_only
      WRITE (i, *) 'Temperature subtraction? ', temp_subtract
      IF (use_expansion_order) THEN
         WRITE (i, '(a, l2, i3)') ' Maximum order of Taylor expansion of exp(iG.u)? ', use_expansion_order, max_taylor
      ELSE
         WRITE (i, '(a,l2)') ' Maximum order of Taylor expansion of exp(iG.u)? ', use_expansion_order
      END IF
      IF (use_expansion_max_error) THEN
         WRITE (i, '(a, l2, f10.7)') ' Maximum error in calculation of exp(iG.u)? ', use_expansion_max_error, expansion_max_error
      ELSE
         WRITE (i, '(a, l2)') ' Maximum error in calculation of exp(iG.u)? ', use_expansion_max_error
      END IF
      WRITE (i, *) 'Compare with input data?', fit_data
      WRITE (i, *) 'Refine intensity scale factor?', refine_scale
      WRITE (i, *) 'Refine flat background?', refine_flat_bgr
      WRITE (i, *) 'Refine linear-in-Q background?', refine_linear_bgr
      WRITE (i, *) 'Output .ppm? ', output_ppm, TRIM(colourmap_str)
      WRITE (i, *) 'Output supercell Bragg intensities? ', output_supercell_bragg
      WRITE (i, '(a)') ' Calculating diffuse scattering in range (reciprocal-lattice units):'
      WRITE (i, 52) ' from HKL_X: ', centre(:) - axis(1, :)
      WRITE (i, 52) ' to   HKL_X: ', centre(:) + axis(1, :)
      WRITE (i, 52) ' from HKL_Y: ', centre(:) - axis(2, :)
      WRITE (i, 52) ' to   HKL_Y: ', centre(:) + axis(2, :)
      WRITE (i, 52) ' from HKL_Z: ', centre(:) - axis(3, :)
      WRITE (i, 52) ' to   HKL_Z: ', centre(:) + axis(3, :)
      IF (i == 7) THEN
         WRITE (i, '(a)') ' Calculating diffuse scattering in range (wave-vector units):'
         WRITE (i, 53) ' from Q_X: ', centre_orth(:) - axis_orth(1, :)
         WRITE (i, 53) ' to   Q_X: ', centre_orth(:) + axis_orth(1, :)
         WRITE (i, 53) ' length:   ', axis_length(1)
         WRITE (i, 53) ' from Q_Y: ', centre_orth(:) - axis_orth(2, :)
         WRITE (i, 53) ' to   Q_Y: ', centre_orth(:) + axis_orth(2, :)
         WRITE (i, 53) ' length:   ', axis_length(2)
         WRITE (i, 53) ' from Q_Z: ', centre_orth(:) - axis_orth(3, :)
         WRITE (i, 53) ' to   Q_Z: ', centre_orth(:) + axis_orth(3, :)
         WRITE (i, 53) ' length:   ', axis_length(3)
         WRITE (i, 56) dash
      END IF
52    FORMAT(A12, 3F10.5)
53    FORMAT(A10, 3F10.5)
   END DO
   WRITE (7, '(a, 3f8.4, 3f8.4)') 'Unit-cell dimensions: ', cell_params, cell_angles
   WRITE (7, '(a, 3i3)') 'Supercell dimensions: ', num_cells
   WRITE (7, 56) dash
   WRITE (7, *) 'Fractional coordinates: '
   DO a = 1, num_atoms_in_cell
      WRITE (7, *) cell_coords(:, a)
   END DO
   IF (.NOT. mag_only) THEN
      WRITE (7, 56) dash
      WRITE (7, *) 'Average occupancy: '
      DO a = 1, num_atoms_in_cell
         WRITE (7, *) (TRIM(element_str(e, a)), conc(e + ielement_in_cell(a)), e=1, num_elements(a))
      END DO

      WRITE (7, 56) dash
      WRITE (7, *) 'X-ray form factors: '
      DO a = 1, num_atoms_in_cell
         WRITE (7, *) (TRIM(element_str(e, a)), NEW_LINE('a'), xff(:, e + ielement_in_cell(a)), e=1, num_elements(a))
      END DO

      WRITE (7, 56) dash
      WRITE (7, *) 'Neutron scattering lengths: '
      DO a = 1, num_atoms_in_cell
         WRITE (7, *) (TRIM(element_str(e, a)), nsl(e + ielement_in_cell(a)), e=1, num_elements(a))
      END DO

      WRITE (7, 56) dash
      WRITE (7, *) 'Electron form factors: '
      DO a = 1, num_atoms_in_cell
         WRITE (7, *) (TRIM(element_str(e, a)), NEW_LINE('a'), eff(:, e + ielement_in_cell(a)), e=1, num_elements(a))
      END DO
   END IF

   CLOSE (7)
   mem_unit = 4D0*num_boxes*num_elements_in_cell*PRODUCT(num_cells)/1024D0**3D0 ! footprint for c_r
   mem_footprint = mem_unit + 3D0*mem_unit ! footprint for c_r and u_r
   IF (use_occ) mem_footprint = mem_footprint + 2D0*mem_unit ! plus a_q
   IF (use_mag) mem_footprint = mem_footprint + 3D0*mem_unit + 6D0*mem_unit ! plus spins, s_q
   IF (mem_footprint > mem_warning) THEN ! if >1GB needed to store c_r and u_r
      WRITE (*, *) ''
      WRITE (*, *) 'WARNING: A large amount of memory is required to store supercells and Fourier transforms.'
      WRITE (*, '(a, f8.4)') 'Memory required (GB): ', mem_footprint
      WRITE (*, *) 'Please check you have sufficient system resources, '
      WRITE (*, *) 'or consider using fewer/smaller supercells.'

   END IF

   ! transform nq and origin
   nq(:) = 2*nq(:) + 1
   num_q = PRODUCT(nq)
   mem_footprint = 4D0*num_q/1024D0**3D0
   IF (mem_footprint > mem_warning) THEN ! if >1GB needed to store all intensities
      WRITE (*, *) ''
      WRITE (*, *) 'WARNING: A large amount of memory is required to store the calculated intensities.'
      WRITE (*, '(a, f8.4)') 'Memory required (GB): ', mem_footprint
      WRITE (*, *) 'Please check you have sufficient system resources, '
      WRITE (*, *) 'or consider using a coarser reciprocal-space grid.'
      WRITE (*, *) ''
   END IF

   origin(:) = centre(:) - (axis(1, :) + axis(2, :) + axis(3, :))
   DO i = 1, 3
      origin_orth(i) = centre_orth(i) - axis_orth(i, i)
      IF (nq(i) > 1) THEN
         separation(i, :) = 2D0*axis(i, :)/(1D0*nq(i) - 1D0) ! interval between points
         spacing_orth(i) = 2D0*axis_length(i)/(1D0*nq(i) - 1D0)
      ELSE
         separation(i, :) = 0D0
         spacing_orth(i) = 0D0
      END IF
   END DO

   !! READ DATA
   IF (fit_data) THEN
      num_q = 0
      ! Find size of each dataset
      data_file = 'scatty_data_01.txt'
      OPEN (11, file=data_file, status='old')
      DO
         READ (11, *, END=20) dummy
         num_q = num_q + 1
      END DO
20    CLOSE (11)

      ALLOCATE (intensity_expt(num_q))
      ALLOCATE (errors_sq(num_q))
      ALLOCATE (qvec_from_iq(3, num_q))
      ALLOCATE (q_mod_all(num_q))
      ALLOCATE (valid(num_q))

      intensity_expt = 0D0
      errors_sq = 0D0
      qvec_from_iq = 0D0
      q_mod_all = 0D0
      valid = .FALSE.
      xtal_i_sq_sum = 0D0

      ! Read in experimental xtal data (format: h,k,l,intensity,error)
      OPEN (11, file=data_file, status='old')
      DO iq = 1, num_q
         READ (11, *) qvec_from_iq(:, iq), intensity_expt(iq), errors_sq(iq)
         IF (errors_sq(iq) < 1D-8) STOP 'DATA error: error bars seem to be <= zero for some XTAL data points'
         valid(iq) = .TRUE.
         q_mod_all(iq) = SQRT(DOT_PRODUCT(qvec_from_iq(:, iq), MATMUL(q_metric, qvec_from_iq(:, iq))))
         errors_sq(iq) = 1D0/(errors_sq(iq)**2D0)
         xtal_i_sq_sum = errors_sq(iq)*intensity_expt(iq)**2D0 + xtal_i_sq_sum
      END DO
      WRITE (*, *) 'Data file ', TRIM(ADJUSTL(data_file)), ' read successfully!'
      WRITE (*, *) num_q, 'data points read'
      CLOSE (11)

      axis2 = 2D0*axis
      ! write(*,*) 'origin', origin
      ! do i=1,3
      ! write(*,*) 'axis',i,axis2(:,i),nq(i)
      ! enddo
      ! stop
      ! END READ DATA
      out_str = TRIM(ADJUSTL(filename_stem))//'_'//TRIM(name)//'_data'
      CALL write_fit(num_q, nq, qvec_from_iq, q_orth, hkl2xyz, valid, origin, axis2, intensity_expt, title, &
      & out_str, colourmap, min, max, output_vtk, output_ascii, output_ppm, suppress_errors)
      DO iq = 1, num_q
         qvec_from_iq(:, iq) = qvec_from_iq(:, iq)*num_cells
      END DO
   ELSE
      ALLOCATE (qvec_from_iq(3, num_q))
      ALLOCATE (valid(num_q))
      valid = .TRUE.
      qvec_from_iq = 0D0
      iq = 0
      DO qz = 1, nq(3)
         DO qy = 1, nq(2)
            DO qx = 1, nq(1)
               iq = iq + 1
               qvec_from_iq(:, iq) = origin(:) + (qx - 1)*separation(1, :) + (qy - 1)*separation(2, :) + (qz - 1)*separation(3, :)
               qvec_from_iq(:, iq) = qvec_from_iq(:, iq)*num_cells
            END DO
         END DO
      END DO

   END IF
   WRITE (i, 56) dash
   WRITE (*, *) 'Identifying symmetry-equivalent peaks...'
   t1 = OMP_GET_WTIME()

   ! NOTE
   ! qx,qy,qz are coordinates along axes of user-specified cut for S(Q)
   ! h,k,l are coordinates along Cartesian axes in reciprocal space
   ! (in non-orthogonal case, these are the crystallographic reciprocal-lattice axes)

   max_bragg(:) = -HUGE(1)
   DO iq = 1, num_q
      q = qvec_from_iq(:, iq)
      IF (n == 0) THEN
         i0 = NINT(q)
      ELSE
         i0(:) = ABS(FLOOR(q)) + n
         dist_from_bragg = MOD(q, 1D0)
         WHERE (dist_from_bragg > 0.5D0) dist_from_bragg = dist_from_bragg - 1D0
         WHERE (dist_from_bragg < -0.5D0) dist_from_bragg = dist_from_bragg + 1D0
         WHERE (ABS(dist_from_bragg) < eps) i0 = ABS(NINT(q)) + n ! Takes account of rounding error
      END IF
      DO isymm = 1, num_symm_ops(laue_class)
         CALL spinsym_int(laue_class, isymm, i0, i1)
         i1 = ABS(i1)
         WHERE (i1 > max_bragg) max_bragg = i1
      END DO
   END DO

   min_bragg(:) = -max_bragg(:) ! this may lead to more Bragg points than necessary being stored, but makes for easier calculation
   num_hkl = PRODUCT(min_bragg + max_bragg + 1)

   mem_footprint = 4D0*PRODUCT(max_bragg - min_bragg + 1)/1024D0**3D0
   IF (mem_footprint > mem_warning) THEN ! if >2GB needed to store all supercell Bragg peaks
      WRITE (*, *) ''
      WRITE (*, *) 'WARNING: A large amount of memory is required to store supercell Bragg intensities.'
      WRITE (*, '(a, f8.4)') 'Memory required (GB): ', mem_footprint
      WRITE (*, *) 'Please check you have sufficient system resources, '
      WRITE (*, *) 'or consider using a smaller reciprocal-space range.'
   END IF

   ! WORK OUT WHICH h,k,l ARE NEEDED TO CALCULATE SPECIFIED CUT
   ! For each Q-point, this finds the (m)^3 surrounding h,k,l which are needed to determine its value, and remembers them
   ! Then, for each h,k,l in the above list, it finds all the symmetry-equivalent hkl and remembers them too.

   ALLOCATE (keep_hkl(min_bragg(1):max_bragg(1), min_bragg(2):max_bragg(2), min_bragg(3):max_bragg(3)))
   ALLOCATE (keep_hkl_symm(min_bragg(1):max_bragg(1), min_bragg(2):max_bragg(2), min_bragg(3):max_bragg(3)))

   ! find which Bragg points are needed to calculate S(Q), neglecting Laue symmetry
   keep_hkl = .FALSE.
   DO iq = 1, num_q
      q = qvec_from_iq(:, iq)
      IF (n == 0) THEN
         i0 = NINT(q)
         keep_hkl(i0(1), i0(2), i0(3)) = .TRUE.
      ELSE
         i0(:) = FLOOR(q) - n
         dist_from_bragg = MOD(q, 1D0)
         WHERE (dist_from_bragg > 0.5D0) dist_from_bragg = dist_from_bragg - 1D0
         WHERE (dist_from_bragg < -0.5D0) dist_from_bragg = dist_from_bragg + 1D0
         WHERE (ABS(dist_from_bragg) < eps) i0 = NINT(q) - n ! Takes account of rounding error
         i_min(:) = i0(:) + 1
         i_max(:) = i0(:) + m
         keep_hkl(i_min(1):i_max(1), i_min(2):i_max(2), i_min(3):i_max(3)) = .TRUE.
      END IF
   END DO

   ! find which Bragg points are equivalent by Laue symmetry to those determined above
   keep_hkl_symm = .FALSE.
   DO l = min_bragg(3), max_bragg(3)
      DO k = min_bragg(2), max_bragg(2)
         DO h = min_bragg(1), max_bragg(1)
            IF (keep_hkl(h, k, l)) THEN
            DO isymm = 1, num_symm_ops(laue_class)
               CALL spinsym_int(laue_class, isymm, [h, k, l], i1)
               ! next 3 lines take into account Friedel's law (-h,-k,-l = h,k,l)
               ! to avoid calculating same intensities twice
               IF (i1(1) < 0) CYCLE
               IF ((i1(1) == 0) .AND. (i1(2) < 0)) CYCLE
               IF ((i1(1) == 0) .AND. (i1(2) == 0) .AND. (i1(3) < 0)) CYCLE
               keep_hkl_symm(i1(1), i1(2), i1(3)) = .TRUE.
            END DO
            END IF
         END DO
      END DO
   END DO

   q_max = 0D0
   DO l = min_bragg(3), max_bragg(3)
      DO k = min_bragg(2), max_bragg(2)
         DO h = min_bragg(1), max_bragg(1)
            IF (.NOT. keep_hkl_symm(h, k, l)) CYCLE
            hkl = [h, k, l]
            hkl = hkl/num_cells
            q_mod = SQRT(DOT_PRODUCT(hkl, MATMUL(q_metric, hkl)))
            IF (q_mod > q_max) q_max = q_mod
         END DO
      END DO
   END DO

   num_hkl = COUNT(keep_hkl_symm) ! count number of points we need to keep

   t2 = OMP_GET_WTIME()
   WRITE (*, 50) ' ...done in CPU time: ', t2 - t1
   WRITE (*, '(a, i8)') ' Number of supercell peaks to be calculated:', num_hkl
   WRITE (*, '(a, i8)') ' Number of reciprocal-space points to be sampled:', num_q
   WRITE (*, 56) dash

   ALLOCATE (c_r(num_cells(1), num_cells(2), num_cells(3), num_elements_in_cell, num_boxes))
   ALLOCATE (u_r(num_cells(1), num_cells(2), num_cells(3), 3, num_elements_in_cell, num_boxes))

   IF (use_mag) THEN
      ALLOCATE (spins(num_cells(1), num_cells(2), num_cells(3), 3, num_elements_in_cell, num_boxes))
      ALLOCATE (tot_moment_sq(num_elements_in_cell))
      spins = 0D0
   END IF
   IF (temp_subtract) tot_moment_sq = 0D0

   WRITE (*, *) 'Reading supercells...'
   WRITE (*, *) 'Entries read are ATOM and SPIN.'

   t1 = OMP_GET_WTIME()

   IF (use_expansion_max_error .AND. use_expansion_order) THEN
      ex = 1D0/(1D0 + 1D0*max_taylor)
      min_u_cutoff = ((expansion_max_error*GAMMA(1D0*max_taylor + 2D0))**ex)/q_max
   ELSE
      min_u_cutoff = HUGE(0D0)
   END IF

   ALLOCATE (max_num_dft(num_elements_in_cell, num_boxes))
   ALLOCATE (u_rms(num_elements_in_cell))
   ALLOCATE (count_rms(num_elements_in_cell))
   max_num_dft = 0
   u_r = 0D0
   c_r = 0
   supercell = 0
   use_disp = .FALSE.
   u_max = 0D0
   u_rms = 0D0
   count_rms = 0
   ! read spin vectors from each spin config (N.B. max. 9999 configs)
   ibox = 0
   DO k = 1, 9999
      atom_file = TRIM(filename_stem)//'_atoms_'//TRIM(itc(k))//'.txt'
      spin_file = TRIM(filename_stem)//'_spins_'//TRIM(itc(k))//'.txt'
      INQUIRE (file=spin_file, exist=spin_file_exists)
      INQUIRE (file=atom_file, exist=atom_file_exists)
      IF (atom_file_exists) THEN
         read_file = atom_file
         file_exists = .TRUE.
      ELSEIF (spin_file_exists) THEN
         read_file = spin_file
         file_exists = .TRUE.
      ELSE
         file_exists = .FALSE.
      END IF
      IF (file_exists) THEN
         ibox = ibox + 1
         OPEN (11, file=read_file, status='old')
         DO
            READ (11, '(a)', END=12) str
            p = INDEX(str, 'ATOM')
            IF (p .NE. 0) THEN
               IF (use_occ) THEN
                  READ (str(p + 4:), *) supercell(:), s_i, elem
                  found_element = .FALSE.
                  DO e = 1, num_elements(supercell(1))
                     IF (TRIM(elem) == TRIM(element_str(e, supercell(1)))) THEN
                        found_element = .TRUE.
                        EXIT
                     END IF
                  END DO
                  IF (.NOT. found_element) THEN
                     WRITE (*, *) 'Supercell error: Element symbol "', TRIM(elem), &
                     & '" not recognised in "', TRIM(read_file), '" for ATOM ', supercell(:)
                     STOP
                  END IF
               ELSE
                  READ (str(p + 4:), *) supercell(:), s_i
                  e = 1
               END IF
               IF (ABS(DOT_PRODUCT(s_i, s_i)) > eps) use_disp = .TRUE.
               IF ((supercell(1) < 1) .OR. (supercell(1) > num_atoms_in_cell)) THEN
                  WRITE (*, *) 'Supercell error: Invalid entry in "', TRIM(read_file), '" for ATOM ', supercell(:)
                  WRITE (*, *) 'Entry 1 should be >= 1 and <= the number of entries for SITE.'
                  STOP
               END IF
               DO j = 2, 4
                  IF ((supercell(j) < 0) .OR. (supercell(j) > num_cells(j - 1) - 1)) THEN
                     WRITE (*, *) 'Supercell error: Invalid entry in "', TRIM(read_file), '" for ATOM ', supercell(:)
                     WRITE (*, *) 'Entries 2:4 should be >= 0 and <= BOX-1.'
                     STOP
                  END IF
               END DO
               i = e + ielement_in_cell(supercell(1))
               u_mod = DOT_PRODUCT(s_i, MATMUL(r_metric, s_i))
               u_mod = SQRT(u_mod)
               IF (clip_displacements .AND. u_mod > clip_dist) THEN
                  s_i = clip_dist*s_i/u_mod
                  u_mod = clip_dist
               END IF
               u_rms(i) = u_rms(i) + u_mod**2D0
               count_rms(i) = count_rms(i) + 1
               IF (u_mod > u_max) u_max = u_mod
               u_r(supercell(2) + 1, supercell(3) + 1, supercell(4) + 1, :, i, ibox) = s_i ! displacement
               c_r(supercell(2) + 1, supercell(3) + 1, supercell(4) + 1, i, ibox) = 1 ! occupancy
               IF ((use_expansion_max_error) .AND. (u_mod .GE. min_u_cutoff)) max_num_dft(i, ibox) = max_num_dft(i, ibox) + 1
            END IF
            p = INDEX(str, 'SPIN')
            IF (p .NE. 0) THEN
               IF (use_occ) THEN
                  READ (str(p + 4:), *) supercell(:), s_i, elem
                  found_element = .FALSE.
                  DO e = 1, num_elements(supercell(1))
                     IF (TRIM(elem) == TRIM(element_str(e, supercell(1)))) THEN
                        found_element = .TRUE.
                        EXIT
                     END IF
                  END DO
                  IF (.NOT. found_element) THEN
                     WRITE (*, *) 'Supercell error: Element symbol "', TRIM(elem), &
                     & '" not recognised in "', TRIM(read_file), '" for SPIN ', supercell(:)
                     STOP
                  END IF
               ELSE
                  READ (str(p + 4:), *) supercell(:), s_i
                  e = 1
               END IF
               IF ((supercell(1) < 1) .OR. (supercell(1) > num_atoms_in_cell)) THEN
                  WRITE (*, *) 'Supercell error: Invalid entry in "', TRIM(read_file), '" for SPIN ', supercell(:)
                  WRITE (*, *) 'Entry 1 should be >= 1 and <= the number of entries for SITE.'
                  STOP
               END IF
               DO j = 2, 4
                  IF ((supercell(j) < 0) .OR. (supercell(j) > num_cells(j - 1) - 1)) THEN
                     WRITE (*, *) 'Supercell error: Invalid entry in "', TRIM(read_file), '" for SPIN ', supercell(:)
                     WRITE (*, *) 'Entries 2:4 should be >= 0 and <= BOX-1.'
                     STOP
                  END IF
               END DO
               i = e + ielement_in_cell(supercell(1))
               s_i = MATMUL(r_orth, s_i/cell_params)
               c_r(supercell(2) + 1, supercell(3) + 1, supercell(4) + 1, i, ibox) = 1
               IF (temp_subtract) tot_moment_sq(i) = DOT_PRODUCT(s_i, s_i) + tot_moment_sq(i)
               IF (use_mag) THEN
                  spins(supercell(2) + 1, supercell(3) + 1, supercell(4) + 1, :, i, ibox) = s_i
               ELSE
                  IF (xtal_radiation == 1) THEN
                     WRITE (*, *) 'Supercell error: Please specify FORM_FACTOR_J0 in "', &
                     & TRIM(read_file), '" to calculate magnetic scattering.'
                     STOP
                  END IF
               END IF
            END IF
         END DO
12       CLOSE (11)
         WRITE (*, *) 'Supercell "', TRIM(read_file), '" read successfully!'
      END IF
   END DO
   IF (temp_subtract) tot_moment_sq = 2D0*tot_moment_sq/(3D0*num_boxes)

   t2 = OMP_GET_WTIME()

   WRITE (*, 50) ' ...done in CPU time: ', t2 - t1
   WRITE (*, 56) dash
   WRITE (*, '(a,l)') ' Occupational disorder present? ', use_occ
   WRITE (*, '(a,l)') ' Displacive disorder present? ', use_disp
   WRITE (*, '(a,l)') ' Magnetic disorder present? ', use_mag

   ! Check for inconsistencies in input files
   ALLOCATE (conc_recalc(num_elements_in_cell))
   DO i = 1, num_elements_in_cell
      conc_recalc(i) = SUM(c_r(:, :, :, i, :))/(1D0*num_boxes*tot_num_cells)
   END DO
   DO a = 1, num_atoms_in_cell
      DO e = 1, num_elements(a)
         i = e + ielement_in_cell(a)
         IF ((.NOT. mag_only) .AND. (.NOT. use_occ) .AND. (ABS(conc_recalc(i) - 1D0) > eps) &
         & .AND. (ABS(conc_recalc(i)) > eps)) THEN
            WRITE (*, *) 'Supercell error: Integer OCC (0 or 1) is specified for all sites, &
            & but supercell contains fractional occupancies.'
            STOP
         END IF
         IF (ABS(conc(i) - conc_recalc(i)) > 0.05D0) THEN
            WRITE (*, *) ''
            WRITE (*, *) 'WARNING: OCC differs from supercell-averaged occupancy by more than 0.05 for site: ', a
            WRITE (*, *) 'OCC: ', TRIM(element_str(e, a)), conc(i),&
            & ' Supercell average: ', TRIM(element_str(e, a)), conc_recalc(i)
            WRITE (*, *) 'Please check that the occupancies are consistent '
            WRITE (*, *) 'up to statistical (finite box-size) variations.'
         END IF
      END DO
   END DO

   IF (use_disp) THEN
      IF ((.NOT. use_expansion_order) .AND. (.NOT. use_expansion_max_error)) THEN
         WRITE (*, *) 'Config error: for displacive disorder, EXPANSION_ORDER (=> 0) &
         & and/or EXPANSION_MAX_ERROR must be specified.'
         STOP
      END IF
      num_q_bins = INT(q_max*nq_bin) + 1

      ALLOCATE (q_max_bin(num_q_bins))
      ALLOCATE (u_min_bin(num_q_bins))
      ALLOCATE (n_taylor_bin(num_q_bins))
      u_rms = SQRT(u_rms/count_rms)
      u_rms_max = MAXVAL(u_rms)/vol**(1D0/3D0)
      ! max. rms displacement as a fraction of the unit-cell length
      ! u_rms_max=1d0
      q_max_bin = 0D0
      u_min_bin = HUGE(0D0)
      n_taylor_bin = max_taylor

      WRITE (*, '(a, f10.6)') ' Maximum displacement (Ang.):  ', u_max
      WRITE (*, '(a, f10.6)') ' Maximum wavevector (inv. Ang.):  ', q_max
      found_n_taylor = .FALSE.
      IF (use_expansion_max_error) THEN

         IF (max_taylor == 0) THEN
            WRITE (*, *) 'Config error: EXPANSION_ORDER must be > 0 to use EXPANSION_MAX_ERROR.'
            STOP
         ELSE
            DO i = 1, max_taylor
               max_err = (u_max*q_max)**(1D0*i + 1D0)/GAMMA(1D0*i + 2D0)
               IF (max_err .LE. expansion_max_error) THEN
                  n_taylor = i
                  found_n_taylor = .TRUE.
                  EXIT
               END IF
            END DO
         END IF
         IF (.NOT. found_n_taylor) n_taylor = max_taylor
         ex = 1D0/(1D0 + 1D0*n_taylor)
         DO iq_bin = 1, num_q_bins
            q_max_bin(iq_bin) = iq_bin/(1D0*nq_bin)
            IF (q_max_bin(iq_bin) > q_max) q_max_bin(iq_bin) = q_max
            DO i = 1, n_taylor
               max_err = (u_max*q_max_bin(iq_bin))**(1D0*i + 1D0)/GAMMA(1D0*i + 2D0)
               IF (max_err .LE. expansion_max_error) THEN
                  n_taylor_bin(iq_bin) = i
                  EXIT
               END IF
            END DO
            ! write(*,*) expansion_max_error,n_taylor,q_max_bin(iq_bin),ex,Gamma(1d0*n_taylor+2d0)
            u_min_bin(iq_bin) = ((expansion_max_error*GAMMA(1D0*n_taylor + 2D0))**ex)/q_max_bin(iq_bin)
            ! Displacement magnitude above which FFT may give inaccurate results for G = q_max_bin(iq_bin)
         END DO
         IF (found_n_taylor) THEN
            use_direct_ft = .FALSE.
            WRITE (*, '(a, f10.7)') ' Maximum error in calculation of exp(iG.u): ', expansion_max_error
            WRITE (*, '(a, i4)') ' Maximum number of terms in Taylor expansion: ', MAXVAL(n_taylor_bin)
            IF (n_taylor > 19) THEN
               WRITE (*, *) ''
               WRITE (*, *) 'NOTE: There is a large number of terms in the Taylor expansion. '
               WRITE (*, *) 'This may result in a slow calculation.'
            END IF
         ELSE

            !  stop 'Direct Fourier calculation is disabled in this version'
            use_direct_ft = .TRUE.
            WRITE (*, *) ''
            WRITE (*, '(a, f10.7)') ' Maximum error in calculation of exp(iG.u): ', expansion_max_error
            WRITE (*, *) ''
            WRITE (*, *) 'NOTE: SCATTY will perform some calculations by direct Fourier '
            WRITE (*, *) 'summation to ensure that error <= EXPANSION_MAX_ERROR.'
            WRITE (*, *) 'This may result in a slower calculation.'
            max_err = expansion_max_error
            n_taylor = max_taylor
            ALLOCATE (b_r_dft_list(MAXVAL(max_num_dft), num_elements_in_cell, num_boxes, num_q_bins))
            ALLOCATE (r_dft_list(3, MAXVAL(max_num_dft), num_elements_in_cell, num_boxes, num_q_bins))
            ALLOCATE (u_r_dft_list(3, MAXVAL(max_num_dft), num_elements_in_cell, num_boxes, num_q_bins))
            ALLOCATE (num_dft(num_elements_in_cell, num_boxes, num_q_bins))
            ALLOCATE (num_hkl_bin(num_q_bins))
            b_r_dft_list = 0D0
            r_dft_list = 0D0
            u_r_dft_list = 0D0
            num_dft = 0
            DO ibox = 1, num_boxes
               DO a = 1, num_atoms_in_cell
                  DO e = 1, num_elements(a)
                     i = e + ielement_in_cell(a)
                     DO z = 1, num_cells(3)
                        DO y = 1, num_cells(2)
                           DO x = 1, num_cells(1)
                              u_mod = SQRT(DOT_PRODUCT(u_r(x, y, z, :, i, ibox), MATMUL(r_metric, u_r(x, y, z, :, i, ibox))))
                              DO iq_bin = 1, num_q_bins
                                 IF (u_mod > u_min_bin(iq_bin)) THEN
                                    num_dft(i, ibox, iq_bin) = num_dft(i, ibox, iq_bin) + 1
                                    IF (num_dft(i, ibox, iq_bin) > max_num_dft(i, ibox)) STOP 'This should not happen'
                                    b_r_dft_list(num_dft(i, ibox, iq_bin), i, ibox, iq_bin) = c_r(x, y, z, i, ibox)/conc(i)
                                    r_dft_list(:, num_dft(i, ibox, iq_bin), i, ibox, iq_bin) = [x - 1, y - 1, z - 1]
                                    u_r_dft_list(:, num_dft(i, ibox, iq_bin), i, ibox, iq_bin) = u_r(x, y, z, :, i, ibox)
                                 END IF
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
            num_hkl_bin = 0
            DO l = min_bragg(3), max_bragg(3)
               DO k = min_bragg(2), max_bragg(2)
                  DO h = min_bragg(1), max_bragg(1)
                     IF (.NOT. keep_hkl_symm(h, k, l)) CYCLE
                     hkl = [h, k, l]
                     hkl = hkl/num_cells
                     q_mod = SQRT(DOT_PRODUCT(hkl, MATMUL(q_metric, hkl)))
                     iq_bin = INT(q_mod*nq_bin) + 1
                     IF (iq_bin > 0) THEN
                        num_hkl_bin(iq_bin) = num_hkl_bin(iq_bin) + 1
                     END IF
                  END DO
               END DO
            END DO
            tot_num_exp = 0D0
            DO iq_bin = 1, num_q_bins
               tot_num_exp = tot_num_exp + num_hkl_bin(iq_bin)*SUM(num_dft(:, :, iq_bin))/(1D0*num_hkl*num_atoms*num_boxes)
            END DO
            WRITE (*, *) ''
            WRITE (*, '(a, f8.4)') ' Fraction of calculation to be performed by direct Fourier summation: ', tot_num_exp
         END IF

      ELSE

         use_direct_ft = .FALSE.
         n_taylor = max_taylor
         n_taylor_bin = max_taylor
         max_err = (u_max*q_max)**(1D0*n_taylor + 1D0)/GAMMA(1D0*n_taylor + 2D0)
         WRITE (*, '(a, i4)') ' Number of terms in Taylor displacement expansion: ', n_taylor
         WRITE (*, '(a, f8.4)') ' Estimated maximum error in exp(iG.u): ', max_err
         IF (n_taylor > 19) THEN
            WRITE (*, *) ''
            WRITE (*, *) 'NOTE: There is a large number of terms in the Taylor expansion.'
            WRITE (*, *) 'This may result in a slower calculation.'
         ELSE IF (n_taylor == 0) THEN
            WRITE (*, *) ''
            WRITE (*, *) 'WARNING: Number of terms in Taylor displacement expansion is zero.'
            WRITE (*, *) 'Scattering from displacive disorder WILL NOT be calculated.'
         END IF
      END IF

      max_num_terms = 0

      ALLOCATE (num_products(n_taylor))

      DO itaylor = 1, n_taylor
         max_num_terms = max_num_terms + (itaylor + 1)*(itaylor + 2)/2
         num_products(itaylor) = (itaylor + 2)*(itaylor + 1)/2
      END DO

      mem_footprint = 2D0*max_num_terms*mem_unit ! this can get very large
      IF (max_err > 0.1D0) THEN
         ex = 1D0/(1D0 + 1D0*n_taylor)
         WRITE (*, *) ''
         WRITE (*, *) 'WARNING: Estimated maximum error in exp(iG.u) > 0.1. '
         WRITE (*, *) 'Results may be very inaccurate!'
         WRITE (*, '(a, f10.7)') 'Maximum error in exp(iG.u) > 0.1 for |G|(A^-1) > ', ((0.1D0*GAMMA(1D0*n_taylor + 2D0))**ex)/u_max
         WRITE (*, *) 'Please increase EXPANSION_ORDER or reduce '
         WRITE (*, *) 'EXPANSION_MAX_ERROR to obtain accurate results.'
      END IF
      IF (mem_footprint > mem_warning) THEN
         WRITE (*, *) ''
         WRITE (*, *) 'WARNING: A large amount of memory is required for Taylor expansion of displacements.'
         WRITE (*, '(a, f8.4)') ' Memory required (GB): ', mem_footprint
         WRITE (*, *) 'Please check you have sufficient system resources, or consider:'
         WRITE (*, *) ' - using fewer supercells;'
         WRITE (*, *) ' - using smaller supercells;'
         WRITE (*, *) ' - for specified EXPANSION_MAX_ERROR, reducing the reciprocal-space range.'
      END IF

      WRITE (*, 56) dash
      WRITE (*, *) 'Calculating products for Taylor expansion...'

      t1 = OMP_GET_WTIME()

      ALLOCATE (pre_taylor_dft(n_taylor))

      DO itaylor = 1, n_taylor
         pre_taylor_dft(itaylor) = (0D0, 1D0)**(1D0*itaylor)/GAMMA(1D0*itaylor + 1D0)
      END DO

      ALLOCATE (itaylor_from_iterm(max_num_terms))
      ALLOCATE (prod(n_taylor, max_num_terms))
      ALLOCATE (nests(n_taylor))
      prod = 0
      itaylor_from_iterm = 0
      nests = 0
      from = 1
      to = 3

      DO itaylor = 1, n_taylor
         itaylor_from_iterm(from:to) = itaylor
         from = to + 1
         to = to + (itaylor + 2)*(itaylor + 3)/2
      END DO

      nests(1) = 1

      DO iterm = 1, max_num_terms
         prod(:, iterm) = nests
         CALL update(1)
      END DO

      DEALLOCATE (nests)
      ALLOCATE (num_terms(num_q_bins))
      num_terms = 0

      DO iq_bin = 1, num_q_bins
         num_terms(iq_bin) = SUM(num_products(1:n_taylor_bin(iq_bin)))
      END DO

      ALLOCATE (pre_taylor(max_num_terms))

      DO iterm = 1, max_num_terms
         count_taylor(:) = 0
         itaylor = itaylor_from_iterm(iterm)
         DO j = 1, itaylor
            DO k = 1, 3
               IF (prod(j, iterm) == k) count_taylor(k) = count_taylor(k) + 1
            END DO
         END DO
         pre_taylor(iterm) = (0D0, 1D0)**(1D0*itaylor)/PRODUCT(GAMMA(1D0*count_taylor + 1D0), mask=count_taylor .NE. 0)
      END DO

      t2 = OMP_GET_WTIME()

      WRITE (*, 50) ' ...done in CPU time: ', t2 - t1
      WRITE (*, '(a, i4)') ' Maximum number of products in Taylor displacement expansion: ', MAXVAL(num_terms)
   ELSE
      n_taylor = 0
      max_num_terms = 0
   END IF

   IF (ALLOCATED(num_hkl_bin)) DEALLOCATE (num_hkl_bin)
   IF (ALLOCATED(max_num_dft)) DEALLOCATE (max_num_dft)
   IF (ALLOCATED(u_min_bin)) DEALLOCATE (u_min_bin)
   IF (ALLOCATED(q_max_bin)) DEALLOCATE (q_max_bin)
   IF (ALLOCATED(u_rms)) DEALLOCATE (u_rms)
   IF (ALLOCATED(count_rms)) DEALLOCATE (count_rms)
   IF (ALLOCATED(conc_recalc)) DEALLOCATE (conc_recalc)

   t1 = OMP_GET_WTIME()

   IF (use_occ) ALLOCATE (a_q(num_elements_in_cell, num_boxes, num_cells(1), num_cells(2), num_cells(3)))
   IF (use_disp) ALLOCATE (au_q(max_num_terms, num_elements_in_cell, num_boxes, num_cells(1), num_cells(2), num_cells(3)))
   IF (use_mag) ALLOCATE (s_q(3, num_elements_in_cell, num_boxes, num_cells(1), num_cells(2), num_cells(3)))

   t0 = OMP_GET_WTIME()

   ! CALCULATING FFTs
   WRITE (*, 56) dash

   !$OMP  PARALLEL DO DEFAULT(SHARED) &
   !$OMP  PRIVATE(a_r, u_r_scale, u_product, u_product_sgl, b_r, i,j,k, iterm, itaylor)
   DO ibox = 1, num_boxes
      WRITE (*, '(a, i5, a, i5, a, i5, a)') ' Calculating FFT for supercell:', ibox, &
      &'    (Thread', omp_get_thread_num() + 1, '    out of', omp_get_num_threads(), ')...'

      IF (use_occ) THEN
         ALLOCATE (a_r(num_cells(1), num_cells(2), num_cells(3)))
         ALLOCATE (b_r(num_cells(1), num_cells(2), num_cells(3)))
      END IF

      IF (use_disp) THEN
         ALLOCATE (u_product(num_cells(1), num_cells(2), num_cells(3)))
         ALLOCATE (u_product_sgl(num_cells(1), num_cells(2), num_cells(3)))
         ALLOCATE (u_r_scale(num_cells(1), num_cells(2), num_cells(3), 3))
      END IF

      IF (use_occ .AND. (.NOT. use_disp)) THEN
         DO i = 1, num_elements_in_cell
            IF (ABS(conc(i)) > eps) THEN
               a_r(:, :, :) = c_r(:, :, :, i, ibox)/conc(i) - 1D0
               a_q(i, ibox, :, :, :) = fft(dcmplx(a_r, 0D0), inv=.TRUE.)
            ELSE
               a_q(i, ibox, :, :, :) = 0D0
            END IF
         END DO

      ELSE IF (use_disp .AND. (.NOT. use_occ)) THEN
         DO i = 1, num_elements_in_cell
            u_r_scale(:, :, :, :) = u_r(:, :, :, :, i, ibox)/u_rms_max
            DO iterm = 1, max_num_terms
               u_product(:, :, :) = 1D0
               itaylor = itaylor_from_iterm(iterm)
               DO j = 1, itaylor
                  k = prod(j, iterm)
                  u_product(:, :, :) = u_r_scale(:, :, :, k)*u_product(:, :, :)
               END DO
               u_product_sgl = u_product
               au_q(iterm, i, ibox, :, :, :) = fft(dcmplx(u_product_sgl, 0D0), inv=.TRUE.)
            END DO
         END DO

      ELSE IF (use_disp .AND. use_occ) THEN
         DO i = 1, num_elements_in_cell
            u_r_scale(:, :, :, :) = u_r(:, :, :, :, i, ibox)/u_rms_max
            IF (ABS(conc(i)) > eps) THEN
               b_r(:, :, :) = c_r(:, :, :, i, ibox)/conc(i)
               a_r(:, :, :) = b_r(:, :, :) - 1D0
               a_q(i, ibox, :, :, :) = fft(dcmplx(a_r, 0D0), inv=.TRUE.)
            ELSE
               b_r(:, :, :) = 0D0
               a_q(i, ibox, :, :, :) = 0D0
            END IF
            DO iterm = 1, max_num_terms
               u_product(:, :, :) = 1D0
               itaylor = itaylor_from_iterm(iterm)
               DO j = 1, itaylor
                  k = prod(j, iterm)
                  u_product(:, :, :) = u_r_scale(:, :, :, k)*u_product(:, :, :)
               END DO
               u_product(:, :, :) = u_product(:, :, :)*b_r(:, :, :)
               u_product_sgl = u_product
               au_q(iterm, i, ibox, :, :, :) = fft(dcmplx(u_product_sgl, 0D0), inv=.TRUE.)
            END DO
         END DO

      END IF
      IF (use_mag) THEN
         DO i = 1, num_elements_in_cell
            DO j = 1, 3
               s_q(j, i, ibox, :, :, :) = fft(dcmplx(spins(:, :, :, j, i, ibox), 0D0), inv=.TRUE.)
            END DO
         END DO
      END IF

      IF (use_occ) THEN
         DEALLOCATE (a_r)
         DEALLOCATE (b_r)
      END IF

      IF (use_disp) THEN
         DEALLOCATE (u_product)
         DEALLOCATE (u_product_sgl)
         DEALLOCATE (u_r_scale)
      END IF

   END DO
   !$OMP END PARALLEL DO

   t2 = OMP_GET_WTIME()
   WRITE (*, 50) ' ...done in CPU time: ', t2 - t1

   IF (ALLOCATED(spins)) DEALLOCATE (spins)
   IF (ALLOCATED(c_r)) DEALLOCATE (c_r)
   IF (ALLOCATED(u_r)) DEALLOCATE (u_r)

   WRITE (*, 56) dash

   WRITE (*, *) 'Calculating intensity at supercell Bragg positions...'
   ALLOCATE (intensity_bragg(min_bragg(1):max_bragg(1), min_bragg(2):max_bragg(2), min_bragg(3):max_bragg(3)))
   intensity_bragg = 0D0
   t1 = OMP_GET_WTIME()
   num_it = num_hkl/25
   progress = 0
   it_count = 0

   !$OMP PARALLEL DEFAULT(SHARED) &
   !$OMP PRIVATE(h,k,l, hkl_fft, bragg, gamma_point, ihkl, hkl, q, q_sq, q_mod, ss, iq_bin, iq_sq, &
   !$OMP a, sf, e, i, j, struct_fac, mag_struct_fac, xff_hkl, mag_ff_hkl, temp_subtract_term, bragg_sf0, &
   !$OMP hkl_product, hkl_scale, iterm, itaylor, u_hkl_product, u_hkl_product_bragg, bragg_sf, ibox, ft, &
   !$OMP idft, q_dot_u, jtaylor, sf2, amplitude, mag_amplitude) FIRSTPRIVATE(it_count)
   ALLOCATE (sf(num_atoms_in_cell))
   ALLOCATE (struct_fac(num_elements_in_cell))
   ALLOCATE (xff_hkl(num_elements_in_cell))
   IF (use_mag) THEN
      ALLOCATE (mag_struct_fac(num_elements_in_cell))
      ALLOCATE (mag_ff_hkl(num_elements_in_cell))
   END IF
   IF (use_disp) THEN
      ALLOCATE (q_dot_u(n_taylor))
      ALLOCATE (hkl_product(max_num_terms))
      ALLOCATE (u_hkl_product(max_num_terms))
      ALLOCATE (u_hkl_product_bragg(max_num_terms))
   END IF

   CALL progress_bar(0)

   !$OMP DO COLLAPSE(3)
   DO l = min_bragg(3), max_bragg(3)
      DO k = min_bragg(2), max_bragg(2)
         DO h = min_bragg(1), max_bragg(1)
            IF (.NOT. keep_hkl_symm(h, k, l)) CYCLE

            it_count = it_count + 1

            IF (it_count == num_it) THEN
               it_count = 0
               !$OMP ATOMIC
               progress = progress + 1
               !$OMP END ATOMIC
               CALL progress_bar(progress)
            END IF

            hkl_fft = MOD([h, k, l], num_cells)
            bragg = .FALSE.
            IF (h == 0 .AND. k == 0 .AND. l == 0) THEN
               gamma_point = .TRUE.
            ELSE
               gamma_point = .FALSE.
            END IF
            IF (hkl_fft(1) == 0 .AND. hkl_fft(2) == 0 .AND. hkl_fft(3) == 0) THEN
               ihkl = [h, k, l]/num_cells
               IF (centring == 1) THEN
                  bragg = .TRUE.
               ELSEIF (centring == 2) THEN
                  IF (MOD(ihkl(1) + ihkl(2) + ihkl(3), 2) == 0) bragg = .TRUE.
               ELSEIF (centring == 3) THEN
                  IF ((MOD(ihkl(1) + ihkl(2), 2) == 0) &
                  & .AND. (MOD(ihkl(2) + ihkl(3), 2) == 0) &
                  & .AND. (MOD(ihkl(1) + ihkl(3), 2) == 0)) &
                  & bragg = .TRUE.
               ELSEIF (centring == 4) THEN
                  IF (MOD(-ihkl(1) + ihkl(2) + ihkl(3), 3) == 0) bragg = .TRUE.
               ELSEIF (centring == 5) THEN
                  IF (MOD(ihkl(1) + ihkl(2), 2) == 0) bragg = .TRUE.
               ELSEIF (centring == 6) THEN
                  IF (MOD(ihkl(2) + ihkl(3), 2) == 0) bragg = .TRUE.
               ELSEIF (centring == 7) THEN
                  IF (MOD(ihkl(1) + ihkl(3), 2) == 0) bragg = .TRUE.
               ELSEIF (centring == 8) THEN
                  IF (MOD(ihkl(1) - ihkl(2), 3) == 0) bragg = .TRUE.
               END IF
            END IF
            WHERE (hkl_fft < 0) hkl_fft = num_cells + hkl_fft
            hkl_fft = hkl_fft + 1
            hkl = [h, k, l]
            hkl = hkl/num_cells
            IF (xtal_radiation == 2 .OR. xtal_radiation == 3 .OR. use_mag .OR. use_disp) THEN
               q = MATMUL(q_orth, hkl)
               q_sq = DOT_PRODUCT(q, q)
               q_mod = SQRT(q_sq)
               ss = q_sq*q_sq_to_ss
               iq_bin = INT(q_mod*nq_bin) + 1
               IF (gamma_point) THEN
                  iq_sq = 0D0
               ELSE
                  iq_sq = 1D0/q_sq
               END IF
            END IF
            hkl = tpi*hkl

            DO a = 1, num_atoms_in_cell
               sf(a) = EXP((0D0, 1D0)*DOT_PRODUCT(hkl, cell_coords(:, a)) - 0.5D0*u_iso(a)*q_sq)
            END DO

            IF ((xtal_radiation == 1) .AND. (.NOT. mag_only)) THEN
               DO a = 1, num_atoms_in_cell
                  DO e = 1, num_elements(a)
                     i = e + ielement_in_cell(a)
                     struct_fac(i) = sf(a)*nsl(i)*conc(i)
                  END DO
               END DO
            ELSE IF (xtal_radiation == 2) THEN
               DO a = 1, num_atoms_in_cell
                  DO e = 1, num_elements(a)
                     i = e + ielement_in_cell(a)
                     xff_hkl(i) = xff(1, i)*EXP(-ss*xff(2, i)) + xff(3, i)*EXP(-ss*xff(4, i)) + &
                     &xff(5, i)*EXP(-ss*xff(6, i)) + xff(7, i)*EXP(-ss*xff(8, i)) + xff(9, i)
                     struct_fac(i) = sf(a)*xff_hkl(i)*conc(i)
                  END DO
               END DO
            ELSE IF (xtal_radiation == 3) THEN
               DO a = 1, num_atoms_in_cell
                  DO e = 1, num_elements(a)
                     i = e + ielement_in_cell(a)
                     xff_hkl(i) = eff(1, i)*EXP(-ss*eff(6, i)) + eff(2, i)*EXP(-ss*eff(7, i)) + &
                     & eff(3, i)*EXP(-ss*eff(8, i)) + eff(4, i)*EXP(-ss*eff(9, i)) &
                     & + eff(5, i)*EXP(-ss*eff(10, i))
                     struct_fac(i) = sf(a)*xff_hkl(i)*conc(i)
                  END DO
               END DO
            END IF

            IF (use_mag) THEN
               DO a = 1, num_atoms_in_cell
                  DO e = 1, num_elements(a)
                     i = e + ielement_in_cell(a)
                     mag_ff_hkl(i) = ff0(1, i)*EXP(-ss*ff0(2, i)) &
                     & + ff0(3, i)*EXP(-ss*ff0(4, i)) + ff0(5, i)*EXP(-ss*ff0(6, i)) + ff0(7, i)
                     mag_ff_hkl(i) = mag_ff_hkl(i) + ss*c2(i)*(ff2(1, i)*EXP(-ss*ff2(2, i)) &
                     & + ff2(3, i)*EXP(-ss*ff2(4, i)) + ff2(5, i)*EXP(-ss*ff2(6, i)) + ff2(7, i))
                     mag_struct_fac(i) = sf(a)*mag_ff_hkl(i)
                  END DO
               END DO
               IF (temp_subtract) THEN
                  temp_subtract_term = 0D0
                  DO a = 1, num_atoms_in_cell
                     DO e = 1, num_elements(a)
                        i = e + ielement_in_cell(a)
                        temp_subtract_term = tot_moment_sq(i)*mag_ff_hkl(i)**2D0 + temp_subtract_term
                     END DO
                  END DO
               END IF
            END IF

            IF (.NOT. mag_only) THEN

               IF (bragg) bragg_sf0 = tot_num_cells*SUM(struct_fac)

               IF (use_disp) THEN
                  hkl_product = 1D0
                  hkl_scale = u_rms_max*hkl
                  DO iterm = 1, num_terms(iq_bin)
                     itaylor = itaylor_from_iterm(iterm)
                     DO j = 1, itaylor
                        i = prod(j, iterm)
                        hkl_product(iterm) = hkl_scale(i)*hkl_product(iterm)
                     END DO
                     u_hkl_product(iterm) = hkl_product(iterm)*pre_taylor(iterm)
                     ! if(bragg.and.(mod(itaylor,2)==1)) then
                     ! u_hkl_product_bragg(iterm)=(0d0,0d0)
                     ! else
                     IF (bragg) u_hkl_product_bragg(iterm) = u_hkl_product(iterm)
                     ! endif
                  END DO
                  IF (bragg) THEN
                     bragg_sf = (0D0, 0D0)
                     DO ibox = 1, num_boxes
                        DO i = 1, num_elements_in_cell
                           ft = (0D0, 0D0)
                           DO iterm = 1, num_terms(iq_bin)
                              ft = ft + u_hkl_product_bragg(iterm)*au_q(iterm, i, ibox, hkl_fft(1), hkl_fft(2), hkl_fft(3))
                           END DO

                           IF (use_direct_ft) THEN
                              DO idft = 1, num_dft(i, ibox, iq_bin)
                                 q_dot_u(1) = DOT_PRODUCT(hkl, u_r_dft_list(:, idft, i, ibox, iq_bin))
                                 DO itaylor = 2, n_taylor, 2
                                    jtaylor = itaylor/2
                                    q_dot_u(itaylor) = q_dot_u(jtaylor)*q_dot_u(jtaylor)
                                    IF (itaylor < n_taylor) q_dot_u(itaylor + 1) = q_dot_u(1)*q_dot_u(itaylor)
                                 END DO
                                 sf2 = b_r_dft_list(idft, i, ibox, iq_bin) &
                                      & *EXP((0D0, 1D0)*DOT_PRODUCT(hkl, r_dft_list(:, idft, i, ibox, iq_bin)))
                                 ft = ft + sf2*(EXP((0D0, 1D0)*q_dot_u(1)) - 1D0 - SUM(pre_taylor_dft*q_dot_u))
                              END DO
                           END IF

                           bragg_sf = bragg_sf + ft*struct_fac(i)
                        END DO
                     END DO
                     bragg_sf = bragg_sf/num_boxes
                  END IF
               END IF

               IF ((.NOT. use_occ) .AND. (.NOT. use_disp)) THEN
                  IF (bragg) THEN

                     IF (.NOT. remove_bragg) THEN
                        intensity_bragg(h, k, l) = intensity_bragg(h, k, l) + num_boxes*CONJG(bragg_sf0)*bragg_sf0
                     END IF

                  END IF
               ELSE IF (use_occ .AND. (.NOT. use_disp)) THEN
                  DO ibox = 1, num_boxes
                     amplitude = (0D0, 0D0)
                     DO i = 1, num_elements_in_cell
                        amplitude = amplitude + a_q(i, ibox, hkl_fft(1), hkl_fft(2), hkl_fft(3))*struct_fac(i)
                     END DO
                     IF (bragg) THEN
                        IF (.NOT. remove_bragg) amplitude = amplitude + bragg_sf0
                     END IF

                     intensity_bragg(h, k, l) = intensity_bragg(h, k, l) + CONJG(amplitude)*amplitude

                  END DO
               ELSE IF (use_disp .AND. (.NOT. use_occ)) THEN
                  DO ibox = 1, num_boxes
                     amplitude = (0D0, 0D0)
                     DO i = 1, num_elements_in_cell

                        ft = (0D0, 0D0)
                        DO iterm = 1, num_terms(iq_bin)
                           ft = ft + u_hkl_product(iterm)*au_q(iterm, i, ibox, hkl_fft(1), hkl_fft(2), hkl_fft(3))
                        END DO

                        IF (use_direct_ft) THEN
                           DO idft = 1, num_dft(i, ibox, iq_bin)
                              q_dot_u(1) = DOT_PRODUCT(hkl, u_r_dft_list(:, idft, i, ibox, iq_bin))
                              DO itaylor = 2, n_taylor, 2
                                 jtaylor = itaylor/2
                                 q_dot_u(itaylor) = q_dot_u(jtaylor)*q_dot_u(jtaylor)
                                 IF (itaylor < n_taylor) q_dot_u(itaylor + 1) = q_dot_u(1)*q_dot_u(itaylor)
                              END DO
                              sf2 = b_r_dft_list(idft, i, ibox, iq_bin) &
                                   & *EXP((0D0, 1D0)*DOT_PRODUCT(hkl, r_dft_list(:, idft, i, ibox, iq_bin)))
                              ft = ft + sf2*(EXP((0D0, 1D0)*q_dot_u(1)) - 1D0 - SUM(pre_taylor_dft*q_dot_u))
                           END DO
                        END IF

                        amplitude = amplitude + ft*struct_fac(i)

                     END DO

                     IF (bragg) THEN
                        IF (remove_bragg) THEN
                           amplitude = amplitude - bragg_sf
                        ELSE
                           amplitude = amplitude + bragg_sf0
                        END IF
                     END IF

                     intensity_bragg(h, k, l) = intensity_bragg(h, k, l) + CONJG(amplitude)*amplitude

                  END DO
               ELSE IF (use_disp .AND. use_occ) THEN
                  DO ibox = 1, num_boxes
                     amplitude = (0D0, 0D0)
                     DO i = 1, num_elements_in_cell
                        ft = a_q(i, ibox, hkl_fft(1), hkl_fft(2), hkl_fft(3))
                        DO iterm = 1, num_terms(iq_bin)
                           ft = ft + u_hkl_product(iterm)*au_q(iterm, i, ibox, hkl_fft(1), hkl_fft(2), hkl_fft(3))
                        END DO

                        IF (use_direct_ft) THEN
                           DO idft = 1, num_dft(i, ibox, iq_bin)
                              q_dot_u(1) = DOT_PRODUCT(hkl, u_r_dft_list(:, idft, i, ibox, iq_bin))
                              DO itaylor = 2, n_taylor, 2
                                 jtaylor = itaylor/2
                                 q_dot_u(itaylor) = q_dot_u(jtaylor)*q_dot_u(jtaylor)
                                 IF (itaylor < n_taylor) q_dot_u(itaylor + 1) = q_dot_u(1)*q_dot_u(itaylor)
                              END DO
                              sf2 = b_r_dft_list(idft, i, ibox, iq_bin) &
                                   &*EXP((0D0, 1D0)*DOT_PRODUCT(hkl, r_dft_list(:, idft, i, ibox, iq_bin)))
                              ft = ft + sf2*(EXP((0D0, 1D0)*q_dot_u(1)) - 1D0 - SUM(pre_taylor_dft*q_dot_u))
                           END DO
                        END IF

                        amplitude = amplitude + ft*struct_fac(i)

                     END DO

                     IF (bragg) THEN
                        IF (remove_bragg) THEN
                           amplitude = amplitude - bragg_sf
                        ELSE
                           amplitude = amplitude + bragg_sf0
                        END IF
                     END IF

                     intensity_bragg(h, k, l) = intensity_bragg(h, k, l) + CONJG(amplitude)*amplitude

                  END DO
               END IF
            END IF

            IF (use_mag) THEN
               DO ibox = 1, num_boxes
                  mag_amplitude = (0D0, 0D0)
                  DO i = 1, num_elements_in_cell
                     mag_amplitude(:) = mag_amplitude(:) &
                     & + s_q(:, i, ibox, hkl_fft(1), hkl_fft(2), hkl_fft(3))*mag_struct_fac(i)
                  END DO
                  IF (gamma_point) THEN
                     mag_amplitude = mag_amplitude*SQRT(2D0/3D0)
                  ELSE
                     mag_amplitude = mag_amplitude - iq_sq*q*DOT_PRODUCT(q, mag_amplitude)
                  END IF
                  IF (temp_subtract) THEN
                     intensity_bragg(h, k, l) = intensity_bragg(h, k, l) &
                     & + mag_prefactor*(DOT_PRODUCT(mag_amplitude, mag_amplitude) - temp_subtract_term)
                  ELSE
                     intensity_bragg(h, k, l) = intensity_bragg(h, k, l) &
                     & + mag_prefactor*DOT_PRODUCT(mag_amplitude, mag_amplitude)
                  END IF
               END DO
            END IF

            IF (isnan(intensity_bragg(h, k, l))) STOP 'Calculation error: intensity is NaN.'

         END DO
      END DO
   END DO
   !$OMP END DO

   DEALLOCATE (sf)
   DEALLOCATE (struct_fac)
   DEALLOCATE (xff_hkl)
   IF (use_mag) THEN
      DEALLOCATE (mag_struct_fac)
      DEALLOCATE (mag_ff_hkl)
   END IF
   IF (use_disp) THEN
      DEALLOCATE (q_dot_u)
      DEALLOCATE (hkl_product)
      DEALLOCATE (u_hkl_product)
      DEALLOCATE (u_hkl_product_bragg)
   END IF

   !$OMP END PARALLEL
   CALL progress_bar(25)

   norm = 1D0/(1D0*num_boxes*num_atoms)
   intensity_bragg = intensity_bragg*norm
   t2 = OMP_GET_WTIME()

   WRITE (*, 50) ' ...done in CPU time: ', t2 - t1

   IF (ALLOCATED(a_q)) DEALLOCATE (a_q)
   IF (ALLOCATED(au_q)) DEALLOCATE (au_q)
   IF (ALLOCATED(s_q)) DEALLOCATE (s_q)
   IF (ALLOCATED(tot_moment_sq)) DEALLOCATE (tot_moment_sq)
   IF (ALLOCATED(itaylor_from_iterm)) DEALLOCATE (itaylor_from_iterm)
   IF (ALLOCATED(prod)) DEALLOCATE (prod)
   IF (ALLOCATED(pre_taylor)) DEALLOCATE (pre_taylor)
   IF (ALLOCATED(pre_taylor_dft)) DEALLOCATE (pre_taylor_dft)
   IF (ALLOCATED(num_terms)) DEALLOCATE (num_terms)
   IF (ALLOCATED(b_r_dft_list)) DEALLOCATE (b_r_dft_list)
   IF (ALLOCATED(r_dft_list)) DEALLOCATE (r_dft_list)
   IF (ALLOCATED(u_r_dft_list)) DEALLOCATE (u_r_dft_list)
   IF (ALLOCATED(num_dft)) DEALLOCATE (num_dft)

   ! SYMMETRISE
   WRITE (*, 56) dash
   WRITE (*, *) 'Applying symmetry...'
   t1 = OMP_GET_WTIME()

   ALLOCATE (intensity_laue(min_bragg(1):max_bragg(1), min_bragg(2):max_bragg(2), min_bragg(3):max_bragg(3)))
   intensity_laue = 0D0

   DO l = min_bragg(3), max_bragg(3)
      DO k = min_bragg(2), max_bragg(2)
         DO h = min_bragg(1), max_bragg(1)
         IF (keep_hkl(h, k, l)) THEN
            DO isymm = 1, num_symm_ops(laue_class)
               CALL spinsym_int(laue_class, isymm, [h, k, l], i1)
               intensity_laue(h, k, l) = intensity_laue(h, k, l) + intensity_bragg(i1(1), i1(2), i1(3))
            END DO
         END IF
         END DO
      END DO
   END DO

   ! multiply by 2 because 1/2 of intensities were not calculated above (since -h,-k,-l = h,k,l )
   intensity_laue = 2.0*intensity_laue/num_symm_ops(laue_class)
   intensity_laue(0, 0, 0) = intensity_bragg(0, 0, 0) ! must reset gamma point to avoid over-counting intensity
   t2 = OMP_GET_WTIME()
   WRITE (*, 50) ' ...done in CPU time: ', t2 - t1

   WRITE (*, 56) dash
   WRITE (*, *) 'Resampling reciprocal-space points...'
   t1 = OMP_GET_WTIME()
   ALLOCATE (intensity(num_q))
   intensity = 0.0

   ! Work out normalisation constants for resampling (should be ~1.0 but not exact)
   !$OMP PARALLEL DEFAULT(SHARED) &
   !$OMP PRIVATE(qz, qy, qx, i, j, k, box, window, q, norm, i0, &
   !$OMP dist_from_bragg, dhkl, weight, qr)
   ALLOCATE (box(m, 3))
   ALLOCATE (window(m, 3))

   !$OMP DO
   DO iq = 1, num_q

      q = qvec_from_iq(:, iq)
      norm = 0D0

      IF (n == 0) THEN
         i0(:) = NINT(q)
         intensity(iq) = intensity_laue(i0(1), i0(2), i0(3))
      ELSE
         i0(:) = FLOOR(q) - n
         dist_from_bragg = MOD(q, 1D0)
         WHERE (dist_from_bragg > 0.5D0) dist_from_bragg = dist_from_bragg - 1D0
         WHERE (dist_from_bragg < -0.5D0) dist_from_bragg = dist_from_bragg + 1D0
         WHERE (ABS(dist_from_bragg) < eps) i0 = NINT(q) - n ! Takes account of rounding error

         DO i = 1, m

            dhkl(:) = q(:) - (i0(:) + i)
            dhkl(:) = tpi*dhkl(:)
            IF (window_type == 1) THEN ! lanczos
               WHERE (ABS(dhkl) > eps)
                  window(i, :) = m*SIN(dhkl(:)/m)/dhkl(:)
               ELSEWHERE
                  window(i, :) = 1D0
               END WHERE
            ELSE IF (window_type == 2) THEN ! cosine
               WHERE (ABS(dhkl) > eps)
                  window(i, :) = COS(dhkl(:)/(2D0*m))
               ELSEWHERE
                  window(i, :) = 1D0
               END WHERE
            ELSE
               STOP 'Resampling error: invalid WINDOW_TYPE.'
            END IF
            IF (sum_type == 1) THEN ! cuboid
               WHERE (ABS(dhkl) > eps)
                  box(i, :) = 2D0*SIN(dhkl(:)*r_cut)/dhkl(:)
               ELSEWHERE
                  box(i, :) = 2D0*r_cut
               END WHERE
            END IF
         END DO

         DO k = 1, m
            DO j = 1, m
               DO i = 1, m
                  IF (sum_type == 1) THEN ! cuboid
                     weight = box(i, 1)*box(j, 2)*box(k, 3)
                  ELSE IF (sum_type == 2) THEN ! spheroid
                     dhkl(:) = q(:) - [i + i0(1), j + i0(2), k + i0(3)]
                     dhkl = MATMUL(q_orth, dhkl/num_cells)
                     qr = mr_cut*SQRT(DOT_PRODUCT(dhkl, dhkl))
                     IF (ABS(qr) > eps) THEN
                        weight = (SIN(qr)/qr - COS(qr))/(qr*qr)
                     ELSE
                        weight = 1D0/3D0
                     END IF
                     weight = weight*sphere_norm
                  ELSE
                     STOP 'Resampling error: invalid SUM_TYPE.'
                  END IF
                  weight = weight*window(i, 1)*window(j, 2)*window(k, 3)
                  intensity(iq) = weight*intensity_laue(i + i0(1), j + i0(2), k + i0(3)) + intensity(iq)
                  norm = weight + norm
               END DO
            END DO
         END DO

         intensity(iq) = ABS(intensity(iq))/norm
      END IF

   END DO
   !$OMP END DO

   DEALLOCATE (box)
   DEALLOCATE (window)

   !$OMP END PARALLEL
   t2 = OMP_GET_WTIME()
   WRITE (*, 50) ' ...done in CPU time: ', t2 - t1
   WRITE (*, 56) dash

   !! WRITE OUTPUT FILES
   out_str = TRIM(filename_stem)//'_'//TRIM(name)//'_sc'
   IF (output_supercell_bragg) THEN

      WRITE (*, *) 'Writing supercell Bragg peaks as .vtk file...'
      OPEN (11, file=TRIM(ADJUSTL(out_str))//'_supercell.vtk', status='replace')
      WRITE (11, '(a)') '# vtk DataFile Version 2.0'
      WRITE (11, '(a)') 'TITLE diffuse scattering - supercell Bragg peaks'
      WRITE (11, '(a)') 'ASCII'
      WRITE (11, '(a)') 'DATASET STRUCTURED_POINTS'
      WRITE (11, '(a,3i6)') ADJUSTL('DIMENSIONS'), (2*max_bragg(i) + 1, i=1, 3)
      WRITE (11, '(a,3f12.6)') 'ORIGIN', min_bragg/(1.0*num_cells)
      WRITE (11, '(a,3f12.6)') 'SPACING', 1.0/(1.0*num_cells)
      WRITE (11, '(a,i12)') 'POINT_DATA', PRODUCT(2*max_bragg + 1)
      WRITE (11, '(a)') 'SCALARS '//TRIM(ADJUSTL(title))//' float'
      WRITE (11, '(a)') 'LOOKUP_TABLE default'
      DO l = min_bragg(3), max_bragg(3)
         DO k = min_bragg(2), max_bragg(2)
            DO h = min_bragg(1), max_bragg(1)
               WRITE (11, '(f20.12)') intensity_bragg(h, k, l)
            END DO
         END DO
      END DO
      CLOSE (11)

   END IF

   DEALLOCATE (intensity_bragg)
   DEALLOCATE (intensity_laue)

   IF (fit_data) THEN

      CALL calc_chi_sq(intensity, intensity_expt, errors_sq, valid, num_q, q_mod_all, refine_scale, &
      & refine_flat_bgr, refine_linear_bgr, chi_sq, scale, flat_bgr, linear_bgr)
      WRITE (*, *) 'CHI^2, R_WP, SCALE, FLAT BGR, LINEAR BGR =', chi_sq, 100D0*SQRT(chi_sq/xtal_i_sq_sum),&
      &scale, flat_bgr, linear_bgr
      intensity = intensity*scale + flat_bgr + q_mod_all*linear_bgr
      WRITE (*, *) ''
      INQUIRE (file=TRIM(TRIM(name)//'_chisq_rwp_scale_fbgr_lbgr.txt'), exist=file_exists)
      IF (file_exists) THEN
         OPEN (11, file=TRIM(name)//'_chisq_rwp_scale_fbgr_lbgr.txt', status='old', position='append', action='write')
      ELSE
         OPEN (11, file=TRIM(name)//'_chisq_rwp_scale_fbgr_lbgr.txt', status='new', action='write')
      END IF
      WRITE (11, *) chi_sq, 100D0*SQRT(chi_sq/xtal_i_sq_sum), scale, flat_bgr, linear_bgr
      CLOSE (11)

   END IF
   axis = 2D0*axis
   DO iq = 1, num_q
      qvec_from_iq(:, iq) = qvec_from_iq(:, iq)/(1D0*num_cells)
   END DO
   IF (output_ascii) THEN
      OPEN (11, file=TRIM(out_str)//'_list.txt', status='replace')
      DO iq = 1, num_q
         WRITE (11, '(4f20.12, f15.5)') qvec_from_iq(:, iq), intensity(iq), 1.0
      END DO
      CLOSE (11)
   END IF

   CALL write_fit(num_q, nq, qvec_from_iq, q_orth, hkl2xyz, valid, origin, axis, intensity, title, &
   & out_str, colourmap, min, max, output_vtk, output_ascii, output_ppm, suppress_errors)

   WRITE (*, *) 'All done!'
   WRITE (*, 50) ' ...Total CPU time: ', t2 - t0

50 FORMAT(a, f8.4)
56 FORMAT(a)

CONTAINS

   RECURSIVE SUBROUTINE update(place_in_product)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: place_in_product
      INTEGER :: looper
      IF (nests(place_in_product) < 3) THEN
         nests(place_in_product) = nests(place_in_product) + 1
         DO looper = 1, place_in_product
            nests(looper) = nests(place_in_product)
         END DO
      ELSE
         IF (.NOT. (place_in_product == n_taylor)) CALL update(place_in_product + 1)
      END IF
   END SUBROUTINE

   SUBROUTINE progress_bar(iter_count)
      IMPLICIT NONE
      INTEGER(kind=4)::iter_count, looper
      CHARACTER(len=40)::bar = "???% |                         |"
      WRITE (unit=bar(1:3), fmt="(i3)") 4*iter_count
      DO looper = 1, iter_count
         bar(6 + looper:6 + looper) = "*"
      END DO
      WRITE (unit=6, fmt="(a1,a, $)") CHAR(13), bar
      RETURN
   END SUBROUTINE progress_bar

   SUBROUTINE read_config_short
      IMPLICIT NONE

      REAL(8), PARAMETER :: eps = 10D0*EPSILON(0.0)
      CHARACTER(200) :: str
      INTEGER :: i, j, k, l, m, n, o, e, nitems

      title = 'unspecified_title'
      num_atoms_in_cell = 0
      num_cells = 0
      cell_params = 0D0
      cell_angles = 0D0
      use_occ = .FALSE.
      use_mag = .FALSE.

      OPEN (11, file=read_file, status='old')
      DO
         READ (11, '(a)', END=22) str
         i = INDEX(str, 'TITLE')
         IF (i .NE. 0) READ (str(i + 5:), *) title
         i = INDEX(str, 'CELL')
         IF (i .NE. 0) READ (str(i + 4:), *) cell_params(:), cell_angles(:)
         i = INDEX(str, 'BOX')
         IF (i .NE. 0) READ (str(i + 3:), *) num_cells(:)
         i = INDEX(str, 'SITE')
         IF (i .NE. 0) num_atoms_in_cell = num_atoms_in_cell + 1
      END DO
22    REWIND (11)
      ALLOCATE (cell_coords(3, num_atoms_in_cell))
      ALLOCATE (u_iso(num_atoms_in_cell))
      ALLOCATE (num_elements(num_atoms_in_cell))
      cell_coords = 0D0
      u_iso = 0D0
      num_elements = 1
      j = 0
      k = 0
      l = 0
      DO
         READ (11, '(a)', END=23) str
         i = INDEX(str, 'SITE')
         IF (i .NE. 0) THEN
            j = j + 1
            READ (str(i + 4:), *) cell_coords(:, j)
         END IF
         i = INDEX(str, 'OCC')
         IF (i .NE. 0) THEN
            k = k + 1
            IF (k > num_atoms_in_cell) STOP 'Supercell error: Too many entries for OCC.'
            num_elements(k) = (nitems(str) - 1)/2
         END IF
         i = INDEX(str, 'UISO')
         IF (i .NE. 0) THEN
            l = l + 1
            IF (k > num_atoms_in_cell) STOP 'Supercell error: Too many entries for UISO.'
            READ (str(i + 4:), *) u_iso(l)
         END IF
      END DO
23    REWIND (11)
      IF (k == 1) num_elements(:) = num_elements(1)
      max_num_elements = MAXVAL(num_elements)
      ALLOCATE (element_str(max_num_elements, num_atoms_in_cell))

      num_elements_in_cell = SUM(num_elements)
      ALLOCATE (ielement_in_cell(num_atoms_in_cell))
      DO a = 1, num_atoms_in_cell
         ielement_in_cell(a) = SUM(num_elements(1:a - 1))
      END DO
      ALLOCATE (xff(9, num_elements_in_cell))
      ALLOCATE (eff(10, num_elements_in_cell))
      ALLOCATE (nsl(num_elements_in_cell))
      ALLOCATE (conc(num_elements_in_cell))
      ALLOCATE (ff0(7, num_elements_in_cell))
      ALLOCATE (ff2(7, num_elements_in_cell))
      ALLOCATE (c2(num_elements_in_cell))
      xff = 0D0
      eff = 0D0
      nsl = 0D0
      conc = 0D0 ! sites not occupied by default
      ff0 = 0D0
      ff2 = 0D0
      c2 = 0D0
      j = 0
      k = 0
      l = 0
      m = 0
      n = 0
      o = 0
      DO
         READ (11, '(a)', END=24) str
         p = INDEX(str, 'OCC')
         IF (p .NE. 0) THEN
            l = l + 1
            IF (l > num_atoms_in_cell) STOP 'Supercell error: Too many entries for OCC.'
            num_elements(l) = (nitems(str) - 1)/2
            ALLOCATE (item_index(nitems(str) + 1))
            CALL index_items(str, nitems(str), item_index)
            DO e = 1, num_elements(l)
               i = e + ielement_in_cell(l)
               READ (str(item_index(2*e):item_index(2*e + 1) - 1), *) element_str(e, l)
               CALL sl_from_element(element_str(e, l), xff(:, i), nsl(i), eff(:, i))
               READ (str(item_index(2*e + 1):item_index(2*e + 2) - 1), *) conc(i)
               IF ((conc(i) > eps) .AND. (conc(i) < 1D0 - eps)) THEN
                  use_occ = .TRUE. ! If fractional occupancy on any site, set use_occ to true
               ELSE IF (conc(i) < -eps) THEN
                  WRITE (*, *) 'Supercell error: OCC < 0 for site ', l
                  STOP
               ELSE IF (conc(i) > 1D0 + eps) THEN
                  WRITE (*, *) 'Supercell error: OCC > 1 for site ', l
                  STOP
               END IF
            END DO
            DEALLOCATE (item_index)
         END IF
         p = INDEX(str, 'FORM_FACTOR_J0')
         IF (p .NE. 0) THEN
            use_mag = .TRUE.
            m = m + 1
            IF (m > num_atoms_in_cell) STOP 'Supercell error: Too many entries for FORM_FACTOR_J0.'
            READ (str(p + 14:), *) (ff0(:, e + ielement_in_cell(m)), e=1, num_elements(m))
         END IF
         p = INDEX(str, 'FORM_FACTOR_J2')
         IF (p .NE. 0) THEN
            n = n + 1
            IF (n > num_atoms_in_cell) STOP 'Supercell error: Too many entries for FORM_FACTOR_J2.'
            READ (str(p + 14:), *) (ff2(:, e + ielement_in_cell(n)), e=1, num_elements(n))
         END IF
         p = INDEX(str, 'C2')
         IF (p .NE. 0) THEN
            o = o + 1
            IF (o > num_atoms_in_cell) STOP 'Supercell error: Too many entries for FORM_FACTOR_J2.'
            READ (str(p + 2:), *) (c2(e + ielement_in_cell(o)), e=1, num_elements(o))
         END IF
      END DO
24    CLOSE (11)

      IF (use_mag) THEN
         IF (xtal_radiation == 2 .OR. xtal_radiation == 3) THEN ! X-ray / electron
            WRITE (*, *) ''
            WRITE (*, *) 'WARNING: Neutron radiation is not specified, but supercell contains magnetic '
            WRITE (*, *) 'information. Magnetic scattering will not be calculated for X-ray or electron radiation.'
            WRITE (*, *) ''
            use_mag = .FALSE.
         END IF
         IF (m == 1) THEN
            WRITE (*, *) ''
            WRITE (*, *) 'NOTE: A single FORM_FACTOR_J0 line is given in "', TRIM(read_file), '".'
            WRITE (*, *) 'SCATTY will assume that all sites in the unit cell'
            WRITE (*, *) 'have the same FORM_FACTOR_J0, FORM_FACTOR_J2, and C2.'

            DO a = 1, num_atoms_in_cell
               DO e = 1, num_elements(a)
                  ff0(:, e + ielement_in_cell(a)) = ff0(:, e + ielement_in_cell(1))
                  ff2(:, e + ielement_in_cell(a)) = ff2(:, e + ielement_in_cell(1))
                  c2(e + ielement_in_cell(a)) = c2(e + ielement_in_cell(1))
               END DO
            END DO
         ELSE IF (m < num_atoms_in_cell) THEN
            WRITE (*, *) 'Supercell error: Please specify FORM_FACTOR_J0 in "', &
            & TRIM(read_file), '" to calculate magnetic neutron scattering.'
            STOP
         END IF
      END IF

      IF (l == 0) THEN
         IF (.NOT. use_mag) THEN
            WRITE (*, *) 'Supercell error: No entries for OCC or FORM_FACTOR_J0 in "', TRIM(read_file), '".'
            WRITE (*, *) 'Please specify OCC and/or FORM_FACTOR_J0 to calculate '
            WRITE (*, *) 'structural and/or magnetic scattering.'

            STOP
         ELSE
            WRITE (*, *) ''
            WRITE (*, *) 'WARNING: OCC not given in "', TRIM(read_file), '".'
            WRITE (*, *) 'SCATTY will only calculate magnetic scattering '
            WRITE (*, *) 'and will assume full occupancy for normalisation.'

            conc = 1D0 ! Assumes full occupancy
            mag_only = .TRUE.
         END IF
      ELSE IF (l == 1) THEN
         WRITE (*, *) ''
         WRITE (*, *) 'NOTE: A single OCC line is given in "', TRIM(read_file), '".'
         WRITE (*, *) 'SCATTY will assume that all sites in the unit cell have the same OCC.'
         DO e = 1, num_elements(1)
            element_str(e, :) = element_str(e, 1)
         END DO
         DO a = 1, num_atoms_in_cell
            DO e = 1, num_elements(a)
               i = e + ielement_in_cell(a)
               CALL sl_from_element(element_str(e, a), xff(:, i), nsl(i), eff(:, i))
               conc(e + ielement_in_cell(a)) = conc(e + ielement_in_cell(1))
               nsl(e + ielement_in_cell(a)) = nsl(e + ielement_in_cell(1))
               xff(:, e + ielement_in_cell(a)) = xff(:, e + ielement_in_cell(1))
               eff(:, e + ielement_in_cell(a)) = eff(:, e + ielement_in_cell(1))
            END DO
         END DO
      ELSE IF (l < num_atoms_in_cell) THEN
         WRITE (*, *) 'Supercell error: Please specify OCC for all sites in "', TRIM(read_file), '".'
         STOP
      END IF
      !if(abs(cell_angles(1)-90.0)<eps.and.abs(cell_angles(2)-90.0)<eps.and.abs(cell_angles(3)-90.0)<eps) then
      ! orthogonal=.true.
      !endif
      IF (num_atoms_in_cell < 1) THEN
         WRITE (*, *) 'Supercell error: Please specify SITE in "', TRIM(read_file), '".'
         STOP
      END IF

      IF (ABS(SUM(cell_params)) < eps .OR. ABS(SUM(cell_angles)) < eps) THEN
         WRITE (*, *) 'Supercell error: Please specify CELL in "', TRIM(read_file), '".'
         WRITE (*, *) 'Format is 3 lattice parameters (a,b,c in Angstrom) + 3 cell angles (alpha,beta,gamma).'
         STOP
      END IF

      DO i = 1, 3
         IF (cell_params(i) < 0D0 .OR. cell_angles(i) < 0D0) THEN
            WRITE (*, *) 'Supercell error: CELL must be positive.'
            STOP
         END IF
      END DO

      tot_num_cells = PRODUCT(num_cells)
      IF (tot_num_cells == 0) THEN
         WRITE (*, *) 'Supercell error: Please specify BOX in "', TRIM(read_file), '".'
         STOP
      END IF

      num_atoms = tot_num_cells*SUM(conc)
      !write(*,*) 'Num atoms:',num_atoms
      IF (ABS(num_atoms) < eps) THEN
         WRITE (*, *) 'Supercell error: It appears OCC = 0 for all sites. Nothing to calculate.'
         STOP
      END IF

   END SUBROUTINE read_config_short

   SUBROUTINE read_spindiff_config

      CHARACTER(128) :: scale_str, flat_bgr_str, linear_bgr_str

      name = 'unspecified'
      flat_bgr_str = 'unspecified'
      linear_bgr_str = 'unspecified'
      scale_str = 'refine'
      colourmap_str = 'default'
      laue_class_str = '-1'
      window_str = 'lanczos'
      sum_str = 'parallel'
      multiplier = 2
      temp_subtract = .FALSE.
      n = 3
      r_max = 0.5 ! default value considers correlations up to half the box size
      max = -1.0
      min = -1.0
      centre = 0D0
      axis = 0D0
      nq = 0
      output_ascii = .TRUE.
      output_vtk = .TRUE.
      output_ppm = .FALSE.
      output_supercell_bragg = .FALSE.
      radiation_str = '?'
      remove_bragg = .FALSE.
      centring_str = 'P'
      centring = 1
      max_taylor = HUGE(1)
      expansion_max_error = HUGE(0D0)
      use_expansion_order = .FALSE.
      use_expansion_max_error = .FALSE.
      mag_only = .FALSE.
      clip_dist = HUGE(0D0)
      clip_displacements = .FALSE.
      max_num_threads = 0
      scale = 0D0
      flat_bgr = 0D0
      linear_bgr = 0D0
      refine_scale = .TRUE.
      refine_flat_bgr = .FALSE.
      refine_linear_bgr = .FALSE.
      hkl2xyz = 0D0
      hkl2xyz = 0D0
      suppress_errors = .FALSE.

      OPEN (11, file='scatty_config.txt', status='old')
      DO
         READ (11, '(a)', END=8) str
         i = INDEX(str, 'PPM_OUTPUT')
         IF (i .NE. 0) output_ppm = .TRUE.
         i = INDEX(str, 'SUPERCELL_BRAGG_OUTPUT')
         IF (i .NE. 0) output_supercell_bragg = .TRUE.
         i = INDEX(str, 'SYMMETRY ')
         IF (i .NE. 0) READ (str(i + 9:), *) laue_class_str
         i = INDEX(str, 'CENTRE ')
         IF (i .NE. 0) READ (str(i + 7:), *) centre(:)
         i = INDEX(str, 'X_AXIS ')
         IF (i .NE. 0) READ (str(i + 7:), *) axis(1, :), nq(1)
         i = INDEX(str, 'Y_AXIS ')
         IF (i .NE. 0) READ (str(i + 7:), *) axis(2, :), nq(2)
         i = INDEX(str, 'Z_AXIS ')
         IF (i .NE. 0) READ (str(i + 7:), *) axis(3, :), nq(3)
         i = INDEX(str, 'HKL_TO_X')
         IF (i .NE. 0) READ (str(i + 8:), *) hkl2xyz(1, :)
         i = INDEX(str, 'HKL_TO_Y')
         IF (i .NE. 0) READ (str(i + 8:), *) hkl2xyz(2, :)
         i = INDEX(str, 'HKL_TO_Z')
         IF (i .NE. 0) READ (str(i + 8:), *) hkl2xyz(3, :)
         i = INDEX(str, 'WINDOW ')
         IF (i .NE. 0) READ (str(i + 7:), *) n
         i = INDEX(str, 'WINDOW_TYPE ')
         IF (i .NE. 0) READ (str(i + 12:), *) window_str
         i = INDEX(str, 'SUM ')
         IF (i .NE. 0) READ (str(i + 4:), *) sum_str
         i = INDEX(str, 'CUTOFF ')
         IF (i .NE. 0) READ (str(i + 7:), *) multiplier
         i = INDEX(str, 'TEMP_SUBTRACT')
         IF (i .NE. 0) temp_subtract = .TRUE.
         i = INDEX(str, 'PPM_RANGE ')
         IF (i .NE. 0) READ (str(i + 10:), *) min, max
         i = INDEX(str, 'PPM_COLOURMAP ')
         IF (i .NE. 0) READ (str(i + 14:), *) colourmap_str
         i = INDEX(str, 'RADIATION ')
         IF (i .NE. 0) READ (str(i + 10:), *) radiation_str
         i = INDEX(str, 'MAG_ONLY')
         IF (i .NE. 0) mag_only = .TRUE.
         i = INDEX(str, 'SUPPRESS_ERRORS')
         IF (i .NE. 0) suppress_errors = .TRUE.
         i = INDEX(str, 'REMOVE_BRAGG ')
         IF (i .NE. 0) THEN
            READ (str(i + 13:), *) centring_str
            remove_bragg = .TRUE.
         END IF
         i = INDEX(str, 'NAME ')
         IF (i .NE. 0) READ (str(i + 5:), *) name
         i = INDEX(str, 'EXPANSION_ORDER ')
         IF (i .NE. 0) THEN
            use_expansion_order = .TRUE.
            READ (str(i + 16:), *) max_taylor
         END IF
         i = INDEX(str, 'EXPANSION_MAX_ERROR ')
         IF (i .NE. 0) THEN
            use_expansion_max_error = .TRUE.
            READ (str(i + 20:), *) expansion_max_error
         END IF
         i = INDEX(str, 'MAX_NUM_THREADS ')
         IF (i .NE. 0) THEN
            READ (str(i + 16:), *) max_num_threads
         END IF
         i = INDEX(str, 'CLIP_DISPLACEMENTS ')
         IF (i .NE. 0) THEN
            clip_displacements = .TRUE.
            READ (str(i + 19:), *) clip_dist
         END IF
         i = INDEX(str, 'SCALE')
         IF (i .NE. 0) THEN
            READ (str(i + 5:), *) scale_str
            IF ((scale_str == 'refine') .OR. (scale_str == 'REFINE')) THEN
               refine_scale = .TRUE.
            ELSE IF (scale_str == 'unspecified') THEN
               WRITE (*, *) 'CONFIG error: Please specify XTAL_SCALE in config file'
               WRITE (*, *) ''
               WRITE (*, *) 'Allowed values are'
               WRITE (*, *) 'either the word REFINE if the intensity scale is to be refined, or'
               WRITE (*, *) 'a real number giving the fixed intensity scale in barns/(st.atom)'
               WRITE (*, *) ''
               STOP
            ELSE
               READ (str(i + 5:), *) scale
               refine_scale = .FALSE.
            END IF
         END IF
         i = INDEX(str, 'FLAT_BACKGROUND')
         IF (i .NE. 0) THEN
            READ (str(i + 15:), *) flat_bgr_str
            IF ((flat_bgr_str == 'refine') .OR. (flat_bgr_str == 'REFINE')) THEN
               refine_flat_bgr = .TRUE.
            ELSE IF (flat_bgr_str == 'unspecified') THEN
               refine_flat_bgr = .FALSE.
               flat_bgr = 0D0
            ELSE
               refine_flat_bgr = .FALSE.
               READ (str(i + 15:), *) flat_bgr
            END IF
         END IF
         i = INDEX(str, 'LINEAR_BACKGROUND')
         IF (i .NE. 0) THEN
            READ (str(i + 17:), *) linear_bgr_str
            IF ((linear_bgr_str == 'refine') .OR. (linear_bgr_str == 'REFINE')) THEN
               refine_linear_bgr = .TRUE.
            ELSE IF (linear_bgr_str == 'unspecified') THEN
               refine_linear_bgr = .FALSE.
               linear_bgr = 0D0
            ELSE
               refine_linear_bgr = .FALSE.
               READ (str(i + 17:), *) linear_bgr
            END IF
         END IF
      END DO
8     CLOSE (11)

      IF (radiation_str == 'n' .OR. radiation_str == 'N') THEN
         xtal_radiation = 1
      ELSE IF (radiation_str == 'x' .OR. radiation_str == 'X') THEN
         xtal_radiation = 2
         IF (mag_only) THEN
            WRITE (*, *) 'Config error: MAG_ONLY is not available for X-ray scattering.'
            STOP
         END IF
         IF (temp_subtract) THEN
            WRITE (*, *) 'Config error: TEMP_SUBTRACT is not available for X-ray scattering.'
            STOP
         END IF
      ELSE IF (radiation_str == 'e' .OR. radiation_str == 'E') THEN
         xtal_radiation = 3
         IF (mag_only) THEN
            WRITE (*, *) 'Config error: MAG_ONLY is not available for electron scattering.'
            STOP
         END IF
         IF (temp_subtract) THEN
            WRITE (*, *) 'Config error: TEMP_SUBTRACT is not available for electron scattering.'
            STOP
         END IF
      ELSE
         WRITE (*, *) 'Config error: Please specify XTAL_RADIATION.'
         WRITE (*, *) 'Options are "N" for neutron, "X" for X-ray, or "E" for electron.'
         STOP
      END IF

      IF (MOD(n*multiplier/2D0, 1D0) > eps) STOP 'Config error: WINDOW/CUTOFF must be a multiple of [integer]/2.'

      ! check for mistakes in program input
      IF (name == 'unspecified_name') THEN
         WRITE (*, *) 'Config error: Please specify NAME.'
         STOP
      END IF

      IF (nq(1) == 0 .AND. nq(2) == 0 .AND. nq(3) == 0) THEN
         WRITE (*, *) 'Config error: Please specify at least one of X_AXIS, Y_AXIS, and Z_AXIS.'
         STOP
      END IF

      IF (max_taylor < 0) THEN
         WRITE (*, *) 'Config error: EXPANSION_ORDER should be => 0.'
         STOP
      END IF

      IF (expansion_max_error < 2D0*EPSILON(0D0)) THEN
         WRITE (*, *) 'Config error: EXPANSION_MAX_ERROR should be > 0.0.'
         STOP
      END IF

      IF (nq(1) < 0 .OR. nq(2) < 0 .OR. nq(3) < 0) THEN
         WRITE (*, *) 'Config error: Number of points for calculation should be > 0.'
         STOP
      END IF

      IF (centring_str == 'P') THEN
         centring = 1
      ELSEIF (centring_str == 'I') THEN
         centring = 2
      ELSEIF (centring_str == 'F') THEN
         centring = 3
      ELSEIF (centring_str == 'R') THEN
         centring = 4
      ELSEIF (centring_str == 'C') THEN
         centring = 5
      ELSEIF (centring_str == 'A') THEN
         centring = 6
      ELSEIF (centring_str == 'B') THEN
         centring = 7
      ELSEIF (centring_str == 'H') THEN
         centring = 8
      ELSE
         WRITE (*, *) 'Config error: CENTRING not recognised. Options are P,I,F,R,C,A,B, or H.'
      END IF

      IF ((colourmap_str == 'default') .OR. (colourmap_str == 'DEFAULT')) THEN
         colourmap = 1
      ELSE IF ((colourmap_str == 'heat') .OR. (colourmap_str == 'HEAT')) THEN
         colourmap = 2
      ELSE IF ((colourmap_str == 'jet') .OR. (colourmap_str == 'JET')) THEN
         colourmap = 3
      ELSE IF ((colourmap_str == 'grey1') .OR. (colourmap_str == 'GREY1')) THEN
         colourmap = 4
      ELSE IF ((colourmap_str == 'grey2') .OR. (colourmap_str == 'GREY2')) THEN
         colourmap = 5
      ELSE IF ((colourmap_str == 'viridis') .OR. (colourmap_str == 'VIRIDIS')) THEN
         colourmap = 6
      ELSE
         WRITE (*, *) 'Config error: COLOURMAP not recognised.'
         STOP
      END IF

      IF (n < 0 .OR. n == 1) THEN
         WRITE (*, *) 'Config error: WINDOW length must be 0, or greater than 1.'
         STOP
      END IF

      IF ((window_str == 'lanczos') .OR. (window_str == 'LANCZOS')) THEN
         window_type = 1
      ELSE IF ((window_str == 'cosine') .OR. (window_str == 'COSINE')) THEN
         window_type = 2
      ELSE
         WRITE (*, *) 'Config error: Please specify WINDOW [n] [type]. Allowed options are: n = integer > 1, type = COS or SINC.'
         STOP
      END IF

      IF ((sum_str == 'parallel') .OR. (window_str == 'PARALLEL')) THEN
         sum_type = 1
      ELSE IF ((sum_str == 'sphere') .OR. (sum_str == 'SPHERE')) THEN
         sum_type = 2
      ELSE
         WRITE (*, *) 'Config error: Please specify a valid SUM. Allowed options are PARALLEL or SPHERE.'
         STOP
      END IF

      IF (laue_class_str == '-1') THEN
         laue_class = 1
      ELSE IF (TRIM(ADJUSTL(laue_class_str)) == '2|m') THEN
         laue_class = 2
      ELSE IF (TRIM(ADJUSTL(laue_class_str)) == 'mmm') THEN
         laue_class = 3
      ELSE IF (TRIM(ADJUSTL(laue_class_str)) == '4|m') THEN
         laue_class = 4
      ELSE IF (TRIM(ADJUSTL(laue_class_str)) == '4|mmm') THEN
         laue_class = 5
      ELSE IF (TRIM(ADJUSTL(laue_class_str)) == '-3') THEN
         laue_class = 6
      ELSE IF (TRIM(ADJUSTL(laue_class_str)) == '-31m') THEN
         laue_class = 7
      ELSE IF (TRIM(ADJUSTL(laue_class_str)) == '-3m1' .OR. TRIM(ADJUSTL(laue_class_str)) == '-3m') THEN
         laue_class = 8
      ELSE IF (TRIM(ADJUSTL(laue_class_str)) == '6|m') THEN
         laue_class = 9
      ELSE IF (TRIM(ADJUSTL(laue_class_str)) == '6|mmm') THEN
         laue_class = 10
      ELSE IF (TRIM(ADJUSTL(laue_class_str)) == 'm-3' .OR. TRIM(ADJUSTL(laue_class_str)) == 'm3') THEN
         laue_class = 11
      ELSE IF (TRIM(ADJUSTL(laue_class_str)) == 'm-3m' .OR. TRIM(ADJUSTL(laue_class_str)) == 'm3m') THEN
         laue_class = 12
      ELSE
         WRITE (*, *) 'Config error: Laue class not recognised.'
         STOP
      END IF

      WRITE (*, *) 'SCATTY configuration file "scatty_config.txt" read successfully!'

   END SUBROUTINE read_spindiff_config

   SUBROUTINE sl_from_element(elem_str, xff, nsl, eff)

      CHARACTER*(*), INTENT(in) :: elem_str ! input
      REAL(8), INTENT(out) :: nsl, xff(9), eff(10)
      CHARACTER(len=4), DIMENSION(212) :: element_table
      INTEGER :: i
      REAL(8) :: nsl_table(212), xff_table(9, 212), eff_table(10, 212)
      LOGICAL :: found_element

      element_table(1) = "H   "
      element_table(2) = "H.  "
      element_table(3) = "D   "
      element_table(4) = "H1- "
      element_table(5) = "He  "
      element_table(6) = "Li  "
      element_table(7) = "Li1+"
      element_table(8) = "Be  "
      element_table(9) = "Be2+"
      element_table(10) = "B   "
      element_table(11) = "C   "
      element_table(12) = "C.  "
      element_table(13) = "N   "
      element_table(14) = "O   "
      element_table(15) = "O1- "
      element_table(16) = "O2- "
      element_table(17) = "F   "
      element_table(18) = "F1- "
      element_table(19) = "Ne  "
      element_table(20) = "Na  "
      element_table(21) = "Na1+"
      element_table(22) = "Mg  "
      element_table(23) = "Mg2+"
      element_table(24) = "Al  "
      element_table(25) = "Al3+"
      element_table(26) = "Si  "
      element_table(27) = "Si. "
      element_table(28) = "Si4+"
      element_table(29) = "P   "
      element_table(30) = "S   "
      element_table(31) = "Cl  "
      element_table(32) = "Cl1-"
      element_table(33) = "Ar  "
      element_table(34) = "K   "
      element_table(35) = "K1+ "
      element_table(36) = "Ca  "
      element_table(37) = "Ca2+"
      element_table(38) = "Sc  "
      element_table(39) = "Sc3+"
      element_table(40) = "Ti  "
      element_table(41) = "Ti2+"
      element_table(42) = "Ti3+"
      element_table(43) = "Ti4+"
      element_table(44) = "V   "
      element_table(45) = "V2+ "
      element_table(46) = "V3+ "
      element_table(47) = "V5+ "
      element_table(48) = "Cr  "
      element_table(49) = "Cr2+"
      element_table(50) = "Cr3+"
      element_table(51) = "Mn  "
      element_table(52) = "Mn2+"
      element_table(53) = "Mn3+"
      element_table(54) = "Mn4+"
      element_table(55) = "Fe  "
      element_table(56) = "Fe2+"
      element_table(57) = "Fe3+"
      element_table(58) = "Co  "
      element_table(59) = "Co2+"
      element_table(60) = "Co3+"
      element_table(61) = "Ni  "
      element_table(62) = "Ni2+"
      element_table(63) = "Ni3+"
      element_table(64) = "Cu  "
      element_table(65) = "Cu1+"
      element_table(66) = "Cu2+"
      element_table(67) = "Zn  "
      element_table(68) = "Zn2+"
      element_table(69) = "Ga  "
      element_table(70) = "Ga3+"
      element_table(71) = "Ge  "
      element_table(72) = "Ge4+"
      element_table(73) = "As  "
      element_table(74) = "Se  "
      element_table(75) = "Br  "
      element_table(76) = "Br1-"
      element_table(77) = "Kr  "
      element_table(78) = "Rb  "
      element_table(79) = "Rb1+"
      element_table(80) = "Sr  "
      element_table(81) = "Sr2+"
      element_table(82) = "Y   "
      element_table(83) = "Y3+ "
      element_table(84) = "Zr  "
      element_table(85) = "Zr4+"
      element_table(86) = "Nb  "
      element_table(87) = "Nb3+"
      element_table(88) = "Nb5+"
      element_table(89) = "Mo  "
      element_table(90) = "Mo3+"
      element_table(91) = "Mo5+"
      element_table(92) = "Mo6+"
      element_table(93) = "Tc  "
      element_table(94) = "Ru  "
      element_table(95) = "Ru3+"
      element_table(96) = "Ru4+"
      element_table(97) = "Rh  "
      element_table(98) = "Rh3+"
      element_table(99) = "Rh4+"
      element_table(100) = "Pd  "
      element_table(101) = "Pd2+"
      element_table(102) = "Pd4+"
      element_table(103) = "Ag  "
      element_table(104) = "Ag1+"
      element_table(105) = "Ag2+"
      element_table(106) = "Cd  "
      element_table(107) = "Cd2+"
      element_table(108) = "In  "
      element_table(109) = "In3+"
      element_table(110) = "Sn  "
      element_table(111) = "Sn2+"
      element_table(112) = "Sn4+"
      element_table(113) = "Sb  "
      element_table(114) = "Sb2+"
      element_table(115) = "Sb5+"
      element_table(116) = "Te  "
      element_table(117) = "I   "
      element_table(118) = "I1- "
      element_table(119) = "Xe  "
      element_table(120) = "Cs  "
      element_table(121) = "Cs1+"
      element_table(122) = "Ba  "
      element_table(123) = "Ba2+"
      element_table(124) = "La  "
      element_table(125) = "La3+"
      element_table(126) = "Ce  "
      element_table(127) = "Ce3+"
      element_table(128) = "Ce4+"
      element_table(129) = "Pr  "
      element_table(130) = "Pr3+"
      element_table(131) = "Pr4+"
      element_table(132) = "Nd  "
      element_table(133) = "Nd3+"
      element_table(134) = "Pm  "
      element_table(135) = "Pm3+"
      element_table(136) = "Sm  "
      element_table(137) = "Sm3+"
      element_table(138) = "Eu  "
      element_table(139) = "Eu2+"
      element_table(140) = "Eu3+"
      element_table(141) = "Gd  "
      element_table(142) = "Gd3+"
      element_table(143) = "Tb  "
      element_table(144) = "Tb3+"
      element_table(145) = "Dy  "
      element_table(146) = "Dy3+"
      element_table(147) = "Ho  "
      element_table(148) = "Ho3+"
      element_table(149) = "Er  "
      element_table(150) = "Er3+"
      element_table(151) = "Tm  "
      element_table(152) = "Tm3+"
      element_table(153) = "Yb  "
      element_table(154) = "Yb2+"
      element_table(155) = "Yb3+"
      element_table(156) = "Lu  "
      element_table(157) = "Lu3+"
      element_table(158) = "Hf  "
      element_table(159) = "Hf4+"
      element_table(160) = "Ta  "
      element_table(161) = "Ta5+"
      element_table(162) = "W   "
      element_table(163) = "W6+ "
      element_table(164) = "Re  "
      element_table(165) = "Os  "
      element_table(166) = "Os4+"
      element_table(167) = "Ir  "
      element_table(168) = "Ir3+"
      element_table(169) = "Ir4+"
      element_table(170) = "Pt  "
      element_table(171) = "Pt2+"
      element_table(172) = "Pt4+"
      element_table(173) = "Au  "
      element_table(174) = "Au1+"
      element_table(175) = "Au3+"
      element_table(176) = "Hg  "
      element_table(177) = "Hg1+"
      element_table(178) = "Hg2+"
      element_table(179) = "Tl  "
      element_table(180) = "Tl1+"
      element_table(181) = "Tl3+"
      element_table(182) = "Pb  "
      element_table(183) = "Pb2+"
      element_table(184) = "Pb4+"
      element_table(185) = "Bi  "
      element_table(186) = "Bi3+"
      element_table(187) = "Bi5+"
      element_table(188) = "Po  "
      element_table(189) = "At  "
      element_table(190) = "Rn  "
      element_table(191) = "Fr  "
      element_table(192) = "Ra  "
      element_table(193) = "Ra2+"
      element_table(194) = "Ac  "
      element_table(195) = "Ac3+"
      element_table(196) = "Th  "
      element_table(197) = "Th4+"
      element_table(198) = "Pa  "
      element_table(199) = "U   "
      element_table(200) = "U3+ "
      element_table(201) = "U4+ "
      element_table(202) = "U6+ "
      element_table(203) = "Np  "
      element_table(204) = "Np3+"
      element_table(205) = "Np4+"
      element_table(206) = "Np6+"
      element_table(207) = "Pu  "
      element_table(208) = "Pu3+"
      element_table(209) = "Pu4+"
      element_table(210) = "Pu6+"
      element_table(211) = "Am  "
      element_table(212) = "Cm  "

      nsl_table(1) = -0.374100D0
      nsl_table(2) = -0.374100D0
      nsl_table(3) = 0.667400D0
      nsl_table(4) = -0.374100D0
      nsl_table(5) = 0.326000D0
      nsl_table(6) = -0.203000D0
      nsl_table(7) = -0.203000D0
      nsl_table(8) = 0.779000D0
      nsl_table(9) = 0.779000D0
      nsl_table(10) = 0.535000D0
      nsl_table(11) = 0.664800D0
      nsl_table(12) = 0.664800D0
      nsl_table(13) = 0.930000D0
      nsl_table(14) = 0.580500D0
      nsl_table(15) = 0.580500D0
      nsl_table(16) = 0.580500D0
      nsl_table(17) = 0.565000D0
      nsl_table(18) = 0.565000D0
      nsl_table(19) = 0.455000D0
      nsl_table(20) = 0.363000D0
      nsl_table(21) = 0.363000D0
      nsl_table(22) = 0.537500D0
      nsl_table(23) = 0.537500D0
      nsl_table(24) = 0.344900D0
      nsl_table(25) = 0.344900D0
      nsl_table(26) = 0.414900D0
      nsl_table(27) = 0.414900D0
      nsl_table(28) = 0.414900D0
      nsl_table(29) = 0.513000D0
      nsl_table(30) = 0.284700D0
      nsl_table(31) = 0.957900D0
      nsl_table(32) = 0.957900D0
      nsl_table(33) = 0.188400D0
      nsl_table(34) = 0.367000D0
      nsl_table(35) = 0.367000D0
      nsl_table(36) = 0.490000D0
      nsl_table(37) = 0.490000D0
      nsl_table(38) = 1.230000D0
      nsl_table(39) = 1.230000D0
      nsl_table(40) = -0.343800D0
      nsl_table(41) = -0.343800D0
      nsl_table(42) = -0.343800D0
      nsl_table(43) = -0.343800D0
      nsl_table(44) = -0.038200D0
      nsl_table(45) = -0.038200D0
      nsl_table(46) = -0.038200D0
      nsl_table(47) = -0.038200D0
      nsl_table(48) = 0.363500D0
      nsl_table(49) = 0.363500D0
      nsl_table(50) = 0.363500D0
      nsl_table(51) = -0.373000D0
      nsl_table(52) = -0.373000D0
      nsl_table(53) = -0.373000D0
      nsl_table(54) = -0.373000D0
      nsl_table(55) = 0.954000D0
      nsl_table(56) = 0.954000D0
      nsl_table(57) = 0.954000D0
      nsl_table(58) = 0.253000D0
      nsl_table(59) = 0.253000D0
      nsl_table(60) = 0.253000D0
      nsl_table(61) = 1.030000D0
      nsl_table(62) = 1.030000D0
      nsl_table(63) = 1.030000D0
      nsl_table(64) = 0.771800D0
      nsl_table(65) = 0.771800D0
      nsl_table(66) = 0.771800D0
      nsl_table(67) = 0.568000D0
      nsl_table(68) = 0.568000D0
      nsl_table(69) = 0.729000D0
      nsl_table(70) = 0.729000D0
      nsl_table(71) = 0.819300D0
      nsl_table(72) = 0.819300D0
      nsl_table(73) = 0.658000D0
      nsl_table(74) = 0.797000D0
      nsl_table(75) = 0.679000D0
      nsl_table(76) = 0.679000D0
      nsl_table(77) = 0.785000D0
      nsl_table(78) = 0.708000D0
      nsl_table(79) = 0.708000D0
      nsl_table(80) = 0.702000D0
      nsl_table(81) = 0.702000D0
      nsl_table(82) = 0.775000D0
      nsl_table(83) = 0.775000D0
      nsl_table(84) = 0.716000D0
      nsl_table(85) = 0.716000D0
      nsl_table(86) = 0.705400D0
      nsl_table(87) = 0.705400D0
      nsl_table(88) = 0.705400D0
      nsl_table(89) = 0.695000D0
      nsl_table(90) = 0.695000D0
      nsl_table(91) = 0.695000D0
      nsl_table(92) = 0.695000D0
      nsl_table(93) = 0.680000D0
      nsl_table(94) = 0.721000D0
      nsl_table(95) = 0.721000D0
      nsl_table(96) = 0.721000D0
      nsl_table(97) = 0.593000D0
      nsl_table(98) = 0.593000D0
      nsl_table(99) = 0.593000D0
      nsl_table(100) = 0.591000D0
      nsl_table(101) = 0.591000D0
      nsl_table(102) = 0.591000D0
      nsl_table(103) = 0.597000D0
      nsl_table(104) = 0.597000D0
      nsl_table(105) = 0.597000D0
      nsl_table(106) = 0.500000D0
      nsl_table(107) = 0.500000D0
      nsl_table(108) = 0.406000D0
      nsl_table(109) = 0.406000D0
      nsl_table(110) = 0.622800D0
      nsl_table(111) = 0.622800D0
      nsl_table(112) = 0.622800D0
      nsl_table(113) = 0.564000D0
      nsl_table(114) = 0.564000D0
      nsl_table(115) = 0.564000D0
      nsl_table(116) = 0.580000D0
      nsl_table(117) = 0.528000D0
      nsl_table(118) = 0.528000D0
      nsl_table(119) = 0.489000D0
      nsl_table(120) = 0.542000D0
      nsl_table(121) = 0.542000D0
      nsl_table(122) = 0.525000D0
      nsl_table(123) = 0.525000D0
      nsl_table(124) = 0.827000D0
      nsl_table(125) = 0.827000D0
      nsl_table(126) = 0.484000D0
      nsl_table(127) = 0.484000D0
      nsl_table(128) = 0.484000D0
      nsl_table(129) = 0.445000D0
      nsl_table(130) = 0.445000D0
      nsl_table(131) = 0.445000D0
      nsl_table(132) = 0.769000D0
      nsl_table(133) = 0.769000D0
      nsl_table(134) = 1.260000D0
      nsl_table(135) = 1.260000D0
      nsl_table(136) = 0.800000D0
      nsl_table(137) = 0.800000D0
      nsl_table(138) = 0.600000D0
      nsl_table(139) = 0.600000D0
      nsl_table(140) = 0.600000D0
      nsl_table(141) = 0.950000D0
      nsl_table(142) = 0.950000D0
      nsl_table(143) = 0.738000D0
      nsl_table(144) = 0.738000D0
      nsl_table(145) = 1.690000D0
      nsl_table(146) = 1.690000D0
      nsl_table(147) = 0.808000D0
      nsl_table(148) = 0.808000D0
      nsl_table(149) = 0.803000D0
      nsl_table(150) = 0.803000D0
      nsl_table(151) = 0.705000D0
      nsl_table(152) = 0.705000D0
      nsl_table(153) = 1.240000D0
      nsl_table(154) = 1.240000D0
      nsl_table(155) = 1.240000D0
      nsl_table(156) = 0.730000D0
      nsl_table(157) = 0.730000D0
      nsl_table(158) = 0.770000D0
      nsl_table(159) = 0.770000D0
      nsl_table(160) = 0.691000D0
      nsl_table(161) = 0.691000D0
      nsl_table(162) = 0.477000D0
      nsl_table(163) = 0.477000D0
      nsl_table(164) = 0.920000D0
      nsl_table(165) = 1.070000D0
      nsl_table(166) = 1.070000D0
      nsl_table(167) = 1.060000D0
      nsl_table(168) = 1.060000D0
      nsl_table(169) = 1.060000D0
      nsl_table(170) = 0.950000D0
      nsl_table(171) = 0.950000D0
      nsl_table(172) = 0.950000D0
      nsl_table(173) = 0.763000D0
      nsl_table(174) = 0.763000D0
      nsl_table(175) = 0.763000D0
      nsl_table(176) = 1.266000D0
      nsl_table(177) = 1.266000D0
      nsl_table(178) = 1.266000D0
      nsl_table(179) = 0.879000D0
      nsl_table(180) = 0.879000D0
      nsl_table(181) = 0.879000D0
      nsl_table(182) = 0.940100D0
      nsl_table(183) = 0.940100D0
      nsl_table(184) = 0.940100D0
      nsl_table(185) = 0.853300D0
      nsl_table(186) = 0.853300D0
      nsl_table(187) = 0.853300D0
      nsl_table(188) = 0.000000D0
      nsl_table(189) = 0.000000D0
      nsl_table(190) = 0.000000D0
      nsl_table(191) = 0.000000D0
      nsl_table(192) = 1.000000D0
      nsl_table(193) = 1.000000D0
      nsl_table(194) = 0.000000D0
      nsl_table(195) = 0.000000D0
      nsl_table(196) = 0.984000D0
      nsl_table(197) = 0.984000D0
      nsl_table(198) = 0.910000D0
      nsl_table(199) = 0.842000D0
      nsl_table(200) = 0.842000D0
      nsl_table(201) = 0.842000D0
      nsl_table(202) = 0.842000D0
      nsl_table(203) = 1.060000D0
      nsl_table(204) = 1.060000D0
      nsl_table(205) = 1.060000D0
      nsl_table(206) = 1.060000D0
      nsl_table(207) = 0.750000D0
      nsl_table(208) = 0.750000D0
      nsl_table(209) = 0.750000D0
      nsl_table(210) = 0.750000D0
      nsl_table(211) = 0.830000D0
      nsl_table(212) = 0.950000D0

      xff_table(:, 1) = [0.493002, 10.510900, 0.322912, 26.125700, 0.140191, 3.142360, 0.040810, 57.799702, 0.003038]
      xff_table(:, 2) = [0.489918, 20.659300, 0.262003, 7.740390, 0.196767, 49.551899, 0.049879, 2.201590, 0.001305]
      xff_table(:, 3) = [0.489918, 20.659300, 0.262003, 7.740390, 0.196767, 49.551899, 0.049879, 2.201590, 0.001305]
      xff_table(:, 4) = [0.897661, 53.136799, 0.565616, 15.187000, 0.415815, 186.576004, 0.116973, 3.567090, 0.002389]
      xff_table(:, 5) = [0.873400, 9.103700, 0.630900, 3.356800, 0.311200, 22.927601, 0.178000, 0.982100, 0.006400]
      xff_table(:, 6) = [1.128200, 3.954600, 0.750800, 1.052400, 0.617500, 85.390503, 0.465300, 168.261002, 0.037700]
      xff_table(:, 7) = [0.696800, 4.623700, 0.788800, 1.955700, 0.341400, 0.631600, 0.156300, 10.095300, 0.016700]
      xff_table(:, 8) = [1.591900, 43.642700, 1.127800, 1.862300, 0.539100, 103.483002, 0.702900, 0.542000, 0.038500]
      xff_table(:, 9) = [6.260300, 0.002700, 0.884900, 0.831300, 0.799300, 2.275800, 0.164700, 5.114600, -6.109200]
      xff_table(:, 10) = [2.054500, 23.218500, 1.332600, 1.021000, 1.097900, 60.349800, 0.706800, 0.140300, -0.193200]
      xff_table(:, 11) = [2.310000, 20.843901, 1.020000, 10.207500, 1.588600, 0.568700, 0.865000, 51.651199, 0.215600]
      xff_table(:, 12) = [2.260690, 22.690701, 1.561650, 0.656665, 1.050750, 9.756180, 0.839259, 55.594898, 0.286977]
      xff_table(:, 13) = [12.212600, 0.005700, 3.132200, 9.893300, 2.012500, 28.997499, 1.166300, 0.582600, -11.529000]
      xff_table(:, 14) = [3.048500, 13.277100, 2.286800, 5.701100, 1.546300, 0.323900, 0.867000, 32.908901, 0.250800]
      xff_table(:, 15) = [4.191600, 12.857300, 1.639690, 4.172360, 1.526730, 47.017899, -20.306999, -0.014040, 21.941200]
      xff_table(:, 16) = [4.191600, 12.857300, 1.639690, 4.172360, 1.526730, 47.017899, -20.306999, -0.014040, 21.941200]
      xff_table(:, 17) = [3.539200, 10.282500, 2.641200, 4.294400, 1.517000, 0.261500, 1.024300, 26.147600, 0.277600]
      xff_table(:, 18) = [3.632200, 5.277560, 3.510570, 14.735300, 1.260640, 0.442258, 0.940706, 47.343700, 0.653396]
      xff_table(:, 19) = [3.955300, 8.404200, 3.112500, 3.426200, 1.454600, 0.230600, 1.125100, 21.718399, 0.351500]
      xff_table(:, 20) = [4.762600, 3.285000, 3.173600, 8.842200, 1.267400, 0.313600, 1.112800, 129.423996, 0.676000]
      xff_table(:, 21) = [3.256500, 2.667100, 3.936200, 6.115300, 1.399800, 0.200100, 1.003200, 14.039000, 0.404000]
      xff_table(:, 22) = [5.420400, 2.827500, 2.173500, 79.261101, 1.226900, 0.380800, 2.307300, 7.193700, 0.858400]
      xff_table(:, 23) = [3.498800, 2.167600, 3.837800, 4.754200, 1.328400, 0.185000, 0.849700, 10.141100, 0.485300]
      xff_table(:, 24) = [6.420200, 3.038700, 1.900200, 0.742600, 1.593600, 31.547199, 1.964600, 85.088600, 1.115100]
      xff_table(:, 25) = [4.174480, 1.938160, 3.387600, 4.145530, 1.202960, 0.228753, 0.528137, 8.285240, 0.706786]
      xff_table(:, 26) = [6.291500, 2.438600, 3.035300, 32.333698, 1.989100, 0.678500, 1.541000, 81.693703, 1.140700]
      xff_table(:, 27) = [5.662690, 2.665200, 3.071640, 38.663399, 2.624460, 0.916946, 1.393200, 93.545799, 1.247070]
      xff_table(:, 28) = [4.439180, 1.641670, 3.203450, 3.437570, 1.194530, 0.214900, 0.416530, 6.653650, 0.746297]
      xff_table(:, 29) = [6.434500, 1.906700, 4.179100, 27.157000, 1.780000, 0.526000, 1.490800, 68.164497, 1.114900]
      xff_table(:, 30) = [6.905300, 1.467900, 5.203400, 22.215099, 1.437900, 0.253600, 1.586300, 56.172001, 0.866900]
      xff_table(:, 31) = [11.460400, 0.010400, 7.196400, 1.166200, 6.255600, 18.519400, 1.645500, 47.778400, -9.557400]
      xff_table(:, 32) = [18.291500, 0.006600, 7.208400, 1.171700, 6.533700, 19.542400, 2.338600, 60.448601, -16.378000]
      xff_table(:, 33) = [7.484500, 0.907200, 6.772300, 14.840700, 0.653900, 43.898300, 1.644200, 33.392899, 1.444500]
      xff_table(:, 34) = [8.218600, 12.794900, 7.439800, 0.774800, 1.051900, 213.186996, 0.865900, 41.684101, 1.422800]
      xff_table(:, 35) = [7.957800, 12.633100, 7.491700, 0.767400, 6.359000, -0.002000, 1.191500, 31.912800, -4.997800]
      xff_table(:, 36) = [8.626600, 10.442100, 7.387300, 0.659900, 1.589900, 85.748398, 1.021100, 178.436996, 1.375100]
      xff_table(:, 37) = [15.634800, -0.007400, 7.951800, 0.608900, 8.437200, 10.311600, 0.853700, 25.990499, -14.875000]
      xff_table(:, 38) = [9.189000, 9.021300, 7.367900, 0.572900, 1.640900, 136.108002, 1.468000, 51.353100, 1.332900]
      xff_table(:, 39) = [13.400800, 0.298540, 8.027300, 7.962900, 1.659430, -0.286040, 1.579360, 16.066200, -6.666700]
      xff_table(:, 40) = [9.759500, 7.850800, 7.355800, 0.500000, 1.699100, 35.633801, 1.902100, 116.105003, 1.280700]
      xff_table(:, 41) = [9.114230, 7.524300, 7.621740, 0.457585, 2.279300, 19.536100, 0.087899, 61.655800, 0.897155]
      xff_table(:, 42) = [17.734400, 0.220610, 8.738160, 7.047160, 5.256910, -0.157620, 1.921340, 15.976800, -14.652000]
      xff_table(:, 43) = [19.511400, 0.178847, 8.234730, 6.670180, 2.013410, -0.292630, 1.520800, 12.946400, -13.280000]
      xff_table(:, 44) = [10.297100, 6.865700, 7.351100, 0.438500, 2.070300, 26.893801, 2.057100, 102.477997, 1.219900]
      xff_table(:, 45) = [10.106000, 6.881800, 7.354100, 0.440900, 2.288400, 20.300400, 0.022300, 115.122002, 1.229800]
      xff_table(:, 46) = [9.431410, 6.395350, 7.741900, 0.383349, 2.153430, 15.190800, 0.016865, 63.969002, 0.656565]
      xff_table(:, 47) = [15.688700, 0.679003, 8.142080, 5.401350, 2.030810, 9.972780, -9.576000, 0.940464, 1.714300]
      xff_table(:, 48) = [10.640600, 6.103800, 7.353700, 0.392000, 3.324000, 20.262600, 1.492200, 98.739899, 1.183200]
      xff_table(:, 49) = [9.540340, 5.660780, 7.750900, 0.344261, 3.582740, 13.307500, 0.509107, 32.422401, 0.616898]
      xff_table(:, 50) = [9.680900, 5.594630, 7.811360, 0.334393, 2.876030, 12.828800, 0.113575, 32.876099, 0.518275]
      xff_table(:, 51) = [11.281900, 5.340900, 7.357300, 0.343200, 3.019300, 17.867399, 2.244100, 83.754303, 1.089600]
      xff_table(:, 52) = [10.806100, 5.279600, 7.362000, 0.343500, 3.526800, 14.343000, 0.218400, 41.323502, 1.087400]
      xff_table(:, 53) = [9.845210, 4.917970, 7.871940, 0.294393, 3.565310, 10.817100, 0.323613, 24.128099, 0.393974]
      xff_table(:, 54) = [9.962530, 4.848500, 7.970570, 0.283303, 2.760670, 10.485200, 0.054447, 27.573000, 0.251877]
      xff_table(:, 55) = [11.769500, 4.761100, 7.357300, 0.307200, 3.522200, 15.353500, 2.304500, 76.880501, 1.036900]
      xff_table(:, 56) = [11.042400, 4.653800, 7.374000, 0.305300, 4.134600, 12.054600, 0.439900, 31.280899, 1.009700]
      xff_table(:, 57) = [11.176400, 4.614700, 7.386300, 0.300500, 3.394800, 11.672900, 0.072400, 38.556599, 0.970700]
      xff_table(:, 58) = [12.284100, 4.279100, 7.340900, 0.278400, 4.003400, 13.535900, 2.348800, 71.169197, 1.011800]
      xff_table(:, 59) = [11.229600, 4.123100, 7.388300, 0.272600, 4.739300, 10.244300, 0.710800, 25.646601, 0.932400]
      xff_table(:, 60) = [10.338000, 3.909690, 7.881730, 0.238668, 4.767950, 8.355830, 0.725591, 18.349100, 0.286667]
      xff_table(:, 61) = [12.837600, 3.878500, 7.292000, 0.256500, 4.443800, 12.176300, 2.380000, 66.342102, 1.034100]
      xff_table(:, 62) = [11.416600, 3.676600, 7.400500, 0.244900, 5.344200, 8.873000, 0.977300, 22.162600, 0.861400]
      xff_table(:, 63) = [10.780600, 3.547700, 7.758680, 0.223140, 5.227460, 7.644680, 0.847114, 16.967300, 0.386044]
      xff_table(:, 64) = [13.338000, 3.582800, 7.167600, 0.247000, 5.615800, 11.396600, 1.673500, 64.812599, 1.191000]
      xff_table(:, 65) = [11.947500, 3.366900, 7.357300, 0.227400, 6.245500, 8.662500, 1.557800, 25.848700, 0.890000]
      xff_table(:, 66) = [11.816800, 3.374840, 7.111810, 0.244078, 5.781350, 7.987600, 1.145230, 19.896999, 1.144310]
      xff_table(:, 67) = [14.074300, 3.265500, 7.031800, 0.233300, 5.165200, 10.316300, 2.410000, 58.709702, 1.304100]
      xff_table(:, 68) = [11.971900, 2.994600, 7.386200, 0.203100, 6.466800, 7.082600, 1.394000, 18.099501, 0.780700]
      xff_table(:, 69) = [15.235400, 3.066900, 6.700600, 0.241200, 4.359100, 10.780500, 2.962300, 61.413502, 1.718900]
      xff_table(:, 70) = [12.692000, 2.812620, 6.698830, 0.227890, 6.066920, 6.364410, 1.006600, 14.412200, 1.535450]
      xff_table(:, 71) = [16.081600, 2.850900, 6.374700, 0.251600, 3.706800, 11.446800, 3.683000, 54.762501, 2.131300]
      xff_table(:, 72) = [12.917200, 2.537180, 6.700030, 0.205855, 6.067910, 5.479130, 0.859041, 11.603000, 1.455720]
      xff_table(:, 73) = [16.672300, 2.634500, 6.070100, 0.264700, 3.431300, 12.947900, 4.277900, 47.797199, 2.531000]
      xff_table(:, 74) = [17.000601, 2.409800, 5.819600, 0.272600, 3.973100, 15.237200, 4.354300, 43.816299, 2.840900]
      xff_table(:, 75) = [17.178900, 2.172300, 5.235800, 16.579599, 5.637700, 0.260900, 3.985100, 41.432800, 2.955700]
      xff_table(:, 76) = [17.171801, 2.205900, 6.333800, 19.334499, 5.575400, 0.287100, 3.727200, 58.153500, 3.177600]
      xff_table(:, 77) = [17.355499, 1.938400, 6.728600, 16.562300, 5.549300, 0.226100, 3.537500, 39.397202, 2.825000]
      xff_table(:, 78) = [17.178400, 1.788800, 9.643500, 17.315100, 5.139900, 0.274800, 1.529200, 164.934006, 3.487300]
      xff_table(:, 79) = [17.581600, 1.713900, 7.659800, 14.795700, 5.898100, 0.160300, 2.781700, 31.208700, 2.078200]
      xff_table(:, 80) = [17.566299, 1.556400, 9.818400, 14.098800, 5.422000, 0.166400, 2.669400, 132.376007, 2.506400]
      xff_table(:, 81) = [18.087400, 1.490700, 8.137300, 12.696300, 2.565400, 24.565100, -34.193001, -0.013800, 41.402500]
      xff_table(:, 82) = [17.775999, 1.402900, 10.294600, 12.800600, 5.726290, 0.125599, 3.265880, 104.353996, 1.912130]
      xff_table(:, 83) = [17.926800, 1.354170, 9.153100, 11.214500, 1.767950, 22.659901, -33.108002, -0.013190, 40.260201]
      xff_table(:, 84) = [17.876499, 1.276180, 10.948000, 11.916000, 5.417320, 0.117622, 3.657210, 87.662697, 2.069290]
      xff_table(:, 85) = [18.166800, 1.214800, 10.056200, 10.148300, 1.011180, 21.605400, -2.647900, -0.102760, 9.414540]
      xff_table(:, 86) = [17.614201, 1.188650, 12.014400, 11.766000, 4.041830, 0.204785, 3.533460, 69.795700, 3.755910]
      xff_table(:, 87) = [19.881201, 0.019175, 18.065300, 1.133050, 11.017700, 10.162100, 1.947150, 28.338900, -12.912000]
      xff_table(:, 88) = [17.916300, 1.124460, 13.341700, 0.028781, 10.799000, 9.282060, 0.337905, 25.722799, -6.393400]
      xff_table(:, 89) = [3.702500, 0.277200, 17.235600, 1.095800, 12.887600, 11.004000, 3.742900, 61.658401, 4.387500]
      xff_table(:, 90) = [21.166401, 0.014734, 18.201700, 1.030310, 11.742300, 9.536590, 2.309510, 26.630699, -14.421000]
      xff_table(:, 91) = [21.014900, 0.014345, 18.099199, 1.022380, 11.463200, 8.788090, 0.740625, 23.345200, -14.316000]
      xff_table(:, 92) = [17.887100, 1.036490, 11.175000, 8.480610, 6.578910, 0.058881, 0.000000, 0.000000, 0.344941]
      xff_table(:, 93) = [19.130100, 0.864132, 11.094800, 8.144870, 4.649010, 21.570700, 2.712630, 86.847198, 5.404280]
      xff_table(:, 94) = [19.267401, 0.808520, 12.918200, 8.434670, 4.863370, 24.799700, 1.567560, 94.292801, 5.378740]
      xff_table(:, 95) = [18.563801, 0.847329, 13.288500, 8.371640, 9.326020, 0.017662, 3.009640, 22.886999, -3.189200]
      xff_table(:, 96) = [18.500299, 0.844582, 13.178700, 8.125340, 4.713040, 0.036495, 2.185350, 20.850401, 1.423570]
      xff_table(:, 97) = [19.295700, 0.751536, 14.350100, 8.217580, 4.734250, 25.874901, 1.289180, 98.606201, 5.328000]
      xff_table(:, 98) = [18.878500, 0.764252, 14.125900, 7.844380, 3.325150, 21.248699, -6.198900, -0.010360, 11.867800]
      xff_table(:, 99) = [18.854500, 0.760825, 13.980600, 7.624360, 2.534640, 19.331699, -5.652600, -0.010200, 11.283500]
      xff_table(:, 100) = [19.331900, 0.698655, 15.501700, 7.989290, 5.295370, 25.205200, 0.605844, 76.898598, 5.265930]
      xff_table(:, 101) = [19.170099, 0.696219, 15.209600, 7.555730, 4.322340, 22.505699, 0.000000, 0.000000, 5.291600]
      xff_table(:, 102) = [19.249300, 0.683839, 14.790000, 7.148330, 2.892890, 17.914400, -7.949200, 0.005127, 13.017400]
      xff_table(:, 103) = [19.280800, 0.644600, 16.688499, 7.472600, 4.804500, 24.660500, 1.046300, 99.815598, 5.179000]
      xff_table(:, 104) = [19.181200, 0.646179, 15.971900, 7.191230, 5.274750, 21.732599, 0.357534, 66.114700, 5.215720]
      xff_table(:, 105) = [19.164301, 0.645643, 16.245600, 7.185440, 4.370900, 21.407200, 0.000000, 0.000000, 5.214040]
      xff_table(:, 106) = [19.221399, 0.594600, 17.644400, 6.908900, 4.461000, 24.700800, 1.602900, 87.482498, 5.069400]
      xff_table(:, 107) = [19.151400, 0.597922, 17.253500, 6.806390, 4.471280, 20.252100, 0.000000, 0.000000, 5.119370]
      xff_table(:, 108) = [19.162399, 0.547600, 18.559601, 6.377600, 4.294800, 25.849899, 2.039600, 92.802902, 4.939100]
      xff_table(:, 109) = [19.104500, 0.551522, 18.110800, 6.324700, 3.788970, 17.359501, 0.000000, 0.000000, 4.996350]
      xff_table(:, 110) = [19.188900, 5.830300, 19.100500, 0.503100, 4.458500, 26.890900, 2.466300, 83.957100, 4.782100]
      xff_table(:, 111) = [19.109400, 0.503600, 19.054800, 5.837800, 4.564800, 23.375200, 0.487000, 62.206100, 4.786100]
      xff_table(:, 112) = [18.933300, 5.764000, 19.713100, 0.465500, 3.418200, 14.004900, 0.019300, -0.758300, 3.918200]
      xff_table(:, 113) = [19.641800, 5.303400, 19.045500, 0.460700, 5.037100, 27.907400, 2.682700, 75.282501, 4.590900]
      xff_table(:, 114) = [18.975500, 0.467196, 18.933001, 5.221260, 5.107890, 19.590200, 0.288753, 55.511299, 4.696260]
      xff_table(:, 115) = [19.868500, 5.448530, 19.030199, 0.467973, 2.412530, 14.125900, 0.000000, 0.000000, 4.692630]
      xff_table(:, 116) = [19.964399, 4.817420, 19.013800, 0.420885, 6.144870, 28.528400, 2.523900, 70.840302, 4.352000]
      xff_table(:, 117) = [20.147200, 4.347000, 18.994900, 0.381400, 7.513800, 27.766001, 2.273500, 66.877602, 4.071200]
      xff_table(:, 118) = [20.233200, 4.357900, 18.997000, 0.381500, 7.806900, 29.525900, 2.886800, 84.930397, 4.071400]
      xff_table(:, 119) = [20.293301, 3.928200, 19.029800, 0.344000, 8.976700, 26.465900, 1.990000, 64.265800, 3.711800]
      xff_table(:, 120) = [20.389200, 3.569000, 19.106199, 0.310700, 10.662000, 24.387899, 1.495300, 213.904007, 3.335200]
      xff_table(:, 121) = [20.352400, 3.552000, 19.127800, 0.308600, 10.282100, 23.712799, 0.961500, 59.456501, 3.279100]
      xff_table(:, 122) = [20.336100, 3.216000, 19.297001, 0.275600, 10.888000, 20.207300, 2.695900, 167.201996, 2.773100]
      xff_table(:, 123) = [20.180700, 3.213670, 19.113600, 0.283310, 10.905400, 20.055799, 0.773634, 51.745998, 3.029020]
      xff_table(:, 124) = [20.577999, 2.948170, 19.599001, 0.244475, 11.372700, 18.772600, 3.287190, 133.123993, 2.146780]
      xff_table(:, 125) = [20.248899, 2.920700, 19.376301, 0.250698, 11.632300, 17.821100, 0.336048, 54.945301, 2.408600]
      xff_table(:, 126) = [21.167101, 2.812190, 19.769501, 0.226836, 11.851300, 17.608299, 3.330490, 127.112999, 1.862640]
      xff_table(:, 127) = [20.803600, 2.776910, 19.559000, 0.231540, 11.936900, 16.540800, 0.612376, 43.169201, 2.090130]
      xff_table(:, 128) = [20.323500, 2.659410, 19.818600, 0.218850, 12.123300, 15.799200, 0.144583, 62.235500, 1.591800]
      xff_table(:, 129) = [22.044001, 2.773930, 19.669701, 0.222087, 12.385600, 16.766899, 2.824280, 143.643997, 2.058300]
      xff_table(:, 130) = [21.372700, 2.645200, 19.749100, 0.214299, 12.132900, 15.323000, 0.975180, 36.406502, 1.771320]
      xff_table(:, 131) = [20.941299, 2.544670, 20.053900, 0.202481, 12.466800, 14.813700, 0.296689, 45.464298, 1.242850]
      xff_table(:, 132) = [22.684500, 2.662480, 19.684700, 0.210628, 12.774000, 15.885000, 2.851370, 137.903000, 1.984860]
      xff_table(:, 133) = [21.961000, 2.527220, 19.933901, 0.199237, 12.120000, 14.178300, 1.510310, 30.871700, 1.475880]
      xff_table(:, 134) = [23.340500, 2.562700, 19.609501, 0.202088, 13.123500, 15.100900, 2.875160, 132.720993, 2.028760]
      xff_table(:, 135) = [22.552700, 2.417400, 20.110800, 0.185769, 12.067100, 13.127500, 2.074920, 27.449100, 1.194990]
      xff_table(:, 136) = [24.004200, 2.472740, 19.425800, 0.196451, 13.439600, 14.399600, 2.896040, 128.007004, 2.209630]
      xff_table(:, 137) = [23.150400, 2.316410, 20.259899, 0.174081, 11.920200, 12.157100, 2.714880, 24.824200, 0.954586]
      xff_table(:, 138) = [24.627399, 2.387900, 19.088600, 0.194200, 13.760300, 13.754600, 2.922700, 123.174004, 2.574500]
      xff_table(:, 139) = [24.006300, 2.277830, 19.950399, 0.173530, 11.803400, 11.609600, 3.872430, 26.515600, 1.363890]
      xff_table(:, 140) = [23.749701, 2.222580, 20.374500, 0.163940, 11.850900, 11.311000, 3.265030, 22.996599, 0.759344]
      xff_table(:, 141) = [25.070900, 2.253410, 19.079800, 0.181951, 13.851800, 12.933100, 3.545450, 101.398003, 2.419600]
      xff_table(:, 142) = [24.346600, 2.135530, 20.420799, 0.155525, 11.870800, 10.578200, 3.714900, 21.702900, 0.645089]
      xff_table(:, 143) = [25.897600, 2.242560, 18.218500, 0.196143, 14.316700, 12.664800, 2.953540, 115.362000, 3.583240]
      xff_table(:, 144) = [24.955900, 2.056010, 20.327101, 0.149525, 12.247100, 10.049900, 3.773000, 21.277300, 0.691967]
      xff_table(:, 145) = [26.507000, 2.180200, 17.638300, 0.202172, 14.559600, 12.189900, 2.965770, 111.874001, 4.297280]
      xff_table(:, 146) = [25.539499, 1.980400, 20.286100, 0.143384, 11.981200, 9.349720, 4.500730, 19.580999, 0.689690]
      xff_table(:, 147) = [26.904900, 2.070510, 17.294001, 0.197940, 14.558300, 11.440700, 3.638370, 92.656601, 4.567960]
      xff_table(:, 148) = [26.129601, 1.910720, 20.099400, 0.139358, 11.978800, 8.800180, 4.936760, 18.590799, 0.852795]
      xff_table(:, 149) = [27.656300, 2.073560, 16.428499, 0.223545, 14.977900, 11.360400, 2.982330, 105.703003, 5.920460]
      xff_table(:, 150) = [26.722000, 1.846590, 19.774799, 0.137290, 12.150600, 8.362250, 5.173790, 17.897400, 1.176130]
      xff_table(:, 151) = [28.181900, 2.028590, 15.885100, 0.238849, 15.154200, 10.997500, 2.987060, 102.960999, 6.756210]
      xff_table(:, 152) = [27.308300, 1.787110, 19.332001, 0.136974, 12.333900, 7.967780, 5.383480, 17.292200, 1.639290]
      xff_table(:, 153) = [28.664101, 1.988900, 15.434500, 0.257119, 15.308700, 10.664700, 2.989630, 100.417000, 7.566720]
      xff_table(:, 154) = [28.120899, 1.785030, 17.681700, 0.159970, 13.333500, 8.183040, 5.146570, 20.389999, 3.709830]
      xff_table(:, 155) = [27.891701, 1.732720, 18.761400, 0.138790, 12.607200, 7.644120, 5.476470, 16.815300, 2.260010]
      xff_table(:, 156) = [28.947599, 1.901820, 15.220800, 9.985190, 15.100000, 0.261033, 3.716010, 84.329803, 7.976280]
      xff_table(:, 157) = [28.462799, 1.682160, 18.121000, 0.142292, 12.842900, 7.337270, 5.594150, 16.353500, 2.975730]
      xff_table(:, 158) = [29.143999, 1.832620, 15.172600, 9.599900, 14.758600, 0.275116, 4.300130, 72.028999, 8.581540]
      xff_table(:, 159) = [28.813101, 1.591360, 18.460100, 0.128903, 12.728500, 6.762320, 5.599270, 14.036600, 2.396990]
      xff_table(:, 160) = [29.202400, 1.773330, 15.229300, 9.370460, 14.513500, 0.295977, 4.764920, 63.364399, 9.243540]
      xff_table(:, 161) = [29.158701, 1.507110, 18.840700, 0.116741, 12.826800, 6.315240, 5.386950, 12.424400, 1.785550]
      xff_table(:, 162) = [29.081800, 1.720290, 15.430000, 9.225900, 14.432700, 0.321703, 5.119820, 57.056000, 9.887500]
      xff_table(:, 163) = [29.493601, 1.427550, 19.376301, 0.104621, 13.054400, 5.936670, 5.064120, 11.197200, 1.010740]
      xff_table(:, 164) = [28.762100, 1.671910, 15.718900, 9.092270, 14.556400, 0.350500, 5.441740, 52.086102, 10.472000]
      xff_table(:, 165) = [28.189400, 1.629030, 16.155001, 8.979480, 14.930500, 0.382661, 5.675890, 48.164700, 11.000500]
      xff_table(:, 166) = [30.419001, 1.371130, 15.263700, 6.847060, 14.745800, 0.165191, 5.067950, 18.003000, 6.498040]
      xff_table(:, 167) = [27.304899, 1.592790, 16.729601, 8.865530, 15.611500, 0.417916, 5.833770, 45.001099, 11.472200]
      xff_table(:, 168) = [30.415600, 1.343230, 15.862000, 7.109090, 13.614500, 0.204633, 5.820080, 20.325399, 8.279030]
      xff_table(:, 169) = [30.705799, 1.309230, 15.551200, 6.719830, 14.232600, 0.167252, 5.536720, 17.491100, 6.968240]
      xff_table(:, 170) = [27.005899, 1.512930, 17.763901, 8.811740, 15.713100, 0.424593, 5.783700, 38.610298, 11.688300]
      xff_table(:, 171) = [29.842899, 1.329270, 16.722401, 7.389790, 13.215300, 0.263297, 6.352340, 22.942600, 9.853290]
      xff_table(:, 172) = [30.961201, 1.248130, 15.982900, 6.608340, 13.734800, 0.168640, 5.920340, 16.939199, 7.395340]
      xff_table(:, 173) = [16.881901, 0.461100, 18.591299, 8.621600, 25.558201, 1.482600, 5.860000, 36.395599, 12.065800]
      xff_table(:, 174) = [28.010900, 1.353210, 17.820400, 7.739500, 14.335900, 0.356752, 6.580770, 26.404301, 11.229900]
      xff_table(:, 175) = [30.688601, 1.219900, 16.902901, 6.828720, 12.780100, 0.212867, 6.523540, 18.659000, 9.096800]
      xff_table(:, 176) = [20.680901, 0.545000, 19.041700, 8.448400, 21.657499, 1.572900, 5.967600, 38.324600, 12.608900]
      xff_table(:, 177) = [25.085300, 1.395070, 18.497299, 7.651050, 16.888300, 0.443378, 6.482160, 28.226200, 12.020500]
      xff_table(:, 178) = [29.564100, 1.211520, 18.059999, 7.056390, 12.837400, 0.284738, 6.899120, 20.748199, 10.626800]
      xff_table(:, 179) = [27.544600, 0.655150, 19.158400, 8.707510, 15.538000, 1.963470, 5.525930, 45.814899, 13.174600]
      xff_table(:, 180) = [21.398500, 1.471100, 20.472300, 0.517394, 18.747801, 7.434630, 6.828470, 28.848200, 12.525800]
      xff_table(:, 181) = [30.869499, 1.100800, 18.384100, 6.538520, 11.932800, 0.219074, 7.005740, 17.211399, 9.802700]
      xff_table(:, 182) = [31.061701, 0.690200, 13.063700, 2.357600, 18.441999, 8.618000, 5.969600, 47.257900, 13.411800]
      xff_table(:, 183) = [21.788601, 1.336600, 19.568199, 0.488383, 19.140600, 6.772700, 7.011070, 23.813200, 12.473400]
      xff_table(:, 184) = [32.124401, 1.005660, 18.800301, 6.109260, 12.017500, 0.147041, 6.968860, 14.714000, 8.084280]
      xff_table(:, 185) = [33.368900, 0.704000, 12.951000, 2.923800, 16.587700, 8.793700, 6.469200, 48.009300, 13.578200]
      xff_table(:, 186) = [21.805300, 1.235600, 19.502600, 6.241490, 19.105301, 0.469999, 7.102950, 20.318501, 12.471100]
      xff_table(:, 187) = [33.536400, 0.916540, 25.094601, 0.039042, 19.249701, 5.714140, 6.915550, 12.828500, -6.799400]
      xff_table(:, 188) = [34.672600, 0.700999, 15.473300, 3.550780, 13.113800, 9.556420, 7.025880, 47.004501, 13.677000]
      xff_table(:, 189) = [35.316299, 0.685870, 19.021099, 3.974580, 9.498870, 11.382400, 7.425180, 45.471500, 13.710800]
      xff_table(:, 190) = [35.563099, 0.663100, 21.281601, 4.069100, 8.003700, 14.042200, 7.443300, 44.247299, 13.690500]
      xff_table(:, 191) = [35.929901, 0.646453, 23.054701, 4.176190, 12.143900, 23.105200, 2.112530, 150.645004, 13.724700]
      xff_table(:, 192) = [35.763000, 0.616341, 22.906401, 3.871350, 12.473900, 19.988701, 3.210970, 142.324997, 13.621100]
      xff_table(:, 193) = [35.215000, 0.604909, 21.670000, 3.576700, 7.913420, 12.601000, 7.650780, 29.843599, 13.543100]
      xff_table(:, 194) = [35.659698, 0.589092, 23.103201, 3.651550, 12.597700, 18.599001, 4.086550, 117.019997, 13.526600]
      xff_table(:, 195) = [35.173599, 0.579689, 22.111200, 3.414370, 8.192160, 12.918700, 7.055450, 25.944300, 13.463700]
      xff_table(:, 196) = [35.564499, 0.563359, 23.421900, 3.462040, 12.747300, 17.830900, 4.807030, 99.172203, 13.431400]
      xff_table(:, 197) = [35.100700, 0.555054, 22.441799, 3.244980, 9.785540, 13.466100, 5.294440, 23.953300, 13.376000]
      xff_table(:, 198) = [35.884701, 0.547751, 23.294800, 3.415190, 14.189100, 16.923500, 4.172870, 105.250999, 13.428700]
      xff_table(:, 199) = [36.022800, 0.529300, 23.412800, 3.325300, 14.949100, 16.092699, 4.188000, 100.612999, 13.396600]
      xff_table(:, 200) = [35.574699, 0.520480, 22.525900, 3.122930, 12.216500, 12.714800, 5.370730, 26.339399, 13.309200]
      xff_table(:, 201) = [35.371498, 0.516598, 22.532600, 3.050530, 12.029100, 12.572300, 4.798400, 23.458200, 13.267100]
      xff_table(:, 202) = [34.850899, 0.507079, 22.758400, 2.890300, 14.009900, 13.176700, 1.214570, 25.201700, 13.166500]
      xff_table(:, 203) = [36.187401, 0.511929, 23.596399, 3.253960, 15.640200, 15.362200, 4.185500, 97.490799, 13.357300]
      xff_table(:, 204) = [35.707401, 0.502322, 22.613001, 3.038070, 12.989800, 12.144900, 5.432270, 25.492800, 13.254400]
      xff_table(:, 205) = [35.510300, 0.498626, 22.578699, 2.966270, 12.776600, 11.948400, 4.921590, 22.750200, 13.211600]
      xff_table(:, 206) = [35.013599, 0.489810, 22.728600, 2.810990, 14.388400, 12.330000, 1.756690, 22.658100, 13.113000]
      xff_table(:, 207) = [36.525398, 0.499384, 23.808300, 3.263710, 16.770700, 14.945500, 3.479470, 105.980003, 13.381200]
      xff_table(:, 208) = [35.840000, 0.484938, 22.716900, 2.961180, 13.580700, 11.533100, 5.660160, 24.399200, 13.199100]
      xff_table(:, 209) = [35.649300, 0.481422, 22.646000, 2.890200, 13.359500, 11.316000, 5.188310, 21.830099, 13.155500]
      xff_table(:, 210) = [35.173599, 0.473204, 22.718100, 2.738480, 14.763500, 11.553000, 2.286780, 20.930300, 13.058200]
      xff_table(:, 211) = [36.670601, 0.483629, 24.099199, 3.206470, 17.341499, 14.313600, 3.493310, 102.273003, 13.359200]
      xff_table(:, 212) = [36.648800, 0.465154, 24.409599, 3.089970, 17.399000, 13.434600, 4.216650, 88.483398, 13.288700]

      eff_table = 0D0
      eff_table(:, 1) = [0.0088, 0.0449, 0.1481, 0.2356, 0.0914, 0.1152, 1.0867, 4.9755, 16.5591, 43.2743]
      eff_table(:, 5) = [0.0084, 0.0443, 0.1314, 0.1671, 0.0666, 0.0596, 0.5360, 2.4274, 7.7852, 0.3126]
      eff_table(:, 6) = [0.0478, 0.2048, 0.5253, 1.5225, 0.9853, 0.2258, 2.1032, 12.9349, 50.7501, 136.6280]
      eff_table(:, 8) = [0.0423, 0.1874, 0.6019, 1.4311, 0.7891, 0.1445, 1.4180, 8.1165, 27.9705, 74.8684]
      eff_table(:, 10) = [0.0436, 0.1898, 0.6788, 1.3273, 0.5544, 0.1207, 1.1595, 6.2474, 21.0460, 59.3619]
      eff_table(:, 11) = [0.0489, 0.2091, 0.7537, 1.1420, 0.3555, 0.1140, 1.0825, 5.4281, 17.8811, 51.1341]
      eff_table(:, 13) = [0.0267, 0.1328, 0.5301, 1.1020, 0.4215, 0.0541, 0.5165, 2.8207, 10.6297, 34.3764]
      eff_table(:, 14) = [0.0365, 0.1729, 0.5805, 0.8814, 0.3121, 0.0652, 0.6184, 2.9449, 9.6298, 8.2194]
      eff_table(:, 17) = [0.0382, 0.1822, 0.5972, 0.7707, 0.2130, 0.0613, 0.5753, 2.6858, 8.8214, 5.6668]
      eff_table(:, 19) = [0.0380, 0.1785, 0.5494, 0.6942, 0.1918, 0.0554, 0.5087, 2.2639, 7.3316, 1.6912]
      eff_table(:, 20) = [0.1260, 0.6442, 0.8893, 1.8197, 1.2988, 0.1684, 1.7150, 8.8386, 50.8265, 147.2073]
      eff_table(:, 22) = [0.1130, 0.5575, 0.9046, 2.1580, 1.4735, 0.1356, 1.3579, 6.9255, 32.3165, 92.1138]
      eff_table(:, 24) = [0.1165, 0.5504, 1.0179, 2.6295, 1.5711, 0.1295, 1.2619, 6.8242, 28.4577, 88.4750]
      eff_table(:, 26) = [0.0567, 0.3365, 0.8104, 2.4960, 2.1186, 0.0582, 0.6155, 3.2522, 16.7929, 57.6767]
      eff_table(:, 29) = [0.1005, 0.4615, 1.0663, 2.5854, 1.2725, 0.0977, 0.9084, 4.9654, 18.5471, 54.3648]
      eff_table(:, 30) = [0.0915, 0.4312, 1.0847, 2.4671, 1.0852, 0.0838, 0.7788, 4.3462, 15.5846, 44.6365]
      eff_table(:, 31) = [0.0799, 0.3891, 1.0037, 2.3332, 1.0507, 0.0694, 0.6443, 3.5351, 12.5058, 35.8633]
      eff_table(:, 33) = [0.1044, 0.4551, 1.4232, 2.1533, 0.4459, 0.0853, 0.7701, 4.4684, 14.5864, 41.2474]
      eff_table(:, 34) = [0.2149, 0.8703, 2.4999, 2.3591, 3.0318, 0.1660, 1.6906, 8.7447, 46.7825, 165.6923]
      eff_table(:, 36) = [0.2355, 0.9916, 2.3959, 3.7252, 2.5647, 0.1742, 1.8329, 8.8407, 47.4583, 134.9613]
      eff_table(:, 38) = [0.4636, 2.0802, 2.9003, 1.4193, 2.4323, 0.3682, 4.0312, 22.6493, 71.8200, 103.3691]
      eff_table(:, 40) = [0.2123, 0.8960, 2.1765, 3.0436, 2.4439, 0.1399, 1.4568, 6.7534, 33.1168, 101.8238]
      eff_table(:, 44) = [0.2369, 1.0774, 2.1894, 3.0825, 1.7190, 0.1505, 1.6392, 7.5691, 36.8741, 107.8517]
      eff_table(:, 48) = [0.1970, 0.8228, 2.0200, 2.1717, 1.7516, 0.1197, 1.1985, 5.4097, 25.2361, 94.4290]
      eff_table(:, 51) = [0.1943, 0.8190, 1.9296, 2.4968, 2.0625, 0.1135, 1.1313, 5.0341, 24.1798, 80.5598]
      eff_table(:, 55) = [0.1929, 0.8239, 1.8689, 2.3694, 1.9060, 0.1087, 1.0806, 4.7637, 22.8500, 76.7309]
      eff_table(:, 58) = [0.2186, 0.9861, 1.8540, 2.3258, 1.4685, 0.1182, 1.2300, 5.4177, 25.7602, 80.8542]
      eff_table(:, 61) = [0.2313, 1.0657, 1.8229, 2.2609, 1.1883, 0.1210, 1.2691, 5.6870, 27.0917, 83.0285]
      eff_table(:, 64) = [0.3501, 1.6558, 1.9582, 0.2134, 1.4109, 0.1867, 1.9917, 11.3396, 53.2619, 63.2520]
      eff_table(:, 67) = [0.1780, 0.8096, 1.6744, 1.9499, 1.4495, 0.0876, 0.8650, 3.8612, 18.8726, 64.7016]
      eff_table(:, 69) = [0.2135, 0.9768, 1.6669, 2.5662, 1.6790, 0.1020, 1.0219, 4.6275, 22.8742, 80.1535]
      eff_table(:, 71) = [0.2135, 0.9761, 1.6555, 2.8938, 1.6356, 0.0989, 0.9845, 4.5527, 21.5563, 70.3903]
      eff_table(:, 73) = [0.2059, 0.9518, 1.6372, 3.0490, 1.4756, 0.0926, 0.9182, 4.3291, 19.2996, 58.9329]
      eff_table(:, 74) = [0.1574, 0.7614, 1.4834, 3.0016, 1.7978, 0.0686, 0.6808, 3.1163, 14.3458, 44.0455]
      eff_table(:, 75) = [0.1899, 0.8983, 1.6358, 3.1845, 1.1518, 0.0810, 0.7957, 3.9054, 15.7701, 45.6124]
      eff_table(:, 77) = [0.1742, 0.8447, 1.5944, 3.1507, 1.1338, 0.0723, 0.7123, 3.5192, 13.7724, 39.1148]
      eff_table(:, 78) = [0.3781, 1.4904, 3.5753, 3.0031, 3.3272, 0.1557, 1.5347, 9.9947, 51.4251, 185.9828]
      eff_table(:, 80) = [0.3723, 1.4598, 3.5124, 4.4612, 3.3031, 0.1480, 1.4643, 9.2320, 49.8807, 148.0937]
      eff_table(:, 82) = [0.3234, 1.2737, 3.2115, 4.0563, 3.7962, 0.1244, 1.1948, 7.2756, 34.1430, 111.2079]
      eff_table(:, 84) = [0.2997, 1.1879, 3.1075, 3.9740, 3.5769, 0.1121, 1.0638, 6.3891, 28.7081, 97.4289]
      eff_table(:, 86) = [0.1680, 0.9370, 2.7300, 3.8150, 3.0053, 0.0597, 0.6524, 4.4317, 19.5540, 85.5011]
      eff_table(:, 89) = [0.3069, 1.1714, 3.2293, 3.4254, 2.1224, 0.1101, 1.0222, 5.9613, 25.1965, 93.5831]
      eff_table(:, 93) = [0.2928, 1.1267, 3.1675, 3.6619, 2.5942, 0.1020, 0.9481, 5.4713, 23.8153, 82.8991]
      eff_table(:, 94) = [0.2604, 1.0442, 3.0761, 3.2175, 1.9448, 0.0887, 0.8240, 4.8278, 19.8977, 80.4566]
      eff_table(:, 97) = [0.2713, 1.0556, 3.1416, 3.0451, 1.7179, 0.0907, 0.8324, 4.7702, 19.7862, 80.2540]
      eff_table(:, 100) = [0.2003, 0.8779, 2.6135, 2.8594, 1.0258, 0.0659, 0.6111, 3.5563, 12.7638, 44.4283]
      eff_table(:, 103) = [0.2739, 1.0503, 3.1564, 2.7543, 1.4328, 0.0881, 0.8028, 4.4451, 18.7011, 79.2633]
      eff_table(:, 106) = [0.3072, 1.1303, 3.2046, 2.9329, 1.6560, 0.0966, 0.8856, 4.6273, 20.6789, 73.4723]
      eff_table(:, 108) = [0.3564, 1.3011, 3.2424, 3.4839, 2.0459, 0.1091, 1.0452, 5.0900, 24.6578, 88.0513]
      eff_table(:, 110) = [0.2966, 1.1157, 3.0973, 3.8156, 2.5281, 0.0896, 0.8268, 4.2242, 20.6900, 71.3399]
      eff_table(:, 113) = [0.2725, 1.0651, 2.9940, 4.0697, 2.5682, 0.0809, 0.7488, 3.8710, 18.8800, 60.6499]
      eff_table(:, 116) = [0.2422, 0.9692, 2.8114, 4.1509, 2.8161, 0.0708, 0.6472, 3.3609, 16.0752, 50.1724]
      eff_table(:, 117) = [0.2617, 1.0325, 2.8097, 4.4809, 2.3190, 0.0749, 0.6914, 3.4634, 16.3603, 48.2522]
      eff_table(:, 119) = [0.2334, 0.9496, 2.6381, 4.4680, 2.5020, 0.0655, 0.6050, 3.0389, 14.0809, 41.0005]
      eff_table(:, 120) = [0.5713, 2.4866, 4.9795, 4.0198, 4.4403, 0.1626, 1.8213, 11.1049, 49.0568, 202.9987]
      eff_table(:, 122) = [0.5229, 2.2874, 4.7243, 5.0807, 5.6389, 0.1434, 1.6019, 9.4511, 42.7685, 148.4969]
      eff_table(:, 124) = [0.5461, 2.3856, 5.0653, 5.7601, 4.0463, 0.1479, 1.6552, 10.0059, 47.3245, 145.8464]
      eff_table(:, 126) = [0.2227, 1.0760, 2.9482, 5.8496, 7.1834, 0.0571, 0.5946, 3.2022, 16.4253, 95.7030]
      eff_table(:, 129) = [0.5237, 2.2913, 4.6161, 4.7233, 4.8173, 0.1360, 1.5068, 8.8213, 41.9536, 141.2424]
      eff_table(:, 132) = [0.5368, 2.3301, 4.6058, 4.6621, 4.4622, 0.1378, 1.5140, 8.8719, 43.5967, 141.8065]
      eff_table(:, 134) = [0.5232, 2.2627, 4.4552, 4.4787, 4.5073, 0.1317, 1.4336, 8.3087, 40.6010, 135.9196]
      eff_table(:, 136) = [0.5162, 2.2302, 4.3449, 4.3598, 4.4292, 0.1279, 1.3811, 7.9629, 39.1213, 132.7846]
      eff_table(:, 138) = [0.5272, 2.2844, 4.3361, 4.3178, 4.0908, 0.1285, 1.3943, 8.1081, 40.9631, 134.1233]
      eff_table(:, 141) = [0.9664, 3.4052, 5.0803, 1.4991, 4.2528, 0.2641, 2.6586, 16.2213, 80.2060, 92.5359]
      eff_table(:, 143) = [0.5110, 2.1570, 4.0308, 3.9936, 4.2466, 0.1210, 1.2704, 7.1368, 35.0354, 123.5062]
      eff_table(:, 145) = [0.4974, 2.1097, 3.8906, 3.8100, 4.3084, 0.1157, 1.2108, 6.7377, 32.4150, 116.9225]
      eff_table(:, 147) = [0.4679, 1.9693, 3.7191, 3.9632, 4.2432, 0.1069, 1.0994, 5.9769, 27.1491, 96.3119]
      eff_table(:, 149) = [0.5034, 2.1088, 3.8232, 3.7299, 3.8963, 0.1141, 1.1769, 6.6087, 33.4332, 116.4913]
      eff_table(:, 151) = [0.4839, 2.0262, 3.6851, 3.5874, 4.0037, 0.1081, 1.1012, 6.1114, 30.3728, 110.5988]
      eff_table(:, 153) = [0.5221, 2.1695, 3.7567, 3.6685, 3.4274, 0.1148, 1.1860, 6.7520, 35.6807, 118.0692]
      eff_table(:, 156) = [0.4680, 1.9466, 3.5428, 3.8490, 3.6594, 0.1015, 1.0195, 5.6058, 27.4899, 95.2846]
      eff_table(:, 158) = [0.4048, 1.7370, 3.3399, 3.9448, 3.7293, 0.0868, 0.8585, 4.6378, 21.6900, 80.2408]
      eff_table(:, 160) = [0.3835, 1.6747, 3.2986, 4.0462, 3.4303, 0.0810, 0.8020, 4.3545, 19.9644, 73.6337]
      eff_table(:, 162) = [0.3661, 1.6191, 3.2455, 4.0856, 3.2064, 0.0761, 0.7543, 4.0952, 18.2886, 68.0967]
      eff_table(:, 164) = [0.3933, 1.6973, 3.4202, 4.1274, 2.6158, 0.0806, 0.7972, 4.4237, 19.5692, 68.7477]
      eff_table(:, 165) = [0.3854, 1.6555, 3.4129, 4.1111, 2.4106, 0.0787, 0.7638, 4.2441, 18.3700, 65.1071]
      eff_table(:, 167) = [0.3510, 1.5620, 3.2946, 4.0615, 2.4382, 0.0706, 0.6904, 3.8266, 16.0812, 58.7638]
      eff_table(:, 170) = [0.3083, 1.4158, 2.9662, 3.9349, 2.1709, 0.0609, 0.5993, 3.1921, 12.5285, 49.7675]
      eff_table(:, 173) = [0.3055, 1.3945, 2.9617, 3.8990, 2.0026, 0.0596, 0.5827, 3.1035, 11.9693, 47.9106]
      eff_table(:, 176) = [0.3593, 1.5736, 3.5237, 3.8109, 1.6953, 0.0694, 0.6758, 3.8457, 15.6203, 56.6614]
      eff_table(:, 179) = [0.3511, 1.5489, 3.5676, 4.0900, 2.5251, 0.0672, 0.6522, 3.7420, 15.9791, 65.1354]
      eff_table(:, 182) = [0.3540, 1.5453, 3.5975, 4.3152, 2.7743, 0.0668, 0.6465, 3.6968, 16.2056, 61.4909]
      eff_table(:, 185) = [0.3530, 1.5258, 3.5815, 4.5532, 3.0714, 0.0661, 0.6324, 3.5906, 15.9962, 57.5760]
      eff_table(:, 188) = [0.3673, 1.5772, 3.7079, 4.8582, 2.8440, 0.0678, 0.6527, 3.7396, 17.0668, 55.9789]
      eff_table(:, 189) = [0.3547, 1.5206, 3.5621, 5.0184, 3.0075, 0.0649, 0.6188, 3.4696, 15.6090, 49.4818]
      eff_table(:, 190) = [0.4586, 1.7781, 3.9877, 5.7273, 1.5460, 0.0831, 0.7840, 4.3599, 20.0128, 62.1535]
      eff_table(:, 191) = [0.8282, 2.9941, 5.6597, 4.9292, 4.2889, 0.1515, 1.6163, 9.7752, 42.8480, 190.7366]
      eff_table(:, 192) = [1.4129, 4.4269, 7.0460, 1.0573, 8.6430, 0.2921, 3.1381, 19.6767, 102.0436, 113.9798]
      eff_table(:, 194) = [0.7169, 2.5710, 5.1791, 6.3484, 5.6474, 0.1263, 1.2900, 7.3686, 32.4490, 118.0558]
      eff_table(:, 196) = [0.6958, 2.4936, 5.1269, 6.6988, 5.0799, 0.1211, 1.2247, 6.9398, 30.0991, 105.1960]
      eff_table(:, 198) = [1.2502, 4.2284, 7.0489, 1.1390, 5.8222, 0.2415, 2.6442, 16.3313, 73.5757, 91.9401]
      eff_table(:, 199) = [0.6410, 2.2643, 4.8713, 5.9287, 5.3935, 0.1097, 1.0644, 5.7907, 25.0261, 101.3899]
      eff_table(:, 203) = [0.6938, 2.4652, 5.1227, 5.5965, 4.8543, 0.1171, 1.1757, 6.4053, 27.5217, 103.0482]
      eff_table(:, 207) = [0.6902, 2.4509, 5.1284, 5.0339, 4.8575, 0.1153, 1.1545, 6.2291, 27.0741, 111.3150]
      eff_table(:, 211) = [0.7577, 2.7264, 5.4184, 4.8198, 4.1013, 0.1257, 1.3044, 7.1035, 32.4649, 118.8647]
      eff_table(:, 212) = [0.7567, 2.7565, 5.4364, 5.1918, 3.5643, 0.1239, 1.2979, 7.0798, 32.7871, 110.1512]

      found_element = .FALSE.
      DO i = 1, 212
         IF (TRIM(elem_str) == element_table(i)) THEN
            xff(:) = xff_table(:, i)
            nsl = nsl_table(i)
            eff(:) = eff_table(:, i)
            found_element = .TRUE.
            !write(*,*) 'found sls for element',trim(element_str),i
         END IF
      END DO

      IF (.NOT. found_element) THEN
         WRITE (*, *) 'Supercell error: Scattering lengths are not known for element symbol ', TRIM(elem_str)
         STOP
      END IF

   END SUBROUTINE sl_from_element

END PROGRAM scatty

SUBROUTINE write_fit(num_q, num_axis, qvec_from_iq_c, q_c2o, hkl2xyz, valid, origin, axis, &
& intensity, title, str, colourmap, min, max, output_vtk, output_ascii, output_ppm, suppress_errors)

   IMPLICIT NONE

   INTEGER(8), INTENT(in) :: num_q
   INTEGER, INTENT(in) :: num_axis(3), colourmap
   LOGICAL, INTENT(in) :: valid(num_q), output_vtk, output_ascii, output_ppm
   REAL(8), INTENT(in) :: qvec_from_iq_c(3, num_q), origin(3), axis(3, 3), q_c2o(3, 3), hkl2xyz(3, 3)
   REAL(8), INTENT(in) :: intensity(num_q)
   LOGICAL, INTENT(in) :: suppress_errors
   REAL, INTENT(inout):: max, min
   CHARACTER(128), INTENT(in) :: title, str
   INTEGER :: grid_count(num_axis(1), num_axis(2), num_axis(3)), qx, qy, qz, dim, dim_axis(3), i, j, width, height
   INTEGER(8) :: iq
   REAL(8) :: i_q_grid(num_axis(1), num_axis(2), num_axis(3)), q(3), axis_length, dq_axis(3), axis_norm(3, 3), origin_orth(3)
   REAL(8) :: axis_orth(3, 3), q_orth(3, 3), dp(3), q_mod_tmp
   REAL(8) :: q_mod(num_axis(1), num_axis(2), num_axis(3))
   REAL :: colourbar(20, 256)
   REAL, ALLOCATABLE :: intensity_img(:, :)
   REAL(8), PARAMETER ::  eps = 10D0*EPSILON(0.0)
   CHARACTER(128) :: colourbar_str
   LOGICAL :: use_hkl2xyz

   dq_axis = 0D0
   axis_norm = 0D0
   colourbar_str = 'colourbar'

   use_hkl2xyz = .FALSE.
   DO i = 1, 3
      DO j = 1, 3
         IF (ABS(hkl2xyz(i, j)) > eps) use_hkl2xyz = .TRUE.
      END DO
   END DO

   IF (use_hkl2xyz) THEN
      WRITE (*, *) 'NOTE: Using user-defined (h,k,l) orthogonalization matrix:'
      WRITE (*, *) 'HKL_TO_X =', hkl2xyz(1, :)
      WRITE (*, *) 'HKL_TO_Y =', hkl2xyz(2, :)
      WRITE (*, *) 'HKL_TO_Z =', hkl2xyz(3, :)
      q_orth = hkl2xyz
      origin_orth = origin(:)
   ELSE
      q_orth = q_c2o
      origin_orth = MATMUL(q_orth, origin(:))
   END IF

   DO i = 1, 3
      IF (use_hkl2xyz) THEN
         axis_orth(i, :) = axis(i, :)
      ELSE
         axis_orth(i, :) = MATMUL(q_orth, axis(i, :))
      END IF
      axis_length = DOT_PRODUCT(axis_orth(i, :), axis_orth(i, :))
      IF (ABS(axis_length) > eps) THEN
         axis_length = SQRT(axis_length)
         dq_axis(i) = axis_length/(1D0*(num_axis(i) - 1))
         axis_norm(i, :) = axis_orth(i, :)/(1D0*(num_axis(i) - 1))
      ELSE
         IF (num_axis(i) .NE. 1) STOP 'WRITE_FIT error: Invalid number of data points along some dimension(s)'
      END IF
      DO j = 1, i - 1
         IF (ABS(DOT_PRODUCT(axis_orth(i, :), axis_orth(j, :))) > eps) THEN
            WRITE (*, *) 'WARNING: Reciprocal lattice vectors do not seem to be orthogonal; dot product = '&
            &, DOT_PRODUCT(axis_orth(i, :), axis_orth(j, :))
            WRITE (*, *) ''
            STOP
         END IF
      END DO
   END DO

   IF (ABS(DOT_PRODUCT(dq_axis, dq_axis)) < eps) RETURN

   i_q_grid = 0D0
   grid_count = 0
   dim = 3
   dim_axis = 1

   DO i = 1, 3
      IF (num_axis(i) == 1) THEN
         dim = dim - 1
         dim_axis(i) = 0
      END IF
      dp(i) = DOT_PRODUCT(axis_norm(i, :), axis_norm(i, :))
   END DO

! rebin data onto a grid for visualisation
   DO iq = 1, num_q
      IF (.NOT. valid(iq)) CYCLE
      q = MATMUL(q_orth, qvec_from_iq_c(:, iq))
      q_mod_tmp = SQRT(DOT_PRODUCT(q, q))
      q = q - origin_orth(:)
      IF (ABS(dp(1)) > eps) THEN
         qx = NINT(DOT_PRODUCT(q(:), axis_norm(1, :))/dp(1)) + 1
      ELSE
         qx = 1
      END IF
      IF ((qx < 1) .OR. (qx > num_axis(1))) THEN
         IF (.NOT. suppress_errors) WRITE (*, *) 'WRITE_FIT warning: X out of bounds', qx, num_axis(1)
         CYCLE
      END IF
      IF (ABS(dp(2)) > eps) THEN
         qy = NINT(DOT_PRODUCT(q(:), axis_norm(2, :))/dp(2)) + 1
      ELSE
         qy = 1
      END IF
      IF ((qy < 1) .OR. (qy > num_axis(2))) THEN
         IF (.NOT. suppress_errors) WRITE (*, *) 'WRITE_FIT warning: Y out of bounds', qy, num_axis(2)
         CYCLE
      END IF
      IF (ABS(dp(3)) > eps) THEN
         qz = NINT(DOT_PRODUCT(q(:), axis_norm(3, :))/dp(3)) + 1
      ELSE
         qz = 1
      END IF
      IF ((qz < 1) .OR. (qz > num_axis(3))) THEN
         IF (.NOT. suppress_errors) WRITE (*, *) 'WRITE_FIT warning: Z out of bounds', qz, num_axis(3)
         CYCLE
      END IF
      q_mod(qx, qy, qz) = q_mod_tmp
      i_q_grid(qx, qy, qz) = intensity(iq) + i_q_grid(qx, qy, qz)
      grid_count(qx, qy, qz) = grid_count(qx, qy, qz) + 1
   END DO

   IF (output_vtk) THEN
      OPEN (11, file=TRIM(ADJUSTL(str))//'.vtk', status='replace')
      WRITE (11, '(a)') '# vtk DataFile Version 2.0'
      WRITE (11, '(a)') 'TITLE diffuse scattering'
      WRITE (11, '(a)') 'ASCII'
      WRITE (11, '(a)') 'DATASET STRUCTURED_POINTS'
      WRITE (11, '(a,3i6)') ADJUSTL('DIMENSIONS'), num_axis(1), num_axis(2), num_axis(3)
      WRITE (11, '(a,3f12.6)') 'ORIGIN', origin_orth(:)
      WRITE (11, '(a,3f12.6)') 'SPACING', dq_axis(1), dq_axis(2), dq_axis(3)
      WRITE (11, '(a,i12)') 'POINT_DATA', PRODUCT(num_axis)
      WRITE (11, '(a)') 'SCALARS '//TRIM(title)//' float'
      WRITE (11, '(a)') 'LOOKUP_TABLE default'
      DO qz = 1, num_axis(3)
         DO qy = 1, num_axis(2)
            DO qx = 1, num_axis(1)
               IF (grid_count(qx, qy, qz) > 0) i_q_grid(qx, qy, qz) = i_q_grid(qx, qy, qz)/grid_count(qx, qy, qz)
               WRITE (11, '(f20.12)') i_q_grid(qx, qy, qz)
               IF (grid_count(qx, qy, qz) > 1) THEN
                  IF (.NOT. suppress_errors) THEN
                     WRITE (*, *) 'WRITE_FIT warning: grid count > 1', qx, qy, qz, grid_count(qx, qy, qz)
                  END IF
               END IF
            END DO
         END DO
      END DO
      CLOSE (11)
   END IF

   IF (dim == 1 .AND. output_ascii) THEN
      OPEN (11, file=TRIM(ADJUSTL(str))//'.txt', status='replace')
      IF (dim_axis(1) == 1) THEN
      DO qx = 1, num_axis(1)
         WRITE (11, '(f20.12)') i_q_grid(qx, 1, 1)
      END DO
      ELSE IF (dim_axis(2) == 1) THEN
      DO qy = 1, num_axis(2)
         WRITE (11, '(f20.12)') i_q_grid(1, qy, 1)
      END DO
      ELSE IF (dim_axis(2) == 1) THEN
      DO qz = 1, num_axis(3)
         WRITE (11, '(f20.12)') i_q_grid(1, 1, qz)
      END DO
      END IF
      CLOSE (11)
   END IF

   IF (dim == 2 .AND. output_ascii) THEN
      OPEN (11, file=TRIM(ADJUSTL(str))//'.txt', status='replace')
      !open(12,file=trim(adjustl(str))//'qmod.txt',status='replace')
      IF ((dim_axis(1) == 1) .AND. (dim_axis(2) == 1)) THEN
         width = num_axis(1)
         height = num_axis(2)
         ALLOCATE (intensity_img(width, height))
         intensity_img = i_q_grid(:, :, 1)
         DO qy = 1, num_axis(2)
            WRITE (11, *) (i_q_grid(qx, qy, 1), qx=1, num_axis(1))
            ! write(12,*) (q_mod(qx,qy,1),qx=1,num_axis(1))
         END DO
      ELSE IF ((dim_axis(1) == 1) .AND. (dim_axis(3) == 1)) THEN
         width = num_axis(1)
         height = num_axis(3)
         ALLOCATE (intensity_img(width, height))
         intensity_img = i_q_grid(:, 1, :)
         DO qz = 1, num_axis(3)
            WRITE (11, *) (i_q_grid(qx, 1, qz), qx=1, num_axis(1))
            ! write(12,*)  (q_mod(qx,1,qz),qx=1,num_axis(1))
         END DO
      ELSE IF ((dim_axis(2) == 1) .AND. (dim_axis(3) == 1)) THEN
         width = num_axis(2)
         height = num_axis(3)
         ALLOCATE (intensity_img(width, height))
         intensity_img = i_q_grid(1, :, :)
         DO qz = 1, num_axis(3)
            WRITE (11, *) (i_q_grid(1, qy, qz), qy=1, num_axis(2))
            ! WRITE(12,*) (q_mod(1,qy,qz),qy=1,num_axis(2))
         END DO
      END IF
      CLOSE (11)
   END IF

   IF (dim == 2 .AND. output_ppm) THEN
      CALL write_ppm(intensity_img, width, height, colourmap, str, max, min)
   END IF

END SUBROUTINE write_fit

SUBROUTINE calc_chi_sq(diffract_calc, diffract_expt, errors_sq, valid, pw_nq, q,&
& refine_scale, refine_flat_bgr, refine_linear_bgr, chi_sq, scale, flat_bgr, linear_bgr)

! Calculates chi squared and scale and background factors from calc. and expt. NS intensities

   INTEGER(8), INTENT(in) :: pw_nq
   LOGICAL, INTENT(in) :: refine_scale, refine_flat_bgr, refine_linear_bgr
   REAL(8), INTENT(in) :: diffract_calc(pw_nq), diffract_expt(pw_nq), errors_sq(pw_nq), q(pw_nq)
   LOGICAL, INTENT(in) :: valid(pw_nq) ! If there are no data for a given Q-point, valid=false
   REAL(8), INTENT(inout) :: scale, flat_bgr, linear_bgr
   REAL(8), INTENT(out) :: chi_sq
   REAL(8) :: sum_ss, sum_fs, sum_qs, sum_s, sum_e, sum_f, sum_qq, sum_qf, sum_q, num, denom
   REAL(8), PARAMETER ::  eps = 10D0*EPSILON(0.0)

   IF (ALL(.NOT. valid) .EQV. .TRUE.) THEN
      scale = 0D0
      flat_bgr = 0D0
      linear_bgr = 0D0
      chi_sq = 0D0
      RETURN
   END IF

   sum_ss = SUM((diffract_calc**2D0)*errors_sq, valid)
   IF ((.NOT. refine_flat_bgr) .AND. (.NOT. refine_linear_bgr)) THEN
      IF (refine_scale) THEN
         sum_fs = SUM(diffract_calc*(diffract_expt - flat_bgr - q*linear_bgr)*errors_sq, valid)
         scale = sum_fs/sum_ss
      END IF
   ELSE IF ((refine_flat_bgr) .AND. (.NOT. refine_linear_bgr)) THEN
      sum_e = SUM(errors_sq, valid)
      sum_f = SUM((diffract_expt - q*linear_bgr)*errors_sq, valid)
      sum_s = SUM(diffract_calc*errors_sq, valid)
      IF (refine_scale) THEN
         sum_fs = SUM(diffract_calc*(diffract_expt - q*linear_bgr)*errors_sq, valid)
         scale = (sum_e*sum_fs - sum_s*sum_f)/(sum_e*sum_ss - sum_s**2D0)
      END IF
      flat_bgr = (sum_f - sum_s*scale)/sum_e
   ELSE IF ((.NOT. refine_flat_bgr) .AND. (refine_linear_bgr)) THEN
      sum_qq = SUM(errors_sq*q**2D0, valid)
      sum_qf = SUM((diffract_expt - flat_bgr)*errors_sq*q, valid)
      sum_qs = SUM(diffract_calc*errors_sq*q, valid)
      IF (refine_scale) THEN
         sum_fs = SUM(diffract_calc*(diffract_expt - flat_bgr)*errors_sq, valid)
         scale = (sum_qq*sum_fs - sum_qs*sum_qf)/(sum_qq*sum_ss - sum_qs**2D0)
      END IF
      linear_bgr = (sum_qf - sum_qs*scale)/sum_qq
   ELSE IF ((refine_flat_bgr) .AND. (refine_linear_bgr)) THEN
      sum_e = SUM(errors_sq, valid)
      sum_f = SUM(diffract_expt*errors_sq, valid)
      sum_s = SUM(diffract_calc*errors_sq, valid)
      sum_qq = SUM(errors_sq*q**2D0, valid)
      sum_q = SUM(errors_sq*q, valid)
      sum_qf = SUM(diffract_expt*errors_sq*q, valid)
      sum_qs = SUM(diffract_calc*errors_sq*q, valid)
      IF (refine_scale) THEN
         sum_fs = SUM(diffract_calc*diffract_expt*errors_sq, valid)
         num = (sum_e*sum_fs - sum_f*sum_s)*(sum_qq*sum_e - sum_q**2D0) - (sum_e*sum_qf - sum_q*sum_f)*(sum_qs*sum_e - sum_s*sum_q)
         denom = (sum_ss*sum_e - sum_s**2D0)*(sum_qq*sum_e - sum_q**2D0) - (sum_qs*sum_e - sum_q*sum_s)**2D0
         scale = num/denom
      END IF
      linear_bgr = ((sum_e*sum_qf - sum_q*sum_f) - scale*(sum_qs*sum_e - sum_q*sum_s))/(sum_qq*sum_e - sum_q**2D0)
      flat_bgr = (sum_f - scale*sum_s - linear_bgr*sum_q)/sum_e
   ELSE
      STOP 'CALC_CHI_SQ error: Unrecognised option'
   END IF

   chi_sq = SUM(errors_sq*(scale*diffract_calc + flat_bgr + q*linear_bgr - diffract_expt)**2D0, valid)

END SUBROUTINE calc_chi_sq

SUBROUTINE index_items(line, num_items, item_index)

   INTEGER, INTENT(in) :: num_items
   CHARACTER, INTENT(in) :: line*(*)
   INTEGER, INTENT(out) :: item_index(num_items + 1)
   LOGICAL back
   INTEGER length

   back = .TRUE.
   length = LEN_TRIM(line)
   k = INDEX(line(1:length), ' ', back)
   IF (k == 0) THEN
      nitems = 0
      RETURN
   END IF

   nitems = 1
   item_index(num_items + 1) = length + 1
   item_index(num_items - nitems + 1) = k
   DO
      ! starting with the right most blank space,
      ! look for the next non-space character down
      ! indicating there is another item in the line
      DO
         IF (k <= 0) EXIT

         IF (line(k:k) == ' ') THEN
            k = k - 1
            CYCLE
         ELSE
            nitems = nitems + 1
            EXIT
         END IF

      END DO

      ! once a non-space character is found,
      ! skip all adjacent non-space character
      DO
         IF (k <= 0) EXIT

         IF (line(k:k) /= ' ') THEN
            k = k - 1
            CYCLE
         END IF

         EXIT

      END DO

      IF (k <= 0) EXIT
      item_index(num_items - nitems + 1) = k

   END DO
END SUBROUTINE index_items

! number of space-separated items in a line
INTEGER FUNCTION nitems(line)
   CHARACTER line*(*)
   LOGICAL back
   INTEGER length

   back = .TRUE.
   length = LEN_TRIM(line)
   k = INDEX(line(1:length), ' ', back)
   IF (k == 0) THEN
      nitems = 0
      RETURN
   END IF

   nitems = 1
   DO
      ! starting with the right most blank space,
      ! look for the next non-space character down
      ! indicating there is another item in the line
      DO
         IF (k <= 0) EXIT

         IF (line(k:k) == ' ') THEN
            k = k - 1
            CYCLE
         ELSE
            nitems = nitems + 1
            EXIT
         END IF

      END DO

      ! once a non-space character is found,
      ! skip all adjacent non-space character
      DO
         IF (k <= 0) EXIT

         IF (line(k:k) /= ' ') THEN
            k = k - 1
            CYCLE
         END IF

         EXIT

      END DO

      IF (k <= 0) EXIT

   END DO
END FUNCTION nitems

SUBROUTINE calc_metric_tensor(cell_params, cell_angles_deg, r_orth, q_orth, r_metric, q_metric, vol)

   ! Calculates metric tensor

   IMPLICIT NONE
   REAL(8), INTENT(in) :: cell_params(3), cell_angles_deg(3)
   REAL(8), INTENT(out) :: r_metric(3, 3), q_metric(3, 3), q_orth(3, 3), r_orth(3, 3)
   REAL(8) :: pi, cos_alpha_rec, r_orth_tmp(3, 3), cell_angles(3), vol

   pi = ACOS(-1D0)
   cell_angles = cell_angles_deg*pi/180D0 ! convert to radians
   cos_alpha_rec = COS(cell_angles(2))*COS(cell_angles(3)) - COS(cell_angles(1))
   cos_alpha_rec = cos_alpha_rec/(SIN(cell_angles(2))*SIN(cell_angles(3)))
   vol = PRODUCT(cell_params)*SQRT(1D0 - (COS(cell_angles(1)))**2D0 - (COS(cell_angles(2)))**2D0 &
   & - (COS(cell_angles(3)))**2D0 + 2D0*COS(cell_angles(1))*COS(cell_angles(2))*COS(cell_angles(3)))

   ! find real-space orthogonalisation matrix
   r_orth = 0D0
   r_orth(1, 1) = cell_params(1)
   r_orth(1, 2) = cell_params(2)*COS(cell_angles(3))
   r_orth(1, 3) = cell_params(3)*COS(cell_angles(2))
   r_orth(2, 2) = cell_params(2)*SIN(cell_angles(3))
   r_orth(2, 3) = -cell_params(3)*SIN(cell_angles(2))*cos_alpha_rec
   r_orth(3, 3) = vol/(cell_params(1)*cell_params(2)*SIN(cell_angles(3)))

   IF (isnan(vol)) STOP 'CALC_METRIC_TENSOR error: invalid cell parameters'

   r_orth_tmp = r_orth
   ! and reciprocal-space orthogonalisation matrix
   CALL inverse(r_orth_tmp, q_orth, 3)
   q_orth = TRANSPOSE(q_orth)
   q_orth = 2D0*pi*q_orth ! Note factor of 2*pi here

   ! find metric tensors
   r_metric = MATMUL(TRANSPOSE(r_orth), r_orth)
   q_metric = MATMUL(TRANSPOSE(q_orth), q_orth)

END SUBROUTINE calc_metric_tensor

SUBROUTINE cross_product(a, b, c)

   IMPLICIT NONE
   REAL(8), INTENT(in) :: a(3)
   REAL(8), INTENT(in) :: b(3)
   REAL(8), INTENT(out) :: c(3)

   c(1) = a(2)*b(3) - a(3)*b(2)
   c(2) = a(3)*b(1) - a(1)*b(3)
   c(3) = a(1)*b(2) - a(2)*b(1)

END SUBROUTINE cross_product

SUBROUTINE calc_max_radius(orth_matrix, n_vec, max_radius)

   IMPLICIT NONE
   INTEGER, INTENT(in) ::  n_vec(3)
   REAL(8), INTENT(in) :: orth_matrix(3, 3)
   REAL(8), INTENT(out) :: max_radius
   REAL(8) :: max_vec(3), cross(3), max_x(3), max_y(3), max_z(3)

   max_x(:) = MATMUL(orth_matrix, [1D0*n_vec(1), 0D0, 0D0])
   max_y(:) = MATMUL(orth_matrix, [0D0, 1D0*n_vec(2), 0D0])
   max_z(:) = MATMUL(orth_matrix, [0D0, 0D0, 1D0*n_vec(3)])
   CALL cross_product(max_y, max_z, cross)
   cross = cross/SQRT(DOT_PRODUCT(cross, cross))
   max_vec(1) = DOT_PRODUCT(max_x, cross)
   CALL cross_product(max_z, max_x, cross)
   cross = cross/SQRT(DOT_PRODUCT(cross, cross))
   max_vec(2) = DOT_PRODUCT(max_y, cross)
   CALL cross_product(max_x, max_y, cross)
   cross = cross/SQRT(DOT_PRODUCT(cross, cross))
   max_vec(3) = DOT_PRODUCT(max_z, cross)

   max_radius = MINVAL(max_vec)

END SUBROUTINE

SUBROUTINE inverse(a, c, n)

! Routine from
! "Introductory Computational Physics" by Andi Klein and Alexander Godunov (CUP, ISBN 0-521-82862-7 (2006))

!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed
! during the calculation
!===========================================================
   IMPLICIT NONE

   INTEGER, INTENT(in) :: n
   REAL(8), INTENT(inout) :: a(n, n)
   REAL(8), INTENT(out) :: c(n, n)
   REAL(8) :: L(n, n), U(n, n), b(n), d(n), x(n), coeff
   INTEGER i, j, k

   ! step 0: initialization for matrices L and U and b
   ! Fortran 90/95 aloows such operations on matrices
   L = 0D0
   U = 0D0
   b = 0D0

   ! step 1: forward elimination
   DO k = 1, n - 1
      DO i = k + 1, n
         coeff = a(i, k)/a(k, k)
         L(i, k) = coeff
         DO j = k + 1, n
            a(i, j) = a(i, j) - coeff*a(k, j)
         END DO
      END DO
   END DO

   ! Step 2: prepare L and U matrices
   ! L matrix is a matrix of the elimination coefficient
   ! + the diagonal elements are 1.0
   DO i = 1, n
      L(i, i) = 1D0
   END DO
   ! U matrix is the upper triangular part of A
   DO j = 1, n
      DO i = 1, j
         U(i, j) = a(i, j)
      END DO
   END DO

   ! Step 3: compute columns of the inverse matrix C
   DO k = 1, n
      b(k) = 1D0
      d(1) = b(1)
      ! Step 3a: Solve Ld=b using the forward substitution
      DO i = 2, n
         d(i) = b(i)
         DO j = 1, i - 1
            d(i) = d(i) - L(i, j)*d(j)
         END DO
      END DO
      ! Step 3b: Solve Ux=d using the back substitution
      x(n) = d(n)/U(n, n)
      DO i = n - 1, 1, -1
         x(i) = d(i)
         DO j = n, i + 1, -1
            x(i) = x(i) - U(i, j)*x(j)
         END DO
         x(i) = x(i)/u(i, i)
      END DO
      ! Step 3c: fill the solutions x(n) into column k of C
      DO i = 1, n
         c(i, k) = x(i)
      END DO
      b(k) = 0D0
   END DO

END SUBROUTINE inverse

FUNCTION itc(i)

   INTEGER :: i
   CHARACTER(4) :: itc
   CHARACTER(8) :: fmt

   IF (i .GE. 0 .AND. i < 100) THEN
      fmt = '(I2.2)' ! an integer of width 2 with zeros at the left
      WRITE (itc, fmt) i ! converting integer to string using 'internal file'
   ELSE IF (i .GE. 100 .AND. i < 1000) THEN
      fmt = '(I3.3)' ! an integer of width 3 with zeros at the left
      WRITE (itc, fmt) i ! converting integer to string using 'internal file'
   ELSE IF (i .GE. 1000 .AND. i < 10000) THEN
      fmt = '(I4.4)' ! an integer of width 4 with zeros at the left
      WRITE (itc, fmt) i ! converting integer to string using 'internal file'
   ELSE
      STOP 'ITC error - cannot have int > 9999.'
   END IF

END FUNCTION itc

SUBROUTINE spinsym_int(laue_class, isym, q, q_new)

   IMPLICIT NONE
   INTEGER, INTENT(in) :: q(3)
   INTEGER, INTENT(out) :: q_new(3)
   INTEGER, INTENT(in) :: isym, laue_class
   INTEGER :: h, k, l, i, h_new, k_new, l_new

   h = q(1)
   k = q(2)
   l = q(3)

   SELECT CASE (laue_class)

   CASE (1) ! -1
      IF (isym .EQ. 1) THEN
         h_new = h
         k_new = k
         l_new = l
      ELSEIF (isym .EQ. 2) THEN
         h_new = -h
         k_new = -k
         l_new = -l
      ELSE
         STOP 'SPINSYM error: incorrect number of symmetry operations'
      END IF

   CASE (2) ! 2/m (b axis unique)
      IF (isym .EQ. 1) THEN
         h_new = h
         k_new = k
         l_new = l
      ELSEIF (isym .EQ. 2) THEN
         h_new = -h
         k_new = k
         l_new = -l
      ELSEIF (isym .EQ. 3) THEN
         h_new = -h
         k_new = -k
         l_new = -l
      ELSEIF (isym .EQ. 4) THEN
         h_new = h
         k_new = -k
         l_new = l
      ELSE
         STOP 'SPINSYM error: incorrect number of symmetry operations'
      END IF

   CASE (3) ! mmm
      IF (isym .EQ. 1) THEN
         h_new = h
         k_new = k
         l_new = l
      ELSEIF (isym .EQ. 2) THEN
         h_new = -h
         k_new = -k
         l_new = l
      ELSEIF (isym .EQ. 3) THEN
         h_new = -h
         k_new = k
         l_new = -l
      ELSEIF (isym .EQ. 4) THEN
         h_new = h
         k_new = -k
         l_new = -l
      ELSEIF (isym .EQ. 5) THEN
         h_new = -h
         k_new = -k
         l_new = -l
      ELSEIF (isym .EQ. 6) THEN
         h_new = h
         k_new = k
         l_new = -l
      ELSEIF (isym .EQ. 7) THEN
         h_new = h
         k_new = -k
         l_new = l
      ELSEIF (isym .EQ. 8) THEN
         h_new = -h
         k_new = k
         l_new = l
      ELSE
         STOP 'SPINSYM error: incorrect number of symmetry operations'
      END IF

   CASE (4) ! 4/m
      IF (isym .EQ. 1) THEN
         h_new = h
         k_new = k
         l_new = l
      ELSEIF (isym .EQ. 2) THEN
         h_new = -h
         k_new = -k
         l_new = l
      ELSEIF (isym .EQ. 3) THEN
         h_new = -k
         k_new = h
         l_new = -l
      ELSEIF (isym .EQ. 4) THEN
         h_new = k
         k_new = -h
         l_new = -l
      ELSEIF (isym .EQ. 5) THEN
         h_new = -h
         k_new = -k
         l_new = -l
      ELSEIF (isym .EQ. 6) THEN
         h_new = h
         k_new = k
         l_new = -l
      ELSEIF (isym .EQ. 7) THEN
         h_new = k
         k_new = -h
         l_new = l
      ELSEIF (isym .EQ. 8) THEN
         h_new = -k
         k_new = h
         l_new = l
      ELSE
         STOP 'SPINSYM error: incorrect number of symmetry operations'
      END IF

   CASE (5) ! 4/mmm
      IF (isym .EQ. 1) THEN
         h_new = h
         k_new = k
         l_new = l
      ELSEIF (isym .EQ. 2) THEN
         h_new = -h
         k_new = -k
         l_new = l
      ELSEIF (isym .EQ. 3) THEN
         h_new = k
         k_new = -h
         l_new = l
      ELSEIF (isym .EQ. 4) THEN
         h_new = -k
         k_new = h
         l_new = l
      ELSEIF (isym .EQ. 5) THEN
         h_new = -h
         k_new = k
         l_new = -l
      ELSEIF (isym .EQ. 6) THEN
         h_new = h
         k_new = -k
         l_new = -l
      ELSEIF (isym .EQ. 7) THEN
         h_new = k
         k_new = h
         l_new = -l
      ELSEIF (isym .EQ. 8) THEN
         h_new = -k
         k_new = -h
         l_new = -l
      ELSEIF (isym .EQ. 9) THEN
         h_new = -h
         k_new = -k
         l_new = -l
      ELSEIF (isym .EQ. 10) THEN
         h_new = h
         k_new = k
         l_new = -l
      ELSEIF (isym .EQ. 11) THEN
         h_new = -k
         k_new = h
         l_new = -l
      ELSEIF (isym .EQ. 12) THEN
         h_new = k
         k_new = -h
         l_new = -l
      ELSEIF (isym .EQ. 13) THEN
         h_new = h
         k_new = -k
         l_new = l
      ELSEIF (isym .EQ. 14) THEN
         h_new = -h
         k_new = k
         l_new = l
      ELSEIF (isym .EQ. 15) THEN
         h_new = -k
         k_new = -h
         l_new = l
      ELSEIF (isym .EQ. 16) THEN
         h_new = k
         k_new = h
         l_new = l
      ELSE
         STOP 'SPINSYM error: incorrect number of symmetry operations'
      END IF

   CASE (6) ! -3; note HEXAGONAL AXES
      i = -h - k
      IF (isym .EQ. 1) THEN
         h_new = h
         k_new = k
         l_new = l
      ELSEIF (isym .EQ. 2) THEN
         h_new = k
         k_new = i
         l_new = l
      ELSEIF (isym .EQ. 3) THEN
         h_new = i
         k_new = h
         l_new = l
      ELSEIF (isym .EQ. 4) THEN
         h_new = -h
         k_new = -k
         l_new = -l
      ELSEIF (isym .EQ. 5) THEN
         h_new = -k
         k_new = -i
         l_new = -l
      ELSEIF (isym .EQ. 6) THEN
         h_new = -i
         k_new = -h
         l_new = -l
      ELSE
         STOP 'SPINSYM error: incorrect number of symmetry operations'
      END IF

   CASE (7) ! -31m ; note HEXAGONAL AXES
      i = -h - k
      IF (isym .EQ. 1) THEN
         h_new = h
         k_new = k
         l_new = l
      ELSEIF (isym .EQ. 2) THEN
         h_new = k
         k_new = i
         l_new = l
      ELSEIF (isym .EQ. 3) THEN
         h_new = i
         k_new = h
         l_new = l
      ELSEIF (isym .EQ. 4) THEN
         h_new = -k
         k_new = -h
         l_new = -l
      ELSEIF (isym .EQ. 5) THEN
         h_new = -h
         k_new = -i
         l_new = -l
      ELSEIF (isym .EQ. 6) THEN
         h_new = -i
         k_new = -k
         l_new = -l
      ELSEIF (isym .EQ. 7) THEN
         h_new = -h
         k_new = -k
         l_new = -l
      ELSEIF (isym .EQ. 8) THEN
         h_new = -k
         k_new = -i
         l_new = -l
      ELSEIF (isym .EQ. 9) THEN
         h_new = -i
         k_new = -h
         l_new = -l
      ELSEIF (isym .EQ. 10) THEN
         h_new = k
         k_new = h
         l_new = l
      ELSEIF (isym .EQ. 11) THEN
         h_new = h
         k_new = i
         l_new = l
      ELSEIF (isym .EQ. 12) THEN
         h_new = i
         k_new = k
         l_new = l
      ELSE
         STOP 'SPINSYM error: incorrect number of symmetry operations'
      END IF

   CASE (8) ! -3m1 or -3m ; note HEXAGONAL AXES
      i = -h - k
      IF (isym .EQ. 1) THEN
         h_new = h
         k_new = k
         l_new = l
      ELSEIF (isym .EQ. 2) THEN
         h_new = k
         k_new = i
         l_new = l
      ELSEIF (isym .EQ. 3) THEN
         h_new = i
         k_new = h
         l_new = l
      ELSEIF (isym .EQ. 4) THEN
         h_new = -k
         k_new = -h
         l_new = l
      ELSEIF (isym .EQ. 5) THEN
         h_new = -h
         k_new = -i
         l_new = l
      ELSEIF (isym .EQ. 6) THEN
         h_new = -i
         k_new = -k
         l_new = l
      ELSEIF (isym .EQ. 7) THEN
         h_new = -h
         k_new = -k
         l_new = -l
      ELSEIF (isym .EQ. 8) THEN
         h_new = -k
         k_new = -i
         l_new = -l
      ELSEIF (isym .EQ. 9) THEN
         h_new = -i
         k_new = -h
         l_new = -l
      ELSEIF (isym .EQ. 10) THEN
         h_new = k
         k_new = h
         l_new = -l
      ELSEIF (isym .EQ. 11) THEN
         h_new = h
         k_new = i
         l_new = -l
      ELSEIF (isym .EQ. 12) THEN
         h_new = i
         k_new = k
         l_new = -l
      ELSE
         STOP 'SPINSYM error: incorrect number of symmetry operations'
      END IF

   CASE (9) ! 6/m
      i = -h - k
      IF (isym .EQ. 1) THEN
         h_new = h
         k_new = k
         l_new = l
      ELSEIF (isym .EQ. 2) THEN
         h_new = k
         k_new = i
         l_new = l
      ELSEIF (isym .EQ. 3) THEN
         h_new = i
         k_new = h
         l_new = l
      ELSEIF (isym .EQ. 4) THEN
         h_new = -h
         k_new = -k
         l_new = l
      ELSEIF (isym .EQ. 5) THEN
         h_new = -k
         k_new = -i
         l_new = l
      ELSEIF (isym .EQ. 6) THEN
         h_new = -i
         k_new = -h
         l_new = l
      ELSEIF (isym .EQ. 7) THEN
         h_new = -h
         k_new = -k
         l_new = -l
      ELSEIF (isym .EQ. 8) THEN
         h_new = -k
         k_new = -i
         l_new = -l
      ELSEIF (isym .EQ. 9) THEN
         h_new = -i
         k_new = -h
         l_new = -l
      ELSEIF (isym .EQ. 10) THEN
         h_new = h
         k_new = k
         l_new = -l
      ELSEIF (isym .EQ. 11) THEN
         h_new = k
         k_new = i
         l_new = -l
      ELSEIF (isym .EQ. 12) THEN
         h_new = i
         k_new = h
         l_new = -l
      ELSE
         STOP 'SPINSYM error: incorrect number of symmetry operations'
      END IF

   CASE (10) ! 6/mmm
      i = -h - k
      IF (isym .EQ. 1) THEN
         h_new = h
         k_new = k
         l_new = l
      ELSEIF (isym .EQ. 2) THEN
         h_new = k
         k_new = i
         l_new = l
      ELSEIF (isym .EQ. 3) THEN
         h_new = i
         k_new = h
         l_new = l
      ELSEIF (isym .EQ. 4) THEN
         h_new = -h
         k_new = -k
         l_new = l
      ELSEIF (isym .EQ. 5) THEN
         h_new = -k
         k_new = -i
         l_new = l
      ELSEIF (isym .EQ. 6) THEN
         h_new = -i
         k_new = -h
         l_new = l
      ELSEIF (isym .EQ. 7) THEN
         h_new = k
         k_new = h
         l_new = -l
      ELSEIF (isym .EQ. 8) THEN
         h_new = h
         k_new = i
         l_new = -l
      ELSEIF (isym .EQ. 9) THEN
         h_new = i
         k_new = k
         l_new = -l
      ELSEIF (isym .EQ. 10) THEN
         h_new = -k
         k_new = -h
         l_new = -l
      ELSEIF (isym .EQ. 11) THEN
         h_new = -h
         k_new = -i
         l_new = -l
      ELSEIF (isym .EQ. 12) THEN
         h_new = -i
         k_new = -k
         l_new = -l
      ELSEIF (isym .EQ. 13) THEN
         h_new = -h
         k_new = -k
         l_new = -l
      ELSEIF (isym .EQ. 14) THEN
         h_new = -k
         k_new = -i
         l_new = -l
      ELSEIF (isym .EQ. 15) THEN
         h_new = -i
         k_new = -h
         l_new = -l
      ELSEIF (isym .EQ. 16) THEN
         h_new = h
         k_new = k
         l_new = -l
      ELSEIF (isym .EQ. 17) THEN
         h_new = k
         k_new = i
         l_new = -l
      ELSEIF (isym .EQ. 18) THEN
         h_new = i
         k_new = h
         l_new = -l
      ELSEIF (isym .EQ. 19) THEN
         h_new = -k
         k_new = -h
         l_new = l
      ELSEIF (isym .EQ. 20) THEN
         h_new = -h
         k_new = -i
         l_new = l
      ELSEIF (isym .EQ. 21) THEN
         h_new = -i
         k_new = -k
         l_new = l
      ELSEIF (isym .EQ. 22) THEN
         h_new = k
         k_new = h
         l_new = l
      ELSEIF (isym .EQ. 23) THEN
         h_new = h
         k_new = i
         l_new = l
      ELSEIF (isym .EQ. 24) THEN
         h_new = i
         k_new = k
         l_new = l
      ELSE
         STOP 'SPINSYM error: incorrect number of symmetry operations'
      END IF

   CASE (11)
      IF (isym .EQ. 1) THEN
         h_new = h
         k_new = k
         l_new = l
      ELSEIF (isym .EQ. 2) THEN
         h_new = -h
         k_new = k
         l_new = l
      ELSEIF (isym .EQ. 3) THEN
         h_new = h
         k_new = -k
         l_new = l
      ELSEIF (isym .EQ. 4) THEN
         h_new = h
         k_new = k
         l_new = -l
      ELSEIF (isym .EQ. 5) THEN
         h_new = -h
         k_new = -k
         l_new = -l
      ELSEIF (isym .EQ. 6) THEN
         h_new = h
         k_new = -k
         l_new = -l
      ELSEIF (isym .EQ. 7) THEN
         h_new = -h
         k_new = k
         l_new = -l
      ELSEIF (isym .EQ. 8) THEN
         h_new = -h
         k_new = -k
         l_new = l
      ELSEIF (isym .EQ. 9) THEN
         h_new = k
         k_new = l
         l_new = h
      ELSEIF (isym .EQ. 10) THEN
         h_new = -k
         k_new = l
         l_new = h
      ELSEIF (isym .EQ. 11) THEN
         h_new = k
         k_new = -l
         l_new = h
      ELSEIF (isym .EQ. 12) THEN
         h_new = k
         k_new = l
         l_new = -h
      ELSEIF (isym .EQ. 13) THEN
         h_new = -k
         k_new = -l
         l_new = -h
      ELSEIF (isym .EQ. 14) THEN
         h_new = k
         k_new = -l
         l_new = -h
      ELSEIF (isym .EQ. 15) THEN
         h_new = -k
         k_new = l
         l_new = -h
      ELSEIF (isym .EQ. 16) THEN
         h_new = -k
         k_new = -l
         l_new = h
      ELSEIF (isym .EQ. 17) THEN
         h_new = l
         k_new = h
         l_new = k
      ELSEIF (isym .EQ. 18) THEN
         h_new = -l
         k_new = h
         l_new = k
      ELSEIF (isym .EQ. 19) THEN
         h_new = l
         k_new = -h
         l_new = k
      ELSEIF (isym .EQ. 20) THEN
         h_new = l
         k_new = h
         l_new = -k
      ELSEIF (isym .EQ. 21) THEN
         h_new = -l
         k_new = -h
         l_new = -k
      ELSEIF (isym .EQ. 22) THEN
         h_new = l
         k_new = -h
         l_new = -k
      ELSEIF (isym .EQ. 23) THEN
         h_new = -l
         k_new = h
         l_new = -k
      ELSEIF (isym .EQ. 24) THEN
         h_new = -l
         k_new = -h
         l_new = k
      ELSE
         STOP 'SPINSYM error: incorrect number of symmetry operations'
      END IF

   CASE (12)

      IF (isym .EQ. 1) THEN
         h_new = h
         k_new = k
         l_new = l
      ELSEIF (isym .EQ. 2) THEN
         h_new = h
         k_new = k
         l_new = -l
      ELSEIF (isym .EQ. 3) THEN
         h_new = h
         k_new = -k
         l_new = l
      ELSEIF (isym .EQ. 4) THEN
         h_new = -h
         k_new = k
         l_new = l
      ELSEIF (isym .EQ. 5) THEN
         h_new = -h
         k_new = -k
         l_new = -l
      ELSEIF (isym .EQ. 6) THEN
         h_new = -h
         k_new = -k
         l_new = l
      ELSEIF (isym .EQ. 7) THEN
         h_new = -h
         k_new = k
         l_new = -l
      ELSEIF (isym .EQ. 8) THEN
         h_new = h
         k_new = -k
         l_new = -l
      ELSEIF (isym .EQ. 9) THEN
         h_new = k
         k_new = l
         l_new = h
      ELSEIF (isym .EQ. 10) THEN
         h_new = k
         k_new = l
         l_new = -h
      ELSEIF (isym .EQ. 11) THEN
         h_new = k
         k_new = -l
         l_new = h
      ELSEIF (isym .EQ. 12) THEN
         h_new = -k
         k_new = l
         l_new = h
      ELSEIF (isym .EQ. 13) THEN
         h_new = -k
         k_new = -l
         l_new = -h
      ELSEIF (isym .EQ. 14) THEN
         h_new = -k
         k_new = -l
         l_new = h
      ELSEIF (isym .EQ. 15) THEN
         h_new = -k
         k_new = l
         l_new = -h
      ELSEIF (isym .EQ. 16) THEN
         h_new = k
         k_new = -l
         l_new = -h
      ELSEIF (isym .EQ. 17) THEN
         h_new = l
         k_new = h
         l_new = k
      ELSEIF (isym .EQ. 18) THEN
         h_new = l
         k_new = h
         l_new = -k
      ELSEIF (isym .EQ. 19) THEN
         h_new = l
         k_new = -h
         l_new = k
      ELSEIF (isym .EQ. 20) THEN
         h_new = -l
         k_new = h
         l_new = k
      ELSEIF (isym .EQ. 21) THEN
         h_new = -l
         k_new = -h
         l_new = -k
      ELSEIF (isym .EQ. 22) THEN
         h_new = -l
         k_new = -h
         l_new = k
      ELSEIF (isym .EQ. 23) THEN
         h_new = -l
         k_new = h
         l_new = -k
      ELSEIF (isym .EQ. 24) THEN
         h_new = l
         k_new = -h
         l_new = -k
      ELSEIF (isym .EQ. 25) THEN
         h_new = h
         k_new = l
         l_new = k
      ELSEIF (isym .EQ. 26) THEN
         h_new = h
         k_new = l
         l_new = -k
      ELSEIF (isym .EQ. 27) THEN
         h_new = h
         k_new = -l
         l_new = k
      ELSEIF (isym .EQ. 28) THEN
         h_new = -h
         k_new = l
         l_new = k
      ELSEIF (isym .EQ. 29) THEN
         h_new = -h
         k_new = -l
         l_new = -k
      ELSEIF (isym .EQ. 30) THEN
         h_new = -h
         k_new = -l
         l_new = k
      ELSEIF (isym .EQ. 31) THEN
         h_new = -h
         k_new = l
         l_new = -k
      ELSEIF (isym .EQ. 32) THEN
         h_new = h
         k_new = -l
         l_new = -k
      ELSEIF (isym .EQ. 33) THEN
         h_new = k
         k_new = h
         l_new = l
      ELSEIF (isym .EQ. 34) THEN
         h_new = k
         k_new = h
         l_new = -l
      ELSEIF (isym .EQ. 35) THEN
         h_new = k
         k_new = -h
         l_new = l
      ELSEIF (isym .EQ. 36) THEN
         h_new = -k
         k_new = h
         l_new = l
      ELSEIF (isym .EQ. 37) THEN
         h_new = -k
         k_new = -h
         l_new = -l
      ELSEIF (isym .EQ. 38) THEN
         h_new = -k
         k_new = -h
         l_new = l
      ELSEIF (isym .EQ. 39) THEN
         h_new = -k
         k_new = h
         l_new = -l
      ELSEIF (isym .EQ. 40) THEN
         h_new = k
         k_new = -h
         l_new = -l
      ELSEIF (isym .EQ. 41) THEN
         h_new = l
         k_new = k
         l_new = h
      ELSEIF (isym .EQ. 42) THEN
         h_new = l
         k_new = k
         l_new = -h
      ELSEIF (isym .EQ. 43) THEN
         h_new = l
         k_new = -k
         l_new = h
      ELSEIF (isym .EQ. 44) THEN
         h_new = -l
         k_new = k
         l_new = h
      ELSEIF (isym .EQ. 45) THEN
         h_new = -l
         k_new = -k
         l_new = -h
      ELSEIF (isym .EQ. 46) THEN
         h_new = -l
         k_new = -k
         l_new = h
      ELSEIF (isym .EQ. 47) THEN
         h_new = -l
         k_new = k
         l_new = -h
      ELSEIF (isym .EQ. 48) THEN
         h_new = l
         k_new = -k
         l_new = -h
      ELSE
         STOP 'SPINSYM error: incorrect number of symmetry operations'
      END IF

   CASE default
      STOP 'SPINSYM error: incorrect value for Laue class'

   END SELECT

   q_new(1) = h_new
   q_new(2) = k_new
   q_new(3) = l_new

END SUBROUTINE

SUBROUTINE write_ppm(intensity_img, width_ppm, height_ppm, colourmap, filename_stem, max, min)
   IMPLICIT NONE

   INTEGER, INTENT(in) :: width_ppm, height_ppm, colourmap
   REAL, INTENT(in) :: intensity_img(width_ppm, height_ppm)
   REAL, INTENT(inout) :: max, min
   CHARACTER(128), INTENT(in) :: filename_stem
   INTEGER :: i, x, y, rgb(3, width_ppm), width2, height2
   REAL :: rowintensity(width_ppm)
   REAL, ALLOCATABLE :: intensity_img2(:, :)
   REAL, PARAMETER :: eps = 10.0*EPSILON(0.0)

   IF ((ABS(max + 1) < eps) .AND. (ABS(min + 1) < eps)) THEN
      max = MAXVAL(intensity_img)
      min = MINVAL(intensity_img)
   END IF

   OPEN (11, file=TRIM(ADJUSTL(filename_stem))//'_ppm_scale.txt', status='replace')
   WRITE (11, *) 'Maximum intensity (for scale bar): ', max
   WRITE (11, *) 'Minimum intensity (for scale bar): ', min
   CLOSE (11)

   DO i = 11, 12
      IF (i == 11) THEN
         OPEN (i, file=TRIM(ADJUSTL(filename_stem))//'.ppm', status='replace')
         width2 = width_ppm
         height2 = height_ppm
         ALLOCATE (intensity_img2(width2, height2))
         intensity_img2 = intensity_img
      ELSE
         OPEN (i, file=TRIM(ADJUSTL(filename_stem))//'_colourbar.ppm', status='replace')
         width2 = 20
         height2 = 256
         DEALLOCATE (intensity_img2)
         ALLOCATE (intensity_img2(width2, height2))
         DO y = 1, height2
            intensity_img2(:, y) = min + (max - min)*(y - 1)/255.0
         END DO
      END IF
      WRITE (i, '(A2,1X,2(I5,1X),A3)') 'P3', width2, height2, '255'
      DO y = height2, 1, -1
         DO x = 1, width2
         IF (ABS(max - min) < eps) THEN
            rowintensity(x) = 0.0
         ELSE
            rowintensity(x) = INT(255*(intensity_img2(x, y) - min)/(max - min))
         END IF
         IF (rowintensity(x) < 0) rowintensity(x) = 0
         IF (rowintensity(x) > 255) rowintensity(x) = 255
         SELECT CASE (colourmap)
         CASE (1)
            IF (rowintensity(x) == 0) rgb(:, x) = [59, 76, 192]
            IF (rowintensity(x) == 1) rgb(:, x) = [60, 78, 194]
            IF (rowintensity(x) == 2) rgb(:, x) = [61, 80, 195]
            IF (rowintensity(x) == 3) rgb(:, x) = [62, 81, 197]
            IF (rowintensity(x) == 4) rgb(:, x) = [63, 83, 198]
            IF (rowintensity(x) == 5) rgb(:, x) = [64, 85, 200]
            IF (rowintensity(x) == 6) rgb(:, x) = [66, 87, 201]
            IF (rowintensity(x) == 7) rgb(:, x) = [67, 88, 203]
            IF (rowintensity(x) == 8) rgb(:, x) = [68, 90, 204]
            IF (rowintensity(x) == 9) rgb(:, x) = [69, 92, 206]
            IF (rowintensity(x) == 10) rgb(:, x) = [70, 93, 207]
            IF (rowintensity(x) == 11) rgb(:, x) = [71, 95, 209]
            IF (rowintensity(x) == 12) rgb(:, x) = [73, 97, 210]
            IF (rowintensity(x) == 13) rgb(:, x) = [74, 99, 211]
            IF (rowintensity(x) == 14) rgb(:, x) = [75, 100, 213]
            IF (rowintensity(x) == 15) rgb(:, x) = [76, 102, 214]
            IF (rowintensity(x) == 16) rgb(:, x) = [77, 104, 215]
            IF (rowintensity(x) == 17) rgb(:, x) = [79, 105, 217]
            IF (rowintensity(x) == 18) rgb(:, x) = [80, 107, 218]
            IF (rowintensity(x) == 19) rgb(:, x) = [81, 109, 219]
            IF (rowintensity(x) == 20) rgb(:, x) = [82, 110, 221]
            IF (rowintensity(x) == 21) rgb(:, x) = [84, 112, 222]
            IF (rowintensity(x) == 22) rgb(:, x) = [85, 114, 223]
            IF (rowintensity(x) == 23) rgb(:, x) = [86, 115, 224]
            IF (rowintensity(x) == 24) rgb(:, x) = [87, 117, 225]
            IF (rowintensity(x) == 25) rgb(:, x) = [89, 119, 226]
            IF (rowintensity(x) == 26) rgb(:, x) = [90, 120, 228]
            IF (rowintensity(x) == 27) rgb(:, x) = [91, 122, 229]
            IF (rowintensity(x) == 28) rgb(:, x) = [93, 123, 230]
            IF (rowintensity(x) == 29) rgb(:, x) = [94, 125, 231]
            IF (rowintensity(x) == 30) rgb(:, x) = [95, 127, 232]
            IF (rowintensity(x) == 31) rgb(:, x) = [96, 128, 233]
            IF (rowintensity(x) == 32) rgb(:, x) = [98, 130, 234]
            IF (rowintensity(x) == 33) rgb(:, x) = [99, 131, 235]
            IF (rowintensity(x) == 34) rgb(:, x) = [100, 133, 236]
            IF (rowintensity(x) == 35) rgb(:, x) = [102, 135, 237]
            IF (rowintensity(x) == 36) rgb(:, x) = [103, 136, 238]
            IF (rowintensity(x) == 37) rgb(:, x) = [104, 138, 239]
            IF (rowintensity(x) == 38) rgb(:, x) = [106, 139, 239]
            IF (rowintensity(x) == 39) rgb(:, x) = [107, 141, 240]
            IF (rowintensity(x) == 40) rgb(:, x) = [108, 142, 241]
            IF (rowintensity(x) == 41) rgb(:, x) = [110, 144, 242]
            IF (rowintensity(x) == 42) rgb(:, x) = [111, 145, 243]
            IF (rowintensity(x) == 43) rgb(:, x) = [112, 147, 243]
            IF (rowintensity(x) == 44) rgb(:, x) = [114, 148, 244]
            IF (rowintensity(x) == 45) rgb(:, x) = [115, 150, 245]
            IF (rowintensity(x) == 46) rgb(:, x) = [116, 151, 246]
            IF (rowintensity(x) == 47) rgb(:, x) = [118, 153, 246]
            IF (rowintensity(x) == 48) rgb(:, x) = [119, 154, 247]
            IF (rowintensity(x) == 49) rgb(:, x) = [120, 156, 247]
            IF (rowintensity(x) == 50) rgb(:, x) = [122, 157, 248]
            IF (rowintensity(x) == 51) rgb(:, x) = [123, 158, 249]
            IF (rowintensity(x) == 52) rgb(:, x) = [124, 160, 249]
            IF (rowintensity(x) == 53) rgb(:, x) = [126, 161, 250]
            IF (rowintensity(x) == 54) rgb(:, x) = [127, 163, 250]
            IF (rowintensity(x) == 55) rgb(:, x) = [129, 164, 251]
            IF (rowintensity(x) == 56) rgb(:, x) = [130, 165, 251]
            IF (rowintensity(x) == 57) rgb(:, x) = [131, 167, 252]
            IF (rowintensity(x) == 58) rgb(:, x) = [133, 168, 252]
            IF (rowintensity(x) == 59) rgb(:, x) = [134, 169, 252]
            IF (rowintensity(x) == 60) rgb(:, x) = [135, 171, 253]
            IF (rowintensity(x) == 61) rgb(:, x) = [137, 172, 253]
            IF (rowintensity(x) == 62) rgb(:, x) = [138, 173, 253]
            IF (rowintensity(x) == 63) rgb(:, x) = [140, 174, 254]
            IF (rowintensity(x) == 64) rgb(:, x) = [141, 176, 254]
            IF (rowintensity(x) == 65) rgb(:, x) = [142, 177, 254]
            IF (rowintensity(x) == 66) rgb(:, x) = [144, 178, 254]
            IF (rowintensity(x) == 67) rgb(:, x) = [145, 179, 254]
            IF (rowintensity(x) == 68) rgb(:, x) = [147, 181, 255]
            IF (rowintensity(x) == 69) rgb(:, x) = [148, 182, 255]
            IF (rowintensity(x) == 70) rgb(:, x) = [149, 183, 255]
            IF (rowintensity(x) == 71) rgb(:, x) = [151, 184, 255]
            IF (rowintensity(x) == 72) rgb(:, x) = [152, 185, 255]
            IF (rowintensity(x) == 73) rgb(:, x) = [153, 186, 255]
            IF (rowintensity(x) == 74) rgb(:, x) = [155, 187, 255]
            IF (rowintensity(x) == 75) rgb(:, x) = [156, 188, 255]
            IF (rowintensity(x) == 76) rgb(:, x) = [158, 190, 255]
            IF (rowintensity(x) == 77) rgb(:, x) = [159, 191, 255]
            IF (rowintensity(x) == 78) rgb(:, x) = [160, 192, 255]
            IF (rowintensity(x) == 79) rgb(:, x) = [162, 193, 255]
            IF (rowintensity(x) == 80) rgb(:, x) = [163, 194, 255]
            IF (rowintensity(x) == 81) rgb(:, x) = [164, 195, 254]
            IF (rowintensity(x) == 82) rgb(:, x) = [166, 196, 254]
            IF (rowintensity(x) == 83) rgb(:, x) = [167, 197, 254]
            IF (rowintensity(x) == 84) rgb(:, x) = [168, 198, 254]
            IF (rowintensity(x) == 85) rgb(:, x) = [170, 199, 253]
            IF (rowintensity(x) == 86) rgb(:, x) = [171, 199, 253]
            IF (rowintensity(x) == 87) rgb(:, x) = [172, 200, 253]
            IF (rowintensity(x) == 88) rgb(:, x) = [174, 201, 253]
            IF (rowintensity(x) == 89) rgb(:, x) = [175, 202, 252]
            IF (rowintensity(x) == 90) rgb(:, x) = [176, 203, 252]
            IF (rowintensity(x) == 91) rgb(:, x) = [178, 204, 251]
            IF (rowintensity(x) == 92) rgb(:, x) = [179, 205, 251]
            IF (rowintensity(x) == 93) rgb(:, x) = [180, 205, 251]
            IF (rowintensity(x) == 94) rgb(:, x) = [182, 206, 250]
            IF (rowintensity(x) == 95) rgb(:, x) = [183, 207, 250]
            IF (rowintensity(x) == 96) rgb(:, x) = [184, 208, 249]
            IF (rowintensity(x) == 97) rgb(:, x) = [185, 208, 248]
            IF (rowintensity(x) == 98) rgb(:, x) = [187, 209, 248]
            IF (rowintensity(x) == 99) rgb(:, x) = [188, 210, 247]
            IF (rowintensity(x) == 100) rgb(:, x) = [189, 210, 247]
            IF (rowintensity(x) == 101) rgb(:, x) = [190, 211, 246]
            IF (rowintensity(x) == 102) rgb(:, x) = [192, 212, 245]
            IF (rowintensity(x) == 103) rgb(:, x) = [193, 212, 245]
            IF (rowintensity(x) == 104) rgb(:, x) = [194, 213, 244]
            IF (rowintensity(x) == 105) rgb(:, x) = [195, 213, 243]
            IF (rowintensity(x) == 106) rgb(:, x) = [197, 214, 243]
            IF (rowintensity(x) == 107) rgb(:, x) = [198, 214, 242]
            IF (rowintensity(x) == 108) rgb(:, x) = [199, 215, 241]
            IF (rowintensity(x) == 109) rgb(:, x) = [200, 215, 240]
            IF (rowintensity(x) == 110) rgb(:, x) = [201, 216, 239]
            IF (rowintensity(x) == 111) rgb(:, x) = [203, 216, 238]
            IF (rowintensity(x) == 112) rgb(:, x) = [204, 217, 238]
            IF (rowintensity(x) == 113) rgb(:, x) = [205, 217, 237]
            IF (rowintensity(x) == 114) rgb(:, x) = [206, 217, 236]
            IF (rowintensity(x) == 115) rgb(:, x) = [207, 218, 235]
            IF (rowintensity(x) == 116) rgb(:, x) = [208, 218, 234]
            IF (rowintensity(x) == 117) rgb(:, x) = [209, 219, 233]
            IF (rowintensity(x) == 118) rgb(:, x) = [210, 219, 232]
            IF (rowintensity(x) == 119) rgb(:, x) = [211, 219, 231]
            IF (rowintensity(x) == 120) rgb(:, x) = [213, 219, 230]
            IF (rowintensity(x) == 121) rgb(:, x) = [214, 220, 229]
            IF (rowintensity(x) == 122) rgb(:, x) = [215, 220, 228]
            IF (rowintensity(x) == 123) rgb(:, x) = [216, 220, 227]
            IF (rowintensity(x) == 124) rgb(:, x) = [217, 220, 225]
            IF (rowintensity(x) == 125) rgb(:, x) = [218, 220, 224]
            IF (rowintensity(x) == 126) rgb(:, x) = [219, 220, 223]
            IF (rowintensity(x) == 127) rgb(:, x) = [220, 221, 222]
            IF (rowintensity(x) == 128) rgb(:, x) = [221, 221, 221]
            IF (rowintensity(x) == 129) rgb(:, x) = [222, 220, 219]
            IF (rowintensity(x) == 130) rgb(:, x) = [223, 220, 218]
            IF (rowintensity(x) == 131) rgb(:, x) = [224, 219, 216]
            IF (rowintensity(x) == 132) rgb(:, x) = [225, 219, 215]
            IF (rowintensity(x) == 133) rgb(:, x) = [226, 218, 214]
            IF (rowintensity(x) == 134) rgb(:, x) = [227, 218, 212]
            IF (rowintensity(x) == 135) rgb(:, x) = [228, 217, 211]
            IF (rowintensity(x) == 136) rgb(:, x) = [229, 216, 209]
            IF (rowintensity(x) == 137) rgb(:, x) = [230, 216, 208]
            IF (rowintensity(x) == 138) rgb(:, x) = [231, 215, 206]
            IF (rowintensity(x) == 139) rgb(:, x) = [232, 215, 205]
            IF (rowintensity(x) == 140) rgb(:, x) = [232, 214, 203]
            IF (rowintensity(x) == 141) rgb(:, x) = [233, 213, 202]
            IF (rowintensity(x) == 142) rgb(:, x) = [234, 212, 200]
            IF (rowintensity(x) == 143) rgb(:, x) = [235, 212, 199]
            IF (rowintensity(x) == 144) rgb(:, x) = [236, 211, 197]
            IF (rowintensity(x) == 145) rgb(:, x) = [236, 210, 196]
            IF (rowintensity(x) == 146) rgb(:, x) = [237, 209, 194]
            IF (rowintensity(x) == 147) rgb(:, x) = [238, 209, 193]
            IF (rowintensity(x) == 148) rgb(:, x) = [238, 208, 191]
            IF (rowintensity(x) == 149) rgb(:, x) = [239, 207, 190]
            IF (rowintensity(x) == 150) rgb(:, x) = [240, 206, 188]
            IF (rowintensity(x) == 151) rgb(:, x) = [240, 205, 187]
            IF (rowintensity(x) == 152) rgb(:, x) = [241, 204, 185]
            IF (rowintensity(x) == 153) rgb(:, x) = [241, 203, 184]
            IF (rowintensity(x) == 154) rgb(:, x) = [242, 202, 182]
            IF (rowintensity(x) == 155) rgb(:, x) = [242, 201, 181]
            IF (rowintensity(x) == 156) rgb(:, x) = [243, 200, 179]
            IF (rowintensity(x) == 157) rgb(:, x) = [243, 199, 178]
            IF (rowintensity(x) == 158) rgb(:, x) = [244, 198, 176]
            IF (rowintensity(x) == 159) rgb(:, x) = [244, 197, 174]
            IF (rowintensity(x) == 160) rgb(:, x) = [245, 196, 173]
            IF (rowintensity(x) == 161) rgb(:, x) = [245, 195, 171]
            IF (rowintensity(x) == 162) rgb(:, x) = [245, 194, 170]
            IF (rowintensity(x) == 163) rgb(:, x) = [245, 193, 168]
            IF (rowintensity(x) == 164) rgb(:, x) = [246, 192, 167]
            IF (rowintensity(x) == 165) rgb(:, x) = [246, 191, 165]
            IF (rowintensity(x) == 166) rgb(:, x) = [246, 190, 163]
            IF (rowintensity(x) == 167) rgb(:, x) = [246, 188, 162]
            IF (rowintensity(x) == 168) rgb(:, x) = [247, 187, 160]
            IF (rowintensity(x) == 169) rgb(:, x) = [247, 186, 159]
            IF (rowintensity(x) == 170) rgb(:, x) = [247, 185, 157]
            IF (rowintensity(x) == 171) rgb(:, x) = [247, 184, 156]
            IF (rowintensity(x) == 172) rgb(:, x) = [247, 182, 154]
            IF (rowintensity(x) == 173) rgb(:, x) = [247, 181, 152]
            IF (rowintensity(x) == 174) rgb(:, x) = [247, 180, 151]
            IF (rowintensity(x) == 175) rgb(:, x) = [247, 178, 149]
            IF (rowintensity(x) == 176) rgb(:, x) = [247, 177, 148]
            IF (rowintensity(x) == 177) rgb(:, x) = [247, 176, 146]
            IF (rowintensity(x) == 178) rgb(:, x) = [247, 174, 145]
            IF (rowintensity(x) == 179) rgb(:, x) = [247, 173, 143]
            IF (rowintensity(x) == 180) rgb(:, x) = [247, 172, 141]
            IF (rowintensity(x) == 181) rgb(:, x) = [247, 170, 140]
            IF (rowintensity(x) == 182) rgb(:, x) = [247, 169, 138]
            IF (rowintensity(x) == 183) rgb(:, x) = [247, 167, 137]
            IF (rowintensity(x) == 184) rgb(:, x) = [247, 166, 135]
            IF (rowintensity(x) == 185) rgb(:, x) = [246, 164, 134]
            IF (rowintensity(x) == 186) rgb(:, x) = [246, 163, 132]
            IF (rowintensity(x) == 187) rgb(:, x) = [246, 161, 131]
            IF (rowintensity(x) == 188) rgb(:, x) = [246, 160, 129]
            IF (rowintensity(x) == 189) rgb(:, x) = [245, 158, 127]
            IF (rowintensity(x) == 190) rgb(:, x) = [245, 157, 126]
            IF (rowintensity(x) == 191) rgb(:, x) = [245, 155, 124]
            IF (rowintensity(x) == 192) rgb(:, x) = [244, 154, 123]
            IF (rowintensity(x) == 193) rgb(:, x) = [244, 152, 121]
            IF (rowintensity(x) == 194) rgb(:, x) = [244, 151, 120]
            IF (rowintensity(x) == 195) rgb(:, x) = [243, 149, 118]
            IF (rowintensity(x) == 196) rgb(:, x) = [243, 147, 117]
            IF (rowintensity(x) == 197) rgb(:, x) = [242, 146, 115]
            IF (rowintensity(x) == 198) rgb(:, x) = [242, 144, 114]
            IF (rowintensity(x) == 199) rgb(:, x) = [241, 142, 112]
            IF (rowintensity(x) == 200) rgb(:, x) = [241, 141, 111]
            IF (rowintensity(x) == 201) rgb(:, x) = [240, 139, 109]
            IF (rowintensity(x) == 202) rgb(:, x) = [240, 137, 108]
            IF (rowintensity(x) == 203) rgb(:, x) = [239, 136, 106]
            IF (rowintensity(x) == 204) rgb(:, x) = [238, 134, 105]
            IF (rowintensity(x) == 205) rgb(:, x) = [238, 132, 103]
            IF (rowintensity(x) == 206) rgb(:, x) = [237, 130, 102]
            IF (rowintensity(x) == 207) rgb(:, x) = [236, 129, 100]
            IF (rowintensity(x) == 208) rgb(:, x) = [236, 127, 99]
            IF (rowintensity(x) == 209) rgb(:, x) = [235, 125, 97]
            IF (rowintensity(x) == 210) rgb(:, x) = [234, 123, 96]
            IF (rowintensity(x) == 211) rgb(:, x) = [233, 121, 95]
            IF (rowintensity(x) == 212) rgb(:, x) = [233, 120, 93]
            IF (rowintensity(x) == 213) rgb(:, x) = [232, 118, 92]
            IF (rowintensity(x) == 214) rgb(:, x) = [231, 116, 90]
            IF (rowintensity(x) == 215) rgb(:, x) = [230, 114, 89]
            IF (rowintensity(x) == 216) rgb(:, x) = [229, 112, 88]
            IF (rowintensity(x) == 217) rgb(:, x) = [228, 110, 86]
            IF (rowintensity(x) == 218) rgb(:, x) = [227, 108, 85]
            IF (rowintensity(x) == 219) rgb(:, x) = [227, 106, 83]
            IF (rowintensity(x) == 220) rgb(:, x) = [226, 104, 82]
            IF (rowintensity(x) == 221) rgb(:, x) = [225, 102, 81]
            IF (rowintensity(x) == 222) rgb(:, x) = [224, 100, 79]
            IF (rowintensity(x) == 223) rgb(:, x) = [223, 98, 78]
            IF (rowintensity(x) == 224) rgb(:, x) = [222, 96, 77]
            IF (rowintensity(x) == 225) rgb(:, x) = [221, 94, 75]
            IF (rowintensity(x) == 226) rgb(:, x) = [220, 92, 74]
            IF (rowintensity(x) == 227) rgb(:, x) = [218, 90, 73]
            IF (rowintensity(x) == 228) rgb(:, x) = [217, 88, 71]
            IF (rowintensity(x) == 229) rgb(:, x) = [216, 86, 70]
            IF (rowintensity(x) == 230) rgb(:, x) = [215, 84, 69]
            IF (rowintensity(x) == 231) rgb(:, x) = [214, 82, 67]
            IF (rowintensity(x) == 232) rgb(:, x) = [213, 80, 66]
            IF (rowintensity(x) == 233) rgb(:, x) = [212, 78, 65]
            IF (rowintensity(x) == 234) rgb(:, x) = [210, 75, 64]
            IF (rowintensity(x) == 235) rgb(:, x) = [209, 73, 62]
            IF (rowintensity(x) == 236) rgb(:, x) = [208, 71, 61]
            IF (rowintensity(x) == 237) rgb(:, x) = [207, 69, 60]
            IF (rowintensity(x) == 238) rgb(:, x) = [205, 66, 59]
            IF (rowintensity(x) == 239) rgb(:, x) = [204, 64, 57]
            IF (rowintensity(x) == 240) rgb(:, x) = [203, 62, 56]
            IF (rowintensity(x) == 241) rgb(:, x) = [202, 59, 55]
            IF (rowintensity(x) == 242) rgb(:, x) = [200, 57, 54]
            IF (rowintensity(x) == 243) rgb(:, x) = [199, 54, 53]
            IF (rowintensity(x) == 244) rgb(:, x) = [198, 51, 52]
            IF (rowintensity(x) == 245) rgb(:, x) = [196, 49, 50]
            IF (rowintensity(x) == 246) rgb(:, x) = [195, 46, 49]
            IF (rowintensity(x) == 247) rgb(:, x) = [193, 43, 48]
            IF (rowintensity(x) == 248) rgb(:, x) = [192, 40, 47]
            IF (rowintensity(x) == 249) rgb(:, x) = [190, 37, 46]
            IF (rowintensity(x) == 250) rgb(:, x) = [189, 34, 45]
            IF (rowintensity(x) == 251) rgb(:, x) = [188, 30, 44]
            IF (rowintensity(x) == 252) rgb(:, x) = [186, 26, 43]
            IF (rowintensity(x) == 253) rgb(:, x) = [185, 22, 41]
            IF (rowintensity(x) == 254) rgb(:, x) = [183, 17, 40]
            IF (rowintensity(x) == 255) rgb(:, x) = [181, 11, 39]
         CASE (2)
            IF (rowintensity(x) < 85) THEN
               rgb(1, x) = INT((255*rowintensity(x))/85)
               rgb(2, x) = 0
               rgb(3, x) = 0
            ELSE IF (rowintensity(x) < 170) THEN
               rgb(1, x) = 255
               rgb(2, x) = INT((255*(rowintensity(x) - 85))/85)
               rgb(3, x) = 0
            ELSE
               rgb(1, x) = 255
               rgb(2, x) = 255
               rgb(3, x) = INT((255*(rowintensity(x) - 170))/85)
            END IF
         CASE (3)
            IF (rowintensity(x) == 0) rgb(:, x) = [0, 0, 131]
            IF (rowintensity(x) == 1) rgb(:, x) = [0, 0, 135]
            IF (rowintensity(x) == 2) rgb(:, x) = [0, 0, 139]
            IF (rowintensity(x) == 3) rgb(:, x) = [0, 0, 143]
            IF (rowintensity(x) == 4) rgb(:, x) = [0, 0, 147]
            IF (rowintensity(x) == 5) rgb(:, x) = [0, 0, 151]
            IF (rowintensity(x) == 6) rgb(:, x) = [0, 0, 155]
            IF (rowintensity(x) == 7) rgb(:, x) = [0, 0, 159]
            IF (rowintensity(x) == 8) rgb(:, x) = [0, 0, 163]
            IF (rowintensity(x) == 9) rgb(:, x) = [0, 0, 167]
            IF (rowintensity(x) == 10) rgb(:, x) = [0, 0, 171]
            IF (rowintensity(x) == 11) rgb(:, x) = [0, 0, 175]
            IF (rowintensity(x) == 12) rgb(:, x) = [0, 0, 179]
            IF (rowintensity(x) == 13) rgb(:, x) = [0, 0, 183]
            IF (rowintensity(x) == 14) rgb(:, x) = [0, 0, 187]
            IF (rowintensity(x) == 15) rgb(:, x) = [0, 0, 191]
            IF (rowintensity(x) == 16) rgb(:, x) = [0, 0, 195]
            IF (rowintensity(x) == 17) rgb(:, x) = [0, 0, 199]
            IF (rowintensity(x) == 18) rgb(:, x) = [0, 0, 203]
            IF (rowintensity(x) == 19) rgb(:, x) = [0, 0, 207]
            IF (rowintensity(x) == 20) rgb(:, x) = [0, 0, 211]
            IF (rowintensity(x) == 21) rgb(:, x) = [0, 0, 215]
            IF (rowintensity(x) == 22) rgb(:, x) = [0, 0, 219]
            IF (rowintensity(x) == 23) rgb(:, x) = [0, 0, 223]
            IF (rowintensity(x) == 24) rgb(:, x) = [0, 0, 227]
            IF (rowintensity(x) == 25) rgb(:, x) = [0, 0, 231]
            IF (rowintensity(x) == 26) rgb(:, x) = [0, 0, 235]
            IF (rowintensity(x) == 27) rgb(:, x) = [0, 0, 239]
            IF (rowintensity(x) == 28) rgb(:, x) = [0, 0, 243]
            IF (rowintensity(x) == 29) rgb(:, x) = [0, 0, 247]
            IF (rowintensity(x) == 30) rgb(:, x) = [0, 0, 251]
            IF (rowintensity(x) == 31) rgb(:, x) = [0, 0, 255]
            IF (rowintensity(x) == 32) rgb(:, x) = [0, 4, 255]
            IF (rowintensity(x) == 33) rgb(:, x) = [0, 8, 255]
            IF (rowintensity(x) == 34) rgb(:, x) = [0, 12, 255]
            IF (rowintensity(x) == 35) rgb(:, x) = [0, 16, 255]
            IF (rowintensity(x) == 36) rgb(:, x) = [0, 20, 255]
            IF (rowintensity(x) == 37) rgb(:, x) = [0, 24, 255]
            IF (rowintensity(x) == 38) rgb(:, x) = [0, 28, 255]
            IF (rowintensity(x) == 39) rgb(:, x) = [0, 32, 255]
            IF (rowintensity(x) == 40) rgb(:, x) = [0, 36, 255]
            IF (rowintensity(x) == 41) rgb(:, x) = [0, 40, 255]
            IF (rowintensity(x) == 42) rgb(:, x) = [0, 44, 255]
            IF (rowintensity(x) == 43) rgb(:, x) = [0, 48, 255]
            IF (rowintensity(x) == 44) rgb(:, x) = [0, 52, 255]
            IF (rowintensity(x) == 45) rgb(:, x) = [0, 56, 255]
            IF (rowintensity(x) == 46) rgb(:, x) = [0, 60, 255]
            IF (rowintensity(x) == 47) rgb(:, x) = [0, 64, 255]
            IF (rowintensity(x) == 48) rgb(:, x) = [0, 68, 255]
            IF (rowintensity(x) == 49) rgb(:, x) = [0, 72, 255]
            IF (rowintensity(x) == 50) rgb(:, x) = [0, 76, 255]
            IF (rowintensity(x) == 51) rgb(:, x) = [0, 80, 255]
            IF (rowintensity(x) == 52) rgb(:, x) = [0, 84, 255]
            IF (rowintensity(x) == 53) rgb(:, x) = [0, 88, 255]
            IF (rowintensity(x) == 54) rgb(:, x) = [0, 92, 255]
            IF (rowintensity(x) == 55) rgb(:, x) = [0, 96, 255]
            IF (rowintensity(x) == 56) rgb(:, x) = [0, 100, 255]
            IF (rowintensity(x) == 57) rgb(:, x) = [0, 104, 255]
            IF (rowintensity(x) == 58) rgb(:, x) = [0, 108, 255]
            IF (rowintensity(x) == 59) rgb(:, x) = [0, 112, 255]
            IF (rowintensity(x) == 60) rgb(:, x) = [0, 116, 255]
            IF (rowintensity(x) == 61) rgb(:, x) = [0, 120, 255]
            IF (rowintensity(x) == 62) rgb(:, x) = [0, 124, 255]
            IF (rowintensity(x) == 63) rgb(:, x) = [0, 128, 255]
            IF (rowintensity(x) == 64) rgb(:, x) = [0, 131, 255]
            IF (rowintensity(x) == 65) rgb(:, x) = [0, 135, 255]
            IF (rowintensity(x) == 66) rgb(:, x) = [0, 139, 255]
            IF (rowintensity(x) == 67) rgb(:, x) = [0, 143, 255]
            IF (rowintensity(x) == 68) rgb(:, x) = [0, 147, 255]
            IF (rowintensity(x) == 69) rgb(:, x) = [0, 151, 255]
            IF (rowintensity(x) == 70) rgb(:, x) = [0, 155, 255]
            IF (rowintensity(x) == 71) rgb(:, x) = [0, 159, 255]
            IF (rowintensity(x) == 72) rgb(:, x) = [0, 163, 255]
            IF (rowintensity(x) == 73) rgb(:, x) = [0, 167, 255]
            IF (rowintensity(x) == 74) rgb(:, x) = [0, 171, 255]
            IF (rowintensity(x) == 75) rgb(:, x) = [0, 175, 255]
            IF (rowintensity(x) == 76) rgb(:, x) = [0, 179, 255]
            IF (rowintensity(x) == 77) rgb(:, x) = [0, 183, 255]
            IF (rowintensity(x) == 78) rgb(:, x) = [0, 187, 255]
            IF (rowintensity(x) == 79) rgb(:, x) = [0, 191, 255]
            IF (rowintensity(x) == 80) rgb(:, x) = [0, 195, 255]
            IF (rowintensity(x) == 81) rgb(:, x) = [0, 199, 255]
            IF (rowintensity(x) == 82) rgb(:, x) = [0, 203, 255]
            IF (rowintensity(x) == 83) rgb(:, x) = [0, 207, 255]
            IF (rowintensity(x) == 84) rgb(:, x) = [0, 211, 255]
            IF (rowintensity(x) == 85) rgb(:, x) = [0, 215, 255]
            IF (rowintensity(x) == 86) rgb(:, x) = [0, 219, 255]
            IF (rowintensity(x) == 87) rgb(:, x) = [0, 223, 255]
            IF (rowintensity(x) == 88) rgb(:, x) = [0, 227, 255]
            IF (rowintensity(x) == 89) rgb(:, x) = [0, 231, 255]
            IF (rowintensity(x) == 90) rgb(:, x) = [0, 235, 255]
            IF (rowintensity(x) == 91) rgb(:, x) = [0, 239, 255]
            IF (rowintensity(x) == 92) rgb(:, x) = [0, 243, 255]
            IF (rowintensity(x) == 93) rgb(:, x) = [0, 247, 255]
            IF (rowintensity(x) == 94) rgb(:, x) = [0, 251, 255]
            IF (rowintensity(x) == 95) rgb(:, x) = [0, 255, 255]
            IF (rowintensity(x) == 96) rgb(:, x) = [4, 255, 251]
            IF (rowintensity(x) == 97) rgb(:, x) = [8, 255, 247]
            IF (rowintensity(x) == 98) rgb(:, x) = [12, 255, 243]
            IF (rowintensity(x) == 99) rgb(:, x) = [16, 255, 239]
            IF (rowintensity(x) == 100) rgb(:, x) = [20, 255, 235]
            IF (rowintensity(x) == 101) rgb(:, x) = [24, 255, 231]
            IF (rowintensity(x) == 102) rgb(:, x) = [28, 255, 227]
            IF (rowintensity(x) == 103) rgb(:, x) = [32, 255, 223]
            IF (rowintensity(x) == 104) rgb(:, x) = [36, 255, 219]
            IF (rowintensity(x) == 105) rgb(:, x) = [40, 255, 215]
            IF (rowintensity(x) == 106) rgb(:, x) = [44, 255, 211]
            IF (rowintensity(x) == 107) rgb(:, x) = [48, 255, 207]
            IF (rowintensity(x) == 108) rgb(:, x) = [52, 255, 203]
            IF (rowintensity(x) == 109) rgb(:, x) = [56, 255, 199]
            IF (rowintensity(x) == 110) rgb(:, x) = [60, 255, 195]
            IF (rowintensity(x) == 111) rgb(:, x) = [64, 255, 191]
            IF (rowintensity(x) == 112) rgb(:, x) = [68, 255, 187]
            IF (rowintensity(x) == 113) rgb(:, x) = [72, 255, 183]
            IF (rowintensity(x) == 114) rgb(:, x) = [76, 255, 179]
            IF (rowintensity(x) == 115) rgb(:, x) = [80, 255, 175]
            IF (rowintensity(x) == 116) rgb(:, x) = [84, 255, 171]
            IF (rowintensity(x) == 117) rgb(:, x) = [88, 255, 167]
            IF (rowintensity(x) == 118) rgb(:, x) = [92, 255, 163]
            IF (rowintensity(x) == 119) rgb(:, x) = [96, 255, 159]
            IF (rowintensity(x) == 120) rgb(:, x) = [100, 255, 155]
            IF (rowintensity(x) == 121) rgb(:, x) = [104, 255, 151]
            IF (rowintensity(x) == 122) rgb(:, x) = [108, 255, 147]
            IF (rowintensity(x) == 123) rgb(:, x) = [112, 255, 143]
            IF (rowintensity(x) == 124) rgb(:, x) = [116, 255, 139]
            IF (rowintensity(x) == 125) rgb(:, x) = [120, 255, 135]
            IF (rowintensity(x) == 126) rgb(:, x) = [124, 255, 131]
            IF (rowintensity(x) == 127) rgb(:, x) = [128, 255, 128]
            IF (rowintensity(x) == 128) rgb(:, x) = [131, 255, 124]
            IF (rowintensity(x) == 129) rgb(:, x) = [135, 255, 120]
            IF (rowintensity(x) == 130) rgb(:, x) = [139, 255, 116]
            IF (rowintensity(x) == 131) rgb(:, x) = [143, 255, 112]
            IF (rowintensity(x) == 132) rgb(:, x) = [147, 255, 108]
            IF (rowintensity(x) == 133) rgb(:, x) = [151, 255, 104]
            IF (rowintensity(x) == 134) rgb(:, x) = [155, 255, 100]
            IF (rowintensity(x) == 135) rgb(:, x) = [159, 255, 96]
            IF (rowintensity(x) == 136) rgb(:, x) = [163, 255, 92]
            IF (rowintensity(x) == 137) rgb(:, x) = [167, 255, 88]
            IF (rowintensity(x) == 138) rgb(:, x) = [171, 255, 84]
            IF (rowintensity(x) == 139) rgb(:, x) = [175, 255, 80]
            IF (rowintensity(x) == 140) rgb(:, x) = [179, 255, 76]
            IF (rowintensity(x) == 141) rgb(:, x) = [183, 255, 72]
            IF (rowintensity(x) == 142) rgb(:, x) = [187, 255, 68]
            IF (rowintensity(x) == 143) rgb(:, x) = [191, 255, 64]
            IF (rowintensity(x) == 144) rgb(:, x) = [195, 255, 60]
            IF (rowintensity(x) == 145) rgb(:, x) = [199, 255, 56]
            IF (rowintensity(x) == 146) rgb(:, x) = [203, 255, 52]
            IF (rowintensity(x) == 147) rgb(:, x) = [207, 255, 48]
            IF (rowintensity(x) == 148) rgb(:, x) = [211, 255, 44]
            IF (rowintensity(x) == 149) rgb(:, x) = [215, 255, 40]
            IF (rowintensity(x) == 150) rgb(:, x) = [219, 255, 36]
            IF (rowintensity(x) == 151) rgb(:, x) = [223, 255, 32]
            IF (rowintensity(x) == 152) rgb(:, x) = [227, 255, 28]
            IF (rowintensity(x) == 153) rgb(:, x) = [231, 255, 24]
            IF (rowintensity(x) == 154) rgb(:, x) = [235, 255, 20]
            IF (rowintensity(x) == 155) rgb(:, x) = [239, 255, 16]
            IF (rowintensity(x) == 156) rgb(:, x) = [243, 255, 12]
            IF (rowintensity(x) == 157) rgb(:, x) = [247, 255, 8]
            IF (rowintensity(x) == 158) rgb(:, x) = [251, 255, 4]
            IF (rowintensity(x) == 159) rgb(:, x) = [255, 255, 0]
            IF (rowintensity(x) == 160) rgb(:, x) = [255, 251, 0]
            IF (rowintensity(x) == 161) rgb(:, x) = [255, 247, 0]
            IF (rowintensity(x) == 162) rgb(:, x) = [255, 243, 0]
            IF (rowintensity(x) == 163) rgb(:, x) = [255, 239, 0]
            IF (rowintensity(x) == 164) rgb(:, x) = [255, 235, 0]
            IF (rowintensity(x) == 165) rgb(:, x) = [255, 231, 0]
            IF (rowintensity(x) == 166) rgb(:, x) = [255, 227, 0]
            IF (rowintensity(x) == 167) rgb(:, x) = [255, 223, 0]
            IF (rowintensity(x) == 168) rgb(:, x) = [255, 219, 0]
            IF (rowintensity(x) == 169) rgb(:, x) = [255, 215, 0]
            IF (rowintensity(x) == 170) rgb(:, x) = [255, 211, 0]
            IF (rowintensity(x) == 171) rgb(:, x) = [255, 207, 0]
            IF (rowintensity(x) == 172) rgb(:, x) = [255, 203, 0]
            IF (rowintensity(x) == 173) rgb(:, x) = [255, 199, 0]
            IF (rowintensity(x) == 174) rgb(:, x) = [255, 195, 0]
            IF (rowintensity(x) == 175) rgb(:, x) = [255, 191, 0]
            IF (rowintensity(x) == 176) rgb(:, x) = [255, 187, 0]
            IF (rowintensity(x) == 177) rgb(:, x) = [255, 183, 0]
            IF (rowintensity(x) == 178) rgb(:, x) = [255, 179, 0]
            IF (rowintensity(x) == 179) rgb(:, x) = [255, 175, 0]
            IF (rowintensity(x) == 180) rgb(:, x) = [255, 171, 0]
            IF (rowintensity(x) == 181) rgb(:, x) = [255, 167, 0]
            IF (rowintensity(x) == 182) rgb(:, x) = [255, 163, 0]
            IF (rowintensity(x) == 183) rgb(:, x) = [255, 159, 0]
            IF (rowintensity(x) == 184) rgb(:, x) = [255, 155, 0]
            IF (rowintensity(x) == 185) rgb(:, x) = [255, 151, 0]
            IF (rowintensity(x) == 186) rgb(:, x) = [255, 147, 0]
            IF (rowintensity(x) == 187) rgb(:, x) = [255, 143, 0]
            IF (rowintensity(x) == 188) rgb(:, x) = [255, 139, 0]
            IF (rowintensity(x) == 189) rgb(:, x) = [255, 135, 0]
            IF (rowintensity(x) == 190) rgb(:, x) = [255, 131, 0]
            IF (rowintensity(x) == 191) rgb(:, x) = [255, 128, 0]
            IF (rowintensity(x) == 192) rgb(:, x) = [255, 124, 0]
            IF (rowintensity(x) == 193) rgb(:, x) = [255, 120, 0]
            IF (rowintensity(x) == 194) rgb(:, x) = [255, 116, 0]
            IF (rowintensity(x) == 195) rgb(:, x) = [255, 112, 0]
            IF (rowintensity(x) == 196) rgb(:, x) = [255, 108, 0]
            IF (rowintensity(x) == 197) rgb(:, x) = [255, 104, 0]
            IF (rowintensity(x) == 198) rgb(:, x) = [255, 100, 0]
            IF (rowintensity(x) == 199) rgb(:, x) = [255, 96, 0]
            IF (rowintensity(x) == 200) rgb(:, x) = [255, 92, 0]
            IF (rowintensity(x) == 201) rgb(:, x) = [255, 88, 0]
            IF (rowintensity(x) == 202) rgb(:, x) = [255, 84, 0]
            IF (rowintensity(x) == 203) rgb(:, x) = [255, 80, 0]
            IF (rowintensity(x) == 204) rgb(:, x) = [255, 76, 0]
            IF (rowintensity(x) == 205) rgb(:, x) = [255, 72, 0]
            IF (rowintensity(x) == 206) rgb(:, x) = [255, 68, 0]
            IF (rowintensity(x) == 207) rgb(:, x) = [255, 64, 0]
            IF (rowintensity(x) == 208) rgb(:, x) = [255, 60, 0]
            IF (rowintensity(x) == 209) rgb(:, x) = [255, 56, 0]
            IF (rowintensity(x) == 210) rgb(:, x) = [255, 52, 0]
            IF (rowintensity(x) == 211) rgb(:, x) = [255, 48, 0]
            IF (rowintensity(x) == 212) rgb(:, x) = [255, 44, 0]
            IF (rowintensity(x) == 213) rgb(:, x) = [255, 40, 0]
            IF (rowintensity(x) == 214) rgb(:, x) = [255, 36, 0]
            IF (rowintensity(x) == 215) rgb(:, x) = [255, 32, 0]
            IF (rowintensity(x) == 216) rgb(:, x) = [255, 28, 0]
            IF (rowintensity(x) == 217) rgb(:, x) = [255, 24, 0]
            IF (rowintensity(x) == 218) rgb(:, x) = [255, 20, 0]
            IF (rowintensity(x) == 219) rgb(:, x) = [255, 16, 0]
            IF (rowintensity(x) == 220) rgb(:, x) = [255, 12, 0]
            IF (rowintensity(x) == 221) rgb(:, x) = [255, 8, 0]
            IF (rowintensity(x) == 222) rgb(:, x) = [255, 4, 0]
            IF (rowintensity(x) == 223) rgb(:, x) = [255, 0, 0]
            IF (rowintensity(x) == 224) rgb(:, x) = [251, 0, 0]
            IF (rowintensity(x) == 225) rgb(:, x) = [247, 0, 0]
            IF (rowintensity(x) == 226) rgb(:, x) = [243, 0, 0]
            IF (rowintensity(x) == 227) rgb(:, x) = [239, 0, 0]
            IF (rowintensity(x) == 228) rgb(:, x) = [235, 0, 0]
            IF (rowintensity(x) == 229) rgb(:, x) = [231, 0, 0]
            IF (rowintensity(x) == 230) rgb(:, x) = [227, 0, 0]
            IF (rowintensity(x) == 231) rgb(:, x) = [223, 0, 0]
            IF (rowintensity(x) == 232) rgb(:, x) = [219, 0, 0]
            IF (rowintensity(x) == 233) rgb(:, x) = [215, 0, 0]
            IF (rowintensity(x) == 234) rgb(:, x) = [211, 0, 0]
            IF (rowintensity(x) == 235) rgb(:, x) = [207, 0, 0]
            IF (rowintensity(x) == 236) rgb(:, x) = [203, 0, 0]
            IF (rowintensity(x) == 237) rgb(:, x) = [199, 0, 0]
            IF (rowintensity(x) == 238) rgb(:, x) = [195, 0, 0]
            IF (rowintensity(x) == 239) rgb(:, x) = [191, 0, 0]
            IF (rowintensity(x) == 240) rgb(:, x) = [187, 0, 0]
            IF (rowintensity(x) == 241) rgb(:, x) = [183, 0, 0]
            IF (rowintensity(x) == 242) rgb(:, x) = [179, 0, 0]
            IF (rowintensity(x) == 243) rgb(:, x) = [175, 0, 0]
            IF (rowintensity(x) == 244) rgb(:, x) = [171, 0, 0]
            IF (rowintensity(x) == 245) rgb(:, x) = [167, 0, 0]
            IF (rowintensity(x) == 246) rgb(:, x) = [163, 0, 0]
            IF (rowintensity(x) == 247) rgb(:, x) = [159, 0, 0]
            IF (rowintensity(x) == 248) rgb(:, x) = [155, 0, 0]
            IF (rowintensity(x) == 249) rgb(:, x) = [151, 0, 0]
            IF (rowintensity(x) == 250) rgb(:, x) = [147, 0, 0]
            IF (rowintensity(x) == 251) rgb(:, x) = [143, 0, 0]
            IF (rowintensity(x) == 252) rgb(:, x) = [139, 0, 0]
            IF (rowintensity(x) == 253) rgb(:, x) = [135, 0, 0]
            IF (rowintensity(x) == 254) rgb(:, x) = [131, 0, 0]
            IF (rowintensity(x) == 255) rgb(:, x) = [128, 0, 0]
         CASE (4)
            rgb(:, x) = INT(rowintensity(x))
         CASE (5)
            rgb(:, x) = INT(255 - rowintensity(x))
         CASE (6)
            IF (rowintensity(x) == 0) rgb(:, x) = [68, 1, 84]
            IF (rowintensity(x) == 1) rgb(:, x) = [69, 2, 86]
            IF (rowintensity(x) == 2) rgb(:, x) = [69, 4, 87]
            IF (rowintensity(x) == 3) rgb(:, x) = [69, 5, 89]
            IF (rowintensity(x) == 4) rgb(:, x) = [70, 7, 90]
            IF (rowintensity(x) == 5) rgb(:, x) = [70, 8, 92]
            IF (rowintensity(x) == 6) rgb(:, x) = [70, 10, 93]
            IF (rowintensity(x) == 7) rgb(:, x) = [71, 11, 95]
            IF (rowintensity(x) == 8) rgb(:, x) = [71, 13, 96]
            IF (rowintensity(x) == 9) rgb(:, x) = [71, 14, 98]
            IF (rowintensity(x) == 10) rgb(:, x) = [71, 16, 99]
            IF (rowintensity(x) == 11) rgb(:, x) = [72, 17, 100]
            IF (rowintensity(x) == 12) rgb(:, x) = [72, 19, 102]
            IF (rowintensity(x) == 13) rgb(:, x) = [72, 20, 103]
            IF (rowintensity(x) == 14) rgb(:, x) = [72, 22, 104]
            IF (rowintensity(x) == 15) rgb(:, x) = [72, 23, 106]
            IF (rowintensity(x) == 16) rgb(:, x) = [72, 24, 107]
            IF (rowintensity(x) == 17) rgb(:, x) = [72, 26, 108]
            IF (rowintensity(x) == 18) rgb(:, x) = [72, 27, 109]
            IF (rowintensity(x) == 19) rgb(:, x) = [72, 28, 110]
            IF (rowintensity(x) == 20) rgb(:, x) = [72, 30, 112]
            IF (rowintensity(x) == 21) rgb(:, x) = [73, 31, 113]
            IF (rowintensity(x) == 22) rgb(:, x) = [72, 32, 114]
            IF (rowintensity(x) == 23) rgb(:, x) = [72, 34, 115]
            IF (rowintensity(x) == 24) rgb(:, x) = [72, 35, 116]
            IF (rowintensity(x) == 25) rgb(:, x) = [72, 36, 117]
            IF (rowintensity(x) == 26) rgb(:, x) = [72, 37, 118]
            IF (rowintensity(x) == 27) rgb(:, x) = [72, 39, 119]
            IF (rowintensity(x) == 28) rgb(:, x) = [72, 40, 120]
            IF (rowintensity(x) == 29) rgb(:, x) = [72, 41, 121]
            IF (rowintensity(x) == 30) rgb(:, x) = [72, 42, 122]
            IF (rowintensity(x) == 31) rgb(:, x) = [72, 44, 123]
            IF (rowintensity(x) == 32) rgb(:, x) = [71, 45, 124]
            IF (rowintensity(x) == 33) rgb(:, x) = [71, 46, 125]
            IF (rowintensity(x) == 34) rgb(:, x) = [71, 47, 125]
            IF (rowintensity(x) == 35) rgb(:, x) = [71, 49, 126]
            IF (rowintensity(x) == 36) rgb(:, x) = [70, 50, 127]
            IF (rowintensity(x) == 37) rgb(:, x) = [70, 51, 128]
            IF (rowintensity(x) == 38) rgb(:, x) = [70, 52, 128]
            IF (rowintensity(x) == 39) rgb(:, x) = [70, 54, 129]
            IF (rowintensity(x) == 40) rgb(:, x) = [69, 55, 130]
            IF (rowintensity(x) == 41) rgb(:, x) = [69, 56, 130]
            IF (rowintensity(x) == 42) rgb(:, x) = [69, 57, 131]
            IF (rowintensity(x) == 43) rgb(:, x) = [68, 58, 132]
            IF (rowintensity(x) == 44) rgb(:, x) = [68, 60, 132]
            IF (rowintensity(x) == 45) rgb(:, x) = [67, 61, 133]
            IF (rowintensity(x) == 46) rgb(:, x) = [67, 62, 133]
            IF (rowintensity(x) == 47) rgb(:, x) = [67, 63, 134]
            IF (rowintensity(x) == 48) rgb(:, x) = [66, 64, 134]
            IF (rowintensity(x) == 49) rgb(:, x) = [66, 66, 135]
            IF (rowintensity(x) == 50) rgb(:, x) = [65, 67, 135]
            IF (rowintensity(x) == 51) rgb(:, x) = [65, 68, 136]
            IF (rowintensity(x) == 52) rgb(:, x) = [65, 69, 136]
            IF (rowintensity(x) == 53) rgb(:, x) = [64, 70, 136]
            IF (rowintensity(x) == 54) rgb(:, x) = [64, 71, 137]
            IF (rowintensity(x) == 55) rgb(:, x) = [63, 73, 137]
            IF (rowintensity(x) == 56) rgb(:, x) = [63, 74, 138]
            IF (rowintensity(x) == 57) rgb(:, x) = [62, 75, 138]
            IF (rowintensity(x) == 58) rgb(:, x) = [62, 76, 138]
            IF (rowintensity(x) == 59) rgb(:, x) = [61, 77, 138]
            IF (rowintensity(x) == 60) rgb(:, x) = [61, 78, 139]
            IF (rowintensity(x) == 61) rgb(:, x) = [60, 79, 139]
            IF (rowintensity(x) == 62) rgb(:, x) = [60, 80, 139]
            IF (rowintensity(x) == 63) rgb(:, x) = [59, 81, 139]
            IF (rowintensity(x) == 64) rgb(:, x) = [59, 83, 140]
            IF (rowintensity(x) == 65) rgb(:, x) = [58, 84, 140]
            IF (rowintensity(x) == 66) rgb(:, x) = [58, 85, 140]
            IF (rowintensity(x) == 67) rgb(:, x) = [57, 86, 140]
            IF (rowintensity(x) == 68) rgb(:, x) = [57, 87, 140]
            IF (rowintensity(x) == 69) rgb(:, x) = [56, 88, 141]
            IF (rowintensity(x) == 70) rgb(:, x) = [56, 89, 141]
            IF (rowintensity(x) == 71) rgb(:, x) = [55, 90, 141]
            IF (rowintensity(x) == 72) rgb(:, x) = [55, 91, 141]
            IF (rowintensity(x) == 73) rgb(:, x) = [54, 92, 141]
            IF (rowintensity(x) == 74) rgb(:, x) = [54, 93, 141]
            IF (rowintensity(x) == 75) rgb(:, x) = [53, 94, 141]
            IF (rowintensity(x) == 76) rgb(:, x) = [53, 95, 142]
            IF (rowintensity(x) == 77) rgb(:, x) = [52, 96, 142]
            IF (rowintensity(x) == 78) rgb(:, x) = [52, 97, 142]
            IF (rowintensity(x) == 79) rgb(:, x) = [52, 98, 142]
            IF (rowintensity(x) == 80) rgb(:, x) = [51, 99, 142]
            IF (rowintensity(x) == 81) rgb(:, x) = [51, 100, 142]
            IF (rowintensity(x) == 82) rgb(:, x) = [50, 101, 142]
            IF (rowintensity(x) == 83) rgb(:, x) = [50, 102, 142]
            IF (rowintensity(x) == 84) rgb(:, x) = [49, 103, 142]
            IF (rowintensity(x) == 85) rgb(:, x) = [49, 104, 142]
            IF (rowintensity(x) == 86) rgb(:, x) = [48, 105, 142]
            IF (rowintensity(x) == 87) rgb(:, x) = [48, 106, 142]
            IF (rowintensity(x) == 88) rgb(:, x) = [48, 107, 143]
            IF (rowintensity(x) == 89) rgb(:, x) = [47, 108, 143]
            IF (rowintensity(x) == 90) rgb(:, x) = [47, 109, 143]
            IF (rowintensity(x) == 91) rgb(:, x) = [46, 110, 143]
            IF (rowintensity(x) == 92) rgb(:, x) = [46, 111, 143]
            IF (rowintensity(x) == 93) rgb(:, x) = [45, 112, 143]
            IF (rowintensity(x) == 94) rgb(:, x) = [45, 113, 143]
            IF (rowintensity(x) == 95) rgb(:, x) = [45, 114, 143]
            IF (rowintensity(x) == 96) rgb(:, x) = [44, 115, 143]
            IF (rowintensity(x) == 97) rgb(:, x) = [44, 116, 143]
            IF (rowintensity(x) == 98) rgb(:, x) = [43, 117, 143]
            IF (rowintensity(x) == 99) rgb(:, x) = [43, 118, 143]
            IF (rowintensity(x) == 100) rgb(:, x) = [43, 119, 143]
            IF (rowintensity(x) == 101) rgb(:, x) = [42, 120, 143]
            IF (rowintensity(x) == 102) rgb(:, x) = [42, 121, 143]
            IF (rowintensity(x) == 103) rgb(:, x) = [42, 122, 143]
            IF (rowintensity(x) == 104) rgb(:, x) = [41, 123, 143]
            IF (rowintensity(x) == 105) rgb(:, x) = [41, 123, 143]
            IF (rowintensity(x) == 106) rgb(:, x) = [40, 124, 143]
            IF (rowintensity(x) == 107) rgb(:, x) = [40, 125, 143]
            IF (rowintensity(x) == 108) rgb(:, x) = [40, 126, 143]
            IF (rowintensity(x) == 109) rgb(:, x) = [39, 127, 143]
            IF (rowintensity(x) == 110) rgb(:, x) = [39, 128, 143]
            IF (rowintensity(x) == 111) rgb(:, x) = [39, 129, 143]
            IF (rowintensity(x) == 112) rgb(:, x) = [38, 130, 143]
            IF (rowintensity(x) == 113) rgb(:, x) = [38, 131, 143]
            IF (rowintensity(x) == 114) rgb(:, x) = [37, 132, 143]
            IF (rowintensity(x) == 115) rgb(:, x) = [37, 133, 142]
            IF (rowintensity(x) == 116) rgb(:, x) = [37, 134, 142]
            IF (rowintensity(x) == 117) rgb(:, x) = [36, 135, 142]
            IF (rowintensity(x) == 118) rgb(:, x) = [36, 136, 142]
            IF (rowintensity(x) == 119) rgb(:, x) = [36, 137, 142]
            IF (rowintensity(x) == 120) rgb(:, x) = [35, 138, 142]
            IF (rowintensity(x) == 121) rgb(:, x) = [35, 139, 142]
            IF (rowintensity(x) == 122) rgb(:, x) = [35, 139, 142]
            IF (rowintensity(x) == 123) rgb(:, x) = [34, 140, 142]
            IF (rowintensity(x) == 124) rgb(:, x) = [34, 141, 142]
            IF (rowintensity(x) == 125) rgb(:, x) = [34, 142, 141]
            IF (rowintensity(x) == 126) rgb(:, x) = [33, 143, 141]
            IF (rowintensity(x) == 127) rgb(:, x) = [33, 144, 141]
            IF (rowintensity(x) == 128) rgb(:, x) = [33, 145, 141]
            IF (rowintensity(x) == 129) rgb(:, x) = [32, 146, 141]
            IF (rowintensity(x) == 130) rgb(:, x) = [32, 147, 141]
            IF (rowintensity(x) == 131) rgb(:, x) = [32, 148, 140]
            IF (rowintensity(x) == 132) rgb(:, x) = [32, 149, 140]
            IF (rowintensity(x) == 133) rgb(:, x) = [31, 150, 140]
            IF (rowintensity(x) == 134) rgb(:, x) = [31, 151, 140]
            IF (rowintensity(x) == 135) rgb(:, x) = [31, 152, 139]
            IF (rowintensity(x) == 136) rgb(:, x) = [31, 153, 139]
            IF (rowintensity(x) == 137) rgb(:, x) = [31, 154, 139]
            IF (rowintensity(x) == 138) rgb(:, x) = [31, 155, 139]
            IF (rowintensity(x) == 139) rgb(:, x) = [31, 156, 138]
            IF (rowintensity(x) == 140) rgb(:, x) = [31, 156, 138]
            IF (rowintensity(x) == 141) rgb(:, x) = [31, 157, 138]
            IF (rowintensity(x) == 142) rgb(:, x) = [31, 158, 137]
            IF (rowintensity(x) == 143) rgb(:, x) = [31, 159, 137]
            IF (rowintensity(x) == 144) rgb(:, x) = [31, 160, 137]
            IF (rowintensity(x) == 145) rgb(:, x) = [31, 161, 136]
            IF (rowintensity(x) == 146) rgb(:, x) = [31, 162, 136]
            IF (rowintensity(x) == 147) rgb(:, x) = [32, 163, 135]
            IF (rowintensity(x) == 148) rgb(:, x) = [32, 164, 135]
            IF (rowintensity(x) == 149) rgb(:, x) = [32, 165, 134]
            IF (rowintensity(x) == 150) rgb(:, x) = [33, 166, 134]
            IF (rowintensity(x) == 151) rgb(:, x) = [33, 167, 134]
            IF (rowintensity(x) == 152) rgb(:, x) = [34, 168, 133]
            IF (rowintensity(x) == 153) rgb(:, x) = [34, 169, 133]
            IF (rowintensity(x) == 154) rgb(:, x) = [35, 170, 132]
            IF (rowintensity(x) == 155) rgb(:, x) = [36, 170, 131]
            IF (rowintensity(x) == 156) rgb(:, x) = [37, 171, 131]
            IF (rowintensity(x) == 157) rgb(:, x) = [38, 172, 130]
            IF (rowintensity(x) == 158) rgb(:, x) = [38, 173, 130]
            IF (rowintensity(x) == 159) rgb(:, x) = [39, 174, 129]
            IF (rowintensity(x) == 160) rgb(:, x) = [40, 175, 128]
            IF (rowintensity(x) == 161) rgb(:, x) = [41, 176, 128]
            IF (rowintensity(x) == 162) rgb(:, x) = [43, 177, 127]
            IF (rowintensity(x) == 163) rgb(:, x) = [44, 178, 126]
            IF (rowintensity(x) == 164) rgb(:, x) = [45, 179, 126]
            IF (rowintensity(x) == 165) rgb(:, x) = [46, 180, 125]
            IF (rowintensity(x) == 166) rgb(:, x) = [48, 180, 124]
            IF (rowintensity(x) == 167) rgb(:, x) = [49, 181, 123]
            IF (rowintensity(x) == 168) rgb(:, x) = [50, 182, 123]
            IF (rowintensity(x) == 169) rgb(:, x) = [52, 183, 122]
            IF (rowintensity(x) == 170) rgb(:, x) = [53, 184, 121]
            IF (rowintensity(x) == 171) rgb(:, x) = [55, 185, 120]
            IF (rowintensity(x) == 172) rgb(:, x) = [56, 186, 119]
            IF (rowintensity(x) == 173) rgb(:, x) = [58, 187, 118]
            IF (rowintensity(x) == 174) rgb(:, x) = [60, 187, 118]
            IF (rowintensity(x) == 175) rgb(:, x) = [61, 188, 117]
            IF (rowintensity(x) == 176) rgb(:, x) = [63, 189, 116]
            IF (rowintensity(x) == 177) rgb(:, x) = [65, 190, 115]
            IF (rowintensity(x) == 178) rgb(:, x) = [67, 191, 114]
            IF (rowintensity(x) == 179) rgb(:, x) = [68, 192, 113]
            IF (rowintensity(x) == 180) rgb(:, x) = [70, 193, 112]
            IF (rowintensity(x) == 181) rgb(:, x) = [72, 193, 111]
            IF (rowintensity(x) == 182) rgb(:, x) = [74, 194, 110]
            IF (rowintensity(x) == 183) rgb(:, x) = [76, 195, 109]
            IF (rowintensity(x) == 184) rgb(:, x) = [78, 196, 108]
            IF (rowintensity(x) == 185) rgb(:, x) = [80, 197, 106]
            IF (rowintensity(x) == 186) rgb(:, x) = [82, 197, 105]
            IF (rowintensity(x) == 187) rgb(:, x) = [84, 198, 104]
            IF (rowintensity(x) == 188) rgb(:, x) = [86, 199, 103]
            IF (rowintensity(x) == 189) rgb(:, x) = [88, 200, 102]
            IF (rowintensity(x) == 190) rgb(:, x) = [90, 200, 101]
            IF (rowintensity(x) == 191) rgb(:, x) = [92, 201, 99]
            IF (rowintensity(x) == 192) rgb(:, x) = [95, 202, 98]
            IF (rowintensity(x) == 193) rgb(:, x) = [97, 203, 97]
            IF (rowintensity(x) == 194) rgb(:, x) = [99, 203, 95]
            IF (rowintensity(x) == 195) rgb(:, x) = [101, 204, 94]
            IF (rowintensity(x) == 196) rgb(:, x) = [103, 205, 93]
            IF (rowintensity(x) == 197) rgb(:, x) = [106, 206, 91]
            IF (rowintensity(x) == 198) rgb(:, x) = [108, 206, 90]
            IF (rowintensity(x) == 199) rgb(:, x) = [110, 207, 89]
            IF (rowintensity(x) == 200) rgb(:, x) = [113, 208, 87]
            IF (rowintensity(x) == 201) rgb(:, x) = [115, 208, 86]
            IF (rowintensity(x) == 202) rgb(:, x) = [117, 209, 84]
            IF (rowintensity(x) == 203) rgb(:, x) = [120, 210, 83]
            IF (rowintensity(x) == 204) rgb(:, x) = [122, 210, 81]
            IF (rowintensity(x) == 205) rgb(:, x) = [125, 211, 80]
            IF (rowintensity(x) == 206) rgb(:, x) = [127, 212, 78]
            IF (rowintensity(x) == 207) rgb(:, x) = [130, 212, 77]
            IF (rowintensity(x) == 208) rgb(:, x) = [132, 213, 75]
            IF (rowintensity(x) == 209) rgb(:, x) = [135, 213, 74]
            IF (rowintensity(x) == 210) rgb(:, x) = [137, 214, 72]
            IF (rowintensity(x) == 211) rgb(:, x) = [140, 215, 71]
            IF (rowintensity(x) == 212) rgb(:, x) = [142, 215, 69]
            IF (rowintensity(x) == 213) rgb(:, x) = [145, 216, 67]
            IF (rowintensity(x) == 214) rgb(:, x) = [147, 216, 66]
            IF (rowintensity(x) == 215) rgb(:, x) = [150, 217, 64]
            IF (rowintensity(x) == 216) rgb(:, x) = [153, 217, 62]
            IF (rowintensity(x) == 217) rgb(:, x) = [155, 218, 61]
            IF (rowintensity(x) == 218) rgb(:, x) = [158, 218, 59]
            IF (rowintensity(x) == 219) rgb(:, x) = [160, 219, 57]
            IF (rowintensity(x) == 220) rgb(:, x) = [163, 219, 55]
            IF (rowintensity(x) == 221) rgb(:, x) = [166, 220, 54]
            IF (rowintensity(x) == 222) rgb(:, x) = [168, 220, 52]
            IF (rowintensity(x) == 223) rgb(:, x) = [171, 221, 50]
            IF (rowintensity(x) == 224) rgb(:, x) = [174, 221, 49]
            IF (rowintensity(x) == 225) rgb(:, x) = [176, 222, 47]
            IF (rowintensity(x) == 226) rgb(:, x) = [179, 222, 45]
            IF (rowintensity(x) == 227) rgb(:, x) = [182, 222, 43]
            IF (rowintensity(x) == 228) rgb(:, x) = [184, 223, 42]
            IF (rowintensity(x) == 229) rgb(:, x) = [187, 223, 40]
            IF (rowintensity(x) == 230) rgb(:, x) = [190, 224, 38]
            IF (rowintensity(x) == 231) rgb(:, x) = [192, 224, 37]
            IF (rowintensity(x) == 232) rgb(:, x) = [195, 224, 35]
            IF (rowintensity(x) == 233) rgb(:, x) = [198, 225, 34]
            IF (rowintensity(x) == 234) rgb(:, x) = [201, 225, 32]
            IF (rowintensity(x) == 235) rgb(:, x) = [203, 225, 31]
            IF (rowintensity(x) == 236) rgb(:, x) = [206, 226, 29]
            IF (rowintensity(x) == 237) rgb(:, x) = [209, 226, 28]
            IF (rowintensity(x) == 238) rgb(:, x) = [211, 226, 27]
            IF (rowintensity(x) == 239) rgb(:, x) = [214, 227, 26]
            IF (rowintensity(x) == 240) rgb(:, x) = [216, 227, 26]
            IF (rowintensity(x) == 241) rgb(:, x) = [219, 227, 25]
            IF (rowintensity(x) == 242) rgb(:, x) = [222, 228, 25]
            IF (rowintensity(x) == 243) rgb(:, x) = [224, 228, 24]
            IF (rowintensity(x) == 244) rgb(:, x) = [227, 228, 24]
            IF (rowintensity(x) == 245) rgb(:, x) = [229, 229, 25]
            IF (rowintensity(x) == 246) rgb(:, x) = [232, 229, 25]
            IF (rowintensity(x) == 247) rgb(:, x) = [235, 229, 26]
            IF (rowintensity(x) == 248) rgb(:, x) = [237, 230, 27]
            IF (rowintensity(x) == 249) rgb(:, x) = [240, 230, 28]
            IF (rowintensity(x) == 250) rgb(:, x) = [242, 230, 29]
            IF (rowintensity(x) == 251) rgb(:, x) = [245, 231, 30]
            IF (rowintensity(x) == 252) rgb(:, x) = [247, 231, 32]
            IF (rowintensity(x) == 253) rgb(:, x) = [249, 231, 33]
            IF (rowintensity(x) == 254) rgb(:, x) = [252, 232, 35]
            IF (rowintensity(x) == 255) rgb(:, x) = [254, 232, 37]

         CASE default
            WRITE (*, *) 'Error writing .ppm file'
            STOP
         END SELECT
         IF (i == 11 .AND. isnan(intensity_img2(x, y))) rgb(:, x) = [255, 255, 255]
         END DO
         WRITE (i, '(4000(3(I3,2X),2X))') (rgb(:, x), x=1, width2)
      END DO
      CLOSE (i)

   END DO

END SUBROUTINE write_ppm
