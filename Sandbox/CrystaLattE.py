#!/usr/bin/env python

import math
import psi4
import re
import os
import sys
from psi4.driver.wrapper_autofrag import auto_fragments

# Imports of outsourced code.
sys.path.insert(0, "Read_CIF")
import Read_CIF

## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Import two functions from the PyCIFRW module.
#from CifFile import CifFile, CifBlock
#
## =============================================================================
#
## Shows an error message and the usage.
#def print_error(msg):
#
#    print('')
#    print('    ERROR:  %s' % msg)
#    print('')
#    print_usage()
#
## =============================================================================
#
## =============================================================================
#
## Converts an "_atom_type_label" into an element name.
#def extract_element(label):
#
#    elem2 = ['He','Li','Be','Ne','Na','Mg','Al','Si','Cl','Ar','Ca','Sc','Ti',
#             'Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',
#             'Rb','Sr','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
#             'Sb','Te','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd',
#             'Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','Re','Os','Ir','Pt',
#             'Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa',
#             'Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr']
#
#    if (label[0:2] in elem2):
#        return label[0:2]
#
#    elem1 = ['H','B','C','N','O','F','P','S','K','V','Y','I','W','U']
#
#    if (label[0] in elem1):
#        return label[0]
#
#    print('WARNING: could not convert "%s" into element name!' % label)
#    return label
#
## =============================================================================
#
## =============================================================================
#
## Basic XYZ format.
#def write_xyz(atoms, box, f):
#
#    # Write the number of atoms.
#    N = len(atoms)
#    f.write('%d\n' % N)
#
#    # Write a comment line with the box size.
#    if (len(box) == 6):
#        # box = (ax,bx,by,cx,cy,cz)
#        ax = box[0]
#        bx = box[1]
#        by = box[2]
#        cx = box[3]
#        cy = box[4]
#        cz = box[5]
#        f.write('Crystal created from CIF file. Box vectors:')
#        f.write(' a= %10.5f %10.5f %10.5f' % (ax, 0.0, 0.0))
#        f.write(' b= %10.5f %10.5f %10.5f' % (bx, by, 0.0))
#        f.write(' c= %10.5f %10.5f %10.5f\n' % (cx, cy, cz))
#    else:
#        # box = (ax,by,cz)
#        f.write('Crystal created from CIF file. Box size:')
#        f.write(' %10.5f %10.5f %10.5f\n' % box)
#
#    # Write atom data (units are Angstroms).
#    # The argument "atoms" has format ('Si', x, y, z) for example
#    for i in range(N):
#        f.write('%-10s %10.6f %10.6f %10.6f\n' % atoms[i])
#
## =============================================================================
#
## =============================================================================
#
## Read CIF file, and extract the necessary info in the form of a dictionary.
## E.g., the value of "_cell_volume" can be found with data['_cell_volume'].
#def read_cif(fNameIn):
#
#    data = {}
#
#
#    # Open the CIF file and read all the lines into a list of strings.
#    try:
#        f = open(fNameIn, 'r')
#        lines = []
#        for line in f:
#            stripped = line.strip()
#            if (len(stripped) > 0):  lines.append(stripped)
#    except:
#        print('Failed to open CIF file "{0}".format(fNameIn)') # Carlos Borca (2017-09-19-1431)
#        sys.exit()
#
#    # Use the CifFile parser to extract the data.  Although there might be
#    # multiple data blocks, we'll only use the first one.
#    cif_file = CifFile(fNameIn)
#
#    for db in cif_file:
#        data_block = db
#        break
#
#
#    try:
#
#        # Extract some parameters, and convert them to floats.
#        data['_cell_length_a']    = float(data_block['_cell_length_a'])
#        data['_cell_length_b']    = float(data_block['_cell_length_b'])
#        data['_cell_length_c']    = float(data_block['_cell_length_c'])
#        data['_cell_angle_alpha'] = float(data_block['_cell_angle_alpha'])
#        data['_cell_angle_beta']  = float(data_block['_cell_angle_beta'])
#        data['_cell_angle_gamma'] = float(data_block['_cell_angle_gamma'])
#        data['_cell_volume']      = float(data_block['_cell_volume'])
#
#        # Get the symbolic operations that define the space group.  In a CIF file
#        # that's the part that looks like:
#        #
#        # loop_
#        # _symmetry_equiv_pos_as_xyz
#        #   'x,y,z'
#        #   'y,x,2/3-z'
#        #   '-y,x-y,2/3+z'
#        #   '-x,-x+y,1/3-z'
#        #   '-x+y,-x,1/3+z'
#        #   'x-y,-y,-z'
#        #
#        # In some cases it's called "_space_group_symop_operation_xyz" apparently?!?!
#        data['_symmetry_equiv_pos_as_xyz'] = []
#
#        try:
#            xyz = data_block["_symmetry_equiv_pos_as_xyz"]
#
#        except KeyError:
#            try:
#                xyz = data_block["_space_group_symop_operation_xyz"]
#            except KeyError:
#                print('Missing item in CIF file: need either \'_symmetry_equiv_pos_as_xyz\' or \'_space_group_symop_operation_xyz\'.') # Carlos Borca (2017-09-19-1424)
#                sys.exit()
#
#
#        # Copy the x,y,z symmetry group operations.  Remove the quotes if there
#        # are any.
#        for op_xyz in xyz:
#
#            if (op_xyz[0] == '\''):
#                data['_symmetry_equiv_pos_as_xyz'].append(op_xyz[1:-1])
#            else:
#                data['_symmetry_equiv_pos_as_xyz'].append(op_xyz)
#
#        # Add x,y,z of the atoms to "data", but make sure to convert
#        # e.g. "0.1549(8)" to "0.1549".
#        data['_atom_site_label'] = data_block['_atom_site_label']
#
#        data['_atom_site_fract_x'] = []
#        for str_x in data_block['_atom_site_fract_x']:
#            data['_atom_site_fract_x'].append( float(str_x.split('(')[0]) )
#
#        data['_atom_site_fract_y'] = []
#        for str_y in data_block['_atom_site_fract_y']:
#            data['_atom_site_fract_y'].append( float(str_y.split('(')[0]) )
#
#        data['_atom_site_fract_z'] = []
#        for str_z in data_block['_atom_site_fract_z']:
#            data['_atom_site_fract_z'].append( float(str_z.split('(')[0]) )
#
#
#    except KeyError as e:
#        print('Error!  Missing item in file.') # Carlos Borca (2017-09-19-1427)
#        print ('e') # Carlos Borca (2017-09-19-1427)
#
#        sys.exit()
#
#
#    #print "ALL DATA:"
#    #print data
#    #print
#
#    # Return the extracted data.
#    return data
#
## =============================================================================
#
## =============================================================================
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Read arguments given.
#
#def read_cif_driver(args):
#
#        # Default settings.
#        fNameIn = ''
#        fNameOut = ''
#        Nx = 1
#        Ny = 1
#        Nz = 1
#        make_rect_box = False
#
#
#        # Read the arguments.  We expect at least 4.
#        if (len(args) <= 4):
#            print_usage()
#
#        i = 1
#        while (i < len(args)):
#
#
#            # Check if the name of the input file was given.
#            if (args[i] == '-i'):
#
#                # Make sure a file name is given.
#                if (i+1 == len(args)):
#                    print_error('no input file name given')
#
#                fNameIn = args[i+1]
#                i = i + 2
#
#
#            # Check if the name of the output file was given.
#            elif (args[i] == '-o'):
#
#                # Make sure a file name is given.
#                if (i+1 == len(args)):
#                    print_error('no output file name given')
#
#                # Check we have a valid file extension.
#                fNameOut = args[i+1]
#                unknown = True
#
#
#                for ext in ['.xyz', '.lammpstrj', '.gro', '.cif']:
#                    if (fNameOut.endswith(ext)):
#                        unknown = False
#
#                if (unknown):
#                    print_error('unknown file extension of output file')
#
#                i = i + 2
#
#
#            # Check if the box size was given.
#            elif (args[i] == '-b'):
#
#                # Make sure 3 integers are given.
#                if (i+3 >= len(args)):
#                    print_error('need 3 integers to indicate box size')
#
#                Nx = int(args[i+1])
#                Ny = int(args[i+2])
#                Nz = int(args[i+3])
#
#                if (Nx == 0  or  Ny == 0  or  Nz == 0):
#                    print_error('box size integers need to be larger than zero')
#
#                i = i + 4
#
#
#            # Check if the final configuration should be in a rectangular shape, or in
#            # the same shape as the unit cell.
#            elif (args[i] == '-r'):
#
#                make_rect_box = True
#                i = i + 1
#
#
#            # Anything else is wrong.
#            else:
#                print_error('invalid argument "%s"' % args[i])
#
#
#        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#        # Read input file.
#
#
#        # Make sure an input file was given.
#        if (fNameIn == ''):
#            print_error('no input file given.  Use:  -i filename')
#
#
#        # Open the CIF file and read the data.
#        data = read_cif(fNameIn)
#
#
#        # Extract lengths and angles from the CIF file.
#        La = float(data['_cell_length_a'])
#        Lb = float(data['_cell_length_b'])
#        Lc = float(data['_cell_length_c'])
#        alpha = radians(float(data['_cell_angle_alpha']))
#        beta = radians(float(data['_cell_angle_beta']))
#        gamma = radians(float(data['_cell_angle_gamma']))
#        volume = float(data['_cell_volume'])
#
#
#        # Extract the symmetry operations.  This will be a list of strings such as:
#        #    ['x,y,z', 'y,x,2/3-z', '-y,x-y,2/3+z', '-x,-x+y,1/3-z', ... ]
#        ops = data['_symmetry_equiv_pos_as_xyz']
#
#        # For proper evaluation, we need to convert "2/3" into "2./3", etc. to prevent
#        # integer division which would turn e.g. 2/3 into 0.
#        for i in range(len(ops)):
#            ops[i] = ops[i].replace("0/", "0./") # also for e.g. 10/9
#            ops[i] = ops[i].replace("1/", "1./")
#            ops[i] = ops[i].replace("2/", "2./")
#            ops[i] = ops[i].replace("3/", "3./")
#            ops[i] = ops[i].replace("4/", "4./")
#            ops[i] = ops[i].replace("5/", "5./")
#            ops[i] = ops[i].replace("6/", "6./")
#            ops[i] = ops[i].replace("7/", "7./")
#            ops[i] = ops[i].replace("8/", "8./")
#            ops[i] = ops[i].replace("9/", "9./")
#        #    ops[i] = ops[i].replace("/", "./")
#
#        # Get the atom labels and coordinates.
#        labels = data['_atom_site_label']
#        fX = [ float(s) for s in data['_atom_site_fract_x'] ]
#        fY = [ float(s) for s in data['_atom_site_fract_y'] ]
#        fZ = [ float(s) for s in data['_atom_site_fract_z'] ]
#
#        # Create a list of 4-tuples, where each tuple is an atom:
#        #   [ ('Si', 0.4697, 0.0, 0.0),  ('O', 0.4135, 0.2669, 0.1191),  ... ]
#        atoms = [ (labels[i], fX[i], fY[i], fZ[i]) for i in range(len(labels)) ]
#
#        # Make sure that all atoms lie within the unit cell.  Also convert names such
#        # as 'Oa1' into 'O'.
#        for i in range(len(atoms)):
#            (name,xn,yn,zn) = atoms[i]
#            xn = (xn + 10.0) % 1.0
#            yn = (yn + 10.0) % 1.0
#            zn = (zn + 10.0) % 1.0
#            name = extract_element(name)
#            atoms[i] = (name,xn,yn,zn)
#
#
#        # Update the user.
#        print('Loaded a CIF file with %d atom coordinates and %d symmetry operations.' % (len(atoms), len(ops))) # Carlos Borca (2017-09-19-1429)
#        print('') # Carlos Borca (2017-12-06-1554)
#
#        # Just for reference, here is a typical example of a CIF file:
#        """
#        _cell_length_a 4.916
#        _cell_length_b 4.916
#        _cell_length_c 5.4054
#        _cell_angle_alpha 90
#        _cell_angle_beta 90
#        _cell_angle_gamma 120
#        _cell_volume 113.131
#        _exptl_crystal_density_diffrn      2.646
#        _symmetry_space_group_name_H-M 'P 32 2 1'
#        loop_
#        _space_group_symop_operation_xyz
#          'x,y,z'
#          'y,x,2/3-z'
#          '-y,x-y,2/3+z'
#          '-x,-x+y,1/3-z'
#          '-x+y,-x,1/3+z'
#          'x-y,-y,-z'
#        loop_
#        _atom_site_label
#        _atom_site_fract_x
#        _atom_site_fract_y
#        _atom_site_fract_z
#        Si   0.46970   0.00000   0.00000
#        O   0.41350   0.26690   0.11910
#        """
#
#
#        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#        # Use symmetry operations to create the unit cell.
#
#
#        # The CIF file consists of a few atom positions plus several "symmetry
#        # operations" that indicate the other atom positions within the unit cell.  So
#        # using these operations, create copies of the atoms until no new copies can be
#        # made.
#
#
#        # Two atoms are on top of each other if they are less than "eps" away.
#        eps = 0.01  # in Angstrom
#
#
#        # For each atom, apply each symmetry operation to create a new atom.
#        imax = len(atoms)
#        i=0
#        while (i < imax):
#
#            label,x,y,z = atoms[i]
#
#            for op in ops:
#
#                # Python is awesome: calling e.g. eval('x,y,1./2+z') will convert the
#                # string into a 3-tuple using the current values for x,y,z!
#                xn,yn,zn = eval(op)
#
#                # Make sure that the new atom lies within the unit cell.
#                xn = (xn + 10.0) % 1.0
#                yn = (yn + 10.0) % 1.0
#                zn = (zn + 10.0) % 1.0
#
#                # Check if the new position is actually new, or the same as a previous
#                # atom.
#                new_atom = True
#                for at in atoms:
#                    if (abs(at[1]-xn) < eps  and  abs(at[2]-yn) < eps  and  abs(at[3]-zn) < eps):
#                        new_atom = False
#
#                        # Check that this is the same atom type.
#                        if (at[0] != label):
#                            print_error('invalid CIF file: atom of type %s overlaps with atom of type %s' % (at[0],label))
#
#                # If the atom is new, add it to the list!
#                if (new_atom):
#                    atoms.append( (label,xn,yn,zn) )  # add a 4-tuple
#
#
#            # Update the loop iterator.
#            i = i + 1
#            imax = len(atoms)
#
#
#        # Sort the atoms according to type alphabetically.
#        atoms = sorted(atoms, key=lambda at: at[0])
#        atoms.reverse()
#
#
#        # Done with creating the unit cell.  Update the user.
#        print('Created a unit cell consisting of %d atoms.' % len(atoms))
#        print('') # Carlos Borca (2017-12-06-1554)
#
#        print('Fractional coordinates:')
#        for atom in atoms:
#            print('%10s  %.3f  %.3f  %.3f' % atom)
#
#        print ('') # Carlos Borca (2017-12-06-1554)
#
#        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#        # Create a larger box made of several unit cells: the super cell.
#
#
#        atomlist = []
#
#        for atom in atoms:
#
#            # Get label and fractional coordinates.
#            label,xf,yf,zf = atom
#
#            for i in range(Nx):
#                    x = i+xf
#
#                    for j in range(Ny):
#                        y = j+yf
#
#                        for k in range(Nz):
#                            z = k+zf
#                            atomlist.append( (label,x,y,z) ) # add 4-tuple
#
#        atoms = atomlist
#
#        # If the user wants us to create a copy of the current CIF file, with
#        # additional atoms, then do that.  Note that the atoms here have *fractional*
#        # coordinates!
#        if (fNameOut.endswith('.cif')):
#
#            write_cif(fNameIn, atoms, fNameOut)
#
#            print('Done writing extended CIF file (%d atoms in total).' % len(atoms))
#            exit(0)
#
#
#        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#        # Convert the fractional coordinates into real coordinates.
#
#
#        # The primitive vectors a,b,c are such that
#        #
#        #   cos(alpha) = b.c / |b||c|
#        #   cos(beta)  = a.c / |a||c|
#        #   cos(gamma) = a.b / |a||b|
#        #
#        # with the convention
#        #
#        #   a = La*xhat
#        #   b = bx*xhat + by*yhat
#        #   c = cx*xhat + cy*yhat + cz*zhat
#        #
#        cosa = cos(alpha)
#        sina = sin(alpha)
#        cosb = cos(beta)
#        sinb = sin(beta)
#        cosg = cos(gamma)
#        sing = sin(gamma)
#
#        cosa2 = cosa * cosa
#        cosb2 = cosb * cosb
#        sing2 = sing * sing
#
#        ax = La
#
#        bx = Lb * cosg
#        by = Lb * sing
#
#        cx = Lc * cosb
#        cy = Lc * (cosa - cosg*cosb) / sing
#        cz = Lc * sqrt( 1 - (cosa2 + cosb2 - 2*cosg*cosb*cosa) / sing2 )
#
#
#        # Use the volume to check if we did the vectors right.
#        V = ax*by*cz
#        if ( abs(V - volume) > 0.1):
#            print_error('volume does not match that calculated from primitive vectors')
#
#
#        # Check if we have a rectangular box.
#        if (bx < eps  and  cx < eps  and cy < eps):
#            make_rect_box = True
#
#
#        # Update the user.
#        print('The primitive unit cell vectors are:')
#        print('   a = [%6.4f, %6.4f, %6.4f]' % (ax,0,0))
#        print('   b = [%6.4f, %6.4f, %6.4f]' % (bx,by,0))
#        print('   c = [%6.4f, %6.4f, %6.4f]' % (cx,cy,cz))
#        print('') # Carlos Borca (2017-12-06-1558)
#        print('This gives a volume of %f A^3.' % V) # Carlos Borca (2017-12-06-1558)
#        print('CIF file indicates it is %f A^3.' % volume) # Carlos Borca (2017-12-06-1558)
#        print('') # Carlos Borca (2017-12-06-1559)
#
#
#        # Determine the box size.
#        Lx = Nx * La
#        Ly = Ny * Lb
#        Lz = Nz * Lc
#
#
#        for i in range(len(atoms)):
#
#            # Get label and fractional coordinates.
#            label,xf,yf,zf = atoms[i]
#
#            xa = xf * ax  # contribution of a-vector to the x-coordinate of this atom
#            #ya = 0       # a-vector has no y-component, so does not affect y of atom
#            #za = 0       # a-vector has no z-component, so does not affect z of atom
#
#            xb = yf * bx  # contribution of b-vector to the x-coordinate of this atom
#            yb = yf * by  # contribution of b-vector to the y-coordinate of this atom
#            #zb = 0       # b-vector has no z-component, so does not affect z of atom
#
#            xc = zf * cx  # contribution of c-vector to the x-coordinate of this atom
#            yc = zf * cy  # contribution of c-vector to the y-coordinate of this atom
#            zc = zf * cz  # contribution of c-vector to the z-coordinate of this atom
#
#            # Add all contributions.
#            xn = xa + xb + xc
#            yn = yb + yc
#            zn = zc
#
#            if (make_rect_box):
#                xn = (xn + Lx) % Lx
#                yn = (yn + Ly) % Ly
#                zn = (zn + Lz) % Lz
#
#            atoms[i] = (label, xn, yn, zn)
#
#
#        # Determine the box-vector.
#        if (make_rect_box):
#            box = (Lx, Ly, Lz)
#        else:
#            box = (Lx, Ly, Lz, Nx*cx, Ny*cy, Nz*cz)
#
#
#        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#        # Create the output file.
#
#        try:
#            fOut = open(fNameOut, 'w')
#
#            if (fNameOut.endswith('.xyz')):
#                write_xyz(atoms, box, fOut)
#
#            elif (fNameOut.endswith('.lammpstrj')):
#                write_lammpstrj(atoms, box, fOut)
#
#            elif (fNameOut.endswith('.gro')):
#                write_gro(atoms, box, fOut)
#
#        except:
#            print_error('Failed to write to output file')
#
#
#        print('Created output file %s (%d atoms in total).' % (fNameOut, len(atoms)))
#        fOut.close()
#
#        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ==================================================================
def psifrg2xyz(p4frag, numfrags, numfatoms, frg_separator):
    """Takes the output of the automatic fragmentation procedure, the total number of fragments,
    the number of atoms in one fragment, and the string that separates the fragments.
    Returns .xyz files for each fragment."""

    # Read the output of the automatic fragmentation.
    with open(p4frag) as cellfrags:
        frags = cellfrags.readlines()

    # Generate .xyz files for each fragment.
    fcounter = 1         # Counter of processed fragments.
    lcounter = 0         # Counter of lines.
    HeaderLine = True    # Flag for first line of a fragment .xyz file.

    for line in frags:
        lcounter += 1

        if line.startswith(frg_separator):
            fcounter += 1
            HeaderLine = True

        else:
            ffnidx = "f" + str(fcounter).zfill(len(str(numfrags))) \
                     + ".xyz"
            with open(ffnidx, "a") as frgxyz: # WARNING: File exists?

                if HeaderLine == True:
                    frgxyz.write(str(numfatoms) + "\n" + "Fragment " # TODO: What if there are more than two molecules in the central unit cell.
                    + str(fcounter) + "\n")
                    HeaderLine = False

                if HeaderLine == False:
                    frgxyz.write(line)
    
    print("%i fragment .xyz files were created." % fcounter)

    return
# ==================================================================

# ==================================================================
def frags_filter(n_atm_frg): # WARNING: The number of atoms in a fragment probably should be passed as a list.
    """It takes the possible numbers of atoms in a fragment and discards all fragments files that contain a smaller number of atoms."""

    # TODO: What if there are more than two types of molecules in the
    #       unit cell? n_atom_frg should be passed as a list with
    #       possible numbers of atoms in a fragment.

    directory = os.getcwd()
    files = os.listdir(directory)
    incmolcounter = 0

    for file in files:

        if re.match('^f[0-9]+.xyz$', file): # Match filenames f???.xyz

            with open(file, 'r') as f:
                lines = f.readlines()

                if len(lines) != n_atm_frg + 2: # Match files with unexpected atoms.
                    print("Expected " + str(n_atm_frg)\
                          + " atoms. Found only " + str(len(lines) - 2)\
                          + ". Removing: " + file)

                    os.remove(os.path.join(directory,file))

                    incmolcounter += 1

    print("\nRemoved %s fragments.\n" % incmolcounter)

    return
# ==================================================================

# ==================================================================
def scellcntr(xyzpattern):
    """Compute the center of coordinates of .xyz files with a matching pattern in the filename.
    Returns a list with fragment names, x, y, and z coordinates in lists of floats and the coordinates of the center."""

    frg = [] # Fragments file names.
    x = [] # Values of x coordinates.
    y = [] # Values of y coordinates.
    z = [] # Values of z coordinates.

    directory = os.getcwd()
    files = os.listdir(directory)

    for file in files:

        if re.match(xyzpattern, file):
            
            with open(file, 'r') as f:
                lines = f.readlines()

                i = 0
                for i in range(len(lines)):

                    if i < 2:
                        continue

                    else:
                        info = lines[i].split()
                        frg.append(file)
                        x.append(float(info[1]))
                        y.append(float(info[2]))
                        z.append(float(info[3]))

    cntr_x = (max(x) - min(x))/2.0
    cntr_y = (max(y) - min(y))/2.0
    cntr_z = (max(z) - min(z))/2.0

    print("Center of the supercell located at:")
    print("x= %10.5f"   % (cntr_x))
    print("y= %10.5f"   % (cntr_y))
    print("z= %10.5f\n" % (cntr_z))

    return frg, x, y, z, cntr_x, cntr_y, cntr_z 
# ==================================================================

# ==================================================================
def posvecmag(xyzpattern):
    """Using the scellcntr() function it computes the magnitude of the position vector of the atoms read by the scellcntr() function.
    Returns the filename and the magnitude of the position vector."""

    frg, x, y, z, cntr_x, cntr_y, cntr_z = scellcntr(xyzpattern)

    r = []

    for j in range(len(x)):
        r.append(math.sqrt((x[j] - cntr_x)**2.0 + (y[j] - cntr_y)**2.0 + (z[j] - cntr_z)**2.0))

    print("Magnitudes of the position vectors calculated.\n")

    return frg, r
# ==================================================================

# ==================================================================
def proximity():
    """Using position vector magnitude, computed in posvecmag(). This function returns a list of fragments filenames
    with fragments containing atoms closest to the center of the supercell first."""

    frg, r = posvecmag("^f[0-9]+.xyz$")

    frgprox = []

    for pair in sorted(enumerate(r), key=lambda item:item[1]):
        if frg[pair[0]] in frgprox:
            continue
        else:
            frgprox.append(frg[pair[0]])

    print("Molecules organized based on the closest atom to the")
    print("center of the supercell.\n")

    return frgprox
# ==================================================================

# ==================================================================
def frgs2mols():
    """Takes fragents and a list of proximity to the center of the supercell, generated by the proximity() function.
    Writes new files with molecules indexed by increasing proximity to the center of the supercell."""

    directory = os.getcwd()
    proxlist = proximity()
    molecule = 0

    for file in proxlist:
    #for file in proxlist[:10]: # NOTE: For tests
        molecule += 1
        molfname = "m" + str(molecule).zfill(len(proxlist[-1]) - 5) + ".xyz"

        with open(file, 'r') as ffile, open(molfname, 'w') as mfile:

            for l in ffile.readlines():

                if l.startswith("Fragment"):
                    newl = "Molecule " + str(molecule) + " (" + l[:-1] + ")\n"
                    mfile.write(newl)

                else:
                    mfile.write(l)

        os.remove(os.path.join(directory, file))

        print("Fragment " + str(file) + " has been rewritten to: " + str(molfname))

    print("")
# ==================================================================

# ==================================================================
def mols2mons():
    """Takes molecules and filters and a list of proximity to the center of the supercell, generated by the proximity() function.
    Writes new files with molecules indexed by increasing proximity to the center of the supercell."""

    directory = os.getcwd()
    molfiles = os.listdir(directory)
    monidx = 0

    for file in molfiles:
        
        if re.match('^m[0-9]+.xyz$', file):
            monidx += 1
            monfname = "1-" + str(file)[1:]

            with open(file, 'r') as moleculef, open(monfname, 'w') as monomerf:

                for l in moleculef.readlines():

                    if l.startswith("Molecule"):
                        newl = "Monomer " + str(monidx) + "\n"
                        monomerf.write(newl)

                    else:
                        monomerf.write(l)

            os.remove(os.path.join(directory, file))

            print("Molecule " + str(file) + " has been rewritten to: " + str(monfname))

    print("")
# ==================================================================

# ==================================================================
def rmin(nmer1, nmer2):
    """Takes two monomer files and returns separation between closest atoms belonging to different monomers."""

    with open(nmer1, 'r') as nmf1, open(nmer2, 'r') as nmf2:
        nm1l = nmf1.readlines()
        nm2l = nmf2.readlines()

        nmf1info = []
        nmf2info = []

        dist = []

        i = 0
        for i in range(len(nm1l)):

            if i < 2:
                continue

            else:
                nmf1atom = nm1l[i].split()
                nmf1info.append(nmf1atom)

        for j in range(len(nm2l)):

            if j < 2:
                continue

            else:
                nmf2atom = nm2l[j].split()
                nmf2info.append(nmf2atom)

        for atom1 in range(len(nmf1info)):

            for atom2 in range(len(nmf2info)):
                x12 = float(nmf1info[atom1][1]) - float(nmf2info[atom2][1])
                y12 = float(nmf1info[atom1][2]) - float(nmf2info[atom2][2])
                z12 = float(nmf1info[atom1][3]) - float(nmf2info[atom2][3])
                r12 = math.sqrt(x12**2.0 + y12**2.0 + z12**2.0)
                dist.append(r12)

    minsep = min(dist)

    return minsep
# ==================================================================

# ==================================================================
def nmerbuilder(nmertype, rcut):
    """Takes a float indicating a cutoff distance in Angstroms.
    Returns dimer files that pass through the filter."""

    if nmertype == "dimers":
        print("Merging monomers with monomers to obtain dimers.\n")
        nmerpatt = "^1-[0-9]+.xyz$"
        ucnmlbls = "Dimer"
        numnmlbl = "2"

    elif nmertype == "trimers":
        print("Merging dimers with monomers to obtain trimers.\n")
        nmerpatt = "^2-[0-9]+.xyz$"
        ucnmlbls = "Trimer"
        numnmlbl = "3"

    if nmertype == "tetramers":
        print("Merging trimers with monomers to obtain tetramers.\n")
        nmerpatt = "^3-[0-9]+.xyz$"
        ucnmlbls = "Tetramer"
        numnmlbl = "4"

    directory = os.getcwd()
    nmerfs = os.listdir(directory)
    nmidx = 0
    dscrd = 0

    for nmer in nmerfs:

        if re.match(nmerpatt, nmer):

            for monomer in nmerfs:

                if re.match('^1-[0-9]+.xyz$', monomer):

                    if nmer < monomer: # WARNING: This may not work for trimers, tetramers, ...
                        nmidx +=1
                        newnm = numnmlbl + "-" + str(nmer)[2:-4] + "+" + str(monomer)[2:]
                        
                        # TODO: Implement monomer-in-the-cell filter for dimers only.

                        mindist = rmin(nmer, monomer) # Separation cutoff filter.

                        if rcut <= mindist:
                            print("%s %i: Discarded. Interaction separation %3.2f greater than cutoff %3.2f" % (ucnmlbls, nmidx, mindist, rcut))
                            dscrd += 1

                        else:

                            with open(nmer, 'r') as mf, open(monomer, 'r') as nf, open(newnm, 'w') as newf:
                                
                                l1idx = 0
                                for l1 in mf.readlines():
                                    l1idx += 1

                                    if l1idx == 1:
                                        newfl1 = str(int(l1) + int(nf.readline())) + "\n"
                                        newf.write(newfl1)

                                    elif l1.startswith("Monomer"):
                                        newl1 = ucnmlbls + " " + str(nmidx) + " (Monomers " + str(nmer)[2:-4] + "+" + str(monomer)[2:-4] + ")\n"
                                        newf.write(newl1)

                                    elif l1.startswith("Dimer"):
                                        newl1 = ucnmlbls + " " + str(nmidx) + " (Monomers " + str(nmer)[2:-4] + "+" + str(monomer)[2:-4] + ")\n"
                                        newf.write(newl1)

                                    elif l1.startswith("Trimer"):
                                        newl1 = ucnmlbls + " " + str(nmidx) + " (Monomers " + str(nmer)[2:-4] + "+" + str(monomer)[2:-4] + ")\n"
                                        newf.write(newl1)
                                    
                                    else:
                                        newf.write(l1)

                                l2idx = 0
                                for l2 in nf.readlines():
                                    l2idx += 1

                                    if l2idx < 2:
                                        continue

                                    else:
                                        newf.write(l2)

                            print(ucnmlbls + " " + str(nmidx) + ": " + str(newnm) + " created from " + str(nmer) + " and " + str(monomer))

    print("")
    print("Discarded %i far-separated n-mers.\n" % dscrd)
# ==================================================================

# ==================================================================
# Main program.
def main():
    "Takes a CIF file and computes the crystal lattice energy using a manybody expansion approach."
    
    print("\n=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n")
    print("                             CrystaLattE                             \n")
    print(" The tool for the automated calculation of crystal lattice energies. \n")
    print("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n")

    # ------------------------------------------------------------------
    # Read a CIF file and generate the unit cell.
    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    # Run the unit cell through auto_fragments() to extract the number
    # of molecules in the central unit cell.
    # UnitCFrags = auto_fragments(molecule = UnitCell)
    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    # Read a CIF file and generates a supercell.
    
    ReadCifIn = 'Benzene-138K.cif'  # CIF input file. # WARNING: Hardcoded for now!
    ReadCifOut = 'Benzene-138K.xyz' # XYZ output file. # WARNING: Hardcoded for now!
    ReadCifA = '5'                  # X replicas. # WARNING: Hardcoded for now!
    ReadCifB = '5'                  # Y replicas. # WARNING: Hardcoded for now!
    ReadCifC = '5'                  # Z replicas. # WARNING: Hardcoded for now!
    
    args = ['', '-i', ReadCifIn, '-o', ReadCifOut, '-b', ReadCifA, ReadCifB, ReadCifC]
    
    print("The following arguments will be passed to the CIF reader script:\n")
    print("./Read_CIF.py" + " ".join(str(arg) for arg in args) + "\n")
    
    Read_CIF.main(args)
    # ------------------------------------------------------------------
    
    # ------------------------------------------------------------------
    # Take the supercell .xyz file
    # And generate .xyz files with all possible monomers.
    print ("")
    
    # Read the lines in the .xyz file recently generated by Read_CIF.py
    with open(ReadCifOut) as fxyz:
        sxyz = fxyz.readlines()
    
    # Generate a SuperCell object.
    #SuperCell = psi4.geometry('\n'.join(sxyz[2:]))
    #SuperCell.update_geometry()
    #print (SuperCell.print_out())
    
    # Generate fragments from SuperCell.
    #CellFrags = auto_fragments(molecule = SuperCell)
    #print (CellFrags.natom())
    #print (CellFrags.nfragments())
    
    # Read the output of the automatic fragmentation.
    
    print("Producing .xyz files for each fragment of the supercell.\n")

    p4frag = "bzfrag.p4" # WARNING: Name of the fragmented super cell file is temporarily hardcoded.
    numfrags = 664       # WARNING: Number of fragments temporarily hardcoded!
    numfatoms = 12       # WARNING: Number of atoms per fragment hardcoded! What if there are two types of molecules?
    frg_separator = "--" # Fragment separator string.

    psifrg2xyz(p4frag, numfrags, numfatoms, frg_separator)

    # Discard fragments that are not a complete molecule.
    
    print("Detecting fragments with incomplete molecules.\n")

    frags_filter(numfatoms)

    # Organize fragments according to their distance to the center of
    # the supercell, produce molecule files.
    
    print("Organizing molecules according to their separation to the")
    print("center of the supercell.\n")

    frgs2mols()
    
    # Create monomer files from fragment files.

    print("Creating monomers from molecule files.\n")

    mols2mons()
    
    # ------------------------------------------------------------------
    
    # ------------------------------------------------------------------
    # Loop through all monomers and generate dimers with all other
    # monomers.
    #
    # Filter dimers that are too distant apart.
    
    rdimer = 10.0

    print("Creating dimers from monomer files.\n")

    nmerbuilder("dimers", rdimer)

    # Filter out and keep count of all non-unique dimers, using the
    # nuclear repulsion energy criteria.
    # ------------------------------------------------------------------
    
    # ------------------------------------------------------------------
    # Loop through all dimers and generate trimers with all other
    # monomers.
    #
    # Filter trimers that are too distant apart.
    #
    # Filter out and keep count of all non-unique trimers, using 
    # ArbAlign.
    # ------------------------------------------------------------------
    
    # .
    # .
    # .
    
    # ------------------------------------------------------------------
    # Run plesantly parallel PSI4 computations on all the final list of 
    # monomers, dimers, trimers, etc.
    #
    # Multiply the resulting energy of each one by the degeneracy factor.
    #
    # Sum results to get a lattice energy.
    # ------------------------------------------------------------------

    print("Execution terminated. Thanks for using CrystaLattE.\n")
# ==================================================================

if __name__ == "__main__":
    main()
