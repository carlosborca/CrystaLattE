#!/usr/bin/env python

#
# @BEGIN LICENSE
#
# CrystaLattE: The tool for the automated calculation of crystal lattice 
# energies.
#
# Copyright (c) 2017-2019 
# Carlos H. Borca
# Brandon W. Bakr
# Lori A. Burns
# C. David Sherrill
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of CrystaLattE.
#
# CrystaLattE is free software; you can redistribute it and/or modify
# it under the tesms of the GNU Lesser General Public License as 
# published by the Free Software Foundation, version 3.
#
# CrystaLattE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public 
# License along with CrystaLattE; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
# 02110-1301 USA.
#
# @END LICENSE
#

# Import standard Python modules.
import itertools
import math
import multiprocessing
import numpy as np
import os
import sys
import shutil
import time

# Import parts of Psi4.
import psi4
#from psi4.driver.qcdb.bfs import BFS

# Import QCEelemental
import qcelemental as qcel

# Import parts of PyCIFRW
from CifFile import CifFile

# ======================================================================
def input_parser(in_f_name):
    """Reads a CrystaLattE input file, sets the arguments for the main()
    function, and calls main() to start the execution of the program.
    
    Arguments:
    <str> in_f_name
        Input filename.
    """

    keywords = {}
    keywords["nmers_up_to"] = 2 
    keywords["cif_a"] = 5
    keywords["cif_b"] = 5
    keywords["cif_c"] = 5
    keywords["bfs_thresh"] = 1.2
    keywords["uniq_filter"] = "ChSEV"
    keywords["r_cut_com"] = 10.0 
    keywords["r_cut_monomer"] = 12.0 
    keywords["r_cut_dimer"] = 10.0 
    keywords["r_cut_trimer"] = 8.0
    keywords["r_cut_tetramer"] = 6.0 
    keywords["r_cut_pentamer"] = 4.0 
    keywords["cle_run_type"] = ["test"]
    keywords["psi4_method"] = "HF/STO-3G"
    keywords["psi4_bsse"] = "cp"
    keywords["psi4_memory"] = "500 MB"
    keywords["verbose"] = 1

    with open(in_f_name, "r") as input_file:
        
        non_blank_lines = []

        for line_inp_f in input_file:

            if line_inp_f == "\n":
                continue

            elif line_inp_f.startswith("#"):
                continue
            
            else:
                non_blank_lines.append(line_inp_f)

        for non_blank_line in non_blank_lines:

            split_line = non_blank_line.split("=")
            keyword_name = split_line[0].lower().strip()

            try:
                keyword_value = split_line[1].split('\n')[0]
            
            except IndexError:
                print("\nERROR: Invalid input file. Check that keywords and values are separated by an equal sign.\n")
                sys.exit()
           
            if keyword_name == "cle_run_type":
                keyword_value = keyword_value.lower()
                keyword_value = keyword_value.replace(" ", "").split("+")
                   
            elif keyword_name in ["uniq_filter", "psi4_bsse", "psi4_memory", "psi4_method", "cif_input", "cif_output"]:
                keyword_value = keyword_value.strip()

            elif keyword_name in ["nmers_up_to", "verbose"]:
                
                try:
                    keyword_value = int(keyword_value)

                except ValueError:
                    print("\nERROR: Invalid input file. Check that given {} is an integer.\n".format(keyword_name))
                    sys.exit()
            
            elif keyword_name in ["cif_a", "cif_b", "cif_c"]:
                
                try:
                    keyword_value = int(keyword_value)
                    
                    if (keyword_value%2) == 0:
                        keyword_value += 1

                except ValueError:
                    print("\nERROR: Invalid input file. Check that given {} is an integer.\n".format(keyword_name))
                    sys.exit()

            else:
                try:
                    keyword_value = float(keyword_value)

                except ValueError:
                    print("\nERROR: Invalid input file. Check that given {} is an float.\n".format(keyword_name))
                    sys.exit()

            keywords[keyword_name] = keyword_value

    # Check CIF is provided.
    if "cif_input" not in keywords.keys():
        print("\nERROR: No CIF given as input.\n")
        sys.exit()

    # Check CIF extension.
    else:
        if not keywords["cif_input"].endswith(".cif"):
            print("\nERROR: Invalid CIF file name. Check that given file ends with a .cif extension.\n")
            sys.exit()

    # Attempt to create a filename for the supercell XYZ file.
    if "cif_output" not in keywords.keys():
    
        # Check proper CIF input filename.
        if keywords["cif_input"].endswith(".cif"):
            # Create supercell filename changing the extension of the input.
            keywords["cif_output"] = keywords["cif_input"][:-4] + ".xyz"

        else:
            print("\nERROR: Invalid CIF file name. Check that given file ends with a .cif extension.\n")
            sys.exit()

    # Print program header.
    if keywords["verbose"] >= 1:

        print_header()

    if keywords["verbose"] >= 2:

        print("CrystaLattE execution setup:\n")
        
        # Get the keys of the keywords dictionary, and put them on a list.
        kw_keys = list(keywords.keys())

        # Sort the list in decreasing order.
        kw_keys.sort()
        
        for kw_key in kw_keys:
            print("  {:15} = {}".format(kw_key, str(keywords[kw_key])))

    if "makefp" in keywords["cle_run_type"]:

        if len(keywords["cle_run_type"]) > 1:

            if "timings" not in keywords["cle_run_type"]:
                print("\nERROR: makefp mode cannot be run at the same time with any other mode, except timings.\n")
                sys.exit()
    
    elif "quiet" in keywords["cle_run_type"]:
        
        if (len(keywords["cle_run_type"]) < 2):
            print("\nERROR: quiet mode be must run together with psi4api or test mode.\n")
            sys.exit()
        
        elif "psithon" in keywords["cle_run_type"]:
            print("\nERROR: quiet and psithon modes cannot be run at the same time.\n")
            sys.exit()

        elif "timings" in keywords["cle_run_type"]:
            print("\nERROR: if running quiet and timings modes together, psi4api mode must be included too.\n")
            sys.exit()

    elif "psi4api" in keywords["cle_run_type"]:

        if "psithon" in keywords["cle_run_type"]:
            print("\nERROR: psi4api and psithon modes cannot be run at the same time.\n")
            sys.exit()

        else:
            print("HELP!! I AM TRAPPED! elif psi4api in keywords[cle_run_type]:")

    if "timings" in keywords["cle_run_type"]:

        if (len(keywords["cle_run_type"]) < 2):
            print("\nERROR: timings mode be must run together with another mode.\n")
            sys.exit()

        else:
        
            func_str = "main("

            for kw_key in kw_keys:
                
                if type(keywords[kw_key]) == str:
                    func_str = func_str + "{}='{}', ".format(kw_key, str(keywords[kw_key]))
                
                else:
                    func_str = func_str + "{}={}, ".format(kw_key, str(keywords[kw_key])) 
            
            func_str = func_str[:-2] + ")"
            
            import cProfile as profile
            profile.run(func_str)

    else:
        main(
            keywords["cif_input"], 
            keywords["cif_output"],
            keywords["cif_a"], 
            keywords["cif_b"],
            keywords["cif_c"], 
            keywords["bfs_thresh"],
            keywords["uniq_filter"], 
            keywords["nmers_up_to"],
            keywords["r_cut_com"], 
            keywords["r_cut_monomer"], 
            keywords["r_cut_dimer"], 
            keywords["r_cut_trimer"], 
            keywords["r_cut_tetramer"], 
            keywords["r_cut_pentamer"], 
            keywords["cle_run_type"],
            keywords["psi4_method"],
            keywords["psi4_bsse"],
            keywords["psi4_memory"],
            keywords["verbose"])

    return keywords
    
# ======================================================================






# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# ======================================================================
def extract_element(label):
    """Converts an "_atom_type_label" into an element name.
    """

    elem2 = ['He','Li','Be','Ne','Na','Mg','Al','Si','Cl','Ar','Ca','Sc','Ti',
             'Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',
             'Rb','Sr','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
             'Sb','Te','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd',
             'Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','Re','Os','Ir','Pt',
             'Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa',
             'Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr']

    if (label[0:2] in elem2):
        return label[0:2]

    elem1 = ['H','B','C','N','O','F','P','S','K','V','Y','I','W','U']

    if (label[0] in elem1):
        return label[0]

    print('WARNING: could not convert "%s" into element name!' % label)
    return label
# ======================================================================


# ======================================================================
def write_xyz(atoms, box, f):
    """Writes a basic XYZ format file.
    """

    # Write the number of atoms.
    N = len(atoms)
    f.write('%d\n' % N)

    # Write a comment line with the box size.
    if (len(box) == 6):
        # box = (ax,bx,by,cx,cy,cz)
        ax = box[0]
        bx = box[1]
        by = box[2]
        cx = box[3]
        cy = box[4]
        cz = box[5]
        f.write('Crystal created from CIF file. Box vectors:')
        f.write(' a= %10.5f %10.5f %10.5f' % (ax, 0.0, 0.0))
        f.write(' b= %10.5f %10.5f %10.5f' % (bx, by, 0.0))
        f.write(' c= %10.5f %10.5f %10.5f\n' % (cx, cy, cz))
    else:
        # box = (ax,by,cz)
        f.write('Crystal created from CIF file. Box size:') 
        f.write(' %10.5f %10.5f %10.5f\n' % box)

    # Write atom data (units are Angstroms).
    # The argument "atoms" has format ('Si', x, y, z) for example
    for i in range(N):
        f.write('%-10s %10.6f %10.6f %10.6f\n' % atoms[i])
# =============================================================================


# =============================================================================
def read_cif(fNameIn):
    """Read CIF file, and extract the necessary info in the form of a 
    dictionary. E.g., the value of "_cell_volume" can be found with 
    data['_cell_volume'].
    """
    data = {}

    # Open the CIF file and read all the lines into a list of strings.
    try:
        f = open(fNameIn, 'r')
        lines = []

        for line in f:
            stripped = line.strip()

            if (len(stripped) > 0):  lines.append(stripped)
    except:
        print("\nERROR: Failed to open CIF file '{0}'".format(fNameIn))
        sys.exit()

    # Use the CifFile parser to extract the data. Although there might
    # be multiple data blocks, we'll only use the first one.
    
    #TODO: Known bug: Sometimes not all the blocks in the CIF files are
    #      read, and the following extraction fails.
    cif_file = CifFile(fNameIn)

    for db in cif_file:
        data_block = db
        break

    try:
        # Extract some parameters, and convert them to floats.
        data['_cell_length_a']    = float(data_block['_cell_length_a'].replace("(", "").replace(")", ""))
        data['_cell_length_b']    = float(data_block['_cell_length_b'].replace("(", "").replace(")", ""))
        data['_cell_length_c']    = float(data_block['_cell_length_c'].replace("(", "").replace(")", ""))
        data['_cell_angle_alpha'] = float(data_block['_cell_angle_alpha'].replace("(", "").replace(")", ""))
        data['_cell_angle_beta']  = float(data_block['_cell_angle_beta'].replace("(", "").replace(")", ""))
        data['_cell_angle_gamma'] = float(data_block['_cell_angle_gamma'].replace("(", "").replace(")", ""))
        data['_cell_volume']      = float(data_block['_cell_volume'].replace("(", "").replace(")", ""))

        # Get the symbolic operations that define the space group. In
        # a CIF file that's the part that looks like:
        
        # loop_
        # _symmetry_equiv_pos_as_xyz
        #   'x,y,z'
        #   'y,x,2/3-z'
        #   '-y,x-y,2/3+z'
        #   '-x,-x+y,1/3-z'
        #   '-x+y,-x,1/3+z'
        #   'x-y,-y,-z'
        
        # In some cases it is called:
        # "_space_group_symop_operation_xyz".
        data['_symmetry_equiv_pos_as_xyz'] = []

        try:
            xyz = data_block["_symmetry_equiv_pos_as_xyz"]

        except KeyError:
        
            try:
                xyz = data_block["_space_group_symop_operation_xyz"]
            
            except KeyError:
                print('\nERROR: Missing item in CIF file: need either \'_symmetry_equiv_pos_as_xyz\' or \'_space_group_symop_operation_xyz\'.')
                sys.exit()

        # Copy the x,y,z symmetry group operations. Remove the quotes 
        # if there are any.
        for op_xyz in xyz:

            if (op_xyz[0] == '\''):
                data['_symmetry_equiv_pos_as_xyz'].append(op_xyz[1:-1])
            
            else:
                data['_symmetry_equiv_pos_as_xyz'].append(op_xyz)


        # Add x,y,z of the atoms to "data", but make sure to convert
        # e.g. "0.1549(8)" to "0.1549".
        data['_atom_site_label'] = data_block['_atom_site_label']

        data['_atom_site_fract_x'] = []
        
        for str_x in data_block['_atom_site_fract_x']:
            data['_atom_site_fract_x'].append( float(str_x.split('(')[0]) )

        data['_atom_site_fract_y'] = []
        
        for str_y in data_block['_atom_site_fract_y']:
            data['_atom_site_fract_y'].append( float(str_y.split('(')[0]) )

        data['_atom_site_fract_z'] = []
        
        for str_z in data_block['_atom_site_fract_z']:
            data['_atom_site_fract_z'].append( float(str_z.split('(')[0]) )

    
    except KeyError as e:
        print('\nERROR: Missing item in file.')
        print(e)
        sys.exit()

    # Return the extracted data.
    return data
# =============================================================================


# =============================================================================
def cif_main(args):

    # Default settings.
    fNameIn = ''
    fNameOut = ''
    Nx = 1
    Ny = 1
    Nz = 1
    make_rect_box = False
    
    
    # Read the arguments. We expect at least 4.
    if (len(args) <= 4):
        print_usage()
    
    i = 1
    while (i < len(args)):
    
        # Check if the name of the input file was given.
        if (args[i] == '-i'):
            
            # Make sure a file name is given.
            if (i+1 == len(args)):
                print('\nERROR: no input file name given')
            
            fNameIn = args[i+1]
            i = i + 2
    
        # Check if the name of the output file was given.
        elif (args[i] == '-o'):
    
            # Make sure a file name is given.
            if (i+1 == len(args)):
                print('\nERROR: no output file name given')
    
            # Check we have a valid file extension.
            fNameOut = args[i+1]
            unknown = True
    
            for ext in ['.xyz', '.lammpstrj', '.gro', '.cif']:
                if (fNameOut.endswith(ext)):
                    unknown = False
    
            if (unknown):
                print('\nERROR: unknown file extension of output file')
    
            i = i + 2
    
        # Check if the box size was given.
        elif (args[i] == '-b'):
    
            # Make sure 3 integers are given.
            if (i+3 >= len(args)):
                print('\nERROR: need 3 integers to indicate box size')
    
            Nx = int(args[i+1])
            Ny = int(args[i+2])
            Nz = int(args[i+3])
    
            if (Nx == 0  or  Ny == 0  or  Nz == 0):
                print('\nERROR: box size integers need to be larger than zero')
    
            i = i + 4
    
        # Check if the final configuration should be in a rectangular
        # shape, or in the same shape as the unit cell.
        elif (args[i] == '-r'):
    
            make_rect_box = True
            i = i + 1
    
    
        # Anything else is wrong.
        else:
            print('\nERROR: invalid argument "%s"' % args[i])
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Read input file.
    
    # Make sure an input file was given.
    if (fNameIn == ''):
        print('\nERROR: no input file given.  Use:  -i filename')
    
    # Open the CIF file and read the data.
    data = read_cif(fNameIn)
    
    # Extract lengths and angles from the CIF file.
    La = float(data['_cell_length_a'])
    Lb = float(data['_cell_length_b'])
    Lc = float(data['_cell_length_c'])
    alpha = math.radians(float(data['_cell_angle_alpha']))
    beta = math.radians(float(data['_cell_angle_beta']))
    gamma = math.radians(float(data['_cell_angle_gamma']))
    volume = float(data['_cell_volume'])
    
    # Extract the symmetry operations.  This will be a list of strings
    # such as:
    #    ['x,y,z', 'y,x,2/3-z', '-y,x-y,2/3+z', '-x,-x+y,1/3-z', ... ]
    ops = data['_symmetry_equiv_pos_as_xyz']
    
    # For proper evaluation, we need to convert "2/3" into "2./3", etc.
    # to prevent integer division which would turn e.g. 2/3 into 0.
    for i in range(len(ops)):
        ops[i] = ops[i].replace("0/", "0./") # also for e.g. 10/9
        ops[i] = ops[i].replace("1/", "1./")
        ops[i] = ops[i].replace("2/", "2./")
        ops[i] = ops[i].replace("3/", "3./")
        ops[i] = ops[i].replace("4/", "4./")
        ops[i] = ops[i].replace("5/", "5./")
        ops[i] = ops[i].replace("6/", "6./")
        ops[i] = ops[i].replace("7/", "7./")
        ops[i] = ops[i].replace("8/", "8./")
        ops[i] = ops[i].replace("9/", "9./")
    #    ops[i] = ops[i].replace("/", "./")
    
    # Get the atom labels and coordinates.
    labels = data['_atom_site_label']
    fX = [ float(s) for s in data['_atom_site_fract_x'] ]
    fY = [ float(s) for s in data['_atom_site_fract_y'] ]
    fZ = [ float(s) for s in data['_atom_site_fract_z'] ]
    
    # Create a list of 4-tuples, where each tuple is an atom:
    #   [ ('Si', 0.4697, 0.0, 0.0),  ('O', 0.4135, 0.2669, 0.1191),  ... ]
    atoms = [ (labels[i], fX[i], fY[i], fZ[i]) for i in range(len(labels)) ]
    
    # Make sure that all atoms lie within the unit cell.  Also convert 
    # names such as 'Oa1' into 'O'.
    for i in range(len(atoms)):
        (name,xn,yn,zn) = atoms[i]
        xn = (xn + 10.0) % 1.0
        yn = (yn + 10.0) % 1.0
        zn = (zn + 10.0) % 1.0
        name = extract_element(name)
        atoms[i] = (name,xn,yn,zn)
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Use symmetry operations to create the unit cell.
    
    # The CIF file consists of a few atom positions plus several "symmetry
    # operations" that indicate the other atom positions within the unit cell.  So
    # using these operations, create copies of the atoms until no new copies can be
    # made.
    
    # Two atoms are on top of each other if they are less than "eps" away.
    eps = 0.01  # in Angstrom
    
    # For each atom, apply each symmetry operation to create a new atom.
    imax = len(atoms)
    i=0
    while (i < imax):
    
        label,x,y,z = atoms[i]
    
        for op in ops:
    
            # Python is awesome: calling e.g. eval('x,y,1./2+z') will convert the
            # string into a 3-tuple using the current values for x,y,z!
            xn,yn,zn = eval(op)
    
            # Make sure that the new atom lies within the unit cell.
            xn = (xn + 10.0) % 1.0
            yn = (yn + 10.0) % 1.0
            zn = (zn + 10.0) % 1.0
    
            # Check if the new position is actually new, or the same as a previous
            # atom.
            new_atom = True
            for at in atoms:
                if (abs(at[1]-xn) < eps  and  abs(at[2]-yn) < eps  and  abs(at[3]-zn) < eps):
                    new_atom = False
    
                    # Check that this is the same atom type.
                    if (at[0] != label):
                        print('\nERROR: invalid CIF file: atom of type %s overlaps with atom of type %s' % (at[0],label))
    
            # If the atom is new, add it to the list!
            if (new_atom):
                atoms.append( (label,xn,yn,zn) )  # add a 4-tuple
    
        # Update the loop iterator.
        i = i + 1
        imax = len(atoms)
    
    # Sort the atoms according to type alphabetically.
    atoms = sorted(atoms, key=lambda at: at[0])
    atoms.reverse()
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Create a larger box made of several unit cells: the super cell.
    
    atomlist = []
    
    for atom in atoms:
    
        # Get label and fractional coordinates.
        label,xf,yf,zf = atom
    
        for i in range(Nx):
                x = i+xf
    
                for j in range(Ny):
                    y = j+yf
    
                    for k in range(Nz):
                        z = k+zf
                        atomlist.append( (label,x,y,z) ) # add 4-tuple
    
    atoms = atomlist
    
    # If the user wants us to create a copy of the current CIF file, with
    # additional atoms, then do that.  Note that the atoms here have *fractional*
    # coordinates!
    if (fNameOut.endswith('.cif')):
        
        write_cif(fNameIn, atoms, fNameOut)
    
        print('\nERROR: Done writing extended CIF file (%d atoms in total).' % len(atoms))
        exit(0)
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Convert the fractional coordinates into real coordinates.
    # The primitive vectors a,b,c are such that 
    #
    #   cos(alpha) = b.c / |b||c|
    #   cos(beta)  = a.c / |a||c|
    #   cos(gamma) = a.b / |a||b|
    #
    # with the convention
    #
    #   a = La*xhat
    #   b = bx*xhat + by*yhat
    #   c = cx*xhat + cy*yhat + cz*zhat
    #
    cosa = math.cos(alpha)
    #sina = math.sin(alpha)
    cosb = math.cos(beta)
    #sinb = math.sin(beta)
    cosg = math.cos(gamma)
    sing = math.sin(gamma)
    
    cosa2 = cosa * cosa
    cosb2 = cosb * cosb
    sing2 = sing * sing
    
    ax = La
    
    bx = Lb * cosg
    by = Lb * sing
    
    cx = Lc * cosb
    cy = Lc * (cosa - cosg*cosb) / sing
    cz = Lc * math.sqrt( 1 - (cosa2 + cosb2 - 2*cosg*cosb*cosa) / sing2 )
    
    # Use the volume to check if we did the vectors right.
    V = ax*by*cz
    
    if ( abs(V - volume) > 0.1):
        print('WARNING: Volume of the unit cell declared in CIF ({:.2f} A^3) is different than the calculated from primitive vectors ({:.2f} A^3).\n'.format(volume, V))
    
    # Check if we have a rectangular box.
    if (bx < eps  and  cx < eps  and cy < eps):
        make_rect_box = True
    
    # Determine the box size.
    Lx = Nx * La
    Ly = Ny * Lb
    Lz = Nz * Lc
    
    for i in range(len(atoms)):
    
        # Get label and fractional coordinates.
        label,xf,yf,zf = atoms[i]
    
        xa = xf * ax  # contribution of a-vector to the x-coordinate of this atom
        #ya = 0       # a-vector has no y-component, so does not affect y of atom
        #za = 0       # a-vector has no z-component, so does not affect z of atom
        
        xb = yf * bx  # contribution of b-vector to the x-coordinate of this atom
        yb = yf * by  # contribution of b-vector to the y-coordinate of this atom
        #zb = 0       # b-vector has no z-component, so does not affect z of atom
    
        xc = zf * cx  # contribution of c-vector to the x-coordinate of this atom
        yc = zf * cy  # contribution of c-vector to the y-coordinate of this atom
        zc = zf * cz  # contribution of c-vector to the z-coordinate of this atom
    
        # Add all contributions.
        xn = xa + xb + xc
        yn = yb + yc
        zn = zc
    
        if (make_rect_box):
            xn = (xn + Lx) % Lx
            yn = (yn + Ly) % Ly
            zn = (zn + Lz) % Lz
    
        atoms[i] = (label, xn, yn, zn)
    
    # Determine the box-vector.
    if (make_rect_box):
        box = (Lx, Ly, Lz)
    
    else:
        box = (Lx, Ly, Lz, Nx*cx, Ny*cy, Nz*cz)
    
    try:
        fOut = open(fNameOut, 'w')
    
        if (fNameOut.endswith('.xyz')):
            write_xyz(atoms, box, fOut)
    
        elif (fNameOut.endswith('.lammpstrj')):
            write_lammpstrj(atoms, box, fOut)
    
        elif (fNameOut.endswith('.gro')):
            write_gro(atoms, box, fOut)
    
    except:
        print('\nERROR: Failed to write to output file')
    
    fOut.close()
    

# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||






# ======================================================================
def cif_driver(cif_input, cif_output, cif_a, cif_b, cif_c, verbose=1):
    """Takes the name of a CIF input file and the name of a .xyz output
    file, as well as the number of replicas of the rectangular cell in
    each direction (A, B, and C). It then calls 
    Read_CIF() and passes that information as arguments to generate an
    .xyz file of the supercell.

    Arguments:
    <str> cif_input
        CIF input filename.
    <str> cif_output
        XYZ output filename.
    <int> cif_a
        Number of replicas of the cartesian unit cell in `a` direction.
    <int> cif_b
        Number of replicas of the cartesian unit cell in `b` direction.
    <int> cif_c
        Number of replicas of the cartesian unit cell in `c` direction.
    <int> verbose
        Adjusts the level of detail of the printouts.
    """

    cif_arguments = ["", "-i", cif_input, "-o", cif_output, "-b", cif_a, cif_b, cif_c, "-r"]

    if verbose >= 2:
        print("\nGenerating the supercell .xyz file.")
        print("\nThe following arguments will be passed to the CIF reader script:")
        print("./Read_CIF.py" + " ".join(str(cif_argument) for cif_argument in cif_arguments) + "\n")
    
    return cif_arguments
# ======================================================================


# ======================================================================
def center_supercell(cif_output, verbose=0):
    """Takes the supercell file produced by Read_CIF and computes the
    center of the supercell coordinates to translate the supercell to
    the origin.

    Arguments:
    <str> cif_output
        Name of the file with the cartesian coordinates of the
        supercell.
    <int> verbose
        Adjusts the level of detail of the printouts.

    Returns:
    <numpy.ndarray> scell_geom_max_coords
        Array of 3 numbers with the maximum value of the coordinates x,
        y, and z of the centered supercell.
    <numpy.ndarray> scell_geom
        3-columns array with the x, y, and z coordinates of the 
        centered supercell.
    <numpy.ndarray> scell_elem
        1-column array with element symbols of the centered supercell.
    """
    
    # Find out if there are duplicate coordinates in the XYZ of the
    # supercell. Generate the indexes of unique atoms to exclude
    # superimposed atoms.
    scell_dupl = np.loadtxt(cif_output, skiprows=2, usecols=(0, 1, 2, 3), dtype="str")
    dirty_cell = [tuple(row) for row in scell_dupl]
    clean_cell, unique_idx = np.unique(dirty_cell, return_index=True, axis=0)

    # Compare the size of the clean (all unique coordinates) supercell
    # with that of the dirty supercell. If different, need to clean up
    # the supercell to avoid superimposed atoms.
    if len(scell_dupl) != len(clean_cell):
        print("WARNING: The supercell contains {} duplicate coordinates. Skipping duplicates.".format(len(scell_dupl) - len(clean_cell)))
    
        # Creates two NumPy arrays: one with the coordinates of atoms in the
        # supercell and other with the element symbols of the atoms in it.
        scell_geom_dupl = np.loadtxt(cif_output, skiprows=2, usecols=(1, 2, 3), dtype=np.float64)
        scell_elem_dupl = np.loadtxt(cif_output, skiprows=2, usecols=(0), dtype="str")

        # Clean up duplicate coordiantes in the new NumPy arrays.
        scell_geom = scell_geom_dupl[unique_idx]
        scell_elem = scell_elem_dupl[unique_idx]

        print("         The number of unique coordinates in the supercell is now: {}\n".format(len(clean_cell)))

    else:
        # Creates two NumPy arrays: one with the coordinates of atoms in the
        # supercell and other with the element symbols of the atoms in it.
        scell_geom = np.loadtxt(cif_output, skiprows=2, usecols=(1, 2, 3), dtype=np.float64)
        scell_elem = np.loadtxt(cif_output, skiprows=2, usecols=(0), dtype="str")

    if verbose >= 2:
        print("Generating monomers for all complete molecules in the supercell:")

    # Distances will be handled in Bohr.
    scell_geom = scell_geom / qcel.constants.bohr2angstroms

    # Calculation of the supercell center as the midpoint of all
    # coordinates.
    scell_cntr = (np.max(scell_geom, axis=0) - np.min(scell_geom, axis=0))/2.0
    
    if verbose >= 2:
        print("\nCurrently, the center of the supercell is located at:")
        print("x = %10.5f" % (scell_cntr[0]))
        print("y = %10.5f" % (scell_cntr[1]))
        print("z = %10.5f" % (scell_cntr[2]))
        print("\nThe supercell coordinates will be translated and centered on the origin.")

    # Translate the supercell to the origin.
    scell_geom -= scell_cntr
    
    # Return a NumPy array with 3 numbers: the maximum of each x, y, and
    # z.
    scell_geom_max_coords = np.max(scell_geom, axis=0)

    return scell_geom_max_coords, scell_geom, scell_elem
# ======================================================================


# ======================================================================
def supercell2monomers(cif_output, r_cut_monomer, bfs_thresh, verbose=1):
    """Takes the supercell cartesian coordinates file produced by
    Read_CIF, and passes it to the `center_supercell()` function which
    translates the supercell to the origin.

    The centered supercell geometries and elements arrays are passed to
    the Breadth-First Search which returns all fragments found in the
    supercell.

    This function also takes a cutoff distance, measured from the
    origin, which is then used to decide if a monomer should be included
    or not based on its proximity to the origin.

    Returns a dictionary with all the fragments (monomers) that are
    located within the cutoff region.
    
    Arguments:
    <str> cif_output
        Name of the file with the cartesian coordinates of the
        supercell.
    <float> r_cut_monomer
        Cutoff value to include fragments found in the supercell in the
        dictionary of N-mers.
    <int> verbose
        Adjusts the level of detail of the printouts.

    Returns:
    <dict> nmers
        A dictionary populated with N-mers (monomers at this time) and
        their corresponding atributes.
    """

    # Centering the supercell.
    scell_geom_max_coords, scell_geom, scell_elem = center_supercell(cif_output, verbose)

    # Check if each of dimensions of the supercell satisfies the 
    # condition of being the twice as long as the cutoff. This helps
    # filtering out fragments that contain incomplete molecules located
    # at the edges of the supercell.
    if (r_cut_monomer / qcel.constants.bohr2angstroms) > np.min(scell_geom_max_coords):
        print("\nERROR: Cutoff (%3.2f A) longer than half the smallest dimension of the supercell (%3.2f A)." \
              % (r_cut_monomer, np.min(scell_geom_max_coords)*qcel.constants.bohr2angstroms))
        print("       Please increase the dimensions of the supercell to at least twice r_cut_monomer or reduce the lenght of the cutoff.\n")
        sys.exit()

    # Start the BFS timer.
    bfs_start_time = time.time()
    
    # Passes the supercell geometry and elements to the breadth-first
    # search algorithm of QCDB to obtain fragments.
    #fragments = BFS(scell_geom, scell_elem, None, bfs_thresh)
    fragments = psi4.driver.qcdb.bfs.BFS(scell_geom, scell_elem, None, bfs_thresh)

    # Stop the BFS timer.
    bfs_stop_time = time.time() - bfs_start_time
    
    # Two lists containing geometries and elements of a fragment.
    frag_geoms = [scell_geom[fr] for fr in fragments]
    frag_elems = [scell_elem[fr] for fr in fragments]

    # List containing the magnitude of the shortest position vector of
    # each fragment.
    frag_r_mins = []

    for frag in frag_geoms:
        # Get the minium of the norm of the position vectors of all
        # atoms in the fragment.
        r_min_atom = np.min(np.linalg.norm(frag, axis=1))
        frag_r_mins.append(r_min_atom)
    
    # Convert the list of minimums to a NumPy array
    frag_r_mins = np.array(frag_r_mins)

    # Array mapping the indices of the fragments sorted by r_min_atom to
    # the order in the frag_r_mins list.
    mapper = frag_r_mins.argsort()

    # A dictionary for storing N-mers will be created. 

    # Then, for each fragment in the frag_geoms list, an index will be 
    # defined based on the mapper list.
    
    # If the shortest position vector is shorter that the monomers 
    # cutoff, a new dictionary will be created to store all the 
    # information corresponding to such N-mer.
    nmers = {}

    for i in range(len(frag_geoms)):
        index = mapper[i]

        # Discard edges of the supercell and keeps the monomers within a
        # sphere of the cutoff radius.
        if frag_r_mins[index] <= (r_cut_monomer / qcel.constants.bohr2angstroms):
            name = "1mer-" + str(i)
            
            # The dictionary of one N-mer.
            nmers[name] = {}
            
            nmers[name]["monomers"] = [i]
            nmers[name]["elem"] = frag_elems[index]
            nmers[name]["coords"] = frag_geoms[index]
            nmers[name]["rmin"] = frag_r_mins[index]
            nmers[name]["delimiters"] = []
            nmers[name]["com"] = center_of_mass(frag_elems[index], frag_geoms[index])
            nmers[name]["priority_min"] = 0.0
            nmers[name]["priority_com"] = 0.0

    total_number_of_monomers = len(nmers.keys())

    if verbose >= 2:
        print("\nThe BFS algorithm found {} monomers in the supercell in {:.2f} s".format(total_number_of_monomers, bfs_stop_time))
    
    return nmers
# ======================================================================


# ======================================================================
def create_nmer(nmers, ref_monomer, other_monomers, verbose=1):
    """Takes a `nmers` dictionary, and two strings with the keys of a
    refrence monomer and other monomer.

    This function will merge both monomers and create a new N-mer.

    It returns a string `nm_new_name` and a `nm_new` dictionary
    containing the information of the newly generated N-mer.

    Arguments:

    Returns:
    <str> nm_new_name
        Name of the new N-mer, typically `n`mer-`i`+`j`+`k`+... where
        `n` is the order of the N-mer and `i`, `j`, ... the indices of
        the monomer that compose it.
    <dict> nm_new
        Dictionary of the new N-mer containing its corresponding
        information.
    """
    
    nm_new = {}

    # Monomers included in the new N-mer.
    nm_new_monomers = []
    nm_new_monomers.append(ref_monomer["monomers"][0])
    nm_new_monomers.extend(other_monomers)
    nm_new["monomers"] = nm_new_monomers
    
    # Elements and coordinates of the atoms in the new N-mer.
    nm_new_elem_arrays = []
    nm_new_coords_arrays = []
    nm_new["atoms_per_monomer"] = []

    for monomer in nm_new_monomers:
        name = "1mer-" + str(monomer)
        nm_new["atoms_per_monomer"].append(len(nmers[name]["elem"]))
        nm_new_elem_arrays.append(nmers[name]["elem"])
        nm_new_coords_arrays.append(nmers[name]["coords"])

    nm_new["elem"] = np.concatenate(nm_new_elem_arrays, axis=0)
    nm_new["coords"] = np.concatenate(nm_new_coords_arrays, axis=0)
    
    # Indices of each monomer in the array of elements and coordinates 
    # (For spliting N-mers into monomers).
    nm_new["delimiters"] = np.cumsum(nm_new["atoms_per_monomer"])
    
    # Norm of the shortest position vector of an atom in the new N-mer.
    # The reference is always closest to the origin so it is always the
    # minimum rmin
    nm_new["rmin"] = ref_monomer["rmin"]

    # Nuclear repulsion energy of the new N-mer.
    nm_new["nre"] = nre(nm_new["elem"], nm_new["coords"])

    # Chemical space eigenvalues of the new N-mer.
    nm_new["chsev"] = chemical_space(nm_new["elem"], nm_new["coords"])

    # Non additive N-body energy of the N-mer.
    nm_new["nambe"] = 0.0

    # Counter of the number of replicas of this N-mer.
    nm_new["replicas"] = 1

    # Contribution of unique structure to crystal lattice energy.
    nm_new["contrib"] = 0.0

    # Shortest separation vector between the atoms of the monomers in
    # the N-mer.
    nm_new["min_monomer_separations"] = []

    # Separation vector between center of mass of the monomers in the
    # N-mer.
    nm_new["com_monomer_separations"] = []

    for a in nm_new_monomers:
        a_name = "1mer-" + str(a)
        
        for b in nm_new_monomers:
            
            if b <= a:
                continue
            
            b_name = "1mer-" + str(b)

            distm, r_min = distance_matrix(nmers[a_name]["coords"], nmers[b_name]["coords"]) 
            nm_new["min_monomer_separations"].append(r_min)
            nm_new["com_monomer_separations"].append(np.linalg.norm(nmers[a_name]["com"] - nmers[b_name]["com"]))
    
    
    # Criterion to launch energy calculations.
    nm_new["priority_min"] = 0.0

    priority_min = 1.0
    
    for rmin in nm_new["min_monomer_separations"]:
        one_over_rmin3 = 1.0/rmin**3
        priority_min *= one_over_rmin3
    
    nm_new["priority_min"] = priority_min

    # Alternative criterion to launch energy calculations.
    nm_new["priority_com"] = 0.0

    priority_com = 1.0
    
    for rcom in nm_new["com_monomer_separations"]:
        one_over_rcom3 = 1.0/rcom**3
        priority_com *= one_over_rcom3
    
    nm_new["priority_com"] = priority_com

    # Key of the new N-mer in the nmers dictionary.
    nm_new_name = str(len(nm_new_monomers)) + "mer-" + "+".join(map(str, nm_new_monomers))

    return nm_new_name, nm_new
# ======================================================================


# ======================================================================
def center_of_mass(elems, geoms):
    """Takes the element symbols and coordinates of a set of atoms and
    computes the center of mass of the molecule.
    
    Arguments:
    <numpy.ndarray> elems
        Array of 1 column with the atomic symbols of each atom.
    <numpy.ndarray> geoms
        Array of 3 columns with the coordinates of the system.

    Returns:
    <numpy.ndarray> com
        Array of 3 numbers with the center of mass of the system.
    """
    
    com = np.array([0.0, 0.0, 0.0])
    total_mass = 0.0

    for at in range(len(geoms)):
        m = qcel.periodictable.to_mass(elems[at])
        com += m*geoms[at]
        total_mass += m
    
    com /= total_mass
    
    return com
# ======================================================================


# ======================================================================
def distance_matrix(a, b):
    """Euclidean distance matrix between rows of arrays `a` and `b`.
    Equivalent to `scipy.spatial.distance.cdist(a, b, 'euclidean')`.
    Returns a.shape[0] x b.shape[0] array. It also returns the shortest
    distance between the atoms of `a` and of `b`.
    """

    assert(a.shape[1] == b.shape[1])
    distm = np.zeros([a.shape[0], b.shape[0]])
    
    for i in range(a.shape[0]):
    
        for j in range(b.shape[0]):
            distm[i, j] = np.linalg.norm(a[i] - b[j])

    r_min = np.min(distm)

    return distm, r_min
# ======================================================================


# ======================================================================
def nre(elem, geom):
    """Takes two Numpy arrays, one with the element symbols and one with
    coordinates of a set of atoms and returns a float number with the 
    computed nuclear repulsion energy.
    """
    
    nre = 0.
    for at1 in range(geom.shape[0]):

        for at2 in range(at1):
            dist = np.linalg.norm(geom[at1] - geom[at2])
            nre += qcel.periodictable.to_Z(elem[at1]) * qcel.periodictable.to_Z(elem[at2]) / dist

    return nre
# ======================================================================


# ======================================================================
def chemical_space(elem, geom):
    """Takes the element symbols and coordinates of a set of atoms and
    computes the chemical space matrix, eigenvalues, and eigenvectors
    of the chemical system and returns a list with the sorted
    eigenvalues.

    Arguments:
    <numpy.ndarray> elems
        Array of 1 column with the atomic symbols of each atom.
    <numpy.ndarray> geoms
        Array of 3 columns with the coordinates of the system.

    Returns:
    <numpy.ndarray> chem_spc_eigen_values
        Array of x numbers with the sorted eigenvalues of the molecular
        system, where x is the number of atoms in the system.
    """

    # Get the number of atoms in the N-mer
    natoms = geom.shape[0]
    
    # Create a NumPy matrix
    M = np.zeros((natoms,natoms))

    # Iterate over atoms
    for i in range(natoms):
   
        # Get the charge of atom i
        charge_i = qcel.periodictable.to_Z(elem[i])
        
        # Fill the diagonal with the special polynomial from of:
        # DOI: 10.1103/PhysRevLett.108.058301     
        M[i,i] = 0.5 * np.power(charge_i, 2.4)
    
        for j in range(i):
            
            # Get the charge of atom j
            charge_j = qcel.periodictable.to_Z(elem[j])
            
            # Compute distance between i and j
            dist = np.linalg.norm(geom[i] - geom[j])
            
            # Compute Coulomb interaction between i and j
            ij_elem = charge_i * charge_j / dist

            M[i,j] = ij_elem
            # Symmetric Matrix 
            M[j,i] = ij_elem

    # Solve the eigenvalue problem
    eigenvalues, eigenvectors = np.linalg.eig(M)

    # Eigenvalues must be in list to use 'sort'
    sorted_eigenvalues = list(eigenvalues)
    
    # Sort the eigenvalues in order of "decreasing absolute value".
    # This first sort is done to guarantee the same sort in the case 
    # that two eigenvalues are the same magnitude but different sign.
    sorted_eigenvalues.sort(key = lambda x: -x)
    sorted_eigenvalues.sort(key = lambda x: -abs(x))
    
    # Cast back to a NumPy array.
    chem_spc_eigen_values = np.array(sorted_eigenvalues)

    return chem_spc_eigen_values
# ======================================================================


# ======================================================================
def build_nmer(nmers, total_monomers, nmer_type, nmer_separation_cutoff, coms_separation_cutoff, uniq_filter, verbose=1):
    """Takes a float indicating a cutoff distance in Angstroms.
    Returns dimer files that pass through the filter.
    """

    # Function reused for different types of N-mers.
    
    build_nmers_start_time = time.time()

    if nmer_type == "dimers":
        #nm_dictname_pattern = "1mer-"
        num_monomers = 2

    elif nmer_type == "trimers":
        #nm_dictname_pattern = "2mer-"
        num_monomers = 3

    elif nmer_type == "tetramers":
        #nm_dictname_pattern = "3mer-"
        num_monomers = 4

    elif nmer_type == "pentamers":
        #nm_dictname_pattern = "4mer-"
        num_monomers = 5
    
    else:
        print("\nERROR: The N-mer type must be defined as 'dimers', 'trimers', 'tetramers', or 'pentamer'.")

    if verbose >= 2:
        print("")

    counter_new_nmers = 0 # Number of new N-mers generated
    counter_dscrd_sep = 0 # N-mers filtered out by atomic separation
    counter_dscrd_com = 0 # N-mers filtered out by COM separation
    counter_dscrd_rep = 0 # N-mers filtered out as a replicas

    new_nmers = {}

    # TODO: Support for crystals with more than one molecule in the
    #       primitive unit cell.
    num_ref_monomers = 1

    for ref_monomer_idx in range(num_ref_monomers):
        monomer_key = "1mer-" + str(ref_monomer_idx)
        ref_monomer = nmers[monomer_key]

        for other_monomers in itertools.combinations(range(ref_monomer_idx+1, total_monomers), num_monomers-1):

            # Start N-mer building timer.
            nmer_start_time = time.time()

            new_nmer_name, new_nmer = create_nmer(nmers, ref_monomer, other_monomers, verbose)

            max_mon_sep = max(new_nmer["min_monomer_separations"])
            max_com_sep = max(new_nmer["com_monomer_separations"])

            if max_mon_sep > (nmer_separation_cutoff / qcel.constants.bohr2angstroms):
                
                if verbose >= 2:
                    
                    # Stop N-mer building timer.
                    nmer_stop_time = time.time() - nmer_start_time

                    print(
                        "{} discarded in {:.2f} s: Maximum separation between closest atoms of different monomers is {:3.2f} A, longer than cutoff {:3.2f} A.".format(
                            new_nmer_name, nmer_stop_time, max_mon_sep*qcel.constants.bohr2angstroms, nmer_separation_cutoff))
                        
                    counter_dscrd_sep += 1
            
            elif max_com_sep > (coms_separation_cutoff / qcel.constants.bohr2angstroms):
                
                if verbose >= 2:
                    
                    # Stop N-mer building timer.
                    nmer_stop_time = time.time() - nmer_start_time
                    print(
                        "{} discarded in {:.2f} s: Maximum separation between closest COMs of different monomers is {:3.2f} A, longer than cutoff {:3.2f} A.".format(
                            new_nmer_name, nmer_stop_time, max_com_sep*qcel.constants.bohr2angstroms, coms_separation_cutoff))

                    counter_dscrd_com += 1

            else:
                found_duplicate = False
                
                #if False:
                if uniq_filter == "Dreamaligner":
                    chemical_space = False
                    dreamliner = True
                
                else:
                    chemical_space = True
                    dreamliner = False

                nre_filter_ran = False
                chsev_filter_ran = False
                rmsd_filter_ran = False

                for kexisting, existing in new_nmers.items():

                    # Nuclear repulsion energy filter.
                    nre_diff = abs(existing["nre"] - new_nmer["nre"])
                    nre_filter_ran = True

                    # If NRE difference is large, this is a new N-mer.
                    # Threfore reset posterior filters ran flags, there
                    # is no need to run further filters than NRE.
                    if nre_diff > 1.e-5:

                        chsev_filter_ran = False
                        rmsd_filter_ran = False                        
                    
                    # The NRE is small, proceed to next filters.
                    else:
                        # Chemical space eigenvalues filter.
                        if chemical_space == True:

                            chsev_diff = np.linalg.norm(existing["chsev"] - new_nmer["chsev"])
                            chsev_filter_ran = True

                            if chsev_diff < 1.e-3:
                                found_duplicate = True

                                if verbose >= 2:

                                    # Stop N-mer building timer.
                                    nmer_stop_time = time.time() - nmer_start_time

                                    print("{} discarded in {:.2f} s: Replica of {}. NRE difference is {:.1e} and ChSEV difference is {:.1e}.".format(
                                        new_nmer_name, nmer_stop_time, kexisting, nre_diff, chsev_diff))

                                existing["replicas"] += 1
                                counter_dscrd_rep += 1

                                break

                        # RMSD filter.
                        if dreamliner == True:
                        
                            # Block B787 printout
                            sys.stdout = open(os.devnull, 'w')
                            
                            # Call the dreamliner from QCDB.
                            rmsd, mill = qcel.molutil.B787(rgeom=existing["coords"], cgeom=new_nmer["coords"], runiq=existing["elem"], cuniq=new_nmer["elem"], run_mirror=True, verbose=2)
                            rmsd_filter_ran = True

                            # Reanable printout
                            sys.stdout = sys.__stdout__

                            if rmsd < 1.e-3:
                                found_duplicate = True

                                if verbose >= 2:
                                    
                                    # Stop N-mer building timer.
                                    nmer_stop_time = time.time() - nmer_start_time

                                    print("{} discarded in {:.2f} s: Replica of {}. NRE difference is {:.1e} and RMSD is {:.1e}.".format(
                                        new_nmer_name, nmer_stop_time, kexisting, nre_diff, rmsd))

                                existing["replicas"] += 1
                                counter_dscrd_rep += 1

                                break

                if not found_duplicate:
                    new_nmers[new_nmer_name] = new_nmer

                    if verbose >= 2:

                        # Stop N-mer building timer.
                        nmer_stop_time = time.time() - nmer_start_time

                        if nre_filter_ran == True:
                            
                            if chsev_filter_ran == True:
                                print("{} generated in {:.2f} s: New N-mer NRE difference is {:.1e}, lowest ChSEV difference found is {:.1e}.".format(
                                    new_nmer_name, nmer_stop_time, nre_diff, chsev_diff))
                        
                            if rmsd_filter_ran == True:
                                print("{} generated in {:.2f} s: New N-mer NRE difference is {:.1e}, lowest RMSD found is {:.1e}.".format(
                                    new_nmer_name, nmer_stop_time, nre_diff, rmsd))

                            if chsev_filter_ran == False and rmsd_filter_ran == False:
                                print("{} generated in {:.2f} s: New N-mer NRE difference is {:.1e}.".format(new_nmer_name, nmer_stop_time, nre_diff))

                        else:
                            print("{} generated in {:.2f} s: New N-mer NRE is {:.12f}.".format(new_nmer_name, nmer_stop_time, new_nmer["nre"]))

                    counter_new_nmers += 1
    
    if verbose >= 2:
        build_nmers_stop_time = time.time() - build_nmers_start_time
        print("\n{} unique {} were found and generated in {:.1f} s.".format(counter_new_nmers, nmer_type, build_nmers_stop_time))

    if verbose >= 2:
        print("\n{} {} did not meet the atomic separation cutoff and were discarded.".format(counter_dscrd_sep, nmer_type))
        print("{} {} did not meet the center of mass separation cutoff and were discarded.".format(counter_dscrd_com, nmer_type))
        print("{} {} were duplicates of another dimer and were discarded.".format(counter_dscrd_rep, nmer_type))
    
    nmers.update(new_nmers)

    return nmers
# ======================================================================


# ======================================================================
def nmer2psiapimol(nmers, keynmer, nmer, verbose=0):
    """Takes the `nmers` dictionary; `keynmer`, the key of a given N-mer of
    the N-mers dictionary; and `nmer`, its corresponding dictionary.
    Returns a string `psi_api_molecule` that defines a molecule in PSI4
    API mode. This string starts with three lines specifying the use of
    atomic units, avoidance of center of mass translation, and
    avoidance of reorientation. The next lines contain the element
    symbol and coordinates, for each atom in the N-mer, separating each
    monomer with lines containing a double hyphen.
    """

    psi_api_molecule = "\nunits = au\n"
    psi_api_molecule += "no_com\n"
    psi_api_molecule += "no_reorient\n"

    for at in range(nmer["coords"].shape[0]):

        if at in nmer["delimiters"]:
            psi_api_molecule += "--\n"

        psi_api_molecule += "  {:6} {:16.8f} {:16.8f} {:16.8f} \n".format(nmer["elem"][at], nmer["coords"][at][0], nmer["coords"][at][1], nmer["coords"][at][2])

    return psi_api_molecule
# ======================================================================


# ======================================================================
def monomer2makefp(cif_output, nmer, verbose=0):
    """.
    """

    makefp_input =  "! CrystaLattE: GAMESS input for generation of EFP parameters\n"
    makefp_input += " $contrl units=angs local=boys runtyp=makefp\n"
    makefp_input += "   mult=1 icharg=0 coord=cart icut=11 $end\n"
    makefp_input += " $system timlim=99999 mwords=200 $end\n"
    makefp_input += " $scf dirscf=.t. soscf=.f. diis=.t. conv=1.0d-06 $end\n"
    makefp_input += " $basis gbasis=n311 ngauss=6 npfunc=2 ndfunc=3 nffunc=1\n"
    makefp_input += "   diffs=.t. diffsp=.t. $end\n"
    makefp_input += " $makefp pol=.t. disp=.t. exrep=.t. chtr=.f. $end\n"
    makefp_input += " $stone\n"
    makefp_input += "   bigexp=0.0\n"
    makefp_input += " $end\n"
    makefp_input += " $data\n"
    makefp_input += "  {}\n".format(cif_output.split(".")[0])
    makefp_input += "c1\n"

    for at in range(nmer["coords"].shape[0]):

        if at in nmer["delimiters"]:
            makefp_input += "--\n"
        
        makefp_input += " {:6} {:4.1f} {:16.8f} {:16.8f} {:16.8f} \n".format(nmer["elem"][at] + str(at + 1), qcel.periodictable.to_Z(nmer["elem"][at]), nmer["coords"][at][0], nmer["coords"][at][1], nmer["coords"][at][2])

    makefp_input += " $end"

    owd = os.getcwd()
    makefp_folder = cif_output[:-4]

    try:
        os.mkdir(makefp_folder)
    
    except FileExistsError:
        pass

    os.chdir(makefp_folder)
    makefp_filename = cif_output.split(".")[0] + "-makefp-lb.inp"

    with open(makefp_filename, "w") as makefp_f:
        
        for line in makefp_input:
            makefp_f.write(line)

    os.chdir(owd)
# ======================================================================


# ======================================================================
def nmer2psithon(cif_output, nmers, keynmer, nmer, rminseps, rcomseps, psi4_method, psi4_bsse, psi4_memory, verbose=0):
    """.
    """
    
    psithon_input =  "# PSI4 file produced by CrystaLattE\n\n"
    psithon_input += "# Psithon input for N-mer:      {}\n".format(keynmer)
    psithon_input += "# Number of replicas:           {}\n".format(nmer["replicas"])
    psithon_input += "# Priority index for input:     {:12.8e}\n".format(nmer["priority_min"])
    psithon_input += "# Minimum monomer separations:  {}\n".format(rminseps.lstrip(" "))
    psithon_input += "# COM priority index for input: {:12.8e}\n".format(nmer["priority_com"])
    psithon_input += "# Minimum COM separations:      {}\n".format(rcomseps.lstrip(" "))
    
    psithon_input += "\nmemory {}\n".format(psi4_memory)
    
    mymol = keynmer.replace("2mer", "Dimer").replace("3mer", "Trimer").replace("4mer", "Tetramer").replace("5mer", "Pentamer").replace("-", "_").replace("+", "_")
    psithon_input += "\nmolecule {} {{\n".format(mymol)
    
    for at in range(nmer["coords"].shape[0]):

        if at in nmer["delimiters"]:
            psithon_input += "--\n"
        
        psithon_input += "  {:6} {:16.8f} {:16.8f} {:16.8f} \n".format(nmer["elem"][at], nmer["coords"][at][0], nmer["coords"][at][1], nmer["coords"][at][2])

    psithon_input += "units = au\n"
    psithon_input += "no_com\n"
    psithon_input += "no_reorient\n"
    
    psithon_input += "}\n"

    psithon_input += "\nset {\n"
    psithon_input += "  e_convergence 9\n"
    psithon_input += "  scf_type df\n"
    psithon_input += "  mp2_type df\n" 
    psithon_input += "  cc_type df\n"
    psithon_input += "  freeze_core true\n"
    psithon_input += "}\n"

    # Hartree-Fock is called with the 'scf' string in Psithon mode.
    if psi4_method.lower().startswith("hf"):
        psithon_method = "scf" + psi4_method[2:]

    else:
        psithon_method = psi4_method

    psithon_input += "\nenergy('{}', bsse_type = '{}')\n".format(psithon_method, psi4_bsse)
    psithon_input += "\n"

    owd = os.getcwd()
    psithon_folder = cif_output[:-4]

    try:
        os.mkdir(psithon_folder)
    
    except FileExistsError:
        pass

    os.chdir(psithon_folder)
    psithon_filename = keynmer + ".in"

    with open(psithon_filename, "w") as psithon_f:
        
        for line in psithon_input:
            psithon_f.write(line)

    os.chdir(owd)
# ======================================================================


# ======================================================================
def psi4api_energies(cif_output, nmers, keynmer, nmer, cpus, cle_run_type, psi4_method, psi4_bsse, psi4_memory, verbose=0):
    """
    Arguments:
    
    """

    # If the output is going to be kept, setup the filename.
    if "quiet" or "test" in cle_run_type:
        psi4.core.be_quiet()
        
    # If the output is not kept, do not print to screen.
    else:
        owd = os.getcwd()

        p4folder = cif_output[:-4]
        
        try:
            os.mkdir(p4folder)
        except FileExistsError:
            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
            print("WARNING: A folder with the same name as the CIF file already exists.")
            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

        os.chdir(p4folder)
        
        p4out = keynmer + ".dat"
        psi4.core.set_output_file(p4out)
        
        os.chdir(owd)

    psi_api_molecule = nmer2psiapimol(nmers, keynmer, nmer, verbose)
    mymol = psi4.geometry(psi_api_molecule)
    
    # Set the number of threads to run Psi4.
    psi4.core.set_num_threads(cpus)
    
    # Set Psi4 memory.
    psi4.set_memory(psi4_memory)
    
    # Set Psi4 block of options.
    psi4.set_options({'scf_type': 'df', 'mp2_type': 'df', 'cc_type': 'df', 'freeze_core': 'true', 'e_convergence': '9'})
    
    # Execute Psi4 energy calculations, unless running on test mode.

    # Example:  psi4.energy('MP2/aug-cc-pV[D,T]Z', molecule=mymol, bsse_type=['cp', 'nocp', 'vmfc'])
    #           psi4.energy('HF/STO-3G', molecule=mymol, bsse_type=['nocp'])
    #           psi4.energy('MP2/aug-cc-pVDZ', molecule=mymol, bsse_type=['vmfc'])
    
    if "test" not in cle_run_type:
        psi4.energy(psi4_method, molecule=mymol, bsse_type=[psi4_bsse])
    
        # This is how you print out all the variables to figure out errors.
        #import pprint
        #pprint.pprint(psi4.core.variables())

        # Get the non-additive n-body contribution, exclusive of all
        # previous-body interactions.
        varstring = "{}-CORRECTED {}-BODY INTERACTION ENERGY".format(psi4_bsse.upper(), str(len(nmer["monomers"])))
    
        n_body_energy = psi4.core.variable(varstring)
    
        if len(nmer["monomers"]) > 2:
            varstring = "{}-CORRECTED {}-BODY INTERACTION ENERGY".format(psi4_bsse.upper(), str(len(nmer["monomers"]) - 1))
            n_minus_1_body_energy = psi4.core.variable(varstring)
            nmer["nambe"] = n_body_energy - n_minus_1_body_energy
        
        else:
            nmer["nambe"] = n_body_energy
# ======================================================================


# ======================================================================
def cle_manager(cif_output, nmers, cle_run_type, psi4_method, psi4_bsse, psi4_memory, verbose=0):
    """Manages which mode of CrystaLattE calculation will be employed.
    
    Arguments:
    <dict> nmers
        Dictionary containing dictionaries for each N-mer in the 
        system.
    <list> cle_run_type
        List of keywords indicating the modes in which CrystaLattE is
        going to be run.
    <str> psi4_method
        Method and basis set for the energy calculation, separated by a
        slash.
    <str> psi4_bsse
        Method for correction of the basis set superposition error.
    <str> psi4_memory
        Memory allocation for PSI4, written as `500 MB`, or `60 GB`.
    <int> verbose
        Adjusts the level of detail of the printouts.

    Returns:
    <float> crystal_lattice_energy
        The value of the accumulated crystal lattice energy in atomic
        units.
    <list> results
        A summary of the results for all N-mers including energies,
        replicas, contributions, cumulative lattice energy, priority,
        and minimum atomic separations; for the print_results function.
    """

    crystal_lattice_energy = 0.0
    results = []
    
    # Get the keys of the N-mers dictionary, and put them on a list.
    nmer_keys = list(nmers.keys())

    # Sort the list in decreasing priority order.
    nmer_keys.sort(key = lambda x: -nmers[x]['priority_min'])
    #nmer_keys.sort(key = lambda x: -nmers[x]['rmin'])

    # The next line was replaced to trigger the calculations in order.
    #for keynmer, nmer in nmers.items():
    for keynmer in nmer_keys:

        nmer = nmers[keynmer]
        # Energies are not calculated for monomers. Rigid body approximation.
        if len(nmer["monomers"]) == 1:
            continue

        # Find out the number of CPUs in the local system.
        cpus = multiprocessing.cpu_count()
        
        # Start wall-clock timer.
        energies_start = time.time()
        
        # Generate a string with an ordered list of minimum separations
        # between atoms belonging to different monomers.
        rminseps = ""
       
        nmer_min_monomer_separations = nmer["min_monomer_separations"]
        nmer_min_monomer_separations.sort()
        
        for r in nmer_min_monomer_separations:
            rminseps += "{:6.3f} ".format(r * qcel.constants.bohr2angstroms)
        
        # Generate a string with an ordered list of minimum separations
        # between the center of masses of the monomers.
        rcomseps = ""
       
        nmer_min_com_separations = nmer["com_monomer_separations"]
        nmer_min_com_separations.sort()
        
        for r in nmer_min_com_separations:
            rcomseps += "{:6.3f} ".format(r * qcel.constants.bohr2angstroms)

        # Produce Psithon inputs
        if "psithon" in cle_run_type:
            nmer2psithon(cif_output, nmers, keynmer, nmer, rminseps, rcomseps, psi4_method, psi4_bsse, psi4_memory, verbose)

        # Run energies in PSI4 API.
        else:
            psi4api_energies(cif_output, nmers, keynmer, nmer, cpus, cle_run_type, psi4_method, psi4_bsse, psi4_memory, verbose)

        nmer["contrib"] = nmer["nambe"] * nmer["replicas"] / float(len(nmer["monomers"]))
        crystal_lattice_energy += nmer["contrib"]

        # Stop wall-clock timer.
        energies_end = time.time()
        
        # Calculate execution time.
        energies_wallclock = energies_end - energies_start
        
        if "test" or "psithon" in cle_run_type:
            nmer_result = "{:26} | {:>12} | {:>4} | {:>12} | {:>13} | {:12.6e} | {}".format(
                    keynmer,
                    "Not Computed", 
                    nmer["replicas"],
                    "Not Computed",
                    "Not Computed",
                    nmer["priority_min"],
                    rminseps)

        else:
            nmer_result = "{:26} | {:>12.8f} | {:>4} | {:>12.8f} | {:>13.8f} | {:12.6e} | {}".format(
                    keynmer, 
                    nmer["nambe"] * qcel.constants.hartree2kcalmol * qcel.constants.cal2J, 
                    nmer["replicas"], 
                    nmer["contrib"] * qcel.constants.hartree2kcalmol * qcel.constants.cal2J,
                    crystal_lattice_energy * qcel.constants.hartree2kcalmol * qcel.constants.cal2J,
                    nmer["priority_min"],
                    rminseps)
        
        results.append(nmer_result)

        if verbose >= 2:

            if "psi4api" in cle_run_type:
                print("{} elapsed {:.2f} s processing on {} threads. Cumulative Lattice Energy = {:9.8f} kJ/mol".format(
                    keynmer, 
                    energies_wallclock, 
                    cpus, 
                    crystal_lattice_energy * qcel.constants.hartree2kcalmol * qcel.constants.cal2J))

            if "psithon" in cle_run_type:
                print("{} written.".format(keynmer))
        
    return crystal_lattice_energy, results
# ======================================================================


# ======================================================================
def print_header():
    """.
    """

    print("")
    print("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n")
    print("                              CrystaLattE                              \n")
    print("  The tool for the automated calculation of crystal lattice energies.  \n")
    print("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n")
# ======================================================================


# ======================================================================
def print_results(results, crystal_lattice_energy, verbose=0):
    """Prints a summary of the energy results at the end of the
    execution.
    
    Arguments:
    <int> verbose
        Adjusts the level of detail of the printouts.
    """

    if verbose >= 1:
        print("Summary of results:")
        print("---------------------------+--------------+------+--------------+---------------+--------------+----------------{}".format("-"*(shutil.get_terminal_size().columns - 112)))
        print("                           | Non-Additive | Num. |        N-mer | Partial Crys. |  Calculation | Minimum Monomer")
        print("N-mer Name                 |    MB Energy | Rep. | Contribution | Lattice Ener. |     Priority | Separations")
        print("                           |     (kJ/mol) |  (#) |     (kJ/mol) |      (kJ/mol) | (Arb. Units) | (A)")
        print("---------------------------+--------------+------+--------------+---------------+--------------+----------------{}".format("-"*(shutil.get_terminal_size().columns - 112)))
        for result in results:
            print(result)
        print("---------------------------+--------------+------+--------------+---------------+--------------+----------------{}\n".format("-"*(shutil.get_terminal_size().columns - 112)))
        #print("Crystal Lattice Energy (Eh)       = {:5.8f}".format(crystal_lattice_energy))
        print("Crystal Lattice Energy (kJ/mol)   = {:9.8f}".format(crystal_lattice_energy * qcel.constants.hartree2kcalmol * qcel.constants.cal2J))
        print("Crystal Lattice Energy (kcal/mol) = {:9.8f}\n".format(crystal_lattice_energy * qcel.constants.hartree2kcalmol))
# ======================================================================


# ======================================================================
def print_end_msg(start, verbose=0):
    """Prints a success message and timing information at the end of the
    execution.
       
    Arguments:
    <int> verbose
        Adjusts the level of detail of the printouts.
    """
    
    if verbose >= 1:
        print("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n")
        print("Execution terminated succesfully.")
        print("Total elapsed wall-clock time: {:.2f} s\n".format(time.time() - start))
        print("Thank you for using CrystaLattE.\n")
        print("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n")

# ======================================================================


# ======================================================================
def main(cif_input, cif_output="sc.xyz", cif_a=5, cif_b=5, cif_c=5, bfs_thresh=1.2, uniq_filter="ChSEV", nmers_up_to=2, r_cut_com=10.0, r_cut_monomer=12.0, r_cut_dimer=10.0, r_cut_trimer=8.0, r_cut_tetramer=6.0, r_cut_pentamer=4.0, cle_run_type=["test"], psi4_method="HF/STO-3G", psi4_bsse="cp", psi4_memory="500 MB", verbose=1):
    """Takes a CIF file and computes the crystal lattice energy using a
    many-body expansion approach.
    """
   
    # Start counting time.
    start = time.time()

    # Read a CIF file and generate the unit cell.
    cif_arguments = cif_driver(cif_input, cif_output, cif_a, cif_b, cif_c, verbose)
    cif_main(cif_arguments)
    
    # Read the output of the automatic fragmentation.
    nmers = supercell2monomers(cif_output, r_cut_monomer, bfs_thresh, verbose)
    total_monomers = len(nmers)

    # If makefp mode requested, produce a MAKEFP input file for gamess
    # and exit.
    if "makefp" in cle_run_type:
        monomer2makefp(cif_output, nmers["1mer-0"], verbose)
        print("\nThe makefp file for {} has been created.\n".format(cif_output[:-4]))
        
        if verbose >= 1:
            print("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n")
            print("Execution terminated succesfully.")
            print("Total elapsed wall-clock time: {:.2f} s\n".format(time.time() - start))
            print("Thank you for using CrystaLattE.\n")
            print("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n")

        sys.exit(0)

    # Loop through all monomers and monomers in the central unit cell
    # to generate dimers with at least one monomer in the central cell.
    # Then loop through all existing N-mers and generate higher-order 
    # N-mers with all monomers.

    if nmers_up_to < 2:
        print("\nERROR: CrystaLattE is designed to use at least dimers.")
        print("       Please use 2 <= nmer_up_to < 5.")
        sys.exit()
    
    if nmers_up_to >= 2:
        
        if verbose >= 2:
            print("\nMerging monomers with monomers to obtain dimers.")
        build_nmer(nmers, total_monomers, "dimers", r_cut_dimer, r_cut_com, uniq_filter, verbose)

    if nmers_up_to >= 3:
        
        if verbose >= 2:
            print("\nMerging dimers with monomers to obtain trimers.")
        build_nmer(nmers, total_monomers, "trimers", r_cut_trimer, r_cut_com, uniq_filter, verbose)

    if nmers_up_to >= 4:
        
        if verbose >= 2:
            print("\nMerging trimers with monomers to obtain tetramers.")
        build_nmer(nmers, total_monomers, "tetramers", r_cut_tetramer, r_cut_com, uniq_filter, verbose)

    if nmers_up_to == 5:
        
        if verbose >= 2:
            print("\nMerging tetramers with monomers to obtain pentamers.")
        build_nmer(nmers, total_monomers, "pentamers", r_cut_pentamer, r_cut_com, uniq_filter, verbose)

    if nmers_up_to > 5:
        print("\nERROR: The current implementation of CrystaLattE is limited to pentamers.")
        print("       Please use 2 <= nmer_up_to < 5.")
        sys.exit()
   
    if verbose >= 2:
        
        if "psi4api" in cle_run_type:
            print ("\nComputing interaction energies of N-mers:")
        
        if "psithon" in cle_run_type:
            print ("\nWriting N-mer coordinates to Psithon input files:")

    crystal_lattice_energy, results = cle_manager(cif_output, nmers, cle_run_type, psi4_method, psi4_bsse, psi4_memory, verbose)
    # ------------------------------------------------------------------

    if verbose >= 2:
        print("")

    # Print the final results.
    print_results(results, crystal_lattice_energy, verbose)
    
    # Print exit message and timings information.
    print_end_msg(start, verbose)

    return nmers, crystal_lattice_energy

# ======================================================================


if __name__ == "__main__":

    # Hard-coded Test
    if "crystalatte.py" in sys.argv[-1]:
        
        print_header()

        main(   cif_input="../Tests/Ammonia/Ammonia.cif",
                cif_output="../Tests/Ammonia/Ammonia.xyz",
                cif_a=3,
                cif_b=3,
                cif_c=3,
                bfs_thresh=1.2,
                uniq_filter="ChSEV",
                nmers_up_to=2,
                r_cut_com=6.5,
                r_cut_monomer=3.5,
                r_cut_dimer=2.6,
                r_cut_trimer=3.7,
                r_cut_tetramer=3.7,
                r_cut_pentamer=6.1,
                cle_run_type=["psi4api"],
                psi4_method="HF/STO-3G",
                psi4_bsse="nocp",
                psi4_memory="500 MB",
                verbose=2)

        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        print("WARNING: No input was provided. The previous execution was just a test.")
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        print("\nCrystaLattE execution command:  ./crystalatte.py YourInput.cle\n")

    # Normal execution using an input file.
    else:
        input_parser(sys.argv[-1])
