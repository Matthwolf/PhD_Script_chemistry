# Import required libraries
import pymatgen as pmg
from pymatgen.io.gaussian import GaussianOutput

# Get the Gaussian output file name
filename = input('Enter the Gaussian output file name (including .out): ')

# Read the Gaussian output file
gaussian_output = GaussianOutput(filename)

# Get the energies, final structure, Mulliken charges, and structures from the Gaussian output
energies = gaussian_output.energies
final_structure = gaussian_output.final_structure
mulliken_charges = gaussian_output.Mulliken_charges
structures = gaussian_output.structures

# Get the basis set and functional from the Gaussian output
basis_set = gaussian_output.basis_set
functional = gaussian_output.functional

# Count the number of different atom types (C, O, H) in the final structure
carbon_count = 0
oxygen_count = 0
hydrogen_count = 0
for site in final_structure:
    if str(site.specie) == "C" : carbon_count = 1
    if str(site.specie) == "O" : oxygen_count = 1
    if str(site.specie) == "H" : hydrogen_count = 1

atom_type_count = carbon_count + oxygen_count + hydrogen_count

# Get the number of atoms in the molecule
atom_count = len(structures[1])

# Set the size of the box
box_size = 10.0

# Set the indices for C, O, and H
carbon_index = 1
oxygen_index = 2
hydrogen_index = 3

#incrementation

i=1

# Create a boxed structure of the molecule
mol_box = final_structure.get_boxed_structure(box_size, box_size, box_size)

# Convert the box size to an integer and a string
box_int = int(box_size)
box_str = str(box_int)

# Create the name of the data file
data_file = "Molecule_" + box_str + ".data"

# Write the data file
with open(data_file, "w") as file:
    file.write("LAMMPS DATA FILE. " + functional + "/" + basis_set + "\n")
    file.write(str(atom_count) + " atoms\n")
    file.write("0 bonds\n")
    file.write("0 angles\n")
    file.write("0 dihedrals\n")
    file.write("0 impropers\n")
    file.write(str(atom_type_count) + " atom types\n")
    file.write("0 bond types\n")
    file.write("0 angle types\n")
    file.write("0 dihedral types\n")
    file.write("0 improper types\n")
    file.write("0.0 " + str(box_size) + " xlo xhi\n")
    file.write("0.0 " + str(box_size) + " ylo yhi\n")
    file.write("0.0 " + str(box_size) + " zlo zhi\n")
    ###WARING THIS VERSION IS FOR MOLECULE CONTAINING C O H keep the same index for atome type and mass
    if carbon_count == 1 : file.write("\n"+"#"+" "+str(carbon_index)+" "+"C")
    if oxygen_count == 1 : file.write("\n"+"#"+" "+str(oxygen_index)+" "+"O")
    if hydrogen_count == 1 : file.write("\n"+"#"+" "+str(hydrogen_index)+" "+"H")
    file.write("\n")
    file.write("\n"+"Masses")
    file.write("\n")
    if carbon_count == 1 : file.write("\n"+str(carbon_index)+" "+"12.010700 # C")
    if oxygen_count == 1 : file.write("\n"+str(oxygen_index)+" "+"15.99940 # O")
    if hydrogen_count == 1 : file.write("\n"+str(hydrogen_index)+" "+"1.007940 # H")
    file.write("\n")
    file.write("\n"+"Atoms # charge")
    file.write("\n")
    for site in mol_box:
        if site.specie.symbol == "C" :
            if mulliken_charges[i][1] > 0:
                file.write("\n"+str(i)+" "+str(carbon_index)+" "+" "+str(mulliken_charges[i][1])+"".join(["%8.3f" % x for x in mol_box.cart_coords[i-1]])+" "+"#"+" "+site.specie.symbol)
            else:
                file.write("\n"+str(i)+" "+str(carbon_index)+" "+str(mulliken_charges[i][1])+"".join(["%8.3f" % x for x in mol_box.cart_coords[i-1]])+" "+"#"+" "+site.specie.symbol)
        if site.specie.symbol == "O" :
            if mulliken_charges[i][1] > 0:
                file.write("\n"+str(i)+" "+str(oxygen_index)+" "+" "+str(mulliken_charges[i][1])+"".join(["%8.3f" % x for x in mol_box.cart_coords[i-1]])+" "+"#"+" "+site.specie.symbol)
            else:
                file.write("\n"+str(i)+" "+str(oxygen_index)+" "+str(mulliken_charges[i][1])+"".join(["%8.3f" % x for x in mol_box.cart_coords[i-1]])+" "+"#"+" "+site.specie.symbol)
        if site.specie.symbol == "H" :
            if mulliken_charges[i][1] > 0:
                file.write("\n"+str(i)+" "+str(hydrogen_index)+" "+" "+str(mulliken_charges[i][1])+"".join(["%8.3f" % x for x in mol_box.cart_coords[i-1]])+" "+"#"+" "+site.specie.symbol)
            else:
                file.write("\n"+str(i)+" "+str(hydrogen_index)+" "+str(mulliken_charges[i][1])+"".join(["%8.3f" % x for x in mol_box.cart_coords[i-1]])+" "+"#"+" "+site.specie.symbol)
        i = i +1
    file.write("\n"+" ")