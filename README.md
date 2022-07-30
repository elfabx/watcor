# WatCor: utility to adjust water sites in a Gromacs .gro file

This utility replaces existing water molecules in a Gromacs coordinate file
(.gro format) with coordinates that match the specified water model.
This may involve adding / removing virtual sites.

The water O atom coordinates are conserved. New H atoms are created so that the
H-O-H angle matches the selected water model and the bisector of this agle
points to the same direction as in the original file. Virtual (or lone pair)
sites are then generated geometrically from the new O and H coordinates.
The aim of this procedure is to keep any H-bonding arrangement in the source
file intact as much as possible when using the selected model.
It will still be necessary to (re)optimise and (re)equilibrate the resulting
structure. It is left to the user to change the water model in the topology
file.

For a normal biomolecular simulation you are unlikely to need this tool.
It is probably easier to first remove all solvent molecules and then solvate
the biomolecule using the new water model. This tool is meant to help those
who model water with a defined structure, such as hydrates or ice.

## Building

From the source directory:
```
mkdir build
cd build
cmake ..
make
```

## Usage

Move the executable to your working directory or intall in the path:

```
./watcor -m model input.gro output.gro
```

If `-m model` is omitted tip3p is assumed. If `output.gro` is not given, the
new coordinate file is written to stdout.

Currently the following models are included:

- tip3p (default)
- tip3p-fb
- spc/e
- spc/fw
- spc/eb
- opc3
- opc
- tip4p
- tip4p-ew
- tip4p-fb
- tip5p
- tip5p-e


