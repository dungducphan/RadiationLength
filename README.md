# RadiationLength
A program to calculate the radiation length of materials. This program used the command line argument parsing tool [`cxxopts`](https://github.com/jarro2783/cxxopts) by author [jarro2783](https://github.com/jarro2783).

```
Usage:
  radLength [OPTION...]

	If --file option is chosen, the program will read the input data file and calculate the radiation length
	of the composite material described by the data file. The data file is an ASCII file, formatted in three
	columns respectively containing the atomic number (Z), the atomic mass (A) and the mass fraction (%).
	An example of the file can be as following:

		1	 1.008		30
		6	12.012		60
		8	16.002		10

	This file describes a material composed of H, C and O with mass fraction of 30%, 60% and 10%, respectively.


OPTIONS:

      --file FILE               Specify material composition data file,
                                --material and --atomicid options will be ignored.
      --material MATERIAL_NAME  Specify a pre-defined material in the
                                software's dictionary.
      --density DENSITY         Specify the density (in g/cm3) of the
                                material.
      --dictionary              Print the software's dictionary of
                                pre-defined materials.
      --atomicid Z,A            Specify atomic number and atomic mass of
                                mono-nucleus material.
      --help                    Print this help.
      
```
