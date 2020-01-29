#include <cxxopts.hpp>
#include <fstream>
#include <iostream>
#include <map>
#include <utility>
#include <cmath>

const double gFineStructureConstant = 0.00729735256;
const double gRadiationLengthConstant = 0.001395852;

typedef unsigned int AtomicNumber_t;
typedef double AtomicMass_t;
typedef double MassFraction_t;
typedef std::map<std::pair<AtomicNumber_t, AtomicMass_t>, MassFraction_t>
    CompositionMap_t;

double L_rad(AtomicNumber_t vZ) {
  double l_rad = 0;

  if (vZ <= 0) {
    std::cout << "Invalid atomic number. Exiting..." << std::endl;
    return 0;
  }

  switch (vZ) {
  case 1: {
    l_rad = 5.31;
    break;
  }
  case 2: {
    l_rad = 4.79;
    break;
  }
  case 3: {
    l_rad = 4.74;
    break;
  }
  case 4: {
    l_rad = 4.71;
    break;
  }
  default: {
    l_rad = log(184.15 * pow(vZ, -1. / 3.));
  }
  }

  return l_rad;
}

double L_rad_primed(AtomicNumber_t vZ) {
  double l_rad_primed = 0;

  if (vZ <= 0) {
    std::cout << "Invalid atomic number. Exiting..." << std::endl;
    return 0;
  }

  switch (vZ) {
  case 1: {
    l_rad_primed = 6.144;
    break;
  }
  case 2: {
    l_rad_primed = 5.621;
    break;
  }
  case 3: {
    l_rad_primed = 5.805;
    break;
  }
  case 4: {
    l_rad_primed = 5.924;
    break;
  }
  default: {
    l_rad_primed = log(1194 * pow(vZ, -2. / 3.));
  }
  }

  return l_rad_primed;
}

double F_z(AtomicNumber_t vZ) {
  double a = gFineStructureConstant * (double)vZ;
  double a_sq = a * a;
  return a_sq * (1. / (1 + a_sq) + 0.20206 - 0.0369 * a_sq +
                 0.0083 * pow(a, 4) - 0.002 * pow(a, 6));
}

double Calc_radiation_length_mononucleus_material(AtomicNumber_t vZ,
                                                  AtomicMass_t vA) {
  double invertedRadLength =
      gRadiationLengthConstant * ((double)(vZ * vZ) * (L_rad(vZ) - F_z(vZ)) +
                                  (double)vZ * L_rad_primed(vZ));
  return vA / invertedRadLength;
}

double Calc_radiation_length_composite_material(
    CompositionMap_t vMaterialCompositionMap) {
  double invertedRadLength = 0.;
  for (auto material : vMaterialCompositionMap) {
    double matRadLength = Calc_radiation_length_mononucleus_material(
        material.first.first, material.first.second);
    invertedRadLength += material.second / matRadLength;
  }

  return 1. / invertedRadLength;
}

double gMaterialDensity = 0.;
std::string gMaterialCompositionDatFile;
std::pair<AtomicNumber_t, AtomicMass_t> gMaterialAtomicID;
std::string gPredefinedMaterialName;

enum EMaterialClass {
  kMonoNucleus,
  kPredefined,
  kComposition,
  kPrintMaterialDictionary,
  kMaterialUndefinedError
} gMatericalClass;

void ArgumentParsing(int argc, char *argv[]) {
  try {
    cxxopts::Options options(argv[0],
                             "A radiation length calculator for materials.\n");
    options.custom_help(
        "[OPTION...]\n"
        "\n\tIf --file option is chosen, the program will read the input data "
        "file and calculate the radiation length"
        "\n\tof the composite material described by the data file. The data "
        "file is an ASCII file, formatted in three"
        "\n\tcolumns respectively containing the atomic number (Z), the atomic "
        "mass (A) and the mass fraction (%). "
        "\n\tAn example of the file can be as following:"
        "\n"
        "\n\t\t1\t 1.008\t\t30"
        "\n\t\t6\t12.012\t\t60"
        "\n\t\t8\t16.002\t\t10"
        "\n"
        "\n\tThis file describes a material composed of H, C and O with mass "
        "fraction of 30%, 60% and 10%, respectively."
        "\n\n\n"
        "OPTIONS:");

    bool apple = false;

    options.add_options()("file",
                          "Specify material composition data file, --material "
                          "and --atomicid options will be ignored.",
                          cxxopts::value<std::string>(), "FILE")(
        "material",
        "Specify a pre-defined material in the software's dictionary.",
        cxxopts::value<std::string>(), "MATERIAL_NAME")(
        "density", "Specify the density (in g/cm3) of the material.",
        cxxopts::value<double>(),
        "DENSITY")("dictionary",
                   "Print the software's dictionary of pre-defined materials.")(
        "atomicid",
        "Specify atomic number and atomic mass of mono-nucleus material.",
        cxxopts::value<std::vector<double>>(),
        "Z,A")("help", "Print this help.")
#ifdef CXXOPTS_USE_UNICODE
        ("unicode", u8"A help option with non-ascii: à. Here the size of the"
                    " string should be correct")
#endif
        ;

    auto result = options.parse(argc, argv);

    if (result.count("help")) {
      std::cout << options.help({""}) << std::endl;
      exit(0);
    }

    if (result.count("density")) {
      gMaterialDensity = result["density"].as<double>();
    }

    if (result.count("file")) {
      gMaterialCompositionDatFile = result["file"].as<std::string>();
      gMatericalClass = kComposition;
    } else {
      if (result.count("atomicid")) {
        auto grepZA = result["atomicid"].as<std::vector<double>>();
        double tmp_z = grepZA.at(0);
        double tmp_a = grepZA.at(1);
        gMaterialAtomicID = std::make_pair((AtomicNumber_t)tmp_z, tmp_a);
        gMatericalClass = kMonoNucleus;
      } else {
        if (result.count("material")) {
          gPredefinedMaterialName = result["material"].as<std::string>();
          gMatericalClass = kPredefined;
        } else {
          if (result.count("dictionary")) {
            gMatericalClass = kPrintMaterialDictionary;
          } else {
            gMatericalClass = kMaterialUndefinedError;
          }
        }
      }
    }

  } catch (const cxxopts::OptionException &e) {
    std::cout << "Error parsing options: " << e.what() << std::endl;
    exit(1);
  }
}

void PrintUndefinedMaterialError() {
  std::cout << std::endl;
  std::cout << "ERROR: " << std::endl;
  std::cout << "\t- Use one of the options -f, -m and -z to specify a material."
            << std::endl;
  std::cout << "\t- Or use option -d to print a list of pre-defined materials."
            << std::endl;
  std::cout << std::endl;
  std::cout << "Exiting..." << std::endl;
}

void PrintMaterialDictionary() {
  std::cout << std::endl;
  std::cout << "\t*** Printing Material Dictionary *** \t" << std::endl;
}

void PrintMononucleusRadiationLength() {
  std::cout << std::endl;
  double radLength = Calc_radiation_length_mononucleus_material(
      gMaterialAtomicID.first, gMaterialAtomicID.second);
  std::cout << "Radiation length of material with Z = " << gMaterialAtomicID.first << ", A = " << gMaterialAtomicID.second << " is " << radLength << " g/cm2. ";
  if (gMaterialDensity != 0.) {
    std::cout << "Value corrected for a density of " << gMaterialDensity << " g/cm3 is " << radLength / gMaterialDensity << " cm.";
  }
  std::cout << std::endl;
}

CompositionMap_t ReadCompositionDatFile(std::string vCompositionDatFile) {
  std::ifstream datfile(vCompositionDatFile.c_str());

  CompositionMap_t compositionMap;
  double tmp_z = 0.;
  double tmp_a = 0.;
  double tmp_f = 0.;
  while (datfile >> tmp_z >> tmp_a >> tmp_f) {
    compositionMap.insert(
        std::make_pair(std::make_pair(tmp_z, tmp_a), tmp_f / 100.));
  }

  datfile.close();

  return compositionMap;
}

void PrintCompositionRadiationLength() {
  std::cout << std::endl;
  std::cout << "Input material composition data file: "
            << gMaterialCompositionDatFile << "." << std::endl;

  CompositionMap_t composition =
      ReadCompositionDatFile(gMaterialCompositionDatFile);
  double radLength = Calc_radiation_length_composite_material(composition);

  std::cout << "Radiation length of \"" << gMaterialCompositionDatFile.c_str() << "\" is " << radLength << " g/cm2. ";
  if (gMaterialDensity != 0.) {
    std::cout << "Value corrected for a density of " << gMaterialDensity << " g/cm3 is " << radLength / gMaterialDensity << " cm.";
  }
  std::cout << std::endl;
  std::cout << std::endl;
}

int main(int argc, char *argv[]) {
  ArgumentParsing(argc, argv);

  switch (gMatericalClass) {
    case kMaterialUndefinedError: {
        PrintUndefinedMaterialError();
        break;
    }
    case kPrintMaterialDictionary: {
        PrintMaterialDictionary();
        break;
    }
    case kMonoNucleus: {
        PrintMononucleusRadiationLength();
        break;
    }
    case kComposition: {
        PrintCompositionRadiationLength();
        break;
    }
    default: {
        break;
    }
  }

  return 0;
}
