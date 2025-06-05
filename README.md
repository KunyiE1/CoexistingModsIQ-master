# Coexisting Modifications Combinations Identification and Quantification
This tool is implemented by changing original topmg which is included in TopPIC. 
For further information and detailed manual, please visit https://www.toppic.org/software/toppic/

## System requirements

* GCC version higher than 5.5.0 for C++14 support
* CMake (>= 3.1)

### Linux (Ubuntu 22.04):

```sh
# install compiling tools
sudo apt install build-essential cmake

# install other dependencies
sudo apt install zlib1g-dev 
sudo apt install libboost-filesystem-dev 
sudo apt install libboost-program-options-dev 
sudo apt install libboost-system-dev 
sudo apt install libboost-thread-dev 
sudo apt install libboost-iostreams-dev 
sudo apt install libboost-chrono-dev 

sudo apt install libxerces-c-dev  
sudo apt install libeigen3-dev 
sudo apt install nlohmann-json3-dev


# install the catch unit test framework (https://github.com/philsquared/Catch)
sudo apt install catch

# Qt5 for GUI
sudo apt install qtbase5-dev

# building
mkdir build
cd build
cmake ..
make -j$(nproc)
make install
```

### Linux (CentOS Stream 8):

```sh
# install Extra Packages for Enterprise Linux (EPEL)
sudo dnf install 'dnf-command(config-manager)'
sudo dnf config-manager --set-enabled powertools
sudo dnf install epel-release 

# install compiling tools
sudo dnf install gcc gcc-c++ make cmake

# install dependencies
sudo dnf install zlib-devel
sudo dnf install boost-devel 
sudo dnf install xerces-c-devel
sudo dnf install eigen3-devel
sudo dnf install json-devel

# Qt5 for GUI
sudo dnf install qt5-qtbase-devel

# building
mkdir build
cd build
cmake ..
make -j$(nproc)
make install
```

The installation time is about 1 minutes. 

## Simple Manual for Getting Started
The input of the program includes:
- **var_mods.txt**: a text file with variable modification informations
- **Spectrum data**: .msalign file ans .feature file
- **ref_peptide.txt**: containing the sequences corresponding to the spectra, should be in the same directory as spectrum data files.

### Run on a toy dataset
After TopPIC is successfully built(or you can only build topmg by run topmg in IDE like Clion)
You can simply run the program on a toy(demo) dataset with command line:
```sh
topmg -i ToyCaseForTesting/var_mods.txt ToyCaseForTesting/sim_ms2.msalign
```
according to the building, the command line could also be:
```sh
bin/topmg -i ToyCaseForTesting/var_mods.txt ToyCaseForTesting/sim_ms2.msalign
```
The toy dataset contains 10 spectra and their corresponding peptide sequence. The expected running time is about 2 seconds. 

The results will be stored in the results.txt as an ouput file.

## Simulator
The python3.9 code for simulator could be found in the folder **"HomMTM PSM simulator"**

packages required:
- numpy
- random
- re
- pickle
- math

just run the python script **RunSimulator.py**, simulated HomMTM PSMs will be generated. 

### Change Settings for Simulator
In **RunSimulator.py**, line 8, **simulation_num = 1000** means 1000 simulated PSMs will be generated at one time.
output directory is set in line 31.

In **Simulator.py**, line 196, **abund_1, abund_2 = 0.2, 0.8** means the abundances of two isoforms in the simulated PSM will be 0.2 and 0.8.

