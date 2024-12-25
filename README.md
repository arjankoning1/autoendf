
# AUTOENDF
AUTOENDF is a collection of bash shell scripts for the automatic testing and processing of ENDF-6 files 
with the checking and processing codes CHECKR, FIZCON, PSYCHE, INTER, PREPRO and NJOY. 
The scripts in AUTOENDF will produce input files, run the codes and automatically analyze the output for errors and warnings. 
This diagnosis is obtained by grepping the error and warning messages as printed by the various codes.

## Documentation and reference
The user manual for AUTOENDF can be found here: [AUTOENDF User Manual (pdf)](https://github.com/arjankoning1/autoendf/blob/main/doc/tools.pdf).
The reference to be used for AUTOENDF is

A.J. Koning, D. Rochman, J.-Ch. Sublet, N. Dzysiuk, M. Fleming, and S. van der Marck, *TENDL: Complete Nuclear Data Library for innovative Nuclear Science and Technology*, Nuclear Data Sheets 155,1 (2019).

## Installation

### Prerequisites:

The following are the prerequisites for using AUTOENDF:
  - git (if the package is downloaded via Github)
  - bash shell environment
  - Installed Fortran codes: CHECKR, FIZCON, PSYCHE, INTER (from BNL checking code suite), PREPRO and NJOY. The executables are expected to have the following location and names:
    - \<your bin directory\>/checkr
    - \<your bin directory\>/fizcon
    - \<your bin directory\>/psyche
    - \<your bin directory\>/inter
    - \<your bin directory\>/recent  (and similar for the other PREPRO codes sigma1, groupie etc.)
    - \<your bin directory\>/xnjoy

For this you will need to install
  - [PREPRO](https://github.com/IAEA-NDS/PREPRO)
  - [BNL ENDF-utility-codes](https://github.com/IAEA-NDS/ENDF-utility-codes)
  - [NJOY](https://github.com/njoy)

### Instructions:

To install AUTOENDF:
```
  - git clone https://github.com/arjankoning1/autoendf.git
  - cd autoendf/bin
  - Edit autobnl, autoprepro, autonjoy and change the Thome and bin variables to match with <your bin directory>
```

## Sample case

A successful installation can be verified by running the sample case. 
```
cd autoendf/sample
./verify
```

For simple ENDF-6 format validation, the most important output is probably given in the *checkr.ers*, *njoy.ers* etc. files, 
which contain a diagnosis of the (usually long) output files of the various checking and processing codes. 
When the *.ers* files are not empty, the associated output files *checkr.out*, *njoy.out*, 
etc. can be consulted for more information on the errors or warnings.

## The AUTOENDF package

The *autoendf/* directory contains the following directories and files:

+ `README.md` is this README file.
+ `LICENSE` is the License file.
+ `bin/` contains the scripts *autobnl*, *autoprepro*, *autonjoy* and, not currently active, *autofudge* and *automcnp*.
+ `doc/` contains the tutorial in pdf format (in the appendices).
+ `sample/` contains the ENDF-6 formatted file *n-Nb093.tendl.2021* to be used as test case.

In total, you will need about 500 Mb of free disk space to run the sample case of AUTOENDF.

## License and Copyright
This software is distributed and copyrighted according to the [LICENSE](LICENSE) file.
