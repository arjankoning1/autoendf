# AUTOENDF

AUTOENDF is a collection of Bash scripts for automatically checking and processing ENDF-6 formatted nuclear data files with several widely used nuclear-data tools. It can run the BNL checking codes, PREPRO, NJOY and FUDGE, generate the required input files, and extract errors and warnings from their output into compact diagnostic files.

The main scripts are:

- `autobnl` for CHECKR, FIZCON, PSYCHE and INTER
- `autoprepro` for PREPRO
- `autonjoy` for NJOY
- `autofudge` for FUDGE

Each script prints its available options when invoked without arguments.

## Documentation and reference

A description of AUTOENDF and its options can be found in the [AUTOENDF tutorial (pdf)](https://github.com/arjankoning1/autoendf/blob/main/doc/tools.pdf).

The reference to be used for AUTOENDF is:

A.J. Koning, D. Rochman, J.-Ch. Sublet, N. Dzysiuk, M. Fleming, and S. van der Marck, *TENDL: Complete Nuclear Data Library for innovative Nuclear Science and Technology*, Nuclear Data Sheets 155, 1 (2019).

## Installation

### Prerequisites

AUTOENDF itself consists of Bash scripts and does not require compilation.

The following are required for full use of AUTOENDF:

- a Bash shell environment and standard Unix command-line utilities
- the BNL ENDF checking codes CHECKR, FIZCON, PSYCHE and INTER
- PREPRO
- NJOY
- FUDGE
- git, only when AUTOENDF is downloaded using `git clone`

The external packages are available from:

- [PREPRO](https://github.com/IAEA-NDS/PREPRO)
- [BNL ENDF utility codes](https://github.com/IAEA-NDS/ENDF-utility-codes)
- [NJOY](https://github.com/njoy)
- [FUDGE](https://github.com/LLNL/fudge)

### Downloads

AUTOENDF can be downloaded in one of the following ways.

#### 1. Latest version without git

Users who do not have git can download a snapshot of the current `main` branch directly from GitHub:

```bash
curl -L \
  -o autoendf-main.tar.gz \
  https://github.com/arjankoning1/autoendf/archive/refs/heads/main.tar.gz

tar zxf autoendf-main.tar.gz
mv autoendf-main autoendf
```

This produces the same `autoendf/` directory structure as the git version, but without the git history.

The downloaded snapshot contains the latest version of the `main` branch at the time of download. To obtain a newer version later, download the snapshot again.

#### 2. Latest version using git

Users with git can clone the repository with

```bash
git clone https://github.com/arjankoning1/autoendf.git
```

The advantage of this method is that the local AUTOENDF installation can subsequently be updated with

```bash
cd autoendf
git pull --ff-only
```

### Runtime configuration

The AUTOENDF scripts are located in `autoendf/bin/`. To run them from anywhere, add this directory to `PATH`, for example:

```bash
export PATH="/path/to/autoendf/bin:$PATH"
```

The external executables can be selected in two ways.

#### 1. Specify the executable directory on the command line

The scripts accept a `-bin` option. For example:

```bash
autobnl    -file myfile.endf -bin /path/to/bnl/bin/
autoprepro -file myfile.endf -bin /path/to/prepro/bin/
autonjoy   -file myfile.endf -bin /path/to/njoy/bin/
autofudge  -file myfile.endf -bin /path/to/fudge/.venv/bin/
```

For NJOY, the default executable name is `xnjoy`; another name can be selected with `-version`.

#### 2. Use AUTOENDF_HOME

If `-bin` is not specified, the scripts derive their default external-code locations from `AUTOENDF_HOME`. If `AUTOENDF_HOME` is not defined, `$HOME` is used.

The default layout is:

```text
$AUTOENDF_HOME/
├── bin/
│   ├── checkr
│   ├── fizcon
│   ├── psyche
│   ├── inter
│   ├── linear
│   ├── recent
│   ├── sigma1
│   ├── ...
│   └── xnjoy
└── tools/
    └── fudge/
        └── .venv/
            └── bin/
                ├── endf2gnds.py
                ├── gnds2endf.py
                └── checkGNDS.py
```

For such an installation, set for example:

```bash
export AUTOENDF_HOME="/path/to/your/nuclear-data-tools"
```

This line and the AUTOENDF `PATH` setting can be added to `~/.zshrc` or `~/.profile`.

It is not necessary to edit `Thome` or `bin` inside the AUTOENDF scripts.

## Running AUTOENDF

Each wrapper requires an ENDF-6 file through the `-file` option. Typical examples are:

```bash
autobnl    -file myfile.endf
autoprepro -file myfile.endf
autonjoy   -file myfile.endf
autofudge  -file myfile.endf
```

The scripts generate the required input files, run the corresponding external programs and, by default, perform a simple diagnosis of their output.

For ENDF-6 format validation, particularly useful files are the diagnostic `*.ers` files such as:

```text
checkr.ers
fizcon.ers
psyche.ers
inter.ers
njoy.ers
fudge.ers
```

A non-empty diagnostic file indicates that AUTOENDF found messages matching its error or warning criteria. The corresponding full output file, such as `checkr.out`, `njoy.out` or `fudge.out`, should then be inspected for details.

## Sample case

A successful installation can be tested with the supplied Nb-93 TENDL-2025 evaluation:

```bash
cd autoendf/samples
./verify
```

The verification script runs:

```text
autobnl
autoprepro
autonjoy
autofudge
```

on:

```text
n-Nb093.tendl.2025
```

The required external codes must be available through the default `AUTOENDF_HOME` layout described above for this unmodified verification script.

## The AUTOENDF package

The `autoendf/` directory contains:

- `README.md` this README file
- `LICENSE` the license file
- `bin/` the AUTOENDF scripts `autobnl`, `autoprepro`, `autonjoy`, `autofudge` and `automcnp`
- `doc/` the AUTOENDF documentation
- `samples/` the Nb-93 ENDF-6 test file and the `verify` script

The `automcnp` script is included in the repository but is not part of the standard `samples/verify` sequence.

Approximately 500 MB of free working space is recommended when running the full sample processing chain.

## License and Copyright

This software is distributed and copyrighted according to the [LICENSE](LICENSE) file.
