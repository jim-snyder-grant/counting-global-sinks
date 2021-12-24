# counting-global-sinks

Code to generate and count digraphs that have a single global sink. 

Built on the nauty library of B D McKay.

## Download and build

The makefile relies on a copy of the nauty library files in a parallel nauty directory.
It also needs data files as input, available from B D McKays data repository. Some of these are very large,
so instad of putting them in source control, there is a shell file, getdata.sh, to download them. 

From an empty working local directory:
```
git clone https://github.com/jim-snyder-grant/counting-global-sinks.git
git clone https://github.com/jim-snyder-grant/nauty.git
cd nauty
./configure
make
cd ../counting-global-sinks
make
./getdata.sh
./gsinks   
```
(./gsinks -help for options)
