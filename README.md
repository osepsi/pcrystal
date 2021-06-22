# Introduction

This package is used to calculate the crystallization of three-dimensional radially growing crystal seeds using a probabilistic method given an initial set of seeds and birth times. As a result the crystallization curve as well as the size distribution is obtained. A constant growth speed is assumed. The core routine is written in Fortran for fast calculation and it is extended in python for easy usage. Please see below the compilation process of the fortran module. Apart of the Fortran module only scientific python packages (numpy, scipy, matplotlib) is needed for efficient usage.


# Compiling gcrystal
## On linux

0. Install gfortran
	`sudo apt install gfortran`

1. Install anaconda
	a. get anaconda `curl -O https://repo.anaconda.com/archive/Anaconda3-2019.07-Linux-x86_64.sh`
	b. install anaconda `bash Anaconda3-2019.07-Linux-x86_64.sh`
	c. Answer "yes" to the question "Do you wish the installer to initialize Anaconda3 by running conda init?"                  

2. create a new environment
	a. conda create --name myenv
	b. conda activate myenv
	c. conda install numpy scipy matplotlib jupyter

3. Now you should have f2py. Check it!
	`which f2py` should give something like user_home_folder/anaconda3/envs/myenv/bin/f2py

4. Compile fcrystal.f90
	a. Go to the directory where fcrystal.f90 resides.
	b. `f2py -c -m fcrystal fcrystal.f90`
After executing this, I have a fresh "fcrystal.cpython-37m-x86_64-linux-gnu.so" file. 
Yours might have slightly different name, but it shoud start with fcrystal and should
be created freshly. Check the creation date with `ls -lh`

5. Test it
	a. open jupyter from the same terminal `jupyter notebook`
 	   In case no broswer opens, open it manually by copying the given web addrress.
	b. open the test.ipynb file from jupyter
	c. Run the first cell.

For later use always open terminal and execute `conda antivate myenv`, then `jupyter notebook`.


---

## On windows

It is strongly advised to use Windows Subsystem for Linux and continue with the Linux version. Otherwise:

1. Install anaconda, and necessary packages (numpy, scipy, matplotlib, jupyter, ...)

2. Install MinGW
	- MinGW-64 if anaconda is 64 bit (https://sourceforge.net/projects/mingw-w64/files/latest/download)
		- Architecture: x86_64
		- Other options remain the default
	- or MinGW32 if anaconda installation is 32 bit (https://sourceforge.net/projects/mingw/files/MinGW/Base/gcc/Version6/)

3. Setup environment variables
	- create "C_INCLUDE_PATH" = path_to_mingw_install\include
		- For example on my computer:
		  C:\mingw-w64\x86_64-8.1.0-posix-seh-rt_v6-rev0\mingw64\include
	- append to "Path" path_to_mingw_install\bin
		- For example on my computer:
		  C:\mingw-w64\x86_64-8.1.0-posix-seh-rt_v6-rev0\mingw64\bin

4. Open Anaconda terminal
	python anaconda_path\Scripts\f2py.py -c --fcompiler=gnu95 --compiler=mingw32 -m fcrystal fcrystal.f90
	On my computer: python c:\Users\Lenovo\Anaconda3\Scripts\f2py.py -c --fcompiler=gnu95 --compiler=mingw32 -m fcrystal fcrystal.f90
