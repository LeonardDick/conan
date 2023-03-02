Installation
============

To use CONAn, simply clone the public repository to your machine using using GitHub.

    $ git clone ...

The code is written in Python. Multiple libraries need to be installed to run the code, which are listed in the requirements.txt file. 
For the installation using pip, run the following command:

    $ pip install -r requirements.txt

Or manually install all packages using the following command:

    $ pip install pandas numpy scipy prettytable matplotlib

For the installation in a new conda environment, run the following commands:

    $ conda create -n conan python=3.8
    $ conda config --env --add channels conda-forge
    $ conda install -n conan --file requirements.txt 
    $ conda activate conan

Now the code is ready to be used. To start the program, simply run the following command:
    
    $ python3.8 CONAn.py
    
The code works by asking questions to the user and prints the results to the terminal. 
A log file will be written called ``conan.log``, which contains everything that is printed to the terminal. 
For certain modules, additional files will be created (e.g. a .csv file for the results of the analysis).

