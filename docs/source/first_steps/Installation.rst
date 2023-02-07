Installation
============

To use CONAn, simply clone the public repository to your machine using using GitHub.
.. code-block::

    $ git clone ...

The code is written in Python 3.8.10., however multiple libraries are relied on to run the code. These are listed in the requirements.txt file. To install these, run the following command::
.. code-block::

    $ pip install -r requirements.txt

You can also manually install all packages with the following command:
.. code-block::

    $ pip install argparse gettext os xml.dom pandas numpy math time scipy gc matplotlib.pyplot prettytable cdist

Now the code is ready to run. To run the code, simply run
.. code-block::
    
    $ python3 CONAn.py
    
The code works by asking questions to the user and prints the results to the terminal. 
A log file will be written called ``conan.log``, which contains everything, that is printed to the terminal. 
For certain modules, additional files will be created (e.g. a .csv file for the results of the analysis).

