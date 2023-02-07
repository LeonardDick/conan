Molecule Identifier
===================
A molecule identifier is included in the program. 
All bonds between atoms are identified via distance criteria. 
The cutoff distances, below which two atoms are considered bonded, depend on the elements. 
All combinations of elements are listed in a library.
The program deduces all molecules from the chosen frame in the trajectory and therefore cannot handle bond breaking or formation in the trajectory.
In case the trajectory is in pdb format, the program directly identiies molecules from the molecule information given in the pdb file. 
The module is therefore much faster using the pdb format.
With every molecule identified, the program printst the number of unique molecules and how many of each are in the system.

.. note::
    The molecule identifier is not yet able to hande PBC. 
    If the molecules have to be identified in a periodic system, the molecules have to be wrapped before.
    The library of cutoff distances is not yet complete, and cutoff distances with elements apart from C, H, N, O have yet to be implemented.



