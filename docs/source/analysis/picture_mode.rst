Picture mode
============
The picture mode is able to generate seperate xyz files, pulled from the trajectory. The xyz files can be used to generate pictures with the help of other programs (e.g. VMD).


Simulation Box Picture
----------------------
Save the xyz structure of the whole simulation box.
The frame of choice is set beforehand by the user when prompted to select a frame with which the carbon structures shall be identified.

Pore Picture
-----------
Save the structure of the pore structure in the simulation box in a seperate xyz file.
Additionally, the center point of the CNT can be included in to the xyz file, labeled as 'X'.

CNT Picture
----------
Save the xyz structure of the CNT (without the outer walls).
Additionally, it is possible to add all atoms to the xyz file, which are contained inside the CNT. 
As some molecules might be partially inside the CNT and partially outside, some molecules might be fragmented.
The other option is to add all molecules to the xyz file, which are completely, or just partially, contained inside the CNT.



