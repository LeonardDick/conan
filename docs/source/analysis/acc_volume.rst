Accessible Volume
=================
The accessible volume is calculated by scanning over the whole trajectory to identify the atom which is furthest displaced from the center line of the CNT, while including the van der Waals radius of the given element.
The accessible volume :math:`V_{acc}` is then calculated by simplifying the CNT as a cylinder, using the following equation:

.. math::

    V_{acc} = \pi*r_{acc}^2*h_{CNT}    

where :math:`r_{acc}` is the radius of the atom which is furthest displaced from the center line of the CNT. :math:`h_{CNT}` is the length of a given CNT.

.. note::

    The accessible volume is therefore dependent on the most displaced atom and not directly derived by the CNTs diameter.
    Results can thus slightly differ for varying simulations with the same CNT, and a large sampling is recommended to yield the best results.
    
The calculated volumes are the groundwork for further analysis regarding e.g. density profiles along the z axis.

