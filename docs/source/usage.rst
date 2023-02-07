Usage
=====

.. _installation:

Installation
------------

To use CONAn, simply clone the public repository to your machine using using GitHub::

    $ git clone (EDIT HERE)

The code is written in Python 3.8.10., however multiple libraries are relied on to run the code. These are listed in the requirements.txt file. To install these, run the following command::

    $ pip install -r requirements.txt

You can also manually install all packages with the following command::

    $ pip install argparse gettext os xml.dom pandas numpy math time scipy gc matplotlib.pyplot prettytable cdist


Creating recipes
----------------

To retrieve a list of random ingredients,
you can use the ``lumache.get_random_ingredients()`` function:

.. autofunction:: lumache.get_random_ingredients

The ``kind`` parameter should be either ``"meat"``, ``"fish"``,
or ``"veggies"``. Otherwise, :py:func:`lumache.get_random_ingredients`
will raise an exception.

.. autoexception:: lumache.InvalidKindError

For example:

>>> import lumache
>>> lumache.get_random_ingredients()
['shells', 'gorgonzola', 'parsley']

