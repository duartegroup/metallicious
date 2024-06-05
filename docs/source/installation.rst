Usage
=====

.. _installation:

Installation
------------

The easiest way to install _metallicious_ is to use conda enviroment and pip:

conda create --name metallicious
conda activate metallicious
conda install rdkit autode psiresp mdanalysis networkx qcelemental==0.25.1 ambertools --channel conda-forge
pip install metallicious

.. code-block:: console

   (.venv) $ pip install lumache

Dependencies
----------------

The core dependencies of _metallicious_:
- [rdkit](https://www.rdkit.org/)
- [networkx](https://networkx.org/)
- [MDAnalysis](https://www.mdanalysis.org/)
- [ParmEd](https://parmed.github.io/ParmEd/html/index.html)

(**Optional**) the parametrization of templates requires:
- [autode](https://github.com/duartegroup/autodE)
- [ORCA](https://orcaforum.kofo.mpg.de/app.php/portal)
- [psiRESP](https://github.com/lilyminium/psiresp)

(**Optional**) simple force-field parametrization with General Amber Force-field (GAFF) requires:
- [ambertools](https://ambermd.org/AmberTools.php)




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

