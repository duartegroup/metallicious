Installation
=====
.. _installation:


The easiest way to install *metallicious* is to use `Anaconda <https://anaconda.org/anaconda/python>`_, as *metallicious* depends on external packages.
Create a new conda environment:

.. code-block:: bash

    conda create --name metallicious
    conda activate metallicious

Install dependencies and the *metallicious*:

.. code-block:: bash

    conda install rdkit parmed autode psiresp mdanalysis networkx qcelemental==0.25.1 ambertools --channel conda-forge
    pip install metallicious

Alternatively, if you do not need to parametrize templates, install core dependencies + antechamber :

.. code-block:: bash

    conda install rdkit mdanalysis networkx ambertools --channel conda-forge
    pip install metallicious


Dependencies
----------------

The core dependencies of *metallicious*:

* `rdkit <https://www.rdkit.org/>`_
* `networkx <https://networkx.org/>`_
* `MDAnalysis <https://www.mdanalysis.org/>`_
* `ParmEd <https://parmed.github.io/ParmEd/html/index.html>`_

(**Optional**) the parametrization of templates requires:

* `autode <https://github.com/duartegroup/autodE>`_
* `ORCA <https://orcaforum.kofo.mpg.de/app.php/portal>`_
* `psiRESP <https://github.com/lilyminium/psiresp>`_

(**Optional**) simple force-field parametrization with General Amber Force-field (GAFF) requires:

* `ambertools <https://ambermd.org/AmberTools.php>`_
