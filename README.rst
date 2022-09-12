=======================
Scipion GROMACS plugin
=======================

In order to use this plug-in, you need to have Scipion3 installed
(https://scipion-em.github.io/docs/docs/scipion-modes/how-to-install.html).

GROMACS version 21.5 is automatically installed.
If you wanted to use any other preinstalled version, you would need to modify the scipion.conf file adding:

GROMACS_HOME=<path_to_gromacs_home_folder>

To install the plugin,  you have to follow the following steps:

**Clone this repository:**

.. code-block::

    git clone https://github.com/scipion-chem/scipion-chem-gromacs.git


**Install the plugin in devel mode**

.. code-block::

    scipion3 installp -p path_to_scipion-chem-gromacs --devel -j <numberOfProcessors>
    
OR
    
**Install the plugin in user mode**

.. code-block::

    scipion3 installp -p path_to_scipion-chem-gromacs
