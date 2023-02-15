=======================
Scipion GROMACS plugin
=======================

**Documentation under development, sorry for the inconvenience**

In order to use this plug-in, you need to have Scipion3 installed
(https://scipion-em.github.io/docs/docs/scipion-modes/how-to-install.html).

GROMACS version 2021.5 is automatically installed.
If you wanted to use any other preinstalled version, you would need to modify the scipion.conf file adding:

GROMACS_HOME=<path_to_gromacs_home_folder>

To install the plugin,  you have to follow the following steps:

1. **Clone this repository:**

.. code-block::

    git clone https://github.com/scipion-chem/scipion-chem-gromacs.git


2. **Switch to the desired branch** (master or devel):

Scipion-chem-gromacs is constantly under development and including new features.
If you want a relatively older an more stable version, use master branch (default).
If you want the latest changes and developments, user devel branch.

.. code-block::

            cd scipion-chem-gromacs
            git checkout devel

3. **Install**:

.. code-block::

    scipion3 installp -p path_to_scipion-chem-gromacs --devel -j <numberOfProcessors>
    
OR
    
**Install the plugin in user mode** (not available yet)

.. code-block::

    scipion3 installp -p path_to_scipion-chem-gromacs
