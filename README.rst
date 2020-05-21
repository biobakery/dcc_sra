##########
DCC -> SRA
##########

This module contains anadama workflows and pipelines for gathering
sequences and metadata from the iHMP OSDF instance and submitting it
to SRA.

.. contents::

________________________________________________

Installation
============

You'll need the following dependencies installed before getting started:

- ``git``  
- Python 2.7
- `setuptools <https://pypi.python.org/pypi/setuptools>`_
- `pip <https://pypi.python.org/pypi/pip/>`_

First, install OSDF-python and cutlass::

  pip install 'git+https://github.com/ihmpdcc/osdf-python.git@master#egg=osdf-python-0.3.1'
  pip install 'git+https://github.com/ihmpdcc/cutlass.git@master#egg=osdf-cutlass-0.8.1'

Then, install AnADAMA::

   pip install 'git+https://bitbucket.org/biobakery/anadama.git@master#egg=anadama-0.0.1'

Finally, install the dcc_sra module::

  pip install 'git+https://github.com/biobakery/dcc_sra.git@master#egg=dcc_sra-0.0.1'


  
Using the pipeline
==================

There are many ways to do the same thing

Use the ``pipeline`` command::

  anadama pipeline dcc_sra


Tired of the pipeline asking you questions? Use options::

  anadama pipeline dcc_sra -o 'serialize.study_id: 700f3b77246a81f33c15af8787319115'

Don't want to write big long shell commands? Use the skeleton::

  anadama skeleton dcc_sra
  # options are now editable as YAML config files
  emacs input/_options/serialize.txt

  anadama run
  

Got lost? Read the help::

  anadama help pipeline dcc_sra
  anadama help pipeline
  anadama help

To read a webpage about AnADAMA, go to the `docs <http://huttenhower.sph.harvard.edu/docs/anadama/index.html>`_.

