.. _quickstart:


Installation
============
First, make sure you have the third party dependencies:

* `samtools <https://www.htslib.org/>`_
* `minimap2 <https://github.com/lh3/minimap2>`_
* `bowtie2 <https://github.com/BenLangmead/bowtie2>`_

You can then install hicstuff latest stable version using pip:

.. code-block:: bash

    pip3 install --user hicstuff

Or the latest development (unstable) version using:

.. code-block:: bash

    pip3 install -e git+https://github.com/koszullab/hicstuff.git@master#egg=hicstuff


Quickstart
==========

The fastest way to generate Hi-C matrices is to use the `hicstuff pipeline` command:

.. code-block:: bash

    hicstuff pipeline -g bt2_index for.fq rev.fq

However, you most likely want to have a look at `the command line help <https://hicstuff.readthedocs.io/en/latest/api/hicstuff.html#hicstuff.commands.Pipeline>`_ to select appropriate options, such as the enzyme used in the experiment. The help can be displayed using:

.. code-block:: bash

    hicstuff pipeline --help

Matrices generated in the default coordinate format can then be visualised using the `view <https://hicstuff.readthedocs.io/en/latest/api/hicstuff.html#hicstuff.commands.View>`_ command, modified using the `rebin <https://hicstuff.readthedocs.io/en/latest/api/hicstuff.html#hicstuff.commands.Rebin>`_ and `convert <https://hicstuff.readthedocs.io/en/latest/api/hicstuff.html#hicstuff.commands.Convert>`_ commands, or used as input for other softwares.
