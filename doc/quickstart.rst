.. _quickstart:


Installation
============
hicstuff needs Python 3.10 or later, on Linux or macOS. On Windows, use WSL2:
pysam, a dependency, publishes no Windows wheels.

.. code-block:: bash

    pip install hicstuff

That alone runs the whole pipeline, from FASTQ files to matrices, with bwa:
``hicstuff pipeline -a bwa -e DpnII -g genome.fa R1.fq.gz R2.fq.gz`` needs no
other program. bwa then runs in-process through
`pybwa <https://github.com/fulcrumgenomics/pybwa>`_, which pip installs with
hicstuff.

Or the latest development (unstable) version:

.. code-block:: bash

    pip install git+https://github.com/koszullab/hicstuff.git@master

Aligners and other programs
---------------------------
``--aligner`` chooses the aligner that maps the reads:

* `bowtie2 <https://github.com/BenLangmead/bowtie2>`_, the default: the
  ``bowtie2`` program, which is not on PyPI (install it with apt, conda or
  brew);
* `bwa <https://github.com/lh3/bwa>`_: the ``bwa`` program if it is on PATH,
  otherwise in-process through pybwa, which comes with hicstuff.
  ``--aligner-backend`` (or ``HICSTUFF_ALIGNER_BACKEND``) chooses;
* `minimap2 <https://github.com/lh3/minimap2>`_: the ``minimap2`` program,
  which is not on PyPI either.

Mapping mode ``5spm`` aligns with bwa only, through pybwa unless
``--aligner-backend binary``. ``hicstuff iteralign`` and ``hicstuff missview``
always run the aligner, and ``hicstuff pipeline`` only to start from FASTQ
files: pipeline runs that start from BAM or pairs files need no aligner,
except ``bowtie2-inspect`` (part of bowtie2) when the genome is given as a
bowtie2 index.

No other program is needed. `bgzip <https://www.htslib.org/>`_ compresses the
final pairs file on several threads when it is installed, in the BGZF format,
which pairix can index; pysam compresses it on one thread otherwise. Pairs are
sorted by GNU, BSD or uutils ``sort`` when one is installed, and in Python on
one thread otherwise (busybox ``sort`` lacks the options hicstuff needs); set
``HICSTUFF_SORT=python`` to always sort them in Python. samtools is not
needed: hicstuff runs the samtools that pysam bundles.

Platforms that need a compiler
------------------------------
pip installs wheels where they exist, and builds the other dependencies from
source, which needs a C compiler, make and the zlib headers. As of November
2026: macOS with Python 3.14 (pybwa has no wheel for it yet), macOS on Intel
(biopython has published no Intel wheel since 1.87: ``pip install
--prefer-binary hicstuff`` installs the last version that has one), Linux with
an older glibc such as RHEL or CentOS 7 (pybwa and pysam among others; conda
is easier there), and musl-based Linux such as Alpine (biopython, and
matplotlib on aarch64).

Licence
-------
hicstuff is distributed under the BSD-3-Clause licence. pybwa's metadata says
MIT, but every pybwa wheel compiles in bwa, which is distributed under the
GPLv3: an installation of hicstuff therefore contains GPLv3 code. This matters
only to whoever redistributes it, in a Docker image or a packed environment
for instance, who must then meet the GPLv3's terms for bwa.


Quickstart
==========

The fastest way to generate Hi-C matrices is to use the `hicstuff pipeline` command:

.. code-block:: bash

    hicstuff pipeline -g bt2_index for.fq rev.fq

However, you most likely want to have a look at `the command line help <https://hicstuff.readthedocs.io/en/latest/api/hicstuff.html#hicstuff.commands.Pipeline>`_ to select appropriate options, such as the enzyme used in the experiment. The help can be displayed using:

.. code-block:: bash

    hicstuff pipeline --help

Matrices generated in the default coordinate format can then be visualised using the `view <https://hicstuff.readthedocs.io/en/latest/api/hicstuff.html#hicstuff.commands.View>`_ command, modified using the `rebin <https://hicstuff.readthedocs.io/en/latest/api/hicstuff.html#hicstuff.commands.Rebin>`_ and `convert <https://hicstuff.readthedocs.io/en/latest/api/hicstuff.html#hicstuff.commands.Convert>`_ commands, or used as input for other softwares.
