Tutorial
========

Introduction
------------

The Taiji pipeline aims at integrating different kinds of high throughput profiling techniques to construct TF regulatory networks and identify key regulators through network analysis.
That being said, one only needs ATAC-seq or DNase-seq data to run the analysis, though the result will be better given other information, especially RNA-seq data.

The data input to the Taiji pipeline are fastq or bam files.
For gene expression data, you have the option to provide an external gene expression table instead of having the pipeline analyze the raw RNA-seq data for you.

How to use
----------

To run the Taiji pipeline, you would need 2 configuration files.

The first configuration file is used to specify the options used by the pipeline.
Please look at the `example configuration file <https://github.com/kaizhang/Taiji/blob/master/example_config.yml>`_ for details.

The second configuration file contains the information about the input data sets.
Take a look at this `example <https://github.com/kaizhang/Taiji/blob/master/example_input.yml>`_.

To run the pipeline, supply the ``taiji`` with the first configuration file: ``taiji run --config example_config.yml`` or ``taiji run --config example_config.yml --remote`` if the program is compiled with ``sge`` flag.

Parallelism
-----------

Taiji supports two levels of parallelism -- node level and workflow level. Node
level parallelism is automatically turned on when compiling with the ``SGE`` flag.
The workflow level parallelism can be turned on using `-N <num_of_process>`.
However, this is only recommended for users with a super computer, as it will
consume a lot of memory.

Auto-recovery
-------------

The pipeline supports auto-recovery, which means you can stop the program at any time and it will resume from the last checkpoint.
The checkpoints are saved in a file called "sciflow.db".
Delete this file if you want a fresh run.

Results
-------

The Taiji pipeline outputs many files, distributed in several directories.

``OUTPUTDIR/Rank/``
^^^^^^^^^^^^^^^^^^^

This is the primary output of Taiji pipeline. It contains the TF ranks under different
conditions / cell-types. Cell-type-specific TF can be identified by looking at
the TF rank dynamics (fold-change) across different cell types.

- `GeneRank_all.tsv`: Ranks for all genes.
- `GeneRank_filtered.tsv`: Genes with lower rank and that are less variate are removed.

``OUTPUTDIR/TFBS/``
^^^^^^^^^^^^^^^^^^^

BED files of predicted TF binding sites in each cell type.

``OUTPUTDIR/Network/``
^^^^^^^^^^^^^^^^^^^^^^

Static gene regulatory network for each cell type.

``OUTPUTDIR/ATAC_Seq/``
^^^^^^^^^^^^^^^^^^^^^^^

BAM, BED, and peak files for ATAC-seq data.

``OUTPUTDIR/RNA_Seq/``
^^^^^^^^^^^^^^^^^^^^^^

BAM, gene and transcript quantification files for RNA-seq data.


Visualize the results
---------------------

We will release the codes of data visualization soon. Stay tuned...
