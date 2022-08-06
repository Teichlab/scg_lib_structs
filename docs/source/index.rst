Single Cell Genomics Library Preprocessing Pipelines
====================================================

As more and more methods are developed, the list in the `scg_lib_structs <https://github.com/Teichlab/scg_lib_structs>`_ GitHub repository grows bigger and bigger. To extract information from the data, the first step is always data preprocessing, a procedure that converts the raw data ( ``fastq`` files ) to some sort of count matrices, such as gene-by-cell or peak-by-cell matrices. If you Google it, you will find some tutorials on this topic, and different experimental methods have quite different preprocessing pipelines. With the development of the computational tools, it is now possible to perform the preprocessing procedures using a unified pipeline (sort of ...). This documentation showcases how to perform the data preprocessing step using just `starsolo <https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md>`_ for all scRNA-seq methods and `chromap <https://github.com/haowenz/chromap>`_ + `MACS <https://github.com/macs3-project/MACS>`_ for all scATAC-seq methods. This documentation does not really provide a ready-to-use pipeline for each method. Instead, it documents the commands used in the pipeline and, more importantly, explains the reason and the rationale behind the chosen parameters of each software. The point here is to showcase how and why we do this so that you can build and customise the pipeline on your own.

.. note::

   - **Feedback needed !!!** This project is still under development and will be updated according to my own time. If you have questions, spot any errors, see something confusing and have suggestions for improvement, please do get in touch by `raising an issue <https://github.com/Teichlab/scg_lib_structs/issues>`_ in the **scg_lib_structs** GitHub repository, or by email: ``chenx9@sustech.edu.cn``.
   - If you go through a few methods, you will find some text is repeating. The reason is that I want to make each method as a self-contained and independent page. I want people to be able to just click a method they are interested in and immediately start reading, without having to read other methods first.

.. tip::

   Make sure you are familiar with different sequencing modes from different Illumina machines by looking at `this page <https://teichlab.github.io/scg_lib_structs/methods_html/Illumina.html>`_.

Required softwares
------------------

The softwares needed for the preprocessing is very standard. All of them are used routinely in genomics and can be installed via ``conda`` or the like, so I'm not going to talk about software installation. I'm stating the version of each software I'm currently using (24-Jul-2022). Also make sure they are executable and in your ``$PATH``.

* General utilities

.. code-block:: text

   curl v7.79.1
   wget v1.20.3
   samtools v1.13
   bedtools v2.30.0
   tabix v0.2.5
   bgzip v0.2.5
   bedClip
   faSize
   bcl2fastq v2.20.0.422 (only if you want to practice generating FastQ)
   sratoolkit v3.0.0 (only for GEO data)

* scRNA-seq

.. code-block:: text

   STAR v2.7.9a

* scATAC-seq

.. code-block:: text

   chromap v0.2.1-r369
   MACS2 v2.2.7.1

Note that ``bedClip`` and ``faSize`` are from the `UCSC genome browser executables <http://hgdownload.soe.ucsc.edu/admin/exe/>`_.

Before you start:
-----------------

.. toctree::
   :maxdepth: 2

   Building_reference.md


Methods:
--------

.. toctree::
   :maxdepth: 2

   Gene_expression.md
   Epigenetics.md
   Multi-omics.md
