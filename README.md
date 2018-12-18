[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/fastphylo/README.html)

# Fastphylo #

Fastphylo is software project containing the implementations of the algorithms "Fast Computation of
Distance Estimators", "Fast Neighbor Joining", and more. See the paper "FastPhylo: Fast tools for phylogenetics" for more information.

The software is licensed under the MIT license.


## Installation ##

There is presently one easy option for installing FastPhylo and that is using BioConda. So,

1. [Set up Bio Conda for your system](https://bioconda.github.io/index.html), including adding the bioconda channel.
2. [Run the install command](https://bioconda.github.io/recipes/fastphylo/README.html): ```conda install fastphylo```

The alternative is to build the software locally from source.

## Source code ##

Source code is available through cloning of this repository -- see GitHub links and resources.
See the file "INSTALL" for compilation and installation instructions!

## Documentation ##

For the moment, we refer to the docs at our old web site: fastphylo.sourceforge.net

### Relevant publications

* Elias & Lagergren, [Fast Computation of Distance Estimators](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-8-89), BMC Bioinformatics, 2007.
* Elias & Lagergren, [Fast Neighbor Joining](https://www.sciencedirect.com/science/article/pii/S0304397508009079?via%3Dihub), Theoretical Computer Science, 2009.
* Khan _et al_, [FastPhylo: Fast tools for phylogenetics](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-334), BMC Bioinformatics, 2013.

-------------------------------------------------------------------------------

The directory layout of this package:

* src/c++       the c++ sources
* src/docbook   the docbook sources for the html documentation on the homepage

