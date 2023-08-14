---
title: 'pvlib python: 2023 update'
tags:
  - Python
  - solar energy
  - photovoltaics
  - renewable energy
authors:
  - name: Kevin S. Anderson
    orcid: 0000-0002-1166-7957
    affiliation: 1
  - name: William F. Holmgren
    orcid: 0000-0001-6218-9767
    affiliation: 2
  - name: Clifford W. Hansen
    orcid: 0000-0002-8620-5378
    affiliation: 1
  - name: Mark A. Mikofski
    orcid: 0000-0001-8001-8582
    affiliation: 2
  - name: Adam R. Jensen
    orcid: 0000-0002-5554-9856
    affiliation: 3
  - name: Anton Driesse
    orcid: 0000-0003-3023-2155
    affiliation: 4
affiliations:
 - name: Sandia National Laboratories
   index: 1
 - name: DNV
   index: 2
 - name: Technical University of Denmark
   index: 3
 - name: PV Performance Labs
   index: 4

date: 20 July 2023
bibliography: paper.bib
---


# Summary

pvlib python is a community-developed, open-source software toolbox
for simulating the performance of solar photovoltaic (PV) energy
systems.  It provides peer-reviewed reference implementations of
over 100 empirical and physics-based models from the scientific literature,
including solar position algorithms, irradiance models, thermal models,
and PV electrical models.  In addition to these individual low-level
model implementations, pvlib python provides high-level constructs
that chain these models together like building blocks to form complete
"weather-to-power" photovoltaic system models.  It also provides functions to
fetch and import a wide variety of weather datasets useful for PV modeling. 

pvlib python has been developed since 2015 and follows modern best
practices for open-source python software, with comprehensive automated
testing, standards-based packaging, and semantic versioning.
Its source code is developed openly on GitHub and releases are
distributed via the Python Package Index (PyPI) and the conda-forge
repository.  pvlib python's source code is made freely available under
the permissive BSD-3 license.

Here we present an update on pvlib python describing capability
and community development since our 2018 publication [@pvlibjoss2018].


# Statement of need

PV performance models are used throughout the field of solar phovoltaics.
The rapid increase in scale, technological diversity, and sophistication
of the global solar energy industry demands correspondingly more
capable models.  Per the United States Department of Energy,
"the importance of accurate modeling is hard to overstate" [@seto2022].

Compared with other PV modeling tools, pvlib python stands out in several
key aspects.  One is its reusable toolbox design, providing the user a
level of flexibility beyond that of other tools.  Rather than organizing
the user interface around pre-built modeling workflows, pvlib python
makes the individual "building blocks" of PV performance models accessible to
the user.  This allows the user to assemble their own model workflows, including
the ability of inserting their own custom modeling steps.  This flexibility
is essential for applications in both academia and industry.

Another key aspect of pvlib python is that it is used via
a general-purpose programming language (Python).  In addition to being more
generally flexible, powerful, and scalable than a traditional
graphical user interface, pvlib python's interface as a programming library
allows it to be combined with other Python packages.  This enables
integration with database query, data manipulation, numerical optimization,
plotting, and reporting packages, to name a few.

A final key aspect of pvlib python is its open peer review approach and
foundation in published scientific research, allowing it to be developed by
a decentralized and diverse community of PV researchers and practitioners
without compromising its focus on transparent and reliable model
implementations.

These core tenets, along with sustained contributions from a passionate and
committed community, have lead to its widespread adoption across the PV
field [@Stein2022].  In support of the claim that pvlib python provides
meaningful value and addresses real needs, we list these quantitative metrics:

1. Its 2018 JOSS publication, at the time of this writing,
ranks 14th by citation count out of the 2000+ papers published by JOSS to date.
2. The Python Package Index (PyPI) classifies pvlib python as "critical project"
due to being in the top 1% of the index's packages by download count.
3. The project's online documentation currently receives over 400,000 page views
per year.
4. pvlib python was found to be the third most used python project
in the broader open-source sustainability software landscape, with the first
two being netCDF4 utilities applicable across many scientific fields [@Augspurger2023].


# Functionality additions

To meet new needs of the PV industry, substantial new functionality has been
added in the roughly five years since the 2018 JOSS publication.

First, several dozen new models have been
implemented, expanding the package's capability in both existing and new
modeling areas and prompting the creation of several
new modules within pvlib python.
Notable additions include methods of fitting empirical PV performance models
to measurements, models for performance loss mechanisms like soiling and snow
coverage, and models for simulating rear-side irradiance (prompted by the recent rapid
increase in bifacial PV deployment).
\autoref{fig:functions-comparison} summarizes the number of models (or functions)
per module for pvlib python versions 0.6.0 and 0.10.1, showing a substantial
capability expansion over the last five years.

![Comparison of public function counts for selected pvlib modules for v0.6.0 and v0.10.1. Some modules are smaller in v0.10.1 due to moving functions to new modules (e.g. from `pvsystem` to `iam`).\label{fig:functions-comparison}](functions_06_010.png)

Second, in addition to the new function-level model implementations,
the package's high-level classes have also been expanded to better support
the complexity of emerging system designs, including systems with differing
mounting or electrical parameters across subsystems and systems that require
custom orientation/steering models.

Third, the creation of `pvlib.iotools`, a sub-package for fetching and parsing
datasets relevant to PV modeling.  These functions provide a standardized
interface for reading data files in various complex data formats, offering
conveniences like optionally standardizing the dataset labels and units
to pvlib's conventions [@Jensen2023].  As of version 0.10.1, `pvlib.iotools` contains
functions to download data from ten online datasets,
plus file reading/parsing functions for a dozen solar resource file formats.

These additions are discussed in more detail in [@pvpmc_2023_update] and [@pvpmc2022_pvlib_update].
Complete descriptions of the changes in each release can be found in the
project's documentation.


# Community growth

It is difficult or impossible to comprehensively describe the community around
open-source projects like pvlib python, but some aspects of it are more easily
quantifiable than others.  Here we examine the community from a few convenient
perspectives, but we emphasize that this section is necessarily a limited view of
the community as a whole.

First, we examine contributors to pvlib python's code repository.  The
project's use of version control software enables easy quantification of
repository additions (to code, documentation, tests, etc) over time.  The
project's repository currently comprises contributions from over 100 people
spanning industry, academia, and government research institutions.
\autoref{fig:community} shows the number of unique repository
contributors over time, demonstrating continued and generally accelerating
attraction of new contributors.

![Total repository contributor count over time (left) and other community size statistics (right).\label{fig:community}](community.png)

However, the project as a whole is the product of not only of those who contribute
code but also those who submit bug reports, propose ideas for new features,
participate in online fora, and support the project in other ways.
Along those lines, two easily tracked metrics are the number of people
registered in the pvlib python online discussion forum and the number of
GitHub "stars" on the pvlib python code repository.  \autoref{fig:community}
shows these counts over time.  Although these numbers
almost certainly substantially underestimate the true size of the pvlib community,
their increase over time indicates continued and accelerating community growth.

In addition to continuous interaction online, community members sometimes
meet in person at user's group and tutorial sessions run by pvlib python
maintainers and community members alike.
To date, these meetings have been held at the IEEE Photovoltaics Specialists
Conference (PVSC), the PVPMC Workshops, and the PyData Global conference.
\autoref{fig:timeline} shows a timeline of these meetings, along with other
notable events in the project's history.

![pvlib python major event timeline: releases (top), community events (middle), and other project milestones (bottom).\label{fig:timeline}](timeline2.png)

Finally, it is worth pointing out that pvlib python contributors and users
are part of a broader community around not just pvlib python itself
but also other members of the pvlib software "family": pvanalytics, a
package for PV data quality assurance and feature recognition
algorithms [@pvpmc2022_pvanalytics_update], and twoaxistracking, a package
for simulating self-shading in arrays of two-axis solar trackers [@Jensen2022].
Moreover, looking beyond pvlib and its affiliated packages, we see that Python
is proving to be the most common programming language for open-source PV
modeling and analysis software in general.  The packages mentioned here
make up one portion of a growing landscape of Python-for-PV projects [@Holmgren2018].


# Acknowledgements

Although much of the development and maintenance of pvlib python is on a
volunteer basis, the project has also been supported by projects with various
funding sources, including:

- The U.S. Department of Energy’s Solar Energy Technology Office, through
  the PV Performance Modeling Collaborative (PVPMC) and other projects
- The Danish Energy Agency through grant no. 64020-1082
- NumFOCUS's Small Development Grant program
- Google's Summer of Code program

pvlib python benefits enormously from building on top of
various high-quality packages that have become de facto standards in the python
ecosystem: numpy [@numpy], pandas [@pandas], scipy [@scipy], and numba [@numba]
for numerics, matplotlib [@matplotlib] for plotting,
sphinx [@sphinx] for documentation, and pytest [@pytest] for automated testing.
The project also benefits from online infrastructure generously provided free
of charge, including GitHub (code development and automated testing) and
ReadTheDocs.org (documentation building and hosting).

This work was supported by the U.S. Department of Energy’s Office of Energy
Efficiency and Renewable Energy (EERE) under the Solar Energy Technologies
Office Award Number 38267. Sandia National Laboratories is a multimission
laboratory managed and operated by National Technology & Engineering Solutions
of Sandia, LLC, a wholly owned subsidiary of Honeywell International Inc., for
the U.S. Department of Energy’s National Nuclear Security Administration under
contract DE-NA0003525. This paper describes objective technical results and
analysis. Any subjective views or opinions that might be expressed in the paper
do not necessarily represent the views of the U.S. Department of Energy or
the United States Government.

# NOTES (to be deleted)

Additional authors? 

- PVPMC?  Josh, Marios, Dan
- Tutorials?  Silvana
- Former pvlib python contributors?  Rob, Tony, Cameron, etc

The future?  Commercial adoption?

# References
