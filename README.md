Analyzing biases in genealogies using demographic microsimulation
================

- <a href="#set-up" id="toc-set-up">Set up</a>

## Introduction

Genealogies are promising sources for addressing many questions in
historical and kinship demography. So far, an incomplete understanding
of the biases that affect their representativeness has hindered their
full exploitation. Here, we conduct a series of experiments on synthetic
populations aimed at understanding how different sources of bias in
ascendant genealogies can affect the accuracy of demographic estimates.
We use the SOCSIM demographic microsimulation program and data for
Sweden from the [Human Fertility Collection
(HFC)](https://www.fertilitydata.org/) (1751-1890), the [Human Fertility
Database (HFD)](https://www.humanfertility.org/) (1891-2022), and the
[Human Mortality Database (HMD)](https://www.mortality.org/)
(1751-2022). We analyze three sources of bias: selection in direct
lineages, incomplete reconstruction of family trees, and missing
information on some subpopulations. We evaluate their effect by
comparing common demographic measures estimated from
‘perfectly-recorded’ and ‘bias-infused’ synthetic populations.

# Set up

To retrieve the data, run the microsimulations and reproduce the results
of the paper (scripts numbered 0 to 6), you need R, some standard
packages, HMDHFDplus and rsocsim.

To write the input rate files (script 0) and to compare the simulation
results with the original input data (script 2), we use data from
[HFC](https://www.fertilitydata.org/),
[HFD](https://www.humanfertility.org/), and
[HMD](https://www.mortality.org/). The last two are retrieved through
the HMDHFDplus package. To gain full access to the databases, you need
to become a registered user, after accepting the user agreement at
[HFD](https://www.humanfertility.org/Account/Auth) and
[HMD](https://mortality.org/Account/Auth). Once you obtain your
credentials and have HMDHFDplus installed, please type your username and
password for HMD and HFD in the corresponding space in scripts 0 and 2.

To find information on how to install rsocsim, please go to
<https://github.com/MPIDR/rsocsim>
