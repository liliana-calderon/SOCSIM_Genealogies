Code for reproducible results for the paper “Analyzing biases in
genealogies using demographic microsimulation”
================

- [Introduction](#introduction)
- [Scripts for reproducible results](#scripts-for-reproducible-results)

# Introduction

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

# Scripts for reproducible results

**NB**: If you want to re-run the analysis from scratch (i.e., retrieve
the data from [HFC](https://www.fertilitydata.org/),
[HFD](https://www.humanfertility.org/), and
[HMD](https://www.mortality.org/)), run the microsimulations and
reproduce the results), you need R, some standard packages, `HMDHFDplus`
and `rsocsim`. Please, follow these steps:

1.  Download this repository to your computer.  
2.  Register at [HFD](https://www.humanfertility.org/) and
    [HMD](https://www.mortality.org/) and install `HMDHFDplus`.  
    To write the input rate files and to compare the simulation results
    with the original input data, we use data from
    [HFC](https://www.fertilitydata.org/),
    [HFD](https://www.humanfertility.org/), and
    [HMD](https://www.mortality.org/). The last two are retrieved
    through the `HMDHFDplus` package. To gain full access to the
    databases, you need to become a registered user at both databases,
    after accepting the user agreement on their websites:
    <https://www.humanfertility.org/Account/Auth> and
    <https://mortality.org/Account/Auth>.
3.  Install `rsocsim`.  
    Please go to <https://github.com/MPIDR/rsocsim> to find the
    instructions on how to install `rsocsim`.  
4.  Open the provided RStudio project.
5.  Run the scripts sequentially in the order suggested by the numbers
    of the files (0-6).
6.  In scripts 0 and 2, please type your username and password for
    [HFD](https://www.humanfertility.org/) and
    [HMD](https://www.mortality.org/) in the corresponding space before
    running the code.  
7.  To get exactly the same simulation results, you should use the same
    randomly generated seeds mentioned in script 1.