# Analyzing biases in genealogies using demographic microsimulation

Genealogies are promising sources for addressing many questions in historical and kinship demography. So far, an incomplete understanding of the biases that affect their representativeness has hindered their full exploitation. Here, we conduct a series of experiments on synthetic populations aimed at understanding how different sources of bias in ascendant genealogies can affect the accuracy of demographic estimates. We use the SOCSIM demographic microsimulation program and data for Sweden from the Human Fertility Collection (HFD) (1751-1890), the Human Fertility Database (HFD) (1891-2022), and the Human Mortality Database (HMD) (1751-2022). We analyze three sources of bias: selection in direct lineages, incomplete reconstruction of family trees, and missing information on some subpopulations. We evaluate their effect by comparing common demographic measures estimated from ‘perfectly-recorded’ and ‘bias-infused’ synthetic populations.


# Set up

To run this code, you need R, some standard packages, HMDHFDplus and rsocsim. 
Please go to <https://github.com/MPIDR/rsocsim>  to find information on how to install rsocsim.

To write the input rate files from HFD and HMD and to compare the simulation results with the original input data, 
we use the HMDHFDplus package. For this, you need to create an account ag 
Please type your username and password for HMD and HFD in the corresponding space in scripts 0 and 1.

To reproduce the analysis, please run the scripts numbered 0 to 6. 
