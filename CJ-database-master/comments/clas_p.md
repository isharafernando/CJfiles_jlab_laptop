# Proton and Deuteron F2 and Reduced Cross Sections from CLAS6

## Reference:
* Proton
    * Data table    : CLAS-NOTE [2003-001](https://arxiv.org/abs/hep-ex/0309052)
    * uncertainties : Table I in Phys. Rev. D 67, 092001
* Deuterium:  
   * Phys.Rev. C73 (2006) 045205
   * CLAS-NOTE [2005-013](https://arxiv.org/abs/hep-ex/0507098v2). 

## Source: 
Emails from Mikhail Osipenko.


## Data files: 
  * F2      proton   : [xlsx](../data/JAM/10057.xlsx), [csv](../data/JAM/csv/10057.csv)   
  * F2      deuteron : [xlsx](../data/JAM/10058.xlsx), [csv](../data/JAM/csv/10058.csv)   
  * sig_r   proton   : [xlsx](../data/JAM/10059.xlsx), [csv](../data/JAM/csv/10059.csv)   
  * sig_r   deuteron : [xlsx](../data/JAM/10060.xlsx), [csv](../data/JAM/csv/10060.csv)   


## Uncertainties:
**Words on systematic uncertainties (from Osipenko):

F2p has even systematic uncertainties split into different sources. Unfortunately there is no "overall normalizzation" part, however because CLAS is not very precise detector the dominant uncertainty is point-to-point correlated (due to wide acceptance and complex detector structure). The overall part is perhaps 1-3%, but we never separated it.

1) eff - efficiency+acceptance, the data are corrected by means of GEANT3 Monte Carlo simulation of CLAS detector, these are rough estimate of possible mismatch in detector description and reality,

2) eep - e+e- contamination correction, i.e. mostly e- coming from pi0 decay,

3) phe - photo-electron correction, means efficiency of Cherenkov counter of CLAS used to identify electrons,

4) rdc - radiative correction (i.e. model used in calculations),

5) mom - momentum correction, small adjustments of electron momentum due to non-ideal knowledge of CLAS magnetic field,

6) rlt - for F2 only, model dependence due to assumed R=sigma_L/sigma_T ratio used in the extraction.

All these corrections were applied to the data in tables. Uncertainties, being an estimate of unknown error, are just rough guess of how much these corrections could be wrong.

All these uncertainties are __point-to-point correlated__, but in unknown way, what means that if eg. the model of radiative correction or R_LT was wrong by 10% all data points will be affected together, the magnitude of effect will be similar to the given uncertainty, but the specific kinematic function could be different.

These uncertainties are independent of each other, you can sum them in quadrature."

*Additional normalization uncertainties:
According to Ref.3 the F2p raw yield is normalized by luminosity then checked against the ep elastic cross section model (see section C-E). We know the ep cross section at 1% level, so it's reasonable to add a 1% nomalization uncertainty to all F2 and reduced cross section data.


