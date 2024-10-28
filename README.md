## About
This repository includes QUODcarb software compatible with MATLAB for 
 calculating marine CO2-system variables in the exactly-determined and
 over-determined cases. QUODcarb, as presented in Fennell & Primeau 
 (2024), stands for Quantifying Uncertainty in an Over-Determined marine
 CARBonate dataset. The calculation also estimates the output 
 uncertainties and requires input standard measurement uncertainties. See 
 Fennell & Primeau (2024) for methods and introduction to the solver. 
 
 We intend for QUODcarb to have the
 same capabilities as CO2SYSv3, in addition to the following improvements:
 - Allows over-determined input
 	- i.e. two, three, four, or five (and beyond) measurements for the same water parcel may be input into a single QUODcarb calculation
 - Allows multiple input and output temperatures and pressures
 	-  i.e. pH at 25 Celsius and pCO2 at 20 Celsius input in tp(1) and tp(2) for a single QUODcarb calculation
 - Uncertainty calculation inherent to the solver
 	- allows for uncertainty estimates in the over-determined cases that are consistent with the uncertainties in the measurements and thermodynamic model (i.e., the parameterized pKs and mass balance totals)
  - Makes use of the redundancy in over-determined datasets to re-estimate the parameterized dissociation constants and mass balance totals
  	-  pK1, pK2, TB etc., in the over-determined cases only
 - Provides a summary of the joint posterior probability distribution for all the parameters used in the calculation in terms of the location of the maximum and the full covariance matrix.

## History
CO2SYS was initially developed by Lewis and Wallace (1998) for MS DOS, 
 later adapted for MS Excel and MATLAB by Pierrot et al. (2006). The code was 
 vectorized, refined, and optimized for computational speed by van Heuven 
 et al. (2011). Options for error propagation were added by Orr et al. 
 (2018). The CO2SYSv3 update by Sharp (2022, 2023) added numerous options 
 and extended capabilities of the MATLAB software. The QUODcarb software
 builds upon these previous versions. See `resources` for a more in-depth description of the updates to the code by the various previous authors.

## Installation and Use
Download the files of interest in the `src` repository and put them together in a directory
 that is accessible by MATLAB. 
To perform CO2-system calculations, use `src/QUODcarb.m` as directed in this 
documentation and in the comments at the top of the file. The `src/QUODcarb.m` file is a standalone function, it does not need the other files to function. Download the example files only if you are interested in trying them yourself or interested in seeing the example syntax. Download the
 `src/compare.m` file only if you would like to compare the QUODcarb output to
 CO2SYSv3’s output (see `src/example4.m` for example use of `src/compare.m`).

## Syntax Example
QUODcarb requires on input an `obs` structure and an `opt` structure with, respectively, the input conditions and measurements and the options for equilibrium constants etc. 
- `obs` structure with temperature, pressure, etc. may be added to in a for loop at each `obs(i).TC` etc.
	- `obs.sal = salinity;`
	- `obs.usal = sal_uncertainty;` % (± sigma)
	- `obs.TC = total_carbon;` % DIC (umol/kg)
	- `obs.uTC = TC_uncertainty;` % (± sigma)
	- `obs.TA = total_alk;`	% (umol/kg)
	- `obs.uTA = TA_uncertainty;` % (± sigma)
	- `obs.tp(1).T = temp;` % deg Celcius
	- `obs.tp(1).uT = temp_uncertainty;` % (± sigma)
	- `obs.tp(1).P = pressure;` % (dbar)
	- `obs.tp(1).uP = press_uncertainty;` % (± sigma)
	- `obs.tp(1).ph = ph_meas;` % any ph scale
	- `obs.tp(1).uph = ph_uncertainty;` % (± sigma)

- `opt` structure with the input options, the equilibrium constant options are the same as for CO2SYSv3
	- `opt.K1K2 = 10;` % Lueker et al., 2000
	- `opt.KSO4 = 1;` % Dickson et al., 1990a
	- `opt.KF = 2;`	% Perez and Fraga, 1987
	- `opt.TB = 1;`	% Lee et al., 2010
	- `opt.phscale = 1;` % (1=tot, 2=free, 3=sws, 4=nbs)
	- `opt.printcsv	= 1;` % (1=on, 0=off)
	- `opt.fname = 'output.csv';` % CSV filename
	- `opt.printmes = 1;` % (1=on, 0=off)
	- `opt.co2press = 1;` % (1=on, 0=off)
	- `opt.Revelle = 1;` % (1=on, 0=off) Calculate Revelle factor?
 - see the top of the `src/QUODcarb.m` file for all possible input options
 - see `resources/QUODcarb_inputs_outputs.pdf` for a list of the possible input and outut parameters and variable names

## Quick Start
See example scripts (`src/example1.m`, `src/example2.m`, `src/example3.m`) showing possible
 inputs, including correct formulation and capitalization. The function
 QUODcarb requires an input observations structure (`obs`) and an input
 options structure (`opt`). All possible inputs and outputs can be found
 in the document `resources/QUODcarb_inputs_outputs`. To see how we parsed the data
 for the paper, take a look at `paper_scripts/driver.m` (can download the 
 GOMECC-3 dataset via `paper_scripts/parse_gomecc3_data.m` to try it out for yourself).

 View tutorial videos on youtube: https://www.youtube.com/playlist?list=PL8Z4MrPvr_ClapjsVT4urcL6KyktVkyff.

## Recreate Fennell & Primeau 2024 Figures
To recreate the figures from the first QUODcarb paper, see the directory `paper_scripts`. Start by running `paper_scripts/parse_gomecc_data.m` to load the GOMECC-3 Dataset to your local computer. Then run `paper_scripts/driver.m` to run QUODcarb through all 1119 datapoints 26 times, one for each possible calculation. This requires making new directories on your device to store the output .mat files. The code takes ~4hrs to run on a remote device and 20GB of storage. Once the 26 .mat files are loaded you may choose any of the plotting scripts to plot the desired figures.

## Tips for Success
- Capitalization DOES matter, please be careful to capitalize and/or use
     lower case exactly as in the examples. 
	- Example: obs.m(1).ph ≠ obs.m(1).pH
 	- See the full list in `resources/QUODcarb_inputs_outputs.pdf`
- As with CO2SYS, QUODcarb requires two separate CO2-system measurements. The possible inputs include: (input at least two)
	- `obs.TC` (± `obs.uTC`)
	- `obs.TA` (± `obs.uTA`)
	- `obs.tp(1).ph` (± `obs.tp(1).uph`)
	- `obs.tp(2).pco2` (± `obs.tp(2).upco2`)
		- or `obs.tp(2).fco2` ± `obs.tp(2).ufco2`
  		- or `obs.tp(2).co2st` ± `obs.tp(2).uco2st`
	- `obs.tp(3).co3` (± `obs.tp(3).uco3`)
 		- or `obs.tp(3).hco3` ± `obs.tp(3).uhco3`
        - pCO2 and CO3 can also go into tp(1) if they are at the same
          temperature
        - This does not show the salinity, nor the tp(1) through tp(3)
          temperatures and pressures, please do not forget them.
- Even when working with just TC and TA, tp(1) still needs at least a
     temperature and pressure.
        - This is the same as CO2SYS, just formatted differently.
- The input to QUODcarb must always start at `obs(1)` or else it does not
 initialize properly, do not try to give it `obs(2)` with an empty `obs(1)`.
- The input uncertainty must be a non-zero number as -log10(0) is undefined and breaks the code.

## Citation
Not published yet, documentation will be updated when available.

## References
Lewis, E., & Wallace, D. W. R. (1998). Program developed for CO2 system calculations. ORNL/CDIAC-105. Carbon Dioxide Information Analysis Center, Oak Ridge National Laboratory, Oak Ridge, TN.

Sharp, J. D., Pierrot, D., Humphreys, M. P., Epitalon, J. M., Orr, J. C., Lewis, E. R., & Wallace, D. W. R. (2023). CO2SYSv3 for MATLAB (Version v3.2.1). http://doi.org/10.5281/zenodo.3950562

Pierrot, D., Lewis, E., & Wallace, D. W. R. (2006). CO2SYS DOS Program developed for CO2 system calculations. ORNL/CDIAC-105. Carbon Dioxide Information Analysis Center, Oak Ridge National Laboratory, US Department of Energy, Oak Ridge, TN.

Van Heuven, S., Pierrot, D., Rae, J. W. N., Lewis, E., & Wallace, D. W. R. (2011). MATLAB Program developed for CO2 system calculations. ORNL/CDIAC-105b. Carbon Dioxide Information Analysis Center, Oak Ridge National Laboratory, Oak Ridge, TN.

Orr. J. C., Epitalon, J. M., Dickson, A. G., & Gattuso, J. P. (2018). Routine uncertainty propagation for the marine carbon dioxide system. Marine Chemistry 207, 84-107.

## References within the code, as taken from CO2SYSv3 code notes
- K1K2 formulations
	- 1 = Roy et al 1993
 		- Roy et al, Marine Chemistry, 44:249-267, 1993
   		- See also: Erratum, Marine Chemistry 45:337, 1994 and Erratum, Marine Chemistry 52:183, 1996
   		- Typo: In the abstract on p. 249, in the eq. for lnK1* the last term should have S raised to the power 1.5
     		- lnK1 = eq. 29 on p. 254 and what they use in their abstract
       		- lnK2 = eq. 30 on p. 254 and what they use in their abstract    
	- 2 = Goyet and Poisson, 1989
 		- Goyet and Poisson, Deep-Sea Research, 36(11):1635-1654, 1989
   		- pK1 = in Table 5 on p. 1652 and what they use in the abstract
		- pK2 = in Table 5 on p. 1652 and what they use in the abstract
	- 3 = Hansson refit by Dickson and Millero
	  	- Hansson (1973) refit by Dickson and Millero 1987 (DM)
		- Dickson and Millero, Deep-Sea Research, 34(10):1733-1743, 1987
			- See also Corrigenda, Deep-Sea Research 36:983, 1989
		- Refit data of Hansson, Deep-Sea Research, 20:461-478, 1973 and Hansson, Acta Chemica Scandanavia, 27:931-944, 1973
		- Typo in DM on p. 1739 in Table 4: the equation for pK2* for Hansson should have a 0.000132 *S^2 instead of a 0.000116 *S^2
	- 4 = Mehrbach refit by Dickson and Millero
		- Mehrbach (1973) refit by Dickson and Millero (1987)
		- Dickson and Millero, Deep-Sea Research, 34(10):1733-1743, 1987
			- See also Corrigenda, Deep-Sea Research , 36:983, 1989
		- Refit data of Mehrbach et al., Limn Oc, 18(6):897-907, 1973
	- 5 = Hansson and Mehrbach refit by Dickson and Millero
		- Hansson (1973) and Mehrbach (1973) refit by Dickson and Millero (1987)
		- Typo in DM on p. 1740 in Table 5: the second equation should be pK2* = not pK1* = 
	- 6 = GEOSECS (i.e., original Mehrbach) (Takahashi et al., 1982)
		- GEOSECS and Peng et al use K1, K2 from Mehrbach et al., Limnology and Oceanography, 18(6):897-907, 1973
			- i.e. these are the original Mehrbach dissociation constants
	- 7 = Peng (i.e., Mehrbach modified)
		- GEOSECS and Peng et al 1987 use K1, K2 from Mehrbach et al., Limnology and Oceanography, 18(6):897-907, 1973
			- i.e. these are the original Mehrbach dissociation constants
	- 8 = Millero 1979 (pure water)
		- Millero, F. J., Geochimica et Cosmochemica Acta 43:1651-1661, 1979
		- K1 from refit data from Harned and Davis, J American Chemical Society, 65:2030,2037, 1943
		- K2 from refit data from Harned and Scholes, J American Chemical Society, 43:1706-1709, 1941
	- 9 = Cai and Wang 1998
		- Cai and Wang, 1998, for estuarine use
		- Data used in this work is from:
			- K1: Mehrbach (1973) for S>15, for S<15 Mook and Keone (1975)
			- K2: Mehrbach (1973) for S>20, for S<20: Edmond and Gieskes (1970)
	- 10 = Lueker et al 2000
		- Lueker, Dickson, and Keeling, Mar. Chem. 70:105-119, 2000
		- This is Mehrbach’s data refit after conversion to the total scale, for comparison with their equilibrator work 
	- 11 = Mojica Prieto and Millero, 2002
		- Mojica Prieto and Millero 2002. Geochim. Et Cosmochim. Acta. 66(14):2529-2540
	- 12 = Millero et al 2002
		- Millero et al., 2002. Deep-Sea Res. I(49):1705-1723
		- Calculated from overdetermined WOCE-era field measurements
	- 13 = Millero et al 2006
		- Millero, Grahan, Huand, Bustos-Serrano, Pierrot. Mar. Chem. 100:80-94, 2006
		- From Millero 2006 work on pK1 and pK2 from titrations in Gulf Stream seawater
	- 14 = Millero 2010
		- Millero, 2010 for estuarine use
		- Marine and Freshwater Research, v.61 p.139-142
  		- Fits through compilation of real seawater titration results: Mehrbach et al. (1973), Mojica-Prieto & Millero (2002), Millero et al. (2006)
	- 15 = Waters, Millero and Woosley 2014
		- Waters, Millero, and Woosley. Mar. Chem., 165, 66-67, 2014
		- Corrigendum to “The free proton concentration scale for seawater pH”
		- Effectively, this is an update of Millero (2010) formulation
	- 16 = Sulpis et al 2020
		- Sulpis, Lauvset, Hagens. Ocean Science, 16:847-862, 2020
		- This study uses overdeterminations of the carbonate system to iteratively fir K1 and K2
	- 17 = Schockman and Byrne 2021
		- Schockman and Byrne. Geochimica et Cosmochimica Acta, 300:231-245, 2021
		- This study uses spectrophotometric pH measurements to determine K1*K2 and presents a new parameterization for K2 based on these determinations
		- K1 is taken from Waters, Millero, Woosley, 2014
- KSO4 formulations
	- 1 = Dickson 1990a
		- Dickson, A. G., J. Chemical Thermodynamics, 22:113-127, 1990
		- TYPO on p.121: the constant e9 should be e8, this is from eqs 22 and 23 on p. 123, and Table 4 p. 121
	- 2 = Khoo et al., 1977
		- Khoo et al, Analytical Chemistry, 49(1):29-34, 1977
		- They find log(beta) = CO2SYS’s pKS
		- Eq 20 on p. 33
	- 3 = Waters and Millero, 2013
		- Waters and Millero, Marine Chemistry, 149:8-22, 2013 with corrections from Waters et al, Marine Chemistry 165:66-67, 2014 (addendum)
- KF formulations
	- 1 = Dickson and Riley 1979
		- Dickson, A. G. and Riley, J. P., Marine Chemistry 7:89-99, 1979
	- 2 = Perez and Frage 1987
		- Perez and Frage 1987
- TB formulations
	- 1 = Uppstrom, 1979
		- Uppstrom, L. Deep-Sea Research 21:161-162, 1974
	- 2 = Lee et al 2010
		- Lee, Kim, Byrne, Millero, Feely, Yong-Ming Liu, 2010. Geochimica et Cosmochimica Acta 74 (6): 1801-1811
	- For GEOSECS cases: 
		- Culkin, F., in Chemical Oceanography, ed. Riley and Skirrow, 1965
		- This is 1% lower than Uppstrom’s value
- Gas Constant = 83.14462618 ml bar-1 K-1 mol-1
	- Recommended by NIST
	- https://physics.nist.gov/cgi-bin/cuu/Value?r
- Calculate TCa (Total Calcium)
	- Riley, J. P. and Tongudui, M., Chemical Geology 2:263-269, 1967
	- For GEOSECS cases (6 & 7):
		- Culkin, F. in Chemical Oceanography, ed. Riley and Skirrow, 1965
		- Quoted in Takahashi et al, GEOSECS Pacific Expedition v. 3, 1982
- Calculate TF (Total Fluoride)
	- Riley, J. P., Deep-Sea Research 12:219-220, 1965
- Calculate TS (Total Sulfide)
	- Morris, A. W., and Riley, J. P., Deep-Sea Research 13:699-705, 1966
- Calculate K0
	- Weiss, R. F., Marine Chemistry 2:203-215, 1974
	- Partial molal volume of CO2 (cm^3/mol) = 32.3, from Weiss 1974 Appendix, paragraph 3
- Calculate IonS
	- From the DOE handbook, Chapter 5, p. 13/22, eq. 7.2.4, 1994
- Calculate fH
	- Takahashi et al, Chapter 3 in GEOSECS Pacific Expedition, v. 3, 1982 (p. 80)
	- GEOSECS case (7):
		- Peng et al, Tellus 39B:439-458, 1987
		- They reference the GEOSECS report, but round the value given there off so that it is about 0.08 (1%) lower. It doesn’t agree with the check value they give on p. 456.
- Calculate KB
	- Dickson, A. G., Deep-Sea Research 37:755-766, 1990 (aka Dickson 1990b)
	- GEOSECS cases (6 & 7):
		- Lyman, John, UCLA Thesis, 1957
		- Fit by Li et al, JGR 74:5507-5525, 1969
		- This is for GEOSECS and Peng et al.
- Calculate KW
	- Millero, Geochimica et Cosmochemica Acta 59:661-677, 1995
		- His check value of 1.6 umol/kg-SW should be 6.2
	- GEOSECS cases (6 & 7):
		- Millero, Geochemica et Comochemica Acta 43:1651-1661, 1979
	- Freshwater case (8):
		- Millero, Geochimica et Cosmochemica Acta 43: 1651-1661, 1979
- Calculate KP1, KP2, KP3 and KSi
	- Yao and Millero, Aquatic Geochemistry 1:53-88, 1995
		- KP1, KP2, KP3 are on the SWS pH scale in mol/kg-sw
		- KSi was given on the SWS scale in molal units
	- GEOSECS case (7):
		- KP1 = 0.02; Peng et al don’t include the contribution from this term, but it is so small it doesn’t contribute. It needs to be kept to that the routinces work ok. 
		- KP2, KP3 from Kester, D. R., and Pytkowicz, R. M., Limnology and Oceanography 12:243-252, 1967
- Calculate KNH
	- Ammonia dissociation constant from Clegg and Whitfield (1995)
- Calculate KNH4
	- First hydrogen sulfide dissociation constant from Millero et a. (1988)
- Correct K1 K2 for pressure
	- From Millero, 1995. They are the same as Millero 1979 and Millero, 1992. They are from data of Culberson and Pytkowicz, 1968.
	- Millero, F. J., Geochemica et Cosmochemica Acta 59: 661-677, 1995.
		- Typo: a factor of 10^3 was left out of the definition of Kappa
		- Typo: the value of R given is incorrect with the wrong units
		- Typo: the values of the a’s for H2S and H20 are from the 1983 values for fresh water
		- Typo: the value of a1 for B(OH)3 should be +0.1622
		- Table 9 on p. 675 has no values for Si.
		- There are a variety of other typos in Table 9 on p. 675
		- There are other typos in the paper, and most of the check values given don’t check.
	- Millero, F. J., Geochemica et Cosmochemica Acta 43:1651-1661, 1979.
		- See Table 5 and eqs. 7, 7a, 7b on pp. 1656-1657
	- Millero, F. J., Sohn, Mary l., Chemical Oceanography, CRC Press, 1992. See Chapter 6.
		- Typo: this chapter has numerous typos (eqs. 36, 52, 56, 65, 72, 79, and 96 have typos).
	- Culberson, C. H., and Pytkowicz, R. M., Limnology and Oceanography 13:403-417, 1968
	- GEOSECS cases (6 & 7):
		- Takahashi et al., GEOSECS Pacific Expedition v.3, 1982 quotes Culberson and Pytkowicz, 1968, but the fits are the same as those in Edmond and Gieskes, GCA, 34:1261-1291, 1970 who in turn quote Li, personal communication
		- Edmond, John M. and Gieskes, J. M. T. M., Geochemica et Cosmochemica Acta, 34:1261-1291, 1970
- Correct KW for pressure
	- Millero, Chapter 43 in Chemical Oceanography, Academic Press, 1983.
	- Note the temperature dependence of KappaK1 and KappaKW for fresh water in Millero, 1983 are the same
	- GEOSECS case (6 & 7):
		- GEOSECS doesn’t include OH term, so this doesn’t matter. Peng et al didn’t include pressure, but here we assume that the KW correction is the same as for the other seawater cases.
		- from Millero, 1983 and his programs CO2ROY(T).BAS
- Correct KF and KS for pressure
	- From Millero, 1995, which is the same as Millero, 1983
- Correct KP1, KP2, KP3, KSi for pressure
	- From Millero, 1995, which are the same as Millero, 1983
	- The only mention of this is Millero, 1995 where it is stated that the values have been estimated from the values of boric acid. However, there is no listing of the values in the table
	- Use the values for KB from above.
- Fugacity constants
	- In previous versions of CO2SYS, the fugacity factor was calculated assuming pressure at one atmosphere, or close to it. Starting with v3.2.1, an option to use in situ pressure is provided
	- Weiss, R. F., Marine Chemistry 2:203-315, 1974
		- They fit the data of Goff and Gratch (1946) with the vapor pressure lowering by sea salt as given by Robinson (1954)
		- This fits the more complicated Goff and Gratch, and Robinson equations from 273 to 313 deg K and 0 to 40 Sali with a standard error of 0.015%, about 5 uatm over this range
		- Goff, J. A., and Gratch, S., Trans. Am. Soc. Heating and Ventilating Engineers 52:95-122, 1946
		- Robinson, J. Marine Bio Assoc. of the U.K. 33:449-455, 1954
			- Eq. 10 on p. 350, this is in atmospheres
	- GEOSECS and Peng (6 & 7) assume pCO2 = fCO2, or FugFac = 1 (called ‘p2f’ in QUODcarb)
 - Initial pH estimates obtained via approach of Munhoven (2013) added by JD Sharp, CO2SYSv3

## All References, Alphabetical

 Broecker, W. S., Spencer, D. W., & Craig, H (1982). Pacific expedition: hydrographic data 1973-1974. International Decade of Ocean Exploration, National Science Foundation.

Cai, W. J., & Wang, Y. (1998). The chemistry, fluxes, and sources of carbon dioxide in the estuarine waters of the Satilla and Altamaha Rivers, Georgia. Limnology and Oceanography, 43(4), 657-668.

Clegg, S. L., & Whitfield, M. (1995). A chemical model of seawater including dissolved ammonia and the stoichiometric dissociation constant of ammonia in estuarine water and seawater from− 2 to 40 C. Geochimica et Cosmochimica Acta, 59(12), 2403-2421.

Culberson, C. H. & Pytkowicz, R. M. (1968). Effect of pressure on carbonic acid, boric acid, and the pH of seawater, Limnology and Oceanography 13:403-417.

Dickson, A. G., & Riley, J. P. (1979). The estimation of acid dissociation constants in seawater media from potentiometric titrations with strong base. I. The ionic product of water—Kw. Marine Chemistry, 7(2), 89-99.

Dickson, A. G., & Millero, F. J. (1987). A comparison of the equilibrium constants for the dissociation of carbonic acid in seawater media. Deep Sea Research Part A. Oceanographic Research Papers, 34(10), 1733-1743.

Dickson, A. G. (1990). Standard potential of the reaction: AgCl (s)+ 12H2 (g)= Ag (s)+ HCl (aq), and and the standard acidity constant of the ion HSO4− in synthetic sea water from 273.15 to 318.15 K. The Journal of Chemical Thermodynamics, 22(2), 113-127.

Dickson, A. G. (1990). Thermodynamics of the dissociation of boric acid in synthetic seawater from 273.15 to 318.15 K. Deep Sea Research Part A. Oceanographic Research Papers, 37(5), 755-766.

DOE (1994) Handbook of methods for the analysis of the various parameters of the carbon dioxide system in sea water; version 2, A. G. Dickson & C. Goyet, eds. ORNL/CDIAC-74.

Edmond, J. M. & Gieskes, J. M. T. M. (1970). The calculation of the degree of seawater with respect to calcium carbonate under in situ conditions, Geochemica et Cosmochemica Acta, 34:1261-1291.
Goff, J. A., & Gratch, S. (1946). Transactions of the American Society of Heating and Ventilating Engineers, 52:95-122.

Goyet, C., & Poisson, A. (1989). New determination of carbonic acid dissociation constants in seawater as a function of temperature and salinity. Deep Sea Research Part A. Oceanographic Research Papers, 36(11), 1635-1654.

Hansson, I. (1973). The determination of dissociation constants of carbonic acid in synthetic sea water in the salinity range of 20–40‰ and temperature range of 5–30◦ C. Acta Chemica Scandanavia, 27, 931-944.

Hansson, I. (1973, May). A new set of acidity constants for carbonic acid and boric acid in sea water. In Deep Sea Research and Oceanographic Abstracts (Vol. 20, No. 5, pp. 461-478). Elsevier.

Harned, H. S., & Davis Jr, R. (1943). The ionization constant of carbonic acid in water and the solubility of carbon dioxide in water and aqueous salt solutions from 0 to 50. Journal of the American Chemical Society, 65(10), 2030-2037.

Harned, H. S., & Scholes Jr, S. R. (1941). The Ionization Constant of HCO3-from 0 to 50. Journal of the American Chemical Society, 63(6), 1706-1709.

Humphreys, M.P., Lewis, E.R., Sharp, J.D., & Pierrot, D. (2022). PyCO2SYS: marine carbonate system calculations in Python. Geoscientific Model Development 15, 15-43.

Kester, D. R., & Pytkowicx, R. M. (1967). Determination of the apparent dissociation constants of phosphoric acid in seawater 1. Limnology and Oceanography, 12(2), 243-252.

Khoo, K. H., Ramette, R. W., Culberson, C. H., & Bates, R. G. (1977). Determination of hydrogen ion concentrations in seawater from 5 to 40. degree. C: standard potentials at salinities from 20 to 45%. Analytical Chemistry, 49(1), 29-34.

Pierrot, D., Lewis, E., & Wallace, D. W. R. (2006). CO2SYS DOS Program developed for CO2 system calculations. ORNL/CDIAC-105. Carbon Dioxide Information Analysis Center, Oak Ridge National Laboratory, US Department of Energy, Oak Ridge, TN.

Prieto, F. J. M., & Millero, F. J. (2002). The values of pK1+ pK2 for the dissociation of carbonic acid in seawater. Geochimica et Cosmochimica Acta, 66(14), 2529-2540.

Lueker, T. J., Dickson, A. G., & Keeling, C. D. (2000). Ocean pCO2 calculated from dissolved inorganic carbon, alkalinity, and equations for K1 and K2: validation based on laboratory measurements of CO2 in gas and seawater at equilibrium. Marine chemistry, 70(1-3), 105-119.

Lee, K., Kim, T. W., Byrne, R. H., Millero, F. J., Feely, R. A., & Liu, Y. M. (2010). The universal ratio of boron to chlorinity for the North Pacific and North Atlantic oceans. Geochimica et Cosmochimica Acta, 74(6), 1801-1811.

Lewis, E., & Wallace, D. W. R. (1998). Program Developed for CO2 System Calculations. ORNL/CDIAC-105. Carbon Dioxide Information Analysis Center, Oak Ridge National Laboratory, Oak Ridge, TN.

Li, Y. H., Takahashi, T., & Broecker, W. S. (1969). Degree of saturation of CaCO3 in the oceans. Journal of Geophysical Research, 74(23), 5507-5525.

Lyman, J. (1957). Buffer mechanism of seawater. Ph.D Thesis, University of California, Los Angeles 196 pp.

Mehrbach, C., Culberson, C. H., Hawley, J. E., & Pytkowicx, R. M. (1973). Measurement of the apparent dissociation constants of carbonic acid in seawater at atmospheric pressure 1. Limnology and oceanography, 18(6), 897-907.

Millero, F. J. (1979). The thermodynamics of the carbon dioxide system in seawater, Geochemica et Cosmochemica Acta 43:1651-1661.

Millero, Frank J. (1983) Influence of pressure on chemical processes in the sea. Chapter 43 in Chemical Oceanography, eds. Riley, J. P. andChester, R., Academic Press.

Millero, F. J., Plese, T., & Fernandez, M. (1988). The dissociation of hydrogen sulfide in seawater 1. Limnology and Oceanography, 33(2), 269-274.

Millero, Frank J., & Sohn, Mary L., Chemical Oceanography, CRC Press (1992). See chapter 6.

Millero, F. J. (1995). Thermodynamics of the carbon dioxide system in the oceans, Geochemica et Cosmochemica Acta 59:661-677.

Millero, F. J., Pierrot, D., Lee, K., Wanninkhof, R., Feely, R., Sabine, C. L., Key, R. M., & Takahashi, T. (2002). Dissociation constants for carbonic acid determined from field measurements. Deep Sea Research Part I: Oceanographic Research Papers, 49(10), 1705-1723.

Millero, F. J., Graham, T. B., Huang, F., Bustos-Serrano, H., & Pierrot, D. (2006). Dissociation constants of carbonic acid in seawater as a function of salinity and temperature. Marine Chemistry, 100(1-2), 80-94.

Millero, F. J. (2010). Carbonate constants for estuarine waters. Marine and Freshwater Research, 61(2), 139-142.

Mojica Prieto, F. J., & Millero, F. J. (2002). The values of pK1+ pK2 for the dissociation of carbonic acid in seawater. Geochimica et Cosmochimica Acta, 66(14), 2529-2540.

Mook, W. G., & Koene, B. K. S. (1975). Chemistry of dissolved inorganic carbon in estuarine and coastal brackish waters. Estuar. Coastal Mar. Sci., 3: 325-336.

Morris, A. W., & Riley, J. P. (1966, August). The bromide/chlorinity and sulphate/chlorinity ratio in sea water. In Deep sea research and oceanographic Abstracts (Vol. 13, No. 4, pp. 699-705). Elsevier.

Munhoven, G. (2013). Mathematics of the total alkalinity–pH equation – pathway to robust and universal solution algorithms: the SolveSAPHE package v1.0.1. Geoscientific Model Development 6, 1367–1388.

Orr, J. C., Epitalon, J. M., & Gattuso, J. P. (2015). Comparison of ten packages that compute ocean carbonate chemistry. Biogeosciences, 12(5), 1483-1510.

Orr, J.C., Epitalon, J.-M., Dickson, A. G., & Gattuso, J.-P. (2018). Routine uncertainty propagation for the marine carbon dioxide system. Marine Chemistry 207, 84-107.

Peng, T. H., Takahashi, T., Broecker, W. S., & Olafsson, J. O. N. (1987). Seasonal variability of carbon dioxide, nutrients and oxygen in the northern North Atlantic surface water: observations and a model. Tellus B: Chemical and Physical Meteorology, 39(5), 439-458.

Perez, F. F., & Fraga, F. (1987). Association constant of fluoride and hydrogen ions in seawater. Marine Chemistry, 21(2), 161-168.

Riley, J. P. (1965). The occurrence of anomalously high fluoride concentrations in the North Atlantic. Deep Sea Research A, 12(2), 219-220.

Robinson, R. A. (1954). The vapour pressure and osmotic equivalence of sea water. Journal of the Marine Biological Association of the United Kingdom, 33(2), 449-455.

Roy, R. N., Roy, L. N., Vogel, K. M., Porter-Moore, C., Pearson, T., Good, C. E., Millero, F. J., & Campbell, D. M. (1993). The dissociation constants of carbonic acid in seawater at salinities 5 to 45 and temperatures 0 to 45 C. Marine Chemistry, 44(2-4), 249-267.

Schockman, K. M., & Byrne, R. H. (2021). Spectrophotometric determination of the bicarbonate dissociation constant in seawater. Geochimica et Cosmochimica Acta, 300, 231-245.

Sharp, J.D., Pierrot, D., Humphreys, M.P., Epitalon, J.-M., Orr, J.C., Lewis, E.R., & Wallace, D.W.R. (2023, Jan. 19). CO2SYSv3 for MATLAB (Version v3.2.1). Zenodo. http://doi.org/10.5281/zenodo.3950562

Skirrow, G., & Riley, J. P. (Eds.). (1965). Chemical oceanography. Academic Press.

Sulpis, O., Lauvset, S. K., & Hagens, M. (2020). Current estimates of K1* and K2* appear inconsistent with measured CO2 system parameters in cold oceanic regions. Ocean Science Discussions, 1-27.

Takahashi, T., Williams, R., & Bos, D. (1982) Chapter 3: Carbonate Chemistry. In W. S. Broecker, D. W. Spencer, H. Craig (Eds.), Pacific Expedition: Hydrographic Data 1973-1974, 77-82. International Decade of Ocean Exploration, National Science Foundation.

Uppström, L. R. (1974). The boron/chlorinity ratio of deep-sea water from the Pacific Ocean. Deep Sea Research A, 21(2), 161-162.

van Heuven, S., Pierrot, D., Rae, J.W.B., Lewis, E., & Wallace, D.W.R. (2011). MATLAB Program Developed for CO2 System Calculations. ORNL/CDIAC-105b. Carbon Dioxide Information Analysis Center, Oak Ridge National Laboratory, Oak Ridge, TN.

Waters, J. F., & Millero, F. J. (2013). The free proton concentration scale for seawater pH. Marine Chemistry, 149, 8-22.

Waters, J., Millero, F. J., & Woosley, R. J. (2014). Corrigendum to “The free proton concentration scale for seawater pH”,[MARCHE: 149 (2013) 8–22]. Mar. Chem, 165, 66-67.

Weiss, R. (1974). Carbon dioxide in water and seawater: the solubility of a non-ideal gas. Marine chemistry, 2(3), 203-215.

Yao, W., & Millero, F. J. (1995). The chemistry of the anoxic waters in the Framvaren Fjord, Norway. Aquatic Geochemistry, 1, 53-88.
![image](https://github.com/user-attachments/assets/4f664959-bbca-46db-a533-df75421d8ee4)



