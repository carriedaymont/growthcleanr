 
 
 *********Modifying infant algo for Aim 2
 *2025-10-02
 *Modified from last 2023 code. See source and mods below. 

 
 *Mods from 2025-10-01: changed to do all CHOP data
 
 *
*****  Change directory
 *

cd "/Users/carriedaymont/Dropbox/Dropbox Documents/Research/__Chris-growthcleanr/Share"

**************************************************************************************
*Step 0: Dataset Preparation
**************************************************************************************

 *
******  Using all of CHOP data
 *
 
 
timer clear // Clears all existing timers
timer on 1 // Starts timer 1

use chopdall, clear


*Dataset specific stuff
*Fix param
replace param="WEIGHTKG" if param=="Weight"
replace param="HEIGHTCM" if param=="Height"
replace param="HEADCM" if param=="Head Circumference"
replace param="GESTAGE" if param=="Gestational Age"

*Fix birth head circumference in CHOP data
replace measurement=(measurement / 0.0283495) * 2.54 if param=="HEADCM" & agedays==0 & dataset=="chop"

*/

**************************************************************************************
*Step 1: Dataset Preparation
**************************************************************************************


******************************       1A       ****************************

*ID variable

*Make id variable
sort subjid param agedays measurement
gen id=_n



******************************       1B       ****************************

*Variables for ageyears and agemonths

*Make agemonths and ageyears
gen ageyears=agedays/365.25
gen agemonths=agedays/30.4375

******************************       1C       ****************************

*Drop those 20 and older
keep if agedays<20*365.25

*Drop if no values
drop if param!="WEIGHTKG" & param!="HEIGHTCM" & param!="HEADCM"
drop if measurement==.

*Drop if HC over 3y
drop if param=="HEADCM" & agedays>3*365.25


*May be Stata-specific, generate param-specific measurement variables

*Make wt ht variables
gen wt=measurement if param=="WEIGHTKG"
gen ht=measurement if param=="HEIGHTCM"
gen hc=measurement if param=="HEADCM"

*Make whichvar variable
gen whichvar=""



save temp_chopsetup, replace

 *
*****		Exporting 1p dataset - will change this when doing whole dataset
 *
export delimited temp_chopsetup_1p.csv, replace
*/

**************************************************************************************
*Step 2a: Z-scores
**************************************************************************************

*Calculate both CDC and WHO z and SD from growth data files and use SD from both

*This foreach collapse in 1{ just facilitates collapsing sections of code to make it easier to work with

foreach collapse in 1 {

use temp_chopsetup, clear

***2aA**  Merge with growthdata files, using 2023-01-27 version
merge m:1 agedays sex using growthfile_cdc_ext, gen(check1)
drop if check1==2
merge m:1 agedays sex using growthfile_who, gen(check2)
drop if check2==2
drop check*

***2aB*** Make modified z-scores (using same method as last time, when they were called SD scores)
foreach c in who cdc {
	foreach p in wt ht hc{
		gen `c'`p'z=(`p'-`c'_`p'_m)/`c'_`p'_csd_pos
		replace `c'`p'z=(`p'-`c'_`p'_m)/`c'_`p'_csd_neg if `p'<`c'_`p'_m
	}
}
 
***2aC***
*For wt and ht smooth from 2y to 4y, simple weighted average
*Note for ht, length/ht transition also happens at 2y
gen whoweight=4-ageyears
gen cdcweight=ageyears-2
foreach p in wt ht {
	gen s`p'z=(who`p'z * whoweight + cdc`p'z * cdcweight)/2
	replace s`p'z=who`p'z if ageyears<=2
	replace s`p'z=cdc`p'z if ageyears>=4
}

*Don't smooth for HC
gen shcz=whohcz

***2aD*** Indicate these are original csd's -- need for later
rename who_* orig_who_* 
rename cdc_* orig_cdc_*

save temp_chopz, replace
export delimited temp_chopz.csv, replace


}


*/

**************************************************************************************
*Step 2b: Corrected z-scores
**************************************************************************************

use temp_chopz, clear

/*
For first weight for a subject that is before 10 completed months AND swtz <-2
*Identify which Fenton M it is closest to, assign that as its "gestational age"
*Determine postmenstrual agedays (assigned ga in days + chron age in days) 
*Determine ga-corrected agedays: agedays - (280-assigned ga)
*If firstwt is higher than median wt for Fenton at 50 weeks, 
	*determine which corrected age would put firstwt at median for WHO
*Until chron age 2y, use ga-corrected z-score
	-Fenton z-score for CGA <=50 weeks
	-Age-corrected WHO for CGA >50 weeks
*From 2-4y chron age, calc ga-corrected z-score and uncorrected z-score and smooth
*For 4y and up, use uncorrected z-score
*Check whether corrected or original first wtz are consistent with other early z
	*If corrected first wtz more consistent, use corrected z's
	*If original first wtz more consistent, first weight was likely in error, use original z
*/

**2bA** Identify first weights in subjects potentially eligible for correction
sort subjid param agedays measurement
by subjid param: gen sn_wt=_n
gen potcorr_wt=sn_wt==1 & swtz<-2 & agemonths<10
bysort subjid: egen potcorr=max(potcorr_wt)

**2bB** Convert to integer grams rounded to nearest 10 to facilitate merging
gen intwt=int(wt*100) * 10
*Smallest wt in Fenton is 510 for week 22, but ht and hc only go to week 23 (wt=570)
*Replace intwt with 570 if intwt is 250 to 560
replace intwt=570 if intwt>=250 & intwt<=560

**2bC** Merge with Fenton to get assigned GA (aga) based on p
merge m:1 intwt sex using fentlms_foraga, gen(check_intwtmerge)
drop if check_intwtmerge==2
*Note not correctable if did not merge (weight too low)
replace potcorr_wt=0 if check_intwtmerge==1
drop potcorr
bysort subjid: egen potcorr=max(potcorr_wt)
gen aga_i=fengadays if potcorr_wt==1
bysort subjid: egen aga=min(aga_i)

**2bD** Determine pma and corrected agedays if <37 weeks (<259 days)
gen pmagedays=aga + agedays if aga<259
gen cagedays=pmagedays-280 if aga<259
replace cagedays=. if ageyears>4
replace aga=. if aga>=259
replace fengadays=pmagedays

**2bE** Merge with Fenton
merge m:1 fengadays sex using fentlms_forz, gen(check_fenlms)
drop if check_fenlms==2
drop check_fenlms


**2bF** Determine Fenton z-scores using pma
*Minimal skew so just use regular LMS method
*Need to make wt into grams for this
gen wtg=wt*1000
rename fen_wt_* fen_wtg_*

foreach p in wtg ht hc {
	gen fen`p'z=(((`p'/fen_`p'_m)^fen_`p'_l)-1) / (fen_`p'_l * fen_`p'_s)
}
rename fenwtgz fenwtz

**2bG** Make ga-corrected z-scores
*Merge with WHO again using corrected age
rename agedays agedaysorig
rename cagedays agedays
merge m:1 agedays sex using growthfile_who, gen(checkwho)
drop if checkwho==2
merge m:1 agedays sex using growthfile_cdc_ext, gen(checkcdc)
drop if checkcdc==2
rename agedays cagedays
rename agedaysorig agedays
rename who_* cwho_*
rename cdc_* ccdc_*
*If agedays >730 but cagedays<=730, then LMS will be ht for agedays but lt for cagedays
*Therefore, add 0.8 to ht if agedays>730 and cagedays<=730 for WHO
*Therefore, add 0.7 to ht if agedays>730 and cagedays<=730 for CDC
gen cwho_cht=ht
replace cwho_cht=ht+0.8 if agedays>730 & cagedays<=730
gen ccdc_cht=ht
replace ccdc_cht=ht+0.7 if agedays>730 & cagedays<=730

foreach c in cwho ccdc {
	foreach p in wt hc{
		gen `c'`p'z=(`p'-`c'_`p'_m)/`c'_`p'_csd_pos
		replace `c'`p'z=(`p'-`c'_`p'_m)/`c'_`p'_csd_neg if `p'<`c'_`p'_m
	}
}

foreach c in cwho ccdc {
	foreach p in ht {
		gen `c'htz=(`c'_cht-`c'_`p'_m)/`c'_`p'_csd_pos
		replace `c'htz=(`c'_cht-`c'_`p'_m)/`c'_`p'_csd_neg if `p'<`c'_`p'_m
	}
}
drop *_cht check*

*For potentially correctable, make corrected z-score = 
	*Fenton z-score when available
	*Corrected WHO z-score from 0-2 when available
foreach p in wt ht hc {
	gen c`p'z=fen`p'z if potcorr==1
	replace c`p'z=cwho`p'z if potcorr==1 & c`p'z==. & ageyears<=2
}


**2bH** Make smoothed corrected z-scores for WHO/CDC from 2-4
*simple weighted average

*(Note do not do this for HC)
gen cwhoweight=4-ageyears
gen ccdcweight=ageyears-2
foreach p in wt ht {
	replace c`p'z=(cwho`p'z * cwhoweight + ccdc`p'z * ccdcweight)/2 if ageyears>2
}
replace chcz=cwhohcz if ageyears>2

**2bI** Make smoothed z-scores between ga-corrected and uncorrected smoothed from 2-4
gen uncorrweight=4-ageyears
gen corrweight=ageyears-2
foreach p in wt ht hc {
	*Smooth
	gen sc`p'z=(s`p'z * uncorrweight + c`p'z * corrweight)/2
	*For <2, use corrected Fenton/WHO for correctable and available, otherwise use regular z-score
	replace sc`p'z=c`p'z if ageyears<=2 & potcorr==1
	replace sc`p'z=s`p'z if c`p'z==.
	*Keep original z if ageyears>=4
	replace sc`p'z=s`p'z if ageyears>=4
	*Keep original z if not eligible for correction
	replace sc`p'z=s`p'z if potcorr==0	
}

**2bJ** Done at different points above, but makes sense to describe/double check that it's done down here


**2bK** Check for whether first weight is likely in error by comparing differences with corrected zs for next 3 weights
gen firstwtage_i=potcorr_wt==1
bysort subjid: egen firstwtage=min(firstwtage_i)
foreach n in 2 3 4{
	foreach p in wt {
		gen s`p'z`n'_i=s`p'z if sn==`n' & ageyears<2 & potcorr==1
		bysort subjid: egen s`p'z`n'=max(s`p'z`n'_i)
		gen sc`p'z`n'_i=sc`p'z if sn==`n' & ageyears<2 & potcorr==1
		bysort subjid: egen sc`p'z`n'=max(sc`p'z`n'_i)
		gen ds`p'z`n'_i=s`p'z-s`p'z`n' if potcorr==1 & sn==1
		bysort subjid: egen ds`p'z`n'=max(ds`p'z`n'_i)		
		gen dsc`p'z`n'_i=sc`p'z-sc`p'z`n' if potcorr==1 & sn==1
		bysort subjid: egen dsc`p'z`n'=max(dsc`p'z`n'_i)
	}
}

*Do not correct if there are other values and if |sum of diffs| is > for corrected than for uncorrected
foreach t in s sc {
	gen abssumdiff_`t'=abs(d`t'wtz2 + d`t'wtz3 + d`t'wtz4)
	replace abssumdiff_`t'=abs(d`t'wtz2 + d`t'wtz3) if abssumdiff_`t'==.
	replace abssumdiff_`t'=abs(d`t'wtz2) if abssumdiff_`t'==.	
}
gen uncorr_i=abssumdiff_sc>abssumdiff_s & abssumdiff_s!=. & potcorr_wt==1
bysort subjid: egen uncorr=max(uncorr_i)

*Replace sc`p'z with s`p'z for these values
foreach p in wt ht hc {
	replace sc`p'z=s`p'z if uncorr==1
}

*Recenter corrected with everyone else in next step

save temp_corrz, replace
export delimited temp_corrz.csv, replace



*/


**************************************************************************************
*Step 3: Recentering (using external file)
**************************************************************************************


foreach collapse in 1 {

use temp_corrz, clear

**3A** Merge with rc data file

merge m:1 agedays sex param using rcfile-2023-08-15, gen(check1)



*Recenter z-scores, making "to be cleaned" variables

*Rename to match Stata code
gen rcwtz = sdmedian if param=="WEIGHTKG"
gen rchtz = sdmedian if param=="HEIGHTCM"
gen rchcz = sdmedian if param=="HEADCM"

foreach p in wt ht hc {
	gen tbc`p'z=s`p'z-rc`p'z
}

**3B** Recenter corrected z-scores
foreach p in wt ht hc {
	gen ctbc`p'z=sc`p'z-rc`p'z
}


**For Carrie's workflow**
*Set up exclusion variable and parameter-specific subject ID
tostring(subjid), gen(subjid_str)
foreach p in wt ht hc{
	gen exc_`p'=0
	replace exc_`p'=1 if `p'==.
	gen subjid_`p'=subjid_str
	replace subjid_`p'="" if exc_`p'!=0
}

lab def exc 0 "Include" 1 "Missing"
lab val exc_* exc

}

save temp_rc, replace




*/

*
*
*
**************************************************************************************
*Testers (for testing, not a step)
**************************************************************************************

*
*
*


foreach collapse in 1 {


*Identify testers for specific things

use temp_rc, clear

*Label truth variable
lab val truth exc

*Identify testers from sample

replace testernote="" if subjid==1


****General
replace testernote="Lots of ht exclusions" if subjid==273608

****Should not be excluded - Weight
replace testernote="Plaus WT - all are plausible, swtz up to 18, changed BIV for wt to 22" if subjid==415033
replace testernote="Plaus WT - all are plausible" if subjid==391947
replace testernote="Plaus WT - all are plausible" if subjid==77925
replace testernote="Plaus WT - all are plausible" if subjid==552576
replace testernote="Plaus WT - all are plausible" if subjid==439679
replace testernote="Plaus WT - all are plausible" if subjid==605828
replace testernote="Plaus WT - all are plausible" if subjid==214765
replace testernote="Plaus WT - all are plausible" if subjid==101376
replace testernote="Plaus WT - all are plausible" if subjid==382250
*replace testernote="" if subjid==1
*replace testernote="" if subjid==1
*replace testernote="" if subjid==1
*replace testernote="" if subjid==1

gen plauswt=subjid==415033 | subjid==391947 | subjid==77925 | subjid==552576 | subjid==439679 | subjid==605828 | ///
	subjid==214765 | subjid==101376 | subjid==382250
replace plauswt=. if wt==.



*****Borderline - weight
replace testernote="Borderline WT - 1501" if subjid==59883
replace testernote="Borderline WT - 741" if subjid==427543
replace testernote="Borderline WT - 69, think is real" if subjid==336691
replace testernote="Borderline WT - 780, 985" if subjid==135373
*replace testernote="" if subjid==1
*replace testernote="" if subjid==1
*replace testernote="" if subjid==1
*replace testernote="" if subjid==1


****Evil Twins
replace testernote="ET WT - multiple wt unit errors in a row" if subjid==415033
replace testernote="ET WT - implausibly high birth weight" if subjid==518862
replace testernote="ET WT - pair of high weights later, ?plausible" if subjid==292515
replace testernote="ET WT - classic wt unit error evil twin" if subjid==98785
replace testernote="ET WT - true wt z-score difference 8 from birth to 4 years" if subjid==160383
replace testernote="ET WT - BW z-score 7, at 2 years was -0.4, extremely unlikely" if subjid==262498
replace testernote="ET WT - true increase from wt z 1.07 at birth to 6 at 552 days" if subjid==269867
replace testernote="ET WT - 4 wt errors in a row" if subjid==703
replace testernote="ET WT - should exclude higher z-scores rather than lower, excludes 7/16" if subjid==415033
replace testernote="ET HT - evil twins" if subjid==234717
replace testernote="ET HT - evil twins may be corrected but should be excluded later" if subjid==341961
replace testernote="ET HT - evil twins may be corrected but should be excluded later" if subjid==517758
replace testernote="Pair HT" if subjid==175941
replace testernote="ET HT - z scores 2.3 (birth), -4 (4 years), -1.2, -1.2; unclear which is right" if subjid==17139
replace testernote="ET HC - z-score of -0.5 that is extreme error" if subjid==552196
replace testernote="ET HC - z-score of -3.8 that is extreme error" if subjid==95308
replace testernote="ET HT - 2nd and 3rd hts probably correct, first height is probably high but not necessarily implaus" if subjid==572912


****EWMA 1
replace testernote="EWMA1 WT - exclude 616" if subjid==104084
replace testernote="EWMA1 WT - exclude 439" if subjid==242948
replace testernote="EWMA1 WT - exclude 439" if subjid==263972
replace testernote="EWMA1 WT - exclude 7" if subjid==620197
replace testernote="EWMA1 WT - exclude 2992" if subjid==468634
replace testernote="EWMA1 WT - exclude 6083 (SDEs, both bad)" if subjid==346417
replace testernote="EWMA1 WT - exclude 283" if subjid==451136
replace testernote="EWMA1 WT - exclude 181" if subjid==10826
replace testernote="EWMA1 WT - exclude 322" if subjid==179911
*replace testernote="" if subjid==1
*replace testernote="" if subjid==1
*replace testernote="" if subjid==1
*replace testernote="" if subjid==1
*replace testernote="" if subjid==1
*replace testernote="" if subjid==1


*****SDE
replace testernote="SDE WT - SDE-EWMA d75" if subjid==61928
*replace testernote="" if subjid==1
*replace testernote="" if subjid==1
*replace testernote="" if subjid==1
*replace testernote="" if subjid==1
*replace testernote="" if subjid==1



******EWMA 2
replace testernote="EWMA2 HC - ??hydro inc from -0.14 to 1.9 in 1 month,  almost 6cm" if subjid==6395308
replace testernote="EWMA2 WT - BW 2.353, swtz inc by 2, gain 32g/d x 124d, high but plaus" if subjid==617824
replace testernote="EWMA2 WT - may be real but probably appropriate to exclude" if subjid==605828
replace testernote="EWMA2 HT - 2 days probably exclude, birth prob wrong too" if subjid==280707
replace testernote="EWMA2 HT - 2 days exclude" if subjid==17704
replace testernote="EWMA2 HT - 3 days exclude" if subjid==1
replace testernote="EWMA2 HT - 3 days exclude" if subjid==105356
replace testernote="EWMA2 HT - 3 days prob exclude" if subjid==419123
replace testernote="EWMA2 HT - 3 days exclude" if subjid==221405
replace testernote="EWMA2 HT - exclude 1st and 3rd measurements" if subjid==76970
replace testernote="EWMA2 HT - 3 days exclude" if subjid==290278
replace testernote="EWMA2 HT - 3 days exclude" if subjid==498586
replace testernote="EWMA2 HT - 3 days exclude" if subjid==370522
replace testernote="EWMA2 HT - 4 days exclude" if subjid==384691
replace testernote="EWMA2 HT - 4 days exclude" if subjid==101206
replace testernote="EWMA2 HT - 4 days exclude" if subjid==470727
replace testernote="EWMA2 HT - 5 days exclude" if subjid==350926
replace testernote="EWMA2 HT - 5 days exclude" if subjid==133295
replace testernote="EWMA2 HT - 17 days exclude (not 5)" if subjid==605128
replace testernote="EWMA2 HT - 6 days exclude" if subjid==636212
replace testernote="EWMA2 HT - 7 days exclude" if subjid==73097
replace testernote="EWMA2 HT - 8 days exclude" if subjid==546898
replace testernote="EWMA2 HT - 11 days exclude" if subjid==437369
replace testernote="EWMA2 HT - 18 days exclude" if subjid==377874
replace testernote="EWMA2 HT - 31 days exclude" if subjid==591554
replace testernote="EWMA2 HT - 34 days exclude" if subjid==602549
replace testernote="EWMA2 HT - birth exclude (not second measurement)" if subjid==28813
replace testernote="EWMA2 HT - 1577 days exclude" if subjid==17139
replace testernote="EWMA2 HT - 368, 517, 655 days exclude" if subjid==399756
replace testernote="EWMA2 HT - birth prob exclude" if subjid==405645
replace testernote="EWMA2 HT - 715 days exclude" if subjid==407136
replace testernote="EWMA2 HT - 17 days exclude" if subjid==411316
replace testernote="EWMA2 HT - 665 days exclude" if subjid==413795
replace testernote="EWMA2 HT - ALL include" if subjid==416385
replace testernote="EWMA2 HT - 3 days exclude" if subjid==419123
replace testernote="EWMA2 HT - birth exclude" if subjid==422233
replace testernote="EWMA2 HT - 61 days exclude" if subjid==432158

replace testernote="EWMA2 WT - 663 days exclude" if subjid==301323


*replace testernote="EWMA2 HT -  days exclude" if subjid==


*Maxdiff
replace testernote="Maxdiff HT" if subjid==556959
replace testernote="Maxdiff HT" if subjid==561636
replace testernote="Maxdiff HT" if subjid==69860
replace testernote="Maxdiff HT" if subjid==233126

******Error load
replace testernote="Error load WT" if subjid==59119
*replace testernote="" if subjid==1
*replace testernote="" if subjid==1
*replace testernote="" if subjid==1
*replace testernote="" if subjid==1


*Tester variable
replace tester=1 if testernote!=""

}

save temp_testers, replace




*/


*
*
*
**************************************************************************************
*Subset (for testing, not a step)
*For smaller datasets this should not remove anyone
**************************************************************************************

*
*
*


*Select subset of data
*Start with random sample
*Can change seed for sample

use temp_testers, clear


*Usually keep those with early/birth meas only but need later measurements to test correction and other things
keep if ind_hc==0 | ind_hc==1 | synthetic==1

sort subjid param agedays measurement


set seed 6
gen ran=runiform()
sort subjid
bysort subjid: gen sn=_n
sort sn ran
gen gnum_i=_n if sn==1
bysort subjid: egen gnum=min(gnum_i)
**# Bookmark #1
********************************************************************************************************************************
keep if gnum<10000000 | tester==1 | synthetic==1

drop ran gnum*


save temp_subset, replace


*Remove files for space
foreach f in chopsetup corrz {
	capture confirm file "temp_`f'.dta"
	if !_rc {  
		rm temp_`f'.dta
	}
}


*/

**************************************************************************************
*Step 4: No need to EWMA
**************************************************************************************

use temp_subset, clear

**4A** Identify clean subjects with no need to EWMA (nnte), will join back up for absolute value checks

foreach p in wt ht hc {
	
	*|Diff| between min and max tbc`p'z < 1
	*Try alternate criteria - difference between min and max <2.5, no adjacents more than one different 
	bysort subjid param: egen maxtbc`p'z=max(tbc`p'z)
	bysort subjid param: egen mintbc`p'z=min(tbc`p'z)
	gen difftbc`p'z=maxtbc`p'z-mintbc`p'z
	
	*Calculate difference between adjacent z-scores AND measurements 
	sort subjid param agedays
	by subjid param: gen diffprevtbc`p'z=tbc`p'z-tbc`p'z[_n-1]

	*No SDEs
	duplicates tag subjid param agedays, gen(dupsde_`p')
	
	*Crude CF - no identical measurements
	duplicates tag subjid param measurement, gen (dupcf_`p')

	*NNTE if all criteria above are met
	gen nnte_`p'_i=tbc`p'z>-3 & tbc`p'z<3 & abs(difftbc`p'z)<2.5 & (abs(diffprevtbc`p'z)<1 | diffprevtbc`p'z==.) & dupsde_`p'==0 & dupcf_`p'==0
	
	*Subject-level parameter-specific variable
	bysort subjid param: egen nnte_`p'=min(nnte_`p'_i)
	replace nnte_`p'=1 if `p'==.
}

*Overall variable indicating does not need to be EWMAed
gen nnte_i=0 if nnte_wt==0 | nnte_ht==0 | nnte_hc==0
replace nnte_i=1 if nnte_i==.
bysort subjid: egen nnte=min(nnte_i)

lab def nnte 0 "Should be EWMAed" 1 "No need to EWMA"
lab val nnte nnte_wt nnte_ht nnte_hc nnte

*drop mintbc*z maxtbc*z difftbc*z diff* dup* nnte_*

**Carrie's workflow to identify these values to be able to addend them later** 


save temp_nnte, replace

use temp_nnte, clear

keep if nnte==1

save temp_nnte_forappend, replace

use temp_nnte, clear

keep if nnte==0

save temp_tbc, replace



rm temp_nnte.dta


*/
**************************************************************************************
*Step 5: Temporary SDEs
**************************************************************************************


use temp_tbc, clear

*Calculate medians for all parameters
foreach p in wt ht hc {
	duplicates tag subjid_`p' param agedays, gen(dup_`p')
	bysort subjid_`p': egen median_`p'z_i=median(tbc`p'z) if exc_`p'==0 & dup_`p'==0
	bysort subjid_`p': egen median_`p'z=min(median_`p'z_i)
}

*Determine differences from medians
local oi="1"
foreach p in wt ht hc {
	gen temp`p'=`oi'
	local oi=`oi'+1
}
foreach p in wt ht hc {
	local o empty
	if temp`p'==1 {
		local o ht	
	}
	if temp`p'==2 {
		local o wt
	}
	if temp`p'==3 {
		local o ht
	}
	gen absd_median_`p'z=abs(tbc`p'z-median_`p'z)
	gen absd_median_`o'z=abs(tbc`p'z-median_`o'z)
	sort subjid_`p' agedays absd_median_`p'z absd_median_`o'z
	by subjid_`p' agedays: gen d_sort_`p'_n=_n if exc_`p'==0
	replace exc_`p'=2 if exc_`p'==0 & d_sort_`p'_n>1
	replace subjid_`p'=subjid_str
	replace subjid_`p'="" if exc_`p'!=0
	drop absd_median_*  d_sort*
}

drop dup_* median*z* temp*
lab def exc 2 "SDE (temp)", add


save temp_tempsde, replace



*/


**************************************************************************************
*Step 6: Carried Forwards
**************************************************************************************


use temp_tempsde, clear

*Note: took out unit error/switch


foreach p in wt ht hc {

**6A**
	*If no duplicates
	sort subjid_`p' agedays `p'
	by subjid_`p': gen d_prev_`p'=measurement-measurement[_n-1]
	replace exc_`p'=3 if d_prev_`p'==0 & exc_`p'==0
	replace subjid_`p'=subjid_str
	replace subjid_`p'="" if exc_`p'!=0
	drop d_prev_`p'
	
**6B**
	*If duplicates
	tostring(subjid), gen(subjid_`p'2)
	replace subjid_`p'="" if exc_`p'==1
	sort subjid_`p'2 agedays exc_`p'
	by subjid_`p'2 agedays: gen tempn=_n if subjid_`p'2!=""
	by subjid_`p'2 agedays: gen temptot=_N if subjid_`p'2!=""
	sort subjid_`p'2 tempn agedays
	by subjid_`p'2 tempn: gen agen_i=_n if tempn==1 & subjid_`p'2!=""
	bysort subjid_`p'2 agedays: egen agen=min(agen_i)
	drop agen_i
	gen agen_prev=agen-1 if agen>1
	sort subjid_`p'2 tempn agen
	by subjid_`p'2 tempn: gen tempn_prev_i=tempn[_n-1]
	bysort subjid_`p'2 agen: egen tempn_prev=min(tempn_prev_i)
	drop tempn_prev_i
	egen maxtempn=max(tempn)
	egen maxagen=max(agen)
	local mtempn=maxtempn[1]
	local magen=maxagen[1]
	forvalues tn=1/`mtempn' {
		gen d_prevn_`tn'=.
	}
	forvalues an=2/`magen' {
		forvalues tn=1/`mtempn' {
			gen tempcomp_i=measurement if agen==(`an'-1) & tempn==`tn' & subjid_`p'2!=""
			bysort subjid_`p'2: egen tempcomp=min(tempcomp_i)
			replace d_prevn_`tn'=abs(measurement-tempcomp) if agen==`an' & subjid_`p'2!=""
			drop tempcomp tempcomp_i
		}
	}
	egen mincomp=rowmin(d_prevn_*) if subjid_`p'2!="" & agen>2
	replace exc_`p'=3 if mincomp==0 & exc_`p'!=1
	replace subjid_`p'=subjid_str
	replace subjid_`p'="" if exc_`p'!=0
	drop subjid_`p'2 mincomp temp* *prev* *agen* *tempn* mincomp
}


**6C** Redo temporary SDEs
*Calculate medians for all parameters
foreach p in wt ht hc {
	duplicates tag subjid_`p' param agedays, gen(dup_`p')
	bysort subjid_`p': egen median_`p'z_i=median(tbc`p'z) if exc_`p'==0 & dup_`p'==0
	bysort subjid_`p': egen median_`p'z=min(median_`p'z_i)
}

*Determine differences from medians
local oi="1"
foreach p in wt ht hc {
	gen temp`p'=`oi'
	local oi=`oi'+1
}
foreach p in wt ht hc {
	local o empty
	if temp`p'==1 {
		local o ht	
	}
	if temp`p'==2 {
		local o wt
	}
	if temp`p'==3 {
		local o ht
	}
	gen absd_median_`p'z=abs(tbc`p'z-median_`p'z)
	gen absd_median_`o'z=abs(tbc`p'z-median_`o'z)
	sort subjid_`p' agedays absd_median_`p'z absd_median_`o'z
	by subjid_`p' agedays: gen d_sort_`p'_n=_n if exc_`p'==0
	replace exc_`p'=2 if exc_`p'==0 & d_sort_`p'_n>1
	replace subjid_`p'=subjid_str
	replace subjid_`p'="" if exc_`p'>0
	drop absd_median_*  d_sort*
}

*Identify different categories of carried forward

**6D**
*Determine if value is in whole pounds
gen lbs=wt*2.20462262
gen wholelbs=abs(mod(lbs, 1)<0.01) & wt!=.
*Determine if ht value is in whole/half inches
gen htin=ht/2.54
gen htwholein=abs(mod(htin, 1)<0.01) & ht!=.
gen hthalfin=abs(mod(htin, 0.5)<0.01) & ht!=.
*Determine if hc value is in whole/half inches
gen hcin=ht/2.54
gen hcwholein=abs(mod(hcin, 1)<0.01) & hc!=.
gen hchalfin=abs(mod(hcin, 0.5)<0.01) & hc!=.
*Determine if wholeimp
gen wholeimp=wholelbs==1 | htwholein==1 | hcwholein==1
*Determine if whole/halfimp
gen wholehalfimp=wholelbs==1 | htwholein==1 | hcwholein==1 | hthalfin==1 | hchalfin==1




foreach p in wt ht hc {
	*Temporarily reinclude cfs in subjid
	replace subjid_`p'=subjid_str if exc_`p'==3
	sort subjid_`p' agedays
	
	**6E** Determine length of CF string
	sort subjid_`p' agedays
	by subjid_`p': gen cfstr_`p'_i=exc_`p'==3
	forvalues n=2/50 {
		local m=`n'-1
		sort subjid_`p' agedays
		by subjid_`p': replace cfstr_`p'_i=exc_`p'==`n' if exc_`p'==3 & cfstr_`p'_i[_n-1]==`m'
	}
	by subjid_`p': egen cfstr_`p'=max(cfstr_`p'_i)
	drop *_i
	
	**6F** *Determine difference between initial value and CF z-score
	*Only go up to length of string
	*Use nonrecentered z-score
	gen zdiff_`p'=.
	sort subjid_`p' agedays
	forvalues n=1/50 {
		sort subjid_`p' agedays
		by subjid_`p': replace zdiff_`p'=s`p'z[_n-`n']-s`p'z if cfstr_`p'>=`n'
	}
	
	*Identify those meeting other criteria
	**6G**
	*Determine if difference between first CF and initial value is below <0.05 and CF string <2
	replace exc_`p'=4 if exc_`p'==3 & cfstr_`p'==1 & zdiff_`p'<0.05
	
	**6H**
	*Determine if difference between first CF and initial value is below 0.1 and whole/half imperial and CF string <2
	replace exc_`p'=5 if (exc_`p'==3 | exc_`p'==4 ) & cfstr_`p'==1 & zdiff_`p'<0.1 & wholehalfimp==1
	
	**6I**
	*Determine if difference between first CF and initial value is below tighter z-limit (0.05) and teen
	replace exc_`p'=6 if exc_`p'==3 & zdiff_`p'<0.05 & ((sex==0 & ageyears>=17) | (sex==1 & ageyears>=16))
	
	**6J**
	*Determine if difference between first CF and initial value is below 0.01 and whole/half imperial and teen
	replace exc_`p'=7 if (exc_`p'==3 | exc_`p'==6) & zdiff_`p'<0.1 & wholehalfimp==1 & ((sex==0 & ageyears>=17) | (sex==1 & ageyears>=16))
	
	*Clean-up*
	*Replace subjid_`p'
	replace subjid_`p'=subjid_str
	replace subjid_`p'="" if exc_`p'==3 | exc_`p'==4 | exc_`p'==5 | exc_`p'==6 | exc_`p'==7
}


*Clean-up*
drop dup_* median*z* temp*
lab def exc 3 "Carried forward" 4 "1 CF deltaZ <0.05" 5 "1 CF imperial deltaZ <0.1" 6 "Teen 2 plus CF deltaZ <0.05" 7 "Teen 2 plus CF imperial deltaZ <0.1", add 


save temp_cf, replace



*Remove files for space
foreach f in tbc {
	capture confirm file "temp_`f'.dta"
	if !_rc {  
		rm temp_`f'.dta
	}
}

*/

**************************************************************************************
*Step 7: BIVs
**************************************************************************************


*Wt limits from do-file Oct 10 2022, in dataset, highest plausible weight=29.5kg, 65 lbs
*HC Limits based on analysis in do-file from Oct 11 2022: in dataset, lowest plausible 18.5, z=-13, highest plausible 63.5, z=11.2  (at birth, highest=42, z=6)


use temp_cf, clear

**7A**
*Identify absolute cutoffs
*Change to BIV exclusion code even if temp SDE or CF

*Min/max weight at birth based on published births 
*Note, even extremely small wts (0.004) are only -8.61 z-score at birth
replace exc_wt=44 if wt<0.2 & exc_wt!=1
replace exc_wt=44 if wt>10.5 & agedays==0 & exc_wt!=1
*Min weight after birth based on rules about who can go home with alowance for significant weight loss -- this is for outpatient data only so very low weight babies would still be in NICU
replace exc_wt=44 if wt<1 & agedays>0 & exc_wt!=1
*Max weight for <2y based on eval (see do-file from Oct 10 2022), highest plaus weight 29.5kg
replace exc_wt=44 if wt>35 & ageyears<2 & exc_wt!=1
*Max weight for all based on published data
replace exc_wt=44 if wt>600 & exc_wt!=1


*Min/max ht based on analysis in do file from 
*Also, 18 is z=-6 for 22 0/7 in Fenton and 65 is z=6 for 40 0/7
replace exc_ht=44 if ht<18 & exc_ht!=1
replace exc_ht=44 if ht>65 & agedays==0 & exc_ht!=1
replace exc_ht=44 if ht>244 & exc_ht!=1


*Min/max HC based on analysis in do file from Oct 11 2022
*Also, 13 is z=-6 for 22 0/7 in Fenton and 
replace exc_hc=44 if hc<13 & exc_hc!=1
replace exc_hc=44 if hc>75 & hc!=. & exc_hc!=1
replace exc_hc=44 if hc>50 & agedays==0 & exc_hc!=1

**7B**
*Z-scores - ***Note, using unrecentered values***
*For weight only do after birth
replace exc_wt=55 if swtz<-25 & ageyears<1 & exc_wt!=1 & exc_wt!=44
replace exc_wt=55 if swtz<-15 & ageyears>=1 & exc_wt!=1 & exc_wt!=44
replace exc_wt=55 if swtz>22 & exc_wt!=1 & exc_wt!=44

*Max z-score for height based on analysis of CHOP data because 15/25 too loose for upper limits
replace exc_ht=55 if shtz<-25 & ageyears<1 & exc_ht!=1 & exc_wt!=44
replace exc_ht=55 if shtz<-15 & ageyears>=1 & exc_ht!=1 & exc_wt!=44
replace exc_ht=55 if shtz>8 & exc_ht!=1 & exc_wt!=44

*HC
replace exc_hc=55 if shcz<-15 & exc_hc!=1 & exc_wt!=44
replace exc_hc=55 if shcz>15 & exc_hc!=1 & exc_wt!=44

foreach p in wt ht hc {
	replace subjid_`p'=subjid_str
	replace subjid_`p'="" if exc_`p'!=0
}

*Redo temporary SDEs
*Calculate medians for all parameters
foreach p in wt ht hc {
	duplicates tag subjid_`p' param agedays, gen(dup_`p')
	bysort subjid_`p': egen median_`p'z_i=median(tbc`p'z) if exc_`p'==0 & dup_`p'==0
	bysort subjid_`p': egen median_`p'z=min(median_`p'z_i)
}

*Determine differences from medians
local oi="1"
foreach p in wt ht hc {
	gen temp`p'=`oi'
	local oi=`oi'+1
}
foreach p in wt ht hc {
	local o empty
	if temp`p'==1 {
		local o ht	
	}
	if temp`p'==2 {
		local o wt
	}
	if temp`p'==3 {
		local o ht
	}
	gen absd_median_`p'z=abs(tbc`p'z-median_`p'z)
	gen absd_median_`o'z=abs(tbc`p'z-median_`o'z)
	sort subjid_`p' agedays absd_median_`p'z absd_median_`o'z
	by subjid_`p' agedays: gen d_sort_`p'_n=_n if exc_`p'==0
	replace exc_`p'=2 if exc_`p'==0 & d_sort_`p'_n>1
	replace subjid_`p'=subjid_str
	replace subjid_`p'="" if exc_`p'!=0
	drop absd_median_*  d_sort*
}

drop dup_* median*z* temp*
lab def exc 44 "Absolute BIV" 55 "Standardized BIV", add 


save temp_biv, replace



*Remove files for space
foreach f in tempsde {
	capture confirm file "temp_`f'.dta"
	if !_rc {  
		rm temp_`f'.dta
	}
}

*/

**************************************************************************************
*Step 9: Evil Twins
**************************************************************************************



/*Overall strategy
-Compare each z to one before/after
-If |delta z| is >5, delete the value further away from the median tbcz for the subject (worked better than ewma or variation)
-If corrected z is available has to be excluded by both corr and uncorr to be exclued
-Different from most steps in that once a |delta z|>5 is identified, at least one value for that subject/parameter for that round will get excluded
-Continue until no big delta z's are identified
-Don't do pairs at this stage
*/

****Adjusting limits
*Weight: one set of unit errors was caught but was just on cusp of z-score difference of 5
*Weight: one true increase of just over 8 (that was over 4 years, may need to stratify at 2 years or so)
*Weight: one true increase of just over 5 (from 1 at birth to 6 at 552 days - maybe stratify at 1 year?)
*Height: one where more extreme value got included, it was probably (although not definitely) the wrong one

use temp_biv, clear


foreach p in wt ht hc {
	
	local i=1
	while `i'>0 {

		gen hinext`p'=0
		gen hiprev`p'=0
		
		**9A**
		*Create variable for difference in z-score
		*Uncorrected
		sort subjid_`p' agedays
		by subjid_`p': gen `p'zdiffnext=tbc`p'z[_n+1] - tbc`p'z
		by subjid_`p': gen `p'zdiffprev=tbc`p'z - tbc`p'z[_n-1]
		
		**9B**
		*Corrz setup
		gen cexc_`p'=exc_`p'
		replace cexc_`p'=20 if ctbc`p'z==.
		gen csubjid_`p'=subjid_`p'
		replace csubjid_`p'="" if ctbc`p'z==.
		
		*Create variable for difference in corrected z-score
		sort csubjid_`p' agedays
		by csubjid_`p': gen `p'czdiffnext=ctbc`p'z[_n+1] - ctbc`p'z
		by csubjid_`p': gen `p'czdiffprev=ctbc`p'z - ctbc`p'z[_n-1]

		**9C**
		*Indicator for Out of Bounds
		gen outbounds_`p'=abs(`p'zdiffnext)>5 & (abs(`p'czdiffnext)>5 | `p'czdiffnext==.) & `p'zdiffnext!=. & hinext`p'==0
		replace outbounds_`p'=1 if abs(`p'zdiffprev)>5 & `p'zdiffprev!=. & (abs(`p'czdiffprev)>5 | `p'czdiffnext==.) & hiprev`p'==0	
		
		**9D**
		*Determine which one is further away from median z for subject
		bysort subjid_`p': egen median`p'z=median(tbc`p'z)
		gen absd`p'z=abs(tbc`p'z-median`p'z)
		gsort subjid_`p' outbounds_`p' -absd`p'z
		by subjid_`p' outbounds_`p': gen outn_`p'=_n if outbounds_`p'==1
		by subjid_`p' outbounds_`p': gen totout_`p'=_N if outbounds_`p'==1
		*Exclude pairs
		bysort subjid_`p': gen tot_`p'=_N if subjid_`p'!=""
		replace exc_`p'=222 if outn_`p'==1 & outbounds_`p'==1 & tot_`p'>2 & exc_`p'==0
		replace subjid_`p'=subjid_str
		replace subjid_`p'="" if exc_`p'!=0
		
		**9E**
		*See if any more rounds required
		*One subject will always have at least two outbounds values
		*See if any have more than 2
		replace totout_`p'=0 if outbounds_`p'==0 | subjid_`p'==""
		egen max_totout_`p'=max(totout_`p')
		local i=(max_totout_`p'[1]-2)
		*Drop remaining 
		drop *diff* outbounds* median* absd* outn_* *totout_`p' *prev hinext* hiprev* tot_ csubjid* cexc*
		*Runs again if remaining outbounds values
	}


}

*Clean up 
foreach p in wt ht hc {
	replace subjid_`p'=subjid_str
	replace subjid_`p'="" if exc_`p'!=0
}

**9F**
*Redo temporary SDEs - use local macro dp instead of p to distinguish it
*Calculate medians for all parameters
foreach dp in wt ht hc {
	duplicates tag subjid_`dp' param agedays, gen(dup_`dp')
	bysort subjid_`dp': egen median_`dp'z_i=median(tbc`dp'z) if exc_`dp'==0 & dup_`dp'==0
	bysort subjid_`dp': egen median_`dp'z=min(median_`dp'z_i)
}


foreach dp in wt ht hc {
	replace whichvar="`dp'"
	if whichvar[1]=="wt" {
		local o ht
	}
	if whichvar[1]=="ht" {
		local o wt
	}
	if whichvar[1]=="hc" {
		local o ht
	}
	gen absd_median_`dp'z=abs(tbc`dp'z-median_`dp'z)
	gen absd_median_`o'z=abs(tbc`dp'z-median_`o'z)
	sort subjid_`dp' agedays absd_median_`dp'z absd_median_`o'z
	by subjid_`dp' agedays: gen d_sort_`dp'_n=_n if exc_`dp'==0
	replace exc_`dp'=2 if exc_`dp'==0 & d_sort_`dp'_n>1
	replace subjid_`dp'=subjid_str
	replace subjid_`dp'="" if exc_`dp'>0
	drop absd_median_*  d_sort*
}
drop dup* *median* max* min* sn* 
				

lab def exc 222 "Evil Twins", add 


save temp_et, replace


foreach f in cf {
	capture confirm file "temp_`f'.dta"
	if !_rc {  
		rm temp_`f'.dta
	}
}

*/

**************************************************************************************
*Step 11: Extreme EWMA
**************************************************************************************


use temp_et, clear

foreach p in wt ht hc {

	 
	**11A**
	*Run as long as there are any meeting EWMA criteria
		*Regular setup
		replace whichvar="`p'"
		local a 5
		local i 1 		
		while `i'>0 {
		
		*Corrz setup
		gen cexc_`p'=exc_`p'
		replace cexc_`p'=20 if ctbc`p'z==.
		gen csubjid_`p'=subjid_`p'
		replace csubjid_`p'="" if ctbc`p'z==.
		
		*Determine exponent based on the max age difference between value and the adjacent values
		sort subjid param agedays
		by subjid param: gen maxagediff_`p'=abs(agedays-agedays[_n-1])
		by subjid param: replace maxagediff_`p'=abs(agedays-agedays[_n+1]) if abs(agedays-agedays[_n+1])>abs(agedays-agedays[_n-1]) & abs(agedays-agedays[_n+1])!=.
		gen expon_`p'=-1.5
		replace expon_`p'=-2.5 if maxagediff_`p'>365.25 & maxagediff_`p'!=.
		replace expon_`p'=-3.5 if maxagediff_`p'>730.5 & maxagediff_`p'!=.
			
		**11B** (in middle of 11A)
		*Calc dEWMA for corrz first
		*Calculate EWMAs/dEWMAs including before/after
		gen str3 ba=""
		gen seldup_`p'=1 if cexc_`p'==0
		replace cexc_`p'=0 if cexc_`p'==2
		replace csubjid_`p'=subjid_str if cexc_`p'==0
		sort csubjid_`p' agedays
		by csubjid_`p': gen tot_`p'=_N
		by csubjid_`p': gen sn_`p'=_n
		summ tot_`p' if cexc_`p'==0
		local m=r(max)
		forvalues n=1/`m' {
			  by csubjid_`p': gen temp_ctbc`p'z_`n'=ctbc`p'z[`n'] if seldup_`p'[`n']==1
			  by csubjid_`p': gen temp_age_`p'_`n'=agedays[`n'] if seldup_`p'[`n']==1
			  gen double temp_ewt_`p'_`n'=abs(`a'+agedays-temp_age_`p'_`n')^expon_`p' if (agedays-temp_age_`p'_`n')!=0 
			  gen double temp_ewtz_`p'_`n'=temp_ctbc`p'z_`n'*temp_ewt_`p'_`n' 
		}
		egen double sumewt_`p'=rowtotal(temp_ewt_`p'_*) if cexc_`p'==0
		egen double sumewtz_`p'=rowtotal(temp_ewtz_`p'_*) if cexc_`p'==0
		gen ewma_`p'=sumewtz_`p'/sumewt_`p' if csubjid_`p'!=""
		gen dewma_`p'=ctbc`p'z-ewma_`p'
		drop temp* sumewt* tot* sn* ba seldup* ewma*
		rename dewma* cdewma*
			
		*Calculate EWMAs/dEWMAs including before/after
		gen str3 ba=""
		gen seldup_`p'=1 if exc_`p'==0
		replace exc_`p'=0 if exc_`p'==2
		replace subjid_`p'=subjid_str if exc_`p'==0
		sort subjid_`p' agedays
		by subjid_`p': gen tot_`p'=_N
		by subjid_`p': gen sn_`p'=_n
		summ tot_`p' if exc_`p'==0
		local m=r(max)
		forvalues n=1/`m' {
			  by subjid_`p': gen temp_tbc`p'z_`n'=tbc`p'z[`n'] if seldup_`p'[`n']==1
			  by subjid_`p': gen temp_age_`p'_`n'=agedays[`n'] if seldup_`p'[`n']==1
			  gen double temp_ewt_`p'_`n'=abs(`a'+agedays-temp_age_`p'_`n')^expon_`p' if (agedays-temp_age_`p'_`n')!=0 
			  gen double temp_ewtz_`p'_`n'=temp_tbc`p'z_`n'*temp_ewt_`p'_`n' 
		}
		egen double sumewt_`p'=rowtotal(temp_ewt_`p'_*) if exc_`p'==0
		egen double sumewtz_`p'=rowtotal(temp_ewtz_`p'_*) if exc_`p'==0
		gen ewma_`p'=sumewtz_`p'/sumewt_`p' if subjid_`p'!=""
		gen dewma_`p'=tbc`p'z-ewma_`p'
		*Before & After
		foreach t in bef aft {
			  replace ba="`t'"
			  gen double sumewt_`p'_`t'=sumewt_`p'
			  gen double sumewtz_`p'_`t'=sumewtz_`p'
			  gen ewma_`p'_`t'=ewma_`p'
			  gen dewma_`p'_`t'=dewma_`p'
			  forvalues d=1/`m' {
					if ba=="bef" {
							local s=`d'+1
					}
					if ba=="aft" {
						local s=`d'-1
					}
					replace sumewt_`p'_`t'=sumewt_`p'-temp_ewt_`p'_`d' if temp_ewt_`p'_`d'!=. & sn_`p'==`s'
					replace sumewtz_`p'_`t'=sumewtz_`p'-temp_ewtz_`p'_`d' if temp_ewtz_`p'_`d'!=. & sn_`p'==`s'
					replace ewma_`p'_`t'=sumewtz_`p'_`t'/sumewt_`p'_`t' if subjid_`p'!="" & sn_`p'==`s'
					replace dewma_`p'_`t'=tbc`p'z-ewma_`p'_`t' if subjid_`p'!="" & sn_`p'==`s'
			  }
		}
		
		*Clean Up
		drop temp*

		**11C**
		*Calculate next/previous z and absolute value of tbcsd of interest
		*Then calc abssum = |z + dewma| - high if value deviant from norms AND individual value, if abnl for pop but nl for kid will be smaller
		sort subjid_`p' agedays
		by subjid_`p': gen d_nextz_`p'=tbc`p'z-tbc`p'z[_n+1]
		by subjid_`p': gen d_prevz_`p'=tbc`p'z-tbc`p'z[_n-1]
		gen abs_tbc`p'z=abs(tbc`p'z)
		bysort subjid_`p': egen max_abs_tbc`p'z=max(abs_tbc`p'z)
		gen abssum_sd_dewma_`p'=abs(tbc`p'z+dewma_`p')

		**11D**
		*Potential exclusions
		*Identify values with |dewma|>3.5 (and |cdewma|>3.5 when relevant)
		gen extreme_exc_`p'=0
		replace extreme_exc_`p'=1 if dewma_`p'>3.5 & (dewma_`p'_bef>3 | dewma_`p'_bef==.) & (dewma_`p'_aft>3|dewma_`p'_aft==.) & (cdewma_`p'>3.5|cdewma_`p'==.) ///
			& tbc`p'z>3.5 & exc_`p'==0 & tot_`p'>2
		replace extreme_exc_`p'=1 if dewma_`p'<-3.5 & (dewma_`p'_bef<-3|dewma_`p'==.) & (dewma_`p'_aft<-3|dewma_`p'_aft==.) & (cdewma_`p'>-3.5|cdewma_`p'==.) /// 
			& tbc`p'z<-3.5 & exc_`p'==0 & tot_`p'>2

		**11E**
		*Of those meeting EWMA criteria, select value with highest abssum and exclude it
		*Do not exclude first/last measurements
		gen firstlast_`p'=sn_`p'==1|sn_`p'==tot_`p'
		gsort subjid_`p' extreme_exc_`p' firstlast_`p' -abssum_sd_dewma_`p'
		by subjid_`p' extreme_exc_`p' firstlast_`p': gen ext_exc_num_`p'=_n if extreme_exc_`p'==1 & firstlast_`p'==0
		replace ext_exc_num_`p'=0 if ext_exc_num_`p'==.
		egen max_ext_exc_num_`p'=max(ext_exc_num_`p')
		replace exc_`p'=333 if ext_exc_num_`p'==1
		replace subjid_`p'=subjid_str
		replace subjid_`p'="" if exc_`p'!=0


		**11G** (Note trigger is calculated before redo temp SDE, but temp SDEs are done before rerunning) 
		*Trigger to run again if additional values meeting EWMA criteria
		local i=max_ext_exc_num_`p'[1]
		drop maxagediff_`p'-max_ext_exc_num_`p'
		*Drop variables
		drop cexc* csubjid*

		**11F**
		*Redo temporary SDEs - use local macro dp instead of p to distinguish it
		*Calculate medians for all parameters
		foreach dp in wt ht hc {
			duplicates tag subjid_`dp' param agedays, gen(dup_`dp')
			bysort subjid_`dp': egen median_`dp'z_i=median(tbc`dp'z) if exc_`dp'==0 & dup_`dp'==0
			bysort subjid_`dp': egen median_`dp'z=min(median_`dp'z_i)
		}
		
		*Determine differences from medians
		local oi="1"
		foreach dp in wt ht hc {
			gen temp`dp'=`oi'
			local oi=`oi'+1
		}
		
		foreach dp in wt ht hc {
			local o empty
			if temp`dp'==1 {
				local o ht	
			}
			if temp`dp'==2 {
				local o wt
			}
			if temp`dp'==3 {
				local o ht
			}
			gen absd_median_`dp'z=abs(tbc`dp'z-median_`p'z)
			gen absd_median_`o'z=abs(tbc`dp'z-median_`o'z)
			sort subjid_`dp' agedays absd_median_`dp'z absd_median_`o'z
			by subjid_`dp' agedays: gen d_sort_`dp'_n=_n if exc_`dp'==0
			replace exc_`dp'=2 if exc_`dp'==0 & d_sort_`dp'_n>1
			replace subjid_`dp'=subjid_str
			replace subjid_`dp'="" if exc_`dp'>0
			drop absd_median_*  d_sort*
		}
		drop dup_* median*z* temp*

	}
}
	
lab def exc 333 "EWMA1 high sum", add


save temp_ewma1, replace



foreach f in biv {
	capture confirm file "temp_`f'.dta"
	if !_rc {  
		rm temp_`f'.dta
	}
}

*/

**************************************************************************************
*Step 13: SDEs
**************************************************************************************



use temp_ewma1, clear

******************************      Set-up      ****************************

*Identify even-numbered days for later
gen evenday=mod(agedays, 2)==0

******************************       13A       ****************************

*Remove identical SDEs first
foreach p in wt ht hc {
	replace exc_`p'=0 if exc_`p'==2
	replace subjid_`p'=subjid_str if exc_`p'==0
	egen identical_`p' = group(subjid_`p' param agedays measurement)
	bysort identical_`p': gen gn_`p'=_n
	replace exc_`p'=770 if exc_`p'==0 & gn>1
	replace subjid_`p'="" if exc_`p'!=0
	drop identical* gn_*
}

******************************       13B       ****************************

*Identify left over duplicate groups
foreach p in wt ht hc {
	duplicates tag subjid_`p' agedays, gen(dup_`p')
	replace dup_`p'=0 if subjid_`p'==""
	bysort subjid_`p' agedays: gen dgn_`p'=_n
	replace dgn_`p'=. if subjid_`p'=="" | dup_`p'==0
	gen dgn_`p'1=dgn_`p'==1
	*Identify days with measurements
	bysort subjid_`p' agedays: gen san_`p'=_n
	gen san_`p'1=san_`p'==1
	*Total measurements on a day
	bysort subjid_`p' agedays: gen satot_`p'=_N
	replace satot_`p'=. if subjid_`p'==""
}	

*Identify same-day similars
foreach p in wt ht hc {
	bysort subjid_`p' agedays: egen min_dup_`p'=min(`p')
	bysort subjid_`p' agedays: egen max_dup_`p'=max(`p')
	gen similar_`p'=0
}	


local oi="1"
foreach dp in wt ht hc {
	gen temp`dp'=`oi'
	local oi=`oi'+1
}
		
foreach p in wt ht hc {
	local o empty
	if temp`p'==1 {
		replace similar_`p'=1 if (max_dup_`p'/min_dup_`p')<1.03 & (max_dup_`p'-min_dup_`p')<2.5
		list subjid_wt exc_wt agedays wt max_dup_wt min_dup_wt similar_wt if subjid==13001009
	}
	if temp`p'==2 {
		replace similar_`p'=1 if (max_dup_`p'-min_dup_`p')<5.081 & min_dup_`p'>=127
		replace similar_`p'=1 if (max_dup_`p'-min_dup_`p')<2.541 & min_dup_`p'<127
		list subjid_ht exc_ht agedays ht max_dup_ht min_dup_ht similar_ht if subjid==13001008
	}
	if temp`p'==3 {
		replace similar_`p'=1 if max_dup_`p'-min_dup_`p'<1.271 & dup_`p'>0
	}
}

drop temp*

******************************       13C       ****************************

*Identify groups of duplicates for which all should be excluded
foreach p in wt ht hc {	  
	*Ratio criteria - add epsilon
		bysort subjid_`p': egen subjdups_`p'=sum(dgn_`p'1)
		bysort subjid_`p': egen subjdays_`p'=sum(san_`p'1)
		replace subjdups_`p'=. if subjid_`p'==""
		replace subjdays_`p'=. if subjid_`p'==""
		gen ratio_`p'_i=subjdups_`p'/subjdays_`p'
		gen crit_ratio_`p'_i=ratio_`p'>(1/3) + 0.0001
		bysort subjid_`p': egen crit_ratio_`p'=min(crit_ratio_`p')

	*Adjacent
		sort subjid_`p' san_`p' agedays
		gen adjacent_`p'_i=1 if san_`p'==1 & dup_`p'>0 & dup_`p'[_n+1]>0 & agedays!=agedays[_n+1] & san_`p'==san_`p'[_n+1]
		replace adjacent_`p'_i=1 if san_`p'==1 & dup_`p'>0 & dup_`p'[_n-1]>0 & agedays!=agedays[_n-1] & san_`p'==san_`p'[_n-1] 
		bysort subjid_`p' agedays: egen adjacent_`p'=max(adjacent_`p'_i)
					  
	*Adjacent and ratio criteria
		gen crit_adjratio_`p'=adjacent_`p'==1 & ratio_`p'>(1/4)


	*Adjacent and first or last criteria
		*Need subjid variable specifically for first of any ageday (dup or not)
		*Variable for which ageday of subject/parameter - forward and reverse
		gen subjid1_`p'=subjid_`p' if san_`p'==1
		sort subjid1_`p' agedays
		by subjid1_`p': gen agedaynum1_`p'=_n
		replace agedaynum1_`p'=. if subjid1_`p'==""
		gsort subjid1_`p' -agedays
		by subjid1_`p': gen ragedaynum1_`p'=_n
		replace ragedaynum1_`p'=. if subjid1_`p'==""

		*Identify max number of first values - only need to do up to 1/3 of that value because otherwise will be excluded by ratio
		summ subjdays_`p', d
		local m=r(max)/3 + 1
		
		*Identify first adjacents 
		sort subjid1_`p' agedaynum1_`p'
		gen first_adj_`p'_i=0
		by subjid1_`p': gen first_adj_`p'_1=adjacent_`p'[1]==1 & adjacent_`p'[2] == 1
		replace first_adj_`p'_i=1 if first_adj_`p'_1==1 & agedaynum1_`p'==1
		forvalues n=2/`m' {
				sort subjid1_`p' agedays
				by subjid1_`p': gen first_adj_`p'_`n'=adjacent_`p'[`n']==1 & adjacent_`p'[`n'-1]==1 & first_adj_`p'_1==1
				replace first_adj_`p'_i=1 if first_adj_`p'_`n'==1 & agedaynum1_`p'==`n'
		}
		bysort subjid_`p' agedays: egen first_adj_`p'=max(first_adj_`p'_i)
		
		*Identify last adjacents 
		sort subjid1_`p' ragedaynum1_`p'
		gen last_adj_`p'_i=0
		by subjid1_`p': gen last_adj_`p'_1=adjacent_`p'[1]==1 & adjacent_`p'[2] == 1
		replace last_adj_`p'_i=1 if last_adj_`p'_1==1 & ragedaynum1_`p'==1
		forvalues n=2/`m' {
				gsort subjid1_`p' -agedays
				by subjid1_`p': gen last_adj_`p'_`n'=adjacent_`p'[`n']==1 & adjacent_`p'[`n'-1]==1 & last_adj_`p'_1==1
				replace last_adj_`p'_i=1 if last_adj_`p'_`n'==1 & ragedaynum1_`p'==`n'
		}
		bysort subjid_`p' agedays: egen last_adj_`p'=max(last_adj_`p'_i)
		
		*Combine first and last adjacents
		gen crit_adjfirstlast_`p'=first_adj_`p'==1 | last_adj_`p'==1

	  
	*Adjacent 3 criteria
		sort subjid_`p' agedays
		gen crit_adjacent3_`p'=1 if dup_`p'>0 & dup_`p'[_n+1]>0 & dup_`p'[_n+2]>0 & agedays==agedays[_n+1] & agedays==agedays[_n+2]
		replace crit_adjacent3_`p'=1 if dup_`p'>0 & dup_`p'[_n-1]>0 & dup_`p'[_n-2]>0 & agedays==agedays[_n-1] & agedays==agedays[_n-2] 
		replace crit_adjacent3_`p'=1 if dup_`p'>0 & dup_`p'[_n-1]>0 & dup_`p'[_n+1]>0 & agedays==agedays[_n+1] & agedays==agedays[_n-1] 
		*list subjid agedays wt dup_wt similar_wt crit_adjacent3_wt if subjid==13001014
	  
	*Exclude if any of above are positive
		replace exc_`p'=771 if crit_ratio_`p'==1 & similar_`p'==0 & exc_`p'==0
		replace exc_`p'=771 if crit_adjfirstlast_`p'==1 & similar_`p'==0 & exc_`p'==0
		*list subjid agedays wt adjacent_wt last* crit_adjfirstlast_wt exc_wt if subjid==13001009
		replace exc_`p'=771 if crit_adjratio_`p'==1 & similar_`p'==0 & exc_`p'==0
		replace exc_`p'=771 if crit_adjacent3_`p'==1 & similar_`p'==0 & exc_`p'==0
		replace subjid_`p'=subjid_str
		replace subjid_`p'="" if exc_`p'!=0 & similar_`p'==0
		*list subjid agedays hc crit*hc exc_hc if subjid==13001010

	  
	*Drop some variables
		drop adjacent*`p' subjdup*`p' subjday*`p' first*`p' last*`p'
}

******************************       13D       ****************************

*Set up for EWMA step

*Update dup variables and identify total number of duplicate groups per subject
foreach p in wt ht hc {	  
	replace dup_`p'=0 if subjid_`p'==""
	drop dgn_`p'*
	bysort subjid_`p' agedays: gen dgn_`p'=_n		
	bysort subjid_`p' agedays: gen dgtot_`p'=_N
	replace dgn_`p'=. if subjid_`p'=="" | dup_`p'==0
	gen dgn_`p'1=dgn_`p'==1
	bysort subjid_`p': egen subjdups_`p'=sum(dgn_`p'1)
	bysort subjid_`p': egen subjdays_`p'=sum(san_`p'1)
	gen nondups_`p'=subjdays_`p'>subjdups_`p' & dup_`p'>0
}

*Generate special exclusion variable that includes non-SDEs and similar SDEs but not non-similar SDEs
foreach p in wt ht hc {
	gen excs_`p'=exc_`p'
	replace excs_`p'=99999 if exc_`p'==0 & dup_`p'>0 & similar_`p'==0
	gen subjids_`p'=subjid_`p'
	replace subjids_`p'="" if excs_`p'!=0
}

******************************       13E       ****************************

*Use EWMA to exclude all but one value from remaining groups of duplicates for which there is at least one other value for the subject
*Include values if they are non SDEs or similar SDEs - do it by making ageweight zero




local a 5
gen str3 ba=""
foreach p in wt ht hc {
	
	*Determine exponent based on min age difference between value and next value
	sort subjid param agedays
	by subjid param: gen maxagediff_`p'=abs(agedays-agedays[_n-1])
	by subjid param: replace maxagediff_`p'=abs(agedays-agedays[_n+1]) if abs(agedays-agedays[_n+1])>abs(agedays-agedays[_n-1]) & abs(agedays-agedays[_n+1])!=.
	gen expon_`p'=-1.5
	replace expon_`p'=-2.5 if maxagediff_`p'>365.25 & maxagediff_`p'!=.
	replace expon_`p'=-3.5 if maxagediff_`p'>730.5 & maxagediff_`p'!=.

	bysort subjid_`p' agedays: gen `p'_sort_n=_n if exc_`p'==0
	bysort subjid_`p': gen tot_`p'=_N
	bysort subjid_`p': gen sn_`p'=_n
	summ tot_`p' if exc_`p'==0
	local m=r(max)
	di "`m'"
	forvalues n=1/`m' {
		by subjid_`p': gen temp_tbc`p'z_`n'=tbc`p'z[`n'] if `p'_sort_n[`n']==1
		by subjid_`p': gen temp_age_`p'_`n'=agedays[`n'] if `p'_sort_n[`n']==1
		gen double temp_ewt_`p'_`n'=abs(`a'+agedays-temp_age_`p'_`n')^expon_`p' if (agedays-temp_age_`p'_`n')!=0 
		by subjid_`p': gen temp_excs_`p'_`n'=excs_`p'[`n']
		replace temp_ewt_`p'_`n'=0 if temp_excs_`p'_`n'!=0
		gen double temp_ewtz_`p'_`n'=temp_tbc`p'z_`n'*temp_ewt_`p'_`n'
		}
	egen double sumewt_`p'=rowtotal(temp_ewt_`p'_*) 
	egen double sumewtz_`p'=rowtotal(temp_ewtz_`p'_*)
	gen ewma_`p'=sumewtz_`p'/sumewt_`p' if subjid_`p'!=""
	gen dewma_`p'=tbc`p'z-ewma_`p'
	gen absdewma_`p'=abs(dewma_`p')
	sort subjid_`p' agedays absdewma_`p'
	by subjid_`p' agedays: gen order_`p'=_n
	*Exclude all values if lowest absdewma is >1
	replace exc_`p'=772 if exc_`p'==0 & dup_`p'>0 & subjdays_`p'>1 & absdewma_`p'>1
	*list subjid agedays wt exc_wt  if subjid==13001009
	*Exclude all but lowest absdewma if there is at least one other value for subject
	replace exc_`p'=773 if exc_`p'==0 & dup_`p'>0 & subjdays_`p'>1 & order_`p'>1
	*list subjid agedays hc shcz absdewma_hc dup_hc subjdays_hc order_hc excs_hc exc_hc if subjid==13001011
	replace subjid_`p'=subjid_str
	replace subjid_`p'="" if exc_`p'!=0
	drop order_`p' excs_`p' subjids_`p' expon_`p'
}

******************************       13F       ****************************


*Use median of other parameter to exclude all but one value from groups of duplicates for which there are no other SDEs for subject
	
*Generate medians for use if only SDEs available for parameter (will be similar)
*Using dp instead of p just to match other similar code that has to distinguish between p and dp
*Determine differences from medians
local oi="1"
foreach p in wt ht hc {
	bysort subjid_`p': egen median_`p'z_i=median(tbc`p'z) if exc_`p'==0
	bysort subjid_`p': egen median_`p'z=min(median_`p'z) 
	gen temp`p'=`oi'
	local oi=`oi'+1
}
	
foreach p in wt ht hc {
	local o empty
	if temp`p'==1 {
		local o ht	
	}
	if temp`p'==2 {
		local o wt
	}
	if temp`p'==3 {
		local o ht
	}
	gen absdiff_median_`p'z=abs(tbc`p'z-median_`o'z)
	sort subjid_`p' agedays absdiff_median_`p'z
	by subjid_`p' agedays: gen order_`p'=_n
	*Exclude all if smallest absolute difference from other median is >2
	replace exc_`p'=772 if exc_`p'==0 & dup_`p'>0 & nondups_`p'==0 & median_`o'z!=. & absdiff_median_`p'z>2
	*For remaining, exclude all but smallest difference
	replace exc_`p'=774 if exc_`p'==0 & dup_`p'>0 & nondups_`p'==0 & median_`o'z!=. & order_`p'>1
	replace subjid_`p'=subjid_str
	replace subjid_`p'="" if exc_`p'!=0
	drop order_`p' temp`p'
	
******************************       13G       ****************************

	*If only one day with SDEs, keep one closest to median of POI if there are 3 or more in SDE group
	*If there are an even number of values in the SDE group and the two middle values are different, they will be equidistant from the median; in that case regain the higher value on an even numbered day and exclude the lower value on an odd numbered day
	bysort subjid_`p': egen median_only`p'z_i=median(tbc`p'z) if exc_`p'==0
	bysort subjid_`p': egen median_only`p'z=min(median_`p'z) 
	gen absdiff_median_only`p'z=abs(tbc`p'z-median_only`p'z)
	*Determine if the lowest absdiff has another equal lowest absdiff
	bysort subjid_`p' agedays: egen mindiff_`p'=min(absdiff_median_only`p'z)
	gen ismindiff_`p'=absdiff_median_only`p'z==mindiff_`p'
	bysort subjid_`p' agedays: egen nummindiff_`p'=sum(ismindiff_`p')
	*If there are not two equal mindiffs
	replace exc_`p'=774 if exc_`p'==0 & dup_`p'>0 & subjdays_`p'==1 & median_`o'z==. & dgtot_`p'>2 & nummindiff_`p'==1 & ismindiff_`p'!=1
	*If there are two equal mindiffs
	sort subjid_`p' agedays absdiff_median_only`p'z measurement
	by subjid_`p' agedays: gen order_`p'=_n
	*Even days - keep lower (first) value (exclude higher/second)
	replace exc_`p'=774 if exc_`p'==0 & dup_`p'>0 & subjdays_`p'==1 & median_`o'z==. & dgtot_`p'>2 & nummindiff_`p'==2 & evenday==1 & order_`p'>1
	*Odd days - keep higher (second) value (exclude lower/first)
	replace exc_`p'=774 if exc_`p'==0 & dup_`p'>0 & subjdays_`p'==1 & median_`o'z==. & dgtot_`p'>2 & nummindiff_`p'==2 & evenday==0 & order_`p'!=2
	replace subjid_`p'=subjid_str
	replace subjid_`p'="" if exc_`p'!=0	
	drop order_`p'
	
******************************       13H       ****************************

	*If only values are SDEs and there are only two on a day, include lower/exclude higher on even-numbered day and include higher/exclude lower value on odd-numbered day
	sort subjid_`p' agedays measurement
	by subjid_`p' agedays: gen order_`p'=_n
	*Even days - keep lower (first) value (exclude higher/second)
	replace exc_`p'=774 if exc_`p'==0 & dup_`p'>0 & nondups_`p'==0 & median_`o'z==. & dgtot_`p'==2 & evenday==1 & order_`p'==2
	*Odd days - keep higher (second) value (exclude lower/first)
	replace exc_`p'=774 if exc_`p'==0 & dup_`p'>0 & nondups_`p'==0 & median_`o'z==. & dgtot_`p'==2 & evenday==0 & order_`p'==1
	replace subjid_`p'=subjid_str
	replace subjid_`p'="" if exc_`p'!=0
}

******************************    Clean-up    ****************************

drop evenday-order_hc


*Do NOT replace exclusion code if previously excluded as extremes as was done before -- want to still identify those as extremes

lab def exc 770 "SDE-Identical", modify
lab def exc 771 "SDE-All-Exclude", modify
lab def exc 772 "SDE-Extreme", modify
lab def exc 773 "SDE-EWMA", modify
lab def exc 774 "SDE-One-Day", modify

save temp_sde, replace



foreach f in et {
	capture confirm file "temp_`f'.dta"
	if !_rc {  
		rm temp_`f'.dta
	}
}

*/

**************************************************************************************
*Step 15: EWMA Moderate 
**************************************************************************************


use temp_sde, clear

**1A**
*Plus/minus measurements

gen wt_plus=wt+0.05*wt
gen wt_minus=wt-0.05*wt
gen ht_plus=ht+1
gen ht_minus=ht-1
gen hc_plus=hc+1
gen hc_minus=hc-1

lab var wt_plus "Weight  +  5%"
lab var wt_minus "Weight - 5%"
lab var ht_plus "Height + 1 cm"
lab var ht_minus "Height - 1 cm"
lab var hc_plus "HC + 0.5 cm"
lab var hc_minus "HC - 0.5 cm"

**15B**
*Make Z scores for plus minus
foreach c in who cdc {
	foreach d in plus minus {
		foreach p in wt ht hc{
			gen `c'`p'z_`d'=(`p'_`d' - orig_`c'_`p'_m) / orig_`c'_`p'_csd_pos
			replace `c'`p'z_`d'=(`p'_`d' - orig_`c'_`p'_m) / orig_`c'_`p'_csd_neg if `p' < orig_`c'_`p'_m
		}
	}
}

*Smooth and recenter
replace whoweight=4-ageyears
replace cdcweight=ageyears-2
foreach d in plus minus {
	foreach p in wt ht hc{
		gen s`p'z_`d'=(who`p'z_`d' * whoweight + cdc`p'z_`d' * cdcweight)/2
		replace s`p'z_`d'=who`p'z_`d' if ageyears<=2
		replace s`p'z_`d'=cdc`p'z_`d' if ageyears>=4
		gen tbc`p'z_`d'=s`p'z_`d'-rc`p'z
	}
}

lab var tbcwtz_plus "Z score for weight + 5%"
lab var tbcwtz_minus "Z score for weight + 5%"
lab var tbchtz_plus "Z score for height + 1 cm"
lab var tbchtz_minus "Z score for height - 1 cm"
lab var tbchcz_plus "Z score for HC + 1 cm"
lab var tbchcz_minus "Z score for HC - 1 cm"

**15C**
*Temporarily exclude birth for HT and HC only
foreach p in ht hc {
	replace exc_`p'=5000 if agedays==0 & exc_`p'==0
	replace subjid_`p'="" if exc_`p'==5000
}


**15D**
*EWMA
local a 5
foreach p in wt ht hc {
	local i 1
	local run 1
	while `i'>0 {	

		*Corrz setup
		gen cexc_`p'=exc_`p'
		replace cexc_`p'=20 if ctbc`p'z==.
		gen csubjid_`p'=subjid_`p'
		replace csubjid_`p'="" if ctbc`p'z==.
		
		*Regular setup
		replace whichvar="`p'"
		local e -1.5
		local a 5
		local i 1 
      *Run as long as there are any meeting EWMA criteria
		while `i'>0 {
			*Determine exponent based on min age difference between value and next value
			sort subjid param agedays
			by subjid param: gen maxagediff_`p'=abs(agedays-agedays[_n-1])
			by subjid param: replace maxagediff_`p'=abs(agedays-agedays[_n+1]) if abs(agedays-agedays[_n+1])>abs(agedays-agedays[_n-1]) & abs(agedays-agedays[_n+1])!=.
			gen expon_`p'=-1.5
			replace expon_`p'=-2.5 if maxagediff_`p'>365.25 & maxagediff_`p'!=.
			replace expon_`p'=-3.5 if maxagediff_`p'>730.5 & maxagediff_`p'!=.
			
			*Calc dEWMA for corrz first
			*Calculate EWMAs/dEWMAs including before/after
			gen str3 ba=""
			gen seldup_`p'=1 if cexc_`p'==0
			replace cexc_`p'=0 if cexc_`p'==2
			replace csubjid_`p'=subjid_str if cexc_`p'==0
			sort csubjid_`p' agedays
			by csubjid_`p': gen tot_`p'=_N
			by csubjid_`p': gen sn_`p'=_n
			summ tot_`p' if cexc_`p'==0
			local m=r(max)
			forvalues n=1/`m' {
				  by csubjid_`p': gen temp_ctbc`p'z_`n'=ctbc`p'z[`n'] if seldup_`p'[`n']==1
                  by csubjid_`p': gen temp_age_`p'_`n'=agedays[`n'] if seldup_`p'[`n']==1
                  gen double temp_ewt_`p'_`n'=abs(`a'+agedays-temp_age_`p'_`n')^expon_`p' if (agedays-temp_age_`p'_`n')!=0 
                  gen double temp_ewtz_`p'_`n'=temp_ctbc`p'z_`n'*temp_ewt_`p'_`n' 
			}
            egen double sumewt_`p'=rowtotal(temp_ewt_`p'_*) if cexc_`p'==0
            egen double sumewtz_`p'=rowtotal(temp_ewtz_`p'_*) if cexc_`p'==0
            gen ewma_`p'=sumewtz_`p'/sumewt_`p' if csubjid_`p'!=""
            gen dewma_`p'=ctbc`p'z-ewma_`p'
			drop temp* sumewt* tot* sn* ba seldup* ewma*
			rename dewma* cdewma*
		
		gen str3 ba=""
		sort subjid_`p' agedays
		by subjid_`p': gen tot_`p'=_N
		by subjid_`p': gen sn_`p'=_n
		summ tot_`p' if exc_`p'==0
		local m=r(max)
		forvalues n=1/`m' {
			  by subjid_`p': gen temp_tbc`p'z_`n'=tbc`p'z[`n'] 
			  by subjid_`p': gen temp_age_`p'_`n'=agedays[`n'] 
			  gen double temp_ewt_`p'_`n'=abs(`a'+agedays-temp_age_`p'_`n')^expon_`p' if sn_`p'!=`n'
			  gen double temp_ewtz_`p'_`n'=temp_tbc`p'z_`n'*temp_ewt_`p'_`n' 
			  }
		egen double sumewt_`p'=rowtotal(temp_ewt_`p'_*) if exc_`p'==0
		egen double sumewtz_`p'=rowtotal(temp_ewtz_`p'_*) if exc_`p'==0
		gen ewma_`p'=sumewtz_`p'/sumewt_`p' if subjid_`p'!=""
		gen dewma_`p'=tbc`p'z-ewma_`p'
		*Before & After
		foreach t in bef aft {
			replace ba="`t'"
			gen double sumewt_`p'_`t'=sumewt_`p'
			gen double sumewtz_`p'_`t'=sumewtz_`p'
			gen ewma_`p'_`t'=ewma_`p'
			gen dewma_`p'_`t'=dewma_`p'
			forvalues d=1/`m' {
				if ba=="bef" {
					  local s=`d'+1
				}
				if ba=="aft" {
					  local s=`d'-1
				}
				replace sumewt_`p'_`t'=sumewt_`p'-temp_ewt_`p'_`d' if temp_ewt_`p'_`d'!=. & sn_`p'==`s'
				replace sumewtz_`p'_`t'=sumewtz_`p'-temp_ewtz_`p'_`d' if temp_ewtz_`p'_`d'!=. & sn_`p'==`s'
				replace ewma_`p'_`t'=sumewtz_`p'_`t'/sumewt_`p'_`t' if subjid_`p'!="" & sn_`p'==`s'
				replace dewma_`p'_`t'=tbc`p'z-ewma_`p'_`t' if subjid_`p'!="" & sn_`p'==`s'
			}
		}     
		drop temp*
		drop sn_`p'
		sort subjid_`p' agedays
		by subjid_`p': gen sn_`p'=_n
		by subjid_`p': gen vistot_`p'=_N
		by subjid_`p': gen d_agedays_prev_`p'=agedays-agedays[_n-1]
		by subjid_`p': gen d_agedays_next_`p'=agedays[_n+1]-agedays

		/*Indicator for whether something is first/final value and how far from prior value
		missing if pair or single
		1 for first value
		2 for last and <2 years after prior
		3 for last and >=2 years after prior
		0 for everything in the middle
		missing for all others
		*/
		gen fl_`p'=0
		replace fl_`p'=1 if sn_`p'==1
		replace fl_`p'=2 if sn_`p'==vistot_`p' * (d_agedays_prev_`p'<(365.25*2))
		replace fl_`p'=3 if sn_`p'==vistot_`p' * (d_agedays_prev_`p'>=(365.25*2))
		replace fl_`p'=. if exc_`p'>0
		replace fl_`p'=. if vistot_`p'<3
		sort subjid_`p' agedays
		
		**15E**
		/*Variables for prior/next and difference between prior/next (current minus prior, current minus next)*/
		by subjid_`p': gen tbc`p'z_prev=tbc`p'z[_n-1]
		by subjid_`p': gen tbc`p'z_next=tbc`p'z[_n+1]
		by subjid_`p': gen d_prevz_`p'=tbc`p'z-tbc`p'z[_n-1]
		by subjid_`p': gen d_nextz_`p'=tbc`p'z-tbc`p'z[_n+1]
		
		/*Variables for difference between plus/minus of current and prior or next*/
		foreach d in plus minus {
			  by subjid_`p': gen d_prevz_`d'_`p'=tbc`p'z_`d'-tbc`p'z[_n-1]
			  by subjid_`p': gen d_nextz_`d'_`p'=tbc`p'z_`d'-tbc`p'z[_n+1]
		}
		gen exc_temp_`p'=.
					
		**15F**
		*Additional criteria for all exclusions - bef/aft>1, difference prior/next>1, and difference plus/minus prior/next>1
		gen addcrithigh=dewma_`p'_bef>1 & dewma_`p'_aft>1 & ///
			((d_nextz_`p'>1 & d_nextz_plus_`p'>1 & d_nextz_minus_`p'>1) | d_nextz_`p'==.) & ///
			((d_prevz_`p'>1 & d_prevz_plus_`p'>1 & d_prevz_minus_`p'>1) | d_prevz_`p'==.)
			  
		gen addcritlow=dewma_`p'_bef<-1 & dewma_`p'_aft<-1 & ///
			((d_nextz_`p'<-1 & d_nextz_plus_`p'<-1 & d_nextz_minus_`p'<-1) | d_nextz_`p'==.) & ///
			((d_prevz_`p'<-1 & d_prevz_plus_`p'<-1 & d_prevz_minus_`p'<-1) | d_prevz_`p'==.)		
		
		**15H** Potential exclusions (G is embedded below)
		
		*For middle values, absdewma>1 
		replace exc_temp_`p'=1008 if fl_`p'==0 & ///
			dewma_`p'>1 & dewma_`p'!=. & (cdewma_`p'>1 | cdewma_`p'==.) & addcrithigh==1
		replace exc_temp_`p'=1008 if fl_`p'==0 & ///
			dewma_`p'<-1 & (cdewma_`p'<-1 | cdewma_`p'==.) & addcritlow==1

		
		*For birth values with next value within one year, absdewma>3
		replace exc_temp_`p'=1018 if agedays==0 & d_agedays_next_`p'<365.25 & ///
			dewma_`p'>3 & dewma_`p'!=. & (cdewma_`p'>3 | cdewma_`p'==.) & addcrithigh==1
		replace exc_temp_`p'=1018 if agedays==0 & d_agedays_next_`p'<365.25 & ///
			dewma_`p'<-3 & (cdewma_`p'<-3 | cdewma_`p'==.) & addcritlow==1			  
		
		*For birth with next value >=1 year, absdewma>4
		replace exc_temp_`p'=1029 if agedays==0 & d_agedays_next_`p'>=365.25 & /// 
			dewma_`p'>4 & dewma_`p'!=. & (cdewma_`p'>4 | cdewma_`p'==.) & addcrithigh==1
		replace exc_temp_`p'=1029 if agedays==0 & d_agedays_next_`p'>=365.25 & ///
			dewma_`p'<-4 & (cdewma_`p'<-4 | cdewma_`p'==.) & addcritlow==1					
		
		*For first non-birth values with next value within one year, absdewma>2
		replace exc_temp_`p'=1028 if fl_`p'==1 & agedays!=0 & d_agedays_next_`p'<365.25 & ///
			dewma_`p'>2 & dewma_`p'!=. & (cdewma_`p'>2 | cdewma_`p'==.) & addcrithigh==1
		replace exc_temp_`p'=1028 if fl_`p'==1 & agedays!=0 & d_agedays_next_`p'<365.25 & ///
			dewma_`p'<-2 & (cdewma_`p'<-2 | cdewma_`p'==.) & addcritlow==1
			  
		*For first non-birth with next value >=1 year, absdewma>3
		replace exc_temp_`p'=1029 if fl_`p'==1 & d_agedays_next_`p'>=365.25 & /// 
			dewma_`p'>3 & dewma_`p'!=. & (cdewma_`p'>3 | cdewma_`p'==.) & addcrithigh==1
		replace exc_temp_`p'=1029 if fl_`p'==1 & d_agedays_next_`p'>=365.25 & ///
			dewma_`p'<-3 & (cdewma_`p'<-3 | cdewma_`p'==.) & addcritlow==1
			  
		*For last values with prior value <2 year and prior |tbcsd|<2, absdewma>2
		replace exc_temp_`p'=1038 if fl_`p'==2 & abs(tbc`p'z_prev)<2 & ///
			dewma_`p'>2 & dewma_`p'!=. & (cdewma_`p'>2 | cdewma_`p'==.) & addcrithigh==1
		replace exc_temp_`p'=1038 if fl_`p'==2 & abs(tbc`p'z_prev)<2 & ///
			dewma_`p'<-2 & (cdewma_`p'<-2 | cdewma_`p'==.) & addcritlow==1
			  
		*For last values with prior value <2 year and prior |tbcsd|>=2, absdewma>tbcsd_prev
		*Do corrz>3
		replace exc_temp_`p'=1039 if fl_`p'==2 & abs(tbc`p'z_prev)>=2 & ///
			dewma_`p'>abs(tbc`p'z_prev) & dewma_`p'!=. & (cdewma_`p'>3 | cdewma_`p'==.) & addcrithigh==1
		replace exc_temp_`p'=1039 if fl_`p'==2 & abs(tbc`p'z_prev)>=2 & ///
			dewma_`p'<(-1 *abs(tbc`p'z_prev)) & (cdewma_`p'<-3 | cdewma_`p'==.) & addcritlow==1
		
		**15G**
		*Determine medians for parameter and other parameters and difference from medians for value of interest
		local oi="1"
		foreach dp in wt ht hc {
			gen temp`dp'=`oi'
			local oi=`oi'+1
			bysort subjid_`dp': egen median_tbc`dp'z_i=median(tbc`dp'z) if exc_`dp'==0
			bysort subjid: egen median_tbc`dp'z=min(median_tbc`dp'z_i) 
		}

		foreach dp in `p' {
			if temp`dp'==1 {
				local o ht	
			}
			
			if temp`dp'==2 {
				local o wt
			}
			
			if temp`dp'==3 {
				local o ht
			}

			gen absd_median_tbc`dp'z=abs(tbc`dp'z-median_tbc`p'z)
			gen absd_median_tbc`o'z=abs(tbc`dp'z-median_tbc`o'z)
		}

		*Also comparison tbcz for other parameter on the same day as measurement of interest
		gen comp_tbc`o'z_i=tbc`o'z if exc_`o'==0
		bysort subjid agedays: egen comp_tbc`o'z=min(comp_tbc`o'z_i)
		
		/*For last values with prior value >=2 year and prior |tbcsd|<2
			-absdewma>3
			-difference between value of interest and other parameter is >4 (or median of other parameter if no value on same day)
		*/
		replace exc_temp_`p'=1048 if fl_`p'==3 & abs(tbc`p'z_prev)<2 & ///
			dewma_`p'>3 & dewma_`p'!=. & (cdewma_`p'>3 | cdewma_`p'==.) & addcrithigh==1 & ///
			(((tbc`p'z-tbc`o'z)>4 & tbc`o'z!=.)|((tbc`p'z-median_tbc`o'z)>4 & median_tbc`o'z!=. & tbc`o'z==.)|median_tbc`o'z==.)
		replace exc_temp_`p'=1048 if fl_`p'==3 & abs(tbc`p'z_prev)<2 & ///
			dewma_`p'<-3 & (cdewma_`p'<-3 | cdewma_`p'==.) & addcritlow==1 & ///
			(((tbc`p'z-tbc`o'z)<-4 & tbc`o'z!=.)|((tbc`p'z-median_tbc`o'z)<-4 & median_tbc`o'z!=. & tbc`o'z==.)|median_tbc`o'z==.)
			  
		/*For last values with prior value >=2 year and prior |tbcsd|>=2
			-absdewma>1+|tbcsd prior value|
			-difference between value of interest and other parameter is >4 (or median of other parameter if no value on same day)
		*/
		replace exc_temp_`p'=1049 if fl_`p'==3 & abs(tbc`p'z_prev)>=2 & ///
			dewma_`p'>(1+abs(tbc`p'z_prev)) & dewma_`p'!=. &(cdewma_`p'>3 | cdewma_`p'==.) & addcrithigh==1 & ///
			(((tbc`p'z-tbc`o'z)>4 & tbc`o'z!=.)|((tbc`p'z-median_tbc`o'z)>4 & median_tbc`o'z!=. & tbc`o'z==.)|median_tbc`o'z==.)
		replace exc_temp_`p'=1049 if fl_`p'==3 & abs(tbc`p'z_prev)>=2 & ///
			dewma_`p'<(-1-abs(tbc`p'z_prev)) & (cdewma_`p'<-3 | cdewma_`p'==.) & addcritlow==1 & ///
			(((tbc`p'z-tbc`o'z)<-4 & tbc`o'z!=.)|((tbc`p'z-median_tbc`o'z)<-4 & median_tbc`o'z!=. & tbc`o'z==.)|median_tbc`o'z==.)
		
		**15I**
		gen abssum_sd_dewma_`p'=abs(tbc`p'z+dewma_`p')
		
		*15J**
		*Identify potential exclusion for each subject/param with highest sum of absolute sd and absolute dewma
		gen exc_ind_`p'=exc_temp_`p'!=.
		gsort subjid_`p' agedays exc_ind_`p' fl_`p' -abssum_sd_dewma_`p'
		by subjid_`p' agedays exc_ind_`p': gen exc_num_`p'=_n if exc_ind_`p'==1
		replace exc_num_`p'=0 if exc_num_`p'==.
		egen max_exc_num_`p'=max(exc_num_`p')
		
		*Exclude identified value, reset other potential exclusions to include, start cycle again if any had >1 potential exclusions
		replace exc_`p'=exc_temp_`p' if exc_num_`p'==1
		replace subjid_`p'=subjid_str
		replace subjid_`p'="" if exc_`p'!=0		
		replace cexc_`p'=exc_`p' if cexc_`p'==0
		replace csubjid_`p'="" if cexc_`p'>0
		local i=max_exc_num_`p'[1]
		di "i at end of run `run'=`i'"
		local run=`run'+1
		*Drop variables
		drop ba-max_exc_num_`p'
		drop cdewma* maxagediff_`p' expon_`p'

    }
}

}

*Clean-up
drop cexc* csubjid* 

foreach p in ht hc {
	replace exc_`p'=0 if exc_`p'==5000
	replace subjid_`p'=subjid_str if exc_`p'==0
}

lab def exc ///
	1008 "EWMA2 middle" ///
	1018 "EWMA2 birth WT" ///
	1019 "EWMA2 birth WT ext" /// 
	1028 "EWMA2 first" ///
	1029 "EWMA2 first ext" /// 
	1038 "EWMA2 last" ///
	1039 "EWMA2 last high" ///
    1048 "EWMA2 last ext" ///
	1049 "EWMA2 last ext high" ///
	, modify


save temp_ewma2, replace



foreach f in ewma1 {
	capture confirm file "temp_`f'.dta"
	if !_rc {  
		rm temp_`f'.dta
	}
}

*/

**************************************************************************************
*Step 16: EWMA Moderate HEIGHT AND HC BIRTH
**************************************************************************************


use temp_ewma2, clear

*16A** Calculate EWMA values
*EWMA - only need to run once
local a 5
foreach p in ht hc {
	
	*Corrz setup
	gen cexc_`p'=exc_`p'
	replace cexc_`p'=20 if ctbc`p'z==.
	gen csubjid_`p'=subjid_`p'
	replace csubjid_`p'="" if ctbc`p'z==.

	*Regular setup
	replace whichvar="`p'"
	local a 5
	local i 1 
	
	*Determine exponent based on min age difference between value and next value
	sort subjid param agedays
	by subjid param: gen maxagediff_`p'=abs(agedays-agedays[_n-1])
	by subjid param: replace maxagediff_`p'=abs(agedays-agedays[_n+1]) if abs(agedays-agedays[_n+1])>abs(agedays-agedays[_n-1]) & abs(agedays-agedays[_n+1])!=.
	gen expon_`p'=-1.5
	replace expon_`p'=-2.5 if maxagediff_`p'>365.25 & maxagediff_`p'!=.
	replace expon_`p'=-3.5 if maxagediff_`p'>730.5 & maxagediff_`p'!=.		
	
	*Calc dEWMA for corrz first
	*Calculate EWMAs/dEWMAs including before/after
	gen str3 ba=""
	gen seldup_`p'=1 if cexc_`p'==0
	replace cexc_`p'=0 if cexc_`p'==2
	replace csubjid_`p'=subjid_str if cexc_`p'==0
	sort csubjid_`p' agedays
	by csubjid_`p': gen tot_`p'=_N
	by csubjid_`p': gen sn_`p'=_n
	summ tot_`p' if cexc_`p'==0
	local m=r(max)
	forvalues n=1/`m' {
		  by csubjid_`p': gen temp_ctbc`p'z_`n'=ctbc`p'z[`n'] if seldup_`p'[`n']==1
		  by csubjid_`p': gen temp_age_`p'_`n'=agedays[`n'] if seldup_`p'[`n']==1
		  gen double temp_ewt_`p'_`n'=abs(`a'+agedays-temp_age_`p'_`n')^expon_`p' if (agedays-temp_age_`p'_`n')!=0 
		  gen double temp_ewtz_`p'_`n'=temp_ctbc`p'z_`n'*temp_ewt_`p'_`n' 
	}
	egen double sumewt_`p'=rowtotal(temp_ewt_`p'_*) if cexc_`p'==0
	egen double sumewtz_`p'=rowtotal(temp_ewtz_`p'_*) if cexc_`p'==0
	gen ewma_`p'=sumewtz_`p'/sumewt_`p' if csubjid_`p'!=""
	gen dewma_`p'=ctbc`p'z-ewma_`p'
	drop temp* sumewt* tot* sn* ba seldup* ewma*
	rename dewma* cdewma*
	
	*Regular EWMA
	gen str3 ba=""
	sort subjid_`p' agedays
	by subjid_`p': gen tot_`p'=_N
	by subjid_`p': gen sn_`p'=_n
	summ tot_`p' if exc_`p'==0
	local m=r(max)
	forvalues n=1/`m' {
		  by subjid_`p': gen temp_tbc`p'z_`n'=tbc`p'z[`n'] 
		  by subjid_`p': gen temp_age_`p'_`n'=agedays[`n'] 
		  gen double temp_ewt_`p'_`n'=abs(`a'+agedays-temp_age_`p'_`n')^expon_`p' if sn_`p'!=`n'
		  gen double temp_ewtz_`p'_`n'=temp_tbc`p'z_`n'*temp_ewt_`p'_`n' 
		  }
	egen double sumewt_`p'=rowtotal(temp_ewt_`p'_*) if exc_`p'==0
	egen double sumewtz_`p'=rowtotal(temp_ewtz_`p'_*) if exc_`p'==0
	gen ewma_`p'=sumewtz_`p'/sumewt_`p' if subjid_`p'!=""
	gen dewma_`p'=tbc`p'z-ewma_`p'
	*Before & After
	foreach t in bef aft {
		replace ba="`t'"
		gen double sumewt_`p'_`t'=sumewt_`p'
		gen double sumewtz_`p'_`t'=sumewtz_`p'
		gen ewma_`p'_`t'=ewma_`p'
		gen dewma_`p'_`t'=dewma_`p'
		forvalues d=1/`m' {
			if ba=="bef" {
				  local s=`d'+1
			}
			if ba=="aft" {
				  local s=`d'-1
			}
			replace sumewt_`p'_`t'=sumewt_`p'-temp_ewt_`p'_`d' if temp_ewt_`p'_`d'!=. & sn_`p'==`s'
			replace sumewtz_`p'_`t'=sumewtz_`p'-temp_ewtz_`p'_`d' if temp_ewtz_`p'_`d'!=. & sn_`p'==`s'
			replace ewma_`p'_`t'=sumewtz_`p'_`t'/sumewt_`p'_`t' if subjid_`p'!="" & sn_`p'==`s'
			replace dewma_`p'_`t'=tbc`p'z-ewma_`p'_`t' if subjid_`p'!="" & sn_`p'==`s'
		}
	}     
	drop temp*
	drop sn_`p'
	sort subjid_`p' agedays
	by subjid_`p': gen sn_`p'=_n
	by subjid_`p': gen vistot_`p'=_N
	by subjid_`p': gen d_agedays_next_`p'=agedays[_n+1]-agedays
	
	**16B**
	/*Variables for next and difference between value and next*/
	by subjid_`p': gen tbc`p'z_next=tbc`p'z[_n+1]
	by subjid_`p': gen d_nextz_`p'=tbc`p'z-tbc`p'z[_n+1]
	
	/*Variables for difference between plus/minus of current and next*/
	foreach d in plus minus {
		  by subjid_`p': gen d_nextz_`d'_`p'=tbc`p'z_`d'-tbc`p'z[_n+1]
	}
	gen exc_temp_`p'=.
			
	**16C**	
	*Additional criteria for all exclusions - bef/aft>1, difference next>1, and difference plus/minus next>1
	gen addcrithigh=dewma_`p'_bef>1 & dewma_`p'_aft>1 & d_nextz_`p'>1 & d_nextz_plus_`p'>1 & d_nextz_minus_`p'>1
		  
	gen addcritlow=dewma_`p'_bef<-1 & dewma_`p'_aft<-1 & d_nextz_`p'<-1 & d_nextz_plus_`p'<-1 & d_nextz_minus_`p'<-1
	
	**16D**
	*Identify exclusions - birth only	
	*For birth values with next value within one year, absdewma>3
	replace exc_`p'=1068 if agedays==0 & d_agedays_next_`p'<365.25 & ///
		dewma_`p'>3 & dewma_`p'!=. & (cdewma_`p'>3 | cdewma_`p'==.) & addcrithigh==1
	replace exc_temp_`p'=1068 if agedays==0 & d_agedays_next_`p'<365.25 & ///
		dewma_`p'<-3 & (cdewma_`p'<-3 | cdewma_`p'==.) & addcritlow==1			  
	
	*For birth with next value >1 year, absdewma>4
	replace exc_`p'=1069 if agedays==0 & d_agedays_next_`p'>365.25 & /// 
		dewma_`p'>4 & dewma_`p'!=. & (cdewma_`p'>4 | cdewma_`p'==.) & addcrithigh==1
	replace exc_`p'=1069 if agedays==0 & d_agedays_next_`p'>365.25 & ///
		dewma_`p'<-4 & (cdewma_`p'<-4 | cdewma_`p'==.) & addcritlow==1			
		
	
	*Clean up
	replace subjid_`p'=subjid_str
	replace subjid_`p'="" if exc_`p'!=0	
	*Drop variables
	drop ba-d_nextz_minus_`p'
	drop cdewma* maxagediff_`p' expon_`p'
	drop cexc* csubjid* addcrit*

}


lab def exc /// 
	1068 "EWMA2 birth HT HC" ///
	1069 "EWMA2 birth HT HC ext" /// 
	, modify

	  
save temp_ewma2b, replace



foreach f in sde {
	capture confirm file "temp_`f'.dta"
	if !_rc {  
		rm temp_`f'.dta
	}
}

*/

**************************************************************************************
*Step 17: Absolute (Raw) Differences
**************************************************************************************


use temp_ewma2b, clear

local a 5
local i 1
foreach p in ht hc {
	replace whichvar="`p'"
	while `i'>0 {
		gen str3 ba=""
		**17A**
		*Identify midpoints and intervals
		sort subjid_`p' agedays
		by subjid_`p': gen d_agedays_`p'=agedays[_n+1]-agedays if subjid_`p'!=""
		by subjid_`p': gen mid_agedays_`p'=0.5*(agedays+agedays[_n+1]) if subjid_`p'!=""
		
		
		**17B**
		gen tanner_months=6+12*(round(mid_agedays_`p'/365.25))
		*Set tanner_months to missing if first age is <30 months
		replace tanner_months=. if agemonths<30
		
		*Set up variables
		gen mindiff_`p'=.
		gen maxdiff_`p'=.
		
		**17C**
		*Merge with Tanner - Tanner values only available for height
		if whichvar[1]=="ht" {
			merge m:1 sex tanner_months using tanner_`p'_vel_rev.dta, nogen
			
			drop if subjid==.

			*Adjust max_`p'_vel because low for older ages
			replace max_`p'_vel=2.54 if max_`p'_vel<2.54
			replace max_`p'_vel=2*2.54 if d_agedays_ht>2*30.4375 & max_`p'_vel<2*2.54
			replace max_`p'_vel=4*2.54 if d_agedays_ht>0.5*365.25 & max_`p'_vel<4*2.54
			replace max_`p'_vel=8*2.54 if d_agedays_ht>365.25 & max_`p'_vel<8*2.54
			*Scale min and max for interval and add extra allowance
			replace mindiff_`p'=0.5*min_`p'_vel*(d_agedays/365.25)^2-3 if d_agedays<365.25
			replace mindiff_`p'=0.5*min_`p'_vel-3 if d_agedays>365.25
			replace maxdiff_`p'=2*max_`p'_vel*(d_agedays/365.25)^1.5+5.5 if d_agedays>365.25
			replace maxdiff_`p'=2*max_`p'_vel*(d_agedays/365.25)^0.33+5.5 if d_agedays<365.25
		}
		***************************************************************************************************************************************ADD

		**17D**

		gen whoagegrp_`p'=round(agedays/30.4375, 1) if agemonths<=24
		*No WHO age group if second age months >24
		sort subjid_`p' agedays
		by subjid_`p': replace whoagegrp_`p'=. if agemonths[_n+1]>24
		*Generate indicator variable for separate intervals
		gen whoinc_age_`p'=.

		
		*WHO values - do analysis belowseparately because 1 month interval not available for HC
		
		*Identify closest WHO interval
		**17E**
		if whichvar[1]=="ht" {
			replace whoinc_age_`p'=1 if d_agedays_`p'>=20 & d_agedays_`p'<46
			replace whoinc_age_`p'=2 if d_agedays_`p'>=46 & d_agedays_`p'<76
			replace whoinc_age_`p'=3 if d_agedays_`p'>=76 & d_agedays_`p'<107
			replace whoinc_age_`p'=4 if d_agedays_`p'>=107 & d_agedays_`p'<153
			replace whoinc_age_`p'=6 if d_agedays_`p'>=153 & d_agedays_`p'<199			
			*If d_agedays under lowest age interval use lowest age interval
			replace whoinc_age_`p'=1 if d_agedays_`p'<20
			*If d_agedays higher than highest interval, use highest interval and replace d_agedays with value one higher than highest age in highest interval
			replace d_agedays_`p'=200 if d_agedays_`p'>199 & d_agedays_`p'!=.												
			replace whoinc_age_`p'=6 if d_agedays_`p'>=200

			
			**17F** 
			*Merge
			merge m:1 sex whoagegrp_`p' using who_`p'_vel_3sd.dta, nogen
			merge m:1 sex whoagegrp_`p' using who_`p'_maxvel_3sd.dta, nogen
			drop if subjid==.
			gen who_mindiff_`p'=.
			gen who_maxdiff_`p'=.
			foreach i in 1 2 3 4 6 {
				replace who_mindiff_`p'=whoinc_`i'_`p' if whoinc_age_`p'==`i'
				replace who_maxdiff_`p'=max_whoinc_`i'_`p' if whoinc_age_`p'==`i'
			}
			
			**17G and H** 
			*Scale WHO based on ratio of d_agedays and age interval
			replace who_mindiff_`p'=who_mindiff_`p'*d_agedays_`p'/(whoinc_age_`p'*30.4375) if d_agedays_`p'<(whoinc_age_`p'*30.4375)
			replace who_maxdiff_`p'=who_maxdiff_`p'*d_agedays_`p'/(whoinc_age_`p'*30.4375) if d_agedays_`p'>(whoinc_age_`p'*30.4375)
			
			*Adjust limits (make mindiff smaller and maxdiff bigger) and add tolerance for measurement error
			*AND choose WHO vs Tanner
			*Use WHO values if WHO is available
			replace mindiff_`p'=0.5*who_mindiff_`p'-3 if who_mindiff_`p'!=. & d_agedays_`p'<(9*30.4375)
			replace maxdiff_`p'=2*who_maxdiff_`p'+3 if who_maxdiff_`p'!=. & d_agedays_`p'<(9*30.4375)
			*If d_agedays >9 months, use Tanner when available, otherwise use WHO
			replace mindiff_`p'=0.5*who_mindiff_`p'-3 if mindiff_`p'==. & who_mindiff_`p'!=.
			replace maxdiff_`p'=2*who_maxdiff_`p'+3 if maxdiff_`p'==. & who_maxdiff_`p'!=.
			*If neither Tanner nor WHO is available, mindiff is measurement error tolerance, and there is no maxdiff
			replace mindiff_`p'=-3 if mindiff_`p'==. & min_`p'_vel==.
			replace maxdiff_`p'=. if who_mindiff_`p'==. & min_`p'_vel==.
			*If birth measurement, make allowance 1.5 cm bigger
			replace mindiff_`p'=mindiff_`p'-1.5 if agedays==0
			replace maxdiff_`p'=maxdiff_`p'+1.5 if agedays==0
			**17I**
			*Generate variable containing min/maxdiff of previous measurement
			foreach m in min max{
				sort subjid_`p' agedays
				by subjid_`p': gen `m'diff_prev_`p'=`m'diff_`p'[_n-1]
			}
		}

		
		**17J**
		if whichvar[1]=="hc" {
			replace whoinc_age_`p'=2 if d_agedays_`p'>=46 & d_agedays_`p'<76
			replace whoinc_age_`p'=3 if d_agedays_`p'>=76 & d_agedays_`p'<107
			replace whoinc_age_`p'=4 if d_agedays_`p'>=107 & d_agedays_`p'<153
			replace whoinc_age_`p'=6 if d_agedays_`p'>=153 & d_agedays_`p'<199 

			**17K**
			merge m:1 sex whoagegrp_`p' using who_`p'_vel_3sd.dta, nogen
			merge m:1 sex whoagegrp_`p' using who_`p'_maxvel_3sd.dta, nogen
			
			drop if subjid==.
			gen who_mindiff_`p'=.
			gen who_maxdiff_`p'=.
			foreach i in 2 3 4 6 {
				replace who_mindiff_`p'=whoinc_`i'_`p' if whoinc_age_`p'==`i'
				replace who_maxdiff_`p'=max_whoinc_`i'_`p' if whoinc_age_`p'==`i'
			}
			
			**17L and M**
			*Scale WHO based on ratio of d_agedays and age interval
			replace who_mindiff_`p'=who_mindiff_`p'*d_agedays_`p'/(whoinc_age_`p'*30.4375) if d_agedays_`p'<(whoinc_age_`p'*30.4375)
			replace who_maxdiff_`p'=who_maxdiff_`p'*d_agedays_`p'/(whoinc_age_`p'*30.4375) if d_agedays_`p'>(whoinc_age_`p'*30.4375)
			
			*Adjust limits and add tolerance for measurement error (don't have to choose between WHO and Tanner)
			replace mindiff_`p'=0.5*who_mindiff_`p'-1.5 if who_mindiff_`p'!=.
			replace maxdiff_`p'=2*who_maxdiff_`p'+1.5 if who_maxdiff_`p'!=. 
			*If WHO is not available, mindiff is measurement error tolerance, and there is no maxdiff
			replace mindiff_`p'=-1.5 if mindiff_`p'==.
			replace maxdiff_`p'=. if who_mindiff_`p'==.
			*If birth measurement, make allowance 0.5 cm bigger
			replace mindiff_`p'=mindiff_`p'-0.5 if agedays==0
			replace maxdiff_`p'=maxdiff_`p'+0.5 if agedays==0
			
			
			**17N**
			*Generate variable containing min/maxdiff of previous measurement
			foreach m in min max{
				sort subjid_`p' agedays
				by subjid_`p': gen `m'diff_prev_`p'=`m'diff_`p'[_n-1]
			}
		}

		**17O** EWMA
		*Determine exponent based on min age difference between value and next value
		sort subjid param agedays
		by subjid param: gen maxagediff_`p'=abs(agedays-agedays[_n-1])
		by subjid param: replace maxagediff_`p'=abs(agedays-agedays[_n+1]) if abs(agedays-agedays[_n+1])>abs(agedays-agedays[_n-1]) & abs(agedays-agedays[_n+1])!=.
		gen expon_`p'=-1.5
		replace expon_`p'=-2.5 if maxagediff_`p'>365.25 & maxagediff_`p'!=.
		replace expon_`p'=-3.5 if maxagediff_`p'>730.5 & maxagediff_`p'!=.
		
		*EWMA
		sort subjid_`p' agedays
		by subjid_`p': gen tot_`p'=_N
		by subjid_`p': gen sn_`p'=_n
		summ tot_`p' if exc_`p'==0
		local m=r(max)
		forvalues n=1/`m' {
			by subjid_`p': gen temp_tbc`p'z_`n'=tbc`p'z[`n'] 
			by subjid_`p': gen temp_age_`p'_`n'=agedays[`n'] 
			gen double temp_ewt_`p'_`n'=abs(`a' + agedays-temp_age_`p'_`n')^expon_`p' if sn_`p'!=`n'
			gen double temp_ewtz_`p'_`n'=temp_tbc`p'z_`n'*temp_ewt_`p'_`n' 
			}
		egen double sumewt_`p'=rowtotal(temp_ewt_`p'_*) if exc_`p'==0
		egen double sumewtz_`p'=rowtotal(temp_ewtz_`p'_*) if exc_`p'==0
		gen ewma_`p'=sumewtz_`p'/sumewt_`p' if subjid_`p'!=""
		gen dewma_`p'=tbc`p'z-ewma_`p'
		gen d_prev_`p'=.
		sort subjid_`p' agedays
		by subjid_`p': replace d_prev_`p'=`p'-`p'[_n-1] if sn_`p'!=1
		gen d_`p'=.
		by subjid_`p': replace d_`p'=d_prev_`p'[_n+1] if sn_`p'!=tot_`p'  
		gen bef_g_aftm1_`p'=.
		gen aft_g_befp1_`p'=.
		*Before & After
		foreach t in bef aft {
			replace ba="`t'"
			gen double sumewt_`p'_`t'=sumewt_`p'
			gen double sumewtz_`p'_`t'=sumewtz_`p'
			gen ewma_`p'_`t'=ewma_`p'
			gen dewma_`p'_`t'=dewma_`p'
			forvalues d=1/`m' {
				  if ba=="bef" local s=`d'+1
				  if ba=="aft" local s=`d'-1
				  replace sumewt_`p'_`t'=sumewt_`p'-temp_ewt_`p'_`d' if temp_ewt_`p'_`d'!=. & sn_`p'==`s'
				  replace sumewtz_`p'_`t'=sumewtz_`p'-temp_ewtz_`p'_`d' if temp_ewtz_`p'_`d'!=. & sn_`p'==`s'
				  replace ewma_`p'_`t'=sumewtz_`p'_`t'/sumewt_`p'_`t' if subjid_`p'!="" & sn_`p'==`s'
				  replace dewma_`p'_`t'=tbc`p'z-ewma_`p'_`t' if subjid_`p'!="" & sn_`p'==`s'
			}
		  }     
		sort subjid_`p' agedays
		
		**17 P Q R S**
		*Identify which has bigger dewma if other value of measurement pair is not taken into account
		gen pair=1 if (d_prev_`p'<mindiff_prev_`p' | d_`p'<mindiff_`p' | d_prev_`p'>maxdiff_prev_`p'  | d_`p'>maxdiff_`p') & exc_`p'==0
		by subjid_`p': replace bef_g_aftm1_`p'=1 if abs(dewma_`p'_bef)>abs(dewma_`p'_aft[_n-1]) & sn_`p'!=1 & pair==1 & pair[_n-1]==1
		by subjid_`p': replace aft_g_befp1_`p'=1 if abs(dewma_`p'_aft)>abs(dewma_`p'_bef[_n+1]) & sn_`p'!=tot_`p' & pair==1 & pair[_n+1]==1

		gen temp_exc_`p'=151 if d_prev_`p'<mindiff_prev_`p' & bef_g_aftm1_`p'==1 & exc_`p'==0 & mindiff_prev_`p'!=.
		gen temp_diff=abs(dewma_`p'_bef) if temp_exc_`p'==151

		replace temp_exc_`p'=152 if d_`p'<mindiff_`p' & aft_g_befp1_`p'==1 & exc_`p'==0 & mindiff_`p'!=.
		replace temp_diff=abs(dewma_`p'_aft) if temp_exc_`p'==152

		replace temp_exc_`p'=161 if d_prev_`p'>maxdiff_prev_`p' & bef_g_aftm1_`p'==1 & exc_`p'==0 & mindiff_prev_`p'!=.
		replace temp_diff=abs(dewma_`p'_bef) if temp_exc_`p'==161

		replace temp_exc_`p'=162 if d_`p'>maxdiff_`p' & aft_g_befp1_`p'==1 & exc_`p'==0 & mindiff_`p'!=.
		replace temp_diff=abs(dewma_`p'_aft) if temp_exc_`p'==162

		
		*Identify more extreme value if only 2 measurements for SP (dewma will be the same)
		gen abstbc`p'z=abs(tbc`p'z)
		replace temp_exc_`p'=153 if d_prev_`p'<mindiff_prev_`p' & abstbc`p'z>abstbc`p'z[_n-1] & tot_`p'==2 & mindiff_prev_`p'!=.     
		replace temp_exc_`p'=154 if d_`p'<mindiff_`p' & abstbc`p'z>abstbc`p'z[_n+1] & tot_`p'==2 & mindiff_`p'!=.
		replace temp_exc_`p'=163 if d_prev_`p'>maxdiff_prev_`p' & abstbc`p'z>abstbc`p'z[_n-1] & tot_`p'==2 & maxdiff_prev_`p'!=.
		replace temp_exc_`p'=164 if d_`p'>maxdiff_`p' & abstbc`p'z>abstbc`p'z[_n+1] & tot_`p'==2 & maxdiff_`p'!=.
		gen temp_bin_exc_`p'=temp_exc_`p'!=.
		
		gsort subjid_`p' -temp_bin_exc_`p' -temp_diff
		
		by subjid_`p': gen dn=_n if temp_bin_exc_`p'==1
		replace exc_`p'=int(temp_exc_`p'/10)+2000 if dn==1
		replace subjid_`p'=subjid_str
		replace subjid_`p'="" if exc_`p'!=0
		
		**17T**
		bysort subjid: egen tot_temp_exc_`p'=total(temp_bin_exc_`p')
		egen max_tot_temp_exc_`p'=max(tot_temp_exc_`p')
		
		**Clean-up**
		local i=max_tot_temp_exc_`p'[1]
		drop ba-max_tot_temp_exc_`p'
	}
	local i 1
}


*Clean-up**
lab def exc 2015 "Min diff" 2016 "Max diff", modify


save temp_abs, replace


foreach f in ewma2 {
	capture confirm file "temp_`f'.dta"
	if !_rc {  
		rm temp_`f'.dta
	}
}

*/

************************************************************************************************
* #19: 1 or 2 measurements
************************************************************************************************


use temp_abs, clear

*Pairs = only 2 nonexcluded values for an SP remain at this step
*Singles = only 1 nonexcluded value for an SP remains at this step

foreach p in wt ht hc {
	sort subjid_`p' agedays
	
	**19A**
	*Determine total # values for parameter and order them
	by subjid_`p': gen vistot_`p'=_N
	by subjid_`p': gen sn_`p'=_n
	
	
	**19B**
	*For pairs, determine: z-score/ageday of other value, diff between values, median of values
	*	and absolute value of z-score
	gen tbc`p'z_other=tbc`p'z[_n+1] if vistot_`p'==2 & sn_`p'==1
	replace tbc`p'z_other=tbc`p'z[_n-1] if vistot_`p'==2 & sn_`p'==2
	gen agedays_`p'_other=agedays[_n+1] if vistot_`p'==2 & sn_`p'==1
	replace agedays_`p'_other=agedays[_n-1] if vistot_`p'==2 & sn_`p'==2
	gen absd_tbc`p'z=abs(tbc`p'z-tbc`p'z_other) if vistot_`p'==2
	gen absd_agedays_`p'=abs(agedays-agedays_`p'_other) if vistot_`p'==2
		gen abs_tbc`p'z=abs(tbc`p'z) if vistot_`p'==2
	*Whether pair/single or not, determine median z-scores (only need for DOP, but will calc for all)
	bysort subjid_`p': egen median_tbc`p'z_i=median(tbc`p'z) if exc_`p'==0
	bysort subjid: egen median_tbc`p'z=min(median_tbc`p'z_i) 
	
	
	**19C**
	*Repeat for corrected z-scores (not ageday values)
	gen ctbc`p'z_other=ctbc`p'z[_n+1] if vistot_`p'==2 & sn_`p'==1
	replace ctbc`p'z_other=ctbc`p'z[_n-1] if vistot_`p'==2 & sn_`p'==2
	gen absd_ctbc`p'z=abs(ctbc`p'z-ctbc`p'z_other) if vistot_`p'==2
	bysort subjid_`p': egen median_ctbc`p'z_i=median(ctbc`p'z) if exc_`p'==0
	bysort subjid: egen median_ctbc`p'z=min(median_ctbc`p'z_i) 
	gen abs_ctbc`p'z=abs(ctbc`p'z) if vistot_`p'==2
}

**19D**
*Generate variables to identify DOP
local oi="1"
foreach p in wt ht hc {
	gen temp`p'=`oi'
	local oi=`oi'+1
}
foreach p in wt ht hc {
	local o empty
	if temp`p'==1 {
		local o ht	
	}
	if temp`p'==2 {
		local o wt
	}
	if temp`p'==3 {
		local o ht
	}

	*When a value for the DOP exists on the same day as the VOI, copmare the DOP and VOI by determining difference in z-score on the same day
	gen abs_d_`p'_`o'=abs(tbc`p'z-tbc`o'z)
	*When there is no value for the DOP on the same day as the VOI, compare the DOP and VOI by determine difference between VOI z-score and the median z for the DOP
	replace abs_d_`p'_`o'=abs(tbc`p'z-median_tbc`o'z) if tbc`o'z==.

	**19E**
	*For pairs, identify the more discrepant value. This is the one with the larger absolute difference between the VOI z-score and DOP (preferentially same-day). If the difference between z-scores are both the same, or they are missing, then identify the value with the larger absolute z-score. 
	sort subjid_`p' abs_d_`p'_`o' abs_tbc`p'z
	by subjid_`p': gen dn_`p'=_n if exc_`p'==0 & vistot_`p'==2
	
	**19F**
	*If the difference between the two z-scores (and corrected if available) for the same SP is >4 and the difference in agedays is >=365.25, exclude the more discrepant value.
	replace exc_`p'=3101 if vistot_`p'==2 & absd_agedays_`p'>=365.25 & absd_tbc`p'z>4 & dn_`p'==2
	
	**19G**
	*If the difference between the two z-scores for the same SP is >2.5 and the difference in agedays is <365.25, exclude the more discrepant value.
	replace exc_`p'=3100 if vistot_`p'==2 & absd_agedays_`p'<365.25 & absd_tbc`p'z>2.5 & dn_`p'==2
	
	**19H**
	*Reidentify how many values there are per SP -- if one measurement was excluded from a pair above, it gets re-evaluated as a single
	drop vistot_`p'
	bysort subjid_`p': gen vistot_`p'=_N if exc_`p'==0
	
	**19I**
	*For singles, exclude if |z|>3 and |diff| between VOI and DOP is >5
	replace exc_`p'=3200 if vistot_`p'==1 & abs(tbc`p'z)>3 & abs(tbc`p'z-tbc`o'z)>5 & tbc`o'z!=. & exc_`p'==0
	replace exc_`p'=3200 if vistot_`p'==1 & abs(tbc`p'z)>3 & tbc`o'z==. & abs(tbc`p'z-median_tbc`o'z)>5 & median_tbc`o'z!=. & exc_`p'==0
	*For singles without any available measurements in the DOP, exclude if the |z| for the VOI is >5
	replace exc_`p'=3200 if vistot_`p'==1 & abs(tbc`p'z)>5 & tbc`o'z==. & median_tbc`o'z==. & exc_`p'==0
	
	replace subjid_`p'=subjid_str
	replace subjid_`p'="" if exc_`p'!=0
	
}

*Clean-up
drop sn_wt-vistot_hc
lab def exc 3101 "2 meas >1 year" 3100 "2 meas <1 year" 3200 "1 meas", add
save temp_1and2, replace

/*																														*******Remove ewma2b later**
foreach f in ewma2b {
	capture confirm file "temp_`f'.dta"
	if !_rc {  
		rm temp_`f'.dta
	}
}
*/

************************************************************************************************
* #21: Error load
************************************************************************************************



use temp_1and2, clear


foreach p in wt ht hc {
	
	*Two binary variables: exc_bin -- excluded for error vs not (either included or excluded for CF or SDE) and inc_bin (included vs not)
	**21A**
	recode exc_`p' 0/1=. 3/7=. 770/774=. nonm=1, gen(exc_bin_`p')
	recode exc_`p' 0=1 nonm=., gen(inc_bin_`p')
	
	**21B**
	bysort subjid: egen tot_exc_`p'=total(exc_bin_`p')
	bysort subjid: egen tot_inc_`p'=total(inc_bin_`p')	
	
	**21C**
	*If ratio of exc_bin/(exc_bin + inc_bin)>0.4 then all values for a parameter get dropped. 	
	replace exc_`p'=4100 if tot_exc_`p'>=2 & tot_exc_`p'/(tot_exc_`p' + tot_inc_`p')>0.4 & exc_`p'==0

	replace subjid_`p'=subjid_str
	replace subjid_`p'="" if exc_`p'!=0
	
}

*Clean-up
drop *_bin_* tot_exc_* tot_inc_*

lab def exc 4100 "Error load", add
save temp_load, replace


*rm temp_abs.dta



************************************************************************************************
* #22: Clean-up (not applicable to R, not an official step)
************************************************************************************************

use temp_load, clear

keep ///
subjid subjid_wt subjid_ht subjid_hc sex agedays param measurement ///
synthetic tester truth testernotes ///
ind_birth ind_hc ind_2y dataset ///
visit_type ga ///
id ///
ageyears agemonths wt ht hc ///
who*z cdc*z s*z tbc*z ctbc*z ///
potcorr* intwt fengadays aga pmagedays cagedays ///
rc*z ///
exc_wt exc_ht exc_hc ///
nnte 



save temp_keep, replace


foreach f in 1and2 {
	capture confirm file "temp_`f'.dta"
	if !_rc {  
		rm temp_`f'.dta
	}
}


************************************************************************************************
* #23: Addend NNTE values
************************************************************************************************

**23A**
use temp_keep, clear

append using temp_nnte_forappend


save temp_clean_unmod_2025-09-30, replace


foreach f in load {
	capture confirm file "temp_`f'.dta"
	if !_rc {  
		rm temp_`f'.dta
	}
}

************************************************************************************************
* #25: Generate dataset for comparison or other use
************************************************************************************************

use temp_clean, clear

save "2025-10-02-gc-Stata-chopdall-results", replace

export delimited using "2025-10-02-gc-Stata-chopdall-results.csv", replace

*Make smaller dataset with fewer variables
/*
keep subjid sex agedays param measurement id

export delimited using "DATE-gc-Stata-results-limitedvars", replace
*/

timer off 1 // Stops timer 1
timer list // Lists the accumulated time for all timers
 



************************************************************************************************
* Source and Modification Info
************************************************************************************************
/*

Based on /Users/carriedaymont/Dropbox/Dropbox Documents/Research/__Chris growthcleanr/gc-revised/Stata/gc-Stata-infant-2025-09-26.do  
 
Running list of modifications to code
N Means this will likely need to be changed in NORC code as well

09-30-2025
- Removed code to create dataset
N - Adjusted rc code to fix bugs
N - Fixed bug in subset data - had to drop ran



*/
