# Per Protocol Analysis Code for the EAGeR Trial

<img src="https://github.com/ainaimi/EAGeR-PerProtocol/blob/main/EAGeR.gif" align="right"
     alt="EAGeR Logo" width="120" height="75">

This repository contains code needed to reproduce analyses used to obtain per
protocol findings from the Effects of Aspirin on Gestation and Reproduction Trial. 

The code was run on a computing cluster via the g_submission shell script. This 
shell script run the R coded that analyses the EAGeR data with g computation algorithm, 
estimated via the monte carlo estimator, and implements the bootstrap to obtain 95% CIs
for all contrasts provided in the manuscript. 

The basic command used to submit the g computation code (gF.R) is as follows:

``` Rscript --no-save --no-restore --verbose ./R/gF.R [ARGUMENTS]  > gF_NaturalCourse.Rout 2>&1 ```

Where the `[ARGUMENTS]` option include values for various objects that enable the R 
program to run. 

The arguments used in the `gF.R` program include the following, presented in the order
they appear in the shell script file:

- `thresh`: number of days per week threshold at which a woman was deemed "adherent". Set at a default of 5/7 = 0.7142857, but impact of alternative values was explored.
- `bootNum`: partial  number of bootstrap resamples used to compute variance of point estimates. This number should be multiplied by `array_num` for the total number of resamples.
- `montecarlo`: size of the monte carlo resample used to implement the g computation algorithm.
- `array_num`: number of arrays over which the submission jobs were spread. This is used to avoid delays in obtaining program outputs. 
- `strat`: used to determine whether all regression models will be stratified by randomization stratum (`strat = 1`) or eligibility stratum (`strat = 2`) 
- `rand`: used in the g computation algorithm to set all individuals to a given randomization stratum.
- `expo`: used in the g computation algorithm to set all individuals to a given adherence level.
- `cens`: used in the g computation algorithm to adjust for withdrawal.
- `int`: used in the g computation algorithm to evaluate impact of initiating adherence at varying weeks post conception.


This work was funded by an Extramural Research Grant from the Eunice Kennedy Shriver National Institute of Child Health and Human Development, grant number R01 HD093602, and by the Intramural Research Program of the Eunice Kennedy Shriver National Institute of Child Health and Human Development, National Institues of Health, Bethesda, Maryland (contract numbers HHSN267200603423, HHSN267200603424, and HHSN267200603426).
