# jpsi_pA
This code computes the correlation of gluon (charged hadron) with J/Psi production at different rapidities, studying the effect of e-by-e fluctuations.

usage:<br>
To generate the MV tables and then run the calculation:<br>
./main <br>
or<br>
./main --readTable 0<br>
To read in the MV tables from file and then do the calculation (faster):<br>
./main --readTable 1<br>

To include fluctuations use option<br>
--fluctuations 1
