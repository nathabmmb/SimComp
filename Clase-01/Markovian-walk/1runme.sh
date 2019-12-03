
idum=1225168
A="1. 2."
B="11. 16."
Temp=2.


# ---------------------------------------------
gfortran -o Markovian Markovian.f90


echo $idum  > Tar
echo $A    >> Tar
echo $B    >> Tar
echo $Temp >> Tar

./Markovian
