

C=0.3         # Carbon coverage [Monolayers]

nx=32         # Ni(100) Surface Dimensions (nx x ny)
ny=32     

Ea0=0.25      # Activ. E. for C diffusion [eV]
Eint=-0.10       # Interaction between C_ad neighbors. 
              # 3|Eint| < Ea

Temp=200.     # [K]

hpic=100      # Take picture every hpic frames 
nrun=200000   # KMC steps

idum=7318792
#------------------------------------------------------------------

echo "$nx  $ny  $C" 	         > Tar
echo "$Ea0  $Eint"       	>> Tar
echo "$nrun  $hpic"     	>> Tar
echo "$Temp  $idum"     	>> Tar
    
gfortran -o kmc kmc.f90 
./kmc
