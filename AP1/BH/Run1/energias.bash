grep -e "-4" BH_Structures/* | sort -k2
rm Energies.dat
grep -e "-4" BH_Structures/* | awk '{print $2}' >> Energies.dat