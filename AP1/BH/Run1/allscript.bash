for j in {1..9}
do

    for i in {0..20}
    do
        rm BH_output.txt
        ../BH_Lammps/bh.x lammps
        # ./energias.bash
        # gnuplot --persist ploteador &
        python notfinder.py
    done
        python notaverager.py