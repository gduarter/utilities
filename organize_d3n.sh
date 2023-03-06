
prot=()
anchors=()
dirs=("01.files" "02.dockprep" "03.spheres" "04.grid" "05.ref_footprint" "06.multigrid" "07.d3n" "08.cartesianmin" "09.rescore" "10.dock")

for val in "${prot[@]}" #PDB code
do
    for elem in "${dirs[@]}" 
    do
        mkdir -p ${root}/${val}/${elem}
    done
done
