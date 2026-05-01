aafile=( \
  "SCA" "SCV" "SCL" "SCI" "SCC" "SCM" "SCS" "SCT" "SCN" "SCQ" \
  "SCF" "SCY" "SCW" "SCP" "GLYD" \
  "SCHE" "SCHD" "SCDN" "SCEN" "SCKN" "SCRN" \
  "SCHP" "SCD" "SCE" "SCCM" "SCYM" "SCK" "SCR" \
  "SCRN-1" "SCW-1" \
  "NONE" \
)

for aa in "${aafile[@]}"
do
  mkdir POPC-$aa
  cd POPC-$aa
  DIR=../../results/POPC-aa/POPC-$aa
  traj=(1 2 3)
  for t in "${traj[@]}"
  do
    mkdir trajectory$t
    cd trajectory$t
    files=""
    for (( i=401; i<=500; i++ ))
    do
      files="$files $DIR/charmm-gui/namd/out$t/section$i.dcd"
    done
    catdcd -o section401-500.dcd -stride 1 $files
    files=""
    for (( i=501; i<=600; i++ ))
    do
      files="$files $DIR/charmm-gui/namd/out$t/section$i.dcd"
    done
    catdcd -o section501-600.dcd -stride 1 $files
    files=""
    for (( i=601; i<=700; i++ ))
    do
      files="$files $DIR/charmm-gui/namd/out$t/section$i.dcd"
    done
    catdcd -o section601-700.dcd -stride 1 $files
    files=""
    for (( i=701; i<=800; i++ ))
    do
      files="$files $DIR/charmm-gui/namd/out$t/section$i.dcd"
    done
    catdcd -o section701-800.dcd -stride 1 $files
    files=""
    for (( i=801; i<=900; i++ ))
    do
      files="$files $DIR/charmm-gui/namd/out$t/section$i.dcd"
    done
    catdcd -o section801-900.dcd -stride 1 $files
    files=""
    for (( i=901; i<=1000; i++ ))
    do
      files="$files $DIR/charmm-gui/namd/out$t/section$i.dcd"
    done
    catdcd -o section901-1000.dcd -stride 1 $files
    files=""
    cd ..
  done
  cd ..
  mv ../results/POPC-aa/POPC-$aa/charmm-gui/namd/step5_input.psf ./POPC-$aa/popc-${aa,,}.psf
done

