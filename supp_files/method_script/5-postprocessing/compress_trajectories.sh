#A FAIRE: "SCY" "SCI" "SCL" "SCF" "SCS" "SCT" "SCKN" "SCRN" "SCDN" "SCEN" "SCQ" "SCN" "SCHE" "SCHD" "SCC" "SCM" "SCYM"
#(-N): "SCHP" "SCK" "SCR" "SCCM" "SCD" "SCE" #dens_prof (en cours) à faire (SCE et SCD à faire)
#traj4:"SCV" "GLYD" "SCA" "SCP" "SCW" #dens_prof à faire
#aafile=("SCHP" "SCK" "SCR" "SCCM" "SCD" "SCE" "SCV" "GLYD" "SCA" "SCP" "SCW")
aafile=("SCY" "SCI" "SCL" "SCF" "SCS" "SCT" "SCKN" "SCRN" "SCDN" "SCEN" "SCQ" "SCN" "SCHE" "SCHD" "SCC" "SCM" "SCYM" "SCHP" "SCK" "SCR" "SCCM" "SCD" "SCE" "SCV" "GLYD" "SCA" "SCP" "SCW" "SCRN-1" "SCW-1" "NONE")
#aafile=("NONE")
#unfinished : SCHD SCI SCKN SCK-N SCL SCM SCN SCQ SCRN SCS SCT SCY SCYM
#aafile=("SCDN" "SCHE" "SCHP" "SCR" "SCD" "SCE" "SCV" "GLYD" "SCA" "SCP" "SCW" "SCCM" "SCEN" "SCF" "NONE" "SCC") #finished

for aa in "${aafile[@]}"
do
  #mkdir POPC-$aa
  cd POPC-$aa
  DIR=/media/bories/Backup/bories/Documents/Travail/results/homoPOPC-aa/homoPOPC-$aa
  traj=(1 2 3)
  if [[ "$aa" == "SCHP" || "$aa" == "SCK" || "$aa" == "SCR" || "$aa" == "SCCM" || "$aa" == "SCD" || "$aa" == "SCE" ]]; then
    DIR=/media/bories/Backup/bories/Documents/Travail/results/homoPOPC-aa/homoPOPC-$aa-N
  elif [[ "$aa" == "GLYD" || "$aa" == "SCV" || "$aa" == "SCA" || "$aa" == "SCP" || "$aa" == "SCW" ]]; then
    traj=(4 2 3)
  elif [[ "$aa" == "SCYM" ]]; then
    traj=(4 5 6)
  fi
  for t in "${traj[@]}"
  do
    #mkdir trajectory$t
    if [[ $t == 4 ]]; then
      #mv trajectory$t trajectory1
      cd trajectory1
    elif [[ $t == 5 ]]; then
      #mv trajectory$t trajectory2
      cd trajectory2
    elif [[ $t == 6 ]]; then
      #mv trajectory$t trajectory3
      cd trajectory3
    else
      cd trajectory$t
    fi
    files=""
    for (( i=401; i<=500; i++ ))
    do
      files="$files $DIR/charmm-gui/namd/out$t/section$i.dcd"
    done
    #catdcd -o section401-500.dcd -stride 1 $files
    files=""
    for (( i=501; i<=600; i++ ))
    do
      files="$files $DIR/charmm-gui/namd/out$t/section$i.dcd"
    done
    #catdcd -o section501-600.dcd -stride 1 $files
    files=""
    for (( i=601; i<=700; i++ ))
    do
      files="$files $DIR/charmm-gui/namd/out$t/section$i.dcd"
    done
    #catdcd -o section601-700.dcd -stride 1 $files
    files=""
    for (( i=701; i<=800; i++ ))
    do
      files="$files $DIR/charmm-gui/namd/out$t/section$i.dcd"
    done
    #catdcd -o section701-800.dcd -stride 1 $files
    files=""
    for (( i=801; i<=900; i++ ))
    do
      files="$files $DIR/charmm-gui/namd/out$t/section$i.dcd"
    done
    #catdcd -o section801-900.dcd -stride 1 $files
    files=""
    for (( i=901; i<=1000; i++ ))
    do
      files="$files $DIR/charmm-gui/namd/out$t/section$i.dcd"
    done
    #catdcd -o section901-1000.dcd -stride 1 $files
    files=""
    #gzip section*dcd
    cd ..
  done
  cd ..
  #cp $DIR/charmm-gui/namd/step5_input.psf ./POPC-$aa/popc-${aa,,}.psf
  mv ./POPC-$aa/step5_input.psf ./POPC-$aa/popc-${aa,,}.psf
done

