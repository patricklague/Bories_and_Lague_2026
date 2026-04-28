#!/bin/bash

# may2022 (already coded)
# center the membrane
# cell dimensions
# membrane thickness (P-P distance, hydrophobic thickness)
# electron densities (water penetration and P-P distance, hydrophobic thickness)

# todo:
# area compressibility modulus
# area per lipid (voronoi tesselation)
# solute distributions (center of mass)
# free energies
# averages and StdDev for bloc analysis (order parameters, electron densities)


# see Pogozheva2022 pour les analyses à réaliser
# electron densities
# order parameters (the 16:0 sn-1 chain of the most abundant lipid in each system)
# distance between lipid head groups, hydrophobic thickness, area per lipid
# depth penetration H2O, sterol tilt angle, area compressibility modulus


# 5) calculate the PMF using the script freeEprofile2.py
#	 python freeEprofile2.py densityProfileSCF-200-300ns.dat 32 35 > freeE-200-300ns.dat
# 5b) calculate normalized distribution functions (brut results before pmf)

# 6) compare with MacCallum values:
#       gnuplot> set xrange [0:]
#       gnuplot> set yrange [-18:8]
#       gnuplot> plot "freeE-200-300ns.dat" u 1:3 w l, "phe-maccallum.dat" u ($1*10.0):2:3 with errorbars


######################
# default values
######################

# directory with the analysis scripts
aa=AMACID
t=NUMTRAJ
nb_aa=NUMAA
DIR=RAWDATA
traj="charmm-gui/namd/out$t"
SCRIPTDIR="."
DATA="$DIR/analyses/traj$t/data"
PSFFILE="$DIR/charmm-gui/namd/step5_input.psf"
PDBFILE="$DIR/charmm-gui/namd/step5_input.pdb"

equils=(600)

plots="1"	# make plots
help="0"	# print help menu

first=1	# first section
last=1000	# last section
stride1=1	# with stride=1, we have 100K frame over the 1000ns

# actions (whole trajectory. timeseries)
center="1"	# get a centered bilayer trajectory
blocs="1"	# get trajectory separated in blocs (401-600, 601-800, 801-1000) / if ="0", the trajectory is taken from 401 to 1000
watdel="1"	# delete water for the dcd file (helpful to save space on disk, and to speed up the analysis) / can't be used to get the density profiles "water" and "total"
cells="1"	# extract cell dimensions along trajectory
thickness="1"	# compute membrane thickness along trajectory

# bloc analysis
densProf="1"	# density profiles
orders="1"	# lipid chain order parameters
solute1=$aa
solute2=""
suppr="1"	# put 1 if you want to suppress big dcd files to free space on disk
    		# or 0 if you want to keep them for further analysis

#########################################################
# do not edit below
#########################################################


############################
# values from command line
############################

for action in "$@"
do
    case "${action}" in
        catdcd) center="1"; echo "Catdcd and center the trajectory for analysis";;
        cell) cells="1" ; echo "Will extract cell dimensions along trajectory";;
        thickness) thickness="1" ; echo "Will compute membrane thickness along trajectory" ;;
        order) orders="1" ; echo "Will compute lipid chain order parameters" ;;
        density) densProf="1" ; echo "Will compute density profiles" ;;
        solutes) solutes="1" ; echo "Will compute the solute probability distributions and PMFs" ;;
        plots) plots="1" ; echo "Will do the plots";;
    esac
done

while getopts f:l:s: flag
do
    case "${flag}" in
        f) first=${OPTARG};;
        l) last=${OPTARG};;
        s) stride=${OPTARG};;
        h) help="1";;
        s1) solute1=${OPTARG};;
        s2) solute2=${OPTARG};;
    esac
done

if [ "$help" != "0" ]; then
  echo "Command line options are (keywords in ()):"
  echo "  (catdcd) the trajectory"
  echo "  Extract (cell) dimensions along trajectory"
  echo "  Compute membrane (thickness) along trajectory"
  echo "  Compute (density) profiles"
  echo "  Compute lipid chain (order) parameters"
  echo "  (solute) density profiles"
  echo "  Do the (plots)"
  echo ""
  echo "Additional arguments for catdcd are (default):"
  echo "  -f (1) : 1st section"
  echo "  -l (2000) : last section"
  echo "  -s (10) : stride value"
  echo ""
  echo "Additional arguments for solutes are (default):"
  echo "  -s1 (SCR)"
  echo "  -ss (SCRN)"
  echo ""
fi

############################
# code for actions
############################

# create directories if new analysis
mkdir -p $DATA
cp $PSFFILE ./bilayer.psf

#
# bilayer centering
#
if [ "$center" != "0" ]; then

  echo "Doing the catdcd file and centering with these options:"
  echo "  First section: $first";
  echo "  Equilibration section: $equil";
  echo "  Last section: $last";
  echo "  Stride: $stride1";

  if [ "$blocs" != "0" ]; then
    files=""
    for (( i=401; i<=600; i++ ))
    do
      files="$files $DIR/$traj/section$i.dcd"
    done

    if [ "$watdel" != "0" ]; then
      cat $SCRIPTDIR/indexNoWat.vmd | sed s=PSFFILE=$PSFFILE= | sed s=PDBFILE=$PDBFILE= | sed s=AAA=$aa= > temp.vmd
      vmd -dispdev text -e temp.vmd  >& /dev/null
      rm temp.vmd
      catdcd -o trajectory.dcd -stride $stride1 -i findexfile.ind $files
      PSFFILE2="./analog-$aa.psf"
      cat $SCRIPTDIR/center.vmd | sed s=PSFFILE=$PSFFILE2= | sed s=DCDFILE="./trajectory.dcd"= > temp.vmd
      vmd -dispdev text -e $SCRIPTDIR/temp.vmd
      mv centered.dcd bloc3.dcd
    else
      catdcd -o trajectory.dcd -stride $stride1 $files
      cat $SCRIPTDIR/center.vmd | sed s=PSFFILE=$PSFFILE= | sed s=DCDFILE="./trajectory.dcd"= > temp.vmd
      vmd -dispdev text -e $SCRIPTDIR/temp.vmd
      mv centered.dcd bloc3.dcd
    fi
      
  
    files=""
    for (( i=601; i<=800; i++ ))
    do
      files="$files $DIR/$traj/section$i.dcd"
    done

    if [ "$watdel" != "0" ]; then
      cat $SCRIPTDIR/indexNoWat.vmd | sed s=PSFFILE=$PSFFILE= | sed s=PDBFILE=$PDBFILE= | sed s=AAA=$aa= > temp.vmd
      vmd -dispdev text -e temp.vmd  >& /dev/null
      rm temp.vmd
      catdcd -o trajectory.dcd -stride $stride1 -i findexfile.ind $files
      PSFFILE2="./analog-$aa.psf"
      cat $SCRIPTDIR/center.vmd | sed s=PSFFILE=$PSFFILE2= | sed s=DCDFILE="./trajectory.dcd"= > temp.vmd
      vmd -dispdev text -e $SCRIPTDIR/temp.vmd
      mv centered.dcd bloc4.dcd
    else
      catdcd -o trajectory.dcd -stride $stride1 $files
      cat $SCRIPTDIR/center.vmd | sed s=PSFFILE=$PSFFILE= | sed s=DCDFILE="./trajectory.dcd"= > temp.vmd
      vmd -dispdev text -e $SCRIPTDIR/temp.vmd
      mv centered.dcd bloc4.dcd
    fi
 
    files=""
    for (( i=801; i<=1000; i++ ))
    do
      files="$files $DIR/$traj/section$i.dcd"
    done

    if [ "$watdel" != "0" ]; then
      cat $SCRIPTDIR/indexNoWat.vmd | sed s=PSFFILE=$PSFFILE= | sed s=PDBFILE=$PDBFILE= | sed s=AAA=$aa= > temp.vmd
      vmd -dispdev text -e temp.vmd  >& /dev/null
      rm temp.vmd
      catdcd -o trajectory.dcd -stride $stride1 -i findexfile.ind $files
      PSFFILE2="./analog-$aa.psf"
      cat $SCRIPTDIR/center.vmd | sed s=PSFFILE=$PSFFILE2= | sed s=DCDFILE="./trajectory.dcd"= > temp.vmd
      vmd -dispdev text -e $SCRIPTDIR/temp.vmd
      mv centered.dcd bloc5.dcd
    else
      catdcd -o trajectory.dcd -stride $stride1 $files
      cat $SCRIPTDIR/center.vmd | sed s=PSFFILE=$PSFFILE= | sed s=DCDFILE="./trajectory.dcd"= > temp.vmd
      vmd -dispdev text -e $SCRIPTDIR/temp.vmd
      mv centered.dcd bloc5.dcd
    fi
    
  else
    for i in "${equils[@]}"
    do
      files=""
      for (( j=($last-$i+1); j<=$last; j++ ))
      do
        files="$files $DIR/$traj/section$j.dcd"
      done
      if [ "$watdel" != "0" ]; then
        cat $SCRIPTDIR/indexNoWat.vmd | sed s=PSFFILE=$PSFFILE= | sed s=PDBFILE=$PDBFILE= | sed s=AAA=$aa= > temp.vmd
        vmd -dispdev text -e temp.vmd  >& /dev/null
    
        rm temp.vmd
        catdcd -o trajectory.dcd -stride $stride1 -i findexfile.ind $files
        PSFFILE2="./analog-$aa.psf"
        cat $SCRIPTDIR/center.vmd | sed s=PSFFILE=$PSFFILE2= | sed s=DCDFILE="./trajectory.dcd"= > temp.vmd
        vmd -dispdev text -e $SCRIPTDIR/temp.vmd
      else
        catdcd -o trajectory.dcd -stride $stride10 $files
        PSFFILE2="./bilayer.psf"
        cat $SCRIPTDIR/center.vmd | sed s=PSFFILE=$PSFFILE2= | sed s=DCDFILE="./trajectory.dcd"= > temp.vmd
        vmd -dispdev text -e $SCRIPTDIR/temp.vmd
      fi
      rm trajectory.dcd temp.vmd
      mv centered.dcd centered$i.dcd
    done
  fi
fi

if [ "$watdel" != "0" ]; then
  PSFFILE2="./analog-$aa.psf"
else
  PSFFILE2="./bilayer.psf"
fi


#
# cell dimensions (whole traj)
#

if [ "$cells" = "1" ]; then

  echo "Extracting cell dimensions..."
  echo "#section    x     y    z" > cell.dat

  for (( i=$first; i<=$last; i++ ))
  do

    # extract cell dimensions (files not distributed but generated during the simulation, see readme.md for further details)
    zcat $DIR/$traj/section$i.xsc.gz | tail -n 1 > temp.dat
    x=`cat temp.dat | awk '{print $2}'`
    y=`cat temp.dat | awk '{print $6}'`
    z=`cat temp.dat | awk '{print $10}'`
    echo "$i   $x   $y   $z" >> cell.dat

  done

  mv cell.dat $DATA
  rm temp.dat
  echo "    done!"
    
fi


#
# membrane thickness (whole traj)
#

if [ "$thickness" = "1" ]; then


  echo "Computing bilayer thickness..."
  cat $SCRIPTDIR/thickness.vmd | sed s=PSFFILE=$PSFFILE2= | sed s=DCDFILE="./trajectory.dcd"= > temp.vmd
  ln -sf centered$i.dcd trajectory.dcd
  vmd -dispdev none -e temp.vmd

  echo "#frame    thickness     bilayer-center" > thickness.dat
  tail -n +2 bilayer-thickness.dat > temp.dat
  cat temp.dat | awk '{print $1, $2, $3}' >> thickness.dat

  if [ "$plots" = "1" ]; then
    gnuplot $SCRIPTDIR/thickness.gnu
  fi

  mv thickness.dat $DATA
  mv thickness.png $PLOT
  rm bilayer-thickness.dat temp.dat temp.vmd trajectory.dcd
  echo "    done!"

fi

#
# density profiles (bloc analysis)
#

if [ "$densProf" != "0" ]; then

  mkdir -p $DATA
  mkdir -p $DATA/densityProfiles

  if [ "$blocs" != "0" ]; then
    for i in 3 4 5
    do
      ln -sf bloc$i.dcd trajectory.dcd
      if [ "$watdel" != "0" ]; then
        cat $SCRIPTDIR/densityProfiles-popc.vmd | sed s=PSFFILE=$PSFFILE2= > temp.vmd
      else
        cat $SCRIPTDIR/densityProfiles-water-total.vmd | sed s=PSFFILE=$PSFFILE2= > temp.vmd
      fi
      vmd -dispdev text -e $SCRIPTDIR/temp.vmd
      rm trajectory.dcd
      if [ "$watdel" != "0" ]; then
        mv profile-carbonyl.dat  $DATA/densityProfiles/profile-carbonyl$i.dat
        mv profile-chains.dat  $DATA/densityProfiles/profile-chains$i.dat
        mv profile-choline.dat  $DATA/densityProfiles/profile-choline$i.dat
        mv profile-popc.dat  $DATA/densityProfiles/profile-popc$i.dat
        mv profile-phosphate.dat  $DATA/densityProfiles/profile-phosphate$i.dat
      else
        mv profile-total.dat  $DATA/densityProfiles/profile-total$i-2.dat
        mv profile-water.dat $DATA/densityProfiles/profile-water$i-2.dat
      fi
    done
  else
    for i in "${equils[@]}"
    do
      ln -sf centered$i.dcd trajectory.dcd
      if [ "$watdel" != "0" ]; then
        cat $SCRIPTDIR/densityProfiles-popc.vmd | sed s=PSFFILE=$PSFFILE2= > temp.vmd
      else
        cat $SCRIPTDIR/densityProfiles-water-total.vmd | sed s=PSFFILE=$PSFFILE2= > temp.vmd
      fi
      vmd -dispdev text -e $SCRIPTDIR/temp.vmd
      rm trajectory.dcd
      if [ "$watdel" != "0" ]; then
        mv profile-carbonyl.dat  $DATA/densityProfiles/profile-carbonyl$i.dat
        mv profile-chains.dat  $DATA/densityProfiles/profile-chains$i.dat
        mv profile-choline.dat  $DATA/densityProfiles/profile-choline$i.dat
        mv profile-popc.dat  $DATA/densityProfiles/profile-popc$i.dat
        mv profile-phosphate.dat  $DATA/densityProfiles/profile-phosphate$i.dat
      else
        mv profile-total.dat  $DATA/densityProfiles/profile-total$i.dat
        mv profile-water.dat $DATA/densityProfiles/profile-water$i.dat
      fi
    done
    rm temp.vmd
  fi
fi


#
# lipid acyl chain order parameters
#

if [ "$orders" != "0" ]; then

  echo "Computing lipid order parameters..."
  mkdir -p $DATA/orderParameters
  cat $SCRIPTDIR/orders-popc.vmd | sed s=PSFFILE=$PSFFILE2= > orders-temp.vmd
  if [ "$blocs" != "0" ]; then
    for i in 3 4 5
    do

      ln -sf bloc$i.dcd trajectory.dcd
      vmd -dispdev text -e $SCRIPTDIR/orders-temp.vmd
      rm trajectory.dcd
      mv orderparameters-chain2.dat_plot\~ $DATA/orderParameters/orderparameters-chain2$i.dat
      mv orderparameters-chain3.dat_plot\~ $DATA/orderParameters/orderparameters-chain3$i.dat

    done
  else
    for i in "${equils[@]}"
    do
    
      ln -sf centered$i.dcd trajectory.dcd
      vmd -dispdev text -e $SCRIPTDIR/orders-temp.vmd
      rm trajectory.dcd
      mv orderparameters-chain2.dat_plot\~ $DATA/orderParameters/orderparameters-chain2$i.dat
      mv orderparameters-chain3.dat_plot\~ $DATA/orderParameters/orderparameters-chain3$i.dat

    done
    rm orders-temp.vmd
    echo "  Done!"
    echo ""
  fi

fi

    
#    
# suppress all files created   
#
    
if [ "$suppr" != "0" ]; then

  echo "Making space in the disk : suppressing extra files..." $solute1 #$solute2
      
  rm $SCRIPTDIR/*.dcd
  rm -r networks/
  rm findexfile.ind analog*
  

  echo "  Done!"
  echo ""

fi

#### Fin de script ####


exit
