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
DIR=RAWDATA
traj="charmm-gui/namd/out$t"
SCRIPTDIR="/media/bories/Backup/bories/Documents/Travail/results/homoPOPC-aa/analyses/scripts2"
DATA="$DIR/analyses/traj$t/data"
PLOT="$DIR/analyses/traj$t/plot"
PSFFILE="$DIR/charmm-gui/namd/step5_input.psf"
PDBFILE="$DIR/charmm-gui/namd/step5_input.pdb"

equils=(600)

plots="1"	# make plots
help="0"	# print help menu

first=1	# first section
equil=600	# equilibration section
last=1000	# last section
stride10=10	# with stride=10, we have 10K frame over the 1000ns
stride1=1	# with stride=1, we have 100K frame over the 1000ns

# actions (whole trajectory. timeseries)
logfile="0"
center="1"	# get a centered bilayer trajectory
watdel="1"	# delete water from the dcd file
memdel="1"     # delete membrane from the dcd file (mainly to free memory space)
cells="0"	# extract cell dimensions along trajectory
thickness="0"	# compute membrane thickness along trajectory

# bloc analysis
densProf="0"	# density profiles
aadensProf="1"	# amino acid density profiles (="0" for the membrane-only system)
orders="0"	# lipid chain order parameters
solutes="0"     # solute COM probability distribution, and PMFs
solute1=$aa
solute2=""
suppr="0"	# put 1 if you want to suppress big dcd files to free space on disk
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
mkdir -p $DATA $PLOT
cp $PSFFILE ./bilayer.psf

#
# creation of the logfile to calculate the energy
#
if [ "$logfile" != "0" ]; then

  gunzip $DIR/$traj/section*out.gz
  
  files=""
  for (( i=$first; i<=$last; i++ ))
  do
    files="$files $DIR/$traj/section$i.out"
  done
  
  tail -n +269 $files > all.out
  gzip $DIR/$traj/section*out

  grep "ENERGY" all.out | grep -vi I  > energy.dat

  mv energy.dat $DATA/

fi


#
# bilayer centering
#
if [ "$center" != "0" ]; then

  echo "Doing the catdcd file and centering with these options:"
  echo "  First section: $first";
  echo "  Equilibration section: $equil";
  echo "  Last section: $last";
  echo "  Stride: $stride1";

  
  for i in "${equils[@]}"
  do
    files=""
    for (( j=($last-$i); j<=$last; j++ ))
    do
      files="$files $DIR/$traj/section$j.dcd"
    done
    if [ "$watdel" != "0" ]; then
      cat $SCRIPTDIR/indexNoWat.vmd | sed s=PSFFILE=$PSFFILE= | sed s=PDBFILE=$PDBFILE= | sed s=AAA=$aa= > temp.vmd
      
      vmd -dispdev text -e temp.vmd  >& /dev/null
    
      rm temp.vmd
      catdcd -o trajectory.dcd -stride $stride1 -i findexfile.ind $files
      #catdcd -o analog-$aa.dcd -stride 1 -i findexfile.ind $SCRIPTDIR/centered.dcd >& /dev/null
      #rm centered.dcd
      #mv analog-$aa.dcd centered.dcd
      #rm findexfile.ind
      PSFFILE2="./analog-$aa.psf"
      PDBFILE2="./analog-$aa.pdb"
      cat $SCRIPTDIR/center.vmd | sed s=PSFFILE=$PSFFILE2= | sed s=DCDFILE="./trajectory.dcd"= > temp.vmd
      vmd -dispdev text -e $SCRIPTDIR/temp.vmd
      rm trajectory.dcd
      if [ "$memdel" != "0" ]; then
        cat $SCRIPTDIR/indexNoMemb.vmd | sed s=PSFFILE=$PSFFILE2= | sed s=PDBFILE=$PDBFILE2= | sed s=AAA=$aa= > temp.vmd
        vmd -dispdev text -e temp.vmd  >& /dev/null
        catdcd -o trajectory.dcd -stride $stride1 -i findexfile.ind centered.dcd
        rm centered.dcd
        mv trajectory.dcd centered.dcd
      fi
    else
      catdcd -o trajectory.dcd -stride $stride10 $files
      PSFFILE2="./bilayer.psf"
      cat $SCRIPTDIR/center.vmd | sed s=PSFFILE=$PSFFILE2= | sed s=DCDFILE="./trajectory.dcd"= > temp.vmd
    vmd -dispdev text -e $SCRIPTDIR/temp.vmd
    fi
    rm trajectory.dcd temp.vmd
    mv centered.dcd centered$i.dcd
  done
  rm temp.vmd
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

    # extract cell dimensions
    zcat $DIR/$traj/section$i.xsc.gz | tail -n 1 > temp.dat
    x=`cat temp.dat | awk '{print $2}'`
    y=`cat temp.dat | awk '{print $6}'`
    z=`cat temp.dat | awk '{print $10}'`
    echo "$i   $x   $y   $z" >> cell.dat

  done

  if [ "$plots" = "1" ]; then
    gnuplot $SCRIPTDIR/cell.gnu
    mv cell.png $PLOT
  fi

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

  echo "#section    thickness     bilayer-center" > thickness.dat
  tail -n +2 bilayer-thickness.dat > temp.dat
  cat temp.dat | awk '{print $1+1, $2, $3}' >> thickness.dat

  if [ "$plots" = "1" ]; then
    gnuplot $SCRIPTDIR/thickness.gnu
  fi

  mv thickness.dat $DATA
  mv thickness.png $PLOT
  rm bilayer-thickness.dat temp.dat temp.vmd
  echo "    done!"

fi

#
# density profiles
#

if [ "$densProf" != "0" ]; then

  mkdir -p $DATA
  mkdir -p $DATA/densityProfiles
  
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

#
# aa density profiles (bloc analysis)
#

if [ "$aadensProf" != "0" ]; then
  
  mkdir -p $DATA
  mkdir -p $DATA/densityProfiles

  cat $SCRIPTDIR/densityProfiles-aa-total.vmd | sed s=AAA=$solute1=  | sed s=PSFFILE=$PSFFILE2= | sed s=DCDFILE="./trajectory.dcd"= > solutes.vmd
  
  for i in "${equils[@]}"
  do
    ln -sf centered$i.dcd trajectory.dcd
    vmd -dispdev text -e $SCRIPTDIR/solutes.vmd
    rm trajectory.dcd
    mv -f profile-aa.dat  $DATA/densityProfiles/profile-aa-$i.dat
  done

  #rm solutes.vmd

fi


#
# lipid acyl chain order parameters
#

if [ "$orders" != "0" ]; then

  echo "Computing lipid order parameters..."
  mkdir -p $DATA/orderParameters
  cat $SCRIPTDIR/orders-popc.vmd | sed s=PSFFILE=$PSFFILE2= > orders-temp.vmd

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
