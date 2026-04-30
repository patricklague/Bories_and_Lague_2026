#! /bin/bash
#SBATCH --account=def-plague
#SBATCH --job-name=BOR
#SBATCH --gpus-per-node=1          # Number of GPU(s) per node
#SBATCH --cpus-per-task=8         # CPU cores/threads
#SBATCH --mem 2048               # memory per node
#SBATCH --time=0-02:00            # time (DD-HH:MM)
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# modules may vary according to the stack, check with "module avail"
module load StdEnv/2020
module load intel/2020.1.217
module load cuda/11.4
module load namd-multicore/2.14
# Number of DCD sections 
# time of trajectory = $nstep * $nsection * timestep
# nsteps is now set in step5 file
NSECTION=1000

# namd command, adjust +p10 to number of available cores as requested to queue system
#NAMD="/home/plague/projects/def-plague/bin/NAMD_Git-2018-06-14_Linux-x86_64-multicore-CUDA/namd2 +p5 +idlepoll"
NAMD="namd2 +p$SLURM_CPUS_PER_TASK +idlepoll"
out=1 #2 and 3 (correspond to each trajectory)

#
# Do not edit below...
#

# need to start the trajectory?
if [ ! -f out$out/section1.dcd ]; then  
    if [ $out -eq 1 ]; then
        #start trajectory
        $NAMD step6.1_equilibration.inp >& step6.1_equilibration.out
        $NAMD step6.2_equilibration.inp >& step6.2_equilibration.out
        $NAMD step6.3_equilibration.inp >& step6.3_equilibration.out
        $NAMD step6.4_equilibration.inp >& step6.4_equilibration.out
        $NAMD step6.5_equilibration.inp >& step6.5_equilibration.out
        $NAMD step6.6_equilibration.inp >& step6.6_equilibration.out
    fi
    
    mkdir -p out$out
    cp step6.6_equilibration.xst out$out/restart.xst
    cp step6.6_equilibration.xsc out$out/restart.xsc
    cp step6.6_equilibration.coor out$out/restart.coor
    cp step6.6_equilibration.vel out$out/restart.vel

fi

#
# restart trajectory
#

# find the number of the last section
SECTION=0
if [ -f out$out/section1.dcd ]; then
    for file in `ls out$out/section*vel.gz`;
        do
        temp=`echo $file | awk -F . '{ print $1 }' | awk -F on '{ print $2 }'`
        if  [ "$temp" -gt "$SECTION" ]; then
            SECTION=$temp
        fi
    done
    echo Section $SECTION will be used to restart the trajectory
fi

# make sure we have the right files for restart
if [ $SECTION -gt 0 ]; then

    rm -f out$out/restart*
    cp -f out$out/section$SECTION.xst.gz out$out/restart.xst.gz
    cp -f out$out/section$SECTION.coor.gz out$out/restart.coor.gz
    cp -f out$out/section$SECTION.vel.gz out$out/restart.vel.gz
    cp -f out$out/section$SECTION.xsc.gz out$out/restart.xsc.gz
    gunzip out$out/restart*
fi
let COUNTER=SECTION+1

# create appropriate input file using step5_production.inp as template
inputname="out$out/restart"
outputname="out$out/section$COUNTER"
sed "s/step6.6_equilibration/out$out\/restart/" step7_production.inp | \
  sed "s/step7_production/out$out\/section$COUNTER/" > step7_run.inp

# run the simulation for 1 nanosecond
$NAMD step7_run.inp > out$out/section$COUNTER.out
gzip -f out$out/section$COUNTER.xst
gzip -f out$out/section$COUNTER.coor
#gzip out$out/section$COUNTER.dcd
gzip -f out$out/section$COUNTER.vel
gzip -f out$out/section$COUNTER.xsc
gzip -f out$out/section$COUNTER.out

#rm out$out/*restart*
if [ $COUNTER -le $NSECTION ]; then
  /opt/software/slurm/bin/sbatch submit-alliancecan.sh > jobid
fi

exit

