#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12     	# number of MPI processes per node
#SBATCH --mem=0                   	# memory; default unit is megabytes, specify 0 for all on node.
#SBATCH --mail-type=ALL   		#email notification for all events
#SBATCH --mail-user=andrew.mcfadden@ucalgary.ca
#SBATCH --time=1-00:00           	# Maximum run time in hh:mm:ss, or d-hh:mm
#SBATCH --partition=parallel		# lattice or single for small jobs (<8cores)/cpu2019 for big jobs
#SBATCH --job-name=rsm_lrr_transInit_PIVJetSim	#Job name
FOAM=`which pimpleCentralFoam`		#Foam Solver

# 2018-05-31 DSP - Example batch job script for a parallel OpenFOAM run on ARC

# Initialize OpenFOAM environment.
module purge
module load gcc/10.2.0;
module load openmpi/4.1.1-gnu
source $HOME/OpenFOAM/OpenFOAM-v2012/etc/bashrc

export FOAM_RUN=$PWD

echo "Working in $PWD"

CORES=$SLURM_NTASKS
echo "Running on $CORES cores."

echo "Make a new decomposeParDict file"
DATE=`date`

cat > system/decomposeParDict <<EOF
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      decomposeParDict;
}

// decomposeParDict created at ${DATE}.

numberOfSubdomains $CORES;

method          scotch;

EOF

echo "Forcing new decomposition"

decomposePar -force

echo "Using mpiexec: `which mpiexec`"

echo "About to run $FOAM at `date`"

mpiexec $FOAM -parallel > ${SLURM_JOB_NAME}_scotch_arc_${CORES}cores_${SLURM_JOB_ID}.out

echo "Finished at `date`"

echo "Running reconstructPar at `date`."
reconstructPar -newTimes
echo "Finished reconstructPar at `date`."
echo "Manually delete processor directories if reconstruction succeeded. "


