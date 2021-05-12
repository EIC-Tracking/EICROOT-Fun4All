# EICROOT-Fun4All

Hi!

I will be updating the scripts that I have been working on, editing the default ones in https://github.com/eic/EicToyModel/tree/master/fun4all_with_eicroot

// the new edits to get various detector subsystems incorporated are as follows;

copied the EicToyModel directory to my $HOME

singularity shell -B /home  -B /cvmfs /cvmfs/eic.opensciencegrid.org/singularity/rhic_sl7_ext.simg

source /cvmfs/eic.opensciencegrid.org/x8664_sl7/opt/fun4all/core/bin/eic_setup.sh -n

source /cvmfs/eic.opensciencegrid.org/x8664_sl7/opt/fun4all/core/bin/setup_local.sh $HOME/myinstall

export LD_LIBRARY_PATH=$HOME/EicToyModel/build/lib:${OPT_SPHENIX}/vgm/lib64:${LD_LIBRARY_PATH}

source /cvmfs/eic.opensciencegrid.org/x8664_sl7/opt/fun4all/core/bin/setup_local.sh $HOME/EicToyModel/fun4all_with_eicroot

export ROOT_INCLUDE_PATH=$HOME/EicToyModel/build/include/etm:${ROOT_INCLUDE_PATH}






