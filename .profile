# ~/.profile: executed by the command interpreter for login shells.
# This file is not read by bash(1), if ~/.bash_profile or ~/.bash_login
# exists.
# see /usr/share/doc/bash/examples/startup-files for examples.
# the files are located in the bash-doc package.

# the default umask is set in /etc/profile; for setting the umask
# for ssh logins, install and configure the libpam-umask package.
#umask 022

# if running bash
if [ -n "$BASH_VERSION" ]; then
    # include .bashrc if it exists
    if [ -f "$HOME/.bashrc" ]; then
	. "$HOME/.bashrc"
    fi
fi

# set PATH so it includes user's private bin if it exists
if [ -d "$HOME/bin" ] ; then
    PATH="$HOME/bin:$PATH"
fi

# set PATH so it includes user's private bin if it exists
if [ -d "$HOME/.local/bin" ] ; then
    PATH="$HOME/.local/bin:$PATH"
fi


## DAVID
#source /usr/local/progs/geant4/geant4.10.07.p02-install/bin/geant4.sh
export PATH="/usr/local/progs/gate/gate-9.1-install/bin:/usr/local/progs/ImageJ:$PATH"
export LD_LIBRARY_PATH="/usr/local/progs/geant4/geant4.10.07.p02-install/lib:/usr/local/progs/gate/gate-9.1-src/libtorch/lib:/usr/local/progs/root/root_v6.24.06/lib:$LD_LIBRARY_PATH"

export G4ENSDFSTATEDATA="/usr/local/progs/geant4/geant4.10.07.p02-install/share/Geant4-10.7.2/data/G4ENSDFSTATE2.3"
export G4PIIDATA="/usr/local/progs/geant4/geant4.10.07.p02-install/share/Geant4-10.7.2/data/G4PII1.3"
export G4INCLDATA="/usr/local/progs/geant4/geant4.10.07.p02-install/share/Geant4-10.7.2/data/G4INCL1.0"
export G4LEDATA="/usr/local/progs/geant4/geant4.10.07.p02-install/share/Geant4-10.7.2/data/G4EMLOW7.13"
export G4PARTICLEXSDATA="/usr/local/progs/geant4/geant4.10.07.p02-install/share/Geant4-10.7.2/data/G4PARTICLEXS3.1.1"
export G4NEUTRONHPDATA="/usr/local/progs/geant4/geant4.10.07.p02-install/share/Geant4-10.7.2/data/G4NDL4.6"
export G4SAIDXSDATA="/usr/local/progs/geant4/geant4.10.07.p02-install/share/Geant4-10.7.2/data/G4SAIDDATA2.0"
export G4REALSURFACEDATA="/usr/local/progs/geant4/geant4.10.07.p02-install/share/Geant4-10.7.2/data/RealSurface2.2"
export G4ABLADATA="/usr/local/progs/geant4/geant4.10.07.p02-install/share/Geant4-10.7.2/data/G4ABLA3.1"
export G4LEVELGAMMADATA="/usr/local/progs/geant4/geant4.10.07.p02-install/share/Geant4-10.7.2/data/PhotonEvaporation5.7"
export G4RADIOACTIVEDATA="/usr/local/progs/geant4/geant4.10.07.p02-install/share/Geant4-10.7.2/data/RadioactiveDecay5.6"

