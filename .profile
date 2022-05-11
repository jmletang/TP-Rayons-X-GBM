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
if test -e /usr/local/progs/geant4/geant4.10.07.p02-install/bin/geant4.sh
then
	source /usr/local/progs/geant4/geant4.10.07.p02-install/bin/geant4.sh
fi
export PATH="/usr/local/progs/gate/gate-9.1-install/bin:/usr/local/progs/ImageJ:$PATH"
export LD_LIBRARY_PATH="/usr/local/progs/gate/gate-9.1-src/libtorch/lib:/usr/local/progs/root/root_v6.24.06/lib:$LD_LIBRARY_PATH"

