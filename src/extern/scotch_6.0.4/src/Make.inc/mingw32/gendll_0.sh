#! /bin/sh
# (C) 2008 Yves Secretan (yves.secretan@ete.inrs.ca)
# This software is governed by the CeCILL-C license under French law
# and abiding by the rules of distribution of free software. You can
# use, modify and/or redistribute the software under the terms of the
# CeCILL-C license as circulated by CEA, CNRS and INRIA at the following
# URL: "http://www.cecill.info".
#
# To be executed in a MSYS window.
#
# This file will link the PT-Scotch DLL. It must be adapted to reflect
# your environment, in particular library path and library name.

export OBJS="../../libscotch/*.o"
export LDFLAGS="--shared -Wl,--allow-multiple-definition"
export PGMFILES="/l"
export LDPATH="-L../../../lib"
export LDLIBS="-lptscotch -lptscotcherr"

gcc --output ../../../lib/libptscotch.dll $LDFLAGS $LDPATH $OBJS $LDLIBS
