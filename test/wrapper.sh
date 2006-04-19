#!/bin/sh

# Run "prove" on all Perl programs

TESTS="matrix.pl unitcell.pl rings.pl smarts.pl aromatest.pl cml.pl test-set.pl"
PROVE=prove

if [ -d ../src/formats/.libs ]; then
    if [ "x${BABEL_LIBDIR}" = "x" ]; then
	BABEL_LIBDIR="`pwd`/../src/formats/.libs:`pwd`/../src/formats/xml/.libs"
	export BABEL_LIBDIR
    fi
fi

prove ${TESTS}
