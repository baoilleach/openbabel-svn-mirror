#!/bin/sh
# test Open Babel with roundtripping if the test suite exists
#   (just un-tar into this folder as "test-set" and it'll run)

echo
if [ -d test-set ]; then
    echo "Testing molecule file roundtripping..."
    echo "This test will take quite some time!"
    echo " and it currently always passes, so check the results file"
    cd test-set
    source test.sh
    cd ..
else
    echo "Roundtrip test set not found. Skipping test."
    exit 77
fi
