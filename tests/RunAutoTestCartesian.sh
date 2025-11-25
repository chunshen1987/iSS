#!/usr/bin/env bash

cp tests/music_input_91 tests/music_input

cp tests/iSS_parameters_ideal.dat tests/iSS_parameters_ideal_fixSeed.dat
# Change the seed in the parameter file (iSS_parameters_ideal_fixSeed.dat) to 42 for this test
sed -i '/^[[:space:]]*randomSeed[[:space:]]*=/s/-1/42/' tests/iSS_parameters_ideal_fixSeed.dat

# Verify the seed was changed; fail early if not.
if ! grep -qE '^[[:space:]]*randomSeed[[:space:]]*=[[:space:]]*42' tests/iSS_parameters_ideal_fixSeed.dat; then
	echo "ERROR: failed to set random seed to 42 in tests/iSS_parameters_ideal_fixSeed.dat" >&2
	echo "Current matching lines:" >&2
	grep -nE 'randomSeed' tests/iSS_parameters_ideal_fixSeed.dat >&2 || true
	exit 1
else
	echo "random seed set to 42 in tests/iSS_parameters_ideal_fixSeed.dat"
fi

# Run the test with fluid cell 5 (time is 1 in Milne coordinates)
./iSS.e tests/iSS_parameters_ideal_fixSeed.dat tests testIdealOneFluidCell5.dat
# Move the test output to a different name with _test1 suffix
mv check_211_spectra.dat check_211_spectra_test1.dat
mv check_2212_spectra.dat check_2212_spectra_test1.dat
mv checkReconstructedTmunu.dat checkReconstructedTmunu_test1.dat

# Run the test with fluid cell 5 in Cartesian coordinates (time is 1)
./iSS.e tests/iSS_parameters_ideal_Cartesian.dat tests testIdealOneFluidCell5.dat
# Move the test output to a different name with _test2 suffix
mv check_211_spectra.dat check_211_spectra_test2.dat
mv check_2212_spectra.dat check_2212_spectra_test2.dat
mv checkReconstructedTmunu.dat checkReconstructedTmunu_test2.dat

# Run the test with fluid cell 6 (time is 1 and z=5 in Cartesian coordinates)
./iSS.e tests/iSS_parameters_ideal_Cartesian.dat tests testIdealOneFluidCell6.dat
# Move the test output to a different name with _test3 suffix
mv check_211_spectra.dat check_211_spectra_test3.dat
mv check_2212_spectra.dat check_2212_spectra_test3.dat
mv checkReconstructedTmunu.dat checkReconstructedTmunu_test3.dat

# Run the test with fluid cell 1 (time is 10 in Milne coordinates)
./iSS.e tests/iSS_parameters_ideal_fixSeed.dat tests testIdealOneFluidCell1.dat
# Move the test output to a different name with _test4 suffix
mv check_211_spectra.dat check_211_spectra_test4.dat
mv check_2212_spectra.dat check_2212_spectra_test4.dat
mv checkReconstructedTmunu.dat checkReconstructedTmunu_test4.dat

python3 tests/TestOutputFilesCartesian.py
STATUS=$?

if [ $STATUS == 0 ]; then
    echo "All tests passed! :)"
	rm -fr ./check*.dat
	rm -fr tests/*_fixSeed.dat
else
    echo "Tests FAILED :("
    exit 1
fi
