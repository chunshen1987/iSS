#!/usr/bin/env bash
set -euo pipefail
curl -fL 'https://www.dropbox.com/s/c4twi65l5ediyjy/NEoS4D_deltafCoeffs.tar.gz?dl=1' -o NEoS4D_deltafCoeffs.tar.gz
tar -xzf NEoS4D_deltafCoeffs.tar.gz
rm NEoS4D_deltafCoeffs.tar.gz
