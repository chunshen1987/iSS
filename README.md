iSS
===

Monte Carlo sampler for particle distribution from Cooper-Frye freeze-out procedure

Unit tests
----------

Build and run the focused C++ tests with CMake and CTest. They use the
bundled doctest header and generate synthetic EoS tables in the build tree.

```sh
cmake -S . -B build -Dunittest=ON
cmake --build build --parallel 2
ctest --test-dir build --output-on-failure
```

The light-fermion momentum diagnostic currently fails against independent
Fermi–Dirac quadrature. Run it directly from the repository root to reproduce
the discrepancy.

```sh
build/src/unittest_numerics.e --test-case='Light fermion momentum*'
```

Resonance decays (`perform_decays = 1`)
---------------------------------------

A decay channel is drawn with the table's branching ratios, normalised over
the channels the mother can reach at its pole mass. A channel that is closed
at the daughters' pole masses is reached by giving its broad daughters (width
above 0.1 MeV) masses from their Breit-Wigner distributions, each no lower
than the mass at which all of its own channels stay open, so that the
daughter keeps the table's branching ratios; a channel that cannot be reached
that way (e.g. a0(980), f0(980) -> K Kbar in `pdg-s95pv1.dat`) is skipped and
its share goes to the others. A particle that is below all of its channels
is kept in the final state, and the number of such particles is reported at
the end of the run; with the shipped tables this does not happen, because a
daughter is never given a mass at which it cannot decay. A negative number
of daughters (how
`pdg-urqmd_v3.3+.dat` and `pdg-s95pv1.dat` mark closed channels) is read as
its absolute value. Two- and three-body channels use the original routines and
random-number sequence; four- and five-body ones are uniform in phase space
(Raubold-Lynch).
