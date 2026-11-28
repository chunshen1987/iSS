iSS
===

Monte Carlo sampler for particle distribution from Cooper-Frye freeze-out procedure

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

Direct UrQMD output
-------------------

Set `use_UrQMD_format = 1` to write sampled events in UrQMD's restart
format to `fort.14` in the working directory. Pass that file to UrQMD as
`ftn10` in place of the file normally produced by `osc2u`. The option takes
priority over OSCAR, gzip, and binary event output and works with
`MC_sampling = 1` through `4`. It holds sampled events in memory while
writing and omits empty samples.

Use `afterburner_type = 1` to select the UrQMD particle table. The writer
reports an error if a sampled PDG ID has no UrQMD mapping. Update the
third line of `iSS_tables/OSCAR_header.txt` for the projectile, target,
reference frame, and lab kinetic energy per nucleon. The shipped header
describes Au+Au at 100 GeV per nucleon in the equal-speed frame. As in
`osc2u`, the impact parameter is zero and the frame setting describes the
sampled momenta. The writer back-propagates the hadrons to time zero and
records their original hadronization times for UrQMD.
