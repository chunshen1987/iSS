# Review of pull request #6: "Decay through every channel of the decay table"

- PR: https://github.com/chunshen1987/iSS/pull/6
- Base: `XSCAPE`. Head: `IFUSP-Heavy-Ions:fix-n-body-decays` (`b0e62e7`)
- Changes `src/particle_decay.{h,cpp}` and `README.md`. It adds 4- and
  5-body decays sampled uniformly in phase space, opens channels that
  are closed at the pole masses through Breit-Wigner daughter masses,
  and gives unreachable channels' branching ratios to the others.
- Status: **changes requested** (review posted 2026-11-26).

The PR is careful otherwise. The 2- and 3-body paths and their
random-number sequence are unchanged, and its phase-space and
conservation checks are thorough. Two Codex findings were checked
against the code and both hold.

## P1: a resonance with no reachable channel hangs the run

When no channel is reachable, `perform_decays` returns the mother
unchanged. It is still `stable == 0`, so all three callers re-queue it
and decay it again, forever:

- `FSSW.cpp` (`decay_one_primary`, about line 1389)
- `FSSW.cpp` (resonance feed-down, about line 2195)
- `emissionfunction.cpp` (resonance feed-down, about line 3990)

`n_unreachable_` is only reported in the destructor, which a hung run
never reaches.

This can happen with the shipped tables. A broad daughter may be given
masses down to `lowest_mass_ + kMassMargin`. Its own channels with broad
daughters only count as reachable when
`M - kMassMargin > threshold_min`, and `threshold_min` adds
`kMassMargin` for each broad daughter. So there is a window just above
the lower bound where none of its channels is reachable. The fallback in
`sample_daughter_masses` puts masses close to that bound when the mother
is just above threshold.

The PR's reachability rules were re-implemented in a script (iSS was
not run) and evaluated for each broad species at
`lowest_mass_ + kMassMargin`. No species is stuck at its pole mass. At
that lowest mass these have no reachable channel:

- `pdg-urqmd_v3.3+.dat`: h1'(1380) (10333), ω(1420) (50223)
- `pdg-s95pv1.dat`: h1(1170), b1(1235)±/0, a1(1260)⁰, ω(1420), φ(1680),
  f2(2010)

Suggested fix: make the daughters' lower mass bound and the
reachability test agree (for example, a `lowest_mass_` that already
includes the margins its channels need). As a safety net, return a
particle with no reachable channel in a form the callers treat as final.

## P2: an exact-threshold 4- or 5-body channel calls `exit(1)`

Channel selection accepts a channel when `M >= threshold_pole`, and
masses are only sampled when `M < threshold_pole`.
`perform_n_body_decay` exits when `M <= m_sum`, and with pole masses
`m_sum == threshold_pole`. A mother exactly at threshold picks the
channel and then stops the run. This is unlikely, but the selection
should be strict (`M > threshold_pole`), matching the decay routine.

## Also noted

The PR reports two channels in `pdg-s95pv1.dat` that do not conserve
charge: η2(1645) (10225) → K⁺ K̄⁰ and → K⁰ K⁻, branching ratio 0.02 each.
The review asked for that to go in a separate PR.
