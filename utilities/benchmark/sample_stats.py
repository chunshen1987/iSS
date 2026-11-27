#!/usr/bin/env python3
"""Summarize iSS binary samples and compare two summaries statistically.

  sample_stats.py summarize particle_samples.bin summary.json
  sample_stats.py compare a.json b.json [--show 211,2212]

The binary layout is the one written by use_binary_format = 1
(src/sample_output.h): per event an int count, then per hadron an int pid
and 9 floats (mass, t, x, y, z, E, px, py, pz).

For each species the summary keeps the mean and its standard error of the
multiplicity per event, <pT>, <pT^2>, v2 at |y| < 1, <y^2>, <E>, <px>,
<py>, <px^2> and v1. The comparison prints z = (a - b)/sqrt(da^2 + db^2)
for each observable, and their maximum and rms. For two independent
samples of the same distribution the rms is about 1.
"""
import argparse
import json
import sys

import numpy as np

RECORD = np.dtype([("pid", "<i4"), ("mass", "<f4"), ("t", "<f4"),
                   ("x", "<f4"), ("y", "<f4"), ("z", "<f4"), ("E", "<f4"),
                   ("px", "<f4"), ("py", "<f4"), ("pz", "<f4")])
SPECIES = [211, -211, 111, 321, -321, 2212, -2212, 3122, 3312, 3334, 333,
           223]


def load(filename):
    raw = np.fromfile(filename, dtype=np.uint8)
    offset, hadrons, events = 0, [], []
    iev = 0
    while offset < raw.size:
        count = int(raw[offset:offset + 4].view("<i4")[0])
        offset += 4
        nbytes = count*RECORD.itemsize
        hadrons.append(raw[offset:offset + nbytes].view(RECORD))
        events.append(np.full(count, iev, dtype=np.int32))
        offset += nbytes
        iev += 1
    return np.concatenate(hadrons), np.concatenate(events), iev


def mean_error(values):
    values = np.asarray(values, dtype=np.float64)
    if values.size == 0:
        return (0., 0.)
    return (float(values.mean()), float(values.std()/np.sqrt(values.size)))


def species_stats(hadrons, events, nev, pid):
    selected = hadrons["pid"] == pid
    h = hadrons[selected]
    multiplicity = np.bincount(events[selected], minlength=nev)
    px = h["px"].astype(np.float64)
    py = h["py"].astype(np.float64)
    pz = h["pz"].astype(np.float64)
    E = h["E"].astype(np.float64)
    pT = np.hypot(px, py)
    phi = np.arctan2(py, px)
    y = 0.5*np.log((E + pz)/(E - pz))
    mid = np.abs(y) < 1
    return {
        "N/ev": mean_error(multiplicity),
        "<pT>": mean_error(pT),
        "<pT^2>": mean_error(pT**2),
        "v2(|y|<1)": mean_error(np.cos(2*phi[mid])),
        "<y^2>": mean_error(y**2),
        "<E>": mean_error(E),
        "<px>": mean_error(px),
        "<py>": mean_error(py),
        "<px^2>": mean_error(px**2),
        "v1": mean_error(np.cos(phi)),
    }


def summarize(args):
    hadrons, events, nev = load(args.samples)
    summary = {str(pid): species_stats(hadrons, events, nev, pid)
               for pid in SPECIES}
    summary["all"] = {"N/ev": mean_error(np.bincount(events, minlength=nev))}
    with open(args.summary, "w") as f:
        json.dump(summary, f)
    print(f"{nev} events, {hadrons.size} hadrons -> {args.summary}")


def compare(args):
    with open(args.a) as f:
        a = json.load(f)
    with open(args.b) as f:
        b = json.load(f)
    show = set(args.show.split(",")) if args.show else set(a)
    zs = []
    for species in a:
        line = [f"{species:>6s}"]
        for key, (ma, da) in a[species].items():
            mb, db = b[species][key]
            # species absent from both samples, e.g. resonances after decays
            if not np.all(np.isfinite([ma, da, mb, db])) or (da == 0
                                                              and db == 0):
                continue
            z = (ma - mb)/np.sqrt(da**2 + db**2)
            zs.append(z)
            line.append(f"{key}={ma:.5g}/{mb:.5g}(z={z:+.1f})")
        if species in show:
            print("  ".join(line))
    zs = np.array(zs)
    print(f"{zs.size} observables: max |z| = {np.abs(zs).max():.2f}, "
          f"rms z = {np.sqrt((zs**2).mean()):.2f}, "
          f"|z| > 3: {(np.abs(zs) > 3).sum()}")


def main():
    parser = argparse.ArgumentParser(
        description="Summarize and compare iSS binary samples.")
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("summarize", help="binary samples -> JSON summary")
    p.add_argument("samples")
    p.add_argument("summary")
    p.set_defaults(func=summarize)
    p = sub.add_parser("compare", help="compare two JSON summaries")
    p.add_argument("a")
    p.add_argument("b")
    p.add_argument("--show", default="",
                   help="comma-separated species to print, default all")
    p.set_defaults(func=compare)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    sys.exit(main())
