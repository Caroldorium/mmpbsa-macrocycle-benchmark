#!/usr/bin/env python3
"""Split the connected edge table into one edge-directory list per family."""
import csv, collections, os

rows = list(csv.DictReader(open("edge_list_connected.csv")))
by = collections.defaultdict(list)
for r in rows:
    fam = r["family"].replace("/", "_").replace(" ", "_")
    by[fam].append(r)

for fam, g in sorted(by.items()):
    path = "list_%s.txt" % fam
    with open(path, "w") as fh:
        for r in g:
            fh.write("%s__%s\n" % (r["a"], r["b"]))
    nconn = sum(1 for r in g if str(r.get("connector", "")).lower() == "true")
    missing = sum(1 for r in g
                  if not os.path.exists("edges/%s__%s/complex.prmtop" % (r["a"], r["b"])))
    print("%-14s %4d edges (%2d connectors)  missing topology: %d"
          % (fam, len(g), nconn, missing))
print("total %d edges" % len(rows))
