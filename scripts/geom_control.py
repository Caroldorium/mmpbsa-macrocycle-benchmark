"""Interface-size controls for Zhao et al. 2025 Dataset I.

Computes, from the deposited crystal coordinates only, the descriptors that an
entropy-free MM/PBSA score is largely a smooth function of: buried solvent
accessible surface area, interface atom and contact counts, and peptide size.
No force field, no minimisation.  The point is to ask whether the reported
cross-system correlation with pKd survives a size control.
"""
import sys, os, glob, csv
import numpy as np

# Bondi radii, the set MM/PBSA's LCPO nonpolar term is built on
RADII = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80,
         "F": 1.47, "CL": 1.75, "BR": 1.85, "I": 1.98, "SE": 1.90,
         "ZN": 1.39, "MG": 1.73, "CA": 2.31, "NA": 2.27, "K": 2.75,
         "FE": 2.00, "MN": 2.00, "CU": 1.40, "H": 1.20}
PROBE = 1.4
N_SPHERE = 240


def _sphere(n):
    i = np.arange(n, dtype=float) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5 ** 0.5) * i
    return np.column_stack([np.cos(theta) * np.sin(phi),
                            np.sin(theta) * np.sin(phi), np.cos(phi)])


SPHERE = _sphere(N_SPHERE)


def read_pdb(path, heavy_only=True):
    xyz, elem, res = [], [], []
    with open(path) as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            alt = line[16]
            if alt not in (" ", "A"):
                continue
            e = line[76:78].strip().upper()
            if not e:
                nm = line[12:16].strip()
                e = nm[0] if nm[0].isalpha() else nm[1]
            if heavy_only and e == "H":
                continue
            if e in ("D",):
                continue
            xyz.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            elem.append(e)
            res.append((line[21], line[22:27].strip(), line[17:20].strip()))
    return np.asarray(xyz, float), elem, res


def sasa(xyz, elem, probe=PROBE):
    """Shrake-Rupley, returns per-atom SASA in A^2."""
    r = np.array([RADII.get(e, 1.70) for e in elem]) + probe
    n = len(xyz)
    out = np.zeros(n)
    # neighbour list via a coarse grid on the max radius
    rmax = r.max()
    cell = 2 * rmax
    keys = np.floor(xyz / cell).astype(int)
    buckets = {}
    for i, k in enumerate(map(tuple, keys)):
        buckets.setdefault(k, []).append(i)
    offs = [(a, b, c) for a in (-1, 0, 1) for b in (-1, 0, 1) for c in (-1, 0, 1)]
    for i in range(n):
        kx, ky, kz = keys[i]
        cand = []
        for da, db, dc in offs:
            cand.extend(buckets.get((kx + da, ky + db, kz + dc), ()))
        cand = np.array([j for j in cand if j != i], dtype=int)
        if len(cand):
            d = np.linalg.norm(xyz[cand] - xyz[i], axis=1)
            keep = cand[d < (r[i] + r[cand])]
        else:
            keep = cand
        pts = xyz[i] + r[i] * SPHERE
        if len(keep):
            dd = np.linalg.norm(pts[:, None, :] - xyz[keep][None, :, :], axis=2)
            buried = (dd < r[keep][None, :]).any(axis=1)
            frac = 1.0 - buried.mean()
        else:
            frac = 1.0
        out[i] = 4 * np.pi * r[i] ** 2 * frac
    return out


POLAR = {"N", "O", "S"}


def descriptors(prot_pdb, pep_pdb):
    px, pe, pr = read_pdb(prot_pdb)
    lx, le, lr = read_pdb(pep_pdb)
    cx = np.vstack([px, lx]); ce = pe + le
    s_p = sasa(px, pe).sum()
    s_l = sasa(lx, le).sum()
    s_c_all = sasa(cx, ce)
    s_c = s_c_all.sum()
    dsasa = s_p + s_l - s_c
    # apolar / polar split of the buried area
    s_p_at = sasa(px, pe); s_l_at = sasa(lx, le)
    free = np.concatenate([s_p_at, s_l_at])
    d_at = free - s_c_all
    ap = np.array([e not in POLAR for e in ce])
    # contacts
    d = np.linalg.norm(px[:, None, :] - lx[None, :, :], axis=2)
    c4 = int((d < 4.0).sum()); c5 = int((d < 5.0).sum()); c6 = int((d < 6.0).sum())
    iface_prot = int((d.min(axis=1) < 4.5).sum())
    iface_pep = int((d.min(axis=0) < 4.5).sum())
    iface_res = len({pr[i][:2] for i in np.where(d.min(axis=1) < 4.5)[0]})
    return dict(
        pep_heavy=len(lx), prot_heavy=len(px),
        pep_res=len({r[:2] for r in lr}),
        dSASA=dsasa, dSASA_apolar=float(d_at[ap].sum()), dSASA_polar=float(d_at[~ap].sum()),
        contacts_4A=c4, contacts_5A=c5, contacts_6A=c6,
        iface_prot_atoms=iface_prot, iface_pep_atoms=iface_pep, iface_prot_res=iface_res,
        pep_polar_frac=float(np.mean([e in POLAR for e in le])),
    )


if __name__ == "__main__":
    root, out = sys.argv[1], sys.argv[2]
    rows = []
    for d in sorted(glob.glob(os.path.join(root, "*"))):
        if not os.path.isdir(d):
            continue
        pdb = os.path.basename(d)
        prot = os.path.join(d, f"{pdb}_protein.pdb")
        pep = os.path.join(d, f"{pdb}_CP.pdb")
        if not (os.path.exists(prot) and os.path.exists(pep)):
            print("MISSING", pdb); continue
        rec = descriptors(prot, pep)
        rec["pdb"] = pdb
        rows.append(rec)
        print(pdb, "dSASA=%.0f" % rec["dSASA"], "c4=%d" % rec["contacts_4A"], flush=True)
    cols = ["pdb"] + [k for k in rows[0] if k != "pdb"]
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols); w.writeheader(); w.writerows(rows)
    print("wrote", out, len(rows))
