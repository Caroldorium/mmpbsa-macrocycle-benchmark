#!/usr/bin/env python3
"""Pipeline schematic, written for readers with no computational background."""
import os, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

OUT = os.path.expanduser('~/Downloads/lpta_methods_20260927/fig1_pipeline.png')
NAVY, CRIM, GREY, INK = '#002A5C', '#B51F58', '#5E6873', '#1F2933'

fig, ax = plt.subplots(figsize=(13.2, 6.4))
ax.set_xlim(0, 13.2); ax.set_ylim(0, 6.4); ax.axis('off')

def box(x, y, w, h, title, lines, edge=NAVY, fill='#F4F7FA'):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.07",
                                linewidth=1.5, edgecolor=edge, facecolor=fill))
    ax.text(x + w/2, y + h - 0.30, title, ha='center', va='top',
            fontsize=10.5, fontweight='bold', color=edge)
    for i, t in enumerate(lines):
        ax.text(x + w/2, y + h - 0.66 - i*0.29, t, ha='center', va='top',
                fontsize=8.6, color=INK)

def arrow(x1, y1, x2, y2):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle='-|>',
                                 mutation_scale=15, linewidth=1.4, color=GREY))

box(0.15, 3.55, 2.5, 2.65, "1  Structure",
    ["Peptide docked or", "modelled in the", "protein's binding site", "", "One 3-D pose"])
box(3.0, 3.55, 2.5, 2.65, "2  Chemistry",
    ["Force field assigns", "charges and bonded", "terms to every atom", "",
     "Non-standard residues", "parameterised here"])
box(5.85, 3.55, 2.5, 2.65, "3  Relax",
    ["Brief energy", "minimisation removes", "bad contacts", "",
     "Seconds, not", "a simulation"])
box(8.7, 3.55, 2.5, 2.65, "4  Score",
    ["Interaction energy", "minus desolvation", "cost (MM/PBSA)", "",
     "One number per", "peptide"])
box(11.55, 3.55, 1.5, 2.65, "5  Rank", ["Sort", "designs", "", "Pick what", "to make"],
    edge=CRIM, fill='#FDF2F6')
for x in (2.72, 5.57, 8.42, 11.27):
    arrow(x, 4.85, x + 0.22, 4.85)

ax.text(6.6, 3.22, "Validation — the part that decides whether any of this is trustworthy",
        ha='center', fontsize=11, fontweight='bold', color=CRIM)

box(0.15, 0.25, 4.0, 2.65, "Benchmark",
    ["234 macrocyclic peptides across 5 proteins", "with measured binding affinities",
     "(Awoonor-Williams et al. 2026)", "",
     "Published FEP+ values give a ceiling", "for what rigorous physics achieves"])
box(4.5, 0.25, 4.1, 2.65, "The trap",
    ["Bigger peptides bind tighter because", "they bury more surface.",
     "Any method looks good if size varies.", "",
     "So we compare only peptides with the", "SAME heavy-atom count"], edge=CRIM, fill='#FDF2F6')
box(8.95, 0.25, 4.1, 2.65, "The measure",
    ["Give it two peptides, ask which binds", "tighter. A coin gets 50%.", "",
     "MM/PBSA: 60.1%   FEP+: 84.1%", "heavy-atom count: no information",
     "(the sizes are equal)"], edge=CRIM, fill='#FDF2F6')
arrow(4.2, 1.55, 4.45, 1.55); arrow(8.65, 1.55, 8.9, 1.55)
ax.add_patch(FancyArrowPatch((9.95, 3.5), (6.6, 3.05), arrowstyle='-|>', mutation_scale=15,
                             linewidth=1.4, color=GREY, connectionstyle="arc3,rad=0.25"))
fig.savefig(OUT, dpi=300, bbox_inches='tight', facecolor='white')
print('wrote', OUT)
