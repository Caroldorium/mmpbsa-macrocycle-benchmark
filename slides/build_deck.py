#!/usr/bin/env python3
"""Rebuild the collaboration-meeting deck around the Zha25 single-structure protocol.

Starts from the existing deck rather than the bare template, so the Brown Lab master, the LPS
pathway slide and the embedded movie survive untouched -- python-pptx cannot copy those cleanly
between files.  Unwanted slides are removed, kept slides retitled, new slides appended, then the
whole set reordered.  Supplementary slides are marked hidden.
"""
import copy, pathlib, re
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from pptx.opc.packuri import PackURI

SRC = '/Users/daddycl/Downloads/202601019-collabMeeting.pptx'
OUT = '/Users/daddycl/Downloads/20261019-collabMeeting-v2.pptx'
REPO = pathlib.Path('/Users/daddycl/Desktop/Lab_computation_results.nosync/projects/'
                    'mmpbsa-macrocycle-benchmark')
FIG = REPO / 'figures'

NAVY = RGBColor(0x00, 0x2A, 0x5C)
INK = RGBColor(0x1F, 0x29, 0x33)
GREY = RGBColor(0x5E, 0x68, 0x73)
ACC = RGBColor(0xB5, 0x1F, 0x58)

SW, SH = Inches(13.333), Inches(7.5)
BODY = (Inches(0.47), Inches(1.21), Inches(12.39), Inches(5.46))

prs = Presentation(SRC)
BLANK = next(l for l in prs.slide_masters[0].slide_layouts if l.name == '1_Blank')
sldIdLst = prs.slides._sldIdLst


# ---------------------------------------------------------------- helpers
def style(run, size=18, color=INK, bold=False, italic=False):
    run.font.name = 'Arial'
    run.font.size = Pt(size)
    run.font.bold = bold
    run.font.italic = italic
    run.font.color.rgb = color


def set_title(slide, text, size=26):
    ph = slide.shapes.title
    tf = ph.text_frame
    tf.clear()
    p = tf.paragraphs[0]
    style(p.add_run(), size)
    p.runs[0].text = text
    style(p.runs[0], size, NAVY, bold=True)
    tf.word_wrap = True
    return ph


def footer(slide, text, size=8):
    tb = slide.shapes.add_textbox(Inches(0.34), Inches(7.02), Inches(12.6), Inches(0.38))
    tf = tb.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]
    r = p.add_run(); r.text = text
    style(r, size, GREY)
    return tb


def drop_body_placeholder(slide):
    for sh in list(slide.placeholders):
        if sh.placeholder_format.idx != 0:
            sh._element.getparent().remove(sh._element)


_next = [0]


def new_slide():
    """Add a blank slide, giving its part a partname above every existing one.

    python-pptx numbers a new slide part by the count of slides, not by the first free index.
    After deleting slides that leaves gaps, so a new part can collide with a kept one and the
    saved package ends up with duplicate zip entries.  Allocating above the current maximum
    avoids it without renaming anything that already exists.
    """
    if not _next[0]:
        used = [int(re.search(r'slide(\d+)\.xml$', sl.part.partname).group(1))
                for sl in prs.slides]
        _next[0] = max(used) + 1
    s = prs.slides.add_slide(BLANK)
    s.part.partname = PackURI('/ppt/slides/slide%d.xml' % _next[0])
    _next[0] += 1
    drop_body_placeholder(s)
    return s


def fit(img, box, top_pad=0.0):
    """Centre an image inside box=(l,t,w,h) in EMU, preserving aspect."""
    from PIL import Image
    iw, ih = Image.open(img).size
    l, t, w, h = box
    t = t + Emu(int(top_pad * 914400))
    h = h - Emu(int(top_pad * 914400))
    scale = min(w / iw, h / ih)
    nw, nh = int(iw * scale), int(ih * scale)
    return int(l + (w - nw) / 2), int(t + (h - nh) / 2), nw, nh


def fig_slide(title, img, note=None, cite=None, title_size=26):
    s = new_slide()
    set_title(s, title, title_size)
    pad = 0.0
    if note:
        tb = s.shapes.add_textbox(Inches(0.47), Inches(1.06), Inches(12.39), Inches(0.44))
        tf = tb.text_frame; tf.word_wrap = True
        r = tf.paragraphs[0].add_run(); r.text = note
        style(r, 13, ACC, bold=True)
        pad = 0.46
    l, t, w, h = fit(str(img), BODY, top_pad=pad)
    s.shapes.add_picture(str(img), l, t, w, h)
    if cite:
        footer(s, cite)
    return s


def bullets(slide, items, left=0.47, top=1.25, width=12.39, height=5.4, base=18, gap=10):
    tb = slide.shapes.add_textbox(Inches(left), Inches(top), Inches(width), Inches(height))
    tf = tb.text_frame; tf.word_wrap = True
    first = True
    for lvl, text, *rest in items:
        sz = rest[0] if rest else (base if lvl == 0 else base - 2)
        p = tf.paragraphs[0] if first else tf.add_paragraph()
        first = False
        p.level = lvl
        p.space_after = Pt(gap if lvl == 0 else gap - 4)
        r = p.add_run(); r.text = text
        style(r, sz, INK if lvl == 0 else GREY, bold=(lvl == 0))
    return tb


def stat(slide, x, y, value, label, w=3.0, vsize=27, lsize=11.5, color=ACC):
    tb = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(1.15))
    tf = tb.text_frame; tf.word_wrap = True
    p = tf.paragraphs[0]; r = p.add_run(); r.text = value
    style(r, vsize, color, bold=True)
    p2 = tf.add_paragraph(); r2 = p2.add_run(); r2.text = label
    style(r2, lsize, GREY)
    return tb


def table(slide, rows, left, top, width, height, fs=12, hdr=True, widths=None):
    nr, nc = len(rows), len(rows[0])
    gf = slide.shapes.add_table(nr, nc, Inches(left), Inches(top), Inches(width), Inches(height))
    tbl = gf.table
    if widths:
        for i, cw in enumerate(widths):
            tbl.columns[i].width = Inches(cw)
    for ri, row in enumerate(rows):
        for ci, val in enumerate(row):
            cell = tbl.cell(ri, ci)
            cell.text = ''
            p = cell.text_frame.paragraphs[0]
            r = p.add_run(); r.text = str(val)
            bold = (hdr and ri == 0) or str(val).startswith('**')
            if r.text.startswith('**'):
                r.text = r.text.strip('*')
            style(r, fs, NAVY if (hdr and ri == 0) else INK, bold=bold)
            if ci > 0:
                p.alignment = PP_ALIGN.RIGHT
    return gf


def section(title, sub=None):
    s = new_slide()
    t = s.shapes.title
    t._element.getparent().remove(t._element)
    tb = s.shapes.add_textbox(Inches(0.8), Inches(2.7), Inches(11.7), Inches(2.0))
    tf = tb.text_frame; tf.word_wrap = True
    p = tf.paragraphs[0]; r = p.add_run(); r.text = title
    style(r, 36, NAVY, bold=True)
    if sub:
        p2 = tf.add_paragraph(); r2 = p2.add_run(); r2.text = sub
        style(r2, 16, GREY)
    return s


# ------------------------------------------------- 1. keep / delete originals
orig = list(prs.slides)
KEEP = {1: 'title', 2: 'lps', 5: 'protocol', 8: 'reagents', 16: 'descriptors',
        19: 'mlchart', 20: 'fullsweep', 21: 'ablation', 22: 'anchors', 23: 'krasmetrics'}
kept = {name: orig[i - 1] for i, name in KEEP.items()}
objective = orig[16].placeholders[1].text_frame.text.strip()

ids = list(sldIdLst)
for i in range(len(orig), 0, -1):
    if i not in KEEP:
        prs.part.drop_rel(ids[i - 1].rId)
        sldIdLst.remove(ids[i - 1])
print(f'kept {len(prs.slides)} original slides')


# ------------------------------------------------- 2. edit the kept slides
# --- title slide
t = kept['title']
for ph in t.placeholders:
    if ph.placeholder_format.idx == 1:
        tf = ph.text_frame; tf.clear()
        p = tf.paragraphs[0]
        r = p.add_run(); r.text = 'Collaboration Meeting  |  19 October 2026'
        style(r, 20, INK)
        p2 = tf.add_paragraph()
        r2 = p2.add_run()
        r2.text = ('Carlos Liu  |  PhD Student  |  Brown Lab  |  '
                   'Leslie Dan Faculty of Pharmacy  |  University of Toronto')
        style(r2, 18, GREY)

# --- protocol slide (was 5): retitle, replace body text, keep the movie
pr_s = kept['protocol']
set_title(pr_s, 'The oracle is a single minimised structure scored in seconds, '
                'not a trajectory', 22)
for sh in list(pr_s.shapes):
    if sh.has_text_frame and not sh.is_placeholder:
        sh._element.getparent().remove(sh._element)
for sh in list(pr_s.placeholders):
    if sh.placeholder_format.idx != 0:
        sh._element.getparent().remove(sh._element)
bullets(pr_s, [
    (0, 'Build and minimise one complex, then score it'),
    (1, 'tleap parameterisation; ncAA partial charges and force-field parameters generated '
        'in-house, geometry and chirality checked'),
    (1, 'OpenMM minimisation, OBC1 — no explicit-solvent trajectory'),
    (1, 'MM/PBSA with an interior dielectric of 2 (pb_2), ff03'),
    (0, 'Cost is the point'),
    (1, 'FEP+ costs GPU-days per compound; this costs 1–10 s — roughly three orders '
        'of magnitude cheaper'),
    (1, 'Section 3.7 shows the expensive alternative buys nothing (slide 9)'),
], left=0.47, top=1.30, width=5.6, base=16, gap=9)
footer(pr_s, 'Zhao et al. 2025 single-structure MM/PBSA protocol, reimplemented · '
             'movie: canonical vs ncAA-swapped complex, illustrative only')

# --- ML chart (was 19)
ml = kept['mlchart']
set_title(ml, 'Morgan + PLEC reaches ρ 0.563 — two-thirds of the FEP+ ceiling', 24)
footer(ml, 'Equal-family mean Spearman ρ over 5 target families, 140 ten-anchor starts, '
           '7 contexts · FEP+ here is 0.838, scored on the held-out query compounds of each '
           'anchor split — not the same estimand as the 0.831 over all 235 compounds used on '
           'slides 6–8, so the two must not be read off one axis')

# --- anchor curve (was 22)
an = kept['anchors']
set_title(an, 'Twenty anchors buys ρ 0.735; five buys 0.626', 24)
for sh in list(an.placeholders):
    if sh.placeholder_format.idx != 0:
        sh._element.getparent().remove(sh._element)
footer(an, 'Random anchor selection, 3-family cohort (6XIC / 8BPV / 8QSO) — a narrower set '
           'than the 5-family benchmark on the preceding slide · maxmin selection on the same '
           'data runs 0.688 → 0.749')

# --- supplementary: retitle the kept ML/MD slides
set_title(kept['descriptors'], 'Four descriptor families tested on the same anchor protocol', 24)
for sh in kept['descriptors'].shapes:
    if sh.has_table:
        for row in sh.table.rows:
            c = row.cells[3]
            if '0.331' in c.text:
                c.text_frame.paragraphs[0].runs[0].text = '0.524–0.540 alone · 0.563 + Morgan'
footer(kept['descriptors'],
       'CORRECTED: the PLEC-alone figure previously shown (0.331) came from a leave-family-out '
       'benchmark over 225 compounds, not the 140-anchor protocol stated above. Under this '
       'protocol PLEC alone is 0.524 (Bayesian linear) or 0.540 (pairwise logistic).')
set_title(kept['fullsweep'], 'Thirteen methods on the same five families — KRAS breaks all of them', 23)
set_title(kept['ablation'], 'Added descriptors do not improve the Morgan + PLEC ranker', 24)
set_title(kept['krasmetrics'], 'On KRAS, an H-bond count out-ranks the score itself', 24)
for _sh in kept['krasmetrics'].shapes:                 # fix a typo carried from the source deck
    if _sh.has_text_frame:
        for _p in _sh.text_frame.paragraphs:
            for _r in _p.runs:
                if 'Prehaps' in _r.text:
                    _r.text = _r.text.replace('Prehaps', 'Perhaps')
footer(kept['krasmetrics'],
       'KRAS same-20, 4 replicas × 10 ns, 2–10 ns window · every metric is tagged '
       'EXPLORATORY_POST_HOC in the source, and a same-week feature-fusion analysis ranks '
       'H-bonds and native contacts last by F-score — treat as a hypothesis, not a result')
set_title(kept['reagents'], 'The earlier three-swap purchase, made before the result reversed', 23)
_t = kept['reagents'].shapes.title          # source box starts at x=0.09 and runs off the edge
_t.left, _t.top, _t.width = 0, 0, Inches(13.333)
for sh in list(kept['reagents'].shapes):
    if sh.has_text_frame and 'Is MD validated' in sh.text_frame.text:
        sh.text_frame.paragraphs[0].runs[0].text = ''
footer(kept['reagents'],
       'S03 = M1 Q6X/QAY/Q2W, on the RF-parent branch — a different series from the '
       'derivative-7 purchase on the previous slide. Reagents were bought when the 3 × 10 ns '
       'swap-minus-parent verdict was INCONCLUSIVE (+1.72 kcal/mol, CI straddles zero); the '
       '4 × 250 ns rerun returns UNFAVORABLE (+6.21, CI [+3.06, +8.86]). The size-matched '
       'control on the previous slide is the check that was not available at the time.')


# ------------------------------------------------- 3. new core slides
new = {}

new['pipeline'] = fig_slide(
    'Our design pipeline, end to end',
    FIG / 'fig1_pipeline.png',
    cite='Generative backbone design → sequence design → ncAA substitution → MM/PBSA triage '
         '→ synthesis and assay')

s = new_slide()
set_title(s, 'Two choices decide whether a binding benchmark means anything', 24)
bullets(s, [
    (0, '1.  Average within a receptor, never across receptors'),
    (1, 'Pooling compounds from different crystal forms rewards a method for separating two '
        'groups, not for ranking inside either'),
    (1, 'On KRAS this reverses the sign: pooled +0.225, macro over forms −0.130 (supplementary)'),
    (0, '2.  Compare compounds of the same size'),
    (1, 'Experimental affinity tracks heavy-atom count, and the score tracks it harder — so an '
        'unmatched comparison partly measures size'),
    (1, 'A heavy-atom count reaches ρ 0.487 on its own. Any method must beat that, not chance'),
    (0, 'Both were fixed before the results were read, and both are stated in the paper'),
], top=1.45, base=17, gap=13)
footer(s, 'Decisive pairs only (|ΔΔG| ≥ 0.5 kcal/mol) · 95% intervals bootstrap over compounds, '
          'not pairs, because pairs drawn from n compounds are not independent')
new['design'] = s

new['pairwise'] = fig_slide(
    'Among equally sized peptides, the score picks the tighter binder above chance',
    FIG / 'size_matched_pairwise_accuracy.png',
    note='Within a three-atom window: 0.624 [0.555, 0.689] against chance 0.500 — '
         'FEP+ reaches 0.839',
    cite='A size descriptor cannot answer at all when the sizes are equal — which is the '
         'point: it is a proxy for buried interface area, not a rival scoring function. '
         'Within three atoms it reaches 0.602, so that is the bar, not chance.',
    title_size=23)

new['parity'] = fig_slide(
    'Per target: ρ +0.545 against an FEP+ ceiling of +0.831',
    FIG / 'zha25_fepplus_parity_5families.png',
    cite='234 compounds, 5 targets · row 2 applies a per-family linear calibration; the fitted '
         'slopes of 0.099–0.208 mean the raw score moves 5–10 kcal/mol per 1 kcal/mol of real '
         'affinity — it ranks, it is not a ΔG')

new['confound'] = fig_slide(
    'Experiment tracks ligand size; the score tracks it harder',
    FIG / 'fig2_size_confound.png',
    note='Experiment ρ −0.43 · MM/PBSA ρ −0.55 — which is why the headline result is '
         'size-matched',
    cite='Mean within-target Spearman across the five families')

new['sampling'] = fig_slide(
    'Four 10 ns replicas buy nothing — the noise exceeds the signal',
    FIG / 'sampling_noise.png',
    note='Replica-to-replica SD 6.81 kcal/mol against a between-compound SD of 5.74',
    cite='This is the measurement behind scoring a single minimised structure instead of a '
         'trajectory: the noise from which replicas you happened to run is larger than the '
         'entire signal distinguishing one compound from another.')

new['rosetta'] = fig_slide(
    'Rosetta cannot serve as the reference standard',
    FIG / 'rosetta_vs_experiment.png',
    note='Macro ρ −0.149 over 239 compounds, wrong-signed on four of five targets',
    cite='A reference that does not itself track experiment cannot be used to decide which '
         'force field or protocol is better — stated in the paper as a caution on the LptA '
         'force-field comparison')

new['so3lr'] = fig_slide(
    'A machine-learned force field lands below a heavy-atom count',
    FIG / 'so3lr_comparison.png',
    note='SO3LR-SF 0.292 · heavy-atom null 0.487 · MM/PBSA 0.545',
    cite='SO3LR-SF has no solvation model — it is a gas-phase interaction energy. Adding the '
         'PB desolvation term to our own gas-phase term moves macro ρ from 0.394 to 0.545, the '
         'largest single lever in the comparison.')

new['ff'] = fig_slide(
    'Force-field choice does not transfer between datasets',
    FIG / 'zha25_lpta_forcefield.png',
    cite='ff03 wins four of five benchmark targets (0.545 vs 0.495 for ff19SB); on the LptA '
         'designs the ordering reverses. The LptA comparison is scored against Rosetta, which '
         'the previous slide shows does not track experiment — so this is reported as a '
         'caution, not a recommendation.')

new['lpta212'] = fig_slide(
    '212 designs scored in three branches — pooling them gives an artefact',
    FIG / 'zha25_lpta_212.png',
    note='Pooled ρ −0.316 is an artefact of three branches occupying different score ranges; '
         'within branch, +0.671 / +0.669 / +0.484',
    cite='20 canonical RFpeptide designs, 96 ncAA swaps on those parents, 96 ncAA swaps on a '
         'thanatin derivative-7 scaffold · no design has an experimental affinity, so this is '
         'agreement between two computational methods, not accuracy',
    title_size=24)

new['lptacanon'] = fig_slide(
    'On the 20 canonical designs, agreement of +0.671 falls to +0.221 under size control',
    FIG / 'zha25_lpta_canonical.png',
    note='The 95% interval [−0.33, +0.69] crosses zero — and panel C shows minimisation depth '
         'is irrelevant here (artefact/signal 0.04)',
    cite='Head-to-tail cyclic, canonical residues, inside the source protocol’s stated '
         'applicability domain',
    title_size=22)

new['swap'] = fig_slide(
    'The same confound sits inside our own swaps',
    FIG / 'zha25_lpta_parent_swap.png',
    note='75 of 96 swaps beat their parent — and bury 70 Å² more interface. '
         'ρ +0.209 falls to +0.053 once buried area is controlled',
    cite='96 ncAA swaps each matched to the canonical parent it was built from, differing at '
         'one residue · the collapse holds in ff03, ff14SB and ff19SB')

s = new_slide()
set_title(s, 'Ten compounds, ranking frozen before synthesis — and what is still unmeasured', 23)
table(s, [
    ['', 'registered design'],
    ['compounds', '10 de novo LptA designs, 2 canonical + 8 ncAA'],
    ['size window', '149–153 heavy atoms — essentially size-matched'],
    ['predicted range', '−104.0 to −69.8 kcal/mol, spanning the range rather than '
                        'picking winners'],
    ['ranking', 'fixed at selection, will not be revised'],
], left=0.6, top=1.28, width=7.3, height=1.9, fs=11.5, widths=[1.7, 5.6])
table(s, [
    ['n', 'power at ρ = 0.30', 'ρ = 0.50', 'ρ = 0.75'],
    ['10', '0.12', '0.28', '0.67'],
    ['20', '0.22', '0.58', '0.96'],
    ['40', '0.43', '0.88', '1.00'],
], left=8.3, top=1.28, width=4.4, height=1.6, fs=10.5)
bullets(s, [
    (0, 'Ten compounds cannot confirm the ranking — 12% power at the measured effect size, '
        'stated in advance rather than discovered afterwards'),
    (0, 'Ten compounds can do two things'),
    (1, 'Detect a catastrophic inversion — 66% power against ρ = −0.75, a live concern '
        'because exactly that happened on the external target'),
    (1, 'Supply ten measured anchors for within-series calibration'),
    (0, 'No LptA design has a measured affinity today. Nothing here is a validation.'),
], top=3.45, base=14, gap=6)
footer(s, 'Registered as a sanity check and anchor acquisition, not confirmatory validation · '
          'confirmatory validation would need roughly 40 compounds')
new['prospective'] = s

s = new_slide()
set_title(s, 'Ten building blocks bought for derivative-7 — and how they rank in our own panel', 21)
rows = [['reagent (ncAA code)', 'site', 'USD', 'rank /96', 'size-controlled']]
rows += [
    ['Fmoc-L-Asp(OFm)-OH  (Q2B)', '3', '130.25', '1', '2'],
    ['Fmoc-L-Cys(STmp)-OH  (Q6W)', '3', '183.00', '3', '**67'],
    ['Fmoc-L-Dap(Ac)-OH  (QF2)', '14', '29.00', '6', '9'],
    ['Fmoc-L-3,3-diphenyl-Ala  (Q8S)', '14', '16.75', '8', '7'],
    ['Fmoc-L-Dab(Ac)-OH  (QEW)', '9', '103.25', '11', '5'],
    ['Fmoc-L-Phe(2,4-F2)-OH  (QC3)', '5', '64.60', '14', '14'],
    ['Fmoc-O-benzyl-L-Tyr  (Z1T)', '5', '20.00', '15', '18'],
    ['Fmoc-L-Phe(2-Me)-OH  (Q7I)', '5', '11.50', '16', '26'],
    ['Fmoc-L-Lys(propionyl)-OH  (Q8P)', '14', '43.00', '23', '11'],
    ['Fmoc-S-4-MeO-Bzl-L-Cys  (Q99)', '5', '2.50', '27', '**54'],
]
table(s, rows, left=0.45, top=1.22, width=8.55, height=4.15, fs=10.5,
      widths=[3.5, 0.75, 1.1, 1.3, 1.9])
stat(s, 9.35, 1.30, '$603.85', 'ten building blocks across four sites (3, 5, 9, 14)', w=3.5)
stat(s, 9.35, 2.62, '12 of 96', 'median predicted rank — 48 if chosen at random', w=3.5)
stat(s, 9.35, 3.94, '8 of 10', 'still inside the top 30 once buried area is removed',
     w=3.5, color=RGBColor(0x2F, 0x6B, 0x3C))
tb = s.shapes.add_textbox(Inches(0.45), Inches(5.52), Inches(12.4), Inches(1.35))
tf = tb.text_frame; tf.word_wrap = True
for i, (txt, col, bold) in enumerate([
        ('The two that move are the ones to watch. Fmoc-L-Cys(STmp)-OH — the most expensive '
         'item at $183 — falls from 3rd to 67th once buried interface area is controlled, and '
         'Fmoc-S-4-MeO-Bzl-L-Cys from 27th to 54th. Their apparent ranking was largely size.',
         INK, False),
        ('The other eight hold their position, so the selection as a whole is not the size '
         'confound of slide 17 reappearing — this is the check being applied before spending, '
         'rather than after.', INK, False)]):
    p_ = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
    p_.space_after = Pt(6)
    r_ = p_.add_run(); r_.text = txt
    style(r_, 12.5, col, bold=bold)
footer(s, 'Zha25 MM/PBSA on the 96 single-swap derivative-7 designs, converged scores · '
          'size control is the rank of each design\'s score after regressing out buried dSASA '
          '· no design in this panel has a measured affinity, so these are predicted ranks, '
          'not predicted potencies')
new['purchase'] = s


# ------------------------------------------------- 4. supplementary slides
supp = []

supp.append(section('Supplementary',
                    'Benchmark detail · rejected approaches · the retired long-MD track · '
                    'ML detail · provenance'))

s = new_slide()
set_title(s, 'Per-target performance on the benchmark cohort', 24)
table(s, [
    ['target', 'n', 'MM/PBSA r', 'R²', 'MM/PBSA ρ', 'FEP+ ρ', 'slope', 'offset'],
    ['Cyclin A/B', '21', '+0.609', '0.371', '+0.591', '+0.791', '0.202', '−59.4'],
    ['JAK2', '48', '+0.588', '0.345', '+0.595', '+0.788', '0.103', '−44.8'],
    ['KRAS', '45', '+0.300', '0.090', '+0.225', '+0.786', '0.125', '−74.8'],
    ['MCL-1', '57', '+0.693', '0.480', '+0.728', '+0.859', '0.099', '−69.0'],
    ['PCSK9', '63', '+0.551', '0.304', '+0.584', '+0.929', '0.208', '−64.9'],
    ['**macro', '**234', '**+0.548', '**0.318', '**+0.545', '**+0.831', '—', '—'],
], left=0.7, top=1.6, width=11.9, height=3.4, fs=13)
footer(s, 'Pearson r and Spearman ρ computed within receptor group · "slope" is the per-target '
          'regression slope of score on experiment; its reciprocal gives the factor by which '
          'ΔΔG spread is inflated')
supp.append(s)

s = new_slide()
set_title(s, 'Pairwise accuracy on decisive pairs, with intervals', 24)
table(s, [
    ['target', 'identical size', '95% CI', 'within 3 atoms', '95% CI'],
    ['Cyclin A/B', '0.533', '[0.053, 1.000]', '0.674', '[0.454, 0.863]'],
    ['JAK2', '0.637', '[0.397, 0.831]', '0.598', '[0.437, 0.732]'],
    ['KRAS', '0.518', '[0.289, 0.727]', '0.496', '[0.339, 0.664]'],
    ['MCL-1', '0.549', '[0.328, 0.757]', '0.594', '[0.444, 0.733]'],
    ['PCSK9', '**0.765', '[0.613, 0.882]', '**0.757', '[0.671, 0.830]'],
    ['**macro, MM/PBSA', '**0.601', '[0.472, 0.720]', '**0.624', '[0.555, 0.689]'],
    ['macro, FEP+', '0.841', '[0.730, 0.912]', '0.839', '[0.782, 0.886]'],
], left=0.9, top=1.5, width=11.5, height=3.8, fs=12.5)
footer(s, 'Chance = 0.500 · decisive pairs are |ΔΔG| ≥ 0.5 kcal/mol · intervals are 95% '
          'bootstrap over COMPOUNDS; a pair-level bootstrap understates them badly '
          '(it gives [0.530, 0.693] where the compound-level interval is [0.472, 0.720]) · '
          'only PCSK9 is individually above chance')
supp.append(s)

s = new_slide()
set_title(s, 'Pooling across crystal forms reverses the KRAS result', 24)
table(s, [
    ['KRAS', 'n', 'MM/PBSA ρ', 'FEP+ ρ'],
    ['within form 7YUZ', '17', '−0.307', '+0.576'],
    ['within form 8JJS', '28', '+0.046', '+0.782'],
    ['**macro over forms', '45', '**−0.130', '+0.679'],
    ['pooled over both', '45', '**+0.225', '+0.786'],
], left=1.6, top=1.45, width=10.0, height=1.9, fs=13)
bullets(s, [
    (0, 'The mechanism is explicit'),
    (1, 'The score is offset by −3.05 kcal/mol between the two crystal forms, in the same '
        'direction as a −1.15 kcal/mol offset in the experimental labels'),
    (1, 'So the pooled +0.225 is credit for separating two groups, not for ranking within '
        'either; FEP+ is far less affected (+0.679 vs +0.786)'),
    (0, 'This happens inside what is conventionally reported as a within-family metric'),
], top=3.65, base=14, gap=8)
supp.append(s)

s = new_slide()
set_title(s, 'Why a single-target number and the cohort number differ', 24)
bullets(s, [
    (0, 'Both measurements are real; they measure different things'),
    (1, 'ρ = −0.101 is the KRAS same-20 panel — one target, twenty ligands'),
    (1, 'ρ = +0.545 is the macro average over five targets and 234 compounds'),
    (0, 'KRAS is the hardest of the five targets for this method'),
    (1, 'Its own macro value over the two crystal forms is −0.130 — within noise of the −0.101'),
    (1, 'Its size-matched pairwise accuracy, 0.496, is the only one of five that sits at chance'),
    (0, 'The cohort number is the one that generalises'),
    (1, 'Four of five targets are positive; PCSK9 reaches 0.757 pairwise accuracy '
        '[0.671, 0.830]'),
    (1, 'A single target cannot establish or refute a method — which is why the benchmark '
        'is five'),
], top=1.45, base=16, gap=11)
footer(s, 'Table 3 and Table 5 of the manuscript · KRAS same-20 figures from '
          'kras-md-metric-correlation-20260917')
supp.append(s)

s = new_slide()
set_title(s, 'The score ranks, but it is not a ΔG', 24)
stat(s, 0.8, 1.45, '0.099 – 0.208', 'per-family regression slope of score on experiment')
stat(s, 5.2, 1.45, '5 – 10 kcal/mol', 'of raw score per 1 kcal/mol of real affinity')
stat(s, 9.4, 1.45, '1.11 vs 0.85', 'macro RMSE after calibration, against FEP+')
bullets(s, [
    (0, 'What this means in practice'),
    (1, 'The raw number cannot be quoted as a binding free energy — its scale is roughly five '
        'to ten times the experimental range'),
    (1, 'A per-family linear calibration makes the RMSE interpretable, but needs measured '
        'anchors in that family to fit — which is why the prospective test is framed as '
        'anchor acquisition as well as a sanity check'),
    (0, 'Use it to order compounds within a series. Do not use it to predict affinity.'),
], top=3.05, base=15, gap=9)
supp.append(s)

s = new_slide()
set_title(s, 'How the intervals are computed', 24)
bullets(s, [
    (0, 'Resample compounds, not pairs'),
    (1, 'Pairs drawn from n compounds share compounds and are not independent; a pair-level '
        'bootstrap badly understates the uncertainty'),
    (1, 'Exact-match band: pair-level gives [0.530, 0.693]; compound-level gives '
        '[0.472, 0.720] — which no longer excludes chance'),
    (0, 'Estimator: macro, not pooled'),
    (1, 'Pairs are formed WITHIN a target, then accuracies averaged across targets with equal '
        'weight'),
    (1, 'Pooling across targets invents cross-receptor comparisons and changes the answer'),
    (0, 'Exclusions, stated once and applied throughout'),
    (1, 'One compound (MCL_1 cpd_55) carried an unrelieved clash, VDW +277,293 kcal/mol — '
        'dropped from correlations, leaving 234'),
    (1, 'Censored labels dropped for pairwise analysis, leaving 220'),
], top=1.4, base=15, gap=9)
footer(s, 'scripts/size_matched_pairwise_accuracy.py in the deposited repository')
supp.append(s)

supp.append(fig_slide(
    'Consensus scoring was tested, not assumed — and it does not help',
    FIG / 'consensus_fails.png',
    cite='A speculative claim that the two methods fail on different compounds was written, '
         'then tested and withdrawn. The rank errors correlate at +0.714, so the consensus '
         '(0.445) loses to MM/PBSA alone (0.545).',
    title_size=24))

supp.append(fig_slide(
    'A sealed external test on a sixth target — everything failed',
    FIG / 'external_12ca5.png',
    cite='17 peptidomimetics vs the anti-haemagglutinin antibody 12CA5, affinities sealed '
         'before prediction · the lesson carried into the paper is that end-point scoring '
         'inherits the quality of its input pose, and pose provenance must be reported '
         'alongside any score',
    title_size=24))

s = new_slide()
set_title(s, 'Why the machine-learned force field underperforms: no desolvation term', 23)
bullets(s, [
    (0, 'Four independent lines point the same way'),
    (1, 'It correlates +0.866 with our gas-phase MM term (VDW+EEL), only +0.610 with full '
        'MM/PBSA, and −0.696 with the PB desolvation penalty'),
    (1, 'Its mean magnitude, −131 kcal/mol, sits in the gas-phase range (−109) not the '
        'solvated one (−73) — against an experimental mean of −9.9'),
    (1, 'Its skill pattern against pocket polarity (r = −0.798) nearly matches our gas-phase '
        'term’s (−0.753): it wins only on MCL-1, the greasiest pocket, and collapses on '
        'PCSK9, the most polar'),
    (1, 'The preprint itself states optimisation "is performed in the gas phase while the '
        'input structures were prepared in solvent"'),
    (0, 'Adding PB desolvation to our gas-phase term moves macro ρ from 0.394 to 0.545'),
    (1, 'That +0.151 is the single largest lever in the whole comparison'),
], top=1.4, base=15, gap=9)
footer(s, 'SO3LR-SF: Agha, Suárez-Dou, Kabylda, Tkatchenko & Volkamer, ChemRxiv '
          '10.26434/chemrxiv.15008482/v1 · our MM/GBSA arm reproduces their reported '
          'comparator (0.478 vs 0.44), so the shortfall is specific to the new method')
supp.append(s)

supp.append(fig_slide(
    'The 10 ns stability result did not survive longer sampling',
    FIG / 'md_legacy_vs_current.png',
    cite='The source repository marks the 10 ns files "historical context [that] must not '
         'replace the current conclusion" (ncflow-3-swap/results/00_READ_ME_FIRST.md) · '
         'kept here as the record of why designs are no longer ranked by long '
         'explicit-solvent simulation',
    title_size=24))

s = new_slide()
set_title(s, 'The two FEP+ baselines in this deck are different estimands', 23)
table(s, [
    ['target', 'held-query splits', 'all 235 compounds'],
    ['Cyclin A/B', '0.792', '0.791'],
    ['JAK2', '0.792', '0.788'],
    ['KRAS', '**0.828', '**0.786'],
    ['MCL-1', '**0.903', '**0.859'],
    ['PCSK9', '**0.874', '**0.929'],
    ['**mean', '**0.838', '**0.831'],
], left=2.9, top=1.32, width=7.5, height=2.45, fs=11.5)
bullets(s, [
    (0, '0.838 — FEP+ on the held-out query compounds of each ten-anchor split (~147), '
        'ρ per context, averaged context → family → families equally'),
    (0, '0.831 — FEP+ over all 235 labelled compounds, ρ per family'),
    (0, 'Both reproduce exactly; neither is stale. Three of five targets differ, so this is a '
        'cohort difference, not rounding. The 0.838 KRAS value rests on a two-compound context'),
    (0, 'Consequence: ML and benchmark ρ values must not be placed on one axis'),
], top=3.95, base=13, gap=7)
supp.append(s)

s = new_slide()
set_title(s, 'Where every number on these slides comes from', 24)
table(s, [
    ['claim', 'source'],
    ['ρ +0.545, +0.831, per-target table', 'manuscript Table 3 · results/zha25_cohort_scores.csv'],
    ['pairwise 0.624 [0.555, 0.689]', 'Table 5 · scripts/size_matched_pairwise_accuracy.py'],
    ['KRAS pooled vs macro', 'Table 4'],
    ['size confound −0.43 / −0.55', 'scripts/fig_size_confound.py'],
    ['replica noise 6.81 vs 5.74', 'kras_all40_md_score_and_stability_metrics.csv'],
    ['Rosetta −0.149', 'results/rosetta_vs_experiment_fepplus.csv'],
    ['SO3LR-SF 0.292', 'so3lr_sf_evaluation/so3lr_sf_scores.csv'],
    ['consensus +0.714 / 0.445', 'so3lr_sf_evaluation/consensus.py'],
    ['12CA5 −0.891', 'results/zha25_12ca5_built_pose_scores.csv'],
    ['LptA 212 / canonical / swaps', 'results/zha25_lpta_*.csv'],
    ['registered ten', 'results/prospective_10_denovo_size_matched.csv'],
], left=0.8, top=1.35, width=11.7, height=4.8, fs=11.5, widths=[4.4, 7.3])
footer(s, 'All of the above are in the public deposit: '
          'github.com/Caroldorium/mmpbsa-macrocycle-benchmark · '
          'every figure has a generating script in scripts/')
supp.append(s)

s = new_slide()
set_title(s, 'Numbers in this deck that come from outside the benchmark repository', 23)
table(s, [
    ['claim', 'file', 'dated'],
    ['Morgan+PLEC 0.563, 3NN 0.452, FEP+ 0.838',
     'anchor-probability-and-activity-cliffs-20260916/…/prediction_metrics_equal_family.csv',
     '16 Sep 2026'],
    ['anchor curve 0.626 → 0.735',
     'transfer-prior-and-anchor-learning-20260916/…/learning_curve_fixed_cohort_summary.csv',
     '16 Sep 2026'],
    ['descriptor families 0.493 / 0.509 / 0.524–0.540 / 0.419',
     'residue-energy-fingerprint-fusion-20260917 · mordred-fingerprint-fusion-20260916',
     '16–17 Sep 2026'],
    ['KRAS nine MD metrics, H-bond 0.541',
     'kras-md-metric-correlation-20260917/…/kras_same20_md_metric_correlations.csv',
     '17 Sep 2026'],
    ['250 ns structural + binding summaries',
     'molecular-dynamics/results/01_… and 02_all_13_binding_energy_summary.tsv',
     '13 Aug 2026'],
    ['legacy 10 ns stability (historical)',
     'untitled folder/ncflow-3-swap/results/10_legacy_10ns_all_systems.tsv',
     '31 Jul 2026'],
], left=0.55, top=1.35, width=12.2, height=4.6, fs=10, widths=[3.5, 6.6, 2.1])
footer(s, 'Morgan + PLEC = 0.563 reproduces bit-identically across five independent reruns '
          'between 16 and 23 Sep 2026 · the legacy 10 ns row is retained only as history')
supp.append(s)


# ------------------------------------------------- 5. order and hide
ORDER = [kept['title'], kept['lps'], new['pipeline'], kept['protocol'],
         new['design'], new['pairwise'], new['parity'], new['confound'],
         new['sampling'], new['rosetta'], new['so3lr'], new['ff'],
         kept['mlchart'], kept['anchors'],
         new['lpta212'], new['lptacanon'], new['swap'],
         new['purchase'], kept['reagents'], new['prospective']]
CORE = len(ORDER)
ORDER += [supp[0], supp[1], supp[2], supp[3], supp[4], supp[5], supp[6],
          supp[7], supp[8], supp[9], supp[10],
          kept['descriptors'], kept['fullsweep'], kept['ablation'], supp[11],
          kept['krasmetrics'], supp[12], supp[13]]

by_id = {}
for sid in list(sldIdLst):
    by_id[prs.part.related_part(sid.rId)] = sid
for sid in list(sldIdLst):
    sldIdLst.remove(sid)
for sl in ORDER:
    sldIdLst.append(by_id[sl.part])

for sl in ORDER[:CORE]:
    if sl._element.get('show') == '0':        # slide 2 ships hidden in the source deck
        del sl._element.attrib['show']
for sl in ORDER[CORE:]:
    sl._element.set('show', '0')

prs.save(OUT)
print(f'wrote {OUT}')
print(f'  {len(prs.slides)} slides: {CORE} core + {len(prs.slides) - CORE} supplementary (hidden)')
