#!/usr/bin/env python3
"""Render MANUSCRIPT.md to DOCX and HTML with figures embedded.

The source names each figure as a backticked filename inside its caption rather than as Markdown
image syntax, so pandoc alone emits a text-only document.  This inserts an image directive after
each caption in a build copy; MANUSCRIPT.md itself is left untouched.
"""
import pathlib, re, subprocess, sys

HERE = pathlib.Path(__file__).resolve().parent.parent
MD   = HERE / 'manuscript' / 'MANUSCRIPT.md'
SEARCH = [HERE / 'manuscript', HERE / 'figures']

def find(name):
    for d in SEARCH:
        if (d / name).exists(): return (d / name).resolve()
    return None

src = MD.read_text(encoding='utf-8')
# only inside the figures section: a caption ending in `something.png`
seen, missing = [], []
def sub(m):
    name = m.group(2)
    p = find(name)
    if p is None:
        missing.append(name); return m.group(0)
    seen.append(name)
    return f"{m.group(1)}`{name}`\n\n![]({p.as_posix()})\n"

out = re.sub(r'(\*\*Figure \d+:.*?)`([A-Za-z0-9_\-]+\.png)`', sub, src, flags=re.S)
if missing:
    print('MISSING FIGURES:', ', '.join(missing), file=sys.stderr); sys.exit(1)

build = HERE / 'manuscript' / '.build.md'
build.write_text(out, encoding='utf-8')
rp = ':'.join(str(d) for d in SEARCH)
for fmt, extra in (('docx', []), ('html', ['--standalone', '--embed-resources'])):
    tgt = MD.with_suffix('.' + fmt)
    subprocess.run(['pandoc', str(build), '-o', str(tgt), f'--resource-path={rp}', *extra], check=True)
    print(f'{tgt.name}: {tgt.stat().st_size/1e6:.2f} MB')
build.unlink()
print(f'embedded {len(seen)} figures: {", ".join(seen)}')
