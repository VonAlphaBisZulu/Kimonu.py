"""Preview build: Obsidian outline (Kimonu.md) -> Science template (.tex/.pdf).

Strips callouts, ==editorial comments== and image embeds; links become plain text (numeric
citations come later); sections from Introduction to Materials and Methods become the main text;
figures are placed at the end at full text width. The preview swaps newtx (not installed here)
for the template's mathptmx fallback and adds amsmath/amssymb, which the template itself says
not to add -- fine for a look, to be revisited for submission.

    python md_to_science_tex.py <Kimonu.md> <template_dir> <figure_dir> <out_dir>
"""
import os
import re
import shutil
import subprocess
import sys

md_path, tpl_dir, fig_dir, out = sys.argv[1:5]
PANDOC = os.environ.get("PANDOC", "/home/sadreyev_lab/pschneider/miniconda3/envs/bioinfo/bin/pandoc")
os.makedirs(out, exist_ok=True)
for f in ("scicite.sty", "sciencemag.bst"):
    shutil.copy(os.path.join(tpl_dir, f), out)
for f in ("Figure1_draft.png", "Figure2_draft.png"):
    shutil.copy(os.path.join(fig_dir, f), out)

raw = open(md_path).read()
md = re.sub(r'(?m)^>.*\n?', '', raw)
md = re.sub(r'==\[.*?\]==|==\*.*?\*==', '', md, flags=re.S)
md = re.sub(r'!\[\[[^\]]*\]\]', '', md)
md = re.sub(r'\[([^\]]+)\]\((https?://[^)]+)\)', r'\1', md)
parts = re.split(r'(?m)^# ', md)
sec = {p.split('\n', 1)[0].strip(): (p.split('\n', 1)[1] if '\n' in p else '') for p in parts[1:]}


def pdc(text, shift=True):
    args = [PANDOC, '-f', 'markdown-auto_identifiers', '-t', 'latex', '--wrap=preserve']
    if shift:
        args.append('--shift-heading-level-by=1')
    return subprocess.run(args, input=text, capture_output=True, text=True, check=True).stdout


abstract = pdc(sec['Abstract'].strip(), False)
body = pdc(sec['Introduction'].strip(), False)
for k in ['Theoretical background', 'Kimonu', 'Results', 'Discussion', 'Materials and Methods']:
    if k in sec:
        body += '\n\\subsection*{' + k + '}\n' + pdc(sec[k].strip())
# pandoc emits numbered \subsection/\subsubsection; the template uses unnumbered headings
body = body.replace('\\subsubsection{', '\\paragraph*{').replace('\\subsection{', '\\subsubsection*{')

tpl = open(os.path.join(tpl_dir, 'science_template.tex')).read()
pre = tpl[:tpl.index('\\begin{document}')]
pre = pre.replace('\\usepackage{newtxtext,newtxmath}',
                  '\\usepackage{mathptmx}\n\\usepackage{amsmath,amssymb}\n'
                  '\\providecommand{\\tightlist}{\\setlength{\\itemsep}{0pt}\\setlength{\\parskip}{0pt}}')
title = re.search(r'\*\*Working title:\*\*\s*(.+)', raw).group(1).strip()

doc = pre + r'''
\begin{document}
\title{\bfseries\boldmath ''' + title + r'''}
\author{Philipp~Schneider$^{1\ast}$\\ \small$^{1}$Affiliation to be added.\\ \small$^\ast$Corresponding author.}
\date{}
\maketitle
\begin{abstract} \bfseries \boldmath
''' + abstract + r'''
\end{abstract}
\noindent
''' + body + r'''
\clearpage
\begin{figure}
\centering
\includegraphics[width=\textwidth]{Figure1_draft}
\caption{\textbf{Null-space readouts on the example network of Langary et al.}
(\textbf{A}) Metabolic network. (\textbf{B}) Complex network with concordance modules (colours) and the kinetic module (bold outline).
(\textbf{C}) Steady-state basis $\mathbf{K}=\ker(\mathbf{N})$ with fully coupled reactions and elementary modes.
(\textbf{D}) Complex activities $\mathbf{M}=\mathbf{A}\mathbf{K}$: zero rows balanced, proportional rows concordant.}
\end{figure}
\begin{figure}
\centering
\includegraphics[width=\textwidth,height=0.8\textheight,keepaspectratio]{Figure2_draft}
\caption{\textbf{Reproduction of key figures across networks.}
(\textbf{A}) Kinetic modules and (\textbf{B}) ACRR pairs on ten networks of increasing size. (\textbf{C}) Runtime against network size for all 33 networks, log-log fits with 95\% confidence bands.}
\end{figure}
\end{document}
'''
open(os.path.join(out, 'kimonu_preview.tex'), 'w').write(doc)
for _ in range(2):
    subprocess.run(['pdflatex', '-interaction=nonstopmode', 'kimonu_preview.tex'], cwd=out, capture_output=True)
print('built', os.path.join(out, 'kimonu_preview.pdf'))
