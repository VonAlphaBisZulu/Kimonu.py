r"""Kimonu.md and Kimonu_Supplement.md -> one LaTeX file in the Science template, in Science Advances form.

Reads the markdown as if all markup were accepted: struck text and ==comments== go, bold proposals become plain.
Citations become numbers from one bibliography for main text and Supplement, references cited only in the
Supplement last. Display equations are numbered, with [ref: ...] placeholders resolved by EQ_LABELS. Figures sit
at first mention with their legends. Availability moves from Materials and Methods to Acknowledgments.

newtx is not installed on this cluster, so the build uses the template's mathptmx fallback, plus amsmath, amssymb
and longtable, which the template asks not to add. Science Advances has no required format at first submission.

    python md_to_sciadv.py <Kimonu.md> <Kimonu_Supplement.md> <template dir> <dir with figure PNGs> <out dir> <kimonu.bib>
PNG names: Figure1.png ... and FigS_<name>.png, rendered from the vault SVGs.
"""
import os, re, shutil, subprocess, sys

md_path, sm_path, tpl, figs, out, bib = sys.argv[1:7]
PANDOC = os.environ.get("PANDOC", os.path.expanduser("~/miniconda3/envs/bioinfo/bin/pandoc"))
os.makedirs(out, exist_ok=True)
for f in ("scicite.sty", "sciencemag.bst"):
    shutil.copy(os.path.join(tpl, f), out)
shutil.copy(bib, os.path.join(out, "kimonu.bib"))

EQ_LABELS = {r"\text{mass balance: linear in }": "lin/log eq", r"\mathbf{v}^\top\ln{\mathbf{x}}": "eq:psi_c/psi_d", r"E(\mathbf{k}) = \{": "full problem"}
UNICODE = {"∅": r"\ensuremath{\varnothing}", "→": r"\ensuremath{\rightarrow}", "⇌": r"\ensuremath{\rightleftharpoons}", "≤": r"\ensuremath{\le}", "≥": r"\ensuremath{\ge}",
           "−": r"\ensuremath{-}", "×": r"\ensuremath{\times}", "⊆": r"\ensuremath{\subseteq}", "⊇": r"\ensuremath{\supseteq}", "≈": r"\ensuremath{\approx}", "∝": r"\ensuremath{\propto}",
           "√": r"\ensuremath{\surd}", "ψ": r"\ensuremath{\psi}", "ε": r"\ensuremath{\varepsilon}", "δ": r"\ensuremath{\delta}", "τ": r"\ensuremath{\tau}", "σ": r"\ensuremath{\sigma}", "Δ": r"\ensuremath{\Delta}",
           "²": r"\textsuperscript{2}", "³": r"\textsuperscript{3}", "⁴": r"\textsuperscript{4}", "⁵": r"\textsuperscript{5}", "⁶": r"\textsuperscript{6}", "⁷": r"\textsuperscript{7}", "⁸": r"\textsuperscript{8}", "⁹": r"\textsuperscript{9}",
           "¹": r"\textsuperscript{1}", "⁰": r"\textsuperscript{0}", "₁": r"\ensuremath{{}_1}", "∘": r"\ensuremath{\circ}", "⊤": r"\ensuremath{\top}", "≠": r"\ensuremath{\ne}", "∈": r"\ensuremath{\in}", "α": r"\ensuremath{\alpha}", "β": r"\ensuremath{\beta}", "γ": r"\ensuremath{\gamma}", "λ": r"\ensuremath{\lambda}", "⇔": r"\ensuremath{\Leftrightarrow}", "·": r"\ensuremath{\cdot}", "⋅": r"\ensuremath{\cdot}", "′": "'", "ᵀ": r"\ensuremath{{}^\top}"}
tokens = {}


def token(tex):
    k = f"ZZTOK{len(tokens)}ZZ"
    tokens[k] = tex
    return k


def accept(t):
    t = re.sub(r"(?m)^>.*\n?", "", t)
    t = re.sub(r"==\[.*?\]==|==\*.*?\*==|==[^=\n]*==", "", t, flags=re.S)
    t = re.sub(r"~~.+?~~", "", t, flags=re.S)
    return re.sub(r"\*\*(.+?)\*\*", lambda m: m.group(0) if m.group(1) in ("KI", "MO", "NU") else m.group(1), t)


def caption(text):
    """(label as written, bold title sentence, rest with bold panel letters). The legend is running text, so it can cross a page."""
    m = re.match(r"\s*(Figure \d+\.|fig\. S\d+\.)\s*", text, flags=re.I)
    label = m.group(1) if m else ""
    text = text.strip()[len(m.group(0).lstrip()):] if m else text.strip()
    m = re.match(r"(.+?\.)\s+(.*)", text, flags=re.S)
    head, rest = (m.group(1), m.group(2)) if m else (text, "")
    rest = re.sub(r"(?<![\w(])\(([A-I](?:\s?[,–]\s?[A-I])*)\)(?=\s)", lambda k: token(r"(\textbf{" + k.group(1) + "})"), rest)
    return label, head, rest


keys = {d.lower(): k for k, d in re.findall(r"@\w+\{([^,]+),[^@]*?doi = \{([^}]+)\}", open(bib).read())}
CITE = r"\[([^\]]+)\]\(https://doi\.org/((?:[^()\s]|\([^()\s]*\))+)\)"


def cites(t):
    miss = {d for _, d in re.findall(CITE, t) if d.lower().rstrip(".") not in keys}
    assert not miss, f"not in the bibliography: {miss}"
    t = re.sub(r"\(((?:\s*" + CITE + r"\s*[;,]?)+)\)", lambda m: token(r"\cite{" + ",".join(keys[d.lower()] for _, d in re.findall(CITE, m.group(1))) + "}"), t)
    t = re.sub(CITE, lambda m: re.sub(r",?\s*\d{4}$", "", m.group(1)) + " " + token(r"\cite{" + keys[m.group(2).lower()] + "}"), t)
    return re.sub(r"\[([^\]]+)\]\((https?://[^)]+)\)", r"\1", t)


eq_no, eq_ref = [0], {}


def equations(t, prefix):
    def one(m):
        body = m.group(1).strip()
        eq_no[0] += 1
        lab = f"eq{prefix}{eq_no[0]}"
        for needle, name in EQ_LABELS.items():
            if needle in body:
                eq_ref[name] = lab
        return token(r"\begin{equation}\label{" + lab + "}\n" + body + "\n\\end{equation}")
    return re.sub(r"\$\$(.*?)\$\$", one, t, flags=re.S)


def figures(t, label):
    """An embed followed by its legend line becomes a figure at that place. Returns text; label(name) gives the PNG, or None to drop."""
    def one(m):
        name = os.path.splitext(os.path.basename(m.group(1).split("|")[0]))[0]
        png = label(name)
        if png is None:
            return ""
        label, head, rest = caption(m.group(2) or "")
        return token("@@FIG:" + png + "@@" + label + " " + head + "@@" + rest + "@@")
    return re.sub(r"!\[\[([^\]]+)\]\]\s*\n([^\n]*)", one, t)


def pdc(t, shift=0):
    a = [PANDOC, "-f", "markdown-auto_identifiers-implicit_figures-blank_before_header+tex_math_dollars", "-t", "latex", "--wrap=preserve", "--no-highlight"]
    if shift:
        a.append(f"--shift-heading-level-by={shift}")
    return subprocess.run(a, input=t, capture_output=True, text=True, check=True).stdout


def finish(x):
    for _ in range(4):
        def fig(m):
            png, head, rest = m.group(1), m.group(2), m.group(3)
            return ("\n\\begin{figure}[H]\n\\centering\n\\includegraphics[width=\\textwidth,height=0.8\\textheight,keepaspectratio]{" + png + "}\n\\end{figure}\n"
                    "\\noindent{\\small\\textbf{" + pdc(head).strip() + "} " + pdc(rest).strip() + "}\n\n")
        x = re.sub(r"@@FIG:([^@]+)@@([^@]*)@@([^@]*)@@", fig, x)
        for k, v in tokens.items():
            x = x.replace(k, v)
    x = re.sub(r"\\texorpdfstring\{([^{}]*)\}\{[^{}]*\}", r"\1", x)
    return x


raw = open(md_path, encoding="utf-8").read()
title = re.search(r"\*\*Title:\*\*\s*(.+)", raw).group(1).strip()
short = re.search(r"\*\*Short title:\*\*\s*(.+)", raw).group(1).strip()
md = accept(raw)
md = re.sub(r"!\[\[figures/Figure1_draft2\.svg[^\]]*\]\]\n?", "", md)
md = cites(md)
md = equations(md, "")
n = {"Figure1_draft": "Figure1", "Figure2_draft": "Figure2", "Figure3_draft": "Figure3"}
md = figures(md, lambda name: n.get(name))
sec = {}
for p in re.split(r"(?m)^# ", md)[1:]:
    h, _, b = p.partition("\n")
    sec.setdefault(h.strip(), b)
for name, lab in eq_ref.items():
    pass


def refs(t):
    t = t.replace(" (Equation [ref: N=YA factorization])", "")
    t = re.sub(r"Equation \[ref: lin/log eq\]", lambda m: token(r"Eq.~\ref{" + eq_ref["lin/log eq"] + "}"), t)
    t = re.sub(r"Equation \[ref: full problem\]", lambda m: token(r"Eq.~\ref{" + eq_ref["full problem"] + "}"), t)
    t = re.sub(r"\[eq:psi_c/psi_d\]", lambda m: token(r"Eq.~\ref{" + eq_ref["eq:psi_c/psi_d"] + "}"), t)
    return t


meth = sec["Materials and Methods"]
avail = re.search(r"\n(Kimonu is available on GitHub.*?)\n", meth + "\n").group(1)
meth = meth.replace(avail, "").replace("### Availability\n", "")
body = "\\section*{Introduction}\n" + pdc(refs(sec["Introduction"]).strip())
body += "\\section*{Results}\n" + pdc(refs(sec["Results"]).strip()).replace("\\subsubsection{", "\\subsubsection*{").replace("\\subsection{", "\\subsection*{")
body += "\\section*{Discussion}\n" + pdc(refs(sec["Discussion"]).strip())
body += "\\section*{Materials and Methods}\n" + pdc(refs(meth).strip()).replace("\\subsubsection{", "\\subsubsection*{").replace("\\subsection{", "\\subsection*{")
abstract = pdc(sec["Abstract"].strip())

# ---- Supplement
sm = open(sm_path, encoding="utf-8").read()
sm = accept(sm)
sm = sm[sm.index("## S1. The example network"):]
sm = cites(sm)
eq_no[0] = 0
sm = equations(sm, "S")
sm = figures(sm, lambda name: name)
sm_tex = pdc(sm.strip()).replace("\\subsubsection{", "\\subsubsection*{").replace("\\subsection{", "\\subsection*{")

tp = open(os.path.join(tpl, "science_template.tex"), encoding="utf-8").read()
pre = tp[:tp.index("\\begin{document}")]
pre = pre.replace("\\usepackage{newtxtext,newtxmath}", "\\usepackage[utf8]{inputenc}\n\\usepackage{mathptmx}\n\\usepackage{amsmath,amssymb}\n\\usepackage{longtable,booktabs,array,calc,float}\n"
                  "\\providecommand{\\tightlist}{\\setlength{\\itemsep}{0pt}\\setlength{\\parskip}{0pt}}\n\\providecommand{\\pandocbounded}[1]{#1}\n\\providecommand{\\texorpdfstring}[2]{#1}")
pre += "".join(f"\\DeclareUnicodeCharacter{{{ord(c):04X}}}{{{v}}}\n" for c, v in UNICODE.items())
pre += "\\allowdisplaybreaks\n\\emergencystretch=2em\n"
doc = pre + r"""
\begin{document}
\renewcommand{\scititle}{""" + title + r"""}
\title{\bfseries\boldmath\scititle}
\author{[Authors and affiliations to be added]\\ \small Short title: """ + short + r"""}
\date{}
\maketitle
\begin{abstract} \bfseries \boldmath
""" + abstract + r"""
\end{abstract}
\noindent
""" + body + r"""
\clearpage
\bibliography{kimonu}
\bibliographystyle{sciencemag}
\section*{Acknowledgments}
\paragraph*{Funding:} [to be added]
\paragraph*{Author contributions:} [to be added, by initials]
\paragraph*{Competing interests:} The authors declare that they have no competing interests. [to be confirmed]
\paragraph*{Data, code and materials availability:}
""" + pdc(cites(accept(avail))) + r"""
\clearpage
\setcounter{figure}{1}\setcounter{table}{0}\setcounter{equation}{0}\setcounter{page}{1}
\renewcommand{\thefigure}{S\arabic{figure}}
\renewcommand{\thetable}{S\arabic{table}}
\renewcommand{\theequation}{S\arabic{equation}}
\renewcommand{\thepage}{S\arabic{page}}
\section*{Supplementary Materials for\\ \scititle}
""" + sm_tex + "\n\\end{document}\n"
doc = finish(doc)
open(os.path.join(out, "kimonu_sciadv.tex"), "w", encoding="utf-8").write(doc)
for f in os.listdir(figs):
    if f.endswith(".png"):
        shutil.copy(os.path.join(figs, f), out)
log = ""
for cmd in ("pdflatex", "bibtex", "pdflatex", "pdflatex"):
    r = subprocess.run([cmd] + (["-interaction=nonstopmode"] if cmd == "pdflatex" else []) + ["kimonu_sciadv" + (".tex" if cmd == "pdflatex" else "")],
                       cwd=out, capture_output=True, errors="replace")
    log = r.stdout
print("\n".join(l for l in log.split("\n") if l.startswith("!") or "Output written" in l or "undefined" in l.lower())[:3000])
