r"""Kimonu.md and Kimonu_Supplement.md -> one LaTeX file (main text, then the Supplementary Materials).

Look of the September preview (kimonu_preview.tex): the original Science template with Times (mathptmx), 12 pt,
1.5 line spacing, 1 inch margins and line numbers, plus microtype for even word spacing. Headings are real
headings: sections, subsections, and the third level as an italic heading on its own line.

The markdown is read as if all markup were accepted: struck text, ==comments== and callouts go, bold proposals
become plain. Citations become numbers from one bibliography (sciencemag.bst). Display equations are numbered and
lose their final punctuation, and a colon that leads into a display equation goes too. Figures sit at first
mention, with the legend as running text below, so a long legend can cross a page.

    python md_to_kimonu_tex.py <Kimonu.md> <Kimonu_Supplement.md> <template dir> <dir with PNGs> <out dir> <kimonu.bib>
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
           "¹": r"\textsuperscript{1}", "⁰": r"\textsuperscript{0}", "⁻": r"\textsuperscript{\ensuremath{-}}", "₁": r"\ensuremath{{}_1}", "∘": r"\ensuremath{\circ}", "⊤": r"\ensuremath{\top}", "≠": r"\ensuremath{\ne}", "∈": r"\ensuremath{\in}",
           "α": r"\ensuremath{\alpha}", "β": r"\ensuremath{\beta}", "γ": r"\ensuremath{\gamma}", "λ": r"\ensuremath{\lambda}", "μ": r"\ensuremath{\mu}", "⇔": r"\ensuremath{\Leftrightarrow}", "·": r"\ensuremath{\cdot}", "⋅": r"\ensuremath{\cdot}", "′": "'", "ᵀ": r"\ensuremath{{}^\top}",
           "≥": r"\ensuremath{\ge}", "∞": r"\ensuremath{\infty}", "±": r"\ensuremath{\pm}", "∼": r"\ensuremath{\sim}", "≡": r"\ensuremath{\equiv}", "→": r"\ensuremath{\rightarrow}"}
tokens = {}


def token(tex):
    k = f"ZZTOK{len(tokens)}ZZ"
    tokens[k] = tex
    return k


def accept(t):
    t = re.sub(r"(?m)^>.*\n?", "", t)
    t = re.sub(r"==\[.*?\]==|==\*.*?\*==|==[^=\n]*==", "", t, flags=re.S)
    t = re.sub(r"~~.+?~~", "", t, flags=re.S)
    t = re.sub(r"(?m)^---\s*$", "", t)                  # separators
    t = re.sub(r"\bFigures (?=\d)", "Figs. ", t)           # Science style in the running text
    t = re.sub(r"(?<!\*\*)\bFigure (?=\d)", "Fig. ", t)  # legend labels (**Figure N.**) stay for caption()
    t = re.sub(r"\b(Figs?\. \d) (?=[A-I]\b)", r"\1", t)    # "Fig. 1 A" -> "Fig. 1A"
    return re.sub(r"(?m)^([^\n-][^\n]*)\n(- )", r"\1\n\n\2", t)   # a list needs a blank line before it


def caption(text):
    """(label, bold title sentence, rest with bold panel letters)."""
    text = text.replace("**", "")
    m = re.match(r"\s*(Figure \d+\.|fig\. S\d+\.)\s*", text, flags=re.I)
    label = m.group(1) if m else ""
    text = text.strip()[len(m.group(0).lstrip()):] if m else text.strip()
    text = re.sub(r"\*\*(\([A-I](?:, [A-I])?\))\*\*", r"\1", text)
    m = re.match(r"(.+?\.)\s+(.*)", text, flags=re.S)
    head, rest = (m.group(1), m.group(2)) if m else (text, "")
    head = head.replace("**", "")
    rest = rest.replace("**", "") if rest.count("**") % 2 == 0 else rest
    rest = re.sub(r"(?<![\w(])\(([A-I](?:\s?[,–]\s?[A-I])*)\)(?=\s)", lambda k: token(r"(\textbf{" + k.group(1) + "})"), rest)
    label = re.sub(r"^Figure (\d+)\.", r"Fig. \1.", label)
    label = re.sub(r"^fig\. (S\d+)\.", r"Fig. \1.", label, flags=re.I)
    return label, head, rest


keys = {d.lower(): k for k, d in re.findall(r"@\w+\{([^,]+),[^@]*?doi = \{([^}]+)\}", open(bib).read())}
CITE = r"\[([^\]]+)\]\(https://doi\.org/((?:[^()\s]|\([^()\s]*\))+)\)"


def key(d):
    d = d.lower().rstrip(".")
    return keys.get(d) or keys[d.replace("(", "%28").replace(")", "%29")]


def cites(t):
    miss = {d for _, d in re.findall(CITE, t) if d.lower().rstrip(".") not in keys and d.lower().replace("(", "%28").replace(")", "%29") not in keys}
    assert not miss, f"not in the bibliography: {miss}"
    t = re.sub(r"\(((?:\s*" + CITE + r"\s*[;,]?)+)\)", lambda m: token(r"\cite{" + ",".join(key(d) for _, d in re.findall(CITE, m.group(1))) + "}"), t)
    t = re.sub(CITE, lambda m: re.sub(r",?\s*\d{4}$", "", m.group(1)) + " " + token(r"\cite{" + key(m.group(2)) + "}"), t)
    return re.sub(r"\[([^\]]+)\]\((https?://[^)]+)\)", lambda m: m.group(1) + " (" + token(r"\url{" + m.group(2) + "}") + ")", t)


eq_no, eq_ref = [0], {}


def equations(t, prefix):
    def one(m):
        body = m.group(1).strip()
        if r"\begin{array}" in body:
            return token("\\[\n\\scriptsize " + body.rstrip(".,;: ") + "\n\\]")     # the matrices of S1, unnumbered
        body = re.sub(r"\s*[.,;:]\s*$", "", body)
        body = re.sub(r"\\tag\{[^}]*\}", "", body)
        eq_no[0] += 1
        lab = f"eq{prefix}{eq_no[0]}"
        for needle, name in EQ_LABELS.items():
            if needle in body:
                eq_ref.setdefault(name, lab)
        return token(r"\begin{equation}\label{" + lab + "}\n" + body + "\n\\end{equation}")
    t = re.sub(r"\$\$(.*?)\$\$", one, t, flags=re.S)
    eqs = "|".join(k for k, v in tokens.items() if v.startswith(("\\begin{equation}", "\\[")))
    return re.sub(r":(\s*)(" + eqs + ")", r"\1\2", t) if eqs else t       # no colon leads into a display equation


def figures(t, label):
    def one(m):
        name = os.path.splitext(os.path.basename(m.group(1).split("|")[0]))[0]
        png = label(name)
        if png is None:
            return ""
        lab, head, rest = caption(m.group(2) or "")
        return "\n\n" + token("@@FIG:" + png + "@@" + lab + " " + head + "@@" + rest + "@@") + "\n\n"
    return re.sub(r"!\[\[([^\]]+)\]\]\s*\n([^\n]*)", one, t)


def pdc(t):
    a = [PANDOC, "-f", "markdown-auto_identifiers-implicit_figures-blank_before_header+tex_math_dollars", "-t", "latex", "--wrap=preserve", "--no-highlight"]
    x = subprocess.run(a, input=t, capture_output=True, text=True, check=True).stdout
    # heading levels: # section, ## subsection, ### italic display heading (\paragraph redefined below)
    x = x.replace("\\subsubsection{", "\\paragraph*{").replace("\\subsection{", "\\subsection*{").replace("\\section{", "\\section*{")
    return x


def finish(x):
    for _ in range(5):
        def fig(m):
            png, head, rest = m.group(1), m.group(2), m.group(3)
            return ("\n\\begin{figure}[H]\n\\centering\n\\includegraphics[width=\\textwidth,height=0.85\\textheight,keepaspectratio]{" + png + "}\n\\end{figure}\n"
                    "\\noindent{\\small\\textbf{" + pdc(head).strip() + "} " + pdc(rest).strip() + "\\par}\n\\medskip\n\n")
        x = re.sub(r"@@FIG:([^@]+)@@([^@]*)@@([^@]*)@@", fig, x)
        for k, v in tokens.items():
            x = x.replace(k, v)
    x = re.sub(r"\\texorpdfstring\{([^{}]*)\}\{[^{}]*\}", r"\1", x)
    def spec(m):      # natural column widths instead of pandoc's equal p-columns
        cols = re.findall(r"\\(raggedright|raggedleft|centering)\\arraybackslash\}p\{", m.group(1))
        return "\\begin{longtable}[]{@{}" + "".join({"raggedright": "l", "raggedleft": "r", "centering": "c"}[c] for c in cols) + "@{}}" if cols else m.group(0)
    x = re.sub(r"\\begin\{longtable\}\[\]\{@\{\}(.*?)@\{\}\}", spec, x, flags=re.S)
    x = re.sub(r"\\begin\{minipage\}\[b\]\{\\linewidth\}\\raggedright\s*(.*?)\s*\\end\{minipage\}", r"\1", x, flags=re.S)
    x = re.sub(r"\\begin\{minipage\}\[b\]\{\\linewidth\}\\raggedleft\s*(.*?)\s*\\end\{minipage\}", r"\1", x, flags=re.S)
    return x


raw = open(md_path, encoding="utf-8").read()
field = lambda name: re.search(r"\*\*" + name + r":\*\*\s*(.+)", raw).group(1).strip()
title, short, affil = field("Title"), field("Short title"), field("Affiliation")
affil_tex = re.sub(r",? *MA \d{5}", ", MA", affil).replace("Research Institute, ", "Research Institute,\\\\ \\small ")
orcid = re.search(r"ORCID (\S+)", field("Author")).group(1)
md = accept(raw)
md = cites(md)
md = equations(md, "")
md = figures(md, lambda name: name)
sec = {}
for p in re.split(r"(?m)^# ", md)[1:]:
    h, _, b = p.partition("\n")
    sec.setdefault(h.strip(), b)


def refs(t):
    t = re.sub(r"\(Equation \[ref: lin/log eq\]\)", lambda m: token(r"(Eq.~\ref{" + eq_ref["lin/log eq"] + "})"), t)
    t = re.sub(r"Equation \[ref: lin/log eq\]", lambda m: token(r"Eq.~\ref{" + eq_ref["lin/log eq"] + "}"), t)
    t = re.sub(r"Equation \[ref: full problem\]", lambda m: token(r"Eq.~\ref{" + eq_ref["full problem"] + "}"), t)
    t = re.sub(r"\[eq:psi_c/psi_d\]", lambda m: token(r"Eq.~\ref{" + eq_ref["eq:psi_c/psi_d"] + "}"), t)
    assert "[ref:" not in t and "[eq:" not in t, re.findall(r"\[(?:ref|eq):[^\]]*\]", t)
    return t


def h(name, text):
    return f"\\section*{{{name}}}\n" + pdc(refs(text).strip()) + "\n"


meth = re.sub(r"(?m)^### ", "## ", sec["Materials and Methods"])
body = h("Introduction", sec["Introduction"]) + h("Results", sec["Results"]) + h("Discussion", sec["Discussion"]) + h("Materials and Methods", meth)
smlist = pdc(sec["Supplementary Materials"].strip())
abstract = pdc(sec["Abstract"].strip())
ack = sec["Acknowledgments"].replace("archived on Zenodo .", "archived on Zenodo (DOI to be added).")
ack = re.sub(r"(?m)^(Funding|Author contributions|Use of AI tools|Competing interests|Data, code and materials availability):\s*", r"\n\n**\1:** ", ack)
ack_tex = pdc(ack.strip())

# ---- Supplement
sm = accept(open(sm_path, encoding="utf-8").read())
sm = sm[sm.index("**This PDF includes**"):]
sm = re.sub(r"\(figures/supp/(DataS\d_[\w]+\.csv)\)", r"(\1)", sm)
sm = cites(sm)
eq_no[0] = 0
sm = equations(sm, "S")
sm = figures(sm, lambda name: name)
sm = re.sub(r"(?m)^## ", "# ", re.sub(r"(?m)^### ", "## ", sm))
sm_tex = pdc(sm.strip())

pre = r"""\documentclass[12pt]{article}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{mathptmx}
\usepackage{amsmath,amssymb}
\usepackage[mathlines]{lineno}
\renewcommand\linenumberfont{\normalfont\small}
\usepackage{microtype}
\usepackage{graphicx}
\usepackage{float}
\usepackage{longtable,booktabs,array,calc}
\usepackage[letterpaper,margin=1in]{geometry}
\linespread{1.5}
\frenchspacing
\renewenvironment{abstract}{\quotation}{\endquotation}
\date{}
\renewcommand\refname{References and Notes}
\usepackage{scicite}
\usepackage{url}
\providecommand{\tightlist}{\setlength{\itemsep}{0pt}\setlength{\parskip}{0pt}}
\providecommand{\pandocbounded}[1]{#1}
\providecommand{\texorpdfstring}[2]{#1}
\makeatletter
% line numbers next to display equations (amsmath environments)
\newcommand*\patchAmsMathEnvironmentForLineno[1]{%
  \expandafter\let\csname old#1\expandafter\endcsname\csname #1\endcsname
  \expandafter\let\csname oldend#1\expandafter\endcsname\csname end#1\endcsname
  \renewenvironment{#1}{\linenomath\csname old#1\endcsname}{\csname oldend#1\endcsname\endlinenomath}}
\newcommand*\patchBothAmsMathEnvironmentsForLineno[1]{\patchAmsMathEnvironmentForLineno{#1}\patchAmsMathEnvironmentForLineno{#1*}}
\AtBeginDocument{\patchBothAmsMathEnvironmentsForLineno{equation}}
% headings: section bold large, subsection bold, third level italic on its own line
\renewcommand\section{\@startsection{section}{1}{\z@}{-3ex \@plus -1ex \@minus -.2ex}{1.2ex \@plus .2ex}{\normalfont\large\bfseries}}
\renewcommand\subsection{\@startsection{subsection}{2}{\z@}{-2.5ex \@plus -1ex \@minus -.2ex}{0.8ex \@plus .2ex}{\normalfont\normalsize\bfseries}}
\renewcommand\paragraph{\@startsection{paragraph}{4}{\z@}{-2ex \@plus -.5ex \@minus -.2ex}{0.5ex \@plus .1ex}{\normalfont\normalsize\itshape}}
\renewcommand\verbatim@font{\normalfont\fontsize{7}{8.5}\selectfont\ttfamily}
\makeatother
\setlength{\LTleft}{0pt}
\usepackage{etoolbox}
\AtBeginEnvironment{longtable}{\scriptsize\linespread{1.1}\selectfont\setlength{\tabcolsep}{3pt}}
"""
pre += "".join(f"\\DeclareUnicodeCharacter{{{ord(c):04X}}}{{{v}}}\n" for c, v in UNICODE.items())
pre += "\\allowdisplaybreaks\n"
doc = pre + r"""
\begin{document}
\linenumbers
\title{\bfseries\boldmath """ + pdc(title).strip() + r"""}
\author{Philipp~C.~Schneider$^{1\ast}$\\[1ex]
\small$^{1}$""" + affil_tex + r"""\\
\small$^\ast$Corresponding author. Email: p.schneider.chuzel@gmail.com\\
\small ORCID """ + orcid + r"""}
\maketitle
\begin{abstract}
""" + abstract + r"""
\end{abstract}
\noindent\textbf{Short title:} """ + short + r"""
\clearpage

""" + body + r"""
\section*{Supplementary Materials}
""" + smlist + r"""
\bibliography{kimonu}
\bibliographystyle{sciencemag}
\section*{Acknowledgments}
""" + ack_tex + r"""
\clearpage
\setcounter{figure}{0}\setcounter{table}{0}\setcounter{equation}{0}\setcounter{page}{1}
\renewcommand{\thefigure}{S\arabic{figure}}
\renewcommand{\thetable}{S\arabic{table}}
\renewcommand{\theequation}{S\arabic{equation}}
\renewcommand{\thepage}{S\arabic{page}}
\begin{center}
{\Large\bfseries Supplementary Materials for}\\[1ex]
{\large\bfseries\boldmath """ + pdc(title).strip() + r"""}\\[1ex]
Philipp~C.~Schneider\\
\small Corresponding author. Email: p.schneider.chuzel@gmail.com
\end{center}
""" + sm_tex + "\n\\end{document}\n"
doc = finish(doc)
assert "@@FIG" not in doc and "ZZTOK" not in doc
open(os.path.join(out, "kimonu.tex"), "w", encoding="utf-8").write(doc)
for f in os.listdir(figs):
    if f.endswith(".png"):
        shutil.copy(os.path.join(figs, f), out)
log = ""
for cmd in ("pdflatex", "bibtex", "pdflatex", "pdflatex"):
    r = subprocess.run([cmd] + (["-interaction=nonstopmode"] if cmd == "pdflatex" else []) + ["kimonu" + (".tex" if cmd == "pdflatex" else "")],
                       cwd=out, capture_output=True, errors="replace")
    log = r.stdout
print("\n".join(l for l in log.split("\n") if l.startswith("!") or "Output written" in l or "undefined" in l.lower() or "Overfull" in l)[:4000])
