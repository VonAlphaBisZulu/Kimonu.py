"""Collect every doi.org link in the outline into a BibTeX file for the Science template.

Entries come from doi.org content negotiation; keys are FirstAuthorYear (a/b suffix on clashes).
Titles are double-braced so the style keeps their capitalisation (METATOOL, COBRA-k).
Existing entries are kept, so hand edits to the .bib survive a rerun.

    python make_bib.py <Kimonu.md> <out.bib>
"""
import re
import sys
import unicodedata
import urllib.request

TEX = {'ü': '{\\"u}', 'ö': '{\\"o}', 'ä': '{\\"a}', 'á': "{\\'a}", 'é': "{\\'e}", 'ñ': '{\\~n}'}

md_path, bib_path = sys.argv[1:3]
dois = list(dict.fromkeys(d.rstrip('.').lower() for d in
                          re.findall(r'https://doi\.org/([^)\s]+)', open(md_path).read())))
try:
    old = open(bib_path).read()
except FileNotFoundError:
    old = ''
have = {d.lower() for d in re.findall(r'(?im)^\s*doi\s*=\s*\{([^}]+)\}', old)}
keys = set(re.findall(r'@\w+\{([^,]+),', old))

out = [old.rstrip()] if old.strip() else []
for doi in dois:
    if doi in have:
        continue
    req = urllib.request.Request('https://doi.org/' + doi, headers={'Accept': 'application/x-bibtex'})
    raw = urllib.request.urlopen(req, timeout=30).read().decode().strip()
    typ = re.match(r'@(\w+)\{', raw).group(1)
    fields = dict((k.lower(), v) for k, v in re.findall(r'(\w+)=\{((?:[^{}]|\{[^{}]*\})*)\}', raw))
    fields.update((k, v) for k, v in re.findall(r'(\w+)=(\d+)\b', raw) if k not in fields)
    surname = fields['author'].split(',')[0].split(' and ')[0].strip()
    base = unicodedata.normalize('NFKD', surname).encode('ascii', 'ignore').decode()
    base = re.sub(r'\W', '', base) + fields.get('year', '')
    key, n = base, 0
    while key in keys:
        n += 1
        key = base + 'abcdefg'[n - 1]
    keys.add(key)
    fields['doi'] = doi
    fields.pop('url', None)
    fields.pop('month', None)
    fields['title'] = '{' + fields['title'].rstrip('.').replace(' – ', ' -- ') + '}'
    if 'journal' not in fields:  # preprint: the template's form is a note naming server and id
        fields.pop('publisher', None)
        fields['note'] = 'bioRxiv ' + doi.split('/', 1)[1]
    for k in ('issn', 'publisher', 'editor'):
        fields.pop(k, None)
    author = fields['author'].replace('√°', 'á').replace('√±', 'ñ')  # Mac-Roman mojibake in some Crossref records
    fields['author'] = ''.join(TEX.get(c, c) for c in author)  # bibtex drops raw non-ASCII
    if 'pages' in fields:
        fields['pages'] = fields['pages'].replace('–', '--')
    body = ',\n'.join(f'  {k.lower()} = {{{v}}}' for k, v in fields.items())
    out.append(f'@{typ.capitalize()}{{{key},\n{body}\n}}')
    print(key, doi)
open(bib_path, 'w').write('\n\n'.join(out) + '\n')
