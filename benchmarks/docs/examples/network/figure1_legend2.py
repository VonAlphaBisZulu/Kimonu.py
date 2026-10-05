"""Figure 1: the two kinetic modules get their own legend entries; module 2 (C+E, A+E) in mint."""
import sys
src, out = sys.argv[1:3]
s = open(src).read()
MINT = "#cdf3df"
def rep(a, b, n=1):
    global s
    assert s.count(a) == n, (s.count(a), a[:60]); s = s.replace(a, b)
# module 2 outline in panel B
rep('<rect x="465.6" y="112.6" width="47" height="58.4" rx="8" ry="8" fill="#c9ecf4"',
    f'<rect x="465.6" y="112.6" width="47" height="58.4" rx="8" ry="8" fill="{MINT}"')
# its Y^T.T chips in panel F
for y in ("373.176", "403.176"):
    rep(f'<rect x="364.346" y="{y}" width="19.624" height="8.6" rx="2" ry="2" fill="#c9ecf4"/>',
        f'<rect x="364.346" y="{y}" width="19.624" height="8.6" rx="2" ry="2" fill="{MINT}"/>')
# legend: numbered entries, module 2 in the free third row of the right column
rep('<tspan x="0" y="0">kinetic module</tspan>', '<tspan x="0" y="0">kinetic module 1</tspan>')
sw = '<rect x="416.462" y="195.838" width="7.853" height="6.425" rx="3.165" ry="3.165" fill="#c9ecf4" stroke="#000" stroke-dasharray="6 2" stroke-miterlimit="10" stroke-width="1.5"/>'
sw2 = (sw.replace('y="195.838"', 'y="209.833"').replace('fill="#c9ecf4"', f'fill="{MINT}"')
       + '<text transform="translate(430.329 215.16)" fill="#1b2430" font-family="ArialMT, Arial" font-size="8"><tspan x="0" y="0">kinetic module 2</tspan></text>')
rep(sw, sw + "\n  " + sw2)
open(out, "w").write(s)
