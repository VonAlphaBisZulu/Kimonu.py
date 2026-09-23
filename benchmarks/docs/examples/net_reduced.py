import math
S={'A':(120,330),'B':(470,100),'C':(870,330),'D':(740,305),'E':(470,640),'F':(575,245),'G':(470,380),'H':(700,540),'I':(250,540)}
# reaction, substrates, products, reversible, stem (J1 -> J2) or None, label (x, y, anchor), bend / fork spread
RX=[('R3',['B'],['C'],0,None,(700,205,'start'),30),('R5',['C'],['D'],0,None,(793,298,'middle'),0),
    ('R6',['D','E'],['F'],0,((655,335),(625,295)),(612,326,'end'),70),
    ('R8',['F'],['B','G'],0,((548,245),(515,245)),(530,232,'middle'),45),
    ('R9',['C','G'],['H'],0,((665,445),(685,490)),(697,462,'start'),70),
    ('R11',['H'],['C','E'],0,((728,572),(748,596)),(762,588,'start'),28),
    ('R12',['A','G'],['I'],0,((280,445),(265,490)),(262,458,'end'),70),
    ('R14',['I'],['A','E'],0,((222,572),(202,596)),(188,588,'end'),28)]
rs=24; T=rs+3
def unit(dx,dy): L=math.hypot(dx,dy) or 1; return dx/L,dy/L
def bez(p0,p1,p2,p3,t):
    u=1-t; return tuple(u**3*a+3*u*u*t*b+3*u*t*t*c+t**3*d for a,b,c,d in zip(p0,p1,p2,p3))
def split(p0,p1,p2,p3,t):
    lerp=lambda a,b:(a[0]+(b[0]-a[0])*t,a[1]+(b[1]-a[1])*t)
    a,b,c=lerp(p0,p1),lerp(p1,p2),lerp(p2,p3); d,e=lerp(a,b),lerp(b,c); f=lerp(d,e)
    return (p0,a,d,f),(f,e,c,p3)
def trim(curve, center, at_end):
    # find t where the curve leaves the node circle
    ts=[i/400 for i in range(401)]
    if at_end:
        t=next(t for t in reversed(ts) if math.dist(bez(*curve,t),center)>=T); return split(*curve,t)[0]
    t=next(t for t in ts if math.dist(bez(*curve,t),center)>=T); return split(*curve,t)[1]
def d_(c): return 'M{:.1f},{:.1f} C{:.1f},{:.1f} {:.1f},{:.1f} {:.1f},{:.1f}'.format(*c[0],*c[1],*c[2],*c[3])
out=[]; lab=[]
for r,subs,prods,rev,stem,L,bend in RX:
    ms=' marker-start="url(#ah)"' if rev else ''
    if stem is None:
        (x1,y1),(x2,y2)=S[subs[0]],S[prods[0]]; nx,ny=unit(-(y2-y1),x2-x1); m=((x1+x2)/2+nx*bend,(y1+y2)/2+ny*bend)
        c=((x1,y1),((x1+2*m[0])/3,(y1+2*m[1])/3),((x2+2*m[0])/3,(y2+2*m[1])/3),(x2,y2))
        c=trim(trim(c,S[subs[0]],False),S[prods[0]],True); out.append(f'<path d="{d_(c)}"{ms} marker-end="url(#ah)"/>')
    else:
        J1,J2=stem; dx,dy=unit(J2[0]-J1[0],J2[1]-J1[1])
        for s in subs:
            P=S[s]; k=min(0.45*math.dist(P,J1),70); ux,uy=unit(J1[0]-P[0],J1[1]-P[1])
            c=(P,(P[0]+ux*k*0.5,P[1]+uy*k*0.5),(J1[0]-dx*k,J1[1]-dy*k),J1)
            out.append(f'<path d="{d_(trim(c,P,False))}"{ms}/>')
        out.append(f'<path d="M{J1[0]},{J1[1]} L{J2[0]},{J2[1]}" class="stem"/>')
        for p in prods:
            P=S[p]; k=min(0.45*math.dist(P,J2),bend); ux,uy=unit(P[0]-J2[0],P[1]-J2[1])
            c=(J2,(J2[0]+dx*k,J2[1]+dy*k),(P[0]-ux*k*0.5,P[1]-uy*k*0.5),P)
            out.append(f'<path d="{d_(trim(c,P,True))}" marker-end="url(#ah)"/>')
    lab.append(f'<text x="{L[0]}" y="{L[1]}" class="rl" text-anchor="{L[2]}">{r}</text>')
nodes=''.join(f'<circle cx="{x}" cy="{y}" r="{rs}" class="{"hub" if m in "EG" else "sp"}"/><text x="{x}" y="{y+6}" class="sl">{m}</text>' for m,(x,y) in S.items())
body='<g class="edges">'+''.join(out)+'</g>'+''.join(lab)+nodes
open('net_reduced.svg','w').write(body)
style='.edges path{fill:none;stroke:#4B5563;stroke-width:1.8}.edges .stem{stroke-width:3}.sp{fill:#E3EEF9;stroke:#377EB8;stroke-width:1.8}.hub{fill:#FDE7D3;stroke:#D9731A;stroke-width:1.8}.sl{font:600 18px monospace;text-anchor:middle}.rl{font:13px monospace;fill:#374151;paint-order:stroke;stroke:white;stroke-width:4px}'
open('langary_network_reduced.svg','w').write(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="60 60 880 610" width="880" height="610"><style>{style}</style><rect x="60" y="60" width="880" height="610" fill="white"/><defs><marker id="ah" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0,1L9,5L0,9z" fill="#4B5563"/></marker></defs>{body}</svg>')
print('ok')
