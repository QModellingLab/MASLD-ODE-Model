#!/usr/bin/env python3
"""
Fig 4 v11 — uniform-width 1×6 stack
  - Intra: individual node-to-node, red
  - Inter: individual node-to-node, purple, via shared right/left bus channel
  - P_*:   individual node-to-P_*, gold dashed, via shared far-right channel
  - Inh:   individual, blue T-bar
  ALL lines connect to specific nodes — no floating.
"""
from collections import defaultdict

DISP = {
    "IL_6":"IL-6","IL_6R":"IL-6R","SOCS3":"SOCS3","TNFa":"TNFα",
    "TNFR1":"TNFR1","INS":"INS","INSR":"INSR","LEP":"LEP",
    "ObR":"ObR","ACDC":"ACDC","adipoR":"AdipoR",
    "IRS_1_2":"IRS-1/2","PI3K":"PI3K","Akt":"Akt","GSK_3":"GSK-3",
    "LXR_a":"LXRα","RXR":"RXR","SREBP_1c":"SREBP-1c",
    "ChREBP":"ChREBP","L_PK":"L-PK",
    "AMPK":"AMPK","p38":"p38","PPAR_a":"PPARα","PPAR_g":"PPARγ",
    "IRE1a":"IRE1α","XBP1":"XBP1",
    "PERK":"PERK","eIF2a":"eIF2α","ATF4":"ATF4","CHOP":"CHOP",
    "CYP2E1":"CYP2E1","CxI":"Cx I","CxII":"Cx II","CxIII":"Cx III","CxIV":"Cx IV",
    "IKKb":"IKKβ","NF_kB":"NF-κB","JNK1_2":"JNK1/2",
    "AP_1":"AP-1","ASK1":"ASK1","TRAF2":"TRAF2",
    "IL_1":"IL-1","IL_8":"IL-8",
    "TGF_b1":"TGF-β1","FasL":"FasL","C_EBPa":"C/EBPα",
    "Cdc42":"Cdc42","Rac1":"Rac1","MLK3":"MLK3",
    "ITCH":"ITCH","Bim":"Bim","Bax":"Bax","Bid":"Bid",
    "CASP8":"CASP8","CASP3":"CASP3","CASP7":"CASP7",
    "Cytc":"Cyt c","Fas":"Fas",
}
MODULES = [
    ("① Adipocytokine / Receptor","#F5F5F5","#9E9E9E","#555"),
    ("② Insulin / PI3K-Akt","#EDE7F6","#7E57C2","#4A148C"),
    ("③ AMPK / Lipid Metabolism","#E8F5E9","#66BB6A","#1B5E20"),
    ("④ ER Stress / Mitochondria","#FFF3E0","#FF9800","#E65100"),
    ("⑤ Inflammation / NF-κB","#FFEBEE","#EF5350","#B71C1C"),
    ("⑥ Apoptosis / Cell Death","#E3F2FD","#42A5F5","#0D47A1"),
]
MOD_NODES = [
    ["IL_6","IL_6R","SOCS3","TNFa","TNFR1","INS","INSR","LEP","ObR","ACDC","adipoR"],
    ["IRS_1_2","PI3K","Akt","GSK_3","LXR_a","RXR","SREBP_1c","ChREBP","L_PK"],
    ["AMPK","p38","PPAR_a","PPAR_g"],
    ["IRE1a","XBP1","PERK","eIF2a","ATF4","CHOP","CYP2E1","CxI","CxII","CxIII","CxIV"],
    ["TRAF2","ASK1","IKKb","NF_kB","JNK1_2","AP_1","IL_1","IL_8","TGF_b1","FasL","C_EBPa"],
    ["Cdc42","Rac1","MLK3","ITCH","Fas","CASP8","Bid","Bim","Bax","Cytc","CASP3","CASP7"],
]
node2mod={}
for mi,ns in enumerate(MOD_NODES):
    for n in ns: node2mod[n]=mi

NW,NH = 58,21; FS=8.5; sp=76; vs=27; GAP=30; MPAD=10; LBL_H=13; SX=20
STACK_W = sp*5 + NW + MPAD*2

pos={}
bx = SX + MPAD
y_cursor = 18

def place_rows(nodes_per_row, y_start):
    y = y_start + LBL_H + MPAD
    for row in nodes_per_row:
        for ci,n in enumerate(row):
            pos[n] = (bx + ci*sp, y)
        y += vs
    return y - vs + NH + MPAD  # module bottom

# M1
m_bots = []
bot = place_rows([
    ["IL_6","IL_6R","SOCS3","TNFa","TNFR1"],
    ["INS","INSR","LEP","ObR"],
    ["ACDC","adipoR"],
], y_cursor)
m_bots.append(bot); m_tops = [y_cursor]; y_cursor = bot + GAP

# M2
m_tops.append(y_cursor)
bot = place_rows([
    ["IRS_1_2","PI3K","Akt","GSK_3"],
    ["LXR_a","RXR","SREBP_1c","ChREBP","L_PK"],
], y_cursor)
m_bots.append(bot); y_cursor = bot + GAP

# M3
m_tops.append(y_cursor)
bot = place_rows([["AMPK","p38","PPAR_a","PPAR_g"]], y_cursor)
m_bots.append(bot); y_cursor = bot + GAP

# M4
m_tops.append(y_cursor)
y4 = y_cursor + LBL_H + MPAD
pos["IRE1a"]=(bx,y4); pos["XBP1"]=(bx+sp,y4); pos["PERK"]=(bx+sp*2,y4); pos["eIF2a"]=(bx+sp*3,y4)
pos["ATF4"]=(bx,y4+vs); pos["CHOP"]=(bx+sp,y4+vs); pos["CYP2E1"]=(bx+sp*2,y4+vs)
cx2=int((STACK_W-MPAD*2-NW)/4)
pos["CxI"]=(bx,y4+vs*2); pos["CxII"]=(bx+cx2,y4+vs*2); pos["CxIII"]=(bx+cx2*2,y4+vs*2); pos["CxIV"]=(bx+cx2*3,y4+vs*2)
bot = y4+vs*2+NH+MPAD; m_bots.append(bot); y_cursor = bot + GAP

# M5
m_tops.append(y_cursor)
bot = place_rows([
    ["TRAF2","ASK1","IKKb","NF_kB"],
    ["JNK1_2","AP_1","IL_1","IL_8"],
    ["FasL","C_EBPa","TGF_b1"],
], y_cursor)
m_bots.append(bot); y_cursor = bot + GAP

# M6
m_tops.append(y_cursor)
bot = place_rows([
    ["Cdc42","Rac1","MLK3","ITCH"],
    ["Fas","CASP8","Bid","Bim"],
    ["Bax","Cytc","CASP3","CASP7"],
], y_cursor)
m_bots.append(bot)

H = int(bot + 38)

# P_*
P_OUT=[("P_Hyperinsulinemia","Hyperinsulinemia"),("P_De_novo_fatty_acid_synthesis","De novo FA synthesis"),
("P_Improvement_of_NAFLD","Improvement of NAFLD"),("P_Development_of_NAFLD","Development of NAFLD"),
("P_Adipogenesis","Adipogenesis"),("P_HCC_proliferation","HCC proliferation"),
("P_Apoptosis","Apoptosis"),("P_Development_of_steatohepatitis","Steatohepatitis"),
("P_Cell_death","Cell death"),("P_Inflammation","Inflammation"),
("P_Fibrosis","Fibrosis"),("P_Hepatocyte_injury","Hepatocyte injury")]
PW,PH=115,21; px=SX+STACK_W+105; pv=30
p_h=(len(P_OUT)-1)*pv+PH; p_top=max(50,(bot-p_h)//2)
for i,(pn,_) in enumerate(P_OUT): pos[pn]=(px,p_top+i*pv)
W=px+PW+15

# ============================================================
ACTS=[("ACDC","adipoR"),("AMPK","P_Improvement_of_NAFLD"),("AMPK","p38"),("AP_1","IL_1"),("AP_1","IL_6"),("AP_1","TNFa"),("ASK1","JNK1_2"),("ATF4","CHOP"),("Bax","Cytc"),("Bid","Bax"),("Bim","Bax"),("CASP3","P_Hepatocyte_injury"),("CASP7","P_Hepatocyte_injury"),("CASP8","Bid"),("CHOP","Bim"),("CYP2E1","FasL"),("CYP2E1","IKKb"),("CYP2E1","IL_8"),("CYP2E1","JNK1_2"),("CYP2E1","TGF_b1"),("CYP2E1","TNFa"),("C_EBPa","P_Adipogenesis"),("Cdc42","MLK3"),("ChREBP","L_PK"),("ChREBP","P_De_novo_fatty_acid_synthesis"),("CxI","FasL"),("CxI","IKKb"),("CxI","IL_8"),("CxI","JNK1_2"),("CxI","TGF_b1"),("CxI","TNFa"),("CxII","FasL"),("CxII","IKKb"),("CxII","IL_8"),("CxII","JNK1_2"),("CxII","TGF_b1"),("CxII","TNFa"),("CxIII","FasL"),("CxIII","IKKb"),("CxIII","IL_8"),("CxIII","JNK1_2"),("CxIII","TGF_b1"),("CxIII","TNFa"),("CxIV","FasL"),("CxIV","IKKb"),("CxIV","IL_8"),("CxIV","JNK1_2"),("CxIV","TGF_b1"),("CxIV","TNFa"),("Cytc","CASP3"),("Cytc","CASP7"),("Fas","CASP8"),("FasL","Fas"),("FasL","P_Cell_death"),("GSK_3","INS"),("GSK_3","P_Hyperinsulinemia"),("IKKb","NF_kB"),("IL_1","P_Development_of_steatohepatitis"),("IL_6","IL_6R"),("IL_6","P_Development_of_steatohepatitis"),("IL_6R","SOCS3"),("IL_8","P_Inflammation"),("INS","INSR"),("INSR","IRS_1_2"),("INSR","LXR_a"),("IRE1a","TRAF2"),("IRE1a","XBP1"),("IRS_1_2","PI3K"),("ITCH","CASP8"),("JNK1_2","AP_1"),("JNK1_2","ITCH"),("JNK1_2","P_Apoptosis"),("JNK1_2","P_HCC_proliferation"),("LEP","ObR"),("LXR_a","RXR"),("MLK3","JNK1_2"),("NF_kB","IL_1"),("NF_kB","IL_6"),("NF_kB","TNFa"),("ObR","AMPK"),("PERK","eIF2a"),("PI3K","Akt"),("PPAR_a","P_Improvement_of_NAFLD"),("PPAR_g","P_Development_of_NAFLD"),("RXR","PPAR_g"),("RXR","SREBP_1c"),("Rac1","MLK3"),("SOCS3","SREBP_1c"),("SREBP_1c","P_De_novo_fatty_acid_synthesis"),("TGF_b1","P_Fibrosis"),("TGF_b1","P_Inflammation"),("TNFR1","NF_kB"),("TNFa","JNK1_2"),("TNFa","P_Cell_death"),("TNFa","P_Development_of_steatohepatitis"),("TNFa","TNFR1"),("TRAF2","ASK1"),("TRAF2","IKKb"),("XBP1","C_EBPa"),("XBP1","P_De_novo_fatty_acid_synthesis"),("adipoR","AMPK"),("eIF2a","ATF4"),("p38","PPAR_a")]
INHS=[("Akt","GSK_3"),("JNK1_2","IRS_1_2"),("SOCS3","IRS_1_2")]

# Classify
intra_a=[]; inter_a=[]; p_a=[]
for s,t in ACTS:
    if t.startswith("P_"): p_a.append((s,t))
    elif node2mod[s]==node2mod[t]: intra_a.append((s,t))
    else: inter_a.append((s,t))
intra_i=[(s,t) for s,t in INHS if node2mod[s]==node2mod[t]]
inter_i=[(s,t) for s,t in INHS if node2mod[s]!=node2mod[t]]
print(f"Intra:{len(intra_a)}+{len(intra_i)}, Inter:{len(inter_a)}+{len(inter_i)}, P_*:{len(p_a)}, Total:{len(ACTS)+len(INHS)}")

# Routing
def nc(n):
    x,y=pos[n]; w=PW if n.startswith("P_") else NW; h=PH if n.startswith("P_") else NH
    return (x+w/2,y+h/2)
def pt(n,side):
    cx,cy=nc(n); hw=(PW if n.startswith("P_") else NW)/2; hh=(PH if n.startswith("P_") else NH)/2
    return {'r':(cx+hw,cy),'l':(cx-hw,cy),'t':(cx,cy-hh),'b':(cx,cy+hh)}[side]
def auto_s(s,t):
    sx,sy=nc(s);tx,ty=nc(t);dx=tx-sx;dy=ty-sy
    if abs(dx)>=abs(dy): return ('r','l') if dx>0 else ('l','r')
    return ('b','t') if dy>0 else ('t','b')
def ortho(s,t,off=0):
    ss,ts=auto_s(s,t); x1,y1=pt(s,ss); x2,y2=pt(t,ts)
    if ss in 'rl' and ts in 'rl':
        mx=(x1+x2)/2+off
        return f"M{x1:.0f},{y1:.0f} L{mx:.0f},{y1:.0f} L{mx:.0f},{y2:.0f} L{x2:.0f},{y2:.0f}"
    elif ss in 'tb' and ts in 'tb':
        my=(y1+y2)/2+off
        return f"M{x1:.0f},{y1:.0f} L{x1:.0f},{my:.0f} L{x2:.0f},{my:.0f} L{x2:.0f},{y2:.0f}"
    else:
        if ss in 'rl': return f"M{x1:.0f},{y1:.0f} L{x2:.0f},{y1:.0f} L{x2:.0f},{y2:.0f}"
        else: return f"M{x1:.0f},{y1:.0f} L{x1:.0f},{y2:.0f} L{x2:.0f},{y2:.0f}"

stack_right = SX + STACK_W

# Inter-module routing: all use RIGHT bus for downward, LEFT bus for upward
# Allocate unique bus-x slot per edge to avoid overlap
inter_down = [(s,t) for s,t in inter_a if nc(t)[1]>nc(s)[1]]
inter_up   = [(s,t) for s,t in inter_a if nc(t)[1]<=nc(s)[1]]
inter_i_down = [(s,t) for s,t in inter_i if nc(t)[1]>nc(s)[1]]
inter_i_up   = [(s,t) for s,t in inter_i if nc(t)[1]<=nc(s)[1]]

# ============================================================
# SVG
# ============================================================
svg=[]
svg.append(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" font-family="Arial, Helvetica, sans-serif">')
svg.append('''<defs>
  <marker id="a" markerWidth="6" markerHeight="5" refX="6" refY="2.5" orient="auto" markerUnits="userSpaceOnUse"><path d="M0,0.5 L5.5,2.5 L0,4.5" fill="none" stroke="#C62828" stroke-width="0.7" stroke-linejoin="round"/></marker>
  <marker id="i" markerWidth="2" markerHeight="9" refX="1" refY="4.5" orient="auto" markerUnits="userSpaceOnUse"><line x1="1" y1="0" x2="1" y2="9" stroke="#1565C0" stroke-width="1.8"/></marker>
  <marker id="x" markerWidth="6" markerHeight="5" refX="6" refY="2.5" orient="auto" markerUnits="userSpaceOnUse"><path d="M0,0.5 L5.5,2.5 L0,4.5" fill="none" stroke="#8E24AA" stroke-width="0.6"/></marker>
  <marker id="xi" markerWidth="2" markerHeight="9" refX="1" refY="4.5" orient="auto" markerUnits="userSpaceOnUse"><line x1="1" y1="0" x2="1" y2="9" stroke="#8E24AA" stroke-width="1.6"/></marker>
  <marker id="p" markerWidth="5" markerHeight="4" refX="5" refY="2" orient="auto" markerUnits="userSpaceOnUse"><path d="M0,0.3 L4.5,2 L0,3.7" fill="none" stroke="#C49000" stroke-width="0.7"/></marker>
</defs>''')
svg.append(f'<rect width="{W}" height="{H}" fill="white"/>')

# Module bgs
for mi,(ml,mf,ms,mt) in enumerate(MODULES):
    svg.append(f'<rect x="{SX}" y="{m_tops[mi]}" width="{STACK_W}" height="{m_bots[mi]-m_tops[mi]}" rx="5" fill="{mf}" stroke="{ms}" stroke-width="1.3" opacity="0.8"/>')
    svg.append(f'<text x="{SX+5}" y="{m_tops[mi]+11}" font-size="7.5" font-weight="bold" fill="{mt}">{ml}</text>')
# P_* bg
px0=px-8;py0=p_top-10-LBL_H;pw0=PW+16;ph0=(len(P_OUT)-1)*pv+PH+22+LBL_H
svg.append(f'<rect x="{px0}" y="{py0}" width="{pw0}" height="{ph0}" rx="5" fill="#FFFDE7" stroke="#C49000" stroke-width="1.3" opacity="0.8"/>')
svg.append(f'<text x="{px0+5}" y="{py0+11}" font-size="7.5" font-weight="bold" fill="#F57F17">Pathway Outputs</text>')

# ===== EDGES =====

# 1) Intra-module (red, node-to-node)
oc={}
for s,t in intra_a:
    k=tuple(sorted([s,t])); oc[k]=oc.get(k,0); o=oc[k]*3; oc[k]+=1
    svg.append(f'<path d="{ortho(s,t,o)}" fill="none" stroke="#C62828" stroke-width="0.65" marker-end="url(#a)" opacity="0.55"/>')
for s,t in intra_i:
    svg.append(f'<path d="{ortho(s,t)}" fill="none" stroke="#1565C0" stroke-width="1.0" marker-end="url(#i)"/>')

# 2) Inter-module downward (purple, right bus, node-to-node)
for ei,(s,t) in enumerate(inter_down):
    bx_i = stack_right + 6 + ei*1.2
    sx,sy = pt(s,'r'); tx,ty = pt(t,'r')
    svg.append(f'<path d="M{sx:.0f},{sy:.0f} L{bx_i:.0f},{sy:.0f} L{bx_i:.0f},{ty:.0f} L{tx:.0f},{ty:.0f}" '
               f'fill="none" stroke="#8E24AA" stroke-width="0.5" marker-end="url(#x)" opacity="0.32"/>')

# 3) Inter-module upward / feedback (purple, left bus, node-to-node)
for ei,(s,t) in enumerate(inter_up):
    lx_i = SX - 5 - ei*1.5
    sx,sy = pt(s,'l'); tx,ty = pt(t,'l')
    svg.append(f'<path d="M{sx:.0f},{sy:.0f} L{lx_i:.0f},{sy:.0f} L{lx_i:.0f},{ty:.0f} L{tx:.0f},{ty:.0f}" '
               f'fill="none" stroke="#8E24AA" stroke-width="0.5" marker-end="url(#x)" opacity="0.32"/>')

# Inter inhibitions
for ei,(s,t) in enumerate(inter_i_down + inter_i_up):
    is_down = nc(t)[1] > nc(s)[1]
    if is_down:
        bx_i = stack_right + 6 + (len(inter_down)+ei)*1.2
        sx,sy=pt(s,'r'); tx,ty=pt(t,'r')
        svg.append(f'<path d="M{sx:.0f},{sy:.0f} L{bx_i:.0f},{sy:.0f} L{bx_i:.0f},{ty:.0f} L{tx:.0f},{ty:.0f}" '
                   f'fill="none" stroke="#8E24AA" stroke-width="0.7" marker-end="url(#xi)" opacity="0.45"/>')
    else:
        lx_i = SX - 5 - (len(inter_up)+ei)*1.5
        sx,sy=pt(s,'l'); tx,ty=pt(t,'l')
        svg.append(f'<path d="M{sx:.0f},{sy:.0f} L{lx_i:.0f},{sy:.0f} L{lx_i:.0f},{ty:.0f} L{tx:.0f},{ty:.0f}" '
                   f'fill="none" stroke="#8E24AA" stroke-width="0.7" marker-end="url(#xi)" opacity="0.45"/>')

# 4) P_* edges (gold dashed, right channel, node-to-P_*)
p_chan_base = stack_right + 6 + (len(inter_down)+len(inter_i))*1.2 + 10
for ei,(s,t) in enumerate(sorted(p_a, key=lambda e:pos[e[1]][1])):
    cx = p_chan_base + ei*1.5
    sx,sy = pt(s,'r'); tx,ty = pt(t,'l')
    svg.append(f'<path d="M{sx:.0f},{sy:.0f} L{cx:.0f},{sy:.0f} L{cx:.0f},{ty:.0f} L{tx:.0f},{ty:.0f}" '
               f'fill="none" stroke="#C49000" stroke-width="0.45" stroke-dasharray="3.5,2" marker-end="url(#p)" opacity="0.40"/>')

# ===== NODES =====
for mi,(ml,mf,ms,mt) in enumerate(MODULES):
    for n in MOD_NODES[mi]:
        x,y=pos[n]; label=DISP.get(n,n)
        svg.append(f'<rect x="{x}" y="{y}" width="{NW}" height="{NH}" rx="3" fill="white" stroke="{ms}" stroke-width="1.2"/>')
        svg.append(f'<text x="{x+NW/2}" y="{y+NH/2+3}" text-anchor="middle" font-size="{FS}" font-weight="bold" fill="#212121">{label}</text>')
for pn,pl in P_OUT:
    x,y=pos[pn]
    svg.append(f'<rect x="{x}" y="{y}" width="{PW}" height="{PH}" rx="3" fill="#FFFDE7" stroke="#C49000" stroke-width="1.0" stroke-dasharray="4,2"/>')
    svg.append(f'<text x="{x+PW/2}" y="{y+PH/2+3}" text-anchor="middle" font-size="7.5" fill="#996515">{pl}</text>')

# Legend
lx,ly=20,H-28
svg.append(f'<rect x="{lx}" y="{ly}" width="540" height="22" rx="3" fill="#FAFAFA" stroke="#BDBDBD" stroke-width="0.7"/>')
for ix,color,sw,mk,dash,op,label in [
    (6,"#C62828","0.65","a","","","Intra-module"),(98,"#1565C0","1.0","i","","","Inhibition"),
    (185,"#8E24AA","0.5","x","","opacity=\"0.35\" ","Inter-module"),(290,"#C49000","0.45","p","3.5,2","","Output readout"),
]:
    x0=lx+ix; da=f'stroke-dasharray="{dash}" ' if dash else ""
    svg.append(f'<line x1="{x0}" y1="{ly+12}" x2="{x0+18}" y2="{ly+12}" stroke="{color}" stroke-width="{sw}" {da}{op}marker-end="url(#{mk})"/>')
    svg.append(f'<text x="{x0+22}" y="{ly+15}" font-size="7" fill="#333">{label}</text>')
svg.append(f'<text x="{W-8}" y="{H-5}" text-anchor="end" font-size="6.5" fill="#AAA">70 nodes · 106 interactions · KEGG hsa04932</text>')
svg.append('</svg>')

with open('/home/claude/fig4_v11.svg','w') as f: f.write('\n'.join(svg))

# Verify: every node has at least 1 path touching it
node_touched=set()
for s,t in ACTS: node_touched.add(s); node_touched.add(t)
for s,t in INHS: node_touched.add(s); node_touched.add(t)
all_nodes=set()
for ns in MOD_NODES: all_nodes.update(ns)
untouched=all_nodes-node_touched
print(f"Canvas: {W}x{H}")
print(f"Paths drawn: {len(intra_a)+len(intra_i)+len(inter_down)+len(inter_up)+len(inter_i)+len(p_a)}")
print(f"Untouched nodes: {len(untouched)} {'— '+str(untouched) if untouched else '✅ ALL connected'}")
