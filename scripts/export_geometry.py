#!/usr/bin/env python3
"""export_geometry.py: write a KiCad board's copper geometry as the JSON that kicad-thermal reads.

    python3 scripts/export_geometry.py BOARD.kicad_pcb OUT.json [--fill] [--stackup STACKUP.json]

Runs in a Python that can `import pcbnew` (KiCad 9; on Linux the distro's python3). The format is the one
`src/geometry.rs` parses; coordinates are in mm in KiCad's frame (y down), and the loader flips y.

  layers      enabled copper layers, top to bottom
  bbox        [x, y, w, h] of the board edges
  outline     the board outline (first outline of KiCad's board polygon)
  copper      {layer: [{net, kind, outline, holes}]}, kind = pad | track | zone | via_pad
  vias        [{x, y, drill, dia, top, bottom, net, kind}], kind = via | pth (plated through-hole pads: the barrel)
  pads        [{ref, num, net, x, y, layers}]   (used to place sources and sinks)
  footprints  [{ref, value, x, y, layer}]
  stackup     [{name, type, thickness_mm, er, tand}] from the board's stackup, or from --stackup
  nets        every net name

Polygons come from KiCad's own geometry engine (maximum chord error 0.005 mm, error inside the shape):
pads use their effective polygon per copper layer, tracks and arcs their outline, zones their filled polygons.
Via pads are copper too: a track that ends inside a via's pad but outside its drill is connected in KiCad, so each
via's pad is added to the copper of its top and bottom layers (kind via_pad). A plated through-hole pad's copper is
in its pad polygons; its barrel is listed under vias with dia = drill + 0.05 mm (a 0.025 mm plated wall).

--fill refills every zone before export (use it when the board file may hold stale fills).
--stackup reads {"stackup": [...], "note": "..."} (or a bare list) and uses it instead of the board's own stackup;
the file's note is kept as stackup_note.
"""
import argparse
import json
import re
import sys

import pcbnew

MAXERR_MM = 0.005
PLATED_WALL_MM = 0.025


def stackup_from_file(path):
    """The board's (stackup ...) section, read from the file text (layer name, type, thickness, er, loss tangent)."""
    txt = open(path, encoding='utf-8').read()
    start = txt.find('(stackup')
    if start < 0:
        return []
    depth, end = 0, start
    for end in range(start, len(txt)):
        depth += {'(': 1, ')': -1}.get(txt[end], 0)
        if depth == 0:
            break
    body = txt[start:end + 1]
    out = []
    # each (layer "NAME" ...) block at depth 1 inside (stackup ...)
    for m in re.finditer(r'\(layer "([^"]+)"', body):
        depth, i = 0, m.start()
        for i in range(m.start(), len(body)):
            depth += {'(': 1, ')': -1}.get(body[i], 0)
            if depth == 0:
                break
        blk = body[m.start():i + 1]
        typ = re.search(r'\(type "([^"]+)"\)', blk)
        if not typ:
            continue
        num = lambda key: (float(x.group(1)) if (x := re.search(r'\(%s ([-\d.eE]+)\)' % key, blk)) else None)
        out.append({'name': m.group(1), 'type': typ.group(1), 'thickness_mm': num('thickness'),
                    'er': num('epsilon_r'), 'tand': num('loss_tangent')})
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('board')
    ap.add_argument('out')
    ap.add_argument('--fill', action='store_true', help='refill all zones before exporting')
    ap.add_argument('--stackup', help='stackup JSON to use instead of the board\'s own')
    a = ap.parse_args()

    b = pcbnew.LoadBoard(a.board)
    if b is None:
        sys.exit(f'pcbnew could not load {a.board}')
    if a.fill:
        pcbnew.ZONE_FILLER(b).Fill(b.Zones())
    mm = pcbnew.ToMM
    maxerr = pcbnew.FromMM(MAXERR_MM)
    cu = [(b.GetLayerName(l), l) for l in b.GetEnabledLayers().CuStack()]
    top, bottom = cu[0][0], cu[-1][0]

    def pts(chain):
        return [[round(mm(chain.CPoint(i).x), 4), round(mm(chain.CPoint(i).y), 4)] for i in range(chain.PointCount())]

    def polys(ps, net, kind):
        out = []
        for i in range(ps.OutlineCount()):
            ring = pts(ps.Outline(i))
            if len(ring) < 3:
                continue
            holes = [h for h in (pts(ps.Hole(i, k)) for k in range(ps.HoleCount(i))) if len(h) >= 3]
            out.append({'net': net, 'kind': kind, 'outline': ring, 'holes': holes})
        return out

    def xy(p):
        return round(mm(p.x), 4), round(mm(p.y), 4)

    copper = {name: [] for name, _ in cu}
    vias, pads, fps = [], [], []
    for fp in b.GetFootprints():
        x, y = xy(fp.GetPosition())
        fps.append({'ref': fp.GetReference(), 'value': fp.GetValue(), 'x': x, 'y': y, 'layer': b.GetLayerName(fp.GetLayer())})
        for pad in fp.Pads():
            net = pad.GetNetname()
            on = [name for name, l in cu if pad.IsOnLayer(l)]
            px, py = xy(pad.GetPosition())
            pads.append({'ref': fp.GetReference(), 'num': pad.GetNumber(), 'net': net, 'x': px, 'y': py, 'layers': on})
            for name, l in cu:
                if pad.IsOnLayer(l):
                    copper[name] += polys(pad.GetEffectivePolygon(l, pcbnew.ERROR_INSIDE), net, 'pad')
            drill = mm(pad.GetDrillSize().x)
            if drill > 0 and pad.GetAttribute() == pcbnew.PAD_ATTRIB_PTH:
                vias.append({'x': px, 'y': py, 'drill': drill, 'dia': round(drill + 2 * PLATED_WALL_MM, 4),
                             'top': top, 'bottom': bottom, 'net': net, 'kind': 'pth'})
    for t in b.GetTracks():
        cls = t.GetClass()
        if cls in ('PCB_TRACK', 'PCB_ARC'):
            ps = pcbnew.SHAPE_POLY_SET()
            t.TransformShapeToPolygon(ps, t.GetLayer(), 0, maxerr, pcbnew.ERROR_INSIDE)
            copper[b.GetLayerName(t.GetLayer())] += polys(ps, t.GetNetname(), 'track')
        elif cls == 'PCB_VIA':
            x, y = xy(t.GetPosition())
            vt, vb = b.GetLayerName(t.TopLayer()), b.GetLayerName(t.BottomLayer())
            vias.append({'x': x, 'y': y, 'drill': mm(t.GetDrillValue()), 'dia': mm(t.GetWidth(pcbnew.F_Cu)),
                         'top': vt, 'bottom': vb, 'net': t.GetNetname(), 'kind': 'via'})
            for name, l in cu:
                if name in (vt, vb):
                    ps = pcbnew.SHAPE_POLY_SET()
                    t.TransformShapeToPolygon(ps, l, 0, maxerr, pcbnew.ERROR_INSIDE)
                    copper[name] += polys(ps, t.GetNetname(), 'via_pad')
    for z in b.Zones():
        if z.GetIsRuleArea():
            continue
        for name, l in cu:
            if z.IsOnLayer(l):
                copper[name] += polys(z.GetFilledPolysList(l), z.GetNetname(), 'zone')

    edge = pcbnew.SHAPE_POLY_SET()
    b.GetBoardPolygonOutlines(edge)           # KiCad 9 signature
    outline = pts(edge.Outline(0)) if edge.OutlineCount() else []
    bb = b.GetBoardEdgesBoundingBox()

    geo = {'board': a.board, 'kicad': pcbnew.GetBuildVersion(), 'layers': [n for n, _ in cu],
           'bbox': [mm(bb.GetX()), mm(bb.GetY()), mm(bb.GetWidth()), mm(bb.GetHeight())], 'outline': outline,
           'copper': copper, 'vias': vias, 'pads': pads, 'footprints': fps,
           'nets': sorted(n.GetNetname() for n in b.GetNetsByNetcode().values())}
    if a.stackup:
        st = json.load(open(a.stackup))
        geo['stackup'] = st['stackup'] if isinstance(st, dict) else st
        if isinstance(st, dict) and st.get('note'):
            geo['stackup_note'] = st['note']
    else:
        geo['stackup'] = stackup_from_file(a.board)
    with open(a.out, 'w') as f:
        json.dump(geo, f, sort_keys=True)
    n_cu = sum(geo['stackup'][i]['type'] == 'copper' for i in range(len(geo['stackup'])))
    print(f"{a.board}: polygons {{{', '.join(f'{k}: {len(v)}' for k, v in copper.items())}}}, "
          f"{sum(v['kind'] == 'via' for v in vias)} vias, {sum(v['kind'] == 'pth' for v in vias)} PTH barrels, "
          f"{len(pads)} pads, {len(fps)} footprints, stackup {len(geo['stackup'])} entries "
          f"({'complete' if n_cu == len(cu) else 'incomplete: the solver will use its default'}) -> {a.out}")


if __name__ == '__main__':
    main()
