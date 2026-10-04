#!/usr/bin/env python3
"""Write the synthetic test boards (kicad_fdtd geometry JSON) and the case list.

    python scripts/make_fixtures.py            # writes tests/fixtures/*.json

The boards are small enough for the suite to run in seconds but exercise everything the
solver touches: 2/4/6-layer stackups, heavy copper, split planes, plane cutouts that
overlap (an invalid polygon shapely repairs), pours, tracks, vias, a non-rectangular
outline with air cells, a board with no stackup (default FR4 fallback), and DC paths
through narrow necks and through vias between layers. Two more are the analytic
cases from kicad_fdtd's own tests (copper plate convection, track IR drop).
"""
import json
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))
FIX = os.path.join(HERE, "..", "tests", "fixtures")
rect = lambda xa, ya, xb, yb: [[xa, ya], [xb, ya], [xb, yb], [xa, yb]]
X0, Y0 = 100.0, 50.0


def stack(layers, diel, cu_t):
    st = []
    for i, l in enumerate(layers):
        st.append({"name": l, "type": "copper", "thickness_mm": cu_t[i], "er": None, "tand": None})
        if i < len(layers) - 1:
            st.append({"name": "d%d" % i, "type": "prepreg", "thickness_mm": diel[i], "er": 4.5, "tand": 0.02})
    return st


def board(W, H, layers, diel, cu_t, n_tracks, n_vias, n_holes, n_parts, pours, neck, seed):
    """Random board with a /PWR pour-neck-pour DC path on the top layer (S.1 -> K.1)."""
    random.seed(seed)
    x0, y0 = X0, Y0
    cu = {l: [] for l in layers}
    inner = layers[1:-1]
    holes = [rect(x0 + x, y0 + y, x0 + x + 0.6, y0 + y + 0.6)
             for x, y in ((random.uniform(2, W - 3), random.uniform(2, H - 3)) for _ in range(n_holes))]
    full = (x0 + 0.5, y0 + 0.5, x0 + W - 0.5, y0 + H - 0.5)
    if inner:
        cu[inner[0]].append({"net": "GND", "kind": "zone", "outline": rect(*full), "holes": holes})
        for i, l in enumerate(inner[1:]):
            half = (x0 + 0.5, y0 + 0.5, x0 + W / 2 - 0.5, y0 + H - 0.5) if i % 2 == 0 else (x0 + W / 2 + 0.5, y0 + 0.5, x0 + W - 0.5, y0 + H - 0.5)
            cu[l].append({"net": "/12V" if i % 2 == 0 else "GND", "kind": "zone", "outline": rect(*half), "holes": holes[:n_holes // 2]})
    else:
        cu[layers[-1]].append({"net": "GND", "kind": "zone", "outline": rect(*full), "holes": holes})
    for l in (layers[0], layers[-1]):
        for _ in range(pours):
            pw, ph = min(25, W / 3), min(20, H / 2)
            xa, ya = random.uniform(1, W - pw - 1), random.uniform(1, H - ph - 1)
            cu[l].append({"net": random.choice(["GND", "/VBAT", "/12V"]), "kind": "zone", "outline": rect(x0 + xa, y0 + ya, x0 + xa + pw, y0 + ya + ph), "holes": []})
        for i in range(n_tracks):
            ln, w = random.uniform(3, min(18, W / 3)), random.choice([0.2, 0.25, 0.4, 1.0])
            xa, ya = random.uniform(1, W - ln - 1), random.uniform(1, H - w - 1)
            cu[l].append({"net": "/S%d" % i, "kind": "track", "outline": rect(x0 + xa, y0 + ya, x0 + xa + ln, y0 + ya + w), "holes": []})
    bw, nw, nl = neck
    yc = y0 + H / 2
    cu[layers[0]] += [{"net": "/PWR", "kind": "zone", "outline": rect(x0 + 2, yc - bw / 2, x0 + 2 + bw, yc + bw / 2), "holes": []},
                      {"net": "/PWR", "kind": "track", "outline": rect(x0 + 2 + bw, yc - nw / 2, x0 + 2 + bw + nl, yc + nw / 2), "holes": []},
                      {"net": "/PWR", "kind": "zone", "outline": rect(x0 + 2 + bw + nl, yc - bw / 2, x0 + 2 + 2 * bw + nl, yc + bw / 2), "holes": []}]
    pads = [{"ref": "S", "num": "1", "net": "/PWR", "x": x0 + 2 + bw * 0.3, "y": yc, "layers": [layers[0]]},
            {"ref": "K", "num": "1", "net": "/PWR", "x": x0 + 2 + bw * 1.7 + nl, "y": yc, "layers": [layers[0]]}]
    fps = [{"ref": "S", "value": "src", "x": x0 + 2 + bw * 0.3, "y": yc, "layer": layers[0]},
           {"ref": "K", "value": "sink", "x": x0 + 2 + bw * 1.7 + nl, "y": yc, "layer": layers[0]}]
    vias = [{"x": x0 + random.uniform(1, W - 1), "y": y0 + random.uniform(1, H - 1), "drill": 0.3, "dia": 0.6,
             "top": layers[0], "bottom": layers[-1], "net": random.choice(["GND", "/12V"]), "kind": "via"} for _ in range(n_vias)]
    for k in range(n_parts):
        xc, yc2 = x0 + random.uniform(5, W - 8), y0 + random.uniform(3, H - 6)
        fps.append({"ref": "U%d" % k, "value": "IC%d" % k, "x": xc, "y": yc2, "layer": layers[0]})
        for p in range(8):
            pads.append({"ref": "U%d" % k, "num": str(p + 1), "net": "/12V", "x": xc + (p % 4) * 1.0, "y": yc2 + (p // 4) * 2.0, "layers": [layers[0]]})
    return {"board": "synthetic", "kicad": None, "layers": layers, "bbox": [x0, y0, W, H], "outline": rect(x0, y0, x0 + W, y0 + H),
            "stackup": stack(layers, diel, cu_t), "copper": cu, "vias": vias, "pads": pads, "footprints": fps,
            "nets": ["GND", "/12V", "/VBAT", "/PWR"] + ["/S%d" % i for i in range(n_tracks)]}


def lshape_vias():
    """L-shaped 2-layer board, no stackup in the file, DC from a top pad to a bottom pad
    through four vias, plus a through-hole connector (pth barrel)."""
    x0, y0, W, H = X0, Y0, 30.0, 30.0
    L = ["F.Cu", "B.Cu"]
    outline = [[x0, y0], [x0 + W, y0], [x0 + W, y0 + 12], [x0 + 12, y0 + 12], [x0 + 12, y0 + H], [x0, y0 + H]]
    cu = {"F.Cu": [{"net": "/VIN", "kind": "zone", "outline": rect(x0 + 1, y0 + 1, x0 + 10, y0 + 10), "holes": []},
                   {"net": "/VIN", "kind": "pad", "outline": rect(x0 + 3, y0 + 3, x0 + 5, y0 + 5), "holes": []},
                   {"net": "GND", "kind": "zone", "outline": rect(x0 + 13, y0 + 1, x0 + 29, y0 + 11), "holes": []},
                   {"net": "/SIG", "kind": "track", "outline": rect(x0 + 2, y0 + 15, x0 + 10, y0 + 15.3), "holes": []}],
          "B.Cu": [{"net": "/VIN", "kind": "zone", "outline": rect(x0 + 6, y0 + 1, x0 + 11, y0 + 28), "holes": [rect(x0 + 8, y0 + 14, x0 + 9, y0 + 16)]},
                   {"net": "GND", "kind": "zone", "outline": rect(x0 + 1, y0 + 13, x0 + 5, y0 + 29), "holes": []}]}
    vias = [{"x": x0 + 7 + 0.8 * i, "y": y0 + 8, "drill": 0.3, "dia": 0.6, "top": "F.Cu", "bottom": "B.Cu", "net": "/VIN", "kind": "via"} for i in range(4)]
    vias.append({"x": x0 + 20, "y": y0 + 6, "drill": 1.0, "dia": 1.05, "top": "F.Cu", "bottom": "B.Cu", "net": "GND", "kind": "pth"})
    pads = [{"ref": "J1", "num": "1", "net": "/VIN", "x": x0 + 4, "y": y0 + 4, "layers": ["F.Cu"]},
            {"ref": "U1", "num": "3", "net": "/VIN", "x": x0 + 8.5, "y": y0 + 26, "layers": ["B.Cu"]},
            {"ref": "U1", "num": "4", "net": "GND", "x": x0 + 3, "y": y0 + 26, "layers": ["B.Cu"]},
            {"ref": "Q1", "num": "1", "net": "GND", "x": x0 + 20, "y": y0 + 6, "layers": ["F.Cu", "B.Cu"]},
            {"ref": "Q1", "num": "2", "net": "GND", "x": x0 + 24, "y": y0 + 6, "layers": ["F.Cu"]},
            {"ref": "H1", "num": "", "net": "", "x": x0 + 26, "y": y0 + 9, "layers": []}]
    fps = [{"ref": "J1", "value": "XT30", "x": x0 + 4, "y": y0 + 4, "layer": "F.Cu"},
           {"ref": "U1", "value": "BUCK", "x": x0 + 6, "y": y0 + 26, "layer": "B.Cu"},
           {"ref": "Q1", "value": "FET", "x": x0 + 22, "y": y0 + 6, "layer": "F.Cu"},
           {"ref": "H1", "value": "MountingHole", "x": x0 + 26, "y": y0 + 9, "layer": "F.Cu"}]
    return {"board": "synthetic", "kicad": None, "layers": L, "bbox": [x0, y0, W, H], "outline": outline, "stackup": [],
            "copper": cu, "vias": vias, "pads": pads, "footprints": fps, "nets": ["/VIN", "GND", "/SIG"]}


def synth(layers, W=0.2, L=12.0):
    """kicad_fdtd tests/test_geometry.py: synth()."""
    x0, y0 = X0, Y0
    copper = {l: [] for l in layers}
    copper[layers[0]].append({"net": "/T", "kind": "track", "outline": rect(x0 + 4, y0 + 4 - W / 2, x0 + 4 + L, y0 + 4 + W / 2), "holes": []})
    copper[layers[0]].append({"net": "/CLK_P", "kind": "track", "outline": rect(x0 + 4, y0 + 2 - W / 2, x0 + 4 + L, y0 + 2 + W / 2), "holes": []})
    copper[layers[0]].append({"net": "/CLK_N", "kind": "track", "outline": rect(x0 + 4, y0 + 2.5 - W / 2, x0 + 4 + L, y0 + 2.5 + W / 2), "holes": []})
    copper[layers[1]].append({"net": "GND", "kind": "zone", "outline": rect(x0, y0, x0 + 20, y0 + 8), "holes": []})
    st = stack(layers, [1.53] * (len(layers) - 1), [0.035] * len(layers))
    return {"layers": layers, "bbox": [x0, y0, 20.0, 8.0], "outline": rect(x0, y0, x0 + 20, y0 + 8), "stackup": st, "copper": copper,
            "vias": [{"x": x0 + 10, "y": y0 + 6, "drill": 0.3, "dia": 0.35, "top": layers[0], "bottom": layers[-1], "net": "GND", "kind": "via"}],
            "pads": [{"ref": "P", "num": "1", "net": "/T", "x": x0 + 4.3, "y": y0 + 4, "layers": [layers[0]]},
                     {"ref": "P", "num": "2", "net": "/T", "x": x0 + 4 + L - 0.3, "y": y0 + 4, "layers": [layers[0]]},
                     {"ref": "J", "num": "1", "net": "/T", "x": x0 + 1, "y": y0 + 1, "layers": layers}],
            "footprints": [{"ref": "P", "value": "port", "x": x0 + 4, "y": y0 + 4, "layer": layers[0]}],
            "nets": ["/T", "/CLK_P", "/CLK_N", "GND"]}


def plate():
    """kicad_fdtd tests/test_thermal.py: test_plate_convection."""
    g = synth(["F.Cu", "B.Cu"])
    full = rect(X0, Y0, X0 + 20, Y0 + 8)
    g["copper"]["F.Cu"] = [{"net": "GND", "kind": "zone", "outline": full, "holes": []}]
    g["copper"]["B.Cu"] = [{"net": "GND", "kind": "zone", "outline": full, "holes": []}]
    return g


L4 = ["F.Cu", "In1.Cu", "In2.Cu", "B.Cu"]
L6 = ["F.Cu", "In1.Cu", "In2.Cu", "In3.Cu", "In4.Cu", "B.Cu"]
BOARDS = {
    "a4_power": board(40, 30, L4, [0.2, 1.065, 0.2], [0.035] * 4, 80, 150, 120, 8, 3, (8, 0.8, 6), seed=1),
    "b2_compact": board(30, 20, ["F.Cu", "B.Cu"], [1.53], [0.035] * 2, 80, 60, 60, 6, 2, (6, 0.6, 5), seed=2),
    "c6_heavy": board(30, 24, L6, [0.1, 0.5, 0.2, 0.5, 0.1], [0.035, 0.07, 0.07, 0.07, 0.07, 0.035], 60, 150, 100, 6, 2, (7, 1.0, 5), seed=3),
    "d2_strip": board(80, 12, ["F.Cu", "B.Cu"], [1.53], [0.07] * 2, 60, 80, 60, 5, 2, (5, 1.0, 30), seed=4),
    "e2_lshape_vias": lshape_vias(),
    "plate": plate(),
    "track": synth(["F.Cu", "B.Cu"]),
}

pw = lambda refs: [[r, round(0.05 + 0.03 * i, 4)] for i, r in enumerate(refs)]
CASES = [
    {"name": "a4_power", "cell": 0.5, "h": 8.0, "tamb": 40.0, "powers": pw(["U%d" % i for i in range(8)]),
     "dc": {"net": "/PWR", "from": ["S", "1"], "to": ["K", "1"], "current": 1.5}},
    {"name": "b2_compact", "cell": 0.25, "h": 8.0, "tamb": 40.0, "powers": pw(["U%d" % i for i in range(6)]),
     "dc": {"net": "/PWR", "from": ["S", "1"], "to": ["K", "1"], "current": 1.0}},
    {"name": "c6_heavy", "cell": 0.5, "h": 8.0, "tamb": 40.0, "powers": pw(["U%d" % i for i in range(6)]),
     "dc": {"net": "/PWR", "from": ["S", "1"], "to": ["K", "1"], "current": 2.0}},
    {"name": "d2_strip", "cell": 0.4, "h": 8.0, "tamb": 40.0, "powers": pw(["U%d" % i for i in range(5)]),
     "dc": {"net": "/PWR", "from": ["S", "1"], "to": ["K", "1"], "current": 6.0}},
    {"name": "e2_lshape_vias", "cell": 0.3, "h": 10.0, "tamb": 35.0, "powers": [["U1", 0.8], ["Q1", 0.5]],
     "dc": {"net": "/VIN", "from": ["J1", "1"], "to": ["U1", "3"], "current": 3.0}},
    {"name": "plate", "cell": 0.5, "h": 10.0, "tamb": 25.0, "powers": [["P", 1.0]], "dc": None},
    {"name": "track", "cell": 0.1, "h": 10.0, "tamb": 25.0, "powers": [],
     "dc": {"net": "/T", "from": ["P", "1"], "to": ["P", "2"], "current": 1.0}},
]

if __name__ == "__main__":
    os.makedirs(FIX, exist_ok=True)
    for name, g in BOARDS.items():
        with open(os.path.join(FIX, name + ".json"), "w") as f:
            json.dump(g, f)
    with open(os.path.join(FIX, "cases.json"), "w") as f:
        json.dump(CASES, f, indent=1)
    print("wrote %d boards and %d cases to %s" % (len(BOARDS), len(CASES), os.path.normpath(FIX)))
