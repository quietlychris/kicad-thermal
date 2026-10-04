"""Synthetic benchmark boards for kicad_fdtd's thermal/DC solver: writes $BENCH_BOARDS (default bench/boards)/<name>.json."""
import json, os, random
OUT = os.environ.get("BENCH_BOARDS") or os.path.join(os.path.dirname(os.path.abspath(__file__)), "boards")
os.makedirs(OUT, exist_ok=True)
rect = lambda xa, ya, xb, yb: [[xa, ya], [xb, ya], [xb, yb], [xa, yb]]

def make(name, W, H, layers, diel, cu_t, n_tracks, n_vias, n_holes, n_parts, pours, neck, seed=1):
    random.seed(seed); x0, y0 = 100.0, 50.0
    cu = {l: [] for l in layers}
    inner = layers[1:-1]
    holes = [rect(x0+x, y0+y, x0+x+0.6, y0+y+0.6) for x, y in ((random.uniform(2, W-3), random.uniform(2, H-3)) for _ in range(n_holes))]
    full = (x0+0.5, y0+0.5, x0+W-0.5, y0+H-0.5)
    if inner:
        cu[inner[0]].append({"net": "GND", "kind": "zone", "outline": rect(*full), "holes": holes})
        for i, l in enumerate(inner[1:]):
            half = (x0+0.5, y0+0.5, x0+W/2-0.5, y0+H-0.5) if i % 2 == 0 else (x0+W/2+0.5, y0+0.5, x0+W-0.5, y0+H-0.5)
            cu[l].append({"net": "/12V" if i % 2 == 0 else "GND", "kind": "zone", "outline": rect(*half), "holes": holes[:n_holes//2]})
    else:   # 2-layer: GND pour on the bottom
        cu[layers[-1]].append({"net": "GND", "kind": "zone", "outline": rect(*full), "holes": holes})
    for l in (layers[0], layers[-1]):
        for i in range(pours):
            pw, ph = min(25, W/3), min(20, H/2)
            xa, ya = random.uniform(1, W-pw-1), random.uniform(1, H-ph-1)
            cu[l].append({"net": random.choice(["GND", "/VBAT", "/12V"]), "kind": "zone", "outline": rect(x0+xa, y0+ya, x0+xa+pw, y0+ya+ph), "holes": []})
        for i in range(n_tracks):
            ln, w = random.uniform(3, min(18, W/3)), random.choice([0.2, 0.25, 0.4, 1.0])
            xa, ya = random.uniform(1, W-ln-1), random.uniform(1, H-w-1)
            cu[l].append({"net": "/S%d" % i, "kind": "track", "outline": rect(x0+xa, y0+ya, x0+xa+ln, y0+ya+w), "holes": []})
    # DC test path on the top layer: pour - neck - pour, net /PWR
    (bw, nw, nl) = neck
    yc = y0 + H/2
    cu[layers[0]] += [{"net": "/PWR", "kind": "zone", "outline": rect(x0+2, yc-bw/2, x0+2+bw, yc+bw/2), "holes": []},
                      {"net": "/PWR", "kind": "track", "outline": rect(x0+2+bw, yc-nw/2, x0+2+bw+nl, yc+nw/2), "holes": []},
                      {"net": "/PWR", "kind": "zone", "outline": rect(x0+2+bw+nl, yc-bw/2, x0+2+2*bw+nl, yc+bw/2), "holes": []}]
    pads = [{"ref": "S", "num": "1", "net": "/PWR", "x": x0+2+bw*0.3, "y": yc, "layers": [layers[0]]},
            {"ref": "K", "num": "1", "net": "/PWR", "x": x0+2+bw*1.7+nl, "y": yc, "layers": [layers[0]]}]
    fps = [{"ref": "S", "value": "src", "x": x0+2+bw*0.3, "y": yc, "layer": layers[0]},
           {"ref": "K", "value": "sink", "x": x0+2+bw*1.7+nl, "y": yc, "layer": layers[0]}]
    vias = [{"x": x0+random.uniform(1, W-1), "y": y0+random.uniform(1, H-1), "drill": 0.3, "dia": 0.6, "top": layers[0], "bottom": layers[-1],
             "net": random.choice(["GND", "/12V"]), "kind": "via"} for _ in range(n_vias)]
    for k in range(n_parts):
        xc, yc2 = x0+random.uniform(5, W-8), y0+random.uniform(3, H-6)
        fps.append({"ref": "U%d" % k, "value": "x", "x": xc, "y": yc2, "layer": layers[0]})
        for p in range(8):
            pads.append({"ref": "U%d" % k, "num": str(p+1), "net": "/12V", "x": xc+(p % 4)*1.0, "y": yc2+(p//4)*2.0, "layers": [layers[0]]})
    st = []
    for i, l in enumerate(layers):
        st.append({"name": l, "type": "copper", "thickness_mm": cu_t[i], "er": None, "tand": None})
        if i < len(layers)-1:
            st.append({"name": "d%d" % i, "type": "prepreg", "thickness_mm": diel[i], "er": 4.5, "tand": 0.02})
    g = {"layers": layers, "bbox": [x0, y0, W, H], "outline": rect(x0, y0, x0+W, y0+H), "stackup": st, "copper": cu, "vias": vias,
         "pads": pads, "footprints": fps, "nets": ["GND", "/12V", "/VBAT", "/PWR"] + ["/S%d" % i for i in range(n_tracks)]}
    json.dump(g, open(os.path.join(OUT, "%s.json" % name), "w"))

L4 = ["F.Cu", "In1.Cu", "In2.Cu", "B.Cu"]
L6 = ["F.Cu", "In1.Cu", "In2.Cu", "In3.Cu", "In4.Cu", "B.Cu"]
make("A_4L_power", 120, 90, L4, [0.2, 1.065, 0.2], [0.035]*4, 300, 600, 400, 15, 6, (20, 1.0, 15))
make("B_2L_compact", 60, 40, ["F.Cu", "B.Cu"], [1.53], [0.035]*2, 250, 200, 150, 10, 3, (10, 0.6, 8), seed=2)
make("C_6L_heavy", 100, 80, L6, [0.1, 0.5, 0.2, 0.5, 0.1], [0.035, 0.07, 0.07, 0.07, 0.07, 0.035], 300, 800, 500, 12, 5, (18, 1.5, 12), seed=3)
make("D_2L_long_strip", 200, 30, ["F.Cu", "B.Cu"], [1.53], [0.07]*2, 200, 300, 200, 8, 4, (14, 2.0, 60), seed=4)
print("ok")
