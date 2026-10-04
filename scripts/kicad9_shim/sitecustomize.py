# KiCad 9 compatibility shim for kicad_fdtd's geometry exporter (commit 74e21eb), loaded only when this folder is on
# PYTHONPATH: kicad_fdtd calls BOARD.GetBoardPolygonOutlines(outline, True), but KiCad 9's signature is
# (SHAPE_POLY_SET&, OUTLINE_ERROR_HANDLER*=None, bool=False, bool=False). This drops a bool second argument.
# kicad_fdtd's own code is unchanged.
try:
    import pcbnew
    _orig = pcbnew.BOARD.GetBoardPolygonOutlines

    def _compat(self, outlines, *args):
        if args and isinstance(args[0], bool):
            args = args[1:]
        return _orig(self, outlines, *args)

    pcbnew.BOARD.GetBoardPolygonOutlines = _compat
except Exception:
    pass
