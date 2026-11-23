"""
D8 hydrology reproducing the ArcGIS Spatial Analyst tools used by the original
ValleyMorph Tool: FlowDirection (NORMAL), FlowAccumulation (FLOAT),
StreamOrder (STRAHLER) and StreamToFeature (SIMPLIFY).

Direction codes follow ArcGIS: E=1, SE=2, S=4, SW=8, W=16, NW=32, N=64, NE=128.
"""

import math
import warnings
from collections import deque

import numpy as np

# (code, drow, dcol) in ArcGIS code order
DIRS = [(1, 0, 1), (2, 1, 1), (4, 1, 0), (8, 1, -1),
        (16, 0, -1), (32, -1, -1), (64, -1, 0), (128, -1, 1)]
CODE_TO_D = {c: (dr, dc) for c, dr, dc in DIRS}


def _shift(a, dr, dc, fill=np.nan):
    """b[r, c] = a[r + dr, c + dc] (fill outside)."""
    out = np.full(a.shape, fill, dtype=a.dtype)
    nr, nc = a.shape
    rs = slice(max(0, -dr), nr - max(0, dr))
    cs = slice(max(0, -dc), nc - max(0, dc))
    rd = slice(max(0, dr), nr - max(0, -dr))
    cd = slice(max(0, dc), nc - max(0, -dc))
    out[rs, cs] = a[rd, cd]
    return out


def fill_sinks(z_in):
    """
    Fills depressions (priority-flood, Wang & Liu 2006 / Barnes et al. 2014):
    every cell is raised to the lowest level at which it can drain to the edge
    of the data (raster boundary or NoData).  Filled areas become flat and are
    routed by flow_direction's flat handling.
    """
    import heapq
    z = z_in.copy()
    valid = ~np.isnan(z)
    nr, nc = z.shape
    done = ~valid
    heap = []
    border = np.zeros((nr, nc), dtype=bool)
    border[0, :] = border[-1, :] = border[:, 0] = border[:, -1] = True
    edge = valid & (border | ~np.stack([~np.isnan(_shift(z, dr, dc)) for _, dr, dc in DIRS]).all(axis=0))
    for r, c in zip(*np.nonzero(edge)):
        heapq.heappush(heap, (z[r, c], int(r), int(c)))
        done[r, c] = True
    while heap:
        h, r, c = heapq.heappop(heap)
        for _, dr, dc in DIRS:
            rr, cc = r + dr, c + dc
            if 0 <= rr < nr and 0 <= cc < nc and not done[rr, cc]:
                done[rr, cc] = True
                if z[rr, cc] < h:
                    z[rr, cc] = h
                heapq.heappush(heap, (z[rr, cc], rr, cc))
    return z


def flow_direction(z_in):
    """
    ArcGIS FlowDirection, force_flow = NORMAL.

    * Single-cell sinks that are not on an edge are filled to their lowest
      neighbour before directions are assigned.
    * Each cell drains to the neighbour with the steepest drop
      (drop = dz / distance, diagonals at sqrt(2)).  Ties are broken by
      enlarging the neighbourhood one step in each tied direction, then by
      ArcGIS code order.
    * Cells on the raster boundary whose steepest inward drop is not positive
      drain out of the surface.  NoData neighbours are otherwise ignored.
    * Flat cells drain toward the nearest cell of the flat that has a lower
      outlet.  ArcGIS 10.1 does not document its own flat rule, so paths
      across large flats can differ from ArcGIS by a few cells.  Flats without
      a lower outlet drain toward NoData where they touch it, and are sinks
      otherwise.
    * Cells next to NoData that are lower than all valid neighbours drain
      toward the NoData.
    """
    z = z_in.copy()
    valid = ~np.isnan(z)
    nr, nc = z.shape
    dist = [math.sqrt(2.0) if dr and dc else 1.0 for _, dr, dc in DIRS]

    def neighbours(zz):
        return [_shift(zz, dr, dc) for _, dr, dc in DIRS]

    nb = neighbours(z)
    nbvalid = np.stack([~np.isnan(a) for a in nb])
    # Edge = the physical edge of the raster.  NoData neighbours are ignored.
    border = np.zeros((nr, nc), dtype=bool)
    border[0, :] = border[-1, :] = border[:, 0] = border[:, -1] = True
    edge = valid & border

    # Fill single-cell sinks: interior cells with eight valid, higher neighbours.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # cells with no valid neighbour
        nbmin = np.nanmin(np.stack(nb), axis=0)
    sink1 = valid & ~edge & nbvalid.all(axis=0) & (z < nbmin)
    z[sink1] = nbmin[sink1]
    nb = neighbours(z)

    drops = np.stack([(z - nb[k]) / dist[k] for k in range(8)])  # NaN where invalid
    with np.errstate(invalid="ignore"):
        maxdrop = np.nanmax(np.where(np.isnan(drops), -np.inf, drops), axis=0)

    fdir = np.zeros((nr, nc), dtype=np.int32)

    # Case A: positive steepest drop.
    caseA = valid & (maxdrop > 0)
    for r, c in zip(*np.nonzero(caseA)):
        m = maxdrop[r, c]
        cand = [k for k in range(8) if drops[k, r, c] == m]
        if len(cand) > 1:
            cand = _break_tie(z, r, c, cand, dist)
        fdir[r, c] = DIRS[cand[0]][0]

    # Case B: edge cells with no positive inward drop drain off the surface.
    caseB = valid & edge & ~caseA
    for r, c in zip(*np.nonzero(caseB)):
        fdir[r, c] = _middle_direction([not (0 <= r + dr < nr and 0 <= c + dc < nc) for _, dr, dc in DIRS])

    # Case C: interior flats drain to the nearest cell of the flat with a lower outlet.
    flat = valid & ~edge & ~caseA & (maxdrop == 0)
    nodata_nb = valid & ~nbvalid.all(axis=0)
    if flat.any():
        _resolve_flats(z, fdir, flat, valid, nodata_nb)

    # Case D: cells next to NoData that are lower than all their valid
    # neighbours (not filled as sinks) drain toward the NoData.
    pit = valid & ~edge & ~caseA & (maxdrop < 0) & nodata_nb
    for r, c in zip(*np.nonzero(pit)):
        fdir[r, c] = _toward_nodata(r, c, nbvalid)
    return fdir


def _toward_nodata(r, c, nbvalid):
    return _middle_direction([not nbvalid[k, r, c] for k in range(8)])


def _middle_direction(open_dirs):
    """Pick a direction among the open (NoData / off-raster) neighbours: the
    middle of the longest run of adjacent open directions, preferring the
    cardinal direction when the run has an even length."""
    ks = [k for k in range(8) if open_dirs[k]]
    if not ks:
        return 0
    if len(ks) == 8:
        return DIRS[0][0]
    # rotate so the sequence starts right after a closed direction
    start = next(k for k in range(8) if not open_dirs[k])
    best, run = [], []
    for i in range(1, 9):
        k = (start + i) % 8
        if open_dirs[k]:
            run.append(k)
        else:
            if len(run) > len(best):
                best = run
            run = []
    if len(run) > len(best):
        best = run
    n = len(best)
    if n % 2:
        k = best[n // 2]
    else:
        a, b = best[n // 2 - 1], best[n // 2]
        k = a if a % 2 == 0 else b  # even index in DIRS = cardinal
    return DIRS[k][0]


def _break_tie(z, r, c, cand, dist):
    nr, nc = z.shape
    best, bestv = [], -np.inf
    for k in cand:
        _, dr, dc = DIRS[k]
        rr, cc = r + 2 * dr, c + 2 * dc
        if 0 <= rr < nr and 0 <= cc < nc and not np.isnan(z[rr, cc]):
            v = (z[r, c] - z[rr, cc]) / (2 * dist[k])
        else:
            v = -np.inf
        if v > bestv:
            best, bestv = [k], v
        elif v == bestv:
            best.append(k)
    return best if best else cand


def _resolve_flats(z, fdir, flat, valid, nodata_nb):
    nr, nc = z.shape
    seen = np.zeros(flat.shape, dtype=bool)

    def bfs(q):
        while q:
            r, c = q.popleft()
            for code, dr, dc in DIRS:
                rr, cc = r + dr, c + dc
                if 0 <= rr < nr and 0 <= cc < nc and flat[rr, cc] and not seen[rr, cc] \
                        and z[rr, cc] == z[r, c]:
                    fdir[rr, cc] = _code_from(-dr, -dc)  # drains back toward (r, c)
                    seen[rr, cc] = True
                    q.append((rr, cc))

    # 1) Seeds: flat cells next to an equal-height cell that already drains downhill.
    q = deque()
    for r, c in zip(*np.nonzero(flat)):
        for code, dr, dc in DIRS:
            rr, cc = r + dr, c + dc
            if 0 <= rr < nr and 0 <= cc < nc and valid[rr, cc] and not flat[rr, cc] \
                    and fdir[rr, cc] != 0 and z[rr, cc] == z[r, c]:
                fdir[r, c] = code
                seen[r, c] = True
                q.append((r, c))
                break
    bfs(q)
    # 2) Flats with no lower outlet drain off the surface where they touch NoData.
    q = deque()
    nbv = np.stack([~np.isnan(_shift(z, dr, dc)) for _, dr, dc in DIRS])
    for r, c in zip(*np.nonzero(flat & ~seen & nodata_nb)):
        fdir[r, c] = _toward_nodata(r, c, nbv)
        seen[r, c] = True
        q.append((r, c))
    bfs(q)


def _code_from(dr, dc):
    for code, a, b in DIRS:
        if a == dr and b == dc:
            return code
    return 0


def downstream_index(fdir, valid):
    """Flat index of the downstream cell, or -1 when flow leaves the surface/stops."""
    nr, nc = fdir.shape
    down = np.full(nr * nc, -1, dtype=np.int64)
    for code, dr, dc in DIRS:
        rs, cs_ = np.nonzero((fdir == code) & valid)
        rr, cc = rs + dr, cs_ + dc
        ok = (rr >= 0) & (rr < nr) & (cc >= 0) & (cc < nc)
        ok[ok] &= valid[rr[ok], cc[ok]]
        down[rs[ok] * nc + cs_[ok]] = rr[ok] * nc + cc[ok]
    return down


def topological_order(down, valid_flat):
    """Cells ordered so every cell comes after all cells draining into it."""
    n = down.size
    indeg = np.bincount(down[down >= 0], minlength=n)
    order = []
    frontier = np.nonzero(valid_flat & (indeg == 0))[0]
    indeg = indeg.copy()
    while frontier.size:
        order.append(frontier)
        d = down[frontier]
        d = d[d >= 0]
        np.subtract.at(indeg, d, 1)
        frontier = np.unique(d[indeg[d] == 0])
    return np.concatenate(order) if order else np.array([], dtype=np.int64)


def flow_accumulation(fdir, valid):
    """ArcGIS FlowAccumulation (FLOAT): number of upstream cells."""
    nr, nc = fdir.shape
    down = downstream_index(fdir, valid)
    vf = valid.ravel()
    order = topological_order(down, vf)
    acc = np.zeros(nr * nc)
    for i in order:
        d = down[i]
        if d >= 0:
            acc[d] += acc[i] + 1
    out = acc.reshape(nr, nc)
    out[~valid] = np.nan
    return out, down, order


def stream_order_strahler(stream, down, order):
    """ArcGIS StreamOrder (STRAHLER) on the stream cells (bool, flat-indexed)."""
    n = down.size
    maxo = np.zeros(n, dtype=np.int32)
    cnt = np.zeros(n, dtype=np.int32)
    so = np.zeros(n, dtype=np.int32)
    for i in order:
        if not stream[i]:
            continue
        so[i] = 1 if maxo[i] == 0 else maxo[i] + (1 if cnt[i] >= 2 else 0)
        d = down[i]
        if d >= 0 and stream[d]:
            if so[i] > maxo[d]:
                maxo[d], cnt[d] = so[i], 1
            elif so[i] == maxo[d]:
                cnt[d] += 1
    return so


def stream_to_feature(grid, stream2d, fdir, simplify_tol=None):
    """
    ArcGIS StreamToFeature.  Lines run downstream between sources, junctions and
    outlets.  A source line starts half a cell upstream of the source cell centre
    and an outlet line ends half a cell downstream of the outlet cell centre.
    Returns a list of (points[(x, y)], from_cell, to_cell).
    """
    nr, nc = stream2d.shape
    cs = grid.cs
    stream = stream2d.ravel()
    valid = ~np.isnan(grid.data)
    down = downstream_index(fdir, valid)
    sdown = np.where((down >= 0) & stream[np.maximum(down, 0)], down, -1)
    sdown[~stream] = -1
    upcount = np.bincount(sdown[sdown >= 0], minlength=nr * nc)

    def xy(i):
        r, c = divmod(int(i), nc)
        return grid.center(r, c)

    def dvec(i):
        r, c = divmod(int(i), nc)
        code = int(fdir[r, c])
        if code in CODE_TO_D:
            dr, dc = CODE_TO_D[code]
            return dc * cs, -dr * cs
        return 0.0, 0.0

    arcs = []
    starts = np.nonzero(stream & ((upcount == 0) | (upcount >= 2)))[0]
    for s in starts:
        x, y = xy(s)
        pts = []
        if upcount[s] == 0:
            vx, vy = dvec(s)
            pts.append((x - 0.5 * vx, y - 0.5 * vy))
        pts.append((x, y))
        i = s
        while True:
            d = sdown[i]
            if d < 0:
                vx, vy = dvec(i)
                if vx or vy:
                    pts.append((pts[-1][0] + 0.5 * vx, pts[-1][1] + 0.5 * vy))
                break
            pts.append(xy(d))
            i = d
            if upcount[d] >= 2:
                break
        if simplify_tol is not None:
            pts = douglas_peucker(pts, simplify_tol)
        if len(pts) >= 2:
            arcs.append((pts, int(s), int(i)))
    return arcs


def douglas_peucker(pts, tol):
    if len(pts) <= 2:
        return pts
    # remove exact collinear first
    keep = [True] * len(pts)
    stack = [(0, len(pts) - 1)]
    while stack:
        a, b = stack.pop()
        ax, ay = pts[a]
        bx, by = pts[b]
        L = math.hypot(bx - ax, by - ay)
        best, bi = -1.0, None
        for k in range(a + 1, b):
            px, py = pts[k]
            if L == 0:
                dd = math.hypot(px - ax, py - ay)
            else:
                dd = abs((bx - ax) * (ay - py) - (ax - px) * (by - ay)) / L
            if dd > best:
                best, bi = dd, k
        if bi is not None and best > tol:
            stack.append((a, bi))
            stack.append((bi, b))
        else:
            for k in range(a + 1, b):
                keep[k] = False
    return [p for p, k in zip(pts, keep) if k]
