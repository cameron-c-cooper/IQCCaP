using System;
using System.Collections.Generic;
using System.Linq.Expressions;
using System.Numerics;
using PicoGK;

namespace IcosahedralQC
{
    public static class IcosahedralQuasicrystal
    {
        // ── Golden ratio and normalisation ────────────────────────────────────
        private static readonly double Tau = (1.0 + Math.Sqrt(5.0)) / 2.0;
        private static readonly double Norm = 1.0 / Math.Sqrt(1.0 + Tau * Tau);

        // ── Paper SI Eq.(1) matrix (un-normalised) ────────────────────────────
        // rows 0-2: physical space;   rows 3-5: perpendicular space
        private static readonly double[,] M = {
            { Tau,  Tau,    0,   -1,    0,    1 },  // row 0  phys-x
            {   0,    0,    1,  Tau,    1,  Tau },  // row 1  phys-y
            {   1,   -1, -Tau,    0,  Tau,    0 },  // row 2  phys-z
            { Tau, -Tau,    1,    0,   -1,    0 },  // row 3  perp-x
            {  -1,   -1,    0, -Tau,    0,  Tau },  // row 4  perp-y
            {   0,    0,  Tau,   -1,  Tau,   -1 },  // row 5  perp-z
        };

        // ── RT window: precomputed face normals and half-widths ───────────────
        private static readonly double[,] FaceN = new double[15, 3];
        private static readonly double[] FaceH = new double[15];

        // ── Normalised perp column vectors (v_c = col c of M_perp * Norm) ────
        private static readonly double[,] Vperp = new double[6, 3];

        // ── 20 combinations of 3 directions from {0..5} for tile enumeration ─
        private static readonly (int di, int dj, int dk)[] Dir3;

        // ── Sphere circumradius^2 (documentation only) ────────────────────────
        private static readonly double Rw2Sphere;

        private static readonly int[] CellSubsets = BuildCellSubsets();

        // ── Auto-shell safety margin (extra integer steps beyond corner reach)─
        // 2 extra shells guarantee the RT-window-clipped node cloud comfortably
        // covers the domain corners even after perp-space rejection.
        private const int ShellMargin = 2;

        // ── Hard cap on auto NShell to keep the 6D sweep tractable ────────────
        // (2*NShell+1)^6 iterations.  16 -> ~1.4e9, the practical ceiling.
        private const int MaxAutoShell = 16;

        static IcosahedralQuasicrystal()
        {
            // ── Perp column vectors ───────────────────────────────────────────
            for (int c = 0; c < 6; c++)
            {
                Vperp[c, 0] = M[3, c] * Norm;
                Vperp[c, 1] = M[4, c] * Norm;
                Vperp[c, 2] = M[5, c] * Norm;
            }

            // ── Sphere circumradius^2 (for documentation) ──────────────────────
            double s = 0;
            for (int r = 3; r <= 5; r++) s += M[r, 0] * M[r, 0];
            Rw2Sphere = s * Norm * Norm;   // = 1.0

            // ── RT face normals and half-widths ───────────────────────────────
            int fi = 0;
            for (int i = 0; i < 6; i++)
                for (int j = i + 1; j < 6; j++)
                {
                    double nx = Vperp[i, 1] * Vperp[j, 2] - Vperp[i, 2] * Vperp[j, 1];
                    double ny = Vperp[i, 2] * Vperp[j, 0] - Vperp[i, 0] * Vperp[j, 2];
                    double nz = Vperp[i, 0] * Vperp[j, 1] - Vperp[i, 1] * Vperp[j, 0];
                    FaceN[fi, 0] = nx;
                    FaceN[fi, 1] = ny;
                    FaceN[fi, 2] = nz;

                    double h = 0;
                    for (int k = 0; k < 6; k++)
                        h += Math.Abs(nx * Vperp[k, 0] + ny * Vperp[k, 1] + nz * Vperp[k, 2]);
                    FaceH[fi] = h * 0.5;
                    fi++;
                }

            // ── Direction triples ─────────────────────────────────────────────
            var list = new List<(int, int, int)>();
            for (int a = 0; a < 6; a++)
                for (int b = a + 1; b < 6; b++)
                    for (int c = b + 1; c < 6; c++)
                        list.Add((a, b, c));
            Dir3 = list.ToArray();   // 20 elements
        }

        /// Minimum 6D search half-width needed for accepted lattice points to
        /// reach the corners of the physical domain cube.
        ///
        /// Reach to a corner = (domainSize/2)*sqrt(3). Each single-axis 6D step
        /// projects to `edgeLength` mm, so the body-diagonal step count is
        /// reach/edgeLength. Round up and add ShellMargin, then cap.
        ///
        /// If `requested > 0`, that value is honoured but a warning is logged
        /// when it is below the computed minimum.
        public static int ResolveShell(float edgeLength, float domainSize, int requested)
        {
            double cornerReach = (domainSize / 2.0) * Math.Sqrt(3.0);
            int needed = (int)Math.Ceiling(cornerReach / edgeLength) + ShellMargin;
            int auto = Math.Min(needed, MaxAutoShell);

            if (requested > 0)
            {
                if (requested < needed)
                    Library.Log(
                        $"      [WARN] NShell={requested} under-fills the domain " +
                        $"(need >= {needed} for a={edgeLength}mm, domain={domainSize}mm). " +
                        $"Corners will be empty; auto beam radius will be inflated. " +
                        $"Pass NShell<=0 to auto-size.");
                return requested;
            }

            if (needed > MaxAutoShell)
                Library.Log(
                    $"      [WARN] Auto NShell capped at {MaxAutoShell} " +
                    $"(domain corners want {needed}). Increase MaxAutoShell or " +
                    $"the domain corners may be slightly under-filled.");

            Library.Log($"      Auto NShell : {auto}  " +
                        $"(corner reach {cornerReach:F1}mm / a {edgeLength}mm + margin {ShellMargin})");
            return auto;
        }

        /// Generate the IQC primal truss.
        public static void GeneratePrimalTruss(
            float edgeLength, float domainSize, int nShell, float windowScale,
            out List<Vector3> nodes, out List<(int A, int B)> edges)
        {
            int shell = ResolveShell(edgeLength, domainSize, nShell);
            BuildPrimal(edgeLength, domainSize, shell, windowScale,
                        out nodes, out edges, out _, out _, out _);
        }

        // WARNING: This fails. Dont use it. Please. For me
        public static void GenerateDualTruss(
            float edgeLength, float domainSize, int nShell, float windowScale,
            out List<Vector3> dualNodes, out List<(int A, int B)> dualEdges)
        {
            int shell = ResolveShell(edgeLength, domainSize, nShell);

            // ── Step 1: build primal + extended lookup ────────────────────────
            BuildPrimal(edgeLength, domainSize, shell, windowScale,
                        out List<Vector3> _primalNodes,
                        out List<(int, int)> _primalEdges,
                        out List<int[]> _coords6D,
                        out Dictionary<string, int> _keyToIndex,
                        out Dictionary<string, int> extIndex);

            Library.Log($"      Extended node set: {extIndex.Count}");

            double scale = GetScale(edgeLength);
            float half = domainSize / 2f;

            // Shrink the dual-node acceptance boundary inward by half an edge length.
            float dualKeep = half + edgeLength;

            Vector3[] G = BuildPhysGenerators(scale);

            // ── Step 2: build coord array for iteration over extIndex ─────────
            var allExtCoords = new int[extIndex.Count][];
            foreach (var kv in extIndex)
            {
                var p = kv.Key.Split(',');
                var c = new int[6];
                for (int i = 0; i < 6; i++) c[i] = int.Parse(p[i]);
                allExtCoords[kv.Value] = c;
            }

            // ── Step 3: enumerate all cells ─────────────────────────────────
            var centroids = new List<Vector3>();
            var cellFaceKeys = new List<List<string>>();
            var seen = new HashSet<string>();
            int n3 = 0, n4 = 0, n5 = 0;

            foreach (int[] v in allExtCoords)
            {
                foreach (int S in CellSubsets)
                {
                    if (!TryBuildCell(v, S, G, extIndex, scale,
                        out Vector3 cen, out var fkeys, out string key)) continue;
                    if (!seen.Add(key)) continue;
                    if (Math.Abs(cen.X) > dualKeep || Math.Abs(cen.Y) > dualKeep || Math.Abs(cen.Z) > dualKeep)
                        continue;
                    centroids.Add(cen);
                    cellFaceKeys.Add(fkeys);
                    switch (BitOperations.PopCount((uint)S))
                    {
                        case 3:
                            n3++;
                            break;
                        case 4:
                            n4++;
                            break;
                        default:
                            n5++;
                            break;
                    }

                }
            }
            Library.Log($"      Cells found: {centroids.Count} (rhombohedra {n3}, dodeca {n4}, icosa {n5}");
            var faceOwners = new Dictionary<string, List<int>>(centroids.Count * 6);
            for (int ci = 0; ci < centroids.Count; ci++)
            {
                foreach (string fk in cellFaceKeys[ci])
                {
                    if (!faceOwners.TryGetValue(fk, out var owners))
                        faceOwners[fk] = owners = new List<int>(2);
                    owners.Add(ci);
                }
            }
            var edgesList = new List<(int, int)>();
            foreach (var owners in faceOwners.Values)
            {
                if (owners.Count == 2) edgesList.Add((owners[0], owners[1]));
                else if (owners.Count > 2)
                    Library.Log($"      [WARN] face shared by {owners.Count} cells");
            }

            // ── Step 6: keep only the largest connected component ─────────────
            (centroids, edgesList) = LargestConnectedComponent(centroids, edgesList);

            dualNodes = centroids;
            dualEdges = edgesList;

            if (dualNodes.Count > 0)
                LogFillExtent("Dual", dualNodes, half);
            if (dualEdges.Count > 0)
                LogEdgeLengthRange("Dual", dualNodes, dualEdges);
            Library.Log($"      Dual nodes   : {dualNodes.Count}");
            Library.Log($"      Dual edges   : {dualEdges.Count}");
        }

        // ─────────────────────────────────────────────────────────────────────
        // CORE PRIMAL BUILDER
        // ─────────────────────────────────────────────────────────────────────

        private static void BuildPrimal(
            float edgeLength,
            float domainSize,
            int nShell,
            float windowScale,
            out List<Vector3> nodes,
            out List<(int A, int B)> edges,
            out List<int[]> coords6D,
            out Dictionary<string, int> keyToIndex,
            out Dictionary<string, int> extIndex)
        {
            double scale = GetScale(edgeLength);
            double half = domainSize / 2.0;
            int extShell = nShell + 2;  // ensures rhombohedron corners are present

            // ── Extended index: window only, no physical domain clip ──────────
            extIndex = new Dictionary<string, int>(capacity: 20000);
            for (int n0 = -extShell; n0 <= extShell; n0++)
                for (int n1 = -extShell; n1 <= extShell; n1++)
                    for (int n2 = -extShell; n2 <= extShell; n2++)
                        for (int n3 = -extShell; n3 <= extShell; n3++)
                            for (int n4 = -extShell; n4 <= extShell; n4++)
                                for (int n5 = -extShell; n5 <= extShell; n5++)
                                {
                                    if (!InRT(n0, n1, n2, n3, n4, n5, windowScale)) continue;
                                    extIndex[$"{n0},{n1},{n2},{n3},{n4},{n5}"] = extIndex.Count;
                                }

            // ── Domain primal: window + physical domain clip ──────────────────
            nodes = new List<Vector3>(extIndex.Count / 2);
            coords6D = new List<int[]>(extIndex.Count / 2);
            keyToIndex = new Dictionary<string, int>(extIndex.Count / 2);

            for (int n0 = -nShell; n0 <= nShell; n0++)
                for (int n1 = -nShell; n1 <= nShell; n1++)
                    for (int n2 = -nShell; n2 <= nShell; n2++)
                        for (int n3 = -nShell; n3 <= nShell; n3++)
                            for (int n4 = -nShell; n4 <= nShell; n4++)
                                for (int n5 = -nShell; n5 <= nShell; n5++)
                                {
                                    if (!InRT(n0, n1, n2, n3, n4, n5, windowScale)) continue;

                                    double rx = PhysCoord(0, n0, n1, n2, n3, n4, n5) * scale;
                                    double ry = PhysCoord(1, n0, n1, n2, n3, n4, n5) * scale;
                                    double rz = PhysCoord(2, n0, n1, n2, n3, n4, n5) * scale;
                                    if (Math.Abs(rx) > half || Math.Abs(ry) > half || Math.Abs(rz) > half) continue;

                                    string key = $"{n0},{n1},{n2},{n3},{n4},{n5}";
                                    keyToIndex[key] = nodes.Count;
                                    nodes.Add(new Vector3((float)rx, (float)ry, (float)rz));
                                    coords6D.Add(new[] { n0, n1, n2, n3, n4, n5 });
                                }

            int interiorCount = nodes.Count;
            Library.Log($"      Primal nodes: {interiorCount}");

            // Fill diagnostic on INTERIOR nodes only (before stubs are appended;
            // stubs lie outside the cube by design and would falsely trip the
            // under-fill warning).
            if (interiorCount > 0)
                LogFillExtent("Primal", nodes.GetRange(0, interiorCount), (float)half);

            // Stub nodes are appended AFTER all interior nodes, so indices
            // [0, interiorCount) are interior and [interiorCount, nodes.Count) are
            // stubs. Stub nodes are intentionally NOT added to keyToIndex, so they
            // never spawn their own edges or get treated as lattice vertices, they
            // exist only to anchor the outward beam, and the voxel clip removes the
            // portion beyond the face.
            edges = new List<(int, int)>();
            var stubIndex = new Dictionary<string, int>(); // outside-node key -> node list index

            for (int i = 0; i < interiorCount; i++)
            {
                int[] ci = coords6D[i];
                for (int d = 0; d < 6; d++)
                    for (int sgn = -1; sgn <= 1; sgn += 2)   // both +1 and -1 neighbours
                    {
                        int[] cj = Sh(ci, d, sgn);
                        string kj = Kstr(cj);

                        if (keyToIndex.TryGetValue(kj, out int j))
                        {
                            // both endpoints interior — add once to avoid duplicates.
                            // The +1 scan from the lower-index node already covers this
                            // pair, so only emit when j > i.
                            if (j > i) edges.Add((i, j));
                            continue;
                        }

                        if (!extIndex.ContainsKey(kj)) continue;

                        if (!stubIndex.TryGetValue(kj, out int sj))
                        {
                            // materialise the outside neighbour as a stub node
                            double sx = PhysCoord(0, cj[0], cj[1], cj[2], cj[3], cj[4], cj[5]) * scale;
                            double sy = PhysCoord(1, cj[0], cj[1], cj[2], cj[3], cj[4], cj[5]) * scale;
                            double sz = PhysCoord(2, cj[0], cj[1], cj[2], cj[3], cj[4], cj[5]) * scale;
                            sj = nodes.Count;
                            nodes.Add(new Vector3((float)sx, (float)sy, (float)sz));
                            coords6D.Add((int[])cj.Clone());
                            stubIndex[kj] = sj;
                        }
                        edges.Add((i, sj)); // interior -> stub (sliced by the clip box)
                    }
            }

            int stubCount = nodes.Count - interiorCount;
            Library.Log($"      Primal edges: {edges.Count}  " +
                        $"(incl. {stubCount} boundary stub nodes for Fig-1f-style faces)");

            // Edge-length range (interior + stub edges are all one lattice step = a)
            if (edges.Count > 0)
                LogEdgeLengthRange("Primal", nodes, edges);
        }

        // ─────────────────────────────────────────────────────────────────────
        // WINDOW TEST — RHOMBIC TRIACONTAHEDRON
        // ─────────────────────────────────────────────────────────────────────

        private static bool InRT(
            int n0, int n1, int n2, int n3, int n4, int n5,
            double windowScale)
        {
            double px = (M[3, 0] * n0 + M[3, 1] * n1 + M[3, 2] * n2 + M[3, 3] * n3 + M[3, 4] * n4 + M[3, 5] * n5) * Norm;
            double py = (M[4, 0] * n0 + M[4, 1] * n1 + M[4, 2] * n2 + M[4, 3] * n3 + M[4, 4] * n4 + M[4, 5] * n5) * Norm;
            double pz = (M[5, 0] * n0 + M[5, 1] * n1 + M[5, 2] * n2 + M[5, 3] * n3 + M[5, 4] * n4 + M[5, 5] * n5) * Norm;

            for (int f = 0; f < 15; f++)
            {
                double dot = FaceN[f, 0] * px + FaceN[f, 1] * py + FaceN[f, 2] * pz;
                if (Math.Abs(dot) >= FaceH[f] * windowScale) return false;
            }
            return true;
        }

        // ─────────────────────────────────────────────────────────────────────
        // MATH HELPERS
        // ─────────────────────────────────────────────────────────────────────

        private static double PhysCoord(int axis,
            int n0, int n1, int n2, int n3, int n4, int n5)
            => (M[axis, 0] * n0 + M[axis, 1] * n1 + M[axis, 2] * n2
               + M[axis, 3] * n3 + M[axis, 4] * n4 + M[axis, 5] * n5) * Norm;

        private static double GetScale(double edgeLength)
        {
            double dx = M[0, 0] * Norm, dy = M[1, 0] * Norm, dz = M[2, 0] * Norm;
            return edgeLength / Math.Sqrt(dx * dx + dy * dy + dz * dz);
        }

        private static Vector3 PhysPos(int[] n, double scale) => new Vector3(
            (float)(PhysCoord(0, n[0], n[1], n[2], n[3], n[4], n[5]) * scale),
            (float)(PhysCoord(1, n[0], n[1], n[2], n[3], n[4], n[5]) * scale),
            (float)(PhysCoord(2, n[0], n[1], n[2], n[3], n[4], n[5]) * scale));

        private static int[] Sh(int[] v, int axis, int delta)
        { int[] w = (int[])v.Clone(); w[axis] += delta; return w; }

        private static string Kstr(int[] v) =>
            $"{v[0]},{v[1]},{v[2]},{v[3]},{v[4]},{v[5]}";

        /// Logs the bounding-box half-extent of the accepted nodes vs the domain
        /// half-extent.  A fill fraction well below 1.0 means the window did not
        /// reach the domain boundary (under-fill) — the usual cause of an
        /// inflated auto beam radius.
        private static void LogFillExtent(string tag, List<Vector3> nodes, float half)
        {
            float mx = 0, my = 0, mz = 0;
            foreach (var n in nodes)
            {
                if (Math.Abs(n.X) > mx) mx = Math.Abs(n.X);
                if (Math.Abs(n.Y) > my) my = Math.Abs(n.Y);
                if (Math.Abs(n.Z) > mz) mz = Math.Abs(n.Z);
            }
            float worst = Math.Min(mx, Math.Min(my, mz)); // tightest axis governs fill
            float frac = half > 0 ? worst / half : 0;
            Library.Log($"      {tag} fill    : extent (±{mx:F1}, ±{my:F1}, ±{mz:F1}) mm " +
                        $"vs domain ±{half:F1} mm  (min-axis fill {frac:P0})");
            if (frac < 0.95f)
                Library.Log($"      [WARN] {tag} under-fills domain (min-axis fill {frac:P0} < 95%). " +
                            $"Increase NShell or use auto-sizing (NShell<=0).");
        }

        /// Logs the min/max bonded-edge length.  For the primal truss both should
        /// equal the tile edge `a`; for the dual truss they vary (beams of
        /// varying lengths, per the paper).
        private static void LogEdgeLengthRange(
            string tag, List<Vector3> nodes, List<(int A, int B)> edges)
        {
            float min = float.MaxValue, max = 0;
            foreach (var (a, b) in edges)
            {
                float d = Vector3.Distance(nodes[a], nodes[b]);
                if (d < min) min = d;
                if (d > max) max = d;
            }
            Library.Log($"      {tag} edge len: {min:F3} – {max:F3} mm");
        }

        // ─────────────────────────────────────────────────────────────────────
        // GRAPH HELPERS
        // ─────────────────────────────────────────────────────────────────────

        private static (List<Vector3> nodes, List<(int A, int B)> edges)
            LargestConnectedComponent(
                List<Vector3> nodes,
                List<(int A, int B)> edges)
        {
            int N = nodes.Count;
            if (N == 0) return (nodes, edges);

            var adj = new List<int>[N];
            for (int i = 0; i < N; i++) adj[i] = new List<int>();
            foreach (var (a, b) in edges)
            { adj[a].Add(b); adj[b].Add(a); }

            var compId = new int[N];
            var compSz = new List<int>();
            Array.Fill(compId, -1);

            for (int start = 0; start < N; start++)
            {
                if (compId[start] >= 0) continue;
                int cid = compSz.Count;
                compSz.Add(0);
                var queue = new Queue<int>();
                queue.Enqueue(start);
                compId[start] = cid;
                while (queue.Count > 0)
                {
                    int cur = queue.Dequeue();
                    compSz[cid]++;
                    foreach (int nb in adj[cur])
                        if (compId[nb] < 0)
                        { compId[nb] = cid; queue.Enqueue(nb); }
                }
            }

            int best = 0;
            for (int c = 1; c < compSz.Count; c++)
                if (compSz[c] > compSz[best]) best = c;

            int removed = N - compSz[best];
            if (compSz.Count > 1)
                Library.Log($"      Connected components: {compSz.Count}  " +
                            $"(largest: {compSz[best]} nodes, removed: {removed} from " +
                            $"{compSz.Count - 1} smaller component(s))");

            var remap = new int[N];
            var newNodes = new List<Vector3>(compSz[best]);
            Array.Fill(remap, -1);

            for (int i = 0; i < N; i++)
                if (compId[i] == best)
                { remap[i] = newNodes.Count; newNodes.Add(nodes[i]); }

            var newEdges = new List<(int, int)>();
            foreach (var (a, b) in edges)
                if (remap[a] >= 0 && remap[b] >= 0)
                    newEdges.Add((remap[a], remap[b]));

            return (newNodes, newEdges);
        }

        // NOTE: this is the minimum VERTEX SPACING among accepted nodes, which in
        // an IQC is the short diagonal of the oblate rhombohedron (~0.56*a), NOT
        // the tile edge length.  Do not use it to validate the scale; use the
        // bonded-edge length range (LogEdgeLengthRange) instead.
        private static float NNDist(List<Vector3> pts)
        {
            float m = float.MaxValue;
            int lim = Math.Min(pts.Count, 600);
            for (int i = 0; i < lim; i++)
                for (int j = i + 1; j < lim; j++)
                { float d = Vector3.DistanceSquared(pts[i], pts[j]); if (d < m) m = d; }
            return (float)Math.Sqrt(m);
        }
        private static int[] BuildCellSubsets()
        {
            var list = new List<int>();
            for (int mask = 1; mask < 64; mask++)
            {
                int pc = BitOperations.PopCount((uint)mask);
                if (pc >= 3 && pc <= 5) list.Add(mask);
            }
            return list.ToArray();
        }
        private static Vector3[] BuildPhysGenerators(double scale)
        {
            var g = new Vector3[6];
            for (int c = 0; c < 6; c++)
            {
                g[c] = new Vector3(
                    (float)(M[0, c] * Norm * scale),
                    (float)(M[1, c] * Norm * scale),
                    (float)(M[2, c] * Norm * scale));
            }
            return g;
        }

        private static int[] AddMask(int[] v, int mask)
        {
            int[] w = (int[])v.Clone();
            int m = mask;
            while (m != 0)
            {
                int b = BitOperations.TrailingZeroCount(m);
                m &= m - 1;
                w[b] += 1;
            }
            return w;
        }
        private static bool TryBuildCell(
            int[] v, int S, Vector3[] G, Dictionary<string, int> extIndex,
            double scale, out Vector3 centroid, out List<string> faceKeys,
            out string dedupKey)
        {
            centroid = default;
            faceKeys = null;
            dedupKey = null;
            if (!extIndex.ContainsKey(Kstr(AddMask(v, S)))) return false;
            int k = BitOperations.PopCount((uint)S);
            int expectedVerts = k * k - k + 2;
            var faces = new List<int[]>(k * (k - 1));
            var vertexMasks = new HashSet<int>();
            int[] sBits = BitsOf(S);
            for (int a = 0; a < sBits.Length; a++)
            {
                for (int b = a + 1; b < sBits.Length; b++)
                {
                    int i = sBits[a], j = sBits[b];
                    int rest = S & ~((1 << i) | (1 << j));
                    Vector3 n = Vector3.Cross(G[i], G[j]);
                    for (int side = -1; side <= 1; side += 2)
                    {
                        int baseMask = 0, r = rest;
                        while (r != 0)
                        {
                            int m = BitOperations.TrailingZeroCount(r);
                            r &= r - 1;
                            if ((Vector3.Dot(n, G[m]) > -0f ? 1 : -1) == side)
                                baseMask |= (1 << m);
                        }
                        int c0 = baseMask,
                            c1 = baseMask | (1 << i),
                            c2 = baseMask | (1 << j),
                            c3 = baseMask | (1 << i) | (1 << j);
                        faces.Add(new[] { c0, c1, c2, c3 });
                        vertexMasks.Add(c0);
                        vertexMasks.Add(c1);
                        vertexMasks.Add(c2);
                        vertexMasks.Add(c3);
                    }
                }
            }
            // degen geometry guard
            if (vertexMasks.Count != expectedVerts) return false;
            // validate every vertex is present, interior subset sum absent
            var maskToId = new Dictionary<int, int>(expectedVerts);
            for (int t = S; ; t = (t - 1) & S)
            {
                bool present = extIndex.TryGetValue(Kstr(AddMask(v, t)), out int id);
                if (vertexMasks.Contains(t))
                {
                    if (!present) return false; // missing hull vertex
                    maskToId[t] = id;
                }
                else if (present) return false; // interior point exists -> subdivided
                if (t == 0) break;
            }
            // centroid = v + 1/2 sum(generators)
            Vector3 cen = PhysPos(v, scale);

            foreach (int i in sBits) cen += 0.5f * G[i];

            centroid = cen;
            faceKeys = new List<String>(faces.Count);
            foreach (var f in faces)
            {
                int[] ids = {
                    maskToId[f[0]],
                    maskToId[f[1]],
                    maskToId[f[2]],
                    maskToId[f[3]],
                };
                Array.Sort(ids);
                faceKeys.Add($"{ids[0]}_{ids[1]}_{ids[2]}_{ids[3]}");
            }
            var vids = new List<int>(maskToId.Values);
            vids.Sort();
            dedupKey = string.Join("|", vids);
            return true;
        }
        private static int[] BitsOf(int mask)
        {
            var list = new List<int>(6);
            while (mask != 0)
            {
                int b = BitOperations.TrailingZeroCount(mask);
                list.Add(b);
                mask &= mask - 1;
            }
            return list.ToArray();
        }
    }
}
