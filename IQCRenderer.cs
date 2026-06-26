using System;
using System.IO;
using System.Numerics;
using System.Collections.Generic;
using PicoGK;

namespace IcosahedralQC
{
    public enum TrussMode { Primal, Dual }

    public class IQCParameters
    {
        public float EdgeLength = 3.9f;

        public int NShell = 0;

        public float WindowScale = 1.0f;

        public TrussMode Mode = TrussMode.Primal;

        public float RelativeDensity = 0.10f;

        public float BeamRadius = 0.41f;

        public float MinBeamRadius = 0.25f;

        public float DomainSize = 57.0f;

        public float SmoothDistMM = 0.15f;

        public int SmoothPasses = 1;

        public string StlFileName;
    }

    public static class IQCRenderer
    {
        public static IQCParameters Params = new IQCParameters();

        public static void Task() => RunTask(Params);

        private static void RunTask(IQCParameters p)
        {
            Library.Log("=== IQC Truss Generator (v6 — auto-shell) ===");
            Library.Log($"  Mode           : {p.Mode}");
            Library.Log($"  Edge length    : {p.EdgeLength} mm");
            Library.Log($"  Domain size    : {p.DomainSize} mm  ({p.DomainSize / p.EdgeLength:F1}a)");
            Library.Log($"  NShell         : {(p.NShell <= 0 ? "auto" : p.NShell.ToString())}");
            Library.Log($"  WindowScale    : {p.WindowScale}");
            Library.Log($"  Rel. density   : {p.RelativeDensity:P0}");

            // ── 1. Generate truss ─────────────────────────────────────────────
            Library.Log("[1/5] Generating truss...");

            List<Vector3> nodes;
            List<(int A, int B)> edges;

            if (p.Mode == TrussMode.Dual)
            {
                IcosahedralQuasicrystal.GenerateDualTruss(
                    edgeLength: p.EdgeLength,
                    domainSize: p.DomainSize,
                    nShell: p.NShell,
                    windowScale: p.WindowScale,
                    dualNodes: out nodes,
                    dualEdges: out edges);
            }
            else if (p.Mode == TrussMode.Primal)
            {
                IcosahedralQuasicrystal.GeneratePrimalTruss(
                    edgeLength: p.EdgeLength,
                    domainSize: p.DomainSize,
                    nShell: p.NShell,
                    windowScale: p.WindowScale,
                    nodes: out nodes,
                    edges: out edges);
            }
            else
            {
                throw new Exception("Unknown TrussMode.");
            }

            Library.Log($"      Nodes: {nodes.Count},  Edges: {edges.Count}");

            if (nodes.Count == 0)
                throw new Exception(
                    "No nodes generated. Check WindowScale; auto-shell should fill the domain.");
            if (edges.Count == 0)
                throw new Exception(
                    "No edges generated. The primal or dual graph is empty.");

            // ── 2. Compute beam radius ────────────────────────────────────────
            Library.Log("[2/5] Computing beam radius...");

            float beamRadius = ComputeBeamRadius(nodes, edges, p);

            // ── 3. Build PicoGK Lattice ───────────────────────────────────────
            Library.Log("[3/5] Building Lattice beams...");
            var lat = new Lattice();

            foreach (var (a, b) in edges)
            {
                if (Vector3.DistanceSquared(nodes[a], nodes[b]) < 1e-6f) continue;
                lat.AddBeam(nodes[a], beamRadius,
                            nodes[b], beamRadius, false);
            }

            float nodeR = beamRadius * 1.15f;
            float half = p.DomainSize / 2f;
            // Node-joint spheres only at REAL lattice vertices inside the domain.
            // Boundary stub nodes sit outside the cube (they only anchor the
            // outward beam, which the clip slices); adding a sphere there would
            // leave a spurious bump just outside a face. The clip would trim most
            // of it, but skipping them is cleaner and faster.
            foreach (var nd in nodes)
                if (Math.Abs(nd.X) <= half && Math.Abs(nd.Y) <= half && Math.Abs(nd.Z) <= half)
                    lat.AddSphere(nd, nodeR);

            // foreach (var (a, b) in edges)
            // {
            //     if (Vector3.DistanceSquared(nodes[a], nodes[b]) < 1e-6f) continue;
            //     lat.AddBeam(nodes[a], beamRadius, nodes[b], beamRadius, false);
            // }

            // ── 4. Voxelise + domain clip ─────────────────────────────────────
            Library.Log("[4/5] Voxelising...");
            var voxTruss = new Voxels(lat);

            // mshCreateCube(vecSize, vecCenter): vecSize is the FULL edge length.
            // For a domain cube of side DomainSize (half-extent DomainSize/2), the
            // full edge length is DomainSize — NOT DomainSize*2. The previous *2
            // made a 114mm box (half-extent 57mm) that clipped nothing, because the
            // primal node clip had already capped everything at ±28.5mm. With the
            // correct 57mm box, this BoolIntersect actually SLICES beams that cross
            // the ±28.5mm faces, leaving the partial-beam stubs seen in the paper's
            // fabricated samples (Figure 1f). It pairs with the boundary-beam change
            // in IcosahedralQuasicrystal.BuildPrimal (keepBoundaryBeams).
            float full = p.DomainSize;
            Vector3 fullVec = new Vector3(full, full, full);
            Vector3 origin = Vector3.Zero;
            Mesh mshBox = Utils.mshCreateCube(fullVec, origin);
            var voxBox = new Voxels(mshBox);
            voxTruss.BoolIntersect(voxBox);
            Mesh topPlate = Utils.mshCreateCube(new Vector3(p.DomainSize, p.DomainSize, 2), new Vector3(0, 0, p.DomainSize / 2 + beamRadius));
            Mesh botPlate = Utils.mshCreateCube(new Vector3(p.DomainSize, p.DomainSize, 2), new Vector3(0, 0, -p.DomainSize / 2 - beamRadius));
            voxTruss.BoolAdd(new Voxels(topPlate));
            voxTruss.BoolAdd(new Voxels(botPlate));

            // ── 5. Smooth ─────────────────────────────────────────────────────
            if (p.SmoothPasses > 0 && p.SmoothDistMM > 0f)
            {
                Library.Log($"[4b/5] Smoothing ({p.SmoothPasses}× @ {p.SmoothDistMM} mm)...");
                for (int i = 0; i < p.SmoothPasses; i++)
                    voxTruss.Smoothen(p.SmoothDistMM);
            }
            else
            {
                Library.Log("[4b/5] Smoothing skipped.");
            }

            // ── 6. Export STL ─────────────────────────────────────────────────
            Library.Log("[5/5] Exporting STL...");
            string stlPath = Path.Combine(Library.strLogFolder, p.StlFileName + ".stl");
            new Mesh(voxTruss).SaveToStlFile(stlPath);
            Library.Log($"      Saved: {stlPath}");

            // ── 7. Viewer ─────────────────────────────────────────────────────
            Library.oViewer().SetGroupMaterial(0,
                p.Mode == TrussMode.Dual ? "CC3333FF" : "3355CCFF",
                0.5f, 0.8f);
            Library.oViewer().Add(voxTruss);

            PrintStats(nodes, edges, beamRadius, p);
            Library.Log("=== Done ===");
        }

        // ─────────────────────────────────────────────────────────────────────

        private static float ComputeBeamRadius(
            List<Vector3> nodes,
            List<(int A, int B)> edges,
            IQCParameters p)
        {
            double totalLength = 0;
            foreach (var (a, b) in edges)
                totalLength += Vector3.Distance(nodes[a], nodes[b]);

            float beamRadius = p.BeamRadius;

            if (p.RelativeDensity > 0f)
            {
                // ρ = π r² L_total / V_domain  →  r = √(ρ V / π L)
                // Idealised cylinder sum: ignores joint/sphere overlap, so the
                // result is a lower bound on the CAD-true radius. Expect a few-%
                // offset from paper Table 2 even with correct domain fill.
                double vDomain = Math.Pow(p.DomainSize, 3.0);
                beamRadius = (float)Math.Sqrt(
                    p.RelativeDensity * vDomain / (Math.PI * totalLength));

                if (beamRadius < p.MinBeamRadius)
                {
                    Library.Log($"      Auto-radius {beamRadius:F3} mm < MinBeamRadius " +
                                $"{p.MinBeamRadius:F3} mm — clamped. Actual ρ will exceed target.");
                    beamRadius = p.MinBeamRadius;
                }
            }

            Library.Log($"      Beam radius : {beamRadius:F3} mm  (dia {beamRadius * 2:F3} mm)");
            Library.Log($"      Total strut length : {totalLength / 1000.0:F2} m");
            return beamRadius;
        }

        // ── Filter primal nodes/edges to strict domain box (unused; kept) ─────
        private static void FilterToDomain(
            ref List<Vector3> nodes,
            ref List<(int A, int B)> edges,
            float domainSize)
        {
            float half = domainSize / 2f;
            var keep = new bool[nodes.Count];
            for (int i = 0; i < nodes.Count; i++)
            {
                var n = nodes[i];
                keep[i] = Math.Abs(n.X) <= half &&
                           Math.Abs(n.Y) <= half &&
                           Math.Abs(n.Z) <= half;
            }

            var newNodes = new List<Vector3>();
            var remap = new int[nodes.Count];
            for (int i = 0; i < nodes.Count; i++)
            {
                remap[i] = keep[i] ? newNodes.Count : -1;
                if (keep[i]) newNodes.Add(nodes[i]);
            }

            var newEdges = new List<(int, int)>();
            foreach (var (a, b) in edges)
            {
                if (remap[a] >= 0 && remap[b] >= 0)
                    newEdges.Add((remap[a], remap[b]));
            }

            nodes = newNodes;
            edges = newEdges;
        }

        // ── Stats ─────────────────────────────────────────────────────────────

        private static void PrintStats(
            List<Vector3> nodes,
            List<(int A, int B)> edges,
            float beamRadius,
            IQCParameters p)
        {
            double totalLen = 0;
            float minEdge = float.MaxValue;
            float maxEdge = 0;
            foreach (var (a, b) in edges)
            {
                float d = Vector3.Distance(nodes[a], nodes[b]);
                totalLen += d;
                if (d < minEdge) minEdge = d;
                if (d > maxEdge) maxEdge = d;
            }

            double vStruts = Math.PI * beamRadius * beamRadius * totalLen;
            double vDomain = Math.Pow(p.DomainSize, 3.0);

            Library.Log("─── Geometry summary ───────────────────────────────");
            Library.Log($"  Mode                 : {p.Mode}");
            Library.Log($"  Nodes                : {nodes.Count}");
            Library.Log($"  Edges                : {edges.Count}");
            Library.Log($"  Beam diameter        : {beamRadius * 2:F3} mm");
            Library.Log($"  Edge length range    : {minEdge:F2} – {maxEdge:F2} mm");
            Library.Log($"  Total strut length   : {totalLen / 1000.0:F2} m");
            Library.Log($"  Achieved rel. density: {vStruts / vDomain:P1} (idealised, no joint overlap)");
            Library.Log("────────────────────────────────────────────────────");
        }
    }
}
