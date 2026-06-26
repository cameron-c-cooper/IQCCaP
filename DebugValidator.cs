// // =============================================================================
// //  DebugValidator.cs  (optional standalone test — no PicoGK required)
// //
// //  Run this first to validate the quasicrystal geometry before committing
// //  to a full voxelisation run. It checks:
// //    • Node and edge counts are in expected ranges
// //    • Edge lengths cluster around EdgeLength
// //    • The dual graph is connected
// //    • Basic icosahedral symmetry (5-fold about the z-axis)
// //
// //  Usage: add this file to the project and call DebugValidator.Run()
// //         from Program.cs instead of IQCRenderer.Task(), or run standalone.
// // =============================================================================

// using System;
// using System.Collections.Generic;
// using System.Numerics;
// using IcosahedralQC;

// namespace IcosahedralQC
// {
//     public static class DebugValidator
//     {
//         public static void Run()
//         {
//             Console.WriteLine("=== IQC Geometry Validator ===\n");

//             float edgeLength = 6.0f;
//             float domain     = 57.0f;
//             int   nShell     = 4;

//             // ── Test 1: Primal truss ──────────────────────────────────────────
//             Console.WriteLine("── Primal IQC truss ────────────────────────────────");
//             IcosahedralQuasicrystal.GenerateTruss(
//                 edgeLength, domain, nShell, 1.0f,
//                 out var primalNodes, out var primalEdges);

//             Console.WriteLine($"  Nodes: {primalNodes.Count}");
//             Console.WriteLine($"  Edges: {primalEdges.Count}");

//             AnalyseEdgeLengths(primalNodes, primalEdges, edgeLength, "Primal");

//             // ── Test 2: Dual truss ────────────────────────────────────────────
//             Console.WriteLine("\n── Dual IQC truss ──────────────────────────────────");
//             IcosahedralQuasicrystal.GenerateDualTruss(
//                 edgeLength, domain, nShell, 1.0f,
//                 out var dualNodes, out var dualEdges);

//             Console.WriteLine($"  Nodes: {dualNodes.Count}");
//             Console.WriteLine($"  Edges: {dualEdges.Count}");

//             AnalyseEdgeLengths(dualNodes, dualEdges, edgeLength, "Dual");

//             // ── Test 3: Connectivity check ────────────────────────────────────
//             Console.WriteLine("\n── Connectivity (BFS from node 0) ─────────────────");
//             var adj = new List<int>[dualNodes.Count];
//             for (int i = 0; i < dualNodes.Count; i++) adj[i] = new List<int>();
//             foreach (var (a, b) in dualEdges)
//             {
//                 adj[a].Add(b);
//                 adj[b].Add(a);
//             }
//             int reached = BFS(0, adj, dualNodes.Count);
//             bool connected = reached == dualNodes.Count;
//             Console.WriteLine($"  Reachable: {reached} / {dualNodes.Count}  →  {(connected ? "CONNECTED ✓" : "DISCONNECTED ✗")}");

//             // ── Test 4: Bounding box ──────────────────────────────────────────
//             Console.WriteLine("\n── Bounding box ────────────────────────────────────");
//             BBox(dualNodes, out var bMin, out var bMax);
//             Console.WriteLine($"  min: ({bMin.X:F1}, {bMin.Y:F1}, {bMin.Z:F1})");
//             Console.WriteLine($"  max: ({bMax.X:F1}, {bMax.Y:F1}, {bMax.Z:F1})");

//             // ── Test 5: Approx 5-fold symmetry count ──────────────────────────
//             Console.WriteLine("\n── 5-fold symmetry test (C5 about Z-axis) ──────────");
//             int symCount = CountApproxC5Nodes(dualNodes, edgeLength * 0.3f);
//             Console.WriteLine($"  Nodes with ~C5-symmetric partner: {symCount}");

//             Console.WriteLine("\n=== Validation complete ===");
//         }

//         // ── Helpers ──────────────────────────────────────────────────────────

//         private static void AnalyseEdgeLengths(
//             List<Vector3> nodes, List<(int A, int B)> edges,
//             float nominalLength, string label)
//         {
//             if (edges.Count == 0) { Console.WriteLine("  No edges!"); return; }

//             float min = float.MaxValue, max = 0, sum = 0;
//             var hist = new int[10]; // bins: [0.5a, 1.5a] in 0.1a steps
//             foreach (var (a, b) in edges)
//             {
//                 float d = Vector3.Distance(nodes[a], nodes[b]);
//                 if (d < min) min = d;
//                 if (d > max) max = d;
//                 sum += d;
//                 int bin = (int)((d / nominalLength - 0.5) * 10);
//                 if (bin >= 0 && bin < 10) hist[bin]++;
//             }
//             float avg = sum / edges.Count;
//             Console.WriteLine($"  {label} edge lengths — min:{min:F2}  avg:{avg:F2}  max:{max:F2}  (nominal={nominalLength:F2} mm)");
//         }

//         private static int BFS(int start, List<int>[] adj, int n)
//         {
//             var visited = new bool[n];
//             var queue   = new Queue<int>();
//             queue.Enqueue(start);
//             visited[start] = true;
//             int count = 0;
//             while (queue.Count > 0)
//             {
//                 int v = queue.Dequeue();
//                 count++;
//                 foreach (int nb in adj[v])
//                     if (!visited[nb]) { visited[nb] = true; queue.Enqueue(nb); }
//             }
//             return count;
//         }

//         private static void BBox(List<Vector3> nodes, out Vector3 bMin, out Vector3 bMax)
//         {
//             bMin = new Vector3(float.MaxValue);
//             bMax = new Vector3(float.MinValue);
//             foreach (var n in nodes)
//             {
//                 bMin = Vector3.Min(bMin, n);
//                 bMax = Vector3.Max(bMax, n);
//             }
//         }

//         private static int CountApproxC5Nodes(List<Vector3> nodes, float tol)
//         {
//             // For each node, rotate 72° about Z and look for a nearby node
//             int count = 0;
//             double a72 = 2 * Math.PI / 5;
//             float cos72 = (float)Math.Cos(a72);
//             float sin72 = (float)Math.Sin(a72);

//             foreach (var n in nodes)
//             {
//                 float rx = n.X * cos72 - n.Y * sin72;
//                 float ry = n.X * sin72 + n.Y * cos72;
//                 // Find closest node to rotated position
//                 float bestD = float.MaxValue;
//                 foreach (var m in nodes)
//                 {
//                     float d = (float)Math.Sqrt((m.X-rx)*(m.X-rx)+(m.Y-ry)*(m.Y-ry)+(m.Z-n.Z)*(m.Z-n.Z));
//                     if (d < bestD) bestD = d;
//                 }
//                 if (bestD < tol) count++;
//             }
//             return count;
//         }
//     }
// }
