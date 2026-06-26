//  PAPER TABLE 2 REFERENCE (Rosa et al. 2024):
//    IQC primal,  a=3.9mm, R=0.37mm, ρ≈10%
//    IQC dual,    a=3.9mm, R=0.41mm, ρ≈10%
//    IQC primal,  a=6.0mm, R=0.57mm, ρ≈10%
//    IQC dual,    a=6.0mm, R=0.63mm, ρ≈10%
//  Note: the auto-computed R from pi * r² L / V is approximate (doesn't account
//  for joint overlap).  The paper used CAD-derived ρ vs R/a curves (Fig.S7).
//  For exact sample geometry, use RelativeDensity=0 and set BeamRadius directly.
//
//  VOXEL SIZE GUIDE:
//    Dual beam dia ≈ 0.82mm at ρ=10%  →  voxel ≤ dia/3 ≈ 0.27mm
//    0.30mm — fast, good for iteration
//    0.20mm — higher quality
//    0.13mm — very high quality, slow (~10× vs 0.30mm)
// =============================================================================

using System;
using System.IO;
using PicoGK;
using IcosahedralQC;

// ── Output folder ─────────────────────────────────────────────────────────────
string strLogFolder = Path.Combine(
    Environment.GetFolderPath(Environment.SpecialFolder.MyDocuments),
    "IQC_Output");
Directory.CreateDirectory(strLogFolder);

// ── Parameters ────────────────────────────────────────────────────────────────
IQCRenderer.Params = new IQCParameters
{
    // Dual or Primal
    // Mode = TrussMode.Dual,
    Mode = TrussMode.Dual,

    // ── Topology ───────────────────────────────────────────────────────────────
    EdgeLength = 6.0f,
    NShell = 6, // 0 triggers autosizing, however 6 fills sufficiently
    WindowScale = 1.0f,

    // ── Beam sizing ─────────────────────────────────────────────────────────────
    // Auto-compute from RelativeDensity.  Paper target: rho = 10%.
    // I have mixed feelings on this. On one hand, from a visual perspective
    // I am nailing the geometry on the head. On the other hand, it also says
    // that the way I am computing the relative density leads to a different
    // relative density with the same beam size. For now, I would recommend
    // just using a fixed beam radius and setting relative density to 0
    // RelativeDensity = 0.10f,
    RelativeDensity = 0.10f,
    BeamRadius = 0.37f,         // only used when RelativeDensity == 0
    MinBeamRadius = 0.15f,         // clamp floor; RT test gives larger R -> loosen this

    // ── Domain ─────────────────────────────────────────────────────────────────
    // Paper is symmetric about center
    DomainSize = 100.0f,

    // ── Post-processing ─────────────────────────────────────────────────────────
    SmoothDistMM = 0.15f,
    SmoothPasses = 0,

    // ── Output ─────────────────────────────────────────────────────────────────
    StlFileName = "IQC_DualTruss_a6mm_rho10_plated_temp",
    // StlFileName = "IQC_DualTruss_a6mm_rho10_plated_fixed_radius",
    // StlFileName = "IQC_DualTruss_a39mm_rho10_plated",
    // StlFileName = "IQC_DualTruss_a39mm_rho10_plated_fixed_radius",
    // StlFileName = "IQC_PrimalTruss_a6mm_rho10_plated",
    // StlFileName = "IQC_PrimalTruss_a6mm_rho10_plated_fixed_radius",
    // StlFileName = "IQC_PrimalTruss_a39mm_rho10_plated",
    // StlFileName = "IQC_PrimalTruss_a39mm_rho10_plated_fixed_radius",
};

// ── Voxel size ────────────────────────────────────────────────────────────────
float fVoxelSizeMM = 0.2f;

// ── Run ───────────────────────────────────────────────────────────────────────
try
{
    PicoGK.Library.Go(
        fVoxelSizeMM,
        IQCRenderer.Task,
        strLogFolder);
}
catch (Exception e)
{
    Console.WriteLine("Task failed: " + e.Message);
    Console.WriteLine(e.ToString());
}
