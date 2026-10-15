"""Text summaries printed after a run.

The wording is part of the command-line contract: tests and the examples
match these lines.
"""

from __future__ import annotations

from collections.abc import Sequence

from marse.core.reactive_transport import ReactiveTransportResult
from marse.core.simulation import BiofilmProfileResult, SimulationResult
from marse.core.well_mixed import WellMixedResult


def _summarise(result: SimulationResult | BiofilmProfileResult) -> None:
    print(f"run_id      {result.manifest.run_id}")
    print(f"experiment  {result.config.experiment_id}")
    print(f"kind        {result.config.kind}")

    if isinstance(result, BiofilmProfileResult):
        outputs = result.manifest.outputs
        settings = result.config.biofilm
        if settings is None:
            raise ValueError("a biofilm profile result is missing its biofilm settings")
        print(f"biofilm     {settings.thickness_um:g} um in {settings.cells} nodes")
        print(f"solved      {outputs['newton_iterations']} Newton iterations")
        print(f"penetration {outputs['penetration_depth_um']:.1f} um")
        print(f"active zone {outputs['active_zone_um']:.1f} um")
        for name, rate in outputs["mean_growth_rate_per_h"].items():
            surface = outputs["surface_growth_rate_per_h"][name]
            share = rate / surface if surface else 0.0
            print(f"  {name:<24} mean {rate:.6g} /h ({share:.1%} of the surface rate)")
        return

    final = result.final_state
    print(f"steps       {final.step} over {final.time_h:g} h")
    print(
        f"substrate   {final.substrate_mm:.6g} mM left of {result.config.substrate.initial_mm:g} mM"
    )
    for name, value in zip(result.organism_names, final.biomass_g_per_l, strict=True):
        print(f"  {name:<24} {value:.6g} g/L")


_BIOFILM_KEYS = ("penetration_depth_um", "active_zone_um", "base_concentration_mm")


def _summarise_well_mixed(result: WellMixedResult) -> None:
    outputs = result.manifest.outputs
    substeps = outputs["substeps"]
    print(f"run_id      {result.manifest.run_id}")
    print(f"experiment  {result.config.experiment_id}")
    print(f"kind        {result.config.kind}")
    print(
        f"steps       {result.config.steps} over {outputs['final_time_h']:g} h: "
        f"{substeps['accepted']} adaptive substeps, {substeps['rejected']} retried, "
        f"{substeps['limited']} limited to keep concentrations positive"
    )
    print("final       mol per m3")
    for name, value in outputs["final_mol_per_m3"].items():
        print(f"  {name:<24} {value:.6g}")
    if "ph" in outputs:
        ph = outputs["ph"]
        print(
            f"pH          {ph['final']:.3f} at the end; lowest {ph['lowest']:.3f} "
            f"at {ph['lowest_at_h']:g} h"
        )
    worst = max(q["largest_relative_residual"] for q in outputs["balance"].values())
    print(f"balance     {_listed(outputs['balance'])} conserved to {worst:.1e} of their totals")


def _listed(quantities: Sequence[str]) -> str:
    """``carbon, nitrogen and electrons``, from the quantities a balance lists."""
    names = list(quantities)
    return names[0] if len(names) == 1 else ", ".join(names[:-1]) + " and " + names[-1]


def _summarise_in_space(result: ReactiveTransportResult) -> None:
    outputs = result.manifest.outputs
    substeps = outputs["substeps"]
    grid = result.config.domain.grid
    size = " x ".join(f"{s:g}" for s in grid.size_um)
    voxels = " x ".join(str(n) for n in grid.shape)
    print(f"run_id      {result.manifest.run_id}")
    print(f"experiment  {result.config.experiment_id}")
    print(f"kind        {result.config.kind}")
    print(f"space       {grid.dimensions}-D, {voxels} voxels of {grid.voxel_um:g} um ({size} um)")
    print(
        f"steps       {result.config.steps} over {outputs['final_time_h']:g} h: "
        f"{substeps['accepted']} implicit substeps, {substeps['rejected']} retried, "
        f"{substeps['limited']} limited to keep concentrations positive"
    )
    faces = ["the film" if "mouth" in outputs else "the air" if "air" in outputs else "the top"]
    if "mouth" in outputs and "air" in outputs:
        faces.append("the air")
    if "surface" in outputs:
        faces.append("the substratum")
    boundary = _listed(faces)
    if "mouth" in outputs and any(intake.retained for intake in result.config.domain.diet):
        boundary += " (food left on the teeth included)"
    print(f"{'final':<24}   mol per m2 of surface   entered through {boundary}")
    for name, value in outputs["totals_mol_per_m2"].items():
        entered = outputs["imported_mol_per_m2"][name]
        print(f"  {name:<24} {value:>12.6g}   {entered:>+12.4g}")
    if "surface" in outputs:
        _summarise_surface(result)
    if "mouth" in outputs:
        _summarise_mouth(result)
    if "spreading" in outputs:
        spreading = outputs["spreading"]
        final = spreading["final"]
        print(
            f"spreading   {spreading['mechanism']}, every {spreading['interval_h']:g} h: "
            f"{spreading['spreads']} spreads, at most {spreading['most_rounds']} round(s) each, "
            f"{spreading['moved_um3_per_um2']:.4g} um3 per um2 moved"
        )
        print(
            f"  biovolume {final['biovolume_um3_per_um2']:.4g} um3 per um2, maximum thickness "
            f"{final['maximum_thickness_um']:g} um, "
            f"fullest voxel {final['largest_volume_fraction']:.6f}"
        )
    if "plaque" in outputs:
        _summarise_plaque(result)
    if "air" in outputs:
        _summarise_air(result)
    if "ph" in outputs:
        ph = outputs["ph"]
        final = ph["final"]
        print(
            f"pH          at the substratum {final['substratum_mean']:.3f} at the end "
            f"(lowest {ph['lowest_at_substratum']:.3f} at {ph['lowest_at_substratum_h']:g} h); "
            f"in the box {final['box_min']:.3f} to {final['box_max']:.3f}"
        )
    worst = max(q["largest_relative_residual"] for q in outputs["balance"].values())
    print(
        f"balance     {_listed(outputs['balance'])}, counting what crossed {boundary}, "
        f"conserved to {worst:.1e}"
    )


def _summarise_mouth(result: ReactiveTransportResult) -> None:
    mouth = result.manifest.outputs["mouth"]
    worst = max(q["largest_relative_residual"] for q in mouth["balance"].values())
    print(
        f"mouth       {mouth['swallows']} swallows in {mouth['spans']} spans; "
        f"{mouth['final_volume_ml']:.3g} mL at the end"
    )
    print(f"  {'saliva at the end':<24} mol per m3")
    for name, value in mouth["final_mol_per_m3"].items():
        print(f"  {name:<24} {value:>12.6g}")
    if result.mouth is not None and "ph" in result.mouth:
        print(f"  {'pH':<24} {result.mouth['ph'][-1]:>12.3f}")
    counted = ["secreted", "swallowed"]
    if "intakes" in mouth:
        counted = ["secreted", "eaten", "swallowed"]
        print(f"  {'taken in':<24} mol per m2 of plaque, from {mouth['intakes']} intake(s)")
        for name, value in mouth["eaten_mol_per_m2"].items():
            if value:
                print(f"  {name:<24} {value:>12.6g}")
    if "expelled_mol_per_m2" in mouth:
        counted.append("expelled")
    if "air" in result.manifest.outputs:
        counted.append("exchanged with the air")
    print(
        f"  box and mouth together, counting what was {_listed(counted)}, conserved to {worst:.1e}"
    )


def _summarise_plaque(result: ReactiveTransportResult) -> None:
    plaque = result.manifest.outputs["plaque"]
    assert result.plaque is not None
    thickness = result.plaque["thickness_um"]
    print(
        f"plaque      {plaque['final_thickness_um']:.4g} um at the end, "
        f"{thickness.min():.4g} to {thickness.max():.4g} um over the run "
        f"(spreading {plaque['spreading']})"
    )
    detached = ", ".join(f"{n} {v:.4g}" for n, v in plaque["detached_mol_per_m2"].items() if v)
    print(f"  detached, worn and over the top, mol per m2: {detached or 'nothing'}")
    if "removed_mol_per_m2" in plaque:
        cleanings = result.manifest.outputs["mouth"].get("cleanings", 0)
        removed = ", ".join(f"{n} {v:.4g}" for n, v in plaque["removed_mol_per_m2"].items() if v)
        print(f"  taken off by {cleanings} cleaning(s), mol per m2: {removed or 'nothing'}")


def _summarise_air(result: ReactiveTransportResult) -> None:
    air = result.manifest.outputs["air"]
    gave = ", ".join(f"{n} {v:.4g}" for n, v in air["exchanged_mol_per_m2"].items())
    into = "the box and the mouth" if "mouth" in result.manifest.outputs else "the box"
    print(f"air         gave {gave} mol per m2 to {into}, less what it took back")


def _summarise_surface(result: ReactiveTransportResult) -> None:
    domain = result.config.domain
    final = result.manifest.outputs["surface"]["final"]
    substratum = domain.substratum
    assert substratum is not None
    species = [s.attached for s in domain.suspension]
    print(f"surface     {domain.surface.conditioning_film if domain.surface else ''}")
    print(f"  {'bound cells per cm2':<22}" + "".join(f"{s:>14}" for s in species) + "   covered")
    for material in substratum.materials:
        cells = "".join(f"{final[f'{material}_{s}_cells_per_cm2']:>14.4g}" for s in species)
        print(f"  {material:<22}{cells}   {final[f'{material}_covered']:>6.1%}")
