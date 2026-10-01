# Surface adhesion and detachment

> This page describes a rule of the version 1 two-dimensional ecosystem engine.
> Cells binding to a substratum in the version 2 engine in space are a
> different model: see [environments.md](environments.md) and
> [theory.md §6.4](theory.md#64-attachment-to-surfaces).

Species can opt into a simple surface interaction with `adhesion_edges`,
choosing any of `top`, `bottom`, `left`, or `right`. `adhesion_per_h` transfers
biomass from adjacent interior cells onto those surfaces. `detachment_per_h`
removes biomass from the selected surface with first-order kinetics.

This is a bounded phenomenological surface model, not a mechanical biofilm
solver. It does not infer attachment forces, matrix rheology, shear stress, or
cell shape. With no adhesion edges or zero rates, the default ecosystem
behavior is unchanged.
