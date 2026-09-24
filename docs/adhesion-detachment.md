# Surface adhesion and detachment

Species can opt into a simple surface interaction with `adhesion_edges`,
choosing any of `top`, `bottom`, `left`, or `right`. `adhesion_per_h` transfers
biomass from adjacent interior cells onto those surfaces. `detachment_per_h`
removes biomass from the selected surface with first-order kinetics.

This is a bounded phenomenological surface model, not a mechanical biofilm
solver. It does not infer attachment forces, matrix rheology, shear stress, or
cell shape. With no adhesion edges or zero rates, the default ecosystem
behavior is unchanged.
