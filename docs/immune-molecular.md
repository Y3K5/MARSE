# Immune and molecular interaction primitives

`marse.immune` provides explicit field-level primitives for immune pressure
and molecular neutralization. An `ImmuneInteraction` maps an effector field to
a species-specific first-order biomass loss using a bounded Hill response.
`MolecularNeutralizer` reduces the effective effector field, representing a
declared competing molecule, inhibitor, or drug interaction.

`ImmuneCellType` and `ImmuneAgent` provide a small game-like action layer:
`move_toward` follows the steepest adjacent biomass gradient, `attack` removes
target biomass with first-order kinetics, and `secrete` adds an effector field.
Each action has explicit rates and energy scaling, so resistance can be
represented by lower susceptibility in `ImmuneInteraction` rather than by an
implicit rule.

These are not agent-based immune cells, receptor kinetics, pharmacokinetics,
or clinical predictions. Fields must be supplied by the caller, units and
sources remain part of the experiment definition, and the output reports
removed biomass explicitly. Agent trajectories, contact mechanics, cytokine
networks, and host-tissue models remain future opt-in layers.
