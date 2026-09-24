# Action and resource economy

`ActionRule` and `ResourceBudget` make the game-like layer explicit. An action
has a name, a resource cost, and optionally a required capability. `try_action`
returns a new budget when the action executes, or a structured reason when it
is blocked by a missing capability or insufficient resource.

This is a bounded accounting layer, not a claim that all microbial or immune
decisions reduce to one energy pool. Experiments can define ATP-like energy,
oxygen, nutrients, drug-resistance capacity, or other measured resources, and
later providers can connect these budgets to spatial fields.
