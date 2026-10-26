# Changelog

## Initial experimental public candidate

Prepared the connected studio in an isolated branch with portable launch/navigation, local coordinates, third-party notices and software checks. Private history and generated illustrations excluded. No engine changes. Publication remains gated by PUBLICATION_REVIEW.md.

## Code-only packaging

Excluded all runtime asset directories and community fields from Git. Added checksum-verified local import and a fail-closed launch check. Source bytes and notices are unchanged. Fresh clones require the separately held exact assets; no public bundle or download is implied.

## CI hygiene correction

Formatted and sorted the three Python support files using the pinned Ruff tool, replaced the test lambda and unittest-style assertions to match repository lint rules, and fixed two JSON final newlines. No runtime asset bytes, scientific data or checks were weakened.

## Main synchronization and Windows timeout diagnosis

Preserved both changelog entries while integrating main. The timed-out Windows job was executing long periodontal sweeps; main already marks these slow and runs them through nightly.yml, while default PR tests exclude slow. No test threshold, job timeout or scientific input was changed for this resolution. Fresh Windows CI is still required after push.
