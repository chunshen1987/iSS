# iSS repository instructions

## Scientific scope

iSS (iSpectraSampler) is scientific computing software for high-energy particle collisions. It samples hadrons from hydrodynamic freeze-out surfaces using the Cooper–Frye prescription and can pass those samples to hadronic afterburners such as UrQMD.

When changing a physics or numerical result, state the relevant assumptions, units, and conventions. Verify the change with the applicable tests or a reproducible comparison. For stochastic comparisons, record the seed or the statistical criteria used to assess agreement.

## Development reports

For every development task, write a Markdown summary report under `doc/` before handing off the work. Describe what changed, why it changed, how it was verified, and any remaining limitations. Include reproduction commands, inputs, and settings when they are relevant to the results. Link to the report in the final response.

## Git commits

Use a concise commit subject and include a brief summary of the changes in the commit body. End agent-authored commits with a `Signed-off-by` trailer using the AI agent's name. For example, Codex signs with `Signed-off-by: Codex`. Use the agent's email in the trailer when one is configured for that agent.
