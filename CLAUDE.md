# PEAR-TREE — working rules for Claude

## No placeholders in commands — EVER

Never hand Jeremy a command containing a placeholder to fill in (`<DONORS>`,
`<SAME AS ...>`, `<PROJECT_ID>`, `/path/to/...`, `...`, etc.). A command with a
placeholder is not runnable and wastes a round-trip.

Instead, before giving the command:
1. **Find the real value** — from the conversation (e.g. an LSF `bacct`/`bjobs`
   stanza that prints the exact `FOFN=`/`OUTDIR=` the job ran with), the repo, a
   config, or a manifest.
2. If the value can't be found but *can be computed*, **compute it inline** inside
   the command (a `$(...)` that derives it from files that exist on the farm), so
   the command is self-contained and deterministic.
3. Only if a value is genuinely unknowable and uncomputable, **ask for it in a
   separate question first**, get the answer, then give the finished command —
   don't smuggle the unknown into the command as a placeholder.

Also give **full absolute paths**, never elisions like `/...` (see the
`full-paths-in-bash-commands` memory).
