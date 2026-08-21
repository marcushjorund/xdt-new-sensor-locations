---
description: "Add a new sensor selection result file to the Shiny app and deploy script"
agent: agent
argument-hint: "Filename in results/ and display label for the UI dropdown"
---

Add a new result configuration. This is an **atomic operation** — all three steps must succeed or none should be applied.

## Steps

1. **Verify** the named file exists in `results/`. If it does not exist, stop and report the missing file — do not proceed with edits.

2. **Add an entry to `configs`** in [shiny/app.R](../../shiny/app.R). Follow the exact pattern of existing entries in the `configs` named list (label → filename string).

3. **Add the filename to `rds_files`** in [deploy.R](../../deploy.R). Add it to the same vector as all other result filenames.

## Rules

- Steps 2 and 3 are inseparable. Never do one without the other.
- Do not change any other code in `app.R` or `deploy.R`.
- After completing, confirm: "Added `<filename>` to configs in app.R and rds_files in deploy.R."
