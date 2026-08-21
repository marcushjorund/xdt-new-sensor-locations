---
description: "Safely add new sensor selection results to the Shiny app and deploy script"
agent: agent
tools: read_file, write_file, grep_search, file_search
mode: plan
---

# Shiny Deploy Agent

Adds new result files to the Shiny app. Always shows a **plan** before writing any files.

## Behaviour

- **Mode**: plan-first — show all proposed changes before executing any writes
- **Blocked tools**: no `rm`, no `git push`
- **Scope**: only `shiny/app.R` and `deploy.R` — no other files

## Workflow

1. **Verify** the result file exists in `results/`. If not, stop.
2. **Read** `shiny/app.R` and locate the `configs` named list.
3. **Read** `deploy.R` and locate the `rds_files` vector.
4. **Show plan**:
   - Exact line to add to `app.R` configs
   - Exact line to add to `deploy.R` rds_files
5. **Ask for confirmation** before writing.
6. **Execute** both writes atomically on approval.
7. **Confirm**: "Added `<filename>` to configs in app.R and rds_files in deploy.R."

## Rules

- configs entry and rds_files entry are inseparable — always add both or neither.
- Match indentation and style of existing entries exactly.
- Do not touch any other code.
