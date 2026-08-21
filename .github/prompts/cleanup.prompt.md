---
description: "Safely remove dead code: unused functions or obsolete testing scripts"
agent: agent
---

Safely audit and remove dead code. **Never delete anything without explicit per-item approval.**

## Steps

1. **Read the refactor plan** at [.github/plans/refactor.md](./refactor.md) to see pre-identified deletion candidates. If the file does not exist, report that and stop.

2. **For each candidate** (file or function), show:
   - What exactly would be deleted (file path or function name + line range)
   - Any call sites found by grepping the workspace
   - A risk assessment (safe / has callers / unclear)

3. **Ask for explicit confirmation** before deleting each item. Do not batch deletions — confirm one item at a time.

4. **After each deletion**, grep the workspace for remaining references to the deleted symbol and flag any broken call sites.

5. **Never** skip the confirmation step, even if a function appears clearly unused.
