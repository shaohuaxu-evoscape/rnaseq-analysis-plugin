# Setup Wizard — Stage 5: Confirm & Generate

Show a full summary table before writing:

```
┌─ Configuration Summary ───────────────────────────────────┐
│
│  Project
│    name       : {name}
│    organism   : {organism}  (taxid={taxid}, KEGG={kegg})
│    strain     : {strain}
│
│  Data
│    gene counts: {path}
│    genome     : {local/remote path}
│
│  Comparison
│    condition1 : {cond1}
│    condition2 : {cond2}
│    timepoints : {list}
│
│  Normalization: {method}  (min_expr={x}, min_samples={n})
│
│  Advanced     : {enabled modules or "none"}
│
│  [1]  Confirm and write configs/analysis_case.yaml   ←
│  [2]  Go back and edit a specific section
└───────────────────────────────────────────────────────────┘
Enter number:
```

If **[2]**, ask which section to go back to:
```
┌─ Go back to ──────────────────────────────────────────────┐
│  [1]  Project Basics          → re-read wizard-1-basics.md
│  [2]  Data Sources            → re-read wizard-2-data.md
│  [3]  Remote Preprocessing    → re-read wizard-3-remote.md
│  [4]  Normalization/Enrichment → re-read wizard-4-analysis.md
└───────────────────────────────────────────────────────────┘
Enter number:
```

On confirmation, write `configs/analysis_case.yaml` and suggest next commands:

```
Configuration saved. Next steps:
  python ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_run.py -c configs/analysis_case.yaml --steps 1a-5b --dry-run   # preview
  python ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_run.py -c configs/analysis_case.yaml --steps 1a-5b             # full run
  python ${CLAUDE_PLUGIN_ROOT}/scripts/rnaseq_run.py --list-steps                                             # list steps
```

See `references/rnaseq-run-pipeline.md` for execution details.
