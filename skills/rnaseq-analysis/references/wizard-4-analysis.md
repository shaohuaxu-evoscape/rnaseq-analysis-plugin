# Setup Wizard — Stage 4: Normalization, Enrichment & Advanced Analysis

## Normalization

```
┌─ Normalization Method ────────────────────────────────────┐
│  [1]  TPM   — Transcripts Per Million (requires GTF)   ← default
│  [2]  CPM   — Counts Per Million (no GTF needed)
│  [3]  FPKM  — Fragments Per Kilobase Per Million
└───────────────────────────────────────────────────────────┘
Enter number:
```

Then show current gene filtering defaults and ask to confirm or edit:
```
Gene filtering defaults:
  min_expr    = 1.0   (minimum normalized expression)
  min_samples = 2     (minimum samples above threshold)

┌──────────────────────────────────────────────────────────┐
│  [1]  Keep defaults   ← recommended
│  [2]  Edit values
└──────────────────────────────────────────────────────────┘
Enter number:
```

Echo:
```
✓ Normalization : {method}
✓ Gene filtering: min_expr={x}, min_samples={n}
```

Update `normalization.method` and `normalization.gene_filtering`.

## Enrichment Analysis

```
┌─ Organism for Enrichment ─────────────────────────────────┐
│  Select organism to auto-fill taxid and KEGG code:
│
│  [1]  Yarrowia lipolytica    taxid=4952   KEGG=yli
│  [2]  Saccharomyces cerevisiae  taxid=4932  KEGG=sce
│  [3]  Escherichia coli K-12  taxid=83333  KEGG=eco
│  [4]  Homo sapiens           taxid=9606   KEGG=hsa
│  [5]  Mus musculus           taxid=10090  KEGG=mmu
│  [6]  Enter manually
│  [7]  Skip (fill in later)
└───────────────────────────────────────────────────────────┘
Enter number:
```

If **[6]**, ask for taxid and KEGG code separately as free text.

Echo:
```
✓ Organism taxid : {taxid}
✓ KEGG code      : {kegg}
✓ Gene ID prefix : {prefix}  (e.g., YALI1_ for Y. lipolytica)
```

Update `project.organism_taxid`, `project.gene_id_prefix`, `differential_analysis.go_enrichment.organism_taxid`, `differential_analysis.kegg_enrichment.organism`.

## Advanced Analysis (Optional)

```
┌─ Advanced Analysis ───────────────────────────────────────┐
│  All modules below are disabled by default.
│  Enter numbers to enable (comma-separated), or 0 to skip:
│
│  [1]  Heterologous gene expression (4c)
│  [2]  Fermentation overview (4d)
│  [3]  Transporter & aminopeptidase analysis (4a)
│  [4]  Temporal causality (4b)
└───────────────────────────────────────────────────────────┘
Enable (e.g., 1,3), or 0 to skip all:
```

Parse the user's reply as a comma-separated list of numbers to enable. If "0" or empty, keep all disabled.

For any module toggled on, ask for required parameters:
- **4c (Heterologous genes)**: ask for gene list as free text (name, function, color per gene)
- **4d (Fermentation overview)**: ask for data file path
- **4b (Temporal causality)**: ask for gene module lists (TORC1, sterol, etc.)

Echo:
```
✓ Heterologous genes  : {enabled/disabled}
✓ Fermentation data   : {enabled/disabled}
✓ Transporter analysis: {enabled/disabled}
✓ Temporal causality  : {enabled/disabled}
```

Update the relevant `advanced_analysis.*.enabled` flags.

## Next

Read `references/wizard-5-confirm.md` and present the configuration summary.
