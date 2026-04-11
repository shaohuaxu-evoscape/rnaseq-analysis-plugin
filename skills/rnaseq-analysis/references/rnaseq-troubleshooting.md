# Troubleshooting

> **Prerequisite:** Read [`../SKILL.md`](../SKILL.md) for pipeline overview and step routing.

Common errors and fixes for the RNA-seq analysis pipeline.

---

## 1. Singular Design Matrix (Step 3a)

**Error message**:
```
ValueError: Differential analysis step 3a design matrix is singular after adding
condition, timepoint, and batch fixed effects. Current batch/condition layout is
not identifiable.
```

**Cause**: Perfect collinearity among condition, timepoint, and batch in the OLS model (legacy fallback). Common scenarios:
- Only 1 shared timepoint (need at least 2)
- Batch and condition are completely confounded (all R1 from batch A, all R2 from batch B, no overlap)

**Fix**:
1. Ensure `experiment.timepoints` contains at least 2 timepoints
2. Check the batch/condition distribution in `sample_manifest.tsv`: ensure they're not completely confounded
3. For single-batch comparisons (most common), the model handles this automatically (batch_id has only one level, producing no extra column)

**Note**: With the DESeq2 Wald test (default), the design is auto-detected and this error is much less likely.

---

## 2. Missing Timepoints (Step 3a)

**Error message**:
```
ValueError: Differential analysis step 3a is missing samples for condition1 at
timepoints: [78]
```

**Cause**: A timepoint specified in `experiment.timepoints` has no samples for one condition.

**Fix**:
- **Option A**: Remove the missing timepoint from `experiment.timepoints`
- **Option B**: Delete the `experiment.timepoints` field entirely to let the pipeline auto-detect all shared timepoints

---

## 3. TPM/FPKM Gene Loss (Step 1a)

**Warning message**:
```
WARNING: N genes missing from GTF -- excluded from TPM
```

**Cause**: Some gene IDs in the count matrix have no exon annotations in the GTF file, so gene length cannot be computed. These genes are completely excluded.

**Impact**: If N is large (> 100), downstream analysis completeness may be affected.

**Fix**:
- Verify GTF matches the reference genome (same version)
- If gene length data is unreliable, switch to CPM: `normalization.method: "cpm"`
- Check whether missing genes are important heterologous/engineered genes (these are typically not in the original GTF)

---

## 4. No Shared Timepoints

**Error message**:
```
ValueError: Target conditions have no shared timepoints
```

**Cause**: The two target conditions have no overlapping timepoints in the source data.

**Fix**:
1. Verify gene_counts.tsv column names follow the format: `{condition}-{timepoint}` (e.g., `R1-18`, `R2-24`)
2. Confirm `condition` in `target_conditions` matches the column prefix
3. For multi-batch merges, verify each batch covers the target timepoints

---

## 5. Preprocessing Not Running

**Problem**: User wants to re-run preprocessing (fastp/HISAT2/htseq-count), but specifying `--steps 0a-0d` gives an error.

**Error message**:
```
ValueError: Preprocessing is now automatic and hidden. Run public analysis steps
only, for example '--steps 1a-2d' or '--steps all'.
```

**Cause**: Steps 0a-0d are hidden and cannot be manually specified. They only trigger automatically when a batch lacks a usable gene_counts.tsv.

**Fix**: To force re-preprocessing:
1. Comment out or remove the `gene_counts` path in the batch config:
   ```yaml
   batches:
     batch1:
       # gene_counts: "..."   # comment out or remove
       raw_fastq_dir: "/path/to/raw/fastq"
   ```
2. Delete the existing gene_counts.tsv at the default output path:
   ```bash
   rm results/shared/<batch_id>/01_preprocessing/gene_counts/gene_counts.tsv
   ```
3. Re-run the pipeline; preprocessing will trigger automatically

---

## 6. Network Request Failure (GO/KEGG Enrichment)

**Error message**:
```
RuntimeError: Failed to fetch https://rest.uniprot.org/... after 3 attempts: ...
```

**Cause**: UniProt or KEGG REST API is unreachable, or rate limiting was triggered.

**Fix**:
1. Check network connectivity
2. Steps 3b/3c cache downloaded data in `cache/` subdirectories. If partial data was downloaded:
   - Use `--resume` to skip completed steps
   - Or re-run the failed step directly
3. If failures persist, check `timeout` and `retries` params in `scripts/net.py`

**Cache locations**:
```bash
ls results/<run_name>/<pair_name>/03_differential_analysis/go_enrichment/cache/
ls results/<run_name>/<pair_name>/03_differential_analysis/kegg_enrichment/cache/
```

---

## 7. Resume Not Working

**Warning message**:
```
WARNING: Existing pipeline_run_manifest.json does not match the current config;
starting a fresh run
```

**Cause**: `--resume` checks config fingerprint (SHA-256 hash). If any field in `analysis_case.yaml` affecting the resolved config has changed, the fingerprint won't match.

**Fix**:
- To resume: revert config changes and re-run
- For a new config: don't use `--resume`, run from scratch (existing outputs will be overwritten)
- Changing `--steps` with the same config is fine; only the underlying config matters

---

## 8. Image Format / Quality Issues

**Problem**: Generated plots are blurry or in the wrong format.

**Fix**: Edit the `plot` section in `configs/analysis_case.yaml`:

```yaml
plot:
  dpi: 300       # Increase to 300 (publication quality)
  format: "png"  # Or "pdf" / "svg" (vector formats)
```

**Notes**:
- DPI 200 is sufficient for screen viewing; DPI 300 for publication/print
- PDF/SVG are vector formats, not limited by DPI, but produce larger files
- Changing `format` affects all steps; already generated plots won't update automatically -- re-run the relevant steps
