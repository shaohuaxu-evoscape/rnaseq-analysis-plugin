# Setup Wizard — Stage 1: Project Basics

Interactive project setup. Automatically triggered when `configs/analysis_case.yaml` does NOT exist.

## Interaction Style

**Use a menu-driven format throughout.** Whenever options can be enumerated, present them as a numbered list. The user selects by entering a number. Free-text input is used only when no finite list applies.

**Do NOT add any section headers, round labels, or step labels before menus or questions.** Just show the menu or question directly.

After each selection, echo confirmation before the next question:
```
✓ {Field} : {value}
```

**Menu format** (use consistently):
```
┌─ Title ───────────────────────────────────────────────────┐
│  [1]  Option A
│  [2]  Option B   ← recommended
│  [3]  Enter custom value
└───────────────────────────────────────────────────────────┘
Enter number:
```

Mark the recommended/default option with `← default` or `← recommended`.

## Project Scaffolding (automatic)

The project folder name was already determined in SKILL.md Case A (`$PROJECT_DIR`). Use that same name here — do NOT re-prompt. Create the directory structure inside it:

```bash
mkdir -p ${PROJECT_DIR}/configs ${PROJECT_DIR}/inputs/ref ${PROJECT_DIR}/results
cp ${CLAUDE_PLUGIN_ROOT}/configs/analysis_case.template.yaml ${PROJECT_DIR}/configs/analysis_case.yaml
```

Then inform the user:

```
Project scaffolded:
  {PROJECT_DIR}/configs/       — config file (we will fill this in together)
  {PROJECT_DIR}/inputs/        — place your data here (gene counts, reference genome, FASTQ)
  {PROJECT_DIR}/results/       — pipeline outputs will go here
```

> All subsequent file operations in this wizard use `{PROJECT_DIR}/` as the base path.

## Project Basics

**Project name (free text):**

```
请输入项目名称（自由文本，例如 "RNA Analysis Testcase"）：
```

**Organism (numbered menu).** Only show AFTER receiving the project name:

```
请选择实验物种：

┌─ Organism ────────────────────────────────────────────────┐
│  [1]  Yarrowia lipolytica
│  [2]  Saccharomyces cerevisiae
│  [3]  Escherichia coli
│  [4]  Homo sapiens
│  [5]  Mus musculus
│  [6]  Other (enter name)
└───────────────────────────────────────────────────────────┘
输入编号：
```

**Strain (free text, optional).** Only show after organism is confirmed:

```
请输入菌株名称（选填，例如 "A316"；直接回车跳过）：
```

Echo confirmation after all three are collected:
```
✓ Project  : {name}
✓ Organism : {organism}
✓ Strain   : {strain or "—"}
```

Update `project.name`, `project.organism`, `project.strain` in the config.

## Next

Read `references/wizard-2-data.md` and continue with Data Sources.
