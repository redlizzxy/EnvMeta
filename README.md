# EnvMeta

**An environmental microbiome metagenomic visualization platform**

[中文版 README](README_CN.md) · [Online demo](https://envmeta-3xjhcu7lv2gkj4pjtk8gsb.streamlit.app/) · [License](LICENSE)

EnvMeta solves the core pain point of environmental-microbiome PhD students: **the sequencing vendor hands you a pile of TSVs and you don't know which file lets you do which analysis**. One-click file recognition + 14 publication-grade plot types + automatic biogeochemical-cycle figure inference + a YAML hypothesis scorer + standalone interactive HTML export. **Fully open-source, free, offline-capable.**

> 🌐 **Try it online (zero install)**: <https://envmeta-3xjhcu7lv2gkj4pjtk8gsb.streamlit.app/>
>
> Click "📦 Load arsenic-slag remediation example data" on the home page → all 14 analyses run in 3 seconds.

---

## Why we built it

| Real pain point | Traditional path | EnvMeta |
|---|---|---|
| Vendor returns a stack of `.tsv`, unclear what each one is for | Write your own parser | **Drop it in → auto-recognized** as 1 of 11 file types |
| Want stacked-bar / PCoA / LEfSe but can't write R | Pay vendor for premium analysis / grind through code | **3 clicks → publication-grade PDF** |
| Need an element-cycle figure to argue Fe-As-S coupling | No public tool exists; people draw it by hand in PowerPoint | **Auto-inferred from KO annotations**: 4 elements × 18 pathways |
| Want to assess whether the data actually supports your mechanistic hypothesis | Subjective narrative | **YAML scorer** + permutation null-p + weight-sensitivity |
| Reviewer wants to reproduce SI | PDF attachment + data tables + email back-and-forth | **400 KB standalone HTML**, reviewer interacts in their browser |

## Selling points (vs. competitors)

| Capability | Krona | Anvi'o | MicrobiomeAnalyst | Vendor cloud platforms | **EnvMeta** |
|---|:-:|:-:|:-:|:-:|:-:|
| Auto-inferred element-cycle figure | ❌ | ❌ | ❌ | ❌ | **✅ unique** |
| Hypothesis scoring + null-p + weight sensitivity | ❌ | ❌ | ❌ | ❌ | **✅ unique** |
| Standalone offline interactive HTML (SI killer) | ❌ | ⚠️ static | ❌ web-only | ❌ | **✅ unique** |
| Fork Bundle (paper-tool binding) | ❌ | ❌ | ❌ | ❌ | **✅ unique** |
| Cross-element coupling (As ↔ H₂S → As₂S₃) | ❌ | ❌ | ❌ | ❌ | **✅ unique** |
| All parameters tweakable | ⚠️ | ✅ | ⚠️ | ❌ locked | ✅ |
| Fully open-source, free | ✅ | ✅ | ✅ | ❌ paid | ✅ |

## Feature matrix (v0.9.0 / hypothesis stress test complete)

| Module | What it covers |
|---|---|
| 📁 File recognition | metadata / abundance (MAG / Taxon layered) / distance / alpha / CheckM / env / KO wide+long / Keystone / MAG taxonomy / Gephi nodes+edges — **11 types** |
| 📊 Reads-based (7 figures) | Taxonomy stacked bar / α-diversity / β-diversity PCoA / RDA ordination / LEfSe / element-cycle gene heatmap / gene log2FC |
| 🧬 MAG-based (5 figures) | MAG quality scatter / MAG abundance heatmap / pathway completeness / element-cycle gene profile / co-occurrence network (Gephi-prep) |
| 🔄 Biogeochemical-cycle figure ⭐ | 4 elements (As/N/S/Fe) × 18 pathways auto-inferred + cross-element chemical coupling + keystone ★ annotation |
| 🧪 Mechanistic-hypothesis YAML scorer | **6 claim types** (pathway_active / **`pathway_inactive`** [v0.9 ⭐ Popperian falsification] / coupling_possible / env_correlation / keystone_in_pathway / group_contrast) + 3 confidence indicators (Fisher permutation p / Saltelli weight-sensitivity / Bradford-Hill required veto) + 9-tier interpretation |
| 📐 Hypothesis writing guide (v0.9 ⭐) | Two-layer template (calibration + stress claims) + pre-registration discipline + pre-prediction template + 6 claim-type selection guide + Bradford-Hill mapping. See [`docs/hypothesis_writing_guide.md`](docs/hypothesis_writing_guide.md) |
| 📦 Fork Bundle | Pack KB + YAML + config + KEGG snapshot → zip; reviewers reproduce in one click |
| 🌐 Standalone interactive HTML | 400 KB single file with D3.js inlined, 4-quadrant force layout + click-through + SVG/JSON export, fully offline |
| 💾 Export center | PNG / PDF / SVG / TIFF 600 dpi / TSV / runnable `.py` reproduction script, batch ZIP |
| 🧭 Beginner-onboarding kit | Data-prep guide (11 upstream tools → EnvMeta mapping) + chart-selection wizard (question → recommended analysis) + "how to read" expander on each of the 14 figures + sample data one-click load |

## Target users

- 🎓 Master's/PhD students in **environmental microbiology / soil remediation / ecology**
- 📚 Labs that send samples to a sequencing vendor and want to run downstream analyses themselves
- 📝 Authors who need an auto-inferred element-cycle figure as a mechanism figure
- 🔬 No programming background required; users with R/Python experience get deeper customization

## Three ways to use it

### 1. Online, zero install (easiest)

<https://envmeta-3xjhcu7lv2gkj4pjtk8gsb.streamlit.app/>

If the page first shows "App is sleeping", click "Wake up" and wait ~30 s.

### 2. Local install (beginner-friendly)

See [docs/install_for_beginners.md](docs/install_for_beginners.md) (Chinese): zero-base 20-min install (Miniconda + clone repo + one-line dependency install).

### 3. Local install (familiar with bioinformatics)

```bash
git clone https://github.com/redlizzxy/EnvMeta.git
cd EnvMeta
conda create -n envmeta python=3.11 -y
conda activate envmeta
pip install -r requirements.txt
streamlit run app.py
```

The browser auto-opens `http://localhost:8501`.

## 📜 Recent releases

Beta phase ships frequent bug fixes / features. Full list in **[CHANGELOG.md](CHANGELOG.md)**.

### v0.9.0 — 2026-05-09 (hypothesis scoring controlled experiment complete + stress-test discrimination evidence) ⭐

**Paper 3 core evidence fully in place.** A controlled experiment over 4 KEGG-curated metagenomic datasets (in-house + Liu 2023 cold seep + Grettenberger 2021 AMD stream + Ayala 2020 pit lake) with all hypothesis YAMLs pre-registered (git timestamp anchored before EnvMeta runs).

- ✨ **New `pathway_inactive` claim type** — the 6th claim type, Popperian falsifiability primary instrument. Evaluates as `satisfied` when n_active_mags == 0 (matches the "should NOT be active" prediction); backward-compatible with all existing YAMLs.
- ✨ **Two-layer hypothesis writing tutorial** — [`docs/hypothesis_writing_guide.md`](docs/hypothesis_writing_guide.md) for users (calibration + stress dual-layer template, pre-registration discipline, pre-prediction template, 6 claim-type guide, Bradford-Hill mapping) + [`paper/hypotheses/HYPOTHESIS_DESIGN_PRINCIPLES.md`](paper/hypotheses/HYPOTHESIS_DESIGN_PRINCIPLES.md) for paper Methods (design philosophy, 4 categories of limitations).
- ✨ **General stress runner CLI**: [`tools/external_benchmarks/run_stress_yaml.py`](tools/external_benchmarks/run_stress_yaml.py) takes any `--dataset` + `--yaml` and runs cycle_diagram + scoring; supports `--all`.
- 📊 **4 KEGG-curated datasets all STRONG (calibration)** + **3 stress tests with score gaps below calibration** (Grettenberger weak 0.250 / Liu suggestive 0.625 / Ayala suggestive 0.455). Cross-topic `arsenate_reduction_should_dominate` rejected in **2/2 non-arsenic datasets** (n=0 active MAGs in both Grettenberger and Ayala) — ironclad evidence that EnvMeta's scoring engine is domain-neutral.
- 📚 **Reference DOI audit** — verified DOIs for all 16 hypothesis claims × 13 review citations; corrected 4 citation errors transparently (Yin 2011 wrong journal; Bothe "2007 FEMS Rev" non-existent → Bothe 2000; Cabrera 2006 wrong journal+topic; Auld 2017 wrong topic → replaced by Dai 2014 PLoS One [primary AMD nifHDK metagenomic evidence] + Méndez-García 2015 Front Microbiol review).
- 🐛 6 hypothesis YAMLs reference metadata corrected (no claim entity changed; pre-registration audit trail preserved in git history).
- 🧪 pytest **297/297 green** (+4 new `pathway_inactive` test cases).

### v0.8.2 — 2026-05-08 (RDA values aligned to R vegan + 11-figure side-by-side validation complete)

- 🐛 Fixed RDA values disagreeing with R `vegan` (skbio normalization difference caused 16-20× inertia bias and ANOVA F/p reversal)
  - Switched to SS-based formulas (vegan-equivalent); after fix all F / r / explained-variance match R to 4 decimal places
- 📚 Completed R/Python side-by-side validation for all 11 figures (`paper/benchmarks/validation/`)
  - 5 figures with exact numerical agreement + 6 with algorithmic equivalence + 11 READMEs + paper citation template
- 📚 Paper 3 (EnvMeta methodology paper) pre-submission task tracker (`paper/manuscript/`)
- 🧪 pytest 293/293 green (no regressions)

### v0.8.1 — 2026-04-21 (first batch of macOS beta feedback fixes)

- 🐛 Fixed interactive HTML losing chemical–pathway links after switching groups (drag/hover broken)
- 🐛 Fixed YAML hypothesis upload on macOS raising `[Errno 63] File name too long`
- 📚 Beginner install guide adds 3 macOS gotchas (conda ToS / protobuf resolver / Streamlit welcome email)
- 📚 README adds "Updating to a new version" section

### v0.8.0 — 2026-04-19 (v0.8 beta, "Sunday Sprint")

- ✨ Beginner-onboarding kit (data-prep guide / chart wizard / 14-figure how-to-read / one-click sample load)
- ✨ Export center: 4-tab unified entry + batch ZIP
- ✨ Interactive HTML export v1.3 (D3 inline + standalone offline)
- ✨ Streamlit Cloud online deployment

---

## 🔄 Updating to a new version (beta phase)

EnvMeta ships **frequent fixes / features** during beta — pulling once a week is recommended.

**Online users** need do nothing — Streamlit Cloud redeploys automatically.

**Local install users**, run in Terminal (Mac) or Anaconda Prompt (Windows):

```bash
# 1) If streamlit is running, Ctrl+C to stop it
# 2) Enter the EnvMeta directory
cd ~/Desktop/EnvMeta              # Mac
# cd %USERPROFILE%\Desktop\EnvMeta  # Windows

# 3) Activate the environment
conda activate envmeta

# 4) Pull latest code (this is the core update command)
git pull origin master

# 5) Only when the terminal mentions requirements changes / startup throws ModuleNotFoundError
pip install -r requirements.txt

# 6) Restart
streamlit run app.py
```

### Common scenarios

**A: `git pull` says `Already up to date.`** — You're current; just start the app.

**B: `git pull` reports `Your local changes ... would be overwritten`** — You've edited files locally. Stash, pull, optionally restore:

```bash
git stash              # stash local changes
git pull origin master
git stash pop          # restore them if you want
```

**C: After update, startup throws `No module named 'xxx'`** — Run `pip install -r requirements.txt` to install new dependencies.

### How to know there is an update

- Check the [GitHub repo home page](https://github.com/redlizzxy/EnvMeta) commit list — new commits = new updates
- Or `git fetch && git log HEAD..origin/master --oneline` to see how many commits behind you are
- Pull first when you hit a bug — it may already be fixed

## Data formats

Read [docs/data_preparation_zh.md](docs/data_preparation_zh.md) (also browsable in-app on the "Data preparation guide" page; currently Chinese, English version planned). It covers output → EnvMeta input mapping for 11 upstream tools: CoverM / HUMAnN3 / eggNOG / DRAM / GTDB-Tk / CheckM2 / KofamScan / QIIME2 / Kraken2+Bracken / MetaPhlAn4 / iCAMP.

## Sample data

`tests/sample_data/` contains a slim version of the arsenic-slag / steel-slag microbial-remediation study (3 groups CK/A/B, 168 MAGs × 10 samples × 57 target KOs). One click on the home page runs all 14 analyses end-to-end.

## Tech stack

- **UI front-end**: Streamlit
- **Publication-grade rendering**: matplotlib + seaborn
- **Statistics**: pandas / scipy / scikit-bio / statsmodels
- **Interactive HTML**: D3.js v7 (inlined; no external deps)
- **Cycle-figure inference**: in-house (`envmeta/geocycle/`) — KEGG-driven KB + permutation tests + sensitivity scan

## Project layout

```
envmeta/
├── app.py                         # Streamlit entry point
├── envmeta/                       # Core package
│   ├── file_manager/              # File recognition (11 types)
│   ├── analysis/                  # 14-figure analysis engine
│   ├── geocycle/                  # Cycle inference + hypothesis scorer + HTML export
│   │   ├── knowledge_base/        # 4 elements × 18 pathways × 57 KOs (KEGG-driven)
│   │   ├── hypothesis.py          # YAML scorer
│   │   └── html_exporter.py       # Standalone interactive HTML
│   ├── help/                      # Beginner-onboarding kit
│   ├── tools/                     # Fork Bundle / KB builder / Gephi prep
│   └── export/                    # PNG/PDF/SVG/TIFF + .py reproduction script
├── docs/
│   ├── data_preparation_zh.md     # Upstream tool → EnvMeta mapping
│   └── install_for_beginners.md   # Beginner install guide (Chinese)
├── paper/
│   ├── bundles/                   # Example paper Fork Bundles
│   ├── benchmarks/                # Validation data + efficiency comparison
│   └── user_study/                # Beta survey design
├── tests/
│   ├── sample_data/               # Slim paper data (one-click load)
│   └── test_*.py                  # 293 cases, all green
└── requirements.txt
```

## Roadmap

- [x] **Phase 0** — Project skeleton + environment + knowledge base
- [x] **Phase 1** — 7 Reads-based figures + base parameter tuning + export
- [x] **Phase 2** — 5 MAG-based figures + code generator
- [x] **Phase 3** — Cycle-figure inference + hypothesis scorer + Fork Bundle + standalone interactive HTML
- [x] **v0.8 Sunday Sprint** — Beginner kit + export center + HTML v1.3
- [x] **English README + LICENSE** (this release)
- [ ] Zenodo DOI (after `v1.0` release)
- [ ] Paper Methods + Results + Discussion drafts
- [ ] Phase 4 — Plugin framework (after paper acceptance)

## Beta tester feedback wanted

EnvMeta is collecting beta feedback for the methodology paper (target: iMeta / Bioinformatics).

If you'd like to try EnvMeta and share feedback, the best paths are:

- **GitHub Issues** (preferred for bug reports / feature requests): <https://github.com/redlizzxy/EnvMeta/issues>
- **Email**: 18872605913@163.com

(There is also an internal Chinese-language beta survey available on the [Chinese README](README_CN.md).)

## Citation

EnvMeta is being prepared as a methodology paper. Until publication, please cite this repository URL. A DOI will be appended here once issued (Zenodo, planned for the `v1.0` release).

## Acknowledgments

EnvMeta originated from a metagenomic study of arsenic-slag / steel-slag microbial remediation. We thank beta-stage testers for their feedback (named individually in the acknowledgments section of the published paper).

## License

[MIT License](LICENSE) — Copyright (c) 2026 redlizzxy and EnvMeta contributors.

## Contact (email preferred)

- GitHub Issues: <https://github.com/redlizzxy/EnvMeta/issues>
- Bug reports / feature requests: 18872605913@163.com
