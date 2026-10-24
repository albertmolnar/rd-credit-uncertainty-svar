# Credit Supply, Uncertainty and U.S. Private R&D Investment

Replication files for:

> Albert Molnár, "Credit Supply, Uncertainty and U.S. Private R&D Investment: Evidence from a Sign-Restricted VAR."
> Keleti Károly Faculty of Business and Management, Óbuda University.

The code estimates a Bayesian structural VAR for the United States (1986Q1–2019Q4). It identifies credit supply, uncertainty and aggregate supply shocks with sign and relative-magnitude restrictions. It then estimates lag-augmented local projections, constructs the R&D capital stock following the BEA R&D Satellite Account (Sliker, 2007), and produces every table and figure in the paper.

## Requirements

- MATLAB R2026a (earlier releases from R2020a should work). No toolboxes are required.
- Run time is about 6–10 minutes on a laptop.
- The random seed is fixed (`cfg.seed` in `code/lib/default_config.m`), so the same MATLAB release reproduces the posterior draws exactly.

## How to run

In MATLAB, change to the `code` folder and run:

```matlab
main                 % data, VAR, identification, local projections, robustness, coverage check
make_paper_outputs   % compact figures and tables used in the paper
export_to_paper      % copies tables/figures into paper/ and builds paper_v2_overleaf.zip
report_numbers       % prints the numbers quoted in the text
```

Optional:

```matlab
rd_capital_analysis  % R&D capital stock: benchmark regression and sensitivity to the depreciation rate
check_deflator       % baseline with nominal instead of deflated R&D (Section III-B)
```

Outputs are written to `code/output/` (tables as `.tex`, figures as `.pdf`, a log file).

To compile the paper, run `export_to_paper`, then compile `paper/paper_v2.tex` with pdfLaTeX and BibTeX, or upload `paper_v2_overleaf.zip` to Overleaf.

## Data

All series are public. The quarterly dataset is built from the raw files in `code/data_raw/`. These are the vintage used in the paper, downloaded on 16 October 2026. Daily and monthly series are averaged to quarters.

| File | Series | Source | Code |
|---|---|---|---|
| `RD_NOM.csv` | Private fixed investment, IPP: research and development (nominal, SAAR) | BEA via FRED | Y006RC1Q027SBEA |
| `GDP.csv` | Gross domestic product (nominal, SAAR) | BEA via FRED | GDP |
| `GDPDEF.csv` | GDP implicit price deflator | BEA via FRED | GDPDEF |
| `GS1.csv` | 1-year Treasury constant maturity yield (daily) | Federal Reserve H.15 via FRED | DGS1 |
| `EPU.csv` | Economic Policy Uncertainty index for the US (daily) | Baker, Bloom and Davis via FRED | USEPUINDXD |
| `EBP.csv` | Excess bond premium (monthly) | Federal Reserve Board, [ebp_csv.csv](https://www.federalreserve.gov/econres/notes/feds-notes/ebp_csv.csv) | – |
| *not included* | Moody's Baa yield minus 10-year Treasury (daily) | Moody's via FRED | BAA10Y |
| *not included* | CBOE Volatility Index (daily) | Cboe via FRED | VIXCLS |

Real GDP is GDP/GDPDEF×100. Real R&D investment is nominal R&D deflated by the GDP deflator, because an R&D price index was not available to the author.

**BAA10Y and VIXCLS are copyrighted** by Moody's and Cboe and are therefore not redistributed here. They are used only in two robustness checks, Table VI rows "Baa–10y spread for EBP" and "VIX for EPU". If the files are missing, `main.m` tries to download them from FRED. If the download fails, those two rows are skipped and everything else is unaffected. You can also download them manually from FRED as CSV or XLSX and save them in `code/data_local/` as `BAA10Y.xlsx` and `VIXCLS.xlsx`.

## Repository structure

```
code/
  main.m                 full pipeline
  make_paper_outputs.m   compact figures/tables for the paper
  export_to_paper.m      copies outputs into paper/ and builds the Overleaf zip
  report_numbers.m       numbers quoted in the text
  rd_capital_analysis.m  R&D capital stock (Sliker, 2007)
  check_deflator.m       nominal-R&D robustness check
  lp_coverage_check.m    Monte Carlo coverage of the local-projection intervals (called by main.m)
  lib/                   functions (VAR estimation, posterior sampling, sign restrictions,
                         local projections, unit-root tests, perpetual inventory, ...)
  data_raw/              raw data (vintage used in the paper)
paper/
  paper_v2.tex, references.bib, IEEEtran.cls
```

## Licence

Code: MIT licence (see `LICENSE`). Data: subject to the terms of the original providers listed above.
