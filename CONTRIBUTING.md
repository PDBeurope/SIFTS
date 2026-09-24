# Contributing to pdbe_sifts

Thank you for your interest in contributing to **pdbe_sifts**, the open-source SIFTS pipeline from EMBL-EBI PDBe!
This document covers everything you need to go from zero to an accepted pull request.

---

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [Ways to Contribute](#ways-to-contribute)
- [Setting Up Your Development Environment](#setting-up-your-development-environment)
- [Making Changes](#making-changes)
- [Code Style and Linting](#code-style-and-linting)
- [Running Tests](#running-tests)
- [Building the Documentation](#building-the-documentation)
- [Submitting a Pull Request](#submitting-a-pull-request)
- [Reporting Bugs](#reporting-bugs)
- [Suggesting Features](#suggesting-features)

---

## Code of Conduct

This project follows the [EMBL-EBI code of conduct](https://www.ebi.ac.uk/about/our-culture/).
Please be respectful and constructive in all interactions.

---

## Ways to Contribute

You do not need to write code to contribute meaningfully. We welcome:

- **Bug reports** -- something broken? Open an issue.
- **Documentation** -- typos, unclear steps, missing examples.
- **Tests** -- additional unit or integration tests are always valuable.
- **Feature requests** -- open an issue to discuss before implementing.
- **Code fixes and features** -- bug fixes, performance improvements, new capabilities.

If you are unsure where to start, look for open issues or check the areas flagged in the README.

---

## Setting Up Your Development Environment

### Prerequisites

The following tools must be installed and available on your `PATH` before you begin:

| Tool | Purpose | Install |
|---|---|---|
| [MMseqs2](https://github.com/soedinglab/MMseqs2) | Fast sequence search | `conda install -c conda-forge mmseqs2` |
| [FASTA36](https://fasta.bioch.virginia.edu/wrpearson/fasta/) (`lalign36`) | Local pairwise alignment | `conda install -c bioconda fasta3` |
| [BLAST+](https://blast.ncbi.nlm.nih.gov/) | Optional alternative to MMseqs2 | `conda install -c bioconda blast` |

**Python ≥ 3.10** and **16 GB RAM minimum** (32 GB+ recommended for large datasets) are also required.

### Option A -- uv (recommended for development)

```bash
# 1. Fork the repository on GitHub, then clone your fork
git clone https://github.com/<your-username>/SIFTS.git
cd SIFTS

# 2. Add the upstream remote so you can stay in sync
git remote add upstream https://github.com/PDBeurope/SIFTS.git

# 3. Create a virtual environment and install all dependencies
uv sync

# 4. Activate the environment
source .venv/bin/activate        # macOS / Linux
# .venv\Scripts\activate         # Windows

# 5. Install the package in editable mode
pip install -e ".[dev]"

# 6. Install pre-commit hooks
pre-commit install
```

### Option B -- micromamba

```bash
micromamba env create -f environment.yml
micromamba activate pdbe_sifts
pip install -e ".[dev]"
pre-commit install
```

---

## Making Changes

1. **Sync with upstream** before starting any new work:

   ```bash
   git fetch upstream
   git checkout master
   git merge upstream/master
   ```

2. **Create a feature branch** from `master`. Use a short, descriptive name:

   ```bash
   git checkout -b fix/short-description       # for bug fixes
   git checkout -b feat/short-description      # for new features
   git checkout -b docs/short-description      # for documentation only
   git checkout -b test/short-description      # for tests only
   ```

3. **Make your changes**, then stage and commit:

   ```bash
   git add <changed-files>
   git commit -m "fix: brief description of what changed"
   ```

   We use [Conventional Commits](https://www.conventionalcommits.org/) style loosely:
   `fix:`, `feat:`, `docs:`, `test:`, `refactor:`, `chore:`.

4. **Push** your branch to your fork:

   ```bash
   git push origin fix/short-description
   ```

---

## Code Style and Linting

This project uses [Ruff](https://docs.astral.sh/ruff/) for both linting and formatting
(configured in `pyproject.toml`). Pre-commit hooks run Ruff automatically on every commit.

To run manually:

```bash
# Lint (and auto-fix where possible)
ruff check --fix .

# Format
ruff format .
```

Key style rules (enforced by Ruff):

- Line length: **80 characters**
- Import ordering: **isort-compatible** (`I` rules)
- Target Python version: **3.10+** idioms (`UP` rules)
- No unused imports or undefined names (`F` rules)

If pre-commit blocks your commit, run the commands above, then `git add` the
auto-fixed files and commit again.

---

## Running Tests

```bash
# Run the full test suite
pytest

# Run a specific test file
pytest tests/test_<module>.py

# Run with verbose output
pytest -v

# Run with coverage report
pytest --cov=pdbe_sifts --cov-report=term-missing
```

Tests live under `tests/`. Please add or update tests when your change affects
existing behaviour or introduces new functionality.

---

## Building the Documentation

Docs are built with [MkDocs](https://www.mkdocs.org/) and
[mkdocstrings](https://mkdocstrings.github.io/):

```bash
# Install docs dependencies
pip install -e ".[docs]"

# Serve locally (auto-reloads on changes)
mkdocs serve

# Build static site
mkdocs build
```

Docstrings in the source code feed directly into the API reference pages, so
improving them also improves the published documentation.

---

## Submitting a Pull Request

1. Ensure all pre-commit hooks pass (`pre-commit run --all-files`).
2. Ensure all tests pass (`pytest`).
3. Open a pull request from your branch to `PDBeurope/SIFTS:master` on GitHub.
4. Fill in the PR description with:
   - **What** the change does.
   - **Why** it is needed (link to a related issue if applicable: `Closes #<number>`).
   - Any **testing** you did beyond the automated suite.
5. A maintainer will review your PR. Please respond to review comments promptly.

Keep PRs **focused** -- one logical change per PR is much easier to review
and merge than a large mixed-purpose change.

---

## Reporting Bugs

Open an issue at <https://github.com/PDBeurope/SIFTS/issues> and include:

- Your operating system and Python version.
- The exact command or code you ran.
- The full error message / traceback.
- The version of `pdbe_sifts` (run `pdbe_sifts --version` or check `pyproject.toml`).

---

## Suggesting Features

Open an issue describing:

- The use case or problem you are trying to solve.
- Your proposed solution, if you have one in mind.
- Any alternatives you considered.

Discussing a feature in an issue first avoids wasted effort and helps the
maintainers give early feedback before you invest time in implementation.

---

## Questions?

For general usage questions, open a GitHub issue with the `question` label.
For anything sensitive, contact the PDBe team directly via
[pdbehelp@ebi.ac.uk](mailto:pdbehelp@ebi.ac.uk).
