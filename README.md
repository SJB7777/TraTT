# DynamicFTT

Transient formation dynamics of acoustic Bragg satellites in dynamical X-ray
diffraction. See `DESIGN.md` for the physics conventions, the analytic ground
truth, and the phase plan.

## Setup

```
uv sync --extra dev
uv run pytest -q          # 12 analytic checks on the solver
```

## Layout

```
src/dftt/        solver.  core.py (Laue IVP + Bragg Riccati), drive.py (phi(zeta;t))
tests/           the analytic ground truth of DESIGN.md, as asserts
explore/         exploration notebooks.  outputs kept locally, stripped on commit
figs/            paper figures only.  fixed names fig1..fig6.  overwritten, never accumulated
data/            computed cubes + runs.csv.  gitignored, regenerate with scripts/
scripts/         only things that write to figs/ or data/
```

## File-naming rule

Figure sprawl comes from putting parameters in filenames. Don't.

- parameters live **inside** the file (`npz` meta, figure title/footer), never in its name
- `data/runs/000042.npz` + one `data/runs.csv` index maps id -> parameters
- a file may enter `figs/` only if it has a paper figure number; everything else
  stays inline in a notebook

## Notebooks

`explore/*.ipynb` are tracked, but `nbstripout` is installed as a git filter
(`.gitattributes`), so **outputs are removed at commit time while staying in your
working copy**. Nothing to remember, nothing to run by hand.

The filter is registered in `.git/config` with an absolute path to `.venv`. If
you move the repo or rebuild the venv, re-run:

```
uv run nbstripout --install --attributes .gitattributes
```

To strip a notebook manually anyway: `uv run nbstripout explore/foo.ipynb`.

Pick the `.venv` kernel in VSCode/Jupyter; `dftt` is installed editable by
`uv sync`, so `import dftt` just works.
