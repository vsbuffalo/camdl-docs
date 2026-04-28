"""Run camdl simulate and parse the result as a polars DataFrame.

Two roles ``camdl simulate`` plays in book chapters:

1. *Pedagogical demo* — the chapter wants to show the literal command
   and its output (stdout + colored stderr). Use ``run_cli`` from
   ``cli.py`` for this; the output is rendered HTML.

2. *Plumbing* — the chapter wants the simulated trajectory as a
   DataFrame to plot, summarize, etc. Use ``run_simulate`` (this
   module). The function parses stdout as TSV and returns
   ``polars.DataFrame``. Pass ``show=True`` to *also* surface the
   command + stderr in the rendered cell (useful when the simulate
   step itself is part of the lesson, not just plumbing).

Both helpers share ``_run_capture`` from ``cli.py`` so subprocess
invocation, env handling, and shell-vs-arglist behavior live in one
place.
"""
from __future__ import annotations

import io
import tempfile
from pathlib import Path
from typing import Union

import polars as pl

from .cli import _run_capture, ansi_to_html, _wrap, truncate


def run_simulate(
    model: Union[str, Path],
    params: Union[str, Path, None] = None,
    *,
    seed: int = 42,
    replicates: Union[int, None] = None,
    backend: str = "chain_binomial",
    dt: float = 1.0,
    scenario: Union[str, None] = None,
    obs: Union[bool, str, Path] = False,
    extra_args: Union[list[str], None] = None,
    cwd: Union[str, Path, None] = None,
    show: bool = False,
    show_max_lines: Union[int, None] = 30,
) -> Union[pl.DataFrame, tuple[pl.DataFrame, pl.DataFrame]]:
    """Run ``camdl simulate`` and return the parsed trajectory.

    Parameters
    ----------
    model
        Path to the ``.camdl`` model file (relative to ``cwd`` if set).
    params
        Path to a parameter TOML (relative to ``cwd`` if set). Optional —
        camdl can also resolve params from the model defaults.
    seed
        RNG seed (passed via ``--seed``).
    replicates
        If set, passes ``--replicates N`` for a multi-rep ensemble.
    backend, dt, scenario
        Standard ``camdl simulate`` flags. ``scenario`` is omitted if
        ``None``.
    obs
        - ``False`` (default): don't request an observation file; return
          only the trajectory.
        - ``True``: write observations to a temp file, parse, and return
          ``(trajectory, observations)`` as a 2-tuple.
        - path-like: write observations to that path; same return shape
          as ``True``.
    extra_args
        Extra CLI args appended after the standard ones (e.g.
        ``["--enable", "vaccination"]``).
    cwd
        Working directory for the subprocess.
    show
        If ``True``, also display the command + stderr inline in the
        rendered cell (HTML output, ANSI-styled). Useful when the
        simulate step is part of the lesson; default ``False`` keeps
        the cell quiet for plumbing-only use.
    show_max_lines
        When ``show=True``, truncate stderr to this many lines so a
        chatty simulator doesn't flood the rendered output.

    Returns
    -------
    polars.DataFrame
        Trajectory (always returned).
    (polars.DataFrame, polars.DataFrame)
        If ``obs`` is truthy: ``(trajectory, observations)``.
    """
    cmd: list[str] = ["camdl", "simulate", str(model)]
    if params is not None:
        cmd.extend(["--params", str(params)])
    cmd.extend(["--seed", str(seed),
                "--backend", backend,
                "--dt", str(dt)])
    if replicates is not None:
        cmd.extend(["--replicates", str(replicates)])
    if scenario is not None:
        cmd.extend(["--scenario", scenario])

    obs_path: Path | None = None
    cleanup_tmp = False
    if obs is True:
        tmpdir = Path(tempfile.mkdtemp(prefix="camdl_sim_"))
        obs_path = tmpdir / "obs.tsv"
        cleanup_tmp = True
    elif isinstance(obs, (str, Path)):
        obs_path = Path(obs)
    if obs_path is not None:
        cmd.extend(["--obs", str(obs_path)])

    if extra_args:
        cmd.extend(extra_args)

    result = _run_capture(cmd, cwd=str(cwd) if cwd else None, check=True)

    if show:
        from IPython.display import display
        echo = " ".join(cmd)
        body_parts: list[str] = []
        if result.stderr.strip():
            stderr_text = truncate(result.stderr.rstrip(),
                                    max_lines=show_max_lines)
            body_parts.append(
                f'<span style="opacity:0.7">{ansi_to_html(stderr_text)}</span>'
            )
        body = "\n".join(body_parts) if body_parts else "(no output)"
        display(_wrap(body, echo=echo, collapse=True))

    # Parse trajectory from stdout, skipping comment lines
    tsv_lines = [l for l in result.stdout.split("\n") if not l.startswith("#")]
    traj = pl.read_csv(io.StringIO("\n".join(tsv_lines)), separator="\t",
                       infer_schema_length=10000)

    if obs_path is not None:
        obs_df = pl.read_csv(obs_path, separator="\t", comment_prefix="#",
                              infer_schema_length=10000)
        if cleanup_tmp:
            obs_path.unlink(missing_ok=True)
            obs_path.parent.rmdir()
        return traj, obs_df

    return traj
