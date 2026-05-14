"""Shared CLI runner for camdl-book chapters.

Runs shell commands, captures stdout/stderr separately, converts ANSI SGR
escape sequences to styled HTML, and returns IPython HTML objects ready
for display in Quarto notebooks.

Usage in a .qmd cell::

    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
    from _lib.cli import run_cli

    run_cli("camdl simulate sir.camdl --seed 42",
            echo="camdl simulate sir.camdl --seed 42")
"""

from __future__ import annotations

import html as _html
import re
import subprocess
from typing import Literal

from IPython.display import HTML

# ── ANSI SGR → CSS ──────────────────────────────────────────────────────
# Mirrors _filters/ansi.lua so rendered HTML matches the Lua filter's
# fenced-block output.

_SGR: dict[str, str | None] = {
    # Reset
    "0":  None,
    # Bold / faint / italic / underline
    "1":  "font-weight:bold",
    "2":  "opacity:0.6",
    "3":  "font-style:italic",
    "4":  "text-decoration:underline",
    # Standard foreground
    "30": "color:#1a1a1a",
    "31": "color:#cc3333",
    "32": "color:#2e8b2e",
    "33": "color:#b8860b",
    "34": "color:#3366cc",
    "35": "color:#9b59b6",
    "36": "color:#1abc9c",
    "37": "color:#cccccc",
    # Bright foreground
    "90": "color:#666666",
    "91": "color:#e74c3c",
    "92": "color:#27ae60",
    "93": "color:#f39c12",
    "94": "color:#5dade2",
    "95": "color:#c39bd3",
    "96": "color:#48c9b0",
    "97": "color:#ffffff",
}

_ESC_RE = re.compile(r"\033\[([0-9;]*)m")


def ansi_to_html(text: str) -> str:
    """Convert ANSI SGR sequences in *text* to inline-styled HTML spans.

    Handles compound codes (e.g. ``\\033[1;31m`` → bold red) and resets.
    All literal ``<``, ``>``, ``&`` are escaped first so the result is
    safe for embedding in ``<pre><code>`` blocks.
    """
    escaped = _html.escape(text, quote=False)
    result: list[str] = []
    open_spans = 0
    pos = 0

    for m in _ESC_RE.finditer(escaped):
        # Emit text before this escape
        if m.start() > pos:
            result.append(escaped[pos:m.start()])

        params = m.group(1)
        if params in ("", "0"):
            # Reset: close all open spans
            result.append("</span>" * open_spans)
            open_spans = 0
        else:
            styles = []
            for code in params.split(";"):
                s = _SGR.get(code)
                if s:
                    styles.append(s)
            if styles:
                result.append(f'<span style="{";".join(styles)}">')
                open_spans += 1

        pos = m.end()

    # Remainder after last escape
    result.append(escaped[pos:])
    result.append("</span>" * open_spans)
    return "".join(result)


def truncate(
    text: str,
    *,
    max_lines: int | None = None,
    head: int | None = None,
    tail: int | None = None,
) -> str:
    """Truncate *text* by line count.

    Three modes (first match wins):

    ``head`` only
        Keep the first *head* lines.
    ``tail`` only
        Keep the last *tail* lines.
    ``max_lines`` (with optional ``tail``)
        If the text exceeds *max_lines*, show the first
        ``max_lines - tail`` lines, an "[…N lines…]" marker,
        then the last *tail* lines (default 6).
    """
    lines = text.split("\n")

    if head is not None and tail is None and max_lines is None:
        if len(lines) <= head:
            return text
        return "\n".join(lines[:head] + [f"  [...{len(lines) - head} more lines...]"])

    if tail is not None and head is None and max_lines is None:
        if len(lines) <= tail:
            return text
        return "\n".join([f"  [...{len(lines) - tail} lines above...]"] + lines[-tail:])

    if max_lines is not None:
        _tail = tail if tail is not None else 6
        if len(lines) <= max_lines:
            return text
        n_head = max_lines - _tail
        omitted = len(lines) - n_head - _tail
        return "\n".join(
            lines[:n_head]
            + [f"  [...{omitted} lines omitted...]"]
            + lines[-_tail:]
        )

    return text


Stream = Literal["stdout", "stderr", "both"]


def _run_capture(
    cmd: str | list[str],
    *,
    cwd: str | None = None,
    env: dict[str, str] | None = None,
    check: bool = False,
) -> subprocess.CompletedProcess:
    """Thin wrapper over ``subprocess.run`` that captures both streams as text.

    Used by both ``run_cli`` (which renders the result as HTML) and
    ``run_simulate`` (which parses stdout as TSV). Centralizes the
    shell vs. arglist + env-merge logic.
    """
    import os as _os
    shell = isinstance(cmd, str)
    run_env = None
    if env:
        run_env = {**_os.environ, **env}
    result = subprocess.run(
        cmd, capture_output=True, text=True, shell=shell, cwd=cwd, env=run_env,
    )
    if check and result.returncode != 0:
        import sys as _sys
        cmd_str = cmd if isinstance(cmd, str) else " ".join(map(str, cmd))
        print(f"command failed (exit {result.returncode}): {cmd_str}", file=_sys.stderr)
        if result.stdout:
            print("--- stdout ---", file=_sys.stderr)
            print(result.stdout, file=_sys.stderr)
        if result.stderr:
            print("--- stderr ---", file=_sys.stderr)
            print(result.stderr, file=_sys.stderr)
        raise subprocess.CalledProcessError(
            result.returncode, cmd, result.stdout, result.stderr,
        )
    return result


def fit_clean(*configs: "str | Path", root: "str | Path" = "results/fits") -> None:
    """Remove stale fit dirs for one or more fit configs.

    Each ``config`` can be a path to a fit toml (``fits/foo.toml``) or
    a bare stem (``foo``). The function deletes every
    ``<root>/<stem>-*`` directory it finds.

    Use case: pinned at the top of a chapter's setup cell so each
    render starts with at most one fit dir per stem. After
    ``fit_clean(...)`` runs, ``Path(root).glob('<stem>-*').next()``
    is unambiguous (one match) for every fit the chapter performs.

    No-op if a stem has no matching dirs. Order-independent. Cheap —
    typical chapter rendering pays a few hundred microseconds.
    """
    import shutil
    from pathlib import Path
    root_p = Path(root)
    for cfg in configs:
        stem = Path(cfg).stem if str(cfg).endswith(".toml") else str(cfg)
        for d in root_p.glob(f"{stem}-*"):
            if d.is_dir():
                shutil.rmtree(d)


def run_cli(
    cmd: str | list[str],
    *,
    echo: str | None = None,
    cwd: str | None = None,
    collapse: bool = True,
    show: Stream = "both",
    max_lines: int | None = 40,
    head: int | None = None,
    tail: int | None = None,
    check: bool = False,
    env: dict[str, str] | None = None,
    stderr_style: str = "opacity:0.7",
    font_size: str | None = None,
) -> HTML:
    """Run a shell command and return styled HTML output.

    Parameters
    ----------
    cmd
        Command string (passed to shell) or arg list.
    echo
        Command text shown as the prompt / summary line.  If *collapse*
        is True, it becomes a ``<details><summary>`` toggle.
    cwd
        Working directory for the subprocess.
    collapse
        Wrap output in a collapsible ``<details>`` element (requires *echo*).
    show
        Which stream(s) to display: ``"stdout"``, ``"stderr"``, or
        ``"both"`` (default).  When ``"both"``, stderr lines are wrapped
        in a dim span (see *stderr_style*) so the two streams are
        visually distinguishable.
    max_lines
        If set, truncate with head + tail (see :func:`truncate`).
        Set to ``None`` to disable truncation.
    head
        Show only the first *head* lines (overrides *max_lines*).
    tail
        Show only the last *tail* lines, or the tail portion of
        *max_lines* truncation.
    check
        Raise ``subprocess.CalledProcessError`` on non-zero exit.
    env
        Extra environment variables (merged with ``os.environ``).
    stderr_style
        CSS applied to stderr lines when ``show="both"``.
    font_size
        Optional CSS font-size applied to the output ``<pre>`` (e.g.
        ``"0.75em"``, ``"11px"``). Default ``None`` leaves theme styling
        untouched. Useful for wide flag listings or tabular output.

    Returns
    -------
    IPython.display.HTML
        Ready to ``display()`` in a Quarto/Jupyter cell.
    """
    result = _run_capture(cmd, cwd=cwd, env=env, check=check)

    # Assemble the raw text according to `show`
    if show == "stdout":
        raw = result.stdout.rstrip()
    elif show == "stderr":
        raw = result.stderr.rstrip()
    else:
        # "both": interleave with stderr visually dimmed
        parts: list[str] = []
        if result.stdout.rstrip():
            parts.append(result.stdout.rstrip())
        if result.stderr.rstrip():
            stderr_html = ansi_to_html(
                truncate(result.stderr.rstrip(),
                         max_lines=max_lines, head=head, tail=tail)
            )
            # For "both", we build the combined HTML manually so we can
            # style stderr differently.
            stdout_text = result.stdout.rstrip()
            stdout_html = ansi_to_html(
                truncate(stdout_text, max_lines=max_lines, head=head, tail=tail)
            ) if stdout_text else ""

            body_parts: list[str] = []
            if stdout_html:
                body_parts.append(stdout_html)
            if stderr_html:
                body_parts.append(
                    f'<span style="{stderr_style}">{stderr_html}</span>'
                )
            body = "\n".join(body_parts)
            return _wrap(body, echo=echo, collapse=collapse, font_size=font_size)

        raw = "\n".join(parts) if parts else ""

    # Single-stream path
    raw = truncate(raw, max_lines=max_lines, head=head, tail=tail)
    body = ansi_to_html(raw)
    return _wrap(body, echo=echo, collapse=collapse, font_size=font_size)


def _wrap(
    body: str,
    *,
    echo: str | None,
    collapse: bool,
    font_size: str | None = None,
) -> HTML:
    """Wrap converted HTML body in the appropriate container."""
    style_attr = f' style="font-size:{font_size}"' if font_size else ""
    if collapse and echo:
        return HTML(
            f'<details class="cli-run">'
            f'<summary><code>$ {_html.escape(echo)}</code></summary>'
            f'<pre class="ansi-output"{style_attr}><code>{body}</code></pre>'
            f'</details>'
        )
    prompt = f'<span style="opacity:0.6">$ {_html.escape(echo)}</span>\n' if echo else ""
    return HTML(f'<pre class="ansi-output"{style_attr}><code>{prompt}{body}</code></pre>')
