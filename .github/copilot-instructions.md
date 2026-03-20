# Copilot Instructions for MSAexplorer

This file defines repository-specific standards for AI-assisted edits.
Use these rules for all code, tests, and documentation changes.

## 1) Core Principles

- Keep changes small, readable, and reviewable.
- Prefer simple and explicit code over clever abstractions.
- Avoid overengineering: introduce new classes/indirection only when they clearly improve maintainability.
- Keep behavior stable unless the task explicitly requests a breaking change.
- Preserve existing architecture and naming patterns.

## 2) Repository Structure and Interaction

This repository contains two Python packages with different roles:

- `msaexplorer/` (core library)
  - Analysis and data model (`explore.py`, `_data_classes.py`)
  - Plotting (`draw.py`)
  - Export helpers (`export.py`)
  - CLI entrypoint (`cli.py`)
- `app_src/` (Shiny app frontend)
  - UI/server glue for interactive exploration (`shiny_user_interface.py`, `shiny_server.py`, `shiny_plots.py`)

How they interact:

- `app_src` should consume stable APIs from `msaexplorer`.
- `msaexplorer` must remain independently usable without the app.
- Feature work starts in `msaexplorer`; `app_src` is adapted afterwards if needed.
- Keep app-only concerns out of core modules.

## 3) Code Cleanliness Standards

- Use single quotes for normal strings: `'text'`.
- Use double quotes only when required (e.g., quote escaping readability) or for docstrings.
- Keep functions focused on one task; split long functions when logic becomes hard to scan.
- Prefer descriptive names over abbreviations.
- Remove dead code, unused imports, and outdated comments.
- Add comments only for non-obvious reasoning, not for trivial operations.
- Keep public APIs typed (input and return types) where practical.
- Reuse existing helpers/dataclasses before adding new structures.

## 4) Naming Conventions

- Functions/methods: `snake_case`
- Variables: `snake_case`
- Classes/dataclasses: `PascalCase`
- Constants: `UPPER_SNAKE_CASE`
- Internal helpers: prefix with `_` (module-private intent)
- Boolean names should read naturally (`include_ambig`, `show_legend`, `is_valid`)

## 5) Docstring Standards

Use concise, practical docstrings.

- Use triple double quotes for docstrings.
- First line: short summary sentence in imperative/present style.
- Describe parameters, return values, and raised exceptions when relevant.
- Keep docstrings short and precise; avoid tutorial-length blocks in function docstrings.
- Document every public class and public function.
- For internal/private helpers, add docstrings when behavior is not obvious.

Recommended structure:

```python
"""
Calculate pairwise distances against the current reference.

:param distance_type: Distance metric key.
:return: Pairwise distances container.
:raises ValueError: If the metric is unsupported.
"""
```

## 6) Documentation Guidelines (pdoc-oriented)

- API docs are generated primarily from docstrings (pdoc).
- Keep module/class/function docstrings accurate and synchronized with behavior.
- Prefer small executable examples over long narrative snippets.
- Any example added to docs (including package `__init__.py` examples) must be tested locally before integration.
- Do not document planned behavior as if already implemented.
- When refactoring signatures or return types, update related docstrings in the same change.

## 7) Testing Requirements

Every code change should be validated before merge.

Minimum expectations:

- Run targeted tests for edited modules.
- Add/adjust tests when behavior, output format, or API contracts change.
- Ensure export/plot/stat outputs are validated with deterministic assertions.

Typical commands:

```bash
pytest -q
```

For local source verification during development:

```bash
PYTHONPATH="/absolute/path/to/MSAexplorer" pytest -q
```

Guidance:

- Prefer explicit regression tests over broad smoke-only checks.
- If a bug is fixed, include a test that fails before and passes after.
- Keep fixtures minimal and representative.

## 8) Change Discipline

Before submitting changes:

- Confirm naming and style consistency.
- Confirm docstrings match implementation.
- Confirm tests pass.
- Confirm no unrelated edits were introduced.

If uncertain, prefer the simpler design.

