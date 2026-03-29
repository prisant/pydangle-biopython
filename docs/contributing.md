# Contributing

Thanks for considering contributing to this project!  Here's how to get set up
and what to keep in mind.

## Development setup

```bash
git clone https://github.com/prisant/pydangle-biopython.git
cd pydangle-biopython

python -m venv .venv
source .venv/bin/activate

make install
```

This installs the package in editable mode with dev dependencies and
configures the git changelog hook.

## Making changes

1. Create a branch from `main`:

    ```bash
    git checkout -b feature/my-change
    ```

2. Make your changes in the `src/pydangle_biopython/` directory.

3. Add or update tests in `tests/`.  Use the shared fixtures from `conftest.py`
   when possible.

4. Run the checks:

    ```bash
    make all    # lint + type check + tests
    ```

5. Commit with a descriptive message:

    ```bash
    git commit -m "feat: add DNA-specific builtins"
    ```

## Commit message conventions

The changelog generator filters commits by prefix.  Use these conventions:

| Prefix     | Purpose                | In changelog? |
|------------|------------------------|---------------|
| `feat:`    | New feature            | Yes           |
| `fix:`     | Bug fix                | Yes           |
| `docs:`    | Documentation update   | Yes           |
| `test:`    | Test additions/changes | Yes           |
| `refactor:`| Code restructuring     | Yes           |
| `chore:`   | Maintenance tasks      | No            |
| `wip:`     | Work in progress       | No            |

## Code style

This project uses:

- **ruff** for linting and formatting (configured in `pyproject.toml`)
- **mypy** with strict mode for type checking
- **pytest** for testing

Run them individually or all at once:

```bash
make lint        # ruff check
make format      # ruff format
make typecheck   # mypy
make test        # pytest
make all         # all of the above (except format)
```

## Project structure

When adding new functionality, follow the existing module pattern:

- **`builtins.py`** — Built-in measurement command definitions
- **`parser.py`** — Command string parsing and translation
- **`measure.py`** — Geometry computation and structure traversal
- **`main.py`** — CLI entry point (keep thin, delegate to measure)
- **`exceptions.py`** — Custom exception classes

## Adding dependencies

Add runtime dependencies to `[project.dependencies]` and dev-only tools
to `[project.optional-dependencies.dev]` in `pyproject.toml`, then run:

```bash
pip install -e ".[dev]"
```

## Documentation

Docs live in `docs/` and are built with MkDocs:

```bash
make docs        # build and serve locally at http://127.0.0.1:8000
```

Update `mkdocs.yml` if you add new pages.

## Releasing

```bash
# Tag the release
git tag -a v0.2.0 -m "Description of release"

# Push the tag
git push origin v0.2.0

# Trigger changelog grouping
git commit --allow-empty -m "chore: start v0.3.0 development"
```
