# Contributing

All contributions are welcome: bug reports, suggestions, documentation and
pull requests, on [GitHub](https://github.com/koszullab/hicstuff).

## Development set-up

```bash
git clone https://github.com/koszullab/hicstuff.git
cd hicstuff
uv sync --extra dev
```

`uv run hicstuff --version` then runs the version in the repository.

## Tests, lint and format

```bash
uv run pytest tests/                 # the tests
uv run ruff check src/               # lint
uv run ruff format src/              # format (line length 100)
uv run pre-commit install            # run ruff and the file checks on every commit
```

Tests that need a program that is not installed (bowtie2, bwa, minimap2) are
skipped. New functions should come with tests in `tests/`, and the lines a
pull request changes must be at least 80% covered: Codecov checks it.

Docstrings follow the [numpydoc standard](https://numpydoc.readthedocs.io/en/latest/format.html):
they are the [Python API reference](api/hicstuff.md) of this site.

## Documentation

This site is built from the Markdown files of `docs/` by
[Zensical](https://zensical.org), and served by GitHub Pages at
<https://koszullab.github.io/hicstuff/>.

```bash
uv sync --extra docs
uv run zensical serve     # live preview at http://localhost:8000
uv run zensical build     # the site, in site/
```

- **Pages** are the Markdown files of `docs/`; `nav` in `zensical.toml` lists
  them in order. A new page goes in both.
- **The Python API reference** is written from the docstrings: each page of
  `docs/api/` holds `::: hicstuff.<module>` directives, which
  [mkdocstrings](https://mkdocstrings.github.io) expands.
- **The CLI reference**, `docs/cli.md`, is written by hand. `tests/test_docs.py`
  fails when a command or an option of the CLI is missing from it, or when it
  documents one the CLI does not have.
- **Diagrams** are [Mermaid](https://mermaid.js.org) code blocks.
- **The look** comes from [chromokit](https://github.com/js2264/chromokit),
  the UI kit of metacooler's documentation. Its files (`overrides/`,
  `docs/stylesheets/chromokit.css`, `docs/javascripts/chromokit.js`) are
  listed in `chromokit.json`: do not edit them, put hicstuff's own styles in
  `docs/stylesheets/site.css`. To take up a newer chromokit, and to check that
  no kit file was edited:

    ```bash
    npx github:js2264/chromokit sync
    npx github:js2264/chromokit check
    ```

Every push to `master` builds the site and deploys it
([docs.yml](https://github.com/koszullab/hicstuff/blob/master/.github/workflows/docs.yml));
pull requests build it, with warnings, without deploying.

## Changelog

Every change that a user could notice gets an entry under `## Unreleased` in
`docs/changelog.md`, in its section: *Breaking changes*, *New*, *Changed* or
*Fixed*.

## Releasing

Releases are cut with [bump-my-version](https://github.com/callowayproject/bump-my-version),
configured in `pyproject.toml`. From a clean, up-to-date `master`:

```bash
uv sync --extra dev                  # once, to install bump-my-version
uv run bump-my-version bump minor    # or patch / major
git push --follow-tags
```

The bump sets the version in `pyproject.toml`, `uv.lock`, the `Dockerfile`
and `zensical.toml` (the version this site shows), turns `## Unreleased` in
`docs/changelog.md` into a dated section for the version, commits as
`chore(release): vX.Y.Z` and creates the `vX.Y.Z` tag. Preview it with
`uv run bump-my-version bump patch --dry-run -vv`.

Pushing the tag starts the
[publish workflow](https://github.com/koszullab/hicstuff/blob/master/.github/workflows/publish.yml):
it checks that the tag matches the version of `pyproject.toml`, publishes the
package to PyPI and pushes the Docker image to `ghcr.io/koszullab/hicstuff`.
The push of the release commit to `master` deploys this site, with the new
version in its header.
