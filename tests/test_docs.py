"""The CLI reference of the documentation (docs/cli.md) follows the CLI.

docs/cli.md is written by hand: one section per command, headed
``### `command` { #command }``, with a table whose first column lists each
option in backticks, short form first (``-g, --genome``). These tests fail
when a command or an option of the CLI is missing from it, or when it
documents one the CLI does not have.
"""

import re
from pathlib import Path

import click
import pytest

from hicstuff.cli import cli

CLI_MD = Path(__file__).resolve().parents[1] / "docs" / "cli.md"

#: The heading of a command's section: ### `pipeline` { #pipeline }
_SECTION = re.compile(r"^### `([a-z0-9_-]+)`", re.MULTILINE)
#: The first cell of an options table row: | `-g, --genome` | ...
_OPTION_CELL = re.compile(r"^\| `(-[^`]+)` \|", re.MULTILINE)


def documented_options():
    """The options documented for each command, by command name."""
    text = CLI_MD.read_text(encoding="utf-8")
    headings = list(_SECTION.finditer(text))
    sections = {}
    for heading, following in zip(headings, headings[1:] + [None]):
        body = text[heading.end() : following.start() if following else len(text)]
        options = set()
        for cell in _OPTION_CELL.findall(body):
            options.update(name.strip() for name in cell.split(","))
        sections[heading.group(1)] = options
    return sections


def cli_options(command):
    """The options of a click command, without --help, nor the hidden ones,
    which the help does not show either: a former name, such as --matfmt."""
    options = set()
    for param in command.params:
        if isinstance(param, click.Option) and not param.hidden:
            options.update(param.opts + param.secondary_opts)
    return options - {"-h", "--help"}


def test_every_command_is_documented():
    ctx = click.Context(cli)
    assert set(documented_options()) == set(cli.list_commands(ctx))


@pytest.mark.parametrize("name", sorted(cli.list_commands(click.Context(cli))))
def test_every_option_is_documented(name):
    command = cli.get_command(click.Context(cli), name)
    documented = documented_options()[name]
    missing = cli_options(command) - documented
    unknown = documented - cli_options(command)
    assert not missing, f"options of `hicstuff {name}` missing from docs/cli.md: {missing}"
    assert not unknown, f"docs/cli.md documents options `hicstuff {name}` lacks: {unknown}"
