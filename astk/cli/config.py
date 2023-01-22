# -*- coding: utf-8 -*-
"""
astk.cli.config
~~~~~~~~~~~~~~~~~
This module provides command line api configure.
"""
from pathlib import Path
from gettext import gettext as _

import click
from click.core import _check_iter
from click.exceptions import BadParameter, UsageError
import configparser

from astk.utils import RunConfigure
from astk.constant import BASE_DIR


parser = configparser.RawConfigParser()
parser.read([BASE_DIR / "cli/aliases.ini"])
ALIASES_DIC = dict(parser.items("aliases"))


class MultiOption(click.Option):
    """Thanks https://stackoverflow.com/questions/48391777/nargs-equivalent-for-options-in-click"""

    def __init__(self, *args, **kwargs):
        self.save_other_options = kwargs.pop('save_other_options', True)
        nargs = kwargs.pop('nargs', -1)
        assert nargs == -1, f'nargs, if set, must be -1 not {nargs}'
        # kwargs["nargs"] = -1
        super(MultiOption, self).__init__(*args, **kwargs)
        self._previous_parser_process = None
        self._eat_all_parser = None

    def add_to_parser(self, parser, ctx):
        
        def parser_process(value, state):
            # method to hook to the parser.process
            done = False
            value = [value]
            if self.save_other_options:
                # grab everything up to the next option
                while state.rargs and not done:
                    for prefix in self._eat_all_parser.prefixes:
                        if state.rargs[0].startswith(prefix):
                            done = True
                    if not done:
                        value.append(state.rargs.pop(0))
            else:
                # grab everything remaining
                value += state.rargs
                state.rargs[:] = []
            value = tuple(value)

            # call the actual process
            self._previous_parser_process(value, state)

        retval = super(MultiOption, self).add_to_parser(parser, ctx)
        for name in self.opts:
            our_parser = parser._long_opt.get(name) or parser._short_opt.get(name)
            if our_parser:
                self._eat_all_parser = our_parser
                self._previous_parser_process = our_parser.process
                our_parser.process = parser_process
                break
        return retval

    def type_cast_value(self, ctx,  value):
        """Convert and validate a value against the option's
        :attr:`type`, :attr:`multiple`, and :attr:`nargs`.
        """
        if value is None:
            return () if self.multiple or self.nargs == -1 else None

        def check_iter(value):
            try:
                return _check_iter(value)
            except TypeError:
                # This should only happen when passing in args manually,
                # the parser should construct an iterable when parsing
                # the command line.
                raise BadParameter(
                    _("Value must be an iterable."), ctx=ctx, param=self
                ) from None

        def convert(value):
            value = tuple(check_iter(value))
            return tuple(self.type(x, self, ctx) for x in value)
        return convert(value)


class AliasedGroup(click.Group):
    """refer to https://click.palletsprojects.com/en/8.1.x/advanced/"""
    def get_command(self, ctx, cmd_name):
        rv = click.Group.get_command(self, ctx, cmd_name)
        if rv is not None:
            return rv
        if cmd_name in ALIASES_DIC:
            actual_cmd = ALIASES_DIC[cmd_name]
            return click.Group.get_command(self, ctx, actual_cmd)
        else:            
            matches = [x for x in self.list_commands(ctx)
                    if x.startswith(cmd_name)]
            if not matches:
                return None
            elif len(matches) == 1:
                return click.Group.get_command(self, ctx, matches[0])
            ctx.fail(f"Too many matches: {', '.join(sorted(matches))}")

    def resolve_command(self, ctx, args):
        # always return the full command name
        _, cmd, args = super().resolve_command(ctx, args)
        return cmd.name, cmd, args


@click.group(
        cls=AliasedGroup, 
        context_settings={
            "help_option_names": ['-h', '--help'],
            "max_content_width": 120})
def cli_fun():
     """
     Welcome to use ASTK!\f
     """


@cli_fun.command(help="ASTK configure setting")
@click.option('-R', '--R', "RPath", type=click.Path(exists=True), 
               required=True, help="R path setting")
def config(*args, **kwargs):

    rc = RunConfigure()
    rc.update(Rscript=str(Path(kwargs['RPath']).with_name("Rscript")))
    rc.save()
