# -*- coding: utf-8 -*-
"""
astk.cli.config
~~~~~~~~~~~~~~~~~
This module provides command line api configure.
"""
import json
from pathlib import Path
from gettext import gettext as _

import click
from click import UsageError
from click.core import _check_iter
from click.exceptions import BadParameter

from astk.utils import RunConfigure


class CustomMultiCommand(click.Group):
    """https://stackoverflow.com/questions/46641928/python-click-multiple-command-names/53144555"""
    ALIASES  = {}
    def command(self, *args, **kwargs):
        """Behaves the same as `click.Group.command()` except if passed
        a list of names, all after the first will be aliases for the first.
        """
        
        def decorator(f):
            if args and args[0]:
                _args = [args[0][0]] + list(args[1:])
                for alias in args[0][1:]:
                    self.ALIASES[args[0][0]] = args[0][1]
                    cmd = super(CustomMultiCommand, self).command(
                        alias, *args[1:], **kwargs)(f)
                    cmd.short_help = "Alias for '{}'".format(_args[0])
            else:
                _args = args
            from click.decorators import command
            cmd = command(*_args, **kwargs)(f)
            self.add_command(cmd)
            return cmd
        return decorator

    def add_command(self, cmd, name=None) -> None:
        """Registers another :class:`Command` with this group.  If the name
        is not provided, the name of the command is used.
        """
        name = name or [cmd.name]
        if name is None:
            raise TypeError("Command has no name.")
        # _check_multicommand(self, name, cmd, register=True)
        for n in name:
            self.commands[n] = cmd

    def resolve_command(self, ctx, args):
        from click.utils import make_str
        from click.parser import split_opt
        cmd_name = make_str(args[0])
        original_cmd_name = cmd_name

        cmd_name = self.ALIASES.get(cmd_name, cmd_name)
        # Get the command
        cmd = self.get_command(ctx, cmd_name)

        # If we can't find the command but there is a normalization
        # function available, we try with that one.
        if cmd is None and ctx.token_normalize_func is not None:
            cmd_name = ctx.token_normalize_func(cmd_name)
            cmd = self.get_command(ctx, cmd_name)

        # If we don't find the command we want to show an error message
        # to the user that it was not provided.  However, there is
        # something else we should do: if the first argument looks like
        # an option we want to kick off parsing again for arguments to
        # resolve things like --help which now should go to the main
        # place.
        if cmd is None and not ctx.resilient_parsing:
            if split_opt(cmd_name)[0]:
                self.parse_args(ctx, ctx.args)
            ctx.fail(("No such command {name!r}.").format(name=original_cmd_name))
        return cmd_name if cmd else None, cmd, args[1:]


class MultiOption(click.Option):
    """Thanks https://stackoverflow.com/questions/48391777/nargs-equivalent-for-options-in-click"""

    def __init__(self, *args, **kwargs):
        self.save_other_options = kwargs.pop('save_other_options', True)
        nargs = kwargs.pop('nargs', -1)
        assert nargs == -1, 'nargs, if set, must be -1 not {}'.format(nargs)
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


@click.command(help="ASTK configure setting")
@click.option('-R', '--R', "RPath", type=click.Path(exists=True), help="R path setting")
def sc_setting(*args, **kwargs):

    rc = RunConfigure()
    rc.update(Rscript=str(Path(kwargs['RPath']).with_name("Rscript")))
    rc.save()
