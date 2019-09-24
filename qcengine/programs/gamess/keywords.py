import collections
import textwrap


def format_keyword(opt, val, lop_off=True):
    """Reformat `val` for option `opt` from python into GAMESS-speak."""

    text = ''

    # Transform booleans into Fortran booleans
    if str(val) == 'True':
        text += '.true.'
    elif str(val) == 'False':
        text += '.false.'

    # No Transform
    else:
        text += str(val).lower()

    if lop_off:
        return opt[7:].lower(), text
    else:
        return opt.lower(), text


def format_keywords(options):
    """From GAMESS-directed, non-default options dictionary `options`, write a GAMESS deck."""

    grouped_options = collections.defaultdict(dict)
    for group_key, val in options.items():
        group, key = group_key.split('__')
        grouped_options[group.lower()][key.lower()] = val

    grouped_lines = {}
    for group, opts in sorted(grouped_options.items()):
        line = []
        line.append(f'${group.lower()}')
        for key, val in sorted(grouped_options[group].items()):
            line.append('='.join(format_keyword(key, val, lop_off=False)))
        line.append('$end\n')
        grouped_lines[group] = textwrap.fill(' '.join(line), initial_indent=' ', subsequent_indent='  ')

    return '\n'.join(grouped_lines.values()) + '\n'
