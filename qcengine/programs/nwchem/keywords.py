import collections


def format_keyword(opt, val, lop_off=True):
    """Function to reformat value *val* for option *opt* from python into nwchem-speak."""

    # Transform string booleans into " "
    if val is True:
        return opt.lower(), ''
    elif val is False:
        return '', ''

    # complete hack
    if opt.upper() == 'MEMORY':
        return opt.lower(), f'{val} byte'

    elif isinstance(val, list):
        text = ' '.join([str(v) for v in val])
    else:
        text = str(val)

    if lop_off:
        return opt[7:].lower(), text
    else:
        return opt.lower(), text


def format_keywords(options):
    """From NWCHEM-directed, non-default options dictionary `options`, write a NWCHEM deck."""

    grouped_options = collections.defaultdict(dict)
    for group_key, val in options.items():
        nesting = group_key.split('__')
        if len(nesting) == 1:
            group, key = 'aaaglobal', nesting[0]
        elif len(nesting) == 2:
            group, key = nesting
        else:
            print(nesting)
            raise ValueError('Nesting!')
        grouped_options[group][key] = val

    grouped_lines = {}
    for group, opts in sorted(grouped_options.items()):
        lines = []
        group_level_lines = []
        for key, val in grouped_options[group].items():
            line = ' '.join(format_keyword(key, val, lop_off=False))
            if group.lower() == 'basis' and any(
                [word in line for word in ['spherical', 'cartesian', 'print', 'noprint', 'rel']]):
                group_level_lines.append(line)
            else:
                lines.append(line)
        if group == 'aaaglobal':
            grouped_lines[group] = '\n'.join(lines) + '\n'
        else:
            grouped_lines[group] = f'{group.lower()} ' + ' '.join(group_level_lines) + '\n  ' + '\n  '.join(
                lines) + '\nend\n'

    return '\n'.join(grouped_lines.values()) + '\n'
