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
    elif isinstance(val, dict):
        text = []
        for k, v in val.items():
            merge = [k]
            merge.extend(str(v) if isinstance(v, (int, float)) else list(map(str,v)))
            text.append(' '.join(merge))
        text = ' '.join(text)
    else:
        text = str(val)

    if lop_off:
        return opt[7:].lower(), text
    else:
        return opt.lower(), text


def format_keywords(options):
    """From NWCHEM-directed, non-default options dictionary `options`, write a NWCHEM deck."""

    def rec_dd():
        return collections.defaultdict(rec_dd)
    grouped_options = rec_dd()

    for group_key, val in options.items():
        nesting = group_key.split('__')
        if len(nesting) == 1:
            key = nesting[0]
            grouped_options['aaaglobal'][key] = val
        elif len(nesting) == 2:
            g1, key = nesting
            grouped_options[g1][key] = val
        elif len(nesting) == 3:
            g1, g2, key = nesting
            grouped_options[g1][g2][key] = val
        else:
            print(nesting)
            raise ValueError('Nesting N!')

    grouped_lines = {}
    for group, opts in sorted(grouped_options.items()):
        lines = []
        group_level_lines = []
        for key, val in grouped_options[group].items():
            if isinstance(val, dict):
                g2_level_lines = []
                g2_level_lines.append(key.lower())
                for k2, v2 in val.items():
                    line2 = ' '.join(format_keyword(k2, v2, lop_off=False))
                    g2_level_lines.append(line2)
                g2_level_lines = ' '.join(g2_level_lines)
                lines.append(g2_level_lines)
            else:
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
