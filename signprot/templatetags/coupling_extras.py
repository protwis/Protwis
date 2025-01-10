from django import template


register = template.Library()

@register.filter
def split(objs, arg):
    delimiter, index = arg.split(',')
    return objs.split(delimiter)[int(index)]

@register.filter
def join_list(objs):
    return '<br>'.join(sorted(objs))

@register.filter
def guidetopharma_reword(objs):
    return objs.replace('GuideToPharma', 'Guide to Pharmacology').replace('-agonist','')

@register.filter
def dict_check(objs, arg):
    if arg in objs:
        if objs[arg]!=None:
            return objs[arg]
        else:
            return '-'
    else:
        return False

@register.filter
def add_quote(objs):
    if objs!='-':
        return "{}'".format(objs)
    else:
        return objs

@register.filter
def zero_reformat(objs):
    if objs in [0.0,0.00,0.000]:
        return 0
    else:
        return objs

@register.filter
def replace(objs, arg):
    if len(arg)==1:
        return objs.replace(arg, '')
    else:
        replace_from = arg[0]
        replace_to = arg[1]
        return objs.replace(replace_from, replace_to)

@register.filter
def pub_out(objs):
    # if objs.authors:
    main_auth = objs.authors.split(',')[0]
    return '{} et al. <i>{}</i> {}'.format(main_auth, objs.journal.name, objs.year)
    # else:
    #     return 'weird {}'.format(objs)

@register.filter
def dict_get(dictionary, key):
    try:
        return dictionary.get(key, '')
    except (TypeError, AttributeError):
        return ''

@register.filter
def remove_prefix(value):
    prefixes = ['Psites_', 'Osites_', 'Plong_', 'Olong_']
    for prefix in prefixes:
        if value.startswith(prefix):
            return value[len(prefix):]
    return value

@register.filter
def join_if_list(value, separator=', '):
    if isinstance(value, list):
        return separator.join(map(str, value))
    return value

@register.filter
def compute_width_min(label, min_width=3):
    length = len(label)
    if length < min_width:
        length = min_width

    # Calculate width: length * character width factor + padding
    # Adjust the character width factor and padding as needed
    width = length * 1 + 1.5  # Character width factor and padding

    # Return width as a string with 'em' unit
    return f'{width}em'

@register.filter
def range_filter(n):
    return range(n)
