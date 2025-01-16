from django import template
from common.definitions import G_PROTEIN_DISPLAY_NAME as g_prot_dict
from common.definitions import ARRESTIN_DISPLAY_NAME as arr_dict

import re

register = template.Library()

@register.filter
def join_attr(obj_list, attr_name, sep=', '):
    return sep.join(getattr(i, attr_name) for i in obj_list)

# @register.filter
# def lineformat ( objs ):
#     elements = [obj.name for obj in objs]
#     if len(elements) > 0:
#         return ", ".join(elements)
#     else:
#         return '-'

@register.filter
def endo_format ( objs ):
    elements = list(set([obj.ligand.name for obj in objs]))
    if len(elements) > 0:
        return ", ".join(elements)
    else:
        return '-'

@register.filter
def create_struct_links ( objs ):
    return ", ".join(["<a href=\"/structure/" + obj.pdb_code.index + "\">" + obj.pdb_code.index + "</a>" for obj in objs])

@register.filter
def dashwhenempty (obj):
    if obj:
        return obj
    else:
        return '-'

# .replace('<sup>','').replace('</sup>','').replace('<sub>','').replace('</sub>','').replace('&alpha;','alpha').replace('&beta;','beta')

@register.filter
def ligandrole ( objs ):
    objs = sorted(objs, key=lambda obj: obj.ligand_role.name)
    elements = [obj.ligand_role.name for obj in objs]
    if len(elements) > 0:
        if len(elements) > 1:
            if elements[0]==elements[1]:
                return elements[0]
            else:
                return "\n".join(elements)
        else:
            return "\n".join(elements)
    else:
        return '-'

@register.filter
def ligandtype ( objs ):
    objs = sorted(objs, key=lambda obj: obj.ligand_role.name)
    elements = [obj.ligand.ligand_type.name for obj in objs]
    if len(elements) > 0:
        return "\n".join(elements)
    else:
        return '-'

@register.filter
def only_gproteins ( objs ):
    elements = [element for obj in objs for element in obj.name.split(',') if re.match(".*G.*", element) and not re.match(".*thase.*|PGS", element)]
    if len(elements) > 0:
        return "\n".join(elements)
    else:
        return '-'

@register.filter
def alpha_display_name ( obj ):
    if obj.family.name != "GPa1":
        return '&alpha;' + obj.family.name[3:].lower()
    else:
        return obj.family.name

@register.filter
def only_one_subunit ( objs, arg ):
    protfams, value = arg.split(',')
    if '-' in protfams:
        protfams = protfams.split('-')
    else:
        protfams = [protfams]
    if "Alpha" in protfams or "Arrestin" in protfams:
        elements = [element for element in objs if element.category in ["G alpha", "Arrestin"]]
    else:
        elements = [element for element in objs if element.wt_protein.family.parent.name in protfams]
    if len(elements) > 0:
        if value=='name':
            if "Alpha" in protfams and elements[0].display_name[0]=='G':
                return '&alpha;'+elements[0].display_name[1:]
            elif elements[0].display_name[0]=='G':
                return elements[0].display_name[1:]
            elif "Arrestin" in protfams:
                return elements[0].wt_protein.name
            elif "Beta" in protfams:
                return elements[0].display_name.split(" ")[0]
        elif value=='id':
            return elements[0].id
        elif value=='entry_name':
            return elements[0].wt_protein.entry_name
        elif value=='species':
            return elements[0].wt_protein.species.common_name
        elif value=='note':
            if elements[0].note:
                return elements[0].note
            else:
                return '-'
        elif value=='family':
            return elements[0].wt_protein.family.parent.name
        elif value=='coverage':
            return elements[0].wt_coverage
        else:
            return '-'
    else:
        return '-'

@register.filter
def only_arrestins ( objs ):
    elements = [element for obj in objs for element in obj.name.split(',') if re.match(".*rrestin.*", element)]
    if len(elements) > 0:
        return "\n".join(elements)
    else:
        return '-'

@register.filter
def only_fusions ( objs ):
    elements = [element for obj in objs for element in obj.name.split(',') if re.match(".*thase.*|PGS|BRIL|.*Lysozyme|.*b562.*|TrxA|Flavodoxin|Rubredoxin|Sialidase|.*Thioredoxin.*|Endolysin|.*cytochrome.*|.*DARPin.*", element)] #not re.match(".*bod.*|.*Ab.*|.*Sign.*|.*G.*|.*restin.*|.*scFv.*|.*Fab.*|.*activity.*|.*RAMP.*|.*peptide.*|.*CD4.*", element) or
    if len(elements) > 0:
        return "\n".join(elements)
    else:
        return '-'

@register.filter
def only_antibodies ( objs ):
    elements = [element for obj in objs for element in obj.name.split(',') if re.match(".*bod.*|.*Ab.*|.*scFv.*|.*Fab.*|.*activity.*|.*RAMP.*|.*GRK.*|Unidentified peptide|.*CD4.*|.*IgG.*|.*NB.*|.*Fv.*", element)]
    if len(elements) > 0:
        return "\n".join(elements)
    else:
        return '-'

@register.filter
def only_other_proteins ( objs ):
    elements = [element for obj in objs for element in obj.name.split(',') if not re.match(".*thase.*|PGS|BRIL|.*Lysozyme|.*b562.*|TrxA|Flavodoxin|Rubredoxin|Sialidase|.*Thioredoxin.*|Endolysin|.*cytochrome.*|.*bod.*|.*Ab.*|.*scFv.*|.*Fab.*|.*activity.*|.*RAMP.*|.*CD4.*|.*IgG.*|.*NB.*|.*Fv.*|Go|Gi1|G11|Gs|Gt|G alpha|Gi2|.*G protein.*", element)]
    if len(elements) > 0:
        return "\n".join(elements)
    else:
        return '-'

@register.filter
def last_author ( objs ):
    if objs:
        return objs.split(',')[-1]
    else:
        return '-'

@register.filter
def cut_at_20 ( objs ):
    return objs[:20]

@register.filter
def get_refined_model_version ( objs ):
    return objs.pdb_data.pdb.split('\n')[0][-10:]

@register.filter
def cut_refined ( objs ):
    return objs.split('_')[0]

@register.filter
def cut_classname ( objs ):
    if objs == "Other GPCRs":
        return "Other"
    else:
        return objs[5:]

@register.filter
def entry_short ( objs ):
    return objs.split("_")[0].upper()

@register.filter
def gprot_short ( objs ):
    return g_prot_dict[objs]

@register.filter
def arr_short ( objs ):
    return arr_dict[objs]

@register.filter
def receptor_short ( objs ):
    if not objs.startswith('mGlu'):
        objs = objs[0].upper()+objs[1:]
    return objs.replace(" receptor","").replace("-adrenoceptor","")
