from common.alignment import Alignment as GenericAlignment
from protein.models import ProteinSegment

from copy import deepcopy


class Alignment(GenericAlignment):
    """A class representing a protein sequence alignment, with or without a reference sequence"""
    ECD_segments = ProteinSegment.objects.filter(name__startswith='ECD').values_list('slug', flat=True)
    sequence_based_segments = ProteinSegment.objects.filter(proteinfamily__in=['Alpha','Arrestin']).values_list('slug', flat=True)

    def merge_generic_numbers(self):
        """Check whether there are many display numbers for each position, and merge them if there are"""
        generic_numbers = deepcopy(self.generic_numbers) # deepcopy is required because the dictionary changes during the loop
        merged_num = False
        for ns, segments in generic_numbers.items():
            for segment, positions in segments.items():
                for pos, dns in positions.items():
                    if not dns:
                        self.generic_numbers[ns][segment][pos] = ""
                    elif len(dns) == 1:
                        self.generic_numbers[ns][segment][pos] = self.format_generic_number(dns[0])
                    else:
                        helix_num = False
                        seq_based_nums = []
                        struct_based_num = False
                        for dn in dns:
                            # split and format the numbers
                            split_dn = dn.split('x')
                            split_dn2 = split_dn[0].split('.')
                            
                            helix_num = split_dn2[0]
                            seq_based_num = split_dn2[1]
                            struct_based_num = split_dn[1]
                            split_seq_dn = split_dn[0].split('.')

                            if not seq_based_nums:
                                seq_based_nums.append(seq_based_num)
                            elif len(seq_based_nums) == 1 and seq_based_num not in seq_based_nums:
                                seq_based_nums.append(seq_based_num)
                            elif len(seq_based_nums) == 2 and seq_based_num not in seq_based_nums:
                                seq_based_nums[1] = seq_based_num
                        
                        # make sure the sequence-based numbers are in the correct order
                        seq_based_nums.sort()

                        # create the merged number, e.g. 1.34-35x35
                        merged_num = str(helix_num) + '.' + ('-').join(seq_based_nums) + 'x' + str(struct_based_num)
                        self.generic_numbers[ns][segment][pos] = self.format_generic_number(merged_num)

    def format_generic_number(self, generic_number):
        """Format a generic number for display in alignment viewer"""
        if not generic_number:
            return ""

        if len(generic_number.split('.'))==3:
            formatted_gn = generic_number
            split_gn = generic_number.split('.')
            subdomain = split_gn[0]
            segment = split_gn[1]
            full_segment = '.'.join([subdomain, segment])
            if full_segment in self.sequence_based_segments:# or segment in self.sequence_based_segments: # Second part is temp fix for arrestins
                seq_gn = split_gn[2]
                if seq_gn.startswith('0'):
                    seq_gn = seq_gn[1:]
                seq_gn = '.'+seq_gn
                struct_gn = '-'
            else:
                seq_gn = '-'
                struct_gn = split_gn[2]
                if struct_gn.startswith('0'):
                    struct_gn = struct_gn[1:] 
                struct_gn = '.'+struct_gn

            formatted_gn = ''
            formatted_gn += '<b>{}</b><br /><b>{}</b><br /><span class="ali-td-generic-num-normal">{}<br /></span>{}'.format(subdomain, segment, seq_gn, struct_gn)
        # elif len(generic_number.split('.'))==2:
        #     formatted_gn = '<br />'
        #     seq_class = 'ali-td-generic-num-normal'
        #     formatted_gn += '<span class="{:s}">{:s}<br /></span>'.format(seq_class, generic_number)
        # ECD format
        elif generic_number.split('x')[0] in self.ECD_segments:
            formatted_gn = '<br />'
            seq_class = 'ali-td-generic-num-normal'
            formatted_gn += '<span class="{:s}">{:s}<br /></span>'.format(seq_class, generic_number)
        elif generic_number.startswith('B.GPS'):
            formatted_gn = generic_number
            split_gn = generic_number.split('.')
            subdomain = split_gn[0]
            segment = 'GPS'
            full_segment = '.'.join([subdomain, segment])
            seq_gn = '-'
            struct_gn =  split_gn[1][-2:]
            formatted_gn = '<b>{}</b><br /><b>{}</b><br /><span class="ali-td-generic-num-normal">{}<br /></span>{}'.format(subdomain, segment, seq_gn, struct_gn)
        else:
            split_gn = generic_number.split("x")
            split_gn_helix = split_gn[0].split('.')
            helix_num = split_gn_helix[0]
            seq_num = split_gn_helix[1]
            str_num = split_gn[1]

            # is's ugly to have HTML here, but adding it in the template would awkward

            # subdomain placeholder
            formatted_gn = '-<br />'

            # helix number
            formatted_gn += '<b>{:s}</b><br />'.format(helix_num)
            
            # sequence-based number
            seq_class = 'ali-td-generic-num-normal'
            # smaller font size if there are more than one sequence-based numbers
            if len(seq_num) > 2:
                seq_class = 'ali-td-generic-num-tiny'
            # color number red if it differs from the structure-based number
            if seq_num != str_num:
                seq_class += ' ali-td-generic-num-red'
            formatted_gn += '<span class="{:s}">.{:s}<br /></span>'.format(seq_class, seq_num)
            
            # structure-based number
            formatted_gn += 'x{:s}'.format(str_num)

        return formatted_gn