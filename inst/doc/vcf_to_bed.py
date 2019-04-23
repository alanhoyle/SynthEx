from __future__ import print_function
import re
import sys

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def get_sample_info(sample_format, format_map):
    sample_info = {}
    fields = sample_format.split(':')
    sample_info['GT'] = fields[format_map.index('GT')]
    sample_info['DP'] = fields[format_map.index('DP')]
    if 'RO' in format_map: 
        sample_info['RO'] = fields[format_map.index('RO')] # i.e. Freebayes
    if 'AD' in format_map :  # i.e. Strelka2
        sample_info['RO'] = fields[format_map.index('AD')].split(',')[0]
    return sample_info

if len(sys.argv) < 3:
    eprint ("ERROR:",sys.argv[0], "requires two numeric parameters at the command line. ")
    eprint ("USAGE: ")
    eprint ("   ",sys.argv[0], "$MIN_READ_DEPTH $MIN_ALLELE_FREQUENCY")
    sys.exit(1)

min_read_depth = float(sys.argv[1])
min_AF = float(sys.argv[2])

pass_count = 0
filter_count = 0

for line in sys.stdin:
    if line.startswith('#'):
        continue
    else:
        line = line.rstrip()
        fields = line.split('\t')
        if len(fields) > 8:

            if not (len(fields[3]) == 1 and len(fields[4]) == 1 and \
                re.search('[TCGA]', fields[3]) and re.search('[TCGA]', fields[4])):
                filter_count += 1
                continue

            format = fields[8]

            format_map = format.split(':')
            sample = fields[9]
            sample_format = get_sample_info(sample, format_map)
            
            GT = sample_format['GT']

            if GT != '0/1':
                filter_count += 1
                continue

            DP = 0
            if 'DP' in sample_format:
                DP = abs(float(sample_format['DP']))

            if DP < min_read_depth:
                filter_count += 1
                continue
            
            RO = 0
            if 'RO' in sample_format:
                RO = abs(float(sample_format['RO']))
            elif 'AD' in sample_format:
                RO =  abs(float(sample_format['AD'].split(',')[0]))
                
            RAF = 0
            if DP > 0:
                RAF = RO/DP
            
            MAF = RAF
            if RAF > 0.5:
                MAF = 1 - RAF

                
            if MAF >= min_AF:
                MAF = round(MAF, 3) 
                pass_count += 1
                # BED is 0-based exclusive, VCF is 1-based, so need to translate
                seq = (str(fields[0]).replace('chr', ''), \
                    str(int(fields[1])-1), \
                        str(int(fields[1])), \
                            str(MAF))
                output = '\t'.join(seq)
                print (output)
            else:
                filter_count += 1

eprint('Passed:',  str(pass_count))
eprint('Filtered: ', str(filter_count))
eprint('Done.')

