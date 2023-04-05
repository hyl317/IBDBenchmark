# after calling genotype prob by bcftools, set heterozygous C/T, A/G site to missing
# usage: python3 removeDeaminationSite.py input.vcf.gz | bgzip -c > output.vcf.gz

import gzip
import sys

args = sys.argv
with gzip.open(args[1]) as f:
    for line in f:
        line = line.decode('utf-8').strip()
        cols = line.split('\t')
        if line.startswith('#'):
            print(line)
        else:
            sample = cols[9]
            if sample == '.':
                print(line)
            else:
                ref, alt = cols[3], cols[4]
                if (ref == 'C' and alt == 'T') or (ref == 'T' and alt == 'C') \
                    or (ref == 'A' and alt == 'G') or (ref == 'G' and alt == 'A'):
                    if sample.split(':')[0] == '0/1':
                        ch, pos, _, ref, alt, *_ = cols
                        print('\t'.join([ch, pos, '.', ref, alt, '.', '.', '.', 'GT', '.']))
                    else:
                        print(line)
                else:
                    print(line)