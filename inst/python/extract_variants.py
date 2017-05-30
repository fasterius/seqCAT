#!/usr/bin/env python3

# Import modules
import argparse
import vcf
import os.path

# Argument parser
parser = argparse.ArgumentParser(epilog='Extracts variant data from a VCF.')
parser.add_argument('input', type=str, help='input VCF file path')
parser.add_argument('sample', type=str, help='sample to extract')
parser.add_argument('output', type=str, help='output file path')
args = parser.parse_args()

# Sample to extract from VCF
sample = args.sample

# Remove output file if already existing
if os.path.isfile(args.output):
    os.remove(args.output)

# Open output file for appending
output_file = open(args.output, 'a')

# Open input VCF file
vcf_reader = vcf.Reader(open(args.input, 'r'))

for record in vcf_reader:
    ref=record.REF
    alt=record.ALT
    id=record.ID
    chrom=record.CHROM
    pos=record.POS

    # Collect genotype info (skip record if any is missing)
    try:
        gt = record.genotype(str(sample))['GT']
        ad1 = record.genotype(str(sample))['AD']
        dp = record.genotype(str(sample))['DP']

    except AttributeError:
        continue

    try:
        ann=record.INFO['ANN']
    except KeyError:
        continue

    qual=record.QUAL
    filt1=record.FILTER

    if filt1:
         filt=filt1[0]
    else:
         filt='None'


    if gt:

        gts = gt.split('/')
        A1 = gts[0]
        A2 = gts[1]

        if A1=='0':
            A1GT=ref
        else:
            if A1=='1':
                A1GT=alt[0]
            else:
                if A1=='2':
                    A1GT=alt[1]
                else:
                    A1GT='NONE'

        if A2=='0':
           A2GT=ref
        else:
            if A2=='1':
                A2GT=alt[0]
            else:
                if A2=='2':
                    A2GT=alt[1]
                else:
                    A2GT='NONE'

        result = ['%s %s %s %s %s %s %s %s %s %s' %
            (line.split('|')[3], line.split('|')[4], line.split('|')[2],
             line.split('|')[1], line.split('|')[5], line.split('|')[6],
             line.split('|')[7], line.split('|')[9], line.split('|')[10],
             line.split('|')[15]) for line in ann]

        result = set(result)
        priority_list = list(reversed(['HIGH', 'MODERATE', 'LOW', 'MODIFIER']))
        max_index = -1

        for line in result:
            gene, ensgid, impact, effect, feature, enst, biotype, nucl, aacid, warnings = line.split(' ')
            max_index = max(max_index, priority_list.index(impact))
        max_impact = priority_list[max_index]

        for line in result:
            gene, ensgid, impact, effect, feature, enst, biotype, nucl, aacid, warnings = line.split(' ')
            if impact == max_impact:
                current = str(chrom) + "\t" + \
                          str(pos) + "\t" + \
                          str(id) + "\t" +  \
                          str(ref) + "\t" + \
                          str(alt) + "\t" + \
                          str(gene) + "\t" + \
                          str(enst) + "\t" + \
                          str(ensgid) + "\t" + \
                          str(impact) + "\t" + \
                          str(effect) + "\t" + \
                          str(feature) + "\t" + \
                          str(biotype) + "\t" + \
                          str(A1GT) + "\t" + \
                          str(A2GT) + "\t" + \
                          str(dp) + "\t" + \
                          str(filt) + "\t" + \
                          str(ad1) + "\t" + \
                          str(nucl) + "\t" + \
                          str(aacid) + "\t" + \
                          str(warnings) + "\n"

                # Write to file
                output_file.write(current)

# Close output file
output_file.close()
