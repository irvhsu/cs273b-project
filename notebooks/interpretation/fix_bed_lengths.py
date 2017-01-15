def main(in_bed):
    for line in in_bed:
        chrom, start, end, a, b, strand = line.strip().split()[:6]
        start = str(int(start) - 145 + 31)
        end   = str(int(end)   + 145)
        print '\t'.join([chrom, start, end, a, b, strand])
        
if __name__ == '__main__':
    import sys
    with open(sys.argv[1], 'r') as f:
        main(f)
