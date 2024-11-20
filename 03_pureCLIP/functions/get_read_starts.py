'''
Usage:
    get_read_starts.py <in_bed> <out_bed>

Description:
    Get a bed file of the read start sites from the alignment bed file
'''


from docopt import docopt

def main(argv=None):
    
    #Parse arguments to get all input file names
    args = docopt(__doc__, argv=argv)

    f = open(args['<in_bed>'], 'r')
    output = open(args['<out_bed>'], 'w')
    for line in f:
        read = line.strip('\n').split('\t')
        strand = read[5]
        if strand == '+':
            start = int(read[1]) -1 
            end = start + 1
        elif strand == '-':
            start = int(read[2])
            end = start + 1
        if start < 0:
            start = 0
            end = 1
        output.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(read[0], str(start), str(end), read[3], read[4], read[5]))

    output.close()

if __name__ == '__main__':
    main()
