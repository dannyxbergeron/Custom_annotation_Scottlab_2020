import sys

in_file = sys.argv[1]
out_file = sys.argv[2]

chr = [str(x) for x in range(1,23)] + ['X', 'Y', 'MT']

with open(in_file, 'r') as in_f:
    with open(out_file, 'w') as out_f:
        for line in in_f.read().splitlines():
            fields = line.split('\t')
            if fields[0] in chr or line.startswith('#!'):
                out_f.write(line)
                out_f.write('\n')
