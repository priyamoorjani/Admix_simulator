#! /usr/bin/env python3
import argparse
import random
from shutil import copyfile

def main():
    parser = argparse.ArgumentParser(description="Replace values in a dataset with missing data into a dataset randomly.")
    parser.add_argument("-r", help="ratio of SNPs to convert to missing data.", dest="ratio", type=float, required=True)
    parser.add_argument("-f", help="file to convert, must be eigenstrat format.", dest="filename", type=str, required=True)
    parser.add_argument("-o", help="target filename for conversion.", dest="outfile", type=str, required=False)
    parser.add_argument("-i", help="designate conversion as in-place.", action="store_true", required=False)
    parser.add_argument("-a", help="model ancient DNA. create pseudo-diploid reads.", action="store_true", required=False)
    parser.add_argument("-s", help="random seed.", dest="seed", type=int, required=False)

    parser.set_defaults(func=replace_data)
    args = parser.parse_args()
    args.func(args)

def replace_data(args):
    if bool(args.outfile) == bool(args.i): # XOR
        assert False, "Must designate one and only one of outfile or in-place modification."

    if bool(args.i):
        args.outfile = args.filename

    if args.seed is not None:
        random.seed(args.seed)

    lines = []
    with open(args.filename) as f:
        for line in f:
            if line != None:
                lines.append(line.strip())

    def replace(c):
        if random.uniform(0, 1) < r:
            return '9'
        if bool(args.a) and c == '1':
            return '0' if random.uniform(0, 1) > 0.5 else '2'
        return c

    r = args.ratio
    for i in range(len(lines)):
        lines[i] = "".join(replace(c) for c in lines[i])

    with open(args.outfile, 'w') as f:
        for line in lines:
            f.write(line + "\n")


if __name__=="__main__":
    main()