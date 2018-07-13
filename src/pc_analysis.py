#! /usr/bin/env python3
import argparse
import matplotlib.pyplot as plt
from collections import defaultdict

def main():
    parser = argparse.ArgumentParser(description="Perform a naive estimation " +
        "of admixture ratio by projecting the admixed population and the " +
        "reference populations onto their first two principal components and " +
        "do a simple variance analysis to determine the estimated ratio and " +
        "standard error.")
    parser.add_argument("-f", help="file to convert, list of eigenvectors " +
        "from smartpca.", dest="filename", type=str, required=True)
    parser.add_argument("-p", help="populations to be processing. Format is:" +
        "'Admixed;Ref1,Ref2'", dest="populations", type=str, required=True)
    parser.add_argument("-o", help="if designated, output a plot of the \
        individuals to the given filename", dest="outfile", type=str, 
        required=False)

    parser.set_defaults(func=pc_analysis)
    args = parser.parse_args()
    args.func(args)




def pc_analysis(args):
    with open(args.filename) as file:
        line = file.readline()
        eig1, eig2 = line.split()[1], line.split()[2]
        vals = defaultdict(lambda: [])
        for line in file:
            parts = line.split()
            val1, val2, indiv = parts[1], parts[2], parts[-1]
            if not vals[indiv]:
                vals[indiv] = []
            vals[indiv].append([val1, val2])

    pops = args.populations.split(';')
    admix = pops[0]
    reference = pops[1].split(',')

    if (args.outfile):
        plt.figure()
        plt.xlabel("Principal Component 1, eigenvalue: " + str(eig1))
        plt.ylabel("Principal Component 2, eigenvalue: " + str(eig2))
        plt.title("Simulation data compared against ancestor populations in PC")
        for key in vals:
            x1, x2 = zip(*vals[key])
            colors = {reference[0]: 'r', reference[1]: 'b', admix: 'g'}
            plt.scatter(x1, x2, c=colors[key], label=key)

        plt.legend()
        plt.savefig(args.outfile)


    x1, x2 = zip(*vals[admix])
    x1 = [float(x) for x in x1]
    x2 = [float(x) for x in x2]
    mean_x1 = sum(x1) / len(x1)
    mean_x2 = sum(x2) / len(x2)
    var_x1 = sum([(x_i - mean_x1)**2 for x_i in x1]) / len(x1)
    var_x2 = sum([(x_i - mean_x2)**2 for x_i in x2]) / len(x2)
    print("Estimated mean and deviation of the simulated data " +
        "on the first principal component:")
    print(mean_x1, "+/-", var_x1 ** 0.5)
    print()
    print("Estimated mean and deviation of the simulated data " + 
        "on the second principal component:")
    print(mean_x2, "+/-", var_x2 ** 0.5)
    print()

    refA1, refA2 = zip(*vals[reference[0]])
    refB1, refB2 = zip(*vals[reference[1]])
    mean_f = sum([float(val) for val in refA1]) / len(refA1)
    mean_y = sum([float(val) for val in refB1]) / len(refB1)
    var_f = sum([(float(x_i) - mean_f)**2 for x_i in refA1]) / len(refA1)
    var_y = sum([(float(x_i) - mean_y)**2 for x_i in refB1]) / len(refB1)


    # Naive estimate of ratio
    ratio = (mean_y - mean_x1) / (mean_y - mean_f)

    # Compute empirical covariances
    cov_x1_y = 0
    for x in x1:
        for y in refB1:
            cov_x1_y += (x - mean_x1) * (float(y) - mean_y)
    cov_x1_y /= (len(x1) * len(refB1))

    cov_x1_f = 0
    for x in x1:
        for f in refA1:
            cov_x1_f += (x - mean_x1) * (float(f) - mean_f)
    cov_x1_f /= (len(x1) * len(refA1))

    cov_f_y = 0
    for f in refA1:
        for y in refB1:
            cov_f_y += (float(f) - mean_f) * (float(y) - mean_y)
    cov_f_y /= (len(refA1) * len(refB1))


    # covariance of numerator and denom
    var_diff = var_y + var_x1 - 2 * cov_x1_y
    var_merged = var_y + var_f - 2 * cov_f_y
    # var a / b = a^2 / b^2 ( var(a) / a^2 + var(b) / b^2 - 2 cov(ab) / ab)
    # cov(y-x, y-f) = var y - cov yx - cov yf + cov fx
    true_var = ratio**2 * ( ( var_diff / (mean_y - mean_x1)**2 ) + 
                           ( var_merged / ( mean_y - mean_f)**2 ) - 
                           2 * ( ( var_y - cov_x1_y - cov_f_y + cov_x1_f ) / ( mean_y * ( mean_y - mean_f ) ) ) )

    print("Naive estimation of admixture ratio based on principal component analysis:")
    print(ratio, "+/-", true_var ** 0.5)

if __name__=="__main__":
    main()
