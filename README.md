# Admixture simulator

This simulator can be used for generating admixed genomes. As input the method takes phased data from two populations and generates admixed individuals for a given time of admixture and proportion of ancestry from each population. Method assumes instantaneous admixture but to generate data for multiple pulses or continuous admixture, one can run the same code in a loop. Details of the method can be found in  <a href="https://priyamoorjani.files.wordpress.com/2018/06/2011_moorjani_african-gene-flow-in-levant_plosgenetics.pdf">Moorjani et al. 2011</a>.  

#### Command line: 
```
./simulation.py -p <parfilename> 
```
#### Input:
This program requires two sets of phased individuals in the EIGENSTRAT format (See https://reich.hms.harvard.edu/software/InputFileFormats). The input phased geno, snp, and ind files must be consistent (same number of SNPs in .snp and .geno file, and there should be two rows for each individual which is correctly indentified in the ind file). For an example, see example/data/CEU.ph*
 
The program also requires a parameter file. See format below.
#### Parameter file arguments:
```
ancestorAgeno:  .phgeno file for the first ancestral population
ancestorAsnp:   .phsnp file for the first ancestral population
ancestorAind:   .phind file for the first ancestral population
ancestorBgeno:  .phgeno file for the second ancestral population
ancestorBsnp:   .phsnp file for the second ancestral population
ancestorBind:   .phind file for the second ancestral population
outputfilename: The filename that outputs will be written to. 
lambda:     	The age of the admixture. Any number greater than 0.
theta:      	The admixture ratio - theta should be a float between 0 and 1, and indicates the proportion of ancestorA in the produced individuals.
n:      	The number of strands to output. Thus, if the output is haploid, the simulator will create n strands, otherwise the simulator will create 2n strands and merge them.
haploidoutput:  If 'False', output normal diploid genotypes, otherwise output phased individuals.
trackancestry:  If 'False', do not track ancestry and output information to <output>.ancestry. Each line of the output file will have the ancestry at a SNP - individuals separated by hyphens ('-').
```

#### Output:
The simulator outputs to <output>.geno, <output>.ind, <output>.snp under normal conditions, and to <output>.phgeno, <output>.phind, <output>.phsnp when haploidoutput is set to True. Additionally, the utility will write to <output>.ancestry if trackancestry is set to True. The output files will be in EIGENSTRAT format. 
The simulated individuals are labeled: NA<number>, and have gender Unknown. The population name of the outputted individuals is always "Simulation". 

#### Requirements:
The number of strands that the simulator creates must be lower than the number of individuals in the smaller of the two parental pools, as we require that at any point in the set of simulated individuals, no two draw from the same haplotype (to reduce inbreeding-like effects). Thus, at some point, it's possible that every strand being simulated draws from one pool, so if there are more strands than individuals in a pool, the simulator cannot continue. Note also that we guarantee that every crossover event includes a derangement - individuals cannot draw from one parent, crossover, and end up drawing from the same individual that they were copying from originally. At the ends of chromosomes, we redraw ancestry completely from random, essentially a crossover event that does not guarantee derangement.

#### Example:
Look in the directory example/ for full details. The command will be:
```
./simulation.py -p parfiles/simulation.par
```

# Simulating missing data

This utility allows one to set a number of genotypes (based on the user defined proportion) as missing data. This is useful for testing the impact of missing data on the inference (particularly helpful in ancient DNA studies). This is implemented by checking the genotype at every SNP for each individual, and choosing to replace the position with 9 with given probability based on the input missing fraction (-r). The analysis can be performed in two modes: 1) Diploid mode: output is left as diploid for non-missing data; 2) Pseudo-haploid mode: output is haploid or missing. This aims to mimic ancient DNA where its hard to call heterozygous sites with limited coverage. Hence typically one randomly samples a read mapping to that site and uses a pseudo-haploid call for the inference. Similarly, in our simulation we randomly sample one allele at each heterozygous site.

#### Command line: 
```
./missingdata.py -r <threshold> -f <input> -o <output> -a -s <seed>
```
#### Input: 
Inputs must be in EIGENSTRAT format, PACKEDANCESTRYMAP is not currently supported (See https://reich.hms.harvard.edu/software/InputFileFormats).

#### Arguments (REQUIRED):
```
    -r <float> : Designates the ratio of data to remove. With probability r, each character in the .geno file will be converted to a '9'
    -f <filename> : Name of the geno file to read and convert.
    -o <outfilename> XOR -i : The name outfilename designates the file to which to write the new output. -i indicates that the modification is meant to be done in-place - identical to having the output filename be the same as the input file name. You may only use one of these two options. 
```

#### Arguments (OPTIONAL):
```
    -s <seed> : Random seed for generation.
    -a : Indicates ancient data conversion. For characters that aren't chosen to be removed by the main part of the utility, if the character is a '1', indicating a heterozygous site, the utility will replace the character with either a '0' or a '2', with even probability of the outcomes.
```

#### Output:
The output file as indicated by the arguments will contain the characters from the input file, with each character being replaced with a '9' with probability r, and with '1's being replaced by a '0' or '2' with probability 0.5 each. 

#### Example:
In the directory example/, run:
```
./missingdata.py -r 0.2 -f family.geno -o family_missingdata.geno -a -s 42
```

This will use data in family.geno (output from simulator/realdata) and output a new pseudo-haploid geno file with missing-data set based on the proportion -r. 

#### Requirements: This utility only uses standard python modules. Developed in Python 3.5.2. 

# Naive PCA-based estimator for admixture proportions

This utility takes the eigenstrat output file from smartpca (https://github.com/DReichLab/EIG/tree/master/EIGENSTRAT) and performs a naive analysis to infer the admixture proportion in simulated individuals. Assuming a two-way admixture, we can perform a PCA between the admixed individuals and the two ancestral populations and then use the projection on PC 1 and 2 to infer the admixture proportion. Assume $m_1, m_2, m_3$ is the mean position on PC1 for ancestral pop1, ancestral pop2 and the admixed pop respectively. We estimate $theta$ as $abs(m_1-m_3)/abs(m_1-m_2). We recognize that this approximation does not hold under complex mixture or when true ancestrals are not available. Though this serves as a handy tool to infer admixture proportions in simulations where the truth is known. We estimate standard errors by using a simple empirical variance calculation.

#### Command line: 
```
./pc_analysis.py -f <input evec file> -p <population list>
```

#### Input:
The input must be the output style of smartpca (*.evec). The first column will be the individual IDs, the second will be the projection of the individual onto the first PC, the third will be the projection of the individual onto the second PC, and the final will be the population title. The first line will contain the eigenvalues (for display purposes only).

#### Arguments:
```
    -f <filename> : This designates the input file to the program, which is itself the output from a smartpca run on the individuals. 
    -p <poplist> : This indicates the populations to read out. The format for this input is "<Admixed population>;<Reference population 1>,<Reference population 2>". 
```

#### Output:
This utility will output to the standard output - with information on the population means and standard errors on the first two principal components, as well as the estimation of admixture ratio based on positions in the first two principal components. 

#### Example:
```
./pc_analysis.py -f "family.evec" -p "Simulation;French,Yoruba"
```

If you want to store the outputs:
```
./pc_analysis.py -f "family.evec" -p "Simulation;French,Yoruba" > pc_analysis.log
```

#### Experimental:
```
    -o <outfilename> : This will make the program try to use matplotlib to plot the individuals on the first two PCs. This is currently experimental, as whether or not this will work depends on the backend that your system's matplotlib uses.
```
#### Requirements:
The package requires Python 3.5.2 and standard packages. For plotting, matplotlib is also required. 

# Support:
Send queries to Neel Alex (salexucb@berkeley.edu) and Priya Moorjani (moorjani@berkeley.edu)
