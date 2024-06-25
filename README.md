# GSAlign_replica
# A replica of GSAlign with changes made for simplicity. Code for API from NCBI and Kegg to extract an organism's information. 
# This information can be parsed using a code in order to run a comparison between the known genes of two species. 
# Another code is in this repository to act as a means to extract the similarity scores from the output.

# Assumptions:
# 1. Mutations occur independently of each other. A mutation in one area would not affect the probability of one occurring somewhere else, even though this happens in vivo.

# 2. Mutation rates are constant, and not influenced by biological or environmental factors.

# 3. The fitness is only dependent on the total gene expression level, and does not change as a result of a change in environment. In addition, the alleles are both neutral in terms of fitness.

# 4. There is no selection pressure on the sequences or the mutations.

# 5. Mating is random and certain phenotypes are not selected for. The population is also hermaphroditic, so biological sex does not play a role in mate selection.

# Data Structures:

# 1. Population array: NumPy array with 2 dimensions where one is the population size and the other is the alleles that the individual has (represented by either 1 or 0)

# 2. Allele frequency array: NumPy array with 2 dimensions where one is the mean frequency of allele 1 and the other is the mean frequency of allele 2 for that generation.

# 3. Expression level/fitness arrays: NumPy arrays to store expression levels and fitness over time for each generation.

# Distributions used:

# 1. Uniform distribution: in cases like the substitution mutation, where it was a random choice of ATGC for base substitutions

# 2. Bernoulli distribution: in binary outcomes like whether or not a mutation will occur (the randomly generated number is either greater than or less than the mutation rate, so the event either occurs or it doesn't)



