author = "Cove Soyars"

#This program takes in a csv of miRNA expression where each column is a case of either skin cancer or AML and each row
#represents an miRNA. It finds, using a t-test, miRNAs that are differentially expressed between the cancers. This code
#is largely inspired by the permutation test code from class. It incorporates the methods for importing the data and
#retrieving the indicies for a given condition (here it is a type of cancer). It also uses similar methods to the t test
#KNN filtering code from class

import csv
from scipy import stats

########################################### DEFINE FUNCTIONS ##########################################################

def read_csv(expressionfile):
    """this function retrieves the miRNA expression from the CSV"""
    table = [] #define list to hold rows, making a 2D list
    with open(expressionfile) as file: #open file
        next(file) #skip header
        reader = csv.reader(file, quoting=csv.QUOTE_NONNUMERIC, #define csv reader
                                quotechar='"')
        for row in reader: #for row in reader,
            if sum(row) != 0: #if row is not all zeros,
                table.append(row) #add row to table
    return table

def miRNA_ID(ID_file):
    """This function reads the ID file, containing the IDs of the miRNAs"""
    IDs = [] #define list to hold IDs
    with open(ID_file) as file:
        reader = csv.reader(file) #define reader
        for ID in reader:
            IDs.append(ID[0]) #append each ID to the list, indexed at 0 because each ID is in a one-item list
    return IDs

def read_headers(expressionfile):
    """This function gets the headers from the CSV"""
    with open(expressionfile) as file: #open file
        headers = next(file).strip().split(',') #headers = first row of file, stripped of formatting characters and split
        #by commas
    return headers

def get_cancer_type(type, header):
    """This function gets the indices associated with the given cancer type from a header"""
    indices = [] #define list to hold indices
    for i, column in enumerate(header): #loop through header, saving index
        if type == column: #if given type of cancer is equal to that column in the header,
            indices.append(i) #append the index
    return indices

def get_differentially_expressed(expression, names, type1, type2):
    """This function returns the miRNAs that have a t-test p value > 0.001"""
    #calculate t-test pvalue for each miRNA, and return it's ID if it passes the p value test
    differentially_expressed = [] #define list to hold differentially expressed miRNAs
    for i in range(0, len(expression)):
        miRNA = expression[i]  # make variable for gene to make list comprehension easier in next line
        pval = stats.ttest_ind([miRNA[ind] for ind in type1], [miRNA[ind] for ind in type2]).pvalue #calculate p value
        if pval <= 0.05/len(expression): #if p value passes test with Bonferroni correction, (see comments below function)
            differentially_expressed.append(names[i]) #add the miRNA's name to the list
    return differentially_expressed

#The before using the Bonferroni correction, using an alpha of 0.05 resulted in 475 miRNAs passing the threshold. This
#is more than the alpha value can support, which is 72.1 (1442 * 0.05). Using the Bonferronni correction's alpha
# (0.05/1442) results in 68 miRNAs passing the threshold, which is much closer to 72.1 and more believable.

########################################## MAIN ######################################################################

expression = read_csv("miRNA_project_data.csv") #read in expression

ID = miRNA_ID('soyars_project_data_rownames.csv') #read in miRNA IDs

header = read_headers("soyars_project_data.csv") #read in header

##define AML and Skin cancer indices
AML = get_cancer_type('AML',header)
Skin = get_cancer_type('Skin', header)

diff = get_differentially_expressed(expression, ID, AML, Skin) #find differentially expressed miRNAs

#Print results
print "The amount of differentially expressed miRNAs is %s" % (len(diff))

#write results to file for use in mirPath:
out_fileH = open("differentially_expressed_miRNAs.txt", 'w')
for miRNA in diff:
    out_fileH.write(miRNA)
    out_fileH.write('\n')
out_fileH.close()

