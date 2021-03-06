import csv as csv
import sys
import getopt

def read_biogrid(filename):
    """
    Reads a BioGRID TAB2 formatted file - extracts the protein-protein
    interactions only. Each interaction is stored twice - as (a,b) and (b,a)
    Returns a set of tuples representing the interactions
    """
    interactions = set()
    with open(filename, "rU") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            if r['Experimental System Type'] == 'physical':
                if (r['Official Symbol Interactor A'] and 
                    r['Official Symbol Interactor B']):
                    interactions.add((r['Official Symbol Interactor A'],
                                      r['Official Symbol Interactor B']))
                    interactions.add((r['Official Symbol Interactor B'],
                                      r['Official Symbol Interactor A']))
    return interactions


def read_interactions(filename, directed=False):
    """
    Reads a tab delimited file of interactions - format is 'GeneA\tGeneB\n'
    and no header is required.

    directed indicates whether the interaction should be considered
    directed or not - e.g. for kinase-substrate interactions. If not directed
    interactions are stored as both (GeneA,GeneB) and (GeneB,GeneA)

    Returns a set of tuples representing the interactions
    """
    interactions = set()
    with open(filename, "rU") as f:
        reader = csv.reader(f, delimiter="\t")
        for r in reader:
            a = r[0]
            b = r[1]
            interactions.add((a, b))
            if not directed:
                interactions.add((b, a))
    return interactions


def annotate_dependencies(
        dependencies_file, output_file, interactions, rowname):
    """
    Annotates a set of associations generated by James Campbell's R Scripts
    according to whether they occur between gene pairs known to interact
    
    Args:
        dependencies_file - associations file from the R script
        output_file - name of the file to output annotated dependencies
        interactions - set of tuples representing known interactions
        rowname - column name to use for the annotations 
    """
    with open(dependencies_file, "rU") as fin, open(output_file, "w") as fout:
        reader = csv.DictReader(fin, delimiter="\t")
        fieldnames = reader.fieldnames + [rowname]
        writer = csv.DictWriter(fout, fieldnames, delimiter="\t")
        writer.writeheader()
        for row in reader:
            geneA = row['marker'].split('_')[0]
            geneB = row['target'].split('_')[0]
            row[rowname] = (geneA, geneB) in interactions
            writer.writerow(row)
    return

# Below is the executable part of the program

help_message = '''
This is a simple script to annotated genetic dependencies according to known
functional relationships. This reads an associations file and outputs it with
an additional column indicating whether each dependency involves a gene pair
with a known functional relationship.

Required parameters are as follows:
-a --associations (filename): the associations file to use (from James Campbells R functions)
-o --output (filename): the file to output to
-i --interactions (filename): the interactions file to use

Optional parameters as follows:
-b --biogrid: indicates that interaction file is in Biogrid Tab2 format, defaults to False
-d --directed: treat interactions as directed (e.g kinase-substrate), defaults to False
-n --name (string): column name for the annotation, defaults to "FunctionalInteraction"
-h --help: displays this message

'''


class Usage(Exception):

    def __init__(self, msg):
        self.msg = msg


def main(argv=None):
    biogrid = False
    directed = False
    rowname = "FunctionalInteraction"
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "ha:o:i:bdn:",
                                       ["help", "associations=", "output=", "interactions=",
                                        "biogrid", "directed","name="])
        except getopt.error as msg:
            raise Usage(msg)

        # option processing
        for option, value in opts:
            if option in ("-h", "--help"):
                raise Usage(help_message)
            elif option in ("-a", "--associations"):
                associations_file = value
            elif option in ("-o", "--output"):
                output_file = value
            elif option in ("-i", "--input"):
                interactions_file = value
            elif option in ("-b", "--biogrid"):
                biogrid = True
            elif option in ("-d", "--directed"):
                directed = True
            elif option in ("-n", "--name"):
                rowname = value
        try:
            interactions_file
            output_file
            associations_file
        except NameError:
            print("Interactions, associations and output filename must all be specified\n")
            raise Usage(help_message)
        else:
            # The runnable part of the program...
            if biogrid:
                interactions = read_biogrid(interactions_file)
            else:
                interactions = read_interactions(interactions_file, directed)
            annotate_dependencies(associations_file, output_file, interactions,
                rowname)

    except Usage as err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, "\t for help use --help"
        return 2

    return 0

if __name__ == "__main__":
    sys.exit(main())
