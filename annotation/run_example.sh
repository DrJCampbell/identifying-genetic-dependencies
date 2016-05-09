#!/bin/bash          
python annotate_dependencies.py -a ../results/kinome_rb_mut_associations.txt -o ../results/annotated_dependencies1.txt -i BIOGRID-GENE-111860-3.4.136.tab2.txt -n PPI -b
python annotate_dependencies.py -a ../results/annotated_dependencies1.txt -o ../results/annotated_dependencies2.txt -i RB1_regulates.txt -n GeneRegulatory -d
