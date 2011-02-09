from itertools import product
import os, os.path, optparse, csv

if __name__ == '__main__':

    parser = optparse.OptionParser("usage: %prog [options]")

    parser.add_option('-o', '--outfile', type = 'string', dest = 'outfile',
						help = 'Output File')
    parser.add_option('-d', '--in-direc', type = 'string', dest = 'indirec',
                        default = None,
                        help = 'Input Directory')
    parser.add_option('-s', '--searchfile', type = 'string', dest = 'searchfile',
						default = None,
						help = 'Search file')
	
    (options, args) = parser.parse_args()

    with open(options.searchfile) as handle:
        search_terms = set()
        for line in handle:
            search_terms.add(line.strip())

    gene_list_files = os.listdir(options.indirec)
    print gene_list_files
    gene_list_dict = {}
    for f in gene_list_files:
        listname = f.split('.')[0]
        with open(os.path.join(options.indirec, f)) as handle:
            gene_list_dict[listname] = [x.strip() for x in handle if len(x.strip())]

    with open(options.outfile, 'w') as handle:
        fields = ('Search-Term', 'List-Name', 'GeneID')
        writer = csv.DictWriter(handle, fields, delimiter = '\t')
        writer.writerow(dict(zip(fields, fields)))
        for listname, search_term in product(gene_list_dict.keys(), search_terms):
            for gene in gene_list_dict[listname]:            
                writer.writerow({'Search-Term': search_term,
                                'List-Name':listname,
                                'GeneID':gene})
                
    
    
