import csv, os, os.path
from itertools import islice
from collections import defaultdict
from mpmath import loggamma
from math import log, exp
from suds.client import Client
import optparse, re, shutil, gzip, urllib2, urllib


def logchoose(ni, ki):
    #n = max(ni, ki)
    #k = min(ni, ki)
    try:
        lgn1 = loggamma(ni+1)
        lgk1 = loggamma(ki+1)
        lgnk1 = loggamma(ni-ki+1)
    except ValueError:
        #print ni,ki
        raise ValueError


    return lgn1 - (lgnk1 + lgk1)

def gauss_hypergeom(X, n, m, N):
    """Returns the probability of drawing X successes of m marked items
     in n draws from a bin of N total items."""

    assert N >= m, 'Number of items %i must be larger than the number of marked items %i' % (N, m)
    assert m >= X, 'Number of marked items %i must be larger than the number of sucesses %i' % (m, X)
    assert n >= X, 'Number of draws %i must be larger than the number of sucesses %i' % (n, X)
    assert N >= n, 'Number of draws %i must be smaller than the total number of items %i' % (n, N)


    r1 = logchoose(m, X)
    try:
        r2 = logchoose(N-m, n-X)
    except ValueError:
        return 0
    r3 = logchoose(N,n)

    return exp(r1 + r2 - r3)

def hypergeo_cdf(X, n, m, N):

    assert N >= m, 'Number of items %i must be larger than the number of marked items %i' % (N, m)
    assert m >= X, 'Number of marked items %i must be larger than the number of sucesses %i' % (m, X)
    assert n >= X, 'Number of draws %i must be larger than the number of sucesses %i' % (n, X)
    assert N >= n, 'Number of draws %i must be smaller than the total number of items %i' % (n, N)
    assert N-m >= n-X, 'There are more failures %i than unmarked items %i' % (N-m, n-X)

    s = 0
    for i in range(1, X+1):
        s += max(gauss_hypergeom(i, n, m, N), 0.0)
    return min(max(s,0.0), 1)


def SearchPUBMED(search_sent, recent_date = None, BLOCK_SIZE = 100000, START = 0):

    POST_URL = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&'
    POST_URL += 'retmax=%i&' % BLOCK_SIZE
    if START > 0:
        POST_URL += 'retstart=%i&' % START

    search_url = POST_URL + '&term=' + urllib.quote_plus(search_sent)
    if recent_date:
        time_delta = datetime.today()-recent_date
        search_url += '&reldate=' + str(time_delta.days)
    
    xml_data = urllib2.urlopen(search_url).read()
    id_list = re.findall('<Id>(\d*)</Id>', xml_data)

    if len(id_list) == BLOCK_SIZE:
        return id_list + SearchPUBMED(search_sent, recent_date = recent_date,
                                      BLOCK_SIZE = BLOCK_SIZE, START = START+BLOCK_SIZE-1)
    else:
        return id_list


if __name__ == '__main__':


    parser = optparse.OptionParser("usage: %prog [options]")
    parser.add_option('-o', '--outfile', type = 'string', dest = 'outfile',
                        help = 'Output File')
    parser.add_option('-s', '--searchfile', type = 'string', dest = 'searchfile',
                        help = 'Search Term File')
    parser.add_option('-c', '--searchcache', type = 'string', dest = 'searchcache',
                        default = 'data/',
                        help = 'Cache file for storing search results')
    parser.add_option('-f', '--forceupdate', dest = 'forceupdate',
                        action = 'store_true', default = False,
                        help = 'Force updating search cache and gen2pubmed files')
    parser.add_option('-n', '--num-genes', type = 'int', dest = 'numtake',
                        default = 200)
	
    (options, args) = parser.parse_args()
    

    gene2pub_file = os.path.join(options.searchcache, 'gene2pubmed')
    if not os.path.exists(gene2pub_file) or options.forceupdate:
        url = 'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz'
        print 'downloading'
        with open(gene2pub_file + '.gz', 'w') as outhandle:
            ihandle = urllib2.urlopen(url)
            outhandle.write(ihandle.read())
        ihandle = gzip.open(gene2pub_file + '.gz')
        print 'unzipping'
        with open(gene2pub_file, 'w') as ohandle:
            ohandle.write(ihandle.read())

    print 'reading gene2pub'
    gene2pub = defaultdict(set)
    all_arts = set()
    with open(gene2pub_file) as handle:
        handle.next()
        for row in csv.DictReader(handle, fieldnames = ('taxid', 'gene', 'pubid'), delimiter = '\t'):
            if row['taxid'] == '9606':
                gene2pub[row['gene']].add(row['pubid'])
                all_arts.add(row['pubid'])
    
    genelist_dict = defaultdict(list)
    check_terms = set()
    check_pairs = set()
    print 'reading search-file'
    with open(options.searchfile) as handle:
        for row in csv.DictReader(handle, delimiter = '\t'):
            check_terms.add(row['Search-Term'])
            check_pairs.add((row['Search-Term'], row['List-Name']))
            genelist_dict[row['List-Name']].append(row['GeneID'])
    
    publist_dict = defaultdict(set)
    for term, genelist in genelist_dict.items():
        for gene in islice(genelist, options.numtake):
            publist_dict[term] |= gene2pub[gene]

    searchcache_file = os.path.join(options.searchcache, 'cachefile.txt')
    search_results = defaultdict(set)
    if os.path.exists(searchcache_file):
        with open(searchcache_file) as handle:
            for line in handle:
                term, pub = handle.strip().split('\t')
                search_results[term].add(pub)

    for term in check_terms:
        if term not in search_results:
            print 'searching for term', term
            search_results[term] = set(SearchPUBMED(term)) & all_arts
            print 'got %i articles' % len(search_results[term])

    print 'processing enrichment'
    with open(options.outfile, 'w') as handle:
        fields = ('SearchName', 'ListName', 'p-value', 'ArticleOverlap',
                    'SearchArticles', 'ListArticles')        
        writer = csv.DictWriter(handle, fields, delimiter = '\t')
        writer.writerow(dict(zip(fields, fields)))
        total_articles = len(all_arts)
        for search_term, listname in check_pairs:
            x_suc = len(publist_dict[listname] & search_results[term])
            m_marked = len(search_results[term])
            n_draws = len(publist_dict[listname])
            p = 1-hypergeo_cdf(x_suc, n_draws, m_marked, total_articles)
            print search_term, listname, p
            writer.writerow({
                    'SearchName':search_term,
                    'ListName':listname,
                    'p-value':p,
                    'ArticleOverlap':x_suc,
                    'SearchArticles':m_marked, 
                    'ListArticles':n_draws
                            })
            
            

