from __future__ import with_statement
from BeautifulSoup import BeautifulStoneSoup
from suds.client import Client
import optparse, os.path, os, re, shutil
from itertools import islice
import logging
from multiprocessing import Pool
from itertools import imap


#create a silly little mock
class fake_memcache():
	class Client():
		def __call__(self, *args, **kwargs):
			return None
		def __init__(self, *args, **kwargs):
			pass
		def set(self, *args, **kwargs):
			return None
		def get(self, *args, **kwargs):
			return None
		


try:
	import memcache
except ImportError:
	memcache = fake_memcache

def nth(iterable, n, default=None):
	"Returns the nth item or a default value"
	return next(islice(iterable, n, None), default)
def take(n, iterable):
	"Return first n items of the iterable as a list"
	return list(islice(iterable, n))

def TermGen(args, resume_line):
	counter = 0
	for f in args:
		with open(f) as handle:
			for line in handle:
				counter += 1
				if counter < resume_line:
					continue
				parts = line.strip().split('\t')
				
				yield parts[0], parts[1], counter

def GetFromCache(line_num, fname):
	with open(fname) as handle:
		it = nth(handle, line_num)
	parts = it.split('\t')
	return parts


def AskWhatzit(client, uni_prot, query):
	
	out = client.service.search(pipelineName = 'whatizitSwissprotDisease',
								query = query+' AND ' + uni_prot, 
								limit = 5000)
	return out
	
def ParseXML(xml_data):
	
	lines = []
	for i, line in enumerate(xml_data.split('\n')):
		if 'MedlineCitationSet'.lower() in line.lower():
			lines.append(i-1)
	if len(lines) < 2:
		return []
	parts = islice(xml_data.split('\n'), lines[0], lines[1])
	ptree = BeautifulStoneSoup(''.join(parts))

	pmids = []
	for cit in ptree.findAll('medlinecitation'):
		pmids.append(cit.pmid.string)
	
	
	return pmids
	

def AddToCache(fname, uni_prot, query, pmids, d):
	with open(fname, 'a') as handle:
		handle.write('\t'.join([uni_prot, query, 
										str(len(pmids)), ','.join(pmids)])+'\n')
	d[(uni_prot, query)] = max(d.keys())+1
	

	
def RetrieveChunk(chunk_tuple):
	
	logging.warning('Starting SOAP')
	url = 'http://www.ebi.ac.uk/webservices/whatizit/ws?wsdl'
	client = Client(url, faults = False, retxml = True)
	
	
	mc = memcache.Client(['127.0.0.1:11211'])
	
	output = []
	for query, uni_prot, cur_line in chunk_tuple:
		logging.warning('Retrieving: %s at line %d' % (uni_prot, cur_line))
		val = mc.get(str(hash(uni_prot+'_'+query)))
		if val:
			logging.warning('Used Cache!!: %s at line %d' % (uni_prot, cur_line))
			output.append(val.split('_'))
			continue
		xmlout = AskWhatzit(client, uni_prot, query)
		pmids = ParseXML(xmlout)
		output.append([uni_prot, query, str(len(pmids)), ','.join(pmids)])
		val = mc.set(str(hash(uni_prot+'_'+query)), 
				'\t'.join([uni_prot, query, 
				str(len(pmids)), ','.join(pmids)]))
		
	return output
	
	
def ChunkGen(size, iterable):
	
	block = take(size, iterable)
	while block:
		yield block
		block = take(size, iterable)
	
	
	
if __name__ == '__main__':
	
	
	parser = optparse.OptionParser("usage: %prog [options]")
	parser.add_option('-o', '--outfile', type = 'string', dest = 'outfile',
						help = 'Output File')
	parser.add_option('-c', '--searchcache', type = 'string', dest = 'searchcache',
						default = 'C:\\publookup\\cachefolder\\',
						help = 'Cache file for storing search results')
	parser.add_option('-w', '--workers', type = 'int', dest = 'numworkers',
						default = 3,
						help = 'Number of workers to use')
	parser.add_option('-f', '--forcenew', action = 'store_true', dest = 'forcenew',
						default = False,
						help = 'Forces the restarting of the retrieval')
	parser.add_option('-u', '--updatecache', action = 'store_true', dest = 'updatecache',
						default = False,
						help = 'Updates the cache file')
	
	if options.updatecache:
		memcache = fake_memcache
	
	(options, args) = parser.parse_args()
	
	line_search_dict = {}

	if options.searchcache and USE_MEM and not options.updatecache:
		mc = memcache.Client(['127.0.0.1:11211'])
		if os.path.isdir(options.searchcache):
			files = map(lambda x: os.path.join(options.searchcache, x), 
						os.listdir(options.searchcache))
		else:
			files = (options.searchcache, )
		for f in files:
			with open(f) as handle:
				for line in handle:
					parts = line.strip().split('\t')
					ind = (parts[0], parts[1])
					mc.set(str(hash('_'.join(ind))), '_'.join(parts)) 
	
	restart_line = 0
	if os.path.exists(options.outfile) and not options.forcenew and not options.updatecache:
		with open(options.outfile) as handle:
			for i,line in enumerate(handle):
				pass
		try:
			restart_line = i
		except:
			pass
	elif os.path.exists(options.outfile) and options.forcenew:
		with open(options.outfile, 'w') as handle:
			pass

	if options.updatecache:
		options.outfile = os.path.join(options.searchcache, options.outfile)
		args = map(lambda x: os.path.join(options.searchcache, x), 
						os.listdir(options.searchcache))
			
	for un, q, counter in TermGen(args, 0):
		pass
	
	pool = Pool(options.numworkers)
	with open(options.outfile, 'a') as handle:
		
		for out in pool.imap(RetrieveChunk, 
				ChunkGen(10, TermGen(args, restart_line)), options.numworkers):
			for sout in out:
				logging.warning('Writting: %s' % sout[0])
				handle.write('\t'.join(map(str, sout))+'\n')
			handle.flush()
			os.fsync(handle.fileno())
		
	
	
	
	
	
	
	
	
	
	
	
	
	