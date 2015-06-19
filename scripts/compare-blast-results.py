from Bio import SearchIO
from sys import argv

def parseBlastOutFile(filename):
	if filename[-3:] == "xml":
		qResultGen = SearchIO.parse(filename, 'blast-xml')
	elif filename[-3:] == "txt":
		qResultGen = SearchIO.parse(filename, 'blast-tab')
	else:
		print("Unrecognized filetype.")
		assert False

	parsed = {qRes.id : qRes for qRes in qResultGen}
	print("Parsed "+filename)

	return parsed

def analyzeBlastOut(gold, comp):
	# make sure we have only use queries which both results have
	intersection = gold.viewkeys() & comp.viewkeys()

	totalMiss = 0
	totalGoldHits = 0
	perfectQueries = 0
	badQueries = 0
	veryBadQueries = 0
	for qId in intersection:
		goldHits = len(gold[qId])
		totalGoldHits += goldHits

		miss = len(gold[qId]) - len(comp[qId])
		assert miss >= 0
		totalMiss += miss
		if miss == 0:
			perfectQueries += 1

		missFrac = float(miss) / goldHits 
		if missFrac > 0.1:
			badQueries += 1
		if missFrac > 0.2:
			veryBadQueries += 1

	print "Queries: {}".format(len(intersection))
	print "Misses: {}".format(totalMiss)
	print "Total hits in gold-standard: {}".format(totalGoldHits)
	print "Overall accuracy: {}".format( 1 - float(totalMiss)/totalGoldHits)
	print "Queries with perfect accuracy: {}".format(  perfectQueries)
	print "Queries with greater than 90% of hits: {}".format(  len(intersection) - badQueries)
	print "Queries with greater than 80% of hits: {}".format(  len(intersection) - veryBadQueries)


def main(args):
	if len(args) != 3:
		print("Usage: [gold-standard] [comparison]")
		return
	else:
		analyzeBlastOut( parseBlastOutFile(args[1]), parseBlastOutFile(args[2]))


if __name__ == '__main__':
	main(argv)