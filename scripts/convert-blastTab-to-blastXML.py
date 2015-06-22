from Bio import SearchIO
from sys import argv




def main(args):
	if len(args) == 2:
		filenameRoot = args[1].split(".")[0]
		filenameXML = filenameRoot + ".xml"
		SearchIO.convert(args[1], 'blast-tab', filenameXML, 'blast-xml')

	elif len(args) == 3:
		filenameRoot = args[1].split(".")[0]
		filenameXML = args[2]
		SearchIO.convert(args[1], 'blast-tab', filenameXML, 'blast-xml')

	else:
		print("Usage: path/to/blast/tabular/file [optional path/for/new/blast/xml/file]")


if __name__ == '__main__':
	main(argv)