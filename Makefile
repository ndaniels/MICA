all: gofmt blosum/blosum.go install

install:
	go install -p 6 . ./compress

blosum/blosum.go:
	scripts/mkBlosum | gofmt > blosum/blosum.go

gofmt:
	gofmt -w *.go compress/*.go
	scripts/colcheck *.go compress/*.go

tags:
	find ./ \( \
			-name '*.go' \
		\) -print0 \
		| xargs -0 gotags > TAGS

loc:
	find ./ -name '*.go' -and -not -name 'blosum.go' -print | sort | xargs wc -l

push:
	git push origin master
	git push tufts master
	git push github master

pushl:
	git push origin linkedlist
	git push tufts linkedlist
	git push github linkedlist

par:
	time GOMAXPROCS=8 ./compress/compress par-coarse.fasta par-coarse.links par-compressed.cbp data/medium.fasta
