all: gofmt blosum/blosum.go install

install:
	go install -p 6 . \
		./cmd/cablastp-compress ./cmd/cablastp-decompress \
		./cmd/cablastp-search ./cmd/cablastp-psisearch \
		./cmd/cablastp-deltasearch

blosum/blosum.go:
	scripts/mkBlosum | gofmt > blosum/blosum.go

gofmt:
	gofmt -w *.go cmd/*/*.go
	scripts/colcheck *.go cmd/*/*.go

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

