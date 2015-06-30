all: gofmt blosum/blosum.go install

install:
	go install -p 6 . \
		./cmd/mica-compress ./cmd/mica-decompress \
		./cmd/mica-search ./cmd/mica-psisearch \
		./cmd/mica-deltasearch ./cmd/mica-xsearch

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

