all: gofmt blosum/blosum.go install

install:
	go install -p 6 .

blosum/blosum.go:
	scripts/mkBlosum | gofmt > blosum/blosum.go

gofmt:
	gofmt -w *.go
	scripts/colcheck *.go

tags:
	find ./ \( \
			-name '*.go' \
		\) -print0 \
		| xargs -0 gotags > TAGS

loc:
	find ./ -name '*.go' -and -not -name 'blosum.go' -print | sort | xargs wc -l

