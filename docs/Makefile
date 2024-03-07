.SUFFIXES: .ipynb .md

MD := $(patsubst %.ipynb,%.md,$(wildcard md/expl/*.ipynb))

.ipynb.md:
	jupyter nbconvert --to markdown $*.ipynb
	perl -pi -e "s/\\$$\\$$/\\$$/g" $*.md
#	[ ! -f $*.md ] || mv $*.md ./md/expl/
#	[ ! -d $*_files ] || cp -rf $*_files ./md/expl/ && rm -rf $*_files

all: $(MD) 
	julia --color=yes make.jl

GSplines.jl:
	git clone -b gh-pages https://github.com/AlgebraicGeometricModeling/GSplines.jl.git

publish:
	cd GSplines.jl; git add docs; git commit -am "update gh-pages"; git push

clean:
	rm ./md/expl/*.md

