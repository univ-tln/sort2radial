BINDIR          = bin
SRCDIR          = src

default all: sort2radial

sort2radial: 
	-cd $(SRCDIR) && make
	mkdir -p $(BINDIR)
	mv  $(SRCDIR)/sort2radial $(BINDIR)/

test: 
	cd example && make

clean:
	-rm -rf $(BINDIR)
	-cd src && make clean
	-cd example && make clean
