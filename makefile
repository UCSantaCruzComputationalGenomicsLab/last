CXXFLAGS = -O3
all:
	@cd src && $(MAKE) CXXFLAGS="$(CXXFLAGS)"

prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
install: all
	mkdir -p $(bindir)
	cp src/last?? src/last-split src/last-merge-batches scripts/* $(bindir)

clean:
	@cd src && $(MAKE) clean

html:
	@cd doc && $(MAKE)

distdir = last-`hg id -n`

RSYNCFLAGS = -rC --exclude 'last??' --exclude last-split --exclude last-merge-batches

dist: log html
	@cd src && $(MAKE) version.hh
	rsync $(RSYNCFLAGS) doc examples makefile scripts src *.txt $(distdir)
	zip -qrm $(distdir) $(distdir)

log:
	hg log --style changelog > ChangeLog.txt
