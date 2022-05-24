FUTHARK?=futhark

.PHONY: test
test: lib/github.com/diku-dk/cpprandom
	$(MAKE) -C lib/github.com/diku-dk/linalg test

.PHONY: doc
doc: lib/github.com/diku-dk/cpprandom
	$(MAKE) -C lib/github.com/diku-dk/linalg doc

.PHONY: clean
clean:
	rm -rf *~
	$(MAKE) -C lib/github.com/diku-dk/linalg clean
	$(MAKE) -C benchmarks clean

.PHONY: realclean
realclean: clean
	rm -rf lib/github.com/diku-dk/cpprandom

lib/github.com/diku-dk/cpprandom:
	$(FUTHARK) pkg sync
