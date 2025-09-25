FUTHARK?=futhark

TESTS=lu_tests.fut linalg_tests.fut qr_tests.fut nmf_tests.fut lup_tests.fut perm_tests.fut

.PHONY: test
test: $(TESTS)
	$(FUTHARK) test --pass-compiler-option=--Werror $^

.PHONY: doc
doc:
	$(FUTHARK) doc -o doc .

.PHONY: clean
clean:
	rm -rf *~ *.c lu_tests linalg_tests nmf_tests qr_tests linalg lup_tests perm_tests nmf lu qr *.actual *.expected doc
