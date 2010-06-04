EXTRA_DIST = NEWS README AUTHORS ChangeLog include

include_HEADERS = include/p3dfft.h include/p3dfft.mod

SUBDIRS = \
	build \
	sample

all-local: 
	echo $(SUCCESS)
clean-local:
	-test -z "include/p3dfft.mod" || rm -f include/p3dfft.mod
