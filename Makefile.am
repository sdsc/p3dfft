EXTRA_DIST = NEWS README AUTHORS ChangeLog P3DFFT_USER_GUIDE_2.4.pdf \
		build/bcomm1.F90 build/bcomm1_trans.F90 build/bcomm2.F90 \
		build/btran.F90 build/fcomm1.F90 build/fcomm2.F90 \
		build/fcomm2_trans.F90 build/ftran.F90 build/ghosts.F90 \
		build/init_plan.F90 build/reorder.F90 build/setup.F90 

include_HEADERS = include/p3dfft.h include/p3dfft.mod

SUBDIRS = \
	build \
	sample

all-local: 
	echo $(SUCCESS)
clean-local:
	-test -z "include/p3dfft.mod" || rm -f include/p3dfft.mod
