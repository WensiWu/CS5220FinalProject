#
# To build with a different compiler / on a different platform, use
#     make PLATFORM=xxx
#
# where xxx is
#     icc = Intel compilers
#     gcc = GNU compilers
#     clang = Clang compiler (OS X default)
#
# Or create a Makefile.in.xxx of your own!
#

PLATFORM=icc
include Makefile.in.$(PLATFORM)

.PHONY: exe exe-vec clean realclean


# === Executables

exe: 3D_geom_nonlin_truss.x

exe: 3D_geom_nonlin_truss_vec.x

3D_geom_nonlin_truss.x: 3D_geom_nonlin_truss.o 
	$(CC) $(OMP_CFLAGS) $^ -o $@

3D_geom_nonlin_truss_vec.x: 3d_geom_nonlin_truss_vec.o
	$(CC) $(OMP_CFLAGS) $^-o $@

3D_geom_nonlin_truss.o: 3D_geom_nonlin_truss.c
	$(CC) -c $(OMP_CFLAGS) $<

3D_geom_nonlin_truss_vec.o: 3D_geom_nonlin_truss_vec.c
	$(CC) -c $(OMP_CFLAGS) $<

%.o: %.c
	$(CC) -c $(CFLAGS) $<


# === Profiling

.PHONY: maqao-cqa maqao-perf scan-build vtune-report

maqao-cqa: 3D_geom_nonlin_truss.x
	( module load maqao ; \
	  maqao cqa ./3D_geom_nonlin_truss.x fct=main uarch=HASWELL > maqao-cqa.report )

maqao-perf: 3D_geom_nonlin_truss.x
	( module load maqao ; \
	  maqao perf ./3D_geom_nonlin_truss.x fct=main uarch=HASWELL > maqao-perf.report )

scan-build:
	( module load llvm-analyzer ; \
	  scan-build -v --use-analyzer=/share/apps/llvm-3.7.0/bin/clang make )

vtune-report:
	amplxe-cl -R hotspots -report-output vtune-report.csv -format csv -csv-delimiter comma


# === Cleanup and tarball

clean:
	rm -f *.o

realclean: clean
	rm -f 3D_geom_nonlin_truss.x
