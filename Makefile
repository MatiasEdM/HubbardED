HUBB.x : main.F90 inout.F90 utilities.F90 lattice.F90 realspaceDET.F90 hopMAP.F90 doubleOCC.F90 hamiltonian.F90 sign_problem.F90 stoq.F90 truncation.F90 wrap_lapack.F90 
	gfortran -o HUBB.x main.F90 inout.F90 utilities.F90 lattice.F90 realspaceDET.F90 hopMAP.F90 doubleOCC.F90 hamiltonian.F90 indxDET.F90 energies.F90 spin.F90 sign_problem.F90 stoq.F90 truncation.F90 fnaBULK.F90 fnaRELAX.F90 fnaUTIL.F90 kspaceHUBB.F90 wrap_lapack.F90 -llapack -lblas

clean :
	rm -rf *out *png

plotwfc :
	python3 plot_wfc.py
