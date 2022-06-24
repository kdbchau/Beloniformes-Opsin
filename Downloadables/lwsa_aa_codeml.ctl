      seqfile = lwsa_msa_1030.phy	** sequence data filename
     treefile = bel_38_tree.phy 	** tree structure file name
      outfile = lwsa_msa_1030_AA_OUT	** main result file name
      
   aaRatefile = ../dat/lg.dat * empirical aa models (.dat file)
        ncatG = 4 * number of rate classes
    fix_alpha = 0 * fix gamma shape parameter
        alpha = 0.5 * initial or fixed alpha
        model = 3 * 2: Empirical; 3: Empirical + F

      seqtype = 3 * 2:AAs; 3:codons->AAs
        icode = 0 * 0:universal code; 1:mammalian mt; etc...
    cleandata = 0 * 1:remove sites with ambiguity data

 RateAncestor = 1 * ancestral states
        getSE = 0 * 1: want S.E.s of estimates
      	noisy = 9 * how much data displayed during run 
	verbose = 2 * how detailed the rst file it
      runmode = 0 * 0: user tree
        clock = 0 * 0:no clock
        ndata = 1 * number of data sets
   Small_Diff = .5e-6 * small value used in difference approximation
  fix_blength = 0 * 0: ignore, -1: random; 1: initial; 2: fixed
       method = 1 * 0: simultaneous; 1: one branch at a time
