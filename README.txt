Before beginning, set up the following directories:

$ ln -s /w/hallc-scshelf2102/nps/$USER/ROOTfiles ROOTfiles
$ mkdir hist
$ mkdir cuts
$ mkdir plots

Then,
1.  hcana [0] .x make_hist_hms_optics.C(1544,kFALSE,kFALSE,kFALSE,-1)
2.  hcana [0] .x set_ytar_delta_cuts.C(1544,-1)
3.  hcana [0] .x make_hist_hms_optics.C(1544,kTRUE,kFALSE,kFALSE,-1)

4a.  hcana [0] .x plot_yfp_cuts.C(1544,-1)
4b.  hcana [0] .x set_ypfp_yfp_cuts.C(1544,-1,5)
4c.  hcana [0] .x plot_yfp_cuts.C(1544,-1)

5a.  hcana [0] .x plot_xfp_cuts.C(1544,-1)
5b.  hcana [0] .x set_xpfp_xfp_cuts.C(1544,-1,5)
5c.  hcana [0] .x plot_xfp_cuts.C(1544,-1)

6. hcana [0] .x make_fit_ntuple.C(1544,-1)

7. hcana [0] .x plot_yptar_diff.C(1544,-1)
8. hcana [0] .x plot_ytar_diff.C(1544,-1)
9. hcana [0] .x plot_xptar_diff.C(1544,-1)

10. hcana [0] .x fit_opt_matrix.C()

--------- With New Elements --------------
copy matrix file and remove everything up to first order line, and use this with the following:

a.  hcana [0] .x make_hist_hms_optics_v2.C(1544,kTRUE,kFALSE,kFALSE,-1)
   a'. Well, you may need to remake the yptar cuts, and make sure you're not overwriting cuts files so change the extension names or something.
b.  hcana [0] .x plot_yfp_cuts.C(1544,-1)
c.  hcana [0] .x plot_xfp_cuts.C(1544,-1)

Things to check - ztar vs ztarg vs ztarcalc within make_hist_hms_optics_v2.C


