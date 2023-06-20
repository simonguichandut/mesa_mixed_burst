This is all of the code necessary to run the [MESA](https://docs.mesastar.org/en/latest/index.html) simulations presented in Guichandut & Cumming 2023: 

[The imprint of convection on Type I X-ray bursts: Pauses in photospheric radius expansion lightcurves](https://arxiv.org/abs/2301.08769)

Designed to run with MESA [version r15140](https://zenodo.org/record/4311514).  We made some additional changes to the source code (see `mesa_source_code_changes/readme`), but these are only necessary if you wish to use our pgstar inlists, and run simulations with the "close gaps" prescriptions (see paper).

The initial model `ns_env.mod` is created in `make_env`, which is a local copy of the MESA test suite problem, with some modifications for a NS with a mass of 1.4Msun, radius 12km, and iron-56 column depth of 1e11 g/cm2.

There are 4 inlists for a single burst. They are located in `base/inlist_folder`:
- `inlist1_accrete`: accrete solar onto the fresh NS until it ignites.
- `inlist2_flash`: follow the thermonuclear runaway and convection until the outer layers reach 90% of the Eddington luminosity.
- `inlist3_relax_tau`: relax the optical depth of the outer boundary to 2/3.
- `inlist4_wind`: run the hydrodynamic wind.

`base/base_inlist` is designed to cycle through each of these inlists. The scripts `base/next_inlist` and `base/prev_inlist` can go through the main inlist file and comment/un-comment the relevant lines.

In our paper, each run starts from the same ignition model, created by the first inlist (`runs/make_ignition_model`). Then, we vary the convective prescription in the second inlist, and the third and fourth are essentially the same.

Here's how to run our main model "Schwarzschild":
1. `cd base && ./mk` (compile mesa)
2. `cp base_inlist_template base_inlist` . Replace all intances of `PATH` with your own path where the code is installed. Do the same with `inlist_folder/inlist1_accrete_template` and `inlist_folder/inlist2_flash_template` (copy to same filenames without _template, and replace PATH).
3. `cd ../runs/Schwarzschild && cp ../../base/base_inlist inlist`
4. `../base/next_inlist` Cycle once so that we start from the second inlist
5. `../base/star` Run MESA
6. Repeat 4&5 two more times to do the rest of the simulation

A template script is included (`template_run.sh`) to do steps 3 through 6 automatically.

Note that these simulations are expensive. Both the flash and wind parts take tens of thousands of models to compute. For convenience, we provide the models at 90% Eddington in the `runs/...` directories.
