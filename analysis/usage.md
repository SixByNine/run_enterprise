# Short commands with defualt settings in pulsar_glitch package:

- Do model selection to evaluate models and find the best model for each glitch
`python /nvme1/yliu/yangliu/run_enterprise/analysis/model_selection.py -p bst_psrn.par -t chp_psrn.tim -u taug -g g --sigma 100 100 100 100 -r 0 1 2`

- Combine the best model for all target glitches and run stride fit
`python /nvme1/yliu/yangliu/run_enterprise/analysis/stride_plots.py -p bst_psrn.par -t chp_psrn.tim -u taug -g g -r r`

- Check results in tempo2
`tempo2 -gr plk -f fnl_psrn.par chp_psrn_all.tim`

- Write notes and comments for pulsar in documentation  
`vim psrn_comments.txt`

- Generate a summary (or a short one without the posteriors plots) latex file for pulsar
`python /nvme1/yliu/yangliu/run_enterprise/analysis/make_glitch_summary.py psrn`
`(python /nvme1/yliu/yangliu/run_enterprise/analysis/concise_glitch_summary.py psrn)`

- Complie the latex file to generate PDF file
`pdflatex -draftmode psrn_sum.tex && pdflatex psrn_sum.tex`

- Complie the short summary file to generate PDF file
`cd ..`
`pdflatex -draftmode psrn_sum.tex && pdflatex psrn_sum.tex`


## Full commands

1. Enterprise fitting:

- Emcee:

`/nvme1/yliu/yangliu/run_enterprise/run_enterprise.py --gl-all --auto-add -t 16 bst_psrn_?.par chp_psrn_?.tim (--truth-file trh_psrn.txt) --outdir psrn_? --emcee -N 8000 --nwalkers 128 --plot-chain --plot-derived -j --red-prior-log -A -8 --tspan-mult 1.1 --glitch-alt-f0 --glitch-alt-f0t 200 --alt-f0t-gltd --glitch-epoch-range 100 --measured-prior --measured-sigma 50 --glitch-td-min 1 --glitch-td-max 3 --glitch-f0d-range 3.0 --glitch-f0-range 0.8 --glitch-f1-range 0.8 --glitch-f2-range 0 --glitch-td-split 1.8 2.8 |& tee opt_psrn_?.txt`

- Dynesty:

`/nvme1/yliu/yangliu/run_enterprise/run_enterprise.py --gl-all --auto-add -t 16 bst_psrn_?.par chp_psrn_?.tim --truth-file trh_psrn.txt --outdir psrn_? --dynesty --nlive 500 --dynesty-plots --plot-derived -j --red-prior-log -A -8 --tspan-mult 1.1 --glitch-alt-f0 --glitch-alt-f0t 200 --alt-f0t-gltd --glitch-epoch-range 100 --measured-prior --measured-sigma 50 --glitch-td-min 1 --glitch-td-max 3 --glitch-f0d-range 3.0 --glitch-f0-range 0.8 --glitch-f1-range 0.8 --glitch-f2-range 0 --glitch-td-split 1.8 2.8 2>&1 | tee opt_psrn_?.txt`

- Multinest:

`/nvme1/yliu/yangliu/run_enterprise/run_enterprise.py --gl-all --auto-add -t 16 bst_psrn_?.par chp_psrn_?.tim --truth-file trh_psrn.txt --outdir psrn_? --multinest --plot-chain --plot-derived -j --red-prior-log -A -8 --tspan-mult 1.1 --glitch-alt-f0 --glitch-alt-f0t 200 --alt-f0t-gltd --glitch-epoch-range 100 --measured-prior --measured-sigma 50 --glitch-td-min 1 --glitch-td-max 3 --glitch-f0d-range 3.0 --glitch-f0-range 0.8 --glitch-f1-range 0.8 --glitch-f2-range 0 --glitch-td-split 1.8 2.8 |& tee opt_psrn_?.txt`

2. Tempo2:

`tempo2 -gr plk -f bst_psrn_?.par.post chp_psrn_?.tim`

3. Stride fitting:

`python stride_plots.py -p fnl_psrn_?.par -t chp_psrn_?.tim -u taug`

4. Make comments:

`vim psrn_comments.txt`

5. Create latex:

`python /nvme1/yliu/yangliu/yang/scripts/make_glitch_summary.py psrn`

6. Convert latex to pdf:

`pdflatex -draftmode psrn_sum.tex && pdflatex psrn_sum.tex`


## Others

- Create new parameter files

`tempo2 - f in.par in.tim -newpar`

- Get the value of parameters from parameter files

`grep GLF2_1 new.par`


# Defualt definitions

## Large glitch:

dF/F>10^-6

## Default GLTD for model selection:

GLTD=50, GLTD2=200, GLTD3=800

## Defualt split for mulitple recoveries:

- 2 recoveries: Split at 10^2.2 (~158)

- 3 recoveries: Split at 10^1.8 and 10^2.8 (~63/630)

### Exponentials of 10:

1.0: 10.000;

1.1: 12.589;

1.2: 15.848;

1.3: 19.952;

1.4: 25.118;

1.5: 31.622;

1.6: 39.810;

1.7: 50.118;

1.8: 63.095;

1.9: 79.432.


# Model Abbreviation

- N: fit GLF0 and GLF1 with tempo2 (for tiny glitches)

- F: 0 recovery

- R: 1 recovery

- D: 2 recoveries

- T: 3 recoveries

- G: extra tiny glitches

- n: do not fit for GLF2

- y: fit for GLF2
