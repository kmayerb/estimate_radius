# estimate_radius

example of more intensive TCR radius estimation 

![image](https://user-images.githubusercontent.com/46639063/160673001-debde1b6-8a68-4194-a143-609f38a1b583.png)


## Version 1.1

## requirements 

Must have tcrdist3, olga, and tcrsampler installed with Britanova or other appropriate data downloaded.

If you've never run tcrsampler before, do this first:

```python
import tcrsampler 
ts = tcrsampler.TCRsampler()
# For Human, Beta 
ts.download_background_file('britanova_human_beta_t_cb.tsv.sampler.tsv.zip')
# For Mouse
ts.download_background_file( 'ruggiero_mouse_sampler.zip')
```

## Run from the commandline

Estimate radius per clone using dense V-J gene matched synthetic neighborhoos

```
arguments:
  -h, --help          show this help message and exit
  -p P                path to input file
  -pout POUT          path to store outputs
  -f F                filename clones file
  -n_cpus N_CPUS      availale cpus to distribute
  -n_runs N_RUNS      number of times to estimate ideal radius
  -organism ORGANISM  human or mouse
  -chains CHAINS      comma separated list of beta beta,alpha or alpha
```

For human, `britanova_human_beta` is currently used for beta chain and `olga_human_alpha` for alpha chain.
```bash
python radii_estimation.py -p data -pout output -f dash_human.csv -n_cpus 2 -n_runs 10 -organism human -chains alpha,beta
```

For mouse, the `ruggiero` dataset is used for the background.
```bash
python radii_estimation.py -p data -pout output -f mouse_test.csv -n_cpus 2 -n_runs 3 -organism mouse -chains alpha,beta
``` 

## Recommedations 

Try this out first on a small file and set n_runs to 3
