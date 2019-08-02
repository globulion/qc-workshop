# Using Psi4 as Python Module

Psi4 can be used directly as a Python module `psi4` or 
by running Psi4 executable that reads input file, written
in Psithon: a Python language extended with some additional
features of Psi4. Some of those features are:
 * setting memory
 * setting options
 * setting molecule
 * setting custom basis set 

Inputting in Psi4 is very straightforward, so here we won't focus on
that. Example of input that runs some few jobs
is given below:

```python
# test molecule
molecule h2o {
0 1
O          -1.519303384584    -0.304272418768    -0.146864073795
H          -0.606285058641    -0.110049229060    -0.005964489169
H          -1.942559797908    -0.201234603396     0.685272417230
---
0 1
O           1.445083057756     0.288290510148     0.134486679373
H           1.993668810692    -0.443686923120    -0.085207965000
H           1.733106492991     1.008615243162    -0.397661731522

units angstrom
}

set 
{
    # ==> General Psi4 Options <== #
    step_type                     nr
    opt_coordinates               cartesian
    g_convergence                 gau_verytight
    full_hess_every               -1
    geom_maxiter                  50
    intrafrag_step_limit          0.0005
    intrafrag_step_limit_max      0.0005
    intrafrag_step_limit_min      0.0005
    dynamic_level                 0
    basis                         6-31G*
    scf_type                      direct
    guess                         core
    e_convergence                 1e-10
    d_convergence                 1e-11
    print                         1
    puream                        True
    freeze_core                   False
    onepdm                        False
    opdm_relax                    False
    cc_type                       df
    df_basis_scf                  aug-cc-pvdz-jkfit
    df_basis_cc                   aug-cc-pvdz-ri
    df_basis_sapt                 aug-cc-pvdz-ri
}

# run optimization
e, w = optimize('scf', molecule=h2o, return_wfn=True)

# compute interaction energy decomposition
energy('sapt0', molecule=h2o, ref_wfn=w)
```

The above task is a Psithon script that performs geometry optimization
of water dimer in ground state by using the HF/6-31G model. After that,
it decomposes the interaction energy by using the symmetry adapted perturbation theory (SAPT) method.

Below, we shall go through a few technical aspects that we will encounter later on during the workshop.

## Molecule object

## Basis set object

## Mints: molecular integrals

## Psi4 drivers: computing wavefunctions and more

## Play around in Python
### 2-index transformation
### 4-index transformation



-----------
[Main Page](https://github.com/globulion/qc-workshop)
