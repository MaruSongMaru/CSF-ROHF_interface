### File description
In `src/`, you can find the following files.
- `GUGA_diag.py`: Evaluates the diagonal H matrix element of a CSF.
- `gen_spinfree_rdm.py`: Generates one- and two-body RDMs of a CSF.
- `IntegralClass.py`: Reads in an FCIDUMP file.
- `runmolcas.py`: Interfaces to run CSF-ROHF using OpenMolcas RASSCF module.

### Example
After you set `molcas_path` in `src/runmolcas.py` line 72 to your Molcas installation directory,
you can run the example input in `example/`.
Go to `example/` and `$ python ../src/runmolcas.py N2 1 1 1 2 2 2`.
`N2` is the Molcas input filenmae without extension, `1 1 1 2 2 2` is the CSF
you use for the ROHF optimization in the step-vector format.
This executes Molcas with `N2.inp` and feed RDMs and the RDM energy to Molcas
for every RASSCF iteration.
