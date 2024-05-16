# CHAMFER AMM15 Configuration Based on CO9_AMM15_P2.0

MY_SRC and EXPREF for 1993 initialisation of
 [CO9_AMM15_P2.0](https://github.com/JMMP-Group/CO_AMM15/tree/CO9). 
The configuration diverges as follows:

- GLOSEA6 lateral forcing for both the atlantic and baltic boundary conditions
- MY_SRC includes Momentum Trend capability with Anthony Wise's atmospheric
  pressure additions

## ARCHER2 Specifics
The simulation has been configured on ARCHER2. BDY, SBC and RST forcing files
are available on at:

```
/work/n01/n01/shared/CO_AMM15/CHAMFER
```

The model is compiled using the ARCHER2 compilation scripts that is packaged
with NEMO (`arch-X86_ARCHER2-Cray.fcm`).

The XIOS executable version is stored on ARCHER2:
```
/work/n01/shared/nemo/XIOS2_Cray/
```

The NEMO version is a non-standard branch that downloaded from the subversion
repository using
```
 svn -r 15194 co https://forge.ipsl.jussieu.fr/nemo/svn/NEMO/branches/UKMO/NEMO_4.0.4_momentum_trends nemo_4.0.4_trd
```

This verion allows the use of momentum trend diagnostics, which have been broken
for all NEMO releases since v4.0
