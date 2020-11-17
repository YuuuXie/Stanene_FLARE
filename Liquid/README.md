Folder `LMP`: molecular dynamics using LAMMPS

- folder `Atoms10k_Temp500K`: production simulation with 10000 atoms at around 500K, simulating the melting process

- `in.lammps`, `data.lammps`: melting process simulation with 64 atoms, to generate an initial configuration for AIMD

Folder `AIMD`: ab-initio molecular dynamics of melted tin

- `init_struc.xyz`: the initial structure obtained from LAMMPS MD

- `md.in`: the input file of AIMD. The cell is set based on experimental density of liquid tin, such that there's no vacuum in the cell.

- The AIMD output is too large, and is uploaded to the Materials Cloud Archive
