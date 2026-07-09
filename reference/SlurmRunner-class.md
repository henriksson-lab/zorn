# SLURM runner

Runner backend that submits Bascet jobs to a SLURM scheduler.

## Slots

- `settings`:

  Reserved settings string.

- `ncpu`:

  Number of CPU cores requested.

- `partition`:

  SLURM partition.

- `account`:

  SLURM account.

- `time`:

  SLURM time limit.

- `prepend`:

  Command prefix.

- `mem`:

  Memory limit string.

- `direct`:

  Whether to wait for completion before returning.

- `verbose`:

  Whether to print additional debug output.

- `deleteScript`:

  Whether generated scripts are deleted after submission.

- `benchmark`:

  Whether benchmark logging is enabled.

- `logTime`:

  Whether task runtime logging is enabled.
