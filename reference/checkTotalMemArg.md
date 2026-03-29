# Given memory amount in mb, format it in a format suitable for Rust parsing

Helper for functions. Set total memory argument based on container and
runner, if not user specified

## Usage

``` r
checkTotalMemArg(totalMem, runner, bascetInstance)
```

## Arguments

- totalMem:

  Total memory to allocate, as a string (e.g. "8g"). If NULL, derived
  from runner

- runner:

  The job runner, used to extract memory settings if totalMem is NULL

- bascetInstance:

  A Bascet instance, used to subtract container memory overhead
