# Pick random KMERs from KMC3 database. The choice is among KMERs within a frequency range

Pick random KMERs from KMC3 database. The choice is among KMERs within a
frequency range

## Usage

``` r
ChooseInformativeKMERs(kmerHist, minFreq = 0.005, maxFreq = 1)
```

## Arguments

- kmerHist:

  Data.frame with KMER frequency

- minFreq:

  Pick KMERs having minimum frequency

- maxFreq:

  Pick KMERs having maxiumum frequency

## Value

List of KMERs
