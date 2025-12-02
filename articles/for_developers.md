# Installation

## Compiling Bascet

Bascet is located [in a separate Github
repository](https://github.com/henriksson-lab/bascet)

You first need to [download rust](https://rustup.rs)

To download Bascet and get the Rust dependencies:

``` console
git clone https://github.com/henriksson-lab/bascet.git   
#or: git clone git@github.com:henriksson-lab/bascet.git   

cd bascet
make install_rust  #this will install rust nightly build -- only need to do this once
```

To compile Bascet, you just run “make” whenever you make a change:

``` console
cd bascet
make
```

## Editing the Bascet code

To modify the code, we recommend
[VSCode](https://code.visualstudio.com/Download) and in particular,
install the [rust code
extensions](https://code.visualstudio.com/docs/languages/rust).

This setup enables you to see compile errors and get hints while coding!

## Accessing a locally compiled Bascet in Zorn

You can do this to link to your own Rust binary:

``` r
bascetInstance.default <- BascetInstance(
  bin="/path/to/bascet/target/debug/bascet",
  tempdir="./"
)
```

If you want access to all the software in Docker/Singularity, you can
also pass an image; while Bascet is present, providing an absolute path
like above will override it.
