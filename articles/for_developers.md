# Installation

## Compiling Bascet

Bascet is located [in a separate Github
repository](https://github.com/henriksson-lab/bascet).

First install Rust from
[rust-lang.org](https://rust-lang.org/tools/install/):

``` bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

Then download bascet

``` console
git clone https://github.com/henriksson-lab/bascet.git
```

and install what is needed for a nightly build:

``` console
cd bascet
make install_rust  #this will install rust nightly build -- only need to do this once
```

To compile Bascet, you just run “make” whenever you make a change:

You have to make a release build before using Bascet from Zorn. Note if
you are editing the source code: Zorn will not ensure that it is
recompiled when you run it! Rerun whenever needed

``` console
make
```

## Editing the Bascet code

To modify the code, we recommend
[VSCode](https://code.visualstudio.com/Download) and in particular,
install the [rust code
extensions](https://code.visualstudio.com/docs/languages/rust).

This setup enables you to see compile errors and get hints while coding!
