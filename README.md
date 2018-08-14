# Linear algebra library for Futhark

A small collection of linear algebra operations written in Futhark.
There's not a lot here yet; please open an issue if you need something
that's missing.

## Installation

```
$ futhark-pkg add github.com/diku-dk/sorts
$ futhark-pkg sync
```

## Usage

```
$ futharki
> import "lib/github.com/diku-dk/sorts/radix_sort"
> radix_sort_int i32.num_bits i32.get_bit [5,7,-1,2,-2]
[-2i32, -1i32, 2i32, 5i32, 7i32]
```

## See also

* https://github.com/athas/distance
* https://github.com/athas/vector
