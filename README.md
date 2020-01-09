# CHTKC

CHTKC is a k-mer counting software written in C language. It is based on a lock-free hash table which uses chaining to resolve collisions.


## Installation

The latest release of CHTKC source code can be downloaded from [github][1].

To compile CHTKC, please ensure `zlib` (to support gzip-compressed inputs) and `cmake` (version 3.0.2 or higher) are installed on the target system.

The source code can be compiled using:

```shell
mkdir build
cd build
cmake ..
make
```

The `build` directory will contain two versions of CHTKC binary files:

* `chtkc`: normal version.
* `chtkco`: optimized version.

## Documentation

The documentation of CHTKC can be found [here][2].

## License

CHTKC is distributed under GNU GPL 3 license.


[1]: https://github.com/wjnjlcn/chtkc/releases
[2]: https://chtkc-doc.readthedocs.io/en/latest/