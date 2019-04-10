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

The `build` directory will contain two versions of CHTKC binary files -- the normal version `chtkc` and the optimized version `chtkco`.


## Note

Under small RAM usage, the optimized version has better performance, because it can allocate more nodes than the normal version.

However, because the optimized version can only use up to 4G (4294967295) nodes, when the available memory is big enough, the normal version may instead allocate more nodes than the optimized the version.

For example, when counting 28-mers, if the RAM usage is over 150 GB, please consider using the normal version.


## License

CHTKC is distributed under GNU GPL 3 license.


[1]: https://github.com/wjnjlcn/chtkc/releases
