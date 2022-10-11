# gtooljl

This package reads and writes Gtool3 format files.
It supports UR4 UR8 MR4 MR8 formats.

# Install
In Julia,
```
julia> ]
pkg> add https://github.com/mchwavy/gtooljl.git
```

# Usage

To read header information,
```
using gtoolgl

dir = "."
varname = "occo2f"

filename = dir * "/" * varname
f = opengtool( filename, "r" )
chead = readchead( f )
```


To read data,

```
using gtoolgl

dir = "."
varname = "occo2f"

nx, ny, nz, nt, lon, lat, dep, tarray, array = readgtool( dir, varname )
```
readgtool uses opengtool.

To write,
```
using gtoolgl

f = opengtool( filename, "w" )

nxyz = nx * ny * nz
chead[64]=@sprintf("%16s", nxyz)

for it=1:12
    arrayxy = array2[:, :, it]
    chead[25] = @sprintf("%16i", tarray[it])
    writegtool( f, chead, "MR8", arrayxy )
end

closegtool( f )
```


[![Build Status](https://travis-ci.com/mchwavy/gtooljl.jl.svg?branch=main)](https://travis-ci.com/mchwavy/gtooljl.jl)
[![Coverage](https://codecov.io/gh/mchwavy/gtooljl.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mchwavy/gtooljl.jl)
[![Coverage](https://coveralls.io/repos/github/mchwavy/gtooljl.jl/badge.svg?branch=main)](https://coveralls.io/github/mchwavy/gtooljl.jl?branch=main)
