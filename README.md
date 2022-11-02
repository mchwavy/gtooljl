# gtooljl

This package is a Julia Language module that reads and writes Gtool3 format files.
It supports UR4 UR8 MR4 MR8 and URY formats.

# Install
In Julia,
```
julia> ]
pkg> add https://github.com/mchwavy/gtooljl.git
```

# Usage

To read header information,
```
using gtooljl

dir = "."
varname = "occo2f"

filename = dir * "/" * varname
f = opengtool( filename, "r" )
chead = readchead( f )
```

To read data, you need specify the directory containing the axis data files in advance in the GTAXDIR environment variable.

```
using gtooljl

dir = "."
varname = "occo2f"
filename = dir * "/" * varname

nx, ny, nz, nt, lon, lat, dep, tarray, array = readgtool( filename )
```
readgtool uses opengtool and readchead internally, and reads axis files.
To get information on how many data are available in the time direction, readgtool uses the ngtstat command, which was developed to handle gtool3 files.
Or you can also tell readgtool how many data to read in the time direction.
```
ntin = 1
nx, ny, nz, nt, lon, lat, dep, tarray, array = readgtool( filename, ntin )
```

To obtain the size of an array in gtool file,
```
nx, ny, nz, nt = gtoolarraysize( filename, 1 )
```

To write,
```
using gtooljl

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
The type of chead is Vector{FString{16}} and chead can be initialized as:
```
chead=Vector{FString{16}}(undef, 64)
```


To convert time data in gtool format file to DateTime, call gtooltime2datetime:
```
date = gtooltime2datetime( tarray[1] )
date = gtooltime2datetime.( tarray )
```
To convert time data in DateTime to gtool format, call datetime2gtooltime:
```
d = DateTime(0001,1,1,0,0,0)
t = datetime2gtooltime( d )
```


[![Build Status](https://travis-ci.com/mchwavy/gtooljl.jl.svg?branch=main)](https://travis-ci.com/mchwavy/gtooljl.jl)
[![Coverage](https://codecov.io/gh/mchwavy/gtooljl.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mchwavy/gtooljl.jl)
[![Coverage](https://coveralls.io/repos/github/mchwavy/gtooljl.jl/badge.svg?branch=main)](https://coveralls.io/github/mchwavy/gtooljl.jl?branch=main)
