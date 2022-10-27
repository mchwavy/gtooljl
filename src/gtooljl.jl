module gtooljl

export opengtool, closegtool, writegtool, readgtool, readaxis, readchead, ipack32len, unpack_bits_from32, unpack_bit, pack_bits_into32, assert, gtooltime2datetime

using FortranFiles
using OffsetArrays
using Printf
using Dates

function assert( in1, in2, process )
    if in1 != in2
        println(in1, " should be ", in2, " in ", process)
        exit(1)
    end
end

function ipack32len( len_orig, nbit )

    bwidth = 32

    len = nbit * floor(Int, len_orig / bwidth ) + floor( Int, (nbit * mod(len_orig, bwidth) + bwidth - 1) / bwidth)
    return len
end

function unpack_bit( ilen, ipack1 )

    ipack0 = OffsetArray(ipack1, 0:length(ipack1)-1)

    nfull = ilen >>> 5

    idata1 = zeros( Int64, ilen )
    idata0 = OffsetArray(idata1, 0:ilen-1)

    for m = 0 : nfull-1
        ival = ipack0[ m ]
        n = m << 5

        idata0[ n +  0 ] = 1 &  (ival >> 31)
        idata0[ n +  1 ] = 1 &  (ival >> 30)
        idata0[ n +  2 ] = 1 &  (ival >> 29)
        idata0[ n +  3 ] = 1 &  (ival >> 28)
        idata0[ n +  4 ] = 1 &  (ival >> 27)
        idata0[ n +  5 ] = 1 &  (ival >> 26)
        idata0[ n +  6 ] = 1 &  (ival >> 25)
        idata0[ n +  7 ] = 1 &  (ival >> 24)
        idata0[ n +  8 ] = 1 &  (ival >> 23)
        idata0[ n +  9 ] = 1 &  (ival >> 22)
        idata0[ n +  10 ] = 1 & ( ival >> 21)
        idata0[ n +  11 ] = 1 & ( ival >> 20)
        idata0[ n +  12 ] = 1 & ( ival >> 19)
        idata0[ n +  13 ] = 1 & ( ival >> 18)
        idata0[ n +  14 ] = 1 & ( ival >> 17)
        idata0[ n +  15 ] = 1 & ( ival >> 16)
        idata0[ n +  16 ] = 1 & ( ival >> 15)
        idata0[ n +  17 ] = 1 & ( ival >> 14)
        idata0[ n +  18 ] = 1 & ( ival >> 13)
        idata0[ n +  19 ] = 1 & ( ival >> 12)
        idata0[ n +  20 ] = 1 & ( ival >> 11)
        idata0[ n +  21 ] = 1 & ( ival >> 10)
        idata0[ n +  22 ] = 1 & ( ival >> 9)
        idata0[ n +  23 ] = 1 & ( ival >> 8)
        idata0[ n +  24 ] = 1 & ( ival >> 7)
        idata0[ n +  25 ] = 1 & ( ival >> 6)
        idata0[ n +  26 ] = 1 & ( ival >> 5)
        idata0[ n +  27 ] = 1 & ( ival >> 4)
        idata0[ n +  28 ] = 1 & ( ival >> 3)
        idata0[ n +  29 ] = 1 & ( ival >> 2)
        idata0[ n +  30 ] = 1 & ( ival >> 1)
        idata0[ n +  31 ] = 1 & ( ival )
    end

    id0size = ((nfull-1) << 5) + 31

    if (ilen & 31) != 0 
        ival = ipack0[ nfull ]
        n = nfull << 5
        for i = 0:ilen-n-1
            idata0[ n + i ] = 1 & (ival >> (31 - i) )
        end
        id0size = ilen - 1
    end

    idata1 = OffsetArrays.no_offset_view(idata0)
    id1size = id0size + 1
    
    return idata1, id1size
end

function pack_bit( idata1, len_orig )

    width = 32
    mask = 1

    ilen2 = floor( Int, len_orig / width )

    packed1 = zeros(Int32, ilen2+1)
    packed0 = OffsetArray(packed1, 0:ilen2)

    idata0 = OffsetArray(idata1, 0:length(idata1)-1)

    ilen = ilen2
    if len_orig > width * ilen2
        ilen = ilen + 1
    end
   
    for i = 0:ilen2-1
        m  = width * i
        tmp = (  ( ( idata0[ m + 0] & mask ) <<  31 )
                      + ( ( idata0[ m + 1] & mask ) <<  30 ) 
                      + ( ( idata0[ m + 2] & mask ) <<  29 )
                      + ( ( idata0[ m + 3] & mask ) <<  28 ) 
                      + ( ( idata0[ m + 4] & mask ) <<  27 )
                      + ( ( idata0[ m + 5] & mask ) <<  26 )
                      + ( ( idata0[ m + 6] & mask ) <<  25 )
                      + ( ( idata0[ m + 7] & mask ) <<  24 )
                      + ( ( idata0[ m + 8] & mask ) <<  23 )
                      + ( ( idata0[ m + 9] & mask ) <<  22 )
                      + ( ( idata0[ m + 10] & mask ) << 21 )
                      + ( ( idata0[ m + 11] & mask ) <<  20 )
                      + ( ( idata0[ m + 12] & mask ) <<  19 )
                      + ( ( idata0[ m + 13] & mask ) <<  18 )
                      + ( ( idata0[ m + 14] & mask ) <<  17 )
                      + ( ( idata0[ m + 15] & mask ) <<  16 )
                      + ( ( idata0[ m + 16] & mask ) <<  15 )
                      + ( ( idata0[ m + 17] & mask ) <<  14 )
                      + ( ( idata0[ m + 18] & mask ) <<  13 )
                      + ( ( idata0[ m + 19] & mask ) <<  12 )
                      + ( ( idata0[ m + 20] & mask ) <<  11 )
                      + ( ( idata0[ m + 21] & mask ) <<  10 )
                      + ( ( idata0[ m + 22] & mask ) <<  9 )
                      + ( ( idata0[ m + 23] & mask ) <<  8 )
                      + ( ( idata0[ m + 24] & mask ) <<  7 )
                      + ( ( idata0[ m + 25] & mask ) <<  6 )
                      + ( ( idata0[ m + 26] & mask ) <<  5 )
                      + ( ( idata0[ m + 27] & mask ) <<  4 )
                      + ( ( idata0[ m + 28] & mask ) <<  3 )
                      + ( ( idata0[ m + 29] & mask ) <<  2 )
                      + ( ( idata0[ m + 30] & mask ) <<  1)
                      + ( ( idata0[ m + 31] & mask )  ) )
        if (tmp >> (32-1)) & 1 == 1
            sign = -1
            packed0[i] = - ( Int( typemax(UInt32) + 1 ) - tmp)
        else
            sign = 1
            packed0[i] = tmp
        end
    end
    
    if ilen > ilen2 
        m = width * ilen2
        tmp = 0
        for j = 0:len_orig-m-1
            value = (idata0[ m+j ] & mask) << (width - j - 1)
            tmp = tmp | value
        end
        if (tmp >> (32-1)) & 1 == 1
            sign = -1
            packed0[ilen2] = - ( Int( typemax(UInt32) + 1 ) - tmp)
        else
            sign = 1
            packed0[ilen2] = tmp
        end
    end

    packed1 = OffsetArrays.no_offset_view(packed0)
    return packed1, ilen
end

function unpack_bits_from32( ilen, ipack1, nbit )
    bwidth = 32
    ipack0 = OffsetArray(ipack1, 0:length(ipack1)-1)
    if nbit == 1
        idata1 = zeros( Int64, ilen )
        idata_tmp, idata_size = unpack_bit( ilen, ipack0 )
        idata1[1:idata_size] = idata_tmp[1:idata_size]
    else
        imask = (1 << nbit ) - 1
#        nfull = Int32(ilen) >> 5
#        idata1 = zeros( Int64, (nfull-1)<<5+31+1 )
#        idata0 = OffsetArray(idata1, 0:((nfull-1)<<5 + 31))
        idata1 = zeros( Int64, ilen )
        idata0 = OffsetArray(idata1, 0:ilen-1 )
        for i = 0: ilen - 1
            i2 = i >> 5
            i3 = i & 31 

            ipos = nbit * i2 + ( ( nbit * i3 ) >> 5 )
            ioff = nbit + ( (nbit * i3) & 31)

            ival = ( ipack0[ipos] << (ioff - bwidth) ) & imask
            if ioff > bwidth
                if 2*bwidth > ioff
                    ival = ival | ( ( ipack0[ipos+1] >>> ( 2 * bwidth - ioff ) ) & imask )
                else
                    ival = ival | ( ( ipack0[ipos+1] << ( ioff-2*bwidth ) ) & imask )
                end
            end

            idata0[i] = ival
        end
        idata1 = OffsetArrays.no_offset_view(idata0)
    end
    return idata1
end

function pack_bits_into32( idata1, len_orig, nbit )

    if len_orig == 0
        println("Do nothing for 0-length input in pack_bits_into32.")
        exit(1)
    end
    
    idata0 = OffsetArray(idata1, 0:length(idata1)-1)

    bwidth = 32

    packed1 = zeros( Int32, len_orig )
    packed0 = OffsetArray(packed1, 0:len_orig-1)
    packed164 = zeros( Int64, len_orig )
    packed064 = OffsetArray(packed164, 0:len_orig - 1)

    if nbit == 1 
        packed1, ilen = pack_bit(idata1, len_orig )
        packed0 = OffsetArray(packed1, 0:length(packed1)-1)
    else
        mask = (1 << nbit) - 1
        ilen = ipack32len( len_orig, nbit )

        for i = 0:ilen-1
            packed064[ i ] = 0
        end
        
        ip = 0
        ioff = 0
        for i = 0:len_orig - 1
            if ioff > bwidth
                ioff = ioff - bwidth
                ip = ip + 1
            end
            
            ival = idata0[i] & mask
            packed064[ ip ] = packed064[ip] | ( ival << (bwidth - ioff - nbit) )
            if ioff + nbit > bwidth
                packed064[ip+1] = packed064[ip+1] | (( ival << (2 * bwidth - ioff - nbit ) ) & typemax(UInt32) )
            end

            ioff = ioff + nbit
        end

        for i = 0:len_orig-1
            if packed064[i] > typemax(Int32) && packed064[i] < typemax(UInt32)
                tmp0 = - (typemax(UInt32) + 1 - packed064[i])
                tmp = Int32(tmp0)
                packed0[i] = tmp
            elseif packed064[i] > typemax(Int32) && packed064[i] >= typemax(UInt32)
                tmp = Int64( typemax(UInt32) ) & packed064[i]
                if (tmp >> (32-1)) & 1 == 1
                    sign = -1
                    packed0[i] = - ( Int( typemax(UInt32) + 1 ) - tmp )
                else
                    sign = 1
                    packed0[i] = tmp
                end
            else
                tmp = typemax(Int32) & packed064[i]

                if (tmp >> (32-1)) & 1 == 1
                    sign = -1
                    packed0[i] = - ( Int( typemax(UInt32) + 1 ) - tmp )
                else
                    sign = 1
                    packed0[i] = tmp
                end
            end
        end

    end
    
    packed1 = OffsetArrays.no_offset_view(packed0)
    return packed1, ilen
end

function writeMR4( f, chead, nmr, gdata, miss )
    gmskd = zeros( Float32, nmr )
    imsk = zeros( Int32, nmr )
    iflag = zeros( Int64, nmr )
 
    nnn = Int32(0)
    for ij = 1:nmr
        if gdata[ij] != miss
            nnn = nnn + 1
            gmskd[ nnn ] = gdata[ ij ]
            iflag[ ij ] = 1
        else
            iflag[ ij ] = 0
        end
    end
    imsk, nmsk = pack_bits_into32( iflag, nmr, 1 )

    write( f, chead )
    outint = Int32(nnn)
    write( f, outint )
    outint = Int32.(imsk[1:nmsk])
    write( f, outint )
    arrayout = Float32.(gmskd[1:nnn])
    write( f, arrayout )
end

function writeMR8( f, chead, nmr, gdata, miss )
    gmskd = zeros( Float64, nmr )
    imsk = zeros( Int32, nmr )
    iflag = zeros( Int64, nmr )
 
    nnn = 0
    for ij = 1:nmr
        if gdata[ij] != miss
            nnn = nnn + 1
            gmskd[ nnn ] = gdata[ ij ]
            iflag[ ij ] = 1
        else
            iflag[ ij ] = 0
        end
    end
    imsk, nmsk = pack_bits_into32( iflag, nmr, 1 )

    write( f, chead )
    outint = Int32(nnn)
    write( f, outint )
    outint = Int32.(imsk[1:nmsk])
    write( f, outint )
    arrayout = Float64.(gmskd[1:nnn])
    write( f, arrayout )
end

#gfxscp( gdata[:, :, k], nxy, 1 << ibit, miss )
function gfxscp( gdatak, inum, nr, vmiss )
    dbl_max = 1.7e+308
    dbl_min = 2.3e-308

    imiss = nr - 1
    
    # << get offset and scale >>
    
    dmin = dbl_max
    dmax = - dbl_max
    dmin = minimum( gdatak[ gdatak .!= vmiss ] )
    dmax = maximum( gdatak[ gdatak .!= vmiss ] )
#    for ij = 1: inum
#        if gdatak[ij] != vmiss
#            dmin = min( gdatak[ ij ], dmin )
#            dmax = max( gdatak[ ij ], dmax )
#        end
#    end

    dx0  = 1 / max( imiss - 1, 1 )
    if dmax > dmin 
        dstp = dmax * dx0 - dmin * dx0
        dstp = max( dstp, dbl_min )
        dstp = min( dstp, dbl_max * dx0 )
        dxr = 1 / dstp
    else
        dstp = 0
        dxr = 0
    end
    if ! (dmin >= 0 || dmax < 0)
        if dmax == 0
            dmin = dmax
            dstp = - dstp
            dxr = - dxr
        else

            amin = abs( dmin )
            amax = abs( dmax )
            if amin < amax 
                xi = amin / amax
            else
                xi = amax / amin
            end
            if ! (xi < 1e-10)
                i0 = floor(Int,  (imiss - 1) / (1 + amax / amin) )

                if ! ( i0 == 0 || i0 == imiss - 1 )
                    
                    astp = real( amin / i0 )
                    if astp == 0 then
                        dstp = amin / i0
                    else
                        dstp = astp
                    end
                    
                    dstp = max( dstp, dbl_min )
                    dstp = min( dstp, dbl_max * dx0 )
                    dxr = 1 / dstp
                    dmin = - dstp * i0
                end
            end
        end
    end

    # << scaling >>
    gdatak1 = reshape( gdatak, inum )
    idata = zeros( Int64, inum )
    for ij = 1:inum
        if gdatak1[ij] != vmiss
            val = min( gdatak1[ ij ] - dmin, dbl_max )
            idata[ ij ] = min( round( Int, val * dxr ), imiss - 1 )
        else
            idata[ ij ] = imiss
        end
    end
    
    return idata, dmin, dstp
end
    
function writeURY( f, chead, nxy, nz, gdata, miss, fmt )

    dma = zeros( Float64, (2, nz) )
    ibit = 0
    ibit = parse( Int, fmt[4:5] )

    if ibit < 1
        ibit = 1
    end
    if ibit > 31
        ibit = 31
    end

    izlen = ipack32len( nxy, ibit )
    iz = zeros( Int32, nz * nxy )

    for k = 1:nz
        idata, dma[1, k], dma[2, k] = gfxscp( gdata[:, :, k], nxy, 1 << ibit, miss )
        iz[1+(k-1)*izlen:(k-1)*izlen+nxy], idummy = pack_bits_into32( idata, nxy, ibit )
    end
        
    write( f, chead )
    write( f, dma[ 1:2, 1:nz ]  )
    write( f, iz[ 1:izlen*nz ] )

end
    
function readMR4(f, nmr, miss)
    varT = fill(0.e0, nmr)
    nmsk = ipack32len( nmr, 1)
    ntmp = read(f, (Int32, 1) )
    nnn = ntmp[1]
    imsk = read(f, (Int32, nmsk) )
    gmskd = read(f, (Float32, nnn))
    iflag = unpack_bits_from32( nmr, imsk, 1 )
    n = 0
    for i = 1: nmr
        if iflag[ i ] == 1
            n = n + 1
            varT[ i ] = Float64( gmskd[ n ] )
        else
            varT[ i ] = miss
        end
    end
    return varT
end

function readMR8(f, nmr, miss)
    varT = fill(0.e0, nmr)
    nmsk = ipack32len( nmr, 1)
    ntmp = read(f, (Int32, 1) )
    nnn = ntmp[1]
    imsk = read(f, (Int32, nmsk) )
    gmskd = read(f, (Float64, nnn))
    iflag = unpack_bits_from32( nmr, imsk, 1 )
    n = 0
    for i = 1: nmr
        if iflag[ i ] == 1
            n = n + 1
            varT[ i ] = gmskd[ n ]
        else
            varT[ i ] = miss
        end
    end
    return varT
end

function readURY( f, nxy, nz, miss, dfmt )

    varT = fill(0.e0, nxy, nz)
    nbit = 0
    nbit = parse( Int, dfmt[4:5] )
    if nbit < 1
        nbit = 1
    elseif nbit > 31
        nbit = 31
    end

    izlen = ipack32len( nxy, nbit ) #  packed length for each z-level

    dma = read( f, (Float64, 2, nz) )
    iz = read( f, (Int32, izlen*nz) )

    imiss = ( 1 << nbit ) - 1

    for k = 1: nz
        idata = unpack_bits_from32( nxy, iz[ 1 + (k-1)*izlen:k*izlen], nbit )

        for i = 1: nxy
            if idata[i] != imiss
                varT[ i, k ] = dma[ 1, k ] + idata[ i ] * dma[ 2, k ]
            else
                varT[ i, k ] = miss
            end
        end
    end

    return varT
end


function readaxis( filename, direction )
    if isfile( "./" * filename )
        AxisFile = "./" * filename
    elseif isfile( ENV["GTAXDIR"] * "/" * filename )
        AxisFile = ENV["GTAXDIR"] * "/" * filename
    else
        println("Axis file ", filename, " not found in GTAXDIR.")
        exit(1)
    end
    fw = opengtool( AxisFile, "r" )
    cheadW = read(fw, (FString{16}, 64))
    dfmtW = lstrip( trimstring( cheadW[38] ), ' ')
    if dfmtW != "UR4" && dfmtW != "UR8" && dfmtW != "MR4" && dfmtW != "MR8" && dfmtW[1:3] != "URY"
        println("Format ", dfmtW,  " is not supported.\n")
        exit(1)
    end
    sc = lstrip( trimstring( cheadW[31] ), ' ')
    if direction == "lon"
        nwread = parse( Int64, sc )
        nw = nwread - 1
    else
        nwread = parse( Int64, sc )
        nw = nwread
    end
    if dfmtW == "UR4"
        axW = read(fw, (Float32, nw))
    elseif dfmtW == "UR8"
        axW = read(fw, (Float64, nw))
    elseif dfmtW == "MR4"
        axW = readMR4(fw, nw, -999.)
    elseif dfmtW == "MR8"
        axW = readMR8(fw, nw, -999.)
    elseif dfmtW[1:3] == "URY"
        axW = readURY(fw, nw, -999.)
    end
    close( fw )
    return axW, nw
end

function readchead( f )
    chead = read(f, (FString{16}, 64) )
    return chead
end

function readgtool( filename, ntin::Int64=0 )

#    filename = dir * "/" * varname

    if ntin == 0
        ngtsString = read( pipeline( `ngtstat $filename` , `tail -n1` ), String )
        cnt = split( ngtsString, r"\s+" )[2]
        nt = parse( Int64, cnt )
    else
        nt = ntin
    end

    f = opengtool( filename, "r" )
    chead = read(f, (FString{16}, 64) )

    dfmt = lstrip( trimstring( chead[38] ), ' ')
    if dfmt != "UR4" && dfmt != "UR8" && dfmt != "MR4" && dfmt != "MR8" && dfmt[1:3] != "URY"
        println("Format ", dfmt,  " is not supported.\n")
        exit(1)
    end

    # X
    sc = lstrip( trimstring( chead[29] ), ' ')
    filename = "GTAXLOC." * sc
    lon, nx = readaxis( filename, "lon" )

    # Y
    sc = lstrip( trimstring( chead[32] ), ' ')
    filename = "GTAXLOC." * sc
    lat, ny = readaxis( filename, "lat" )

    # Z
    sc = lstrip( trimstring( chead[35] ), ' ')
    filename = "GTAXLOC." * sc
    dep, nz = readaxis( filename, "dep" )

    unit = lstrip( trimstring( chead[16] ), ' ')
    sc = lstrip( trimstring( chead[39] ), ' ')
    miss = parse( Float64, sc)

    if nz == 1
        varXY = fill(0.e0, nx, ny)
        varXYT = fill(0.e0, nx, ny, nt)
    else
        varXYZ = fill(0.e0, nx, ny, nz)
        varXYZT = fill(0.e0, nx, ny, nz, nt)
    end
    tarray = fill(0.e0, nt)
    
    rewind( f )

    for it=1:nt
        chead = read(f, (FString{16}, 64))
        sc = lstrip( trimstring( chead[25] ), ' ')
        tarray[it] = parse( Int64, sc )
        if dfmt == "UR4"
            if nz == 1
                varXY = read(f, (Float32, (nx, ny)))
            else
                varXYZ = read(f, (Float32, (nx, ny, nz) ) )
            end
        elseif dfmt == "UR8"
            if nz == 1
                varXY = read(f, (Float64, (nx, ny)))
            else
                varXYZ = read(f, (Float64, (nx, ny, nz) ) )
            end
        elseif dfmt[1:3] == "URY"
            nxy = nx*ny
            varT1 = readURY(f, nxy, nz, miss, dfmt)
            if nz == 1
                varXY = collect( reshape( varT1, nx, ny ))
            else
                varXYZ = collect( reshape( varT1, nx, ny, nz ))
            end
        elseif dfmt == "MR4"
            nmr = nx*ny*nz
            varT1 = readMR4(f, nmr, miss)
            if nz == 1
                varXY = collect( reshape( varT1, nx, ny ))
            else
                varXYZ = collect( reshape( varT1, nx, ny, nz ))
            end
        elseif dfmt == "MR8"
            nmr = nx*ny*nz
            varT1 = readMR8(f, nmr, miss)
            if nz == 1
                varXY = collect( reshape( varT1, nx, ny ))
            else
                varXYZ = collect( reshape( varT1, nx, ny, nz ))
            end
        end
        if nz == 1
            varXYT[:, :, it]   = varXY[:, :]
        else
            varXYZT[:, :, :, it]   = varXYZ[:, :, :]
        end
    end

    if nz==1
        varXYT = replace!( varXYT, miss => NaN )
        varReturn = varXYT
    else
        varXYZT = replace!( varXYZT, miss => NaN )
        varReturn = varXYZT
    end

    close( f )
    return nx, ny, nz, nt, lon, lat, dep, tarray, varReturn

end 

function writegtool( f, chead, fmt, arraytmp )

    sc = lstrip( trimstring( chead[39] ), ' ')
    miss = parse( Float64, sc)
    arraytmp = replace!( arraytmp, NaN => miss )
    if fmt == "UR4" 
        chead[38]=@sprintf("%-16s", "UR4")
        arrayout = Float32.(arraytmp)
        write(f, chead)
        write(f, arrayout)
    elseif fmt == "UR8"
        chead[38]=@sprintf("%-16s", "UR8")
        arrayout = Float64.(arraytmp)
        write(f, chead)
        write(f, arrayout)
    elseif fmt == "MR4"
        chead[38]=@sprintf("%-16s", "MR4")
        nmr = length( arraytmp )
        writeMR4( f, chead, nmr, arraytmp, miss )
    elseif fmt == "MR8"
        chead[38]=@sprintf("%-16s", "MR8")
        nmr = length( arraytmp )
        writeMR8( f, chead, nmr, arraytmp, miss )
    elseif fmt[1:3] == "URY"
        chead[38]=@sprintf("%-16s", fmt)
        nxy = size( arraytmp, 1 ) * size( arraytmp, 2)
        nz = size( arraytmp, 3 )
        writeURY( f, chead, nxy, nz, arraytmp, miss, fmt )
    else
        println("Output data format ", fmt, " is not supported.")
        exit(1)
    end

end 

function opengtool( filename, mode )
    f = FortranFile( filename, mode, access="sequential", convert="big-endian" )
end

function closegtool( f )
    close( f )
end

function gtooltime2datetime( time )
    # Gtoole time is in hour since 0000-01-01T00:00:00
    # Julia time is in day since 0000-12-31T00:00:00 (Yr 0000 has 366 days)
    date = rata2datetime(floor(Int64, (time-365*24)/24)) + Hour(mod((time-365*24)/24, 1)*24)
    return date
end

end # module gtool

#######
