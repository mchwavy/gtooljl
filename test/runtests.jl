using gtooljl
using Test


# test1

@testset "test1" begin
    
    nn = 9
    idata = zeros( Int64, nn )
    packed = zeros( Int32, nn )

    nbit = 12
    idata[1] = 273  # Z'111'
    idata[2] = 546  # Z'222'
    idata[3] = 819  # Z'333'
    idata[4] = 1092 # Z'444'
    idata[5] = 1365 # Z'555'
    idata[6] = 1638 # Z'666'
    idata[7] = 1911 # Z'777'
    idata[8] = 2184 # Z'888'
    idata[9] = 1911 # Z'777'

    packed, ilen = pack_bits_into32( idata, 9, nbit )

    @test ilen ==  4
    @test packed[1] == 286401075
    @test packed[2] == 876893526
    @test packed[3] == 1719105672
    @test packed[4] == 2003828736

    #     12-bit test(2)
    idata[1] = 4095  # Z'fff'
    idata[2] = 3822  # Z'eee'
    idata[3] = 3549  # Z'ddd'
    idata[4] = 3276  # Z'ccc'
    idata[5] = 3003  # Z'bbb'
    idata[6] = 2730  # Z'aaa'
    idata[7] = 2457  # Z'999'
    idata[8] = 2184  # Z'888'
    idata[9] = 4095  # Z'fff'

    packed, ilen = pack_bits_into32( idata, 9, nbit )

    @test ilen == 4
    @test packed[1] == -69923
    @test packed[2] == -590562374
    @test packed[3] == -1432774520
    @test packed[4] == -1048576

    #     12-bit test(3)
    idata[1] = 18    # Z'012'
    idata[2] = 837   # Z'345'
    idata[3] = 1656  # Z'678'
    idata[4] = 2475  # Z'9ab'
    idata[5] = 3294  # Z'cde'
    idata[6] = 3841  # Z'f01'

    packed, ilen = pack_bits_into32( idata, 6, nbit )

    @test ilen == 3
    @test packed[1] == 19088743
    @test packed[2] == -1985229329
    @test packed[3] == 16777216

    # test for nbit==1

    nbit = 1
    idata[1] = 0
    idata[2] = 1
    idata[3] = 1
    idata[4] = 1
    idata[5] = 0
    idata[6] = 1
    idata[7] = 1
    idata[8] = 0
    idata[9] = 1

    packed, ilen = pack_bits_into32( idata, 8, nbit )

    @test ilen == 1
    @test packed[1] == 1979711488

    packed, ilen = pack_bits_into32( idata, 9, nbit )
    @test ilen == 1
    @test packed[1] == 1988100096

end
#######

# test2()
@testset "test2" begin
    
    # test for nbit = 1
    nbit = 1

    nn = 5000

    idata00 = zeros(Int64, nn)

    imask = (1 <<  nbit) - 1
    for i = 1: nn
        idata00[i] = (i - 1) & imask
    end

    idata01, ilen = pack_bits_into32( idata00, nn, nbit )

    idata02 =  unpack_bits_from32( nn, idata01, nbit )

    for i = 1:nn
        @test idata00[i] == idata02[i]
    end
end
  
####

# function test3()

#######
###
### This doesn't work for Julia.
###
#######

#    idata00 = zeros(Int64, 1)
#    idata01 = zeros(Int32, 1)
#
#    idata00[1] = 999           # dummy
#    idata01[1] = -1
#    idata01, ilen = pack_bits_into32( idata00, 0, 12 )    # 0-length input => do nothing

#    assert(ilen,  0, "test3")
#    assert(idata01[1], -1, "test3 ")

#    println("Done.")

#end

#######

# test4()
@testset "test4" begin
    
    nn = 329

    nbit = 1

    # test(1)
    idata = fill(1, nn)
    packed = fill(-1, nn)
    packed, ilen = pack_bits_into32( idata, nn, nbit )
    
    @test ilen == 11
    for i = 1:10
        @test packed[i] ==  -1
    end
    @test packed[11] == -8388608

    # test(2)
    idata = fill(0, nn)
    packed = fill(-1, nn)

    packed, ilen = pack_bits_into32( idata, nn, nbit )

    @test ilen == 11
    for i = 1:11
        @test packed[i] == 0
    end

    # test(3)
    idata = fill(0, nn)
    for i = 32 : 32 : nn
        idata[i] = 1
    end 
    packed = fill(-1, nn)
    packed, ilen = pack_bits_into32( idata, nn, nbit )

    @test ilen == 11
    for i = 1:10
        @test packed[i] == 1
    end
    @test packed[11] == 0

end

#######

@testset "test5" begin
    # test5()

    ipack = zeros(Int64, 4)
    idata = zeros(Int64, 100)

    ipack[1]=1
    for i=1:32
        idata_tmp, idata_size = unpack_bit( 32, ipack )
        idata[1:idata_size] = idata_tmp[1:idata_size]
        for m = 1:32
            if m + i == 33
                @test idata[m] == 1
            else
                @test idata[m] == 0
            end
        end
        ipack[1] = 2 * ipack[1]
    end

    idata = zeros(Int64, 100)
    for i = 1:100
        idata[i] = -999
    end

    ipack[1] = 1
    ipack[2] = -1
    ipack[3] = 3
    ipack[4] = -1

    idata_tmp, idata_size = unpack_bit( 99, ipack )
    idata[1:idata_size] = idata_tmp[1:idata_size]

    @test idata[100] == -999
    for i = 1:31
        @test idata[i] == 0
    end
    
    for i = 32:64
        @test idata[i] == 1
    end 
    for i = 65:94
        @test idata[i] == 0
    end
    for i = 95:99
        @test idata[i] ==  1
    end

end

#@testset "gtool.jl" begin
    
#    @test testgtool.test1 == 0

#end
