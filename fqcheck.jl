#!/usr/bin/julia
using Libz

# fqcheck
function fqcheck(infile::ASCIIString, read_length::Int64)
    if endswith(infile, ".gz")
        hand = ZlibInflateInputStream(open(infile, "r"), reset_on_end=true)
    else
        hand = open(infile)
    end
    outArray::Array{Int64,2} = zero(Array{Int64,2}(read_length,48))
    orderDict::Dict{Char,Int64} = Dict{Char,Int64}('A'=>1, 'T'=>2, 'C'=>3, 'G'=>4, 'N'=>5)
    for line1::ASCIIString in eachline(hand)
        line2::ASCIIString = readline(hand)[1:end-1]
        _::ASCIIString, line4::ASCIIString = readline(hand), readline(hand)
        len::Int64 = length(line2)
        for i in 1:len
            base::Char = line2[i]
            baseInArray::Int64 = orderDict[base]
            outArray[i, baseInArray]::Int64 += 1
            # phred 33
            qual::Int64 = Int(line4[i])-33
            qualInArray::Int64 = qual+6
            outArray[i, qualInArray]::Int64 += 1
        end
    end
    return outArray
end

function ration(a::Int64,t::Int64,c::Int64,g::Int64,n::Int64)
    s::Int64 = a+t+c+g+n
    ra::Float64 = a/s
    rt::Float64 = t/s
    rc::Float64 = c/s
    rg::Float64 = g/s
    rn::Float64 = n/s
    return(ra, rt, rc, rg, rn)
end

function printArray(inArray::Array{Int64,2}, outfile::ASCIIString)
    f::IOStream = open(outfile, "w")
    all_base = sum(inArray[1:end, 6:end])
    q20 = sum(inArray[1:end, 26:end])/all_base
    q30 = sum(inArray[1:end, 36:end])/all_base
    print(f, "#base: ", all_base, '\n')
    print(f, "#Q20: ", q20*100, '\n')
    print(f, "#Q30: ", q30*100, '\n')

    print(f, "pos A T C G N")
    for qual in 0:42
        print(f, " ", qual)
    end
    println(f)
    reads_len::Int64 = size(inArray, 1)
    for pos in 1:reads_len
        A::Int64 = inArray[pos,1]
        T::Int64 = inArray[pos,2]
        C::Int64 = inArray[pos,3]
        G::Int64 = inArray[pos,4]
        N::Int64 = inArray[pos,5]
        if A+T+C+G+N == 0 ## 
            return 0
        end
        ra::Float64, rt::Float64, rc::Float64, rg::Float64, rn::Float64 = ration(A, T, C, G, N)
        print(f, pos, " ", ra, " ", rt, " ", rc, " ", rg, " ",rn)
        for q_plus in 6:48
            qual_count::Int64 = resultArray[pos,q_plus]
            print(f, " ", qual_count)
        end
        println(f)
    end
    close(f)
end

# main
    if length(ARGS) != 1
        info("\n\tUsage:\tjulia <reads.fastq.gz>")
        exit(-1)
    end
    infile = abspath(ARGS[1])
    name = split(infile, '/')[end]
    name1 = split(name, '.')[1]
    outfile = string(name1, ".fqcheck")
    reads_len = 151 # the longest 
    resultArray = fqcheck(infile, reads_len)
    printArray(resultArray, outfile)
