using CSV
using DataFrames
using Plots


"""
    read_img(img_file::String, bigendian=true)

Read an image from a binary file.
Returns the expected shape and pixels in all
dimensions from the file header as well as the
image info in a nD array.
"""
function read_img(img_file::String, bigendian=true)

    bin_info = read(img_file)
    metadata = bin_info[1:18]
    img_bin  = bin_info[19:end]

    img_pix  = reshape(metadata[1:6], 2, :)
    img_pix  = reinterpret(reshape, Int16  , img_pix )
    img_size = reshape(metadata[7:end], 4, :)
    img_size = reinterpret(reshape, Float32, img_size)

    img_data = reshape(img_bin, 4, :)
    img_data = reinterpret(reshape, Float32, img_data)
    if bigendian
        img_pix  = ntoh.(img_pix )
        img_size = ntoh.(img_size)
        img_data = ntoh.(img_data)
    end

    return img_pix, img_size, reshape(img_data, Int.(img_pix)...)
end


"""
    write_img(img_file::String,
              npixel::Tuple{Int64, Int64, Int64},
              img_size::Tuple{Float32, Float32, Float32},
              img_data::Array{Float32, 3})
Write array as a binary image file.
"""
function write_img(img_file::String                          ,
                   npixel  ::Tuple{UInt16 , UInt16 , UInt16 },
                   img_size::Tuple{Float32, Float32, Float32},
                   img_data::Array{Float32, 3}               )
    open(img_file, "w") do img_out
        write.(Ref(img_out), hton.(npixel))
        write.(Ref(img_out), hton.(img_size))
        write.(Ref(img_out), hton.(img_data))
    end;
end


function split_parse(row::String)::Vector{Float32}
    parse.(Float32, split(row))
end


"""
    read_foms(file_name::String)::DataFrame

Read a single image file figure of merit output.
"""
function read_foms(file_name::String)::DataFrame
    raw_input = CSV.read(file_name, DataFrame)

    #Header read as single String
    all_names = names(raw_input)[1]
    col_names = Symbol.(lstrip.(rstrip.(split(all_names, "   "))))
    return combine(raw_input, Symbol(all_names) => ByRow(split_parse) => col_names)
end


"""
    read_allfoms(file_name::String)::DataFrame

Read a file with figures of merit from multiple
iterations of mlem to a DataFrame.
"""
function read_allfoms(file_name::String)::DataFrame
    (col_types, diams, vals) = open(file_name) do fin
        fom_vals  = Vector{Float32}[]
        fom_types = String[]
        diams     = Float32[]
        for (i, line) in enumerate(eachline(fin))
            if i == 1
                fom_types = split(line)
            elseif i == 3
                diams = unique(split(line))
            elseif !in((2, 4))(i)
                push!(fom_vals, parse.(Float32, split(line)))
            end
        end
        return fom_types, diams, fom_vals
      end
      # Some gymnastics to get things in the right order.
      col_types = [col_types[1], col_types[2]*" "*col_types[3], col_types[4]]
      col_names = vcat("iteration", [typ*diam for typ in col_types for diam in diams])
      return DataFrame(transpose(hcat(vals...)), col_names)
end

# Version with CSV that sometimes doesn't work!!
# function read_allfoms(file_name::String)::DataFrame
#     raw_input = CSV.read(file_name, DataFrame, header=[1, 2], skipto=4)

#     function stitch_colnames(in_header::String)::Vector{Symbol}
#         nms      = unique(lstrip.(rstrip.(split(in_header, "  "))))
#         is_fom   = isnothing.(tryparse.(Float32, nms[2:end]))
#         foms     = nms[2:end][is_fom]
#         rads     = repeat(nms[2:end][.!is_fom], 1, length(foms))
#         all_cols = [.*(fom, vrad) for (fom, vrad) in zip(foms, eachcol(rads))]
#         return Symbol.(vcat("iteration", all_cols...))
#     end

#     header_str = names(raw_input)[1]
#     col_names  = stitch_colnames(header_str)
#     return combine(raw_input, Symbol(header_str) => ByRow(split_parse) => col_names)
# end
