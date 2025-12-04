using NCDatasets
using Images
using FileIO
using Statistics
using Printf
using Glob

# Start timer
start_time = time()

# Configuration
file_pattern = "/mnt/nfs_Vol3T/received/hvs-2/E2H-MTG-1/2025/11/29/W_XX*0072_00[0-9][0-9].nc"
files = sort(glob(file_pattern))
pop!(files)  # Remove last file if needed
output_file = "fci_composite.png"
resolution = 11136  # 1km resolution

# Initialize output arrays
composite_R = zeros(Float32, resolution, resolution)
composite_G = zeros(Float32, resolution, resolution)
composite_B = zeros(Float32, resolution, resolution)
total_subsolar_lat = zeros(Float32, resolution, resolution)
total_subsolar_lon = zeros(Float32, resolution, resolution)

# Process each file
for file_path in files
    println("Processing: $file_path")
    
    try
        ds = NCDataset(file_path, "r")
        
        index_offset = ds["index_offset"][:][1]
        println("index_offset: $index_offset")
        
        # Read radiance data
        radiance_B = ds["data"]["vis_04"]["measured"]["effective_radiance"][:, :]
        radiance_G = ds["data"]["vis_05"]["measured"]["effective_radiance"][:, :]
        radiance_R = ds["data"]["vis_06"]["measured"]["effective_radiance"][:, :]
        
        # Get index_map
        index_map = UInt16.(ds["data"]["vis_04"]["measured"]["index_map"][:, :])
        
        # Read position information
        start_row = Int(ds["data"]["vis_06"]["measured"]["start_position_row"][:][1])
        start_col = Int(ds["data"]["vis_06"]["measured"]["start_position_column"][:][1])
        end_row = Int(ds["data"]["vis_06"]["measured"]["end_position_row"][:][1])
        end_col = Int(ds["data"]["vis_06"]["measured"]["end_position_column"][:][1])
        
        # Get geometric parameter vectors
        subsolar_lat = ds["state"]["celestial"]["subsolar_latitude"][:]
        subsolar_lon = ds["state"]["celestial"]["subsolar_longitude"][:]
        
        # Place segment in composite (convert to 0-based indexing)
        composite_R[start_row:end_row-1, start_col:end_col-1] .= radiance_R
        composite_G[start_row:end_row-1, start_col:end_col-1] .= radiance_G
        composite_B[start_row:end_row-1, start_col:end_col-1] .= radiance_B
        
        # Adjust index_map by subtracting index_offset
        index_map_adjusted = similar(index_map)
        for i in eachindex(index_map)
            if index_map[i] == 65535
                index_map_adjusted[i] = 65535
            else
                index_map_adjusted[i] = index_map[i] - UInt16(index_offset)
            end
        end
        
        # Clamp indices and extract values
        indices_flat = vec(index_map_adjusted)
        valid_indices = clamp.(indices_flat, 1, length(subsolar_lat))
        
        # Handle fill values
        subsolar_lat_values = similar(index_map_adjusted, Float32)
        subsolar_lon_values = similar(index_map_adjusted, Float32)
        
        for i in eachindex(index_map_adjusted)
            if index_map_adjusted[i] == 65535
                subsolar_lat_values[i] = 0.0f0
                subsolar_lon_values[i] = 0.0f0
            else
                subsolar_lat_values[i] = subsolar_lat[valid_indices[i]]
                subsolar_lon_values[i] = subsolar_lon[valid_indices[i]]
            end
        end
        
        total_subsolar_lat[start_row:end_row-1, start_col:end_col-1] .= subsolar_lat_values
        total_subsolar_lon[start_row:end_row-1, start_col:end_col-1] .= subsolar_lon_values
        
        close(ds)
        
    catch e
        println("Error processing $file_path: $e")
        continue
    end
end

# Normalization
println("composite_R raw range min: $(minimum(composite_R))")
println("composite_R raw range max: $(maximum(composite_R))")
println("composite_R unique values: $(length(unique(composite_R)))")

valid_data_R = composite_R[composite_R .!= 65535.0]
valid_data_G = composite_G[composite_G .!= 65535.0]
valid_data_B = composite_B[composite_B .!= 65535.0]

println("Valid data R len: $(length(valid_data_R))")

# Normalize each channel
function normalize_channel(data)
    valid_data = data[data .!= 65535.0]
    if length(valid_data) > 0
        vmin = quantile(valid_data, 0.0)
        vmax = quantile(valid_data, 1.0)
        normalized = similar(data)
        for i in eachindex(data)
            if data[i] == 65535.0
                normalized[i] = 0.0f0
            else
                normalized[i] = clamp((data[i] - vmin) / (vmax - vmin) * 255, 0, 255)
            end
        end
        return normalized
    else
        return data
    end
end

composite_normalized_R = normalize_channel(composite_R)
composite_normalized_G = normalize_channel(composite_G)
composite_normalized_B = normalize_channel(composite_B)

# Convert to uint8
composite_uint8_R = UInt8.(composite_normalized_R)
composite_uint8_G = UInt8.(composite_normalized_G)
composite_uint8_B = UInt8.(composite_normalized_B)

# Create alpha channel
alpha_channel = [composite_R[i, j] == 65535.0 ? 0x00 : 0xff for i in axes(composite_R, 1), j in axes(composite_R, 2)]

# Create RGBA image
rgb_img = RGB.(composite_uint8_R, composite_uint8_G, composite_uint8_B)
rgba_img = RGBA.(rgb_img, alpha_channel)

# Flip and save
rgba_flipped = reverse(rgba_img; dims=1)
save(output_file, rgba_flipped)

println("Composite image saved to: $output_file")
println("Image shape: $(size(composite_uint8_R))")
println("Value range: $(minimum(composite_uint8_R)) - $(maximum(composite_uint8_R))")

end_time = time()
elapsed_time = end_time - start_time
@printf "Script execution time: %.2f seconds\n" elapsed_time
@printf "Script execution time: %.2f minutes\n" elapsed_time / 60