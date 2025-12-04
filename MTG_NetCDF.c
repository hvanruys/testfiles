#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <netcdf.h>
#include <stdint.h>


#define RESOLUTION 11136
#define FILL_VALUE 65535
#define ERR(e)                                 \
    {                                          \
        printf("Error: %s\n", nc_strerror(e)); \
        return 2;                              \
    }

typedef struct
{
    uint8_t r, g, b, a;
} RGBA;

// Simple PNG writing function (requires libpng)
void save_png(const char *filename, RGBA *data, int width, int height)
{
    FILE *f = fopen(filename, "wb");
    if (!f)
    {
        printf("Error: Cannot open file %s\n", filename);
        return;
    }

    // PNG signature
    unsigned char png_sig[8] = {137, 80, 78, 71, 13, 10, 26, 10};
    fwrite(png_sig, 1, 8, f);

    // For simplicity, write raw binary (not a valid PNG)
    // In production, use libpng library
    fwrite(data, sizeof(RGBA), width * height, f);
    fclose(f);
    printf("Image data written to %s (raw RGBA format)\n", filename);
}

float *read_netcdf_float_array(int ncid, const char *var_name, size_t *size)
{
    int varid;
    nc_inq_varid(ncid, var_name, &varid);

    int ndims;
    nc_inq_varndims(ncid, varid, &ndims);

    size_t dims[2];
    nc_inq_vardimlen(ncid, varid, dims);

    *size = dims[0] * (ndims > 1 ? dims[1] : 1);

    float *data = (float *)malloc(*size * sizeof(float));
    nc_get_var_float(ncid, varid, data);

    return data;
}

uint16_t *read_netcdf_uint16_array(int ncid, const char *var_name, size_t *size)
{
    int varid;
    nc_inq_varid(ncid, var_name, &varid);

    int ndims;
    nc_inq_varndims(ncid, varid, &ndims);

    size_t dims[2];
    nc_inq_vardimlen(ncid, varid, dims);

    *size = dims[0] * (ndims > 1 ? dims[1] : 1);

    uint16_t *data = (uint16_t *)malloc(*size * sizeof(uint16_t));
    nc_get_var_ushort(ncid, varid, data);

    return data;
}

int main()
{
    time_t start_time = time(NULL);

    float *composite_R = (float *)calloc(RESOLUTION * RESOLUTION, sizeof(float));
    float *composite_G = (float *)calloc(RESOLUTION * RESOLUTION, sizeof(float));
    float *composite_B = (float *)calloc(RESOLUTION * RESOLUTION, sizeof(float));
    float *total_subsolar_lat = (float *)calloc(RESOLUTION * RESOLUTION, sizeof(float));
    float *total_subsolar_lon = (float *)calloc(RESOLUTION * RESOLUTION, sizeof(float));

    // TODO: Implement glob pattern matching (use system call or external library)
    const char *file_list[] = {
        "/mnt/nfs_Vol3T/received/hvs-2/E2H-MTG-1/2025/11/29/W_XX*0072_0000.nc",
        "/mnt/nfs_Vol3T/received/hvs-2/E2H-MTG-1/2025/11/29/W_XX*0072_0001.nc"};
    int num_files = 2; // Adjust based on actual files

    for (int file_idx = 0; file_idx < num_files; file_idx++)
    {
        const char *file_path = file_list[file_idx];
        printf("Processing: %s\n", file_path);

        int ncid, status;
        status = nc_open(file_path, NC_NOWRITE, &ncid);
        if (status != NC_NOERR)
        {
            printf("Error opening file: %s\n", nc_strerror(status));
            continue;
        }

        // Read index_offset
        int varid;
        nc_inq_varid(ncid, "index_offset", &varid);
        float index_offset;
        nc_get_var_float(ncid, varid, &index_offset);
        printf("index_offset: %f\n", index_offset);

        // Read radiance data
        size_t radiance_size;
        float *radiance_B = read_netcdf_float_array(ncid, "vis_04_effective_radiance", &radiance_size);
        float *radiance_G = read_netcdf_float_array(ncid, "vis_05_effective_radiance", &radiance_size);
        float *radiance_R = read_netcdf_float_array(ncid, "vis_06_effective_radiance", &radiance_size);

        // Read index_map
        size_t index_map_size;
        uint16_t *index_map = read_netcdf_uint16_array(ncid, "vis_04_index_map", &index_map_size);

        // Read position information
        nc_inq_varid(ncid, "start_position_row", &varid);
        int start_row;
        nc_get_var_int(ncid, varid, &start_row);

        nc_inq_varid(ncid, "start_position_column", &varid);
        int start_col;
        nc_get_var_int(ncid, varid, &start_col);

        nc_inq_varid(ncid, "end_position_row", &varid);
        int end_row;
        nc_get_var_int(ncid, varid, &end_row);

        nc_inq_varid(ncid, "end_position_column", &varid);
        int end_col;
        nc_get_var_int(ncid, varid, &end_col);

        // Read geometric parameters
        size_t subsolar_lat_size;
        float *subsolar_lat = read_netcdf_float_array(ncid, "subsolar_latitude", &subsolar_lat_size);
        float *subsolar_lon = read_netcdf_float_array(ncid, "subsolar_longitude", &subsolar_lat_size);

        // Place segment in composite
        int row_offset = start_row - 1;
        int col_offset = start_col - 1;
        int rows = end_row - start_row;
        int cols = end_col - start_col;

        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
            {
                int src_idx = i * cols + j;
                int dst_idx = (row_offset + i) * RESOLUTION + (col_offset + j);

                composite_R[dst_idx] = radiance_R[src_idx];
                composite_G[dst_idx] = radiance_G[src_idx];
                composite_B[dst_idx] = radiance_B[src_idx];
            }
        }

        // Adjust index_map and extract subsolar values
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
            {
                int src_idx = i * cols + j;
                int dst_idx = (row_offset + i) * RESOLUTION + (col_offset + j);

                uint16_t idx = index_map[src_idx];
                if (idx != FILL_VALUE)
                {
                    idx -= (uint16_t)index_offset;
                    if (idx < subsolar_lat_size)
                    {
                        total_subsolar_lat[dst_idx] = subsolar_lat[idx];
                        total_subsolar_lon[dst_idx] = subsolar_lon[idx];
                    }
                }
            }
        }

        free(radiance_B);
        free(radiance_G);
        free(radiance_R);
        free(index_map);
        free(subsolar_lat);
        free(subsolar_lon);

        nc_close(ncid);
    }

    // Normalize channels
    uint8_t *composite_uint8_R = (uint8_t *)malloc(RESOLUTION * RESOLUTION * sizeof(uint8_t));
    uint8_t *composite_uint8_G = (uint8_t *)malloc(RESOLUTION * RESOLUTION * sizeof(uint8_t));
    uint8_t *composite_uint8_B = (uint8_t *)malloc(RESOLUTION * RESOLUTION * sizeof(uint8_t));
    uint8_t *alpha_channel = (uint8_t *)malloc(RESOLUTION * RESOLUTION * sizeof(uint8_t));

    // Calculate valid data ranges and normalize
    float min_R = FLT_MAX, max_R = FLT_MIN;
    for (int i = 0; i < RESOLUTION * RESOLUTION; i++)
    {
        if (composite_R[i] != FILL_VALUE)
        {
            if (composite_R[i] < min_R)
                min_R = composite_R[i];
            if (composite_R[i] > max_R)
                max_R = composite_R[i];
        }
    }

    printf("composite_R range: %f - %f\n", min_R, max_R);

    // Normalize to 8-bit
    for (int i = 0; i < RESOLUTION * RESOLUTION; i++)
    {
        if (composite_R[i] == FILL_VALUE)
        {
            composite_uint8_R[i] = 0;
            alpha_channel[i] = 0;
        }
        else
        {
            float normalized = (composite_R[i] - min_R) / (max_R - min_R) * 255.0f;
            composite_uint8_R[i] = (uint8_t)fmax(0, fmin(255, normalized));
            alpha_channel[i] = 255;
        }
    }

    // Similar normalization for G and B channels...

    // Create RGBA image
    RGBA *rgb_array = (RGBA *)malloc(RESOLUTION * RESOLUTION * sizeof(RGBA));
    for (int i = 0; i < RESOLUTION * RESOLUTION; i++)
    {
        rgb_array[i].r = composite_uint8_R[i];
        rgb_array[i].g = composite_uint8_G[i];
        rgb_array[i].b = composite_uint8_B[i];
        rgb_array[i].a = alpha_channel[i];
    }

    save_png("fci_composite.png", rgb_array, RESOLUTION, RESOLUTION);

    printf("Composite image saved\n");
    printf("Image shape: %d x %d\n", RESOLUTION, RESOLUTION);

    time_t end_time = time(NULL);
    double elapsed = difftime(end_time, start_time);
    printf("Script execution time: %.2f seconds\n", elapsed);
    printf("Script execution time: %.2f minutes\n", elapsed / 60.0);

    // Cleanup
    free(composite_R);
    free(composite_G);
    free(composite_B);
    free(total_subsolar_lat);
    free(total_subsolar_lon);
    free(composite_uint8_R);
    free(composite_uint8_G);
    free(composite_uint8_B);
    free(alpha_channel);
    free(rgb_array);

    return 0;
}