import numpy as np

def fci_latlon_to_grid(lat, lon, ssd=1.0):
    """
    Convert latitude/longitude to FCI Level 1c reference grid row/column.
    
    Parameters:
    -----------
    lat : float or array-like
        Latitude in degrees
    lon : float or array-like
        Longitude in degrees
    ssd : float
        Spatial Sampling Distance in km (0.5, 1.0, or 2.0)
    
    Returns:
    --------
    row, column : float or array
        Row and column positions in the reference grid (1-indexed, can be fractional)
        Returns NaN if the point is not visible from the satellite
    
    Notes:
    ------
    Based on MTG FCI L1 Product User Guide sections 5.2 and 5.3
    This is the forward projection from geographic to satellite viewing angles.
    """
    
    # Grid parameters for different SSDs
    grid_params = {
        0.5: {
            'lambda_0': 8.9142405037 * np.pi / 180,
            'phi_0': -8.9142405037 * np.pi / 180,
            'azimuth_sampling': 0.000800524494 * np.pi / 180,
            'elevation_sampling': 0.000800524494 * np.pi / 180
        },
        1.0: {
            'lambda_0': 8.9138402398 * np.pi / 180,
            'phi_0': -8.9138402398 * np.pi / 180,
            'azimuth_sampling': 0.001601048988 * np.pi / 180,
            'elevation_sampling': 0.001601048988 * np.pi / 180
        },
        2.0: {
            'lambda_0': 8.9130397083 * np.pi / 180,
            'phi_0': -8.9130397083 * np.pi / 180,
            'azimuth_sampling': 0.003202097973 * np.pi / 180,
            'elevation_sampling': 0.003202097973 * np.pi / 180
        }
    }
    
    if ssd not in grid_params:
        raise ValueError("SSD must be 0.5, 1.0, or 2.0 km")
    
    params = grid_params[ssd]
    
    # Earth ellipsoid parameters
    r_eq = 6378.137  # km
    r_pol = r_eq * (1 - 1/298.257223563)  # km
    h = 35786.4 + r_eq  # km - geostationary radius
    
    # Sub-satellite point longitude
    lambda_D = 0.0 * np.pi / 180  # radians
    
    # Convert lat/lon to radians
    lat_rad = np.asarray(lat) * np.pi / 180
    lon_rad = np.asarray(lon) * np.pi / 180
    
    # Forward projection: geographic coordinates to viewing angles
    # Step 1: Convert geodetic latitude to geocentric coordinates
    cos_lat = np.cos(lat_rad)
    sin_lat = np.sin(lat_rad)
    cos_lon = np.cos(lon_rad)
    sin_lon = np.sin(lon_rad)
    
    # Geocentric latitude
    c = 1.0 / np.sqrt(cos_lat**2 + (r_pol/r_eq)**2 * sin_lat**2)
    
    # Position on Earth's surface in geocentric coordinates
    x = r_eq * c * cos_lat * cos_lon
    y = r_eq * c * cos_lat * sin_lon
    z = r_pol * c * (r_pol/r_eq) * sin_lat
    
    # Step 2: Calculate satellite position
    x_sat = h
    y_sat = 0.0
    z_sat = 0.0
    
    # Step 3: Vector from satellite to point on Earth
    dx = x - x_sat
    dy = y
    dz = z
    
    # Step 4: Check visibility
    # Point must be on the Earth-facing hemisphere
    distance = np.sqrt(dx**2 + dy**2 + dz**2)
    
    # For a point to be visible, the angle between the satellite-to-point vector
    # and the satellite-to-Earth-center vector must be less than the horizon angle
    cos_angle = -dx / distance  # dot product with normalized vectors
    horizon_angle = np.arcsin(r_eq / h)
    visible = cos_angle > np.cos(np.pi/2 + horizon_angle)
    
    # Step 5: Calculate viewing angles
    # phi_s: elevation angle (North-South)
    # lambda_s: azimuth angle (East-West, but note the sign convention!)
    
    # Distance in the equatorial plane
    rho_c = np.sqrt(dx**2 + dy**2)
    
    # Elevation angle
    phi_s = np.arctan(dz / rho_c)
    
    # Azimuth angle (note: increases from East to West in the viewing convention)
    lambda_s = np.arctan(dy / dx)
    
    # Convert viewing angles to row and column
    # Note: equations from Section 5.2, solved for row and column
    column = (params['lambda_0'] - lambda_s) / params['azimuth_sampling'] + 1
    row = (phi_s - params['phi_0']) / params['elevation_sampling'] + 1
    
    # Set invalid points to NaN
    row = np.where(visible, row, np.nan)
    column = np.where(visible, column, np.nan)
    
    return row, column

def fci_grid_to_latlon(row, column, ssd=1.0):
    """
    Convert FCI Level 1c reference grid row/column to latitude/longitude.
    
    Parameters:
    -----------
    row : int or array-like
        Row number(s) in the reference grid (1-indexed, starting from south)
    column : int or array-like
        Column number(s) in the reference grid (1-indexed, starting from west)
    ssd : float
        Spatial Sampling Distance in km (0.5, 1.0, or 2.0)
    
    Returns:
    --------
    lon, lat : float or array
        Longitude and latitude in degrees (NaN for space pixels)
    earth_mask : bool or array
        True for pixels viewing Earth, False for deep space pixels
    
    Notes:
    ------
    Based on MTG FCI L1 Product User Guide sections 5.2 and 5.3
    """
    
    # Grid parameters for different SSDs (from Table 3 in the document)
    grid_params = {
        0.5: {
            'lambda_0': 8.9142405037 * np.pi / 180,  # radians
            'phi_0': -8.9142405037 * np.pi / 180,    # radians
            'azimuth_sampling': 0.000800524494 * np.pi / 180,  # radians
            'elevation_sampling': 0.000800524494 * np.pi / 180  # radians
        },
        1.0: {
            'lambda_0': 8.9138402398 * np.pi / 180,
            'phi_0': -8.9138402398 * np.pi / 180,
            'azimuth_sampling': 0.001601048988 * np.pi / 180,
            'elevation_sampling': 0.001601048988 * np.pi / 180
        },
        2.0: {
            'lambda_0': 8.9130397083 * np.pi / 180,
            'phi_0': -8.9130397083 * np.pi / 180,
            'azimuth_sampling': 0.003202097973 * np.pi / 180,
            'elevation_sampling': 0.003202097973 * np.pi / 180
        }
    }
    
    if ssd not in grid_params:
        raise ValueError("SSD must be 0.5, 1.0, or 2.0 km")
    
    params = grid_params[ssd]
    
    # Earth ellipsoid parameters
    r_eq = 6378.137  # km - equatorial radius
    r_pol = r_eq * (1 - 1/298.257223563)  # km - polar radius
    h = 35786.4 + r_eq  # km - geostationary radius (altitude + equatorial radius)
    
    # Sub-satellite point longitude (degrees)
    lambda_D = 0.0  # degrees
    
    # Calculate viewing angles (lambda_s, phi_s) from row and column
    # Note: rows and columns are 1-indexed
    lambda_s = params['lambda_0'] - (column - 1) * params['azimuth_sampling']
    phi_s = params['phi_0'] + (row - 1) * params['elevation_sampling']
    
    # Normalized Geostationary Projection (inverse projection)
    # Calculate intermediate values
    cos_lambda_s = np.cos(lambda_s)
    sin_lambda_s = np.sin(lambda_s)
    cos_phi_s = np.cos(phi_s)
    sin_phi_s = np.sin(phi_s)
    s5 = h**2 - r_eq**2
    s4 = (r_eq / r_pol)**2

    s_d_2 = (h * cos_lambda_s * cos_phi_s)**2 - s5 * (cos_phi_s**2 + s4 * sin_phi_s**2)
    s_d = np.sqrt(s_d_2)

    # Calculate s_n
    s_n = (h * cos_lambda_s * cos_phi_s - s_d ) / (cos_phi_s**2 + s4 * sin_phi_s**2)
    # Calculate s1, s2, s3
    s1 = h - s_n * cos_lambda_s * cos_phi_s
    s2 = - s_n * sin_lambda_s * cos_phi_s
    s3 = s_n * sin_phi_s
    
    # Calculate s_xy and s4
    s_xy = np.sqrt(s1**2 + s2**2)
    
    # Check if the line of sight intersects Earth
    # The discriminant s_d must be non-negative for Earth intersection
    # Create mask for valid Earth pixels
    earth_mask = s_d >= 0
    
    # Initialize output arrays with NaN
    lat_deg = np.full_like(lambda_s, np.nan, dtype=float)
    lon_deg = np.full_like(lambda_s, np.nan, dtype=float)
    
    # Only calculate lat/lon for pixels that see Earth
    if np.any(earth_mask):
        # Calculate latitude and longitude for Earth pixels
        lat = np.arctan((s3 / s_xy) * s4)
        lon = np.arctan(s2 / s1) + lambda_D * np.pi / 180
        
        # Convert to degrees
        lat_deg = np.where(earth_mask, lat * 180 / np.pi, np.nan)
        lon_deg = np.where(earth_mask, lon * 180 / np.pi, np.nan)
    
    return lon_deg, lat_deg, earth_mask


# Example usage
if __name__ == "__main__":
    # Example: Check corner points (these should be in space!)
    print("Corner points (should be in deep space):")
    corners = [(1, 1), (1, 11136), (11136, 1), (11136, 11136)]
    for row, col in corners:
        lon, lat, is_earth = fci_grid_to_latlon(row, col, ssd=1.0)
        status = "Earth" if is_earth else "Space"
        print(f"Row {row:5d}, Column {col:5d} -> {status:5s}: Lon: {lon:8.4f}°, Lat: {lat:8.4f}°")
    
    # Example: Center point (sub-satellite point, should be on Earth)
    print("\nCenter point (sub-satellite, should be on Earth):")
    row, col = 5568, 5568
    lon, lat, is_earth = fci_grid_to_latlon(row, col, ssd=1.0)
    status = "Earth" if is_earth else "Space"
    print(f"Row {row}, Column {col} -> {status}: Lon: {lon:.4f}°, Lat: {lat:.4f}°")
    
    # Example: Calculate for multiple pixels
    print("\nMultiple pixels (mix of Earth and space):")
    rows = np.array([1, 3000, 5568, 8000, 11136])
    cols = np.array([1, 3000, 5568, 8000, 11136])
    lons, lats, earth_masks = fci_grid_to_latlon(rows, cols, ssd=1.0)
    
    for r, c, lon, lat, is_earth in zip(rows, cols, lons, lats, earth_masks):
        status = "Earth" if is_earth else "Space"
        print(f"Row {r:5d}, Column {c:5d} -> {status:5s}: Lon: {lon:8.4f}°, Lat: {lat:8.4f}°")
    
    # Example: Create a grid and count Earth vs Space pixels
    print("\nFull disc analysis:")
    row_grid, col_grid = np.meshgrid(
        np.arange(1, 11137, 100),  # Sample every 100 rows
        np.arange(1, 11137, 100),  # Sample every 100 columns
        indexing='ij'
    )
    lon_grid, lat_grid, earth_grid = fci_grid_to_latlon(row_grid, col_grid, ssd=1.0)
    
    n_earth = np.sum(earth_grid)
    n_space = np.sum(~earth_grid)
    total = earth_grid.size
    
    print(f"Sampled grid shape: {lon_grid.shape}")
    print(f"Earth pixels: {n_earth} ({100*n_earth/total:.1f}%)")
    print(f"Space pixels: {n_space} ({100*n_space/total:.1f}%)")
    print(f"Lon range (Earth only): {np.nanmin(lon_grid):.4f}° to {np.nanmax(lon_grid):.4f}°")
    print(f"Lat range (Earth only): {np.nanmin(lat_grid):.4f}° to {np.nanmax(lat_grid):.4f}°")
    print("")

    # Your database of locations
    lats = np.array([50.85, 51.51, 48.86])  # Brussels, London, Paris
    lons = np.array([4.35, -0.13, 2.35])

    # Convert to FCI grid coordinates
    rows, cols = fci_latlon_to_grid(lats, lons, ssd=1.0)
    # Now you can extract data from FCI images at these locations
    # Note: round to nearest integer for pixel indexing
    row_indices = np.round(rows).astype(int)
    col_indices = np.round(cols).astype(int)


    lon_new_grid, lat_new_grid, earth_new_grid = fci_grid_to_latlon(row_indices, col_indices, ssd=1.0)

    for lat, lon, row, col, is_earth in zip(lat_new_grid, lon_new_grid, row_indices, col_indices, earth_new_grid):
        status = "Earth" if is_earth else "Space"
        print(f"Location (Lat: {lat:.2f}°, Lon: {lon:.2f}°) -> Row: {row}, Column: {col} -> {status}")
    #rows, cols = fci_latlon_to_grid(0.0, 0.0, ssd=1.0)
    #print(f"\nEquator and Prime Meridian (0°, 0°) -> Row: {rows}, Column: {cols}")  
    p_row, p_col = fci_latlon_to_grid(50.85, 4.35, ssd=1.0)
    print(f"\nBrussels (50.85°N, 4.35°E) -> Row: {p_row}, Column: {p_col}")
    lon_bxl, lat_bxl, is_earth_bxl = fci_grid_to_latlon(np.round(p_row).astype(int), np.round(p_col).astype(int), ssd=1.0)
    status_bxl = "Earth" if is_earth_bxl else "Space"
    print(f"Back conversion -> Lon: {lon_bxl:.4f}°, Lat: {lat_bxl:.4f}° -> {status_bxl}")
    