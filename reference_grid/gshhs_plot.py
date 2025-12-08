import struct
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import cartopy.crs as ccrs

class GSHHSReader:
    """Reader for GSHHS (Global Self-consistent Hierarchical High-resolution Geography) binary files."""
    
    def __init__(self, filepath):
        self.filepath = filepath
        self.polygons = []
        
    def read(self):
        """Read all polygons from the GSHHS binary file."""
        self.polygons = []
        
        with open(self.filepath, 'rb') as f:
            while True:
                # Read polygon header (44 bytes)
                header_data = f.read(44)
                if len(header_data) < 44:
                    break
                
                # Parse header (big-endian format)
                header = struct.unpack('>iiiiiiiiiii', header_data)
                
                polygon_id = header[0]
                n_points = header[1]
                flag = header[2]
                west = header[3] * 1e-6   # Convert to degrees
                east = header[4] * 1e-6
                south = header[5] * 1e-6
                north = header[6] * 1e-6
                area = header[7]
                area_full = header[8]
                container = header[9]
                ancestor = header[10]
                
                # Level: 1=land, 2=lake, 3=island_in_lake, 4=pond_in_island
                level = (flag >> 8) & 255
                
                # Read polygon points
                points_data = f.read(n_points * 8)
                if len(points_data) < n_points * 8:
                    break
                
                # Parse points (longitude, latitude pairs)
                lons = []
                lats = []
                for i in range(n_points):
                    lon, lat = struct.unpack('>ii', points_data[i*8:(i+1)*8])
                    lons.append(lon * 1e-6)  # Convert to degrees
                    lats.append(lat * 1e-6)
                
                polygon = {
                    'id': polygon_id,
                    'n_points': n_points,
                    'level': level,
                    'bounds': (west, east, south, north),
                    'area': area,
                    'lons': np.array(lons),
                    'lats': np.array(lats)
                }
                
                self.polygons.append(polygon)
        
        return self.polygons
    
    def filter_by_level(self, level):
        """Filter polygons by level (1=land, 2=lake, 3=island, 4=pond)."""
        return [p for p in self.polygons if p['level'] == level]
    
    def filter_by_bounds(self, west, east, south, north):
        """Filter polygons that intersect with given bounds."""
        filtered = []
        for p in self.polygons:
            pw, pe, ps, pn = p['bounds']
            if not (pe < west or pw > east or pn < south or ps > north):
                filtered.append(p)
        return filtered


def plot_gshhs(polygons, title='GSHHS Coastlines', figsize=(15, 10), 
               color_map=None, show_levels=True, central_lat=0, central_lon=0, 
               distance=35786400, max_jump=10):
    """
    Plot GSHHS polygons using general perspective projection.
    
    Parameters:
    -----------
    polygons : list
        List of polygon dictionaries from GSHHSReader
    title : str
        Plot title
    figsize : tuple
        Figure size (width, height)
    color_map : dict
        Dictionary mapping level to color. If None, uses default colors.
    show_levels : bool
        Whether to show different levels in different colors
    central_lat : float
        Central latitude for perspective projection (degrees)
    central_lon : float
        Central longitude for perspective projection (degrees)
    distance : float
        Distance from center of Earth to observation point (in meters)
        Default: 35786400 m (geostationary altitude)
    max_jump : float
        Maximum allowed longitude jump in degrees before breaking line (default: 10)
    """
    if color_map is None:
        color_map = {
            1: 'brown',      # Land
            2: 'lightblue',  # Lake
            3: 'brown',      # Island in lake
            4: 'lightblue'   # Pond in island
        }
    
    # Create figure with general perspective projection
    proj = ccrs.Geostationary(central_longitude=central_lon, 
                              satellite_height=distance)
    fig, ax = plt.subplots(figsize=figsize, 
                            subplot_kw=dict(projection=proj))
    
    for polygon in polygons:
        level = polygon['level']
        color = color_map.get(level, 'gray') if show_levels else 'black'
        
        lons = polygon['lons']
        lats = polygon['lats']
        
        # Break polygon into segments to avoid long lines across dateline
        lon_diffs = np.abs(np.diff(lons))
        breaks = np.where(lon_diffs > max_jump)[0]
        
        if len(breaks) == 0:
            # No discontinuities, plot entire polygon
            ax.plot(lons, lats, color=color, linewidth=0.5,
                   transform=ccrs.PlateCarree())
        else:
            # Plot segments between breaks
            start = 0
            for break_idx in breaks:
                end = break_idx + 1
                ax.plot(lons[start:end], lats[start:end], 
                       color=color, linewidth=0.5,
                       transform=ccrs.PlateCarree())
                start = end
            # Plot final segment
            ax.plot(lons[start:], lats[start:], 
                   color=color, linewidth=0.5,
                   transform=ccrs.PlateCarree())
    
    ax.coastlines(resolution='50m', linewidth=0.5, alpha=0.5)
    ax.gridlines(draw_labels=False, alpha=0.3)
    ax.set_title(title)
    
    if show_levels:
        # Create legend
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='brown', label='Land/Islands'),
            Patch(facecolor='lightblue', label='Lakes/Ponds')
        ]
        ax.legend(handles=legend_elements, loc='upper right')
    
    plt.tight_layout()
    return fig, ax


# Example usage
if __name__ == '__main__':
    # Example 1: Read and plot entire file
    print("Reading GSHHS file...")
    
    # Replace with your GSHHS file path
    # Download from: https://www.ngdc.noaa.gov/mgg/shorelines/gshhs.html
    filepath = 'gshhs_c.b'  # 'c' = crude resolution
    
    try:
        reader = GSHHSReader(filepath)
        polygons = reader.read()
        
        print(f"Read {len(polygons)} polygons")
        
        # Plot all polygons
        fig, ax = plot_gshhs(polygons, title='GSHHS World Coastlines (Crude)')
        #plt.savefig('gshhs_world.png', dpi=150)
        plt.show()
        
        # Example 2: Plot only land masses (level 1)
        #land_polygons = reader.filter_by_level(1)
        #print(f"\nLand polygons: {len(land_polygons)}")
        
        #fig, ax = plot_gshhs(land_polygons, 
        #                    title='GSHHS Land Masses Only',
        #                    show_levels=False)
        #plt.savefig('gshhs_land.png', dpi=150)
        #plt.show()
        
        # Example 3: Plot specific region (e.g., Mediterranean)
        #med_polygons = reader.filter_by_bounds(
        #    west=-10, east=40, south=30, north=50
        #)
        #print(f"\nMediterranean region polygons: {len(med_polygons)}")
        
        #fig, ax = plot_gshhs(med_polygons, 
        #                    title='Mediterranean Region')
        #ax.set_xlim(-10, 40)
        #ax.set_ylim(30, 50)
        #plt.savefig('gshhs_mediterranean.png', dpi=150)
        #plt.show()
        
        # Print some statistics
        print("\nPolygon statistics:")
        print(f"Total points: {sum(p['n_points'] for p in polygons)}")
        print(f"Level distribution:")
        for level in range(1, 5):
            count = len(reader.filter_by_level(level))
            level_names = {1: 'Land', 2: 'Lake', 3: 'Island', 4: 'Pond'}
            print(f"  Level {level} ({level_names[level]}): {count}")
        
    except FileNotFoundError:
        print(f"Error: File '{filepath}' not found.")
        print("\nTo use this script:")
        print("1. Download GSHHS files from:")
        print("   https://www.ngdc.noaa.gov/mgg/shorelines/gshhs.html")
        print("2. Extract the .b files (e.g., gshhs_c.b, gshhs_l.b, etc.)")
        print("3. Update the 'filepath' variable with your file location")
        print("\nResolutions available:")
        print("  gshhs_c.b - Crude (lowest resolution)")
        print("  gshhs_l.b - Low")
        print("  gshhs_i.b - Intermediate")
        print("  gshhs_h.b - High")
        print("  gshhs_f.b - Full (highest resolution)")