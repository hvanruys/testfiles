"""
Rayleigh correction for FCI MTG-I1 imagery
Python implementation with NumPy
"""

import numpy as np
from typing import Tuple, Optional
from dataclasses import dataclass

# Constants for Rayleigh scattering
EARTH_SUN_DISTANCE = 1.0  # AU, adjust seasonally if needed
RAYLEIGH_PHASE_FUNCTION = 0.9596  # (1 + cos²θ) normalized


@dataclass
class CorrectionParams:
    """Parameters for Rayleigh correction"""
    rayleigh_optical_depth: float
    surface_pressure: float = 1013.25  # hPa
    band_number: int = 1


class RayleighCorrector:
    """Rayleigh correction for FCI MTG-I1 imagery"""
    
    @staticmethod
    def get_rayleigh_optical_depth(band_number: int) -> float:
        """
        Get Rayleigh optical depth for FCI bands
        
        Args:
            band_number: FCI band number (1-6 for VIS/NIR)
            
        Returns:
            Rayleigh optical depth
        """
        optical_depths = {
            1: 0.2174,  # VIS 0.4 (444 nm)
            2: 0.1121,  # VIS 0.5 (510 nm)
            3: 0.0492,  # VIS 0.6 (640 nm)
            4: 0.0318,  # VIS 0.8 (865 nm)
            5: 0.0087,  # NIR 1.6 (1640 nm)
            6: 0.0028,  # NIR 2.2 (2250 nm)
        }
        return optical_depths.get(band_number, 0.0)
    
    @staticmethod
    def correct(
        reflectance: np.ndarray,
        solar_zenith: np.ndarray,
        view_zenith: np.ndarray,
        relative_azimuth: np.ndarray,
        params: CorrectionParams
    ) -> np.ndarray:
        """
        Perform Rayleigh correction on reflectance data
        
        Args:
            reflectance: Input TOA reflectance values
            solar_zenith: Solar zenith angles in degrees
            view_zenith: Viewing zenith angles in degrees
            relative_azimuth: Relative azimuth angles in degrees
            params: Correction parameters
            
        Returns:
            Rayleigh-corrected surface reflectance
        """
        # Validate input shapes
        if not (reflectance.shape == solar_zenith.shape == 
                view_zenith.shape == relative_azimuth.shape):
            raise ValueError("Input array shapes must match")
        
        # Adjust optical depth for actual surface pressure
        tau_r = params.rayleigh_optical_depth * (params.surface_pressure / 1013.25)
        
        # Convert angles to radians
        theta_s = np.radians(solar_zenith)
        theta_v = np.radians(view_zenith)
        phi = np.radians(relative_azimuth)
        
        # Cosines
        mu_s = np.cos(theta_s)
        mu_v = np.cos(theta_v)
        
        # Create output array
        corrected = np.full_like(reflectance, -999.0)
        
        # Valid data mask
        valid = (mu_s > 0) & (mu_v > 0) & (reflectance >= 0)
        
        # Process only valid pixels
        if np.any(valid):
            # Scattering angle
            cos_scatter = (-mu_s[valid] * mu_v[valid] + 
                          np.sin(theta_s[valid]) * np.sin(theta_v[valid]) * 
                          np.cos(phi[valid]))
            
            # Rayleigh phase function: P(Θ) = 3/4 * (1 + cos²Θ)
            phase = 0.75 * (1.0 + cos_scatter**2)
            
            # Atmospheric transmission
            T_s = np.exp(-tau_r / mu_s[valid])  # Sun to surface
            T_v = np.exp(-tau_r / mu_v[valid])  # Surface to sensor
            T_total = T_s * T_v
            
            # Rayleigh reflectance (path radiance contribution)
            rho_ray = ((tau_r * phase) / (4.0 * np.pi * mu_s[valid] * mu_v[valid]) * 
                      (1.0 - np.exp(-tau_r * (1.0/mu_s[valid] + 1.0/mu_v[valid]))))
            
            # Correct for multiple scattering (simple approximation)
            rho_multiple = rho_ray / (1.0 - 0.5 * tau_r)
            
            # Remove Rayleigh contribution and account for transmission
            rho_surface = (reflectance[valid] - rho_multiple) / T_total
            
            # Ensure non-negative values
            corrected[valid] = np.maximum(0.0, rho_surface)
        
        return corrected
    
    @staticmethod
    def correct_image(
        image: np.ndarray,
        solar_zenith: np.ndarray,
        view_zenith: np.ndarray,
        relative_azimuth: np.ndarray,
        params: CorrectionParams
    ) -> np.ndarray:
        """
        Correct a 2D image array
        
        Args:
            image: Input image (height, width) with TOA reflectance
            solar_zenith: Solar zenith angles (height, width) in degrees
            view_zenith: Viewing zenith angles (height, width) in degrees
            relative_azimuth: Relative azimuth angles (height, width) in degrees
            params: Correction parameters
            
        Returns:
            Corrected image array
        """
        # Flatten arrays for processing
        reflectance_flat = image.ravel()
        solar_zenith_flat = solar_zenith.ravel()
        view_zenith_flat = view_zenith.ravel()
        relative_azimuth_flat = relative_azimuth.ravel()
        
        # Perform correction
        corrected_flat = RayleighCorrector.correct(
            reflectance_flat,
            solar_zenith_flat,
            view_zenith_flat,
            relative_azimuth_flat,
            params
        )
        
        # Reshape back to image dimensions
        return corrected_flat.reshape(image.shape)


def example_usage():
    """Example usage of the RayleighCorrector"""
    
    # Create sample data (1000x1000 image)
    height, width = 1000, 1000
    
    # Generate synthetic TOA reflectance
    reflectance = np.random.rand(height, width) * 0.5
    
    # Generate angle arrays (these should come from your MTG-I1 data)
    solar_zenith = np.random.rand(height, width) * 60  # 0-60 degrees
    view_zenith = np.random.rand(height, width) * 45   # 0-45 degrees
    relative_azimuth = np.random.rand(height, width) * 180  # 0-180 degrees
    
    # Set correction parameters for band 3 (VIS 0.6)
    params = CorrectionParams(
        band_number=3,
        rayleigh_optical_depth=RayleighCorrector.get_rayleigh_optical_depth(3),
        surface_pressure=1013.25
    )
    
    # Perform correction
    print(f"Processing {height}x{width} image for band {params.band_number}")
    corrected = RayleighCorrector.correct_image(
        reflectance,
        solar_zenith,
        view_zenith,
        relative_azimuth,
        params
    )
    
    print(f"Correction complete!")
    print(f"Input reflectance range: [{reflectance.min():.4f}, {reflectance.max():.4f}]")
    print(f"Corrected reflectance range: [{corrected[corrected >= 0].min():.4f}, "
          f"{corrected[corrected >= 0].max():.4f}]")
    print(f"Invalid pixels: {np.sum(corrected < 0)}")
    
    return corrected


# Example with real-world workflow
def process_fci_band(
    reflectance_file: str,
    angle_data: dict,
    band_number: int,
    output_file: Optional[str] = None
) -> np.ndarray:
    """
    Process a single FCI band with Rayleigh correction
    
    Args:
        reflectance_file: Path to reflectance data file
        angle_data: Dictionary containing 'solar_zenith', 'view_zenith', 
                   'relative_azimuth' arrays
        band_number: FCI band number
        output_file: Optional path to save corrected data
        
    Returns:
        Corrected reflectance array
    """
    # Load reflectance data (adjust based on your file format)
    # reflectance = np.load(reflectance_file)  # For .npy files
    # or use netCDF4, h5py, etc. for HDF5/NetCDF files
    
    # Set correction parameters
    params = CorrectionParams(
        band_number=band_number,
        rayleigh_optical_depth=RayleighCorrector.get_rayleigh_optical_depth(band_number),
        surface_pressure=1013.25  # Adjust if you have surface pressure data
    )
    
    # Perform correction
    # corrected = RayleighCorrector.correct_image(
    #     reflectance,
    #     angle_data['solar_zenith'],
    #     angle_data['view_zenith'],
    #     angle_data['relative_azimuth'],
    #     params
    # )
    
    # Save if output file specified
    # if output_file:
    #     np.save(output_file, corrected)
    
    # return corrected
    pass


if __name__ == "__main__":
    # Run example
    corrected_image = example_usage()