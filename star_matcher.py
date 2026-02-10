"""
Star Matching Algorithm for Astrometry
Digantara Assessment - Image Processing Engineer

This implements a robust star matching and plate solving algorithm using:
- Triangle pattern matching with scale/rotation invariant features
- KD-tree for efficient pattern retrieval
- RANSAC-based verification
- Non-linear least squares optimization for orientation estimation
"""

import numpy as np
import pandas as pd
from scipy.spatial import KDTree
from scipy.optimize import least_squares
from itertools import combinations
import matplotlib.pyplot as plt
from typing import Tuple, Dict, List, Optional
import json
import warnings
warnings.filterwarnings('ignore')


class StarMatcher:
    """
    Robust star matching algorithm for astrometry and plate solving.
    """
    
    def __init__(self, focal_length_mm: float, pixel_pitch_um: float, 
                 image_size: Tuple[int, int], ra_deg: float, dec_deg: float,
                 fov_x_deg: float, fov_y_deg: float):
        """
        Initialize star matcher with sensor parameters.
        
        Parameters:
        -----------
        focal_length_mm : float
            Focal length in millimeters
        pixel_pitch_um : float
            Pixel pitch in micrometers
        image_size : tuple
            Image dimensions (width, height) in pixels
        ra_deg, dec_deg : float
            Approximate field center (RA, DEC) in degrees
        fov_x_deg, fov_y_deg : float
            Field of view in degrees
        """
        self.focal_length_mm = focal_length_mm
        self.pixel_pitch_um = pixel_pitch_um
        self.image_width, self.image_height = image_size
        self.ra_deg = ra_deg
        self.dec_deg = dec_deg
        self.fov_x_deg = fov_x_deg
        self.fov_y_deg = fov_y_deg
        
        # Calculate plate scale (arcsec/pixel)
        # Formula: scale = (pixel_size / focal_length) * 206265
        self.plate_scale = (pixel_pitch_um / focal_length_mm) * 206265 / 1000
        
        print(f"Sensor Configuration:")
        print(f"  Image size: {self.image_width} x {self.image_height} pixels")
        print(f"  Focal length: {self.focal_length_mm} mm")
        print(f"  Pixel pitch: {self.pixel_pitch_um} µm")
        print(f"  Plate scale: {self.plate_scale:.4f} arcsec/pixel")
        print(f"  Field of View: {self.fov_x_deg:.4f}° x {self.fov_y_deg:.4f}°")
        
        # Storage for catalog and results
        self.catalog = None
        self.matched_image_coords = None
        self.matched_catalog_coords = None
        
    def generate_synthetic_catalog(self, n_stars: int = 500) -> pd.DataFrame:
        """
        Generate synthetic star catalog for the field.
        
        In production, this would query GAIA DR3 or Hipparcos using:
        - astroquery.vizier for catalog access
        - Cone search around field center
        - Magnitude filtering
        
        For this assessment, we simulate realistic catalog data.
        """
        print(f"\n[A1] Generating synthetic catalog...")
        print(f"     Field center: RA={self.ra_deg:.4f}°, DEC={self.dec_deg:.4f}°")
        
        # Generate stars with margin for field rotation/offset
        margin_factor = 1.8
        ra_min = self.ra_deg - self.fov_x_deg * margin_factor / 2
        ra_max = self.ra_deg + self.fov_x_deg * margin_factor / 2
        dec_min = self.dec_deg - self.fov_y_deg * margin_factor / 2
        dec_max = self.dec_deg + self.fov_y_deg * margin_factor / 2
        
        # Handle RA wrap-around at 0°/360°
        if ra_min < 0:
            ra_min += 360
        if ra_max > 360:
            ra_max -= 360
            
        # Generate star positions
        np.random.seed(42)  # Reproducible results
        ra_stars = np.random.uniform(ra_min, ra_max, n_stars)
        dec_stars = np.random.uniform(dec_min, dec_max, n_stars)
        
        # Realistic magnitude distribution (exponential, brighter stars rarer)
        magnitudes = np.random.exponential(scale=2.5, size=n_stars) + 4.0
        magnitudes = np.clip(magnitudes, 4.0, 13.0)
        
        catalog = pd.DataFrame({
            'RA': ra_stars,
            'DEC': dec_stars,
            'Magnitude': magnitudes
        })
        
        # Sort by magnitude (brightest first)
        catalog = catalog.sort_values('Magnitude').reset_index(drop=True)
        
        print(f"     Generated {len(catalog)} catalog stars")
        print(f"     Magnitude range: {catalog['Magnitude'].min():.1f} - {catalog['Magnitude'].max():.1f}")
        
        self.catalog = catalog
        return catalog
    
    def create_triangle_patterns(self, coords: np.ndarray, 
                                 n_patterns: int = 200,
                                 name: str = "patterns") -> List[Dict]:
        """
        Create triangle patterns for star matching.
        
        Triangle patterns use scale and rotation invariant features:
        - Side length ratios (invariant to scaling)
        - Angular features (invariant to rotation)
        
        This is more robust than direct coordinate matching.
        
        Parameters:
        -----------
        coords : np.ndarray
            Star coordinates [N x 2]
        n_patterns : int
            Number of patterns to create
        name : str
            Identifier for logging
            
        Returns:
        --------
        patterns : list of dict
            Each pattern contains triangle geometry and invariant features
        """
        n_stars = len(coords)
        patterns = []
        
        print(f"     Creating {name}...")
        print(f"     Using {n_stars} stars")
        
        # Use brightest/central stars for pattern creation
        max_pattern_stars = min(30, n_stars)
        
        for i in range(max_pattern_stars):
            # Find nearest neighbors (Euclidean distance)
            distances = np.linalg.norm(coords - coords[i], axis=1)
            nearest = np.argsort(distances)[1:12]  # Exclude self, get 11 nearest
            
            # Create triangles from star i and pairs of nearest neighbors
            for j, k in combinations(nearest[:9], 2):
                p1, p2, p3 = coords[i], coords[j], coords[k]
                
                # Calculate side lengths
                d12 = np.linalg.norm(p1 - p2)
                d23 = np.linalg.norm(p2 - p3)
                d31 = np.linalg.norm(p3 - p1)
                
                # Skip degenerate triangles
                if min(d12, d23, d31) < 1e-6:
                    continue
                
                # Sort sides for canonical representation
                sides = sorted([d12, d23, d31])
                
                # Scale-invariant features: side ratios
                ratio1 = sides[1] / sides[2]
                ratio2 = sides[0] / sides[2]
                
                # Additional feature: largest angle cosine (rotation invariant)
                # Cosine rule: cos(C) = (a² + b² - c²) / (2ab)
                cos_largest = ((sides[0]**2 + sides[1]**2 - sides[2]**2) / 
                              (2 * sides[0] * sides[1] + 1e-10))
                cos_largest = np.clip(cos_largest, -1, 1)
                
                patterns.append({
                    'indices': (i, j, k),
                    'coords': np.array([p1, p2, p3]),
                    'sides': sides,
                    'ratio1': ratio1,
                    'ratio2': ratio2,
                    'cos_angle': cos_largest,
                    'perimeter': sum(sides)
                })
                
                if len(patterns) >= n_patterns:
                    break
            
            if len(patterns) >= n_patterns:
                break
        
        print(f"     Created {len(patterns)} triangle patterns")
        return patterns
    
    def match_patterns(self, image_patterns: List[Dict], 
                      catalog_patterns: List[Dict],
                      tolerance: float = 0.06) -> List[Dict]:
        """
        Match triangle patterns using KD-tree for efficient retrieval.
        
        Parameters:
        -----------
        image_patterns : list
            Patterns from image centroids
        catalog_patterns : list
            Patterns from catalog stars
        tolerance : float
            Feature distance threshold
            
        Returns:
        --------
        matches : list
            Pattern correspondences with quality scores
        """
        if not catalog_patterns:
            return []
        
        print(f"     Matching {len(image_patterns)} image patterns against {len(catalog_patterns)} catalog patterns")
        
        # Build feature vectors (ratio1, ratio2)
        catalog_features = np.array([[p['ratio1'], p['ratio2']] 
                                     for p in catalog_patterns])
        
        # Build KD-tree for O(log n) queries
        tree = KDTree(catalog_features)
        
        matches = []
        for img_pattern in image_patterns:
            img_feature = np.array([img_pattern['ratio1'], img_pattern['ratio2']])
            
            # Query k nearest neighbors
            distances, indices = tree.query(img_feature, k=min(5, len(catalog_features)))
            
            # Accept matches within tolerance
            for dist, idx in zip(distances, indices):
                if dist < tolerance:
                    matches.append({
                        'image_pattern': img_pattern,
                        'catalog_pattern': catalog_patterns[idx],
                        'distance': dist,
                        'quality': np.exp(-dist * 10)  # Quality score
                    })
        
        print(f"     Found {len(matches)} pattern matches")
        return matches
    
    def verify_matches(self, matches: List[Dict], 
                      centroid_data: pd.DataFrame,
                      min_votes: int = 2) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
        """
        [A4] Verify star matches using voting mechanism.
        
        Each triangle match provides 3 point correspondences.
        Stars appearing in multiple consistent triangles are verified.
        This is similar to RANSAC consensus.
        
        Parameters:
        -----------
        matches : list
            Pattern matches from matching step
        centroid_data : DataFrame
            Image centroid data
        min_votes : int
            Minimum votes required for verification
            
        Returns:
        --------
        image_coords, catalog_coords : np.ndarray or None
            Verified correspondences
        """
        print(f"\n[A4] Verifying matches with voting...")
        
        # Collect votes: image_idx -> list of (catalog_idx, quality)
        votes = {}
        
        for match in matches:
            img_indices = match['image_pattern']['indices']
            cat_indices = match['catalog_pattern']['indices']
            quality = match['quality']
            
            # Each triangle vertex is a potential correspondence
            for img_idx, cat_idx in zip(img_indices, cat_indices):
                if img_idx not in votes:
                    votes[img_idx] = []
                votes[img_idx].append((cat_idx, quality))
        
        # Verify by consensus voting
        verified = {}
        vote_counts = {}
        
        for img_idx, candidates in votes.items():
            if len(candidates) >= min_votes:
                # Find most common catalog match (weighted by quality)
                cat_indices = [c[0] for c in candidates]
                qualities = [c[1] for c in candidates]
                
                unique_cats, inverse = np.unique(cat_indices, return_inverse=True)
                
                # Sum quality scores for each unique catalog star
                quality_sums = np.zeros(len(unique_cats))
                for i, cat in enumerate(unique_cats):
                    mask = (inverse == i)
                    quality_sums[i] = sum(q for j, q in enumerate(qualities) if mask[j])
                
                # Select best match
                best_idx = np.argmax(quality_sums)
                verified[img_idx] = unique_cats[best_idx]
                vote_counts[img_idx] = sum(inverse == best_idx)
        
        if not verified:
            print(f"     ERROR: No verified matches!")
            return None, None
        
        print(f"     Verified {len(verified)} star matches from {len(votes)} candidates")
        
        if len(verified) < 4:
            print(f"     WARNING: Only {len(verified)} matches (4+ recommended for robust astrometry)")
        
        # Extract coordinates
        image_indices = list(verified.keys())
        catalog_indices = list(verified.values())
        
        centroid_array = centroid_data[['x_centroid', 'y_centroid']].values
        image_coords = centroid_array[image_indices]
        catalog_coords = self.catalog[['RA', 'DEC']].values[catalog_indices]
        
        # Print statistics
        votes_list = [vote_counts[idx] for idx in image_indices]
        print(f"     Vote statistics: min={min(votes_list)}, max={max(votes_list)}, mean={np.mean(votes_list):.1f}")
        
        return image_coords, catalog_coords
    
    def project_catalog_to_image(self, catalog_coords: np.ndarray,
                                 ra_center: float, dec_center: float, 
                                 roll: float) -> np.ndarray:
        """
        Project sky coordinates to image plane using gnomonic projection.
        
        The gnomonic (tangent plane) projection is standard in astrometry:
        - Projects celestial sphere onto tangent plane
        - Great circles become straight lines
        - Preserves angles at projection center
        
        Transform sequence:
        1. Spherical → Standard coordinates (ξ, η) via gnomonic projection
        2. Apply field rotation
        3. Scale to pixels using plate scale
        4. Add image center offset
        
        Parameters:
        -----------
        catalog_coords : np.ndarray
            Sky coordinates [N x 2]: [[RA, DEC], ...]
        ra_center, dec_center : float
            Field center in degrees
        roll : float
            Field rotation in degrees
            
        Returns:
        --------
        image_coords : np.ndarray
            Pixel coordinates [N x 2]: [[x, y], ...]
        """
        # Convert to radians
        ra = np.radians(catalog_coords[:, 0])
        dec = np.radians(catalog_coords[:, 1])
        ra0 = np.radians(ra_center)
        dec0 = np.radians(dec_center)
        roll_rad = np.radians(roll)
        
        # Handle RA wrap-around
        dra = ra - ra0
        dra = np.where(dra > np.pi, dra - 2*np.pi, dra)
        dra = np.where(dra < -np.pi, dra + 2*np.pi, dra)
        
        # Gnomonic projection
        cos_dec = np.cos(dec)
        sin_dec = np.sin(dec)
        cos_dec0 = np.cos(dec0)
        sin_dec0 = np.sin(dec0)
        cos_dra = np.cos(dra)
        
        denominator = sin_dec0 * sin_dec + cos_dec0 * cos_dec * cos_dra
        
        # Standard coordinates (radians)
        xi = cos_dec * np.sin(dra) / denominator
        eta = (cos_dec0 * sin_dec - sin_dec0 * cos_dec * cos_dra) / denominator
        
        # Convert to arcseconds
        xi_arcsec = np.degrees(xi) * 3600
        eta_arcsec = np.degrees(eta) * 3600
        
        # Apply rotation
        cos_r = np.cos(roll_rad)
        sin_r = np.sin(roll_rad)
        
        xi_rot = xi_arcsec * cos_r - eta_arcsec * sin_r
        eta_rot = xi_arcsec * sin_r + eta_arcsec * cos_r
        
        # Convert to pixels
        x_rel = xi_rot / self.plate_scale
        y_rel = eta_rot / self.plate_scale
        
        # Image center offset
        cx = self.image_width / 2
        cy = self.image_height / 2
        
        x = cx + x_rel
        y = cy - y_rel  # Y increases downward
        
        return np.column_stack([x, y])
    
    def estimate_orientation(self, image_coords: np.ndarray, 
                           catalog_coords: np.ndarray) -> Dict:
        """
        Estimate optimal orientation (RA, DEC, ROLL) using least squares.
        
        This solves the plate solving problem:
        Given matched stars in image and catalog, find the transformation
        that minimizes reprojection error.
        
        Method: Non-linear least squares (Levenberg-Marquardt)
        Parameters: [RA_center, DEC_center, ROLL]
        Objective: Minimize sum of squared residuals
        
        Parameters:
        -----------
        image_coords : np.ndarray
            Observed image coordinates [N x 2]
        catalog_coords : np.ndarray
            Catalog coordinates [N x 2]
            
        Returns:
        --------
        orientation : dict
            Optimized parameters and error metrics
        """
        print(f"\n[A5] Estimating optimal orientation...")
        print(f"      Using {len(image_coords)} matched stars")
        
        # Initial guess
        x0 = np.array([self.ra_deg, self.dec_deg, 0.0])
        
        def residual_fn(params):
            """Compute reprojection residuals"""
            ra, dec, roll = params
            projected = self.project_catalog_to_image(catalog_coords, ra, dec, roll)
            return (image_coords - projected).flatten()
        
        # Optimize
        print(f"      Running non-linear least squares optimization...")
        result = least_squares(
            residual_fn, x0,
            method='lm',  # Levenberg-Marquardt
            max_nfev=1000,
            ftol=1e-8,
            xtol=1e-8
        )
        
        # Extract results
        ra_opt, dec_opt, roll_opt = result.x
        
        # Normalize angles
        ra_opt = ra_opt % 360
        roll_opt = roll_opt % 360
        if roll_opt > 180:
            roll_opt -= 360
        
        # Calculate errors
        projected = self.project_catalog_to_image(catalog_coords, ra_opt, dec_opt, roll_opt)
        residuals = image_coords - projected
        residual_norms = np.linalg.norm(residuals, axis=1)
        
        # Calculate transformation matrix
        cd_matrix = self.calculate_cd_matrix(ra_opt, dec_opt, roll_opt)
        
        orientation = {
            'RA': ra_opt,
            'DEC': dec_opt,
            'ROLL': roll_opt,
            'cd_matrix': cd_matrix,
            'n_matches': len(image_coords),
            'mean_residual_px': np.mean(residual_norms),
            'median_residual_px': np.median(residual_norms),
            'std_residual_px': np.std(residual_norms),
            'rms_residual_px': np.sqrt(np.mean(residual_norms**2)),
            'max_residual_px': np.max(residual_norms),
            'optimization_success': result.success,
            'cost': result.cost
        }
        
        # Convert to angular units
        orientation['rms_residual_arcsec'] = orientation['rms_residual_px'] * self.plate_scale
        
        print(f"      Optimization {'succeeded' if result.success else 'failed'}")
        print(f"      Final cost: {result.cost:.2e}")
        
        return orientation
    
    def calculate_cd_matrix(self, ra_center: float, dec_center: float, 
                           roll: float) -> np.ndarray:
        """
        Calculate WCS CD transformation matrix.
        
        The CD matrix maps pixel offsets to sky coordinate offsets:
        [dRA * cos(DEC)]  =  [CD1_1  CD1_2] [dx]
        [dDEC           ]     [CD2_1  CD2_2] [dy]
        
        Returns:
        --------
        cd_matrix : np.ndarray
            2x2 CD matrix in degrees/pixel
        """
        scale = self.plate_scale / 3600.0  # degrees/pixel
        roll_rad = np.radians(roll)
        
        cd = np.array([
            [scale * np.cos(roll_rad), -scale * np.sin(roll_rad)],
            [scale * np.sin(roll_rad),  scale * np.cos(roll_rad)]
        ])
        
        return cd
    
    def save_matched_data(self, centroid_data: pd.DataFrame,
                         image_coords: np.ndarray,
                         catalog_coords: np.ndarray,
                         output_path: str):
        """
        [A5] Save matched star data to CSV.
        
        Output format:
        x_centroid, y_centroid, brightness, RA, DEC, Magnitude
        """
        matched_records = []
        
        centroid_array = centroid_data[['x_centroid', 'y_centroid']].values
        
        for i, img_coord in enumerate(image_coords):
            # Find corresponding centroid
            distances = np.linalg.norm(centroid_array - img_coord, axis=1)
            idx = np.argmin(distances)
            
            # Find catalog star
            cat_coord = catalog_coords[i]
            cat_distances = np.linalg.norm(
                self.catalog[['RA', 'DEC']].values - cat_coord, axis=1
            )
            cat_idx = np.argmin(cat_distances)
            
            matched_records.append({
                'x_centroid': centroid_data.iloc[idx]['x_centroid'],
                'y_centroid': centroid_data.iloc[idx]['y_centroid'],
                'brightness': centroid_data.iloc[idx]['brightness'],
                'RA': catalog_coords[i, 0],
                'DEC': catalog_coords[i, 1],
                'Magnitude': self.catalog.iloc[cat_idx]['Magnitude']
            })
        
        df = pd.DataFrame(matched_records)
        df.to_csv(output_path, index=False)
        print(f"\n      Saved matched data: {output_path}")
        print(f"      Format: x_centroid, y_centroid, brightness, RA, DEC, Magnitude")
    
    def plot_results(self, centroid_data: pd.DataFrame,
                    image_coords: np.ndarray,
                    catalog_coords: np.ndarray,
                    orientation: Dict,
                    output_path: str):
        """
        Create visualization of matching results.
        """
        fig, axes = plt.subplots(1, 2, figsize=(16, 7))
        
        # Plot 1: Matched stars overlay
        ax1 = axes[0]
        
        ax1.scatter(centroid_data['x_centroid'], centroid_data['y_centroid'],
                   c='lightblue', s=20, alpha=0.4, label='All centroids', zorder=1)
        
        ax1.scatter(image_coords[:, 0], image_coords[:, 1],
                   c='red', s=80, marker='x', linewidths=2, 
                   label=f'Matched ({len(image_coords)})', zorder=2)
        
        # Label first few matches
        for i in range(min(8, len(image_coords))):
            ax1.annotate(str(i+1), (image_coords[i, 0], image_coords[i, 1]),
                        xytext=(5, 5), textcoords='offset points',
                        fontsize=8, color='darkred', fontweight='bold')
        
        ax1.set_xlabel('X (pixels)', fontsize=12, fontweight='bold')
        ax1.set_ylabel('Y (pixels)', fontsize=12, fontweight='bold')
        ax1.set_title('Star Centroids with Matches', fontsize=14, fontweight='bold')
        ax1.legend(fontsize=10)
        ax1.invert_yaxis()
        ax1.grid(True, alpha=0.3, linestyle='--')
        ax1.set_aspect('equal', adjustable='box')
        
        # Plot 2: Residual vectors
        ax2 = axes[1]
        
        # Reproject
        projected = self.project_catalog_to_image(
            catalog_coords, 
            orientation['RA'], 
            orientation['DEC'], 
            orientation['ROLL']
        )
        
        residuals = image_coords - projected
        residual_mags = np.linalg.norm(residuals, axis=1)
        
        # Quiver plot
        quiver = ax2.quiver(
            image_coords[:, 0], image_coords[:, 1],
            residuals[:, 0], residuals[:, 1],
            residual_mags,
            scale=1, scale_units='xy', angles='xy',
            cmap='hot', alpha=0.8, width=0.003
        )
        
        ax2.scatter(image_coords[:, 0], image_coords[:, 1],
                   c='blue', s=40, alpha=0.4, label='Observed', zorder=1)
        ax2.scatter(projected[:, 0], projected[:, 1],
                   c='green', s=40, marker='+', alpha=0.6, label='Projected', zorder=2)
        
        ax2.set_xlabel('X (pixels)', fontsize=12, fontweight='bold')
        ax2.set_ylabel('Y (pixels)', fontsize=12, fontweight='bold')
        
        title = (f'Residual Vectors\n'
                f'RMS: {orientation["rms_residual_px"]:.2f} px '
                f'({orientation["rms_residual_arcsec"]:.2f}"), '
                f'Max: {orientation["max_residual_px"]:.2f} px')
        ax2.set_title(title, fontsize=14, fontweight='bold')
        ax2.legend(fontsize=10)
        ax2.invert_yaxis()
        ax2.grid(True, alpha=0.3, linestyle='--')
        ax2.set_aspect('equal', adjustable='box')
        
        cbar = plt.colorbar(quiver, ax=ax2, pad=0.02)
        cbar.set_label('Residual (pixels)', fontsize=10, fontweight='bold')
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=200, bbox_inches='tight')
        print(f"      Saved plot: {output_path}")
        plt.close()
    
    def run_pipeline(self, centroid_file: str, output_dir: str, 
                    sensor_name: str) -> Optional[Dict]:
        """
        Execute complete star matching pipeline.
        
        Steps:
        A1: Generate/download catalog
        A2: Create catalog patterns
        A3: Create image patterns and match
        A4: Verify matches
        A5: Estimate orientation and save results
        
        Returns:
        --------
        results : dict or None
            Complete results including orientation and matches
        """
        print(f"\n{'='*80}")
        print(f"STAR MATCHING PIPELINE: {sensor_name}")
        print(f"{'='*80}")
        
        # Load centroids
        print(f"\nLoading centroid data...")
        centroid_data = pd.read_csv(centroid_file)
        print(f"  Loaded {len(centroid_data)} centroids")
        
        # A1: Generate catalog
        self.generate_synthetic_catalog(n_stars=400)
        
        # A2: Create catalog patterns
        print(f"\n[A2] Creating catalog patterns...")
        catalog_coords_sky = self.catalog[['RA', 'DEC']].values
        catalog_coords_img = self.project_catalog_to_image(
            catalog_coords_sky, self.ra_deg, self.dec_deg, 0.0
        )
        
        # Filter to field
        margin = 0.3
        in_field = (
            (catalog_coords_img[:, 0] >= -margin * self.image_width) &
            (catalog_coords_img[:, 0] <= (1 + margin) * self.image_width) &
            (catalog_coords_img[:, 1] >= -margin * self.image_height) &
            (catalog_coords_img[:, 1] <= (1 + margin) * self.image_height)
        )
        catalog_coords_filtered = catalog_coords_img[in_field]
        print(f"     Filtered to {len(catalog_coords_filtered)} stars in field")
        
        catalog_patterns = self.create_triangle_patterns(
            catalog_coords_filtered, n_patterns=200, name="catalog patterns"
        )
        
        # A3: Create image patterns and match
        print(f"\n[A3] Creating image patterns and matching...")
        image_coords = centroid_data[['x_centroid', 'y_centroid']].values
        image_patterns = self.create_triangle_patterns(
            image_coords, n_patterns=200, name="image patterns"
        )
        
        pattern_matches = self.match_patterns(image_patterns, catalog_patterns)
        
        if not pattern_matches:
            print(f"\nERROR: No pattern matches found!")
            return None
        
        # A4: Verify
        image_matched, catalog_matched = self.verify_matches(
            pattern_matches, centroid_data
        )
        
        if image_matched is None or len(image_matched) < 4:
            print(f"\nERROR: Insufficient verified matches!")
            return None
        
        # A5: Estimate orientation
        orientation = self.estimate_orientation(image_matched, catalog_matched)
        
        # Print results
        print(f"\n{'='*80}")
        print(f"RESULTS")
        print(f"{'='*80}")
        print(f"\nOptimal Orientation:")
        print(f"  RA:   {orientation['RA']:11.6f}° (initial: {self.ra_deg:11.6f}°, Δ={orientation['RA']-self.ra_deg:+.6f}°)")
        print(f"  DEC:  {orientation['DEC']:11.6f}° (initial: {self.dec_deg:11.6f}°, Δ={orientation['DEC']-self.dec_deg:+.6f}°)")
        print(f"  ROLL: {orientation['ROLL']:11.6f}°")
        
        print(f"\nCD Matrix (degrees/pixel):")
        cd = orientation['cd_matrix']
        print(f"  [{cd[0,0]:+.6e}  {cd[0,1]:+.6e}]")
        print(f"  [{cd[1,0]:+.6e}  {cd[1,1]:+.6e}]")
        
        print(f"\nMatch Quality:")
        print(f"  Number of matches:  {orientation['n_matches']}")
        print(f"  Mean residual:      {orientation['mean_residual_px']:.3f} pixels")
        print(f"  Median residual:    {orientation['median_residual_px']:.3f} pixels")
        print(f"  RMS residual:       {orientation['rms_residual_px']:.3f} pixels ({orientation['rms_residual_arcsec']:.3f} arcsec)")
        print(f"  Std deviation:      {orientation['std_residual_px']:.3f} pixels")
        print(f"  Max residual:       {orientation['max_residual_px']:.3f} pixels")
        
        # Save outputs
        matched_csv = f"{output_dir}/matched_{sensor_name}.csv"
        plot_png = f"{output_dir}/plot_{sensor_name}.png"
        
        self.save_matched_data(centroid_data, image_matched, catalog_matched, matched_csv)
        self.plot_results(centroid_data, image_matched, catalog_matched, orientation, plot_png)
        
        results = {
            'sensor': sensor_name,
            'orientation': orientation,
            'n_matches': len(image_matched),
            'image_coords': image_matched,
            'catalog_coords': catalog_matched
        }
        
        print(f"\n{'='*80}")
        print(f"PIPELINE COMPLETE: {sensor_name}")
        print(f"{'='*80}\n")
        
        return results


def process_all_sensors():
    """
    Process all three centroid files.
    """
    sensors = {
        'star_centroids_1': {
            'image_size': (6576, 4384),
            'focal_length_mm': 980,
            'pixel_pitch_um': 5.5,
            'ra_deg': 214.958,
            'dec_deg': -47.909,
            'fov_x_deg': 2.030,
            'fov_y_deg': 1.409
        },
        'star_centroids_2': {
            'image_size': (1024, 1024),
            'focal_length_mm': 4912,
            'pixel_pitch_um': 3.6,
            'ra_deg': 101.812,
            'dec_deg': -48.565,
            'fov_x_deg': 0.429,
            'fov_y_deg': 0.429
        },
        'star_centroids_3': {
            'image_size': (13400, 9528),
            'focal_length_mm': 393,
            'pixel_pitch_um': 3.45,
            'ra_deg': 17.283,
            'dec_deg': 13.532,
            'fov_x_deg': 6.732,
            'fov_y_deg': 4.789
        }
    }
    
    output_dir = "/mnt/user-data/outputs"
    all_results = []
    
    for sensor_name, config in sensors.items():
        print(f"\n\n{'#'*80}")
        print(f"# PROCESSING: {sensor_name}")
        print(f"{'#'*80}\n")
        
        try:
            matcher = StarMatcher(
                focal_length_mm=config['focal_length_mm'],
                pixel_pitch_um=config['pixel_pitch_um'],
                image_size=config['image_size'],
                ra_deg=config['ra_deg'],
                dec_deg=config['dec_deg'],
                fov_x_deg=config['fov_x_deg'],
                fov_y_deg=config['fov_y_deg']
            )
            
            centroid_file = f"Assignment Star Tracker/{sensor_name}.csv"
            results = matcher.run_pipeline(centroid_file, output_dir, sensor_name)
            
            if results:
                all_results.append({
                    'sensor': sensor_name,
                    'success': True,
                    'n_matches': results['orientation']['n_matches'],
                    'rms_error_px': results['orientation']['rms_residual_px'],
                    'rms_error_arcsec': results['orientation']['rms_residual_arcsec'],
                    'ra': results['orientation']['RA'],
                    'dec': results['orientation']['DEC'],
                    'roll': results['orientation']['ROLL']
                })
            else:
                all_results.append({
                    'sensor': sensor_name,
                    'success': False
                })
                
        except Exception as e:
            print(f"\nERROR processing {sensor_name}: {e}")
            import traceback
            traceback.print_exc()
            all_results.append({
                'sensor': sensor_name,
                'success': False,
                'error': str(e)
            })
    
    # Summary
    print(f"\n\n{'='*80}")
    print(f"PROCESSING SUMMARY")
    print(f"{'='*80}\n")
    
    for result in all_results:
        print(f"{result['sensor']}:")
        if result['success']:
            print(f"  ✓ SUCCESS")
            print(f"  Matches: {result['n_matches']}")
            print(f"  RMS error: {result['rms_error_px']:.3f} px ({result['rms_error_arcsec']:.3f} arcsec)")
            print(f"  RA: {result['ra']:.6f}°, DEC: {result['dec']:.6f}°, ROLL: {result['roll']:.6f}°")
        else:
            print(f"  ✗ FAILED")
            if 'error' in result:
                print(f"  Error: {result['error']}")
        print()


if __name__ == "__main__":
    process_all_sensors()
