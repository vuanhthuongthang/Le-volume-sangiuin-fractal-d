import numpy as np
from scipy.ndimage import distance_transform_edt
from skimage.morphology import skeletonize
from scipy.optimize import curve_fit
from scipy.integrate import quad

class FractalVascularVolume:
    def __init__(self, binary_volume: np.ndarray, voxel_spacing: float):
        """
        Initialise l'analyseur de volume vasculaire.
        
        :param binary_volume: Masque 3D numpy (booléen ou 0/1) du réseau scanné.
        :param voxel_spacing: Résolution isométrique du voxel en mm (ex: 0.5 mm).
        """
        self.volume = binary_volume.astype(bool)
        self.voxel_spacing = voxel_spacing
        self.voxel_volume = voxel_spacing ** 3
        
        self.skeleton = None
        self.radii_mm = None
        self.fractal_dimension = None
        self.fit_params = None
        
    def extract_skeleton_and_radii(self):
        """
        Étape 1 & 2 : Squelettisation 3D et calcul des rayons locaux (Transformée de Distance).
        """
        print("Calcul de la transformée de distance (EDT)...")
        # L'EDT donne la distance euclidienne de chaque pixel vasculaire au bord le plus proche
        edt = distance_transform_edt(self.volume)
        
        print("Extraction du squelette topologique 3D...")
        # Squelettisation 3D préservant la topologie
        self.skeleton = skeletonize(self.volume)
        
        # Extraction des rayons aux points exacts du squelette, convertis en mm
        radii_voxels = edt[self.skeleton]
        self.radii_mm = radii_voxels * self.voxel_spacing
        
        print(f"Squelette extrait : {len(self.radii_mm)} segments locaux trouvés.")
        return self.radii_mm

    def fit_power_law(self, bins=40):
        """
        Étape 3 : Fitter la loi de puissance V(r) proportionnelle à r^(3-d).
        """
        if self.radii_mm is None:
            self.extract_skeleton_and_radii()
            
        # Création de l'histogramme des rayons du squelette
        # N(r) est le nombre de voxels du squelette ayant un rayon r (représente la longueur L(r))
        counts, bin_edges = np.histogram(self.radii_mm, bins=bins)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        dr = bin_edges[1] - bin_edges[0]
        
        # Le volume d'un segment cylindrique local est V = pi * r^2 * L
        # L(r) = counts * voxel_spacing
        lengths = counts * self.voxel_spacing
        volumes = np.pi * (bin_centers ** 2) * lengths
        
        # Densité volumique v(r) = dV / dr
        v_r = volumes / dr
        
        # Filtrage des bins vides pour le passage au logarithme
        valid = v_r > 0
        r_valid = bin_centers[valid]
        v_valid = v_r[valid]
        
        # Régression linéaire dans l'espace Log-Log : ln(v(r)) = ln(C) + p * ln(r)
        # où p = 3 - d
        def log_model(r_log, ln_C, p):
            return ln_C + p * r_log
            
        popt, _ = curve_fit(log_model, np.log(r_valid), np.log(v_valid))
        
        ln_C, p = popt
        self.fit_params = {'C': np.exp(ln_C), 'p': p}
        self.fractal_dimension = 3 - p
        
        print(f"Régression terminée : Dimension fractale locale d = {self.fractal_dimension:.3f}")
        return self.fit_params

    def calculate_total_volume(self, r_min=0.008, r_max_scanner=0.5):
        """
        Étape 4 : Intégration analytique pour l'extrapolation microvasculaire.
        
        :param r_min: Rayon capillaire limite (par défaut 8 µm = 0.008 mm).
        :param r_max_scanner: Limite de résolution de certitude du scanner (en mm).
        """
        if self.fit_params is None:
            self.fit_power_law()
            
        C = self.fit_params['C']
        p = self.fit_params['p']
        
        # Densité volumique ajustée v(r) = C * r^p
        def volume_density(r):
            return C * (r ** p)
            
        # Intégration analytique sur la partie manquante (Microcirculation)
        v_micro, error = quad(volume_density, r_min, r_max_scanner)
        
        # Volume macroscopique directement compté depuis le masque
        v_macro_scan = np.sum(self.volume) * self.voxel_volume
        
        v_total_eff = v_macro_scan + v_micro
        
        print("\n--- Bilan Volumique (mm³) ---")
        print(f"Volume Macro scanné (V_scan) : {v_macro_scan:.2f} mm³")
        print(f"Volume Micro extrapolé       : {v_micro:.2f} mm³")
        print(f"Volume Sanguin Total Effectif: {v_total_eff:.2f} mm³")
        
        return v_total_eff

# ==========================================
# Bloc d'exécution et de test (Génération d'un volume factice)
# ==========================================
if __name__ == "__main__":
    # Génération d'un tenseur 3D binaire factice (un simple tube cylindrique pour le test)
    # Dans un cas réel, on chargerait un NIfTI ou DICOM : volume = nibabel.load('scan.nii.gz').get_fdata() > seuil
    print("Génération d'un volume vasculaire 3D factice (100x100x100)...")
    grid = np.mgrid[-50:50, -50:50, -50:50]
    # Création d'une branche principale (cylindre)
    cylinder = (grid[0]**2 + grid[1]**2) <= 15**2
    # Ajout d'une bifurcation plus petite
    branch = ( (grid[1]-20)**2 + grid[2]**2 <= 7**2 ) & (grid[0] > 0)
    mock_volume = cylinder | branch

    # Instanciation de l'analyseur (résolution de 0.5 mm par voxel)
    analyzer = FractalVascularVolume(mock_volume, voxel_spacing=0.5)
    
    # Exécution du pipeline complet
    analyzer.extract_skeleton_and_radii()
    analyzer.fit_power_law()
    analyzer.calculate_total_volume(r_min=0.008, r_max_scanner=0.5)
