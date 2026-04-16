# baseline_config.py

from pathlib import Path

# Paths
PHI_DIR = Path("PHI")
HMI_DIR = Path("HMI")
OUT_DIR = Path("baseline_outputs")

# Dates to include
ONLY_DATES = {"20221027", "20221028"}

# Matching
MAX_TIME_DIFF_SEC = 600  # 10 minutes

# Smooth PHI-HMI blend
R_INNER = 0.70
R_OUTER = 0.90

# LOS -> Br approximation
MU_MIN = 0.40
ALPHA = 0.80

# Disk usage
DISK_FRACTION = 0.98

# Carrington-style binning
NLAT = 180
NLON = 360