from pathlib import Path

file_dir = Path(__file__).parent
project_root = file_dir.parent

# ********************************** PATHS ******************************************* #
OUTPUT_DIR = project_root / "output"
EXOMOL_DATA_DIR = None

# ******************************* PROCESSING ***************************************** #
# chunk size for .states files: approx 1,000,000 per 1GB of RAM
STATES_CHUNK_SIZE = 1_000_000
# chunk size for .trans files: roughly 10,000,000 per 1GB of RAM
TRANS_CHUNK_SIZE = 10_000_000

# ****************************** LOCAL CONFIG **************************************** #
# load the local config:
from .config_local import *
