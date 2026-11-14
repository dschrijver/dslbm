import h5py
import os
from pathlib import Path
from natsort import natsorted

def read_data():
    script_path = str(os.path.dirname(os.path.abspath(__file__)))
    data_path = script_path.rsplit("/", 1)[0]
    data_path = Path(data_path)
    data_files = data_path.glob("data_*.h5")
    data_files = natsorted(data_files, key=lambda x: x.name)
    return data_files