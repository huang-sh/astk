import hashlib
import pickle
from pathlib import Path

from astk.types import *
from .AS_event import make_events
from  .gtf_parse import construct_genome
from astk.suppa.lib.diff_tools import multiple_conditions_analysis as mca
from astk.constant import AS_TYPE


def check_gtf_used(gtf):
    gtf = Path(gtf)
    gtf_hash = hashlib.blake2b()
    gtf_size = gtf.stat().st_size
    gtf_hash.update(str(gtf_size).encode('utf-8'))
    f = gtf.open("rb")
    for _ in range(10):
        gtf_hash.update(f.readline().strip())
    hash_val = gtf_hash.hexdigest()
    return hash_val


def gtf_parse_cache(gtf):
    home = Path.home()
    astk_cache_dir = home / f".astk"
    astk_cache_dir.mkdir(exist_ok=True)
    gtf_hash = check_gtf_used(gtf)
    gtf_cache_path = astk_cache_dir / "gtfParse"
    gtf_cache_path.mkdir(exist_ok=True)
    pkl = (gtf_cache_path / gtf_hash).with_suffix(".pkl")
    return pkl
    

def generate_events(gtf, event_types, output, split):
    import time
    T1 = time.time()

    T21 =time.time()
    print('程序运行时间:%s秒' % ((T21 - T1)))

    genome = construct_genome(gtf)
    T2 =time.time()
    print('程序运行时间:%s秒' % ((T2 - T21)))
    
    if event_types == "ALL":
        event_types =  ['SE', "A5", "A3", "MX", "RI", 'AF', 'AL']
    else:
        event_types = [event_types]

    make_events(output, genome, event_types, split)

    T3 =time.time()
    print('程序运行时间:%s秒' % ((T3 - T1)))


FilePath = TypeVar('FilePath', str, Path)


def diff_splice(psi_files: Sequence[FilePath],
                exp_files: Sequence[FilePath],
                reference: FilePath,
                output: FilePath,
                method: str
    ):
    mca(method, psi_files, exp_files, reference, 1000, 0, 
        False, True, 0.05, True, False, False, 0, 0, str(output))                 
