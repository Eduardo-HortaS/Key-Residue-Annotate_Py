import os
import subprocess
import sys
from joblib import Parallel, delayed

# Para chamar um script no bash
subprocess.run(["python", "prepare_fasta_per_domain_tentativa_ruim_parallel_inside_script.py", "--per_dom_json", "per_dom.json", "--resource_dir", "resource_dir", "--output_dir", "output_dir", "--log", "log.log"], check=True)

# Para paralelizar
Parallel(n_jobs=4)(delayed(subprocess.run)(["python", "prepare_fasta_per_domain_tentativa_ruim_parallel_inside_script.py", "--per_dom_json", "per_dom.json", "--resource_dir", "resource_dir", "--output_dir", "output_dir", "--log", "log.log"]) for i in range(4))

# Ver depois