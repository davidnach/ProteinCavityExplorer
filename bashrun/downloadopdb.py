import requests
import sys
import os.path
from subprocess import Popen

save_path = './pdbs'

pid = sys.argv[1]
completename = os.path.join(save_path, pid +'.pdb')
url = 'https://files.rcsb.org/view/'+ pid +'.pdb'
r = requests.get(url, allow_redirects=True)

if (r.status_code == 200) :
    file1 = open(os.path.join(save_path, pid +'.pdb'), 'wb').write(r.content)

