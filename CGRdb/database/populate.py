from .substance import Substance
from .reaction import Reaction
from .config import db
from CGRtools import SMILESRead
from pony.orm import commit, db_session
from CGRtools import smiles
from multiprocessing import Queue, Process
import zipfile
from tqdm import tqdm
from CGRtools.exceptions import InvalidAromaticRing, IncorrectSmiles, ValenceError, MappingError, EmptyMolecule


def upload_smi(filename, provider='postgres', user='postgres', host='localhost', password="example", database='test',
               port=5432, workers=10, reaction=False, verbose=True):

    def worker(q, reaction=reaction):
        db.bind(provider=provider, user=user, host=host, password=password, database=database,
                port=port)
        db.generate_mapping()
        while True:
            data = q.get()
            if data is None:
                break
            for i in data:
                mol, label = i
                mol = smiles(mol)
                try:
                    mol.canonicalize()
                except (InvalidAromaticRing, IncorrectSmiles, ValenceError, MappingError, EmptyMolecule):
                    if verbose:
                        print(mol, "standardization_failed")
                    continue
                if reaction:
                    for _ in range(10):
                        try:
                            with db_session():
                                Reaction(mol)
                                break
                        except Exception:
                            continue
                    else:
                        if verbose:
                            print(mol, "upload_failed")
                else:
                    subs = [(x, None) for x in mol.split()]
                    for _ in range(10):
                        try:
                            with db_session():
                                Substance(subs)
                                break
                        except Exception:
                            continue
                    else:
                        if verbose:
                            print(mol, "upload_failed")
        print("finished")

    q = Queue(maxsize=workers*2)
    num_workers = workers
    pr = [Process(target=worker,
                  args=[q]) for _ in range(num_workers)]
    [p.start() for p in pr]
    if filename.endswith(".zip"):
        with zipfile.ZipFile(filename, 'r') as zip_ref:
            zip_ref.extractall(filename[:-4])
        filename = filename[:-4]
    elif not filename.endswith(".smi"):
        raise ("""file with smiles should have .smi extension, it also can be zipped (.smi.zip). 
                        It should have following structure in each row: SMILES(required) id(optional) 
                        """)
    with open(filename) as f:
        tmp = []
        for n, row in tqdm(enumerate(f)):
            columns = row.rstrip('\n').lstrip().split(" ")
            smi = columns[0]
            idx = None
            if len(columns) > 1:
                if columns[1].startswith("|"):
                    smi += columns[1]
                    if len(columns) > 2:
                        idx = columns[2]
                else:
                    idx = columns[1]
            if not idx:
                idx = n
            tmp.append((smi, idx))
            if len(tmp) == 1000:
                q.put(tmp)
                tmp = []
        else:
            q.put(tmp)
        for i in range(num_workers):
            q.put(None)

