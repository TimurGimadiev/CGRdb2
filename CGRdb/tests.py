from CGRdb.database import db
from pony.orm import db_session
from collections import namedtuple, defaultdict
from .database import Molecule
from time import time
from CGRtools import smiles

Result = namedtuple("Result", ["smiles", "spent_time", "results_requested", "results_found"])

Report = namedtuple("Report", ["request_type", "results", "average_time", "min", "max"])

mol_queries = {"benzene": "c1ccccc1",
               "atropine": "CN1[C@H]2CC[C@@H]1C[C@H](OC(=O)C(CO)c1ccccc1)C2",
               "aspirin": "CC(=O)Oc1ccccc1C(=O)O",
               "ethanol": "CCO", "paracetamol": "CC(=O)Nc1ccc(O)cc1",
               "haloperidol": "O=C(CCCN1CCC(O)(c2ccc(Cl)cc2)CC1)c1ccc(F)cc1",
               "cinnoline": "c1ccc2nnccc2c1",
               "glycine": "NCC(=O)O",
               "carbenicillin": "CC1(C)S[C@@H]2[C@H](NC(=O)C(C(=O)O)c3ccccc3)C(=O)N2[C@H]1C(=O)O",
               "phenol": "Oc1ccccc1"}

def spent_time(start):
    return time()-start

class Test:
    def __init__(self, provider='postgres', user='postgres', host='localhost', password="example", database='test',
                 port=5432, test_sequences=(1, 10, 100, 1000, -1)):
        db.bind(provider=provider, user=user, host=host, password=password, database=database,
                port=port)
        db.generate_mapping(create_tables=False)
        self.test_sequences = test_sequences

    def test_mols(self, smi: list):
        """
        will test speed of requests:
            - get first
            - get first 10
            - get first 100
            - get first 1000
            - get all
        :return:
        """
        results = defaultdict(list)
        for i in smi:
            mol = smiles(i)
            for n in self.test_sequences:
                results["similarity_ordered"].append(self.query_mol(mol, n=n, ordered=True))
                results["substructure_ordered"].append(self.query_mol(mol, n=n, similarity=False, ordered=True))
                results["similarity_unordered"].append(self.query_mol(mol, n=n, ordered=False))
                results["substructure_unordered"].append(self.query_mol(mol, n=n, similarity=False, ordered=False))
                print(f"{n} sequence finished for {i}")
        return results

    def query_mol(self, mol, n=1, similarity=True, ordered=True, no_graph=True):
        with db_session():
            results_found = 0
            if similarity:
                start = time()
                res = Molecule.similars(mol, ordered=ordered, no_graph=no_graph)
            else:
                start = time()
                res = Molecule.substructres(mol, ordered=ordered, no_graph=no_graph)
            for i, _ in enumerate(res, 1):
                results_found += 1
                if i == n:
                    break
            t = spent_time(start)
        return Result(smiles=str(mol), spent_time=t, results_requested=n, results_found=results_found)
