from CGRdb.database import db
from pony.orm import db_session
from collections import namedtuple, defaultdict
from .database import Molecule, CGR
from time import time
from CGRtools import smiles

Result = namedtuple("Result", ["smiles", "spent_time", "results_requested", "results_found"])

Report = namedtuple("Report", ["request_type", "results", "average_time", "min", "max"])

mol_queries = {"benzene": "c1ccccc1",
               "atropine": "CN1[C@H]2CC[C@@H]1C[C@H](OC(=O)C(CO)c1ccccc1)C2",
               "aspirin": "CC(=O)Oc1ccccc1C(=O)O",
               "ethanol": "CCO",
               "paracetamol": "CC(=O)Nc1ccc(O)cc1",
               "haloperidol": "O=C(CCCN1CCC(O)(c2ccc(Cl)cc2)CC1)c1ccc(F)cc1",
               "cinnoline": "c1ccc2nnccc2c1",
               "glycine": "NCC(=O)O",
               "carbenicillin": "CC1(C)S[C@@H]2[C@H](NC(=O)C(C(=O)O)c3ccccc3)C(=O)N2[C@H]1C(=O)O",
               "phenol": "Oc1ccccc1"}

cgr_queries = {1: 'NCCS[.>-]C(c1scnc1)[->.]O'}

reaction_queries = {1:'[CH2:5]([OH:11])[c:6]1[s:10][cH:9][n:8][cH:7]1.[SH:4][CH2:3][CH2:2][NH2:1]>>[NH2:1][CH2:2][CH2:3][S:4][CH2:5][c:6]1[s:10][cH:9][n:8][cH:7]1'}

def spent_time(start):
    return time()-start

class Test:
    def __init__(self, provider='postgres', user='postgres', host='localhost', password="example", database='test',
                 port=5432, test_sequences=(1, 10, 100, 1000, -1)):
        db.bind(provider=provider, user=user, host=host, password=password, database=database,
                port=port)
        db.generate_mapping(create_tables=False)
        self.test_sequences = test_sequences

    def test_mols(self, smi: list, cgr=False):
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
        if cgr:
            query = self.query_cgr
        else:
            query = self.query_mol
        for i in smi:
            mol = smiles(i)
            for n in self.test_sequences:
                with db_session:
                    results["similarity_ordered"].append(query(mol, n=n, ordered=True))
                with db_session:
                    results["substructure_ordered"].append(query(mol, n=n, similarity=False, ordered=True))
                with db_session:
                    results["similarity_unordered"].append(query(mol, n=n, ordered=False))
                with db_session:
                    results["substructure_unordered"].append(query(mol, n=n, similarity=False, ordered=False))
                print(f"{n} sequence finished for {i}")
        return results

    def query_mol(self, mol, n=1, similarity=True, ordered=True):
        with db_session():
            results_found = 0
            if similarity:
                start = time()
                res = Molecule.similars(mol, ordered=ordered)
            else:
                start = time()
                res = Molecule.substructres(mol, ordered=ordered)
            for i, r in enumerate(res, 1):
                r[1].unload()
                results_found += 1
                if i == n:
                    break
            t = spent_time(start)
        return Result(smiles=str(mol), spent_time=t, results_requested=n, results_found=results_found)

    def query_cgr(self, mol, n=1, similarity=True, ordered=True):
        with db_session():
            results_found = 0
            if similarity:
                start = time()
                res = CGR.similars(mol, ordered=ordered)
            else:
                start = time()
                res = CGR.substructres(mol, ordered=ordered)
            for i, r in enumerate(res, 1):
                r[0].unload()
                results_found += 1
                if i == n:
                    break
            t = spent_time(start)
        return Result(smiles=str(mol), spent_time=t, results_requested=n, results_found=results_found)
