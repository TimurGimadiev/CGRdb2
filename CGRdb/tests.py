from CGRdb.database import db
from pony.orm import db_session
from collections import namedtuple, defaultdict
from .database import Molecule, CGR, Reaction
from time import time
from CGRtools import smiles, ReactionContainer, MoleculeContainer
from typing import Optional

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

cgr_queries = {
1: 'NCCS[.>-]C(c1scnc1)[->.]O',
2: 'OB(O)[->.]c1(ccccc1)[.>-]c2([->.]Br)cc3C=C(C=O)COc3cc2',
3: 'CC1(C)N[->.]C1[.>-]SCCO',
4: 'CC(C)([.>-]OCC)[=>-]CCCC(C)CCO[->.]C(C)=O'
               }

reaction_queries = {
1:'[CH2:5]([OH:11])[c:6]1[s:10][cH:9][n:8][cH:7]1.[SH:4][CH2:3][CH2:2][NH2:1]>>[NH2:1][CH2:2][CH2:3][S:4][CH2:5][c:6]1[s:10][cH:9][n:8][cH:7]1',
2:'[O:17]1[CH2:18][C:3](=[CH:4][c:5]2[c:16]1[cH:15][cH:14][c:7]([Br:19])[cH:6]2)[CH:2]=[O:1].[cH:9]1[cH:10][cH:11][cH:12][cH:13][c:8]1[B:21]([OH:20])[OH:22]>>[cH:10]1[cH:9][c:8]([cH:13][cH:12][cH:11]1)-[c:7]2[cH:14][cH:15][c:16]3[O:17][CH2:18][C:3](=[CH:4][c:5]3[cH:6]2)[CH:2]=[O:1]',
3:'[CH2:8]([OH:9])[CH2:7][SH:6].[CH3:1][C:2]1([CH3:3])[CH2:5][NH:4]1>>[CH3:1][C:2]([CH3:3])([NH2:4])[CH2:5][S:6][CH2:7][CH2:8][OH:9]',
4: '[CH2:2]([OH:3])[CH3:1].[CH3:5][C:4]([CH3:6])=[CH:7][CH2:8][CH2:9][CH:10]([CH3:11])[CH2:12][CH2:13][O:14][C:16]([CH3:15])=[O:17]>>[CH3:5][C:4]([CH3:6])([O:3][CH2:2][CH3:1])[CH2:7][CH2:8][CH2:9][CH:10]([CH2:12][CH2:13][OH:14])[CH3:11]'
                    }

def spent_time(start):
    return time()-start

class Test:
    def __init__(self, provider='postgres', user='postgres', host='localhost', password="example", database='test',
                 port=5432, test_sequences=(1, 10, 100, 1000, -1)):
        db.bind(provider=provider, user=user, host=host, password=password, database=database,
                port=port)
        db.generate_mapping(create_tables=False)
        self.test_sequences = test_sequences

    def test_mols(self, smi: Optional[list] = None, cgr=False):
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
        if smi is None:
            if cgr:
                smi = cgr_queries.items()
            else:
                smi = mol_queries.items()
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

    def query_cgr(self, mol: MoleculeContainer, n=1, similarity=True, ordered=True):
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

    def test_reactions(self, smi: Optional[list] = None):
        """
        will test speed of requests:
            - get first
            - get first 10
            - get first 100
            - get first 1000
            - get all
        :return:
        """
        if smi is None:
            smi = reaction_queries.values()
        results = defaultdict(list)
        for i in smi:
            mol = smiles(i)
            for n in self.test_sequences:
                with db_session:
                    results["similarity_ordered"].append(self.query_reaction(mol, n=n, ordered=True))
                with db_session:
                    results["substructure_ordered"].append(self.query_reaction(mol, n=n,
                                                                               similarity=False, ordered=True))
                with db_session:
                    results["similarity_ordered_checkRC"].append(self.query_reaction(mol, n=n, ordered=True,
                                                                                     mapping=True))
                with db_session:
                    results["substructure_ordered_checkRC"].append(self.query_reaction(mol, n=n, similarity=False,
                                                                                       ordered=True, mapping=True))
                with db_session:
                    results["similarity_unordered"].append(self.query_reaction(mol, n=n, ordered=False))
                with db_session:
                    results["substructure_unordered"].append(self.query_reaction(mol, n=n,
                                                                                 similarity=False, ordered=False))
                with db_session:
                    results["similarity_unordered_checkRC"].append(self.query_reaction(mol, n=n, ordered=False, mapping=True))
                with db_session:
                    results["substructure_unordered_checkRC"].append(self.query_reaction(mol, n=n, similarity=False,
                                                                                 ordered=False, mapping=True))
                print(f"{n} sequence finished for {i}")
        return results

    def query_reaction(self, reaction: ReactionContainer, n=1, similarity: bool = True, ordered: bool = True,
                       fix_roles: bool = True, mapping: bool = False, request_only: bool = False):
        with db_session():
            results_found = 0
            if similarity:
                start = time()
                res = Reaction.similars(reaction, ordered=ordered, fix_roles=fix_roles,
                                        mapping=mapping, request_only=request_only)
            else:
                start = time()
                res = Reaction.substructures(reaction, ordered=ordered, fix_roles=fix_roles,
                                             mapping=mapping, request_only=request_only)
            for i, r in enumerate(res, 1):
                #r[0].unload()
                results_found += 1
                if i == n:
                    break
            t = spent_time(start)
        return Result(smiles=str(reaction), spent_time=t, results_requested=n, results_found=results_found)
