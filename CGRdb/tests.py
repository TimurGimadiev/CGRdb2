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
3: 'CN[.>-]C(C[.>-]C1OCCc2c3c(cccc3)sc12)[->.]Cl',
4: 'c12C(=O)N(C(c3cccc(ccc1)c23)=O)[.>-]C([.>-]N4CCN(CC4)c5ccccc5)[=>.]O',
5: 'O=C(C)CCC([.>-]n1ccnc1)[->.]Br'
               }

reaction_queries = {
1:'[CH2:5]([OH:11])[c:6]1[s:10][cH:9][n:8][cH:7]1.[SH:4][CH2:3][CH2:2][NH2:1]>>[NH2:1][CH2:2][CH2:3][S:4][CH2:5][c:6]1[s:10][cH:9][n:8][cH:7]1',
2:'[O:17]1[CH2:18][C:3](=[CH:4][c:5]2[c:16]1[cH:15][cH:14][c:7]([Br:19])[cH:6]2)[CH:2]=[O:1].[cH:9]1[cH:10][cH:11][cH:12][cH:13][c:8]1[B:21]([OH:20])[OH:22]>>[cH:10]1[cH:9][c:8]([cH:13][cH:12][cH:11]1)-[c:7]2[cH:14][cH:15][c:16]3[O:17][CH2:18][C:3](=[CH:4][c:5]3[cH:6]2)[CH:2]=[O:1]',
3:'[CH2:3]([CH3:4])[Cl:18].[NH2:2][CH3:1].[c:10]12[c:9]([c:17]3[c:12]([cH:13][cH:14][cH:15][cH:16]3)[s:11]1)[CH2:8][CH2:7][O:6][CH2:5]2>>[c:12]12[c:17]([c:9]3[c:10]([CH:5]([CH2:4][CH2:3][NH:2][CH3:1])[O:6][CH2:7][CH2:8]3)[s:11]1)[cH:16][cH:15][cH:14][cH:13]2',
4: '[O:29]=[CH2:16].[cH:22]1[cH:23][cH:24][cH:25][cH:26][c:21]1[N:20]2[CH2:27][CH2:28][NH:17][CH2:18][CH2:19]2.[cH:4]1[cH:5][cH:6][c:7]2[cH:8][cH:9][cH:10][c:11]3[C:13](=[O:14])[NH:15][C:2](=[O:1])[c:3]1[c:12]23>>[c:3]12[C:2](=[O:1])[N:15]([C:13]([c:11]3[cH:10][cH:9][cH:8][c:7]([cH:6][cH:5][cH:4]1)[c:12]23)=[O:14])[CH2:16][N:17]4[CH2:18][CH2:19][N:20]([CH2:27][CH2:28]4)[c:21]5[cH:26][cH:25][cH:24][cH:23][cH:22]5',
5: '[CH2:4]([CH2:5][CH2:6][Br:12])[C:2](=[O:3])[CH3:1].[n:7]1[cH:8][cH:9][nH:10][cH:11]1>>[n:7]1([cH:11][n:10][cH:9][cH:8]1)[CH2:6][CH2:5][CH2:4][C:2](=[O:3])[CH3:1]'

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
