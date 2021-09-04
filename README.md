### CGRdb

Chemical cartridge for reactions and molecules.

### Install Postgres
Use docker-compose.yaml to set up service or do installation by yourself.
Password, user and port can be changed in docker-compose.yaml

    docker-compose up -d

### Library installation
Set up your virtual environment, activate it, go into CGRdb directory and do:

    pip install -U -e .

### Usage
Please use jupyter notebooks from examples folder to upload data and test functionality, other
documentation in progress.

### Data
All data required for examples are in dataset folder.

### Reproducibility of results from paper
Results can be reproduced with Test class from CGRdb.tests.
It contains two methods Test.test_reactions() and Test.test_mols(),
for testing CGR results use Test.test_mols(cgr=True)
