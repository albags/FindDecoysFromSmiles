# FindDecoysFromSmiles

Find decoys from canonical smiles using Fingerprints and Tanimoto Index.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 

### Prerequisites

* [Anaconda Python 2.7](https://www.anaconda.com/download/)
* [Rdkit](http://www.rdkit.org/docs/Install.html)

### Installing & running
Clone the project on your local machine

`$ git clone https://github.com/albags/FindDecoysFromSmiles.git`

`$ cd FindDecoysFromSmiles/`

Go to your terminal and run the python file.

For knowing the arguments
`$ python Molecule_decoys.py --help`

Exemples of arguments:
`$ python Molecule_decoys.py -s 'C1=CC=CN=C1' -i smile -I 'c1cccnc1'`

`$ python Molecule_decoys.py -s 'C1=CC=CN=C1' -i txt -I 110000smiles.txt -o txt -O output.txt`

## Licence
ProjectTargetPredictionApp is released under the [MIT License](LICENSE).
