import sys
from rdkit import Chem
from decimal import Decimal
from rdkit.Chem import rdMolDescriptors
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
import argparse

#default values:
HBD_t = 1
MW_t = 25
RB_t = 1
HBA_t = 2
ClogP_t = float(Decimal(1))#1.5
tanimoto_t = Decimal('0.75')

class Molecule(object):    
	def __init__(self, mol_input, input_type = "smile"):
		mol = Chem.MolFromSmiles(mol_input)
		self.HBD = self.compute_HBD(mol)
		self.Wt = self.compute_Wt(mol)
		self.NumRotatableBonds = self.compute_NumRotatableBonds(mol)
		self.HBA = self.compute_HBA(mol)
		self.logP = self.compute_logP(mol)
		self.FingerPrint = self.compute_fingerprint(mol)
		
	def compute_HBD(self, mol_input):
		return rdMolDescriptors.CalcNumHBD(mol_input)
	
	def compute_Wt(self, mol_input):
		return rdMolDescriptors.CalcExactMolWt(mol_input)

	def compute_NumRotatableBonds(self, mol_input):
		return rdMolDescriptors.CalcNumRotatableBonds(mol_input)

	def compute_HBA(self, mol_input):
		return rdMolDescriptors.CalcNumHBA(mol_input)

	def compute_logP(self, mol_input):
		logP, mr = rdMolDescriptors.CalcCrippenDescriptors(mol_input)
		return logP
	
	def compute_fingerprint(self, mol_input):
		return FingerprintMols.FingerprintMol(mol_input)
	
	def compute_similarity(self, mol_input):
		return DataStructs.FingerprintSimilarity(self.FingerPrint, mol_input.FingerPrint, metric=DataStructs.TanimotoSimilarity)
			
	def is_decoy(self, mol_input):
		if self.HBD - HBD_t <= mol_input.HBD <= self.HBD + HBD_t:
			if self.Wt - MW_t <= mol_input.Wt <= self.Wt + MW_t:
				if self.NumRotatableBonds - RB_t <= mol_input.NumRotatableBonds <= self.NumRotatableBonds + RB_t:
					if self.HBA - HBA_t <= mol_input.HBA <= self.HBA + HBA_t :
						if self.logP - ClogP_t <= mol_input.logP <= self.logP + ClogP_t:
							if  self.compute_similarity(mol_input) <= tanimoto_t:
								return True
		return False

	def find_decoys(self, input_file, output_file, output_type ="sdf", input_type ="sdf"):
		lista = []
		if input_type == "sdf":
			suppl = Chem.SDMolSupplier(input_file)
			for mol in suppl:
				try: 
					smi = mol.GetProp("SMILES")
					m = Molecule(smi, input_type="smile")
					if self.is_decoy(m):
						lista.append(mol)
				except KeyError:
					pass
				sys.stdout.flush()
		elif input_type == "txt":
			f = open(input_file, "r")
			line = f.readline()
			while line:
				smi = line.split(",")
				try: 
					m = Molecule(smi[1], input_type="smile")
					if self.is_decoy(m):
						lista.append(line)
				except Exception as e:
					pass
				line = f.readline()
				sys.stdout.flush()
			f.close()
		if output_type == "sdf":
			self.write_output(lista, output_file, output_type=input_type)
		elif output_type == "txt":
			self.write_output(lista, output_file, output_type=input_type)
	
	def write_output(self, lista, output_file, output_type="sdf"):
		if output_type == "sdf":
			w = Chem.SDWriter(output_file)
			for mol in lista:
				w.write(mol)
		elif output_type == "txt":
			f = open(output_file, "w")
			for mol in lista:
				f.write(mol)
			f.close()
		print "FET!"

#Program
def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument("-s", "--smile", required=True, help="reference molecule")
	parser.add_argument("-i", "--type_input", choices=['smile', 'sdf', 'txt'], help="input type", required=True,)
	parser.add_argument("-I", "--_input", help="name input which should contain the data""data input", required=True)
	parser.add_argument("-o", "--type_output", choices=['smile', 'sdf', 'txt'], help="output type", required=False, default=None)
	parser.add_argument("-O", "--_output", help="name output which will contain the data""data output", required=False, default=None)
	opts = parser.parse_args()

	if opts.type_input == "smile":
		return opts
		
	elif not (opts._input) or not (opts._output):
		parser.error("Incorrect number of arguments, required an input file and an output file")
		
	elif not (opts._input.endswith(opts.type_input)) or not (opts._output.endswith(opts.type_output)):
		if not (opts._input.endswith(opts.type_input)) and not (opts._output.endswith(opts.type_output)):
			parser.error("Incorrect type input and type output")
		elif not (opts._input.endswith(opts.type_input)):
			parser.error("Incorrect type input")
		else:
			parser.error("Incorrect type output")
			
	else:
		return opts

def main():
	opts = parse_args()
	#.txt files
	if opts.type_input == "txt":
		smile = Molecule(opts.smile)
		smile.find_decoys(opts._input, opts._output, output_type=opts.type_output, input_type="txt")
	#.sdf files
	elif opts.type_input == "sdf":
		smile = Molecule(opts.smile)
		smile.find_decoys(opts._input, opts._output, output_type=opts.type_output, input_type="sdf")
	#smile to compare
	elif opts.type_input == "smile":
		smile = Molecule(opts.smile)
		_input = Molecule(opts._input)
		print smile.is_decoy(_input)

if __name__ == "__main__":
	main()


