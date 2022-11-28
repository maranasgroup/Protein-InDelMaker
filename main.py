print('-'*200)
print('***** Protein-InDelMaker : Developed at Prof. Costas Maranas lab at The Pennsylvania State University by Veda Sheersh Boorla. This package lets you make specific insertions,\
deletions,mutations or any combinations of these in any order.\
The package makes use of PyRosetta - version 4 *****')
print('-'*200)
print('***** Always examine the output structures to see if they make physical sense. Known bugs and fixes are given in README file.  *****')
print('***** Package author: Veda Sheersh Boorla\
. Any comments and bugs can be reported to mailforveda@gmail.com *****')
print('-'*200)

import os
import argparse
from core import indel

def interactive_main():
	print('-'*100)
	run_name = input('Enter a name for storing your results (Results are stored in ./results/name/)\nName can contain \
numbers,letters,underscore only\n:')
	print('-'*100)
	results_path = './results/{}'.format(run_name)
	os.system('mkdir ./results/{}'.format(run_name))
	os.system('mkdir ./results/{}/temp'.format(run_name))
	print('-'*100)
	pdb_path = input('Enter name of .pdb file\nThe name should contain fullpath\
to the file\n:')
	print('-'*100)
	input_path = input('Enter name of .input file\nThe name should contain \
fullpath to the file\n:')
	print('-'*100)
	ligand = input('Do you have a ligand in place? Enter yes or no:\n:')
	print('-'*100)
	if ligand=='yes':
		ligand=True
		ligand_params = input('Enter ligand params file:?\n:')
		print('-'*100)
		partners = input('Enter ligand partners:\n:')
		print('-'*100)
	else:
		ligand = False
		ligand_params=[]
		partners = ''
		
	close_cycles = input('Enter number of trials for loop closure\n\
Default is 20\n:')
	print('-'*100)
	refine_cycles = input('Enter number of refinement cycles after loop closure\n\
Default is 1\n:')
	print('-'*100)
	relax_cycles = input('Enter number of relax cycles at the end of each simulation\n\
Default is 3\n:')
	print('-'*100)
	debug = input('Run in debug mode?\n\
Default is yes (enter yes or no)\n:')
	if debug.lower()=='yes' or debug.lower()=='y':
		debug = True
	elif debug.lower()=='no' or debug.lower()=='n':
		debug = False
	print('-'*100)


	indel_mover = indel.InDelMut(temp_path='{}/temp'.format(results_path),
					results_path = results_path, pdb_path=pdb_path,
					input_path=input_path,close_cycles = int(close_cycles), 
					refine_cycles = int(refine_cycles),relax_cycles=int(relax_cycles),
					DEBUG=debug,watch=True,ligand=ligand,ligand_params=[ligand_params],partners=partners)
	indel_mover.run()
	return indel_mover

def main(args):
	
	run_name = args.run_name
	pdb_path = args.pdb_path
	input_path = args.input_path
	ligand = args.score_with_ligand
	close_cycles = args.n_loop_closure_cycles
	refine_cycles = args.n_loop_closure_refine_cycles
	relax_cycles = args.n_relax_cycles
	debug = args.debug_mode
	
	results_path = os.path.abspath('./results/{}'.format(run_name))
	if os.path.exists(results_path):
		print(f'{results_path} already exists..this will overwrite previous any previous files, quit now if you donot want that')
	else:
		os.mkdir(results_path)
	try:
		os.mkdir(os.path.abspath('./results/{}/temp'.format(run_name)))
	except:
		pass
	
	pdb_path = os.path.abspath('{}'.format(pdb_path))
	input_path = os.path.abspath('{}'.format(input_path))
	if ligand:
		ligand_params = args.ligand_params_path
		partners = args.ligand_partners
	else:
		ligand = False
		ligand_params=[]
		partners = ''

	indel_mover = indel.InDelMut(temp_path='{}/temp'.format(results_path),
					results_path = results_path, pdb_path=pdb_path,
					input_path=input_path,close_cycles = int(close_cycles), 
					refine_cycles = int(refine_cycles),relax_cycles=int(relax_cycles),
					DEBUG=debug,watch=True,ligand=ligand,ligand_params=[ligand_params],partners=partners)
	indel_mover.run()
	return indel_mover

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Protein-InDelMaker")

	parser.add_argument("run_name", help="'Enter a {name} for storing your results, {name} cannot have spaces (Results are stored in directory ./results/{name})", type=str)
	parser.add_argument("pdb_path", help="'Enter path to input .pdb file", type=str)
	parser.add_argument("input_path", help="'Enter path to input .input file", type=str)
	parser.add_argument("--score_with_ligand", help="Is there a ligand present in the pdb file for scoring?", action="store_true", default=False)
	parser.add_argument("--ligand_params_path", help="params file in Rosetta format for ligand", type=str)
	parser.add_argument("--ligand_partners", help="partners string in Rosetta format for scoring (eg., A_B)", type=str)
	parser.add_argument("--n_loop_closure_cycles", help="number of cycles for loop closure", type=int,default=20)
	parser.add_argument("--n_loop_closure_refine_cycles", help="number of cycles for refinment after loop closure", type=int,default=1)
	parser.add_argument("--n_relax_cycles", help="number of cycles for relax", type=int,default=3)
	parser.add_argument("--debug_mode", help="Run in debug mode?", action="store_true", default=False)

	args = parser.parse_args()
	
	main(args)
