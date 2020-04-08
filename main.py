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

if __name__ == '__main__':
	interactive_main()
