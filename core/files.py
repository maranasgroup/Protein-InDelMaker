from core.objects import *

class PDBFile(object):
    """The fundamental structure file - PDB format"""
    def __init__(self, fullpath):
        """
        Args:-
            fullpath: `str` of full path leading to the file
        Returns:-
            None
        """
        super(PDBFile, self).__init__()
        self.path = fullpath
        foos = fullpath.split('/')
        _, filename = foos[:-1],foos[-1]

        try:
            assert(self._test_extension('.pdb'))
            self.filename = filename
            self.pdb = self._get_pdbobj()

        except AssertionError:
            print('Provided file at {} does not have .pdb extension!'\
                  .format(fullpath))
            self.filename = filename
            self.pdb = None

    def _test_extension(self,ext='.pdb'):
        """Simple extension test"""
        if self.path.strip()[-4:]!=ext:
            return False
        else:
            return True

    def _get_pdbobj(self):
        return PDBStructure(self.path)

    def _is_res_list_present(self,res_list,chain):
        '''
        Check if all resnums of res_list are present in `chain` of pdb
        Args:
            'res_list': 'list', list of residue numbers
            'chain': 'str', single letter name of chain
        Returns:
            True: if all of res_list present
            False: otherwise
        '''

        c = self.pdb.get_chain(chain)
        count = 0

        if c:
            residues = map(int,c.residues.keys())
            for resnum in res_list:
                if resnum in residues:
                    count+=1
                else:
                    continue
            if count==len(res_list):
                return True
            else:
                return False
        else:
            return False

    def _is_res_range_present(self,start,end,chain):
        res_list = range(start,end+1)
        return self._is_res_list_present(res_list,chain)

    def dump_pdb(self, outname, chains = []):
        '''
        Dumps pdb in current state
        Args:
            'outname': 'str', /path/name for output
            'chains': 'list', list of chain names to output (Outputs all
                                 chains by default)
        Returns:
            None
        '''
        f = open(outname,'w')
        for c,lines in self.pdb.pdb_lines().items():
            if chains:
                if not c in chains:
                    continue
                else:
                    for line in lines:
                        f.write(line)
            else:
                for line in lines:
                    f.write(line)
        f.close()
        print('INFO:PDBFile:dump_pdb: Saved pdb as {}'.format(outname))

    def delete_res_list(self,reslist,chain):
        if self._is_res_list_present(reslist,chain):
            for res in reslist:
                self.pdb.delete_res(res,chain)
        else:
            print('ERROR:PDBFile.delete_res_list: Unable to delete! \
                One or more residues not found!')

    def delete_res_range(self,start,end,chain):
        if self._is_res_range_present(start,end,chain):
            for res in range(start,end+1):
                self.pdb.delete_res(res,chain)
        else:
            print('ERROR:PDBFile.delete_res_range: Unable to delete! \
                One or more residues not found!')


class InputFile(object):
    """docstring for InputFile"""
    def __init__(self, filepath):
        super(InputFile, self).__init__()
        self.path = filepath

        self._load()

    def _load(self):
        # Each line of input begins with either of
        HEADERS = ['INSERT','DELETE','MUTATE']

        f_handle = open(self.path)
        lines = []

        for line in f_handle:
            head, chain, tasks = line.strip().split(' ')
            if head.strip() in HEADERS:
                chain = chain.strip().split('_')[-1]
                tasks = tasks.strip().split(',')
                lines.append((head,chain,tasks))

        self.lines = lines
