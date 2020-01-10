import copy

import sys

import math

import numpy as np

class OrderedDict(dict):
    """
    Description:
        A wrapper object around dict which stores additional info. of the order
        in which (key,value) pairs are added to it.
    """
    def __init__(self):
        self.data = {}
        self.order = []

    def __setitem__(self, key, item):
        self.data[key] = item
        self.order.append(key)

    def __delitem__(self, key):
        del self.data[key]
        self.order.remove(key)

    def __iter__(self):
        return self.order.__iter__()

    def next(self):
        """
        Note:-
            This 'self' is the return value of __iter__
        """
        for key in self:
            return self.data[key]

    def __repr__(self): return repr(self.data)
    __hash__ = None # Avoid Py3k warning
    def __len__(self): return len(self.data)
    def __getitem__(self, key):
        if key in self.data:
            return self.data[key]
        if hasattr(self.__class__, "__missing__"):
            return self.__class__.__missing__(self, key)
        raise KeyError(key)

    def insert_after_index(self,index,key,value):
        print('Pushed key {} after {}'.format(key,index+1))
        self.order.insert(index+1,key)
        self.data[key] = value

    def get_by_index(self,index): return self.data[self.order[index]]
    def keys(self): return [key for key in self.order]
    def items(self): return [(key,self.data[key]) for key in self.order]
    def values(self): return [self.data[key] for key in self.order]
    def has_key(self, key): return key in self.order
    def __contains__(self, key): return key in self.order


class vector3d :
    def __init__( self, x=0, y=0, z=0 ) :
        self.x_ = round(x,3)
        self.y_ = round(y,3)
        self.z_ = round(z,3)
    def __add__( self, other ) :
        v = vector3d()
        v.x_ = self.x_ + other.x_
        v.y_ = self.y_ + other.y_
        v.z_ = self.z_ + other.z_
        return v
    def __sub__( self, other ) :
        v = vector3d()
        v.x_ = self.x_ - other.x_
        v.y_ = self.y_ - other.y_
        v.z_ = self.z_ - other.z_
        return v
    def __mul__( self, other ) :
        dotprod = 0
        dotprod += self.x_ * other.x_
        dotprod += self.y_ * other.y_
        dotprod += self.z_ * other.z_
        return dotprod
    def __str__( self ) :
        me = "("
        me += str( self.x_ ) + ", "
        me += str( self.y_ ) + ", "
        me += str( self.z_ ) + ")"
        return me

    def format_for_pdb( self ):
        rounded = (round(self.x,3),round(self.y,3),round(self.z,3))
        str_rounded = map(str,rounded)
        f_str = '{:>8}'
        return f_str.format(str_rounded[0])+\
                f_str.format(str_rounded[1])+\
                    f_str.format(str_rounded[2])

    def xyz(self):
        return (self.x(), self.y(), self.z())

    def x(self) :
        return self.x_
    def y(self) :
        return self.y_
    def z(self) :
        return self.z_
    def copy(self, other ) :
        self.x_ = other.x_;
        self.y_ = other.y_;
        self.z_ = other.z_;
    def set( self, x, y, z ):
        self.x_ = x
        self.y_ = y
        self.z_ = z
    def set_x( self, val ) :
        self.x_ = val
    def set_y( self, val ) :
        self.y_ = val
    def set_z( self, val ) :
        self.z_ = val
    def distance_squared( self, other ) :
        return ( self.x_ - other.x_ ) * (self.x_ - other.x_ ) + \
               ( self.y_ - other.y_ ) * (self.y_ - other.y_ ) + \
               ( self.z_ - other.z_ ) * (self.z_ - other.z_ )
    def distance( self, other ) :
        return math.sqrt( self.distance_squared(other) )
    def normalize( self ) :
        len = self.length()
        if len == 0.0 :
            return
        invlength = 1 / len
        self.scale( invlength )
        return self
    def scale( self, stretch ) :
        self.x_ *= stretch
        self.y_ *= stretch
        self.z_ *= stretch
        return self
    def length( self ) :
        return math.sqrt( self * self )
    def max( self, other ) :
        if self.x_ < other.x_ : self.x_ = other.x_
        if self.y_ < other.y_ : self.y_ = other.y_
        if self.z_ < other.z_ : self.z_ = other.z_
    def min( self, other ) :
        if self.x_ > other.x_ : self.x_ = other.x_
        if self.y_ > other.y_ : self.y_ = other.y_
        if self.z_ > other.z_ : self.z_ = other.z_

    def cross(self,other):
        return vector3d(*(np.cross(np.array(self.x(),self.y(),self.z())\
            ,np.array(other.x(),other.y(),other.z()))))

    def angle(self,other):
        return 180*np.arccos(\
            self.__mul__(other)/self.length()/other.length()\
            )/np.pi

def is_integer( s ) :
    try :
        int(s)
        return True
    except :
        return False

class Atom :
    def __init__( self ) :
        self.xyz = vector3d()
        self.name = ""
        self.pdb_index = 0
        self.occupancy = 1.0
        self.bfactor = 0
        self.het = False
        self.elem = ""

    def translate( self, xt, yt, zt):
        x,y,z = self.xyz.xyz()
        self.xyz = vector3d(x-xt,y-yt,z-zt)

    def rotate(self, angle_deg, axis):
        assert(axis.lower() in ['x','y','z'])

        r = Rotator.from_euler(axis.lower(), angle_deg, degrees=True)
        vec = self.xyz.xyz()
        self.xyz = vector3d(*r.apply(vec))

class Residue :
    def __init__( self ) :
        self.atoms = OrderedDict()
        self.stripped_atoms = OrderedDict()
        self.stripped_atmap = {}
        self.resstring = "" # residue index + insertion code -- a string
        self.resname = ""
        self.insertion_code = ""
        self.chain = None #pointer to the containing chain

    def is_Nterm(self):
        # Am I the first residue in Chain?

        # there got to be a better way to do this?
        ind = self.chain.residues.values().index(self)
        if ind==0:
            return True
        else:
            return False

    def is_Cterm(self):
        # Am I the first residue in Chain?

        ind = self.chain.residues.values().index(self)
        if ind==len(self.chain.residues):
            return True
        else:
            return False

    def add_atom( self, atom ) :
        self.atoms[atom.name] = atom
        self.stripped_atoms[ atom.name.strip() ] = atom
        self.stripped_atmap[ atom.name.strip() ] = atom

    def atom( self, atname ) :
        if atname in self.atmap :
            return self.atmap[ atname ]
        elif atname.strip() in self.stripped_atmap :
            return self.stripped_atmap[ atname.strip() ]
        else:
            print("Error in looking up atom",atname,"in residue", self.resname)
            sys.exit(1)

    def get_atom( self, atname ) :
        if atname in self.atoms :
            return self.atoms[ atname ]
        elif atname.strip() in self.stripped_atoms :
            return self.stripped_atoms[ atname.strip() ]
        else:
            print("Error in looking up atom",atname,"in residue", self.resname)
            sys.exit(1)

    def has_atom( self, atname ) :
        return (atname in self.atmap) or (atname.strip() in self.stripped_atmap)

    def claim( self, containing_chain ) :
        self.chain = containing_chain

    def resid( self ) :
        assert( self.chain )
        return self.chain.chain_name + " " + self.resstring

    def get_atom_vectors(self, atoms='ALL'):
        vectors = []
        if atoms=='ALL':
            atoms = self.atoms.keys()
        for a in atoms:
            atom = self.get_atom(str(a))
            vec = atom.xyz
            vectors.append(np.array([vec.x(),vec.y(),vec.z()]))
        return vectors

    def translate(self, xt, yt, zt, atoms='ALL'):
        if atoms=='ALL':
            atoms = self.atoms.keys()
        for a in atoms:
            atom = self.get_atom(str(a))
            atom.translate(xt,yt,zt)

    def rotate(self, angle_deg, axis, atoms='ALL'):
        assert(axis.lower() in ['x','y','z'])

        if atoms=='ALL':
            atoms = self.atoms.keys()
        for at in atoms:
            atom = self.get_atom(at)
            atom.rotate(angle_deg,axis.lower())

class Chain:
    def __init__( self ) :
        self.chain_name = ""
        self.residues = OrderedDict()

    def __len__(self): return len(self.residues)

    def add_residue( self, residue ):
        self.residues[residue.resstring] = residue
        residue.claim( self )

    def insert_residue(self, indexpos, resstring, residue):
        try:
            assert(isinstance(residue,Residue))
        except AssertionError:
            print('ERROR:Chain:insert_residue: Cannot insert object of \
                {} type into Chain'.format(type(residue)))
            raise(AssertionError)
            return None

        try:
            index = self.residues.order.index(str(indexpos))
        except ValueError:
            print('ERROR:Chain.insert_residue: Failed to insert at position {}\
                 - position does not exist!'.format(pos))
            return None

        print('INFO:Chain.insert_residue: Inserting residue {} \
            after position {}'.format(resstring,index))
        self.residues.insert_after_index(index, resstring, residue)

    def replace_residue( self, newresidue ):
        # in place replacement; keep the original location of the
        # residue in the self.residues array
        copyres = copy.copy( newresidue )
        copyres.claim( self )
        if newresidue.resstring not in self.resmap :
            print("Could not replace residue", newresidue.resstring)
            print(len( self.resmap ))
        assert( newresidue.resstring in self.resmap )
        for i in xrange(len( self.residues )) :
            if self.residues[ i ].resstring == newresidue.resstring :
                self.residues[ i ] = copyres
        self.resmap[ newresidue.resstring ] = copyres

    def get_residue( self, resstring ) :
        return self.residues[ resstring ]

    def get_atom_vectors(self, residues='ALL'):
        vectors = []
        if residues=='ALL':
            residues = self.residues.keys()
        for res in residues:
            residue = self.get_residue(str(res))
            vectors+=residue.get_atom_vectors()
        return vectors

    def translate(self, xt, yt, zt, residues='ALL'):
        if residues=='ALL':
            residues = self.residues.keys()
        for r in residues:
            residue = self.get_residue(r)
            residue.translate(xt,yt,zt)

    def rotate(self, angle_deg, axis, residues='ALL'):
        assert(axis.lower() in ['x','y','z'])

        if residues=='ALL':
            residues = self.residues.keys()
        for res in residues:
            residue = self.get_residue(res)
            residue.rotate(angle_deg,axis.lower())

    def get_residues(self, residues='ALL'):
        if residues=='ALL':
            return self.residues

class PDBStructure :
    def __init__( self, pdbfile=None ):
        self.chains = OrderedDict()
        self.chain_name = ""
        if pdbfile:
            self.read_from_lines( open( pdbfile ).readlines() )

    def info( self ):
        '''
        Display basic info
        '''
        for c,chain in self.chains.items():
            resfirst = chain.residues.get_by_index(0)
            reslast = chain.residues.get_by_index(-1)

            print('Chain {} : {} to {} = {} residues'.\
                format(c, resfirst.resstring+'_'+resfirst.resname, \
                    reslast.resstring+'_'+reslast.resname, len(chain)))

    def insert_pdb(self, pdb, pos, c):
        '''
        Insert another PDBStructure instance after position 'pos' in chain 'c'
        '''
        try:
            assert(isinstance(pdb,PDBStructure))
        except AssertionError:
            raise(AssertionError)
            print('ERROR:PDBStructure:insert_pdb: Cannot insert object of \
                {} type into pdb'.format(type(pdb)))
            return None

        chain = self.get_chain(c)
        res = chain.get_residue(str(pos))

        respos2 = 1
        # for every chain of pdb
        for c2, chain2 in pdb.chains.items():
            # insert every residue
            for r2, res2 in chain2.residues.items():
                resstring2 = '{}I'.format(respos2)
                res2.resstring = resstring2
                chain.insert_residue(pos,resstring2,res2)
                respos2+=1
                pos = resstring2

    def delete_res(self, resnum, chain_name):
        chain = self.get_chain(chain_name)
        if str(resnum) in chain.residues:
            del(chain.residues[str(resnum)])
            print('INFO:PDBStructure:delete_res:Deleted residue {} from chain{}\
            '.format(resnum,chain_name))

    def get_chain(self, chain_name):
        try:
            chain = self.chains[chain_name]
            return chain
        except KeyError:
            print('ERROR:PDBStructure.get_chain(chain_name={}): \
            Chain does not exist!'.format(chain_name))
            return None

    def get_residue( self, chnm, resstring ) :
        return self.chains[ chnm ].get_residue( resstring )

    def get_residues(self, chains='ALL'):
        residues = OrderedDict()
        if chains=='ALL':
            chains = self.chains.keys()
        for c in chains:
            chain = self.get_chain(c)
            res = chain.get_residues()
            residues[c] = res
        return residues

    def get_atom_vectors(self, chains='ALL'):
        vectors = []
        if chains=='ALL':
            chains = self.chains.keys()
        for c in chains:
            chain = self.get_chain(c)
            vecs = chain.get_atom_vectors()
            vectors+=vecs

        return np.array(vectors)

    def add_chain( self, chain ) :
        self.chains[chain.chain_name] = chain

    def read_from_lines( self, lines ):
        last_chain = Chain()
        last_residue = Residue()
        for line in lines :
            if line[0:4] == "ATOM" or line[0:6] == "HETATM" :
                chnm = self.chain_name_from_pdbline( line )
                if last_chain.chain_name != "" and chnm!=last_chain.chain_name :
                    last_chain.add_residue( last_residue )
                    last_residue = Residue()
                    if last_chain.chain_name not in self.chains :
                        self.add_chain( last_chain )
                    if chnm not in self.chains :
                        last_chain = Chain()
                    else :
                        # this chain has already been added,
                        # but now we have more residues
                        last_chain = self.chains[ chnm ]
                if last_chain.chain_name == "":
                    last_chain.chain_name = chnm
                resstring = self.resstring_from_pdbline( line )
                if last_residue.resname!= "" and \
                                            last_residue.resstring!=resstring :
                    last_chain.add_residue( last_residue )
                    last_residue = Residue()
                if last_residue.resname == "" :
                    last_residue.resname = self.resname_from_pdbline( line )
                    last_residue.resstring = resstring
                    #print "Read residue", last_residue.resname,
                    # last_chain.chain_name, last_residue.resstring
                atom = Atom()
                atom.xyz = self.xyz_from_pdbline( line )
                atom.name = self.atname_from_pdbline( line )
                atom.pdb_index = self.atnum_from_pdbline( line )
                atom.occupancy = self.occupancy_from_pdbline( line )
                atom.bfactor = self.bfactor_from_pdbline( line )
                atom.element = self.element_from_pdbline( line )
                atom.het = line[0:6] == "HETATM"
                last_residue.add_atom( atom )
        if last_residue.resname != "":
            last_chain.add_residue( last_residue )
        if last_chain.chain_name != "" and \
                                    last_chain.chain_name not in self.chains :
            self.add_chain( last_chain )
    def pdb_atname_range( self ) :     return ( 12, 16 )
    def pdb_atnum_range( self ) :      return (  6, 11 )
    def pdb_resname_range( self ) :    return ( 17, 20 )
    def pdb_chain_name_range( self ) : return ( 21, 22 )
    def pdb_resstring_range( self ) :  return ( 22, 27 )
    def pdb_xcoord_range( self ) :     return ( 30, 38 )
    def pdb_ycoord_range( self ) :     return ( 38, 46 )
    def pdb_zcoord_range( self ) :     return ( 46, 54 )
    def pdb_occupancy_range( self ) :  return ( 56, 60 )
    def pdb_bfactor_range( self ) :    return ( 61, 66 )
    def element_range( self ) :        return ( 76, 78 )

    def atname_from_pdbline( self, line ) :
        return line[ self.pdb_atname_range()[0]:self.pdb_atname_range()[1] ]
    def atnum_from_pdbline( self, line ):
        return line[ self.pdb_atnum_range()[0]:self.pdb_atnum_range()[1] ]
    def resname_from_pdbline( self, line ):
        return line[ self.pdb_resname_range()[0]:self.pdb_resname_range()[1] ]
    def chain_name_from_pdbline( self, line ):
        return line[ self.pdb_chain_name_range()[0]:self.pdb_chain_name_range()[1] ]
    def resstring_from_pdbline( self, line ):
        return line[ self.pdb_resstring_range()[0]:self.pdb_resstring_range()[1] ].strip()
    def xyz_from_pdbline( self, line ):
        if len( line ) < 50 :
            return None
        xstr = line[ self.pdb_xcoord_range()[0]:self.pdb_xcoord_range()[1] ]
        ystr = line[ self.pdb_ycoord_range()[0]:self.pdb_ycoord_range()[1] ]
        zstr = line[ self.pdb_zcoord_range()[0]:self.pdb_zcoord_range()[1] ]
        return vector3d( float(xstr), float(ystr),float(zstr))
    def occupancy_from_pdbline( self, line ):
        return 1.0 if len(line) < self.pdb_occupancy_range()[1] else float( line[ self.pdb_occupancy_range()[0]:self.pdb_occupancy_range()[1] ] )
    def bfactor_from_pdbline( self, line ):
        return 0.0 if len(line) < self.pdb_bfactor_range()[1] else float( line[ self.pdb_bfactor_range()[0]:self.pdb_bfactor_range()[1] ] )
    def element_from_pdbline( self, line ) :
        therange = self.element_range()
        return "" if len(line) < therange[1] else line[ therange[0]:therange[1] ]
    def pdb_lines( self ) :
        lines_all = OrderedDict()
        count_atoms = 0
        for c,chain in self.chains.items() :
            lines = []
            for res in chain.residues.values() :
                for atom in res.atoms.values() :
                    line = " "*80
                    count_atoms += 1
                    if atom.het :
                        line = line[:0] + "HETATM" + line[6:]
                    else :
                        line = line[:0] + "ATOM" + line[4:]
                    line = line[:self.pdb_atname_range()[0]] + atom.name + line[self.pdb_atname_range()[1]:]
                    line = line[:self.pdb_atnum_range()[0]] + ( "%5d" % count_atoms ) + line[self.pdb_atnum_range()[1]:]
                    line = line[:self.pdb_resname_range()[0]] + res.resname  + line[self.pdb_resname_range()[1]:]
                    line = line[:self.pdb_chain_name_range()[0]] + chain.chain_name + line[self.pdb_chain_name_range()[1]:]
                    if is_integer( res.resstring ) :
                        line = line[:self.pdb_resstring_range()[0]] + ("%5s" % (res.resstring+" ") ) + line[self.pdb_resstring_range()[1]:]
                    else :
                        line = line[:self.pdb_resstring_range()[0]] + ("%5s" % (res.resstring    ) ) + line[self.pdb_resstring_range()[1]:]
                    line = line[:self.pdb_xcoord_range()[0]] + ( " %7.3f" % atom.xyz.x() ) + line[self.pdb_xcoord_range()[1]:]
                    line = line[:self.pdb_ycoord_range()[0]] + ( " %7.3f" % atom.xyz.y() ) + line[self.pdb_ycoord_range()[1]:]
                    line = line[:self.pdb_zcoord_range()[0]] + ( " %7.3f" % atom.xyz.z() ) + line[self.pdb_zcoord_range()[1]:]
                    line = line[:self.pdb_occupancy_range()[0]] + ( "%4.2f" % atom.occupancy ) + line[self.pdb_occupancy_range()[1]:]
                    line = line[:self.pdb_bfactor_range()[0]] + ( "%5.2f" % atom.bfactor ) + line[self.pdb_bfactor_range()[1]:]
                    line = line[:self.element_range()[0]] + ( "%2s" % atom.element ) + line[self.element_range()[1]:]

                    lines.append( line + "\n" )
            lines.append( "TER\n")
            lines_all[c] = lines
        return lines_all

    def translate(self, xt, yt, zt, chains='ALL'):
        if chains=='ALL':
            chains = self.chains.keys()
        for c in chains:
            chain = self.get_chain(c)
            chain.translate(xt,yt,zt)

    def rotate(self, angle_deg, axis, chains='ALL'):
        assert(axis.lower() in ['x','y','z'])

        if chains=='ALL':
            chains = self.chains.keys()
        for c in chains:
            chain = self.get_chain(c)
            chain.rotate(angle_deg,axis.lower())
