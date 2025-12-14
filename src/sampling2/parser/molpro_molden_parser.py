#! /usr/bin/env python2

# mndo molden format

import sys
import os
import re
import math
from element import element


class molpro_molden_parser():
    """ parser every format of log file """
    def __init__(self, config = {}):
        """
        init. variables
        """
        self.files = {'log': 'molden.dat',  'newlog': 'force.dat'}
        self.freq = {}
        self.dim = {'n_atom': 10, 'n_mode': 24}

        if config != {}:
            self.files['log'] = config['log']
            self.files['newlog'] = config['newlog']            
            self.dim['n_atom'] = config['n_atom']
            self.dim['n_mode'] = config['n_mode']
            # self.read_log()
        
        return
   
    def __read_molden_log_geom(self):
        """ read gaussian geom data """      
        filename = self.files['log']
        fp = open(filename, "r")

        line = "EMPTY"
        
        # pat = re.compile("\[FR-COORD\]\s+(\d+)")
        # pat0 = re.compile("\[NATOM\]\s+(\S+)")
        pat = r'[FR-COORD]'
        # pat0 = r'[NATOM]'

        pos = 0
        # atom_pos = 0

        while line != "":
            line = fp.readline()
            # m = pat.search(line)
            # m0 = pat0.search(line)            
            if pat in line:
                # string = m.group(1)
                # self.dim['n_atom'] = int(fp.readline().strip().split()[-1])
                print "The number of atoms are: ", self.dim['n_atom']
                # atom_pos = fp.tell()
                pos = fp.tell()

        # read geom
        if pos != 0:
            print "Find Geometry"
            fp.seek(pos)
        else:
            print "FAILED TO READ GAUSSIAN FREQ. CHECK LOG FILE"
            exit(1)

        n_atom = self.dim['n_atom']
        geom = [{} for i in xrange(n_atom)]
        # start to read geom
        for i in xrange(n_atom):
            coord = [0.0 for j in xrange(3)]
            line = fp.readline()
            record = line.split()
            atom_label = record[0]
            # atom_num = int(record[2])
            coord[0] = float(record[1])
            coord[1] = float(record[2])
            coord[2] = float(record[3])
            atom = {'atom_label': atom_label, 'coord': coord}
            geom[i] = atom
            # geom[i]['atom_number'] = atom_num
            
        # fp.seek(atom_pos)
        e = element()
        for i in xrange(n_atom):
        #     line = fp.readline()
        #     record = line.split()
        #     atom_number = int(record[2])
            e.atom_label = geom[i]['atom_label']
            geom[i]['atom_number'] = e.get_data()
            e.atom_label = ''
 
        # print(geom)

        self.freq['geom'] = geom
        self.freq['n_atom'] = n_atom
        
        # dump xyz file for check
        file_xyz = open("check.xyz", "w")
        n_atom = self.dim['n_atom']
        print >>file_xyz, "%10d\n" % n_atom
        b2a = 0.529177
        for atom in geom:
            print >>file_xyz, "%10s%15.8f%15.8f%15.8f" % (atom['atom_label'],
                                                          atom['coord'][0]*b2a,
                                                          atom['coord'][1]*b2a,
                                                          atom['coord'][2]*b2a)
        file_xyz.close()

        
        fp.close()
        
        return



    def __read_molden_log_freq(self):
        """
        read in one block and assign variables
        """
        filename = self.files['log'] 
        fp = open(filename, "r")        
        n_atom = self.dim['n_atom']
        self.freq['freq'] = []
        line = fp.readline()
        
        # pat_freq = re.compile("\[FREQ\]\s+(\d+)")
        pat_freq = r'[FREQ]'
        while line != "":
            line = fp.readline()            
            # m1 = pat_freq.search(line)                
            if pat_freq in line:
                 break
            
        n_mode = self.dim['n_mode']
        mass = [1.0 for x in range(n_mode)]
        self.freq['mass'] = mass

        # skip the 0.0 frequencies, there are 3n modes 
        for i in range(6): 
            fp.readline() 
        
        for i in xrange(n_mode):
            line = fp.readline()
            string = line.strip()
            self.freq['freq'].append(float(string))
            
        # pat_rmass = r'[RMASS]'
        # while line != "":
        #     line = fp.readline() 
        #     if pat_rmass in line:
        #         break 

        # for i in xrange(n_mode):
        #     line = fp.readline() 
        #     string = line.strip() 
        #     self.freq['mass'].append(float(string))

        return

 
    def __read_molden_log_vib(self):
        """ read in every normal mode """
        filename = self.files['log']
        n_atom = self.dim['n_atom']
        fp = open(filename, "r")
        # pat1 = re.compile("\[FR-NORM-COORD\]\s+(\d+)")
        pat1 = r'[FR-NORM-COORD]'
        line = "EMPTY"
        while line != "":
            line = fp.readline()
            # m = pat1.search(line)
            # if m is not None:
            if pat1 in line:
                print "Find the normal modes"
                break
 
        # read data now
        n_mode = self.dim['n_mode']     
        normal_mode = [[0.0 for i in xrange(n_atom*3)] for j in xrange(n_mode)]

        for i in range(6 * (1 + n_atom)):
            fp.readline()

        # Write down the final transfermation matrix S(n_mode, n_atom)        
        for i in xrange(n_mode):
            line = fp.readline()
            for j in xrange(n_atom):
                line = fp.readline()
                record = line.split()
                normal_mode[i][j*3+0] = float(record[0])
                normal_mode[i][j*3+1] = float(record[1])
                normal_mode[i][j*3+2] = float(record[2])
                
        self.freq['coor_vib'] = normal_mode
        self.freq['n_mode'] = n_mode
        fp.close()
        
        return


    def dump_molden_log_freq(self):
        """
        write down freq related data in a specific format.
        """
        filename = self.files['newlog']
        fp = open(filename, "w")
        n_atom = self.freq['n_atom']
        n_mode = self.freq['n_mode']
        geom = self.freq['geom']
        print >>fp, "ATOMIC GEOMETRYS(AU)"
        # geometry
        print >>fp, "%10d" % (n_atom)
        for i in xrange(n_atom):
            atom = geom[i]
            print >>fp, "%5s%10d" % (atom['atom_label'], atom['atom_number']),
            for j in xrange(3):
                print >>fp, "%15.8f" % atom['coord'][j],
            print >>fp, ""
        # normal modes
        print >>fp, "%10d%10d" % (n_mode, n_atom)
        freq = self.freq['freq']
        mass = self.freq['mass']
        coor_vib = self.freq['coor_vib']
        print >>fp, "FREQUENCY(cm**-1)"
        for i in xrange(n_mode):
            print >>fp, "%15.8f" % (freq[i])
#            if (i != 0 and i % 5 == 4) or (i+1 == n_mode):
#                print >>fp, ""
                
        print >>fp, "REDUCED MASS(AMU)"        
        for i in xrange(n_mode):
            print >>fp, "%15.8f" % (mass[i])
#            if (i != 0 and i % 5 == 4) or (i+1 == n_mode):
#                print >>fp, ""
        print >>fp, "NORMAL MODE"
        for i in xrange(n_mode):
            this_mode = coor_vib[i]
            for j in xrange(n_atom):
                print >>fp, "%15.8f%15.8f%15.8f" \
                % (this_mode[j*3+0], this_mode[j*3+1], this_mode[j*3+2])
        fp.close()
        
        return
    
    def read_log(self):
        """
        read molden log file
        """ 
        filename = self.files['log']
        file_out = open(filename, "r")
        # read molden geometry
 
        self.__read_molden_log_geom()
        self.__read_molden_log_freq()
        self.__read_molden_log_vib()
        self.dump_molden_log_freq()

        return
    
### main program
if __name__ == "__main__":
    mol = molden_parser()
    mol.read_log()
 


    

