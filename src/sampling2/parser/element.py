from os import path 

class element:
    def __init__(self, atom_label_list=[], atom_num_list=[], atom_label='', atom_num=0):
        self.atom_label_dic = {}
        self.atom_num_dic = {} 

        self.atom_label_list = atom_label_list
        self.atom_num_list = atom_num_list
        self.atom_label = atom_label
        self.atom_num = atom_num 

        self.read_elements()

        return

    def read_elements(self):
        eF = open(path.split(path.realpath(__file__))[0] + '/' + 'elements.dat')
        eF.readline() 
        eF.readline()
        eF.readline() 
        for line in eF:
            line = line.strip().split() 
            if line == []:
                break 
            atom_label = line[0]
            atom_num = int(line[3])
            atom_mass = float(line[-1])
            self.atom_label_dic.update({atom_label.upper(): atom_num})
            self.atom_num_dic.update({atom_num: atom_label})

        return 

    def get_data(self):
        if self.atom_label_list != []:
            return [self.atom_label_dic[x.upper()] for x in self.atom_label_list]
        
        elif self.atom_num_list != []:
            return [self.atom_num_dic[x] for x in self.atom_num_list]
        
        elif self.atom_num != 0:
            return self.atom_num_dic[self.atom_num]
    
        elif self.atom_label != '':
            return self.atom_label_dic[self.atom_label.upper()]

        else:
            print('no data is processed')
            exit()
        
