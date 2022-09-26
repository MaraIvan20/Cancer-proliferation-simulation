

class Cell:

    def __init__(self,alive=False):
        self.alive=alive

    def __str__(self):
        if self.alive==False:
            return '.'
        else:
            return 'O'  


    def is_alive(self):
        return self.alive


    def update_cell(self,neigh_matrix):
        n_alive=0
        for i in range(3):
            for j in range (3):
                if neigh_matrix[i][j].is_alive()==True:
                    n_alive+=1
                    
        if self.is_alive()==True:
            if n_alive>=5 or n_alive<=2: #5 and 2 respectively, because we include the cell itself, which is alive
                self.alive=False
        else:
            if n_alive==3: 
                self.alive=True



                

class Cancer(Cell):

    def __str__(self):
        if self.alive==False:
            return '.'
        else:
            return 'X'


    def update_cell(self,neigh_matrix):
        n_alive=0
        for i in range(3):
            for j in range (3):
                if neigh_matrix[i][j].is_alive()==True:
                    n_alive+=1
                    
        if self.is_alive()==True:
            if n_alive>=6 or n_alive<=2: #6 and 2 respectively, because we include the cell itself, which is alive
                self.alive=False
        else:
            if n_alive==3:
                self.alive=True

    


class Tissue:

    def __init__(self,rows=1,cols=1,CellType=Cell):
        self.rows=rows
        self.cols=cols
        self.CellType=CellType  

        self.matrix=list()

        for i in range (self.rows):
            self.matrix.append([])
            for j in range (self.cols):
                self.matrix[i].append(self.CellType(False)) #initialise the cells as being not alive, although for cells of class Cell or Cancer this is the default
            #self.matrix will be a nested list, ie a list that has self.rows number of lists, each one with self.cols number of elements


        #initially, I used the code below to construct self.matrix  
        #the method below was far more efficient in terms of time, but I modified it due to errors
        #writing instance.matrix[2][0].alive=True in script.py would change the entire row 2 instead of a single element
                
        #for i in range (self.rows):
            #self.matrix.append([self.CellType(False)]*self.cols)
            

    def __str__(self):
        
        matrix_str=''

        for i in range (self.rows):
            for j in range (self.cols):
                matrix_str+=str(self.matrix[i][j])  #this will call the __str__() methods defined in the cell type classes
            matrix_str+='\n' #after every row we need a new line character

        return matrix_str


    def __getitem__(self,key):
        #print('getitem is called')
        return self.matrix[key]
    
    def __setitem__(self,key,value):
        #print('setitem is called')
        self.matrix[key]=value  

        #obs: with __getitem__ and __setitem__ defined like this, when we do tissue[1]=[cellsim.Cell(True)]*9, __setitem__ is called, while when we do tissue[1][2]=cellsim.Cell(True), __getitem__ is called
        


    def seed_from_matrix(self,eq_matrix):
        
        self.rows=len(eq_matrix)
        if self.rows>0:
            self.cols=len(eq_matrix[0])  
        else:
            self.cols=0

        self.matrix=list()

        if self.rows>0 and self.cols>0:
            self.CellType=type(eq_matrix[0][0])
            for i in range (self.rows):
                self.matrix.append([])
                for j in range (self.cols):
                    self.matrix[i].append(eq_matrix[i][j])
                    
        elif self.rows>0:  #because we may have empty rows, case in which we want self.matrix to still be a 2D array
            self.CellType=None
            self.matrix=[[]]*self.rows
            
        else:  #self.matrix will be a simple empty list if we don't have rows in our original argument list eq_matrix
            self.CellType=None
            



    def seed_from_file(self,filename,CellType=Cell):
        
        fin=open(filename,'r')
        self.matrix=[]
        self.CellType=CellType
        line=fin.readline().strip()
        self.cols=len(line)
        self.rows=0
        while line:
            self.rows+=1
            self.matrix.append([])
            for ch in line:
                if ch=='.':
                    self.matrix[self.rows-1].append(CellType(False))#append the rows with objects of type CellType that are dead, ie have the attribute alive False
                    #for self.rows number of lists, we will have lists from index 0 to index (self.rows-1)
                else: 
                    self.matrix[self.rows-1].append(CellType(True))
            line=fin.readline().strip()
            
        fin.close()



    def seed_random(self,confluency,CellType=Cell):

        from random import choices

        self.CellType=CellType
        bool_list=[True,False]
        results=choices(bool_list,weights=[confluency, 1-confluency],k=self.rows*self.cols) #create a list with k boolean values, where True has a proportion of confluency out of 1, so False has a proportion of (1-confluency)
        for i in range(self.rows):
            for j in range (self.cols):
                self.matrix[i][j]=self.CellType(results[i*self.cols+j]) #append self.matrix with elements from the list with boolean values in order; the element on row i and column j in self.matrix will be the (i*self.cols+j)th element in results




    def next_state(self):
        from copy import deepcopy

        #create a matrix by putting margins of dead cell to the initial self.matrix, in order to account for the neighbours of the cells that are on the margins in tissue:
        #the extended_matrix will have a dimension of (self.rows+2)x(self.cols+2)
        extended_matrix=[[self.CellType(False)]*(self.cols+2)]  #the first row contains only dead cells
        for i in self.matrix:
            extended_matrix.append([self.CellType(False)]+deepcopy(i)+[self.CellType(False)]) #the middle rows begin and end with a dead cell and have rows from self.matrix in between
            #it is essential to work with a deepcopy of each row, because the rows in self.matrix will be changing as we update each cell in self.matrix, causing an erroneus result otherwise      
        extended_matrix.append([self.CellType(False)]*(self.cols+2)) #the last row contains only dead cells

        #update each element of self.matrix based on its neighbours:
        for i in range(self.rows):
            for j in range (self.cols):
                neigh_matrix=[[extended_matrix[i][x] for x in [j,j+1,j+2]],[extended_matrix[i+1][x] for x in [j,j+1,j+2]],[extended_matrix[i+2][x] for x in [j,j+1,j+2]]] #create the matrix of neighbours for each cell
                #due to the extension of the self.matrix in extended_matrix, self.matrix[i][j] will be extended_matrix[i+1][j+1]
                self.matrix[i][j].update_cell(neigh_matrix)  #update each cell
                     
