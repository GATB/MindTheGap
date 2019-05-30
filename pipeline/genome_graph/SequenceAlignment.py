''' ***********************************************
    
    Class(es) to perform pairwise sequence comparisons/alignments
    Author - Claire Lemaitre
    
    Usage:
    nw = NeedlemanWunsch(sequence1, sequence2, 10, -5, -5)
    id = nw.getIdentity()
    
    *********************************************** '''


class NeedlemanWunsch:
    ''' Pairwise Global Alignment with the Needleman and Wunsch algorithm
        stores a matrix |S|x|T| (|S|+1 lines and |T|+1columns), sequences S and T and the score system (match, mismatch, gap)
        defines Needleman and Wunsch global alignment functions 
        and computes identity ratio of the best alignment'''

    def __init__(self, S, T, match, mismatch, gap):
        ''' defines and stores initial values
            initiates the first row and first column for Global Alignment
            '''
        self.S = S
        self.T = T
        self.lS = len(S)
        self.lT = len(T)
        self.gap = gap
        self.match = match
        self.mismatch = mismatch
        self.matrix = [[0 for j in range(len(T)+1)] for i in range(len(S)+1)]
        for i in range(1,self.lS+1):
            self.matrix[i][0]=self.matrix[i-1][0]+self.gap
        for j in range(1,self.lT+1):
            self.matrix[0][j]=self.matrix[0][j-1]+self.gap

    def getIdentity(self):
        ''' main function '''
        self.fillMatrix()
        return self.backTrackAndIdentity()


    def score(self,letter1,letter2):
        if letter1 == letter2:
            return self.match
        else:
            return self.mismatch


    def fillMatrix(self):
        '''Fills the matrix'''
        for i in range(1,self.lS+1):
            for j in range(1,self.lT+1):
                self.matrix[i][j] = max(self.matrix[i-1][j] + self.gap,
                                        self.matrix[i][j-1] + self.gap,
                                        self.matrix[i-1][j-1] + self.score(self.S[i-1],self.T[j-1]) )

    def backTrackAndIdentity(self):
        ''' backtracks the filled matrix 
            and returns the Identity ratio of the best alignment
            id = nb_match / max(len(S),len(T))
            '''
        nb_match = 0
        # alignment_length = 0 # not used, usefull if one prefers id = nb_match/align_len
        i=self.lS
        j=self.lT
        while i>0 and j>0:
            if (self.matrix[i][j] == self.matrix[i-1][j-1] + self.score(self.S[i-1],self.T[j-1])):
                if self.S[i-1] == self.T[j-1]:
                    nb_match += 1
                i -= 1
                j -= 1
            elif self.matrix[i][j] == self.matrix[i-1][j] + self.gap:
                i -=1
            else:
                j -=1
            # alignment_length += 1
        '''
        usefull only if one prefers id = nb_match/align_len
        while i>0:
            i -=1
            alignment_length += 1
        while j>0:
            j -=1
            alignment_length += 1
            '''

        # print(f"nb_match={nb_match},l1={self.lS},l2={self.lT}")
        return nb_match/max(self.lS,self.lT)


