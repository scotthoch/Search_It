#Scott Hoch June 2013
'''
TODO
-require that a word is shorter than the line you are looking at before searching in wrap mode (no double letters)
-search double letters in word list using regex and find those locations in word search before brute force search
-learn cython and implement KMP search to replace find() method
-Start code with os.path.getsize(filepath) to see if the wordsearch is bigger than 25% of RAM.  If it is, read in the
entire word list and then load the wordsearch in chunks that overlap enough to find the longest unfound word in
the overlap
'''
#Word search code
#Required to take an input file containing the dimensions of the word search grid,
#the letters in the word search grid listed by row,
#a command 'WRAP', or 'NO_WRAP' indicating if words can be wrapped across the grid,
#the number of words to be found,
#and a list of the words to be found.
#The output should contain the start and end coordinates of all of the words that are found in the
#word search, or the text 'NOT FOUND' if a word is not found.

import re as re

############################
######## define functions ########
############################
def print_found_words(word_list, word_dict):
        print word_dict
        for word in word_list:
            if type(word_dict['%s'%word][0])==str:
                print word_dict['%s'%word]
            else:
                print word_dict['%s'%word][0], word_dict['%s'%word][1]

def join_list(List, delimiter, wrap):
    '''turn a 2 D list of strings into a single string'''
    multiplier=1+wrap
    return "%s"%delimiter.join(["".join(List[x]*multiplier) for x in range(len(List))])

def word_search(letter_string, word_to_find, direction, diag=0):
    '''letter_string -> a string that contains a collection of letters to search
       word_to_find -> a string to search for in letter_string
       direction -> 1 for forwards, -1 for backwards
       if diag =1, return the number of commas behind the found word
       if word is found returns  start, and stop position of word, else returns a, 0'''
    #look for word in forward direction.  SInce mostly single searches, this is faster than regex
    start=letter_string.find('%s' %word_to_find[::direction])
    #if  word was found get its start and stop position
    if start>=0:
        comma_count=letter_string[0:start].count(',')
        start=start-comma_count#get rid of commas
        #forward spelling case
        if direction ==1:
            stop=start+len(word_to_find)-1
        #backward spelling case
        else:
            stop=start#reverse order when spelled backwards!
            start=stop+len(word_to_find)-1
        if diag==0:
            return start, stop
        else:
            return start, stop, comma_count
    else:
        if diag==0:
            return 'a', 0
        else:
            return 'a', 0, 0

def num_to_row_col(index, dim_mat, transpose,wrap):
    wrapper=1+wrap
    if transpose:
        index_row = index%dim_mat[0]
        index_col= index/dim_mat[0]/wrapper
    else:
        index_row = index/dim_mat[1]/wrapper
        index_col = index%dim_mat[1]
    return index_row, index_col

def remove_letters(letter_array, start=[0,0], stop=[0,0]):
    '''letter_array ->  list of lists containing the puzzle
       start, stop -> starting and stopping columns to remove letters from
       row -> row to remove letters from
       returns letter_array with letters [start:stop] replaced by '~' '''
    #make sure start and stop are in the correct order for the range function
    drow=stop[0]-start[0]
    dcol=stop[1]-start[1]
    #for a horizontal word
    if (dcol!=0 and drow==0):
        a, b = min(start[1],stop[1]), max(start[1],stop[1])
        for i in range(a,b+1):
            letter_array[start[0]][i]='~'
    #for a vertical word
    elif(drow!=0 and dcol==0):
        a, b = min(start[0],stop[0]), max(start[0],stop[0])
        for i in range(a,b+1):
            letter_array[i][start[1]]='~'
    #for a diagonal word
    elif (drow>0 and dcol>0) or (drow<0 and dcol<0) :
        a, b = min(start[0],stop[0]), min(start[1],stop[1])
        for i in range(0, abs(dcol)+1):
            letter_array[a+i][b+i]='~'
    return letter_array

def subtract_main_diags(start, stop, diag, main_length, wrap):
    wrapper = wrap+1
    let_to_skip=diag*main_length*wrapper
    return start - let_to_skip , stop - let_to_skip

def coordinate_to_grid(start, stop, diag, mis_match, match_side=1):
    if match_side*mis_match >=0: #more columns than rows or a square matrix
        return [start,start+diag],[stop, stop+diag]
    else:
        return [start+diag,start],[stop+diag,stop]

def fix_transpose_diag(start,stop,num_col):
    '''takes the start and stop coordinates and returns
        the start and stop coordinates with the mirror operation undone'''
    start[1]=num_col-1-start[1]
    stop[1]=num_col-1-stop[1]
    return start, stop

def get_diags(letter_array, dim_mat):
    '''Returns a matrix where the rows are continuous, non overlapping
       diagonals with negative slope from letter_array.
       starts at the bottom left of the matrix and works its way around'''
    #check to see if the matrix is square
    mis_match=dim_mat[1]-dim_mat[0]
    main_diags=abs(mis_match)
    diags_to_find=sum(dim_mat)-1
    #get the main diagonal
    diag_mat=[]
    main_length=min(dim_mat)
    diag_mat.append([letter_array[n][n] for n in range(0,main_length)])
    diags_to_find-=1
    if mis_match>=0: #for a matrix that has more columns than rows, or square
        for i in range(1,main_diags+1): #The number of main (longest) diagonals will be the dimension mismatch plus 1
            diag_mat.append([letter_array[n][n+i] for n in range(0,main_length)])
            diags_to_find-=1
        for j in range(1, int(0.5*diags_to_find)):
            diag_mat.append([letter_array[n][n+main_diags+j] for n in range(0,main_length-j)])
            diag_mat.append([letter_array[n+j][n] for n in range(0,main_length-j)])

    else: #for a matrix with more rows than columns,
        for i in range(1,abs(main_diags)+1):
            diag_mat.append([letter_array[n+i][n] for n in range(0,main_length)])
            diags_to_find-=1
        for j in range(1, int(0.5*diags_to_find)):
            diag_mat.append([letter_array[n+main_diags+j][n] for n in range(0,main_length-j)])
            diag_mat.append([letter_array[n][n+j] for n in range(0,main_length-j)])

    return diag_mat, mis_match, main_length

def row_search(letter_array, word_list, dim_mat, wrap, transpose):
    ''' letter_array ->  list of lists containing the word search grid
        word_list -> dictionary form {words to find : location if already found, or NOT_FOUND}
        dim_mat -> [rows, colums] dimension of letter_array
        transpose=1 if you want to search for words in the columns instead of rows
        searches for words row and column wise spelled both forwards and backwards
        returns letter_array with used letters replaced by '~', and updated word_list'''
    #turn nested lists into single line for faster searching, but keep original in memory for eliminating used letters
    if transpose==False:
        letter_array_holder=join_list(letter_array,',', wrap)
    else: #transpose matrix to look in columns
        letter_array_holder=join_list(map(list, zip(*letter_array)),',',wrap)
    #Get all the words in the dictionary for a search
    for word, found in word_list.items():
        #only search for a word if it has not been found yet
        if found=='NOT_FOUND':

            #check forward spelling case
            start, stop = word_search(letter_array_holder, word, 1)
            if type(start)==str:  #if word was not found check to see if it is spelled backwards
                start, stop = word_search(letter_array_holder, word, -1)
            if type(start)==int: #if word was found forward or backwards
                start_row, start_col = num_to_row_col(start, dim_mat, transpose, wrap)
                stop_row, stop_col = num_to_row_col(stop, dim_mat, transpose, wrap)
                word_list[word]=((start_row,start_col),(stop_row, stop_col))
                letter_array=remove_letters(letter_array,[start_row, start_col],[stop_row, stop_col] )

    return letter_array, word_list

def diag_search(letter_array, word_list, dim_mat, wrap, slope ):
    '''send in an array and this searches for words along the diagonals
       with slope \pm 1 '''
    #transform the array such that the negative slope diagonals become rows
    if slope==-1:
        diag_mat, mis_match, main_length = get_diags([w[::-1] for w in letter_array], dim_mat)
    else:
        diag_mat, mis_match, main_length = get_diags(letter_array, dim_mat)
    diag_lookup = [len(row) for row in diag_mat]
    diag_mat=join_list(diag_mat,',', wrap)
    print diag_mat
    num_main_diag=abs(mis_match)+1
    wrapper = 1+wrap
    for word, found in word_list.items():
    #only search for a word if it has not been found yet
        if found=='NOT_FOUND':
            #check forward spelling case
            start, stop, diag = word_search(diag_mat, word, 1, 1)
            if type(start)==str: #check backward spelling if not found
                start, stop, diag = word_search(diag_mat, word, -1, 1)
            if type(start)==int: # if we have found the word in the puzzle
                if diag==0: #on main diagonal
                    start, stop = start%main_length, stop%main_length
                    start, stop = coordinate_to_grid(start, stop, diag, mis_match)
                elif diag<num_main_diag+1: # if on one of the main diagonals, or the next diagonal
                    start, stop = subtract_main_diags(start, stop, diag, main_length, wrap)
                    start, stop = start%diag_lookup[diag], stop%diag_lookup[diag]
                    start, stop = coordinate_to_grid(start, stop, diag, mis_match)
                else:
                    #subtract numbers as though all diagonals were the same length
                    start, stop =subtract_main_diags(start, stop, diag, main_length, wrap)
                    #add in the extra numbers positions that are missing.  two cases for even and odd value of diag
                    diag_past_main=diag-num_main_diag # How many incomplete diagonals do we need to correct for
                    if diag_past_main%2==0: #same side as mis_match
                        num_shorter=int(0.5*(2+(diag_past_main)))-1 #prepare for a sum of even numbers
                        start, stop = wrapper*num_shorter*(num_shorter+1)+start, wrapper*num_shorter*(num_shorter+1)+stop
                        offset=num_main_diag+(diag_past_main/2)
                        start, stop = start%diag_lookup[diag], stop%diag_lookup[diag]
                        start, stop = coordinate_to_grid(start, stop, offset, mis_match, 1)
                    if diag_past_main%2==1: #opposite side as mis_match
                        num_shorter=int(0.5*(2+(diag_past_main-1)))-1 # prepare for the sum of even numbers
                        start = wrapper*(num_shorter*(num_shorter+1)+(diag_past_main/2)+1)+start
                        stop = wrapper* (num_shorter*(num_shorter+1)+(diag_past_main/2)+1)+stop
                        offset=((diag_past_main+1)/2)
                        start, stop = start%diag_lookup[diag], stop%diag_lookup[diag]
                        start, stop = coordinate_to_grid(start, stop, offset,mis_match, -1)

                if slope==-1:
                    start, stop = fix_transpose_diag(start,stop,dim_mat[1])
                word_list[word] = ((start[0], start[1]),(stop[0], stop[1]))
                letter_array=remove_letters(letter_array,[start[0],start[1]],[stop[0],stop[1]])
    return letter_array, word_list

def check_all(letter_array, word_dict, dim_mat, wrap):
    '''send in an array of the letters that is the word search field, letter_array,
    a dictionary with keys being the words to find, and values default to not found, word_dict,
    along with the dimensions of the search field, dim_mat'''
    #first check rows
    letter_array, word_dict = row_search(letter_array, word_dict, dim_mat, wrap, transpose=False)
    letter_array, word_dict = row_search(letter_array, word_dict, dim_mat, wrap,  transpose=True)
    letter_array, word_dict = diag_search(letter_array, word_dict, dim_mat, wrap, slope=1)
    letter_array, word_dict = diag_search(letter_array, word_dict, dim_mat, wrap, slope= -1)
    return word_dict


############################
########## Main code ##########
############################

def main():
    global puzzle
    #get file path from user and read in all of the text
    fil='/Users/trogdor/Documents/personal/PR/code/Factual/wordInput.txt'#str(raw_input('Please type the path of the file to load ' ))
    f=open(fil,'r')
    raw_dat=f.read()

    #Convert string to an array dealing with new lines and blank spaces on every OS I know of
    raw_dat = filter(None,re.split('\r| |\n|\rn',raw_dat))

    #Create puzzle in a list of lists
    dim_mat=[int(raw_dat[0]),int(raw_dat[1])]
    puzzle = [list(raw_dat[x]) for x in range(2,dim_mat[0]+2)]
    #See if we are dealing with a wrapped function
    wrap=raw_dat[2+dim_mat[0]].lower()=='WRAP'.lower()
    #load the words into a dictionary with the default assumpion that they have not been found yet
    #I have to load them into a separate list so that I can print their locations in the correct order at the end of the code
    word_list=raw_dat[4+dim_mat[0]::]
    words_found=dict((word,'NOT_FOUND') for word in word_list)

    #complete the word search depending on the wrap criteria
    words_found=check_all(puzzle, words_found, dim_mat,wrap)
    print_found_words(word_list, words_found)

