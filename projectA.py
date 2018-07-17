"""
Artificial Intelligence Part A 
Student Name and Number: 1. Huang Jiang Ng (789425)
                         2. Xiang Yao Ng (779576)

"""

from builtins import int
from multiprocessing.dummy import list
from queue import PriorityQueue
from copy import copy, deepcopy
from argparse import OPTIONAL
from ipaddress import _split_optional_netmask


# Move Phase

###########################################################################


def possible_moves(board, piece_pos):
    """
     possible_moves is used to calculate possible moves for given piece
     input: a 2D list board configuration and a list of coordinates for 
             a piece
     output: a list of list of possible moves position
    """

    row = piece_pos[ROW]
    col = piece_pos[COL]
    possible_moves_position = []
    
    # calculate east direction
    if (col+1 < MAX_COL):
        if (board[row][col+1] == SPACE):
            possible_moves_position.append([row, col+1])
        
        else:
            if (col+2 < MAX_COL and board[row][col+2] == SPACE):
                possible_moves_position.append([row, col+2])
    
    # calculate west direction
    if (col-1 >= 0):        
        if (board[row][col-1] == SPACE):
            possible_moves_position.append([row, col-1])
        
        else:
            if (col-2 >= 0 and board[row][col-2] == SPACE):
                possible_moves_position.append([row, col-2])        
    
    # calculate south direction
    if (row+1 < MAX_ROW):
        if (board[row+1][col] == SPACE):
            possible_moves_position.append([row+1, col])
        
        else:
            if (row+2 < MAX_COL and board[row+2][col] == SPACE):
                possible_moves_position.append([row+2, col])    
    
    # calculate north direction
    if (row-1 >= 0):
        if (board[row-1][col] == SPACE):
            possible_moves_position.append([row-1, col])
        
        else:
            if (row-2 >= 0 and board[row-2][col] == SPACE):
                possible_moves_position.append([row-2, col])    
    
    return possible_moves_position

# massacre phase

#########################################################################


def is_black_eliminated(eliminate, board_pos):
    """
        is_black_eliminated checks whether the black piece in the 
        board in the function elimination has been eliminated before 
        checking whether the white piece should be eliminated
        input: eliminate, board_pos
        output: bool
    """
    
    for eliminated in eliminate:
        # check whether the black position on board has been eliminated
        if (eliminated == board_pos):
            return True
        
    return False


def elimination(board, piece_pos, prev_pos):
    """
        elimination defines whether any pieces on the board will be 
        eliminated if a piece is moved to its current position
        input: 2D list board configuration, 1D list of coordinates piece
             position
        output: a list of list of piece positions to be eliminated (eliminate)
        
    """
    
    eliminate = []
    new_board = move_piece(board, piece_pos, prev_pos)
    
    # check whether the current white piece is surrounding black color piece
    # then check whether another white color piece is surrounding the
    # black color piece
    
    # check east
    
    if (piece_pos[COL] + 1 < MAX_COL and 
            new_board[piece_pos[ROW]][piece_pos[COL] + 1] == BLACK):
        if (piece_pos[COL] + 2 < MAX_COL and 
                new_board[piece_pos[ROW]][piece_pos[COL] + 2] in (WHITE, CORNER)):
            eliminate.append([piece_pos[ROW], piece_pos[COL] + 1])
            
    # check west
    if (piece_pos[COL] - 1 >= MIN_COL and 
            new_board[piece_pos[ROW]][piece_pos[COL] - 1] == BLACK):
        if (piece_pos[COL] - 2 > MIN_COL and 
                new_board[piece_pos[ROW]][piece_pos[COL] - 2] in (WHITE, CORNER)):
            eliminate.append([piece_pos[ROW], piece_pos[COL] - 1])    

    # check south
    if (piece_pos[ROW] + 1 < MAX_ROW and 
            new_board[piece_pos[ROW] + 1][piece_pos[COL]] == BLACK):
        if (piece_pos[ROW] + 2 < MAX_ROW and 
                new_board[piece_pos[ROW] + 2][piece_pos[COL]] in (WHITE, CORNER)):
            eliminate.append([piece_pos[ROW] + 1, piece_pos[COL]])       
            
    # check north
    if (piece_pos[ROW] - 1 >= MIN_ROW and 
            new_board[piece_pos[ROW] - 1][piece_pos[COL]] == BLACK):
        if (piece_pos[ROW] - 2 >= MIN_ROW and 
                new_board[piece_pos[ROW] - 2][piece_pos[COL]] in (WHITE, CORNER)):
            eliminate.append([piece_pos[ROW] - 1, piece_pos[COL]])    
  
    # now we check whether the current white piece is instead being
    # surrounded by black pieces
    
    # check east-west
    if ((piece_pos[COL] + 1 < MAX_COL) and 
        (new_board[piece_pos[ROW]][piece_pos[COL] + 1] in (BLACK, CORNER))
        and (piece_pos[COL] - 1 >= MIN_COL) and 
            (new_board[piece_pos[ROW]][piece_pos[COL] - 1] in (BLACK, CORNER))):
        east_pos = [piece_pos[ROW], piece_pos[COL] + 1]
        west_pos = [piece_pos[ROW], piece_pos[COL] - 1]
        if ((is_black_eliminated(eliminate, east_pos) is False) and
                (is_black_eliminated(eliminate, west_pos) is False)):                                   
            eliminate.append([piece_pos[ROW], piece_pos[COL]])
            
    # check north-south
    if ((piece_pos[ROW] + 1 < MAX_ROW) and 
        (new_board[piece_pos[ROW] + 1][piece_pos[COL]] in (BLACK, CORNER))
        and (piece_pos[ROW] - 1 >= MIN_ROW) and 
            (new_board[piece_pos[ROW] - 1][piece_pos[COL]] in
                (BLACK, CORNER))):
        south_pos = [piece_pos[ROW] + 1, piece_pos[COL]]
        north_pos = [piece_pos[ROW] - 1, piece_pos[COL]]
        if ((is_black_eliminated(eliminate, south_pos) is False) and
                (is_black_eliminated(eliminate, north_pos) is False)):   
            eliminate.append([piece_pos[ROW], piece_pos[COL]])       
             
    return eliminate

def manhattan_dist(white_coord, black_coord):
    '''
        Find the manhattan distance between 2 coordinates 
        input: 2 lists of coordinates
        output: int
    '''
    row_diff = abs(white_coord[ROW] - black_coord[ROW])
    col_diff = abs(white_coord[COL] - black_coord[COL])
    return row_diff + col_diff


def find_closest_whites(black_id, white_dict, black_dict):
    '''
        i) Find the closest white pieces to this black piece according 
        to their manhattan distances
        ii) Find the total manhattan distances between the closest 
        white pieces and this black piece
        input: A list of coordinates of a black piece, and the 2 
                dictionaries of pieces
        output: void
    '''
    if (black_dict[black_id][IS_ALIVE] == DEAD):
        return None
    
    black_coord = black_dict[black_id][COORD]
    # Contain the ids of the closest and second closest whites
    whites_assigned = [None, None] 
    # Initialise the search for the 2 whites pieces
    closest_dist = MAX_DIST
    second_closest_dict = MAX_DIST
    
    for (key, white_stats) in white_dict.items():
        if (white_stats[IS_ALIVE] == DEAD):
            continue 
        else:
            dist = manhattan_dist(white_stats[COORD], black_coord)
            if (dist < closest_dist):
                second_closest_dict = closest_dist
                whites_assigned[SECOND] = whites_assigned[CLOSEST_WHITE]
                closest_dist = dist
                whites_assigned[CLOSEST_WHITE] = key
            elif ((closest_dist < dist and dist < second_closest_dict) 
                  or closest_dist == dist):
                second_closest_dict = dist
                whites_assigned[SECOND] = key
                
    black_dict[black_id][WHITES_ASSIGNED] = whites_assigned
    
    # Calculate the TOTAL manhattan distances if 2 closest whites are found        
    if (bool(whites_assigned[CLOSEST_WHITE]) and bool(whites_assigned[SECOND])):
        black_dict[black_id][TOTAL_DIST] = closest_dist + second_closest_dict
    # Take the single manhanttan distance if there's only 1 white piece on
    # board
    elif (bool(whites_assigned[CLOSEST_WHITE])):
        black_dict[black_id][TOTAL_DIST] = closest_dist
    # No white pieces on board at all
    else:
        black_dict[black_id][TOTAL_DIST] = None
    
    
def next_kill(board, black_dict):
    '''
        Check which black piece to be eliminated next based on the 
        location of closest white pieces. The lower the total manhattan 
        distance from whites to the black, the higher the priority for 
        the black.
        input: black dictionary
        output: black id to be killed next 
    '''
    # Sort total distance
    # Only include black pieces that are alive
    dist_dict = {}
    for (key, black_stats) in black_dict.items():
        if (black_stats[IS_ALIVE]):
            dist_dict[key] = black_stats[TOTAL_DIST]
    
    # End Massacre() if all black pieces have been eliminated
    if (bool(dist_dict) is False):
        return None
    
    # Check if the black piece of the highest priority can be eliminated
    # normally or by corner
    for key in sorted(dist_dict, key=dist_dict.get):
        elimination_method = can_eliminate(board, key, black_dict)
        if (elimination_method == NORMAL_ELIMNATION):
            return key
        elif (elimination_method == CORNER_ELIMINATION):
            # Delete the second closest white piece attached to this black
            black_dict[key][WHITES_ASSIGNED][SECOND] = None
            return key
            
    # End Massacre() when the available black pieces on board can't be
    # eliminated by any way
    return None


# Need to check if the black piece is still alive before calling this function
def can_eliminate(board, black_id, black_dict):
    '''
        Check whether adjacent squares of this black piece on the same
        row or column are empty, and whether this black piece has 2 
        white pieces assigned to eliminate it. Also check if this
        black piece can be eliminated by any corner. 
        input: 2-d array board, the black piece's id, black dictionary
        output: Black piece can be eliminated normally or by any corner, 
                or can't be eliminated yet
                (Condition represented by integers)
    '''
    black_coord = black_dict[black_id][COORD]
    # Calculate the number of white pieces assigned to this black piece
    whites_assigned = black_dict[black_id][WHITES_ASSIGNED]
    num_of_whites = FULL
    for white in whites_assigned:
        if (bool(white) is False):
            num_of_whites -= 1
    
    # check east west direction
    if ((black_coord[COL] + 1 < MAX_COL) and (black_coord[COL] - 1 >= 0)):
        if ((board[black_coord[ROW]][black_coord[COL] + 1] in 
             (WHITE, SPACE)) and 
                (board[black_coord[ROW]][black_coord[COL] - 1] in 
                    (WHITE, SPACE)) and (num_of_whites == FULL)):
            return NORMAL_ELIMNATION
        elif (board[black_coord[ROW]][black_coord[COL] + 1] == CORNER or
              board[black_coord[ROW]][black_coord[COL] - 1] == CORNER and
              num_of_whites >= 1):
            return CORNER_ELIMINATION

    # check north-south direction
    
    if ((black_coord[ROW] + 1 < MAX_ROW) and (black_coord[ROW] - 1 >= 0)):
        
        if ((board[black_coord[ROW] + 1][black_coord[COL]] in 
             (WHITE, SPACE)) and 
                (board[black_coord[ROW] - 1][black_coord[COL]] in 
                    (WHITE, SPACE)) and (num_of_whites == FULL)):
            return NORMAL_ELIMNATION
        elif (board[black_coord[ROW] + 1][black_coord[COL]] == CORNER or 
              board[black_coord[ROW] - 1][black_coord[COL]] == CORNER and 
              num_of_whites >= 1):
            return CORNER_ELIMINATION  
    return CANT_BE_ELIMINATED
   
                    
def total_hdist(black_id, path, black_dict):
    ''' 
        The total heuristic distance for A* is defined as the sum of 
        each white piece's manhattan distance to the closest adjacent
        space of the target black piece. There can be at most 2 white
        pieces assigned to eliminate the target black piece depending
        on the elimination method. 
        input: black_id, path, black_dict
        output: total heuristic distance to this black piece

    '''
    black_coord = black_dict[black_id][COORD]
    white1 = black_dict[black_id][WHITES_ASSIGNED][CLOSEST_WHITE]
    # white2 is None if this is a corner elimination
    white2 = black_dict[black_id][WHITES_ASSIGNED][SECOND]
    dist_array = []
    
    if (white2 is not None):
        # calculate distance of east west direction and store each
        # distance into array
        for i in (-1, 1):
            adj_coord = [black_coord[ROW], black_coord[COL] + i]
            dist1 = manhattan_dist(path[CLOSEST_WHITE][LAST], adj_coord)
            dist2 = manhattan_dist(path[SECOND][LAST], adj_coord)
            dist_array.append([dist1, dist2])
                
        # check north-south direction
        for i in (-1, 1):
            adj_coord = [black_coord[ROW] + i, black_coord[COL]]
            dist1 = manhattan_dist(path[CLOSEST_WHITE][LAST], adj_coord)
            dist2 = manhattan_dist(path[SECOND][LAST], adj_coord)
            dist_array.append([dist1, dist2])
            
        min_white1 = MAX_DIST
        min_white2 = MAX_DIST
        for [x, y] in dist_array:
            if (x < min_white1):
                min_white1 = x
        
            if (y < min_white2):
                min_white2 = y
    
        return min_white1 + min_white2
    
    else:
        # calculate distance of east west direction and store each
        # distance into array
        for i in (-1, 1):
            adj_coord = [black_coord[ROW], black_coord[COL] + i]
            dist1 = manhattan_dist(path[CLOSEST_WHITE][LAST], adj_coord)
            dist_array.append(dist1)
                
        # check north-south direction
        for i in (-1, 1):
            adj_coord = [black_coord[ROW] + i, black_coord[COL]]
            dist1 = manhattan_dist(path[CLOSEST_WHITE][LAST], adj_coord)
            dist_array.append(dist1)
            
        min_white = MAX_DIST
        for x in dist_array:
            if (x < min_white):
                min_white = x
        return min_white


def remove_piece(board, piece, black_dict, white_dict):
    """
        remove_piece eliminates pieces on the board when running Astar
        input: board, piece_id to be eliminated, black_dict, white_dict
        output: new_board
    """
    
    new_board = deepcopy(board)
    
    if ('W' in piece):
        piece_pos = white_dict[piece][COORD]
        new_board[piece_pos[ROW]][piece_pos[COL]] = SPACE
        
    elif ('B' in piece):
        piece_pos = black_dict[piece][COORD]
        new_board[piece_pos[ROW]][piece_pos[COL]] = SPACE
    
    return new_board   


def move_piece(board, curr_pos, prev_pos):
    """ 
        move_piece updates the current position of the piece
        on the board
        input: board, curr_pos, prev_pos
        output: new_board
    """
    new_board = deepcopy(board)
    new_board[prev_pos[ROW]][prev_pos[COL]] = SPACE
    new_board[curr_pos[ROW]][curr_pos[COL]] = WHITE
    
    return new_board


def recalculate_path(path, white_piece, curr_pos):
    """ 
        recalculate_path updates the current path of the white piece
        on the board
        input: path, which white piece, curr_pos of the white piece
        output: new_path
    """
    new_path = deepcopy(path)
    new_path[white_piece].append(curr_pos)
    
    return new_path


def retrieve_piece_type(board, coord):  
    return board[coord[ROW]][coord[COL]]


def debug_board(board):
    print("board positions: ")    
    for y in range(8):
        for z in range(8):
            print(board[y][z] + " ", end='')
        print() 

    
def find_blackid(black_pos, black_dict):
    """
        find_blackid searches for black pieces in the 
        dictionary
    """
    
    for (key, stats) in black_dict.items():
        if (stats[COORD] == black_pos):
            return key
    
    return None
    
    
def Astar(board, white_dict, black_dict, black_id):   
    """
        Astar is the main program running the massacre phase
        It checks all possible moves for the white pieces
        to reach black piece and find the smallest possible
        sequence of moves
        
    """ 
    q = PriorityQueue()
    # first initialize the queue with current board position
    path = [[], []]
    white1 = black_dict[black_id][WHITES_ASSIGNED][CLOSEST_WHITE]
    white1_pos = white_dict[white1][COORD]
    path[CLOSEST_WHITE].append(white1_pos)
    
    # check whether only 1 white piece is needed for elimination
    is_corner_elimination = \
        not bool(black_dict[black_id][WHITES_ASSIGNED][SECOND])
    if (not is_corner_elimination):
        white2 = black_dict[black_id][WHITES_ASSIGNED][SECOND]
        white2_pos = white_dict[white2][COORD]
        path[SECOND].append(white2_pos)
        
    h = total_hdist(black_id, path, black_dict)
    f = h
    
    seq_path = []
    initial_node = [f, h, board, path, seq_path]
    q.put(initial_node)
    
    # for each move of the white pieces, perform A*    
    while not q.empty():
        curr_node = q.get()
        white1_moves = possible_moves(board, curr_node[PATH][CLOSEST_WHITE][LAST])
        seq_path = curr_node[SEQ_PATH]
        
        for move in white1_moves:
            new_seq = deepcopy(seq_path)
            # goal test
            eliminated = elimination(curr_node[BOARD], move, curr_node[PATH][CLOSEST_WHITE][LAST])
            if (bool(eliminated)):
                # case 1: eliminate target black piece (and other black pieces)
                # case 2: eliminate other black piece(s) instead of the target
                # case 3: white piece gets eliminated without eliminating the target (bad move)
                target_eliminated = False
                white_eliminated = False
                # Store other black pieces' coords that can be eliminated 
                optional_elimination = [] 
                
                for i in range(len(eliminated)):
                    # goal state has been reached
                    if (black_dict[black_id][COORD] == eliminated[i]):
                        new_board = move_piece(curr_node[BOARD], move, 
                                               curr_node[PATH][CLOSEST_WHITE][LAST])
                        new_board = remove_piece(new_board, black_id, 
                                                 black_dict, white_dict)
                        new_path = recalculate_path(curr_node[PATH], CLOSEST_WHITE, move)
                        new_seq.append([curr_node[PATH][CLOSEST_WHITE][LAST], move])
                        black_dict[black_id][IS_ALIVE] = False
                        target_eliminated = True
                        
                    elif (move == eliminated[i]):
                        # If this moving white is the only eliminated piece, 
                        # Astar() won't expand it 
                        white_eliminated = True
                    
                    elif (retrieve_piece_type(curr_node[BOARD], 
                                              eliminated[i]) == BLACK):
                        optional_elimination.append(eliminated[i])
                
                # Explore alternative goal state
                if (optional_elimination and not target_eliminated):
                    # Eliminate this coincidental black piece if the white 
                    # piece in action is not sacrificed
                    if (not white_eliminated):
                        new_board = move_piece(curr_node[BOARD], move,
                                               curr_node[PATH][CLOSEST_WHITE][LAST])
                        new_path = recalculate_path(curr_node[PATH], CLOSEST_WHITE, move)
                        new_seq.append([curr_node[PATH][CLOSEST_WHITE][LAST], move])

                        for opt in optional_elimination:
                            other_id = find_blackid(opt, black_dict)
                            black_dict[other_id][IS_ALIVE] = False
                            new_board = remove_piece(new_board, other_id,
                                                     black_dict, white_dict)
                            
                        # Eliminate the target black piece next time
                        return (new_board, new_path, new_seq)
                        
                if (target_eliminated):
                    if (white_eliminated):
                        white_dict[white1][IS_ALIVE] = False
                        new_board = remove_piece(new_board, white1,
                                                 black_dict, white_dict)
                    
                    # If optional_elimination is not empty, remove 
                    # the additional eliminated black pieces
                    for opt in optional_elimination:
                        other_id = find_blackid(opt, black_dict)
                        black_dict[other_id][IS_ALIVE] = False
                        new_board = remove_piece(new_board, other_id, 
                                                 black_dict, white_dict)
                        
                    return (new_board, new_path, new_seq)        
                        
            else:
                new_board = move_piece(curr_node[BOARD], move, 
                                       curr_node[PATH][CLOSEST_WHITE][LAST])
                new_path = recalculate_path(curr_node[PATH], CLOSEST_WHITE, move)
                depth = len(new_path[CLOSEST_WHITE]) - 1 
                new_h = total_hdist(black_id, new_path, black_dict)
                new_seq.append([curr_node[PATH][CLOSEST_WHITE][LAST], move])
                q.put([depth + new_h, new_h, new_board, new_path, new_seq])
        
        #Explore the second white piece assigned if it's a normal elimination
        #Structure is the same as above
        if (not is_corner_elimination):
            white2_moves = possible_moves(curr_node[BOARD], curr_node[PATH][SECOND][LAST])
            for move in white2_moves:  
                eliminated = elimination(curr_node[BOARD], move,
                                         curr_node[PATH][SECOND][LAST])
                new_seq = deepcopy(seq_path)
                
                if (bool(eliminated)):
                    target_eliminated = False
                    white_eliminated = False
                    optional_elimination = [] 
                    
                    for i in range(len(eliminated)):
                        if (black_dict[black_id][COORD] == eliminated[i]):
                            new_board = move_piece(curr_node[BOARD], move, 
                                                   curr_node[PATH][SECOND][LAST])
                            new_board = remove_piece(new_board, black_id, 
                                                     black_dict, white_dict)
                            new_path = recalculate_path(curr_node[PATH], SECOND, move)
                            
                            new_seq.append([curr_node[PATH][SECOND][LAST], move])
                            black_dict[black_id][IS_ALIVE] = False
                            target_eliminated = True
                            
                        elif (move == eliminated[i]):
                            white_eliminated = True
                    
                        elif (retrieve_piece_type(curr_node[BOARD],
                                                  eliminated[i]) == BLACK):
                            optional_elimination.append(eliminated[i])
                    
                    if (optional_elimination and not target_eliminated):
                        if (not white_eliminated):
                            new_board = move_piece(curr_node[BOARD], move, 
                                                   curr_node[PATH][SECOND][LAST])
                            new_path = recalculate_path(curr_node[PATH], SECOND, move)
                            new_seq.append([curr_node[PATH][SECOND][LAST], move])
                            
                            for opt in optional_elimination:
                                other_id = find_blackid(opt, black_dict)
                                black_dict[other_id][IS_ALIVE] = False
                                new_board = \
                                    remove_piece(new_board, other_id, 
                                                 black_dict, white_dict)
                                
                            return (new_board, new_path, new_seq)
                            
                    if (target_eliminated):
                        if (white_eliminated):
                            white_dict[white2][IS_ALIVE] = False
                            new_board = \
                                remove_piece(new_board, white2, 
                                             black_dict, white_dict)
                        
                        for opt in optional_elimination:
                            other_id = find_blackid(opt, black_dict)
                            black_dict[other_id][IS_ALIVE] = False
                            new_board = remove_piece(new_board, opt, 
                                                     black_dict, white_dict)
                            
                        return (new_board, new_path, new_seq)        
                else:
                    
                    new_board = move_piece(curr_node[BOARD], move, 
                                           curr_node[PATH][SECOND][LAST])
                    new_path = recalculate_path(curr_node[PATH], SECOND, move)
                    depth = len(new_path[SECOND]) - 1
                    new_h = total_hdist(black_id, new_path, black_dict)
                    new_seq.append([curr_node[PATH][SECOND][LAST], move])
                    q.put([depth + new_h, new_h, new_board, new_path,
                           new_seq])
     

def moves_phase(board, white_dict, black_dict):
    """
        moves_phase is used to run moves and calculate the possible moves
        by black and white pieces
    """
    num_moves = 0
    
    # implement move phase
    for key in white_dict:
        moves = len(possible_moves(board, white_dict[key][COORD]))
        num_moves += moves
    
    print(num_moves)
    
    num_moves = 0
    for key in black_dict:
        moves = len(possible_moves(board, black_dict[key][COORD]))
        num_moves += moves
        
    print(num_moves)


def print_seq(seq_path):
    """
        helper function for massacre_phase
    """
    
    new_seq = deepcopy(seq_path)
    # switch the order of the coordinates from (row, col) to (col, row)
    for seq in new_seq:
        for coord in seq:
            temp = coord[0]
            coord[0] = coord[1]
            coord[1] = temp
            
    for seq in new_seq:
        print(str(tuple(seq[PREV])) + " -> " + str(tuple(seq[NEXT])))


def massacre_phase(board, white_dict, black_dict):  
    """
        massacre_phase runs massacre and returns a printing of 
        sequence of paths
    """
    
    while (True):
        for (key, stats) in black_dict.items():
            if (stats[IS_ALIVE]):
                find_closest_whites(key, white_dict, black_dict)
        
        curr_black = next_kill(board, black_dict)
        
        #All black pieces are eliminated, or can't find ways to eliminate remaining pieces 
        if (curr_black is None):
            break
        
        (board, path, seq_path) = Astar(board, white_dict, black_dict, 
                                        curr_black)
        
        white1 = black_dict[curr_black][WHITES_ASSIGNED][CLOSEST_WHITE]
        white_dict[white1][COORD] = path[CLOSEST_WHITE][LAST]
        
        is_corner_elimination = \
            not bool(black_dict[curr_black][WHITES_ASSIGNED][SECOND])
        if (not is_corner_elimination):
            white2 = black_dict[curr_black][WHITES_ASSIGNED][SECOND]
            white_dict[white2][COORD] = path[SECOND][LAST]
            
        print_seq(seq_path)
        
    
def dict_init(piece_id, my_dict, row, col):
    """
        dict_init is used to initialize dictionaries containing information 
        on the pieces
    """
    
    my_dict[piece_id] = [None]*DICT_SIZE
    my_dict[piece_id][COORD] = [row, int(col/2)]
    my_dict[piece_id][IS_ALIVE] = STILL_ALIVE
  
  
def main(text):
    
    board = [["-" for x in range(8)] for y in range(8)]
    # white_dict stores information for white pieces
    # key: piece id (e.g. W1)
    # values: [[own coordinates], is_alive]
    white_dict = {}
    # black_dict stores information for black pieces
    # key: piece id (e.g. B2)
    # values: [[own coordinates], is_alive, [flanking white ids]]
    black_dict = {}
    
    white_id_iter = iter(range(1, MAX_PIECE))
    black_id_iter = iter(range(1, MAX_PIECE))
    # read first 8 lines of inputs for board configuration
    for row in range(0, MAX_ROW):
        
        for col in range(0, 15):
            if (text[row][col] == WHITE):
                board[row][int(col/2)] = WHITE
                piece_id = 'W' + str(next(white_id_iter))
                dict_init(piece_id, white_dict, row, col)
                
            elif (text[row][col] == BLACK):
                board[row][int(col/2)] = BLACK
                piece_id = 'B' + str(next(black_id_iter))
                dict_init(piece_id, black_dict, row, col)
                
            elif (text[row][col] == CORNER):
                board[row][int(col/2)] = CORNER
                
    # implement moves/massacre phase
    if (text[LAST] == "Moves"):
        moves_phase(board, white_dict, black_dict)
    elif (text[LAST] == "Massacre"):
        massacre_phase(board, white_dict, black_dict)
 
#############################################################################


WHITE = 'O'
BLACK = '@'
CORNER = 'X'
SPACE = '-'
MAX_PIECE = 20

ROW = 0     # index for coordinates
COL = 1     # index for coordinates

# indices for white_dict/black_dict
COORD = 0                   
IS_ALIVE = 1               
WHITES_ASSIGNED = 2                   
CLOSEST_WHITE = 0       
SECOND = 1  
TOTAL_DIST = 3             

# represent whether a piece is eliminated or not
STILL_ALIVE = True
DEAD = False

MAX_DIST = 20   # an arbitrary large number to initialise distance-related variables
FULL = 2        # if a black piece has 2 white pieces assigned to eliminate it
MAX_ROW = 8     # max number of rows
MAX_COL = 8     # max number of columns
MIN_ROW = 0     # min number of rows
MIN_COL = 0     # min number of columns
DICT_SIZE = 4   # for dict_init()

# to represent elimination method for first_kill() and can_eliminate()
NORMAL_ELIMNATION = 1
CORNER_ELIMINATION = 2
CANT_BE_ELIMINATED = 0

text = [None] * 9

# for indexing in Astar()
BOARD = 2
PATH = 3
SEQ_PATH = 4
LAST = -1

# for indexing seq_path
PREV = 0 
NEXT = 1

for row in range(0, 9):
    text[row] = input()

main(text)
