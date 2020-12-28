import SudokuBoard
import Variable
import Domain
import Trail
import Constraint
import ConstraintNetwork
import time
import random
from collections import defaultdict

class BTSolver:

    # ==================================================================
    # Constructors
    # ==================================================================

    def __init__ ( self, gb, trail, val_sh, var_sh, cc ):
        self.network = ConstraintNetwork.ConstraintNetwork(gb)
        self.hassolution = False
        self.gameboard = gb
        self.trail = trail

        self.varHeuristics = var_sh
        self.valHeuristics = val_sh
        self.cChecks = cc

    # ==================================================================
    # Consistency Checks
    # ==================================================================

    # Basic consistency check, no propagation done
    def assignmentsCheck ( self ):
        for c in self.network.getConstraints():
            if not c.isConsistent():
                return False
        return True

    """
        Part 1 TODO: Implement the Forward Checking Heuristic

        This function will do both Constraint Propagation and check
        the consistency of the network

        (1) If a variable is assigned then eliminate that value from
            the square's neighbors.

        Note: remember to trail.push variables before you assign them
        Return: a tuple of a dictionary and a bool. The dictionary contains all MODIFIED variables, mapped to their MODIFIED domain.
                The bool is true if assignment is consistent, false otherwise.
    """
    def forwardChecking ( self ):
        d = {}
        b = True
        for av in self.network.variables:
            if av.isAssigned():
                for n in self.network.getNeighborsOfVariable(av):
                    if n.isChangeable and not n.isAssigned() and n.getDomain().contains(av.getAssignment()):
                        self.trail.push(n)
                        n.removeValueFromDomain(av.getAssignment())
                        d[n] = n.getDomain()
                        if n.domain.size() == 0:
                            b = False
                    '''
                    if n.domain.size() == 1:
                        self.trail.push(n)
                        assignedVars.append(n)
                    #'''
        #b = self.assignmentsCheck()
        return (d,b)

    # =================================================================
	# Arc Consistency
	# =================================================================
    def arcConsistency( self ):
        assignedVars = []
        for c in self.network.constraints:
            for v in c.vars:
                if v.isAssigned():
                    assignedVars.append(v)
        while len(assignedVars) != 0:
            av = assignedVars.pop(0)
            for neighbor in self.network.getNeighborsOfVariable(av):
                if neighbor.isChangeable and not neighbor.isAssigned() and neighbor.getDomain().contains(av.getAssignment()):
                    neighbor.removeValueFromDomain(av.getAssignment())
                    if neighbor.domain.size() == 1:
                        neighbor.assignValue(neighbor.domain.values[0])
                        assignedVars.append(neighbor)

    
    """
        Part 2 TODO: Implement both of Norvig's Heuristics

        This function will do both Constraint Propagation and check
        the consistency of the network

        (1) If a variable is assigned then eliminate that value from
            the square's neighbors.

        (2) If a constraint has only one possible place for a value
            then put the value there.

        Note: remember to trail.push variables before you assign them
        Return: a pair of a dictionary and a bool. The dictionary contains all variables 
		        that were ASSIGNED during the whole NorvigCheck propagation, and mapped to the values that they were assigned.
                The bool is true if assignment is consistent, false otherwise.
    """
    def norvigCheck ( self ): 
        d = {}
        b = True
        #'''
        self.forwardChecking()
        rows = defaultdict(list)
        cols = defaultdict(list)
        blocks = defaultdict(list)
        units = [rows, cols, blocks]
        for v in self.network.variables:
            rows[v.row].append(v)
            cols[v.col].append(v)
            blocks[v.block].append(v)
        for unit in units:
            for u in unit.values():
                counter = [0 for i in range(self.gameboard.N)]
                for i in range(self.gameboard.N):
                    for val in u[i].getDomain().values:
                        counter[val-1] += 1
                for j in range(self.gameboard.N):
                    if counter[j] == 1:
                        for x in range(9):
                            if j+1 in u[x].getDomain().values:
                                variable = u[x]
                                d[variable] = j+1
                                self.trail.push(variable)
                                variable.assignValue(j+1)
                                self.forwardChecking()
                                break
        b = self.assignmentsCheck()
        #'''
        return (d,b)

    """
         Optional TODO: Implement your own advanced Constraint Propagation

         Completing the three tourn heuristic will automatically enter
         your program into a tournament.
     """
    def getTournCC ( self ):
        return False

    # ==================================================================
    # Variable Selectors
    # ==================================================================

    # Basic variable selector, returns first unassigned variable
    def getfirstUnassignedVariable ( self ):
        for v in self.network.variables:
            if not v.isAssigned():
                return v

        # Everything is assigned
        return None

    """
        Part 1 TODO: Implement the Minimum Remaining Value Heuristic

        Return: The unassigned variable with the smallest domain
    """
    def getMRV ( self ):
        minD = 1000000 if len(self.network.variables) != 0 else 0
        for v in self.network.variables:
            d = v.getDomain().size()
            if not v.isAssigned() and d < minD:
                minD = int(d)
        for v in self.network.variables:
            if not v.isAssigned() and minD == v.getDomain().size():
                return v
        return None

    """
        Part 2 TODO: Implement the Degree Heuristic

        Return: The unassigned variable with the most unassigned neighbors
    """
    def getDegree ( self ):
        minVal = -1
        tup = ()
        for v in self.network.variables:
            count = 0
            for n in self.network.getNeighborsOfVariable(v):
                if not n.isAssigned():
                    count += 1
            if count > minVal:
                minVal = count
                tup = (v, minVal)
        return tup[0]

    """
        Part 2 TODO: Implement the Minimum Remaining Value Heuristic
                       with Degree Heuristic as a Tie Breaker

        Return: The unassigned variable with the smallest domain and affecting the  most unassigned neighbors.
                If there are multiple variables that have the same smallest domain with the same number of unassigned neighbors, add them to the list of Variables.
                If there is only one variable, return the list of size 1 containing that variable.
    """
    def MRVwithTieBreaker ( self ):
        variables = []
        d1 = defaultdict(int)
        # MRV
        minD = 1000000 if len(self.network.variables) != 0 else 0
        for v in self.network.variables:
            d = v.getDomain().size()
            if not v.isAssigned() and d < minD:
                minD = int(d)
        for v in self.network.variables:
            if not v.isAssigned() and minD == v.getDomain().size():
                variables.append(v)
        # degree heuristic
        maxN = -1
        for v in variables:
            for n in self.network.getNeighborsOfVariable(v):
                if not n.isAssigned():
                    d1[v] += 1
            maxN = d1[v] if d1[v] > maxN else maxN
        return [x for x, y in d1.items() if y==maxN]

    """
         Optional TODO: Implement your own advanced Variable Heuristic

         Completing the three tourn heuristic will automatically enter
         your program into a tournament.
     """
    def getTournVar ( self ):
        return None

    # ==================================================================
    # Value Selectors
    # ==================================================================

    # Default Value Ordering
    def getValuesInOrder ( self, v ):
        values = v.domain.values
        return sorted( values )

    """
        Part 1 TODO: Implement the Least Constraining Value Heuristic

        The Least constraining value is the one that will knock the least
        values out of it's neighbors domain.

        Return: A list of v's domain sorted by the LCV heuristic
                The LCV is first and the MCV is last
    """
    def getValuesLCVOrder ( self, v ):
        d = defaultdict(int)
        neighbors = self.network.getNeighborsOfVariable(v)
        for val in v.domain.values:
            for n in neighbors:
                if n.getDomain().contains(val) and not n.isAssigned():
                    d[val] += 1
        values = [x for x, y in sorted(d.items(), key = lambda item: item[1])]
        return values
    """
         Optional TODO: Implement your own advanced Value Heuristic

         Completing the three tourn heuristic will automatically enter
         your program into a tournament.
     """
    def getTournVal ( self, v ):
        return None

    # ==================================================================
    # Engine Functions
    # ==================================================================

    def solve ( self, time_left=600):
        if time_left <= 60:
            return -1

        start_time = time.time()
        if self.hassolution:
            return 0

        # Variable Selection
        v = self.selectNextVariable()

        # check if the assigment is complete
        if ( v == None ):
            # Success
            self.hassolution = True
            return 0

        # Attempt to assign a value
        for i in self.getNextValues( v ):

            # Store place in trail and push variable's state on trail
            self.trail.placeTrailMarker()
            self.trail.push( v )

            # Assign the value
            v.assignValue( i )

            # Propagate constraints, check consistency, recur
            if self.checkConsistency():
                elapsed_time = time.time() - start_time 
                new_start_time = time_left - elapsed_time
                if self.solve(time_left=new_start_time) == -1:
                    return -1
                
            # If this assignment succeeded, return
            if self.hassolution:
                return 0

            # Otherwise backtrack
            self.trail.undo()
        
        return 0

    def checkConsistency ( self ):
        if self.cChecks == "forwardChecking":
            return self.forwardChecking()[1]

        if self.cChecks == "norvigCheck":
            return self.norvigCheck()[1]

        if self.cChecks == "tournCC":
            return self.getTournCC()

        else:
            return self.assignmentsCheck()

    def selectNextVariable ( self ):
        if self.varHeuristics == "MinimumRemainingValue":
            return self.getMRV()

        if self.varHeuristics == "Degree":
            return self.getDegree()

        if self.varHeuristics == "MRVwithTieBreaker":
            return self.MRVwithTieBreaker()[0]

        if self.varHeuristics == "tournVar":
            return self.getTournVar()

        else:
            return self.getfirstUnassignedVariable()

    def getNextValues ( self, v ):
        if self.valHeuristics == "LeastConstrainingValue":
            return self.getValuesLCVOrder( v )

        if self.valHeuristics == "tournVal":
            return self.getTournVal( v )

        else:
            return self.getValuesInOrder( v )

    def getSolution ( self ):
        return self.network.toSudokuBoard(self.gameboard.p, self.gameboard.q)
