# This is a file that contains modularized elements of a stress fiber (sf) contraction simulation.

from dolfin import *
#import mshr
import ufl
import meshio
import numpy as np

'''
	This class overloads the dolfin "Expression" function.
	It defines the initial orientation of the sf distribution.

	Considerations: how do I allow for different distributions to be defined?
	answer: I should probably add an input for distribution type, and then it can be chosen. I think.
	but I would have to alter the class itself every time I want to add a new distribution; I guess that's ok.
'''

class FiberOrientation(UserExpression):
        '''
            User-defined stress fiber orientation in the cell
            in the undeformed configuration
        '''

        # still need to modularize this better/make it easier to adapt/use.

        #current expression: non-continuous, +z orientation when z-coordinate > 35, +y otherwise
        def eval(self, value, x):

            # undeformed; makes it easier.
            value[0] = 0
            value[1] = 0
            value[2] = 1

            

        def value_shape(self):
            return (3,)

# this represents the scalar contractile strength field that is to be applied to our mesh            
class ContractileStrength(UserExpression):
    
    # first spatially defined sf contraction will just say that the contractile strength is large far from nucleus
    # will it also apply to the gel? I'll understand this more by understanding the math more; but need to 

    # overloading the __init__ function in order to allow time to be tracked.
    def __init__(self, t,**kwargs):

        '''
        Old code/method -- no idea if it's faster or slower; both work.

        # call super class to pass other objects defined here
        #super(ContractileStrength,self).__init__(**kwargs)
        #self.t = t
        #self._element = kwargs["element"] --this line is garbage/not needed
        '''

        self.t = t
        super().__init__(**kwargs)

    
    def eval(self, value, x):
        
        # define contractile strenght spatially
        #value = self.t*5 -- old method that wasn't working; "value" is vector valued, not an int.
        value[0] = self.t*5

    def value_shape(self):
        return (1,)

'''

	Notes, and previous definitions of fiber orientation


	
            x_c = x[0]
            y_c = x[1]
            z_c = x[2]        
            value[0] = x_c
            value[1] = y_c
            value[2] = z_c

            #normalizing/turning into unit vector
            for i in range(0,3):
                value[i] = value[i]/sqrt(x_c**2 + y_c**2 + z_c**2) #be careful of when this is small/division by zero.
			 

			
            # notes: (THESE NOTES ARE IMPORTANT)

            # defining contractile strength as a function of position

            # michael mentioned using spherical harmonics to specify sf,
            # general ideas: map geometry of cell to sphere, on the sphere we 
            # prescribe the fiber directions and then map it back to the cell.
            # (try to meet with Michael about this to understand exactly what he is asking)

            # can try something more complicated with real cell geometry. 
            # method: can we combine sparse info about where fiber SHOULD go to better describe the distribution

            #1. how do we know direction? (have some idea around protrusions)
            #2. how do we describe continuous distribution?

            # think of ways to specify density of sf here
'''