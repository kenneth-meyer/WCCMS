# This is a file that contains modularized elements of a stress fiber (sf) contraction simulation.

from dolfin import *


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
        #current expression: non-continuous, +z orientation when z-coordinate > 35, +y otherwise
        def eval(self, value, x):

            # undeformed; makes it easier.
            value[0] = 0
            value[1] = 0
            value[2] = 1

            '''
            x_c = x[0]
            y_c = x[1]
            z_c = x[2]        
            value[0] = x_c
            value[1] = y_c
            value[2] = z_c

            #normalizing/turning into unit vector
            for i in range(0,3):
                value[i] = value[i]/sqrt(x_c**2 + y_c**2 + z_c**2) #be careful of when this is small/division by zero.
			''' 
            # notes:

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


        def value_shape(self):
            return (3,)

# this represents the scalar contractile strength field that is to be applied to our mesh            
class ContractileStrength(UserExpression): # will it recognize this as a Expression function in dolfin?

    #t = 0. #time-dependent as well...
    
    # first spatially defined sf contraction will just say that the contractile strength is large far from nucleus
    # will it also apply to the gel? I'll understand this more by understanding the math more; but need to 

    #t = 0.

    # overloading the __init__ function in order to allow time to be tracked.
    def __init__(self, **kwargs):
         self._t = kwargs["t"]
         self._element = kwargs["element"]
         #self.element = element


    def eval(self, value, x):
        
        # define contractile strenght spatially

        # start off with changing it with z position or something.

        '''
        x_c = x[0]
        y_c = x[1]
        z_c = x[2]


        # contractile strength increases as we move further away from the nucleus.
        value = t*sqrt(x_c**2 + y_c**2 + z_c**2)/100 #divide by 100 because I don't want it to be too large.
        '''
        
        value = t*5 # t should be advanced; it is a field.


    def value_shape(self):
        return (1,)
        # not sure how to define the shape/how to access the shape of the stress field.
