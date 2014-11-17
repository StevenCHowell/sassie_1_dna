import sh_box_of_balls
import numpy
from numpy import linalg 
import random

def evaluate_rotation_matrix(R,this_theta):
	pi = numpy.pi
	print
	print '###################'
	theta = this_theta * 180.0/pi
	condition_number = linalg.cond(R)
	smallest_singular_value = linalg.cond(R,-2)
	two_norm = linalg.cond(R,2)
	print 'theta = ',theta,'\n',R,'\ncondition number = ',condition_number
	print 'largest_singluar_value = ',two_norm
	print 'smallest_singular_value = ',smallest_singular_value
	print '###################'
	print 

	return

if __name__ == '__main__':

	pi = numpy.pi

	input_vector = [0.0,1.0,0.0]

	initial_theta = -160.0 * pi/180.0
	delta_theta = -5.0 * pi/180.0

	for i in xrange(1000):
		x = random.random()
		y = random.random()
		z = random.random()
		input_vector = [x,y,z]
		this_theta = initial_theta + i*delta_theta
		R = box_of_balls.get_rotation_matrix(input_vector,this_theta)
		evaluate_rotation_matrix(R,this_theta)

